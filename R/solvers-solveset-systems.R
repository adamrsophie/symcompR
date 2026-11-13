# Port of the systems part of sympy/solvers/solveset.py: linear_eq_to_matrix,
# linsolve, substitution and nonlinsolve, with linear_coeffs from
# sympy/core/exprtools.py and _linsolve from sympy/polys/matrices/linsolve.py.

# ---------------------------------------------------------------------------
# Helpers

# Python's type name as SymPy's messages print it
.sset_pytype <- function(o) {
  if (is.numeric(o)) return(if (all(o == round(o))) "<class 'int'>" else "<class 'float'>")
  if (is.character(o)) return("<class 'str'>")
  if (!is_sym(o)) return("<class 'list'>")
  switch(o$type,
         Symbol = "<class 'sympy.core.symbol.Symbol'>",
         Dummy = "<class 'sympy.core.symbol.Dummy'>",
         Integer = "<class 'sympy.core.numbers.Integer'>",
         Rational = "<class 'sympy.core.numbers.Rational'>",
         Add = "<class 'sympy.core.add.Add'>",
         Mul = "<class 'sympy.core.mul.Mul'>",
         Pow = "<class 'sympy.core.power.Pow'>",
         Function = o$name,
         paste0("<class '", o$type, "'>"))
}

# a dict given as a list of (key, value) pairs or as a named list
.sset_as_dict <- function(r, syms) {
  if (is.list(r) && identical(names(r), c("k", "v"))) return(r)
  d <- .sd_new()
  if (!length(r)) return(d)
  nms <- names(r)
  for (i in seq_along(r)) {
    if (!is.null(nms) && nzchar(nms[i])) {
      j <- match(nms[i], vapply(syms, function(s) s$name, ""))
      k <- if (is.na(j)) Symbol(nms[i]) else syms[[j]]
      d <- .sd_set(d, k, as_sym(r[[i]]))
    } else d <- .sd_set(d, as_sym(r[[i]][[1L]]), as_sym(r[[i]][[2L]]))
  }
  d
}

# Python dict equality
.sset_dict_eq <- function(a, b) length(a$k) == length(b$k) && identical(.sd_ident(a), .sd_ident(b))
.sset_dict_in <- function(d, lst) any(vapply(lst, function(z) .sset_dict_eq(d, z), NA))

# Basic.dummy_eq: equality up to the names of bound and dummy variables
.sset_alpha_eq <- function(a, b) {
  if (identical(a$key, b$key)) return(TRUE)
  if (!is_sym(a) || !is_sym(b) || a$type != b$type || length(a$args) != length(b$args)) return(FALSE)
  if (!length(a$args)) return(FALSE)
  if (a$type == "Lambda") {
    va <- a$args[[1L]]$args; vb <- b$args[[1L]]$args
    if (length(va) != length(vb)) return(FALSE)
    b2 <- xreplace(b$args[[2L]], Map(list, vb, va))
    return(.sset_alpha_eq(a$args[[2L]], b2))
  }
  if (a$type == "ConditionSet") {
    b2 <- xreplace(b$args[[2L]], list(list(b$args[[1L]], a$args[[1L]])))
    return(.sset_alpha_eq(a$args[[2L]], b2) && .sset_alpha_eq(a$args[[3L]], b$args[[3L]]))
  }
  all(mapply(.sset_alpha_eq, a$args, b$args))
}
.sset_free_dummies <- function(e) Filter(function(s) s$type == "Dummy", free_symbols(e))
.sset_dummy_eq_any <- function(a, b) {
  fd <- .sset_free_dummies(a)
  if (length(fd) != 1L) return(.sset_alpha_eq(a, b))
  fs <- free_symbols(b)
  if (length(fs) != 1L) return(.sset_alpha_eq(a, b))
  tmp <- Dummy()
  .sset_alpha_eq(xreplace(a, list(list(fd[[1L]], tmp))), xreplace(b, list(list(fs[[1L]], tmp))))
}

.sset_eq_expr <- function(e) {
  e <- as_sym(e)
  if (e$type == "Equality") .sub(e$args[[1L]], e$args[[2L]]) else e
}
.sset_seq <- function(x) if (inherits(x, "symmatrix")) x$data else if (is_sym(x) && x$type == "Tuple") x$args else x

# ---------------------------------------------------------------------------
# linear_coeffs and linear_eq_to_matrix

# _linear_eq_to_dict: list(coefficient dicts, constant terms)
.sset_linear_eq_to_dict <- function(eqs, syms) {
  coeffs <- list(); ind <- list()
  for (e in eqs) {
    if (e$type == "Equality") {
      l <- .pv_lin_eq2dict(e$args[[1L]], syms)
      rr <- .pv_lin_eq2dict(e$args[[2L]], syms)
      tk <- l[[2L]]$keys; tv <- l[[2L]]$values
      for (i in seq_along(rr[[2L]]$keys)) {
        k <- rr[[2L]]$keys[[i]]
        j <- which(vapply(tk, function(z) identical(z$key, k$key), TRUE))
        if (length(j)) tv[[j]] <- .sub(tv[[j]], rr[[2L]]$values[[i]])
        else { tk[[length(tk) + 1L]] <- k; tv[[length(tv) + 1L]] <- .neg(rr[[2L]]$values[[i]]) }
      }
      keep <- !vapply(tv, is_zero_num, TRUE)
      coeffs[[length(coeffs) + 1L]] <- .pydict(tk[keep], tv[keep])
      ind[[length(ind) + 1L]] <- .sub(l[[1L]], rr[[1L]])
    } else {
      r <- .pv_lin_eq2dict(e, syms)
      coeffs[[length(coeffs) + 1L]] <- r[[2L]]
      ind[[length(ind) + 1L]] <- r[[1L]]
    }
  }
  list(coeffs, ind)
}

#' Coefficients of a linear expression
#'
#' `linear_coeffs(eq, ...)` returns the coefficients of the symbols `...` in
#' the linear expression `eq` followed by the constant term, as SymPy's
#' `linear_coeffs()`; with `dict = TRUE` it returns a dictionary keyed by the
#' symbols present, with key `1` for the constant. Nonlinear terms raise a
#' `NonlinearError`.
#'
#' @param eq An expression.
#' @param ... The symbols.
#' @param dict Return a dictionary.
#' @return A list, or a dictionary.
#' @examples
#' x <- Symbol("x"); y <- Symbol("y")
#' linear_coeffs(3*x + 2*y - 1, x, y)
#' @export
linear_coeffs <- function(eq, ..., dict = FALSE) {
  eq <- as_sym(eq)
  syms <- list(...)
  if (length(syms) == 1L && is.list(syms[[1L]]) && !is_sym(syms[[1L]]))
    .sset_err("ValueError", "expecting unpacked symbols, *syms")
  syms <- lapply(syms, as_sym)
  if (anyDuplicated(vapply(syms, .k, ""))) .sset_err("ValueError", "duplicate symbols given")
  r <- .sset_linear_eq_to_dict(list(eq), syms)
  d <- r[[1L]][[1L]]; c <- r[[2L]][[1L]]
  if (dict) {
    if (!is_zero_num(c)) d <- .pydict(c(d$keys, list(S_One)), c(d$values, list(c)))
    return(d)
  }
  rv <- rep(list(S_Zero), length(syms) + 1L)
  rv[[length(rv)]] <- c
  keys <- vapply(d$keys, .k, "")
  for (i in seq_along(syms)) {
    j <- match(syms[[i]]$key, keys)
    if (!is.na(j)) rv[[i]] <- d$values[[j]]
  }
  rv
}

#' Matrix form of a linear system
#'
#' `linear_eq_to_matrix(equations, ...)` writes the linear equations (or
#' expressions equated to zero) as `A x = b` for the given symbols and returns
#' the pair `(A, b)`, as SymPy's `linear_eq_to_matrix()`. The symbols may be
#' given one by one or as a list; nonlinear terms raise a `NonlinearError`.
#'
#' @param equations A list of equations or expressions, one of them, or a
#'   matrix of expressions.
#' @param ... The symbols, or one list of them.
#' @return A tuple of the matrices `A` and `b`.
#' @examples
#' x <- Symbol("x"); y <- Symbol("y")
#' linear_eq_to_matrix(list(x + 2*y - 3, 3*x - y), list(x, y))
#' @export
linear_eq_to_matrix <- function(equations, ...) {
  symbols <- list(...)
  if (!length(symbols))
    .sset_err("ValueError", .filldedent("Symbols must be given, for which coefficients are to be found."))
  if (is.list(symbols[[1L]]) && !is_sym(symbols[[1L]])) symbols <- symbols[[1L]]
  else if (is_sym(symbols[[1L]]) && symbols[[1L]]$type == "Tuple") symbols <- symbols[[1L]]$args
  symbols <- lapply(symbols, as_sym)
  if (anyDuplicated(vapply(symbols, .k, ""))) .sset_err("ValueError", "Symbols must be unique")
  if (inherits(equations, "symmatrix")) equations <- equations$data
  else if (is_sym(equations) && equations$type == "Tuple") equations <- equations$args
  else if (is_sym(equations) || !is.list(equations)) {
    e <- as_sym(equations)
    if (!(.sset_is_expr(e) || e$type == "Equality"))
      .sset_err("ValueError", .filldedent("Equation(s) must be given as a sequence, Expr, Eq or Matrix."))
    equations <- list(e)
  }
  equations <- lapply(equations, as_sym)
  r <- .sset_linear_eq_to_dict(equations, symbols)
  eq <- r[[1L]]; cs <- r[[2L]]
  n <- length(eq); m <- length(symbols)
  keys <- vapply(symbols, .k, "")
  A <- rep(list(S_Zero), n * m)
  for (row in seq_len(n)) {
    d <- eq[[row]]
    for (i in seq_along(d$keys)) A[[(row - 1L) * m + match(d$keys[[i]]$key, keys)]] <- d$values[[i]]
  }
  .pytuple(list(.new_mat(A, n, m), .new_mat(lapply(cs, .neg), n, 1L)))
}

# ---------------------------------------------------------------------------
# linsolve

# sympy.polys.matrices.linsolve._linsolve: a dict of the solution, or NULL
.sset__linsolve <- function(eqs, syms) {
  nsyms <- length(syms)
  r <- .sset_linear_eq_to_dict(eqs, syms)
  eqsdict <- r[[1L]]; const <- r[[2L]]
  # sympy_dict_to_dm
  elems <- c(const, unlist(lapply(eqsdict, function(d) d$values), recursive = FALSE))
  elems <- elems[!duplicated(vapply(elems, .k, ""))]
  cd <- construct_domain(elems, field = TRUE, extension = TRUE)
  K <- cd[[1L]]
  conv <- cd[[2L]]
  if (!is.list(conv) || inherits(conv, "sym") || length(elems) == 1L && !is.list(conv)) conv <- list(conv)
  ekeys <- vapply(elems, .k, "")
  emap <- function(e) conv[[match(e$key, ekeys)]]
  skeys <- vapply(syms, .k, "")
  rows <- list()
  for (i in seq_along(eqsdict)) {
    d <- eqsdict[[i]]; rhs <- const[[i]]
    row <- list()
    for (j in seq_along(d$keys)) row[[.sk(match(d$keys[[j]]$key, skeys) - 1L)]] <- emap(d$values[[j]])
    if (!is_zero_num(rhs)) row[[.sk(nsyms)]] <- -emap(rhs)
    if (length(row)) rows[[length(rows) + 1L]] <- row
  }
  names(rows) <- .sk(seq_along(rows) - 1L)
  if (!length(rows)) {
    sol <- .sd_new()
    for (s in syms) sol <- .sd_set(sol, s, s)
    return(sol)
  }
  # ponytail: sdm_irref also for float matrices; SymPy uses ddm_rref there
  rr <- .sdm_irref(rows)
  Arref <- rr$rref; pivots <- rr$pivots; nzcols <- rr$nzc
  if (length(pivots) && pivots[length(pivots)] == nsyms) return(NULL)
  P <- .sdm_particular_from_rref(Arref, nsyms + 1L, pivots)
  ns <- if (!length(pivots)) list(K = lapply(seq_len(nsyms) - 1L, function(i) { z <- list(); z[[.sk(i)]] <- .de_one(K); z }),
                                  nonpivots = seq_len(nsyms) - 1L)
        else .sdm_nullspace_from_rref_core(Arref, .de_one(K), nsyms, pivots, nzcols)
  terms <- list()
  add_term <- function(i, t) {
    k <- .sk(i)
    terms[[k]] <<- c(terms[[k]], list(t))
  }
  for (k in names(P)) add_term(as.integer(k), .de_to_sympy(K, P[[k]]))
  for (q in seq_along(ns$nonpivots)) {
    sym <- syms[[ns$nonpivots[q] + 1L]]
    Vi <- ns$K[[q]]
    for (k in names(Vi)) add_term(as.integer(k), .mul(list(sym, .de_to_sympy(K, Vi[[k]]))))
  }
  sol <- .sd_new()
  for (i in seq_len(nsyms)) {
    t <- terms[[.sk(i - 1L)]]
    sol <- .sd_set(sol, syms[[i]], if (is.null(t)) S_Zero else .add(t))
  }
  sol
}

#' Solve a linear system as a set
#'
#' `linsolve(system, ...)` solves a system of linear equations and returns
#' the `FiniteSet` of one `Tuple` (or the empty set), as SymPy's `linsolve()`.
#' Free unknowns stay as parameters. The system may be a list of equations or
#' expressions, a list `(A, b)` of matrices, or an augmented matrix; the
#' symbols follow one by one or as a list, and for matrix input they default
#' to `tau0, tau1, ...`.
#'
#' @param system The linear system.
#' @param ... The symbols, or one list of them.
#' @return A set of tuples.
#' @examples
#' x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
#' linsolve(list(x + y + z - 1, x + y + 2*z - 3), list(x, y, z))
#' linsolve(Matrix(list(list(1, 2, 3), list(4, 5, 6))))
#' @export
linsolve <- function(system, ...) {
  symbols <- list(...)
  if (!is_sym(system) && !length(system) && !inherits(system, "symmatrix")) return(S_EmptySet)
  if (length(symbols) && is.list(symbols[[1L]]) && !is_sym(symbols[[1L]]) && !inherits(symbols[[1L]], "symmatrix"))
    symbols <- symbols[[1L]]
  else if (length(symbols) && is_sym(symbols[[1L]]) && symbols[[1L]]$type == "Tuple") symbols <- symbols[[1L]]$args
  dup_msg <- "duplicate symbols given"
  A <- NULL; b <- NULL
  if (!inherits(system, "symmatrix") && ((is.list(system) && !is_sym(system)) || (is_sym(system) && system$type == "Tuple"))) {
    system <- .sset_seq(system)
    if (length(system) == 2L && inherits(system[[1L]], "symmatrix")) { A <- system[[1L]]; b <- system[[2L]] }
    if (!inherits(system[[1L]], "symmatrix")) {
      if (!length(symbols))
        .sset_err("ValueError", .filldedent("When passing a system of equations, the explicit symbols for which a solution is being sought must be given as a sequence, too."))
      symbols <- lapply(symbols, as_sym)
      if (anyDuplicated(vapply(symbols, .k, ""))) .sset_err("ValueError", dup_msg)
      eqs <- lapply(system, as_sym)
      sol <- .sset__linsolve(eqs, symbols)
      if (is.null(sol)) return(S_EmptySet)
      return(FiniteSet(do.call(Tuple, lapply(symbols, function(s) { v <- .sd_get(sol, s); if (is.null(v)) s else v }))))
    }
  } else if (inherits(system, "symmatrix") && !(length(symbols) && inherits(symbols[[1L]], "symmatrix"))) {
    nc <- system$ncol
    A <- .from_rows(lapply(.mrows(system), function(r) r[-nc]))
    b <- .from_rows(lapply(.mrows(system), function(r) r[nc]))
  }
  if (is.null(b)) .sset_err("ValueError", "Invalid arguments")
  gen <- NULL
  if (!length(symbols)) {
    symbols <- lapply(seq_len(A$ncol), function(i) Dummy())
    names_ab <- unique(unlist(lapply(c(A$data, b$data), function(e) vapply(free_symbols(e), function(s) s$name, ""))))
    nm <- "tau"
    strip <- function(s) sub("[0-9]+$", "", s)
    while (any(nm == strip(names_ab))) nm <- .sset_numbered_incr(nm)
    gen <- nm
  } else symbols <- lapply(symbols, as_sym)
  eqs <- list()
  rows <- .mrows(A)
  for (i in seq_along(rows)) {
    terms <- list()
    for (j in seq_along(rows[[i]])) if (!is_zero_num(rows[[i]][[j]])) terms[[length(terms) + 1L]] <- .mul(list(rows[[i]][[j]], symbols[[j]]))
    terms[[length(terms) + 1L]] <- .neg(b$data[[i]])
    eqs[[length(eqs) + 1L]] <- .add(terms)
  }
  polys <- lapply(eqs, function(e) .sp_as_poly(e, symbols))
  if (any(vapply(polys, is.null, NA))) .sset_err("NonlinearError", "Nonlinear term encountered in solve_lin_sys")
  sol <- .sp_solve_lin_sys(polys, symbols)
  if (is.null(sol)) return(S_EmptySet)
  tup <- do.call(Tuple, lapply(symbols, function(s) { v <- .sd_get(sol, s); if (is.null(v)) s else v }))
  if (!is.null(gen)) {
    fk <- vapply(free_symbols(tup), .k, "")
    rep <- list(); i <- 0L
    for (s in symbols) if (s$key %in% fk) { rep[[length(rep) + 1L]] <- list(s, Symbol(paste0(gen, i))); i <- i + 1L }
    if (length(rep)) tup <- subs(tup, rep)
  }
  FiniteSet(tup)
}

# numbered_string_incr from uniquely_named_symbol
.sset_numbered_incr <- function(s) {
  m <- regmatches(s, regexpr("[0-9]*$", s))
  base <- substr(s, 1L, nchar(s) - nchar(m))
  paste0(base, if (nzchar(m)) as.integer(m) + 1L else 0L)
}

# ---------------------------------------------------------------------------
# substitution

.sset_return_conditionset <- function(eqs, symbols, exclude = list()) {
  conds <- c(lapply(eqs, function(lhs) Eq(lhs, 0)), lapply(exclude, function(d) Ne(d, 0)))
  ConditionSet(do.call(Tuple, symbols), do.call(And, conds),
               do.call(ProductSet, rep(list(Complexes), length(symbols))))
}

.sset_process_exclusions_for_branch <- function(branch, exclude, symbols) {
  complements <- .sd_new()
  conditions <- list()
  isexpr <- vapply(branch$v, function(v) !is_set(v), NA)
  expr_branch <- list(k = branch$k[isexpr], v = branch$v[isexpr])
  set_valued <- vapply(branch$k[!isexpr], .k, "")
  skeys <- vapply(symbols, .k, "")
  for (denominator in exclude) {
    transformed <- cancel(.sd_subs(denominator, expr_branch))
    z <- .sset_fact(transformed, "zero")
    if (isTRUE(z)) return(NULL)
    if (identical(z, FALSE)) next
    fs <- free_symbols(transformed)
    fk <- vapply(fs, .k, "")
    unresolved <- fs[fk %in% skeys]
    if (length(unresolved) == 1L && !any(fk %in% set_valued)) {
      parameter <- unresolved[[1L]]
      bv <- .sd_get(branch, parameter)
      if (is.null(bv) || .sset_eq(bv, parameter)) {
        forbidden <- if (isTRUE(.sset_fact(parameter, "real"))) solveset_real(transformed, parameter)
                     else solveset_complex(transformed, parameter)
        if (!.has_type(forbidden, "ConditionSet")) {
          old <- .sd_get(complements, parameter)
          complements <- .sd_set(complements, parameter, Union(if (is.null(old)) S_EmptySet else old, forbidden))
          next
        }
      }
    }
    conditions[[length(conditions) + 1L]] <- Ne(transformed, 0)
  }
  list(complements, conditions)
}

.sset_branch_as_conditionset <- function(branch, symbols, intersections, complements, conditions) {
  relations <- conditions
  values <- list()
  for (symbol in symbols) {
    value <- .sd_get(branch, symbol)
    if (is.null(value)) value <- symbol
    values[[length(values) + 1L]] <- value
    if (is_set(value)) {
      relations[[length(relations) + 1L]] <- Contains(symbol, value)
      restricted_value <- symbol
    } else {
      if (!.sset_eq(value, symbol)) relations[[length(relations) + 1L]] <- Eq(symbol, value)
      restricted_value <- value
    }
    iv <- .sd_get(intersections, symbol)
    if (!is.null(iv)) relations[[length(relations) + 1L]] <- .contains(iv, restricted_value)
    cv <- .sd_get(complements, symbol)
    if (!is.null(cv)) relations[[length(relations) + 1L]] <- Not(.contains(cv, restricted_value))
  }
  domain <- do.call(ProductSet, lapply(symbols, function(s) if (isTRUE(.sset_fact(s, "real"))) Reals else Complexes))
  condition <- if (length(relations)) do.call(And, relations) else S_true
  if (all(vapply(values, function(v) !is_set(v) && .expr_is_number(v), NA))) {
    domain <- Intersection(domain, FiniteSet(do.call(Tuple, values)))
    condition <- subs(condition, lapply(seq_along(symbols), function(i) list(symbols[[i]], values[[i]])), simultaneous = TRUE)
  }
  ConditionSet(do.call(Tuple, symbols), condition, domain)
}

#' Solve a system of equations by substitution
#'
#' `substitution(system, symbols, result, known_symbols, exclude,
#' all_symbols)` solves the equations of `system` one after the other for the
#' unknowns `symbols`, substituting the values found so far, as SymPy's
#' `substitution()`. It returns a `FiniteSet` of `Tuple`s ordered like
#' `all_symbols` (by default `symbols`); unknowns left free stay symbolic, and
#' periodic solutions are `ImageSet`s.
#'
#' @param system A list of equations or expressions.
#' @param symbols A list of the unknowns.
#' @param result A list of known partial solutions, each a list of
#'   `list(symbol, value)` pairs or a named list.
#' @param known_symbols Symbols whose values are known.
#' @param exclude Expressions that must not vanish, such as denominators.
#' @param all_symbols The symbols to report, in order.
#' @return A set of tuples, or a `ConditionSet`.
#' @examples
#' x <- Symbol("x"); y <- Symbol("y")
#' substitution(list(x + y), list(x), list(list(list(y, 1))), list(y), list(), list(x, y))
#' substitution(list(x^2 - 1, y - x), list(x, y))
#' @export
substitution <- function(system, symbols, result = list(list()), known_symbols = list(),
                         exclude = list(), all_symbols = NULL) {
  if (missing(symbols)) .sset_err("TypeError", "substitution() missing 1 required positional argument: 'symbols'")
  system <- .sset_seq(system)
  if (is_sym(system)) system <- list(system)
  if (!length(system)) return(S_EmptySet)
  system <- lapply(system, .sset_eq_expr)
  if (!length(symbols) || identical(symbols, 0) || identical(symbols, FALSE))
    .sset_err("ValueError", .filldedent("Symbols must be given, for which solution of the system is to be found."))
  if (!is.list(symbols) && !(is_sym(symbols) && symbols$type == "Tuple") || (is_sym(symbols) && symbols$type != "Tuple"))
    .sset_err("TypeError", .filldedent(sprintf("symbols should be given as a sequence, e.g. a list.Not type %s: %s",
                                               .sset_pytype(symbols), if (is_sym(symbols)) .str(symbols) else format(symbols))))
  symbols <- lapply(.sset_seq(symbols), as_sym)
  if (!(symbols[[1L]]$type %in% c("Symbol", "Dummy")))
    .sset_err("ValueError", .filldedent(sprintf("Iterable of symbols must be given as second argument, not type %s: %s",
                                                .sset_pytype(symbols[[1L]]), .str(symbols[[1L]]))))
  all_symbols <- if (is.null(all_symbols)) symbols else lapply(.sset_seq(all_symbols), as_sym)
  known_symbols <- lapply(.sset_seq(known_symbols), as_sym)
  exclude <- lapply(.sset_seq(exclude), as_sym)
  dsyms <- c(all_symbols, symbols, known_symbols)
  result <- lapply(result, .sset_as_dict, syms = dsyms)
  .sset_substitution(system, symbols, result, known_symbols, exclude, all_symbols)
}

.sset_substitution <- function(system, symbols, result, known_symbols, exclude, all_symbols) {
  old_result <- result
  st <- new.env()
  st$complements <- .sd_new()
  st$intersections <- .sd_new()
  total_conditionset <- -1L
  total_solveset_call <- -1L
  kkeys <- vapply(known_symbols, .k, "")
  akeys <- vapply(all_symbols, .k, "")
  unsolved_syms <- function(eq, sort = FALSE) {
    free <- free_symbols(eq)
    u <- Filter(function(s) !(s$key %in% kkeys) && s$key %in% akeys, free)
    if (sort) u <- .sp_sorted(u)
    u
  }
  eqs_in_better_order <- if (length(system) < 2L) system else
    .key_sort(lapply(system, function(e) list(length(unsolved_syms(e)), .nodes(e), .sort_key(e))), system)

  add_intersection_complement <- function(result, intersection_dict, complement_dict) {
    final_result <- list()
    for (res in result) {
      res_copy <- res
      for (i in seq_along(res$k)) {
        key_res <- res$k[[i]]; value_res <- res$v[[i]]
        intersect_set <- .sd_get(intersection_dict, key_res)
        complement_set <- .sd_get(complement_dict, key_res)
        truthy <- function(s) !is.null(s) && !.sset_is_empty(s)
        if (truthy(intersect_set) || truthy(complement_set)) {
          new_value <- FiniteSet(value_res)
          if (truthy(intersect_set) && !.sset_set_is(intersect_set, "Complexes"))
            new_value <- Intersection(new_value, intersect_set)
          if (truthy(complement_set)) new_value <- Complement(new_value, complement_set)
          if (.sset_is_empty(new_value)) { res_copy <- NULL; break }
          if (.sset_is_fs(new_value) && length(new_value$args) == 1L) res_copy <- .sd_set(res_copy, key_res, new_value$args[[1L]])
          else res_copy <- .sd_set(res_copy, key_res, new_value)
        }
      }
      if (!is.null(res_copy)) final_result[[length(final_result) + 1L]] <- res_copy
    }
    final_result
  }

  extract_main_soln <- function(sym, sol, soln_imageset) {
    if (sol$type == "ConditionSet") sol <- sol$args[[3L]]
    if (sol$type == "Complement") {
      st$complements <- .sd_set(st$complements, sym, sol$args[[2L]])
      sol <- sol$args[[1L]]
    }
    if (sol$type == "Union") {
      sol_args <- .set_print_order(sol$args)
      sol <- S_EmptySet
      for (a2 in sol_args) sol <- Union(sol, if (.sset_is_fs(a2)) a2 else FiniteSet(a2))
    }
    if (sol$type == "Intersection") {
      a1 <- sol$args[[1L]]
      if (!(.sset_set_is(a1, "Reals") || .sset_set_is(a1, "Complexes")))
        st$intersections <- .sd_set(st$intersections, sym, a1)
      sol <- sol$args[[2L]]
    }
    if (sol$type == "ImageSet") {
      expr2 <- sol$args[[1L]]$args[[2L]]
      soln_imageset <- .sd_set(soln_imageset, expr2, sol)
      sol <- FiniteSet(expr2)
    }
    if (!.sset_is_fs(sol) && !.sset_is_empty(sol)) sol <- FiniteSet(sol)
    list(sol, soln_imageset)
  }

  restore_imgset <- function(rnew, original_imageset) {
    for (i in seq_along(rnew$k)) {
      img <- .sd_get(original_imageset, rnew$k[[i]])
      if (!is.null(img)) rnew$v[[i]] <- img
    }
    rnew
  }

  append_eq <- function(eq, res, n = NULL) {
    u <- Dummy("u")
    if (!is.null(n)) eq <- subs(eq, n, S_Zero)
    satisfy <- if (.is_true(eq)) TRUE else if (.is_false(eq)) FALSE else .sp_checksol(u, .sd1(u, eq), list(minimal = TRUE))
    if (identical(satisfy, FALSE)) list(res = .sd_new(), delete = TRUE, appended = FALSE)
    else list(res = res, delete = FALSE, appended = TRUE)
  }

  # _append_new_soln: list(newresult, delete_soln)
  append_new_soln <- function(rnew, sym, sol, imgset_yes, soln_imageset, original_imageset, newresult, eq = NULL) {
    delete_soln <- FALSE
    appended <- FALSE
    push <- function(nr, r) { if (!.sset_dict_in(r, nr)) nr[[length(nr) + 1L]] <- r; nr }
    if (!is.null(imgset_yes)) {
      local_n <- imgset_yes[[1L]]
      base <- imgset_yes[[2L]]
      if (!is.null(sym) && !is.null(sol) && !is_zero_num(sol)) {
        dl <- atoms(sol, "Dummy")
        lam <- Lambda(local_n, if (length(dl)) subs(sol, lapply(dl, function(d) list(d, local_n))) else sol)
        rnew <- .sd_set(rnew, sym, ImageSet(lam, base))
      }
      if (!is.null(eq)) {
        ae <- append_eq(eq, rnew, local_n)
        rnew <- ae$res; delete_soln <- ae$delete
        if (ae$appended) { newresult[[length(newresult) + 1L]] <- rnew; appended <- TRUE }
      }
    } else if (!is.null(eq)) {
      ae <- append_eq(eq, rnew)
      rnew <- ae$res; delete_soln <- ae$delete
      if (ae$appended) { newresult[[length(newresult) + 1L]] <- rnew; appended <- TRUE }
    } else if (!is.null(.sd_get(soln_imageset, sol))) {
      rnew <- .sd_set(rnew, sym, .sd_get(soln_imageset, sol))
      rnew <- restore_imgset(rnew, original_imageset)
      newresult <- push(newresult, rnew)
      appended <- TRUE
    } else {
      newresult[[length(newresult) + 1L]] <- rnew
      appended <- TRUE
    }
    rnew <- restore_imgset(rnew, original_imageset)
    # the appended dict is the same object as rnew in SymPy
    if (appended) newresult[[length(newresult)]] <- rnew
    else newresult <- push(newresult, rnew)
    list(newresult, delete_soln)
  }

  # list.remove with Python's iteration semantics: drop the first equal dict
  remove_res <- function(lst, res) {
    j <- which(vapply(lst, function(z) .sset_dict_eq(z, res), NA))
    if (length(j)) lst[-j[1L]] else lst
  }

  solve_using_known_values <- function(result, real) {
    solver <- if (real) solveset_real else solveset_complex
    soln_imageset <- .sd_new()
    total_solvest_call <- 0L
    total_conditionst <- 0L
    input <- result
    aliased <- TRUE
    for (index in seq_along(eqs_in_better_order)) {
      eq <- eqs_in_better_order[[index]]
      newresult <- list()
      imgset_yes <- NULL
      idx <- 1L
      while (idx <= length(result)) {
        res <- result[[idx]]
        pos <- idx
        idx <- idx + 1L
        original_imageset <- .sd_new()
        got_symbol <- list()
        for (i in seq_along(res$k)) {
          v <- res$v[[i]]
          if (is_sym(v) && v$type == "ImageSet") {
            res$v[[i]] <- v$args[[1L]]$args[[2L]]
            original_imageset <- .sd_set(original_imageset, res$k[[i]], v)
            dn <- atoms(v$args[[1L]]$args[[2L]], "Dummy")
            imgset_yes <- list(dn[[1L]], v$args[[2L]])
          }
        }
        result[[pos]] <- res
        if (aliased) input <- result
        eq2 <- expand(.sd_subs(eq, res))
        if (!is.null(imgset_yes) && !.has1(eq2, imgset_yes[[1L]])) imgset_yes <- NULL
        us <- unsolved_syms(eq2, sort = TRUE)
        if (!length(us)) {
          if (length(res$k)) {
            r <- append_new_soln(res, NULL, NULL, imgset_yes, soln_imageset, original_imageset, newresult, eq2)
            newresult <- r[[1L]]
            if (r[[2L]]) {
              result <- remove_res(result, res)
              if (aliased) input <- result
            }
          }
          next
        }
        dep <- .sp_indep(eq2, us)
        if ((.sset_has_type(dep[[1L]], "Abs") || .sset_has_type(dep[[2L]], "Abs")) && !real)
          .sset_err("NotImplementedError", "nonlinsolve cannot solve equations with Abs in the complex domain")
        soln_imageset <- .sd_new()
        for (sym in us) {
          not_solvable <- FALSE
          soln <- tryCatch({
            soln <- solver(eq2, sym)
            total_solvest_call <- total_solvest_call + 1L
            soln_new <- S_EmptySet
            if (soln$type == "Complement") {
              st$complements <- .sd_set(st$complements, sym, soln$args[[2L]])
              soln <- soln$args[[1L]]
            }
            if (soln$type == "Intersection") {
              if (!.sset_set_is(soln$args[[1L]], "Reals")) st$intersections <- .sd_set(st$intersections, sym, soln$args[[1L]])
              soln_new <- Union(soln_new, soln$args[[2L]])
            }
            if (!.sset_is_empty(soln_new)) soln <- soln_new
            if (index > 1L && real && !(soln$type %in% c("ImageSet", "ConditionSet")))
              soln <- Union(soln, solveset_complex(eq2, sym))
            if (!(soln$type %in% c("FiniteSet", "ImageSet", "ConditionSet", "Union", "EmptySet")))
              .sset_err("NotImplementedError", sprintf("nonlinsolve cannot handle solution of type %s for symbol %s. Got %s from equation %s.",
                                                       soln$type, .str(sym), .str(soln), .str(eq2)))
            soln
          }, ValueError = function(e) NULL)
          if (is.null(soln)) next
          if (soln$type == "ConditionSet") {
            if (.sset_set_is(soln$args[[3L]], "Reals") || .sset_set_is(soln$args[[3L]], "Complexes")) {
              soln <- S_EmptySet
              not_solvable <- TRUE
              total_conditionst <- total_conditionst + 1L
            } else soln <- soln$args[[3L]]
          }
          if (!.sset_is_empty(soln)) {
            r <- extract_main_soln(sym, soln, soln_imageset)
            soln <- r[[1L]]; soln_imageset <- r[[2L]]
          }
          for (sol in .sset_elems(soln)) {
            r <- extract_main_soln(sym, sol, soln_imageset)
            sol <- r[[1L]]; soln_imageset <- r[[2L]]
            if (.sset_is_empty(sol)) next
            sol <- .sset_elems(sol)[[1L]]
            if (length(got_symbol) && any(vapply(got_symbol, function(g) .has1(sol, g), NA))) next
            rnew <- res
            if (!is_set(sol)) for (i in seq_along(res$k)) {
              v <- res$v[[i]]
              if (!is_set(v)) rnew$v[[i]] <- subs(v, sym, sol)
            }
            img <- .sd_get(soln_imageset, sol)
            if (!is.null(img)) {
              lam <- img$args[[1L]]
              rnew <- .sd_set(rnew, sym, do.call(apply_lambda, c(list(lam), rep(list(S_Zero), length(lam$args[[1L]]$args)))))
            } else rnew <- .sd_set(rnew, sym, sol)
            r <- append_new_soln(rnew, sym, sol, imgset_yes, soln_imageset, original_imageset, newresult)
            newresult <- r[[1L]]
            if (r[[2L]]) {
              result <- remove_res(result, res)
              if (aliased) input <- result
            }
          }
          if (!not_solvable) got_symbol[[length(got_symbol) + 1L]] <- sym
        }
      }
      if (length(newresult)) { result <- newresult; aliased <- FALSE }
    }
    list(result, total_solvest_call, total_conditionst, input)
  }

  rr <- solve_using_known_values(old_result, TRUE)
  new_result_real <- rr[[1L]]
  rc <- solve_using_known_values(rr[[4L]], FALSE)
  new_result_complex <- rc[[1L]]
  total_conditionset <- total_conditionset + rr[[3L]] + rc[[3L]]
  total_solveset_call <- total_solveset_call + rr[[2L]] + rc[[2L]]
  if (total_conditionset == total_solveset_call && total_solveset_call != -1L)
    return(.sset_return_conditionset(eqs_in_better_order, all_symbols, exclude))

  filtered_complex <- list()
  for (i in new_result_complex) {
    dup <- FALSE
    for (j in new_result_real) {
      if (!setequal(vapply(i$k, .k, ""), vapply(j$k, .k, "")) || length(i$k) != length(j$k)) next
      if (all(mapply(.sset_dummy_eq_any, i$v, j$v))) { dup <- TRUE; break }
    }
    if (!dup) filtered_complex[[length(filtered_complex) + 1L]] <- i
  }
  result <- c(new_result_real, filtered_complex)

  result_all_variables <- list()
  result_infinite <- list()
  for (res in result) {
    if (!length(res$k)) next
    if (length(res$k) < length(all_symbols)) {
      for (s in all_symbols) if (is.null(.sd_get(res, s))) res <- .sd_set(res, s, s)
      result_infinite[[length(result_infinite) + 1L]] <- res
    }
    if (!.sset_dict_in(res, result_all_variables)) result_all_variables[[length(result_all_variables) + 1L]] <- res
  }
  if (length(result_infinite)) result_all_variables <- result_infinite
  processed <- list()
  for (res in result_all_variables) {
    pr <- .sset_process_exclusions_for_branch(res, exclude, all_symbols)
    if (is.null(pr)) next
    branch_complements <- pr[[1L]]; conditions <- pr[[2L]]
    for (i in seq_along(st$complements$k)) {
      symbol <- st$complements$k[[i]]
      old <- .sd_get(branch_complements, symbol)
      branch_complements <- .sd_set(branch_complements, symbol, Union(if (is.null(old)) S_EmptySet else old, st$complements$v[[i]]))
    }
    restricted <- add_intersection_complement(list(res), st$intersections, branch_complements)
    for (branch in restricted) processed[[length(processed) + 1L]] <- list(res, branch_complements, conditions, branch)
  }
  if (any(vapply(processed, function(p) length(p[[3L]]) > 0L, NA))) {
    return(.sset_union(lapply(processed, function(p)
      .sset_branch_as_conditionset(p[[1L]], all_symbols, st$intersections, p[[2L]], p[[3L]]))))
  }
  FiniteSet(lapply(processed, function(p) do.call(Tuple, lapply(all_symbols, function(s) {
    v <- .sd_get(p[[4L]], s)
    if (is.null(v)) .sset_err("KeyError", .str(s))
    v
  }))))
}

# ---------------------------------------------------------------------------
# nonlinsolve

.sset_solveset_work <- function(system, symbols) {
  soln <- solveset(system[[1L]], symbols[[1L]])
  if (.sset_is_empty(soln)) return(S_EmptySet)
  if (.sset_is_fs(soln)) return(FiniteSet(lapply(soln$args, function(s) Tuple(s))))
  FiniteSet(Tuple(soln))
}

.sset_handle_positive_dimensional <- function(polys, symbols, denominators) {
  syms <- .sp_sorted(symbols)
  basis <- .sp_groebner(lapply(polys, .sp_poly_expr), syms)
  .sset_substitution(basis, symbols, list(.sd_new()), list(), denominators, symbols)
}

# leading monomials in lex order cover every generator by a pure power
.sset_zero_dimensional <- function(basis, symbols) {
  n <- length(symbols)
  ex <- integer(n)
  for (b in basis) {
    p <- .sp_as_poly(b, symbols)
    if (is.null(p) || !length(p$terms)) next
    mons <- lapply(p$terms, `[[`, "monom")
    lm <- mons[[1L]]
    for (m in mons[-1L]) {
      d <- m - lm
      nz <- which(d != 0L)
      if (length(nz) && d[nz[1L]] > 0L) lm <- m
    }
    if (sum(lm > 0L) == 1L) ex <- ex + lm
  }
  all(ex > 0L)
}

.sset_separate_poly_nonpoly <- function(system, symbols) {
  polys <- list(); polys_expr <- list(); nonpolys <- list(); unrad_changed <- list()
  denominators <- list()
  poly <- NULL
  for (eq in system) {
    for (d in .sp_simple_dens(eq, symbols))
      if (!any(vapply(denominators, function(z) .sset_eq(z, d), NA))) denominators[[length(denominators) + 1L]] <- d
    if (eq$type == "Equality") eq <- .sub(eq$args[[1L]], eq$args[[2L]])
    wr <- .sp_unrad(simplify(eq), symbols)
    if (!is.null(wr)) {
      unrad_changed[[length(unrad_changed) + 1L]] <- eq
      if (!length(wr[[2L]])) eq <- wr[[1L]]
    }
    if (.sset_is_expr(eq)) {
      eq <- .as_numer_denom(eq)[[1L]]
      poly <- .sp_as_poly(eq, symbols)
    } else if (is_zero_num(eq)) next
    if (!is.null(poly)) {
      polys[[length(polys) + 1L]] <- poly
      polys_expr[[length(polys_expr) + 1L]] <- .sp_poly_expr(poly)
    } else nonpolys[[length(nonpolys) + 1L]] <- eq
  }
  list(polys = polys, polys_expr = polys_expr, nonpolys = nonpolys, denominators = denominators, unrad_changed = unrad_changed)
}

.sset_has_float <- function(e) {
  if (e$type == "Float") return(TRUE)
  for (a in e$args) if (.sset_has_float(a)) return(TRUE)
  FALSE
}

# _handle_poly: list(poly_sol, poly_eqs)
.sset_handle_poly <- function(polys, symbols) {
  exprs <- lapply(polys, .sp_poly_expr)
  inexact <- any(vapply(exprs, .sset_has_float, NA))
  if (inexact) exprs <- lapply(exprs, function(e) nsimplify(N(e), rational = TRUE))
  # ponytail: the lex basis directly; SymPy's grevlex basis gives the same
  # inconsistency and dimension tests and fglm gives this lex basis
  basis <- .sp_groebner(exprs, symbols)
  if (any(vapply(basis, function(b) is_exact_number(b) && !is_zero_num(b), NA)))
    return(list(list(), list()))
  if (.sset_zero_dimensional(basis, symbols)) {
    if (inexact) basis <- lapply(basis, nfloat)
    result <- tryCatch(.sp_solve_poly_system(lapply(basis, function(b) .sp_as_poly(b, symbols)), symbols, strict = TRUE),
                       error = function(e) if (grepl("UnsolvableFactorError", conditionMessage(e)) || inherits(e, "UnsolvableFactorError")) "unsolvable" else stop(e))
    if (identical(result, "unsolvable")) return(list(list(.sd_new()), basis))
    poly_sol <- lapply(result, function(res) {
      d <- .sd_new()
      for (i in seq_along(symbols)) d <- .sd_set(d, symbols[[i]], res[[i]])
      d
    })
    return(list(poly_sol, list()))
  }
  if (inexact) basis <- lapply(basis, nfloat)
  list(list(.sd_new()), basis)
}

#' Solve a system of nonlinear equations as a set
#'
#' `nonlinsolve(system, ...)` solves a system of equations for the given
#' symbols and returns a `FiniteSet` of `Tuple`s, as SymPy's `nonlinsolve()`:
#' polynomial systems go through a Groebner basis (positive-dimensional ones
#' leave free unknowns as parameters), other equations are solved by
#' [substitution()], with complex and periodic (`ImageSet`) solutions, and a
#' `ConditionSet` when nothing can be solved.
#'
#' @param system A list of equations or expressions.
#' @param ... The symbols, or one list of them.
#' @return A set of tuples, or a `ConditionSet`.
#' @examples
#' x <- Symbol("x"); y <- Symbol("y")
#' nonlinsolve(list(x^2 + y^2 - 1, x - y), list(x, y))
#' nonlinsolve(list(x*y - 1, x - 2), list(x, y))
#' @export
nonlinsolve <- function(system, ...) {
  system <- .sset_seq(system)
  if (is_sym(system)) system <- list(system)
  if (!length(system)) return(S_EmptySet)
  symbols <- list(...)
  if (!length(symbols))
    .sset_err("ValueError", .filldedent("Symbols must be given, for which solution of the system is to be found."))
  if (is.list(symbols[[1L]]) && !is_sym(symbols[[1L]])) symbols <- symbols[[1L]]
  else if (is_sym(symbols[[1L]]) && symbols[[1L]]$type == "Tuple") symbols <- symbols[[1L]]$args
  if (!length(symbols))
    .sset_err("IndexError", .filldedent("Symbols must be given, for which solution of the system is to be found."))
  symbols <- lapply(symbols, as_sym)
  system <- lapply(system, as_sym)
  rc <- .sp_recast(system, symbols)
  if (length(rc$restore)) {
    soln <- nonlinsolve(rc$f, rc$symbols)
    return(subs(soln, rc$restore, simultaneous = TRUE))
  }
  if (length(system) == 1L && length(symbols) == 1L) return(.sset_solveset_work(system, symbols))
  sp <- .sset_separate_poly_nonpoly(system, symbols)
  poly_eqs <- list()
  poly_sol <- list(.sd_new())
  if (length(sp$polys)) {
    hp <- .sset_handle_poly(sp$polys, symbols)
    poly_sol <- hp[[1L]]; poly_eqs <- hp[[2L]]
    if (length(poly_sol) && length(poly_sol[[1L]]$k)) {
      fsk <- function(lst) unique(unlist(lapply(lst, function(e) vapply(free_symbols(e), .k, ""))))
      poly_syms <- fsk(sp$polys_expr)
      unrad_syms <- fsk(sp$unrad_changed)
      if (length(sp$unrad_changed) && setequal(unrad_syms, poly_syms))
        poly_sol <- Filter(function(sol) isTRUE(.sp_checksol(sp$unrad_changed, sol)), poly_sol)
    }
  }
  remaining <- c(poly_eqs, sp$nonpolys)
  to_tuple <- function(sol) do.call(Tuple, lapply(symbols, function(s) { v <- .sd_get(sol, s); if (is.null(v)) s else v }))
  if (!length(remaining)) {
    if (length(sp$denominators)) {
      branches <- lapply(poly_sol, function(sol) {
        d <- .sd_new()
        for (s in symbols) { v <- .sd_get(sol, s); d <- .sd_set(d, s, if (is.null(v)) s else v) }
        d
      })
      return(.sset_substitution(list(S_Zero), symbols, branches, list(), sp$denominators, symbols))
    }
    return(FiniteSet(lapply(poly_sol, to_tuple)))
  }
  subs_res <- .sset_substitution(remaining, symbols, poly_sol, list(), sp$denominators, symbols)
  if (!.sset_is_fs(subs_res)) return(subs_res)
  if (length(sp$unrad_changed)) {
    result <- lapply(subs_res$args, function(sol) {
      d <- .sd_new()
      for (i in seq_along(symbols)) d <- .sd_set(d, symbols[[i]], sol$args[[i]])
      d
    })
    correct <- Filter(function(sol) !identical(.sp_checksol(sp$unrad_changed, sol), FALSE), result)
    return(FiniteSet(lapply(correct, to_tuple)))
  }
  subs_res
}

# The named-list solution of a linear system (NULL when inconsistent), for
# internal callers that read the values by symbol name.
.linsolve_named <- function(eqs, syms) {
  eqs <- lapply(eqs, function(e) expand(.eq_to_expr(e)))
  syms <- .as_gens(list(syms))
  if (!all(vapply(eqs, .is_linear_in, TRUE, syms = syms))) stop("system is not linear")
  .solve_linear_system(eqs, syms)
}
