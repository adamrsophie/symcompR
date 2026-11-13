# Arbitrary constants in ODE solutions, after sympy/solvers/ode/ode.py:
# constantsimp, constant_renumber and their helpers.

# Basic.count(query): the number of subtrees equal to query
.cs_count <- function(e, q) {
  n <- 0L
  walk <- function(z) {
    if (identical(z$key, q$key)) n <<- n + 1L
    for (a in z$args) walk(a)
  }
  walk(e)
  n
}

.cs_subset <- function(fs, Cs) all(vapply(fs, function(s) .sp_in(s, Cs), TRUE))

.cs_get_constant_subexpressions <- function(expr, Cs) {
  Ces <- list()
  walk <- function(e) {
    fs <- free_symbols(e)
    if (length(fs) && .cs_subset(fs, Cs)) {
      Ces[[length(Ces) + 1L]] <<- e
      return(invisible())
    }
    if (e$type == "Function" && e$name == "exp") e <- expand_mul(e)
    if (e$type %in% c("Add", "Mul")) {
      tr <- Filter(function(i) .cs_subset(free_symbols(i), Cs), e$args)
      if (length(tr) > 1L) {
        x <- if (e$type == "Add") .add(tr) else .mul(tr)
        if (!.sp_is_number(x)) Ces[[length(Ces) + 1L]] <<- x
      }
    } else if (e$type == "Integral") {
      if (.cs_subset(free_symbols(e), Cs) && all(vapply(.integral_limits(e), function(l) length(l) == 3L, TRUE)))
        Ces[[length(Ces) + 1L]] <<- e
    }
    for (a in e$args) walk(a)
  }
  walk(expr)
  Ces
}

.integral_limits <- function(e) lapply(e$args[-1L], function(l) if (l$type == "Tuple") l$args else list(l))

.cs_sort_str <- function(xs) xs[order(vapply(xs, .str, ""), method = "radix")]

.cs_remove_linear_redundancies <- function(expr, Cs) {
  cnts <- vapply(Cs, function(i) .cs_count(expr, i), 0L)
  Cs <- Cs[cnts > 0L]
  cnts <- cnts[cnts > 0L]
  linear <- function(e) {
    if (e$type == "Add") {
      xs <- list()
      for (j in seq_along(Cs)) {
        i <- Cs[[j]]
        if (.cs_count(e, i) == cnts[j] && is_zero_num(.diff1(.diff1(e, i), i))) xs[[length(xs) + 1L]] <- i
      }
      dk <- character(); dv <- list()
      for (x in xs) {
        y <- .diff1(e, x)
        k <- y$key
        if (!(k %in% dk)) { dk <- c(dk, k); dv[[k]] <- list() }
        dv[[k]] <- c(dv[[k]], list(x))
      }
      for (k in dk) {
        if (length(dv[[k]]) > 1L) {
          srt <- .cs_sort_str(dv[[k]])
          for (x in srt[-1L]) e <- subs(e, x, S_Zero)
        }
      }
    }
    e
  }
  walk <- function(e) {
    if (length(e$args)) e <- sym_rebuild(e, lapply(e$args, walk))
    linear(e)
  }
  if (e_is_eq <- (expr$type == "Equality")) {
    lhs <- walk(expr$args[[1L]]); rhs <- walk(expr$args[[2L]])
    f <- function(i) is_number_type(i) || .sp_in(i, Cs)
    if (lhs$type == "Symbol" && .sp_in(lhs, Cs)) { tmp <- lhs; lhs <- rhs; rhs <- tmp }
    atomic <- function(z) !length(z$args)
    if (lhs$type %in% c("Add", "Symbol") && rhs$type %in% c("Add", "Symbol")) {
      la <- if (atomic(lhs)) list(lhs) else lhs$args
      ra <- if (atomic(rhs)) list(rhs) else rhs$args
      lt <- Filter(f, la); lf <- Filter(Negate(f), la)
      rt <- Filter(f, ra); rf <- Filter(Negate(f), ra)
      if (!length(lt)) lt <- list(S_Zero)
      if (!length(lf)) lf <- list(S_Zero)
      if (!length(rt)) rt <- list(S_Zero)
      if (!length(rf)) rf <- list(S_Zero)
      lhs <- .sub(.add(lf), .add(rf))
      rhs <- .sub(.add(rt), .add(lt))
    } else if (lhs$type %in% c("Mul", "Symbol") && rhs$type %in% c("Mul", "Symbol")) {
      la <- if (atomic(lhs)) list(lhs) else lhs$args
      lt <- Filter(f, la); lf <- Filter(Negate(f), la)
      if (length(lt)) {
        if (!length(lf)) lf <- list(S_One)
        lhs <- .mul(lf)
        rhs <- .div(rhs, .mul(lt))
      }
    }
    return(Eq(lhs, rhs))
  }
  walk(expr)
}

# Poly(f, expand=False) over the generators found in f, as terms_gcd uses it
.cs_terms_gcd <- function(f, clear = TRUE, deep = FALSE) {
  orig <- f
  if (f$type == "Equality") return(Eq(.cs_terms_gcd(f$args[[1L]], clear, deep), .cs_terms_gcd(f$args[[2L]], clear, deep)))
  if (!.is_expr_node(f) || !length(f$args)) return(orig)
  if (deep) {
    new <- sym_rebuild(f, lapply(f$args, .cs_terms_gcd, clear = clear, deep = TRUE))
    return(.cs_terms_gcd(new, clear = clear, deep = FALSE))
  }
  gens <- .find_gens(list(f))
  if (!length(gens)) return(orig)
  p <- .sp_as_poly_noexpand(f, gens)
  if (is.null(p) || !length(p$terms)) return(orig)
  if (!all(vapply(p$terms, function(t) is_exact_number(t$coeff), TRUE))) return(orig)
  J <- Reduce(pmin, lapply(p$terms, `[[`, "monom"))
  vals <- lapply(p$terms, function(t) gmp::as.bigq(t$coeff$value))
  den <- Reduce(function(a, b) gmp::lcm.bigz(a, b), lapply(vals, gmp::denominator))
  nums <- lapply(vals, function(v) gmp::numerator(v * den))
  cont <- abs(Reduce(function(a, b) gmp::gcd.bigz(a, b), nums))
  coeff <- make_rational(gmp::as.bigq(cont, den))
  rest <- .add(lapply(seq_along(p$terms), function(i) {
    t <- p$terms[[i]]
    m <- t$monom - J
    .mul(c(list(make_rational(gmp::as.bigq(nums[[i]], cont))),
           Map(function(g, k) .pow(g, .int(k)), gens[m > 0L], m[m > 0L])))
  }))
  term <- .mul(Map(function(g, k) .pow(g, .int(k)), gens[J > 0L], J[J > 0L]))
  if (is_one_num(coeff)) {
    coeff <- S_One
    if (is_one_num(term)) return(orig)
  }
  if (clear) return(.keep_coeff(coeff, .mul(list(term, rest))))
  kc <- as_coeff_mul_split(.keep_coeff(coeff, rest, clear = FALSE))
  .keep_coeff(kc[[1L]], .mul(list(term, kc[[2L]])), clear = FALSE)
}

# polynomial terms of f in gens without expanding f first
.sp_as_poly_noexpand <- function(f, gens) {
  n <- length(gens)
  gk <- vapply(gens, .k, "")
  terms <- list(); keys <- character()
  for (t in .aargs(f)) {
    coeff <- list(); mon <- integer(n)
    for (fac in .margs(t)) {
      if (is_number_type(fac)) { coeff[[length(coeff) + 1L]] <- fac; next }
      dp <- .decompose_power(fac)
      base <- dp[[1L]]; ex <- dp[[2L]]
      if (ex < 0) { ex <- -ex; base <- .pow(base, S_NegativeOne) }
      j <- match(base$key, gk)
      if (is.na(j)) return(NULL)
      mon[j] <- mon[j] + as.integer(ex)
    }
    k <- paste(mon, collapse = ",")
    i <- match(k, keys)
    cf <- .mul(coeff)
    if (is.na(i)) {
      keys <- c(keys, k)
      terms[[length(terms) + 1L]] <- list(monom = mon, coeff = cf)
    } else terms[[i]]$coeff <- .add(list(terms[[i]]$coeff, cf))
  }
  list(terms = Filter(function(t) !is_zero_num(t$coeff), terms), gens = gens)
}

#' Simplify arbitrary constants in an ODE solution
#'
#' `constantsimp(expr, constants)` absorbs subexpressions made only of the
#' arbitrary constants into those constants and removes redundant ones, as
#' SymPy's `constantsimp()` does for the solutions returned by [dsolve()].
#' `constant_renumber(expr)` renumbers constants named `C1`, `C2`, ... in
#' printing order.
#'
#' @param expr An expression or equation.
#' @param constants A list of the constant symbols.
#' @param variables Symbols that are not constants; by default the symbols
#'   named `C` followed by digits are the constants.
#' @param newconstants Symbols to use as the new constants.
#' @return An expression or equation.
#' @examples
#' x <- Symbol("x"); C1 <- Symbol("C1"); C2 <- Symbol("C2")
#' constantsimp(2 * C1 * x, list(C1))
#' constantsimp(C1 + C2 + x, list(C1, C2))
#' constant_renumber(C2 * x + C1)
#' @export
constantsimp <- function(expr, constants) {
  expr <- as_sym(expr)
  Cs <- lapply(constants, as_sym)
  orig <- expr
  for (xe in .cs_get_constant_subexpressions(expr, Cs)) {
    xes <- free_symbols(xe)
    if (!length(xes)) next
    if (all(vapply(xes, function(c) .cs_count(expr, c) == .cs_count(xe, c), TRUE))) {
      xes <- .cs_sort_str(xes)
      expr <- subs(expr, xe, xes[[1L]])
    }
  }
  cs <- tryCatch(cse(expr), error = function(e) NULL)
  if (!is.null(cs)) {
    commons <- rev(cs$replacements)
    rexpr <- cs$reduced[[1L]]
    for (s in commons) {
      csyms <- Filter(function(z) z$type == "Symbol", atoms(s[[2L]]))
      ratoms <- Filter(function(z) z$type == "Symbol", atoms(rexpr))
      others <- Filter(function(ex) !identical(ex[[1L]]$key, s[[1L]]$key), commons)
      if (length(csyms) == 1L && .sp_in(csyms[[1L]], Cs) && !.sp_in(csyms[[1L]], ratoms) &&
          !any(vapply(others, function(ex) .has1(ex[[2L]], csyms[[1L]]) || identical(ex[[1L]]$key, csyms[[1L]]$key), TRUE))) {
        rexpr <- subs(rexpr, s[[1L]], csyms[[1L]])
      } else {
        rexpr <- subs(rexpr, s[[1L]], s[[2L]])
      }
    }
    expr <- rexpr
  }
  expr <- .cs_remove_linear_redundancies(expr, Cs)
  new_expr <- .cs_terms_gcd(expr, clear = FALSE, deep = TRUE)
  if (new_expr$type == "Mul") {
    infac <- FALSE
    asfac <- FALSE
    for (m in new_expr$args) {
      if (m$type == "Function" && m$name == "exp") {
        asfac <- TRUE
      } else if (m$type == "Add") {
        infac <- any(vapply(m$args, function(t) any(vapply(.margs(t), function(fi) fi$type == "Function" && fi$name == "exp", TRUE)), TRUE))
      }
      if (asfac && infac) {
        new_expr <- expr
        break
      }
    }
  }
  expr <- new_expr
  if (!identical(orig$key, expr$key)) return(constantsimp(expr, Cs))
  expr
}

#' @rdname constantsimp
#' @export
constant_renumber <- function(expr, variables = NULL, newconstants = NULL) {
  if (is.list(expr) && !is_sym(expr)) {
    r <- constant_renumber(new_sym("Tuple", lapply(expr, as_sym)), variables, newconstants)
    return(r$args)
  }
  expr <- as_sym(expr)
  if (!is.null(variables)) {
    variables <- lapply(variables, as_sym)
    constantsymbols <- Filter(function(s) !.sp_in(s, variables), free_symbols(expr))
  } else {
    variables <- list()
    constantsymbols <- Filter(function(s) grepl("^C[0-9]+$", s$name), free_symbols(expr))
  }
  found <- list()
  C_1 <- lapply(constantsymbols, function(ci) list(ci, S_One))
  skey <- function(arg) .sort_key(if (length(C_1)) subs(arg, C_1) else arg)
  rec <- function(e) {
    if (e$type == "Tuple") return(new_sym("Tuple", lapply(e$args, rec)))
    if (e$type == "Equality") return(Eq(rec(e$args[[1L]]), rec(e$args[[2L]])))
    if (!(e$type %in% c("Mul", "Add", "Pow")) && e$type != "Function" && !.has_free(e, constantsymbols)) return(e)
    if (e$type == "Piecewise") return(e)
    if (.sp_in(e, constantsymbols)) {
      if (!.sp_in(e, found)) found[[length(found) + 1L]] <<- e
      return(e)
    }
    if (e$type %in% c("Function", "Pow")) return(sym_rebuild(e, lapply(e$args, rec)))
    args <- e$args
    if (length(args) > 1L) args <- .key_sort(lapply(args, skey), args)
    sym_rebuild(e, lapply(args, rec))
  }
  expr <- rec(expr)
  found <- Filter(function(c) !.sp_in(c, variables), found)
  if (!length(found)) return(expr)
  if (is.null(newconstants)) {
    vn <- vapply(variables, function(v) if (is.null(v$name)) "" else v$name, "")
    newc <- list()
    i <- 1L
    while (length(newc) < length(found)) {
      nm <- paste0("C", i)
      if (!(nm %in% vn)) newc[[length(newc) + 1L]] <- Symbol(nm)
      i <- i + 1L
    }
  } else {
    newc <- Filter(function(s) !.sp_in(s, variables), lapply(newconstants, as_sym))[seq_along(found)]
  }
  subs(expr, Map(list, found, newc), simultaneous = TRUE)
}
