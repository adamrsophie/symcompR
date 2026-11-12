# Further solvers after sympy/solvers: decompogen.py, ode/ode.py
# (homogeneous_order), solvers.py (solve_linear_system_LU), polysys.py
# (solve_triangulated, factor_system), inequalities.py (rational and
# absolute value inequalities) and solveset.py (substitution).

# Dense nested coefficient list of e in gens, as Poly.rep.to_list().
.dense_rep <- function(e, gens) {
  if (!length(gens)) {
    v <- num_value_dbl(e)
    return(if (is.na(v)) .sort_key(e) else v)
  }
  e <- expand(e)
  if (is_zero_num(e)) return(list())
  g <- gens[[1L]]
  d <- degree(e, g)
  lapply(d:0, function(k) {
    ck <- coeff(e, g, k)
    if (length(gens) == 1L) .dense_rep(ck, list()) else .dense_rep(ck, gens[-1L])
  })
}

.poly_sort_key <- function(e, gens) {
  list(lapply(gens, function(g) if (.has1(e, g)) degree(e, g) else 0), .dense_rep(e, gens))
}

# Systems sorted as SymPy _sort_systems: each system by descending poly key,
# then the systems by descending zipped keys.
.sort_systems <- function(systems, gens) {
  systems <- lapply(systems, function(s) rev(.key_sort(lapply(s, .poly_sort_key, gens), s)))
  keys <- lapply(systems, function(s) {
    k <- lapply(s, .poly_sort_key, gens)
    list(lapply(k, `[[`, 1L), lapply(k, `[[`, 2L))
  })
  rev(.key_sort(keys, systems))
}

# Minimal sets of factors, one from each list, as SymPy _factor_sets.
.factor_sets <- function(sets) {
  if (!length(sets)) return(list(character()))
  res <- list()
  rec <- function(remaining, cur) {
    if (!length(remaining)) {
      res[[length(res) + 1L]] <<- sort(cur)
      return(invisible())
    }
    i <- which.min(lengths(remaining))
    nxt <- remaining[[i]]
    rest <- remaining[-i]
    for (f in nxt) rec(Filter(function(s) !(f %in% s), rest), c(cur, f))
  }
  rec(sets, character())
  res <- res[!duplicated(vapply(res, paste, "", collapse = "\r"))]
  Filter(function(s1) !any(vapply(res, function(s2) length(s2) < length(s1) && all(s2 %in% s1), TRUE)), res)
}

#' More equation solvers
#'
#' `decompogen(f, x)` splits `f` into a list of functions whose composition
#' is `f`. `homogeneous_order(eq, ...)` is the order of homogeneity of `eq`
#' in the given symbols or function applications, or `NULL`.
#' `solve_linear_system_LU(M, syms)` solves an augmented matrix with LU
#' decomposition. `solve_triangulated(polys, gens)` solves a zero
#' dimensional polynomial system through a lex Groebner basis and returns
#' its rational solutions as lists ordered like `gens`. `factor_system(eqs,
#' gens)` splits a polynomial system into irreducible subsystems: `list()`
#' means no solution and `list(list())` means every value solves it.
#'
#' @param f,eq An expression.
#' @param x A symbol.
#' @param ... Symbols or function applications.
#' @param M An augmented matrix.
#' @param syms A list of symbols.
#' @param polys,eqs Lists of polynomial expressions, each equal to zero.
#' @param gens Generators; all free symbols by default for `factor_system()`.
#' @return A list, a number or `NULL`.
#' @examples
#' x <- Symbol("x"); y <- Symbol("y")
#' decompogen(sin(x)^2 + sin(x) + 1, x)
#' homogeneous_order(x^2 + x*y, x, y)
#' solve_triangulated(list(x^2 + y^2 - 1, x - y + 1), list(x, y))
#' factor_system(list(x^2 - 1, y - 1))
#' @name more-solvers
NULL

#' @rdname more-solvers
#' @export
decompogen <- function(f, x) {
  f <- as_sym(f)
  x <- as_sym(x)
  if (is_relational(f)) stop("expecting an expression, got a relational", call. = FALSE)
  if (!.has1(f, x)) return(list(f))
  if (f$type == "Function" && f$name %in% c("Max", "Min")) {
    args <- f$args
    d0 <- NULL
    d <- list(x)
    for (i in seq_along(args)) {
      if (!.has1(args[[i]], x)) next
      d <- decompogen(args[[i]], x)
      if (length(d) == 1L) d <- c(list(x), d)
      tail_keys <- vapply(d[-1L], .k, "")
      if (is.null(d0)) d0 <- d[-1L]
      else if (!identical(tail_keys, vapply(d0, .k, ""))) {
        d <- list(x)
        break
      }
      args[[i]] <- d[[1L]]
    }
    if (identical(d[[1L]]$key, x$key)) return(list(f))
    return(c(list(fcall(f$name, args)), d0))
  }
  if (f$type %in% c("Function", "Pow")) {
    arg <- if (f$type == "Pow" && is_const(f$args[[1L]], "E")) f$args[[2L]] else f$args[[1L]]
    if (identical(arg$key, x$key)) return(list(f))
    return(c(list(subs(f, arg, x)), decompogen(arg, x)))
  }
  gens <- Filter(function(g) .has1(g, x), Poly(f)$gens)
  if (length(gens) == 1L && !identical(gens[[1L]]$key, x$key))
    return(c(list(subs(f, gens[[1L]], x)), decompogen(gens[[1L]], x)))
  if (length(free_symbols(f)) > 1L) return(list(f))
  tryCatch(sym_decompose(f, x), error = function(e) list(f))
}

#' @rdname more-solvers
#' @export
homogeneous_order <- function(eq, ...) {
  syms <- lapply(list(...), as_sym)
  if (!length(syms)) stop("homogeneous_order: no symbols were given.", call. = FALSE)
  eq <- as_sym(eq)
  if (.has_type(eq, c("Order", "Derivative"))) return(NULL)
  if (!length(free_symbols(eq))) return(S_Zero)
  keep <- list()
  for (s in syms) {
    if (s$type %in% c("Symbol", "Dummy")) {
      keep[[length(keep) + 1L]] <- s
      next
    }
    ak <- vapply(s$args, .k, "")
    if (!all(ak %in% vapply(syms, .k, ""))) return(NULL)
    d <- Dummy("d")
    eq <- subs(eq, s, d)
    keep[[length(keep) + 1L]] <- d
  }
  syms <- keep[!duplicated(vapply(keep, .k, ""))]
  if (!any(vapply(syms, function(s) .has1(eq, s), TRUE))) return(NULL)
  if (eq$type == "Function") {
    o <- do.call(homogeneous_order, c(list(eq$args[[1L]]), syms))
    return(if (!is.null(o) && is_zero_num(o)) S_Zero else NULL)
  }
  t <- Dummy("t", positive = TRUE)
  et <- xreplace(eq, lapply(syms, function(s) list(s, .mul(list(t, s)))))
  sv <- separatevars(et, list(t), dict = TRUE)
  if (is.null(sv)) return(NULL)
  tpart <- sv[[t$name]]
  if (identical(tpart$key, S_One$key)) return(S_Zero)
  facs <- if (tpart$type == "Mul") tpart$args else list(tpart)
  d <- .mul(Filter(function(fc) .has1(fc, t), facs))
  be <- .base_exp(d)
  if (identical(be[[1L]]$key, t$key)) be[[2L]] else NULL
}

#' @rdname more-solvers
#' @export
solve_linear_system_LU <- function(M, syms) {
  M <- Matrix(M)
  if (M$nrow != M$ncol - 1L) stop("Rows should be equal to columns - 1", call. = FALSE)
  syms <- .as_gens(list(syms))
  rows <- .mrows(M)
  n <- M$nrow
  A <- Matrix(lapply(rows, function(r) r[seq_len(n)]))
  b <- lapply(rows, function(r) r[[n + 1L]])
  lu <- LUdecomposition(A)
  L <- .mrows(lu$L)
  U <- .mrows(lu$U)
  pb <- b[lu$perm]
  y <- vector("list", n)
  for (i in seq_len(n)) {
    s <- pb[[i]]
    if (i > 1L) for (j in seq_len(i - 1L)) s <- .sub(s, .mul(list(L[[i]][[j]], y[[j]])))
    y[[i]] <- s
  }
  xs <- vector("list", n)
  for (i in rev(seq_len(n))) {
    if (is_zero_num(cancel(U[[i]][[i]]))) stop("Matrix is singular", call. = FALSE)
    s <- y[[i]]
    if (i < n) for (j in (i + 1L):n) s <- .sub(s, .mul(list(U[[i]][[j]], xs[[j]])))
    xs[[i]] <- cancel(.div(s, U[[i]][[i]]))
  }
  stats::setNames(xs, vapply(syms, function(s) s$name, ""))
}

#' @rdname more-solvers
#' @export
solve_triangulated <- function(polys, gens) {
  gens <- .as_gens(list(gens))
  n <- length(gens)
  G <- rev(as.list(do.call(groebner, c(list(polys), gens))))
  # roots in the ground domain: rational functions of the parameters with
  # rational coefficients, as SymPy's ground_roots()
  in_ground <- function(r) {
    if (is_exact_number(r)) return(TRUE)
    if (any(vapply(gens, function(g) .has1(r, g), TRUE))) return(FALSE)
    nd <- fraction(together(r))
    ps <- free_symbols(r)
    all(vapply(nd, function(z) !is.null(expr_to_qpoly(expand(z), ps, strict = TRUE)), TRUE))
  }
  rat_roots <- function(p, v) {
    rs <- Filter(in_ground, roots(p, v))
    rs[!duplicated(vapply(rs, .k, ""))]
  }
  sols <- lapply(rat_roots(G[[1L]], gens[[n]]), list)
  G <- G[-1L]
  if (n > 1L) for (k in (n - 1L):1L) {
    var <- gens[[k]]
    vars <- gens[(k + 1L):n]
    allowed <- vapply(c(list(var), vars), .k, "")
    new <- list()
    for (vals in sols) {
      H <- list()
      for (g in G) {
        fk <- vapply(Filter(function(s) any(vapply(gens, function(v) identical(v$key, s$key), TRUE)),
                            free_symbols(g)), .k, "")
        if (!all(fk %in% allowed) || !.has1(g, var)) next
        h <- expand(xreplace(g, Map(list, vars, vals)))
        if (degree(g, var) == degree(h, var)) H[[length(H) + 1L]] <- h
      }
      p <- H[[which.min(vapply(H, function(h) degree(h, var), 0))]]
      for (z in rat_roots(p, var)) new[[length(new) + 1L]] <- c(list(z), vals)
    }
    sols <- new
  }
  .key_sort(lapply(sols, function(s) lapply(s, .sort_key)), sols)
}

#' @rdname more-solvers
#' @export
factor_system <- function(eqs, gens = NULL) {
  if (!is.list(eqs) || inherits(eqs, "sym")) eqs <- list(eqs)
  eqs <- lapply(eqs, function(e) expand(.eq_to_expr(e)))
  fs <- unique_syms(lapply(eqs, free_symbols))
  gens <- if (is.null(gens)) fs else .as_gens(list(gens))
  if (!length(gens)) return(if (all(vapply(eqs, is_zero_num, TRUE))) list(list()) else list())
  params <- Filter(function(s) !any(vapply(gens, function(g) identical(g$key, s$key), TRUE)), fs)
  pool <- list()
  sets <- list()
  for (e in eqs) {
    if (is_zero_num(e)) next
    fl <- do.call(factor_list, c(list(e), gens, params))
    facs <- lapply(fl$factors, `[[`, 1L)
    dep <- vapply(facs, function(f) any(vapply(gens, function(g) .has1(f, g), TRUE)), TRUE)
    pf <- facs[dep]
    if (any(!dep)) pf <- c(pf, list(.mul(facs[!dep])))
    else if (!length(pf)) return(list())
    for (f in pf) pool[[f$key]] <- f
    sets[[length(sets) + 1L]] <- vapply(pf, .k, "")
  }
  if (!length(sets)) return(list(list()))
  out <- lapply(.factor_sets(sets), function(ks) unname(pool[ks]))
  generic <- Filter(function(s) all(vapply(s, function(f) any(vapply(gens, function(g) .has1(f, g), TRUE)), TRUE)), out)
  .sort_systems(generic, gens)
}

# ---------------------------------------------------------------------------
# Inequalities

# Solution set of one inequality, as a subset of the reals.
.ineq_set <- function(e, rel, x) {
  e <- as_sym(e)
  if (!.has1(e, x)) return(if (.is_true(Rel(e, 0, rel))) Reals else S_EmptySet)
  solve_univariate_inequality(Rel(e, 0, rel), x, relational = FALSE)
}

.rat_part <- function(p) as_sym(if (inherits(p, "Poly")) as_expr(p) else p)

#' Rational and absolute value inequalities
#'
#' `solve_rational_inequalities(eqs)` solves `eqs`, a list of lists whose
#' inner entries `list(list(n, d), rel)` mean `n/d rel 0`; inner lists are
#' joined by `And` and the outer list by `Or`. `reduce_abs_inequality(expr,
#' rel, x)` reduces `expr rel 0` with nested absolute values, and
#' `reduce_abs_inequalities(exprs, x)` a list of `list(expr, rel)` pairs.
#'
#' @param eqs Nested list of inequalities.
#' @param expr An expression.
#' @param rel A relation: `"<"`, `"<="`, `">"`, `">="`, `"=="` or `"!="`.
#' @param x The real variable.
#' @param exprs A list of `list(expr, rel)` pairs.
#' @return A set, or a relational condition on `x`.
#' @examples
#' x <- Symbol("x")
#' solve_rational_inequalities(list(list(list(list(x - 1, x + 1), "<"))))
#' reduce_abs_inequality(Abs(x - 5) - 3, "<", x)
#' @name abs-inequalities
NULL

#' @rdname abs-inequalities
#' @export
solve_rational_inequalities <- function(eqs) {
  result <- S_EmptySet
  for (grp in eqs) {
    if (!length(grp)) next
    sol <- Reals
    for (q in grp) {
      n <- .rat_part(q[[1L]][[1L]])
      d <- .rat_part(q[[1L]][[2L]])
      rel <- q[[2L]]
      x <- c(free_symbols(n), free_symbols(d))
      x <- if (length(x)) x[[1L]] else Symbol("x")
      s <- .ineq_set(expand(.mul(list(n, d))), rel, x)
      if (.has1(d, x)) s <- Complement(s, .ineq_set(d, "==", x))
      sol <- Intersection(sol, s)
      if (identical(sol$key, S_EmptySet$key)) break
    }
    result <- Union(result, sol)
  }
  result
}

# SymPy _bottom_up_scan: every sign case of the absolute values, as
# (expression, conditions) pairs.
.abs_cases <- function(e) {
  if (e$type %in% c("Add", "Mul")) {
    out <- NULL
    for (a in e$args) {
      cs <- .abs_cases(a)
      if (is.null(out)) {
        out <- cs
        next
      }
      nxt <- list()
      for (p in out) for (q in cs)
        nxt[[length(nxt) + 1L]] <- list(e = if (e$type == "Add") .add(list(p$e, q$e)) else .mul(list(p$e, q$e)),
                                       conds = c(p$conds, q$conds))
      out <- nxt
    }
    return(out)
  }
  if (e$type == "Pow") {
    n <- e$args[[2L]]
    if (n$type != "Integer") stop("Only Integer Powers are allowed on Abs.", call. = FALSE)
    return(lapply(.abs_cases(e$args[[1L]]), function(p) list(e = .pow(p$e, n), conds = p$conds)))
  }
  if (e$type == "Function" && e$name == "Abs") {
    out <- list()
    for (p in .abs_cases(e$args[[1L]])) {
      out[[length(out) + 1L]] <- list(e = p$e, conds = c(p$conds, list(list(p$e, ">="))))
      out[[length(out) + 1L]] <- list(e = .neg(p$e), conds = c(p$conds, list(list(p$e, "<"))))
    }
    return(out)
  }
  list(list(e = e, conds = list()))
}

#' @rdname abs-inequalities
#' @export
reduce_abs_inequality <- function(expr, rel, x) {
  x <- as_sym(x)
  if (identical(.facts(x)[["extended_real"]], FALSE))
    stop("Cannot solve inequalities with absolute values containing non-real variables.", call. = FALSE)
  flip <- c("<" = ">", "<=" = ">=")
  sol <- S_EmptySet
  for (cs in .abs_cases(as_sym(expr))) {
    main <- if (rel %in% names(flip)) list(.neg(cs$e), flip[[rel]]) else list(cs$e, rel)
    s <- Reals
    for (q in c(list(main), cs$conds)) s <- Intersection(s, .ineq_set(q[[1L]], q[[2L]], x))
    sol <- Union(sol, s)
  }
  as_relational(sol, x)
}

#' @rdname abs-inequalities
#' @export
reduce_abs_inequalities <- function(exprs, x) {
  do.call(And, lapply(exprs, function(q) reduce_abs_inequality(q[[1L]], q[[2L]], x)))
}

#' Solve a system by substitution
#'
#' `substitution(system, symbols, result, known_symbols, exclude,
#' all_symbols)` solves `system` for `symbols` after substituting each
#' partial solution in `result`, and returns a `FiniteSet` of `Tuple`s
#' ordered like `all_symbols` (default `symbols`). Unknowns left free stay
#' symbolic. This is SymPy's `substitution()` interface over [solve()].
#'
#' @param system A list of equations.
#' @param symbols Unknowns.
#' @param result A list of named lists of known values.
#' @param known_symbols,exclude Ignored; kept for SymPy compatibility.
#' @param all_symbols Symbols to report, in order.
#' @return A `FiniteSet`.
#' @examples
#' x <- Symbol("x"); y <- Symbol("y")
#' substitution(list(x + y), list(x), list(list(y = 1)), list(y), list(), list(x, y))
#' @export
substitution <- function(system, symbols, result = list(list()), known_symbols = list(),
                         exclude = list(), all_symbols = NULL) {
  if (!is.list(system) || inherits(system, "sym")) system <- list(system)
  symbols <- .as_gens(list(symbols))
  all_symbols <- if (is.null(all_symbols)) symbols else .as_gens(list(all_symbols))
  out <- list()
  for (r in result) {
    eqs <- lapply(system, function(e) if (length(r)) subs(.eq_to_expr(e), r) else .eq_to_expr(e))
    eqs <- Filter(function(e) !is_zero_num(e), eqs)
    sols <- if (length(eqs) == 1L && length(symbols) == 1L) {
      # one unknown: real solutions first, then the complex families, as
      # SymPy's substitution() does with solveset()
      parts <- function(S) if (S$type == "FiniteSet") S$args else if (S$type == "Union") do.call(c, lapply(S$args, parts))
                           else if (S$type == "Intersection" && any(vapply(S$args, function(a) a$type == "FiniteSet", TRUE)))
                             parts(Filter(function(a) a$type == "FiniteSet", S$args)[[1L]])
                           else if (S$type == "EmptySet") list() else list(S)
      els <- c(parts(solveset(eqs[[1L]], symbols[[1L]], Reals)), parts(solveset(eqs[[1L]], symbols[[1L]], Complexes)))
      # the two solveset() calls use different Dummy indices
      els <- els[!duplicated(vapply(els, .str, ""))]
      lapply(els, function(v) stats::setNames(list(v), symbols[[1L]]$name))
    } else solve(eqs, symbols, dict = TRUE)
    for (s in sols) {
      vals <- c(r, s)
      # exclusions (mostly denominators) must not vanish at a solution
      if (length(exclude) && !any(vapply(vals, function(v) is_sym(v) && is_set(v), TRUE)) &&
          any(vapply(exclude, function(ex) is_zero_num(simplify(subs(as_sym(ex), vals))), TRUE))) next
      out[[length(out) + 1L]] <- do.call(Tuple, lapply(all_symbols, function(v) {
        z <- vals[[v$name]]
        if (is.null(z)) return(v)
        z <- as_sym(z)
        # only free symbols are replaced (CRootOf keeps its variable)
        if (is_set(z)) {
          for (ex in exclude) {
            bad <- tryCatch(solveset(as_sym(ex), v), error = function(e) NULL)
            if (!is.null(bad) && bad$type == "FiniteSet") z <- Complement(z, bad)
          }
          return(z)
        }
        if (any(vapply(free_symbols(z), function(q) !is.null(vals[[q$name]]), TRUE))) subs(z, vals) else z
      }))
    }
  }
  if (!length(out)) return(S_EmptySet)
  FiniteSet(out)
}

# ---------------------------------------------------------------------------
# Linear recurrences with polynomial coefficients (sympy/solvers/recurr.py)

.rs_deg <- function(e, n) if (is_zero_num(e)) -Inf else degree(e, n)

# Nonnegative integer roots of p in x.
# Nonnegative integer roots of p in x, as bigz.
.nni_roots <- function(p, x) {
  p <- expand(p)
  if (is_zero_num(p) || !.has1(p, x)) return(gmp::as.bigz(integer()))
  rs <- Filter(function(r) r$type == "Integer" && r$value >= 0, roots(p, x))
  if (!length(rs)) return(gmp::as.bigz(integer()))
  gmp::as.bigz(vapply(rs, function(r) as.character(gmp::numerator(gmp::as.bigq(r$value))), ""))
}
.small_int <- function(z) {
  if (abs(z) > .Machine$integer.max) stop("recurrence degree bound is too large", call. = FALSE)
  as.integer(as.character(z))
}

# k-th forward difference of p at a, as SymPy's _delta.
.rs_delta <- function(p, k, n, a) {
  B <- S_One
  D <- subs(p, n, .add(list(a, .int(k))))
  if (k >= 1L) for (i in seq_len(k)) {
    B <- .mul(list(B, Rational(i - k - 1L, i)))
    D <- .add(list(D, .mul(list(B, subs(p, n, .add(list(a, .int(k - i))))))))
  }
  D
}

.rsolve_poly <- function(coeffs, f, n, shift = 0L) {
  f <- as_sym(f)
  if (!is_polynomial(f, n)) return(NULL)
  homogeneous <- is_zero_num(expand(f))
  coeffs <- lapply(coeffs, function(c) expand(as_sym(c)))
  r <- length(coeffs) - 1L
  polys <- lapply(0:r, function(i)
    expand(.add(lapply(i:r, function(j) .mul(list(coeffs[[j + 1L]], .int(choose(j, i))))))))
  degs <- vapply(polys, .rs_deg, 0, n = n)
  lcs <- Map(function(p, d) if (is.finite(d)) coeff(p, n, d) else S_Zero, polys, degs)
  d <- max(degs)
  b <- max(degs - 0:r)
  x <- Dummy("x")
  ff <- function(i) .mul(lapply(seq_len(i) - 1L, function(k) .sub(x, .int(k))))
  degree_poly <- .add(lapply(which(degs - 0:r == b), function(i) .mul(list(lcs[[i]], ff(i - 1L)))))
  nni <- .nni_roots(degree_poly, x)
  N <- if (length(nni)) .small_int(max(nni)) else integer()
  N <- if (homogeneous) c(N, -b - 1) else c(N, degree(f, n) - b, -b - 1)
  N <- max(N)
  if (N < 0) return(if (homogeneous) list(result = S_Zero, C = list()) else NULL)
  if (N <= r) {
    C <- lapply(0:N, function(i) Symbol(paste0("C", i + shift)))
    y <- .add(Map(function(c, i) .mul(list(c, .pow(n, .int(i)))), C, 0:N))
    E <- .add(lapply(0:r, function(i) .mul(list(coeffs[[i + 1L]], subs(y, n, .add(list(n, .int(i))))))))
    sol <- solve_undetermined_coeffs(.sub(E, f), C, n)
    if (is.null(sol)) return(NULL)
    C0 <- C
    C <- Filter(function(c) is.null(sol[[c$name]]), C)
    result <- if (length(sol)) xreplace(y, unname(Map(function(v, k) list(Symbol(k), v), sol, names(sol)))) else y
  } else {
    A <- r
    U <- N + A + b + 1L
    nr <- .nni_roots(polys[[r + 1L]], n)
    a <- if (length(nr)) make_rational(gmp::as.bigq(max(nr) + 1)) else S_Zero
    alpha <- list()
    for (i in (-A):d) {
      I <- rep(list(S_One), d + 1L)
      if (d >= 1L) for (k in seq_len(d))
        I[[k + 1L]] <- .div(.mul(list(I[[k]], .add(list(x, .int(i - k + 1L))))), .int(k))
      s <- list()
      for (j in 0:A) for (k in 0:d) {
        m <- i + j
        if (m < 0 || m > k) next
        s[[length(s) + 1L]] <- .mul(list(I[[k + 1L]], .int(choose(k, m)), .rs_delta(polys[[j + 1L]], k, n, a)))
      }
      alpha[[as.character(i)]] <- .add(s)
    }
    al <- function(i, v) subs(alpha[[as.character(i)]], x, .int(v))
    V <- lapply(seq_len(U) - 1L, function(i) lapply(seq_len(A) - 1L, function(j) if (i == j) S_One else S_Zero))
    G <- rep(list(S_Zero), U)
    for (i in A:(U - 1L)) {
      v <- rep(list(S_Zero), A)
      g <- S_Zero
      for (k in seq_len(A + b)) {
        if (i - k < 0) break
        B <- al(k - A, i - k)
        for (j in seq_len(A)) v[[j]] <- .add(list(v[[j]], .mul(list(B, V[[i - k + 1L]][[j]]))))
        if (!homogeneous) g <- .add(list(g, .mul(list(B, G[[i - k + 1L]]))))
      }
      den <- al(-A, i)
      for (j in seq_len(A)) V[[i + 1L]][[j]] <- .neg(.div(v[[j]], den))
      if (!homogeneous) G[[i + 1L]] <- .div(.sub(.rs_delta(f, i - A, n, a), g), den)
    }
    P <- rep(list(S_One), U)
    if (U > 1L) for (i in 2:U)
      P[[i]] <- expand(.div(.mul(list(P[[i - 1L]], .sub(.sub(n, a), .int(i - 2L)))), .int(i - 1L)))
    Q <- lapply(seq_len(A), function(j) .add(lapply(seq_len(U), function(i) expand(.mul(list(V[[i]][[j]], P[[i]]))))))
    h <- if (!homogeneous) .add(lapply(seq_len(U), function(i) expand(.mul(list(G[[i]], P[[i]]))))) else S_Zero
    C <- lapply(seq_len(A) - 1L, function(i) Symbol(paste0("C", i + shift)))
    gf <- function(i) .add(Map(function(c, q) .mul(list(c, .rs_delta(q, i, n, a))), C, Q))
    rng <- if (N + 1L <= U - 1L) (N + 1L):(U - 1L) else integer()
    E <- lapply(rng, function(i) expand(if (homogeneous) gf(i) else .add(list(gf(i), .rs_delta(h, i, n, a)))))
    sol <- list()
    if (length(E)) {
      sol <- .solve_linear_system(E, C)
      if (!length(sol)) return(if (homogeneous) list(result = S_Zero, C = list()) else NULL)
    }
    result <- h
    C0 <- C
    keep <- list()
    for (k in seq_along(C)) {
      c <- C[[k]]
      if (!is.null(sol[[c$name]])) s <- .mul(list(sol[[c$name]], Q[[k]]))
      else {
        s <- .mul(list(c, Q[[k]]))
        keep[[length(keep) + 1L]] <- c
      }
      result <- .add(list(result, expand(s)))
    }
    C <- keep
  }
  if (!identical(vapply(C, .k, ""), vapply(C0, .k, "")[seq_along(C)])) {
    result <- xreplace(result, Map(list, C, C0[seq_along(C)]))
    C <- C0[seq_along(C)]
  }
  list(result = result, C = C)
}

.rsolve_ratio <- function(coeffs, f, n) {
  f <- as_sym(f)
  if (!is_polynomial(f, n)) return(NULL)
  coeffs <- lapply(coeffs, as_sym)
  r <- length(coeffs) - 1L
  A <- expand(subs(coeffs[[r + 1L]], n, .sub(n, .int(r))))
  B <- coeffs[[1L]]
  h <- Dummy("h")
  res <- resultant(A, subs(B, n, .add(list(n, h))), n)
  if (!is_polynomial(res, h)) {
    nd <- fraction(together(res))
    res <- quo(nd[[1L]], nd[[2L]], h)
  }
  nni <- .nni_roots(res, h)
  if (!length(nni)) return(.rsolve_poly(coeffs, f, n))
  C <- S_One
  for (i in .small_int(max(nni)):0) {
    d <- gcd(A, subs(B, n, .add(list(n, .int(i)))), n)
    A <- quo(A, d, n)
    B <- quo(B, subs(d, n, .sub(n, .int(i))), n)
    C <- .mul(c(list(C), lapply(0:i, function(j) subs(d, n, .sub(n, .int(j))))))
  }
  denoms <- lapply(0:r, function(i) subs(C, n, .add(list(n, .int(i)))))
  numers <- vector("list", r + 1L)
  for (i in seq_len(r + 1L)) {
    g <- gcd(coeffs[[i]], denoms[[i]], n)
    numers[[i]] <- quo(coeffs[[i]], g, n)
    denoms[[i]] <- quo(denoms[[i]], g, n)
  }
  numers <- lapply(seq_len(r + 1L), function(i) .mul(c(list(numers[[i]]), denoms[-i])))
  res <- .rsolve_poly(numers, .mul(c(list(f), denoms)), n)
  if (is.null(res)) return(NULL)
  list(result = simplify(.div(res$result, C)), C = res$C)
}

.is_hypergeometric <- function(e, n) !is.null(hypersimp(e, n))

.rsolve_hyper <- function(coeffs, f, n) {
  coeffs <- lapply(coeffs, as_sym)
  f <- as_sym(f)
  r <- length(coeffs) - 1L
  kernel <- list()
  result <- S_Zero
  if (!is_zero_num(f)) {
    if (f$type == "Add") {
      reps <- list()
      sums <- list()
      fe <- expand(f)
      for (g in if (fe$type == "Add") fe$args else list(fe)) {
        if (!.is_hypergeometric(g, n)) return(NULL)
        hit <- 0L
        for (k in seq_along(reps)) if (hypersimilar(g, reps[[k]], n)) {
          hit <- k
          break
        }
        if (hit) sums[[hit]] <- .add(list(sums[[hit]], g))
        else {
          reps[[length(reps) + 1L]] <- g
          sums[[length(sums) + 1L]] <- S_Zero
        }
      }
      inhom <- Map(function(g, s) .add(list(g, s)), reps, sums)
    } else if (.is_hypergeometric(f, n)) {
      inhom <- list(f)
    } else return(NULL)
    for (i in seq_along(inhom)) {
      cf <- S_One
      polys <- coeffs
      denoms <- rep(list(S_One), r + 1L)
      s <- hypersimp(inhom[[i]], n)
      if (is.null(s)) return(NULL)
      for (j in seq_len(r)) {
        cf <- .mul(list(cf, subs(s, n, .add(list(n, .int(j - 1L))))))
        nd <- fraction(together(cf))
        polys[[j + 1L]] <- .mul(list(polys[[j + 1L]], nd[[1L]]))
        denoms[[j + 1L]] <- nd[[2L]]
      }
      polys <- lapply(seq_len(r + 1L), function(j) .mul(c(list(polys[[j]]), denoms[-j])))
      R <- .rsolve_ratio(polys, .mul(denoms), n)
      if (!is.null(R)) {
        Rv <- R$result
        if (length(R$C)) Rv <- xreplace(Rv, lapply(R$C, function(c) list(c, S_Zero)))
      } else {
        R2 <- .rsolve_poly(polys, .mul(denoms), n)
        Rv <- if (is.null(R2)) NULL else R2$result
      }
      if (is.null(Rv) || is_zero_num(Rv)) return(NULL)
      inhom[[i]] <- .mul(list(inhom[[i]], Rv))
      result <- simplify(.add(inhom))
    }
  }
  Z <- Dummy("Z")
  ur <- function(p, v = n) {
    rs <- roots(p, v)
    rs[!duplicated(vapply(rs, .k, ""))]
  }
  p_f <- ur(coeffs[[1L]])
  q_f <- ur(subs(coeffs[[r + 1L]], n, .add(list(n, .int(1L - r)))))
  mid <- list(list(S_One, S_One))
  for (p in p_f) for (q in q_f) {
    if (p$type == "Integer" && q$type == "Integer" && p$value <= q$value) next
    mid[[length(mid) + 1L]] <- list(.sub(n, p), .sub(n, q))
  }
  factors <- c(lapply(p_f, function(p) list(.sub(n, p), S_One)), mid,
               lapply(q_f, function(q) list(S_One, .sub(n, q))))
  for (AB in factors) {
    A <- AB[[1L]]
    B <- AB[[2L]]
    D <- .mul(list(A, subs(B, n, .add(list(n, .int(r - 1L))))))
    pl <- list()
    degs <- c()
    for (i in 0:r) {
      a <- .mul(lapply(seq_len(i) - 1L, function(j) subs(A, n, .add(list(n, .int(j))))))
      bb <- .mul(lapply(if (i < r) i:(r - 1L) else integer(), function(j) subs(B, n, .add(list(n, .int(j))))))
      pl[[i + 1L]] <- expand(quo(.mul(list(coeffs[[i + 1L]], a, bb)), D, n))
      if (!is_zero_num(pl[[i + 1L]])) degs <- c(degs, degree(pl[[i + 1L]], n))
    }
    if (!length(degs)) return(NULL)
    d <- max(degs)
    poly <- .add(lapply(0:r, function(i) .mul(list(coeff(pl[[i + 1L]], n, d), .pow(Z, .int(i))))))
    if (!.has1(poly, Z)) next
    for (z in ur(poly, Z)) {
      if (is_zero_num(z)) next
      rc <- lapply(0:r, function(i) .mul(list(pl[[i + 1L]], .pow(z, .int(i)))))
      if (d == 0 && !is_zero_num(simplify(.add(lapply(seq_len(r), function(j) .mul(list(rc[[j + 1L]], .int(j)))))))) {
        sol <- list(Symbol("C0"))
      } else {
        rp <- .rsolve_poly(rc, S_Zero, n, 0L)
        sol <- lapply(rp$C, function(s) coeff(expand(rp$result), s))
      }
      for (Cs in sol) {
        ratio <- simplify(.div(.mul(list(z, A, subs(Cs, n, .add(list(n, S_One))))), .mul(list(B, Cs))))
        n0 <- S_Zero
        for (nr in ur(fraction(together(ratio))[[2L]])) {
          if (.has1(nr, S_I)) return(NULL)
          if (.is_true(Lt(n0, .add(list(nr, S_One))))) n0 <- .add(list(nr, S_One))
        }
        K <- product(ratio, list(n, n0, .sub(n, S_One)))
        if (.has_type(K, c("Function")) && any(c("factorial", "FallingFactorial", "RisingFactorial") %in% .fn_names(K)))
          K <- simplify(K)
        if (!is_zero_num(casoratian(c(kernel, list(K)), n, zero = FALSE))) kernel[[length(kernel) + 1L]] <- K
      }
    }
  }
  kernel <- .key_sort(lapply(kernel, .sort_key), kernel)
  syms <- lapply(seq_along(kernel) - 1L, function(i) Symbol(paste0("C", i)))
  for (i in seq_along(kernel)) result <- .add(list(result, .mul(list(syms[[i]], kernel[[i]]))))
  list(result = result, C = syms)
}

.fn_names <- function(e) {
  out <- if (e$type == "Function") e$name else character()
  for (a in e$args) out <- c(out, .fn_names(a))
  unique(out)
}

#' Building blocks of the recurrence solver
#'
#' The recurrence is `sum(coeffs[[i]] * y(n + i - 1)) = f`, with
#' polynomial coefficients in `n`. `rsolve_poly()` finds its polynomial
#' solutions (Abramov's algorithm), `rsolve_ratio()` its rational solutions
#' and `rsolve_hyper()` its hypergeometric solutions (Petkovsek's
#' algorithm). Free constants are `C0`, `C1`, ... The result is `NULL` when
#' there is no solution of that kind.
#'
#' @param coeffs A list of coefficients, lowest shift first.
#' @param f The inhomogeneous term.
#' @param n The index symbol.
#' @param shift First index for the constant names (`rsolve_poly()` only).
#' @return An expression or `NULL`.
#' @examples
#' n <- Symbol("n", integer = TRUE)
#' rsolve_poly(list(-1, 1), 4*n^3, n)
#' rsolve_ratio(list(-2*n^3 + n^2 + 2*n - 1, 2*n^3 + n^2 - 6*n,
#'                   -2*n^3 - 11*n^2 - 18*n - 9, 2*n^3 + 13*n^2 + 22*n + 8), 0, n)
#' rsolve_hyper(list(-1, -1, 1), 0, n)
#' @name rsolve-parts
NULL

#' @rdname rsolve-parts
#' @export
rsolve_poly <- function(coeffs, f, n, shift = 0L) {
  r <- .rsolve_poly(coeffs, f, as_sym(n), as.integer(shift))
  if (is.null(r)) NULL else r$result
}

#' @rdname rsolve-parts
#' @export
rsolve_ratio <- function(coeffs, f, n) {
  r <- .rsolve_ratio(coeffs, f, as_sym(n))
  if (is.null(r)) NULL else r$result
}

#' @rdname rsolve-parts
#' @export
rsolve_hyper <- function(coeffs, f, n) {
  r <- .rsolve_hyper(coeffs, f, as_sym(n))
  if (is.null(r)) NULL else r$result
}
