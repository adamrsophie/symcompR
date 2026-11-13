# Equation solving, after sympy/solvers/solvers.py, polysys.py and solveset.py.

.eq_to_expr <- function(e) {
  e <- as_sym(e)
  if (e$type == "Equality") return(.sub(e$args[[1L]], e$args[[2L]]))
  if (e$type == "BooleanAtom") return(if (e$name == "true") S_Zero else S_One)
  e
}

#' Solve equations
#'
#' `solve(f, x)` solves `f = 0` (or an [Eq()]) for `x` following SymPy's
#' `solve()`: a single equation in one symbol gives a list of solutions, a
#' linear system a named list (a dict), other systems a list of solution
#' tuples or, with `dict = TRUE`, of named lists. Extra symbols may follow
#' `b`, and SymPy's flags (`dict`, `set`, `check`, `simplify`, `rational`,
#' `exclude`, `implicit`, `manual`, `cubics`, `quartics`, `quintics`, ...)
#' are passed by name. Solutions are simplified, checked and sorted as in
#' SymPy. Inequalities go to [reduce_inequalities()]. This is a method for
#' the base R generic [base::solve()].
#'
#' @param a An equation, expression, list of them, or a character string.
#' @param b A symbol or list of symbols; by default all free symbols.
#' @param ... More symbols, then named flags as in SymPy.
#' @return A list of solutions, a named list, or a list of tuples.
#' @examples
#' x <- Symbol("x"); y <- Symbol("y")
#' solve(x^2 - 4, x)
#' solve(Eq(exp(x), 3), x)
#' solve(list(x + y - 3, x - y - 1), list(x, y))
#' solve(list(x^2 + y^2 - 5, x - y + 1), list(x, y))
#' @export
solve.sym <- function(a, b, ...) .solve_front(a, if (missing(b)) list() else list(b), list(...))

#' @rdname solve.sym
#' @export
solve.character <- function(a, b, ...) {
  e <- sympify(a)
  if (missing(b)) return(solve(e, ...))
  if (is.character(b)) b <- lapply(.expand_symbol_names(b), function(n) .resolve_name(e, n))
  solve(e, b, ...)
}

#' @rdname solve.sym
#' @export
solve.Poly <- function(a, b, ...) .solve_front(a, if (missing(b)) list() else list(b), list(...))

#' @rdname solve.sym
#' @export
solve.list <- function(a, b, ...) {
  a <- lapply(a, function(e) if (is.character(e)) sympify(e) else e)
  .solve_front(a, if (missing(b)) list() else list(b), list(...))
}

.solve_front <- function(f, b, dots) {
  nm <- names(dots)
  if (is.null(nm)) nm <- rep("", length(dots))
  .sp_solve(f, c(b, dots[!nzchar(nm)]), dots[nzchar(nm)])
}

unique_syms <- function(lsts) {
  all <- unlist(lsts, recursive = FALSE)
  if (!length(all)) return(list())
  all[!duplicated(vapply(all, .k, ""))]
}

# ---------------------------------------------------------------------------
# One equation, one unknown

.unique_sols <- function(sols) {
  if (!length(sols)) return(list())
  sols[!duplicated(vapply(sols, .k, ""))]
}

.solve_univariate <- function(f, x, check = TRUE) {
  f <- as_sym(f)
  if (!.has1(f, x)) {
    if (is_zero_num(f)) stop("equation is an identity; every value is a solution")
    return(list())
  }
  nd <- fraction(together(f))
  num <- nd[[1L]]
  den <- nd[[2L]]
  exact <- expand(num)
  exact <- is_polynomial(exact, x)
  sols <- .solve_numer(num, x)
  if (is.null(sols)) stop("could not solve ", .str(f), " for ", .str(x))
  sols <- .unique_sols(sols)
  if (check && length(sols)) {
    sols <- Filter(function(s) {
      if (.has1(den, x) && .vanishes(den, x, s)) return(FALSE)
      if (!.fits_assumptions(s, x)) return(FALSE)
      if (exact) return(TRUE)
      .residual_ok(num, x, s)
    }, sols)
  }
  sols
}

# A candidate contradicting an assumption on the unknown (positive, real,
# integer, ...) is dropped, as SymPy's solve() does.
.fits_assumptions <- function(s, x) {
  af <- if (is_symbol_like(x)) x$extra$facts else NULL
  if (is.null(af)) return(TRUE)
  fs <- .facts(s)
  for (k in c("positive", "negative", "nonnegative", "nonpositive", "real", "integer", "rational", "zero", "nonzero",
              "even", "odd", "prime", "composite", "irrational", "finite")) {
    want <- af[[k]]
    if (is.null(want) || is.na(want)) next
    if (isTRUE(want) && identical(fs[[k]], FALSE)) return(FALSE)
    if (identical(want, FALSE) && isTRUE(fs[[k]])) return(FALSE)
  }
  TRUE
}

.neval_safe <- function(s) tryCatch(.neval(s, list()), error = function(e) 1)

# Does expression e vanish at x = s? Exact zero test first, numeric second.
.vanishes <- function(e, x, s) {
  d <- tryCatch(subs(e, x, s), error = function(err) NULL)
  if (is.null(d)) return(FALSE)
  if (is_zero_num(d)) return(TRUE)
  if (length(free_symbols(d))) return(FALSE)
  z <- tryCatch(simplify(d), error = function(err) d)
  is_zero_num(z)
}

# Accept s as a root of num unless a numeric residual, scaled by the size
# of the terms, clearly does not vanish and simplification cannot prove zero.
.residual_ok <- function(num, x, s) {
  v <- tryCatch(subs(num, x, s), error = function(e) NULL)
  if (is.null(v)) return(TRUE)
  if (is_nan(v) || is_inf_const(v)) return(FALSE)
  if (length(free_symbols(v))) return(TRUE)
  terms <- if (v$type == "Add") v$args else list(v)
  vals <- tryCatch(vapply(terms, function(t) as.complex(.neval(t, list())), 0i), error = function(e) NULL)
  if (is.null(vals)) return(TRUE)
  scale <- sum(base::Mod(vals))
  if (!is.finite(scale)) return(TRUE)
  if (base::Mod(sum(vals)) <= 1e-8 * max(1, scale)) return(TRUE)
  is_zero_num(tryCatch(simplify(v), error = function(e) v))
}

# Solve num = 0 for x; NULL when no method applies.
.solve_numer <- function(num, x) {
  num <- expand(num)
  if (!.has1(num, x)) return(list())
  if (is_polynomial(num, x)) {
    r <- roots(num, x, rootof = TRUE)
    deg <- length(.coeff_list(num, x)) - 1L
    return(as.list(r))
  }
  # factor products
  fac <- factor_sym(num)
  if (fac$type == "Mul" || fac$type == "Pow") {
    parts <- if (fac$type == "Mul") fac$args else list(fac)
    parts <- Filter(function(p) .has1(p, x), parts)
    if (length(parts) > 1L || (length(parts) == 1L && parts[[1L]]$type == "Pow")) {
      out <- list()
      for (p in parts) {
        if (p$type == "Pow" && isTRUE(ask_fact(p$args[[2L]], "positive"))) p <- p$args[[1L]]
        s <- .solve_numer(p, x)
        if (is.null(s)) return(NULL)
        out <- c(out, s)
      }
      return(out)
    }
  }
  # polynomial in a single generator g(x)
  r <- .solve_via_generator(num, x)
  if (!is.null(r)) return(r)
  inv <- .invert(num, list(S_Zero), x)
  if (!is.null(inv)) return(inv)
  NULL
}

.solve_via_generator <- function(num, x) {
  num <- .rewrite_exp_powers(num, x)
  gens <- .find_gens(list(num))
  dep <- Filter(function(g) .has1(g, x), gens)
  if (length(dep) != 1L) return(NULL)
  g <- dep[[1L]]
  if (identical(g$key, x$key)) return(NULL)
  u <- Dummy("u")
  pu <- subs(num, g, u)
  if (.has1(pu, x) || !is_polynomial(pu, u)) return(NULL)
  us <- roots(pu, u, rootof = TRUE)
  out <- list()
  for (val in us) {
    s <- .invert(.sub(g, val), list(S_Zero), x)
    if (is.null(s)) return(NULL)
    out <- c(out, s)
  }
  out
}

# exp(2*x) -> exp(x)^2 so that a single generator appears.
.rewrite_exp_powers <- function(e, x) {
  exps <- atoms(e, "Function")
  exps <- Filter(function(f) f$name == "exp" && .has1(f, x), exps)
  if (length(exps) < 2L) return(e)
  args <- lapply(exps, function(f) f$args[[1L]])
  base <- args[[1L]]
  ks <- lapply(args, function(a) cancel(.div(a, base)))
  if (!all(vapply(ks, is_exact_number, TRUE))) return(e)
  den <- Reduce(gmp::lcm.bigz, lapply(ks, function(k) gmp::denominator(k$value)))
  num_g <- Reduce(gmp::gcd.bigz, lapply(ks, function(k) gmp::numerator(k$value * den)))
  unit <- .mul(list(base, make_rational(gmp::as.bigq(num_g, den))))
  g <- fcall("exp", list(unit))
  rules <- lapply(seq_along(exps), function(i) {
    k <- cancel(.div(args[[i]], unit))
    list(exps[[i]], .pow(g, k))
  })
  xreplace(e, rules)
}

# Invert lhs = rhs for x by peeling off invertible operations.
.invert <- function(lhs, rhss, x) {
  for (iter in 1:50) {
    if (identical(lhs$key, x$key)) return(lapply(rhss, function(r) simplify_if_cheap(r)))
    t <- lhs$type
    if (t == "Add") {
      dep <- Filter(function(a) .has1(a, x), lhs$args)
      indep <- Filter(function(a) !.has1(a, x), lhs$args)
      if (length(dep) != 1L) return(NULL)
      ind <- .add(indep)
      rhss <- lapply(rhss, function(r) .sub(r, ind))
      lhs <- dep[[1L]]
      next
    }
    if (t == "Mul") {
      dep <- Filter(function(a) .has1(a, x), lhs$args)
      indep <- Filter(function(a) !.has1(a, x), lhs$args)
      if (length(dep) != 1L) return(NULL)
      ind <- .mul(indep)
      rhss <- lapply(rhss, function(r) .div(r, ind))
      lhs <- dep[[1L]]
      next
    }
    if (t == "Pow") {
      b <- lhs$args[[1L]]
      e <- lhs$args[[2L]]
      if (.has1(b, x) && !.has1(e, x)) {
        if (e$type == "Integer") {
          n <- as.integer(as.character(e$value))
          new <- list()
          for (r in rhss) {
            if (n > 0) {
              y <- Dummy("y")
              new <- c(new, as.list(roots(.sub(.pow(y, e), r), y)))
            } else {
              new[[length(new) + 1L]] <- .pow(r, .inv(e))
              if (n %% 2L == 0L) new[[length(new) + 1L]] <- .neg(.pow(r, .inv(e)))
            }
          }
          rhss <- new
        } else {
          rhss <- lapply(rhss, function(r) .pow(r, .inv(e)))
        }
        lhs <- b
        next
      }
      if (!.has1(b, x) && .has1(e, x)) {
        rhss <- lapply(rhss, function(r) .div(fcall("log", list(r)), fcall("log", list(b))))
        lhs <- e
        next
      }
      return(NULL)
    }
    if (t == "Function") {
      if (length(lhs$args) != 1L) return(NULL)
      a <- lhs$args[[1L]]
      nm <- lhs$name
      new <- switch(EXPR = nm,
        exp = lapply(rhss, function(r) fcall("log", list(r))),
        log = lapply(rhss, function(r) fcall("exp", list(r))),
        sin = unlist(lapply(rhss, function(r) list(fcall("asin", list(r)), .sub(S_Pi, fcall("asin", list(r))))), recursive = FALSE),
        cos = unlist(lapply(rhss, function(r) list(fcall("acos", list(r)), .sub(.mul(list(S_Two, S_Pi)), fcall("acos", list(r))))), recursive = FALSE),
        tan = lapply(rhss, function(r) fcall("atan", list(r))),
        cot = lapply(rhss, function(r) fcall("acot", list(r))),
        sinh = unlist(lapply(rhss, function(r) list(fcall("asinh", list(r)), .sub(.mul(list(S_I, S_Pi)), fcall("asinh", list(r))))), recursive = FALSE),
        cosh = unlist(lapply(rhss, function(r) list(fcall("acosh", list(r)), .neg(fcall("acosh", list(r))))), recursive = FALSE),
        tanh = lapply(rhss, function(r) fcall("atanh", list(r))),
        asin = lapply(rhss, function(r) fcall("sin", list(r))),
        acos = lapply(rhss, function(r) fcall("cos", list(r))),
        atan = lapply(rhss, function(r) fcall("tan", list(r))),
        asinh = lapply(rhss, function(r) fcall("sinh", list(r))),
        acosh = lapply(rhss, function(r) fcall("cosh", list(r))),
        atanh = lapply(rhss, function(r) fcall("tanh", list(r))),
        LambertW = lapply(rhss, function(r) .mul(list(r, fcall("exp", list(r))))),
        erf = lapply(rhss, function(r) fcall("erfinv", list(r))),
        erfc = lapply(rhss, function(r) fcall("erfcinv", list(r))),
        erfinv = lapply(rhss, function(r) fcall("erf", list(r))),
        erfcinv = lapply(rhss, function(r) fcall("erfc", list(r))),
        Abs = if (isTRUE(ask_fact(a, "real")) || TRUE)
          unlist(lapply(rhss, function(r) list(r, .neg(r))), recursive = FALSE) else NULL,
        NULL)
      if (is.null(new)) return(NULL)
      rhss <- new
      lhs <- a
      next
    }
    return(NULL)
  }
  NULL
}

simplify_if_cheap <- function(e) {
  if (.cops(e) > 60L) return(e)
  tryCatch(simplify(e), error = function(err) e)
}

# ---------------------------------------------------------------------------
# Systems

.is_linear_in <- function(e, syms) {
  e <- expand(e)
  terms <- if (e$type == "Add") e$args else list(e)
  for (t in terms) {
    d <- 0L
    for (s in syms) {
      ds <- .degree_in(t, s)
      if (is.na(ds)) return(FALSE)
      d <- d + ds
    }
    if (d > 1L) return(FALSE)
  }
  TRUE
}

.solve_system <- function(eqs, syms, check = TRUE) {
  nds <- lapply(eqs, function(e) fraction(together(e)))
  dens <- Filter(function(d) !is_number_type(d), lapply(nds, `[[`, 2L))
  sols <- .solve_system_numer(lapply(nds, function(nd) expand(nd[[1L]])), syms)
  if (!check || !length(dens)) return(sols)
  Filter(function(sol) {
    for (d in dens) {
      v <- tryCatch(simplify(subs(d, sol)), error = function(e) NULL)
      if (!is.null(v) && is_zero_num(v)) return(FALSE)
    }
    TRUE
  }, sols)
}

.solve_system_numer <- function(eqs, syms) {
  if (all(vapply(eqs, .is_linear_in, TRUE, syms = syms))) {
    r <- .solve_linear_system(eqs, syms)
    return(if (is.null(r)) list() else list(r))
  }
  if (all(vapply(eqs, function(e) is_polynomial(e, syms), TRUE))) {
    r <- .solve_poly_system_fast(eqs, syms)
    if (!is.null(r)) return(r)
    return(.solve_poly_system(eqs, syms))
  }
  .solve_by_substitution(eqs, syms)
}

.solve_linear_system <- function(eqs, syms) {
  n <- length(syms)
  rows <- lapply(eqs, function(eq) {
    coefs <- lapply(syms, function(s) .linear_coeff(eq, s))
    rest <- expand(.sub(eq, .add(lapply(seq_len(n), function(i) .mul(list(coefs[[i]], syms[[i]]))))))
    c(coefs, list(.neg(rest)))
  })
  sol <- .gauss_solve(rows, n)
  if (is.null(sol)) return(NULL)
  piv <- attr(sol, "pivots")
  free <- setdiff(seq_len(n), piv)
  if (length(free)) return(.solve_linear_param(eqs, syms, piv, free))
  out <- list()
  for (k in piv) out[[syms[[k]]$name]] <- simplify_if_cheap(sol[[k]])
  out
}

# Parametric linear solution: solve for pivot variables treating free
# variables as parameters.
.solve_linear_param <- function(eqs, syms, piv, free) {
  psyms <- syms[piv]
  n <- length(psyms)
  rows <- lapply(eqs, function(eq) {
    coefs <- lapply(psyms, function(s) .linear_coeff(eq, s))
    rest <- expand(.sub(eq, .add(lapply(seq_len(n), function(i) .mul(list(coefs[[i]], psyms[[i]]))))))
    c(coefs, list(.neg(rest)))
  })
  sol <- .gauss_solve(rows, n)
  out <- list()
  for (k in seq_len(n)) out[[psyms[[k]]$name]] <- simplify_if_cheap(sol[[k]])
  out
}

# Polynomial systems: propagate univariate and linear equations, branch on
# roots, fall back to a lex Groebner basis.
.solve_poly_system <- function(eqs, syms, depth = 0L) {
  if (depth > 200L) stop("polynomial system too deep")
  eqs <- lapply(eqs, expand)
  eqs <- Filter(function(e) !is_zero_num(e), eqs)
  if (any(vapply(eqs, function(e) !length(free_symbols(e)) || !any(vapply(syms, function(s) .has1(e, s), TRUE)), TRUE))) {
    # equations without unknowns must vanish identically
    bad <- Filter(function(e) !any(vapply(syms, function(s) .has1(e, s), TRUE)), eqs)
    if (any(vapply(bad, function(e) !is_zero_num(simplify(e)), TRUE))) return(list())
    eqs <- Filter(function(e) any(vapply(syms, function(s) .has1(e, s), TRUE)), eqs)
  }
  if (!length(eqs)) return(list(list()))
  sk <- vapply(syms, .k, "")
  used <- lapply(eqs, function(e) which(vapply(syms, function(s) .has1(e, s), TRUE)))
  nvars <- vapply(used, length, 0L)
  # 1. univariate equations: choose the one of lowest degree
  uni <- which(nvars == 1L)
  if (length(uni)) {
    degs <- vapply(uni, function(i) length(.coeff_list(eqs[[i]], syms[[used[[i]]]])) - 1L, 0L)
    i <- uni[which.min(degs)]
    v <- syms[[used[[i]]]]
    vals <- as.list(roots(eqs[[i]], v, rootof = TRUE))
    if (!length(vals) && min(degs) > 0) return(list())
    out <- list()
    rest_syms <- syms[-used[[i]]]
    for (val in vals) {
      sub_eqs <- lapply(eqs[-i], function(e) expand(subs(e, v, val)))
      rec <- if (length(rest_syms)) .solve_poly_system(sub_eqs, rest_syms, depth + 1L) else {
        if (all(vapply(sub_eqs, function(e) is_zero_num(simplify(e)), TRUE))) list(list()) else list()
      }
      for (s in rec) {
        s2 <- s
        s2[[v$name]] <- val
        out[[length(out) + 1L]] <- s2
      }
    }
    return(out)
  }
  # 2. linear in some variable with a constant nonzero coefficient
  best <- NULL
  for (i in seq_along(eqs)) {
    e <- eqs[[i]]
    for (j in used[[i]]) {
      s <- syms[[j]]
      cl <- tryCatch(.coeff_list(e, s), error = function(err) NULL)
      if (is.null(cl) || length(cl) != 2L) next
      if (!is_number_type(cl[[1L]]) || is_zero_num(cl[[1L]])) next
      size <- .cops(e)
      if (is.null(best) || size < best$size) best <- list(i = i, j = j, size = size, cl = cl)
    }
  }
  if (!is.null(best)) {
    s <- syms[[best$j]]
    val <- expand(.div(.neg(best$cl[[2L]]), best$cl[[1L]]))
    sub_eqs <- lapply(eqs[-best$i], function(e) expand(subs(e, s, val)))
    rest_syms <- syms[-best$j]
    rec <- if (length(rest_syms)) .solve_poly_system(sub_eqs, rest_syms, depth + 1L) else {
      if (all(vapply(sub_eqs, is_zero_num, TRUE))) list(list()) else list()
    }
    out <- list()
    for (sol in rec) {
      v <- val
      for (nm in names(sol)) v <- subs(v, .resolve_name(v, nm), sol[[nm]])
      sol[[s$name]] <- expand(v)
      out[[length(out) + 1L]] <- sol
    }
    return(out)
  }
  # 3. Groebner basis in lex order, then recurse on the triangular basis
  G <- groebner(eqs, syms, order = "lex")
  if (length(G) == 1L && is_number_type(G[[1L]]) && !is_zero_num(G[[1L]])) return(list())
  Gl <- lapply(G, identity)
  if (identical(vapply(Gl, .k, ""), vapply(eqs, .k, ""))) {
    stop("could not solve the polynomial system (infinitely many solutions or unsupported)")
  }
  .solve_poly_system(Gl, syms, depth + 1L)
}

.solve_by_substitution <- function(eqs, syms) {
  # pick an equation with a single unknown
  for (i in seq_along(eqs)) {
    us <- Filter(function(s) .has1(eqs[[i]], s), syms)
    if (length(us) == 1L) {
      v <- us[[1L]]
      vals <- .solve_univariate(eqs[[i]], v)
      out <- list()
      rest <- syms[!vapply(syms, function(s) identical(s$key, v$key), TRUE)]
      for (val in vals) {
        sub_eqs <- lapply(eqs[-i], function(e) subs(e, v, val))
        # with no unknowns left, every other equation must vanish
        if (!length(rest) && !all(vapply(sub_eqs, function(e) is_zero_num(simplify(e)), TRUE))) next
        rec <- if (length(rest)) .solve_system(Filter(function(e) !is_zero_num(e), sub_eqs), rest) else list(list())
        for (s in rec) {
          s[[v$name]] <- val
          out[[length(out) + 1L]] <- s
        }
      }
      return(out)
    }
  }
  # solve one equation for one unknown and substitute
  for (i in seq_along(eqs)) {
    for (v in syms) {
      if (!.has1(eqs[[i]], v)) next
      vals <- tryCatch(.solve_univariate(eqs[[i]], v), error = function(e) NULL)
      if (is.null(vals) || !length(vals)) next
      rest <- syms[!vapply(syms, function(s) identical(s$key, v$key), TRUE)]
      out <- list()
      for (val in vals) {
        sub_eqs <- Filter(function(e) !is_zero_num(e), lapply(eqs[-i], function(e) simplify(subs(e, v, val))))
        rec <- if (length(sub_eqs)) .solve_system(sub_eqs, rest) else list(list())
        for (s in rec) {
          vv <- val
          for (nm in names(s)) vv <- subs(vv, .resolve_name(vv, nm), s[[nm]])
          s[[v$name]] <- simplify(vv)
          out[[length(out) + 1L]] <- s
        }
      }
      return(out)
    }
  }
  stop("could not solve the system")
}

#' Numerical solution of equations
#'
#' `nsolve(f, x, x0)` finds a root near `x0` by Newton's method (secant
#' fallback); for systems give lists of equations, symbols and starting
#' values.
#'
#' @param f An equation or list of equations.
#' @param x A symbol or list of symbols.
#' @param x0 Starting value(s).
#' @param tol Tolerance.
#' @param maxit Maximum iterations.
#' @return A Float or a named list of Floats.
#' @examples
#' x <- Symbol("x")
#' nsolve(cos(x) - x, x, 1)
#' @export
nsolve <- function(f, x, x0, tol = 1e-14, maxit = 100) {
  if (is.list(f) && !inherits(f, "sym")) {
    fs <- lapply(f, .eq_to_expr)
    xs <- .as_gens(list(x))
    J <- lapply(fs, function(fi) lapply(xs, function(v) .diff1(fi, v)))
    cur <- as.complex(unlist(lapply(x0, function(v) if (is_sym(v)) .neval(v, list()) else v)))
    nm <- vapply(xs, function(v) v$name, "")
    ok <- FALSE
    for (it in seq_len(maxit)) {
      env <- stats::setNames(as.list(cur), nm)
      Fv <- vapply(fs, function(fi) as.complex(.neval(fi, env)), 0i)
      Jv <- matrix(vapply(unlist(J, recursive = FALSE), function(d) as.complex(.neval(d, env)), 0i),
                   length(fs), length(xs), byrow = TRUE)
      step <- base::solve(Jv, Fv)
      cur <- cur - step
      if (max(base::Mod(step)) < tol * max(1, max(base::Mod(cur)))) {
        ok <- TRUE
        break
      }
    }
    # as mpmath's findroot(verify = TRUE): the squared residual norm must
    # not exceed tol
    if (ok) {
      env <- stats::setNames(as.list(cur), nm)
      Fv <- vapply(fs, function(fi) as.complex(.neval(fi, env)), 0i)
      ok <- all(is.finite(Fv)) && sum(base::Mod(Fv)^2) <= tol
    }
    if (!ok) stop("nsolve() did not converge in ", maxit, " iterations; try another start point", call. = FALSE)
    return(stats::setNames(lapply(cur, function(v) .round_sym(.cplx_clean(v), 17)), nm))
  }
  fe <- .eq_to_expr(f)
  x <- as_sym(x)
  df <- .diff1(fe, x)
  cur <- if (is.complex(x0)) x0 else if (is_sym(x0)) .neval(x0, list()) else as.numeric(x0)
  ok <- FALSE
  for (it in seq_len(maxit)) {
    env <- stats::setNames(list(cur), x$name)
    fv <- .neval(fe, env)
    dv <- .neval(df, env)
    if (!is.finite(dv) || dv == 0) break
    step <- fv / dv
    cur <- cur - step
    if (base::Mod(step) < tol * max(1, base::Mod(cur))) {
      ok <- TRUE
      break
    }
  }
  # a small step is not enough: as mpmath's findroot(verify = TRUE), the
  # squared residual must not exceed tol
  if (ok) {
    r <- base::Mod(.neval(fe, stats::setNames(list(cur), x$name)))
    ok <- is.finite(r) && r^2 <= tol
  }
  if (!ok) stop("nsolve() did not converge in ", maxit, " iterations; try another start point", call. = FALSE)
  .round_sym(.cplx_clean(cur), 17)
}

# Fast path for polynomial systems with rational coefficients: works on
# qpolys, substitutes rational roots with vectorized evaluation and linear
# eliminations with polynomial substitution. Returns NULL when coefficients
# are not rational, so the expression level solver takes over.
.solve_poly_system_fast <- function(eqs, syms) {
  cv <- tryCatch(.to_qpolys(eqs, syms), error = function(e) NULL)
  if (is.null(cv) || length(cv$gens) != length(syms)) return(NULL)
  n <- length(syms)
  res <- .qp_system_rec(cv$polys, n, depth = 0L)
  if (is.null(res)) return(NULL)
  lapply(res, function(sol) {
    out <- list()
    for (j in seq_len(n)) {
      v <- sol[[j]]
      out[[syms[[j]]$name]] <- if (is.null(v)) syms[[j]] else v
    }
    out
  })
}

# Each solution is a list of length n holding sym values (NULL = free).
.qp_system_rec <- function(P, n, depth) {
  if (depth > 500L) return(NULL)
  P <- Filter(Negate(qp_is_zero), P)
  if (any(vapply(P, qp_is_const, TRUE))) return(list())
  if (!length(P)) return(list(vector("list", n)))
  P <- P[!duplicated(vapply(P, function(p) paste(.mono_keys(p$e), as.character(p$c), collapse = ";"), ""))]
  used <- lapply(P, qp_vars_used)
  nv <- vapply(used, length, 0L)
  uni <- which(nv == 1L)
  if (length(uni)) {
    degs <- vapply(uni, function(i) qp_degree(P[[i]], used[[i]]), 0)
    i <- uni[which.min(degs)]
    j <- used[[i]]
    up <- up_from_qpoly(P[[i]], j)
    # gcd with every other univariate polynomial in the same variable
    for (k in uni) if (k != i && identical(used[[k]], j)) up <- up_gcd(up, up_from_qpoly(P[[k]], j))
    if (up_deg(up) <= 0) return(list())
    zl <- .zp_factor_list(gmp::numerator(up * gmp::as.bigq(.lcm_den(up))))
    rest <- P[-i]
    out <- list()
    for (pr in zl$factors) {
      g <- pr[[1L]]
      if (length(g) == 2L) {
        val <- gmp::as.bigq(-g[2L], g[1L])
        sub <- lapply(rest, qp_eval_var, i = j, v = val)
        rec <- .qp_system_rec(sub, n, depth + 1L)
        if (is.null(rec)) return(NULL)
        for (sol in rec) {
          sol[j] <- list(make_rational(val))
          out[[length(out) + 1L]] <- sol
        }
      } else {
        return(NULL)
      }
    }
    return(out)
  }
  best <- NULL
  for (i in seq_along(P)) {
    p <- P[[i]]
    for (j in used[[i]]) {
      if (qp_degree(p, j) != 1) next
      w <- which(p$e[, j] == 1L)
      lc <- qp_norm(p$e[w, , drop = FALSE] - matrix(as.integer(seq_len(n) == j), length(w), n, byrow = TRUE), p$c[w])
      if (!qp_is_const(lc)) next
      # replacement degree and variable count, to keep substitutions small
      others <- setdiff(used[[i]], j)
      rdeg <- max(0, rowSums(p$e[p$e[, j] == 0L, , drop = FALSE]))
      size <- length(p$c) + 10 * length(others) + 100 * max(0, rdeg - 1)
      if (is.null(best) || size < best$size) best <- list(i = i, j = j, size = size, lc = qp_lc(lc),
                                                          nother = length(others), rdeg = rdeg)
    }
  }
  if (!is.null(best) && (best$nother <= 1L && best$rdeg <= 1)) {
    p <- P[[best$i]]
    j <- best$j
    w <- which(p$e[, j] == 1L)
    restp <- qp(p$e[-w, , drop = FALSE], p$c[-w])
    restp <- qp_norm(restp$e, restp$c)
    val <- qp_scale(restp, -1 / best$lc)
    sub <- lapply(P[-best$i], qp_subs_var, i = j, q = val)
    rec <- .qp_system_rec(sub, n, depth + 1L)
    if (is.null(rec)) return(NULL)
    out <- list()
    for (sol in rec) {
      v <- val
      for (k in qp_vars_used(v)) {
        if (!is.null(sol[[k]]) && is_exact_number(sol[[k]])) v <- qp_eval_var(v, k, sol[[k]]$value)
      }
      if (qp_is_const(v)) {
        sol[j] <- list(make_rational(if (qp_is_zero(v)) gmp::as.bigq(0) else v$c[1L]))
      } else {
        return(NULL)
      }
      out[[length(out) + 1L]] <- sol
    }
    return(out)
  }
  if (getOption("symcompR.trace_solve", FALSE)) {
    message("solve: Groebner on ", length(P), " polynomials, ", sum(vapply(P, function(p) length(p$c), 0L)), " terms")
  }
  key_of <- function(L) sort(vapply(L, function(p) {
    p <- qp_monic(qp_reorder(p, "lex"))
    paste(.mono_keys(p$e), as.character(p$c), collapse = ";")
  }, ""))
  for (ord in c("grevlex", "lex")) {
    G <- .groebner_qp(P, ord)
    G <- lapply(G, qp_reorder, order = "lex")
    if (length(G) == 1L && qp_is_const(G[[1L]])) return(list())
    if (!identical(key_of(G), key_of(P))) return(.qp_system_rec(G, n, depth + 1L))
  }
  NULL
}
