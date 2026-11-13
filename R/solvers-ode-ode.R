# The ordinary differential equation driver, after sympy/solvers/ode/ode.py:
# dsolve, classify_ode, solve_ics, odesimp, ode_sol_simplicity, the handling
# of initial conditions and of redundant solutions, the first order power
# series solution and checkinfsol.

# ode.allhints: the order in which classify_ode() lists the hints
.ODE_ALLHINTS <- c(
  "factorable", "nth_algebraic", "separable", "1st_exact", "1st_linear", "Bernoulli",
  "1st_rational_riccati", "Riccati_special_minus2", "1st_homogeneous_coeff_best",
  "1st_homogeneous_coeff_subs_indep_div_dep", "1st_homogeneous_coeff_subs_dep_div_indep",
  "almost_linear", "linear_coefficients", "separable_reduced", "1st_power_series", "lie_group",
  "nth_linear_constant_coeff_homogeneous", "nth_linear_euler_eq_homogeneous",
  "nth_linear_constant_coeff_undetermined_coefficients",
  "nth_linear_euler_eq_nonhomogeneous_undetermined_coefficients",
  "nth_linear_constant_coeff_variation_of_parameters",
  "nth_linear_euler_eq_nonhomogeneous_variation_of_parameters", "Liouville", "2nd_linear_airy",
  "2nd_linear_bessel", "2nd_linear_bessel_transform", "2nd_hypergeometric", "2nd_hypergeometric_Integral",
  "nth_order_reducible", "2nd_power_series_ordinary", "2nd_power_series_regular",
  "nth_algebraic_Integral", "separable_Integral", "1st_exact_Integral", "1st_linear_Integral",
  "Bernoulli_Integral", "1st_homogeneous_coeff_subs_indep_div_dep_Integral",
  "1st_homogeneous_coeff_subs_dep_div_indep_Integral", "almost_linear_Integral",
  "linear_coefficients_Integral", "separable_reduced_Integral",
  "nth_linear_constant_coeff_variation_of_parameters_Integral",
  "nth_linear_euler_eq_nonhomogeneous_variation_of_parameters_Integral", "Liouville_Integral",
  "2nd_nonlinear_autonomous_conserved", "2nd_nonlinear_autonomous_conserved_Integral")

# ---------------------------------------------------------------------------
# Numbered constants

.ode_const_exclude <- function(eqs) {
  keys <- character()
  for (e in eqs) {
    if (!is_sym(e)) next
    keys <- c(keys, vapply(free_symbols(e), .k, ""))
    for (fn in atoms(e, "Function")) keys <- c(keys, Symbol(fn$name)$key)
  }
  keys
}

# iter_numbered_constants(eq, start, prefix), taken num at a time
.ode_numbered_constants <- function(eqs, num, start = 1L, prefix = "C", exclude = NULL) {
  excl <- c(.ode_const_exclude(eqs), exclude)
  out <- list()
  i <- start
  while (length(out) < num) {
    s <- Symbol(paste0(prefix, i))
    i <- i + 1L
    if (!(s$key %in% excl)) out[[length(out) + 1L]] <- s
  }
  out
}

.ode_as_eq_list <- function(eq) {
  if (is_sym(eq)) {
    if (eq$type == "Tuple") return(eq$args)
    return(list(eq))
  }
  if (is.list(eq)) return(eq)
  .ode_stop("ValueError", paste0("Expected Expr or iterable but got ", format(eq)))
}

# ode.get_numbered_constants: one constant, or a tuple of num of them
.ode_get_numbered_constants <- function(eq, num = 1L, start = 1L, prefix = "C") {
  cs <- .ode_numbered_constants(.ode_as_eq_list(eq), num, start, prefix)
  if (num == 1L) cs[[1L]] else .pytuple(cs)
}

# ---------------------------------------------------------------------------
# Initial conditions: SymPy's {f(x0): value, f(x).diff(x, n).subs(x, x0): value}
# as a list of list(key, value). The package's earlier form list(point,
# order, value) is accepted too.

.ode_norm_ics <- function(ics, func) {
  if (is.null(ics)) return(NULL)
  if (inherits(ics, "pydict")) return(lapply(seq_along(ics$keys), function(i) list(as_sym(ics$keys[[i]]), as_sym(ics$values[[i]]))))
  if (is_sym(ics)) .ode_stop("ValueError", "Enter boundary conditions of the form ics={f(point): value, f(x).diff(x, order).subs(x, point): value}")
  lapply(ics, function(p) {
    if (length(p) == 3L) {
      x <- func$args[[1L]]
      pt <- as_sym(p[[1L]]); ord <- as.integer(as.numeric(p[[2L]])); val <- as_sym(p[[3L]])
      key <- if (ord == 0L) xreplace(func, list(list(x, pt))) else Subs(.ode_d(func, x, ord), x, pt)
      return(list(key, val))
    }
    list(as_sym(p[[1L]]), as_sym(p[[2L]]))
  })
}

.ode_ics_legacy <- function(ics) !is.null(ics) && !inherits(ics, "pydict") && is.list(ics) && any(vapply(ics, length, 0L) == 3L)

# ---------------------------------------------------------------------------
# dsolve

#' Solve ordinary differential equations
#'
#' `dsolve(eq, func)` solves the ordinary differential equation `eq` (an
#' expression taken equal to zero, or an [Eq()]) for the unknown function
#' `func`, as SymPy's `dsolve()`. The hints of [classify_ode()] are tried in
#' SymPy's order; `hint` names one of them, or the meta hints `"default"`,
#' `"all"`, `"all_Integral"` and `"best"`. A list of equations is solved as
#' a system of ODEs.
#'
#' @param eq The equation, or a list of equations for a system.
#' @param func The unknown function such as `f(x)`, or a list of functions
#'   for a system; detected from the equation when `NULL`.
#' @param hint The solution method.
#' @param simplify Whether to simplify the solution with `odesimp()`.
#' @param ics Initial conditions: a list of pairs `list(f(x0), value)` or
#'   `list(Subs(Derivative(f(x), x, n), x, x0), value)`; the triple form
#'   `list(x0, n, value)` is accepted too.
#' @param xi,eta Infinitesimals for the `lie_group` hint.
#' @param x0 The point of expansion for the power series hints.
#' @param n The number of power series terms.
#' @param ... Further options passed on to the hints.
#' @return An equation, or a list of equations; a named list for the `"all"`
#'   hints.
#' @examples
#' x <- Symbol("x"); f <- Function("f")
#' dsolve(Derivative(f(x), x, 2) - 3 * Derivative(f(x), x) + 2 * f(x), f(x))
#' dsolve(Derivative(f(x), x) - f(x), f(x), ics = list(list(f(0), 1)))
#' @export
dsolve <- function(eq, func = NULL, hint = "default", simplify = TRUE, ics = NULL, xi = NULL, eta = NULL,
                   x0 = 0, n = 6, ...) {
  kw <- list(...)
  if (is.list(eq) && !is_sym(eq)) return(.ode_dsolve_iterable(lapply(eq, as_sym), func, ics))
  eq <- as_sym(eq)
  if (eq$type == "Tuple") return(.ode_dsolve_iterable(eq$args, func, ics))
  if (!is.null(func)) func <- as_sym(func)
  if (is.character(func)) func <- NULL
  if (.ode_ics_legacy(ics) && is.null(func)) func <- .ode_preprocess(.eq_to_expr(eq))[[2L]]
  ics <- .ode_norm_ics(ics, func)
  given_hint <- hint
  hints <- .ode_desolve(eq, func, hint = hint, simplify = TRUE, ics = ics,
                        prep = !isFALSE(kw$prep), kw = list(xi = xi, eta = eta, x0 = as_sym(x0), n = n))
  if (!is.null(hints$eq)) eq <- hints$eq
  hints$eq <- NULL
  all_ <- isTRUE(hints$all)
  hints$all <- NULL
  if (all_) {
    retdict <- list()
    failed <- list()
    gethints <- .ode_classify(eq, NULL, dict = TRUE, kw = list(hint = "all"))
    orderedhints <- gethints$ordered_hints
    last <- NULL
    for (h in names(hints)) {
      last <- h
      rv <- tryCatch(.ode_helper_simplify(eq, h, hints[[h]], simplify), NotImplementedError = function(e) e)
      if (inherits(rv, "condition")) failed[[h]] <- rv else retdict[[h]] <- rv
    }
    func <- hints[[last]]$func
    simp <- vapply(retdict, function(s) .ode_sol_simplicity(s, func, trysolving = !simplify), 0)
    best <- retdict[[which.min(simp)]]
    if (given_hint == "best") return(best)
    keys <- names(retdict)
    retdict$best <- best
    for (i in orderedhints) {
      if (i %in% keys && .ode_same(best, retdict[[i]])) {
        retdict$best_hint <- i
        break
      }
    }
    retdict$default <- gethints$default
    retdict$order <- gethints$order
    for (h in names(failed)) retdict[[h]] <- failed[[h]]
    return(.pydict(as.list(names(retdict)), unname(retdict)))
  }
  .ode_helper_simplify(eq, hints$hint, hints, simplify, ics = ics)
}

# structural equality of solutions, lists of solutions included
.ode_same <- function(a, b) {
  if (is.null(a) || is.null(b)) return(FALSE)
  if (is_sym(a) && is_sym(b)) return(identical(a$key, b$key))
  if (is_sym(a) || is_sym(b)) return(FALSE)
  length(a) == length(b) && all(vapply(seq_along(a), function(i) .ode_same(a[[i]], b[[i]]), TRUE))
}

.ode_sol_list <- function(s) if (is_sym(s)) list(s) else s

# ode._helper_simplify
.ode_helper_simplify <- function(eq, hint, match, simplify = TRUE, ics = NULL) {
  r <- match
  func <- r$func
  order <- r$order
  m <- r[[hint]]
  free <- free_symbols(eq)
  fkeys <- vapply(free, .k, "")
  cons <- function(s) Filter(function(c) !(c$key %in% fkeys), free_symbols(s))
  is_solver <- inherits(m, "ode_solver")
  base <- sub("_Integral$", "", hint)
  if (simplify) {
    sols <- if (is_solver) .ode_get_general_solution(m) else .ode_series_solver(base)(eq, func, order, m)
    if (!is_sym(sols)) {
      rv <- list()
      for (s in sols) {
        simp <- .ode_odesimp(eq, s, func, hint)
        rv <- c(rv, .ode_sol_list(simp))
      }
    } else rv <- .ode_odesimp(eq, sols, func, hint)
  } else {
    exprs <- if (is_solver) .ode_get_general_solution(m, simplify = FALSE) else {
      m$simplify <- FALSE
      .ode_series_solver(base)(eq, func, order, m)
    }
    rv <- if (!is_sym(exprs)) lapply(exprs, .ode_handle_Integral, func = func, hint = hint) else .ode_handle_Integral(exprs, func, hint)
  }
  if (!is_sym(rv)) {
    if (simplify) rv <- .ode_remove_redundant_solutions(eq, rv, order, func$args[[1L]])
    if (length(rv) == 1L) rv <- rv[[1L]]
  }
  if (length(ics) && !grepl("power_series", hint, fixed = TRUE)) {
    rv1 <- list()
    for (s in .ode_sol_list(rv)) {
      sc <- tryCatch(.ode_solve_ics(list(s), list(func), cons(s), ics), ValueError = function(e) NULL)
      if (is.null(sc)) next
      rv1[[length(rv1) + 1L]] <- .ode_subs_dict(s, sc)
    }
    rv_c <- .ode_get_constant_solutions(eq, func, ics)
    if (length(rv_c)) {
      rv1 <- Filter(function(s) !(s$type == "BooleanAtom" && s$name == "false"), rv1)
      for (s in rv_c) if (!any(vapply(rv1, function(z) identical(z$key, s$key), TRUE))) rv1[[length(rv1) + 1L]] <- s
    }
    if (simplify && length(rv1) > 1L) rv1 <- .ode_remove_redundant_solutions(eq, rv1, order, func$args[[1L]])
    if (length(rv1) == 1L) return(rv1[[1L]])
    rv <- rv1
  }
  rv
}

.ode_subs_dict <- function(e, d) {
  if (!length(d$keys)) return(e)
  subs(e, lapply(seq_along(d$keys), function(i) list(d$keys[[i]], d$values[[i]])))
}

# ode.solve_ics: the values of the constants for the initial conditions,
# a dictionary keyed by the constants
.ode_solve_ics <- function(sols, funcs, constants, ics) {
  ics <- .ode_norm_ics(ics, funcs[[1L]])
  x <- funcs[[1L]]$args[[1L]]
  diff_sols <- list()
  subs_sols <- list()
  for (p in ics) {
    funcarg <- p[[1L]]
    value <- p[[2L]]
    if (.is_undefined_application(funcarg)) {
      x0 <- funcarg$args[[1L]]
      matching_func <- Filter(function(f) f$name == funcarg$name, funcs)[[1L]]
      S <- sols
    } else if (funcarg$type %in% c("Subs", "Derivative")) {
      if (funcarg$type == "Subs") funcarg <- doit(funcarg)
      if (funcarg$type == "Subs") {
        deriv <- funcarg$args[[1L]]
        x0 <- .subs_parts(funcarg)[[2L]][[1L]]
        variables <- .ode_deriv_vars(deriv)
        matching_func <- deriv
      } else {
        deriv <- funcarg
        x0 <- .ode_deriv_vars(funcarg)[[1L]]
        variables <- rep(list(x), length(.ode_deriv_vars(funcarg)))
        matching_func <- subs(deriv, x0, x)
      }
      fname <- deriv$args[[1L]]$name
      for (sol in sols) {
        if (length(Filter(function(a) a$name == fname, atoms(sol, "Function")))) {
          spec <- lapply(variables, function(v) list(v, 1L))
          dl <- .diff_spec_eval(sol$args[[1L]], spec)
          dr <- .diff_spec_eval(sol$args[[2L]], spec)
          diff_sols[[length(diff_sols) + 1L]] <- Eq(dl, dr)
        }
      }
      S <- diff_sols
    } else .ode_nie("Unrecognized initial condition")
    for (sol in S) {
      if (.has1(sol, matching_func)) {
        sol2 <- subs(sol, x, x0)
        sol2 <- subs(sol2, funcarg, value)
        if (sol2$type != "BooleanAtom" || !length(subs_sols)) {
          subs_sols <- Filter(function(s) s$type != "BooleanAtom", subs_sols)
          subs_sols[[length(subs_sols) + 1L]] <- sol2
        }
      }
    }
  }
  solved <- tryCatch(.sp_solve(subs_sols, list(constants), list(dict = TRUE, `_raw` = TRUE)),
                     NotImplementedError = function(e) list())
  if (isTRUE(solved)) .ode_stop("ValueError", "Initial conditions did not produce any solutions for constants. Perhaps they are degenerate.")
  if (!length(solved)) .ode_stop("ValueError", "Couldn't solve for initial conditions")
  if (length(solved) > 1L) .ode_nie("Initial conditions produced too many solutions for constants")
  .pydict(solved[[1L]]$k, solved[[1L]]$v)
}

# ode._get_constant_solutions
.ode_get_constant_solutions <- function(eq, func, ics) {
  ics <- .ode_norm_ics(ics, func)
  x <- func$args[[1L]]
  const_val <- NULL
  for (p in ics) {
    funcarg <- p[[1L]]
    value <- p[[2L]]
    if (funcarg$type == "Subs") {
      funcarg <- doit(funcarg)
      if (funcarg$type == "Subs") funcarg <- funcarg$args[[1L]]
    }
    if (.is_undefined_application(funcarg)) {
      if (func$name != funcarg$name) next
      if (is.null(const_val)) const_val <- value
      else if (!identical(const_val$key, value$key)) return(list())
    } else if (funcarg$type == "Derivative") {
      if (func$name != funcarg$args[[1L]]$name) next
      if (!is_zero_num(value)) return(list())
    } else .ode_nie("Unrecognized initial condition")
  }
  eq <- .eq_to_expr(eq)
  eq <- .map_bottom_up(eq, function(e) if (e$type == "Derivative" && identical(e$args[[1L]]$key, func$key)) S_Zero else e)
  if (is_zero_num(eq)) {
    if (is.null(const_val)) return(list(Eq(func, .ode_get_numbered_constants(eq, 1L))))
    return(list(Eq(func, const_val)))
  }
  result <- list()
  for (sol in .sp_solve(eq, list(func), list())) {
    if (!.has1(sol, x) && (is.null(const_val) || identical(sol$key, const_val$key))) result[[length(result) + 1L]] <- Eq(func, sol)
  }
  result
}

# ---------------------------------------------------------------------------
# classify_ode

#' Classify an ordinary differential equation
#'
#' `classify_ode(eq, func)` returns the hints, in SymPy's order, of the
#' methods [dsolve()] can use to solve the ODE `eq`, as SymPy's
#' `classify_ode()`. With `dict = TRUE` it returns a dictionary of the
#' matching hints with the keys `order`, `default` and `ordered_hints`.
#'
#' @param eq The equation.
#' @param func The unknown function; detected when `NULL`.
#' @param dict Whether to return the dictionary of matches.
#' @param ics Initial conditions, as for [dsolve()].
#' @param prep Whether to preprocess the equation.
#' @param xi,eta Infinitesimals for the `lie_group` hint.
#' @param n The number of power series terms.
#' @param ... Further options: `hint` and `x0`.
#' @return A character vector of hints, or a dictionary.
#' @examples
#' x <- Symbol("x"); f <- Function("f")
#' classify_ode(Derivative(f(x), x) - f(x), f(x))
#' @export
classify_ode <- function(eq, func = NULL, dict = FALSE, ics = NULL, prep = TRUE, xi = NULL, eta = NULL, n = NULL, ...) {
  eq <- as_sym(eq)
  if (!is.null(func)) func <- as_sym(func)
  if (.ode_ics_legacy(ics) && is.null(func)) func <- .ode_preprocess(.eq_to_expr(eq))[[2L]]
  ics <- .ode_norm_ics(ics, func)
  kw <- list(...)
  if (!is.null(kw$x0)) kw$x0 <- as_sym(kw$x0)
  r <- .ode_classify(eq, func, dict = dict, ics = ics, prep = prep, xi = xi, eta = eta, n = n, kw = kw)
  if (!dict) return(r)
  .ode_hints_dict(r)
}

# the R list of matches as a Python dictionary
.ode_hints_dict <- function(r) {
  vals <- lapply(r, function(v) if (is.list(v) && !is_sym(v) && !inherits(v, c("ode_solver", "pytuple")) && !is.null(names(v)))
    .pydict(as.list(names(v)), unname(v)) else v)
  .pydict(as.list(names(r)), unname(vals))
}

.ode_has_bad <- function(e) .has1(e, S_Infinity) || .has1(e, S_ComplexInfinity) || .has1(e, S_NaN) || .has1(e, S_NegativeInfinity)

.ode_classify <- function(eq, func = NULL, dict = FALSE, ics = NULL, prep = TRUE, xi = NULL, eta = NULL, n = NULL, kw = list()) {
  if (!is.null(func) && length(func$args) != 1L)
    .ode_stop("ValueError", paste0("dsolve() and classify_ode() only work with functions of one variable, not ", .str(func)))
  if (eq$type == "Equality") eq <- .sub(eq$args[[1L]], eq$args[[2L]])
  eq_orig <- eq
  if (prep || is.null(func)) {
    pp <- .ode_preprocess(eq, func)
    eq <- pp[[1L]]
    if (is.null(func)) func <- pp[[2L]]
  }
  x <- func$args[[1L]]
  fname <- func$name
  fx <- func
  y <- Dummy("y")
  terms <- if (is.null(n)) 5L else n
  order <- .ode_order(eq, fx)
  matching_hints <- list(order = order)
  df <- .ode_d(fx, x)
  d2 <- .ode_d(fx, x, 2L)
  a <- Wild("a", exclude = list(fx))
  d <- Wild("d", exclude = list(df, d2))
  e <- Wild("e", exclude = list(df))
  c1 <- Wild("c1", exclude = list(x))
  a3 <- Wild("a3", exclude = list(fx, df, d2))
  b3 <- Wild("b3", exclude = list(fx, df, d2))
  c3 <- Wild("c3", exclude = list(fx, df, d2))
  boundary <- list()
  C1 <- Symbol("C1")
  for (p in ics) {
    funcarg <- p[[1L]]
    val <- p[[2L]]
    if (funcarg$type %in% c("Subs", "Derivative")) {
      if (funcarg$type == "Subs") {
        deriv <- funcarg$args[[1L]]
        sp <- .subs_parts(funcarg)
        old <- sp[[1L]][[1L]]
        new <- sp[[2L]][[1L]]
      } else {
        deriv <- funcarg
        old <- x
        new <- .ode_deriv_vars(funcarg)[[1L]]
      }
      dv <- if (deriv$type == "Derivative") .ode_deriv_vars(deriv) else list()
      if (deriv$type == "Derivative" && .is_undefined_application(deriv$args[[1L]]) && deriv$args[[1L]]$name == fname &&
          length(deriv$args[[1L]]$args) == 1L && identical(old$key, x$key) && !.has1(new, x) &&
          all(vapply(dv, function(i) identical(i$key, dv[[1L]]$key), TRUE)) && !.has1(val, x)) {
        dorder <- .ode_order(deriv, x)
        temp <- paste0("f", dorder)
        boundary[[temp]] <- new
        boundary[[paste0(temp, "val")]] <- val
      } else .ode_stop("ValueError", "Invalid boundary conditions for Derivatives")
    } else if (.is_undefined_application(funcarg)) {
      if (funcarg$name == fname && length(funcarg$args) == 1L && !.has1(funcarg$args[[1L]], x) && !.has1(val, x)) {
        boundary$f0 <- funcarg$args[[1L]]
        boundary$f0val <- val
      } else .ode_stop("ValueError", "Invalid boundary conditions for Function")
    } else .ode_stop("ValueError", "Enter boundary conditions of the form ics={f(point): value, f(x).diff(x, order).subs(x, point): value}")
  }
  ode <- .ode_problem(eq_orig, func, x, prep = prep, params = list(xi = xi, eta = eta))
  user_hint <- if (is.null(kw$hint)) "default" else kw$hint
  early_exit <- user_hint == "default"
  user_hint <- sub("_Integral$", "", user_hint)
  user_map <- .ODE_SOLVER_MAP
  if (!(user_hint %in% c("default", "all", "all_Integral", "best")) && user_hint %in% names(.ODE_SOLVER_MAP))
    user_map <- .ODE_SOLVER_MAP[user_hint]
  for (hint in names(user_map)) {
    solver <- .ode_solver_new(user_map[[hint]], ode)
    if (.ode_matches(solver)) {
      matching_hints[[hint]] <- solver
      if (isTRUE(solver$def$has_integral)) matching_hints[[paste0(hint, "_Integral")]] <- solver
      if (dict && early_exit) {
        matching_hints$default <- hint
        return(matching_hints)
      }
    }
  }
  eq <- expand(eq)
  reduced_eq <- NULL
  if (eq$type == "Add") {
    deriv_coef <- coeff(eq, .ode_d(fx, x, order))
    if (!is_one_num(deriv_coef) && !is_zero_num(deriv_coef)) {
      r <- sym_match(deriv_coef, a * .pow(fx, c1))
      if (!is.null(r) && .ode_true(r$c1)) {
        den <- .pow(fx, r$c1)
        reduced_eq <- .add(lapply(eq$args, function(arg) .div(arg, den)))
      }
    }
  }
  if (is.null(reduced_eq) || is_zero_num(reduced_eq)) reduced_eq <- eq
  if (order == 1L) {
    r <- sym_match(collect(eq, df, exact = TRUE), d + e * df)
    if (length(r)) {
      rd <- subs(r$d %||% S_Zero, fx, y)
      re <- subs(r$e %||% S_Zero, fx, y)
      point <- if (is.null(boundary$f0)) S_Zero else boundary$f0
      value <- if (is.null(boundary$f0val)) C1 else boundary$f0val
      check <- cancel(.div(rd, re))
      check1 <- subs(check, list(list(x, point), list(y, value)))
      if (!.ode_has_bad(check1)) {
        check2 <- subs(.diff1(check1, x), list(list(x, point), list(y, value)))
        if (!.ode_has_bad(check2))
          matching_hints[["1st_power_series"]] <- list(d = rd, e = re, y = y, terms = terms, f0 = point, f0val = value)
      }
    }
  } else if (order == 2L) {
    deq <- a3 * d2 + b3 * df + c3 * fx
    r <- sym_match(collect(reduced_eq, list(d2, df, fx)), deq)
    ordinary <- FALSE
    if (length(r)) {
      if (!all(vapply(r, function(v) is_polynomial(v), TRUE))) {
        nd <- .as_numer_denom(reduced_eq)
        reduced_eq <- expand(nd[[1L]])
        r <- sym_match(collect(reduced_eq, list(d2, df, fx)), deq)
      }
    }
    if (length(r) && !is_zero_num(r$a3 %||% S_Zero)) {
      ra3 <- r$a3; rb3 <- r$b3 %||% S_Zero; rc3 <- r$c3 %||% S_Zero
      p <- cancel(.div(rb3, ra3))
      q <- cancel(.div(rc3, ra3))
      point <- if (is.null(kw$x0)) S_Zero else kw$x0
      check <- subs(p, x, point)
      if (!.ode_has_bad(check)) {
        check <- subs(q, x, point)
        if (!.ode_has_bad(check)) {
          ordinary <- TRUE
          matching_hints[["2nd_power_series_ordinary"]] <- list(a3 = ra3, b3 = rb3, c3 = rc3, x0 = point, terms = terms)
        }
      }
      if (!ordinary) {
        p <- cancel(.mul(list(.sub(x, point), p)))
        check <- subs(p, x, point)
        if (!.ode_has_bad(check)) {
          q <- cancel(.mul(list(.pow(.sub(x, point), S_Two), q)))
          check <- subs(q, x, point)
          if (!.ode_has_bad(check))
            matching_hints[["2nd_power_series_regular"]] <- list(p = p, q = q, x0 = point, terms = terms)
        }
      }
    }
  }
  retlist <- .ODE_ALLHINTS[.ODE_ALLHINTS %in% names(matching_hints)]
  if (dict) {
    matching_hints["default"] <- list(if (length(retlist)) retlist[[1L]] else NULL)
    matching_hints$ordered_hints <- .pytuple(as.list(retlist))
    return(matching_hints)
  }
  retlist
}

# ---------------------------------------------------------------------------
# odesimp and friends

.ode_eq_map <- function(eq, fn) {
  if (eq$type == "Equality") return(Eq(fn(eq$args[[1L]]), fn(eq$args[[2L]])))
  fn(eq)
}

# ode._handle_Integral
.ode_handle_Integral <- function(expr, func, hint) {
  if (hint == "nth_linear_constant_coeff_homogeneous") return(expr)
  if (!endsWith(hint, "_Integral")) return(doit(expr))
  expr
}

# ode.odesimp
.ode_odesimp <- function(ode, eq, func, hint) {
  x <- func$args[[1L]]
  C1 <- .ode_get_numbered_constants(eq, 1L)
  okeys <- vapply(free_symbols(ode), .k, "")
  constants <- Filter(function(s) !(s$key %in% okeys), free_symbols(eq))
  eq <- .ode_handle_Integral(eq, func, hint)
  if (startsWith(hint, "nth_linear_euler_eq_nonhomogeneous")) eq <- simplify(eq)
  if (eq$type != "Equality") .ode_stop("TypeError", "eq should be an instance of Equality")
  if (length(constants)) {
    reps <- lapply(constants, function(i) list(i, Dummy(i$name, nonzero = TRUE)))
    eq <- xreplace(eq, reps)
    eq <- xreplace(eq, lapply(reps, rev))
  }
  eq <- constantsimp(eq, constants)
  if (eq$type != "Equality") eqs <- list(eq)
  else {
    lhs <- eq$args[[1L]]
    rhs <- eq$args[[2L]]
    if (identical(rhs$key, func$key) && !.has1(lhs, func)) {
      eq <- Eq(rhs, lhs)
      lhs <- eq$args[[1L]]; rhs <- eq$args[[2L]]
    }
    if (eq$type == "Equality" && identical(lhs$key, func$key) && !.has1(rhs, func)) eqs <- list(eq)
    else {
      floats <- length(Filter(function(z) z$type == "Float", atoms(eq))) > 0L
      flags <- list(force = TRUE)
      if (floats) flags$rational <- FALSE
      eqsol <- tryCatch({
        s <- .sp_solve(eq, list(func), flags)
        if (!length(s)) NULL else s
      }, NotImplementedError = function(err) NULL, PolynomialError = function(err) NULL)
      if (is.null(eqsol)) eqs <- list(eq)
      else {
        .expand <- function(expr) {
          nd <- .as_numer_denom(expr)
          if (nd[[2L]]$type == "Add") return(expr)
          powsimp(expand(expr), combine = "exp", deep = TRUE)
        }
        eqs <- lapply(eqsol, function(t) Eq(func, .expand(t)))
      }
      if (startsWith(hint, "1st_homogeneous_coeff")) {
        eqs <- lapply(eqs, function(eqi) {
          newi <- .ode_eq_map(eqi, function(z) logcombine(z, force = TRUE))
          if (newi$type == "Equality" && newi$args[[1L]]$type == "Function" && newi$args[[1L]]$name == "log" &&
              is_zero_num(newi$args[[2L]]))
            newi <- Eq(.div(newi$args[[1L]]$args[[1L]], C1), C1)
          newi
        })
      }
    }
  }
  ofree <- free_symbols(ode)
  eqs <- lapply(eqs, function(eqi) {
    eqi <- constantsimp(eqi, constants)
    constant_renumber(eqi, ofree)
  })
  if (length(eqs) == 1L) return(eqs[[1L]])
  eqs
}

# ode.ode_sol_simplicity
.ode_sol_simplicity <- function(sol, func, trysolving = TRUE) {
  if (!is_sym(sol)) {
    for (i in sol) if (.ode_sol_simplicity(i, func, trysolving) == Inf) return(Inf)
    return(nchar(paste0("[", paste(vapply(sol, .sp_py_str, ""), collapse = ", "), "]")))
  }
  if (.has_type(sol, "Integral")) return(Inf)
  lhs <- sol$args[[1L]]
  rhs <- sol$args[[2L]]
  if ((identical(lhs$key, func$key) && !.has1(rhs, func)) || (identical(rhs$key, func$key) && !.has1(lhs, func))) return(-2)
  if (trysolving) {
    ok <- tryCatch(length(.sp_solve(sol, list(func), list())) > 0L, NotImplementedError = function(e) FALSE)
    if (ok) return(-1)
  }
  nchar(.sp_py_str(sol))
}

# ode._extract_funcs
.ode_extract_funcs <- function(eqs) {
  funcs <- list()
  for (eq in eqs) {
    derivs <- Filter(function(z) z$type == "Derivative", .ode_preorder(eq))
    for (d in derivs) funcs <- c(funcs, .ode_undef_atoms(d))
  }
  if (length(funcs)) funcs[!duplicated(vapply(funcs, .k, ""))] else funcs
}

.ode_preorder <- function(e) {
  out <- list()
  walk <- function(z) {
    out[[length(out) + 1L]] <<- z
    for (a in .sympy_args(z)) walk(a)
  }
  walk(e)
  out
}

# ode._remove_redundant_solutions
.ode_remove_redundant_solutions <- function(eq, solns, order, var) {
  uniq <- list()
  for (s1 in solns) {
    keep <- TRUE
    for (s2 in uniq) {
      if (.ode_is_special_case_of(s1, s2, eq, order, var)) {
        keep <- FALSE
        break
      } else if (.ode_is_special_case_of(s2, s1, eq, order, var)) {
        uniq <- Filter(function(z) !identical(z$key, s2$key), uniq)
      }
    }
    if (keep) uniq[[length(uniq) + 1L]] <- s1
  }
  uniq
}

.ode_is_special_case_of <- function(soln1, soln2, eq, order, var) {
  soln1 <- .sub(soln1$args[[2L]], soln1$args[[1L]])
  soln2 <- .sub(soln2$args[[2L]], soln2$args[[1L]])
  h1 <- .has_type(soln1, "Order")
  h2 <- .has_type(soln2, "Order")
  if (h1 && h2) {
    if (identical(getO(soln1)$key, getO(soln2)$key)) {
      soln1 <- removeO(soln1)
      soln2 <- removeO(soln2)
    } else return(FALSE)
  } else if (h1 || h2) return(FALSE)
  ekeys <- vapply(free_symbols(eq), .k, "")
  constants1 <- Filter(function(s) !(s$key %in% ekeys), free_symbols(soln1))
  constants2 <- Filter(function(s) !(s$key %in% ekeys), free_symbols(soln2))
  new1 <- .ode_numbered_constants(list(soln1, soln2), length(constants1))
  for (i in seq_along(constants1)) soln1 <- subs(soln1, constants1[[i]], new1[[i]])
  lhs <- soln1
  rhs <- soln2
  eqns <- list(Eq(lhs, rhs))
  for (k in seq_len(order - 1L)) {
    lhs <- .diff1(lhs, var)
    rhs <- .diff1(rhs, var)
    eqns[[length(eqns) + 1L]] <- Eq(lhs, rhs)
  }
  if (any(vapply(eqns, function(z) z$type == "BooleanAtom" && z$name == "false", TRUE))) return(FALSE)
  eqns <- Filter(function(z) !(z$type == "BooleanAtom" && z$name == "true"), eqns)
  csols <- tryCatch(.sp_solve(eqns, list(constants2), list(dict = TRUE, `_raw` = TRUE)),
                    NotImplementedError = function(e) NULL)
  if (is.null(csols)) return(FALSE)
  for (cs in csols) {
    for (z in eqns) {
      z <- .sub(z$args[[2L]], z$args[[1L]])
      if (!isTRUE(.sp_checksol(z, cs, list()))) return(FALSE)
    }
  }
  for (cs in csols) if (!any(vapply(cs$v, function(c) .has1(c, var), TRUE))) return(TRUE)
  FALSE
}

# ---------------------------------------------------------------------------
# Power series solutions

.ode_series_solver <- function(base) {
  switch(base,
    "1st_power_series" = .ode_1st_power_series,
    "2nd_power_series_ordinary" = .ode_2nd_power_series_ordinary,
    "2nd_power_series_regular" = .ode_2nd_power_series_regular,
    .ode_nie(paste0("no solver for the hint ", base)))
}

# ode.ode_1st_power_series
.ode_1st_power_series <- function(eq, func, order, match) {
  x <- func$args[[1L]]
  y <- match$y
  h <- .neg(.div(match$d, match$e))
  point <- match$f0
  value <- match$f0val
  terms <- match$terms
  F <- h
  if (!.ode_true(h)) return(Eq(func, value))
  series <- value
  at <- function(z) subs(z, list(list(x, point), list(y, value)))
  if (terms > 1L) {
    hc <- at(h)
    if (.has1(hc, S_Infinity) || .has1(hc, S_NaN) || .has1(hc, S_ComplexInfinity)) return(Eq(func, S_Infinity))
    else if (.ode_true(hc)) series <- .add(list(series, .mul(list(hc, .sub(x, point)))))
  }
  for (factcount in seq(2L, length.out = max(0L, terms - 2L))) {
    Fnew <- .add(list(.diff1(F, x), .mul(list(.diff1(F, y), h))))
    Fnewc <- at(Fnew)
    if (.ode_has_bad(Fnewc)) return(Eq(func, S_Infinity))
    series <- .add(list(series, .div(.mul(list(Fnewc, .pow(.sub(x, point), .int(factcount)))), fcall("factorial", list(.int(factcount))))))
    F <- Fnew
  }
  series <- .add(list(series, O(.pow(x, .int(terms)))))
  Eq(func, series)
}

# the constant term of an expression in n (as_independent(n)[0])
.ode_n_offset <- function(e, n) .sp_indep(e, list(n), as_add = TRUE)[[1L]]

# ode.ode_2nd_power_series_ordinary
.ode_2nd_power_series_ordinary <- function(eq, func, order, match) {
  x <- func$args[[1L]]
  C <- .ode_get_numbered_constants(eq, num = 2L)
  C0 <- C[[1L]]; C1 <- C[[2L]]
  n <- Dummy("n", integer = TRUE)
  s <- Wild("s")
  k <- Wild("k", exclude = list(x))
  x0 <- match$x0
  terms <- as.integer(match$terms)
  p <- match$a3; q <- match$b3; r <- match$c3
  recurr <- function(a) fcall("r", list(a))
  sk <- list(); sv <- list()
  sset <- function(key, v) {
    j <- match(key$key, vapply(sk, .k, ""))
    if (is.na(j)) { sk[[length(sk) + 1L]] <<- key; sv[[length(sv) + 1L]] <<- v } else sv[[j]] <<- v
  }
  coefflist <- list(list(recurr(n), r), list(.mul(list(n, recurr(n))), q),
                    list(.mul(list(n, .sub(n, S_One), recurr(n))), p))
  for (index in 0:2) {
    coeff <- coefflist[[index + 1L]]
    if (!.ode_true(coeff[[2L]])) next
    f2 <- powsimp(expand(subs(.mul(list(coeff[[2L]], .pow(.sub(x, x0), .sub(n, .int(index))))), x, .add(list(x, x0)))))
    addargs <- if (f2$type == "Add") .sympy_args(f2) else list(f2)
    for (arg in addargs) {
      powm <- sym_match(arg, .mul(list(s, .pow(x, k))))
      term <- .mul(list(coeff[[1L]], powm$s))
      if (!(powm$k$type %in% c("Symbol", "Dummy"))) term <- subs(term, n, .sub(n, .ode_n_offset(powm$k, n)))
      startind <- subs(powm$k, n, .int(index))
      if (.ode_true(startind)) {
        for (i in rev(seq_len(as.integer(.ode_num(startind))) - 1L)) {
          if (!.ode_true(subs(term, n, .int(i)))) sset(term, .int(i))
          else {
            sset(term, .int(i + 1L))
            break
          }
        }
      } else sset(term, S_Zero)
    }
  }
  teq <- S_Zero
  req <- .add(sk)
  if (any(vapply(sv, .ode_true, TRUE))) {
    maxval <- sv[[which.max(vapply(sv, .ode_num, 0))]]
    for (i in seq_along(sk)) {
      val <- sv[[i]]
      if (!identical(val$key, maxval$key))
        for (jj in seq_len(.ode_num(maxval) - .ode_num(val))) teq <- .add(list(teq, subs(sk[[i]], n, val)))
    }
  }
  fk <- list(); fv <- list()
  fset <- function(key, v) {
    j <- match(key$key, vapply(fk, .k, ""))
    if (is.na(j)) { fk[[length(fk) + 1L]] <<- key; fv[[length(fv) + 1L]] <<- v } else fv[[j]] <<- v
  }
  undef <- function(e) Filter(.is_undefined_application, atoms(e, "Function"))
  argnum <- function(f) .ode_num(f$args[[1L]])
  if (.ode_true(teq)) {
    fargs <- undef(teq)
    if (length(fargs) == 1L) fset(fargs[[1L]], S_Zero)
    else {
      maxf <- fargs[[which.max(vapply(fargs, argnum, 0))]]
      sol <- .sp_solve(teq, list(maxf), list())
      fset(maxf, sol[[1L]])
    }
  }
  fargs <- undef(req)
  off <- vapply(fargs, function(f) .ode_num(.ode_n_offset(f$args[[1L]], n)), 0)
  maxf <- fargs[[which.max(off)]]
  minf <- fargs[[which.min(off)]]
  startiter <- if (minf$args[[1L]]$type %in% c("Symbol", "Dummy")) 0L else -as.integer(.ode_num(.ode_n_offset(minf$args[[1L]], n)))
  lhs <- maxf
  rhs <- .sp_solve(req, list(maxf), list())[[1L]]
  tcounter <- sum(vapply(fv, .ode_true, TRUE))
  if (tcounter < terms - 3L) for (it in seq(tcounter, terms - 4L)) {
    check <- subs(rhs, n, .int(startiter))
    nlhs <- subs(lhs, n, .int(startiter))
    nrhs <- .ode_subs_unordered(check, fk, fv)
    fset(nlhs, nrhs)
    startiter <- startiter + 1L
  }
  series <- .add(list(C0, .mul(list(C1, .sub(x, x0)))))
  for (i in seq_along(fk)) {
    if (.ode_true(fv[[i]])) {
      fact <- fk[[i]]$args[[1L]]
      series <- .add(list(series, .mul(list(subs(fv[[i]], list(list(recurr(S_Zero), C0), list(recurr(S_One), C1))),
                                             .pow(.sub(x, x0), fact)))))
    }
  }
  series <- .add(list(collect(expand_mul(series), list(C0, C1)), O(.pow(x, .int(terms)))))
  Eq(func, series)
}

# ode.ode_2nd_power_series_regular
.ode_2nd_power_series_regular <- function(eq, func, order, match) {
  x <- func$args[[1L]]
  C <- .ode_get_numbered_constants(eq, num = 2L)
  C0 <- C[[1L]]; C1 <- C[[2L]]
  m <- Dummy("m")
  x0 <- match$x0
  terms <- as.integer(match$terms)
  p <- match$p; q <- match$q
  indicial <- list()
  for (term in list(p, q)) {
    if (!.has1(term, x)) indicial[[length(indicial) + 1L]] <- term
    else {
      term <- series(term, x = x, x0 = x0, n = 1L)
      if (term$type == "Order") indicial[[length(indicial) + 1L]] <- S_Zero
      else for (arg in .sympy_args(term)) if (!.has1(arg, x)) {
        indicial[[length(indicial) + 1L]] <- arg
        break
      }
    }
  }
  p0 <- indicial[[1L]]; q0 <- indicial[[2L]]
  sollist <- .sp_solve(.add(list(.mul(list(m, .sub(m, S_One))), .mul(list(m, p0)), q0)), list(m), list())
  if (!length(sollist) || !all(vapply(sollist, function(z) isTRUE(.sp_fact(z, "real")), TRUE))) return(NULL)
  serdict2 <- list()
  if (length(sollist) == 1L) {
    m1 <- m2 <- sollist[[1L]]
    if (.ode_num(.sub(.sub(.int(terms), m1), S_One)) <= 0) return(Eq(func, O(.int(terms))))
    serdict1 <- .ode_frobenius(.sub(.sub(.int(terms), m1), S_One), m1, p0, q0, p, q, x0, x, C0)
  } else {
    m1 <- sollist[[1L]]; m2 <- sollist[[2L]]
    if (.ode_num(m1) < .ode_num(m2)) { tmp <- m1; m1 <- m2; m2 <- tmp }
    serdict1 <- .ode_frobenius(.sub(.sub(.int(terms), m1), S_One), m1, p0, q0, p, q, x0, x, C0)
    n2 <- .sub(.sub(.int(terms), m2), S_One)
    serdict2 <- if (!isTRUE(.sp_fact(.sub(m1, m2), "integer"))) .ode_frobenius(n2, m2, p0, q0, p, q, x0, x, C1)
                else .ode_frobenius(n2, m2, p0, q0, p, q, x0, x, C1, check = m1)
  }
  if (!is.list(serdict1) || !length(serdict1)) return(NULL)
  power_of <- function(nm) as.integer(substring(nm, 2L))
  fs1 <- C0
  for (nm in names(serdict1)) fs1 <- .add(list(fs1, .mul(list(serdict1[[nm]], .pow(.sub(x, x0), .int(power_of(nm)))))))
  fs1 <- .mul(list(.pow(.sub(x, x0), m1), fs1))
  fs2 <- S_Zero
  if (is.list(serdict2) && length(serdict2)) {
    for (nm in names(serdict2)) fs2 <- .add(list(fs2, .mul(list(serdict2[[nm]], .pow(.sub(x, x0), .int(power_of(nm)))))))
    fs2 <- .add(list(fs2, C1))
    fs2 <- .mul(list(.pow(.sub(x, x0), m2), fs2))
  }
  Eq(func, .add(list(collect(.add(list(fs1, fs2)), list(C0, C1)), O(.pow(x, .int(terms))))))
}

# ode._frobenius: the coefficients C1, ..., Cn of a Frobenius series, as a
# list named by them, or FALSE
.ode_frobenius <- function(n, m, p0, q0, p, q, x0, x, c, check = NULL) {
  n <- as.integer(trunc(.ode_num(n)))
  m2 <- check
  d <- Dummy("d")
  serlist <- list()
  for (ser in list(p, q)) {
    if (isTRUE(is_polynomial(ser, x)) && degree(Poly(ser, x)) <= n) {
      if (.ode_true(x0)) ser <- subs(ser, x, .add(list(x, x0)))
      e <- expand(ser)
    } else e <- expand(removeO(series(ser, x = x, x0 = x0, n = n + 1L)))
    serlist[[length(serlist) + 1L]] <- lapply(0:n, function(i) coeff(e, x, i))
  }
  ps <- serlist[[1L]]; qs <- serlist[[2L]]
  indicial <- .add(list(.mul(list(d, .sub(d, S_One))), .mul(list(d, p0)), q0))
  frob <- list()
  if (n >= 1L) for (i in seq_len(n)) {
    num <- .mul(list(c, .add(list(.mul(list(m, ps[[i + 1L]])), qs[[i + 1L]]))))
    if (i > 1L) for (j in seq_len(i - 1L))
      num <- .add(list(num, .mul(list(frob[[paste0("C", j)]], .add(list(.mul(list(.add(list(m, .int(j))), ps[[i - j + 1L]])), qs[[i - j + 1L]]))))))
    if (!is.null(m2) && identical(.int(i)$key, .sub(m2, m)$key)) {
      if (.ode_true(num)) return(FALSE)
      frob[[paste0("C", i)]] <- S_Zero
    } else frob[[paste0("C", i)]] <- .neg(.div(num, subs(indicial, d, .add(list(m, .int(i))))))
  }
  frob
}

# ---------------------------------------------------------------------------
# checkinfsol

#' Check infinitesimals of an ODE
#'
#' `checkinfsol(eq, infinitesimals, func)` checks whether each set of
#' infinitesimals (a named list with the entries `xi` and `eta`, or a list
#' of such lists) satisfies the determining PDE of the first order ODE `eq`,
#' as SymPy's `checkinfsol()`.
#'
#' @param eq The ODE.
#' @param infinitesimals The infinitesimals, as returned by
#'   [infinitesimals()].
#' @param func The unknown function; detected when `NULL`.
#' @param order The order of the ODE; computed when `NULL`.
#' @return A list of pairs `(TRUE, 0)` or `(FALSE, residual)`.
#' @examples
#' x <- Symbol("x"); f <- Function("f")
#' eq <- Derivative(f(x), x) - x^2
#' checkinfsol(eq, infinitesimals(eq, f(x)), f(x))
#' @export
checkinfsol <- function(eq, infinitesimals, func = NULL, order = NULL) {
  eq <- as_sym(eq)
  if (eq$type == "Equality") eq <- .sub(eq$args[[1L]], eq$args[[2L]])
  if (is.null(func)) {
    pp <- .ode_preprocess(eq)
    eq <- pp[[1L]]
    func <- pp[[2L]]
  }
  func <- as_sym(func)
  variables <- func$args
  if (length(variables) != 1L) .ode_stop("ValueError", "ODE's have only one independent variable")
  x <- variables[[1L]]
  if (is.null(order)) order <- .ode_order(eq, func)
  if (order != 1L) .ode_nie("Lie groups solver has been implemented only for first order differential equations")
  df <- .ode_d(func, x)
  a <- Wild("a", exclude = list(df))
  b <- Wild("b", exclude = list(df))
  m <- sym_match(collect(expand(eq), df), a * df + b)
  if (!is.null(m)) h <- .neg(simplify(.div(m$b %||% S_Zero, m$a %||% S_Zero)))
  else {
    sol <- tryCatch(.sp_solve(eq, list(df), list()), NotImplementedError = function(e)
      .ode_nie("Infinitesimals for the first order ODE could not be found"))
    h <- sol[[1L]]
  }
  y <- Dummy("y")
  h <- subs(h, func, y)
  xi <- fcall("xi", list(x, y))
  eta <- fcall("eta", list(x, y))
  dxi <- fcall("xi", list(x, func))
  deta <- fcall("eta", list(x, func))
  pde <- .add(list(.diff1(eta, x), .mul(list(.sub(.diff1(eta, y), .diff1(xi, x)), h)),
                   .neg(.mul(list(.diff1(xi, y), .pow(h, S_Two)))), .neg(.mul(list(xi, .diff1(h, x)))),
                   .neg(.mul(list(eta, .diff1(h, y))))))
  infs <- .ode_inf_list(infinitesimals)
  lapply(infs, function(s) {
    tsol <- list(list(xi, subs(as_sym(.ode_inf_get(s, dxi, "xi")), func, y)),
                 list(eta, subs(as_sym(.ode_inf_get(s, deta, "eta")), func, y)))
    r <- simplify(doit(subs(pde, tsol)))
    if (.ode_true(r)) .pytuple(list(FALSE, subs(r, y, func))) else .pytuple(list(TRUE, S_Zero))
  })
}

# infinitesimals come as dictionaries keyed by xi(x, f(x)) and eta(x, f(x)),
# or as named lists with the entries xi and eta
.ode_inf_list <- function(infs) {
  if (.ode_inf_one(infs)) return(list(.ode_inf_dict(infs)))
  lapply(infs, .ode_inf_dict)
}
# one dictionary: a pydict, a named list, or a list of (key, value) pairs
.ode_inf_one <- function(s) {
  inherits(s, "pydict") || (is.list(s) && !is.null(names(s)) && all(c("xi", "eta") %in% names(s))) ||
    (is.list(s) && length(s) && all(vapply(s, function(p) is.list(p) && !is_sym(p) && length(p) == 2L && is_sym(p[[1L]]), TRUE)))
}
.ode_inf_dict <- function(s) {
  if (inherits(s, "pydict") || !is.null(names(s))) return(s)
  .pydict(lapply(s, `[[`, 1L), lapply(s, function(p) as_sym(p[[2L]])))
}
.ode_inf_get <- function(s, key, nm) {
  if (inherits(s, "pydict")) {
    for (i in seq_along(s$keys)) if (identical(as_sym(s$keys[[i]])$key, key$key)) return(s$values[[i]])
    for (i in seq_along(s$keys)) if (is_sym(s$keys[[i]]) && s$keys[[i]]$type == "Function" && s$keys[[i]]$name == nm) return(s$values[[i]])
    stop(.ode_cond("KeyError", .str(key)))
  }
  s[[nm]]
}
