# Ordinary differential equations and recurrences, after sympy/solvers/ode
# and sympy/solvers/recurr.py.

# Integration constants C1, C2, ... start after any C<k> already in the
# equation, so they never capture a parameter of the input.
.ode_env <- new.env()
.ode_env$shift <- 0L
.C <- function(i) Symbol(paste0("C", i + .ode_env$shift))
.is_gen_const <- function(s) grepl("^C[0-9]+$", s$name) && as.integer(substring(s$name, 2L)) > .ode_env$shift

# Replace f(x) and its derivatives by plain symbols. Returns the rewritten
# equation and the list of symbols y0, y1, ... up to the order.
.ode_symbols <- function(eq, fx, x) {
  order <- ode_order(eq, fx)
  ys <- lapply(0:order, function(k) Dummy(paste0("y", k)))
  walk <- function(z) {
    if (z$type == "Derivative" && identical(z$args[[1L]]$key, fx$key)) {
      n <- as.integer(as.numeric(z$args[[2L]]$args[[2L]]))
      return(ys[[n + 1L]])
    }
    if (identical(z$key, fx$key)) return(ys[[1L]])
    if (!length(z$args)) return(z)
    sym_rebuild(z, lapply(z$args, walk))
  }
  list(eq = walk(eq), ys = ys, order = order)
}

# sympy.solvers.ode.ode.allhints
.ODE_ALLHINTS <- c("factorable", "nth_algebraic", "separable", "1st_exact", "1st_linear",
                   "Bernoulli", "1st_rational_riccati", "Riccati_special_minus2",
                   "1st_homogeneous_coeff_best", "1st_homogeneous_coeff_subs_indep_div_dep",
                   "1st_homogeneous_coeff_subs_dep_div_indep", "almost_linear",
                   "linear_coefficients", "separable_reduced", "1st_power_series", "lie_group",
                   "nth_linear_constant_coeff_homogeneous", "nth_linear_euler_eq_homogeneous",
                   "nth_linear_constant_coeff_undetermined_coefficients",
                   "nth_linear_euler_eq_nonhomogeneous_undetermined_coefficients",
                   "nth_linear_constant_coeff_variation_of_parameters",
                   "nth_linear_euler_eq_nonhomogeneous_variation_of_parameters", "Liouville",
                   "2nd_linear_airy", "2nd_linear_bessel", "2nd_linear_bessel_transform",
                   "2nd_hypergeometric", "2nd_hypergeometric_Integral", "nth_order_reducible",
                   "2nd_power_series_ordinary", "2nd_power_series_regular",
                   "nth_algebraic_Integral", "separable_Integral", "1st_exact_Integral",
                   "1st_linear_Integral", "Bernoulli_Integral",
                   "1st_homogeneous_coeff_subs_indep_div_dep_Integral",
                   "1st_homogeneous_coeff_subs_dep_div_indep_Integral", "almost_linear_Integral",
                   "linear_coefficients_Integral", "separable_reduced_Integral",
                   "nth_linear_constant_coeff_variation_of_parameters_Integral",
                   "nth_linear_euler_eq_nonhomogeneous_variation_of_parameters_Integral",
                   "Liouville_Integral", "2nd_nonlinear_autonomous_conserved",
                   "2nd_nonlinear_autonomous_conserved_Integral")

# _preprocess: the one undefined function applied inside the derivatives
.ode_detect_func <- function(e) {
  fs <- list()
  walk <- function(z) {
    if (z$type == "Derivative") for (f in Filter(.is_undefined_application, atoms(z, "Function"))) fs[[f$key]] <<- f
    for (a in z$args) walk(a)
  }
  walk(e)
  if (length(fs) != 1L)
    stop(.sympy_condition("ValueError", paste0("The function cannot be automatically detected for ", .sp_py_str(e), ".")))
  fs[[1L]]
}

# the checks _desolve makes before solving
.ode_check <- function(e, fx, hint) {
  if (length(fx$args) != 1L)
    stop(.sympy_condition("ValueError", paste0("dsolve() and classify_ode() only work with functions of one variable, not ", .sp_py_str(fx))))
  if (!identical(hint, "default") && !(hint %in% c(.ODE_ALLHINTS, "all", "all_Integral", "best")))
    stop(.sympy_condition("ValueError", paste0("Hint not recognized: ", hint)))
  if (ode_order(e, fx) == 0L)
    stop(.sympy_condition("ValueError", paste0(.sp_py_str(e), " is not a solvable differential equation in ", .sp_py_str(fx))))
}

# the checks classify_ode makes of ics, given as list(list(key, value), ...)
.ode_check_ics <- function(ics, fx) {
  x <- fx$args[[1L]]
  is_f <- function(z) z$type == "Function" && identical(z$name, fx$name) && length(z$args) == 1L
  for (p in ics) {
    key <- as_sym(p[[1L]])
    val <- as_sym(p[[2L]])
    if (key$type %in% c("Subs", "Derivative")) {
      if (key$type == "Subs") {
        vp <- .subs_parts(key)
        deriv <- key$args[[1L]]; old <- vp[[1L]][[1L]]; new <- vp[[2L]][[1L]]
      } else {
        deriv <- key; old <- x; new <- key$args[[2L]]$args[[1L]]
      }
      ok <- deriv$type == "Derivative" && is_f(deriv$args[[1L]]) && identical(old$key, x$key) && !.has1(new, x) &&
        !.has1(val, x)
      if (ok) {
        dv <- vapply(deriv$args[-1L], function(t) if (t$type == "Tuple") t$args[[1L]]$key else t$key, "")
        ok <- all(dv == dv[1L])
      }
      if (!ok) stop(.sympy_condition("ValueError", "Invalid boundary conditions for Derivatives"))
    } else if (.is_undefined_application(key)) {
      if (!(is_f(key) && !.has1(key$args[[1L]], x) && !.has1(val, x)))
        stop(.sympy_condition("ValueError", "Invalid boundary conditions for Function"))
    } else {
      stop(.sympy_condition("ValueError", paste0("Enter boundary conditions of the form ics={f(point): value, ",
                                                 "f(x).diff(x, order).subs(x, point): value}")))
    }
  }
}

#' Solve ordinary differential equations
#'
#' `dsolve(eq, fx)` solves an ordinary differential equation for the unknown
#' function application `fx`, such as `f(x)`, returning an [Eq()]. Supported
#' classes: separable and linear first order equations, Bernoulli and exact
#' equations, linear equations with constant coefficients of any order (with
#' undetermined coefficients for polynomial, exponential and sinusoidal
#' forcing) and Cauchy-Euler equations. Arbitrary constants are `C1`, `C2`,
#' ... `ics` applies initial conditions, given as a list of
#' `list(point, order, value)`. `ode_order()` returns the order and
#' `classify_ode()` the matching classes.
#'
#' @param eq The equation or expression equal to zero.
#' @param fx The unknown, for example `f(x)`; found from the derivatives in
#'   `eq` when omitted.
#' @param ics Optional initial conditions: for `dsolve()` a list of
#'   `list(point, order, value)`; for `classify_ode()` a list of
#'   `list(key, value)` pairs, where a key is `f(point)` or a derivative of
#'   `f` at a point, checked as SymPy checks them.
#' @param hint A SymPy hint name; any recognized hint is accepted, and the
#'   solution comes from the classes above.
#' @param dict For `classify_ode()`, return the order, the default hint and
#'   the ordered hints as a list.
#' @param ... Further SymPy options, ignored.
#' @return An [Eq()], or a list of them for several solutions.
#' @examples
#' x <- Symbol("x"); f <- Function("f")
#' dsolve(Derivative(f(x), x) - f(x), f(x))
#' dsolve(Derivative(f(x), x, 2) + f(x), f(x))
#' @export
dsolve <- function(eq, fx = NULL, ics = NULL, hint = "default", ...) {
  e <- .eq_to_expr(eq)
  # dsolve(eq, f(x), "hint") as SymPy's positional form
  if (is.character(ics)) { hint <- ics; ics <- NULL }
  fx <- if (is.null(fx)) .ode_detect_func(e) else as_sym(fx)
  .ode_check(e, fx, hint)
  x <- fx$args[[1L]]
  e <- doit(e)
  used <- Filter(function(s) grepl("^C[0-9]+$", s$name), free_symbols(e))
  old <- .ode_env$shift
  .ode_env$shift <- if (length(used)) max(vapply(used, function(s) as.integer(substring(s$name, 2L)), 0L)) else 0L
  on.exit(.ode_env$shift <- old)
  sv <- .ode_symbols(e, fx, x)
  sol <- .ode_dispatch(sv, fx, x)
  if (is.null(sol)) stop("dsolve() has no method for this equation")
  if (!is.null(ics)) sol <- .apply_ics(sol, fx, x, ics)
  sol
}

#' @rdname dsolve
#' @export
ode_order <- function(eq, fx) {
  fx <- as_sym(fx)
  # sympy.solvers.deutils.ode_order: a derivative of another expression adds
  # its own order to the highest order among that expression's arguments
  rec <- function(z) {
    if (!.has1(z, fx)) return(0L)
    if (z$type == "Derivative") {
      nv <- sum(vapply(z$args[-1L], function(v) as.integer(as.numeric(v$args[[2L]])), 0L))
      if (identical(z$args[[1L]]$key, fx$key)) return(nv)
      a <- z$args[[1L]]$args
      return(nv + if (length(a)) max(vapply(a, rec, 0L)) else 0L)
    }
    if (length(z$args)) max(vapply(z$args, rec, 0L)) else 0L
  }
  rec(as_sym(eq))
}

.linear_ode_coeffs <- function(sv, x) {
  e <- expand(sv$eq)
  ys <- sv$ys
  cs <- lapply(ys, function(y) .linear_coeff(e, y))
  rest <- expand(.sub(e, .add(Map(function(c, y) .mul(list(c, y)), cs, ys))))
  dep <- function(z) any(vapply(ys, function(y) .has1(z, y), TRUE))
  if (any(vapply(cs, dep, TRUE)) || dep(rest)) return(NULL)
  list(c = cs, r = .neg(rest))
}

#' @rdname dsolve
#' @export
classify_ode <- function(eq, fx = NULL, dict = FALSE, ics = NULL, ...) {
  if (!is.null(fx) && length(as_sym(fx)$args) != 1L)
    stop(.sympy_condition("ValueError", paste0("dsolve() and classify_ode() only work with functions of one variable, not ",
                                               .sp_py_str(as_sym(fx)))))
  e <- doit(.eq_to_expr(eq))
  fx <- if (is.null(fx)) .ode_detect_func(e) else as_sym(fx)
  if (!is.null(ics)) .ode_check_ics(ics, fx)
  if (dict) {
    hs <- classify_ode(e, fx)
    return(list(order = ode_order(e, fx), default = if (length(hs)) hs[[1L]], ordered_hints = hs))
  }
  x <- fx$args[[1L]]
  sv <- .ode_symbols(e, fx, x)
  out <- character()
  lin <- .linear_ode_coeffs(sv, x)
  if (!is.null(lin)) {
    const <- all(vapply(lin$c, function(c) !.has1(c, x), TRUE))
    if (sv$order == 1L) out <- c(out, "1st_linear")
    if (const) out <- c(out, if (is_zero_num(lin$r)) "nth_linear_constant_coeff_homogeneous"
                             else "nth_linear_constant_coeff_undetermined_coefficients")
    if (.is_euler(lin, x)) out <- c(out, "nth_linear_euler_eq_homogeneous")
  }
  if (sv$order == 1L && !is.null(.separate(sv, x))) out <- c(out, "separable")
  out
}

.ode_dispatch <- function(sv, fx, x) {
  lin <- .linear_ode_coeffs(sv, x)
  if (!is.null(lin)) {
    const <- all(vapply(lin$c, function(c) !.has1(c, x), TRUE))
    if (const) {
      r <- .ode_const_coeff(lin, fx, x)
      if (!is.null(r)) return(r)
    }
    if (.is_euler(lin, x)) {
      r <- .ode_euler(lin, fx, x)
      if (!is.null(r)) return(r)
    }
    if (sv$order == 1L) {
      r <- .ode_first_linear(lin, fx, x)
      if (!is.null(r)) return(r)
    }
  }
  if (sv$order == 1L) {
    r <- .ode_bernoulli(sv, fx, x)
    if (!is.null(r)) return(r)
    r <- .ode_separable(sv, fx, x)
    if (!is.null(r)) return(r)
    r <- .ode_exact(sv, fx, x)
    if (!is.null(r)) return(r)
  }
  NULL
}

# a(x) y' + b(x) y = r(x)
.ode_first_linear <- function(lin, fx, x) {
  a <- lin$c[[2L]]; b <- lin$c[[1L]]; r <- lin$r
  P <- simplify(.div(b, a))
  Q <- simplify(.div(r, a))
  IP <- integrate(P, x)
  if (IP$type == "Integral") return(NULL)
  mu <- simplify(fcall("exp", list(IP)))
  I2 <- integrate(simplify(.mul(list(mu, Q))), x)
  if (I2$type == "Integral") return(NULL)
  Eq(fx, simplify(.div(.add(list(.C(1), I2)), mu)), evaluate = FALSE)
}

# y' = g(x) h(y): returns list(g, h) or NULL.
.separate <- function(sv, x) {
  y <- sv$ys[[1L]]; yp <- sv$ys[[2L]]
  e <- expand(sv$eq)
  a <- .linear_coeff(e, yp)
  if (.has1(a, yp)) return(NULL)
  rest <- expand(.sub(e, .mul(list(a, yp))))
  if (.has1(rest, yp)) return(NULL)
  rhs <- factor_sym(simplify(.neg(.div(rest, a))))
  facs <- if (rhs$type == "Mul") rhs$args else list(rhs)
  g <- list(); h <- list()
  for (f in facs) {
    hx <- .has1(f, x); hy <- .has1(f, y)
    if (hx && hy) return(NULL)
    if (hy) h[[length(h) + 1L]] <- f else g[[length(g) + 1L]] <- f
  }
  list(g = .mul(g), h = .mul(h))
}

.ode_separable <- function(sv, fx, x) {
  s <- .separate(sv, x)
  if (is.null(s)) return(NULL)
  y <- sv$ys[[1L]]
  H <- integrate(.inv(s$h), y)
  G <- integrate(s$g, x)
  if (H$type == "Integral" || G$type == "Integral") return(NULL)
  eq <- .sub(H, .add(list(G, .C(1))))
  sols <- tryCatch(.solve_univariate(eq, y, check = FALSE), error = function(err) NULL)
  if (is.null(sols) || !length(sols)) return(Eq(subs(H, y, fx), .add(list(G, .C(1))), evaluate = FALSE))
  out <- lapply(sols, function(v) Eq(fx, simplify_if_cheap(v), evaluate = FALSE))
  if (length(out) == 1L) out[[1L]] else out
}

# y' + P y = Q y^n
.ode_bernoulli <- function(sv, fx, x) {
  y <- sv$ys[[1L]]; yp <- sv$ys[[2L]]
  e <- expand(sv$eq)
  a <- .linear_coeff(e, yp)
  if (.has1(a, y) || .has1(a, yp)) return(NULL)
  rest <- expand(.div(.sub(e, .mul(list(a, yp))), a))
  terms <- if (rest$type == "Add") rest$args else list(rest)
  P <- list(); Q <- list(); n <- NULL
  for (t in terms) {
    d <- .degree_in(t, y)
    if (is.na(d) || d == 0L) return(NULL)
    if (d == 1L) {
      P[[length(P) + 1L]] <- .div(t, y)
    } else {
      if (!is.null(n) && n != d) return(NULL)
      n <- d
      Q[[length(Q) + 1L]] <- .neg(.div(t, .pow(y, .int(d))))
    }
  }
  if (is.null(n) || !length(P)) return(NULL)
  P <- .add(P); Q <- .add(Q)
  # v = y^(1-n): v' + (1-n) P v = (1-n) Q
  k <- .int(1L - n)
  lin <- list(c = list(.mul(list(k, P)), S_One), r = .mul(list(k, Q)))
  v <- Dummy("v")
  r <- .ode_first_linear(lin, v, x)
  if (is.null(r)) return(NULL)
  sol <- .pow(r$args[[2L]], .inv(k))
  Eq(fx, simplify_if_cheap(sol), evaluate = FALSE)
}

# M(x, y) + N(x, y) y' = 0 with dM/dy = dN/dx
.ode_exact <- function(sv, fx, x) {
  y <- sv$ys[[1L]]; yp <- sv$ys[[2L]]
  e <- expand(sv$eq)
  N <- .linear_coeff(e, yp)
  if (.has1(N, yp)) return(NULL)
  M <- expand(.sub(e, .mul(list(N, yp))))
  if (!is_zero_num(simplify(.sub(.diff1(M, y), .diff1(N, x))))) return(NULL)
  F <- integrate(M, x)
  if (F$type == "Integral") return(NULL)
  rest <- simplify(.sub(N, .diff1(F, y)))
  if (.has1(rest, x)) return(NULL)
  G <- integrate(rest, y)
  if (G$type == "Integral") return(NULL)
  Eq(subs(.add(list(F, G)), y, fx), .C(1), evaluate = FALSE)
}

# Roots of the characteristic polynomial, returned as real exponent and
# oscillation pairs together with the basis functions.
.char_basis <- function(cs, x) {
  m <- Dummy("m")
  n <- length(cs) - 1L
  poly <- .add(lapply(seq_along(cs), function(k) .mul(list(cs[[k]], .pow(m, .int(k - 1L))))))
  rs <- .roots_legacy(poly, m)  # SymPy uses rootof(): real roots sorted first
  if (sum(attr(rs, "multiplicity")) != n) return(NULL)
  list(roots = rs, mult = attr(rs, "multiplicity"))
}

.ode_homog_basis <- function(rs, mult, x, var = x, real_coeffs = TRUE) {
  basis <- list()
  done <- character()
  for (i in seq_along(rs)) {
    r <- rs[[i]]
    if (r$key %in% done) next
    if (!real_coeffs) {
      for (k in 0:(mult[i] - 1L)) basis[[length(basis) + 1L]] <- .mul(list(.pow(var, .int(k)), fcall("exp", list(.mul(list(r, x))))))
      done <- c(done, r$key)
      next
    }
    ri <- .re_im(r)
    a <- ri[[1L]]; b <- ri[[2L]]
    for (k in 0:(mult[i] - 1L)) {
      pw <- .pow(var, .int(k))
      if (is_zero_num(b)) {
        basis[[length(basis) + 1L]] <- .mul(list(pw, fcall("exp", list(.mul(list(a, x))))))
      } else {
        b <- fcall("Abs", list(b))
        basis[[length(basis) + 1L]] <- .mul(list(pw, fcall("exp", list(.mul(list(a, x)))), fcall("sin", list(.mul(list(b, x))))))
        basis[[length(basis) + 1L]] <- .mul(list(pw, fcall("exp", list(.mul(list(a, x)))), fcall("cos", list(.mul(list(b, x))))))
      }
    }
    done <- c(done, r$key)
    if (!is_zero_num(ri[[2L]])) done <- c(done, fcall("conjugate", list(r))$key)
  }
  basis
}

.ode_const_coeff <- function(lin, fx, x) {
  ch <- .char_basis(lin$c, x)
  if (is.null(ch)) return(NULL)
  real_coeffs <- all(vapply(lin$c, function(c) isTRUE(ask_fact(c, "real")), TRUE))
  basis <- .ode_homog_basis(ch$roots, ch$mult, x, real_coeffs = real_coeffs)
  yh <- .add(lapply(seq_along(basis), function(i) .mul(list(.C(i), basis[[i]]))))
  yp <- S_Zero
  if (!is_zero_num(lin$r)) {
    yp <- .undetermined(lin, x, ch)
    if (is.null(yp)) return(NULL)
  }
  Eq(fx, .add(list(yh, yp)), evaluate = FALSE)
}

# Particular solution by undetermined coefficients.
.undetermined <- function(lin, x, ch) {
  r <- expand(lin$r)
  terms <- if (r$type == "Add") r$args else list(r)
  groups <- list()
  for (t in terms) {
    facs <- if (t$type == "Mul") t$args else list(t)
    ae <- S_Zero; bw <- S_Zero; deg <- 0L
    for (f in facs) {
      if (!.has1(f, x)) next
      d <- .degree_in(f, x)
      if (!is.na(d)) {
        deg <- deg + d
        next
      }
      if (f$type == "Function" && f$name == "exp") {
        cl <- tryCatch(.coeff_list(f$args[[1L]], x), error = function(e) NULL)
        if (is.null(cl) || length(cl) != 2L || !is_zero_num(cl[[2L]])) return(NULL)
        ae <- cl[[1L]]
        next
      }
      if (f$type == "Function" && f$name %in% c("sin", "cos")) {
        cl <- tryCatch(.coeff_list(f$args[[1L]], x), error = function(e) NULL)
        if (is.null(cl) || length(cl) != 2L || !is_zero_num(cl[[2L]])) return(NULL)
        bw <- fcall("Abs", list(cl[[1L]]))
        next
      }
      return(NULL)
    }
    key <- paste(ae$key, bw$key)
    if (is.null(groups[[key]])) groups[[key]] <- list(a = ae, b = bw, deg = deg)
    else groups[[key]]$deg <- max(groups[[key]]$deg, deg)
  }
  unknowns <- list()
  trial <- list()
  k <- 0L
  fresh <- function() {
    k <<- k + 1L
    u <- Dummy(paste0("A", k))
    unknowns[[length(unknowns) + 1L]] <<- u
    u
  }
  for (g in groups) {
    root <- expand(.add(list(g$a, .mul(list(S_I, g$b)))))
    s <- 0L
    for (i in seq_along(ch$roots)) if (is_zero_num(simplify(.sub(ch$roots[[i]], root)))) s <- ch$mult[i]
    ex <- fcall("exp", list(.mul(list(g$a, x))))
    for (d in 0:g$deg) {
      pw <- .pow(x, .int(d + s))
      if (is_zero_num(g$b)) {
        trial[[length(trial) + 1L]] <- .mul(list(fresh(), pw, ex))
      } else {
        trial[[length(trial) + 1L]] <- .mul(list(fresh(), pw, ex, fcall("cos", list(.mul(list(g$b, x))))))
        trial[[length(trial) + 1L]] <- .mul(list(fresh(), pw, ex, fcall("sin", list(.mul(list(g$b, x))))))
      }
    }
  }
  Y <- .add(trial)
  L <- list()
  cur <- Y
  for (i in seq_along(lin$c)) {
    L[[i]] <- .mul(list(lin$c[[i]], cur))
    cur <- .diff1(cur, x)
  }
  res <- expand(.sub(.add(L), lin$r))
  # coefficients of independent functions of x must vanish
  gens <- Filter(function(g) .has1(g, x), .find_gens(list(res)))
  gens <- Filter(function(g) !any(vapply(unknowns, function(u) identical(u$key, g$key), TRUE)), gens)
  dmap <- lapply(gens, function(g) Dummy("g"))
  rs <- xreplace(res, Map(list, gens, dmap))
  terms <- if (rs$type == "Add") rs$args else list(rs)
  grp <- list()
  for (t in terms) {
    mono <- S_One
    facs <- if (t$type == "Mul") t$args else list(t)
    for (f in facs) {
      b <- .base_exp(f)[[1L]]
      if (any(vapply(dmap, function(d) identical(d$key, b$key), TRUE))) mono <- .mul(list(mono, f))
    }
    grp[[mono$key]] <- c(grp[[mono$key]], list(.div(t, mono)))
  }
  eqs <- lapply(grp, .add)
  sol <- .linsolve_exprs(unname(eqs), unknowns)
  if (is.null(sol)) return(NULL)
  xreplace(Y, sol)
}

.is_euler <- function(lin, x) {
  if (!is_zero_num(lin$r)) return(FALSE)
  for (k in seq_along(lin$c)) {
    q <- simplify(.div(lin$c[[k]], .pow(x, .int(k - 1L))))
    if (.has1(q, x)) return(FALSE)
  }
  length(lin$c) > 2L
}

# a_n x^n y^(n) + ... + a_0 y = 0
.ode_euler <- function(lin, fx, x) {
  n <- length(lin$c) - 1L
  m <- Dummy("m")
  terms <- list()
  for (k in 0:n) {
    a <- simplify(.div(lin$c[[k + 1L]], .pow(x, .int(k))))
    fall <- if (k == 0L) S_One else .mul(lapply(0:(k - 1L), function(j) .sub(m, .int(j))))
    terms[[length(terms) + 1L]] <- .mul(list(a, fall))
  }
  rs <- .roots_legacy(expand(.add(terms)), m)
  if (sum(attr(rs, "multiplicity")) != n) return(NULL)
  mult <- attr(rs, "multiplicity")
  basis <- list()
  lx <- fcall("log", list(x))
  for (i in seq_along(rs)) {
    ri <- .re_im(rs[[i]])
    for (k in 0:(mult[i] - 1L)) {
      pw <- .pow(lx, .int(k))
      if (is_zero_num(ri[[2L]])) {
        basis[[length(basis) + 1L]] <- .mul(list(pw, .pow(x, ri[[1L]])))
      } else if (isTRUE(ask_fact(ri[[2L]], "positive"))) {
        basis[[length(basis) + 1L]] <- .mul(list(pw, .pow(x, ri[[1L]]), fcall("sin", list(.mul(list(ri[[2L]], lx))))))
        basis[[length(basis) + 1L]] <- .mul(list(pw, .pow(x, ri[[1L]]), fcall("cos", list(.mul(list(ri[[2L]], lx))))))
      }
    }
  }
  Eq(fx, .add(lapply(seq_along(basis), function(i) .mul(list(.C(i), basis[[i]])))), evaluate = FALSE)
}

# Solve linear equations for the constants that they determine; constants
# left free stay symbolic.
.solve_constants <- function(eqs, cs) {
  eqs <- lapply(eqs, expand)
  lin <- all(vapply(eqs, function(e) .is_linear_in(e, cs), TRUE))
  if (!lin) {
    # nonlinear in the constants: general solver, first verified branch
    sols <- tryCatch(solve(eqs, cs, dict = TRUE), error = function(e) list())
    for (sl in sols) {
      if (!is.list(sl) || is.null(names(sl))) next
      res <- vapply(eqs, function(e) is_zero_num(simplify(subs(e, sl))), TRUE)
      if (all(res)) return(lapply(Filter(function(c) !is.null(sl[[c$name]]), cs), function(c) list(c, sl[[c$name]])))
    }
    stop("initial conditions are inconsistent")
  }
  s <- .solve_linear_system(eqs, cs)
  if (is.null(s)) stop("initial conditions are inconsistent")
  lapply(Filter(function(c) !is.null(s[[c$name]]), cs), function(c) list(c, s[[c$name]]))
}

.apply_ics <- function(sol, fx, x, ics) {
  if (is.list(sol) && !inherits(sol, "sym")) {
    for (cand in sol) {
      r <- tryCatch(.apply_ics(cand, fx, x, ics), error = function(e) NULL)
      if (!is.null(r)) return(r)
    }
    stop("no solution branch satisfies the initial conditions")
  }
  if (!identical(sol$args[[1L]]$key, fx$key)) return(.apply_ics_implicit(sol, fx, x, ics))
  rhs <- sol$args[[2L]]
  cs <- Filter(.is_gen_const, free_symbols(rhs))
  eqs <- lapply(ics, function(ic) {
    d <- rhs
    for (k in seq_len(ic[[2L]])) d <- .diff1(d, x)
    .sub(subs(d, x, as_sym(ic[[1L]])), as_sym(ic[[3L]]))
  })
  s <- .solve_constants(eqs, cs)
  Eq(fx, simplify_if_cheap(xreplace(rhs, s)), evaluate = FALSE)
}

# Implicit solution F(x, f(x)) = C1 with f(x0) = y0.
.apply_ics_implicit <- function(sol, fx, x, ics) {
  if (length(ics) != 1L || ics[[1L]][[2L]] != 0L) stop("implicit solutions accept one value condition")
  x0 <- as_sym(ics[[1L]][[1L]]); y0 <- as_sym(ics[[1L]][[3L]])
  F <- .sub(sol$args[[1L]], sol$args[[2L]])
  y <- Dummy("y")
  Fy <- subs(F, fx, y)
  cs <- Filter(.is_gen_const, free_symbols(Fy))
  rule <- .solve_constants(list(subs(subs(Fy, x, x0), y, y0)), cs)
  G <- xreplace(Fy, rule)
  ys <- tryCatch(.solve_univariate(G, y, check = FALSE), error = function(e) list())
  for (v in ys) {
    at <- tryCatch(simplify(.sub(subs(v, x, x0), y0)), error = function(e) NULL)
    if (!is.null(at) && is_zero_num(at)) return(Eq(fx, simplify_if_cheap(v), evaluate = FALSE))
  }
  Eq(subs(G, y, fx), S_Zero, evaluate = FALSE)
}

# The constant coefficient recurrence solver used by fps(); rsolve() is in
# solvers-recurr.R.
.rsolve_cc <- function(eq, yn, n, init = NULL) {
  e <- expand(.eq_to_expr(eq))
  yn <- as_sym(yn)
  n <- as_sym(n)
  name <- yn$name
  shifts <- list()
  walk <- function(z) {
    if (z$type == "Function" && z$name == name) {
      d <- expand(.sub(z$args[[1L]], n))
      if (d$type != "Integer") stop("rsolve() needs integer shifts of the index")
      shifts[[z$key]] <<- list(node = z, k = .int_val(d))
    }
    for (a in z$args) walk(a)
  }
  walk(e)
  shifts <- unname(shifts)
  if (!length(shifts)) stop("the unknown does not appear in the recurrence")
  ks <- vapply(shifts, `[[`, 0, "k")
  lo <- min(ks); hi <- max(ks)
  ds <- lapply(shifts, function(s) Dummy("s"))
  er <- xreplace(e, Map(function(s, d) list(s$node, d), shifts, ds))
  coefs <- rep(list(S_Zero), hi - lo + 1)
  for (i in seq_along(shifts)) coefs[[ks[i] - lo + 1]] <- .linear_coeff(er, ds[[i]])
  if (any(vapply(coefs, function(c) .has1(c, n), TRUE))) {
    rest <- expand(.sub(er, .add(Map(function(c, d) .mul(list(c, d)), lapply(seq_along(shifts), function(i) coefs[[ks[i] - lo + 1]]), ds))))
    return(.rsolve_general(coefs, lo, rest, n, init))
  }
  rest <- expand(.sub(er, .add(Map(function(c, d) .mul(list(c, d)), lapply(seq_along(shifts), function(i) coefs[[ks[i] - lo + 1]]), ds))))
  m <- Dummy("m")
  poly <- .add(lapply(seq_along(coefs), function(k) .mul(list(coefs[[k]], .pow(m, .int(k - 1L))))))
  rs <- .roots_legacy(poly, m)
  mult <- attr(rs, "multiplicity")
  if (sum(mult) != length(coefs) - 1L) stop("could not factor the characteristic polynomial")
  terms <- list()
  idx <- 0L
  cs <- list()
  for (i in seq_along(rs)) for (j in 0:(mult[i] - 1L)) {
    c <- Symbol(paste0("C", idx))
    cs[[length(cs) + 1L]] <- c
    terms[[length(terms) + 1L]] <- .mul(list(c, .pow(n, .int(j)), .pow(rs[[i]], n)))
    idx <- idx + 1L
  }
  yh <- .add(terms)
  yp <- S_Zero
  if (!is_zero_num(rest)) yp <- .rsolve_particular(coefs, lo, .neg(rest), n, rs, mult)
  sol <- .add(list(yh, yp))
  if (!is.null(init)) {
    eqs <- lapply(init, function(p) .sub(subs(sol, n, as_sym(p[[1L]])), as_sym(p[[2L]])))
    # inconsistent initial values: no solution, as SymPy's None
    k <- tryCatch(.solve_constants(eqs, cs), error = function(e) NULL)
    if (is.null(k)) return(NULL)
    sol <- simplify_if_cheap(xreplace(sol, k))
  }
  sol
}

# SymPy's rsolve() driver for polynomial or rational coefficients: clear
# denominators, shift the lowest term to y(n) and call rsolve_hyper().
.rsolve_general <- function(coefs, lo, rest, n, init) {
  coefs <- lapply(coefs, simplify)
  common <- S_One
  for (c in coefs) if (!is_polynomial(c, n)) common <- lcm(common, fraction(together(c))[[2L]], n)
  ind <- fraction(together(rest))
  if (is_polynomial(ind[[2L]], n)) common <- lcm(common, ind[[2L]], n)
  if (!identical(common$key, S_One$key)) {
    coefs <- lapply(coefs, function(c) {
      nd <- fraction(together(c))
      .mul(list(nd[[1L]], quo(common, nd[[2L]], n)))
    })
    rest <- .mul(list(ind[[1L]], quo(common, ind[[2L]], n)))
  }
  sh <- .add(list(n, .int(-lo)))
  coefs <- lapply(coefs, function(c) expand(subs(c, n, sh)))
  rest <- expand(subs(rest, n, sh))
  r <- .rsolve_hyper(coefs, .neg(rest), n)
  if (is.null(r)) return(NULL)
  sol <- r$result
  if (length(r$C) && length(init)) {
    eqs <- lapply(init, function(p) .sub(subs(sol, n, as_sym(p[[1L]])), as_sym(p[[2L]])))
    cs <- tryCatch(.solve_constants(eqs, r$C), error = function(e) NULL)
    if (is.null(cs)) return(NULL)
    sol <- simplify_if_cheap(xreplace(sol, cs))
  }
  sol
}

.rsolve_particular <- function(coefs, lo, r, n, rs, mult) {
  terms <- if (r$type == "Add") r$args else list(r)
  groups <- list()
  for (t in terms) {
    base <- S_One; deg <- 0L
    facs <- if (t$type == "Mul") t$args else list(t)
    for (f in facs) {
      if (!.has1(f, n)) next
      d <- .degree_in(f, n)
      if (!is.na(d)) {
        deg <- deg + d
      } else if (f$type == "Pow" && !.has1(f$args[[1L]], n) && identical(f$args[[2L]]$key, n$key)) {
        base <- f$args[[1L]]
      } else stop("rsolve() cannot handle this inhomogeneous term")
    }
    k <- base$key
    groups[[k]] <- list(base = base, deg = max(deg, if (is.null(groups[[k]])) 0L else groups[[k]]$deg))
  }
  unknowns <- list()
  trial <- list()
  for (g in groups) {
    s <- 0L
    for (i in seq_along(rs)) if (identical(rs[[i]]$key, g$base$key)) s <- mult[i]
    for (d in 0:g$deg) {
      u <- Dummy(paste0("A", length(unknowns)))
      unknowns[[length(unknowns) + 1L]] <- u
      trial[[length(trial) + 1L]] <- .mul(list(u, .pow(n, .int(d + s)), .pow(g$base, n)))
    }
  }
  Y <- .add(trial)
  L <- .add(lapply(seq_along(coefs), function(k) .mul(list(coefs[[k]], subs(Y, n, .add(list(n, .int(lo + k - 1L))))))))
  res <- expand(.sub(L, r))
  res <- expand_power_exp(res)
  gens <- .find_gens(list(res))
  pw <- Filter(function(g) g$type == "Pow" && identical(g$args[[2L]]$key, n$key), gens)
  bases <- unique(c(list(n), pw))
  dm <- lapply(pw, function(g) Dummy("g"))
  rs2 <- xreplace(res, Map(list, pw, dm))
  eqs <- list()
  for (t in if (rs2$type == "Add") rs2$args else list(rs2)) {
    mono <- S_One
    for (f in if (t$type == "Mul") t$args else list(t)) {
      b <- .base_exp(f)[[1L]]
      if (identical(b$key, n$key) || any(vapply(dm, function(d) identical(d$key, b$key), TRUE))) mono <- .mul(list(mono, f))
    }
    eqs[[mono$key]] <- c(eqs[[mono$key]], list(.div(t, mono)))
  }
  sol <- .linsolve_exprs(unname(lapply(eqs, .add)), unknowns)
  if (is.null(sol)) stop("could not find a particular solution")
  xreplace(Y, sol)
}
