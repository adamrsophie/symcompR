# Regions and integrals over them, after sympy/vector/parametricregion.py,
# implicitregion.py and integrals.py.

#' Parametric regions
#'
#' `ParametricRegion(definition, ...)` describes a region by its parametric
#' coordinates; each further argument is a parameter, or `list(parameter,
#' lower, upper)` giving its bounds. `region$definition`,
#' `region$parameters`, `region$limits` and `region$dimensions` read it.
#' `parametric_region_list()` converts a geometric entity (point, segment,
#' curve, circle, ellipse, polygon) or an [ImplicitRegion()] into a list of
#' parametric regions.
#'
#' @param definition An expression or a list of expressions.
#' @param ... Parameters, or `list(parameter, lower, upper)`.
#' @param reg A geometric entity or an implicit region.
#' @param parameter The name of the parameter (or names, for implicit
#'   regions).
#' @return A `ParametricRegion`, or a list of them.
#' @examples
#' t <- Symbol("t")
#' p <- ParametricRegion(list(t, t^2), list(t, 0, 1))
#' p$limits
#' parametric_region_list(Segment(Point(0, 0), Point(1, 1)))
#' @export
ParametricRegion <- function(definition, ...) {
  bounds <- list(...)
  for (bd in bounds) {
    if (is.list(bd) && !is_sym(bd) || (is_sym(bd) && bd$type == "Tuple")) {
      n <- if (is_sym(bd)) length(bd$args) else length(bd)
      if (n != 3L) .vec_stop("ValueError", "Tuple should be in the form (parameter, lowerbound, upperbound)")
    }
  }
  bargs <- lapply(bounds, function(bd) if (is.list(bd) && !is_sym(bd)) do.call(Tuple, lapply(bd, as_sym)) else as_sym(bd))
  def <- if (is.list(definition) && !is_sym(definition)) definition else if (is_sym(definition) && definition$type == "Tuple") definition$args else list(definition)
  o <- new_sym("ParametricRegion", c(list(do.call(Tuple, lapply(def, as_sym))), bargs))
  class(o) <- c("ParametricRegion", "sym")
  o
}

register_rebuilder("ParametricRegion", function(e, args) do.call(ParametricRegion, c(list(args[[1L]]), args[-1L])))

.pr_parameters <- function(reg) lapply(reg$args[-1L], function(b) if (b$type == "Tuple") b$args[[1L]] else b)
.pr_limits <- function(reg) {
  keys <- list(); vals <- list()
  for (b in reg$args[-1L]) if (b$type == "Tuple") {
    i <- match(b$args[[1L]]$key, vapply(keys, .k, ""))
    if (is.na(i)) { keys[[length(keys) + 1L]] <- b$args[[1L]]; vals[[length(vals) + 1L]] <- .pytuple(b$args[2:3]) }
    else vals[[i]] <- .pytuple(b$args[2:3])
  }
  .pydict(keys, vals)
}
.pr_limit <- function(lims, p) lims$values[[match(p$key, vapply(lims$keys, .k, ""))]]

#' @export
`$.ParametricRegion` <- function(x, name) {
  if (name %in% c("type", "args", "name", "value", "extra", "key")) return(.subset2(x, name))
  switch(name,
         definition = .subset2(x, "args")[[1L]],
         parameters = .pytuple(.pr_parameters(x)),
         limits = .pr_limits(x),
         dimensions = length(.pr_limits(x)$keys),
         .vec_stop("AttributeError", paste0("'ParametricRegion' object has no attribute '", name, "'")))
}

# ---------------------------------------------------------------------------
# Integrals over parametric regions

#' Integrals of fields over regions
#'
#' `ParametricIntegral(field, region)` integrates a scalar or vector field
#' over a [ParametricRegion()]: line integrals for one parameter, surface
#' (flux) integrals for two and volume integrals for three. When the
#' integral cannot be computed it stays unevaluated. `vector_integrate()`
#' also accepts implicit regions and geometric entities, and otherwise falls
#' back to [integrate()].
#'
#' @param field A scalar or vector field.
#' @param parametricregion,... A region, or integration limits.
#' @return An expression.
#' @examples
#' C <- CoordSys3D("C")
#' t <- Symbol("t")
#' ParametricIntegral(C$x * C$i + C$y * C$j, ParametricRegion(list(t, t), list(t, 0, 1)))
#' @export
ParametricIntegral <- function(field, parametricregion) {
  field <- as_sym(field)
  reg <- parametricregion
  cs <- .get_coord_systems(field)
  if (!length(cs)) coord_sys <- CoordSys3D("C")
  else if (length(cs) > 1L) .vec_stop("ValueError", "")
  else coord_sys <- cs[[1L]]
  lims <- .pr_limits(reg)
  dims <- length(lims$keys)
  if (dims == 0L) return(S_Zero)
  bv <- .cs(coord_sys)$base_vectors
  bs <- .cs(coord_sys)$base_scalars
  def <- reg$args[[1L]]$args
  pf <- field
  r <- .VZERO
  for (n in seq_along(def)) r <- .bd_add2(r, .mul_any(bv[[n]], def[[n]]))
  if (length(cs)) for (n in seq_along(def)) pf <- subs(pf, bs[[n]], def[[n]])
  dif <- function(v, p) if (.is_bd(v)) .bd_map(v, function(m) .diff1(m, p)) else .diff1(v, p)
  if (dims == 1L) {
    parameter <- .pr_parameters(reg)[[1L]]
    r_diff <- dif(r, parameter)
    lu <- .pr_limit(lims, parameter)
    integrand <- if (.is_vec(pf)) simplify(.vec_dot_method(r_diff, pf)) else simplify(.mul(list(magnitude(r_diff), pf)))
    result <- integrate(integrand, list(parameter, lu[[1L]], lu[[2L]]))
  } else if (dims == 2L) {
    uv <- .bounds_case(.pr_parameters(reg), lims)
    u <- uv[[1L]]; v <- uv[[2L]]
    normal_vector <- simplify(.vec_cross_method(dif(r, u), dif(r, v)))
    integrand <- if (.is_vec(pf)) .vec_dot_method(pf, normal_vector) else .mul(list(pf, magnitude(normal_vector)))
    integrand <- simplify(integrand)
    lu <- .pr_limit(lims, u); lv <- .pr_limit(lims, v)
    result <- integrate(integrand, list(u, lu[[1L]], lu[[2L]]), list(v, lv[[1L]], lv[[2L]]))
  } else {
    variables <- .bounds_case(.pr_parameters(reg), lims)
    coeff <- det(jacobian(def, variables))
    integrand <- simplify(.mul_any(pf, coeff))
    l <- lapply(variables, function(var) { lv <- .pr_limit(lims, var); list(var, lv[[1L]], lv[[2L]]) })
    result <- do.call(integrate, c(list(integrand), l))
  }
  if (!(is_sym(result) && result$type == "Integral")) return(result)
  new_sym("ParametricIntegral", list(field, reg))
}

# ParametricIntegral._bounds_case: order the parameters so that a parameter
# whose bounds use another one is integrated first
.bounds_case <- function(parameters, limits) {
  V <- limits$keys
  E <- list()
  for (p in V) {
    lu <- .pr_limit(limits, p)
    for (q in V) if (!identical(p$key, q$key) && (.contains_key(lu[[1L]], q$key) || .contains_key(lu[[2L]], q$key)))
      E[[length(E) + 1L]] <- list(p, q)
  }
  if (!length(E)) return(parameters)
  .vtopological_sort(V, E)
}

# sympy.utilities.iterables.topological_sort with key=default_sort_key
.vtopological_sort <- function(V, E) {
  L <- list()
  targets <- vapply(E, function(ed) ed[[2L]]$key, "")
  S <- Filter(function(v) !(v$key %in% targets), V)
  S <- .key_sort(lapply(S, .sort_key), S)
  S <- rev(S)
  kgt <- function(a, b) .key_cmp(.sort_key(a), .sort_key(b)) > 0
  while (length(S)) {
    node <- S[[length(S)]]
    S <- S[-length(S)]
    L[[length(L) + 1L]] <- node
    for (ed in E) {
      if (!identical(ed[[1L]]$key, node$key)) next
      j <- which(vapply(E, function(z) identical(z[[1L]]$key, ed[[1L]]$key) && identical(z[[2L]]$key, ed[[2L]]$key), TRUE))[1L]
      E <- E[-j]
      v <- ed[[2L]]
      if (any(vapply(E, function(z) identical(z[[2L]]$key, v$key), TRUE))) next
      pos <- NA_integer_
      for (s in seq_along(S)) if (kgt(v, S[[s]])) { pos <- s; break }
      S <- if (is.na(pos)) c(S, list(v)) else append(S, list(v), after = pos - 1L)
    }
  }
  if (length(E)) .vec_stop("ValueError", "cycle detected")
  L
}

#' @rdname ParametricIntegral
#' @export
vector_integrate <- function(field, ...) {
  region <- list(...)
  if (length(region) == 1L) {
    r1 <- region[[1L]]
    if (inherits(r1, "ParametricRegion")) return(ParametricIntegral(field, r1))
    if (inherits(r1, "ImplicitRegion")) return(vector_integrate(field, parametric_region_list(r1)[[1L]]))
    if (inherits(r1, "GeoEntity")) {
      result <- S_Zero
      for (reg in parametric_region_list(r1)) result <- .add_any(result, vector_integrate(field, reg))
      return(result)
    }
  }
  do.call(integrate, c(list(field), region))
}

#' @rdname ParametricRegion
#' @export
parametric_region_list <- function(reg, parameter = "t") {
  if (inherits(reg, "GeoPoint")) return(list(ParametricRegion(reg$coords)))
  if (inherits(reg, "GeoCurve")) {
    return(list(ParametricRegion(reg$functions, reg$limits)))
  }
  if (inherits(reg, c("GeoEllipse", "GeoCircle"))) {
    t <- Symbol(parameter, real = TRUE)
    c0 <- reg$center$coords
    hr <- if (inherits(reg, "GeoCircle")) reg$radius else reg$hradius
    vr <- if (inherits(reg, "GeoCircle")) reg$radius else reg$vradius
    def <- list(.add(list(c0[[1L]], .mul(list(hr, fcall("cos", list(t)))))), .add(list(c0[[2L]], .mul(list(vr, fcall("sin", list(t)))))))
    return(list(ParametricRegion(def, list(t, S_Zero, .mul(list(S_Two, S_Pi))))))
  }
  if (inherits(reg, "GeoSegment")) {
    t <- Symbol(parameter, real = TRUE)
    p1 <- reg$p1$coords; p2 <- reg$p2$coords
    def <- Map(function(a, b) simplify(.add(list(a, .mul(list(t, .sub(b, a)))))), p1, p2)
    bounds <- NULL
    for (n in seq_len(3L)) {
      # solve() of an identity is [] in SymPy
      sv <- function(eq) if (is_zero_num(eq)) list() else solve(eq, t)
      lo <- sv(.sub(def[[n]], p1[[n]]))
      hi <- sv(.sub(def[[n]], p2[[n]]))
      if (length(lo) == 1L && length(hi) == 1L) { bounds <- list(t, lo[[1L]], hi[[1L]]); break }
    }
    return(list(ParametricRegion(def, bounds)))
  }
  if (inherits(reg, "GeoPolygon")) {
    v <- reg$vertices
    n <- length(v)
    return(lapply(seq_len(n), function(m) parametric_region_list(Segment(v[[m]], v[[if (m == n) 1L else m + 1L]]), parameter)[[1L]]))
  }
  if (inherits(reg, "ImplicitRegion")) {
    if (identical(parameter, "t")) parameter <- c("t", "s")
    definition <- rational_parametrization(reg, parameter)
    bounds <- list()
    for (n in seq_len(length(reg$args[[1L]]$args) - 1L)) {
      p <- Symbol(parameter[[n]], real = TRUE)
      definition <- lapply(definition, function(el) trigsimp(subs(el, p, fcall("tan", list(.div(p, S_Two))))))
      bounds[[length(bounds) + 1L]] <- list(p, S_Zero, .mul(list(S_Two, S_Pi)))
    }
    return(list(do.call(ParametricRegion, c(list(definition), bounds))))
  }
  .vec_stop("ValueError", "SymPy cannot determine parametric representation of the region.")
}

# ---------------------------------------------------------------------------
# Implicit regions

#' Implicit regions
#'
#' `ImplicitRegion(variables, equation)` is the region where `equation`
#' (or an [Eq()]) vanishes. `region$variables`, `region$equation` and
#' `implicit_degree(region)` read it; `regular_point()`,
#' `singular_points()`, `multiplicity()` (via the region method below) and
#' `rational_parametrization()` port the methods of SymPy's
#' `ImplicitRegion`, and `conic_coeff()` gives the coefficients of a conic.
#'
#' @param variables A list of symbols.
#' @param equation An expression or equation.
#' @param region An implicit region.
#' @param parameters Parameter names.
#' @param reg_point A regular point to use for a conic.
#' @return An `ImplicitRegion`, expression, list or set.
#' @examples
#' x <- Symbol("x"); y <- Symbol("y")
#' r <- ImplicitRegion(list(x, y), x^2 + y^2 - 4)
#' regular_point(r)
#' rational_parametrization(r)
#' @export
ImplicitRegion <- function(variables, equation) {
  variables <- if (is_sym(variables) && variables$type == "Tuple") variables else do.call(Tuple, lapply(variables, as_sym))
  equation <- as_sym(equation)
  if (equation$type == "Equality") equation <- .sub(equation$args[[1L]], equation$args[[2L]])
  o <- new_sym("ImplicitRegion", list(variables, equation))
  class(o) <- c("ImplicitRegion", "sym")
  o
}

register_rebuilder("ImplicitRegion", function(e, args) ImplicitRegion(args[[1L]], args[[2L]]))

#' @export
`$.ImplicitRegion` <- function(x, name) {
  if (name %in% c("type", "args", "name", "value", "extra", "key")) return(.subset2(x, name))
  switch(name,
         variables = .subset2(x, "args")[[1L]],
         equation = .subset2(x, "args")[[2L]],
         degree = implicit_degree(x),
         .vec_stop("AttributeError", paste0("'ImplicitRegion' object has no attribute '", name, "'")))
}

#' @rdname ImplicitRegion
#' @export
implicit_degree <- function(region) .int(total_degree(region$args[[2L]]))

.ir_vars <- function(region) region$args[[1L]]$args

#' @rdname ImplicitRegion
#' @export
conic_coeff <- function(variables, equation) {
  if (total_degree(equation) != 2) .vec_stop("ValueError", "")
  variables <- if (is_sym(variables) && variables$type == "Tuple") variables$args else variables
  x <- variables[[1L]]; y <- variables[[2L]]
  equation <- expand(equation)
  a <- coeff(equation, .sq(x))
  b <- coeff(equation, .mul(list(x, y)))
  cc <- coeff(equation, .sq(y))
  d <- coeff(coeff(equation, x, 1), y, 0)
  e <- coeff(coeff(equation, y, 1), x, 0)
  f <- coeff(coeff(equation, x, 0), y, 0)
  .pytuple(list(a, b, cc, d, e, f))
}

#' @rdname ImplicitRegion
#' @export
regular_point <- function(region) {
  equation <- region$args[[2L]]
  vars <- .ir_vars(region)
  if (length(vars) == 1L) return(.pytuple(list(.set_elements(solveset(equation, vars[[1L]], domain = Reals))[[1L]])))
  if (length(vars) == 2L && implicit_degree(region)$value == 2) {
    cf <- conic_coeff(vars, equation)
    a <- cf[[1L]]; b <- cf[[2L]]; cc <- cf[[3L]]
    if (identical(.sq(b)$key, .mul(list(.int(4), a, cc))$key)) return(.pytuple(do.call(.regular_point_parabola, cf)))
    return(.pytuple(do.call(.regular_point_ellipse, cf)))
  }
  if (length(vars) == 3L) {
    for (xr in -10:9) for (yr in -10:9) {
      e2 <- subs(equation, list(list(vars[[1L]], .int(xr)), list(vars[[2L]], .int(yr))))
      sols <- solveset(e2, vars[[3L]], domain = Reals)
      if (!identical(sols$key, S_EmptySet$key))
        return(.pytuple(list(.int(xr), .int(yr), .set_elements(solveset(e2, vars[[3L]]))[[1L]])))
    }
  }
  .vec_stop("NotImplementedError", "")
}

.set_elements <- function(s) if (s$type == "FiniteSet") .sympy_set_order(s) else list()
.sympy_set_order <- function(s) s$args

.regular_point_parabola <- function(a, b, c, d, e, f) {
  nz <- function(z) !is_zero_num(z)
  ok <- (nz(a) || nz(d)) && (nz(c) || nz(e)) && identical(.sq(b)$key, .mul(list(.int(4), a, c))$key) && (nz(a) || nz(c))
  if (!ok) .vec_stop("ValueError", "Rational Point on the conic does not exist")
  if (nz(a)) {
    d_dash <- .sub(.mul(list(.int(4), a, e)), .mul(list(S_Two, b, d)))
    f_dash <- .sub(.mul(list(.int(4), a, f)), .sq(d))
    if (nz(d_dash)) {
      y_reg <- .neg(.div(f_dash, d_dash))
      x_reg <- .neg(.div(.add(list(d, .mul(list(b, y_reg)))), .mul(list(S_Two, a))))
    } else ok <- FALSE
  } else if (nz(c)) {
    d_dash <- .sub(.mul(list(.int(4), c, d)), .mul(list(S_Two, b, e)))
    f_dash <- .sub(.mul(list(.int(4), c, f)), .sq(e))
    if (nz(d_dash)) {
      x_reg <- .neg(.div(f_dash, d_dash))
      y_reg <- .neg(.div(.add(list(e, .mul(list(b, x_reg)))), .mul(list(S_Two, c))))
    } else ok <- FALSE
  }
  if (!ok) .vec_stop("ValueError", "Rational Point on the conic does not exist")
  list(x_reg, y_reg)
}

.regular_point_ellipse <- function(a, b, c, d, e, f) {
  D <- .sub(.mul(list(.int(4), a, c)), .sq(b))
  if (is_zero_num(D)) .vec_stop("ValueError", "Rational Point on the conic does not exist")
  if (is_zero_num(a) && is_zero_num(c)) {
    K <- S_NegativeOne
    L <- .mul(list(.int(4), .sub(.mul(list(d, e)), .mul(list(b, f)))))
  } else if (!is_zero_num(c)) {
    K <- D
    L <- .add(list(.mul(list(.int(4), .sq(c), .sq(d))), .mul(list(.int(-4), b, c, d, e)), .mul(list(.int(4), a, c, .sq(e))),
                   .mul(list(.int(4), .sq(b), c, f)), .mul(list(.int(-16), a, .sq(c), f))))
  } else {
    K <- D
    L <- .add(list(.mul(list(.int(4), .sq(a), .sq(e))), .mul(list(.int(-4), b, a, d, e)), .mul(list(.int(4), .sq(b), a, f))))
  }
  ok <- !is_zero_num(L) && !(num_sign(K) > 0 && num_sign(L) < 0)
  if (!ok) .vec_stop("ValueError", "Rational Point on the conic does not exist")
  Kq <- .nt_limit_den(K); Lq <- .nt_limit_den(L)
  k1 <- gmp::numerator(Kq); k2 <- gmp::denominator(Kq)
  l1 <- gmp::numerator(Lq); l2 <- gmp::denominator(Lq)
  g <- gmp::gcd(k2, l2)
  a1 <- gmp::as.bigq(l2 * k2, g); b1 <- gmp::as.bigq(k1 * l2, g); c1 <- -gmp::as.bigq(l1 * k2, g)
  sq_part <- function(v) {
    cv <- core(gmp::as.bigz(abs(v)), 2)
    cv <- if (is_sym(cv)) gmp::numerator(cv$value) else gmp::as.bigz(cv)
    gmp::as.bigq(sign(as.numeric(v)) * cv)
  }
  a2 <- sq_part(a1); b2 <- sq_part(b1); c2 <- sq_part(c1)
  r1 <- .sqrt(make_rational(a1 / a2)); r2 <- .sqrt(make_rational(b1 / b2)); r3 <- .sqrt(make_rational(c1 / c2))
  qg <- function(p, q) gmp::as.bigq(gmp::gcd(gmp::numerator(p), gmp::numerator(q)))
  g <- qg(qg(a2, b2), c2)
  a2 <- a2 / g; b2 <- b2 / g; c2 <- c2 / g
  g1 <- qg(a2, b2); a2 <- a2 / g1; b2 <- b2 / g1; c2 <- c2 * g1
  g2 <- qg(a2, c2); a2 <- a2 / g2; b2 <- b2 * g2; c2 <- c2 / g2
  g3 <- qg(b2, c2); a2 <- a2 * g3; b2 <- b2 / g3; c2 <- c2 / g3
  X <- Symbol("x"); Y <- Symbol("y"); Z <- Symbol("z")
  eq <- .add(list(.mul(list(make_rational(a2), .sq(X))), .mul(list(make_rational(b2), .sq(Y))), .mul(list(make_rational(c2), .sq(Z)))))
  sols <- diophantine(eq)
  sols <- if (is_sym(sols)) .set_elements(sols) else sols
  if (!length(sols)) .vec_stop("ValueError", "Rational Point on the conic does not exist")
  flag <- FALSE
  xs <- list(X, Y, Z)
  for (sol in sols) {
    sl <- if (is_sym(sol)) sol$args else sol
    syms <- free_symbols(do.call(Tuple, sl))
    rep <- lapply(syms, function(s) list(s, .int(3)))
    setrep <- function(sym, val) c(Filter(function(r) !identical(r[[1L]]$key, sym$key), rep), list(list(sym, val)))
    sol_z <- sl[[3L]]
    if (is_zero_num(sol_z)) { flag <- TRUE; next }
    if (is_number_type(sol_z)) next
    sz <- free_symbols(sol_z)
    if (length(sz) == 1L) {
      p <- sz[[1L]]
      rep <- setrep(p, .first_integer_not_root(sol_z, p))
    }
    if (length(sz) == 2L) {
      p <- sz[[1L]]; q <- sz[[2L]]
      # S.Integers iterates 0, 1, -1, 2, -2, ...
      for (k in .integers_seq(50L)) {
        qv <- .first_integer_not_root(subs(sol_z, p, k), q)
        if (!is.null(qv)) { rep <- setrep(p, k); rep <- setrep(q, qv); break }
      }
    }
    xs <- if (length(syms)) lapply(sl, function(s) .subs_dict(s, rep)) else sl
    flag <- FALSE
    break
  }
  if (flag) .vec_stop("ValueError", "Rational Point on the conic does not exist")
  x <- .div(.mul(list(xs[[1L]], make_rational(g3))), r1)
  y <- .div(.mul(list(xs[[2L]], make_rational(g2))), r2)
  z <- .div(.mul(list(xs[[3L]], make_rational(g1))), r3)
  x <- .div(x, z); y <- .div(y, z)
  if (is_zero_num(a) && is_zero_num(c)) {
    x_reg <- .div(.sub(.add(list(x, y)), .mul(list(S_Two, e))), .mul(list(S_Two, b)))
    y_reg <- .div(.sub(.sub(x, y), .mul(list(S_Two, d))), .mul(list(S_Two, b)))
  } else if (!is_zero_num(c)) {
    x_reg <- .div(.add(list(x, .mul(list(.int(-2), d, c)), .mul(list(b, e)))), K)
    y_reg <- .div(.sub(.sub(y, .mul(list(b, x_reg))), e), .mul(list(S_Two, c)))
  } else {
    y_reg <- .div(.add(list(x, .mul(list(.int(-2), e, a)), .mul(list(b, d)))), K)
    x_reg <- .div(.sub(.sub(y, .mul(list(b, y_reg))), d), .mul(list(S_Two, a)))
  }
  list(x_reg, y_reg)
}

.integers_seq <- function(n) lapply(c(0L, as.vector(rbind(seq_len(n), -seq_len(n)))), .int)

# next(iter(Complement(S.Integers, solveset(Eq(e, 0), p, S.Integers)))), or
# NULL when every integer is a root
.first_integer_not_root <- function(e, p) {
  e <- expand(e)
  if (is_zero_num(e)) return(NULL)
  roots <- if (.has1(e, p)) Filter(function(r) is_sym(r) && r$type == "Integer", solve(e, p)) else list()
  rk <- vapply(roots, .k, "")
  for (k in .integers_seq(1000L)) if (!(k$key %in% rk)) return(k)
  NULL
}

# Rational(v).limit_denominator(10**12) for an exact number
.nt_limit_den <- function(v) {
  q <- if (v$type %in% c("Integer", "Rational")) v$value else gmp::as.bigq(num_value_dbl(v))
  q
}

#' @rdname ImplicitRegion
#' @export
singular_points <- function(region) {
  eq <- region$args[[2L]]
  vars <- .ir_vars(region)
  eqs <- c(list(eq), lapply(vars, function(v) .diff1(eq, v)))
  sols <- solve(eqs, vars)
  if (!length(sols)) return(S_EmptySet)
  do.call(FiniteSet, lapply(sols, function(s) do.call(Tuple, lapply(vars, function(v) if (is.null(s[[v$name]])) v else s[[v$name]]))))
}

.ir_multiplicity <- function(region, point) {
  if (inherits(point, "GeoPoint")) point <- point$coords
  if (is_sym(point) && point$type == "Tuple") point <- point$args
  point <- lapply(point, as_sym)
  me <- region$args[[2L]]
  vars <- .ir_vars(region)
  for (n in seq_along(vars)) me <- subs(me, vars[[n]], .add(list(vars[[n]], point[[n]])))
  me <- expand(me)
  if (length(me$args)) {
    m <- min(vapply(me$args, function(tm) as.numeric(total_degree(tm)), 0))
  } else {
    m <- as.numeric(total_degree(me))
  }
  .int(m)
}

#' @rdname ImplicitRegion
#' @export
rational_parametrization <- function(region, parameters = c("t", "s"), reg_point = NULL) {
  equation <- region$args[[2L]]
  degree <- implicit_degree(region)$value
  vars <- .ir_vars(region)
  if (degree == 1) {
    if (length(vars) == 1L) return(.pytuple(list(equation)))
    if (length(vars) == 2L) {
      y_par <- .set_elements(solveset(equation, vars[[2L]]))[[1L]]
      return(.pytuple(list(vars[[1L]], y_par)))
    }
    .vec_stop("NotImplementedError", "")
  }
  point <- list()
  sp <- NULL
  singular <- function() { if (is.null(sp)) sp <<- singular_points(region); sp }
  if (degree == 2) {
    if (!is.null(reg_point)) point <- reg_point
    else if (length(singular()$args)) point <- singular()$args[[1L]]$args
    else point <- regular_point(region)
  }
  if (length(singular()$args)) {
    for (spoint in singular()$args) {
      pt <- spoint$args
      syms <- free_symbols(spoint)
      if (length(syms)) pt <- lapply(pt, function(s) subs(s, lapply(syms, function(z) list(z, S_Two))))
      if (.ir_multiplicity(region, pt)$value == degree - 1) { point <- pt; break }
    }
  }
  if (!length(point)) .vec_stop("NotImplementedError", "")
  me <- equation
  for (n in seq_along(vars)) me <- subs(me, vars[[n]], .add(list(vars[[n]], point[[n]])))
  me <- expand(me)
  hn <- S_Zero; hn_1 <- S_Zero
  for (tm in me$args) {
    if (total_degree(tm) == degree) hn <- .add(list(hn, tm)) else hn_1 <- .add(list(hn_1, tm))
  }
  hn_1 <- .neg(hn_1)
  parameters <- as.character(parameters)
  if (length(vars) == 2L) {
    p1 <- parameters[[1L]]
    s <- Symbol(if (p1 == "s") "s_" else "s", real = TRUE)
    t <- Symbol(p1, real = TRUE)
    hn <- .subs_dict(hn, list(list(vars[[1L]], s), list(vars[[2L]], t)))
    hn_1 <- .subs_dict(hn_1, list(list(vars[[1L]], s), list(vars[[2L]], t)))
    x_par <- .add(list(subs(.mul(list(s, .div(hn_1, hn))), s, S_One), point[[1L]]))
    y_par <- .add(list(subs(.mul(list(t, .div(hn_1, hn))), s, S_One), point[[2L]]))
    return(.pytuple(list(x_par, y_par)))
  }
  if (length(vars) == 3L) {
    rr <- Symbol(if ("r" %in% parameters) "r_" else "r", real = TRUE)
    s <- Symbol(parameters[[2L]], real = TRUE)
    t <- Symbol(parameters[[1L]], real = TRUE)
    rep <- list(list(vars[[1L]], rr), list(vars[[2L]], s), list(vars[[3L]], t))
    hn <- .subs_dict(hn, rep); hn_1 <- .subs_dict(hn_1, rep)
    q <- .div(hn_1, hn)
    return(.pytuple(list(.add(list(subs(.mul(list(rr, q)), rr, S_One), point[[1L]])),
                         .add(list(subs(.mul(list(s, q)), rr, S_One), point[[2L]])),
                         .add(list(subs(.mul(list(t, q)), rr, S_One), point[[3L]])))))
  }
  .vec_stop("NotImplementedError", "")
}
