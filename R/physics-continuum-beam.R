# Beams, after sympy/physics/continuum_mechanics/beam.py (Beam and
# Beam3D). The loads are singularity functions of the coordinate along the
# beam; the reactions come from the boundary conditions by linsolve().

#' Beams
#'
#' `Beam(length, elastic_modulus, second_moment)` is a beam bent in a plane.
#' `$apply_load(value, start, order, end)` adds a load,
#' `$apply_support(loc, type)` a `"pin"`, `"roller"` or `"fixed"` support
#' and `$apply_rotation_hinge()` and `$apply_sliding_hinge()` hinges;
#' `$solve_for_reaction_loads()` finds the reactions. `$shear_force()`,
#' `$bending_moment()`, `$slope()`, `$deflection()` and `$shear_stress()`
#' give the curves along the beam, `$max_shear_force()`, `$max_bmoment()`,
#' `$max_deflection()` and `$point_cflexure()` their extremes, and
#' `$solve_for_ild_reactions()`, `$solve_for_ild_shear()` and
#' `$solve_for_ild_moment()` the influence lines. The boundary conditions are
#' `$bc_deflection`, `$bc_slope`, `$bc_bending_moment` and
#' `$bc_shear_force`. `$join()` joins two beams. The `$plot_*()` methods draw
#' base R plots. `Beam3D(length, elastic_modulus, shear_modulus,
#' second_moment, area)` is a beam loaded in three dimensions. After
#' `sympy.physics.continuum_mechanics.beam`.
#'
#' @param length,elastic_modulus The length and the elastic modulus.
#' @param second_moment The second moment of area, or a geometry entity
#'   giving the cross section.
#' @param area The cross section area.
#' @param variable The coordinate along the beam.
#' @param base_char The base name of the integration constants.
#' @param ild_variable The position of the moving load of the influence
#'   lines.
#' @param shear_modulus The shear modulus.
#' @return A `Beam` or `Beam3D` object.
#' @examples
#' E <- Symbol("E"); I <- Symbol("I")
#' b <- Beam(30, E, I)
#' R1 <- b$apply_support(10, "roller"); R2 <- b$apply_support(30, "roller")
#' b$apply_load(-8, 0, -1)
#' b$solve_for_reaction_loads(R1, R2)
#' b$reaction_loads
#' b$bending_moment()
#' @export
Beam <- function(length, elastic_modulus, second_moment, area = Symbol("A"), variable = Symbol("x"),
                 base_char = "C", ild_variable = Symbol("a")) {
  self <- .pyenv_new("Beam")
  .beam_init(self, length, elastic_modulus, second_moment, area, variable, base_char, ild_variable)
  self
}

.beam_init <- function(self, length, elastic_modulus, second_moment, area, variable, base_char, ild_variable) {
  self$length <- length
  self$elastic_modulus <- elastic_modulus
  if (.beam_is_geo(second_moment)) {
    self$cross_section <- second_moment
  } else {
    .cm_set(self, "_cross_section", NULL)
    self$second_moment <- second_moment
  }
  self$variable <- variable
  .cm_set(self, "ild_variable", ild_variable)
  .cm_set(self, "_base_char", base_char)
  for (f in c("_bc_deflection", "_bc_slope", "_bc_bending_moment", "_bc_shear_force", "_applied_supports",
              "_applied_support_symbols", "_applied_rotation_hinges", "_applied_sliding_hinges",
              "_rotation_hinge_symbols", "_sliding_hinge_symbols", "_support_as_loads", "_applied_loads"))
    .cm_set(self, f, list())
  .cm_set(self, "_load", S_Zero)
  self$area <- area
  empty <- .pydict(list(), list())
  for (f in c("_reaction_loads", "_ild_reactions", "_rotation_jumps", "_deflection_jumps",
              "_ild_rotations_jumps", "_ild_deflection_jumps")) .cm_set(self, f, empty)
  .cm_set(self, "_ild_shear", S_Zero)
  .cm_set(self, "_ild_moment", S_Zero)
  .cm_set(self, "_original_load", S_Zero)
  .cm_set(self, "_joined_beam", FALSE)
  invisible(self)
}

.beam_is_geo <- function(s) inherits(s, c("GeoEntity", "GeometryEntity"))

.beam_dict <- function(keys, values) .pydict(keys, values)

# the points of a circle or polygon as Point objects, as SymPy's geometry
# constructors store them
.cm_geo_points <- function(s) {
  pt <- function(p) if (inherits(p, "GeoPoint")) p else Point(unclass(p))
  if (inherits(s, "GeoCircle")) s$center <- pt(s$center)
  else if (inherits(s, "GeoPolygon")) s$vertices <- lapply(s$vertices, pt)
  s
}

# second_moment_of_area() of a circle (Ellipse) or a polygon about its
# centroid: (I_xx, I_yy, I_xy)
.cm_second_moment_of_area <- function(s) {
  if (inherits(s, "GeoCircle")) {
    hr <- s$radius; vr <- Abs(s$radius)
    return(list(S_Pi * hr * vr^3 / 4L, S_Pi * hr^3 * vr / 4L, S_Zero))
  }
  if (!inherits(s, "GeoPolygon"))
    .physics_error("AttributeError", "object has no attribute 'second_moment_of_area'")
  v <- s$vertices; n <- length(v)
  I_xx <- S_Zero; I_yy <- S_Zero; I_xy <- S_Zero
  for (i in seq_len(n)) {
    p1 <- v[[if (i == 1L) n else i - 1L]]$coords; p2 <- v[[i]]$coords
    x1 <- p1[[1L]]; y1 <- p1[[2L]]; x2 <- p2[[1L]]; y2 <- p2[[2L]]
    w <- x1 * y2 - x2 * y1
    I_xx <- I_xx + (y1^2 + y1 * y2 + y2^2) * w
    I_yy <- I_yy + (x1^2 + x1 * x2 + x2^2) * w
    I_xy <- I_xy + (x1 * y2 + 2L * x1 * y1 + 2L * x2 * y2 + x2 * y1) * w
  }
  A <- area(s); cen <- centroid(s)$coords
  c_x <- cen[[1L]]; c_y <- cen[[2L]]
  list(I_xx / 12L - A * c_y^2, I_yy / 12L - A * c_x^2, I_xy / 24L - A * (c_x * c_y))
}

.beam_dict_get <- function(d, k) {
  i <- match(as_sym(k)$key, vapply(d$keys, function(z) as_sym(z)$key, ""))
  if (is.na(i)) .physics_error("KeyError", .str(k))
  d$values[[i]]
}

.beam_add_load <- function(self, term) {
  .cm_set(self, "_load", .cm_get(self, "_load") + term)
  .cm_set(self, "_original_load", .cm_get(self, "_original_load") + term)
}

.beam_apply_load <- function(self, value, start, order, end = NULL) {
  x <- self$variable
  value <- as_sym(value); start <- as_sym(start)
  .cm_push(self, "_applied_loads", .pytuple(list(value, start, order, end)))
  .beam_add_load(self, value * SingularityFunction(x, start, order))
  if (!is.null(end)) {
    if (isTRUE(as_sym(order) < 0))
      .physics_error("ValueError", "If 'end' is provided the 'order' of the load cannot be negative, i.e. 'end' is only valid for distributed loads.")
    .beam_add_load(self, -.cm_taylor(x, value, start, as_sym(order), as_sym(end)))
  }
  invisible(NULL)
}

.beam_remove_load <- function(self, value, start, order, end = NULL) {
  x <- self$variable
  value <- as_sym(value); start <- as_sym(start)
  loads <- .cm_get(self, "_applied_loads")
  i <- which(vapply(loads, function(t) .cm_tuple_eq(unclass(t), list(value, start, order, end)), TRUE))
  if (!length(i)) .physics_error("ValueError", "No such load distribution exists on the beam object.")
  .beam_add_load(self, -value * SingularityFunction(x, start, order))
  .cm_set(self, "_applied_loads", loads[-i[[1L]]])
  if (!is.null(end)) {
    if (isTRUE(as_sym(order) < 0))
      .physics_error("ValueError", "If 'end' is provided the 'order' of the load cannot be negative, i.e. 'end' is only valid for distributed loads.")
    .beam_add_load(self, .cm_taylor(x, value, start, as_sym(order), as_sym(end)))
  }
  invisible(NULL)
}

.beam_apply_support <- function(self, loc, type = "fixed") {
  loc <- as_sym(loc)
  if (!type %in% c("fixed", "roller", "pin"))
    .physics_error("ValueError", "Invalid support type. Choose from 'pin', 'roller', or 'fixed'.")
  .cm_push(self, "_applied_supports", .pytuple(list(loc, type)))
  R <- Symbol(paste0("R_", .str(loc)))
  .cm_push(self, "_applied_support_symbols", R)
  if (type %in% c("pin", "roller")) {
    self$apply_load(R, loc, -1L)
    .cm_push(self, "_bc_deflection", .pytuple(list(loc, 0L)))
  } else {
    M <- Symbol(paste0("M_", .str(loc)))
    .cm_push(self, "_applied_support_symbols", M)
    self$apply_load(R, loc, -1L)
    self$apply_load(M, loc, -2L)
    .cm_push(self, "_bc_deflection", .pytuple(list(loc, 0L)))
    .cm_push(self, "_bc_slope", .pytuple(list(loc, 0L)))
    .cm_push(self, "_support_as_loads", .pytuple(list(M, loc, -2L, NULL)))
  }
  .cm_push(self, "_support_as_loads", .pytuple(list(R, loc, -1L, NULL)))
  if (type %in% c("pin", "roller")) R else .pytuple(list(R, M))
}

.beam_get_I <- function(self, loc) {
  I <- self$second_moment
  if (!(is_sym(I) && I$type == "Piecewise")) return(I)
  for (p in .pw_pairs(I)) if (isTRUE(loc <= p[[2L]]$args[[2L]])) return(p[[1L]])
  NULL
}

.beam_hinge <- function(self, loc, prefix, order, locs, syms, bc) {
  loc <- as_sym(loc)
  jump <- Symbol(paste0(prefix, .str(loc)))
  .cm_push(self, locs, loc)
  .cm_push(self, syms, jump)
  self$apply_load(self$elastic_modulus * .beam_get_I(self, loc) * jump, loc, order)
  .cm_push(self, bc, .pytuple(list(loc, 0L)))
  jump
}

# the terms of eqs.args without an infinite number among their arguments
.beam_finite_terms <- function(eqs) {
  inf <- function(n) is_sym(n) && isTRUE(ask_fact(n, "infinite"))
  keep <- Filter(function(a) !any(vapply(a$args, inf, TRUE)), eqs$args)
  if (!length(keep)) return(S_Zero)
  .add(keep)
}

.beam_unique <- function(l) l[!duplicated(vapply(l, function(z) z$key, ""))]

.beam_solve <- function(self, ...) {
  reactions <- list(...)
  if (anyDuplicated(vapply(reactions, function(r) r$key, "")))
    .physics_error("ValueError", "Duplicate Symbols passed into solve_for_reaction_loads()")
  x <- self$variable; l <- self$length
  C3 <- Symbol("C3"); C4 <- Symbol("C4")
  rotation_jumps <- .cm_get(self, "_rotation_hinge_symbols")
  deflection_jumps <- .cm_get(self, "_sliding_hinge_symbols")
  sf <- self$shear_force(); bm <- self$bending_moment()
  eqs <- list(limit(sf, x, l), limit(bm, x, l))
  for (bc in .cm_get(self, "_bc_shear_force")) eqs[[length(eqs) + 1L]] <- .beam_finite_terms(subs(sf, x, bc[[1L]]) - bc[[2L]])
  for (bc in .cm_get(self, "_bc_bending_moment")) eqs[[length(eqs) + 1L]] <- .beam_finite_terms(subs(bm, x, bc[[1L]]) - bc[[2L]])
  slope_curve <- integrate(bm, x) + C3
  for (bc in .cm_get(self, "_bc_slope")) eqs[[length(eqs) + 1L]] <- subs(slope_curve, x, bc[[1L]]) - bc[[2L]]
  deflection_curve <- integrate(slope_curve, x) + C4
  for (bc in .cm_get(self, "_bc_deflection")) eqs[[length(eqs) + 1L]] <- subs(deflection_curve, x, bc[[1L]]) - bc[[2L]]
  total <- .beam_unique(c(.cm_get(self, "_applied_support_symbols"), reactions))
  sol <- .cm_linsolve(eqs, c(list(C3, C4), total, rotation_jumps, deflection_jumps))
  nt <- length(total); nr <- length(rotation_jumps)
  .cm_set(self, "_reaction_loads", .pydict(total, sol[2L + seq_len(nt)]))
  .cm_set(self, "_rotation_jumps", .pydict(rotation_jumps, sol[2L + nt + seq_len(nr)]))
  .cm_set(self, "_deflection_jumps", .pydict(deflection_jumps, sol[2L + nt + nr + seq_along(deflection_jumps)]))
  ld <- .cm_get(self, "_load")
  for (d in c("_reaction_loads", "_rotation_jumps", "_deflection_jumps")) ld <- .cm_dict_subs(ld, .cm_get(self, d))
  .cm_set(self, "_load", ld)
  invisible(NULL)
}

# sorted(): Python compares with __lt__, which raises TypeError when the
# relational stays undecided
.cm_sorted <- function(l) {
  lt <- function(a, b) {
    r <- as_sym(a) < as_sym(b)
    if (!is.logical(r))
      .physics_error("TypeError", paste0("cannot determine truth value of Relational: ", sstr(r)))
    r
  }
  out <- list()
  for (v in l) {
    k <- length(out)
    while (k >= 1L && lt(v, out[[k]])) k <- k - 1L
    out <- append(out, list(v), after = k)
  }
  out
}

# the singularity points of a curve, as max_shear_force() reads them
.beam_singularities <- function(curve) {
  terms <- if (curve$type == "Add") curve$args else curve$args
  pts <- lapply(terms, function(t) {
    # term.args[-1]: SymPy orders the SingularityFunction last in a Mul
    if (t$type == "Mul") {
      sf <- Filter(function(a) identical(a$name, "SingularityFunction"), t$args)
      t <- if (length(sf)) sf[[1L]] else t$args[[length(t$args)]]
    }
    t$args[[2L]]
  })
  pts <- .beam_unique(pts)
  .cm_sorted(pts)
}

.beam_max_curve <- function(self, curve, slope) {
  x <- self$variable
  sing <- .beam_singularities(curve)
  intervals <- list(); values <- list()
  for (i in seq_along(sing)) {
    s <- sing[[i]]
    if (is_zero_num(s)) next
    prev <- sing[[if (i == 1L) length(sing) else i - 1L]]
    # solve() of the slope between two singularity points; SymPy cannot
    # solve a slope that vanishes on the whole interval
    p <- .cm_on_interval(slope, x, prev)
    if (!is_zero_num(p)) {
      points <- .cm_roots_in(p, x, prev, s, closed = FALSE)
      val <- lapply(points, function(pt) Abs(subs(curve, x, pt)))
      points <- c(points, list(prev, s))
      val <- c(val, list(Abs(limit(curve, x, prev, dir = "+")), Abs(limit(curve, x, s, dir = "-"))))
      k <- .beam_argmax(val)
      values[[length(values) + 1L]] <- val[[k]]
      intervals[[length(intervals) + 1L]] <- points[[k]]
    } else {
      init <- limit(curve, x, prev, dir = "+"); fin <- limit(curve, x, s, dir = "-")
      mid <- subs(curve, x, (prev + s) / 2L)
      if (identical(mid$key, ((init + fin) / 2L)$key) && !identical(init$key, fin$key)) {
        values <- c(values, list(init, fin)); intervals <- c(intervals, list(prev, s))
      } else {
        values[[length(values) + 1L]] <- fin
        intervals[[length(intervals) + 1L]] <- Interval(prev, s)
      }
    }
  }
  values <- lapply(values, Abs)
  k <- .beam_argmax(values)
  .pytuple(list(intervals[[k]], values[[k]]))
}

# the SingularityFunction applications in e
.cm_sf_atoms <- function(e) Filter(function(z) z$name == "SingularityFunction", atoms(e, "Function"))

# an expression of singularity functions on the interval that starts at lo
.cm_on_interval <- function(e, x, lo) {
  e <- as_sym(e)
  for (sf in .cm_sf_atoms(e)) {
    a <- sf$args[[2L]]; n <- sf$args[[3L]]
    v <- if (isTRUE(n < 0) || !isTRUE(a <= lo)) S_Zero else (sf$args[[1L]] - a)^n
    e <- xreplace(e, list(list(sf, v)))
  }
  expand(e)
}

# the real roots of p in (lo, hi), or in [lo, hi) when closed
.cm_roots_in <- function(p, x, lo, hi, closed = FALSE) {
  r <- solve(p, x)
  if (is.null(r)) return(list())
  if (!is.list(r) || is_sym(r)) r <- list(r)
  Filter(function(z) isTRUE(ask_fact(z, "real")) && (if (closed) isTRUE(z >= lo) else isTRUE(z > lo)) && isTRUE(z < hi), r)
}

# max() and list.index() of comparable values: the first maximum
.beam_argmax <- function(vals) {
  k <- 1L
  for (i in seq_along(vals)[-1L]) if (isTRUE(vals[[i]] > vals[[k]])) k <- i
  k
}

.beam_point_cflexure <- function(self) {
  x <- self$variable; l <- self$length
  bm <- self$bending_moment()
  terms <- list()
  for (term in (if (bm$type == "Add") bm$args else bm$args)) {
    sfs <- .cm_sf_atoms(term)
    if (!length(sfs)) terms[[length(terms) + 1L]] <- term
    else for (sf in sfs) if (isTRUE(sf$args[[3L]] >= 0)) terms[[length(terms) + 1L]] <- term
  }
  nsb <- if (length(terms)) .add(terms) else S_Zero
  locs <- .beam_unique(lapply(.cm_sf_atoms(nsb), function(t) t$args[[2L]]))
  locs <- locs[order(vapply(locs, function(p) as.numeric(N(p)), 0))]
  pts <- c(list(S_Zero), Filter(function(p) !is_zero_num(p) && !identical(p$key, l$key), locs), list(l))
  roots <- list()
  for (i in seq_len(length(pts) - 1L)) {
    st <- pts[[i]]; en <- pts[[i + 1L]]
    if (identical(st$key, en$key)) next
    iv <- Interval(st, en)
    roots <- c(roots, list(st, en))
    # solveset() rewrites the singularity functions as Piecewise; inside the
    # interval that is the polynomial the active terms give
    sol <- solveset(.cm_on_interval(nsb, x, st), x, iv)
    if (sol$type == "FiniteSet") roots <- c(roots, sol$args)
    else if (sol$type == "Interval") roots <- c(roots, list(sol$args[[1L]], sol$args[[2L]]))
  }
  roots <- .beam_unique(roots)
  roots <- .cm_sorted(roots)
  if (length(roots) <= 1L) return(list())
  out <- list()
  for (i in seq_along(roots)) {
    p <- roots[[i]]
    lm <- if (i == 1L) p / 2L else (roots[[i - 1L]] + p) / 2L
    rm <- if (i == length(roots)) (p + self$length) / 2L else (p + roots[[i + 1L]]) / 2L
    prod <- subs(bm, x, lm) * subs(bm, x, rm)
    if (isTRUE(ask_fact(prod, "negative"))) out[[length(out) + 1L]] <- p
  }
  out
}

# the slope and deflection of a joined beam with a piecewise second moment
.beam_joined <- function(self, sign, what) {
  x <- self$variable; E <- self$elastic_modulus
  prs <- .pw_pairs(self$second_moment)
  prev_slope <- S_Zero; prev_def <- S_Zero; prev_end <- S_Zero; out <- S_Zero
  for (i in seq_along(prs)) {
    if (i > 1L) prev_end <- prs[[i - 1L]][[2L]]$args[[2L]]
    end_i <- prs[[i]][[2L]]$args[[2L]]
    sv <- sign * S_One / E * integrate(self$bending_moment() / prs[[i]][[1L]], list(x, prev_end, x))
    if (what == "slope") {
      seg <- prev_slope + sv
    } else {
      dv <- integrate(prev_slope + sv, list(x, prev_end, x))
      seg <- prev_def + dv
      prev_def <- subs(dv, x, end_i)
    }
    out <- out + seg * SingularityFunction(x, prev_end, 0L)
    if (i != length(prs)) out <- out - seg * SingularityFunction(x, end_i, 0L)
    prev_slope <- subs(sv, x, end_i)
  }
  out
}

.beam_slope <- function(self) {
  x <- self$variable; E <- self$elastic_modulus; I <- self$second_moment
  if (!length(.cm_get(self, "_bc_slope"))) return(diff(self$deflection(), x))
  if (is_sym(I) && I$type == "Piecewise" && isTRUE(.cm_get(self, "_joined_beam"))) return(.beam_joined(self, -1L, "slope"))
  C3 <- Symbol("C3")
  sc <- -integrate(S_One / (E * I) * self$bending_moment(), x) + C3
  eqs <- lapply(.cm_get(self, "_bc_slope"), function(bc) subs(sc, x, bc[[1L]]) - bc[[2L]])
  subs(sc, C3, .cm_linsolve(eqs, list(C3))[[1L]])
}

.beam_deflection <- function(self) {
  x <- self$variable; E <- self$elastic_modulus; I <- self$second_moment
  bd <- .cm_get(self, "_bc_deflection"); bs <- .cm_get(self, "_bc_slope")
  joined <- is_sym(I) && I$type == "Piecewise" && isTRUE(.cm_get(self, "_joined_beam"))
  bc <- .cm_get(self, "_base_char")
  if (!length(bd) && !length(bs)) {
    if (joined) return(.beam_joined(self, -1L, "deflection"))
    return(S_One / (E * I) * integrate(-integrate(self$bending_moment(), x), x) +
             Symbol(paste0(bc, "3")) * x + Symbol(paste0(bc, "4")))
  }
  if (!length(bd)) return(integrate(self$slope(), x) + Symbol(paste0(bc, "4")))
  if (!length(bs)) {
    if (joined) return(.beam_joined(self, -1L, "deflection"))
    C3 <- Symbol(paste0(bc, "3")); C4 <- Symbol(paste0(bc, "4"))
    sc <- -integrate(self$bending_moment(), x) + C3
    dc <- integrate(sc, x) + C4
    eqs <- lapply(bd, function(b) subs(dc, x, b[[1L]]) - b[[2L]])
    k <- .cm_linsolve(eqs, list(C3, C4))
    dc <- subs(subs(dc, C3, k[[1L]]), C4, k[[2L]])
    return(S_One / (E * I) * dc)
  }
  if (joined) return(.beam_joined(self, 1L, "deflection"))
  C4 <- Symbol("C4")
  dc <- integrate(self$slope(), x) + C4
  eqs <- lapply(bd, function(b) subs(dc, x, b[[1L]]) - b[[2L]])
  subs(dc, C4, .cm_linsolve(eqs, list(C4))[[1L]])
}

.beam_max_deflection <- function(self) {
  x <- self$variable; l <- self$length
  sl <- self$slope()
  bps <- lapply(.cm_sf_atoms(sl), function(t) t$args[[2L]])
  bps <- Filter(function(p) isTRUE(p > 0) && isTRUE(p < l), bps)
  bps <- .beam_unique(c(list(S_Zero), bps, list(l)))
  bps <- .cm_sorted(bps)
  points <- list()
  for (i in seq_len(length(bps) - 1L)) {
    p <- .cm_on_interval(sl, x, bps[[i]])
    if (is_zero_num(p))
      .physics_error("NotImplementedError", "An expression is already zero in a region: solve can only represent discrete, not interval, solutions.")
    points <- c(points, .cm_roots_in(p, x, bps[[i]], bps[[i + 1L]], closed = i > 1L))
  }
  dc <- self$deflection()
  defl <- lapply(points, function(p) Abs(subs(dc, x, p)))
  if (!length(defl)) return(NULL)
  k <- .beam_argmax(defl)
  .pytuple(list(points[[k]], defl[[k]]))
}

.beam_join <- function(self, beam, via = "fixed") {
  x <- self$variable; E <- self$elastic_modulus
  new_length <- self$length + beam$length
  if (!identical(as_sym(self$elastic_modulus)$key, as_sym(beam$elastic_modulus)$key))
    .physics_error("NotImplementedError", "Joining beams with different Elastic modulus is not implemented.")
  I <- if (!identical(self$second_moment$key, beam$second_moment$key))
    Piecewise(list(self$second_moment, x <= self$length), list(beam$second_moment, x <= new_length))
  else self$second_moment
  if (!via %in% c("fixed", "hinge")) .physics_error("ValueError", "Invalid joining method. Choose from 'fixed' or 'hinge'.")
  nb <- Beam(new_length, E, I, x)
  .cm_set(nb, "_joined_beam", TRUE)
  if (via == "hinge") nb$apply_rotation_hinge(self$length)
  nb
}

# _solve_for_ild_equations
.beam_ild_eqs <- function(self, value) {
  x <- self$variable; a <- .cm_get(self, "ild_variable")
  load <- .cm_get(self, "_load") + value * SingularityFunction(x, a, -1L)
  sf <- -integrate(load, x)
  list(sf, integrate(sf, x))
}

.beam_ild_reactions <- function(self, value, ...) {
  reactions <- list(...)
  value <- as_sym(value)
  e <- .beam_ild_eqs(self, value)
  x <- self$variable; l <- self$length; a <- .cm_get(self, "ild_variable")
  rj <- .cm_get(self, "_rotation_hinge_symbols"); dj <- .cm_get(self, "_sliding_hinge_symbols")
  C3 <- Symbol("C3"); C4 <- Symbol("C4")
  SF <- function(p, q, r) SingularityFunction(p, q, r)
  eqs <- list(limit(e[[1L]], x, l) - value * (SF(a, 0L, 0L) - SF(a, l, 0L)),
              limit(e[[2L]], x, l) - value * (l * SF(a, 0L, 0L) - SF(a, 0L, 1L) + SF(a, l, 1L)))
  for (bc in .cm_get(self, "_bc_shear_force")) {
    pos <- bc[[1L]]
    eqs[[length(eqs) + 1L]] <- .beam_finite_terms(subs(self$shear_force(), x, pos) - bc[[2L]]) -
      value * (SF(-a, -as_sym(pos), 0L) - SF(-a, 0L, 0L))
  }
  for (bc in .cm_get(self, "_bc_bending_moment")) {
    pos <- bc[[1L]]
    eqs[[length(eqs) + 1L]] <- .beam_finite_terms(subs(self$bending_moment(), x, pos) - bc[[2L]]) -
      value * (pos * SF(a, 0L, 0L) - SF(a, 0L, 1L) + SF(a, pos, 1L))
  }
  sc <- integrate(e[[2L]], x) + C3
  for (bc in .cm_get(self, "_bc_slope"))
    eqs[[length(eqs) + 1L]] <- subs(sc, x, bc[[1L]]) - bc[[2L]] + value * (SF(-a, 0L, 1L) + bc[[1L]] * SF(-a, 0L, 0L))^2 / 2L
  dc <- integrate(sc, x) + C4
  for (bc in .cm_get(self, "_bc_deflection"))
    eqs[[length(eqs) + 1L]] <- subs(dc, x, bc[[1L]]) - bc[[2L]] + value * (SF(-a, 0L, 1L) + bc[[1L]] * SF(-a, 0L, 0L))^3 / 6L
  sol <- .cm_linsolve(eqs, c(list(C3, C4), reactions, rj, dj))
  nr <- length(reactions); nj <- length(rj)
  .cm_set(self, "_ild_reactions", .pydict(reactions, sol[2L + seq_len(nr)]))
  .cm_set(self, "_ild_rotations_jumps", .pydict(rj, sol[2L + nr + seq_len(nj)]))
  .cm_set(self, "_ild_deflection_jumps", .pydict(dj, sol[2L + nr + nj + seq_along(dj)]))
  invisible(NULL)
}

.beam_ild_shear <- function(self, distance, value, ...) {
  x <- self$variable; l <- self$length; a <- .cm_get(self, "ild_variable")
  value <- as_sym(value); distance <- as_sym(distance)
  sf <- .beam_ild_eqs(self, value)[[1L]]
  c1 <- value - limit(sf, x, distance)
  c2 <- (limit(sf, x, l) - limit(sf, x, distance)) - value
  d <- .cm_get(self, "_ild_reactions")
  for (r in list(...)) { v <- .beam_dict_get(d, r); c1 <- subs(c1, r, v); c2 <- subs(c2, r, v) }
  .cm_set(self, "_ild_shear", c1 - (c1 - c2) * SingularityFunction(a, distance, 0L) -
            value * SingularityFunction(-a, 0L, 0L) + value * SingularityFunction(a, l, 0L))
  invisible(NULL)
}

.beam_ild_moment <- function(self, distance, value, ...) {
  x <- self$variable; l <- self$length; a <- .cm_get(self, "ild_variable")
  value <- as_sym(value); distance <- as_sym(distance)
  m <- .beam_ild_eqs(self, value)[[2L]]
  SF <- function(p, q, r) SingularityFunction(p, q, r)
  c1 <- value * (distance * SF(a, 0L, 0L) - SF(a, 0L, 1L) + SF(a, distance, 1L)) - limit(m, x, distance)
  c2 <- limit(m, x, l) - limit(m, x, distance) - value * (l * SF(a, 0L, 0L) - SF(a, 0L, 1L) + SF(a, l, 1L))
  d <- .cm_get(self, "_ild_reactions")
  for (r in list(...)) { v <- .beam_dict_get(d, r); c1 <- subs(c1, r, v); c2 <- subs(c2, r, v) }
  .cm_set(self, "_ild_moment", c1 - (c1 - c2) * SingularityFunction(a, distance, 0L))
  invisible(NULL)
}

# the values of the free symbols other than the variable must be given
.beam_plot_subs <- function(e, var, subs) {
  subs <- subs %||% .pydict(list(), list())
  for (s in free_symbols(e)) {
    if (identical(s$key, var$key)) next
    if (!any(vapply(subs$keys, function(k) identical(as_sym(k)$key, s$key), TRUE)))
      .physics_error("ValueError", sprintf("Value of %s was not passed.", .str(s)))
  }
  subs
}

.beam_plot <- function(self, e, title, ylab, col, subs = NULL) {
  subs <- .beam_plot_subs(e, self$variable, subs)
  len <- .cm_dict_subs(self$length, subs)
  .cm_plot(.cm_dict_subs(e, subs), self$variable, len, title, ylab, col)
}

local({
  sym_setter <- function(field) function(self, v) .cm_set(self, field, as_sym(v))
  bc_prop <- function(field) function(self) .cm_get(self, field)
  bc_set <- function(field) function(self, v) .cm_set(self, field, as.list(v))
  .pyclass_def("Beam",
    props = list(
      reaction_loads = function(self) .cm_get(self, "_reaction_loads"),
      rotation_jumps = function(self) .cm_get(self, "_rotation_jumps"),
      deflection_jumps = function(self) .cm_get(self, "_deflection_jumps"),
      ild_shear = function(self) .cm_get(self, "_ild_shear"),
      ild_reactions = function(self) .cm_get(self, "_ild_reactions"),
      ild_rotation_jumps = function(self) .cm_get(self, "_ild_rotations_jumps"),
      ild_deflection_jumps = function(self) .cm_get(self, "_ild_deflection_jumps"),
      ild_moment = function(self) .cm_get(self, "_ild_moment"),
      length = function(self) .cm_get(self, "_length"),
      area = function(self) .cm_get(self, "_area"),
      variable = function(self) .cm_get(self, "_variable"),
      elastic_modulus = function(self) .cm_get(self, "_elastic_modulus"),
      second_moment = function(self) .cm_get(self, "_second_moment"),
      cross_section = function(self) .cm_get(self, "_cross_section"),
      boundary_conditions = function(self)
        .pydict(list("deflection", "slope", "bending_moment", "shear_force"),
                lapply(c("_bc_deflection", "_bc_slope", "_bc_bending_moment", "_bc_shear_force"), function(f) .cm_get(self, f))),
      bc_shear_force = bc_prop("_bc_shear_force"), bc_bending_moment = bc_prop("_bc_bending_moment"),
      bc_slope = bc_prop("_bc_slope"), bc_deflection = bc_prop("_bc_deflection"),
      load = function(self) .cm_get(self, "_load"),
      applied_loads = function(self) .cm_get(self, "_applied_loads")),
    setters = list(
      length = sym_setter("_length"), area = sym_setter("_area"), elastic_modulus = sym_setter("_elastic_modulus"),
      variable = function(self, v) {
        if (!(is_sym(v) && v$type == "Symbol")) .physics_error("TypeError", "The variable should be a Symbol object.")
        .cm_set(self, "_variable", v)
      },
      second_moment = function(self, i) {
        .cm_set(self, "_cross_section", NULL)
        if (.beam_is_geo(i)) .physics_error("ValueError", "To update cross-section geometry use `cross_section` attribute")
        .cm_set(self, "_second_moment", as_sym(i))
      },
      cross_section = function(self, s) {
        if (!is.null(s)) {
          s <- .cm_geo_points(s)
          .cm_set(self, "_second_moment", .cm_second_moment_of_area(s)[[1L]])
        }
        .cm_set(self, "_cross_section", s)
      },
      bc_shear_force = bc_set("_bc_shear_force"), bc_bending_moment = bc_set("_bc_bending_moment"),
      bc_slope = bc_set("_bc_slope"), bc_deflection = bc_set("_bc_deflection")),
    methods = list(
      `__str__` = function(self) {
        cs <- .cm_get(self, "_cross_section")
        sprintf("Beam(%s, %s, %s)", sstr(self$length), sstr(self$elastic_modulus),
                if (!is.null(cs)) paste(format(cs), collapse = "") else sstr(self$second_moment))
      },
      join = .beam_join,
      apply_support = .beam_apply_support,
      apply_rotation_hinge = function(self, loc)
        .beam_hinge(self, loc, "P_", -3L, "_applied_rotation_hinges", "_rotation_hinge_symbols", "_bc_bending_moment"),
      apply_sliding_hinge = function(self, loc)
        .beam_hinge(self, loc, "W_", -4L, "_applied_sliding_hinges", "_sliding_hinge_symbols", "_bc_shear_force"),
      apply_load = .beam_apply_load,
      remove_load = .beam_remove_load,
      solve_for_reaction_loads = .beam_solve,
      shear_force = function(self) -integrate(self$load, self$variable),
      max_shear_force = function(self) .beam_max_curve(self, self$shear_force(), .cm_get(self, "_load")),
      bending_moment = function(self) integrate(self$shear_force(), self$variable),
      max_bmoment = function(self) .beam_max_curve(self, self$bending_moment(), self$shear_force()),
      point_cflexure = .beam_point_cflexure,
      slope = .beam_slope,
      deflection = .beam_deflection,
      max_deflection = .beam_max_deflection,
      shear_stress = function(self) self$shear_force() / .cm_get(self, "_area"),
      solve_for_ild_reactions = .beam_ild_reactions,
      solve_for_ild_shear = .beam_ild_shear,
      solve_for_ild_moment = .beam_ild_moment,
      plot_shear_stress = function(self, subs = NULL) .beam_plot(self, self$shear_stress(), "Shear Stress", "tau", "red", subs),
      plot_shear_force = function(self, subs = NULL) .beam_plot(self, self$shear_force(), "Shear Force", "V", "green", subs),
      plot_bending_moment = function(self, subs = NULL) .beam_plot(self, self$bending_moment(), "Bending Moment", "M", "blue", subs),
      plot_slope = function(self, subs = NULL) .beam_plot(self, self$slope(), "Slope", "theta", "magenta", subs),
      plot_deflection = function(self, subs = NULL) .beam_plot(self, self$deflection(), "Deflection", "delta", "red", subs),
      plot_loading_results = function(self, subs = NULL) {
        .beam_plot_subs(self$deflection(), self$variable, subs)
        op <- graphics::par(mfrow = c(4L, 1L)); on.exit(graphics::par(op))
        self$plot_shear_force(subs); self$plot_bending_moment(subs); self$plot_slope(subs); self$plot_deflection(subs)
      },
      plot_ild_reactions = function(self, subs = NULL) {
        d <- .cm_get(self, "_ild_reactions")
        if (!length(d$keys)) .physics_error("ValueError", "I.L.D. reaction equations not found. Please use solve_for_ild_reactions() to generate the I.L.D. reaction equations.")
        a <- .cm_get(self, "ild_variable")
        for (v in d$values) .beam_plot_subs(v, a, subs)
        .beam_plot_subs(self$length, a, subs)
        for (i in seq_along(d$keys)) .cm_plot(.cm_dict_subs(d$values[[i]], subs %||% .pydict(list(), list())), a,
                                              .cm_dict_subs(self$length, subs %||% .pydict(list(), list())),
                                              "I.L.D. for Reactions", .str(d$keys[[i]]), "blue")
      },
      plot_ild_shear = function(self, subs = NULL) {
        e <- .cm_get(self, "_ild_shear")
        if (is_zero_num(e)) .physics_error("ValueError", "I.L.D. shear equation not found. Please use solve_for_ild_shear() to generate the I.L.D. shear equations.")
        a <- .cm_get(self, "ild_variable")
        subs <- .beam_plot_subs(e, a, subs); .beam_plot_subs(self$length, a, subs)
        .cm_plot(.cm_dict_subs(e, subs), a, .cm_dict_subs(self$length, subs), "I.L.D. for Shear", "V", "blue")
      },
      plot_ild_moment = function(self, subs = NULL) {
        e <- .cm_get(self, "_ild_moment")
        if (is_zero_num(e)) .physics_error("ValueError", "I.L.D. moment equation not found. Please use solve_for_ild_moment() to generate the I.L.D. moment equations.")
        a <- .cm_get(self, "ild_variable")
        subs <- .beam_plot_subs(e, a, subs); .beam_plot_subs(self$length, a, subs)
        .cm_plot(.cm_dict_subs(e, subs), a, .cm_dict_subs(self$length, subs), "I.L.D. for Moment", "M", "blue")
      }))
})

#' @rdname Beam
#' @export
Beam3D <- function(length, elastic_modulus, shear_modulus, second_moment, area, variable = Symbol("x")) {
  self <- .pyenv_new("Beam3D")
  # SymPy passes the variable to Beam() in the place of the area
  .beam_init(self, length, elastic_modulus, second_moment, variable, Symbol("x"), "C", Symbol("a"))
  self$shear_modulus <- shear_modulus
  self$area <- area
  for (f in c("_load_vector", "_moment_load_vector", "_load_Singularity", "_slope", "_deflection"))
    .cm_set(self, f, list(0L, 0L, 0L))
  .cm_set(self, "_torsion_moment", .pd_empty())
  .cm_set(self, "_angular_deflection", 0L)
  self
}

# apply_load() compares the direction as given; the other methods lower it first
.b3_dir <- function(dir) switch(dir, x = 1L, y = 2L, 3L)

.b3_add <- function(self, field, k, v) {
  l <- .cm_get(self, field)
  l[[k]] <- as_sym(l[[k]]) + v
  .cm_set(self, field, l)
}

.b3_apply_load <- function(self, value, start, order, dir = "y") {
  x <- self$variable
  value <- as_sym(value); start <- as_sym(start); order <- as_sym(order)
  k <- .b3_dir(dir)
  if (!.cm_same(order, -1L)) .b3_add(self, "_load_vector", k, value)
  .b3_add(self, "_load_Singularity", k, value * SingularityFunction(x, start, order))
  invisible(NULL)
}

.b3_apply_moment_load <- function(self, value, start, order, dir = "y") {
  x <- self$variable
  value <- as_sym(value); start <- as_sym(start); order <- as_sym(order)
  k <- .b3_dir(dir)
  if (!.cm_same(order, -2L)) {
    .b3_add(self, "_moment_load_vector", k, value)
  } else if (k == 1L) {
    tm <- .cm_get(self, "_torsion_moment")
    .cm_set(self, "_torsion_moment", .pd_set(tm, start, if (.pd_has(tm, start)) .pd_get(tm, start) + value else value))
  }
  # SymPy adds the moments of every direction to the x singularity load
  .b3_add(self, "_load_Singularity", 1L, value * SingularityFunction(x, start, order))
  invisible(NULL)
}

.b3_apply_support <- function(self, loc, type = "fixed") {
  if (!type %in% c("pin", "roller", "fixed"))
    .physics_error("ValueError", "Invalid support type. Choose from 'pin', 'roller', or 'fixed'.")
  nm <- if (is_sym(loc)) .str(loc) else .pyfmt(loc)
  R <- Symbol(paste0("R_", nm))
  zero3 <- list(0L, 0L, 0L)
  if (type %in% c("pin", "roller")) {
    .cm_set(self, "_reaction_loads", .pd_set(.cm_get(self, "_reaction_loads"), R, R))
    .cm_push(self, "_bc_deflection", .pytuple(list(loc, zero3)))
  } else {
    M <- Symbol(paste0("M_", nm))
    .cm_set(self, "_reaction_loads", .pd_set(.cm_get(self, "_reaction_loads"), R, list(R, M)))
    .cm_push(self, "_bc_deflection", .pytuple(list(loc, zero3)))
    .cm_push(self, "_bc_slope", .pytuple(list(loc, zero3)))
  }
  invisible(NULL)
}

.b3_solve <- function(self, ...) {
  reaction <- list(...)
  x <- self$variable; l <- self$length
  q <- lapply(.cm_get(self, "_load_Singularity"), as_sym)
  shear <- lapply(q, function(ld) integrate(ld, x))
  moment <- lapply(shear, function(s) integrate(s, x))
  for (i in 1:3) {
    react <- Filter(function(r) has(shear[[i]], r) || has(moment[[i]], r), reaction)
    if (!length(react)) next
    sol <- .cm_linsolve(list(limit(shear[[i]], x, l), limit(moment[[i]], x, l)), react)
    rl <- .cm_get(self, "_reaction_loads")
    for (k in seq_along(react)) {
      if (.pd_has(rl, react[[k]])) {
        old <- .pd_get(rl, react[[k]])
        if (!(is_sym(old) && .cm_same(old, sol[[k]])))
          .physics_error("ValueError", sprintf("Ambiguous solution for %s in different directions.", .str(react[[k]])))
      }
    }
    for (k in seq_along(react)) rl <- .pd_set(rl, react[[k]], sol[[k]])
    .cm_set(self, "_reaction_loads", rl)
  }
  invisible(NULL)
}

.b3_shear_force <- function(self) {
  x <- self$variable
  lapply(.cm_get(self, "_load_vector"), function(q) integrate(-as_sym(q), x))
}

.b3_bending_moment <- function(self) {
  x <- self$variable
  m <- lapply(.cm_get(self, "_moment_load_vector"), as_sym)
  sh <- self$shear_force()
  list(integrate(-m[[1L]], x), integrate(-m[[2L]] + sh[[3L]], x), integrate(-m[[3L]] - sh[[2L]], x))
}

.b3_polar_moment <- function(self) {
  I <- self$second_moment
  if (is_sym(I)) return(2L * I)
  .add(lapply(I, as_sym))
}

.b3_solve_for_torsion <- function(self) {
  x <- self$variable
  tm <- .cm_get(self, "_torsion_moment")
  total <- .add(lapply(tm$values, as_sym))
  pts <- tm$keys
  td <- Piecewise(list(total, x <= pts[[1L]]), list(0L, x >= pts[[1L]]))
  for (i in seq_along(pts)[-1L]) {
    total <- total - tm$values[[i - 1L]]
    td <- td + Piecewise(list(0L, x <= pts[[i - 1L]]), list(total, x <= pts[[i]]), list(0L, x >= pts[[i]]))
  }
  .cm_set(self, "_angular_deflection", integrate(td, x) / (self$shear_modulus * self$polar_moment()))
  invisible(NULL)
}

# the solutions of the beam equations, as dsolve() gives them, with the
# integration constants fixed by the end conditions
.b3_solve_slope_deflection <- function(self) {
  x <- self$variable; l <- self$length
  E <- self$elastic_modulus; G <- self$shear_modulus
  I <- self$second_moment
  if (is_sym(I)) { Iy <- I; Iz <- I } else { Iy <- as_sym(I[[1L]]); Iz <- as_sym(I[[2L]]) }
  A <- .cm_get(self, "_area")
  ld <- lapply(.cm_get(self, "_load_vector"), as_sym)
  mo <- lapply(.cm_get(self, "_moment_load_vector"), as_sym)
  C1 <- Symbol("C1"); C2 <- Symbol("C2"); Ci <- Symbol("C_i")
  int2 <- function(f) expand(integrate(integrate(f, x), x))
  ends <- function(e, a, b) {
    k <- .cm_linsolve(list(subs(e, x, 0L), subs(e, x, l)), list(a, b))
    subs(subs(e, a, k[[1L]]), b, k[[2L]])
  }
  ends_k <- function(e, a, b) .cm_linsolve(list(subs(e, x, 0L), subs(e, x, l)), list(a, b))
  defl <- .cm_get(self, "_deflection"); slope <- .cm_get(self, "_slope")
  def_x <- ends(int2(-ld[[1L]] / (E * A)) + C1 + C2 * x, C1, C2)
  defl[[1L]] <- def_x
  slope[[1L]] <- diff(def_x, x)
  slope_z <- ends(int2(-(integrate(-ld[[2L]], x) + Ci + mo[[3L]]) / (E * Iz)) + C1 + C2 * x, C1, C2)
  def_y <- expand(integrate((-(ld[[2L]] * x) + Ci + G * A * slope_z) / (G * A), x)) + C1
  k <- ends_k(def_y, C1, Ci)
  defl[[2L]] <- subs(subs(def_y, C1, k[[1L]]), Ci, k[[2L]])
  slope[[3L]] <- subs(slope_z, Ci, k[[2L]])
  slope_y <- ends(int2(-(integrate(ld[[3L]], x) - Ci + mo[[2L]]) / (E * Iy)) + C1 + C2 * x, C1, C2)
  def_z <- expand(integrate((-(ld[[3L]] * x) + Ci - G * A * slope_y) / (G * A), x)) + C1
  k <- ends_k(def_z, C1, Ci)
  defl[[3L]] <- subs(subs(def_z, C1, k[[1L]]), Ci, k[[2L]])
  slope[[2L]] <- subs(slope_y, Ci, k[[2L]])
  .cm_set(self, "_deflection", defl)
  .cm_set(self, "_slope", slope)
  invisible(NULL)
}

.b3_max <- function(self, curve, slope) {
  x <- self$variable; l <- self$length
  curve <- as_sym(curve)
  if (is_zero_num(curve)) return(.pytuple(list(0L, 0L)))
  slope <- expand(as_sym(slope))
  if (is_zero_num(slope))
    .physics_error("NotImplementedError", "An expression is already zero in a region: solve can only represent discrete, not interval, solutions.")
  points <- c(.cm_roots_in(slope, x, S_Zero, l), list(S_Zero, l))
  vals <- lapply(points, function(p) Abs(subs(curve, x, p)))
  k <- .beam_argmax(vals)
  .pytuple(list(points[[k]], vals[[k]]))
}

.b3_plot <- function(self, what, dir = "all", subs = NULL) {
  dirs <- if (tolower(dir) %in% c("x", "y", "z")) tolower(dir) else c("x", "y", "z")
  curves <- switch(what, shear_force = self$shear_force(), bending_moment = self$bending_moment(),
                   slope = self$slope(), shear_stress = self$shear_stress(), self$deflection())
  title <- switch(what, shear_force = "Shear Force", bending_moment = "Bending Moment", slope = "Slope",
                  shear_stress = "Shear stress", "Deflection")
  if (length(dirs) > 1L) { op <- graphics::par(mfrow = c(3L, 1L)); on.exit(graphics::par(op)) }
  for (d in dirs) {
    e <- as_sym(curves[[.b3_dir(d)]])
    sb <- .beam_plot_subs(e, self$variable, subs)
    .cm_plot(.cm_dict_subs(e, sb), self$variable, .cm_dict_subs(self$length, sb),
             sprintf("%s along %s direction", title, d), sprintf("%s(%s)", substr(what, 1L, 1L), d),
             c(x = "red", y = "green", z = "blue")[[d]])
  }
  invisible(NULL)
}

local({
  .pyclass_def("Beam3D", parent = "Beam",
    props = list(
      shear_modulus = function(self) .cm_get(self, "_shear_modulus"),
      second_moment = function(self) .cm_get(self, "_second_moment"),
      area = function(self) .cm_get(self, "_area"),
      load_vector = function(self) .cm_get(self, "_load_vector"),
      moment_load_vector = function(self) .cm_get(self, "_moment_load_vector")),
    setters = list(
      shear_modulus = function(self, e) .cm_set(self, "_shear_modulus", as_sym(e)),
      second_moment = function(self, i) .cm_set(self, "_second_moment", if (is.list(i) && !is_sym(i)) lapply(i, as_sym) else as_sym(i)),
      area = function(self, a) .cm_set(self, "_area", as_sym(a))),
    methods = list(
      polar_moment = .b3_polar_moment,
      apply_load = .b3_apply_load,
      apply_moment_load = .b3_apply_moment_load,
      apply_support = .b3_apply_support,
      solve_for_reaction_loads = .b3_solve,
      shear_force = .b3_shear_force,
      axial_force = function(self) self$shear_force()[[1L]],
      shear_stress = function(self) lapply(self$shear_force(), function(s) s / .cm_get(self, "_area")),
      axial_stress = function(self) self$axial_force() / .cm_get(self, "_area"),
      bending_moment = .b3_bending_moment,
      torsional_moment = function(self) self$bending_moment()[[1L]],
      solve_for_torsion = .b3_solve_for_torsion,
      solve_slope_deflection = .b3_solve_slope_deflection,
      slope = function(self) .cm_get(self, "_slope"),
      deflection = function(self) .cm_get(self, "_deflection"),
      angular_deflection = function(self) .cm_get(self, "_angular_deflection"),
      max_shear_force = function(self)
        lapply(1:3, function(d) .b3_max(self, self$shear_force()[[d]], .cm_get(self, "_load_vector")[[d]])),
      max_bending_moment = function(self)
        lapply(1:3, function(d) .b3_max(self, self$bending_moment()[[d]], self$shear_force()[[d]])),
      max_bmoment = function(self) self$max_bending_moment(),
      max_deflection = function(self)
        lapply(1:3, function(d) .b3_max(self, self$deflection()[[d]], self$slope()[[d]])),
      plot_shear_force = function(self, dir = "all", subs = NULL) .b3_plot(self, "shear_force", dir, subs),
      plot_bending_moment = function(self, dir = "all", subs = NULL) .b3_plot(self, "bending_moment", dir, subs),
      plot_slope = function(self, dir = "all", subs = NULL) .b3_plot(self, "slope", dir, subs),
      plot_deflection = function(self, dir = "all", subs = NULL) .b3_plot(self, "deflection", dir, subs),
      plot_shear_stress = function(self, dir = "all", subs = NULL) .b3_plot(self, "shear_stress", dir, subs)))
})
