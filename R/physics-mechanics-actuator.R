# Pathways, wrapping geometries and actuators, after
# sympy/physics/mechanics/pathway.py, wrapping_geometry.py and actuator.py.
#
# These are mutable environments whose attributes, as in SymPy, can be set
# once: a second assignment raises an AttributeError.

.pm_new_obj <- function(classes) {
  self <- new.env(parent = emptyenv())
  class(self) <- c(classes, "physics_obj")
  assign("_uid", .pv_next_uid(), envir = self)
  self
}

.pm_has <- function(self, nm) exists(nm, envir = self, inherits = FALSE)

.pm_immutable <- function(self, field, attr, value) {
  if (.pm_has(self, field))
    .pv_err("AttributeError", sprintf("Can't set attribute `%s` to %s as it is immutable.", attr, .pv_objrepr(value)))
}

.pm_sympify_strict <- function(x) {
  if (is.character(x)) .pv_err("SympifyError", sprintf("SympifyError: %s", .pv_objrepr(x)))
  if (is.null(x)) .pv_err("SympifyError", "SympifyError: None")
  as_sym(x)
}

.PM_POINT_CLASS <- "<class 'sympy.physics.vector.point.Point'>"

#' @export
format.PathwayBase <- function(x, ...) .pv_objrepr(x)
#' @export
format.WrappingGeometryBase <- function(x, ...) .pv_objrepr(x)
#' @export
format.ActuatorBase <- function(x, ...) .pv_objrepr(x)

# ---------------------------------------------------------------------------
# pathway.py

.pm_point_pair_relative_position <- function(point_1, point_2) .pv_pos_from(point_2, point_1)
.pm_point_pair_length <- function(point_1, point_2) .pv_magnitude(.pm_point_pair_relative_position(point_1, point_2))
.pm_point_pair_extension_velocity <- function(point_1, point_2) .pv_diff(.pm_point_pair_length(point_1, point_2), .pv_t())

.pm_set_attachments <- function(self, attachments, min3 = FALSE) {
  tup <- .pytuple(attachments)
  .pm_immutable(self, "_attachments", "attachments", tup)
  n <- length(attachments)
  if (!min3 && n != 2L)
    .pv_err("ValueError", sprintf("Value %s passed to `attachments` was an iterable of length %d, must be an iterable of length 2.",
                                  .pv_objrepr(tup), n))
  if (min3 && n <= 2L)
    .pv_err("ValueError", sprintf("Value %s passed to `attachments` was an iterable of length %d, must be an iterable of length 3 or greater.",
                                  .pv_objrepr(tup), n))
  for (i in seq_len(n)) {
    point <- attachments[[i]]
    if (!.pv_is_point(point))
      .pv_err("TypeError", sprintf("Value %s passed to `attachments` at index %d was of type %s, must be %s.",
                                   .pv_objrepr(point), i - 1L, .pv_pyclass(point), .PM_POINT_CLASS))
  }
  self[["_attachments"]] <- tup
}

.pvR_PathwayBase <- function(x)
  sprintf("%s(%s)", class(x)[1L], paste(vapply(x[["_attachments"]], .pv_str, ""), collapse = ", "))
.pvP_PathwayBase_attachments <- function(self) self[["_attachments"]]
.pvS_PathwayBase_attachments <- function(self, value) .pm_set_attachments(self, .pm_iterable(value))
.pvS_ObstacleSetPathway_attachments <- function(self, value) .pm_set_attachments(self, .pm_iterable(value), TRUE)

#' Pathways of actuators
#'
#' `LinearPathway()` is the straight line between two points,
#' `ObstacleSetPathway()` the piecewise straight path through three or more
#' points and `WrappingPathway()` the geodesic between two points over a
#' wrapping geometry, as in SymPy. Each exposes `attachments`, `length`,
#' `extension_velocity` and `to_loads(force)` through `$`.
#'
#' @param ... The attachment points.
#' @param attachment_1,attachment_2 The end points of a wrapping pathway.
#' @param geometry A wrapping geometry, see [WrappingSphere()].
#' @return A pathway object.
#' @examples
#' N <- ReferenceFrame("N")
#' pA <- PhysPoint("pA"); pB <- PhysPoint("pB")
#' pB$set_pos(pA, Symbol("q") * N$x)
#' path <- LinearPathway(pA, pB)
#' path$length
#' @export
LinearPathway <- function(...) {
  self <- .pm_new_obj(c("LinearPathway", "PathwayBase"))
  .pm_set_attachments(self, list(...))
  self
}

.pvP_LinearPathway_length <- function(self) do.call(.pm_point_pair_length, unclass(self[["_attachments"]]))
.pvP_LinearPathway_extension_velocity <- function(self)
  do.call(.pm_point_pair_extension_velocity, unclass(self[["_attachments"]]))
.pvM_LinearPathway_to_loads <- function(self, force) {
  a <- self[["_attachments"]]
  relative_position <- .pm_point_pair_relative_position(a[[1L]], a[[2L]])
  len <- .pvP_LinearPathway_length(self)
  force <- as_sym(force)
  list(Force(a[[1L]], .pv_vdiv(.pv_vmul(relative_position, .neg(force)), len)),
       Force(a[[length(a)]], .pv_vdiv(.pv_vmul(relative_position, force), len)))
}

#' @rdname LinearPathway
#' @export
ObstacleSetPathway <- function(...) {
  self <- .pm_new_obj(c("ObstacleSetPathway", "PathwayBase"))
  .pm_set_attachments(self, list(...), TRUE)
  self
}

.pm_pairs <- function(a) lapply(seq_len(length(a) - 1L), function(i) list(a[[i]], a[[i + 1L]]))

.pvP_ObstacleSetPathway_length <- function(self) {
  out <- S_Zero
  for (p in .pm_pairs(self[["_attachments"]])) out <- .add(list(out, .pm_point_pair_length(p[[1L]], p[[2L]])))
  out
}
.pvP_ObstacleSetPathway_extension_velocity <- function(self) {
  out <- S_Zero
  for (p in .pm_pairs(self[["_attachments"]])) out <- .add(list(out, .pm_point_pair_extension_velocity(p[[1L]], p[[2L]])))
  out
}
.pvM_ObstacleSetPathway_to_loads <- function(self, force) {
  force <- as_sym(force)
  loads <- list()
  for (p in .pm_pairs(self[["_attachments"]])) {
    relative_position <- .pm_point_pair_relative_position(p[[1L]], p[[2L]])
    len <- .pm_point_pair_length(p[[1L]], p[[2L]])
    loads <- c(loads, list(Force(p[[1L]], .pv_vdiv(.pv_vmul(relative_position, .neg(force)), len)),
                           Force(p[[2L]], .pv_vdiv(.pv_vmul(relative_position, force), len))))
  }
  loads
}

#' @rdname LinearPathway
#' @export
WrappingPathway <- function(attachment_1, attachment_2, geometry) {
  self <- .pm_new_obj(c("WrappingPathway", "PathwayBase"))
  .pm_set_attachments(self, list(attachment_1, attachment_2))
  .pvS_WrappingPathway_geometry(self, geometry)
  self
}

.PM_WRAP_CLASS <- "<class 'sympy.physics.mechanics.wrapping_geometry.WrappingGeometryBase'>"

.pvP_WrappingPathway_geometry <- function(self) self[["_geometry"]]
.pvS_WrappingPathway_geometry <- function(self, geometry) {
  .pm_immutable(self, "_geometry", "geometry", geometry)
  if (!inherits(geometry, "WrappingGeometryBase"))
    .pv_err("TypeError", sprintf("Value %s passed to `geometry` was of type %s, must be %s.",
                                 .pv_objrepr(geometry), .pv_pyclass(geometry), .PM_WRAP_CLASS))
  self[["_geometry"]] <- geometry
}
.pvP_WrappingPathway_length <- function(self) {
  a <- self[["_attachments"]]
  self[["_geometry"]]$geodesic_length(a[[1L]], a[[2L]])
}
.pvP_WrappingPathway_extension_velocity <- function(self) .pv_diff(self$length, .pv_t())
.pvM_WrappingPathway_to_loads <- function(self, force) {
  force <- as_sym(force)
  a <- self[["_attachments"]]
  pA <- a[[1L]]; pB <- a[[2L]]
  pO <- self[["_geometry"]]$point
  ends <- self[["_geometry"]]$geodesic_end_vectors(pA, pB)
  pO_force <- .pv_vneg(.pv_vadd(ends[[1L]], ends[[2L]]))
  list(Force(pA, .pv_vmul(ends[[1L]], force)), Force(pB, .pv_vmul(ends[[2L]], force)),
       Force(pO, .pv_vmul(pO_force, force)))
}
.pvR_WrappingPathway <- function(x)
  sprintf("%s(%s, geometry=%s)", class(x)[1L], paste(vapply(x[["_attachments"]], .pv_str, ""), collapse = ", "),
          .pv_objrepr(x[["_geometry"]]))

# ---------------------------------------------------------------------------
# wrapping_geometry.py

#' Wrapping geometries
#'
#' `WrappingSphere()`, `WrappingCylinder()` and `WrappingCone()` are the
#' surfaces over which a [WrappingPathway()] wraps, as in SymPy. Each
#' exposes `point`, `point_on_surface(point)`, `geodesic_length(p1, p2)`
#' and `geodesic_end_vectors(p1, p2)` through `$`.
#'
#' @param radius The radius.
#' @param point The center of a sphere, or a point on the axis of a
#'   cylinder.
#' @param axis The axis vector.
#' @param alpha The half angle of a cone.
#' @param apex The apex of a cone.
#' @return A wrapping geometry.
#' @examples
#' N <- ReferenceFrame("N")
#' O <- PhysPoint("O")
#' WrappingCylinder(Symbol("r"), O, N$x)
#' @export
WrappingSphere <- function(radius, point) {
  self <- .pm_new_obj(c("WrappingSphere", "WrappingGeometryBase"))
  self[["_radius"]] <- radius
  self[["_point"]] <- point
  self
}

.pvR_WrappingGeometryBase <- function(x) sprintf("%s()", class(x)[1L])
.pvR_WrappingSphere <- function(x) sprintf("%s(radius=%s, point=%s)", class(x)[1L], .pv_str(x[["_radius"]]), .pv_str(x[["_point"]]))
.pvP_WrappingSphere_radius <- function(self) self[["_radius"]]
.pvS_WrappingSphere_radius <- function(self, radius) self[["_radius"]] <- radius
.pvP_WrappingSphere_point <- function(self) self[["_point"]]
.pvS_WrappingSphere_point <- function(self, point) self[["_point"]] <- point

.pm_radius_sq <- function(point_vector) {
  if (.pv_is_vec(point_vector)) return(.pv_vdot(point_vector, point_vector))
  .pow(as_sym(point_vector), S_Two)
}

.pvM_WrappingSphere_point_on_surface <- function(self, point) {
  point_radius_squared <- .pm_radius_sq(.pv_pos_from(point, self[["_point"]]))
  Eq(point_radius_squared, .pow(as_sym(self[["_radius"]]), S_Two), evaluate = FALSE)
}
.pvM_WrappingSphere_geodesic_length <- function(self, point_1, point_2) {
  v1 <- .pv_normalize(.pv_pos_from(point_1, self[["_point"]]))
  v2 <- .pv_normalize(.pv_pos_from(point_2, self[["_point"]]))
  central_angle <- fcall("acos", list(.pv_vdot(v2, v1)))
  .mul(list(as_sym(self[["_radius"]]), central_angle))
}
.pvM_WrappingSphere_geodesic_end_vectors <- function(self, point_1, point_2) {
  pA <- point_1; pB <- point_2; pO <- self[["_point"]]
  pA_vec <- .pv_pos_from(pA, pO)
  pB_vec <- .pv_pos_from(pB, pO)
  if (isTRUE(.pv_vec_eq(.pv_vcross(pA_vec, pB_vec), 0)))
    .pv_err("ValueError", sprintf(paste("Can't compute geodesic end vectors for the pair of points %s and %s on a",
                                        "sphere %s as they are diametrically opposed, thus the geodesic is not defined."),
                                  .pv_str(pA), .pv_str(pB), .pv_objrepr(self)))
  .pytuple(list(.pv_normalize(.pv_vcross(.pv_vcross(pA_vec, .pv_pos_from(pB, pA)), pA_vec)),
                .pv_normalize(.pv_vcross(.pv_vcross(pB_vec, .pv_pos_from(pA, pB)), pB_vec))))
}

#' @rdname WrappingSphere
#' @export
WrappingCylinder <- function(radius, point, axis) {
  self <- .pm_new_obj(c("WrappingCylinder", "WrappingGeometryBase"))
  self[["_radius"]] <- radius
  self[["_point"]] <- point
  self[["_axis"]] <- .pv_normalize(axis)
  self
}

.pvR_WrappingCylinder <- function(x) sprintf("%s(radius=%s, point=%s, axis=%s)", class(x)[1L], .pv_str(x[["_radius"]]),
                                             .pv_str(x[["_point"]]), .pv_str(x[["_axis"]]))
.pvP_WrappingCylinder_radius <- function(self) self[["_radius"]]
.pvS_WrappingCylinder_radius <- function(self, radius) self[["_radius"]] <- radius
.pvP_WrappingCylinder_point <- function(self) self[["_point"]]
.pvS_WrappingCylinder_point <- function(self, point) self[["_point"]] <- point
.pvP_WrappingCylinder_axis <- function(self) self[["_axis"]]
.pvS_WrappingCylinder_axis <- function(self, axis) self[["_axis"]] <- .pv_normalize(axis)

.pm_vscale <- function(v, s) .pv_vmul(v, as_sym(s))

.pvM_WrappingCylinder_point_on_surface <- function(self, point) {
  ax <- self[["_axis"]]
  relative_position <- .pv_pos_from(point, self[["_point"]])
  parallel <- .pm_vscale(ax, .pv_vdot(relative_position, ax))
  point_vector <- .pv_vsub(relative_position, parallel)
  Eq(trigsimp(.pm_radius_sq(point_vector)), .pow(as_sym(self[["_radius"]]), S_Two), evaluate = FALSE)
}

.pm_perp_unit <- function(rel, ax) .pv_normalize(.pv_vsub(rel, .pm_vscale(ax, .pv_vdot(rel, ax))))

.pvM_WrappingCylinder_geodesic_length <- function(self, point_1, point_2) {
  ax <- self[["_axis"]]
  relative_position <- .pv_pos_from(point_2, point_1)
  parallel_length <- .pv_vdot(relative_position, ax)
  p1 <- .pm_perp_unit(.pv_pos_from(point_1, self[["_point"]]), ax)
  p2 <- .pm_perp_unit(.pv_pos_from(point_2, self[["_point"]]), ax)
  central_angle <- .pm_directional_atan(cancel(.pv_vdot(.pv_vcross(p1, p2), ax)), cancel(.pv_vdot(p1, p2)))
  planar_arc_length <- .mul(list(as_sym(self[["_radius"]]), central_angle))
  .sqrt(.add(list(.pow(parallel_length, S_Two), .pow(planar_arc_length, S_Two))))
}

.pvM_WrappingCylinder_geodesic_end_vectors <- function(self, point_1, point_2) {
  ax <- self[["_axis"]]
  r1 <- .pv_pos_from(point_1, self[["_point"]])
  r2 <- .pv_pos_from(point_2, self[["_point"]])
  if (isTRUE(.pv_vec_eq(r1, r2)))
    .pv_err("ValueError", sprintf("Cannot compute geodesic end vectors for coincident points %s and %s as no geodesic exists.",
                                  .pv_str(point_1), .pv_str(point_2)))
  n1 <- .pv_vsub(r1, .pm_vscale(ax, .pv_vdot(r1, ax)))
  n2 <- .pv_vsub(r2, .pm_vscale(ax, .pv_vdot(r2, ax)))
  if (isTRUE(.pv_vec_eq(n1, n2))) {
    perp1 <- .pv_vec0()
    perp2 <- .pv_vec0()
  } else {
    perp1 <- .pv_normalize(.pv_vcross(ax, n1))
    perp2 <- .pv_vneg(.pv_normalize(.pv_vcross(ax, n2)))
  }
  geodesic_length <- .pvM_WrappingCylinder_geodesic_length(self, point_1, point_2)
  parallel_length <- .pv_vdot(.pv_pos_from(point_2, point_1), ax)
  planar_arc_length <- .sqrt(.sub(.pow(geodesic_length, S_Two), .pow(parallel_length, S_Two)))
  v1 <- .pv_normalize(.pv_vadd(.pm_vscale(perp1, planar_arc_length), .pm_vscale(ax, parallel_length)))
  v2 <- .pv_normalize(.pv_vsub(.pm_vscale(perp2, planar_arc_length), .pm_vscale(ax, parallel_length)))
  .pytuple(list(v1, v2))
}

#' @rdname WrappingSphere
#' @export
WrappingCone <- function(alpha, apex, axis) {
  alpha <- as_sym(alpha)
  if (.expr_is_number(alpha)) {
    if (identical(alpha$key, S_Zero$key)) .pv_err("ValueError", "Cone angle alpha must be positive.")
    if (identical(alpha$key, .div(S_Pi, S_Two)$key)) .pv_err("ValueError", "Cone angle alpha must be less than pi/2.")
  } else if (identical(ask_fact(alpha, "real"), FALSE) || identical(ask_fact(alpha, "positive"), FALSE)) {
    .pv_err("ValueError", "Cone angle alpha must be real and positive.")
  }
  if (!.pv_is_point(apex)) .pv_err("TypeError", "The 'apex' must be a Point object.")
  if (!.pv_is_vec(axis)) .pv_err("TypeError", "The 'axis' must be a Vector object.")
  self <- .pm_new_obj(c("WrappingCone", "WrappingGeometryBase"))
  self[["_alpha"]] <- alpha
  self[["_apex"]] <- apex
  self[["_axis"]] <- .pv_normalize(axis)
  self
}

.pvR_WrappingCone <- function(x) sprintf("%s(alpha=%s, apex=%s, axis=%s)", class(x)[1L], .pv_str(x[["_alpha"]]),
                                         .pv_str(x[["_apex"]]), .pv_str(x[["_axis"]]))
.pvP_WrappingCone_point <- function(self) self[["_apex"]]
.pvP_WrappingCone_alpha <- function(self) self[["_alpha"]]
.pvS_WrappingCone_alpha <- function(self, alpha) self[["_alpha"]] <- alpha
.pvP_WrappingCone_apex <- function(self) self[["_apex"]]
.pvS_WrappingCone_apex <- function(self, apex) self[["_apex"]] <- apex
.pvP_WrappingCone_axis <- function(self) self[["_axis"]]
.pvS_WrappingCone_axis <- function(self, axis) self[["_axis"]] <- .pv_normalize(axis)

.pvM_WrappingCone_point_on_surface <- function(self, point) {
  ax <- self[["_axis"]]
  position <- .pv_pos_from(point, self[["_apex"]])
  axis_component <- .pm_vscale(ax, .pv_vdot(position, ax))
  radial_component <- .pv_vsub(position, axis_component)
  lhs <- .pv_vdot(radial_component, radial_component)
  rhs <- .mul(list(.pv_vdot(axis_component, axis_component), .pow(fcall("tan", list(self[["_alpha"]])), S_Two)))
  Eq(lhs, rhs, evaluate = FALSE)
}

.pm_cone_central <- function(n1, n2, ax) {
  central <- .pm_directional_atan(cancel(.pv_vdot(.pv_vcross(n1, n2), ax)), cancel(.pv_vdot(n1, n2)))
  Piecewise(list(central, central <= S_Pi), list(.sub(.mul(list(S_Two, S_Pi)), central), TRUE))
}

.pvM_WrappingCone_geodesic_length <- function(self, point_1, point_2) {
  ax <- self[["_axis"]]; alpha <- self[["_alpha"]]
  pos1 <- .pv_pos_from(point_1, self[["_apex"]])
  pos2 <- .pv_pos_from(point_2, self[["_apex"]])
  z1 <- .pv_vdot(pos1, ax); z2 <- .pv_vdot(pos2, ax)
  s1 <- .div(z1, cos(alpha)); s2 <- .div(z2, cos(alpha))
  n1 <- .pv_normalize(.pv_vsub(pos1, .pm_vscale(ax, z1)))
  n2 <- .pv_normalize(.pv_vsub(pos2, .pm_vscale(ax, z2)))
  delta_u <- .mul(list(.pm_cone_central(n1, n2, ax), sin(alpha)))
  .sqrt(.add(list(.pow(s1, S_Two), .pow(s2, S_Two), .neg(.mul(list(S_Two, s1, s2, cos(delta_u)))))))
}

.pvM_WrappingCone_geodesic_end_vectors <- function(self, point_1, point_2) {
  ax <- self[["_axis"]]; alpha <- self[["_alpha"]]
  pos1 <- .pv_pos_from(point_1, self[["_apex"]])
  pos2 <- .pv_pos_from(point_2, self[["_apex"]])
  if (isTRUE(.pv_vec_eq(pos1, pos2)))
    .pv_err("ValueError", sprintf("No unique geodesic exists for coincident points %s and %s.", .pv_str(point_1), .pv_str(point_2)))
  if (identical(.pv_magnitude(pos1)$key, S_Zero$key)) {
    v <- .pv_normalize(pos2)
    return(.pytuple(list(v, .pv_vneg(v))))
  }
  if (identical(.pv_magnitude(pos2)$key, S_Zero$key)) {
    v <- .pv_normalize(pos1)
    return(.pytuple(list(.pv_vneg(v), v)))
  }
  z1 <- .pv_vdot(pos1, ax); z2 <- .pv_vdot(pos2, ax)
  s1 <- .div(z1, cos(alpha)); s2 <- .div(z2, cos(alpha))
  L <- .pvM_WrappingCone_geodesic_length(self, point_1, point_2)
  if (identical(L$key, S_Zero$key)) return(.pytuple(list(.pv_vec0(), .pv_vec0())))
  n1 <- .pv_normalize(.pv_vsub(pos1, .pm_vscale(ax, z1)))
  n2 <- .pv_normalize(.pv_vsub(pos2, .pm_vscale(ax, z2)))
  delta_u <- .mul(list(.pm_cone_central(n1, n2, ax), sin(alpha)))
  g1 <- .pv_normalize(pos1); c1 <- .pv_vcross(ax, n1)
  g2 <- .pv_normalize(pos2); c2 <- .pv_vcross(ax, n2)
  v1 <- .pv_vadd(.pm_vscale(g1, .div(.sub(.mul(list(s2, cos(delta_u))), s1), L)),
                 .pm_vscale(c1, .div(.mul(list(s2, sin(delta_u))), L)))
  v2 <- .pv_vadd(.pm_vscale(g2, .div(.sub(.mul(list(s1, cos(delta_u))), s2), L)),
                 .pm_vscale(c2, .div(.mul(list(.neg(s1), sin(delta_u))), L)))
  .pytuple(list(v1, v2))
}

.pm_directional_atan <- function(numerator, denominator) {
  numerator <- as_sym(numerator); denominator <- as_sym(denominator)
  nn <- .expr_is_number(numerator); dn <- .expr_is_number(denominator)
  if (nn && dn) {
    angle <- atan2(numerator, denominator)
    if (isTRUE(.pm_lt0(angle))) angle <- .add(list(angle, .mul(list(S_Two, S_Pi))))
    return(angle)
  }
  if (nn) .pv_err("NotImplementedError", sprintf(paste("Cannot compute a directional atan when the numerator %s is",
                                                     "numeric and the denominator %s is symbolic."),
                                               .pv_str(numerator), .pv_str(denominator)))
  if (dn) .pv_err("NotImplementedError", sprintf(paste("Cannot compute a directional atan when the numerator %s is",
                                                     "symbolic and the denominator %s is numeric."),
                                               .pv_str(numerator), .pv_str(denominator)))
  ratio <- as_sym(trigsimp(.div(numerator, denominator)))
  if (ratio$type == "Function" && identical(ratio$name, "tan")) return(ratio$args[[1L]])
  if (ratio$type == "Mul" && identical(ratio$args[[1L]]$key, S_NegativeOne$key) &&
      ratio$args[[2L]]$type == "Function" && identical(ratio$args[[2L]]$name, "tan"))
    return(.sub(.mul(list(S_Two, S_Pi)), ratio$args[[2L]]$args[[1L]]))
  .pv_err("NotImplementedError", sprintf("Cannot compute a directional atan for the value %s.", .pv_str(ratio)))
}

# ---------------------------------------------------------------------------
# actuator.py

.PM_PATHWAY_CLASS <- "<class 'sympy.physics.mechanics.pathway.PathwayBase'>"

.pvM_ActuatorBase___repr__ <- function(self) .pv_objrepr(self)
.pvR_ActuatorBase <- function(x) sprintf("%s()", class(x)[1L])

#' Actuators
#'
#' The actuators of SymPy's `sympy.physics.mechanics`: `ForceActuator()`
#' applies a force along a pathway, `LinearSpring()`, `LinearDamper()`,
#' `DuffingSpring()` and `CoulombKineticFriction()` compute that force from
#' the pathway's length or extension velocity, and `TorqueActuator()`
#' applies a torque about an axis between two frames.
#' `TorqueActuator_at_pin_joint()` is SymPy's
#' `TorqueActuator.at_pin_joint()`. Each exposes `to_loads()` through `$`.
#'
#' @param force The force magnitude.
#' @param pathway A pathway, see [LinearPathway()].
#' @param stiffness,damping Spring stiffness and damping coefficient.
#' @param equilibrium_length Spring length at rest.
#' @param linear_stiffness,nonlinear_stiffness Coefficients of a Duffing
#'   spring.
#' @param torque The torque magnitude.
#' @param axis The torque axis.
#' @param target_frame,reaction_frame Frames (or rigid bodies) the torque
#'   and its reaction act on.
#' @param pin_joint A [PinJoint()].
#' @param mu_k,mu_s Kinetic and static friction coefficients.
#' @param f_n The normal force.
#' @param v_s The Stribeck velocity.
#' @param sigma The viscous friction coefficient.
#' @return An actuator object.
#' @examples
#' N <- ReferenceFrame("N")
#' pA <- PhysPoint("pA"); pB <- PhysPoint("pB")
#' pB$set_pos(pA, Symbol("q") * N$x)
#' spring <- LinearSpring(Symbol("k"), LinearPathway(pA, pB))
#' spring$force
#' spring$to_loads()
#' @export
ForceActuator <- function(force, pathway) {
  self <- .pm_new_obj(c("ForceActuator", "ActuatorBase"))
  .pvS_ForceActuator_force(self, force)
  .pvS_ForceActuator_pathway(self, pathway)
  self
}

.pvP_ForceActuator_force <- function(self) self[["_force"]]
.pvS_ForceActuator_force <- function(self, force) {
  .pm_immutable(self, "_force", "force", force)
  self[["_force"]] <- .pm_sympify_strict(force)
}
.pvP_ForceActuator_pathway <- function(self) self[["_pathway"]]
.pvS_ForceActuator_pathway <- function(self, pathway) {
  .pm_immutable(self, "_pathway", "pathway", pathway)
  if (!inherits(pathway, "PathwayBase"))
    .pv_err("TypeError", sprintf("Value %s passed to `pathway` was of type %s, must be %s.",
                                 .pv_objrepr(pathway), .pv_pyclass(pathway), .PM_PATHWAY_CLASS))
  self[["_pathway"]] <- pathway
}
.pvM_ForceActuator_to_loads <- function(self) self[["_pathway"]]$to_loads(self$force)
.pvR_ForceActuator <- function(x) sprintf("%s(%s, %s)", class(x)[1L], .pv_str(x$force), .pv_objrepr(x[["_pathway"]]))

.pm_computed_force <- function(self, force) .pv_err("AttributeError", "Can't set computed attribute `force`.")

.pm_strict_setter <- function(field, attr) function(self, value) {
  .pm_immutable(self, field, attr, value)
  self[[field]] <- .pm_sympify_strict(value)
}

#' @rdname ForceActuator
#' @export
LinearSpring <- function(stiffness, pathway, equilibrium_length = 0) {
  self <- .pm_new_obj(c("LinearSpring", "ForceActuator", "ActuatorBase"))
  .pvS_LinearSpring_stiffness(self, stiffness)
  .pvS_ForceActuator_pathway(self, pathway)
  .pvS_LinearSpring_equilibrium_length(self, equilibrium_length)
  self
}
.pvP_LinearSpring_force <- function(self)
  .neg(.mul(list(self[["_stiffness"]], .sub(self[["_pathway"]]$length, self[["_equilibrium_length"]]))))
.pvS_LinearSpring_force <- .pm_computed_force
.pvP_LinearSpring_stiffness <- function(self) self[["_stiffness"]]
.pvS_LinearSpring_stiffness <- .pm_strict_setter("_stiffness", "stiffness")
.pvP_LinearSpring_equilibrium_length <- function(self) self[["_equilibrium_length"]]
.pvS_LinearSpring_equilibrium_length <- .pm_strict_setter("_equilibrium_length", "equilibrium_length")
.pvR_LinearSpring <- function(x) {
  s <- sprintf("%s(%s, %s", class(x)[1L], .pv_str(x[["_stiffness"]]), .pv_objrepr(x[["_pathway"]]))
  if (identical(x[["_equilibrium_length"]]$key, S_Zero$key)) paste0(s, ")")
  else sprintf("%s, equilibrium_length=%s)", s, .pv_str(x[["_equilibrium_length"]]))
}

#' @rdname ForceActuator
#' @export
LinearDamper <- function(damping, pathway) {
  self <- .pm_new_obj(c("LinearDamper", "ForceActuator", "ActuatorBase"))
  .pvS_LinearDamper_damping(self, damping)
  .pvS_ForceActuator_pathway(self, pathway)
  self
}
.pvP_LinearDamper_force <- function(self) .neg(.mul(list(self[["_damping"]], self[["_pathway"]]$extension_velocity)))
.pvS_LinearDamper_force <- .pm_computed_force
.pvP_LinearDamper_damping <- function(self) self[["_damping"]]
.pvS_LinearDamper_damping <- .pm_strict_setter("_damping", "damping")
.pvR_LinearDamper <- function(x) sprintf("%s(%s, %s)", class(x)[1L], .pv_str(x[["_damping"]]), .pv_objrepr(x[["_pathway"]]))

#' @rdname ForceActuator
#' @export
TorqueActuator <- function(torque, axis, target_frame, reaction_frame = NULL) {
  self <- .pm_new_obj(c("TorqueActuator", "ActuatorBase"))
  .pvS_TorqueActuator_torque(self, torque)
  .pvS_TorqueActuator_axis(self, axis)
  .pvS_TorqueActuator_target_frame(self, target_frame)
  .pvS_TorqueActuator_reaction_frame(self, reaction_frame)
  self
}

#' @rdname ForceActuator
#' @export
TorqueActuator_at_pin_joint <- function(torque, pin_joint) {
  if (!inherits(pin_joint, "PinJoint"))
    .pv_err("TypeError", sprintf("Value %s passed to `pin_joint` was of type %s, must be %s.", .pv_objrepr(pin_joint),
                                 .pv_pyclass(pin_joint), "<class 'sympy.physics.mechanics.joint.PinJoint'>"))
  TorqueActuator(torque, pin_joint$joint_axis, pin_joint$child_interframe, pin_joint$parent_interframe)
}
.pvM_TorqueActuator_at_pin_joint <- function(self, torque, pin_joint) TorqueActuator_at_pin_joint(torque, pin_joint)

.PM_VECTOR_CLASS <- "<class 'sympy.physics.vector.vector.Vector'>"

.pvP_TorqueActuator_torque <- function(self) self[["_torque"]]
.pvS_TorqueActuator_torque <- .pm_strict_setter("_torque", "torque")
.pvP_TorqueActuator_axis <- function(self) self[["_axis"]]
.pvS_TorqueActuator_axis <- function(self, axis) {
  .pm_immutable(self, "_axis", "axis", axis)
  if (!.pv_is_vec(axis))
    .pv_err("TypeError", sprintf("Value %s passed to `axis` was of type %s, must be %s.", .pv_objrepr(axis),
                                 .pv_pyclass(axis), .PM_VECTOR_CLASS))
  self[["_axis"]] <- axis
}
.pvP_TorqueActuator_target_frame <- function(self) self[["_target_frame"]]
.pvS_TorqueActuator_target_frame <- function(self, target_frame) {
  .pm_immutable(self, "_target_frame", "target_frame", target_frame)
  if (inherits(target_frame, "RigidBody")) target_frame <- target_frame$frame
  else if (!.pv_is_frame(target_frame))
    .pv_err("TypeError", sprintf("Value %s passed to `target_frame` was of type %s, must be %s.", .pv_objrepr(target_frame),
                                 .pv_pyclass(target_frame), .PV_FRAME_CLASS))
  self[["_target_frame"]] <- target_frame
}
.pvP_TorqueActuator_reaction_frame <- function(self) self[["_reaction_frame"]]
.pvS_TorqueActuator_reaction_frame <- function(self, reaction_frame) {
  if (.pm_has(self, "_reaction_frame")) .pm_immutable(self, "_reaction_frame", "reaction_frame", reaction_frame)
  if (inherits(reaction_frame, "RigidBody")) reaction_frame <- reaction_frame$frame
  else if (!.pv_is_frame(reaction_frame) && !is.null(reaction_frame))
    .pv_err("TypeError", sprintf("Value %s passed to `reaction_frame` was of type %s, must be %s.", .pv_objrepr(reaction_frame),
                                 .pv_pyclass(reaction_frame), .PV_FRAME_CLASS))
  assign("_reaction_frame", reaction_frame, envir = self)
}
.pvM_TorqueActuator_to_loads <- function(self) {
  vec <- .pv_vmul(self[["_axis"]], self[["_torque"]])
  loads <- list(Torque(self[["_target_frame"]], vec))
  if (!is.null(self[["_reaction_frame"]])) loads[[2L]] <- Torque(self[["_reaction_frame"]], .pv_vneg(vec))
  loads
}
.pvR_TorqueActuator <- function(x) {
  s <- sprintf("%s(%s, axis=%s, target_frame=%s", class(x)[1L], .pv_str(x[["_torque"]]), .pv_str(x[["_axis"]]),
               .pv_str(x[["_target_frame"]]))
  if (!is.null(x[["_reaction_frame"]])) sprintf("%s, reaction_frame=%s)", s, .pv_str(x[["_reaction_frame"]]))
  else paste0(s, ")")
}

#' @rdname ForceActuator
#' @export
DuffingSpring <- function(linear_stiffness, nonlinear_stiffness, pathway, equilibrium_length = 0) {
  self <- .pm_new_obj(c("DuffingSpring", "ForceActuator", "ActuatorBase"))
  .pvS_DuffingSpring_linear_stiffness(self, .pm_sympify_strict(linear_stiffness))
  .pvS_DuffingSpring_nonlinear_stiffness(self, .pm_sympify_strict(nonlinear_stiffness))
  .pvS_DuffingSpring_equilibrium_length(self, .pm_sympify_strict(equilibrium_length))
  if (!inherits(pathway, "PathwayBase")) .pv_err("TypeError", "pathway must be an instance of PathwayBase.")
  self[["_pathway"]] <- pathway
  self
}
.pvP_DuffingSpring_linear_stiffness <- function(self) self[["_linear_stiffness"]]
.pvS_DuffingSpring_linear_stiffness <- .pm_strict_setter("_linear_stiffness", "linear_stiffness")
.pvP_DuffingSpring_nonlinear_stiffness <- function(self) self[["_nonlinear_stiffness"]]
.pvS_DuffingSpring_nonlinear_stiffness <- .pm_strict_setter("_nonlinear_stiffness", "nonlinear_stiffness")
.pvP_DuffingSpring_equilibrium_length <- function(self) self[["_equilibrium_length"]]
.pvS_DuffingSpring_equilibrium_length <- .pm_strict_setter("_equilibrium_length", "equilibrium_length")
.pvP_DuffingSpring_force <- function(self) {
  displacement <- .sub(self[["_pathway"]]$length, self[["_equilibrium_length"]])
  .sub(.neg(.mul(list(self[["_linear_stiffness"]], displacement))),
       .mul(list(self[["_nonlinear_stiffness"]], .pow(displacement, as_sym(3L)))))
}
.pvS_DuffingSpring_force <- .pm_strict_setter("_force", "force")
.pvR_DuffingSpring <- function(x)
  sprintf("%s(%s, %s, %s, equilibrium_length=%s)", class(x)[1L], .pv_str(x[["_linear_stiffness"]]),
          .pv_str(x[["_nonlinear_stiffness"]]), .pv_objrepr(x[["_pathway"]]), .pv_str(x[["_equilibrium_length"]]))

#' @rdname ForceActuator
#' @export
CoulombKineticFriction <- function(mu_k, f_n, pathway, v_s = NULL, sigma = NULL, mu_s = NULL) {
  self <- .pm_new_obj(c("CoulombKineticFriction", "ForceActuator", "ActuatorBase"))
  self[["_mu_k"]] <- if (!is.null(mu_k)) .pm_sympify_strict(mu_k) else S_One
  self[["_mu_s"]] <- if (!is.null(mu_s)) .pm_sympify_strict(mu_s) else self[["_mu_k"]]
  self[["_f_n"]] <- .pm_sympify_strict(f_n)
  self[["_sigma"]] <- if (!is.null(sigma)) .pm_sympify_strict(sigma) else S_Zero
  self[["_v_s"]] <- if (!is.null(v_s)) .pm_sympify_strict(v_s) else 0.01
  .pvS_ForceActuator_pathway(self, pathway)
  self
}
.pvP_CoulombKineticFriction_mu_k <- function(self) self[["_mu_k"]]
.pvP_CoulombKineticFriction_mu_s <- function(self) self[["_mu_s"]]
.pvP_CoulombKineticFriction_f_n <- function(self) self[["_f_n"]]
.pvP_CoulombKineticFriction_sigma <- function(self) self[["_sigma"]]
.pvP_CoulombKineticFriction_v_s <- function(self) self[["_v_s"]]
.pvP_CoulombKineticFriction_force <- function(self) {
  v <- self[["_pathway"]]$extension_velocity
  f_c <- .mul(list(self[["_mu_k"]], self[["_f_n"]]))
  f_max <- .mul(list(self[["_mu_s"]], self[["_f_n"]]))
  stribeck_term <- .mul(list(.sub(f_max, f_c), exp(.neg(.pow(.div(v, as_sym(self[["_v_s"]])), S_Two)))))
  viscous_term <- .mul(list(self[["_sigma"]], v))
  .sub(.mul(list(.add(list(f_c, stribeck_term)), .neg(fcall("sign", list(v))))), viscous_term)
}
.pvS_CoulombKineticFriction_force <- .pm_computed_force
.pvR_CoulombKineticFriction <- function(x)
  sprintf("%s(%s, %s %s, %s, %s, %s)", class(x)[1L], .pv_str(x[["_mu_k"]]), .pv_str(x[["_mu_s"]]), .pv_str(x[["_f_n"]]),
          .pv_objrepr(x[["_pathway"]]), .pv_str(x[["_v_s"]]), .pv_str(x[["_sigma"]]))
