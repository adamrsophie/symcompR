# Bodies, inertia and loads of classical mechanics, after
# sympy/physics/mechanics/body_base.py, particle.py, rigidbody.py,
# inertia.py and loads.py.
#
# Particles and rigid bodies are mutable environments, as frames and points
# are. Inertia, Force and Torque are SymPy named tuples; here they are
# classed lists with named fields, so `I$dyadic`, `I[[1]]` and `length(I)`
# behave as in SymPy.

.pm_is_body <- function(x) inherits(x, "BodyBase")

# Python's repr() of a value
.pv_objrepr <- function(x) {
  if (is.character(x) && length(x) == 1L) return(paste0("'", x, "'"))
  if (inherits(x, "physics_obj")) {
    for (cl in class(x)) {
      f <- get0(paste0(".pvR_", cl), envir = .pv_ns, inherits = FALSE)
      if (!is.null(f)) return(f(x))
    }
    return(format(x))
  }
  if (inherits(x, "pytuple")) {
    s <- vapply(x, .pv_objrepr, "")
    return(paste0("(", paste(s, collapse = ", "), if (length(s) == 1L) ",", ")"))
  }
  if (is.list(x) && !is_sym(x) && !inherits(x, "symmatrix"))
    return(paste0("[", paste(vapply(x, .pv_objrepr, ""), collapse = ", "), "]"))
  .pv_str(x)
}

.pvM_physics_obj___repr__ <- function(self) .pv_objrepr(self)
.pvM_physics_obj___str__ <- function(self) format(self)

# Named tuples: a namedtuple of Python is a tuple, so these errors follow
# the operators SymPy forbids on Inertia, Force and Torque.
.pm_namedtuple <- function(cls, fields) structure(fields, class = c(cls, "pynamedtuple", "physics_obj"))

#' @export
Ops.pynamedtuple <- function(e1, e2) {
  self <- if (inherits(e1, "pynamedtuple")) e1 else e2
  other <- if (inherits(e1, "pynamedtuple")) e2 else e1
  if (!(.Generic %in% c("+", "*", "==", "!="))) .pv_err("TypeError", sprintf("unsupported operand type(s) for %s", .Generic))
  if (.Generic %in% c("==", "!=")) {
    eq <- .pm_tuple_eq(e1, e2)
    return(if (.Generic == "==") eq else !eq)
  }
  .pv_err("TypeError", sprintf("unsupported operand type(s) for %s: '%s' and '%s'", .Generic,
                               class(self)[1L], .pv_tname(other)))
}

.pm_py_eq <- function(a, b) {
  if (inherits(a, c("PhysVector", "PhysDyadic")) || inherits(b, c("PhysVector", "PhysDyadic"))) return(isTRUE(a == b))
  if (is.environment(a) || is.environment(b)) return(identical(a, b))
  if (is_sym(a) || is_sym(b)) return(is_sym(a) && is_sym(b) && identical(a$key, b$key))
  identical(a, b)
}

.pm_tuple_eq <- function(a, b) {
  if (!is.list(a) || !is.list(b) || length(a) != length(b)) return(FALSE)
  all(vapply(seq_along(a), function(i) .pm_py_eq(a[[i]], b[[i]]), TRUE))
}

#' @export
format.pynamedtuple <- function(x, ...) .pv_objrepr(x)

# ---------------------------------------------------------------------------
# inertia.py

#' Inertia dyadics and the Inertia tuple
#'
#' `inertia()` builds the inertia dyadic with the given moments and products
#' of inertia in `frame`, as SymPy's `inertia()`. `inertia_of_point_mass()`
#' is the inertia dyadic of a point mass at `pos_vec`. `Inertia()` pairs an
#' inertia dyadic with its reference point, as SymPy's `Inertia` named
#' tuple; `Inertia_from_inertia_scalars()` is SymPy's
#' `Inertia.from_inertia_scalars()`.
#'
#' @param frame A [ReferenceFrame()].
#' @param ixx,iyy,izz Moments of inertia.
#' @param ixy,iyz,izx Products of inertia.
#' @param mass Mass of the point.
#' @param pos_vec Position vector of the point from the reference point.
#' @param dyadic An inertia dyadic.
#' @param point The reference point, a [PhysPoint()].
#' @return A `PhysDyadic`, or an `Inertia` tuple.
#' @examples
#' N <- ReferenceFrame("N")
#' inertia(N, 1, 2, 3)
#' O <- PhysPoint("O")
#' Inertia(inertia(N, 1, 2, 3), O)
#' @export
inertia <- function(frame, ixx, iyy, izz, ixy = 0, iyz = 0, izx = 0) {
  if (!.pv_is_frame(frame)) .pv_err("TypeError", "Need to define the inertia in a frame")
  ixx <- as_sym(ixx); iyy <- as_sym(iyy); izz <- as_sym(izz)
  ixy <- as_sym(ixy); iyz <- as_sym(iyz); izx <- as_sym(izx)
  x <- frame$x; y <- frame$y; z <- frame$z
  ixx * .pv_outer(x, x) + ixy * .pv_outer(x, y) +
    izx * .pv_outer(x, z) + ixy * .pv_outer(y, x) +
    iyy * .pv_outer(y, y) + iyz * .pv_outer(y, z) +
    izx * .pv_outer(z, x) + iyz * .pv_outer(z, y) +
    izz * .pv_outer(z, z)
}

#' @rdname inertia
#' @export
inertia_of_point_mass <- function(mass, pos_vec, frame) {
  x <- frame$x; y <- frame$y; z <- frame$z
  as_sym(mass) * ((.pv_outer(x, x) + .pv_outer(y, y) + .pv_outer(z, z)) *
                    (.pv_vdot(pos_vec, pos_vec)) - .pv_outer(pos_vec, pos_vec))
}

#' @rdname inertia
#' @export
Inertia <- function(dyadic, point) {
  if (.pv_is_point(dyadic) && .pv_is_dyad(point)) {
    tmp <- point
    point <- dyadic
    dyadic <- tmp
  }
  if (!.pv_is_point(point)) .pv_err("TypeError", "Reference point should be of type Point")
  if (!.pv_is_dyad(dyadic)) .pv_err("TypeError", "Inertia value should be expressed as a Dyadic")
  .pm_namedtuple("Inertia", list(dyadic = dyadic, point = point))
}

#' @rdname inertia
#' @export
Inertia_from_inertia_scalars <- function(point, frame, ixx, iyy, izz, ixy = 0, iyz = 0, izx = 0)
  Inertia(inertia(frame, ixx, iyy, izz, ixy, iyz, izx), point)

.pvR_Inertia <- function(x) sprintf("Inertia(dyadic=%s, point=%s)", .pv_objrepr(x$dyadic), .pv_objrepr(x$point))
.pvM_Inertia_from_inertia_scalars <- function(self, ...) Inertia_from_inertia_scalars(...)

# ---------------------------------------------------------------------------
# body_base.py, particle.py and rigidbody.py

.pm_body_init <- function(self, name, masscenter, mass) {
  if (!(is.character(name) && length(name) == 1L)) .pv_err("TypeError", "Supply a valid name.")
  assign("_uid", .pv_next_uid(), envir = self)
  self[["_name"]] <- name
  if (is.null(mass)) mass <- Symbol(paste0(name, "_mass"))
  if (is.null(masscenter)) masscenter <- PhysPoint(paste0(name, "_masscenter"))
  .pvS_BodyBase_mass(self, mass)
  .pvS_BodyBase_masscenter(self, masscenter)
  .pvS_BodyBase_potential_energy(self, 0)
  self[["points"]] <- list()
  invisible(self)
}

#' @export
format.BodyBase <- function(x, ...) x[["_name"]]

.pvR_BodyBase <- function(x) sprintf("%s(%s, masscenter=%s, mass=%s)", class(x)[1L], .pv_objrepr(x[["_name"]]),
                                     .pv_objrepr(x[["_masscenter"]]), .pv_objrepr(x[["_mass"]]))

.pvP_BodyBase_name <- function(self) self[["_name"]]
.pvP_BodyBase_masscenter <- function(self) self[["_masscenter"]]
.pvS_BodyBase_masscenter <- function(self, point) {
  if (!.pv_is_point(point)) .pv_err("TypeError", "The body's center of mass must be a Point object.")
  self[["_masscenter"]] <- point
}
.pvP_BodyBase_mass <- function(self) self[["_mass"]]
.pvS_BodyBase_mass <- function(self, mass) self[["_mass"]] <- as_sym(mass)
.pvP_BodyBase_potential_energy <- function(self) self[["_potential_energy"]]
.pvS_BodyBase_potential_energy <- function(self, scalar) self[["_potential_energy"]] <- as_sym(scalar)
.pvM_BodyBase_set_potential_energy <- function(self, scalar) {
  warning("set_potential_energy() is deprecated; assign the potential_energy attribute instead.", call. = FALSE)
  .pvS_BodyBase_potential_energy(self, scalar)
  invisible(NULL)
}

#' Particles and rigid bodies
#'
#' `Particle()` and `RigidBody()` create the bodies of SymPy's
#' `sympy.physics.mechanics`. A particle has a mass and a point; a rigid
#' body adds a reference frame and an inertia. Attributes and methods
#' (`mass`, `masscenter`, `frame`, `inertia`, `central_inertia`,
#' `potential_energy`, `linear_momentum()`, `angular_momentum()`,
#' `kinetic_energy()` and `parallel_axis()`) are reached with `$`.
#'
#' @param name Body name.
#' @param point,masscenter The mass center, a [PhysPoint()]; created when
#'   `NULL`.
#' @param mass The mass; a symbol `<name>_mass` when `NULL`.
#' @param frame The body frame; a new [ReferenceFrame()] when `NULL`.
#' @param inertia An [Inertia()] tuple or a list of a dyadic and a point;
#'   built from symbolic moments of inertia when `NULL`.
#' @return A `Particle` or `RigidBody`.
#' @examples
#' N <- ReferenceFrame("N")
#' P <- PhysPoint("P")
#' P$set_vel(N, 3 * N$x)
#' p <- Particle("p", P, 2)
#' p$linear_momentum(N)
#' p$kinetic_energy(N)
#' B <- RigidBody("B")
#' B$central_inertia
#' @export
Particle <- function(name, point = NULL, mass = NULL) {
  self <- new.env(parent = emptyenv())
  class(self) <- c("Particle", "BodyBase", "physics_obj")
  .pm_body_init(self, name, point, mass)
  self
}

.pvP_Particle_point <- function(self) self[["_masscenter"]]
.pvS_Particle_point <- function(self, point) .pvS_BodyBase_masscenter(self, point)

.pvM_Particle_linear_momentum <- function(self, frame) self[["_mass"]] * .pv_point_vel(self[["_masscenter"]], frame)

.pvM_Particle_angular_momentum <- function(self, point, frame)
  .pv_cross_dispatch(.pv_pos_from(self[["_masscenter"]], point),
                     self[["_mass"]] * .pv_point_vel(self[["_masscenter"]], frame))

.pvM_Particle_kinetic_energy <- function(self, frame) {
  v <- .pv_point_vel(self[["_masscenter"]], frame)
  S_Half * self[["_mass"]] * .pv_dot_dispatch(v, v)
}

.pvM_Particle_parallel_axis <- function(self, point, frame)
  inertia_of_point_mass(self[["_mass"]], .pv_pos_from(self[["_masscenter"]], point), frame)

#' @rdname Particle
#' @export
RigidBody <- function(name, masscenter = NULL, frame = NULL, mass = NULL, inertia = NULL) {
  self <- new.env(parent = emptyenv())
  class(self) <- c("RigidBody", "BodyBase", "physics_obj")
  .pm_body_init(self, name, masscenter, mass)
  if (is.null(frame)) frame <- ReferenceFrame(paste0(name, "_frame"))
  .pvS_RigidBody_frame(self, frame)
  if (is.null(inertia)) {
    s <- function(n) Symbol(paste0(name, "_", n))
    inertia <- Inertia_from_inertia_scalars(self[["_masscenter"]], self[["_frame"]],
                                            s("ixx"), s("iyy"), s("izz"), s("ixy"), s("iyz"), s("izx"))
  }
  .pvS_RigidBody_inertia(self, inertia)
  self
}

.pvR_RigidBody <- function(x) sprintf("%s(%s, masscenter=%s, frame=%s, mass=%s, inertia=%s)", class(x)[1L],
                                      .pv_objrepr(x[["_name"]]), .pv_objrepr(x[["_masscenter"]]),
                                      .pv_objrepr(x[["_frame"]]), .pv_objrepr(x[["_mass"]]),
                                      .pv_objrepr(x[["_inertia"]]))

.pvP_RigidBody_frame <- function(self) self[["_frame"]]
.pvS_RigidBody_frame <- function(self, F) {
  if (!.pv_is_frame(F)) .pv_err("TypeError", "RigidBody frame must be a ReferenceFrame object.")
  self[["_frame"]] <- F
}
.pvP_RigidBody_x <- function(self) self[["_frame"]]$x
.pvP_RigidBody_y <- function(self) self[["_frame"]]$y
.pvP_RigidBody_z <- function(self) self[["_frame"]]$z
.pvP_RigidBody_inertia <- function(self) self[["_inertia"]]
.pvS_RigidBody_inertia <- function(self, I) {
  ok <- is.list(I) && !inherits(I, c("PhysVector", "PhysDyadic")) && length(I) == 2L &&
    .pv_is_dyad(I[[1L]]) && .pv_is_point(I[[2L]])
  if (!ok) .pv_err("TypeError", "RigidBody inertia must be a tuple of the form (Dyadic, Point).")
  self[["_inertia"]] <- Inertia(I[[1L]], I[[2L]])
  I_Ss_O <- inertia_of_point_mass(self[["_mass"]], .pv_pos_from(self[["_masscenter"]], I[[2L]]), self[["_frame"]])
  self[["_central_inertia"]] <- I[[1L]] - I_Ss_O
}
.pvP_RigidBody_central_inertia <- function(self) self[["_central_inertia"]]
.pvS_RigidBody_central_inertia <- function(self, I) {
  if (!.pv_is_dyad(I)) .pv_err("TypeError", "RigidBody inertia must be a Dyadic object.")
  .pvS_RigidBody_inertia(self, Inertia(I, self[["_masscenter"]]))
}

.pvM_RigidBody_linear_momentum <- function(self, frame) self[["_mass"]] * .pv_point_vel(self[["_masscenter"]], frame)

.pvM_RigidBody_angular_momentum <- function(self, point, frame) {
  I <- self[["_central_inertia"]]
  w <- .pv_ang_vel_in(self[["_frame"]], frame)
  m <- self[["_mass"]]
  r <- .pv_pos_from(self[["_masscenter"]], point)
  v <- .pv_point_vel(self[["_masscenter"]], frame)
  .pv_ddot(I, w) + .pv_vcross(r, m * v)
}

.pvM_RigidBody_kinetic_energy <- function(self, frame) {
  w <- .pv_ang_vel_in(self[["_frame"]], frame)
  rotational_KE <- S_Half * .pv_dot_dispatch(w, .pv_dot_dispatch(self[["_central_inertia"]], w))
  v <- .pv_point_vel(self[["_masscenter"]], frame)
  translational_KE <- S_Half * self[["_mass"]] * .pv_dot_dispatch(v, v)
  rotational_KE + translational_KE
}

.pvM_RigidBody_parallel_axis <- function(self, point, frame = NULL) {
  if (is.null(frame)) frame <- self[["_frame"]]
  self[["_central_inertia"]] + inertia_of_point_mass(self[["_mass"]], .pv_pos_from(self[["_masscenter"]], point), frame)
}

# ---------------------------------------------------------------------------
# loads.py

#' Loads: forces and torques
#'
#' `Force()` is a force applied at a point and `Torque()` a torque applied
#' to a frame, as SymPy's named tuples of the same names. A body may stand
#' for its mass center (forces) or its frame (torques). `gravity()` returns
#' the gravitational forces on the given bodies.
#'
#' @param point A [PhysPoint()] or a body.
#' @param force,torque A vector.
#' @param frame A [ReferenceFrame()] or a rigid body.
#' @param acceleration The gravitational acceleration vector.
#' @param ... Bodies.
#' @return A `Force` or `Torque`; a list of forces for `gravity()`.
#' @examples
#' N <- ReferenceFrame("N")
#' P <- PhysPoint("P")
#' Force(P, 3 * N$x)
#' Torque(N, 2 * N$z)
#' gravity(-9.81 * N$y, Particle("p", P, 1))
#' @export
Force <- function(point, force) {
  if (.pm_is_body(point)) point <- point[["_masscenter"]]
  if (!.pv_is_point(point)) .pv_err("TypeError", "Force location should be a Point.")
  if (!.pv_is_vec(force)) .pv_err("TypeError", "Force vector should be a Vector.")
  .pm_namedtuple(c("Force", "LoadBase"), list(location = point, vector = force))
}

.pvR_Force <- function(x) sprintf("Force(point=%s, force=%s)", .pv_str(x$location), .pv_str(x$vector))
.pvP_Force_point <- function(self) .subset2(self, "location")
.pvP_Force_force <- function(self) .subset2(self, "vector")

#' @rdname Force
#' @export
Torque <- function(frame, torque) {
  if (.pm_is_body(frame)) frame <- frame$frame
  if (!.pv_is_frame(frame)) .pv_err("TypeError", "Torque location should be a ReferenceFrame.")
  if (!.pv_is_vec(torque)) .pv_err("TypeError", "Torque vector should be a Vector.")
  .pm_namedtuple(c("Torque", "LoadBase"), list(location = frame, vector = torque))
}

.pvR_Torque <- function(x) sprintf("Torque(frame=%s, torque=%s)", .pv_str(x$location), .pv_str(x$vector))
.pvP_Torque_frame <- function(self) .subset2(self, "location")
.pvP_Torque_torque <- function(self) .subset2(self, "vector")

#' @rdname Force
#' @export
gravity <- function(acceleration, ...) {
  gravity_force <- list()
  for (body in list(...)) {
    if (!.pm_is_body(body)) .pv_err("TypeError", sprintf("%s is not a body type", .pv_pyclass(body)))
    gravity_force[[length(gravity_force) + 1L]] <- Force(body[["_masscenter"]], body[["_mass"]] * acceleration)
  }
  gravity_force
}

.pm_is_tuple <- function(x) inherits(x, "pytuple")

.pm_parse_load <- function(load) {
  if (inherits(load, "LoadBase")) return(load)
  if (.pm_is_tuple(load)) {
    if (length(load) != 2L) .pv_err("ValueError", sprintf("Load %s should have a length of 2.", .pv_objrepr(load)))
    if (.pv_is_point(load[[1L]])) return(Force(load[[1L]], load[[2L]]))
    if (.pv_is_frame(load[[1L]])) return(Torque(load[[1L]], load[[2L]]))
    .pv_err("ValueError", sprintf("Load not recognized. The load location %s should either be a Point or a ReferenceFrame.",
                                  .pv_str(load[[1L]])))
  }
  .pv_err("TypeError", sprintf("Load type %s not recognized as a load. It should be a Force, Torque or tuple.",
                               .pv_pyclass(load)))
}
