# Joints and the deprecated Body and JointsMethod classes, after
# sympy/physics/mechanics/joint.py, body.py and jointsmethod.py.
#
# A joint is a mutable environment of class c(<JointType>, "Joint",
# "physics_obj"). The joint types supply the hooks `_generate_coordinates`,
# `_generate_speeds`, `_orient_frames`, `_set_angular_velocity` and
# `_set_linear_velocity` as `.pvM_<JointType>_<hook>` methods, which the
# shared constructor calls in SymPy's order.

.pm_is_joint <- function(x) inherits(x, "Joint")

.pm_hook <- function(self, name, ...) .pv_lookup(self, ".pvM_", name)(self, ...)

# sympy.utilities.iterables.iterable() for joint arguments
.pm_iter_list <- function(x) {
  if (is.null(x)) return(list())
  if (inherits(x, "symmatrix")) return(x$data)
  if (is.list(x) && !is_sym(x) && !inherits(x, "physics_obj")) return(unclass(x))
  NULL
}

.pm_is_dyn_coord <- function(v) is_sym(v) && (v$type == "Derivative" || .is_undefined_application(v))

#' Joints between bodies
#'
#' The joints of SymPy's `sympy.physics.mechanics`: `PinJoint()` (one
#' rotation about `joint_axis`), `PrismaticJoint()` (one translation along
#' `joint_axis`), `CylindricalJoint()` (a rotation and a translation along
#' the same axis), `PlanarJoint()` (a rotation about the parent
#' interframe's x axis and two translations in its y z plane),
#' `SphericalJoint()` (three rotations) and `WeldJoint()` (no motion).
#' Each joint orients the child frame and sets the velocities of the
#' attachment points when it is created. `Joint()` is the abstract base and
#' cannot be instantiated. Attributes (`coordinates`, `speeds`, `kdes`,
#' `parent_point`, `child_point`, `parent_interframe`,
#' `child_interframe`, `joint_axis` and others) are reached with `$`.
#'
#' @param name Joint name.
#' @param parent,child The parent and child bodies.
#' @param coordinates,speeds Generalized coordinates and speeds; created
#'   with [dynamicsymbols()] when `NULL`.
#' @param parent_point,child_point Attachment points, as points or
#'   position vectors from the mass centers.
#' @param parent_interframe,child_interframe Intermediate frames, as frames
#'   or vectors to align with the x axis.
#' @param parent_axis,child_axis,parent_joint_pos,child_joint_pos
#'   Deprecated arguments of SymPy.
#' @param joint_axis The axis of rotation or translation.
#' @param rotation_coordinate,translation_coordinate,rotation_speed,translation_speed
#'   Coordinates and speeds of a cylindrical or planar joint.
#' @param planar_coordinates,planar_speeds The two planar coordinates and
#'   speeds of a planar joint.
#' @param rot_type,amounts,rot_order Orientation of a spherical joint.
#' @return A joint object.
#' @examples
#' P <- RigidBody("P")
#' C <- RigidBody("C")
#' J <- PinJoint("J", P, C)
#' J$coordinates
#' J$kdes
#' C$frame$dcm(P$frame)
#' @export
Joint <- function(name, parent, child, coordinates = NULL, speeds = NULL, parent_point = NULL,
                  child_point = NULL, parent_interframe = NULL, child_interframe = NULL,
                  parent_axis = NULL, child_axis = NULL, parent_joint_pos = NULL,
                  child_joint_pos = NULL) {
  .pv_err("TypeError", paste0("Can't instantiate abstract class Joint without an implementation for abstract ",
                              "methods '_generate_coordinates', '_generate_speeds', '_orient_frames', ",
                              "'_set_angular_velocity', '_set_linear_velocity'"))
}

.pm_joint_new <- function(cls, name, parent, child, coordinates, speeds, parent_point, child_point,
                          parent_interframe, child_interframe, parent_axis = NULL, child_axis = NULL,
                          parent_joint_pos = NULL, child_joint_pos = NULL, fields = list()) {
  self <- new.env(parent = emptyenv())
  class(self) <- c(cls, "Joint", "physics_obj")
  assign("_uid", .pv_next_uid(), envir = self)
  for (nm in names(fields)) assign(nm, fields[[nm]], envir = self)
  if (!(is.character(name) && length(name) == 1L)) .pv_err("TypeError", "Supply a valid name.")
  self[["_name"]] <- name
  if (!.pm_is_body(parent)) .pv_err("TypeError", "Parent must be a body.")
  self[["_parent"]] <- parent
  if (!.pm_is_body(child)) .pv_err("TypeError", "Child must be a body.")
  self[["_child"]] <- child
  if (!is.null(parent_axis) || !is.null(child_axis)) {
    warning("The parent_axis and child_axis arguments of Joint are deprecated; use parent_interframe and child_interframe.",
            call. = FALSE)
    if (is.null(parent_interframe)) parent_interframe <- parent_axis
    if (is.null(child_interframe)) child_interframe <- child_axis
  }
  self[["_parent_frame"]] <- if (inherits(parent, "RigidBody")) parent$frame
                             else if (.pv_is_frame(parent_interframe)) parent_interframe
                             else ReferenceFrame(sprintf("%s_%s_frame", name, parent$name))
  self[["_child_frame"]] <- if (inherits(child, "RigidBody")) child$frame
                            else if (.pv_is_frame(child_interframe)) child_interframe
                            else ReferenceFrame(sprintf("%s_%s_frame", name, child$name))
  self[["_parent_interframe"]] <- .pm_locate_joint_frame(self, parent, parent_interframe, self[["_parent_frame"]])
  self[["_child_interframe"]] <- .pm_locate_joint_frame(self, child, child_interframe, self[["_child_frame"]])
  self[["_parent_axis"]] <- .pm_joint_axis(parent_axis, self[["_parent_frame"]])
  self[["_child_axis"]] <- .pm_joint_axis(child_axis, self[["_child_frame"]])
  if (!is.null(parent_joint_pos) || !is.null(child_joint_pos)) {
    warning("The parent_joint_pos and child_joint_pos arguments of Joint are deprecated; use parent_point and child_point.",
            call. = FALSE)
    if (is.null(parent_point)) parent_point <- parent_joint_pos
    if (is.null(child_point)) child_point <- child_joint_pos
  }
  self[["_parent_point"]] <- .pm_locate_joint_pos(self, parent, parent_point, self[["_parent_frame"]])
  self[["_child_point"]] <- .pm_locate_joint_pos(self, child, child_point, self[["_child_frame"]])
  self[["_coordinates"]] <- .pm_hook(self, "_generate_coordinates", coordinates)
  self[["_speeds"]] <- .pm_hook(self, "_generate_speeds", speeds)
  .pm_validate_coordinates(self[["_coordinates"]], self[["_speeds"]])
  self[["_kdes"]] <- .pvM_Joint__generate_kdes(self)
  .pm_hook(self, "_orient_frames")
  .pm_hook(self, "_set_angular_velocity")
  .pm_hook(self, "_set_linear_velocity")
  self
}

#' @export
format.Joint <- function(x, ...) {
  cl <- class(x)[1L]
  if (cl %in% c("PinJoint", "PrismaticJoint", "CylindricalJoint", "PlanarJoint", "SphericalJoint", "WeldJoint"))
    return(sprintf("%s: %s  parent: %s  child: %s", cl, x[["_name"]], format(x[["_parent"]]), format(x[["_child"]])))
  x[["_name"]]
}
.pvR_Joint <- function(x) format(x)

.pvP_Joint_name <- function(self) self[["_name"]]
.pvP_Joint_parent <- function(self) self[["_parent"]]
.pvP_Joint_child <- function(self) self[["_child"]]
.pvP_Joint_coordinates <- function(self) self[["_coordinates"]]
.pvP_Joint_speeds <- function(self) self[["_speeds"]]
.pvP_Joint_kdes <- function(self) self[["_kdes"]]
.pvP_Joint_parent_axis <- function(self) self[["_parent_axis"]]
.pvP_Joint_child_axis <- function(self) self[["_child_axis"]]
.pvP_Joint_parent_point <- function(self) self[["_parent_point"]]
.pvP_Joint_child_point <- function(self) self[["_child_point"]]
.pvP_Joint_parent_interframe <- function(self) self[["_parent_interframe"]]
.pvP_Joint_child_interframe <- function(self) self[["_child_interframe"]]

.pvM_Joint__to_vector <- function(self, matrix, frame) PhysVector(list(list(matrix, frame)))

# Joint._axis(ax, *frames)
.pm_joint_axis <- function(ax, ...) {
  frames <- list(...)
  if (is.null(ax)) return(frames[[1L]]$x)
  if (!.pv_is_vec(ax)) .pv_err("TypeError", "Axis must be a Vector.")
  ref_frame <- NULL
  for (frame in frames) {
    ok <- tryCatch({ .pv_to_matrix(ax, frame); TRUE }, ValueError = function(e) FALSE)
    if (ok) {
      ref_frame <- frame
      break
    }
  }
  if (is.null(ref_frame)) .pv_err("ValueError", "Axis cannot be expressed in one of the body's frames.")
  if (!isTRUE(.pv_vec_eq(.pv_time_derivative(ax, ref_frame), 0)))
    .pv_err("ValueError", "Axis cannot be time-varying when viewed from the associated body.")
  ax
}
.pvM_Joint__axis <- function(self, ax, ...) .pm_joint_axis(ax, ...)

.pm_ne0 <- function(e) !identical(as_sym(e)$key, S_Zero$key)

.pm_choose_rotation_axis <- function(frame, axis) {
  components <- .pv_to_matrix(axis, frame)$data
  x <- components[[1L]]; y <- components[[2L]]; z <- components[[3L]]
  if (.pm_ne0(x)) {
    if (.pm_ne0(y) && .pm_ne0(z)) return(.pv_vcross(axis, frame$x))
    if (.pm_ne0(z)) return(frame$y)
    return(frame$z)
  }
  if (.pm_ne0(y)) return(frame$x)
  frame$y
}
.pvM_Joint__choose_rotation_axis <- function(self, frame, axis) .pm_choose_rotation_axis(frame, axis)

.pm_create_aligned_interframe <- function(frame, align_axis, frame_axis = NULL, frame_name = NULL) {
  if (.pm_is_body(frame)) frame <- frame$frame
  if (is.null(frame_axis)) frame_axis <- frame$x
  if (is.null(frame_name)) {
    nm <- frame$name
    frame_name <- if (endsWith(nm, "_frame")) paste0(substr(nm, 1L, nchar(nm) - 6L), "_int_frame")
                  else paste0(nm, "_int_frame")
  }
  angle <- .pv_angle_between(frame_axis, align_axis)
  rotation_axis <- .pv_vcross(frame_axis, align_axis)
  if (isTRUE(.pv_vec_eq(rotation_axis, 0)) && !.pm_ne0(angle)) return(frame)
  if (identical(angle$key, S_Pi$key)) rotation_axis <- .pm_choose_rotation_axis(frame, align_axis)
  int_frame <- ReferenceFrame(frame_name)
  .pv_orient_axis(int_frame, frame, rotation_axis, angle)
  int_frame$set_ang_vel(frame, 0 * rotation_axis)
  int_frame
}
.pvM_Joint__create_aligned_interframe <- function(self, frame, align_axis, frame_axis = NULL, frame_name = NULL)
  .pm_create_aligned_interframe(frame, align_axis, frame_axis, frame_name)

.pvM_Joint__generate_kdes <- function(self) {
  t <- .pv_t()
  q <- self[["_coordinates"]]$data
  u <- self[["_speeds"]]$data
  .pm_colmat(lapply(seq_along(q), function(i) .add(list(.neg(.pv_diff(q[[i]], t)), u[[i]]))))
}

.pm_locate_joint_pos <- function(self, body, joint_pos, body_frame = NULL) {
  if (is.null(body_frame)) body_frame <- body$frame
  if (is.null(joint_pos)) return(body$masscenter)
  if (!(.pv_is_point(joint_pos) || .pv_is_vec(joint_pos)))
    .pv_err("TypeError", "Attachment point must be a Point or Vector.")
  if (.pv_is_vec(joint_pos))
    joint_pos <- body$masscenter$locatenew(sprintf("%s_%s_joint", self[["_name"]], body$name), joint_pos)
  if (!isTRUE(.pv_vec_eq(.pv_time_derivative(.pv_pos_from(joint_pos, body$masscenter), body_frame), 0)))
    .pv_err("ValueError", "Attachment point must be fixed to the associated body.")
  joint_pos
}
.pvM_Joint__locate_joint_pos <- function(self, body, joint_pos, body_frame = NULL)
  .pm_locate_joint_pos(self, body, joint_pos, body_frame)

.pm_locate_joint_frame <- function(self, body, interframe, body_frame = NULL) {
  if (is.null(body_frame)) body_frame <- body$frame
  if (is.null(interframe)) return(body_frame)
  if (.pv_is_vec(interframe)) {
    interframe <- .pm_create_aligned_interframe(body_frame, interframe,
                                                frame_name = sprintf("%s_%s_int_frame", self[["_name"]], body$name))
  } else if (!.pv_is_frame(interframe)) .pv_err("TypeError", "Interframe must be a ReferenceFrame.")
  if (!isTRUE(.pv_vec_eq(.pv_ang_vel_in(interframe, body_frame), 0)))
    .pv_err("ValueError", sprintf("Interframe %s is not fixed to body %s.", format(interframe), format(body)))
  body$masscenter$set_vel(interframe, 0)
  interframe
}
.pvM_Joint__locate_joint_frame <- function(self, body, interframe, body_frame = NULL)
  .pm_locate_joint_frame(self, body, interframe, body_frame)

.pvM_Joint__fill_coordinate_list <- function(self, coordinates, n_coords, label = "q", offset = 0,
                                             number_single = FALSE) {
  create_symbol <- function(number) {
    if (n_coords == 1 && !isTRUE(number_single)) return(dynamicsymbols(sprintf("%s_%s", label, self[["_name"]])))
    dynamicsymbols(sprintf("%s%d_%s", label, as.integer(number), self[["_name"]]))
  }
  name <- if (identical(label, "q")) "generalized coordinate" else "generalized speed"
  generated <- list()
  coords <- .pm_iter_list(coordinates)
  if (is.null(coords)) coords <- list(coordinates)
  if (!(length(coords) == 0L || length(coords) == n_coords))
    .pv_err("ValueError", sprintf("Expected %d %ss, instead got %d %ss.", as.integer(n_coords), name, length(coords), name))
  for (i in seq_along(coords)) {
    coord <- coords[[i]]
    if (is.null(coord)) generated[[length(generated) + 1L]] <- create_symbol(i - 1L + offset)
    else if (.pm_is_dyn_coord(coord)) generated[[length(generated) + 1L]] <- coord
    else .pv_err("TypeError", sprintf("The %s %s should have been a dynamicsymbol.", name, .pv_str(coord)))
  }
  rng <- seq_len(max(0L, n_coords - length(coords))) - 1L + length(coords) + offset
  for (i in rng) generated[[length(generated) + 1L]] <- create_symbol(i)
  .pm_colmat(generated)
}

.pm_qi <- function(self, i) self[["_coordinates"]]$data[[i]]
.pm_ui <- function(self, i) self[["_speeds"]]$data[[i]]

# ---------------------------------------------------------------------------
# PinJoint

#' @rdname Joint
#' @export
PinJoint <- function(name, parent, child, coordinates = NULL, speeds = NULL, parent_point = NULL,
                     child_point = NULL, parent_interframe = NULL, child_interframe = NULL,
                     parent_axis = NULL, child_axis = NULL, joint_axis = NULL,
                     parent_joint_pos = NULL, child_joint_pos = NULL)
  .pm_joint_new("PinJoint", name, parent, child, coordinates, speeds, parent_point, child_point,
                parent_interframe, child_interframe, parent_axis, child_axis, parent_joint_pos,
                child_joint_pos, fields = list(`_joint_axis` = joint_axis))

.pvP_PinJoint_joint_axis <- function(self) self[["_joint_axis"]]
.pvM_PinJoint__generate_coordinates <- function(self, coordinate) .pvM_Joint__fill_coordinate_list(self, coordinate, 1, "q")
.pvM_PinJoint__generate_speeds <- function(self, speed) .pvM_Joint__fill_coordinate_list(self, speed, 1, "u")
.pvM_PinJoint__orient_frames <- function(self) {
  self[["_joint_axis"]] <- .pm_joint_axis(self[["_joint_axis"]], self[["_parent_interframe"]])
  .pv_orient_axis(self[["_child_interframe"]], self[["_parent_interframe"]], self[["_joint_axis"]], .pm_qi(self, 1L))
}
.pvM_PinJoint__set_angular_velocity <- function(self)
  self[["_child_interframe"]]$set_ang_vel(self[["_parent_interframe"]], .pm_ui(self, 1L) * .pv_normalize(self[["_joint_axis"]]))
.pvM_PinJoint__set_linear_velocity <- function(self) {
  self[["_child_point"]]$set_pos(self[["_parent_point"]], 0)
  self[["_parent_point"]]$set_vel(self[["_parent_frame"]], 0)
  self[["_child_point"]]$set_vel(self[["_child_frame"]], 0)
  self[["_child"]]$masscenter$v2pt_theory(self[["_parent_point"]], self[["_parent_frame"]], self[["_child_frame"]])
}

# ---------------------------------------------------------------------------
# PrismaticJoint

#' @rdname Joint
#' @export
PrismaticJoint <- function(name, parent, child, coordinates = NULL, speeds = NULL, parent_point = NULL,
                           child_point = NULL, parent_interframe = NULL, child_interframe = NULL,
                           parent_axis = NULL, child_axis = NULL, joint_axis = NULL,
                           parent_joint_pos = NULL, child_joint_pos = NULL)
  .pm_joint_new("PrismaticJoint", name, parent, child, coordinates, speeds, parent_point, child_point,
                parent_interframe, child_interframe, parent_axis, child_axis, parent_joint_pos,
                child_joint_pos, fields = list(`_joint_axis` = joint_axis))

.pvP_PrismaticJoint_joint_axis <- function(self) self[["_joint_axis"]]
.pvM_PrismaticJoint__generate_coordinates <- function(self, coordinate) .pvM_Joint__fill_coordinate_list(self, coordinate, 1, "q")
.pvM_PrismaticJoint__generate_speeds <- function(self, speed) .pvM_Joint__fill_coordinate_list(self, speed, 1, "u")
.pvM_PrismaticJoint__orient_frames <- function(self) {
  self[["_joint_axis"]] <- .pm_joint_axis(self[["_joint_axis"]], self[["_parent_interframe"]])
  .pv_orient_axis(self[["_child_interframe"]], self[["_parent_interframe"]], self[["_joint_axis"]], S_Zero)
}
.pvM_PrismaticJoint__set_angular_velocity <- function(self)
  self[["_child_interframe"]]$set_ang_vel(self[["_parent_interframe"]], 0)
.pvM_PrismaticJoint__set_linear_velocity <- function(self) {
  axis <- .pv_normalize(self[["_joint_axis"]])
  self[["_child_point"]]$set_pos(self[["_parent_point"]], .pm_qi(self, 1L) * axis)
  self[["_parent_point"]]$set_vel(self[["_parent_frame"]], 0)
  self[["_child_point"]]$set_vel(self[["_child_frame"]], 0)
  self[["_child_point"]]$set_vel(self[["_parent_frame"]], .pm_ui(self, 1L) * axis)
  self[["_child"]]$masscenter$set_vel(self[["_parent_frame"]], .pm_ui(self, 1L) * axis)
}

# ---------------------------------------------------------------------------
# CylindricalJoint

#' @rdname Joint
#' @export
CylindricalJoint <- function(name, parent, child, rotation_coordinate = NULL, translation_coordinate = NULL,
                             rotation_speed = NULL, translation_speed = NULL, parent_point = NULL,
                             child_point = NULL, parent_interframe = NULL, child_interframe = NULL,
                             joint_axis = NULL)
  .pm_joint_new("CylindricalJoint", name, parent, child, list(rotation_coordinate, translation_coordinate),
                list(rotation_speed, translation_speed), parent_point, child_point,
                parent_interframe, child_interframe, fields = list(`_joint_axis` = joint_axis))

.pvP_CylindricalJoint_joint_axis <- function(self) self[["_joint_axis"]]
.pvP_CylindricalJoint_rotation_coordinate <- function(self) .pm_qi(self, 1L)
.pvP_CylindricalJoint_translation_coordinate <- function(self) .pm_qi(self, 2L)
.pvP_CylindricalJoint_rotation_speed <- function(self) .pm_ui(self, 1L)
.pvP_CylindricalJoint_translation_speed <- function(self) .pm_ui(self, 2L)
.pvM_CylindricalJoint__generate_coordinates <- function(self, coordinates) .pvM_Joint__fill_coordinate_list(self, coordinates, 2, "q")
.pvM_CylindricalJoint__generate_speeds <- function(self, speeds) .pvM_Joint__fill_coordinate_list(self, speeds, 2, "u")
.pvM_CylindricalJoint__orient_frames <- function(self) {
  self[["_joint_axis"]] <- .pm_joint_axis(self[["_joint_axis"]], self[["_parent_interframe"]])
  .pv_orient_axis(self[["_child_interframe"]], self[["_parent_interframe"]], self[["_joint_axis"]], .pm_qi(self, 1L))
}
.pvM_CylindricalJoint__set_angular_velocity <- function(self)
  self[["_child_interframe"]]$set_ang_vel(self[["_parent_interframe"]], .pm_ui(self, 1L) * .pv_normalize(self[["_joint_axis"]]))
.pvM_CylindricalJoint__set_linear_velocity <- function(self) {
  self[["_child_point"]]$set_pos(self[["_parent_point"]], .pm_qi(self, 2L) * .pv_normalize(self[["_joint_axis"]]))
  self[["_parent_point"]]$set_vel(self[["_parent_frame"]], 0)
  self[["_child_point"]]$set_vel(self[["_child_frame"]], 0)
  self[["_child_point"]]$set_vel(self[["_parent_frame"]], .pm_ui(self, 2L) * .pv_normalize(self[["_joint_axis"]]))
  self[["_child"]]$masscenter$v2pt_theory(self[["_child_point"]], self[["_parent_frame"]], self[["_child_interframe"]])
}

# ---------------------------------------------------------------------------
# PlanarJoint

#' @rdname Joint
#' @export
PlanarJoint <- function(name, parent, child, rotation_coordinate = NULL, planar_coordinates = NULL,
                        rotation_speed = NULL, planar_speeds = NULL, parent_point = NULL, child_point = NULL,
                        parent_interframe = NULL, child_interframe = NULL)
  .pm_joint_new("PlanarJoint", name, parent, child, list(rotation_coordinate, planar_coordinates),
                list(rotation_speed, planar_speeds), parent_point, child_point,
                parent_interframe, child_interframe)

.pvP_PlanarJoint_rotation_coordinate <- function(self) .pm_qi(self, 1L)
.pvP_PlanarJoint_planar_coordinates <- function(self) .pm_rows(self[["_coordinates"]], 1L, self[["_coordinates"]]$nrow)
.pvP_PlanarJoint_rotation_speed <- function(self) .pm_ui(self, 1L)
.pvP_PlanarJoint_planar_speeds <- function(self) .pm_rows(self[["_speeds"]], 1L, self[["_speeds"]]$nrow)
.pvP_PlanarJoint_rotation_axis <- function(self) self[["_parent_interframe"]]$x
.pvP_PlanarJoint_planar_vectors <- function(self) list(self[["_parent_interframe"]]$y, self[["_parent_interframe"]]$z)
.pvM_PlanarJoint__generate_coordinates <- function(self, coordinates) {
  rotation <- .pvM_Joint__fill_coordinate_list(self, coordinates[[1L]], 1, "q", number_single = TRUE)
  planar <- .pvM_Joint__fill_coordinate_list(self, coordinates[[2L]], 2, "q", 1)
  .pm_cjoin(rotation, planar)
}
.pvM_PlanarJoint__generate_speeds <- function(self, speeds) {
  rotation <- .pvM_Joint__fill_coordinate_list(self, speeds[[1L]], 1, "u", number_single = TRUE)
  planar <- .pvM_Joint__fill_coordinate_list(self, speeds[[2L]], 2, "u", 1)
  .pm_cjoin(rotation, planar)
}
.pvM_PlanarJoint__orient_frames <- function(self)
  .pv_orient_axis(self[["_child_interframe"]], self[["_parent_interframe"]], self$rotation_axis, .pm_qi(self, 1L))
.pvM_PlanarJoint__set_angular_velocity <- function(self)
  self[["_child_interframe"]]$set_ang_vel(self[["_parent_interframe"]], .pm_ui(self, 1L) * self$rotation_axis)
.pvM_PlanarJoint__set_linear_velocity <- function(self) {
  pv <- self$planar_vectors
  self[["_child_point"]]$set_pos(self[["_parent_point"]], .pm_qi(self, 2L) * pv[[1L]] + .pm_qi(self, 3L) * pv[[2L]])
  self[["_parent_point"]]$set_vel(self[["_parent_interframe"]], 0)
  self[["_child_point"]]$set_vel(self[["_child_interframe"]], 0)
  self[["_child_point"]]$set_vel(self[["_parent_frame"]], .pm_ui(self, 2L) * pv[[1L]] + .pm_ui(self, 3L) * pv[[2L]])
  self[["_child"]]$masscenter$v2pt_theory(self[["_child_point"]], self[["_parent_frame"]], self[["_child_frame"]])
}

# ---------------------------------------------------------------------------
# SphericalJoint

#' @rdname Joint
#' @export
SphericalJoint <- function(name, parent, child, coordinates = NULL, speeds = NULL, parent_point = NULL,
                           child_point = NULL, parent_interframe = NULL, child_interframe = NULL,
                           rot_type = "BODY", amounts = NULL, rot_order = 123)
  .pm_joint_new("SphericalJoint", name, parent, child, coordinates, speeds, parent_point, child_point,
                parent_interframe, child_interframe,
                fields = list(`_rot_type` = rot_type, `_amounts` = amounts, `_rot_order` = rot_order))

.pvM_SphericalJoint__generate_coordinates <- function(self, coordinates) .pvM_Joint__fill_coordinate_list(self, coordinates, 3, "q")
.pvM_SphericalJoint__generate_speeds <- function(self, speeds)
  .pvM_Joint__fill_coordinate_list(self, speeds, .pm_len(self[["_coordinates"]]), "u")
.pvM_SphericalJoint__orient_frames <- function(self) {
  rot_type <- self[["_rot_type"]]
  if (!(toupper(rot_type) %in% c("BODY", "SPACE")))
    .pv_err("NotImplementedError", sprintf("Rotation type \"%s\" is not implemented. Implemented rotation types are: ('BODY', 'SPACE')",
                                           rot_type))
  amounts <- if (is.null(self[["_amounts"]])) self[["_coordinates"]] else self[["_amounts"]]
  .pv_orient_dispatch(self[["_child_interframe"]], self[["_parent_interframe"]], rot_type, amounts, self[["_rot_order"]])
}
.pvM_SphericalJoint__set_angular_velocity <- function(self) {
  t <- .pv_t()
  q <- self[["_coordinates"]]$data
  u <- self[["_speeds"]]$data
  n <- min(length(q), length(u))
  rule <- .pydict(lapply(seq_len(n), function(i) .pv_diff(q[[i]], t)), u[seq_len(n)])
  vel <- .pvM_PhysVector_xreplace(.pv_ang_vel_in(self[["_child_interframe"]], self[["_parent_interframe"]]), rule)
  self[["_child_interframe"]]$set_ang_vel(self[["_parent_interframe"]], vel)
}
.pvM_SphericalJoint__set_linear_velocity <- function(self) {
  self[["_child_point"]]$set_pos(self[["_parent_point"]], 0)
  self[["_parent_point"]]$set_vel(self[["_parent_frame"]], 0)
  self[["_child_point"]]$set_vel(self[["_child_frame"]], 0)
  self[["_child"]]$masscenter$v2pt_theory(self[["_parent_point"]], self[["_parent_frame"]], self[["_child_frame"]])
}

# ---------------------------------------------------------------------------
# WeldJoint

#' @rdname Joint
#' @export
WeldJoint <- function(name, parent, child, parent_point = NULL, child_point = NULL,
                      parent_interframe = NULL, child_interframe = NULL) {
  self <- .pm_joint_new("WeldJoint", name, parent, child, list(), list(), parent_point, child_point,
                        parent_interframe, child_interframe)
  self[["_kdes"]] <- .pm_empty(0L, 1L)
  self
}

.pvM_WeldJoint__generate_coordinates <- function(self, coordinate) .pm_empty()
.pvM_WeldJoint__generate_speeds <- function(self, speed) .pm_empty()
.pvM_WeldJoint__orient_frames <- function(self)
  .pv_orient_axis(self[["_child_interframe"]], self[["_parent_interframe"]], self[["_parent_interframe"]]$x, S_Zero)
.pvM_WeldJoint__set_angular_velocity <- function(self)
  self[["_child_interframe"]]$set_ang_vel(self[["_parent_interframe"]], 0)
.pvM_WeldJoint__set_linear_velocity <- function(self) {
  self[["_child_point"]]$set_pos(self[["_parent_point"]], 0)
  self[["_parent_point"]]$set_vel(self[["_parent_frame"]], 0)
  self[["_child_point"]]$set_vel(self[["_child_frame"]], 0)
  self[["_child"]]$masscenter$set_vel(self[["_parent_frame"]], 0)
}

# ---------------------------------------------------------------------------
# body.py: the deprecated Body class, both a rigid body and a particle.

#' Bodies of the deprecated joints framework
#'
#' `Body()` is SymPy's deprecated `Body`: a rigid body, or a particle when
#' only a mass is given, that also stores the loads applied to it. Methods
#' `apply_force()`, `apply_torque()`, `clear_loads()`, `remove_load()`,
#' `masscenter_vel()`, `ang_vel_in()`, `dcm()`, `kinetic_energy()` and
#' `parallel_axis()` and the attributes `loads` and `is_rigidbody` are
#' reached with `$`. Creating a body emits a deprecation warning, as in
#' SymPy; new code should use [RigidBody()] or [Particle()].
#'
#' @param name Body name.
#' @param masscenter The mass center; created when `NULL`.
#' @param mass The mass; a symbol when `NULL`.
#' @param frame The body frame; created when `NULL`.
#' @param central_inertia The central inertia dyadic.
#' @return A `Body`.
#' @examples
#' B <- suppressWarnings(Body("B"))
#' B$apply_force(B$x)
#' B$loads
#' @export
Body <- function(name, masscenter = NULL, mass = NULL, frame = NULL, central_inertia = NULL) {
  warning("Body is deprecated; use RigidBody or Particle instead.", call. = FALSE)
  self <- new.env(parent = emptyenv())
  class(self) <- c("Body", "RigidBody", "Particle", "BodyBase", "physics_obj")
  self[["_loads"]] <- list()
  if (is.null(frame)) frame <- ReferenceFrame(paste0(name, "_frame"))
  if (is.null(masscenter)) masscenter <- PhysPoint(paste0(name, "_masscenter"))
  if (is.null(central_inertia) && is.null(mass)) {
    s <- function(n) Symbol(paste0(name, "_", n))
    inert <- Inertia_from_inertia_scalars(masscenter, frame, s("ixx"), s("iyy"), s("izz"), s("ixy"), s("iyz"), s("izx"))
  } else {
    inert <- .pytuple(list(central_inertia, masscenter))
  }
  particle <- is.null(central_inertia) && !is.null(mass)
  if (is.null(mass)) mass <- Symbol(paste0(name, "_mass"))
  masscenter$set_vel(frame, 0)
  .pm_body_init(self, name, masscenter, mass)
  .pvS_RigidBody_frame(self, frame)
  if (particle) {
    self[["_central_inertia"]] <- .pv_dyad0()
  } else {
    .pvS_RigidBody_inertia(self, inert)
  }
  self
}

.pvR_Body <- function(x) {
  if (.pvP_Body_is_rigidbody(x)) return(.pvR_RigidBody(x))
  .pvR_BodyBase(x)
}

.pvP_Body_loads <- function(self) self[["_loads"]]
.pvP_Body_is_rigidbody <- function(self) exists("_inertia", envir = self, inherits = FALSE)
.pvP_Body_inertia <- function(self) {
  if (.pvP_Body_is_rigidbody(self)) return(self[["_inertia"]])
  .pytuple(list(self[["_central_inertia"]], self[["_masscenter"]]))
}
.pvS_Body_inertia <- function(self, I) .pvS_RigidBody_inertia(self, I)

.pm_body_frame_of <- function(body) {
  if (.pv_is_frame(body)) return(body)
  if (inherits(body, "Body")) return(body$frame)
  .pv_err("UnboundLocalError", "cannot access local variable 'frame' where it is not associated with a value")
}

.pvM_Body_kinetic_energy <- function(self, frame) {
  if (.pvP_Body_is_rigidbody(self))
    return(RigidBody(self$name, self$masscenter, self$frame, self$mass,
                     .pytuple(list(self$central_inertia, self$masscenter)))$kinetic_energy(frame))
  Particle(self$name, self$masscenter, self$mass)$kinetic_energy(frame)
}

.pm_load_has <- function(load, obj) any(vapply(load, function(z) identical(z, obj), TRUE))

.pvM_Body_apply_force <- function(self, force, point = NULL, reaction_body = NULL, reaction_point = NULL) {
  if (!.pv_is_point(point)) {
    if (is.null(point)) point <- self$masscenter
    else .pv_err("TypeError", "Force must be applied to a point on the body.")
  }
  if (!.pv_is_vec(force)) .pv_err("TypeError", "Force must be a vector.")
  if (!is.null(reaction_body)) reaction_body$apply_force(-force, point = reaction_point)
  loads <- self[["_loads"]]
  for (i in seq_along(loads)) {
    if (.pm_load_has(loads[[i]], point)) {
      force <- force + loads[[i]][[2L]]
      loads[[i]] <- NULL
      break
    }
  }
  loads[[length(loads) + 1L]] <- .pytuple(list(point, force))
  self[["_loads"]] <- loads
  invisible(NULL)
}

.pvM_Body_apply_torque <- function(self, torque, reaction_body = NULL) {
  if (!.pv_is_vec(torque)) .pv_err("TypeError", "A Vector must be supplied to add torque.")
  if (!is.null(reaction_body)) reaction_body$apply_torque(-torque)
  loads <- self[["_loads"]]
  frame <- self$frame
  for (i in seq_along(loads)) {
    if (.pm_load_has(loads[[i]], frame)) {
      torque <- torque + loads[[i]][[2L]]
      loads[[i]] <- NULL
      break
    }
  }
  loads[[length(loads) + 1L]] <- .pytuple(list(frame, torque))
  self[["_loads"]] <- loads
  invisible(NULL)
}

.pvM_Body_clear_loads <- function(self) {
  self[["_loads"]] <- list()
  invisible(NULL)
}

.pvM_Body_remove_load <- function(self, about = NULL) {
  if (!is.null(about)) {
    if (!.pv_is_point(about)) .pv_err("TypeError", "Load is applied about Point or ReferenceFrame.")
  } else about <- self$frame
  loads <- self[["_loads"]]
  for (i in seq_along(loads)) {
    if (.pm_load_has(loads[[i]], about)) {
      loads[[i]] <- NULL
      break
    }
  }
  self[["_loads"]] <- loads
  invisible(NULL)
}

.pvM_Body_masscenter_vel <- function(self, body) .pv_point_vel(self$masscenter, .pm_body_frame_of(body))
.pvM_Body_ang_vel_in <- function(self, body) .pv_ang_vel_in(self$frame, .pm_body_frame_of(body))
.pvM_Body_dcm <- function(self, body) .pv_dcm(self$frame, .pm_body_frame_of(body))

.pvM_Body_parallel_axis <- function(self, point, frame = NULL) {
  if (.pvP_Body_is_rigidbody(self)) return(.pvM_RigidBody_parallel_axis(self, point, frame))
  if (is.null(frame)) .pv_err("AttributeError", "'NoneType' object has no attribute 'x'")
  .pvM_Particle_parallel_axis(self, point, frame)
}

# ---------------------------------------------------------------------------
# jointsmethod.py

#' Equations of motion of a set of joints
#'
#' `JointsMethod()` is SymPy's deprecated `JointsMethod`: it collects the
#' bodies, loads, coordinates, speeds and kinematic equations of the given
#' joints, and `form_eoms(method)` (reached with `$`) forms the equations
#' of motion with [KanesMethod()] (the default) or [LagrangesMethod()].
#' New code should use [System()].
#'
#' @param newtonion The Newtonian frame, or a body whose frame is used.
#' @param ... Joints.
#' @return A `JointsMethod` object.
#' @examples
#' P <- RigidBody("P")
#' C <- RigidBody("C")
#' J <- PinJoint("J", P, C)
#' method <- suppressWarnings(JointsMethod(P, J))
#' method$q
#' @export
JointsMethod <- function(newtonion, ...) {
  warning("JointsMethod is deprecated; use System instead.", call. = FALSE)
  self <- new.env(parent = emptyenv())
  class(self) <- c("JointsMethod", "MethodBase", "physics_obj")
  assign("_uid", .pv_next_uid(), envir = self)
  self[["_frame"]] <- if (.pm_is_body(newtonion)) newtonion$frame else newtonion
  self[["_joints"]] <- .pytuple(list(...))
  bodies <- list()
  for (joint in self[["_joints"]]) {
    for (b in list(joint$child, joint$parent))
      if (!any(vapply(bodies, function(z) identical(z, b), TRUE))) bodies[[length(bodies) + 1L]] <- b
  }
  self[["_bodies"]] <- bodies
  loads <- list()
  for (body in bodies) if (inherits(body, "Body")) loads <- c(loads, body$loads)
  self[["_loads"]] <- loads
  gen <- function(attr, msg) {
    out <- list()
    for (joint in self[["_joints"]]) for (c in joint[[attr]]$data) {
      if (any(vapply(out, function(z) identical(z$key, c$key), TRUE))) .pv_err("ValueError", msg)
      out[[length(out) + 1L]] <- c
    }
    .pm_colmat(out)
  }
  self[["_q"]] <- gen("_coordinates", "Coordinates of joints should be unique.")
  self[["_u"]] <- gen("_speeds", "Speeds of joints should be unique.")
  kd <- .pm_empty(0L, 1L)
  for (joint in self[["_joints"]]) kd <- .pm_cjoin(kd, joint$kdes)
  self[["_kdes"]] <- kd
  self[["_method"]] <- NULL
  self
}

.pvP_JointsMethod_bodies <- function(self) self[["_bodies"]]
.pvP_JointsMethod_loads <- function(self) self[["_loads"]]
.pvP_JointsMethod_q <- function(self) self[["_q"]]
.pvP_JointsMethod_u <- function(self) self[["_u"]]
.pvP_JointsMethod_kdes <- function(self) self[["_kdes"]]
.pvP_JointsMethod_method <- function(self) self[["_method"]]
.pvP_JointsMethod_frame <- function(self) self[["_frame"]]
.pvP_JointsMethod_forcing_full <- function(self) self[["_method"]]$forcing_full
.pvP_JointsMethod_mass_matrix_full <- function(self) self[["_method"]]$mass_matrix_full
.pvP_JointsMethod_mass_matrix <- function(self) self[["_method"]]$mass_matrix
.pvP_JointsMethod_forcing <- function(self) self[["_method"]]$forcing
.pvP_JointsMethod_holonomic_constraints <- function(self) .pm_empty()
.pvP_JointsMethod_nonholonomic_constraints <- function(self) .pm_empty()
.pvP_JointsMethod_velocity_constraints <- function(self) .pm_empty()
.pvP_JointsMethod_acceleration_constraints <- function(self) .pm_empty()

.pvM_JointsMethod__convert_bodies <- function(self) {
  lapply(self[["_bodies"]], function(body) {
    if (!inherits(body, "Body")) return(body)
    if (.pvP_Body_is_rigidbody(body)) {
      rb <- RigidBody(body$name, body$masscenter, body$frame, body$mass, .pytuple(list(body$central_inertia, body$masscenter)))
    } else {
      rb <- Particle(body$name, body$masscenter, body$mass)
    }
    rb$potential_energy <- body$potential_energy
    rb
  })
}

.pvM_JointsMethod_form_eoms <- function(self, method = KanesMethod) {
  bodylist <- .pvM_JointsMethod__convert_bodies(self)
  if (identical(method, LagrangesMethod)) {
    L <- do.call(Lagrangian, c(list(self[["_frame"]]), bodylist))
    self[["_method"]] <- method(L, self[["_q"]], self[["_loads"]], bodylist, self[["_frame"]])
  } else {
    self[["_method"]] <- method(self[["_frame"]], q_ind = self[["_q"]], u_ind = self[["_u"]], kd_eqs = self[["_kdes"]],
                                forcelist = self[["_loads"]], bodies = bodylist)
  }
  self[["_method"]]$`_form_eoms`()
}
.pvM_JointsMethod__form_eoms <- function(self) .pvM_JointsMethod_form_eoms(self)

.pvM_JointsMethod_rhs <- function(self, inv_method = NULL) self[["_method"]]$rhs(inv_method = inv_method)
