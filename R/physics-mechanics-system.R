# The System and SymbolicSystem classes, after
# sympy/physics/mechanics/system.py.

# System._objects_to_list
.pm_objects_to_list <- function(lst) {
  if (inherits(lst, "symmatrix")) return(lst$data)
  if (inherits(lst, "pynamedtuple")) return(unclass(lst))
  if (is.list(lst) && !is_sym(lst) && !inherits(lst, "physics_obj")) return(unclass(lst))
  list(lst)
}

.pm_obj_in <- function(x, xs) any(vapply(xs, function(z) .pm_py_eq(z, x), TRUE))

.pm_obj_set_str <- function(xs) paste0("{", paste(vapply(xs, .pv_objrepr, ""), collapse = ", "), "}")

.PM_EXPECTED <- list(BodyBase = function(x) .pm_is_body(x), Joint = function(x) .pm_is_joint(x),
                     ActuatorBase = function(x) inherits(x, "ActuatorBase"), Basic = function(x) is_sym(x))

# System._check_objects
.pm_check_objects <- function(objects, obj_lst, expected_type, obj_name, type_name) {
  seen <- unclass(obj_lst)
  duplicates <- list()
  wrong_types <- list()
  ok <- .PM_EXPECTED[[expected_type]]
  for (obj in objects) {
    if (!ok(obj) && !.pm_obj_in(obj, wrong_types)) wrong_types[[length(wrong_types) + 1L]] <- obj
    if (.pm_obj_in(obj, seen)) {
      if (!.pm_obj_in(obj, duplicates)) duplicates[[length(duplicates) + 1L]] <- obj
    } else seen[[length(seen) + 1L]] <- obj
  }
  if (length(wrong_types))
    .pv_err("TypeError", sprintf("%s %s are not %s.", obj_name, .pm_obj_set_str(wrong_types), type_name))
  if (length(duplicates))
    .pv_err("ValueError", sprintf("%s %s have already been added to the system.", obj_name, .pm_obj_set_str(duplicates)))
  invisible(NULL)
}

.pm_colvec <- function(xs) .new_mat(lapply(xs, as_sym), length(xs), 1L)

# System._parse_coordinates
.pm_parse_coordinates <- function(self, new_coords, independent, old_coords_ind, old_coords_dep,
                                  coord_type = "coordinates") {
  coords_ind <- .pm_iterable(old_coords_ind)
  coords_dep <- .pm_iterable(old_coords_dep)
  new_coords <- .pm_objects_to_list(new_coords)
  if (!is.list(independent) && length(independent) <= 1L) independent <- rep(list(independent), length(new_coords))
  independent <- as.list(independent)
  for (i in seq_len(min(length(new_coords), length(independent)))) {
    if (isTRUE(as.logical(independent[[i]]))) coords_ind[[length(coords_ind) + 1L]] <- new_coords[[i]]
    else coords_dep[[length(coords_dep) + 1L]] <- new_coords[[i]]
  }
  current <- list(coordinates = c(self[["_q_ind"]]$data, self[["_q_dep"]]$data),
                  speeds = c(self[["_u_ind"]]$data, self[["_u_dep"]]$data),
                  u_auxiliary = self[["_u_aux"]]$data)
  current[[coord_type]] <- c(coords_ind, coords_dep)
  .pm_validate_coordinates(current$coordinates, current$speeds, u_auxiliary = current$u_auxiliary)
  list(.pm_colvec(coords_ind), .pm_colvec(coords_dep))
}

# System._parse_expressions
.pm_parse_expressions <- function(new_expressions, old_expressions, name, check_negatives = FALSE) {
  old_expressions <- .pm_iterable(old_expressions)
  new_expressions <- .pm_objects_to_list(new_expressions)
  check_exprs <- if (check_negatives) c(old_expressions, lapply(old_expressions, .neg)) else old_expressions
  .pm_check_objects(new_expressions, check_exprs, "Basic", name, "expressions")
  for (expr in new_expressions)
    if (identical(expr$key, S_Zero$key)) .pv_err("ValueError", sprintf("Parsed %s are zero.", name))
  .pm_colvec(c(old_expressions, new_expressions))
}

#' Multibody systems
#'
#' `System()` is SymPy's `System` class: a container of bodies, joints,
#' loads, actuators, generalized coordinates and speeds, kinematic
#' differential equations and constraints, which forms the equations of
#' motion with `form_eoms()` using [KanesMethod()] (the default) or
#' [LagrangesMethod()]. `System_from_newtonian()` is SymPy's
#' `System.from_newtonian()`. Attributes, setters (`$<-`) and the methods
#' `add_coordinates()`, `add_speeds()`, `add_auxiliary_speeds()`,
#' `add_kdes()`, `add_holonomic_constraints()`,
#' `add_nonholonomic_constraints()`, `add_bodies()`, `add_loads()`,
#' `apply_uniform_gravity()`, `add_actuators()`, `add_joints()`,
#' `get_body()`, `get_joint()`, `form_eoms()`, `rhs()` and
#' `validate_system()` are reached with `$`. Methods that take a variable
#' number of objects in SymPy take them through `...`; `independent` is a
#' named argument.
#'
#' @param frame The inertial [ReferenceFrame()]; `inertial_frame` when
#'   `NULL`.
#' @param fixed_point The fixed [PhysPoint()]; `inertial_point` when `NULL`.
#' @param newtonian A rigid body whose frame and mass center are used.
#' @return A `System` object.
#' @examples
#' sys <- System()
#' q <- dynamicsymbols("q"); u <- dynamicsymbols("u")
#' sys$add_coordinates(q)
#' sys$add_speeds(u)
#' sys$add_kdes(u - dynamicsymbols("q", 1))
#' sys$q
#' sys$kdes
#' @export
System <- function(frame = NULL, fixed_point = NULL) {
  self <- new.env(parent = emptyenv())
  class(self) <- c("System", "MethodBase", "physics_obj")
  assign("_uid", .pv_next_uid(), envir = self)
  if (is.null(frame)) frame <- ReferenceFrame("inertial_frame")
  else if (!.pv_is_frame(frame)) .pv_err("TypeError", "Frame must be an instance of ReferenceFrame.")
  self[["_frame"]] <- frame
  if (is.null(fixed_point)) fixed_point <- PhysPoint("inertial_point")
  else if (!.pv_is_point(fixed_point)) .pv_err("TypeError", "Fixed point must be an instance of Point.")
  self[["_fixed_point"]] <- fixed_point
  fixed_point$set_vel(frame, 0)
  for (nm in c("_q_ind", "_q_dep", "_u_ind", "_u_dep", "_u_aux", "_kdes", "_hol_coneqs", "_nonhol_coneqs"))
    self[[nm]] <- .pm_empty(0L, 1L)
  self[["_vel_constrs"]] <- NULL
  self[["_bodies"]] <- list()
  self[["_joints"]] <- list()
  self[["_loads"]] <- list()
  self[["_actuators"]] <- list()
  self[["_eom_method"]] <- NULL
  self
}

#' @rdname System
#' @export
System_from_newtonian <- function(newtonian) {
  if (inherits(newtonian, "Particle"))
    .pv_err("TypeError", "A Particle has no frame so cannot act as the Newtonian.")
  system <- System(frame = newtonian$frame, fixed_point = newtonian$masscenter)
  .pvM_System_add_bodies(system, newtonian)
  system
}
.pvM_System_from_newtonian <- function(self, newtonian) System_from_newtonian(newtonian)

.pm_reset <- function(self) self[["_eom_method"]] <- NULL

.pvP_System_fixed_point <- function(self) self[["_fixed_point"]]
.pvP_System_frame <- function(self) self[["_frame"]]
.pvP_System_x <- function(self) self[["_frame"]]$x
.pvP_System_y <- function(self) self[["_frame"]]$y
.pvP_System_z <- function(self) self[["_frame"]]$z

.pvP_System_bodies <- function(self) .pytuple(self[["_bodies"]])
.pvS_System_bodies <- function(self, bodies) {
  .pm_reset(self)
  bodies <- .pm_objects_to_list(bodies)
  .pm_check_objects(bodies, list(), "BodyBase", "Bodies", "bodies")
  self[["_bodies"]] <- bodies
}
.pvP_System_joints <- function(self) .pytuple(self[["_joints"]])
.pvS_System_joints <- function(self, joints) {
  .pm_reset(self)
  joints <- .pm_objects_to_list(joints)
  .pm_check_objects(joints, list(), "Joint", "Joints", "joints")
  self[["_joints"]] <- list()
  do.call(.pvM_System_add_joints, c(list(self), joints))
}
.pvP_System_loads <- function(self) .pytuple(self[["_loads"]])
.pvS_System_loads <- function(self, loads) {
  .pm_reset(self)
  loads <- .pm_objects_to_list(loads)
  self[["_loads"]] <- lapply(loads, .pm_parse_load)
}
.pvP_System_actuators <- function(self) .pytuple(self[["_actuators"]])
.pvS_System_actuators <- function(self, actuators) {
  .pm_reset(self)
  actuators <- .pm_objects_to_list(actuators)
  .pm_check_objects(actuators, list(), "ActuatorBase", "Actuators", "actuators")
  self[["_actuators"]] <- actuators
}

.pvP_System_q <- function(self) .pm_cjoin(self[["_q_ind"]], self[["_q_dep"]])
.pvP_System_u <- function(self) .pm_cjoin(self[["_u_ind"]], self[["_u_dep"]])
.pvP_System_q_ind <- function(self) self[["_q_ind"]]
.pvP_System_q_dep <- function(self) self[["_q_dep"]]
.pvP_System_u_ind <- function(self) self[["_u_ind"]]
.pvP_System_u_dep <- function(self) self[["_u_dep"]]
.pvP_System_u_aux <- function(self) self[["_u_aux"]]

.pm_set_pair <- function(self, a, b, r) {
  self[[a]] <- r[[1L]]
  self[[b]] <- r[[2L]]
}
.pvS_System_q_ind <- function(self, q_ind) {
  .pm_reset(self)
  .pm_set_pair(self, "_q_ind", "_q_dep", .pm_parse_coordinates(self, .pm_objects_to_list(q_ind), TRUE, list(),
                                                               self[["_q_dep"]], "coordinates"))
}
.pvS_System_q_dep <- function(self, q_dep) {
  .pm_reset(self)
  .pm_set_pair(self, "_q_ind", "_q_dep", .pm_parse_coordinates(self, .pm_objects_to_list(q_dep), FALSE,
                                                               self[["_q_ind"]], list(), "coordinates"))
}
.pvS_System_u_ind <- function(self, u_ind) {
  .pm_reset(self)
  .pm_set_pair(self, "_u_ind", "_u_dep", .pm_parse_coordinates(self, .pm_objects_to_list(u_ind), TRUE, list(),
                                                               self[["_u_dep"]], "speeds"))
}
.pvS_System_u_dep <- function(self, u_dep) {
  .pm_reset(self)
  .pm_set_pair(self, "_u_ind", "_u_dep", .pm_parse_coordinates(self, .pm_objects_to_list(u_dep), FALSE,
                                                               self[["_u_ind"]], list(), "speeds"))
}
.pvS_System_u_aux <- function(self, u_aux) {
  .pm_reset(self)
  self[["_u_aux"]] <- .pm_parse_coordinates(self, .pm_objects_to_list(u_aux), TRUE, list(), list(), "u_auxiliary")[[1L]]
}

.pvP_System_kdes <- function(self) self[["_kdes"]]
.pvS_System_kdes <- function(self, kdes) {
  .pm_reset(self)
  self[["_kdes"]] <- .pm_parse_expressions(.pm_objects_to_list(kdes), list(), "kinematic differential equations")
}
.pvP_System_holonomic_constraints <- function(self) self[["_hol_coneqs"]]
.pvS_System_holonomic_constraints <- function(self, constraints) {
  .pm_reset(self)
  self[["_hol_coneqs"]] <- .pm_parse_expressions(.pm_objects_to_list(constraints), list(), "holonomic constraints")
}
.pvP_System_nonholonomic_constraints <- function(self) self[["_nonhol_coneqs"]]
.pvS_System_nonholonomic_constraints <- function(self, constraints) {
  .pm_reset(self)
  self[["_nonhol_coneqs"]] <- .pm_parse_expressions(.pm_objects_to_list(constraints), list(), "nonholonomic constraints")
}
.pvP_System_velocity_constraints <- function(self) {
  if (is.null(self[["_vel_constrs"]]))
    return(.pm_cjoin(.pv_mdiff(self[["_hol_coneqs"]], .pv_t()), self[["_nonhol_coneqs"]]))
  self[["_vel_constrs"]]
}
.pvS_System_velocity_constraints <- function(self, constraints) {
  .pm_reset(self)
  if (is.null(constraints)) {
    self[["_vel_constrs"]] <- NULL
    return(invisible(NULL))
  }
  self[["_vel_constrs"]] <- .pm_parse_expressions(.pm_objects_to_list(constraints), list(), "velocity constraints")
}
.pvP_System_acceleration_constraints <- function(self) .pv_mdiff(self$velocity_constraints, .pv_t())
.pvP_System_eom_method <- function(self) self[["_eom_method"]]

.pm_independent <- function(independent) if (is.list(independent)) independent else as.list(independent)

.pvM_System_add_coordinates <- function(self, ..., independent = TRUE) {
  .pm_reset(self)
  .pm_set_pair(self, "_q_ind", "_q_dep", .pm_parse_coordinates(self, list(...), independent, self[["_q_ind"]],
                                                               self[["_q_dep"]], "coordinates"))
  invisible(NULL)
}
.pvM_System_add_speeds <- function(self, ..., independent = TRUE) {
  .pm_reset(self)
  .pm_set_pair(self, "_u_ind", "_u_dep", .pm_parse_coordinates(self, list(...), independent, self[["_u_ind"]],
                                                               self[["_u_dep"]], "speeds"))
  invisible(NULL)
}
.pvM_System_add_auxiliary_speeds <- function(self, ...) {
  .pm_reset(self)
  self[["_u_aux"]] <- .pm_parse_coordinates(self, list(...), TRUE, self[["_u_aux"]], list(), "u_auxiliary")[[1L]]
  invisible(NULL)
}
.pvM_System_add_kdes <- function(self, ...) {
  .pm_reset(self)
  self[["_kdes"]] <- .pm_parse_expressions(list(...), self[["_kdes"]], "kinematic differential equations",
                                           check_negatives = TRUE)
  invisible(NULL)
}
.pvM_System_add_holonomic_constraints <- function(self, ...) {
  .pm_reset(self)
  self[["_hol_coneqs"]] <- .pm_parse_expressions(list(...), self[["_hol_coneqs"]], "holonomic constraints",
                                                 check_negatives = TRUE)
  invisible(NULL)
}
.pvM_System_add_nonholonomic_constraints <- function(self, ...) {
  .pm_reset(self)
  self[["_nonhol_coneqs"]] <- .pm_parse_expressions(list(...), self[["_nonhol_coneqs"]], "nonholonomic constraints",
                                                    check_negatives = TRUE)
  invisible(NULL)
}
.pvM_System_add_bodies <- function(self, ...) {
  .pm_reset(self)
  bodies <- list(...)
  .pm_check_objects(bodies, self[["_bodies"]], "BodyBase", "Bodies", "bodies")
  self[["_bodies"]] <- c(self[["_bodies"]], bodies)
  invisible(NULL)
}
.pvM_System_add_loads <- function(self, ...) {
  .pm_reset(self)
  loads <- lapply(list(...), .pm_parse_load)
  self[["_loads"]] <- c(self[["_loads"]], loads)
  invisible(NULL)
}
.pvM_System_apply_uniform_gravity <- function(self, acceleration) {
  .pm_reset(self)
  do.call(.pvM_System_add_loads, c(list(self), do.call(gravity, c(list(acceleration), self[["_bodies"]]))))
  invisible(NULL)
}
.pvM_System_add_actuators <- function(self, ...) {
  .pm_reset(self)
  actuators <- list(...)
  .pm_check_objects(actuators, self[["_actuators"]], "ActuatorBase", "Actuators", "actuators")
  self[["_actuators"]] <- c(self[["_actuators"]], actuators)
  invisible(NULL)
}

# OrderedSet update and difference, keyed by sympy keys or object identity
.pm_oset_add <- function(set, xs) {
  for (x in xs) if (!.pm_obj_in(x, set)) set[[length(set) + 1L]] <- x
  set
}
.pm_oset_diff <- function(set, other) Filter(function(x) !.pm_obj_in(x, other), set)

.pvM_System_add_joints <- function(self, ...) {
  .pm_reset(self)
  joints <- list(...)
  .pm_check_objects(joints, self[["_joints"]], "Joint", "Joints", "joints")
  self[["_joints"]] <- c(self[["_joints"]], joints)
  coordinates <- list(); speeds <- list(); kdes <- list(); bodies <- list()
  for (joint in joints) {
    coordinates <- .pm_oset_add(coordinates, joint$coordinates$data)
    speeds <- .pm_oset_add(speeds, joint$speeds$data)
    kdes <- .pm_oset_add(kdes, joint$kdes$data)
    bodies <- .pm_oset_add(bodies, list(joint$parent, joint$child))
  }
  coordinates <- .pm_oset_diff(coordinates, self$q$data)
  speeds <- .pm_oset_diff(speeds, self$u$data)
  kdes <- .pm_oset_diff(kdes, c(self[["_kdes"]]$data, lapply(self[["_kdes"]]$data, .neg)))
  bodies <- .pm_oset_diff(bodies, self[["_bodies"]])
  do.call(.pvM_System_add_coordinates, c(list(self), coordinates))
  do.call(.pvM_System_add_speeds, c(list(self), speeds))
  do.call(.pvM_System_add_kdes, c(list(self), Filter(function(k) !identical(k$key, S_Zero$key), kdes)))
  do.call(.pvM_System_add_bodies, c(list(self), bodies))
  invisible(NULL)
}

.pvM_System_get_body <- function(self, name) {
  for (body in self[["_bodies"]]) if (identical(body$name, name)) return(body)
  NULL
}
.pvM_System_get_joint <- function(self, name) {
  for (joint in self[["_joints"]]) if (identical(joint$name, name)) return(joint)
  NULL
}

.pvM_System__form_eoms <- function(self) .pvM_System_form_eoms(self)

.pm_kw_set_str <- function(ks) paste0("{", paste0("'", ks, "'", collapse = ", "), "}")

.pvM_System_form_eoms <- function(self, eom_method = KanesMethod, ...) {
  kwargs <- list(...)
  loads <- self[["_loads"]]
  for (act in self[["_actuators"]]) loads <- c(loads, act$to_loads())
  loads <- if (length(loads)) .pytuple(loads) else NULL
  if (identical(eom_method, KanesMethod)) {
    disallowed <- c("frame", "q_ind", "u_ind", "kd_eqs", "q_dependent", "u_dependent", "u_auxiliary",
                    "configuration_constraints", "velocity_constraints", "forcelist", "bodies")
    wrong <- intersect(names(kwargs), disallowed)
    if (length(wrong))
      .pv_err("ValueError", sprintf("The following keyword arguments are not allowed to be overwritten in KanesMethod: %s.",
                                    .pm_kw_set_str(wrong)))
    base <- list(frame = self[["_frame"]], q_ind = self[["_q_ind"]], u_ind = self[["_u_ind"]],
                 kd_eqs = self[["_kdes"]], q_dependent = self[["_q_dep"]], u_dependent = self[["_u_dep"]],
                 configuration_constraints = self[["_hol_coneqs"]], velocity_constraints = self$velocity_constraints,
                 u_auxiliary = self[["_u_aux"]], forcelist = loads, bodies = self$bodies, explicit_kinematics = FALSE)
    base[names(kwargs)] <- kwargs
    self[["_eom_method"]] <- do.call(eom_method, base)
  } else if (identical(eom_method, LagrangesMethod)) {
    disallowed <- c("frame", "qs", "forcelist", "bodies", "hol_coneqs", "nonhol_coneqs", "Lagrangian")
    wrong <- intersect(names(kwargs), disallowed)
    if (length(wrong))
      .pv_err("ValueError", sprintf("The following keyword arguments are not allowed to be overwritten in LagrangesMethod: %s.",
                                    .pm_kw_set_str(wrong)))
    base <- list(frame = self[["_frame"]], qs = self$q, forcelist = loads, bodies = self$bodies,
                 hol_coneqs = self[["_hol_coneqs"]], nonhol_coneqs = self[["_nonhol_coneqs"]])
    base[names(kwargs)] <- kwargs
    if (is.null(base$Lagrangian)) base$Lagrangian <- do.call(Lagrangian, c(list(base$frame), unclass(base$bodies)))
    self[["_eom_method"]] <- do.call(eom_method, base)
  } else {
    .pv_err("NotImplementedError", "The equations of motion method has not been implemented.")
  }
  self[["_eom_method"]]$`_form_eoms`()
}

.pvM_System_rhs <- function(self, inv_method = NULL, ...) self[["_eom_method"]]$rhs(inv_method = inv_method, ...)
.pvP_System_mass_matrix <- function(self) self[["_eom_method"]]$mass_matrix
.pvP_System_mass_matrix_full <- function(self) self[["_eom_method"]]$mass_matrix_full
.pvP_System_forcing <- function(self) self[["_eom_method"]]$forcing
.pvP_System_forcing_full <- function(self) self[["_eom_method"]]$forcing_full

.pvM_System_validate_system <- function(self, eom_method = KanesMethod, check_duplicates = FALSE) {
  msgs <- character()
  n_hc <- self[["_hol_coneqs"]]$nrow
  n_vc <- self$velocity_constraints$nrow
  n_q_dep <- self[["_q_dep"]]$nrow
  n_u_dep <- self[["_u_dep"]]$nrow
  q_set <- .pm_oset_add(list(), self$q$data)
  u_set <- .pm_oset_add(list(), self$u$data)
  n_q <- length(q_set)
  n_u <- length(u_set)
  if (n_q_dep != n_hc)
    msgs <- c(msgs, .pv_fill(sprintf("The number of dependent generalized coordinates %d should be equal to the number of holonomic constraints %d.",
                                     n_q_dep, n_hc)))
  missing_q <- list()
  for (joint in self[["_joints"]]) missing_q <- .pm_oset_add(missing_q, .pm_oset_diff(joint$coordinates$data, q_set))
  if (length(missing_q))
    msgs <- c(msgs, .pv_fill(sprintf("The generalized coordinates %s used in joints are not added to the system.",
                                     .pm_set_str(missing_q))))
  if (identical(eom_method, KanesMethod)) {
    n_kdes <- .pm_len(self[["_kdes"]])
    missing_kdes <- list(); missing_u <- list()
    kd_all <- c(self[["_kdes"]]$data, lapply(self[["_kdes"]]$data, .neg))
    for (joint in self[["_joints"]]) {
      missing_u <- .pm_oset_add(missing_u, .pm_oset_diff(joint$speeds$data, u_set))
      missing_kdes <- .pm_oset_add(missing_kdes, .pm_oset_diff(joint$kdes$data, kd_all))
    }
    if (length(missing_u))
      msgs <- c(msgs, .pv_fill(sprintf("The generalized speeds %s used in joints are not added to the system.",
                                       .pm_set_str(missing_u))))
    if (length(missing_kdes))
      msgs <- c(msgs, .pv_fill(sprintf("The kinematic differential equations %s used in joints are not added to the system.",
                                       .pm_set_str(missing_kdes))))
    if (n_u_dep != n_vc)
      msgs <- c(msgs, .pv_fill(sprintf("The number of dependent generalized speeds %d should be equal to the number of velocity constraints %d.",
                                       n_u_dep, n_vc)))
    if (n_q > n_u)
      msgs <- c(msgs, .pv_fill(sprintf("The number of generalized coordinates %d should be less than or equal to the number of generalized speeds %d.",
                                       n_q, n_u)))
    if (n_u != n_kdes)
      msgs <- c(msgs, .pv_fill(sprintf("The number of generalized speeds %d should be equal to the number of kinematic differential equations %d.",
                                       n_u, n_kdes)))
  } else if (identical(eom_method, LagrangesMethod)) {
    qd <- .pv_mdiff(self$q, .pv_t())$data
    not_qdots <- .pm_oset_diff(.pm_oset_add(list(), self$u$data), qd)
    for (joint in self[["_joints"]]) not_qdots <- .pm_oset_add(not_qdots, .pm_oset_diff(joint$speeds$data, qd))
    if (length(not_qdots))
      msgs <- c(msgs, .pv_fill(sprintf(paste("The generalized speeds %s are not supported by this method. Only derivatives",
                                             "of the generalized coordinates are supported. If these symbols are used in",
                                             "your expressions, then this will result in wrong equations of motion."),
                                       .pm_set_str(not_qdots))))
    if (.pm_truthy(self[["_u_aux"]]))
      msgs <- c(msgs, .pv_fill(sprintf(paste("This method does not support auxiliary speeds. If these symbols are used in",
                                             "your expressions, then this will result in wrong equations of motion. The",
                                             "auxiliary speeds are %s."), .pv_str(self[["_u_aux"]]))))
  } else {
    .pv_err("NotImplementedError", "The equations of motion method has not been implemented.")
  }
  if (isTRUE(check_duplicates)) {
    checks <- list(list("generalized coordinates", self$q$data), list("generalized speeds", self$u$data),
                   list("auxiliary speeds", self[["_u_aux"]]$data), list("bodies", self[["_bodies"]]),
                   list("joints", self[["_joints"]]))
    for (chk in checks) {
      seen <- list(); duplicates <- list()
      for (x in chk[[2L]]) {
        if (.pm_obj_in(x, seen)) duplicates <- .pm_oset_add(duplicates, list(x))
        else seen[[length(seen) + 1L]] <- x
      }
      if (length(duplicates))
        msgs <- c(msgs, .pv_fill(sprintf("The %s %s exist multiple times within the system.", chk[[1L]],
                                         .pm_obj_set_str(duplicates))))
    }
  }
  if (length(msgs)) .pv_err("ValueError", paste(msgs, collapse = "\n"))
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# SymbolicSystem

#' Symbolic equations of motion
#'
#' `SymbolicSystem()` stores equations of motion in one of SymPy's three
#' forms (explicit, implicit with a combined mass matrix, or implicit with
#' separate kinematic and dynamic parts), as SymPy's `SymbolicSystem`.
#' Attributes such as `coordinates`, `speeds`, `states`,
#' `comb_explicit_rhs`, `comb_implicit_mat` and `comb_implicit_rhs`, and
#' the methods `compute_explicit_form()`, `dynamic_symbols()` and
#' `constant_symbols()`, are reached with `$`.
#'
#' @param coord_states Coordinates, or all states when `speeds` is `NULL`.
#' @param right_hand_side The right hand side of the equations.
#' @param speeds Generalized speeds.
#' @param mass_matrix Mass matrix of the implicit forms.
#' @param coordinate_derivatives Explicit kinematic equations.
#' @param alg_con Indices of algebraic constraint equations.
#' @param output_eqns A named list of output equations.
#' @param coord_idxs,speed_idxs Zero based indices of the coordinates and
#'   speeds in `coord_states`.
#' @param bodies,loads Bodies and loads of the system.
#' @return A `SymbolicSystem` object.
#' @examples
#' x <- dynamicsymbols("x"); v <- dynamicsymbols("v")
#' s <- SymbolicSystem(list(x, v), Matrix(list(v, -x)))
#' s$states
#' s$comb_explicit_rhs
#' @export
SymbolicSystem <- function(coord_states, right_hand_side, speeds = NULL, mass_matrix = NULL,
                           coordinate_derivatives = NULL, alg_con = NULL, output_eqns = list(),
                           coord_idxs = NULL, speed_idxs = NULL, bodies = NULL, loads = NULL) {
  self <- new.env(parent = emptyenv())
  class(self) <- c("SymbolicSystem", "physics_obj")
  assign("_uid", .pv_next_uid(), envir = self)
  if (is.null(speeds)) {
    self[["_states"]] <- .pm_colmat(coord_states)
    cs <- .pm_iterable(coord_states)
    self[["_coordinates"]] <- if (is.null(coord_idxs)) NULL else .pm_colmat(lapply(coord_idxs, function(i) cs[[i + 1L]]))
    self[["_speeds"]] <- if (is.null(speed_idxs)) NULL else .pm_colmat(lapply(speed_idxs, function(i) cs[[i + 1L]]))
  } else {
    self[["_coordinates"]] <- .pm_colmat(coord_states)
    self[["_speeds"]] <- .pm_colmat(speeds)
    self[["_states"]] <- .pm_cjoin(self[["_coordinates"]], self[["_speeds"]])
  }
  if (!is.null(coordinate_derivatives)) {
    self[["_kin_explicit_rhs"]] <- coordinate_derivatives
    self[["_dyn_implicit_rhs"]] <- right_hand_side
    self[["_dyn_implicit_mat"]] <- mass_matrix
  } else if (!is.null(mass_matrix)) {
    self[["_comb_implicit_rhs"]] <- right_hand_side
    self[["_comb_implicit_mat"]] <- mass_matrix
  } else {
    self[["_comb_explicit_rhs"]] <- right_hand_side
  }
  if (!is.null(alg_con) && !is.null(coordinate_derivatives))
    alg_con <- lapply(alg_con, function(i) i + .pm_len(coordinate_derivatives))
  self[["_alg_con"]] <- alg_con
  self[["output_eqns"]] <- output_eqns
  if (!is.null(bodies) && !inherits(bodies, "pytuple")) bodies <- .pytuple(.pm_iterable(bodies))
  if (!is.null(loads) && !inherits(loads, "pytuple")) loads <- .pytuple(.pm_iterable(loads))
  self[["_bodies"]] <- bodies
  self[["_loads"]] <- loads
  self
}

#' @export
format.SymbolicSystem <- function(x, ...) "<SymbolicSystem object>"

.pm_ss_get <- function(self, nm) if (exists(nm, envir = self, inherits = FALSE)) self[[nm]] else NULL
.pm_ss_need <- function(self, nm, msg) {
  v <- .pm_ss_get(self, nm)
  if (is.null(v)) .pv_err("AttributeError", msg)
  v
}

.pvP_SymbolicSystem_coordinates <- function(self) .pm_ss_need(self, "_coordinates", "The coordinates were not specified.")
.pvP_SymbolicSystem_speeds <- function(self) .pm_ss_need(self, "_speeds", "The speeds were not specified.")
.pvP_SymbolicSystem_states <- function(self) self[["_states"]]
.pvP_SymbolicSystem_alg_con <- function(self) self[["_alg_con"]]
.pvP_SymbolicSystem_dyn_implicit_mat <- function(self)
  .pm_ss_need(self, "_dyn_implicit_mat", "dyn_implicit_mat is not specified for equations of motion form [1] or [2].")
.pvP_SymbolicSystem_dyn_implicit_rhs <- function(self)
  .pm_ss_need(self, "_dyn_implicit_rhs", "dyn_implicit_rhs is not specified for equations of motion form [1] or [2].")
.pvP_SymbolicSystem_comb_implicit_mat <- function(self) {
  v <- .pm_ss_get(self, "_comb_implicit_mat")
  if (!is.null(v)) return(v)
  dm <- .pm_ss_get(self, "_dyn_implicit_mat")
  if (is.null(dm)) .pv_err("AttributeError", "comb_implicit_mat is not specified for equations of motion form [1].")
  num_kin_eqns <- .pm_len(self[["_kin_explicit_rhs"]])
  num_dyn_eqns <- .pm_len(self[["_dyn_implicit_rhs"]])
  inter1 <- .pm_rjoin(.mech_eye(num_kin_eqns), .mech_zeros(num_kin_eqns, num_dyn_eqns))
  inter2 <- .pm_rjoin(.mech_zeros(num_dyn_eqns, num_kin_eqns), dm)
  self[["_comb_implicit_mat"]] <- .pm_cjoin(inter1, inter2)
  self[["_comb_implicit_mat"]]
}
.pvP_SymbolicSystem_comb_implicit_rhs <- function(self) {
  v <- .pm_ss_get(self, "_comb_implicit_rhs")
  if (!is.null(v)) return(v)
  dr <- .pm_ss_get(self, "_dyn_implicit_rhs")
  if (is.null(dr)) .pv_err("AttributeError", "comb_implicit_mat is not specified for equations of motion in form [1].")
  self[["_comb_implicit_rhs"]] <- .pm_cjoin(self[["_kin_explicit_rhs"]], dr)
  self[["_comb_implicit_rhs"]]
}
.pvM_SymbolicSystem_compute_explicit_form <- function(self) {
  if (!is.null(.pm_ss_get(self, "_comb_explicit_rhs"))) .pv_err("AttributeError", "comb_explicit_rhs is already formed.")
  inter1 <- tryCatch(.pvP_SymbolicSystem_kin_explicit_rhs(self), AttributeError = function(e) NULL)
  if (!is.null(inter1)) {
    inter2 <- .pm_LUsolve(self[["_dyn_implicit_mat"]], self[["_dyn_implicit_rhs"]])
    out <- .pm_cjoin(inter1, inter2)
  } else {
    out <- .pm_LUsolve(self[["_comb_implicit_mat"]], self[["_comb_implicit_rhs"]])
  }
  self[["_comb_explicit_rhs"]] <- out
  invisible(NULL)
}
.pvP_SymbolicSystem_comb_explicit_rhs <- function(self)
  .pm_ss_need(self, "_comb_explicit_rhs",
              "Please run .combute_explicit_form before attempting to access comb_explicit_rhs.")
.pvP_SymbolicSystem_kin_explicit_rhs <- function(self)
  .pm_ss_need(self, "_kin_explicit_rhs", "kin_explicit_rhs is not specified for equations of motion form [1] or [2].")

.pm_ss_eom_exprs <- function(self) {
  if (is.null(.pm_ss_get(self, "_comb_explicit_rhs")))
    return(c(self$comb_implicit_mat$data, self$comb_implicit_rhs$data))
  self[["_comb_explicit_rhs"]]$data
}

.pvM_SymbolicSystem_dynamic_symbols <- function(self) {
  out <- list()
  for (expr in .pm_ss_eom_exprs(self)) out <- .pm_oset_add(out, find_dynamicsymbols(expr))
  out <- .pm_oset_add(out, self[["_states"]]$data)
  .pytuple(out)
}
.pvM_SymbolicSystem_constant_symbols <- function(self) {
  out <- list()
  for (expr in .pm_ss_eom_exprs(self)) out <- .pm_oset_add(out, free_symbols(as_sym(expr)))
  t <- .pv_t()
  if (!.pm_obj_in(t, out)) .pv_err("KeyError", "t")
  .pytuple(Filter(function(z) !identical(z$key, t$key), out))
}
.pvP_SymbolicSystem_bodies <- function(self) .pm_ss_need(self, "_bodies", "bodies were not specified for the system.")
.pvP_SymbolicSystem_loads <- function(self) .pm_ss_need(self, "_loads", "loads were not specified for the system.")
