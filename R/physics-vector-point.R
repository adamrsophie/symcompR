# Points of classical mechanics, after sympy/physics/vector/point.py.

#' Points of classical mechanics
#'
#' `PhysPoint()` creates a point of `sympy.physics.vector` (SymPy's
#' `Point`; the name differs from the geometry `Point()` of this package).
#' Points are mutable: `set_pos`, `set_vel` and `set_acc` record position,
#' velocity and acceleration, `pos_from`, `vel` and `acc` query them, and
#' `locatenew`, `v1pt_theory`, `v2pt_theory`, `a1pt_theory`,
#' `a2pt_theory` and `partial_velocity` follow SymPy. All are reached with
#' `$`.
#'
#' @param name Point name.
#' @return A `PhysPoint`.
#' @examples
#' N <- ReferenceFrame("N")
#' O <- PhysPoint("O")
#' q <- dynamicsymbols("q")
#' P <- O$locatenew("P", q * N$x)
#' O$set_vel(N, 0)
#' P$vel(N)
#' @export
PhysPoint <- function(name) {
  p <- new.env(parent = emptyenv())
  assign("_uid", .pv_next_uid(), envir = p)
  p$name <- name
  p[["_pos_dict"]] <- list()
  p[["_vel_dict"]] <- list()
  p[["_acc_dict"]] <- list()
  class(p) <- c("PhysPoint", "physics_obj")
  p
}

#' @export
format.PhysPoint <- function(x, ...) .pv_str(get("name", envir = x))

.pv_point_dict <- function(p, num) get(c("_pos_dict", "_vel_dict", "_acc_dict")[num + 1L], envir = p)

.pv_check_point <- function(other) if (!.pv_is_point(other)) .pv_err("TypeError", "A Point must be supplied")

.pv_pdict_list <- function(self, other, num) {
  p <- .pv_connect_path(self, other, function(pt) .pv_point_dict(pt, num))
  if (!is.null(p)) return(p)
  .pv_err("ValueError", paste0("No Connecting Path found between ", .pv_str(other$name), " and ", .pv_str(self$name)))
}
.pvM_PhysPoint__pdict_list <- .pv_pdict_list

.pv_pos_from <- function(self, otherpoint) {
  outvec <- .pv_vec0()
  plist <- .pv_pdict_list(self, otherpoint, 0L)
  for (i in seq_len(length(plist) - 1L))
    outvec <- .pv_vadd(outvec, .pvod_get(plist[[i]][["_pos_dict"]], plist[[i + 1L]]))
  outvec
}
.pvM_PhysPoint_pos_from <- .pv_pos_from

.pvM_PhysPoint_set_acc <- function(self, frame, value) {
  if (.pv_is0(value)) value <- .pv_vec0()
  value <- .pv_check_vector(value)
  .pv_check_frame(frame)
  self[["_acc_dict"]] <- .pvod_set(self[["_acc_dict"]], frame, value)
  invisible(NULL)
}

.pvM_PhysPoint_set_pos <- function(self, otherpoint, value) {
  if (.pv_is0(value)) value <- .pv_vec0()
  value <- .pv_check_vector(value)
  .pv_check_point(otherpoint)
  self[["_pos_dict"]] <- .pvod_set(self[["_pos_dict"]], otherpoint, value)
  otherpoint[["_pos_dict"]] <- .pvod_set(otherpoint[["_pos_dict"]], self, .pv_vneg(value))
  invisible(NULL)
}

.pvM_PhysPoint_set_vel <- function(self, frame, value) {
  if (.pv_is0(value)) value <- .pv_vec0()
  value <- .pv_check_vector(value)
  .pv_check_frame(frame)
  self[["_vel_dict"]] <- .pvod_set(self[["_vel_dict"]], frame, value)
  invisible(NULL)
}

.pvM_PhysPoint_locatenew <- function(self, name, value) {
  if (!is.character(name)) .pv_err("TypeError", "Must supply a valid name")
  if (.pv_is0(value)) value <- .pv_vec0()
  value <- .pv_check_vector(value)
  p <- PhysPoint(name)
  .pvM_PhysPoint_set_pos(p, self, value)
  .pvM_PhysPoint_set_pos(self, p, .pv_vneg(value))
  p
}

.pv_point_vel <- function(self, frame) {
  .pv_check_frame(frame)
  vd <- self[["_vel_dict"]]
  if (!.pvod_has(vd, frame)) {
    valid <- FALSE
    is_cyclic <- FALSE
    visited <- list()
    queue <- list(self)
    candidates <- list()
    while (length(queue)) {
      node <- queue[[1L]]
      queue <- queue[-1L]
      if (.pv_in(node, visited)) next
      visited[[length(visited) + 1L]] <- node
      for (ent in node[["_pos_dict"]]) {
        neighbor <- ent$key
        if (.pv_in(neighbor, visited)) next
        ok <- tryCatch({ .pv_express(ent$val, frame); TRUE },
                       ValueError = function(e) FALSE)
        if (!ok) next
        if (.pv_in(neighbor, queue)) is_cyclic <- TRUE
        if (!.pvod_has(neighbor[["_vel_dict"]], frame)) {
          queue[[length(queue) + 1L]] <- neighbor
          next
        }
        nv <- .pvod_get(neighbor[["_vel_dict"]], frame)
        candidates[[length(candidates) + 1L]] <- neighbor
        if (!valid) {
          .pvM_PhysPoint_set_vel(self, frame, .pv_vadd(.pv_time_derivative(.pv_pos_from(self, neighbor), frame), nv))
          valid <- TRUE
        }
      }
    }
    if (is_cyclic) warning(.pv_fill("Kinematic loops are defined among the positions of points. This is likely not desired and may cause errors in your calculations."), call. = FALSE)
    if (length(candidates) > 1L)
      warning(.pv_fill(sprintf("Velocity of %s automatically calculated based on point %s but it is also possible from points(s): %s. Velocities from these points are not necessarily the same. This may cause errors in your calculations.",
                               self$name, candidates[[1L]]$name, .pv_str(candidates[-1L]))), call. = FALSE)
    if (valid) return(.pvod_get(self[["_vel_dict"]], frame))
    .pv_err("ValueError", .pv_fill(sprintf("Velocity of point %s has not been defined in ReferenceFrame %s.", self$name, frame$name)))
  }
  .pvod_get(vd, frame)
}
.pvM_PhysPoint_vel <- .pv_point_vel

.pv_point_acc <- function(self, frame) {
  .pv_check_frame(frame)
  if (!.pvod_has(self[["_acc_dict"]], frame)) {
    if (!.pv_vec_eq(.pv_point_vel(self, frame), 0))
      return(.pv_time_derivative(.pvod_get(self[["_vel_dict"]], frame), frame))
    return(.pv_vec0())
  }
  .pvod_get(self[["_acc_dict"]], frame)
}
.pvM_PhysPoint_acc <- .pv_point_acc

.pvM_PhysPoint_a1pt_theory <- function(self, otherpoint, outframe, interframe) {
  .pv_check_frame(outframe)
  .pv_check_frame(interframe)
  .pv_check_point(otherpoint)
  dist <- .pv_pos_from(self, otherpoint)
  v <- .pv_point_vel(self, interframe)
  a1 <- .pv_point_acc(otherpoint, outframe)
  a2 <- .pv_point_acc(self, interframe)
  omega <- .pv_ang_vel_in(interframe, outframe)
  alpha <- .pv_ang_acc_in(interframe, outframe)
  val <- .pv_vadd(.pv_vadd(.pv_vadd(.pv_vadd(a2, .pv_vmul(.pv_vcross(omega, v), 2)), a1),
                           .pv_vcross(alpha, dist)), .pv_vcross(omega, .pv_vcross(omega, dist)))
  .pvM_PhysPoint_set_acc(self, outframe, val)
  .pv_point_acc(self, outframe)
}

.pvM_PhysPoint_a2pt_theory <- function(self, otherpoint, outframe, fixedframe) {
  .pv_check_frame(outframe)
  .pv_check_frame(fixedframe)
  .pv_check_point(otherpoint)
  dist <- .pv_pos_from(self, otherpoint)
  a <- .pv_point_acc(otherpoint, outframe)
  omega <- .pv_ang_vel_in(fixedframe, outframe)
  alpha <- .pv_ang_acc_in(fixedframe, outframe)
  .pvM_PhysPoint_set_acc(self, outframe, .pv_vadd(.pv_vadd(a, .pv_vcross(alpha, dist)), .pv_vcross(omega, .pv_vcross(omega, dist))))
  .pv_point_acc(self, outframe)
}

.pvM_PhysPoint_v1pt_theory <- function(self, otherpoint, outframe, interframe) {
  .pv_check_frame(outframe)
  .pv_check_frame(interframe)
  .pv_check_point(otherpoint)
  dist <- .pv_pos_from(self, otherpoint)
  v1 <- .pv_point_vel(self, interframe)
  v2 <- .pv_point_vel(otherpoint, outframe)
  omega <- .pv_ang_vel_in(interframe, outframe)
  .pvM_PhysPoint_set_vel(self, outframe, .pv_vadd(.pv_vadd(v1, v2), .pv_vcross(omega, dist)))
  .pv_point_vel(self, outframe)
}

.pvM_PhysPoint_v2pt_theory <- function(self, otherpoint, outframe, fixedframe) {
  .pv_check_frame(outframe)
  .pv_check_frame(fixedframe)
  .pv_check_point(otherpoint)
  dist <- .pv_pos_from(self, otherpoint)
  v <- .pv_point_vel(otherpoint, outframe)
  omega <- .pv_ang_vel_in(fixedframe, outframe)
  .pvM_PhysPoint_set_vel(self, outframe, .pv_vadd(v, .pv_vcross(omega, dist)))
  .pv_point_vel(self, outframe)
}

.pvM_PhysPoint_partial_velocity <- function(self, frame, ...) {
  vel <- .pv_point_vel(self, frame)
  partials <- .pv_partial_velocity(list(vel), list(...), frame)[[1L]]
  if (length(partials) == 1L) partials[[1L]] else .pytuple(partials)
}

.pvP_PhysPoint__pos_dict <- function(self) .od_pydict(self[["_pos_dict"]])
.pvP_PhysPoint__vel_dict <- function(self) .od_pydict(self[["_vel_dict"]])
.pvP_PhysPoint__acc_dict <- function(self) .od_pydict(self[["_acc_dict"]])
# the dictionaries are plain attributes in SymPy, so callers may replace them
.pv_od_from_pydict <- function(d) {
  out <- list()
  if (inherits(d, "pydict")) for (i in seq_along(d$keys)) out <- .pvod_set(out, d$keys[[i]], d$values[[i]])
  out
}
.pvS_PhysPoint__pos_dict <- function(self, value) self[["_pos_dict"]] <- .pv_od_from_pydict(value)
.pvS_PhysPoint__vel_dict <- function(self, value) self[["_vel_dict"]] <- .pv_od_from_pydict(value)
.pvS_PhysPoint__acc_dict <- function(self, value) self[["_acc_dict"]] <- .pv_od_from_pydict(value)
