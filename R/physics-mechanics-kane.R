# Kane's method, after sympy/physics/mechanics/method.py and kane.py, with
# the small matrix helpers the equation of motion classes share.

# ---------------------------------------------------------------------------
# Matrix helpers: SymPy's Matrix() of a list, col_join, row_join, jacobian,
# negation and xreplace, building entries as SymPy does.

.pm_empty <- function(n = 0L, m = 0L) .new_mat(list(), as.integer(n), as.integer(m))
.pm_len <- function(M) if (is.null(M)) 0L else if (inherits(M, "symmatrix")) length(M$data) else length(M)
.pm_truthy <- function(M) !is.null(M) && .pm_len(M) > 0L

# Matrix(x) if x else Matrix()
.pm_none_handler <- function(x) {
  if (inherits(x, "symmatrix")) return(if (length(x$data)) x else .pm_empty())
  if (!.pm_truthy(x)) return(.pm_empty())
  .pm_colmat(x)
}

# Matrix(list): a column of the entries; a list of matrices is stacked
.pm_colmat <- function(x) {
  if (inherits(x, "symmatrix")) return(x)
  if (is_sym(x) || is.numeric(x)) return(.pv_colv(list(x)))
  x <- unclass(x)
  if (!length(x)) return(.pm_empty())
  if (all(vapply(x, function(z) inherits(z, "symmatrix"), TRUE))) return(.pm_vstack(x))
  if (all(vapply(x, function(z) is.list(z) && !is_sym(z) && !inherits(z, "symmatrix"), TRUE)))
    return(.from_rows(lapply(x, function(r) lapply(r, as_sym))))
  .pv_colv(x)
}

.pm_vstack <- function(mats) {
  mats <- Filter(function(M) M$nrow > 0L, lapply(mats, .pm_colmat))
  if (!length(mats)) return(.pm_empty())
  Reduce(.pm_cjoin, mats)
}

.pm_cjoin <- function(A, B) {
  if (A$nrow == 0L && A$ncol != B$ncol) A <- .pm_empty(0L, B$ncol)
  if (A$ncol != B$ncol)
    .pv_err("ShapeError", sprintf("The matrices have incompatible number of columns (%d and %d)", A$ncol, B$ncol))
  .new_mat(c(A$data, B$data), A$nrow + B$nrow, B$ncol)
}

.pm_rjoin <- function(A, B) {
  if (A$ncol == 0L && A$nrow != B$nrow) A <- .pm_empty(B$nrow, 0L)
  if (A$nrow != B$nrow)
    .pv_err("ShapeError", sprintf("The matrices have incompatible number of rows (%d and %d)", A$nrow, B$nrow))
  d <- list()
  for (i in seq_len(A$nrow)) {
    if (A$ncol) d <- c(d, A$data[(i - 1L) * A$ncol + seq_len(A$ncol)])
    if (B$ncol) d <- c(d, B$data[(i - 1L) * B$ncol + seq_len(B$ncol)])
  }
  .new_mat(d, A$nrow, A$ncol + B$ncol)
}

.pm_mneg <- function(M) .pv_mmap(M, .neg)
.mech_zeros <- function(n, m = 1L) .pv_zeros(as.integer(n), as.integer(m))
.mech_eye <- function(n) {
  d <- rep(list(S_Zero), n * n)
  for (i in seq_len(n)) d[[(i - 1L) * n + i]] <- S_One
  .new_mat(d, as.integer(n), as.integer(n))
}

# rows lo..hi (zero based, hi exclusive) and all columns
.pm_rows <- function(M, lo, hi) {
  hi <- min(hi, M$nrow)
  if (hi <= lo) return(.pm_empty(0L, M$ncol))
  .new_mat(M$data[(lo * M$ncol + 1L):(hi * M$ncol)], hi - lo, M$ncol)
}
# columns lo..hi (zero based, hi exclusive) and all rows
.pm_cols <- function(M, lo, hi) {
  hi <- min(hi, M$ncol)
  if (hi <= lo) return(.pm_empty(M$nrow, 0L))
  d <- list()
  for (i in seq_len(M$nrow)) d <- c(d, M$data[(i - 1L) * M$ncol + seq.int(lo + 1L, hi)])
  .new_mat(d, M$nrow, hi - lo)
}

.pm_jacobian <- function(M, X) {
  X <- .pm_colmat(X)
  if (M$nrow == 1L) m <- M$ncol else if (M$ncol == 1L) m <- M$nrow
  else .pv_err("TypeError", "``self`` must be a row or a column matrix")
  if (X$ncol == 1L) n <- X$nrow else if (X$nrow == 1L) n <- X$ncol
  else .pv_err("TypeError", "X must be a row or a column matrix")
  d <- vector("list", m * n)
  for (j in seq_len(m)) for (i in seq_len(n)) d[[(j - 1L) * n + i]] <- .pv_diff(M$data[[j]], X$data[[i]])
  .new_mat(d, as.integer(m), as.integer(n))
}

# dict.fromkeys(keys, 0) as a pydict
.pm_zero_dict <- function(keys) {
  keys <- if (inherits(keys, "symmatrix")) keys$data else unclass(keys)
  .pydict(lapply(keys, as_sym), rep(list(S_Zero), length(keys)))
}

.pm_dict_pairs <- function(d) Map(list, d$keys, d$values)
.pm_mxreplace <- function(M, d) {
  if (!length(d$keys)) return(M)
  pairs <- .pm_dict_pairs(d)
  .pv_mmap(M, function(z) xreplace(z, pairs))
}
.pm_xreplace <- function(e, d) if (!length(d$keys)) as_sym(e) else xreplace(as_sym(e), .pm_dict_pairs(d))

.pm_dict_get <- function(d, k) {
  i <- match(as_sym(k)$key, vapply(d$keys, .k, ""))
  if (is.na(i)) .pv_err("KeyError", .pv_str(k))
  d$values[[i]]
}

.pm_mdiff_t <- function(M) .pv_mmap(M, function(z) .pv_diff(z, .pv_t()))

.pm_mat_str_list <- function(xs) paste0("[", paste(vapply(xs, .pv_str, ""), collapse = ", "), "]")

# ---------------------------------------------------------------------------
# MethodBase

.pvP_MethodBase_constraints_jacobian <- function(self) .pv_linear_eq_to_matrix(self$velocity_constraints$data, self$u$data)[[1L]]

.pvM_MethodBase_rhs <- function(self, inv_method = NULL, ...) {
  if (is.null(inv_method)) return(.pm_LUsolve(self$mass_matrix_full, self$forcing_full))
  .pv_matmul(inverse(self$mass_matrix_full), self$forcing_full)
}

#' @export
format.MethodBase <- function(x, ...) sprintf("<%s object>", class(x)[1L])

# ---------------------------------------------------------------------------
# KanesMethod

#' Kane's method
#'
#' `KanesMethod()` forms the equations of motion by Kane's method, as
#' SymPy's `KanesMethod`. Call `kanes_equations(bodies, loads)` (reached
#' with `$`) to form Fr and Fr*, then read `mass_matrix`, `forcing`,
#' `mass_matrix_full`, `forcing_full`, `rhs()`, `kindiffdict()` or
#' `linearize()`.
#'
#' @param frame The inertial [ReferenceFrame()].
#' @param q_ind,u_ind Independent generalized coordinates and speeds.
#' @param kd_eqs Kinematic differential equations.
#' @param q_dependent,u_dependent Dependent coordinates and speeds.
#' @param configuration_constraints,velocity_constraints,acceleration_constraints,nonholonomic_constraints
#'   Constraint equations.
#' @param u_auxiliary Auxiliary speeds.
#' @param bodies,forcelist Bodies and loads.
#' @param explicit_kinematics Solve the kinematic equations explicitly.
#' @param kd_eqs_solver,constraint_solver Linear solvers: `"LU"` or a
#'   function of `A` and `b`.
#' @return A `KanesMethod` object.
#' @examples
#' q <- dynamicsymbols("q"); u <- dynamicsymbols("u")
#' N <- ReferenceFrame("N")
#' P <- PhysPoint("P")
#' P$set_vel(N, u * N$x)
#' km <- KanesMethod(N, list(q), list(u), kd_eqs = list(dynamicsymbols("q", 1) - u))
#' invisible(km$kanes_equations(list(Particle("p", P, Symbol("m"))),
#'                              list(Force(P, -Symbol("k") * q * N$x))))
#' km$mass_matrix
#' km$forcing
#' @export
KanesMethod <- function(frame, q_ind, u_ind, kd_eqs = NULL, q_dependent = NULL,
                        configuration_constraints = NULL, u_dependent = NULL,
                        velocity_constraints = NULL, acceleration_constraints = NULL,
                        u_auxiliary = NULL, bodies = NULL, forcelist = NULL,
                        explicit_kinematics = TRUE, kd_eqs_solver = "LU",
                        constraint_solver = "LU", nonholonomic_constraints = NULL) {
  self <- new.env(parent = emptyenv())
  class(self) <- c("KanesMethod", "MethodBase", "physics_obj")
  assign("_uid", .pv_next_uid(), envir = self)
  if (!.pm_truthy(q_ind)) {
    q_ind <- list(dynamicsymbols("dummy_q"))
    kd_eqs <- list(dynamicsymbols("dummy_kd"))
  }
  if (!.pv_is_frame(frame)) .pv_err("TypeError", "An inertial ReferenceFrame must be supplied")
  self[["_inertial"]] <- frame
  self[["_fr"]] <- NULL
  self[["_frstar"]] <- NULL
  self[["_forcelist"]] <- forcelist
  self[["_bodylist"]] <- bodies
  self[["explicit_kinematics"]] <- explicit_kinematics
  self[["_constraint_solver"]] <- constraint_solver
  .pm_km_init_vectors(self, q_ind, q_dependent, u_ind, u_dependent, u_auxiliary)
  .pm_validate_coordinates(self[["_q"]], self[["_u"]])
  .pm_km_init_kindiffeq(self, kd_eqs, kd_eqs_solver)
  .pm_km_init_constraints(self, configuration_constraints, velocity_constraints,
                          acceleration_constraints, nonholonomic_constraints,
                          linear_solver = constraint_solver)
  self
}

.pm_is_iterable <- function(x) !is_sym(x) && !is.numeric(x) && !inherits(x, "physics_obj") && (is.list(x) || inherits(x, "symmatrix"))

.pm_km_init_vectors <- function(self, q_ind, q_dep, u_ind, u_dep, u_aux) {
  t <- .pv_t()
  q_dep <- .pm_none_handler(q_dep)
  if (!.pm_is_iterable(q_ind)) .pv_err("TypeError", "Generalized coordinates must be an iterable.")
  q_ind <- .pm_colmat(q_ind)
  self[["_qdep"]] <- q_dep
  self[["_q"]] <- .pm_vstack(list(q_ind, q_dep))
  self[["_qdot"]] <- .pm_mdiff_t(self[["_q"]])
  u_dep <- .pm_none_handler(u_dep)
  if (!.pm_is_iterable(u_ind)) .pv_err("TypeError", "Generalized speeds must be an iterable.")
  u_ind <- .pm_colmat(u_ind)
  self[["_udep"]] <- u_dep
  self[["_u"]] <- .pm_vstack(list(u_ind, u_dep))
  self[["_udot"]] <- .pm_mdiff_t(self[["_u"]])
  self[["_uaux"]] <- .pm_none_handler(u_aux)
  invisible(self)
}

.pm_km_init_constraints <- function(self, config, vel, acc, nonholonomic, linear_solver = "LU") {
  linear_solver <- .pm_parse_linear_solver(linear_solver)
  t <- .pv_t()
  num_speeds <- .pm_len(self[["_u"]])
  num_dep_speeds <- .pm_len(self[["_udep"]])
  num_dof <- num_speeds - num_dep_speeds
  self[["_f_h"]] <- .pm_none_handler(config)
  if (.pm_len(self[["_qdep"]]) != .pm_len(self[["_f_h"]]))
    .pv_err("ValueError", "There must be an equal number of dependent coordinates and configuration constraints.")
  if (!is.null(nonholonomic) && !is.null(vel))
    .pv_err("ValueError", "Only one of the kwargs velocity_constraints or nonholonomic_constraints can be supplied at a time.")
  qmap <- self[["_qdot_u_map"]]
  self[["_nonholonomic_constraints"]] <- msubs(.pm_none_handler(nonholonomic), qmap)
  vel <- .pm_none_handler(vel)
  hol <- self[["_f_h"]]
  nonhol <- self[["_nonholonomic_constraints"]]
  if (.pm_truthy(vel)) {
    self[["_velocity_constraints"]] <- msubs(vel, qmap)
  } else if (.pm_truthy(hol) && .pm_truthy(nonhol)) {
    con_diff <- msubs(.pm_mdiff_t(hol), qmap)
    self[["_velocity_constraints"]] <- .pm_cjoin(con_diff, nonhol)
  } else if (!.pm_truthy(hol) && .pm_truthy(nonhol)) {
    self[["_velocity_constraints"]] <- nonhol
  } else if (.pm_truthy(hol) && !.pm_truthy(nonhol)) {
    self[["_velocity_constraints"]] <- msubs(.pm_mdiff_t(hol), qmap)
  } else {
    self[["_velocity_constraints"]] <- .pm_empty()
  }
  acc <- .pm_none_handler(acc)
  if (.pm_len(self[["_velocity_constraints"]]) != num_dep_speeds)
    .pv_err("ValueError", "There must be an equal number of dependent speeds and velocity constraints.")
  if (.pm_truthy(acc) && .pm_len(acc) != num_dep_speeds)
    .pv_err("ValueError", "There must be an equal number of dependent speeds and acceleration constraints.")
  if (.pm_truthy(self[["_velocity_constraints"]])) {
    lm <- .pv_linear_eq_to_matrix(self[["_velocity_constraints"]]$data, self[["_u"]]$data)
    self[["_k_nh"]] <- lm[[1L]]
    self[["_f_nh"]] <- .pm_mneg(lm[[2L]])
    if (!.pm_truthy(acc)) {
      f_dnh <- .pv_madd(.pv_matmul(.pm_mdiff_t(self[["_k_nh"]]), self[["_u"]]), .pm_mdiff_t(self[["_f_nh"]]))
      f_dnh <- msubs(f_dnh, qmap)
      self[["_f_dnh"]] <- f_dnh
      self[["_k_dnh"]] <- self[["_k_nh"]]
      self[["_acceleration_constraints"]] <- .pv_madd(.pv_matmul(self[["_k_dnh"]], self[["_udot"]]), self[["_f_dnh"]])
    } else {
      qdd_repl <- .pydict(lapply(qmap$keys, function(qi) .pv_diff(qi, t)),
                          lapply(qmap$values, function(expr) msubs(.pv_diff(expr, t), qmap)))
      acc <- msubs(acc, qdd_repl, qmap)
      self[["_acceleration_constraints"]] <- acc
      lm <- .pv_linear_eq_to_matrix(acc$data, self[["_udot"]]$data)
      self[["_k_dnh"]] <- lm[[1L]]
      self[["_f_dnh"]] <- .pm_mneg(lm[[2L]])
    }
    self[["_B_ind"]] <- .pm_cols(self[["_k_nh"]], 0L, num_dof)
    self[["_B_dep"]] <- .pm_cols(self[["_k_nh"]], num_dof, num_speeds)
    self[["_Ars"]] <- tryCatch(.pm_mneg(linear_solver(self[["_B_dep"]], self[["_B_ind"]])),
                               NonInvertibleMatrixError = function(e)
                                 .pv_err("NonInvertibleMatrixError", "Ars is singular, maybe due to a wrong set of dependent speeds. Try a different set of dependent speeds."))
  } else {
    for (nm in c("_f_nh", "_k_nh", "_f_dnh", "_k_dnh", "_acceleration_constraints", "_B_ind", "_B_dep", "_Ars"))
      self[[nm]] <- .pm_empty()
  }
  invisible(self)
}

.pvP_KanesMethod_constraints_jacobian <- function(self) .pm_rjoin(self[["_B_ind"]], self[["_B_dep"]])

.pm_km_init_kindiffeq <- function(self, kdeqs, linear_solver = "LU") {
  linear_solver <- .pm_parse_linear_solver(linear_solver)
  if (.pm_truthy(kdeqs)) {
    if (.pm_len(self[["_q"]]) != .pm_len(kdeqs))
      .pv_err("ValueError", "There must be an equal number of kinematic differential equations and coordinates.")
    u <- self[["_u"]]
    qdot <- self[["_qdot"]]
    kdeqs <- .pm_colmat(kdeqs)
    u_zero <- .pm_zero_dict(u)
    uaux_zero <- .pm_zero_dict(self[["_uaux"]])
    qdot_zero <- .pm_zero_dict(qdot)
    k_ku <- .pm_jacobian(kdeqs, u)
    k_kqdot <- .pm_jacobian(kdeqs, qdot)
    f_k <- .pm_mxreplace(.pm_mxreplace(kdeqs, u_zero), qdot_zero)
    dy_syms <- find_dynamicsymbols(.pm_rjoin(.pm_rjoin(k_ku, k_kqdot), f_k))
    dk <- vapply(dy_syms, .k, "")
    nonlin_vars <- Filter(function(v) v$key %in% dk, c(u$data, qdot$data))
    if (length(nonlin_vars))
      .pv_err("ValueError", sprintf("The provided kinematic differential equations are nonlinear in %s. They must be linear in the generalized speeds and derivatives of the generalized coordinates.",
                                    .pm_mat_str_list(nonlin_vars)))
    self[["_f_k_implicit"]] <- .pm_mxreplace(f_k, uaux_zero)
    self[["_k_ku_implicit"]] <- .pm_mxreplace(k_ku, uaux_zero)
    self[["_k_kqdot_implicit"]] <- k_kqdot
    f_k_explicit <- linear_solver(k_kqdot, f_k)
    k_ku_explicit <- linear_solver(k_kqdot, k_ku)
    vals <- .pm_mneg(.pv_madd(.pv_matmul(k_ku_explicit, u), f_k_explicit))
    self[["_qdot_u_map"]] <- .pydict(qdot$data, vals$data)
    self[["_f_k"]] <- .pm_mxreplace(f_k_explicit, uaux_zero)
    self[["_k_ku"]] <- .pm_mxreplace(k_ku_explicit, uaux_zero)
    self[["_k_kqdot"]] <- .mech_eye(.pm_len(qdot))
  } else {
    self[["_qdot_u_map"]] <- .pydict(list(), list())
    self[["_f_k_implicit"]] <- self[["_f_k"]] <- .pm_empty()
    self[["_k_ku_implicit"]] <- self[["_k_ku"]] <- .pm_empty()
    self[["_k_kqdot_implicit"]] <- self[["_k_kqdot"]] <- .pm_empty()
  }
  invisible(self)
}

.pm_msubs_vec <- function(v, d) if (length(d$keys)) msubs(v, d) else v

.pm_km_form_fr <- function(self, fl) {
  if (!is.null(fl) && (.pm_len(fl) == 0L || !.pm_is_iterable(fl)))
    .pv_err("ValueError", "Force pairs must be supplied in an non-empty iterable or None.")
  N <- self[["_inertial"]]
  qmap <- self[["_qdot_u_map"]]
  vf <- .pm_f_list_parser(fl, N)
  vel_list <- lapply(vf[[1L]], function(i) msubs(i, qmap))
  f_list <- lapply(vf[[2L]], function(i) msubs(i, qmap))
  o <- .pm_len(self[["_u"]])
  b <- length(f_list)
  FR <- .mech_zeros(o, 1L)
  partials <- .pv_partial_velocity(vel_list, self[["_u"]], N)
  for (i in seq_len(o))
    FR$data[[i]] <- Reduce(function(acc, j) acc + .pv_dot_dispatch(partials[[j]][[i]], f_list[[j]]), seq_len(b), S_Zero)
  if (.pm_truthy(self[["_udep"]])) {
    p <- o - .pm_len(self[["_udep"]])
    FRtilde <- .pm_rows(FR, 0L, p)
    FRold <- .pm_rows(FR, p, o)
    FRtilde <- .pv_madd(FRtilde, .pv_matmul(.pv_T(self[["_Ars"]]), FRold))
    FR <- FRtilde
  }
  self[["_forcelist"]] <- fl
  self[["_fr"]] <- FR
  FR
}

.pm_km_form_frstar <- function(self, bl) {
  if (!.pm_is_iterable(bl)) .pv_err("TypeError", "Bodies must be supplied in an iterable.")
  t <- .pv_t()
  N <- self[["_inertial"]]
  qmap <- self[["_qdot_u_map"]]
  udot_zero <- .pm_zero_dict(self[["_udot"]])
  uaux_zero <- .pm_zero_dict(self[["_uaux"]])
  uauxdot <- lapply(self[["_uaux"]]$data, function(i) .pv_diff(i, t))
  uauxdot_zero <- .pm_zero_dict(uauxdot)
  q_ddot_u_map <- .pydict(c(lapply(qmap$keys, function(k) .pv_diff(k, t)), qmap$keys),
                          c(lapply(qmap$values, function(v) .pm_xreplace(.pv_diff(v, t), qmap)), qmap$values))
  get_partial_velocity <- function(body) {
    if (inherits(body, "RigidBody")) vlist <- list(.pv_point_vel(body$masscenter, N), .pv_ang_vel_in(body$frame, N))
    else if (inherits(body, "Particle")) vlist <- list(.pv_point_vel(body$point, N))
    else .pv_err("TypeError", "The body list may only contain either RigidBody or Particle as list elements.")
    v <- lapply(vlist, function(vel) msubs(vel, qmap))
    .pv_partial_velocity(v, self[["_u"]], N)
  }
  bl <- unclass(bl)
  partials <- lapply(bl, get_partial_velocity)
  o <- .pm_len(self[["_u"]])
  MM <- .mech_zeros(o, o)
  nonMM <- .mech_zeros(o, 1L)
  zero_uaux <- function(expr) msubs(expr, uaux_zero)
  zero_udot_uaux <- function(expr) msubs(msubs(expr, udot_zero), uaux_zero)
  mm_add <- function(j, k, v) MM$data[[(j - 1L) * o + k]] <<- MM$data[[(j - 1L) * o + k]] + v
  for (i in seq_along(bl)) {
    body <- bl[[i]]
    if (inherits(body, "RigidBody")) {
      M <- zero_uaux(body$mass)
      I <- zero_uaux(body$central_inertia)
      vel <- zero_uaux(.pv_point_vel(body$masscenter, N))
      omega <- zero_uaux(.pv_ang_vel_in(body$frame, N))
      acc <- zero_udot_uaux(.pv_point_acc(body$masscenter, N))
      inertial_force <- .pv_diff(M, t) * vel + M * acc
      inertial_torque <- zero_uaux((.pv_ddot(.pv_time_derivative(I, body$frame), omega)) +
                                     msubs(.pv_ddot(I, .pvM_ReferenceFrame_ang_acc_in(body$frame, N)), udot_zero) +
                                     (.pv_vcross(omega, .pv_ddot(I, omega))))
      for (j in seq_len(o)) {
        tmp_vel <- zero_uaux(partials[[i]][[1L]][[j]])
        tmp_ang <- zero_uaux(.pv_ddot(I, partials[[i]][[2L]][[j]]))
        for (k in seq_len(o)) {
          mm_add(j, k, M * .pv_vdot(tmp_vel, partials[[i]][[1L]][[k]]))
          mm_add(j, k, .pv_vdot(tmp_ang, partials[[i]][[2L]][[k]]))
        }
        nonMM$data[[j]] <- nonMM$data[[j]] + .pv_vdot(inertial_force, partials[[i]][[1L]][[j]])
        nonMM$data[[j]] <- nonMM$data[[j]] + .pv_vdot(inertial_torque, partials[[i]][[2L]][[j]])
      }
    } else {
      M <- zero_uaux(body$mass)
      vel <- zero_uaux(.pv_point_vel(body$point, N))
      acc <- zero_udot_uaux(.pv_point_acc(body$point, N))
      inertial_force <- .pv_diff(M, t) * vel + M * acc
      for (j in seq_len(o)) {
        temp <- zero_uaux(partials[[i]][[1L]][[j]])
        for (k in seq_len(o)) mm_add(j, k, M * .pv_vdot(temp, partials[[i]][[1L]][[k]]))
        nonMM$data[[j]] <- nonMM$data[[j]] + .pv_vdot(inertial_force, partials[[i]][[1L]][[j]])
      }
    }
  }
  MM <- zero_uaux(msubs(MM, q_ddot_u_map))
  nonMM <- msubs(msubs(nonMM, q_ddot_u_map), udot_zero, uauxdot_zero, uaux_zero)
  fr_star <- .pm_mneg(.pv_madd(.pv_matmul(MM, msubs(self[["_udot"]], uauxdot_zero)), nonMM))
  if (.pm_truthy(self[["_udep"]])) {
    p <- o - .pm_len(self[["_udep"]])
    fr_star_ind <- .pm_rows(fr_star, 0L, p)
    fr_star_dep <- .pm_rows(fr_star, p, o)
    ArsT <- .pv_T(self[["_Ars"]])
    fr_star <- .pv_madd(fr_star_ind, .pv_matmul(ArsT, fr_star_dep))
    MMi <- .pm_rows(MM, 0L, p)
    MMd <- .pm_rows(MM, p, o)
    MM <- .pv_madd(MMi, .pv_matmul(ArsT, MMd))
    nonMM <- .pv_madd(.pm_rows(nonMM, 0L, p), .pv_matmul(ArsT, .pm_rows(nonMM, p, o)))
  }
  self[["_bodylist"]] <- bl
  self[["_frstar"]] <- fr_star
  self[["_k_d"]] <- MM
  self[["_f_d"]] <- .pm_mneg(.pv_msub(self[["_fr"]], nonMM))
  fr_star
}

.pvM_KanesMethod__form_fr <- function(self, fl) .pm_km_form_fr(self, fl)
.pvM_KanesMethod__form_frstar <- function(self, bl) .pm_km_form_frstar(self, bl)

.pvM_KanesMethod_to_linearizer <- function(self, linear_solver = "LU") {
  if (is.null(self[["_fr"]]) || is.null(self[["_frstar"]])) .pv_err("ValueError", "Need to compute Fr, Fr* first.")
  t <- .pv_t()
  f_c <- self[["_f_h"]]
  f_v <- if (.pm_truthy(self[["_f_nh"]]) && .pm_truthy(self[["_k_nh"]]))
    .pv_madd(self[["_f_nh"]], .pv_matmul(self[["_k_nh"]], self[["_u"]])) else .pm_empty()
  f_a <- if (.pm_truthy(self[["_f_dnh"]]) && .pm_truthy(self[["_k_dnh"]]))
    .pv_madd(self[["_f_dnh"]], .pv_matmul(self[["_k_dnh"]], self[["_udot"]])) else .pm_empty()
  u_zero <- .pm_zero_dict(self[["_u"]])
  ud_zero <- .pm_zero_dict(self[["_udot"]])
  qd_zero <- .pm_zero_dict(self[["_qdot"]])
  qd_u_zero <- .pm_zero_dict(.pm_vstack(list(self[["_qdot"]], self[["_u"]])))
  f_0 <- .pv_madd(msubs(self[["_f_k"]], u_zero), .pv_matmul(self[["_k_kqdot"]], self[["_qdot"]]))
  f_1 <- .pv_madd(msubs(self[["_f_k"]], qd_zero), .pv_matmul(self[["_k_ku"]], self[["_u"]]))
  f_2 <- msubs(self[["_frstar"]], qd_u_zero)
  f_3 <- .pv_madd(msubs(self[["_frstar"]], ud_zero), self[["_fr"]])
  f_4 <- .mech_zeros(.pm_len(f_2), 1L)
  q <- self[["_q"]]
  u <- self[["_u"]]
  q_i <- if (.pm_truthy(self[["_qdep"]])) .pm_rows(q, 0L, q$nrow - .pm_len(self[["_qdep"]])) else q
  q_d <- self[["_qdep"]]
  u_i <- if (.pm_truthy(self[["_udep"]])) .pm_rows(u, 0L, u$nrow - .pm_len(self[["_udep"]])) else u
  u_d <- self[["_udep"]]
  uaux <- self[["_uaux"]]
  uauxdot <- .pm_mdiff_t(uaux)
  uaux_zero <- .pm_zero_dict(.pm_vstack(list(uaux, uauxdot)))
  sym_list <- .pm_vstack(list(q, self[["_qdot"]], u, self[["_udot"]], uaux, uauxdot))$data
  for (M in list(self[["_k_kqdot"]], self[["_k_ku"]], self[["_f_k"]], self[["_k_dnh"]], self[["_f_dnh"]], self[["_k_d"]]))
    if (length(find_dynamicsymbols(M, sym_list)))
      .pv_err("ValueError", paste0("Cannot have dynamicsymbols outside dynamic", strrep(" ", 30L), "forcing vector."))
  r <- unclass(find_dynamicsymbols(msubs(self[["_f_d"]], uaux_zero), sym_list))
  r <- .pm_sort_default(r)
  rk <- vapply(r, .k, "")
  for (i in r) if (.pv_diff(i, t)$key %in% rk)
    .pv_err("ValueError", paste0("Cannot have derivatives of specified", strrep(" ", 34L), "quantities when linearizing forcing terms."))
  Linearizer(f_0, f_1, f_2, f_3, f_4, f_c, f_v, f_a, q, u, q_i, q_d, u_i, u_d, r, linear_solver = linear_solver)
}

# list.sort(key=default_sort_key)
.pm_sort_default <- function(xs) {
  if (length(xs) < 2L) return(xs)
  .key_sort(lapply(xs, .sort_key), xs)
}

.pvM_KanesMethod_linearize <- function(self, ..., new_method = NULL, linear_solver = "LU") {
  linearizer <- .pvM_KanesMethod_to_linearizer(self, linear_solver = linear_solver)
  result <- linearizer$linearize(...)
  .pytuple(c(unclass(result), list(linearizer$r)))
}

.pvM_KanesMethod_kanes_equations <- function(self, bodies = NULL, loads = NULL) {
  if (is.null(bodies)) bodies <- self[["_bodylist"]]
  if (is.null(loads) && !is.null(self[["_forcelist"]])) loads <- self[["_forcelist"]]
  if (is.list(loads) && !inherits(loads, "pytuple") && !length(loads)) loads <- NULL
  if (!.pm_truthy(self[["_k_kqdot"]]))
    .pv_err("AttributeError", "Create an instance of KanesMethod with kinematic differential equations to use this method.")
  fr <- .pm_km_form_fr(self, loads)
  frstar <- .pm_km_form_frstar(self, bodies)
  if (.pm_truthy(self[["_uaux"]])) {
    if (!.pm_truthy(self[["_udep"]])) {
      km <- KanesMethod(self[["_inertial"]], self[["_q"]], self[["_uaux"]], u_auxiliary = self[["_uaux"]],
                        constraint_solver = self[["_constraint_solver"]])
    } else {
      km <- KanesMethod(self[["_inertial"]], self[["_q"]], self[["_uaux"]], u_auxiliary = self[["_uaux"]],
                        u_dependent = self[["_udep"]],
                        velocity_constraints = .pv_madd(.pv_matmul(self[["_k_nh"]], self[["_u"]]), self[["_f_nh"]]),
                        acceleration_constraints = .pv_madd(.pv_matmul(self[["_k_dnh"]], self[["_udot"]]), self[["_f_dnh"]]),
                        constraint_solver = self[["_constraint_solver"]])
    }
    km[["_qdot_u_map"]] <- self[["_qdot_u_map"]]
    self[["_km"]] <- km
    fraux <- .pm_km_form_fr(km, loads)
    frstaraux <- .pm_km_form_frstar(km, bodies)
    self[["_aux_eq"]] <- .pv_madd(fraux, frstaraux)
    self[["_fr"]] <- .pm_cjoin(fr, fraux)
    self[["_frstar"]] <- .pm_cjoin(frstar, frstaraux)
  }
  .pytuple(list(self[["_fr"]], self[["_frstar"]]))
}

.pvM_KanesMethod__form_eoms <- function(self) {
  r <- .pvM_KanesMethod_kanes_equations(self, self$bodies, self$loads)
  .pv_madd(r[[1L]], r[[2L]])
}

.pvM_KanesMethod_rhs <- function(self, inv_method = NULL) {
  nq <- .pm_len(self[["_q"]])
  rhs <- .mech_zeros(nq + .pm_len(self[["_u"]]), 1L)
  kdes <- .pvM_KanesMethod_kindiffdict(self)
  for (i in seq_len(nq)) rhs$data[[i]] <- .pm_dict_get(kdes, .pv_diff(self[["_q"]]$data[[i]], .pv_t()))
  sol <- if (is.null(inv_method)) .pm_LUsolve(self$mass_matrix, self$forcing)
         else .pv_matmul(inverse(self$mass_matrix), self$forcing)
  for (i in seq_along(sol$data)) rhs$data[[nq + i]] <- sol$data[[i]]
  rhs
}

.pvM_KanesMethod_kindiffdict <- function(self) {
  if (!length(self[["_qdot_u_map"]]$keys))
    .pv_err("AttributeError", "Create an instance of KanesMethod with kinematic differential equations to use this method.")
  self[["_qdot_u_map"]]
}

.pm_need_fr <- function(self) {
  if (!.pm_truthy(self[["_fr"]]) || !.pm_truthy(self[["_frstar"]])) .pv_err("ValueError", "Need to compute Fr, Fr* first.")
}

.pvP_KanesMethod_auxiliary_eqs <- function(self) {
  .pm_need_fr(self)
  if (!.pm_truthy(self[["_uaux"]])) .pv_err("ValueError", "No auxiliary speeds have been declared.")
  self[["_aux_eq"]]
}

.pvP_KanesMethod_mass_matrix_kin <- function(self)
  if (isTRUE(self[["explicit_kinematics"]])) self[["_k_kqdot"]] else self[["_k_kqdot_implicit"]]

.pvP_KanesMethod_forcing_kin <- function(self) {
  if (isTRUE(self[["explicit_kinematics"]]))
    .pm_mneg(.pv_madd(.pv_matmul(self[["_k_ku"]], self[["_u"]]), self[["_f_k"]]))
  else .pm_mneg(.pv_madd(.pv_matmul(self[["_k_ku_implicit"]], self[["_u"]]), self[["_f_k_implicit"]]))
}

.pvP_KanesMethod_mass_matrix <- function(self) {
  .pm_need_fr(self)
  .pm_vstack(list(self[["_k_d"]], self[["_k_dnh"]]))
}

.pvP_KanesMethod_forcing <- function(self) {
  .pm_need_fr(self)
  .pm_mneg(.pm_vstack(list(self[["_f_d"]], self[["_f_dnh"]])))
}

.pvP_KanesMethod_mass_matrix_full <- function(self) {
  .pm_need_fr(self)
  o <- .pm_len(self[["_u"]])
  n <- .pm_len(self[["_q"]])
  .pm_cjoin(.pm_rjoin(self$mass_matrix_kin, .mech_zeros(n, o)), .pm_rjoin(.mech_zeros(o, n), self$mass_matrix))
}

.pvP_KanesMethod_forcing_full <- function(self) .pm_vstack(list(self$forcing_kin, self$forcing))

.pvP_KanesMethod_q <- function(self) self[["_q"]]
.pvP_KanesMethod_u <- function(self) self[["_u"]]
.pvP_KanesMethod_bodylist <- function(self) {
  warning("'bodylist' is deprecated, use 'bodies' instead.", call. = FALSE)
  self[["_bodylist"]]
}
.pvP_KanesMethod_forcelist <- function(self) {
  warning("'forcelist' is deprecated, use 'loads' instead.", call. = FALSE)
  self[["_forcelist"]]
}
.pvP_KanesMethod_bodies <- function(self) self[["_bodylist"]]
.pvP_KanesMethod_loads <- function(self) self[["_forcelist"]]
.pvP_KanesMethod_frame <- function(self) self[["_inertial"]]
.pvP_KanesMethod_holonomic_constraints <- function(self) self[["_f_h"]]
.pvP_KanesMethod_nonholonomic_constraints <- function(self) self[["_nonholonomic_constraints"]]
.pvP_KanesMethod_velocity_constraints <- function(self) self[["_velocity_constraints"]]
.pvP_KanesMethod_acceleration_constraints <- function(self) self[["_acceleration_constraints"]]
