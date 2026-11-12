# Lagrange's method, after sympy/physics/mechanics/lagrange.py.

#' Lagrange's method
#'
#' `LagrangesMethod()` forms the equations of motion from a Lagrangian, as
#' SymPy's `LagrangesMethod`. Call `form_lagranges_equations()` (reached
#' with `$`), then read `mass_matrix`, `forcing`, `mass_matrix_full`,
#' `forcing_full`, `rhs()`, `solve_multipliers()` or `linearize()`.
#'
#' @param Lagrangian The Lagrangian, a scalar expression.
#' @param qs Generalized coordinates.
#' @param forcelist Loads.
#' @param bodies Bodies.
#' @param frame The inertial [ReferenceFrame()], needed with loads.
#' @param hol_coneqs,nonhol_coneqs Holonomic and nonholonomic constraints.
#' @return A `LagrangesMethod` object.
#' @examples
#' q <- dynamicsymbols("q")
#' qd <- dynamicsymbols("q", 1)
#' m <- Symbol("m"); k <- Symbol("k")
#' lm <- LagrangesMethod(m * qd^2 / 2 - k * q^2 / 2, list(q))
#' lm$form_lagranges_equations()
#' @export
LagrangesMethod <- function(Lagrangian, qs, forcelist = NULL, bodies = NULL, frame = NULL,
                            hol_coneqs = NULL, nonhol_coneqs = NULL) {
  self <- new.env(parent = emptyenv())
  class(self) <- c("LagrangesMethod", "MethodBase", "physics_obj")
  assign("_uid", .pv_next_uid(), envir = self)
  t <- .pv_t()
  self[["_L"]] <- .pv_colv(list(as_sym(Lagrangian)))
  self[["eom"]] <- NULL
  for (nm in c("_m_cd", "_m_d", "_f_cd", "_f_d", "lam_coeffs")) self[[nm]] <- .pm_empty()
  if (!is.null(forcelist)) {
    if (!.pm_is_iterable(forcelist)) .pv_err("TypeError", "Force pairs must be supplied in an iterable.")
    self[["_forcelist"]] <- unclass(forcelist)
  } else self[["_forcelist"]] <- list()
  if (is.null(bodies)) self[["_bodies"]] <- list()
  else {
    if (!.pm_is_iterable(bodies)) .pv_err("TypeError", "bodies should be an iterable.")
    self[["_bodies"]] <- unclass(bodies)
  }
  if (!is.null(frame) && !.pv_is_frame(frame)) .pv_err("TypeError", "frame must be a valid ReferenceFrame")
  self[["inertial"]] <- frame
  for (nm in c("lam_vec", "_term1", "_term2", "_term3", "_term4")) self[[nm]] <- .pm_empty()
  if (!.pm_is_iterable(qs)) .pv_err("TypeError", "Generalized coordinates must be an iterable")
  self[["_q"]] <- .pm_colmat(qs)
  self[["_qdots"]] <- .pm_mdiff_t(self[["_q"]])
  self[["_qdoubledots"]] <- .pm_mdiff_t(self[["_qdots"]])
  .pm_validate_coordinates(self[["_q"]])
  hol_coneqs <- .pm_none_handler(hol_coneqs)
  nonhol_coneqs <- .pm_none_handler(nonhol_coneqs)
  self[["coneqs"]] <- .pm_vstack(list(.pm_mdiff_t(hol_coneqs), nonhol_coneqs))
  self[["_hol_coneqs"]] <- hol_coneqs
  self[["_nonhol_coneqs"]] <- nonhol_coneqs
  self
}

.pvM_LagrangesMethod_form_lagranges_equations <- function(self) {
  t <- .pv_t()
  qds <- self[["_qdots"]]
  qdd_zero <- .pm_zero_dict(self[["_qdoubledots"]])
  n <- .pm_len(self[["_q"]])
  self[["_term1"]] <- .pv_T(.pm_mdiff_t(.pm_jacobian(self[["_L"]], qds)))
  self[["_term2"]] <- .pv_T(.pm_jacobian(self[["_L"]], self[["_q"]]))
  if (.pm_truthy(self[["coneqs"]])) {
    coneqs <- self[["coneqs"]]
    m <- .pm_len(coneqs)
    self[["lam_vec"]] <- .pm_colmat(dynamicsymbols(paste0("lam1:", m + 1L)))
    self[["lam_coeffs"]] <- .pv_linear_eq_to_matrix(.pm_mneg(coneqs)$data, qds$data)[[1L]]
    self[["_term3"]] <- .pv_matmul(.pv_T(self[["lam_coeffs"]]), self[["lam_vec"]])
    self[["_m_cd"]] <- .pm_mneg(self[["lam_coeffs"]])
    diffconeqs <- .pm_mdiff_t(coneqs)
    self[["_diffconeqs"]] <- diffconeqs
    self[["_f_cd"]] <- .pm_mneg(.pv_msubs(diffconeqs, qdd_zero))
  } else self[["_term3"]] <- .mech_zeros(n, 1L)
  loads <- self[["_forcelist"]]
  if (.pm_truthy(loads)) {
    N <- self[["inertial"]]
    term4 <- .mech_zeros(n, 1L)
    for (i in seq_len(n)) {
      qd <- qds$data[[i]]
      vf <- .pm_f_list_parser(loads, N)
      term4$data[[i]] <- Reduce(function(acc, j) acc + .pv_dot_dispatch(.pv_vdiff(vf[[1L]][[j]], qd, N), vf[[2L]][[j]]),
                                seq_along(vf[[1L]]), S_Zero)
    }
    self[["_term4"]] <- term4
  } else self[["_term4"]] <- .mech_zeros(n, 1L)
  without_lam <- .pv_msub(.pv_msub(self[["_term1"]], self[["_term2"]]), self[["_term4"]])
  lm <- .pv_linear_eq_to_matrix(without_lam$data, self[["_qdoubledots"]]$data)
  self[["_m_d"]] <- lm[[1L]]
  self[["_f_d"]] <- lm[[2L]]
  self[["eom"]] <- .pv_msub(without_lam, self[["_term3"]])
  self[["eom"]]
}

.pvM_LagrangesMethod__form_eoms <- function(self) .pvM_LagrangesMethod_form_lagranges_equations(self)

.pm_need_eom <- function(self) if (is.null(self[["eom"]])) .pv_err("ValueError", "Need to compute the equations of motion first")

.pvP_LagrangesMethod_mass_matrix <- function(self) {
  .pm_need_eom(self)
  if (.pm_truthy(self[["coneqs"]])) .pm_rjoin(self[["_m_d"]], .pv_T(self[["lam_coeffs"]])) else self[["_m_d"]]
}

.pvP_LagrangesMethod_mass_matrix_full <- function(self) {
  .pm_need_eom(self)
  n <- .pm_len(self[["_q"]])
  m <- .pm_len(self[["coneqs"]])
  row1 <- .pm_rjoin(.mech_eye(n), .mech_zeros(n, n + m))
  row2 <- .pm_rjoin(.mech_zeros(n, n), self$mass_matrix)
  if (.pm_truthy(self[["coneqs"]])) {
    row3 <- .pm_rjoin(.pm_rjoin(.mech_zeros(m, n), self[["_m_cd"]]), .mech_zeros(m, m))
    return(.pm_cjoin(.pm_cjoin(row1, row2), row3))
  }
  .pm_cjoin(row1, row2)
}

.pvP_LagrangesMethod_forcing <- function(self) {
  .pm_need_eom(self)
  self[["_f_d"]]
}

.pvP_LagrangesMethod_forcing_full <- function(self) {
  .pm_need_eom(self)
  if (.pm_truthy(self[["coneqs"]])) return(.pm_cjoin(.pm_cjoin(self[["_qdots"]], self$forcing), self[["_f_cd"]]))
  .pm_cjoin(self[["_qdots"]], self$forcing)
}

.pvM_LagrangesMethod_to_linearizer <- function(self, q_ind = NULL, qd_ind = NULL, q_dep = NULL, qd_dep = NULL,
                                               linear_solver = "LU") {
  t <- .pv_t()
  q <- self[["_q"]]
  u <- self[["_qdots"]]
  ud <- .pm_mdiff_t(u)
  lams <- self[["lam_vec"]]
  q_i <- .pm_none_handler(q_ind)
  q_d <- .pm_none_handler(q_dep)
  u_i <- .pm_none_handler(qd_ind)
  u_d <- .pm_none_handler(qd_dep)
  f_c <- self[["_hol_coneqs"]]
  f_v <- self[["coneqs"]]
  f_a <- .pm_mdiff_t(f_v)
  f_0 <- u
  f_1 <- .pm_mneg(u)
  f_2 <- self[["_term1"]]
  f_3 <- .pm_mneg(.pv_madd(self[["_term2"]], self[["_term4"]]))
  f_4 <- .pm_mneg(self[["_term3"]])
  if (.pm_len(q_d) != .pm_len(f_c) || .pm_len(u_d) != .pm_len(f_v))
    .pv_err("ValueError", sprintf("Must supply %d dependent coordinates, and %d dependent speeds", .pm_len(f_c), .pm_len(f_v)))
  keys <- function(M) unique(vapply(M$data, .k, ""))
  if (!setequal(keys(.pm_vstack(list(q_i, q_d))), keys(q)))
    .pv_err("ValueError", "Must partition q into q_ind and q_dep, with no extra or missing symbols.")
  if (!setequal(keys(.pm_vstack(list(u_i, u_d))), keys(u)))
    .pv_err("ValueError", "Must partition qd into qd_ind and qd_dep, with no extra or missing symbols.")
  insyms <- .pm_vstack(list(q, u, ud, lams))$data
  r <- .pm_sort_default(unclass(find_dynamicsymbols(f_3, insyms)))
  rk <- vapply(r, .k, "")
  for (i in r) if (.pv_diff(i, t)$key %in% rk)
    .pv_err("ValueError", paste0("Cannot have derivatives of specified", strrep(" ", 34L), "quantities when linearizing forcing terms."))
  Linearizer(f_0, f_1, f_2, f_3, f_4, f_c, f_v, f_a, q, u, q_i, q_d, u_i, u_d, r, lams, linear_solver = linear_solver)
}

.pvM_LagrangesMethod_linearize <- function(self, q_ind = NULL, qd_ind = NULL, q_dep = NULL, qd_dep = NULL,
                                           linear_solver = "LU", ...) {
  linearizer <- .pvM_LagrangesMethod_to_linearizer(self, q_ind, qd_ind, q_dep, qd_dep, linear_solver = linear_solver)
  result <- linearizer$linearize(...)
  .pytuple(c(unclass(result), list(linearizer$r)))
}

.pvM_LagrangesMethod_solve_multipliers <- function(self, op_point = NULL, sol_type = "dict") {
  k <- .pm_len(self[["lam_vec"]])
  if (k == 0L) .pv_err("ValueError", "System has no lagrange multipliers to solve for.")
  if (!is.null(op_point) && !is.list(op_point))
    .pv_err("TypeError", "op_point must be either a dictionary or an iterable of dictionaries.")
  opd <- .pm_op_point_dict(op_point)
  mass_matrix <- .pm_cjoin(self$mass_matrix, .pm_rjoin(.pm_mneg(self[["lam_coeffs"]]), .mech_zeros(k, k)))
  force_matrix <- .pm_cjoin(self$forcing, self[["_f_cd"]])
  mass_matrix <- msubs(mass_matrix, opd)
  force_matrix <- msubs(force_matrix, opd)
  sol <- .pm_LUsolve(mass_matrix, .pm_mneg(force_matrix))$data
  sol_list <- sol[(length(sol) - k + 1L):length(sol)]
  if (identical(sol_type, "dict")) return(.pydict(self[["lam_vec"]]$data, sol_list))
  if (identical(sol_type, "Matrix")) return(.pm_colmat(sol_list))
  .pv_err("ValueError", sprintf("Unknown sol_type %s.", sol_type))
}

.pvP_LagrangesMethod_q <- function(self) self[["_q"]]
.pvP_LagrangesMethod_u <- function(self) self[["_qdots"]]
.pvP_LagrangesMethod_velocity_constraints <- function(self) self[["coneqs"]]
.pvP_LagrangesMethod_acceleration_constraints <- function(self) {
  if (exists("_diffconeqs", envir = self, inherits = FALSE)) return(self[["_diffconeqs"]])
  .pm_mdiff_t(self[["coneqs"]])
}
.pvP_LagrangesMethod_constraints_jacobian <- function(self) .pm_mneg(self[["lam_coeffs"]])
.pvP_LagrangesMethod_bodies <- function(self) self[["_bodies"]]
.pvP_LagrangesMethod_forcelist <- function(self) {
  warning("'forcelist' is deprecated, use 'loads' instead.", call. = FALSE)
  self[["_forcelist"]]
}
.pvP_LagrangesMethod_loads <- function(self) self[["_forcelist"]]
.pvP_LagrangesMethod_frame <- function(self) self[["inertial"]]
.pvP_LagrangesMethod_holonomic_constraints <- function(self) self[["_hol_coneqs"]]
.pvP_LagrangesMethod_nonholonomic_constraints <- function(self) self[["_nonhol_coneqs"]]
