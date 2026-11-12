# Musculotendon models, after sympy/physics/biomechanics/musculotendon.py.

.BM_FORMULATIONS <- c("RIGID_TENDON", "FIBER_LENGTH_EXPLICIT", "TENDON_FORCE_EXPLICIT",
                      "FIBER_LENGTH_IMPLICIT", "TENDON_FORCE_IMPLICIT")

#' Musculotendon formulations
#'
#' The integer enumeration SymPy calls `MusculotendonFormulation`: `0` is
#' `RIGID_TENDON`, `1` `FIBER_LENGTH_EXPLICIT`, `2` `TENDON_FORCE_EXPLICIT`,
#' `3` `FIBER_LENGTH_IMPLICIT` and `4` `TENDON_FORCE_IMPLICIT`. A member is
#' selected by its value or by its name.
#'
#' @param value An integer from 0 to 4, or a member name.
#' @param x A `MusculotendonFormulation` value.
#' @param ... Ignored.
#' @return A `MusculotendonFormulation` value, an integer that prints as
#'   its value.
#' @examples
#' MusculotendonFormulation("RIGID_TENDON")
#' MusculotendonFormulation(2) == 2
#' @export
MusculotendonFormulation <- function(value) {
  if (inherits(value, "MusculotendonFormulation")) return(value)
  if (is.character(value) && length(value) == 1L) {
    v <- match(value, .BM_FORMULATIONS) - 1L
    if (is.na(v)) .pv_err("KeyError", sprintf("'%s'", value))
  } else {
    ok <- (is.numeric(value) || is_sym(value)) && length(value) == 1L
    v <- if (ok) suppressWarnings(as.numeric(if (is_sym(value)) .neval(value, list()) else value)) else NA
    if (is.na(v) || !(v %in% 0:4)) .pv_err("ValueError", sprintf("%s is not a valid MusculotendonFormulation", .pv_objrepr(value)))
    v <- as.integer(v)
  }
  structure(v, class = "MusculotendonFormulation")
}

#' @rdname MusculotendonFormulation
#' @export
format.MusculotendonFormulation <- function(x, ...) as.character(unclass(x))

#' @rdname MusculotendonFormulation
#' @export
print.MusculotendonFormulation <- function(x, ...) {
  cat(sprintf("<MusculotendonFormulation.%s: %d>\n", .BM_FORMULATIONS[unclass(x) + 1L], unclass(x)))
  invisible(x)
}

.bm_formulation_value <- function(x) {
  if (inherits(x, "MusculotendonFormulation")) return(as.integer(unclass(x)))
  if (is.numeric(x) && length(x) == 1L && x %in% 0:4) return(as.integer(x))
  if (is_sym(x) && x$type == "Integer") {
    v <- as.integer(.neval(x, list()))
    if (v %in% 0:4) return(v)
  }
  NA_integer_
}

.bm_symbols <- function(spec) lapply(.pv_expand_names(spec), Symbol)

#' Musculotendon actuators
#'
#' Musculotendon models after SymPy's `MusculotendonDeGroote2016`: a force
#' actuator along a pathway whose force follows from the activation
#' dynamics, the characteristic curves of De Groote et al. (2016) and the
#' chosen [MusculotendonFormulation()]. Constants left as `NULL` become
#' symbols named after the model. The rigid tendon formulation and the
#' explicit fiber length and tendon force formulations are available; the
#' implicit ones raise `NotImplementedError` as in SymPy.
#' `MusculotendonDeGroote2016_with_defaults()` uses the published curve
#' constants and default fiber velocity, pennation and damping values.
#' `MusculotendonBase()` is the abstract base class and raises an error.
#'
#' The objects expose SymPy's attributes through `$`, among them `force`,
#' `pathway`, `activation_dynamics`, `state_vars` (`x`), `input_vars`
#' (`r`), `constants` (`p`), `M`, `F`, `rhs()` and `to_loads()`.
#'
#' @param name Name of the musculotendon.
#' @param pathway A pathway, such as [LinearPathway()].
#' @param activation_dynamics An activation model, such as
#'   [FirstOrderActivationDeGroote2016()].
#' @param musculotendon_dynamics A [MusculotendonFormulation()] value.
#' @param tendon_slack_length,peak_isometric_force,optimal_fiber_length,maximal_fiber_velocity,optimal_pennation_angle,fiber_damping_coefficient
#'   Model constants.
#' @param with_defaults Use the published curve constants.
#' @param ... Ignored.
#' @return A musculotendon object.
#' @examples
#' N <- ReferenceFrame("N")
#' O <- PhysPoint("O"); I <- PhysPoint("I")
#' I$set_pos(O, dynamicsymbols("q") * N$x)
#' muscle <- MusculotendonDeGroote2016("m", LinearPathway(O, I), FirstOrderActivationDeGroote2016("m"))
#' muscle$x
#' muscle$constants
#' @name musculotendon
NULL

#' @rdname musculotendon
#' @export
MusculotendonBase <- function(...) .bm_abstract("MusculotendonBase", "curves")

#' @rdname musculotendon
#' @export
MusculotendonDeGroote2016 <- function(name, pathway, activation_dynamics, musculotendon_dynamics = 0L,
                                      tendon_slack_length = NULL, peak_isometric_force = NULL,
                                      optimal_fiber_length = NULL, maximal_fiber_velocity = NULL,
                                      optimal_pennation_angle = NULL, fiber_damping_coefficient = NULL,
                                      with_defaults = FALSE) {
  self <- .pm_new_obj(c("MusculotendonDeGroote2016", "MusculotendonBase", "ForceActuator", "ActuatorBase"))
  .bm_musculotendon_init(self, name, pathway, activation_dynamics, musculotendon_dynamics, tendon_slack_length,
                         peak_isometric_force, optimal_fiber_length, maximal_fiber_velocity,
                         optimal_pennation_angle, fiber_damping_coefficient, with_defaults)
  self
}

#' @rdname musculotendon
#' @export
MusculotendonDeGroote2016_with_defaults <- function(name, pathway, activation_dynamics, musculotendon_dynamics = 0L,
                                                    tendon_slack_length = NULL, peak_isometric_force = NULL,
                                                    optimal_fiber_length = NULL,
                                                    maximal_fiber_velocity = Float("10.0"),
                                                    optimal_pennation_angle = Float("0.0"),
                                                    fiber_damping_coefficient = Float("0.1"))
  MusculotendonDeGroote2016(name, pathway, activation_dynamics = activation_dynamics,
                            musculotendon_dynamics = musculotendon_dynamics, tendon_slack_length = tendon_slack_length,
                            peak_isometric_force = peak_isometric_force, optimal_fiber_length = optimal_fiber_length,
                            maximal_fiber_velocity = maximal_fiber_velocity,
                            optimal_pennation_angle = optimal_pennation_angle,
                            fiber_damping_coefficient = fiber_damping_coefficient, with_defaults = TRUE)

.pvP_MusculotendonDeGroote2016_curves <- function(self) .bm_degroote_curves()

.bm_musculotendon_init <- function(self, name, pathway, activation_dynamics, musculotendon_dynamics,
                                   tendon_slack_length, peak_isometric_force, optimal_fiber_length,
                                   maximal_fiber_velocity, optimal_pennation_angle, fiber_damping_coefficient,
                                   with_defaults) {
  .bm_set_name(self, name)
  .pvS_ForceActuator_force(self, Symbol("F"))
  .pvS_ForceActuator_pathway(self, pathway)
  if (!inherits(activation_dynamics, "ActivationBase"))
    .pv_err("TypeError", sprintf("Can't set attribute `activation_dynamics` to %s as it must be of type `ActivationBase`, not %s.",
                                 .pv_str(activation_dynamics), .pv_pyclass(activation_dynamics)))
  self[["_activation_dynamics"]] <- activation_dynamics
  self[["_child_objects"]] <- list(activation_dynamics)
  nm <- self$name
  self[["_l_T_slack"]] <- if (!is.null(tendon_slack_length)) tendon_slack_length else Symbol(paste0("l_T_slack_", nm))
  self[["_F_M_max"]] <- if (!is.null(peak_isometric_force)) peak_isometric_force else Symbol(paste0("F_M_max_", nm))
  self[["_l_M_opt"]] <- if (!is.null(optimal_fiber_length)) optimal_fiber_length else Symbol(paste0("l_M_opt_", nm))
  self[["_v_M_max"]] <- if (!is.null(maximal_fiber_velocity)) maximal_fiber_velocity else Symbol(paste0("v_M_max_", nm))
  self[["_alpha_opt"]] <- if (!is.null(optimal_pennation_angle)) optimal_pennation_angle else Symbol(paste0("alpha_opt_", nm))
  self[["_beta"]] <- if (!is.null(fiber_damping_coefficient)) fiber_damping_coefficient else Symbol(paste0("beta_", nm))
  self[["_with_defaults"]] <- isTRUE(with_defaults)
  v <- .bm_formulation_value(musculotendon_dynamics)
  if (is.na(v))
    .pv_err("TypeError", sprintf("Musculotendon dynamics %s passed to `musculotendon_dynamics` was of type %s, must be <enum 'MusculotendonFormulation'>.",
                                 .pv_objrepr(musculotendon_dynamics), .pv_pyclass(musculotendon_dynamics)))
  switch(v + 1L,
         .bm_rigid_tendon(self),
         .bm_fiber_length_explicit(self),
         .bm_tendon_force_explicit(self),
         .pv_err("NotImplementedError", ""),
         .pv_err("NotImplementedError", ""))
  self[["_musculotendon_dynamics"]] <- musculotendon_dynamics
  self[["_force"]] <- -self[["_F_T"]]
  invisible(self)
}

# the constants used in the expressions, as sym
.bm_k <- function(self, field) as_sym(self[[field]])

.bm_curve_call <- function(self, curve, x, spec) {
  if (self[["_with_defaults"]]) return(list(f = .bm_with_defaults(curve, x), k = list()))
  k <- .bm_symbols(sprintf(spec, self$name))
  list(f = do.call(curve, c(list(x), k)), k = k)
}

.bm_rigid_tendon <- function(self) {
  curves <- self$curves
  pw <- self[["_pathway"]]
  l_MT <- pw$length
  v_MT <- pw$extension_velocity
  l_T_slack <- .bm_k(self, "_l_T_slack"); l_M_opt <- .bm_k(self, "_l_M_opt"); alpha_opt <- .bm_k(self, "_alpha_opt")
  self[["_l_MT"]] <- l_MT
  self[["_v_MT"]] <- v_MT
  self[["_l_T"]] <- l_T_slack
  self[["_l_T_tilde"]] <- S_One
  l_M <- sqrt((l_MT - l_T_slack)^2 + (l_M_opt * sin(alpha_opt))^2)
  self[["_l_M"]] <- l_M
  l_M_tilde <- l_M / l_M_opt
  self[["_l_M_tilde"]] <- l_M_tilde
  v_M <- v_MT * (l_MT - l_T_slack) / l_M
  self[["_v_M"]] <- v_M
  v_M_tilde <- v_M / .bm_k(self, "_v_M_max")
  self[["_v_M_tilde"]] <- v_M_tilde
  fl_T <- .bm_curve_call(self, curves$tendon_force_length, S_One, "c_0:4_fl_T_%s")
  fl_M_pas <- .bm_curve_call(self, curves$fiber_force_length_passive, l_M_tilde, "c_0:2_fl_M_pas_%s")
  fl_M_act <- .bm_curve_call(self, curves$fiber_force_length_active, l_M_tilde, "c_0:12_fl_M_act_%s")
  fv_M <- .bm_curve_call(self, curves$fiber_force_velocity, v_M_tilde, "c_0:4_fv_M_%s")
  self[["_fl_T"]] <- fl_T$f; self[["_fl_M_pas"]] <- fl_M_pas$f
  self[["_fl_M_act"]] <- fl_M_act$f; self[["_fv_M"]] <- fv_M$f
  F_M_tilde <- self$a * fl_M_act$f * fv_M$f + fl_M_pas$f + .bm_k(self, "_beta") * v_M_tilde
  self[["_F_M_tilde"]] <- F_M_tilde
  self[["_F_T_tilde"]] <- F_M_tilde
  F_M <- F_M_tilde * .bm_k(self, "_F_M_max")
  self[["_F_M"]] <- F_M
  self[["_cos_alpha"]] <- cos(alpha_opt)
  self[["_F_T"]] <- F_M * self[["_cos_alpha"]]
  self[["_state_vars"]] <- zeros(0L, 1L)
  self[["_input_vars"]] <- zeros(0L, 1L)
  self[["_state_eqns"]] <- zeros(0L, 1L)
  self[["_curve_constants"]] <- .bm_colmat(c(fl_T$k, fl_M_pas$k, fl_M_act$k, fv_M$k))
}

.bm_fiber_length_explicit <- function(self) {
  curves <- self$curves
  pw <- self[["_pathway"]]
  l_M_tilde <- dynamicsymbols(paste0("l_M_tilde_", self$name))
  self[["_l_M_tilde"]] <- l_M_tilde
  l_MT <- pw$length
  v_MT <- pw$extension_velocity
  self[["_l_MT"]] <- l_MT
  self[["_v_MT"]] <- v_MT
  l_M_opt <- .bm_k(self, "_l_M_opt")
  l_M <- l_M_tilde * l_M_opt
  self[["_l_M"]] <- l_M
  l_T <- l_MT - sqrt(l_M^2 - (l_M_opt * sin(.bm_k(self, "_alpha_opt")))^2)
  self[["_l_T"]] <- l_T
  l_T_tilde <- l_T / .bm_k(self, "_l_T_slack")
  self[["_l_T_tilde"]] <- l_T_tilde
  cos_alpha <- (l_MT - l_T) / l_M
  self[["_cos_alpha"]] <- cos_alpha
  fl_T <- .bm_curve_call(self, curves$tendon_force_length, l_T_tilde, "c_0:4_fl_T_%s")
  fl_M_pas <- .bm_curve_call(self, curves$fiber_force_length_passive, l_M_tilde, "c_0:2_fl_M_pas_%s")
  fl_M_act <- .bm_curve_call(self, curves$fiber_force_length_active, l_M_tilde, "c_0:12_fl_M_act_%s")
  self[["_fl_T"]] <- fl_T$f; self[["_fl_M_pas"]] <- fl_M_pas$f; self[["_fl_M_act"]] <- fl_M_act$f
  F_T_tilde <- fl_T$f
  self[["_F_T_tilde"]] <- F_T_tilde
  F_T <- F_T_tilde * .bm_k(self, "_F_M_max")
  self[["_F_T"]] <- F_T
  F_M <- F_T / cos_alpha
  self[["_F_M"]] <- F_M
  F_M_tilde <- F_M / .bm_k(self, "_F_M_max")
  self[["_F_M_tilde"]] <- F_M_tilde
  fv_M <- (F_M_tilde - fl_M_pas$f) / (self$a * fl_M_act$f)
  self[["_fv_M"]] <- fv_M
  v_M_tilde <- .bm_curve_call(self, curves$fiber_force_velocity_inverse, fv_M, "c_0:4_fv_M_%s")
  self[["_v_M_tilde"]] <- v_M_tilde$f
  dl_M_tilde_dt <- .bm_k(self, "_v_M_max") / l_M_opt * v_M_tilde$f
  self[["_dl_M_tilde_dt"]] <- dl_M_tilde_dt
  self[["_state_vars"]] <- .pv_colv(list(l_M_tilde))
  self[["_input_vars"]] <- zeros(0L, 1L)
  self[["_state_eqns"]] <- .pv_colv(list(dl_M_tilde_dt))
  self[["_curve_constants"]] <- .bm_colmat(c(fl_T$k, fl_M_pas$k, fl_M_act$k, v_M_tilde$k))
}

.bm_tendon_force_explicit <- function(self) {
  curves <- self$curves
  pw <- self[["_pathway"]]
  t <- .pv_t()
  F_T_tilde <- dynamicsymbols(paste0("F_T_tilde_", self$name))
  self[["_F_T_tilde"]] <- F_T_tilde
  l_MT <- pw$length
  v_MT <- pw$extension_velocity
  self[["_l_MT"]] <- l_MT
  self[["_v_MT"]] <- v_MT
  self[["_fl_T"]] <- F_T_tilde
  fl_T_inv <- .bm_curve_call(self, curves$tendon_force_length_inverse, F_T_tilde, "c_0:4_fl_T_%s")
  self[["_fl_T_inv"]] <- fl_T_inv$f
  l_T_tilde <- fl_T_inv$f
  self[["_l_T_tilde"]] <- l_T_tilde
  l_T <- l_T_tilde * .bm_k(self, "_l_T_slack")
  self[["_l_T"]] <- l_T
  l_M_opt <- .bm_k(self, "_l_M_opt")
  l_M <- sqrt((l_MT - l_T)^2 + (l_M_opt * sin(.bm_k(self, "_alpha_opt")))^2)
  self[["_l_M"]] <- l_M
  l_M_tilde <- l_M / l_M_opt
  self[["_l_M_tilde"]] <- l_M_tilde
  fl_M_pas <- .bm_curve_call(self, curves$fiber_force_length_passive, l_M_tilde, "c_0:2_fl_M_pas_%s")
  fl_M_act <- .bm_curve_call(self, curves$fiber_force_length_active, l_M_tilde, "c_0:12_fl_M_act_%s")
  self[["_fl_M_pas"]] <- fl_M_pas$f; self[["_fl_M_act"]] <- fl_M_act$f
  cos_alpha <- (l_MT - l_T) / l_M
  self[["_cos_alpha"]] <- cos_alpha
  F_T <- F_T_tilde * .bm_k(self, "_F_M_max")
  self[["_F_T"]] <- F_T
  F_M <- F_T / cos_alpha
  self[["_F_M"]] <- F_M
  F_M_tilde <- F_M / .bm_k(self, "_F_M_max")
  self[["_F_M_tilde"]] <- F_M_tilde
  fv_M <- (F_M_tilde - fl_M_pas$f) / (self$a * fl_M_act$f)
  self[["_fv_M"]] <- fv_M
  fv_M_inv <- .bm_curve_call(self, curves$fiber_force_velocity_inverse, fv_M, "c_0:4_fv_M_%s")
  self[["_fv_M_inv"]] <- fv_M_inv$f
  v_M_tilde <- fv_M_inv$f
  self[["_v_M_tilde"]] <- v_M_tilde
  v_M <- v_M_tilde * .bm_k(self, "_v_M_max")
  self[["_v_M"]] <- v_M
  v_T <- v_MT - v_M / cos_alpha
  self[["_v_T"]] <- v_T
  v_T_tilde <- v_T / .bm_k(self, "_l_T_slack")
  self[["_v_T_tilde"]] <- v_T_tilde
  fl_T <- if (self[["_with_defaults"]]) .bm_with_defaults(curves$tendon_force_length, l_T_tilde)
          else do.call(curves$tendon_force_length, c(list(l_T_tilde), fl_T_inv$k))
  self[["_fl_T"]] <- fl_T
  dF <- .pv_subs(.pv_diff(fl_T, t), .pydict(list(.pv_diff(l_T_tilde, t)), list(v_T_tilde)))
  self[["_dF_T_tilde_dt"]] <- dF
  self[["_state_vars"]] <- .pv_colv(list(F_T_tilde))
  self[["_input_vars"]] <- zeros(0L, 1L)
  self[["_state_eqns"]] <- .pv_colv(list(dF))
  self[["_curve_constants"]] <- .bm_colmat(c(fl_T_inv$k, fl_M_pas$k, fl_M_act$k, fv_M_inv$k))
}

.pvP_MusculotendonBase_tendon_slack_length <- function(self) self[["_l_T_slack"]]
.pvP_MusculotendonBase_l_T_slack <- function(self) self[["_l_T_slack"]]
.pvP_MusculotendonBase_peak_isometric_force <- function(self) self[["_F_M_max"]]
.pvP_MusculotendonBase_F_M_max <- function(self) self[["_F_M_max"]]
.pvP_MusculotendonBase_optimal_fiber_length <- function(self) self[["_l_M_opt"]]
.pvP_MusculotendonBase_l_M_opt <- function(self) self[["_l_M_opt"]]
.pvP_MusculotendonBase_maximal_fiber_velocity <- function(self) self[["_v_M_max"]]
.pvP_MusculotendonBase_v_M_max <- function(self) self[["_v_M_max"]]
.pvP_MusculotendonBase_optimal_pennation_angle <- function(self) self[["_alpha_opt"]]
.pvP_MusculotendonBase_alpha_opt <- function(self) self[["_alpha_opt"]]
.pvP_MusculotendonBase_fiber_damping_coefficient <- function(self) self[["_beta"]]
.pvP_MusculotendonBase_beta <- function(self) self[["_beta"]]
.pvP_MusculotendonBase_activation_dynamics <- function(self) self[["_activation_dynamics"]]
.pvP_MusculotendonBase_excitation <- function(self) self[["_activation_dynamics"]][["_e"]]
.pvP_MusculotendonBase_e <- function(self) self[["_activation_dynamics"]][["_e"]]
.pvP_MusculotendonBase_activation <- function(self) self[["_activation_dynamics"]][["_a"]]
.pvP_MusculotendonBase_a <- function(self) self[["_activation_dynamics"]][["_a"]]
.pvP_MusculotendonBase_musculotendon_dynamics <- function(self) self[["_musculotendon_dynamics"]]

.bm_children <- function(self, attr) lapply(self[["_child_objects"]], function(ch) .pv_getattr(ch, attr))

.pvP_MusculotendonBase_state_vars <- function(self) .bm_vstack(c(list(self[["_state_vars"]]), .bm_children(self, "state_vars")))
.pvP_MusculotendonBase_x <- .pvP_MusculotendonBase_state_vars
.pvP_MusculotendonBase_input_vars <- function(self) .bm_vstack(c(list(self[["_input_vars"]]), .bm_children(self, "input_vars")))
.pvP_MusculotendonBase_r <- .pvP_MusculotendonBase_input_vars

.pvP_MusculotendonBase_constants <- function(self) {
  cs <- list(self[["_l_T_slack"]], self[["_F_M_max"]], self[["_l_M_opt"]], self[["_v_M_max"]],
             self[["_alpha_opt"]], self[["_beta"]])
  .bm_vstack(c(list(.bm_symbolic_constants(cs)), .bm_children(self, "constants"), list(self[["_curve_constants"]])))
}
.pvP_MusculotendonBase_p <- .pvP_MusculotendonBase_constants

.pvP_MusculotendonBase_M <- function(self)
  .bm_block_diag(c(list(eye(self[["_state_vars"]]$nrow)), .bm_children(self, "M")))
.pvP_MusculotendonBase_F <- function(self) .bm_vstack(c(list(self[["_state_eqns"]]), .bm_children(self, "F")))

.pvM_MusculotendonBase_rhs <- function(self) {
  v <- .bm_formulation_value(self[["_musculotendon_dynamics"]])
  child_rhs <- lapply(self[["_child_objects"]], function(ch) ch$rhs())
  if (v == 0L) return(.bm_vstack(child_rhs))
  if (v %in% c(1L, 2L)) return(.bm_vstack(c(list(self[["_state_eqns"]]), child_rhs)))
  solve(self$M, self$F)
}

# Matrix.vstack, keeping the column count of empty blocks
.bm_vstack <- function(mats) {
  if (!length(mats)) return(zeros(0L, 0L))
  ncol <- mats[[1L]]$ncol
  for (M in mats) if (M$nrow > 0L || M$ncol > 0L) { ncol <- M$ncol; break }
  for (M in mats)
    if (M$ncol != ncol && M$nrow > 0L)
      .pv_err("ShapeError", sprintf("The matrices have incompatible number of columns (%d and %d)", ncol, M$ncol))
  data <- unlist(lapply(mats, function(M) M$data), recursive = FALSE)
  .new_mat(if (is.null(data)) list() else data, sum(vapply(mats, function(M) M$nrow, 0L)), ncol)
}

# sympy.matrices.dense.diag of matrix blocks
.bm_block_diag <- function(mats) {
  nr <- sum(vapply(mats, function(M) M$nrow, 0L)); nc <- sum(vapply(mats, function(M) M$ncol, 0L))
  R <- zeros(nr, nc)
  r0 <- 0L; c0 <- 0L
  for (M in mats) {
    for (i in seq_len(M$nrow)) for (j in seq_len(M$ncol))
      R$data[[(r0 + i - 1L) * nc + c0 + j]] <- M$data[[(i - 1L) * M$ncol + j]]
    r0 <- r0 + M$nrow; c0 <- c0 + M$ncol
  }
  R
}

.pvR_MusculotendonBase <- function(x)
  sprintf(paste0("%s(%s, pathway=%s, activation_dynamics=%s, musculotendon_dynamics=%s, tendon_slack_length=%s, ",
                 "peak_isometric_force=%s, optimal_fiber_length=%s, maximal_fiber_velocity=%s, ",
                 "optimal_pennation_angle=%s, fiber_damping_coefficient=%s)"),
          class(x)[1L], .pv_objrepr(x$name), .pv_objrepr(x$pathway), .pv_objrepr(x$activation_dynamics),
          .pv_str(x$musculotendon_dynamics), .pv_objrepr(x[["_l_T_slack"]]), .pv_objrepr(x[["_F_M_max"]]),
          .pv_objrepr(x[["_l_M_opt"]]), .pv_objrepr(x[["_v_M_max"]]), .pv_objrepr(x[["_alpha_opt"]]),
          .pv_objrepr(x[["_beta"]]))

#' @export
format.MusculotendonBase <- function(x, ...) .pv_str(x$force)
