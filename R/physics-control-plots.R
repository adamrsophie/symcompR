# Numerical data and plots of SISO systems, after
# sympy/physics/control/control_plots.py. The plots use base R graphics.

.cp_check_system <- function(system) {
  if (!inherits(system, "SISOLinearTimeInvariant"))
    .physics_error("NotImplementedError", "Only SISO LTI systems are currently supported.")
  e <- system$to_expr()
  if (length(free_symbols(e)) > 1L)
    .physics_error("ValueError", "Extra degree of freedom found. Make sure that there are no free symbols in the dynamical system other than the variable of Laplace transform.")
  if (.lti_has_fn(e, "exp")) .physics_error("NotImplementedError", "Time delay terms are not supported.")
  invisible(NULL)
}

# a number as an R double, or complex when its imaginary part is not zero
.cp_value <- function(z) {
  v <- as.complex(N(z))
  if (Im(v) == 0) Re(v) else v
}

# _poly_roots: exact roots over ZZ and QQ, numerical roots otherwise
.cp_poly_roots <- function(f, x) {
  cs <- all_coeffs(Poly(f, x))
  if (length(cs) <= 1L) return(list())
  r <- if (.rc_domain(cs) %in% c("ZZ", "QQ")) all_roots(f, x) else nroots(f, x)
  lapply(r, .cp_value)
}

# the complex coefficients of a polynomial, highest degree first
.cp_coeffs <- function(f, x) vapply(all_coeffs(Poly(expand(f), x)), function(c) as.complex(N(c)), 0i)

.cp_horner <- function(cf, z) {
  out <- rep(0i, length(z))
  for (c in cf) out <- out * z + c
  out
}

# the inverse Laplace transform of num/den at times t, from the partial
# fractions of the strictly proper part over the poles
.cp_inverse_laplace <- function(num, den, t) {
  num <- num / den[[1L]]; den <- den / den[[1L]]
  n <- length(den) - 1L
  if (length(num) > n) {
    rem <- num
    while (length(rem) > n) {
      q <- rem[[1L]]
      rem <- rem[-1L] - q * c(den[-1L], rep(0, length(rem) - length(den)))
    }
    num <- rem
  }
  if (n == 0L) return(rep(0, length(t)))
  poles <- polyroot(rev(den))
  groups <- list()
  for (p in poles) {
    hit <- FALSE
    for (g in seq_along(groups)) {
      if (Mod(groups[[g]]$a - p) < 1e-5 * max(1, Mod(p))) {
        groups[[g]]$m <- groups[[g]]$m + 1L
        groups[[g]]$sum <- groups[[g]]$sum + p
        hit <- TRUE
        break
      }
    }
    if (!hit) groups[[length(groups) + 1L]] <- list(a = p, m = 1L, sum = p)
  }
  groups <- lapply(groups, function(g) list(a = g$sum / g$m, m = g$m))
  terms <- do.call(rbind, lapply(groups, function(g) cbind(a = rep(g$a, g$m), k = seq_len(g$m))))
  rad <- 1 + 2 * max(Mod(poles))
  z <- rad * exp(2i * pi * (seq_len(2L * n) - 0.5) / (2L * n))
  Fz <- .cp_horner(num, z) / .cp_horner(den, z)
  A <- sapply(seq_len(nrow(terms)), function(c) 1 / (z - terms[c, "a"])^Re(terms[c, "k"]))
  cf <- qr.solve(A, Fz)
  y <- rep(0i, length(t))
  for (c in seq_len(nrow(terms))) {
    k <- Re(terms[c, "k"])
    y <- y + cf[[c]] * t^(k - 1) * exp(terms[c, "a"] * t) / factorial(k - 1)
  }
  Re(y)
}

.cp_response <- function(system, extra_power, scale, lower_limit, upper_limit, n) {
  .cp_check_system(system)
  tf <- system$doit()
  num <- .cp_coeffs(tf$num, tf$var) * as.complex(N(as_sym(scale)))
  den <- c(.cp_coeffs(tf$den, tf$var), rep(0, extra_power))
  x <- seq(as.numeric(lower_limit), as.numeric(upper_limit), length.out = n)
  list(x, .cp_inverse_laplace(num, den, x))
}

.cp_lower <- function(lower_limit)
  if (isTRUE(as.numeric(lower_limit) < 0)) .physics_error("ValueError", "Lower limit of time must be greater than or equal to zero.")

#' Numerical data and plots of SISO systems
#'
#' `pole_zero_numerical_data(system)` gives the zeros and the poles of a
#' SISO system. `impulse_response_numerical_data()`,
#' `step_response_numerical_data()` and `ramp_response_numerical_data()`
#' sample the impulse, unit step and ramp responses on `n` equally spaced
#' times; `bode_magnitude_numerical_data()` and
#' `bode_phase_numerical_data()` sample the magnitude (dB) and the phase on
#' `n` logarithmically spaced frequencies. `nyquist_plot_expr()` and
#' `nichols_plot_expr()` give the expressions of the Nyquist and Nichols
#' plots. The `*_plot()` functions draw the plots with base R graphics.
#' After `sympy.physics.control.control_plots`.
#'
#' @param system A SISO system.
#' @param prec The precision of the evaluation (kept for compatibility).
#' @param lower_limit,upper_limit The time range.
#' @param slope The slope of the ramp.
#' @param n The number of points.
#' @param adaptive Ignored; the points are always equally spaced.
#' @param initial_exp,final_exp The frequency range as powers of ten.
#' @param freq_unit `"rad/sec"` or `"Hz"`.
#' @param phase_unit `"rad"` or `"deg"`.
#' @param phase_unwrap Whether to unwrap the phase.
#' @param initial_omega,final_omega The frequency range of the Nyquist and
#'   Nichols plots.
#' @param color,pole_color,zero_color The colors.
#' @param pole_markersize,zero_markersize The marker sizes.
#' @param grid,show_axes,show Whether to draw a grid, the axes, and the plot.
#' @param ... Further arguments.
#' @return A list of the x and y values (`pole_zero_numerical_data()`: of the
#'   zeros and poles), a list of expressions, or, for the plots, the data
#'   invisibly.
#' @examples
#' s <- Symbol("s")
#' G <- TransferFunction(1, s^2 + s + 1, s)
#' pole_zero_numerical_data(G)
#' step_response_numerical_data(G, n = 5)
#' @export
pole_zero_numerical_data <- function(system) {
  .cp_check_system(system)
  tf <- system$doit()
  .pytuple(list(.cp_poly_roots(tf$num, tf$var), .cp_poly_roots(tf$den, tf$var)))
}

#' @rdname pole_zero_numerical_data
#' @export
impulse_response_numerical_data <- function(system, prec = 8L, lower_limit = 0L, upper_limit = 10L, n = 1000L, adaptive = FALSE, ...) {
  .cp_lower(lower_limit)
  .cp_response(system, 0L, 1L, lower_limit, upper_limit, n)
}

#' @rdname pole_zero_numerical_data
#' @export
step_response_numerical_data <- function(system, prec = 8L, lower_limit = 0L, upper_limit = 10L, n = 1000L, adaptive = FALSE, ...) {
  .cp_lower(lower_limit)
  .cp_response(system, 1L, 1L, lower_limit, upper_limit, n)
}

#' @rdname pole_zero_numerical_data
#' @export
ramp_response_numerical_data <- function(system, slope = 1L, prec = 8L, lower_limit = 0L, upper_limit = 10L, n = 1000L, adaptive = FALSE, ...) {
  if (isTRUE(as.numeric(slope) < 0)) .physics_error("ValueError", "Slope must be greater than or equal to zero.")
  .cp_lower(lower_limit)
  .cp_response(system, 2L, slope, lower_limit, upper_limit, n)
}

.cp_freq <- function(system, initial_exp, final_exp, freq_unit, n) {
  .cp_check_system(system)
  if (!freq_unit %in% c("rad/sec", "Hz")) .physics_error("ValueError", "Only \"rad/sec\" and \"Hz\" are accepted frequency units.")
  e <- system$to_expr()
  nd <- .as_numer_denom(e)
  x <- 10^seq(as.numeric(initial_exp), as.numeric(final_exp), length.out = n)
  z <- if (freq_unit == "Hz") 2i * pi * x else 1i * x
  list(x = x, H = .cp_horner(.cp_coeffs(nd[[1L]], system$var), z) / .cp_horner(.cp_coeffs(nd[[2L]], system$var), z))
}

#' @rdname pole_zero_numerical_data
#' @export
bode_magnitude_numerical_data <- function(system, initial_exp = -5L, final_exp = 5L, freq_unit = "rad/sec", n = 1000L, ...) {
  f <- .cp_freq(system, initial_exp, final_exp, freq_unit, n)
  list(f$x, 20 * log10(Mod(f$H)))
}

#' @rdname pole_zero_numerical_data
#' @export
bode_phase_numerical_data <- function(system, initial_exp = -5L, final_exp = 5L, freq_unit = "rad/sec",
                                      phase_unit = "rad", phase_unwrap = TRUE, n = 1000L, ...) {
  .cp_check_system(system)
  if (!freq_unit %in% c("rad/sec", "Hz")) .physics_error("ValueError", "Only \"rad/sec\" and \"Hz\" are accepted frequency units.")
  if (!phase_unit %in% c("rad", "deg")) .physics_error("ValueError", "Only \"rad\" and \"deg\" are accepted phase units.")
  f <- .cp_freq(system, initial_exp, final_exp, freq_unit, n)
  y <- Arg(f$H)
  if (phase_unit == "deg") y <- y * 180 / pi
  if (isTRUE(phase_unwrap)) {
    half <- if (phase_unit == "rad") pi else 180
    unit <- 2 * half
    if (length(y) > 1L) for (i in 2:length(y)) {
      d <- y[[i]] - y[[i - 1L]]
      if (d > half) y[[i]] <- y[[i]] - unit else if (d < -half) y[[i]] <- y[[i]] + unit
    }
  }
  list(f$x, y)
}

#' @rdname pole_zero_numerical_data
#' @export
nyquist_plot_expr <- function(system) {
  w <- Dummy("w", real = TRUE)
  w_expr <- subs(system$to_expr(), system$var, S_I * w)
  ri <- .as_real_imag(w_expr)
  .pytuple(list(ri[[1L]], ri[[2L]], w))
}

#' @rdname pole_zero_numerical_data
#' @export
nichols_plot_expr <- function(system) {
  w <- Dummy("w", real = TRUE)
  H <- subs(system$to_expr(), system$var, S_I * w)
  .pytuple(list(20L * log(Abs(H), 10L), arg(H) * 180L / S_Pi, w))
}

.cp_frame <- function(xlab, ylab, main, grid, show_axes, ...) {
  if (isTRUE(grid)) graphics::grid()
  if (isTRUE(show_axes)) graphics::abline(h = 0, v = 0, col = "black")
  graphics::title(main = main, xlab = xlab, ylab = ylab)
}

.cp_name <- function(system) .lti_str(system)

#' @rdname pole_zero_numerical_data
#' @export
pole_zero_plot <- function(system, pole_color = "blue", pole_markersize = 10, zero_color = "orange",
                           zero_markersize = 7, grid = TRUE, show_axes = TRUE, show = TRUE, ...) {
  d <- pole_zero_numerical_data(system)
  z <- as.complex(unlist(d[[1L]])); p <- as.complex(unlist(d[[2L]]))
  allv <- c(z, p, 0)
  graphics::plot(Re(allv), Im(allv), type = "n", xlab = "", ylab = "")
  graphics::points(Re(p), Im(p), pch = 4, cex = pole_markersize / 7, col = pole_color)
  graphics::points(Re(z), Im(z), pch = 1, cex = zero_markersize / 7, col = zero_color)
  .cp_frame("Real Axis", "Imaginary Axis", paste("Poles and Zeros of", .cp_name(system)), grid, show_axes)
  invisible(d)
}

.cp_line_plot <- function(d, color, xlab, ylab, main, grid, show_axes, log = "") {
  graphics::plot(d[[1L]], d[[2L]], type = "l", col = color, xlab = "", ylab = "", log = log)
  .cp_frame(xlab, ylab, main, grid, show_axes)
  invisible(d)
}

#' @rdname pole_zero_numerical_data
#' @export
step_response_plot <- function(system, color = "blue", prec = 8L, lower_limit = 0L, upper_limit = 10L,
                               show_axes = FALSE, grid = TRUE, show = TRUE, ...) {
  d <- step_response_numerical_data(system, prec = prec, lower_limit = lower_limit, upper_limit = upper_limit, ...)
  .cp_line_plot(d, color, "Time (s)", "Amplitude", paste("Unit Step Response of", .cp_name(system)), grid, show_axes)
}

#' @rdname pole_zero_numerical_data
#' @export
impulse_response_plot <- function(system, color = "blue", prec = 8L, lower_limit = 0L, upper_limit = 10L,
                                  show_axes = FALSE, grid = TRUE, show = TRUE, ...) {
  d <- impulse_response_numerical_data(system, prec = prec, lower_limit = lower_limit, upper_limit = upper_limit, ...)
  .cp_line_plot(d, color, "Time (s)", "Amplitude", paste("Impulse Response of", .cp_name(system)), grid, show_axes)
}

#' @rdname pole_zero_numerical_data
#' @export
ramp_response_plot <- function(system, slope = 1L, color = "blue", prec = 8L, lower_limit = 0L, upper_limit = 10L,
                               show_axes = FALSE, grid = TRUE, show = TRUE, ...) {
  d <- ramp_response_numerical_data(system, slope = slope, prec = prec, lower_limit = lower_limit, upper_limit = upper_limit, ...)
  .cp_line_plot(d, color, "Time (s)", "Amplitude",
                sprintf("Ramp Response of %s [Slope = %s]", .cp_name(system), format(slope)), grid, show_axes)
}

#' @rdname pole_zero_numerical_data
#' @export
bode_magnitude_plot <- function(system, initial_exp = -5L, final_exp = 5L, color = "blue", show_axes = FALSE,
                                grid = TRUE, show = TRUE, freq_unit = "rad/sec", ...) {
  d <- bode_magnitude_numerical_data(system, initial_exp = initial_exp, final_exp = final_exp, freq_unit = freq_unit)
  .cp_line_plot(d, color, sprintf("Frequency (%s) [Log Scale]", freq_unit), "Magnitude (dB)",
                paste("Bode Plot (Magnitude) of", .cp_name(system)), grid, show_axes, log = "x")
}

#' @rdname pole_zero_numerical_data
#' @export
bode_phase_plot <- function(system, initial_exp = -5L, final_exp = 5L, color = "blue", show_axes = FALSE,
                            grid = TRUE, show = TRUE, freq_unit = "rad/sec", phase_unit = "rad", phase_unwrap = TRUE, ...) {
  d <- bode_phase_numerical_data(system, initial_exp = initial_exp, final_exp = final_exp, freq_unit = freq_unit,
                                 phase_unit = phase_unit, phase_unwrap = phase_unwrap)
  .cp_line_plot(d, color, sprintf("Frequency (%s) [Log Scale]", freq_unit), sprintf("Phase (%s)", phase_unit),
                paste("Bode Plot (Phase) of", .cp_name(system)), grid, show_axes, log = "x")
}

#' @rdname pole_zero_numerical_data
#' @export
bode_plot <- function(system, initial_exp = -5L, final_exp = 5L, grid = TRUE, show_axes = FALSE, show = TRUE,
                      freq_unit = "rad/sec", phase_unit = "rad", phase_unwrap = TRUE, ...) {
  mag <- bode_magnitude_numerical_data(system, initial_exp = initial_exp, final_exp = final_exp, freq_unit = freq_unit)
  ph <- bode_phase_numerical_data(system, initial_exp = initial_exp, final_exp = final_exp, freq_unit = freq_unit,
                                  phase_unit = phase_unit, phase_unwrap = phase_unwrap)
  op <- graphics::par(mfrow = c(2L, 1L))
  on.exit(graphics::par(op))
  .cp_line_plot(mag, "blue", "", "Magnitude (dB)", paste("Bode Plot of", .cp_name(system)), grid, show_axes, log = "x")
  .cp_line_plot(ph, "blue", sprintf("Frequency (%s) [Log Scale]", freq_unit), sprintf("Phase (%s)", phase_unit),
                "", grid, show_axes, log = "x")
  invisible(list(mag, ph))
}

.cp_parametric <- function(fx, fy, w, initial_omega, final_omega, n = 1000L) {
  ws <- seq(as.numeric(initial_omega), as.numeric(final_omega), length.out = n)
  ev <- function(e) vapply(ws, function(v) as.numeric(N(subs(e, w, v))), 0)
  list(ev(fx), ev(fy))
}

#' @rdname pole_zero_numerical_data
#' @export
nyquist_plot <- function(system, initial_omega = 0.01, final_omega = 100, show = TRUE, color = "blue", n = 1000L, ...) {
  .cp_check_system(system)
  e <- nyquist_plot_expr(system)
  d <- .cp_parametric(e[[1L]], e[[2L]], e[[3L]], initial_omega, final_omega, n)
  graphics::plot(d[[1L]], d[[2L]], type = "l", col = color, xlab = "Real Axis", ylab = "Imaginary Axis",
                 ylim = range(c(d[[2L]], -d[[2L]])), main = paste("Nyquist Plot of", .cp_name(system)))
  graphics::lines(d[[1L]], -d[[2L]], col = color)
  invisible(d)
}

#' @rdname pole_zero_numerical_data
#' @export
nichols_plot <- function(system, initial_omega = 0.01, final_omega = 100, show = TRUE, color = "blue", n = 1000L, ...) {
  .cp_check_system(system)
  e <- nichols_plot_expr(system)
  d <- .cp_parametric(e[[2L]], e[[1L]], e[[3L]], initial_omega, final_omega, n)
  graphics::plot(d[[1L]], d[[2L]], type = "l", col = color, xlab = "Phase [deg]", ylab = "Magnitude [dB]",
                 main = paste("Nichols Plot of", .cp_name(system)))
  invisible(d)
}
