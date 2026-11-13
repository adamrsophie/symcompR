# Quantum mechanical model systems, after sympy/physics/hydrogen.py,
# sho.py, qho_1d.py and pring.py, and the reduced Planck constant of
# sympy/physics/quantum/constants.py that the last two use.

# HBar: a positive real NumberSymbol that prints as hbar.
.HBAR <- new_sym("HBar", name = "hbar", key = "HBar", extra = list(class_key = list(2, 0, "HBar")))
.physics_hbar <- function() .HBAR
register_str_printer("HBar", function(e) "hbar")
register_latex_printer("HBar", function(e) "\\hbar")
.facts_hooks$HBar <- function(e) c(commutative = TRUE, real = TRUE, positive = TRUE, negative = FALSE,
                                   nonzero = TRUE, zero = FALSE, finite = TRUE, irrational = TRUE,
                                   rational = FALSE, integer = FALSE, complex = TRUE, extended_real = TRUE,
                                   nonnegative = TRUE, nonpositive = FALSE, extended_positive = TRUE,
                                   extended_negative = FALSE, extended_nonzero = TRUE, infinite = FALSE,
                                   imaginary = FALSE, hermitian = TRUE)
register_neval("HBar", function(e, env) 1.05457162e-34)

# ---------------------------------------------------------------------------
# Hydrogen atom

#' The hydrogen atom
#'
#' Wavefunctions and energies of the nonrelativistic and relativistic
#' hydrogen atom in atomic units, after `sympy.physics.hydrogen`.
#' `hydrogen_R_nl()` is the radial wavefunction, `hydrogen_Psi_nlm()` the
#' full wavefunction, `hydrogen_E_nl()` the nonrelativistic energy and
#' `hydrogen_E_nl_dirac()` the relativistic energy from the Dirac equation.
#' `hydrogen_Y_lm()` and `hydrogen_Z_lm()` are the expanded complex and real
#' spherical harmonics.
#'
#' @param n The principal quantum number.
#' @param l The angular momentum quantum number.
#' @param m The magnetic quantum number.
#' @param r The radial coordinate.
#' @param phi,theta The azimuthal and polar angles.
#' @param Z The atomic number (1 for hydrogen).
#' @param spin_up `TRUE` for the electron spin up state.
#' @param c The speed of light in atomic units.
#' @return A `sym` expression.
#' @examples
#' r <- Symbol("r")
#' hydrogen_R_nl(1, 0, r)
#' hydrogen_E_nl(2)
#' hydrogen_E_nl_dirac(2, 0, Z = 1, c = 1)
#' @export
hydrogen_R_nl <- function(n, l, r, Z = 1) {
  n <- as_sym(n); l <- as_sym(l); r <- as_sym(r); Z <- as_sym(Z)
  n_r <- n - l - 1
  a <- 1 / Z
  r0 <- 2 * r / (n * a)
  C <- sqrt((S(2) / (n * a))^3 * fcall("factorial", list(n_r)) / (2 * n * fcall("factorial", list(n + l))))
  C * r0^l * expand(fcall("assoc_laguerre", list(n_r, 2 * l + 1, r0))) * exp(-r0 / 2)
}

#' @rdname hydrogen_R_nl
#' @export
hydrogen_Y_lm <- function(l, m, phi, theta) {
  expand(fcall("Ynm", lapply(list(l, m, theta, phi), as_sym)), func = TRUE)
}

#' @rdname hydrogen_R_nl
#' @export
hydrogen_Z_lm <- function(l, m, phi, theta) {
  expand(Znm(l, m, theta, phi), func = TRUE)
}

#' @rdname hydrogen_R_nl
#' @export
hydrogen_Psi_nlm <- function(n, l, m, r, phi, theta, Z = 1) {
  n <- as_sym(n); l <- as_sym(l); m <- as_sym(m); r <- as_sym(r)
  phi <- as_sym(phi); theta <- as_sym(theta); Z <- as_sym(Z)
  if (isTRUE(ask_fact(n, "integer")) && isTRUE(n < 1))
    .physics_error("ValueError", "'n' must be positive integer")
  if (isTRUE(ask_fact(l, "integer")) && !isTRUE(n > l))
    .physics_error("ValueError", "'n' must be greater than 'l'")
  if (isTRUE(ask_fact(m, "integer")) && !isTRUE(abs(m) <= l))
    .physics_error("ValueError", "|'m'| must be less or equal 'l'")
  hydrogen_R_nl(n, l, r, Z) * expand(fcall("Ynm", list(l, m, theta, phi)), func = TRUE)
}

#' @rdname hydrogen_R_nl
#' @export
hydrogen_E_nl <- function(n, Z = 1) {
  n <- as_sym(n); Z <- as_sym(Z)
  if (isTRUE(ask_fact(n, "integer")) && isTRUE(n < 1))
    .physics_error("ValueError", "'n' must be positive integer")
  -Z^2 / (2 * n^2)
}

#' @rdname hydrogen_R_nl
#' @export
hydrogen_E_nl_dirac <- function(n, l, spin_up = TRUE, Z = 1, c = Float("137.035999037")) {
  n <- as_sym(n); l <- as_sym(l); Z <- as_sym(Z); c <- as_sym(c)
  if (!isTRUE(l >= 0)) .physics_error("ValueError", "'l' must be positive or zero")
  if (!isTRUE(n > l)) .physics_error("ValueError", "'n' must be greater than 'l'")
  if (identical(l$key, S_Zero$key) && isFALSE(spin_up)) .physics_error("ValueError", "Spin must be up for l==0.")
  skappa <- if (isTRUE(spin_up)) -l - 1 else -l
  beta <- sqrt(skappa^2 - Z^2 / c^2)
  c^2 / sqrt(1 + Z^2 / (n + skappa + beta)^2 / c^2) - c^2
}

# A SymPy exception class as an R condition.
.physics_error <- function(cls, msg) {
  stop(structure(class = c(cls, "error", "condition"), list(message = msg, call = NULL)))
}

# ---------------------------------------------------------------------------
# Three dimensional isotropic harmonic oscillator

#' The isotropic harmonic oscillator
#'
#' `sho_R_nl()` is the radial wavefunction of the three dimensional
#' isotropic harmonic oscillator and `sho_E_nl()` its energy, after
#' `sympy.physics.sho`.
#'
#' @param n The nodal quantum number.
#' @param l The angular momentum quantum number.
#' @param nu The mass times the angular frequency over twice hbar.
#' @param r The radial coordinate.
#' @param hw The energy quantum, hbar times omega.
#' @return A `sym` expression.
#' @examples
#' r <- Symbol("r"); nu <- Symbol("nu")
#' sho_R_nl(0, 0, 1, r)
#' sho_E_nl(1, 2, Symbol("hw"))
#' @export
sho_R_nl <- function(n, l, nu, r) {
  n <- as_sym(n); l <- as_sym(l); nu <- as_sym(nu); r <- as_sym(r)
  n <- n + 1
  C <- sqrt(((2 * nu)^(l + Rational(3, 2)) * 2^(n + l + 1) * fcall("factorial", list(n - 1))) /
              (sqrt(S_Pi) * fcall("factorial2", list(2 * n + 2 * l - 1))))
  C * r^l * exp(-nu * r^2) * fcall("assoc_laguerre", list(n - 1, l + S_Half, 2 * nu * r^2))
}

#' @rdname sho_R_nl
#' @export
sho_E_nl <- function(n, l, hw) {
  (2 * as_sym(n) + l + Rational(3, 2)) * hw
}

# ---------------------------------------------------------------------------
# One dimensional harmonic oscillator

#' The one dimensional quantum harmonic oscillator
#'
#' `qho_psi_n()` is the wavefunction, `qho_E_n()` the energy and
#' `qho_coherent_state()` the coefficient of the coherent state, after
#' `sympy.physics.qho_1d`. The reduced Planck constant appears as `hbar`.
#'
#' @param n The quantum number.
#' @param x The coordinate.
#' @param m The mass.
#' @param omega The angular frequency.
#' @param alpha The eigenvalue of the annihilation operator.
#' @return A `sym` expression.
#' @examples
#' x <- Symbol("x")
#' qho_psi_n(0, x, 1, 1)
#' qho_E_n(1, Symbol("omega"))
#' @export
qho_psi_n <- function(n, x, m, omega) {
  n <- as_sym(n); x <- as_sym(x); m <- as_sym(m); omega <- as_sym(omega)
  nu <- m * omega / .HBAR
  C <- (nu / S_Pi)^Rational(1, 4) * sqrt(1 / (2^n * fcall("factorial", list(n))))
  C * exp(-nu * x^2 / 2) * fcall("hermite", list(n, sqrt(nu) * x))
}

#' @rdname qho_psi_n
#' @export
qho_E_n <- function(n, omega) {
  .HBAR * omega * (as_sym(n) + S_Half)
}

#' @rdname qho_psi_n
#' @export
qho_coherent_state <- function(n, alpha) {
  alpha <- as_sym(alpha)
  exp(-Abs(alpha)^2 / 2) * (alpha^n) / sqrt(fcall("factorial", list(as_sym(n))))
}

# ---------------------------------------------------------------------------
# Particle on a ring

#' The particle on a ring
#'
#' `pring_wavefunction()` is the normalized wavefunction and
#' `pring_energy()` the energy of a particle on a ring, after
#' `sympy.physics.pring`.
#'
#' @param n The quantum number.
#' @param x The angle.
#' @param m The mass.
#' @param r The radius of the ring.
#' @return A `sym` expression.
#' @examples
#' x <- Symbol("x")
#' pring_wavefunction(1, x)
#' pring_energy(1, Symbol("m"), Symbol("r"))
#' @export
pring_wavefunction <- function(n, x) {
  n <- as_sym(n); x <- as_sym(x)
  exp(n * S_I * x) / sqrt(2 * S_Pi)
}

#' @rdname pring_wavefunction
#' @export
pring_energy <- function(n, m, r) {
  n <- as_sym(n); m <- as_sym(m); r <- as_sym(r)
  if (isTRUE(ask_fact(n, "integer"))) return((n^2 * .HBAR^2) / (2 * m * r^2))
  .physics_error("ValueError", "'n' must be integer")
}
