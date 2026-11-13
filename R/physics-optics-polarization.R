# Jones and Mueller calculus of polarized light, after
# sympy/physics/optics/polarization.py.

#' Jones and Stokes vectors and polarization optics
#'
#' `jones_vector(psi, chi)` and `stokes_vector(psi, chi, p, I)` are the
#' Jones and Stokes vectors of light with polarization ellipse angles `psi`
#' and `chi`, degree of polarization `p` and intensity `I`;
#' `jones_2_stokes(e)` converts a Jones vector to a Stokes vector.
#' `linear_polarizer()`, `phase_retarder()`, `half_wave_retarder()`,
#' `quarter_wave_retarder()`, `transmissive_filter()` and
#' `reflective_filter()` are Jones matrices of optical elements;
#' `mueller_matrix(J)` is the Mueller matrix of a Jones matrix and
#' `polarizing_beam_splitter()` the Jones matrix of a polarizing beam
#' splitter. After `sympy.physics.optics.polarization`.
#'
#' @param psi,chi The angles of the polarization ellipse.
#' @param p The degree of polarization.
#' @param I The intensity.
#' @param e A Jones vector.
#' @param theta The angle of the fast axis.
#' @param delta The phase difference.
#' @param T,R The transmittance and reflectance.
#' @param J A Jones matrix.
#' @param Tp,Rs,Ts,Rp The transmissivities and reflectivities of the two
#'   polarizations.
#' @param phia,phib The phase differences of the two reflections.
#' @return A `symmatrix`.
#' @examples
#' jones_vector(0, 0)
#' linear_polarizer(Symbol("theta"))
#' mueller_matrix(linear_polarizer(0))
#' @export
jones_vector <- function(psi, chi) {
  psi <- as_sym(psi); chi <- as_sym(chi)
  Matrix(list(list(-S_I * sin(chi) * sin(psi) + cos(chi) * cos(psi)),
              list(S_I * sin(chi) * cos(psi) + sin(psi) * cos(chi))))
}

#' @rdname jones_vector
#' @export
stokes_vector <- function(psi, chi, p = 1L, I = 1L) {
  psi <- as_sym(psi); chi <- as_sym(chi); p <- as_sym(p); I <- as_sym(I)
  S0 <- I
  S1 <- I * p * cos(2 * psi) * cos(2 * chi)
  S2 <- I * p * sin(2 * psi) * cos(2 * chi)
  S3 <- I * p * sin(2 * chi)
  Matrix(list(list(S0), list(S1), list(S2), list(S3)))
}

#' @rdname jones_vector
#' @export
jones_2_stokes <- function(e) {
  ex <- e$data[[1L]]; ey <- e$data[[2L]]
  Matrix(list(list(Abs(ex)^2 + Abs(ey)^2),
              list(Abs(ex)^2 - Abs(ey)^2),
              list(2 * re(ex * conjugate(ey))),
              list(-2 * im(ex * conjugate(ey)))))
}

#' @rdname jones_vector
#' @export
linear_polarizer <- function(theta = 0L) {
  theta <- as_sym(theta)
  Matrix(list(list(cos(theta)^2, sin(theta) * cos(theta)),
              list(sin(theta) * cos(theta), sin(theta)^2)))
}

#' @rdname jones_vector
#' @export
phase_retarder <- function(theta = 0L, delta = 0L) {
  theta <- as_sym(theta); delta <- as_sym(delta)
  R <- Matrix(list(list(cos(theta)^2 + exp(S_I * delta) * sin(theta)^2,
                        (1 - exp(S_I * delta)) * cos(theta) * sin(theta)),
                   list((1 - exp(S_I * delta)) * cos(theta) * sin(theta),
                        sin(theta)^2 + exp(S_I * delta) * cos(theta)^2)))
  R * exp(-S_I * delta / 2)
}

#' @rdname jones_vector
#' @export
half_wave_retarder <- function(theta) phase_retarder(theta, S_Pi)

#' @rdname jones_vector
#' @export
quarter_wave_retarder <- function(theta) phase_retarder(theta, S_Pi / 2)

#' @rdname jones_vector
#' @export
transmissive_filter <- function(T) {
  T <- as_sym(T)
  Matrix(list(list(sqrt(T), 0L), list(0L, sqrt(T))))
}

#' @rdname jones_vector
#' @export
reflective_filter <- function(R) {
  R <- as_sym(R)
  Matrix(list(list(sqrt(R), 0L), list(0L, -sqrt(R))))
}

#' @rdname jones_vector
#' @export
mueller_matrix <- function(J) {
  A <- Matrix(list(list(1L, 0L, 0L, 1L), list(1L, 0L, 0L, -1L), list(0L, 1L, 1L, 0L), list(0L, -S_I, S_I, 0L)))
  simplify(A %*% .optics_kron(J, .mat_map(J, function(z) conjugate(z))) %*% inverse(A))
}

# TensorProduct of two explicit matrices: the Kronecker product
.optics_kron <- function(A, B) {
  nr <- A$nrow * B$nrow; nc <- A$ncol * B$ncol
  Matrix(nr, nc, function(i, j) {
    .mul(list(.mget(A, (i - 1L) %/% B$nrow + 1L, (j - 1L) %/% B$ncol + 1L),
              .mget(B, (i - 1L) %% B$nrow + 1L, (j - 1L) %% B$ncol + 1L)))
  })
}

#' @rdname jones_vector
#' @export
polarizing_beam_splitter <- function(Tp = 1L, Rs = 1L, Ts = 0L, Rp = 0L, phia = 0L, phib = 0L) {
  Tp <- as_sym(Tp); Rs <- as_sym(Rs); Ts <- as_sym(Ts); Rp <- as_sym(Rp); phia <- as_sym(phia); phib <- as_sym(phib)
  Matrix(list(list(sqrt(Tp), 0L, S_I * sqrt(Rp), 0L),
              list(0L, sqrt(Ts), 0L, -S_I * sqrt(Rs) * exp(S_I * phia)),
              list(S_I * sqrt(Rp), 0L, sqrt(Tp), 0L),
              list(0L, -S_I * sqrt(Rs) * exp(S_I * phib), 0L, sqrt(Ts))))
}
