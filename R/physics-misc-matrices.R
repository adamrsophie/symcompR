# Known matrices of physics, after sympy/physics/matrices.py.

.phys_mat <- function(rows) Matrix(lapply(rows, function(r) lapply(r, as_sym)))

#' Matrices of physics
#'
#' `msigma(i)` is the Pauli matrix with index 1, 2 or 3; `mgamma(mu, lower)`
#' the Dirac gamma matrix with index 0, 1, 2, 3 or 5 (upper index by
#' default; `lower = TRUE` lowers it with the Minkowski metric);
#' `minkowski_tensor` the Minkowski metric `diag(1, -1, -1, -1)`;
#' `pat_matrix(m, dx, dy, dz)` the parallel axis theorem matrix that
#' translates the inertia matrix of a point mass `m` by `(dx, dy, dz)`; and
#' `mdft(n)` the discrete Fourier transform matrix (deprecated in SymPy in
#' favor of `DFT(n)`). After `sympy.physics.matrices`.
#'
#' @param i The Pauli index.
#' @param mu The Dirac index.
#' @param lower Whether to lower the index.
#' @param m A mass.
#' @param dx,dy,dz The displacement.
#' @param n The size of the transform.
#' @return A `symmatrix`.
#' @examples
#' msigma(1) %*% msigma(2)
#' mgamma(5)
#' minkowski_tensor
#' pat_matrix(2, 1, 0, 0)
#' @export
msigma <- function(i) {
  mat <- if (identical(i, 1) || identical(i, 1L)) list(list(0, 1), list(1, 0))
         else if (identical(i, 2) || identical(i, 2L)) list(list(0, -S_I), list(S_I, 0))
         else if (identical(i, 3) || identical(i, 3L)) list(list(1, 0), list(0, -1))
         else NULL
  if (is_sym(i) && i$type == "Integer") return(msigma(as.integer(i$value)))
  if (is.null(mat)) .physics_error("IndexError", "Invalid Pauli index")
  .phys_mat(mat)
}

#' @rdname msigma
#' @export
pat_matrix <- function(m, dx, dy, dz) {
  dx <- as_sym(dx); dy <- as_sym(dy); dz <- as_sym(dz)
  dxdy <- -dx * dy; dydz <- -dy * dz; dzdx <- -dz * dx
  dxdx <- dx^2; dydy <- dy^2; dzdz <- dz^2
  mat <- list(list(dydy + dzdz, dxdy, dzdx),
              list(dxdy, dxdx + dzdz, dydz),
              list(dzdx, dydz, dydy + dxdx))
  .phys_mat(mat) * as_sym(m)
}

#' @rdname msigma
#' @export
mgamma <- function(mu, lower = FALSE) {
  if (is_sym(mu) && mu$type == "Integer") mu <- as.integer(mu$value)
  if (!(is.numeric(mu) && length(mu) == 1L && mu %in% c(0, 1, 2, 3, 5)))
    .physics_error("IndexError", "Invalid Dirac index")
  I <- S_I
  mat <- switch(as.character(mu),
    "0" = list(list(1, 0, 0, 0), list(0, 1, 0, 0), list(0, 0, -1, 0), list(0, 0, 0, -1)),
    "1" = list(list(0, 0, 0, 1), list(0, 0, 1, 0), list(0, -1, 0, 0), list(-1, 0, 0, 0)),
    "2" = list(list(0, 0, 0, -I), list(0, 0, I, 0), list(0, I, 0, 0), list(-I, 0, 0, 0)),
    "3" = list(list(0, 0, 1, 0), list(0, 0, 0, -1), list(-1, 0, 0, 0), list(0, 1, 0, 0)),
    "5" = list(list(0, 0, 1, 0), list(0, 0, 0, 1), list(1, 0, 0, 0), list(0, 1, 0, 0)))
  m <- .phys_mat(mat)
  if (isTRUE(lower) && mu %in% c(1, 2, 3, 5)) m <- -m
  m
}

#' @rdname msigma
#' @format NULL
#' @export
minkowski_tensor <- .phys_mat(list(list(1, 0, 0, 0), list(0, -1, 0, 0), list(0, 0, -1, 0), list(0, 0, 0, -1)))

#' @rdname msigma
#' @export
mdft <- function(n) {
  warning("The sympy.physics.matrices.mdft method is deprecated. Use sympy.DFT(n).as_explicit() instead.", call. = FALSE)
  n <- as.integer(n)
  w <- fcall("exp", list(.div(.mul(list(.int(-2), S_Pi, S_I)), .int(n))))
  Matrix(n, n, function(i, j) .div(.pow(w, .int((i - 1L) * (j - 1L))), .sqrt(.int(n))))
}
