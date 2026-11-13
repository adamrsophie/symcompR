# The Routh array of a polynomial, after
# sympy/physics/control/routh_table.py.

#' Routh-Hurwitz table
#'
#' `RouthHurwitz(polynomial, var)` is the Routh array of a polynomial: a
#' matrix with `degree + 1` rows whose first column decides the stability
#' of the polynomial. A row of zeros is replaced by the derivative of its
#' auxiliary polynomial and a zero in the first column is handled by the
#' shifting method. The result is a `symmatrix` with the extra fields
#' `zero_row_case`, `zero_col_infos` and `auxiliary_polynomials`. After
#' `sympy.physics.control.routh_table`.
#'
#' @param polynomial A polynomial expression.
#' @param var Its variable.
#' @return A `RouthHurwitz` matrix.
#' @examples
#' s <- Symbol("s")
#' RouthHurwitz(s^4 + s^3 + 3*s^2 + 3*s + 3, s)
#' @export
RouthHurwitz <- function(polynomial, var) {
  if (!(is_sym(var) && var$type == "Symbol")) .physics_error("ValueError", "var must be a Symbol")
  P <- Poly(polynomial, var)
  coeffs <- all_coeffs(P)
  deg <- length(coeffs) - 1L
  nr <- deg + 1L; nc <- deg %/% 2L + 1L
  d <- rep(list(S_Zero), nr * nc)
  at <- function(i, j) d[[i * nc + j + 1L]]
  put <- function(i, j, v) d[[i * nc + j + 1L]] <<- as_sym(v)
  is0 <- function(e) identical(e$key, S_Zero$key)
  zero_row_case <- FALSE
  zero_col_infos <- list()
  aux_degrees <- integer()
  special <- function(i) {
    active <- nc - i %/% 2L
    if (all(vapply(seq_len(active) - 1L, function(j) is0(at(i, j)), TRUE))) {
      zero_row_case <<- TRUE
      aux <- deg - i + 1L
      aux_degrees <<- c(aux_degrees, aux)
      for (j in seq_len(nc) - 1L) {
        cdeg <- aux - 2L * j
        if (cdeg < 0L) next
        put(i, j, .mul(list(at(i - 1L, j), .int(cdeg))))
      }
      return(invisible(NULL))
    }
    if (is0(at(i, 0L))) {
      n_zeros <- 0L
      for (j in seq_len(nc) - 1L) { if (!is0(at(i, j))) break; n_zeros <- n_zeros + 1L }
      zero_col_infos[[length(zero_col_infos) + 1L]] <<- .pytuple(list(i, n_zeros))
      if (n_zeros < active) {
        vals <- lapply(seq(n_zeros, active - 1L), function(j) at(i, j))
        sgn <- if (n_zeros %% 2L == 0L) S_One else S_NegativeOne
        for (k in seq_along(vals) - 1L) put(i, k, .add(list(at(i, k), .mul(list(sgn, vals[[k + 1L]])))))
      }
    }
    invisible(NULL)
  }
  if (deg < 1L) {
    put(0L, 0L, coeffs[[1L]])
  } else {
    row <- 0L; col <- 0L
    for (cf in coeffs) {
      put(row, col, cf)
      row <- (row + 1L) %% 2L
      col <- col + 1L - row
    }
    if (is0(at(1L, 0L))) special(1L)
    if (nr > 2L) for (i in 2:(nr - 1L)) {
      active <- nc - i %/% 2L
      for (j in seq_len(active) - 1L) {
        num <- .sub(.mul(list(at(i - 1L, 0L), at(i - 2L, j + 1L))), .mul(list(at(i - 2L, 0L), at(i - 1L, j + 1L))))
        put(i, j, .div(num, at(i - 1L, 0L)))
      }
      special(i)
    }
  }
  M <- .new_mat(d, nr, nc)
  aux_polys <- if (!zero_row_case) NULL else lapply(aux_degrees, function(ad) {
    r <- deg - ad
    ex <- seq(ad, 0L, by = -2L)
    Poly(.add(lapply(seq_along(ex), function(k) .mul(list(at(r, k - 1L), .pow(var, .int(ex[[k]])))))), var)
  })
  M$zero_row_case <- zero_row_case
  M$zero_col_infos <- zero_col_infos
  M$auxiliary_polynomials <- aux_polys
  class(M) <- c("RouthHurwitz", class(M))
  M
}
