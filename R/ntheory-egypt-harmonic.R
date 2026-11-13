# The harmonic prefix of Egyptian fractions, after egypt_harmonic() in
# sympy/ntheory/egyptian_fraction.py.

#' Harmonic prefix of an Egyptian fraction
#'
#' `egypt_harmonic(r)` takes the terms `1, 1/2, 1/3, ...` of the harmonic
#' series while their sum stays at most `r`, and returns the list of the
#' denominators taken and the remainder, as SymPy's `egypt_harmonic()`
#' (the first step of [egyptian_fraction()] for rationals above 1).
#'
#' @param r A rational number.
#' @return A list of a list of integers and a rational.
#' @examples
#' egypt_harmonic(Rational(5, 2))
#' egypt_harmonic(Rational(1, 3))
#' @export
egypt_harmonic <- function(r) {
  r <- as_sym(r)
  if (!.is_rat(r)) stop(.sympy_condition("TypeError", "expected a Rational"))
  q <- gmp::as.bigq(r$value)
  rv <- list()
  d <- gmp::as.bigz(1)
  acc <- gmp::as.bigq(0)
  while (acc + 1 / gmp::as.bigq(d) <= q) {
    acc <- acc + 1 / gmp::as.bigq(d)
    rv[[length(rv) + 1L]] <- .int(d)
    d <- d + 1
  }
  list(rv, make_rational(q - acc))
}
