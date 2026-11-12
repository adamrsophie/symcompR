# taylor_term() and pade(): helpers over the series machinery of
# R/series-driver.R.

#' Taylor term
#'
#' `taylor_term(e, x, k, x0)` is the k-th Taylor term of `e` about `x0`.
#'
#' @param e An expression.
#' @param x The variable.
#' @param k Index of the term.
#' @param x0 The expansion point.
#' @return An expression.
#' @export
taylor_term <- function(e, x, k, x0 = 0) {
  e <- as_sym(e); x <- as_sym(x); x0 <- as_sym(x0)
  d <- e
  for (j in seq_len(k)) d <- .diff1(d, x)
  .mul(list(.div(subs(d, x, x0), make_rational(gmp::as.bigq(gmp::factorialZ(k)))), .pow(.sub(x, x0), .int(k))))
}

#' Pade approximant
#'
#' `pade(e, x, m, n)` returns the `[m/n]` Pade approximant of `e` at zero.
#'
#' @param e An expression.
#' @param x The variable.
#' @param m,n Degrees of numerator and denominator.
#' @return A rational function.
#' @export
pade <- function(e, x, m, n) {
  e <- as_sym(e); x <- as_sym(x)
  s <- expand(removeO(series(e, x, 0, m + n + 1)))
  cvec <- rep(list(S_Zero), m + n + 1L)
  for (t in (if (s$type == "Add") s$args else list(s))) {
    ce <- .so_as_coeff_exponent(t, x)
    k <- if (ce[[2L]]$type == "Integer") as.integer(num_value_dbl(ce[[2L]])) else -1L
    if (k >= 0L && k <= m + n) cvec[[k + 1L]] <- .add(list(cvec[[k + 1L]], ce[[1L]]))
  }
  qs <- lapply(seq_len(n), function(i) Dummy(paste0("q", i)))
  eqs <- lapply((m + 1):(m + n), function(k) .add(c(list(cvec[[k + 1L]]), lapply(seq_len(n), function(i)
    if (k - i >= 0) .mul(list(qs[[i]], cvec[[k - i + 1L]])) else S_Zero))))
  sol <- .linsolve_exprs(eqs, qs)
  qv <- lapply(sol, `[[`, 2L)
  Q <- .add(c(list(S_One), lapply(seq_len(n), function(i) .mul(list(qv[[i]], .pow(x, .int(i)))))))
  P <- .add(lapply(0:m, function(k) .mul(list(.add(c(list(cvec[[k + 1L]]), lapply(seq_len(min(k, n)), function(i)
    .mul(list(qv[[i]], cvec[[k - i + 1L]]))))), .pow(x, .int(k))))))
  .div(P, Q)
}

