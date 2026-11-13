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
#' `pade(e, x, m, n)` returns the `[m/n]` Pade approximant of `e` at `x0`,
#' computed from the extended Euclidean algorithm as in SymPy's
#' `pade_approximant()`.
#'
#' @param e An expression.
#' @param x The variable.
#' @param m,n Degrees of numerator and denominator.
#' @param x0 The expansion point.
#' @return A rational function.
#' @examples
#' x <- Symbol("x")
#' pade(exp(x), x, 2, 2)
#' @export
pade <- function(e, x, m, n = m, x0 = 0) {
  e <- as_sym(e); x <- as_sym(x); x0 <- as_sym(x0)
  m <- as.integer(m); n <- as.integer(n)
  if (m < 0L || n < 0L) stop(.sympy_condition("ValueError", sprintf("m=%d and n=%d must be non-negative integers", m, n)))
  # sympy.series.approximants.pade_approximant
  s <- expand(removeO(series(subs(e, x, .add(list(x, x0))), x, 0, m + n + 1L)))
  cf <- list(); ex <- integer()
  for (t in (if (s$type == "Add") s$args else list(s))) {
    ce <- .so_as_coeff_exponent(t, x)
    if (ce[[2L]]$type != "Integer" || .has1(ce[[1L]], x))
      stop(.sympy_condition("ValueError", sprintf("f=%s does not have a laurent expansion in x around x0.", .str(e))))
    ex <- c(ex, as.integer(num_value_dbl(ce[[2L]]))); cf <- c(cf, list(ce[[1L]]))
  }
  hp <- max(0L, -min(c(0L, ex)))
  if (n < hp) stop(.sympy_condition("ValueError", sprintf(
    "Expression f=%s has pole of order %d, which is greater than %d, the requested order of the denominator.", .str(e), hp, n)))
  ex <- ex + hp
  deg <- max(c(0L, ex))
  f <- rep(list(S_Zero), deg + 1L)
  for (i in seq_along(ex)) f[[deg - ex[i] + 1L]] <- .op_K(.add(list(f[[deg - ex[i] + 1L]], cf[[i]])))
  nd <- .pade_gcdex(.op_strip(f), m, n - hp)
  r <- .mul(list(.pow(x, .int(-hp)), .div(.op_expr(nd[[1L]], x), .op_expr(nd[[2L]], x))))
  if (is_zero_num(x0)) r else subs(r, x, .sub(x, x0))
}

# Dense coefficient lists, highest degree first; list() is the zero polynomial.
.op_strip <- function(f) {
  while (length(f) && is_zero_num(f[[1L]])) f <- f[-1L]
  f
}
.op_divmod <- function(f, g) {
  f <- .op_strip(f); g <- .op_strip(g)
  if (length(f) < length(g)) return(list(list(), f))
  q <- rep(list(S_Zero), length(f) - length(g) + 1L)
  for (i in seq_along(q)) {
    c <- .op_K(.div(f[[i]], g[[1L]]))
    q[[i]] <- c
    for (j in seq_along(g)) f[[i + j - 1L]] <- .op_K(.sub(f[[i + j - 1L]], .mul(list(c, g[[j]]))))
  }
  list(.op_strip(q), .op_strip(f))
}

# pade_approximant_gcdex(): the (r, t) pairs of gcdex_steps(x^(m+n+1), f)
.pade_gcdex <- function(f, m, n) {
  order <- m + n
  if (any(!vapply(f, is_zero_num, TRUE))) {
    nz <- which(!vapply(f, is_zero_num, TRUE))
    em <- length(f) - max(nz)
    if (em > m) stop(.sympy_condition("ValueError", sprintf(
      "polynomial has zero of order %d, which is greater than %d, the requested order of the numerator.", em, m)))
  }
  if (length(f) - 1L > order) f <- .op_strip(f[(length(f) - order):length(f)])
  r1 <- c(list(S_One), rep(list(S_Zero), order + 1L)); r2 <- f
  t1 <- list(); t2 <- list(S_One)
  out <- list(r2, t2)
  for (k in seq_len(order + 2L)) {
    if (k > 1L) {
      if (length(r2) - 1L < m) return(out)
      out <- list(r2, t2)
    }
    qr <- .op_divmod(r1, r2)
    if (!length(qr[[2L]])) break
    r1 <- r2; r2 <- qr[[2L]]
    qt <- if (length(qr[[1L]]) && length(t2)) .op_mul(qr[[1L]], t2) else list()
    nt <- .op_strip(.op_sub(t1, qt))
    t1 <- t2; t2 <- nt
  }
  out
}

