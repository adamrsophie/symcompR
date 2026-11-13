# Gaussian quadrature nodes and weights, after sympy/integrals/quadrature.py.
# Polynomials are gmp rational coefficient vectors, lowest degree first.

.qp_q <- function(x) gmp::as.bigq(x)

.qp_add <- function(a, b) {
  n <- max(length(a), length(b))
  a <- c(a, .qp_q(rep(0, n - length(a)))); b <- c(b, .qp_q(rep(0, n - length(b))))
  a + b
}

.qp_mul <- function(a, b) {
  out <- .qp_q(rep(0, length(a) + length(b) - 1L))
  for (i in seq_along(a)) for (j in seq_along(b)) out[i + j - 1L] <- out[i + j - 1L] + a[i] * b[j]
  out
}

.qp_eval <- function(p, r) {
  v <- .qp_q(0)
  for (i in rev(seq_along(p))) v <- v * r + p[i]
  v
}

.qp_diff <- function(p) if (length(p) < 2L) .qp_q(0) else p[-1L] * .qp_q(seq_len(length(p) - 1L))

# generalized binomial coefficient with a rational top
.qp_binom <- function(t, k) {
  v <- .qp_q(1)
  if (k > 0) for (i in 0:(k - 1)) v <- v * (t - i) / (i + 1)
  v
}

.qp_legendre <- function(n) {
  p0 <- .qp_q(1); if (n == 0) return(p0)
  p1 <- .qp_q(c(0, 1))
  if (n >= 2) for (k in 1:(n - 1)) {
    p2 <- .qp_add(.qp_mul(.qp_q(c(0, 2 * k + 1)), p1), -k * p0) / (k + 1)
    p0 <- p1; p1 <- p2
  }
  p1
}

.qp_hermite <- function(n) {
  p0 <- .qp_q(1); if (n == 0) return(p0)
  p1 <- .qp_q(c(0, 2))
  if (n >= 2) for (k in 1:(n - 1)) {
    p2 <- .qp_add(.qp_mul(.qp_q(c(0, 2)), p1), -2 * k * p0)
    p0 <- p1; p1 <- p2
  }
  p1
}

.qp_laguerre <- function(n, a = .qp_q(0)) {
  a <- .qp_q(a)
  out <- .qp_q(rep(0, n + 1))
  for (k in 0:n) {
    v <- .qp_q(1)
    if (k < n) for (j in (k + 1):n) v <- v * (a + j)
    out[k + 1] <- (-1)^k * v / (gmp::factorialZ(n - k) * gmp::factorialZ(k))
  }
  out
}

.qp_jacobi <- function(n, a, b) {
  a <- .qp_q(a); b <- .qp_q(b)
  out <- .qp_q(0)
  xm <- .qp_q(c(-1, 1)) / 2; xp <- .qp_q(c(1, 1)) / 2
  for (s in 0:n) {
    t <- .qp_q(.qp_binom(n + a, n - s) * .qp_binom(n + b, s))
    for (i in seq_len(s)) t <- .qp_mul(t, xm)
    for (i in seq_len(n - s)) t <- .qp_mul(t, xp)
    out <- .qp_add(out, t)
  }
  out
}

# the real roots of p in increasing order, as rationals within 10^-prec;
# exact when rational (a root with a small denominator that is exact)
.qp_real_roots <- function(p, digits) {
  while (length(p) > 1L && p[length(p)] == 0) p <- p[-length(p)]
  deg <- length(p) - 1L
  if (deg < 1L) return(list())
  dp <- .qp_diff(p)
  z <- polyroot(as.numeric(p))
  z <- sort(Re(z[abs(Im(z)) < 1e-6 * pmax(1, abs(z))]))
  scale <- gmp::as.bigz(10)^(digits + 40)
  lapply(z, function(x0) {
    r <- gmp::as.bigq(gmp::as.bigz(round(x0 * 2^40)), gmp::as.bigz(2)^40)
    for (it in 1:200) {
      rn <- r - .qp_eval(p, r) / .qp_eval(dp, r)
      rn <- gmp::as.bigq((gmp::numerator(rn) * scale) %/% gmp::denominator(rn), scale)
      done <- abs(rn - r) * scale <= 1
      r <- rn
      if (done) break
    }
    # a rational root with a small denominator is exact
    cf <- .qp_small_rational(r)
    if (!is.null(cf) && .qp_eval(p, cf) == 0) cf else r
  })
}

# the best rational approximation with denominator below 10^6, from the
# continued fraction expansion
.qp_small_rational <- function(r) {
  h0 <- gmp::as.bigz(0); h1 <- gmp::as.bigz(1); k0 <- gmp::as.bigz(1); k1 <- gmp::as.bigz(0)
  x <- r
  best <- NULL
  for (i in 1:40) {
    a <- gmp::numerator(x) %/% gmp::denominator(x)
    h2 <- a * h1 + h0; k2 <- a * k1 + k0
    if (k2 > 1000000) break
    best <- gmp::as.bigq(h2, k2)
    h0 <- h1; h1 <- h2; k0 <- k1; k1 <- k2
    fr <- x - a
    if (fr == 0) break
    x <- 1 / fr
  }
  best
}

.qp_n <- function(e, d) N(e, d)
.qp_r <- function(q) make_rational(gmp::as.bigq(q))

#' Gaussian quadrature
#'
#' Nodes and weights of Gaussian quadrature rules to `n_digits` significant
#' digits, as the functions of `sympy.integrals.quadrature`:
#' `gauss_legendre()` on `[-1, 1]`, `gauss_laguerre()` with weight
#' `exp(-x)` on `[0, oo)`, `gauss_hermite()` with weight `exp(-x^2)`,
#' `gauss_gen_laguerre()` with weight `x^alpha exp(-x)`,
#' `gauss_chebyshev_t()` and `gauss_chebyshev_u()` of the first and second
#' kind, `gauss_jacobi()` with weight `(1 - x)^alpha (1 + x)^beta` and
#' `gauss_lobatto()`, whose nodes include the end points. Each returns
#' `list(nodes, weights)`. The nodes are the real roots of the orthogonal
#' polynomial, refined with exact rational arithmetic far below the last
#' digit before rounding.
#'
#' @param n The number of nodes.
#' @param n_digits Significant digits.
#' @param alpha,beta Parameters of the weight function (rational).
#' @return A list of two lists of numbers.
#' @examples
#' gauss_legendre(3, 10)
#' gauss_chebyshev_t(2, 10)
#' @export
gauss_legendre <- function(n, n_digits) {
  p <- .qp_legendre(n); pd <- .qp_diff(p)
  rs <- .qp_real_roots(p, n_digits)
  list(lapply(rs, function(r) .qp_n(.qp_r(r), n_digits)),
       lapply(rs, function(r) .qp_n(.qp_r(2 / ((1 - r^2) * .qp_eval(pd, r)^2)), n_digits)))
}

#' @rdname gauss_legendre
#' @export
gauss_laguerre <- function(n, n_digits) {
  p <- .qp_laguerre(n); p1 <- .qp_laguerre(n + 1)
  rs <- .qp_real_roots(p, n_digits)
  list(lapply(rs, function(r) .qp_n(.qp_r(r), n_digits)),
       lapply(rs, function(r) .qp_n(.qp_r(r / ((n + 1)^2 * .qp_eval(p1, r)^2)), n_digits)))
}

#' @rdname gauss_legendre
#' @export
gauss_hermite <- function(n, n_digits) {
  p <- .qp_hermite(n); p1 <- .qp_hermite(n - 1)
  rs <- .qp_real_roots(p, n_digits)
  list(lapply(rs, function(r) .qp_n(.qp_r(r), n_digits)),
       lapply(rs, function(r) .qp_n(.mul(list(.qp_r(gmp::as.bigq(2)^(n - 1) * gmp::factorialZ(n) / (n^2 * .qp_eval(p1, r)^2)),
                                              .sqrt(S_Pi))), n_digits)))
}

#' @rdname gauss_legendre
#' @export
gauss_gen_laguerre <- function(n, alpha, n_digits) {
  a <- .qp_alpha(alpha)
  p <- .qp_laguerre(n, a); p1 <- .qp_laguerre(n - 1, a); p2 <- .qp_laguerre(n - 1, a + 1)
  g <- fcall("gamma", list(.add(list(.qp_r(a), .int(n)))))
  rs <- .qp_real_roots(p, n_digits)
  list(lapply(rs, function(r) .qp_n(.qp_r(r), n_digits)),
       lapply(rs, function(r) .qp_n(.mul(list(g, .qp_r(1 / (n * gmp::factorialZ(n - 1) * .qp_eval(p1, r) * .qp_eval(p2, r))))), n_digits)))
}

.qp_alpha <- function(a) {
  a <- as_sym(a)
  if (!.is_rat(a)) stop(.sympy_condition("NotImplementedError", "alpha and beta must be rational numbers here"))
  gmp::as.bigq(a$value)
}

#' @rdname gauss_legendre
#' @export
gauss_chebyshev_t <- function(n, n_digits) {
  list(lapply(seq_len(n), function(i) .qp_n(cos(.mul(list(.rat(2 * i - 1, 2 * n), S_Pi))), n_digits)),
       lapply(seq_len(n), function(i) .qp_n(.div(S_Pi, .int(n)), n_digits)))
}

#' @rdname gauss_legendre
#' @export
gauss_chebyshev_u <- function(n, n_digits) {
  list(lapply(seq_len(n), function(i) .qp_n(cos(.mul(list(.rat(i, n + 1), S_Pi))), n_digits)),
       lapply(seq_len(n), function(i) .qp_n(.mul(list(.div(S_Pi, .int(n + 1)), .pow(sin(.mul(list(.rat(i, n + 1), S_Pi))), S_Two))), n_digits)))
}

#' @rdname gauss_legendre
#' @export
gauss_jacobi <- function(n, alpha, beta, n_digits) {
  a <- .qp_alpha(alpha); b <- .qp_alpha(beta)
  p <- .qp_jacobi(n, a, b); pd <- .qp_diff(p); pn <- .qp_jacobi(n + 1, a, b)
  G <- .mul(list(fcall("gamma", list(.qp_r(n + a + 1))), fcall("gamma", list(.qp_r(n + b + 1))),
                 .pow(fcall("gamma", list(.qp_r(n + a + b + 1))), S_NegativeOne), .pow(S_Two, .qp_r(a + b))))
  k <- -(2 * n + a + b + 2) / (n + a + b + 1) / gmp::factorialZ(n + 1)
  rs <- .qp_real_roots(p, n_digits)
  list(lapply(rs, function(r) .qp_n(.qp_r(r), n_digits)),
       lapply(rs, function(r) .qp_n(.mul(list(G, .qp_r(k / (.qp_eval(pd, r) * .qp_eval(pn, r))))), n_digits)))
}

#' @rdname gauss_legendre
#' @export
gauss_lobatto <- function(n, n_digits) {
  p <- .qp_legendre(n - 1); pd <- .qp_diff(p)
  rs <- .qp_real_roots(pd, n_digits)
  xi <- lapply(rs, function(r) .qp_n(.qp_r(r), n_digits))
  w <- lapply(rs, function(r) .qp_n(.qp_r(2 / (n * (n - 1) * .qp_eval(p, r)^2)), n_digits))
  we <- .qp_n(.qp_r(gmp::as.bigq(2, n * (n - 1))), n_digits)
  list(c(list(S_NegativeOne), xi, list(S_One)), c(list(we), w, list(we)))
}
