# Univariate power series rings over ZZ and QQ, a port of
# sympy/polys/series/ringpython.py (PythonPowerSeriesRingZZ and
# PythonPowerSeriesRingQQ). A series is the pair (coeffs, prec): coeffs is a
# dense coefficient list, highest degree first and without leading zeros,
# and prec is the order of the error term O(x**prec), or NULL (None) for an
# exact polynomial. The ring operations truncate at the ring precision.

# dense coefficient helpers (gmp vectors, highest degree first)
.us_K <- function(dom, v) if (dom == "ZZ") gmp::as.bigz(v) else gmp::as.bigq(v)
.us_empty <- function(dom) .us_K(dom, integer(0))
.us_strip <- function(f) {
  i <- 1L
  while (i <= length(f) && f[i] == 0) i <- i + 1L
  if (i > length(f)) f[0] else f[i:length(f)]
}
.us_deg <- function(f) if (length(f)) length(f) - 1L else -Inf
.us_trunc <- function(f, n) {
  if (n <= 0) return(f[0])
  if (length(f) > n) .us_strip(f[(length(f) - n + 1L):length(f)]) else f
}
.us_add <- function(f, g) {
  n <- max(length(f), length(g))
  if (!n) return(f)
  a <- .us_lift(f, n); b <- .us_lift(g, n)
  .us_strip(a + b)
}
.us_lift <- function(f, n) {
  z <- if (inherits(f, "bigq")) gmp::as.bigq(rep(0, n)) else gmp::as.bigz(rep(0, n))
  if (length(f)) z[(n - length(f) + 1L):n] <- f
  z
}
.us_neg <- function(f) if (length(f)) -f else f
.us_sub <- function(f, g) .us_add(f, .us_neg(g))
.us_mul <- function(f, g) {
  if (!length(f) || !length(g)) return(f[0])
  n <- length(f) + length(g) - 1L
  r <- .us_lift(f[0], n)
  if (inherits(g, "bigq") && !inherits(r, "bigq")) r <- gmp::as.bigq(r)
  for (i in seq_along(f)) for (j in seq_along(g)) r[i + j - 1L] <- r[i + j - 1L] + f[i] * g[j]
  .us_strip(r)
}
.us_series_mul <- function(f, g, n) .us_trunc(.us_mul(.us_trunc(f, n), .us_trunc(g, n)), n)
.us_series_pow <- function(f, k, n) {
  r <- .us_one_like(f)
  b <- .us_trunc(f, n)
  while (k > 0) {
    if (k %% 2 == 1) r <- .us_series_mul(r, b, n)
    k <- k %/% 2
    if (k > 0) b <- .us_series_mul(b, b, n)
  }
  .us_trunc(r, n)
}
.us_one_like <- function(f) if (inherits(f, "bigq")) gmp::as.bigq(1) else gmp::as.bigz(1)
.us_mul_ground <- function(f, c) if (c == 0 || !length(f)) f[0] else .us_strip(f * c)
.us_add_ground <- function(f, c) {
  if (!length(f)) return(.us_strip(c(f[0], c)))
  f[length(f)] <- f[length(f)] + c
  .us_strip(f)
}
.us_diff <- function(f) {
  n <- length(f) - 1L
  if (n < 1L) return(f[0])
  .us_strip(f[seq_len(n)] * gmp::as.bigz(n:1))
}
.us_integrate <- function(f) {
  n <- length(f)
  if (!n) return(f[0])
  .us_strip(c(gmp::as.bigq(f) / gmp::as.bigq(n:1), gmp::as.bigq(0)))
}
.us_rshift <- function(f, n) if (n >= length(f)) f[0] else f[seq_len(length(f) - n)]
# dup_exquo: exact division or ExactQuotientFailed
.us_exquo <- function(f, g, dom) {
  r <- f; q <- f[0]
  dq <- length(f) - length(g) + 1L
  if (dq < 1L) { if (length(f)) stop(ExactQuotientFailed("no exact quotient")); return(f) }
  q <- .us_lift(f[0], dq)
  if (dom == "QQ") q <- gmp::as.bigq(q)
  r <- .us_lift(f, length(f))
  if (dom == "QQ") r <- gmp::as.bigq(r)
  for (i in seq_len(dq)) {
    c <- r[i] / g[1L]
    if (dom == "ZZ") { if (gmp::denominator(gmp::as.bigq(c)) != 1) stop(ExactQuotientFailed("no exact quotient")); c <- gmp::as.bigz(c) }
    q[i] <- c
    r[i:(i + length(g) - 1L)] <- r[i:(i + length(g) - 1L)] - c * g
  }
  if (any(r != 0)) stop(ExactQuotientFailed("no exact quotient"))
  .us_strip(q)
}
# dup_revert: 1/f mod x**n (the constant term is a unit)
.us_revert <- function(f, n, dom) {
  a <- rev(f)
  b <- .us_lift(f[0], n)
  if (dom == "QQ") b <- gmp::as.bigq(b)
  inv0 <- if (dom == "ZZ") gmp::as.bigz(1 / gmp::as.bigq(a[1L])) else 1 / a[1L]
  if (n >= 1L) b[1L] <- inv0
  if (n >= 2L) for (k in 2:n) {
    s <- b[1L] * 0
    for (j in 2:min(k, length(a))) if (min(k, length(a)) >= 2L) s <- s + a[j] * b[k - j + 1L]
    b[k] <- -inv0 * s
  }
  .us_strip(rev(b))
}
.us_compose <- function(f, g, n) {
  r <- f[0]
  for (i in seq_along(f)) r <- .us_trunc(.us_add_ground(.us_series_mul(r, g, n), f[i]), n)
  r
}
# dup_series_reversion: the compositional inverse g of f mod x**n
.us_reversion <- function(f, n, dom) {
  if (n <= 1L) return(f[0])
  a1 <- rev(f)[2L]
  inv1 <- if (dom == "ZZ") gmp::as.bigz(1 / gmp::as.bigq(a1)) else 1 / a1
  g <- c(inv1, inv1 * 0)
  if (n >= 3L) for (k in 2:(n - 1L)) {
    h <- rev(.us_lift(.us_compose(f, g, k + 1L), k + 1L))
    g <- .us_sub(c(-inv1 * h[k + 1L], .us_lift(g[0], k)), .us_neg(.us_lift(g, k + 1L)))
    g <- .us_strip(g)
  }
  .us_trunc(g, n)
}

# _useries and the operations of ringpython.py
.uss <- function(coeffs, prec) .pytuple(list(coeffs, prec))
.us_c <- function(s) s[[1L]]
.us_p <- function(s) s[[2L]]
.us_new <- function(coeffs, prec, dom, rp) {
  deg <- .us_deg(coeffs)
  if (is.null(prec)) {
    if (deg < rp) return(.uss(coeffs, NULL))
    prec <- rp
  }
  prec <- min(prec, rp)
  if (deg < prec) return(.uss(coeffs, prec))
  .uss(.us_trunc(coeffs, prec), prec)
}
.us_valuation <- function(s) {
  cf <- .us_c(s)
  if (!length(cf)) return(-1L)
  i <- length(cf)
  while (cf[i] == 0) i <- i - 1L
  length(cf) - i
}
.us_unify <- function(s1, s2, rp) {
  c1 <- .us_c(s1); p1 <- .us_p(s1); c2 <- .us_c(s2); p2 <- .us_p(s2)
  if (!is.null(p1) && !is.null(p2) && p1 == p2 && p1 == rp) return(list(c1, c2, rp))
  if (!is.null(p1) && !is.null(p2) && p1 == p2 && p1 <= rp) return(list(c1, c2, p1))
  u <- if (is.null(p1) && is.null(p2)) NULL else if (is.null(p1)) p2 else if (is.null(p2)) p1 else min(p1, p2)
  u <- if (is.null(u)) rp else min(u, rp)
  if (.us_deg(c1) >= u) c1 <- .us_trunc(c1, u)
  if (.us_deg(c2) >= u) c2 <- .us_trunc(c2, u)
  list(c1, c2, u)
}
.us_same <- function(a, b) length(a) == length(b) && (!length(a) || all(a == b))
.us_equality <- function(s1, s2, rp) {
  if (is.null(.us_p(s1)) && is.null(.us_p(s2))) return(.us_same(.us_c(s1), .us_c(s2)))
  u <- .us_unify(s1, s2, rp)
  if (!.us_same(u[[1L]], u[[2L]])) return(FALSE)
  NA
}
.us_equal_repr <- function(s1, s2) .us_same(.us_c(s1), .us_c(s2)) && identical(.us_p(s1), .us_p(s2))
.us_pos <- function(s, rp) {
  cf <- .us_c(s); prec <- .us_p(s); deg <- .us_deg(cf)
  if (is.null(prec)) {
    if (deg >= rp) { cf <- .us_trunc(cf, rp); prec <- rp }
  } else {
    prec <- min(prec, rp)
    if (deg >= prec) cf <- .us_trunc(cf, prec)
  }
  .uss(cf, prec)
}
.us_negs <- function(s, rp) .us_pos(.uss(.us_neg(.us_c(s)), .us_p(s)), rp)
.us_addsub <- function(s1, s2, rp, op) {
  md <- max(.us_deg(.us_c(s1)), .us_deg(.us_c(s2)))
  if (is.null(.us_p(s1)) && is.null(.us_p(s2)) && md < rp) return(.uss(op(.us_c(s1), .us_c(s2)), NULL))
  u <- .us_unify(s1, s2, rp)
  .uss(op(u[[1L]], u[[2L]]), u[[3L]])
}
.us_add_ground_s <- function(s, n, dom, rp) {
  if (n == 0) return(s)
  .us_new(.us_add_ground(.us_c(s), n), .us_p(s), dom, rp)
}
.us_rsub_ground_s <- function(s, n, dom, rp) {
  if (n == 0) return(s)
  .us_new(.us_add_ground(.us_neg(.us_c(s)), n), .us_p(s), dom, rp)
}
.us_muls <- function(s1, s2, rp) {
  if (is.null(.us_p(s1)) && is.null(.us_p(s2))) {
    if (.us_deg(.us_c(s1)) + .us_deg(.us_c(s2)) < rp) return(.uss(.us_mul(.us_c(s1), .us_c(s2)), NULL))
  }
  u <- .us_unify(s1, s2, rp)
  .uss(.us_series_mul(u[[1L]], u[[2L]], u[[3L]]), u[[3L]])
}
.us_mul_ground_s <- function(s, n, dom, rp) {
  if (n == 0) return(.uss(.us_c(s)[0], .us_p(s)))
  if (n == 1) return(s)
  if (n == -1) return(.us_negs(s, rp))
  .us_new(.us_mul_ground(.us_c(s), n), .us_p(s), dom, rp)
}
.us_div_ground_s <- function(s, n, dom, rp) {
  if (n == 0) stop(.sympy_condition("ZeroDivisionError", "Division by zero in power series"))
  if (n == 1) return(s)
  if (n == -1) return(.us_negs(s, rp))
  .us_new(.us_mul_ground(.us_c(s), 1 / gmp::as.bigq(n)), .us_p(s), dom, rp)
}
.us_square_s <- function(s, rp) {
  prec <- .us_p(s)
  if (is.null(prec)) {
    if (.us_deg(.us_c(s)) * 2 < rp) return(.uss(.us_series_pow(.us_c(s), 2, rp), NULL))
    prec <- rp
  }
  .uss(.us_series_pow(.us_c(s), 2, prec), prec)
}
.us_is_unit <- function(c, dom) if (dom == "ZZ") abs(c) == 1 else c != 0
.us_inverse_s <- function(s, dom, rp) {
  cf <- .us_c(s); prec <- .us_p(s)
  if (!length(cf) || !.us_is_unit(cf[length(cf)], dom))
    stop(.sympy_condition("NotReversible", "Series inverse requires the constant term to be a unit"))
  if (is.null(prec)) {
    if (length(cf) == 1L) return(.uss(.us_revert(cf, 1L, dom), NULL))
    prec <- rp
  } else prec <- min(prec, rp)
  .uss(.us_trunc(.us_revert(cf, prec, dom), prec), prec)
}
.us_div_direct <- function(s1, s2, dom, rp) {
  if (!.us_is_unit(.us_c(s2)[length(.us_c(s2))], dom)) stop(.sympy_condition("ValueError", "Trailing coefficient of the divisor must be a unit"))
  mp <- .us_unify(s1, s2, rp)[[3L]]
  .us_muls(s1, .us_inverse_s(s2, dom, mp), rp)
}
.us_div_s <- function(s1, s2, dom, rp) {
  c1 <- .us_c(s1); p1 <- .us_p(s1); c2 <- .us_c(s2); p2 <- .us_p(s2)
  if (!length(c2)) stop(.sympy_condition("ZeroDivisionError", "Series division by zero"))
  v1 <- .us_valuation(s1); v2 <- .us_valuation(s2)
  if (v1 < v2) stop(.sympy_condition("ValueError", "quotient would not be a power series"))
  if (is.null(p1) && is.null(p2)) {
    q <- tryCatch(.us_exquo(c1, c2, dom), ExactQuotientFailed = function(e) NULL)
    if (!is.null(q)) return(.uss(q, NULL))
  }
  if (v2 == 0L) return(.us_div_direct(s1, s2, dom, rp))
  p1 <- if (!is.null(p1)) p1 - v2 else rp
  p2 <- if (!is.null(p2)) p2 - v2 else rp
  .us_div_direct(.us_new(.us_rshift(c1, v2), p1, dom, rp), .us_new(.us_rshift(c2, v2), p2, dom, rp), dom, rp)
}
.us_pow_s <- function(s, n, dom, rp) {
  if (n < 0) {
    s <- .us_pow_s(s, -n, dom, rp)
    return(tryCatch(.us_inverse_s(s, dom, rp),
                    NotReversible = function(e) stop(.sympy_condition("ValueError", "Result would not be a power series"))))
  }
  cf <- .us_c(s); prec <- .us_p(s)
  if (n == 0) return(.uss(.us_one_like(cf), prec))
  if (is.null(prec)) {
    if (.us_deg(cf) * n < rp) return(.uss(.us_series_pow(cf, n, rp), NULL))
    prec <- rp
  } else prec <- min(prec, rp)
  .uss(.us_series_pow(cf, n, prec), prec)
}
.us_truncate_s <- function(s, n) {
  if (n < 0) stop(.sympy_condition("ValueError", "Truncation precision must be non-negative"))
  if (.us_deg(.us_c(s)) < n) return(s)
  .uss(.us_trunc(.us_c(s), n), n)
}
.us_compose_s <- function(s1, s2, dom, rp) {
  c1 <- .us_c(s1); c2 <- .us_c(s2)
  if (length(c2) && c2[length(c2)] != 0)
    stop(.sympy_condition("ValueError", "Series composition requires the constant term of the second series to be zero."))
  if (is.null(.us_p(s1)) && is.null(.us_p(s2))) {
    comp <- .us_compose(c1, c2, rp)
    return(.uss(comp, if (.us_deg(c1) * .us_deg(c2) < rp) NULL else rp))
  }
  u <- .us_unify(s1, s2, rp)
  .uss(.us_compose(u[[1L]], u[[2L]], u[[3L]]), u[[3L]])
}
.us_reversion_s <- function(s, dom, rp) {
  cf <- .us_c(s); prec <- .us_p(s)
  if (!length(cf) || cf[length(cf)] != 0)
    stop(.sympy_condition("NotReversible", "Series compositional inverse requires the constant term to be zero."))
  if (length(cf) >= 2L && !.us_is_unit(cf[length(cf) - 1L], dom))
    stop(.sympy_condition("NotReversible", "Series compositional inverse requires the linear term to be unit."))
  prec <- if (is.null(prec)) rp else min(prec, rp)
  .uss(.us_reversion(cf, prec, dom), prec)
}
.us_derivative_s <- function(s, dom, rp) {
  prec <- .us_p(s)
  if (!is.null(prec) && prec != 0) prec <- prec - 1L
  .us_new(.us_diff(.us_c(s)), prec, dom, rp)
}
.us_integrate_s <- function(s, dom, rp) {
  prec <- .us_p(s)
  if (!is.null(prec) && prec != 0) prec <- prec + 1L
  .us_new(.us_integrate(.us_c(s)), prec, dom, rp)
}
# exsqrt over QQ: the exact square root or NULL
.us_exsqrt <- function(c) {
  q <- gmp::as.bigq(c)
  if (q < 0) return(NULL)
  n <- gmp::numerator(q); d <- gmp::denominator(q)
  rn <- .isqrt(n); rd <- .isqrt(d)
  if (rn * rn != n || rd * rd != d) return(NULL)
  gmp::as.bigq(rn, rd)
}
# the square root of a series with constant term one, mod x**n
.us_sqrt1 <- function(f, n) {
  a <- rev(.us_lift(gmp::as.bigq(f), max(n, length(f))))
  b <- gmp::as.bigq(rep(0, n))
  if (n >= 1L) b[1L] <- 1
  if (n >= 2L) for (k in 2:n) {
    s <- a[k]
    if (k >= 3L) for (j in 2:(k - 1L)) s <- s - b[j] * b[k - j + 1L]
    b[k] <- s / 2
  }
  .us_strip(rev(b))
}
.us_sqrt_s <- function(s, dom, rp) {
  cf <- .us_c(s); prec <- .us_p(s)
  if (!length(cf)) return(s)
  c0 <- cf[length(cf)]
  if (c0 == 0) stop(.sympy_condition("ValueError", "Cannot compute square root of a series with zero constant term."))
  r0 <- .us_exsqrt(c0)
  if (is.null(r0)) stop(.sympy_condition("ValueError", "Constant term is not a perfect square in the domain."))
  ns <- .us_div_ground_s(s, c0, dom, rp)
  if (is.null(prec)) {
    ds <- .us_deg(cf)
    if (ds %% 2 == 0) {
      root <- .us_sqrt1(.us_c(ns), ds %/% 2 + 1L)
      sq <- .us_square_s(.uss(root, ds + 1L), ds + 1L)
      if (.us_same(.us_c(sq), .us_c(ns))) return(.us_mul_ground_s(.uss(root, NULL), r0, dom, rp))
    }
  }
  ap <- if (!is.null(prec)) min(prec, rp) else rp
  sq <- .us_mul_ground_s(.uss(.us_sqrt1(.us_c(ns), ap), ap), r0, dom, ap)
  .uss(.us_c(sq), ap)
}
.us_check0 <- function(s, msg) {
  cf <- .us_c(s)
  if (length(cf) && cf[length(cf)] != 0) stop(.sympy_condition("ValueError", msg))
}
.us_prec0 <- function(s, rp) if (is.null(.us_p(s))) rp else min(.us_p(s), rp)
# compose the series sum(c[k] x**k) (ascending c) with s
.us_fcompose <- function(c, s, prec, dom, rp) .us_compose_s(.uss(.us_strip(rev(gmp::as.bigq(c))), prec), .uss(.us_c(s), prec), dom, rp)
.us_coef_series <- function(kind, n) {
  c <- gmp::as.bigq(rep(0, max(n, 1L)))
  fact <- function(k) gmp::factorialZ(k)
  for (k in 0:(n - 1L)) {
    v <- switch(kind,
      log = if (k == 0) 0 else gmp::as.bigq((-1)^(k - 1), k),
      exp = 1 / gmp::as.bigq(fact(k)),
      sin = if (k %% 2 == 1) gmp::as.bigq((-1)^((k - 1) / 2)) / gmp::as.bigq(fact(k)) else 0,
      sinh = if (k %% 2 == 1) 1 / gmp::as.bigq(fact(k)) else 0,
      cos = if (k %% 2 == 0) gmp::as.bigq((-1)^(k / 2)) / gmp::as.bigq(fact(k)) else 0,
      cosh = if (k %% 2 == 0) 1 / gmp::as.bigq(fact(k)) else 0)
    c[k + 1L] <- v
  }
  c[seq_len(n)]
}
.us_elem <- function(kind, s, dom, rp) {
  cf <- .us_c(s)
  msgs <- c(exp = "Exponential", sin = "Sine", sinh = "Hyperbolic sine", cos = "Cosine", cosh = "Hyperbolic cosine",
            atan = "Arctangent", atanh = "Hyperbolic Arctangent", asin = "Hyperbolic arcsine", asinh = "Hyperbolic Arcsine",
            tan = "Tangent", tanh = "Hyperbolic Tangent")
  if (!length(cf)) return(if (kind %in% c("exp", "cos", "cosh")) .uss(gmp::as.bigq(1), .us_p(s)) else s)
  .us_check0(s, paste0(msgs[[kind]], " requires the constant term of the input series to be zero."))
  prec <- .us_prec0(s, rp)
  if (kind %in% c("exp", "sin", "sinh", "cos", "cosh")) return(.us_fcompose(.us_coef_series(kind, prec), s, prec, dom, rp))
  s <- .uss(cf, prec)
  if (kind %in% c("atan", "atanh")) {
    ds <- .us_derivative_s(s, dom, rp)
    s2 <- .us_square_s(s, rp)
    q <- if (kind == "atan") .us_add_ground_s(s2, 1, dom, rp) else .us_rsub_ground_s(s2, 1, dom, rp)
    return(.us_integrate_s(.us_muls(ds, .us_inverse_s(q, dom, rp), rp), dom, rp))
  }
  if (kind %in% c("asin", "asinh")) {
    c <- gmp::as.bigq(c(0, 1, 0))
    k <- 3L
    while (k < prec) {
      t <- gmp::as.bigq((k - 2)^2) * c[length(c) - 1L] / gmp::as.bigq(k * (k - 1))
      c <- c(c, if (kind == "asin") t else -t, gmp::as.bigq(0))
      k <- k + 2L
    }
    return(.us_compose_s(.uss(.us_strip(rev(c)), prec), s, dom, rp))
  }
  # tan and tanh: Newton's iteration converges to the unique series
  sn <- .us_coef_series(if (kind == "tan") "sin" else "sinh", prec)
  cs <- .us_coef_series(if (kind == "tan") "cos" else "cosh", prec)
  tq <- .us_series_mul(.us_strip(rev(sn)), .us_revert(.us_strip(rev(cs)), prec, "QQ"), prec)
  .us_compose_s(.uss(tq, prec), s, dom, rp)
}
.us_log_s <- function(s, dom, rp) {
  cf <- .us_c(s)
  if (!length(cf)) stop(.sympy_condition("ValueError", "Logarithm of a zero series is undefined."))
  if (cf[length(cf)] != 1) stop(.sympy_condition("ValueError", "Logarithm of a power series requires the constant term to be one."))
  if (length(cf) == 1L) return(.uss(cf[0], .us_p(s)))
  prec <- .us_prec0(s, rp)
  s <- .us_add_ground_s(.uss(cf, prec), -1, dom, rp)
  .us_compose_s(.uss(.us_strip(rev(.us_coef_series("log", prec))), prec), s, dom, rp)
}
.us_pretty <- function(s, sym = "x", ascending = TRUE) {
  cf <- .us_c(s); prec <- .us_p(s)
  if (!is.null(prec) && prec == 0) return(paste0("O(", sym, "**", prec, ")"))
  deg <- length(cf) - 1L
  terms <- character()
  for (i in seq_along(cf)) {
    c <- cf[i]; e <- deg - i + 1L
    if (c == 0) next
    cs <- if (inherits(c, "bigq")) .de_str(c) else as.character(c)
    t <- if (e == 0L) cs else {
      base <- if (e == 1L) sym else paste0(sym, "**", e)
      if (c == 1) base else if (c == -1) paste0("-", base) else paste0(cs, "*", base)
    }
    terms <- c(terms, t)
  }
  if (ascending) terms <- rev(terms)
  poly <- gsub("+ -", "- ", paste(terms, collapse = " + "), fixed = TRUE)
  if (!nzchar(poly)) poly <- "0"
  if (poly == "0") return(if (!is.null(prec)) paste0("O(", sym, "**", prec, ")") else "0")
  if (!is.null(prec)) paste0(poly, " + O(", sym, "**", prec, ")") else poly
}
.us_in <- function(dom, coeffs) {
  v <- lapply(coeffs, function(c) {
    if (inherits(c, c("bigz", "bigq"))) return(.us_K(dom, c))
    if (is.numeric(c) && length(c) == 2L && dom == "QQ") return(gmp::as.bigq(c[1L], c[2L]))
    if (is.numeric(c) && length(c) == 1L && c == round(c)) return(.us_K(dom, c))
    stop(.sympy_condition("TypeError", "Unsupported coefficient type"))
  })
  if (!length(v)) .us_empty(dom) else do.call(c, v)
}
.us_from_list <- function(dom, coeffs, prec, rp) {
  cf <- .us_strip(rev(.us_in(dom, coeffs)))
  n <- length(.us_in(dom, coeffs))
  if (is.null(prec)) {
    if (n <= rp) return(.uss(cf, NULL))
    prec <- rp
  }
  if (n > prec) cf <- .us_trunc(cf, prec)
  .uss(cf, prec)
}
.us_arg <- function(dom, s) {
  cf <- .us_c(s)
  cf <- if (is.list(cf)) { if (length(cf)) do.call(c, cf) else .us_empty(dom) } else cf
  .uss(.us_K(dom, cf), .us_p(s))
}
.us_out <- function(s) {
  cf <- .us_c(s)
  .uss(lapply(seq_along(cf), function(i) cf[i]), if (is.null(.us_p(s))) NULL else as.integer(.us_p(s)))
}

#' Power series rings over the integers and the rationals
#'
#' `PythonPowerSeriesRingZZ(prec)` and `PythonPowerSeriesRingQQ(prec)` are
#' univariate power series rings with coefficients in ZZ or QQ, as SymPy's
#' classes of the same names. A series is the pair `(coeffs, prec)`: the
#' coefficients from the highest degree down, and the order of the error
#' term `O(x**prec)` or `NULL` for an exact polynomial. Call the ring on a
#' vector of coefficients (lowest degree first) to make a series, and use
#' its methods with `$`: `add`, `subtract`, `multiply`, `divide`, `pow_int`,
#' `square`, `compose`, `inverse`, `reversion`, `truncate`,
#' `differentiate`, `pretty`, `to_list`, and over QQ `integrate`, `sqrt`,
#' `log`, `exp`, `sin`, `cos`, `tan`, `atan`, `hypot` and the hyperbolic functions.
#'
#' @param prec The ring precision, a nonnegative integer (default 6).
#' @param x A power series ring.
#' @param name The method name.
#' @param ... Unused.
#' @return A power series ring. Called on coefficients it returns a series.
#' @examples
#' R <- PythonPowerSeriesRingZZ()
#' s <- R(c(1, 2, 3))
#' R$pretty(R$pow_int(s, 2))
#' R$pretty(R$inverse(R(c(1, 1))))
#' Q <- PythonPowerSeriesRingQQ(4)
#' Q$pretty(Q$exp(Q$gen))
#' @name power_series_rings
NULL

.us_ring <- function(dom, prec) {
  if (prec < 0) stop(.sympy_condition("ValueError", "Power series precision must be non-negative"))
  prec <- as.integer(prec)
  rp <- prec
  f <- function(coeffs, prec = NULL) .us_out(.us_from_list(dom, coeffs, prec, rp))
  structure(f, class = "PythonPowerSeriesRing", dom = dom, prec = prec)
}
#' @rdname power_series_rings
#' @export
PythonPowerSeriesRingZZ <- function(prec = 6L) .us_ring("ZZ", prec)
#' @rdname power_series_rings
#' @export
PythonPowerSeriesRingQQ <- function(prec = 6L) .us_ring("QQ", prec)
#' @rdname power_series_rings
#' @export
`$.PythonPowerSeriesRing` <- function(x, name) {
  dom <- attr(x, "dom"); rp <- attr(x, "prec")
  A <- function(s) .us_arg(dom, s)
  O <- .us_out
  g <- function(c) .us_K(dom, c)
  fn <- function(kind) function(s) O(.us_elem(kind, A(s), dom, rp))
  qq_only <- c("integrate", "sqrt", "log", "log1p", "exp", "expm1", "sin", "sinh", "cos", "cosh", "tan", "tanh",
               "atan", "atanh", "asin", "asinh", "hypot")
  if (dom == "ZZ" && name %in% qq_only) name <- paste0("\001", name)
  switch(name,
    domain = if (dom == "ZZ") ZZ else QQ,
    prec = rp,
    one = O(if (rp == 0L) .uss(.us_empty(dom), 0L) else .uss(g(1), NULL)),
    zero = O(if (rp == 0L) .uss(.us_empty(dom), 0L) else .uss(.us_empty(dom), NULL)),
    gen = O(if (rp < 2L) .uss(.us_empty(dom), rp) else .uss(g(c(1, 0)), NULL)),
    pretty = function(s, symbol = "x", ascending = TRUE) .us_pretty(A(s), symbol, ascending),
    print = function(s, symbol = "x", ascending = TRUE) { cat(.us_pretty(A(s), symbol, ascending), "\n", sep = ""); invisible(NULL) },
    from_list = function(coeffs, prec = NULL) O(.us_from_list(dom, coeffs, prec, rp)),
    from_element = function(s) O(.us_new(.us_c(A(s)), .us_p(A(s)), dom, rp)),
    to_list = function(s) { cf <- rev(.us_c(A(s))); lapply(seq_along(cf), function(i) cf[i]) },
    to_dense = function(s) { cf <- .us_c(A(s)); lapply(seq_along(cf), function(i) cf[i]) },
    series_prec = function(s) .us_p(s),
    equal = function(s1, s2) .us_equality(A(s1), A(s2), rp),
    equal_repr = function(s1, s2) .us_equal_repr(A(s1), A(s2)),
    is_ground = function(s) if (rp == 0L) NA else length(.us_c(A(s))) <= 1L,
    constant_coefficient = function(s) { cf <- .us_c(A(s)); if (length(cf)) cf[length(cf)] else g(0) },
    positive = function(s) O(.us_pos(A(s), rp)),
    negative = function(s) O(.us_negs(A(s), rp)),
    add = function(s1, s2) O(.us_addsub(A(s1), A(s2), rp, .us_add)),
    add_ground = function(s, n) O(.us_add_ground_s(A(s), g(n), dom, rp)),
    subtract = function(s1, s2) O(.us_addsub(A(s1), A(s2), rp, .us_sub)),
    subtract_ground = function(s, n) O(if (n == 0) A(s) else .us_new(.us_add_ground(.us_c(A(s)), -g(n)), .us_p(A(s)), dom, rp)),
    rsubtract_ground = function(s, n) O(.us_rsub_ground_s(A(s), g(n), dom, rp)),
    multiply = function(s1, s2) O(.us_muls(A(s1), A(s2), rp)),
    multiply_ground = function(s, n) O(.us_mul_ground_s(A(s), g(n), dom, rp)),
    divide = function(s1, s2) O(.us_div_s(A(s1), A(s2), dom, rp)),
    pow_int = function(s, n) O(.us_pow_s(A(s), n, dom, rp)),
    square = function(s) O(.us_muls(A(s), A(s), rp)),
    compose = function(s1, s2) O(.us_compose_s(A(s1), A(s2), dom, rp)),
    inverse = function(s) O(.us_inverse_s(A(s), dom, rp)),
    reversion = function(s) O(.us_reversion_s(A(s), dom, rp)),
    truncate = function(s, n) O(.us_truncate_s(A(s), n)),
    differentiate = function(s) O(.us_derivative_s(A(s), dom, rp)),
    integrate = function(s) O(.us_integrate_s(A(s), dom, rp)),
    sqrt = function(s) O(.us_sqrt_s(A(s), dom, rp)),
    log = function(s) O(.us_log_s(A(s), dom, rp)),
    log1p = function(s) { s <- A(s); .us_check0(s, "Log1p requires the constant term of the input series to be zero."); O(.us_log_s(.us_add_ground_s(s, 1, dom, rp), dom, rp)) },
    exp = fn("exp"),
    expm1 = function(s) O(.us_add_ground_s(.us_elem("exp", A(s), dom, rp), -1, dom, rp)),
    sin = fn("sin"), sinh = fn("sinh"), cos = fn("cos"), cosh = fn("cosh"),
    tan = fn("tan"), tanh = fn("tanh"), atan = fn("atan"), atanh = fn("atanh"),
    asin = fn("asin"), asinh = fn("asinh"),
    hypot = function(s1, s2) O(.us_sqrt_s(.us_addsub(.us_muls(A(s1), A(s1), rp), .us_muls(A(s2), A(s2), rp), rp, .us_add), dom, rp)),
    stop(.sympy_condition("AttributeError", paste0("'PythonPowerSeriesRing", dom, "' object has no attribute '", sub("^\001", "", name), "'"))))
}
#' @rdname power_series_rings
#' @export
format.PythonPowerSeriesRing <- function(x, ...) paste0("Python Power Series Ring over ", attr(x, "dom"), " with precision ", attr(x, "prec"))
#' @rdname power_series_rings
#' @export
print.PythonPowerSeriesRing <- function(x, ...) { cat(format(x), "\n", sep = ""); invisible(x) }
