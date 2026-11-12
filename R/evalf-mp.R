# Arbitrary precision binary floating point on gmp integers, after mpmath's
# libmpf and libelefun (the engine behind SymPy's evalf). A number is
# list(m = bigz mantissa, e = exponent) with value m * 2^e; results are
# rounded to p bits, half to even. Complex numbers are list(re, im).
# Elementary functions work in fixed point (bigz scaled by 2^wp) with guard
# bits and argument reduction.

.Z0 <- gmp::as.bigz(0)
.Z1 <- gmp::as.bigz(1)
.mpf <- function(m, e) list(m = gmp::as.bigz(m), e = e)
.mp_zero <- list(m = gmp::as.bigz(0), e = 0)
.mp_is_zero <- function(a) a$m == 0
.pow2 <- function(k) gmp::as.bigz(2)^k

# round(m / 2^s) to nearest, ties to even
.shr_round <- function(m, s) {
  if (s <= 0) return(m * .pow2(-s))
  d <- .pow2(s)
  q <- gmp::divq.bigz(m, d)
  r2 <- 2 * (m - q * d)
  if (r2 > d || (r2 == d && gmp::mod.bigz(q, 2) == 1)) q <- q + 1
  q
}
.mp_round <- function(m, e, p) {
  m <- gmp::as.bigz(m)
  if (m == 0) return(.mp_zero)
  b <- .bitlen(m)
  if (b > p) {
    s <- b - p
    m <- .shr_round(m, s); e <- e + s
    if (.bitlen(m) > p) { m <- .shr_round(m, 1); e <- e + 1 }
  }
  # strip trailing zero bits, as mpmath's normalized form
  if (m != 0) {
    tz <- 0
    while (gmp::mod.bigz(m, 2) == 0) { m <- gmp::divq.bigz(m, 2); tz <- tz + 1 }
    e <- e + tz
  }
  list(m = m, e = e)
}

.mp_from_q <- function(q, p) {
  q <- gmp::as.bigq(q)
  n <- gmp::numerator(q); d <- gmp::denominator(q)
  if (n == 0) return(.mp_zero)
  k <- p + .bitlen(d) - .bitlen(n) + 4
  num <- if (k >= 0) n * .pow2(k) else n
  den <- if (k >= 0) d else d * .pow2(-k)
  qq <- gmp::divq.bigz(num, den); r <- num - qq * den
  .mp_round(2 * qq + (r != 0), -k - 1, p)
}
.mp_from_dbl <- function(x) {
  if (x == 0) return(.mp_zero)
  ex <- floor(log2(abs(x))) - 52
  m <- round(x / 2^ex)
  while (abs(m) >= 2^53) { ex <- ex + 1; m <- round(x / 2^ex) }
  .mp_round(gmp::as.bigz(m), ex, 53)
}
.mp_to_q <- function(a) if (a$e >= 0) gmp::as.bigq(a$m * .pow2(a$e)) else gmp::as.bigq(a$m, .pow2(-a$e))
.mp_to_dbl <- function(a) {
  if (a$m == 0) return(0)
  b <- .bitlen(a$m)
  if (b > 60) { s <- b - 60; return(as.numeric(gmp::divq.bigz(a$m, .pow2(s))) * 2^(a$e + s)) }
  as.numeric(a$m) * 2^a$e
}
.mp_neg <- function(a) list(m = -a$m, e = a$e)
.mp_sign <- function(a) if (a$m > 0) 1 else if (a$m < 0) -1 else 0
.mp_mag <- function(a) if (a$m == 0) -Inf else .bitlen(a$m) + a$e

.mp_add <- function(a, b, p) {
  if (a$m == 0) return(.mp_round(b$m, b$e, p))
  if (b$m == 0) return(.mp_round(a$m, a$e, p))
  if (a$e < b$e) { t <- a; a <- b; b <- t }
  d <- a$e - b$e
  # b lies far below the rounding position of a: it only decides rounding
  if (.mp_mag(a) - .mp_mag(b) > p + 4) {
    k <- p + 4 + max(0, p - .bitlen(a$m))
    return(.mp_round(a$m * .pow2(k) + sign(as.numeric(b$m)), a$e - k, p))
  }
  .mp_round(a$m * .pow2(d) + b$m, b$e, p)
}
.mp_sub <- function(a, b, p) .mp_add(a, .mp_neg(b), p)
.mp_mul <- function(a, b, p) .mp_round(a$m * b$m, a$e + b$e, p)
.mp_div <- function(a, b, p) {
  if (b$m == 0) stop("division by zero in multiprecision arithmetic", call. = FALSE)
  if (a$m == 0) return(.mp_zero)
  num <- a$m; den <- b$m
  if (den < 0) { num <- -num; den <- -den }
  k <- max(0, p + .bitlen(den) - .bitlen(num) + 4)
  num <- num * .pow2(k)
  q <- gmp::divq.bigz(num, den); r <- num - q * den
  .mp_round(2 * q + (r != 0), a$e - b$e - k - 1, p)
}
.mp_cmp <- function(a, b) .mp_sign(.mp_sub(a, b, max(.bitlen(a$m), .bitlen(b$m)) + abs(a$e - b$e) + 2))
.mp_ldexp <- function(a, k) list(m = a$m, e = a$e + k)
.mp_isqrt <- function(n) {
  if (n < 2) return(n)
  x <- .pow2(ceiling(.bitlen(n) / 2))
  repeat { y <- gmp::divq.bigz(x + gmp::divq.bigz(n, x), 2); if (y >= x) return(x); x <- y }
}
.mp_sqrt <- function(a, p) {
  if (a$m == 0) return(.mp_zero)
  if (a$m < 0) stop("square root of a negative number", call. = FALSE)
  k <- max(0, 2 * (p + 4) - .bitlen(a$m))
  if ((a$e - k) %% 2 != 0) k <- k + 1
  n <- a$m * .pow2(k)
  r <- .mp_isqrt(n)
  .mp_round(2 * r + (r * r != n), (a$e - k) / 2 - 1, p)
}
.mp_floor <- function(a) if (a$e >= 0) a$m * .pow2(a$e) else gmp::divq.bigz(a$m, .pow2(-a$e))

# division rounding toward zero, so series terms of either sign reach 0
.tdiv <- function(a, b) {
  a <- gmp::as.bigz(a); b <- gmp::as.bigz(b)
  q <- gmp::divq.bigz(abs(a), abs(b))
  if ((a < 0) != (b < 0)) -q else q
}

# fixed point: bigz x with value x / 2^wp
.fxp <- function(a, wp) {
  s <- a$e + wp
  if (s >= 0) a$m * .pow2(s) else gmp::divq.bigz(a$m, .pow2(-s))
}
.fx_mp <- function(x, wp, p) .mp_round(x, -wp, p)

# constants by series, cached per precision
.mp_cache <- new.env()
.fx_atan_inv <- function(n, wp) {
  # atan(1/n) for integer n, by its Taylor series
  one <- .pow2(wp); n2 <- gmp::as.bigz(n)^2
  term <- gmp::divq.bigz(one, gmp::as.bigz(n)); s <- term; k <- 1
  repeat {
    term <- gmp::divq.bigz(term, n2)
    if (term == 0) break
    t <- gmp::divq.bigz(term, 2 * k + 1)
    s <- if (k %% 2 == 1) s - t else s + t
    k <- k + 1
  }
  s
}
.fx_pi <- function(wp) {
  key <- paste0("pi", wp)
  if (!is.null(v <- .mp_cache[[key]])) return(v)
  w <- wp + 10
  v <- gmp::divq.bigz(16 * .fx_atan_inv(5, w) - 4 * .fx_atan_inv(239, w), .pow2(10))
  assign(key, v, envir = .mp_cache)
  v
}
.fx_ln2 <- function(wp) {
  key <- paste0("ln2", wp)
  if (!is.null(v <- .mp_cache[[key]])) return(v)
  # log(2) = 2*atanh(1/3) = 2 * sum 1/((2k+1) 3^(2k+1))
  w <- wp + 10; one <- .pow2(w)
  term <- gmp::divq.bigz(one, 3); s <- term; k <- 1
  repeat {
    term <- gmp::divq.bigz(term, 9)
    if (term == 0) break
    s <- s + gmp::divq.bigz(term, 2 * k + 1); k <- k + 1
  }
  v <- gmp::divq.bigz(2 * s, .pow2(10))
  assign(key, v, envir = .mp_cache)
  v
}
.mp_pi <- function(p) .fx_mp(.fx_pi(p + 10), p + 10, p)
.mp_ln2 <- function(p) .fx_mp(.fx_ln2(p + 10), p + 10, p)

# exp(x) for real x
.mp_exp <- function(a, p) {
  if (a$m == 0) return(.mpf(1, 0))
  mag <- .mp_mag(a)
  wp <- p + 20 + max(0, mag)
  # x = n*log(2) + r, |r| <= log(2)/2
  ln2 <- .fx_mp(.fx_ln2(wp), wp, wp)
  n <- .mp_floor(.mp_add(.mp_div(a, ln2, wp), .mpf(1, -1), wp))
  r <- .mp_sub(a, .mp_mul(.mp_round(n, 0, wp), ln2, wp), wp)
  # r / 2^k, Taylor series, then k squarings
  k <- max(0, floor(sqrt(wp) / 2))
  w2 <- wp + k + 10
  x <- .tdiv(.fxp(r, w2), .pow2(k))
  one <- .pow2(w2)
  s <- one; term <- one; j <- 1
  repeat {
    term <- .tdiv(.tdiv(term * x, one), gmp::as.bigz(j))
    if (term == 0) break
    s <- s + term; j <- j + 1
  }
  for (i in seq_len(k)) s <- .tdiv(s * s, one)
  .mp_ldexp(.fx_mp(s, w2, p + 4), as.numeric(n)) |> (\(z) .mp_round(z$m, z$e, p))()
}

# log(x) for real x > 0: x = y * 2^n with y near 1, log(y) = 2*atanh((y-1)/(y+1))
.mp_log <- function(a, p, extra = 0) {
  if (a$m <= 0) stop("log of a nonpositive number in real arithmetic", call. = FALSE)
  wp <- p + 20 + extra
  b <- .bitlen(a$m)
  n <- b + a$e
  y <- list(m = a$m, e = -b)                  # y in [1/2, 1)
  if (.mp_cmp(y, .mp_sqrt(.mpf(1, -1), wp)) < 0) { y <- .mp_ldexp(y, 1); n <- n - 1 }
  z <- .mp_div(.mp_sub(y, .mpf(1, 0), wp), .mp_add(y, .mpf(1, 0), wp), wp)
  w2 <- wp + 10; one <- .pow2(w2)
  zf <- .fxp(z, w2); z2 <- .tdiv(zf * zf, one)
  s <- zf; term <- zf; k <- 1
  repeat {
    term <- .tdiv(term * z2, one)
    if (term == 0) break
    s <- s + .tdiv(term, 2 * k + 1); k <- k + 1
  }
  ly <- .fx_mp(2 * s, w2, wp)
  r <- .mp_add(ly, .mp_mul(.mp_round(gmp::as.bigz(n), 0, wp), .fx_mp(.fx_ln2(wp), wp, wp), wp), p)
  # log of a number near 1: redo with the bits lost to cancellation
  if (extra == 0 && r$m != 0 && .mp_mag(r) < -10) return(.mp_log(a, p, extra = 10 - .mp_mag(r)))
  r
}

# sin and cos of real x, by reduction modulo pi/2
.mp_sincos <- function(a, p, extra = 0) {
  if (a$m == 0) return(list(sin = .mp_zero, cos = .mpf(1, 0)))
  # |x| below the working precision: sin(x) = x and cos(x) = 1 to p bits
  if (.mp_mag(a) < -(p + 10)) return(list(sin = .mp_round(a$m, a$e, p), cos = .mpf(1, 0)))
  mag <- max(0, .mp_mag(a))
  wp <- p + 20 + mag + extra + max(0, -.mp_mag(a))
  pi2 <- .mp_ldexp(.fx_mp(.fx_pi(wp), wp, wp), -1)
  n <- .mp_floor(.mp_add(.mp_div(a, pi2, wp), .mpf(1, -1), wp))
  r <- .mp_sub(a, .mp_mul(.mp_round(n, 0, wp), pi2, wp), wp)
  # x near a multiple of pi/2: redo with the bits lost to cancellation
  if (extra == 0 && n != 0 && r$m != 0 && .mp_mag(r) < -10) return(.mp_sincos(a, p, extra = 10 - .mp_mag(r)))
  w2 <- wp + 10; one <- .pow2(w2)
  x <- .fxp(r, w2); x2 <- .tdiv(x * x, one)
  s <- x; c <- one; ts <- x; tc <- one; k <- 1
  repeat {
    ts <- -.tdiv(.tdiv(ts * x2, one), gmp::as.bigz((2 * k) * (2 * k + 1)))
    tc <- -.tdiv(.tdiv(tc * x2, one), gmp::as.bigz((2 * k - 1) * (2 * k)))
    if (ts == 0 && tc == 0) break
    s <- s + ts; c <- c + tc; k <- k + 1
  }
  S <- .fx_mp(s, w2, p); C <- .fx_mp(c, w2, p)
  switch(as.integer(gmp::mod.bigz(n, 4)) + 1L,
         list(sin = S, cos = C), list(sin = C, cos = .mp_neg(S)),
         list(sin = .mp_neg(S), cos = .mp_neg(C)), list(sin = .mp_neg(C), cos = S))
}

# atan(x) for real x
.mp_atan <- function(a, p) {
  if (a$m == 0) return(.mp_zero)
  wp <- p + 20 + max(0, -.mp_mag(a))
  if (.mp_cmp(list(m = abs(a$m), e = a$e), .mpf(1, 0)) > 0) {
    half_pi <- .mp_ldexp(.mp_pi(wp), -1)
    r <- .mp_sub(half_pi, .mp_atan(.mp_div(.mpf(1, 0), list(m = abs(a$m), e = a$e), wp), wp), p)
    return(if (a$m < 0) .mp_neg(r) else r)
  }
  # halve the argument: atan(x) = 2*atan(x / (1 + sqrt(1 + x^2)))
  x <- a; h <- 0
  while (.mp_mag(x) > -8) {
    x <- .mp_div(x, .mp_add(.mpf(1, 0), .mp_sqrt(.mp_add(.mpf(1, 0), .mp_mul(x, x, wp), wp), wp), wp), wp)
    h <- h + 1
  }
  w2 <- wp + 10; one <- .pow2(w2)
  xf <- .fxp(x, w2); x2 <- .tdiv(xf * xf, one)
  s <- xf; term <- xf; k <- 1
  repeat {
    term <- -.tdiv(term * x2, one)
    if (term == 0) break
    s <- s + .tdiv(term, 2 * k + 1); k <- k + 1
  }
  .mp_ldexp(.fx_mp(s, w2, p), h)
}

.mp_pow_int <- function(a, n, p) {
  n <- gmp::as.bigz(n)
  if (n == 0) return(.mpf(1, 0))
  wp <- p + 10 + .bitlen(n)
  r <- .mpf(1, 0); b <- a; k <- abs(n)
  while (k > 0) {
    if (gmp::mod.bigz(k, 2) == 1) r <- .mp_mul(r, b, wp)
    k <- gmp::divq.bigz(k, 2)
    if (k > 0) b <- .mp_mul(b, b, wp)
  }
  if (n < 0) r <- .mp_div(.mpf(1, 0), r, wp)
  .mp_round(r$m, r$e, p)
}

# decimal string with dps significant digits, after mpmath's to_str
.mp_to_str <- function(a, dps, strip_zeros = TRUE, min_fixed = NULL, max_fixed = NULL) {
  # mpmath's to_str: zero is always 0.0
  if (a$m == 0) return("0.0")
  if (is.null(min_fixed)) min_fixed <- min(-(dps %/% 3), -5)
  if (is.null(max_fixed)) max_fixed <- dps
  neg <- a$m < 0
  q <- abs(.mp_to_q(a))
  # exponent estimate, corrected after rounding to dps digits
  ex <- floor((.mp_mag(a) - 1) * log10(2))
  repeat {
    sc <- q * gmp::as.bigq(10)^(dps - 1 - ex)
    n <- gmp::numerator(sc); d <- gmp::denominator(sc)
    digits <- gmp::divq.bigz(2 * n + d, 2 * d)
    nd <- nchar(as.character(digits))
    if (nd > dps) { ex <- ex + 1; next }
    if (nd < dps) { ex <- ex - 1; next }
    break
  }
  ds <- as.character(digits)
  fixed <- min_fixed < ex && ex < max_fixed
  if (fixed) {
    if (ex >= 0) {
      ip <- substr(ds, 1L, ex + 1L); fp <- substring(ds, ex + 2L)
      if (nchar(ip) < ex + 1L) ip <- paste0(ip, strrep("0", ex + 1L - nchar(ip)))
    } else {
      ip <- "0"; fp <- paste0(strrep("0", -ex - 1L), ds)
    }
    if (strip_zeros) fp <- sub("0+$", "", fp)
    if (!nzchar(fp)) fp <- "0"
    s <- paste0(ip, ".", fp)
  } else {
    fp <- substring(ds, 2L)
    if (strip_zeros) fp <- sub("0+$", "", fp)
    if (!nzchar(fp)) fp <- "0"
    s <- paste0(substr(ds, 1L, 1L), ".", fp, "e", if (ex >= 0) "+" else "-", abs(ex))
  }
  if (neg) paste0("-", s) else s
}

.dps_to_prec <- function(dps) max(1, round((dps + 1) * 3.3219280948873626))
.prec_to_dps <- function(prec) max(1, round(prec / 3.3219280948873626) - 1)

# ---------------------------------------------------------------------------
# Complex values and expression evaluation

.mpc <- function(re, im = .mp_zero) list(re = re, im = im)
.mpc_real <- function(z) z$im$m == 0
.mpc_add <- function(a, b, p) .mpc(.mp_add(a$re, b$re, p), .mp_add(a$im, b$im, p))
.mpc_mul <- function(a, b, p) {
  if (.mpc_real(a) && .mpc_real(b)) return(.mpc(.mp_mul(a$re, b$re, p)))
  .mpc(.mp_sub(.mp_mul(a$re, b$re, p + 10), .mp_mul(a$im, b$im, p + 10), p),
       .mp_add(.mp_mul(a$re, b$im, p + 10), .mp_mul(a$im, b$re, p + 10), p))
}
.mpc_div <- function(a, b, p) {
  if (.mpc_real(b)) return(.mpc(.mp_div(a$re, b$re, p), .mp_div(a$im, b$re, p)))
  w <- p + 10
  n <- .mp_add(.mp_mul(b$re, b$re, w), .mp_mul(b$im, b$im, w), w)
  .mpc(.mp_div(.mp_add(.mp_mul(a$re, b$re, w), .mp_mul(a$im, b$im, w), w), n, p),
       .mp_div(.mp_sub(.mp_mul(a$im, b$re, w), .mp_mul(a$re, b$im, w), w), n, p))
}
.mp_atan2 <- function(y, x, p) {
  if (x$m == 0) {
    if (y$m == 0) return(.mp_zero)
    h <- .mp_ldexp(.mp_pi(p), -1)
    return(if (y$m > 0) h else .mp_neg(h))
  }
  t <- .mp_atan(.mp_div(y, x, p + 10), p + 10)
  if (x$m > 0) return(.mp_round(t$m, t$e, p))
  pi <- .mp_pi(p + 10)
  if (y$m >= 0) .mp_add(t, pi, p) else .mp_sub(t, pi, p)
}
.mpc_abs <- function(z, p) .mp_sqrt(.mp_add(.mp_mul(z$re, z$re, p + 10), .mp_mul(z$im, z$im, p + 10), p + 10), p)
.mpc_exp <- function(z, p) {
  r <- .mp_exp(z$re, p + 10)
  if (.mpc_real(z)) return(.mpc(.mp_round(r$m, r$e, p)))
  sc <- .mp_sincos(z$im, p + 10)
  .mpc(.mp_mul(r, sc$cos, p), .mp_mul(r, sc$sin, p))
}
.mpc_log <- function(z, p) {
  if (.mpc_real(z) && z$re$m > 0) return(.mpc(.mp_log(z$re, p)))
  if (z$re$m == 0 && z$im$m == 0) stop("log of zero", call. = FALSE)
  .mpc(.mp_log(.mpc_abs(z, p + 10), p), .mp_atan2(z$im, z$re, p))
}
.mpc_sqrt <- function(z, p) {
  if (.mpc_real(z)) {
    if (z$re$m >= 0) return(.mpc(.mp_sqrt(z$re, p)))
    return(.mpc(.mp_zero, .mp_sqrt(.mp_neg(z$re), p)))
  }
  w <- p + 10
  r <- .mpc_abs(z, w)
  a <- .mp_sqrt(.mp_ldexp(.mp_add(r, z$re, w), -1), w)
  b <- .mp_sqrt(.mp_ldexp(.mp_sub(r, z$re, w), -1), w)
  if (z$im$m < 0) b <- .mp_neg(b)
  .mpc(.mp_round(a$m, a$e, p), .mp_round(b$m, b$e, p))
}
.mpc_pow <- function(b, x, p, xq = NULL) {
  # a real float exponent that is an integer or 1/2, as mpc_pow_mpf
  if (is.null(xq) && x$im$m == 0) {
    if (x$re$m == 0) xq <- gmp::as.bigq(0)
    else if (x$re$e >= 0) xq <- gmp::as.bigq(x$re$m * .pow2(x$re$e))
    else if (x$re$e == -1) xq <- gmp::as.bigq(x$re$m, 2)
  }
  # half-integers: sqrt(b)^m
  if (!is.null(xq) && gmp::denominator(xq) == 2 && gmp::numerator(xq) != 1)
    return(.mpc_pow(.mpc_sqrt(b, p + 10 + .bitlen(gmp::numerator(xq))), .mpc(.mpf(gmp::numerator(xq), 0)), p, gmp::as.bigq(gmp::numerator(xq))))
  # an integer exponent: repeated squaring; otherwise exp(x*log(b))
  if (!is.null(xq) && gmp::denominator(xq) == 1) {
    n <- gmp::numerator(xq)
    if (.mpc_real(b)) return(.mpc(.mp_pow_int(b$re, n, p)))
    w <- p + 10 + .bitlen(n); r <- .mpc(.mpf(1, 0)); bb <- b; k <- abs(n)
    while (k > 0) { if (gmp::mod.bigz(k, 2) == 1) r <- .mpc_mul(r, bb, w); k <- gmp::divq.bigz(k, 2); if (k > 0) bb <- .mpc_mul(bb, bb, w) }
    return(if (n < 0) .mpc_div(.mpc(.mpf(1, 0)), r, p) else r)
  }
  if (!is.null(xq) && gmp::denominator(xq) == 2 && gmp::numerator(xq) == 1) return(.mpc_sqrt(b, p))
  if (b$re$m == 0 && b$im$m == 0) return(.mpc(.mp_zero))
  w <- p + 20 + max(0, .mp_mag(x$re))
  .mpc_exp(.mpc_mul(x, .mpc_log(b, w), w), p)
}
.mpc_sin <- function(z, p) {
  sc <- .mp_sincos(z$re, p + 10)
  if (.mpc_real(z)) return(.mpc(.mp_round(sc$sin$m, sc$sin$e, p)))
  ch <- .mpc_cosh_sinh(z$im, p + 10)
  .mpc(.mp_mul(sc$sin, ch$cosh, p), .mp_mul(sc$cos, ch$sinh, p))
}
.mpc_cos <- function(z, p) {
  sc <- .mp_sincos(z$re, p + 10)
  if (.mpc_real(z)) return(.mpc(.mp_round(sc$cos$m, sc$cos$e, p)))
  ch <- .mpc_cosh_sinh(z$im, p + 10)
  .mpc(.mp_mul(sc$cos, ch$cosh, p), .mp_neg(.mp_mul(sc$sin, ch$sinh, p)))
}
.mpc_cosh_sinh <- function(x, p) {
  e1 <- .mp_exp(x, p + 10); e2 <- .mp_div(.mpf(1, 0), e1, p + 10)
  list(cosh = .mp_ldexp(.mp_add(e1, e2, p), -1), sinh = .mp_ldexp(.mp_sub(e1, e2, p), -1))
}
.mpc_I <- .mpc(.mp_zero, .mpf(1, 0))
.mpc_neg <- function(z) .mpc(.mp_neg(z$re), .mp_neg(z$im))
.mpc_one <- .mpc(.mpf(1, 0))
.mpc_sub <- function(a, b, p) .mpc_add(a, .mpc_neg(b), p)

# Elementary functions on complex arguments (log formulas off the real axis).
.mpc_fun <- function(nm, z, p) {
  mz <- suppressWarnings(max(.mp_mag(z$re), .mp_mag(z$im)))
  if (is.finite(mz) && mz < -4 && !(nm %in% c("Abs", "re", "im", "conjugate", "sign", "floor", "ceiling", "arg"))) {
    r <- .mpc_fun_raw(nm, z, p + 10 - mz)
    if (is.null(r)) return(NULL)
    return(.mpc(.mp_round(r$re$m, r$re$e, p), .mp_round(r$im$m, r$im$e, p)))
  }
  .mpc_fun_raw(nm, z, p)
}
.mpc_fun_raw <- function(nm, z, p) {
  w <- p + 20
  one <- .mpc_one; I <- .mpc_I
  real <- .mpc_real(z)
  inside <- function(lo, hi) real && .mp_cmp(z$re, lo) >= 0 && .mp_cmp(z$re, hi) <= 0
  m1 <- .mpf(-1, 0); p1 <- .mpf(1, 0)
  switch(EXPR = nm,
    exp = .mpc_exp(z, p), log = .mpc_log(z, p), sqrt = .mpc_sqrt(z, p),
    sin = .mpc_sin(z, p), cos = .mpc_cos(z, p),
    tan = .mpc_div(.mpc_sin(z, w), .mpc_cos(z, w), p), cot = .mpc_div(.mpc_cos(z, w), .mpc_sin(z, w), p),
    sec = .mpc_div(one, .mpc_cos(z, w), p), csc = .mpc_div(one, .mpc_sin(z, w), p),
    sinh = { e1 <- .mpc_exp(z, w); e2 <- .mpc_div(one, e1, w); r <- .mpc_sub(e1, e2, w); .mpc(.mp_ldexp(r$re, -1), .mp_ldexp(r$im, -1)) },
    cosh = { e1 <- .mpc_exp(z, w); e2 <- .mpc_div(one, e1, w); r <- .mpc_add(e1, e2, w); .mpc(.mp_ldexp(r$re, -1), .mp_ldexp(r$im, -1)) },
    tanh = { e2 <- .mpc_exp(.mpc(.mp_ldexp(z$re, 1), .mp_ldexp(z$im, 1)), w); .mpc_div(.mpc_sub(e2, one, w), .mpc_add(e2, one, w), p) },
    coth = { e2 <- .mpc_exp(.mpc(.mp_ldexp(z$re, 1), .mp_ldexp(z$im, 1)), w); .mpc_div(.mpc_add(e2, one, w), .mpc_sub(e2, one, w), p) },
    sech = .mpc_div(one, .mpc_fun("cosh", z, w), p), csch = .mpc_div(one, .mpc_fun("sinh", z, w), p),
    atan = if (real) .mpc(.mp_atan(z$re, p)) else {
      # atan(z) = I/2*(log(1 - I*z) - log(1 + I*z))
      iz <- .mpc_mul(I, z, w)
      d <- .mpc_sub(.mpc_log(.mpc_sub(one, iz, w), w), .mpc_log(.mpc_add(one, iz, w), w), w)
      r <- .mpc_mul(I, d, w); .mpc(.mp_ldexp(r$re, -1), .mp_ldexp(r$im, -1))
    },
    acot = if (real) {
      if (z$re$m == 0) .mpc(.mp_ldexp(.mp_pi(p), -1)) else .mpc(.mp_atan(.mp_div(p1, z$re, w), p))
    } else .mpc_fun("atan", .mpc_div(one, z, w), p),
    # mpmath on the real axis beyond [-1, 1]: asin(x) = sign(x)*(pi/2 - I*acosh(|x|))
    asin = if (real && !inside(m1, p1)) {
      h <- .mp_ldexp(.mp_pi(p), -1)
      ac <- .mpc_fun("acosh", .mpc(list(m = abs(z$re$m), e = z$re$e)), p)$re
      if (z$re$m > 0) .mpc(h, .mp_neg(ac)) else .mpc(.mp_neg(h), ac)
    } else if (inside(m1, p1)) {
      if (.mp_cmp(list(m = abs(z$re$m), e = z$re$e), p1) == 0)
        .mpc(.mp_ldexp(if (z$re$m > 0) .mp_pi(p) else .mp_neg(.mp_pi(p)), -1))
      else .mpc(.mp_atan(.mp_div(z$re, .mp_sqrt(.mp_sub(p1, .mp_mul(z$re, z$re, w), w), w), w), p))
    } else {
      # -I*log(I*z + sqrt(1 - z^2))
      r <- .mpc_log(.mpc_add(.mpc_mul(I, z, w), .mpc_sqrt(.mpc_sub(one, .mpc_mul(z, z, w), w), w), w), w)
      .mpc(.mp_round(r$im$m, r$im$e, p), .mp_neg(.mp_round(r$re$m, r$re$e, p)))
    },
    # mpmath on the real axis: acos(x) = I*acosh(x) for x > 1 and
    # pi - I*acosh(-x) for x < -1
    acos = if (real && !inside(m1, p1)) {
      ac <- .mpc_fun("acosh", .mpc(list(m = abs(z$re$m), e = z$re$e)), p)$re
      if (z$re$m > 0) .mpc(.mp_zero, ac) else .mpc(.mp_pi(p), .mp_neg(ac))
    } else .mpc_sub(.mpc(.mp_ldexp(.mp_pi(w), -1)), .mpc_fun("asin", z, w), p),
    asec = .mpc_fun("acos", .mpc_div(one, z, w), p), acsc = .mpc_fun("asin", .mpc_div(one, z, w), p),
    asinh = .mpc_log(.mpc_add(z, .mpc_sqrt(.mpc_add(.mpc_mul(z, z, w), one, w), w), w), p),
    acosh = .mpc_log(.mpc_add(z, .mpc_mul(.mpc_sqrt(.mpc_add(z, one, w), w), .mpc_sqrt(.mpc_sub(z, one, w), w), w), w), p),
    atanh = { r <- .mpc_sub(.mpc_log(.mpc_add(one, z, w), w), .mpc_log(.mpc_sub(one, z, w), w), w); .mpc(.mp_ldexp(r$re, -1), .mp_ldexp(r$im, -1)) },
    acoth = .mpc_fun("atanh", .mpc_div(one, z, w), p),
    Abs = .mpc(.mpc_abs(z, p)), re = .mpc(z$re), im = .mpc(z$im),
    arg = .mpc(.mp_atan2(z$im, z$re, p)), conjugate = .mpc(z$re, .mp_neg(z$im)),
    sign = if (real) .mpc(.mpf(.mp_sign(z$re), 0)) else .mpc_div(z, .mpc(.mpc_abs(z, w)), p),
    floor = if (real) .mpc(.mp_round(.mp_floor(z$re), 0, p)) else NULL,
    ceiling = if (real) .mpc(.mp_neg(.mp_round(.mp_floor(.mp_neg(z$re)), 0, p))) else NULL,
    gamma = if (real) .mpc(.mp_gamma(z$re, p)) else NULL,
    loggamma = if (real && z$re$m > 0) .mpc(.mp_loggamma(z$re, p)) else NULL,
    factorial = if (real) .mpc(.mp_gamma(.mp_add(z$re, p1, w), p)) else NULL,
    erf = if (real) .mpc(.mp_erf(z$re, p)) else NULL,
    erfc = if (real) .mpc(.mp_erfc(z$re, p)) else NULL,
    zeta = if (real) .mpc(.mp_zeta(z$re, p)) else NULL,
    NULL)
}

.mp_const <- function(nm, p) {
  switch(EXPR = nm,
    pi = .mp_pi(p), E = .mp_exp(.mpf(1, 0), p),
    GoldenRatio = .mp_ldexp(.mp_add(.mpf(1, 0), .mp_sqrt(.mpf(5, 0), p + 10), p), -1),
    EulerGamma = .mp_euler(p), Catalan = .mp_catalan(p),
    NULL)
}

# Evaluate e at working precision p; NULL when a part has no
# arbitrary precision evaluator.
.mpeval <- function(e, p, env = list()) {
  t <- e$type
  r <- switch(t,
    Integer = , Rational = .mpc(.mp_from_q(e$value, p)),
    Float = .mpc(.float_mp(e)),
    AlgebraicNumber = .mpeval(.an_value(e), p, env),
    Symbol = , Dummy = { v <- env[[e$name]]; if (is.null(v)) NULL else .mpeval(as_sym(v), p, env) },
    Constant = if (e$name == "I") .mpc_I else { c <- .mp_const(e$name, p); if (is.null(c)) NULL else .mpc(c) },
    Add = {
      # evalf_add: a sum that cancels below the target accuracy is redone with
      # the lost bits added to the working precision, up to a cap; agreeing
      # results at two precisions cannot certify it, since both may lose the
      # same bits.
      cap <- p + 2000
      inexact <- any(vapply(e$args, function(a) a$type != "Integer", NA))
      wp <- p
      repeat {
        vs <- vector("list", length(e$args))
        for (i in seq_along(e$args)) {
          v <- .mpeval(e$args[[i]], wp, env)
          if (is.null(v)) return(NULL)
          vs[[i]] <- v
        }
        s <- .mpc(.mp_zero)
        for (v in vs) s <- .mpc_add(s, v, wp)
        lost <- 0
        for (part in c("re", "im")) {
          mx <- max(vapply(vs, function(v) .mp_mag(v[[part]]), 0))
          if (is.finite(mx)) {
            r <- .mp_mag(s[[part]])
            lost <- max(lost, if (is.finite(r)) mx - r else Inf)
          }
        }
        if (!inexact || wp - lost >= p || wp >= cap) break
        wp <- min(cap, max(wp + 32, if (is.finite(lost)) p + lost + 16 else cap))
      }
      .mpc(.mp_round(s$re$m, s$re$e, p), .mp_round(s$im$m, s$im$e, p))
    },
    Mul = {
      s <- .mpc_one
      for (a in e$args) { v <- .mpeval(a, p, env); if (is.null(v)) return(NULL); s <- .mpc_mul(s, v, p) }
      s
    },
    Pow = {
      b <- .mpeval(e$args[[1L]], p, env); x <- e$args[[2L]]
      if (is.null(b)) return(NULL)
      xq <- if (x$type %in% c("Integer", "Rational")) gmp::as.bigq(x$value) else NULL
      xv <- .mpeval(x, p, env); if (is.null(xv)) return(NULL)
      .mpc_pow(b, xv, p, xq)
    },
    Function = {
      if (length(e$args) != 1L) return(NULL)
      a <- .mpeval(e$args[[1L]], p, env)
      if (is.null(a)) return(NULL)
      if (e$name %in% c("floor", "ceiling", "frac") && .mpc_real(a)) {
        r <- .mp_near_integer(e$name, e$args[[1L]], a$re, p, env)
        if (!is.null(r)) return(r)
      }
      .mpc_fun(e$name, a, p)
    },
    NULL)
  r
}

# floor, ceiling and frac of an argument within rounding error of an
# integer k, after get_integer_part: the sign of arg - k at higher precision
# decides; a difference that stays zero means arg = k.
.mp_near_integer <- function(nm, arg, a, p, env) {
  k <- .mp_floor(.mp_add(a, .mpf(1, -1), p + 10))
  kq <- gmp::as.bigq(k)
  d0 <- .mp_sub(a, .mp_round(k, 0, p + 10), p + 10)
  if (d0$m != 0 && .mp_mag(d0) >= max(0, .mp_mag(a)) - p + 8) return(NULL)
  dz <- tryCatch(.mpeval_adaptive(.sub(arg, make_rational(kq)), 20L, env, maxprec = max(2000, 8 * p)),
                 error = function(err) NULL)
  s <- if (is.null(dz) || dz$re$m == 0) 0 else if (dz$re$m > 0) 1 else -1
  kk <- .mp_round(k, 0, p)
  switch(EXPR = nm,
    floor = .mpc(if (s < 0) .mp_sub(kk, .mpf(1, 0), p) else kk),
    ceiling = .mpc(if (s > 0) .mp_add(kk, .mpf(1, 0), p) else kk),
    frac = if (s == 0) .mpc(.mp_zero) else if (s > 0) .mpc(.mp_round(dz$re$m, dz$re$e, p)) else .mpc(.mp_add(.mpf(1, 0), dz$re, p)))
}

# Adaptive evaluation to prec bits: repeat at higher working precision until
# two evaluations agree to prec bits, as mpmath's autoprec.
.mpeval_adaptive <- function(e, prec, env = list(), maxprec = NULL) {
  if (is.null(maxprec)) maxprec <- max(2000, 8 * prec)
  wp <- prec + 30
  prev <- .mpeval(e, wp, env)
  if (is.null(prev)) return(NULL)
  repeat {
    wp2 <- min(maxprec, 2 * wp + 30)
    cur <- .mpeval(e, wp2, env)
    if (is.null(cur)) return(NULL)
    ok <- function(a, b) {
      if (a$m == 0 && b$m == 0) return(TRUE)
      d <- .mp_sub(a, b, wp2 + 10)
      d$m == 0 || .mp_mag(d) < max(.mp_mag(a), .mp_mag(b)) - prec - 4
    }
    conv <- ok(prev$re, cur$re) && ok(prev$im, cur$im)
    # a zero may come from information lost to rounding: keep going
    if (conv && cur$re$m == 0 && cur$im$m == 0 && !is_exact_number(e) && wp2 < maxprec) conv <- FALSE
    if (conv || wp2 >= maxprec) {
      r <- .mpc(.mp_round(cur$re$m, cur$re$e, prec), .mp_round(cur$im$m, cur$im$e, prec))
      attr(r, "converged") <- conv
      return(r)
    }
    prev <- cur; wp <- wp2
  }
}

# ---------------------------------------------------------------------------
# Special functions and constants

.mp_bernoulli <- function(n) {
  key <- paste0("B", n)
  if (!is.null(v <- .mp_cache[[key]])) return(v)
  v <- gmp::as.bigq(bernoulli(.int(n))$value)
  assign(key, v, envir = .mp_cache)
  v
}

# loggamma(x) for real x > 0: shift to z >= wp/4, Stirling's series, undo
.mp_loggamma <- function(a, p) {
  wp <- p + 30
  x <- a
  z0 <- max(10, ceiling(wp / 4))
  shift <- max(0, ceiling(z0 - .mp_to_dbl(x)))
  z <- .mp_add(x, .mpf(shift, 0), wp)
  lz <- .mp_log(z, wp)
  s <- .mp_sub(.mp_mul(.mp_sub(z, .mpf(1, -1), wp), lz, wp), z, wp)
  s <- .mp_add(s, .mp_ldexp(.mp_log(.mp_ldexp(.mp_pi(wp), 1), wp), -1), wp)
  zinv <- .mp_div(.mpf(1, 0), z, wp); z2 <- .mp_mul(zinv, zinv, wp); zp <- zinv
  for (k in 1:2000) {
    t <- .mp_mul(.mp_from_q(.mp_bernoulli(2 * k) / (2 * k * (2 * k - 1)), wp), zp, wp)
    if (t$m == 0 || .mp_mag(t) < .mp_mag(s) - wp) break
    s <- .mp_add(s, t, wp)
    zp <- .mp_mul(zp, z2, wp)
  }
  if (shift > 0) {
    pr <- .mpf(1, 0)
    for (j in 0:(shift - 1)) pr <- .mp_mul(pr, .mp_add(x, .mpf(j, 0), wp), wp)
    s <- .mp_sub(s, .mp_log(list(m = abs(pr$m), e = pr$e), wp), wp)
  }
  .mp_round(s$m, s$e, p)
}
.mp_gamma <- function(a, p) {
  wp <- p + 30
  if (a$e >= 0 || .mp_cmp(.mp_round(.mp_floor(a), 0, wp), a) == 0) {
    n <- .mp_floor(a)
    if (.mp_cmp(.mp_round(n, 0, wp), a) == 0) {
      if (n <= 0) stop("gamma has a pole at a nonpositive integer", call. = FALSE)
      if (n < 3000) return(.mp_round(gmp::factorialZ(as.numeric(n) - 1), 0, p))
    }
  }
  if (.mp_cmp(a, .mpf(1, -1)) < 0) {
    # reflection: gamma(x) = pi / (sin(pi*x) * gamma(1 - x))
    s <- .mp_sincos(.mp_mul(.mp_pi(wp), a, wp), wp)$sin
    g <- .mp_gamma(.mp_sub(.mpf(1, 0), a, wp), wp)
    return(.mp_div(.mp_pi(wp), .mp_mul(s, g, wp), p))
  }
  .mp_exp(.mp_loggamma(a, wp), p)
}

# erfc(x): 1 - erf(x) with the bits lost to cancellation for moderate x,
# the asymptotic series for large x (terms decrease up to n = x^2)
.mp_erfc <- function(a, p) {
  if (a$m <= 0) return(.mp_sub(.mpf(1, 0), .mp_erf(a, p + 10), p))
  x <- .mp_to_dbl(a)
  lost <- ceiling(x * x * 1.4427) + 10
  if (lost < 4 * p + 2000) return(.mp_sub(.mpf(1, 0), .mp_erf(a, p + lost), p))
  wp <- p + 30
  x2 <- .mp_mul(a, a, wp)
  inv2x2 <- .mp_div(.mpf(1, 0), .mp_ldexp(x2, 1), wp)
  s <- .mpf(1, 0); term <- .mpf(1, 0); n <- 1
  repeat {
    term <- .mp_neg(.mp_mul(.mp_mul(term, .mpf(2 * n - 1, 0), wp), inv2x2, wp))
    if (term$m == 0 || .mp_mag(term) < -wp) break
    s <- .mp_add(s, term, wp); n <- n + 1
    if (n > x * x) break
  }
  f <- .mp_div(.mp_exp(.mp_neg(x2), wp), .mp_mul(a, .mp_sqrt(.mp_pi(wp), wp), wp), wp)
  .mp_mul(f, s, p)
}

.mp_erf <- function(a, p) {
  if (a$m == 0) return(.mp_zero)
  x2 <- .mp_to_dbl(a)^2
  if (x2 > 0.7 * p + 10) return(.mpf(if (a$m > 0) 1 else -1, 0))
  wp <- p + 30 + ceiling(x2 * 1.4427)
  w2 <- wp + 10; one <- .pow2(w2)
  xf <- .fxp(a, w2); x2f <- .tdiv(xf * xf, one)
  term <- xf; s <- xf; n <- 1
  repeat {
    term <- .tdiv(-.tdiv(term * x2f, one), gmp::as.bigz(n))
    if (term == 0) break
    s <- s + .tdiv(term, gmp::as.bigz(2 * n + 1)); n <- n + 1
  }
  r <- .fx_mp(2 * s, w2, wp)
  .mp_div(r, .mp_sqrt(.mp_pi(wp), wp), p)
}

# zeta(s) for real s != 1: Borwein's algorithm for s >= 1/2, the
# functional equation below
.mp_zeta <- function(a, p) {
  wp <- p + 30
  if (a$m == 0) return(.mpf(-1, -1))
  if (.mp_cmp(a, .mpf(1, 0)) == 0) stop("zeta has a pole at 1", call. = FALSE)
  if (.mp_cmp(a, .mpf(1, -1)) < 0) {
    one_s <- .mp_sub(.mpf(1, 0), a, wp)
    if (a$e >= 0 && a$m < 0 && gmp::mod.bigz(.mp_floor(a), 2) == 0) return(.mp_zero)
    f <- .mp_mul(.mp_exp(.mp_mul(a, .mp_log(.mpf(2, 0), wp), wp), wp),
                 .mp_exp(.mp_mul(.mp_sub(a, .mpf(1, 0), wp), .mp_log(.mp_pi(wp), wp), wp), wp), wp)
    f <- .mp_mul(f, .mp_sincos(.mp_ldexp(.mp_mul(.mp_pi(wp), a, wp), -1), wp)$sin, wp)
    f <- .mp_mul(f, .mp_gamma(one_s, wp), wp)
    return(.mp_mul(f, .mp_zeta(one_s, wp), p))
  }
  n <- ceiling(wp * 0.3934) + 10
  # d_k = n * sum_{i <= k} (n + i - 1)! 4^i / ((n - i)! (2i)!)
  d <- vector("list", n + 1L); acc <- gmp::as.bigq(0)
  for (i in 0:n) {
    acc <- acc + gmp::as.bigq(gmp::factorialZ(n + i - 1) * gmp::as.bigz(4)^i, gmp::factorialZ(n - i) * gmp::factorialZ(2 * i))
    d[[i + 1L]] <- n * acc
  }
  dn <- d[[n + 1L]]
  s <- .mp_zero
  for (k in 0:(n - 1)) {
    c <- .mp_from_q(d[[k + 1L]] - dn, wp)
    pw <- .mp_exp(.mp_mul(a, .mp_log(.mpf(k + 1, 0), wp), wp), wp)
    t <- .mp_div(c, pw, wp)
    s <- if (k %% 2 == 0) .mp_add(s, t, wp) else .mp_sub(s, t, wp)
  }
  den <- .mp_mul(.mp_from_q(dn, wp), .mp_sub(.mpf(1, 0), .mp_exp(.mp_mul(.mp_sub(.mpf(1, 0), a, wp), .mp_log(.mpf(2, 0), wp), wp), wp), wp), wp)
  .mp_neg(.mp_div(s, den, p))
}

# Euler's constant by the Brent-McMillan algorithm
.mp_euler <- function(p) {
  key <- paste0("gamma", p)
  if (!is.null(v <- .mp_cache[[key]])) return(v)
  wp <- p + 30
  N <- ceiling(wp * 0.1733) + 2
  w2 <- wp + 20; one <- .pow2(w2)
  A <- -.fxp(.mp_log(.mpf(N, 0), w2), w2); B <- one
  U <- A; V <- B; N2 <- gmp::as.bigz(N)^2
  for (k in 1:ceiling(3.6 * N + 10)) {
    B <- .tdiv(B * N2, gmp::as.bigz(k)^2)
    A <- .tdiv(.tdiv(A * N2, gmp::as.bigz(k)) + B, gmp::as.bigz(k))
    U <- U + A; V <- V + B
    if (A == 0 && B == 0) break
  }
  v <- .mp_div(.fx_mp(U, w2, wp), .fx_mp(V, w2, wp), p)
  assign(key, v, envir = .mp_cache)
  v
}

# Catalan's constant: pi/8*log(2 + sqrt(3)) + 3/8*sum (k!)^2/((2k)! (2k+1)^2)
.mp_catalan <- function(p) {
  wp <- p + 30
  w2 <- wp + 10; one <- .pow2(w2)
  t <- one; s <- one; k <- 1
  repeat {
    t <- .tdiv(t * gmp::as.bigz(k), gmp::as.bigz(2 * (2 * k - 1)))
    if (t == 0) break
    s <- s + .tdiv(t, gmp::as.bigz(2 * k + 1)^2); k <- k + 1
  }
  sm <- .fx_mp(3 * s, w2, wp)
  l <- .mp_log(.mp_add(.mpf(2, 0), .mp_sqrt(.mpf(3, 0), wp), wp), wp)
  .mp_ldexp(.mp_add(.mp_mul(.mp_pi(wp), l, wp), sm, p), -3)
}

# ---------------------------------------------------------------------------
# Floats of arbitrary precision

# The binary value and precision of a Float node; plain Floats are doubles.
.float_mp <- function(e) if (!is.null(e$extra$mp)) e$extra$mp else .mp_from_dbl(e$value)
.float_prec <- function(e) if (!is.null(e$extra$prec)) e$extra$prec else 53
.mp_float <- function(a, prec) {
  if (prec == 53) {
    v <- .mp_to_dbl(a)
    # a double that underflowed or is subnormal keeps the binary value
    if (is.finite(v) && (a$m == 0 || abs(v) >= 2.2250738585072014e-308)) return(.float_raw(v))
  }
  new_sym("Float", value = .mp_to_dbl(a), extra = list(mp = a, prec = prec),
          key = paste0("FM", as.character(a$m), "e", a$e, "p", prec))
}
.is_mp_float <- function(e) e$type == "Float" && !is.null(e$extra$mp)

# Float arithmetic at the larger precision, as SymPy's Float
.num_mp_op <- function(op, a, b) {
  p <- max(if (a$type == "Float") .float_prec(a) else 0, if (b$type == "Float") .float_prec(b) else 0)
  ma <- if (a$type == "Float") .float_mp(a) else .mp_from_q(a$value, p + 10)
  mb <- if (b$type == "Float") .float_mp(b) else .mp_from_q(b$value, p + 10)
  r <- switch(op, add = .mp_add(ma, mb, p), mul = .mp_mul(ma, mb, p), div = .mp_div(ma, mb, p))
  .mp_float(r, p)
}

# str(Float): mpmath's to_str at the precision's digits; trailing zeros are
# stripped when the Float is printed inside another expression
.fmt_float_sym <- function(e, strip = TRUE) {
  v <- e$value
  if (is.nan(v)) return("nan")
  if (is.null(e$extra$mp) && is.infinite(v)) return(if (v > 0) "oo" else "-oo")
  prec <- .float_prec(e)
  dps <- if (prec < 5) 0 else .prec_to_dps(prec)
  s <- .mp_to_str(.float_mp(e), max(dps, 1), strip_zeros = strip)
  if (startsWith(s, "-.0")) s <- paste0("-0.", substring(s, 4L)) else if (startsWith(s, ".0")) s <- paste0("0.", substring(s, 3L))
  s
}

# Sign of a real number expression certified by evaluation to increasing
# precision, as Expr._eval_is_extended_positive_negative: 1, -1, or NA when
# the value is complex, zero, or has no significant digits.
.mp_num_sign <- function(e) {
  z <- tryCatch(.mpeval_adaptive(e, 20L, list(), maxprec = 4000), error = function(err) NULL)
  if (is.null(z) || !isTRUE(attr(z, "converged")) || z$im$m != 0 || z$re$m == 0) return(NA_real_)
  if (z$re$m > 0) 1 else -1
}

# An mpc value as a Float, an imaginary Float multiple of I, or their sum.
.mpc_sym <- function(r, prec, zero_int = FALSE) {
  # EvalfMixin.evalf: a zero result is the integer 0
  if (zero_int && r$re$m == 0 && r$im$m == 0) return(S_Zero)
  re <- .mp_float(r$re, prec)
  if (r$im$m == 0) return(re)
  im <- .mp_float(r$im, prec)
  if (r$re$m == 0) return(.mul(list(im, S_I)))
  .add(list(re, .mul(list(im, S_I))))
}

# The precision of the Floats in a list of numbers (53 when there are none).
.floats_prec <- function(args) {
  ps <- vapply(args, function(a) if (a$type == "Float") .float_prec(a) else 0, 0)
  if (any(ps > 0)) max(ps) else 53
}
