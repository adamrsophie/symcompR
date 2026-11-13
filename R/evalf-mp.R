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
  if (gmp::mod.bigz(m, 2) == 0) {
    tz <- .v2(abs(m))
    m <- gmp::divq.bigz(m, .pow2(tz)); e <- e + tz
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
  m <- a$m; e <- a$e
  if (b > 60) { s <- b - 60; m <- gmp::divq.bigz(m, .pow2(s)); e <- e + s }
  .ldexp_dbl(as.numeric(m), e)
}
# x * 2^e in steps, so that a subnormal result is not lost to an
# intermediate underflow of 2^e
.ldexp_dbl <- function(x, e) {
  while (e > 1000) { x <- x * 2^1000; e <- e - 1000 }
  while (e < -1000 && x != 0) { x <- x * 2^-1000; e <- e + 1000 }
  x * 2^e
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
# A constant in fixed point is computed once at a precision a little above
# the largest one requested so far and shifted down for lower precisions,
# as mpmath's constant_memo.
.fx_const_memo <- function(name, wp, compute) {
  memo <- .mp_cache[[name]]
  if (is.null(memo) || memo$prec < wp) {
    w <- as.integer(wp * 1.05 + 10)
    memo <- list(prec = w, v = compute(w))
    assign(name, memo, envir = .mp_cache)
  }
  gmp::divq.bigz(memo$v, .pow2(memo$prec - wp))
}
.fx_pi <- function(wp) .fx_const_memo("pi", wp, function(w)
  gmp::divq.bigz(16 * .fx_atan_inv(5, w + 10) - 4 * .fx_atan_inv(239, w + 10), .pow2(10)))
.fx_ln2 <- function(wp) .fx_const_memo("ln2", wp, function(w) {
  # log(2) = 2*atanh(1/3) = 2 * sum 1/((2k+1) 3^(2k+1))
  w <- w + 10; one <- .pow2(w)
  term <- gmp::divq.bigz(one, 3); s <- term; k <- 1
  repeat {
    term <- gmp::divq.bigz(term, 9)
    if (term == 0) break
    s <- s + gmp::divq.bigz(term, 2 * k + 1); k <- k + 1
  }
  gmp::divq.bigz(2 * s, .pow2(10))
})
.mp_pi <- function(p) .fx_mp(.fx_pi(p + 20), p + 20, p)
.mp_ln2 <- function(p) .fx_mp(.fx_ln2(p + 20), p + 20, p)

# Fixed point helpers: x / 2^wp with signed floor shifts, as Python's >>
.fx_shift <- function(x, k) if (k >= 0) x * .pow2(k) else gmp::divq.bigz(x, .pow2(-k))
.to_fixed <- function(a, wp) .fx_shift(a$m, a$e + wp)

# exp(x / 2^prec) * 2^prec for a fixed point x: mpmath's exp_basecase, the
# even and odd Taylor sums of x / 2^r followed by r squarings, or
# exponential_series above EXP_COSH_CUTOFF (600 bits with the Python backend)
.fx_exp_basecase <- function(x, prec) {
  if (prec > 600) return(.fx_exponential_series(x, prec))
  r <- as.integer(floor(sqrt(prec)))
  prec <- prec + r
  s0 <- s1 <- .pow2(prec)
  k <- 2
  a <- x2 <- .fx_shift(x * x, -prec)
  while (a != 0) {
    a <- gmp::divq.bigz(a, k); s0 <- s0 + a; k <- k + 1
    a <- gmp::divq.bigz(a, k); s1 <- s1 + a; k <- k + 1
    a <- .fx_shift(a * x2, -prec)
  }
  s1 <- .fx_shift(s1 * x, -prec)
  s <- s0 + s1
  for (i in seq_len(r)) s <- .fx_shift(s * s, -prec)
  .fx_shift(s, -r)
}
# exponential_series(x, prec, type = 0): exp(x) from cosh(x) by
# exp = cosh + sqrt(cosh^2 - 1), with the argument halved r times
.fx_exponential_series <- function(x, prec) {
  sign <- x < 0
  if (sign) x <- -x
  r <- as.integer(floor(0.5 * sqrt(prec)))
  xmag <- .bitlen(x) - prec
  r <- max(0, xmag + r)
  extra <- 10 + 2 * max(r, -xmag)
  wp <- prec + extra
  x <- .fx_shift(x, extra - r)
  one <- .pow2(wp)
  if (prec < 1500) {
    x2 <- a <- .fx_shift(x * x, -wp)
    x4 <- .fx_shift(x2 * x2, -wp)
    s0 <- s1 <- gmp::as.bigz(0)
    k <- 2
    while (a != 0) {
      a <- gmp::divq.bigz(a, (k - 1) * k); s0 <- s0 + a; k <- k + 2
      a <- gmp::divq.bigz(a, (k - 1) * k); s1 <- s1 + a; k <- k + 2
      a <- .fx_shift(a * x4, -wp)
    }
    s1 <- .fx_shift(x2 * s1, -wp)
    c <- s1 + s0 + one
  } else {
    u <- as.integer(floor(0.3 * prec^0.35))
    x2 <- a <- .fx_shift(x * x, -wp)
    xpowers <- list(one, x2)
    for (i in seq_len(u - 1L)) xpowers[[length(xpowers) + 1L]] <- .fx_shift(xpowers[[length(xpowers)]] * x2, -wp)
    sums <- rep(list(gmp::as.bigz(0)), u)
    k <- 2
    while (a != 0) {
      for (i in seq_len(u)) {
        a <- gmp::divq.bigz(a, (k - 1) * k)
        sums[[i]] <- sums[[i]] + a
        k <- k + 2
      }
      a <- .fx_shift(a * xpowers[[length(xpowers)]], -wp)
    }
    for (i in seq_len(u)[-1L]) sums[[i]] <- .fx_shift(sums[[i]] * xpowers[[i]], -wp)
    c <- Reduce(`+`, sums) + one
  }
  s <- .mp_isqrt(c * c - .pow2(wp) * one)
  v <- if (sign) c - s else c + s
  for (i in seq_len(r)) v <- .fx_shift(v * v, -wp)
  .fx_shift(v, -extra)
}

# exp(x) for real x, after mpmath's mpf_exp: x = n*log(2) + t with
# 0 <= t < log(2) in fixed point, the exponent n kept as an exact integer
.mp_exp <- function(a, p) {
  if (a$m == 0) return(.mpf(1, 0))
  mag <- .mp_mag(a)
  wp <- p + 14
  if (mag < -wp) return(.mpf(1, 0))
  if (mag > 1) {
    wpmod <- wp + mag
    t <- .to_fixed(a, wpmod)
    lg2 <- .fx_ln2(wpmod)
    n <- gmp::divq.bigz(t, lg2)
    t <- .fx_shift(t - n * lg2, -mag)
  } else {
    t <- .to_fixed(a, wp)
    n <- gmp::as.bigz(0)
  }
  # exponents are doubles, exact below 2^53. A huge result cannot be held;
  # a tiny one (a quadrature node far out) only enters sums as negligible,
  # so its exponent may be rounded. .mpc_sym refuses such a final value.
  if (n > 2^52) stop("exponent overflow in multiprecision exp", call. = FALSE)
  .mp_round(.fx_exp_basecase(t, wp), as.numeric(n) - wp, p)
}

# log(x) for real x > 0: x = y * 2^n with y near 1, log(y) = 2*atanh(z)
# with z = (y-1)/(y+1); the fixed point series is scaled to the size of z,
# so x near 1 keeps its relative accuracy
.mp_log <- function(a, p, extra = 0) {
  if (a$m <= 0) stop("log of a nonpositive number in real arithmetic", call. = FALSE)
  wp <- p + 20 + extra
  b <- .bitlen(a$m)
  n <- b + a$e
  y <- list(m = a$m, e = -b)                  # y in [1/2, 1)
  if (.mp_cmp(y, .mp_sqrt(.mpf(1, -1), wp)) < 0) { y <- .mp_ldexp(y, 1); n <- n - 1 }
  num <- .mp_add(y, .mpf(-1, 0), .bitlen(y$m) + abs(y$e) + 2)
  if (num$m == 0) {
    ly <- .mp_zero
  } else {
    z <- .mp_div(num, .mp_add(y, .mpf(1, 0), wp), wp)
    w2 <- wp + 10 + max(0, -.mp_mag(z)); one <- .pow2(w2)
    zf <- .to_fixed(z, w2); z2 <- gmp::divq.bigz(zf * zf, one)
    s <- zf; term <- zf; k <- 1
    repeat {
      term <- .tdiv(term * z2, one)
      if (term == 0) break
      s <- s + .tdiv(term, 2 * k + 1); k <- k + 1
    }
    ly <- .fx_mp(2 * s, w2, wp)
  }
  if (n == 0) return(.mp_round(ly$m, ly$e, p))
  .mp_add(ly, .mp_mul(.mp_round(gmp::as.bigz(n), 0, wp), .fx_mp(.fx_ln2(wp), wp, wp), wp), p)
}

# cos(t / 2^wp), sin(t / 2^wp) in fixed point for 0 <= t <= pi/2: Taylor
# series with guard bits
.fx_cos_sin_basecase <- function(t, wp) {
  w2 <- wp + 12; one <- .pow2(w2)
  x <- .fx_shift(t, 12); x2 <- gmp::divq.bigz(x * x, one)
  s <- x; c <- one; ts <- x; tc <- one; k <- 1
  repeat {
    ts <- -.tdiv(.tdiv(ts * x2, one), gmp::as.bigz((2 * k) * (2 * k + 1)))
    tc <- -.tdiv(.tdiv(tc * x2, one), gmp::as.bigz((2 * k - 1) * (2 * k)))
    if (ts == 0 && tc == 0) break
    s <- s + ts; c <- c + tc; k <- k + 1
  }
  list(c = .fx_shift(c, -12), s = .fx_shift(s, -12))
}

# cos and sin of x, or of pi*x, after mpmath's mpf_cos_sin: reduction modulo
# pi/2 with the working precision raised until the reduced argument is
# resolved (mod_pi2); for pi*x the nearest half integer is subtracted
# exactly, so integers and half integers give exact zeros.
.mp_cos_sin <- function(a, p, pi = FALSE) {
  if (a$m == 0) return(list(cos = .mpf(1, 0), sin = .mp_zero))
  sign <- a$m < 0
  man <- abs(a$m); ex <- a$e
  mag <- .bitlen(man) + ex
  wp <- p + 10
  if (mag < 0 && mag < -wp) {
    x <- if (pi) .mp_mul(a, .mp_pi(wp), wp + 4) else a
    return(list(cos = .mpf(1, 0), sin = .mp_round(x$m, x$e, p)))
  }
  if (pi) {
    if (ex >= -1) {
      if (ex == -1) {
        odd <- gmp::mod.bigz(gmp::divq.bigz(man, 2), 2) == 1
        return(list(cos = .mp_zero, sin = .mpf(if (xor(odd, sign)) -1 else 1, 0)))
      }
      if (ex == 0) return(list(cos = .mpf(-1, 0), sin = .mp_zero))
      return(list(cos = .mpf(1, 0), sin = .mp_zero))
    }
    n <- gmp::divq.bigz(gmp::divq.bigz(man, .pow2(-ex - 2)) + 1, 2)
    man <- man - n * .pow2(-ex - 1)
    mag2 <- .bitlen(abs(man)) + ex
    wp <- p + 10 - mag2
    t <- .fx_shift(man, ex + wp)
    t <- gmp::divq.bigz(t * .fx_pi(wp), .pow2(wp))
    # the reduced argument may be negative: cos is even, sin odd
    tneg <- t < 0
    cs <- .fx_cos_sin_basecase(abs(t), wp)
    if (tneg) cs$s <- -cs$s
  } else {
    if (mag > 0) {
      i <- 0
      repeat {
        cp <- 20 * 2^i
        wpmod <- wp + mag + cp
        pi2 <- .fx_pi(wpmod - 1)
        pi4 <- gmp::divq.bigz(pi2, 2)
        t <- .fx_shift(man, wpmod + ex)
        n <- gmp::divq.bigz(t, pi2); y <- t - n * pi2
        small <- if (y > pi4) pi2 - y else y
        if (gmp::divq.bigz(small, .pow2(wp + mag - 10)) != 0) {
          t <- gmp::divq.bigz(y, .pow2(mag))
          wp <- wpmod - mag
          break
        }
        i <- i + 1
      }
    } else {
      wp <- wp - mag
      t <- .fx_shift(man, ex + wp)
      n <- gmp::as.bigz(0)
    }
    cs <- .fx_cos_sin_basecase(t, wp)
  }
  c <- cs$c; s <- cs$s
  m <- as.integer(gmp::mod.bigz(n, 4))
  if (m == 1) { tmp <- c; c <- -s; s <- tmp }
  else if (m == 2) { c <- -c; s <- -s }
  else if (m == 3) { tmp <- c; c <- s; s <- -tmp }
  if (sign) s <- -s
  list(cos = .mp_round(c, -wp, p), sin = .mp_round(s, -wp, p))
}
.mp_sincos <- function(a, p) .mp_cos_sin(a, p)

# ---- complex helpers (libmpc) ----
.mp_xmul <- function(a, b) list(m = a$m * b$m, e = a$e + b$e)
.mpc_mul <- function(a, b, p) {
  if (.mpc_real(a) && .mpc_real(b)) return(.mpc(.mp_mul(a$re, b$re, p)))
  # exact products, one rounding, as mpc_mul
  .mpc(.mp_sub(.mp_xmul(a$re, b$re), .mp_xmul(a$im, b$im), p),
       .mp_add(.mp_xmul(a$re, b$im), .mp_xmul(a$im, b$re), p))
}
.mpc_div <- function(a, b, p) {
  if (.mpc_real(b)) return(.mpc(.mp_div(a$re, b$re, p), .mp_div(a$im, b$re, p)))
  w <- p + 10
  n <- .mp_add(.mp_xmul(b$re, b$re), .mp_xmul(b$im, b$im), w)
  .mpc(.mp_div(.mp_add(.mp_xmul(a$re, b$re), .mp_xmul(a$im, b$im), w), n, p),
       .mp_div(.mp_sub(.mp_xmul(a$im, b$re), .mp_xmul(a$re, b$im), w), n, p))
}
.mpc_abs <- function(z, p) {
  if (z$im$m == 0) return(.mp_round(abs(z$re$m), z$re$e, p))
  if (z$re$m == 0) return(.mp_round(abs(z$im$m), z$im$e, p))
  .mp_sqrt(.mp_add(.mp_xmul(z$re, z$re), .mp_xmul(z$im, z$im), p + 4), p)
}
# log(sqrt(a^2 + b^2)) without cancellation near the unit circle
.mp_log_hypot <- function(a, b, p) {
  if (b$m == 0) { t <- a; a <- b; b <- t }
  if (a$m == 0) {
    if (b$m == 0) stop("log of zero", call. = FALSE)
    return(.mp_log(list(m = abs(b$m), e = b$e), p))
  }
  a2 <- .mp_xmul(a, a); b2 <- .mp_xmul(b, b)
  extra <- 20
  h2 <- .mp_add(a2, b2, p + extra)
  cancelled <- .mp_add(h2, .mpf(-1, 0), 10)
  if (cancelled$m == 0 || .mp_mag(cancelled) < -extra %/% 2)
    h2 <- .mp_add(a2, b2, p + extra - min(a2$e, b2$e))
  .mp_ldexp(.mp_log(h2, p), -1)
}
.mpc_log <- function(z, p) {
  if (.mpc_real(z) && z$re$m > 0) return(.mpc(.mp_log(z$re, p)))
  if (z$re$m == 0 && z$im$m == 0) stop("log of zero", call. = FALSE)
  .mpc(.mp_log_hypot(z$re, z$im, p), .mp_atan2(z$im, z$re, p))
}
# principal square root, as mpc_sqrt: the smaller part by division
.mpc_sqrt <- function(z, p) {
  if (.mpc_real(z)) {
    if (z$re$m >= 0) return(.mpc(.mp_sqrt(z$re, p)))
    return(.mpc(.mp_zero, .mp_sqrt(.mp_neg(z$re), p)))
  }
  a <- z$re; b <- z$im
  wp <- p + 20
  if (a$m >= 0) {
    t <- .mp_add(.mpc_abs(z, wp), a, wp)
    re <- .mp_sqrt(.mp_ldexp(t, -1), p)
    im <- .mp_div(b, .mp_sqrt(.mp_ldexp(t, 1), wp), p)
  } else {
    t <- .mp_sub(.mpc_abs(z, wp), a, wp)
    im <- .mp_sqrt(.mp_ldexp(t, -1), p)
    re <- .mp_div(b, .mp_sqrt(.mp_ldexp(t, 1), wp), p)
    if (b$m < 0) { re <- .mp_neg(re); im <- .mp_neg(im) }
  }
  .mpc(re, im)
}
# cosh and sinh of real x, after mpf_cosh_sinh: extra precision for tiny x
.mpc_cosh_sinh <- function(x, p) {
  if (x$m == 0) return(list(cosh = .mpf(1, 0), sinh = .mp_zero))
  sign <- x$m < 0
  mag <- .mp_mag(x)
  wp <- p + 14
  if (mag < -4) {
    if (mag < -wp) return(list(cosh = .mpf(1, 0), sinh = .mp_round(x$m, x$e, p)))
    wp <- wp - mag
  }
  ax <- list(m = abs(x$m), e = x$e)
  if (mag > 10 && 3 * 2^(mag - 1) > wp) {
    c <- .mp_ldexp(.mp_exp(ax, p), -1)
    return(list(cosh = c, sinh = if (sign) .mp_neg(c) else c))
  }
  if (mag > 1) {
    wpmod <- wp + mag
    t <- .to_fixed(ax, wpmod)
    lg2 <- .fx_ln2(wpmod)
    n <- gmp::divq.bigz(t, lg2)
    t <- .fx_shift(t - n * lg2, -mag)
    n <- as.numeric(n)
  } else {
    t <- .to_fixed(ax, wp); n <- 0
  }
  ea <- .fx_exp_basecase(t, wp)
  eb <- gmp::divq.bigz(.pow2(2 * wp), ea)
  ebs <- gmp::divq.bigz(eb, .pow2(2 * n))
  ch <- ea + ebs; sh <- ea - ebs
  if (sign) sh <- -sh
  list(cosh = .mp_round(ch, n - wp - 1, p), sinh = .mp_round(sh, n - wp - 1, p))
}
.mpc_exp <- function(z, p) {
  if (z$im$m == 0) return(.mpc(.mp_exp(z$re, p)))
  if (z$re$m == 0) { cs <- .mp_cos_sin(z$im, p); return(.mpc(cs$cos, cs$sin)) }
  mg <- .mp_exp(z$re, p + 4)
  cs <- .mp_cos_sin(z$im, p + 4)
  .mpc(.mp_mul(mg, cs$cos, p), .mp_mul(mg, cs$sin, p))
}
.mpc_sin <- function(z, p) {
  if (z$im$m == 0) return(.mpc(.mp_cos_sin(z$re, p)$sin))
  if (z$re$m == 0) return(.mpc(.mp_zero, .mpc_cosh_sinh(z$im, p)$sinh))
  wp <- p + 6
  cs <- .mp_cos_sin(z$re, wp); ch <- .mpc_cosh_sinh(z$im, wp)
  .mpc(.mp_mul(cs$sin, ch$cosh, p), .mp_mul(cs$cos, ch$sinh, p))
}
.mpc_cos <- function(z, p) {
  if (z$im$m == 0) return(.mpc(.mp_cos_sin(z$re, p)$cos))
  if (z$re$m == 0) return(.mpc(.mpc_cosh_sinh(z$im, p)$cosh))
  wp <- p + 6
  cs <- .mp_cos_sin(z$re, wp); ch <- .mpc_cosh_sinh(z$im, wp)
  .mpc(.mp_mul(cs$cos, ch$cosh, p), .mp_neg(.mp_mul(cs$sin, ch$sinh, p)))
}
# cos(pi*z) and sin(pi*z), as mpc_cos_sin_pi
.mpc_cos_sin_pi <- function(z, p) {
  if (z$im$m == 0) { cs <- .mp_cos_sin(z$re, p, pi = TRUE); return(list(cos = .mpc(cs$cos), sin = .mpc(cs$sin))) }
  b <- .mp_mul(z$im, .mp_pi(p + 5), p + 5)
  if (z$re$m == 0) { ch <- .mpc_cosh_sinh(b, p); return(list(cos = .mpc(ch$cosh), sin = .mpc(.mp_zero, ch$sinh))) }
  wp <- p + 6
  cs <- .mp_cos_sin(z$re, wp, pi = TRUE); ch <- .mpc_cosh_sinh(b, wp)
  list(cos = .mpc(.mp_mul(cs$cos, ch$cosh, p), .mp_neg(.mp_mul(cs$sin, ch$sinh, p))),
       sin = .mpc(.mp_mul(cs$sin, ch$cosh, p), .mp_mul(cs$cos, ch$sinh, p)))
}
# (a + b i)^n exactly for n >= 0
.complex_int_pow <- function(a, b, n) {
  wre <- gmp::as.bigz(1); wim <- gmp::as.bigz(0)
  n <- gmp::as.bigz(n)
  while (n > 0) {
    if (gmp::mod.bigz(n, 2) == 1) {
      t <- wre * a - wim * b; wim <- wim * a + wre * b; wre <- t
      n <- n - 1
    }
    t <- a * a - b * b; b <- 2 * a * b; a <- t
    n <- gmp::divq.bigz(n, 2)
  }
  list(re = wre, im = wim)
}
.mpc_reciprocal <- function(z, p) {
  m <- .mp_add(.mp_xmul(z$re, z$re), .mp_xmul(z$im, z$im), p + 10)
  .mpc(.mp_div(z$re, m, p), .mp_neg(.mp_div(z$im, m, p)))
}
.mpc_square <- function(z, p) {
  .mpc(.mp_sub(.mp_xmul(z$re, z$re), .mp_xmul(z$im, z$im), p), .mp_ldexp(.mp_mul(z$re, z$im, p), 1))
}
# z^n for an integer n, as mpc_pow_int
.mpc_pow_int <- function(z, n, p) {
  n <- gmp::as.bigz(n)
  if (z$im$m == 0) return(.mpc(.mp_pow_int(z$re, n, p)))
  if (z$re$m == 0) {
    v <- .mp_pow_int(z$im, n, p)
    return(switch(as.integer(gmp::mod.bigz(n, 4)) + 1L, .mpc(v), .mpc(.mp_zero, v), .mpc(.mp_neg(v)), .mpc(.mp_zero, .mp_neg(v))))
  }
  if (n == 0) return(.mpc_one)
  if (n == 1) return(.mpc(.mp_round(z$re$m, z$re$e, p), .mp_round(z$im$m, z$im$e, p)))
  if (n == 2) return(.mpc_square(z, p))
  if (n == -1) return(.mpc_reciprocal(z, p))
  if (n < 0) return(.mpc_reciprocal(.mpc_pow_int(z, -n, p + 4), p))
  a <- z$re; b <- z$im
  de <- a$e - b$e
  exact_size <- as.numeric(n) * (abs(de) + max(.bitlen(a$m), .bitlen(b$m)))
  if (exact_size < 10000) {
    am <- a$m; bm <- b$m; ae <- a$e; be <- b$e
    if (de > 0) { am <- am * .pow2(de); ae <- be } else { bm <- bm * .pow2(-de); be <- ae }
    r <- .complex_int_pow(am, bm, n)
    nn <- as.numeric(n)
    return(.mpc(.mp_round(r$re, nn * ae, p), .mp_round(r$im, nn * be, p)))
  }
  l <- .mpc_log(z, p + 10)
  .mpc_exp(.mpc(.mp_mul(l$re, .mpf(n, 0), p + 10), .mp_mul(l$im, .mpf(n, 0), p + 10)), p)
}

# x^n for an integer n, as mpf_pow_int: exact when small, otherwise binary
# powering with truncation at p + 4*bitlen(n) + 4 bits
.mp_pow_int <- function(a, n, p) {
  n <- gmp::as.bigz(n)
  if (n == 0) return(.mpf(1, 0))
  if (a$m == 0) {
    if (n < 0) stop("division by zero in multiprecision arithmetic", call. = FALSE)
    return(.mp_zero)
  }
  if (n == 1) return(.mp_round(a$m, a$e, p))
  if (n == -1) return(.mp_div(.mpf(1, 0), a, p))
  if (n < 0) return(.mp_div(.mpf(1, 0), .mp_pow_int(a, -n, p + 5), p))
  neg <- a$m < 0 && gmp::mod.bigz(n, 2) == 1
  man <- abs(a$m); ex <- a$e
  nn <- as.numeric(n)
  bc <- .bitlen(man)
  if (man == 1) return(list(m = if (neg) gmp::as.bigz(-1) else gmp::as.bigz(1), e = ex * nn))
  if (bc * nn < 1000) {
    r <- .mp_round(man^n, ex * nn, p)
    if (neg) r$m <- -r$m
    return(r)
  }
  workprec <- p + 4 * .bitlen(n) + 4
  pm <- gmp::as.bigz(1); pe <- 0
  repeat {
    if (gmp::mod.bigz(n, 2) == 1) {
      pm <- pm * man; pe <- pe + ex
      b <- .bitlen(pm)
      if (b > workprec) { pm <- gmp::divq.bigz(pm, .pow2(b - workprec)); pe <- pe + b - workprec }
      n <- n - 1
      if (n == 0) break
    }
    man <- man * man; ex <- ex + ex
    b <- .bitlen(man)
    if (b > workprec) { man <- gmp::divq.bigz(man, .pow2(b - workprec)); ex <- ex + b - workprec }
    n <- gmp::divq.bigz(n, 2)
  }
  .mp_round(if (neg) -pm else pm, pe, p)
}

# b^x after mpmath's mpf_pow, mpc_pow_mpf and mpc_pow; xq is an exact
# rational exponent when the caller knows one
.mpc_pow <- function(b, x, p, xq = NULL) {
  if (x$im$m == 0) {
    if (is.null(xq)) {
      t <- x$re
      if (t$m == 0) return(.mpc_one)
      t <- .mp_round(t$m, t$e, .bitlen(t$m))
      if (t$e >= 0) xq <- gmp::as.bigq(t$m * .pow2(t$e))
      else if (t$e == -1) xq <- gmp::as.bigq(t$m, 2)
    }
    if (!is.null(xq)) {
      num <- gmp::numerator(xq); den <- gmp::denominator(xq)
      if (den == 1) return(.mpc_pow_int(b, num, p))
      if (den == 2) {
        if (.mpc_real(b) && b$re$m >= 0) {
          if (num == 1) return(.mpc_sqrt(b, p))
          if (num == -1) return(.mpc(.mp_div(.mpf(1, 0), .mp_sqrt(b$re, p + 10), p)))
          return(.mpc(.mp_pow_int(.mp_sqrt(b$re, p + 10), num, p)))
        }
        return(.mpc_pow_int(.mpc_sqrt(b, p + 10), num, p))
      }
    }
    if (b$re$m == 0 && b$im$m == 0) return(.mpc(.mp_zero))
    if (.mpc_real(b) && b$re$m > 0) {
      if (.mp_cmp(b$re, .mpf(1, 0)) == 0) return(.mpc_one)
      return(.mpc(.mp_exp(.mp_xmul(x$re, .mp_log(b$re, p + 10)), p)))
    }
    l <- .mpc_log(b, p + 10)
    return(.mpc_exp(.mpc(.mp_mul(l$re, x$re, p + 10), .mp_mul(l$im, x$re, p + 10)), p))
  }
  if (b$re$m == 0 && b$im$m == 0) return(.mpc(.mp_zero))
  .mpc_exp(.mpc_mul(.mpc_log(b, p + 10), x, p + 10), p)
}

# ---- gamma functions ----
# Stirling's series for log gamma(w), real (mpf) or complex (mpc) w with |w|
# large enough that the terms fall below 2^-wp relative to the sum before
# they start to grow; an expansion that does not converge is an error, never
# a silently truncated value.
.mp_stirling_terms <- function(w, wp, cx) {
  if (cx) {
    lw <- .mpc_log(w, wp)
    half <- .mpc(.mpf(1, -1))
    s <- .mpc_sub(.mpc_mul(.mpc_sub(w, half, wp), lw, wp), w, wp)
    s <- .mpc_add(s, .mpc(.mp_half_log_2pi(wp)), wp)
    winv <- .mpc_div(.mpc_one, w, wp); w2 <- .mpc_mul(winv, winv, wp)
    mag <- function(v) suppressWarnings(max(.mp_mag(v$re), .mp_mag(v$im)))
  } else {
    lw <- .mp_log(w, wp)
    s <- .mp_sub(.mp_mul(.mp_sub(w, .mpf(1, -1), wp), lw, wp), w, wp)
    s <- .mp_add(s, .mp_half_log_2pi(wp), wp)
    winv <- .mp_div(.mpf(1, 0), w, wp); w2 <- .mp_mul(winv, winv, wp)
    mag <- .mp_mag
  }
  wpow <- winv
  smag <- mag(s)
  if (!is.finite(smag)) smag <- 0
  prev <- Inf
  k <- 1
  repeat {
    c <- .mp_stirling_coef(k, wp)
    t <- if (cx) .mpc(.mp_mul(c, wpow$re, wp), .mp_mul(c, wpow$im, wp)) else .mp_mul(c, wpow, wp)
    tm <- mag(t)
    if (!is.finite(tm) || tm < smag - wp - 4) break
    if (tm > prev) stop("Stirling series for the gamma function did not converge", call. = FALSE)
    s <- if (cx) .mpc_add(s, t, wp) else .mp_add(s, t, wp)
    prev <- tm
    wpow <- if (cx) .mpc_mul(wpow, w2, wp) else .mp_mul(wpow, w2, wp)
    k <- k + 1
  }
  s
}
# the series constants at wp bits, cached: log(2*pi)/2 and B_2k/(2k(2k-1))
.mp_half_log_2pi <- function(wp) {
  key <- paste0("hl2pi_", wp)
  if (!is.null(v <- .mp_cache[[key]])) return(v)
  v <- .mp_ldexp(.mp_log(.mp_ldexp(.mp_pi(wp), 1), wp), -1)
  assign(key, v, envir = .mp_cache)
  v
}
.mp_stirling_coef <- function(k, wp) {
  key <- paste0("stc", k, "_", wp)
  if (!is.null(v <- .mp_cache[[key]])) return(v)
  v <- .mp_from_q(.mp_bernoulli(2 * k) / (2 * k * (2 * k - 1)), wp)
  assign(key, v, envir = .mp_cache)
  v
}
# bits needed above p so that exp(log gamma) keeps p bits: |log gamma(z)|
# is about |z| log|z|
.mp_lgamma_extra <- function(mag) {
  if (!is.finite(mag) || mag <= 0) return(24)
  24 + mag + ceiling(log2(mag + 1))
}
.mp_exact_sub <- function(a, b) .mp_add(a, .mp_neg(b), max(.bitlen(a$m), .bitlen(b$m)) + abs(a$e - b$e) + 4)
.mp_is_int <- function(a) a$m == 0 || .mp_round(a$m, a$e, .bitlen(a$m))$e >= 0

# gamma(x) (rgamma = FALSE) or 1/gamma(x) (rgamma = TRUE) for real x, with
# one final rounding to p bits
.mp_gamma_r <- function(a, p, rgamma = FALSE) {
  if (.mp_is_int(a)) {
    n <- .mp_floor(a)
    if (n <= 0) {
      if (rgamma) return(.mp_zero)
      stop(.sympy_condition("ValueError", "gamma function pole"))
    }
    if (n < 3000) {
      f <- gmp::factorialZ(as.numeric(n) - 1)
      return(if (rgamma) .mp_div(.mpf(1, 0), list(m = f, e = 0), p) else .mp_round(f, 0, p))
    }
  }
  if (.mp_cmp(a, .mpf(1, -1)) < 0) {
    # reflection: gamma(x) = pi / (sin(pi x) gamma(1 - x))
    wp <- p + 15
    s <- .mp_cos_sin(a, wp, pi = TRUE)$sin
    g <- .mp_gamma_r(.mp_exact_sub(.mpf(1, 0), a), wp)
    if (rgamma) return(.mp_div(.mp_mul(s, g, wp), .mp_pi(wp), p))
    return(.mp_div(.mp_pi(wp), .mp_mul(s, g, wp), p))
  }
  mg <- .mp_mag(a)
  wp <- p + .mp_lgamma_extra(mg)
  z0 <- max(10, ceiling(wp / 4))
  shift <- 0
  if (mg <= 30) shift <- max(0, ceiling(z0 - .mp_to_dbl(a)))
  w <- if (shift > 0) .mp_add(a, .mpf(shift, 0), wp) else a
  lg <- .mp_stirling_terms(w, wp, FALSE)
  pr <- NULL
  if (shift > 0) {
    pr <- .mpf(1, 0)
    for (j in 0:(shift - 1)) pr <- .mp_mul(pr, .mp_add(a, .mpf(j, 0), wp + 10), wp + 10)
  }
  if (rgamma) {
    r <- .mp_exp(.mp_neg(lg), wp)
    return(if (is.null(pr)) .mp_round(r$m, r$e, p) else .mp_mul(r, pr, p))
  }
  r <- .mp_exp(lg, wp)
  if (is.null(pr)) .mp_round(r$m, r$e, p) else .mp_div(r, pr, p)
}
.mp_gamma <- function(a, p) .mp_gamma_r(a, p)
.mp_rgamma <- function(a, p) .mp_gamma_r(a, p, TRUE)

# complex gamma and its reciprocal
.mpc_gamma_r <- function(z, p, rgamma = FALSE) {
  if (z$im$m == 0) return(.mpc(.mp_gamma_r(z$re, p, rgamma)))
  if (.mp_cmp(z$re, .mpf(1, -1)) < 0) {
    wp <- p + 15
    s <- .mpc_cos_sin_pi(z, wp)$sin
    g <- .mpc_gamma_r(.mpc(.mp_exact_sub(.mpf(1, 0), z$re), .mp_neg(z$im)), wp)
    pi <- .mpc(.mp_pi(wp))
    if (rgamma) return(.mpc_div(.mpc_mul(s, g, wp), pi, p))
    return(.mpc_div(pi, .mpc_mul(s, g, wp), p))
  }
  mg <- max(.mp_mag(z$re), .mp_mag(z$im))
  wp <- p + .mp_lgamma_extra(mg)
  z0 <- max(10, ceiling(wp / 4))
  shift <- 0
  if (mg <= 30) {
    zr <- .mp_to_dbl(z$re); zi <- .mp_to_dbl(z$im)
    if (zr * zr + zi * zi < z0 * z0) shift <- max(0, ceiling(sqrt(max(0, z0 * z0 - zi * zi)) - zr))
  }
  w <- if (shift > 0) .mpc_add(z, .mpc(.mpf(shift, 0)), wp) else z
  lg <- .mp_stirling_terms(w, wp, TRUE)
  pr <- NULL
  if (shift > 0) {
    pr <- .mpc_one
    for (j in 0:(shift - 1)) pr <- .mpc_mul(pr, .mpc_add(z, .mpc(.mpf(j, 0)), wp + 10), wp + 10)
  }
  if (rgamma) {
    r <- .mpc_exp(.mpc_neg(lg), wp)
    return(if (is.null(pr)) .mpc(.mp_round(r$re$m, r$re$e, p), .mp_round(r$im$m, r$im$e, p)) else .mpc_mul(r, pr, p))
  }
  r <- .mpc_exp(lg, wp)
  if (is.null(pr)) .mpc(.mp_round(r$re$m, r$re$e, p), .mp_round(r$im$m, r$im$e, p)) else .mpc_div(r, pr, p)
}
.mpc_gamma <- function(z, p) .mpc_gamma_r(z, p)

# log gamma(x) for real x > 0, redone at higher precision when the result
# cancels (x near 1 or 2)
.mp_loggamma <- function(a, p) {
  extra <- 0
  repeat {
    mg <- .mp_mag(a)
    wp <- p + .mp_lgamma_extra(mg) + extra
    z0 <- max(10, ceiling(wp / 4))
    shift <- if (mg <= 30) max(0, ceiling(z0 - .mp_to_dbl(a))) else 0
    w <- if (shift > 0) .mp_add(a, .mpf(shift, 0), wp) else a
    lg <- .mp_stirling_terms(w, wp, FALSE)
    if (shift > 0) {
      pr <- .mpf(1, 0)
      for (j in 0:(shift - 1)) pr <- .mp_mul(pr, .mp_add(a, .mpf(j, 0), wp + 10), wp + 10)
      lp <- .mp_log(pr, wp)
      r <- .mp_sub(lg, lp, wp)
      if (r$m == 0) {
        if (extra > 4 * p + 200) return(.mp_zero)
        extra <- 2 * extra + p; next
      }
      lost <- max(.mp_mag(lg), .mp_mag(lp)) - .mp_mag(r)
      if (lost > 10 + extra) { extra <- lost + 10; next }
      return(.mp_round(r$m, r$e, p))
    }
    return(.mp_round(lg$m, lg$e, p))
  }
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


# the nd significant decimal digits of a (rounded half up) and the decimal
# exponent of the first one. Exact through rationals while those stay
# small; a huge or tiny exponent scales by a power of ten in binary
# arithmetic with a guard of 64 bits.
.mp_digits_exp <- function(a, nd) {
  ex <- floor((.mp_mag(a) - 1) * log10(2))
  exact <- abs(a$e) + .bitlen(a$m) < 20000 && abs(ex) < 6000
  wp <- ceiling(nd * log2(10)) + .bitlen(gmp::as.bigz(abs(ex) + 1)) + 64
  repeat {
    k <- nd - 1 - ex
    if (exact) {
      sc <- abs(.mp_to_q(a)) * gmp::as.bigq(10)^k
      n <- gmp::numerator(sc); d <- gmp::denominator(sc)
      digits <- gmp::divq.bigz(2 * n + d, 2 * d)
    } else {
      t <- .mp_mul(list(m = abs(a$m), e = a$e), .mp_pow_int(.mpf(5, 1), k, wp), wp)
      t <- .mp_add(t, .mpf(1, -1), wp)
      digits <- if (t$e >= 0) t$m * .pow2(t$e) else gmp::divq.bigz(t$m, .pow2(-t$e))
    }
    ns <- nchar(as.character(digits))
    if (ns > nd) { ex <- ex + 1; next }
    if (ns < nd) { ex <- ex - 1; next }
    break
  }
  list(digits = as.character(digits), ex = ex)
}

# decimal string with dps significant digits, after mpmath's to_str
.mp_to_str <- function(a, dps, strip_zeros = TRUE, min_fixed = NULL, max_fixed = NULL) {
  # mpmath's to_str: zero is always 0.0
  if (a$m == 0) return("0.0")
  # no digits (a Float of precision below 5): ".0" and the exponent,
  # rounded up when the leading digit is 5 or more
  if (dps == 0) {
    r <- .mp_digits_exp(a, 3)
    ex <- r$ex + (substr(r$digits, 1L, 1L) >= "5")
    return(paste0(if (a$m < 0) "-" else "", ".0e", if (ex >= 0) "+" else "-", format(abs(ex), scientific = FALSE)))
  }
  if (is.null(min_fixed)) min_fixed <- min(-(dps %/% 3), -5)
  if (is.null(max_fixed)) max_fixed <- dps
  neg <- a$m < 0
  r <- .mp_digits_exp(a, dps)
  ex <- r$ex; digits <- r$digits
  ds <- digits
  fixed <- min_fixed < ex && ex < max_fixed
  if (fixed) {
    if (ex >= 0) {
      ip <- substr(ds, 1L, ex + 1L); fp <- substring(ds, ex + 2L)
      if (nchar(ip) < ex + 1L) ip <- paste0(ip, strrep("0", ex + 1L - nchar(ip)))
    } else {
      ip <- "0"; fp <- paste0(strrep("0", -ex - 1L), ds)
    }
    if (strip_zeros) fp <- sub("0+$", "", fp)
    # a stripped empty fraction becomes .0; an unstripped one stays empty
    if (strip_zeros && !nzchar(fp)) fp <- "0"
    s <- paste0(ip, ".", fp)
  } else {
    fp <- substring(ds, 2L)
    if (strip_zeros) fp <- sub("0+$", "", fp)
    if (strip_zeros && !nzchar(fp)) fp <- "0"
    s <- paste0(substr(ds, 1L, 1L), ".", fp, "e", if (ex >= 0) "+" else "-", format(abs(ex), scientific = FALSE))
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
    gamma = if (real) .mpc(.mp_gamma(z$re, p)) else .mpc_gamma(z, p),
    loggamma = if (real && z$re$m > 0) .mpc(.mp_loggamma(z$re, p)) else NULL,
    factorial = if (real) .mpc(.mp_gamma(.mp_add(z$re, p1, w), p)) else .mpc_gamma(.mpc_add(z, one, w), p),
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
      r <- .mp_sympy_fun(e, p, env)
      if (!isFALSE(r)) return(r)
      h <- .FUNS[[e$name]]$evalf
      if (!is.null(h)) return(.mpeval(h(e), p, env))
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
  if (!is.null(r <- .mp_sympy_const(e, prec))) return(r)
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
.mp_euler <- function(p) .fx_mp(.fx_euler(p + 20), p + 20, p)
.fx_euler <- function(wp) .fx_const_memo("euler", wp, function(w) {
  wp <- w + 10
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
  gmp::divq.bigz(U * .pow2(w), V)
})

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
  s <- .mp_to_str(.float_mp(e), dps, strip_zeros = strip)
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
  if (abs(r$re$e) > 2^52 || abs(r$im$e) > 2^52)
    stop("exponent beyond the exact range of multiprecision numbers", call. = FALSE)
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
