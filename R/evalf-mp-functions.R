# Special functions as mpmath evaluates them (functions/bessel.py,
# functions/expintegrals.py and the fixed point routines of libmp/libhyper.py),
# on the mpmath number model of evalf-mp-ctx.R and the hypergeometric layer of
# evalf-mp-hyper.R. Each function follows its mpmath definition line by line;
# the SymPy entry points are registered at the end of this file.

# ---------------------------------------------------------------------------
# libmp helpers

# from_man_exp(m, e) without a precision: the exact normalized value
.lm <- function(m, e, p = Inf) .mp_round(m, e, p)
# Python's x // y (floor division) on big integers
.fdiv <- function(x, y) gmp::divq.bigz(gmp::as.bigz(x), gmp::as.bigz(y))
# sqrt_fixed(x, prec) = isqrt(x << prec)
.sqrt_fixed <- function(x, prec) .mp_isqrt(x * .pow2(prec))
.mp_one <- list(m = gmp::as.bigz(1), e = 0)

# ctx.fmul(a, b, prec = p): a product rounded to p bits, exact for p = 0
.hx_fmul <- function(a, b, prec = .MPX$prec) {
  old <- .MPX$prec
  on.exit(.MPX$prec <- old)
  .MPX$prec <- if (prec <= 0) Inf else prec
  .hx_mul(.mv_convert(a), .mv_convert(b))
}

# An exact rational hypergeometric parameter, as mpmath's (p, q) tuples
.mv_tuple <- function(p, q) .mv_new("tuple", q = gmp::as.bigq(p, q))

# ctx.nthroot(x, n) for an integer n >= 2: the root to prec + 20 bits, then
# rounded (a correctly rounded value; mpmath's Newton iteration is not ported)
.hx_nthroot <- function(x, n) {
  x <- .mv_convert(x); p <- .MPX$prec
  if (x$k == "mpf" && x$re$m >= 0) {
    if (x$re$m == 0) return(x)
    a <- x$re; w <- p + 20
    # floor(2^w * a^(1/n)) by integer n-th root
    k <- n * w - a$e
    sh <- k - (k %% n) + n
    m <- if (a$e + sh >= 0) a$m * .pow2(a$e + sh) else gmp::divq.bigz(a$m, .pow2(-(a$e + sh)))
    r <- .int_nthroot(m, n)
    return(.mv_mpf(.mp_round(2 * r + 1, -(sh %/% n) - 1, p)))
  }
  z <- .mv_as_mpc(x)
  w <- p + 20
  l <- .mpc_log(.mpc(z$re, z$im), w)
  r <- .mpc_exp(.mpc(.mp_div(l$re, .mpf(n, 0), w), .mp_div(l$im, .mpf(n, 0), w)), w)
  .mv_mpc(.mp_round(r$re$m, r$re$e, p), .mp_round(r$im$m, r$im$e, p))
}
.int_nthroot <- function(m, n) {
  if (m < 2) return(m)
  x <- .pow2(ceiling(.bitlen(m) / n) + 1)
  repeat {
    y <- gmp::divq.bigz((n - 1) * x + gmp::divq.bigz(m, x^(n - 1)), n)
    if (y >= x) break
    x <- y
  }
  while (x^n > m) x <- x - 1
  while ((x + 1)^n <= m) x <- x + 1
  x
}

# c_memo: a constant cached per context with the precision it was computed at
.SF_MEMO <- new.env(parent = emptyenv())
.hx_c_memo <- function(name, f) {
  key <- paste0(if (is.null(.MPX$ctx)) "mp" else .MPX$ctx, ":", name)
  prec <- .MPX$prec
  hit <- .SF_MEMO[[key]]
  if (!is.null(hit) && hit$p >= prec) return(.mv_pos(hit$v))
  v <- f()
  assign(key, list(p = prec, v = v), envir = .SF_MEMO)
  v
}
.mp_with_ctx <- function(ctx, expr) {
  old <- .MPX$ctx
  .MPX$ctx <- ctx
  on.exit(.MPX$ctx <- old)
  force(expr)
}

# ---------------------------------------------------------------------------
# Bessel functions (functions/bessel.py)

# mpf_besseljn and mpc_besseljn: J_n for an integer n in fixed point
.lm_besseljn <- function(n, x, prec) {
  negate <- n < 0 && gmp::mod.bigz(n, 2) == 1
  n <- abs(n)
  nn <- as.numeric(n)
  if (x$k == "mpf") {
    a <- x$re
    prec <- prec + 50
    mag <- if (a$m == 0) 0 else .mp_mag(a)
    wp <- prec + 20 + nn * .bitlen(n)
    if (mag < 0) wp <- wp - nn * mag
    xf <- .hx_fix(a, wp)
    x2 <- .hx_shr(xf^2, wp)
    if (nn == 0) {
      s <- t <- .pow2(wp)
    } else {
      s <- t <- .hx_shr(.fdiv(xf^nn, gmp::factorialZ(nn)), (nn - 1) * wp + nn)
    }
    k <- 1
    while (t != 0) {
      t <- .hx_shr(.fdiv(t * x2, -4 * k * (k + nn)), wp)
      s <- s + t
      k <- k + 1
    }
    if (negate) s <- -s
    # rounded to prec + 50 bits, as mpf_besseljn returns it
    return(.mv_mpf(.lm(s, -wp, prec)))
  }
  zre <- x$re; zim <- x$im
  origprec <- prec
  mz <- function(a) if (a$m == 0) 0 else .mp_mag(a)
  mag <- max(mz(zre), mz(zim))
  prec <- prec + 20 + nn * .bitlen(n) + abs(mag)
  if (mag < 0) prec <- prec - nn * mag
  zre <- .hx_fix(zre, prec); zim <- .hx_fix(zim, prec)
  z2re <- .hx_shr(zre^2 - zim^2, prec)
  z2im <- .hx_shr(zre * zim, prec - 1)
  if (nn == 0) {
    sre <- tre <- .pow2(prec); sim <- tim <- gmp::as.bigz(0)
  } else {
    ri <- .complex_int_pow(zre, zim, nn)
    sre <- tre <- .hx_shr(.fdiv(ri[[1]], gmp::factorialZ(nn)), (nn - 1) * prec + nn)
    sim <- tim <- .hx_shr(.fdiv(ri[[2]], gmp::factorialZ(nn)), (nn - 1) * prec + nn)
  }
  k <- 1
  while (abs(tre) + abs(tim) > 3) {
    pk <- -4 * k * (k + nn)
    t1 <- tre * z2re - tim * z2im; t2 <- tim * z2re + tre * z2im
    tre <- .hx_shr(.fdiv(t1, pk), prec)
    tim <- .hx_shr(.fdiv(t2, pk), prec)
    sre <- sre + tre; sim <- sim + tim
    k <- k + 1
  }
  if (negate) { sre <- -sre; sim <- -sim }
  .mv_mpc(.lm(sre, -prec, origprec), .lm(sim, -prec, origprec))
}

.hx_int_of <- function(x) {
  x <- .mv(x)
  if (x$k == "int") return(gmp::numerator(x$q))
  v <- .mp_norm(.mv_convert(.hx_re(x))$re)
  .mp_to_int(v)
}

.sf_besselj <- function(n, z, kw = list()) {
  n <- .mv(n)
  if (n$k == "int") {
    n_isint <- TRUE
  } else {
    n <- .mv_convert(n)
    n_isint <- .hx_isint(n)
    if (n_isint) n <- .mv_int(.hx_int_of(n))
  }
  if (n_isint && gmp::numerator(n$q) < 0) {
    ni <- gmp::numerator(n$q)
    return(.hx_mul(if (gmp::mod.bigz(ni, 2) == 1) -1 else 1, .sf_besselj(.mv_int(-ni), z, kw)))
  }
  z <- .mv_convert(z)
  M <- .hx_mag(z)
  if (n_isint && abs(M) < 10 && abs(gmp::numerator(n$q)) < 20)
    return(.lm_besseljn(gmp::numerator(n$q), z, .MPX$prec))
  if (.hx_iszero(z)) {
    if (.hx_iszero(n)) v <- .hx_add(.hx_add(.hx_one(), n), z)
    else if (.mv_cmp(.hx_re(n), 0) > 0) v <- .hx_mul(n, z)
    else v <- .hx_add(.hx_add(.mv_inf(), z), n)
  } else {
    orig <- .MPX$prec
    on.exit(.MPX$prec <- orig)
    .MPX$prec <- .MPX$prec + min(3 * abs(M), .MPX$prec)
    w <- .hx_fmul(z, 0.5, 0)
    h <- function(n) {
      r <- .hx_neg_exact(.hx_fmul(w, w, max(0, .MPX$prec + M)))
      list(.hx_term(list(w), list(n), list(), list(.hx_add(n, 1)), list(), list(.hx_add(n, 1)), r))
    }
    v <- .hx_hypercomb(h, list(n), kw = kw)
    .MPX$prec <- orig
  }
  .mv_pos(v)
}

.sf_besseli <- function(n, z, kw = list()) {
  n <- .mv_convert(n); z <- .mv_convert(z)
  if (!.hx_iszero(n) && .hx_isnpint(n)) return(.sf_besseli(.hx_neg(n), z, kw))
  if (.hx_iszero(z)) {
    if (.hx_iszero(n)) return(.hx_add(.hx_add(1, n), z))
    if (.hx_isint(n)) return(.hx_mul(0, .hx_add(n, z)))
    r <- .hx_re(n)
    if (.hx_iszero(r)) return(.hx_mul(.mv_nan(), .hx_add(n, z)))
    if (.mv_cmp(r, 0) > 0) return(.hx_mul(0, .hx_add(n, z)))
    return(.hx_add(.mv_inf(), .hx_add(n, z)))
  }
  M <- .hx_mag(z)
  h <- function(n) {
    w <- .hx_fmul(z, 0.5, 0)
    r <- .hx_fmul(w, w, max(0, .MPX$prec + M))
    list(.hx_term(list(w), list(n), list(), list(.hx_add(n, 1)), list(), list(.hx_add(n, 1)), r))
  }
  .hx_hypercomb(h, list(n), kw = kw)
}

.sf_bessely <- function(n, z, kw = list()) .hx_wrapped(function(n, z) {
  if (.hx_iszero(z)) {
    if (.hx_iszero(n)) return(.hx_add(.mv_ninf(), .hx_add(n, z)))
    if (!.hx_iszero(.hx_im(n))) return(.hx_mul(.mv_nan(), .hx_add(n, z)))
    r <- .hx_re(n)
    q <- .hx_add(n, 0.5)
    if (.hx_isint(q)) {
      if (.mv_cmp(n, 0) > 0) return(.hx_add(.mv_ninf(), .hx_add(n, z)))
      return(.hx_mul(0, .hx_add(n, z)))
    }
    if (.mv_cmp(r, 0) < 0 && gmp::mod.bigz(.mp_floor(.mv_convert(q)$re), 2) == 1)
      return(.hx_add(.mv_inf(), .hx_add(n, z)))
    return(.hx_add(.mv_ninf(), .hx_add(n, z)))
  }
  .MPX$prec <- .MPX$prec + 10
  d <- .hx_nint_distance(n)$d
  if (d < -.MPX$prec) {
    h <- .mv_pos(.hx_eps())
    .MPX$prec <- .MPX$prec * 2
    n <- .hx_add(n, h)
  } else if (d < 0) {
    .MPX$prec <- .MPX$prec - d
  }
  cs <- .hx_cospi_sinpi(n)
  .hx_div(.hx_sub(.hx_mul(.sf_besselj(n, z, kw), cs$cos), .sf_besselj(.hx_neg(n), z, kw)), cs$sin)
}, n, z)

.sf_besselk <- function(n, z, kw = list()) .hx_wrapped(function(n, z) {
  if (.hx_iszero(z)) return(.mv_inf())
  M <- .hx_mag(z)
  if (M < 1) {
    h <- function(n) {
      r <- .hx_pow(.hx_div(z, 2), 2)
      list(.hx_term(list(z, 2), list(.hx_neg(n), .hx_sub(n, 1)), list(n), list(), list(), list(.hx_sub(1, n)), r),
           .hx_term(list(z, 2), list(n, .hx_sub(.hx_neg(n), 1)), list(.hx_neg(n)), list(), list(), list(.hx_add(1, n)), r))
    }
  } else {
    if (is.finite(M)) .MPX$prec <- .MPX$prec + M
    h <- function(n) list(.hx_term(list(.hx_div(.hx_pi(), 2), z, .hx_exp(.hx_neg(z))), list(0.5, -0.5, 1),
                                   list(), list(), list(.hx_add(n, 0.5), .hx_sub(0.5, n)), list(),
                                   .hx_div(-1, .hx_mul(2, z))))
  }
  .hx_hypercomb(h, list(n), kw = kw)
}, n, z)

.sf_hankel <- function(n, x, sgn, kw = list()) .hx_wrapped(function(n, x) {
  jj <- if (sgn > 0) .hx_j() else .hx_neg(.hx_j())
  .hx_sum_accurately(function() list(.sf_besselj(n, x, kw), .hx_mul(jj, .sf_bessely(n, x, kw))))
}, n, x)

# ---------------------------------------------------------------------------
# Airy functions (functions/bessel.py), the function and its first derivative

.sf_airy_c <- function(name) .hx_c_memo(name, function() switch(name,
  ai1 = .hx_div(1, .hx_mul(.hx_nthroot(9, 3), .hx_gamma(.hx_div(.mv_mpf(.mpf(2, 0)), 3)))),
  ai2 = .hx_div(-1, .hx_mul(.hx_nthroot(3, 3), .hx_gamma(.hx_div(.mv_mpf(.mpf(1, 0)), 3)))),
  bi1 = .hx_div(1, .hx_mul(.hx_nthroot(3, 6), .hx_gamma(.hx_div(.mv_mpf(.mpf(2, 0)), 3)))),
  bi2 = .hx_div(.hx_nthroot(3, 6), .hx_gamma(.hx_div(.mv_mpf(.mpf(1, 0)), 3)))))

.sf_airy <- function(z, which, derivative = 0, kw = list()) {
  z <- .mv_convert(z)
  if (.mv_sp(z)) {
    if (derivative == 0) {
      if (which == "ai" && z$sp %in% c("inf", "-inf")) return(.hx_div(1, z))
      if (which == "bi" && z$sp == "inf") return(z)
      if (which == "bi" && z$sp == "-inf") return(.hx_div(1, z))
    }
    .hx_valueerror(sprintf("essential singularity of %s(z)", if (which == "ai") "Ai" else "Bi"))
  }
  extraprec <- if (!.hx_iszero(z)) max(0, as.integer(1.5 * .hx_mag(z))) else 0
  big <- function() {
    .MPX$prec <- .MPX$prec + extraprec
    w <- .hx_pow(z, 1.5); r <- .hx_div(-0.75, w); u <- .hx_div(.hx_mul(-2, w), 3)
    .MPX$prec <- .MPX$prec - extraprec
    list(r = r, u = u)
  }
  small_w <- function() {
    .MPX$prec <- .MPX$prec + extraprec
    w <- .hx_div(.hx_pow(z, 3), 9)
    .MPX$prec <- .MPX$prec - extraprec
    w
  }
  if (which == "ai" && derivative == 0) {
    h <- function() {
      if (.mv_cmp(.hx_re(z), 4) > 0) {
        b <- big()
        C <- .hx_div(.hx_exp(b$u), .hx_mul(.hx_mul(2, .hx_sqrt(.hx_pi())), .hx_nthroot(z, 4)))
        return(list(.hx_term(list(C), list(1), list(), list(), list(.mv_tuple(1, 6), .mv_tuple(5, 6)), list(), b$r)))
      }
      w <- small_w()
      list(.hx_term(list(.sf_airy_c("ai1")), list(1), list(), list(), list(), list(.mv_q(2, 3)), w),
           .hx_term(list(.hx_mul(z, .sf_airy_c("ai2"))), list(1), list(), list(), list(), list(.mv_q(4, 3)), w))
    }
  } else if (which == "ai") {
    h <- function() {
      if (.mv_cmp(.hx_re(z), 4) > 0) {
        b <- big()
        C <- .hx_mul(.hx_div(.hx_neg(.hx_exp(b$u)), .hx_mul(2, .hx_sqrt(.hx_pi()))), .hx_nthroot(z, 4))
        return(list(.hx_term(list(C), list(1), list(), list(), list(.mv_tuple(-1, 6), .mv_tuple(7, 6)), list(), b$r)))
      }
      w <- small_w()
      list(.hx_term(list(.hx_mul(.sf_airy_c("ai1"), 0.5), z), list(1, 2), list(), list(), list(), list(.mv_q(5, 3)), w),
           .hx_term(list(.sf_airy_c("ai2")), list(1), list(), list(), list(), list(.mv_q(1, 3)), w))
    }
  } else if (derivative == 0) {
    h <- function() {
      w <- small_w()
      list(.hx_term(list(.sf_airy_c("bi1")), list(1), list(), list(), list(), list(.mv_q(2, 3)), w),
           .hx_term(list(.hx_mul(z, .sf_airy_c("bi2"))), list(1), list(), list(), list(), list(.mv_q(4, 3)), w))
    }
  } else {
    h <- function() {
      w <- small_w()
      list(.hx_term(list(.hx_mul(.sf_airy_c("bi1"), 0.5), z), list(1, 2), list(), list(), list(), list(.mv_q(5, 3)), w),
           .hx_term(list(.sf_airy_c("bi2")), list(1), list(), list(), list(), list(.mv_q(1, 3)), w))
    }
  }
  .hx_hypercomb(h, list(), kw = kw)
}

# ---------------------------------------------------------------------------
# Error functions (libmp/libhyper.py and functions/expintegrals.py)

.lm_erf <- function(x, prec) {
  if (x$m == 0) return(.mp_zero)
  size <- .mp_mag(x)
  # erf(x) = 1 to within the precision
  if (size > 3 && 2 * (size - 1) + 0.528766 > log2(prec))
    return(.mpf(if (x$m < 0) -1 else 1, 0))
  # erf(x) ~ 2*x/sqrt(pi) close to 0
  if (size < -prec) {
    c <- .mp_sqrt(.mp_pi(prec + 20), prec + 20)
    return(.mp_div(.mp_ldexp(x, 1), c, prec))
  }
  wp <- prec + abs(size) + 25
  t <- abs(.hx_fix(x, wp))
  t2 <- .hx_shr(t * t, wp)
  s <- t; term <- 12345; k <- 1
  while (term != 0) {
    t <- .fdiv(.hx_shr(t * t2, wp), k)
    term <- .fdiv(t, 2 * k + 1)
    if (k %% 2 == 1) s <- s - term else s <- s + term
    k <- k + 1
  }
  s <- .fdiv(s * .pow2(wp + 1), .sqrt_fixed(.fx_pi(wp), wp))
  if (x$m < 0) s <- -s
  .lm(s, -wp, prec)
}

.lm_erfc <- function(x, prec) {
  if (x$m == 0) return(.mp_one)
  wp <- prec + 20
  mag <- .mp_mag(x)
  wp <- wp + max(0, 2 * mag)
  regular_erf <- x$m < 0 || mag < 2
  check <- function() { n <- .mp_to_int(x); as.numeric(n)^2 * 1.44 > wp }
  if (regular_erf || !check()) {
    if (regular_erf) return(.mp_sub(.mp_one, .lm_erf(x, prec + 10), prec))
    n <- as.numeric(.mp_to_int(x)) + 1
    return(.mp_sub(.mp_one, .lm_erf(x, prec + floor(n^2 * 1.44) + 10), prec))
  }
  s <- term <- .pow2(wp)
  term_prev <- gmp::as.bigz(0)
  t <- .hx_shr(2 * .hx_fix(x, wp)^2, wp)
  k <- 1
  repeat {
    term <- .fdiv(term * (2 * k - 1) * .pow2(wp), t)
    if ((k > 4 && term > term_prev) || term == 0) break
    if (k %% 2 == 1) s <- s - term else s <- s + term
    term_prev <- term
    k <- k + 1
  }
  s <- .fdiv(s * .pow2(wp), .sqrt_fixed(.fx_pi(wp), wp))
  s <- .lm(s, -wp, wp)
  z <- .mp_exp(.mp_round(-.mp_mul(x, x, wp)$m, .mp_mul(x, x, wp)$e, wp), wp)
  .mp_div(.mp_mul(z, s, wp), x, prec)
}

# square_exp_arg: z**2 (or 1/z**2) to 4*prec + 20 bits, times mult exactly
.sf_square_exp_arg <- function(z, mult = 1, reciprocal = FALSE) {
  prec <- .MPX$prec * 4 + 20
  z2 <- .hx_fmul(z, z, prec)
  if (reciprocal) {
    old <- .MPX$prec; .MPX$prec <- prec
    z2 <- .hx_div(.hx_one(), z2)
    .MPX$prec <- old
  }
  if (mult != 1) z2 <- .hx_fmul(z2, mult, 0)
  z2
}

.sf_erf_complex <- function(z) .hx_wrapped(function(z) {
  re_z <- .hx_re(z)
  if (.mv_cmp(re_z, 2) > 0) {
    v <- .sf_erf_complex(.hx_neg_exact(z))
    return(.hx_neg_exact(v))
  } else if (.mv_cmp(re_z, -2) < 0) {
    v <- .hx_sub(.sf_erfc_complex(.hx_neg_exact(z)), 1)
  } else {
    z2 <- .sf_square_exp_arg(z, -1)
    v <- .hx_mul(.hx_mul(.hx_div(2, .hx_sqrt(.hx_pi())), z),
                 .hx_hyper(list(.mv_tuple(1, 2)), list(.mv_tuple(3, 2)), z2))
  }
  if (.hx_iszero(re_z)) v <- .hx_mul(.hx_im(v), .hx_j())
  v
}, z)

.sf_erfc_complex <- function(z) .hx_wrapped(function(z) {
  re_z <- .hx_re(z)
  if (.mv_cmp(re_z, 2) > 0) {
    z2 <- .sf_square_exp_arg(z)
    nz2 <- .hx_neg_exact(z2)
    v <- .hx_mul(.hx_div(.hx_exp(nz2), .hx_sqrt(.hx_pi())), .sf_hyperu(.mv_tuple(1, 2), .mv_tuple(1, 2), z2))
  } else {
    v <- .hx_sub(1, .sf_erf_complex(z))
  }
  if (.hx_iszero(re_z)) v <- .hx_add(1, .hx_mul(.hx_im(v), .hx_j()))
  v
}, z)

.sf_erf <- function(z) {
  z <- .mv_convert(z)
  if (z$k == "mpf" && !.mv_sp(z)) return(.mv_mpf(.lm_erf(z$re, .MPX$prec)))
  if (z$k == "mpc" && !.mv_sp(z) && z$im$m == 0) return(.mv_mpc(.lm_erf(z$re, .MPX$prec), .mp_zero))
  .sf_erf_complex(z)
}
.sf_erfc <- function(z) {
  z <- .mv_convert(z)
  if (z$k == "mpf" && !.mv_sp(z)) return(.mv_mpf(.lm_erfc(z$re, .MPX$prec)))
  if (z$k == "mpc" && !.mv_sp(z) && z$im$m == 0) return(.mv_mpc(.lm_erfc(z$re, .MPX$prec), .mp_zero))
  .sf_erfc_complex(z)
}
.sf_erfi <- function(z) .hx_wrapped(function(z) {
  if (.hx_iszero(z)) return(z)
  z2 <- .sf_square_exp_arg(z)
  v <- .hx_mul(.hx_mul(.hx_div(2, .hx_sqrt(.hx_pi())), z), .hx_hyper(list(.mv_tuple(1, 2)), list(.mv_tuple(3, 2)), z2))
  if (.hx_iszero(.hx_re(z))) v <- .hx_mul(.hx_im(v), .hx_j())
  v
}, z)

# hyperu(a, b, z): the 2F0 series, else the connection formula
.sf_hyperu <- function(a, b, z, kw = list()) {
  pa <- .hx_param(a); pb <- .hx_param(b)
  a <- pa$v; b <- pb$v
  z <- .mv_convert(z)
  if (.hx_iszero(z)) {
    if (.mv_cmp(.hx_re(b), 1) <= 0) return(.hx_gammaprod(list(.hx_sub(1, b)), list(.hx_add(.hx_sub(a, b), 1))))
    return(.hx_add(.mv_inf(), z))
  }
  pbb <- .hx_param(.hx_sub(.hx_add(1, a), b))
  orig <- .MPX$prec
  v <- tryCatch({
    .MPX$prec <- .MPX$prec + 10
    v <- .hx_hypsum(2, 0, c(pa$t, pbb$t), list(a, pbb$v), .hx_div(-1, z), list(maxterms = .MPX$prec))
    .hx_div(v, .hx_pow(z, a))
  }, NoConvergence = function(err) NULL, finally = .MPX$prec <- orig)
  if (!is.null(v)) return(v)
  h <- function(a, b) {
    w <- .hx_sinpi(b)
    list(.hx_term(list(.hx_pi(), w), list(1, -1), list(), list(.hx_add(.hx_sub(a, b), 1), b), list(a), list(b), z),
         .hx_term(list(.hx_neg(.hx_pi()), w, z), list(1, -1, .hx_sub(1, b)), list(), list(a, .hx_sub(2, b)),
                  list(.hx_add(.hx_sub(a, b), 1)), list(.hx_sub(2, b)), z))
  }
  .hx_hypercomb(h, list(a, b), kw = kw)
}

# ---------------------------------------------------------------------------
# Exponential integrals (libmp/libhyper.py)

.lm_ei_taylor <- function(x, prec) {
  s <- t <- x; k <- 2
  while (t != 0) {
    t <- .fdiv(.hx_shr(t * x, prec), k)
    s <- s + .fdiv(t, k)
    k <- k + 1
  }
  s
}
.lm_complex_ei_taylor <- function(zre, zim, prec) {
  sre <- tre <- zre; sim <- tim <- zim; k <- 2
  while (abs(tre) + abs(tim) > 5) {
    a <- .hx_shr(.fdiv(tre * zre - tim * zim, k), prec)
    b <- .hx_shr(.fdiv(tre * zim + tim * zre, k), prec)
    tre <- a; tim <- b
    sre <- sre + .fdiv(tre, k); sim <- sim + .fdiv(tim, k)
    k <- k + 1
  }
  list(sre, sim)
}
.lm_ei_asymptotic <- function(x, prec) {
  one <- .pow2(prec)
  x <- t <- .fdiv(one * .pow2(prec), x)
  s <- one + x; k <- 2
  while (t != 0) {
    t <- .hx_shr(k * t * x, prec)
    s <- s + t; k <- k + 1
  }
  s
}
.lm_complex_ei_asymptotic <- function(zre, zim, prec) {
  one <- .pow2(prec)
  M <- .hx_shr(zim * zim + zre * zre, prec)
  xre <- tre <- .fdiv(zre * .pow2(prec), M)
  xim <- tim <- .fdiv(-zim * .pow2(prec), M)
  sre <- one + xre; sim <- xim; k <- 2
  while (abs(tre) + abs(tim) > 1000) {
    a <- .hx_shr((tre * xre - tim * xim) * k, prec)
    b <- .hx_shr((tre * xim + tim * xre) * k, prec)
    tre <- a; tim <- b
    sre <- sre + tre; sim <- sim + tim
    k <- k + 1
    if (k > prec) .hx_noconv()
  }
  list(sre, sim)
}

# mpf_ei(x) (e1 = TRUE: E1(x) = -Ei(-x)); NULL stands for ComplexResult
.lm_ei <- function(x, prec, e1 = FALSE) {
  if (e1) x <- .mp_neg(x)
  if (e1 && x$m >= 0) {
    if (x$m == 0) return("inf")
    return(NULL)
  }
  if (x$m == 0) return(if (e1) "inf" else "-inf")
  xabs <- list(m = abs(x$m), e = x$e)
  xmag <- .mp_mag(x)
  wp <- prec + 20
  can_use_asymp <- xmag > wp
  if (!can_use_asymp) {
    xabsint <- .mp_to_int(xabs)
    can_use_asymp <- xabsint > as.integer(wp * 0.693) + 10
  }
  if (can_use_asymp) {
    v <- if (xmag > wp) .mp_one else .lm(.lm_ei_asymptotic(.hx_fix(x, wp), wp), -wp)
    v <- .mp_mul(v, .mp_exp(x, wp), wp)
    v <- .mp_div(v, x, prec)
  } else {
    wp <- wp + 2 * as.numeric(.mp_to_int(xabs))
    u <- .hx_fix(x, wp)
    v <- .lm_ei_taylor(u, wp) + .fx_euler(wp)
    t1 <- .lm(v, -wp)
    t2 <- .mp_log(xabs, wp)
    v <- .mp_add(t1, t2, prec)
  }
  if (e1) v <- .mp_neg(v)
  v
}

.lm_ei_c <- function(z, prec, e1 = FALSE) {
  if (e1) z <- .mpc_neg(z)
  a <- z$re; b <- z$im
  if (b$m == 0) {
    if (e1) {
      x <- .lm_ei(a, prec)
      x <- if (is.character(x)) x else .mp_neg(x)
      y <- if (a$m >= 0) .mp_neg(.mp_pi(prec)) else .mp_zero
      return(list(x, y))
    }
    return(list(.lm_ei(a, prec), .mp_zero))
  }
  wp <- prec + 40
  mz <- function(a) if (a$m == 0) 0 else .mp_mag(a)
  zmag <- max(mz(a), mz(b))
  can_use_asymp <- zmag > wp
  if (!can_use_asymp) {
    zabsint <- abs(.mp_to_int(a)) + abs(.mp_to_int(b))
    can_use_asymp <- zabsint > as.integer(wp * 0.693) + 20
  }
  if (can_use_asymp) {
    v <- tryCatch({
      if (zmag > wp) {
        v <- .mpc(.mp_one)
      } else {
        r <- .lm_complex_ei_asymptotic(.hx_fix(a, wp), .hx_fix(b, wp), wp)
        v <- .mpc(.lm(r[[1]], -wp), .lm(r[[2]], -wp))
      }
      v <- .mpc_mul(v, .mpc_exp(z, wp), wp)
      v <- .mpc_div(v, z, wp)
      if (e1) {
        v <- .mpc(.mp_round(-v$re$m, v$re$e, prec), .mp_round(-v$im$m, v$im$e, prec))
      } else {
        v <- .mpc(.mp_pos(v$re, prec),
                  if (b$m < 0) .mp_sub(v$im, .mp_pi(wp), prec) else .mp_add(v$im, .mp_pi(wp), prec))
      }
      v
    }, NoConvergence = function(err) NULL)
    if (!is.null(v)) return(list(v$re, v$im))
  }
  wp <- wp + 2 * as.numeric(.mp_to_int(.mpc_abs(z, 5)))
  r <- .lm_complex_ei_taylor(.hx_fix(a, wp), .hx_fix(b, wp), wp)
  v <- .mpc(.lm(r[[1]] + .fx_euler(wp), -wp), .lm(r[[2]], -wp))
  u <- if (e1) .mpc_log(.mpc_neg(z), wp) else .mpc_log(z, wp)
  v <- .mpc_add(v, u, prec)
  if (e1) v <- .mpc_neg(v)
  list(v$re, v$im)
}

.lm_special <- function(v) if (is.character(v)) .mv_special(v) else .mv_mpf(v)
.lm_cpair <- function(r) {
  if (is.character(r[[1]])) return(.mv_special(r[[1]], TRUE))
  .mv_mpc(r[[1]], r[[2]])
}

# ctx.ei and ctx.e1 (_wrap_libmp_function: a real argument falls back to the
# complex routine on ComplexResult)
.sf_ei <- function(z, e1 = FALSE) {
  z <- .mv_convert(z); p <- .MPX$prec
  if (.mv_sp(z)) {
    if (z$sp == "nan") return(z)
    pos <- z$sp == "inf"
    if (e1) return(if (pos) .hx_zero() else .mv_special("nan"))
    return(if (pos) .mv_inf() else .hx_zero())
  }
  if (z$k == "mpf") {
    r <- .lm_ei(z$re, p, e1)
    if (!is.null(r)) return(.lm_special(r))
    return(.lm_cpair(.lm_ei_c(.mpc(z$re), p, e1)))
  }
  .lm_cpair(.lm_ei_c(.mpc(z$re, z$im), p, e1))
}

# mpf_expint(n, x): E_n(x) for an integer n and real x (gamma = TRUE gives
# the upper incomplete gamma function Gamma(n, x)); NULL for NotImplementedError
.lm_expint <- function(n, x, prec, gamma = FALSE) {
  if (x$m == 0) {
    if (gamma) {
      if (n <= 0) return(list("inf", NULL))
      return(list(.mp_round(gmp::factorialZ(n - 1), 0, prec), NULL))
    }
    if (n > 1) return(list(.mp_from_q(gmp::as.bigq(1, n - 1), prec), NULL))
    return(list("inf", NULL))
  }
  n_orig <- n
  if (gamma) n <- 1 - n
  wp <- prec + 20
  xmag <- .mp_mag(x)
  if (xmag < -10) return(NULL)
  nmag <- .bitlen(abs(gmp::as.bigz(n)))
  have_imag <- n > 0 && x$m < 0
  negx <- .mp_neg(x)
  if (n == 0 || 2 * nmag - xmag < -wp) {
    v <- .mp_exp(negx, wp)
    re <- if (gamma) .mp_mul(v, .mp_pow_int(x, n_orig - 1, wp), prec) else .mp_div(v, x, prec)
  } else {
    can_use_asymptotic_series <- -3 * wp < n && n <= 0
    if (!can_use_asymptotic_series) {
      xi <- as.numeric(abs(.mp_to_int(x)))
      m <- min(max(1, xi - n), 2 * wp)
      siz <- -n * nmag + (m + n) * .bitlen(gmp::as.bigz(m + n)) - m * xmag - (144 * m) %/% 100
      tol <- -wp - 10
      can_use_asymptotic_series <- siz < tol
    }
    if (can_use_asymptotic_series) {
      r <- .fdiv(-.pow2(wp + wp), .hx_fix(x, wp))
      m <- n
      t <- r * m
      s <- .pow2(wp)
      while (m != 0 && t != 0) {
        s <- s + t
        m <- m + 1
        t <- .hx_shr(m * r * t, wp)
      }
      v <- .mp_exp(negx, wp)
      v <- if (gamma) .mp_mul(v, .mp_pow_int(x, n_orig - 1, wp), wp) else .mp_div(v, x, wp)
      re <- .mp_mul(v, .lm(s, -wp), prec)
    } else if (n == 1) {
      re <- .lm_ei(negx, prec)
      if (is.character(re)) return(NULL)
      re <- .mp_neg(re)
    } else if (n > 0 && n < 3 * wp) {
      T1 <- .lm_ei(negx, wp)
      if (is.character(T1)) return(NULL)
      T1 <- .mp_neg(T1)
      if (gamma) {
        if (n_orig %% 2 == 1) T1 <- .mp_neg(T1)
      } else {
        T1 <- .mp_mul(T1, .mp_pow_int(negx, n - 1, wp), wp)
      }
      r <- t <- .hx_fix(x, wp)
      facs <- rep(list(gmp::as.bigz(1)), n - 1)
      if (n > 2) for (k in 1:(n - 2)) facs[[k + 1]] <- facs[[k]] * k
      facs <- rev(facs)
      s <- facs[[1]] * .pow2(wp)
      if (n > 2) for (k in 1:(n - 2)) {
        if (k %% 2 == 1) s <- s - facs[[k + 1]] * t else s <- s + facs[[k + 1]] * t
        t <- .hx_shr(t * r, wp)
      }
      T2 <- .lm(s, -wp, wp)
      T2 <- .mp_mul(T2, .mp_exp(negx, wp), Inf)
      if (gamma) T2 <- .mp_mul(T2, .mp_pow_int(x, n_orig, wp), wp)
      R <- .mp_add(T1, T2, Inf)
      re <- .mp_div(R, .mp_round(gmp::factorialZ(n - 1), 0, Inf), prec)
    } else {
      return(NULL)
    }
  }
  if (have_imag) {
    M <- list(m = -gmp::factorialZ(n - 1), e = 0)
    if (gamma) {
      im <- .mp_div(.mp_pi(wp), M, prec)
      if (n_orig %% 2 == 1) im <- .mp_neg(im)
    } else {
      im <- .mp_div(.mp_mul(.mp_pi(wp), .mp_pow_int(negx, n_orig - 1, wp), wp), M, prec)
    }
    return(list(re, im))
  }
  list(re, NULL)
}
.lm_expint_value <- function(r) {
  if (is.character(r[[1]])) return(.mv_special(r[[1]]))
  if (is.null(r[[2]])) .mv_mpf(r[[1]]) else .mv_mpc(r[[1]], r[[2]])
}

# ctx._gamma_upper_int and ctx._expint_int; NULL for NotImplementedError
.sf_gamma_upper_int <- function(n, z) {
  n <- as.numeric(.hx_int_of(.hx_re(n)))
  if (n == 0) return(.sf_ei(z, e1 = TRUE))
  z <- .mv_convert(z)
  if (z$k != "mpf" || .mv_sp(z)) return(NULL)
  r <- .lm_expint(n, z$re, .MPX$prec, gamma = TRUE)
  if (is.null(r)) NULL else .lm_expint_value(r)
}
.sf_expint_int <- function(n, z) {
  n <- as.numeric(.hx_int_of(n))
  if (n == 1) return(.sf_e1(z))
  z <- .mv_convert(z)
  if (z$k != "mpf" || .mv_sp(z)) return(NULL)
  r <- .lm_expint(n, z$re, .MPX$prec)
  if (is.null(r)) NULL else .lm_expint_value(r)
}

.sf_e1 <- function(z) {
  z <- .mv_convert(z)
  # ctx._e1 never raises NotImplementedError for mpf and mpc arguments
  .sf_ei(z, e1 = TRUE)
}

.sf_expint <- function(n, z) .hx_wrapped(function(n, z) {
  if (.hx_isint(n) && .hx_is_real_type(z)) {
    v <- .sf_expint_int(n, z)
    if (!is.null(v)) return(v)
  }
  if (.hx_isnan(n) || .hx_isnan(z)) return(.hx_mul(z, n))
  if (.mv_sp(z) && z$sp == "inf") return(.hx_div(1, z))
  if (.hx_iszero(z)) {
    if (.mv_cmp(.hx_re(n), 1) <= 0) return(.mv_inf())
    return(.hx_div(.hx_one(), .hx_sub(n, 1)))
  }
  if (.hx_numeq(n, 0)) return(.hx_div(.hx_exp(.hx_neg(z)), z))
  if (.hx_numeq(n, -1)) return(.hx_div(.hx_mul(.hx_exp(.hx_neg(z)), .hx_add(z, 1)), .hx_pow(z, 2)))
  .hx_mul(.hx_pow(z, .hx_sub(n, 1)), .sf_upper_gamma(.hx_sub(1, n), z))
}, n, z)

.sf_li <- function(z) .hx_wrapped(function(z) {
  if (.hx_iszero(z)) return(z)
  if (.hx_numeq(z, 1)) return(.mv_ninf())
  .sf_ei(.hx_log(z))
}, z)

# Ci and Si (libmp/libhyper.py)
.lm_ci_si_taylor <- function(x, wp, which = 0) {
  x <- .hx_fix(x, wp)
  x2 <- .hx_shr(-(x * x), wp)
  if (which == 0) { s <- gmp::as.bigz(0); t <- .pow2(wp); k <- 2 }
  else { s <- x; t <- x; k <- 3 }
  while (t != 0) {
    t <- .hx_shr(.fdiv(t * x2, k * (k - 1)), wp)
    s <- s + .fdiv(t, k)
    k <- k + 2
  }
  .lm(s, -wp)
}
.lm_mpc_ci_si_taylor <- function(re, im, wp, which = 0) {
  if (re$m != 0) mag <- .mp_mag(re) else if (im$m != 0) mag <- .mp_mag(im)
  if (im$m != 0) mag <- max(mag, .mp_mag(im))
  if (mag > 2 || mag < -wp) return(NULL)
  wp <- wp + (2 - mag)
  zre <- .hx_fix(re, wp); zim <- .hx_fix(im, wp)
  z2re <- .hx_shr(zim * zim - zre * zre, wp)
  z2im <- .hx_shr(-2 * zre * zim, wp)
  if (which == 0) { sre <- gmp::as.bigz(0); sim <- gmp::as.bigz(0); tre <- .pow2(wp); tim <- gmp::as.bigz(0); k <- 2 }
  else { sre <- zre; sim <- zim; tre <- zre; tim <- zim; k <- 3 }
  while (max(abs(tre), abs(tim)) > 2) {
    f <- k * (k - 1)
    a <- .hx_shr(.fdiv(tre * z2re - tim * z2im, f), wp)
    b <- .hx_shr(.fdiv(tre * z2im + tim * z2re, f), wp)
    tre <- a; tim <- b
    sre <- sre + .fdiv(tre, k); sim <- sim + .fdiv(tim, k)
    k <- k + 2
  }
  .mpc(.lm(sre, -wp), .lm(sim, -wp))
}
.lm_ci_si <- function(x, prec, which = 2) {
  wp <- prec + 20
  ci <- NULL; si <- NULL
  if (x$m == 0) return(list(ci = "-inf", si = .mp_zero))
  mag <- .mp_mag(x)
  if (mag < -wp) {
    if (which != 0) si <- .mp_pos(x, prec)
    if (which != 1) ci <- .mp_add(.mp_euler(wp), .mp_log(list(m = abs(x$m), e = x$e), wp), prec)
    return(list(ci = ci, si = si))
  } else if (mag > wp) {
    if (which != 0) si <- .mp_ldexp(if (x$m < 0) .mp_neg(.mp_pi(prec)) else .mp_pi(prec), -1)
    if (which != 1) ci <- .mp_div(.mp_cos_sin(x, wp)$sin, x, prec)
    return(list(ci = ci, si = si))
  } else {
    wp <- wp + abs(mag)
  }
  asymptotic <- mag - 1 > log2(wp)
  if (!asymptotic) {
    if (which != 0) si <- .mp_pos(.lm_ci_si_taylor(x, wp, 1), prec)
    if (which != 1) {
      ci <- .lm_ci_si_taylor(x, wp, 0)
      ci <- .mp_add(ci, .mp_euler(wp), wp)
      ci <- .mp_add(ci, .mp_log(list(m = abs(x$m), e = x$e), wp), prec)
    }
    return(list(ci = ci, si = si))
  }
  neg <- x$m < 0
  x <- list(m = abs(x$m), e = x$e)
  xf <- .hx_fix(x, wp)
  xr <- .fdiv(.pow2(2 * wp), xf)
  s1 <- .pow2(wp); s2 <- xr; t <- xr; k <- 2
  while (t != 0) {
    t <- -t
    t <- .hx_shr(t * xr * k, wp); k <- k + 1
    s1 <- s1 + t
    t <- .hx_shr(t * xr * k, wp); k <- k + 1
    s2 <- s2 + t
  }
  s1 <- .mp_div(.lm(s1, -wp), x, wp)
  s2 <- .mp_div(.lm(s2, -wp), x, wp)
  cs <- .mp_cos_sin(x, wp)
  if (which != 0) {
    si <- .mp_add(.mp_mul(cs$cos, s1, Inf), .mp_mul(cs$sin, s2, Inf), wp)
    si <- .mp_sub(.mp_ldexp(.mp_pi(wp), -1), si, wp)
    if (neg) si <- .mp_neg(si)
    si <- .mp_pos(si, prec)
  }
  if (which != 1) ci <- .mp_sub(.mp_mul(cs$sin, s1, Inf), .mp_mul(cs$cos, s2, Inf), prec)
  list(ci = ci, si = si)
}

# ctx.ci, ctx.si: the libmp routines, else the generic formulas
.sf_ci <- function(z) {
  z <- .mv_convert(z); p <- .MPX$prec
  r <- NULL
  if (!.mv_sp(z)) {
    if (z$k == "mpf" || z$im$m == 0) {
      v <- .lm_ci_si(z$re, p, 0)$ci
      v <- if (is.character(v)) .mv_special(v) else .mv_mpf(v)
      if (z$re$m < 0) return(.mv_mpc(v$re, .mp_pi(p)))
      return(if (z$k == "mpc") .mv_mpc(v$re, .mp_zero) else v)
    }
    wp <- p + 20
    t <- .lm_mpc_ci_si_taylor(z$re, z$im, wp, 0)
    if (!is.null(t)) {
      t <- .mpc(.mp_add(t$re, .mp_euler(wp), wp), t$im)
      ci <- .mpc_add(t, .mpc_log(.mpc(z$re, z$im), wp), p)
      return(.mv_mpc(ci$re, ci$im))
    }
  }
  .sf_ci_generic(z)
}
.sf_ci_generic <- function(z) .hx_wrapped(function(z) {
  if (.hx_isinf(z)) return(if (z$sp == "inf") .hx_zero() else .hx_mul(.hx_pi(), .mv_mpc(.mp_zero, .mpf(1, 0))))
  jz <- .hx_fmul(.hx_j(), z, 0)
  njz <- .hx_neg_exact(jz)
  v <- .hx_mul(0.5, .hx_add(.sf_ei(jz), .sf_ei(njz)))
  zreal <- .hx_re(z); zimag <- .hx_im(z)
  pij <- function(f) .hx_mul(.hx_pi(), .mv_mpc(.mp_zero, .mp_from_dbl(f)))
  if (.hx_iszero(zreal)) {
    if (.mv_cmp(zimag, 0) > 0) v <- .hx_add(v, pij(0.5))
    if (.mv_cmp(zimag, 0) < 0) v <- .hx_sub(v, pij(0.5))
  }
  if (.mv_cmp(zreal, 0) < 0) {
    if (.mv_cmp(zimag, 0) >= 0) v <- .hx_add(v, pij(1))
    if (.mv_cmp(zimag, 0) < 0) v <- .hx_sub(v, pij(1))
  }
  if (.hx_is_real_type(z) && .mv_cmp(zreal, 0) > 0) v <- .hx_re(v)
  v
}, z)

.sf_si <- function(z) {
  z <- .mv_convert(z); p <- .MPX$prec
  if (!.mv_sp(z)) {
    if (z$k == "mpf") return(.mv_mpf(.lm_ci_si(z$re, p, 1)$si))
    if (z$im$m == 0) return(.mv_mpc(.lm_ci_si(z$re, p, 1)$si, .mp_zero))
    t <- .lm_mpc_ci_si_taylor(z$re, z$im, p + 20, 1)
    if (!is.null(t)) return(.mv_mpc(.mp_pos(t$re, p), .mp_pos(t$im, p)))
  }
  .sf_si_generic(z)
}
.sf_si_generic <- function(z) .hx_wrapped(function(z) {
  if (.hx_isinf(z)) return(.hx_mul(if (z$sp == "inf") 0.5 else -0.5, .hx_pi()))
  if (.hx_mag(z) >= -1) {
    jz <- .hx_fmul(.hx_j(), z, 0)
    njz <- .hx_neg_exact(jz)
    v <- .hx_mul(.mv_mpc(.mp_zero, .mpf(-1, -1)), .hx_sub(.sf_ei(jz), .sf_ei(njz)))
    zreal <- .hx_re(z)
    if (.mv_cmp(zreal, 0) > 0) v <- .hx_sub(v, .hx_mul(0.5, .hx_pi()))
    if (.mv_cmp(zreal, 0) < 0) v <- .hx_add(v, .hx_mul(0.5, .hx_pi()))
    if (.hx_is_real_type(z)) v <- .hx_re(v)
    return(v)
  }
  .hx_mul(z, .hx_hyper(list(.mv_tuple(1, 2)), list(.mv_tuple(3, 2), .mv_tuple(3, 2)), .hx_mul(.hx_mul(-0.25, z), z)))
}, z)

.sf_chi <- function(z) .hx_wrapped(function(z) {
  nz <- .hx_fmul(z, -1, 0)
  v <- .hx_mul(0.5, .hx_add(.sf_ei(z), .sf_ei(nz)))
  zreal <- .hx_re(z); zimag <- .hx_im(z)
  pij <- function(f) .hx_mul(.hx_pi(), .mv_mpc(.mp_zero, .mp_from_dbl(f)))
  if (.mv_cmp(zimag, 0) > 0) v <- .hx_add(v, pij(0.5))
  else if (.mv_cmp(zimag, 0) < 0) v <- .hx_sub(v, pij(0.5))
  else if (.mv_cmp(zreal, 0) < 0) v <- .hx_add(v, pij(1))
  v
}, z)

.sf_shi <- function(z) .hx_wrapped(function(z) {
  if (.hx_mag(z) >= -1) {
    nz <- .hx_fmul(z, -1, 0)
    v <- .hx_mul(0.5, .hx_sub(.sf_ei(z), .sf_ei(nz)))
    zimag <- .hx_im(z)
    if (.mv_cmp(zimag, 0) > 0) v <- .hx_sub(v, .hx_mul(.mv_mpc(.mp_zero, .mpf(1, -1)), .hx_pi()))
    if (.mv_cmp(zimag, 0) < 0) v <- .hx_add(v, .hx_mul(.mv_mpc(.mp_zero, .mpf(1, -1)), .hx_pi()))
    return(v)
  }
  .hx_mul(z, .hx_hyper(list(.mv_tuple(1, 2)), list(.mv_tuple(3, 2), .mv_tuple(3, 2)), .hx_mul(.hx_mul(0.25, z), z)))
}, z)

.sf_fresnel <- function(z, which) .hx_wrapped(function(z) {
  if (.mv_sp(z) && z$sp %in% c("inf", "-inf")) return(.mv_mpf(.mpf(if (z$sp == "inf") 1 else -1, -1)))
  w <- .hx_div(.hx_mul(.hx_neg(.hx_pow(.hx_pi(), 2)), .hx_pow(z, 4)), 16)
  if (which == "s")
    return(.hx_mul(.hx_div(.hx_mul(.hx_pi(), .hx_pow(z, 3)), 6),
                   .hx_hyper(list(.mv_tuple(3, 4)), list(.mv_tuple(3, 2), .mv_tuple(7, 4)), w)))
  .hx_mul(z, .hx_hyper(list(.mv_tuple(1, 4)), list(.mv_tuple(1, 2), .mv_tuple(5, 4)), w))
}, z)

# ---------------------------------------------------------------------------
# Incomplete gamma functions (functions/expintegrals.py)

.sf_gammainc <- function(z, a = 0, b = NULL, regularized = FALSE) {
  z <- .mv_convert(z)
  if (is.null(a)) { a <- .hx_zero(); lower_modified <- FALSE }
  else { a <- .mv_convert(a); lower_modified <- !.hx_iszero(a) }
  if (is.null(b)) { b <- .mv_inf(); upper_modified <- FALSE }
  else { b <- .mv_convert(b); upper_modified <- !(.mv_sp(b) && b$sp == "inf") }
  if (!(upper_modified || lower_modified)) {
    if (regularized) {
      r <- .mv_cmp(.hx_re(z), 0)
      return(if (r < 0) .mv_inf() else if (r > 0) .hx_one() else .mv_nan())
    }
    return(.hx_gamma(z))
  }
  if (.hx_numeq(a, b)) return(.hx_zero())
  if (.mv_cmp(.hx_re(a), .hx_re(b)) > 0) return(.hx_neg(.sf_gammainc(z, b, a, regularized)))
  if (upper_modified && lower_modified) return(.mv_pos(.sf_gamma3(z, a, b, regularized)))
  if (lower_modified) return(.sf_upper_gamma(z, a, regularized))
  .sf_lower_gamma(z, b, regularized)
}

.sf_lower_gamma <- function(z, b, regularized = FALSE) {
  z <- .mv_convert(z); b <- .mv_convert(b)
  if (.hx_isnpint(z)) return(.mv_inf())
  G <- if (regularized) list(z) else list()
  negb <- .hx_neg_exact(b)
  h <- function(z) list(.hx_term(list(.hx_exp(negb), b, z), list(1, z, -1), list(), G, list(1), list(.hx_add(1, z)), b))
  .hx_hypercomb(h, list(z))
}

.sf_upper_gamma <- function(z, a, regularized = FALSE) {
  z <- .mv_convert(z); a <- .mv_convert(a)
  if (.hx_isint(z)) {
    if (regularized) {
      if (.hx_isnpint(z)) return(.hx_zero())
      orig <- .MPX$prec
      .MPX$prec <- .MPX$prec + 10
      v <- .sf_gamma_upper_int(z, a)
      if (!is.null(v)) { v <- .hx_div(v, .hx_gamma(z)); .MPX$prec <- orig; return(v) }
      .MPX$prec <- orig
    } else {
      v <- .sf_gamma_upper_int(z, a)
      if (!is.null(v)) return(v)
    }
  }
  if (.hx_numeq(z, 2) && .hx_numeq(a, -1)) return(.hx_mul(.hx_add(z, a), 0))
  if (.hx_numeq(z, 3) && (.hx_numeq(a, .mv_mpc(.mpf(-1, 0), .mpf(-1, 0))) || .hx_numeq(a, .mv_mpc(.mpf(-1, 0), .mpf(1, 0)))))
    return(.hx_mul(.hx_add(z, a), 0))
  nega <- .hx_neg_exact(a)
  G <- if (regularized) list(z) else list()
  v <- tryCatch({
    h <- function(z) {
      r <- .hx_sub(z, 1)
      list(.hx_term(list(.hx_exp(nega), a), list(1, r), list(), G, list(1, .hx_neg(r)), list(), .hx_div(1, nega)))
    }
    .hx_hypercomb(h, list(z), kw = list(force_series = TRUE))
  }, NoConvergence = function(err) NULL)
  if (!is.null(v)) return(v)
  h <- function(z) list(.hx_term(list(), list(1, .hx_sub(z, 1)), list(z), G, list(), list(), 0),
                        .hx_term(list(.hx_neg(.hx_exp(nega)), a, z), list(1, z, -1), list(), G, list(1), list(.hx_add(1, z)), a))
  .hx_hypercomb(h, list(z))
}

.sf_gamma3 <- function(z, a, b, regularized = FALSE) {
  pole <- .hx_isnpint(z)
  if (regularized && pole) return(.hx_zero())
  .MPX$prec <- .MPX$prec + 15
  on.exit(.MPX$prec <- .MPX$prec - 15)
  T1 <- .sf_gammainc(z, a, regularized = regularized)
  T2 <- .sf_gammainc(z, b, regularized = regularized)
  R <- .hx_sub(T1, T2)
  if (.hx_mag(R) - max(.hx_mag(T1), .hx_mag(T2)) > -10) return(R)
  if (!pole) {
    T1 <- .sf_gammainc(z, 0, b, regularized = regularized)
    T2 <- .sf_gammainc(z, 0, a, regularized = regularized)
    return(.hx_sub(T1, T2))
  }
  .mp_error("NotImplementedError", "")
}

# ---------------------------------------------------------------------------
# Orthogonal polynomials at general degree (functions/orthogonal.py)

.sf_hyp2f1 <- function(a, b, c, z, kw = list()) .hx_hyper(list(a, b), list(c), z, kw)

.sf_legendre <- function(n, x) .hx_wrapped(function(n, x) {
  if (.hx_isint(n)) {
    n <- .hx_int_of(n)
    if (gmp::mod.bigz(n + (n < 0), 2) == 1) {
      if (.hx_iszero(x)) return(x)
      mag <- .hx_mag(x)
      if (mag < -2 * .MPX$prec - 10) return(x)
      if (mag < -5) .MPX$prec <- .MPX$prec - mag
    }
    n <- .mv_int(n)
  }
  .sf_hyp2f1(.hx_neg(n), .hx_add(n, 1), 1, .hx_div(.hx_sub(1, x), 2))
}, n, x)

.sf_legenp <- function(n, m, z) {
  n <- .mv_convert(n); m <- .mv_convert(m)
  if (.hx_iszero(m)) return(.sf_legendre(n, z))
  h <- function(n, m) {
    g <- .hx_mul(m, 0.5)
    list(.hx_term(list(.hx_add(1, z), .hx_sub(1, z)), list(g, .hx_neg(g)), list(), list(.hx_sub(1, m)),
                  list(.hx_neg(n), .hx_add(n, 1)), list(.hx_sub(1, m)), .hx_mul(0.5, .hx_sub(1, z))))
  }
  .hx_hypercomb(h, list(n, m))
}

.sf_chebyt <- function(n, x) .hx_wrapped(function(n, x) {
  if (.hx_iszero(x) && .hx_isint(n) && gmp::mod.bigz(.hx_int_of(.hx_re(n)), 2) == 1) return(.hx_mul(x, 0))
  .sf_hyp2f1(.hx_neg(n), n, .mv_tuple(1, 2), .hx_div(.hx_sub(1, x), 2), list(force_series = TRUE))
}, n, x)

.sf_chebyu <- function(n, x) .hx_wrapped(function(n, x) {
  if (.hx_iszero(x) && .hx_isint(n) && gmp::mod.bigz(.hx_int_of(.hx_re(n)), 2) == 1) return(.hx_mul(x, 0))
  .hx_mul(.hx_add(n, 1), .sf_hyp2f1(.hx_neg(n), .hx_add(n, 2), .mv_tuple(3, 2), .hx_div(.hx_sub(1, x), 2), list(force_series = TRUE)))
}, n, x)

.sf_hermite <- function(n, z) .hx_hypercomb(function() {
  pn <- .hx_param(n); n <- pn$v
  z <- .mv_convert(z)
  q <- .mv_q(-1, 2)
  if (.hx_iszero(z))
    return(list(.hx_term(list(2, .hx_pi()), list(n, 0.5), list(), list(.hx_mul(q, .hx_sub(n, 1))), list(), list(), 0)))
  can_use_2f0 <- .hx_isnpint(.hx_neg(n)) || .mv_cmp(.hx_re(z), 0) > 0 ||
    (.hx_iszero(.hx_re(z)) && .mv_cmp(.hx_im(z), 0) > 0)
  expprec <- .MPX$prec * 4 + 20
  w <- z
  w2 <- .hx_fmul(w, w, expprec)
  old <- .MPX$prec; .MPX$prec <- expprec
  rw2 <- .hx_div(1, w2)
  .MPX$prec <- old
  nrw2 <- .hx_neg_exact(rw2)
  nw <- .hx_neg_exact(w)
  if (can_use_2f0)
    return(list(.hx_term(list(2, w), list(n, n), list(), list(), list(.hx_mul(q, n), .hx_mul(q, .hx_sub(n, 1))), list(), nrw2)))
  list(.hx_term(list(2, nw), list(n, n), list(), list(), list(.hx_mul(q, n), .hx_mul(q, .hx_sub(n, 1))), list(), nrw2),
       .hx_term(list(2, .hx_pi(), nw), list(.hx_add(n, 2), 0.5, 1), list(), list(.hx_mul(q, n)),
                list(.hx_mul(q, .hx_sub(n, 1))), list(.hx_sub(1, q)), w2))
}, list())

.sf_gegenbauer <- function(n, a, z) .hx_wrapped(function(n, a, z) {
  if (.hx_isnpint(a)) return(.hx_mul(0, .hx_add(z, n)))
  term <- function(n, a) {
    a2 <- .hx_mul(2, a)
    list(.hx_term(list(), list(), list(.hx_add(n, a2)), list(.hx_add(n, 1), a2), list(.hx_neg(n), .hx_add(n, a2)),
                  list(.hx_add(a, 0.5)), .hx_mul(0.5, .hx_sub(1, z))))
  }
  if (.hx_isnpint(.hx_add(a, 0.5))) {
    if (.hx_isnpint(.hx_add(n, 1))) .mp_error("NotImplementedError", "Gegenbauer function with two limits")
    return(.hx_hypercomb(function(a) term(n, a), list(a)))
  }
  .hx_hypercomb(function(n) term(n, a), list(n))
}, n, a, z)

.sf_jacobi <- function(n, a, b, x) .hx_wrapped(function(n, a, b, x) {
  if (!.hx_isnpint(a)) {
    h <- function(n) list(.hx_term(list(), list(), list(.hx_add(.hx_add(a, n), 1)), list(.hx_add(n, 1), .hx_add(a, 1)),
                                   list(.hx_neg(n), .hx_add(.hx_add(.hx_add(a, b), n), 1)), list(.hx_add(a, 1)),
                                   .hx_mul(.hx_sub(1, x), 0.5)))
    return(.hx_hypercomb(h, list(n)))
  }
  if (!.hx_isint(b)) {
    h2 <- function(n, a) list(.hx_term(list(), list(), list(.hx_neg(b)), list(.hx_add(n, 1), .hx_sub(.hx_neg(b), n)),
                                      list(.hx_neg(n), .hx_add(.hx_add(.hx_add(a, b), n), 1)), list(.hx_add(b, 1)),
                                      .hx_mul(.hx_add(x, 1), 0.5)))
    return(.hx_hypercomb(h2, list(n, a)))
  }
  .mp_error("NotImplementedError", "jacobi with a nonpositive integer a and an integer b")
}, n, a, b, x)

# ---------------------------------------------------------------------------
# Elliptic integrals (libmp/libhyper.py and functions/elliptic.py)

# agm_fixed: a, b fixed point numbers close to unit magnitude
.lm_agm_fixed <- function(a, b, prec) {
  i <- 0
  repeat {
    anew <- .hx_shr(a + b, 1)
    if (i > 4 && abs(a - anew) < 8) return(a)
    b <- .mp_isqrt(a * b)
    a <- anew
    i <- i + 1
  }
}
.lm_agm <- function(a, b, prec) {
  if (a$m == 0 || b$m == 0) return(.mp_zero)
  wp <- prec + 20
  amag <- .mp_mag(a); bmag <- .mp_mag(b)
  abs_mag_delta <- abs(amag - bmag)
  if (abs_mag_delta > 10) {
    while (abs_mag_delta > 10) {
      a1 <- .mp_ldexp(.mp_add(a, b, wp), -1)
      b <- .mp_sqrt(.mp_mul(a, b, wp), wp)
      a <- a1
      abs_mag_delta <- abs_mag_delta %/% 2
    }
    amag <- .mp_mag(a); bmag <- .mp_mag(b)
  }
  min_mag <- min(amag, bmag); max_mag <- max(amag, bmag)
  n <- 0
  if (min_mag < -8) n <- -min_mag else if (max_mag > 20) n <- -max_mag
  if (n != 0) { a <- .mp_ldexp(a, n); b <- .mp_ldexp(b, n) }
  g <- .lm_agm_fixed(.hx_fix(a, wp), .hx_fix(b, wp), wp)
  .lm(g, -wp - n, prec)
}
.lm_mpc_agm <- function(a, b, prec) {
  if ((a$re$m == 0 && a$im$m == 0) || (b$re$m == 0 && b$im$m == 0)) return(.mpc(.mp_zero))
  wp <- prec + 20
  eps <- .mpf(1, -wp + 10)
  repeat {
    s <- .mpc_add(a, b, wp)
    a1 <- .mpc(.mp_ldexp(s$re, -1), .mp_ldexp(s$im, -1))
    b1 <- .mpc_sqrt(.mpc_mul(a, b, wp), wp)
    a <- a1; b <- b1
    size <- max(.mp_to_dbl(.mpc_abs(a, 10)), .mp_to_dbl(.mpc_abs(b, 10)))
    sz <- if (.mp_to_dbl(.mpc_abs(a, 10)) >= .mp_to_dbl(.mpc_abs(b, 10))) .mpc_abs(a, 10) else .mpc_abs(b, 10)
    err <- .mpc_abs(.mpc_sub(a, b, 10), 10)
    if (sz$m == 0 || .mp_cmp(err, .mp_mul(eps, sz, Inf)) < 0) return(a)
  }
}
.lm_ellipk <- function(x, prec) {
  if (x$m == 0) return(.mp_ldexp(.mp_pi(prec), -1))
  if (.mp_cmp(x, .mp_one) == 0) return("inf")
  wp <- prec + 15
  a <- .mp_sqrt(.mp_sub(.mp_one, x, wp), wp)
  v <- .lm_agm(.mp_one, a, wp)
  .mp_ldexp(.mp_div(.mp_pi(wp), v, prec), -1)
}
.lm_mpc_ellipk <- function(z, prec) {
  if (z$im$m == 0 && .mp_cmp(z$re, .mp_one) <= 0) return(.mpc(.lm_ellipk(z$re, prec)))
  wp <- prec + 15
  a <- .mpc_sqrt(.mpc_sub(.mpc_one, z, wp), wp)
  v <- .lm_mpc_agm(.mpc_one, a, wp)
  r <- .mpc_div(.mpc(.mp_pi(wp)), v, prec)
  .mpc(.mp_ldexp(r$re, -1), .mp_ldexp(r$im, -1))
}
.lm_ellipe <- function(x, prec) {
  if (x$m == 0) return(.mp_ldexp(.mp_pi(prec), -1))
  if (.mp_cmp(x, .mp_one) == 0) return(.mp_one)
  wp <- prec + 20
  mag <- .mp_mag(x)
  if (mag < -wp) return(.mp_ldexp(.mp_pi(prec), -1))
  p <- max(mag, 0) - wp
  h <- .mpf(1, p)
  K <- .lm_ellipk(x, 2 * wp)
  Kh <- .lm_ellipk(.mp_sub(x, h, Inf), 2 * wp)
  Kdiff <- .mp_ldexp(.mp_sub(K, Kh, Inf), -p)
  t <- .mp_sub(.mp_one, x, Inf)
  b <- .mp_mul(Kdiff, .mp_ldexp(x, 1), wp)
  .mp_mul(t, .mp_add(K, b, Inf), prec)
}
.lm_mpc_ellipe <- function(z, prec) {
  if (z$im$m == 0 && .mp_cmp(z$re, .mp_one) <= 0) return(.mpc(.lm_ellipe(z$re, prec)))
  wp <- prec + 15
  mag <- .mpc_abs(z, 1)
  p <- max(.mp_mag(mag), 0) - wp
  h <- .mpf(1, p)
  K <- .lm_mpc_ellipk(z, 2 * wp)
  Kh <- .lm_mpc_ellipk(.mpc_add(z, .mpc(h), 2 * wp), 2 * wp)
  Kdiff <- .mpc_sub(Kh, K, wp)
  Kdiff <- .mpc(.mp_ldexp(Kdiff$re, -p), .mp_ldexp(Kdiff$im, -p))
  t <- .mpc_sub(.mpc_one, z, wp)
  b <- .mpc_mul(Kdiff, .mpc(.mp_ldexp(z$re, 1), .mp_ldexp(z$im, 1)), wp)
  .mpc_mul(t, .mpc_add(K, b, wp), prec)
}
# ctx.ellipk and ctx._ellipe: the real routine, the complex one on ComplexResult
.sf_ellipk <- function(m) {
  m <- .mv_convert(m); p <- .MPX$prec
  if (m$k == "mpf" && .mp_cmp(m$re, .mp_one) <= 0) {
    v <- .lm_ellipk(m$re, p)
    return(if (is.character(v)) .mv_special(v) else .mv_mpf(v))
  }
  if (m$k == "mpf" && .mp_cmp(m$re, .mp_one) == 0) return(.mv_inf())
  r <- .lm_mpc_ellipk(.mpc(m$re, if (m$k == "mpc") m$im else .mp_zero), p)
  .mv_mpc(r$re, r$im)
}
.sf_ellipe1 <- function(m) {
  m <- .mv_convert(m); p <- .MPX$prec
  if (m$k == "mpf" && .mp_cmp(m$re, .mp_one) <= 0) return(.mv_mpf(.lm_ellipe(m$re, p)))
  r <- .lm_mpc_ellipe(.mpc(m$re, if (m$k == "mpc") m$im else .mp_zero), p)
  .mv_mpc(r$re, r$im)
}

# ctx.root(x, n) for n = -6 (RF_calc, RJ_calc)
.hx_root <- function(x, n) {
  if (n > 0) return(.hx_nthroot(x, n))
  old <- .MPX$prec; .MPX$prec <- old + 5
  r <- .hx_nthroot(x, -n)
  .MPX$prec <- old
  .hx_div(1, r)
}
.hx_acos_r <- function(a) .mv_from_mpc(.mpc_fun("acos", .mv_to_mpc(a), .MPX$prec))
.hx_acosh_r <- function(a) .mv_from_mpc(.mpc_fun("acosh", .mv_to_mpc(a), .MPX$prec))
.mv_to_mpc <- function(a) { a <- .mv_convert(a); .mpc(a$re, if (a$k == "mpc") a$im else .mp_zero) }
.mv_from_mpc <- function(z) if (z$im$m == 0) .mv_mpf(z$re) else .mv_mpc(z$re, z$im)
.hx_maxabs <- function(xs) { m <- NULL; for (x in xs) { a <- .hx_abs(x); if (is.null(m) || .mv_cmp(a, m) > 0) m <- a }; m }

.sf_rc_calc <- function(x, y, r, pv = TRUE) {
  if (pv && .hx_iszero(.hx_im(y)) && .mv_cmp(.hx_re(y), 0) < 0)
    return(.hx_mul(.hx_sqrt(.hx_div(x, .hx_sub(x, y))), .sf_rc_calc(.hx_sub(x, y), .hx_neg(y), r)))
  if (.hx_numeq(x, y)) return(.hx_div(1, .hx_sqrt(x)))
  extraprec <- 2 * max(0, -.hx_mag(.hx_sub(x, y)) + .hx_mag(x))
  .MPX$prec <- .MPX$prec + extraprec
  if (.hx_is_real_type(x) && .hx_is_real_type(y)) {
    x <- .hx_re(x); y <- .hx_re(y)
    a <- .hx_sqrt(.hx_div(x, y))
    if (.mv_cmp(x, y) < 0) v <- .hx_div(.hx_acos_r(a), .hx_sqrt(.hx_sub(y, x)))
    else v <- .hx_div(.hx_acosh_r(a), .hx_sqrt(.hx_sub(x, y)))
  } else {
    sx <- .hx_sqrt(x); sy <- .hx_sqrt(y)
    v <- .hx_div(.hx_acos_r(.hx_div(sx, sy)), .hx_mul(.hx_sqrt(.hx_sub(1, .hx_div(x, y))), sy))
  }
  .MPX$prec <- .MPX$prec - extraprec
  v
}
.sf_rf_calc <- function(x, y, z, r) {
  if (.hx_numeq(y, z)) return(.sf_rc_calc(x, y, r))
  if (.hx_numeq(x, z)) return(.sf_rc_calc(y, x, r))
  if (.hx_numeq(x, y)) return(.sf_rc_calc(z, x, r))
  xm <- x; ym <- y; zm <- z
  A0 <- Am <- .hx_div(.hx_add(.hx_add(x, y), z), 3)
  Q <- .hx_mul(.hx_root(.hx_mul(3, r), -6), .hx_maxabs(list(.hx_sub(A0, x), .hx_sub(A0, y), .hx_sub(A0, z))))
  g <- .mv_mpf(.mpf(1, -2)); pow4 <- .hx_one()
  repeat {
    xs <- .hx_sqrt(xm); ys <- .hx_sqrt(ym); zs <- .hx_sqrt(zm)
    lm <- .hx_add(.hx_add(.hx_mul(xs, ys), .hx_mul(xs, zs)), .hx_mul(ys, zs))
    Am1 <- .hx_mul(.hx_add(Am, lm), g)
    xm <- .hx_mul(.hx_add(xm, lm), g); ym <- .hx_mul(.hx_add(ym, lm), g); zm <- .hx_mul(.hx_add(zm, lm), g)
    if (.mv_cmp(.hx_mul(pow4, Q), .hx_abs(Am)) < 0) break
    Am <- Am1
    pow4 <- .hx_mul(pow4, g)
  }
  t <- .hx_div(pow4, Am)
  X <- .hx_mul(.hx_sub(A0, x), t); Y <- .hx_mul(.hx_sub(A0, y), t)
  Z <- .hx_sub(.hx_neg(X), Y)
  E2 <- .hx_sub(.hx_mul(X, Y), .hx_pow(Z, 2)); E3 <- .hx_mul(.hx_mul(X, Y), Z)
  s <- .hx_sub(.hx_add(.hx_add(.hx_sub(9240, .hx_mul(924, E2)), .hx_mul(385, .hx_pow(E2, 2))), .hx_mul(660, E3)),
               .hx_mul(.hx_mul(630, E2), E3))
  .hx_div(.hx_mul(.hx_power(Am, -0.5), s), 9240)
}
.sf_rj_calc <- function(x, y, z, p, r) {
  if (.hx_iszero(p)) return(.mv_inf())
  if (.hx_iszero(x) + .hx_iszero(y) + .hx_iszero(z) > 1) return(.mv_inf())
  re <- function(t) .hx_re(t); im <- function(t) .hx_im(t)
  ok <- .mv_cmp(re(x), 0) >= 0 && .mv_cmp(re(y), 0) >= 0 && .mv_cmp(re(z), 0) >= 0 && .mv_cmp(re(p), 0) > 0
  if (!ok && (.hx_numeq(x, p) || .hx_numeq(y, p) || .hx_numeq(z, p))) ok <- TRUE
  if (!ok && (!.hx_iszero(im(p)) || .mv_cmp(re(p), 0) >= 0)) {
    rr <- function(a, b, c) .hx_iszero(im(a)) && .mv_cmp(re(a), 0) >= 0 && .hx_numeq(.hx_conj(b), c)
    if (rr(x, y, z) || rr(y, x, z) || rr(z, x, y)) ok <- TRUE
  }
  # outside the region: integrate from 0 to N and move the arguments by N
  initial_integral <- .hx_zero()
  if (!ok) {
    lo <- re(x)
    for (t in list(y, z, p)) if (.mv_cmp(re(t), lo) < 0) lo <- re(t)
    N <- .hx_add(.mv_mpf(.mp_round(-.mp_floor(.mv_convert(lo)$re), 0, Inf)), 1)
    args4 <- list(x, y, z, p)
    upper <- function(t) .mv_cmp(im(t), 0) >= 0 || .mv_cmp(re(t), 0) > 0
    if (all(vapply(args4, upper, TRUE))) {
      margin <- .hx_j()
    } else if (all(vapply(args4, function(t) .mv_cmp(im(t), 0) < 0 || .mv_cmp(re(t), 0) > 0, TRUE))) {
      margin <- .hx_neg(.hx_j())
    } else {
      # through the upper half plane, below any parameter in the lower one
      margin <- .mv_int(1)
      for (t in args4) {
        if (upper(t)) next
        cand <- .hx_mul(.hx_abs(im(t)), 0.5)
        if (.mv_cmp(cand, margin) < 0) margin <- cand
      }
      margin <- .hx_mul(margin, .hx_j())
    }
    N <- .hx_add(N, margin)
    F <- function(t) .hx_div(1, .hx_mul(.hx_mul(.hx_mul(.hx_sqrt(.hx_add(t, x)), .hx_sqrt(.hx_add(t, y))),
                                                .hx_sqrt(.hx_add(t, z))), .hx_add(t, p)))
    initial_integral <- .hx_mul(1.5, .mp_quadsubdiv(F, list(.hx_zero(), N)))
    x <- .hx_add(x, N); y <- .hx_add(y, N); z <- .hx_add(z, N); p <- .hx_add(p, N)
  }
  xm <- x; ym <- y; zm <- z; pm <- p
  A0 <- Am <- .hx_div(.hx_add(.hx_add(.hx_add(x, y), z), .hx_mul(2, p)), 5)
  delta <- .hx_mul(.hx_mul(.hx_sub(p, x), .hx_sub(p, y)), .hx_sub(p, z))
  Q <- .hx_mul(.hx_root(.hx_mul(0.25, r), -6), .hx_maxabs(list(.hx_sub(A0, x), .hx_sub(A0, y), .hx_sub(A0, z), .hx_sub(A0, p))))
  g <- .mv_mpf(.mpf(1, -2)); pow4 <- .hx_one(); S <- 0
  repeat {
    sx <- .hx_sqrt(xm); sy <- .hx_sqrt(ym); sz <- .hx_sqrt(zm); sp <- .hx_sqrt(pm)
    lm <- .hx_add(.hx_add(.hx_mul(sx, sy), .hx_mul(sx, sz)), .hx_mul(sy, sz))
    Am1 <- .hx_mul(.hx_add(Am, lm), g)
    xm <- .hx_mul(.hx_add(xm, lm), g); ym <- .hx_mul(.hx_add(ym, lm), g)
    zm <- .hx_mul(.hx_add(zm, lm), g); pm <- .hx_mul(.hx_add(pm, lm), g)
    dm <- .hx_mul(.hx_mul(.hx_add(sp, sx), .hx_add(sp, sy)), .hx_add(sp, sz))
    em <- .hx_div(.hx_mul(delta, .hx_pow(pow4, 3)), .hx_pow(dm, 2))
    if (.mv_cmp(.hx_mul(pow4, Q), .hx_abs(Am)) < 0) break
    T <- .hx_div(.hx_mul(.sf_rc_calc(.hx_one(), .hx_add(.hx_one(), em), r), pow4), dm)
    S <- .hx_add(S, T)
    pow4 <- .hx_mul(pow4, g)
    Am <- Am1
  }
  t <- .hx_div(pow4, Am)
  X <- .hx_mul(.hx_sub(A0, x), t); Y <- .hx_mul(.hx_sub(A0, y), t); Z <- .hx_mul(.hx_sub(A0, z), t)
  P <- .hx_div(.hx_sub(.hx_sub(.hx_neg(X), Y), Z), 2)
  XY <- .hx_mul(X, Y); XYZ <- .hx_mul(XY, Z)
  E2 <- .hx_sub(.hx_add(.hx_add(XY, .hx_mul(X, Z)), .hx_mul(Y, Z)), .hx_mul(3, .hx_pow(P, 2)))
  E3 <- .hx_add(.hx_add(XYZ, .hx_mul(.hx_mul(2, E2), P)), .hx_mul(4, .hx_pow(P, 3)))
  E4 <- .hx_mul(.hx_add(.hx_add(.hx_mul(2, XYZ), .hx_mul(E2, P)), .hx_mul(3, .hx_pow(P, 3))), P)
  E5 <- .hx_mul(XYZ, .hx_pow(P, 2))
  Pn <- .hx_sub(24024, .hx_mul(5148, E2))
  Pn <- .hx_add(Pn, .hx_mul(2457, .hx_pow(E2, 2)))
  Pn <- .hx_add(Pn, .hx_mul(4004, E3))
  Pn <- .hx_sub(Pn, .hx_mul(.hx_mul(4158, E2), E3))
  Pn <- .hx_sub(Pn, .hx_mul(3276, E4))
  Pn <- .hx_add(Pn, .hx_mul(2772, E5))
  v1 <- .hx_div(.hx_mul(.hx_mul(pow4, .hx_power(Am, -1.5)), Pn), 24024)
  v2 <- .hx_mul(6, S)
  .hx_add(.hx_add(initial_integral, v1), v2)
}
.sf_carlson <- function(f) {
  prec <- .MPX$prec
  on.exit(.MPX$prec <- prec)
  .MPX$prec <- prec + 20
  tol <- .hx_mul(.hx_eps(), 2^10)
  v <- f(tol)
  .MPX$prec <- prec
  .mv_pos(v)
}
.sf_elliprf <- function(x, y, z) { x <- .mv_convert(x); y <- .mv_convert(y); z <- .mv_convert(z); .sf_carlson(function(r) .sf_rf_calc(x, y, z, r)) }
.sf_elliprj <- function(x, y, z, p) {
  x <- .mv_convert(x); y <- .mv_convert(y); z <- .mv_convert(z); p <- .mv_convert(p)
  .sf_carlson(function(r) .sf_rj_calc(x, y, z, p, r))
}
.sf_elliprd <- function(x, y, z) .sf_elliprj(x, y, z, z)

.hx_nint <- function(x) {
  x <- .mv_convert(x)
  .mv_mpf(.mp_round(.mp_floor(.mp_add(x$re, .mpf(1, -1), Inf)), 0, Inf))
}
.sf_ellip_reduce <- function(z, m, complete_fn) {
  x <- .hx_re(z)
  .MPX$prec <- .MPX$prec + max(0, .hx_mag(x))
  pi <- .mv_pos(.hx_pi())
  away <- .mv_cmp(.hx_abs(x), .hx_div(pi, 2)) > 0
  list(x = x, pi = pi, away = away)
}
.sf_ellipf <- function(phi, m) .hx_wrapped(function(z, m) {
  if (.hx_iszero(m)) return(.hx_add(z, m))
  if (.hx_iszero(z)) return(.hx_mul(z, m))
  rd <- .sf_ellip_reduce(z, m)
  if (.hx_numeq(m, 1) && rd$away) return(.mv_inf())
  if (rd$away) {
    d <- .hx_nint(.hx_div(rd$x, rd$pi))
    z <- .hx_sub(z, .hx_mul(rd$pi, d))
    P <- .hx_mul(.hx_mul(2, d), .sf_ellipk(m))
  } else P <- 0
  cs <- .hx_cos_sin(z)
  .hx_add(.hx_mul(cs$sin, .sf_elliprf(.hx_pow(cs$cos, 2), .hx_sub(1, .hx_mul(m, .hx_pow(cs$sin, 2))), 1)), P)
}, phi, m)
.sf_ellipe <- function(...) {
  args <- list(...)
  if (length(args) == 1L) return(.hx_wrapped(function(m) .sf_ellipe1(m), args[[1L]]))
  .hx_wrapped(function(z, m) {
    if (.hx_iszero(m)) return(.hx_add(z, m))
    if (.hx_iszero(z)) return(.hx_mul(z, m))
    rd <- .sf_ellip_reduce(z, m)
    if (rd$away) {
      d <- .hx_nint(.hx_div(rd$x, rd$pi))
      z <- .hx_sub(z, .hx_mul(rd$pi, d))
      P <- .hx_mul(.hx_mul(2, d), .sf_ellipe(m))
    } else P <- 0
    terms <- function() {
      cs <- .hx_cos_sin(z)
      x <- .hx_pow(cs$cos, 2); y <- .hx_sub(1, .hx_mul(m, .hx_pow(cs$sin, 2)))
      RF <- .sf_elliprf(x, y, 1); RD <- .sf_elliprd(x, y, 1)
      list(.hx_mul(cs$sin, RF), .hx_div(.hx_mul(.hx_mul(.hx_neg(m), .hx_pow(cs$sin, 3)), RD), 3))
    }
    .hx_add(.hx_sum_accurately(terms), P)
  }, args[[1L]], args[[2L]])
}
.sf_ellippi <- function(...) {
  args <- list(...)
  complete <- length(args) == 2L
  f <- function(n, z, m) {
    if (complete) {
      if (.hx_numeq(m, 1)) return(.hx_div(.mv_ninf(), .hx_sign(.hx_sub(n, 1))))
      away <- FALSE
    } else {
      rd <- .sf_ellip_reduce(z, m); away <- rd$away
    }
    if (away) {
      d <- .hx_nint(.hx_div(rd$x, rd$pi))
      z <- .hx_sub(z, .hx_mul(rd$pi, d))
      P <- .hx_mul(.hx_mul(2, d), .sf_ellippi(n, m))
    } else P <- 0
    terms <- function() {
      if (complete) { c <- .hx_zero(); s <- .hx_one() }
      else { cs <- .hx_cos_sin(z); c <- cs$cos; s <- cs$sin }
      x <- .hx_pow(c, 2); y <- .hx_sub(1, .hx_mul(m, .hx_pow(s, 2)))
      RF <- .sf_elliprf(x, y, 1)
      RJ <- .sf_elliprj(x, y, 1, .hx_sub(1, .hx_mul(n, .hx_pow(s, 2))))
      list(.hx_mul(s, RF), .hx_div(.hx_mul(.hx_mul(n, .hx_pow(s, 3)), RJ), 3))
    }
    .hx_add(.hx_sum_accurately(terms), P)
  }
  if (complete) return(.hx_wrapped(function(n, m) f(n, .hx_div(.hx_pi(), 2), m), args[[1L]], args[[2L]]))
  .hx_wrapped(f, args[[1L]], args[[2L]], args[[3L]])
}

# ---------------------------------------------------------------------------
# Lambert W (mpmath's lambertw): a starting value from the series at the
# branch point, the asymptotic expansion or a double precision fit, then
# Halley iteration at prec + 20 + mag(k or 1) bits

.sf_lambertw_special <- function(z, k) {
  if (.hx_iszero(z)) return(if (k == 0) z else .hx_add(.mv_ninf(), z))
  if (identical(z$sp, "inf")) return(if (k == 0) z else .hx_add(z, .hx_mul(.hx_mul(2 * k, .hx_pi()), .hx_j())))
  if (identical(z$sp, "-inf")) return(.hx_add(.hx_neg(z), .hx_mul(.hx_mul(2 * k + 1, .hx_pi()), .hx_j())))
  .hx_log(z)
}

# _lambertw_approx_hybrid: Python float or complex arithmetic
.sf_lambertw_approx <- function(z, k) {
  out <- function(v, cx) if (cx) .mv_mpc(.mp_from_dbl(Re(v)), .mp_from_dbl(Im(v))) else .mv_float(Re(v))
  imag_sign <- 0
  x <- .mv_dbl(.hx_re(z)); yv <- .hx_im(z)
  if (!.hx_iszero(yv)) imag_sign <- if (.mv_cmp(yv, 0) < 0) -1 else 1
  y <- .mv_dbl(yv)
  zc <- complex(real = x, imaginary = y); cx <- TRUE
  r <- -0.367879441171442
  if (k == 0) {
    if (-4.0 < y && y < 4.0 && -1.0 < x && x < 2.5) {
      if (imag_sign != 0) {
        if (y > 1.00) return(out(complex(real = 0.876, imaginary = 0.645) + complex(real = 0.118, imaginary = -0.174) * (zc - complex(real = 0.75, imaginary = 2.5)), TRUE))
        if (y > 0.25) return(out(complex(real = 0.505, imaginary = 0.204) + complex(real = 0.375, imaginary = -0.132) * (zc - complex(real = 0.75, imaginary = 0.5)), TRUE))
        if (y < -1.00) return(out(complex(real = 0.876, imaginary = -0.645) + complex(real = 0.118, imaginary = 0.174) * (zc - complex(real = 0.75, imaginary = -2.5)), TRUE))
        if (y < -0.25) return(out(complex(real = 0.505, imaginary = -0.204) + complex(real = 0.375, imaginary = 0.132) * (zc - complex(real = 0.75, imaginary = -0.5)), TRUE))
      }
      if (x < -0.5) {
        if (imag_sign >= 0) return(out(complex(real = -0.318, imaginary = 1.34) + complex(real = -0.697, imaginary = -0.593) * (zc + 1), TRUE))
        return(out(complex(real = -0.318, imaginary = -1.34) + complex(real = -0.697, imaginary = 0.593) * (zc + 1), TRUE))
      }
      if (imag_sign == 0 && x > r) { zc <- x; cx <- FALSE }
      if (x < -0.2) return(out(-1 + 2.33164398159712 * (zc - r)^0.5 - 1.81218788563936 * (zc - r), cx))
      if (x < 0.5) return(out(zc, cx))
      return(out(0.2 + 0.3 * zc, cx))
    }
    if (imag_sign == 0 && x > 0.0) { L1 <- log(x); L2 <- log(L1); cx <- FALSE }
    else { L1 <- log(zc); L2 <- log(L1) }
  } else {
    if (imag_sign == 0 && r < x && x < 0.0) { zc <- x; cx <- FALSE }
    if (imag_sign >= 0 && y < 0.1 && -0.6 < x && x < -0.2)
      return(out(-1 - 2.33164398159712 * (zc - r)^0.5 - 1.81218788563936 * (zc - r), cx))
    if (imag_sign == 0 && -0.2 <= x && x < 0.0) { L1 <- log(-x); return(out(L1 - log(-L1), FALSE)) }
    cx <- TRUE
    L1 <- if (imag_sign == -1 && y == 0 && x < 0.0) log(zc) - complex(imaginary = 3.1415926535897932)
          else log(zc) - complex(imaginary = 6.2831853071795865)
    L2 <- log(L1)
  }
  out(L1 - L2 + L2 / L1 + L2 * (L2 - 2) / (2 * L1^2), cx)
}

# _lambertw_series; it may raise the context precision for the caller
.sf_lambertw_series <- function(z, k, tol) {
  magz <- .hx_mag(z)
  if (-10 < magz && magz < 900 && -1000 < k && k < 1000) {
    if (magz < 1 && .mv_cmp(.hx_abs(.hx_add(z, 0.36787944117144)), 0.05) < 0) {
      imz <- .hx_im(z)
      if (k == 0 || (k == -1 && .mv_cmp(imz, 0) >= 0) || (k == 1 && .mv_cmp(imz, 0) < 0)) {
        delta <- .hx_sum_accurately(function() list(z, .hx_exp(.mv_int(-1))))
        cancellation <- -.hx_mag(delta)
        .MPX$prec <- .MPX$prec + cancellation
        p <- .hx_sqrt(.hx_mul(2, .hx_add(.hx_mul(.hx_exp(.hx_one()), z), 1)))
        .MPX$prec <- .MPX$prec - cancellation
        u <- list(.mv_mpf(.mpf(-1, 0)), .mv_mpf(.mpf(1, 0)))   # u[[l + 1]] is u_l
        a <- list(.mv_mpf(.mpf(1, 1)), .mv_mpf(.mpf(-1, 0)))
        if (k != 0) p <- .hx_neg(p)
        s <- .hx_zero()
        for (l in 0:(max(2, cancellation) - 1)) {
          if (l + 1 > length(u)) {
            a[[l + 1]] <- .hx_fsum(lapply(seq_len(max(0, l - 2)) + 1, function(j) .hx_mul(u[[j + 1]], u[[l + 2 - j]])))
            u[[l + 1]] <- .hx_sub(.hx_sub(.hx_div(.hx_mul(l - 1, .hx_add(.hx_div(u[[l - 1]], 2), .hx_div(a[[l - 1]], 4))), l + 1),
                                          .hx_div(a[[l + 1]], 2)), .hx_div(u[[l]], l + 1))
          }
          term <- .hx_mul(u[[l + 1]], .hx_pow(p, l))
          s <- .hx_add(s, term)
          if (.hx_mag(term) < -tol) return(list(s, TRUE))
        }
        .MPX$prec <- .MPX$prec + cancellation %/% 2
        return(list(s, FALSE))
      }
    }
    if (k == 0 || k == -1) return(list(.sf_lambertw_approx(z, k), FALSE))
  }
  if (k == 0) {
    if (magz < -1) return(list(.hx_mul(z, .hx_sub(1, z)), FALSE))
    L1 <- .hx_log(z); L2 <- .hx_log(L1)
  } else if (k == -1 && .hx_iszero(.hx_im(z)) && .mv_cmp(.hx_re(z), -0.36787944117144) > 0 && .mv_cmp(.hx_re(z), 0) < 0) {
    L1 <- .hx_log(.hx_neg(z))
    return(list(.hx_sub(L1, .hx_log(.hx_neg(L1))), FALSE))
  } else {
    L1 <- .hx_add(.hx_log(z), .hx_mul(.hx_mul(.mv_mpc(.mp_zero, .mpf(1, 1)), .hx_pi()), k))
    L2 <- .hx_log(L1)
  }
  list(.hx_add(.hx_add(.hx_sub(L1, L2), .hx_div(L2, L1)), .hx_div(.hx_mul(L2, .hx_sub(L2, 2)), .hx_mul(2, .hx_pow(L1, 2)))), FALSE)
}

.sf_lambertw <- function(z, k = 0) {
  z <- .mv_convert(z)
  k <- as.integer(.mv_dbl(k))
  if (.hx_iszero(z) || .mv_sp(z)) return(.sf_lambertw_special(z, k))
  prec <- .MPX$prec
  on.exit(.MPX$prec <- prec)
  .MPX$prec <- .MPX$prec + 20 + .hx_mag(if (k != 0) k else 1)
  wp <- .MPX$prec
  tol <- wp - 5
  sr <- .sf_lambertw_series(z, k, tol)
  w <- sr[[1L]]
  if (!sr[[2L]]) {
    two <- .mv_mpf(.mpf(1, 1))
    for (i in 0:99) {
      ew <- .hx_exp(w)
      wew <- .hx_mul(w, ew)
      wewz <- .hx_sub(wew, z)
      wn <- .hx_sub(w, .hx_div(wewz, .hx_sub(.hx_add(wew, ew),
                                             .hx_div(.hx_mul(.hx_add(w, two), wewz), .hx_add(.hx_mul(two, w), two)))))
      done <- .hx_mag(.hx_sub(wn, w)) <= .hx_mag(wn) - tol
      w <- wn
      if (done) break
    }
  }
  .MPX$prec <- prec
  .mv_pos(w)
}
