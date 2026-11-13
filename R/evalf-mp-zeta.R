# mpmath's zeta functions as SymPy's evalf calls them: zeta(s) and the
# Hurwitz zeta(s, a) (mpmath/functions/zeta.py over libmp/gammazeta.py),
# altzeta, bernpoly, harmonic and polylog, and SymPy's dirichlet_eta, whose
# _eval_evalf evaluates its rewrite in terms of zeta in Float arithmetic.
# The Riemann-Siegel formula (used by mpmath only for |Im s| > 500 * prec) is
# not ported: such s go to the Euler-Maclaurin sum.

# ---------------------------------------------------------------------------
# libmp level (mpf = list(m, e), mpc = list(re, im))

# borwein_coefficients(n), cached
.zt_borwein <- function(n) {
  key <- paste0("borwein", n)
  if (!is.null(v <- .mp_cache[[key]])) return(v)
  ds <- vector("list", n + 1L)
  d <- gmp::as.bigz(1); s <- gmp::as.bigz(1); ds[[1L]] <- s
  for (i in seq_len(n)) {
    d <- d * 4 * gmp::as.bigz(n + i - 1) * gmp::as.bigz(n - i + 1)
    d <- gmp::divq.bigz(d, gmp::as.bigz(2 * i) * gmp::as.bigz(2 * i - 1))
    s <- s + d
    ds[[i + 1L]] <- s
  }
  assign(key, ds, envir = .mp_cache)
  ds
}

# log_int_fixed, exp_fixed and cos_sin_fixed: fixed point values with wp
# fractional bits, computed through the mpf functions with ten guard bits;
# log_int_fixed caches its values, as mpmath does
.zt_log_int <- function(k, wp) {
  key <- paste0("logint", k, "_", wp)
  if (!is.null(v <- .mp_cache[[key]])) return(v)
  v <- .to_fixed(.mp_log(.mp_round(k, 0, Inf), wp + 10), wp)
  assign(key, v, envir = .mp_cache)
  v
}
.zt_exp_fixed <- function(u, wp) .to_fixed(.mp_exp(.mp_round(u, -wp, Inf), wp + 10), wp)
.zt_cos_sin_fixed <- function(t, wp) {
  cs <- .mp_sincos(.mp_round(t, -wp, Inf), wp + 10)
  list(.to_fixed(cs$cos, wp), .to_fixed(cs$sin, wp))
}

# list_primes(n): the primes up to n
.zt_primes <- function(n) {
  if (n < 2) return(integer(0))
  sv <- rep(TRUE, n); sv[1L] <- FALSE
  i <- 2L
  while (i * i <= n) {
    if (sv[i]) sv[seq(i * i, n, by = i)] <- FALSE
    i <- i + 1L
  }
  which(sv)
}

# mpf_pow(s, t, prec) with t normalized first
.zt_mpf_pow <- function(s, t, p) {
  t <- .mp_norm(t)
  if (t$m == 0) return(.mpf(1, 0))
  .mp_pow_real(s, t, p)
}
# mpc_pow(z, w, prec)
.zt_mpc_pow <- function(z, w, p) {
  if (w$im$m == 0) return(.mpc_pow_mpf(z, w$re, p))
  .mpc_exp(.mpc_mul(.mpc_log(z, p + 10), w, p + 10), p)
}

# mpf_zeta_int: zeta(s) for an integer s, at prec + 20 bits on some paths
.zt_zeta_int <- function(s, prec) {
  s <- as.numeric(s)
  wp <- prec + 20
  if (s < 2) {
    if (s == 1) .mp_error("ValueError", "zeta(1) pole")
    if (s == 0) return(.mpf(-1, -1))
    return(.mp_div(.mp_from_q(.mp_bernoulli(-s + 1), wp), .mp_round(s - 1, 0, Inf), prec))
  }
  if (s >= wp) return(.mp_round(1, 0, prec))
  if (s >= wp * 0.431) {
    one <- .pow2(wp)
    t <- one + .pow2(wp - s) + gmp::divq.bigz(one, gmp::as.bigz(3)^s) + .pow2(max(0, wp - s * 2))
    return(.mp_round(t, -wp, prec))
  }
  m <- wp / (s - 1) + 1
  if (m < 30) {
    needed <- as.integer(2^m + 1)
    if (needed < as.integer(wp / 2.54 + 5) / 10) {
      t <- .mpf(1, 0)
      for (k in .zt_primes(needed)) {
        powprec <- as.integer(wp - s * log(k) / log(2))
        if (powprec < 2) break
        a <- .mp_sub(.mpf(1, 0), .mp_pow_int(.mp_round(k, 0, Inf), -s, powprec), wp)
        t <- .mp_mul(t, a, wp)
      }
      return(.mp_div(.mpf(1, 0), t, wp))
    }
  }
  n <- as.integer(wp / 2.54 + 5)
  d <- .zt_borwein(n); dn <- d[[n + 1L]]
  t <- gmp::as.bigz(0); one <- .pow2(wp)
  for (k in 0:(n - 1L)) {
    c <- d[[k + 1L]] - dn
    if (k %% 2L == 1L) c <- -c
    t <- t + gmp::divq.bigz(c * one, gmp::as.bigz(k + 1)^s)
  }
  t <- gmp::divq.bigz(t * one, -dn)
  t <- gmp::divq.bigz(t * one, one - .pow2(wp + 1 - s))
  .mp_round(t, -wp - wp, prec)
}

# mpf_zeta(s, prec, alt): zeta(s), or the alternating zeta (Dirichlet eta)
.zt_mpf_zeta <- function(s, prec, alt = FALSE) {
  one <- .mpf(1, 0)
  if (s$m == 0) return(if (alt) .mpf(1, -1) else .mpf(-1, -1))
  s <- .mp_norm(s)
  wp <- prec + 20
  sgn <- s$m < 0; ex <- s$e; bc <- .bitlen(abs(s$m))
  # first term vanishes?
  if (!sgn && ex + bc > log(wp) / log(2) + 2) return(.mp_round(1, 0, prec))
  if (ex >= 0) {
    n <- s$m * .pow2(ex)
    if (alt) {
      if (n == 1) return(.mp_ln2(prec))
      z <- .zt_zeta_int(n, wp)
      q <- .mp_sub(one, .zt_mpf_pow(.mpf(1, 1), .mp_sub(one, s, wp), wp), wp)
      return(.mp_mul(z, q, prec))
    }
    return(.zt_zeta_int(n, prec))
  }
  # negative s: the reflection formula
  if (sgn) {
    if (alt) {
      q <- .mp_sub(one, .zt_mpf_pow(.mpf(1, 1), .mp_sub(one, s, wp), wp), wp)
      return(.mp_mul(.zt_mpf_zeta(s, wp), q, prec))
    }
    y <- .mp_sub(one, s, 10 * wp)
    a <- .mp_gamma(y, wp)
    b <- .zt_mpf_zeta(y, wp)
    c <- .mp_cos_sin(.mp_ldexp(s, -1), wp, pi = TRUE)$sin
    wp2 <- wp + max(0, ex + bc)
    pi <- .mp_pi(wp + wp2)
    d <- .mp_div(.zt_mpf_pow(.mp_ldexp(pi, 1), s, wp2), pi, wp2)
    return(.mp_mul(a, .mp_mul(b, .mp_mul(c, d, wp), wp), prec))
  }
  # near the pole
  r <- .mp_sub(one, s, wp)
  pole_dist <- -2 * .mp_mag(r)
  if (pole_dist > wp) {
    if (alt) return(.mp_ln2(prec))
    q <- .mp_neg(.mp_div(one, r, wp))
    return(.mp_add(q, .mp_euler(wp), prec))
  }
  wp <- wp + max(0, pole_dist)
  # Borwein's algorithm
  n <- as.integer(wp / 2.54 + 5)
  d <- .zt_borwein(n); dn <- d[[n + 1L]]
  t <- gmp::as.bigz(0)
  sf <- .to_fixed(s, wp)
  for (k in 0:(n - 1L)) {
    u <- .fx_shift(-sf * .zt_log_int(k + 1, wp), -wp)
    w <- (d[[k + 1L]] - dn) * .zt_exp_fixed(u, wp)
    t <- if (k %% 2L == 1L) t - w else t + w
  }
  t <- gmp::divq.bigz(t, -dn)
  t <- .mp_round(t, -wp, wp)
  if (alt) return(.mp_round(t$m, t$e, prec))
  q <- .mp_sub(one, .zt_mpf_pow(.mpf(1, 1), .mp_sub(one, s, wp), wp), wp)
  .mp_div(t, q, prec)
}

# mpc_zeta(s, prec, alt, force); NotImplementedError for |s| > prec
.zt_mpc_zeta <- function(s, prec, alt = FALSE, force = FALSE) {
  re <- s$re; im <- s$im
  if (im$m == 0) return(.mpc(.zt_mpf_zeta(re, prec, alt)))
  if (!force && .mp_cmp(.mpc_abs(s, 10), .mp_round(prec, 0, Inf)) > 0) .mp_error("NotImplementedError", "")
  wp <- prec + 20
  one <- .mpf(1, 0); mone <- .mpc(one)
  # near the pole
  r <- .mpc_sub(mone, s, wp)
  pole_dist <- -2 * .mp_mag(.mpc_abs(r, 10))
  if (pole_dist > wp) {
    if (alt) {
      q <- .mp_ln2(wp)
      y <- .mp_mul(q, .mp_euler(wp), wp)
      g <- .mp_sub(y, .mp_ldexp(.mp_mul(q, q, wp), -1), Inf)
      g <- .mp_neg(g)
      z <- .mpc(.mp_add(.mp_mul(r$re, g, wp), q, wp), .mp_mul(r$im, g, wp))
      return(.mpc(.mp_round(z$re$m, z$re$e, prec), .mp_round(z$im$m, z$im$e, prec)))
    }
    q <- .mpc_neg(.mpc_div(mone, r, wp))
    return(.mpc(.mp_add(q$re, .mp_euler(wp), prec), .mp_round(q$im$m, q$im$e, prec)))
  }
  wp <- wp + max(0, pole_dist)
  # the reflection formula
  if (re$m < 0) {
    two <- .mpc(.mpf(1, 1))
    if (alt) {
      q <- .mpc_sub(mone, .zt_mpc_pow(two, .mpc_sub(mone, s, wp), wp), wp)
      return(.mpc_mul(.zt_mpc_zeta(s, wp), q, prec))
    }
    y <- .mpc_sub(mone, s, 10 * wp)
    a <- .mpc_gamma(y, wp)
    b <- .zt_mpc_zeta(y, wp)
    c <- .mpc_cos_sin_pi(.mpc(.mp_ldexp(re, -1), .mp_ldexp(im, -1)), wp)$sin
    mag <- max(.mp_mag(re), .mp_mag(im))
    wp2 <- wp + max(0, mag)
    pi <- .mp_pi(wp + wp2)
    d <- .zt_mpc_pow(.mpc(.mp_ldexp(pi, 1)), s, wp2)
    d <- .mpc(.mp_div(d$re, pi, wp2), .mp_div(d$im, pi, wp2))
    return(.mpc_mul(a, .mpc_mul(b, .mpc_mul(c, d, wp), wp), prec))
  }
  # Borwein's algorithm
  n <- as.integer(wp / 2.54 + 5)
  n <- n + as.integer(0.9 * abs(trunc(.mp_to_dbl(im))))
  d <- .zt_borwein(n); dn <- d[[n + 1L]]
  ref <- .to_fixed(re, wp); imf <- .to_fixed(im, wp)
  tre <- gmp::as.bigz(0); tim <- gmp::as.bigz(0)
  one_2wp <- .pow2(2 * wp)
  critical_line <- .mp_cmp(re, .mpf(1, -1)) == 0
  for (k in 0:(n - 1L)) {
    lg <- .zt_log_int(k + 1, wp)
    w <- if (critical_line) gmp::divq.bigz(one_2wp, .mp_isqrt(gmp::as.bigz(k + 1) * one_2wp))
         else .zt_exp_fixed(.fx_shift(-ref * lg, -wp), wp)
    w <- if (k %% 2L == 1L) w * (dn - d[[k + 1L]]) else w * (d[[k + 1L]] - dn)
    cs <- .zt_cos_sin_fixed(.fx_shift(-imf * lg, -wp), wp)
    tre <- tre + .fx_shift(w * cs[[1L]], -wp)
    tim <- tim + .fx_shift(w * cs[[2L]], -wp)
  }
  tre <- gmp::divq.bigz(tre, -dn); tim <- gmp::divq.bigz(tim, -dn)
  tt <- .mpc(.mp_round(tre, -wp, wp), .mp_round(tim, -wp, wp))
  if (alt) return(.mpc(.mp_round(tt$re$m, tt$re$e, prec), .mp_round(tt$im$m, tt$im$e, prec)))
  q <- .mpc_sub(mone, .zt_mpc_pow(.mpc(.mpf(1, 1)), r, wp), wp)
  .mpc_div(tt, q, prec)
}

# mpf_psi0: the digamma function of a real argument
.zt_psi0 <- function(x, prec) {
  wp <- prec + 10
  x <- .mp_norm(x)
  sgn <- x$m < 0; ex <- x$e; bc <- .bitlen(abs(x$m))
  if (x$m == 0 || (ex >= 0 && sgn)) .mp_error("ValueError", "polygamma pole")
  one_f <- .mpf(1, 0)
  # near 0, fixed point arithmetic becomes bad
  if (ex + bc < -5) {
    v <- .zt_psi0(.mp_add(x, one_f, prec), prec)
    return(.mp_sub(v, .mp_div(one_f, x, wp), prec))
  }
  # reflection formula
  if (sgn && ex + bc > 3) {
    cs <- .mp_cos_sin(x, wp, pi = TRUE)
    q <- .mp_mul(.mp_div(cs$cos, cs$sin, wp), .mp_pi(wp), wp)
    p <- .zt_psi0(.mp_sub(one_f, x, wp), wp)
    return(.mp_sub(p, q, prec))
  }
  # the logarithmic term is accurate enough
  if (!sgn && bc + ex > wp) return(.mp_log(.mp_sub(x, one_f, wp), prec))
  # initial recurrence to obtain a large enough x
  m <- as.numeric(if (ex >= 0) x$m * .pow2(ex) else sign(x$m) * gmp::divq.bigz(abs(x$m), .pow2(-ex)))
  n <- as.integer(0.11 * wp) + 2
  s <- gmp::as.bigz(0)
  xf <- .to_fixed(x, wp)
  one <- .pow2(wp)
  if (m < n) for (k in seq(m, n - 1)) {
    s <- s - gmp::divq.bigz(one * one, xf)
    xf <- xf + one
  }
  xf <- xf - one
  # logarithmic term
  s <- s + .to_fixed(.mp_log(.mp_round(xf, -wp, wp), wp), wp)
  # endpoint term in the Euler-Maclaurin expansion
  s <- s + gmp::divq.bigz(one * one, 2 * xf)
  # Euler-Maclaurin remainder sum
  x2 <- .fx_shift(xf * xf, -wp)
  t <- one; prev <- gmp::as.bigz(0); k <- 1
  repeat {
    t <- .fx_shift(t * x2, -wp)
    b <- .mp_from_q(.mp_bernoulli(2 * k), wp)
    offset <- b$e + 2 * wp
    bman <- abs(b$m)
    term <- gmp::divq.bigz(.fx_shift(bman, offset), t * (2 * k))
    s <- if (k %% 2 == 1) s - term else s + term
    if (k > 2 && term >= prev) break
    prev <- term
    k <- k + 1
  }
  .mp_round(s, -wp, wp)
}

# ---------------------------------------------------------------------------
# context level (mpv values at .MPX$prec)

# ctx.bernoulli (mpmath's B_1 = -1/2)
.zt_bernoulli <- function(n) {
  if (n == 1) return(.mv_mpf(.mpf(-1, -1)))
  .mv_mpf(.mp_from_q(.mp_bernoulli(n), .MPX$prec))
}

# ctx.harmonic for a real argument
.zt_harmonic <- function(x) {
  x <- .mv_convert(x)
  if (.hx_iszero(x) || .mv_sp(x)) return(x)
  prec <- .MPX$prec
  a <- .zt_psi0(.mp_add(.mpf(1, 0), x$re, prec + 5), prec)
  .mv_mpf(.mp_add(a, .mp_euler(prec + 5), prec))
}

# ctx._zeta and ctx._altzeta (_wrap_libmp_function over mpf_zeta, mpc_zeta)
.zt_lib_zeta <- function(x, alt = FALSE) {
  x <- .mv_convert(x)
  prec <- .MPX$prec
  if (x$k == "mpf") {
    if (.mv_sp(x)) return(if (x$sp == "inf") .hx_one() else .mv_nan())
    return(.mv_mpf(.zt_mpf_zeta(x$re, prec, alt)))
  }
  if (.mv_sp(x)) return(.mv_nan(TRUE))
  r <- .zt_mpc_zeta(.mpc(x$re, x$im), prec, alt)
  .mv_mpc(r$re, r$im)
}
.zt_not_implemented <- function(expr) tryCatch(expr, NotImplementedError = function(e) NULL)

# ctx.powm1(x, y): x^y - 1
.zt_powm1 <- function(x, y) {
  w <- .hx_sub(.hx_pow(x, y), .hx_one())
  M <- .hx_mag(w)
  if (M > -8) return(w)
  if (.hx_iszero(w) && (.hx_iszero(y) || (any(vapply(list(1, -1), function(v) isTRUE(.mv_eq(x, v)), TRUE)) && .hx_isint(y))))
    return(w)
  lnx <- .hx_log(x)
  if (.hx_mag(y) + .hx_mag(lnx) < -.MPX$prec) {
    ly <- .hx_mul(lnx, y)
    return(.hx_add(ly, .hx_div(.hx_pow(ly, 2), 2)))
  }
  .hx_sum_accurately(function() list(.hx_pow(x, y), .mv_int(-1)), 1)
}

# ctx.altzeta
.zt_altzeta <- function(s) {
  r <- .zt_not_implemented(.zt_lib_zeta(s, alt = TRUE))
  if (!is.null(r)) return(r)
  prec <- .MPX$prec
  on.exit(.MPX$prec <- prec)
  .MPX$prec <- prec + 10
  s <- .mv_convert(s)
  r <- if (isTRUE(.mv_eq(s, 1))) .hx_add(.mv_mpf(.mp_ln2(.MPX$prec)), .hx_mul(0, s))
       else .hx_mul(.hx_neg(.zt_powm1(.mv_int(2), .hx_sub(1, s))), .zt_zeta(s))
  .MPX$prec <- prec
  .mv_pos(r)
}

# ctx.bernpoly(n, z)
.zt_bernpoly <- function(n, z) {
  n <- as.integer(n)
  if (n < 0) .mp_error("ValueError", "Bernoulli polynomials only defined for n >= 0")
  z <- .mv(z)
  if (.hx_isinf(z)) return(.hx_pow(z, n))
  if (.hx_isnan(z)) return(z)
  if (n <= 3L) {
    if (n == 0L) return(.hx_pow(z, .mv_int(0)))
    if (n == 1L) return(.hx_sub(z, 0.5))
    if (n == 2L) return(.hx_div(.hx_add(.hx_mul(.hx_mul(6, z), .hx_sub(z, 1)), 1), 6))
    return(.hx_mul(z, .hx_add(.hx_mul(z, .hx_sub(z, 1.5)), 0.5)))
  }
  if (isTRUE(.mv_eq(z, 0)) || isTRUE(.mv_eq(z, 1))) return(.zt_bernoulli(n))
  if (isTRUE(.mv_eq(z, 0.5))) return(.hx_mul(.hx_sub(.hx_ldexp(.mv_int(1), 1L - n), 1), .zt_bernoulli(n)))
  if (.mv_cmp(.hx_abs(z), 2) > 0) {
    terms <- function() {
      out <- list(); t <- .hx_one(); out[[1L]] <- t
      r <- .hx_div(.hx_one(), z)
      t <- .hx_mul(.hx_mul(t, n), r)
      out[[2L]] <- .hx_div(.hx_neg(t), 2)
      k <- 2L
      while (k <= n) {
        t <- .hx_mul(.hx_div(.hx_mul(t, n + 1L - k), k), r)
        if (k %% 2L == 0L) out[[length(out) + 1L]] <- .hx_mul(t, .zt_bernoulli(k))
        k <- k + 1L
      }
      out
    }
    return(.hx_mul(.hx_sum_accurately(terms), .hx_pow(z, n)))
  }
  terms <- function() {
    out <- list(.zt_bernoulli(n)); t <- .hx_one(); k <- 1L
    while (k < n - 1L) {
      t <- .hx_mul(.hx_div(.hx_mul(t, n + 1L - k), k), z)
      m <- n - k
      if (m %% 2L == 0L) out[[length(out) + 1L]] <- .hx_mul(t, .zt_bernoulli(m))
      k <- k + 1L
    }
    t <- .hx_mul(.hx_div(.hx_mul(t, 2), n - 1L), z)
    out[[length(out) + 1L]] <- .hx_div(.hx_neg(t), 2)
    t <- .hx_mul(.hx_div(t, n), z)
    out[[length(out) + 1L]] <- t
    out
  }
  .hx_sum_accurately(terms)
}

# ctx.zeta(s, a, derivative)
.zt_zeta <- function(s, a = .mv_int(1), d = 0) {
  a1 <- isTRUE(.mv_eq(a, 1))
  if (a1 && d == 0) {
    r <- .zt_not_implemented(.zt_lib_zeta(s))
    if (!is.null(r)) return(r)
  }
  s <- .mv_convert(s)
  prec <- .MPX$prec
  if (.hx_iszero(s) && d == 0) return(.hx_sub(.mv_mpf(.mpf(1, -1)), .hx_param(a)$v))
  if (isTRUE(.mv_eq(s, 1))) return(.mv_inf())
  abss <- .hx_abs(s)
  if (.hx_isinf(abss)) {
    if (isTRUE(.mv_eq(.hx_re(s), .mv_inf()))) return(if (d == 0) .hx_one() else .hx_zero())
    return(.hx_mul(s, 0))
  }
  if (.hx_isnan(abss)) return(.hx_div(1, s))
  if (a1 && d == 0 && .mv_cmp(.hx_re(s), 2 * prec) > 0) return(.hx_add(.hx_one(), .hx_pow(.mv_mpf(.mpf(1, 1)), .hx_neg(s))))
  .mv_pos(.zt_hurwitz(s, a, d))
}

# ctx._hurwitz
.zt_hurwitz <- function(s, a = .mv_int(1), d = 0) {
  prec <- .MPX$prec
  on.exit(.MPX$prec <- prec)
  extraprec <- 10
  .MPX$prec <- .MPX$prec + extraprec
  pa <- .hx_param(a); a <- pa$v; atype <- pa$t
  if (.mv_cmp(.hx_re(a), 0) < 0 && .hx_isnpint(a)) .mp_error("ValueError", "Hurwitz zeta complex infinity")
  if (.mv_cmp(.hx_re(s), 0) < 0) {
    r <- .zt_not_implemented(.zt_hurwitz_reflection(s, a, d, atype))
    if (!is.null(r)) return(r)
  }
  repeat {
    .MPX$prec <- prec + extraprec
    tt <- .zt_hurwitz_em(s, a, d, prec + 10)
    T12 <- .hx_add(tt[[1L]], tt[[2L]])
    cancellation <- .hx_mag(tt[[1L]]) - .hx_mag(T12)
    if (!is.nan(cancellation) && cancellation < extraprec) return(T12)
    extraprec <- max(2 * extraprec, min(cancellation + 5, 100 * prec))
    if (extraprec > 100 * prec) .mp_error("NoConvergence", "zeta: too much cancellation")
  }
}

.zt_hurwitz_reflection <- function(s, a, d, atype) {
  if (d != 0) .mp_error("NotImplementedError", "")
  res <- .hx_re(s)
  negs <- .hx_neg(s)
  # integer reflection formula
  if (.hx_isnpint(s)) {
    n <- as.numeric(.hx_int_value(res))
    if (n <= 0) return(.hx_div(.zt_bernpoly(1 - n, a), n - 1))
  }
  if (!(atype %in% c("Q", "Z"))) .mp_error("NotImplementedError", "")
  t <- .hx_sub(1, s)
  # a standardized to 0 < a <= 1
  v <- .mv_int(0); shift <- 0
  b <- .mv_convert(a)
  while (.mv_cmp(.hx_re(b), 1) > 0) {
    b <- .hx_sub(b, 1); v <- .hx_sub(v, .hx_pow(b, negs)); shift <- shift - 1
  }
  while (.mv_cmp(.hx_re(b), 0) <= 0) {
    v <- .hx_add(v, .hx_pow(b, negs)); b <- .hx_add(b, 1); shift <- shift + 1
  }
  # rational reflection formula
  p <- gmp::numerator(a$q); q <- gmp::denominator(a$q)
  p <- p + shift * q
  q <- as.integer(as.numeric(q))
  g <- .hx_fsum(lapply(seq_len(q), function(k) {
    arg <- .hx_sub(.hx_div(t, 2), .hx_mul(2 * k, b))
    .hx_mul(.hx_cospi(arg), .zt_hurwitz(t, .mv_new("tuple", q = gmp::as.bigq(k, q))))
  }))
  g <- .hx_mul(g, .hx_div(.hx_mul(2, .hx_gamma(t)), .hx_pow(.hx_mul(.hx_mul(2, .hx_pi()), q), t)))
  .hx_add(v, g)
}

# _hurwitz_em for derivative 0: the truncated L-series and the
# Euler-Maclaurin tail
.zt_hurwitz_em <- function(s, a, d, prec) {
  # derivatives: the general recurrence (evalf-eval.R)
  if (d != 0) return(.ev_hurwitz_em_d(s, a, d, prec))
  a <- .mv_convert(a)
  tol <- -prec
  M1 <- 0; M2 <- prec %/% 3; N <- M2
  lsum <- .mv_int(0)
  if (.hx_isint(s)) s <- .mv_int(.hx_int_value(.hx_re(s)))
  s1 <- .hx_sub(s, 1)
  repeat {
    l <- .zt_zetasum(s, .hx_add(M1, a), M2 - M1 - 1)
    lsum <- .hx_add(lsum, l)
    M2a <- .hx_add(M2, a)
    logM2ad <- .hx_pow(.hx_log(M2a), .mv_int(0))
    rM2a <- .hx_div(1, M2a)
    M2as <- .hx_pow(M2a, .hx_neg(s))
    tailsum <- .hx_div(1, .hx_mul(s1, .hx_pow(M2a, s1)))
    tailsum <- .hx_add(tailsum, .hx_mul(.hx_mul(0.5, logM2ad), M2as))
    U <- .mv_int(1)
    r <- M2as
    fact <- gmp::as.bigz(2)
    for (j in seq_len(N)) {
      j2 <- 2 * j
      upds <- if (j == 1) 1 else c(j2 - 2, j2 - 1)
      for (m in upds) {
        U <- .hx_mul(.hx_sub(1 - m, s), U)
        r <- .hx_mul(r, rM2a)
      }
      # fdot(U, logs) with logs = [1]
      t <- .hx_div(.hx_mul(.hx_mul(.mv_pos(.mv_convert(U)), r), .zt_bernoulli(j2)), .mv_int(-fact))
      tailsum <- .hx_add(tailsum, t)
      if (.hx_mag(t) < tol) return(list(lsum, .hx_mul(1, tailsum)))
      fact <- fact * (j2 + 1) * (j2 + 2)
    }
    M1 <- M2; M2 <- M2 * 2
    if (.mv_cmp(.hx_re(s), 0) < 0) N <- N + N %/% 2
  }
}

# ctx._zetasum for derivative 0 without reflection: sum (a + k)^-s, k = 0..n
.zt_zetasum <- function(s, a, n) {
  if (.mv_cmp(.hx_abs(.hx_re(s)), 0.5 * .MPX$prec) < 0) {
    r <- .zt_zetasum_fast(s, a, n)
    if (!is.null(r)) return(r)
  }
  negs <- .hx_neg_exact(.mv_convert(s))
  .hx_fsum(lapply(seq_len(n + 1) - 1, function(k) .hx_pow(.hx_add(a, k), negs)))
}

# _zetasum_fast (mpc_zetasum): complex s and integer a, in fixed point
.zt_zetasum_fast <- function(s, a, n) {
  if (!(.hx_isint(a) && .mv(s)$k == "mpc")) return(NULL)
  a <- as.numeric(.hx_int_value(a))
  prec <- .MPX$prec; wp <- prec + 10
  critical_line <- .mp_cmp(s$re, .mpf(1, -1)) == 0
  sre <- .to_fixed(s$re, wp); sim <- .to_fixed(s$im, wp)
  one_2wp <- .pow2(2 * wp)
  xre <- gmp::as.bigz(0); xim <- gmp::as.bigz(0)
  for (w in seq(a, a + n)) {
    lg <- .zt_log_int(w, wp)
    cs <- .zt_cos_sin_fixed(.fx_shift(-sim * lg, -wp), wp)
    u <- if (critical_line) gmp::divq.bigz(one_2wp, .mp_isqrt(gmp::as.bigz(w) * one_2wp))
         else .zt_exp_fixed(.fx_shift(-sre * lg, -wp), wp)
    xre <- xre + .fx_shift(u * cs[[1L]], -wp)
    xim <- xim + .fx_shift(u * cs[[2L]], -wp)
  }
  .mv_mpc(.mp_round(xre, -wp, prec), .mp_round(xim, -wp, prec))
}

# ---------------------------------------------------------------------------
# polylog (defun_wrapped: ten extra bits, the result rounded)

.zt_polylog_series <- function(s, z) {
  tol <- .mv_pos(.hx_eps())
  l <- .hx_zero(); k <- 1; zk <- z
  repeat {
    term <- .hx_div(zk, .hx_pow(.mv_int(k), s))
    l <- .hx_add(l, term)
    if (isTRUE(.mv_cmp(.hx_abs(term), tol) < 0)) break
    zk <- .hx_mul(zk, z)
    k <- k + 1
  }
  l
}

.zt_twopij <- function() .hx_mul(.mv_mpc(.mp_zero, .mpf(1, 1)), .hx_pi())

.zt_polylog_continuation <- function(n, z) {
  if (n < 0) return(.hx_mul(z, 0))
  if (.hx_is_real_type(z) && .hx_isinf(z) && n > 0)
    return(if (.mv_cmp(z, 0) < 0) .mv_ninf() else .mv_new("mpc", re = .mp_zero, im = .mp_zero, sp = "nan"))
  twopij <- .zt_twopij()
  a <- .hx_mul(.hx_div(.hx_neg(.hx_pow(twopij, n)), .hx_factorial(n)),
               .zt_bernpoly(n, .hx_div(.hx_log(z), twopij)))
  if (.hx_is_real_type(z) && .mv_cmp(z, 0) < 0) a <- .hx_re(a)
  imz <- .hx_im(z)
  if (.mv_cmp(imz, 0) < 0 || (.hx_iszero(imz) && .mv_cmp(.hx_re(z), 1) >= 0))
    a <- .hx_sub(a, .hx_div(.hx_mul(twopij, .hx_pow(.hx_log(z), n - 1)), .hx_factorial(n - 1)))
  a
}

.zt_polylog_unitcircle <- function(n, z) {
  tol <- .mv_pos(.hx_eps())
  if (n > 1) {
    l <- .hx_zero(); logz <- .hx_log(z); logmz <- .hx_one(); m <- 0
    repeat {
      if (n - m != 1) {
        term <- .hx_div(.hx_mul(.zt_zeta(.mv_int(n - m)), logmz), .hx_factorial(m))
        if (!.hx_iszero(term) && .mv_cmp(.hx_abs(term), tol) < 0) break
        l <- .hx_add(l, term)
      }
      logmz <- .hx_mul(logmz, logz)
      m <- m + 1
    }
    l <- .hx_add(l, .hx_mul(.hx_div(.hx_pow(.hx_log(z), n - 1), .hx_factorial(n - 1)),
                            .hx_sub(.zt_harmonic(.mv_int(n - 1)), .hx_log(.hx_neg(.hx_log(z))))))
  } else if (n < 1) {
    l <- .hx_mul(.hx_factorial(-n), .hx_pow(.hx_neg(.hx_log(z)), n - 1))
    logz <- .hx_log(z); logkz <- .hx_one(); k <- 0
    repeat {
      b <- .zt_bernoulli(k - n + 1)
      if (!.hx_iszero(b)) {
        term <- .hx_div(.hx_mul(b, logkz), .hx_mul(.hx_factorial(k), k - n + 1))
        if (.mv_cmp(.hx_abs(term), tol) < 0) break
        l <- .hx_sub(l, term)
      }
      logkz <- .hx_mul(logkz, logz)
      k <- k + 1
    }
  } else .mp_error("ValueError", "")
  if (.hx_is_real_type(z) && .mv_cmp(z, 0) < 0) l <- .hx_re(l)
  l
}

.zt_polylog_general <- function(s, z) {
  v <- .hx_zero()
  u <- .hx_log(z)
  if (!isTRUE(.mv_cmp(.hx_abs(u), 5) < 0)) {
    j <- .hx_j()
    v <- .hx_sub(1, s)
    y <- .hx_div(.hx_log(.hx_neg(z)), .hx_mul(.hx_mul(2, .hx_pi()), j))
    A <- .hx_mul(.hx_pow(j, v), .zt_zeta(v, .hx_add(0.5, y)))
    B <- .hx_mul(.hx_pow(j, .hx_neg(v)), .zt_zeta(v, .hx_sub(0.5, y)))
    return(.hx_div(.hx_mul(.hx_gamma(v), .hx_add(A, B)), .hx_pow(.hx_mul(2, .hx_pi()), v)))
  }
  t <- .mv_int(1); k <- 0
  prec <- .MPX$prec
  on.exit(.MPX$prec <- prec)
  if (.hx_isfinite(s)) .MPX$prec <- .MPX$prec + max(0, -.hx_nint_distance(s)[[2L]])
  repeat {
    term <- .hx_mul(.zt_zeta(.hx_sub(s, k)), t)
    if (!isTRUE(.mv_cmp(.hx_abs(term), .hx_eps()) >= 0)) break
    v <- .hx_add(v, term)
    k <- k + 1
    t <- .hx_mul(t, u)
    t <- .hx_div(t, k)
  }
  .hx_add(.hx_mul(.hx_gamma(.hx_sub(1, s)), .hx_pow(.hx_neg(u), .hx_sub(s, 1))), v)
}

.zt_polylog <- function(s, z) {
  prec <- .MPX$prec
  on.exit(.MPX$prec <- prec)
  .MPX$prec <- prec + 10
  s <- .mv_convert(s); z <- .mv_convert(z)
  r <- if (isTRUE(.mv_eq(z, 1))) .zt_zeta(s)
  else if (isTRUE(.mv_eq(z, -1))) .hx_neg(.zt_altzeta(s))
  else if (isTRUE(.mv_eq(s, 0))) .hx_div(z, .hx_sub(1, z))
  else if (isTRUE(.mv_eq(s, 1))) .hx_neg(.hx_log(.hx_sub(1, z)))
  else if (isTRUE(.mv_eq(s, -1))) .hx_div(z, .hx_pow(.hx_sub(1, z), 2))
  else {
    az <- .hx_abs(z)
    sint <- .hx_isint(s)
    if (isTRUE(.mv_cmp(az, 0.75) <= 0) || (!sint && isTRUE(.mv_cmp(az, 0.9) < 0))) .zt_polylog_series(s, z)
    else if (isTRUE(.mv_cmp(az, 1.4) >= 0) && sint)
      .hx_add(.hx_mul(.hx_pow(.mv_int(-1), .hx_add(s, 1)), .zt_polylog_series(s, .hx_div(1, z))),
              .zt_polylog_continuation(as.numeric(.hx_int_value(.hx_re(s))), z))
    else if (.hx_isnan(z)) {
      if (.hx_is_real_type(z) && .hx_isnpint(s)) .mv_nan() else .mv_nan(TRUE)
    }
    else if (sint) .zt_polylog_unitcircle(as.numeric(.hx_int_value(.hx_re(s))), z)
    else .zt_polylog_general(s, z)
  }
  .MPX$prec <- prec
  .mv_pos(r)
}

# ---------------------------------------------------------------------------
# SymPy's dirichlet_eta(s)._eval_evalf(prec) for a rational s:
# rewrite(zeta)._eval_evalf(prec), the product (1 - 2**(1 - s))*zeta(s)
# in Float arithmetic. SymPy writes 2**x as 2**n * 2**r (n = floor(x),
# 0 < r < 1); the Floats are 1, 2**n * mpf_pow(2, r) and zeta(s) at prec.
.zt_eta_eval <- function(e, P, env) {
  if (length(e$args) != 1L) return(NULL)
  s <- e$args[[1L]]
  if (!(s$type %in% c("Integer", "Rational"))) {
    # a complex s: SymPy's Float expression keeps 2.0**(c*I) factors that
    # as_real_imag and evalf resolve; here the value is mpmath's altzeta(s)
    # at prec, which agrees with it to the printed digits
    a <- .mp_sympy_arg(s, P + 5, env)
    if (is.null(a)) return(NULL)
    return(.mp_workprec(P, .zt_altzeta(a)))
  }
  x <- 1 - gmp::as.bigq(s$value)
  if (x == 0) return(NULL)
  zv <- .MP_SYMPY[["zeta"]](fcall("zeta", list(s)), P, env)
  if (is.null(zv)) return(NULL)
  zv <- .mv(zv)
  if (.mv_sp(zv)) return(NULL)
  zv <- .mv_convert(zv)
  if (zv$k == "mpc" && zv$im$m != 0) return(NULL)
  Z <- .mp_round(zv$re$m, zv$re$e, P)
  if (gmp::denominator(x) == 1) {
    A <- .mp_from_q(1 - gmp::as.bigq(2)^as.numeric(x), P)
  } else {
    n <- gmp::numerator(x) %/% gmp::denominator(x)
    pw <- .zt_mpf_pow(.mpf(1, 1), .mp_from_q(x - n, P), P)
    A <- .mp_sub(.mpf(1, 0), .mp_ldexp(pw, as.numeric(n)), P)
  }
  .mv_mpf(.mp_mul(A, Z, P))
}

.mp_reg("zeta", .mp_generic(function(s, a = .mv_int(1)) .zt_zeta(s, a), 1:2))
.mp_reg("polylog", .mp_generic(function(s, z) .zt_polylog(s, z), 2L))
.mp_reg("dirichlet_eta", .zt_eta_eval)
