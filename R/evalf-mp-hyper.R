# Multiprecision hypergeometric functions after mpmath: the series summator
# of libmp/libhyper.py (make_hyp_summator), ctx.hypsum (ctx_mp.py), and
# hypercomb, hyper with its special cases and meijerg from
# functions/hypergeometric.py. Values are the mpmath numbers of
# evalf-mp-ctx.R, so integer and rational parameters stay exact as in mpmath.

.hx_default_maxprec <- function(p) as.integer(floor(1000 * p^0.25 + 4 * p))
.hx_kw <- function(kw, name, default = NULL) if (is.null(kw[[name]])) default else kw[[name]]

# ---------------------------------------------------------------------------
# _convert_param: a value and its type, Z (int), Q (MPQ with a power of two
# denominator at most 16), R (mpf) or C (mpc)

.hx_param <- function(x) {
  x <- .mv(x)
  if (x$k == "int") return(list(v = x, t = "Z"))
  # an exact (p, q) tuple
  if (x$k == "tuple") {
    if (gmp::denominator(x$q) == 1) return(list(v = .mv_int(gmp::numerator(x$q)), t = "Z"))
    return(list(v = .mv_new("mpq", q = x$q), t = "Q"))
  }
  if (x$k %in% c("mpq", "float")) x <- .mv_convert(x)
  if (.mv_sp(x)) .hx_valueerror(sprintf("mantissa and exponent are undefined for %s", x$sp))
  if (x$k == "mpc" && x$im$m != 0) return(list(v = x, t = "C"))
  v <- .mp_norm(x$re)
  if (v$m != 0) {
    if (v$e >= -4) {
      if (v$e >= 0) return(list(v = .mv_int(v$m * .pow2(v$e)), t = "Z"))
      return(list(v = .mv_new("mpq", q = gmp::as.bigq(v$m, .pow2(-v$e))), t = "Q"))
    }
    return(list(v = .mv_mpf(v), t = "R"))
  }
  list(v = .mv_int(0), t = "Z")
}
.hx_param_eq <- function(a, b) a$t == b$t && .mv_eq(a$v, b$v)

# ---------------------------------------------------------------------------
# make_hyp_summator: the series in fixed point with wp fractional bits

# Python's x >> k (floor) and x << k
.hx_shr <- function(x, k) if (k >= 0) gmp::divq.bigz(x, .pow2(k)) else x * .pow2(-k)
.hx_fix <- function(a, wp) .hx_shr(a$m, -(a$e + wp))

.hx_summator <- function(p, q, coeffs, types, z, prec, wp, epsshift, magcheck, maxterms) {
  zc <- z$k == "mpc"
  have_complex <- any(types == "C") || zc
  MAX <- if (is.null(maxterms)) wp * 100 else maxterms
  HIGH <- .pow2(epsshift); LOW <- -HIGH
  one <- .pow2(wp)
  SRE <- PRE <- one
  SIM <- PIM <- gmp::as.bigz(0)
  ZRE <- .hx_fix(z$re, wp)
  ZIM <- if (zc) .hx_fix(z$im, wp) else NULL
  n_par <- length(coeffs)
  isa <- seq_len(n_par) <= p
  INT <- vector("list", n_par); RP <- INT; RQ <- INT; REAL <- INT; CRE <- INT; CIM <- INT
  for (i in seq_len(n_par)) {
    c <- coeffs[[i]]
    switch(types[i],
      Z = INT[[i]] <- gmp::numerator(c$q),
      Q = { RP[[i]] <- gmp::numerator(c$q); RQ[[i]] <- gmp::denominator(c$q) },
      R = REAL[[i]] <- .hx_fix(c$re, wp),
      C = { CRE[[i]] <- .hx_fix(c$re, wp); CIM[[i]] <- .hx_fix(c$im, wp) })
  }
  idx <- function(tt, a) which(types == tt & isa == a)
  aint <- idx("Z", TRUE); bint <- idx("Z", FALSE)
  arat <- idx("Q", TRUE); brat <- idx("Q", FALSE)
  areal <- idx("R", TRUE); breal <- idx("R", FALSE)
  acomplex <- idx("C", TRUE); bcomplex <- idx("C", FALSE)
  nc <- min(length(areal), length(breal))
  ncnum <- areal[seq_along(areal) > nc]; ncden <- breal[seq_along(breal) > nc]
  has_mul <- length(aint) + length(arat) + length(brat) > 0
  mkeys <- names(magcheck)
  n <- 0
  while (n < 1e8 - 1) {
    n <- n + 1
    if (length(mkeys)) {
      key <- format(n, scientific = FALSE)
      if (key %in% mkeys) {
        pm <- .bitlen(abs(PRE))
        if (have_complex) pm <- max(pm, .bitlen(abs(PIM)))
        magcheck[[key]] <- wp - pm
      }
    }
    mul <- gmp::as.bigz(1)
    for (i in aint) mul <- mul * INT[[i]]
    for (i in arat) mul <- mul * RP[[i]]
    for (i in brat) mul <- mul * RQ[[i]]
    div <- gmp::as.bigz(1)
    for (i in bint) div <- div * INT[[i]]
    for (i in brat) div <- div * RP[[i]]
    for (i in arat) div <- div * RQ[[i]]
    div <- div * n
    if (div == 0) {
      if (has_mul && mul == 0) break
      .mp_error("ZeroDivisionError", "")
    }
    if (have_complex) {
      for (k in seq_len(nc)) PRE <- gmp::divq.bigz(PRE * REAL[[areal[k]]], REAL[[breal[k]]])
      for (i in ncnum) PRE <- .hx_shr(PRE * REAL[[i]], wp)
      for (i in ncden) PRE <- gmp::divq.bigz(PRE * one, REAL[[i]])
      for (k in seq_len(nc)) PIM <- gmp::divq.bigz(PIM * REAL[[areal[k]]], REAL[[breal[k]]])
      for (i in ncnum) PIM <- .hx_shr(PIM * REAL[[i]], wp)
      for (i in ncden) PIM <- gmp::divq.bigz(PIM * one, REAL[[i]])
      if (has_mul) {
        if (zc) {
          a <- gmp::divq.bigz(mul * (PRE * ZRE - PIM * ZIM), div)
          b <- gmp::divq.bigz(mul * (PIM * ZRE + PRE * ZIM), div)
          PRE <- .hx_shr(a, wp); PIM <- .hx_shr(b, wp)
        } else {
          PRE <- gmp::divq.bigz(.hx_shr(mul * PRE * ZRE, wp), div)
          PIM <- gmp::divq.bigz(.hx_shr(mul * PIM * ZRE, wp), div)
        }
      } else {
        if (zc) {
          a <- gmp::divq.bigz(PRE * ZRE - PIM * ZIM, div)
          b <- gmp::divq.bigz(PIM * ZRE + PRE * ZIM, div)
          PRE <- .hx_shr(a, wp); PIM <- .hx_shr(b, wp)
        } else {
          PRE <- gmp::divq.bigz(.hx_shr(PRE * ZRE, wp), div)
          PIM <- gmp::divq.bigz(.hx_shr(PIM * ZRE, wp), div)
        }
      }
      for (i in acomplex) {
        a <- PRE * CRE[[i]] - PIM * CIM[[i]]
        b <- PIM * CRE[[i]] + PRE * CIM[[i]]
        PRE <- .hx_shr(a, wp); PIM <- .hx_shr(b, wp)
      }
      for (i in bcomplex) {
        mg <- CRE[[i]] * CRE[[i]] + CIM[[i]] * CIM[[i]]
        re <- PRE * CRE[[i]] + PIM * CIM[[i]]
        im <- PIM * CRE[[i]] - PRE * CIM[[i]]
        PRE <- gmp::divq.bigz(re * one, mg)
        PIM <- gmp::divq.bigz(im * one, mg)
      }
      SRE <- SRE + PRE
      SIM <- SIM + PIM
      if (PRE < HIGH && PRE > LOW && PIM < HIGH && PIM > LOW) break
    } else {
      for (k in seq_len(nc)) PRE <- gmp::divq.bigz(PRE * REAL[[areal[k]]], REAL[[breal[k]]])
      for (i in ncnum) PRE <- .hx_shr(PRE * REAL[[i]], wp)
      for (i in ncden) PRE <- gmp::divq.bigz(PRE * one, REAL[[i]])
      PRE <- if (has_mul) gmp::divq.bigz(.hx_shr(PRE * mul * ZRE, wp), div)
             else gmp::divq.bigz(.hx_shr(PRE * ZRE, wp), div)
      SRE <- SRE + PRE
      if (PRE < HIGH && PRE > LOW) break
    }
    if (n > MAX) .hx_noconv("Hypergeometric series converges too slowly. Try increasing maxterms.")
    for (i in c(aint, bint)) INT[[i]] <- INT[[i]] + 1
    for (i in c(arat, brat)) RP[[i]] <- RP[[i]] + RQ[[i]]
    for (i in c(areal, breal)) REAL[[i]] <- REAL[[i]] + one
    for (i in c(acomplex, bcomplex)) CRE[[i]] <- CRE[[i]] + one
  }
  a <- .mp_round(SRE, -wp, prec)
  if (have_complex) {
    b <- .mp_round(SIM, -wp, prec)
    magn <- if (SRE != 0) { if (SIM != 0) max(.mp_mag(a), .mp_mag(b)) else .mp_mag(a) }
            else if (SIM != 0) .mp_mag(b) else -wp + 1
    return(list(v = .mv_mpc(a, b), mag = magn, magcheck = magcheck))
  }
  magn <- if (SRE != 0) .mp_mag(a) else -wp + 1
  list(v = .mv_mpf(a), mag = magn, magcheck = magcheck)
}

.hx_hypsum_msg <- function(prec, wp) sprintf(paste0("hypsum() failed to converge to the requested %d bits of accuracy\n",
  "using a working precision of %d bits. Try with a higher maxprec,\n",
  "maxterms, or set zeroprec."), as.integer(prec), as.integer(wp))

# ctx.hypsum(p, q, flags, coeffs, z, accurate_small=True, **kwargs)
.hx_hypsum <- function(p, q, types, coeffs, z, kw = list()) {
  coeffs <- lapply(coeffs, .mv)
  z <- .mv_convert(z)
  accurate_small <- .hx_kw(kw, "accurate_small", TRUE)
  for (i in seq_along(coeffs)) {
    if (i <= p || types[i] != "Z") next
    c <- coeffs[[i]]
    if (c <= 0) {
      ok <- FALSE
      for (ii in seq_len(p)) {
        cc <- coeffs[[ii]]
        if (types[ii] == "Z" && cc <= 0 && c <= cc) ok <- TRUE
      }
      if (!ok) .hx_zerodiv("pole in hypergeometric series")
    }
  }
  if (.hx_isinf(z)) {
    nn <- -1L; best <- 0
    for (i in seq_len(p)) {
      if (types[i] == "Z" && coeffs[[i]] < 0) {
        if (nn < 0 || coeffs[[i]] < best) { nn <- i; best <- coeffs[[i]] }
      }
    }
    if (nn >= 0) {
      n <- -gmp::numerator(coeffs[[nn]]$q)
      t <- .hx_pow(z, .mv_int(n))
      k <- 0
      while (k < n) {
        for (i in seq_len(p)) t <- .hx_mul(t, .hx_add(coeffs[[i]], k))
        for (i in p + seq_len(q)) t <- .hx_div(t, .hx_add(coeffs[[i]], k))
        t <- .hx_div(t, k + 1)
        k <- k + 1
      }
      return(t)
    }
  }
  prec <- .MPX$prec
  maxprec <- .hx_kw(kw, "maxprec", .hx_default_maxprec(prec))
  extraprec <- 50
  epsshift <- 25
  magnitude_check <- list()
  max_total_jump <- 0
  for (i in seq_along(coeffs)) {
    if (types[i] == "Z") next
    nd <- .hx_nint_distance(coeffs[[i]])
    n <- -nd$n
    d <- -nd$d
    if (i > p && n >= 0 && d > 4) {
      key <- as.character(n)
      magnitude_check[[key]] <- (if (is.null(magnitude_check[[key]])) 0 else magnitude_check[[key]]) + d
      extraprec <- max(extraprec, d - prec + 60)
    }
    max_total_jump <- max_total_jump + abs(d)
  }
  repeat {
    if (extraprec > maxprec) .hx_noconv(.hx_hypsum_msg(prec, prec + extraprec))
    wp <- prec + extraprec
    mag_dict <- if (length(magnitude_check)) stats::setNames(as.list(rep(NA_real_, length(magnitude_check))), names(magnitude_check)) else list()
    r <- .hx_summator(p, q, coeffs, types, z, prec, wp, epsshift, mag_dict, kw$maxterms)
    cancel <- -r$mag
    jumps_resolved <- TRUE
    if (extraprec < max_total_jump) {
      for (v in r$magcheck) if (is.na(v) || v < prec) { jumps_resolved <- FALSE; break }
    }
    accurate <- cancel < extraprec - 25 - 5 || !accurate_small
    if (jumps_resolved) {
      if (accurate) break
      zeroprec <- kw$zeroprec
      if (!is.null(zeroprec) && cancel > zeroprec) {
        return(if (r$v$k == "mpc") .mv_mpc(.mp_zero, .mp_zero) else .hx_zero())
      }
    }
    extraprec <- extraprec * 2
    epsshift <- epsshift + 5
    extraprec <- extraprec + 5
  }
  r$v
}

# ---------------------------------------------------------------------------
# hypercomb: a sum of terms, each a product of powers, gamma ratios and a
# hypergeometric series; parameters are perturbed past gamma poles

.hx_term <- function(w, c, alpha, beta, a, b, z)
  list(w = lapply(w, .mv), c = lapply(c, .mv), alpha = lapply(alpha, .mv), beta = lapply(beta, .mv),
       a = lapply(a, .mv), b = lapply(b, .mv), z = .mv(z))

.hx_check_need_perturb <- function(terms, prec, discard_known_zeros) {
  perturb <- recompute <- FALSE
  extraprec <- 0
  discard <- integer()
  for (ti in seq_along(terms)) {
    term <- terms[[ti]]
    have_singular_nongamma_weight <- FALSE
    for (k in seq_along(term$w)) {
      w <- term$w[[k]]
      if (.hx_iszero(w)) {
        ck <- term$c[[k]]
        if (.hx_re(.mv_convert(ck)) <= 0 && !.hx_iszero(ck)) {
          perturb <- recompute <- TRUE
          have_singular_nongamma_weight <- TRUE
        }
      }
    }
    pole_count <- c(0, 0, 0)
    groups <- list(term$alpha, term$beta, term$b)
    for (di in 1:3) {
      for (x in groups[[di]]) {
        nd <- .hx_nint_distance(x)
        if (nd$n > 0) next
        if (nd$d == -Inf) {
          ok <- FALSE
          if (di == 3) {
            for (u in term$a) if (.hx_isnpint(u) && u >= .mv_int(nd$n)) { ok <- TRUE; break }
          }
          if (ok) next
          pole_count[di] <- pole_count[di] + 1
        } else if (nd$d < -4) {
          extraprec <- extraprec - nd$d
          recompute <- TRUE
        }
      }
    }
    if (discard_known_zeros && pole_count[2] > pole_count[1] + pole_count[3] && !have_singular_nongamma_weight)
      discard <- c(discard, ti)
    else if (sum(pole_count)) perturb <- recompute <- TRUE
  }
  list(perturb = perturb, recompute = recompute, extraprec = extraprec, discard = discard)
}

.hx_hypercomb_msg <- function(orig, prec) sprintf(paste0("\nhypercomb() failed to converge to the requested %d bits of accuracy\n",
  "using a working precision of %d bits. The function value may be zero or\n",
  "infinite; try passing zeroprec=N or infprec=M to bound finite values between\n",
  "2^(-N) and 2^M. Otherwise try a higher maxprec or maxterms.\n"), as.integer(orig), as.integer(prec))

# fn takes the parameters as separate arguments and returns a list of
# terms made by .hx_term
.hx_hypercomb <- function(fn, params = list(), discard_known_zeros = TRUE, kw = list()) {
  orig <- .MPX$prec
  on.exit(.MPX$prec <- orig)
  sumvalue <- .hx_zero()
  orig_params <- lapply(params, .mv)
  maxprec <- .hx_kw(kw, "maxprec", .hx_default_maxprec(orig))
  kw$maxprec <- maxprec
  zeroprec <- kw$zeroprec
  infprec <- kw$infprec
  perturbed_reference_value <- NULL
  hextra <- 0
  repeat {
    .MPX$prec <- .MPX$prec + 10
    if (.MPX$prec > maxprec) .hx_valueerror(.hx_hypercomb_msg(orig, .MPX$prec))
    orig2 <- .MPX$prec
    params <- orig_params
    terms <- do.call(fn, params)
    chk <- .hx_check_need_perturb(terms, orig, discard_known_zeros)
    .MPX$prec <- .MPX$prec + chk$extraprec
    if (chk$perturb) {
      hmag <- if (!is.null(kw$hmag)) kw$hmag else orig + 10 + hextra
      h <- .hx_ldexp(.hx_one(), -hmag)
      .MPX$prec <- orig2 + 10 + hmag + 10
      for (k in seq_along(params)) {
        params[[k]] <- .hx_add(params[[k]], h)
        h <- .hx_add(h, .hx_div(h, k))
      }
    }
    if (chk$recompute) terms <- do.call(fn, params)
    if (discard_known_zeros && length(chk$discard)) terms <- terms[-chk$discard]
    if (!length(terms)) return(.hx_zero())
    evaluated <- lapply(terms, function(t) {
      .hx_fprod(c(list(.hx_hyper(t$a, t$b, t$z, kw)),
                  lapply(t$alpha, .hx_gamma), lapply(t$beta, .hx_rgamma),
                  Map(.hx_power, t$w, t$c)))
    })
    if (length(terms) == 1L && !chk$perturb) { sumvalue <- evaluated[[1L]]; break }
    sumvalue <- .hx_fsum(evaluated)
    mags <- vapply(evaluated, function(x) as.numeric(.hx_mag(x)), 0)
    max_magnitude <- max(mags)
    sum_magnitude <- as.numeric(.hx_mag(sumvalue))
    cancellation <- max_magnitude - sum_magnitude
    precision_ok <- isTRUE(cancellation < .MPX$prec - orig)
    zero_ok <- !is.null(zeroprec) && max_magnitude - .MPX$prec < -zeroprec
    inf_ok <- !is.null(infprec) && max_magnitude > infprec
    if ((precision_ok && !chk$perturb) || is.nan(cancellation)) break
    if (precision_ok) {
      if (is.null(perturbed_reference_value)) {
        hextra <- hextra + 20
        perturbed_reference_value <- sumvalue
        next
      } else if (.hx_mag(.hx_sub(sumvalue, perturbed_reference_value)) <= .hx_mag(sumvalue) - orig) {
        break
      } else if (zero_ok) {
        sumvalue <- .hx_zero(); break
      } else if (inf_ok) {
        sumvalue <- .mv_inf(); break
      } else if (!is.null(kw$hmag)) {
        break
      } else {
        hextra <- hextra * 2
        perturbed_reference_value <- sumvalue
      }
    } else {
      increment <- min(max(cancellation, orig %/% 2), max(chk$extraprec, orig))
      .MPX$prec <- .MPX$prec + increment
      next
    }
  }
  .MPX$prec <- orig
  .mv_pos(sumvalue)
}

# ctx.power(x, y) = convert(x) ** convert(y)
.hx_power <- function(x, y) .hx_pow(.mv_convert(x), .mv_convert(y))

# ---------------------------------------------------------------------------
# hyper and its special cases

.hx_hyper <- function(a_s, b_s, z, kw = list()) {
  z <- .mv_convert(z)
  if (.hx_isnan(z)) return(.mv_nan())
  ap <- lapply(a_s, .hx_param); bq <- lapply(b_s, .hx_param)
  p <- length(ap); q <- length(bq)
  if (!isFALSE(kw$eliminate)) {
    elim_nonpositive <- isTRUE(kw$eliminate_all)
    i <- 1L
    while (i <= q && length(ap)) {
      b <- bq[[i]]
      hit <- which(vapply(ap, .hx_param_eq, NA, b))
      if (length(hit) && (elim_nonpositive || !.hx_isnpint(b$v))) {
        ap <- ap[-hit[1L]]
        hb <- which(vapply(bq, .hx_param_eq, NA, b))
        bq <- bq[-hb[1L]]
        p <- p - 1L; q <- q - 1L
      } else i <- i + 1L
    }
  }
  if (p == 0) {
    if (q == 1) return(.hx_hyp0f1(bq, z, kw))
    if (q == 0) return(.hx_exp(z))
  } else if (p == 1) {
    if (q == 1) return(.hx_hyp1f1(ap, bq, z, kw))
    if (q == 2) return(.hx_hyp1f2(ap, bq, z, kw))
    if (q == 0) return(.hx_hyp1f0(ap[[1L]]$v, z))
  } else if (p == 2) {
    if (q == 1) return(.hx_hyp2f1(ap, bq, z, kw))
    if (q == 2) return(.hx_hyp2f2(ap, bq, z, kw))
    if (q == 3) return(.hx_hyp2f3(ap, bq, z, kw))
    if (q == 0) return(.hx_hyp2f0(ap, bq, z, kw))
  } else if (p == q + 1) {
    return(.hx_hypq1fq(p, q, ap, bq, z, kw))
  } else if (p > q + 1 && !isTRUE(kw$force_series)) {
    return(.hx_hyp_borel(p, q, ap, bq, z, kw))
  }
  pq <- c(ap, bq)
  .hx_hypsum(p, q, vapply(pq, `[[`, "", "t"), lapply(pq, `[[`, "v"), z, kw)
}
.hx_vals <- function(ps) lapply(ps, `[[`, "v")
.hx_types <- function(ps) vapply(ps, `[[`, "", "t")

# defun_wrapped: arguments converted, ten guard bits, the result rounded
.hx_wrapped <- function(f, ...) {
  args <- lapply(list(...), .mv_convert)
  prec <- .MPX$prec
  .MPX$prec <- prec + 10
  r <- tryCatch(do.call(f, args), finally = .MPX$prec <- prec)
  .mv_pos(r)
}

.hx_hyp1f0 <- function(a, z) .hx_wrapped(function(a, z) .hx_pow(.hx_sub(1, z), .hx_neg(a)), a, z)

.hx_hyp0f1 <- function(bq, z, kw) {
  b <- bq[[1L]]$v; btype <- bq[[1L]]$t
  magz <- if (!.hx_iszero(z)) .hx_mag(z) else 0
  if (magz >= 8 && !isTRUE(kw$force_series)) {
    orig <- .MPX$prec
    v <- tryCatch({
      .MPX$prec <- .MPX$prec + 12 + magz %/% 2
      h <- function() {
        w <- .hx_sqrt(.hx_neg(z))
        jw <- .hx_mul(.hx_j(), w)
        u <- .hx_div(1, .hx_mul(4, jw))
        c <- .hx_sub(.mv_q(1, 2), b)
        E <- .hx_exp(.hx_mul(2, jw))
        list(.hx_term(list(.hx_neg(jw), E), list(c, -1), list(), list(), list(.hx_sub(b, .mv_q(1, 2)), .hx_sub(.mv_q(3, 2), b)), list(), .hx_neg(u)),
             .hx_term(list(jw, E), list(c, 1), list(), list(), list(.hx_sub(b, .mv_q(1, 2)), .hx_sub(.mv_q(3, 2), b)), list(), u))
      }
      v <- .hx_hypercomb(h, list(), kw = list(force_series = TRUE))
      .hx_mul(.hx_div(.hx_gamma(b), .hx_mul(2, .hx_sqrt(.hx_pi()))), v)
    }, NoConvergence = function(err) NULL, finally = .MPX$prec <- orig)
    if (!is.null(v)) {
      if (.hx_is_real_type(b) && .hx_is_real_type(z)) v <- .hx_re(v)
      return(.mv_pos(v))
    }
  }
  .hx_hypsum(0, 1, btype, list(b), z, kw)
}

.hx_hyp1f1 <- function(ap, bq, z, kw) {
  a <- ap[[1L]]$v; b <- bq[[1L]]$v
  if (.hx_iszero(z)) return(.hx_add(.hx_one(), z))
  magz <- .hx_mag(z)
  if (magz >= 7 && !(.hx_isint(a) && .hx_re(.mv_convert(a)) <= 0)) {
    if (.hx_isinf(z) && .hx_sign(a) == 1 && .hx_sign(b) == 1 && .hx_sign(z) == 1) return(.mv_inf())
    if (is.infinite(magz)) magz <- 0
    orig <- .MPX$prec
    v <- tryCatch({
      .MPX$prec <- orig + magz
      sector <- .hx_im0(z) < 0
      h <- function(a, b) {
        E <- if (sector) .hx_expjpi(.hx_neg_exact(a)) else .hx_expjpi(a)
        rz <- .hx_div(1, z)
        list(.hx_term(list(E, z), list(1, .hx_neg(a)), list(b), list(.hx_sub(b, a)), list(a, .hx_sub(.hx_add(1, a), b)), list(), .hx_neg(rz)),
             .hx_term(list(.hx_exp(z), z), list(1, .hx_sub(a, b)), list(b), list(a), list(.hx_sub(b, a), .hx_sub(1, a)), list(), rz))
      }
      v <- .hx_hypercomb(h, list(a, b), kw = list(force_series = TRUE))
      if (.hx_is_real_type(a) && .hx_is_real_type(b) && .hx_is_real_type(z)) v <- .hx_re(v)
      .mv_pos(v)
    }, NoConvergence = function(err) NULL, finally = .MPX$prec <- orig)
    if (!is.null(v)) return(v)
  }
  .hx_hypsum(1, 1, c(ap[[1L]]$t, bq[[1L]]$t), list(a, b), z, kw)
}

.hx_hyp2f1_gosper <- function(a, b, c, z, kw) {
  .z <- z
  orig <- .MPX$prec
  maxprec <- .hx_kw(kw, "maxprec", 100 * orig)
  extra <- 10
  repeat {
    .MPX$prec <- orig + extra
    z <- .mv_convert(.z)
    d <- .hx_zero(); e <- .hx_one(); f <- .hx_zero()
    k <- 0
    abz <- .hx_mul(.hx_mul(a, b), z)
    ch <- .hx_mul(c, .mv_q(1, 2))
    c1h <- .hx_mul(.hx_add(c, 1), .mv_q(1, 2))
    nz <- .hx_sub(1, z)
    g <- .hx_div(z, nz)
    abg <- .hx_mul(.hx_mul(a, b), g)
    cba <- .hx_sub(.hx_sub(c, b), a)
    z2 <- .hx_sub(z, 2)
    tol <- -.MPX$prec - 10
    maxmag <- -Inf
    repeat {
      kch <- .hx_add(k, ch)
      kakbz <- .hx_div(.hx_mul(.hx_mul(.hx_add(k, a), .hx_add(k, b)), z),
                       .hx_mul(.hx_mul(.hx_mul(4, k + 1), kch), .hx_add(k, c1h)))
      d1 <- .hx_mul(kakbz, .hx_sub(e, .hx_mul(.hx_mul(.hx_add(k, cba), d), g)))
      e1 <- .hx_mul(kakbz, .hx_add(.hx_mul(d, abg), .hx_mul(.hx_add(k, c), e)))
      inner <- .hx_sub(.hx_mul(k, .hx_sub(.hx_add(.hx_mul(cba, z), .hx_mul(k, z2)), c)), abz)
      ft <- .hx_div(.hx_mul(d, inner), .hx_mul(.hx_mul(2, kch), nz))
      f1 <- .hx_sub(.hx_add(f, e), ft)
      maxmag <- max(maxmag, .hx_mag(f1))
      if (.hx_mag(.hx_sub(f1, f)) < tol) break
      d <- d1; e <- e1; f <- f1
      k <- k + 1
    }
    cancellation <- maxmag - .hx_mag(f1)
    if (cancellation < extra) break
    extra <- extra + cancellation
    if (extra > maxprec) .hx_noconv("")
  }
  f1
}

.hx_hyp2f1 <- function(ap, bq, z, kw) {
  a <- ap[[1L]]$v; b <- ap[[2L]]$v; c <- bq[[1L]]$v
  isintle0 <- function(x) .hx_isint(x) && x <= 0
  if (z == 1) {
    convergent <- .hx_re(.mv_convert(.hx_sub(.hx_sub(c, a), b))) > 0
    finite <- isintle0(a) || isintle0(b)
    zerodiv <- .hx_isint(c) && c <= 0 &&
      !((.hx_isint(a) && c <= a && a <= 0) || (.hx_isint(b) && c <= b && b <= 0))
    if ((convergent || finite) && !zerodiv)
      return(.hx_gammaprod(list(c, .hx_sub(.hx_sub(c, a), b)), list(.hx_sub(c, a), .hx_sub(c, b)), infsign = TRUE))
    return(.hx_mul(.hx_hyper(list(a, b), list(c), .hx_sub(1, .hx_mul(.hx_eps(), 2))), .mv_inf()))
  }
  if (.hx_iszero(z)) {
    if (!.hx_iszero(c) || .hx_iszero(a) || .hx_iszero(b)) return(.hx_add(1, z))
    return(.mv_nan())
  }
  if (.hx_isint(c) && c <= 0) {
    if (!((.hx_isint(a) && c <= a && a <= 0) || (.hx_isint(b) && c <= b && b <= 0))) return(.mv_inf())
  }
  absz <- .hx_abs(z)
  if (.hx_isfinite(z) && (absz <= 0.8 || (.hx_isint(a) && -1000 <= a && a <= 0) ||
                          (.hx_isint(b) && -1000 <= b && b <= 0))) {
    r <- tryCatch(.hx_hypsum(2, 1, c(ap[[1L]]$t, ap[[2L]]$t, bq[[1L]]$t), list(a, b, c), z, kw),
                  NoConvergence = function(err) { if (isTRUE(kw$force_series)) stop(err); NULL })
    if (!is.null(r)) return(r)
  }
  orig <- .MPX$prec
  v <- tryCatch({
    .MPX$prec <- orig + 10
    if (absz >= 1.3) {
      h <- function(a, b) {
        t <- .hx_sub(.mv_q(1), c); ab <- .hx_sub(a, b); rz <- .hx_div(1, z)
        list(.hx_term(list(.hx_neg(z)), list(.hx_neg(a)), list(c, .hx_neg(ab)), list(b, .hx_sub(c, a)), list(a, .hx_add(t, a)), list(.hx_add(.mv_q(1), ab)), rz),
             .hx_term(list(.hx_neg(z)), list(.hx_neg(b)), list(c, ab), list(a, .hx_sub(c, b)), list(b, .hx_add(t, b)), list(.hx_sub(.mv_q(1), ab)), rz))
      }
      .hx_hypercomb(h, list(a, b), kw = kw)
    } else if (.hx_abs(.hx_sub(1, z)) <= 0.75) {
      h <- function(a, b) {
        t <- .hx_sub(.hx_sub(c, a), b); ca <- .hx_sub(c, a); cb <- .hx_sub(c, b); rz <- .hx_sub(1, z)
        list(.hx_term(list(), list(), list(c, t), list(ca, cb), list(a, b), list(.hx_sub(1, t)), rz),
             .hx_term(list(rz), list(t), list(c, .hx_sub(.hx_add(a, b), c)), list(a, b), list(ca, cb), list(.hx_add(1, t)), rz))
      }
      .hx_hypercomb(h, list(a, b), kw = kw)
    } else if (.hx_abs(.hx_div(z, .hx_sub(z, 1))) <= 0.75) {
      .hx_div(.hx_hyper(list(a, .hx_sub(c, b)), list(c), .hx_div(z, .hx_sub(z, 1))), .hx_pow(.hx_sub(1, z), a))
    } else {
      .hx_hyp2f1_gosper(a, b, c, z, kw)
    }
  }, finally = .MPX$prec <- orig)
  .mv_pos(v)
}

.hx_hypq1fq <- function(p, q, ap, bq, z, kw) {
  a_s <- .hx_vals(ap); b_s <- .hx_vals(bq)
  types <- c(.hx_types(ap), .hx_types(bq))
  absz <- .hx_abs(z)
  ispoly <- FALSE
  for (a in a_s) if (.hx_isint(a) && a <= 0) { ispoly <- TRUE; break }
  if (absz < 1 || ispoly) {
    r <- tryCatch(.hx_hypsum(p, q, types, c(a_s, b_s), z, kw),
                  NoConvergence = function(err) { if (absz > 1.1 || ispoly) stop(err); NULL })
    if (!is.null(r)) return(r)
  }
  if (z == 1) {
    S <- .hx_re(.mv_convert(.hx_sub(.hx_sum(b_s), .hx_sum(a_s))))
    if (S <= 0) return(.hx_mul(.hx_hyper(a_s, b_s, 0.9, kw), .mv_inf()))
  }
  if (p == 3 && q == 2 && .hx_abs(.hx_sub(z, 1)) < 0.05) {
    a1 <- a_s[[1L]]; a2 <- a_s[[2L]]; a3 <- a_s[[3L]]
    b1 <- b_s[[1L]]; b2 <- b_s[[2L]]
    u <- .hx_sub(.hx_add(b1, b2), a3)
    initial <- .hx_gammaprod(list(.hx_sub(b2, a3), .hx_sub(b1, a3), a1, a2), list(.hx_sub(b2, a3), .hx_sub(b1, a3), 1, u))
    cache <- new.env(parent = emptyenv())
    cache[["0"]] <- initial
    term <- function(k) {
      u <- .hx_add(.hx_sub(.hx_add(b1, b2), a3), k)
      # the cache is a dict keyed by the integer k
      kn <- .hx_int_value(k)
      if (!isTRUE(.mv_eq(k, .mv_int(kn)))) .mp_error("KeyError", format(.mv(k)))
      key <- as.character(kn)
      if (!is.null(cache[[key]])) {
        t <- cache[[key]]
      } else {
        t <- cache[[as.character(kn - 1)]]
        if (is.null(t)) .mp_error("KeyError", as.character(kn - 1))
        t <- .hx_mul(t, .hx_mul(.hx_sub(.hx_sub(.hx_add(b1, k), a3), 1), .hx_sub(.hx_sub(.hx_add(b2, k), a3), 1)))
        t <- .hx_div(t, .hx_mul(k, .hx_sub(u, 1)))
        cache[[key]] <- t
      }
      .hx_mul(t, .hx_hyper(list(a1, a2), list(u), z))
    }
    S <- tryCatch(.mp_nsum(term, list(0, .mv_inf()), strict = .hx_kw(kw, "strict", TRUE)),
                  NoConvergence = function(err) NULL)
    if (!is.null(S)) return(.hx_mul(S, .hx_gammaprod(list(b1, b2), list(a1, a2, a3))))
  }
  if (absz < 1.1 && .hx_re(z) <= 1) {
    cache <- new.env(parent = emptyenv())
    cache[["0"]] <- .hx_one()
    term <- function(kk) {
      kk <- .mv(kk)
      k <- .hx_int_value(kk)
      if (!(kk == .mv_int(k))) {
        t <- .hx_div(.hx_pow(z, .mv_convert(kk)), .hx_factorial(kk))
        for (a in a_s) t <- .hx_mul(t, .mp_rf(a, kk))
        for (b in b_s) t <- .hx_div(t, .mp_rf(b, kk))
        return(t)
      }
      key <- as.character(k)
      if (!is.null(cache[[key]])) return(cache[[key]])
      t <- term(k - 1)
      m <- k - 1
      for (j in seq_len(p)) t <- .hx_mul(t, .hx_add(a_s[[j]], m))
      for (j in seq_len(q)) t <- .hx_div(t, .hx_add(b_s[[j]], m))
      t <- .hx_mul(t, z)
      t <- .hx_div(t, k)
      cache[[key]] <- t
      t
    }
    sum_method <- .hx_kw(kw, "sum_method", "r+s+e")
    r <- tryCatch(.mp_nsum(term, list(0, .mv_inf()), strict = .hx_kw(kw, "strict", TRUE),
                           method = gsub("e", "", sum_method, fixed = TRUE)),
                  NoConvergence = function(err) { if (!grepl("e", sum_method, fixed = TRUE)) stop(err); NULL })
    if (!is.null(r)) return(r)
    log_diffs <- function(k0) {
      b2 <- c(b_s, list(1))
      i <- -1
      function() {
        i <<- i + 1
        if (i == 0) {
          s <- .hx_zero()
          v <- .hx_sum(lapply(a_s, function(a) .hx_loggamma(.hx_add(a, k0))))
          v <- .hx_sub(v, .hx_sum(lapply(b2, function(b) .hx_loggamma(.hx_add(b, k0)))))
          return(.hx_add(v, .hx_mul(k0, .hx_log(z))))
        }
        j <- i - 1
        v <- .hx_sub(.hx_sum(lapply(a_s, function(a) .mp_psi(j, .hx_add(a, k0)))),
                     .hx_sum(lapply(b2, function(b) .mp_psi(j, .hx_add(b, k0)))))
        if (j == 0) v <- .hx_add(v, .hx_log(z))
        v
      }
    }
    hyper_diffs <- function(k0) {
      C <- .hx_gammaprod(b_s, a_s)
      gen <- .mp_diffs_exp(log_diffs(k0))
      function() .hx_mul(C, gen())
    }
    tol <- .hx_div(.hx_eps(), 1024)
    prec <- .MPX$prec
    v <- tryCatch({
      trunc <- 50 * .mp_dps(.MPX$prec)
      .MPX$prec <- .MPX$prec + 20
      out <- NULL
      for (i in 0:4) {
        head <- .hx_fsum(lapply(seq_len(trunc) - 1, term))
        te <- .mp_sumem(term, list(trunc, .mv_inf()), tol = tol, adiffs = hyper_diffs(trunc), error = TRUE, fast_abort = TRUE)
        if (te$err < tol) { out <- .hx_add(head, te$value); break }
        trunc <- trunc * 2
        .MPX$prec <- .MPX$prec + .MPX$prec %/% 2
        if (i == 4) .hx_noconv("Euler-Maclaurin summation did not converge")
      }
      out
    }, finally = .MPX$prec <- prec)
    return(.mv_pos(v))
  }
  h <- function(...) {
    args <- list(...)
    a_s <- args[seq_len(p)]; b_s <- args[p + seq_len(q)]
    recz <- .hx_div(.hx_one(), z)
    negz <- .hx_neg_exact(z)
    lapply(seq_len(q + 1), function(k) {
      ak <- a_s[[k]]
      others <- setdiff(seq_len(q + 1), k)
      .hx_term(list(negz), list(.hx_neg(ak)),
               c(b_s, list(ak), lapply(others, function(j) .hx_sub(a_s[[j]], ak))),
               c(a_s, lapply(seq_len(q), function(j) .hx_sub(b_s[[j]], ak))),
               c(list(ak), lapply(seq_len(q), function(j) .hx_add(.hx_sub(ak, b_s[[j]]), 1))),
               lapply(others, function(j) .hx_add(.hx_sub(1, a_s[[j]]), ak)),
               recz)
    })
  }
  .hx_hypercomb(h, c(a_s, b_s), kw = kw)
}

.hx_hyp_borel <- function(p, q, ap, bq, z, kw) {
  a_s <- .hx_vals(ap); b_s <- .hx_vals(bq)
  types <- c(.hx_types(ap), .hx_types(bq))
  kw$maxterms <- .hx_kw(kw, "maxterms", .MPX$prec)
  r <- tryCatch(.hx_hypsum(p, q, types, c(a_s, b_s), z, kw), NoConvergence = function(err) NULL)
  if (!is.null(r)) return(r)
  prec <- .MPX$prec
  s <- tryCatch({
    tol <- .hx_kw(kw, "asymp_tol", .hx_div(.hx_eps(), 4))
    .MPX$prec <- .MPX$prec + 10
    t <- .hx_one(); s <- .hx_one(); out <- NULL
    k <- 1
    while (k < .MPX$prec) {
      for (a in a_s) t <- .hx_mul(t, .hx_add(a, k - 1))
      for (b in b_s) t <- .hx_div(t, .hx_add(b, k - 1))
      t <- .hx_mul(t, z)
      t <- .hx_div(t, k)
      s <- .hx_add(s, t)
      if (.hx_abs(t) <= tol) { out <- s; break }
      k <- k + 1
    }
    out
  }, finally = .MPX$prec <- prec)
  if (!is.null(s)) return(s)
  if (p <= q + 3) {
    contour <- kw$contour
    if (is.null(contour)) {
      if (.hx_arg(z) < 0.25) {
        u <- .hx_div(z, if (1 >= .hx_abs(z)) 1 else .hx_abs(z))
        j2 <- .mv_mpc(.mp_zero, .mpf(1, 1))
        if (.hx_arg(z) >= 0) contour <- list(0, j2, .hx_div(.hx_add(j2, 2), u), .hx_div(2, u), .mv_inf())
        else contour <- list(0, .hx_neg(j2), .hx_div(.hx_add(.hx_neg(j2), 2), u), .hx_div(2, u), .mv_inf())
      } else contour <- list(0, .mv_inf())
    }
    g <- function(t) .hx_mul(.hx_exp(.hx_neg(t)), .hx_hyper(a_s, c(b_s, list(1)), .hx_mul(t, z)))
    qr <- .mp_quad(g, contour, error = TRUE)
    if (.hx_abs(qr$err) <= .hx_mul(.hx_mul(.hx_abs(qr$value), .hx_eps()), 8)) return(qr$value)
  }
  .hx_noconv("")
}

.hx_hyp2f2 <- function(ap, bq, z, kw) {
  a1 <- ap[[1L]]$v; a2 <- ap[[2L]]$v; b1 <- bq[[1L]]$v; b2 <- bq[[2L]]$v
  absz <- .hx_abs(z)
  magz <- .hx_mag(z)
  orig <- .MPX$prec
  asymp_extraprec <- magz
  can_use_asymptotic <- !isTRUE(kw$force_series) && .hx_mag(absz) > 3
  if (can_use_asymptotic) {
    v <- tryCatch({
      .MPX$prec <- .MPX$prec + asymp_extraprec
      h <- function(a1, a2, b1, b2) {
        X <- .hx_sub(.hx_sub(.hx_add(a1, a2), b1), b2)
        A2 <- .hx_add(a1, a2)
        B2 <- .hx_add(b1, b2)
        cc <- list()
        cc[["0"]] <- .hx_one()
        cc[["1"]] <- .hx_sub(.hx_add(.hx_mul(.hx_sub(A2, 1), X), .hx_mul(b1, b2)), .hx_mul(a1, a2))
        s1 <- 0; k <- 0; tprev <- 0
        repeat {
          key <- as.character(k)
          if (is.null(cc[[key]])) {
            uu1 <- .hx_sub(1, B2)
            uu1 <- .hx_add(uu1, .hx_mul(2, a1))
            uu1 <- .hx_add(uu1, .hx_pow(a1, 2))
            uu1 <- .hx_add(uu1, .hx_mul(2, a2))
            uu1 <- .hx_add(uu1, .hx_pow(a2, 2))
            uu1 <- .hx_sub(uu1, .hx_mul(A2, B2))
            uu1 <- .hx_add(uu1, .hx_mul(a1, a2))
            uu1 <- .hx_add(uu1, .hx_mul(b1, b2))
            uu1 <- .hx_add(uu1, .hx_mul(.hx_sub(.hx_mul(2, B2), .hx_mul(3, .hx_add(A2, 1))), k))
            uu1 <- .hx_add(uu1, 2 * k^2)
            uu2 <- .hx_mul(.hx_mul(.hx_sub(.hx_add(.hx_sub(k, A2), b1), 1), .hx_sub(.hx_add(.hx_sub(k, A2), b2), 1)),
                           .hx_sub(.hx_sub(k, X), 2))
            cc[[key]] <- .hx_mul(.hx_div(.hx_one(), k),
                                 .hx_sub(.hx_mul(uu1, cc[[as.character(k - 1)]]), .hx_mul(uu2, cc[[as.character(k - 2)]])))
          }
          t1 <- .hx_mul(cc[[key]], .hx_pow(z, -k))
          if (.hx_abs(t1) < .hx_mul(0.1, .hx_eps())) break
          if (k > 5 && .hx_div(.hx_abs(tprev), .hx_abs(t1)) < 1.5) .hx_noconv("")
          s1 <- .hx_add(s1, t1)
          tprev <- t1
          k <- k + 1
        }
        S <- .hx_mul(.hx_exp(z), s1)
        list(.hx_term(list(z, S), list(X, 1), list(b1, b2), list(a1, a2), list(), list(), 0),
             .hx_term(list(.hx_neg(z)), list(.hx_neg(a1)), list(b1, b2, .hx_sub(a2, a1)), list(a2, .hx_sub(b1, a1), .hx_sub(b2, a1)),
                      list(a1, .hx_add(.hx_sub(a1, b1), 1), .hx_add(.hx_sub(a1, b2), 1)), list(.hx_add(.hx_sub(a1, a2), 1)), .hx_div(-1, z)),
             .hx_term(list(.hx_neg(z)), list(.hx_neg(a2)), list(b1, b2, .hx_sub(a1, a2)), list(a1, .hx_sub(b1, a2), .hx_sub(b2, a2)),
                      list(a2, .hx_add(.hx_sub(a2, b1), 1), .hx_add(.hx_sub(a2, b2), 1)), list(.hx_add(.hx_add(.hx_neg(a1), a2), 1)), .hx_div(-1, z)))
      }
      v <- .hx_hypercomb(h, list(a1, a2, b1, b2), kw = list(force_series = TRUE, maxterms = 4 * .MPX$prec))
      if (all(vapply(list(a1, a2, b1, b2, z), .hx_is_real_type, NA))) v <- .hx_re(v)
      v
    }, NoConvergence = function(err) NULL, finally = .MPX$prec <- orig)
    if (!is.null(v)) return(v)
  }
  .hx_hypsum(2, 2, c(ap[[1L]]$t, ap[[2L]]$t, bq[[1L]]$t, bq[[2L]]$t), list(a1, a2, b1, b2), z, kw)
}

# the oscillating sum shared by the asymptotic 1F2 and 2F3 expansions
.hx_osc_sum <- function(cc_next, X, z, power) {
  s1 <- 0; s2 <- 0; k <- 0; tprev <- 0
  repeat {
    ck <- cc_next(k)
    w <- .hx_mul(ck, power(.hx_neg(z), -0.5 * k))
    t1 <- .hx_mul(.hx_mul(.hx_pow(.hx_neg(.hx_j()), k), .hx_pow(.hx_convert(2), -k)), w)
    t2 <- .hx_mul(.hx_mul(.hx_pow(.hx_j(), k), .hx_pow(.hx_convert(2), -k)), w)
    if (.hx_abs(t1) < .hx_mul(0.1, .hx_eps())) break
    if (k > 5 && .hx_div(.hx_abs(tprev), .hx_abs(t1)) < 1.5) .hx_noconv("")
    s1 <- .hx_add(s1, t1)
    s2 <- .hx_add(s2, t2)
    tprev <- t1
    k <- k + 1
  }
  .hx_add(.hx_mul(.hx_expj(.hx_add(.hx_mul(.hx_pi(), X), .hx_mul(2, .hx_sqrt(.hx_neg(z))))), s1),
          .hx_mul(.hx_expj(.hx_neg(.hx_add(.hx_mul(.hx_pi(), X), .hx_mul(2, .hx_sqrt(.hx_neg(z)))))), s2))
}

.hx_hyp1f2 <- function(ap, bq, z, kw) {
  a1 <- ap[[1L]]$v; b1 <- bq[[1L]]$v; b2 <- bq[[2L]]$v
  absz <- .hx_abs(z)
  magz <- .hx_mag(z)
  orig <- .MPX$prec
  asymp_extraprec <- if (.hx_iszero(z)) 0 else magz %/% 2
  can_use_asymptotic <- !isTRUE(kw$force_series) && .hx_mag(absz) > 19 && .hx_sqrt(absz) > 1.5 * orig
  if (can_use_asymptotic) {
    v <- tryCatch({
      .MPX$prec <- .MPX$prec + asymp_extraprec
      h <- function(a1, b1, b2) {
        Q <- function(n, d) .mv_q(n, d)
        X <- .hx_mul(Q(1, 2), .hx_add(.hx_sub(.hx_sub(a1, b1), b2), Q(1, 2)))
        cc <- list()
        cc[["0"]] <- .hx_one()
        cc[["1"]] <- .hx_mul(2, .hx_sub(.hx_add(.hx_mul(.hx_mul(Q(1, 4), .hx_sub(.hx_add(.hx_add(.hx_mul(3, a1), b1), b2), 2)), .hx_sub(.hx_sub(a1, b1), b2)), .hx_mul(b1, b2)), Q(3, 16)))
        base <- .hx_sub(.hx_add(.hx_mul(b1, b2), .hx_mul(.hx_mul(Q(1, 4), .hx_sub(.hx_sub(a1, b1), b2)), .hx_sub(.hx_add(.hx_add(.hx_mul(3, a1), b1), b2), 2))), Q(3, 16))
        inner <- .hx_mul(.hx_mul(.hx_mul(-16, .hx_sub(.hx_mul(2, a1), 3)), b1), b2)
        inner <- .hx_add(inner, .hx_mul(.hx_mul(4, .hx_sub(.hx_sub(a1, b1), b2)),
                                        .hx_sub(.hx_add(.hx_add(.hx_add(.hx_mul(-8, .hx_pow(a1, 2)), .hx_mul(11, a1)), b1), b2), 2)))
        inner <- .hx_sub(inner, 3)
        cc[["2"]] <- .hx_add(.hx_mul(2, .hx_pow(base, 2)), .hx_mul(Q(1, 16), inner))
        cc_next <- function(k) {
          key <- as.character(k)
          if (!is.null(cc[[key]])) return(cc[[key]])
          uu1 <- .hx_add(.hx_mul(3, k^2), .hx_mul(.hx_sub(.hx_add(.hx_add(.hx_mul(-6, a1), .hx_mul(2, b1)), .hx_mul(2, b2)), 4), k))
          uu1 <- .hx_add(uu1, .hx_mul(3, .hx_pow(a1, 2)))
          uu1 <- .hx_sub(uu1, .hx_pow(.hx_sub(b1, b2), 2))
          uu1 <- .hx_sub(uu1, .hx_mul(.hx_mul(2, a1), .hx_sub(.hx_add(b1, b2), 2)))
          uu1 <- .hx_add(uu1, Q(1, 4))
          uu2 <- .hx_mul(.hx_mul(.hx_sub(.hx_sub(.hx_add(.hx_sub(k, a1), b1), b2), Q(1, 2)),
                                 .hx_sub(.hx_add(.hx_sub(.hx_sub(k, a1), b1), b2), Q(1, 2))),
                         .hx_sub(.hx_add(.hx_add(.hx_sub(k, a1), b1), b2), Q(5, 2)))
          v <- .hx_mul(.hx_div(.hx_one(), 2 * k), .hx_sub(.hx_mul(uu1, cc[[as.character(k - 1)]]), .hx_mul(uu2, cc[[as.character(k - 2)]])))
          cc[[key]] <<- v
          v
        }
        S <- .hx_osc_sum(cc_next, X, z, .hx_pow)
        list(.hx_term(list(.hx_mul(0.5, S), .hx_pi(), .hx_neg(z)), list(1, -0.5, X), list(b1, b2), list(a1), list(), list(), 0),
             .hx_term(list(.hx_neg(z)), list(.hx_neg(a1)), list(b1, b2), list(.hx_sub(b1, a1), .hx_sub(b2, a1)),
                      list(a1, .hx_add(.hx_sub(a1, b1), 1), .hx_add(.hx_sub(a1, b2), 1)), list(), .hx_div(1, z)))
      }
      v <- .hx_hypercomb(h, list(a1, b1, b2), kw = list(force_series = TRUE, maxterms = 4 * .MPX$prec))
      if (all(vapply(list(a1, b1, b2, z), .hx_is_real_type, NA))) v <- .hx_re(v)
      v
    }, NoConvergence = function(err) NULL, finally = .MPX$prec <- orig)
    if (!is.null(v)) return(v)
  }
  .hx_hypsum(1, 2, c(ap[[1L]]$t, bq[[1L]]$t, bq[[2L]]$t), list(a1, b1, b2), z, kw)
}

.hx_hyp2f3 <- function(ap, bq, z, kw) {
  a1 <- ap[[1L]]$v; a2 <- ap[[2L]]$v; b1 <- bq[[1L]]$v; b2 <- bq[[2L]]$v; b3 <- bq[[3L]]$v
  absz <- .hx_abs(z)
  magz <- .hx_mag(z)
  asymp_extraprec <- if (.hx_iszero(z)) 0 else magz %/% 2
  orig <- .MPX$prec
  can_use_asymptotic <- !isTRUE(kw$force_series) && .hx_mag(absz) > 19 && .hx_sqrt(absz) > 1.5 * orig
  if (can_use_asymptotic) {
    v <- tryCatch({
      .MPX$prec <- .MPX$prec + asymp_extraprec
      h <- function(a1, a2, b1, b2, b3) {
        Q <- function(n, d) .mv_q(n, d)
        X <- .hx_mul(Q(1, 2), .hx_add(.hx_sub(.hx_sub(.hx_sub(.hx_add(a1, a2), b1), b2), b3), Q(1, 2)))
        A2 <- .hx_add(a1, a2)
        B3 <- .hx_add(.hx_add(b1, b2), b3)
        A <- .hx_mul(a1, a2)
        B <- .hx_add(.hx_add(.hx_mul(b1, b2), .hx_mul(b3, b2)), .hx_mul(b1, b3))
        R <- .hx_mul(.hx_mul(b1, b2), b3)
        cc <- list()
        cc[["0"]] <- .hx_one()
        cc[["1"]] <- .hx_mul(2, .hx_sub(.hx_add(.hx_sub(B, A), .hx_mul(.hx_mul(Q(1, 4), .hx_sub(.hx_add(.hx_mul(3, A2), B3), 2)), .hx_sub(A2, B3))), Q(3, 16)))
        inner <- .hx_mul(.hx_mul(-16, .hx_sub(.hx_mul(2, A2), 3)), .hx_sub(B, A))
        inner <- .hx_add(inner, .hx_mul(32, R))
        inner <- .hx_add(inner, .hx_mul(.hx_mul(4, .hx_sub(.hx_add(.hx_add(.hx_add(.hx_mul(-8, .hx_pow(A2, 2)), .hx_mul(11, A2)), .hx_mul(8, A)), B3), 2)), .hx_sub(A2, B3)))
        inner <- .hx_sub(inner, 3)
        cc[["2"]] <- .hx_add(.hx_mul(Q(1, 2), .hx_pow(cc[["1"]], 2)), .hx_mul(Q(1, 16), inner))
        cc_next <- function(k) {
          key <- as.character(k)
          if (!is.null(cc[[key]])) return(cc[[key]])
          X2 <- .hx_mul(2, X)
          uu1 <- .hx_mul(.hx_mul(.hx_mul(.hx_sub(.hx_sub(k, X2), 3), .hx_sub(.hx_sub(.hx_sub(k, X2), .hx_mul(2, b1)), 1)),
                                 .hx_sub(.hx_sub(.hx_sub(k, X2), .hx_mul(2, b2)), 1)), .hx_sub(.hx_sub(.hx_sub(k, X2), .hx_mul(2, b3)), 1))
          k1 <- k - 1
          uu2 <- .hx_sub(.hx_mul(4, k1^3), .hx_mul(.hx_mul(6, .hx_add(.hx_mul(4, X), B3)), k1^2))
          uu2 <- .hx_add(uu2, .hx_mul(.hx_mul(2, .hx_sub(.hx_add(.hx_add(.hx_add(.hx_mul(24, .hx_pow(X, 2)), .hx_mul(.hx_mul(12, B3), X)), .hx_mul(4, B)), B3), 1)), k1))
          uu2 <- .hx_sub(uu2, .hx_mul(32, .hx_pow(X, 3)))
          uu2 <- .hx_sub(uu2, .hx_mul(.hx_mul(24, B3), .hx_pow(X, 2)))
          uu2 <- .hx_sub(uu2, .hx_mul(4, B))
          uu2 <- .hx_sub(uu2, .hx_mul(8, R))
          uu2 <- .hx_sub(uu2, .hx_mul(.hx_mul(4, .hx_sub(.hx_add(.hx_mul(4, B), B3), 1)), X))
          uu2 <- .hx_add(uu2, .hx_mul(2, B3))
          uu2 <- .hx_sub(uu2, 1)
          uu3 <- .hx_add(.hx_add(5 * k1^2, .hx_mul(.hx_mul(2, .hx_add(.hx_sub(.hx_add(.hx_mul(-10, X), A2), .hx_mul(3, B3)), 3)), k1)), .hx_mul(2, cc[["1"]]))
          v <- .hx_mul(.hx_div(.hx_one(), 2 * k),
                       .hx_add(.hx_sub(.hx_mul(uu1, cc[[as.character(k - 3)]]), .hx_mul(uu2, cc[[as.character(k - 2)]])), .hx_mul(uu3, cc[[as.character(k - 1)]])))
          cc[[key]] <<- v
          v
        }
        S <- .hx_osc_sum(cc_next, X, z, .hx_power)
        rz <- .hx_div(1, z)
        list(.hx_term(list(.hx_mul(0.5, S), .hx_pi(), .hx_neg(z)), list(1, -0.5, X), list(b1, b2, b3), list(a1, a2), list(), list(), 0),
             .hx_term(list(.hx_neg(z)), list(.hx_neg(a1)), list(b1, b2, b3, .hx_sub(a2, a1)), list(a2, .hx_sub(b1, a1), .hx_sub(b2, a1), .hx_sub(b3, a1)),
                      list(a1, .hx_add(.hx_sub(a1, b1), 1), .hx_add(.hx_sub(a1, b2), 1), .hx_add(.hx_sub(a1, b3), 1)), list(.hx_add(.hx_sub(a1, a2), 1)), rz),
             .hx_term(list(.hx_neg(z)), list(.hx_neg(a2)), list(b1, b2, b3, .hx_sub(a1, a2)), list(a1, .hx_sub(b1, a2), .hx_sub(b2, a2), .hx_sub(b3, a2)),
                      list(a2, .hx_add(.hx_sub(a2, b1), 1), .hx_add(.hx_sub(a2, b2), 1), .hx_add(.hx_sub(a2, b3), 1)), list(.hx_add(.hx_add(.hx_neg(a1), a2), 1)), rz))
      }
      v <- .hx_hypercomb(h, list(a1, a2, b1, b2, b3), kw = list(force_series = TRUE, maxterms = 4 * .MPX$prec))
      if (all(vapply(list(a1, a2, b1, b2, b3, z), .hx_is_real_type, NA))) v <- .hx_re(v)
      v
    }, NoConvergence = function(err) NULL, finally = .MPX$prec <- orig)
    if (!is.null(v)) return(v)
  }
  .hx_hypsum(2, 3, c(ap[[1L]]$t, ap[[2L]]$t, bq[[1L]]$t, bq[[2L]]$t, bq[[3L]]$t), list(a1, a2, b1, b2, b3), z, kw)
}

.hx_hyp2f0 <- function(ap, bq, z, kw) {
  a <- ap[[1L]]$v; b <- ap[[2L]]$v
  kwb <- kw
  kwb$maxterms <- .hx_kw(kwb, "maxterms", .MPX$prec)
  r <- tryCatch(.hx_hypsum(2, 0, c(ap[[1L]]$t, ap[[2L]]$t), list(a, b), z, kwb),
                NoConvergence = function(err) { if (isTRUE(kw$force_series)) stop(err); NULL })
  if (!is.null(r)) return(r)
  h <- function(a, b) {
    w <- .hx_sinpi(b)
    rz <- .hx_div(-1, z)
    list(.hx_term(list(.hx_pi(), w, rz), list(1, -1, a), list(), list(.hx_add(.hx_sub(a, b), 1), b), list(a), list(b), rz),
         .hx_term(list(.hx_neg(.hx_pi()), w, rz), list(1, -1, .hx_sub(.hx_add(1, a), b)), list(), list(a, .hx_sub(2, b)),
                  list(.hx_add(.hx_sub(a, b), 1)), list(.hx_sub(2, b)), rz))
  }
  .hx_hypercomb(h, list(a, .hx_sub(.hx_add(1, a), b)), kw = kw)
}

# meijerg(a_s, b_s, z, r=1, series=None)
.hx_meijerg <- function(an, ap, bm, bq, z, r = 1, series = NULL, kw = list()) {
  n <- length(an); p <- n + length(ap); m <- length(bm); q <- m + length(bq)
  a <- lapply(c(an, ap), .mv_convert)
  b <- lapply(c(bm, bq), .mv_convert)
  z <- .mv_convert(z)
  r <- .mv(r)
  if (is.null(series)) {
    if (p < q) series <- 1
    if (p > q) series <- 2
    if (p == q) series <- if (m + n == p && .hx_abs(z) > 1) 2 else 1
  }
  if (series == 1) {
    h <- function(...) {
      args <- list(...)
      a <- args[seq_len(p)]; b <- args[p + seq_len(q)]
      lapply(seq_len(m), function(k) {
        bk <- b[[k]]
        gn <- c(lapply(setdiff(seq_len(m), k), function(j) .hx_sub(b[[j]], bk)),
                lapply(seq_len(n), function(j) .hx_add(.hx_sub(1, a[[j]]), bk)))
        gd <- c(lapply(n + seq_len(p - n), function(j) .hx_sub(a[[j]], bk)),
                lapply(m + seq_len(q - m), function(j) .hx_add(.hx_sub(1, b[[j]]), bk)))
        hn <- lapply(seq_len(p), function(j) .hx_add(.hx_sub(1, a[[j]]), bk))
        hd <- lapply(setdiff(seq_len(q), k), function(j) .hx_add(.hx_sub(1, b[[j]]), bk))
        hz <- .hx_mul(.hx_pow(.hx_neg(.hx_one()), p - m - n), .hx_pow(z, .hx_div(.hx_one(), r)))
        .hx_term(list(z), list(.hx_div(bk, r)), gn, gd, hn, hd, hz)
      })
    }
  } else {
    h <- function(...) {
      args <- list(...)
      a <- args[seq_len(p)]; b <- args[p + seq_len(q)]
      lapply(seq_len(n), function(k) {
        ak <- a[[k]]
        ex <- if (r == 1) .hx_sub(ak, 1) else .hx_div(.hx_sub(ak, 1), .mv_convert(r))
        gn <- c(lapply(setdiff(seq_len(n), k), function(j) .hx_sub(ak, a[[j]])),
                lapply(seq_len(m), function(j) .hx_add(.hx_sub(1, ak), b[[j]])))
        gd <- c(lapply(m + seq_len(q - m), function(j) .hx_sub(ak, b[[j]])),
                lapply(n + seq_len(p - n), function(j) .hx_add(.hx_sub(1, ak), a[[j]])))
        hn <- lapply(seq_len(q), function(j) .hx_add(.hx_sub(1, ak), b[[j]]))
        hd <- lapply(setdiff(seq_len(p), k), function(j) .hx_sub(.hx_add(1, a[[j]]), ak))
        hz <- .hx_div(.hx_pow(.hx_neg(.hx_one()), q - m - n), .hx_pow(z, .hx_div(.hx_one(), r)))
        .hx_term(list(z), list(ex), gn, gd, hn, hd, hz)
      })
    }
  }
  .hx_hypercomb(h, c(a, b), kw = kw)
}
