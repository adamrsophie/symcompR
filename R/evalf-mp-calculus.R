# mpmath's series acceleration, Euler-Maclaurin summation and tanh-sinh
# quadrature (calculus/extrapolation.py, calculus/differentiation.py,
# calculus/quadrature.py) on the values of evalf-mp-ctx.R, with the helpers
# they call: digamma and polygamma (libmp/gammazeta.py), log gamma, rising
# factorials and factorials.

# prec_to_dps
.mp_dps <- function(p) max(1L, as.integer(round(p / 3.3219280948873626) - 1))

# ctx.factorial: gamma(x + 1) with x + 1 exact
.hx_factorial <- function(a) {
  a <- .mv_convert(a)
  if (a$k == "mpf") {
    x1 <- .mp_add(a$re, .mpf(1, 0), max(.bitlen(a$re$m), 1) + abs(a$re$e) + 4)
    if (.mp_is_int(x1) && x1$m <= 0) .hx_valueerror("gamma function pole")
    return(.mv_mpf(.mp_gamma_r(x1, .MPX$prec)))
  }
  .hx_gamma(.hx_add(a, 1))
}

# rf(x, n) = gammaprod([x + n], [x]) with x + n at twice the precision
.mp_rf <- function(x, n) {
  p <- .MPX$prec
  .MPX$prec <- 2 * p
  xn <- .hx_add(x, n)
  .MPX$prec <- p
  .hx_gammaprod(list(xn), list(x))
}

# ctx.sign
.hx_csign <- function(x) {
  x <- .mv_convert(x)
  if (.hx_iszero(x) || .hx_isnan(x)) return(x)
  if (x$k != "mpc") return(if (x > 0) .hx_one() else .hx_neg(.hx_one()))
  .hx_div(x, .hx_abs(x))
}

# ctx.loggamma: real log gamma for x > 0, else the principal branch from the
# complex Stirling series
.hx_loggamma <- function(x) {
  x <- .mv_convert(x)
  p <- .MPX$prec
  if (x$k == "mpf" && x$re$m > 0) return(.mv_mpf(.mp_loggamma(x$re, p)))
  if (.hx_isnpint(x)) return(.mv_inf())
  z <- .mv_as_mpc(x)
  mg <- max(.mp_mag(z$re), .mp_mag(z$im))
  wp <- p + 20 + max(0, mg)
  shift <- 0
  zr <- .mp_to_dbl(z$re); zi <- .mp_to_dbl(z$im)
  z0 <- max(10, wp / 4)
  if (zr * zr + zi * zi < z0 * z0 || zr < 0) shift <- max(0, ceiling(z0 - zr))
  w <- .mpc_add(.mpc(z$re, z$im), .mpc(.mpf(shift, 0)), wp)
  s <- .mp_stirling_terms(w, wp, TRUE)
  if (shift > 0) {
    zz <- .mpc(z$re, z$im)
    for (j in 0:(shift - 1)) s <- .mpc_sub(s, .mpc_log(.mpc_add(zz, .mpc(.mpf(j, 0)), wp), wp), wp)
  }
  .mv_mpc(.mp_round(s$re$m, s$re$e, p), .mp_round(s$im$m, s$im$e, p))
}

# ctx.psi(m, z): mpf_psi0, mpc_psi0 and mpc_psi
.mp_psi <- function(m, z) {
  z <- .mv_convert(z)
  m <- as.integer(m)
  p <- .MPX$prec
  if (z$k == "mpf") {
    if (m == 0) return(.mv_mpf(.zt_psi0(z$re, p)))
    r <- .mpc_psi(m, .mpc(z$re), p)
    return(.mv_mpf(r$re))
  }
  r <- if (m == 0) .mpc_psi0(.mpc(z$re, z$im), p) else .mpc_psi(m, .mpc(z$re, z$im), p)
  .mv_mpc(r$re, r$im)
}
.mp_to_int <- function(a) if (a$e >= 0) a$m * .pow2(a$e) else { q <- gmp::divq.bigz(abs(a$m), .pow2(-a$e)); if (a$m < 0) -q else q }
.mpc_psi0 <- function(z, prec) {
  if (z$im$m == 0) return(.mpc(.zt_psi0(z$re, prec)))
  wp <- prec + 20
  re <- z$re
  if (re$m < 0 && .mp_mag(re) > 3) {
    cs <- .mpc_cos_sin_pi(z, wp)
    q <- .mpc_mul(.mpc_div(cs$cos, cs$sin, wp), .mpc(.mp_pi(wp)), wp)
    pp <- .mpc_psi0(.mpc_sub(.mpc_one, z, wp), wp)
    return(.mpc_sub(pp, q, prec))
  }
  if (re$m > 0 && .mp_mag(re) > wp) return(.mpc_log(.mpc_sub(z, .mpc_one, wp), prec))
  w <- if (re$m == 0) gmp::as.bigz(0) else .mp_to_int(re)
  n <- as.integer(0.11 * wp) + 2
  s <- .mpc(.mp_zero)
  if (w < n) {
    k <- w
    while (k < n) {
      s <- .mpc_sub(s, .mpc_div(.mpc_one, z, wp), wp)
      z <- .mpc(.mp_add(z$re, .mpf(1, 0), wp), z$im)
      k <- k + 1
    }
  }
  z <- .mpc_sub(z, .mpc_one, wp)
  s <- .mpc_add(s, .mpc_log(z, wp), wp)
  s <- .mpc_add(s, .mpc_div(.mpc(.mpf(1, -1)), z, wp), wp)
  z2 <- .mpc_mul(z, z, wp)
  t <- .mpc_one
  szprev <- .mp_zero
  eps <- .mpf(1, -wp + 2)
  k <- 1
  repeat {
    t <- .mpc_mul(t, z2, wp)
    bern <- .mp_from_q(.mp_bernoulli(2 * k), wp)
    t2 <- .mpc(.mp_mul(t$re, .mpf(2 * k, 0), wp), .mp_mul(t$im, .mpf(2 * k, 0), wp))
    m2 <- .mp_add(.mp_xmul(t2$re, t2$re), .mp_xmul(t2$im, t2$im), wp + 10)
    term <- .mpc(.mp_div(.mp_xmul(t2$re, bern), m2, wp), .mp_div(.mp_neg(.mp_xmul(t2$im, bern)), m2, wp))
    s <- .mpc_sub(s, term, wp)
    szterm <- .mpc_abs(term, 10)
    if (k > 2 && (.mp_cmp(szterm, eps) <= 0 || .mp_cmp(szprev, szterm) <= 0)) break
    szprev <- szterm
    k <- k + 1
  }
  s
}
.mpc_psi <- function(m, z, prec) {
  if (m == 0) return(.mpc_psi0(z, prec))
  wp <- prec + 20
  w <- if (z$re$m == 0) gmp::as.bigz(0) else .mp_to_int(z$re)
  n <- as.integer(0.4 * wp + 4 * m)
  s <- .mpc(.mp_zero)
  if (w < n) {
    k <- w
    while (k < n) {
      s <- .mpc_add(s, .mpc_pow_int(z, -m - 1, wp), wp)
      z <- .mpc(.mp_add(z$re, .mpf(1, 0), wp), z$im)
      k <- k + 1
    }
  }
  zm <- .mpc_pow_int(z, -m, wp)
  z2 <- .mpc_pow_int(z, -2, wp)
  s <- .mpc_add(s, .mpc(.mp_div(zm$re, .mpf(m, 0), wp), .mp_div(zm$im, .mpf(m, 0), wp)), wp)
  zz <- .mpc_div(zm, z, wp)
  s <- .mpc_add(s, .mpc(.mp_mul(zz$re, .mpf(1, -1), wp), .mp_mul(zz$im, .mpf(1, -1), wp)), wp)
  a <- gmp::as.bigz(m + 1); b <- gmp::as.bigz(2)
  k <- 1
  magn <- .mpc_abs(s, 10)
  eps <- .mpf(1, .mp_mag(magn) - wp + 2)
  repeat {
    zm <- .mpc_mul(zm, z2, wp)
    bern <- .mp_from_q(.mp_bernoulli(2 * k), wp)
    scal <- .mp_round(bern$m * a, bern$e, wp)
    scal <- .mp_div(scal, .mpf(b, 0), wp)
    term <- .mpc(.mp_mul(zm$re, scal, wp), .mp_mul(zm$im, scal, wp))
    s <- .mpc_add(s, term, wp)
    szterm <- .mpc_abs(term, 10)
    if (k > 2 && .mp_cmp(szterm, eps) <= 0) break
    a <- a * (m + 2 * k) * (m + 2 * k + 1)
    b <- b * (2 * k + 1) * (2 * k + 2)
    k <- k + 1
  }
  g <- .mp_gamma_r(.mpf(m + 1, 0), wp)
  v <- .mpc(.mp_mul(s$re, g, prec), .mp_mul(s$im, g, prec))
  if (m %% 2 == 0) v <- .mpc(.mp_neg(v$re), .mp_neg(v$im))
  v
}

# ---------------------------------------------------------------------------
# Extrapolation (calculus/extrapolation.py)

.mp_richardson <- function(seq) {
  if (length(seq) < 3) .hx_valueerror("seq should be of minimum length 3")
  L <- length(seq)
  if (!.mv_eq(.hx_csign(.hx_sub(seq[[L]], seq[[L - 1]])), .hx_csign(.hx_sub(seq[[L - 1]], seq[[L - 2]]))))
    seq <- seq[seq(1, L, by = 2)]
  N <- length(seq) %/% 2 - 1
  s <- .hx_zero()
  c <- .hx_div(.mv_int((-1)^N * gmp::pow.bigz(N, N)), .mv_convert(.mv_int(gmp::factorialZ(N))))
  maxc <- .mv_int(1)
  for (k in 0:N) {
    s <- .hx_add(s, .hx_mul(c, seq[[N + k + 1]]))
    ac <- .hx_abs(c)
    maxc <- if (maxc > ac) maxc else ac
    c <- .hx_mul(c, .hx_mul(k - N, .hx_pow(.mv_convert(k + N + 1), N)))
    c <- .hx_div(c, .hx_mul(1 + k, .hx_pow(.mv_convert(k + N), N)))
  }
  list(value = s, maxc = maxc)
}

.mp_shanks <- function(seq, table = list()) {
  if (length(seq) < 2) .hx_valueerror("seq should be of minimum length 2")
  START <- length(table)
  STOP <- length(seq) - 1
  if (STOP %% 2 == 1) STOP <- STOP - 1
  one <- .hx_one()
  eps <- .mv_pos(.hx_eps())
  i <- START
  while (i < STOP) {
    row <- list()
    for (j in 0:i) {
      if (j == 0) {
        a <- 0; b <- .hx_sub(seq[[i + 2]], seq[[i + 1]])
      } else {
        a <- if (j == 1) seq[[i + 1]] else table[[i]][[j - 1]]
        b <- .hx_sub(row[[j]], table[[i]][[j]])
      }
      # mpmath draws a random multiple of eps here (randomized=True)
      if (.hx_iszero(b)) b <- .hx_mul(1, eps)
      row[[j + 1L]] <- .hx_add(a, .hx_div(one, b))
    }
    table[[i + 1L]] <- row
    i <- i + 1
  }
  table
}

.mp_adaptive_extrapolation <- function(update, kw) {
  tol <- if (!is.null(kw$tol)) kw$tol else .hx_div(.hx_eps(), 2^10)
  maxterms <- if (!is.null(kw$maxterms)) kw$maxterms else .mp_dps(.MPX$prec) * 10
  method <- strsplit(if (is.null(kw$method)) "r+s" else kw$method, "+", fixed = TRUE)[[1]]
  strict <- isTRUE(kw$strict)
  if ("d" %in% method || "direct" %in% method) {
    TRY_RICHARDSON <- TRY_SHANKS <- FALSE
  } else {
    TRY_RICHARDSON <- "r" %in% method || "richardson" %in% method
    TRY_SHANKS <- "s" %in% method || "shanks" %in% method
  }
  last_richardson_value <- 0
  shanks_table <- list()
  index <- 0
  step <- 0
  partial <- list()
  best <- .hx_zero()
  orig <- .MPX$prec
  on.exit(.MPX$prec <- orig)
  .MPX$prec <- if (TRY_RICHARDSON || TRY_SHANKS) (.MPX$prec + 10) * 4 else .MPX$prec + 30
  repeat {
    if (index >= maxterms) break
    step <- step + 10
    partial <- update(partial, index:(index + step - 1))
    index <- index + step
    best <- partial[[length(partial)]]
    error <- .hx_abs(.hx_sub(best, partial[[length(partial) - 1]]))
    if (error <= tol) return(best)
    if (TRY_RICHARDSON) {
      rr <- .mp_richardson(partial)
      richardson_error <- .hx_abs(.hx_sub(rr$value, last_richardson_value))
      if (richardson_error <= tol) return(rr$value)
      last_richardson_value <- rr$value
      if (.hx_mul(.hx_eps(), rr$maxc) > tol) TRY_RICHARDSON <- FALSE
      if (richardson_error < error) { error <- richardson_error; best <- rr$value }
    }
    if (TRY_SHANKS) {
      shanks_table <- .mp_shanks(partial, shanks_table)
      row <- shanks_table[[length(shanks_table)]]
      if (length(row) == 2) {
        est1 <- row[[2]]; shanks_error <- 0; maxc <- NULL
      } else {
        est1 <- row[[length(row)]]; maxc <- .hx_abs(row[[length(row) - 1]]); est2 <- row[[length(row) - 2]]
        shanks_error <- .hx_abs(.hx_sub(est1, est2))
      }
      if (shanks_error <= tol) return(est1)
      if (!is.null(maxc) && .hx_mul(.hx_eps(), maxc) > tol) TRY_SHANKS <- FALSE
      if (shanks_error < error) { error <- shanks_error; best <- est1 }
    }
  }
  if (strict) .hx_noconv("")
  best
}

# nsum(f, [a, inf]) for one interval with finite a
.mp_nsum <- function(f, interval, ...) {
  kw <- list(...)
  a <- interval[[1L]]
  g <- function(k) f(.hx_add(k, a))
  update <- function(partial, indices) {
    psum <- if (length(partial)) partial[[length(partial)]] else .hx_zero()
    for (k in indices) {
      psum <- .hx_add(psum, g(.mv_convert(k)))
      partial[[length(partial) + 1L]] <- psum
    }
    partial
  }
  .mv_pos(.mp_adaptive_extrapolation(update, kw))
}

# ---------------------------------------------------------------------------
# Derivatives of exp(f) from those of f (diffs_exp with dpoly)

.MP_DPOLY <- new.env(parent = emptyenv())
.mp_dpoly <- function(n) {
  key <- as.character(n)
  if (!is.null(.MP_DPOLY[[key]])) return(.MP_DPOLY[[key]])
  if (n == 0) { r <- list(keys = list(0L), counts = list(gmp::as.bigz(1))); assign("0", r, envir = .MP_DPOLY); return(r) }
  R <- .mp_dpoly(n - 1)
  rk <- lapply(R$keys, function(c) c(c, 0L))
  rc <- R$counts
  ks <- list(); cs <- list(); idx <- new.env(parent = emptyenv())
  put <- function(k, v) {
    s <- paste(k, collapse = ",")
    i <- idx[[s]]
    if (is.null(i)) {
      ks[[length(ks) + 1L]] <<- k; cs[[length(cs) + 1L]] <<- v
      assign(s, length(ks), envir = idx)
    } else cs[[i]] <<- cs[[i]] + v
  }
  for (i in seq_along(rk)) {
    pw <- rk[[i]]
    pw1 <- pw; pw1[1] <- pw1[1] + 1L
    put(pw1, rc[[i]])
  }
  for (i in seq_along(rk)) {
    pw <- rk[[i]]
    if (!sum(pw)) next
    for (k in seq_along(pw)) {
      pk <- pw[k]
      if (pk) {
        pw2 <- pw; pw2[k] <- pk - 1L; pw2[k + 1] <- pw[k + 1] + 1L
        put(pw2, pk * rc[[i]])
      }
    }
  }
  r <- list(keys = ks, counts = cs)
  assign(key, r, envir = .MP_DPOLY)
  r
}
# fdiffs is a generator function returning the next derivative on each call
.mp_diffs_exp <- function(fdiffs) {
  cache <- list()
  fn <- function(k) {
    while (length(cache) <= k) cache[[length(cache) + 1L]] <<- fdiffs()
    cache[[k + 1L]]
  }
  f0 <- NULL
  i <- -1
  function() {
    i <<- i + 1
    if (i == 0) { f0 <<- .hx_exp(fn(0)); return(f0) }
    s <- .hx_zero()
    dp <- .mp_dpoly(i)
    for (j in seq_along(dp$keys)) {
      pw <- dp$keys[[j]]
      fs <- list()
      for (k in seq_along(pw)) if (pw[k]) fs[[length(fs) + 1L]] <- .hx_pow(fn(k), pw[k])
      s <- .hx_add(s, .hx_mul(.mv_int(dp$counts[[j]]), .hx_fprod(fs)))
    }
    .hx_mul(s, f0)
  }
}

# sumem(f, [a, inf], tol, adiffs=..., error=True, _fast_abort=...)
.mp_sumem <- function(f, interval, tol = NULL, reject = 10, adiffs, error = FALSE, fast_abort = FALSE) {
  if (is.null(tol)) tol <- .mv_pos(.hx_eps())
  a <- .mv_convert(interval[[1L]])
  err <- .hx_zero()
  prev <- 0
  orig <- .MPX$prec
  on.exit(.MPX$prec <- orig)
  .MPX$prec <- .MPX$prec + 10
  s <- .hx_zero()
  for (k in 0:9999) {
    da <- adiffs()
    if (k %% 2 == 1) {
      term <- .hx_div(.hx_mul(.hx_sub(0, da), .mv_mpf(.mp_from_q(.mp_bernoulli(k + 1), .MPX$prec))), .hx_factorial(k + 1))
      mag <- .hx_abs(term)
      if (k > 4 && mag < tol) {
        s <- .hx_add(s, term); break
      } else if (k > 4 && .hx_div(.hx_abs(prev), mag) < reject) {
        err <- .hx_add(err, mag)
        if (fast_abort) return(list(value = s, err = err))
        break
      } else s <- .hx_add(s, term)
      prev <- term
    }
  }
  s <- .hx_add(s, .hx_div(f(a), 2))
  qr <- .mp_quad(f, interval, error = TRUE)
  s <- .hx_add(s, qr$value)
  err <- .hx_add(err, qr$err)
  list(value = s, err = err)
}

# ---------------------------------------------------------------------------
# Tanh-sinh quadrature (QuadratureRule and TanhSinh)

.MP_TS_NODES <- new.env(parent = emptyenv())
.mp_ts_calc_nodes <- function(degree, prec) {
  key <- paste(degree, prec)
  if (!is.null(.MP_TS_NODES[[key]])) return(.MP_TS_NODES[[key]])
  old <- .MPX$prec
  on.exit(.MPX$prec <- old)
  .MPX$prec <- .MPX$prec + 20
  nodes <- list()
  tol <- .hx_ldexp(1, -prec - 10)
  pi4 <- .hx_div(.hx_pi(), 4)
  t0 <- .hx_ldexp(1, -degree)
  if (degree == 1) {
    nodes[[1L]] <- list(x = .hx_zero(), w = .hx_div(.hx_pi(), 2))
    h <- t0
  } else h <- .hx_mul(t0, 2)
  expt0 <- .hx_exp(t0)
  a <- .hx_mul(pi4, expt0)
  b <- .hx_div(pi4, expt0)
  udelta <- .hx_exp(h)
  urdelta <- .hx_div(1, udelta)
  for (k in 0:(20 * 2^degree)) {
    c <- .hx_exp(.hx_sub(a, b))
    d <- .hx_div(1, c)
    co <- .hx_div(.hx_add(c, d), 2)
    si <- .hx_div(.hx_sub(c, d), 2)
    x <- .hx_div(si, co)
    w <- .hx_div(.hx_add(a, b), .hx_pow(co, 2))
    diff <- .hx_abs(.hx_sub(x, 1))
    if (diff <= tol) break
    nodes[[length(nodes) + 1L]] <- list(x = x, w = w)
    nodes[[length(nodes) + 1L]] <- list(x = .hx_neg(x), w = w)
    a <- .hx_mul(a, udelta)
    b <- .hx_mul(b, urdelta)
  }
  assign(key, nodes, envir = .MP_TS_NODES)
  nodes
}
.mp_ts_transform <- function(nodes, a, b) {
  a <- .mv_convert(a); b <- .mv_convert(b)
  one <- .hx_one()
  if (.mv_eq(a, -1) && .mv_eq(b, 1)) return(nodes)
  half <- .mv_convert(0.5)
  if (.hx_isinf(a) || .hx_isinf(b)) {
    if (identical(a$sp, "-inf") && identical(b$sp, "inf")) {
      return(lapply(nodes, function(nd) {
        x2 <- .hx_mul(nd$x, nd$x); px1 <- .hx_sub(one, x2); spx1 <- .hx_pow(px1, .hx_neg(half))
        list(x = .hx_mul(nd$x, spx1), w = .hx_mul(nd$w, .hx_div(spx1, px1)))
      }))
    }
    if (identical(a$sp, "-inf")) {
      b1 <- .hx_add(b, 1)
      return(lapply(nodes, function(nd) {
        u <- .hx_div(2, .hx_add(nd$x, one))
        list(x = .hx_sub(b1, u), w = .hx_mul(nd$w, .hx_mul(half, .hx_pow(u, 2))))
      }))
    }
    if (identical(b$sp, "inf")) {
      a1 <- .hx_sub(a, 1)
      return(lapply(nodes, function(nd) {
        u <- .hx_div(2, .hx_add(nd$x, one))
        list(x = .hx_add(a1, u), w = .hx_mul(nd$w, .hx_mul(half, .hx_pow(u, 2))))
      }))
    }
    return(lapply(.mp_ts_transform(nodes, b, a), function(nd) list(x = nd$x, w = .hx_neg(nd$w))))
  }
  C <- .hx_div(.hx_sub(b, a), 2)
  D <- .hx_div(.hx_add(b, a), 2)
  lapply(nodes, function(nd) list(x = .hx_add(D, .hx_mul(C, nd$x)), w = .hx_mul(C, nd$w)))
}
# fdot over (w, f(x)): exact products, one rounding per part
.mp_fdot <- function(pairs) {
  real <- list(); imag <- list(); cx <- FALSE
  for (pr in pairs) {
    a <- .mv_convert(pr[[1L]]); b <- .mv_convert(pr[[2L]])
    if (a$k == "mpf" && b$k == "mpf") { real[[length(real) + 1L]] <- .mp_xmul(a$re, b$re); next }
    cx <- TRUE
    if (a$k == "mpf") {
      real[[length(real) + 1L]] <- .mp_xmul(a$re, b$re); imag[[length(imag) + 1L]] <- .mp_xmul(a$re, b$im)
    } else if (b$k == "mpf") {
      real[[length(real) + 1L]] <- .mp_xmul(a$re, b$re); imag[[length(imag) + 1L]] <- .mp_xmul(a$im, b$re)
    } else {
      real[[length(real) + 1L]] <- .mp_xmul(a$re, b$re); real[[length(real) + 1L]] <- .mp_neg(.mp_xmul(a$im, b$im))
      imag[[length(imag) + 1L]] <- .mp_xmul(a$re, b$im); imag[[length(imag) + 1L]] <- .mp_xmul(a$im, b$re)
    }
  }
  s <- .mp_sum(real, .MPX$prec)
  if (cx) .mv_mpc(s, .mp_sum(imag, .MPX$prec)) else .mv_mpf(s)
}
.mp_ts_estimate_error <- function(results, prec, epsilon) {
  n <- length(results)
  if (n == 2) return(.hx_abs(.hx_sub(results[[1L]], results[[2L]])))
  if (.mv_eq(results[[n]], results[[n - 1]]) && .mv_eq(results[[n - 1]], results[[n - 2]])) return(.hx_zero())
  D1 <- .hx_div(.hx_log(.hx_abs(.hx_sub(results[[n]], results[[n - 1]]))), .hx_log(.mv_convert(10)))
  D2 <- .hx_div(.hx_log(.hx_abs(.hx_sub(results[[n]], results[[n - 2]]))), .hx_log(.mv_convert(10)))
  D3 <- -prec
  m <- .hx_div(.hx_pow(D1, 2), D2)
  mx <- m
  if (.hx_mul(2, D1) > mx) mx <- .hx_mul(2, D1)
  if (.mv(D3) > mx) mx <- .mv(D3)
  D4 <- if (mx < 0) mx else .mv_int(0)
  d4 <- .mv_dbl(D4)
  .hx_pow(.mv_convert(10), .mv_int(trunc(d4)))
}
.mp_quad <- function(f, points, error = FALSE) {
  prec <- .MPX$prec
  epsilon <- .hx_div(.hx_eps(), 8)
  m <- as.integer(4 + max(0, log2(prec / 30))) + 2
  old <- .MPX$prec
  on.exit(.MPX$prec <- old)
  .MPX$prec <- .MPX$prec + 20
  I <- .hx_zero(); total_err <- .hx_zero()
  for (i in seq_len(length(points) - 1L)) {
    a <- points[[i]]; b <- points[[i + 1L]]
    if (.mv_eq(.mv_convert(a), .mv_convert(b))) next
    results <- list()
    err <- .hx_zero()
    for (degree in seq_len(m)) {
      .MPX$prec <- prec + 20
      nodes <- .mp_ts_transform(.mp_ts_calc_nodes(degree, prec), a, b)
      .MPX$prec <- prec + 20
      h <- .hx_pow(.mv_convert(2), -degree)
      S <- if (length(results)) .hx_div(results[[length(results)]], .hx_mul(h, 2)) else .hx_zero()
      S <- .hx_add(S, .mp_fdot(lapply(nodes, function(nd) list(nd$w, f(nd$x)))))
      results[[length(results) + 1L]] <- .hx_mul(h, S)
      if (degree > 1) {
        err <- .mp_ts_estimate_error(results, prec, epsilon)
        if (err <= epsilon) break
      }
    }
    I <- .hx_add(I, results[[length(results)]])
    total_err <- .hx_add(total_err, err)
  }
  .MPX$prec <- old
  if (error) return(list(value = .mv_pos(I), err = total_err))
  .mv_pos(I)
}

# mpmath's quadsubdiv: quad on each interval of a queue, bisecting an
# interval while its error estimate is not below tol
.mp_quadsubdiv <- function(f, interval, tol = NULL, maxintervals = NULL) {
  queue <- list()
  for (i in seq_len(length(interval) - 1L)) queue[[length(queue) + 1L]] <- list(interval[[i]], interval[[i + 1L]])
  total <- .hx_zero(); total_error <- .hx_zero()
  if (is.null(maxintervals)) maxintervals <- 10 * .MPX$prec
  count <- 0
  if (is.null(tol)) tol <- .mv_pos(.hx_eps())
  orig <- .MPX$prec
  on.exit(.MPX$prec <- orig)
  .MPX$prec <- .MPX$prec + 5
  while (length(queue)) {
    ab <- queue[[length(queue)]]; queue[[length(queue)]] <- NULL
    a <- ab[[1L]]; b <- ab[[2L]]
    r <- .mp_quad(f, list(a, b), error = TRUE)
    if (.mv_cmp(r$err, tol) < 0 || count > maxintervals) {
      total <- .hx_add(total, r$value); total_error <- .hx_add(total_error, r$err)
    } else {
      count <- count + 1
      ainf <- identical(.mv(a)$sp, "-inf"); binf <- identical(.mv(b)$sp, "inf")
      m <- if (ainf && binf) .hx_zero()
           else if (ainf) { u <- .hx_sub(b, 1); v <- .hx_mul(2, b); if (.mv_cmp(v, u) < 0) v else u }
           else if (binf) { u <- .hx_add(a, 1); v <- .hx_mul(2, a); if (.mv_cmp(v, u) > 0) v else u }
           else .hx_add(a, .hx_div(.hx_sub(b, a), 2))
      queue[[length(queue) + 1L]] <- list(a, m)
      queue[[length(queue) + 1L]] <- list(m, b)
    }
  }
  .MPX$prec <- orig
  .mv_pos(total)
}
