# Multivariate factorization over the integers by Wang's EEZ algorithm,
# after dmp_zz_wang and its helpers in sympy/polys/factortools.py, with the
# dense recursive representation of sympy/polys/densebasic.py: a polynomial
# of level u > 0 is a list of level u - 1 polynomials (leading first), a
# level 0 polynomial is a bigz vector (leading first).

.dz <- function(x) gmp::as.bigz(x)
.dz0 <- function() gmp::as.bigz(integer(0))

# densebasic
.dm_zero <- function(u) if (u == 0L) .dz0() else list()
.dm_one <- function(u) if (u == 0L) .dz(1L) else list(.dm_one(u - 1L))
.dm_ground <- function(c, u) if (u == 0L) { if (c == 0) .dz0() else .dz(c) } else { if (c == 0) list() else list(.dm_ground(c, u - 1L)) }
.dm_zero_p <- function(f, u) length(f) == 0L
.dm_strip <- function(f, u) {
  if (u == 0L) return(.zp_trim(f))
  i <- 1L
  while (i <= length(f) && .dm_zero_p(f[[i]], u - 1L)) i <- i + 1L
  if (i > length(f)) list() else f[i:length(f)]
}
.dm_degree <- function(f, u) if (.dm_zero_p(f, u)) -Inf else length(f) - 1L
.dm_degree_in <- function(f, j, u) {
  if (j == 0L) return(.dm_degree(f, u))
  if (.dm_zero_p(f, u)) return(-Inf)
  max(vapply(f, function(c) .dm_degree_in(c, j - 1L, u - 1L), 0))
}
.dm_degree_list <- function(f, u) vapply(0:u, function(j) .dm_degree_in(f, j, u), 0)
.dm_LC <- function(f, u) if (.dm_zero_p(f, u)) .dm_zero(u - 1L) else f[[1L]]
.dm_ground_LC <- function(f, u) { while (u > 0L) { f <- .dm_LC(f, u); u <- u - 1L }; if (length(f)) f[1L] else .dz(0L) }
.dm_ground_coeffs <- function(f, u) if (u == 0L) f else do.call(c, c(list(.dz0()), lapply(f, .dm_ground_coeffs, u = u - 1L)))
.dm_nest <- function(f, l) { for (i in seq_len(l)) f <- list(f); f }
.dm_raise <- function(f, l, u) {
  if (l == 0L) return(f)
  if (u == 0L) {
    if (!length(f)) return(list())
    return(lapply(seq_along(f), function(i) .dm_ground(f[i], l - 1L)))
  }
  lapply(f, .dm_raise, l = l, u = u - 1L)
}

# densearith
.dm_pad <- function(f, n, u) if (length(f) >= n) f else c(rep(list(.dm_zero(u - 1L)), n - length(f)), f)
.dm_add <- function(f, g, u) {
  if (u == 0L) return(.zp_add(f, g))
  if (!length(f)) return(g)
  if (!length(g)) return(f)
  n <- max(length(f), length(g))
  f <- .dm_pad(f, n, u); g <- .dm_pad(g, n, u)
  .dm_strip(lapply(seq_len(n), function(i) .dm_add(f[[i]], g[[i]], u - 1L)), u)
}
.dm_neg <- function(f, u) if (u == 0L) -f else lapply(f, .dm_neg, u = u - 1L)
.dm_sub <- function(f, g, u) .dm_add(f, .dm_neg(g, u), u)
.dm_mul <- function(f, g, u) {
  if (u == 0L) return(.zp_mul(f, g))
  if (!length(f) || !length(g)) return(list())
  df <- length(f); dg <- length(g)
  out <- rep(list(.dm_zero(u - 1L)), df + dg - 1L)
  for (i in seq_len(df)) {
    if (.dm_zero_p(f[[i]], u - 1L)) next
    for (j in seq_len(dg)) out[[i + j - 1L]] <- .dm_add(out[[i + j - 1L]], .dm_mul(f[[i]], g[[j]], u - 1L), u - 1L)
  }
  .dm_strip(out, u)
}
.dm_mul_ground <- function(f, c, u) {
  if (u == 0L) return(.zp_trim(f * c))
  .dm_strip(lapply(f, .dm_mul_ground, c = c, u = u - 1L), u)
}
.dm_quo_ground <- function(f, c, u) {
  if (u == 0L) return(.zp_trim(gmp::divq.bigz(f, c)))
  .dm_strip(lapply(f, .dm_quo_ground, c = c, u = u - 1L), u)
}
.dm_pow <- function(f, n, u) { r <- .dm_one(u); for (i in seq_len(n)) r <- .dm_mul(r, f, u); r }
.dm_sub_mul <- function(f, g, h, u) .dm_sub(f, .dm_mul(g, h, u), u)
.dm_add_mul <- function(f, g, h, u) .dm_add(f, .dm_mul(g, h, u), u)
.dm_expand <- function(F, u) { r <- .dm_one(u); for (f in F) r <- .dm_mul(r, f, u); r }
.dm_max_norm <- function(f, u) { cs <- .dm_ground_coeffs(f, u); if (length(cs)) max(abs(cs)) else .dz(0L) }
.dm_ground_trunc <- function(f, p, u) {
  if (u == 0L) return(.zp_trunc(f, p))
  .dm_strip(lapply(f, .dm_ground_trunc, p = p, u = u - 1L), u)
}
# exact division in Z[x_0, ..., x_u] (dmp_exquo / dmp_quo for exact cases)
.dm_quo <- function(f, g, u) {
  if (u == 0L) {
    dv <- .zp_div_exact(f, g)
    return(dv)
  }
  q <- list(); r <- f
  dg <- .dm_degree(g, u); lcg <- .dm_LC(g, u)
  while (!.dm_zero_p(r, u) && .dm_degree(r, u) >= dg) {
    lcr <- .dm_LC(r, u)
    c <- .dm_quo(lcr, lcg, u - 1L)
    j <- .dm_degree(r, u) - dg
    term <- c(list(c), rep(list(.dm_zero(u - 1L)), j))
    q <- .dm_add(q, term, u)
    r <- .dm_sub(r, .dm_mul(term, g, u), u)
  }
  q
}
# exact univariate integer division (the quotient of dup_div when the
# division is exact)
.zp_div_exact <- function(a, b) {
  a <- .zp_trim(a); b <- .zp_trim(b)
  db <- length(b) - 1L
  if (length(a) - 1L < db) return(.dz0())
  r <- a; nq <- length(a) - db
  q <- .dz(rep(0L, nq))
  lc <- b[1L]
  for (i in seq_len(nq)) {
    cf <- gmp::divq.bigz(r[i], lc)
    q[i] <- cf
    if (cf != 0) { idx <- i:(i + db); r[idx] <- r[idx] - b * cf }
  }
  .zp_trim(q)
}

# dmp_eval_in, dmp_eval_tail, dmp_diff_eval_in
.dup_eval <- function(f, a) { r <- .dz(0L); for (i in seq_along(f)) r <- r * a + f[i]; r }
.dm_eval <- function(f, a, u) {
  # evaluate the main variable at a
  if (u == 0L) return(.dup_eval(f, a))
  r <- .dm_zero(u - 1L)
  for (c in f) r <- .dm_add(.dm_mul_ground(r, a, u - 1L), c, u - 1L)
  r
}
.dm_eval_in <- function(f, a, j, u) {
  if (j == 0L) return(.dm_eval(f, a, u))
  r <- lapply(f, .dm_eval_in, a = a, j = j - 1L, u = u - 1L)
  if (u - 1L == 0L) return(.zp_trim(do.call(c, c(list(.dz0()), lapply(r, function(z) if (length(z)) z else .dz(0L))))))
  .dm_strip(r, u - 1L)
}
.dm_eval_tail <- function(f, A, u) {
  # dmp_eval_tail with _rec_eval_tail
  if (!length(A)) return(f)
  if (.dm_zero_p(f, u)) return(.dm_zero(u - length(A)))
  nA <- length(A)
  rec <- function(g, i) {
    if (i == u) return(.dup_eval(g, A[[nA]]))
    h <- lapply(g, rec, i = i + 1L)
    if (i < u - nA) return(h)
    if (i == u - nA) return(.zp_trim(do.call(c, c(list(.dz0()), lapply(h, function(z) if (length(z)) z else .dz(0L))))))
    .dup_eval(do.call(c, c(list(.dz0()), lapply(h, function(z) if (length(z)) z else .dz(0L)))), A[[nA - u + i]])
  }
  e <- rec(f, 0L)
  if (u == nA - 1L) return(e)
  .dm_strip(e, u - nA)
}
.dm_diff <- function(f, m, u) {
  # m-th derivative in the main variable
  n <- .dm_degree(f, u)
  if (n < m) return(.dm_zero(u))
  out <- list()
  deg <- n
  res <- vector("list", 0)
  for (i in seq_len(length(f) - m)) {
    k <- deg - (i - 1L)
    fac <- .dz(1L); for (t in 0:(m - 1L)) fac <- fac * (k - t)
    res[[i]] <- if (u == 0L) f[i] * fac else .dm_mul_ground(f[[i]], fac, u - 1L)
  }
  if (u == 0L) .zp_trim(do.call(c, c(list(.dz0()), res))) else .dm_strip(res, u)
}
.dm_diff_in <- function(f, m, j, u) {
  if (j == 0L) return(.dm_diff(f, m, u))
  .dm_strip(lapply(f, .dm_diff_in, m = m, j = j - 1L, u = u - 1L), u)
}
.dm_diff_eval_in <- function(f, m, a, j, u) .dm_eval_in(.dm_diff_in(f, m, j, u), a, j, u)

# content and primitive part over the ground domain
.dm_ground_content <- function(f, u) { cs <- .dm_ground_coeffs(f, u); cs <- cs[cs != 0]; if (!length(cs)) return(.dz(0L)); g <- abs(cs[1L]); for (c in as.list(cs[-1L])) g <- gmp::gcd.bigz(g, c); g }
.dm_ground_primitive <- function(f, u) {
  c <- .dm_ground_content(f, u)
  if (c == 0 || c == 1) return(list(c, f))
  list(c, .dm_quo_ground(f, c, u))
}

# GF(p) helpers with bigz coefficients for the Diophantine equations
.gz_trim <- function(a) .zp_trim(a)
.gz_mod <- function(a, p) .gz_trim(gmp::mod.bigz(a, p))
.gz_add <- function(a, b, p) .gz_mod(.zp_add(a, b), p)
.gz_sub <- function(a, b, p) .gz_mod(.zp_add(a, -b), p)
.gz_mul <- function(a, b, p) .gz_mod(.zp_mul(a, b), p)
.gz_divmod <- function(a, b, p) {
  a <- .gz_mod(a, p); b <- .gz_mod(b, p)
  db <- length(b) - 1L
  if (length(a) - 1L < db) return(list(q = .dz0(), r = a))
  inv <- gmp::inv.bigz(b[1L], p)
  r <- a; nq <- length(a) - db
  q <- .dz(rep(0L, nq))
  for (i in seq_len(nq)) {
    cf <- gmp::mod.bigz(r[i] * inv, p)
    q[i] <- cf
    if (cf != 0) { idx <- i:(i + db); r[idx] <- gmp::mod.bigz(r[idx] - b * cf, p) }
  }
  rr <- if (db > 0L) r[(length(r) - db + 1L):length(r)] else .dz0()
  list(q = .gz_trim(q), r = .gz_trim(rr))
}
.gz_monic <- function(f, p) {
  if (!length(f)) return(list(.dz(0L), f))
  lc <- gmp::mod.bigz(f[1L], p)
  list(lc, .gz_mod(f * gmp::inv.bigz(lc, p), p))
}
.gz_gcdex <- function(f, g, p) {
  # gf_gcdex: s*f + t*g = h with h monic
  if (!length(f) && !length(g)) return(list(.dz(1L), .dz0(), .dz0()))
  m0 <- .gz_monic(f, p); m1 <- .gz_monic(g, p)
  p0 <- m0[[1L]]; r0 <- m0[[2L]]; p1 <- m1[[1L]]; r1 <- m1[[2L]]
  if (!length(f)) return(list(.dz0(), .dz(gmp::inv.bigz(p1, p)), r1))
  if (!length(g)) return(list(.dz(gmp::inv.bigz(p0, p)), .dz0(), r0))
  s0 <- .dz(gmp::inv.bigz(p0, p)); s1 <- .dz0()
  t0 <- .dz0(); t1 <- .dz(gmp::inv.bigz(p1, p))
  repeat {
    qr <- .gz_divmod(r0, r1, p)
    Q <- qr$q; R <- qr$r
    if (!length(R)) break
    mr <- .gz_monic(R, p)
    lc <- mr[[1L]]; r0 <- r1; r1 <- mr[[2L]]
    inv <- gmp::inv.bigz(lc, p)
    s <- .gz_sub(s0, .gz_mul(s1, Q, p), p)
    t <- .gz_sub(t0, .gz_mul(t1, Q, p), p)
    s0 <- s1; s1 <- .gz_mod(s * inv, p)
    t0 <- t1; t1 <- .gz_mod(t * inv, p)
  }
  list(s1, t1, r1)
}
.gz_lshift <- function(f, n) if (!length(f)) f else c(f, .dz(rep(0L, n)))
.gz_to_int <- function(f, p) .zp_trunc(f, p)

# dup_zz_diophantine and dmp_zz_diophantine
.dup_zz_diophantine <- function(F, m, p) {
  if (length(F) == 2L) {
    f <- .gz_mod(F[[1L]], p); g <- .gz_mod(F[[2L]], p)
    gx <- .gz_gcdex(g, f, p)
    s <- .gz_lshift(gx[[1L]], m); t <- .gz_lshift(gx[[2L]], m)
    qr <- .gz_divmod(s, f, p)
    s <- qr$r
    t <- .gz_add(t, .gz_mul(qr$q, g, p), p)
    return(list(.gz_to_int(s, p), .gz_to_int(t, p)))
  }
  G <- list(F[[length(F)]])
  if (length(F) > 2L) for (f in rev(F[2:(length(F) - 1L)])) G <- c(list(.zp_mul(f, G[[1L]])), G)
  S <- list(); Tl <- list(.dz(1L))
  for (i in seq_along(G)) {
    r <- .dm_zz_diophantine(list(G[[i]], F[[i]]), Tl[[length(Tl)]], list(), 0L, p, 0L)
    Tl[[length(Tl) + 1L]] <- r[[1L]]; S[[length(S) + 1L]] <- r[[2L]]
  }
  S <- c(S, list(Tl[[length(Tl)]]))
  lapply(seq_along(F), function(i) .gz_to_int(.gz_divmod(.gz_lshift(.gz_mod(S[[i]], p), m), .gz_mod(F[[i]], p), p)$r, p))
}
.dm_zz_diophantine <- function(F, c, A, d, p, u) {
  if (!length(A)) {
    S <- rep(list(.dz0()), length(F))
    n <- length(c) - 1L
    for (i in seq_along(c)) {
      cf <- c[i]
      if (cf == 0) next
      Tt <- .dup_zz_diophantine(F, n - (i - 1L), p)
      for (j in seq_along(S)) S[[j]] <- .zp_trunc(.zp_add(S[[j]], .zp_trim(Tt[[j]] * cf)), p)
    }
    return(S)
  }
  n <- length(A)
  e <- .dm_expand(F, u)
  a <- A[[n]]; A <- A[-n]
  B <- lapply(F, function(f) .dm_quo(e, f, u))
  G <- lapply(F, function(f) .dm_eval_in(f, a, n, u))
  C <- .dm_eval_in(c, a, n, u)
  v <- u - 1L
  S <- .dm_zz_diophantine(G, C, A, d, p, v)
  S <- lapply(S, .dm_raise, l = 1L, u = v)
  for (i in seq_along(S)) c <- .dm_sub_mul(c, S[[i]], B[[i]], u)
  c <- .dm_ground_trunc(c, p, u)
  m1 <- .dz(c(1L, 0L)); m1[2L] <- -a
  m <- .dm_nest(m1, n)
  M <- .dm_one(n)
  for (k in seq_len(d) - 1L) {
    if (.dm_zero_p(c, u)) break
    M <- .dm_mul(M, m, u)
    C <- .dm_diff_eval_in(c, k + 1L, a, n, u)
    if (!.dm_zero_p(C, v)) {
      C <- .dm_quo_ground(C, gmp::factorialZ(k + 1L), v)
      Tt <- .dm_zz_diophantine(G, C, A, d, p, v)
      Tt <- lapply(Tt, function(t) .dm_mul(.dm_raise(t, 1L, v), M, u))
      for (i in seq_along(S)) S[[i]] <- .dm_add(S[[i]], Tt[[i]], u)
      for (i in seq_along(Tt)) c <- .dm_sub_mul(c, Tt[[i]], B[[i]], u)
      c <- .dm_ground_trunc(c, p, u)
    }
  }
  lapply(S, .dm_ground_trunc, p = p, u = u)
}

# Wang's EEZ algorithm
.dup_sqf_p <- function(f) {
  if (length(f) <= 1L) return(TRUE)
  df <- .zp_trim(f[-length(f)] * gmp::as.bigz((length(f) - 1L):1))
  g <- up_gcd(gmp::as.bigq(f), gmp::as.bigq(df))
  length(g) <= 1L
}
.dup_primitive <- function(f) { pr <- .zp_primitive(f); list(pr$content, pr$poly) }
.wang_eval_failed <- function() stop(structure(class = c("EvaluationFailed", "error", "condition"), list(message = "no luck", call = NULL)))
.wang_extraneous <- function() stop(structure(class = c("ExtraneousFactors", "error", "condition"), list(message = "extraneous factors", call = NULL)))

.dm_zz_wang_non_divisors <- function(E, cs, ct) {
  result <- list(cs * ct)
  for (q in E) {
    q <- abs(q)
    for (r in rev(result)) {
      while (r != 1) { r <- gmp::gcd.bigz(r, q); q <- gmp::divq.bigz(q, r) }
      if (q == 1) return(NULL)
    }
    result[[length(result) + 1L]] <- q
  }
  result[-1L]
}
.dm_zz_wang_test_points <- function(f, Tl, ct, A, u) {
  lt <- .dm_eval_tail(.dm_LC(f, u), A, u - 1L)
  if (length(lt) == 0L || all(lt == 0)) .wang_eval_failed()
  g <- .dm_eval_tail(f, A, u)
  if (!.dup_sqf_p(g)) .wang_eval_failed()
  ch <- .dup_primitive(g); c <- ch[[1L]]; h <- ch[[2L]]
  if (h[1L] < 0) { c <- -c; h <- -h }
  v <- u - 1L
  E <- lapply(Tl, function(tk) .dm_eval_tail(tk[[1L]], A, v))
  E <- lapply(E, function(z) if (length(z)) z[1L] else .dz(0L))
  D <- .dm_zz_wang_non_divisors(E, c, ct)
  if (is.null(D)) .wang_eval_failed()
  list(c, h, E)
}
.dm_zz_wang_lead_coeffs <- function(f, Tl, cs, E, H, A, u) {
  v <- u - 1L
  C <- list(); J <- rep(0L, length(E))
  for (h in H) {
    c <- .dm_one(v)
    d <- h[1L] * cs
    for (i in rev(seq_along(E))) {
      k <- 0L; e <- E[[i]]; tt <- Tl[[i]][[1L]]
      while (gmp::mod.bigz(d, e) == 0) { d <- gmp::divq.bigz(d, e); k <- k + 1L }
      if (k != 0L) { c <- .dm_mul(c, .dm_pow(tt, k, v), v); J[i] <- 1L }
    }
    C[[length(C) + 1L]] <- c
  }
  if (!all(J == 1L)) .wang_extraneous()
  CC <- list(); HH <- list()
  for (i in seq_along(C)) {
    c <- C[[i]]; h <- H[[i]]
    d <- .dm_eval_tail(c, A, v); d <- if (length(d)) d[1L] else .dz(0L)
    lc <- h[1L]
    if (cs == 1) cc <- gmp::divq.bigz(lc, d)
    else {
      g <- gmp::gcd.bigz(lc, d)
      d <- gmp::divq.bigz(d, g); cc <- gmp::divq.bigz(lc, g)
      h <- .zp_trim(h * d); cs <- gmp::divq.bigz(cs, d)
    }
    CC[[length(CC) + 1L]] <- .dm_mul_ground(c, cc, v)
    HH[[length(HH) + 1L]] <- h
  }
  if (cs == 1) return(list(f, HH, CC))
  CCC <- lapply(CC, .dm_mul_ground, c = cs, u = v)
  HHH <- lapply(HH, function(h) .zp_trim(h * cs))
  f <- .dm_mul_ground(f, cs^(length(H) - 1L), u)
  list(f, HHH, CCC)
}
.dm_zz_wang_hensel_lifting <- function(f, H, LC, A, p, u) {
  S <- list(f); n <- length(A); v <- u - 1L
  if (n > 1L) {
    rA <- rev(A[-1L])
    for (i in seq_along(rA)) {
      s <- .dm_eval_in(S[[1L]], rA[[i]], n - (i - 1L), u - (i - 1L))
      S <- c(list(.dm_ground_trunc(s, p, v - (i - 1L))), S)
    }
  }
  d <- max(.dm_degree_list(f, u)[-1L])
  for (jj in seq_len(n)) {
    j <- jj + 1L
    s <- S[[jj]]; a <- A[[jj]]
    G <- H; w <- j - 1L
    I <- if (j - 2L > 0L) A[seq_len(j - 2L)] else list()
    Jt <- if (j - 1L < n) A[j:n] else list()
    for (i in seq_along(H)) {
      lc <- .dm_ground_trunc(.dm_eval_tail(LC[[i]], Jt, v), p, w - 1L)
      h <- H[[i]]
      rest <- if (length(h) > 1L) .dm_raise(h[-1L], 1L, w - 1L) else list()
      H[[i]] <- c(list(lc), rest)
    }
    m1 <- .dz(c(1L, 0L)); m1[2L] <- -a
    m <- .dm_nest(m1, w)
    M <- .dm_one(w)
    c <- .dm_sub(s, .dm_expand(H, w), w)
    dj <- .dm_degree_in(s, w, w)
    for (k in seq_len(max(dj, 0)) - 1L) {
      if (.dm_zero_p(c, w)) break
      M <- .dm_mul(M, m, w)
      C <- .dm_diff_eval_in(c, k + 1L, a, w, w)
      if (!.dm_zero_p(C, w - 1L)) {
        C <- .dm_quo_ground(C, gmp::factorialZ(k + 1L), w - 1L)
        Tt <- .dm_zz_diophantine(G, C, I, d, p, w - 1L)
        for (i in seq_along(H)) {
          h <- .dm_add_mul(H[[i]], .dm_raise(Tt[[i]], 1L, w - 1L), M, w)
          H[[i]] <- .dm_ground_trunc(h, p, w)
        }
        h <- .dm_sub(s, .dm_expand(H, w), w)
        c <- .dm_ground_trunc(h, p, w)
      }
    }
  }
  if (!.dm_equal(.dm_expand(H, u), f, u)) .wang_extraneous()
  H
}
.dm_equal <- function(f, g, u) {
  if (u == 0L) return(length(f) == length(g) && all(f == g))
  length(f) == length(g) && all(vapply(seq_along(f), function(i) .dm_equal(f[[i]], g[[i]], u - 1L), NA))
}

# dmp_zz_wang; fac_lc(g) factors a level u - 1 leading coefficient into
# list(content, list(list(factor, multiplicity), ...))
.dm_zz_wang <- function(f, u, fac_lc, mod = NULL) {
  ft <- fac_lc(.dm_LC(f, u), u - 1L)
  ct <- ft[[1L]]; Tl <- ft[[2L]]
  b <- .dm_zz_mignotte_bound(f, u)
  p <- gmp::nextprime(b)
  if (is.null(mod)) mod <- if (u == 1L) 2L else 1L
  history <- character(); configs <- list(); A <- rep(list(.dz(0L)), u); r <- NULL
  tp <- tryCatch(.dm_zz_wang_test_points(f, Tl, ct, A, u), EvaluationFailed = function(err) NULL)
  if (!is.null(tp)) {
    H <- .zp_wang_factor_sqf(tp[[2L]])
    r <- length(H)
    if (r == 1L) return(list(f))
    configs <- list(list(tp[[2L]], tp[[1L]], tp[[3L]], H, A))
  }
  nconf <- 3L; ntries <- 5L; mstep <- 2L
  while (length(configs) < nconf) {
    full <- FALSE
    for (try in seq_len(ntries)) {
      A <- lapply(seq_len(u), function(i) .dz(.wang_randint(-mod, mod)))
      key <- paste(vapply(A, as.character, ""), collapse = ",")
      if (key %in% history) next
      history <- c(history, key)
      tp <- tryCatch(.dm_zz_wang_test_points(f, Tl, ct, A, u), EvaluationFailed = function(err) NULL)
      if (is.null(tp)) next
      H <- .zp_wang_factor_sqf(tp[[2L]])
      rr <- length(H)
      if (!is.null(r)) {
        if (rr != r) { if (rr < r) { configs <- list(); r <- rr } else next }
      } else r <- rr
      if (r == 1L) return(list(f))
      configs[[length(configs) + 1L]] <- list(tp[[2L]], tp[[1L]], tp[[3L]], H, A)
      if (length(configs) == nconf) { full <- TRUE; break }
    }
    if (!full) mod <- mod + mstep
  }
  norms <- vapply(configs, function(cf) as.character(max(abs(cf[[1L]]))), "")
  s_arg <- 1L; best <- gmp::as.bigz(norms[1L])
  for (i in seq_along(configs)[-1L]) { nv <- gmp::as.bigz(norms[i]); if (nv < best) { best <- nv; s_arg <- i } }
  cf <- configs[[s_arg]]
  cs <- cf[[2L]]; E <- cf[[3L]]; H <- cf[[4L]]; A <- cf[[5L]]
  res <- tryCatch({
    lcr <- .dm_zz_wang_lead_coeffs(f, Tl, cs, E, H, A, u)
    .dm_zz_wang_hensel_lifting(lcr[[1L]], lcr[[2L]], lcr[[3L]], A, p, u)
  }, ExtraneousFactors = function(err) NULL)
  if (is.null(res)) return(.dm_zz_wang(f, u, fac_lc, mod + 1L))
  lapply(res, function(g) {
    g <- .dm_ground_primitive(g, u)[[2L]]
    if (.dm_ground_LC(g, u) < 0) g <- .dm_neg(g, u)
    g
  })
}
.dm_zz_mignotte_bound <- function(f, u) {
  a <- .dm_max_norm(f, u)
  b <- abs(.dm_ground_LC(f, u))
  n <- sum(.dm_degree_list(f, u))
  gmp::as.bigz(floor(sqrt(n + 1))) * gmp::as.bigz(2)^n * a * b
}
.zp_wang_factor_sqf <- function(s) {
  g <- .dup_primitive(s)[[2L]]
  if (g[1L] < 0) g <- -g
  if (length(g) <= 1L) return(list())
  if (length(g) == 2L) return(list(g))
  .zp_zassenhaus(g)
}

# a deterministic generator for the evaluation points (the factorization
# does not depend on them; R's global random state is left alone)
.wang_state <- new.env()
.wang_randint <- function(lo, hi) {
  st <- .wang_state$x %||% 12345
  st <- (st * 1103515245 + 12345) %% 2147483648
  .wang_state$x <- st
  lo + floor(st / 2147483648 * (hi - lo + 1))
}

# qpoly <-> dense recursive polynomial in the variables ord (main first)
.qp_to_dm <- function(f, ord) {
  e <- f$e[, ord, drop = FALSE]
  cf <- gmp::numerator(f$c)
  u <- length(ord) - 1L
  rec <- function(rows, k) {
    if (!length(rows)) return(.dm_zero(u - k + 1L))
    dk <- e[rows, k]
    D <- max(dk)
    if (k == length(ord)) {
      out <- gmp::as.bigz(rep(0L, D + 1L))
      for (r in rows) out[D - e[r, k] + 1L] <- out[D - e[r, k] + 1L] + cf[r]
      return(.zp_trim(out))
    }
    .dm_strip(lapply(D:0, function(d) rec(rows[dk == d], k + 1L)), u - k + 1L)
  }
  rec(seq_len(nrow(e)), 1L)
}
.dm_to_qp <- function(g, ord, n) {
  u <- length(ord) - 1L
  exps <- list(); cs <- list()
  rec <- function(h, k, pre) {
    if (k == length(ord)) {
      d <- length(h) - 1L
      for (i in seq_along(h)) if (h[i] != 0) { ex <- integer(n); ex[ord] <- c(pre, d - i + 1L); exps[[length(exps) + 1L]] <<- ex; cs[[length(cs) + 1L]] <<- h[i] }
      return(invisible())
    }
    d <- length(h) - 1L
    for (i in seq_along(h)) rec(h[[i]], k + 1L, c(pre, d - i + 1L))
  }
  rec(g, 1L, integer(0))
  if (!length(exps)) return(qp_const(0, n))
  qp_norm(do.call(rbind, exps), gmp::as.bigq(do.call(c, cs)))
}

# factor a primitive square-free integer qpoly with v as main variable
.wang_factor <- function(h, v) {
  n <- qp_nvars(h)
  others <- setdiff(qp_vars_used(h), v)
  ord <- c(v, others)
  u <- length(ord) - 1L
  f <- .qp_to_dm(h, ord)
  sub_ord <- others
  fac_lc <- function(lc, lu) {
    q <- .dm_to_qp(lc, sub_ord, n)
    fl <- qp_factor_list(q)
    Tl <- list()
    for (pr in fl$factors) Tl[[length(Tl) + 1L]] <- list(.qp_to_dm(pr[[1L]], sub_ord), pr[[2L]])
    list(gmp::numerator(fl$coef), Tl)
  }
  H <- .dm_zz_wang(f, u, fac_lc)
  lapply(H, .dm_to_qp, ord = ord, n = n)
}
