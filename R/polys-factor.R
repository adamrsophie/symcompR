# Factorization over the integers, after sympy/polys/factortools.py
# (Zassenhaus with Hensel lifting) and sympy/polys/euclidtools.py.
# Integer polynomials are bigz vectors, highest degree first.

.zp_trim <- function(a) {
  if (!length(a)) return(a)
  nz <- which(!(a == 0))
  if (!length(nz)) return(gmp::as.bigz(integer(0)))
  a[nz[1L]:length(a)]
}

.zp_sub <- function(a, b) .zp_add(a, -b)

.zp_trunc <- function(a, m) {
  if (!length(a)) return(a)
  r <- gmp::mod.bigz(a, m)
  half <- gmp::divq.bigz(m, 2)
  big <- r > half
  if (any(big)) r[big] <- r[big] - m
  .zp_trim(r)
}

# Division by a monic integer polynomial.
.zp_div_monic <- function(a, b) {
  db <- length(b) - 1L
  if (length(a) - 1L < db) return(list(q = gmp::as.bigz(integer(0)), r = a))
  r <- a
  nq <- length(a) - db
  q <- gmp::as.bigz(rep(0, nq))
  for (i in seq_len(nq)) {
    cf <- r[i]
    q[i] <- cf
    if (cf != 0) {
      idx <- i:(i + db)
      r[idx] <- r[idx] - b * cf
    }
  }
  rr <- if (db > 0L) r[(length(r) - db + 1L):length(r)] else gmp::as.bigz(integer(0))
  list(q = .zp_trim(q), r = .zp_trim(rr))
}

.zp_from_gf <- function(a, p) {
  v <- ifelse(a > p / 2, a - p, a)
  gmp::as.bigz(v)
}

.zp_max_norm <- function(a) max(abs(a))
.zp_l1_norm <- function(a) sum(abs(a))

.zp_content <- function(a) {
  g <- abs(a[1L])
  if (length(a) > 1L) for (k in 2:length(a)) {
    g <- gmp::gcd.bigz(g, a[k])
    if (g == 1) break
  }
  g
}

.zp_primitive <- function(a) {
  a <- .zp_trim(a)
  if (!length(a)) return(list(content = gmp::as.bigz(0), poly = a))
  g <- .zp_content(a)
  if (a[1L] < 0) g <- -g
  list(content = g, poly = gmp::divq.bigz(a, g))
}

# Exact division over the integers (b need not be monic); NULL if inexact.
.zp_exquo <- function(a, b) {
  qa <- gmp::as.bigq(a)
  qb <- gmp::as.bigq(b)
  d <- up_divmod(qa, qb)
  if (length(d$r)) return(NULL)
  if (any(gmp::denominator(d$q) != 1)) return(NULL)
  gmp::numerator(d$q)
}

.zp_hensel_step <- function(m, f, g, h, s, t) {
  M <- m * m
  e <- .zp_trunc(.zp_sub(f, .zp_mul(g, h)), M)
  dv <- .zp_div_monic(.zp_mul(s, e), h)
  q <- .zp_trunc(dv$q, M)
  r <- .zp_trunc(dv$r, M)
  u <- .zp_add(.zp_mul(t, e), .zp_mul(q, g))
  G <- .zp_trunc(.zp_add(g, u), M)
  H <- .zp_trunc(.zp_add(h, r), M)
  u <- .zp_add(.zp_mul(s, G), .zp_mul(t, H))
  b <- .zp_trunc(.zp_sub(u, gmp::as.bigz(1)), M)
  dv <- .zp_div_monic(.zp_mul(s, b), H)
  cc <- .zp_trunc(dv$q, M)
  d <- .zp_trunc(dv$r, M)
  u <- .zp_add(.zp_mul(t, b), .zp_mul(cc, G))
  Sn <- .zp_trunc(.zp_sub(s, d), M)
  Tn <- .zp_trunc(.zp_sub(t, u), M)
  list(G = G, H = H, S = Sn, T = Tn)
}

.zp_hensel_lift <- function(p, f, flist, l) {
  r <- length(flist)
  lc <- f[1L]
  pl <- gmp::as.bigz(p)^l
  if (r == 1L) {
    inv <- gmp::inv.bigz(gmp::mod.bigz(lc, pl), pl)
    return(list(.zp_trunc(f * inv, pl)))
  }
  m <- gmp::as.bigz(p)
  k <- r %/% 2L
  d <- ceiling(log2(l))
  g <- .gf_from_bigz(lc, p)
  for (fi in flist[seq_len(k)]) g <- .gf_mul(g, .gf_from_bigz(fi, p), p)
  h <- .gf_from_bigz(flist[[k + 1L]], p)
  if (k + 2L <= r) for (fi in flist[(k + 2L):r]) h <- .gf_mul(h, .gf_from_bigz(fi, p), p)
  ex <- .gf_gcdex(g, h, p)
  g <- .zp_from_gf(g, p)
  h <- .zp_from_gf(h, p)
  s <- .zp_from_gf(ex$s, p)
  t <- .zp_from_gf(ex$t, p)
  for (i in seq_len(d)) {
    st <- .zp_hensel_step(m, f, g, h, s, t)
    g <- st$G; h <- st$H; s <- st$S; t <- st$T
    m <- m * m
  }
  c(.zp_hensel_lift(p, g, flist[seq_len(k)], l),
    .zp_hensel_lift(p, h, flist[(k + 1L):r], l))
}

.combn_list <- function(n, k) {
  if (k > n) return(list())
  cm <- utils::combn(n, k)
  lapply(seq_len(ncol(cm)), function(j) cm[, j])
}

# Zassenhaus factorization of a square-free primitive integer polynomial.
.zp_zassenhaus <- function(f) {
  n <- length(f) - 1L
  if (n <= 1L) return(list(f))
  fc <- f[length(f)]
  A <- .zp_max_norm(f)
  b <- f[1L]
  B <- gmp::as.bigz(ceiling(sqrt(n + 1) * 1024)) * gmp::as.bigz(2)^n * A * abs(b) / 1024
  B <- gmp::as.bigz(B) + 1
  cands <- list()
  for (px in .small_primes[.small_primes > 2]) {
    if (gmp::mod.bigz(b, px) == 0) next
    Fp <- .gf_from_bigz(f, px)
    if (length(Fp) - 1L != n) next
    if (!.gf_sqf_p(Fp, px)) next
    fs <- .gf_factor_sqf(Fp, px)
    cands[[length(cands) + 1L]] <- list(px, fs)
    if (length(fs) < 15L || length(cands) > 4L) break
  }
  best <- cands[[which.min(vapply(cands, function(cd) length(cd[[2L]]), 0L))]]
  p <- best[[1L]]
  fsqf <- best[[2L]]
  if (length(fsqf) == 1L) return(list(f))
  l <- ceiling(log(2 * as.double(B) + 1, p))
  if (!is.finite(l)) l <- ceiling(as.double(log2(2 * B + 1)) / log2(p))
  modular <- lapply(fsqf, .zp_from_gf, p = p)
  g <- .zp_hensel_lift(p, f, modular, l)
  Tset <- seq_along(g)
  factors <- list()
  s <- 1L
  pl <- gmp::as.bigz(p)^l
  half <- gmp::divq.bigz(pl, 2)
  while (2L * s <= length(Tset)) {
    found <- FALSE
    for (Sidx in .combn_list(length(Tset), s)) {
      S <- Tset[Sidx]
      if (b == 1) {
        q <- gmp::as.bigz(1)
        for (i in S) q <- q * g[[i]][length(g[[i]])]
        q <- gmp::mod.bigz(q, pl)
        if (q > half) q <- q - pl
        if (q != 0 && gmp::mod.bigz(fc, q) != 0) next
        G <- gmp::as.bigz(1)
        for (i in S) G <- .zp_mul(G, g[[i]])
        G <- .zp_trunc(G, pl)
      } else {
        G <- b
        for (i in S) G <- .zp_mul(G, g[[i]])
        G <- .zp_trunc(G, pl)
        G <- .zp_primitive(G)$poly
        q <- G[length(G)]
        if (q != 0 && gmp::mod.bigz(fc, q) != 0) next
      }
      H <- b
      TS <- setdiff(Tset, S)
      for (i in TS) H <- .zp_mul(H, g[[i]])
      H <- .zp_trunc(H, pl)
      if (.zp_l1_norm(G) * .zp_l1_norm(H) <= B) {
        Tset <- TS
        G <- .zp_primitive(G)$poly
        f <- .zp_primitive(H)$poly
        factors[[length(factors) + 1L]] <- G
        b <- f[1L]
        fc <- f[length(f)]
        found <- TRUE
        break
      }
    }
    if (!found) s <- s + 1L
  }
  c(factors, list(f))
}

# Full factorization of an integer polynomial: list(content, factors) where
# factors is a list of list(bigz poly, multiplicity).
.zp_factor_list <- function(f) {
  f <- .zp_trim(gmp::as.bigz(f))
  pr <- .zp_primitive(f)
  cont <- pr$content
  f <- pr$poly
  if (length(f) <= 1L) return(list(content = cont * if (length(f)) f[1L] else 1, factors = list()))
  out <- list()
  k <- 0L
  while (f[length(f)] == 0) {
    f <- f[-length(f)]
    k <- k + 1L
  }
  if (k) out[[length(out) + 1L]] <- list(gmp::as.bigz(c(1, 0)), k)
  if (length(f) > 1L) {
    sq <- up_sqf_list(gmp::as.bigq(f))
    for (pair in sq) {
      g <- gmp::numerator(pair[[1L]] * gmp::as.bigq(.lcm_den(pair[[1L]])))
      g <- .zp_primitive(g)$poly
      fs <- if (length(g) <= 2L) list(g) else .zp_zassenhaus(g)
      for (h in fs) out[[length(out) + 1L]] <- list(.zp_primitive(h)$poly, pair[[2L]])
    }
    # Fix the overall sign and unit content.
    prod_all <- gmp::as.bigz(1)
    for (pr2 in out) for (r in seq_len(pr2[[2L]])) prod_all <- .zp_mul(prod_all, pr2[[1L]])
    ratio <- gmp::as.bigq(f[1L]) / gmp::as.bigq(prod_all[1L])
    cont <- cont * gmp::numerator(ratio)
  }
  ord <- order(vapply(out, function(pp) length(pp[[1L]]), 0L),
               vapply(out, function(pp) paste(as.character(pp[[1L]]), collapse = ","), ""))
  list(content = cont, factors = out[ord])
}

.lcm_den <- function(q) {
  d <- gmp::denominator(q)
  l <- d[1L]
  if (length(d) > 1L) for (k in 2:length(d)) l <- gmp::lcm.bigz(l, d[k])
  l
}

# ---------------------------------------------------------------------------
# Heuristic GCD (sympy dmp_zz_heu_gcd) for integer qpolys, with a
# primitive PRS fallback.

.qp_int <- function(p) {
  # integer qpoly with positive leading coefficient, and its content
  pr <- qp_primitive(p)
  pr
}

.qp_max_norm <- function(p) if (qp_is_zero(p)) gmp::as.bigz(0) else max(abs(gmp::numerator(p$c)))

.qp_interpolate <- function(h, x, i, n) {
  # h: integer qpoly without variable i; write coefficients in base x
  out <- qp_zero(n)
  k <- 0L
  half <- gmp::divq.bigz(x, 2)
  while (!qp_is_zero(h)) {
    num <- gmp::numerator(h$c)
    g <- gmp::mod.bigz(num, x)
    big <- g > half
    if (any(big)) g[big] <- g[big] - x
    gpoly <- qp_norm(h$e, gmp::as.bigq(g))
    if (!qp_is_zero(gpoly)) {
      ge <- gpoly$e
      ge[, i] <- ge[, i] + k
      out <- qp_add(out, qp(ge, gpoly$c))
    }
    h <- qp_scale(qp_sub(h, gpoly), gmp::as.bigq(1, x))
    k <- k + 1L
  }
  out
}

.heu_gcd <- function(f, g, depth = 0L) {
  n <- qp_nvars(f)
  if (qp_is_zero(f)) return(qp_primitive(g)$poly)
  if (qp_is_zero(g)) return(qp_primitive(f)$poly)
  vf <- qp_vars_used(f)
  vg <- qp_vars_used(g)
  cf <- qp_content(f)
  cg <- qp_content(g)
  gc <- gmp::as.bigq(gmp::gcd.bigz(gmp::numerator(cf), gmp::numerator(cg)),
                     gmp::lcm.bigz(gmp::denominator(cf), gmp::denominator(cg)))
  if (!length(vf) || !length(vg)) return(qp_const(gc, n))
  f <- qp_primitive(f)$poly
  g <- qp_primitive(g)$poly
  common <- intersect(vf, vg)
  if (!length(common)) {
    # gcd of polynomials in disjoint variables is the gcd of contents in
    # each variable; fall back to the content-based PRS.
    return(qp_scale(.prs_gcd(f, g), gc))
  }
  i <- common[1L]
  fn <- .qp_max_norm(f)
  gn <- .qp_max_norm(g)
  B <- 2 * min(fn, gn) + 29
  lcf <- abs(gmp::numerator(.lc_in_var(f, i)))
  lcg <- abs(gmp::numerator(.lc_in_var(g, i)))
  x <- max(min(B, 99 * gmp::as.bigz(floor(sqrt(as.double(B))))),
           2 * min(gmp::divq.bigz(fn, max(lcf)), gmp::divq.bigz(gn, max(lcg))) + 4)
  for (attempt in 1:6) {
    ff <- qp_eval_var(f, i, gmp::as.bigq(x))
    gg <- qp_eval_var(g, i, gmp::as.bigq(x))
    if (!qp_is_zero(ff) && !qp_is_zero(gg)) {
      h <- if (!length(qp_vars_used(ff)) && !length(qp_vars_used(gg))) {
        qp_const(gmp::gcd.bigz(gmp::numerator(ff$c[1L]), gmp::numerator(gg$c[1L])), n)
      } else {
        .heu_gcd(ff, gg, depth + 1L)
      }
      hh <- qp_primitive(.qp_interpolate(h, x, i, n))$poly
      if (!qp_is_zero(hh)) {
        if (!is.null(qp_exquo(f, hh)) && !is.null(qp_exquo(g, hh))) return(qp_scale(hh, gc))
      }
    }
    x <- gmp::divq.bigz(73794 * x * gmp::as.bigz(floor(sqrt(sqrt(as.double(x))))), 27011)
  }
  qp_scale(.prs_gcd(f, g), gc)
}

.lc_in_var <- function(p, i) {
  d <- max(p$e[, i])
  w <- which(p$e[, i] == d)
  p$c[w]
}

# Coefficients of p as a polynomial in variable i: list indexed by degree+1.
.coeffs_in_var <- function(p, i) {
  d <- max(p$e[, i])
  out <- vector("list", d + 1L)
  for (k in 0:d) {
    w <- which(p$e[, i] == k)
    e <- p$e[w, , drop = FALSE]
    e[, i] <- 0L
    out[[k + 1L]] <- qp_norm(e, p$c[w])
  }
  out
}

.gcd_list <- function(ps) {
  ps <- Filter(Negate(qp_is_zero), ps)
  if (!length(ps)) return(NULL)
  g <- ps[[1L]]
  if (length(ps) > 1L) for (q in ps[-1L]) {
    g <- qp_gcd(g, q)
    if (qp_is_const(g)) break
  }
  g
}

# Content of p with respect to variable i (gcd of its coefficients).
qp_content_in_var <- function(p, i) {
  g <- .gcd_list(.coeffs_in_var(p, i))
  if (is.null(g)) return(qp_zero(qp_nvars(p)))
  qp_primitive(g)$poly
}

# Primitive PRS gcd in the first common variable, recursive on contents.
.prs_gcd <- function(f, g) {
  n <- qp_nvars(f)
  vf <- qp_vars_used(f)
  vg <- qp_vars_used(g)
  if (!length(vf) || !length(vg)) return(qp_const(1, n))
  common <- intersect(vf, vg)
  if (!length(common)) {
    i <- vf[1L]
    cf <- qp_content_in_var(f, i)
    return(.prs_gcd(cf, g))
  }
  i <- common[1L]
  cf <- qp_content_in_var(f, i)
  cg <- qp_content_in_var(g, i)
  cgcd <- if (qp_is_const(cf) || qp_is_const(cg)) qp_const(1, n) else qp_gcd(cf, cg)
  a <- qp_exquo(f, cf)
  b <- qp_exquo(g, cg)
  if (qp_degree(a, i) < qp_degree(b, i)) {
    tmp <- a; a <- b; b <- tmp
  }
  while (!qp_is_zero(b) && qp_degree(b, i) > 0) {
    r <- .qp_prem(a, b, i)
    a <- b
    if (qp_is_zero(r)) {
      b <- r
      break
    }
    b <- qp_exquo(r, qp_content_in_var(r, i))
    if (is.null(b)) stop("internal error in PRS gcd")
  }
  res <- if (qp_is_zero(b)) a else qp_const(1, n)
  res <- qp_exquo(res, qp_content_in_var(res, i))
  qp_primitive(qp_mul(res, cgcd))$poly
}

# Pseudo-remainder of a by b in variable i.
.qp_prem <- function(a, b, i) {
  da <- qp_degree(a, i)
  db <- qp_degree(b, i)
  n <- qp_nvars(a)
  lcb <- .lc_poly_in_var(b, i)
  r <- a
  dr <- da
  k <- da - db + 1
  while (!qp_is_zero(r) && qp_degree(r, i) >= db) {
    dr <- qp_degree(r, i)
    lcr <- .lc_poly_in_var(r, i)
    mono <- integer(n)
    mono[i] <- dr - db
    r <- qp_sub(qp_mul(r, lcb), qp_mul(qp_mul_term(b, mono, gmp::as.bigq(1)), lcr))
    k <- k - 1
  }
  if (k > 0) r <- qp_mul(r, qp_pow(lcb, k))
  r
}

.lc_poly_in_var <- function(p, i) {
  d <- max(p$e[, i])
  w <- which(p$e[, i] == d)
  e <- p$e[w, , drop = FALSE]
  e[, i] <- 0L
  qp_norm(e, p$c[w])
}

#' @noRd
qp_gcd <- function(f, g) {
  if (qp_is_zero(f) && qp_is_zero(g)) return(f)
  if (qp_is_zero(f)) return(qp_primitive(g)$poly)
  if (qp_is_zero(g)) return(qp_primitive(f)$poly)
  vf <- qp_vars_used(f)
  vg <- qp_vars_used(g)
  if (length(vf) == 1L && length(vg) == 1L && vf == vg) {
    i <- vf
    a <- up_from_qpoly(f, i)
    b <- up_from_qpoly(g, i)
    h <- up_primitive(up_gcd(a, b))$poly
    return(up_to_qpoly(h, i, qp_nvars(f)))
  }
  r <- .heu_gcd(f, g)
  qp_primitive(r)$poly
}

# ---------------------------------------------------------------------------
# Multivariate factorization: contents, square-free parts, then Kronecker
# substitution with recombination.
# ponytail: Kronecker substitution is exponential in the number of
# univariate factors; Wang's EEZ algorithm is the upgrade path.

# Square-free decomposition of a primitive multivariate qpoly in var i.
.qp_sqf_list_var <- function(f, i) {
  out <- list()
  fp <- qp_diff(f, i)
  if (qp_is_zero(fp)) return(list(list(f, 1L)))
  g <- qp_gcd(f, fp)
  if (qp_is_const(g)) return(list(list(f, 1L)))
  w <- qp_exquo(f, g)
  y <- qp_exquo(fp, g)
  k <- 1L
  repeat {
    z <- qp_sub(y, qp_diff(w, i))
    if (qp_is_zero(z)) {
      if (!qp_is_const(w)) out[[length(out) + 1L]] <- list(qp_primitive(w)$poly, k)
      break
    }
    h <- qp_gcd(w, z)
    if (!qp_is_const(h)) out[[length(out) + 1L]] <- list(qp_primitive(h)$poly, k)
    w <- qp_exquo(w, h)
    y <- qp_exquo(z, h)
    k <- k + 1L
    if (qp_is_const(w)) break
  }
  out
}

.kronecker_factor <- function(f) {
  n <- qp_nvars(f)
  vars <- qp_vars_used(f)
  if (length(vars) == 1L) {
    i <- vars
    zl <- .zp_factor_list(gmp::numerator(up_from_qpoly(f, i)))
    res <- list()
    for (pr in zl$factors) for (r in seq_len(pr[[2L]]))
      res[[length(res) + 1L]] <- up_to_qpoly(gmp::as.bigq(pr[[1L]]), i, n)
    return(res)
  }
  degs <- vapply(vars, function(v) qp_degree(f, v), 0)
  D <- numeric(length(vars))
  D[1L] <- 1
  if (length(vars) > 1L) for (k in 2:length(vars)) D[k] <- D[k - 1L] * (degs[k - 1L] + 1)
  total <- sum(degs * D)
  if (total > 400) return(list(f))
  # univariate image
  ue <- as.vector(f$e[, vars, drop = FALSE] %*% D)
  t1 <- vars[1L]
  em <- matrix(0L, length(ue), n)
  em[, t1] <- as.integer(ue)
  F <- qp_norm(em, f$c)
  zl <- .zp_factor_list(gmp::numerator(up_from_qpoly(F, t1)))
  ufs <- list()
  for (pr in zl$factors) for (r in seq_len(pr[[2L]])) ufs[[length(ufs) + 1L]] <- pr[[1L]]
  if (length(ufs) <= 1L) return(list(f))
  if (length(ufs) > 16L) return(list(f))
  decode <- function(u) {
    # bigz univariate -> multivariate qpoly via mixed radix digits
    d <- length(u) - 1L
    nzk <- which(!(u == 0))
    e <- matrix(0L, length(nzk), n)
    cs <- u[nzk]
    for (r in seq_along(nzk)) {
      pw <- d - nzk[r] + 1
      for (k in length(vars):1) {
        dk <- pw %/% D[k]
        if (dk > degs[k]) return(NULL)
        e[r, vars[k]] <- as.integer(dk)
        pw <- pw - dk * D[k]
      }
    }
    qp_norm(e, gmp::as.bigq(cs))
  }
  factors <- list()
  remaining <- seq_along(ufs)
  cur <- f
  s <- 1L
  while (2L * s <= length(remaining)) {
    found <- FALSE
    for (Sidx in .combn_list(length(remaining), s)) {
      S <- remaining[Sidx]
      G <- gmp::as.bigz(1)
      for (j in S) G <- .zp_mul(G, ufs[[j]])
      for (sg in c(1, -1)) {
        cand <- decode(G * sg)
        if (is.null(cand) || qp_is_const(cand)) next
        q <- qp_exquo(cur, cand)
        if (!is.null(q)) {
          factors[[length(factors) + 1L]] <- qp_primitive(cand)$poly
          cur <- q
          remaining <- setdiff(remaining, S)
          found <- TRUE
          break
        }
      }
      if (found) break
    }
    if (!found) s <- s + 1L
  }
  if (!qp_is_const(cur)) factors[[length(factors) + 1L]] <- qp_primitive(cur)$poly
  factors
}

# Factor a qpoly over QQ: list(coef = bigq, factors = list(list(qpoly, mult))).
qp_factor_list <- function(f) {
  n <- qp_nvars(f)
  if (qp_is_zero(f)) return(list(coef = gmp::as.bigq(0), factors = list()))
  pr <- qp_primitive(f)
  coef <- pr$content
  f <- pr$poly
  if (qp_is_const(f)) return(list(coef = coef * qp_lc(f), factors = list()))
  out <- list()
  # monomial content
  mins <- apply(f$e, 2L, min)
  if (any(mins > 0L)) {
    for (j in which(mins > 0L)) out[[length(out) + 1L]] <- list(qp_var(j, n), mins[j])
    f <- qp(f$e - matrix(mins, nrow(f$e), n, byrow = TRUE), f$c)
    f <- qp_norm(f$e, f$c)
  }
  pending <- list(f)
  while (length(pending)) {
    g <- pending[[1L]]
    pending <- pending[-1L]
    if (qp_is_const(g)) next
    vars <- qp_vars_used(g)
    split_done <- FALSE
    if (length(vars) > 1L) {
      for (v in vars) {
        ct <- qp_content_in_var(g, v)
        if (!qp_is_const(ct)) {
          pending <- c(pending, list(ct, qp_exquo(g, ct)))
          split_done <- TRUE
          break
        }
      }
    }
    if (split_done) next
    v <- vars[1L]
    sq <- .qp_sqf_list_var(g, v)
    for (pr2 in sq) {
      h <- pr2[[1L]]
      fs <- if (qp_total_degree(h) <= 1) list(h)
            else if (length(qp_vars_used(h)) > 1L) .wang_factor(qp_primitive(h)$poly, v)
            else .kronecker_factor(h)
      for (fh in fs) out[[length(out) + 1L]] <- list(qp_primitive(fh)$poly, pr2[[2L]])
    }
  }
  # merge repeated factors
  keys <- vapply(out, function(pr2) paste(.mono_keys(pr2[[1L]]$e), as.character(pr2[[1L]]$c), collapse = ";"), "")
  merged <- list()
  for (u in unique(keys)) {
    w <- which(keys == u)
    merged[[length(merged) + 1L]] <- list(out[[w[1L]]][[1L]], sum(vapply(out[w], `[[`, 0, 2L)))
  }
  prod_all <- qp_const(1, n)
  for (pr2 in merged) prod_all <- qp_mul(prod_all, qp_pow(pr2[[1L]], pr2[[2L]]))
  coef <- coef * (qp_lc(pr$poly) / qp_lc(prod_all))
  list(coef = coef, factors = merged)
}
