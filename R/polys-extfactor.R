# Factorization over algebraic number fields and the Gaussian domains, after
# sympy/polys/factortools.py (dup_ext_factor, dmp_ext_factor,
# dmp_trial_division, dmp_zz_i_factor, dmp_qq_i_factor, dmp_factor_list),
# sympy/polys/sqfreetools.py (dmp_sqf_norm, dmp_norm, dmp_sqf_part,
# dmp_sqf_p, _dmp_check_degrees) and the gcd SymPy runs over QQ<a> in
# sympy/polys/euclidtools.py (dmp_inner_gcd with dmp_multi_deflate, then
# dmp_ff_prs_gcd on the subresultant PRS).
#
# Dense recursive representation of densebasic.py with coefficients in
# K = Q(theta) (K from .nf_field()): an element is a bigq vector reduced
# modulo the minimal polynomial of theta (highest power first, zero is the
# empty vector, as the ANP rep); a level 0 polynomial is a list of elements
# (the .nfp_* representation), a level u polynomial a list of level u - 1
# polynomials; leading coefficient first, zero is list().
#
# The norm Res_a(m(a), g(a, X)) of dmp_sqf_norm is the determinant of the
# multiplication by g on the power basis of K (Bareiss over Q[X]); with m
# monic both are prod g(a_i, X) exactly, so r is SymPy's r. The norm is
# factored over QQ by qp_factor_list (Wang EEZ in R/polys-wang.R).

# ---------------------------------------------------------------------------
# elements of K

# K of .nf_field() with the reductions of theta^d, ..., theta^(2d - 2) as
# vectors of length d (highest power first), for products without division
.kd_field <- function(exts) {
  K <- .nf_field(exts)
  d <- K$d
  K$red <- lapply(seq_len(max(d - 1L, 0L)), function(j) {
    r <- .nfe_red(gmp::as.bigq(c(1L, integer(d + j - 1L))), K)
    c(gmp::as.bigq(integer(d - length(r))), r)
  })
  K
}
.ke_add <- function(a, b) if (!length(a)) b else if (!length(b)) a else up_add(a, b)
.ke_mul <- function(a, b, K) {
  la <- length(a); lb <- length(b)
  if (!la || !lb) return(.bq0())
  if (la == 1L) return(b * a)
  if (lb == 1L) return(a * b)
  p <- up_mul(a, b)
  n <- length(p); d <- K$d
  if (n <= d) return(p)
  # fold theta^k, k >= d, back with the table of .kd_field
  low <- p[(n - d + 1L):n]
  for (i in seq_len(n - d)) if (p[i] != 0) low <- low + p[i] * K$red[[n - i - d + 1L]]
  up_trim(low)
}
# inverse: 1/a for rationals, the conjugate over the norm in quadratic
# fields, else extended Euclid modulo the minimal polynomial
.ke_inv <- function(a, K) {
  if (length(a) == 1L) return(1 / a)
  if (K$d == 2L) {
    m1 <- K$mod[2L]; m0 <- K$mod[3L]
    a1 <- a[1L]; a0 <- a[2L]
    N <- a0 * a0 - m1 * a0 * a1 + m0 * a1 * a1
    return(up_trim(c(-a1, a0 - m1 * a1) / N))
  }
  .nfe_inv(a, K)
}
.ke_quo <- function(a, b, K) if (!length(a)) a else .ke_mul(a, .ke_inv(b, K), K)
.ke_pow <- function(a, n, K) { r <- gmp::as.bigq(1L); for (i in seq_len(n)) r <- .ke_mul(r, a, K); r }
.ke_is_one <- function(a) length(a) == 1L && a == 1
# AlgebraicField.is_negative: the sign of the leading coefficient of the rep
.ke_is_negative <- function(a) length(a) > 0L && a[1L] < 0
# multiples of the generator: s*K.unit
.ke_unit <- function(s) if (s == 0) .bq0() else gmp::as.bigq(c(s, 0L))

# ---------------------------------------------------------------------------
# densebasic and densearith over K

.kd_czero <- function(u) if (u == 0L) .bq0() else list()
.kd_one <- function(u) if (u == 0L) list(gmp::as.bigq(1L)) else list(.kd_one(u - 1L))
.kd_one_p <- function(f, u) {
  if (length(f) != 1L) return(FALSE)
  if (u == 0L) .ke_is_one(f[[1L]]) else .kd_one_p(f[[1L]], u - 1L)
}
.kd_ground <- function(c, u) {
  if (!length(c)) return(list())
  if (u == 0L) list(c) else list(.kd_ground(c, u - 1L))
}
.kd_strip <- function(f, u) {
  n <- length(f); i <- 1L
  while (i <= n && !length(f[[i]])) i <- i + 1L
  if (i == 1L) f else if (i > n) list() else f[i:n]
}
.kd_degree <- function(f) if (!length(f)) -Inf else length(f) - 1L
.kd_degree_in <- function(f, j, u) {
  if (j == 0L) return(.kd_degree(f))
  if (!length(f)) return(-Inf)
  max(vapply(f, .kd_degree_in, 0, j = j - 1L, u = u - 1L))
}
.kd_degree_list <- function(f, u) vapply(0:u, function(j) .kd_degree_in(f, j, u), 0)
.kd_ground_LC <- function(f, u) {
  while (u > 0L) { if (!length(f)) return(.bq0()); f <- f[[1L]]; u <- u - 1L }
  if (length(f)) f[[1L]] else .bq0()
}

.kd_add <- function(f, g, u) {
  lf <- length(f); lg <- length(g)
  if (!lf) return(g)
  if (!lg) return(f)
  if (lf < lg) { tmp <- f; f <- g; g <- tmp; lf <- length(f); lg <- length(g) }
  off <- lf - lg
  for (i in seq_len(lg)) f[[off + i]] <- if (u == 0L) .ke_add(f[[off + i]], g[[i]]) else .kd_add(f[[off + i]], g[[i]], u - 1L)
  if (off == 0L) .kd_strip(f, u) else f
}
.kd_neg <- function(f, u) if (u == 0L) lapply(f, function(a) -a) else lapply(f, .kd_neg, u = u - 1L)
.kd_sub <- function(f, g, u) .kd_add(f, .kd_neg(g, u), u)
.kd_mul <- function(f, g, u, K) {
  lf <- length(f); lg <- length(g)
  if (!lf || !lg) return(list())
  out <- rep(list(.kd_czero(u)), lf + lg - 1L)
  for (i in seq_len(lf)) {
    fi <- f[[i]]
    if (!length(fi)) next
    for (j in seq_len(lg)) {
      gj <- g[[j]]
      if (!length(gj)) next
      k <- i + j - 1L
      out[[k]] <- if (u == 0L) .ke_add(out[[k]], .ke_mul(fi, gj, K)) else .kd_add(out[[k]], .kd_mul(fi, gj, u - 1L, K), u - 1L)
    }
  }
  .kd_strip(out, u)
}
.kd_mul_ground <- function(f, c, u, K) {
  if (!length(c)) return(list())
  if (u == 0L) lapply(f, .ke_mul, b = c, K = K) else lapply(f, .kd_mul_ground, c = c, u = u - 1L, K = K)
}
# dmp_mul_term: f*c*x_0^j with c a coefficient (an element when u = 0)
.kd_mul_term <- function(f, c, j, u, K) {
  if (!length(f) || !length(c)) return(list())
  h <- if (u == 0L) lapply(f, .ke_mul, b = c, K = K) else lapply(f, .kd_mul, g = c, u = u - 1L, K = K)
  if (j > 0L) c(h, rep(list(.kd_czero(u)), j)) else h
}
# dmp_add_term: f + c*x_0^i
.kd_add_term <- function(f, c, i, u) {
  if (!length(c)) return(f)
  n <- length(f)
  if (i >= n) return(c(list(c), rep(list(.kd_czero(u)), i - n), f))
  m <- n - i
  f[[m]] <- if (u == 0L) .ke_add(f[[m]], c) else .kd_add(f[[m]], c, u - 1L)
  if (m == 1L) .kd_strip(f, u) else f
}
.kd_pow <- function(f, n, u, K) { r <- .kd_one(u); for (i in seq_len(n)) r <- .kd_mul(r, f, u, K); r }
.kd_ground_monic <- function(f, u, K) {
  if (!length(f)) return(f)
  lc <- .kd_ground_LC(f, u)
  if (.ke_is_one(lc)) f else .kd_mul_ground(f, .ke_inv(lc, K), u, K)
}

# dmp_ff_div: division with remainder over K
.kd_div <- function(f, g, u, K) {
  dg <- .kd_degree(g)
  if (dg < 0) stop("polynomial division by zero", call. = FALSE)
  q <- list(); r <- f; dr <- .kd_degree(f)
  if (dr < dg) return(list(q = q, r = r))
  lc_g <- g[[1L]]
  if (u == 0L) inv <- .ke_inv(lc_g, K)
  repeat {
    if (u == 0L) {
      c <- .ke_mul(r[[1L]], inv, K)
    } else {
      cr <- .kd_div(r[[1L]], lc_g, u - 1L, K)
      if (length(cr$r)) break
      c <- cr$q
    }
    j <- dr - dg
    q <- .kd_add_term(q, c, j, u)
    r <- .kd_sub(r, .kd_mul_term(g, c, j, u, K), u)
    odr <- dr; dr <- .kd_degree(r)
    if (dr < dg) break
    if (!(dr < odr)) stop("polynomial division failed", call. = FALSE)
  }
  list(q = q, r = r)
}
.kd_quo <- function(f, g, u, K) .kd_div(f, g, u, K)$q

# dmp_prem
.kd_prem <- function(f, g, u, K) {
  dg <- .kd_degree(g)
  if (dg < 0) stop("polynomial division by zero", call. = FALSE)
  r <- f; dr <- .kd_degree(f)
  if (dr < dg) return(r)
  N <- dr - dg + 1L
  lc_g <- g[[1L]]
  repeat {
    lc_r <- r[[1L]]
    j <- dr - dg; N <- N - 1L
    r <- .kd_sub(.kd_mul_term(r, lc_g, 0L, u, K), .kd_mul_term(g, lc_r, j, u, K), u)
    odr <- dr; dr <- .kd_degree(r)
    if (dr < dg) break
    if (!(dr < odr)) stop("polynomial division failed", call. = FALSE)
  }
  cc <- if (u == 0L) .ke_pow(lc_g, N, K) else .kd_pow(lc_g, N, u - 1L, K)
  .kd_mul_term(r, cc, 0L, u, K)
}

# arithmetic of the coefficients of a level u polynomial
.kd_cops <- function(u, K) {
  if (u == 0L) return(list(one = gmp::as.bigq(1L), neg = function(a) -a,
                           pow = function(a, n) .ke_pow(a, n, K), mul = function(a, b) .ke_mul(a, b, K),
                           quo = function(a, b) .ke_quo(a, b, K)))
  v <- u - 1L
  list(one = .kd_one(v), neg = function(a) .kd_neg(a, v), pow = function(a, n) .kd_pow(a, n, v, K),
       mul = function(a, b) .kd_mul(a, b, v, K), quo = function(a, b) .kd_quo(a, b, v, K))
}

# dup_inner_subresultants and dmp_inner_subresultants
.kd_inner_subresultants <- function(f, g, u, K) {
  op <- .kd_cops(u, K)
  n <- .kd_degree(f); m <- .kd_degree(g)
  if (n < m) { tmp <- f; f <- g; g <- tmp; tmp <- n; n <- m; m <- tmp }
  if (!length(f)) return(list(R = list(), S = list()))
  if (!length(g)) return(list(R = list(f), S = list(op$one)))
  R <- list(f, g); d <- n - m
  b <- op$pow(op$neg(op$one), d + 1L)
  h <- .kd_mul_term(.kd_prem(f, g, u, K), b, 0L, u, K)
  lc <- g[[1L]]
  cc <- op$pow(lc, d)
  S <- list(op$one, cc)
  cc <- op$neg(cc)
  while (length(h)) {
    k <- .kd_degree(h)
    R[[length(R) + 1L]] <- h
    f <- g; g <- h; d <- m - k; m <- k
    b <- op$mul(op$neg(lc), op$pow(cc, d))
    h <- lapply(.kd_prem(f, g, u, K), op$quo, b = b)
    lc <- g[[1L]]
    cc <- if (d > 1) op$quo(op$pow(op$neg(lc), d), op$pow(cc, d - 1L)) else op$neg(lc)
    S[[length(S) + 1L]] <- op$neg(cc)
  }
  list(R = R, S = S)
}

# dmp_content and dmp_primitive (u >= 1)
.kd_content <- function(f, u, K) {
  v <- u - 1L
  if (!length(f)) return(list())
  cont <- f[[1L]]
  for (c in f[-1L]) {
    cont <- .kd_gcd(cont, c, v, K)
    if (.kd_one_p(cont, v)) break
  }
  if (.ke_is_negative(.kd_ground_LC(cont, v))) .kd_neg(cont, v) else cont
}
.kd_primitive <- function(f, u, K) {
  cont <- .kd_content(f, u, K); v <- u - 1L
  if (!length(f) || .kd_one_p(cont, v)) return(list(cont = cont, f = f))
  list(cont = cont, f = lapply(f, .kd_quo, g = cont, u = v, K = K))
}

# dup_ff_prs_gcd and dmp_ff_prs_gcd with _dmp_ff_trivial_gcd and
# _dmp_simplify_gcd (USE_SIMPLIFY_GCD)
.kd_ff_prs_gcd <- function(f, g, u, K) {
  zf <- !length(f); zg <- !length(g)
  if (zf && zg) return(list(h = list(), cff = list(), cfg = list()))
  if (zf) return(list(h = .kd_ground_monic(g, u, K), cff = list(), cfg = .kd_ground(.kd_ground_LC(g, u), u)))
  if (zg) return(list(h = .kd_ground_monic(f, u, K), cff = .kd_ground(.kd_ground_LC(f, u), u), cfg = list()))
  if (u == 0L) {
    R <- .kd_inner_subresultants(f, g, 0L, K)$R
    h <- .kd_ground_monic(R[[length(R)]], 0L, K)
    return(list(h = h, cff = .kd_quo(f, h, 0L, K), cfg = .kd_quo(g, h, 0L, K)))
  }
  df <- .kd_degree(f); dg <- .kd_degree(g); v <- u - 1L
  if (!(df > 0 && dg > 0)) {
    if (df == 0 && dg == 0) { F <- f[[1L]]; G <- g[[1L]] }
    else if (df == 0) { F <- f[[1L]]; G <- .kd_content(g, u, K) }
    else { F <- .kd_content(f, u, K); G <- g[[1L]] }
    h <- .kd_gcd(F, G, v, K)
    return(list(h = list(h), cff = lapply(f, .kd_quo, g = h, u = v, K = K), cfg = lapply(g, .kd_quo, g = h, u = v, K = K)))
  }
  pf <- .kd_primitive(f, u, K); pg <- .kd_primitive(g, u, K)
  R <- .kd_inner_subresultants(pf$f, pg$f, u, K)$R
  cc <- .kd_ff_prs_gcd(pf$cont, pg$cont, v, K)$h
  h <- .kd_primitive(R[[length(R)]], u, K)$f
  h <- .kd_ground_monic(.kd_mul_term(h, cc, 0L, u, K), u, K)
  list(h = h, cff = .kd_quo(f, h, u, K), cfg = .kd_quo(g, h, u, K))
}

# dmp_inner_gcd: dmp_multi_deflate, the gcd, dmp_inflate
.kd_inner_gcd <- function(f, g, u, K) {
  if (u == 0L) return(.kd_ff_prs_gcd(f, g, 0L, K))
  tf <- .kd_terms(f, u); tg <- .kd_terms(g, u)
  B <- integer(u + 1L)
  for (t in c(tf, tg)) for (i in seq_len(u + 1L)) B[i] <- .igcd2(B[i], t[[1L]][i])
  B[B == 0L] <- 1L
  if (all(B == 1L)) return(.kd_ff_prs_gcd(f, g, u, K))
  defl <- function(tl) .kd_from_terms(lapply(tl, function(t) list(t[[1L]] %/% B, t[[2L]])), u)
  infl <- function(h) .kd_from_terms(lapply(.kd_terms(h, u), function(t) list(t[[1L]] * B, t[[2L]])), u)
  r <- .kd_ff_prs_gcd(defl(tf), defl(tg), u, K)
  lapply(r, infl)
}
.kd_gcd <- function(f, g, u, K) .kd_inner_gcd(f, g, u, K)$h

# terms (exponent vector, element) and back (dmp_to_dict, dmp_from_dict)
.kd_terms <- function(f, u) {
  out <- list()
  rec <- function(h, k, pre) {
    d <- length(h) - 1L
    for (i in seq_along(h)) {
      c <- h[[i]]
      if (!length(c)) next
      e <- c(pre, d - i + 1L)
      if (k == u) out[[length(out) + 1L]] <<- list(e, c) else rec(c, k + 1L, e)
    }
  }
  rec(f, 0L, integer(0))
  out
}
.kd_from_terms <- function(terms, u) {
  if (!length(terms)) return(list())
  e1 <- vapply(terms, function(t) as.integer(t[[1L]][1L]), 0L)
  D <- max(e1)
  if (u == 0L) {
    out <- rep(list(.bq0()), D + 1L)
    for (t in terms) { i <- D - t[[1L]][1L] + 1L; out[[i]] <- .ke_add(out[[i]], t[[2L]]) }
    return(.kd_strip(out, 0L))
  }
  out <- lapply(D:0, function(d) .kd_from_terms(lapply(terms[e1 == d], function(t) list(t[[1L]][-1L], t[[2L]])), u - 1L))
  .kd_strip(out, u)
}
.kd_from_qp <- function(p, u) {
  if (qp_is_zero(p)) return(list())
  .kd_from_terms(lapply(seq_along(p$c), function(i) list(p$e[i, ], p$c[i])), u)
}

# dmp_diff_in and dmp_shift
.kd_mul_int <- function(f, k, u) if (u == 0L) lapply(f, function(a) if (length(a)) a * k else a) else lapply(f, .kd_mul_int, k = k, u = u - 1L)
.kd_diff <- function(f, u) {
  n <- length(f) - 1L
  if (n < 1L) return(list())
  .kd_strip(lapply(seq_len(n), function(i) .kd_mul_int(list(f[[i]]), n - i + 1L, u)[[1L]]), u)
}
.kd_diff_in <- function(f, j, u) {
  if (j == 0L) return(.kd_diff(f, u))
  .kd_strip(lapply(f, .kd_diff_in, j = j - 1L, u = u - 1L), u)
}
.kd_shift <- function(f, a, u, K) {
  if (!length(f)) return(f)
  a0 <- a[[1L]]
  if (u > 0L && any(vapply(a[-1L], length, 0L) > 0L)) f <- lapply(f, .kd_shift, a = a[-1L], u = u - 1L, K = K)
  n <- length(f) - 1L
  if (length(a0) && n > 0L) {
    for (i in n:1) for (j in 0:(i - 1L)) {
      f[[j + 2L]] <- if (u == 0L) .ke_add(f[[j + 2L]], .ke_mul(f[[j + 1L]], a0, K))
                     else .kd_add(f[[j + 2L]], .kd_mul_ground(f[[j + 1L]], a0, u - 1L, K), u - 1L)
    }
  }
  .kd_strip(f, u)
}

# dmp_terms_gcd
.kd_terms_gcd <- function(f, u) {
  tl <- .kd_terms(f, u)
  if (!length(tl)) return(list(J = integer(u + 1L), f = f))
  G <- Reduce(pmin, lapply(tl, `[[`, 1L))
  if (all(G == 0L)) return(list(J = G, f = f))
  list(J = G, f = .kd_from_terms(lapply(tl, function(t) list(t[[1L]] - G, t[[2L]])), u))
}

# ---------------------------------------------------------------------------
# sqfreetools

# dup_sqf_part and dmp_sqf_part over K
.kd_sqf_part <- function(f, u, K) {
  if (!length(f)) return(f)
  if (.ke_is_negative(.kd_ground_LC(f, u))) f <- .kd_neg(f, u)
  g <- if (u == 0L) .kd_gcd(f, .kd_diff(f, 0L), 0L, K) else {
    gg <- f
    for (i in 0:u) gg <- .kd_gcd(gg, .kd_diff_in(f, i, u), u, K)
    gg
  }
  .kd_ground_monic(.kd_quo(f, g, u, K), u, K)
}

# dmp_norm: Res_a(m(a), f(a, X)) as a qpoly in the u + 1 generators
.kd_norm <- function(f, u, K) {
  d <- K$d; n <- u + 1L
  tl <- .kd_terms(f, u)
  if (!length(tl)) return(qp_zero(n))
  E <- do.call(rbind, lapply(tl, `[[`, 1L))
  coord <- function(a, i) { L <- length(a); if (i < L) a[L - i] else gmp::as.bigq(0L) }
  G <- lapply(0:(d - 1L), function(i) qp_norm(E, do.call(c, lapply(tl, function(t) coord(t[[2L]], i)))))
  P <- lapply(0:(2L * d - 2L), function(k) .nfe_red(gmp::as.bigq(c(1L, integer(k))), K))
  M <- matrix(list(), d, d)
  for (r in 0:(d - 1L)) for (cc in 0:(d - 1L)) {
    acc <- qp_zero(n)
    for (i in 0:(d - 1L)) {
      w <- coord(P[[i + cc + 1L]], r)
      if (w != 0 && !qp_is_zero(G[[i + 1L]])) acc <- qp_add(acc, qp_scale(G[[i + 1L]], w))
    }
    M[[r + 1L, cc + 1L]] <- acc
  }
  .qp_det(M, n)
}
# fraction-free (Bareiss) determinant of a matrix of qpolys
.qp_det <- function(M, n) {
  d <- nrow(M)
  if (d == 1L) return(M[[1L, 1L]])
  sgn <- 1L; prev <- NULL
  for (k in seq_len(d - 1L)) {
    if (qp_is_zero(M[[k, k]])) {
      piv <- which(vapply((k + 1L):d, function(r) !qp_is_zero(M[[r, k]]), TRUE))
      if (!length(piv)) return(qp_zero(n))
      r <- k + piv[1L]; tmp <- M[k, ]; M[k, ] <- M[r, ]; M[r, ] <- tmp; sgn <- -sgn
    }
    for (i in (k + 1L):d) for (j in (k + 1L):d) {
      v <- qp_sub(qp_mul(M[[k, k]], M[[i, j]]), qp_mul(M[[i, k]], M[[k, j]]))
      M[[i, j]] <- if (is.null(prev)) v else qp_exquo(v, prev)
    }
    prev <- M[[k, k]]
  }
  if (sgn < 0L) qp_neg(M[[d, d]]) else M[[d, d]]
}
# dmp_sqf_p over QQ
.qp_sqf_p <- function(r) {
  if (qp_is_zero(r)) return(TRUE)
  for (i in seq_len(qp_nvars(r))) {
    fp <- qp_diff(r, i)
    if (qp_is_zero(fp)) next
    if (qp_degree(qp_gcd(r, fp), i) != 0) return(FALSE)
  }
  TRUE
}
# dup_sqf_norm and dmp_sqf_norm with the shifts of _dmp_sqf_norm_shifts:
# list(s, g = f(X - s*a), r = Norm(g) as a qpoly)
.kd_sqf_norm <- function(f, u, K) {
  if (u == 0L) {
    s <- 0L
    repeat {
      r <- .kd_norm(f, 0L, K)
      if (.qp_sqf_p(r)) return(list(s = s, g = f, r = r))
      f <- .kd_shift(f, list(.ke_unit(-1L)), 0L, K); s <- s + 1L
    }
  }
  n <- u + 1L
  try_s <- function(sv) {
    g <- if (all(sv == 0L)) f else .kd_shift(f, lapply(sv, function(si) .ke_unit(-si)), u, K)
    r <- .kd_norm(g, u, K)
    if (.qp_sqf_p(r)) list(s = sv, g = g, r = r) else NULL
  }
  res <- try_s(integer(n))
  if (!is.null(res)) return(res)
  dl <- .kd_degree_list(f, u)
  for (i in order(dl)) if (dl[i] > 0) {
    sv <- integer(n); sv[i] <- 1L
    res <- try_s(sv)
    if (!is.null(res)) return(res)
  }
  j <- 0L
  repeat { j <- j + 1L; res <- try_s(rep(j, n)); if (!is.null(res)) return(res) }
}

# ---------------------------------------------------------------------------
# factortools

# _sort_factors: (len(f), multiplicity, rep); the coefficients compare as
# ANP reps over QQ<a> and as [y, x] over ZZ_I and QQ_I
.kd_cmp_el <- function(a, b, gauss = FALSE) {
  if (gauss) { a <- .kd_yx(a); b <- .kd_yx(b) }
  for (i in seq_len(min(length(a), length(b)))) if (a[i] != b[i]) return(if (a[i] < b[i]) -1L else 1L)
  as.integer(sign(length(a) - length(b)))
}
.kd_yx <- function(a) switch(length(a) + 1L, gmp::as.bigq(c(0L, 0L)), c(gmp::as.bigq(0L), a), a)
.kd_cmp_poly <- function(f, g, u, gauss = FALSE) {
  for (i in seq_len(min(length(f), length(g)))) {
    r <- if (u == 0L) .kd_cmp_el(f[[i]], g[[i]], gauss) else .kd_cmp_poly(f[[i]], g[[i]], u - 1L, gauss)
    if (r != 0L) return(r)
  }
  as.integer(sign(length(f) - length(g)))
}
.kd_sort_factors <- function(fs, u, gauss = FALSE) {
  key <- function(a, b) {
    if (length(a[[1L]]) != length(b[[1L]])) return(length(a[[1L]]) - length(b[[1L]]))
    if (a[[2L]] != b[[2L]]) return(a[[2L]] - b[[2L]])
    .kd_cmp_poly(a[[1L]], b[[1L]], u, gauss)
  }
  out <- list()
  for (f in fs) {
    i <- length(out) + 1L
    while (i > 1L && key(out[[i - 1L]], f) > 0) i <- i - 1L
    out <- append(out, list(f), i - 1L)
  }
  out
}

# dmp_trial_division
.kd_trial_division <- function(f, factors, u, K) {
  result <- list()
  for (fac in factors) {
    k <- 0L
    repeat {
      d <- .kd_div(f, fac, u, K)
      if (length(d$r)) break
      f <- d$q; k <- k + 1L
    }
    if (k == 0L) stop("trial division failed", call. = FALSE)
    result[[length(result) + 1L]] <- list(fac, k)
  }
  .kd_sort_factors(result, u)
}
# _dup_check_degrees and _dmp_check_degrees
.kd_check_degrees <- function(f, u, result) {
  degs <- numeric(u + 1L)
  for (fk in result) degs <- degs + fk[[2L]] * .kd_degree_list(fk[[1L]], u)
  if (!isTRUE(all(degs == .kd_degree_list(f, u)))) stop("assertion failed: degrees of the computed factorization", call. = FALSE)
  invisible(TRUE)
}

# dup_ext_factor and dmp_ext_factor: list(lc, factors = list(list(h, k)))
.kd_ext_factor <- function(f, u, K) {
  if (u == 0L) {
    n <- .kd_degree(f); lc <- .kd_ground_LC(f, 0L)
    f <- .kd_ground_monic(f, 0L, K)
    if (n <= 0) return(list(lc = lc, factors = list()))
    if (n == 1) return(list(lc = lc, factors = list(list(f, 1L))))
  } else {
    lc <- .kd_ground_LC(f, u)
    f <- .kd_ground_monic(f, u, K)
    if (all(.kd_degree_list(f, u) <= 0)) return(list(lc = lc, factors = list()))
  }
  F <- f
  f <- .kd_sqf_part(f, u, K)
  sn <- .kd_sqf_norm(f, u, K)
  fl <- qp_factor_list(sn$r)$factors
  if (length(fl) <= 1L) {
    if (u == 0L) return(list(lc = lc, factors = list(list(f, as.integer(n %/% .kd_degree(f))))))
    factors <- list(f)
  } else {
    g <- sn$g
    a <- lapply(sn$s, .ke_unit)
    factors <- lapply(fl, function(fk) {
      ig <- .kd_inner_gcd(.kd_from_qp(fk[[1L]], u), g, u, K)
      g <<- ig$cfg
      .kd_shift(ig$h, a, u, K)
    })
  }
  result <- .kd_trial_division(F, factors, u, K)
  .kd_check_degrees(F, u, result)
  list(lc = lc, factors = result)
}

# Gaussian integers of an element (x + y*I as c(x, y)) and back
.kd_el_gi <- function(a) switch(length(a) + 1L, gmp::as.bigq(c(0L, 0L)), c(a, gmp::as.bigq(0L)), c(a[2L], a[1L]))
.kd_gi_el <- function(p) if (p[2L] != 0) gmp::as.bigq(c(p[2L], p[1L])) else if (p[1L] != 0) gmp::as.bigq(p[1L]) else .bq0()
.kd_els <- function(f, u) if (u == 0L) f else do.call(c, c(list(list()), lapply(f, .kd_els, u = u - 1L)))
# dmp_ground_primitive over ZZ_I: the gcd of the coefficients normalized by
# canonical_unit (dup_content stops at one)
.kd_zzi_primitive <- function(f, u, K) {
  cont <- gmp::as.bigq(c(0L, 0L))
  for (a in .kd_els(f, u)) {
    cont <- .gi_gcd(cont, .kd_el_gi(a))
    if (cont[1L] == 1 && cont[2L] == 0) break
  }
  ce <- .kd_gi_el(cont)
  if (!length(f) || .ke_is_one(ce)) return(list(cont = ce, f = f))
  list(cont = ce, f = .kd_mul_ground(f, .ke_inv(ce, K), u, K))
}
# dmp_zz_i_factor after the factorization over QQ<I>: each monic factor is
# cleared of denominators and made ZZ_I primitive, the leading coefficient
# takes the contents
.kd_zz_i_post <- function(coeff, factors, u, K) {
  out <- list()
  for (fk in factors) {
    fac <- fk[[1L]]; i <- fk[[2L]]
    den <- gmp::as.bigz(1L)
    for (a in .kd_els(fac, u)) { ds <- gmp::denominator(a); for (k in seq_along(ds)) den <- gmp::lcm.bigz(den, ds[k]) }
    num <- if (den == 1) fac else .kd_mul_ground(fac, gmp::as.bigq(den), u, K)
    pp <- .kd_zzi_primitive(num, u, K)
    coeff <- .ke_quo(.ke_mul(coeff, .ke_pow(pp$cont, i, K), K), gmp::as.bigq(den)^i, K)
    out[[length(out) + 1L]] <- list(pp$f, i)
  }
  list(coeff = coeff, factors = out)
}

# dmp_factor_list over QQ<a> (dom = "alg"), QQ_I or ZZ_I: list(coeff, factors)
.kd_factor_list <- function(f, u, K, dom = "alg") {
  tg <- .kd_terms_gcd(f, u); J <- tg$J; f <- tg$f
  # dmp_ground_primitive: the content is one over QQ<a> (Field.gcd without a
  # ring), cancels against the leading coefficient over QQ_I
  cont <- gmp::as.bigq(1L)
  if (dom == "ZZ_I") { pp <- .kd_zzi_primitive(f, u, K); cont <- pp$cont; f <- pp$f }
  r <- .kd_ext_factor(f, u, K)
  coeff <- r$lc; factors <- r$factors
  if (dom == "ZZ_I") { z <- .kd_zz_i_post(coeff, factors, u, K); coeff <- z$coeff; factors <- z$factors }
  for (i in seq_along(J) - 1L) {
    j <- J[u + 1L - i]
    if (!j) next
    e <- integer(u + 1L); e[u + 1L - i] <- 1L
    factors <- c(list(list(.kd_from_terms(list(list(e, gmp::as.bigq(1L))), u), as.integer(j))), factors)
  }
  list(coeff = .ke_mul(coeff, cont, K), factors = .kd_sort_factors(factors, u, dom != "alg"))
}

# ---------------------------------------------------------------------------
# expressions

.kd_from_expr <- function(e, gens, K) {
  e <- expand(as_sym(e))
  # generators that are not symbols (f(t), its derivatives) become symbols,
  # so that a symbol generator inside them does not hide the polynomial
  nonsym <- which(!vapply(gens, function(g) g$type %in% c("Symbol", "Dummy"), TRUE))
  if (length(nonsym)) {
    ds <- lapply(nonsym, function(i) Dummy())
    e <- xreplace(e, lapply(seq_along(nonsym), function(k) list(gens[[nonsym[k]]], ds[[k]])))
    gens[nonsym] <- ds
  }
  terms <- if (e$type == "Add") e$args else list(e)
  if (length(terms) == 1L && is_zero_num(terms[[1L]])) return(list())
  # the monomial in the generators of each term and its coefficient
  cs <- list(); ex <- list()
  for (t in terms) {
    ev <- integer(length(gens)); c <- t
    for (i in seq_along(gens)) {
      d <- .degree_in(c, gens[[i]])
      if (is.na(d)) stop("not a polynomial in ", .str(gens[[i]]), ": ", .str(e), call. = FALSE)
      if (d > 0L) { ev[i] <- d; c <- .div(c, .pow(gens[[i]], .int(d))) }
    }
    k <- paste(ev, collapse = ",")
    if (is.null(cs[[k]])) { cs[[k]] <- c; ex[[k]] <- ev } else cs[[k]] <- .add(list(cs[[k]], c))
  }
  .kd_from_terms(lapply(names(cs), function(k) list(ex[[k]], .nfe_from_expr(cs[[k]], K))), length(gens) - 1L)
}
# Poly.as_expr: the sum of coeff*x^i*y^j*...
.kd_to_expr <- function(f, gens, u, K) {
  tl <- .kd_terms(f, u)
  if (!length(tl)) return(S_Zero)
  .add(lapply(tl, function(t) {
    e <- t[[1L]]
    .mul(c(list(.nfe_to_expr(t[[2L]], K)), lapply(which(e > 0L), function(i) .pow(gens[[i]], .int(e[i])))))
  }))
}
.kd_gens <- function(e, gens) {
  gs <- .gens_of(list(e), gens)
  Filter(function(g) !(g$type == "Pow" && is_exact_number(g$args[[1L]])) && !is_const(g, "I"), gs)
}

# sqf_part over QQ<a> or QQ_I (dmp_sqf_part): the monic square-free part.
.field_sqf_part <- function(f, Fd, gens) {
  e <- expand(as_sym(f))
  K <- .kd_field(Fd$exts)
  gs <- .kd_gens(e, gens)
  if (!length(gs)) return(e)
  u <- length(gs) - 1L
  .kd_to_expr(.kd_sqf_part(.kd_from_expr(e, gs, K), u, K), gs, u, K)
}

# sqf_norm(f, extension): list(s, g, r) of Poly.sqf_norm, s a list of Integers.
.ext_sqf_norm <- function(f, extension) {
  f <- expand(as_sym(f))
  K <- .kd_field(if (is_sym(extension)) list(extension) else lapply(as.list(extension), as_sym))
  gs <- .poly_sort_gens(free_symbols(f))
  u <- length(gs) - 1L
  sn <- .kd_sqf_norm(.kd_from_expr(f, gs, K), u, K)
  list(lapply(sn$s, .int), .kd_to_expr(sn$g, gs, u, K), qpoly_to_expr(sn$r, gs))
}
