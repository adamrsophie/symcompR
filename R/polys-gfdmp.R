# Dense recursive polynomials over GF(p), after the generic algorithms SymPy
# runs on the domain FF(p) (DMP_Python): sympy/polys/densebasic.py,
# densearith.py, densetools.py and euclidtools.py. A level u polynomial is a
# list of level u - 1 polynomials, leading coefficient first, zero is
# list(); a level -1 element is a bigz scalar in [0, p). The modulus p is a
# bigz. Over a field SymPy's gcd is ground monic and unique, so the
# deflation and simplification shortcuts of dmp_inner_gcd are not needed to
# reproduce its results.

# ---------------------------------------------------------------------------
# elements and basic structure

.gd_z <- function(v) if (v < 0L) gmp::as.bigz(0L) else list()
.gd_isz <- function(a, v) if (v < 0L) a == 0 else !length(a)
.gd_one <- function(v) if (v < 0L) gmp::as.bigz(1L) else list(.gd_one(v - 1L))
.gd_ground <- function(c, v) {
  # the scalar c as a level v polynomial
  if (v < 0L) return(c)
  if (c == 0) return(list())
  list(.gd_ground(c, v - 1L))
}
.gd_strip <- function(f, u) {
  while (length(f) && .gd_isz(f[[1L]], u - 1L)) f <- f[-1L]
  f
}
.gd_deg <- function(f) if (!length(f)) -Inf else length(f) - 1L
.gd_LC <- function(f, u) if (length(f)) f[[1L]] else .gd_z(u - 1L)
.gd_ground_LC <- function(f, u) {
  if (u < 0L) return(f)
  if (!length(f)) return(gmp::as.bigz(0L))
  .gd_ground_LC(f[[1L]], u - 1L)
}
.gd_one_p <- function(f, u) {
  if (u < 0L) return(f == 1)
  length(f) == 1L && .gd_one_p(f[[1L]], u - 1L)
}

# ---------------------------------------------------------------------------
# arithmetic (p is threaded through every call)

.gd_add <- function(f, g, u, p) {
  if (u < 0L) return((f + g) %% p)
  n <- max(length(f), length(g))
  if (!n) return(list())
  z <- .gd_z(u - 1L)
  f <- c(rep(list(z), n - length(f)), f); g <- c(rep(list(z), n - length(g)), g)
  .gd_strip(Map(function(a, b) .gd_add(a, b, u - 1L, p), f, g), u)
}
.gd_neg <- function(f, u, p) {
  if (u < 0L) return((-f) %% p)
  lapply(f, .gd_neg, u = u - 1L, p = p)
}
.gd_sub <- function(f, g, u, p) .gd_add(f, .gd_neg(g, u, p), u, p)
.gd_mul <- function(f, g, u, p) {
  if (u < 0L) return((f * g) %% p)
  if (!length(f) || !length(g)) return(list())
  n <- length(f) + length(g) - 1L
  out <- rep(list(.gd_z(u - 1L)), n)
  for (i in seq_along(f)) {
    if (.gd_isz(f[[i]], u - 1L)) next
    for (j in seq_along(g)) {
      if (.gd_isz(g[[j]], u - 1L)) next
      out[[i + j - 1L]] <- .gd_add(out[[i + j - 1L]], .gd_mul(f[[i]], g[[j]], u - 1L, p), u - 1L, p)
    }
  }
  .gd_strip(out, u)
}
# f times a level u - 1 coefficient c, times x^j
.gd_mul_term <- function(f, c, j, u, p) {
  if (!length(f) || .gd_isz(c, u - 1L)) return(list())
  c(lapply(f, .gd_mul, g = c, u = u - 1L, p = p), rep(list(.gd_z(u - 1L)), j))
}
.gd_add_term <- function(f, c, j, u, p) .gd_add(f, .gd_mul_term(.gd_one(u), c, j, u, p), u, p)
.gd_mul_scalar <- function(f, c, u, p) {
  if (u < 0L) return((f * c) %% p)
  .gd_strip(lapply(f, .gd_mul_scalar, c = c, u = u - 1L, p = p), u)
}
.gd_pow <- function(f, n, u, p) {
  r <- .gd_one(u)
  for (i in seq_len(n)) r <- .gd_mul(r, f, u, p)
  r
}
.gd_inv <- function(a, p) {
  r <- suppressWarnings(gmp::inv.bigz(a %% p, p))
  if (is.na(r)) stop("NotInvertible: zero divisor", call. = FALSE)
  r
}
.gd_ground_monic <- function(f, u, p) {
  if (u < 0L) return(if (f == 0) f else gmp::as.bigz(1L))
  if (!length(f)) return(f)
  lc <- .gd_ground_LC(f, u)
  if (lc == 1) f else .gd_mul_scalar(f, .gd_inv(lc, p), u, p)
}
.gd_diff <- function(f, u, p) {
  n <- length(f) - 1L
  if (n <= 0L) return(list())
  .gd_strip(Map(function(a, k) .gd_mul_scalar(a, gmp::as.bigz(k) %% p, u - 1L, p), f[-length(f)], n:1), u)
}

# ---------------------------------------------------------------------------
# division

# dmp_ff_div (dup_ff_div for u = 0)
.gd_div <- function(f, g, u, p) {
  if (u < 0L) {
    if (g == 0) stop("ZeroDivisionError: polynomial division", call. = FALSE)
    return(list(q = (f * .gd_inv(g, p)) %% p, r = gmp::as.bigz(0L)))
  }
  dg <- .gd_deg(g)
  if (dg < 0) stop("ZeroDivisionError: polynomial division", call. = FALSE)
  q <- list(); r <- f; dr <- .gd_deg(f)
  if (dr < dg) return(list(q = q, r = r))
  lc_g <- .gd_LC(g, u)
  repeat {
    lc_r <- .gd_LC(r, u)
    cr <- .gd_div(lc_r, lc_g, u - 1L, p)
    if (!.gd_isz(cr$r, u - 1L)) break
    j <- dr - dg
    q <- .gd_add_term(q, cr$q, j, u, p)
    r <- .gd_sub(r, .gd_mul_term(g, cr$q, j, u, p), u, p)
    pdr <- dr; dr <- .gd_deg(r)
    if (dr < dg) break
    if (!(dr < pdr)) stop("PolynomialDivisionFailed", call. = FALSE)
  }
  list(q = q, r = r)
}
.gd_quo <- function(f, g, u, p) .gd_div(f, g, u, p)$q
.gd_rem <- function(f, g, u, p) .gd_div(f, g, u, p)$r
.gd_exquo <- function(f, g, u, p, fe = NULL, ge = NULL) {
  d <- .gd_div(f, g, u, p)
  if (!.gd_isz(d$r, u)) stop("ExactQuotientFailed: ", if (is.null(ge)) "g" else ge, " does not divide ",
                             if (is.null(fe)) "f" else fe, call. = FALSE)
  d$q
}

# dmp_pdiv and dmp_prem
.gd_pdiv <- function(f, g, u, p, rem_only = FALSE) {
  dg <- .gd_deg(g)
  if (dg < 0) stop("ZeroDivisionError: polynomial division", call. = FALSE)
  q <- list(); r <- f; dr <- .gd_deg(f)
  if (dr < dg) return(list(q = q, r = r))
  N <- dr - dg + 1
  lc_g <- .gd_LC(g, u)
  repeat {
    lc_r <- .gd_LC(r, u)
    j <- dr - dg; N <- N - 1
    if (!rem_only) q <- .gd_add_term(.gd_mul_term(q, lc_g, 0L, u, p), lc_r, j, u, p)
    r <- .gd_sub(.gd_mul_term(r, lc_g, 0L, u, p), .gd_mul_term(g, lc_r, j, u, p), u, p)
    pdr <- dr; dr <- .gd_deg(r)
    if (dr < dg) break
    if (!(dr < pdr)) stop("PolynomialDivisionFailed", call. = FALSE)
  }
  cN <- .gd_pow(lc_g, N, u - 1L, p)
  list(q = .gd_mul_term(q, cN, 0L, u, p), r = .gd_mul_term(r, cN, 0L, u, p))
}
.gd_prem <- function(f, g, u, p) .gd_pdiv(f, g, u, p, rem_only = TRUE)$r

# ---------------------------------------------------------------------------
# subresultants, resultants, discriminants

.gd_inner_subresultants <- function(f, g, u, p) {
  v <- u - 1L
  n <- .gd_deg(f); m <- .gd_deg(g)
  if (n < m) { t <- f; f <- g; g <- t; t <- n; n <- m; m <- t }
  if (!length(f)) return(list(R = list(), S = list()))
  if (!length(g)) return(list(R = list(f), S = list(.gd_one(v))))
  R <- list(f, g)
  d <- n - m
  b <- .gd_pow(.gd_ground(p - 1, v), d + 1, v, p)
  h <- .gd_mul_term(.gd_prem(f, g, u, p), b, 0L, u, p)
  lc <- .gd_LC(g, u)
  cc <- .gd_pow(lc, d, v, p)
  S <- list(.gd_one(v), cc)
  cc <- .gd_neg(cc, v, p)
  while (length(h)) {
    k <- .gd_deg(h)
    R[[length(R) + 1L]] <- h
    f <- g; g <- h; d <- m - k; m <- k
    b <- .gd_mul(.gd_neg(lc, v, p), .gd_pow(cc, d, v, p), v, p)
    h <- .gd_prem(f, g, u, p)
    h <- .gd_strip(lapply(h, function(ch) .gd_quo(ch, b, v, p)), u)
    lc <- .gd_LC(g, u)
    if (d > 1) cc <- .gd_quo(.gd_pow(.gd_neg(lc, v, p), d, v, p), .gd_pow(cc, d - 1, v, p), v, p)
    else cc <- .gd_neg(lc, v, p)
    S[[length(S) + 1L]] <- .gd_neg(cc, v, p)
  }
  list(R = R, S = S)
}
.gd_subresultants <- function(f, g, u, p) .gd_inner_subresultants(f, g, u, p)$R
.gd_resultant <- function(f, g, u, p) {
  if (!length(f) || !length(g)) return(.gd_z(u - 1L))
  rs <- .gd_inner_subresultants(f, g, u, p)
  if (.gd_deg(rs$R[[length(rs$R)]]) > 0) return(.gd_z(u - 1L))
  rs$S[[length(rs$S)]]
}
.gd_discriminant <- function(f, u, p) {
  d <- .gd_deg(f); v <- u - 1L
  if (d <= 0) return(.gd_z(v))
  s <- if (((d * (d - 1)) %/% 2) %% 2 == 0) gmp::as.bigz(1L) else p - 1
  cc <- .gd_mul_scalar(.gd_LC(f, u), s, v, p)
  .gd_quo(.gd_resultant(f, .gd_diff(f, u, p), u, p), cc, v, p)
}

# ---------------------------------------------------------------------------
# gcd, lcm, content

.gd_trivial_gcd <- function(f, g, u, p) {
  zf <- !length(f); zg <- !length(g)
  if (zf && zg) return(list(h = list(), cff = list(), cfg = list()))
  if (zf) return(list(h = .gd_ground_monic(g, u, p), cff = list(), cfg = .gd_ground(.gd_ground_LC(g, u), u)))
  if (zg) return(list(h = .gd_ground_monic(f, u, p), cff = .gd_ground(.gd_ground_LC(f, u), u), cfg = list()))
  NULL
}
# dmp_ff_prs_gcd: list(h, cff, cfg)
.gd_cofactors <- function(f, g, u, p) {
  if (u < 0L) return(list(h = gmp::as.bigz(1L), cff = f, cfg = g))
  tr <- .gd_trivial_gcd(f, g, u, p)
  if (!is.null(tr)) return(tr)
  if (u == 0L) {
    sr <- .gd_subresultants(f, g, u, p)
    h <- .gd_ground_monic(sr[[length(sr)]], u, p)
  } else {
    pf <- .gd_primitive(f, u, p); pg <- .gd_primitive(g, u, p)
    sr <- .gd_subresultants(pf$f, pg$f, u, p)
    cc <- .gd_cofactors(pf$cont, pg$cont, u - 1L, p)$h
    h <- .gd_primitive(sr[[length(sr)]], u, p)$f
    h <- .gd_ground_monic(.gd_mul_term(h, cc, 0L, u, p), u, p)
  }
  list(h = h, cff = .gd_quo(f, h, u, p), cfg = .gd_quo(g, h, u, p))
}
.gd_gcd <- function(f, g, u, p) .gd_cofactors(f, g, u, p)$h
.gd_lcm <- function(f, g, u, p) {
  if (!length(f) || !length(g)) return(list())
  .gd_ground_monic(.gd_quo(.gd_mul(f, g, u, p), .gd_gcd(f, g, u, p), u, p), u, p)
}
# dmp_content and dmp_primitive in the main variable (u >= 1)
.gd_primitive <- function(f, u, p) {
  v <- u - 1L
  cont <- .gd_LC(f, u)
  if (!length(f)) return(list(cont = cont, f = f))
  for (cc in f[-1L]) {
    cont <- .gd_gcd(cont, cc, v, p)
    if (.gd_one_p(cont, v)) break
  }
  if (.gd_one_p(cont, v)) return(list(cont = cont, f = f))
  list(cont = cont, f = lapply(f, .gd_quo, g = cont, u = v, p = p))
}

# ---------------------------------------------------------------------------
# univariate: extended Euclid, composition and decomposition

.gd_half_gcdex <- function(f, g, p) {
  a <- list(gmp::as.bigz(1L)); b <- list()
  while (length(g)) {
    d <- .gd_div(f, g, 0L, p)
    f <- g; g <- d$r
    t <- .gd_sub(a, .gd_mul(d$q, b, 0L, p), 0L, p); a <- b; b <- t
  }
  a <- .gd_mul_scalar(a, .gd_inv(.gd_LC(f, 0L), p), 0L, p)
  list(s = a, h = .gd_ground_monic(f, 0L, p))
}
.gd_gcdex <- function(f, g, p) {
  hg <- .gd_half_gcdex(f, g, p)
  Fp <- .gd_sub(hg$h, .gd_mul(hg$s, f, 0L, p), 0L, p)
  list(s = hg$s, t = .gd_quo(Fp, g, 0L, p), h = hg$h)
}
.gd_invert <- function(f, g, p) {
  hg <- .gd_half_gcdex(f, g, p)
  if (!.gd_one_p(hg$h, 0L)) stop("NotInvertible: zero divisor", call. = FALSE)
  .gd_rem(hg$s, g, 0L, p)
}
# dmp_compose: f(g) by Horner's rule in the main variable
.gd_compose <- function(f, g, u, p) {
  if (!length(f)) return(list())
  h <- list(f[[1L]])
  for (cc in f[-1L]) h <- .gd_add(.gd_mul(h, g, u, p), .gd_ground_poly(cc, u), u, p)
  h
}
# a level u - 1 coefficient as a constant level u polynomial
.gd_ground_poly <- function(c, u) if (.gd_isz(c, u - 1L)) list() else list(c)

# dup_decompose and its helpers _dup_right_decompose, _dup_left_decompose
.gd_decompose <- function(f, p) {
  F <- list()
  repeat {
    res <- .gd_decompose_step(f, p)
    if (is.null(res)) break
    f <- res$g
    F <- c(list(res$h), F)
  }
  c(list(f), F)
}
.gd_decompose_step <- function(f, p) {
  df <- length(f) - 1L
  if (df < 3L) return(NULL)
  for (s in 2:(df - 1L)) {
    if (df %% s != 0L) next
    h <- .gd_right_decompose(f, s, p)
    g <- .gd_left_decompose(f, h, p)
    if (!is.null(g)) return(list(g = g, h = h))
  }
  NULL
}
.gd_coef <- function(f, i) {
  # the coefficient of x^i of a univariate f
  n <- length(f) - 1L
  if (i > n || i < 0L) gmp::as.bigz(0L) else f[[n - i + 1L]]
}
# _dup_right_decompose: g is kept by degree, g[[k + 1]] the coefficient of x^k
.gd_right_decompose <- function(f, s, p) {
  n <- length(f) - 1L
  lc <- .gd_LC(f, 0L)
  r <- n %/% s
  g <- rep(list(gmp::as.bigz(0L)), s + 1L)
  g[[s + 1L]] <- gmp::as.bigz(1L)
  for (i in seq_len(s - 1L)) {
    cf <- gmp::as.bigz(0L)
    for (j in 0:(i - 1L)) {
      fc <- .gd_coef(f, n + j - i)
      if (fc == 0) next
      cf <- (cf + gmp::as.bigz(i - r * j) %% p * fc %% p * g[[s - j + 1L]]) %% p
    }
    g[[s - i + 1L]] <- (cf * .gd_inv(gmp::as.bigz(i * r) %% p * lc %% p, p)) %% p
  }
  .gd_strip(rev(g), 0L)
}
.gd_left_decompose <- function(f, h, p) {
  g <- list(); i <- 0L
  while (length(f)) {
    d <- .gd_div(f, h, 0L, p)
    if (.gd_deg(d$r) > 0) return(NULL)
    lc <- .gd_LC(d$r, 0L)
    if (lc != 0) g <- .gd_add_term(g, lc, i, 0L, p)
    f <- d$q; i <- i + 1L
  }
  g
}

# ---------------------------------------------------------------------------
# conversion between expressions and level u polynomials

.gd_from_expr <- function(e, gens, p) {
  q <- expr_to_qpoly(expand(as_sym(e)), gens, strict = TRUE)
  u <- length(gens) - 1L
  if (qp_is_zero(q)) return(list())
  cs <- .gf_from_q(lapply(seq_along(q$c), function(i) q$c[i]), gmp::as.bigz(p))
  cs <- gmp::as.bigz(cs)
  E <- matrix(q$e, ncol = length(gens))
  out <- list()
  for (i in seq_len(nrow(E))) if (cs[i] != 0) out <- .gd_add(out, .gd_monom(E[i, ], cs[i], u), u, p)
  out
}
.gd_monom <- function(ex, c, u) {
  if (u < 0L) return(c)
  inner <- .gd_monom(ex[-1L], c, u - 1L)
  c(list(inner), rep(list(.gd_z(u - 1L)), ex[[1L]]))
}
.gd_scalar_expr <- function(a, p, symmetric) {
  if (symmetric && a > p %/% 2) a <- a - p
  .int(a)
}
.gd_to_expr <- function(f, gens, u, p, symmetric = TRUE) {
  if (u < 0L) return(.gd_scalar_expr(f, p, symmetric))
  n <- length(f) - 1L
  if (n < 0L) return(S_Zero)
  x <- gens[[1L]]
  .add(lapply(seq_along(f), function(i) {
    if (.gd_isz(f[[i]], u - 1L)) return(S_Zero)
    .mul(list(.gd_to_expr(f[[i]], gens[-1L], u - 1L, p, symmetric), .pow(x, .int(n - i + 1L))))
  }))
}

# ---------------------------------------------------------------------------
# polytools functions with modulus = p (or Poly arguments over GF(p))

# the modulus of the options or of Poly arguments: list(p, symmetric), or NULL
.gf_modulus_of <- function(polys, opts) {
  m <- opts$modulus; sym <- !isFALSE(opts$symmetric)
  if (is.null(m)) for (P in polys) if (inherits(P, "Poly") && !is.null(P$modulus)) {
    m <- P$modulus
    if (isFALSE(P$symmetric)) sym <- FALSE
    break
  }
  if (is.null(m)) return(NULL)
  ms <- as_sym(m)
  if (ms$type != "Integer" || ms$value <= 0)
    stop("OptionError: 'modulus' must a positive integer, got ", .str(ms), call. = FALSE)
  list(p = gmp::as.bigz(as.character(ms$value)), symmetric = sym)
}

# CoercionFailed of SymPy when the generators leave other symbols in the
# coefficients: the first such coefficient, highest degree first
.gf_check_coeffs <- function(e, gs) {
  keys <- vapply(gs, `[[`, "", "key")
  if (all(vapply(free_symbols(e), function(s) s$key %in% keys, NA))) return(invisible())
  walk <- function(z, k) {
    if (k > length(gs)) return(if (is_number_type(z)) NULL else z)
    for (cf in .coeff_list(z, gs[[k]])) if (!is_zero_num(cf)) { r <- walk(cf, k + 1L); if (!is.null(r)) return(r) }
    NULL
  }
  bad <- walk(e, 1L)
  stop("CoercionFailed: expected an integer, got ", .str(if (is.null(bad)) e else bad), call. = FALSE)
}

# op on the polynomials (expressions or Polys) over GF(p); NULL when no
# modulus applies, so the caller goes on over its usual domain
.gf_hook <- function(op, polys, dots) {
  sp <- .split_poly_opts(dots)
  M <- .gf_modulus_of(polys, sp$opts)
  if (is.null(M)) return(NULL)
  p <- M$p; sym <- M$symmetric
  pl <- Filter(function(z) inherits(z, "Poly"), polys)
  es <- lapply(polys, function(z) expand(as_sym(if (inherits(z, "Poly")) as_expr(z) else z)))
  gens <- .as_gens(sp$gens)
  if (is.null(gens) && length(pl)) gens <- pl[[1L]]$gens
  gs <- .gens_of(es, gens)
  if (!length(gs)) return(NULL)
  for (e in es) .gf_check_coeffs(e, gs)
  u <- length(gs) - 1L
  F <- lapply(es, .gd_from_expr, gens = gs, p = p)
  ex <- function(f, v = u) {
    e <- if (v < 0L) .gd_scalar_expr(f, p, sym) else .gd_to_expr(f, gs[(u - v + 1L):(u + 1L)], v, p, sym)
    if (length(pl) && v == u) do.call(Poly, c(list(e), gs, list(modulus = .int(p), symmetric = sym))) else e
  }
  univariate <- function() if (u > 0L) stop("ValueError: univariate polynomial expected", call. = FALSE)
  not_divide <- function() stop("ExactQuotientFailed: ", .str(.gd_to_expr(F[[2L]], gs, u, p, sym)), " does not divide ",
                                .str(.gd_to_expr(F[[1L]], gs, u, p, sym)), call. = FALSE)
  f <- F[[1L]]; g <- if (length(F) > 1L) F[[2L]] else NULL
  switch(op,
    gcd = ex(.gd_gcd(f, g, u, p)),
    lcm = ex(.gd_lcm(f, g, u, p)),
    cofactors = { r <- .gd_cofactors(f, g, u, p); list(h = ex(r$h), cff = ex(r$cff), cfg = ex(r$cfg)) },
    div = { d <- .gd_div(f, g, u, p); list(q = ex(d$q), r = ex(d$r)) },
    quo = ex(.gd_quo(f, g, u, p)),
    rem = ex(.gd_rem(f, g, u, p)),
    exquo = { d <- .gd_div(f, g, u, p); if (length(d$r)) not_divide(); ex(d$q) },
    pdiv = { d <- .gd_pdiv(f, g, u, p); list(q = ex(d$q), r = ex(d$r)) },
    prem = ex(.gd_prem(f, g, u, p)),
    pquo = ex(.gd_pdiv(f, g, u, p)$q),
    pexquo = { d <- .gd_pdiv(f, g, u, p); if (length(d$r)) not_divide(); ex(d$q) },
    gcdex = { univariate(); r <- .gd_gcdex(f, g, p); list(s = ex(r$s), t = ex(r$t), h = ex(r$h)) },
    half_gcdex = { univariate(); r <- .gd_half_gcdex(f, g, p); list(s = ex(r$s), h = ex(r$h)) },
    invert = { univariate(); ex(.gd_invert(f, g, p)) },
    resultant = { r <- .gd_resultant(f, g, u, p); ex(r, u - 1L) },
    discriminant = { r <- .gd_discriminant(f, u, p); ex(r, u - 1L) },
    subresultants = lapply(.gd_subresultants(f, g, u, p), ex),
    monic = ex(.gd_ground_monic(f, u, p)),
    content = if (length(f)) S_One else S_Zero,
    primitive = list(content = if (length(f)) S_One else S_Zero, primitive = ex(f)),
    compose = ex(.gd_compose(f, g, u, p)),
    decompose = { univariate(); lapply(.gd_decompose(f, p), ex) },
    cancel = { r <- .gd_cofactors(f, g, u, p); list(ex(r$cff), ex(r$cfg)) },
    terms_gcd = {
      if (!length(f)) return(S_Zero)
      J <- .gd_terms_min(f, u)
      if (all(J == 0)) return(as_sym(polys[[1L]]))
      .mul(list(.monomial_vec_expr(J, gs), ex(.gd_terms_shift(f, J, u))))
    },
    sqf_part = {
      if (u == 0L) return(NULL)
      ex(Reduce(function(a, b) .gd_mul(a, b, u, p), lapply(.gd_factor_list(f, u, p)$factors, `[[`, 1L), .gd_one(u - 1L)))
    },
    stop("unknown operation ", op, call. = FALSE))
}

# dmp_terms_gcd: the smallest exponent of each variable over the terms, and
# f divided by that monomial
.gd_terms_min <- function(f, u) {
  if (u < 0L) return(integer())
  n <- length(f) - 1L
  nz <- which(!vapply(f, .gd_isz, NA, v = u - 1L))
  j0 <- n - max(nz) + 1L
  rest <- Reduce(function(a, b) pmin(a, b), lapply(f[nz], .gd_terms_min, u = u - 1L))
  c(j0, if (is.null(rest)) integer() else rest)
}
.gd_terms_shift <- function(f, J, u) {
  if (u < 0L) return(f)
  f <- f[seq_len(length(f) - J[[1L]])]
  lapply(f, function(cc) if (.gd_isz(cc, u - 1L)) cc else .gd_terms_shift(cc, J[-1L], u - 1L))
}

# ---------------------------------------------------------------------------
# Multivariate factorization over GF(p), which SymPy does not implement
# (NotImplementedError). Kronecker's method: with D above every partial
# degree, x_i -> x^(D^(u - i)) maps f to a univariate polynomial whose
# factorization over GF(p) is complete; every factor g of f maps to the
# product of a sub-multiset of those factors (up to a unit), and the map is
# undone exactly on polynomials of partial degrees below D. The smallest
# sub-multiset whose preimage divides f gives an irreducible factor.
# ponytail: the recombination is exponential in the number of univariate
# factors; replace by Hensel lifting (Bernardin, Wang) if that ever matters.

.gd_terms <- function(f, u) {
  if (u < 0L) return(list(list(e = integer(), c = f)))
  n <- length(f) - 1L
  out <- list()
  for (i in seq_along(f)) if (!.gd_isz(f[[i]], u - 1L))
    for (t in .gd_terms(f[[i]], u - 1L)) out[[length(out) + 1L]] <- list(e = c(n - i + 1L, t$e), c = t$c)
  out
}
.gd_kron <- function(f, u, D, p) {
  ts <- .gd_terms(f, u)
  E <- vapply(ts, function(t) sum(t$e * D^(u:0)), 0)
  v <- gmp::as.bigz(rep(0L, max(E) + 1L))
  for (k in seq_along(ts)) v[max(E) - E[k] + 1L] <- ts[[k]]$c
  if (.gf_big(p)) v else as.numeric(v)
}
.gd_unkron <- function(v, u, D, p) {
  n <- length(v) - 1L
  out <- list()
  for (j in seq_along(v)) if (v[j] != 0) {
    E <- n - j + 1
    e <- integer(u + 1L)
    for (i in (u + 1L):1L) { e[i] <- as.integer(E %% D); E <- E %/% D }
    out <- .gd_add(out, .gd_monom(e, gmp::as.bigz(v[j]), u), u, p)
  }
  out
}

# list(lc, factors = list(list(g, k))) with monic irreducible g
.gd_factor_list <- function(f, u, p) {
  P <- gmp::as.bigz(p); pu <- if (.gf_big(p)) P else as.numeric(P)
  if (!length(f)) return(list(lc = P * 0, factors = list()))
  lc <- .gd_ground_LC(f, u)
  f <- .gd_ground_monic(f, u, P)
  out <- list()
  repeat {
    ts <- .gd_terms(f, u)
    D <- max(unlist(lapply(ts, `[[`, "e"))) + 1
    if (D == 1) break
    uf <- .gf_factor(.gd_kron(f, u, D, pu), pu)$factors
    items <- unlist(lapply(uf, function(gk) rep(list(gk[[1L]]), gk[[2L]])), recursive = FALSE)
    m <- length(items)
    h <- NULL
    for (size in seq_len(m %/% 2L)) {
      cmb <- utils::combn(m, size, simplify = FALSE)
      seen <- character()
      for (S in cmb) {
        key <- paste(vapply(items[S], function(a) paste(as.character(a), collapse = ","), ""), collapse = ";")
        if (key %in% seen) next
        seen <- c(seen, key)
        cand <- .gd_ground_monic(.gd_unkron(Reduce(function(a, b) .gf_mul(a, b, pu), items[S]), u, D, P), u, P)
        if (all(unlist(lapply(.gd_terms(cand, u), `[[`, "e")) == 0)) next
        if (!length(.gd_div(f, cand, u, P)$r)) { h <- cand; break }
      }
      if (!is.null(h)) break
    }
    if (is.null(h)) h <- f
    k <- 0L
    repeat {
      d <- .gd_div(f, h, u, P)
      if (length(d$r)) break
      f <- d$q; k <- k + 1L
    }
    out[[length(out) + 1L]] <- list(h, k)
  }
  list(lc = lc, factors = out)
}
