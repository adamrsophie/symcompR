# Dense univariate polynomials over a SymPy domain, after
# sympy/polys/densebasic.py, densearith.py, densetools.py, euclidtools.py and
# sqfreetools.py. A dpoly is list(c = coefficients as expressions, leading
# first, x = the generator, K = the domain from construct_domain()). Each
# coefficient is kept in the canonical form of its domain, so the
# normalizations of SymPy's algorithms (monic gcds over fields, primitive
# gcds with a positive leading coefficient over rings) carry over.

.dk <- function(K) .dom(K)$kind
.dk_field <- function(K) .dk(K) %in% c("QQ", "RR", "CC", "QQ_I", "FracField", "EX", "GF", "AlgField")
.dk_gens <- function(K) .dom(K)$gens
.dk_ground <- function(K) .dom(K)$dom

# the domain of the coefficients of e in x (Poly(e, x) with field = TRUE when asked)
.dp_domain <- function(coeffs, field = FALSE) {
  K <- construct_domain(coeffs)[[1L]]
  if (field && !.dk_field(K)) K <- .dom_get_field(K)
  K
}
.dom_get_field <- function(K) {
  switch(.dk(K), ZZ = QQ, ZZ_I = QQ_I,
         PolyRing = do.call(frac_field, c(list(if (.dk(.dk_ground(K)) == "QQ") ZZ else .dk_ground(K)), .dk_gens(K))),
         K)
}
.dom_get_ring <- function(K) {
  switch(.dk(K), QQ = ZZ, QQ_I = ZZ_I,
         FracField = do.call(poly_ring, c(list(.dk_ground(K)), .dk_gens(K))),
         K)
}

# canonical form of an element
.dn <- function(K, a) {
  switch(.dk(K),
    ZZ = , QQ = , ZZ_I = , QQ_I = a,
    PolyRing = expand(a),
    FracField = cancel(a),
    .ex_simp(a))
}

.d_zero <- function(a) is_zero_num(a)
.d_one <- function(a) is_one_num(a)

# lex leading coefficient of a polynomial expression in gens
.lex_lc <- function(a, gens) {
  if (!length(gens)) return(a)
  if (.contains_key(a, S_I$key)) {
    iI <- .I_dummy()
    return(expand(subs(.lex_lc(expand(subs(a, S_I, iI)), c(gens, list(iI))), iI, S_I)))
  }
  p <- expr_to_qpoly(a, gens, strict = FALSE)
  if (is.null(p) || qp_is_zero(p)) return(a)
  make_rational(qp_lc(p))
}

.dom_is_negative <- function(K, a) {
  switch(.dk(K),
    ZZ = , QQ = is_number_type(a) && num_sign(a) < 0,
    # GaussianDomain.is_negative is False for every element
    ZZ_I = , QQ_I = FALSE,
    PolyRing = { lc <- .lex_lc(a, .dk_gens(K)); is_number_type(lc) && num_sign(lc) < 0 },
    FracField = { nd <- .as_numer_denom(a); lc <- .lex_lc(expand(nd[[1L]]), .dk_gens(K)); is_number_type(lc) && num_sign(lc) < 0 },
    { cm <- as_coeff_mul_split(a)[[1L]]; is_number_type(cm) && num_sign(cm) < 0 })
}

# Gaussian integers as (x, y) with x + y*I
.gi_parts <- function(e) {
  ri <- .as_real_imag(as_sym(e))
  c(gmp::as.bigq(ri[[1L]]$value), gmp::as.bigq(ri[[2L]]$value))
}
.gi_make <- function(x, y) .add(list(make_rational(gmp::as.bigq(x)), .mul(list(make_rational(gmp::as.bigq(y)), S_I))))
.gq_denom <- function(e) { v <- .gi_parts(e); gmp::lcm.bigz(gmp::denominator(v[1]), gmp::denominator(v[2])) }
# GaussianDomain.gcd and canonical_unit on expressions, by the pair
# arithmetic of R/polys-gaussian.R
.zzi_gcd <- function(a, b) {
  g <- .gi_gcd(.gi_parts(a), .gi_parts(b))
  .gi_make(g[1], g[2])
}
.zzi_canonical_unit <- function(a) { u <- .gi_canonical_unit(.gi_parts(a)); .gi_make(u[1], u[2]) }
# K.canonical_unit
.dom_canonical_unit <- function(K, a) {
  if (.dk(K) %in% c("ZZ_I", "QQ_I")) return(.zzi_canonical_unit(a))
  if (.dom_is_negative(K, a)) S_NegativeOne else S_One
}
.dom_is_nonnegative <- function(K, a) if (.dk(K) %in% c("ZZ_I", "QQ_I")) FALSE else !.dom_is_negative(K, a)

# K.gcd for elements
.dom_gcd <- function(K, a, b) {
  k <- .dk(K)
  if (k %in% c("ZZ", "QQ")) {
    if (.d_zero(a)) return(Abs(b)); if (.d_zero(b)) return(Abs(a))
    qa <- gmp::as.bigq(a$value); qb <- gmp::as.bigq(b$value)
    n <- gmp::gcd.bigz(gmp::numerator(qa), gmp::numerator(qb)); d <- gmp::lcm.bigz(gmp::denominator(qa), gmp::denominator(qb))
    return(make_rational(gmp::as.bigq(n, d)))
  }
  if (k == "ZZ_I") return(.zzi_gcd(a, b))
  if (k == "QQ_I") {
    # Field.gcd: the ZZ_I gcd of the numerators over the lcm of the denominators
    if (.d_zero(a) && .d_zero(b)) return(S_Zero)
    da <- .gq_denom(a); db <- .gq_denom(b)
    p <- .zzi_gcd(expand(.mul(list(a, .zs(da)))), expand(.mul(list(b, .zs(db)))))
    return(.div(p, .zs(gmp::lcm.bigz(da, db))))
  }
  if (k == "PolyRing") return(.mpoly_gcd(a, b, .dk_gens(K), .dk(.dk_ground(K)) == "QQ"))
  if (k == "FracField") {
    na <- .as_numer_denom(a); nb <- .as_numer_denom(b)
    gens <- .dk_gens(K)
    p <- .mpoly_gcd(expand(na[[1L]]), expand(nb[[1L]]), gens, FALSE)
    q <- .mpoly_lcm(expand(na[[2L]]), expand(nb[[2L]]), gens)
    return(cancel(.div(p, q)))
  }
  S_One
}
.dom_lcm <- function(K, a, b) .dn(K, .div(.mul(list(a, b)), .dom_gcd(K, a, b)))

# multivariate gcd over ZZ[gens] (content included, positive lex leading
# coefficient) or QQ[gens] (monic)
.I_dummy <- function() { if (is.null(.mi_state$Idummy)) .mi_state$Idummy <- Dummy("ImaginaryUnit"); .mi_state$Idummy }
.mpoly_gcd <- function(a, b, gens, qq) {
  if (.d_zero(a) && .d_zero(b)) return(S_Zero)
  if (.contains_key(a, S_I$key) || .contains_key(b, S_I$key)) {
    # ponytail: I as a generator, so Gaussian common factors such as 1 + I
    # are not found; the ZZ_I and QQ_I gcds of SymPy normalize by quadrant
    iI <- .I_dummy()
    h <- .mpoly_gcd(expand(subs(a, S_I, iI)), expand(subs(b, S_I, iI)), c(gens, list(iI)), qq)
    return(expand(subs(h, iI, S_I)))
  }
  pa <- expr_to_qpoly(a, gens, strict = FALSE); pb <- expr_to_qpoly(b, gens, strict = FALSE)
  if (qp_is_zero(pa)) h <- pb else if (qp_is_zero(pb)) h <- pa else h <- qp_gcd(pa, pb)
  if (qq) return(qpoly_to_expr(qp_monic(h), gens))
  hp <- qp_primitive(h)$poly
  if (qp_lc(hp) < 0) hp <- qp_neg(hp)
  ca <- if (qp_is_zero(pa)) gmp::as.bigz(0) else .bq_int_content(pa)
  cb <- if (qp_is_zero(pb)) gmp::as.bigz(0) else .bq_int_content(pb)
  qpoly_to_expr(qp_scale(hp, gmp::as.bigq(gmp::gcd.bigz(ca, cb))), gens)
}
.mpoly_lcm <- function(a, b, gens) {
  g <- .mpoly_gcd(a, b, gens, FALSE)
  r <- cancel(.div(.mul(list(a, b)), g))
  lc <- .lex_lc(r, gens)
  if (is_number_type(lc) && num_sign(lc) < 0) r <- expand(.neg(r))
  expand(r)
}

# K.quo (exact quotient in rings, division in fields)
.dom_quo <- function(K, a, b) {
  k <- .dk(K)
  if (k == "ZZ") return(.zs(gmp::numerator(gmp::as.bigq(a$value)) %/% gmp::numerator(gmp::as.bigq(b$value))))
  if (k == "PolyRing") return(expand(cancel(.div(a, b))))
  .dn(K, .div(a, b))
}

# ---------------------------------------------------------------------------
# dense polynomials

dpoly <- function(c, x, K) structure(list(c = .dp_strip(c), x = x, K = K), class = "dpoly")
.dp_strip <- function(c) { while (length(c) && .d_zero(c[[1L]])) c <- c[-1L]; c }

# Poly(e, x, domain = K); K from the coefficients when NULL
.dp_from <- function(e, x, K = NULL, field = FALSE) {
  cl <- .coeff_list(expand(as_sym(e)), x)
  if (is.null(K)) K <- .dp_domain(if (length(cl)) cl else list(S_Zero), field)
  dpoly(lapply(cl, function(z) .dn(K, z)), x, K)
}
.dp_expr <- function(p) {
  n <- length(p$c)
  if (!n) return(S_Zero)
  .add(lapply(seq_len(n), function(i) .mul(list(p$c[[i]], .pow(p$x, .int(n - i))))))
}
.dp_deg <- function(p) if (length(p$c)) length(p$c) - 1L else -Inf
.dp_LC <- function(p) if (length(p$c)) p$c[[1L]] else S_Zero
.dp_TC <- function(p) if (length(p$c)) p$c[[length(p$c)]] else S_Zero
.dp_zero <- function(p) !length(p$c)
.dp_set <- function(p, c) dpoly(c, p$x, p$K)
.dp_const <- function(p, a) .dp_set(p, list(a))
.dp_one <- function(p) .dp_const(p, S_One)
.dp_convert <- function(p, K) dpoly(lapply(p$c, function(z) .dn(K, z)), p$x, K)

.dp_add <- function(f, g) {
  a <- f$c; b <- g$c; K <- f$K
  la <- length(a); lb <- length(b)
  if (la < lb) { a <- c(rep(list(S_Zero), lb - la), a) } else if (lb < la) { b <- c(rep(list(S_Zero), la - lb), b) }
  .dp_set(f, Map(function(u, v) .dn(K, .add(list(u, v))), a, b))
}
.dp_neg <- function(f) .dp_set(f, lapply(f$c, function(z) .dn(f$K, .neg(z))))
.dp_sub <- function(f, g) .dp_add(f, .dp_neg(g))
.dp_mul <- function(f, g) {
  if (.dp_zero(f) || .dp_zero(g)) return(.dp_set(f, list()))
  a <- f$c; b <- g$c; K <- f$K
  n <- length(a) + length(b) - 1L
  out <- vector("list", n)
  for (k in seq_len(n)) {
    terms <- list()
    for (i in max(1L, k - length(b) + 1L):min(k, length(a))) terms[[length(terms) + 1L]] <- .mul(list(a[[i]], b[[k - i + 1L]]))
    out[[k]] <- .dn(K, .add(terms))
  }
  .dp_set(f, out)
}
.dp_mul_ground <- function(f, c) .dp_set(f, lapply(f$c, function(z) .dn(f$K, .mul(list(z, c)))))
.dp_quo_ground <- function(f, c) .dp_set(f, lapply(f$c, function(z) .dom_quo(f$K, z, c)))
.dp_mul_term <- function(f, c, j) if (.dp_zero(f) || .d_zero(c)) .dp_set(f, list()) else
  .dp_set(f, c(lapply(f$c, function(z) .dn(f$K, .mul(list(z, c)))), rep(list(S_Zero), j)))
.dp_add_term <- function(f, c, j) .dp_add(f, .dp_set(f, c(list(c), rep(list(S_Zero), j))))
.dp_pow <- function(f, n) { r <- .dp_one(f); for (i in seq_len(n)) r <- .dp_mul(r, f); r }
.dp_diff <- function(f) {
  n <- length(f$c) - 1L
  if (n < 1L) return(.dp_set(f, list()))
  .dp_set(f, lapply(seq_len(n), function(i) .dn(f$K, .mul(list(.int(n - i + 1L), f$c[[i]])))))
}
.dp_eval <- function(f, a) .dn(f$K, .add(lapply(seq_along(f$c), function(i) .mul(list(f$c[[i]], .pow(a, .int(length(f$c) - i)))))))
.dp_monic <- function(f) {
  if (.dp_zero(f)) return(f)
  lc <- .dp_LC(f)
  if (.d_one(lc)) return(f)
  if (.dk_field(f$K)) .dp_set(f, lapply(f$c, function(z) .dn(f$K, .div(z, lc)))) else .dp_quo_ground(f, lc)
}

# dup_content and dup_primitive
.dp_content <- function(f) {
  if (.dp_zero(f)) return(S_Zero)
  K <- f$K; cont <- S_Zero
  for (c in f$c) {
    cont <- .dom_gcd(K, cont, c)
    if (.dk(K) != "QQ" && .d_one(cont)) break
  }
  cont
}
.dp_primitive <- function(f) {
  if (.dp_zero(f)) return(list(S_Zero, f))
  cont <- .dp_content(f)
  if (.d_one(cont)) return(list(cont, f))
  list(cont, .dp_quo_ground(f, cont))
}

# dup_ff_div / dup_rr_div
.dp_div <- function(f, g) {
  if (.dp_zero(g)) stop("polynomial division by zero")
  K <- f$K
  df <- .dp_deg(f); dg <- .dp_deg(g)
  q <- .dp_set(f, list()); r <- f
  if (df < dg) return(list(q, r))
  lcg <- .dp_LC(g); field <- .dk_field(K)
  dr <- df
  repeat {
    lcr <- .dp_LC(r)
    if (!field) {
      cc <- cancel(.div(lcr, lcg))
      if (!.dom_divides(K, lcg, lcr)) break
      c <- .dom_quo(K, lcr, lcg)
    } else c <- .dn(K, .div(lcr, lcg))
    j <- dr - dg
    q <- .dp_add_term(q, c, j)
    r <- .dp_sub(r, .dp_mul_term(g, c, j))
    odr <- dr; dr <- .dp_deg(r)
    # EX coefficients whose canonical forms do not cancel structurally: the
    # leading term is zero by construction
    while (dr == odr && .dk(K) == "EX" && .dp_numzero(.dp_LC(r))) { r <- .dp_set(r, r$c[-1L]); dr <- .dp_deg(r) }
    if (dr < dg) break
    if (dr == odr && .dk(K) %in% c("RR", "CC")) {
      # dup_ff_div: remove the leading term left by rounding error
      r <- .dp_set(r, r$c[-1L]); dr <- .dp_deg(r)
      if (dr < dg) break
    } else if (!(dr < odr)) stop("polynomial division failed")
  }
  list(q, r)
}
.dp_numzero <- function(a) {
  v <- tryCatch(.neval(a, list()), error = function(e) NA)
  !is.na(v) && is.finite(Mod(v)) && Mod(v) < 1e-9
}
.dom_divides <- function(K, a, b) {
  k <- .dk(K)
  if (k == "ZZ") return(gmp::mod.bigz(gmp::numerator(gmp::as.bigq(b$value)), gmp::numerator(gmp::as.bigq(a$value))) == 0)
  if (k == "PolyRing") { q <- cancel(.div(b, a)); return(isTRUE(is_polynomial(q)) && !.has_den_syms(q, .dk_gens(K))) }
  TRUE
}
.has_den_syms <- function(q, gens) { d <- .as_numer_denom(q)[[2L]]; any(vapply(gens, function(g) .has1(d, g), NA)) }
.dp_quo <- function(f, g) .dp_div(f, g)[[1L]]
.dp_rem <- function(f, g) .dp_div(f, g)[[2L]]

.dp_prem <- function(f, g) {
  df <- .dp_deg(f); dg <- .dp_deg(g)
  if (.dp_zero(g)) stop("polynomial division by zero")
  r <- f; dr <- df
  if (df < dg) return(r)
  N <- df - dg + 1L
  lcg <- .dp_LC(g)
  repeat {
    lcr <- .dp_LC(r)
    j <- dr - dg; N <- N - 1L
    R <- .dp_mul_ground(r, lcg)
    G <- .dp_mul_term(g, lcr, j)
    r <- .dp_sub(R, G)
    odr <- dr; dr <- .dp_deg(r)
    if (dr < dg) break
    if (!(dr < odr)) stop("polynomial division failed")
  }
  .dp_mul_ground(r, .pow(lcg, .int(N)))
}

# dup_inner_subresultants: list(R, S)
.dp_inner_subresultants <- function(f, g) {
  K <- f$K
  n <- .dp_deg(f); m <- .dp_deg(g)
  if (n < m) { t <- f; f <- g; g <- t; t <- n; n <- m; m <- t }
  if (.dp_zero(f)) return(list(list(), list()))
  if (.dp_zero(g)) return(list(list(f), list(S_One)))
  R <- list(f, g)
  d <- n - m
  b <- .pow(S_NegativeOne, .int(d + 1L))
  h <- .dp_mul_ground(.dp_prem(f, g), b)
  lc <- .dp_LC(g)
  c <- .dn(K, .pow(lc, .int(d)))
  S <- list(S_One, c)
  c <- .dn(K, .neg(c))
  while (!.dp_zero(h)) {
    k <- .dp_deg(h)
    R[[length(R) + 1L]] <- h
    f <- g; g <- h; d <- m - k; m <- k
    b <- .dn(K, .mul(list(.neg(lc), .pow(c, .int(d)))))
    h <- .dp_quo_ground(.dp_prem(f, g), b)
    lc <- .dp_LC(g)
    if (d > 1L) c <- .dom_quo(K, .dn(K, .pow(.neg(lc), .int(d))), .dn(K, .pow(c, .int(d - 1L))))
    else c <- .dn(K, .neg(lc))
    S[[length(S) + 1L]] <- .dn(K, .neg(c))
  }
  list(R, S)
}
.dp_subresultants <- function(f, g) .dp_inner_subresultants(f, g)[[1L]]
# dup_prs_resultant: list(res, R)
.dp_prs_resultant <- function(f, g) {
  if (.dp_zero(f) || .dp_zero(g)) return(list(S_Zero, list()))
  rs <- .dp_inner_subresultants(f, g)
  R <- rs[[1L]]; S <- rs[[2L]]
  if (.dp_deg(R[[length(R)]]) > 0) return(list(S_Zero, R))
  list(S[[length(S)]], R)
}

# dup_inner_gcd: list(h, cff, cfg)
.dp_inner_gcd <- function(f, g) {
  K <- f$K
  # dup_inner_gcd: over RR the gcd is computed exactly over QQ
  if (.dk(K) == "RR") {
    toq <- function(p) dpoly(lapply(p$c, function(z) if (z$type == "Float") make_rational(gmp::as.bigq(z$value)) else z), p$x, QQ)
    tor <- function(p) dpoly(lapply(p$c, function(z) Float(as.double(num_value_dbl(z)))), p$x, K)
    r <- .dp_inner_gcd(toq(f), toq(g))
    return(lapply(r, tor))
  }
  if (.dk_field(K)) {
    if (.dp_zero(f) && .dp_zero(g)) return(list(f, f, f))
    if (.dp_zero(f)) { h <- .dp_monic(g); return(list(h, .dp_set(f, list()), .dp_const(f, .dn(K, .inv(.dp_LC(g)))))) }
    if (.dp_zero(g)) { h <- .dp_monic(f); return(list(h, .dp_const(f, .dn(K, .inv(.dp_LC(f)))), .dp_set(f, list()))) }
    h <- .dp_monic(.dp_subresultants(f, g)[[length(.dp_subresultants(f, g))]])
    return(list(h, .dp_quo(f, h), .dp_quo(g, h)))
  }
  # dup_rr_prs_gcd
  if (.dp_zero(f) && .dp_zero(g)) return(list(f, f, f))
  if (.dp_zero(f)) {
    if (.dom_is_nonnegative(K, .dp_LC(g))) return(list(g, .dp_set(f, list()), .dp_one(f)))
    return(list(.dp_neg(g), .dp_set(f, list()), .dp_const(f, S_NegativeOne)))
  }
  if (.dp_zero(g)) {
    if (.dom_is_nonnegative(K, .dp_LC(f))) return(list(f, .dp_one(f), .dp_set(f, list())))
    return(list(.dp_neg(f), .dp_const(f, S_NegativeOne), .dp_set(f, list())))
  }
  pf <- .dp_primitive(f); pg <- .dp_primitive(g)
  c <- .dom_gcd(K, pf[[1L]], pg[[1L]])
  sr <- .dp_subresultants(pf[[2L]], pg[[2L]])
  h <- .dp_primitive(sr[[length(sr)]])[[2L]]
  c <- .dn(K, .mul(list(c, .dom_canonical_unit(K, .dp_LC(h)))))
  h <- .dp_mul_ground(h, c)
  list(h, .dp_quo(f, h), .dp_quo(g, h))
}
.dp_gcd <- function(f, g) .dp_inner_gcd(f, g)[[1L]]

# dup_sqf_list: list(coeff, list(list(factor, k), ...))
.dp_sqf_list <- function(f, all = FALSE) {
  K <- f$K
  if (.dk_field(K)) { coeff <- .dp_LC(f); f <- .dp_monic(f) }
  else {
    pr <- .dp_primitive(f); coeff <- pr[[1L]]; f <- pr[[2L]]
    if (.dom_is_negative(K, .dp_LC(f))) { f <- .dp_neg(f); coeff <- .dn(K, .neg(coeff)) }
  }
  if (.dp_deg(f) <= 0) return(list(coeff, list()))
  result <- list(); i <- 1L
  h <- .dp_diff(f)
  ig <- .dp_inner_gcd(f, h); p <- ig[[2L]]; q <- ig[[3L]]
  repeat {
    d <- .dp_diff(p)
    h <- .dp_sub(q, d)
    if (.dp_zero(h)) { result[[length(result) + 1L]] <- list(p, i); break }
    ig <- .dp_inner_gcd(p, h); g <- ig[[1L]]; p <- ig[[2L]]; q <- ig[[3L]]
    if (all || .dp_deg(g) > 0) result[[length(result) + 1L]] <- list(g, i)
    i <- i + 1L
  }
  list(coeff, result)
}

# dup_half_gcdex, dup_gcdex, dup_invert over a field
.dp_half_gcdex <- function(f, g) {
  a <- .dp_one(f); b <- .dp_set(f, list())
  while (!.dp_zero(g)) {
    qr <- .dp_div(f, g)
    f <- g; g <- qr[[2L]]
    nb <- .dp_sub(a, .dp_mul(qr[[1L]], b)); a <- b; b <- nb
  }
  a <- .dp_set(a, lapply(a$c, function(z) .dn(f$K, .div(z, .dp_LC(f)))))
  list(a, .dp_monic(f))
}
.dp_gcdex <- function(f, g) {
  sh <- .dp_half_gcdex(f, g)
  F <- .dp_sub(sh[[2L]], .dp_mul(sh[[1L]], f))
  list(sh[[1L]], .dp_quo(F, g), sh[[2L]])
}
.dp_invert <- function(f, g) {
  sh <- .dp_half_gcdex(f, g)
  if (.dp_deg(sh[[2L]]) == 0 && .d_one(.dp_LC(sh[[2L]]))) return(.dp_rem(sh[[1L]], g))
  stop("zero divisor")
}

# Poly.to_field / retract-like conversion for the auto=True Poly methods
.dp_to_field <- function(f) if (.dk_field(f$K)) f else .dp_convert(f, .dom_get_field(f$K))
.dp_to_ring <- function(f) {
  K <- f$K
  if (!.dk_field(K) || .dk(K) %in% c("EX", "RR", "CC", "GF")) return(f)
  R <- .dom_get_ring(K)
  if (identical(.dk(R), .dk(K))) return(f)
  ok <- all(vapply(f$c, function(z) {
    if (.dk(R) == "ZZ") return(z$type == "Integer")
    nd <- .as_numer_denom(z); !any(vapply(.dk_gens(K), function(g) .has1(nd[[2L]], g), NA)) && is_number_type(nd[[2L]]) && is_one_num(nd[[2L]])
  }, NA))
  if (ok) .dp_convert(f, R) else f
}

# ---------------------------------------------------------------------------
# Domain.unify and unify_composite (sympy/polys/domains/domain.py)

.dk_same <- function(K0, K1) identical(format(K0), format(K1))
.dk_comp <- function(K) .dk(K) %in% c("PolyRing", "FracField")

.unify_gens <- function(f_gens, g_gens) {
  fk <- vapply(f_gens, `[[`, "", "key"); gk <- vapply(g_gens, `[[`, "", "key")
  if (identical(fk, gk)) return(f_gens)
  common <- f_gens[fk %in% gk]
  ck <- vapply(common, `[[`, "", "key")
  k <- 0L
  for (i in seq_along(g_gens)) if (gk[i] %in% ck) { k <- k + 1L; g_gens[[i]] <- common[[k]] }
  gk <- vapply(g_gens, `[[`, "", "key")
  gens <- list()
  for (gen in common) {
    i <- match(gen$key, vapply(f_gens, `[[`, "", "key"))
    gens <- c(gens, f_gens[seq_len(i - 1L)]); f_gens <- f_gens[-seq_len(i)]
    j <- match(gen$key, vapply(g_gens, `[[`, "", "key"))
    gens <- c(gens, g_gens[seq_len(j - 1L)]); g_gens <- g_gens[-seq_len(j)]
    gens <- c(gens, list(gen))
  }
  c(gens, f_gens, g_gens)
}

.dom_unify <- function(K0, K1) {
  if (.dk_same(K0, K1)) return(K0)
  k0 <- .dk(K0); k1 <- .dk(K1)
  if (k0 == "EX") return(K0)
  if (k1 == "EX") return(K1)
  if (.dk_comp(K0) || .dk_comp(K1)) return(.dom_unify_composite(K0, K1))
  if (k1 == "CC") { t <- K0; K0 <- K1; K1 <- t; k0 <- .dk(K0); k1 <- .dk(K1) }
  if (k0 == "CC") return(K0)
  if (k1 == "RR") { t <- K0; K0 <- K1; K1 <- t; k0 <- .dk(K0); k1 <- .dk(K1) }
  if (k0 == "RR") return(if (k1 %in% c("ZZ_I", "QQ_I")) CC else K0)
  if (k1 == "AlgField") { t <- K0; K0 <- K1; K1 <- t; k0 <- .dk(K0); k1 <- .dk(K1) }
  if (k0 == "AlgField") return(K0)
  if (k0 == "QQ_I") return(K0)
  if (k1 == "QQ_I") return(K1)
  if (k0 == "ZZ_I") return(if (k1 == "QQ") QQ_I else K0)
  if (k1 == "ZZ_I") return(if (k0 == "QQ") QQ_I else K1)
  if (k0 == "QQ") return(K0)
  if (k1 == "QQ") return(K1)
  if (k0 == "ZZ") return(K0)
  if (k1 == "ZZ") return(K1)
  EX
}

.dom_unify_composite <- function(K0, K1) {
  g0 <- if (.dk_comp(K0)) .dk_ground(K0) else K0
  g1 <- if (.dk_comp(K1)) .dk_ground(K1) else K1
  s0 <- if (.dk_comp(K0)) .dk_gens(K0) else list()
  s1 <- if (.dk_comp(K1)) .dk_gens(K1) else list()
  domain <- .dom_unify(g0, g1)
  symbols <- .unify_gens(s0, s1)
  k0 <- .dk(K0); k1 <- .dk(K1)
  if (((k0 == "FracField" && k1 == "PolyRing") || (k1 == "FracField" && k0 == "PolyRing")) &&
      (!.dk_field(g0) || !.dk_field(g1)) && .dk_field(domain) && .dk(domain) %in% c("QQ", "QQ_I"))
    domain <- .dom_get_ring(domain)
  cls <- if (.dk_comp(K0) && (!.dk_comp(K1) || k0 == "FracField" || k1 == "PolyRing")) k0 else k1
  if (cls == "PolyRing") do.call(poly_ring, c(list(domain), symbols)) else do.call(frac_field, c(list(domain), symbols))
}
