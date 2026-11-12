# Univariate polynomials over the Gaussian domains ZZ_I and QQ_I, after
# sympy/polys/domains/gaussiandomains.py and the dense gcd and cancel
# routines (dup_rr_prs_gcd, dmp_cancel, PolyElement.cancel). A polynomial is
# list(re, im) of bigq vectors, leading coefficient first.

.gp <- function(re, im) list(re = gmp::as.bigq(re), im = gmp::as.bigq(im))
.gp_zero <- function() .gp(numeric(0), numeric(0))
.gp_len <- function(p) length(p$re)
.gp_is_zero <- function(p) !length(p$re)
.gp_strip <- function(p) {
  n <- length(p$re); k <- 1L
  while (k <= n && p$re[k] == 0 && p$im[k] == 0) k <- k + 1L
  if (k > n) .gp_zero() else if (k == 1L) p else .gp(p$re[k:n], p$im[k:n])
}

# complex scalars as c(re, im) bigq pairs
.gc_mul <- function(a, b) c(a[1] * b[1] - a[2] * b[2], a[1] * b[2] + a[2] * b[1])
.gc_div <- function(a, b) {
  d <- b[1]^2 + b[2]^2
  c((a[1] * b[1] + a[2] * b[2]) / d, (a[2] * b[1] - a[1] * b[2]) / d)
}
.gp_lc <- function(p) c(p$re[1], p$im[1])
.gp_scale <- function(p, s) .gp(p$re * s[1] - p$im * s[2], p$re * s[2] + p$im * s[1])
.gp_sub <- function(f, g) {
  nf <- .gp_len(f); ng <- .gp_len(g); n <- max(nf, ng)
  pad <- function(v, m) c(gmp::as.bigq(rep(0, n - m)), v)
  .gp_strip(.gp(pad(f$re, nf) - pad(g$re, ng), pad(f$im, nf) - pad(g$im, ng)))
}
.gp_shift <- function(p, k) .gp(c(p$re, gmp::as.bigq(rep(0, k))), c(p$im, gmp::as.bigq(rep(0, k))))

# division over QQ_I
.gp_divmod <- function(f, g) {
  dg <- .gp_len(g) - 1L
  q <- .gp(gmp::as.bigq(rep(0, max(.gp_len(f) - dg, 0L))), gmp::as.bigq(rep(0, max(.gp_len(f) - dg, 0L))))
  r <- f; lg <- .gp_lc(g)
  while (!.gp_is_zero(r) && .gp_len(r) - 1L >= dg) {
    j <- .gp_len(r) - 1L - dg
    c <- .gc_div(.gp_lc(r), lg)
    i <- .gp_len(q) - j
    q$re[i] <- c[1]; q$im[i] <- c[2]
    r <- .gp_sub(r, .gp_shift(.gp_scale(g, c), j))
  }
  list(q = .gp_strip(q), r = r)
}
.gp_monic_gcd <- function(f, g) {
  a <- f; b <- g
  while (!.gp_is_zero(b)) { r <- .gp_divmod(a, b)$r; a <- b; b <- r }
  .gp_scale(a, .gc_div(c(gmp::as.bigq(1), gmp::as.bigq(0)), .gp_lc(a)))
}

# Gaussian integers: quadrant normalization and the Euclidean gcd
.gi_quadrant <- function(a) {
  if (a[2] > 0) return(if (a[1] > 0) 0L else 1L)
  if (a[2] < 0) return(if (a[1] < 0) 2L else 3L)
  if (a[1] >= 0) 0L else 2L
}
# units[-quadrant] of (1, I, -1, -I)
.gi_canonical_unit <- function(a) switch(.gi_quadrant(a) + 1L,
  c(1, 0), c(0, -1), c(-1, 0), c(0, 1))
.gi_normalize <- function(a) .gc_mul(a, gmp::as.bigq(.gi_canonical_unit(a)))
.gi_mod <- function(a, b) {
  a <- gmp::as.bigz(a); b <- gmp::as.bigz(b)
  num <- c(a[1] * b[1] + a[2] * b[2], a[2] * b[1] - a[1] * b[2])
  cc <- b[1]^2 + b[2]^2
  q <- c((2 * num[1] + cc) %/% (2 * cc), (2 * num[2] + cc) %/% (2 * cc))
  qb <- c(q[1] * b[1] - q[2] * b[2], q[1] * b[2] + q[2] * b[1])
  gmp::as.bigq(c(a[1] - qb[1], a[2] - qb[2]))
}
.gi_gcd <- function(a, b) {
  while (!(b[1] == 0 && b[2] == 0)) { r <- .gi_mod(a, b); a <- b; b <- r }
  .gi_normalize(a)
}
.gp_content <- function(p) {
  cont <- gmp::as.bigq(c(0, 0))
  for (i in seq_len(.gp_len(p))) {
    cont <- .gi_gcd(cont, c(p$re[i], p$im[i]))
    if (cont[1] == 1 && cont[2] == 0) break
  }
  cont
}
.gp_den <- function(p) Reduce(gmp::lcm.bigz, as.list(c(gmp::denominator(p$re), gmp::denominator(p$im))), gmp::as.bigz(1))
.gp_is_field <- function(...) any(vapply(list(...), function(p) .gp_den(p) != 1, TRUE))

# cofactors over ZZ_I (dup_rr_prs_gcd): h normalized by the canonical unit of
# its leading coefficient times the gcd of the contents
.gp_cofactors <- function(f, g) {
  fc <- .gp_content(f); gc <- .gp_content(g)
  cg <- .gi_gcd(fc, gc)
  hm <- .gp_monic_gcd(.gp_scale(f, .gc_div(gmp::as.bigq(c(1, 0)), fc)), .gp_scale(g, .gc_div(gmp::as.bigq(c(1, 0)), gc)))
  hp <- .gp_scale(hm, gmp::as.bigq(c(.gp_den(hm), 0)))
  hp <- .gp_scale(hp, .gc_div(gmp::as.bigq(c(1, 0)), .gp_content(hp)))
  h <- .gp_scale(hp, .gc_mul(cg, gmp::as.bigq(.gi_canonical_unit(.gp_lc(hp)))))
  list(h = h, p = .gp_divmod(f, h)$q, q = .gp_divmod(g, h)$q)
}

# a Gaussian rational a + b*I as c(a, b) strings, or NULL
.gauss_num <- function(c) {
  rat <- function(t) t$type %in% c("Integer", "Rational")
  if (rat(c)) return(c(as.character(c$value), "0"))
  ca <- .coeff_add_split(c); a <- ca[[1L]]; t <- ca[[2L]]
  if (!rat(a)) return(NULL)
  if (is_const(t, "I")) return(c(as.character(a$value), "1"))
  sp <- as_coeff_mul_split(t)
  if (is_const(sp[[2L]], "I") && rat(sp[[1L]])) return(c(as.character(a$value), as.character(sp[[1L]]$value)))
  NULL
}
.gp_from_coeffs <- function(cl) {
  re <- character(length(cl)); im <- character(length(cl))
  for (i in seq_along(cl)) {
    pc <- .gauss_num(cl[[i]])
    if (is.null(pc)) return(NULL)
    re[i] <- pc[1]; im[i] <- pc[2]
  }
  .gp_strip(.gp(gmp::as.bigq(re), gmp::as.bigq(im)))
}
.gc_expr <- function(a) .add(list(make_rational(a[1]), .mul(list(make_rational(a[2]), S_I))))
.gp_as_expr <- function(p, x) {
  d <- .gp_len(p) - 1L
  .add(lapply(seq_len(.gp_len(p)), function(k)
    .mul(list(.gc_expr(c(p$re[k], p$im[k])), .pow(x, .int(d - k + 1L))))))
}
# coefficient lists in x with Gaussian rational coefficients and some I
.gauss_polys <- function(n, d, x) {
  cl <- tryCatch(list(.coeff_list(n, x), .coeff_list(d, x)), error = function(err) NULL)
  if (is.null(cl)) return(NULL)
  if (!any(vapply(c(cl[[1L]], cl[[2L]]), function(c) .has1(c, S_I), TRUE))) return(NULL)
  f <- .gp_from_coeffs(cl[[1L]]); g <- .gp_from_coeffs(cl[[2L]])
  if (is.null(f) || is.null(g)) return(NULL)
  list(f, g)
}

# PolyElement.cancel over ZZ_I or QQ_I; returns the expression p/q
.gauss_cancel_expr <- function(f, g, x) {
  if (.gp_is_zero(f)) return(S_Zero)
  field <- .gp_is_field(f, g)
  if (field) {
    cq <- .gp_den(f); cp <- .gp_den(g)
    f <- .gp_scale(f, gmp::as.bigq(c(cq, 0))); g <- .gp_scale(g, gmp::as.bigq(c(cp, 0)))
  }
  cf <- .gp_cofactors(f, g); p <- cf$p; q <- cf$q
  if (field) {
    k <- gmp::gcd.bigz(cp, cq)
    p <- .gp_scale(p, gmp::as.bigq(c(cp / k, 0))); q <- .gp_scale(q, gmp::as.bigq(c(cq / k, 0)))
  }
  u <- gmp::as.bigq(.gi_canonical_unit(.gp_lc(q)))
  p <- .gp_scale(p, u); q <- .gp_scale(q, u)
  .div(.gp_as_expr(p, x), .gp_as_expr(q, x))
}

# Poly.cancel(include = False) over ZZ_I or QQ_I: list(c, p, q)
.gauss_poly_cancel <- function(f, g, x) {
  if (.gp_is_zero(f)) return(list(S_One, S_Zero, S_One))
  field <- .gp_is_field(f, g)
  c <- S_One
  if (field) {
    cq <- .gp_den(f); cp <- .gp_den(g)
    f <- .gp_scale(f, gmp::as.bigq(c(cq, 0))); g <- .gp_scale(g, gmp::as.bigq(c(cp, 0)))
    k <- gmp::gcd.bigz(cp, cq)
    c <- make_rational(gmp::as.bigq(cp / k, cq / k))
  }
  cf <- .gp_cofactors(f, g)
  list(c, .gp_as_expr(cf$p, x), .gp_as_expr(cf$q, x))
}
