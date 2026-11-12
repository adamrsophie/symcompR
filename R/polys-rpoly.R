# Univariate Poly objects with SymPy's domain semantics, used by the ports
# of the integration algorithms (risch.py, rde.py, prde.py). A P object is a
# dpoly (R/polys-dpoly.R): coefficients in a domain, a generator, and the
# operations of sympy.polys.polytools.Poly: binary operations unify the
# domains, division, gcdex and invert convert a ring to its field and
# retract the result (auto = True), and the domain of a new Poly comes from
# construct_domain of its coefficients.

.P <- function(e, x, domain = NULL, field = FALSE) {
  if (inherits(e, "dpoly")) {
    if (identical(e$x$key, as_sym(x)$key) && is.null(domain) && !field) return(e)
    e <- .dp_expr(e)
  }
  x <- as_sym(x)
  e <- as_sym(e)
  if (!.is_comm(e)) stop(structure(class = c("PolynomialError", "error", "condition"), list(message = "non-commutative expressions are not supported", call = NULL)))
  cl <- .coeff_list(expand(e), x)
  if (is.null(domain)) domain <- .dp_domain(if (length(cl)) cl else list(S_Zero), field)
  else if (field && !.dk_field(domain)) domain <- .dom_get_field(domain)
  dpoly(lapply(cl, function(z) .dn(domain, z)), x, domain)
}
# Poly(e, x) or NULL when e is not a polynomial in x (as_poly)
P_or_null <- function(e, x, domain = NULL, field = FALSE) {
  e <- as_sym(if (inherits(e, "dpoly")) .dp_expr(e) else e)
  if (!.is_comm(e) || !isTRUE(is_polynomial(e, as_sym(x)))) return(NULL)
  tryCatch(.P(e, x, domain, field), error = function(err) NULL)
}

P_expr <- function(p) .dp_expr(p)
P_deg <- function(p) .dp_deg(p)
P_LC <- function(p) .dp_LC(p)
P_TC <- function(p) .dp_TC(p)
P_is_zero <- function(p) .dp_zero(p)
P_is_one <- function(p) length(p$c) == 1L && is_one_num(p$c[[1L]])
P_is_ground <- function(p) length(p$c) <= 1L
P_is_monomial <- function(p) sum(!vapply(p$c, is_zero_num, NA)) <= 1L
P_has <- function(p, x) .has1(P_expr(p), as_sym(x))
# Poly.nth(k): the coefficient of x^k
P_nth <- function(p, k) { n <- P_deg(p); if (k > n || k < 0) S_Zero else p$c[[n - k + 1L]] }
# Poly.coeffs(): the nonzero coefficients, leading first
P_coeffs <- function(p) Filter(function(z) !is_zero_num(z), p$c)
P_monoms <- function(p) { n <- P_deg(p); (n:0)[!vapply(p$c, is_zero_num, NA)] }
P_all_coeffs <- function(p) if (length(p$c)) p$c else list(S_Zero)

# common domain of two polys (Poly._unify); a plain expression becomes a
# Poly in the same generator
.P_unify <- function(f, g) {
  if (!inherits(g, "dpoly")) g <- .P(as_sym(g), f$x)
  if (!inherits(f, "dpoly")) f <- .P(as_sym(f), g$x)
  K <- .dom_unify(f$K, g$K)
  list(.dp_convert(f, K), .dp_convert(g, K))
}
P_add <- function(f, g) { u <- .P_unify(f, g); .dp_add(u[[1L]], u[[2L]]) }
P_sub <- function(f, g) { u <- .P_unify(f, g); .dp_sub(u[[1L]], u[[2L]]) }
P_mul <- function(f, g) { u <- .P_unify(f, g); .dp_mul(u[[1L]], u[[2L]]) }
P_neg <- function(f) .dp_neg(f)
P_pow <- function(f, n) .dp_pow(f, n)
P_mul_ground <- function(f, c) .dp_mul_ground(f, as_sym(c))
P_quo_ground <- function(f, c) .dp_quo_ground(f, as_sym(c))
P_exquo_ground <- function(f, c) .dp_quo_ground(f, as_sym(c))
P_diff <- function(f) .dp_diff(f)
P_to_field <- function(f) .dp_to_field(f)
P_eval <- function(f, a) .dp_eval(f, as_sym(a))

# Poly.div with auto = True: rings go to their field, the results retract
.P_auto <- function(f, g, op) {
  u <- .P_unify(f, g); f <- u[[1L]]; g <- u[[2L]]
  retract <- FALSE
  if (!.dk_field(f$K)) { f <- .dp_to_field(f); g <- .dp_to_field(g); retract <- TRUE }
  r <- op(f, g)
  if (retract) r <- lapply(r, .dp_to_ring)
  r
}
P_div <- function(f, g) .P_auto(f, g, function(a, b) .dp_div(a, b))
P_quo <- function(f, g) P_div(f, g)[[1L]]
P_rem <- function(f, g) P_div(f, g)[[2L]]
P_exquo <- function(f, g) {
  qr <- P_div(f, g)
  if (!P_is_zero(qr[[2L]])) stop("ExactQuotientFailed", call. = FALSE)
  qr[[1L]]
}
P_gcd <- function(f, g) { u <- .P_unify(f, g); .dp_gcd(u[[1L]], u[[2L]]) }
P_cofactors <- function(f, g) { u <- .P_unify(f, g); .dp_inner_gcd(u[[1L]], u[[2L]]) }
P_lcm <- function(f, g) {
  u <- .P_unify(f, g); f <- u[[1L]]; g <- u[[2L]]
  if (.dp_zero(f) || .dp_zero(g)) return(.dp_set(f, list()))
  h <- .dp_quo(.dp_mul(f, g), .dp_gcd(f, g))
  if (.dk_field(f$K)) .dp_monic(h) else { c <- .dp_LC(h); if (.dom_is_negative(f$K, c)) .dp_neg(h) else h }
}
P_monic <- function(f) .dp_monic(if (.dk_field(f$K)) f else .dp_to_field(f))
P_primitive <- function(f) .dp_primitive(f)
P_half_gcdex <- function(f, g) .P_auto(f, g, function(a, b) .dp_half_gcdex(a, b))
P_gcdex <- function(f, g) .P_auto(f, g, function(a, b) .dp_gcdex(a, b))
P_invert <- function(f, g) .P_auto(f, g, function(a, b) list(.dp_invert(a, b)))[[1L]]

# Poly.cancel(g, include = True): list(p, q)
P_cancel <- function(f, g, include = TRUE) {
  u <- .P_unify(f, g); F <- u[[1L]]; G <- u[[2L]]
  K <- F$K
  if (.dk_field(K) && .dk(K) %in% c("QQ", "QQ_I", "FracField")) {
    cc <- .dp_cancel(F, G)
    if (!include) return(cc)
    nd <- .as_numer_denom(cc[[1L]])
    return(list(.dp_mul_ground(cc[[2L]], nd[[1L]]), .dp_mul_ground(cc[[3L]], nd[[2L]])))
  }
  ig <- .dp_inner_gcd(F, G)
  p <- ig[[2L]]; q <- ig[[3L]]
  cp <- S_One
  pn <- .dom_is_negative(K, .dp_LC(p)); qn <- .dom_is_negative(K, .dp_LC(q))
  if (pn && qn) { p <- .dp_neg(p); q <- .dp_neg(q) }
  else if (pn) { cp <- S_NegativeOne; p <- .dp_neg(p) }
  else if (qn) { cp <- S_NegativeOne; q <- .dp_neg(q) }
  if (!include) return(list(cp, p, q))
  list(.dp_mul_ground(p, cp), q)
}

# Poly.resultant(g, includePRS = True): list(res expression, list of polys)
P_resultant <- function(f, g, includePRS = FALSE) {
  u <- .P_unify(f, g)
  r <- .dp_prs_resultant(u[[1L]], u[[2L]])
  if (includePRS) list(r[[1L]], r[[2L]]) else r[[1L]]
}

# Poly.sqf_list and sqf_list_include
P_sqf_list <- function(f, all = FALSE) .dp_sqf_list(f, all)
P_sqf_list_include <- function(f, all = FALSE) {
  sl <- .dp_sqf_list(f, all)
  coeff <- sl[[1L]]; factors <- sl[[2L]]
  if (length(factors) && factors[[1L]][[2L]] == 1L)
    return(c(list(list(.dp_mul_ground(factors[[1L]][[1L]], coeff), 1L)), factors[-1L]))
  c(list(list(.dp_const(f, coeff), 1L)), factors)
}

# Poly from the monomial dict {k: c}
P_from_dict <- function(d, x, domain = NULL) {
  ks <- as.integer(names(d))
  if (!length(ks)) return(.P(S_Zero, x, domain))
  e <- .add(lapply(seq_along(d), function(i) .mul(list(d[[i]], .pow(x, .int(ks[i]))))))
  .P(e, x, domain)
}

format.dpoly <- function(x, ...) paste0("Poly(", .str(P_expr(x)), ", ", .str(x$x), ", domain='", format(x$K), "')")
print.dpoly <- function(x, ...) { cat(format(x), "\n"); invisible(x) }
