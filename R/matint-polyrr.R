# Poly over RR (sympy/polys/domains/realfield.py): an expression with Float
# coefficients keeps its domain RR, or RR[...] when other symbols enter the
# coefficients, and the coefficients come back as Floats. The rep stays
# over QQ, so only the conversion back to expressions changes.

# Mark a Poly built without domain options from an expression with Floats
.matint_poly_rr <- function(P, expr) {
  fl <- .float_atoms(expr)
  if (length(fl)) P$rr_prec <- max(vapply(fl, .float_prec, 0))
  P
}

# RR.to_sympy (or the PolyRing over RR) for the coefficient of each monomial
.matint_rr_coeffs <- function(p, u, keys) {
  ng <- .poly_ngens(p)
  rep <- p$rep
  extra <- p$gens[-seq_len(ng)]
  lapply(u, function(k) {
    w <- which(keys == k)
    sub <- qp(rep$e[w, -seq_len(ng), drop = FALSE], rep$c[w])
    if (!length(extra)) return(.qp_float_expr(qp_norm(sub$e, sub$c), list(), p$rr_prec))
    sub <- qp_norm(sub$e, sub$c)
    if (qp_is_zero(sub)) S_Zero else .qp_float_expr(sub, extra, p$rr_prec)
  })
}

# The zero of the ground domain: 0.0 over RR, 0 over a polynomial ring
.matint_rr_zero <- function(p) {
  if (is.null(p$rr_prec) || .poly_ngens(p) < length(p$gens)) return(S_Zero)
  .mp_float(.mp_from_q(gmp::as.bigq(0), p$rr_prec), p$rr_prec)
}

# Poly.as_expr() over RR: the RR coefficients times their monomials
.matint_rr_as_expr <- function(p) {
  t <- .poly_terms_expr(p)
  if (!length(t$coeffs)) return(S_Zero)
  .add(lapply(seq_along(t$coeffs), function(i) .mul(list(t$coeffs[[i]], .monomial_vec_expr(t$monoms[i, ], p$gens)))))
}
