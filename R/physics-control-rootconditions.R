# The Routh-Hurwitz and Schur conditions of a polynomial over its domain,
# after sympy/polys/rootconditions.py (Poly.hurwitz_conditions and
# Poly.schur_conditions), which the stability tests of the control module
# use. The coefficient domain decides the algorithm: fraction free over the
# integers, a polynomial ring or a fraction field (after clearing
# denominators), and the EXRAW variants otherwise.

# The kind of the coefficient domain of Poly(f, x): "ZZ", "QQ", "RR",
# "poly" (a polynomial ring), "frac" (a fraction field), "EX" or "complex".
.rc_domain <- function(coeffs) {
  if (any(vapply(coeffs, function(c) .has1(c, S_I), TRUE))) return("complex")
  types <- vapply(coeffs, function(c) c$type, "")
  if (all(types == "Integer")) return("ZZ")
  if (all(types %in% c("Integer", "Rational"))) return("QQ")
  if (all(types %in% c("Integer", "Rational", "Float"))) return("RR")
  gens <- unique(unlist(lapply(coeffs, function(c) vapply(free_symbols(c), function(s) s$key, ""))))
  if (!length(gens)) return("EX")
  is_poly <- function(e) {
    if (e$type %in% c("Integer", "Rational", "Float", "Symbol", "Dummy")) return(TRUE)
    if (e$type %in% c("Add", "Mul")) return(all(vapply(e$args, is_poly, TRUE)))
    if (e$type == "Pow") return(e$args[[2L]]$type == "Integer" && isTRUE(e$args[[2L]]$value >= 0) && is_poly(e$args[[1L]]))
    FALSE
  }
  nds <- lapply(coeffs, function(c) .as_numer_denom(c))
  if (!all(vapply(nds, function(nd) is_poly(nd[[1L]]) && is_poly(nd[[2L]]), TRUE))) return("EX")
  if (all(vapply(nds, function(nd) !length(free_symbols(nd[[2L]])), TRUE))) return("poly")
  "frac"
}

.rc_float_to_rational <- function(c) {
  if (c$type != "Float") return(c)
  .limit_den_bigq(.py_float_to_ratio(as.numeric(N(c))), max(gmp::as.bigz(2)^53 %/% 200, 99))
}

# _dup_routh_hurwitz_fraction_free over a ring given by its operations
.rc_fraction_free <- function(p, R) {
  if (length(p) == 1L) return(list())
  if (length(p) == 2L) return(list(R$mul(p[[1L]], p[[2L]])))
  if (length(p) == 3L) return(list(R$mul(p[[1L]], p[[2L]]), R$mul(p[[1L]], p[[3L]])))
  LC <- p[[1L]]; TC <- p[[length(p)]]
  monic <- R$is_one(LC)
  is_LC_negative <- R$is_negative(LC)
  if (!is.na(is_LC_negative)) {
    for (coeff in list(TC, p[[2L]])) {
      neg <- R$is_negative(coeff)
      if (!is.na(neg) && (R$is_zero(coeff) || xor(is_LC_negative, neg))) return(list(S_NegativeOne))
    }
  }
  p1s <- list(p[[2L]])
  while (length(p) > 3L) {
    qs <- lapply(p[-1L], function(qi) R$mul(p[[2L]], qi))
    if (length(qs) > 2L) for (i in seq(1L, length(qs) - 2L, by = 2L)) qs[[i + 1L]] <- R$sub(qs[[i + 1L]], R$mul(p[[i + 3L]], p[[1L]]))
    p <- qs
    p1 <- p[[2L]]
    if (R$is_zero(p1)) return(list(S_NegativeOne))
    n1 <- length(p1s)
    if (n1 >= 2L) p1 <- R$exquo(p1, p1s[[n1 - 1L]])
    if (n1 >= 3L) {
      p1 <- R$exquo(p1, p1s[[n1 - 2L]])
      p1 <- R$exquo(p1, p1s[[n1 - 2L]])
      p <- lapply(p, function(c) R$exquo(c, p1s[[n1 - 2L]]))
      p <- lapply(p, function(c) R$exquo(c, p1s[[n1 - 2L]]))
    }
    negp1 <- R$is_negative(p1)
    if (!is.na(negp1)) {
      if (length(p1s) %% 2L == 0L) {
        if (!is.na(is_LC_negative) && xor(is_LC_negative, negp1)) return(list(S_NegativeOne))
      } else if (negp1) {
        return(list(S_NegativeOne))
      }
    }
    p1s[[length(p1s) + 1L]] <- p1
  }
  if (!monic) for (i in seq(1L, length(p1s), by = 2L)) p1s[[i]] <- R$mul(p1s[[i]], LC)
  p1s[[length(p1s) + 1L]] <- R$mul(LC, TC)
  p1s
}

.rc_ring_zz <- list(
  mul = function(a, b) .mul(list(a, b)), sub = function(a, b) .sub(a, b),
  exquo = function(a, b) .div(a, b),
  is_zero = function(a) is_zero_num(a), is_one = function(a) identical(a$key, S_One$key),
  is_negative = function(a) isTRUE(a < 0))

.rc_ring_poly <- list(
  mul = function(a, b) expand(.mul(list(a, b))), sub = function(a, b) expand(.sub(a, b)),
  exquo = function(a, b) expand(cancel(.div(a, b))),
  is_zero = function(a) is_zero_num(expand(a)), is_one = function(a) identical(expand(a)$key, S_One$key),
  is_negative = function(a) if (length(free_symbols(a))) NA else isTRUE(a < 0))

# dup_clear_denoms over the integers or a polynomial ring
.rc_clear_denoms <- function(p) {
  dens <- lapply(p, function(c) .as_numer_denom(c)[[2L]])
  common <- Reduce(function(a, b) lcm(a, b), dens, S_One)
  lapply(p, function(c) expand(cancel(.mul(list(c, common)))))
}

.rc_routh_hurwitz <- function(p, kind) {
  if (!length(p)) stop("zero polynomial")
  if (kind == "complex") .physics_error("NotImplementedError", "Routh-Hurwitz is not implemented for complex domains")
  if (kind == "RR") {
    pz <- .rc_clear_denoms(lapply(p, .rc_float_to_rational))
    return(lapply(.rc_fraction_free(pz, .rc_ring_zz), function(c) Float(c)))
  }
  if (kind == "QQ") return(.rc_fraction_free(.rc_clear_denoms(p), .rc_ring_zz))
  if (kind == "ZZ") return(.rc_fraction_free(p, .rc_ring_zz))
  if (kind == "poly") return(.rc_fraction_free(lapply(p, expand), .rc_ring_poly))
  if (kind == "frac") return(.rc_fraction_free(.rc_clear_denoms(p), .rc_ring_poly))
  .rc_exraw(p)
}

.rc_exraw <- function(p) {
  if (all(vapply(p, function(c) c$type %in% c("Integer", "Rational", "Float"), TRUE))) return(.rc_exraw_div(p))
  .rc_exraw_no_div(p, list(), list())
}

.rc_exraw_div <- function(p) {
  if (length(p) < 2L) return(list())
  if (isTRUE(ask_fact(.mul(list(p[[1L]], p[[2L]])), "nonpositive"))) return(list(S_NegativeOne))
  qs <- p
  for (i in seq(2L, length(p), by = 2L)) qs[[i - 1L]] <- .sub(qs[[i - 1L]], .div(.mul(list(p[[i]], p[[1L]])), p[[2L]]))
  qs <- qs[-1L]
  c(list(.mul(list(p[[1L]], p[[2L]]))), .rc_exraw_div(qs))
}

.rc_exraw_no_div <- function(p, previous_cond, nonzeros) {
  if (length(p) < 2L) return(lapply(nonzeros, function(z) .pow(z, S_Two)))
  qs <- lapply(p, function(qi) .mul(list(p[[2L]], qi)))
  for (i in seq(2L, length(p), by = 2L)) qs[[i - 1L]] <- .sub(qs[[i - 1L]], .mul(list(p[[i]], p[[1L]])))
  qs <- qs[-1L]
  r <- .rc_clear_cond(.mul(list(p[[1L]], p[[2L]])), previous_cond, nonzeros)
  if (isTRUE(ask_fact(r$cond, "nonpositive"))) return(list(S_NegativeOne))
  c(list(r$cond), .rc_exraw_no_div(qs, r$previous_cond, r$nonzeros))
}

.rc_set_union <- function(a, b) {
  ka <- vapply(a, .k, "")
  c(a, b[!vapply(b, .k, "") %in% ka])
}

.rc_clear_cond <- function(cond, previous_cond, nonzeros) {
  f1 <- .rc_build_factors(cond, previous_cond)
  nonzeros <- .rc_set_union(nonzeros, f1$nonzeros)
  new_cond <- .mul(c(list(as_sym(f1$sign)), f1$factors))
  f2 <- .rc_build_factors(factor_terms(new_cond), previous_cond)
  nonzeros <- .rc_set_union(nonzeros, f2$nonzeros)
  previous_cond[[length(previous_cond) + 1L]] <- list(sign = f1$sign, factors = f1$factors)
  list(cond = .mul(c(list(as_sym(f2$sign)), f2$factors)), previous_cond = previous_cond, nonzeros = nonzeros)
}

.rc_build_factors <- function(cond, previous_cond) {
  cond <- signsimp(cond)
  sign <- 1L
  if (could_extract_minus_sign(cond)) { cond <- .neg(cond); sign <- -1L }
  factors <- list(); nonzeros <- list()
  pd <- .rc_powers_dict(cond)
  for (i in seq_along(pd$bases)) {
    b <- pd$bases[[i]]; e <- pd$exps[[i]]
    if (e$type != "Integer") factors[[length(factors) + 1L]] <- .pow(b, e)
    else if (gmp::numerator(e$value) %% 2L != 0L) factors[[length(factors) + 1L]] <- b
    else nonzeros[[length(nonzeros) + 1L]] <- b
  }
  fk <- function(l) vapply(l, .k, "")
  for (pc in previous_cond) {
    pk <- fk(pc$factors)
    if (all(pk %in% fk(factors))) {
      sign <- sign * pc$sign
      factors <- factors[!fk(factors) %in% pk]
    }
    nonzeros <- nonzeros[!fk(nonzeros) %in% pk]
  }
  list(sign = sign, factors = factors, nonzeros = nonzeros)
}

# Expr.as_powers_dict of a product
.rc_powers_dict <- function(e) {
  fs <- if (e$type == "Mul") e$args else list(e)
  bases <- list(); exps <- list()
  for (f in fs) {
    be <- if (f$type == "Pow") list(f$args[[1L]], f$args[[2L]]) else list(f, S_One)
    k <- match(be[[1L]]$key, vapply(bases, .k, ""))
    if (is.na(k)) { bases[[length(bases) + 1L]] <- be[[1L]]; exps[[length(exps) + 1L]] <- be[[2L]] }
    else exps[[k]] <- .add(list(exps[[k]], be[[2L]]))
  }
  list(bases = bases, exps = exps)
}

.rc_coeffs <- function(f, x, fast, kind = NULL) {
  cs <- all_coeffs(Poly(f, x))
  list(coeffs = cs, kind = if (isTRUE(fast)) "EX" else kind %||% .rc_domain(cs))
}

# Poly(f, x).hurwitz_conditions(); `kind` overrides the domain of the
# coefficients
.lti_hurwitz <- function(f, x, fast = FALSE, kind = NULL) {
  v <- .rc_coeffs(f, x, fast, kind)
  .rc_routh_hurwitz(v$coeffs, v$kind)
}

# Poly(f, x).schur_conditions()
.lti_schur <- function(f, x, fast = FALSE, kind = NULL) {
  v <- .rc_coeffs(f, x, fast, kind)
  if (!length(v$coeffs)) stop("zero polynomial")
  if (v$kind == "complex") .physics_error("NotImplementedError", "Schur conditions is not implemented for complex domains")
  cs <- v$coeffs
  kind <- v$kind
  if (kind == "RR") { cs <- lapply(cs, .rc_float_to_rational); kind <- "QQ" }
  fm1 <- .add(lapply(seq_along(cs), function(i) .mul(list(cs[[i]], .pow(S_NegativeOne, .int(length(cs) - i))))))
  if (is_zero_num(expand(fm1))) return(list(S_NegativeOne))
  n <- length(cs) - 1L
  ft <- .add(lapply(seq_along(cs), function(i) .mul(list(cs[[i]], .pow(.add(list(x, S_One)), .int(n - i + 1L)),
                                                     .pow(.sub(S_One, x), .int(i - 1L))))))
  ft <- all_coeffs(Poly(expand(ft), x))
  if (kind %in% c("poly", "frac")) ft <- lapply(ft, function(c) expand(cancel(c)))
  r <- .rc_routh_hurwitz(ft, kind)
  if (v$kind == "RR") r <- lapply(r, function(c) Float(c))
  r
}
