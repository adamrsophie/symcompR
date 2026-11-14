# Elements of the polynomial domains as values, after sympy/polys/domains
# (the dtype of each domain), sympy/polys/rings.py (PolyElement) and
# sympy/polys/fields.py (FracElement). These are the entries of
# DomainMatrix (R/polys-domainmatrix.R) and the coefficients of power series
# (R/polys-series.R). Elements are: bigz for ZZ, bigq for QQ, residues
# (sym_modint) for GF(p), Gaussian pairs for ZZ_I and QQ_I, doubles and
# complex numbers for RR and CC, wrapped expressions for EX and EXRAW, and
# PolyElement and FracElement objects for K[X] and K(X) over ZZ or QQ, whose
# terms are kept as a qpoly in lex order (R/polys-core.R).
# Limitation: PolyElement and FracElement take ZZ or QQ coefficients only;
# other grounds raise NotImplementedError.

.de_kind <- function(K) {
  d <- .dom(K)
  if (d$kind == "EX" && d$name == "EXRAW") "EXRAW" else d$kind
}
.de_is_field <- function(K) .de_kind(K) %in% c("QQ", "RR", "CC", "QQ_I", "FracField", "EX", "EXRAW", "GF", "Algebraic")
.de_has_assoc_ring <- function(K) .de_kind(K) %in% c("QQ", "QQ_I", "FracField")
.de_has_assoc_field <- function(K) .de_kind(K) %in% c("ZZ", "ZZ_I", "PolyRing")
.de_same <- function(K0, K1) identical(format(K0), format(K1))

.de_get_field <- function(K) {
  switch(.de_kind(K), ZZ = QQ, ZZ_I = QQ_I,
         PolyRing = do.call(frac_field, c(list(.dk_ground(K)), .dk_gens(K))),
         { if (!.de_is_field(K)) stop(DomainError(paste0("there is no field associated with ", format(K)))); K })
}
.de_get_ring <- function(K) {
  switch(.de_kind(K), QQ = ZZ, QQ_I = ZZ_I,
         FracField = do.call(poly_ring, c(list(.dk_ground(K)), .dk_gens(K))),
         { if (.de_is_field(K)) stop(DomainError(paste0("there is no ring associated with ", format(K)))); K })
}

# ---------------------------------------------------------------------------
# EX elements

.exel <- function(e, raw = FALSE) structure(list(ex = as_sym(e), raw = raw), class = "EXElement")
# Expression.simplify: cancel then expand
.exel_s <- function(e, raw) if (raw) .exel(e, TRUE) else .exel(.ex_simp(e), FALSE)
#' @export
Ops.EXElement <- function(e1, e2) {
  raw <- if (inherits(e1, "EXElement")) e1$raw else e2$raw
  ex <- function(z) if (inherits(z, "EXElement")) z$ex else .de_value_sym(z)
  if (missing(e2)) return(if (.Generic == "-") .exel(.neg(e1$ex), raw) else e1)
  a <- ex(e1); b <- ex(e2)
  z <- function(v) identical(v$key, S_Zero$key)
  switch(.Generic,
    `+` = if (!raw && z(b)) .exel(a, raw) else if (!raw && z(a)) .exel(b, raw) else .exel_s(.add(list(a, b)), raw),
    `-` = if (!raw && z(b)) .exel(a, raw) else if (!raw && z(a)) .exel(.neg(b), raw) else .exel_s(.sub(a, b), raw),
    `*` = if (!raw && (z(a) || z(b))) .exel(S_Zero, raw) else if (!raw && is_number_type(a) && is_number_type(b)) .exel(.mul(list(a, b)), raw)
          else .exel_s(.mul(list(a, b)), raw),
    `/` = .exel_s(.div(a, b), raw), `^` = .exel_s(.pow(a, b), raw),
    `==` = identical(a$key, b$key), `!=` = !identical(a$key, b$key),
    stop("operator ", .Generic, " is not defined for EX elements", call. = FALSE))
}
#' @export
chooseOpsMethod.EXElement <- function(x, y, mx, my, cl, reverse) TRUE
#' @export
format.EXElement <- function(x, ...) paste0("EX(", sstr(x$ex), ")")
#' @export
print.EXElement <- function(x, ...) { cat(format(x), "\n", sep = ""); invisible(x) }

.de_value_sym <- function(z) {
  if (is_sym(z)) return(z)
  if (inherits(z, "bigz") || inherits(z, "bigq")) return(make_rational(gmp::as.bigq(z)))
  as_sym(z)
}

# ---------------------------------------------------------------------------
# PolyElement: an element of K[X] with K = ZZ or QQ

.pe <- function(p, K) structure(list(p = p, K = K), class = "PolyElement")
.pe_check_ground <- function(K) {
  g <- .dk(.dk_ground(K))
  if (!(g %in% c("ZZ", "QQ"))) stop(.sympy_condition("NotImplementedError", paste0("polynomial rings over ", g, " are not supported")))
  g
}
.pe_n <- function(K) length(.dk_gens(K))
.pe_from_sym <- function(K, e) {
  g <- .pe_check_ground(K)
  e <- as_sym(e)
  p <- expr_to_qpoly(e, .dk_gens(K), strict = TRUE)
  if (is.null(p) || (g == "ZZ" && any(gmp::denominator(p$c) != 1)))
    stop(CoercionFailed(paste0("expected an element of ", format(K), ", got ", sstr(e))))
  .pe(p, K)
}
.pe_const <- function(K, v) .pe(qp_const(gmp::as.bigq(v), .pe_n(K)), K)
.de_pe_sym <- function(a) qpoly_to_expr(a$p, .dk_gens(a$K))
.pe_ground_ok <- function(K, p) .dk(.dk_ground(K)) != "ZZ" || all(gmp::denominator(p$c) == 1)

#' @export
Ops.PolyElement <- function(e1, e2) {
  K <- if (inherits(e1, "PolyElement")) e1$K else e2$K
  P <- function(z) if (inherits(z, "PolyElement")) z$p else .de_convert(K, z)$p
  if (missing(e2)) return(if (.Generic == "-") .pe(qp_neg(e1$p), K) else e1)
  if (.Generic == "^") return(.pe(qp_pow(e1$p, as.integer(as.character(e2))), K))
  a <- P(e1); b <- P(e2)
  switch(.Generic,
    `+` = .pe(qp_add(a, b), K), `-` = .pe(qp_sub(a, b), K), `*` = .pe(qp_mul(a, b), K),
    `/` = .pe_exquo(.pe(a, K), .pe(b, K)),
    `==` = qp_equal(a, b), `!=` = !qp_equal(a, b),
    stop("operator ", .Generic, " is not defined for polynomial ring elements", call. = FALSE))
}
#' @export
chooseOpsMethod.PolyElement <- function(x, y, mx, my, cl, reverse) TRUE

.pe_exquo <- function(a, b) {
  if (qp_is_zero(b$p)) stop(.sympy_condition("ZeroDivisionError", "polynomial division"))
  q <- qp_exquo(a$p, b$p)
  if (is.null(q) || !.pe_ground_ok(a$K, q)) stop(ExactQuotientFailed(paste0(format(b), " does not divide ", format(a))))
  .pe(q, a$K)
}

# PolyElement.str
.de_pe_str <- function(a) {
  p <- a$p; K <- a$K
  if (qp_is_zero(p)) return("0")
  gens <- vapply(.dk_gens(K), function(g) { s <- .str(g); if (is_symbol_like(g)) s else paste0("(", s, ")") }, "")
  parts <- character()
  for (i in seq_along(p$c)) {
    c1 <- p$c[i]; ex <- p$e[i, ]
    neg <- c1 < 0
    parts <- c(parts, if (neg) " - " else " + ")
    if (neg) c1 <- -c1
    if (all(ex == 0L)) sc <- .q_str(c1)
    else sc <- if (c1 == 1) "" else .q_str(c1)
    sx <- character()
    for (j in seq_along(ex)) if (ex[j]) sx <- c(sx, if (ex[j] != 1L) paste0(gens[j], "**", ex[j]) else gens[j])
    if (nzchar(sc)) sx <- c(sc, sx)
    parts <- c(parts, paste(sx, collapse = "*"))
  }
  head <- parts[1L]; parts <- parts[-1L]
  if (head == " - ") parts <- c("-", parts)
  paste(parts, collapse = "")
}
.q_str <- function(q) {
  q <- gmp::as.bigq(q)
  if (gmp::denominator(q) == 1) as.character(gmp::numerator(q)) else paste0(gmp::numerator(q), "/", gmp::denominator(q))
}
#' @export
format.PolyElement <- function(x, ...) gsub("**", "^", .de_pe_str(x), fixed = TRUE)
#' @export
print.PolyElement <- function(x, ...) { cat(format(x), "\n", sep = ""); invisible(x) }
# precedence of a PolyElement (precedence_PolyElement)
.pe_prec <- function(a) {
  p <- a$p
  if (length(p$c) == 1L && p$c[1L] == 1 && sum(p$e[1L, ]) == 1L) return(1000L)
  if (qp_is_const(p)) return(if (length(p$c) && (p$c[1L] < 0 || gmp::denominator(p$c[1L]) != 1)) 1000L else 1000L)
  if (length(p$c) == 1L) return(50L)
  40L
}

# ---------------------------------------------------------------------------
# FracElement: an element of K(X), numerator and denominator in K[X]

.fe <- function(n, d, K) structure(list(n = n, d = d, K = K), class = "FracElement")
.fe_ring <- function(K) .de_get_ring(K)
# FracField.new: PolyElement.cancel of the numerator and denominator
.fe_new <- function(n, d, K) {
  if (qp_is_zero(d)) stop(.sympy_condition("ZeroDivisionError", "rational function division by zero"))
  g <- .pe_check_ground(K)
  if (qp_is_zero(n)) return(.fe(n, qp_const(1, ncol(n$e)), K))
  pq <- .pe_cancel(n, d, g == "QQ")
  .fe(pq[[1L]], pq[[2L]], K)
}
# PolyElement.cancel over ZZ (field = FALSE) or QQ (field = TRUE)
.pe_cancel <- function(p, q, field) {
  if (!field) {
    h <- qp_gcd(p, q)
    hp <- qp_primitive(h)$poly
    if (qp_lc(hp) < 0) hp <- qp_neg(hp)
    k <- gmp::gcd.bigz(.bq_int_content(p), .bq_int_content(q))
    hh <- qp_scale(hp, gmp::as.bigq(k))
    p <- qp_exquo(p, hh); q <- qp_exquo(q, hh)
  } else {
    r <- .qp_cancel(p, q)
    p <- r[[1L]]; q <- r[[2L]]
  }
  if (qp_lc(q) < 0) { p <- qp_neg(p); q <- qp_neg(q) }
  list(p, q)
}
.fe_from_sym <- function(K, e) {
  g <- .pe_check_ground(K)
  nd <- .as_numer_denom(together(as_sym(e)))
  gens <- .dk_gens(K)
  n <- expr_to_qpoly(nd[[1L]], gens, strict = TRUE); d <- expr_to_qpoly(nd[[2L]], gens, strict = TRUE)
  if (is.null(n) || is.null(d)) stop(CoercionFailed(paste0("expected an element of ", format(K), ", got ", sstr(e))))
  if (g == "ZZ") {
    # clear the rational denominators into the fraction
    l1 <- .bq_lcm_den(n); l2 <- .bq_lcm_den(d)
    n <- qp_scale(n, gmp::as.bigq(l1 * l2)); d <- qp_scale(d, gmp::as.bigq(l1 * l2))
  }
  .fe_new(n, d, K)
}
.fe_to_sym <- function(a) .div(qpoly_to_expr(a$n, .dk_gens(a$K)), qpoly_to_expr(a$d, .dk_gens(a$K)))

#' @export
Ops.FracElement <- function(e1, e2) {
  K <- if (inherits(e1, "FracElement")) e1$K else e2$K
  F <- function(z) if (inherits(z, "FracElement")) z else .de_convert(K, z)
  if (missing(e2)) return(if (.Generic == "-") .fe(qp_neg(e1$n), e1$d, K) else e1)
  if (.Generic == "^") {
    k <- as.integer(as.character(e2))
    if (k >= 0L) return(.fe(qp_pow(e1$n, k), qp_pow(e1$d, k), K))
    if (qp_is_zero(e1$n)) stop(.sympy_condition("ZeroDivisionError", "rational function division by zero"))
    return(.fe_new(qp_pow(e1$d, -k), qp_pow(e1$n, -k), K))
  }
  a <- F(e1); b <- F(e2)
  switch(.Generic,
    `+` = .fe_new(qp_add(qp_mul(a$n, b$d), qp_mul(b$n, a$d)), qp_mul(a$d, b$d), K),
    `-` = .fe_new(qp_sub(qp_mul(a$n, b$d), qp_mul(b$n, a$d)), qp_mul(a$d, b$d), K),
    `*` = .fe_new(qp_mul(a$n, b$n), qp_mul(a$d, b$d), K),
    `/` = .fe_new(qp_mul(a$n, b$d), qp_mul(a$d, b$n), K),
    `==` = qp_equal(a$n, b$n) && qp_equal(a$d, b$d), `!=` = !(qp_equal(a$n, b$n) && qp_equal(a$d, b$d)),
    stop("operator ", .Generic, " is not defined for rational function field elements", call. = FALSE))
}
#' @export
chooseOpsMethod.FracElement <- function(x, y, mx, my, cl, reverse) TRUE

# StrPrinter._print_FracElement
.fe_str <- function(a) {
  R <- .fe_ring(a$K)
  num <- .pe(a$n, R); den <- .pe(a$d, R)
  if (qp_equal(a$d, qp_const(1, ncol(a$d$e)))) return(.de_pe_str(num))
  ns <- .de_pe_str(num); ds <- .de_pe_str(den)
  if (.pe_prec(num) < 50L) ns <- paste0("(", ns, ")")
  if (.pe_prec(den) < 1000L) ds <- paste0("(", ds, ")")
  paste0(ns, "/", ds)
}
#' @export
format.FracElement <- function(x, ...) gsub("**", "^", .fe_str(x), fixed = TRUE)
#' @export
print.FracElement <- function(x, ...) { cat(format(x), "\n", sep = ""); invisible(x) }

# ---------------------------------------------------------------------------
# Domain operations on elements

.de_zero <- function(K) .de_convert(K, 0L)
.de_one <- function(K) .de_convert(K, 1L)

# K.convert(a): a number, an expression or an element of another domain
.de_convert <- function(K, a, K0 = NULL) {
  k <- .de_kind(K)
  if (inherits(a, "DomainScalar")) a <- a$element
  # elements already of K
  if (k == "ZZ" && inherits(a, "bigz")) return(a)
  if (k == "QQ" && inherits(a, "bigq")) return(a)
  if (k == "PolyRing" && inherits(a, "PolyElement") && .de_same(a$K, K)) return(a)
  if (k == "FracField" && inherits(a, "FracElement") && .de_same(a$K, K)) return(a)
  if (k %in% c("EX", "EXRAW") && inherits(a, "EXElement")) return(if (k == "EXRAW") .exel(a$ex, TRUE) else a)
  if (k == "GF" && inherits(a, "sym_modint")) return(a)
  # elements of other domains go through their value
  if (inherits(a, "PolyElement")) {
    if (k == "FracField") return(.fe_new(a$p, qp_const(1, ncol(a$p$e)), K))
    if (k == "PolyRing") return(.pe_from_sym(K, .de_pe_sym(a)))
    a <- .de_pe_sym(a)
  } else if (inherits(a, "FracElement")) {
    a <- .fe_to_sym(a)
  } else if (inherits(a, "EXElement")) {
    a <- a$ex
  } else if (inherits(a, "sym_modint")) {
    a <- a$val
  } else if (inherits(a, "sym_gaussian")) {
    a <- .gauss_sym(a)
  } else if (inherits(a, "ANP")) {
    a <- .anp_to_sym(a)
  }
  switch(EXPR = k,
    ZZ = {
      if (is.numeric(a) && !is.integer(a) && a != round(a)) stop(CoercionFailed(paste0("expected an integer, got ", a)))
      if (is_sym(a) && a$type != "Integer") stop(CoercionFailed(paste0("expected an integer, got ", sstr(a))))
      .dom_convert(.dom(ZZ), a)
    },
    QQ = {
      if (is_sym(a) && !(a$type %in% c("Integer", "Rational"))) stop(CoercionFailed(paste0("expected `Rational` object, got ", sstr(a))))
      gmp::as.bigq(.dom_convert(.dom(QQ), a))
    },
    GF = .dom_convert(.dom(K), a),
    RR = , CC = .dom_convert(.dom(K), a),
    ZZ_I = , QQ_I = .dom_convert(.dom(K), a),
    EX = .exel(.de_value_sym(a)),
    EXRAW = .exel(.de_value_sym(a), TRUE),
    PolyRing = if (is_sym(a)) .pe_from_sym(K, a) else .pe_const(K, .as_q(a)),
    FracField = if (is_sym(a)) .fe_from_sym(K, a) else { n <- .pe_n(K); .fe_new(qp_const(.as_q(a), n), qp_const(1, n), K) },
    Algebraic = .anp_from_sym(K, .de_value_sym(a)),
    stop(CoercionFailed(paste0("cannot convert to ", format(K)))))
}

.de_to_sympy <- function(K, a) {
  if (inherits(a, "PolyElement")) return(.de_pe_sym(a))
  if (inherits(a, "FracElement")) return(.fe_to_sym(a))
  if (inherits(a, "EXElement")) return(a$ex)
  if (inherits(a, "ANP")) return(.anp_to_sym(a))
  if (is.numeric(a) && !is.integer(a) && .de_kind(K) == "RR") return(Float(a))
  to_sympy(K, a)
}

.de_is_zero <- function(a) {
  if (inherits(a, "bigz") || inherits(a, "bigq")) return(a == 0)
  if (inherits(a, "PolyElement")) return(qp_is_zero(a$p))
  if (inherits(a, "FracElement")) return(qp_is_zero(a$n))
  if (inherits(a, "EXElement")) return(identical(a$ex$key, S_Zero$key))
  if (inherits(a, "sym_modint")) return(a$val == 0)
  if (inherits(a, "sym_gaussian")) return(a$re == 0 && a$im == 0)
  if (inherits(a, "ANP")) return(!length(a$rep))
  if (is_sym(a)) return(is_zero_num(a))
  a == 0
}
.de_nz <- function(a) !.de_is_zero(a)
.de_is_one <- function(K, a) .de_eq(a, .de_one(K))
.de_eq <- function(a, b) {
  if (inherits(a, "PolyElement") && inherits(b, "PolyElement")) return(qp_equal(a$p, b$p))
  if (inherits(a, "FracElement") && inherits(b, "FracElement")) return(qp_equal(a$n, b$n) && qp_equal(a$d, b$d))
  if (inherits(a, "EXElement") && inherits(b, "EXElement")) return(identical(a$ex$key, b$ex$key))
  if (inherits(a, "ANP") && inherits(b, "ANP")) return(length(a$rep) == length(b$rep) && all(a$rep == b$rep))
  isTRUE(a == b)
}

# K.exquo, K.quo, K.rem, K.div and K.inv
.de_exquo <- function(K, a, b) {
  k <- .de_kind(K)
  if (k == "ZZ") {
    if (b == 0) stop(.sympy_condition("ZeroDivisionError", "integer modulo by zero"))
    if (a %% b != 0) stop(ExactQuotientFailed(paste0(.de_str(b), " does not divide ", .de_str(a), " in ZZ")))
    return(a %/% b)
  }
  if (k == "PolyRing") return(.pe_exquo(a, b))
  if (k == "ZZ_I") {
    q <- a / b
    if (gmp::denominator(q$re) != 1 || gmp::denominator(q$im) != 1)
      stop(ExactQuotientFailed(paste0(format(b), " does not divide ", format(a), " in ZZ_I")))
    return(q)
  }
  if (.de_is_zero(b)) stop(.sympy_condition("ZeroDivisionError", "division by zero"))
  a / b
}
.de_quo <- function(K, a, b) {
  k <- .de_kind(K)
  if (k == "ZZ") return(a %/% b)
  if (k == "PolyRing") {
    pr <- qp_div_multi(a$p, list(b$p))
    return(.pe(pr$q[[1L]], K))
  }
  if (.de_is_zero(b)) stop(.sympy_condition("ZeroDivisionError", "division by zero"))
  a / b
}
.de_rem <- function(K, a, b) {
  k <- .de_kind(K)
  if (k == "ZZ") return(a %% b)
  if (.de_is_field(K)) return(.de_zero(K))
  stop(.sympy_condition("NotImplementedError", "rem"))
}
.de_inv <- function(K, a) .de_exquo(K, .de_one(K), a)

.de_is_negative <- function(K, a) {
  switch(EXPR = .de_kind(K),
    ZZ = , QQ = , RR = a < 0,
    PolyRing = !qp_is_zero(a$p) && qp_lc(a$p) < 0,
    FracField = !qp_is_zero(a$n) && qp_lc(a$n) < 0,
    EX = , EXRAW = { cm <- as_coeff_mul_split(a$ex)[[1L]]; is_number_type(cm) && isTRUE(num_sign(cm) < 0) },
    FALSE)
}
.de_abs <- function(K, a) if (.de_is_negative(K, a)) -a else a

# K.gcd, K.lcm and K.cofactors
.de_gcd <- function(K, a, b) {
  switch(.de_kind(K),
    ZZ = gmp::gcd.bigz(a, b),
    QQ = {
      if (a == 0) return(abs(b)); if (b == 0) return(abs(a))
      gmp::as.bigq(gmp::gcd.bigz(gmp::numerator(a), gmp::numerator(b)), gmp::lcm.bigz(gmp::denominator(a), gmp::denominator(b)))
    },
    PolyRing = {
      h <- qp_gcd(a$p, b$p)
      if (.dk(.dk_ground(K)) == "ZZ") {
        hp <- qp_primitive(h)$poly
        if (!qp_is_zero(hp) && qp_lc(hp) < 0) hp <- qp_neg(hp)
        ca <- if (qp_is_zero(a$p)) gmp::as.bigz(0) else .bq_int_content(a$p)
        cb <- if (qp_is_zero(b$p)) gmp::as.bigz(0) else .bq_int_content(b$p)
        h <- qp_scale(hp, gmp::as.bigq(gmp::gcd.bigz(ca, cb)))
      } else h <- qp_monic(h)
      .pe(h, K)
    },
    if (.de_is_field(K)) .de_one(K) else stop(.sympy_condition("NotImplementedError", "gcd")))
}
.de_lcm <- function(K, a, b) {
  switch(.de_kind(K),
    ZZ = gmp::lcm.bigz(a, b),
    QQ = gmp::as.bigq(gmp::lcm.bigz(gmp::numerator(a), gmp::numerator(b)), gmp::gcd.bigz(gmp::denominator(a), gmp::denominator(b))),
    PolyRing = {
      if (.de_is_zero(a) || .de_is_zero(b)) return(.de_zero(K))
      r <- .de_exquo(K, a * b, .de_gcd(K, a, b))
      if (.de_is_negative(K, r)) -r else r
    },
    a * b)
}
.de_cofactors <- function(K, a, b) {
  h <- .de_gcd(K, a, b)
  if (.de_is_zero(h)) return(list(h, .de_zero(K), .de_zero(K)))
  list(h, .de_exquo(K, a, h), .de_exquo(K, b, h))
}
.de_numer <- function(K, a) {
  switch(.de_kind(K), QQ = gmp::as.bigz(gmp::numerator(a)), FracField = .pe(a$n, .de_get_ring(K)), a)
}
.de_denom <- function(K, a) {
  switch(.de_kind(K), QQ = gmp::as.bigz(gmp::denominator(a)), FracField = .pe(a$d, .de_get_ring(K)), .de_one(K))
}

# str() and repr() of an element, as Python prints them
.de_str <- function(a) {
  if (inherits(a, "bigz")) return(as.character(a))
  if (inherits(a, "bigq")) return(.q_str(a))
  if (inherits(a, "PolyElement")) return(.de_pe_str(a))
  if (inherits(a, "FracElement")) return(.fe_str(a))
  if (inherits(a, "EXElement")) return(if (a$raw) gsub("^", "**", sstr(a$ex), fixed = TRUE) else paste0("EX(", gsub("^", "**", sstr(a$ex), fixed = TRUE), ")"))
  if (inherits(a, "sym_modint")) return(paste0(a$val, " mod ", a$mod))
  if (inherits(a, "sym_gaussian")) return(gsub("^", "**", .str(.gauss_sym(a)), fixed = TRUE))
  if (inherits(a, "ANP")) return(.anp_repr(a))
  if (is.logical(a)) return(if (isTRUE(a)) "True" else "False")
  if (is.double(a)) return(.py_mpf_str(a))
  if (is.complex(a)) return(paste0(.py_mpf_str(Re(a)), if (Im(a) < 0) " - " else " + ", .py_mpf_str(abs(Im(a))), "*I"))
  if (is_sym(a)) return(gsub("^", "**", sstr(a), fixed = TRUE))
  as.character(a)
}
.de_repr <- function(a) {
  if (inherits(a, "bigq")) return(paste0("MPQ(", gmp::numerator(a), ",", gmp::denominator(a), ")"))
  if (inherits(a, "sym_modint")) return(paste0("SymmetricModularIntegerMod", a$mod, "(", a$val, ")"))
  if (inherits(a, "sym_gaussian")) return(paste0(a$kind, "(", .q_str(a$re), ", ", .q_str(a$im), ")"))
  if (is.double(a)) return(paste0("mpf('", .py_mpf_repr(a), "')"))
  .de_str(a)
}
# mpmath's str (15 digits) and repr (17 digits) of a double
.py_mpf_str <- function(x, digits = 15L) {
  if (x == 0) return("0.0")
  s <- formatC(x, digits = digits, format = "g", flag = "")
  s <- trimws(s)
  if (grepl("e", s)) {
    m <- sub("e.*$", "", s); e <- as.integer(sub("^.*e", "", s))
    if (!grepl("\\.", m)) m <- paste0(m, ".0")
    return(paste0(m, "e", if (e > 0) "+" else "-", abs(e)))
  }
  if (!grepl("\\.", s)) s <- paste0(s, ".0")
  s
}
.py_mpf_repr <- function(x) .py_mpf_str(x, 17L)
.modint_sym <- function(a) { v <- a$val; if (v > a$mod %/% 2) v <- v - a$mod; as.character(v) }

#' @export
format.ANP <- function(x, ...) .anp_repr(x)
#' @export
print.ANP <- function(x, ...) { cat(format(x), "\n", sep = ""); invisible(x) }

# ---------------------------------------------------------------------------
# ANP: an element of QQ<a> as its coordinates (leading first) modulo the
# minimal polynomial of the primitive element

.anp_cache <- new.env(parent = emptyenv())
.anp_mod <- function(K) {
  d <- .dom(K)
  hit <- .anp_cache[[d$name]]
  if (!is.null(hit)) return(hit)
  ext <- if (length(d$ext) == 1L) d$ext[[1L]] else .add(d$ext)
  X <- Dummy("x")
  mp <- minimal_polynomial(ext, X)
  v <- vapply(.coeff_list(expand(mp), X), function(z) as.character(gmp::as.bigq(z$value)), "")
  assign(d$name, v, envir = .anp_cache)
  v
}
.anp_modq <- function(K) {
  m <- gmp::as.bigq(.anp_mod(K))
  m / m[1L]
}
.anp <- function(rep, K) {
  rep <- gmp::as.bigq(rep)
  while (length(rep) && rep[1L] == 0) rep <- rep[-1L]
  structure(list(rep = rep, K = K), class = "ANP")
}
.anp_gen <- function(K) { d <- .dom(K); if (length(d$ext) == 1L) d$ext[[1L]] else .add(d$ext) }
.anp_to_sym <- function(a) {
  n <- length(a$rep)
  if (!n) return(S_Zero)
  g <- .anp_gen(a$K)
  expand(.add(lapply(seq_len(n), function(i) .mul(list(make_rational(a$rep[i]), .pow(g, .int(n - i)))))))
}
# the coordinates of an expression of the field: a polynomial in the
# primitive element, matched by solving over the power basis numerically
# then exactly (the expression must lie in QQ(a))
.anp_from_sym <- function(K, e) {
  e <- expand(as_sym(e))
  if (e$type %in% c("Integer", "Rational")) return(.anp(e$value, K))
  g <- .anp_gen(K)
  p <- tryCatch(expr_to_qpoly(e, list(g), strict = TRUE), error = function(err) NULL)
  if (is.null(p)) stop(CoercionFailed(paste0("expected an element of ", format(K), ", got ", sstr(e))))
  deg <- max(p$e[, 1L], 0L)
  v <- rep(gmp::as.bigq(0), deg + 1L)
  for (i in seq_along(p$c)) v[deg - p$e[i, 1L] + 1L] <- v[deg - p$e[i, 1L] + 1L] + p$c[i]
  .anp(.anp_reduce(v, .anp_modq(K)), K)
}
.anp_reduce <- function(v, m) {
  v <- gmp::as.bigq(v); n <- length(m) - 1L
  while (length(v) > n) {
    lc <- v[1L]
    if (lc != 0) for (j in seq_along(m)) v[j] <- v[j] - lc * m[j]
    v <- v[-1L]
  }
  while (length(v) && v[1L] == 0) v <- v[-1L]
  v
}
.anp_mulv <- function(a, b) {
  if (!length(a) || !length(b)) return(gmp::as.bigq(integer(0)))
  r <- rep(gmp::as.bigq(0), length(a) + length(b) - 1L)
  for (i in seq_along(a)) for (j in seq_along(b)) r[i + j - 1L] <- r[i + j - 1L] + a[i] * b[j]
  r
}
.anp_addv <- function(a, b) {
  n <- max(length(a), length(b))
  pad <- function(v) c(rep(gmp::as.bigq(0), n - length(v)), v)
  if (!n) return(gmp::as.bigq(integer(0)))
  pad(gmp::as.bigq(a)) + pad(gmp::as.bigq(b))
}
.anp_inv <- function(a) {
  m <- .anp_modq(a$K)
  X <- Dummy("x")
  f <- .add(lapply(seq_along(a$rep), function(i) .mul(list(make_rational(a$rep[i]), .pow(X, .int(length(a$rep) - i))))))
  g <- .add(lapply(seq_along(m), function(i) .mul(list(make_rational(m[i]), .pow(X, .int(length(m) - i))))))
  s <- invert(f, g, X)
  cl <- .coeff_list(expand(s), X)
  .anp(vapply(cl, function(z) as.character(gmp::as.bigq(z$value)), ""), a$K)
}
#' @export
Ops.ANP <- function(e1, e2) {
  K <- if (inherits(e1, "ANP")) e1$K else e2$K
  A <- function(z) if (inherits(z, "ANP")) z else .anp(.as_q(z), K)
  if (missing(e2)) return(if (.Generic == "-") .anp(-e1$rep, K) else e1)
  if (.Generic == "^") {
    k <- as.integer(as.character(e2)); r <- .anp(1, K); b <- if (k < 0) .anp_inv(e1) else e1
    for (i in seq_len(abs(k))) r <- r * b
    return(r)
  }
  a <- A(e1); b <- A(e2); m <- .anp_modq(K)
  switch(.Generic,
    `+` = .anp(.anp_addv(a$rep, b$rep), K), `-` = .anp(.anp_addv(a$rep, -b$rep), K),
    `*` = .anp(.anp_reduce(.anp_mulv(a$rep, b$rep), m), K),
    `/` = { if (!length(b$rep)) stop(.sympy_condition("ZeroDivisionError", "division by zero")); a * .anp_inv(b) },
    `==` = .de_eq(a, b), `!=` = !.de_eq(a, b),
    stop("operator ", .Generic, " is not defined for algebraic numbers", call. = FALSE))
}
#' @export
chooseOpsMethod.ANP <- function(x, y, mx, my, cl, reverse) TRUE
.anp_repr <- function(a) {
  mq <- function(v) paste0("[", paste(vapply(seq_along(v), function(i) paste0("MPQ(", gmp::numerator(v[i]), ",", gmp::denominator(v[i]), ")"), ""), collapse = ", "), "]")
  paste0("ANP(", mq(a$rep), ", ", mq(.anp_modq(a$K)), ", QQ)")
}
