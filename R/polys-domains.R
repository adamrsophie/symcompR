# Polynomial domains, after sympy/polys/domains, sympy/polys/constructor.py,
# sympy/polys/rings.py, sympy/polys/fields.py and sympy/polys/monomials.py.
# A domain is a callable object: ZZ(3), QQ(1, 2), GF(7)(10), ZZ[x](x^2 + 1).
# Elements are native where R has a type (bigz, bigq, double, complex),
# classed objects for modular and Gaussian numbers, Poly objects in
# polynomial rings, and expressions in fraction fields and EX.
# ponytail: elements of algebraic fields are kept as expressions rather than
# as ANP coordinate vectors.

.new_domain <- function(kind, name, ...) {
  d <- list(kind = kind, name = name, ...)
  f <- function(...) .dom_convert(d, ...)
  structure(f, class = c("sym_domain", "function"), dom = d)
}
.dom <- function(K) attr(K, "dom")

#' Polynomial domains
#'
#' SymPy's coefficient domains. `ZZ`, `QQ`, `RR`, `CC`, `EX`, `EXRAW`,
#' `ZZ_I` (Gaussian integers) and `QQ_I` (Gaussian rationals) are domain
#' objects; `GF(p)` (alias `FF`) is the finite field of `p` elements,
#' `poly_ring(K, ...)` (also `K[x, y]` and `PolynomialRing()`) a polynomial
#' ring, `frac_field(K, ...)` (`FractionField()`) a field of rational
#' functions, and `algebraic_field(K, ...)` (`AlgebraicField()`) an algebraic
#' number field such as `QQ<sqrt(2)>`; `CyclotomicField(n)` is `QQ<zeta_n>`.
#' Calling a domain converts its arguments to an element: `ZZ(3)` is a
#' `bigz`, `QQ(1, 2)` a `bigq`, `GF(7)(10)` the residue `3 mod 7`, and
#' `ZZ[x](x^2 + 1)` a [Poly()]. `to_sympy()` and `from_sympy()` convert
#' between elements and expressions (with one argument, `to_sympy(m)` is the
#' quantum matrix helper of [matrix_tensor_product()]), and `Domain(K)`
#' tests for a domain.
#' The names `IntegerRing()`, `RationalField()`, `RealField()`,
#' `ComplexField()`, `ExpressionDomain()`, `FiniteField()` and the Python
#' and GMPY variants of SymPy (`ZZ_python`, `QQ_gmpy`, `FF_python`, ...)
#' return the same domains.
#'
#' @param p A prime modulus, or the numerator for `PythonRational()`.
#' @param q The denominator for `PythonRational()`.
#' @param K A domain.
#' @param ... Generators, extensions or values.
#' @param prec Precision in bits (kept for SymPy compatibility).
#' @param n The order of the root of unity.
#' @param a An element or expression.
#' @param x A domain.
#' @return A domain, an element, or an expression.
#' @examples
#' ZZ(3); QQ(1, 2); GF(7)(10)
#' x <- Symbol("x")
#' ZZ[x]
#' frac_field(QQ, x)
#' algebraic_field(QQ, sqrt(Integer(2)))
#' GF(7)(3) * GF(7)(5)
#' @name domains
NULL

#' @rdname domains
#' @export
ZZ <- .new_domain("ZZ", "ZZ")
#' @rdname domains
#' @export
QQ <- .new_domain("QQ", "QQ")
#' @rdname domains
#' @export
RR <- .new_domain("RR", "RR")
#' @rdname domains
#' @export
CC <- .new_domain("CC", "CC")
#' @rdname domains
#' @export
EX <- .new_domain("EX", "EX")
#' @rdname domains
#' @export
EXRAW <- .new_domain("EX", "EXRAW")
#' @rdname domains
#' @export
ZZ_I <- .new_domain("ZZ_I", "ZZ_I")
#' @rdname domains
#' @export
QQ_I <- .new_domain("QQ_I", "QQ_I")
#' @rdname domains
#' @export
GF <- function(p) {
  p <- gmp::as.bigz(if (is_sym(p)) .Z(p) else p)
  if (!isTRUE(isprime(p))) stop(DomainError(paste0("modulus must be prime, got ", p)))
  .new_domain("GF", paste0("GF(", p, ")"), mod = p)
}
#' @rdname domains
#' @export
FF <- GF
#' @rdname domains
#' @export
poly_ring <- function(K, ...) {
  gens <- lapply(list(...), as_sym)
  .new_domain("PolyRing", paste0(format(K), "[", paste(vapply(gens, .str, ""), collapse = ","), "]"), dom = K, gens = gens)
}
#' @rdname domains
#' @export
frac_field <- function(K, ...) {
  gens <- lapply(list(...), as_sym)
  .new_domain("FracField", paste0(format(K), "(", paste(vapply(gens, .str, ""), collapse = ","), ")"), dom = K, gens = gens)
}
#' @rdname domains
#' @export
algebraic_field <- function(K, ...) {
  ext <- lapply(list(...), as_sym)
  .new_domain("Algebraic", paste0("QQ<", .str(.add(ext)), ">"), dom = K, ext = ext)
}
#' @rdname domains
#' @export
CyclotomicField <- function(n) {
  z <- fcall("exp", list(.div(.mul(list(S_Two, S_Pi, S_I)), as_sym(n))))
  .new_domain("Algebraic", paste0("QQ<zeta", n, ">"), dom = QQ, ext = list(z))
}

#' @export
`[.sym_domain` <- function(x, ...) poly_ring(x, ...)
#' @export
format.sym_domain <- function(x, ...) .dom(x)$name
#' @export
print.sym_domain <- function(x, ...) { cat(format(x), "\n", sep = ""); invisible(x) }
#' @export
Ops.sym_domain <- function(e1, e2) {
  eq <- inherits(e1, "sym_domain") && inherits(e2, "sym_domain") && identical(format(e1), format(e2))
  switch(.Generic, `==` = eq, `!=` = !eq, stop("operator ", .Generic, " is not defined for domains", call. = FALSE))
}

#' @rdname domains
#' @export
Domain <- function(x) inherits(x, "sym_domain")
#' @rdname domains
#' @export
IntegerRing <- function() ZZ
#' @rdname domains
#' @export
RationalField <- function() QQ
#' @rdname domains
#' @export
RealField <- function(prec = 53) RR
#' @rdname domains
#' @export
ComplexField <- function(prec = 53) CC
#' @rdname domains
#' @export
ExpressionDomain <- function() EX
#' @rdname domains
#' @export
FiniteField <- function(p) GF(p)
#' @rdname domains
#' @export
PolynomialRing <- function(K, ...) poly_ring(K, ...)
#' @rdname domains
#' @export
FractionField <- function(K, ...) frac_field(K, ...)
#' @rdname domains
#' @export
AlgebraicField <- function(K, ...) algebraic_field(K, ...)
#' @rdname domains
#' @export
PythonIntegerRing <- IntegerRing
#' @rdname domains
#' @export
GMPYIntegerRing <- IntegerRing
#' @rdname domains
#' @export
ZZ_python <- IntegerRing
#' @rdname domains
#' @export
ZZ_gmpy <- IntegerRing
#' @rdname domains
#' @export
GMPYRationalField <- RationalField
#' @rdname domains
#' @export
QQ_python <- RationalField
#' @rdname domains
#' @export
QQ_gmpy <- RationalField
#' @rdname domains
#' @export
PythonFiniteField <- FiniteField
#' @rdname domains
#' @export
GMPYFiniteField <- FiniteField
#' @rdname domains
#' @export
FF_python <- FiniteField
#' @rdname domains
#' @export
FF_gmpy <- FiniteField

#' @rdname domains
#' @export
PythonRational <- function(p, q = 1) gmp::as.bigq(p, q)

.as_q <- function(a) {
  if (is_sym(a)) {
    if (a$type %in% c("Integer", "Rational")) return(gmp::as.bigq(a$value))
    stop(CoercionFailed(paste0("cannot convert ", .str(a), " to a rational number")))
  }
  if (inherits(a, "bigz") || inherits(a, "bigq")) return(gmp::as.bigq(a))
  if (is.numeric(a) && a == round(a)) return(gmp::as.bigq(gmp::as.bigz(a)))
  if (is.numeric(a)) return(gmp::as.bigq(a))
  stop(CoercionFailed(paste0("cannot convert ", format(a), " to a rational number")))
}

.dom_convert <- function(d, ...) {
  a <- list(...)
  switch(EXPR = d$kind,
    ZZ = {
      q <- .as_q(a[[1L]])
      if (gmp::denominator(q) != 1) stop(CoercionFailed(paste0("cannot convert ", q, " to ZZ")))
      gmp::as.bigz(gmp::numerator(q))
    },
    QQ = if (length(a) == 2L) .as_q(a[[1L]]) / .as_q(a[[2L]]) else .as_q(a[[1L]]),
    RR = as.numeric(if (is_sym(a[[1L]])) num_value_dbl(N(a[[1L]])) else a[[1L]]),
    CC = if (length(a) == 2L) complex(real = as.numeric(a[[1L]]), imaginary = as.numeric(a[[2L]]))
         else as.complex(if (is_sym(a[[1L]])) .neval(a[[1L]], list()) else a[[1L]]),
    EX = as_sym(a[[1L]]),
    ZZ_I = , QQ_I = {
      if (length(a) == 2L) { re <- .as_q(a[[1L]]); im <- .as_q(a[[2L]]) }
      else { ri <- .re_im(as_sym(a[[1L]])); re <- .as_q(ri[[1L]]); im <- .as_q(ri[[2L]]) }
      if (d$kind == "ZZ_I" && (gmp::denominator(re) != 1 || gmp::denominator(im) != 1))
        stop(CoercionFailed("Gaussian integers need integer parts"))
      structure(list(re = re, im = im, kind = d$kind), class = "sym_gaussian")
    },
    GF = .modint(.as_q(a[[1L]]), d$mod),
    PolyRing = do.call(Poly, c(list(as_sym(a[[1L]])), d$gens)),
    FracField = cancel(as_sym(a[[1L]])),
    Algebraic = expand(radsimp(as_sym(a[[1L]]))))
}

.modint <- function(q, p) {
  v <- gmp::numerator(q) * gmp::inv.bigz(gmp::denominator(q), p)
  structure(list(val = gmp::mod.bigz(v, p), mod = p), class = "sym_modint")
}
#' @export
format.sym_modint <- function(x, ...) paste0(as.character(x$val), " mod ", as.character(x$mod))
#' @export
print.sym_modint <- function(x, ...) { cat(format(x), "\n", sep = ""); invisible(x) }
#' @export
Ops.sym_modint <- function(e1, e2) {
  p <- if (inherits(e1, "sym_modint")) e1$mod else e2$mod
  v <- function(z) if (inherits(z, "sym_modint")) z$val else gmp::mod.bigz(gmp::as.bigz(if (is_sym(z)) .Z(z) else z), p)
  if (missing(e2)) return(if (.Generic == "-") .modint(gmp::as.bigq(-e1$val), p) else e1)
  if (.Generic == "^") {
    k <- gmp::as.bigz(if (is_sym(e2)) .Z(e2) else e2)
    b <- if (k < 0) gmp::inv.bigz(v(e1), p) else v(e1)
    return(.modint(gmp::as.bigq(gmp::powm(b, abs(k), p)), p))
  }
  a <- v(e1); b <- v(e2)
  switch(.Generic,
    `+` = .modint(gmp::as.bigq(a + b), p), `-` = .modint(gmp::as.bigq(a - b), p),
    `*` = .modint(gmp::as.bigq(a * b), p), `/` = .modint(gmp::as.bigq(a * gmp::inv.bigz(b, p)), p),
    `==` = a == b, `!=` = a != b,
    stop("operator ", .Generic, " is not defined for GF elements", call. = FALSE))
}
#' @export
chooseOpsMethod.sym_modint <- function(x, y, mx, my, cl, reverse) TRUE

.gauss_sym <- function(z) .add(list(make_rational(z$re), .mul(list(make_rational(z$im), S_I))))
#' @export
format.sym_gaussian <- function(x, ...) .str(.gauss_sym(x))
#' @export
print.sym_gaussian <- function(x, ...) { cat(format(x), "\n", sep = ""); invisible(x) }
#' @export
Ops.sym_gaussian <- function(e1, e2) {
  kind <- if (inherits(e1, "sym_gaussian")) e1$kind else e2$kind
  g <- function(z) if (inherits(z, "sym_gaussian")) z else .dom_convert(list(kind = "QQ_I"), z)
  mk <- function(re, im) structure(list(re = re, im = im, kind = kind), class = "sym_gaussian")
  if (missing(e2)) return(if (.Generic == "-") mk(-e1$re, -e1$im) else e1)
  a <- g(e1)
  if (.Generic == "^") {
    k <- as.integer(if (is_sym(e2)) .Z(e2) else e2)
    r <- mk(gmp::as.bigq(1), gmp::as.bigq(0))
    for (i in seq_len(abs(k))) r <- r * a
    return(if (k < 0) mk(gmp::as.bigq(1), gmp::as.bigq(0)) / r else r)
  }
  b <- g(e2)
  switch(.Generic,
    `+` = mk(a$re + b$re, a$im + b$im), `-` = mk(a$re - b$re, a$im - b$im),
    `*` = mk(a$re * b$re - a$im * b$im, a$re * b$im + a$im * b$re),
    `/` = { n <- b$re^2 + b$im^2; mk((a$re * b$re + a$im * b$im) / n, (a$im * b$re - a$re * b$im) / n) },
    `==` = a$re == b$re && a$im == b$im, `!=` = !(a$re == b$re && a$im == b$im),
    stop("operator ", .Generic, " is not defined for Gaussian numbers", call. = FALSE))
}
#' @export
chooseOpsMethod.sym_gaussian <- function(x, y, mx, my, cl, reverse) TRUE

#' @rdname domains
#' @export
to_sympy <- function(K, a) {
  # one argument: sympy.physics.quantum.matrixutils.to_sympy(m)
  if (missing(a)) return(.qto_sympy(K))
  if (inherits(a, "sym_modint")) {
    # SymPy's symmetric representation of residues
    v <- a$val; if (v > a$mod %/% 2) v <- v - a$mod
    return(make_rational(gmp::as.bigq(v)))
  }
  if (inherits(a, "sym_gaussian")) return(.gauss_sym(a))
  if (inherits(a, "Poly")) return(as_expr(a))
  if (inherits(a, "bigz") || inherits(a, "bigq")) return(make_rational(gmp::as.bigq(a)))
  as_sym(a)
}
#' @rdname domains
#' @export
from_sympy <- function(K, a) K(as_sym(a))

#' Constructing the minimal domain
#'
#' `construct_domain(obj)` finds the smallest domain containing the
#' coefficients `obj` (an expression or a list of them) and converts them,
#' as SymPy's `construct_domain()`: `ZZ` or `QQ` for rational numbers, `RR`
#' or `CC` for floats, `ZZ_I` or `QQ_I` for Gaussian rationals, a polynomial
#' ring or a fraction field over those for expressions in symbols, and `EX`
#' otherwise; with `extension = TRUE` a single radical gives an algebraic
#' field.
#'
#' @param obj An expression or a list of expressions.
#' @param extension Use algebraic fields for radicals.
#' @param field Return a field (`QQ` rather than `ZZ`, a fraction field
#'   rather than a polynomial ring).
#' @return A list of the domain and the converted elements.
#' @examples
#' construct_domain(list(1, Rational(1, 2)))
#' x <- Symbol("x")
#' construct_domain(list(x, 1))[[1L]]
#' construct_domain(list(1/x))[[1L]]
#' @export
construct_domain <- function(obj, extension = FALSE, field = FALSE) {
  single <- !is.list(obj) || is_sym(obj)
  es <- lapply(if (single) list(obj) else obj, as_sym)
  conv <- function(K) { out <- lapply(es, K); if (single) out[[1L]] else out }
  K <- .cd_simple(es, extension)
  if (is.null(K)) K <- .cd_composite(es, field)
  if (is.null(K) || isFALSE(K)) K <- EX
  list(K, conv(K))
}

# pure_complex: list(re, im) for a + b*I with numeric a, b and b != 0
.cd_pure_complex <- function(z) {
  if (!.contains_key(z, S_I$key) || length(free_symbols(z))) return(NULL)
  ri <- .re_im(z)
  ok <- function(v) v$type %in% c("Integer", "Rational", "Float")
  if (ok(ri[[1L]]) && ok(ri[[2L]]) && !is_zero_num(ri[[2L]])) ri else NULL
}

# _construct_simple: a domain, FALSE for EX, or NULL to try composite
.cd_simple <- function(es, extension) {
  rationals <- FALSE; floats <- FALSE; complexes <- FALSE; algebraics <- FALSE
  for (z in es) {
    if (z$type %in% c("Integer", "Rational")) { if (z$type == "Rational") rationals <- TRUE }
    else if (z$type == "Float") floats <- TRUE
    else {
      pc <- .cd_pure_complex(z)
      if (!is.null(pc)) {
        complexes <- TRUE
        if (pc[[1L]]$type == "Float" || pc[[2L]]$type == "Float") floats <- TRUE
        else if (pc[[1L]]$type == "Rational" || pc[[2L]]$type == "Rational") rationals <- TRUE
      } else if (isTRUE(extension) && !length(free_symbols(z)) && isTRUE(ask_fact(z, "algebraic"))) algebraics <- TRUE
      else return(NULL)
    }
  }
  if (algebraics) {
    if (floats) return(FALSE)
    rad <- unique(lapply(Filter(function(p) p$type == "Pow" && p$args[[2L]]$type == "Rational",
                                do.call(c, lapply(es, function(z) atoms(z, "Pow")))), function(p) p))
    if (length(rad) == 1L) return(algebraic_field(QQ, rad[[1L]]))
    return(FALSE)
  }
  if (floats && complexes) return(CC)
  if (floats) return(RR)
  if (complexes) return(if (rationals) QQ_I else ZZ_I)
  if (rationals) return(QQ)
  ZZ
}

# _construct_composite
.cd_composite <- function(es, field = FALSE) {
  nds <- lapply(es, .as_numer_denom)
  parts <- lapply(c(lapply(nds, `[[`, 1L), lapply(nds, `[[`, 2L)), expand)
  gens <- Filter(function(g) !identical(g$key, S_I$key), .find_gens(parts))
  if (!length(gens)) return(NULL)
  if (any(vapply(gens, function(g) !length(free_symbols(g)) && isTRUE(ask_fact(g, "algebraic")), NA))) return(NULL)
  seen <- character()
  for (g in gens) {
    sy <- vapply(free_symbols(g), `[[`, "", "key")
    if (length(intersect(seen, sy))) return(NULL)
    seen <- c(seen, sy)
  }
  hasI <- any(vapply(parts, function(z) .contains_key(z, S_I$key), NA))
  allg <- if (hasI) c(gens, list(S_I)) else gens
  qps <- lapply(parts, function(z) tryCatch(expr_to_qpoly(z, allg, strict = FALSE), error = function(e) NULL))
  if (any(vapply(qps, is.null, NA))) return(NULL)
  k <- length(es)
  n <- length(gens)
  dens <- qps[(k + 1L):(2L * k)]
  isconst <- function(p) qp_is_zero(p) || (length(p$c) == 1L && all(p$e[1L, seq_len(n)] == 0L))
  fractions <- isTRUE(field) || any(!vapply(dens, isconst, NA))
  complexes <- hasI && any(vapply(qps, function(p) !qp_is_zero(p) && any(p$e[, n + 1L] > 0L), NA))
  rationals <- FALSE
  if (!fractions) {
    for (i in seq_len(k)) {
      dv <- if (qp_is_zero(dens[[i]])) gmp::as.bigq(1) else dens[[i]]$c[1L]
      if (any(gmp::denominator(qps[[i]]$c / dv) != 1)) rationals <- TRUE
    }
  } else if (any(vapply(qps, function(p) any(gmp::denominator(p$c) != 1), NA))) rationals <- TRUE
  floats <- any(vapply(es, function(z) .has_type(z, "Float"), NA))
  ground <- if (floats && complexes) CC else if (floats) RR else if (complexes) (if (rationals) QQ_I else ZZ_I) else if (rationals) QQ else ZZ
  if (!fractions) do.call(poly_ring, c(list(ground), gens))
  else do.call(frac_field, c(list(ground), gens))
}

#' Sparse polynomial rings and rational function fields
#'
#' `ring(symbols, domain)` returns a polynomial ring over `domain` followed
#' by its generators, as SymPy's `ring()`; `xring()` returns the ring and a
#' list of the generators, `vring()` assigns the generators in `envir` and
#' returns the ring, and `sring(exprs)` builds the ring of a list of
#' expressions and returns it with their elements. `field()`, `xfield()`,
#' `vfield()` and `sfield()` do the same for fields of rational functions.
#' Ring elements are polynomials kept expanded; field elements are rational
#' functions kept in lowest terms.
#'
#' @param symbols Names as a comma separated string, or a list of symbols.
#' @param domain The coefficient domain.
#' @param order The monomial order (only `"lex"`).
#' @param exprs An expression or a list of them.
#' @param envir The environment receiving the generators.
#' @param e1,e2,x Ring or field elements.
#' @param y,mx,my,cl,reverse Arguments of the `chooseOpsMethod()` generic.
#' @param ... Unused.
#' @return A list with the ring or field and the generators or elements.
#' @examples
#' r <- ring("a,b", ZZ); R <- r[[1L]]; a <- r[[2L]]; b <- r[[3L]]
#' R
#' (a + b)^2
#' f <- field("u,v", QQ); u <- f[[2L]]; v <- f[[3L]]
#' (u^2 - v^2)/(u - v)
#' @name rings
NULL

.ring_syms <- function(symbols) {
  if (is.character(symbols)) return(lapply(trimws(strsplit(symbols, ",", fixed = TRUE)[[1L]]), Symbol))
  lapply(if (is_sym(symbols)) list(symbols) else symbols, as_sym)
}
.new_ring <- function(kind, gens, domain) {
  what <- if (kind == "ring") "Polynomial ring" else "Rational function field"
  structure(list(kind = kind, gens = gens, domain = domain,
                 name = paste0(what, " in ", paste(vapply(gens, .str, ""), collapse = ", "), " over ", format(domain), " with lex order")),
            class = "sym_ring")
}
.relem <- function(R, e) {
  e <- if (R$kind == "ring") expand(e) else cancel(together(e))
  structure(list(expr = e, ring = R), class = "sym_ring_element")
}
#' @rdname rings
#' @export
format.sym_ring <- function(x, ...) x$name
#' @rdname rings
#' @export
print.sym_ring <- function(x, ...) { cat(format(x), "\n", sep = ""); invisible(x) }
#' @rdname rings
#' @export
format.sym_ring_element <- function(x, ...) .relem_str(x)
#' @rdname rings
#' @export
print.sym_ring_element <- function(x, ...) { cat(format(x), "\n", sep = ""); invisible(x) }
#' @rdname rings
#' @export
Ops.sym_ring_element <- function(e1, e2) {
  R <- if (inherits(e1, "sym_ring_element")) e1$ring else e2$ring
  ex <- function(z) if (inherits(z, "sym_ring_element")) z$expr else as_sym(z)
  if (missing(e2)) return(if (.Generic == "-") .relem(R, .neg(e1$expr)) else e1)
  a <- ex(e1); b <- ex(e2)
  if (.Generic %in% c("==", "!=")) { eq <- is_zero_num(expand(.sub(a, b))); return(if (.Generic == "==") eq else !eq) }
  r <- switch(.Generic, `+` = .add(list(a, b)), `-` = .sub(a, b), `*` = .mul(list(a, b)),
              `^` = {
                if (is_number_type(b) && isTRUE(num_sign(b) < 0) && is_zero_num(expand(a))) stop(.sympy_condition("ZeroDivisionError", ""))
                .pow(a, b)
              },
              `/` = {
                q <- .relem_div(R, a, b)
                if (!is.null(q)) return(.relem(R, q))
                q <- .div(a, b)
                if (R$kind == "ring") {
                  qc <- cancel(q)
                  if (!is_one_num(.as_numer_denom(qc)[[2L]])) stop(ExactQuotientFailed(paste0(.str(b), " does not divide ", .str(a))))
                  qc
                } else q
              },
              stop("operator ", .Generic, " is not defined for ring elements", call. = FALSE))
  .relem(R, r)
}
#' @rdname rings
#' @export
chooseOpsMethod.sym_ring_element <- function(x, y, mx, my, cl, reverse) TRUE
#' @rdname rings
#' @export
ring <- function(symbols, domain = ZZ, order = "lex") {
  gens <- .ring_syms(symbols); R <- .new_ring("ring", gens, domain)
  c(list(R), lapply(gens, function(g) .relem(R, g)))
}
#' @rdname rings
#' @export
xring <- function(symbols, domain = ZZ, order = "lex") { r <- ring(symbols, domain, order); list(r[[1L]], r[-1L]) }
#' @rdname rings
#' @export
vring <- function(symbols, domain = ZZ, order = "lex", envir = parent.frame()) {
  r <- ring(symbols, domain, order)
  for (g in r[-1L]) assign(g$expr$name, g, envir = envir)
  r[[1L]]
}
#' @rdname rings
#' @export
sring <- function(exprs, ...) {
  cd <- construct_domain(exprs)
  K <- cd[[1L]]; d <- .dom(K)
  gens <- if (d$kind %in% c("PolyRing", "FracField")) d$gens else .poly_sort_gens(unique_syms(lapply(if (is.list(exprs)) exprs else list(exprs), function(z) free_symbols(as_sym(z)))))
  base <- if (d$kind %in% c("PolyRing", "FracField")) d$dom else K
  R <- .new_ring(if (d$kind == "FracField") "field" else "ring", gens, base)
  els <- lapply(if (is.list(exprs) && !is_sym(exprs)) exprs else list(exprs), function(z) .relem(R, as_sym(z)))
  list(R, if (is.list(exprs) && !is_sym(exprs)) els else els[[1L]])
}
#' @rdname rings
#' @export
field <- function(symbols, domain = ZZ, order = "lex") {
  gens <- .ring_syms(symbols); R <- .new_ring("field", gens, domain)
  c(list(R), lapply(gens, function(g) .relem(R, g)))
}
#' @rdname rings
#' @export
xfield <- function(symbols, domain = ZZ, order = "lex") { r <- field(symbols, domain, order); list(r[[1L]], r[-1L]) }
#' @rdname rings
#' @export
vfield <- function(symbols, domain = ZZ, order = "lex", envir = parent.frame()) {
  r <- field(symbols, domain, order)
  for (g in r[-1L]) assign(g$expr$name, g, envir = envir)
  r[[1L]]
}
#' @rdname rings
#' @export
sfield <- function(exprs, ...) {
  s <- sring(exprs)
  R <- s[[1L]]; R$kind <- "field"
  R$name <- sub("^Polynomial ring", "Rational function field", R$name)
  els <- if (inherits(s[[2L]], "sym_ring_element")) .relem(R, s[[2L]]$expr) else lapply(s[[2L]], function(z) .relem(R, z$expr))
  list(R, els)
}

#' Monomials
#'
#' `Monomial(exponents, gens)` is a monomial given by its exponent vector, as
#' SymPy's `Monomial`. `*`, `/` and `^` combine monomials, and
#' `as_expr(m, ...)` writes it in the generators.
#'
#' @param exponents Nonnegative integer exponents.
#' @param gens Optional generators.
#' @param x,e1,e2 Monomials.
#' @param ... Generators for `as_expr()`.
#' @return A monomial, or an expression from `as_expr()`.
#' @examples
#' m <- Monomial(c(1, 2))
#' m
#' m * Monomial(c(2, 1))
#' @export
Monomial <- function(exponents, gens = NULL) {
  ex <- as.integer(unlist(exponents))
  if (any(ex < 0L)) stop("a monomial has nonnegative exponents", call. = FALSE)
  structure(list(exponents = ex, gens = gens), class = "sym_monomial")
}
#' @rdname Monomial
#' @export
format.sym_monomial <- function(x, ...) paste0("Monomial((", paste(x$exponents, collapse = ", "), if (length(x$exponents) == 1L) ",))" else "))")
#' @rdname Monomial
#' @export
print.sym_monomial <- function(x, ...) { cat(format(x), "\n", sep = ""); invisible(x) }
#' @rdname Monomial
#' @export
Ops.sym_monomial <- function(e1, e2) {
  switch(.Generic,
    `*` = Monomial(e1$exponents + e2$exponents, e1$gens),
    `/` = Monomial(e1$exponents - e2$exponents, e1$gens),
    `^` = Monomial(e1$exponents * as.integer(e2), e1$gens),
    `==` = identical(e1$exponents, e2$exponents), `!=` = !identical(e1$exponents, e2$exponents),
    stop("operator ", .Generic, " is not defined for monomials", call. = FALSE))
}
.monomial_expr <- function(m, gens) {
  if (!length(gens)) gens <- m$gens
  if (length(gens) != length(m$exponents)) stop("the number of generators must match the exponents", call. = FALSE)
  .mul(lapply(seq_along(gens), function(i) .pow(as_sym(gens[[i]]), .int(m$exponents[i]))))
}

#' Polynomial options
#'
#' `Options(gens, args)` processes the options of polynomial functions, as
#' SymPy's `Options`: generators, `domain` (a domain or its name), `field`,
#' `greedy`, `expand`, `extension`, `modulus` and the others, with SymPy's
#' defaults.
#'
#' @param gens Generators.
#' @param args A named list of options.
#' @return A list of options.
#' @examples
#' Options(list(Symbol("x")), list(domain = "ZZ"))$domain
#' @export
Options <- function(gens = list(), args = list()) {
  allowed <- c("domain", "field", "greedy", "expand", "extension", "modulus", "symmetric", "strict", "auto",
               "frac", "formal", "polys", "include", "all", "gaussian", "order", "split", "sort", "wrt", "series",
               "gen", "symbols", "method")
  bad <- setdiff(names(args), allowed)
  if (length(bad)) stop(OptionError(paste0("'", bad[1L], "' is not a valid option")))
  opt <- list(gens = lapply(gens, as_sym), expand = TRUE, greedy = TRUE, auto = TRUE, strict = FALSE, frac = FALSE,
              formal = FALSE, polys = NULL, include = FALSE, all = FALSE, order = "lex")
  for (k in names(args)) opt[[k]] <- args[[k]]
  if (is.character(opt$domain)) {
    d <- opt$domain
    opt$domain <- switch(EXPR = d, ZZ = ZZ, QQ = QQ, RR = RR, CC = CC, EX = EX, EXRAW = EXRAW, ZZ_I = ZZ_I, QQ_I = QQ_I,
                         if (grepl("^(GF|FF)\\(([0-9]+)\\)$", d)) GF(as.numeric(sub("^(GF|FF)\\(([0-9]+)\\)$", "\\2", d)))
                         else stop(OptionError(paste0("expected a valid domain specification, got ", d))))
  }
  if (!is.null(opt$domain) && !is.null(opt$field) && isTRUE(opt$field))
    stop(OptionError("'field' option is not allowed together with 'domain'"))
  if (!is.null(opt$modulus) && !is.null(opt$domain)) stop(OptionError("'modulus' option is not allowed together with 'domain'"))
  if (!is.null(opt$modulus)) opt$domain <- GF(opt$modulus)
  structure(opt, class = "sym_poly_options")
}
