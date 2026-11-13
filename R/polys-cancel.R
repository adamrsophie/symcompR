# cancel() of a pair (p, q), Poly.cancel() and the str() of the elements of
# ring() and field(), ports of sympy/polys/polytools.py (cancel,
# Poly.cancel), sympy/polys/densetools.py (dmp_cancel) and the StrPrinter
# methods for PolyElement and FracElement.

# cancel((p, q), *gens): sring((p, q)) then PolyElement.cancel; returns the
# tuple (1, P, Q), with Polys when both inputs are Polys or polys = TRUE
.cancel_tuple <- function(f, gens, polys = NULL, opts = list()) {
  if (length(f) != 2L)
    stop(.sympy_condition("ValueError", paste0("unexpected argument: (", paste(vapply(f, function(z) .str(as_sym(z)), ""), collapse = ", "), ")")))
  p <- f[[1L]]; q <- f[[2L]]
  r <- .cancel_opts(p, q, opts)
  if (!is.null(r)) return(.pytuple(c(list(S_One), r)))
  pgens <- NULL; dom <- NULL
  if (inherits(p, "Poly") && inherits(q, "Poly")) {
    pgens <- p$gens[seq_len(.poly_ngens(p))]
    dom <- .poly_unified_dom(p, q)
    if (is.null(polys)) polys <- TRUE
  }
  p <- if (inherits(p, "Poly")) as_expr(p) else signsimp(as_sym(p))
  q <- if (inherits(q, "Poly")) as_expr(q) else signsimp(as_sym(q))
  g <- if (length(gens)) gens else .find_gens(list(expand(p), expand(q)))
  if (!length(g)) return(.pytuple(list(S_One, p, q)))
  cv <- tryCatch(.to_qpolys(list(expand(p), expand(q)), g), error = function(e) NULL)
  if (is.null(cv)) return(.pytuple(list(S_One, p, q)))
  pq <- .qp_cancel(cv$polys[[1L]], cv$polys[[2L]])
  P <- qpoly_to_expr(pq[[1L]], cv$gens); Q <- qpoly_to_expr(pq[[2L]], cv$gens)
  if (!isTRUE(polys)) return(.pytuple(list(S_One, P, Q)))
  pg <- if (!is.null(pgens)) pgens else if (length(gens)) gens else cv$gens
  mk <- function(e) { R <- do.call(Poly, c(list(e), pg)); if (!is.null(dom)) R$dom_str <- dom; R }
  .pytuple(list(S_One, mk(P), mk(Q)))
}

# the options of cancel() that change the domain: greedy = FALSE and
# extension = TRUE (the numbers are coefficients) and modulus; NULL when
# they are absent or p and q are not univariate
.cancel_opts <- function(p, q, opts) {
  if (inherits(p, "Poly")) p <- as_expr(p)
  if (inherits(q, "Poly")) q <- as_expr(q)
  if (!is.null(opts$modulus)) return(.cancel_mod(p, q, opts$modulus))
  if (isFALSE(opts$greedy) || isTRUE(opts$extension)) return(.cancel_ex(p, q, isTRUE(opts$extension)))
  NULL
}
# cancel(f, ...) for an expression with those options
.cancel_expr_opts <- function(e, opts) {
  if (is.null(opts$modulus) && !isFALSE(opts$greedy) && !isTRUE(opts$extension)) return(NULL)
  nd <- .as_numer_denom(factor_terms(signsimp(e), radical = TRUE))
  r <- .cancel_opts(nd[[1L]], nd[[2L]], opts)
  if (is.null(r)) NULL else .div(r[[1L]], r[[2L]])
}

# dmp_cancel over GF(p) (a field with no associated ring), coefficients in
# the symmetric representation
.cancel_mod <- function(p, q, m) {
  p <- expand(as_sym(p)); q <- expand(as_sym(q))
  g <- .find_gens(list(p, q))
  if (length(g) != 1L) return(NULL)
  x <- g[[1L]]
  gfv <- function(e) .gf_trim(vapply(.coeff_list(e, x), function(c) {
    v <- .as_q(c)
    (as.numeric(gmp::mod.bigz(gmp::numerator(v), m)) * .gf_inv(as.numeric(gmp::mod.bigz(gmp::denominator(v), m)), m)) %% m
  }, 0))
  a <- gfv(p); b <- gfv(q)
  if (!length(b)) return(NULL)
  h <- .gf_gcd(a, b, m)
  list(.gf_to_expr(.gf_divmod(a, h, m)$q, x, m), .gf_to_expr(.gf_divmod(b, h, m)$q, x, m))
}

# Fraction field domains of Polys. A Poly over ZZ(a) keeps its coefficients
# in the usual representation (a reciprocal such as 1/(a + 1) is an extra
# generator) and records the domain in dom_str; Poly.cancel works over the
# fraction field itself.
.poly_extra_gens <- function(P) { ng <- .poly_ngens(P); if (ng < length(P$gens)) P$gens[-seq_len(ng)] else list() }
.is_reciprocal <- function(g) g$type == "Pow" && is_number_type(g$args[[2L]]) && isTRUE(num_sign(g$args[[2L]]) < 0)
.syms_of <- function(es) {
  s <- unlist(lapply(es, free_symbols), recursive = FALSE)
  s <- s[!duplicated(vapply(s, .k, ""))]
  s[order(vapply(s, .str, ""))]
}
# a Poly whose extra generators include a reciprocal has rational function
# coefficients: its domain is ZZ(symbols), as construct_domain finds
.poly_frac_tag <- function(P) {
  ex <- .poly_extra_gens(P)
  if (!length(ex) || !any(vapply(ex, .is_reciprocal, NA))) return(P)
  P$dom_str <- paste0("ZZ(", paste(vapply(.syms_of(ex), .str, ""), collapse = ","), ")")
  P
}
# Poly(expr, ..., domain = "ZZ[a]"), "ZZ(a)" or QQ: the Poly labelled with
# that domain; NULL for the other options
.poly_with_domain <- function(expr, pa) {
  d <- pa$opts$domain
  if (length(pa$opts) != 1L || is.null(d)) return(NULL)
  if (identical(if (is.character(d)) d else tryCatch(format(d), error = function(e) ""), "QQ")) {
    # integer coefficients over QQ keep the field
    P <- do.call(Poly, c(list(expr), pa$gens))
    if (.poly_ngens(P) == length(P$gens) && is.null(P$dom_str)) P$dom_str <- "QQ"
    return(P)
  }
  # domain = EX: the coefficients stay expressions; the label is kept
  if (identical(d, "EX")) {
    P <- do.call(Poly, c(list(expr), pa$gens))
    P$dom_str <- "EX"
    return(P)
  }
  if (!is.character(d) || !grepl("^(ZZ|QQ)[[(].*[])]$", d)) return(NULL)
  P <- do.call(Poly, c(list(expr), pa$gens))
  P$dom_str <- gsub(" ", "", d, fixed = TRUE)
  P
}
# the unified fraction field of f and g, list(label, K), or NULL when
# neither is over a fraction field
.poly_frac_dom <- function(f, g) {
  parts <- function(P) {
    d <- P$dom_str
    if (is.null(d)) return(NULL)
    m <- regmatches(d, regexec("^(ZZ|QQ)([[(])(.*)[])]$", d))[[1L]]
    if (!length(m)) return(NULL)
    list(ground = m[2L], frac = m[3L] == "(", syms = lapply(strsplit(m[4L], ",", fixed = TRUE)[[1L]], sympify))
  }
  pf <- parts(f); pg <- parts(g)
  if (!isTRUE(pf$frac) && !isTRUE(pg$frac)) return(NULL)
  syms <- .syms_of(c(pf$syms, pg$syms, .poly_extra_gens(f), .poly_extra_gens(g)))
  gr <- if ("QQ" %in% c(pf$ground, pg$ground)) QQ else ZZ
  list(label = paste0(format(gr), "(", paste(vapply(syms, .str, ""), collapse = ","), ")"),
       K = do.call(frac_field, c(list(gr), syms)))
}

# the domain of the result of Poly.cancel: the unified domain of f and g
.poly_unified_dom <- function(f, g) {
  ex <- function(P) { ng <- .poly_ngens(P); if (ng < length(P$gens)) P$gens[-seq_len(ng)] else list() }
  gr <- if ("QQ" %in% c(.domain_of(f$rep), .domain_of(g$rep), f$dom_str, g$dom_str)) "QQ" else "ZZ"
  e <- c(ex(f), ex(g)); e <- e[!duplicated(vapply(e, .k, ""))]
  if (length(e)) return(paste0(gr, "[", paste(vapply(e, .str, ""), collapse = ","), "]"))
  if (gr == "QQ") "QQ" else NULL
}

# cancel() with greedy = FALSE or extension = TRUE: the numbers are
# coefficients, so p and q are polynomials in their one symbol over EX
# (dmp_cancel over a field with no associated ring); NULL when there are
# other generators
.cancel_ex <- function(p, q, extension = FALSE) {
  p <- expand(as_sym(p)); q <- expand(as_sym(q))
  g <- Filter(function(z) !.expr_is_number(z), .find_gens(list(p, q)))
  if (length(g) != 1L) return(NULL)
  Fp <- .dp_from(p, g[[1L]], EX); Gp <- .dp_from(q, g[[1L]], EX)
  r <- .dp_inner_gcd(Fp, Gp)
  P <- r[[2L]]; Q <- r[[3L]]
  if (.dom_is_negative(EX, .dp_LC(P)) && .dom_is_negative(EX, .dp_LC(Q))) { P <- .dp_neg(P); Q <- .dp_neg(Q) }
  if (!extension) return(list(.dp_expr(P), .dp_expr(Q)))
  # QQ<a>: each coefficient in its canonical form
  alg <- function(D) .add(lapply(seq_along(D$c), function(i)
    .mul(list(expand(radsimp(.gd_out(D$c[[i]]))), .pow(D$x, .int(length(D$c) - i))))))
  list(alg(P), alg(Q))
}

# Poly.cancel(g, include): dmp_cancel over the unified domain; the result
# keeps the domain of the inputs
.poly_cancel_method <- function(f, g, include = FALSE) {
  ng <- .poly_ngens(f)
  gens <- f$gens[seq_len(ng)]
  dom <- .poly_unified_dom(f, g)
  mk <- function(e) {
    if (!is.null(f$modulus)) return(do.call(Poly, c(list(e), gens, list(modulus = f$modulus))))
    P <- do.call(Poly, c(list(e), gens))
    if (!is.null(dom)) P$dom_str <- dom
    P
  }
  if (!is.null(f$modulus)) {
    r <- .cancel_mod(as_expr(f), as_expr(g), f$modulus)
    if (!is.null(r)) return(.pytuple(list(S_One, mk(r[[1L]]), mk(r[[2L]]))))
  }
  fd <- if (ng == 1L) .poly_frac_dom(f, g)
  if (!is.null(fd)) {
    # dmp_cancel over the fraction field: denominators cleared, the gcd
    # taken over its ring
    cc <- .dp_cancel(.dp_from(as_expr(f), gens[[1L]], fd$K), .dp_from(as_expr(g), gens[[1L]], fd$K))
    mf <- function(D) { P <- do.call(Poly, c(list(.dp_expr(D)), gens)); P$dom_str <- fd$label; P }
    if (!include) return(.pytuple(list(.gd_out(cc[[1L]]), mf(cc[[2L]]), mf(cc[[3L]]))))
    nd <- .as_numer_denom(.gd_out(cc[[1L]]))
    return(.pytuple(list(mf(.dp_mul_ground(cc[[2L]], nd[[1L]])), mf(.dp_mul_ground(cc[[3L]], nd[[2L]])))))
  }
  r <- .poly_cancel_in(as_expr(f), as_expr(g), gens[[1L]])
  if (!include) return(.pytuple(list(r[[1L]], mk(r[[2L]]), mk(r[[3L]]))))
  nd <- .as_numer_denom(r[[1L]])
  .pytuple(list(mk(expand(.mul(list(nd[[1L]], r[[2L]])))), mk(expand(.mul(list(nd[[2L]], r[[3L]]))))))
}

#' Methods of Poly objects
#'
#' `p$cancel(g, include = FALSE)` cancels the common factors of the
#' polynomials `p` and `g`, as SymPy's `Poly.cancel()`: it returns the
#' constant factor and the two reduced polynomials, or with
#' `include = TRUE` the two polynomials with the constant multiplied in.
#' `p$factor_list()` is SymPy's `Poly.factor_list()`: the constant and a
#' list of (factor, multiplicity) pairs with the factors as polynomials.
#' `p$domain` is the coefficient domain. Other names give the components of
#' the object.
#'
#' @param x A polynomial made by [Poly()].
#' @param name A method or component name.
#' @return The method as a function, or the component.
#' @examples
#' x <- Symbol("x")
#' Poly(4*x^2 - 4, x)$cancel(Poly(2*x - 2, x))
#' @export
`$.Poly` <- function(x, name) {
  v <- unclass(x)[[name, exact = FALSE]]
  if (!is.null(v) || name %in% names(unclass(x))) return(v)
  switch(name,
    cancel = function(g, include = FALSE) .poly_cancel_method(x, g, include),
    factor_list = function() .poly_factor_list_method(x),
    domain = {
      d <- .subset2(x, "dom_str")
      if (is.null(d)) {
        d <- .domain_of(.subset2(x, "rep"))
        ng <- .poly_ngens(x)
        if (ng < length(.subset2(x, "gens"))) d <- paste0(d, "[", paste(vapply(.subset2(x, "gens")[-seq_len(ng)], .str, ""), collapse = ","), "]")
      }
      .sympify_domain(d)
    },
    NULL)
}
# a domain from its name, as Poly(..., domain = name)
.sympify_domain <- function(d) {
  if (d %in% c("ZZ", "QQ", "RR", "CC", "EX", "EXRAW")) return(get(d))
  m <- regmatches(d, regexec("^(ZZ|QQ)([[(])(.*)[])]$", d))[[1L]]
  if (length(m)) return(do.call(if (m[3L] == "(") frac_field else poly_ring, c(list(get(m[2L])), lapply(strsplit(m[4L], ",", fixed = TRUE)[[1L]], sympify))))
  d
}

# str() of an element of ring() or field(): StrPrinter._print_PolyElement and
# _print_FracElement
.relem_str <- function(x) {
  R <- x$ring
  K <- if (R$kind == "ring") do.call(poly_ring, c(list(R$domain), R$gens)) else do.call(frac_field, c(list(R$domain), R$gens))
  el <- tryCatch(.de_convert(K, x$expr), error = function(e) NULL)
  if (is.null(el)) return(.str(x$expr))
  gsub("**", "^", .de_str(el), fixed = TRUE)
}
# the quotient a/b of ring or field elements, or NULL when it falls to the
# expression arithmetic
.relem_div <- function(R, a, b) {
  if (is_zero_num(expand(b))) stop(.sympy_condition("ZeroDivisionError", ""))
  if (R$kind != "ring") return(NULL)
  q <- cancel(.div(a, b))
  den <- .as_numer_denom(q)[[2L]]
  if (.dk(R$domain) == "QQ" && is_number_type(den)) return(q)
  NULL
}
