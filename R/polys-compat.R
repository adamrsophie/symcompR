# Polynomial functions brought in line with SymPy: PurePoly (sympy/polys/
# polytools.py), Poly.all_roots with CRootOf (sympy/polys/rootoftools.py),
# interpolating_poly in Lagrange form (sympy/polys/polyfuncs.py), Groebner
# bases sorted by the monomial order (sympy/polys/polytools.py, groebner),
# and gcd and sqf_part of Polys over GF(p).

# gcd(f, g) when f or g is a Poly over GF(p): a Poly with the same modulus,
# as Poly.gcd returns; NULL otherwise
.poly_mod_gcd <- function(f, g, dots) {
  P <- if (inherits(f, "Poly") && !is.null(f$modulus)) f else if (inherits(g, "Poly") && !is.null(g$modulus)) g
  if (is.null(P)) return(NULL)
  pa <- .split_poly_opts(dots)
  gens <- if (length(pa$gens)) pa$gens else P$gens[seq_len(.poly_ngens(P))]
  ex <- function(z) if (inherits(z, "Poly")) as_expr(z) else z
  h <- do.call(gcd, c(list(ex(f), ex(g)), gens, list(modulus = P$modulus)))
  do.call(Poly, c(list(h), gens, list(modulus = P$modulus)))
}

# sqf_part of a Poly: a Poly in the same generators (and modulus)
.sqf_part_poly <- function(f, dots) {
  pa <- .split_poly_opts(dots)
  gens <- f$gens[seq_len(.poly_ngens(f))]
  m <- if (!is.null(pa$opts$modulus)) pa$opts$modulus else f$modulus
  opts <- pa$opts
  if (!is.null(m)) opts$modulus <- m
  r <- do.call(sqf_part, c(list(as_expr(f)), gens, opts))
  do.call(Poly, c(list(r), gens, if (!is.null(m)) list(modulus = m)))
}

# the product of the square-free factors, expanded, with the coefficients in
# the symmetric range under a modulus
.sqf_part_expr <- function(r, pa) {
  r <- expand(r)
  if (is.null(pa$opts$modulus)) return(r)
  as_expr(do.call(Poly, c(list(r), pa$gens, list(modulus = pa$opts$modulus))))
}

# interpolate(): data given as a list of (x, y) pairs rather than y values
.interp_is_pairs <- function(points) {
  p <- points[[1L]]
  length(p) == 2L && ((is.atomic(p) && !is.object(p)) || (is.list(p) && !is_sym(p)))
}

# interpolating_poly(n, x, X, Y): the Lagrange form, unexpanded
.lagrange_poly <- function(n, x, X, Y) {
  terms <- lapply(seq_len(n), function(i) {
    o <- setdiff(seq_len(n), i)
    numer <- .mul(lapply(o, function(j) .sub(x, X[[j]])))
    denom <- .mul(lapply(o, function(j) .sub(X[[i]], X[[j]])))
    .mul(list(.div(numer, denom), Y[[i]]))
  })
  .add(terms)
}

# the reduced Groebner basis sorted by leading monomial, largest first in
# the requested order
.gb_sorted <- function(G, order) {
  if (length(G) < 2L) return(G)
  lms <- do.call(rbind, lapply(G, function(p) qp_lm(qp_reorder(p, order))))
  G[.mono_order_index(lms, order)]
}

# all_roots(f): the roots over the irreducible factors of f in ZZ[x], with
# multiplicities. Linear and quadratic factors and binomials are solved by
# radicals; every other factor gives CRootOf objects. Real roots come first
# in increasing order; the complex ones follow factor by factor, the factors
# in SymPy's ordered() order, each by real then imaginary part. NULL when f
# is not a univariate polynomial with rational coefficients.
.all_roots_rootof <- function(f, x) {
  if (inherits(f, "Poly")) {
    if (is.null(x)) x <- f$gens[[1L]]
    f <- as_expr(f)
  }
  f <- as_sym(f)
  if (is.null(x)) {
    fs <- free_symbols(f)
    if (length(fs) != 1L) return(NULL)
    x <- fs[[1L]]
  }
  cf <- tryCatch(.coeff_list(f, as_sym(x)), error = function(e) NULL)
  if (length(cf) < 2L || !all(vapply(cf, is_exact_number, TRUE))) return(NULL)
  q <- do.call(c, lapply(cf, `[[`, "value"))
  zl <- .zp_factor_list(gmp::numerator(q * gmp::as.bigq(.lcm_den(q))))
  fx <- list()
  rts <- list()
  mk <- character()
  mv <- integer()
  for (pr in zl$factors) {
    g <- lapply(as.list(pr[[1L]]), function(v) make_rational(gmp::as.bigq(v)))
    n <- length(g) - 1L
    binom <- n > 2L && all(vapply(g[-c(1L, n + 1L)], is_zero_num, TRUE))
    rs <- if (n <= 2L || binom) .roots_of_coeffs(g) else {
      xx <- Dummy("x")
      pe <- .coeffs_to_expr(g, xx)
      lapply(seq_len(n) - 1L, function(i) CRootOf(pe, xx, i))
    }
    fx[[length(fx) + 1L]] <- .coeffs_to_expr(g, Symbol("x"))
    rts[[length(rts) + 1L]] <- rs
    mk <- c(mk, vapply(rs, .k, ""))
    mv <- c(mv, rep(as.integer(pr[[2L]]), length(rs)))
  }
  isreal <- function(r) {
    v <- tryCatch(as.complex(.neval(r, list())), error = function(e) NA)
    !is.na(v) && abs(Im(v)) < 1e-12 * max(1, abs(Re(v)))
  }
  out <- .sort_roots(Filter(isreal, unlist(rts, recursive = FALSE)))
  o <- if (length(fx) > 1L) match(vapply(.bool_print_order(fx), .k, ""), vapply(fx, .k, "")) else seq_along(fx)
  for (k in o) {
    cr <- Filter(Negate(isreal), rts[[k]])
    # CRootOf indices already follow SymPy's order
    if (!all(vapply(cr, function(r) r$type == "CRootOf", NA))) cr <- .sort_roots(cr)
    out <- c(out, cr)
  }
  out <- unlist(lapply(out, function(r) rep(list(r), mv[match(.k(r), mk)])), recursive = FALSE)
  if (is.null(out)) list() else out
}

# PurePoly(f, ...): a Poly whose generators are the placeholders _0, _1, ...
# (skipping names the coefficients use), so polynomials that differ only in
# their generators compare equal
.purepoly <- function(f, ...) {
  P <- if (inherits(f, "Poly")) f else Poly(f, ...)
  ng <- .poly_ngens(P)
  used <- unlist(lapply(P$gens[-seq_len(ng)], function(g) vapply(free_symbols(g), .str, "")))
  ph <- character()
  k <- 0L
  while (length(ph) < ng) {
    s <- paste0("_", k)
    if (!(s %in% used)) ph <- c(ph, s)
    k <- k + 1L
  }
  P <- unclass(P)
  P$gens[seq_len(ng)] <- lapply(ph, Symbol)
  structure(P, class = c("PurePoly", "Poly"))
}

#' @export
print.PurePoly <- function(x, ...) {
  out <- utils::capture.output(NextMethod())
  cat(sub("^Poly\\(", "PurePoly(", out), sep = "\n")
  invisible(x)
}

#' @export
Ops.PurePoly <- function(e1, e2) {
  r <- NextMethod()
  if (inherits(r, "Poly") && !inherits(r, "PurePoly")) class(r) <- c("PurePoly", "Poly")
  r
}

#' @export
chooseOpsMethod.PurePoly <- function(x, y, mx, my, cl, reverse) TRUE

# p$factor_list(): Poly.factor_list, the constant and the irreducible factors
# as Polys in the generators (and modulus) of p
.poly_factor_list_method <- function(P) {
  gens <- .subset2(P, "gens")[seq_len(.poly_ngens(P))]
  m <- .subset2(P, "modulus")
  opt <- if (!is.null(m)) list(modulus = m) else if (identical(.subset2(P, "dom_str"), "QQ")) list(domain = "QQ")
  fl <- do.call(factor_list, c(list(as_expr(P)), gens, if (!is.null(m)) opt))
  .pytuple(list(fl$coeff, lapply(fl$factors, function(f) .pytuple(list(do.call(Poly, c(list(f[[1L]]), gens, opt)), f[[2L]])))))
}
