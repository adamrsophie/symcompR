# Rational function manipulation: fraction, together, cancel, apart, after
# sympy/simplify/radsimp.py, sympy/polys/rationaltools.py and
# sympy/polys/partfrac.py.

#' Numerator and denominator
#'
#' `fraction(e)` returns `list(numer, denom)` by collecting factors with
#' negative exponents, like SymPy's `fraction()`; it does not combine sums.
#' `numer()` and `denom()` return one part.
#'
#' @param e An expression.
#' @param exact If `TRUE`, only powers with a negative numeric exponent go to
#'   the denominator and the parts are not re-evaluated.
#' @return A list of two expressions.
#' @examples
#' x <- Symbol("x"); y <- Symbol("y")
#' fraction(x / y)
#' @export
fraction <- function(e, exact = FALSE) {
  e <- as_sym(e)
  num <- list(); den <- list()
  for (term in .margs(e)) {
    if (.is_comm(term) && (term$type == "Pow" || (term$type == "Function" && term$name == "exp"))) {
      be <- if (term$type == "Pow") term$args else list(S_E, term$args[[1L]])
      b <- be[[1L]]; ex <- be[[2L]]
      if (exact) {
        if (!length(free_symbols(ex)) && could_extract_minus_sign(ex)) {
          den[[length(den) + 1L]] <- if (is_negone_num(ex)) b else new_sym("Pow", list(b, .neg(ex)))
        } else num[[length(num) + 1L]] <- term
      } else if (isTRUE(ask_fact(ex, "negative"))) {
        den[[length(den) + 1L]] <- if (is_negone_num(ex)) b else .pow(b, .neg(ex))
      } else if (isTRUE(ask_fact(ex, "positive"))) {
        num[[length(num) + 1L]] <- term
      } else if (ex$type == "Mul") {
        nd <- .as_numer_denom(term)
        if (!is_one_num(nd[[1L]])) num[[length(num) + 1L]] <- nd[[1L]]
        den[[length(den) + 1L]] <- nd[[2L]]
      } else num[[length(num) + 1L]] <- term
    } else if (term$type == "Rational") {
      if (gmp::numerator(term$value) != 1) num[[length(num) + 1L]] <- .zs(gmp::numerator(term$value))
      den[[length(den) + 1L]] <- .zs(gmp::denominator(term$value))
    } else num[[length(num) + 1L]] <- term
  }
  # Mul(*numer, evaluate=not exact)
  build <- function(l) {
    if (!exact) return(.mul(l))
    if (!length(l)) S_One else if (length(l) == 1L) l[[1L]] else new_sym("Mul", l)
  }
  list(build(num), build(den))
}

#' @rdname fraction
#' @export
numer <- function(e) fraction(e)[[1L]]

#' @rdname fraction
#' @export
denom <- function(e) fraction(e)[[2L]]

#' Combine a sum of fractions
#'
#' `together()` rewrites a sum of rational terms over a common denominator
#' without expanding, like SymPy's `together()`.
#'
#' @param e An expression.
#' @param deep Also combine inside function arguments.
#' @param fraction Put the result over a common denominator.
#' @return An expression.
#' @examples
#' x <- Symbol("x"); y <- Symbol("y")
#' together(1/x + 1/y)
#' @export
together <- function(e, deep = FALSE, fraction = TRUE) {
  if (inherits(e, "symmatrix")) return(.mat_map(e, function(z) together(z, deep, fraction)))
  .together_rec(as_sym(e), deep, fraction)
}

# SymPy together._together
.together_rec <- function(e, deep, fraction) {
  if (!length(e$args) || (e$type == "Function" && !deep)) return(e)
  rec <- function(a) .together_rec(a, deep, fraction)
  if (e$type == "Add") return(gcd_terms(lapply(e$args, rec), fraction = fraction))
  if (e$type == "Pow") return(.pow(rec(e$args[[1L]]), if (deep) rec(e$args[[2L]]) else e$args[[2L]]))
  sym_rebuild(e, lapply(e$args, rec))
}

#' Cancel common factors of a rational function
#'
#' Puts `e` in the form `p/q` with `p` and `q` expanded polynomials with no
#' common factor, like SymPy's `cancel()`.
#'
#' @param e An expression.
#' @param ... Generators.
#' @return An expression.
#' @examples
#' x <- Symbol("x")
#' cancel((x^2 - 1) / (x - 1))
#' @export
cancel <- function(e, ...) {
  if (is.list(e) && !is_sym(e) && !inherits(e, c("symmatrix", "Poly"))) return(.cancel_tuple(e, .as_gens(Filter(Negate(is.logical), list(...))), list(...)$polys, list(...)))
  if (is_sym(e) && !is.null(r <- .cancel_expr_opts(e, list(...)))) return(r)
  if (inherits(e, "symmatrix")) return(.mat_map(e, cancel))
  .cancel(as_sym(e), .as_gens(list(...)))
}

# SymPy cancel(): signsimp, factor_terms(radical = TRUE), then the numerator
# and denominator cancelled as polynomials (PolyElement.cancel)
.cancel <- function(f, gens = NULL, signsimp_ = TRUE) {
  if (signsimp_) f <- signsimp(f)
  if (is_number_type(f) || is_relational(f) || !.is_expr_node(f)) return(f)
  f <- factor_terms(f, radical = TRUE)
  nd <- .as_numer_denom(f)
  if (.has_type(f, "Piecewise") || !.is_comm(f)) return(.cancel_fallback(f))
  p <- expand(nd[[1L]]); q <- expand(nd[[2L]])
  if (is.null(gens)) {
    gens <- .find_gens(list(p, q))
    # I is a coefficient (the ZZ_I and QQ_I domains), not a generator
    gi <- vapply(gens, function(g) identical(g$key, S_I$key), TRUE)
    if (any(gi)) {
      gens <- gens[!gi]
      if (!length(gens)) return(expand(f))
      if (length(gens) == 1L) {
        gp <- .gauss_polys(p, q, gens[[1L]])
        if (!is.null(gp)) return(.gauss_cancel_expr(gp[[1L]], gp[[2L]], gens[[1L]]))
      }
      gens <- c(gens, list(S_I))
    }
  }
  if (!length(gens)) return(expand(f))
  cv <- tryCatch(.to_qpolys(list(p, q), gens), error = function(err) NULL)
  if (is.null(cv)) return(.cancel_fallback(f))
  real <- .has_type(nd[[1L]], "Float") || .has_type(nd[[2L]], "Float")
  pq <- .qp_cancel(cv$polys[[1L]], cv$polys[[2L]], field = real)
  if (real) {
    # the RR domain: coefficients come back as Floats
    fl <- function(P) .add(lapply(seq_along(P$c), function(i)
      .mul(c(list(N(make_rational(P$c[i]))), lapply(which(P$e[i, ] > 0L), function(j) .pow(cv$gens[[j]], .int(P$e[i, j])))))))
    return(.div(fl(pq[[1L]]), fl(pq[[2L]])))
  }
  ii <- which(vapply(cv$gens, function(g) identical(g$key, S_I$key), TRUE))
  if (length(ii)) return(.div(.qp_gauss_expr(pq[[1L]], cv$gens, ii), .qp_gauss_expr(pq[[2L]], cv$gens, ii)))
  .div(qpoly_to_expr(pq[[1L]], cv$gens), qpoly_to_expr(pq[[2L]], cv$gens))
}

# A qpoly with I among its generators read back over QQ_I: powers of I
# reduced by I^2 = -1 and the terms of each monomial in the other generators
# joined into one Gaussian coefficient, as PolyElement.as_expr does there.
# ponytail: the gcd still treats I as a variable, so a common Gaussian factor
# of several terms is not cancelled; port PolyRing over ZZ_I if that matters
.qp_gauss_expr <- function(P, gens, ii) {
  if (qp_is_zero(P)) return(S_Zero)
  E <- P$e[, -ii, drop = FALSE]
  k <- P$e[, ii] %% 4L
  key <- apply(E, 1L, paste, collapse = ",")
  terms <- list()
  for (u in unique(key)) {
    w <- which(key == u)
    re <- gmp::as.bigq(0); im <- gmp::as.bigq(0)
    for (i in w) {
      v <- P$c[i]
      if (k[i] == 0L) re <- re + v else if (k[i] == 1L) im <- im + v
      else if (k[i] == 2L) re <- re - v else im <- im - v
    }
    if (re == 0 && im == 0) next
    cf <- .add(list(make_rational(re), .mul(list(make_rational(im), S_I))))
    mono <- lapply(which(E[w[1L], ] > 0L), function(j) .pow(gens[-ii][[j]], .int(E[w[1L], j])))
    terms[[length(terms) + 1L]] <- .mul(c(list(cf), mono))
  }
  .add(terms)
}

# the noncommutative and Piecewise branch of cancel()
.cancel_fallback <- function(f) {
  if (f$type %in% c("Add", "Mul")) {
    keep <- vapply(f$args, function(a) .is_comm(a) && !.has_type(a, "Piecewise"), TRUE)
    op <- if (f$type == "Add") .add else .mul
    return(op(c(list(.cancel(op(f$args[keep]))), lapply(f$args[!keep], .cancel))))
  }
  if (!length(f$args)) return(f)
  sym_rebuild(f, lapply(f$args, function(a) if (.is_expr_node(a) && !(a$type %in% c("Tuple", "BooleanAtom"))) .cancel(a) else a))
}

# Poly(n, z).cancel(Poly(d, z), include = False): the common factor removed,
# the content and the signs of the leading coefficients in z collected in a
# rational c; returns list(c, numerator, denominator) as expressions
.poly_cancel_in <- function(n, d, z, joint = FALSE) {
  n <- expand(n); d <- expand(d)
  if (is_zero_num(n)) return(list(S_One, S_Zero, S_One))
  gp <- .gauss_polys(n, d, z)
  if (!is.null(gp)) return(.gauss_poly_cancel(gp[[1L]], gp[[2L]], z))
  cl <- tryCatch(list(.coeff_list(n, z), .coeff_list(d, z)), error = function(err) NULL)
  if (!is.null(cl)) {
    ex <- if (joint) .uses_ex_domain(c(cl[[1L]], cl[[2L]])) else .uses_ex_domain(cl[[1L]]) || .uses_ex_domain(cl[[2L]])
    if (ex) return(.ex_cancel(cl[[1L]], cl[[2L]], n, d, z))
  }
  gens <- c(list(z), Filter(function(g) !identical(g$key, z$key), .find_gens(list(n, d))))
  cv <- tryCatch(.to_qpolys(list(n, d), gens), error = function(err) NULL)
  if (is.null(cv)) return(list(S_One, n, d))
  p <- cv$polys[[1L]]; q <- cv$polys[[2L]]
  cq <- .bq_lcm_den(p); cp <- .bq_lcm_den(q)
  f <- qp_scale(p, gmp::as.bigq(cq)); g <- qp_scale(q, gmp::as.bigq(cp))
  h <- qp_gcd(f, g)
  if (!qp_is_const(h)) { h <- qp_primitive(h)$poly; f <- qp_exquo(f, h); g <- qp_exquo(g, h) }
  k <- gmp::gcd.bigz(.bq_int_content(f), .bq_int_content(g))
  if (k != 1) { f <- qp_scale(f, gmp::as.bigq(1, k)); g <- qp_scale(g, gmp::as.bigq(1, k)) }
  G <- gmp::gcd.bigz(cp, cq)
  c <- gmp::as.bigq(cp / G) / gmp::as.bigq(cq / G)
  pn <- qp_lc(f) < 0; qn <- qp_lc(g) < 0
  if (pn && qn) { f <- qp_neg(f); g <- qp_neg(g) }
  else if (pn) { c <- -c; f <- qp_neg(f) }
  else if (qn) { c <- -c; g <- qp_neg(g) }
  list(make_rational(c), .poly_as_expr_in(f, cv$gens), .poly_as_expr_in(g, cv$gens))
}

# construct_domain falls back to EX: an irrational algebraic number among
# numeric coefficients, or composite generators that are algebraic numbers or
# share free symbols (there could be relations between them)
.uses_ex_domain <- function(coeffs) {
  coeffs <- Filter(function(c) !is_zero_num(c), coeffs)
  if (!length(coeffs)) return(FALSE)
  if (all(vapply(coeffs, function(c) !length(free_symbols(c)), TRUE))) {
    simple <- vapply(coeffs, function(c) {
      if (c$type %in% c("Integer", "Rational", "Float")) return("rat")
      if (!is.null(.gauss_num(c))) return("rat")
      if (isTRUE(ask_fact(c, "algebraic"))) return("alg")
      "other"
    }, "")
    if (!any(simple == "other")) return(any(simple == "alg"))
  }
  nds <- unlist(lapply(coeffs, .as_numer_denom), recursive = FALSE)
  gens <- .find_gens(nds)
  if (!length(gens)) return(FALSE)
  if (any(vapply(gens, function(g) !length(free_symbols(g)) && isTRUE(ask_fact(g, "algebraic")), TRUE))) return(TRUE)
  seen <- character()
  for (g in gens) {
    s <- vapply(free_symbols(g), `[[`, "", "key")
    if (length(intersect(seen, s))) return(TRUE)
    seen <- c(seen, s)
  }
  FALSE
}

# arithmetic in EX: every result is cancel()ed and expanded (mul only)
.ex_simp <- function(x) expand(cancel(x), deep = FALSE, mul = TRUE, multinomial = FALSE,
                               power_exp = FALSE, power_base = FALSE, log = FALSE)
.exl_strip <- function(p) { while (length(p) && is_zero_num(p[[1L]])) p <- p[-1L]; p }
.exl_sub <- function(f, g) {
  if (!length(f)) return(lapply(g, function(c) .ex_simp(.neg(c))))
  if (!length(g)) return(f)
  df <- length(f); dg <- length(g)
  if (df == dg) return(.exl_strip(Map(function(a, b) .ex_simp(.sub(a, b)), f, g)))
  k <- abs(df - dg)
  if (df > dg) { h <- f[seq_len(k)]; f <- f[-seq_len(k)] }
  else { h <- lapply(g[seq_len(k)], function(c) .ex_simp(.neg(c))); g <- g[-seq_len(k)] }
  c(h, Map(function(a, b) .ex_simp(.sub(a, b)), f, g))
}
.exl_mul_term <- function(f, c, i) if (is_zero_num(c) || !length(f)) list() else
  c(lapply(f, function(cf) .ex_simp(.mul(list(cf, c)))), rep(list(S_Zero), i))
.exl_add_term <- function(f, c, i) {
  if (is_zero_num(c)) return(f)
  n <- length(f); m <- n - 1L
  if (i == m) return(.exl_strip(c(list(.ex_simp(.add(list(f[[1L]], c)))), f[-1L])))
  if (i >= n) return(c(list(c), rep(list(S_Zero), i - n), f))
  f[[n - i]] <- .ex_simp(.add(list(f[[n - i]], c)))
  f
}
# dup_ff_div over EX
.exl_div <- function(f, g) {
  df <- length(f) - 1L; dg <- length(g) - 1L
  q <- list(); r <- f; dr <- df
  if (df < dg) return(list(q = q, r = r))
  lcg <- g[[1L]]
  repeat {
    c <- .ex_simp(.div(r[[1L]], lcg))
    j <- dr - dg
    q <- .exl_add_term(q, c, j)
    r <- .exl_sub(r, .exl_mul_term(g, c, j))
    odr <- dr; dr <- length(r) - 1L
    if (dr < dg) break
    if (!(dr < odr)) stop("polynomial division failed over EX")
  }
  list(q = q, r = r)
}
.exl_as_expr <- function(p, z) {
  d <- length(p) - 1L
  .add(lapply(seq_along(p), function(k) .mul(list(p[[k]], .pow(z, .int(d - k + 1L))))))
}

# dmp_cancel over EX (a field with no associated ring): divide by the monic
# gcd in z with EX arithmetic, then fix the signs of the leading coefficients
.ex_cancel <- function(fc, gc, n, d, z) {
  h <- list(S_One)
  if (length(fc) > 1L && length(gc) > 1L) {
    G <- tryCatch(gcd(n, d, z), error = function(err) S_One)
    hc <- tryCatch(.coeff_list(G, z), error = function(err) list(S_One))
    if (length(hc) > 1L) h <- lapply(hc, function(c) .ex_simp(.div(c, hc[[1L]])))
  }
  p <- .exl_div(fc, h)$q; q <- .exl_div(gc, h)$q
  cp <- S_One
  pn <- could_extract_minus_sign(p[[1L]]); qn <- could_extract_minus_sign(q[[1L]])
  neg <- function(l) lapply(l, function(c) .ex_simp(.neg(c)))
  if (pn && qn) { p <- neg(p); q <- neg(q) }
  else if (pn) { cp <- S_NegativeOne; p <- neg(p) }
  else if (qn) { cp <- S_NegativeOne; q <- neg(q) }
  list(cp, .exl_as_expr(p, z), .exl_as_expr(q, z))
}

# Poly(..., z).as_expr() over a polynomial ring: sum of coefficient*z^k with
# each coefficient a polynomial in the other generators
.poly_as_expr_in <- function(p, gens) {
  if (qp_is_zero(p)) return(S_Zero)
  ks <- sort(unique(p$e[, 1L]), decreasing = TRUE)
  rest <- gens[-1L]
  .add(lapply(ks, function(k) {
    rows <- which(p$e[, 1L] == k)
    sub <- qp(p$e[rows, -1L, drop = FALSE], p$c[rows])
    .mul(list(qpoly_to_expr(sub, rest), .pow(gens[[1L]], .int(k))))
  }))
}

.bq_lcm_den <- function(p) {
  if (qp_is_zero(p)) return(gmp::as.bigz(1))
  Reduce(gmp::lcm.bigz, as.list(gmp::denominator(p$c)), gmp::as.bigz(1))
}
.bq_int_content <- function(p) {
  if (qp_is_zero(p)) return(gmp::as.bigz(0))
  Reduce(gmp::gcd.bigz, as.list(abs(gmp::numerator(p$c))))
}

# PolyElement.cancel over ZZ or QQ: clear denominators, divide by the gcd over
# ZZ (integer content included), and make the denominator's leading
# coefficient positive
.qp_cancel <- function(p, q, field = FALSE) {
  if (qp_is_zero(p)) return(list(p, qp_const(1, qp_nvars(q))))
  if (field) {
    # a field with no associated ring (RR): divide by the monic gcd
    h <- qp_monic(qp_gcd(p, q))
    if (!qp_is_const(h)) { p <- qp_exquo(p, h); q <- qp_exquo(q, h) }
    if (qp_lc(q) < 0) { p <- qp_neg(p); q <- qp_neg(q) }
    return(list(p, q))
  }
  cq <- .bq_lcm_den(p); cp <- .bq_lcm_den(q)
  f <- qp_scale(p, gmp::as.bigq(cq)); g <- qp_scale(q, gmp::as.bigq(cp))
  h <- qp_gcd(f, g)
  if (!qp_is_const(h)) {
    h <- qp_primitive(h)$poly
    f <- qp_exquo(f, h); g <- qp_exquo(g, h)
  }
  k <- gmp::gcd.bigz(.bq_int_content(f), .bq_int_content(g))
  if (k != 1) { f <- qp_scale(f, gmp::as.bigq(1, k)); g <- qp_scale(g, gmp::as.bigq(1, k)) }
  G <- gmp::gcd.bigz(cp, cq)
  f <- qp_scale(f, gmp::as.bigq(cp / G)); g <- qp_scale(g, gmp::as.bigq(cq / G))
  if (qp_lc(g) < 0) { f <- qp_neg(f); g <- qp_neg(g) }
  list(f, g)
}

#' Partial fraction decomposition
#'
#' `apart(e, x)` decomposes a rational function of `x` into partial fractions
#' over the rationals, using undetermined coefficients as SymPy does.
#' With `full = TRUE` the decomposition runs over the roots of the
#' denominator factors (Bronstein's full decomposition), leaving `RootSum`
#' terms for irreducible factors. `apart_list()` returns that decomposition in
#' SymPy's structured form: the common factor, the polynomial part as a
#' [Poly()], and for each term the root polynomial in a dummy `w`, the
#' numerator and denominator as functions of a root, and the exponent.
#' `assemble_partfrac_list()` rebuilds the expression from that list.
#'
#' @param e A rational function.
#' @param x The variable.
#' @param full Full decomposition over the roots of the denominator.
#' @param partial_list A result of `apart_list()`.
#' @return An expression.
#' @examples
#' x <- Symbol("x")
#' apart(1 / (x^2 - 1), x)
#' apart((x^2 + 1) / (x^3 - x), x)
#' @export
apart <- function(e, x = NULL, full = FALSE) {
  e <- as_sym(e)
  if (is.null(x)) {
    fs <- free_symbols(e)
    if (length(fs) != 1L) stop("apart() needs the variable for multivariate input")
    x <- fs[[1L]]
  }
  x <- as_sym(x)
  if (isTRUE(full)) {
    pl <- apart_list(e, x)
    if (is_sym(pl)) return(pl)
    full_terms <- assemble_partfrac_list(list(S_One, pl[[2L]], pl[[3L]]))
    poly <- as_expr(pl[[2L]])
    rest <- .sub(full_terms, poly)
    terms <- lapply(if (rest$type == "Add") rest$args else list(rest), function(t) if (.has_type(t, "RootSum")) t else sym_factor(t))
    return(.mul(list(pl[[1L]], .add(c(list(poly), terms)))))
  }
  if (!length(e$args)) return(e)
  # P.cancel(Q), P.div(Q), then partial fractions of the remainder, each term
  # factored, all times the common factor
  nd <- .as_numer_denom(e)
  cc <- .poly_cancel_in(nd[[1L]], nd[[2L]], x, joint = TRUE)
  common <- cc[[1L]]; P <- cc[[2L]]; Q <- cc[[3L]]
  r <- .apart_terms(P, Q, x)
  terms <- lapply(r$terms, function(t) if (.has_type(t, "RootSum")) t else sym_factor(t))
  .mul(list(common, .add(c(list(r$poly), terms))))
}

#' @rdname apart
#' @export
apart_list <- function(e, x = NULL) {
  e <- as_sym(e)
  if (!length(e$args)) return(e)
  nd <- .as_numer_denom(e)
  gens <- if (is.null(x)) .gens_of(nd, list()) else list(as_sym(x))
  if (length(gens) != 1L) stop("multivariate partial fraction decomposition", call. = FALSE)
  x <- gens[[1L]]
  P <- expand(nd[[1L]]); Q <- expand(nd[[2L]])
  # Poly.cancel: common factor and coprime numerator and denominator
  g <- gcd(P, Q, x)
  P <- cancel(.div(P, g)); Q <- cancel(.div(Q, g))
  # dmp_cancel: over QQ the denominators are cleared into the common factor
  # (over ZZ nothing is taken out), and the denominator gets a positive
  # leading coefficient
  dP <- .apart_den(P, x); dQ <- .apart_den(Q, x)
  common <- .div(dQ, dP)
  P <- expand(.mul(list(P, dP))); Q <- expand(.mul(list(Q, dQ)))
  lp <- could_extract_minus_sign(LC(P, x)); lq <- could_extract_minus_sign(LC(Q, x))
  if (lq) { Q <- .neg(Q); if (lp) P <- .neg(P) else common <- .neg(common) }
  else if (lp) { P <- .neg(P); common <- .neg(common) }
  dv <- div(P, Q, x)
  polypart <- Poly(dv[[1L]], x)
  P <- expand(dv[[2L]])
  list(common, polypart, .apart_full_decomposition(P, Q, x))
}

# The lcm of the denominators of the rational coefficients of f in x.
.apart_den <- function(f, x) {
  cs <- Filter(function(c) c$type == "Rational", .coeff_list(f, x))
  if (!length(cs)) return(S_One)
  make_rational(gmp::as.bigq(Reduce(gmp::lcm.bigz, lapply(cs, function(c) gmp::denominator(gmp::as.bigq(c$value))))))
}

# apart_list_full_decomposition (Bronstein's algorithm)
.apart_full_decomposition <- function(P, Q, x) {
  f <- .div(P, Q)
  u <- Function("u")(x)
  a <- Dummy("a"); w <- Dummy("w")
  U <- list(); partial <- list()
  sl <- tryCatch(sqf_list(Q, x), error = function(err) {
    # coefficients with parameters: a squarefree denominator is its own list
    if (degree(gcd(Q, diff(Q, x), x), x) > 0) stop(err)
    list(coeff = S_One, factors = list(list(Q, 1L)))
  })
  # sqf_list_include(all = TRUE): every multiplicity, the coefficient in the first factor
  mx <- if (length(sl$factors)) max(vapply(sl$factors, function(p) as.integer(p[[2L]]), 0L)) else 1L
  fl <- lapply(seq_len(mx), function(k) {
    hit <- Filter(function(p) as.integer(p[[2L]]) == k, sl$factors)
    list(if (length(hit)) hit[[1L]][[1L]] else S_One, k)
  })
  fl[[1L]][[1L]] <- expand(.mul(list(sl$coeff, fl[[1L]][[1L]])))
  for (dn in fl) {
    d <- dn[[1L]]; n <- dn[[2L]]
    b <- d
    U[[length(U) + 1L]] <- if (n == 1L) u else diff(u, x, n - 1L)
    h <- .div(cancel(.mul(list(f, .pow(b, .int(n))))), .pow(u, .int(n)))
    H <- list(h); subsl <- list()
    for (j in seq_len(n - 1L)) H[[j + 1L]] <- .div(diff(H[[j]], x), .int(j))
    for (j in seq_len(n)) subsl[[j]] <- list(U[[j]], .div(diff(b, x, j), .int(j)))
    for (j in 0:(n - 1L)) {
      nd <- .as_numer_denom(cancel(H[[j + 1L]]))
      Pj <- nd[[1L]]; Qj <- nd[[2L]]
      for (i in 0:j) Pj <- subs(Pj, subsl[[j - i + 1L]][[1L]], subsl[[j - i + 1L]][[2L]])
      Qj <- subs(Qj, subsl[[1L]][[1L]], subsl[[1L]][[2L]])
      Pj <- expand(Pj); Qj <- expand(Qj)
      G <- gcd(Pj, d, x)
      D <- quo(d, G, x)
      hg <- half_gcdex(Qj, D, x)
      bb <- rem(expand(.mul(list(Pj, quo(hg$s, hg$h, x)))), D, x)
      Dw <- Poly(subs(D, x, w), w)
      numer <- Lambda(a, subs(bb, x, a))
      denom <- Lambda(a, .sub(x, a))
      partial[[length(partial) + 1L]] <- list(Dw, numer, denom, .int(n - j))
    }
  }
  partial
}

#' @rdname apart
#' @export
assemble_partfrac_list <- function(partial_list) {
  common <- as_sym(partial_list[[1L]])
  pfd <- as_expr(partial_list[[2L]])
  for (t in partial_list[[3L]]) {
    r <- t[[1L]]; nf <- t[[2L]]; df <- t[[3L]]; ex <- t[[4L]]
    if (inherits(r, "Poly")) {
      an <- nf$args[[1L]]$args[[1L]]; nu <- nf$args[[2L]]
      ad <- df$args[[1L]]$args[[1L]]; de <- subs(df$args[[2L]], ad, an)
      w <- r$gens[[1L]]
      pfd <- .add(list(pfd, RootSum(as_expr(r), Lambda(an, .div(nu, .pow(de, ex))), w, auto = FALSE)))
    } else {
      for (root in r) pfd <- .add(list(pfd, .div(apply_lambda(nf, root), .pow(apply_lambda(df, root), ex))))
    }
  }
  .mul(list(common, pfd))
}

# .apart_terms(P, Q, x) for P/Q, or .apart_terms(f, x) for a rational function f
.apart_terms <- function(P, Q, x) {
  if (missing(x)) {
    x <- Q
    nd <- .as_numer_denom(cancel(P)); P <- nd[[1L]]; Q <- nd[[2L]]
  }
  P <- expand(P); Q <- expand(Q)
  # P.div(Q) with a constant Q: the quotient as Poly.as_expr, grouped by powers of x
  if (!.has1(Q, x)) return(list(poly = .coeffs_to_expr(lapply(.coeff_list(P, x), function(c) expand(.div(c, Q))), x), terms = list()))
  pc <- .coeff_list(P, x)
  qc <- .coeff_list(Q, x)
  # polynomial division with expression coefficients
  dv <- .expr_poly_divmod(pc, qc)
  poly_part <- .coeffs_to_expr(dv$q, x)
  rem <- dv$r
  if (!length(rem) || all(vapply(rem, is_zero_num, TRUE))) return(list(poly = poly_part, terms = list()))
  # Q.degree() <= 1: the remainder over Q is the partial fraction
  if (length(qc) <= 2L) return(list(poly = poly_part, terms = list(.div(.coeffs_to_expr(rem, x), Q))))
  fl <- factor_list(Q, x)
  lcq <- fl$coeff
  facs <- fl$factors
  # Unknown numerators: for factor f^k of degree d, k numerators of degree d-1.
  unknowns <- list()
  terms <- list()
  idx <- 0L
  for (pr in facs) {
    f <- pr[[1L]]
    k <- pr[[2L]]
    d <- length(.coeff_list(f, x)) - 1L
    for (j in seq_len(k)) {
      cs <- lapply(seq_len(d), function(t) {
        idx <<- idx + 1L
        Dummy(paste0("A", idx))
      })
      unknowns <- c(unknowns, cs)
      numer_j <- .add(lapply(seq_len(d), function(t) .mul(list(cs[[t]], .pow(x, .int(d - t))))))
      terms[[length(terms) + 1L]] <- list(numer = numer_j, f = f, j = j)
    }
  }
  # sum numer_j * Q / (lc * f^j) must equal remainder polynomial
  R <- .coeffs_to_expr(rem, x)
  lhs <- .add(lapply(terms, function(t) {
    .mul(list(t$numer, cancel(.div(Q, .mul(list(lcq, .pow(t$f, .int(t$j))))))))
  }))
  eqs_expr <- .coeff_list(expand(.sub(lhs, .div(R, lcq))), x)
  sol <- .linsolve_exprs(eqs_expr, unknowns)
  if (is.null(sol)) stop("apart() could not solve for the partial fraction coefficients")
  out <- lapply(terms, function(t) {
    num <- xreplace(t$numer, sol)
    .div(cancel(num), .pow(t$f, .int(t$j)))
  })
  out <- Filter(function(z) !is_zero_num(z), out)
  list(poly = poly_part, terms = out)
}

.coeffs_to_expr <- function(cl, x) {
  if (!length(cl)) return(S_Zero)
  d <- length(cl) - 1L
  .add(lapply(seq_along(cl), function(k) .mul(list(cl[[k]], .pow(x, .int(d - k + 1L))))))
}

# Division of dense coefficient lists with expression coefficients.
.expr_poly_divmod <- function(a, b) {
  db <- length(b) - 1L
  da <- length(a) - 1L
  if (da < db) return(list(q = list(), r = a))
  r <- a
  q <- rep(list(S_Zero), da - db + 1L)
  lb <- b[[1L]]
  for (i in seq_along(q)) {
    cf <- cancel(.div(r[[i]], lb))
    q[[i]] <- cf
    if (!is_zero_num(cf)) {
      for (j in 0:db) r[[i + j]] <- expand(.sub(r[[i + j]], .mul(list(b[[j + 1L]], cf))))
    }
  }
  rr <- if (db > 0L) r[(length(r) - db + 1L):length(r)] else list()
  while (length(rr) && is_zero_num(cancel(rr[[1L]]))) rr <- rr[-1L]
  list(q = q, r = rr)
}

# Solve a linear system given as expressions equal to zero; returns a rule
# list of list(unknown, value) pairs, or NULL when inconsistent.
.linsolve_exprs <- function(eqs, unknowns) {
  n <- length(unknowns)
  rows <- list()
  for (eq in eqs) {
    eq <- expand(eq)
    if (is_zero_num(eq)) next
    coefs <- lapply(unknowns, function(u) .linear_coeff(eq, u))
    rest <- expand(.sub(eq, .add(lapply(seq_len(n), function(i) .mul(list(coefs[[i]], unknowns[[i]]))))))
    rows[[length(rows) + 1L]] <- c(coefs, list(.neg(rest)))
  }
  sol <- .gauss_solve(rows, n)
  if (is.null(sol)) return(NULL)
  lapply(seq_len(n), function(i) list(unknowns[[i]], sol[[i]]))
}

# Coefficient of u in a linear expression eq.
.linear_coeff <- function(eq, u) {
  terms <- if (eq$type == "Add") eq$args else list(eq)
  out <- list()
  for (t in terms) {
    if (.has1(t, u)) out[[length(out) + 1L]] <- .div(t, u)
  }
  .add(out)
}

# Gauss-Jordan elimination on an augmented matrix given as rows of
# expressions. Free variables are set to zero. Returns NULL if inconsistent.
.gauss_solve <- function(rows, n) {
  m <- length(rows)
  if (!m) return(rep(list(S_Zero), n))
  A <- rows
  piv_cols <- integer()
  r <- 1L
  for (c in seq_len(n)) {
    if (r > m) break
    pr <- NA_integer_
    for (i in r:m) {
      v <- cancel(A[[i]][[c]])
      A[[i]][[c]] <- v
      if (!is_zero_num(v)) {
        pr <- i
        break
      }
    }
    if (is.na(pr)) next
    if (pr != r) {
      tmp <- A[[r]]; A[[r]] <- A[[pr]]; A[[pr]] <- tmp
    }
    pv <- A[[r]][[c]]
    A[[r]] <- lapply(A[[r]], function(z) cancel(.div(z, pv)))
    for (i in seq_len(m)) {
      if (i == r) next
      fct <- A[[i]][[c]]
      if (is_zero_num(fct)) next
      A[[i]] <- lapply(seq_len(n + 1L), function(j) cancel(.sub(A[[i]][[j]], .mul(list(fct, A[[r]][[j]])))))
    }
    piv_cols <- c(piv_cols, c)
    r <- r + 1L
  }
  if (r <= m) for (i in r:m) if (!is_zero_num(cancel(A[[i]][[n + 1L]]))) return(NULL)
  sol <- rep(list(S_Zero), n)
  for (k in seq_along(piv_cols)) sol[[piv_cols[k]]] <- A[[k]][[n + 1L]]
  attr(sol, "pivots") <- piv_cols
  sol
}
