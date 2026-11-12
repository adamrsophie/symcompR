# Further polynomial functions of sympy/polys: functional decomposition,
# greatest factorial factorization, subresultants, symmetric reduction,
# rational interpolation, minimal polynomials of algebraic numbers,
# Swinnerton-Dyer polynomials, root powers, real root isolation, rational
# roots and the Routh-Hurwitz and Schur stability conditions.

.px <- function(f, x) {
  f <- as_sym(if (inherits(f, "Poly")) as_expr(f) else f)
  if (is.null(x)) {
    fs <- free_symbols(f)
    if (length(fs) != 1L) stop("give the generator explicitly")
    x <- fs[[1L]]
  }
  list(f = expand(f), x = as_sym(x))
}

#' More polynomial functions
#'
#' `sym_decompose(f)` (SymPy's `decompose()`; the name clashes with [stats::decompose()]) writes `f` as a composition `g1(g2(...))` and returns the
#' list of components. `gff_list(f)` is the greatest factorial factorization,
#' a list of `list(factor, k)`. `subresultants(f, g)` returns the subresultant
#' polynomial remainder sequence. `sym_trunc(f, p)` reduces coefficients to
#' symmetric residues modulo `p`. `symmetrize(f, ...)` writes a symmetric
#' polynomial in elementary symmetric polynomials and returns
#' `list(symmetric, remainder)` (with `formal = TRUE` also the symbols `s1`,
#' `s2`, ... and their definitions). `rational_interpolate(data, degnum, x)`
#' fits a rational function with numerator degree `degnum` through the
#' points. `minimal_polynomial(e, x)` (alias `minpoly()`) is the minimal
#' polynomial of an algebraic number built from rationals, `I`, radicals,
#' sums, products and `cos`/`sin` of rational multiples of `pi`.
#' `swinnerton_dyer_poly(n, x)` has the roots `+-sqrt(2) +- ... +- sqrt(p_n)`.
#' `nth_power_roots_poly(f, n)` has the `n`-th powers of the roots of `f`.
#' `intervals(f)` isolates the real roots in rational intervals (refined to
#' width `eps`, within `inf` and `sup`) by continued fractions as SymPy,
#' `refine_root(f, s, t, eps, steps)` narrows one, and `ground_roots(f)` gives the
#' rational roots with multiplicities. `hurwitz_conditions(f, x)` and
#' `schur_conditions(f, x)` are the stability conditions: expressions that
#' must all be positive. `poly_from_expr(e)` returns `list(Poly, gens)`.
#'
#' @param f,g Polynomials.
#' @param x The generator.
#' @param p Modulus.
#' @param ... Generators.
#' @param formal Return symbols for the elementary polynomials.
#' @param data A list of `c(x, y)` points.
#' @param degnum Numerator degree.
#' @param e An algebraic number.
#' @param n Integer.
#' @param s,t Interval endpoints.
#' @param eps Width of the refined interval.
#' @param inf,sup Bounds of the search interval.
#' @param steps Number of refinement steps.
#' @param fast Use fast continued fraction steps.
#' @param sqf Treat `f` as square-free and return bare intervals.
#' @return Expressions, lists or vectors.
#' @examples
#' x <- Symbol("x"); y <- Symbol("y")
#' sym_decompose(x^4 + 2*x^3 - x - 1)
#' subresultants(x^2 + 2*x + 1, x^2 - 1)
#' symmetrize(x^2 + y^2)
#' minimal_polynomial(sqrt(S(2)) + sqrt(S(3)), x)
#' intervals(x^2 - 3)
#' @name polys-extra
NULL

# Dense coefficient lists (highest first) of exact symbols.
.cl <- function(f, x) .ex_trim(.coeff_list(f, x))
.deg <- function(c) length(c) - 1L

.cl_eval_poly <- function(c, h) {
  # Horner with polynomial h, as coefficient lists
  out <- list()
  for (a in c) out <- .ex_sub(.ex_mul(out, h), list(.neg(a)))
  out
}

.right_decompose <- function(f, s) {
  n <- .deg(f)
  lc <- f[[1L]]
  fc <- function(d) if (d < 0 || d > n) S_Zero else f[[n - d + 1L]]
  g <- list(); g[[as.character(s)]] <- S_One
  r <- n %/% s
  for (i in seq_len(s - 1L)) {
    co <- S_Zero
    for (j in 0:(i - 1L)) {
      gj <- g[[as.character(s - j)]]
      if (is.null(gj)) next
      co <- .add(list(co, .mul(list(.int(i - r * j), fc(n + j - i), gj))))
    }
    g[[as.character(s - i)]] <- cancel(.div(co, .mul(list(.int(i * r), lc))))
  }
  lapply(s:0, function(d) { v <- g[[as.character(d)]]; if (is.null(v)) S_Zero else v })
}

.left_decompose <- function(f, h) {
  g <- list()
  while (length(f)) {
    d <- .ex_divmod(f, h)
    if (.deg(d$r) > 0L) return(NULL)
    g <- c(list(if (length(d$r)) d$r[[1L]] else S_Zero), g)
    f <- d$q
  }
  .ex_trim(g)
}

#' @rdname polys-extra
#' @export
sym_decompose <- function(f, x = NULL) {
  v <- .px(f, x); x <- v$x
  fl <- .cl(v$f, x)
  comps <- list()
  repeat {
    df <- .deg(fl)
    found <- FALSE
    for (s in seq_len(max(0L, df - 2L)) + 1L) {
      if (df %% s != 0L) next
      h <- .right_decompose(fl, s)
      g <- .left_decompose(fl, h)
      if (!is.null(g)) {
        comps <- c(list(.ex_to_expr(h, x)), comps)
        fl <- g
        found <- TRUE
        break
      }
    }
    if (!found) break
  }
  c(list(.ex_to_expr(fl, x)), comps)
}

#' @rdname polys-extra
#' @export
gff_list <- function(f, x = NULL) {
  v <- .px(f, x); x <- v$x
  rec <- function(f) {
    f <- expand(.div(f, LC(f, x)))
    if (degree(f, x) <= 0) return(list())
    g <- gcd(f, expand(subs(f, x, .add(list(x, S_One)))))
    H <- rec(g)
    for (i in seq_along(H)) {
      g <- expand(.mul(list(g, subs(H[[i]][[1L]], x, .sub(x, .int(H[[i]][[2L]]))))))
      H[[i]][[2L]] <- H[[i]][[2L]] + 1L
    }
    q <- cancel(.div(f, g))
    if (degree(q, x) <= 0) return(H)
    c(list(list(q, 1L)), H)
  }
  rec(v$f)
}

#' @rdname polys-extra
#' @export
subresultants <- function(f, g, x = NULL) {
  f <- as_sym(f); g <- as_sym(g)
  if (is.null(x)) x <- unique_syms(list(free_symbols(f), free_symbols(g)))[[1L]]
  x <- as_sym(x)
  f <- expand(f); g <- expand(g)
  n <- degree(f, x); m <- degree(g, x)
  if (n < m) { t <- f; f <- g; g <- t; t <- n; n <- m; m <- t }
  if (is_zero_num(f)) return(list())
  if (is_zero_num(g)) return(list(f))
  R <- list(f, g)
  d <- n - m
  b <- .int((-1)^(d + 1))
  h <- expand(.mul(list(prem(f, g, x), b)))
  lc <- LC(g, x)
  cc <- .pow(lc, .int(d))
  cc <- .neg(cc)
  while (!is_zero_num(h)) {
    k <- degree(h, x)
    R[[length(R) + 1L]] <- h
    f <- g; g <- h; d <- m - k; m <- k
    b <- .mul(list(.neg(lc), .pow(cc, .int(d))))
    h <- expand(cancel(.div(prem(f, g, x), b)))
    lc <- LC(g, x)
    cc <- if (d > 1) cancel(.div(.pow(.neg(lc), .int(d)), .pow(cc, .int(d - 1)))) else .neg(lc)
  }
  R
}

#' @rdname polys-extra
#' @export
sym_trunc <- function(f, p, x = NULL) {
  v <- .px(f, x); x <- v$x
  p <- as_sym(p)
  cl <- .coeff_list(v$f, x)
  r <- lapply(cl, function(c) {
    m <- .sub(c, .mul(list(p, fcall("floor", list(.div(c, p))))))
    if (isTRUE(ask_fact(.sub(.mul(list(S_Two, m)), p), "positive"))) .sub(m, p) else m
  })
  expand(.ex_to_expr(r, x))
}

.elementary <- function(vars) {
  n <- length(vars)
  lapply(seq_len(n), function(k) .add(lapply(combinations(seq_len(n), k), function(ix) .mul(vars[ix]))))
}

#' @rdname polys-extra
#' @export
symmetrize <- function(f, ..., formal = FALSE) {
  f <- expand(as_sym(f))
  # sring: the generators are inferred as for Poly (sin(x), not x)
  vars <- .gens_of(list(f), list(...))
  n <- length(vars)
  if (!n) return(if (formal) list(f, S_Zero, list()) else list(f, S_Zero))
  el <- .elementary(vars)
  ss <- lapply(seq_len(n), function(i) Symbol(paste0("s", i)))
  sym_part <- list()
  rest <- f
  weights <- n:1
  # SymPy's smp_symmetrize: among the terms with non-increasing exponents,
  # take the one of largest weighted height (lex order breaks ties)
  repeat {
    if (is_zero_num(rest)) break
    ms <- lapply(do.call(monoms, c(list(rest), vars)), as.integer)
    cs <- do.call(coeffs, c(list(rest), vars))
    o <- do.call(order, c(lapply(seq_len(n), function(i) -vapply(ms, `[`, 0L, i))))
    best <- NULL; height <- -1
    for (i in o) {
      ex <- ms[[i]]
      if (n > 1L && any(diff(ex) > 0)) next
      h <- max(weights * ex)
      if (h > height) { height <- h; best <- i }
    }
    if (is.null(best)) break
    ex <- ms[[best]]; c1 <- cs[[best]]
    e <- c(ex[-n] - ex[-1L], ex[n])
    sym_part[[length(sym_part) + 1L]] <- .mul(c(list(c1), lapply(seq_len(n), function(i) .pow(ss[[i]], .int(e[i])))))
    rest <- expand(.sub(rest, .mul(c(list(c1), lapply(seq_len(n), function(i) .pow(el[[i]], .int(e[i])))))))
  }
  sp <- .add(sym_part)
  if (formal) return(list(sp, rest, lapply(seq_len(n), function(i) list(ss[[i]], el[[i]]))))
  list(xreplace(sp, lapply(seq_len(n), function(i) list(ss[[i]], el[[i]]))), rest)
}

#' @rdname polys-extra
#' @export
rational_interpolate <- function(data, degnum, x = Symbol("x")) {
  xd <- lapply(data, function(p) as_sym(p[[1L]])); yd <- lapply(data, function(p) as_sym(p[[2L]]))
  N <- length(xd)
  k <- N - degnum - 1L
  if (k < 0) stop("too few values for the required degree")
  rows <- lapply(seq_len(N), function(i) {
    c(lapply(0:degnum, function(j) .pow(xd[[i]], .int(j))),
      lapply(0:k, function(j) .neg(.mul(list(yd[[i]], .pow(xd[[i]], .int(j)))))))
  })
  ns <- nullspace(Matrix(rows))
  if (!length(ns)) stop("no rational interpolant")
  r <- ns[[1L]]$data
  num <- .add(lapply(0:degnum, function(j) .mul(list(r[[j + 1L]], .pow(x, .int(j))))))
  den <- .add(lapply(0:k, function(j) .mul(list(r[[degnum + 2L + j]], .pow(x, .int(j))))))
  cancel(.div(num, den))
}

# ---------------------------------------------------------------------------
# Minimal polynomials of algebraic numbers

.annihilator <- function(e, x) {
  y <- Dummy("y")
  if (is_exact_number(e)) return(.sub(x, e))
  if (is_const(e, "I")) return(.add(list(.sq(x), S_One)))
  t <- e$type
  if (t == "Add" || t == "Mul") {
    a <- .annihilator(e$args[[1L]], y)
    rest <- if (length(e$args) == 2L) e$args[[2L]] else new_sym(t, e$args[-1L])
    b <- .annihilator(rest, x)
    if (t == "Add") return(expand(resultant(a, expand(subs(b, x, .sub(x, y))), y)))
    db <- degree(b, x)
    hb <- expand(.mul(list(.pow(y, .int(db)), subs(b, x, .div(x, y)))))
    return(expand(resultant(a, cancel(hb), y)))
  }
  if (t == "Pow" && e$args[[2L]]$type %in% c("Integer", "Rational")) {
    q <- e$args[[2L]]$value
    p <- gmp::numerator(q); d <- gmp::denominator(q)
    a <- .annihilator(e$args[[1L]], y)
    if (p < 0) {
      a <- expand(.mul(list(.pow(y, .int(degree(a, y))), subs(a, y, .inv(y)))))
      p <- -p
    }
    return(expand(resultant(a, .sub(.pow(x, .int(d)), .pow(y, .int(p))), y)))
  }
  if (t == "Function" && e$name %in% c("cos", "sin")) {
    r <- cancel(.div(e$args[[1L]], S_Pi))
    if (!is_exact_number(r)) stop("not an algebraic number: ", .str(e))
    if (e$name == "sin") r <- .sub(S_Half, r)
    p <- gmp::numerator(r$value); q <- gmp::denominator(r$value)
    # T_q(cos(p*pi/q)) = cos(p*pi) = (-1)^p
    return(expand(.sub(chebyshevt(as.numeric(q), x), .int((-1)^as.numeric(p %% 2)))))
  }
  stop("not an algebraic number: ", .str(e))
}

#' @rdname polys-extra
#' @export
minimal_polynomial <- function(e, x = Symbol("x")) {
  if (is_sym(e) && e$type == "AlgebraicNumber") return(minpoly_of_element(e, x))
  e <- as_sym(e); x <- as_sym(x)
  A <- .annihilator(e, x)
  v <- as.complex(.neval(e, list()))
  fl <- factor_list(A, x)$factors
  best <- NULL; bd <- Inf
  for (fm in fl) {
    fac <- fm[[1L]]
    if (degree(fac, x) < 1) next
    d <- tryCatch(base::Mod(as.complex(.neval(fac, stats::setNames(list(v), x$name)))), error = function(err) Inf)
    sc <- tryCatch(sum(base::Mod(vapply(.coeff_list(fac, x), function(c) as.complex(.neval(c, list())), 0i))) *
                     max(1, base::Mod(v))^degree(fac, x), error = function(err) 1)
    if (d / sc < bd) { bd <- d / sc; best <- fac }
  }
  if (is.null(best)) stop("no minimal polynomial found")
  # primitive integer polynomial with a positive leading coefficient
  cl <- .coeff_list(best, x)
  den <- Reduce(function(a, b) a * b / gmp::gcd.bigz(a, b), lapply(cl, function(c) gmp::denominator(c$value)))
  cl <- lapply(cl, function(c) c$value * den)
  g <- Reduce(gmp::gcd.bigz, lapply(cl, function(c) abs(gmp::numerator(c))))
  sg <- if (cl[[1L]] < 0) -1 else 1
  expand(.ex_to_expr(lapply(cl, function(c) make_rational(c / g * sg)), x))
}

#' @rdname polys-extra
#' @export
minpoly <- minimal_polynomial

#' @rdname polys-extra
#' @export
swinnerton_dyer_poly <- function(n, x = Symbol("x")) {
  x <- as_sym(x); n <- as.integer(n)
  if (n < 1L) stop("n must be positive")
  ps <- vapply(seq_len(n), function(i) as.numeric(.prime_nth(i)), 0)
  y <- Dummy("y")
  P <- .sub(.sq(x), .int(ps[1L]))
  for (p in ps[-1L]) P <- expand(resultant(expand(subs(P, x, .sub(x, y))), .sub(.sq(y), .int(p)), y))
  P
}

#' @rdname polys-extra
#' @export
nth_power_roots_poly <- function(f, n, x = NULL) {
  v <- .px(f, x); x <- v$x
  y <- Dummy("y")
  r <- expand(resultant(subs(v$f, x, y), .sub(x, .pow(y, .int(n))), y))
  lt <- .pow(LC(v$f, x), .int(n))
  expand(.mul(list(.div(lt, LC(r, x)), r)))
}

#' @rdname polys-extra
#' @export
ground_roots <- function(f, x = NULL) {
  v <- .px(f, x); x <- v$x
  rs <- roots(v$f, x)
  m <- attr(rs, "multiplicity")
  keep <- vapply(rs, is_exact_number, TRUE)
  out <- stats::setNames(m[keep], vapply(rs[keep], .str, ""))
  out[order(vapply(rs[keep], function(r) as.numeric(r$value), 0))]
}

#' @rdname polys-extra
#' @export
intervals <- function(f, x = NULL, eps = NULL, inf = NULL, sup = NULL, fast = FALSE, sqf = FALSE) {
  v <- .px(f, x); x <- v$x
  q <- function(z) if (is.null(z)) NULL else gmp::as.bigq(as_sym(z)$value)
  eps <- if (is.null(eps)) NULL else .q_exact(eps)
  rs <- .zp_isolate_real_roots(.zp_from_expr(v$f, x), eps = eps, inf = q(inf), sup = q(sup), fast = fast, sqf = sqf)
  lapply(rs, function(r) {
    iv <- list(make_rational(r[[1L]]), make_rational(r[[2L]]))
    if (sqf) iv else list(iv, r[[3L]])
  })
}

#' @rdname polys-extra
#' @export
refine_root <- function(f, s, t, eps = NULL, steps = NULL, fast = FALSE, x = NULL) {
  v <- .px(f, x); x <- v$x
  # Poly.refine_root: one step when neither eps nor steps is given
  if (is.null(eps) && is.null(steps)) steps <- 1L
  r <- .zp_refine(.zp_from_expr(v$f, x), .q_exact(s), .q_exact(t),
                  eps = if (is.null(eps)) NULL else .q_exact(eps), steps = steps, fast = fast)
  list(make_rational(r[[1L]]), make_rational(r[[2L]]))
}

# An exact rational from a number, a string or a Float (QQ.convert).
.q_exact <- function(z) {
  z <- as_sym(z)
  if (z$type == "Float") return(.mp_to_q(.float_mp(z)))
  gmp::as.bigq(z$value)
}

# ---------------------------------------------------------------------------
# Stability conditions (sympy/polys/rootconditions.py)

.hurwitz_ff <- function(p) {
  isneg <- function(c) if (is_number_type(c)) num_sign(c) < 0 else NA
  if (length(p) == 1L) return(list())
  if (length(p) == 2L) return(list(.mul(list(p[[1L]], p[[2L]]))))
  if (length(p) == 3L) return(list(expand(.mul(list(p[[1L]], p[[2L]]))), expand(.mul(list(p[[1L]], p[[3L]])))))
  LCo <- p[[1L]]; TC <- p[[length(p)]]
  lcneg <- isneg(LCo)
  if (!is.na(lcneg)) for (co in list(TC, p[[2L]])) {
    ng <- isneg(co)
    if (!is.na(ng) && (is_zero_num(co) || xor(lcneg, ng))) return(list(S_NegativeOne))
  }
  p1s <- list(p[[2L]])
  while (length(p) > 3L) {
    qs <- lapply(p[-1L], function(qi) expand(.mul(list(p[[2L]], qi))))
    for (i in seq(2L, length(qs) - 1L, by = 2L)) if (i <= length(qs) - 1L) qs[[i]] <- expand(.sub(qs[[i]], .mul(list(p[[i + 2L]], p[[1L]]))))
    p <- qs
    p1 <- p[[2L]]
    if (is_zero_num(p1)) return(list(S_NegativeOne))
    np <- length(p1s)
    if (np >= 2L) p1 <- cancel(.div(p1, p1s[[np - 1L]]))
    if (np >= 3L) {
      p1 <- cancel(.div(p1, p1s[[np - 2L]])); p1 <- cancel(.div(p1, p1s[[np - 2L]]))
      p <- lapply(p, function(z) cancel(.div(.div(z, p1s[[np - 2L]]), p1s[[np - 2L]])))
    }
    ng <- isneg(p1)
    if (!is.na(ng)) {
      if (np %% 2L == 0L) { if (!is.na(lcneg) && xor(lcneg, ng)) return(list(S_NegativeOne)) }
      else if (ng) return(list(S_NegativeOne))
    }
    p1s[[length(p1s) + 1L]] <- p1
  }
  if (!is_one_num(LCo)) for (i in seq(1L, length(p1s), by = 2L)) p1s[[i]] <- expand(.mul(list(p1s[[i]], LCo)))
  p1s[[length(p1s) + 1L]] <- expand(.mul(list(LCo, TC)))
  lapply(p1s, expand)
}

.cleared_coeffs <- function(f, x) {
  cl <- .coeff_list(f, x)
  den <- S_One
  for (c in cl) { d <- fraction(together(c))[[2L]]; den <- sym_lcm(den, d) }
  lapply(cl, function(c) cancel(.mul(list(c, den))))
}

#' @rdname polys-extra
#' @export
hurwitz_conditions <- function(f, x = NULL) {
  v <- .px(f, x)
  if (is_zero_num(v$f)) stop("zero polynomial")
  .hurwitz_ff(.cleared_coeffs(v$f, v$x))
}

#' @rdname polys-extra
#' @export
schur_conditions <- function(f, x = NULL) {
  v <- .px(f, x); x <- v$x
  if (is_zero_num(v$f)) stop("zero polynomial")
  if (is_zero_num(simplify(subs(v$f, x, S_NegativeOne)))) return(list(S_NegativeOne))
  n <- degree(v$f, x)
  ft <- expand(cancel(.mul(list(.pow(.sub(S_One, x), .int(n)), subs(v$f, x, .div(.add(list(x, S_One)), .sub(S_One, x)))))))
  .hurwitz_ff(.cleared_coeffs(ft, x))
}

#' @rdname polys-extra
#' @export
poly_from_expr <- function(e, ...) {
  P <- Poly(e, ...)
  list(P, list(gens = P$gens))
}

symcompR_count_roots <- function(f, x, a, b) count_roots(f, x, a, b)
