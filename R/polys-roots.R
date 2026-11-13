# Roots of univariate polynomials, after sympy/polys/polyroots.py and
# sympy/polys/rootoftools.py.

.roots_quadratic <- function(cf) {
  a <- cf[[1L]]; b <- cf[[2L]]; c <- cf[[3L]]
  if (is_zero_num(c)) return(list(S_Zero, .neg(cancel(.div(b, a)))))
  if (is_zero_num(b)) {
    r <- .sqrt(cancel(.div(.neg(c), a)))
    return(list(.neg(r), r))
  }
  d <- expand(.sub(.sq(b), .mul(list(.int(4), a, c))))
  A <- .mul(list(S_Two, a))
  B <- cancel(.div(.neg(b), A))
  D <- .div(.sqrt(d), A)
  if (is_zero_num(d)) return(list(B, B))
  list(expand(.sub(B, D)), expand(.add(list(B, D))))
}

.roots_cubic <- function(cf) {
  a0 <- cf[[1L]]
  cf <- lapply(cf, function(z) cancel(.div(z, a0)))
  a <- cf[[2L]]; b <- cf[[3L]]; c <- cf[[4L]]
  if (is_zero_num(c)) {
    rr <- .roots_quadratic(list(S_One, a, b))
    return(list(rr[[1L]], S_Zero, rr[[2L]]))
  }
  p <- cancel(.sub(b, .div(.sq(a), .int(3))))
  q <- cancel(.add(list(c, .neg(.div(.mul(list(a, b)), .int(3))), .mul(list(.rat(2, 27), .pow(a, .int(3)))))))
  pon3 <- .div(p, .int(3))
  aon3 <- .div(a, .int(3))
  u1 <- NULL
  if (is_zero_num(p)) {
    if (is_zero_num(q)) return(rep(list(.neg(aon3)), 3))
    u1 <- if (isTRUE(ask_fact(q, "positive"))) .neg(real_root(q, 3)) else real_root(.neg(q), 3)
  } else if (is_zero_num(q)) {
    rr <- .roots_quadratic(list(S_One, S_Zero, p))
    return(lapply(list(rr[[1L]], S_Zero, rr[[2L]]), function(t) .sub(t, aon3)))
  } else if (isTRUE(ask_fact(q, "real")) && isTRUE(ask_fact(q, "negative"))) {
    u1 <- .neg(real_root(.add(list(.div(.neg(q), S_Two), .sqrt(.add(list(.div(.sq(q), .int(4)), .pow(pon3, .int(3))))))), 3))
  }
  coeff <- .mul(list(S_I, .sqrt(.int(3)), S_Half))
  if (is.null(u1)) {
    us <- list(S_One, .add(list(.rat(-1, 2), coeff)), .sub(.rat(-1, 2), coeff))
    D0 <- .sub(.sq(a), .mul(list(.int(3), b)))
    D1 <- .add(list(.mul(list(S_Two, .pow(a, .int(3)))), .mul(list(.int(-9), a, b)), .mul(list(.int(27), c))))
    C <- .pow(.div(.add(list(D1, .sqrt(.sub(.sq(D1), .mul(list(.int(4), .pow(D0, .int(3)))))))), S_Two), .rat(1, 3))
    return(lapply(us, function(uk) expand(.neg(.div(.add(list(a, .mul(list(uk, C)), .div(D0, .mul(list(C, uk))))), .int(3))))))
  }
  u2 <- .mul(list(u1, .add(list(.rat(-1, 2), coeff))))
  u3 <- .mul(list(u1, .sub(.rat(-1, 2), coeff)))
  if (is_zero_num(p)) return(lapply(list(u1, u2, u3), function(u) .sub(u, aon3)))
  lapply(list(u1, u2, u3), function(u) expand(.add(list(.neg(u), .div(pon3, u), .neg(aon3)))))
}

.roots_quartic <- function(cf) {
  a0 <- cf[[1L]]
  cf <- lapply(cf, function(z) cancel(.div(z, a0)))
  a <- cf[[2L]]; b <- cf[[3L]]; c <- cf[[4L]]; d <- cf[[5L]]
  # depressed quartic y^4 + p y^2 + q y + r with x = y - a/4
  p <- cancel(.sub(b, .mul(list(.rat(3, 8), .sq(a)))))
  q <- cancel(.add(list(c, .neg(.mul(list(S_Half, a, b))), .mul(list(.rat(1, 8), .pow(a, .int(3)))))))
  r <- cancel(.add(list(d, .neg(.mul(list(.rat(1, 4), a, c))), .mul(list(.rat(1, 16), .sq(a), b)),
                        .mul(list(.rat(-3, 256), .pow(a, .int(4)))))))
  shift <- .neg(.div(a, .int(4)))
  if (is_zero_num(q)) {
    zs <- .roots_quadratic(list(S_One, p, r))
    out <- list()
    for (z in zs) {
      s <- .sqrt(z)
      out <- c(out, list(.add(list(.neg(s), shift)), .add(list(s, shift))))
    }
    return(lapply(out, expand))
  }
  # resolvent cubic 8m^3 + 8p m^2 + (2p^2 - 8r) m - q^2 = 0
  res <- list(.int(8), .mul(list(.int(8), p)), .sub(.mul(list(S_Two, .sq(p))), .mul(list(.int(8), r))), .neg(.sq(q)))
  ms <- .roots_of_coeffs(res)
  ms <- Filter(function(m) !is_zero_num(m), ms)
  rat <- Filter(is_exact_number, ms)
  m <- if (length(rat)) rat[[1L]] else ms[[1L]]
  s1 <- .sqrt(.mul(list(S_Two, m)))
  base <- .add(list(.div(p, S_Two), m))
  t <- .div(q, .mul(list(S_Two, s1)))
  q1 <- .roots_quadratic(list(S_One, .neg(s1), .add(list(base, t))))
  q2 <- .roots_quadratic(list(S_One, s1, .sub(base, t)))
  lapply(c(q1, q2), function(y) expand(.add(list(y, shift))))
}

.root_of_unity <- function(k, n) {
  ang <- .mul(list(.rat(2 * k, n), S_Pi))
  c1 <- fcall("cos", list(ang))
  s1 <- fcall("sin", list(ang))
  if (c1$type == "Function" || s1$type == "Function") return(fcall("exp", list(.mul(list(S_I, ang)))))
  expand(.add(list(c1, .mul(list(S_I, s1)))))
}

.roots_binomial <- function(cf) {
  n <- length(cf) - 1L
  a <- cf[[1L]]
  b <- cf[[length(cf)]]
  alpha <- cancel(.div(.neg(b), a))
  out <- list()
  if (is_number_type(alpha) && num_sign(alpha) < 0) {
    mag <- .pow(.neg(alpha), .rat(1, n))
    for (k in 0:(n - 1L)) {
      ang <- .mul(list(.rat(2 * k + 1, n), S_Pi))
      c1 <- fcall("cos", list(ang)); s1 <- fcall("sin", list(ang))
      z <- if (c1$type == "Function" || s1$type == "Function") fcall("exp", list(.mul(list(S_I, ang))))
           else .add(list(c1, .mul(list(S_I, s1))))
      out[[length(out) + 1L]] <- expand(.mul(list(mag, z)))
    }
  } else {
    base <- .pow(alpha, .rat(1, n))
    for (k in 0:(n - 1L)) out[[length(out) + 1L]] <- expand(.mul(list(base, .root_of_unity(k, n))))
  }
  out
}

# Roots of an irreducible (or unfactorable) polynomial given by coefficients.
.roots_of_coeffs <- function(cf) {
  while (length(cf) && is_zero_num(cf[[1L]])) cf <- cf[-1L]
  n <- length(cf) - 1L
  if (n <= 0L) return(list())
  if (n == 1L) return(list(cancel(.div(.neg(cf[[2L]]), cf[[1L]]))))
  inner_zero <- all(vapply(cf[-c(1L, length(cf))], is_zero_num, TRUE))
  if (inner_zero && n > 2L) return(.roots_binomial(cf))
  if (n == 2L) return(.roots_quadratic(cf))
  if (n == 3L) return(.roots_cubic(cf))
  if (n == 4L) return(.roots_quartic(cf))
  # x^n structure in a smaller power: f(x) = g(x^k)
  NULL
}

.sort_roots <- function(rs) {
  if (length(rs) < 2L) return(rs)
  vals <- lapply(rs, function(r) tryCatch(as.complex(.neval(r, list())), error = function(e) NA))
  re <- vapply(vals, function(v) if (is.na(v)) Inf else Re(v), 0)
  im <- vapply(vals, function(v) if (is.na(v)) Inf else Im(v), 0)
  isreal <- abs(im) < 1e-12 * pmax(1, abs(re))
  rs[order(!isreal, re, im)]
}

#' Roots of a polynomial
#'
#' `roots(f, x)` returns the roots of a univariate polynomial, solved by
#' factorization over the rationals and by radical formulas up to degree
#' four (and for binomials), like SymPy's `roots()`. The result is a list of
#' distinct roots with a `multiplicity` attribute; with `multiple = TRUE`
#' each root is repeated by its multiplicity. Factors that are not solvable
#' by these formulas are left out unless `rootof = TRUE`, in which case they
#' appear as [CRootOf()] objects.
#'
#' @param f A polynomial expression or [Poly()].
#' @param x The variable; when `k` is missing, `x` is the index and the
#'   variable is the single free symbol of `f`.
#' @param multiple Repeat roots by multiplicity.
#' @param rootof Include `CRootOf` roots for unsolvable factors.
#' @param filter Keep only `"real"`, `"complex"`, `"rational"`, or all roots
#'   (`NULL`).
#' @param cubics,quartics,quintics Use the formulas for degrees three to five.
#' @param trig Use the trigonometric form for cubics with three real roots.
#' @param predicate A function selecting roots.
#' @param strict Signal an error when some roots are not found.
#' @param auto Convert ring coefficients to their field.
#' @return A list of expressions.
#' @examples
#' x <- Symbol("x")
#' roots(x^2 - 3*x + 2, x)
#' roots(x^3 - 1, x)
#' roots((x - 1)^2 * (x + 3), x)
#' @export
roots <- function(f, x = NULL, multiple = FALSE, rootof = FALSE, filter = NULL, cubics = TRUE, trig = FALSE,
                  quartics = TRUE, quintics = FALSE, predicate = NULL, strict = FALSE, auto = TRUE) {
  if (inherits(f, "CartanType")) return(.cartan_num_roots(f))
  if (!rootof) return(.roots_sympy(f, x, auto = auto, cubics = cubics, trig = trig, quartics = quartics, quintics = quintics,
                                   multiple = multiple, filter = filter, predicate = predicate, strict = strict))
  .roots_legacy(f, x, multiple, rootof, filter)
}

.roots_legacy <- function(f, x = NULL, multiple = FALSE, rootof = FALSE, filter = NULL) {
  if (inherits(f, "Poly")) {
    if (is.null(x)) x <- f$gens[[1L]]
    f <- as_expr(f)
  }
  f <- as_sym(f)
  if (is.null(x)) {
    fs <- free_symbols(f)
    if (length(fs) != 1L) stop("roots() needs the variable")
    x <- fs[[1L]]
  }
  x <- as_sym(x)
  cf <- .coeff_list(f, x)
  vals <- list()
  mult <- integer()
  add_root <- function(r, m) {
    k <- match(r$key, vapply(vals, .k, ""))
    if (is.na(k)) {
      vals[[length(vals) + 1L]] <<- r
      mult[length(mult) + 1L] <<- as.integer(m)
    } else {
      mult[k] <<- mult[k] + as.integer(m)
    }
  }
  if (!length(cf) || length(cf) == 1L) return(structure(list(), multiplicity = integer()))
  numeric_cf <- all(vapply(cf, is_exact_number, TRUE))
  if (numeric_cf) {
    q <- do.call(c, lapply(cf, `[[`, "value"))
    den <- .lcm_den(q)
    zl <- .zp_factor_list(gmp::numerator(q * gmp::as.bigq(den)))
    for (pr in zl$factors) {
      g <- lapply(as.list(pr[[1L]]), function(v) make_rational(gmp::as.bigq(v)))
      rs <- .roots_of_coeffs(g)
      if (is.null(rs)) rs <- .roots_by_decompose(g, rootof)
      for (r in rs) add_root(r, pr[[2L]])
    }
  } else {
    fl <- tryCatch(factor_list(f), error = function(e) NULL)
    facs <- if (is.null(fl)) list(list(f, 1L)) else fl$factors
    for (pr in facs) {
      if (!.has1(pr[[1L]], x)) next
      g <- .coeff_list(pr[[1L]], x)
      rs <- .roots_of_coeffs(g)
      if (is.null(rs)) rs <- list()
      for (r in rs) add_root(r, pr[[2L]])
    }
  }
  ord <- if (length(vals) > 1L) {
    o <- .sort_roots(vals)
    match(vapply(o, .k, ""), vapply(vals, .k, ""))
  } else seq_along(vals)
  vals <- vals[ord]
  mult <- mult[ord]
  if (!is.null(filter)) {
    # roots(filter=...): SymPy's 'Z', 'Q', 'R', 'I' and 'C' handlers, and the
    # long names of this port
    keep <- vapply(vals, function(r) switch(filter,
                                            Z = r$type == "Integer",
                                            Q = r$type %in% c("Integer", "Rational"),
                                            R = all(vapply(.as_numer_denom(r), function(a) isTRUE(ask_fact(a, "real")), NA)),
                                            I = isTRUE(ask_fact(r, "imaginary")),
                                            real = isTRUE(ask_fact(r, "real")),
                                            rational = isTRUE(ask_fact(r, "rational")),
                                            complex = , C = TRUE,
                                            stop("Invalid filter: ", filter, call. = FALSE)), TRUE)
    vals <- vals[keep]
    mult <- mult[keep]
  }
  if (multiple) return(unlist(Map(function(r, m) rep(list(r), m), vals, mult), recursive = FALSE))
  structure(vals, multiplicity = mult)
}

# Roots of an irreducible integer polynomial of degree >= 5: try
# f(x) = g(x^k), otherwise CRootOf placeholders.
.roots_by_decompose <- function(g, rootof) {
  n <- length(g) - 1L
  nz <- which(!vapply(g, is_zero_num, TRUE))
  pows <- n - (nz - 1L)
  k <- Reduce(function(a, b) { while (b) { t <- a %% b; a <- b; b <- t }; a }, pows)
  if (k > 1L) {
    h <- g[seq(1L, length(g), by = k)]
    inner <- .roots_of_coeffs(h)
    if (!is.null(inner)) {
      out <- list()
      for (r in inner) out <- c(out, .roots_binomial(c(list(S_One), rep(list(S_Zero), k - 1L), list(.neg(r)))))
      return(out)
    }
  }
  if (!rootof) return(list())
  xx <- Dummy("x")
  pe <- .coeffs_to_expr(g, xx)
  lapply(seq_len(n) - 1L, function(i) CRootOf(pe, xx, i))
}

#' Indexed complex roots of polynomials
#'
#' `CRootOf(f, x, k)` is the `k`-th root (0 based) of the polynomial `f` in
#' `x`, ordered with real roots first in increasing order and then complex
#' roots by real and imaginary part, as in SymPy. It evaluates numerically
#' with [N()].
#'
#' @param f A polynomial.
#' @param x The variable; when `k` is missing, `x` is the index and the
#'   variable is the single free symbol of `f`.
#' @param k Index of the root, from 0.
#' @param radicals Accepted for SymPy compatibility; the root is always kept
#'   as a CRootOf.
#' @return A `sym` CRootOf object.
#' @export
CRootOf <- function(f, x, k, radicals = TRUE) {
  f <- as_sym(f)
  # CRootOf(f, k) for a univariate f, as in SymPy
  if (missing(k)) {
    k <- x
    fs <- free_symbols(f)
    if (length(fs) != 1L) stop("CRootOf(f, k) needs a univariate polynomial", call. = FALSE)
    x <- fs[[1L]]
  }
  x <- as_sym(x)
  if (is_sym(k)) k <- as.integer(as.character(k$value))
  xx <- Symbol("x")
  pe <- expand(subs(f, x, xx))
  new_sym("CRootOf", list(pe, .int(k)))
}

register_str_printer("CRootOf", function(e) paste0("CRootOf(", .str(e$args[[1L]]), ", ", .str(e$args[[2L]]), ")"))
register_bound_vars("CRootOf", function(e) "S1:x")
register_neval("CRootOf", function(e, env) {
  rs <- .numeric_roots(e$args[[1L]], Symbol("x"))
  rs[[as.integer(as.character(e$args[[2L]]$value)) + 1L]]
})
.facts_hooks$CRootOf <- function(e) {
  v <- .neval(e, list())
  out <- c(commutative = TRUE, algebraic = TRUE, finite = TRUE)
  if (abs(Im(v)) < 1e-12) out["real"] <- TRUE else out["real"] <- FALSE
  out
}

.numeric_roots <- function(f, x, sqf = TRUE) {
  if (sqf) {
    # repeated roots: solve each square-free factor, then repeat
    sl <- tryCatch(sqf_list(f, x), error = function(e) NULL)
    if (!is.null(sl) && length(sl$factors) && (length(sl$factors) > 1L || sl$factors[[1L]][[2L]] > 1)) {
      rs <- unlist(lapply(sl$factors, function(fm) rep(.numeric_roots(fm[[1L]], x, FALSE), fm[[2L]])), recursive = FALSE)
      v <- vapply(rs, function(z) as.complex(z), 0i)
      o <- order(Im(v) != 0, Re(v), abs(Im(v)), Im(v))
      return(rs[o])
    }
  }
  cf <- .coeff_list(f, x)
  v <- vapply(cf, function(z) as.complex(.neval(z, list())), 0i)
  if (all(Im(v) == 0)) v <- Re(v)
  r <- polyroot(rev(v))
  # Newton polishing
  dv <- if (length(v) > 1L) v[-length(v)] * ((length(v) - 1L):1) else 0
  ev <- function(c, z) { s <- 0; for (a in c) s <- s * z + a; s }
  for (it in 1:5) {
    d <- ev(dv, r)
    # leave repeated roots (zero derivative) as polyroot() found them
    st <- ifelse(base::Mod(d) > 0, ev(v, r) / d, 0)
    st[!is.finite(base::Mod(st))] <- 0
    r <- r - st
  }
  re <- Re(r)
  im <- Im(r)
  # polishing noise only: a tiny imaginary part relative to the root itself
  im[abs(im) < 1e-10 * abs(re)] <- 0
  re[abs(re) < 1e-14 * abs(im)] <- 0
  if (is.double(v)) {
    near <- im != 0 & abs(im) < 1e-6 * pmax(1, abs(re))
    for (i in which(near)) if (abs(ev(v, re[i])) <= 1e-9 * ev(abs(v), abs(re[i]))) im[i] <- 0
  }
  r <- complex(real = re, imaginary = im)
  real <- im == 0
  # rounded keys: conjugates whose parts differ only by polishing noise sort
  # with the negative imaginary part first, as SymPy orders CRootOf
  o <- order(!real, signif(re, 12), signif(abs(im), 12), im)
  lapply(r[o], function(z) if (Im(z) == 0) Re(z) else z)
}

#' Numerical and real roots
#'
#' `nroots(f, x)` returns all roots as floating point numbers (via
#' [polyroot()] with Newton polishing). `real_roots(f, x)` returns the real
#' roots with multiplicity, exact where possible. `count_roots(f, x, a, b)`
#' counts real roots in `[a, b]` with a Sturm sequence. `sturm()` returns the
#' Sturm sequence.
#'
#' @param f A polynomial.
#' @param x The variable; when `k` is missing, `x` is the index and the
#'   variable is the single free symbol of `f`.
#' @param a,b Interval ends (default the whole line).
#' @return A list of expressions, or an integer.
#' @examples
#' x <- Symbol("x")
#' nroots(x^5 - x + 1, x)
#' count_roots(x^3 - 2*x, x)
#' @export
nroots <- function(f, x = NULL) {
  f <- as_sym(f)
  if (is.null(x)) x <- free_symbols(f)[[1L]]
  lapply(.numeric_roots(f, as_sym(x)), function(v) .round_sym(v, 15))
}

#' @rdname nroots
#' @export
real_roots <- function(f, x = NULL) {
  f <- as_sym(f)
  if (is.null(x)) x <- free_symbols(f)[[1L]]
  x <- as_sym(x)
  r <- roots(f, x, rootof = TRUE)
  m <- attr(r, "multiplicity")
  keep <- vapply(r, function(z) isTRUE(ask_fact(z, "real")), TRUE)
  out <- unlist(Map(function(z, k) rep(list(z), k), r[keep], m[keep]), recursive = FALSE)
  if (is.null(out)) list() else out
}

#' @rdname nroots
#' @export
sturm <- function(f, x = NULL) {
  f <- as_sym(f)
  if (is.null(x)) x <- free_symbols(f)[[1L]]
  x <- as_sym(x)
  p <- up_from_qpoly(expr_to_qpoly(f, list(x)), 1L)
  seqs <- .sturm_seq(p)
  lapply(seqs, up_to_expr, x = x)
}

.sturm_seq <- function(p) {
  out <- list(p, up_diff(p))
  repeat {
    n <- length(out)
    if (!length(out[[n]])) {
      out <- out[-n]
      break
    }
    r <- up_divmod(out[[n - 1L]], out[[n]])$r
    if (!length(r)) break
    out[[n + 1L]] <- -r
  }
  out
}

#' @rdname nroots
#' @export
count_roots <- function(f, x = NULL, a = -Inf, b = Inf) {
  f <- as_sym(f)
  if (is.null(x)) x <- free_symbols(f)[[1L]]
  x <- as_sym(x)
  p <- up_from_qpoly(expr_to_qpoly(f, list(x)), 1L)
  sq <- up_sqf_list(p)
  sf <- gmp::as.bigq(1)
  for (pr in sq) sf <- up_mul(sf, pr[[1L]])
  seqs <- .sturm_seq(sf)
  changes <- function(at) {
    s <- vapply(seqs, function(q) {
      if (is.infinite(at)) {
        d <- up_deg(q)
        return((if (q[1L] > 0) 1 else -1) * if (at < 0 && d %% 2L == 1L) -1 else 1)
      }
      v <- up_eval(q, gmp::as.bigq(at))
      if (v > 0) 1 else if (v < 0) -1 else 0
    }, 0)
    s <- s[s != 0]
    sum(diff(s) != 0)
  }
  av <- if (is.infinite(a)) a else as_bigq(as_sym(a))
  bv <- if (is.infinite(b)) b else as_bigq(as_sym(b))
  n <- changes(av) - changes(bv)
  if (!is.infinite(a) && up_eval(sf, av) == 0) n <- n + 1L
  as.integer(n)
}
