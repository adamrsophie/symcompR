# User-facing polynomial functions, after sympy/polys/polytools.py.

.as_gens <- function(gens) {
  if (!length(gens)) return(NULL)
  out <- list()
  for (g in gens) {
    if (is.list(g) && !inherits(g, "sym")) out <- c(out, lapply(g, as_sym))
    else out[[length(out) + 1L]] <- as_sym(g)
  }
  out
}

.gens_of <- function(exprs, gens) {
  gens <- .as_gens(gens)
  if (!is.null(gens)) return(gens)
  .find_gens(lapply(exprs, function(z) expand(as_sym(z))))
}

.domain_of <- function(p) {
  if (qp_is_zero(p)) return("ZZ")
  if (all(gmp::denominator(p$c) == 1)) "ZZ" else "QQ"
}

#' Polynomial objects
#'
#' `Poly(expr, x, y, ...)` represents `expr` as a polynomial in the given
#' generators (by default, every symbol and non polynomial subexpression).
#' Coefficients are exact rationals; symbols that are not generators become
#' part of the coefficients when generators are given.
#'
#' @param expr An expression.
#' @param ... Generators.
#' @param p A `Poly`.
#' @return A `Poly` object.
#' @examples
#' x <- Symbol("x"); y <- Symbol("y")
#' p <- Poly(x^2 + 2 * x * y + 1, x, y)
#' p
#' degree(p, x)
#' @export
Poly <- function(expr, ...) {
  if (inherits(expr, "Poly")) return(expr)
  expr <- as_sym(expr)
  pa <- .split_poly_opts(list(...))
  if (length(pa$opts)) {
    if (!is.null(P <- .poly_with_domain(expr, pa))) return(P)
    Fd <- .poly_field(pa$opts, list(expr))
    P <- do.call(Poly, c(list(expr), pa$gens))
    if (Fd$kind != "GF") return(P)
    # coefficients in GF(p), stored in the symmetric representation
    p <- Fd$p
    cs <- vapply(seq_along(P$rep$c), function(i) {
      q <- P$rep$c[i]
      v <- (as.numeric(gmp::mod.bigz(gmp::numerator(q), p)) * .gf_inv(as.numeric(gmp::mod.bigz(gmp::denominator(q), p)), p)) %% p
      .gf_int(v, p, Fd$symmetric)
    }, 0)
    P$rep <- qp_norm(P$rep$e, gmp::as.bigq(cs))
    P$modulus <- p
    return(P)
  }
  gens <- .gens_of(list(expr), list(...))
  rep <- expr_to_qpoly(expr, gens, strict = TRUE)
  if (is.null(rep)) {
    return(.poly_symbolic(expr, gens))
  }
  structure(list(rep = rep, gens = gens, coeff_expr = NULL), class = "Poly")
}

# Poly whose coefficients involve symbols other than the generators: stored
# over all symbols, with the generators listed first.
.poly_symbolic <- function(expr, gens) {
  all_gens <- .find_gens(list(expand(expr)))
  gk <- vapply(gens, .k, "")
  extra <- Filter(function(g) !(g$key %in% gk), all_gens)
  full <- c(gens, extra)
  rep <- expr_to_qpoly(expr, full, strict = FALSE)
  .poly_frac_tag(structure(list(rep = rep, gens = full, ngens = length(gens)), class = "Poly"))
}

.poly_ngens <- function(p) if (is.null(p$ngens)) length(p$gens) else p$ngens

#' @rdname Poly
#' @export
as_expr <- function(p, ...) {
  if (inherits(p, "Poly")) return(qpoly_to_expr(p$rep, p$gens))
  if (is_sym(p) && p$type == "AlgebraicNumber") return(.an_value(p))
  if (inherits(p, "sym_monomial")) return(.monomial_expr(p, list(...)))
  as_sym(p)
}

#' @export
print.Poly <- function(x, ...) {
  ng <- .poly_ngens(x)
  gens <- vapply(x$gens[seq_len(ng)], .str, "")
  dom <- if (!is.null(x$dom_str)) x$dom_str else if (ng < length(x$gens)) paste0(.domain_of(x$rep), "[", paste(vapply(x$gens[-seq_len(ng)], .str, ""),
                                                  collapse = ","), "]")
         else .domain_of(x$rep)
  tail <- if (!is.null(x$modulus)) paste0("modulus=", x$modulus) else paste0("domain='", dom, "'")
  cat("Poly(", .poly_terms_str(x), ", ", paste(gens, collapse = ", "), ", ", tail, ")\n", sep = "")
  invisible(x)
}

# SymPy's StrPrinter._print_Poly: terms in lex order, each monomial written
# in generator order.
.poly_terms_str <- function(x) {
  ms <- monoms(x); cs <- coeffs(x)
  if (!length(ms)) return("0")
  gs <- vapply(x$gens[seq_len(.poly_ngens(x))], function(g) {
    s <- .str(g)
    if (g$type %in% c("Symbol", "Dummy", "Wild", "Integer", "Constant", "Float")) s else paste0("(", s, ")")
  }, "")
  terms <- character()
  for (i in seq_along(ms)) {
    e <- ms[[i]]; c1 <- cs[[i]]
    mon <- paste(unlist(lapply(seq_along(e), function(j)
      if (e[j] == 1) gs[j] else if (e[j] > 1) paste0(gs[j], "^", e[j]))), collapse = "*")
    if (nzchar(mon) && c1$type != "Add") {
      if (is_one_num(c1)) { terms <- c(terms, "+", mon); next }
      if (identical(c1$key, S_NegativeOne$key)) { terms <- c(terms, "-", mon); next }
    }
    sc <- if (c1$type == "Add" && nzchar(mon)) paste0("(", .str(c1), ")") else .str(c1)
    st <- if (nzchar(mon)) paste0(sc, "*", mon) else sc
    terms <- if (startsWith(st, "-")) c(terms, "-", substring(st, 2L)) else c(terms, "+", st)
  }
  if (terms[1L] == "-") terms[2L] <- paste0("-", terms[2L])
  paste(terms[-1L], collapse = " ")
}

#' @export
format.Poly <- function(x, ...) .str(as_expr(x))

#' @export
Ops.Poly <- function(e1, e2) {
  # a generator that is not a symbol (exp(1/2)) stays independent: the
  # representations are combined directly, as Poly arithmetic does in SymPy
  if (!missing(e2) && inherits(e1, "Poly") && is.null(e1$modulus) && any(vapply(e1$gens, function(z) z$type != "Symbol", NA))) {
    same <- inherits(e2, "Poly") && is.null(e2$modulus) && .poly_ngens(e1) == .poly_ngens(e2) &&
      identical(vapply(e1$gens, .k, ""), vapply(e2$gens, .k, ""))
    op <- if (same) switch(.Generic, `+` = qp_add, `-` = qp_sub, `*` = qp_mul)
          else if (.Generic == "^" && is.numeric(e2) && length(e2) == 1L && e2 >= 0 && e2 == round(e2)) function(a, b) qp_pow(a, as.integer(e2))
    if (!is.null(op)) { out <- unclass(e1); out$rep <- op(e1$rep, if (same) e2$rep); return(structure(out, class = "Poly")) }
  }
  g <- if (inherits(e1, "Poly")) e1$gens else e2$gens
  a <- if (inherits(e1, "Poly")) as_expr(e1) else as_sym(e1)
  if (missing(e2)) return(Poly(get(.Generic)(a), g))
  b <- if (inherits(e2, "Poly")) as_expr(e2) else as_sym(e2)
  r <- get(.Generic)(a, b)
  if (.Generic %in% c("==", "!=")) return(r)
  if (.Generic == "/") {
    q <- qp_exquo(Poly(a, g)$rep, Poly(b, g)$rep)
    if (!is.null(q)) return(structure(list(rep = q, gens = g), class = "Poly"))
    return(r)
  }
  do.call(Poly, c(list(r), g))
}

.poly_args <- function(f, gens) {
  if (inherits(f, "Poly")) return(f)
  do.call(Poly, c(list(as_sym(f)), gens))
}

.gen_index <- function(p, gen) {
  if (is.null(gen)) return(1L)
  if (is.numeric(gen)) return(as.integer(gen))
  gen <- as_sym(gen)
  j <- match(gen$key, vapply(p$gens, .k, ""))
  if (is.na(j)) return(NA_integer_)
  j
}

#' Degrees and coefficients of polynomials
#'
#' `degree(f, x)` is the degree in `x` (`-Inf` for zero), `total_degree()`
#' the total degree, `degree_list()` the degree in each generator. `LC()`,
#' `LM()` and `LT()` give the leading coefficient, monomial and term in lex
#' order; `coeffs()` the nonzero coefficients; `all_coeffs()` the dense
#' coefficient list of a univariate polynomial (highest degree first);
#' `monoms()` the exponent vectors; `coeff()` the coefficient of a monomial.
#'
#' @param f An expression or [Poly()].
#' @param gen A generator.
#' @param ... Generators.
#' @param x A monomial such as `x^2*y`, for `coeff()`.
#' @param n Power of `x`, for `coeff()`.
#' @return A number, expression or list.
#' @examples
#' x <- Symbol("x"); y <- Symbol("y")
#' degree(x^3 + x*y^4, x)
#' all_coeffs(Poly(2*x^3 - x + 5, x))
#' coeff(3*x^2*y + x, x, 2)
#' @export
degree <- function(f, gen = NULL) {
  if (!inherits(f, "Poly")) {
    f <- as_sym(f)
    if (is_zero_num(f)) return(-Inf)
    if (is.null(gen)) {
      fs <- free_symbols(f)
      if (!length(fs)) return(0)
      if (length(fs) > 1L) stop("specify the generator for a multivariate expression")
      gen <- fs[[1L]]
    }
    p <- .poly_args(f, list(gen))
  } else {
    p <- f
  }
  if (qp_is_zero(p$rep)) return(-Inf)
  j <- .gen_index(p, gen)
  if (is.na(j)) return(0)
  qp_degree(p$rep, j)
}

#' @rdname degree
#' @export
total_degree <- function(f, ...) {
  p <- .poly_args(f, list(...))
  ng <- .poly_ngens(p)
  if (qp_is_zero(p$rep)) return(-Inf)
  max(rowSums(p$rep$e[, seq_len(ng), drop = FALSE]))
}

#' @rdname degree
#' @export
degree_list <- function(f, ...) {
  p <- .poly_args(f, list(...))
  vapply(seq_len(.poly_ngens(p)), function(j) qp_degree(p$rep, j), 0)
}

# Coefficients as expressions in the generators of a Poly (other symbols go
# into the coefficients). Returns list(monoms = matrix, coeffs = list).
.poly_terms_expr <- function(p) {
  ng <- .poly_ngens(p)
  rep <- p$rep
  if (qp_is_zero(rep)) return(list(monoms = matrix(0L, 0L, ng), coeffs = list()))
  me <- rep$e[, seq_len(ng), drop = FALSE]
  keys <- .mono_keys(me)
  u <- unique(keys)
  monoms <- me[match(u, keys), , drop = FALSE]
  o <- .mono_order_index(monoms, "lex")
  u <- u[o]
  monoms <- monoms[o, , drop = FALSE]
  extra <- p$gens[-seq_len(ng)]
  coeffs <- lapply(u, function(k) {
    w <- which(keys == k)
    sub <- qp(rep$e[w, -seq_len(ng), drop = FALSE], rep$c[w])
    if (!length(extra)) return(make_rational(sum(rep$c[w])))
    qpoly_to_expr(qp_norm(sub$e, sub$c), extra)
  })
  list(monoms = monoms, coeffs = coeffs)
}

#' @rdname degree
#' @export
LC <- function(f, ...) {
  p <- .poly_args(f, list(...))
  t <- .poly_terms_expr(p)
  if (!length(t$coeffs)) return(S_Zero)
  t$coeffs[[1L]]
}

#' @rdname degree
#' @export
LM <- function(f, ...) {
  p <- .poly_args(f, list(...))
  t <- .poly_terms_expr(p)
  if (!length(t$coeffs)) return(S_Zero)
  .monomial_vec_expr(t$monoms[1L, ], p$gens)
}

#' @rdname degree
#' @export
LT <- function(f, ...) .mul(list(LC(f, ...), LM(f, ...)))

.monomial_vec_expr <- function(m, gens) {
  .mul(lapply(which(m > 0L), function(j) .pow(gens[[j]], .int(m[j]))))
}

#' @rdname degree
#' @export
coeffs <- function(f, ...) {
  if (is_sym(f) && f$type == "AlgebraicNumber")
    return(if (length(f$extra$rep)) lapply(as.list(f$extra$rep), make_rational) else list(S_Zero))
  .poly_terms_expr(.poly_args(f, list(...)))$coeffs
}

#' @rdname degree
#' @export
monoms <- function(f, ...) {
  t <- .poly_terms_expr(.poly_args(f, list(...)))
  lapply(seq_len(nrow(t$monoms)), function(i) t$monoms[i, ])
}

#' @rdname degree
#' @export
poly_terms <- function(f, ...) {
  t <- .poly_terms_expr(.poly_args(f, list(...)))
  lapply(seq_along(t$coeffs), function(i) list(monom = t$monoms[i, ], coeff = t$coeffs[[i]]))
}

#' @rdname degree
#' @export
all_coeffs <- function(f, ...) {
  p <- .poly_args(f, list(...))
  if (.poly_ngens(p) != 1L) stop("all_coeffs() needs a univariate polynomial")
  t <- .poly_terms_expr(p)
  if (!length(t$coeffs)) return(list(S_Zero))
  d <- max(t$monoms[, 1L])
  out <- rep(list(S_Zero), d + 1L)
  for (i in seq_along(t$coeffs)) out[[d - t$monoms[i, 1L] + 1L]] <- t$coeffs[[i]]
  out
}

#' @rdname degree
#' @export
coeff <- function(f, x, n = 1) {
  f <- as_sym(f)
  if (.is_seq(f)) return(.seq_coeff(f, x))
  x <- as_sym(x)
  # special targets, in SymPy's order: zero, the expression itself, one
  if (is_zero_num(x)) return(S_Zero)
  if (identical(x$key, f$key)) return(if (n == 1) S_One else S_Zero)
  terms <- if (f$type == "Add") f$args else list(f)
  if (identical(x$key, S_One$key))
    return(.add(Filter(function(t) identical(as_coeff_mul_split(t)[[1L]]$key, S_One$key), terms)))
  if (identical(as_sym(n)$key, S_Zero$key)) {
    if (x$type == "Add" && f$type == "Add") {
      c <- coeff(f, x)
      if (is_zero_num(c)) return(S_Zero)
      return(.sub(f, .add(lapply(if (c$type == "Add") c$args else list(c), function(a) .mul(list(a, x))))))
    }
    return(.add(Filter(function(t) !.has1(t, x), terms)))
  }
  if (n != 1) x <- .pow(x, as_sym(n))
  if (f$type == "Add" && !.is_comm(f) && !.is_comm(x)) {
    # SymPy first keeps the terms having every factor of x as a factor
    margs <- function(t) vapply(if (t$type == "Mul") t$args else list(t), .k, "")
    xk <- margs(x)
    ncs <- function(t) vapply(Filter(Negate(.is_comm), if (t$type == "Mul") t$args else list(t)), .k, "")
    xn <- ncs(x)
    contig <- function(t) {
      tn <- ncs(t)
      if (!length(xn)) return(TRUE)
      if (length(xn) > length(tn)) return(FALSE)
      any(vapply(seq_len(length(tn) - length(xn) + 1L), function(i) all(tn[i:(i + length(xn) - 1L)] == xn), TRUE))
    }
    terms <- Filter(function(t) all(xk %in% margs(t)) && contig(t), terms)
    if (!length(terms)) return(S_Zero)
    f <- .add(terms)
  }
  .coeff_match(f, x, terms)
}

# SymPy Expr.coeff matching: a term contributes when the factors of x are
# among its factors; noncommutative factors must appear in order.
.coeff_match <- function(f, x, terms) {
  cnc <- function(t) {
    fs <- if (t$type == "Mul") t$args else list(t)
    # args_cnc splits a negative number into -1 and its magnitude
    if (is_number_type(fs[[1L]]) && num_sign(fs[[1L]]) < 0 && !identical(fs[[1L]]$key, S_NegativeOne$key))
      fs <- c(list(S_NegativeOne, .neg(fs[[1L]])), fs[-1L])
    comm <- vapply(fs, .is_comm, TRUE)
    list(c = fs[comm], nc = fs[!comm])
  }
  self_c <- .is_comm(f); x_c <- .is_comm(x)
  if (self_c && !x_c) return(S_Zero)
  xs <- cnc(x)
  xk <- vapply(xs$c, .k, "")
  one_c <- self_c || x_c
  co <- list()
  for (a in terms) {
    m <- cnc(a)
    mk <- vapply(m$c, .k, "")
    if (length(xk) > length(mk) || !all(xk %in% mk)) next
    resid <- m$c[!(mk %in% xk)]
    co[[length(co) + 1L]] <- list(r = resid, nc = m$nc)
  }
  if (!length(co)) return(S_Zero)
  if (one_c) return(.add(lapply(co, function(z) .mul(c(z$r, z$nc)))))
  nx <- vapply(xs$nc, .k, "")
  # SymPy's default (left) coefficient takes the last occurrence of x
  find <- function(l) {
    lk <- vapply(l, .k, "")
    if (!length(nx) || length(nx) > length(lk)) return(NULL)
    for (i in rev(seq_len(length(lk) - length(nx) + 1L))) if (all(lk[i:(i + length(nx) - 1L)] == nx)) return(i)
    NULL
  }
  nck <- vapply(co, function(z) paste(vapply(z$nc, .k, ""), collapse = "|"), "")
  if (all(nck == nck[1L])) {
    ii <- find(co[[1L]]$nc)
    if (!is.null(ii)) return(.mul(c(list(.add(lapply(co, function(z) .mul(z$r)))), co[[1L]]$nc[seq_len(ii - 1L)])))
  }
  common <- function(ls) Reduce(function(a, b) {
    n <- min(length(a), length(b))
    k <- 0L
    while (k < n && identical(a[[k + 1L]]$key, b[[k + 1L]]$key)) k <- k + 1L
    a[seq_len(k)]
  }, ls)
  beg <- common(lapply(co, `[[`, "nc"))
  if (length(beg)) {
    ii <- find(beg)
    if (!is.null(ii)) {
      gk <- Reduce(intersect, lapply(co, function(z) vapply(z$r, .k, "")))
      g <- Filter(function(f) f$key %in% gk, co[[1L]]$r)
      return(.mul(c(g, beg[seq_len(ii - 1L)])))
    }
  }
  end <- rev(common(lapply(co, function(z) rev(z$nc))))
  if (length(end)) {
    ii <- find(end)
    if (!is.null(ii)) {
      cut <- -length(end) + ii - 1L
      return(.add(lapply(co, function(z) {
        keep <- if (cut < 0L) length(z$nc) + cut else 0L
        .mul(c(z$r, z$nc[seq_len(max(0L, keep))]))
      })))
    }
  }
  hits <- Filter(function(z) !is.null(find(z$nc)), co)
  if (length(hits) == 1L) {
    z <- hits[[1L]]
    return(.mul(c(z$r, z$nc[seq_len(find(z$nc) - 1L)])))
  }
  S_Zero
}

# ---------------------------------------------------------------------------
# Division

.common_polys <- function(f, g, gens) {
  if (inherits(f, "Poly")) gens <- c(gens, f$gens)
  if (inherits(g, "Poly")) gens <- c(gens, g$gens)
  fe <- if (inherits(f, "Poly")) as_expr(f) else as_sym(f)
  ge <- if (inherits(g, "Poly")) as_expr(g) else as_sym(g)
  gens <- .as_gens(gens)
  if (!is.null(gens)) gens <- gens[!duplicated(vapply(gens, .k, ""))]
  cv <- .to_qpolys(list(fe, ge), gens)
  list(f = cv$polys[[1L]], g = cv$polys[[2L]], gens = cv$gens,
       poly_out = inherits(f, "Poly") || inherits(g, "Poly"))
}

.wrap_out <- function(rep, gens, as_poly) {
  if (as_poly) return(structure(list(rep = rep, gens = gens), class = "Poly"))
  qpoly_to_expr(rep, gens)
}

# Univariate arithmetic in x over expression coefficients, for calls that
# name the generator while other symbols act as coefficients (the
# fraction field QQ(y) of SymPy).
.ex_coeff_mode <- function(f, g, gens) {
  if (length(gens) != 1L || inherits(f, "Poly") || inherits(g, "Poly")) return(NULL)
  x <- tryCatch(as_sym(gens[[1L]]), error = function(e) NULL)
  if (is.null(x) || !is_symbol_like(x)) return(NULL)
  others <- setdiff(names(c(.free_sym_list(as_sym(f)), .free_sym_list(as_sym(g)))), x$key)
  if (length(others)) return(x)
  # constant coefficients that are not rational (sqrt(2), pi, floats)
  cf <- tryCatch(c(.coeff_list(as_sym(f), x), .coeff_list(as_sym(g), x)), error = function(e) NULL)
  if (!is.null(cf) && !all(vapply(cf, function(c) c$type %in% c("Integer", "Rational"), TRUE))) return(x)
  NULL
}
.ex_trim <- function(c) {
  c <- lapply(c, cancel)
  zero <- function(z) is_zero_num(z) || (z$type == "Float" && z$value == 0)
  while (length(c) && zero(c[[1L]])) c <- c[-1L]
  c
}
.ex_to_expr <- function(c, x) {
  n <- length(c)
  if (!n) return(S_Zero)
  .add(lapply(seq_len(n), function(i) .mul(list(c[[i]], .pow(x, .int(n - i))))))
}
.ex_sub <- function(a, b) {
  n <- max(length(a), length(b))
  pad <- function(v) c(rep(list(S_Zero), n - length(v)), v)
  .ex_trim(Map(.sub, pad(a), pad(b)))
}
.ex_mul <- function(a, b) {
  if (!length(a) || !length(b)) return(list())
  out <- rep(list(S_Zero), length(a) + length(b) - 1L)
  for (i in seq_along(a)) for (j in seq_along(b))
    out[[i + j - 1L]] <- .add(list(out[[i + j - 1L]], .mul(list(a[[i]], b[[j]]))))
  .ex_trim(out)
}
.ex_divmod <- function(a, b) {
  a <- .ex_trim(a); b <- .ex_trim(b)
  if (length(a) < length(b)) return(list(q = list(), r = a))
  nq <- length(a) - length(b) + 1L
  q <- rep(list(S_Zero), nq)
  while (length(a) && length(a) >= length(b)) {
    sh <- length(a) - length(b)
    k <- cancel(.div(a[[1L]], b[[1L]]))
    q[[nq - sh]] <- k
    na <- .ex_sub(a, c(lapply(b, function(z) .mul(list(k, z))), rep(list(S_Zero), sh)))
    # the leading term must go, even when rounding leaves a residue
    if (length(na) >= length(a)) na <- na[-1L]
    a <- na
  }
  list(q = .ex_trim(q), r = a)
}

#' Polynomial division
#'
#' `div(f, g)` returns `list(q, r)` with `f = q*g + r`, using the first
#' generator as main variable over the rationals. `quo()` and `rem()` return
#' one part; `exquo()` fails unless the division is exact. `pdiv()`,
#' `prem()` and `pquo()` perform pseudo division. `reduced()` divides by a
#' list of polynomials.
#'
#' @param f,g Expressions or [Poly()] objects.
#' @param ... Generators.
#' @return Expressions, or `Poly` objects when inputs are `Poly`.
#' @examples
#' x <- Symbol("x")
#' div(x^3 + 2*x + 1, x - 1)
#' rem(x^4 + 1, x^2 + 1)
#' @export
div <- function(f, g, ...) {
  pa <- .split_poly_opts(list(...))
  if (length(pa$opts)) {
    Fd <- .poly_field(pa$opts, list(f, g))
    if (Fd$kind != "QQ") return(.field_binary(f, g, Fd, pa$gens, "div"))
    return(do.call(div, c(list(f, g), pa$gens)))
  }
  x <- .ex_coeff_mode(f, g, list(...))
  if (!is.null(x)) {
    d <- .ex_divmod(.coeff_list(as_sym(f), x), .coeff_list(as_sym(g), x))
    return(list(q = expand(.ex_to_expr(d$q, x)), r = expand(.ex_to_expr(d$r, x))))
  }
  cp <- .common_polys(f, g, list(...))
  if (qp_is_zero(cp$g)) stop("polynomial division by zero")
  vars_g <- qp_vars_used(cp$g)
  vars_f <- qp_vars_used(cp$f)
  if (length(union(vars_f, vars_g)) <= 1L) {
    i <- if (length(vars_f)) vars_f[1L] else if (length(vars_g)) vars_g[1L] else 1L
    d <- up_divmod(up_from_qpoly(cp$f, i), up_from_qpoly(cp$g, i))
    n <- length(cp$gens)
    return(list(q = .wrap_out(up_to_qpoly(d$q, i, n), cp$gens, cp$poly_out),
                r = .wrap_out(up_to_qpoly(d$r, i, n), cp$gens, cp$poly_out)))
  }
  d <- qp_div_multi(cp$f, list(cp$g))
  list(q = .wrap_out(d$q[[1L]], cp$gens, cp$poly_out),
       r = .wrap_out(d$r, cp$gens, cp$poly_out))
}

#' @rdname div
#' @export
quo <- function(f, g, ...) div(f, g, ...)$q

#' @rdname div
#' @export
rem <- function(f, g, ...) div(f, g, ...)$r

#' @rdname div
#' @export
exquo <- function(f, g, ...) {
  d <- div(f, g, ...)
  r <- if (inherits(d$r, "Poly")) as_expr(d$r) else d$r
  if (!is_zero_num(r)) stop(ExactQuotientFailed(paste0(.str(as_expr(g)), " does not divide ", .str(as_expr(f)))))
  d$q
}

#' @rdname div
#' @export
pdiv <- function(f, g, ...) {
  cp <- .common_polys(f, g, list(...))
  i <- 1L
  dg <- qp_degree(cp$g, i)
  df <- qp_degree(cp$f, i)
  lc <- .lc_poly_in_var(cp$g, i)
  k <- max(df - dg + 1, 0)
  scaled <- qp_mul(cp$f, qp_pow(lc, k))
  d <- qp_div_multi(scaled, list(cp$g))
  list(q = .wrap_out(d$q[[1L]], cp$gens, cp$poly_out), r = .wrap_out(d$r, cp$gens, cp$poly_out))
}

#' @rdname div
#' @export
prem <- function(f, g, ...) pdiv(f, g, ...)$r

#' @rdname div
#' @export
pquo <- function(f, g, ...) pdiv(f, g, ...)$q

#' @rdname div
#' @param G A list of polynomials.
#' @param order Monomial order: `"lex"`, `"grlex"` or `"grevlex"`.
#' @export
reduced <- function(f, G, ..., order = "lex") {
  f <- as_sym(f)
  G <- lapply(G, as_sym)
  gens <- .gens_of(c(list(f), G), list(...))
  old <- .ORDER$current
  .ORDER$current <- order
  on.exit(.ORDER$current <- old)
  cv <- .to_qpolys(c(list(f), G), gens)
  ps <- lapply(cv$polys, qp_reorder, order = order)
  d <- qp_div_multi(ps[[1L]], ps[-1L], order)
  list(q = lapply(d$q, qpoly_to_expr, gens = cv$gens), r = qpoly_to_expr(d$r, cv$gens))
}

# ---------------------------------------------------------------------------
# GCD, LCM and friends

.num_gcd <- function(a, b) {
  a <- as_sym(a)
  b <- as_sym(b)
  if (a$type == "Integer" && b$type == "Integer") {
    return(make_rational(gmp::as.bigq(gmp::gcd.bigz(gmp::numerator(a$value), gmp::numerator(b$value)))))
  }
  pa <- a$value
  pb <- b$value
  make_rational(gmp::as.bigq(gmp::gcd.bigz(gmp::numerator(pa), gmp::numerator(pb)),
                             gmp::lcm.bigz(gmp::denominator(pa), gmp::denominator(pb))))
}

#' Greatest common divisor and least common multiple
#'
#' Works for integers, rationals and polynomials in any number of variables
#' (heuristic GCD with a primitive PRS fallback, as in SymPy). The result is
#' normalized to a positive leading coefficient with integer content.
#' `gcd_list()` and `lcm_list()` take a list. `gcdex()` returns `s`, `t` and
#' `h` with `s*f + t*g = h` for univariate polynomials; `invert()` computes
#' the inverse of `f` modulo `g`.
#'
#' @param f,g Expressions or [Poly()] objects.
#' @param ... Generators.
#' @param seq A list of expressions.
#' @return An expression.
#' @examples
#' x <- Symbol("x"); y <- Symbol("y")
#' gcd(x^2 - 1, x^2 + 2*x + 1)
#' gcd(x^2*y - y, x*y^2 + y^2)
#' sym_lcm(S(4), S(6))
#' @export
gcd <- function(f, g, ...) {
  fs <- as_sym(if (inherits(f, "Poly")) as_expr(f) else f)
  gs <- as_sym(if (inherits(g, "Poly")) as_expr(g) else g)
  if (is_exact_number(fs) && is_exact_number(gs)) return(.num_gcd(fs, gs))
  if (!is.null(P <- .poly_mod_gcd(f, g, list(...)))) return(P)
  pa <- .split_poly_opts(list(...))
  if (length(pa$opts)) {
    Fd <- .poly_field(pa$opts, list(fs, gs))
    if (Fd$kind != "QQ") return(.field_binary(fs, gs, Fd, pa$gens, "gcd"))
    return(do.call(gcd, c(list(f, g), pa$gens)))
  }
  # SymPy's domains: over RR each Float becomes a nearby rational
  # (RealField.to_rational), the GCD is computed over QQ, and its
  # coefficients come back as Floats; over QQ the GCD is monic; over ZZ it
  # keeps the integer content.
  fl <- if (inherits(f, "Poly") || inherits(g, "Poly")) list() else c(.float_atoms(fs), .float_atoms(gs))
  if (length(fl)) {
    prec <- max(vapply(fl, .float_prec, 0))
    rule <- lapply(fl, function(v) list(v, make_rational(.rr_to_rational(.mp_to_q(.float_mp(v)), prec))))
    f <- xreplace(fs, rule); g <- xreplace(gs, rule)
  }
  cp <- .common_polys(f, g, list(...))
  h <- qp_gcd(cp$f, cp$g)
  zz <- all(gmp::denominator(c(cp$f$c, cp$g$c, gmp::as.bigq(1))) == 1)
  if (length(fl) || !zz) {
    h <- if (qp_is_zero(cp$f) && qp_is_zero(cp$g)) qp_zero(length(cp$gens)) else qp_monic(h)
    if (length(fl)) return(.qp_float_expr(h, cp$gens, prec))
  } else if (qp_is_const(h)) {
    cf <- qp_content(cp$f)
    cg <- qp_content(cp$g)
    h <- qp_const(gmp::as.bigq(gmp::gcd.bigz(gmp::numerator(cf), gmp::numerator(cg)),
                               gmp::lcm.bigz(gmp::denominator(cf), gmp::denominator(cg))),
                  length(cp$gens))
    if (qp_is_zero(cp$f) && qp_is_zero(cp$g)) h <- qp_zero(length(cp$gens))
  } else {
    cf <- qp_content(cp$f)
    cg <- qp_content(cp$g)
    cgc <- gmp::as.bigq(gmp::gcd.bigz(gmp::numerator(cf), gmp::numerator(cg)),
                        gmp::lcm.bigz(gmp::denominator(cf), gmp::denominator(cg)))
    if (qp_is_zero(cp$f)) cgc <- cg
    if (qp_is_zero(cp$g)) cgc <- cf
    h <- qp_scale(qp_primitive(h)$poly, cgc)
  }
  .wrap_out(h, cp$gens, cp$poly_out)
}

# Float atoms of an expression.
.float_atoms <- function(e) {
  if (e$type == "Float") return(list(e))
  out <- list()
  for (a in e$args) if (is_sym(a)) out <- c(out, .float_atoms(a))
  out[!duplicated(vapply(out, `[[`, "", "key"))]
}

# RealField.to_rational(limit = TRUE): the exact rational q of a Float, or
# the closest rational with denominator at most max(2^prec %/% 200, 99).
.rr_to_rational <- function(q, prec = 53) {
  maxd <- max(gmp::as.bigz(2)^prec %/% 200, gmp::as.bigz(99))
  p <- gmp::numerator(q); d <- gmp::denominator(q)
  if (d <= maxd) return(q)
  p0 <- gmp::as.bigz(0); q0 <- gmp::as.bigz(1); p1 <- gmp::as.bigz(1); q1 <- gmp::as.bigz(0)
  n <- p
  repeat {
    a <- gmp::divq.bigz(n, d)
    q2 <- q0 + a * q1
    if (q2 > maxd) break
    t0 <- p1; t1 <- q1
    p1 <- p0 + a * p1; q1 <- q2; p0 <- t0; q0 <- t1
    t <- d; d <- n - a * d; n <- t
  }
  k <- gmp::divq.bigz(maxd - q0, q1)
  b1 <- gmp::as.bigq(p0 + k * p1, q0 + k * q1)
  b2 <- gmp::as.bigq(p1, q1)
  if (b2 == 0 || b1 == 0) return(q)
  if (abs(b2 - q) <= abs(b1 - q)) b2 else b1
}

# A qpoly over QQ as an expression with Float coefficients (RR.as_expr).
.qp_float_expr <- function(p, gens, prec) {
  if (qp_is_zero(p)) return(.mp_float(.mp_from_q(gmp::as.bigq(0), prec), prec))
  terms <- lapply(seq_along(p$c), function(i) {
    parts <- list(.mp_float(.mp_from_q(p$c[i], prec), prec))
    for (j in which(p$e[i, ] > 0L)) parts[[length(parts) + 1L]] <- .pow(gens[[j]], .int(p$e[i, j]))
    .mul(parts)
  })
  .add(terms)
}

#' @noRd
lcm <- function(f, g, ...) {
  if (missing(f) || missing(g)) return(.base_call(quote(graphics::lcm)))
  fs <- as_sym(if (inherits(f, "Poly")) as_expr(f) else f)
  gs <- as_sym(if (inherits(g, "Poly")) as_expr(g) else g)
  if (is_exact_number(fs) && is_exact_number(gs)) {
    if (is_zero_num(fs) || is_zero_num(gs)) return(S_Zero)
    return(make_rational(abs(fs$value * gs$value)) / .num_gcd(fs, gs))
  }
  pa <- .split_poly_opts(list(...))
  if (length(pa$opts)) {
    Fd <- .poly_field(pa$opts, list(fs, gs))
    if (Fd$kind != "QQ") return(.field_binary(fs, gs, Fd, pa$gens, "lcm"))
    return(do.call(lcm, c(list(f, g), pa$gens)))
  }
  cp <- .common_polys(f, g, list(...))
  h <- qp_gcd(cp$f, cp$g)
  l <- qp_exquo(qp_mul(qp_primitive(cp$f)$poly, qp_primitive(cp$g)$poly), qp_primitive(h)$poly)
  # over QQ SymPy's LCM is monic
  if (!all(gmp::denominator(c(cp$f$c, cp$g$c, gmp::as.bigq(1))) == 1))
    return(.wrap_out(if (qp_is_zero(cp$f) || qp_is_zero(cp$g)) qp_zero(length(cp$gens)) else qp_monic(l),
                     cp$gens, cp$poly_out))
  cf <- qp_content(cp$f)
  cg <- qp_content(cp$g)
  cl <- gmp::as.bigq(gmp::lcm.bigz(gmp::numerator(cf), gmp::numerator(cg)),
                     gmp::gcd.bigz(gmp::denominator(cf), gmp::denominator(cg)))
  .wrap_out(qp_scale(qp_primitive(l)$poly, cl), cp$gens, cp$poly_out)
}

#' @rdname gcd
#' @export
gcd_list <- function(seq, ...) Reduce(function(a, b) gcd(a, b, ...), seq)

#' @rdname gcd
#' @export
lcm_list <- function(seq, ...) Reduce(function(a, b) lcm(a, b, ...), seq)

#' @rdname gcd
#' @export
gcdex <- function(f, g, ...) {
  x <- .ex_coeff_mode(f, g, list(...))
  if (!is.null(x)) {
    a <- .ex_trim(.coeff_list(as_sym(f), x)); b <- .ex_trim(.coeff_list(as_sym(g), x))
    s0 <- list(S_One); s1 <- list(); t0 <- list(); t1 <- list(S_One)
    while (length(b)) {
      d <- .ex_divmod(a, b)
      a <- b; b <- d$r
      tmp <- .ex_sub(s0, .ex_mul(d$q, s1)); s0 <- s1; s1 <- tmp
      tmp <- .ex_sub(t0, .ex_mul(d$q, t1)); t0 <- t1; t1 <- tmp
    }
    lc <- a[[1L]]
    nz <- function(c) .ex_to_expr(lapply(c, function(z) cancel(.div(z, lc))), x)
    return(list(s = nz(s0), t = nz(t0), h = nz(a)))
  }
  cp <- .common_polys(f, g, list(...))
  i <- 1L
  a <- up_from_qpoly(cp$f, i)
  b <- up_from_qpoly(cp$g, i)
  s0 <- gmp::as.bigq(1); s1 <- .bq0()
  t0 <- .bq0(); t1 <- gmp::as.bigq(1)
  while (length(b)) {
    d <- up_divmod(a, b)
    a <- b
    b <- d$r
    tmp <- up_sub(s0, up_mul(d$q, s1)); s0 <- s1; s1 <- tmp
    tmp <- up_sub(t0, up_mul(d$q, t1)); t0 <- t1; t1 <- tmp
  }
  lc <- a[1L]
  x <- cp$gens[[i]]
  list(s = up_to_expr(s0 / lc, x), t = up_to_expr(t0 / lc, x), h = up_to_expr(a / lc, x))
}

#' @rdname gcd
#' @export
invert <- function(f, g, ...) {
  fs <- as_sym(f)
  gs <- as_sym(g)
  if (fs$type == "Integer" && gs$type == "Integer") {
    return(make_rational(gmp::as.bigq(gmp::inv.bigz(gmp::numerator(fs$value), gmp::numerator(gs$value)))))
  }
  r <- gcdex(f, g, ...)
  if (!is_one_num(r$h)) stop("polynomial is not invertible modulo g")
  r$s
}

#' Content and primitive part
#'
#' `primitive(f)` returns `list(content, primitive)`, and `terms_gcd()` pulls
#' out the greatest common monomial and numeric factor.
#'
#' @param f An expression.
#' @param ... Generators.
#' @return A list or an expression.
#' @export
primitive <- function(f, ...) {
  f <- as_sym(f)
  gens <- .gens_of(list(f), list(...))
  p <- expr_to_qpoly(f, gens, strict = FALSE)
  pr <- qp_primitive(p)
  list(content = make_rational(pr$content), primitive = qpoly_to_expr(pr$poly, gens))
}

#' @rdname primitive
#' @export
content <- function(f, ...) primitive(f, ...)$content

#' @rdname primitive
#' @export
terms_gcd <- function(f, ...) {
  f <- as_sym(f)
  gens <- .gens_of(list(f), list(...))
  p <- expr_to_qpoly(f, gens, strict = FALSE)
  if (qp_is_zero(p)) return(S_Zero)
  pr <- qp_primitive(p)
  mins <- apply(pr$poly$e, 2L, min)
  q <- qp_norm(pr$poly$e - matrix(mins, nrow(pr$poly$e), length(mins), byrow = TRUE), pr$poly$c)
  .mul(list(make_rational(pr$content), .monomial_vec_expr(mins, gens), qpoly_to_expr(q, gens)))
}

# ---------------------------------------------------------------------------
# Square-free decomposition and factorization

#' @noRd
factor <- function(x, ...) {
  if (missing(x)) return(.base_call(quote(base::factor)))
  if (is_sym(x) && .is_bd(x)) return(.bd_factor(x, ...))
  if (is_sym(x) || inherits(x, c("Poly", "symmatrix")) ||
      (is.list(x) && length(x) && all(vapply(x, is_sym, TRUE)))) {
    return(factor_sym(x, ...))
  }
  base::factor(x, ...)
}

#' @rdname sym_factor
#' @export
factor_sym <- function(f, ..., deep = FALSE) {
  if (inherits(f, "pyobj")) return(.pyobj_factor(f, ..., deep = deep))
  if (inherits(f, "symmatrix")) return(.mat_map(f, function(z) factor_sym(z, ..., deep = deep)))
  if (is.list(f) && !inherits(f, "sym") && !inherits(f, "Poly")) return(lapply(f, factor_sym, ..., deep = deep))
  if (inherits(f, "Poly")) f <- as_expr(f)
  f <- as_sym(f)
  if (is_relational(f)) return(sym_rebuild(f, lapply(f$args, factor_sym, ..., deep = deep)))
  # SymPy's factor() hands noncommutative input to factor_nc()
  if (!.is_comm(f)) return(factor_nc(f))
  pa <- .split_poly_opts(list(...))
  if (length(pa$opts)) {
    Fd <- .poly_field(pa$opts, list(f))
    if (Fd$kind != "QQ") return(.field_factor(f, Fd, pa$gens, "factor"))
    return(do.call(factor_sym, c(list(f), pa$gens, list(deep = deep))))
  }
  if (deep) {
    # factor(deep = True): bottom up, then the products and sums once more
    gens <- list(...)
    try_factor <- function(z) { fac <- do.call(factor_sym, c(list(z), gens)); if (fac$type %in% c("Mul", "Pow")) fac else z }
    f <- .bottom_up_args(f, try_factor)
    partials <- list()
    for (p in unique(c(atoms(f, "Mul"), atoms(f, "Add")))) {
      fac <- do.call(factor_sym, c(list(p), gens))
      if (fac$type %in% c("Mul", "Pow") && !identical(fac$key, p$key)) partials[[length(partials) + 1L]] <- list(p, fac)
    }
    return(if (length(partials)) xreplace(f, partials) else f)
  }
  .symbolic_factor(f, list(...))
}

# _symbolic_factor and _symbolic_factor_list of sympy/polys/polytools.py:
# every factor of together(expr) is factored on its own, numeric bases with
# symbolic exponents are kept, and the content is kept apart (_keep_coeff)
.symbolic_factor <- function(expr, gens, Fd = NULL) {
  if (!.is_expr_node(expr)) return(sym_rebuild(expr, lapply(expr$args, .symbolic_factor, gens = gens, Fd = Fd)))
  r <- .symbolic_factor_list(together(expr), gens, Fd)
  # a numerator expanding to zero is a factor 0 in SymPy (PolificationFailed)
  if (is_zero_num(r$coef)) return(S_Zero)
  .keep_coeff(r$coef, .mul(lapply(r$factors, function(p) .pow(p[[1L]], p[[2L]]))))
}

.symbolic_factor_list <- function(expr, gens, Fd = NULL) {
  coeff <- S_One; factors <- list()
  # the third element, when present, carries the domain and generators of a
  # factor (Polys over ZZ_I or QQ_I) for _sorted_factors
  addf <- function(b, e, m = NULL) factors[[length(factors) + 1L]] <<- list(b, e, m)
  for (arg in .margs(expr)) {
    if (is_number_type(arg) || .pure_complex(arg)) { coeff <- .mul(list(coeff, arg)); next }
    if (arg$type == "Pow" && !is_const(arg$args[[1L]], "E")) {
      base <- arg$args[[1L]]; ex <- arg$args[[2L]]
      if (is_number_type(base) && is_number_type(ex)) { coeff <- .mul(list(coeff, arg)); next }
      if (is_number_type(base)) { addf(base, ex); next }
    } else { base <- arg; ex <- S_One }
    fp <- .factor_parts(base, gens, Fd)
    if (is.null(fp)) { addf(base, ex); next }
    c_ <- fp$coef; fs <- fp$factors
    if (!(is_one_num(c_) && c_$type == "Integer")) {
      if (ex$type == "Integer") coeff <- .mul(list(coeff, .pow(c_, ex)))
      else if (isTRUE(ask_fact(c_, "positive"))) addf(c_, ex)
      else fs[[length(fs) + 1L]] <- list(c_, S_One)
    }
    qm <- function(q) if (length(q) >= 3L) q[[3L]] else NULL
    if (is_one_num(ex) && ex$type == "Integer") { for (q in fs) addf(q[[1L]], q[[2L]], qm(q)) }
    else if (isTRUE(ask_fact(ex, "integer"))) { for (q in fs) addf(q[[1L]], .mul(list(q[[2L]], ex)), qm(q)) }
    else {
      other <- list()
      for (q in fs) if (isTRUE(ask_fact(q[[1L]], "positive"))) addf(q[[1L]], .mul(list(q[[2L]], ex))) else other[[length(other) + 1L]] <- q
      addf(.mul(lapply(other, function(q) .pow(q[[1L]], q[[2L]]))), ex)
    }
  }
  # duplicates collected
  kb <- list(); ke <- list(); km <- list()
  for (q in factors) {
    k <- q[[1L]]$key
    if (is.null(kb[[k]])) { kb[[k]] <- q[[1L]]; ke[[k]] <- q[[2L]]; km[k] <- list(q[[3L]]) } else ke[[k]] <- .add(list(ke[[k]], q[[2L]]))
  }
  list(coef = coeff, factors = lapply(names(kb), function(k) if (is.null(km[[k]])) list(kb[[k]], ke[[k]]) else list(kb[[k]], ke[[k]], km[[k]])))
}

# pure_complex: a + b*I with numbers a and b, b nonzero
.pure_complex <- function(v) {
  ca <- .coeff_add_split(v); t <- ca[[2L]]
  if (is_zero_num(t)) return(FALSE)
  sp <- as_coeff_mul_split(t)
  is_const(sp[[2L]], "I") && is_number_type(sp[[1L]]) && is_number_type(ca[[1L]])
}

# Coefficient and factors of e, or NULL when e is not a polynomial; with Fd
# the factorization over that field.
.factor_parts <- function(e, gens, Fd = NULL) {
  e <- expand(e)
  if (is_number_type(e)) return(list(coef = e, factors = list()))
  gens2 <- .gens_of(list(e), gens)
  if (!length(gens2)) return(NULL)
  # factor(..., extension = , gaussian = , domain = ): over the field asked for
  if (!is.null(Fd)) {
    r <- .field_factor_list(e, Fd, gens, "factor")
    return(list(coef = r$coeff, factors = lapply(r$factors, function(p) list(p[[1L]], as_sym(p[[2L]])))))
  }
  # dmp_factor_list over RR: factor exactly over QQ, then each factor is
  # divided by its max norm and converted back
  # generators keep their Floats (log(x - 68.5) stays a generator); the
  # domain is RR only for Float coefficients
  fg <- if (.has_type(e, "Float")) Filter(function(g) !is_number_type(g) && .has_type(g, "Float"), gens2) else list()
  fd <- lapply(fg, function(g) Dummy())
  back <- function(z) if (length(fg)) xreplace(z, lapply(seq_along(fg), function(k) list(fd[[k]], fg[[k]]))) else z
  e0 <- if (length(fg)) xreplace(e, lapply(seq_along(fg), function(k) list(fg[[k]], fd[[k]]))) else e
  if (.has_type(e0, "Float")) {
    ex <- .map_bottom_up(e0, function(z) if (z$type == "Float") make_rational(gmp::as.bigq(z$value)) else z)
    gx <- Filter(function(g) !is_number_type(g), .gens_of(list(ex), gens))
    q <- tryCatch(expr_to_qpoly(ex, gx, strict = TRUE), error = function(err) NULL)
    if (!is.null(q) && length(gx)) {
      fl <- qp_factor_list(q)
      coef <- fl$coef; fs <- list()
      for (pr in fl$factors) {
        mn <- max(abs(pr[[1L]]$c))
        fe <- .add(lapply(seq_along(pr[[1L]]$c), function(i)
          .mul(list(Float(as.double(pr[[1L]]$c[i] / mn)), .monomial_vec_expr(pr[[1L]]$e[i, ], gx)))))
        fs[[length(fs) + 1L]] <- list(back(fe), .int(pr[[2L]]))
        coef <- coef * mn^pr[[2L]]
      }
      return(list(coef = Float(as.double(coef)), factors = fs))
    }
  }
  # ZZ_I and QQ_I: coefficients a + b*I with rational a and b
  if (.contains_key(e, S_I$key) && !any(vapply(gens, function(g) is_const(as_sym(g), "I"), NA))) {
    gx <- Filter(function(g) !is_const(g, "I"), gens2)
    iv <- Dummy("i")
    q <- if (length(gx)) tryCatch(expr_to_qpoly(xreplace(e, list(list(S_I, iv))), c(gx, list(iv)), strict = TRUE), error = function(err) NULL) else NULL
    if (!is.null(q)) {
      # construct_domain: ZZ_I for Gaussian integer coefficients, else QQ_I
      Fd <- list(kind = "NF", exts = list(S_I), gauss = if (all(gmp::denominator(q$c) == 1)) "ZZ_I" else "QQ_I")
      r <- .field_factor_list(e, Fd, gens, "factor")
      # the factors stay Polys over ZZ_I or QQ_I in gx for _sorted_factors
      meta <- list(dom = Fd$gauss, gens = gx)
      return(list(coef = r$coeff, factors = lapply(r$factors, function(p) list(p[[1L]], as_sym(p[[2L]]), meta))))
    }
  }
  p <- tryCatch(expr_to_qpoly(e, gens2, strict = FALSE), error = function(err) NULL)
  if (is.null(p) && length(gens)) {
    # other generators form the ground domain (ZZ[y]): factor over all of
    # them; the factors free of the given generators belong to the coefficient
    gk <- vapply(gens2, `[[`, "", "key")
    allg <- c(gens2, Filter(function(g) !(g$key %in% gk), .find_gens(list(e))))
    q <- tryCatch(expr_to_qpoly(e, allg, strict = FALSE), error = function(err) NULL)
    if (is.null(q)) return(NULL)
    fl <- qp_factor_list(q)
    coef <- make_rational(fl$coef); fs <- list()
    for (pr in fl$factors) {
      fe <- qpoly_to_expr(pr[[1L]], allg)
      if (any(vapply(gens2, function(g) .has1(fe, g), NA))) fs[[length(fs) + 1L]] <- list(fe, .int(pr[[2L]]))
      else coef <- .mul(list(coef, .pow(fe, .int(pr[[2L]]))))
    }
    return(list(coef = coef, factors = fs))
  }
  if (is.null(p)) return(NULL)
  fl <- qp_factor_list(p)
  list(coef = make_rational(fl$coef), factors = lapply(fl$factors, function(pr) list(qpoly_to_expr(pr[[1L]], gens2), .int(pr[[2L]]))))
}

.map_args_deep <- function(e, fn) {
  if (!length(e$args)) return(e)
  if (e$type == "Function") return(sym_rebuild(e, lapply(e$args, fn)))
  sym_rebuild(e, lapply(e$args, function(a) .map_args_deep(a, fn)))
}

.factor_expr <- function(e, gens, sqf_only = FALSE) {
  e <- expand(e)
  if (is_number_type(e) || !length(free_symbols(e)) && !.has_type(e, "Function")) {
    if (is_number_type(e)) return(e)
  }
  gens2 <- .gens_of(list(e), gens)
  if (!length(gens2)) return(e)
  p <- expr_to_qpoly(e, gens2, strict = FALSE)
  fl <- if (sqf_only) .qp_sqf_full(p) else qp_factor_list(p)
  parts <- list(make_rational(fl$coef))
  for (pr in fl$factors) parts[[length(parts) + 1L]] <- .pow(.sign_normalize(qpoly_to_expr(pr[[1L]], gens2)), .int(pr[[2L]]))
  .mul(parts)
}

.sign_normalize <- function(e) e

.qp_sqf_full <- function(p) {
  n <- qp_nvars(p)
  pr <- qp_primitive(p)
  f <- pr$poly
  out <- list()
  mins <- apply(f$e, 2L, min)
  if (length(f$c) && any(mins > 0L)) {
    for (j in which(mins > 0L)) out[[length(out) + 1L]] <- list(qp_var(j, n), mins[j])
    f <- qp_norm(f$e - matrix(mins, nrow(f$e), n, byrow = TRUE), f$c)
  }
  pending <- list(f)
  while (length(pending)) {
    g <- pending[[1L]]
    pending <- pending[-1L]
    if (qp_is_const(g)) next
    v <- qp_vars_used(g)[1L]
    ct <- qp_content_in_var(g, v)
    if (!qp_is_const(ct)) {
      pending <- c(pending, list(ct))
      g <- qp_exquo(g, ct)
    }
    out <- c(out, .qp_sqf_list_var(g, v))
  }
  if (length(out) > 1L) {
    keys <- vapply(out, function(pr) paste(.mono_keys(pr[[1L]]$e), as.character(pr[[1L]]$c), collapse = ";"), "")
    merged <- list()
    for (u in unique(keys)) {
      w <- which(keys == u)
      merged[[length(merged) + 1L]] <- list(out[[w[1L]]][[1L]], sum(vapply(out[w], `[[`, 0, 2L)))
    }
    out <- merged
  }
  prod_all <- qp_const(1, n)
  for (q in out) prod_all <- qp_mul(prod_all, qp_pow(q[[1L]], q[[2L]]))
  coef <- pr$content * (if (length(pr$poly$c)) qp_lc(pr$poly) / qp_lc(prod_all) else 1)
  list(coef = coef, factors = out)
}

#' @rdname sym_factor
#' @export
factor_list <- function(f, ...) {
  f <- as_sym(if (inherits(f, "Poly")) as_expr(f) else f)
  pa <- .split_poly_opts(list(...))
  if (length(pa$opts)) {
    Fd <- .poly_field(pa$opts, list(f))
    if (Fd$kind != "QQ") return(.field_factor_list(f, Fd, pa$gens, "factor"))
    return(do.call(factor_list, c(list(f), pa$gens)))
  }
  # _generic_factor_list: the numerator and denominator of together(f),
  # each through _symbolic_factor_list
  egens <- pa$gens
  nd <- .as_numer_denom(together(f))
  cp <- .symbolic_factor_list(nd[[1L]], egens)
  cq <- .symbolic_factor_list(nd[[2L]], egens)
  if (length(cq$factors)) stop(structure(class = c("PolynomialError", "error", "condition"), list(message = paste("a polynomial expected, got", .str(f)), call = NULL)))
  fs <- .sympy_sorted_factors(cp$factors, "factor", if (length(egens)) egens else NULL)
  # integer multiplicities as R integers, symbolic ones as expressions
  fs <- lapply(fs, function(p) list(p[[1L]], if (is_sym(p[[2L]]) && p[[2L]]$type == "Integer") as.integer(as.character(p[[2L]]$value)) else p[[2L]]))
  return(list(coeff = .div(cp$coef, cq$coef), factors = fs))
  gens <- .gens_of(list(f), list(...))
  p <- tryCatch(expr_to_qpoly(f, gens, strict = FALSE), error = function(err) NULL)
  if (is.null(p) && length(list(...))) {
    # other symbols form the ground domain (ZZ[a]): factor over all of them
    # and move the factors free of the generators into the coefficient
    gk <- vapply(gens, `[[`, "", "key")
    allg <- c(gens, Filter(function(g) !(g$key %in% gk), .find_gens(list(expand(f)))))
    fl <- qp_factor_list(expr_to_qpoly(f, allg, strict = FALSE))
    coeff <- make_rational(fl$coef); facs <- list()
    for (pr in fl$factors) {
      fe <- qpoly_to_expr(pr[[1L]], allg)
      if (any(vapply(gens, function(g) .has1(fe, g), TRUE))) facs[[length(facs) + 1L]] <- list(fe, pr[[2L]])
      else coeff <- .mul(list(coeff, .pow(fe, .int(pr[[2L]]))))
    }
    return(list(coeff = coeff, factors = .sympy_sorted_factors(facs, "factor", gens)))
  }
  if (is.null(p)) p <- expr_to_qpoly(f, gens, strict = FALSE)
  fl <- qp_factor_list(p)
  facs <- lapply(fl$factors, function(pr) list(qpoly_to_expr(pr[[1L]], gens), pr[[2L]]))
  list(coeff = make_rational(fl$coef),
       factors = .sympy_sorted_factors(facs, "factor", if (length(list(...))) gens else NULL))
}

#' @rdname sym_factor
#' @export
sqf_list <- function(f, ...) {
  f <- as_sym(if (inherits(f, "Poly")) as_expr(f) else f)
  pa <- .split_poly_opts(list(...))
  if (length(pa$opts)) {
    Fd <- .poly_field(pa$opts, list(f))
    if (Fd$kind != "QQ") return(.field_factor_list(f, Fd, pa$gens, "sqf"))
    return(do.call(sqf_list, c(list(f), pa$gens)))
  }
  gens <- .gens_of(list(f), list(...))
  p <- expr_to_qpoly(f, gens, strict = FALSE)
  fl <- .qp_sqf_full(p)
  # dmp_sqf_list: one factor per multiplicity, by increasing multiplicity
  ks <- vapply(fl$factors, function(pr) as.numeric(pr[[2L]]), 0)
  facs <- lapply(sort(unique(ks)), function(k) {
    fk <- fl$factors[ks == k]
    pk <- Reduce(qp_mul, lapply(fk, `[[`, 1L))
    list(qpoly_to_expr(pk, gens), k)
  })
  list(coeff = make_rational(fl$coef), factors = facs)
}

#' @rdname sym_factor
#' @export
sqf <- function(f, ...) {
  f <- as_sym(f)
  pa <- .split_poly_opts(list(...))
  if (length(pa$opts)) {
    Fd <- .poly_field(pa$opts, list(f))
    if (Fd$kind != "QQ") return(.field_factor(f, Fd, pa$gens, "sqf"))
    return(do.call(sqf, c(list(f), pa$gens)))
  }
  nd <- fraction(together(f))
  .div(.factor_expr(nd[[1L]], list(...), sqf_only = TRUE), .factor_expr(nd[[2L]], list(...), sqf_only = TRUE))
}

#' @rdname sym_factor
#' @export
sqf_part <- function(f, ...) {
  if (inherits(f, "Poly")) return(.sqf_part_poly(f, list(...)))
  fe <- as_sym(if (inherits(f, "Poly")) as_expr(f) else f)
  if (is_zero_num(fe)) return(S_Zero)
  pa <- .split_poly_opts(list(...))
  if (length(pa$opts)) {
    Fd <- .poly_field(pa$opts, list(fe))
    if (Fd$kind == "NF" && !identical(Fd$gauss, "ZZ_I")) return(.field_sqf_part(fe, Fd, pa$gens))
  }
  sl <- sqf_list(f, ...)
  .sqf_part_expr(.mul(lapply(sl$factors, `[[`, 1L)), pa)
}

# ---------------------------------------------------------------------------
# Resultants and discriminants

#' Resultant and discriminant
#'
#' `resultant(f, g, x)` is the resultant with respect to `x`, computed as the
#' determinant of the Sylvester matrix. `discriminant(f, x)` is
#' `(-1)^(n(n-1)/2) * resultant(f, f') / lc(f)`.
#'
#' @param f,g Expressions.
#' @param x The variable.
#' @return An expression.
#' @examples
#' x <- Symbol("x"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c")
#' discriminant(a * x^2 + b * x + c, x)
#' resultant(x^2 - 1, x - 2, x)
#' @export
resultant <- function(f, g, x = NULL) {
  f <- as_sym(f)
  g <- as_sym(g)
  if (is.null(x)) x <- free_symbols(.add(list(f, g)))[[1L]]
  x <- as_sym(x)
  a <- .coeff_list(f, x)
  b <- .coeff_list(g, x)
  m <- length(a) - 1L
  n <- length(b) - 1L
  if (m < 0 || n < 0) return(S_Zero)
  if (m == 0) return(.pow(a[[1L]], .int(n)))
  if (n == 0) return(.pow(b[[1L]], .int(m)))
  N <- m + n
  rows <- vector("list", N)
  for (i in seq_len(n)) {
    r <- rep(list(S_Zero), N)
    r[i:(i + m)] <- a
    rows[[i]] <- r
  }
  for (i in seq_len(m)) {
    r <- rep(list(S_Zero), N)
    r[i:(i + n)] <- b
    rows[[n + i]] <- r
  }
  expand(.det_bareiss(rows))
}

#' @rdname resultant
#' @export
discriminant <- function(f, x = NULL) {
  f <- as_sym(f)
  if (is.null(x)) x <- free_symbols(f)[[1L]]
  x <- as_sym(x)
  a <- .coeff_list(f, x)
  n <- length(a) - 1L
  if (n < 1L) return(S_Zero)
  r <- resultant(f, .diff1(f, x), x)
  s <- if ((n * (n - 1L) / 2L) %% 2L == 1L) S_NegativeOne else S_One
  cancel(.mul(list(s, r, .inv(a[[1L]]))))
}

# Dense coefficient list (highest degree first) of f as a polynomial in x,
# with coefficients as expressions.
.coeff_list <- function(f, x) {
  f <- expand(as_sym(f))
  terms <- if (f$type == "Add") f$args else list(f)
  degs <- integer(length(terms))
  cfs <- vector("list", length(terms))
  for (i in seq_along(terms)) {
    t <- terms[[i]]
    d <- .degree_in(t, x)
    if (is.na(d)) stop("not a polynomial in ", .str(x), ": ", .str(f))
    degs[i] <- d
    cfs[[i]] <- if (d == 0L) t else .div(t, .pow(x, .int(d)))
  }
  if (!length(terms) || (length(terms) == 1L && is_zero_num(terms[[1L]]))) return(list())
  D <- max(degs)
  out <- rep(list(S_Zero), D + 1L)
  for (i in seq_along(terms)) out[[D - degs[i] + 1L]] <- .add(list(out[[D - degs[i] + 1L]], cfs[[i]]))
  out
}

# Degree of a single term in x, NA if the term is not polynomial in x.
.degree_in <- function(t, x) {
  if (!.has1(t, x)) return(0L)
  if (identical(t$key, x$key)) return(1L)
  if (t$type == "Pow" && identical(t$args[[1L]]$key, x$key) && t$args[[2L]]$type == "Integer" &&
      t$args[[2L]]$value > 0) return(as.integer(as.character(t$args[[2L]]$value)))
  if (t$type == "Mul") {
    d <- 0L
    for (f in t$args) {
      df <- .degree_in(f, x)
      if (is.na(df)) return(NA_integer_)
      d <- d + df
    }
    return(d)
  }
  NA_integer_
}

#' Is an expression a polynomial?
#'
#' @param f An expression.
#' @param ... Generators; by default all free symbols.
#' @return `TRUE` or `FALSE`.
#' @export
is_polynomial <- function(f, ...) {
  f <- expand(as_sym(f))
  gens <- .as_gens(list(...))
  if (is.null(gens)) gens <- free_symbols(f)
  terms <- if (f$type == "Add") f$args else list(f)
  for (t in terms) for (g in gens) if (is.na(.degree_in(t, g))) return(FALSE)
  TRUE
}

# ---------------------------------------------------------------------------
# Composition, Horner form, interpolation

#' Polynomial composition, Horner form and interpolation
#'
#' `compose(f, g, x)` computes `f(g)`. `horner(f, x)` rewrites in Horner
#' form. `interpolate(points, x)` returns the interpolating polynomial for a
#' list of `c(x, y)` pairs, or for values at `x = 1, 2, ...`.
#'
#' @param f,g Expressions.
#' @param x The variable.
#' @param points A list of `c(x, y)` pairs, or a vector of values.
#' @return An expression.
#' @examples
#' x <- Symbol("x")
#' horner(x^3 + 2*x^2 + 3*x + 4, x)
#' interpolate(list(c(1, 1), c(2, 4), c(3, 9)), x)
#' @export
compose <- function(f, g, x = NULL) {
  f <- as_sym(f)
  if (is.null(x)) x <- free_symbols(f)[[1L]]
  expand(subs(f, as_sym(x), as_sym(g)))
}

#' @rdname compose
#' @export
horner <- function(f, x = NULL) {
  f <- as_sym(f)
  if (is.null(x)) x <- free_symbols(f)[[1L]]
  x <- as_sym(x)
  cl <- .coeff_list(f, x)
  if (!length(cl)) return(S_Zero)
  acc <- cl[[1L]]
  if (length(cl) > 1L) for (k in 2:length(cl)) acc <- .add(list(.mul(list(acc, x)), cl[[k]]))
  acc
}

#' @rdname compose
#' @export
interpolate <- function(points, x) {
  x <- as_sym(x)
  if (!is.list(points) || !.interp_is_pairs(points)) points <- lapply(seq_along(points), function(i) list(i, points[[i]]))
  xs <- lapply(points, function(p) as_sym(p[[1L]]))
  ys <- lapply(points, function(p) as_sym(p[[2L]]))
  n <- length(xs)
  terms <- vector("list", n)
  for (i in seq_len(n)) {
    num <- list(ys[[i]])
    for (j in seq_len(n)) if (j != i) num[[length(num) + 1L]] <- .div(.sub(x, xs[[j]]), .sub(xs[[i]], xs[[j]]))
    terms[[i]] <- .mul(num)
  }
  expand(.add(terms))
}

# _sorted_factors of polytools: each factor as its own Poly (its own
# generators unless gens are given), ordered by (len(rep), len(gens), exp,
# str(domain), rep) for factor and (exp, len(rep), ...) for sqf
.sympy_sorted_factors <- function(factors, method = "factor", gens = NULL) {
  if (length(factors) < 2L) return(factors)
  info <- lapply(factors, function(pr) {
    fe <- expand(as_sym(pr[[1L]]))
    if (length(pr) >= 3L && !is.null(pr[[3L]])) {
      # a Poly over ZZ_I or QQ_I: dense in its generators, Gaussian coefficients
      m <- pr[[3L]]
      rep <- .expr_dense_rep(fe, m$gens)$rep
      return(list(lenrep = length(rep), ngens = length(m$gens), exp = as_sym(pr[[2L]]), dom = m$dom, rep = rep, gauss = TRUE))
    }
    g <- if (is.null(gens)) .poly_sort_gens(.find_gens(list(fe))) else gens
    p <- tryCatch(expr_to_qpoly(fe, g, strict = FALSE), error = function(e) NULL)
    if (is.null(p) && length(g)) {
      # Poly(fe, *gens) over a ground domain: the dense representation in the
      # generators with the coefficients as ground elements
      cr <- .expr_dense_rep(fe, g)
      rep <- cr$rep; nd <- format(construct_domain(if (length(cr$coeffs)) cr$coeffs else list(S_Zero))[[1L]])
    } else {
      rep <- if (is.null(p) || !length(g)) list(fe) else .qp_dense_rep(p, 1L)
      nd <- if (is.null(p)) "EX" else {
        cs <- if (qp_is_zero(p)) list(S_Zero) else lapply(seq_along(p$c), function(i) make_rational(p$c[i]))
        format(construct_domain(cs)[[1L]])
      }
    }
    ex <- as_sym(pr[[2L]])
    list(lenrep = length(rep), ngens = length(g), exp = ex, dom = nd, rep = rep)
  })
  cmp1 <- function(u, v) {
    if (is_sym(u) || is_sym(v)) {
      u <- as_sym(u); v <- as_sym(v)
      if (identical(u$key, v$key)) return(0)
      if (.is_true(Lt(u, v))) return(-1)
      if (.is_true(Gt(u, v))) return(1)
      return(0)
    }
    if (identical(u, v)) 0 else if (u < v) -1 else 1
  }
  less <- function(a, b) {
    ka <- if (method == "sqf") list(a$exp, a$lenrep, a$ngens, a$dom) else list(a$lenrep, a$ngens, a$exp, a$dom)
    kb <- if (method == "sqf") list(b$exp, b$lenrep, b$ngens, b$dom) else list(b$lenrep, b$ngens, b$exp, b$dom)
    for (i in seq_along(ka)) {
      c <- if (is.character(ka[[i]])) (if (ka[[i]] == kb[[i]]) 0 else if (.cmp_codepoints(ka[[i]], kb[[i]]) < 0) -1 else 1) else cmp1(ka[[i]], kb[[i]])
      if (c != 0) return(c < 0)
    }
    if (isTRUE(a$gauss) && isTRUE(b$gauss)) return(.rep_cmp_gauss(a$rep, b$rep) < 0)
    .rep_cmp(a$rep, b$rep) < 0
  }
  ord <- seq_along(factors)
  for (i in seq_along(ord)[-1L]) {
    j <- i
    while (j > 1L && less(info[[ord[j]]], info[[ord[j - 1L]]])) { t <- ord[j]; ord[j] <- ord[j - 1L]; ord[j - 1L] <- t; j <- j - 1L }
  }
  factors[ord]
}

# the dense recursive representation of a qpoly from variable i on
.qp_dense_rep <- function(p, i) {
  n <- qp_nvars(p)
  if (i > n) return(if (qp_is_zero(p)) gmp::as.bigq(0) else p$c[1L])
  if (qp_is_zero(p)) return(list())
  d <- max(p$e[, i])
  lapply(d:0, function(k) {
    rows <- which(p$e[, i] == k)
    sub <- qp(p$e[rows, , drop = FALSE], p$c[rows])
    sub$e[, i] <- 0L
    r <- .qp_dense_rep(sub, i + 1L)
    if (i == n) { if (!length(rows)) gmp::as.bigq(0) else p$c[rows[1L]] } else r
  })
}
# dense reps over ZZ_I or QQ_I: GaussianElement.__lt__ compares [y, x]
.rep_cmp_gauss <- function(a, b) {
  if (is_sym(a) && is_sym(b)) {
    ya <- as_sym(im(a)); yb <- as_sym(im(b))
    if (!identical(ya$key, yb$key)) return(if (.is_true(Lt(ya, yb))) -1 else 1)
    xa <- as_sym(re(a)); xb <- as_sym(re(b))
    if (!identical(xa$key, xb$key)) return(if (.is_true(Lt(xa, xb))) -1 else 1)
    return(0)
  }
  if (!is.list(a) || !is.list(b)) return(0)
  n <- min(length(a), length(b))
  for (i in seq_len(n)) { c <- .rep_cmp_gauss(a[[i]], b[[i]]); if (c != 0) return(c) }
  length(a) - length(b)
}
.rep_cmp <- function(a, b) {
  # ground elements: EX elements compare by sort_key, PolyElements by
  # (number of terms, terms)
  if (is_sym(a) && is_sym(b)) return(.ground_elem_cmp(a, b))
  if (is_sym(a) || is_sym(b)) {
    a <- if (is_sym(a)) a else make_rational(gmp::as.bigq(a)); b <- if (is_sym(b)) b else make_rational(gmp::as.bigq(b))
    return(.ground_elem_cmp(a, b))
  }
  if (!is.list(a) && !is.list(b)) {
    if (!length(a) || !length(b)) return(length(a) - length(b))
    return(if (a < b) -1 else if (a > b) 1 else 0)
  }
  if (!is.list(a) || !is.list(b)) return(0)
  n <- min(length(a), length(b))
  for (i in seq_len(n)) { c <- .rep_cmp(a[[i]], b[[i]]); if (c != 0) return(c) }
  if (length(a) < length(b)) -1 else if (length(a) > length(b)) 1 else 0
}

# Poly(e, *gens).rep for a ground domain holding the other generators: the
# dense representation with ground coefficients as expressions
.expr_dense_rep <- function(e, gens) {
  gk <- vapply(gens, `[[`, "", "key")
  allg <- c(gens, Filter(function(g) !(g$key %in% gk), .find_gens(list(e))))
  q <- tryCatch(expr_to_qpoly(e, allg, strict = FALSE), error = function(err) NULL)
  if (is.null(q)) return(list(rep = list(e), coeffs = list(e)))
  ng <- length(gens); coeffs <- list()
  rec <- function(rows, i) {
    if (i > ng) {
      sub <- qp(q$e[rows, , drop = FALSE], q$c[rows]); sub$e[, seq_len(ng)] <- 0L
      ce <- qpoly_to_expr(qp_norm(sub$e, sub$c), allg)
      coeffs[[length(coeffs) + 1L]] <<- ce
      return(ce)
    }
    if (!length(rows)) return(list())
    d <- max(q$e[rows, i])
    lapply(d:0, function(k) { r2 <- rows[q$e[rows, i] == k]; if (!length(r2) && i == ng) S_Zero else rec(r2, i + 1L) })
  }
  list(rep = rec(seq_len(nrow(q$e)), 1L), coeffs = coeffs)
}
.ground_elem_cmp <- function(a, b) {
  if (identical(a$key, b$key)) return(0)
  if (is_number_type(a) && is_number_type(b)) return(if (.is_true(Lt(a, b))) -1 else 1)
  fa <- free_symbols(a); fb <- free_symbols(b)
  if (length(fa) || length(fb)) {
    g <- .poly_sort_gens(.find_gens(list(a, b)))
    pa <- tryCatch(expr_to_qpoly(expand(a), g, strict = FALSE), error = function(err) NULL)
    pb <- tryCatch(expr_to_qpoly(expand(b), g, strict = FALSE), error = function(err) NULL)
    if (!is.null(pa) && !is.null(pb)) {
      # PolyElement.sort_key: (len(self), self.terms())
      na <- length(pa$c); nb <- length(pb$c)
      if (na != nb) return(if (na < nb) -1 else 1)
      oa <- .qp_terms_lex(pa); ob <- .qp_terms_lex(pb)
      for (i in seq_len(na)) {
        ma <- oa$e[i, ]; mb <- ob$e[i, ]
        if (!identical(ma, mb)) { j <- which(ma != mb)[1L]; return(if (ma[j] < mb[j]) -1 else 1) }
        if (oa$c[i] != ob$c[i]) return(if (oa$c[i] < ob$c[i]) -1 else 1)
      }
      return(0)
    }
  }
  .key_cmp(.sort_key(a), .sort_key(b))
}
.qp_terms_lex <- function(p) {
  o <- do.call(order, c(lapply(seq_len(ncol(p$e)), function(j) -p$e[, j]), list(method = "radix")))
  list(e = p$e[o, , drop = FALSE], c = p$c[o])
}
