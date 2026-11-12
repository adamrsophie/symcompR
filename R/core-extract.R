# Extraction methods of sympy/core/expr.py: extract_multiplicatively,
# extract_additively, as_coefficient, extract_branch_factor, and
# Add.primitive of sympy/core/add.py.

#' Extracting factors and terms
#'
#' `extract_multiplicatively(e, c)` writes `e` as `c` times something "in a
#' nice way" (keeping the properties of the arguments) and returns that
#' something, or `NULL`, as SymPy's `Expr.extract_multiplicatively()`.
#' `extract_additively(e, c)` does the same for a term. `as_coefficient(e,
#' expr)` is the coefficient of `expr` in `e` when `e` is a multiple of
#' `expr`, else `NULL`. `extract_branch_factor(e)` writes `e` as
#' `exp_polar(2*pi*I*n)*z` and returns `list(z, n)`. `add_primitive(e)`
#' returns the rational content of a sum and the primitive remainder, as
#' SymPy's `Add.primitive()`.
#'
#' @param e An expression.
#' @param c,expr The factor or term to extract.
#' @param allow_half Also extract `exp_polar(I*pi)`, as half a branch.
#' @return An expression, `NULL`, or a list.
#' @examples
#' x <- Symbol("x"); y <- Symbol("y")
#' extract_multiplicatively(2 * x, 2)
#' extract_multiplicatively((x * y)^3, x^2 * y)
#' extract_additively(x + 3, 2)
#' as_coefficient(2 * sym_I * x, sym_I)
#' extract_branch_factor(exp_polar(3 * sym_pi * sym_I + x))
#' @name extract-methods
NULL

#' @rdname extract-methods
#' @export
add_primitive <- function(e) {
  e <- as_sym(e)
  if (e$type != "Add") {
    # Expr.primitive: the positive rational coefficient
    sp <- as_coeff_mul_split(e)
    if (!(sp[[1L]]$type %in% c("Integer", "Rational"))) return(list(S_One, e))
    if (num_sign(sp[[1L]]) < 0) return(list(.neg(sp[[1L]]), .neg(sp[[2L]])))
    return(list(sp[[1L]], sp[[2L]]))
  }
  terms <- lapply(e$args, function(a) {
    sp <- as_coeff_mul_split(a)
    if (!(sp[[1L]]$type %in% c("Integer", "Rational"))) return(list(q = gmp::as.bigq(1), m = a))
    list(q = gmp::as.bigq(sp[[1L]]$value), m = sp[[2L]])
  })
  num <- Reduce(gmp::gcd.bigz, lapply(terms, function(t) abs(gmp::numerator(t$q))), gmp::as.bigz(0))
  den <- Reduce(gmp::lcm.bigz, lapply(terms, function(t) gmp::denominator(t$q)), gmp::as.bigz(1))
  if (num == 1 && den == 1) return(list(S_One, e))
  r <- gmp::as.bigq(num, den)
  newt <- lapply(terms, function(t) .keep_coeff(make_rational(t$q / r), t$m))
  list(make_rational(r), .add(newt))
}

#' @rdname extract-methods
#' @export
extract_multiplicatively <- function(e, c) {
  e <- as_sym(e); c <- as_sym(c)
  if (is_nan(e)) return(NULL)
  if (is_one_num(c)) return(e)
  if (identical(c$key, e$key)) return(S_One)
  if (c$type == "Add") {
    pr <- add_primitive(c)
    if (!is_one_num(pr[[1L]])) c <- new_sym("Mul", list(pr[[1L]], pr[[2L]]))
  }
  if (c$type == "Mul") {
    a <- c$args[[1L]]; b <- if (length(c$args) == 2L) c$args[[2L]] else new_sym("Mul", c$args[-1L])
    x <- extract_multiplicatively(e, a)
    return(if (is.null(x)) NULL else extract_multiplicatively(x, b))
  }
  quotient <- .div(e, c)
  fa <- function(z, nm) isTRUE(ask_fact(z, nm))
  if (is_number_type(e) || is_inf_const(e)) {
    if (is_const(e, "oo")) return(if (fa(c, "positive")) S_Infinity else NULL)
    if (is_const(e, "-oo")) return(if (fa(c, "negative")) S_Infinity else if (fa(c, "positive")) S_NegativeInfinity else NULL)
    if (is_const(e, "zoo")) return(if (identical(ask_fact(c, "zero"), FALSE) || (is_number_type(c) && !is_zero_num(c))) S_ComplexInfinity else NULL)
    if (e$type == "Integer") {
      if (quotient$type != "Integer") return(NULL)
      if (num_sign(e) > 0 && num_sign(quotient) < 0) return(NULL)
      return(quotient)
    }
    if (e$type %in% c("Integer", "Rational")) {
      if (!(quotient$type %in% c("Integer", "Rational"))) return(NULL)
      if (num_sign(e) > 0 && num_sign(quotient) < 0) return(NULL)
      return(quotient)
    }
    if (e$type == "Float") {
      if (quotient$type != "Float") return(NULL)
      if (num_sign(e) > 0 && num_sign(quotient) < 0) return(NULL)
      return(quotient)
    }
    return(NULL)
  }
  if ((e$type == "Constant" && !(e$name %in% c("oo", "-oo", "zoo", "nan"))) || is_symbol_like(e)) {
    if (quotient$type == "Mul" && length(quotient$args) == 2L) {
      q1 <- quotient$args[[1L]]
      if (q1$type == "Integer" && num_sign(q1) > 0 && identical(quotient$args[[2L]]$key, e$key)) return(quotient)
    } else if (quotient$type == "Integer" && is_number_type(c)) return(quotient)
    return(NULL)
  }
  if (e$type == "Add") {
    pr <- add_primitive(e); cs <- pr[[1L]]; ps <- pr[[2L]]
    if (is_number_type(c) && !identical(c$key, S_NegativeOne$key)) {
      if (!is_one_num(cs)) {
        xc <- if (num_sign(c) < 0) { r <- extract_multiplicatively(cs, .neg(c)); if (is.null(r)) NULL else .neg(r) }
              else extract_multiplicatively(cs, c)
        if (!is.null(xc)) return(.mul(list(xc, ps)))
      }
      return(NULL)
    }
    if (identical(c$key, ps$key)) return(cs)
    newargs <- list()
    for (a in (if (ps$type == "Add") ps$args else list(ps))) {
      na <- extract_multiplicatively(a, c)
      if (is.null(na)) return(NULL)
      newargs[[length(newargs) + 1L]] <- na
    }
    if (!is_one_num(cs)) return(.add(lapply(newargs, function(t) .mul(list(cs, t)))))
    return(.add(newargs))
  }
  if (e$type == "Mul") {
    args <- e$args
    for (i in seq_along(args)) {
      na <- extract_multiplicatively(args[[i]], c)
      if (!is.null(na)) { args[[i]] <- na; return(.mul(args)) }
    }
    return(NULL)
  }
  is_expf <- function(z) z$type == "Function" && z$name == "exp"
  if (e$type == "Pow" || is_expf(e)) {
    be <- function(z) if (is_expf(z)) list(S_E, z$args[[1L]]) else if (z$type == "Pow") list(z$args[[1L]], z$args[[2L]]) else list(z, S_One)
    s <- be(e); cb <- be(c)
    mk <- function(b, ex) if (identical(b$key, S_E$key)) fcall("exp", list(ex)) else .pow(b, ex)
    if (identical(cb[[1L]]$key, s[[1L]]$key)) {
      ne <- extract_additively(s[[2L]], cb[[2L]])
      if (!is.null(ne)) return(mk(s[[1L]], ne))
    } else if (identical(c$key, s[[1L]]$key)) {
      ne <- extract_additively(s[[2L]], S_One)
      if (!is.null(ne)) return(mk(s[[1L]], ne))
    }
  }
  NULL
}

# _corem: the coefficient of c in eq and the remaining terms
.corem <- function(eq, c) {
  co <- list(); non <- list()
  for (t in (if (eq$type == "Add") eq$args else list(eq))) {
    ci <- .expr_coeff(t, c)
    if (is_zero_num(ci)) non[[length(non) + 1L]] <- t else co[[length(co) + 1L]] <- ci
  }
  list(.add(co), .add(non))
}

# Expr.coeff(x) for one term or a sum (no expansion)
.expr_coeff <- function(e, x) {
  if (e$type == "Add") return(.add(lapply(e$args, .expr_coeff, x = x)))
  if (identical(e$key, x$key)) return(S_One)
  if (e$type == "Mul") {
    xs <- if (x$type == "Mul") x$args else list(x)
    fs <- e$args; keys <- vapply(fs, .k, "")
    rest <- fs
    for (xi in xs) {
      j <- match(xi$key, vapply(rest, .k, ""))
      if (is.na(j)) return(S_Zero)
      rest <- rest[-j]
    }
    r <- .mul(rest)
    if (any(vapply(xs, function(xi) .has1(r, xi), TRUE))) return(S_Zero)
    return(r)
  }
  S_Zero
}

#' @rdname extract-methods
#' @export
extract_additively <- function(e, c) {
  e <- as_sym(e); c <- as_sym(c)
  if (is_nan(e)) return(NULL)
  if (is_zero_num(c)) return(e)
  if (identical(c$key, e$key)) return(S_Zero)
  if (is_zero_num(e)) return(NULL)
  if (is_number_type(e)) {
    if (!is_number_type(c)) return(NULL)
    d <- .sub(e, c)
    co <- num_sign(e); ds <- num_sign(d)
    if ((co > 0 && ds >= 0 && num_cmp(d, e) < 0) || (co < 0 && ds <= 0 && num_cmp(d, e) > 0)) return(d)
    return(NULL)
  }
  coeff_add <- function(z) {
    if (z$type == "Add" && is_number_type(z$args[[1L]])) list(z$args[[1L]], .add(z$args[-1L]))
    else if (is_number_type(z)) list(z, S_Zero) else list(S_Zero, z)
  }
  if (is_number_type(c)) {
    ct <- coeff_add(e)
    xa <- extract_additively(ct[[1L]], c)
    if (is.null(xa)) return(NULL)
    return(.add(list(xa, ct[[2L]])))
  }
  if (c$type == "Add" && is_number_type(c$args[[1L]])) {
    co <- .expr_coeff(e, c)
    x1 <- extract_additively(co, S_One)
    xa0 <- if (is.null(x1)) S_Zero else .mul(list(x1, c))
    if (!is_zero_num(xa0)) {
      d <- .sub(e, .mul(list(co, c)))
      dd <- extract_additively(d, c)
      r <- .add(list(xa0, if (is.null(dd)) d else dd))
      return(if (is_zero_num(r)) NULL else r)
    }
    hc <- coeff_add(c); hs <- coeff_add(e)
    xa <- extract_additively(hs[[1L]], hc[[1L]])
    if (is.null(xa)) return(NULL)
    xa2 <- extract_additively(hs[[2L]], hc[[2L]])
    if (is.null(xa2)) return(NULL)
    return(.add(list(xa, xa2)))
  }
  cr <- .corem(e, c)
  x1 <- extract_additively(cr[[1L]], S_One)
  xa0 <- if (is.null(x1)) S_Zero else .mul(list(x1, c))
  if (!is_zero_num(xa0)) {
    dd <- extract_additively(cr[[2L]], c)
    r <- .add(list(xa0, if (is.null(dd)) cr[[2L]] else dd))
    return(if (is_zero_num(r)) NULL else r)
  }
  coeffs <- list()
  self <- e
  for (a in (if (c$type == "Add") c$args else list(c))) {
    sp <- as_coeff_mul_split(a); ac <- sp[[1L]]; at <- sp[[2L]]
    co <- .expr_coeff(self, at)
    if (is_zero_num(co)) return(NULL)
    ct <- coeff_add(co)
    xa <- extract_additively(ct[[1L]], ac)
    if (is.null(xa)) return(NULL)
    self <- .sub(self, .mul(list(co, at)))
    coeffs[[length(coeffs) + 1L]] <- .mul(list(.add(list(ct[[2L]], xa)), at))
  }
  .add(c(coeffs, list(self)))
}

#' @rdname extract-methods
#' @export
as_coefficient <- function(e, expr) {
  r <- extract_multiplicatively(e, expr)
  if (!is.null(r) && !is_zero_num(r) && !.has1(r, as_sym(expr))) r else NULL
}

#' @rdname extract-methods
#' @export
extract_branch_factor <- function(e, allow_half = FALSE) {
  e <- as_sym(e)
  n <- S_Zero; res <- S_One; exps <- list()
  for (a in (if (e$type == "Mul") e$args else list(e))) {
    if (a$type == "Function" && a$name == "exp_polar") exps[[length(exps) + 1L]] <- a$args[[1L]]
    else res <- .mul(list(res, a))
  }
  piimult <- S_Zero; extras <- list()
  ipi <- .mul(list(S_Pi, S_I))
  while (length(exps)) {
    ex <- exps[[length(exps)]]; exps <- exps[-length(exps)]
    if (ex$type == "Add") { exps <- c(exps, ex$args); next }
    if (ex$type == "Mul") {
      cf <- as_coefficient(ex, ipi)
      if (!is.null(cf)) { piimult <- .add(list(piimult, cf)); next }
    }
    extras[[length(extras) + 1L]] <- ex
  }
  if (!length(free_symbols(piimult))) { coeff <- piimult; tail <- list() }
  else {
    parts <- if (piimult$type == "Add") piimult$args else list(piimult)
    coeff <- .add(Filter(function(t) !length(free_symbols(t)), parts))
    tail <- Filter(function(t) length(free_symbols(t)) > 0L, parts)
  }
  branch <- .mul(list(fcall("ceiling", list(.sub(.div(coeff, S_Two), S_Half))), S_Two))
  n <- .add(list(n, .div(branch, S_Two)))
  cc <- .sub(coeff, branch)
  if (allow_half) {
    nc <- extract_additively(cc, S_One)
    if (!is.null(nc)) { n <- .add(list(n, S_Half)); cc <- nc }
  }
  newexp <- .add(c(list(.mul(list(ipi, .add(c(list(cc), tail))))), extras))
  if (!is_zero_num(newexp)) res <- .mul(list(res, fcall("exp_polar", list(newexp))))
  list(res, n)
}
