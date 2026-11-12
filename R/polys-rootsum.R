# Sums of a function over the roots of a polynomial, after RootSum in
# sympy/polys/rootoftools.py.

.rootsum_w <- Dummy("w")

#' Sums over the roots of a polynomial
#'
#' `RootSum(expr, func, x)` is the sum of `func(r)` over the roots `r` of the
#' polynomial `expr` in `x`, counted with multiplicity, as SymPy's
#' `RootSum`. Linear factors are summed at once; when `func` is a rational
#' function and `auto` is `TRUE` the sum is computed exactly from the
#' coefficients (Viete's formulas), otherwise an unevaluated `RootSum`
#' remains. With `quadratic = TRUE` quadratic factors are summed over their
#' explicit roots. [doit()] sums over explicit roots when [roots()] finds
#' them all, and [N()] sums over numerical roots.
#'
#' @param expr A univariate polynomial.
#' @param func A [Lambda()] of one variable, an R function of one argument,
#'   or `NULL` for the identity.
#' @param x The variable of `expr`; by default its only free symbol.
#' @param auto Evaluate sums of rational functions.
#' @param quadratic Sum over the roots of quadratic factors.
#' @return An expression.
#' @examples
#' x <- Symbol("x")
#' RootSum(x^3 + x + 1, Lambda(x, x^2))
#' RootSum(x^3 + x + 1, Lambda(x, 1/(x + 2)))
#' r <- RootSum(x^5 + x + 3, Lambda(x, log(x + 1)))
#' r
#' N(r)
#' @export
RootSum <- function(expr, func = NULL, x = NULL, auto = TRUE, quadratic = FALSE) {
  expr <- as_sym(expr)
  if (is.null(x)) {
    fs <- free_symbols(expr)
    if (length(fs) != 1L) stop("only univariate polynomials are allowed", call. = FALSE)
    x <- fs[[1L]]
  }
  x <- as_sym(x)
  if (is_zero_num(expand(expr))) return(S_NaN)
  w <- .rootsum_w
  pre <- .preprocess_roots(subs(expr, x, w), w)
  poly <- pre$poly
  fexpr <- if (is.null(func)) w
           else if (is_sym(func) && func$type == "Lambda") {
             if (length(func$args[[1L]]$args) != 1L) stop("expected a univariate function, got ", .str(func), call. = FALSE)
             subs(func$args[[2L]], func$args[[1L]]$args[[1L]], w)
           } else if (is.function(func)) as_sym(func(w))
           else stop("expected a univariate function, got ", format(func), call. = FALSE)
  if (!is_one_num(pre$coeff)) fexpr <- subs(fexpr, w, .mul(list(pre$coeff, w)))
  deg <- .int(degree(poly, w))
  if (!.has1(fexpr, w)) return(.mul(list(deg, fexpr)))
  add_const <- S_Zero; mul_const <- S_One
  if (fexpr$type == "Add") {
    add_const <- .add(Filter(function(a) !.has1(a, w), fexpr$args))
    fexpr <- .add(Filter(function(a) .has1(a, w), fexpr$args))
  }
  if (fexpr$type == "Mul") {
    mul_const <- .mul(Filter(function(a) !.has1(a, w), fexpr$args))
    fexpr <- .mul(Filter(function(a) .has1(a, w), fexpr$args))
  }
  rational <- .is_ratfun(fexpr, w)
  # _pure_factors: factors over the coefficient ring, whose generators are
  # factored too
  fl <- factor_list(poly)
  fl$factors <- Filter(function(fk) .has1(fk[[1L]], w), fl$factors)
  # normalized as polynomials in w: positive leading coefficient in lex order
  fl$factors <- lapply(fl$factors, function(fk) {
    gs <- c(list(w), Filter(function(z) !identical(z$key, w$key), .gens_of(list(fk[[1L]]), list())))
    lc <- coeffs(do.call(Poly, c(list(fk[[1L]]), gs)))[[1L]]
    if (num_sign(lc) < 0) fk[[1L]] <- expand(.neg(fk[[1L]]))
    fk
  })
  terms <- lapply(fl$factors, function(fk) {
    f <- fk[[1L]]; k <- .int(fk[[2L]])
    d <- degree(f, w)
    term <- if (d == 1L) {
      cl <- .coeff_list(f, w)  # descending: c1*w + c0
      subs(fexpr, w, .neg(.div(cl[[2L]], cl[[1L]])))
    } else if (isTRUE(quadratic) && d == 2L) {
      .add(lapply(roots(f, w, multiple = TRUE), function(r) subs(fexpr, w, r)))
    } else if (!rational || !isTRUE(auto)) {
      .rootsum_node(f, fexpr, auto)
    } else {
      .rootsum_rational(f, fexpr)
    }
    .mul(list(k, term))
  })
  .add(list(.mul(list(mul_const, .add(terms))), .mul(list(deg, add_const))))
}

.rootsum_node <- function(poly, fexpr, auto = TRUE)
  new_sym("RootSum", list(poly, Lambda(.rootsum_w, fexpr)), extra = list(auto = isTRUE(auto)))

# preprocess_roots: clear denominators, take the primitive part, scale away
# symbolic coefficients (x^2 - a^2 has the roots of x^2 - 1 times a) and
# integer bases (x^5 + 512*x + 1024 at x = 4*y). Returns the scale and the
# polynomial.
.preprocess_roots <- function(f, w) {
  coeff <- S_One
  f <- expand(f)
  cl <- .coeff_list(f, w)
  den <- Reduce(function(u, v) lcm(u, v), lapply(cl, function(c) .as_numer_denom(together(c))[[2L]]))
  cl <- lapply(cl, function(c) cancel(.mul(list(c, den))))
  g <- Reduce(function(u, v) gcd(u, v), Filter(Negate(is_zero_num), cl))
  cl <- lapply(cl, function(c) cancel(.div(c, g)))
  n <- length(cl) - 1L
  poly <- expand(.add(lapply(seq_along(cl), function(i) .mul(list(cl[[i]], .pow(w, .int(n - i + 1L)))))))
  gens <- Filter(function(z) !identical(z$key, w$key), .gens_of(list(poly), list()))
  is_term <- function(c) { c <- expand(c); is_zero_num(c) || c$type != "Add" }
  if (length(gens) && all(vapply(cl, is_term, TRUE))) {
    P <- do.call(Poly, c(list(poly, w), gens))
    ms <- do.call(rbind, lapply(monoms(P), as.integer))
    base <- ms[, 1L]
    for (k in seq_along(gens)) {
      strip <- ms[, k + 1L]
      rev_ <- strip[1L] < strip[length(strip)]
      if (rev_) strip <- rev(strip)
      ratio <- NULL; ok <- TRUE
      for (i in seq_along(base)) {
        a <- base[i]; b <- strip[i]
        if (!a && !b) next
        if (!a || !b || b %% a != 0) { ok <- FALSE; break }
        r <- b %/% a
        if (is.null(ratio)) ratio <- r else if (ratio != r) { ok <- FALSE; break }
      }
      if (ok && !is.null(ratio)) {
        if (rev_) ratio <- -ratio
        poly <- expand(subs(poly, gens[[k]], S_One))
        coeff <- .mul(list(coeff, .pow(gens[[k]], .int(-ratio))))
      }
    }
  }
  cl <- .coeff_list(poly, w)
  if (all(vapply(cl, function(c) c$type == "Integer", TRUE))) {
    basis <- .integer_basis(lapply(cl, function(c) gmp::as.bigz(c$value)))
    if (!is.null(basis)) {
      n <- length(cl) - 1L
      poly <- .add(lapply(seq_along(cl), function(i) {
        k <- n - i + 1L
        .mul(list(.int(gmp::divq.bigz(gmp::as.bigz(cl[[i]]$value), basis^(n - k))), .pow(w, .int(k))))
      }))
      coeff <- .mul(list(coeff, .int(basis)))
    }
  }
  list(coeff = coeff, poly = poly)
}

# _integer_basis: the largest d with p(d*y) = m*q(y) and q integral and
# smaller; cl are the coefficients, highest degree first.
.integer_basis <- function(cl) {
  n <- length(cl) - 1L
  nz <- which(vapply(cl, function(c) c != 0, TRUE))
  monoms <- n - (nz - 1L)
  coeffs <- lapply(cl[nz], abs)
  if (!(coeffs[[1L]] < coeffs[[length(coeffs)]])) return(NULL)
  coeffs <- rev(coeffs)
  m0 <- monoms[1L]
  monoms <- m0 - rev(monoms)
  monoms <- monoms[-length(monoms)]; coeffs <- coeffs[-length(coeffs)]
  if (length(monoms) == 1L) {
    rr <- .iroot(coeffs[[1L]], monoms[1L])
    return(if (rr$exact) rr$root else NULL)
  }
  g <- Reduce(gmp::gcd.bigz, coeffs)
  ds <- rev(gmp::as.bigz(as.character(divisors(Integer(g)))))
  ds <- ds[ds != 1]
  for (i in seq_along(ds)) {
    d <- ds[i]
    if (all(vapply(seq_along(monoms), function(j) gmp::mod.bigz(coeffs[[j]], d^monoms[j]) == 0, TRUE))) return(d)
  }
  NULL
}

# SymPy's Expr.is_rational_function(v)
.is_ratfun <- function(e, v) {
  if (!.has1(e, v) || identical(e$key, v$key)) return(TRUE)
  if (e$type %in% c("Add", "Mul")) return(all(vapply(e$args, .is_ratfun, TRUE, v = v)))
  if (e$type == "Pow") return(e$args[[2L]]$type == "Integer" && .is_ratfun(e$args[[1L]], v))
  FALSE
}

# RootSum._rational_case: symmetric reduction of the sum over formal roots.
.rootsum_rational <- function(poly, fexpr) {
  w <- .rootsum_w
  n <- as.integer(degree(poly, w))
  rs <- lapply(seq_len(n) - 1L, function(i) Dummy(paste0("r", i)))
  f <- .add(lapply(rs, function(r) subs(fexpr, w, r)))
  nd <- .as_numer_denom(together(f))
  vi <- viete(poly, w, rs)
  red <- function(z) {
    s <- do.call(symmetrize, c(list(expand(z)), rs, list(formal = TRUE)))
    if (!is_zero_num(s[[2L]])) stop("the sum over the roots is not symmetric", call. = FALSE)
    xreplace(s[[1L]], lapply(seq_len(n), function(k) list(s[[3L]][[k]][[1L]], vi[[k]][[2L]])))
  }
  sym_factor(.div(red(nd[[1L]]), red(nd[[2L]])))
}

.rootsum_new <- function(poly, fexpr, auto) {
  w <- .rootsum_w
  if (!.has1(fexpr, w)) return(fexpr)
  if (!.is_ratfun(fexpr, w) || !isTRUE(auto)) .rootsum_node(poly, fexpr, auto) else .rootsum_rational(poly, fexpr)
}

register_bound_vars("RootSum", function(e) .rootsum_w$key)
register_rebuilder("RootSum", function(e, args)
  RootSum(args[[1L]], Lambda(.rootsum_w, args[[2L]]$args[[2L]]), .rootsum_w, auto = isTRUE(e$extra$auto)))
register_diff_hook("RootSum", function(e, x)
  .rootsum_new(e$args[[1L]], .diff1(e$args[[2L]]$args[[2L]], x), e$extra$auto))
register_doit("RootSum", function(e) {
  w <- .rootsum_w
  rts <- roots(e$args[[1L]], w, multiple = TRUE)
  if (length(rts) < degree(e$args[[1L]], w)) return(e)
  .add(lapply(rts, function(r) subs(e$args[[2L]]$args[[2L]], w, r)))
})
register_neval("RootSum", function(e, env) {
  rts <- nroots(e$args[[1L]])
  v <- sum(vapply(rts, function(r) as.complex(.neval(subs(e$args[[2L]]$args[[2L]], .rootsum_w, r), env)), 0i))
  # drop round-off in a component that cancels, relative to the sum
  if (abs(Re(v)) < 1e-12 * Mod(v)) v <- complex(real = 0, imaginary = Im(v))
  if (abs(Im(v)) < 1e-12 * Mod(v)) v <- Re(v)
  v
})
register_str_printer("RootSum", function(e) {
  ws <- Symbol("w")
  p <- .str(xreplace(e$args[[1L]], list(list(.rootsum_w, ws))))
  fe <- e$args[[2L]]$args[[2L]]
  if (identical(fe$key, .rootsum_w$key)) return(paste0("RootSum(", p, ")"))
  paste0("RootSum(", p, ", Lambda(w, ", .str(xreplace(fe, list(list(.rootsum_w, ws)))), "))")
})
register_latex_printer("RootSum", function(e) {
  fe <- e$args[[2L]]$args[[2L]]
  args <- .latex(e$args[[1L]])
  if (!identical(fe$key, .rootsum_w$key)) args <- paste0(args, ", ", .latex(e$args[[2L]]))
  paste0("\\operatorname{RootSum} {\\left(", args, "\\right)}")
})
