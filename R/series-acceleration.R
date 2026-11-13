# Convergence acceleration, finite differences and Pade approximants by the
# extended Euclidean algorithm, after sympy/series/acceleration.py,
# sympy/series/kauers.py and sympy/series/approximants.py.

#' Convergence acceleration
#'
#' `richardson(A, k, n, N)` estimates the limit of the sequence `A(k)` as
#' `k` goes to infinity by Richardson extrapolation of order `N` from the
#' terms `A(n), ..., A(n + N)`; `shanks(A, k, n, m)` applies the Shanks
#' transformation `m` times, as SymPy's functions of the same names.
#'
#' @param A An expression in `k`.
#' @param k The index symbol.
#' @param n Starting index.
#' @param N Order of the extrapolation.
#' @param m Number of applications.
#' @return An expression.
#' @examples
#' k <- Symbol("k")
#' richardson((1 + 1/k)^k, k, 10, 3)
#' shanks(1 - 1/k, k, 4)
#' @export
richardson <- function(A, k, n, N) {
  A <- as_sym(A); n <- as.integer(n); N <- as.integer(N)
  s <- S_Zero
  for (j in 0:N)
    s <- s + .div(.mul(list(doit(subs(A, k, .int(n + j))), .pow(.int(n + j), .int(N)), .pow(S_NegativeOne, .int(j + N)))),
                  .mul(list(fcall("factorial", list(.int(j))), fcall("factorial", list(.int(N - j))))))
  s
}

#' @rdname richardson
#' @export
shanks <- function(A, k, n, m = 1) {
  A <- as_sym(A); n <- as.integer(n); m <- as.integer(m)
  table <- lapply(0:(n + m + 1L), function(j) doit(subs(A, k, .int(j))))
  table2 <- table
  if (m >= 1L) for (i in 1:m) {
    for (j in i:(n + m)) {
      x <- table[[j]]; y <- table[[j + 1L]]; z <- table[[j + 2L]]
      table2[[j + 1L]] <- .div(.sub(.mul(list(z, x)), .pow(y, S_Two)), .sub(.add(list(z, x)), .mul(list(S_Two, y))))
    }
    table <- table2
  }
  table[[n + 1L]]
}

#' Finite differences
#'
#' `finite_diff(expression, variable, increment)` is the forward difference
#' `f(x + h) - f(x)` of the expanded expression, and `finite_diff_kauers(sum)`
#' the difference of a [Sum()] with respect to its upper limits, as in
#' `sympy.series.kauers`.
#'
#' @param expression An expression.
#' @param variable A symbol.
#' @param increment The step.
#' @param sum A `Sum`.
#' @return An expression.
#' @examples
#' x <- Symbol("x")
#' finite_diff(x^2 + 2*x + 1, x)
#' n <- Symbol("n")
#' finite_diff_kauers(Sum(x^2, list(x, 1, n)))
#' @export
finite_diff <- function(expression, variable, increment = 1) {
  if (!is_sym(expression))
    stop(.sympy_condition("AttributeError", sprintf("'%s' object has no attribute 'expand'", if (is.integer(expression) || (is.numeric(expression) && expression == round(expression))) "int" else "float")))
  expression <- expand(expression)
  expression2 <- expand(subs(expression, variable, as_sym(variable) + increment))
  expression2 - expression
}

#' @rdname finite_diff
#' @export
finite_diff_kauers <- function(sum) {
  fn <- sum$args[[1L]]
  for (l in sum$args[-1L]) fn <- subs(fn, l$args[[1L]], .add(list(l$args[[length(l$args)]], S_One)))
  fn
}

#' Pade approximants by the extended Euclidean algorithm
#'
#' `pade_approximants_gcdex(f, order)` lists the pairs `list(numerator,
#' denominator)` of the Pade approximants of the polynomial `f` (a [Poly()]
#' or an expression in one variable) truncated at degree `order`, from the
#' steps of [gcdex_steps()] on `x^(order + 1)` and `f`.
#' `pade_approximant_gcdex(f, m, n)` returns the approximant with numerator
#' degree at most `m` and denominator degree at most `n` (`n = m` by
#' default), as SymPy's functions of the same names.
#'
#' @param f A univariate polynomial.
#' @param order The truncation order.
#' @param m,n Degrees of the numerator and the denominator.
#' @return A list of `Poly` pairs, or one pair.
#' @examples
#' x <- Symbol("x")
#' f <- 1 + x + x^2/2 + x^3/6 + x^4/24
#' pade_approximant_gcdex(f, 2)
#' @export
pade_approximants_gcdex <- function(f, order) {
  if (order < 0) stop(.sympy_condition("ValueError", sprintf("order=%d must be a non-negative integer.", as.integer(order))))
  P <- .pade_poly(f)
  x <- P$gens[[1L]]
  cl <- .ex_trim(.coeff_list(as_expr(P), x))
  if (length(cl) - 1L > order) cl <- cl[(length(cl) - order):length(cl)]
  e <- .ex_to_expr(.ex_trim(cl), x)
  lapply(gcdex_steps(.pow(x, .int(order + 1L)), e, gens = x),
         function(st) list(Poly(st[[3L]], x, domain = "QQ"), Poly(st[[2L]], x, domain = "QQ")))
}

.pade_poly <- function(f) {
  P <- if (inherits(f, "Poly")) f else Poly(f)
  if (.poly_ngens(P) != 1L)
    stop(.sympy_condition("ValueError", sprintf("f=%s must be a univariate polynomial.", .pade_poly_str(P))))
  P
}

# the repr of a Poly over QQ, as in SymPy's messages
.pade_poly_str <- function(P) {
  paste0("Poly(", .poly_terms_str(P), ", ", paste(vapply(P$gens, sstr, ""), collapse = ", "), ", domain='QQ')")
}

#' @rdname pade_approximants_gcdex
#' @export
pade_approximant_gcdex <- function(f, m, n = NULL) {
  if (is.null(n)) n <- m
  if (m < 0 || n < 0)
    stop(.sympy_condition("ValueError", sprintf("m=%d and n=%d must be non-negative integers.", as.integer(m), as.integer(n))))
  P <- .pade_poly(f)
  x <- P$gens[[1L]]
  cl <- .ex_trim(.coeff_list(as_expr(P), x))
  nz <- which(!vapply(cl, is_zero_num, TRUE))
  min_degree <- if (length(nz)) length(cl) - max(nz) else 0L
  if (min_degree > m)
    stop(.sympy_condition("ValueError", sprintf(
      "polynomial f=%s has zero of order %d, which is greater than %d, the requested order of the numerator.",
      .pade_poly_str(P), min_degree, as.integer(m))))
  ap <- pade_approximants_gcdex(P, order = m + n)
  res <- ap[[1L]]
  for (nx in ap[-1L]) {
    if (degree(nx[[1L]]) < m) return(res)
    res <- nx
  }
  res
}

#' @rdname pade_approximants_gcdex
#' @param x The variable.
#' @param x0 The expansion point.
#' @export
pade_approximant <- function(f, x, x0 = 0, m = 6, n = NULL) {
  if (is.null(n)) n <- m
  if (m < 0 || n < 0)
    stop(.sympy_condition("ValueError", sprintf("m=%d and n=%d must be non-negative integers", as.integer(m), as.integer(n))))
  f <- as_sym(f); x <- as_sym(x); x0 <- as_sym(x0)
  s <- expand(removeO(series(subs(f, x, .add(list(x, x0))), x, 0, m + n + 1)))
  terms <- if (s$type == "Add") s$args else list(s)
  ks <- integer(); cs <- list()
  for (t in terms) {
    ce <- .so_as_coeff_exponent(t, x)
    if (ce[[2L]]$type != "Integer" || .has1(ce[[1L]], x))
      stop(.sympy_condition("ValueError", sprintf("f=%s does not have a laurent expansion in x around x0.", sstr(f))))
    ks <- c(ks, as.integer(num_value_dbl(ce[[2L]]))); cs[[length(cs) + 1L]] <- ce[[1L]]
  }
  h <- max(0L, -min(c(0L, ks)))
  if (n < h)
    stop(.sympy_condition("ValueError", sprintf(
      "Expression f=%s has pole of order %d, which is greater than %d, the requested order of the denominator.", sstr(f), h, as.integer(n))))
  e <- .add(c(list(S_Zero), lapply(seq_along(ks), function(i) .mul(list(cs[[i]], .pow(x, .int(ks[i] + h)))))))
  nd <- pade_approximant_gcdex(Poly(e, x), m, n - h)
  r <- .mul(list(.pow(x, .int(-h)), .div(as_expr(nd[[1L]]), as_expr(nd[[2L]]))))
  subs(r, x, .sub(x, x0))
}
