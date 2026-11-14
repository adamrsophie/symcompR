# Integral.as_sum() and Integral.principal_value() of
# sympy/integrals/integrals.py.

.integral_check <- function(I) {
  I <- as_sym(I)
  if (I$type != "Integral") stop("TypeError: expected an Integral, got ", .str(I), call. = FALSE)
  I
}

#' Riemann sums and principal values of integrals
#'
#' `integral_as_sum(I, n, method)` approximates the definite integral `I`
#' by a sum over `n` subintervals (a positive integer, or by default a
#' positive integer symbol) with the `"midpoint"`, `"left"`, `"right"` or
#' `"trapezoid"` rule, as SymPy's `Integral.as_sum()`.
#' `principal_value(I)` is the Cauchy principal value of `I`, as SymPy's
#' `Integral.principal_value()`.
#'
#' @param I An unevaluated `Integral()` over one variable with limits.
#' @param n The number of subintervals.
#' @param method The rule.
#' @param evaluate If `FALSE`, return the unevaluated `Sum`.
#' @param ... Passed to `sym_integrate()`.
#' @return An expression.
#' @examples
#' x <- Symbol("x")
#' integral_as_sum(Integral(x^2, list(x, 0, 2)), 4)
#' integral_as_sum(Integral(x^2, list(x, 0, 2)), 4, method = "trapezoid")
#' principal_value(Integral(1/x, list(x, -1, 1)))
#' @export
integral_as_sum <- function(I, n = NULL, method = "midpoint", evaluate = TRUE) {
  I <- .integral_check(I)
  lims <- I$args[-1L]
  if (length(lims) > 1L) stop("NotImplementedError: Multidimensional midpoint rule not implemented yet", call. = FALSE)
  lim <- lims[[1L]]$args
  if (length(lim) != 3L || isFALSE(ask_fact(lim[[2L]], "finite")) || isFALSE(ask_fact(lim[[3L]], "finite")))
    stop("ValueError: Expecting a definite integral over a finite interval.", call. = FALSE)
  n <- if (is.null(n)) Dummy("n", integer = TRUE, positive = TRUE) else as_sym(n)
  if (isFALSE(ask_fact(n, "positive")) || isFALSE(ask_fact(n, "integer")) || isFALSE(ask_fact(n, "finite")))
    stop("ValueError: n must be a positive integer, got ", .str(n), call. = FALSE)
  x <- lim[[1L]]; a <- lim[[2L]]; b <- lim[[3L]]
  dx <- .div(.sub(b, a), n)
  k <- Dummy("k", integer = TRUE, positive = TRUE)
  f <- I$args[[1L]]
  at <- function(p) subs(f, x, p)
  kdx <- .mul(list(k, dx))
  result <- switch(method,
    left = .mul(list(dx, Sum(at(.add(list(a, .mul(list(.sub(k, S_One), dx))))), list(k, 1, n)))),
    right = .mul(list(dx, Sum(at(.add(list(a, kdx))), list(k, 1, n)))),
    midpoint = .mul(list(dx, Sum(at(.sub(.add(list(a, kdx)), .div(dx, S_Two))), list(k, 1, n)))),
    trapezoid = .mul(list(dx, .add(list(.div(.add(list(at(a), at(b))), S_Two),
                                        Sum(at(.add(list(a, kdx))), list(k, 1, .sub(n, S_One))))))),
    stop("ValueError: Unknown method ", method, call. = FALSE))
  if (evaluate) doit(result) else result
}

#' @rdname integral_as_sum
#' @export
principal_value <- function(I, ...) {
  I <- .integral_check(I)
  lims <- I$args[-1L]
  if (length(lims) != 1L || length(lims[[1L]]$args) != 3L)
    stop("ValueError: You need to insert a variable, lower_limit, and upper_limit correctly to calculate cauchy's principal value", call. = FALSE)
  x <- lims[[1L]]$args[[1L]]; a <- lims[[1L]]$args[[2L]]; b <- lims[[1L]]$args[[3L]]
  if (!(.is_comparable(a) && .is_comparable(b) && .is_true(Le(a, b))))
    stop("ValueError: The lower_limit must be smaller than or equal to the upper_limit to calculate cauchy's principal value. Also, a and b need to be comparable.", call. = FALSE)
  if (identical(a$key, b$key)) return(S_Zero)
  r <- Dummy("r")
  f <- I$args[[1L]]
  sing <- singularities(f, x)
  sing <- if (is_sym(sing) && sing$type == "FiniteSet") sing$args else list()
  sing <- Filter(function(s) .is_comparable(s) && .is_true(Le(a, s)) && .is_true(Le(s, b)), sing)
  for (s in sing) if (identical(s$key, a$key) || identical(s$key, b$key))
    stop("ValueError: The principal value is not defined in the given interval due to singularity at ", .str(s), ".", call. = FALSE)
  F <- sym_integrate(f, x, ...)
  if (.has_type(F, "Integral")) return(I)
  val <- if (is_const(a, "-oo") && is_const(b, "oo")) limit(.sub(F, subs(F, x, .neg(x))), x, S_Infinity)
         else .sub(limit(F, x, b, "-"), limit(F, x, a, "+"))
  for (s in sing) val <- .add(list(val, limit(.sub(subs(F, x, .sub(s, r)), subs(F, x, .add(list(s, r)))), r, S_Zero, "+")))
  val
}

#' Change of variables in an integral
#'
#' `integral_transform(I, x, u)` performs the substitution `x = u` in the
#' integral `I`, as SymPy's `Integral.transform()`: when `x` is a symbol,
#' `u` is an expression in a new variable (give `list(u, uvar)` when it has
#' several free symbols); when `u` is a symbol, `x` is an expression in the
#' integration variable that becomes `u`. Limits are mapped, and swapped
#' with a sign change when they come out in decreasing order.
#'
#' @param I An unevaluated `Integral()`.
#' @param x An integration variable, or an expression in one.
#' @param u The expression replacing it, or `list(expr, symbol)`.
#' @return An `Integral`.
#' @examples
#' x <- Symbol("x"); u <- Symbol("u")
#' integral_transform(Integral(x, list(x, 0, 1)), x, 2*u)
#' integral_transform(Integral(x * exp(x^2), list(x, 0, 1)), x^2, u)
#' @export
integral_transform <- function(I, x, u) {
  I <- .integral_check(I)
  x <- as_sym(x)
  d <- Dummy("d")
  lims <- lapply(I$args[-1L], function(l) l$args)
  vars <- lapply(lims, `[[`, 1L)
  vkeys <- vapply(vars, `[[`, "", "key")
  xfree <- Filter(function(s) s$key %in% vkeys, free_symbols(x))
  if (length(xfree) > 1L)
    stop("ValueError: F(x) can only contain one of: ", paste(vapply(vars, .str, ""), collapse = ", "), call. = FALSE)
  xvar <- if (length(xfree)) xfree[[1L]] else d
  if (!xvar$key %in% vkeys) return(I)
  if (is.list(u) && !is_sym(u)) {
    uvar <- as_sym(u[[2L]]); u <- as_sym(u[[1L]])
    if (!uvar$key %in% vapply(free_symbols(u), `[[`, "", "key"))
      stop("ValueError: Expecting a tuple (expr, symbol) where symbol identified a free symbol in expr, but symbol is not in expr's free symbols.", call. = FALSE)
    if (uvar$type != "Symbol") stop("ValueError: Expecting a tuple (expr, symbol) but didn't get a symbol; got ", .str(uvar), call. = FALSE)
  } else {
    u <- as_sym(u)
    uf <- free_symbols(u)
    if (!length(uf)) stop("ValueError: f(u) cannot be a constant", call. = FALSE)
    if (length(uf) > 1L)
      stop("ValueError: When f(u) has more than one free symbol, the one replacing x must be identified: pass f(u) as (f(u), u)", call. = FALSE)
    uvar <- uf[[1L]]
  }
  if (x$type == "Symbol" && u$type == "Symbol") return(xreplace(I, list(list(x, u))))
  if (x$type != "Symbol" && u$type != "Symbol") stop("ValueError: either x or u must be a symbol", call. = FALSE)
  if (identical(uvar$key, xvar$key))
    return(xreplace(integral_transform(I, x, list(subs(u, uvar, d), d)), list(list(d, uvar))))
  if (uvar$key %in% vkeys)
    stop("ValueError: u must contain the same variable as in x or a variable that is not already an integration variable", call. = FALSE)
  if (x$type != "Symbol") {
    F <- list(subs(x, xvar, d))
    soln <- solve(.sub(u, x), xvar, check = FALSE)
    if (!length(soln)) stop("ValueError: no solution for solve(F(x) - f(u), x)", call. = FALSE)
    f <- lapply(soln, function(fi) subs(fi, uvar, d))
  } else {
    f <- list(subs(u, uvar, d))
    pz <- posify(.sub(u, x))
    puvar <- xreplace(uvar, lapply(pz$reps, rev))
    soln <- lapply(solve(pz$expr, puvar), function(s) xreplace(s, pz$reps))
    if (!length(soln)) stop("ValueError: no solution for solve(F(x) - f(u), u)", call. = FALSE)
    F <- lapply(soln, function(fi) subs(fi, xvar, d))
  }
  fun <- I$args[[1L]]
  nf <- lapply(f, function(fi) subs(.mul(list(subs(fun, xvar, fi), .diff1(fi, d))), d, uvar))
  nf <- nf[!duplicated(vapply(nf, `[[`, "", "key"))]
  if (length(nf) > 1L) stop("ValueError: The mapping between F(x) and f(u) did not give a unique integrand.", call. = FALSE)
  newfunc <- nf[[1L]]
  calc1 <- function(Fi, a, b) {
    wok <- subs(Fi, d, a)
    if (is_nan(wok) || (isFALSE(ask_fact(wok, "finite")) && isTRUE(ask_fact(a, "finite"))))
      return(limit(.mul(list(fcall("sign", list(b)), Fi)), d, a))
    wok
  }
  calc <- function(a, b) {
    av <- lapply(F, calc1, a = a, b = b)
    av <- av[!duplicated(vapply(av, `[[`, "", "key"))]
    if (length(av) > 1L) stop("ValueError: The mapping between F(x) and f(u) did not give a unique limit.", call. = FALSE)
    av[[1L]]
  }
  newlims <- list()
  for (xab in lims) {
    if (identical(xab[[1L]]$key, xvar$key)) {
      if (length(xab) == 3L) {
        a <- calc(xab[[2L]], xab[[3L]]); b <- calc(xab[[3L]], xab[[2L]])
        if (.is_true(Gt(.sub(a, b), S_Zero))) { tmp <- a; a <- b; b <- tmp; newfunc <- .neg(newfunc) }
        newlims[[length(newlims) + 1L]] <- list(uvar, a, b)
      } else if (length(xab) == 2L) {
        newlims[[length(newlims) + 1L]] <- list(uvar, calc(xab[[2L]], S_One))
      } else newlims[[length(newlims) + 1L]] <- list(uvar)
    } else newlims[[length(newlims) + 1L]] <- xab
  }
  do.call(Integral, c(list(newfunc), newlims))
}
