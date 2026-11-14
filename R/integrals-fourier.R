# Fourier, sine, cosine and Hankel transforms, after
# sympy/integrals/transforms.py: IntegralTransform.doit with
# _fourier_transform, _sine_cosine_transform and _hankel_transform.

# IntegralTransform.doit: compute directly, else transform the terms of a
# sum, else pull out the factors free of the variable
.it_doit <- function(type, f, x, k, rest, compute, simplify = TRUE, needeval = FALSE) {
  T <- NULL
  if (!.applied_undef_in(f, x)) T <- tryCatch(compute(f, x, k, simplify), IntegralTransformError = function(e) NULL)
  fn <- if (f$type == "Add") f else expand_mul(f)
  if (!is.null(T)) return(T)
  if (fn$type == "Add") {
    res <- lapply(.sympy_args(fn), function(a) .it_doit(type, a, x, k, rest, compute, simplify, needeval))
    ress <- lapply(res, function(r) if (is_sym(r)) r else r[[1L]])
    extra <- Filter(Negate(is.null), lapply(res, function(r) if (is_sym(r) || length(r) < 2L) NULL else r[[2L]]))
    tot <- if (isTRUE(simplify)) simplify(.add(ress)) else .add(ress)
    if (!length(extra)) return(tot)
    cond <- do.call(And, extra)
    if (!identical(cond$key, S_false$key)) return(list(tot, cond))
  }
  if (needeval) .it_error(.it_name(type), "needeval")
  cm <- .as_coeff_mul_dep(fn, x)
  .mul(list(cm[[1L]], do.call(.make_transform, c(list(type, .mul(cm[[2L]]), x, k), rest))))
}
.it_name <- function(type) switch(type, FourierTransform = "Fourier", InverseFourierTransform = "Inverse Fourier",
  SineTransform = "Sine", InverseSineTransform = "Inverse Sine", CosineTransform = "Cosine",
  InverseCosineTransform = "Inverse Cosine", HankelTransform = "Hankel", InverseHankelTransform = "Inverse Hankel", type)

# the definite integral's first Piecewise branch, as the transforms read it
.it_from_piecewise <- function(F, f, name, simplify) {
  if (F$type != "Piecewise") .it_error(name, "could not compute integral")
  b <- F$args[[1L]]
  Fv <- b$args[[1L]]; cond <- b$args[[2L]]
  if (.has_type(Fv, "Integral")) .it_error(name, "integral in unexpected form")
  list(.tr_simplify(Fv, simplify), cond)
}

.fourier_core <- function(a, b, name) function(f, x, k, simplify) {
  F <- integrate(.mul(list(a, f, fcall("exp", list(.mul(list(b, S_I, x, k)))))), list(x, S_NegativeInfinity, S_Infinity))
  if (!.has_type(F, "Integral")) return(list(.tr_simplify(F, simplify), S_true))
  intf <- integrate(f, list(x, S_NegativeInfinity, S_Infinity))
  if (any(vapply(list(S_NegativeInfinity, S_Infinity, S_NaN), function(z) identical(intf$key, z$key), NA)) || .has_type(intf, "Integral"))
    .it_error(name, "function not integrable on real axis")
  .it_from_piecewise(F, f, name, simplify)
}

.sine_cosine_core <- function(kern, name) function(f, x, k, simplify) {
  a <- .div(.sqrt(S_Two), .sqrt(S_Pi))
  F <- integrate(.mul(list(a, f, fcall(kern, list(.mul(list(x, k)))))), list(x, S_Zero, S_Infinity))
  if (!.has_type(F, "Integral")) return(list(.tr_simplify(F, simplify), S_true))
  .it_from_piecewise(F, f, name, simplify)
}

.hankel_core <- function(nu, name) function(f, r, k, simplify) {
  F <- integrate(.mul(list(f, fcall("besselj", list(nu, .mul(list(k, r)))), r)), list(r, S_Zero, S_Infinity))
  if (!.has_type(F, "Integral")) return(list(.tr_simplify(F, simplify), S_true))
  .it_from_piecewise(F, f, name, simplify)
}

.it_result <- function(r, noconds) if (is_sym(r) || !noconds) r else r[[1L]]

#' Fourier, sine, cosine and Hankel transforms
#'
#' Ports of SymPy's `fourier_transform`, `inverse_fourier_transform`,
#' `sine_transform`, `inverse_sine_transform`, `cosine_transform`,
#' `inverse_cosine_transform`, `hankel_transform` and
#' `inverse_hankel_transform` (sympy/integrals/transforms.py), with SymPy's
#' conventions: the Fourier kernel is `exp(-2*pi*I*x*k)`, the sine and cosine
#' transforms carry `sqrt(2/pi)`, and the Hankel transform integrates
#' `f*besselj(nu, k*r)*r` over `(0, oo)`. When the integral is not found the
#' unevaluated transform is returned.
#'
#' @param f,F Expressions.
#' @param x,k,r Variables.
#' @param nu Order of the Hankel transform.
#' @param simplify Simplify the result.
#' @param noconds Drop the convergence condition.
#' @param needeval Signal an error instead of returning an unevaluated
#'   transform.
#' @return The transform, or `list(F, cond)` when `noconds = FALSE`.
#' @examples
#' x <- Symbol("x"); k <- Symbol("k")
#' \donttest{
#' fourier_transform(exp(-x^2), x, k)
#' }
#' @export
fourier_transform <- function(f, x, k, simplify = TRUE, noconds = TRUE, needeval = FALSE) {
  f <- as_sym(f); x <- as_sym(x); k <- as_sym(k)
  .it_result(.it_doit("FourierTransform", f, x, k, list(), .fourier_core(S_One, .mul(list(.int(-2), S_Pi)), "Fourier"), simplify, needeval), noconds)
}
#' @rdname fourier_transform
#' @export
inverse_fourier_transform <- function(F, k, x, simplify = TRUE, noconds = TRUE, needeval = FALSE) {
  F <- as_sym(F); x <- as_sym(x); k <- as_sym(k)
  .it_result(.it_doit("InverseFourierTransform", F, k, x, list(), .fourier_core(S_One, .mul(list(S_Two, S_Pi)), "Inverse Fourier"), simplify, needeval), noconds)
}
#' @rdname fourier_transform
#' @export
sine_transform <- function(f, x, k, simplify = TRUE, noconds = TRUE, needeval = FALSE) {
  f <- as_sym(f); x <- as_sym(x); k <- as_sym(k)
  .it_result(.it_doit("SineTransform", f, x, k, list(), .sine_cosine_core("sin", "Sine"), simplify, needeval), noconds)
}
#' @rdname fourier_transform
#' @export
inverse_sine_transform <- function(F, k, x, simplify = TRUE, noconds = TRUE, needeval = FALSE) {
  F <- as_sym(F); x <- as_sym(x); k <- as_sym(k)
  .it_result(.it_doit("InverseSineTransform", F, k, x, list(), .sine_cosine_core("sin", "Inverse Sine"), simplify, needeval), noconds)
}
#' @rdname fourier_transform
#' @export
cosine_transform <- function(f, x, k, simplify = TRUE, noconds = TRUE, needeval = FALSE) {
  f <- as_sym(f); x <- as_sym(x); k <- as_sym(k)
  .it_result(.it_doit("CosineTransform", f, x, k, list(), .sine_cosine_core("cos", "Cosine"), simplify, needeval), noconds)
}
#' @rdname fourier_transform
#' @export
inverse_cosine_transform <- function(F, k, x, simplify = TRUE, noconds = TRUE, needeval = FALSE) {
  F <- as_sym(F); x <- as_sym(x); k <- as_sym(k)
  .it_result(.it_doit("InverseCosineTransform", F, k, x, list(), .sine_cosine_core("cos", "Inverse Cosine"), simplify, needeval), noconds)
}
# HankelTypeTransform.doit computes directly (no linearity step)
#' @rdname fourier_transform
#' @export
hankel_transform <- function(f, r, k, nu, simplify = TRUE, noconds = TRUE) {
  f <- as_sym(f); r <- as_sym(r); k <- as_sym(k); nu <- as_sym(nu)
  .it_result(.hankel_core(nu, "Hankel")(f, r, k, simplify), noconds)
}
#' @rdname fourier_transform
#' @export
inverse_hankel_transform <- function(F, k, r, nu, simplify = TRUE, noconds = TRUE) {
  F <- as_sym(F); r <- as_sym(r); k <- as_sym(k); nu <- as_sym(nu)
  .it_result(.hankel_core(nu, "Inverse Hankel")(F, k, r, simplify), noconds)
}
