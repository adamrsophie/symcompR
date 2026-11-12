# Unevaluated integral transforms and Laplace helpers, after
# sympy/integrals/transforms.py, laplace.py and singularityfunctions.py.

.TRANSFORMS <- list(
  LaplaceTransform = list(fn = function(f, x, y, ...) laplace_transform(f, x, y, noconds = TRUE), tex = "\\mathcal{L}"),
  InverseLaplaceTransform = list(fn = function(f, x, y, plane = NULL) {
    if (!is.null(plane) && plane$type == "Dummy" && identical(plane$name, "None")) plane <- NULL
    inverse_laplace_transform(f, x, y, plane)
  }, tex = "\\mathcal{L}^{-1}"),
  FourierTransform = list(fn = function(f, x, y, ...) fourier_transform(f, x, y), tex = "\\mathcal{F}"),
  InverseFourierTransform = list(fn = function(f, x, y, ...) inverse_fourier_transform(f, x, y), tex = "\\mathcal{F}^{-1}"),
  MellinTransform = list(fn = function(f, x, y, ...) {
    r <- mellin_transform(f, x, y)
    if (is_sym(r)) r else Tuple(r[[1L]], Tuple(r[[2L]][[1L]], r[[2L]][[2L]]), r[[3L]])
  }, tex = "\\mathcal{M}"),
  InverseMellinTransform = list(fn = function(f, x, y, a, b) {
    nn <- function(z) if (z$type == "Dummy" && identical(z$name, "None")) NULL else z
    inverse_mellin_transform(f, x, y, list(nn(a), nn(b)))
  }, tex = "\\mathcal{M}^{-1}"),
  SineTransform = list(fn = function(f, x, y, ...) sine_transform(f, x, y), tex = "\\mathcal{SIN}"),
  InverseSineTransform = list(fn = function(f, x, y, ...) inverse_sine_transform(f, x, y), tex = "\\mathcal{SIN}^{-1}"),
  CosineTransform = list(fn = function(f, x, y, ...) cosine_transform(f, x, y), tex = "\\mathcal{COS}"),
  InverseCosineTransform = list(fn = function(f, x, y, ...) inverse_cosine_transform(f, x, y), tex = "\\mathcal{COS}^{-1}"),
  HankelTransform = list(fn = function(f, x, y, nu) hankel_transform(f, x, y, nu), tex = "\\mathcal{H}"),
  InverseHankelTransform = list(fn = function(f, x, y, nu) inverse_hankel_transform(f, x, y, nu), tex = "\\mathcal{H}^{-1}"))

.make_transform <- function(type, ...) new_sym(type, lapply(list(...), as_sym))

#' Unevaluated integral transforms
#'
#' `LaplaceTransform(f, t, s)`, `FourierTransform(f, x, k)`,
#' `MellinTransform(f, x, s)`, `SineTransform()`, `CosineTransform()`,
#' `HankelTransform(f, r, k, nu)` and their inverses are unevaluated
#' transform objects; [doit()] computes them with [laplace_transform()],
#' [fourier_transform()] and the other transform functions.
#' `InverseLaplaceTransform()` takes an optional `plane` (ignored), and
#' `InverseMellinTransform()` the strip `a`, `b`.
#'
#' @param f The function.
#' @param t,s,x,k,r Variables: the original first, the transformed second.
#' @param nu The Bessel order.
#' @param plane Ignored.
#' @param a,b The fundamental strip.
#' @return A transform object.
#' @examples
#' t <- Symbol("t", positive = TRUE); s <- Symbol("s")
#' L <- LaplaceTransform(t^2, t, s)
#' doit(L)
#' latex(L)
#' @name integral-transforms
NULL

#' @rdname integral-transforms
#' @export
LaplaceTransform <- function(f, t, s) .make_transform("LaplaceTransform", f, t, s)
#' @rdname integral-transforms
#' @export
InverseLaplaceTransform <- function(f, s, t, plane = NULL)
  .make_transform("InverseLaplaceTransform", f, s, t, if (is.null(plane)) .ilt_none() else plane)
# InverseLaplaceTransform._none_sentinel
.ilt_none <- function() { if (is.null(.mi_state$ilt_none)) .mi_state$ilt_none <- Dummy("None"); .mi_state$ilt_none }
#' @rdname integral-transforms
#' @export
FourierTransform <- function(f, x, k) .make_transform("FourierTransform", f, x, k)
#' @rdname integral-transforms
#' @export
InverseFourierTransform <- function(f, k, x) .make_transform("InverseFourierTransform", f, k, x)
#' @rdname integral-transforms
#' @export
MellinTransform <- function(f, x, s) .make_transform("MellinTransform", f, x, s)
#' @rdname integral-transforms
#' @export
InverseMellinTransform <- function(f, s, x, a, b) .make_transform("InverseMellinTransform", f, s, x, a, b)
#' @rdname integral-transforms
#' @export
SineTransform <- function(f, x, k) .make_transform("SineTransform", f, x, k)
#' @rdname integral-transforms
#' @export
InverseSineTransform <- function(f, k, x) .make_transform("InverseSineTransform", f, k, x)
#' @rdname integral-transforms
#' @export
CosineTransform <- function(f, x, k) .make_transform("CosineTransform", f, x, k)
#' @rdname integral-transforms
#' @export
InverseCosineTransform <- function(f, k, x) .make_transform("InverseCosineTransform", f, k, x)
#' @rdname integral-transforms
#' @export
HankelTransform <- function(f, r, k, nu) .make_transform("HankelTransform", f, r, k, nu)
#' @rdname integral-transforms
#' @export
InverseHankelTransform <- function(f, k, r, nu) .make_transform("InverseHankelTransform", f, k, r, nu)

for (.t in names(.TRANSFORMS)) {
  register_doit(.t, local({
    tt <- .t
    function(e) do.call(.TRANSFORMS[[tt]]$fn, e$args)
  }))
  register_str_printer(.t, local({
    tt <- .t
    function(e) paste0(tt, "(", paste(vapply(e$args, .str, ""), collapse = ", "), ")")
  }))
  register_latex_printer(.t, local({
    tt <- .t
    function(e) paste0(.TRANSFORMS[[tt]]$tex, "_{", .latex(e$args[[2L]]), "}\\left[", .latex(e$args[[1L]]),
                       "\\right]\\left(", .latex(e$args[[3L]]), "\\right)")
  }))
}
rm(.t)

# Laplace transform of an undefined function, its derivatives, sums and
# constant multiples; NULL when another method is needed.
.laplace_rules <- function(f, t, s) {
  if (!.has1(f, t)) return(.div(f, s))
  if (f$type == "Add") {
    parts <- lapply(f$args, .laplace_rules, t = t, s = s)
    if (any(vapply(parts, is.null, TRUE))) return(NULL)
    return(.add(parts))
  }
  if (f$type == "Mul") {
    dep <- vapply(f$args, function(a) .has1(a, t), TRUE)
    if (any(!dep) && sum(dep) >= 1L) {
      r <- .laplace_rules(.mul(f$args[dep]), t, s)
      return(if (is.null(r)) NULL else .mul(c(f$args[!dep], list(r))))
    }
    return(NULL)
  }
  if (.is_undefined_application(f) && length(f$args) == 1L && identical(f$args[[1L]]$key, t$key))
    return(LaplaceTransform(f, t, s))
  if (f$type == "Derivative" && .is_undefined_application(f$args[[1L]]) && length(f$args) == 2L &&
      identical(f$args[[2L]]$args[[1L]]$key, t$key) && length(f$args[[1L]]$args) == 1L) {
    g <- f$args[[1L]]
    n <- .int_val(f$args[[2L]]$args[[2L]])
    out <- list(.mul(list(.pow(s, .int(n)), LaplaceTransform(g, t, s))))
    for (k in 0:(n - 1L)) {
      init <- if (k == 0L) subs(g, t, S_Zero) else Subs(Derivative(g, list(t, k)), t, S_Zero)
      out[[length(out) + 1L]] <- .neg(.mul(list(.pow(s, .int(n - 1L - k)), init)))
    }
    return(.add(out))
  }
  NULL
}

#' Laplace helpers and singularity function integration
#'
#' `laplace_correspondence(f, fdict)` replaces
#' `LaplaceTransform(f(t), t, s)` by `F(s)` (and the inverse transform back)
#' for the pairs `list(list(f, F))`. `laplace_initial_conds(f, t, fdict)`
#' substitutes initial values `list(list(f, c(f0, f1, ...)))` for `f(0)`
#' and the derivatives of `f` at 0. `singularityintegrate(f, x)` integrates
#' singularity functions, or returns `NULL`.
#'
#' @param f An expression.
#' @param t,x Variables.
#' @param fdict Pairs of functions, or of a function and its initial values.
#' @return An expression or `NULL`.
#' @examples
#' x <- Symbol("x"); a <- Symbol("a")
#' singularityintegrate(SingularityFunction(x, a, 3), x)
#' @name transform-helpers
NULL

#' @rdname transform-helpers
#' @export
laplace_correspondence <- function(f, fdict) {
  f <- as_sym(f)
  .map_bottom_up(f, function(z) {
    for (pr in fdict) {
      lo <- as_sym(pr[[1L]](Symbol("x"))); up <- as_sym(pr[[2L]](Symbol("x")))
      # the function must be applied to the transform variable itself
      plain <- function(g, nm) g$type == "Function" && g$name == nm && length(g$args) == 1L && identical(g$args[[1L]]$key, z$args[[2L]]$key)
      if (z$type == "LaplaceTransform" && plain(z$args[[1L]], lo$name)) return(pr[[2L]](z$args[[3L]]))
      if (z$type == "InverseLaplaceTransform" && plain(z$args[[1L]], up$name)) return(pr[[1L]](z$args[[3L]]))
    }
    z
  })
}

#' @rdname transform-helpers
#' @export
laplace_initial_conds <- function(f, t, fdict) {
  f <- as_sym(f); t <- as_sym(t)
  for (pr in fdict) {
    g <- pr[[1L]](t)
    ics <- lapply(pr[[2L]], as_sym)
    for (k in seq_along(ics) - 1L) {
      target <- if (k == 0L) subs(g, t, S_Zero) else Subs(Derivative(g, list(t, k)), t, S_Zero)
      f <- xreplace(f, list(list(target, ics[[k + 1L]])))
    }
  }
  f
}

#' @rdname transform-helpers
#' @export
singularityintegrate <- function(f, x) {
  f <- as_sym(f); x <- as_sym(x)
  isf <- function(z) z$type == "Function" && z$name == "SingularityFunction"
  if (!.has_function(f, "SingularityFunction")) return(NULL)
  SF <- function(a, n) fcall("SingularityFunction", list(x, a, as_sym(n)))
  if (isf(f)) {
    a <- f$args[[2L]]; n <- f$args[[3L]]
    arg <- f$args[[1L]]
    if (!identical(simplify(.diff1(arg, x))$key, S_One$key)) return(NULL)
    SF <- function(a, n) fcall("SingularityFunction", list(arg, a, as_sym(n)))
    if (isTRUE(.facts(n)[["nonnegative"]])) return(.div(SF(a, .add(list(n, S_One))), .add(list(n, S_One))))
    if (n$type == "Integer" && .int_val(n) %in% c(-1, -2)) return(SF(a, .add(list(n, S_One))))
    return(NULL)
  }
  if (f$type == "Mul") {
    s <- Filter(isf, f$args)
    if (length(s) != 1L) return(NULL)
    g <- .mul(Filter(function(z) !isf(z), f$args))
    a <- s[[1L]]$args[[2L]]; n <- s[[1L]]$args[[3L]]
    if (!(n$type == "Integer") || !identical(s[[1L]]$args[[1L]]$key, x$key)) return(NULL)
    ga <- subs(g, x, a)
    if (.int_val(n) == -1L) return(.mul(list(ga, SF(a, 0))))
    if (.int_val(n) == -2L) return(.sub(.mul(list(ga, SF(a, -1))), .mul(list(subs(.diff1(g, x), x, a), SF(a, 0)))))
  }
  NULL
}

