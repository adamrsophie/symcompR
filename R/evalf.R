# Numerical evaluation: N() runs the port of sympy/core/evalf.py
# (evalf-core.R, evalf-eval.R); .neval evaluates in IEEE doubles for the
# functions that have no arbitrary precision evaluator.

.CONST_VALUES <- c(pi = pi, E = exp(1), EulerGamma = 0.57721566490153286061,
                   Catalan = 0.91596559417721901505,
                   GoldenRatio = 1.6180339887498948482,
                   TribonacciConstant = 1.8392867552141611326)


.cplx_clean <- function(v) {
  if (is.complex(v) && !is.na(v) && Im(v) == 0) return(Re(v))
  v
}

# Evaluate to an R double or complex number. env maps symbol names to values.
.neval <- function(e, env) {
  t <- e$type
  switch(t,
         Integer = , Rational = as.double(e$value),
         Float = e$value,
         Symbol = , Dummy = {
           v <- env[[e$name]]
           if (is.null(v)) stop("no numeric value for symbol ", e$name)
           if (is_sym(v)) v <- .neval(v, env)
           v
         },
         Constant = switch(EXPR = e$name,
                           I = 1i, oo = Inf, `-oo` = -Inf,
                           zoo = complex(real = Inf, imaginary = Inf),
                           nan = NaN,
                           .CONST_VALUES[[e$name]]),
         BooleanAtom = as.numeric(e$name == "true"),
         Add = {
           s <- 0
           for (a in e$args) s <- s + .neval(a, env)
           .cplx_clean(s)
         },
         Mul = {
           s <- 1
           for (a in e$args) s <- s * .neval(a, env)
           .cplx_clean(s)
         },
         Pow = {
           b <- .neval(e$args[[1L]], env)
           x <- .neval(e$args[[2L]], env)
           if (!is.complex(b) && !is.complex(x)) {
             if (b < 0 && x != round(x)) return(.cplx_clean(as.complex(b)^x))
             if (b == 0 && x < 0) return(complex(real = Inf, imaginary = Inf))
             return(b^x)
           }
           .cplx_clean(as.complex(b)^x)
         },
         Function = {
           d <- .FUNS[[e$name]]
           if (is.null(d) || is.null(d$numeric))
             stop("cannot evaluate ", e$name, "() numerically")
           vals <- lapply(e$args, .neval, env)
           .cplx_clean(do.call(d$numeric, vals))
         },
         {
           h <- .neval_hooks[[t]]
           if (is.null(h)) stop("cannot evaluate ", t, " numerically")
           h(e, env)
         })
}

#' Numerical evaluation
#'
#' `N()` (alias `evalf()`) evaluates an expression to a floating point
#' number, as SymPy's `evalf`: arbitrary precision binary arithmetic on gmp
#' integers that tracks the accuracy of the real and imaginary parts and
#' raises the working precision where terms cancel. Functions without an
#' arbitrary precision evaluator fall back to R double precision, capped at
#' 15 digits.
#'
#' @param e A `sym` expression.
#' @param digits Significant digits.
#' @param subs Named list of values for symbols.
#' @param maxn Maximum working precision in digits.
#' @param chop `TRUE` to replace tiny real or imaginary parts by exact
#'   zeros, or a number giving the size of a tiny part.
#' @param strict `TRUE` to raise an error when a part cannot be evaluated
#'   to full accuracy within `maxn`.
#' @param quad Quadrature method for integrals (`NULL` for tanh-sinh).
#' @param verbose `TRUE` to print the intermediate results.
#' @return A `sym` Float, a complex `a + b*I` with Float parts, or the
#'   expression with numeric subexpressions evaluated when free symbols
#'   remain.
#' @examples
#' N(sym_pi)
#' N(sqrt(S(2)), 5)
#' x <- Symbol("x")
#' N(sin(x), subs = list(x = 1))
#' N(parse_expr("(pi + pi*I)**2"), chop = TRUE)
#' @export
N <- function(e, digits = 15, subs = list(), maxn = 100, chop = FALSE, strict = FALSE,
              quad = NULL, verbose = FALSE) {
  if (inherits(e, "pyobj")) return(.pyobj_evalf(e, digits))
  if (inherits(e, "symmatrix")) return(.mat_map(e, function(z) N(z, digits, subs, maxn, chop, strict, quad, verbose)))
  if (is.list(e) && !inherits(e, "sym"))
    return(lapply(e, N, digits = digits, subs = subs, maxn = maxn, chop = chop, strict = strict, quad = quad, verbose = verbose))
  e <- as_sym(e)
  if (.is_bd(e)) return(.bd_evalf(e, digits = digits, subs = subs))
  if (is_inf_const(e) || is_nan(e)) return(e)
  # physical quantities, dimensions and prefixes stay symbolic
  if (.has_type(e, c("Quantity", "Dimension", "Prefix"))) return(.units_evalf(e, digits, subs))
  if (!is.null(subs) && !is.list(subs))
    stop(.ev_condition("TypeError", "subs must be given as a dictionary"))
  subs <- lapply(subs, as_sym)
  r <- .ev_N(e, digits, if (length(subs)) subs else NULL, maxn, chop, strict, quad, verbose)
  if (!inherits(r, "ev_unevaluated")) return(r)
  # a number without an arbitrary precision evaluator: R doubles
  v <- r[[1L]]
  if (length(free_symbols(v)) || .has_undef(v)) return(v)
  out <- tryCatch(.neval(v, list()), error = function(err) NULL)
  if (is.null(out)) return(v)
  if (digits > 15) warning("N(): a part of the expression has no arbitrary precision evaluator; digits capped at 15")
  .round_sym(out, min(digits, 15))
}

#' @rdname N
#' @export
evalf <- N

.round_sym <- function(v, digits) {
  if (is.complex(v)) {
    re <- signif(Re(v), digits)
    im <- signif(Im(v), digits)
    if (im == 0) return(Float(re))
    if (re == 0) return(.mul(list(Float(im), S_I)))
    return(.add(list(Float(re), .mul(list(Float(im), S_I)))))
  }
  if (is.nan(v)) return(S_NaN)
  if (is.infinite(v)) return(if (v > 0) S_Infinity else S_NegativeInfinity)
  Float(signif(v, digits))
}

.has_undef <- function(e) {
  if (.is_undefined_application(e)) return(TRUE)
  for (a in e$args) if (is_sym(a) && .has_undef(a)) return(TRUE)
  FALSE
}

#' Convert a numeric expression to an R value
#'
#' @param x A `sym` expression without free symbols.
#' @param ... Ignored.
#' @return A double or complex number.
#' @export
as.numeric.sym <- function(x, ...) {
  v <- .neval(x, list())
  if (is.complex(v)) stop("expression is complex; use as.complex()")
  v
}

#' @rdname as.numeric.sym
#' @export
as.double.sym <- as.numeric.sym

#' @rdname as.numeric.sym
#' @export
as.complex.sym <- function(x, ...) as.complex(.neval(x, list()))

#' @export
as.integer.sym <- function(x, ...) {
  if (x$type == "Integer") return(as.integer(as.character(x$value)))
  as.integer(.neval(x, list()))
}

#' @export
as.logical.sym <- function(x, ...) {
  if (x$type == "BooleanAtom") return(x$name == "true")
  stop("expression has no definite truth value: ", .str(x))
}

#' Exact value of a symbolic integer or rational as gmp big number
#'
#' @param x A `sym` Integer or Rational.
#' @return A `bigz` or `bigq` value.
#' @export
as_bigq <- function(x) {
  x <- as_sym(x)
  if (!is_exact_number(x)) stop("not an exact number: ", .str(x))
  x$value
}

#' @rdname as_bigq
#' @export
as_bigz <- function(x) {
  x <- as_sym(x)
  if (x$type != "Integer") stop("not an integer: ", .str(x))
  gmp::numerator(x$value)
}
