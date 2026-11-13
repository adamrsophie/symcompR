# Numerical evaluation in double precision (sympy/core/evalf.py works at
# arbitrary precision; here precision is limited to IEEE doubles).

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
#' number, optionally substituting values for symbols first. Evaluation is
#' in arbitrary precision binary arithmetic on gmp integers, as SymPy's
#' evalf (mpmath): the working precision grows until `digits` significant
#' digits are certain. Functions without an arbitrary precision evaluator
#' fall back to R double precision, capped at 15 digits.
#'
#' @param e A `sym` expression.
#' @param digits Significant digits.
#' @param subs Named list of values for symbols.
#' @return A `sym` Float, a complex `a + b*I` with Float parts, or the
#'   expression with numeric subexpressions evaluated when free symbols
#'   remain.
#' @examples
#' N(sym_pi)
#' N(sqrt(S(2)), 5)
#' x <- Symbol("x")
#' N(sin(x), subs = list(x = 1))
#' @export
N <- function(e, digits = 15, subs = list()) {
  if (inherits(e, "pyobj")) return(.pyobj_evalf(e, digits))
  if (inherits(e, "symmatrix")) return(.mat_map(e, function(z) N(z, digits, subs)))
  if (is.list(e) && !inherits(e, "sym")) return(lapply(e, N, digits = digits, subs = subs))
  e <- as_sym(e)
  if (.is_bd(e)) return(.bd_evalf(e, digits = digits, subs = subs))
  if (is_inf_const(e) || is_nan(e)) return(e)
  # physical quantities, dimensions and prefixes stay symbolic
  if (.has_type(e, c("Quantity", "Dimension", "Prefix"))) return(.units_evalf(e, digits, subs))
  # symbolic substitutions first; only numbers stay in the environment
  symb <- vapply(subs, function(v) is_sym(v) && length(free_symbols(v)) > 0L, TRUE)
  if (any(symb)) {
    e <- subs(e, lapply(names(subs)[symb], function(nm) list(Symbol(nm), subs[[nm]])))
    subs <- subs[!symb]
  }
  env <- lapply(subs, function(v) if (is_sym(v)) v else v)
  # arbitrary precision (evalf-mp.R) when every part has an evaluator
  fs0 <- free_symbols(e)
  # applied undefined functions stay as they are, arguments included
  if (!.has_undef(e) && (!length(fs0) || all(vapply(fs0, function(s) !is.null(env[[s$name]]), TRUE)))) {
    prec <- .dps_to_prec(digits)
    r <- tryCatch(.mpeval_top(e, prec, lapply(env, as_sym)), error = .mp_error_handler)
    if (is_sym(r)) return(r)
    if (!is.null(r)) return(.mpc_sym(r, prec, zero_int = TRUE))
  } else {
    # free symbols remain: every numeric part at the requested precision
    r <- tryCatch(.partial_evalf_mp(e, lapply(env, as_sym), .dps_to_prec(digits)), error = function(err) NULL)
    if (!is.null(r)) return(r)
  }
  if (digits > 15) warning("N(): a part of the expression has no arbitrary precision evaluator; digits capped at 15")
  digits <- min(digits, 15)
  if (length(subs)) {
    fs <- free_symbols(e)
    left <- Filter(function(s) is.null(env[[s$name]]), fs)
  } else {
    left <- free_symbols(e)
  }
  if (length(left)) return(.partial_evalf(e, env, digits))
  v <- .neval(e, env)
  .round_sym(v, digits)
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

# evalf of an expression with free symbols: numeric subexpressions become
# Floats of prec bits; integer exponents stay exact. NULL when a numeric part
# has no arbitrary precision evaluator.
.partial_evalf_mp <- function(e, env, prec) {
  fail <- FALSE
  rec <- function(z) {
    if (fail || .is_undefined_application(z)) return(z)
    fs <- free_symbols(z)
    if (!.has_undef(z) && (!length(fs) || all(vapply(fs, function(s) !is.null(env[[s$name]]), TRUE)))) {
      if (z$type == "Constant" && z$name %in% c(.INF_NAMES, "nan")) return(z)
      if (z$type == "Constant" && z$name == "I") return(z)
      r <- .mpeval_adaptive(z, prec, env)
      if (is.null(r)) { fail <<- TRUE; return(z) }
      return(.mpc_sym(r, prec, zero_int = TRUE))
    }
    if (!length(z$args)) return(z)
    # a function's own evaluator (SymPy's _eval_evalf)
    if (z$type == "Function" && !is.null(.FUNS[[z$name]]$evalf)) return(rec(.FUNS[[z$name]]$evalf(z)))
    if (z$type == "Pow" && z$args[[2L]]$type == "Integer") return(.pow(rec(z$args[[1L]]), z$args[[2L]]))
    # Mul._eval_evalf: a -1 coefficient stays exact
    if (z$type == "Mul") {
      m1 <- vapply(z$args, function(a) identical(a$key, S_NegativeOne$key), TRUE)
      if (any(m1)) return(.neg(rec(.mul(z$args[!m1]))))
    }
    sym_rebuild(z, lapply(z$args, rec))
  }
  r <- rec(e)
  if (fail) NULL else r
}

.has_undef <- function(e) {
  if (.is_undefined_application(e)) return(TRUE)
  for (a in e$args) if (is_sym(a) && .has_undef(a)) return(TRUE)
  FALSE
}

.partial_evalf <- function(e, env, digits) {
  if (.is_undefined_application(e)) return(e)
  if (!length(e$args)) {
    if (e$type %in% c("Symbol", "Dummy")) {
      v <- env[[e$name]]
      if (is.null(v)) return(e)
      return(.round_sym(.neval(e, env), digits))
    }
    if (is_exact_number(e) || e$type == "Constant") {
      if (e$type == "Constant" && e$name %in% c(.INF_NAMES, "nan", "I")) return(e)
      return(.round_sym(.neval(e, env), digits))
    }
    return(e)
  }
  if (e$type == "Pow" && e$args[[2L]]$type == "Integer") {
    return(.pow(.partial_evalf(e$args[[1L]], env, digits), e$args[[2L]]))
  }
  args <- lapply(e$args, .partial_evalf, env, digits)
  sym_rebuild(e, args)
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
