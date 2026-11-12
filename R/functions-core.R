# Function registry. Each entry may define:
#   eval(args)          special values; return NULL to stay unevaluated
#   diff(args, i)       partial derivative with respect to argument i
#   numeric             R function used by N() on double or complex values
#   facts(args, afacts) assumption facts for the call
#   latex               LaTeX name, or a function(args_latex) -> string
#   nargs               allowed argument counts

.NAN_FUNS <- c("exp", "log", "sin", "cos", "tan", "cot", "sec", "csc", "asin", "acos", "atan", "acot", "asec", "acsc",
               "sinh", "cosh", "tanh", "coth", "sech", "csch", "asinh", "acosh", "atanh", "acoth", "gamma", "loggamma",
               "erf", "erfc", "Abs", "sign", "floor", "ceiling")
.ZOO_VALUES <- c(
  lapply(stats::setNames(nm = c("sin", "cos", "tan", "cot", "sec", "csc", "sinc", "sinh", "cosh", "tanh", "coth", "sech", "csch", "exp")),
         function(n) function() S_NaN),
  lapply(stats::setNames(nm = c("asin", "acos", "asinh", "acosh", "log")), function(n) function() S_ComplexInfinity),
  list(acot = function() S_Zero, acoth = function() S_Zero,
       atanh = function() .mul(list(S_I, AccumBounds(.neg(.div(S_Pi, S_Two)), .div(S_Pi, S_Two)))))
)
fcall <- function(name, args) {
  # SymPy's elementary functions return nan at nan
  if (length(args) == 1L && args[[1L]]$type == "Constant" && args[[1L]]$name == "nan" && name %in% .NAN_FUNS) return(S_NaN)
  # the eval methods at complex infinity
  if (length(args) == 1L && args[[1L]]$type == "Constant" && args[[1L]]$name == "zoo" && name %in% names(.ZOO_VALUES))
    return(.ZOO_VALUES[[name]]())
  if (length(args) == 1L && name %in% .AB_FUNS) {
    a1 <- args[[1L]]
    if (a1$type == "AccumulationBounds" || (a1$type == "Constant" && a1$name %in% c("oo", "-oo") && name != "exp" && name != "log")) {
      r <- .ab_fcall(name, a1)
      if (!is.null(r)) return(r)
    }
  }
  d <- .FUNS[[name]]
  if (!is.null(d)) {
    if (!is.null(d$nargs) && !(length(args) %in% d$nargs))
      stop(name, "() takes ", paste(d$nargs, collapse = " or "),
           " argument(s), got ", length(args), call. = FALSE)
    if (!is.null(d$eval)) {
      r <- d$eval(args)
      if (!is.null(r)) return(r)
    }
    # Function.__new__: arguments that are all Floats or complex Floats
    # evaluate the call numerically
    if (length(args) && all(vapply(args, .should_evalf, TRUE)) && !all(vapply(args, function(a) a$type == "Float", TRUE))) {
      r <- .mp_float_call(name, args, .floats_prec(lapply(args, function(a) if (a$type == "Add") a$args[[1L]] else a)))
      if (!is.null(r)) return(r)
    }
  }
  new_sym("Function", args, name = name)
}

# Build an unevaluated function node without running eval.
fnode <- function(name, args) new_sym("Function", args, name = name)

#' Undefined functions
#'
#' `Function("f")` returns an R function that builds unevaluated applications
#' such as `f(x)`. Derivatives of these stay as [Derivative()] objects, which
#' is how SymPy represents unknown functions in differential equations.
#'
#' @param name Function name.
#' @return An R function of any number of arguments.
#' @examples
#' f <- Function("f")
#' x <- Symbol("x")
#' f(x) + f(x)
#' @export
Function <- function(name) {
  force(name)
  fn <- function(...) fcall(name, lapply(list(...), as_sym))
  structure(fn, class = c("sym_function", "function"), sym_name = name)
}

#' @export
print.sym_function <- function(x, ...) {
  cat(attr(x, "sym_name"), "\n")
  invisible(x)
}

# ---------------------------------------------------------------------------
# Small arithmetic helpers used throughout the package.

.neg <- function(e) .mul(list(S_NegativeOne, e))
.sub <- function(a, b) .add(list(a, .neg(b)))
.div <- function(a, b) {
  if (a$type == "AccumulationBounds" || b$type == "AccumulationBounds") {
    r <- .ab_div(a, b)
    if (!is.null(r)) return(r)
  }
  # Expr.__truediv__: 1/b is the power itself, unflattened
  if (a$type == "Integer" && a$value == 1) return(.pow(b, S_NegativeOne))
  .mul(list(a, .pow(b, S_NegativeOne)))
}
.inv <- function(a) .pow(a, S_NegativeOne)
.int <- function(n) make_rational(gmp::as.bigq(if (is.character(n)) .str_bigz(n) else gmp::as.bigz(n)))
.rat <- function(p, q) make_rational(gmp::as.bigq(p, q))
.sqrt <- function(e) .pow(e, S_Half)
.sq <- function(e) .pow(e, S_Two)

.key_lt <- function(a, b) order(c(a, b), method = "radix")[1L] == 1L &&
  a != b

#' Whether a leading minus sign can be pulled out of an expression
#'
#' Mirrors SymPy's `could_extract_minus_sign()`, which decides canonical forms
#' such as `sin(-x) = -sin(x)`.
#'
#' @param e A `sym` expression.
#' @return `TRUE` or `FALSE`.
#' @export
could_extract_minus_sign <- function(e) {
  e <- as_sym(e)
  t <- e$type
  if (t == "Integer" || t == "Rational" || t == "Float") return(num_sign(e) < 0)
  if (t == "Constant") return(e$name == "-oo")
  if (t == "Mul") {
    a1 <- e$args[[1L]]
    return((is_number_type(a1) && num_sign(a1) < 0) || is_const(a1, "-oo"))
  }
  if (t == "Add") {
    neg <- sum(vapply(e$args, could_extract_minus_sign, TRUE))
    pos <- length(e$args) - neg
    if (neg != pos) return(neg > pos)
    # SymPy: prefer the sign whose default_sort_key is larger
    return(.key_cmp(.sort_key(e), .sort_key(.neg(e))) < 0)
  }
  FALSE
}

# Rational coefficient r of r*pi, or NULL.
.pi_coeff <- function(a) {
  if (is_const(a, "pi")) return(S_One)
  if (is_zero_num(a)) return(S_Zero)
  if (a$type == "Mul" && length(a$args) == 2L && is_exact_number(a$args[[1L]]) &&
      is_const(a$args[[2L]], "pi")) return(a$args[[1L]])
  NULL
}

# Coefficient of pi in one term (Expr.coeff(pi)), or NULL.
.coeff_of_pi <- function(t) {
  if (is_const(t, "pi")) return(S_One)
  if (t$type == "Mul") {
    w <- which(vapply(t$args, function(f) is_const(f, "pi"), TRUE))
    if (length(w) == 1L) return(.mul(t$args[-w]))
  }
  NULL
}

# SymPy's _pi_coeff: for pi or a Mul containing pi, the coefficient of pi
# with its numeric factor reduced mod 2 when the rest is an integer; NULL
# otherwise.
.pi_coeff_sym <- function(a) {
  if (is_const(a, "pi")) return(S_One)
  if (is_zero_num(a)) return(S_Zero)
  if (a$type != "Mul") return(NULL)
  cx <- .coeff_of_pi(a)
  if (is.null(cx)) return(NULL)
  sp <- as_coeff_mul_split(cx)
  c <- sp[[1L]]; x <- sp[[2L]]
  if (c$type == "Float") {
    # exact binary fractions become Rationals
    v <- c$value; f <- abs(v) %% 1
    if (f != 0) {
      m <- 2^(-round(log2(f))); cm <- v * m
      if (trunc(cm) == cm) { c <- Rational(trunc(cm), m); cx <- .mul(list(c, x)) }
    } else { c <- Integer(trunc(v)); cx <- .mul(list(c, x)) }
  }
  if (isTRUE(ask_fact(x, "integer")) && is_exact_number(c)) {
    c2 <- .mod2(gmp::as.bigq(c$value))
    if (c2 == 1) return(x)
    if (c2 == 0) return(if (!is.na(ask_fact(x, "even"))) S_Zero else S_Two)
    return(.mul(list(make_rational(c2), x)))
  }
  cx
}

# SymPy's _peeloff_pi for an Add: list(rest + m1*pi, m2) where the pi
# coefficient (the rational coefficients of pi, symbolic ones included) is
# m1 + m2 with m1 = Mod(coefficient, 1/2); NULL when m2 is zero.
.peeloff_pi <- function(a) {
  ks <- list(); rest <- list()
  for (t in a$args) {
    k <- .coeff_of_pi(t)
    if (!is.null(k) && !is_zero_num(k) && isTRUE(ask_fact(k, "rational"))) ks[[length(ks) + 1L]] <- k
    else rest[[length(rest) + 1L]] <- t
  }
  if (!length(ks)) return(NULL)
  pc <- .add(ks)
  if (is_zero_num(pc)) return(NULL)
  m1 <- Mod(pc, S_Half)
  m2 <- .sub(pc, m1)
  ok <- isTRUE(ask_fact(m2, "integer")) ||
    (isTRUE(ask_fact(.mul(list(S_Two, m2)), "integer")) && identical(ask_fact(m2, "even"), FALSE))
  if (!ok || is_zero_num(m2)) return(NULL)
  list(.add(c(rest, list(.mul(list(m1, S_Pi))))), m2)
}

# If a = c*I with c not containing I, return c; else NULL.
.imag_coeff <- function(a) {
  if (is_const(a, "I")) return(S_One)
  if (a$type == "Mul") {
    w <- which(vapply(a$args, function(f) is_const(f, "I"), TRUE))
    if (length(w) == 1L) return(.mul(a$args[-w]))
  }
  NULL
}

.is_float_arg <- function(a) a$type == "Float"

.has_float <- function(e) {
  if (e$type == "Float") return(TRUE)
  if (!length(e$args)) return(FALSE)
  any(vapply(e$args, .has_float, TRUE))
}

# Numeric evaluation of f(args) when every argument is a Float, or a number
# mixed with Floats.
.float_eval <- function(name, args) {
  if (!any(vapply(args, function(a) a$type == "Float", TRUE))) return(NULL)
  if (!all(vapply(args, function(a) length(free_symbols(a)) == 0L, TRUE)))
    return(NULL)
  # Floats above double precision: evaluate at their precision
  p <- .floats_prec(args)
  if (p != 53 || any(vapply(args, function(a) a$type == "Float" && !is.null(a$extra$mp), TRUE))) {
    r <- .mp_float_call(name, args, p)
    if (!is.null(r)) return(r)
  }
  v <- tryCatch(.neval(fnode(name, args), list()), error = function(e) NULL)
  if (is.null(v) || is.na(v)) return(NULL)
  # base R's complex functions lose digits (acos(3.6) gets a real part
  # 2e-8): a complex value comes from the mp evaluator, as mpmath's
  if (is.complex(v) && Im(v) != 0) {
    r <- .mp_float_call(name, args, p)
    if (!is.null(r)) return(r)
  }
  .complex_to_sym(v)
}

.mp_float_call <- function(name, args, p) {
  r <- tryCatch(.mpeval(fnode(name, args), p + 10, list()), error = function(e) NULL)
  if (is.null(r)) return(NULL)
  .mpc_sym(.mpc(.mp_round(r$re$m, r$re$e, p), .mp_round(r$im$m, r$im$e, p)), p)
}

# Function._should_evalf: a Float, or a + b*I with number parts one of
# which is a Float
.should_evalf <- function(a) {
  if (a$type == "Float") return(TRUE)
  if (a$type != "Add" || length(a$args) != 2L) return(FALSE)
  re <- a$args[[1L]]; im <- a$args[[2L]]
  if (!is_number_type(re) || im$type != "Mul" || length(im$args) != 2L ||
      !is_number_type(im$args[[1L]]) || !is_const(im$args[[2L]], "I")) return(FALSE)
  re$type == "Float" || im$args[[1L]]$type == "Float"
}

.complex_to_sym <- function(v) {
  if (is.complex(v)) {
    if (Im(v) == 0) return(Float(Re(v)))
    return(.add(list(Float(Re(v)), .mul(list(Float(Im(v)), S_I)))))
  }
  if (is.nan(v)) return(S_NaN)
  if (is.infinite(v)) return(if (v > 0) S_Infinity else S_NegativeInfinity)
  Float(v)
}

.odd_eval <- function(name, a) {
  if (could_extract_minus_sign(a)) return(.neg(fcall(name, list(.neg(a)))))
  NULL
}

.even_eval <- function(name, a) {
  if (could_extract_minus_sign(a)) return(fcall(name, list(.neg(a))))
  NULL
}
