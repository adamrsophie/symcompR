# refine(), after sympy/assumptions/refine.py, with the _eval_refine methods
# of Symbol, Pow, exp, Boolean and Predicate.

#' Simplify an expression under assumptions
#'
#' `refine(expr, assumptions)` rewrites `expr` using what [ask()] can decide
#' under the assumptions: `Abs(x)` becomes `x` when `x` is known to be
#' nonnegative, `(-1)^(x + 1)` becomes `-(-1)^x` for integer `x`, a
#' predicate or relation that the assumptions decide becomes `TRUE` or
#' `FALSE`, and so on for powers, `exp()`, `re()`, `im()`, `arg()`,
#' `sign()`, `atan2()`, `sin()`, `cos()`, `floor()`, `ceiling()`,
#' `Heaviside()` and elements of symmetric matrices. The arguments are
#' refined first.
#'
#' @param expr An expression.
#' @param assumptions Assumptions such as `Q$positive(x)`, or `TRUE` for
#'   none.
#' @param ... Further assumptions, joined to `assumptions` with `&`.
#' @return The refined expression.
#' @examples
#' x <- Symbol("x")
#' refine(sqrt(x^2), Q$real(x))
#' refine(sqrt(x^2), Q$positive(x))
#' refine(Q$real(x), Q$positive(x))
#' @export
refine <- function(expr, assumptions = TRUE, ...) {
  extra <- list(...)
  if (length(extra)) assumptions <- do.call(And, c(list(.ask_bool(assumptions)), lapply(extra, .ask_bool)))
  .refine(expr, assumptions)
}

.refine <- function(expr, assumptions) {
  if (inherits(expr, "sym_predicate_fn")) expr <- attr(expr, "predicate")
  if (!is_sym(expr)) return(expr)
  if (length(expr$args)) {
    args <- lapply(expr$args, .refine, assumptions = assumptions)
    # expr.func(*args) gives an equal object when no argument changed
    same <- all(mapply(function(a, b) is_sym(a) && is_sym(b) && .keq(a, b), args, expr$args))
    if (!same) expr <- sym_rebuild(expr, args)
  }
  if (!is_sym(expr)) return(expr)
  ref <- .eval_refine(expr, assumptions)
  if (!is.null(ref)) return(ref)
  h <- .REFINE_HANDLERS[[.refine_name(expr)]]
  if (is.null(h)) return(expr)
  new <- h(expr, assumptions)
  if (is.null(new) || (is_sym(new) && .keq(new, expr))) return(expr)
  if (!is_sym(new) || !.is_expr(new)) return(new)
  .refine(new, assumptions)
}

.refine_name <- function(e) if (e$type == "Function") e$name else e$type

.ask_true <- function(p, assumptions) isTRUE(ask(p, assumptions))

# _eval_refine
.eval_refine <- function(expr, assumptions) {
  t <- expr$type
  if (t %in% c("Symbol", "Dummy", "Wild", "Predicate", "UndefinedPredicate")) return(expr)
  if (t == "Pow") {
    be <- .refine_base_exp(expr)
    b <- be[[1L]]; e <- be[[2L]]
    if (isTRUE(ask(.QA("integer", e), assumptions)) && could_extract_minus_sign(b)) {
      if (isTRUE(ask(.QA("even", e), assumptions))) return(.pow(.neg(b), e))
      if (isTRUE(ask(.QA("odd", e), assumptions))) return(.neg(.pow(.neg(b), e)))
    }
    return(NULL)
  }
  if (.is_fn(expr, "exp")) return(.exp_eval_refine(expr))
  if (.is_boolean_kind(expr)) {
    ret <- ask(expr, assumptions)
    if (isTRUE(ret)) return(S_true)
    if (identical(ret, FALSE)) return(S_false)
  }
  NULL
}

# Pow.as_base_exp
.refine_base_exp <- function(e) {
  b <- e$args[[1L]]; x <- e$args[[2L]]
  if (b$type == "Rational") {
    q <- gmp::as.bigq(b$value)
    if (gmp::numerator(q) == 1) return(list(make_rational(gmp::as.bigq(gmp::denominator(q))), .neg(x)))
  }
  list(b, x)
}

# exp._eval_refine (which asks without the assumptions)
.exp_eval_refine <- function(expr) {
  arg <- expr$args[[1L]]
  if (arg$type != "Mul") return(NULL)
  ioo <- .mul(list(S_I, S_Infinity))
  if (.keq(arg, ioo) || .keq(arg, .neg(ioo))) return(S_NaN)
  coeff <- as_coefficient(arg, .mul(list(S_Pi, S_I)))
  if (is.null(coeff) || .py_falsy(coeff)) return(NULL)
  if (isTRUE(ask(.QA("integer", .mul(list(.int(2L), coeff)))))) {
    if (isTRUE(ask(.QA("even", coeff)))) return(S_One)
    if (isTRUE(ask(.QA("odd", coeff)))) return(S_NegativeOne)
    if (isTRUE(ask(.QA("even", .add(list(coeff, S_Half)))))) return(.neg(S_I))
    if (isTRUE(ask(.QA("odd", .add(list(coeff, S_Half)))))) return(S_I)
  }
  NULL
}

# Python's % on numbers: the remainder with the sign of the divisor
.py_mod2 <- function(a) {
  q <- gmp::as.bigq(a$value)
  r <- q - 2 * gmp::as.bigq(.refine_bigq_floor(q / 2))
  make_rational(r)
}
.refine_bigq_floor <- function(q) {
  n <- gmp::numerator(q); d <- gmp::denominator(q)
  f <- n %/% d
  if (f * d > n) f <- f - 1
  f
}

.refine_abs <- function(expr, assumptions) {
  arg <- expr$args[[1L]]
  if (.ask_true(.QA("real", arg), assumptions) && isTRUE(.fuzzy_not(ask(.QA("negative", arg), assumptions))))
    return(arg)
  if (.ask_true(.QA("negative", arg), assumptions)) return(.neg(arg))
  if (arg$type == "Mul") {
    r <- lapply(.sympy_args(arg), function(a) .refine(fcall("Abs", list(a)), assumptions))
    non_abs <- list(); in_abs <- list()
    for (i in r) {
      if (.is_fn(i, "Abs")) in_abs[[length(in_abs) + 1L]] <- i$args[[1L]]
      else non_abs[[length(non_abs) + 1L]] <- i
    }
    return(.mul(list(.mul(non_abs), fcall("Abs", list(.mul(in_abs))))))
  }
  NULL
}

.refine_pow <- function(expr, assumptions) {
  base <- expr$args[[1L]]; ex <- expr$args[[2L]]
  if (.is_fn(base, "Abs")) {
    arg <- base$args[[1L]]
    if (.ask_true(.QA("even", ex), assumptions)) {
      if (.ask_true(.QA("real", arg), assumptions)) return(.pow(arg, ex))
      if (.ask_true(.QA("imaginary", arg), assumptions))
        return(.mul(list(.pow(S_NegativeOne, .div(ex, .int(2L))), .pow(arg, ex))))
    }
  }
  if (.ask_true(.QA("real", base), assumptions)) {
    if (.is_number(base)) {
      if (.ask_true(.QA("even", ex), assumptions)) return(.pow(fcall("Abs", list(base)), ex))
      if (.ask_true(.QA("odd", ex), assumptions))
        return(.mul(list(fcall("sign", list(base)), .pow(fcall("Abs", list(base)), ex))))
    }
    if (ex$type %in% c("Integer", "Rational") && base$type == "Pow")
      return(.pow(fcall("Abs", list(base$args[[1L]])), .mul(list(base$args[[2L]], ex))))
    if (.keq(base, S_NegativeOne) && ex$type == "Add") {
      old <- expr
      args <- .sympy_args(ex)
      num <- vapply(args, function(a) a$type %in% .NUMBER_TYPES, TRUE)
      coeff <- if (any(num)) args[[which(num)[1L]]] else S_Zero
      terms <- args[!num]
      initial <- length(terms)
      even_terms <- logical(length(terms)); odd_terms <- logical(length(terms))
      for (k in seq_along(terms)) {
        if (.ask_true(.QA("even", terms[[k]]), assumptions)) even_terms[k] <- TRUE
        else if (.ask_true(.QA("odd", terms[[k]]), assumptions)) odd_terms[k] <- TRUE
      }
      new_coeff <- if (sum(odd_terms) %% 2L == 1L) .py_mod2(.add(list(coeff, S_One))) else .py_mod2(coeff)
      terms <- terms[!even_terms & !odd_terms]
      if (!.keq(new_coeff, coeff) || length(terms) < initial) {
        terms[[length(terms) + 1L]] <- new_coeff
        expr <- .pow(base, .add(terms))
      }
      e2 <- .mul(list(.int(2L), .pexp_attr(expr)))
      if (.ask_true(.QA("even", e2), assumptions) && could_extract_minus_sign(e2))
        e2 <- .mul(list(e2, expr$args[[1L]]))
      if (e2$type == "Add") {
        ea <- .sympy_args(e2)
        i <- ea[[1L]]
        p <- .add(ea[-1L])
        if (p$type == "Pow" && .keq(p$args[[1L]], S_NegativeOne)) {
          if (.ask_true(.QA("integer", p$args[[2L]]), assumptions)) {
            i <- .div(.add(list(i, S_One)), .int(2L))
            if (.ask_true(.QA("even", i), assumptions)) return(.pow(expr$args[[1L]], p$args[[2L]]))
            if (.ask_true(.QA("odd", i), assumptions)) return(.pow(expr$args[[1L]], .add(list(p$args[[2L]], S_One))))
            return(.pow(expr$args[[1L]], .add(list(p$args[[2L]], i))))
          }
        }
      }
      if (!.keq(old, expr)) return(expr)
    }
  }
  NULL
}

# expr.exp on what may no longer be a power
.pexp_attr <- function(e) {
  if (e$type != "Pow") .ask_attribute_error(e, "exp")
  e$args[[2L]]
}
.ask_attribute_error <- function(e, attr)
  stop(.sympy_condition("AttributeError", sprintf("'%s' object has no attribute '%s'", .ask_py_class(e), attr)))

.refine_exp <- function(expr, assumptions) {
  coeff <- as_coefficient(expr$args[[1L]], .mul(list(S_Pi, S_I)))
  if (is.null(coeff)) return(expr)
  terms <- if (coeff$type == "Add") .sympy_args(coeff) else list(coeff)
  int_terms <- list(); rem_terms <- list()
  for (t in terms) {
    if (.ask_true(.QA("integer", t), assumptions)) int_terms[[length(int_terms) + 1L]] <- t
    else rem_terms[[length(rem_terms) + 1L]] <- t
  }
  if (!length(int_terms)) return(expr)
  phase <- fcall("exp", list(.mul(list(S_Pi, S_I, .add(rem_terms)))))
  .mul(list(.pow(S_NegativeOne, .add(int_terms)), phase))
}

.refine_atan2 <- function(expr, assumptions) {
  y <- expr$args[[1L]]; x <- expr$args[[2L]]
  at <- function() fcall("atan", list(.div(y, x)))
  if (.ask_true(And(.QA("real", y), .QA("positive", x)), assumptions)) return(at())
  if (.ask_true(And(.QA("negative", y), .QA("negative", x)), assumptions)) return(.sub(at(), S_Pi))
  if (.ask_true(And(.QA("positive", y), .QA("negative", x)), assumptions)) return(.add(list(at(), S_Pi)))
  if (.ask_true(And(.QA("zero", y), .QA("negative", x)), assumptions)) return(S_Pi)
  if (.ask_true(And(.QA("positive", y), .QA("zero", x)), assumptions)) return(.div(S_Pi, .int(2L)))
  if (.ask_true(And(.QA("negative", y), .QA("zero", x)), assumptions)) return(.neg(.div(S_Pi, .int(2L))))
  if (.ask_true(And(.QA("zero", y), .QA("zero", x)), assumptions)) return(S_NaN)
  expr
}

.refine_reim <- function(expr, assumptions) {
  expanded <- expand(expr, complex = TRUE)
  if (!.keq(expanded, expr)) {
    refined <- .refine(expanded, assumptions)
    if (!.keq(refined, expanded)) return(refined)
  }
  NULL
}

.refine_re <- function(expr, assumptions) {
  arg <- expr$args[[1L]]
  if (.ask_true(.QA("real", arg), assumptions)) return(arg)
  if (.ask_true(.QA("imaginary", arg), assumptions)) return(S_Zero)
  .refine_reim(expr, assumptions)
}

.refine_im <- function(expr, assumptions) {
  arg <- expr$args[[1L]]
  if (.ask_true(.QA("real", arg), assumptions)) return(S_Zero)
  if (.ask_true(.QA("imaginary", arg), assumptions)) return(.mul(list(.neg(S_I), arg)))
  .refine_reim(expr, assumptions)
}

.refine_arg <- function(expr, assumptions) {
  rg <- expr$args[[1L]]
  if (.ask_true(.QA("positive", rg), assumptions)) return(S_Zero)
  if (.ask_true(.QA("negative", rg), assumptions)) return(S_Pi)
  NULL
}

# refine_sign (which asks Q.real and Q.imaginary without the assumptions)
.refine_sign <- function(expr, assumptions) {
  arg <- expr$args[[1L]]
  if (.ask_true(.QA("zero", arg), assumptions)) return(S_Zero)
  if (isTRUE(ask(.QA("real", arg)))) {
    if (.ask_true(.QA("positive", arg), assumptions)) return(S_One)
    if (.ask_true(.QA("negative", arg), assumptions)) return(S_NegativeOne)
  }
  if (isTRUE(ask(.QA("imaginary", arg)))) {
    arg_im <- .as_real_imag(arg)[[2L]]
    if (.ask_true(.QA("positive", arg_im), assumptions)) return(S_I)
    if (.ask_true(.QA("negative", arg_im), assumptions)) return(.neg(S_I))
  }
  expr
}

.refine_matrixelement <- function(expr, assumptions) {
  matrix <- expr$args[[1L]]; i <- expr$args[[2L]]; j <- expr$args[[3L]]
  if (.ask_true(.QA("symmetric", matrix), assumptions)) {
    if (could_extract_minus_sign(.sub(i, j))) return(expr)
    return(MatrixElement(matrix, j, i))
  }
  NULL
}

.refine_sin_cos <- function(expr, assumptions) {
  if (!.is_fn(expr, c("sin", "cos"))) .ask_type_error("refine_sin_cos expects a sin or cos function.")
  arg <- expr$args[[1L]]
  is_sin <- expr$name == "sin"
  if (.ask_true(.QA("infinite", arg), assumptions) && .ask_true(.QA("extended_real", arg), assumptions))
    return(AccumBounds(-1, 1))
  if (.ask_true(.QA("zero", arg), assumptions)) return(if (is_sin) S_Zero else S_One)
  half_coeffs <- list(); remaining <- list()
  terms <- if (arg$type == "Add") .sympy_args(arg) else list(arg)
  for (term in terms) {
    cp <- coeff(term, S_Pi)
    if (!.py_falsy(cp) && .ask_true(.QA("integer", .mul(list(.int(2L), cp))), assumptions))
      half_coeffs[[length(half_coeffs) + 1L]] <- .mul(list(.int(2L), cp))
    else remaining[[length(remaining) + 1L]] <- term
  }
  if (!length(half_coeffs)) return(expr)
  known <- S_Zero; unknown <- S_Zero; known_even <- TRUE
  for (cf in half_coeffs) {
    ev <- ask(.QA("even", cf), assumptions)
    if (is.na(ev)) unknown <- .add(list(unknown, cf))
    else {
      known <- .add(list(known, cf))
      known_even <- known_even == ev
    }
  }
  if (.keq(known, S_Zero)) return(expr)
  if (is_sin) {
    k <- .sub(known, S_One)
    k_even <- !known_even
  } else {
    k <- known
    k_even <- known_even
  }
  rem <- .add(list(.add(remaining), .div(.mul(list(unknown, S_Pi)), .int(2L))))
  if (k_even) {
    pw <- .pow(S_NegativeOne, .div(k, .int(2L)))
    rp <- if (pw$type == "Pow") .refine_pow(pw, assumptions) else NULL
    return(.mul(list(if (is.null(rp)) pw else rp, fcall("cos", list(rem)))))
  }
  pw <- .pow(S_NegativeOne, .div(.add(list(k, S_One)), .int(2L)))
  rp <- if (pw$type == "Pow") .refine_pow(pw, assumptions) else NULL
  .mul(list(if (is.null(rp)) pw else rp, fcall("sin", list(rem))))
}

.refine_heaviside <- function(expr, assumptions) {
  arg <- expr$args[[1L]]
  h0 <- if (length(expr$args) > 1L) expr$args[[2L]] else S_Half
  if (.ask_true(.QA("positive", arg), assumptions)) return(S_One)
  if (.ask_true(.QA("negative", arg), assumptions)) return(S_Zero)
  if (.ask_true(.QA("zero", arg), assumptions)) return(h0)
  expr
}

.refine_floor_ceiling <- function(expr, assumptions) {
  arg <- expr$args[[1L]]
  if (.ask_true(.QA("integer", arg), assumptions) || .ask_true(.QA("infinite", arg), assumptions)) return(arg)
  if (arg$type == "Add") {
    g <- list(); ng <- list()
    for (term in .sympy_args(arg)) {
      if (.ask_true(.QA("integer", term), assumptions) || .is_fn(term, c("floor", "ceiling")))
        g[[length(g) + 1L]] <- term
      else ng[[length(ng) + 1L]] <- term
    }
    return(.add(list(.add(g), fcall(expr$name, list(.add(ng))))))
  }
  expr
}

.REFINE_HANDLERS <- list(
  Abs = .refine_abs, Pow = .refine_pow, atan2 = .refine_atan2, re = .refine_re, im = .refine_im,
  arg = .refine_arg, sign = .refine_sign, MatrixElement = .refine_matrixelement,
  cos = .refine_sin_cos, sin = .refine_sin_cos, exp = .refine_exp, Heaviside = .refine_heaviside,
  floor = .refine_floor_ceiling, ceiling = .refine_floor_ceiling)
