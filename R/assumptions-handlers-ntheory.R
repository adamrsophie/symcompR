# Handlers of sympy/assumptions/handlers/ntheory.py: prime, composite, even
# and odd.

.ask_has_float <- function(e) is_sym(e) && (e$type == "Float" || any(vapply(e$args, .ask_has_float, TRUE)))

.prime_number <- function(expr, assumptions) {
  exact <- !.ask_has_float(expr)
  i <- .py_round_int(expr)
  if (is.null(i)) return(FALSE)
  if (identical(.num_equals0(.sub(expr, make_rational(gmp::as.bigq(i)))), FALSE)) return(FALSE)
  if (exact) return(isTRUE(isprime(i)))
  NA
}

.even_number <- function(expr, assumptions) {
  if (expr$type == "Float") {
    q <- .mp_to_q(.float_mp(expr))
    if (gmp::denominator(q) == 1) return(NA)
    return(FALSE)
  }
  i <- .py_round_int(expr)
  if (is.null(i)) return(FALSE)
  if (!isTRUE(.num_equals0(.sub(expr, make_rational(gmp::as.bigq(i)))))) return(FALSE)
  gmp::mod.bigz(i, 2) == 0
}

# ---------------------------------------------------------------------------
# prime

.ask_register("prime", "Expr", .old_fact_handler("prime"))

.ask_register("prime", "Basic", function(expr, assumptions) {
  if (.is_number(expr)) return(.prime_number(expr, assumptions))
  NA
})

.ask_register("prime", "Mul", function(expr, assumptions) {
  if (.is_number(expr)) return(.prime_number(expr, assumptions))
  args <- .sargs(expr)
  for (arg in args) if (!isTRUE(.askr("integer", arg, assumptions))) return(NA)
  for (arg in args) if (.is_number(arg) && isTRUE(.fact(arg, "composite"))) return(FALSE)
  NA
})

.ask_register("prime", "Pow", function(expr, assumptions) {
  if (.is_number(expr)) return(.prime_number(expr, assumptions))
  b <- .pbase(expr); x <- .pexp(expr)
  if (isTRUE(.askr("integer", x, assumptions)) && isTRUE(.askr("integer", b, assumptions))) {
    prime_base <- .askr("prime", b, assumptions)
    if (identical(prime_base, FALSE)) return(FALSE)
    is_exp_one <- .ask_recursive(.apply_pred(.qp("eq"), list(x, S_One)), assumptions)
    if (identical(is_exp_one, FALSE)) return(FALSE)
    if (isTRUE(prime_base) && isTRUE(is_exp_one)) return(TRUE)
  }
  NA
})

.ask_register("prime", "Integer", function(expr, assumptions) isTRUE(isprime(gmp::as.bigz(expr$value))))
.ask_register("prime", c("Rational", "Infinity", "NegativeInfinity", "ImaginaryUnit"), function(expr, assumptions) FALSE)
.ask_register("prime", c("Float", "NumberSymbol"), function(expr, assumptions) .prime_number(expr, assumptions))
.ask_register("prime", "NaN", function(expr, assumptions) NA)

# ---------------------------------------------------------------------------
# composite

.ask_register("composite", "Expr", .old_fact_handler("composite"))

.ask_register("composite", "Basic", function(expr, assumptions) {
  positive <- .askr("positive", expr, assumptions)
  if (isTRUE(positive)) {
    integer <- .askr("integer", expr, assumptions)
    if (isTRUE(integer)) {
      prime <- .askr("prime", expr, assumptions)
      if (is.na(prime)) return(NA)
      is_one <- .ask_recursive(.apply_pred(.qp("eq"), list(expr, S_One)), assumptions)
      if (isTRUE(is_one)) return(FALSE)
      if (is.na(is_one)) return(NA)
      return(!prime)
    }
    return(integer)
  }
  positive
})

# ---------------------------------------------------------------------------
# even

.ask_register("even", "Expr", .old_fact_handler("even"))

.ask_register("even", "Basic", function(expr, assumptions) {
  if (.is_number(expr)) return(.even_number(expr, assumptions))
  NA
})

.ask_register("even", "Mul", function(expr, assumptions) {
  if (.is_number(expr)) return(.even_number(expr, assumptions))
  even <- FALSE; odd <- 0L; irrational <- FALSE; acc <- S_One
  args <- .sargs(expr)
  for (arg in args) {
    if (isTRUE(.askr("integer", arg, assumptions))) {
      if (isTRUE(.askr("even", arg, assumptions))) even <- TRUE
      else if (isTRUE(.askr("odd", arg, assumptions))) odd <- odd + 1L
      else if (!even && !.keq(acc, S_One)) {
        if (isTRUE(.askr("odd", .add(list(acc, arg)), assumptions))) even <- TRUE
      }
    } else if (isTRUE(.askr("irrational", arg, assumptions))) {
      if (irrational) return(NA)
      irrational <- TRUE
    } else {
      return(NA)
    }
    acc <- arg
  }
  if (irrational) return(FALSE)
  if (even) return(TRUE)
  if (odd == length(args)) return(FALSE)
  NA
})

.ask_register("even", "Add", function(expr, assumptions) {
  if (.is_number(expr)) return(.even_number(expr, assumptions))
  result <- TRUE
  for (arg in .sargs(expr)) {
    if (isTRUE(.askr("even", arg, assumptions))) next
    else if (isTRUE(.askr("odd", arg, assumptions))) result <- !result
    else return(NA)
  }
  result
})

.ask_register("even", "Pow", function(expr, assumptions) {
  if (.is_number(expr)) return(.even_number(expr, assumptions))
  b <- .pbase(expr); x <- .pexp(expr)
  if (isTRUE(.askr("integer", x, assumptions))) {
    if (isTRUE(.askr("positive", x, assumptions))) return(.askr("even", b, assumptions))
    else if (isTRUE(.ask_recursive(And(.QN("negative", x), .QA("odd", b)), assumptions))) return(FALSE)
    else if (.keq(b, S_NegativeOne)) return(FALSE)
  }
  NA
})

.ask_register("even", "Integer", function(expr, assumptions) gmp::mod.bigz(gmp::as.bigz(expr$value), 2) == 0)
.ask_register("even", c("Rational", "Infinity", "NegativeInfinity", "ImaginaryUnit"), function(expr, assumptions) FALSE)
.ask_register("even", "NumberSymbol", function(expr, assumptions) .even_number(expr, assumptions))

.ask_register("even", c("Abs", "re"), function(expr, assumptions) {
  if (isTRUE(.askr("real", expr$args[[1L]], assumptions))) return(.askr("even", expr$args[[1L]], assumptions))
  NA
})

.ask_register("even", "im", function(expr, assumptions) {
  if (isTRUE(.askr("real", expr$args[[1L]], assumptions))) return(TRUE)
  NA
})

.ask_register("even", "NaN", function(expr, assumptions) NA)

# ---------------------------------------------------------------------------
# odd

.ask_register("odd", "Expr", .old_fact_handler("odd"))

.ask_register("odd", "Basic", function(expr, assumptions) {
  integer <- .askr("integer", expr, assumptions)
  if (isTRUE(integer)) {
    even <- .askr("even", expr, assumptions)
    if (is.na(even)) return(NA)
    return(!even)
  }
  integer
})
