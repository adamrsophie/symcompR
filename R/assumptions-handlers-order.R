# Handlers of sympy/assumptions/handlers/order.py: negative, nonnegative,
# nonzero, zero, nonpositive, positive and their extended forms.

# A literal zero, which is what Python's "not x" tests on a number.
.py_falsy <- function(e) {
  if (!is_sym(e)) return(isTRUE(e == 0))
  if (is_exact_number(e)) return(e$value == 0)
  if (e$type == "Float") return(.float_sign(e) == 0)
  FALSE
}

.float_sign <- function(e) {
  m <- .float_mp(e)$m
  if (m > 0) 1 else if (m < 0) -1 else 0
}

# The sign of x.evalf(2) for a real number x, or NA when that evaluation
# has no significant digits (SymPy's _prec == 1).
.ev2_sign <- function(r) {
  if (!is_sym(r)) return(NA)
  if (is_exact_number(r)) return(if (r$value > 0) 1 else if (r$value < 0) -1 else 0)
  if (r$type == "Float") return(.float_sign(r))
  if (r$type == "Constant" && r$name %in% c("oo", "-oo")) return(if (r$name == "oo") 1 else -1)
  s <- .mp_num_sign(r)
  if (is.na(s)) NA else s
}

.is_nan_sym <- function(e) is_sym(e) && e$type == "Constant" && identical(e$name, "nan")

# _NegativePredicate_number and _PositivePredicate_number
.sign_number <- function(expr, negative) {
  ri <- .as_real_imag(expr)
  r <- ri[[1L]]; i <- ri[[2L]]
  if (negative && (.is_nan_sym(r) || .is_nan_sym(i))) return(NA)
  test <- function(s) if (negative) s < 0 else s > 0
  if (.py_falsy(i)) {
    s <- .ev2_sign(r)
    if (!is.na(s)) return(test(s))
  } else {
    si <- .ev2_sign(i)
    if (!is.na(si)) {
      if (si != 0) return(FALSE)
      s <- .ev2_sign(r)
      if (!is.na(s)) return(test(s))
    }
  }
  NA
}
.neg_number <- function(expr, assumptions) .sign_number(expr, TRUE)
.pos_number <- function(expr, assumptions) .sign_number(expr, FALSE)

# ---------------------------------------------------------------------------
# negative

.ask_register("negative", "Basic", function(expr, assumptions) {
  if (.is_number(expr)) return(.neg_number(expr, assumptions))
  NA
})
.ask_register("negative", "Expr", .old_fact_handler("negative"))

.ask_register("negative", "Add", function(expr, assumptions) {
  if (.is_number(expr)) return(.neg_number(expr, assumptions))
  r <- .askr("real", expr, assumptions)
  if (!isTRUE(r)) return(r)
  nonpos <- 0L
  args <- .sargs(expr)
  for (arg in args) {
    if (!isTRUE(.askr("negative", arg, assumptions))) {
      if (identical(.askr("positive", arg, assumptions), FALSE)) nonpos <- nonpos + 1L
      else return(NA)
    }
  }
  if (nonpos < length(args)) return(TRUE)
  NA
})

.ask_register("negative", "Mul", function(expr, assumptions) {
  if (.is_number(expr)) return(.neg_number(expr, assumptions))
  result <- NA
  for (arg in .sargs(expr)) {
    if (is.na(result)) result <- FALSE
    if (isTRUE(.askr("negative", arg, assumptions))) result <- !result
    else if (isTRUE(.askr("positive", arg, assumptions))) NULL
    else return(NA)
  }
  result
})

.ask_register("negative", "Pow", function(expr, assumptions) {
  b <- .pbase(expr); x <- .pexp(expr)
  if (.is_E(b)) {
    if (isTRUE(.askr("real", x, assumptions))) return(FALSE)
    return(NA)
  }
  if (.is_number(expr)) return(.neg_number(expr, assumptions))
  if (isTRUE(.askr("real", b, assumptions))) {
    if (isTRUE(.askr("positive", b, assumptions)))
      if (isTRUE(.askr("real", x, assumptions))) return(FALSE)
    if (isTRUE(.askr("even", x, assumptions))) return(FALSE)
    if (isTRUE(.askr("odd", x, assumptions))) return(.askr("negative", b, assumptions))
  }
  NA
})

.ask_register("negative", c("Abs", "ImaginaryUnit"), function(expr, assumptions) FALSE)

.ask_register("negative", "exp", function(expr, assumptions) {
  if (isTRUE(.askr("real", expr$args[[1L]], assumptions))) return(FALSE)
  .MDNI
})

# ---------------------------------------------------------------------------
# nonnegative

.ask_register("nonnegative", "Basic", function(expr, assumptions) {
  if (.is_number(expr)) {
    notnegative <- .fuzzy_not(.neg_number(expr, assumptions))
    if (isTRUE(notnegative)) return(.askr("real", expr, assumptions))
    return(notnegative)
  }
  NA
})
.ask_register("nonnegative", "Expr", .old_fact_handler("nonnegative"))

# ---------------------------------------------------------------------------
# nonzero

.ask_register("nonzero", "Expr", .old_fact_handler("nonzero"))

.ask_register("nonzero", "Basic", function(expr, assumptions) {
  if (identical(.askr("real", expr, assumptions), FALSE)) return(FALSE)
  if (.is_number(expr)) {
    parts <- .as_real_imag(expr)
    nonz <- function(p) {
      if (.py_falsy(p)) return(FALSE)
      s <- .ev2_sign(p)
      if (is.na(s)) NA else s != 0
    }
    return(.fuzzy_or_lazy(2L, function(k) nonz(parts[[k]])))
  }
  NA
})

.ask_register("nonzero", "Add", function(expr, assumptions) {
  args <- .sargs(expr)
  if (.py_all(length(args), function(i) .askr("positive", args[[i]], assumptions)) ||
      .py_all(length(args), function(i) .askr("negative", args[[i]], assumptions))) return(TRUE)
  NA
})

.ask_register("nonzero", "Mul", function(expr, assumptions) {
  for (arg in .sargs(expr)) {
    result <- .askr("nonzero", arg, assumptions)
    if (isTRUE(result)) next
    return(result)
  }
  TRUE
})

.ask_register("nonzero", "Pow", function(expr, assumptions) {
  .fuzzy_and_vals(list(.fuzzy_or_vals(list(.askr("nonzero", .pbase(expr), assumptions),
                                            .askr("zero", .pexp(expr), assumptions))),
                       .askr("real", expr, assumptions)))
})

.ask_register("nonzero", "Abs", function(expr, assumptions) .askr("nonzero", expr$args[[1L]], assumptions))
.ask_register("nonzero", "NaN", function(expr, assumptions) NA)

# ---------------------------------------------------------------------------
# zero

.ask_register("zero", "Expr", .old_fact_handler("zero"))

.ask_register("zero", "Basic", function(expr, assumptions) {
  .fuzzy_and_vals(list(.fuzzy_not(.askr("nonzero", expr, assumptions)), .askr("real", expr, assumptions)))
})

.ask_register("zero", "Mul", function(expr, assumptions) {
  args <- .sargs(expr)
  .fuzzy_or_lazy(length(args), function(i) .askr("zero", args[[i]], assumptions))
})

# ---------------------------------------------------------------------------
# nonpositive

.ask_register("nonpositive", "Expr", .old_fact_handler("nonpositive"))

.ask_register("nonpositive", "Basic", function(expr, assumptions) {
  if (.is_number(expr)) {
    notpositive <- .fuzzy_not(.pos_number(expr, assumptions))
    if (isTRUE(notpositive)) return(.askr("real", expr, assumptions))
    return(notpositive)
  }
  NA
})

# ---------------------------------------------------------------------------
# positive

.ask_register("positive", "Expr", .old_fact_handler("positive"))

.ask_register("positive", "Basic", function(expr, assumptions) {
  if (.is_number(expr)) return(.pos_number(expr, assumptions))
  NA
})

.ask_register("positive", "Mul", function(expr, assumptions) {
  if (.is_number(expr)) return(.pos_number(expr, assumptions))
  result <- TRUE
  for (arg in .sargs(expr)) {
    if (isTRUE(.askr("positive", arg, assumptions))) next
    else if (isTRUE(.askr("negative", arg, assumptions))) result <- !result
    else return(NA)
  }
  result
})

.ask_register("positive", "Add", function(expr, assumptions) {
  if (.is_number(expr)) return(.pos_number(expr, assumptions))
  r <- .askr("real", expr, assumptions)
  if (!isTRUE(r)) return(r)
  nonneg <- 0L
  args <- .sargs(expr)
  for (arg in args) {
    if (!isTRUE(.askr("positive", arg, assumptions))) {
      if (identical(.askr("negative", arg, assumptions), FALSE)) nonneg <- nonneg + 1L
      else return(NA)
    }
  }
  if (nonneg < length(args)) return(TRUE)
  NA
})

.ask_register("positive", "Pow", function(expr, assumptions) {
  b <- .pbase(expr); x <- .pexp(expr)
  if (.is_E(b)) {
    if (isTRUE(.askr("real", x, assumptions))) return(TRUE)
    if (isTRUE(.askr("imaginary", x, assumptions)))
      return(.askr("even", .div(x, .mul(list(S_I, S_Pi))), assumptions))
    return(NA)
  }
  if (.is_number(expr)) return(.pos_number(expr, assumptions))
  if (isTRUE(.askr("even", x, assumptions)) && isTRUE(.askr("real", b, assumptions))) {
    zero_base <- .askr("zero", b, assumptions)
    if (isTRUE(zero_base) && isTRUE(.askr("positive", x, assumptions))) return(FALSE)
    else if (identical(zero_base, FALSE)) return(TRUE)
  }
  if (isTRUE(.askr("positive", b, assumptions)))
    if (isTRUE(.askr("real", x, assumptions))) return(TRUE)
  if (isTRUE(.askr("negative", b, assumptions)))
    if (isTRUE(.askr("odd", x, assumptions))) return(FALSE)
  NA
})

.ask_register("positive", "exp", function(expr, assumptions) {
  x <- expr$args[[1L]]
  if (isTRUE(.askr("real", x, assumptions))) return(TRUE)
  if (isTRUE(.askr("imaginary", x, assumptions)))
    return(.askr("even", .div(x, .mul(list(S_I, S_Pi))), assumptions))
  NA
})

.ask_register("positive", "log", function(expr, assumptions) {
  a <- expr$args[[1L]]
  r <- .askr("real", a, assumptions)
  if (!isTRUE(r)) return(r)
  if (isTRUE(.askr("positive", .sub(a, S_One), assumptions))) return(TRUE)
  if (isTRUE(.askr("negative", .sub(a, S_One), assumptions))) return(FALSE)
  NA
})

.ask_register("positive", "factorial", function(expr, assumptions) {
  x <- expr$args[[1L]]
  if (isTRUE(.ask_recursive(And(.QA("integer", x), .QA("positive", x)), assumptions))) return(TRUE)
  NA
})

.ask_register("positive", "ImaginaryUnit", function(expr, assumptions) FALSE)
.ask_register("positive", "Abs", function(expr, assumptions) .askr("nonzero", expr, assumptions))

.ask_register("positive", c("Trace", "Determinant"), function(expr, assumptions) {
  if (isTRUE(.askr("positive_definite", expr$args[[1L]], assumptions))) return(TRUE)
  NA
})

.ask_register("positive", "MatrixElement", function(expr, assumptions) {
  if (.keq(expr$args[[2L]], expr$args[[3L]]) && isTRUE(.askr("positive_definite", expr$args[[1L]], assumptions)))
    return(TRUE)
  NA
})

.ask_register("positive", "atan", function(expr, assumptions) .askr("positive", expr$args[[1L]], assumptions))

.ask_register("positive", "asin", function(expr, assumptions) {
  x <- expr$args[[1L]]
  if (isTRUE(.ask_recursive(And(.QA("positive", x), .QA("nonpositive", .sub(x, S_One))), assumptions))) return(TRUE)
  if (isTRUE(.ask_recursive(And(.QA("negative", x), .QA("nonnegative", .add(list(x, S_One)))), assumptions))) return(FALSE)
  NA
})

.ask_register("positive", "acos", function(expr, assumptions) {
  x <- expr$args[[1L]]
  if (isTRUE(.ask_recursive(And(.QA("nonpositive", .sub(x, S_One)), .QA("nonnegative", .add(list(x, S_One)))), assumptions)))
    return(TRUE)
  NA
})

.ask_register("positive", "acot", function(expr, assumptions) .askr("real", expr$args[[1L]], assumptions))
.ask_register("positive", "NaN", function(expr, assumptions) NA)

# ---------------------------------------------------------------------------
# extended predicates

.ask_register("extended_negative", "object", function(expr, assumptions)
  .ask_recursive(Or(.QA("negative", expr), .QA("negative_infinite", expr)), assumptions))
.ask_register("extended_positive", "object", function(expr, assumptions)
  .ask_recursive(Or(.QA("positive", expr), .QA("positive_infinite", expr)), assumptions))
.ask_register("extended_nonzero", "object", function(expr, assumptions)
  .ask_recursive(Or(.QA("negative_infinite", expr), .QA("negative", expr), .QA("positive", expr),
                    .QA("positive_infinite", expr)), assumptions))
.ask_register("extended_nonpositive", "object", function(expr, assumptions)
  .ask_recursive(Or(.QA("negative_infinite", expr), .QA("negative", expr), .QA("zero", expr)), assumptions))
.ask_register("extended_nonnegative", "object", function(expr, assumptions)
  .ask_recursive(Or(.QA("zero", expr), .QA("positive", expr), .QA("positive_infinite", expr)), assumptions))
