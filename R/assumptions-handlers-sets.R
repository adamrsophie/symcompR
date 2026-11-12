# Handlers of sympy/assumptions/handlers/sets.py: integer, rational,
# irrational, real, extended_real, hermitian, complex, imaginary,
# antihermitian, algebraic and transcendental.

# int(expr.round()) for a number: a bigz, or NULL where Python raises a
# TypeError (a complex value, an infinity, nan, or no number at all).
.py_round_int <- function(expr) {
  if (!is_sym(expr)) return(NULL)
  rq <- function(q) {
    f <- gmp::as.bigz(floor(q))
    d <- q - f
    if (d > gmp::as.bigq(1, 2) || (d == gmp::as.bigq(1, 2) && gmp::mod.bigz(f, 2) == 1)) f + 1 else f
  }
  if (expr$type == "Integer") return(gmp::as.bigz(expr$value))
  if (expr$type == "Rational") return(rq(expr$value))
  if (expr$type == "Float") return(rq(.mp_to_q(.float_mp(expr))))
  if (!.is_number(expr) || (expr$type == "Constant" && expr$name %in% c("oo", "-oo", "zoo", "nan"))) return(NULL)
  z <- tryCatch(.mpeval_adaptive(expr, 64L, list(), maxprec = 2000), error = function(e) NULL)
  if (is.null(z) || z$im$m != 0) return(NULL)
  rq(.mp_to_q(z$re))
}

# (d).equals(0) for a number d: TRUE, FALSE or NA.
.num_equals0 <- function(d) {
  if (!is_sym(d)) return(isTRUE(d == 0))
  if (is_exact_number(d)) return(d$value == 0)
  if (d$type == "Float") return(.float_sign(d) == 0)
  z <- tryCatch(.mpeval_adaptive(d, 30L, list(), maxprec = 2000), error = function(e) NULL)
  if (!is.null(z) && isTRUE(attr(z, "converged")) && (z$re$m != 0 || z$im$m != 0)) return(FALSE)
  s <- tryCatch(simplify(d), error = function(e) d)
  if (is_sym(s) && is_zero_num(s)) return(TRUE)
  NA
}

.int_number <- function(expr, assumptions) {
  i <- .py_round_int(expr)
  if (is.null(i)) return(FALSE)
  isTRUE(.num_equals0(.sub(expr, make_rational(gmp::as.bigq(i)))))
}

# ---------------------------------------------------------------------------
# integer

.ask_register("integer", c("int", "Integer"), function(expr, assumptions) TRUE)
.ask_register("integer", c("Exp1", "GoldenRatio", "ImaginaryUnit", "Infinity", "NegativeInfinity", "Pi",
                           "Rational", "TribonacciConstant"), function(expr, assumptions) FALSE)
.ask_register("integer", "Expr", .old_fact_handler("integer"))

.ask_register("integer", "Add", function(expr, assumptions) {
  if (.is_number(expr)) return(.int_number(expr, assumptions))
  .test_closed_group(expr, assumptions, "integer")
})

.ask_register("integer", "Pow", function(expr, assumptions) {
  b <- .pbase(expr); x <- .pexp(expr)
  if (.is_number(expr)) return(.int_number(expr, assumptions))
  if (isTRUE(.ask_all(list(.QN("zero", b), .QA("finite", b), .QA("zero", x)), assumptions))) return(TRUE)
  if (isTRUE(.ask_all(list(.QA("integer", b), .QA("integer", x)), assumptions))) {
    neg_exp <- .askr("negative", x, assumptions)
    if (identical(neg_exp, FALSE)) return(TRUE)
    base_pm_1 <- .ask_any(list(.QA("zero", .sub(b, S_One)), .QA("zero", .add(list(b, S_One)))), assumptions)
    if (isTRUE(base_pm_1)) return(TRUE)
    if (isTRUE(neg_exp) && identical(base_pm_1, FALSE)) return(FALSE)
  }
  NA
})

.ask_register("integer", "Mul", function(expr, assumptions) {
  if (.is_number(expr)) return(.int_number(expr, assumptions))
  output <- TRUE
  for (arg in .sargs(expr)) {
    if (!isTRUE(.askr("integer", arg, assumptions))) {
      if (is_exact_number(arg)) {
        if (gmp::denominator(arg$value) == 2) return(.askr("even", .mul(list(.int(2L), expr)), assumptions))
        return(NA)
      } else if (isTRUE(.askr("irrational", arg, assumptions))) {
        if (output) output <- FALSE else return(NA)
      } else {
        return(NA)
      }
    }
  }
  output
})

.ask_register("integer", "Abs", function(expr, assumptions) {
  if (isTRUE(.askr("integer", expr$args[[1L]], assumptions))) return(TRUE)
  NA
})

.ask_register("integer", c("Determinant", "MatrixElement", "Trace"), function(expr, assumptions)
  .askr("integer_elements", expr$args[[1L]], assumptions))

# ---------------------------------------------------------------------------
# rational

.ask_register("rational", "Rational", function(expr, assumptions) TRUE)
.ask_register("rational", "Float", function(expr, assumptions) NA)
.ask_register("rational", c("Exp1", "GoldenRatio", "ImaginaryUnit", "Infinity", "NegativeInfinity", "Pi",
                            "TribonacciConstant"), function(expr, assumptions) FALSE)
.ask_register("rational", "Expr", .old_fact_handler("rational"))

.ask_register("rational", c("Add", "Mul"), function(expr, assumptions) {
  if (.is_number(expr) && !.py_falsy(.as_real_imag(expr)[[2L]])) return(FALSE)
  .test_closed_group(expr, assumptions, "rational")
})

.ask_register("rational", "Pow", function(expr, assumptions) {
  b <- .pbase(expr); x <- .pexp(expr)
  if (.is_E(b)) {
    if (isTRUE(.askr("rational", x, assumptions))) return(.askr("zero", x, assumptions))
    return(NA)
  }
  is_exp_integer <- .askr("integer", x, assumptions)
  if (isTRUE(is_exp_integer)) {
    is_base_rational <- .askr("rational", b, assumptions)
    if (isTRUE(is_base_rational)) {
      is_base_zero <- .askr("zero", b, assumptions)
      if (identical(is_base_zero, FALSE)) return(TRUE)
      if (isTRUE(is_base_zero) && isTRUE(.askr("positive", x, assumptions))) return(TRUE)
    }
    if (identical(.askr("algebraic", b, assumptions), FALSE)) return(.askr("zero", x, assumptions))
    if (isTRUE(.askr("irrational", b, assumptions)) &&
        isTRUE(.ask_recursive(.apply_pred(.qp("eq"), list(x, S_NegativeOne)), assumptions))) return(FALSE)
    return(NA)
  } else if (isTRUE(.askr("rational", x, assumptions))) {
    if (isTRUE(.askr("prime", b, assumptions)) && identical(is_exp_integer, FALSE)) return(FALSE)
    if (isTRUE(.askr("zero", b, assumptions)) && isTRUE(.askr("positive", x, assumptions))) return(TRUE)
    if (isTRUE(.ask_recursive(.apply_pred(.qp("eq"), list(b, S_One)), assumptions))) return(TRUE)
  }
  NA
})

.ask_register("rational", c("asin", "atan", "cos", "sin", "tan"), function(expr, assumptions) {
  x <- expr$args[[1L]]
  if (isTRUE(.askr("rational", x, assumptions))) return(.ask_recursive(.QN("nonzero", x), assumptions))
  NA
})

.ask_register("rational", "exp", function(expr, assumptions) {
  x <- expr$args[[1L]]
  if (isTRUE(.askr("rational", x, assumptions))) return(.ask_recursive(.QN("nonzero", x), assumptions))
  NA
})

.ask_register("rational", c("acot", "cot"), function(expr, assumptions) {
  if (isTRUE(.askr("rational", expr$args[[1L]], assumptions))) return(FALSE)
  NA
})

.ask_register("rational", c("acos", "log"), function(expr, assumptions) {
  x <- expr$args[[1L]]
  if (isTRUE(.askr("rational", x, assumptions))) return(.ask_recursive(.QN("nonzero", .sub(x, S_One)), assumptions))
  NA
})

# ---------------------------------------------------------------------------
# irrational

.ask_register("irrational", "Expr", .old_fact_handler("irrational"))

.ask_register("irrational", "Basic", function(expr, assumptions) {
  real <- .askr("real", expr, assumptions)
  if (isTRUE(real)) {
    rational <- .askr("rational", expr, assumptions)
    if (is.na(rational)) return(NA)
    return(!rational)
  }
  real
})

# ---------------------------------------------------------------------------
# real

.real_number <- function(expr, assumptions) {
  s <- .ev2_sign(.as_real_imag(expr)[[2L]])
  if (is.na(s)) NA else s == 0
}

.ask_register("real", c("Abs", "Exp1", "Float", "GoldenRatio", "im", "Pi", "Rational", "re",
                        "TribonacciConstant"), function(expr, assumptions) TRUE)
.ask_register("real", c("ImaginaryUnit", "Infinity", "NegativeInfinity"), function(expr, assumptions) FALSE)
.ask_register("real", "Expr", .old_fact_handler("real"))

.ask_register("real", "Add", function(expr, assumptions) {
  if (.is_number(expr)) return(.real_number(expr, assumptions))
  .test_closed_group(expr, assumptions, "real")
})

.ask_register("real", "Mul", function(expr, assumptions) {
  if (.is_number(expr)) return(.real_number(expr, assumptions))
  result <- TRUE
  for (arg in .sargs(expr)) {
    if (isTRUE(.askr("real", arg, assumptions))) next
    else if (isTRUE(.askr("imaginary", arg, assumptions))) result <- !result
    else return(NA)
  }
  if (result) return(TRUE)
  .askr("zero", expr, assumptions)
})

.ask_register("real", "Pow", function(expr, assumptions) {
  b <- .pbase(expr); x <- .pexp(expr)
  if (.is_number(expr)) return(.real_number(expr, assumptions))
  if (.is_E(b))
    return(.ask_recursive(Or(.QA("integer", .div(.div(x, S_I), S_Pi)), .QA("real", x)), assumptions))
  if (.is_fn(b, "exp") || (b$type == "Pow" && .is_E(.pbase(b)))) {
    bx <- .pexp(b)
    if (isTRUE(.askr("imaginary", bx, assumptions)))
      if (isTRUE(.askr("imaginary", x, assumptions))) return(TRUE)
    i <- .div(.div(bx, S_I), S_Pi)
    if (isTRUE(.askr("integer", .mul(list(.int(2L), i)), assumptions)))
      return(.askr("real", .pow(.pow(S_NegativeOne, i), x), assumptions))
    return(NA)
  }
  if (isTRUE(.askr("imaginary", b, assumptions))) {
    if (isTRUE(.askr("integer", x, assumptions))) {
      odd <- .askr("odd", x, assumptions)
      if (!is.na(odd)) return(!odd)
      return(NA)
    }
  }
  if (isTRUE(.askr("imaginary", x, assumptions))) {
    imlog <- .askr("imaginary", fcall("log", list(b)), assumptions)
    if (!is.na(imlog)) return(imlog)
  }
  if (isTRUE(.askr("real", b, assumptions))) {
    if (isTRUE(.askr("real", x, assumptions))) {
      if (is_exact_number(x) && isTRUE(.askr("even", .int(gmp::denominator(x$value)), assumptions)))
        return(.askr("nonnegative", b, assumptions))
      base_is_zero <- .askr("zero", b, assumptions)
      if (isTRUE(base_is_zero) || isTRUE(.askr("integer", x, assumptions))) {
        cinf <- .fuzzy_and_vals(list(base_is_zero, .askr("negative", x, assumptions)))
        return(.fuzzy_not(cinf))
      } else if (isTRUE(.askr("positive", b, assumptions))) {
        return(TRUE)
      }
    }
  }
  NA
})

.ask_register("real", c("cos", "sin"), function(expr, assumptions) {
  if (isTRUE(.askr("real", expr$args[[1L]], assumptions))) return(TRUE)
  NA
})

.ask_register("real", "exp", function(expr, assumptions) {
  x <- expr$args[[1L]]
  .ask_recursive(Or(.QA("integer", .div(.div(x, S_I), S_Pi)), .QA("real", x)), assumptions)
})

.ask_register("real", "log", function(expr, assumptions) .askr("positive", expr$args[[1L]], assumptions))

.ask_register("real", c("Determinant", "MatrixElement", "Trace"), function(expr, assumptions)
  .askr("real_elements", expr$args[[1L]], assumptions))

# ---------------------------------------------------------------------------
# extended_real

.ask_register("extended_real", "object", function(expr, assumptions)
  .ask_recursive(Or(.QA("negative_infinite", expr), .QA("negative", expr), .QA("zero", expr),
                    .QA("positive", expr), .QA("positive_infinite", expr)), assumptions))
.ask_register("extended_real", c("Infinity", "NegativeInfinity"), function(expr, assumptions) TRUE)
.ask_register("extended_real", c("Add", "Mul", "Pow"), function(expr, assumptions)
  .test_closed_group(expr, assumptions, "extended_real"))

# ---------------------------------------------------------------------------
# hermitian

.is_matrixbase <- function(e) is_sym(e) && e$type == "ImmutableMatrix"

# The entry (i, j), one based, of an ImmutableMatrix node.
.imm_get <- function(e, i, j) e$args[[3L]]$args[[(i - 1L) * .int_val(e$args[[2L]]) + j]]
.imm_dims <- function(e) c(.int_val(e$args[[1L]]), .int_val(e$args[[2L]]))

.ask_register("hermitian", "object", function(expr, assumptions) {
  if (.is_matrixbase(expr)) return(NA)
  .askr("real", expr, assumptions)
})

.ask_register("hermitian", "Add", function(expr, assumptions) {
  if (.is_number(expr)) return(.MDNI)
  .test_closed_group(expr, assumptions, "hermitian")
})

# shared by the hermitian and antihermitian Mul handlers
.herm_mul <- function(expr, assumptions, result) {
  nccount <- 0L
  for (arg in .sargs(expr)) {
    if (isTRUE(.askr("antihermitian", arg, assumptions))) result <- !result
    else if (!isTRUE(.askr("hermitian", arg, assumptions))) return(NA)
    if (isTRUE(.ask_recursive(.QN("commutative", arg), assumptions))) {
      nccount <- nccount + 1L
      if (nccount > 1L) return(NA)
    }
  }
  result
}

.ask_register("hermitian", "Mul", function(expr, assumptions) {
  if (.is_number(expr)) return(.MDNI)
  .herm_mul(expr, assumptions, TRUE)
})

.ask_register("hermitian", "Pow", function(expr, assumptions) {
  if (.is_number(expr)) return(.MDNI)
  b <- .pbase(expr); x <- .pexp(expr)
  if (.is_E(b)) {
    if (isTRUE(.askr("hermitian", x, assumptions))) return(TRUE)
    return(.MDNI)
  }
  if (isTRUE(.askr("hermitian", b, assumptions)))
    if (isTRUE(.askr("integer", x, assumptions))) return(TRUE)
  .MDNI
})

.ask_register("hermitian", c("cos", "sin", "exp"), function(expr, assumptions) {
  if (isTRUE(.askr("hermitian", expr$args[[1L]], assumptions))) return(TRUE)
  .MDNI
})

.herm_matrix <- function(mat, neg) {
  d <- .imm_dims(mat)
  ret_val <- TRUE
  for (i in seq_len(d[1L])) for (j in seq.int(i, length.out = max(0L, d[2L] - i + 1L))) {
    cj <- fcall("conjugate", list(.imm_get(mat, j, i)))
    cond <- Eq(.imm_get(mat, i, j), if (neg) .neg(cj) else cj)
    if (.is_false(cond)) return(FALSE)
    if (!.is_true(cond)) ret_val <- NA
  }
  if (is.na(ret_val)) return(.MDNI)
  ret_val
}

.ask_register("hermitian", "MatrixBase", function(mat, assumptions) .herm_matrix(mat, FALSE))

# ---------------------------------------------------------------------------
# complex

.ask_register("complex", c("Abs", "cos", "exp", "im", "ImaginaryUnit", "log", "Number", "NumberSymbol",
                           "re", "sin"), function(expr, assumptions) TRUE)
.ask_register("complex", c("Infinity", "NegativeInfinity"), function(expr, assumptions) FALSE)
.ask_register("complex", "Expr", .old_fact_handler("complex"))
.ask_register("complex", c("Add", "Mul"), function(expr, assumptions) .test_closed_group(expr, assumptions, "complex"))
.ask_register("complex", "Pow", function(expr, assumptions) {
  if (.is_E(.pbase(expr))) return(TRUE)
  .test_closed_group(expr, assumptions, "complex")
})
.ask_register("complex", c("Determinant", "MatrixElement", "Trace"), function(expr, assumptions)
  .askr("complex_elements", expr$args[[1L]], assumptions))
.ask_register("complex", "NaN", function(expr, assumptions) NA)

# ---------------------------------------------------------------------------
# imaginary

.imaginary_number <- function(expr, assumptions) {
  s <- .ev2_sign(.as_real_imag(expr)[[1L]])
  if (is.na(s)) NA else s == 0
}

.ask_register("imaginary", "ImaginaryUnit", function(expr, assumptions) TRUE)
.ask_register("imaginary", "Expr", .old_fact_handler("imaginary"))

.ask_register("imaginary", "Add", function(expr, assumptions) {
  if (.is_number(expr)) return(.imaginary_number(expr, assumptions))
  reals <- 0L
  args <- .sargs(expr)
  for (arg in args) {
    if (isTRUE(.askr("imaginary", arg, assumptions))) next
    else if (isTRUE(.askr("real", arg, assumptions))) reals <- reals + 1L
    else return(NA)
  }
  if (reals == 0L) return(TRUE)
  if (reals %in% c(1L, length(args))) return(FALSE)
  NA
})

.ask_register("imaginary", "Mul", function(expr, assumptions) {
  if (.is_number(expr)) return(.imaginary_number(expr, assumptions))
  result <- FALSE
  args <- .sargs(expr)
  for (arg in args) {
    if (isTRUE(.askr("imaginary", arg, assumptions))) result <- !result
    else if (!isTRUE(.askr("real", arg, assumptions))) return(NA)
  }
  if (length(args) == 0L) return(FALSE)
  result
})

.ask_register("imaginary", "Pow", function(expr, assumptions) {
  b <- .pbase(expr); x <- .pexp(expr)
  if (.is_number(expr)) return(.imaginary_number(expr, assumptions))
  if (.is_E(b)) {
    a <- .div(.div(x, S_I), S_Pi)
    return(.ask_recursive(And(.QA("integer", .mul(list(.int(2L), a))), .QN("integer", a)), assumptions))
  }
  if (.is_fn(b, "exp") || (b$type == "Pow" && .is_E(.pbase(b)))) {
    bx <- .pexp(b)
    if (isTRUE(.askr("imaginary", bx, assumptions))) {
      if (isTRUE(.askr("imaginary", x, assumptions))) return(FALSE)
      i <- .div(.div(bx, S_I), S_Pi)
      if (isTRUE(.askr("integer", .mul(list(.int(2L), i)), assumptions)))
        return(.askr("imaginary", .pow(.pow(S_NegativeOne, i), x), assumptions))
    }
  }
  if (isTRUE(.askr("imaginary", b, assumptions))) {
    if (isTRUE(.askr("integer", x, assumptions))) {
      odd <- .askr("odd", x, assumptions)
      if (!is.na(odd)) return(odd)
      return(NA)
    }
  }
  if (isTRUE(.askr("imaginary", x, assumptions))) {
    imlog <- .askr("imaginary", fcall("log", list(b)), assumptions)
    if (!is.na(imlog)) return(FALSE)
  }
  if (isTRUE(.ask_recursive(And(.QA("real", b), .QA("real", x)), assumptions))) {
    if (isTRUE(.askr("positive", b, assumptions))) return(FALSE)
    rat <- .askr("rational", x, assumptions)
    if (!isTRUE(rat)) return(rat)
    if (isTRUE(.askr("integer", x, assumptions))) return(FALSE)
    half <- .askr("integer", .mul(list(.int(2L), x)), assumptions)
    if (isTRUE(half)) return(.askr("negative", b, assumptions))
    return(half)
  }
  NA
})

.ask_register("imaginary", "log", function(expr, assumptions) {
  a <- expr$args[[1L]]
  if (isTRUE(.askr("real", a, assumptions))) {
    if (isTRUE(.askr("positive", a, assumptions))) return(FALSE)
    return(NA)
  }
  if (.is_fn(a, "exp") || (a$type == "Pow" && .is_E(.pbase(a)))) {
    if (.pexp(a)$key %in% c(S_I$key, .neg(S_I)$key)) return(TRUE)
  }
  im <- .askr("imaginary", a, assumptions)
  if (identical(im, FALSE)) return(FALSE)
  NA
})

.ask_register("imaginary", "exp", function(expr, assumptions) {
  a <- .div(.div(expr$args[[1L]], S_I), S_Pi)
  .ask_recursive(And(.QA("integer", .mul(list(.int(2L), a))), .QN("integer", a)), assumptions)
})

.ask_register("imaginary", c("Number", "NumberSymbol"), function(expr, assumptions)
  !.py_falsy(.as_real_imag(expr)[[2L]]))

.ask_register("imaginary", "NaN", function(expr, assumptions) NA)

# ---------------------------------------------------------------------------
# antihermitian

.ask_register("antihermitian", "object", function(expr, assumptions) {
  if (.is_matrixbase(expr)) return(NA)
  if (isTRUE(.askr("zero", expr, assumptions))) return(TRUE)
  .askr("imaginary", expr, assumptions)
})

.ask_register("antihermitian", "Add", function(expr, assumptions) {
  if (.is_number(expr)) return(.MDNI)
  .test_closed_group(expr, assumptions, "antihermitian")
})

.ask_register("antihermitian", "Mul", function(expr, assumptions) {
  if (.is_number(expr)) return(.MDNI)
  .herm_mul(expr, assumptions, FALSE)
})

.ask_register("antihermitian", "Pow", function(expr, assumptions) {
  if (.is_number(expr)) return(.MDNI)
  b <- .pbase(expr); x <- .pexp(expr)
  if (isTRUE(.askr("hermitian", b, assumptions))) {
    if (isTRUE(.askr("integer", x, assumptions))) return(FALSE)
  } else if (isTRUE(.askr("antihermitian", b, assumptions))) {
    if (isTRUE(.askr("even", x, assumptions))) return(FALSE)
    else if (isTRUE(.askr("odd", x, assumptions))) return(TRUE)
  }
  .MDNI
})

.ask_register("antihermitian", "MatrixBase", function(mat, assumptions) .herm_matrix(mat, TRUE))

# ---------------------------------------------------------------------------
# algebraic

.ask_register("algebraic", c("AlgebraicNumber", "GoldenRatio", "ImaginaryUnit", "TribonacciConstant"),
              function(expr, assumptions) TRUE)
.ask_register("algebraic", "Float", function(expr, assumptions) NA)
.ask_register("algebraic", c("ComplexInfinity", "Exp1", "Infinity", "NegativeInfinity", "Pi"),
              function(expr, assumptions) FALSE)
.ask_register("algebraic", c("Add", "Mul"), function(expr, assumptions) .test_closed_group(expr, assumptions, "algebraic"))

.ask_register("algebraic", "Pow", function(expr, assumptions) {
  b <- .pbase(expr); x <- .pexp(expr)
  if (.is_E(b)) {
    if (isTRUE(.askr("algebraic", x, assumptions))) return(.ask_recursive(.QN("nonzero", x), assumptions))
    return(NA)
  }
  if (.keq(b, S_Pi)) {
    if (isTRUE(.askr("integer", x, assumptions)) && isTRUE(.askr("positive", x, assumptions))) return(FALSE)
    return(NA)
  }
  exp_rational <- .askr("rational", x, assumptions)
  base_algebraic <- .askr("algebraic", b, assumptions)
  exp_algebraic <- .askr("algebraic", x, assumptions)
  if (isTRUE(base_algebraic) && isTRUE(exp_algebraic)) {
    if (isTRUE(exp_rational)) return(TRUE)
    if (isTRUE(.ask_recursive(And(.apply_pred(.qp("ne"), list(b, S_Zero)), .apply_pred(.qp("ne"), list(b, S_One))),
                              assumptions)) && identical(exp_rational, FALSE)) return(FALSE)
  }
  exp_integer <- .askr("integer", x, assumptions)
  if (identical(base_algebraic, FALSE) && isTRUE(exp_integer)) {
    if (.rel_truth(Gt(x, S_Zero))) return(FALSE)
  }
  NA
})

.ask_register("algebraic", "Rational", function(expr, assumptions) TRUE)

.ask_register("algebraic", c("asin", "atan", "cos", "sin", "tan", "exp"), function(expr, assumptions) {
  x <- expr$args[[1L]]
  if (isTRUE(.askr("algebraic", x, assumptions))) return(.ask_recursive(.QN("nonzero", x), assumptions))
  NA
})

.ask_register("algebraic", c("acot", "cot"), function(expr, assumptions) {
  if (isTRUE(.askr("algebraic", expr$args[[1L]], assumptions))) return(FALSE)
  NA
})

.ask_register("algebraic", c("acos", "log"), function(expr, assumptions) {
  x <- expr$args[[1L]]
  if (isTRUE(.askr("algebraic", x, assumptions))) return(.ask_recursive(.QN("nonzero", .sub(x, S_One)), assumptions))
  NA
})

# ---------------------------------------------------------------------------
# transcendental

.ask_register("transcendental", "Expr", function(expr, assumptions) {
  ret <- .fact(expr, "transcendental")
  if (!is.na(ret)) return(ret)
  is_complex <- .askr("complex", expr, assumptions)
  if (isTRUE(is_complex)) return(.fuzzy_not(.askr("algebraic", expr, assumptions)))
  is_complex
})
