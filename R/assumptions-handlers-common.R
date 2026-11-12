# Handlers of sympy/assumptions/handlers/common.py and calculus.py: the
# commutative and is_true predicates, and finite, infinite,
# positive_infinite and negative_infinite.

# Arguments in SymPy's order, base and exponent of a power or exp().
.sargs <- function(e) if (e$type %in% c("Add", "Mul")) .sympy_args(e) else e$args
.pbase <- function(e) e$args[[1L]]
.pexp <- function(e) if (e$type == "Pow") e$args[[2L]] else e$args[[1L]]
.is_E <- function(e) is_sym(e) && e$type == "Constant" && identical(e$name, "E")
.is_fn <- function(e, nm) is_sym(e) && e$type == "Function" && e$name %in% nm
.keq <- function(a, b) identical(a$key, b$key)
.is_symbol_node <- function(e) is_sym(e) && e$type %in% c("Symbol", "Dummy", "Wild")

# Python's "if expr > 0:" on a relation: a TypeError when it does not
# evaluate to true or false.
.rel_truth <- function(r) {
  if (.is_true(r)) return(TRUE)
  if (.is_false(r)) return(FALSE)
  .ask_type_error(paste0("cannot determine truth value of Relational: ", .str(r)))
}

# ---------------------------------------------------------------------------
# commutative

.ask_register("commutative", "Symbol", function(expr, assumptions) {
  keys <- vapply(.ask_conjuncts(assumptions), .k, "")
  c0 <- .fact(expr, "commutative")
  if (!is.na(c0)) return(c0 && !(.QN("commutative", expr)$key %in% keys))
  if (.QA("commutative", expr)$key %in% keys) return(TRUE)
  if (.QN("commutative", expr)$key %in% keys) return(FALSE)
  TRUE
})

.ask_register("commutative", "Basic", function(expr, assumptions) {
  for (arg in .sargs(expr)) if (!isTRUE(.askr("commutative", arg, assumptions))) return(FALSE)
  TRUE
})

.ask_register("commutative", c("Number", "NaN"), function(expr, assumptions) TRUE)

# ---------------------------------------------------------------------------
# is_true

.ask_register("is_true", "bool", function(expr, assumptions) expr)
.ask_register("is_true", "BooleanTrue", function(expr, assumptions) TRUE)
.ask_register("is_true", "BooleanFalse", function(expr, assumptions) FALSE)
.ask_register("is_true", "AppliedPredicate", function(expr, assumptions) .ask_recursive(expr, assumptions))

.ask_register("is_true", "Not", function(expr, assumptions) {
  arg <- expr$args[[1L]]
  if (.is_symbol_node(arg)) return(NA)
  value <- .ask_recursive(arg, assumptions)
  if (is.na(value)) NA else !value
})

.ask_register("is_true", "Or", function(expr, assumptions) {
  result <- FALSE
  for (arg in expr$args) {
    p <- .ask_recursive(arg, assumptions)
    if (isTRUE(p)) return(TRUE)
    if (is.na(p)) result <- NA
  }
  result
})

.ask_register("is_true", "And", function(expr, assumptions) {
  result <- TRUE
  for (arg in expr$args) {
    p <- .ask_recursive(arg, assumptions)
    if (identical(p, FALSE)) return(FALSE)
    if (is.na(p)) result <- NA
  }
  result
})

.ask_register("is_true", "Implies", function(expr, assumptions) {
  .ask_recursive(Or(Not(expr$args[[1L]]), expr$args[[2L]]), assumptions)
})

.ask_register("is_true", "Equivalent", function(expr, assumptions) {
  if (length(expr$args) != 2L)
    .ask_value_error(sprintf("too many values to unpack (expected 2, got %d)", length(expr$args)))
  pt <- .ask_recursive(expr$args[[1L]], assumptions)
  if (is.na(pt)) return(NA)
  qt <- .ask_recursive(expr$args[[2L]], assumptions)
  if (is.na(qt)) return(NA)
  pt == qt
})

# ---------------------------------------------------------------------------
# finite

.ask_register("finite", "Symbol", function(expr, assumptions) {
  f <- .fact(expr, "finite")
  if (!is.na(f)) return(f)
  if (.in_conjuncts(.QA("finite", expr), assumptions)) return(TRUE)
  NA
})

# Python's comparisons between -1, True, False and None in the Add handler
.ask_register("finite", "Add", function(expr, assumptions) {
  sign <- -1
  result <- TRUE
  for (arg in .sargs(expr)) {
    bounded <- .askr("finite", arg, assumptions)
    if (isTRUE(bounded)) next
    s <- .askr("extended_positive", arg, assumptions)
    sign_m1 <- identical(sign, -1)
    s_ne_sign <- sign_m1 || !identical(s, sign)
    if ((!sign_m1 && s_ne_sign) || (is.na(s) && (is.na(bounded) || (!sign_m1 && is.na(sign))))) return(NA)
    sign <- s
    if (!identical(result, FALSE)) result <- bounded
  }
  result
})

.ask_register("finite", "Mul", function(expr, assumptions) {
  result <- TRUE
  possible_zero <- FALSE
  for (arg in .sargs(expr)) {
    bounded <- .askr("finite", arg, assumptions)
    if (isTRUE(bounded)) {
      if (!identical(.askr("zero", arg, assumptions), FALSE)) {
        if (identical(result, FALSE)) return(NA)
        possible_zero <- TRUE
      }
    } else if (is.na(bounded)) {
      if (is.na(result)) return(NA)
      if (is.na(.askr("extended_nonzero", arg, assumptions))) return(NA)
      if (!identical(result, FALSE)) result <- NA
    } else {
      if (possible_zero) return(NA)
      result <- FALSE
    }
  }
  result
})

.ask_register("finite", "Pow", function(expr, assumptions) {
  b <- .pbase(expr); x <- .pexp(expr)
  if (.is_E(b)) return(.askr("finite", x, assumptions))
  base_bounded <- .askr("finite", b, assumptions)
  exp_bounded <- .askr("finite", x, assumptions)
  if (is.na(base_bounded) && is.na(exp_bounded)) return(NA)
  if (identical(base_bounded, FALSE) && isTRUE(.askr("extended_nonzero", x, assumptions))) return(FALSE)
  if (isTRUE(base_bounded) && isTRUE(exp_bounded)) {
    is_base_zero <- .askr("zero", b, assumptions)
    is_exp_negative <- .askr("negative", x, assumptions)
    if (isTRUE(is_base_zero) && isTRUE(is_exp_negative)) return(FALSE)
    if (!identical(is_base_zero, FALSE) && !identical(is_exp_negative, FALSE)) return(NA)
    return(TRUE)
  }
  ab <- Abs(b)
  if (.is_true(Le(ab, S_One)) && isTRUE(.askr("extended_positive", x, assumptions))) return(TRUE)
  if (.is_true(Ge(ab, S_One)) && isTRUE(.askr("extended_negative", x, assumptions))) return(TRUE)
  if (.is_true(Ge(ab, S_One)) && identical(exp_bounded, FALSE)) return(FALSE)
  NA
})

.ask_register("finite", "exp", function(expr, assumptions) .askr("finite", expr$args[[1L]], assumptions))

.ask_register("finite", "log", function(expr, assumptions) {
  if (isTRUE(.askr("infinite", expr$args[[1L]], assumptions))) return(FALSE)
  .ask_recursive(.QN("zero", expr$args[[1L]]), assumptions)
})

.ask_register("finite", c("cos", "sin", "Number", "Pi", "Exp1", "GoldenRatio", "TribonacciConstant",
                          "ImaginaryUnit", "sign"), function(expr, assumptions) TRUE)
.ask_register("finite", c("ComplexInfinity", "Infinity", "NegativeInfinity"), function(expr, assumptions) FALSE)
.ask_register("finite", "NaN", function(expr, assumptions) NA)

# ---------------------------------------------------------------------------
# infinite, positive_infinite, negative_infinite

.ask_register("infinite", "Expr", function(expr, assumptions) {
  is_finite <- .askr("finite", expr, assumptions)
  if (is.na(is_finite)) NA else !is_finite
})

.ask_register("positive_infinite", "Infinity", function(expr, assumptions) TRUE)
.ask_register("positive_infinite", c("NegativeInfinity", "ComplexInfinity"), function(expr, assumptions) FALSE)
.ask_register("positive_infinite", "Expr", function(expr, assumptions) {
  if (isTRUE(.askr("complex", expr, assumptions))) return(FALSE)
  if (isTRUE(.askr("finite", expr, assumptions))) return(FALSE)
  NA
})

.ask_register("negative_infinite", "NegativeInfinity", function(expr, assumptions) TRUE)
.ask_register("negative_infinite", c("Infinity", "ComplexInfinity"), function(expr, assumptions) FALSE)
.ask_register("negative_infinite", "Expr", function(expr, assumptions) {
  if (isTRUE(.askr("complex", expr, assumptions))) return(FALSE)
  if (isTRUE(.askr("finite", expr, assumptions))) return(FALSE)
  NA
})
