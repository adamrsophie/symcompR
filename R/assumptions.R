# Three valued assumptions (TRUE, FALSE, NA for unknown), after
# sympy/core/assumptions.py and sympy/core/facts.py.

.FACTS <- .ASSUMPTION_NAMES

# The deduction of sympy/core/facts.py (FactKB.deduce_all_facts) with the
# tables SymPy compiles from the rules of sympy/core/assumptions.py:
# alpha implications closed under transitivity, and beta rules (a
# conjunction implying a fact) with the facts that may trigger them.
.FR_FACTS <- c("algebraic", "commutative", "complex", "composite", "even", "extended_negative", "extended_nonnegative", "extended_nonpositive", "extended_nonzero", "extended_positive", "extended_real", "finite", "imaginary", "infinite", "integer", "irrational", "negative", "noninteger", "nonnegative", "nonpositive", "nonzero", "odd", "positive", "prime", "rational", "real", "transcendental", "zero")
.FR_IMPL <- list(
  "algebraic=F" = c("composite=F", "even=F", "integer=F", "odd=F", "prime=F", "rational=F", "zero=F"),
  "algebraic=T" = c("commutative=T", "complex=T", "finite=T", "infinite=F", "transcendental=F"),
  "commutative=F" = c("algebraic=F", "complex=F", "composite=F", "even=F", "extended_negative=F", "extended_nonnegative=F", "extended_nonpositive=F", "extended_nonzero=F", "extended_positive=F", "extended_real=F", "imaginary=F", "integer=F", "irrational=F", "negative=F", "noninteger=F", "nonnegative=F", "nonpositive=F", "nonzero=F", "odd=F", "positive=F", "prime=F", "rational=F", "real=F", "transcendental=F", "zero=F"),
  "commutative=T" = c(),
  "complex=F" = c("algebraic=F", "composite=F", "even=F", "imaginary=F", "integer=F", "irrational=F", "negative=F", "nonnegative=F", "nonpositive=F", "nonzero=F", "odd=F", "positive=F", "prime=F", "rational=F", "real=F", "transcendental=F", "zero=F"),
  "complex=T" = c("commutative=T", "finite=T", "infinite=F"),
  "composite=F" = c(),
  "composite=T" = c("algebraic=T", "commutative=T", "complex=T", "extended_negative=F", "extended_nonnegative=T", "extended_nonpositive=F", "extended_nonzero=T", "extended_positive=T", "extended_real=T", "finite=T", "imaginary=F", "infinite=F", "integer=T", "irrational=F", "negative=F", "noninteger=F", "nonnegative=T", "nonpositive=F", "nonzero=T", "positive=T", "prime=F", "rational=T", "real=T", "transcendental=F", "zero=F"),
  "even=F" = c("zero=F"),
  "even=T" = c("algebraic=T", "commutative=T", "complex=T", "extended_real=T", "finite=T", "imaginary=F", "infinite=F", "integer=T", "irrational=F", "noninteger=F", "odd=F", "rational=T", "real=T", "transcendental=F"),
  "extended_negative=F" = c("negative=F"),
  "extended_negative=T" = c("commutative=T", "composite=F", "extended_nonnegative=F", "extended_nonpositive=T", "extended_nonzero=T", "extended_positive=F", "extended_real=T", "imaginary=F", "nonnegative=F", "positive=F", "prime=F", "zero=F"),
  "extended_nonnegative=F" = c("composite=F", "extended_positive=F", "nonnegative=F", "positive=F", "prime=F", "zero=F"),
  "extended_nonnegative=T" = c("commutative=T", "extended_negative=F", "extended_real=T", "imaginary=F", "negative=F"),
  "extended_nonpositive=F" = c("extended_negative=F", "negative=F", "nonpositive=F", "zero=F"),
  "extended_nonpositive=T" = c("commutative=T", "composite=F", "extended_positive=F", "extended_real=T", "imaginary=F", "positive=F", "prime=F"),
  "extended_nonzero=F" = c("composite=F", "extended_negative=F", "extended_positive=F", "negative=F", "nonzero=F", "positive=F", "prime=F"),
  "extended_nonzero=T" = c("commutative=T", "extended_real=T", "imaginary=F", "zero=F"),
  "extended_positive=F" = c("composite=F", "positive=F", "prime=F"),
  "extended_positive=T" = c("commutative=T", "extended_negative=F", "extended_nonnegative=T", "extended_nonpositive=F", "extended_nonzero=T", "extended_real=T", "imaginary=F", "negative=F", "nonpositive=F", "zero=F"),
  "extended_real=F" = c("composite=F", "even=F", "extended_negative=F", "extended_nonnegative=F", "extended_nonpositive=F", "extended_nonzero=F", "extended_positive=F", "integer=F", "irrational=F", "negative=F", "noninteger=F", "nonnegative=F", "nonpositive=F", "nonzero=F", "odd=F", "positive=F", "prime=F", "rational=F", "real=F", "zero=F"),
  "extended_real=T" = c("commutative=T", "imaginary=F"),
  "finite=F" = c("algebraic=F", "complex=F", "composite=F", "even=F", "imaginary=F", "infinite=T", "integer=F", "irrational=F", "negative=F", "nonnegative=F", "nonpositive=F", "nonzero=F", "odd=F", "positive=F", "prime=F", "rational=F", "real=F", "transcendental=F", "zero=F"),
  "finite=T" = c("infinite=F"),
  "imaginary=F" = c(),
  "imaginary=T" = c("commutative=T", "complex=T", "composite=F", "even=F", "extended_negative=F", "extended_nonnegative=F", "extended_nonpositive=F", "extended_nonzero=F", "extended_positive=F", "extended_real=F", "finite=T", "infinite=F", "integer=F", "irrational=F", "negative=F", "noninteger=F", "nonnegative=F", "nonpositive=F", "nonzero=F", "odd=F", "positive=F", "prime=F", "rational=F", "real=F", "zero=F"),
  "infinite=F" = c("finite=T"),
  "infinite=T" = c("algebraic=F", "complex=F", "composite=F", "even=F", "finite=F", "imaginary=F", "integer=F", "irrational=F", "negative=F", "nonnegative=F", "nonpositive=F", "nonzero=F", "odd=F", "positive=F", "prime=F", "rational=F", "real=F", "transcendental=F", "zero=F"),
  "integer=F" = c("composite=F", "even=F", "odd=F", "prime=F", "zero=F"),
  "integer=T" = c("algebraic=T", "commutative=T", "complex=T", "extended_real=T", "finite=T", "imaginary=F", "infinite=F", "irrational=F", "noninteger=F", "rational=T", "real=T", "transcendental=F"),
  "irrational=F" = c(),
  "irrational=T" = c("commutative=T", "complex=T", "composite=F", "even=F", "extended_nonzero=T", "extended_real=T", "finite=T", "imaginary=F", "infinite=F", "integer=F", "noninteger=T", "nonzero=T", "odd=F", "prime=F", "rational=F", "real=T", "zero=F"),
  "negative=F" = c(),
  "negative=T" = c("commutative=T", "complex=T", "composite=F", "extended_negative=T", "extended_nonnegative=F", "extended_nonpositive=T", "extended_nonzero=T", "extended_positive=F", "extended_real=T", "finite=T", "imaginary=F", "infinite=F", "nonnegative=F", "nonpositive=T", "nonzero=T", "positive=F", "prime=F", "real=T", "zero=F"),
  "noninteger=F" = c(),
  "noninteger=T" = c("commutative=T", "composite=F", "even=F", "extended_nonzero=T", "extended_real=T", "imaginary=F", "integer=F", "odd=F", "prime=F", "zero=F"),
  "nonnegative=F" = c("composite=F", "positive=F", "prime=F", "zero=F"),
  "nonnegative=T" = c("commutative=T", "complex=T", "extended_negative=F", "extended_nonnegative=T", "extended_real=T", "finite=T", "imaginary=F", "infinite=F", "negative=F", "real=T"),
  "nonpositive=F" = c("negative=F", "zero=F"),
  "nonpositive=T" = c("commutative=T", "complex=T", "composite=F", "extended_nonpositive=T", "extended_positive=F", "extended_real=T", "finite=T", "imaginary=F", "infinite=F", "positive=F", "prime=F", "real=T"),
  "nonzero=F" = c("composite=F", "negative=F", "positive=F", "prime=F"),
  "nonzero=T" = c("commutative=T", "complex=T", "extended_nonzero=T", "extended_real=T", "finite=T", "imaginary=F", "infinite=F", "real=T", "zero=F"),
  "odd=F" = c(),
  "odd=T" = c("algebraic=T", "commutative=T", "complex=T", "even=F", "extended_nonzero=T", "extended_real=T", "finite=T", "imaginary=F", "infinite=F", "integer=T", "irrational=F", "noninteger=F", "nonzero=T", "rational=T", "real=T", "transcendental=F", "zero=F"),
  "polar=T" = c(),
  "positive=F" = c("composite=F", "prime=F"),
  "positive=T" = c("commutative=T", "complex=T", "extended_negative=F", "extended_nonnegative=T", "extended_nonpositive=F", "extended_nonzero=T", "extended_positive=T", "extended_real=T", "finite=T", "imaginary=F", "infinite=F", "negative=F", "nonnegative=T", "nonpositive=F", "nonzero=T", "real=T", "zero=F"),
  "prime=F" = c(),
  "prime=T" = c("algebraic=T", "commutative=T", "complex=T", "composite=F", "extended_negative=F", "extended_nonnegative=T", "extended_nonpositive=F", "extended_nonzero=T", "extended_positive=T", "extended_real=T", "finite=T", "imaginary=F", "infinite=F", "integer=T", "irrational=F", "negative=F", "noninteger=F", "nonnegative=T", "nonpositive=F", "nonzero=T", "positive=T", "rational=T", "real=T", "transcendental=F", "zero=F"),
  "rational=F" = c("composite=F", "even=F", "integer=F", "odd=F", "prime=F", "zero=F"),
  "rational=T" = c("algebraic=T", "commutative=T", "complex=T", "extended_real=T", "finite=T", "imaginary=F", "infinite=F", "irrational=F", "real=T", "transcendental=F"),
  "real=F" = c("composite=F", "even=F", "integer=F", "irrational=F", "negative=F", "nonnegative=F", "nonpositive=F", "nonzero=F", "odd=F", "positive=F", "prime=F", "rational=F", "zero=F"),
  "real=T" = c("commutative=T", "complex=T", "extended_real=T", "finite=T", "imaginary=F", "infinite=F"),
  "transcendental=F" = c(),
  "transcendental=T" = c("algebraic=F", "commutative=T", "complex=T", "composite=F", "even=F", "finite=T", "infinite=F", "integer=F", "odd=F", "prime=F", "rational=F", "zero=F"),
  "zero=F" = c(),
  "zero=T" = c("algebraic=T", "commutative=T", "complex=T", "composite=F", "even=T", "extended_negative=F", "extended_nonnegative=T", "extended_nonpositive=T", "extended_nonzero=F", "extended_positive=F", "extended_real=T", "finite=T", "imaginary=F", "infinite=F", "integer=T", "irrational=F", "negative=F", "noninteger=F", "nonnegative=T", "nonpositive=T", "nonzero=F", "odd=F", "positive=F", "prime=F", "rational=T", "real=T", "transcendental=F")
)
.FR_BETA <- list(
  list(c("even=T", "positive=T", "prime=F"), "composite=T"),
  list(c("composite=F", "positive=T", "prime=F"), "even=F"),
  list(c("integer=T", "odd=F"), "even=T"),
  list(c("extended_positive=F", "extended_real=T", "zero=F"), "extended_negative=T"),
  list(c("extended_nonpositive=T", "extended_nonzero=T"), "extended_negative=T"),
  list(c("extended_negative=F", "extended_real=T"), "extended_nonnegative=T"),
  list(c("extended_positive=F", "extended_real=T"), "extended_nonpositive=T"),
  list(c("extended_real=T", "zero=F"), "extended_nonzero=T"),
  list(c("extended_negative=F", "extended_real=T", "zero=F"), "extended_positive=T"),
  list(c("extended_nonnegative=T", "extended_nonzero=T"), "extended_positive=T"),
  list(c("infinite=F", "real=F"), "extended_real=F"),
  list(c("extended_negative=F", "extended_positive=F", "zero=F"), "extended_real=F"),
  list(c("extended_real=T", "real=F"), "infinite=T"),
  list(c("rational=F", "real=T"), "irrational=T"),
  list(c("positive=F", "real=T", "zero=F"), "negative=T"),
  list(c("nonpositive=T", "nonzero=T"), "negative=T"),
  list(c("extended_negative=T", "finite=T"), "negative=T"),
  list(c("extended_real=T", "integer=F"), "noninteger=T"),
  list(c("negative=F", "real=T"), "nonnegative=T"),
  list(c("extended_nonnegative=T", "finite=T"), "nonnegative=T"),
  list(c("positive=F", "real=T"), "nonpositive=T"),
  list(c("extended_nonpositive=T", "finite=T"), "nonpositive=T"),
  list(c("extended_nonzero=T", "finite=T"), "nonzero=T"),
  list(c("even=F", "integer=T"), "odd=T"),
  list(c("composite=F", "even=T", "prime=F"), "positive=F"),
  list(c("negative=F", "real=T", "zero=F"), "positive=T"),
  list(c("nonnegative=T", "nonzero=T"), "positive=T"),
  list(c("extended_positive=T", "finite=T"), "positive=T"),
  list(c("composite=F", "even=T", "positive=T"), "prime=T"),
  list(c("negative=F", "positive=F", "zero=F"), "real=F"),
  list(c("extended_real=T", "infinite=F"), "real=T"),
  list(c("extended_real=T", "finite=T"), "real=T"),
  list(c("algebraic=F", "complex=T"), "transcendental=T"),
  list(c("extended_negative=F", "extended_positive=F", "extended_real=T"), "zero=T"),
  list(c("negative=F", "positive=F", "real=T"), "zero=T"),
  list(c("extended_nonnegative=T", "extended_nonpositive=T"), "zero=T"),
  list(c("nonnegative=T", "nonpositive=T"), "zero=T")
)
.FR_TRIG <- list(
  "algebraic=F" = c(33L, 12L, 4L, 9L, 30L, 15L, 26L, 14L, 18L, 8L),
  "algebraic=T" = c(11L, 31L, 32L, 28L, 17L, 22L, 20L, 23L),
  "commutative=F" = c(),
  "commutative=T" = c(),
  "complex=F" = c(11L, 13L, 12L, 4L, 9L, 18L, 8L),
  "complex=T" = c(33L, 11L, 31L, 32L, 28L, 17L, 22L, 20L, 23L),
  "composite=F" = c(2L, 29L, 25L),
  "composite=T" = c(24L, 3L),
  "even=F" = c(24L, 12L, 4L, 9L, 30L, 15L, 26L, 8L),
  "even=T" = c(4L, 34L, 9L, 7L, 6L, 15L, 35L, 26L, 21L, 19L, 28L, 17L, 22L, 20L, 23L, 1L, 29L, 25L, 8L),
  "extended_negative=F" = c(12L, 34L, 9L, 6L, 30L, 35L, 26L, 19L),
  "extended_negative=T" = c(31L, 13L, 32L, 30L, 15L, 21L, 17L, 22L, 23L, 18L),
  "extended_nonnegative=F" = c(12L, 4L, 7L, 30L, 15L, 21L, 8L),
  "extended_nonnegative=T" = c(31L, 13L, 32L, 34L, 9L, 10L, 7L, 30L, 35L, 26L, 19L, 20L, 36L, 18L, 8L),
  "extended_nonpositive=F" = c(12L, 9L, 6L, 30L, 26L, 19L, 8L),
  "extended_nonpositive=T" = c(31L, 13L, 32L, 4L, 34L, 5L, 6L, 30L, 15L, 35L, 21L, 22L, 36L, 18L, 8L),
  "extended_nonzero=F" = c(12L, 34L, 7L, 6L, 30L, 35L, 21L, 19L),
  "extended_nonzero=T" = c(31L, 13L, 32L, 4L, 9L, 5L, 10L, 7L, 6L, 30L, 15L, 26L, 23L, 18L),
  "extended_positive=F" = c(12L, 4L, 34L, 7L, 30L, 15L, 35L, 21L),
  "extended_positive=T" = c(31L, 13L, 32L, 30L, 26L, 19L, 28L, 20L, 23L, 18L),
  "extended_real=F" = c(),
  "extended_real=T" = c(31L, 13L, 32L, 4L, 34L, 9L, 7L, 6L, 18L, 8L),
  "finite=F" = c(12L, 4L, 9L, 18L, 8L),
  "finite=T" = c(11L, 31L, 32L, 28L, 17L, 22L, 20L, 23L),
  "imaginary=T" = c(33L),
  "infinite=F" = c(11L, 31L, 32L, 28L, 17L, 22L, 20L, 23L),
  "infinite=T" = c(12L, 4L, 9L, 18L, 8L),
  "integer=F" = c(12L, 4L, 9L, 30L, 15L, 26L, 18L, 8L),
  "integer=T" = c(24L, 3L, 4L, 34L, 9L, 7L, 6L, 15L, 35L, 26L, 21L, 19L, 28L, 17L, 22L, 20L, 23L, 8L),
  "irrational=T" = c(33L, 4L, 9L, 5L, 10L, 7L, 6L, 15L, 26L, 16L, 27L, 21L, 19L, 28L, 17L, 22L, 20L),
  "negative=F" = c(30L, 35L, 26L, 19L),
  "negative=T" = c(33L, 14L, 18L),
  "noninteger=T" = c(31L, 13L, 32L, 4L, 9L, 5L, 10L, 7L, 6L, 30L, 15L, 26L, 23L),
  "nonnegative=F" = c(12L, 4L, 9L, 30L, 15L, 21L, 8L),
  "nonnegative=T" = c(33L, 34L, 9L, 10L, 7L, 35L, 26L, 27L, 21L, 28L, 22L, 23L, 36L, 37L, 14L, 18L, 8L),
  "nonpositive=F" = c(12L, 4L, 9L, 30L, 26L, 19L, 8L),
  "nonpositive=T" = c(33L, 4L, 34L, 5L, 6L, 15L, 35L, 16L, 19L, 17L, 20L, 23L, 36L, 37L, 14L, 18L, 8L),
  "nonzero=F" = c(30L, 35L, 21L, 19L),
  "nonzero=T" = c(33L, 4L, 9L, 5L, 10L, 7L, 6L, 15L, 26L, 16L, 27L, 21L, 19L, 28L, 17L, 22L, 20L, 14L, 18L),
  "odd=F" = c(3L),
  "odd=T" = c(4L, 9L, 5L, 10L, 7L, 6L, 15L, 26L, 16L, 27L, 21L, 19L, 28L, 17L, 22L, 20L),
  "polar=T" = c(),
  "positive=F" = c(30L, 15L, 35L, 21L),
  "positive=T" = c(33L, 1L, 2L, 29L, 14L, 18L),
  "prime=F" = c(1L, 2L, 25L),
  "prime=T" = c(24L, 3L),
  "rational=F" = c(12L, 4L, 9L, 30L, 15L, 26L, 14L, 18L, 8L),
  "rational=T" = c(4L, 34L, 9L, 7L, 6L, 15L, 35L, 26L, 21L, 19L, 28L, 17L, 22L, 20L, 23L, 18L, 8L),
  "real=F" = c(11L, 13L, 12L, 4L, 9L, 18L, 8L),
  "real=T" = c(33L, 4L, 34L, 9L, 7L, 6L, 15L, 35L, 26L, 21L, 19L, 28L, 17L, 22L, 20L, 23L, 14L, 18L, 8L),
  "transcendental=T" = c(11L, 31L, 32L, 12L, 4L, 9L, 30L, 15L, 26L, 28L, 17L, 22L, 20L, 23L, 14L, 18L, 8L),
  "zero=F" = c(12L, 4L, 9L, 30L, 15L, 26L, 8L),
  "zero=T" = c()
)

.fr_parse <- function(l) list(sub("=.$", "", l), endsWith(l, "T"))
.FR_IMPL_P <- lapply(.FR_IMPL, function(v) { v <- as.character(v); list(sub("=.$", "", v), endsWith(v, "T")) })
.FR_BETA_P <- lapply(.FR_BETA, function(r) list(sub("=.$", "", r[[1L]]), endsWith(r[[1L]], "T"), r[[2L]]))

# Closure of a named logical vector of facts. Facts outside SymPy's rules
# (hermitian, antihermitian, polar) are kept; real implies hermitian and
# imaginary antihermitian.
.close_facts <- function(v, strict = TRUE) {
  full <- stats::setNames(rep(NA, length(.FACTS)), .FACTS)
  given <- v[!is.na(v) & names(v) %in% .FACTS]
  tell <- function(k, val) {
    cur <- full[[k]]
    if (!is.na(cur)) {
      if (cur == val) return(FALSE)
      if (strict) stop("inconsistent assumptions: ", k)
      return(FALSE)
    }
    full[[k]] <<- val
    TRUE
  }
  if (identical(unname(given["hermitian"]), FALSE) && is.na(given["real"])) given[["real"]] <- FALSE
  if (identical(unname(given["antihermitian"]), FALSE) && is.na(given["imaginary"])) given[["imaginary"]] <- FALSE
  facts <- paste0(names(given), ifelse(given, "=T", "=F"))
  while (length(facts)) {
    may <- integer()
    for (lv in facts) {
      kv <- .fr_parse(lv)
      if (!tell(kv[[1L]], kv[[2L]])) next
      imp <- .FR_IMPL_P[[lv]]
      if (!is.null(imp)) for (i in seq_along(imp[[1L]])) tell(imp[[1L]][i], imp[[2L]][i])
      may <- union(may, .FR_TRIG[[lv]])
    }
    facts <- character()
    for (bi in may) {
      br <- .FR_BETA_P[[bi]]
      cur <- full[br[[1L]]]
      if (!anyNA(cur) && all(cur == br[[2L]])) facts <- c(facts, br[[3L]])
    }
  }
  if (isTRUE(full[["real"]]) && is.na(full[["hermitian"]])) full[["hermitian"]] <- TRUE
  if (isTRUE(full[["imaginary"]]) && is.na(full[["antihermitian"]])) full[["antihermitian"]] <- TRUE
  full
}

.symbol_facts <- function(given) {
  if (!length(given)) {
    f <- stats::setNames(rep(NA, length(.FACTS)), .FACTS)
    f[["commutative"]] <- TRUE
    return(f)
  }
  if (is.na(given["commutative"])) given["commutative"] <- TRUE
  .close_facts(given)
}


#' Query assumptions about an expression
#'
#' `ask_fact()` returns `TRUE`, `FALSE` or `NA` (unknown) for a single
#' property, following SymPy's `expr.is_<property>` attributes. The helpers
#' `is_positive()`, `is_real()` and friends are shortcuts.
#'
#' @param e A `sym` expression or a value convertible with [as_sym()].
#' @param fact One of the assumption names listed in [Symbol()].
#' @return `TRUE`, `FALSE` or `NA`.
#' @examples
#' x <- Symbol("x", positive = TRUE)
#' ask_fact(x^2 + 1, "positive")
#' is_integer(Symbol("n", integer = TRUE) * 2)
#' @export
ask_fact <- function(e, fact) {
  e <- as_sym(e)
  # is_polar: exp_polar, polar_lift, polar symbols and their products and powers
  if (identical(fact, "polar") && !(e$type %in% c("Symbol", "Dummy", "Wild")) && .is_polar(e)) return(TRUE)
  f <- .facts(e)
  r <- f[fact]
  if (is.null(r) || length(r) != 1L) stop("unknown fact: ", fact)
  unname(r)
}

.facts <- function(e) {
  t <- e$type
  if (t == "Symbol" || t == "Dummy" || t == "Wild") return(e$extra$facts)
  k <- e$key
  cacheable <- nchar(k) < 9000L
  if (cacheable) {
    hit <- .fact_cache[[k]]
    if (!is.null(hit)) return(hit)
  }
  v <- switch(t,
              Integer = , Rational = .number_facts(e),
              Float = .float_facts(e),
              Constant = .constant_facts(e$name),
              Add = .add_facts(e),
              Mul = .mul_facts(e),
              Pow = .pow_facts(e),
              Function = .function_facts(e),
              BooleanAtom = c(commutative = TRUE),
              .other_facts(e))
  full <- tryCatch(.close_facts(v), error = function(err) .close_facts(v, strict = FALSE))
  if (cacheable) {
    if (length(.fact_cache) > 50000L) rm(list = ls(.fact_cache, all.names = TRUE), envir = .fact_cache)
    assign(k, full, envir = .fact_cache)
  }
  full
}

.other_facts <- function(e) {
  r <- .facts_hooks[[e$type]]
  if (!is.null(r)) return(r(e))
  c(commutative = TRUE)
}


.number_facts <- function(e) {
  v <- e$value
  s <- if (v > 0) 1 else if (v < 0) -1 else 0
  isint <- e$type == "Integer"
  f <- c(commutative = TRUE, rational = TRUE, real = TRUE, finite = TRUE,
         algebraic = TRUE, integer = isint, positive = s > 0,
         negative = s < 0, zero = s == 0, irrational = FALSE,
         imaginary = FALSE, transcendental = FALSE, infinite = FALSE)
  if (isint) {
    ev <- gmp::mod.bigz(gmp::numerator(v), 2) == 0
    f["even"] <- ev
    f["odd"] <- !ev
    n <- gmp::numerator(v)
    if (n > 1) {
      pr <- gmp::isprime(n) > 0
      f["prime"] <- pr
      f["composite"] <- !pr
    } else {
      f["prime"] <- FALSE
      f["composite"] <- FALSE
    }
  } else {
    f["even"] <- FALSE
    f["odd"] <- FALSE
    f["prime"] <- FALSE
    f["composite"] <- FALSE
  }
  f
}

.float_facts <- function(e) {
  v <- e$value
  if (is.nan(v)) return(c(commutative = TRUE))
  if (is.infinite(v)) return(c(commutative = TRUE, finite = FALSE,
                               infinite = TRUE, extended_real = TRUE))
  c(commutative = TRUE, real = TRUE, finite = TRUE, positive = v > 0,
    negative = v < 0, zero = v == 0, imaginary = FALSE, infinite = FALSE)
}

.constant_facts <- function(nm) {
  switch(EXPR = nm,
         pi = , E = c(commutative = TRUE, positive = TRUE, irrational = TRUE,
                      transcendental = TRUE),
         EulerGamma = , Catalan = c(commutative = TRUE, positive = TRUE),
         GoldenRatio = , TribonacciConstant = c(commutative = TRUE,
                                                positive = TRUE,
                                                irrational = TRUE,
                                                algebraic = TRUE),
         I = c(commutative = TRUE, imaginary = TRUE, algebraic = TRUE,
               finite = TRUE),
         oo = c(commutative = TRUE, finite = FALSE, infinite = TRUE,
                extended_real = TRUE, real = FALSE, complex = FALSE,
                extended_positive = TRUE),
         `-oo` = c(commutative = TRUE, finite = FALSE, infinite = TRUE,
                   extended_real = TRUE, real = FALSE, complex = FALSE,
                   extended_negative = TRUE),
         zoo = c(commutative = TRUE, finite = FALSE, infinite = TRUE,
                 real = FALSE, complex = FALSE, extended_real = FALSE),
         nan = c(commutative = TRUE),
         c(commutative = TRUE))
}

.args_facts <- function(e) lapply(e$args, .facts)
.all_true <- function(fs, nm) all(vapply(fs, function(f) isTRUE(f[[nm]]), TRUE))
.any_true <- function(fs, nm) any(vapply(fs, function(f) isTRUE(f[[nm]]), TRUE))
.count_false <- function(fs, nm) sum(vapply(fs, function(f) identical(f[[nm]], FALSE), TRUE))

.add_facts <- function(e) {
  fs <- .args_facts(e)
  out <- c(commutative = .all_true(fs, "commutative"))
  if (.all_true(fs, "finite")) out["finite"] <- TRUE
  else if (sum(vapply(fs, function(f) isTRUE(f[["infinite"]]), TRUE)) == 1L &&
           .count_false(fs, "finite") == 1L) out["finite"] <- FALSE
  for (nm in c("real", "integer", "rational", "algebraic", "complex")) {
    if (.all_true(fs, nm)) out[nm] <- TRUE
    else if (.count_false(fs, nm) == 1L && sum(vapply(fs, function(f)
      isTRUE(f[[nm]]), TRUE)) == length(fs) - 1L &&
      .all_true(fs, "finite")) {
      bad <- fs[[which(vapply(fs, function(f) identical(f[[nm]], FALSE), TRUE))]]
      if (nm == "real" && isTRUE(bad[["complex"]])) out[nm] <- FALSE
      if (nm %in% c("integer", "rational") && isTRUE(bad[["real"]]) &&
          isTRUE(bad[["finite"]])) out[nm] <- FALSE
    }
  }
  if (.all_true(fs, "nonnegative") && .any_true(fs, "positive")) out["positive"] <- TRUE
  if (.all_true(fs, "nonpositive") && .any_true(fs, "negative")) out["negative"] <- TRUE
  if (.all_true(fs, "nonnegative")) out["nonnegative"] <- TRUE
  if (.all_true(fs, "nonpositive")) out["nonpositive"] <- TRUE
  if (isTRUE(out["integer"])) {
    par <- vapply(fs, function(f) if (isTRUE(f[["even"]])) 0L
                  else if (isTRUE(f[["odd"]])) 1L else NA_integer_, 0L)
    if (!anyNA(par)) {
      odd <- sum(par) %% 2L == 1L
      out["odd"] <- odd
      out["even"] <- !odd
    }
  }
  if (.all_true(fs, "real") && .all_true(fs, "finite") &&
      is.na(out["positive"]) && is.na(out["negative"])) {
    sg <- .numeric_sign(e)
    # no significant digit: SymPy proves zero by the minimal polynomial
    if (is.na(sg) && isTRUE(.root_bound_zero(e))) sg <- 0
    if (!is.na(sg)) {
      out["positive"] <- sg > 0
      out["negative"] <- sg < 0
      out["zero"] <- sg == 0
    }
  }
  # Add._eval_is_extended_positive and friends: the smallest value of the
  # non-constant part (_monotonic_sign) plus the constant
  if (is.na(out["positive"]) || is.na(out["negative"]) || is.na(out["nonnegative"]) || is.na(out["nonpositive"])) {
    ca <- .coeff_add_split(e); c <- ca[[1L]]; a <- ca[[2L]]
    if (!is_zero_num(c) && length(free_symbols(e))) {
      v <- .monotonic_sign(a)
      if (!is.null(v)) {
        sv <- .add(list(v, c))
        if (!identical(sv$key, e$key)) {
          if (isTRUE(ask_fact(a, "extended_nonnegative"))) {
            if (isTRUE(ask_fact(sv, "extended_positive"))) out["positive"] <- TRUE
            if (isTRUE(ask_fact(sv, "extended_nonnegative"))) out["nonnegative"] <- TRUE
          }
          if (isTRUE(ask_fact(a, "extended_nonpositive"))) {
            if (isTRUE(ask_fact(sv, "extended_negative"))) out["negative"] <- TRUE
            if (isTRUE(ask_fact(sv, "extended_nonpositive"))) out["nonpositive"] <- TRUE
          }
        }
      }
    }
  }
  if (isTRUE(out["rational"]) || is.na(out["rational"])) {
    nirr <- sum(vapply(fs, function(f) isTRUE(f[["irrational"]]), TRUE))
    if (nirr == 1L && sum(vapply(fs, function(f) isTRUE(f[["rational"]]), TRUE)) ==
        length(fs) - 1L) out["irrational"] <- TRUE
  }
  out
}

# exprtools._monotonic_sign: the smallest (or largest) value an expression
# can take, judged from the assumptions of its symbols, or NULL.
# ponytail: the univariate polynomial and exp/log branches return NULL
.MS_EPS <- Dummy("_eps", positive = TRUE)
.monotonic_sign <- function(e) {
  if (!isTRUE(ask_fact(e, "extended_real"))) return(NULL)
  ne <- .neg(e)
  if (ne$type %in% c("Symbol", "Dummy")) {
    rv <- .monotonic_sign(ne)
    return(if (is.null(rv)) NULL else .neg(rv))
  }
  if (e$type != "Add" && !length(free_symbols(.as_numer_denom(e)[[2L]]))) {
    f <- function(nm) ask_fact(e, nm)
    if (isTRUE(f("prime"))) return(if (isTRUE(f("odd"))) .int(3L) else S_Two)
    if (isTRUE(f("composite"))) return(if (isTRUE(f("odd"))) .int(9L) else .int(4L))
    if (isTRUE(f("positive"))) {
      if (isTRUE(f("even"))) return(if (isFALSE(f("prime"))) .int(4L) else S_Two)
      if (isTRUE(f("integer"))) return(S_One)
      return(.MS_EPS)
    }
    if (isTRUE(f("extended_negative"))) {
      if (isTRUE(f("even"))) return(.int(-2L))
      if (isTRUE(f("integer"))) return(S_NegativeOne)
      return(.neg(.MS_EPS))
    }
    if (isTRUE(f("zero")) || isTRUE(f("extended_nonpositive")) || isTRUE(f("extended_nonnegative"))) return(S_Zero)
    return(NULL)
  }
  if (length(free_symbols(e)) == 1L) return(NULL)
  ca <- .coeff_add_split(e); c <- ca[[1L]]; a <- ca[[2L]]
  v <- NULL
  pows <- Filter(function(p) length(free_symbols(p)), atoms(a, "Pow"))
  if (!isTRUE(is_polynomial(a))) {
    nd <- .as_numer_denom(a)
    if (!(!length(free_symbols(nd[[1L]])) || !length(free_symbols(nd[[2L]])))) return(NULL)
    if (a$type %in% c("Mul", "Pow") && isTRUE(ask_fact(a, "rational")) &&
        all(vapply(atoms(a, "Pow"), function(p) p$args[[2L]]$type == "Integer", TRUE)) &&
        (isTRUE(ask_fact(a, "positive")) || isTRUE(ask_fact(a, "negative")))) {
      v <- S_One
      for (ai in .margs(a)) {
        if (!length(free_symbols(ai))) { v <- .mul(list(v, ai)); next }
        reps <- list()
        for (x in free_symbols(ai)) {
          mx <- .monotonic_sign(x)
          if (is.null(mx)) return(NULL)
          reps[[length(reps) + 1L]] <- list(x, mx)
        }
        v <- .mul(list(v, subs(ai, reps)))
      }
    }
  } else if ((!is_zero_num(c) || a$type == "Add") && !length(pows)) {
    signed <- if (!is_zero_num(c)) isTRUE(ask_fact(a, "nonpositive")) || isTRUE(ask_fact(a, "nonnegative"))
              else { ts <- .aargs(a); all(vapply(ts, function(t) isTRUE(ask_fact(t, "nonnegative")), TRUE)) ||
                     all(vapply(ts, function(t) isTRUE(ask_fact(t, "nonpositive")), TRUE)) }
    if (signed) {
      reps <- list()
      for (x in free_symbols(a)) {
        mx <- .monotonic_sign(x)
        if (is.null(mx)) return(NULL)
        reps[[length(reps) + 1L]] <- list(x, mx)
      }
      v <- xreplace(a, reps)
    }
  }
  if (is.null(v)) return(NULL)
  rv <- .add(list(v, c))
  if (isTRUE(ask_fact(a, "nonnegative"))) {
    if (isTRUE(ask_fact(rv, "positive"))) {
      rv <- subs(rv, .MS_EPS, S_Zero)
      return(if (!is_zero_num(rv)) rv else Dummy("pos", positive = TRUE))
    }
    if (isTRUE(ask_fact(rv, "zero"))) return(Dummy("nneg", nonnegative = TRUE))
  }
  if (isTRUE(ask_fact(a, "nonpositive"))) {
    if (isTRUE(ask_fact(rv, "negative"))) {
      rv <- subs(rv, .MS_EPS, S_Zero)
      return(if (!is_zero_num(rv)) rv else Dummy("neg", negative = TRUE))
    }
    if (isTRUE(ask_fact(rv, "zero"))) return(Dummy("npos", nonpositive = TRUE))
  }
  NULL
}

# Sign of an expression without free symbols, as Add._eval_is_extended_positive
# on a number: evalf(2) at increasing precision (see .mp_num_sign), so a
# cancellation anywhere inside the expression cannot flip the sign.
.numeric_sign <- function(e) {
  if (length(free_symbols(e))) return(NA)
  .mp_num_sign(e)
}

# Expr._eval_is_extended_positive_negative when evaluation finds no
# significant digit: SymPy proves zero for a function free algebraic number
# by its minimal polynomial being the variable; R's minimal_polynomial is far
# too slow for that (degree 81 for nested cube roots), so zero is proved
# instead by the BFMSS root bound, which a nonzero value lies above.
#
# The BFMSS separation bound for a number built from rationals, I, sums,
# products and rational powers: with u and l bounding the conjugates of
# the algebraic integers U and L in e = U/L and D a degree bound (the
# product of the distinct root indices), a nonzero e has |e| >= 1/(u^(D-1)*l).
# TRUE when evaluation proves e = 0, NA when the bound does not apply.
# ponytail: gives up when the bound passes 20000 bits.
.root_bound_zero <- function(e) {
  rads <- character(); qs <- numeric(); hasI <- FALSE; scale <- 0
  lg <- function(n) if (n == 0) 0 else as.numeric(log2(abs(gmp::as.bigz(n))))
  ladd <- function(a, b) max(a, b) + log2(1 + 2^(min(a, b) - max(a, b)))
  rec <- function(z) {
    t <- z$type
    r <- if (t %in% c("Integer", "Rational")) c(lg(gmp::numerator(z$value)), lg(gmp::denominator(z$value)))
    else if (is_const(z, "I")) { hasI <<- TRUE; c(0, 0) }
    else if (t == "Add" || t == "Mul") {
      acc <- rec(z$args[[1L]])
      for (b in z$args[-1L]) {
        rb <- rec(b)
        acc <- if (t == "Add") c(ladd(acc[1L] + rb[2L], rb[1L] + acc[2L]), acc[2L] + rb[2L]) else acc + rb
      }
      acc
    } else if (t == "Pow" && z$args[[2L]]$type %in% c("Integer", "Rational")) {
      q <- z$args[[2L]]$value
      b <- rec(z$args[[1L]])
      if (gmp::numerator(q) < 0) b <- rev(b)
      b <- b * abs(as.numeric(gmp::numerator(q)))
      k <- as.numeric(gmp::denominator(q))
      if (k > 1) {
        rk <- paste0(z$args[[1L]]$key, "|", k)
        if (!(rk %in% rads)) { rads <<- c(rads, rk); qs <<- c(qs, k) }
        # one radical: (U/L)^(1/k) = (U*L^(k-1))^(1/k)/L or U/(U^(k-1)*L)^(1/k)
        b <- if (b[1L] >= b[2L]) c((b[1L] + (k - 1) * b[2L]) / k, b[2L]) else c(b[1L], ((k - 1) * b[1L] + b[2L]) / k)
      }
      b
    } else stop("not algebraic")
    scale <<- max(scale, sum(r))
    r
  }
  r <- tryCatch(rec(e), error = function(err) NULL)
  if (is.null(r)) return(NA)
  D <- prod(qs) * (if (hasI) 2 else 1)
  bound <- (D - 1) * r[1L] + r[2L]
  if (!is.finite(bound) || bound > 20000) return(NA)
  v <- tryCatch(.mpeval(e, ceiling(bound + scale) + 64L, list()), error = function(err) NULL)
  if (is.null(v)) return(NA)
  max(.mp_mag(v$re), .mp_mag(v$im)) < -bound
}

.mul_facts <- function(e) {
  fs <- .args_facts(e)
  out <- c(commutative = .all_true(fs, "commutative"))
  allfin <- .all_true(fs, "finite")
  if (allfin) out["finite"] <- TRUE
  # Mul._eval_is_infinite: an infinite factor and no factor that may be zero
  if (.any_true(fs, "infinite") && all(vapply(fs, function(f) isTRUE(f[["infinite"]]) ||
      identical(f[["zero"]], FALSE) || isTRUE(f[["nonzero"]]), TRUE))) {
    out["infinite"] <- TRUE; out["finite"] <- FALSE
  }
  if (.any_true(fs, "zero") && allfin) out["zero"] <- TRUE
  else if (.all_true(fs, "nonzero") && allfin) out["zero"] <- FALSE
  else if (all(vapply(fs, function(f) identical(f[["zero"]], FALSE), TRUE)) && allfin)
    out["zero"] <- FALSE
  if (.all_true(fs, "real")) out["real"] <- TRUE
  if (.all_true(fs, "integer")) out["integer"] <- TRUE
  if (.all_true(fs, "rational")) out["rational"] <- TRUE
  if (.all_true(fs, "algebraic")) out["algebraic"] <- TRUE
  if (.all_true(fs, "complex")) out["complex"] <- TRUE
  nim <- sum(vapply(fs, function(f) isTRUE(f[["imaginary"]]), TRUE))
  nreal <- sum(vapply(fs, function(f) isTRUE(f[["real"]]), TRUE))
  if (nim + nreal == length(fs) && nim > 0) {
    if (nim %% 2L == 0L) out["real"] <- TRUE
    else if (all(vapply(fs, function(f) isTRUE(f[["nonzero"]]) ||
                        isTRUE(f[["imaginary"]]), TRUE))) {
      out["imaginary"] <- TRUE
    }
  }
  sg <- vapply(fs, function(f) if (isTRUE(f[["positive"]])) 1L
               else if (isTRUE(f[["negative"]])) -1L else NA_integer_, 0L)
  if (!anyNA(sg)) {
    p <- prod(sg)
    out["positive"] <- p > 0
    out["negative"] <- p < 0
  } else if (.all_true(fs, "real")) {
    nn <- vapply(fs, function(f) if (isTRUE(f[["nonnegative"]])) 1L
                 else if (isTRUE(f[["nonpositive"]])) -1L else NA_integer_, 0L)
    if (!anyNA(nn)) {
      if (prod(nn) > 0) out["nonnegative"] <- TRUE else out["nonpositive"] <- TRUE
    }
  }
  if (isTRUE(out["integer"])) {
    if (.any_true(fs, "even")) out["even"] <- TRUE
    else if (.all_true(fs, "odd")) out["odd"] <- TRUE
  }
  nrat <- sum(vapply(fs, function(f) isTRUE(f[["rational"]]) &&
                       isTRUE(f[["nonzero"]]), TRUE))
  nirr <- sum(vapply(fs, function(f) isTRUE(f[["irrational"]]), TRUE))
  if (nirr == 1L && nrat == length(fs) - 1L) out["irrational"] <- TRUE
  # Mul._eval_is_extended_real, _eval_is_extended_positive and
  # _eval_is_extended_negative, for factors that may be infinite
  set <- function(nm, v) if (!is.na(v) && is.na(out[nm] %||% NA)) out[nm] <<- v
  set("extended_real", .mul_real_imag(fs, TRUE))
  set("extended_positive", .mul_pos_neg(fs, 1L))
  if (is_const(e$args[[1L]], "NegativeOne") || (e$args[[1L]]$type == "Integer" && e$args[[1L]]$value == -1)) {
    rest <- fs[-1L]
    set("extended_negative", if (length(rest) == 1L) rest[[1L]][["extended_positive"]] else .mul_pos_neg(rest, 1L))
  } else set("extended_negative", .mul_pos_neg(fs, -1L))
  out
}

# Mul._eval_real_imag
.mul_real_imag <- function(fs, real) {
  zero <- FALSE
  t_not <- NULL
  allfin <- function() all(vapply(fs, function(f) isTRUE(f[["finite"]]), NA))
  for (f in fs) {
    ci <- if (isTRUE(f[["complex"]])) TRUE else f[["infinite"]]
    if (identical(ci, FALSE) && isFALSE(f[["extended_real"]])) return(FALSE)
    if (isTRUE(f[["imaginary"]])) {
      real <- !real
    } else if (isTRUE(f[["extended_real"]])) {
      if (!isTRUE(zero)) {
        z <- f[["zero"]]
        if (!isTRUE(z) && identical(zero, FALSE)) zero <- z
        else if (isTRUE(z)) return(if (allfin()) TRUE else NA)
      }
    } else if (isFALSE(f[["extended_real"]]) || isFALSE(f[["imaginary"]])) {
      if (!is.null(t_not)) return(NA)
      t_not <- f
    } else return(NA)
  }
  if (!is.null(t_not)) {
    if (isFALSE(t_not[["extended_real"]]) && real) return(zero)
    if (isFALSE(t_not[["imaginary"]]) && !real) return(zero)
    return(NA)
  }
  if (identical(zero, FALSE)) return(real)
  if (real) return(real)
  NA
}

# Mul._eval_pos_neg
.mul_pos_neg <- function(fs, sign) {
  saw_NON <- saw_NOT <- FALSE
  for (f in fs) {
    if (isTRUE(f[["extended_positive"]])) next
    if (isTRUE(f[["extended_negative"]])) sign <- -sign
    else if (isTRUE(f[["zero"]])) return(if (all(vapply(fs, function(g) isTRUE(g[["finite"]]), NA))) FALSE else NA)
    else if (isTRUE(f[["extended_nonpositive"]])) { sign <- -sign; saw_NON <- TRUE }
    else if (isTRUE(f[["extended_nonnegative"]])) saw_NON <- TRUE
    else if (isFALSE(f[["positive"]])) { sign <- -sign; if (saw_NOT) return(NA); saw_NOT <- TRUE }
    else if (isFALSE(f[["negative"]])) { if (saw_NOT) return(NA); saw_NOT <- TRUE }
    else return(NA)
  }
  if (sign == 1L && !saw_NON && !saw_NOT) return(TRUE)
  if (sign < 0L) return(FALSE)
  NA
}

.pow_facts <- function(e) {
  b <- e$args[[1L]]
  x <- e$args[[2L]]
  fb <- .facts(b)
  fx <- .facts(x)
  out <- c(commutative = isTRUE(fb[["commutative"]]) && isTRUE(fx[["commutative"]]))
  xint <- isTRUE(fx[["integer"]])
  xeven <- isTRUE(fx[["even"]])
  xodd <- isTRUE(fx[["odd"]])
  if (isTRUE(fb[["positive"]]) && isTRUE(fx[["real"]])) {
    out["positive"] <- TRUE
  }
  if (isTRUE(fb[["real"]]) && xint) {
    if (isTRUE(fb[["nonzero"]]) || isTRUE(fx[["nonnegative"]])) out["real"] <- TRUE
    if (xeven && isTRUE(fb[["nonzero"]])) out["positive"] <- TRUE
    if (xeven && (isTRUE(fb[["nonzero"]]) || isTRUE(fx[["nonnegative"]])))
      out["nonnegative"] <- TRUE
    if (xodd && isTRUE(fb[["negative"]])) out["negative"] <- TRUE
    if (xodd && isTRUE(fb[["nonnegative"]])) out["nonnegative"] <- TRUE
  }
  if (isTRUE(fb[["nonnegative"]]) && isTRUE(fx[["real"]]) &&
      isTRUE(fx[["positive"]])) out["nonnegative"] <- TRUE
  if (isTRUE(fb[["negative"]]) && x$type == "Rational" &&
      gmp::mod.bigz(gmp::denominator(x$value), 2) == 0) out["real"] <- FALSE
  if (isTRUE(fb[["integer"]]) && xint && isTRUE(fx[["nonnegative"]])) {
    out["integer"] <- TRUE
    if (isTRUE(fx[["positive"]])) {
      if (isTRUE(fb[["even"]])) out["even"] <- TRUE
      if (isTRUE(fb[["odd"]])) out["odd"] <- TRUE
    }
  }
  if (isTRUE(fb[["rational"]]) && xint &&
      (isTRUE(fb[["nonzero"]]) || isTRUE(fx[["nonnegative"]]))) out["rational"] <- TRUE
  if (is_exact_number(b) && b$value > 0 && x$type == "Rational") {
    out["algebraic"] <- TRUE
    if (abs(gmp::numerator(b$value)) < gmp::as.bigz("1000000000000000000000000000000") &&
        gmp::denominator(b$value) < gmp::as.bigz("1000000000000000000000000000000"))
      out["irrational"] <- TRUE
  }
  if (isTRUE(fb[["algebraic"]]) && isTRUE(fx[["rational"]]) &&
      (isTRUE(fb[["nonzero"]]) || isTRUE(fx[["positive"]]))) out["algebraic"] <- TRUE
  if (isTRUE(fb[["finite"]]) && isTRUE(fx[["finite"]]) &&
      (isTRUE(fb[["nonzero"]]) || isTRUE(fx[["nonnegative"]]))) {
    out["finite"] <- TRUE
    if (isTRUE(fb[["complex"]]) && isTRUE(fx[["complex"]])) out["complex"] <- TRUE
  }
  if (isTRUE(fb[["nonzero"]]) && isTRUE(fb[["finite"]]) && isTRUE(fx[["finite"]]))
    out["zero"] <- FALSE
  if (isTRUE(fb[["zero"]]) && isTRUE(fx[["positive"]])) out["zero"] <- TRUE
  set <- function(nm, v) if (length(v) == 1L && !is.na(v) && is.na(out[nm] %||% NA)) out[nm] <<- v
  set("extended_positive", .pow_ext_positive(b, x, fb, fx))
  set("extended_negative", .pow_ext_negative(b, x, fb, fx))
  set("extended_real", .pow_ext_real(b, x, fb, fx))
  out
}

# Pow._eval_is_extended_positive
.pow_ext_positive <- function(b, x, fb, fx) {
  if (identical(b$key, x$key)) return(if (isTRUE(fb[["extended_nonnegative"]])) TRUE else NA)
  if (isTRUE(fb[["positive"]])) return(if (isTRUE(fx[["real"]])) TRUE else NA)
  if (isTRUE(fb[["extended_negative"]])) return(if (isTRUE(fx[["even"]])) TRUE else if (isTRUE(fx[["odd"]])) FALSE else NA)
  if (isTRUE(fb[["zero"]])) return(if (isTRUE(fx[["extended_real"]])) fx[["zero"]] else NA)
  if (isTRUE(fb[["extended_nonpositive"]])) return(if (isTRUE(fx[["odd"]])) FALSE else NA)
  if (isTRUE(fb[["imaginary"]]) && x$type == "Integer") return(gmp::mod.bigz(gmp::numerator(x$value), 4) == 0)
  NA
}
# Pow._eval_is_extended_negative
.pow_ext_negative <- function(b, x, fb, fx) {
  if (identical(x$key, S_Half$key) && (isTRUE(fb[["complex"]]) || isTRUE(fb[["extended_real"]]))) return(FALSE)
  if (isTRUE(fb[["extended_negative"]])) {
    if (isTRUE(fx[["odd"]]) && isTRUE(fb[["finite"]])) return(TRUE)
    if (isTRUE(fx[["even"]])) return(FALSE)
  } else if (isTRUE(fb[["extended_positive"]]) || isTRUE(fb[["zero"]])) {
    if (isTRUE(fx[["extended_real"]])) return(FALSE)
  } else if (isTRUE(fb[["extended_nonnegative"]])) {
    if (isTRUE(fx[["extended_nonnegative"]])) return(FALSE)
  } else if (isTRUE(fb[["extended_nonpositive"]]) || isTRUE(fb[["extended_real"]])) {
    if (isTRUE(fx[["even"]])) return(FALSE)
  }
  NA
}
# Pow._eval_is_extended_real, the real base and exponent cases
.pow_ext_real <- function(b, x, fb, fx) {
  if (is_const(b, "E")) return(if (isTRUE(fx[["extended_real"]])) TRUE else NA)
  real_b <- fb[["extended_real"]]; real_e <- fx[["extended_real"]]
  if (is.na(real_b) || is.na(real_e)) return(NA)
  if (real_b && real_e) {
    if (isTRUE(fb[["extended_positive"]])) return(TRUE)
    if (isTRUE(fb[["extended_nonnegative"]]) && isTRUE(fx[["extended_nonnegative"]])) return(TRUE)
    if (isTRUE(fx[["integer"]]) && isTRUE(fb[["extended_nonzero"]])) return(TRUE)
    if (isTRUE(fx[["integer"]]) && isTRUE(fx[["nonnegative"]])) return(TRUE)
    if (isTRUE(fb[["extended_negative"]]) && x$type %in% c("Integer", "Rational")) return(FALSE)
  }
  if (real_e && isTRUE(fx[["extended_negative"]]) && isFALSE(fb[["zero"]]))
    return(ask_fact(.pow(b, .neg(x)), "extended_real"))
  NA
}

.function_facts <- function(e) {
  d <- .FUNS[[e$name]]
  if (!is.null(d) && !is.null(d$facts)) {
    r <- d$facts(e$args, lapply(e$args, .facts))
    if (is.null(r)) r <- logical()
    # Function._eval_is_commutative: all arguments commutative
    r["commutative"] <- .all_true(lapply(e$args, .facts), "commutative")
  } else r <- c(commutative = .all_true(lapply(e$args, .facts), "commutative"))
  # Expr._eval_is_extended_positive_negative: a number decides its sign by
  # evaluation (False for both when it is not real)
  if ((is.na(r["positive"] %||% NA) || is.na(r["negative"] %||% NA)) && .expr_is_number(e)) {
    z <- tryCatch(.mpeval_adaptive(e, 20L, list(), maxprec = 4000), error = function(err) NULL)
    if (!is.null(z) && isTRUE(attr(z, "converged"))) {
      if (z$im$m != 0) { r["extended_positive"] <- FALSE; r["extended_negative"] <- FALSE }
      else if (z$re$m != 0) {
        if (is.na(r["positive"] %||% NA)) r["extended_positive"] <- z$re$m > 0
        if (is.na(r["negative"] %||% NA)) r["extended_negative"] <- z$re$m < 0
      }
    }
  }
  r
}

#' Assumption shortcuts
#'
#' Each returns `TRUE`, `FALSE` or `NA` (unknown), like the SymPy attributes
#' `is_positive`, `is_real` and so on.
#'
#' @param e A `sym` expression.
#' @return `TRUE`, `FALSE` or `NA`.
#' @name assumption-queries
NULL

#' @rdname assumption-queries
#' @export
is_positive <- function(e) ask_fact(e, "positive")
#' @rdname assumption-queries
#' @export
is_negative <- function(e) ask_fact(e, "negative")
#' @rdname assumption-queries
#' @export
is_zero <- function(e) ask_fact(e, "zero")
#' @rdname assumption-queries
#' @export
is_nonzero <- function(e) ask_fact(e, "nonzero")
#' @rdname assumption-queries
#' @export
is_nonnegative <- function(e) ask_fact(e, "nonnegative")
#' @rdname assumption-queries
#' @export
is_nonpositive <- function(e) ask_fact(e, "nonpositive")
#' @rdname assumption-queries
#' @export
is_real <- function(e) ask_fact(e, "real")
#' @rdname assumption-queries
#' @export
is_complex <- function(e) ask_fact(e, "complex")
#' @rdname assumption-queries
#' @export
is_integer <- function(e) ask_fact(e, "integer")
#' @rdname assumption-queries
#' @export
is_rational <- function(e) ask_fact(e, "rational")
#' @rdname assumption-queries
#' @export
is_irrational <- function(e) ask_fact(e, "irrational")
#' @rdname assumption-queries
#' @export
is_even <- function(e) ask_fact(e, "even")
#' @rdname assumption-queries
#' @export
is_odd <- function(e) ask_fact(e, "odd")
#' @rdname assumption-queries
#' @export
is_prime_expr <- function(e) ask_fact(e, "prime")
#' @rdname assumption-queries
#' @export
is_finite <- function(e) ask_fact(e, "finite")
#' @rdname assumption-queries
#' @export
is_infinite <- function(e) ask_fact(e, "infinite")
#' @rdname assumption-queries
#' @export
is_imaginary <- function(e) ask_fact(e, "imaginary")
#' @rdname assumption-queries
#' @export
is_commutative <- function(e) ask_fact(e, "commutative")
#' @rdname assumption-queries
#' @export
is_algebraic <- function(e) ask_fact(e, "algebraic")
#' @rdname assumption-queries
#' @export
is_transcendental <- function(e) ask_fact(e, "transcendental")

#' Assumptions attached to a symbol
#'
#' @param x A symbol.
#' @return A named logical vector of every known fact.
#' @export
assumptions <- function(x) {
  x <- as_sym(x)
  f <- .facts(x)
  f[!is.na(f)]
}
