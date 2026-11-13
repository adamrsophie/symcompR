# Ranges, least common multiples of periods and monotonicity, after
# function_range() and lcim() in sympy/calculus/util.py and
# monotonicity_helper() in sympy/calculus/singularities.py.

# the members of a set as a list (Python iteration over a set)
.cu_iter <- function(s) {
  if (s$type == "EmptySet") return(list())
  if (s$type == "FiniteSet") return(s$args)
  if (s$type == "Union") return(s$args)
  NULL
}

#' Range of a function
#'
#' `function_range(f, symbol, domain)` returns the range of the real
#' function `f` of `symbol` over the real set `domain`, as SymPy's
#' `function_range()`: the function is followed on each interval of its
#' continuous domain through its limits at open ends, its values at closed
#' ends and its critical points.
#'
#' @param f An expression.
#' @param symbol The variable.
#' @param domain A subset of the reals.
#' @return A set.
#' @examples
#' x <- Symbol("x")
#' function_range(sin(x), x, Interval(0, 2 * S_Pi))
#' function_range(1 / x, x, Reals)
#' @export
function_range <- function(f, symbol, domain) {
  f <- as_sym(f); symbol <- as_sym(symbol); domain <- as_sym(domain)
  if (domain$type == "EmptySet") return(S_EmptySet)
  if (f$type == "Function" && length(f$args) == 1L && !identical(f$args[[1L]]$key, symbol$key)) {
    dec <- decompogen(f, symbol)
    if (length(dec) > 1L) {
      cur <- function_range(dec[[length(dec)]], symbol, domain)
      for (fn in rev(dec[-length(dec)])) cur <- function_range(fn, symbol, cur)
      return(cur)
    }
  }
  # periodicity() is zero when f simplifies to a constant in symbol
  if (!(symbol$key %in% vapply(free_symbols(simplify(f)), .k, ""))) return(FiniteSet(expand(f)))
  period <- periodicity(f, symbol)
  if (!is.null(period) && is_zero_num(period)) return(FiniteSet(expand(f)))
  if (!is.null(period)) {
    infinite_iv <- function(s) s$type == "Interval" || .set_cls(s) == "Reals"
    if (infinite_iv(domain)) {
      iv <- .iv(domain)
      if (isTRUE(ask_fact(.sub(iv$a, iv$b), "infinite"))) domain <- Interval(0, period)
    } else if (domain$type == "Union") {
      for (sd in domain$args) if (infinite_iv(sd)) {
        iv <- .iv(sd)
        if (isTRUE(ask_fact(.sub(iv$a, iv$b), "infinite"))) domain <- Interval(0, period)
      }
    }
  }
  intervals <- continuous_domain(f, symbol, domain)
  range_int <- S_EmptySet
  icls <- .set_cls(intervals)
  interval_iter <- if (icls %in% c("Interval", "FiniteSet", "Reals")) list(intervals)
                   else if (icls == "Union") intervals$args
                   else stop(.sympy_condition("NotImplementedError", "Unable to find range for the given domain."))
  for (interval in interval_iter) {
    cls <- .set_cls(interval)
    if (cls == "FiniteSet") {
      for (singleton in interval$args)
        if (isTRUE(.py_in(domain, singleton))) range_int <- Union(range_int, FiniteSet(subs(f, symbol, singleton)))
    } else if (cls %in% c("Interval", "Reals")) {
      iv <- .iv(interval)
      vals <- S_EmptySet
      critical_values <- S_EmptySet
      bounds <- list(list(iv$lo, iv$a, "+"), list(iv$ro, iv$b, "-"))
      for (bd in bounds) {
        if (bd[[1L]]) {
          critical_values <- Union(critical_values, FiniteSet(limit(f, symbol, bd[[2L]], bd[[3L]])))
          vals <- Union(vals, critical_values)
        } else vals <- Union(vals, FiniteSet(subs(f, symbol, bd[[2L]])))
      }
      if (identical(set_inf(vals)$key, S_NegativeInfinity$key) && identical(set_sup(vals)$key, S_Infinity$key)) {
        range_int <- Union(range_int, Reals)
        next
      }
      critical_points <- solveset(.diff1(f, symbol), symbol, interval)
      if (critical_points$type == "ImageSet")
        stop(.sympy_condition("NotImplementedError", sprintf("Infinite number of critical points for %s", sstr(f))))
      cp <- .cu_iter(critical_points)
      if (is.null(cp) || critical_points$type == "Union")
        stop(.sympy_condition("NotImplementedError", sprintf("Unable to find critical points for %s", sstr(f))))
      for (pt in cp) vals <- Union(vals, FiniteSet(subs(f, symbol, pt)))
      left_open <- FALSE; right_open <- FALSE
      if (critical_values$type != "EmptySet") {
        if (.py_eq(set_inf(critical_values), set_inf(vals))) left_open <- TRUE
        if (.py_eq(set_sup(critical_values), set_sup(vals))) right_open <- TRUE
      }
      range_int <- Union(range_int, Interval(set_inf(vals), set_sup(vals), left_open, right_open))
    } else stop(.sympy_condition("NotImplementedError", "Unable to find range for the given domain."))
  }
  excluded <- Complement(domain, intervals)
  if (excluded$type == "FiniteSet") {
    for (pt in excluded$args) {
      val <- tryCatch(subs(f, symbol, pt), error = function(e) NULL)
      if (!is.null(val) && isTRUE(ask_fact(val, "real"))) range_int <- Union(range_int, FiniteSet(val))
    }
  }
  range_int
}

#' Least common multiple of periods
#'
#' `lcim(numbers)` is the least common multiple of a list of rational
#' numbers, or of irrational numbers with the same irrational factor (such
#' as multiples of `pi`), and `NULL` otherwise, as SymPy's `lcim()`.
#'
#' @param numbers A list of numbers.
#' @return A number or `NULL`.
#' @examples
#' lcim(list(Rational(1, 2), 2, 3))
#' lcim(list(S_Pi, S_Pi / 2, S_Pi / 4))
#' @export
lcim <- function(numbers) {
  numbers <- lapply(as.list(numbers), as_sym)
  result <- NULL
  if (all(vapply(numbers, function(n) isTRUE(ask_fact(n, "irrational")), TRUE))) {
    factors_num <- lapply(numbers, function(n) .as_coeff_mul1(factor_sym(n)))
    term <- factors_num[[1L]][[2L]]
    if (all(vapply(factors_num, function(cf) identical(cf[[2L]]$key, term$key), TRUE)))
      result <- .mul(list(lcm_list(lapply(factors_num, `[[`, 1L)), term))
  } else if (all(vapply(numbers, function(n) isTRUE(ask_fact(n, "rational")), TRUE))) {
    result <- lcm_list(numbers)
  }
  result
}

# Expr.as_coeff_Mul(): the leading rational coefficient and the rest
.as_coeff_mul1 <- function(e) {
  if (e$type %in% c("Integer", "Rational", "Float")) return(list(e, S_One))
  if (e$type == "Mul" && e$args[[1L]]$type %in% c("Integer", "Rational", "Float")) {
    rest <- e$args[-1L]
    return(list(e$args[[1L]], if (length(rest) == 1L) rest[[1L]] else new_sym("Mul", rest)))
  }
  list(S_One, e)
}

#' Monotonicity on an interval
#'
#' `monotonicity_helper(expression, predicate, interval, symbol)` tells
#' whether `predicate(derivative)` holds on all of `interval`, after
#' checking that `expression` has no singularity inside it, as SymPy's
#' helper behind `is_increasing()` and its relatives. `predicate` maps the
#' derivative to a relational, such as `function(d) Ge(d, 0)`.
#'
#' @param expression An expression.
#' @param predicate A function of the derivative returning a relational.
#' @param interval A set.
#' @param symbol The variable, or `NULL` for the only free symbol.
#' @return `TRUE`, `FALSE` or `NA` when undecided.
#' @examples
#' x <- Symbol("x")
#' monotonicity_helper(x^3, function(d) Ge(d, 0))
#' @export
monotonicity_helper <- function(expression, predicate, interval = Reals, symbol = NULL) {
  expression <- as_sym(expression); interval <- as_sym(interval)
  free <- free_symbols(expression)
  if (is.null(symbol) && length(free) > 1L)
    stop(.sympy_condition("NotImplementedError", "The function has not yet been implemented for all multivariate expressions."))
  variable <- if (!is.null(symbol)) as_sym(symbol) else if (length(free)) free[[1L]] else Symbol("x")
  sings <- tryCatch(singularities(expression, variable, interval), error = function(e) {
    if (inherits(e, "NotImplementedError")) NULL else stop(e)
  })
  if (!is.null(sings) && isTRUE(is_subset(interval, Reals))) {
    interior_sings <- Intersection(.cu_interior(interval), sings)
    if (interior_sings$type != "EmptySet") return(FALSE)
  }
  derivative <- .diff1(expression, variable)
  predicate_interval <- solveset(predicate(derivative), variable, Reals)
  is_subset(interval, predicate_interval)
}

# Set.interior for intervals, unions of intervals and finite sets
.cu_interior <- function(s) {
  cls <- .set_cls(s)
  if (cls == "Interval") { iv <- .iv(s); return(Interval(iv$a, iv$b, TRUE, TRUE)) }
  if (cls == "Union") return(do.call(Union, lapply(s$args, .cu_interior)))
  if (cls == "FiniteSet") return(S_EmptySet)
  Complement(s, boundary(s))
}
