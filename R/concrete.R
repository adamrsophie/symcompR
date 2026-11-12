# Summation and products, after sympy/concrete (summations.py, products.py,
# gosper.py).

#' Symbolic sums and products
#'
#' `summation(f, list(i, a, b))` computes `sum(f, i = a..b)` in closed form
#' with SymPy's `Sum.doit()`: direct sums of short ranges, telescoping,
#' Faulhaber polynomials, geometric sums (a Piecewise for ratio 1), Gosper's
#' algorithm, hypergeometric series through [hyperexpand()] (a Piecewise
#' with the convergence condition), sums of rational functions over infinite
#' ranges by residues, `zeta` sums and Kronecker deltas
#' (`deltasummation()`). `gosper_sum(f, k)` is the indefinite Gosper sum.
#' `product(f, list(i, a, b))` is SymPy's `Product.doit()`: direct products of
#' short ranges, rising factorials over the roots of polynomial factors, sums
#' in exponents, and Kronecker deltas (`deltaproduct()`). `Sum()` and `Product()`
#' build unevaluated objects; use [doit()]. `is_convergent()` applies
#' standard convergence tests.
#'
#' @param f The summand or factor.
#' @param ... Limits `list(i, a, b)`.
#' @return An expression.
#' @examples
#' i <- Symbol("i", integer = TRUE); n <- Symbol("n", integer = TRUE)
#' summation(i^2, list(i, 1, n))
#' summation(1 / i^2, list(i, 1, oo))
#' product(i, list(i, 1, n))
#' @export
summation <- function(f, ...) .sum_doit(Sum(f, ...), deep = FALSE)

#' @rdname summation
#' @export
Sum <- function(f, ...) {
  f <- as_sym(f)
  lims <- lapply(list(...), function(l) lapply(l, as_sym))
  .sum_new(f, lims)
}

# AddWithLimits.__new__ through _common_new: an equation sums both sides,
# nested sums are denested (inner limits first), and Piecewise parts that
# contain a summation variable are folded to the top
.sum_new <- function(f, lims, type = "Sum") {
  if (is_nan(f)) return(S_NaN)
  if (f$type == "Equality") return(Eq(.sum_new(f$args[[1L]], lims, type), .sum_new(f$args[[2L]], lims, type)))
  while (f$type == type) {
    lims <- c(lapply(f$args[-1L], function(t) t$args), lims)
    f <- f$args[[1L]]
  }
  if (.has_type(f, "Piecewise")) {
    vars <- lapply(lims, `[[`, 1L)
    reps <- list()
    for (p in atoms(f, "Piecewise")) {
      if (!any(vapply(vars, function(v) .has1(p, v), TRUE))) reps[[length(reps) + 1L]] <- list(p, Dummy())
    }
    f <- xreplace(f, reps)
    f <- piecewise_fold(f)
    f <- xreplace(f, lapply(reps, function(r) list(r[[2L]], r[[1L]])))
  }
  new_sym(type, c(list(f), lapply(lims, function(l) do.call(Tuple, l))))
}

#' @rdname summation
#' @export
product <- function(f, ...) {
  p <- Product(f, ...)
  if (p$type == "Product") .prod_doit(p, deep = FALSE) else p
}

#' @rdname summation
#' @export
Product <- function(f, ...) .sum_new(as_sym(f), lapply(list(...), function(l) lapply(l, as_sym)), "Product")

for (.t in c("Sum", "Product")) {
  register_bound_vars(.t, function(e) vapply(e$args[-1L], function(t) t$args[[1L]]$key, ""))
  register_str_printer(.t, local({
    tt <- .t
    function(e) paste0(tt, "(", .str(e$args[[1L]]), ", ", paste(vapply(e$args[-1L], .str, ""), collapse = ", "), ")")
  }))
}
rm(.t)
register_doit("Sum", function(e, deep = TRUE) .sum_doit(e, deep))
register_doit("Product", function(e, deep = TRUE) .prod_doit(e, deep))

.poly_deg <- function(p, i) {
  p <- expand(p)
  if (is_zero_num(p)) return(-1L)
  length(.coeff_list(p, i)) - 1L
}

#' @rdname summation
#' @export
gosper_sum <- function(f, ...) {
  f <- as_sym(f)
  k <- list(...)[[1L]]
  # gosper_sum(f, k) is the indefinite sum, gosper_sum(f, list(k, a, b)) the definite one
  if (is_sym(k)) return(.gosper_sum_f(f, k))
  lim <- lapply(k, as_sym)
  .gosper_sum_f(f, lim[[1L]], lim[[2L]], lim[[3L]])
}

