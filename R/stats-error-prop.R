# Error propagation (sympy/stats/error_prop.py) and the public RandomSymbol
# constructor.

#' Random symbols
#'
#' `RandomSymbol(symbol, pspace)` wraps a symbol as a random variable of the
#' probability space `pspace`; without a space it belongs to the abstract
#' space `PSpace()`, as in SymPy.
#'
#' @param symbol A symbol or a name.
#' @param pspace A probability space, or `NULL`.
#' @return A random symbol.
#' @examples
#' a <- RandomSymbol(Symbol("x"))
#' Variance(a)
#' @export
RandomSymbol <- function(symbol, pspace = NULL) .RandomSymbol(symbol, pspace)

.arg0_or_var <- function(v) if (length(v$args) > 0L) v$args[[1L]] else v

#' Propagation of variance
#'
#' Symbolically propagates the variance of the symbols in `expr` (SymPy's
#' `variance_prop()`); every symbol is treated as a random variable unless
#' listed in `consts`.
#'
#' @param expr An expression.
#' @param consts Symbols that are constants.
#' @param include_covar Add covariance terms.
#' @return An expression in `Variance()` and `Covariance()` terms.
#' @examples
#' x <- Symbol("x"); y <- Symbol("y")
#' variance_prop(x * y)
#' @export
variance_prop <- function(expr, consts = list(), include_covar = FALSE) {
  expr <- as_sym(expr)
  if (is_sym(consts)) consts <- consts$args
  args <- expr$args
  if (!length(args) || .is_rs(expr)) {
    if (any(vapply(consts, function(z) identical(as_sym(z)$key, expr$key), NA))) return(S_Zero)
    if (is_random(expr)) return(doit(Variance(expr)))
    if (expr$type == "Symbol") return(doit(Variance(RandomSymbol(expr))))
    return(S_Zero)
  }
  var_args <- lapply(args, variance_prop, consts = consts, include_covar = include_covar)
  pairs <- function(n) if (n < 2L) list() else unlist(lapply(1:(n - 1L), function(i) lapply((i + 1L):n, function(j) c(i, j))), recursive = FALSE)
  if (expr$type == "Add") {
    var_expr <- .add(var_args)
    if (include_covar) {
      terms <- lapply(pairs(length(var_args)), function(ij)
        .mul(list(as_sym(2), expand(Covariance(.arg0_or_var(var_args[[ij[1]]]), .arg0_or_var(var_args[[ij[2]]]))))))
      var_expr <- .add(list(var_expr, .add(terms)))
    }
  } else if (expr$type == "Mul") {
    terms <- lapply(seq_along(args), function(i) var_args[[i]] / args[[i]]^2)
    var_expr <- simplify(expr^2 * .add(terms))
    if (include_covar) {
      terms <- lapply(pairs(length(args)), function(ij)
        .mul(list(as_sym(2), expand(Covariance(.arg0_or_var(var_args[[ij[1]]]), .arg0_or_var(var_args[[ij[2]]]))))) /
          (args[[ij[1]]] * args[[ij[2]]]))
      var_expr <- .add(list(var_expr, .add(terms)))
    }
  } else if (expr$type == "Pow") {
    b <- args[[2L]]
    var_expr <- simplify(var_args[[1L]] * (expr * b / args[[1L]])^2)
  } else if (expr$type == "Function" && identical(expr$name, "exp")) {
    var_expr <- simplify(var_args[[1L]] * expr^2)
  } else {
    var_expr <- Variance(expr)
  }
  var_expr
}
