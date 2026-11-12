# SymPy names that clash with functions of base R's default packages are
# exported with a sym_ prefix, so attaching symcompR masks nothing.

#' @rdname elementary-functions
#' @export
sym_atan2 <- atan2

#' @rdname special-functions
#' @export
sym_beta <- beta

#' Determinant, inverse, trace and friends
#'
#' `sym_det()` computes the determinant (Bareiss elimination by default, or
#' `"berkowitz"` or `"laplace"`). `inverse()` (also `solve(M)`) inverts by Gauss-Jordan
#' elimination; `adjugate()`, `cofactor()`, `minor()`, `trace_matrix()`,
#' `rank_matrix()`, `rref()`, `nullspace()`, `columnspace()` and
#' `rowspace()` follow SymPy. Given a characteristic curve such as
#' [TendonForceLengthDeGroote2016()], `inverse()` returns the constructor of
#' its inverse curve, as SymPy's `Function.inverse()`.
#'
#' @param x,M A `symmatrix`.
#' @param method Algorithm name.
#' @param ... Further arguments, such as `method`.
#' @param i,j Row and column.
#' @return An expression, matrix or list.
#' @examples
#' a <- Symbol("a")
#' M <- Matrix(list(list(a, 1), list(1, a)))
#' sym_det(M)
#' inverse(M)
#' rref(Matrix(list(list(1, 2, 3), list(4, 5, 6))))
#' @export
sym_det <- det

#' Factor polynomials
#'
#' `sym_factor()` factors a polynomial or rational function into irreducibles over
#' the rationals, in any number of variables, like SymPy's `sym_factor()`.
#' `factor_list()` returns `list(coeff, factors)` where `factors` is a list of
#' `list(factor, multiplicity)`. `sqf()`, `sqf_list()` and `sqf_part()` do
#' the same for square-free decomposition.
#'
#' @param x An expression.
#' @param f An expression.
#' @param ... Generators.
#' @param deep Factor inside function arguments.
#' @return An expression or a list.
#' @examples
#' x <- Symbol("x"); y <- Symbol("y")
#' sym_factor(x^4 - 1)
#' sym_factor(x^2 * y - y^3)
#' factor_list(2 * x^3 - 2 * x)
#' sqf((x - 1)^2 * (x + 2))
#' @export
sym_factor <- function(x, ...) factor_sym(x, ...)

#' @rdname special-functions
#' @export
sym_factorial <- factorial

#' @rdname elementary-functions
#' @export
sym_Mod <- Mod

#' @rdname core-helpers
#' @export
sym_ordered <- ordered

#' @rdname special-functions
#' @export
sym_binomial <- binomial

#' Integrate an expression
#'
#' `sym_integrate(f, x)` returns an antiderivative; `sym_integrate(f, list(x, a, b))`
#' a definite integral, evaluated with limits at the endpoints. Several
#' variables integrate in turn. Rules cover polynomials, rational functions
#' (partial fractions with logarithms and arctangents), exponentials,
#' logarithms, trigonometric and hyperbolic functions and their powers,
#' inverse trigonometric functions, substitution, integration by parts and
#' a table of special forms (`erf`, `Si`, `Ei`, `li`, ...). When no rule
#' applies an unevaluated `Integral` is returned. For an R function as first
#' argument this calls [stats::integrate()].
#'
#' @param f An expression, or an R function for [stats::integrate()].
#' @param ... Variables or `list(var, lower, upper)` specifications, or
#'   arguments for [stats::integrate()].
#' @param meijerg,risch,heurisch,manual `TRUE` to use only that algorithm,
#'   `FALSE` to skip it, `NULL` (the default) to try it among the others,
#'   as SymPy's hints.
#' @param conds How conditions of convergence are returned: `"piecewise"`,
#'   `"separate"` or `"none"`.
#' @return An expression.
#' @examples
#' x <- Symbol("x")
#' sym_integrate(x^2 * exp(x), x)
#' sym_integrate(1 / (x^2 + 1), list(x, 0, 1))
#' sym_integrate(exp(-x^2), list(x, -oo, oo))
#' @export
sym_integrate <- integrate

#' @rdname poly-aliases
#' @export
sym_poly <- poly

#' @rdname function-aliases
#' @export
sym_rf <- rf

#' @rdname Symbol
#' @export
sym_var <- var

#' @rdname gcd
#' @export
sym_lcm <- lcm

#' @rdname Symbol
#' @export
sym_symbols <- symbols

# Dummy._count after `import sympy` (SymPy 1.15): the first Dummy() of a
# session is _Dummy_34, as in SymPy. The numbering drifts wherever symcompR
# creates internal dummies that SymPy does not, or the reverse.
.dummy_counter$n <- 34L
