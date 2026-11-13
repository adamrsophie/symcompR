# Solver helpers

`solve_linear_system(M, syms, ...)` solves an augmented matrix; further
symbols may follow `syms`. It returns `NULL` when there is no solution.
`solve_linear(lhs, rhs, symbols, exclude)` solves one equation for one
of `symbols` when it is linear in it, and returns the pair
`(symbol, value)`, `(0, 1)` when nothing is free, or
`(numerator, denominator)` otherwise. `checksol(f, symbol, sol)` tests a
solution by substitution and returns `TRUE`, `FALSE` or `NA` when
undecided. `solve_undetermined_coeffs(eq, coeffs, ...)` makes `eq`
vanish identically in the symbols `...`, as SymPy: a named list for one
solution, a list of them for several, and `NULL` when no solution
exists. `solve_poly_system(eqs, ...)` solves a zero dimensional
polynomial system in the generators `...` and returns a list of solution
tuples; `det_quick()` and `inv_quick()` are
[`sym_det()`](https://adamrsophie.github.io/symcompR/reference/sym_det.md)
and
[`inverse()`](https://adamrsophie.github.io/symcompR/reference/sym_det.md).

## Usage

``` r
solve_linear_system(M, syms, ...)

solve_linear(lhs, rhs = 0, symbols = list(), exclude = list(), x = NULL)

checksol(f, symbol, sol = NULL, ...)

solve_undetermined_coeffs(eq, coeffs, ..., dict = FALSE, set = FALSE)

solve_poly_system(eqs, ..., strict = FALSE)

det_quick(M)

inv_quick(M)
```

## Arguments

- M:

  An augmented matrix.

- syms, coeffs:

  Lists of symbols.

- ...:

  Further symbols for `solve_linear_system()`,
  `solve_undetermined_coeffs()`, or flags such as `numerical = FALSE`
  for `checksol()`.

- lhs, rhs:

  Sides of an equation.

- symbols, exclude:

  Symbols to solve for and symbols to skip.

- x:

  A variable; for `solve_linear()` an alias of `symbols`.

- f, eq:

  An equation or expression.

- symbol:

  A symbol, or a named list (or list of `list(symbol, value)` pairs)
  giving the whole solution.

- sol:

  The value of `symbol`.

- dict, set:

  Output forms of `solve_undetermined_coeffs()`, as for
  [`solve()`](https://rdrr.io/r/base/solve.html).

- eqs:

  A list of equations or expressions.

- strict:

  Kept for SymPy compatibility; roots that cannot be found are left out.

## Value

A list, expression, logical or set.

## Examples

``` r
x <- Symbol("x"); y <- Symbol("y")
checksol(x^2 - 4, list(x = 2))
#> [1] TRUE
```
