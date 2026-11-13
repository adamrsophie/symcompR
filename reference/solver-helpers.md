# Solver helpers

`linear_eq_to_matrix(eqs, syms)` returns `list(A, b)` with `A x = b`.
`solve_linear_system(M, syms, ...)` solves an augmented matrix; further
symbols may follow `syms`. It returns `NULL` when there is no solution.
`solve_linear(lhs, rhs, symbols, exclude)` solves one equation for one
of `symbols` when it is linear in it, and returns the pair
`(symbol, value)`, `(0, 1)` when nothing is free, or
`(numerator, denominator)` otherwise. `checksol(f, symbol, sol)` tests a
solution by substitution and returns `TRUE`, `FALSE` or `NA` when
undecided; `checkodesol(eq, sol, fx)` checks an ODE solution and returns
the pair `(TRUE, 0)` or `(FALSE, residual)`.
`solve_undetermined_coeffs(eq, coeffs, ...)` makes `eq` vanish
identically in the symbols `...`, as SymPy: a named list for one
solution, a list of them for several, and `NULL` when no solution
exists. `nonlinsolve()` solves a system; `solve_poly_system(eqs, ...)`
solves a zero dimensional polynomial system in the generators `...` and
returns a list of solution tuples; `det_quick()` and `inv_quick()` are
[`sym_det()`](https://adamrsophie.github.io/symcompR/reference/sym_det.md)
and
[`inverse()`](https://adamrsophie.github.io/symcompR/reference/sym_det.md).

## Usage

``` r
linear_eq_to_matrix(eqs, syms)

solve_linear_system(M, syms, ...)

solve_linear(lhs, rhs = 0, symbols = list(), exclude = list(), x = NULL)

checksol(f, symbol, sol = NULL, ...)

checkodesol(eq, sol, fx = NULL)

solve_undetermined_coeffs(eq, coeffs, ..., dict = FALSE, set = FALSE)

nonlinsolve(eqs, syms, ...)

solve_poly_system(eqs, ..., strict = FALSE)

det_quick(M)

inv_quick(M)
```

## Arguments

- eqs:

  A list of equations or expressions.

- syms, coeffs:

  Lists of symbols.

- M:

  An augmented matrix.

- ...:

  Further symbols for `solve_linear_system()`,
  `solve_undetermined_coeffs()` or `nonlinsolve()`, or flags such as
  `numerical = FALSE` for `checksol()`.

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

  The value of `symbol`, or an
  [`Eq()`](https://adamrsophie.github.io/symcompR/reference/Eq.md) for
  `checkodesol()`.

- fx:

  Unknown function application; found from `eq` when omitted.

- dict, set:

  Output forms of `solve_undetermined_coeffs()`, as for
  [`solve()`](https://rdrr.io/r/base/solve.html).

- strict:

  Kept for SymPy compatibility; roots that cannot be found are left out.

## Value

A list, expression, logical or set.

## Examples

``` r
x <- Symbol("x"); y <- Symbol("y")
linear_eq_to_matrix(list(x + 2*y - 3, 3*x - y), list(x, y))
#> $A
#> Matrix(2, 2)
#> [1, 2 ]
#> [3, -1]
#> 
#> $b
#> Matrix(2, 1)
#> [3]
#> [0]
#> 
checksol(x^2 - 4, list(x = 2))
#> [1] TRUE
```
