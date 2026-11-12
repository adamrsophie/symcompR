# Solver helpers

`linear_eq_to_matrix(eqs, syms)` returns `list(A, b)` with `A x = b`.
`solve_linear_system(M, syms)` solves an augmented matrix.
`solve_linear(lhs, rhs, x)` solves one linear equation for `x`.
`checksol(f, sol)` tests a solution by substitution, `checkodesol()` an
ODE solution. `solve_undetermined_coeffs(eq, coeffs, x)` makes `eq`
vanish identically in `x`. `reduce_inequalities()` and
`solve_poly_inequality()` solve inequalities in one variable.
`nonlinsolve()` and `solve_poly_system()` solve systems; `det_quick()`
and `inv_quick()` are
[`sym_det()`](https://adamrsophie.github.io/symcompR/reference/sym_det.md)
and
[`inverse()`](https://adamrsophie.github.io/symcompR/reference/sym_det.md).

## Usage

``` r
linear_eq_to_matrix(eqs, syms)

solve_linear_system(M, syms)

solve_linear(lhs, rhs = 0, x = NULL)

checksol(f, sol)

checkodesol(eq, sol, fx)

solve_undetermined_coeffs(eq, coeffs, x)

reduce_inequalities(ineqs, x = NULL)

solve_poly_inequality(poly, rel = ">")

nonlinsolve(eqs, syms)

solve_poly_system(eqs, syms)

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

- lhs, rhs:

  Sides of an equation.

- x:

  A variable.

- f, eq:

  An equation or expression.

- sol:

  A named list of values, or an
  [`Eq()`](https://adamrsophie.github.io/symcompR/reference/Eq.md) for
  `checkodesol()`.

- fx:

  Unknown function application.

- ineqs:

  An inequality or list of inequalities.

- poly:

  A polynomial.

- rel:

  Relational operator such as `">"`.

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
