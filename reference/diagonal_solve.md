# Linear solvers of a matrix

The solvers of SymPy's `Matrix`: `diagonal_solve()`,
`lower_triangular_solve()`, `upper_triangular_solve()`,
`cholesky_solve()`, `LDLsolve()`, `LUsolve()`, `QRsolve()`,
`gauss_jordan_solve()` (returns `list(sol, params)`, plus the free
variable indices with `freevar = TRUE`), `pinv_solve()`,
`cramer_solve()`, `matrix_solve()` (SymPy's `Matrix.solve(rhs, method)`)
and `solve_least_squares()`.

## Usage

``` r
diagonal_solve(M, rhs)

lower_triangular_solve(M, rhs)

upper_triangular_solve(M, rhs)

cholesky_solve(M, rhs)

LDLsolve(M, rhs)

LUsolve(M, rhs, iszerofunc = NULL)

QRsolve(M, B)

gauss_jordan_solve(M, B, freevar = FALSE)

pinv_solve(M, B, arbitrary_matrix = NULL)

cramer_solve(M, rhs, det_method = "laplace")

matrix_solve(M, rhs, method = "GJ")

solve_least_squares(M, rhs, method = "CH")
```

## Arguments

- M:

  The coefficient `symmatrix`.

- rhs, B:

  The right hand side.

- iszerofunc:

  The zero test.

- freevar:

  Also return the indices of the free variables (1 based).

- arbitrary_matrix:

  For `pinv_solve()`: the arbitrary part; dummies by default.

- det_method:

  For `cramer_solve()`: `"laplace"`, `"bird"`, another
  [`det()`](https://rdrr.io/r/base/det.html) method or a function.

- method:

  The method: `"GJ"`, `"LU"`, `"CH"`, `"QR"`, `"LDL"`, `"PINV"`,
  `"CRAMER"`, or an
  [`inverse()`](https://adamrsophie.github.io/symcompR/reference/sym_det.md)
  method.

## Value

A `symmatrix` (a list for `gauss_jordan_solve()`).

## Examples

``` r
A <- Matrix(list(list(2, 1), list(1, 3)))
b <- Matrix(list(1, 2))
LUsolve(A, b)
#> Matrix(2, 1)
#> [1/5]
#> [3/5]
cholesky_solve(A, b)
#> Matrix(2, 1)
#> [1/5]
#> [3/5]
gauss_jordan_solve(Matrix(list(list(1, 2), list(2, 4))), Matrix(list(3, 6)))
#> [[1]]
#> Matrix(2, 1)
#> [3 - 2*tau0]
#> [tau0      ]
#> 
#> [[2]]
#> Matrix(1, 1)
#> [tau0]
#> 
```
