# Matrix logarithm and analytic functions of a matrix

`matrix_log(M)` is SymPy's `Matrix.log()`: the logarithm of each Jordan
block, transformed back (by default the result is `cancel`led).
`analytic_func(M, f, x)` evaluates the analytic function `f(x)` at the
matrix `M` by interpolating on the eigenvalues, as
`Matrix.analytic_func()`.

## Usage

``` r
matrix_log(M, simplify = cancel)

analytic_func(M, f, x)
```

## Arguments

- M:

  A square `symmatrix`.

- simplify:

  A function applied to `M` and to the result, or `FALSE`.

- f:

  An expression in `x`.

- x:

  The symbol of `f`.

## Value

A `symmatrix`.

## Examples

``` r
matrix_log(Matrix(list(list(1, 1), list(0, 1))))
#> Matrix(2, 2)
#> [0, 1]
#> [0, 0]
x <- Symbol("x")
analytic_func(Matrix(list(list(2, 1), list(0, 2))), exp(x), x)
#> Matrix(2, 2)
#> [exp(2), exp(2)]
#> [0     , exp(2)]
```
