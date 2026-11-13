# Matrix distributions

Random matrices of `sympy.stats`: `MatrixGamma()`, `Wishart()`,
`MatrixNormal()` and `MatrixStudentT()`. Matrix arguments are explicit
matrices, nested lists (rows) or matrix symbols. `rv_density(M)` is the
distribution, which is called on a matrix with
[`apply_lambda()`](https://adamrsophie.github.io/symcompR/reference/Lambda.md).

## Usage

``` r
MatrixGamma(symbol, alpha, beta, scale_matrix)

Wishart(symbol, n, scale_matrix)

MatrixNormal(symbol, location_matrix, scale_matrix_1, scale_matrix_2)

MatrixStudentT(symbol, nu, location_matrix, scale_matrix_1, scale_matrix_2)
```

## Arguments

- symbol:

  The name of the random matrix.

- alpha, beta:

  Shape and scale parameters.

- scale_matrix, scale_matrix_1, scale_matrix_2:

  Scale matrices.

- n:

  Degrees of freedom of the Wishart distribution.

- location_matrix:

  The location (mean) matrix.

- nu:

  Degrees of freedom of the matrix Student t distribution.

## Value

A random matrix symbol.

## Examples

``` r
M <- MatrixGamma("M", 1, 2, list(list(1, 0), list(0, 1)))
rv_attr(rv_attr(rv_attr(M, "pspace"), "distribution"), "set")
#> MatrixSet(2, 2, Reals)
```
