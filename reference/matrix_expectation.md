# Matrix expectation, variance and cross covariance

`ExpectationMatrix()`, `VarianceMatrix()` and `CrossCovarianceMatrix()`
of `sympy.stats`: the unevaluated expectation of a random matrix, and
the variance and cross covariance matrices of random vectors.
[`Expectation()`](https://adamrsophie.github.io/symcompR/reference/symbolic_probability.md),
[`Variance()`](https://adamrsophie.github.io/symcompR/reference/symbolic_probability.md)
and
[`Covariance()`](https://adamrsophie.github.io/symcompR/reference/symbolic_probability.md)
build them for matrix arguments;
[`expand()`](https://adamrsophie.github.io/symcompR/reference/expand.md)
applies linearity.

## Usage

``` r
ExpectationMatrix(expr, condition = NULL)

VarianceMatrix(arg, condition = NULL)

CrossCovarianceMatrix(arg1, arg2, condition = NULL)
```

## Arguments

- expr, arg, arg1, arg2:

  Matrix expressions (vectors for the variance and cross covariance).

- condition:

  An optional condition.

## Value

A matrix expression.

## Examples

``` r
k <- Symbol("k", integer = TRUE)
A <- MatrixSymbol("A", k, k)
X <- RandomMatrixSymbol("X", k, 1)
expand(ExpectationMatrix(A * X))
#> A*ExpectationMatrix(X)
expand(VarianceMatrix(A * X))
#> A*VarianceMatrix(X)*A.T
```
