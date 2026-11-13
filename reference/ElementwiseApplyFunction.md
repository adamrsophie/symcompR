# Elementwise application of a function to a matrix expression

`ElementwiseApplyFunction(fun, expr)` applies a scalar function to every
entry of the matrix `expr`, as SymPy's `expr.applyfunc(fun)` does for a
matrix expression. `fun` is an R function of one argument or a
[`Lambda()`](https://adamrsophie.github.io/symcompR/reference/Lambda.md);
it is kept as a `Lambda`. The result prints as
`Lambda(_d, f(_d)).(expr)`.

## Usage

``` r
ElementwiseApplyFunction(fun, expr)
```

## Arguments

- fun:

  A function of one argument or a `Lambda`.

- expr:

  A matrix expression or an explicit matrix.

## Value

A `sym` matrix expression.

## Examples

``` r
k <- Symbol("k")
X <- MatrixSymbol("X", k, k)
ElementwiseApplyFunction(sin, X)
#> Lambda(_d, sin(_d)).(X)
```
