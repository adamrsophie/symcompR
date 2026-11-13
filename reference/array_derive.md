# Derivatives of array and matrix expressions

`array_derive(expr, x)` differentiates an array expression by an array
or matrix symbol, giving an array expression whose first axes are those
of `x`. `matrix_derive(expr, x)` differentiates a matrix expression and
converts the result back to a matrix expression where possible. These
follow `sympy.tensor.array.expressions.arrayexpr_derivatives`;
[`diff()`](https://rdrr.io/r/base/diff.html) uses them for matrix
expressions.

## Usage

``` r
array_derive(expr, x)

matrix_derive(expr, x)

ArrayDerivative(expr, ..., evaluate = FALSE)
```

## Arguments

- expr:

  An array or matrix expression.

- x:

  A matrix or array symbol.

- ...:

  Variables of differentiation, each optionally followed by a count.

- evaluate:

  If `TRUE`, compute the derivative; by default `ArrayDerivative()`
  stays unevaluated, as in SymPy.

## Value

An array or matrix expression.

## Examples

``` r
k <- Symbol("k")
X <- MatrixSymbol("X", k, k); A <- MatrixSymbol("A", k, k)
array_derive(convert_matrix_to_array(A * X), X)
#> PermuteDims(ArrayContraction(ArrayTensorProduct(I, I, A), (1, 5)), (2 3))
matrix_derive(Trace(A * X), X)
#> A.T
```
