# Collect sums of tensor products

Rewrites every
[`ArrayAdd()`](https://adamrsophie.github.io/symcompR/reference/array-expressions.md)
in `expr` whose terms are tensor products that differ in at most one
factor, so that `c1*(A x B) + c2*(C x B)` becomes `(c1*A + c2*C) x B`.

## Usage

``` r
collect_tensor_products(expr)
```

## Arguments

- expr:

  An array expression.

## Value

The rewritten expression.

## Examples

``` r
k <- Symbol("k")
A <- MatrixSymbol("A", k, k); B <- MatrixSymbol("B", k, k); C <- MatrixSymbol("C", k, k)
collect_tensor_products(ArrayAdd(ArrayTensorProduct(A, B), ArrayTensorProduct(A, C)))
#> ArrayTensorProduct(A, B + C)
```
