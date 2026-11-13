# Kronecker products of matrices

`kronecker_product(...)` is the Kronecker (tensor) product of matrices
or matrix expressions, in canonical form: scalar factors come out in
front, nested products flatten and explicit matrices multiply out.
`combine_kronecker(expr)` merges sums, products and powers of Kronecker
products inside a matrix expression.

## Usage

``` r
kronecker_product(...)

combine_kronecker(expr)
```

## Arguments

- ...:

  Matrices or matrix expressions.

- expr:

  A matrix expression.

## Value

A matrix expression or an explicit matrix.

## Examples

``` r
n <- Symbol("n"); m <- Symbol("m")
A <- MatrixSymbol("A", n, m); B <- MatrixSymbol("B", n, m)
kronecker_product(2*A, B)
#> 2*(KroneckerProduct(A, B))
combine_kronecker(kronecker_product(A, B) + kronecker_product(B, B))
#> KroneckerProduct(A + B, B)
```
