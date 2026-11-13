# Hadamard (elementwise) products and powers

`hadamard_product(...)` multiplies matrices elementwise and puts the
result in canonical form; `hadamard_power(base, exp)` raises a matrix
elementwise to a power. With explicit matrices only, the result is an
explicit matrix.

## Usage

``` r
hadamard_product(...)

hadamard_power(base, exp)
```

## Arguments

- ...:

  Matrices or matrix expressions.

- base:

  A matrix expression or a scalar.

- exp:

  The exponent, a scalar.

## Value

A matrix expression or an explicit matrix.

## Examples

``` r
k <- Symbol("k")
A <- MatrixSymbol("A", k, k); B <- MatrixSymbol("B", k, k)
hadamard_product(B, A, A)
#> HadamardProduct(HadamardPower(A, 2), B)
hadamard_power(A, 2)
#> HadamardPower(A, 2)
```
