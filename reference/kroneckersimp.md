# Simplify products of Kronecker deltas

`kroneckersimp(expr)` sets to zero the products of `KroneckerDelta`s
that cannot both be one, such as
`KroneckerDelta(0, j)*KroneckerDelta(1, j)`.

## Usage

``` r
kroneckersimp(expr)
```

## Arguments

- expr:

  An expression.

## Value

An expression.

## Examples

``` r
j <- Symbol("j")
kroneckersimp(1 + KroneckerDelta(0, j)*KroneckerDelta(1, j))
#> 1
```
