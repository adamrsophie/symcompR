# Nest the permutation of a PermuteDims

Moves the axis permutation of a `PermuteDims` expression inside its
argument where possible, as SymPy's `nest_permutation()`.

## Usage

``` r
nest_permutation(expr)
```

## Arguments

- expr:

  An array expression.

## Value

An array expression.

## Examples

``` r
M <- ArraySymbol("M", c(3, 3)); N <- ArraySymbol("N", c(3, 3))
nest_permutation(PermuteDims(ArrayTensorProduct(M, N), c(1, 0, 3, 2)))
#> PermuteDims(ArrayTensorProduct(M, N), (0 1)(2 3))
```
