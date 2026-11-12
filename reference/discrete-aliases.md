# Discrete transform aliases

`ifft()` is
[`sym_ifft()`](https://adamrsophie.github.io/symcompR/reference/convolution.md)
for lists of symbolic or exact values. `fwht()` and `ifwht()` are the
Walsh-Hadamard transform and its inverse. `covering_product()` and
`intersecting_product()` multiply sequences indexed by bitmasks,
combining over unions and intersections of the index sets.

## Usage

``` r
ifft(seq)

fwht(seq)

ifwht(seq)

covering_product(a, b)

intersecting_product(a, b)
```

## Arguments

- seq, a, b:

  Sequences.

## Value

A list of expressions.

## Examples

``` r
covering_product(list(1, 2), list(3, 4))
#> [[1]]
#> 3
#> 
#> [[2]]
#> 18
#> 
```
