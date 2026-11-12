# Egyptian fractions

`egyptian_fraction(r)` writes a positive rational as a sum of distinct
unit fractions and returns their denominators. A rational above 1 first
takes the harmonic series `1 + 1/2 + 1/3 + ...` as far as it fits. The
`algorithm` for the rest is `"Greedy"`, `"Graham Jewett"`,
`"Takenouchi"` or `"Golomb"`.

## Usage

``` r
egyptian_fraction(r, algorithm = "Greedy")
```

## Arguments

- r:

  A positive rational, or `c(p, q)`.

- algorithm:

  Algorithm name.

## Value

A vector of denominators.

## Examples

``` r
egyptian_fraction(Rational(3, 7))
#> [1]   3  11 231
egyptian_fraction(c(3, 7), "Graham Jewett")
#> [1]    7    8    9   56   57   72 3192
egyptian_fraction(Rational(7, 3))
#> [1]  1  2  3  4  5 20
```
