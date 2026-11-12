# Multiplicities of factors

`multiplicity(p, n)` is the largest `m` with `p^m` dividing `n`;
rational `p` and `n` (as `sym` rationals or `gmp` big rationals) give
the exponent of `p` in `n`, possibly negative, and an unevaluated
`sym_factorial(k)` for `n` is handled by
`multiplicity_in_factorial(p, k)`, the exponent of `p` in `k!`.

## Usage

``` r
multiplicity(p, n)

multiplicity_in_factorial(p, n)
```

## Arguments

- p:

  Base: an integer greater than 1, or a rational. For an
  [`ImplicitRegion()`](https://adamrsophie.github.io/symcompR/reference/ImplicitRegion.md),
  the multiplicity of the point `n` on the region.

- n:

  An integer or rational, or a point of an implicit region.

## Value

An integer.

## Examples

``` r
multiplicity(5, 8000)
#> [1] 3
multiplicity(S(3), Rational(1, 9))
#> [1] -2
multiplicity_in_factorial(2, 100)
#> [1] 97
```
