# Handmade distributions

`ContinuousDistributionHandmade()`, `DiscreteDistributionHandmade()` and
`FiniteDistributionHandmade()` of `sympy.stats`: distributions given
directly by a density. The density of a continuous or discrete one is a
`Lambda`; a finite one takes a `Dict` (or a named list) from values to
probabilities. Read `"pdf"`, `"set"`, `"dict"` or `"pmf"` with
[`rv_attr()`](https://adamrsophie.github.io/symcompR/reference/rv_attr.md).

## Usage

``` r
ContinuousDistributionHandmade(pdf, set = NULL)

DiscreteDistributionHandmade(pdf, set = Integers)

FiniteDistributionHandmade(density)
```

## Arguments

- pdf:

  A `Lambda` giving the density.

- set:

  The support; by default the real line, or the integers for a discrete
  distribution.

- density:

  A `Dict` or a named list of probabilities.

## Value

A distribution.

## Examples

``` r
x <- Symbol("x")
CH <- ContinuousDistributionHandmade(Lambda(x, 2 * x), Interval(0, 1))
rv_attr(CH, "pdf", x)
#> 2*x
DH <- DiscreteDistributionHandmade(Lambda(x, 2^(-x)), Range(1, oo))
rv_attr(DH, "pdf", 3)
#> 1/8
FH <- FiniteDistributionHandmade(Dict(list(1, 2), list(Rational(1, 2), Rational(1, 2))))
rv_attr(FH, "set")
#> {1, 2}
```
