# Numerator and denominator

`fraction(e)` returns `list(numer, denom)` by collecting factors with
negative exponents, like SymPy's `fraction()`; it does not combine sums.
`numer()` and `denom()` return one part.

## Usage

``` r
fraction(e, exact = FALSE)

numer(e)

denom(e)
```

## Arguments

- e:

  An expression.

- exact:

  If `TRUE`, only powers with a negative numeric exponent go to the
  denominator and the parts are not re-evaluated.

## Value

A list of two expressions.

## Examples

``` r
x <- Symbol("x"); y <- Symbol("y")
fraction(x / y)
#> [[1]]
#> x
#> 
#> [[2]]
#> y
#> 
```
