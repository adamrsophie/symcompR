# Harmonic prefix of an Egyptian fraction

`egypt_harmonic(r)` takes the terms `1, 1/2, 1/3, ...` of the harmonic
series while their sum stays at most `r`, and returns the list of the
denominators taken and the remainder, as SymPy's `egypt_harmonic()` (the
first step of
[`egyptian_fraction()`](https://adamrsophie.github.io/symcompR/reference/egyptian_fraction.md)
for rationals above 1).

## Usage

``` r
egypt_harmonic(r)
```

## Arguments

- r:

  A rational number.

## Value

A list of a list of integers and a rational.

## Examples

``` r
egypt_harmonic(Rational(5, 2))
#> [[1]]
#> [[1]][[1]]
#> 1
#> 
#> [[1]][[2]]
#> 2
#> 
#> [[1]][[3]]
#> 3
#> 
#> [[1]][[4]]
#> 4
#> 
#> [[1]][[5]]
#> 5
#> 
#> [[1]][[6]]
#> 6
#> 
#> 
#> [[2]]
#> 1/20
#> 
egypt_harmonic(Rational(1, 3))
#> [[1]]
#> list()
#> 
#> [[2]]
#> 1/3
#> 
```
