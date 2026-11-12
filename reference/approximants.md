# Pade approximants of a series

`approximants(l, X, simplify)` returns the successive Pade approximants
of the series whose coefficients are `l`, computed with the continued
fraction method of SymPy's `approximants()` generator.

## Usage

``` r
approximants(l, X = Symbol("x"), simplify = FALSE)
```

## Arguments

- l:

  A list of coefficients.

- X:

  The variable of the approximants.

- simplify:

  Simplify each approximant.

## Value

A list of rational functions.

## Examples

``` r
approximants(list(1, 1, 2, 5, 14, 42, 132))
#> [[1]]
#> 1
#> 
#> [[2]]
#> -1/(x - 1)
#> 
#> [[3]]
#> (x - 1)/(2*x - 1)
#> 
#> [[4]]
#> (1 - 2*x)/(x^2 - 3*x + 1)
#> 
#> [[5]]
#> (x^2 - 3*x + 1)/(3*x^2 - 4*x + 1)
#> 
#> [[6]]
#> (-3*x^2 + 4*x - 1)/(x^3 - 6*x^2 + 5*x - 1)
#> 
#> [[7]]
#> (x^3 - 6*x^2 + 5*x - 1)/(4*x^3 - 10*x^2 + 6*x - 1)
#> 
```
