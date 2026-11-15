# ODE solver hints

`allhints` lists every hint
[`classify_ode()`](https://adamrsophie.github.io/symcompR/reference/classify_ode.md)
can return, in the order it tries them, as SymPy's `ode.allhints`.

## Usage

``` r
allhints
```

## Format

A character vector.

## Examples

``` r
head(allhints)
#> [1] "factorable"    "nth_algebraic" "separable"     "1st_exact"    
#> [5] "1st_linear"    "Bernoulli"    
```
