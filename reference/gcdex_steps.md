# Steps of the extended Euclidean algorithm

`gcdex_steps(f, g, ...)` lists the intermediate triples `(s, t, r)` with
`f*s + g*t = r` of the extended Euclidean algorithm, as SymPy's
`gcdex_steps()`; the last `r` is a greatest common divisor of `f` and
`g`. Multivariate polynomials need the generator, given in `...` or as
`gens`. Integers use the algorithm on the integers.

## Usage

``` r
gcdex_steps(f, g, ..., gens = NULL)
```

## Arguments

- f, g:

  Expressions.

- ...:

  The generator.

- gens:

  The generator.

## Value

A list of lists `list(s, t, r)`.

## Examples

``` r
x <- Symbol("x")
gcdex_steps(x^4 - 2*x^3 - 6*x^2 + 12*x + 15, x^3 + x^2 - 4*x - 4)
#> [[1]]
#> [[1]][[1]]
#> 0
#> 
#> [[1]][[2]]
#> 1
#> 
#> [[1]][[3]]
#> x^3 + x^2 - 4*x - 4
#> 
#> 
#> [[2]]
#> [[2]][[1]]
#> 1
#> 
#> [[2]][[2]]
#> 3 - x
#> 
#> [[2]][[3]]
#> x^2 + 4*x + 3
#> 
#> 
#> [[3]]
#> [[3]][[1]]
#> 3 - x
#> 
#> [[3]][[2]]
#> x^2 - 6*x + 10
#> 
#> [[3]][[3]]
#> 5*x + 5
#> 
#> 
```
