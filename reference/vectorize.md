# Vectorize a function over nested lists

SymPy's `vectorize` decorator: `vectorize(0, 1)(f)` returns a function
that maps `f` over the elements of nested lists (or matrices) given as
its first and second arguments (positions count from 0; names select
named arguments), producing nested lists, the outer one following the
first vectorized argument. Without positions every argument is
vectorized.

## Usage

``` r
vectorize(...)
```

## Arguments

- ...:

  Argument positions (from 0) or names.

## Value

A function that takes a function and returns the vectorized function.

## Examples

``` r
x <- Symbol("x"); y <- Symbol("y")
vdiff <- vectorize(0, 1)(function(f, v) diff(f, v))
vdiff(list(x^2, x * y), list(x, y))
#> [[1]]
#> [[1]][[1]]
#> 2*x
#> 
#> [[1]][[2]]
#> 0
#> 
#> 
#> [[2]]
#> [[2]][[1]]
#> y
#> 
#> [[2]][[2]]
#> x
#> 
#> 
```
