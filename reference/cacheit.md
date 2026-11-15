# Cache a function

SymPy's `cacheit` decorator: `cacheit(f)` returns a function that
remembers its results by the arguments, in the cache that
[`clear_cache()`](https://adamrsophie.github.io/symcompR/reference/clear_cache.md)
empties.

## Usage

``` r
cacheit(f)
```

## Arguments

- f:

  A function.

## Value

The caching function.

## Examples

``` r
slow_square <- function(x) { Sys.sleep(0.1); x^2 }
fast <- cacheit(slow_square)
fast(Symbol("x")); fast(Symbol("x"))
#> x^2
#> x^2
```
