# Raise an error on 0/0

SymPy's `seterr(divide)`: with `divide = TRUE`, `Rational(0, 0)` raises
an error instead of returning `nan`. Changing the setting clears the
cache.

## Usage

``` r
seterr(divide = FALSE)
```

## Arguments

- divide:

  `TRUE` to raise, `FALSE` to return `nan`.

## Value

The previous setting, invisibly.

## Examples

``` r
seterr(TRUE)
tryCatch(Rational(0, 0), error = conditionMessage)
#> [1] "ValueError: Indeterminate 0/0"
seterr(FALSE)
Rational(0, 0)
#> nan
```
