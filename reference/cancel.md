# Cancel common factors of a rational function

Puts `e` in the form `p/q` with `p` and `q` expanded polynomials with no
common factor, like SymPy's `cancel()`.

## Usage

``` r
cancel(e, ...)
```

## Arguments

- e:

  An expression.

- ...:

  Generators.

## Value

An expression.

## Examples

``` r
x <- Symbol("x")
cancel((x^2 - 1) / (x - 1))
#> x + 1
```
