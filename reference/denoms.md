# Denominators of an expression

`denoms(e)` returns the factors of the denominators of all
subexpressions of `e`, as SymPy's `denoms()`; with `symbols` only those
that contain one of them.

## Usage

``` r
denoms(e, ...)
```

## Arguments

- e:

  An expression.

- ...:

  Symbols.

## Value

A list of expressions.

## Examples

``` r
x <- Symbol("x")
denoms(1/x + 1/(x + 1))
#> [[1]]
#> x + 1
#> 
#> [[2]]
#> x
#> 
```
