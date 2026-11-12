# Combine a sum of fractions

`together()` rewrites a sum of rational terms over a common denominator
without expanding, like SymPy's `together()`.

## Usage

``` r
together(e, deep = FALSE, fraction = TRUE)
```

## Arguments

- e:

  An expression.

- deep:

  Also combine inside function arguments.

- fraction:

  Put the result over a common denominator.

## Value

An expression.

## Examples

``` r
x <- Symbol("x"); y <- Symbol("y")
together(1/x + 1/y)
#> (x + y)/(x*y)
```
