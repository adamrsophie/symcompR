# Unevaluated derivative

`Derivative(f, x, ...)` represents a derivative without computing it.
Use [`doit()`](https://adamrsophie.github.io/symcompR/reference/doit.md)
to evaluate.

## Usage

``` r
Derivative(f, ...)
```

## Arguments

- f:

  Expression.

- ...:

  Variables and counts, as for
  [`diff()`](https://rdrr.io/r/base/diff.html).

## Value

A `sym` Derivative.
