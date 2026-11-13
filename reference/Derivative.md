# Unevaluated derivative

`Derivative(f, x, ...)` represents a derivative without computing it.
Use [`doit()`](https://adamrsophie.github.io/symcompR/reference/doit.md)
to evaluate.

## Usage

``` r
Derivative(f, ..., evaluate = FALSE)
```

## Arguments

- f:

  Expression.

- ...:

  Variables and counts, as for
  [`diff()`](https://rdrr.io/r/base/diff.html).

- evaluate:

  If `TRUE`, compute the derivative, as
  [`diff()`](https://rdrr.io/r/base/diff.html) does.

## Value

A `sym` Derivative.
