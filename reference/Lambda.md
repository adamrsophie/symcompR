# Lambda functions

`Lambda(vars, expr)` is an anonymous symbolic function; apply it with
`apply_lambda()` or with `(...)` via
[`as.function()`](https://rdrr.io/r/base/as.function.html).

## Usage

``` r
Lambda(vars, expr)

apply_lambda(f, ...)
```

## Arguments

- vars:

  A symbol or list of symbols.

- expr:

  The body.

- f:

  A Lambda.

- ...:

  Arguments.

## Value

A `sym` Lambda.
