# Evaluate unevaluated objects

Evaluates
[`Derivative()`](https://adamrsophie.github.io/symcompR/reference/Derivative.md),
`Integral`, `Sum`, `Product` and `Limit` objects inside an expression,
like SymPy's `doit()`.

## Usage

``` r
doit(e, deep = TRUE, ...)
```

## Arguments

- e:

  A `sym` expression.

- deep:

  Evaluate nested objects too.

- ...:

  Hints passed on to the objects, such as `evaluate = FALSE` for the
  characteristic curves of
  [`TendonForceLengthDeGroote2016()`](https://adamrsophie.github.io/symcompR/reference/characteristic-curves.md).

## Value

A `sym` expression.
