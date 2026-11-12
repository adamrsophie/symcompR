# Evaluate unevaluated objects

Evaluates
[`Derivative()`](https://adamrsophie.github.io/symcompR/reference/Derivative.md),
`Integral`, `Sum`, `Product` and `Limit` objects inside an expression,
like SymPy's `doit()`.

## Usage

``` r
doit(e, deep = TRUE)
```

## Arguments

- e:

  A `sym` expression.

- deep:

  Evaluate nested objects too.

## Value

A `sym` expression.
