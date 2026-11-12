# Differentiate an expression

`diff(f, x)` differentiates with respect to `x`; `diff(f, x, 3)` gives
the third derivative and `diff(f, x, y)` a mixed partial derivative.
Without variables, `f` must have exactly one free symbol. Derivatives of
unknown functions stay as
[`Derivative()`](https://adamrsophie.github.io/symcompR/reference/Derivative.md)
objects. This is a method for the base R generic
[`base::diff()`](https://rdrr.io/r/base/diff.html).

## Usage

``` r
# S3 method for class 'sym'
diff(x, ..., evaluate = TRUE)
```

## Arguments

- x:

  A `sym` expression.

- ...:

  Variables, optionally each followed by a count.

- evaluate:

  If `FALSE`, return an unevaluated
  [`Derivative()`](https://adamrsophie.github.io/symcompR/reference/Derivative.md).

## Value

A `sym` expression.

## Examples

``` r
x <- Symbol("x"); y <- Symbol("y")
diff(sin(x) * exp(x), x)
#> exp(x)*sin(x) + exp(x)*cos(x)
diff(x^3 * y^2, x, 2, y)
#> 12*x*y
f <- Function("f")
diff(f(x)^2, x)
#> 2*f(x)*Derivative(f(x), x)
```
