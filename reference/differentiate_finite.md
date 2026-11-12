# Finite difference approximation of a derivative

`differentiate_finite(expr, ...)` replaces the derivative of `expr` with
respect to the variables in `...` by a finite difference of `expr`
itself, so the product rule does not hold for it (as SymPy's
`differentiate_finite()`). `points` is the step, or a list of points.

## Usage

``` r
differentiate_finite(expr, ..., points = 1, x0 = NULL, wrt = NULL)
```

## Arguments

- expr:

  An expression.

- ...:

  Variables, each optionally followed by an order.

- points:

  A step or a list of points.

- x0:

  The point at which the derivative is approximated.

- wrt:

  Unused; kept for SymPy compatibility.

## Value

An expression.

## Examples

``` r
x <- Symbol("x"); h <- Symbol("h"); f <- Function("f"); g <- Function("g")
differentiate_finite(f(x)*g(x), x, points = list(x - h, x + h))
#> -f(-h + x)*g(-h + x)/(2*h) + f(h + x)*g(h + x)/(2*h)
```
