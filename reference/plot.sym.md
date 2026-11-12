# Plot symbolic expressions with base graphics

`plot(expr, range)` draws a function of one variable, with `range` given
as `list(x, a, b)` (default `-10` to `10`). `plot_parametric()`,
`plot3d()` (surface) and `plot_implicit()` (zero contour) are the
analogues of SymPy's plotting functions. Extra arguments go to the base
graphics functions.

## Usage

``` r
# S3 method for class 'sym'
plot(x, range = NULL, ..., n = 400)

plot_parametric(fx, fy, range, ..., n = 400)

plot3d(e, range1, range2, ..., n = 50)

plot_implicit(e, range1, range2, ..., n = 200)
```

## Arguments

- x:

  A `sym` expression, for
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html).

- range:

  `list(var, a, b)`.

- ...:

  Graphical parameters.

- n:

  Number of sample points.

- fx, fy:

  Coordinate expressions.

- e:

  An expression.

- range1, range2:

  Ranges for two variables.

## Value

The sampled data, invisibly.

## Examples

``` r
x <- Symbol("x")
if (interactive()) plot(sin(x) / x, list(x, -10, 10))
```
