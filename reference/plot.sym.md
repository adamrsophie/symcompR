# Plot symbolic expressions with base graphics

`plot(expr, range)` draws a function of one variable, with `range` given
as `list(x, a, b)` (default `-10` to `10`). `plot_parametric()`,
`plot3d()` (surface), `plot_implicit()` (zero contour),
`plot_contour()`, `plot3d_parametric_line()` and
`plot3d_parametric_surface()` are the analogues of SymPy's plotting
functions. Extra arguments go to the base graphics functions.

## Usage

``` r
# S3 method for class 'sym'
plot(x, range = NULL, ..., n = 400)

plot_parametric(fx, fy, range, ..., n = 400)

plot3d(e, range1, range2, ..., n = 50)

plot_implicit(e, range1, range2, ..., n = 200)

plot_contour(e, range1, range2, ..., levels = NULL, n = 100)

plot3d_parametric_line(fx, fy, fz, range, ..., n = 300)

plot3d_parametric_surface(fx, fy, fz, range1, range2, ..., n = 30)

PlotGrid(nrows, ncols, ...)
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

- levels:

  Contour levels; by default chosen by
  [`graphics::contour()`](https://rdrr.io/r/graphics/contour.html).

- fz:

  The third coordinate expression.

- nrows, ncols:

  The grid shape.

## Value

The sampled data, invisibly.

## Details

`PlotGrid(nrows, ncols, ...)` draws several plots in one figure; each
argument in `...` is a function of no arguments that draws one plot, as
SymPy's `PlotGrid` takes plot objects.

## Examples

``` r
x <- Symbol("x")
if (interactive()) plot(sin(x) / x, list(x, -10, 10))
```
