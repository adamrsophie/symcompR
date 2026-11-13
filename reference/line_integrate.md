# Line integrals

`line_integrate(f, curve, vars)` integrates `f` along a parametric
curve, a
[`Curve()`](https://adamrsophie.github.io/symcompR/reference/geometry-3d.md)
or a list `list(param, a, b, x(param), y(param), ...)`.

## Usage

``` r
line_integrate(f, curve, vars)
```

## Arguments

- f:

  An expression.

- curve:

  A
  [`Curve()`](https://adamrsophie.github.io/symcompR/reference/geometry-3d.md),
  or a list `list(param, a, b, coord1, coord2, ...)`.

- vars:

  The coordinate variables.

## Value

An expression.

## Examples

``` r
t <- Symbol("t")
line_integrate(1, list(t, 0, 1, t, t), list(Symbol("x"), Symbol("y")))
#> sqrt(2)
line_integrate(Symbol("x"), Curve(list(t, t), list(t, 0, 1)), list(Symbol("x"), Symbol("y")))
#> sqrt(2)/2
```
