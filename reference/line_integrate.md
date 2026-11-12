# Line integrals

`line_integrate(f, curve, vars)` integrates `f` along a parametric curve
given as `list(param, a, b, x(param), y(param), ...)`.

## Usage

``` r
line_integrate(f, curve, vars)
```

## Arguments

- f:

  An expression.

- curve:

  A list `list(param, a, b, coord1, coord2, ...)`.

- vars:

  The coordinate variables.

## Value

An expression.

## Examples

``` r
t <- Symbol("t")
line_integrate(1, list(t, 0, 1, t, t), list(Symbol("x"), Symbol("y")))
#> sqrt(2)
```
