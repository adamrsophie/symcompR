# Parametric regions

`ParametricRegion(definition, ...)` describes a region by its parametric
coordinates; each further argument is a parameter, or
`list(parameter, lower, upper)` giving its bounds. `region$definition`,
`region$parameters`, `region$limits` and `region$dimensions` read it.
`parametric_region_list()` converts a geometric entity (point, segment,
curve, circle, ellipse, polygon) or an
[`ImplicitRegion()`](https://adamrsophie.github.io/symcompR/reference/ImplicitRegion.md)
into a list of parametric regions.

## Usage

``` r
ParametricRegion(definition, ...)

parametric_region_list(reg, parameter = "t")
```

## Arguments

- definition:

  An expression or a list of expressions.

- ...:

  Parameters, or `list(parameter, lower, upper)`.

- reg:

  A geometric entity or an implicit region.

- parameter:

  The name of the parameter (or names, for implicit regions).

## Value

A `ParametricRegion`, or a list of them.

## Examples

``` r
t <- Symbol("t")
p <- ParametricRegion(list(t, t^2), list(t, 0, 1))
p$limits
#> {t: (0, 1)} 
parametric_region_list(Segment(Point(0, 0), Point(1, 1)))
#> [[1]]
#> ParametricRegion((t, t), (t, 0, 1))
#> 
```
