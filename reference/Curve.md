# Curves and parabolas

`Curve()` is a parametric curve `list(x(t), y(t))` over `list(t, a, b)`;
`Parabola()` the parabola with a focus and a directrix. The fields are
`$functions`, `$limits` and `$parameter`, and `$focus` and `$directrix`.
The functions are SymPy's methods: `focal_length()`, `p_parameter()`,
`vertex()`, `axis_of_symmetry()`, and the shared
[`equation()`](https://adamrsophie.github.io/symcompR/reference/geometry-shared.md),
[`intersection()`](https://adamrsophie.github.io/symcompR/reference/geometry-shared.md),
[`arbitrary_point()`](https://adamrsophie.github.io/symcompR/reference/geometry-shared.md)
and
[`geo_length()`](https://adamrsophie.github.io/symcompR/reference/geometry-shared.md).

## Usage

``` r
Curve(functions, limits)

Parabola(focus = NULL, directrix = NULL)

axis_of_symmetry(e)

focal_length(e)

p_parameter(e)

vertex(e)
```

## Arguments

- functions:

  A list of the coordinate functions.

- limits:

  A list `list(t, a, b)`.

- focus:

  The focus.

- directrix:

  The directrix, a line.

- e:

  A curve or a parabola.

## Value

An entity, an expression or a list.

## Examples

``` r
t <- Symbol("t")
C <- Curve(list(t, t^2), list(t, 0, 1))
arbitrary_point(C)
#> Point2D(t, t^2) 
P <- Parabola(Point(0, 0), Line(Point(5, 8), Point(7, 8)))
vertex(P)
#> Point2D(0, 4) 
equation(P)
#> -x^2 - 16*y + 64
```
