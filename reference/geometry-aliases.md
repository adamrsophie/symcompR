# Geometry aliases and helpers

`Point2D()`, `Line2D()` and `Segment2D()` are the 2D constructors.
`rad()` and `deg()` convert angles. `RegularPolygon(center, r, n)` gives
the vertices of a regular polygon as a
[`Polygon()`](https://adamrsophie.github.io/symcompR/reference/Point.md).
`convex_hull()` returns the hull of points (Andrew's monotone chain),
and `closest_points()` and `farthest_points()` the pairs at minimal and
maximal distance.

## Usage

``` r
Point2D(...)

Line2D(p, q)

Segment2D(p, q)

rad(x)

deg(x)

RegularPolygon(center, r, n)

convex_hull(...)

closest_points(...)

farthest_points(...)
```

## Arguments

- ...:

  Coordinates or points.

- p, q:

  Points.

- x:

  An angle.

- center, r, n:

  Regular polygon data.

## Value

A geometric entity, expression or list.

## Examples

``` r
convex_hull(Point(0, 0), Point(1, 1), Point(2, 0), Point(1, 0.5), Point(1, -1))
#> Polygon(Point2D(0, 0), Point2D(1, -1), Point2D(2, 0), Point2D(1, 1)) 
```
