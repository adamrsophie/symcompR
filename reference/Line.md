# Lines, rays and segments

`Line()`, `Ray()` and `Segment()` build linear entities through two
points, as SymPy does; the dimension of the points gives `Line2D`,
`Line3D` and so on. A 2D line may also be given by a point and `slope`,
or by a linear equation in `x` and `y`; a ray by a point and `angle`; a
3D line or ray by a point and `direction_ratio`. A segment between equal
points is that point.

## Usage

``` r
Line(p1, p2 = NULL, slope = NULL, direction_ratio = NULL, x = NULL, y = NULL)

Line2D(p1, pt = NULL, slope = NULL)

Line3D(p1, pt = NULL, direction_ratio = NULL)

Ray(p1, p2 = NULL, angle = NULL, direction_ratio = NULL)

Ray2D(p1, pt = NULL, angle = NULL)

Ray3D(p1, pt = NULL, direction_ratio = NULL)

Segment(p1, p2)

Segment2D(p1, p2)

Segment3D(p1, p2)

direction(l)

smallest_angle_between(l1, l2)

are_concurrent(...)

perpendicular_line(l, p)

parallel_line(l, p)

perpendicular_segment(l, p)

perpendicular_bisector(l, p = NULL)

slope(l)

geo_coefficients(l)

xdirection(l)

ydirection(l)

zdirection(l)

closing_angle(r1, r2)

bisectors(l, l2 = NULL)

bounds(l)
```

## Arguments

- p1:

  A point, a linear entity, or (for `Line()`) an equation.

- p2, pt:

  A second point.

- slope:

  The slope of a 2D line.

- direction_ratio:

  The direction of a 3D line or ray.

- x, y:

  Names or symbols of the variables of an equation.

- angle:

  The angle of a 2D ray.

- l, l1, l2, r1, r2:

  Linear entities.

- ...:

  Lines.

- p:

  A point.

## Value

An entity, an expression, a logical value or a list.

## Details

The other functions are the SymPy methods: `direction()`, `slope()`,
`geo_coefficients()` (SymPy's `coefficients`), `perpendicular_line()`,
`parallel_line()`, `perpendicular_segment()`,
`perpendicular_bisector()`, `smallest_angle_between()`,
`are_concurrent()`, `bisectors()`, `xdirection()`, `ydirection()`,
`zdirection()` and `closing_angle()`. The fields `$p1`, `$p2` and, for a
ray, `$source` hold the defining points.

## Examples

``` r
l <- Line(Point(0, 0), Point(1, 1))
slope(l)
#> 1
equation(l)
#> -x + y
intersection(l, Line(Point(0, 2), Point(2, 0)))
#> [[1]]
#> Point2D(1, 1) 
#> 
Segment(Point(1, 1), Point(1, 1))
#> Point2D(1, 1) 
Line(Point(0, 0), slope = 2)
#> Line2D(Point2D(0, 0), Point2D(1, 2)) 
```
