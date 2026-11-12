# Plane geometry

`Point()`, `Line()`, `Segment()`, `Circle()`, `Polygon()` and
`Triangle()` build geometric entities with exact coordinates.
`distance()` measures between points, `midpoint()`, `area()`,
`perimeter()`, `centroid()`, `slope()`, `is_parallel()`,
`is_perpendicular()`, `is_collinear()`, `intersection()`,
`line_equation()` and `circle_equation()` follow SymPy.

## Usage

``` r
Point(...)

distance(p, q)

midpoint(p, q = NULL)

Line(p, q)

Segment(p, q)

slope(l)

line_equation(l, x = Symbol("x"), y = Symbol("y"))

is_parallel(l, m)

is_perpendicular(l, m)

is_collinear(...)

Circle(center, radius)

circle_equation(l, x = Symbol("x"), y = Symbol("y"))

Polygon(...)

Triangle(a, b, P)

area(p)

perimeter(p)

centroid(p)

intersection(a, b)
```

## Arguments

- ...:

  Coordinates for `Point()`; vertices for `Polygon()`.

- p, q, a, b, P:

  Points or entities.

- l, m:

  Lines.

- x, y:

  Symbols for equations.

- center, radius:

  Circle data.

## Value

A geometric entity, expression or list.

## Examples

``` r
A <- Point(0, 0); B <- Point(3, 4)
distance(A, B)
#> 5
area(Triangle(Point(0, 0), Point(4, 0), Point(0, 3)))
#> 6
intersection(Line(Point(0, 0), Point(1, 1)), Line(Point(0, 2), Point(2, 0)))
#> [[1]]
#> Point2D(1, 1) 
#> 
```
