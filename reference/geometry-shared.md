# Measurements and relations of geometric entities

These functions are SymPy's entity methods and properties that several
classes share: `distance()`, `intersection()` (of any number of
entities, as `sympy.geometry.intersection`), `area()`, `centroid()` (of
points, segments or polygons, as `sympy.geometry.centroid`),
`geo_length()` (SymPy's `length`), `equation()`, `projection()`,
`angle_between()`, `is_parallel()`, `is_perpendicular()`,
`normal_vector()`, `arbitrary_point()`, `parameter_value()`,
`plot_interval()`, `random_point()`, `encloses()`, `equals()`,
`are_similar()`, `are_coplanar()`, `convex_hull()`, `closest_points()`,
`farthest_points()`, `line_equation()` and `circle_equation()`.

## Usage

``` r
distance(p, q)

intersection(..., pairwise = FALSE)

area(e)

centroid(...)

geo_length(e)

equation(e, x = NULL, y = NULL, z = NULL, slope = NULL)

line_equation(e, x = "x", y = "y")

circle_equation(e, x = "x", y = "y")

projection(e, pt, scalar = FALSE)

angle_between(e1, e2)

is_parallel(e1, e2)

is_perpendicular(e1, e2)

normal_vector(e)

arbitrary_point(e, parameter = "t", v = NULL)

parameter_value(e, other, t, v = NULL)

plot_interval(e, parameter = "t")

random_point(e, seed = NULL)

encloses(e, o)

equals(a, b)

are_similar(e1, e2)

are_coplanar(...)

convex_hull(..., polygon = TRUE)

closest_points(...)

farthest_points(...)
```

## Arguments

- ...:

  Entities or points.

- pairwise:

  Intersect every pair of the entities.

- x, y, z:

  Names or symbols of the variables.

- slope:

  The slope of the major axis of a rotated ellipse.

- pt:

  A point.

- scalar:

  For a vector projection, return the scalar component.

- parameter:

  Parameter name or symbol.

- v:

  The second parameter of a plane.

- other:

  A point.

- t:

  The parameter symbol.

- seed:

  A seed for the random point.

- a, b, e, e1, e2, p, q, o:

  Entities.

- polygon:

  Return a polygon (or segment or point) rather than the upper and lower
  hulls.

## Value

An entity, an expression, a logical value or a list.

## Examples

``` r
intersection(Circle(Point(0, 0), 1), Line(Point(0, 0), Point(1, 0)))
#> [[1]]
#> Point2D(-1, 0) 
#> 
#> [[2]]
#> Point2D(1, 0) 
#> 
area(Circle(Point(0, 0), 2))
#> 4*pi
centroid(Point(0, 0), Point(2, 0), Point(1, 3))
#> Point2D(1, 1) 
convex_hull(Point(0, 0), Point(1, 1), Point(2, 0), Point(1, 0.5))
#> Triangle(Point2D(0, 0), Point2D(2, 0), Point2D(1, 1)) 
```
