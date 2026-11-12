# Geometry in three dimensions, rays, planes and conics

`Point3D()`, `Line3D()`, `Segment3D()`, `Ray()` (`Ray2D()`, `Ray3D()`),
`Plane()`, `Ellipse()`, `Parabola()` and `Curve()` extend the plane
geometry of
[`Point()`](https://adamrsophie.github.io/symcompR/reference/Point.md).
[`distance()`](https://adamrsophie.github.io/symcompR/reference/Point.md),
[`intersection()`](https://adamrsophie.github.io/symcompR/reference/Point.md),
[`is_parallel()`](https://adamrsophie.github.io/symcompR/reference/Point.md),
[`area()`](https://adamrsophie.github.io/symcompR/reference/Point.md)
and
[`contains()`](https://adamrsophie.github.io/symcompR/reference/sets.md)
accept them. The accessors follow SymPy's properties:
`direction_ratio()`, `direction_cosine()`, `arbitrary_point()`,
`equation()`, `normal_vector()`, `projection()`, `angle_between()`,
`circumference()`, `eccentricity()`, `foci()`, `major()`,
[`minor()`](https://adamrsophie.github.io/symcompR/reference/sym_det.md)
(an ellipse's minor radius; for a matrix the
[`sym_det()`](https://adamrsophie.github.io/symcompR/reference/sym_det.md)
minor), `focus_distance()`, `encloses_point()`, `focal_length()`,
`vertex()`, `p_parameter()`, `axis_of_symmetry()`, `translate()` and
`geo_length()` (SymPy's `length`). `are_similar()` tests similarity.

## Usage

``` r
Point3D(...)

Line3D(p, q = NULL, direction_ratio = NULL)

Segment3D(p, q)

Ray(p, q)

Ray2D(p, q)

Ray3D(p, q)

direction_ratio(e)

direction_cosine(e)

arbitrary_point(e, t = Symbol("t"))

Plane(p1, p2 = NULL, p3 = NULL, normal_vector = NULL)

normal_vector(e)

equation(e, x = Symbol("x"), y = Symbol("y"), z = Symbol("z"))

projection(e, pt, scalar = FALSE)

angle_between(e1, e2)

Ellipse(center, hradius, vradius)

major(e)

focus_distance(e)

eccentricity(e)

foci(e)

circumference(e)

encloses_point(e, pt)

Parabola(focus, directrix)

focal_length(e)

p_parameter(e)

vertex(e)

axis_of_symmetry(e)

Curve(functions, limits)

translate(e, x = 0, y = 0, z = 0)

geo_length(e)

are_similar(e1, e2)
```

## Arguments

- ...:

  Coordinates, or points.

- p, q:

  Points.

- direction_ratio:

  A direction instead of the second point.

- e, e1, e2:

  Geometric entities.

- t:

  A parameter symbol.

- p1, p2, p3:

  Points of a plane.

- normal_vector:

  The normal of a plane.

- x, y, z:

  Symbols for equations, or offsets for `translate()`.

- pt:

  A point, or for `projection()` of a
  [`CoordSys3D()`](https://adamrsophie.github.io/symcompR/reference/CoordSys3D.md)
  vector the vector to project (`Vector.projection`).

- scalar:

  Return the scalar projection factor instead.

- center, hradius, vradius:

  Ellipse data.

- focus, directrix:

  Parabola data (an axis parallel directrix).

- functions, limits:

  Curve coordinates and `list(t, a, b)`.

## Value

An entity, an expression, a list or a logical.

## Examples

``` r
P <- Point3D(1, 2, 3); Q <- Point3D(4, 6, 8)
L <- Line3D(P, Q)
direction_ratio(L)
#> [[1]]
#> 3
#> 
#> [[2]]
#> 4
#> 
#> [[3]]
#> 5
#> 
distance(L, Point3D(0, 0, 0))
#> 2*sqrt(3)/5
pl <- Plane(Point3D(0, 0, 0), normal_vector = c(1, 1, 1))
intersection(pl, Line3D(Point3D(0, 0, 5), Point3D(1, 1, 5)))
#> [[1]]
#> Point3D(-5/2, -5/2, 5) 
#> 
E <- Ellipse(Point(0, 0), 3, 2)
foci(E)
#> [[1]]
#> Point2D(-sqrt(5), 0) 
#> 
#> [[2]]
#> Point2D(sqrt(5), 0) 
#> 
```
