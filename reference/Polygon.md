# Polygons and triangles

`Polygon()` builds a polygon from its vertices as SymPy does: repeated
and collinear vertices are removed, so that three vertices give a
`Triangle`, two a segment and one a point; `n` gives a regular polygon.
`Triangle()` also accepts `sss`, `asa` or `sas` (side lengths and angles
in degrees). `RegularPolygon()` takes a center, a radius, the number of
sides and a rotation. The vertices are in `$vertices`.

## Usage

``` r
Polygon(..., n = 0)

Triangle(...)

RegularPolygon(c, r, n, rot = 0)

angles(p)

perimeter(p)

first_moment_of_area(p, point = NULL)

sides(p)

cut_section(p, line)

apothem(p)

interior_angle(p)

exterior_angle(p)

rotation(p)

spin(p, angle)

is_similar(t1, t2)

is_equilateral(t)

is_isosceles(t)

is_scalene(t)

is_right(t)

altitudes(t)

orthocenter(t)

circumcenter(t)

circumradius(t)

circumcircle(t)

incenter(t)

inradius(t)

incircle(t)

exradii(t)

excenters(t)

medians(t)

medial(t)

nine_point_circle(t)

eulerline(t)
```

## Arguments

- ...:

  Vertices; for `Triangle()` also `sss`, `asa` or `sas`.

- n:

  The number of sides of a regular polygon.

- c, r:

  The center and radius of a regular polygon.

- rot, angle:

  Rotation angles.

- p, t, t1, t2:

  Polygons, triangles or other entities.

- point:

  A reference point.

- line:

  A line.

## Value

An entity, an expression, a logical value or a list.

## Details

The functions are SymPy's methods and properties: `sides()`, `angles()`,
`perimeter()`,
[`is_convex()`](https://adamrsophie.github.io/symcompR/reference/calculus-utils.md),
[`encloses()`](https://adamrsophie.github.io/symcompR/reference/geometry-shared.md),
`cut_section()`, `first_moment_of_area()`, the triangle centers
`circumcenter()`, `incenter()`, `orthocenter()`, `excenters()` and the
related `circumradius()`, `circumcircle()`, `inradius()`, `incircle()`,
`exradii()`, `altitudes()`, `medians()`, `medial()`,
`nine_point_circle()`, `eulerline()`, `is_right()`, `is_equilateral()`,
`is_isosceles()`, `is_scalene()`, `is_similar()`, and for regular
polygons `apothem()`, `interior_angle()`, `exterior_angle()`,
`rotation()` and `spin()` (which returns the spun polygon). Results
keyed by vertices or sides are named lists.

## Examples

``` r
p <- Polygon(Point(0, 0), Point(1, 0), Point(1, 1), Point(0, 1))
area(p)
#> 1
perimeter(p)
#> 4
t <- Triangle(Point(0, 0), Point(4, 0), Point(0, 3))
circumcenter(t)
#> Point2D(2, 3/2) 
incenter(t)
#> Point2D(1, 1) 
Polygon(Point(0, 0), Point(1, 0), Point(2, 0), Point(1, 1))
#> Triangle(Point2D(0, 0), Point2D(2, 0), Point2D(1, 1)) 
```
