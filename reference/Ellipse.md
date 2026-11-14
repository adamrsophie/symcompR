# Ellipses and circles

`Ellipse()` builds an ellipse from its center and two of `hradius`,
`vradius` and `eccentricity`; equal radii give a `Circle`, a zero radius
a segment. `Circle()` takes a center and a radius, three points, or an
equation in `x` and `y`. The fields `$center`, `$hradius`, `$vradius`
and, for a circle, `$radius` hold the defining data.

## Usage

``` r
Ellipse(center = NULL, hradius = NULL, vradius = NULL, eccentricity = NULL)

Circle(..., x = "x", y = "y", evaluate = TRUE)

major(e)

foci(e)

focus_distance(e)

eccentricity(e)

circumference(e)

periapsis(e)

apoapsis(e)

semilatus_rectum(e)

encloses_point(e, p)

is_tangent(e, o)

tangent_lines(e, p)

normal_lines(e, p, ...)

evolute(e, x = "x", y = "y")

auxiliary_circle(e)

director_circle(e)

second_moment_of_area(e, point = NULL)

polar_second_moment_of_area(e)

section_modulus(e, point = NULL)
```

## Arguments

- center:

  The center.

- hradius, vradius:

  The horizontal and vertical radii.

- eccentricity:

  The eccentricity.

- ...:

  For `Circle()`: a center and a radius, three points, or an equation;
  for `normal_lines()`, `prec`.

- x, y:

  Names or symbols of the variables of an equation.

- evaluate:

  Replace floats in the points by rationals.

- e:

  An ellipse or circle.

- o, p:

  An entity or a point.

- point:

  A reference point.

## Value

An entity, an expression, a logical value or a list.

## Details

The functions are SymPy's methods and properties: `major()`,
[`minor()`](https://adamrsophie.github.io/symcompR/reference/sym_det.md),
`foci()`, `focus_distance()`, `eccentricity()`, `circumference()`,
`periapsis()`, `apoapsis()`, `semilatus_rectum()`, `encloses_point()`,
`is_tangent()`, `tangent_lines()`, `normal_lines()`, `evolute()`,
`auxiliary_circle()`, `director_circle()`, `second_moment_of_area()`,
`polar_second_moment_of_area()` and `section_modulus()`.

## Examples

``` r
e <- Ellipse(Point(0, 0), 3, 2)
foci(e)
#> [[1]]
#> Point2D(-sqrt(5), 0) 
#> 
#> [[2]]
#> Point2D(sqrt(5), 0) 
#> 
eccentricity(e)
#> sqrt(5)/3
equation(Circle(Point(1, 2), 3))
#> (x - 1)^2 + (y - 2)^2 - 9
Circle(Point(0, 0), Point(1, 0), Point(0, 1))
#> Circle(Point2D(1/2, 1/2), sqrt(2)/2) 
tangent_lines(Circle(Point(0, 0), 1), Point(2, 0))
#> [[1]]
#> Line2D(Point2D(2, 0), Point2D(1/2, -sqrt(3)/2)) 
#> 
#> [[2]]
#> Line2D(Point2D(2, 0), Point2D(1/2, sqrt(3)/2)) 
#> 
```
