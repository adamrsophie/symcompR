# Geometry

This vignette covers `sympy.geometry` (exact Euclidean geometry in two
and three dimensions) and quaternions from `sympy.algebras`. Vector
calculus, tensors, differential geometry, Lie algebras, and categories
have their own vignettes.

``` r

library(symcompR)
x <- Symbol("x"); y <- Symbol("y")
t <- Symbol("t"); theta <- Symbol("theta")
```

## Plane geometry

Points, lines, segments, polygons, circles, ellipses, and parabolas are
exact: coordinates are symbolic numbers, so distances and areas come out
as closed forms.

``` r

A <- Point(0, 0); B <- Point(3, 4)
distance(A, B)
#> 5
midpoint(Point(1, 2), Point(3, 6))
#> Point2D(2, 4)
translate(Point(1, 2), 1, 1)
#> Point2D(2, 3)
L <- Line(Point(0, 0), Point(1, 2))
slope(L)
#> 2
line_equation(L)
#> -2*x + y
is_perpendicular(Line(Point(0, 0), Point(1, 1)), Line(Point(0, 2), Point(2, 0)))
#> [1] TRUE
is_collinear(Point(0, 0), Point(1, 1), Point(2, 2))
#> [1] TRUE
intersection(Line(Point(0, 0), Point(1, 1)), Line(Point(0, 2), Point(2, 0)))
#> [[1]]
#> Point2D(1, 1)
```

Triangles and polygons:

``` r

T3 <- Triangle(Point(0, 0), Point(4, 0), Point(0, 3))
area(T3)
#> 6
perimeter(T3)
#> 12
centroid(T3)
#> Point2D(4/3, 1)
area(Polygon(Point(0, 0), Point(2, 0), Point(2, 1), Point(0, 1)))
#> 2
area(RegularPolygon(Point(0, 0), 1, 6))
#> 3*sqrt(3)/2
are_similar(Triangle(Point(0, 0), Point(1, 0), Point(0, 1)),
            Triangle(Point(0, 0), Point(2, 0), Point(0, 2)))
#> [1] TRUE
convex_hull(Point(0, 0), Point(1, 1), Point(2, 0), Point(1, Rational(1, 2)), Point(1, -1))
#> Polygon(Point2D(0, 0), Point2D(1, -1), Point2D(2, 0), Point2D(1, 1))
```

Conic sections:

``` r

C <- Circle(Point(0, 0), 5)
circle_equation(C)
#> x^2 + y^2 - 25
area(C)
#> 25*pi
intersection(C, Line(Point(0, 0), Point(1, 0)))
#> [[1]]
#> Point2D(-5, 0) 
#> 
#> [[2]]
#> Point2D(5, 0)
E <- Ellipse(Point(0, 0), 3, 2)
foci(E)
#> [[1]]
#> Point2D(-sqrt(5), 0) 
#> 
#> [[2]]
#> Point2D(sqrt(5), 0)
eccentricity(E)
#> sqrt(5)/3
major(E)
#> 3
area(E)
#> 6*pi
circumference(E)
#> 12*elliptic_e(5/9)
P <- Parabola(Point(0, 1), Line(Point(-1, -1), Point(1, -1)))
vertex(P)
#> Point2D(0, 0)
equation(P)
#> -x^2 + 4*y
```

The circumference of an ellipse is a complete elliptic integral of the
second kind. Segments have a length and a parametrization:

``` r

S <- Segment(Point(0, 0), Point(3, 4))
geo_length(S)
#> 5
arbitrary_point(S, t)
#> Point2D(3*t, 4*t)
```

## Solid geometry

``` r

L3 <- Line3D(Point3D(1, 2, 3), Point3D(4, 6, 8))
direction_ratio(L3)
#> [[1]]
#> 3
#> 
#> [[2]]
#> 4
#> 
#> [[3]]
#> 5
distance(L3, Point3D(0, 0, 0))
#> 2*sqrt(3)/5
pl <- Plane(Point3D(0, 0, 0), normal_vector = c(1, 1, 1))
equation(pl)
#> x + y + z
intersection(pl, Line3D(Point3D(0, 0, 5), Point3D(1, 1, 5)))
#> [[1]]
#> Point3D(-5/2, -5/2, 5)
normal_vector(Plane(Point3D(1, 0, 0), Point3D(0, 1, 0), Point3D(0, 0, 1)))
#> [[1]]
#> 1
#> 
#> [[2]]
#> 1
#> 
#> [[3]]
#> 1
distance(Plane(Point3D(0, 0, 0), normal_vector = c(0, 0, 1)), Point3D(1, 2, 3))
#> 3
angle_between(Plane(Point3D(0, 0, 0), normal_vector = c(1, 0, 0)),
              Plane(Point3D(0, 0, 0), normal_vector = c(1, 1, 0)))
#> pi/4
```

### Rays, curves, and more properties

[`Ray()`](https://adamrsophie.github.io/symcompR/reference/geometry-3d.md)
starts at a point and passes through another.
[`Ray2D()`](https://adamrsophie.github.io/symcompR/reference/geometry-3d.md),
[`Ray3D()`](https://adamrsophie.github.io/symcompR/reference/geometry-3d.md),
[`Point2D()`](https://adamrsophie.github.io/symcompR/reference/geometry-aliases.md),
[`Line2D()`](https://adamrsophie.github.io/symcompR/reference/geometry-aliases.md),
[`Segment2D()`](https://adamrsophie.github.io/symcompR/reference/geometry-aliases.md),
and
[`Segment3D()`](https://adamrsophie.github.io/symcompR/reference/geometry-3d.md)
are the dimension specific names.

``` r

r <- Ray(Point(0, 0), Point(1, 1))
contains(r, Point(2, 2))
#> [1] TRUE
distance(Point(2, 1), Ray2D(Point(0, 0), Point(1, 0)))
#> 1
distance(Point3D(0, 0, 0), Ray3D(Point3D(0, 0, 2), Point3D(1, 0, 5)))
#> 2
intersection(Ray2D(Point(0, 0), Point(1, 0)), Ray2D(Point(3, 0), Point(2, 0)))
#> [[1]]
#> Segment2D(Point2D(0, 0), Point2D(3, 0))
Point2D(1, 2)
#> Point2D(1, 2)
Line2D(Point(0, 0), Point(1, 1))
#> Line2D(Point2D(0, 0), Point2D(1, 1))
geo_length(Segment2D(Point(0, 0), Point(1, 1)))
#> sqrt(2)
geo_length(Segment3D(Point3D(0, 0, 0), Point3D(3, 4, 5)))
#> 5*sqrt(2)
direction_cosine(Line3D(Point3D(0, 0, 0), Point3D(3, 4, 5)))
#> [[1]]
#> 3*sqrt(2)/10
#> 
#> [[2]]
#> 2*sqrt(2)/5
#> 
#> [[3]]
#> sqrt(2)/2
projection(Line(Point(0, 0), Point(1, 1)), Point(2, 0))
#> Point()
is_parallel(Line(Point(0, 0), Point(1, 1)), Line(Point(0, 1), Point(1, 2)))
#> [1] TRUE
closest_points(Point(0, 0), Point(3, 0), Point(1, 1))
#> [[1]]
#> [[1]][[1]]
#> Point2D(0, 0) 
#> 
#> [[1]][[2]]
#> Point2D(1, 1)
farthest_points(Point(0, 0), Point(3, 0), Point(1, 1))
#> [[1]]
#> [[1]][[1]]
#> Point2D(0, 0) 
#> 
#> [[1]][[2]]
#> Point2D(3, 0)
```

Ellipse and parabola properties:

``` r

encloses_point(E, Point(1, 1))
#> [1] TRUE
focus_distance(E)
#> sqrt(5)
P2 <- Parabola(Point(0, 3), Line(Point(-1, 1), Point(1, 1)))
axis_of_symmetry(P2)
#> Line2D(Point2D(0, 3), Point2D(0, 4))
focal_length(P2)
#> 1
p_parameter(P2)
#> 1
```

[`Curve()`](https://adamrsophie.github.io/symcompR/reference/geometry-3d.md)
is a parametric curve:

``` r

cv <- Curve(list(t, t^2), list(t, 0, 1))
arbitrary_point(cv)
#> Point2D(t, t^2)
```

## Quaternions

``` r

q1 <- Quaternion(1, 2, 3, 4)
q1 + Quaternion(5, 6, 7, 8)
#> 6 + 8*i + 10*j + 12*k
q1 * quat_inverse(q1)
#> 1 + 0*i + 0*j + 0*k
q1 * quat_conjugate(q1)
#> 30 + 0*i + 0*j + 0*k
quat_norm(q1)
#> sqrt(30)
```

A unit quaternion represents a rotation:

``` r

trigsimp(quat_to_rotation_matrix(Quaternion(cos(theta / 2), 0, 0, sin(theta / 2))))
#> Matrix(3, 3)
#> [cos(theta), -sin(theta), 0]
#> [sin(theta), cos(theta) , 0]
#> [0         , 0          , 1]
```
