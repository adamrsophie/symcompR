# Planes

`Plane()` builds a plane through three points or through a point with a
`normal_vector`, as SymPy does. The point and the normal vector are in
`$p1` and `$normal_vector`. The functions are SymPy's methods:
`is_coplanar()`, `parallel_plane()`, `perpendicular_plane()`,
`projection_line()` and the shared
[`distance()`](https://adamrsophie.github.io/symcompR/reference/geometry-shared.md),
[`intersection()`](https://adamrsophie.github.io/symcompR/reference/geometry-shared.md),
[`projection()`](https://adamrsophie.github.io/symcompR/reference/geometry-shared.md),
[`angle_between()`](https://adamrsophie.github.io/symcompR/reference/geometry-shared.md),
[`is_parallel()`](https://adamrsophie.github.io/symcompR/reference/geometry-shared.md),
[`is_perpendicular()`](https://adamrsophie.github.io/symcompR/reference/geometry-shared.md),
[`perpendicular_line()`](https://adamrsophie.github.io/symcompR/reference/Line.md)
and
[`arbitrary_point()`](https://adamrsophie.github.io/symcompR/reference/geometry-shared.md).

## Usage

``` r
Plane(p1, a = NULL, b = NULL, normal_vector = NULL, evaluate = TRUE)

is_coplanar(e, o)

parallel_plane(pl, pt)

perpendicular_plane(pl, ...)

projection_line(pl, line)
```

## Arguments

- p1:

  A point of the plane.

- a, b:

  Two more points, when `normal_vector` is not given.

- normal_vector:

  The normal vector, a list or vector of length 3.

- evaluate:

  Make a point of the normal vector (rationalizing floats).

- e, pl:

  A plane.

- o, pt, line:

  An entity, a point or a line.

- ...:

  Points (for `perpendicular_plane()`) or planes.

## Value

An entity, an expression, a logical value or a list.

## Examples

``` r
p <- Plane(Point3D(1, 1, 1), normal_vector = c(1, 4, 7))
equation(p)
#> x + 4*y + 7*z - 12
distance(p, Point3D(1, 2, 3))
#> 3*sqrt(66)/11
intersection(Plane(Point3D(0, 0, 0), normal_vector = c(1, 0, 0)),
             Plane(Point3D(0, 0, 0), normal_vector = c(0, 1, 0)))
#> [[1]]
#> Line3D(Point3D(0, 0, 0), Point3D(0, 0, 1)) 
#> 
```
