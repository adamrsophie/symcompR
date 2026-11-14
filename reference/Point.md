# Points

`Point()` builds a point of any dimension, as SymPy's `Point`: two
coordinates give a `Point2D`, three a `Point3D`. Float coordinates are
replaced by rationals unless `evaluate = FALSE`; `dim` pads with zeros
or drops trailing zero coordinates. Points support `+`, `-`, `*` and `/`
(each coordinate is simplified, as in SymPy) and `==`, which compares
the coordinates structurally.

## Usage

``` r
Point(..., evaluate = TRUE, dim = NULL, on_morph = "ignore")

Point2D(..., evaluate = TRUE)

Point3D(..., evaluate = TRUE)

point_x(p)

point_y(p)

point_z(p)

ambient_dimension(p)

affine_rank(...)

is_collinear(...)

are_collinear(...)

is_concyclic(...)

is_scalar_multiple(p, q)

midpoint(p, q = NULL)

origin(p)

orthogonal_direction(p)

project(a, b)

taxicab_distance(p, q)

canberra_distance(p, q)

unit(p)

geo_transform(p, matrix)

direction_ratio(p, q = NULL)

direction_cosine(p, q = NULL)
```

## Arguments

- ...:

  Coordinates, a list of coordinates, or a point; for `affine_rank()`,
  `is_collinear()`, `are_collinear()` and `is_concyclic()`, points.

- evaluate:

  Replace floats by rationals.

- dim:

  The dimension of the point.

- on_morph:

  What to do when `dim` changes the dimension: `"ignore"`, `"warn"` or
  `"error"`.

- p, q, a, b:

  Points.

- matrix:

  A 3 by 3 (or, in 3D, 4 by 4) transformation matrix.

## Value

A point, an expression, a logical value or a list.

## Details

The functions below are SymPy's `Point` methods and properties:
`point_x()`, `point_y()` and `point_z()` (the `x`, `y` and `z`
properties), `ambient_dimension()`, `affine_rank()`, `is_collinear()`,
`are_collinear()`, `is_concyclic()`, `is_scalar_multiple()`,
`midpoint()`, `origin()`, `orthogonal_direction()`, `project()`,
`taxicab_distance()`, `canberra_distance()`, `unit()`,
`direction_ratio()`, `direction_cosine()` and `geo_transform()` (SymPy's
`transform`).

## Examples

``` r
A <- Point(0, 0); B <- Point(3, 4)
distance(A, B)
#> 5
A + B
#> Point2D(3, 4) 
Point(0.5, 0.25)
#> Point2D(1/2, 1/4) 
is_collinear(Point(0, 0), Point(1, 1), Point(2, 2))
#> [1] TRUE
midpoint(A, B)
#> Point2D(3/2, 2) 
```
