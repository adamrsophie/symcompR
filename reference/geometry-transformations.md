# Transformations of geometric entities

`geo_rotate()`, `geo_scale()`, `geo_reflect()` and `translate()` are
SymPy's `rotate`, `scale`, `reflect` and `translate` methods of points,
lines, polygons, ellipses and curves;
[`geo_transform()`](https://adamrsophie.github.io/symcompR/reference/Point.md)
applies a transformation matrix to a point. `geo_rotation_matrix()`,
`geo_scale_matrix()` and `geo_translation_matrix()` are the 3 by 3
matrices of `sympy.geometry.entity.rotate`, `scale` and `translate`.

## Usage

``` r
geo_rotation_matrix(th)

geo_translation_matrix(x, y)

geo_scale_matrix(x, y, pt = NULL)

geo_rotate(e, angle, pt = NULL)

geo_scale(e, x = 1, y = 1, pt = NULL, z = 1)

geo_reflect(e, line)

translate(e, x = 0, y = 0, z = 0)
```

## Arguments

- th:

  The angle of the rotation matrix.

- x, y, z:

  The scale factors or translation offsets.

- pt:

  The center of the rotation or scaling.

- e:

  An entity.

- angle:

  The rotation angle.

- line:

  The mirror line.

## Value

The transformed entity.

## Examples

``` r
geo_rotate(Point(1, 0), pi / 2)
#> Point2D(612323399573677/10000000000000000000000000000000, 1) 
geo_scale(Circle(Point(0, 0), 1), 2, 3)
#> Ellipse(Point2D(0, 0), 2, 3) 
geo_reflect(Point(1, 2), Line(Point(0, 0), Point(1, 1)))
#> Point2D(2, 1) 
translate(Point(1, 2), 3, 4)
#> Point2D(4, 6) 
geo_transform(Point(1, 0), geo_rotation_matrix(pi / 2))
#> Point2D(612323399573677/10000000000000000000000000000000, 1) 
```
