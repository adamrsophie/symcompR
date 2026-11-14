# Operators on geometric entities

`==` and `!=` compare entities as SymPy does: points, lines and circles
structurally, polygons up to rotation and reversal of their vertices.
Points also support `+`, `-`, `*`, `/` and unary minus.

## Usage

``` r
# S3 method for class 'GeoEntity'
Ops(e1, e2)
```

## Arguments

- e1, e2:

  Entities, or a point and a scalar or a coordinate list.

## Value

An entity or a logical value.

## Examples

``` r
Point(1, 2) + Point(3, 4)
#> Point2D(4, 6) 
Point(1, 2) * 3
#> Point2D(3, 6) 
Line(Point(0, 0), Point(1, 1)) == Line(Point(0, 0), Point(2, 2))
#> [1] FALSE
```
