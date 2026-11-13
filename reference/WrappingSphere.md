# Wrapping geometries

`WrappingSphere()`, `WrappingCylinder()` and `WrappingCone()` are the
surfaces over which a
[`WrappingPathway()`](https://adamrsophie.github.io/symcompR/reference/LinearPathway.md)
wraps, as in SymPy. Each exposes `point`, `point_on_surface(point)`,
`geodesic_length(p1, p2)` and `geodesic_end_vectors(p1, p2)` through
`$`.

## Usage

``` r
WrappingSphere(radius, point)

WrappingCylinder(radius, point, axis)

WrappingCone(alpha, apex, axis)
```

## Arguments

- radius:

  The radius.

- point:

  The center of a sphere, or a point on the axis of a cylinder.

- axis:

  The axis vector.

- alpha:

  The half angle of a cone.

- apex:

  The apex of a cone.

## Value

A wrapping geometry.

## Examples

``` r
N <- ReferenceFrame("N")
O <- PhysPoint("O")
WrappingCylinder(Symbol("r"), O, N$x)
#> WrappingCylinder(radius=r, point=O, axis=N.x)
```
