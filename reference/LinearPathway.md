# Pathways of actuators

`LinearPathway()` is the straight line between two points,
`ObstacleSetPathway()` the piecewise straight path through three or more
points and `WrappingPathway()` the geodesic between two points over a
wrapping geometry, as in SymPy. Each exposes `attachments`, `length`,
`extension_velocity` and `to_loads(force)` through `$`.

## Usage

``` r
LinearPathway(...)

ObstacleSetPathway(...)

WrappingPathway(attachment_1, attachment_2, geometry)
```

## Arguments

- ...:

  The attachment points.

- attachment_1, attachment_2:

  The end points of a wrapping pathway.

- geometry:

  A wrapping geometry, see
  [`WrappingSphere()`](https://adamrsophie.github.io/symcompR/reference/WrappingSphere.md).

## Value

A pathway object.

## Examples

``` r
N <- ReferenceFrame("N")
pA <- PhysPoint("pA"); pB <- PhysPoint("pB")
pB$set_pos(pA, Symbol("q") * N$x)
path <- LinearPathway(pA, pB)
path$length
#> sqrt(q^2)
```
