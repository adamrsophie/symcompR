# Loads: forces and torques

`Force()` is a force applied at a point and `Torque()` a torque applied
to a frame, as SymPy's named tuples of the same names. A body may stand
for its mass center (forces) or its frame (torques). `gravity()` returns
the gravitational forces on the given bodies.

## Usage

``` r
Force(point, force)

Torque(frame, torque)

gravity(acceleration, ...)
```

## Arguments

- point:

  A
  [`PhysPoint()`](https://adamrsophie.github.io/symcompR/reference/PhysPoint.md)
  or a body.

- force, torque:

  A vector.

- frame:

  A
  [`ReferenceFrame()`](https://adamrsophie.github.io/symcompR/reference/ReferenceFrame.md)
  or a rigid body.

- acceleration:

  The gravitational acceleration vector.

- ...:

  Bodies.

## Value

A `Force` or `Torque`; a list of forces for `gravity()`.

## Examples

``` r
N <- ReferenceFrame("N")
P <- PhysPoint("P")
Force(P, 3 * N$x)
#> Force(point=P, force=3*N.x)
Torque(N, 2 * N$z)
#> Torque(frame=N, torque=2*N.z)
gravity(-9.81 * N$y, Particle("p", P, 1))
#> [[1]]
#> Force(point=P, force=- 9.81*N.y)
#> 
```
