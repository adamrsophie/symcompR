# Momentum, energy and the Lagrangian of a system of bodies

`linear_momentum()`, `angular_momentum()` (about `point`),
`kinetic_energy()` and `potential_energy()` sum the contributions of the
given particles and rigid bodies, as the SymPy functions of the same
names. `Lagrangian()` is the kinetic minus the potential energy and
`center_of_mass()` the position of the mass center from `point`.

## Usage

``` r
linear_momentum(frame, ...)

angular_momentum(point, frame, ...)

kinetic_energy(frame, ...)

potential_energy(...)

center_of_mass(point, ...)

Lagrangian(frame, ...)
```

## Arguments

- frame:

  The
  [`ReferenceFrame()`](https://adamrsophie.github.io/symcompR/reference/ReferenceFrame.md)
  in which the quantity is measured.

- ...:

  Particles and rigid bodies.

- point:

  A
  [`PhysPoint()`](https://adamrsophie.github.io/symcompR/reference/PhysPoint.md).

## Value

A `PhysVector` or a scalar expression.

## Examples

``` r
N <- ReferenceFrame("N")
P <- PhysPoint("P")
P$set_vel(N, 10 * N$x)
p <- Particle("p", P, 2)
linear_momentum(N, p)
#> 20*N.x
kinetic_energy(N, p)
#> 100
```
