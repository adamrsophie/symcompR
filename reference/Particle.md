# Particles and rigid bodies

`Particle()` and `RigidBody()` create the bodies of SymPy's
`sympy.physics.mechanics`. A particle has a mass and a point; a rigid
body adds a reference frame and an inertia. Attributes and methods
(`mass`, `masscenter`, `frame`, `inertia`, `central_inertia`,
`potential_energy`,
[`linear_momentum()`](https://adamrsophie.github.io/symcompR/reference/linear_momentum.md),
[`angular_momentum()`](https://adamrsophie.github.io/symcompR/reference/linear_momentum.md),
[`kinetic_energy()`](https://adamrsophie.github.io/symcompR/reference/linear_momentum.md)
and `parallel_axis()`) are reached with `$`.

## Usage

``` r
Particle(name, point = NULL, mass = NULL)

RigidBody(name, masscenter = NULL, frame = NULL, mass = NULL, inertia = NULL)
```

## Arguments

- name:

  Body name.

- point, masscenter:

  The mass center, a
  [`PhysPoint()`](https://adamrsophie.github.io/symcompR/reference/PhysPoint.md);
  created when `NULL`.

- mass:

  The mass; a symbol `<name>_mass` when `NULL`.

- frame:

  The body frame; a new
  [`ReferenceFrame()`](https://adamrsophie.github.io/symcompR/reference/ReferenceFrame.md)
  when `NULL`.

- inertia:

  An
  [`Inertia()`](https://adamrsophie.github.io/symcompR/reference/inertia.md)
  tuple or a list of a dyadic and a point; built from symbolic moments
  of inertia when `NULL`.

## Value

A `Particle` or `RigidBody`.

## Examples

``` r
N <- ReferenceFrame("N")
P <- PhysPoint("P")
P$set_vel(N, 3 * N$x)
p <- Particle("p", P, 2)
p$linear_momentum(N)
#> 6*N.x
p$kinetic_energy(N)
#> 9
B <- RigidBody("B")
B$central_inertia
#> B_ixx*(B_frame.x|B_frame.x) + B_ixy*(B_frame.x|B_frame.y) + B_izx*(B_frame.x|B_frame.z) + B_ixy*(B_frame.y|B_frame.x) + B_iyy*(B_frame.y|B_frame.y) + B_iyz*(B_frame.y|B_frame.z) + B_izx*(B_frame.z|B_frame.x) + B_iyz*(B_frame.z|B_frame.y) + B_izz*(B_frame.z|B_frame.z)
```
