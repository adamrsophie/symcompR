# Bodies of the deprecated joints framework

`Body()` is SymPy's deprecated `Body`: a rigid body, or a particle when
only a mass is given, that also stores the loads applied to it. Methods
`apply_force()`, `apply_torque()`, `clear_loads()`, `remove_load()`,
`masscenter_vel()`, `ang_vel_in()`, `dcm()`,
[`kinetic_energy()`](https://adamrsophie.github.io/symcompR/reference/linear_momentum.md)
and `parallel_axis()` and the attributes `loads` and `is_rigidbody` are
reached with `$`. Creating a body emits a deprecation warning, as in
SymPy; new code should use
[`RigidBody()`](https://adamrsophie.github.io/symcompR/reference/Particle.md)
or
[`Particle()`](https://adamrsophie.github.io/symcompR/reference/Particle.md).

## Usage

``` r
Body(
  name,
  masscenter = NULL,
  mass = NULL,
  frame = NULL,
  central_inertia = NULL
)
```

## Arguments

- name:

  Body name.

- masscenter:

  The mass center; created when `NULL`.

- mass:

  The mass; a symbol when `NULL`.

- frame:

  The body frame; created when `NULL`.

- central_inertia:

  The central inertia dyadic.

## Value

A `Body`.

## Examples

``` r
B <- suppressWarnings(Body("B"))
B$apply_force(B$x)
B$loads
#> [[1]]
#> (B_masscenter, B_frame.x) 
#> 
```
