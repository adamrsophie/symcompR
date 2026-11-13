# Multibody systems

`System()` is SymPy's `System` class: a container of bodies, joints,
loads, actuators, generalized coordinates and speeds, kinematic
differential equations and constraints, which forms the equations of
motion with `form_eoms()` using
[`KanesMethod()`](https://adamrsophie.github.io/symcompR/reference/KanesMethod.md)
(the default) or
[`LagrangesMethod()`](https://adamrsophie.github.io/symcompR/reference/LagrangesMethod.md).
`System_from_newtonian()` is SymPy's `System.from_newtonian()`.
Attributes, setters (`$<-`) and the methods `add_coordinates()`,
`add_speeds()`, `add_auxiliary_speeds()`, `add_kdes()`,
`add_holonomic_constraints()`, `add_nonholonomic_constraints()`,
`add_bodies()`, `add_loads()`, `apply_uniform_gravity()`,
`add_actuators()`, `add_joints()`, `get_body()`, `get_joint()`,
`form_eoms()`,
[`rhs()`](https://adamrsophie.github.io/symcompR/reference/Eq.md) and
`validate_system()` are reached with `$`. Methods that take a variable
number of objects in SymPy take them through `...`; `independent` is a
named argument.

## Usage

``` r
System(frame = NULL, fixed_point = NULL)

System_from_newtonian(newtonian)
```

## Arguments

- frame:

  The inertial
  [`ReferenceFrame()`](https://adamrsophie.github.io/symcompR/reference/ReferenceFrame.md);
  `inertial_frame` when `NULL`.

- fixed_point:

  The fixed
  [`PhysPoint()`](https://adamrsophie.github.io/symcompR/reference/PhysPoint.md);
  `inertial_point` when `NULL`.

- newtonian:

  A rigid body whose frame and mass center are used.

## Value

A `System` object.

## Examples

``` r
sys <- System()
q <- dynamicsymbols("q"); u <- dynamicsymbols("u")
sys$add_coordinates(q)
sys$add_speeds(u)
sys$add_kdes(u - dynamicsymbols("q", 1))
sys$q
#> Matrix(1, 1)
#> [q(t)]
sys$kdes
#> Matrix(1, 1)
#> [u(t) - Derivative(q(t), t)]
```
