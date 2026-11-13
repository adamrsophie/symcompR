# Actuators

The actuators of SymPy's `sympy.physics.mechanics`: `ForceActuator()`
applies a force along a pathway, `LinearSpring()`, `LinearDamper()`,
`DuffingSpring()` and `CoulombKineticFriction()` compute that force from
the pathway's length or extension velocity, and `TorqueActuator()`
applies a torque about an axis between two frames.
`TorqueActuator_at_pin_joint()` is SymPy's
`TorqueActuator.at_pin_joint()`. Each exposes `to_loads()` through `$`.

## Usage

``` r
ForceActuator(force, pathway)

LinearSpring(stiffness, pathway, equilibrium_length = 0)

LinearDamper(damping, pathway)

TorqueActuator(torque, axis, target_frame, reaction_frame = NULL)

TorqueActuator_at_pin_joint(torque, pin_joint)

DuffingSpring(
  linear_stiffness,
  nonlinear_stiffness,
  pathway,
  equilibrium_length = 0
)

CoulombKineticFriction(
  mu_k,
  f_n,
  pathway,
  v_s = NULL,
  sigma = NULL,
  mu_s = NULL
)
```

## Arguments

- force:

  The force magnitude.

- pathway:

  A pathway, see
  [`LinearPathway()`](https://adamrsophie.github.io/symcompR/reference/LinearPathway.md).

- stiffness, damping:

  Spring stiffness and damping coefficient.

- equilibrium_length:

  Spring length at rest.

- torque:

  The torque magnitude.

- axis:

  The torque axis.

- target_frame, reaction_frame:

  Frames (or rigid bodies) the torque and its reaction act on.

- pin_joint:

  A
  [`PinJoint()`](https://adamrsophie.github.io/symcompR/reference/Joint.md).

- linear_stiffness, nonlinear_stiffness:

  Coefficients of a Duffing spring.

- mu_k, mu_s:

  Kinetic and static friction coefficients.

- f_n:

  The normal force.

- v_s:

  The Stribeck velocity.

- sigma:

  The viscous friction coefficient.

## Value

An actuator object.

## Examples

``` r
N <- ReferenceFrame("N")
pA <- PhysPoint("pA"); pB <- PhysPoint("pB")
pB$set_pos(pA, Symbol("q") * N$x)
spring <- LinearSpring(Symbol("k"), LinearPathway(pA, pB))
spring$force
#> -k*sqrt(q^2)
spring$to_loads()
#> [[1]]
#> Force(point=pA, force=k*q*N.x)
#> 
#> [[2]]
#> Force(point=pB, force=- k*q*N.x)
#> 
```
