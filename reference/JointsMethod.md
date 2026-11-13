# Equations of motion of a set of joints

`JointsMethod()` is SymPy's deprecated `JointsMethod`: it collects the
bodies, loads, coordinates, speeds and kinematic equations of the given
joints, and `form_eoms(method)` (reached with `$`) forms the equations
of motion with
[`KanesMethod()`](https://adamrsophie.github.io/symcompR/reference/KanesMethod.md)
(the default) or
[`LagrangesMethod()`](https://adamrsophie.github.io/symcompR/reference/LagrangesMethod.md).
New code should use
[`System()`](https://adamrsophie.github.io/symcompR/reference/System.md).

## Usage

``` r
JointsMethod(newtonion, ...)
```

## Arguments

- newtonion:

  The Newtonian frame, or a body whose frame is used.

- ...:

  Joints.

## Value

A `JointsMethod` object.

## Examples

``` r
P <- RigidBody("P")
C <- RigidBody("C")
J <- PinJoint("J", P, C)
method <- suppressWarnings(JointsMethod(P, J))
method$q
#> Matrix(1, 1)
#> [q_J(t)]
```
