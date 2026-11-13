# Example multibody models

`multi_mass_spring_damper()` builds the equations of motion of `n`
masses hanging in series from springs and dampers, and
`n_link_pendulum_on_cart()` those of an `n` link pendulum on a sliding
cart, as SymPy's `sympy.physics.mechanics.models`. Both return a
[`KanesMethod()`](https://adamrsophie.github.io/symcompR/reference/KanesMethod.md)
object whose `kanes_equations()` have been formed.

## Usage

``` r
multi_mass_spring_damper(
  n = 1,
  apply_gravity = FALSE,
  apply_external_forces = FALSE
)

n_link_pendulum_on_cart(n = 1, cart_force = TRUE, joint_torques = FALSE)
```

## Arguments

- n:

  Number of masses or links.

- apply_gravity:

  Add gravity to each mass.

- apply_external_forces:

  Add a specified force to each mass.

- cart_force:

  Apply a specified force to the cart.

- joint_torques:

  Apply specified torques at the joints.

## Value

A `KanesMethod` object.

## Examples

``` r
kane <- multi_mass_spring_damper(1)
kane$mass_matrix
#> Matrix(1, 1)
#> [m0]
kane$forcing
#> Matrix(1, 1)
#> [-c0*v0(t) - k0*x0(t)]
```
