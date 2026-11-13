# Arches

`Arch(left_support, right_support, crown_x = , crown_y = )` is a three
hinged parabolic arch between two supports, each given as `list(x, y)`.
The crown abscissa is needed; without `crown_y` the parabola is fitted
through the supports.
`$apply_load(order, label, start, mag, end, angle)` adds a point load
(`order = -1`, with the direction `angle` in degrees) or a vertical load
distributed from `start` to `end` (`order = 0`), and `$remove_load()`
removes one. `$change_support_type()` makes a support a `"roller"` or a
`"hinge"`, `$add_member(y)` adds a tie at height `y`, and `$solve()`
finds the reaction forces, after which `$bending_moment_at()`,
`$shear_force_at()` and `$axial_force_at()` evaluate the internal
forces. The properties are `$get_shape_eqn`, `$get_loads`, `$supports`,
`$left_support`, `$right_support` and `$reaction_force`. `$draw()` draws
the arch with base R graphics. After
`sympy.physics.continuum_mechanics.arch`.

## Usage

``` r
Arch(left_support, right_support, ...)
```

## Arguments

- left_support, right_support:

  The supports, as `list(x, y)`.

- ...:

  `crown_x` and `crown_y`, the position of the crown hinge.

## Value

An `Arch` object.

## Examples

``` r
a <- Arch(list(0, 0), list(10, 0), crown_x = 5, crown_y = 5)
a$get_shape_eqn
#> 5 - (x - 5)^2/5
a$apply_load(-1, "C", 3, 10, angle = 270)
a$solve()
a$reaction_force
#> {R_A_x: 3, R_A_y: 7, R_B_x: -3, R_B_y: 3} 
```
