# Beams

`Beam(length, elastic_modulus, second_moment)` is a beam bent in a
plane. `$apply_load(value, start, order, end)` adds a load,
`$apply_support(loc, type)` a `"pin"`, `"roller"` or `"fixed"` support
and `$apply_rotation_hinge()` and `$apply_sliding_hinge()` hinges;
`$solve_for_reaction_loads()` finds the reactions. `$shear_force()`,
`$bending_moment()`, `$slope()`, `$deflection()` and `$shear_stress()`
give the curves along the beam, `$max_shear_force()`, `$max_bmoment()`,
`$max_deflection()` and `$point_cflexure()` their extremes, and
`$solve_for_ild_reactions()`, `$solve_for_ild_shear()` and
`$solve_for_ild_moment()` the influence lines. The boundary conditions
are `$bc_deflection`, `$bc_slope`, `$bc_bending_moment` and
`$bc_shear_force`. `$join()` joins two beams. The `$plot_*()` methods
draw base R plots.
`Beam3D(length, elastic_modulus, shear_modulus, second_moment, area)` is
a beam loaded in three dimensions. After
`sympy.physics.continuum_mechanics.beam`.

## Usage

``` r
Beam(
  length,
  elastic_modulus,
  second_moment,
  area = Symbol("A"),
  variable = Symbol("x"),
  base_char = "C",
  ild_variable = Symbol("a")
)

Beam3D(
  length,
  elastic_modulus,
  shear_modulus,
  second_moment,
  area,
  variable = Symbol("x")
)
```

## Arguments

- length, elastic_modulus:

  The length and the elastic modulus.

- second_moment:

  The second moment of area, or a geometry entity giving the cross
  section.

- area:

  The cross section area.

- variable:

  The coordinate along the beam.

- base_char:

  The base name of the integration constants.

- ild_variable:

  The position of the moving load of the influence lines.

- shear_modulus:

  The shear modulus.

## Value

A `Beam` or `Beam3D` object.

## Examples

``` r
E <- Symbol("E"); I <- Symbol("I")
b <- Beam(30, E, I)
R1 <- b$apply_support(10, "roller"); R2 <- b$apply_support(30, "roller")
b$apply_load(-8, 0, -1)
b$solve_for_reaction_loads(R1, R2)
b$reaction_loads
#> {R_10: 12, R_30: -4} 
b$bending_moment()
#> 8*SingularityFunction(x, 0, 1) - 12*SingularityFunction(x, 10, 1) + 4*SingularityFunction(x, 30, 1)
```
