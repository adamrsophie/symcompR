# Columns under axial loads

`Column(length, elastic_modulus, area)` is a column loaded along its
axis. `$apply_support(loc)` adds a support,
`$apply_load(value, start, order, end)` a load in singularity function
form and `$apply_telescope_hinge(loc)` a telescope hinge;
`$solve_for_reaction_loads()` finds the reactions, after which
`$axial_force()` and `$extension()` give the axial force and the
extension along the column. The properties are `$length`,
`$elastic_modulus`, `$area`, `$variable`, `$load`, `$applied_loads`,
`$reaction_loads` and `$hinge_extensions`. `$plot_axial_force()` and
`$plot_extension()` draw base R plots. After
`sympy.physics.continuum_mechanics.column`.

## Usage

``` r
Column(length, elastic_modulus, area, variable = Symbol("x"), base_char = "C")
```

## Arguments

- length, elastic_modulus, area:

  The length, the elastic modulus and the cross section area.

- variable:

  The coordinate along the column.

- base_char:

  The base name of the integration constants.

## Value

A `Column` object.

## Examples

``` r
E <- Symbol("E"); A <- Symbol("A")
c <- Column(10, E, A)
c$apply_support(0)
c$apply_load(-1, 5, -1)
c$solve_for_reaction_loads()
c$reaction_loads
#> {R_0: 1} 
```
