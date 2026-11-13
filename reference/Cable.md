# Cables

`Cable(support_1, support_2)` is a cable hung between two supports, each
given as `list(label, x, y)`.
`$apply_load(-1, list(label, x, y, magnitude, direction))` adds a point
load with the direction in degrees and
`$apply_load(0, list(label, magnitude))` a load distributed along the
horizontal; `$remove_loads()` removes loads by label. `$solve()` finds
the tension and the reaction loads; with a distributed load it needs the
abscissa of the lowest point of the cable. `$tension_at(x)` evaluates
the tension of a cable under a distributed load. The properties are
`$supports`, `$left_support`, `$right_support`, `$loads`,
`$loads_position`, `$length`, `$reaction_loads` and `$tension`.
`$change_support()` and `$apply_length()` change the supports and the
length, and `$draw()` and `$plot_tension()` draw base R plots. After
`sympy.physics.continuum_mechanics.cable`.

## Usage

``` r
Cable(support_1, support_2)
```

## Arguments

- support_1, support_2:

  The supports, as `list(label, x, y)`.

## Value

A `Cable` object.

## Examples

``` r
c <- Cable(list("A", 0, 10), list("B", 10, 10))
c$apply_load(-1, list("X", 5, 5, 10, 270))
c$solve()
c$reaction_loads
#> {R_A_x: -5, R_A_y: 5, R_B_x: 5, R_B_y: 5} 
```
