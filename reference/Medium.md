# An optical medium

`Medium(name, permittivity, permeability, n)` is a medium given by any
two of its permittivity, permeability and refractive index `n`, or the
vacuum when none is given. Its properties are `$speed`,
`$refractive_index` and `$n`, and for a medium given by its permittivity
and permeability also `$permittivity`, `$permeability` and
`$intrinsic_impedance`. After `sympy.physics.optics.medium`.

## Usage

``` r
Medium(name, permittivity = NULL, permeability = NULL, n = NULL)
```

## Arguments

- name:

  The name of the medium.

- permittivity, permeability:

  The electric permittivity and magnetic permeability, with units.

- n:

  The refractive index.

## Value

A `Medium` object (of class `MediumPP` or `MediumN`).

## Examples

``` r
m <- Medium("m", n = 2)
m$speed
#> 149896229*meter/second
Medium("vacuum")$refractive_index
#> 1
```
