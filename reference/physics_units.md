# Predefined physical units

Returns the environment that holds the dimensions, prefixes, units,
physical constants and unit systems of `sympy.physics.units`, under
their SymPy names and aliases (`meter`, `m`, `kilometer`,
`speed_of_light`, `c`, `length`, `kilo`, ...), together with the systems
`SI`, `MKS`, `MKSA`, `natural` and `cgs_gauss`, the dimension systems
`dimsys_SI`, `dimsys_MKSA`, `dimsys_length_weight_time`, `dimsys_cgs`
and `dimsys_default`, and the prefix lists `PREFIXES` and
`BIN_PREFIXES`. Since several names (`length`, `time`, `c`, `t`, `T`,
`F`, `S`, `N`) would mask R or package functions, they are not exported
one by one; use `u <- physics_units()` and `u$meter`, or
`with(physics_units(), ...)`.

## Usage

``` r
physics_units()
```

## Value

An environment.

## Examples

``` r
u <- physics_units()
5 * u$kilometer
#> 5*kilometer
with(physics_units(), convert_to(mile, kilometer))
#> 25146*kilometer/15625
```
