# Musculotendon formulations

The integer enumeration SymPy calls `MusculotendonFormulation`: `0` is
`RIGID_TENDON`, `1` `FIBER_LENGTH_EXPLICIT`, `2`
`TENDON_FORCE_EXPLICIT`, `3` `FIBER_LENGTH_IMPLICIT` and `4`
`TENDON_FORCE_IMPLICIT`. A member is selected by its value or by its
name.

## Usage

``` r
MusculotendonFormulation(value)

# S3 method for class 'MusculotendonFormulation'
format(x, ...)

# S3 method for class 'MusculotendonFormulation'
print(x, ...)
```

## Arguments

- value:

  An integer from 0 to 4, or a member name.

- x:

  A `MusculotendonFormulation` value.

- ...:

  Ignored.

## Value

A `MusculotendonFormulation` value, an integer that prints as its value.

## Examples

``` r
MusculotendonFormulation("RIGID_TENDON")
#> <MusculotendonFormulation.RIGID_TENDON: 0>
MusculotendonFormulation(2) == 2
#> [1] TRUE
```
