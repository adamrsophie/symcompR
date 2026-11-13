# Collection of characteristic curves

A frozen record of the seven characteristic curve constructors a
musculotendon model uses, as SymPy's `CharacteristicCurveCollection`
dataclass. Fields are read with `$` and cannot be reassigned.

## Usage

``` r
CharacteristicCurveCollection(
  tendon_force_length,
  tendon_force_length_inverse,
  fiber_force_length_passive,
  fiber_force_length_passive_inverse,
  fiber_force_length_active,
  fiber_force_velocity,
  fiber_force_velocity_inverse
)
```

## Arguments

- tendon_force_length, tendon_force_length_inverse,
  fiber_force_length_passive, fiber_force_length_passive_inverse,
  fiber_force_length_active, fiber_force_velocity,
  fiber_force_velocity_inverse:

  Curve constructors, such as
  [`TendonForceLengthDeGroote2016()`](https://adamrsophie.github.io/symcompR/reference/characteristic-curves.md).

## Value

A `CharacteristicCurveCollection` object.

## Examples

``` r
curves <- CharacteristicCurveCollection(
  TendonForceLengthDeGroote2016, TendonForceLengthInverseDeGroote2016,
  FiberForceLengthPassiveDeGroote2016, FiberForceLengthPassiveInverseDeGroote2016,
  FiberForceLengthActiveDeGroote2016, FiberForceVelocityDeGroote2016,
  FiberForceVelocityInverseDeGroote2016)
curves$tendon_force_length(Symbol("l"), 1, 2, 3, 4)
#> TendonForceLengthDeGroote2016(l, 1, 2, 3, 4)
```
