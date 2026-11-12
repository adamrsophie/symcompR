# Dimensions and scale factors of quantities in a system

The `_QuantityMapper` methods shared by
[`DimensionSystem()`](https://adamrsophie.github.io/symcompR/reference/DimensionSystem.md)
and
[`UnitSystem()`](https://adamrsophie.github.io/symcompR/reference/UnitSystem.md).
`set_quantity_dimension()` and `set_quantity_scale_factor()` record a
quantity's dimension and scale factor in one system (prefixes and
quantities inside the scale factor are replaced by their scale factors
at once); `get_quantity_dimension()` and `get_quantity_scale_factor()`
look them up, falling back to the global definitions made with
[`set_global_dimension()`](https://adamrsophie.github.io/symcompR/reference/Quantity.md)
and
[`set_global_relative_scale_factor()`](https://adamrsophie.github.io/symcompR/reference/Quantity.md).

## Usage

``` r
set_quantity_dimension(system, quantity, dimension)

set_quantity_scale_factor(system, quantity, scale_factor)

get_quantity_dimension(system, unit)

get_quantity_scale_factor(system, unit)
```

## Arguments

- system:

  A `DimensionSystem` or `UnitSystem`.

- quantity, unit:

  A quantity (or, for the scale factor, a product or power of
  quantities).

- dimension:

  A `Dimension` or `1`.

- scale_factor:

  An expression.

## Value

`NULL` invisibly for the setters; a `Dimension` or expression.

## Examples

``` r
u <- physics_units()
get_quantity_dimension(u$SI, u$newton)
#> Dimension(force, F)
get_quantity_scale_factor(u$SI, u$kilometer * u$gram)
#> 1000
```
