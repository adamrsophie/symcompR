# Unit systems

`UnitSystem(base_units, units, name, descr, dimension_system, derived_units)`
is SymPy's coherent set of units on top of a
[`DimensionSystem()`](https://adamrsophie.github.io/symcompR/reference/DimensionSystem.md).
Creating one registers it under its name for `get_unit_system()`; the
predefined systems are `"SI"`, `"MKS"`, `"MKSA"`, `"Natural system"` and
`"cgs_gauss"`. `get_dimensional_expr()` gives the dimension of an
expression as an expression in the dimension names, and
`collect_factor_and_dimension()` (SymPy's
`_collect_factor_and_dimension`) returns
`list(scale_factor, dimension)`, checking that added terms and function
arguments agree in dimension. `get_units_non_prefixed()` lists the units
without a prefix that are not physical constants.

## Usage

``` r
UnitSystem(
  base_units,
  units = list(),
  name = "",
  descr = "",
  dimension_system = NULL,
  derived_units = list()
)

# S3 method for class 'UnitSystem'
extend(
  x,
  base,
  units = list(),
  name = "",
  description = "",
  dimension_system = NULL,
  derived_units = list(),
  ...
)

get_dimension_system(us)

get_unit_system(unit_system)

get_default_unit_system()

derived_units(us)

get_dimensional_expr(us, expr)

collect_factor_and_dimension(us, expr)

get_units_non_prefixed(us)
```

## Arguments

- base_units:

  A list of base quantities.

- units, base:

  A list of further quantities.

- name:

  The name of the system.

- descr, description:

  A description.

- dimension_system:

  A `DimensionSystem`.

- derived_units:

  Preferred units of dimensions, as a list of `list(dimension, unit)`
  pairs or a `units_dict`.

- ...:

  Unused.

- unit_system, us, x:

  A `UnitSystem` or its name.

- expr:

  An expression.

## Value

A `UnitSystem`, an expression, a list or a `units_dict`.

## Examples

``` r
u <- physics_units()
SI <- get_unit_system("SI")
get_dimensional_expr(SI, u$newton / u$meter)
#> force/length
collect_factor_and_dimension(SI, u$meter + u$kilometer)
#> [[1]]
#> 1001
#> 
#> [[2]]
#> Dimension(length, L)
#> 
```
