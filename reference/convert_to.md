# Convert and simplify expressions with units

`convert_to(expr, target_units, unit_system)` rewrites the units and
quantities of `expr` as products of powers of `target_units` (one unit,
or a list of units) whenever the dimensions allow it, as SymPy's
`convert_to()`. `quantity_simplify()` replaces prefixes by their values
and units of the same dimension by one of them; with
`across_dimensions = TRUE` it also converts to the preferred unit of the
resulting dimension in `unit_system`. `check_dimensions()` returns
`expr` when the terms of every sum agree in dimension, and raises a
`ValueError` otherwise. `find_unit()` lists the names of the predefined
units that contain a string, or that share the dimension of a unit,
dimension or expression.

## Usage

``` r
convert_to(expr, target_units, unit_system = "SI")

quantity_simplify(expr, across_dimensions = FALSE, unit_system = NULL)

check_dimensions(expr, unit_system = "SI")

find_unit(quantity, unit_system = "SI")
```

## Arguments

- expr:

  An expression.

- target_units:

  A unit, a list of units, or a `Tuple`.

- unit_system:

  A `UnitSystem` or its name.

- across_dimensions:

  Combine units of different dimensions.

- quantity:

  A string, quantity, dimension or expression.

## Value

An expression, or a character vector for `find_unit()`.

## Examples

``` r
u <- physics_units()
convert_to(u$mile, u$kilometer)
#> 25146*kilometer/15625
convert_to(u$speed_of_light, list(u$meter, u$second))
#> 299792458*meter/second
quantity_simplify(u$kilo * u$foot * u$inch)
#> 250*foot^2/3
find_unit("volt")
#> [1] "volt"           "volts"          "electronvolt"   "electronvolts" 
#> [5] "planck_voltage"
```
