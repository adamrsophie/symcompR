# Physical dimensions

`Dimension(name, symbol)` is SymPy's `Dimension`: a name such as
`"length"` (or any expression) and an optional symbol such as `"L"`.
Dimensions multiply, divide and take powers to form new dimensions, as
in `sympy.physics.units`; adding equal dimensions gives the dimension
again and adding a quantity is an error. Equality ignores the symbol.
`dimension_name()` and `dimension_symbol()` return the parts, and
`has_integer_powers()` tells whether every exponent of the dependencies
in a dimension system is an integer.

## Usage

``` r
Dimension(name, symbol = NULL)

dimension_name(d)

dimension_symbol(d)

has_integer_powers(d, dim_sys)
```

## Arguments

- name:

  A string or expression.

- symbol:

  `NULL`, a string or a symbol.

- d:

  A `Dimension`.

- dim_sys:

  A `DimensionSystem`.

## Value

A `sym` Dimension, an expression, or `TRUE`/`FALSE`.

## Examples

``` r
u <- physics_units()
L <- Dimension("length", "L")
L / Dimension("time")
#> Dimension(length/time)
has_integer_powers(sqrt(u$length), u$dimsys_SI)
#> [1] FALSE
```
