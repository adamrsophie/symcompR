# Physical quantities

`Quantity(name, abbrev)` is SymPy's physical quantity: a unit of
measure, a constant or a generic quantity. It is a positive real atom
that prints as its name and takes part in arithmetic like a symbol
(`Unit` is an alias, as in `sympy.physics.units`). `PhysicalConstant()`
marks a physical constant such as `speed_of_light`. A quantity has no
scale factor or dimension of its own: set them globally with
`set_global_relative_scale_factor()` and `set_global_dimension()`, or
per system with
[`set_quantity_scale_factor()`](https://adamrsophie.github.io/symcompR/reference/set_quantity_dimension.md)
and
[`set_quantity_dimension()`](https://adamrsophie.github.io/symcompR/reference/set_quantity_dimension.md).
`quantity_dimension()` and
[`scale_factor()`](https://adamrsophie.github.io/symcompR/reference/Prefix.md)
read them in the SI system. The predefined units and constants are in
[`physics_units()`](https://adamrsophie.github.io/symcompR/reference/physics_units.md).

## Usage

``` r
Quantity(
  name,
  abbrev = NULL,
  latex_repr = NULL,
  pretty_unicode_repr = NULL,
  pretty_ascii_repr = NULL,
  mathml_presentation_repr = NULL,
  is_prefixed = FALSE,
  ...
)

PhysicalConstant(
  name,
  abbrev = NULL,
  latex_repr = NULL,
  pretty_unicode_repr = NULL,
  pretty_ascii_repr = NULL,
  mathml_presentation_repr = NULL,
  is_prefixed = FALSE,
  ...
)

Unit(
  name,
  abbrev = NULL,
  latex_repr = NULL,
  pretty_unicode_repr = NULL,
  pretty_ascii_repr = NULL,
  mathml_presentation_repr = NULL,
  is_prefixed = FALSE,
  ...
)

set_global_dimension(q, dimension)

set_global_relative_scale_factor(q, scale_factor, reference_quantity)

quantity_name(q)

quantity_abbrev(q)

quantity_dimension(q)

is_prefixed(q)

is_physical_constant(q)
```

## Arguments

- name:

  The name, a string or symbol.

- abbrev:

  The abbreviation; the name by default.

- latex_repr, pretty_unicode_repr, pretty_ascii_repr,
  mathml_presentation_repr:

  Optional printed forms.

- is_prefixed:

  Whether the quantity carries a prefix.

- ...:

  Ignored assumptions.

- q:

  A quantity.

- dimension:

  A `Dimension`.

- scale_factor:

  A number or expression, possibly with prefixes.

- reference_quantity:

  The quantity the scale factor refers to.

## Value

A `sym` Quantity, an expression, `TRUE`/`FALSE`, or `NULL` invisibly for
the setters.

## Examples

``` r
u <- physics_units()
q <- Quantity("q1")
set_global_relative_scale_factor(q, 5000, u$meter)
convert_to(q, u$meter)
#> 5000*meter
scale_factor(u$kilometer)
#> 1000
quantity_dimension(u$newton)
#> Dimension(force, F)
```
