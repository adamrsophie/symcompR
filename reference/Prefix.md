# Unit prefixes

`Prefix(name, abbrev, exponent, base)` is SymPy's prefix with scale
factor `base^exponent`, such as `kilo` (`10^3`) or `kibi` (`2^10`). A
prefix times a quantity gives the scaled quantity, and products and
quotients of prefixes simplify to a prefix or a number, as in
`sympy.physics.units`. The predefined prefixes are in
[`physics_units()`](https://adamrsophie.github.io/symcompR/reference/physics_units.md),
with the dictionaries `PREFIXES` and `BIN_PREFIXES`.
`prefix_unit(unit, prefixes)` returns the prefixed units of `unit` for a
named list of prefixes. `scale_factor()` gives the scale factor of a
prefix or of a quantity in the SI system.

## Usage

``` r
Prefix(name, abbrev, exponent, base = 10, latex_repr = NULL)

scale_factor(x)

prefix_unit(unit, prefixes)
```

## Arguments

- name:

  The name.

- abbrev:

  The abbreviation.

- exponent:

  The exponent of the base.

- base:

  The base, 10 by default.

- latex_repr:

  Optional LaTeX form.

- x:

  A prefix or quantity.

- unit:

  A quantity.

- prefixes:

  A named list of prefixes.

## Value

A `sym` Prefix, a list of quantities, or an expression.

## Examples

``` r
u <- physics_units()
u$kilo * u$meter
#> 1000*meter
u$kilo * u$kilo
#> M
prefix_unit(u$meter, u$PREFIXES[c("m", "c")])
#> [[1]]
#> millimeter
#> 
#> [[2]]
#> centimeter
#> 
```
