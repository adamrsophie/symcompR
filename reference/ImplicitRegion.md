# Implicit regions

`ImplicitRegion(variables, equation)` is the region where `equation` (or
an [`Eq()`](https://adamrsophie.github.io/symcompR/reference/Eq.md))
vanishes. `region$variables`, `region$equation` and
`implicit_degree(region)` read it; `regular_point()`,
`singular_points()`,
[`multiplicity()`](https://adamrsophie.github.io/symcompR/reference/multiplicity.md)
(via the region method below) and `rational_parametrization()` port the
methods of SymPy's `ImplicitRegion`, and `conic_coeff()` gives the
coefficients of a conic.

## Usage

``` r
ImplicitRegion(variables, equation)

implicit_degree(region)

conic_coeff(variables, equation)

regular_point(region)

singular_points(region)

rational_parametrization(region, parameters = c("t", "s"), reg_point = NULL)
```

## Arguments

- variables:

  A list of symbols.

- equation:

  An expression or equation.

- region:

  An implicit region.

- parameters:

  Parameter names.

- reg_point:

  A regular point to use for a conic.

## Value

An `ImplicitRegion`, expression, list or set.

## Examples

``` r
x <- Symbol("x"); y <- Symbol("y")
r <- ImplicitRegion(list(x, y), x^2 + y^2 - 4)
regular_point(r)
#> (0, -2) 
rational_parametrization(r)
#> (4*t/(t^2 + 1), 4*t^2/(t^2 + 1) - 2) 
```
