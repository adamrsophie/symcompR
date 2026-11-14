# Find simple exact forms of numbers

`nsimplify(e)` replaces a number, or the floating point numbers in an
expression, by a simple exact form, as SymPy's `nsimplify()`: the value
is evaluated to 30 digits and identified (with mpmath's
[`identify()`](https://rdrr.io/r/graphics/identify.html)) as a rational,
a quadratic irrational, or a simple formula in the given `constants`;
otherwise its Floats become rationals.

## Usage

``` r
nsimplify(
  e,
  constants = list(),
  tolerance = NULL,
  full = FALSE,
  rational = NULL,
  rational_conversion = "base10"
)
```

## Arguments

- e:

  A number or expression.

- constants:

  A list of constants (such as `sym_pi`) to try.

- tolerance:

  The tolerance of the identification, or `NULL` for `10^-d` with `d`
  the smallest precision (at most 15 digits) of the Floats in `e`; with
  `rational = TRUE`, the tolerance of the rational approximation.

- full:

  Try harder, returning the shortest formula found.

- rational:

  `TRUE` converts the Floats to rationals only; `FALSE` returns the
  identified value or `e` itself; `NULL` falls back to rationals when no
  formula is found.

- rational_conversion:

  `"base10"` converts a Float through its decimal digits, `"exact"`
  through its exact binary value.

## Value

An expression.

## Examples

``` r
nsimplify(0.333333333333)
#> 1/3
nsimplify(1.414213562373095)
#> sqrt(2)
nsimplify(6.283185307179586, constants = list(sym_pi))
#> 2*pi
nsimplify(sym_pi, tolerance = 0.01)
#> 22/7
```
