# Partial fraction decomposition

`apart(e, x)` decomposes a rational function of `x` into partial
fractions over the rationals, using undetermined coefficients as SymPy
does. With `full = TRUE` the decomposition runs over the roots of the
denominator factors (Bronstein's full decomposition), leaving `RootSum`
terms for irreducible factors. `apart_list()` returns that decomposition
in SymPy's structured form: the common factor, the polynomial part as a
[`Poly()`](https://adamrsophie.github.io/symcompR/reference/Poly.md),
and for each term the root polynomial in a dummy `w`, the numerator and
denominator as functions of a root, and the exponent.
`assemble_partfrac_list()` rebuilds the expression from that list.

## Usage

``` r
apart(e, x = NULL, full = FALSE)

apart_list(e, x = NULL)

assemble_partfrac_list(partial_list)
```

## Arguments

- e:

  A rational function.

- x:

  The variable.

- full:

  Full decomposition over the roots of the denominator.

- partial_list:

  A result of `apart_list()`.

## Value

An expression.

## Examples

``` r
x <- Symbol("x")
apart(1 / (x^2 - 1), x)
#> -1/(2*(x + 1)) + 1/(2*(x - 1))
apart((x^2 + 1) / (x^3 - x), x)
#> 1/(x + 1) + 1/(x - 1) - 1/x
```
