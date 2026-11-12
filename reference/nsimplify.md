# Find simple exact forms of numbers

`nsimplify(x)` replaces floating point numbers by nearby rationals,
square roots of rationals, or rational multiples of constants such as
`pi` and `E` (or of the constants listed in `constants`).

## Usage

``` r
nsimplify(e, constants = list(), tolerance = 1e-10, rational = FALSE)
```

## Arguments

- e:

  A number or expression.

- constants:

  A list of constants to try.

- tolerance:

  Relative tolerance.

- rational:

  If `TRUE`, only produce rationals.

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
```
