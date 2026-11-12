# Mathematical constants and special values

`pi`, Euler's number `E`, the imaginary unit `I`, positive and negative
infinity `oo` and `-oo`, complex infinity `zoo`, not a number `nan`, and
the constants `EulerGamma`, `Catalan`, `GoldenRatio` and
`TribonacciConstant`. They are exported as objects of the same names
prefixed with `S_` (for example `S_Pi`) and, for convenience, as
`sym_pi`, `sym_E`, `sym_I`, `oo`, `zoo` and `nan`.

## Usage

``` r
sym_pi

sym_E

sym_I

oo

zoo

nan

EulerGamma

Catalan

GoldenRatio

TribonacciConstant

S_Pi

S_E

S_I

S_Infinity

S_NegativeInfinity

S_ComplexInfinity

S_NaN

S_Zero

S_One

S_Half

S_NegativeOne

S_true

S_false

S_EmptySet
```

## Format

`sym` objects.

## Examples

``` r
sym_pi + 1
#> 1 + pi
sym_I^2
#> -1
oo + 1
#> oo
```
