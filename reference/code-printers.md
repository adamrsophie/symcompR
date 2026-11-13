# Code generation

Print an expression as C (`ccode`), Fortran (`fcode`), JavaScript
(`jscode`), Octave/MATLAB (`octave_code`), R (`rcode`) or Python
(`pycode`) source.

## Usage

``` r
ccode(e, assign_to = NULL, ...)

fcode(e, assign_to = NULL, ...)

jscode(e, assign_to = NULL, ...)

octave_code(e, assign_to = NULL, ...)

rcode(e, assign_to = NULL, ...)

pycode(e, assign_to = NULL, ...)
```

## Arguments

- e:

  A `sym` expression.

- assign_to:

  Optional name of a variable to assign.

- ...:

  Printer settings, as in SymPy (for example `standard`, `precision`,
  `human`, `contract` or `user_functions`).

## Value

A character scalar.

## Examples

``` r
x <- Symbol("x")
ccode(x^2 + sin(x) / 3)
#> [1] "pow(x, 2) + (1.0/3.0)*sin(x)"
rcode(sqrt(x + 1), assign_to = "y")
#> [1] "y = sqrt(x + 1);"
```
