# Interval arithmetic for plotting

SymPy's `sympy.plotting.intervalmath`. `interval(a, b)` is an interval
of doubles with a validity flag: `TRUE` when a function result is
continuous and in the domain, `FALSE` when the argument is outside the
domain, and `NA` (SymPy's `None`) when it is partly outside or the
function is not continuous there. Arithmetic, powers and the elementary
functions (`exp`, `log`, `log10`, `sqrt`, `abs`, `sin`, `cos`, `tan`,
`asin`, `acos`, `atan`, `sinh`, `cosh`, `tanh`, `asinh`, `acosh`,
`atanh`, `floor`, `ceiling`) work on intervals; `imin()`, `imax()` and
`ceil()` are SymPy's names. Comparisons return an `intervalMembership`
pair (truth, validity), combined with `&`, `|`,
[`xor()`](https://rdrr.io/r/base/Logic.html) and `!`.

## Usage

``` r
interval(a, b = NULL, is_valid = TRUE)

imin(...)

imax(...)

ceil(x)

intervalMembership(a, b)
```

## Arguments

- a, b:

  The bounds; one value gives a point interval.

- is_valid:

  The validity flag.

- ...:

  Intervals or numbers.

- x:

  An interval or a number.

## Value

An interval, or an `intervalMembership`.

## Examples

``` r
i <- interval(1, 2)
i + 3
#> interval(4.000000, 5.000000) 
i * interval(-1, 2)
#> interval(-2.000000, 4.000000) 
sqrt(interval(-1, 4))
#> interval(-inf, inf) 
sin(interval(0, 3))
#> interval(0.000000, 1.000000) 
imax(interval(0, 1), interval(-1, 3))
#> interval(0.000000, 3.000000) 
i > 0
#> intervalMembership(True, True) 
i$mid
#> [1] 1.5
```
