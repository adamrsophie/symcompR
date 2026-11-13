# Solve equations as sets over the reals or the complexes

`solveset_real(f, symbol)` is `solveset(f, symbol, Reals)` and
`solveset_complex(f, symbol)` is `solveset(f, symbol, Complexes)`, as
SymPy's functions of the same names.

## Usage

``` r
solveset_real(f, symbol)

solveset_complex(f, symbol)
```

## Arguments

- f:

  An equation, expression or inequality.

- symbol:

  The variable.

## Value

A set.

## Examples

``` r
x <- Symbol("x")
solveset_real(x^2 + 1, x)
#> EmptySet
solveset_complex(x^2 + 1, x)
#> {-I, I}
```
