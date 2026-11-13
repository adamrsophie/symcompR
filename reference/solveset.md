# Solve an equation or inequality as a set

`solveset(f, symbol, domain)` returns the set of values of `symbol` in
`domain` (`Complexes` by default) for which `f = 0` holds, or for which
a relational `f` is true, as SymPy's `solveset()`. Periodic solutions
are `ImageSet`s over the integers, and equations it cannot solve give a
`ConditionSet`. `solveset_real()` and `solveset_complex()` fix the
domain.

## Usage

``` r
solveset(f, symbol = NULL, domain = Complexes)

solveset_real(f, symbol)

solveset_complex(f, symbol)
```

## Arguments

- f:

  An expression, equation or inequality.

- symbol:

  The variable; may be omitted when `f` has one free symbol.

- domain:

  The domain, a set.

## Value

A set.

## Examples

``` r
x <- Symbol("x")
solveset(x^2 - 1, x)
#> {-1, 1}
solveset(exp(x) - 1, x)
#> ImageSet(Lambda(_n, 2*_n*I*pi), Integers)
solveset(cos(x), x, Reals)
#> Union(ImageSet(Lambda(_n, 2*_n*pi + pi/2), Integers), ImageSet(Lambda(_n, 2*_n*pi + 3*pi/2), Integers))
solveset(x^2 > 1, x, Reals)
#> Union(Interval.open(-oo, -1), Interval.open(1, oo))
```
