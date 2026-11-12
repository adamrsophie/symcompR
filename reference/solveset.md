# Solve equations as sets

`solveset(f, x, domain)` returns the solution set of `f = 0` (or of an
inequality) in `domain` (`Complexes` by default), like SymPy's
`solveset()`. Periodic trigonometric equations over the reals give
`ImageSet` unions; unsolvable equations give a `ConditionSet`.

## Usage

``` r
solveset(f, x = NULL, domain = Complexes)
```

## Arguments

- f:

  An equation, expression or inequality.

- x:

  The variable.

- domain:

  `Complexes`, `Reals` or another set.

## Value

A set.

## Examples

``` r
x <- Symbol("x")
solveset(x^2 - 1, x, Reals)
#> {-1, 1}
solveset(x^2 > 1, x, Reals)
#> Union(Interval.open(-oo, -1), Interval.open(1, oo))
solveset(sin(x), x, Reals)
#> Union(ImageSet(Lambda(_n, 2*_n*pi), Integers), ImageSet(Lambda(_n, 2*_n*pi + pi), Integers))
```
