# Solve a system by substitution

`substitution(system, symbols, result, known_symbols, exclude, all_symbols)`
solves `system` for `symbols` after substituting each partial solution
in `result`, and returns a `FiniteSet` of `Tuple`s ordered like
`all_symbols` (default `symbols`). Unknowns left free stay symbolic.
This is SymPy's `substitution()` interface over
[`solve()`](https://rdrr.io/r/base/solve.html).

## Usage

``` r
substitution(
  system,
  symbols,
  result = list(list()),
  known_symbols = list(),
  exclude = list(),
  all_symbols = NULL
)
```

## Arguments

- system:

  A list of equations.

- symbols:

  Unknowns.

- result:

  A list of named lists of known values.

- known_symbols, exclude:

  Ignored; kept for SymPy compatibility.

- all_symbols:

  Symbols to report, in order.

## Value

A `FiniteSet`.

## Examples

``` r
x <- Symbol("x"); y <- Symbol("y")
substitution(list(x + y), list(x), list(list(y = 1)), list(y), list(), list(x, y))
#> {(-1, 1)}
```
