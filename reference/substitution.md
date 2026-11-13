# Solve a system of equations by substitution

`substitution(system, symbols, result, known_symbols, exclude, all_symbols)`
solves the equations of `system` one after the other for the unknowns
`symbols`, substituting the values found so far, as SymPy's
`substitution()`. It returns a `FiniteSet` of `Tuple`s ordered like
`all_symbols` (by default `symbols`); unknowns left free stay symbolic,
and periodic solutions are `ImageSet`s.

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

  A list of equations or expressions.

- symbols:

  A list of the unknowns.

- result:

  A list of known partial solutions, each a list of
  `list(symbol, value)` pairs or a named list.

- known_symbols:

  Symbols whose values are known.

- exclude:

  Expressions that must not vanish, such as denominators.

- all_symbols:

  The symbols to report, in order.

## Value

A set of tuples, or a `ConditionSet`.

## Examples

``` r
x <- Symbol("x"); y <- Symbol("y")
substitution(list(x + y), list(x), list(list(list(y, 1))), list(y), list(), list(x, y))
#> {(-1, 1)}
substitution(list(x^2 - 1, y - x), list(x, y))
#> {(-1, -1), (1, 1)}
```
