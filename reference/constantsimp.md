# Simplify arbitrary constants in an ODE solution

`constantsimp(expr, constants)` absorbs subexpressions made only of the
arbitrary constants into those constants and removes redundant ones, as
SymPy's `constantsimp()` does for the solutions returned by
[`dsolve()`](https://adamrsophie.github.io/symcompR/reference/dsolve.md).
`constant_renumber(expr)` renumbers constants named `C1`, `C2`, ... in
printing order.

## Usage

``` r
constantsimp(expr, constants)

constant_renumber(expr, variables = NULL, newconstants = NULL)
```

## Arguments

- expr:

  An expression or equation.

- constants:

  A list of the constant symbols.

- variables:

  Symbols that are not constants; by default the symbols named `C`
  followed by digits are the constants.

- newconstants:

  Symbols to use as the new constants.

## Value

An expression or equation.

## Examples

``` r
x <- Symbol("x"); C1 <- Symbol("C1"); C2 <- Symbol("C2")
constantsimp(2 * C1 * x, list(C1))
#> C1*x
constantsimp(C1 + C2 + x, list(C1, C2))
#> C1 + x
constant_renumber(C2 * x + C1)
#> C1 + C2*x
```
