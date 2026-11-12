# Simplify boolean expressions

Minimizes a boolean expression into the shortest equivalent `"cnf"` or
`"dnf"` form with the Quine-McCluskey algorithm (up to `max_vars`
variables), like SymPy's `simplify_logic()`.

## Usage

``` r
simplify_logic(e, form = NULL, max_vars = 10)
```

## Arguments

- e:

  A boolean expression.

- form:

  `"cnf"`, `"dnf"`, or `NULL` to pick the shorter.

- max_vars:

  Largest number of variables to minimize.

## Value

A boolean expression.

## Examples

``` r
a <- Symbol("a"); b <- Symbol("b")
simplify_logic(Or(And(a, b), And(a, Not(b))))
#> a
```
