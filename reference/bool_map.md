# Map variables between equivalent boolean expressions

`bool_map(b1, b2)` looks for a one to one renaming of the variables of
`b1` that makes it equivalent to `b2`. It returns `list(expr, mapping)`
with the simplified `b1` and a named list from `b1` variables to `b2`
variables, or `FALSE`.

## Usage

``` r
bool_map(b1, b2)
```

## Arguments

- b1, b2:

  Boolean expressions.

## Value

A list or `FALSE`.

## Examples

``` r
a <- Symbol("a"); b <- Symbol("b"); x <- Symbol("x"); y <- Symbol("y")
bool_map(Or(And(a, b), Not(a)), Or(And(x, y), Not(x)))
#> $expr
#> b | ~a
#> 
#> $mapping
#> $mapping$a
#> x
#> 
#> $mapping$b
#> y
#> 
#> 
```
