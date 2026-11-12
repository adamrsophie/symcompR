# Solve linear systems

`linsolve(eqs, syms)` solves a linear system, returning a named list
with free variables left as parameters, or `NULL` if inconsistent.

## Usage

``` r
linsolve(eqs, syms)
```

## Arguments

- eqs:

  A list of linear equations or expressions.

- syms:

  A list of symbols.

## Value

A named list, or `NULL`.

## Examples

``` r
x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
linsolve(list(x + y + z - 1, x + y + 2*z - 3), list(x, y, z))
#> $x
#> -y - 1
#> 
#> $z
#> 2
#> 
```
