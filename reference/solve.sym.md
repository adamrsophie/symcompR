# Solve equations

`solve(f, x)` solves `f = 0` (or an
[`Eq()`](https://adamrsophie.github.io/symcompR/reference/Eq.md)) for
`x` following SymPy's [`solve()`](https://rdrr.io/r/base/solve.html): a
single equation in one symbol gives a list of solutions, a linear system
a named list (a dict), other systems a list of solution tuples or, with
`dict = TRUE`, of named lists. Extra symbols may follow `b`, and SymPy's
flags (`dict`, `set`, `check`, `simplify`, `rational`, `exclude`,
`implicit`, `manual`, `cubics`, `quartics`, `quintics`, ...) are passed
by name. Solutions are simplified, checked and sorted as in SymPy.
Inequalities go to
[`reduce_inequalities()`](https://adamrsophie.github.io/symcompR/reference/inequality-solvers.md).
This is a method for the base R generic
[`base::solve()`](https://rdrr.io/r/base/solve.html).

## Usage

``` r
# S3 method for class 'sym'
solve(a, b, ...)

# S3 method for class 'character'
solve(a, b, ...)

# S3 method for class 'Poly'
solve(a, b, ...)

# S3 method for class 'list'
solve(a, b, ...)
```

## Arguments

- a:

  An equation, expression, list of them, or a character string.

- b:

  A symbol or list of symbols; by default all free symbols.

- ...:

  More symbols, then named flags as in SymPy.

## Value

A list of solutions, a named list, or a list of tuples.

## Examples

``` r
x <- Symbol("x"); y <- Symbol("y")
solve(x^2 - 4, x)
#> [[1]]
#> -2
#> 
#> [[2]]
#> 2
#> 
solve(Eq(exp(x), 3), x)
#> [[1]]
#> log(3)
#> 
solve(list(x + y - 3, x - y - 1), list(x, y))
#> $x
#> 2
#> 
#> $y
#> 1
#> 
solve(list(x^2 + y^2 - 5, x - y + 1), list(x, y))
#> [[1]]
#> (-2, -1) 
#> 
#> [[2]]
#> (1, 2) 
#> 
```
