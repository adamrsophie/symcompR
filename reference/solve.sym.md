# Solve equations

`solve(f, x)` solves `f = 0` (or an
[`Eq()`](https://adamrsophie.github.io/symcompR/reference/Eq.md)) for
`x`, returning a list of solutions. With a list of equations and a list
of symbols it solves the system and returns a list of solutions, each a
named list mapping symbol names to values. Polynomial equations are
solved exactly (by factorization, radicals or
[`CRootOf()`](https://adamrsophie.github.io/symcompR/reference/CRootOf.md)),
simple transcendental equations by inversion, linear systems by
Gauss-Jordan elimination (free parameters stay symbolic) and polynomial
systems by propagation and Groebner bases. Inequalities are passed to
[`solve_univariate_inequality()`](https://adamrsophie.github.io/symcompR/reference/solve_univariate_inequality.md).
This is a method for the base R generic
[`base::solve()`](https://rdrr.io/r/base/solve.html).

## Usage

``` r
# S3 method for class 'sym'
solve(a, b, ..., dict = FALSE, check = TRUE)

# S3 method for class 'character'
solve(a, b, ...)

# S3 method for class 'list'
solve(a, b, ..., dict = TRUE, check = TRUE)
```

## Arguments

- a:

  An equation, expression, list of them, or a character string.

- b:

  A symbol or list of symbols; by default all free symbols.

- ...:

  Ignored.

- dict:

  For a single equation, return named lists as for systems.

- check:

  Drop solutions that make a denominator vanish.

## Value

A list of solutions.

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
#> [[1]]
#> [[1]]$x
#> 2
#> 
#> [[1]]$y
#> 1
#> 
#> 
solve(list(x^2 + y^2 - 5, x - y + 1), list(x, y))
#> [[1]]
#> [[1]]$x
#> 1
#> 
#> [[1]]$y
#> 2
#> 
#> 
#> [[2]]
#> [[2]]$x
#> -2
#> 
#> [[2]]$y
#> -1
#> 
#> 
```
