# Solve a system of nonlinear equations as a set

`nonlinsolve(system, ...)` solves a system of equations for the given
symbols and returns a `FiniteSet` of `Tuple`s, as SymPy's
`nonlinsolve()`: polynomial systems go through a Groebner basis
(positive-dimensional ones leave free unknowns as parameters), other
equations are solved by
[`substitution()`](https://adamrsophie.github.io/symcompR/reference/substitution.md),
with complex and periodic (`ImageSet`) solutions, and a `ConditionSet`
when nothing can be solved.

## Usage

``` r
nonlinsolve(system, ...)
```

## Arguments

- system:

  A list of equations or expressions.

- ...:

  The symbols, or one list of them.

## Value

A set of tuples, or a `ConditionSet`.

## Examples

``` r
x <- Symbol("x"); y <- Symbol("y")
nonlinsolve(list(x^2 + y^2 - 1, x - y), list(x, y))
#> {(-sqrt(2)/2, -sqrt(2)/2), (sqrt(2)/2, sqrt(2)/2)}
nonlinsolve(list(x*y - 1, x - 2), list(x, y))
#> {(2, 1/2)}
```
