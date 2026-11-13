# Solve a linear system as a set

`linsolve(system, ...)` solves a system of linear equations and returns
the `FiniteSet` of one `Tuple` (or the empty set), as SymPy's
`linsolve()`. Free unknowns stay as parameters. The system may be a list
of equations or expressions, a list `(A, b)` of matrices, or an
augmented matrix; the symbols follow one by one or as a list, and for
matrix input they default to `tau0, tau1, ...`.

## Usage

``` r
linsolve(system, ...)
```

## Arguments

- system:

  The linear system.

- ...:

  The symbols, or one list of them.

## Value

A set of tuples.

## Examples

``` r
x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
linsolve(list(x + y + z - 1, x + y + 2*z - 3), list(x, y, z))
#> {(-y - 1, y, 2)}
linsolve(Matrix(list(list(1, 2, 3), list(4, 5, 6))))
#> {(-1, 2)}
```
