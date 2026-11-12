# Numerical solution of equations

`nsolve(f, x, x0)` finds a root near `x0` by Newton's method (secant
fallback); for systems give lists of equations, symbols and starting
values.

## Usage

``` r
nsolve(f, x, x0, tol = 1e-14, maxit = 100)
```

## Arguments

- f:

  An equation or list of equations.

- x:

  A symbol or list of symbols.

- x0:

  Starting value(s).

- tol:

  Tolerance.

- maxit:

  Maximum iterations.

## Value

A Float or a named list of Floats.

## Examples

``` r
x <- Symbol("x")
nsolve(cos(x) - x, x, 1)
#> 0.739085133215161
```
