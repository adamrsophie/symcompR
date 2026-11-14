# Numerical solution of equations

`nsolve(f, x, x0)` finds a root near `x0` as SymPy's `nsolve()`:
mpmath's secant method for one equation and its damped Newton method for
systems, at `prec` significant digits plus a guard of 20 bits, with the
root verified and rounded to `prec` digits. For systems give lists of
equations, symbols and starting values.

## Usage

``` r
nsolve(
  f,
  x,
  x0,
  tol = NULL,
  maxit = NULL,
  prec = 15,
  maxsteps = NULL,
  verify = TRUE,
  dict = FALSE,
  solver = NULL
)
```

## Arguments

- f:

  An equation or list of equations.

- x:

  A symbol or list of symbols.

- x0:

  Starting value(s).

- tol:

  Tolerance, or `NULL` for mpmath's default of `2^10` times the working
  machine epsilon.

- maxit:

  Maximum iterations, or `NULL` for the solver's default (30 secant
  steps, 10 Newton steps).

- prec:

  Number of significant digits.

- maxsteps:

  SymPy's name for `maxit`.

- verify:

  If `FALSE`, return the last iterate without checking that it is a root
  within the tolerance.

- dict:

  If `TRUE`, return a list holding one named list of the solution, as
  SymPy's `dict=True`.

- solver:

  The mpmath solver: `"secant"` (one equation; give `x0` as two values
  for the two starting points) or `"mdnewton"` (systems).

## Value

A Float, or a named list of Floats (a list of one with `dict`).

## Examples

``` r
x <- Symbol("x")
nsolve(cos(x) - x, x, 1)
#> 0.739085133215161
nsolve(x^2 - 2, x, 1, prec = 30)
#> 1.41421356237309504880168872421
nsolve(x^2 - 2, x, c(1, 2))
#> 1.41421356237310
```
