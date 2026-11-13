# Check solutions of ODEs

`checkodesol(ode, sol, func)` substitutes the solution `sol` into the
ODE `ode` and returns `(TRUE, 0)` when it is a solution, or
`(FALSE, residual)`, as SymPy's `checkodesol()`. A list of solutions
gives a list of results, and a list of equations is checked as a system
with
[`checksysodesol()`](https://adamrsophie.github.io/symcompR/reference/checksysodesol.md).

## Usage

``` r
checkodesol(ode, sol, func = NULL, order = "auto", solve_for_func = TRUE)
```

## Arguments

- ode:

  The ODE, an expression or an equation.

- sol:

  A solution, as an equation or the value of the function, or a list of
  them.

- func:

  The unknown function; detected when `NULL`.

- order:

  The order of the ODE, or `"auto"`.

- solve_for_func:

  Whether to try solving an implicit solution for the function first.

## Value

A tuple `(TRUE, 0)` or `(FALSE, residual)`, or a list of them.

## Examples

``` r
x <- Symbol("x"); f <- Function("f")
checkodesol(Derivative(f(x), x) - f(x), Eq(f(x), exp(x)), f(x))
#> (True, 0) 
```
