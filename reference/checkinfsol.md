# Check infinitesimals of an ODE

`checkinfsol(eq, infinitesimals, func)` checks whether each set of
infinitesimals (a named list with the entries `xi` and `eta`, or a list
of such lists) satisfies the determining PDE of the first order ODE
`eq`, as SymPy's `checkinfsol()`.

## Usage

``` r
checkinfsol(eq, infinitesimals, func = NULL, order = NULL)
```

## Arguments

- eq:

  The ODE.

- infinitesimals:

  The infinitesimals, as returned by
  [`infinitesimals()`](https://adamrsophie.github.io/symcompR/reference/infinitesimals.md).

- func:

  The unknown function; detected when `NULL`.

- order:

  The order of the ODE; computed when `NULL`.

## Value

A list of pairs `(TRUE, 0)` or `(FALSE, residual)`.

## Examples

``` r
x <- Symbol("x"); f <- Function("f")
eq <- Derivative(f(x), x) - x^2
checkinfsol(eq, infinitesimals(eq, f(x)), f(x))
#> [[1]]
#> (True, 0) 
#> 
```
