# Solve linear recurrences

`rsolve(eq, yn, n)` solves a linear recurrence with constant
coefficients for the sequence term `yn`, such as `y(n)`, with polynomial
or geometric inhomogeneous terms. Arbitrary constants are `C0`, `C1`,
... `init` is a list of `list(index, value)` pairs.

## Usage

``` r
rsolve(eq, yn, n, init = NULL)
```

## Arguments

- eq:

  The recurrence, as an expression equal to zero.

- yn:

  The unknown term, for example `y(n)`.

- n:

  The index symbol.

- init:

  Optional initial values.

## Value

An expression.

## Examples

``` r
n <- Symbol("n", integer = TRUE); y <- Function("y")
rsolve(y(n + 2) - y(n + 1) - y(n), y(n), n,
       init = list(list(0, 0), list(1, 1)))
#> (1/2)^n*sqrt(5)*(-(1 - sqrt(5))^n + (1 + sqrt(5))^n)/5
```
