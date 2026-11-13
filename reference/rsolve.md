# Solve linear recurrences

`rsolve(f, y)` solves the linear recurrence `f = 0` (or an
[`Eq()`](https://adamrsophie.github.io/symcompR/reference/Eq.md)) with
polynomial or rational coefficients for the sequence term `y`, such as
`y(n)`, as SymPy's `rsolve()`. The inhomogeneous part must be a sum of
hypergeometric terms. Arbitrary constants are `C0`, `C1`, ... `init`
gives initial values as a list of `list(key, value)` pairs, where a key
is an index or a term such as `y(0)`, or as a list of values for the
indices `0, 1, ...`. The result is `NULL` when no solution is found.

## Usage

``` r
rsolve(f, y, init = NULL, ...)
```

## Arguments

- f:

  The recurrence: an expression equal to zero or an equation.

- y:

  The unknown term, for example `y(n)`.

- init:

  Optional initial values.

- ...:

  For compatibility with the older calling form
  `rsolve(eq, y(n), n, init = ...)`: `init` may be given here by name.

## Value

An expression or `NULL`.

## Examples

``` r
n <- Symbol("n", integer = TRUE); y <- Function("y")
rsolve(y(n + 2) - y(n + 1) - y(n), y(n), list(list(y(0), 0), list(y(1), 1)))
#> -sqrt(5)*(1/2 - sqrt(5)/2)^n/5 + sqrt(5)*(1/2 + sqrt(5)/2)^n/5
rsolve(y(n + 1) - 2 * y(n), y(n))
#> 2^n*C0
```
