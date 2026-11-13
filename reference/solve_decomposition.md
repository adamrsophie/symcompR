# Solve a composite equation by decomposition

`solve_decomposition(f, symbol, domain)` writes `f` as a composition of
simpler functions with
[`decompogen()`](https://adamrsophie.github.io/symcompR/reference/more-solvers.md)
and solves the equations from the outside in, as SymPy's
`solve_decomposition()`.

## Usage

``` r
solve_decomposition(f, symbol, domain)
```

## Arguments

- f:

  An expression.

- symbol:

  The variable.

- domain:

  The domain, a set.

## Value

A set.

## Examples

``` r
x <- Symbol("x")
solve_decomposition(exp(2*x) - 3*exp(x) + 2, x, Reals)
#> {0, log(2)}
```
