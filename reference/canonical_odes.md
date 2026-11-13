# Canonical form of a system of ODEs

`canonical_odes(eqs, funcs, t)` solves the system `eqs` for the highest
derivative of each function, as SymPy's `canonical_odes()`; there is one
system for each solution.

## Usage

``` r
canonical_odes(eqs, funcs, t)
```

## Arguments

- eqs:

  A list of ODEs.

- funcs:

  A list of the unknown functions.

- t:

  The independent variable.

## Value

A list of systems, each a list of equations.

## Examples

``` r
x <- Symbol("x"); f <- Function("f"); g <- Function("g")
canonical_odes(list(Eq(Derivative(f(x), x), g(x)), Eq(Derivative(g(x), x) + f(x), 0)),
               list(f(x), g(x)), x)
#> [[1]]
#> [[1]][[1]]
#> Eq(Derivative(f(x), x), g(x))
#> 
#> [[1]][[2]]
#> Eq(Derivative(g(x), x), -f(x))
#> 
#> 
```
