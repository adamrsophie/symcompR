# Check solutions of systems of ODEs

`checksysodesol(eqs, sols)` substitutes the solutions `sols` (a list of
equations, one per function) into the system `eqs`, as SymPy's
`checksysodesol()`, and returns `(TRUE, residuals)` when all residuals
vanish, `(FALSE, residuals)` otherwise.

## Usage

``` r
checksysodesol(eqs, sols, func = NULL)
```

## Arguments

- eqs:

  A list of ODEs.

- sols:

  A list of solutions.

- func:

  The unknown functions; detected when `NULL`.

## Value

A tuple of a logical and the list of residuals.

## Examples

``` r
t <- Symbol("t"); f <- Function("f"); g <- Function("g")
checksysodesol(list(Eq(Derivative(f(t), t), g(t)), Eq(Derivative(g(t), t), -f(t))),
               list(Eq(f(t), sin(t)), Eq(g(t), cos(t))))
#> (True, [0, 0]) 
```
