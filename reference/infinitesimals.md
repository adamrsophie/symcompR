# Infinitesimals of a first order ODE

`infinitesimals(eq, func)` finds infinitesimals `xi(x, f(x))` and
`eta(x, f(x))` of the Lie point symmetries of the first order ODE `eq`,
as SymPy's `infinitesimals()`: the heuristics `abaco1_simple`,
`abaco1_product`, `abaco2_similar`, `abaco2_unique_unknown`,
`abaco2_unique_general`, `linear`, `function_sum`, `bivariate` and `chi`
are tried in this order.

## Usage

``` r
infinitesimals(eq, func = NULL, order = NULL, hint = "default", match = NULL)
```

## Arguments

- eq:

  The ODE, an expression or an equation.

- func:

  The unknown function; detected when `NULL`.

- order:

  The order of the ODE; computed when `NULL`.

- hint:

  `"default"` for the first heuristic that succeeds, `"all"` for the
  infinitesimals of every heuristic, or the name of one heuristic.

- match:

  Internal: the right hand side `h` and the variable `y` used by the
  `lie_group` solver.

## Value

A list of dictionaries keyed by `xi(x, f(x))` and `eta(x, f(x))`.

## Examples

``` r
x <- Symbol("x"); f <- Function("f")
infinitesimals(Eq(Derivative(f(x), x), x^2 * f(x)), f(x))
#> [[1]]
#> {xi(x, f(x)): 0, eta(x, f(x)): exp(x^3/3)} 
#> 
```
