# Classify a system of ODEs

`classify_sysode(eq, funcs)` returns the properties of the system of
ODEs `eq` that dsolve() uses to pick a solver, as SymPy's
`classify_sysode()`: the number of equations, the equations as
expressions, the functions, their orders, the coefficients of the
functions and their derivatives, the linearity, and the type of the
system.

## Usage

``` r
classify_sysode(eq, funcs = NULL)
```

## Arguments

- eq:

  A list of ODEs.

- funcs:

  A list of the unknown functions, or `NULL` to find them.

## Value

A dictionary keyed by the names of the properties.

## Examples

``` r
t <- Symbol("t"); x <- Function("x"); y <- Function("y")
classify_sysode(list(Eq(Derivative(x(t), t), y(t)^3 * x(t)), Eq(Derivative(y(t), t), y(t)^5)))
#> {no_of_equation: 2, eq: [-x(t)*y(t)^3 + Derivative(x(t), t), -y(t)^5 + Derivative(y(t), t)], func: [x(t), y(t)], order: {x(t): 1, y(t): 1}, func_coeff: {(0, x(t), 0): -y(t)^3, (0, x(t), 1): 1, (0, y(t), 0): 0, (0, y(t), 1): 0, (1, x(t), 0): 0, (1, x(t), 1): 0, (1, y(t), 0): 0, (1, y(t), 1): 1}, is_linear: False, type_of_equation: type1} 
```
