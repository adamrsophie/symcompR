# Matrix form of a system of linear ODEs

`linear_ode_to_matrix(eqs, funcs, t, order)` writes the system of linear
ODEs `eqs` of the given order as
`A_n X^(n) = A_(n-1) X^(n-1) + ... + A_0 X + b`, as SymPy's
`linear_ode_to_matrix()`.

## Usage

``` r
linear_ode_to_matrix(eqs, funcs, t, order)
```

## Arguments

- eqs:

  A list of linear ODEs.

- funcs:

  A list of the unknown functions.

- t:

  The independent variable.

- order:

  The order of the system.

## Value

A tuple of the list of matrices `A_n, ..., A_0` and the vector `b`.

## Examples

``` r
t <- Symbol("t"); f <- Function("f"); g <- Function("g")
linear_ode_to_matrix(list(Eq(Derivative(f(t), t), g(t)), Eq(Derivative(g(t), t), f(t))),
                     list(f(t), g(t)), t, 1)
#> ([Matrix([[1, 0], [0, 1]]), Matrix([[0, 1], [1, 0]])], Matrix([[0], [0]])) 
```
