# Order of a differential equation

`ode_order(expr, func)` returns the order of the differential equation
`expr` in the function `func`, as SymPy's `ode_order()`: a derivative of
an expression adds its own order to the highest order among the
arguments of that expression.

## Usage

``` r
ode_order(expr, func)
```

## Arguments

- expr:

  An expression or equation.

- func:

  The function, such as `f(x)`.

## Value

An integer.

## Examples

``` r
x <- Symbol("x"); f <- Function("f")
ode_order(Derivative(f(x), x, 3) + f(x), f(x))
#> [1] 3
```
