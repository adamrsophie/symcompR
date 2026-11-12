# Solve ordinary differential equations

`dsolve(eq, fx)` solves an ordinary differential equation for the
unknown function application `fx`, such as `f(x)`, returning an
[`Eq()`](https://adamrsophie.github.io/symcompR/reference/Eq.md).
Supported classes: separable and linear first order equations, Bernoulli
and exact equations, linear equations with constant coefficients of any
order (with undetermined coefficients for polynomial, exponential and
sinusoidal forcing) and Cauchy-Euler equations. Arbitrary constants are
`C1`, `C2`, ... `ics` applies initial conditions, given as a list of
`list(point, order, value)`. `ode_order()` returns the order and
`classify_ode()` the matching classes.

## Usage

``` r
dsolve(eq, fx, ics = NULL)

ode_order(eq, fx)

classify_ode(eq, fx)
```

## Arguments

- eq:

  The equation or expression equal to zero.

- fx:

  The unknown, for example `f(x)`.

- ics:

  Optional initial conditions.

## Value

An [`Eq()`](https://adamrsophie.github.io/symcompR/reference/Eq.md), or
a list of them for several solutions.

## Examples

``` r
x <- Symbol("x"); f <- Function("f")
dsolve(Derivative(f(x), x) - f(x), f(x))
#> Eq(f(x), C1*exp(x))
dsolve(Derivative(f(x), x, 2) + f(x), f(x))
#> Eq(f(x), C1*sin(x) + C2*cos(x))
```
