# Simplify Bessel functions

`besselsimp(expr)` simplifies Bessel functions as SymPy's
`besselsimp()`: factors of `I` and `-1` (and their polar forms) are
taken out of the arguments of `besselj` and `besseli`, functions of
half-integer order are rewritten with trigonometric functions, functions
of integer order above one with functions of lower order, the recurrence
`J(a, x) = 2*(a + 1)/x*J(a + 1, x) - J(a + 2, x)` removes redundant
terms, and a changed expression is factored.

## Usage

``` r
besselsimp(expr)
```

## Arguments

- expr:

  An expression.

## Value

An expression.

## Examples

``` r
z <- Symbol("z"); nu <- Symbol("nu")
besselsimp(besselj(nu, z * polar_lift(-1)))
#> exp(I*pi*nu)*besselj(nu, z)
besselsimp(besseli(Rational(-1, 2), z))
#> sqrt(2)*cosh(z)/(sqrt(pi)*sqrt(z))
besselsimp(z * besseli(0, z) + z * besseli(2, z) / 2 + besseli(1, z))
#> 3*z*besseli(0, z)/2
```
