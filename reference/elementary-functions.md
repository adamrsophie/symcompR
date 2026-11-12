# Elementary and complex functions

Symbolic versions of trigonometric, hyperbolic, exponential, complex and
piecewise elementary functions. The R generics
[`sin()`](https://rdrr.io/r/base/Trig.html),
[`cos()`](https://rdrr.io/r/base/Trig.html),
[`exp()`](https://rdrr.io/r/base/Log.html),
[`log()`](https://rdrr.io/r/base/Log.html),
[`sqrt()`](https://rdrr.io/r/base/MathFun.html),
[`abs()`](https://rdrr.io/r/base/MathFun.html),
[`Re()`](https://rdrr.io/r/base/complex.html),
[`Im()`](https://rdrr.io/r/base/complex.html), `sym_Mod()`,
[`Arg()`](https://rdrr.io/r/base/complex.html),
[`Conj()`](https://rdrr.io/r/base/complex.html),
[`floor()`](https://rdrr.io/r/base/Round.html),
[`ceiling()`](https://rdrr.io/r/base/Round.html),
[`max()`](https://rdrr.io/r/base/Extremes.html) and
[`min()`](https://rdrr.io/r/base/Extremes.html) also work on `sym`
objects.

## Usage

``` r
cot(x)

sec(x)

csc(x)

acot(x)

asec(x)

acsc(x)

coth(x)

sech(x)

csch(x)

acoth(x)

asech(x)

acsch(x)

sinc(x)

Abs(x)

re(x)

im(x)

conjugate(x)

arg(x)

frac(x)

Heaviside(x, h0 = NULL)

DiracDelta(x, k = 0)

KroneckerDelta(x, y)

LambertW(x, k = 0)

Max(...)

Min(...)

root(x, n, k = 0)

cbrt(x)

real_root(x, n)

sym_atan2(y, x)

sym_Mod(x, y, ...)
```

## Arguments

- x, y:

  Expressions.

- h0:

  Value of the Heaviside function at zero.

- k:

  Branch index or derivative order.

- ...:

  Expressions.

- n:

  Root index.

## Value

A `sym` expression.

## Examples

``` r
x <- Symbol("x")
cot(sym_pi / 4)
#> 1
Abs(-3 * x)
#> 3*Abs(x)
Max(1, 2, x)
#> Max(2, x)
```
