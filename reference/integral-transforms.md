# Unevaluated integral transforms

`LaplaceTransform(f, t, s)`, `FourierTransform(f, x, k)`,
`MellinTransform(f, x, s)`, `SineTransform()`, `CosineTransform()`,
`HankelTransform(f, r, k, nu)` and their inverses are unevaluated
transform objects;
[`doit()`](https://adamrsophie.github.io/symcompR/reference/doit.md)
computes them with
[`laplace_transform()`](https://adamrsophie.github.io/symcompR/reference/laplace_transform.md),
[`fourier_transform()`](https://adamrsophie.github.io/symcompR/reference/fourier_transform.md)
and the other transform functions. `InverseLaplaceTransform()` takes an
optional `plane` (ignored), and `InverseMellinTransform()` the strip
`a`, `b`.

## Usage

``` r
LaplaceTransform(f, t, s)

InverseLaplaceTransform(f, s, t, plane = NULL)

FourierTransform(f, x, k)

InverseFourierTransform(f, k, x)

MellinTransform(f, x, s)

InverseMellinTransform(f, s, x, a, b)

SineTransform(f, x, k)

InverseSineTransform(f, k, x)

CosineTransform(f, x, k)

InverseCosineTransform(f, k, x)

HankelTransform(f, r, k, nu)

InverseHankelTransform(f, k, r, nu)
```

## Arguments

- f:

  The function.

- t, s, x, k, r:

  Variables: the original first, the transformed second.

- plane:

  Ignored.

- a, b:

  The fundamental strip.

- nu:

  The Bessel order.

## Value

A transform object.

## Examples

``` r
t <- Symbol("t", positive = TRUE); s <- Symbol("s")
L <- LaplaceTransform(t^2, t, s)
doit(L)
#> 2/s^3
latex(L)
#> [1] "\\mathcal{L}_{t}\\left[t^{2}\\right]\\left(s\\right)"
```
