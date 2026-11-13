# Dot product of rotational gradients of spherical harmonics

`dot_rot_grad_Ynm(j, p, l, m, theta, phi)` returns the dot product of
the rotational gradients of two spherical harmonics as an unevaluated
sum over Wigner 3j symbols;
[`doit()`](https://adamrsophie.github.io/symcompR/reference/doit.md)
evaluates it. After `sympy.physics.wigner.dot_rot_grad_Ynm`.

## Usage

``` r
dot_rot_grad_Ynm(j, p, l, m, theta, phi)
```

## Arguments

- j, p:

  Degree and order of the first spherical harmonic.

- l, m:

  Degree and order of the second spherical harmonic.

- theta, phi:

  The angles.

## Value

A `sym` expression.

## Examples

``` r
doit(dot_rot_grad_Ynm(3, 2, 2, 0, Symbol("theta"), Symbol("phi")))
#> 3*sqrt(55)*Ynm(5, 2, theta, phi)/(11*sqrt(pi))
```
