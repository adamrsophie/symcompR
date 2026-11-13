# Clebsch-Gordan coefficients and Wigner symbols as expressions

`CG(j1, m1, j2, m2, j3, m3)` is the unevaluated Clebsch-Gordan
coefficient and `Wigner3j()`, `Wigner6j()` and `Wigner9j()` the Wigner
symbols of `sympy.physics.quantum.cg`;
[`doit()`](https://adamrsophie.github.io/symcompR/reference/doit.md)
evaluates them for numeric arguments. `cg_simp()` simplifies sums and
products of Clebsch-Gordan coefficients with the Varshalovich
identities.

## Usage

``` r
Wigner3j(j1, m1, j2, m2, j3, m3)

CG(j1, m1, j2, m2, j3, m3)

Wigner6j(j1, j2, j12, j3, j, j23)

Wigner9j(j1, j2, j12, j3, j4, j34, j13, j24, j)

cg_simp(e)
```

## Arguments

- j1, m1, j2, m2, j3, m3, j12, j23, j, j4, j34, j13, j24:

  Angular momenta and projections.

- e:

  An expression.

## Value

An expression.

## Examples

``` r
cg <- CG(1, 0, 1, 0, 2, 0)
doit(cg)
#> sqrt(6)/3
doit(Wigner3j(2, 0, 6, 0, 4, 0))
#> sqrt(715)/143
```
