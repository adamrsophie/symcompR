# Function aliases and small special functions

`sym_rf()` and `ff()` are
[`RisingFactorial()`](https://adamrsophie.github.io/symcompR/reference/special-functions.md)
and
[`FallingFactorial()`](https://adamrsophie.github.io/symcompR/reference/special-functions.md),
`ln()` is the natural logarithm, `E1(z)` the exponential integral
`expint(1, z)`, `Li(x)` the offset logarithmic integral `li(x) - li(2)`,
`erf2(x, y)` is `erf(y) - erf(x)`, `Id(x)` the identity, and
`LeviCivita()` (alias `Eijk()`) the permutation symbol of integer
arguments, extended as in SymPy by
`prod((a[j] - a[i]) / sym_factorial(i))`. `chebyshevt_root(n, k)` and
`chebyshevu_root(n, k)` give the `k`-th root (from 0) of the Chebyshev
polynomials. `transpose()` and `adjoint()` act on matrices.

## Usage

``` r
ff(x, k)

ln(x)

E1(z)

Li(x)

erf2(x, y)

Id(x)

LeviCivita(...)

Eijk(...)

chebyshevt_root(n, k)

chebyshevu_root(n, k)

transpose(M)

adjoint(M)

sym_rf(x, k, ...)
```

## Arguments

- x, y, z, k, n:

  Arguments.

- ...:

  Integer arguments.

- M:

  A matrix.

## Value

An expression or matrix.

## Examples

``` r
sym_rf(Symbol("x"), 3)
#> x*(x + 1)*(x + 2)
LeviCivita(1, 2, 3)
#> 1
chebyshevt_root(3, 0)
#> sqrt(3)/2
```
