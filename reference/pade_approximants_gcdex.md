# Pade approximants by the extended Euclidean algorithm

`pade_approximants_gcdex(f, order)` lists the pairs
`list(numerator, denominator)` of the Pade approximants of the
polynomial `f` (a
[`Poly()`](https://adamrsophie.github.io/symcompR/reference/Poly.md) or
an expression in one variable) truncated at degree `order`, from the
steps of
[`gcdex_steps()`](https://adamrsophie.github.io/symcompR/reference/gcdex_steps.md)
on `x^(order + 1)` and `f`. `pade_approximant_gcdex(f, m, n)` returns
the approximant with numerator degree at most `m` and denominator degree
at most `n` (`n = m` by default), as SymPy's functions of the same
names.

## Usage

``` r
pade_approximants_gcdex(f, order)

pade_approximant_gcdex(f, m, n = NULL)

pade_approximant(f, x, x0 = 0, m = 6, n = NULL)
```

## Arguments

- f:

  A univariate polynomial.

- order:

  The truncation order.

- m, n:

  Degrees of the numerator and the denominator.

- x:

  The variable.

- x0:

  The expansion point.

## Value

A list of `Poly` pairs, or one pair.

## Examples

``` r
x <- Symbol("x")
f <- 1 + x + x^2/2 + x^3/6 + x^4/24
pade_approximant_gcdex(f, 2)
#> [[1]]
#> Poly(1/4*x^2 + 3/2*x + 3, x, domain='QQ')
#> 
#> [[2]]
#> Poly(1/4*x^2 - 3/2*x + 3, x, domain='QQ')
#> 
```
