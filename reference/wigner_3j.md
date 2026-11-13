# Wigner 3j, 6j and 9j symbols and related coefficients

Exact values of the Wigner 3j, 6j and 9j symbols, the Clebsch-Gordan
coefficients, the Racah symbol, the Gaunt coefficient and its real
counterpart, after `sympy.physics.wigner`. Arguments are integers or
half integers (numbers or exact rationals); with `prec` the result is
evaluated numerically to that many digits.

## Usage

``` r
wigner_3j(j_1, j_2, j_3, m_1, m_2, m_3, prec = NULL)

clebsch_gordan(j_1, j_2, j_3, m_1, m_2, m_3, prec = NULL)

wigner_6j(j_1, j_2, j_3, j_4, j_5, j_6, prec = NULL)

racah(aa, bb, cc, dd, ee, ff, prec = NULL)

wigner_9j(j_1, j_2, j_3, j_4, j_5, j_6, j_7, j_8, j_9, prec = NULL)

gaunt(l_1, l_2, l_3, m_1, m_2, m_3, prec = NULL)

real_gaunt(l_1, l_2, l_3, mu_1, mu_2, mu_3, prec = NULL)
```

## Arguments

- j_1, j_2, j_3, j_4, j_5, j_6, j_7, j_8, j_9:

  Angular momenta.

- m_1, m_2, m_3:

  Projections.

- prec:

  Optional number of digits for a numerical result.

- aa, bb, cc, dd, ee, ff:

  Arguments of the Racah symbol.

- l_1, l_2, l_3:

  Degrees of the spherical harmonics.

- mu_1, mu_2, mu_3:

  Orders of the real spherical harmonics.

## Value

A `sym` number.

## Examples

``` r
wigner_3j(2, 6, 4, 0, 0, 0)
#> sqrt(715)/143
clebsch_gordan(Rational(3, 2), S_Half, 2, Rational(3, 2), S_Half, 2)
#> 1
wigner_6j(3, 3, 3, 3, 3, 3)
#> -1/14
wigner_9j(1, 1, 1, 1, 1, 1, 1, 1, 0, prec = 20)
#> 0.055555555555555555556
racah(3, 3, 3, 3, 3, 3)
#> -1/14
gaunt(1, 0, 1, 1, 0, -1)
#> -1/(2*sqrt(pi))
real_gaunt(2, 2, 2, 0, 2, 2)
#> -sqrt(5)/(7*sqrt(pi))
```
