# Tschirnhausen transformation of a polynomial

Given a monic irreducible polynomial `T` over the integers, find another
such polynomial `U` defining the same number field, together with the
polynomial `A` that maps the roots of `T` to the roots of `U`. This is
Algorithm 6.3.4 of Cohen, as in SymPy's
`sympy.polys.numberfields.galoisgroups.tschirnhausen_transformation`.

## Usage

``` r
tschirnhausen_transformation(
  T,
  max_coeff = 10,
  max_tries = 30,
  history = NULL,
  fixed_order = TRUE
)
```

## Arguments

- T:

  A monic `Poly` over the integers.

- max_coeff:

  Bound on the coefficients of random candidates.

- max_tries:

  Number of candidates to try.

- history:

  `NULL` or an environment used as a set of coefficient vectors that `U`
  must avoid; `T` is added to it.

- fixed_order:

  Use the deterministic order.

## Value

A list `(A, U)` of `Poly` objects; `U` is in a dummy variable `X`.

## Details

With `fixed_order = TRUE` the candidates `A` are tried in SymPy's fixed
order, from small degrees and coefficients upward; otherwise they are
drawn at random from the Python compatible generator.

## Examples

``` r
x <- Symbol("x")
tschirnhausen_transformation(Poly(x^2 - 2, x))
#> [[1]]
#> Poly(x^2 + x + 1, x, domain='ZZ')
#> 
#> [[2]]
#> Poly(_X^2 - 6*_X + 7, _X, domain='ZZ')
#> 
```
