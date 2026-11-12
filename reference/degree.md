# Degrees and coefficients of polynomials

`degree(f, x)` is the degree in `x` (`-Inf` for zero), `total_degree()`
the total degree, `degree_list()` the degree in each generator. `LC()`,
`LM()` and `LT()` give the leading coefficient, monomial and term in lex
order; `coeffs()` the nonzero coefficients; `all_coeffs()` the dense
coefficient list of a univariate polynomial (highest degree first);
`monoms()` the exponent vectors; `coeff()` the coefficient of a
monomial.

## Usage

``` r
degree(f, gen = NULL)

total_degree(f, ...)

degree_list(f, ...)

LC(f, ...)

LM(f, ...)

LT(f, ...)

coeffs(f, ...)

monoms(f, ...)

poly_terms(f, ...)

all_coeffs(f, ...)

coeff(f, x, n = 1)
```

## Arguments

- f:

  An expression or
  [`Poly()`](https://adamrsophie.github.io/symcompR/reference/Poly.md).

- gen:

  A generator.

- ...:

  Generators.

- x:

  A monomial such as `x^2*y`, for `coeff()`.

- n:

  Power of `x`, for `coeff()`.

## Value

A number, expression or list.

## Examples

``` r
x <- Symbol("x"); y <- Symbol("y")
degree(x^3 + x*y^4, x)
#> [1] 3
all_coeffs(Poly(2*x^3 - x + 5, x))
#> [[1]]
#> 2
#> 
#> [[2]]
#> 0
#> 
#> [[3]]
#> -1
#> 
#> [[4]]
#> 5
#> 
coeff(3*x^2*y + x, x, 2)
#> 3*y
```
