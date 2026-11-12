# Polynomial congruences

`quadratic_congruence(a, b, c, n)` solves `a x^2 + b x + c = 0 (mod n)`.
`polynomial_congruence(expr, m)` solves a univariate polynomial with
integer coefficients modulo `m`. `binomial_mod(n, m, k)` is
`choose(n, m) %% k`, by Granville's extension of Lucas' theorem.

## Usage

``` r
quadratic_congruence(a, b, c, n)

polynomial_congruence(expr, m)

binomial_mod(n, m, k)
```

## Arguments

- a, b, c:

  Coefficients.

- n, m, k:

  Integers; `n` and `m` of `quadratic_congruence()` and
  `polynomial_congruence()` are moduli.

- expr:

  A `sym` polynomial or a string.

## Value

A vector of solutions in increasing order, or an integer.

## Examples

``` r
quadratic_congruence(2, 5, 3, 7)
#> [1] 2 6
polynomial_congruence("x^6 - 2*x^5 - 35", 6125)
#> [1] 3257
binomial_mod(10, 2, 6)
#> [1] 3
```
