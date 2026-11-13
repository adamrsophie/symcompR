# Appell sequences

The Bernoulli, central Bernoulli, Genocchi, Euler and Andre polynomials,
as SymPy's `bernoulli_poly()`, `bernoulli_c_poly()`, `genocchi_poly()`,
`euler_poly()` and `andre_poly()`. With `x = NULL` the variable is the
symbol `_0`, as SymPy writes it; with `polys = TRUE` the result is a
[`Poly()`](https://adamrsophie.github.io/symcompR/reference/Poly.md).

## Usage

``` r
bernoulli_poly(n, x = NULL, polys = FALSE)

bernoulli_c_poly(n, x = NULL, polys = FALSE)

genocchi_poly(n, x = NULL, polys = FALSE)

euler_poly(n, x = NULL, polys = FALSE)

andre_poly(n, x = NULL, polys = FALSE)
```

## Arguments

- n:

  The index, a nonnegative integer.

- x:

  The variable, or `NULL`.

- polys:

  Return a `Poly`.

## Value

An expression or a `Poly`.

## Examples

``` r
x <- Symbol("x")
bernoulli_poly(4, x)
#> x^4 - 2*x^3 + x^2 - 1/30
bernoulli_c_poly(4, x)
#> x^4 - 2*x^2 + 7/15
genocchi_poly(5, x)
#> -5*x^4 + 10*x^3 - 5*x
euler_poly(3, x)
#> x^3 - 3*x^2/2 + 1/4
andre_poly(4, x)
#> x^4 - 6*x^2 + 5
```
