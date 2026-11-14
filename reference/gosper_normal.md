# Gosper's algorithm helpers

`gosper_normal(f, g, n)` computes the Gosper normal form `(A, B, C)` of
the rational function `f/g` in `n`, and `gosper_term(f, n)` a
hypergeometric term `g` with `g(n + 1) - g(n) = f(n)`, or `NULL`, as
SymPy's functions of these names.

## Usage

``` r
gosper_normal(f, g, n, polys = TRUE)

gosper_term(f, n)
```

## Arguments

- f, g:

  Polynomials in `n` (for `gosper_term()`, `f` a hypergeometric term).

- n:

  The variable.

- polys:

  If `TRUE`, return `Poly` objects.

## Value

A list of three polynomials, or an expression or `NULL`.

## Examples

``` r
n <- Symbol("n")
gosper_normal(4*n + 5, 2*(4*n + 1)*(2*n + 3), n)
#> [[1]]
#> Poly(1/4, n, domain='QQ')
#> 
#> [[2]]
#> Poly(n + 3/2, n, domain='QQ')
#> 
#> [[3]]
#> Poly(n + 1/4, n, domain='QQ')
#> 
gosper_term(factorial(n), n)
#> NULL
```
