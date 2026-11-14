# Building blocks of the recurrence solver

The recurrence is `sum(coeffs[[i]] * y(n + i - 1)) = f`, with polynomial
coefficients in `n`. `rsolve_poly()` finds its polynomial solutions
(Abramov's algorithm), `rsolve_ratio()` its rational solutions and
`rsolve_hyper()` its hypergeometric solutions (Petkovsek's algorithm).
Free constants are `C0`, `C1`, ... The result is `NULL` when there is no
solution of that kind.

## Usage

``` r
rsolve_poly(coeffs, f, n, shift = 0L)

rsolve_ratio(coeffs, f, n)

rsolve_hyper(coeffs, f, n)
```

## Arguments

- coeffs:

  A list of coefficients, lowest shift first.

- f:

  The inhomogeneous term.

- n:

  The index symbol.

- shift:

  First index for the constant names (`rsolve_poly()` only).

## Value

An expression or `NULL`.

## Examples

``` r
n <- Symbol("n", integer = TRUE)
rsolve_poly(list(-1, 1), 4*n^3, n)
#> C0 + n^4 - 2*n^3 + n^2
# \donttest{
rsolve_ratio(list(-2*n^3 + n^2 + 2*n - 1, 2*n^3 + n^2 - 6*n,
                  -2*n^3 - 11*n^2 - 18*n - 9, 2*n^3 + 13*n^2 + 22*n + 8), 0, n)
#> C0*(2*n - 3)/(2*(n^2 - 1))
rsolve_hyper(list(-1, -1, 1), 0, n)
#> C0*(1/2 - sqrt(5)/2)^n + C1*(1/2 + sqrt(5)/2)^n
# }
```
