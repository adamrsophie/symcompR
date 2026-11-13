# Convergence acceleration

`richardson(A, k, n, N)` estimates the limit of the sequence `A(k)` as
`k` goes to infinity by Richardson extrapolation of order `N` from the
terms `A(n), ..., A(n + N)`; `shanks(A, k, n, m)` applies the Shanks
transformation `m` times, as SymPy's functions of the same names.

## Usage

``` r
richardson(A, k, n, N)

shanks(A, k, n, m = 1)
```

## Arguments

- A:

  An expression in `k`.

- k:

  The index symbol.

- n:

  Starting index.

- N:

  Order of the extrapolation.

- m:

  Number of applications.

## Value

An expression.

## Examples

``` r
k <- Symbol("k")
richardson((1 + 1/k)^k, k, 10, 3)
#> 64760917612327084418953891207076443/23824634392675154069619793920000000
shanks(1 - 1/k, k, 4)
#> 7/8
```
