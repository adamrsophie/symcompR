# Number theoretic functions

Symbolic functions that evaluate on integers and stay unevaluated on
symbols, using their assumptions where possible (`totient(p)` is `p - 1`
for a prime symbol `p`): Euler's `totient(n)`, Carmichael's
`reduced_totient(n)`, `divisor_sigma(n, k)` (the sum of the `k`-th
powers of the divisors) and its unitary analogue `udivisor_sigma(n, k)`,
`primenu(n)` and `primeomega(n)` (distinct and repeated prime factors),
the Moebius function `mobius(n)`, the Legendre, Jacobi and Kronecker
symbols, and the prime counting function `primepi(n)`. `totient()` and
`reduced_totient()` also accept a factorization as returned by
[`factorint()`](https://adamrsophie.github.io/symcompR/reference/factorint.md).

## Usage

``` r
totient(n)

reduced_totient(n)

divisor_sigma(n, k = NULL)

udivisor_sigma(n, k = NULL)

primenu(n)

primeomega(n)

mobius(n)

legendre_symbol(a, p)

jacobi_symbol(m, n)

kronecker_symbol(a, n)

primepi(n)
```

## Arguments

- n, a, m, p:

  Arguments.

- k:

  Power; `NULL` (the default) means 1.

## Value

A `sym` expression.

## Examples

``` r
totient(25)
#> 20
totient(Symbol("p", prime = TRUE))
#> p - 1
divisor_sigma(18, 2)
#> 455
mobius(30)
#> -1
legendre_symbol(3, 7)
#> -1
jacobi_symbol(45, 77)
#> -1
primepi(100)
#> 25
```
