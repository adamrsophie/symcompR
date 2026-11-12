# Quadratic sieve

`qs_factor(N, prime_bound, M)` factors `N` with the self-initializing
quadratic sieve, using the primes below `prime_bound` as the factor base
and a sieve interval of `2 * M`; it returns a named vector of exponents
whose names are the factors found (not always prime). `qs()` returns the
factors alone. Random choices come from a stream private to the call,
seeded by `seed`.

## Usage

``` r
qs_factor(N, prime_bound, M, ERROR_TERM = 25, seed = 1234)

qs(N, prime_bound, M, ERROR_TERM = 25, seed = 1234)
```

## Arguments

- N:

  An integer greater than 1.

- prime_bound:

  Bound for the factor base.

- M:

  Half the sieve interval.

- ERROR_TERM:

  Slack in the sieve threshold.

- seed:

  Seed.

## Value

A named vector, or a vector of factors.

## Examples

``` r
qs(10009202107, 100, 2000)
#> [1] 100043 100049
```
