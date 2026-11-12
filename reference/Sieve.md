# Sieve of Eratosthenes

`Sieve()` creates a list of primes, totients and Moebius values that
grows on demand, and `sieve` is the global instance the number theory
functions share. `sieve$extend(n)` grows the primes to at least `n`, and
`sieve$extend_to_no(i)` to at least `i` primes.
`sieve$primerange(a, b)`, `sieve$totientrange(a, b)` and
`sieve$mobiusrange(a, b)` return values for the half open range.
`sieve$search(n)` returns the indices of the primes bounding `n`, equal
if `n` is prime, and `sieve$contains(n)` tests membership. `sieve[i]` is
the `i`-th prime, counting from 1.

## Usage

``` r
Sieve(sieve_interval = 1e+06)

sieve
```

## Format

`sieve` is the shared `Sieve` instance.

## Arguments

- sieve_interval:

  Largest number of odd integers sieved in one block.

## Value

A `Sieve` environment.

## Examples

``` r
sieve[1:5]
#> [1]  2  3  5  7 11
sieve$primerange(10, 30)
#> [1] 11 13 17 19 23 29
sieve$totientrange(7, 12)
#> [1]  6  4  6  4 10
sieve$search(25)
#> [1]  9 10
```
