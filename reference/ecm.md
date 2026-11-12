# Lenstra's elliptic curve factoring

`ecm(n)` returns the sorted prime factors of `n`, after trial division
by the primes below 100000. `B1` and `B2` are the even stage 1 and stage
2 bounds and `max_curve` the number of curves tried per factor. Random
curves come from a stream private to the call, seeded by `seed`.

## Usage

``` r
ecm(n, B1 = 10000, B2 = 1e+05, max_curve = 200, seed = 1234)
```

## Arguments

- n:

  An integer.

- B1, B2:

  Even stage bounds.

- max_curve:

  Maximum number of curves.

- seed:

  Seed for the curves.

## Value

A vector of distinct primes.

## Examples

``` r
ecm("25645121643901801")
#> [1]    5394769 4753701529
```
