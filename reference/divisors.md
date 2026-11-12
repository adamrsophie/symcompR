# Divisors

`divisors(n)` lists the positive divisors of `n` in increasing order,
and `proper_divisors(n)` leaves out `n`. `divisor_count(n, modulus)`
counts the divisors that are multiples of `modulus`. `udivisors()` lists
the unitary divisors `d` with `gcd(d, n / d) = 1`, and `antidivisors()`
the integers that divide `n` as unevenly as possible. The `*_count()`
functions count them.

## Usage

``` r
divisors(n, proper = FALSE)

proper_divisors(n)

divisor_count(n, modulus = 1, proper = FALSE)

proper_divisor_count(n, modulus = 1)

udivisors(n)

udivisor_count(n)

antidivisors(n)

antidivisor_count(n)
```

## Arguments

- n:

  An integer.

- proper:

  Leave out `n`.

- modulus:

  Count only multiples of this.

## Value

A vector of divisors, or a count.

## Examples

``` r
divisors(24)
#> [1]  1  2  3  4  6  8 12 24
proper_divisors(24)
#> [1]  1  2  3  4  6  8 12
divisor_count(24, 2)
#> [1] 6
udivisors(60)
#> [1]  1  3  4  5 12 15 20 60
antidivisors(24)
#> [1]  7 16
```
