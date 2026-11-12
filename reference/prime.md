# Generating primes

`prime(nth)` is the `nth` prime, `nextprime(n, ith)` the `ith` prime
after `n`, `prevprime(n)` the largest prime below `n`,
`primerange(a, b)` the primes in `[a, b)` (with one argument, the primes
below `a`), and `randprime(a, b)` a random prime in `[a, b)` using the R
random number generator, so
[`set.seed()`](https://rdrr.io/r/base/Random.html) applies.
`primorial(n, nth)` is the product of the first `n` primes, or with
`nth = FALSE` of the primes up to `n`. `composite(nth)` is the `nth`
composite and `compositepi(n)` counts the composites up to `n`.
`cycle_length(f, x0)` runs Brent's cycle finder on the sequence
`x0, f(x0), ...`.

## Usage

``` r
prime(nth)

nextprime(n, ith = 1)

prevprime(n)

primerange(a, b = NULL)

randprime(a, b)

primorial(n, nth = TRUE)

cycle_length(f, x0, nmax = NULL, values = FALSE)

composite(nth)

compositepi(n)
```

## Arguments

- nth, ith:

  Positive index.

- n, a, b:

  Integers. `n`, `a` and `b` of `prevprime()`, `primerange()` and
  `randprime()` may be numbers or `gmp` big integers; see
  [`isprime()`](https://adamrsophie.github.io/symcompR/reference/isprime.md).

- f:

  A function of one argument.

- x0:

  Starting value.

- nmax:

  Stop after this many iterations; `NULL` runs to the cycle.

- values:

  If `TRUE`, return the list of visited values instead.

## Value

An integer or vector, as described above. `cycle_length()` returns
`list(lambda, mu)` with `mu` missing when `nmax` stopped the search, or
the list of values. `randprime()` returns `NULL` when `a >= b`.

## Details

Integer results are plain doubles while they are exact (below `2^53`)
and `gmp` big integers above that.

## Examples

``` r
prime(25)
#> [1] 97
nextprime(100, 2)
#> [1] 103
prevprime(100)
#> [1] 97
primerange(2, 30)
#>  [1]  2  3  5  7 11 13 17 19 23 29
primorial(5)
#> [1] 2310
composite(10)
#> [1] 18
cycle_length(function(x) (x^2 + 1) %% 7, 2)
#> $lambda
#> [1] 1
#> 
#> $mu
#> [1] 1
#> 
```
