# Perfect powers

`perfect_power(n)` returns `c(b, e)` with `b^e == n` and `e > 1` as
large as possible, or `FALSE`. With `big = FALSE` the smallest prime
exponent is used. `candidates` restricts the exponents tried;
`factor = FALSE` skips the trial division by small primes. A rational
`n` returns `list(b, e)` with `b` a `sym` rational.

## Usage

``` r
perfect_power(n, candidates = NULL, big = TRUE, factor = TRUE)
```

## Arguments

- n:

  An integer or rational.

- candidates:

  Exponents to try, or `NULL`.

- big:

  Prefer the largest exponent.

- factor:

  Try small prime factors first.

## Value

`c(b, e)`, `list(b, e)` for rationals, or `FALSE`.

## Examples

``` r
perfect_power(16)
#> [1] 2 4
perfect_power(16, big = FALSE)
#> [1] 4 2
perfect_power(Rational(4, 9))
#> [[1]]
#> 2/3
#> 
#> [[2]]
#> [1] 2
#> 
perfect_power(-8)
#> [1] -2  3
```
