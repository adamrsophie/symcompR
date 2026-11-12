# Perfect, abundant, deficient, amicable and Carmichael numbers

`abundance(n)` is `sigma(n) - 2n`. `is_perfect()`, `is_abundant()` and
`is_deficient()` test its sign; `is_amicable(m, n)` tests an amicable
pair and `is_carmichael(n)` a Carmichael number.
`find_carmichael_numbers_in_range(x, y)` lists the Carmichael numbers in
`[x, y)` and `find_first_n_carmichaels(n)` the first `n`.
`mersenne_prime_exponent(nth)` is the exponent of the `nth` Mersenne
prime, from the list of the 52 known ones.

## Usage

``` r
mersenne_prime_exponent(nth)

is_perfect(n)

abundance(n)

is_abundant(n)

is_deficient(n)

is_amicable(m, n)

is_carmichael(n)

find_carmichael_numbers_in_range(x, y)

find_first_n_carmichaels(n)
```

## Arguments

- nth:

  Index from 1 to 52.

- n, m:

  Integers.

- x, y:

  Range bounds, `0 <= x <= y`.

## Value

A logical, an integer or a vector.

## Examples

``` r
is_perfect(28)
#> [1] TRUE
abundance(12)
#> [1] 4
is_amicable(220, 284)
#> [1] TRUE
find_first_n_carmichaels(3)
#> [1]  561 1105 1729
mersenne_prime_exponent(10)
#> [1] 89
```
