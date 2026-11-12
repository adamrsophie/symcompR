# Primality tests

`isprime()` tests one integer using GMP (Baillie-PSW plus Miller-Rabin
rounds), so its answer is certain below `2^64` and a probable prime
above, as in SymPy. `is_gaussian_prime()` tests a Gaussian integer. The
remaining functions are the individual pseudoprime tests: Fermat, Euler
and Euler-Jacobi to base `a`, Miller-Rabin `mr()` over `bases`, the
Selfridge Lucas tests, `proth_test()` for Proth numbers and
`is_mersenne_prime()`.
[`is_square()`](https://adamrsophie.github.io/symcompR/reference/dot.md)
tests for a perfect square when given an integer.

## Usage

``` r
isprime(n)

is_gaussian_prime(num)

is_fermat_pseudoprime(n, a)

is_euler_pseudoprime(n, a)

is_euler_jacobi_pseudoprime(n, a)

mr(n, bases)

is_lucas_prp(n)

is_strong_lucas_prp(n)

is_extra_strong_lucas_prp(n)

proth_test(n)

is_mersenne_prime(n)
```

## Arguments

- n:

  An integer.

- num:

  A Gaussian integer: an R complex number, a `sym` expression or an
  integer.

- a:

  Base of the test.

- bases:

  A vector or list of integer bases.

## Value

`TRUE` or `FALSE`.

## Details

Integer arguments may be numbers, `gmp` big integers, decimal strings or
`sym` integers. Doubles at or above `2^53` are refused, because they may
have lost digits already; pass a string or big integer.

## Examples

``` r
isprime(97)
#> [1] TRUE
isprime("170141183460469231731687303715884105727")
#> [1] TRUE
is_gaussian_prime(complex(real = 2, imaginary = 3))
#> [1] TRUE
mr(1373651, c(2, 3))
#> [1] FALSE
```
