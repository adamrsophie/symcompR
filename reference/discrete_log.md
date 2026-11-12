# Discrete logarithm

`discrete_log(n, a, b)` finds `x` with `b^x = a (mod n)`, choosing among
trial multiplication, baby step giant step, Pollard's rho, index
calculus and Pohlig-Hellman as SymPy does. `order` is the order of `b`
if known, and `prime_order` states whether it is prime.

## Usage

``` r
discrete_log(n, a, b, order = NULL, prime_order = NULL)
```

## Arguments

- n:

  Modulus.

- a:

  Target.

- b:

  Base.

- order:

  Order of `b`, or `NULL`.

- prime_order:

  Whether `order` is prime, or `NULL`.

## Value

An integer.

## Examples

``` r
discrete_log(41, 15, 7)
#> [1] 3
discrete_log(587, 2^9, 2)
#> [1] 9
```
