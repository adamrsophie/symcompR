# Pollard's rho and p - 1 factoring

`pollard_rho(n)` searches for a nontrivial factor of `n` with the
sequence `x -> x^2 + a` started at `s`, or with the function `F`,
retrying with random parameters. `pollard_pm1(n, B)` uses Pollard's p -
1 method with smoothness bound `B`. Both return a factor or `NULL`. The
random retries use a stream private to the call, seeded by `seed`, so
the results are reproducible but differ from SymPy's.

## Usage

``` r
pollard_rho(
  n,
  s = 2,
  a = 1,
  retries = 5,
  seed = 1234,
  max_steps = NULL,
  F = NULL
)

pollard_pm1(n, B = 10, a = 2, retries = 0, seed = 1234)
```

## Arguments

- n:

  An integer greater than 4 (`pollard_rho()`) or 3 (`pollard_pm1()`).

- s:

  Starting value.

- a:

  Constant of the default polynomial, or base for p - 1.

- retries:

  Number of retries with new random parameters.

- seed:

  Seed for the retries.

- max_steps:

  Maximum number of steps per attempt.

- F:

  A function of one big integer to use instead of `x^2 + a`.

- B:

  Smoothness bound, greater than 2.

## Value

A factor, or `NULL`.

## Examples

``` r
pollard_rho(8051)
#> [1] 97
pollard_pm1(21477639576571, B = 4000)
#> [1] 4410317
```
