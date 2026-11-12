# Integer factorization

`factorint(n)` returns the prime factorization of `n` as a named vector
of exponents whose names are the primes, in increasing order; `-1`
appears for negative `n` and `0` for zero. The method follows SymPy:
trial division, a perfect power test, Fermat's method, then Pollard p -
1 and rho alongside trial division, and finally Lenstra's elliptic curve
method. With no `limit`, numbers below `10^28` go to the compiled
Pollard rho of GMP, which gives the same unique result faster.

## Usage

``` r
factorint(
  n,
  limit = NULL,
  use_trial = TRUE,
  use_rho = TRUE,
  use_pm1 = TRUE,
  use_ecm = TRUE,
  verbose = FALSE,
  visual = NULL,
  multiple = FALSE
)

factorrat(
  rat,
  limit = NULL,
  use_trial = TRUE,
  use_rho = TRUE,
  use_pm1 = TRUE,
  verbose = FALSE,
  visual = NULL,
  multiple = FALSE
)

primefactors(n, limit = NULL, verbose = FALSE, ...)
```

## Arguments

- n:

  An integer, a named factor vector, a `sym` product, or a
  [`sym_factorial()`](https://adamrsophie.github.io/symcompR/reference/special-functions.md)
  expression.

- limit:

  Trial division bound, or `NULL`.

- use_trial, use_rho, use_pm1, use_ecm:

  Method switches.

- verbose:

  Report progress with
  [`message()`](https://rdrr.io/r/base/message.html).

- visual:

  `TRUE` for a `sym` product, `FALSE` for the vector, `NULL` for the
  default.

- multiple:

  List factors with repetition.

- rat:

  A rational number.

- ...:

  Passed to `factorint()`.

## Value

A named numeric vector, a vector of factors, or a `sym` product.

## Details

`limit` stops trial division at that bound (and disables the elliptic
curve method), so the last factor may be composite. The `use_*` switches
choose methods. `multiple = TRUE` lists the factors with repetition. A
named vector of factors (keys need not be prime) or an unevaluated
product is refactored, and returned as an unevaluated `sym` product
unless `visual = FALSE`; `visual = TRUE` gives that product for any
input. An unevaluated `sym_factorial(k)` is factored with Legendre's
formula.

`factorrat(r)` factors a rational, with negative exponents for the
denominator. `primefactors(n)` lists the distinct primes.

## Examples

``` r
factorint(2000)
#> 2 5 
#> 4 3 
factorint(-60)
#> -1  2  3  5 
#>  1  2  1  1 
factorint(2^20 + 1, multiple = TRUE)
#> [1]    17 61681
factorint(c("4" = 2, "12" = 3))
#> 2^10*3^3
factorint(100, visual = TRUE)
#> 2^2*5^2
factorrat(Rational(8, 9))
#>  2  3 
#>  3 -2 
primefactors(6936)
#> [1]  2  3 17
```
