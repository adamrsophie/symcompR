# Greatest common divisor and least common multiple

Works for integers, rationals and polynomials in any number of variables
(heuristic GCD with a primitive PRS fallback, as in SymPy). The result
is normalized to a positive leading coefficient with integer content.
`gcd_list()` and `lcm_list()` take a list. `gcdex()` returns `s`, `t`
and `h` with `s*f + t*g = h` for univariate polynomials; `invert()`
computes the inverse of `f` modulo `g`.

## Usage

``` r
gcd(f, g, ...)

gcd_list(seq, ...)

lcm_list(seq, ...)

gcdex(f, g, ...)

invert(f, g, ...)

sym_lcm(f, g, ...)
```

## Arguments

- f, g:

  Expressions or
  [`Poly()`](https://adamrsophie.github.io/symcompR/reference/Poly.md)
  objects.

- ...:

  Generators.

- seq:

  A list of expressions.

## Value

An expression.

## Examples

``` r
x <- Symbol("x"); y <- Symbol("y")
gcd(x^2 - 1, x^2 + 2*x + 1)
#> x + 1
gcd(x^2*y - y, x*y^2 + y^2)
#> x*y + y
sym_lcm(S(4), S(6))
#> 12
```
