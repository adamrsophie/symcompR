# Polynomial division

`div(f, g)` returns `list(q, r)` with `f = q*g + r`, using the first
generator as main variable over the rationals. `quo()` and `rem()`
return one part; `exquo()` fails unless the division is exact. `pdiv()`,
`prem()` and `pquo()` perform pseudo division. `reduced()` divides by a
list of polynomials.

## Usage

``` r
div(f, g, ...)

quo(f, g, ...)

rem(f, g, ...)

exquo(f, g, ...)

pdiv(f, g, ...)

prem(f, g, ...)

pquo(f, g, ...)

reduced(f, G, ..., order = "lex")
```

## Arguments

- f, g:

  Expressions or
  [`Poly()`](https://adamrsophie.github.io/symcompR/reference/Poly.md)
  objects.

- ...:

  Generators.

- G:

  A list of polynomials.

- order:

  Monomial order: `"lex"`, `"grlex"` or `"grevlex"`.

## Value

Expressions, or `Poly` objects when inputs are `Poly`.

## Examples

``` r
x <- Symbol("x")
div(x^3 + 2*x + 1, x - 1)
#> $q
#> x^2 + x + 3
#> 
#> $r
#> 4
#> 
rem(x^4 + 1, x^2 + 1)
#> 2
```
