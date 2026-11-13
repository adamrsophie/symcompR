# Algebraic normal form

`to_anf(expr, deep)` rewrites a boolean expression as an exclusive or of
conjunctions of symbols (with `True` kept as a term), `is_anf(expr)`
tests that form, `anf_coeffs(truthvalues)` gives the coefficients of the
form from a truth table, and `ANFform(variables, truthvalues)` builds
the expression, as in `sympy.logic.boolalg`.
`bool_minterm(k, variables)`, `bool_maxterm()` and `bool_monomial()`
give the `k`th minterm, maxterm and monomial, where `k` is an integer or
a vector of bits.

## Usage

``` r
to_anf(expr, deep = TRUE)

is_anf(expr)

anf_coeffs(truthvalues)

ANFform(variables, truthvalues)

bool_minterm(k, variables)

bool_maxterm(k, variables)

bool_monomial(k, variables)
```

## Arguments

- expr:

  A boolean expression.

- deep:

  Convert the arguments too.

- truthvalues:

  A vector of 0 and 1 of length a power of two.

- variables:

  A list of symbols.

- k:

  An integer or a vector of bits.

## Value

A boolean expression, a logical value or an integer vector.

## Examples

``` r
x <- Symbol("x"); y <- Symbol("y")
to_anf(Or(x, y))
#> x ^ y ^ (x & y)
is_anf(Xor(x, And(x, y)))
#> [1] TRUE
anf_coeffs(c(1, 1, 1, 0))
#> [1] 1 0 0 1
ANFform(list(x, y), c(1, 1, 1, 0))
#> True ^ (x & y)
bool_minterm(3, list(x, y))
#> x & y
```
