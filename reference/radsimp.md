# Rationalize denominators and denest square roots

`radsimp()` removes square roots from denominators as SymPy's
`radsimp()`: at most `max_terms` radical terms are handled, and with
`symbolic = FALSE` denominators with symbols are left alone.
`rad_rationalize(num, den)` rationalizes `num/den` and returns the new
numerator and denominator. `collect_sqrt()` collects terms with common
square roots. `split_surds()` splits a sum of surds by their gcd.
`sqrtdenest()` rewrites `sqrt(a + b*sqrt(c))` as a sum of square roots
when possible. `rationalize()` is an alias of `radsimp()`.

## Usage

``` r
radsimp(e, symbolic = TRUE, max_terms = 4L)

rationalize(e, symbolic = TRUE, max_terms = 4L)

rad_rationalize(num, den)

split_surds(expr)

collect_sqrt(expr, evaluate = TRUE)

sqrtdenest(e)
```

## Arguments

- e, expr:

  An expression.

- symbolic:

  Also handle denominators with symbols.

- max_terms:

  The largest number of radical terms handled.

- num, den:

  Numerator and denominator.

- evaluate:

  For `collect_sqrt()`: `FALSE` returns the terms and the number of
  radical terms instead of the sum.

## Value

An expression, or a list.

## Examples

``` r
radsimp(1 / (1 + sqrt(S(2))))
#> -1 + sqrt(2)
sqrtdenest(sqrt(5 + 2 * sqrt(S(6))))
#> sqrt(2) + sqrt(3)
```
