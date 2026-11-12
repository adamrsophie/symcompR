# Continued fractions

`continued_fraction(a)` returns the terms of a rational or quadratic
irrational number. A rational gives a vector of terms; a quadratic
irrational gives a list whose last element is the repeating block, so
`continued_fraction(sqrt(S(7)))` is `list(2, c(1, 1, 1, 4))`. Wherever
terms are accepted, a list means that its last element repeats.
`continued_fraction_periodic(p, q, d, s)` does the same for
`(p + s * sqrt(d)) / q`. `continued_fraction_reduce(cf)` turns terms
back into the number. `continued_fraction_iterator(x, n)` returns up to
`n` terms of any real `x`, exactly for rationals and in double precision
otherwise (reliable for roughly 10 terms).
`continued_fraction_convergents(cf, n)` returns the convergents of the
terms, the first `n` of them when the fraction repeats.

## Usage

``` r
continued_fraction(a)

continued_fraction_periodic(p, q, d = 0, s = 1)

continued_fraction_reduce(cf)

continued_fraction_iterator(x, n = NULL)

continued_fraction_convergents(cf, n = NULL)
```

## Arguments

- a:

  A `sym` rational or quadratic irrational, or something
  [`as_sym()`](https://adamrsophie.github.io/symcompR/reference/as_sym.md)
  accepts.

- p, q, d, s:

  Integers with `d >= 0` and `q != 0`.

- cf:

  Terms, as returned above.

- x:

  A real number.

- n:

  Maximum number of terms or convergents.

## Value

A vector or list of terms, or a `sym` number, or a list of `sym`
convergents.

## Examples

``` r
continued_fraction(Rational(415, 93))
#> [1] 4 2 6 7
continued_fraction(sqrt(S(7)))
#> [[1]]
#> [1] 2
#> 
#> [[2]]
#> [1] 1 1 1 4
#> 
continued_fraction_periodic(3, 2, 7)
#> [[1]]
#> [1] 2
#> 
#> [[2]]
#> [1] 1 4 1 1
#> 
continued_fraction_reduce(list(1, 2))
#> sqrt(2)
continued_fraction_reduce(list(1, c(2)))
#> sqrt(2)
continued_fraction_iterator(sym_pi, 6)
#> [1]   3   7  15   1 292   1
continued_fraction_convergents(c(0, 2, 1, 2))
#> [[1]]
#> 0
#> 
#> [[2]]
#> 1/2
#> 
#> [[3]]
#> 1/3
#> 
#> [[4]]
#> 3/8
#> 
```
