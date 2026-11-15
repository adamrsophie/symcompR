# Content and primitive part

`primitive(f)` returns `list(content, primitive)`, and `terms_gcd()`
pulls out the greatest common monomial and the content of the
coefficients, as SymPy's `terms_gcd()`: an equation is handled side by
side, `deep = TRUE` works bottom up through the arguments,
`clear = FALSE` keeps a rational factor inside a sum when that leaves
integer coefficients, and `expand = FALSE` reads the expression without
expanding it.

## Usage

``` r
primitive(f, ...)

content(f, ...)

terms_gcd(f, ..., deep = FALSE, clear = TRUE, expand = TRUE)
```

## Arguments

- f:

  An expression.

- ...:

  Generators.

- deep:

  Apply to every argument first.

- clear:

  Pull a rational content out of a sum.

- expand:

  Expand before reading the terms.

## Value

A list or an expression.
