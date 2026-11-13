# Substitution for mechanics expressions

`msubs()` substitutes into expressions, matrices, vectors and dyadics as
SymPy's `msubs()`: a structural replacement that leaves derivatives
whole unless they appear as keys. With `smart = TRUE` tangents are
rewritten and denominators that vanish are simplified first.

## Usage

``` r
msubs(expr, ..., smart = FALSE)
```

## Arguments

- expr:

  An expression, matrix, vector or dyadic.

- ...:

  Substitution dictionaries: lists of `list(old, new)` pairs.

- smart:

  Use the smart substitution.

## Value

The substituted object.

## Examples

``` r
q <- dynamicsymbols("q")
msubs(q^2 + dynamicsymbols("q", 1), list(list(q, 2)))
#> Derivative(q(t), t) + 4
```
