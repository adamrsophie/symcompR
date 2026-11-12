# Pattern matching

`sym_match(e, pattern)` matches an expression against a pattern
containing
[`Wild()`](https://adamrsophie.github.io/symcompR/reference/Dummy.md)
symbols and returns a named list of bindings, or `NULL`. Sums and
products match up to commutativity, and a wild may take the identity.
`sym_replace(e, pattern, value)` rewrites every subexpression that
matches, bottom up; `value` is an expression in the same wild symbols or
a function of the named bindings.

## Usage

``` r
sym_match(e, pattern, old = FALSE)

sym_replace(e, pattern, value, exact = NULL)
```

## Arguments

- e:

  An expression.

- pattern:

  A pattern containing wild symbols.

- old:

  Use the old matching that skips the operation count check.

- value:

  The replacement.

- exact:

  `sym_replace()` replaces only when no wild matched zero; by default
  when the pattern has more than one wild.

## Value

A named list or `NULL`; `sym_replace()` returns an expression.

## Examples

``` r
x <- Symbol("x")
a <- Wild("a", exclude = list(x)); b <- Wild("b", exclude = list(x))
sym_match(3 * x + 2, a * x + b)
#> $a
#> 3
#> 
#> $b
#> 2
#> 
sym_replace(sin(x)^2 + 1, sin(Wild("u")), cos(Wild("u")))
#> cos(x)^2 + 1
```
