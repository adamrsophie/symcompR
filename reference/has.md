# Test whether an expression contains a subexpression

Test whether an expression contains a subexpression

## Usage

``` r
has(e, ...)
```

## Arguments

- e:

  A `sym` expression.

- ...:

  Subexpressions to look for. A symbol matches when it is free in `e`; a
  function such as `sin` or a
  [`Function()`](https://adamrsophie.github.io/symcompR/reference/Function.md)
  class matches any of its applications.

## Value

`TRUE` or `FALSE`.
