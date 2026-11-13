# The symbol of a literal

`literal_symbol(literal)` is the symbol of a literal (a symbol or its
negation), as SymPy's `literal_symbol()`; `TRUE` and `FALSE` are
returned unchanged and anything else is an error.

## Usage

``` r
literal_symbol(literal)
```

## Arguments

- literal:

  A literal.

## Value

A symbol or a logical value.

## Examples

``` r
A <- Symbol("A")
literal_symbol(Not(A))
#> A
```
