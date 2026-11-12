# Query assumptions about an expression

`ask_fact()` returns `TRUE`, `FALSE` or `NA` (unknown) for a single
property, following SymPy's `expr.is_<property>` attributes. The helpers
[`is_positive()`](https://adamrsophie.github.io/symcompR/reference/assumption-queries.md),
[`is_real()`](https://adamrsophie.github.io/symcompR/reference/assumption-queries.md)
and friends are shortcuts.

## Usage

``` r
ask_fact(e, fact)
```

## Arguments

- e:

  A `sym` expression or a value convertible with
  [`as_sym()`](https://adamrsophie.github.io/symcompR/reference/as_sym.md).

- fact:

  One of the assumption names listed in
  [`Symbol()`](https://adamrsophie.github.io/symcompR/reference/Symbol.md).

## Value

`TRUE`, `FALSE` or `NA`.

## Examples

``` r
x <- Symbol("x", positive = TRUE)
ask_fact(x^2 + 1, "positive")
#> [1] TRUE
is_integer(Symbol("n", integer = TRUE) * 2)
#> [1] TRUE
```
