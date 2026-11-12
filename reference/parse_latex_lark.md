# Parse LaTeX with the lark grammar

The lark backend of
[`parse_latex()`](https://adamrsophie.github.io/symcompR/reference/parse_latex.md),
with the readings of constants and functions given in `locals`.

## Usage

``` r
parse_latex_lark(s, locals = NULL)
```

## Arguments

- s:

  A character string holding LaTeX.

- locals:

  A named list: names are commands as written in `s` (the constants
  `\\pi`, `\\mathrm{e}`, `\\infty`, `\\imaginaryunit`, the trigonometric
  and hyperbolic functions, `\\exp`, `\\ln`, `\\log`, `\\lg`, or
  `\\operatorname{name}`); values are expressions for constants and
  functions (such as
  [`Function()`](https://adamrsophie.github.io/symcompR/reference/Function.md)
  results) for functions.

## Value

A `sym` expression, a matrix, or a `lark_tree` of the readings.

## Examples

``` r
parse_latex_lark("2\\pi r")
#> 2*pi*r
parse_latex_lark("2\\pi r", locals = list("\\pi" = Symbol("pi")))
#> 2*pi*r
```
