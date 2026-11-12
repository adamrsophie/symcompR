# All the readings of an ambiguous LaTeX string

`parse_latex_lark_all()` returns the list of the expressions `s` may
mean under the lark grammar, with one element when `s` is not ambiguous.
When the same ambiguity occurs several times, the readings choose the
same way at every place where it occurs, as in SymPy's
`parse_latex_lark_all`.

## Usage

``` r
parse_latex_lark_all(s, locals = NULL)
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

A list of expressions.

## Examples

``` r
parse_latex_lark_all("\\sin xy")
#> [[1]]
#> sin(x*y)
#> 
#> [[2]]
#> y*sin(x)
#> 
```
