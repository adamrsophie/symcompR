# Translate a symbol name to LaTeX

Converts Greek letter names, other special names and trailing modifiers
such as `hat` or `dot` to LaTeX, as SymPy's
[`translate()`](https://adamrsophie.github.io/symcompR/reference/geometry-transformations.md).

## Usage

``` r
latex_translate(s)
```

## Arguments

- s:

  A character scalar.

## Value

The LaTeX string.

## Examples

``` r
latex_translate("alphahatdotprime")
#> [1] "{\\dot{\\hat{\\alpha}}}'"
```
