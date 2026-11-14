# Escape a string for LaTeX

Escapes the characters that LaTeX treats specially, so the string prints
as plain text, as SymPy's `latex_escape()`.

## Usage

``` r
latex_escape(s)
```

## Arguments

- s:

  A character scalar.

## Value

The escaped string.

## Examples

``` r
latex_escape("50% of $x_1$")
#> [1] "50\\% of \\$x\\_1\\$"
```
