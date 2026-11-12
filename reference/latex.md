# LaTeX output

LaTeX output

## Usage

``` r
latex(e)
```

## Arguments

- e:

  A `sym` expression or matrix.

## Value

A character scalar.

## Examples

``` r
x <- Symbol("x")
latex(x^2 / 2 + sin(x))
#> [1] "\\frac{x^{2}}{2} + \\sin{\\left(x \\right)}"
```
