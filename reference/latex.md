# LaTeX output

Converts an expression to its LaTeX form, as SymPy's `latex()`.

## Usage

``` r
latex(e, ...)
```

## Arguments

- e:

  A `sym` expression, matrix, `Poly` or a plain R value.

- ...:

  Printer settings, as in SymPy: `mul_symbol`, `fold_short_frac`,
  `long_frac_ratio`, `fold_frac_powers`, `fold_func_brackets`,
  `root_notation`, `inv_trig_style`, `ln_notation`, `mode`, `itex`,
  `mat_str`, `mat_delim`, `symbol_names`, `order`, `full_prec`, `min`,
  `max`, `imaginary_unit`, `diff_operator`, `gothic_re_im`,
  `decimal_separator`, `parenthesize_super`, `mat_symbol_style`,
  `adjoint_style`, `disable_split_super_sub`, `perm_cyclic`,
  `open_delim` and `closed_delim`.

## Value

A character scalar.

## Examples

``` r
x <- Symbol("x")
latex(x^2 / 2 + sin(x))
#> [1] "\\frac{x^{2}}{2} + \\sin{\\left(x \\right)}"
latex(sin(x)^2)
#> [1] "\\sin^{2}{\\left(x \\right)}"
latex(2 * x / 3, fold_short_frac = TRUE)
#> [1] "2 x / 3"
```
