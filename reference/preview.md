# Render an expression with LaTeX

SymPy's `preview()`: typeset `expr` (an expression, or a LaTeX string)
in a standalone LaTeX document and compile it with
[`tools::texi2pdf()`](https://rdrr.io/r/tools/texi2dvi.html).
`output = "pdf"` keeps the PDF; `"png"` converts it with `pdftoppm` when
that program is available; `"tex"` only writes the source. With
`viewer = "file"` the result goes to `filename`; otherwise it is opened
with the system viewer in an interactive session. A LaTeX installation
is needed for `"pdf"` and `"png"`.

## Usage

``` r
preview(
  expr,
  output = "pdf",
  viewer = NULL,
  euler = TRUE,
  packages = character(),
  filename = NULL,
  preamble = NULL,
  extra_preamble = NULL,
  fontsize = NULL,
  ...
)
```

## Arguments

- expr:

  An expression or a LaTeX string.

- output:

  `"pdf"`, `"png"` or `"tex"`.

- viewer:

  `NULL`, or `"file"`.

- euler:

  Whether to use the Euler math fonts.

- packages:

  Extra LaTeX packages.

- filename:

  The output file for `viewer = "file"`.

- preamble:

  A full custom preamble.

- extra_preamble:

  Text added to the default preamble.

- fontsize:

  The font size, a number of points or a string.

- ...:

  Settings for
  [`latex()`](https://adamrsophie.github.io/symcompR/reference/latex.md).

## Value

The path of the output file, invisibly.

## Examples

``` r
x <- Symbol("x")
f <- preview(x^2 + 1, output = "tex", viewer = "file", filename = tempfile(fileext = ".tex"))
cat(readLines(f), sep = "\n")
#> \documentclass[varwidth,12pt]{standalone}
#> 
#> \usepackage{amsmath}
#> \usepackage{amsfonts}
#> \usepackage{euler}
#> \begin{document}
#> 
#> $\displaystyle x^{2} + 1$
#> 
#> \end{document}
```
