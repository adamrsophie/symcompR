# Tables of expressions

`TableForm(data, ...)` arranges expressions in a table, as SymPy's
`TableForm`. `data` is a list of rows, an R matrix, or an explicit
symbolic matrix; short rows are padded.
[`format()`](https://rdrr.io/r/base/format.html) and
[`print()`](https://rdrr.io/r/base/print.html) give the text table,
[`latex()`](https://adamrsophie.github.io/symcompR/reference/latex.md) a
`tabular` environment, and
[`as.matrix()`](https://rdrr.io/r/base/matrix.html) the entries as a
symbolic matrix.

## Usage

``` r
TableForm(
  data,
  headings = NULL,
  alignments = "l",
  formats = NULL,
  pad = NULL,
  wipe_zeros = TRUE
)

# S3 method for class 'TableForm'
format(x, ...)

# S3 method for class 'TableForm'
print(x, ...)

# S3 method for class 'TableForm'
as.matrix(x, ...)
```

## Arguments

- data:

  A list of rows, a matrix, or a symbolic matrix.

- headings:

  `NULL`, `"automatic"`, or a list of row and column headings, either of
  which may be `NULL` or `"automatic"`.

- alignments:

  `"l"`, `"r"`, `"c"` (or `"left"`, `"right"`, `"center"`, `"<"`, `">"`,
  `"^"`) for every column, or one per column, with an optional first one
  for the row headings.

- formats:

  `NULL`, or one entry per column for LaTeX output: `NULL`, a
  [`sprintf()`](https://rdrr.io/r/base/sprintf.html) format, or a
  function of the entry and its zero based row and column that returns a
  string (or `NULL` for the default).

- pad:

  The entry used for missing values and short rows.

- wipe_zeros:

  Print zero entries as blanks.

- x:

  A `TableForm`.

- ...:

  Unused.

## Value

A `TableForm` object; [`format()`](https://rdrr.io/r/base/format.html)
returns a string.

## Examples

``` r
x <- Symbol("x")
TableForm(list(list(5, 7), list(4, 2), list(10, 3)), headings = "automatic")
#>   | 1  2
#> --------
#> 1 | 5  7
#> 2 | 4  2
#> 3 | 10 3
latex(TableForm(list(list(x^2, 1/x)), headings = list(NULL, c("a", "b"))))
#> [1] "\\begin{tabular}{l l}\na & b \\\\\n\\hline\n$x^{2}$ & $\\frac{1}{x}$ \\\\\n\\end{tabular}"
```
