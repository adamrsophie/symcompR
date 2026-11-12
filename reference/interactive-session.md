# Interactive sessions

`init_printing()` chooses how symbolic expressions print at the console,
as SymPy's `init_printing()`: as two dimensional pretty printing
(`pretty_print = TRUE`, the default), as LaTeX (`use_latex = TRUE`), or
as plain strings (`pretty_print = FALSE`). `init_session()` also defines
the symbols `x`, `y`, `z`, `t`, the integer symbols `k`, `m`, `n` and
the functions `f`, `g`, `h` in `envir`, as SymPy's session does.
`interactive_traversal(expr)` walks an expression tree, asking at each
node which argument to enter, and returns the chosen subexpression.

## Usage

``` r
init_printing(pretty_print = TRUE, use_unicode = NULL, use_latex = FALSE)

init_session(
  pretty_print = TRUE,
  use_unicode = NULL,
  quiet = FALSE,
  envir = parent.frame()
)

interactive_traversal(expr)
```

## Arguments

- pretty_print:

  Pretty print expressions.

- use_unicode:

  Unicode characters in pretty printing (`NULL` keeps the current
  choice).

- use_latex:

  Print LaTeX.

- quiet:

  Do not print the banner.

- envir:

  The environment receiving the symbols.

- expr:

  An expression.

## Value

`init_printing()` and `init_session()` return `NULL` invisibly;
`interactive_traversal()` returns an expression.

## Examples

``` r
init_printing(pretty_print = FALSE)
interactive_traversal(Symbol("x") + 1)
#> x + 1
```
