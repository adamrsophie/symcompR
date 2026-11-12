# Parse a LaTeX string into an expression

`parse_latex()` converts a LaTeX math string to an expression with the
grammar and transformer of SymPy's lark LaTeX backend: an Earley parser
reads the whole grammar, including implicit products, fractions, roots,
integrals, sums, limits, derivatives in Leibniz notation, matrices and
relations. An input with several readings that give different
expressions returns all of them in a `lark_tree` with data `"_ambig"`,
printed as SymPy prints it.

## Usage

``` r
parse_latex(s, strict = FALSE, backend = "lark")
```

## Arguments

- s:

  A character string holding LaTeX.

- strict:

  Ignored; kept for the signature of SymPy's `parse_latex()`.

- backend:

  Only `"lark"` is available.

## Value

A `sym` expression, a matrix, or a `lark_tree` of the readings.

## Examples

``` r
parse_latex("\\frac{1 + x}{\\sqrt{y}}")
#> (x + 1)/sqrt(y)
parse_latex("\\int_0^1 x^2 dx")
#> Integral(x^2, (x, 0, 1))
parse_latex("x^2 + 2x + 1 = 0")
#> Eq(x^2 + 2*x + 1, 0)
```
