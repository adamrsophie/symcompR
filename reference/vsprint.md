# Printing for classical mechanics

`vsprint()` is SymPy's `vsprint()`: the string form of `expr` in which a
function of time `q(t)` prints as `q` and its derivatives as `q'`,
`q''`. `vsstrrepr()` is the same with strings quoted, `vprint()` prints
it, `vlatex()` gives LaTeX with dots for time derivatives and
`vpprint()` pretty prints. `init_vprinting()` makes these the default
printers for the session, as SymPy's `init_vprinting()`.
`VectorLatexPrinter()` is the printer behind `vlatex()`, with SymPy's
`doprint()` method.

## Usage

``` r
mprint(expr, ...)

msprint(expr, ...)

mpprint(expr, use_unicode = NULL, ...)

mlatex(expr, ...)

mechanics_printing(...)

vsprint(expr, ...)

vsstrrepr(expr, ...)

vprint(expr, ...)

vlatex(expr, ...)

VectorLatexPrinter(...)

vpprint(expr, use_unicode = NULL, ...)

init_vprinting(...)
```

## Arguments

- expr:

  An expression, vector, dyadic or matrix.

- ...:

  Ignored printer settings, accepted for compatibility.

- use_unicode:

  Unicode or ASCII output for `vpprint()`.

## Value

A character scalar; `vpprint()` returns it invisibly after printing,
`vprint()` and `init_vprinting()` return `NULL` invisibly and
`VectorLatexPrinter()` returns a printer object.

## Examples

``` r
q <- dynamicsymbols("q")
N <- ReferenceFrame("N")
vsprint(dynamicsymbols("q", 1) * N$x)
#> [1] "q'*N.x"
vlatex(dynamicsymbols("q", 2))
#> [1] "\\ddot{q}"
VectorLatexPrinter()$doprint(dynamicsymbols("q", 1))
#> [1] "\\dot{q}"
```
