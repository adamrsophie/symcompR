# Lambda code printers

Printers after `sympy.printing.lambdarepr`. `lambdarepr()` prints an
expression as the code
[`lambdify()`](https://adamrsophie.github.io/symcompR/reference/lambdify.md)
evaluates, with
[`LambdaPrinter()`](https://adamrsophie.github.io/symcompR/reference/code-printer-objects.md).
`NumExprPrinter()` prints code for the numexpr library: the expression
is wrapped in `numexpr.evaluate('...', truediv=True)`, and functions
numexpr lacks are an error. `IntervalPrinter()` prints numbers as mpmath
`mpi` intervals.

## Usage

``` r
lambdarepr(expr, ...)

NumExprPrinter(settings = NULL)

IntervalPrinter(settings = NULL)
```

## Arguments

- expr:

  An expression.

- ...:

  Printer settings.

- settings:

  A named list of printer settings.

## Value

A character scalar of code, or a printer.

## Examples

``` r
x <- Symbol("x")
y <- Symbol("y")
lambdarepr(Piecewise(list(1, x < 1), list(2, TRUE)))
#> [1] "((1) if (x < 1) else (2))"
printer_doprint(NumExprPrinter(), sin(x) + cos(y))
#> [1] "numexpr.evaluate('sin(x) + cos(y)', truediv=True)"
printer_doprint(IntervalPrinter(), x/2 + 1)
#> [1] "(mpi('1/2'))*x + mpi('1')"
```
