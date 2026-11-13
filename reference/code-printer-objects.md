# Code printer objects

Printer objects after SymPy's code printer classes. `C89CodePrinter()`,
`C99CodePrinter()` and `C11CodePrinter()` print C; `CXX98CodePrinter()`,
`CXX11CodePrinter()` and `CXX17CodePrinter()` print C++;
`FCodePrinter()` prints Fortran; `PythonCodePrinter()`,
`MpmathPrinter()`, `NumPyPrinter()` and `SciPyPrinter()` print Python;
`JuliaCodePrinter()`, `OctaveCodePrinter()`, `RustCodePrinter()` and
`JavascriptCodePrinter()` print those languages. `printer_doprint(p, e)`
prints an expression with a printer, as the `doprint` method does; the
printer keeps what the printed code needs (`headers`, `libraries`,
`macros`, `module_imports`), which
[`ast_attr()`](https://adamrsophie.github.io/symcompR/reference/ast-attributes.md)
reads.

## Usage

``` r
printer_doprint(p, e, assign_to = NULL)

C89CodePrinter(settings = NULL)

C99CodePrinter(settings = NULL)

C11CodePrinter(settings = NULL)

CXX98CodePrinter(settings = NULL)

CXX11CodePrinter(settings = NULL)

CXX17CodePrinter(settings = NULL)

FCodePrinter(settings = NULL)

JavascriptCodePrinter(settings = NULL)

JuliaCodePrinter(settings = NULL)

OctaveCodePrinter(settings = NULL)

PythonCodePrinter(settings = NULL)

MpmathPrinter(settings = NULL)

NumPyPrinter(settings = NULL)

SciPyPrinter(settings = NULL)

LambdaPrinter(settings = NULL)

RCodePrinter(settings = NULL)

RustCodePrinter(settings = NULL)
```

## Arguments

- p:

  A printer.

- e:

  An expression or codegen node.

- assign_to:

  Optional variable to assign to.

- settings:

  A named list of printer settings, such as `list(type_aliases = ...)`
  or `list(standard = 95)`.

## Value

A printer, or for `printer_doprint()` a character scalar.

## Examples

``` r
x <- Symbol("x")
p <- C99CodePrinter()
printer_doprint(p, sym_expm1(x))
#> [1] "expm1(x)"
printer_doprint(FCodePrinter(list(source_format = "free")), x^2)
#> [1] "x**2"
```
