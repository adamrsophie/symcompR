# Source files and modules from abstract syntax trees

`render_as_source_file(content)` prints `content` as C source preceded
by the `#include` lines it needs (SymPy's
`cutils.render_as_source_file`). `py_render_as_module(content)` prints
it as a Python module preceded by its imports (SymPy's
`pyutils.render_as_module`), and `f_render_as_module(definitions, name)`
prints Fortran definitions as a module (SymPy's
`futils.render_as_module`).

## Usage

``` r
render_as_source_file(content, Printer = C99CodePrinter, settings = NULL)

py_render_as_module(content, standard = "python3")

f_render_as_module(
  definitions,
  name,
  declarations = list(),
  printer_settings = NULL
)
```

## Arguments

- content:

  An expression or abstract syntax tree.

- Printer:

  A printer constructor such as
  [`C99CodePrinter()`](https://adamrsophie.github.io/symcompR/reference/code-printer-objects.md),
  or a printer class name.

- settings:

  `NULL`, or a list of printer settings.

- standard:

  The Python standard.

- definitions:

  A list of definitions, such as
  [`FunctionDefinition()`](https://adamrsophie.github.io/symcompR/reference/codegen-ast.md)
  nodes.

- name:

  The name of the module.

- declarations:

  A list of declarations.

- printer_settings:

  `NULL` (Fortran 2003 in free form), or a list of printer settings.

## Value

A character scalar.

## Examples

``` r
x <- Symbol("x", real = TRUE)
cat(render_as_source_file(Assignment(x, sin(x))), "\n")
#> #include <math.h>
#> 
#> x = sin(x); 
cat(py_render_as_module(Print(list(x), "%f")), "\n")
#> 
#> 
#> print("%f" % (x), end="") 
```
