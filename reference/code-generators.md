# Code generators

Code generators write routines as source files, as the classes of
`sympy.utilities.codegen`: `CCodeGen()`, `C89CodeGen()`, `C99CodeGen()`,
`FCodeGen()`, `JuliaCodeGen()`, `OctaveCodeGen()` and `RustCodeGen()`.
Their methods (`routine`, `write`, `dump_c`, `dump_h`, `dump_f95`,
`dump_jl`, `dump_m`, `dump_rs`, `get_prototype`, ...) are read with
[`ast_attr()`](https://adamrsophie.github.io/symcompR/reference/ast-attributes.md);
the dump methods write to a connection. `get_code_generator()` returns
the generator for a language, and `codegen()` writes the source files of
named expressions.

## Usage

``` r
CCodeGen(
  project = "project",
  printer = NULL,
  preprocessor_statements = NULL,
  cse = FALSE
)

C89CodeGen(
  project = "project",
  printer = NULL,
  preprocessor_statements = NULL,
  cse = FALSE
)

C99CodeGen(
  project = "project",
  printer = NULL,
  preprocessor_statements = NULL,
  cse = FALSE
)

FCodeGen(project = "project", printer = NULL)

JuliaCodeGen(project = "project", printer = NULL)

OctaveCodeGen(project = "project", printer = NULL)

RustCodeGen(project = "project", printer = NULL)

get_code_generator(language, project = NULL, standard = NULL, printer = NULL)

codegen(
  name_expr,
  language = NULL,
  prefix = NULL,
  project = "project",
  to_files = FALSE,
  header = TRUE,
  empty = TRUE,
  argument_sequence = NULL,
  global_vars = NULL,
  standard = NULL,
  code_gen = NULL,
  printer = NULL
)
```

## Arguments

- project:

  The project name.

- printer:

  A code printer object, or `NULL` for the default.

- preprocessor_statements:

  C preprocessor lines.

- cse:

  Whether to eliminate common subexpressions.

- language:

  "C", "C89", "C99", "F95", "JULIA", "OCTAVE" or "RUST".

- standard:

  The C standard, for language "C".

- name_expr:

  A `list(name, expr)` pair, or a list of them.

- prefix:

  The file name prefix.

- to_files:

  Whether to write files instead of returning them.

- header, empty:

  Whether to write the header comment and empty lines.

- argument_sequence, global_vars:

  As in
  [`make_routine()`](https://adamrsophie.github.io/symcompR/reference/Routine.md).

- code_gen:

  A code generator, instead of `language`.

## Value

A code generator; `codegen()` returns a list of `(file name, contents)`
pairs.

## Examples

``` r
x <- Symbol("x"); y <- Symbol("y")
res <- codegen(list("f", x * y), "C99", header = FALSE, empty = FALSE)
cat(res[[1]][[2]])
#> #include "f.h"
#> #include <math.h>
#> double f(double x, double y) {
#>    double f_result;
#>    f_result = x*y;
#>    return f_result;
#> }
```
