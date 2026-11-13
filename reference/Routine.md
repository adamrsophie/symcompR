# Routines for code generation

`Routine()` holds the arguments, results and variables of one routine,
as SymPy's `sympy.utilities.codegen.Routine`. `make_routine()` builds
one from an expression, as SymPy's `make_routine()`.

## Usage

``` r
Routine(name, arguments, results, local_vars, global_vars)

make_routine(
  name,
  expr,
  argument_sequence = NULL,
  global_vars = NULL,
  language = "F95"
)
```

## Arguments

- name:

  The routine name.

- arguments, results:

  Lists of arguments and results.

- local_vars, global_vars:

  Lists of local and global variables.

- expr:

  An expression, an equation, or a list of them.

- argument_sequence:

  The order of the arguments, or `NULL`.

- language:

  The language whose conventions to follow.

## Value

A routine (class `Routine`); its fields are read with
[`ast_attr()`](https://adamrsophie.github.io/symcompR/reference/ast-attributes.md).

## Examples

``` r
x <- Symbol("x"); y <- Symbol("y")
r <- make_routine("f", x^2 + y)
ast_attr(r, "arguments")
#> [[1]]
#> InputArgument(x) 
#> 
#> [[2]]
#> InputArgument(y) 
#> 
```
