# Data types and arguments of code generation routines

`DataType()` holds the names of a data type in the target languages, and
`get_default_datatype()` chooses one for an expression, as in
`sympy.utilities.codegen`. `InputArgument()`, `OutputArgument()`,
`InOutArgument()` and `Result()` describe the arguments and results of a
routine made by
[`make_routine()`](https://adamrsophie.github.io/symcompR/reference/Routine.md);
their fields are read with
[`ast_attr()`](https://adamrsophie.github.io/symcompR/reference/ast-attributes.md).

## Usage

``` r
DataType(cname, fname, pyname, jlname, octname, rsname)

get_default_datatype(expr, complex_allowed = NULL)

InputArgument(name, datatype = NULL, dimensions = NULL, precision = NULL)

OutputArgument(
  name,
  result_var,
  expr,
  datatype = NULL,
  dimensions = NULL,
  precision = NULL
)

InOutArgument(
  name,
  result_var,
  expr,
  datatype = NULL,
  dimensions = NULL,
  precision = NULL
)

Result(
  expr,
  name = NULL,
  result_var = NULL,
  datatype = NULL,
  dimensions = NULL,
  precision = NULL
)
```

## Arguments

- cname, fname, pyname, jlname, octname, rsname:

  The type names in C, Fortran, Python, Julia, Octave and Rust.

- expr:

  An expression.

- complex_allowed:

  Whether a complex type is allowed. `NULL` reads the option
  `symcompR.complex_allowed` (default `FALSE`), which plays the part of
  SymPy's module flag `COMPLEX_ALLOWED`.

- name:

  A symbol (or, for `Result()`, a name).

- datatype:

  A `DataType()`.

- dimensions:

  A list of `list(lower, upper)` bounds.

- precision:

  The precision, or `NULL`.

- result_var:

  The variable or element the result is stored in.

## Value

An object of the corresponding class.

## Examples

``` r
x <- Symbol("x")
get_default_datatype(x)
#> <DataType double> 
InputArgument(x)
#> InputArgument(x) 
```
