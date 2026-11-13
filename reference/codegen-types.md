# Codegen type helpers

`type_from_expr(e)` deduces a codegen type from an expression, as
SymPy's `Type.from_expr`. `type_cast_check(t, value)` casts `value` to
type `t` and checks that the value did not change by more than the
tolerance and that it lies in the type's range;
`type_cast_nocheck(t, value)` casts without the checks.
`variable_deduced(symbol, value)` builds a
[`Variable()`](https://adamrsophie.github.io/symcompR/reference/codegen-ast.md)
whose type is deduced from the symbol or the value, and
`as_Declaration(v, ...)` turns a variable into a
[`Declaration()`](https://adamrsophie.github.io/symcompR/reference/codegen-ast.md)
with some fields replaced.

## Usage

``` r
type_from_expr(e)

type_cast_nocheck(t, value)

type_cast_check(t, value, rtol = NULL, atol = 0)

variable_deduced(symbol, value = NULL, attrs = list(), cast_check = TRUE)

pointer_deduced(symbol, value = NULL, attrs = list(), cast_check = TRUE)

as_Declaration(v, ...)

from_FunctionDefinition(func_def)

from_FunctionPrototype(func_proto, body)
```

## Arguments

- e, value:

  An expression or number.

- t:

  A codegen type such as `real` or `float32`.

- rtol, atol:

  Relative and absolute tolerances.

- symbol:

  A symbol.

- attrs:

  Attributes of the variable.

- cast_check:

  Check the cast of `value`.

- v:

  A
  [`Variable()`](https://adamrsophie.github.io/symcompR/reference/codegen-ast.md).

- ...:

  Fields to replace.

- func_def:

  A
  [`FunctionDefinition()`](https://adamrsophie.github.io/symcompR/reference/codegen-ast.md).

- func_proto:

  A
  [`FunctionPrototype()`](https://adamrsophie.github.io/symcompR/reference/codegen-ast.md).

- body:

  The body statements.

## Value

A type, a number, or a node.

## Examples

``` r
type_from_expr(Symbol("n", integer = TRUE))
#> integer
type_cast_check(int8, 127)
#> 127
variable_deduced(Symbol("x", real = TRUE))
#> Variable(x, type=real)
```
