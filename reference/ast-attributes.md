# Fields of code generation nodes

`ast_attr(node, name)` returns a field or derived attribute of a codegen
node, as the Python attribute `node.name`: for example
`ast_attr(v, "symbol")`, `ast_attr(block, "left_hand_sides")` or
`ast_attr(float64, "eps")`. `ast_kwargs(node)` returns the fields as a
named list (`Token.kwargs()`), and `ast_func(node)` the constructor, so
`do.call(ast_func(node), sym_args(node))` rebuilds the node.
`ast_getitem(p, i)` indexes a
[`Pointer()`](https://adamrsophie.github.io/symcompR/reference/codegen-ast.md),
giving an
[`Element()`](https://adamrsophie.github.io/symcompR/reference/codegen-ast.md).

## Usage

``` r
ast_attr(node, name)

ast_kwargs(node, exclude = character(), apply = NULL)

ast_func(node)

ast_getitem(p, key)

attr_params(node, looking_for)
```

## Arguments

- node:

  A codegen node.

- name:

  The attribute name.

- exclude:

  Field names to leave out.

- apply:

  Optional function applied to every value.

- p:

  A
  [`Pointer()`](https://adamrsophie.github.io/symcompR/reference/codegen-ast.md).

- key:

  An index or a list of indices.

- looking_for:

  The attribute name.

## Value

The attribute value, a list, or a function.

## Examples

``` r
x <- Symbol("x")
v <- Variable(x, type = real)
ast_attr(v, "type")
#> real
ast_kwargs(v)
#> $symbol
#> x
#> 
#> $type
#> real
#> 
#> $value
#> NoneToken()
#> 
#> $attrs
#> ()
#> 
```
