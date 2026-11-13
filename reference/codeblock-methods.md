# CodeBlock methods

`codeblock_topological_sort(assignments)` orders a list of
[`Assignment()`](https://adamrsophie.github.io/symcompR/reference/codegen-ast.md)
nodes so that every variable is assigned before it is used, keeping the
given order where it is free, as SymPy's `CodeBlock.topological_sort`.
`codeblock_cse(block)` runs common subexpression elimination over the
right hand sides of a block of assignments.

## Usage

``` r
codeblock_topological_sort(assignments)

codeblock_cse(
  block,
  symbols = NULL,
  optimizations = NULL,
  postprocess = NULL,
  order = "canonical"
)
```

## Arguments

- assignments:

  A list of
  [`Assignment()`](https://adamrsophie.github.io/symcompR/reference/codegen-ast.md)
  nodes.

- block:

  A
  [`CodeBlock()`](https://adamrsophie.github.io/symcompR/reference/codegen-ast.md)
  of assignments.

- symbols:

  Optional list of symbols for the new variables.

- optimizations, postprocess:

  Passed on to
  [`cse()`](https://adamrsophie.github.io/symcompR/reference/cse.md).

- order:

  Ordering of the subexpressions.

## Value

A
[`CodeBlock()`](https://adamrsophie.github.io/symcompR/reference/codegen-ast.md).

## Examples

``` r
x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
codeblock_topological_sort(list(Assignment(x, y + z), Assignment(z, 1), Assignment(y, 2)))
#> CodeBlock(
#> Assignment(z, 1),
#> Assignment(y, 2),
#> Assignment(x, y + z)
#> )
```
