# Parse a string into an expression with SymPy's parser

`parse_expr()` is SymPy's `parse_expr`: the string is tokenized as
Python, rewritten by the token `transformations` (see
[standard_transformations](https://adamrsophie.github.io/symcompR/reference/parse_expr_transformations.md))
and evaluated with SymPy's names. Unknown names become symbols or
undefined functions, numbers become exact integers and floats, `x!` is a
factorial and `0.[3]` a repeating decimal. With `evaluate = FALSE` the
operators build unevaluated expressions.
[`sympify()`](https://adamrsophie.github.io/symcompR/reference/sympify.md)
is unchanged and does not use this parser.

## Usage

``` r
parse_expr(
  s,
  local_dict = NULL,
  transformations = standard_transformations,
  global_dict = NULL,
  evaluate = TRUE
)
```

## Arguments

- s:

  A character string.

- local_dict:

  A named list of names to use first.

- transformations:

  A list of token transformations, or `"all"` or `"implicit"` for
  SymPy's groups of them.

- global_dict:

  A named list to use instead of SymPy's names, or `NULL`.

- evaluate:

  If `FALSE`, arithmetic and the elementary functions are left
  unevaluated.

## Value

The value of the expression: usually a `sym`, but Python numbers,
strings, tuples and lists are returned as `pe_int`, `pe_float`,
`pe_str`, `pe_tuple` and `pe_list` objects.

## Examples

``` r
parse_expr("x**2 + 3*x - 1/2")
#> x^2 + 3*x - 1/2
parse_expr("2x sin x", transformations = "implicit")
#> 2*x*sin(x)
parse_expr("0.[3] + x!")
#> factorial(x) + 1/3
parse_expr("1 + 2", evaluate = FALSE)
#> 1 + 2
```
