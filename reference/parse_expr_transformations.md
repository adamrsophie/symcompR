# Token transformations of parse_expr

The transformations of SymPy's
[`parse_expr()`](https://adamrsophie.github.io/symcompR/reference/parse_expr.md),
each a function of the token list, the local namespace and the global
namespace that returns the transformed token list.
`standard_transformations` is the default list: `lambda_notation`,
`auto_symbol`, `repeated_decimals`, `auto_number` and
`factorial_notation`. `convert_xor` reads `^` as a power,
`implicit_multiplication` and `implicit_application` add the omitted `*`
and parentheses, `split_symbols` splits names like `xy`,
`function_exponentiation` reads `sin**2(x)`, `t_rationalize` (SymPy's
`rationalize`) makes floats exact and `convert_equals_signs` reads `=`
as [`Eq()`](https://adamrsophie.github.io/symcompR/reference/Eq.md).
`split_symbols_custom(predicate)` builds a symbol splitter.

## Usage

``` r
function_exponentiation(tokens, local_dict, global_dict)

split_symbols_custom(predicate)

split_symbols(tokens, local_dict, global_dict)

implicit_multiplication(tokens, local_dict, global_dict)

implicit_application(tokens, local_dict, global_dict)

implicit_multiplication_application(tokens, local_dict, global_dict)

auto_symbol(tokens, local_dict, global_dict)

lambda_notation(tokens, local_dict, global_dict)

factorial_notation(tokens, local_dict, global_dict)

convert_xor(tokens, local_dict, global_dict)

repeated_decimals(tokens, local_dict, global_dict)

auto_number(tokens, local_dict, global_dict)

t_rationalize(tokens, local_dict, global_dict)

convert_equals_signs(tokens, local_dict, global_dict)

standard_transformations
```

## Arguments

- tokens:

  A list of tokens, each `c(type, value)`.

- local_dict:

  An environment holding the local names.

- global_dict:

  A named list of global names, or `NULL` for SymPy's.

- predicate:

  A function of a name returning `TRUE` when it is split.

## Value

The transformed token list; `split_symbols_custom()` returns a
transformation.

## Examples

``` r
parse_expr("2x", transformations = c(standard_transformations, implicit_multiplication))
#> 2*x
parse_expr("x^2", transformations = c(standard_transformations, convert_xor))
#> x^2
```
