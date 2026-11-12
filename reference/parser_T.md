# Select parse_expr transformations by number

The `T[...]` shortcut of SymPy's parser: the transformations at the
given 0-based positions in SymPy's numbered list (0 `lambda_notation`, 1
`auto_symbol`, 2 `repeated_decimals`, 3 `auto_number`, 4
`factorial_notation`, 5 `implicit_multiplication_application`, 6
`convert_xor`, 7 `implicit_application`, 8 `implicit_multiplication`, 9
`convert_equals_signs`, 10 `function_exponentiation`, 11
`t_rationalize`). Negative positions count from the end, as in Python.

## Usage

``` r
parser_T(...)
```

## Arguments

- ...:

  Integer positions or vectors of positions.

## Value

A list of transformations.

## Examples

``` r
parse_expr("2x", transformations = parser_T(0:4, 8))
#> 2*x
```
