# Printing helpers

`print_latex()`, `print_ccode()`, `print_fcode()`, `print_rcode()` and
`print_jscode()` print the corresponding strings. `julia_code()`,
`cxxcode()`, `mathematica_code()` and `maple_code()` are further code
printers. `print_tree()` shows the expression tree.

## Usage

``` r
julia_code(e, assign_to = NULL, ...)

cxxcode(e, assign_to = NULL, ...)

maple_code(e, assign_to = NULL, ...)

mathematica_code(e)

print_latex(e)

print_ccode(e, assign_to = NULL)

print_fcode(e, assign_to = NULL)

print_rcode(e, assign_to = NULL)

print_jscode(e, assign_to = NULL)

print_tree(e)
```

## Arguments

- e:

  An expression.

- assign_to:

  Optional variable name.

- ...:

  Printer settings, as in SymPy (for example `standard`).

## Value

The string, invisibly for the `print_` functions.

## Examples

``` r
x <- Symbol("x")
mathematica_code(sin(x)^2 + exp(x))
#> [1] "Exp[x] + Sin[x]^2"
print_tree(x^2 + 1)
#> Add
#>   |-Integer: 1
#>   +-Pow
#>     |-Symbol: x
#>     +-Integer: 2
```
