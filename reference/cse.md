# Common subexpression elimination

`cse(exprs)` finds repeated subexpressions and returns a list with
`replacements` (a list of `list(symbol, expression)` pairs) and
`reduced` (the rewritten expressions), like SymPy's `cse()`.

## Usage

``` r
cse(exprs, prefix = "x")
```

## Arguments

- exprs:

  An expression or list of expressions.

- prefix:

  Prefix of the new symbols.

## Value

A list.

## Examples

``` r
x <- Symbol("x"); y <- Symbol("y")
cse(list((x + y)^2 + sin(x + y), (x + y) * 3))
#> $replacements
#> $replacements[[1]]
#> $replacements[[1]][[1]]
#> x0
#> 
#> $replacements[[1]][[2]]
#> x + y
#> 
#> 
#> 
#> $reduced
#> $reduced[[1]]
#> x0^2 + sin(x0)
#> 
#> $reduced[[2]]
#> 3*x + 3*y
#> 
#> 
```
