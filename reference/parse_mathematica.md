# Parse a Mathematica expression

`parse_mathematica()` translates a string of Wolfram Mathematica code
into an expression, as SymPy's `parse_mathematica` does: standard input
form and FullForm are both read, known heads such as `Sin`, `Log`,
`Plus` or `Power` become their SymPy counterparts, and any other head
becomes an undefined function. `parse_mathematica_to_fullformlist()`
stops at the FullForm, returned as nested lists of strings.

## Usage

``` r
parse_mathematica(s)

parse_mathematica_to_fullformlist(s)
```

## Arguments

- s:

  A character string of Mathematica code.

## Value

For `parse_mathematica()` an expression; for
`parse_mathematica_to_fullformlist()` a string or a nested list.

## Examples

``` r
parse_mathematica("Sin[x]^2 Tan[y]")
#> sin(x)^2*tan(y)
parse_mathematica("Log[2, 8] + 1/2")
#> 7/2
parse_mathematica_to_fullformlist("x*(a + b)")
#> [[1]]
#> [1] "Times"
#> 
#> [[2]]
#> [1] "x"
#> 
#> [[3]]
#> [[3]][[1]]
#> [1] "Plus"
#> 
#> [[3]][[2]]
#> [1] "a"
#> 
#> [[3]][[3]]
#> [1] "b"
#> 
#> 
```
