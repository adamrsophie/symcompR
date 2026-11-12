# Free symbols of an expression

Returns the symbols that appear unbound in `e`, sorted by name.
Variables of integration, summation and lambda arguments are bound and
not included.

## Usage

``` r
free_symbols(e)
```

## Arguments

- e:

  A `sym` expression.

## Value

A list of symbols.

## Examples

``` r
x <- Symbol("x"); y <- Symbol("y")
free_symbols(x^2 + y)
#> [[1]]
#> x
#> 
#> [[2]]
#> y
#> 
```
