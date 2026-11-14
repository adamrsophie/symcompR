# Assumptions shared by expressions

`common_assumptions(exprs, check)` returns the assumptions (among those
named in `check`, by default all of SymPy's) whose values are known and
equal for every expression, as SymPy's `common_assumptions()`.

## Usage

``` r
common_assumptions(exprs, check = NULL)
```

## Arguments

- exprs:

  A list of expressions.

- check:

  Optional names of assumptions to consider.

## Value

A named logical vector, sorted by name.

## Examples

``` r
common_assumptions(list(-4, 0, sqrt(Integer(2)), 2, sym_pi, oo))
#>   commutative     composite extended_real     imaginary           odd 
#>          TRUE         FALSE          TRUE         FALSE         FALSE 
common_assumptions(list(0, 1, 2), c("positive", "integer"))
#> integer 
#>    TRUE 
```
