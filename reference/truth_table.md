# Truth table of a boolean expression

Truth table of a boolean expression

## Usage

``` r
truth_table(e, vars = NULL)
```

## Arguments

- e:

  A boolean expression.

- vars:

  Variables, by default all atoms in `e`.

## Value

A data frame with one logical column per variable and a `value` column.

## Examples

``` r
a <- Symbol("a"); b <- Symbol("b")
truth_table(Implies(a, b))
#>       a     b value
#> 1 FALSE FALSE  TRUE
#> 2 FALSE  TRUE  TRUE
#> 3  TRUE FALSE FALSE
#> 4  TRUE  TRUE  TRUE
```
