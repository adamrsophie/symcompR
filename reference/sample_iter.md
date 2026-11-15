# Draw samples as expressions

SymPy's `sample_iter()`: `numsamples` values of `expr`, each an exact
expression, drawn subject to `condition`. R has no iterators, so the
values come back as a list.

## Usage

``` r
sample_iter(expr, condition = NULL, numsamples = 1, seed = NULL, ...)
```

## Arguments

- expr:

  An expression with random variables.

- condition:

  An optional condition on them.

- numsamples:

  The number of samples.

- seed:

  An optional random seed.

- ...:

  Ignored (SymPy's `library` and `size`).

## Value

A list of expressions.

## Examples

``` r
X <- Die("X", 6)
sample_iter(X, numsamples = 3, seed = 1)
#> [[1]]
#> 3
#> 
#> [[2]]
#> 4
#> 
#> [[3]]
#> 5
#> 
```
