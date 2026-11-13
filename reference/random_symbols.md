# Random symbols in an expression

`random_symbols(expr)` lists the random symbols of `expr`, sorted by
name; `is_random(expr)` tells whether there are any; `pspace(expr)` is
the probability space they live on. These follow `sympy.stats`.

## Usage

``` r
random_symbols(expr)

is_random(expr)

pspace(expr)
```

## Arguments

- expr:

  An expression.

## Value

A list of random symbols, a logical, or a probability space.

## Examples

``` r
X <- Normal("X", 0, 1); Y <- Normal("Y", 0, 1)
random_symbols(2 * X + Y)
#> [[1]]
#> X
#> 
#> [[2]]
#> Y
#> 
is_random(X^2)
#> [1] TRUE
pspace(X)
#> SingleContinuousPSpace(X, NormalDistribution(0, 1))
```
