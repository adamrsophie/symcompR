# Factoring with noncommutative symbols

`factor_nc(e)` factors an expression that contains noncommutative
symbols, as SymPy's `factor_nc()`: common commutative factors and common
noncommutative prefixes and suffixes come out, and the remaining sum is
factored with the noncommutative symbols kept in order.

## Usage

``` r
factor_nc(e)
```

## Arguments

- e:

  An expression.

## Value

An expression.

## Examples

``` r
x <- Symbol("x")
A <- Symbol("A", commutative = FALSE); B <- Symbol("B", commutative = FALSE)
factor_nc(x^2 + 2*A*x + A^2)
#> (x + A)^2
factor_nc(expand((x + A)*(x + B)))
#> (x + A)*(x + B)
factor_nc(A*B + A*Symbol("C", commutative = FALSE))
#> A*(B + C)
```
