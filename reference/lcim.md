# Least common multiple of periods

`lcim(numbers)` is the least common multiple of a list of rational
numbers, or of irrational numbers with the same irrational factor (such
as multiples of `pi`), and `NULL` otherwise, as SymPy's `lcim()`.

## Usage

``` r
lcim(numbers)
```

## Arguments

- numbers:

  A list of numbers.

## Value

A number or `NULL`.

## Examples

``` r
lcim(list(Rational(1, 2), 2, 3))
#> 6
lcim(list(S_Pi, S_Pi / 2, S_Pi / 4))
#> pi
```
