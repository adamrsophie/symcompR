# Reduce a system of rational inequalities

`reduce_rational_inequalities(exprs, gen)` solves a system of rational
inequalities in one variable, as SymPy's function of the same name.
`exprs` is a list of lists (an Or of Ands) of relationals or
`list(expr, op)` pairs.

## Usage

``` r
reduce_rational_inequalities(exprs, gen, relational = TRUE)
```

## Arguments

- exprs:

  A list of lists of relationals.

- gen:

  The variable.

- relational:

  Return a relational instead of a set.

## Value

A relational or a set.

## Examples

``` r
x <- Symbol("x")
reduce_rational_inequalities(list(list(x^2 < 4)), x)
#> (-2 < x) & (x < 2)
```
