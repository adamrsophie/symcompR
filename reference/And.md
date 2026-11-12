# Boolean operators

`And()`, `Or()`, `Not()`, `Xor()`, `Nand()`, `Nor()`, `Implies()`,
`Equivalent()` and `ITE()` build logical expressions. Symbols serve as
boolean variables, and `&`, `|` and `!` work on `sym` objects.

## Usage

``` r
And(...)

Or(...)

Not(a)

Xor(...)

Nand(...)

Nor(...)

Implies(a, b)

Equivalent(...)

ITE(a, b, c)
```

## Arguments

- ...:

  Boolean expressions.

- a, b, c:

  Boolean expressions.

## Value

A boolean `sym` expression.

## Examples

``` r
x <- Symbol("x"); y <- Symbol("y")
And(x, Or(y, Not(x)))
#> x & (y | ~x)
x & !x
#> False
```
