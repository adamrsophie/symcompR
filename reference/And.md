# Boolean operators

`And()`, `Or()`, `Not()`, `Xor()`, `Nand()`, `Nor()`, `Implies()`,
`Equivalent()` and `ITE()` build logical expressions. Symbols serve as
boolean variables, and `&`, `|` and `!` work on `sym` objects.

## Usage

``` r
And(...)

Or(...)

Not(a)

Xor(..., remove_true = TRUE)

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

- remove_true:

  For `Xor()`: with `FALSE` a true argument is kept as an argument
  instead of negating the result, as SymPy's `remove_true`.

## Value

A boolean `sym` expression.

## Examples

``` r
x <- Symbol("x"); y <- Symbol("y")
And(x, Or(y, Not(x)))
#> x & (y | ~x)
x & !x
#> x & ~x
```
