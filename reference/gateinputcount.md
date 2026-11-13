# Number of gate inputs

`gateinputcount(expr)` counts the inputs of the logic gates (`And`,
`Or`, `Xor`, `Not` and `ITE`) in a boolean expression, as SymPy's
`gateinputcount()`.

## Usage

``` r
gateinputcount(expr)
```

## Arguments

- expr:

  A boolean expression.

## Value

An integer.

## Examples

``` r
a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c")
gateinputcount(Or(a, And(b, c)))
#> [1] 4
```
