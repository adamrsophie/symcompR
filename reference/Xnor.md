# Xnor and Exclusive

`Xnor(...)` is the negation of
[`Xor()`](https://adamrsophie.github.io/symcompR/reference/And.md) and
`Exclusive(...)` is true when at most one argument is true, as in SymPy,
where both evaluate to `Not` and `And` expressions.

## Usage

``` r
Xnor(...)

Exclusive(...)
```

## Arguments

- ...:

  Boolean expressions.

## Value

A boolean expression.

## Examples

``` r
A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C")
Xnor(A, B)
#> ~(A ^ B)
Exclusive(A, B, C)
#> ~(A & B) & ~(A & C) & ~(B & C)
```
