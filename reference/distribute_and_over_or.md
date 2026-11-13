# Distribution laws

`distribute_and_over_or(expr)` gives a conjunction of disjunctions and
`distribute_xor_over_and(expr)` an exclusive or of conjunctions, as the
SymPy functions of the same names
([`distribute_or_over_and()`](https://adamrsophie.github.io/symcompR/reference/piecewise_simplify.md)
gives a disjunction of conjunctions).

## Usage

``` r
distribute_and_over_or(expr)

distribute_xor_over_and(expr)
```

## Arguments

- expr:

  A boolean expression.

## Value

A boolean expression.

## Examples

``` r
A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C")
distribute_and_over_or(Or(And(A, B), C))
#> (A | C) & (B | C)
distribute_xor_over_and(And(A, Xor(B, C)))
#> (A & B) ^ (A & C)
```
