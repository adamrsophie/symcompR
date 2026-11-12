# Sum of products and product of sums forms

`SOPform(vars, minterms)` and `POSform(vars, minterms)` build minimal
disjunctive and conjunctive forms from the rows (as integers, most
significant variable first) where the function is true.

## Usage

``` r
SOPform(vars, minterms, dontcares = integer())

POSform(vars, minterms, dontcares = integer())
```

## Arguments

- vars:

  A list of symbols.

- minterms:

  Integer row indices where the function is true.

- dontcares:

  Rows whose value does not matter.

## Value

A boolean expression.

## Examples

``` r
w <- Symbol("w"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
SOPform(list(w, x, y, z), c(1, 3, 7, 11, 15), dontcares = c(0, 2, 5))
#> (y & z) | (z & ~w)
```
