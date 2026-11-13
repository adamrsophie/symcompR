# Sum of products and product of sums forms

`SOPform(vars, minterms)` and `POSform(vars, minterms)` build minimal
disjunctive and conjunctive forms from the rows where the function is
true, with SymPy's Quine-McCluskey reduction. A row is an integer (most
significant variable first), a vector of bits, or a named list giving
the values of some variables (the others take every value).

## Usage

``` r
SOPform(vars, minterms, dontcares = integer())

POSform(vars, minterms, dontcares = integer())
```

## Arguments

- vars:

  A list of symbols.

- minterms:

  Rows where the function is true: an integer vector, or a list of
  integers, bit vectors and named lists.

- dontcares:

  Rows whose value does not matter, in the same forms.

## Value

A boolean expression.

## Examples

``` r
w <- Symbol("w"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
SOPform(list(w, x, y, z), c(1, 3, 7, 11, 15), dontcares = c(0, 2, 5))
#> (y & z) | (~w & ~x)
```
