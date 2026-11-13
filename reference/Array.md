# N-dimensional symbolic arrays

`Array(data, shape)` stores expressions in row-major order.
`tensorproduct()`, `tensorcontraction()`, `permutedims()` and
`derive_by_array()` follow `sympy.tensor.array`. `arr_get()` reads one
entry using 1 based indices.

## Usage

``` r
Array(data, shape = NULL)

arr_get(a, idx)

tensorproduct(...)

tensorcontraction(a, axes)

permutedims(a, axes)

derive_by_array(a, x)

tensordiagonal(a, ...)
```

## Arguments

- data:

  A list of entries (or a nested list with `shape` omitted).

- shape:

  An integer vector.

- a, x:

  Arrays or expressions.

- idx:

  Indices.

- ...:

  Arrays; for `tensordiagonal()`, groups of 1 based axes to diagonalize.

- axes:

  A pair of 1 based axes to contract, or a permutation.

## Value

A `SymArray`, or an expression.

## Examples

``` r
x <- Symbol("x"); y <- Symbol("y")
derive_by_array(x^2 * y, list(x, y))
#> [2*x*y, x^2]
tensorcontraction(Array(list(1, 2, 3, 4), c(2, 2)), c(1, 2))
#> 5
```
