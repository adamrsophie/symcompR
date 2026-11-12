# Cycles

`Cycle(...)` is SymPy's `Cycle`: `Cycle(1, 2)` (or `Cycle(c(1, 2))`)
maps 1 to 2 and 2 to 1. `cycle_mul(c, ...)` is SymPy's call `c(...)`:
the product of `c` and the cycle given by `...`, applying `c` first, so
`cycle_mul(Cycle(1, 2), 2, 3)` is SymPy's `Cycle(1, 2)(2, 3)`.
`cycle_as_list()` is SymPy's `Cycle.list()` and `cycle_size()` its
`size`. A `Cycle` converts to a permutation with
[`Permutation()`](https://adamrsophie.github.io/symcompR/reference/Permutation.md).

## Usage

``` r
Cycle(...)

cycle_mul(c, ...)

cycle_size(c)

cycle_as_list(c, size = NULL)
```

## Arguments

- ...:

  Integers of one cycle, a `Cycle` or a `Permutation`.

- c:

  A `Cycle`.

- size:

  Optional size of the list.

## Value

A `Cycle`, an integer vector or an integer.

## Examples

``` r
Cycle(1, 2)
#> (1 2)
cycle_mul(Cycle(1, 2), 2, 3)
#> (1 3 2)
cycle_as_list(Cycle(1, 2), 4)
#> [1] 0 2 1 3
```
