# Low index subgroups

Finds the subgroups of index at most `N` of the `FpGroup` `G` that
contain the words `Y` (each a generator or inverse of a generator),
returning their coset tables, as SymPy's `low_index_subgroups()`.

## Usage

``` r
low_index_subgroups(G, N, Y = list())

first_in_class(C, Y = list())
```

## Arguments

- G:

  An `FpGroup`.

- N:

  The largest index.

- Y:

  A list of words.

- C:

  A coset table.

## Value

A list of coset tables.

## Examples

``` r
g <- free_group("x, y")
x <- g[[2]]; y <- g[[3]]
f <- FpGroup(g[[1]], list(x^2, y^3, (x * y)^4))
L <- low_index_subgroups(f, 4)
fp_str(lapply(L, function(C) C$table))
#> [1] "[[[0, 0, 0, 0]], [[0, 0, 1, 2], [1, 1, 2, 0], [3, 3, 0, 1], [2, 2, 3, 3]], [[0, 0, 1, 2], [2, 2, 2, 0], [1, 1, 0, 1]], [[1, 1, 0, 0], [0, 0, 1, 1]]]"
```
