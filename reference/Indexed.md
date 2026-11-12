# Indexed objects

`Idx(name, range)` is an integer index symbol, as SymPy's `Idx`; `range`
is an upper bound `n` (giving `0` to `n - 1`) or a pair of bounds.
`IndexedBase(name, shape)` returns a function building indexed objects
`A(i, j)` (SymPy's `A[i, j]`), and `Indexed(base, ...)` builds one
directly. `get_indices(e)` returns the free indices of an expression and
its symmetries; repeated indices in a product are summed (Einstein
convention), and the terms of a sum must have the same free indices.
`get_contraction_structure(e)` lists the summations: index tuples with
the terms summed over them, `NULL` for terms without summation, and
nested structures for subexpressions, as SymPy's dictionary.

## Usage

``` r
Idx(name, range = NULL)

IndexedBase(name, shape = NULL)

Indexed(base, ...)

get_indices(e)

get_contraction_structure(e)

format_contraction_structure(cs)
```

## Arguments

- name:

  The name.

- range:

  An upper bound or a pair of bounds.

- shape:

  Optional dimensions.

- base:

  An `IndexedBase` function or a name.

- ...:

  Indices.

- e:

  An expression.

- cs:

  A contraction structure.

## Value

`Idx()` and `Indexed()` return expressions; `IndexedBase()` a function;
`get_indices()` a list of the free indices and the symmetries;
`get_contraction_structure()` a list of `key`, `value` pairs.

## Examples

``` r
i <- Idx("i"); j <- Idx("j")
A <- IndexedBase("A"); x <- IndexedBase("x"); y <- IndexedBase("y")
A(i, j)
#> A[i, j]
get_indices(A(i, j) * x(j))$indices
#> [[1]]
#> i
#> 
format_contraction_structure(get_contraction_structure(x(i)*y(i) + A(j, j)))
#> [1] "{(i,): {x[i]*y[i]}, (j,): {A[j, j]}}"
```
