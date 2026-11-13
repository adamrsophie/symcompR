# Replace tensors by arrays

Port of SymPy's `TensExpr.replace_with_arrays()`. `replacement_dict` is
a list of pairs `list(key, array)`, where `key` is a tensor such as
`A(i, j)` (its indices fix the meaning of the array axes) or a
`TensorIndexType` (its covariant metric). The result is an
[`Array()`](https://adamrsophie.github.io/symcompR/reference/Array.md)
with axes ordered as `indices` (default: the free indices of `expr`), or
an expression when no free index is left.

## Usage

``` r
replace_with_arrays(expr, replacement_dict, indices = NULL)
```

## Arguments

- expr:

  A tensor expression.

- replacement_dict:

  List of `list(key, array)` pairs.

- indices:

  Optional list of `TensorIndex` giving the axis order.

## Value

A `SymArray` or a `sym`.

## Examples

``` r
L <- TensorIndexType("L", dim = 2)
i <- TensorIndex("i", L); j <- TensorIndex("j", L)
A <- TensorHead("A", list(L))
replace_with_arrays(A(i) * A(j), list(list(A(i), Array(list(1, 2)))), list(i, j))
#> [[1, 2], [2, 4]]
```
