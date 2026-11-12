# Tensor heads

Port of SymPy's `TensorHead` and `tensor_heads()`. A tensor head has a
name, the index types of its slots, a monoterm symmetry (default: none)
and a commutation group `comm` (see
[TensorManager](https://adamrsophie.github.io/symcompR/reference/TensorManager.md)).
The head is a function: `A(i, -j)` builds the tensor with indices `i`
(contravariant) and `-j` (covariant); indices may also be strings such
as `"-j"`. `A$name`, `A$index_types`, `A$symmetry`, `A$rank` and
`A$comm` give its properties.

## Usage

``` r
TensorHead(name, index_types, symmetry = NULL, comm = 0)

tensor_heads(s, index_types, symmetry = NULL, comm = 0)
```

## Arguments

- name:

  Name; for `tensor_heads()` a string of names as in
  [`symbols()`](https://rdrr.io/r/graphics/symbols.html).

- index_types:

  List of `TensorIndexType` objects.

- symmetry:

  A `TensorSymmetry`, or `NULL` for no symmetry.

- comm:

  Commutation group label.

- s:

  A string of names.

## Value

A `TensorHead` function, or a list of them.

## Examples

``` r
L <- TensorIndexType("L", dummy_name = "L")
i <- TensorIndex("i", L); j <- TensorIndex("j", L)
A <- TensorHead("A", list(L, L), TensorSymmetry.fully_symmetric(2))
A(i, -i)
#> A(L_0, -L_0) 
A
#> A(L,L) 
```
