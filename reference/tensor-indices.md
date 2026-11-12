# Tensor index types, indices and symmetries

Ports of SymPy's `TensorIndexType`, `TensorIndex`, `tensor_indices()`
and `TensorSymmetry`. An index type has a name, a dummy name (the head
of generated dummy index names such as `L_0`), a dimension and a metric;
`L$metric`, `L$delta` and `L$epsilon` are its metric, Kronecker delta
and Levi-Civita tensor heads. A `TensorIndex` is contravariant (`is_up`)
or covariant; unary minus flips it. `TensorSymmetry(base, gens)` is a
monoterm slot symmetry given by a BSGS;
`TensorSymmetry.fully_symmetric(n)` (antisymmetric for `n < 0`),
`TensorSymmetry.direct_product(...)`, `TensorSymmetry.riemann()` and
`TensorSymmetry.no_symmetry(n)` build the usual ones.

## Usage

``` r
TensorIndexType(
  name,
  dummy_name = NULL,
  dim = NULL,
  eps_dim = NULL,
  metric_symmetry = 1,
  metric_name = "metric"
)

TensorIndex(name, tensor_index_type, is_up = TRUE)

tensor_indices(s, typ)

TensorSymmetry(...)

TensorSymmetry.fully_symmetric(rank)

TensorSymmetry.direct_product(...)

TensorSymmetry.riemann()

TensorSymmetry.no_symmetry(rank)
```

## Arguments

- name:

  Name (a string).

- dummy_name:

  Head of dummy index names (default: first letter).

- dim:

  Dimension (default: symbol `dim_<dummy_name>`).

- eps_dim:

  Dimension of the epsilon tensor (default `dim`).

- metric_symmetry:

  `1` symmetric, `-1` antisymmetric, `0` no symmetry, `NULL` no metric.

- metric_name:

  Name of the metric tensor.

- tensor_index_type, typ:

  A `TensorIndexType`.

- is_up:

  `TRUE` for a contravariant index.

- s:

  A string of index names, as in
  [`symbols()`](https://rdrr.io/r/graphics/symbols.html).

- ...:

  For `TensorSymmetry()`: a BSGS `list(base, gens)` or `base` and
  `gens`; for `TensorSymmetry.direct_product()`: integers, positive for
  symmetric and negative for antisymmetric groups of indices.

- rank:

  Number of indices; negative for antisymmetric.

## Value

A `TensorIndexType`, `TensorIndex`, list of indices or `TensorSymmetry`.

## Examples

``` r
L <- TensorIndexType("L", dummy_name = "L")
idx <- tensor_indices("i, j", L)
-idx$j
#> -j 
L$metric
#> metric(L,L) 
TensorSymmetry.fully_symmetric(-2)
#> TensorSymmetry((0,), ((0 1)(2 3),)) 
TensorSymmetry.riemann()
#> TensorSymmetry((0, 2), ((0 1)(4 5), (2 3)(4 5), (5)(0 2)(1 3))) 
```
