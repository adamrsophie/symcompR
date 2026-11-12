# Tensor expressions

Ports of SymPy's tensor expression functions. Products and sums of
tensors are built with `*`, `+` and `-` (or `TensMul()`, `TensAdd()`);
repeated indices, one contravariant and one covariant, are contracted
and renamed as dummies (`L_0`, `L_1`, ...). `canon_bp()` brings an
expression to canonical form with the Butler-Portugal algorithm;
`contract_metric(t, g)` contracts the metric `g`;
`contract_delta(t, delta)` contracts a Kronecker delta;
`substitute_indices(t, ...)` replaces free indices given as
`list(old, new)` pairs; `get_free_indices()` and `get_dummy_indices()`
list indices; `riemann_cyclic()` applies the cyclic identity of the
Riemann tensor; `tensor_mul()` multiplies its arguments.

## Usage

``` r
canon_bp(p)

contract_metric(t, g)

contract_delta(t, delta)

substitute_indices(t, ...)

get_free_indices(t)

get_dummy_indices(t)

TensMul(..., is_canon_bp = FALSE)

TensAdd(...)

tensor_mul(...)

riemann_cyclic_replace(t_r)

riemann_cyclic(t2)
```

## Arguments

- t, p:

  A tensor expression (other values are returned unchanged).

- g, delta:

  A metric or delta tensor head.

- ...:

  For `substitute_indices()`, pairs `list(old, new)`; for
  `tensor_mul()`, `TensMul()` and `TensAdd()`, the factors or terms.

- is_canon_bp:

  Whether the product is already canonical.

- t_r:

  A Riemann tensor expression.

- t2:

  A tensor expression.

## Value

A tensor expression or a `sym`; `get_free_indices()` and
`get_dummy_indices()` return lists of `TensorIndex`.

## Examples

``` r
L <- TensorIndexType("L", dummy_name = "L")
idx <- tensor_indices("i, j", L)
A <- TensorHead("A", list(L, L), TensorSymmetry.fully_symmetric(-2))
B <- TensorHead("B", list(L, L), TensorSymmetry.fully_symmetric(2))
canon_bp(A(idx$i, idx$j) * B(-idx$j, -idx$i))
#> 0
p <- TensorHead("p", list(L))
contract_metric(L$metric(idx$i, idx$j) * p(-idx$j), L$metric)
#> p(i) 
```
