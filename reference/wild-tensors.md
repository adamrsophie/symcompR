# Wild tensors and tensor pattern matching

Ports of SymPy's `WildTensorIndex`, `WildTensorHead` and `WildTensor`. A
`WildTensorIndex` matches any index of its index type, of the same
variance unless `ignore_updown = TRUE`. A `WildTensorHead` is called
like a tensor head; the result, a `WildTensor`, matches any tensor
expression whose free indices match its indices (in any order when
`unordered_indices = TRUE`), and with no indices it matches any tensor
expression.
[`sym_match()`](https://adamrsophie.github.io/symcompR/reference/sym_match.md)
matches a tensor expression against a pattern and
[`sym_replace()`](https://adamrsophie.github.io/symcompR/reference/sym_match.md)
replaces the matches, as `expr.match()` and `expr.replace()` do in
SymPy; a matched wild tensor head stands for the matched expression with
its free indices renamed (SymPy's `_WildTensExpr`).

## Usage

``` r
WildTensorIndex(name, tensor_index_type, is_up = TRUE, ignore_updown = FALSE)

WildTensorHead(
  name,
  index_types = NULL,
  symmetry = NULL,
  comm = 0,
  unordered_indices = FALSE
)
```

## Arguments

- name:

  The name of the index or head.

- tensor_index_type:

  A `TensorIndexType`.

- is_up:

  Whether the index is contravariant.

- ignore_updown:

  Whether the index matches indices of either variance.

- index_types, symmetry, comm:

  As for
  [`TensorHead()`](https://adamrsophie.github.io/symcompR/reference/TensorHead.md);
  the index types of a wild tensor are those of its indices, and wild
  heads have no symmetry.

- unordered_indices:

  Whether the indices of the wild tensor match in any order.

## Value

A wild index or a wild tensor head.

## Examples

``` r
R3 <- TensorIndexType("R3", dim = 3)
idx <- tensor_indices("p q", R3)
p <- idx$p; q <- idx$q
a <- WildTensorIndex("a", R3, ignore_updown = TRUE)
K <- TensorHead("K", list(R3)); V <- TensorHead("V", list(R3))
W <- WildTensorHead("W", unordered_indices = TRUE)
sym_match(K(p) * V(q), W(a, p))
#> $a
#> q 
#> 
#> $`W(R3,R3)`
#> _WildTensExpr(K(p)*V(q)) 
#> 
sym_replace(K(p) * V(-p), W(a) * V(-a), 1)
#> 1
sym_replace(K(p) * V(q), W(q) * K(p), W(p) * W(q))
#> V(p)*V(q) 
```
