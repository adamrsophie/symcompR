# Arrays and tensors

``` r

library(symcompR)
x <- Symbol("x"); y <- Symbol("y")
```

This vignette covers `sympy.tensor.array` (N-dimensional arrays),
`sympy.tensor.tensor` (abstract index notation), `sympy.tensor.indexed`,
and the canonicalization algorithms of `sympy.combinatorics.tensor_can`.

## Arrays

`Array(data, shape)` is an N-dimensional array of expressions, filled in
row major order as in SymPy.
[`arr_get()`](https://adamrsophie.github.io/symcompR/reference/Array.md)
reads an entry, with indices starting at 1 as usual in R. The helper
below collects the entries of a one or two dimensional array into an R
character vector or matrix for display.

``` r

entries <- function(a) {
  s <- shape(a)
  if (length(s) == 1) return(vapply(seq_len(s), function(i) sstr(arr_get(a, i)), ""))
  outer(seq_len(s[1]), seq_len(s[2]),
        Vectorize(function(i, j) sstr(arr_get(a, c(i, j)))))
}
M <- Array(list(1, 2, 3, 4), c(2, 2))
shape(M)
#> [1] 2 2
entries(M)
#>      [,1] [,2]
#> [1,] "1"  "2" 
#> [2,] "3"  "4"
arr_get(M, c(1, 2))
#> 2
entries(tensorproduct(Array(list(1, 2)), Array(list(x, y))))
#>      [,1]  [,2] 
#> [1,] "x"   "y"  
#> [2,] "2*x" "2*y"
tensorcontraction(M, c(1, 2))
#> 5
entries(permutedims(M, c(2, 1)))
#>      [,1] [,2]
#> [1,] "1"  "3" 
#> [2,] "2"  "4"
```

[`derive_by_array()`](https://adamrsophie.github.io/symcompR/reference/Array.md)
differentiates with respect to an array of variables, so applying it
twice gives the Hessian.

``` r

g <- derive_by_array(x^2 * y, list(x, y))
entries(g)
#> [1] "2*x*y" "x^2"
entries(derive_by_array(g, list(x, y)))
#>      [,1]  [,2] 
#> [1,] "2*y" "2*x"
#> [2,] "2*x" "0"
```

[`tensordiagonal()`](https://adamrsophie.github.io/symcompR/reference/Array.md)
takes the diagonal over pairs of axes, and
[`NDimArray()`](https://adamrsophie.github.io/symcompR/reference/array-aliases.md)
is another name for
[`Array()`](https://adamrsophie.github.io/symcompR/reference/Array.md).

``` r

C4 <- Array(as.list(1:16), c(2, 2, 2, 2))
entries(tensordiagonal(Array(list(1, 2, 3, 4), c(2, 2)), c(1, 2)))
#> [1] "1" "4"
shape(tensordiagonal(C4, c(1, 2), c(3, 4)))
#> [1] 2 2
shape(NDimArray(list(1, 2, 3, 4, 5, 6), c(2, 3)))
#> [1] 2 3
```

### Transforms over subset lattices

[`walsh_hadamard_transform()`](https://adamrsophie.github.io/symcompR/reference/convolution.md)
and the Mobius transforms of `sympy.discrete.transforms` work on
sequences indexed by subsets:

``` r

walsh_hadamard_transform(list(1, 2, 3, 4))
#> [[1]]
#> 10
#> 
#> [[2]]
#> -2
#> 
#> [[3]]
#> -4
#> 
#> [[4]]
#> 0
mobius_transform(list(1, 2, 3, 4))
#> [[1]]
#> 1
#> 
#> [[2]]
#> 3
#> 
#> [[3]]
#> 4
#> 
#> [[4]]
#> 10
inverse_mobius_transform(list(1, 3, 4, 10))
#> [[1]]
#> 1
#> 
#> [[2]]
#> 2
#> 
#> [[3]]
#> 3
#> 
#> [[4]]
#> 4
```

### Contraction structure of indexed expressions

[`get_contraction_structure()`](https://adamrsophie.github.io/symcompR/reference/Indexed.md)
finds the summation indices of an expression built from `Indexed`
objects, and
[`format_contraction_structure()`](https://adamrsophie.github.io/symcompR/reference/Indexed.md)
makes the result readable.

``` r

i <- Idx("i"); j <- Idx("j")
A <- IndexedBase("A"); B <- IndexedBase("B")
format_contraction_structure(get_contraction_structure(A(i, j) * B(j)))
#> [1] "{(j,): {A[i, j]*B[j]}}"
format_contraction_structure(get_contraction_structure(A(i, i)))
#> [1] "{(i,): {A[i, i]}}"
Indexed(A, i)
#> A[i]
```

## Tensor index notation

[`TensorIndexType()`](https://adamrsophie.github.io/symcompR/reference/tensor-indices.md)
declares a kind of index (with a metric), and
[`TensorHead()`](https://adamrsophie.github.io/symcompR/reference/TensorHead.md)
declares a tensor with given symmetries. Expressions follow the Einstein
summation convention: a repeated upper and lower index is summed, and a
minus sign marks a lower index.
[`canon_bp()`](https://adamrsophie.github.io/symcompR/reference/tensor-expressions.md)
brings an expression to canonical form, which detects symmetry
cancellations.

``` r

L <- TensorIndexType("L", dummy_name = "L")
idx <- tensor_indices("i, j", L)
i <- idx$i; j <- idx$j
A <- TensorHead("A", list(L, L), TensorSymmetry.fully_symmetric(-2))
B <- TensorHead("B", list(L, L), TensorSymmetry.fully_symmetric(2))
p <- TensorHead("p", list(L))
canon_bp(A(i, j) * B(-j, -i))
#> 0
canon_bp(B(i, j) - B(j, i))
#> 0
contract_metric(L$metric(i, j) * p(-j), L$metric)
#> p(i)
```

The contraction of an antisymmetric tensor with a symmetric one
vanishes, and contracting with the metric raises the index of `p`.

### Building tensor expressions

[`TensorIndex()`](https://adamrsophie.github.io/symcompR/reference/tensor-indices.md)
creates a single index, and
[`tensor_heads()`](https://adamrsophie.github.io/symcompR/reference/TensorHead.md)
several tensor heads at once.
[`TensMul()`](https://adamrsophie.github.io/symcompR/reference/tensor-expressions.md),
[`TensAdd()`](https://adamrsophie.github.io/symcompR/reference/tensor-expressions.md),
and
[`tensor_mul()`](https://adamrsophie.github.io/symcompR/reference/tensor-expressions.md)
build products and sums explicitly.
[`get_free_indices()`](https://adamrsophie.github.io/symcompR/reference/tensor-expressions.md)
and
[`get_dummy_indices()`](https://adamrsophie.github.io/symcompR/reference/tensor-expressions.md)
report the free and contracted indices.

``` r

k <- TensorIndex("k", L)
k
#> k
hp <- tensor_heads("q, r", list(L))
q <- hp[[1]]; r <- hp[[2]]
e <- tensor_mul(q(i), r(-i), p(j))
e
#> q(L_0)*r(-L_0)*p(j)
get_free_indices(e)
#> [[1]]
#> j
get_dummy_indices(e)
#> [[1]]
#> L_0 
#> 
#> [[2]]
#> -L_0
TensAdd(q(i), r(i))
#> q(i) + r(i)
doit(TensMul(2, q(i)))
#> 2*q(i)
substitute_indices(q(i) * r(j), list(i, k))
#> q(k)*r(j)
```

### Symmetries

[`TensorSymmetry()`](https://adamrsophie.github.io/symcompR/reference/tensor-indices.md)
describes the symmetry of a tensor through a base and strong generating
set.
[`TensorSymmetry.no_symmetry()`](https://adamrsophie.github.io/symcompR/reference/tensor-indices.md),
[`TensorSymmetry.direct_product()`](https://adamrsophie.github.io/symcompR/reference/tensor-indices.md),
and
[`TensorSymmetry.riemann()`](https://adamrsophie.github.io/symcompR/reference/tensor-indices.md)
cover the common cases.

``` r

TensorSymmetry.no_symmetry(2)
#> TensorSymmetry((), ((3),))
TensorSymmetry.direct_product(2, -2)
#> TensorSymmetry((0, 2), ((5)(0 1), (2 3)(4 5)))
Rm <- TensorHead("R", rep(list(L), 4), TensorSymmetry.riemann())
TensorSymmetry(get_symmetric_group_sgs(2))
#> TensorSymmetry((0,), ((3)(0 1),))
```

The Riemann tensor satisfies the cyclic identity R\_{ijkl} + R\_{iklj} +
R\_{iljk} = 0, which
[`riemann_cyclic_replace()`](https://adamrsophie.github.io/symcompR/reference/tensor-expressions.md)
and
[`riemann_cyclic()`](https://adamrsophie.github.io/symcompR/reference/tensor-expressions.md)
use to rewrite it:

``` r

idx4 <- tensor_indices("m0, m1, m2, m3", L)
riemann_cyclic_replace(Rm(idx4$m0, idx4$m2, idx4$m1, idx4$m3))
#> -(1/3)*R(m0, m3, m2, m1) + (1/3)*R(m0, m1, m2, m3) + (2/3)*R(m0, m2, m1, m3)
riemann_cyclic(Rm(i, j, k, idx4$m0) * (Rm(-i, -j, -k, -idx4$m0) -
                                     2 * Rm(-i, -k, -j, -idx4$m0)))
#> 0
```

### Kronecker delta

[`contract_delta()`](https://adamrsophie.github.io/symcompR/reference/tensor-expressions.md)
contracts the Kronecker delta of an index type:

``` r

Ln <- TensorIndexType("Ln", dim = 3)
d <- tensor_indices("a, b", Ln)
Ln$delta
#> KD(Ln,Ln)
contract_delta(Ln$delta(d$a, -d$b) * Ln$delta(d$b, -d$a), Ln$delta)
#> 3
```

### Partial derivatives

[`PartialDerivative()`](https://adamrsophie.github.io/symcompR/reference/PartialDerivative.md)
is the derivative of a tensor expression with respect to a tensor.
[`pd_perform_derivative()`](https://adamrsophie.github.io/symcompR/reference/PartialDerivative.md)
evaluates it and
[`pd_expand()`](https://adamrsophie.github.io/symcompR/reference/PartialDerivative.md)
distributes it over sums.

``` r

Af <- TensorHead("Af", list(L))
pd <- PartialDerivative(Af(i), Af(j))
get_free_indices(pd)
#> [[1]]
#> i 
#> 
#> [[2]]
#> -j
pd_perform_derivative(pd)
#> KD(i, -j)
pd_expand(PartialDerivative(Af(i) + q(i), Af(j)))
#> PartialDerivative(Af(i), Af(j)) + PartialDerivative(q(i), Af(j))
```

### Replacing tensors with arrays

[`replace_with_arrays()`](https://adamrsophie.github.io/symcompR/reference/replace_with_arrays.md)
evaluates an expression once the components of each tensor and of the
metric are known. Here H^{i}{}\_{j} has the given components and the
metric is \mathrm{diag}(1, -1).

``` r

H <- TensorHead("H", list(L, L))
repl <- list(list(H(i, -j), Array(list(list(1, 2), list(3, 4)))),
             list(L, Array(list(list(1, 0), list(0, -1)))))
entries(replace_with_arrays(H(i, j), repl))
#>      [,1] [,2]
#> [1,] "1"  "-2"
#> [2,] "3"  "-4"
entries(replace_with_arrays(H(i, j), repl, list(i, -j)))
#>      [,1] [,2]
#> [1,] "1"  "2" 
#> [2,] "3"  "4"
```

`TensorManager` holds the commutation rules between tensor groups:

``` r

class(TensorManager)
#> [1] "environment"
```

## Canonicalization of tensor expressions

The `tensor_can` functions implement the Butler-Portugal algorithm that
[`canon_bp()`](https://adamrsophie.github.io/symcompR/reference/tensor-expressions.md)
uses. A tensor expression is encoded as a permutation of index slots,
and the symmetries as base and strong generating sets.
[`get_symmetric_group_sgs()`](https://adamrsophie.github.io/symcompR/reference/tensor-bsgs.md)
gives the BSGS of a symmetric or antisymmetric tensor, and
`riemann_bsgs` that of the Riemann tensor.

``` r

sgs2 <- get_symmetric_group_sgs(2)
sgs2
#> [[1]]
#> [1] 0
#> 
#> [[2]]
#> [[2]][[1]]
#> (3)(0 1)
get_symmetric_group_sgs(3, antisym = TRUE)
#> [[1]]
#> [1] 0 1
#> 
#> [[2]]
#> [[2]][[1]]
#> (0 1)(3 4)
#> 
#> [[2]][[2]]
#> (1 2)(3 4)
riemann_bsgs
#> [[1]]
#> [1] 0 2
#> 
#> [[2]]
#> [[2]][[1]]
#> (0 1)(4 5)
#> 
#> [[2]][[2]]
#> (2 3)(4 5)
#> 
#> [[2]][[3]]
#> (5)(0 2)(1 3)
```

`canonicalize(g, dummies, msym, ...)` returns the canonical form of the
permutation `g`. The example is A^{d_0 d_1} B\_{d_1 d_0} written with
slot order (2, 1, 3, 0) and two signature slots, for symmetric A and B.

``` r

P <- function(...) Permutation(c(...))
canonicalize(P(2, 1, 3, 0, 4, 5), 0:3, 0, list(sgs2[[1]], sgs2[[2]], 1, 0),
             list(sgs2[[1]], sgs2[[2]], 1, 0))
#> [1] 0 2 1 3 4 5
canonical_free(integer(), sgs2[[2]], P(1, 0, 2, 3), 2)
#> [1] 1 0 2 3
```

The helpers build the BSGS of products and of the dummy index
symmetries:

``` r

tensor_gens(sgs2[[1]], sgs2[[2]], list(0, 1), 0)
#> [[1]]
#> [1] 6
#> 
#> [[2]]
#> integer(0)
#> 
#> [[3]]
#> [[3]][[1]]
#> (5)
gens_products(list(sgs2[[1]], sgs2[[2]], list(0, 1), 0))
#> [[1]]
#> [1] 6
#> 
#> [[2]]
#> integer(0)
#> 
#> [[3]]
#> [[3]][[1]]
#> (5)
bsgs_direct_product(sgs2[[1]], sgs2[[2]], sgs2[[1]], sgs2[[2]])
#> [[1]]
#> [1] 0 2
#> 
#> [[2]]
#> [[2]][[1]]
#> (5)(0 1)
#> 
#> [[2]][[2]]
#> (5)(2 3)
dummy_sgs(0:3, 0, 4)
#> [[1]]
#> [1] 1 0 2 3 4 5
#> 
#> [[2]]
#> [1] 0 1 3 2 4 5
#> 
#> [[3]]
#> [1] 2 3 0 1 4 5
perm_af_direct_product(list(c(1, 0, 2, 3), c(0, 1, 3, 2)), list(c(1, 0)), 0)
#> [[1]]
#> [1] 1 0 2 3 4 5
#> 
#> [[2]]
#> [1] 0 1 3 2 4 5
#> 
#> [[3]]
#> [1] 0 1 2 3 5 4
get_minimal_bsgs(c(1, 0), list(P(0, 2, 1, 3, 4), P(1, 0, 2, 3, 4)))
#> NULL
```

[`get_transversals()`](https://adamrsophie.github.io/symcompR/reference/tensor-bsgs.md)
and
[`transversal2coset()`](https://adamrsophie.github.io/symcompR/reference/transversal2coset.md)
list the coset representatives along a base, and
[`double_coset_can_rep()`](https://adamrsophie.github.io/symcompR/reference/double_coset_can_rep.md)
finds the canonical representative of a double coset.

``` r

sgs3 <- get_symmetric_group_sgs(3)
tr <- get_transversals(sgs3[[1]], sgs3[[2]])
length(tr)
#> [1] 2
transversal2coset(5, sgs3[[1]], tr)
#> [[1]]
#> [[1]][[1]]
#> [1] 0 1 2 3 4
#> 
#> [[1]][[2]]
#> [1] 1 0 2 3 4
#> 
#> [[1]][[3]]
#> [1] 2 0 1 3 4
#> 
#> 
#> [[2]]
#> [[2]][[1]]
#> [1] 0 1 2 3 4
#> 
#> [[2]][[2]]
#> [1] 0 2 1 3 4
dg <- dummy_sgs(0:5, 0, 6)
double_coset_can_rep(list(0:5), list(0), c(0, 2), dg,
                     get_transversals(c(0, 2), dg), P(4, 1, 3, 0, 5, 2, 6, 7))
#> [1] 0 2 1 4 3 5 6 7
```

## Array expressions

The array expressions of `sympy.tensor.array.expressions` describe
operations on arrays whose entries are not given.
[`ArraySymbol()`](https://adamrsophie.github.io/symcompR/reference/array-expressions.md)
is an array of a given shape,
[`ArrayTensorProduct()`](https://adamrsophie.github.io/symcompR/reference/array-expressions.md)
the tensor product,
[`ArrayContraction()`](https://adamrsophie.github.io/symcompR/reference/array-expressions.md)
contracts pairs of axes,
[`ArrayDiagonal()`](https://adamrsophie.github.io/symcompR/reference/array-expressions.md)
takes a diagonal over axes,
[`PermuteDims()`](https://adamrsophie.github.io/symcompR/reference/array-expressions.md)
permutes the axes, and
[`ArrayAdd()`](https://adamrsophie.github.io/symcompR/reference/array-expressions.md)
adds arrays. Axes are numbered from zero, as in SymPy.
[`get_shape()`](https://adamrsophie.github.io/symcompR/reference/get_shape.md)
returns the shape as an R list, shown here as a
[`Tuple()`](https://adamrsophie.github.io/symcompR/reference/Tuple.md).

``` r

AA <- ArraySymbol("A", c(3, 3)); BB <- ArraySymbol("B", c(3, 3))
tp <- ArrayTensorProduct(AA, BB)
tp
#> ArrayTensorProduct(A, B)
do.call(Tuple, get_shape(tp))
#> (3, 3, 3, 3)
get_rank(tp)
#> [1] 4
get_ndim(tp)
#> [1] 4
ArrayContraction(tp, c(1, 2))
#> ArrayContraction(ArrayTensorProduct(A, B), (1, 2))
ArrayDiagonal(tp, c(1, 2))
#> ArrayDiagonal(ArrayTensorProduct(A, B), (1, 2))
PermuteDims(AA, c(1, 0))
#> PermuteDims(A, (0 1))
ArrayAdd(AA, BB)
#> ArrayAdd(A, B)
nest_permutation(PermuteDims(tp, c(1, 0, 3, 2)))
#> PermuteDims(ArrayTensorProduct(A, B), (0 1)(2 3))
```

[`ArrayElement()`](https://adamrsophie.github.io/symcompR/reference/array-expressions.md)
is a single entry,
[`ZeroArray()`](https://adamrsophie.github.io/symcompR/reference/array-expressions.md)
and
[`OneArray()`](https://adamrsophie.github.io/symcompR/reference/array-expressions.md)
are the arrays of zeros and ones,
[`Reshape()`](https://adamrsophie.github.io/symcompR/reference/array-expressions.md)
changes the shape,
[`ArrayElementwiseApplyFunc()`](https://adamrsophie.github.io/symcompR/reference/array-expressions.md)
applies a function to every entry, and
[`ArraySum()`](https://adamrsophie.github.io/symcompR/reference/array-expressions.md)
is a sum of array expressions.

``` r

ArrayElement(AA, c(0, 1))
#> A[0, 1]
ZeroArray(2, 3)
#> ZeroArray(2, 3)
OneArray(2, 3)
#> OneArray(2, 3)
as_explicit(OneArray(2, 2))
#> [[1, 1], [1, 1]]
Reshape(AA, 9)
#> Reshape(A, (9,))
do.call(Tuple, get_shape(Reshape(AA, 9)))
#> (9,)
ArrayElementwiseApplyFunc(sin, AA)
#> ArrayElementwiseApplyFunc(Lambda(_d, sin(_d)), A)
i <- Symbol("i"); j <- Symbol("j"); d <- Symbol("d")
ArraySum(i * AA, list(i, 1, 3))
#> Sum(i*A, (i, 1, 3))
doit(ArraySum(i * AA, list(i, 1, 3)))
#> ArrayTensorProduct(6, A)
```

[`ArrayComprehension()`](https://adamrsophie.github.io/symcompR/reference/ArrayComprehension.md)
builds an array from an expression in its index variables;
[`ArrayComprehensionMap()`](https://adamrsophie.github.io/symcompR/reference/ArrayComprehension.md)
takes an R function of the index instead.

``` r

ac <- ArrayComprehension(i * j, list(i, 1, 3), list(j, 2, 4))
do.call(Tuple, get_shape(ac))
#> (3, 3)
doit(ac)
#> [[2, 3, 4], [4, 6, 8], [6, 9, 12]]
acm <- ArrayComprehensionMap(function(i) i + 1, list(i, 1, 5))
acm
#> ArrayComprehensionMap((i, 1, 5))
doit(acm)
#> [2, 3, 4, 5, 6]
```

Matrix expressions convert to array expressions and back, and so do
indexed sums.
[`collect_tensor_products()`](https://adamrsophie.github.io/symcompR/reference/collect_tensor_products.md)
collects sums of tensor products with a common factor.

``` r

k <- Symbol("k")
MP <- MatrixSymbol("P", k, k); MQ <- MatrixSymbol("Q", k, k)
MX <- MatrixSymbol("X", k, k)
cg <- convert_matrix_to_array(MP * t(MQ))
cg
#> ArrayContraction(ArrayTensorProduct(P, Q), (1, 3))
convert_array_to_matrix(cg)
#> P*Q.T
convert_matrix_to_array(Trace(MP * MQ))
#> ArrayContraction(ArrayTensorProduct(P, Q), (0, 3), (1, 2))
convert_indexed_to_array(Sum(MatrixElement(MP, i, d) * MatrixElement(MQ, d, j),
                             list(d, 0, k - 1)), list(i, j))
#> ArrayContraction(ArrayTensorProduct(P, Q), (1, 2))
convert_array_to_indexed(ArrayContraction(ArrayTensorProduct(MP, MQ), c(1, 2)), list(i, j))
#> Sum(P[i, _d0]*Q[_d0, j], (_d0, 0, k - 1))
collect_tensor_products(ArrayAdd(ArrayTensorProduct(MP, MQ), ArrayTensorProduct(MP, MX)))
#> ArrayTensorProduct(P, Q + X)
```

[`array_derive()`](https://adamrsophie.github.io/symcompR/reference/array_derive.md)
differentiates an array expression with respect to a matrix or array
symbol, and
[`ArrayDerivative()`](https://adamrsophie.github.io/symcompR/reference/array_derive.md)
is the unevaluated derivative of an array.

``` r

array_derive(convert_matrix_to_array(MP * MX), MX)
#> PermuteDims(ArrayContraction(ArrayTensorProduct(I, I, P), (1, 5)), (2 3))
ArrayDerivative(MX, MX)
#> Derivative(X, X)
ArrayDerivative(MX, MX, evaluate = TRUE)
#> PermuteDims(ArrayTensorProduct(I, I), (3)(1 2))
tryCatch(ArrayDerivative(ArraySymbol("A", list(2, 2)), ArraySymbol("A", list(2, 2))),
         error = conditionMessage)
#> [1] "Can't calculate derivative wrt A."
ArrayDerivative(Array(list(x^2, x * y)), x)
#> Derivative([x^2, x*y], x)
doit(ArrayDerivative(Array(list(x^2, x * y)), x))
#> [2*x, y]
```
