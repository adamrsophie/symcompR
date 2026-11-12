# Partial derivative of tensor expressions

Port of SymPy's `PartialDerivative` (sympy/tensor/toperators.py): the
derivative of a tensor expression `expr` with respect to tensors or
symbols. The indices of a tensor variable enter with opposite valence;
repeated indices are contracted. Nested derivatives are flattened.
`pd_perform_derivative()` evaluates the derivative with Kronecker deltas
and metrics; `pd_expand()` applies the sum and product rules.

## Usage

``` r
PartialDerivative(expr, ...)

pd_expand(e)

pd_perform_derivative(e)
```

## Arguments

- expr:

  A tensor expression or a scalar expression.

- ...:

  Variables: tensors or symbols.

- e:

  A `PartialDerivative` or tensor expression.

## Value

A `PartialDerivative`, tensor expression or `sym`.

## Examples

``` r
L <- TensorIndexType("L")
i <- TensorIndex("i", L); j <- TensorIndex("j", L)
A <- TensorHead("A", list(L))
PartialDerivative(A(i), A(j))
#> PartialDerivative(A(i), A(j)) 
pd_perform_derivative(PartialDerivative(A(i), A(j)))
#> KD(i, -j) 
```
