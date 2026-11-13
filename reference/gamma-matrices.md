# Gamma matrices

`GammaMatrix` is the tensor head of the Dirac gamma matrices, with one
index of type `LorentzIndex` (dimension 4, dummy name `L`), as in
`sympy.physics.hep.gamma_matrices`. `kahane_simplify()` removes the
contracted gamma matrices of a product by Kahane's algorithm,
`gamma_trace()` is the trace of a product of gamma matrices,
`simplify_gamma_expression()` simplifies the gamma matrices of a tensor
product, `simplify_gpgp()` contracts products of gamma matrices with
momenta and `extract_type_tens()` splits a product into the tensors of
one head and the rest.

## Usage

``` r
LorentzIndex

GammaMatrix(...)

extract_type_tens(expression, component)

simplify_gamma_expression(expression)

simplify_gpgp(ex, sort = TRUE)

gamma_trace(t)

kahane_simplify(expression)
```

## Arguments

- ...:

  Indices of the gamma matrix.

- expression, t, ex:

  Tensor expressions.

- component:

  A tensor head.

- sort:

  Whether to sort the components first.

## Value

A tensor expression, a matrix or, for `extract_type_tens()`, a list of
two expressions.

## Examples

``` r
i0 <- TensorIndex("i0", LorentzIndex); i1 <- TensorIndex("i1", LorentzIndex)
kahane_simplify(GammaMatrix(i0) * GammaMatrix(i1) * GammaMatrix(-i0))
#> -2*GammaMatrix(i1) 
gamma_trace(GammaMatrix(i0) * GammaMatrix(i1))
#> 4*metric(i0, i1) 
```
