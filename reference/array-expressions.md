# Array expressions

Symbolic arrays and the operations on them, following
`sympy.tensor.array.expressions`. `ArraySymbol()` is an array of a given
shape, `ArrayElement()` one of its entries, `ZeroArray()` and
`OneArray()` arrays of zeros and ones. `ArrayTensorProduct()`,
`ArrayAdd()`, `PermuteDims()`, `ArrayDiagonal()`, `ArrayContraction()`,
`Reshape()`, `ArrayElementwiseApplyFunc()` and `ArraySum()` build
expressions; they keep their arguments as given unless
`canonicalize = TRUE`, as in SymPy. Axes are numbered from 0, as in
SymPy.

## Usage

``` r
ArraySymbol(name, shape)

ArrayElement(name, indices)

ZeroArray(...)

OneArray(...)

ArrayTensorProduct(..., canonicalize = FALSE)

ArrayAdd(..., canonicalize = FALSE)

PermuteDims(
  expr,
  permutation = NULL,
  index_order_old = NULL,
  index_order_new = NULL,
  canonicalize = FALSE
)

ArrayDiagonal(expr, ..., canonicalize = FALSE, allow_trivial_diags = FALSE)

ArrayElementwiseApplyFunc(fun, element)

ArrayContraction(expr, ..., canonicalize = FALSE)

Reshape(expr, shape)

ArraySum(f, ...)
```

## Arguments

- name:

  A name (character) or an expression.

- shape:

  A vector or list of dimensions.

- indices:

  A vector or list of indices.

- ...:

  Dimensions, arguments, or axis groups (integer vectors).

- canonicalize:

  Put the result in canonical form.

- expr:

  An array or matrix expression.

- permutation:

  A permutation: a `Permutation`, an array form or a list of cycles.

- index_order_old, index_order_new:

  Index orders defining the permutation, as an alternative to
  `permutation`.

- allow_trivial_diags:

  Allow diagonal groups of one axis.

- fun:

  A function or `Lambda` applied to every element.

- element:

  An array expression.

- f:

  An array expression to sum.

## Value

A `sym` array expression.

## Examples

``` r
A <- ArraySymbol("A", c(3, 3))
B <- ArraySymbol("B", c(3, 3))
ArrayContraction(ArrayTensorProduct(A, B), c(1, 2))
#> ArrayContraction(ArrayTensorProduct(A, B), (1, 2))
PermuteDims(A, c(1, 0))
#> PermuteDims(A, (0 1))
```
