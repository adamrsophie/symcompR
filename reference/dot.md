# Matrix utilities

`dot()` and `cross()` for vectors, `norm_matrix()` (Frobenius or vector
2-norm), `applyfunc()` elementwise, `is_symmetric()`, `is_square()`,
`is_zero_matrix()`, `hadamard_product()`, `kronecker_product()`, `vec()`
and `pinv()` (Moore-Penrose pseudo inverse for full rank matrices).
Given an integer rather than a matrix, `is_square()` tests for a perfect
square, as `sympy.ntheory.primetest.is_square` does.

## Usage

``` r
dot(a, b, doit = FALSE)

cross(a, b, doit = FALSE)

norm_matrix(M)

applyfunc(M, f)

is_square(M)

is_symmetric(M)

is_zero_matrix(M)

hadamard_product(a, b)

kronecker_product(a, b)

vec(M)

pinv(M)
```

## Arguments

- a, b:

  Vectors as matrices, or
  [`CoordSys3D()`](https://adamrsophie.github.io/symcompR/reference/CoordSys3D.md)
  vectors, dyadics and the
  [`Del()`](https://adamrsophie.github.io/symcompR/reference/Del.md)
  operator, for which the methods of `sympy.vector` apply.

- doit:

  For the del operator: evaluate the derivatives.

- M:

  A matrix.

- f:

  A function applied to each entry.

## Value

An expression or matrix.
