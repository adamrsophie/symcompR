# Matrices of physics

`msigma(i)` is the Pauli matrix with index 1, 2 or 3;
`mgamma(mu, lower)` the Dirac gamma matrix with index 0, 1, 2, 3 or 5
(upper index by default; `lower = TRUE` lowers it with the Minkowski
metric); `minkowski_tensor` the Minkowski metric `diag(1, -1, -1, -1)`;
`pat_matrix(m, dx, dy, dz)` the parallel axis theorem matrix that
translates the inertia matrix of a point mass `m` by `(dx, dy, dz)`; and
`mdft(n)` the discrete Fourier transform matrix (deprecated in SymPy in
favor of `DFT(n)`). After `sympy.physics.matrices`.

## Usage

``` r
msigma(i)

pat_matrix(m, dx, dy, dz)

mgamma(mu, lower = FALSE)

minkowski_tensor

mdft(n)
```

## Arguments

- i:

  The Pauli index.

- m:

  A mass.

- dx, dy, dz:

  The displacement.

- mu:

  The Dirac index.

- lower:

  Whether to lower the index.

- n:

  The size of the transform.

## Value

A `symmatrix`.

## Examples

``` r
msigma(1) %*% msigma(2)
#> Matrix(2, 2)
#> [I, 0 ]
#> [0, -I]
mgamma(5)
#> Matrix(4, 4)
#> [0, 0, 1, 0]
#> [0, 0, 0, 1]
#> [1, 0, 0, 0]
#> [0, 1, 0, 0]
minkowski_tensor
#> Matrix(4, 4)
#> [1, 0 , 0 , 0 ]
#> [0, -1, 0 , 0 ]
#> [0, 0 , -1, 0 ]
#> [0, 0 , 0 , -1]
pat_matrix(2, 1, 0, 0)
#> Matrix(3, 3)
#> [0, 0, 0]
#> [0, 2, 0]
#> [0, 0, 2]
```
