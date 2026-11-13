# Wigner small d and D matrices

`wigner_d_small(J, beta)` is the small Wigner d matrix of a rotation by
`beta` about the y axis for angular momentum `J`;
`wigner_d(J, alpha, beta, gamma)` the Wigner D matrix of the rotation
with Euler angles `alpha`, `beta`, `gamma`. After
`sympy.physics.wigner`.

## Usage

``` r
wigner_d_small(J, beta)

wigner_d(J, alpha, beta, gamma)
```

## Arguments

- J:

  An integer or half integer.

- alpha, beta, gamma:

  Euler angles.

## Value

A `symmatrix`.

## Examples

``` r
wigner_d_small(S_Half, S_Pi / 2)
#> Matrix(2, 2)
#> [sqrt(2)/2 , sqrt(2)/2]
#> [-sqrt(2)/2, sqrt(2)/2]
wigner_d(S_Half, Symbol("a"), Symbol("b"), Symbol("c"))
#> Matrix(2, 2)
#> [exp(I*a/2)*exp(I*c/2)*cos(b/2)  , exp(I*a/2)*exp(-I*c/2)*sin(b/2) ]
#> [-exp(-I*a/2)*exp(I*c/2)*sin(b/2), exp(-I*a/2)*exp(-I*c/2)*cos(b/2)]
```
