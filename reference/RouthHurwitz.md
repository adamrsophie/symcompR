# Routh-Hurwitz table

`RouthHurwitz(polynomial, var)` is the Routh array of a polynomial: a
matrix with `degree + 1` rows whose first column decides the stability
of the polynomial. A row of zeros is replaced by the derivative of its
auxiliary polynomial and a zero in the first column is handled by the
shifting method. The result is a `symmatrix` with the extra fields
`zero_row_case`, `zero_col_infos` and `auxiliary_polynomials`. After
`sympy.physics.control.routh_table`.

## Usage

``` r
RouthHurwitz(polynomial, var)
```

## Arguments

- polynomial:

  A polynomial expression.

- var:

  Its variable.

## Value

A `RouthHurwitz` matrix.

## Examples

``` r
s <- Symbol("s")
RouthHurwitz(s^4 + s^3 + 3*s^2 + 3*s + 3, s)
#> Matrix(5, 3)
#> [1 , 3, 3]
#> [1 , 3, 0]
#> [-3, 3, 0]
#> [4 , 0, 0]
#> [3 , 0, 0]
```
