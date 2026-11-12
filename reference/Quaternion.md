# Quaternions

`Quaternion(a, b, c, d)` is `a + b*i + c*j + d*k`. `+`, `-`, `*`
(Hamilton product) and `/` by a scalar work; `quat_conjugate()`,
`quat_norm()`, `quat_inverse()` and `quat_to_rotation_matrix()` follow
SymPy.

## Usage

``` r
Quaternion(a = 0, b = 0, c = 0, d = 0)

quat_conjugate(q)

quat_norm(q)

quat_inverse(q)

quat_to_rotation_matrix(q)
```

## Arguments

- a, b, c, d:

  Components.

- q:

  A quaternion.

## Value

A `Quaternion`, expression or matrix.

## Examples

``` r
Quaternion(1, 2, 3, 4) * Quaternion(5, 6, 7, 8)
#> -60 + 12*i + 30*j + 24*k 
```
