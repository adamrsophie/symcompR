# Quaternions

`Quaternion(a, b, c, d)` is `a + b*i + c*j + d*k`. `+`, `-`, `*`
(Hamilton product), `/` and `^` with an integer power work; over the
real field a complex scalar `z` acts as the quaternion
`re(z) + im(z)*i`, as in SymPy. `quat_conjugate()`, `quat_norm()`,
`quat_inverse()` and `quat_to_rotation_matrix()` follow SymPy.

## Usage

``` r
Quaternion(a = 0, b = 0, c = 0, d = 0, real_field = TRUE)

quat_conjugate(q)

quat_norm(q)

quat_inverse(q)

quat_to_rotation_matrix(q)
```

## Arguments

- a, b, c, d:

  Components.

- real_field:

  `TRUE` for a quaternion over the reals, `FALSE` for components that
  are complex.

- q:

  A quaternion.

## Value

A `Quaternion`, expression or matrix.

## Examples

``` r
Quaternion(1, 2, 3, 4) * Quaternion(5, 6, 7, 8)
#> (-60) + 12*i + 30*j + 24*k 
Quaternion(1, 2, 3, 4)^2
#> (-28) + 4*i + 6*j + 8*k 
```
