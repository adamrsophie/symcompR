# Linearizer of equations of motion

`Linearizer()` holds the general form of a system's equations of motion
and linearizes them about an operating point, as SymPy's `Linearizer`.
It is usually made by the `to_linearizer()` method of
[`KanesMethod()`](https://adamrsophie.github.io/symcompR/reference/KanesMethod.md)
or
[`LagrangesMethod()`](https://adamrsophie.github.io/symcompR/reference/LagrangesMethod.md);
`linearize(op_point, A_and_B, simplify)` (reached with `$`) returns `M`,
`A` and `B`, or `A` and `B`.

## Usage

``` r
Linearizer(
  f_0,
  f_1,
  f_2,
  f_3,
  f_4,
  f_c,
  f_v,
  f_a,
  q,
  u,
  q_i = NULL,
  q_d = NULL,
  u_i = NULL,
  u_d = NULL,
  r = NULL,
  lams = NULL,
  linear_solver = "LU"
)
```

## Arguments

- f_0, f_1, f_2, f_3, f_4, f_c, f_v, f_a:

  The matrices of the general form.

- q, u:

  Generalized coordinates and speeds.

- q_i, q_d, u_i, u_d:

  Independent and dependent coordinates and speeds.

- r:

  Input variables.

- lams:

  Lagrange multipliers.

- linear_solver:

  `"LU"` or a function of `A` and `b`.

## Value

A `Linearizer` object.

## Examples

``` r
q <- dynamicsymbols("q"); u <- dynamicsymbols("u")
qd <- dynamicsymbols("q", 1); ud <- dynamicsymbols("u", 1)
k <- Symbol("k")
lin <- Linearizer(list(qd), list(-u), list(ud), list(k * q), list(0),
                  list(), list(), list(), list(q), list(u),
                  q_i = list(q), u_i = list(u))
lin$linearize(A_and_B = TRUE)
#> (Matrix([[0, 1], [-k, 0]]), Matrix(0, 0, [])) 
```
