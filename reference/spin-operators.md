# Spin operators, rotations and Wigner D functions

The angular momentum operators of `sympy.physics.quantum.spin`:
`JxOp()`, `JyOp()`, `JzOp()`, `J2Op()`, the ladder operators `JplusOp()`
and `JminusOp()` (each takes a name; the objects `Jx`, `Jy`, `Jz`, `J2`,
`Jplus` and `Jminus` are these operators named `J`, as in SymPy), the
rotation `Rotation(alpha, beta, gamma)` and the Wigner D function
`WignerD(j, m, mp, alpha, beta, gamma)` (evaluated by
[`doit()`](https://adamrsophie.github.io/symcompR/reference/doit.md)).
`Rotation_D()` and `Rotation_d()` are SymPy's `Rotation.D` and
`Rotation.d`.

## Usage

``` r
JxOp(...)

JyOp(...)

JzOp(...)

J2Op(...)

JplusOp(...)

JminusOp(...)

Jx

Jy

Jz

J2

Jplus

Jminus

Rotation(...)

Rotation_D(j, m, mp, alpha, beta, gamma)

Rotation_d(j, m, mp, beta)

WignerD(j, m, mp, alpha, beta, gamma, evaluate = FALSE)
```

## Arguments

- ...:

  The name of the operator, or the Euler angles.

- j, m, mp:

  Angular momentum quantum numbers.

- alpha, beta, gamma:

  Euler angles.

- evaluate:

  Evaluate the Wigner D function.

## Value

A quantum operator or an expression.

## Examples

``` r
qapply(Jz * JzKet(1, 1))
#> hbar*|1,1>
represent(JzOp("J"), j = 1)
#> Matrix(3, 3)
#> [hbar, 0, 0    ]
#> [0   , 0, 0    ]
#> [0   , 0, -hbar]
doit(WignerD(1, 1, 0, 0, sym_pi / 2, 0))
#> -sqrt(2)/2
```
