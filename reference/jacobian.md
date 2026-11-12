# Jacobian, Hessian, gradient and Wronskian

Jacobian, Hessian, gradient and Wronskian

## Usage

``` r
jacobian(f, vars)

hessian(f, vars)

gradient(f, vars, doit = TRUE)

wronskian(f, x)
```

## Arguments

- f:

  An expression or list of expressions.

- vars:

  A list of variables. For `gradient()` without `vars`, the gradient of
  a
  [`CoordSys3D()`](https://adamrsophie.github.io/symcompR/reference/CoordSys3D.md)
  scalar or vector field (`sympy.vector`).

- doit:

  Evaluate the derivatives (coordinate system fields only).

- x:

  The variable for the Wronskian.

## Value

A symbolic matrix (see
[`Matrix()`](https://adamrsophie.github.io/symcompR/reference/Matrix.md)).

## Examples

``` r
x <- Symbol("x"); y <- Symbol("y")
jacobian(list(x^2 * y, 5 * x + sin(y)), list(x, y))
#> Matrix(2, 2)
#> [2*x*y, x^2   ]
#> [5    , cos(y)]
hessian(x^3 + x * y^2, list(x, y))
#> Matrix(2, 2)
#> [6*x, 2*y]
#> [2*y, 2*x]
```
