# Solve a system of first order linear ODEs

`linodesolve(A, t, b)` solves `X' = A X + b` for the vector of unknown
functions `X`, as SymPy's `linodesolve()`: constant coefficient systems
by the Jordan form of `A`, and systems with a commutative antiderivative
`B` of `A` by `exp(B)`. The arbitrary constants are Dummy symbols.

## Usage

``` r
linodesolve(A, t, b = NULL, B = NULL, type = "auto", doit = FALSE, tau = NULL)
```

## Arguments

- A:

  The coefficient matrix.

- t:

  The independent variable.

- b:

  The vector of non-homogeneous terms, or `NULL`.

- B:

  A commutative antiderivative of `A`, or `NULL`.

- type:

  One of `"auto"` and `"type1"` to `"type6"`.

- doit:

  Whether to evaluate the integrals of the solution.

- tau:

  The substitution of the independent variable for the types `type5` and
  `type6`.

## Value

A list of the components of the solution.

## Examples

``` r
t <- Symbol("t")
linodesolve(Matrix(list(list(1, 0), list(0, 2))), t)
#> [[1]]
#> _Dummy_1162*exp(t)
#> 
#> [[2]]
#> _Dummy_1163*exp(2*t)
#> 
```
