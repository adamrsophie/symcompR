# Hermitian conjugate

`Dagger(arg)` takes the Hermitian conjugate of a quantum object, an
operator expression, a matrix or a scalar, as SymPy's
`sympy.physics.quantum.Dagger`: kets become bras, products reverse their
order and numbers are conjugated. With `evaluate = FALSE` the conjugate
stays unevaluated.

## Usage

``` r
Dagger(arg, evaluate = TRUE)
```

## Arguments

- arg:

  An expression or a matrix.

- evaluate:

  Evaluate the conjugate.

## Value

An expression or a matrix.

## Examples

``` r
A <- Operator("A"); B <- Operator("B")
Dagger(A * B)
#> Dagger(B)*Dagger(A)
Dagger(Ket("psi"))
#> <psi|
Dagger(sym_I * A)
#> -I*Dagger(A)
```
