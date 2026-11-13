# Quantum operators

The operators of `sympy.physics.quantum.operator`: a general
`Operator(label)` (`O` by default), `HermitianOperator()`,
`UnitaryOperator()`, the deprecated `IdentityOperator(N)`, the outer
product `OuterProduct(ket, bra)` (also `ket * bra`) and
`DifferentialOperator(expr, f)`, which applies `expr` with `f` replaced
by a
[`Wavefunction()`](https://adamrsophie.github.io/symcompR/reference/Wavefunction.md)
under
[`qapply()`](https://adamrsophie.github.io/symcompR/reference/qapply.md).
Operators do not commute; `qcall(A, "inverse")` is `A^-1`.

## Usage

``` r
Operator(...)

HermitianOperator(...)

UnitaryOperator(...)

IdentityOperator(...)

OuterProduct(...)

DifferentialOperator(...)
```

## Arguments

- ...:

  Labels or arguments.

## Value

An operator.

## Examples

``` r
A <- Operator("A"); B <- Operator("B")
A * B - B * A
#> A*B - B*A
Dagger(UnitaryOperator("U")) * UnitaryOperator("U")
#> 1
OuterProduct(Ket("a"), Bra("b"))
#> |a><b|
```
