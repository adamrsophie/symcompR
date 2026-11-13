# Quantum states

Kets and bras of `sympy.physics.quantum.state`. `Ket(...)` and
`Bra(...)` take labels (strings become symbols); with no label the label
is `psi`. `TimeDepKet()` and `TimeDepBra()` take a time as the last
label; `OrthogonalKet()` and `OrthogonalBra()` are orthogonal for
different labels. `QExpr()`, `StateBase()`, `KetBase()`, `BraBase()`,
`State()` and `TimeDepState()` are the base classes. `Bra * Ket` is an
[`InnerProduct()`](https://adamrsophie.github.io/symcompR/reference/InnerProduct.md),
`Ket * Bra` an
[`OuterProduct()`](https://adamrsophie.github.io/symcompR/reference/quantum-operators.md),
and
[`Dagger()`](https://adamrsophie.github.io/symcompR/reference/Dagger.md)
turns one into the other.

## Usage

``` r
QExpr(...)

StateBase(...)

KetBase(...)

BraBase(...)

State(...)

Ket(...)

Bra(...)

TimeDepState(...)

TimeDepKet(...)

TimeDepBra(...)

OrthogonalState(...)

OrthogonalKet(...)

OrthogonalBra(...)
```

## Arguments

- ...:

  Labels.

## Value

A quantum state.

## Examples

``` r
k <- Ket("psi"); b <- Bra("phi")
b * k
#> <phi|psi>
k * b
#> |psi><phi|
Dagger(k)
#> <psi|
TimeDepKet("psi", "t")
#> |psi;t>
```
