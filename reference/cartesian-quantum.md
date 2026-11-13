# Cartesian position and momentum

The operators and states of `sympy.physics.quantum.cartesian`: the
position operators `XOp()`, `YOp()`, `ZOp()` and the momentum operator
`PxOp()`, with their eigenstates `XKet()`, `XBra()`, `PxKet()`,
`PxBra()` and the three dimensional `PositionKet3D()` and
`PositionBra3D()`. Inner products give delta functions and plane waves;
[`represent()`](https://adamrsophie.github.io/symcompR/reference/represent.md)
gives the differential operators.

## Usage

``` r
XOp(...)

YOp(...)

ZOp(...)

PxOp(...)

XKet(...)

XBra(...)

PositionState3D(...)

PositionKet3D(...)

PositionBra3D(...)

PxKet(...)

PxBra(...)
```

## Arguments

- ...:

  Labels.

## Value

A quantum operator or state.

## Examples

``` r
x <- Symbol("x"); px <- Symbol("px")
qapply(XOp() * XKet(x))
#> x*|x>
XBra(x) * PxKet(px)
#> <x|px>
```
