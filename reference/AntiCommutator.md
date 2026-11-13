# Anticommutators

`AntiCommutator(A, B)` is the unevaluated `{A, B} = A*B + B*A` of
`sympy.physics.quantum.anticommutator`; it is `2*A*B` when one argument
commutes and its arguments are put in canonical order.
[`doit()`](https://adamrsophie.github.io/symcompR/reference/doit.md)
evaluates it. The commutator `[A, B]` is
[`Commutator()`](https://adamrsophie.github.io/symcompR/reference/diffgeom.md).

## Usage

``` r
AntiCommutator(A, B)
```

## Arguments

- A, B:

  Operators or expressions.

## Value

An expression.

## Examples

``` r
A <- Operator("A"); B <- Operator("B")
AntiCommutator(B, A)
#> {A,B}
doit(AntiCommutator(A, B))
#> A*B + B*A
```
