# One-dimensional harmonic oscillator

The ladder operators `RaisingOp()` and `LoweringOp()`, the number
operator `NumberOp()`, the Hamiltonian `Hamiltonian()` and the states
`SHOKet(n)` and `SHOBra(n)` of `sympy.physics.quantum.sho1d`. Each
operator takes one label.
[`represent()`](https://adamrsophie.github.io/symcompR/reference/represent.md)
gives matrices of size `ndim` (4 by default) in the number basis, and
[`rewrite()`](https://adamrsophie.github.io/symcompR/reference/rewrite.md)
writes the operators in terms of `"a"`, `"xp"`, `"N"` or `"H"`.

## Usage

``` r
RaisingOp(...)

LoweringOp(...)

NumberOp(...)

Hamiltonian(...)

SHOKet(n)

SHOBra(n)
```

## Arguments

- ...:

  The label.

- n:

  The number of quanta.

## Value

A quantum operator or state.

## Examples

``` r
ad <- RaisingOp("a"); a <- LoweringOp("a")
qapply(ad * SHOKet(1))
#> sqrt(2)*|2>
represent(NumberOp("N"), ndim = 3)
#> Matrix(3, 3)
#> [0, 0, 0]
#> [0, 1, 0]
#> [0, 0, 2]
doit(Commutator(a, ad))
#> 1
```
