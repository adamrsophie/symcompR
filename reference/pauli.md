# Pauli operators

The Pauli operators `SigmaX()`, `SigmaY()`, `SigmaZ()` and the ladder
operators `SigmaMinus()` and `SigmaPlus()` of
`sympy.physics.quantum.pauli`, labelled by an optional name, and the
states `SigmaZKet(n)` and `SigmaZBra(n)` with `n` 0 or 1.
`qsimplify_pauli()` simplifies products of Pauli operators with the same
name.

## Usage

``` r
SigmaX(...)

SigmaY(...)

SigmaZ(...)

SigmaMinus(...)

SigmaPlus(...)

SigmaZKet(n)

SigmaZBra(n)

qsimplify_pauli(e)
```

## Arguments

- ...:

  The optional name of the operator.

- n:

  0 or 1.

- e:

  An expression.

## Value

A quantum operator, state or expression.

## Examples

``` r
sx <- SigmaX(1); sy <- SigmaY(1)
doit(Commutator(sx, sy))
#> 2*I*SigmaZ()
qsimplify_pauli(sx * sy)
#> I*SigmaZ()
qapply(sx * SigmaZKet(0))
#> |1>
```
