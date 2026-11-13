# Quantum Fourier transform

`QFT(start, finish)` and `IQFT(start, finish)` are the quantum Fourier
transform and its inverse on the qubits `start` to `finish - 1`, and
`RkGate(target, k)` (also `Rk()`) the phase gate they decompose into, as
in `sympy.physics.quantum.qft`.

## Usage

``` r
QFT(...)

IQFT(...)

RkGate(...)

Rk(...)
```

## Arguments

- ...:

  The qubit range, or the target and `k`.

## Value

A gate.

## Examples

``` r
qcall(QFT(0, 3), "decompose")
#> SWAP(0,2)*H(0)*C((0),S(1))*H(1)*C((0),T(2))*C((1),S(2))*H(2)
represent(QFT(0, 2), nqubits = 2)
#> Matrix(4, 4)
#> [1/2, 1/2 , 1/2 , 1/2 ]
#> [1/2, I/2 , -1/2, -I/2]
#> [1/2, -1/2, 1/2 , -1/2]
#> [1/2, -I/2, -1/2, I/2 ]
RkGate(1, 2)
#> S(1)
```
