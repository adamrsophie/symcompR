# Wavefunctions

`Wavefunction(expr, ...)` is a function of the coordinates given after
the expression, each a symbol or `list(x, lower, upper)`, as SymPy's
`sympy.physics.quantum.state.Wavefunction`. `qcall(f, "__call__", ...)`
evaluates it (zero outside the limits), `qattr(f, "norm")`,
`qcall(f, "normalize")` and `qcall(f, "prob")` normalize it.

## Usage

``` r
Wavefunction(...)
```

## Arguments

- ...:

  The expression, then the coordinates.

## Value

A `Wavefunction`.

## Examples

``` r
x <- Symbol("x"); L <- Symbol("L", positive = TRUE)
f <- Wavefunction(sqrt(2 / L) * sin(sym_pi * x / L), list(x, 0, L))
qattr(f, "norm")
#> 1
```
