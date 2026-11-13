# Shor's algorithm

`shor(N)` factors `N` with the period finding circuit of
`sympy.physics.quantum.shor`; `period_find()` runs that circuit,
`getr()` and `ratioize()` turn a measurement into a period through its
continued fraction. The modular exponentiation gate `CMod` is not
completed in SymPy either, so `shor()` and `period_find()` raise a
`NotImplementedError`.

## Usage

``` r
CMod(...)

shor(N)

getr(x, y, N)

ratioize(list, N)

period_find(a, N)
```

## Arguments

- ...:

  The arguments of the `CMod` gate.

- N:

  The number to factor, or the bound on the period.

- x, y:

  The numerator and denominator of the measured fraction.

- list:

  The continued fraction terms.

- a:

  The base.

## Value

A factor, a period or an expression.

## Examples

``` r
getr(513, 1024, 10)
#> 2
```
