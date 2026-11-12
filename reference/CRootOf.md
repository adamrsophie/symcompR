# Indexed complex roots of polynomials

`CRootOf(f, x, k)` is the `k`-th root (0 based) of the polynomial `f` in
`x`, ordered with real roots first in increasing order and then complex
roots by real and imaginary part, as in SymPy. It evaluates numerically
with [`N()`](https://adamrsophie.github.io/symcompR/reference/N.md).

## Usage

``` r
CRootOf(f, x, k)
```

## Arguments

- f:

  A polynomial.

- x:

  The variable.

- k:

  Index of the root, from 0.

## Value

A `sym` CRootOf object.
