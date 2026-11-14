# Indexed complex roots of polynomials

`CRootOf(f, x, k)` is the `k`-th root (0 based) of the polynomial `f` in
`x`, ordered with real roots first in increasing order and then complex
roots by real and imaginary part, as in SymPy. It evaluates numerically
with [`N()`](https://adamrsophie.github.io/symcompR/reference/N.md).

## Usage

``` r
CRootOf(f, x, k, radicals = FALSE)
```

## Arguments

- f:

  A polynomial.

- x:

  The variable; when `k` is missing, `x` is the index and the variable
  is the single free symbol of `f`.

- k:

  Index of the root, from 0.

- radicals:

  As in SymPy: when `TRUE`, roots of quadratic and binomial factors are
  given by radicals.

## Value

A `sym` CRootOf object, or the root itself when it is rational (or a
radical, with `radicals = TRUE`).
