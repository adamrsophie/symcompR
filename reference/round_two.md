# Integral bases, prime decomposition and valuations in number fields

`round_two(T)` computes the maximal order of the number field defined by
the irreducible polynomial `T` with the Round Two algorithm, returning
the order as a submodule of the power basis and the field discriminant.
`prime_decomp(p, T)` factors the rational prime `p` into prime ideals
and `prime_valuation(I, P)` is the valuation of an ideal at a prime
ideal, as SymPy's functions of the same names.

## Usage

``` r
round_two(T)

prime_decomp(p, T)

prime_valuation(I, P)
```

## Arguments

- T:

  An irreducible polynomial over ZZ or QQ.

- p:

  A rational prime.

- I:

  A submodule (ideal) of the maximal order.

- P:

  A prime ideal from `prime_decomp()`.

## Value

For `round_two()` a list `(ZK, dK)`; for `prime_decomp()` a list of
prime ideals; for `prime_valuation()` an integer.

## Examples

``` r
x <- Symbol("x")
round_two(x^3 + x^2 - 2*x + 8)
#> (Submodule[[2, 0, 0], [0, 2, 0], [0, 1, 1]]/2, -503) 
prime_decomp(2, x^3 + x^2 - 2*x + 8)
#> [[ (2, x^2 + 1) e=1, f=1 ], [ (2, (x^2 + 3*x + 2)/2) e=1, f=1 ], [ (2, (3*x^2 + 3*x)/2) e=1, f=1 ]] 
```
