# Hilbert spaces

The Hilbert spaces of `sympy.physics.quantum.hilbert`: the abstract
`HilbertSpace()`, the finite or infinite `ComplexSpace(dimension)`, the
square integrable functions `L2(interval)` on an
[`Interval()`](https://adamrsophie.github.io/symcompR/reference/sets.md),
the Fock space `FockSpace()`, and the tensor products, direct sums and
tensor powers built with `*`, `+` and `^` or with the constructors.
`hs_contains(space, obj)` is SymPy's `obj in space`.

## Usage

``` r
HilbertSpace()

ComplexSpace(dimension)

L2(interval)

FockSpace()

TensorProductHilbertSpace(...)

DirectSumHilbertSpace(...)

TensorPowerHilbertSpace(base, exp)

hs_contains(space, obj)
```

## Arguments

- dimension:

  A positive integer, `oo` or a symbol.

- interval:

  An
  [`Interval()`](https://adamrsophie.github.io/symcompR/reference/sets.md).

- ...:

  Hilbert spaces.

- base, exp:

  A Hilbert space and its power.

- space:

  A Hilbert space.

- obj:

  A quantum object.

## Value

A Hilbert space, or `TRUE`/`FALSE` for `hs_contains()`.

## Examples

``` r
ComplexSpace(2) * FockSpace()
#> C(2)*F
ComplexSpace(2)^3
#> C(2)^3
qattr(ComplexSpace(2) + FockSpace(), "dimension")
#> oo
```
