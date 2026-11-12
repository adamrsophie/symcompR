# Ordinals

`Ordinal(...)` is an ordinal number in Cantor normal form, the sum of
terms `OmegaPower(exp, mult)` (that is `w^exp * mult`) with decreasing
exponents, as SymPy's `Ordinal`. `ord0` is the ordinal zero and
`ord_omega` is omega (SymPy's `omega`). `+`, `*` and comparisons follow
ordinal arithmetic, so `1 + ord_omega == ord_omega`; `ord_omega^a`
raises omega to an ordinal or integer. `is_successor_ordinal()`,
`is_limit_ordinal()` and `ordinal_degree()` are SymPy's properties.

## Usage

``` r
OmegaPower(exp, mult)

Ordinal(...)

ord0

ord_omega

is_successor_ordinal(o)

is_limit_ordinal(o)

ordinal_degree(o)

# S3 method for class 'sym_ordinal'
Ops(e1, e2)

# S3 method for class 'sym_ordinal'
chooseOpsMethod(x, y, mx, my, cl, reverse)

# S3 method for class 'sym_omega_power'
Ops(e1, e2)

# S3 method for class 'sym_ordinal'
format(x, ...)

# S3 method for class 'sym_ordinal'
print(x, ...)

# S3 method for class 'sym_omega_power'
format(x, ...)

# S3 method for class 'sym_omega_power'
print(x, ...)
```

## Arguments

- exp:

  An ordinal or a nonnegative integer.

- mult:

  A positive integer.

- ...:

  `OmegaPower` terms, with decreasing exponents.

- o:

  An ordinal.

- e1, e2:

  Ordinals or integers.

- x:

  An ordinal or omega power.

- y, mx, my, cl, reverse:

  Arguments of the
  [`chooseOpsMethod()`](https://rdrr.io/r/base/chooseOpsMethod.html)
  generic.

## Value

An ordinal, an omega power, or a logical.

## Examples

``` r
w <- ord_omega
(w + 3) * (w + 2)
#> w^2 + w*2 + 3
1 + w == w
#> [1] TRUE
Ordinal(OmegaPower(3, 2), OmegaPower(1, 4), OmegaPower(0, 7))
#> w^3*2 + w*4 + 7
```
