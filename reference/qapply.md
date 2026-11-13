# Apply operators to states

`qapply(e)` applies the operators in `e` to the states they act on, as
SymPy's `sympy.physics.quantum.qapply`: products are evaluated from the
right, commutators and tensor products are expanded, and inner products
are evaluated unless `ip_doit = FALSE`. With `dagger = TRUE` operators
are also applied to bras on their left; `sum_doit = TRUE` evaluates
sums.

## Usage

``` r
qapply(e, dagger = FALSE, sum_doit = FALSE, ip_doit = TRUE, ...)
```

## Arguments

- e:

  An expression of operators and states.

- dagger:

  Apply operators to the left too.

- sum_doit:

  Evaluate `Sum` objects.

- ip_doit:

  Evaluate inner products.

- ...:

  Further options for the operators.

## Value

An expression.

## Examples

``` r
qapply(OuterProduct(Ket("a"), Bra("b")) * Ket("c"))
#> <b|c>*|a>
```
