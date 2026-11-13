# Normal ordering of bosonic and fermionic operators

`normal_ordered_form()` moves creation operators to the left of
annihilation operators using the commutation relations, and
`normal_order()` does the same without the commutator terms, as in
`sympy.physics.quantum.operatorordering`.

## Usage

``` r
normal_ordered_form(
  expr,
  independent = FALSE,
  recursive_limit = 10,
  .recursive_depth = 0L
)

normal_order(expr, recursive_limit = 10, .recursive_depth = 0L)
```

## Arguments

- expr:

  An expression of
  [`BosonOp()`](https://adamrsophie.github.io/symcompR/reference/boson.md)
  and
  [`FermionOp()`](https://adamrsophie.github.io/symcompR/reference/fermion.md)
  operators.

- independent:

  Treat operators with different names as independent.

- recursive_limit:

  The maximum number of recursions.

- .recursive_depth:

  The current recursion depth.

## Value

An expression.

## Examples

``` r
a <- BosonOp("a")
normal_ordered_form(a * Dagger(a))
#> 1 + Dagger(a)*a
normal_order(a * Dagger(a))
#> Dagger(a)*a
```
