# Normal ordering, contractions and Wick's theorem

`NO()` is the normal ordering of a product of fermion operators.
`contraction()` is the contraction of two fermion operators, `wicks()`
expands a product of fermion operators by Wick's theorem,
`evaluate_deltas()` removes the Kronecker deltas of a product by
substitution, `substitute_dummies()` renames the dummy indices of an
expression into a canonical form, `PermutationOperator()` is the
operator exchanging two indices and `simplify_index_permutations()`
gathers the terms related by such permutations, as
`sympy.physics.secondquant`.

## Usage

``` r
NO(arg)

contraction(a, b)

evaluate_deltas(e)

substitute_dummies(expr, new_indices = FALSE, pretty_indices = list())

wicks(
  e,
  simplify_kronecker_deltas = FALSE,
  expand = TRUE,
  simplify_dummies = FALSE,
  keep_only_fully_contracted = FALSE
)

PermutationOperator(i, j)

simplify_index_permutations(expr, permutation_operators)
```

## Arguments

- arg, expr, e, a, b:

  Expressions or operators.

- new_indices:

  Whether `substitute_dummies()` makes new dummies.

- pretty_indices:

  A list with the letters `above`, `below` and `general` for the new
  dummies.

- simplify_kronecker_deltas, expand, simplify_dummies,
  keep_only_fully_contracted:

  Options of `wicks()`.

- i, j:

  Indices.

- permutation_operators:

  A list of `PermutationOperator()`s.

## Value

An expression.

## Examples

``` r
i <- Symbol("i", below_fermi = TRUE); a <- Symbol("a", above_fermi = TRUE)
p <- Symbol("p"); q <- Symbol("q")
NO(sym_F(i) * Fd(a))
#> -NO(CreateFermion(a)*AnnihilateFermion(i))
contraction(Fd(i), sym_F(i))
#> 1
wicks(sym_F(p) * Fd(q))
#> KroneckerDelta(_a, q)*KroneckerDelta(p, q) - NO(CreateFermion(q)*AnnihilateFermion(p))
```
