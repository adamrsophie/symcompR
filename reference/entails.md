# Entailment and knowledge bases

`entails(expr, formula_set)` tells whether the formulas imply `expr`.
`PropKB()` is a propositional knowledge base: `kb$tell(sentence)` adds
the clauses of a sentence, `kb$ask(query)` asks whether they entail the
query, `kb$retract(sentence)` removes clauses and `kb$clauses()` lists
them in SymPy's order. `KB()` is the abstract base, whose methods signal
`NotImplementedError`. These follow `sympy.logic.inference`.

## Usage

``` r
entails(expr, formula_set = NULL)

KB(sentence = NULL)

PropKB(sentence = NULL)
```

## Arguments

- expr:

  A boolean expression.

- formula_set:

  A list of boolean expressions.

- sentence:

  A sentence to tell first.

## Value

`TRUE` or `FALSE`; a knowledge base.

## Examples

``` r
A <- Symbol("A"); B <- Symbol("B")
entails(B, list(Implies(A, B), A))
#> [1] TRUE
kb <- PropKB()
kb$tell(Implies(A, B)); kb$tell(A)
kb$ask(B)
#> [1] TRUE
```
