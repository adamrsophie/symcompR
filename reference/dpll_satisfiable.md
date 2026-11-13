# The DPLL algorithm

The building blocks of the Davis Putnam Logemann Loveland algorithm in
`sympy.logic.algorithms.dpll`: `dpll_satisfiable(expr)` returns a model
(a named list) or `FALSE`; `dpll(clauses, symbols, model)` searches from
a partial model; `find_pure_symbol()`, `find_unit_clause()` and
`unit_propagate()` are its steps. The `_int_repr` variants work on
clauses given as integer vectors (positive and negative literals), and
`to_int_repr(clauses, symbols)` encodes clauses so; `pl_true_int_repr()`
evaluates one such clause.

## Usage

``` r
dpll_satisfiable(expr)

dpll(clauses, symbols, model)

find_pure_symbol(symbols, unknown_clauses)

find_unit_clause(clauses, model)

unit_propagate(clauses, symbol)

to_int_repr(clauses, symbols)

pl_true_int_repr(clause, model = list())

unit_propagate_int_repr(clauses, s)

find_pure_symbol_int_repr(symbols, unknown_clauses)

find_unit_clause_int_repr(clauses, model)

dpll_int_repr(clauses, symbols, model)
```

## Arguments

- expr:

  A boolean expression.

- clauses:

  A list of clauses.

- symbols:

  A list of symbols, or an integer vector.

- model:

  A named list of logical values (keyed by symbol names, or by integers
  for the `_int_repr` functions).

- unknown_clauses:

  A list of clauses.

- symbol:

  A literal.

- clause:

  An integer clause.

- s:

  An integer literal.

## Value

A model, `FALSE`, a pair `list(symbol, value)` or a list of clauses.

## Examples

``` r
A <- Symbol("A"); B <- Symbol("B")
dpll_satisfiable(And(A, Not(B)))
#> $B
#> [1] FALSE
#> 
#> $A
#> [1] TRUE
#> 
find_pure_symbol(list(A, B), list(Or(A, Not(B)), Or(Not(B), A)))
#> [[1]]
#> A
#> 
#> [[2]]
#> [1] TRUE
#> 
unit_propagate(list(Or(A, B), Or(Not(A), B)), A)
#> [[1]]
#> B
#> 
```
