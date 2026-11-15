# Binary relation predicates

The classes of SymPy's `Q.eq`, `Q.ne`, `Q.gt`, `Q.ge`, `Q.lt` and
`Q.le`: each constructor returns that predicate, as calling the class
does in SymPy. `AppliedBinaryRelation(p, lhs, rhs)` applies one,
`binrel_properties(p)` gives SymPy's `is_reflexive`, `is_symmetric`,
`reversed` and `negated`, and `BinaryRelation(name)` returns the
predicate of that name.

## Usage

``` r
BinaryRelation(name)

EqualityPredicate()

UnequalityPredicate()

StrictGreaterThanPredicate()

GreaterThanPredicate()

StrictLessThanPredicate()

LessThanPredicate()

AppliedBinaryRelation(p, lhs, rhs)

binrel_properties(p)
```

## Arguments

- name:

  One of `"eq"`, `"ne"`, `"gt"`, `"ge"`, `"lt"`, `"le"`.

- p:

  A binary relation predicate.

- lhs, rhs:

  Expressions.

## Value

A predicate, an applied predicate, or a list.

## Examples

``` r
x <- Symbol("x")
AppliedBinaryRelation(StrictLessThanPredicate(), x, 2)
#> Q.lt(x, 2)
binrel_properties(Q$gt)
#> $is_reflexive
#> [1] FALSE
#> 
#> $is_symmetric
#> [1] FALSE
#> 
#> $reversed
#> Q.lt
#> 
#> $negated
#> Q.le
#> 
```
