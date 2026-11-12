# Symbolic equations and inequalities

`Eq(a, b)` is the equation `a = b`, `Ne()` is `a != b`, and `Lt()`,
`Le()`, `Gt()`, `Ge()` are the strict and non strict inequalities. When
the truth value is decidable they return `S_true` or `S_false`. Use
`lhs()` and `rhs()` to take them apart.

## Usage

``` r
Eq(a, b = 0, evaluate = TRUE)

Ne(a, b, evaluate = TRUE)

Lt(a, b, evaluate = TRUE)

Le(a, b, evaluate = TRUE)

Gt(a, b, evaluate = TRUE)

Ge(a, b, evaluate = TRUE)

lhs(e)

rhs(e)
```

## Arguments

- a, b:

  Expressions.

- evaluate:

  If `FALSE`, keep the relation unevaluated.

- e:

  A relational expression.

## Value

A relational `sym`, or a boolean atom.

## Examples

``` r
x <- Symbol("x")
Eq(x^2, 4)
#> Eq(x^2, 4)
Lt(1, 2)
#> True
```
