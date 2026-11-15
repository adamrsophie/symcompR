# Symbolic equations and inequalities

`Eq(a, b)` is the equation `a = b`, `Ne()` is `a != b`, and `Lt()`,
`Le()`, `Gt()`, `Ge()` are the strict and non strict inequalities. When
the truth value is decidable they return `S_true` or `S_false`. Use
`lhs()` and `rhs()` to take them apart.

## Usage

``` r
Eq(a, b = 0, evaluate = .ev_default())

Ne(a, b, evaluate = .ev_default())

Lt(a, b, evaluate = .ev_default())

Le(a, b, evaluate = .ev_default())

Gt(a, b, evaluate = .ev_default())

Ge(a, b, evaluate = .ev_default())

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
