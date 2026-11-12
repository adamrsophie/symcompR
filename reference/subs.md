# Substitute values or expressions

`subs()` replaces subexpressions and re-evaluates the result, like
SymPy's `expr.subs()`. Substitutions can be given as
`subs(e, old, new)`, as a named list `subs(e, list(x = 1, y = 2))`
(names refer to free symbols), as a list of pairs
`list(list(old, new), ...)`, or as equations. Pairs are applied in
sequence unless `simultaneous = TRUE`. `xreplace()` performs exact
structural replacement without algebraic matching.

## Usage

``` r
subs(e, ..., simultaneous = FALSE)

xreplace(e, rule)
```

## Arguments

- e:

  A `sym` expression.

- ...:

  Substitutions.

- simultaneous:

  Apply all substitutions at once.

- rule:

  A named list or list of pairs mapping old to new expressions.

## Value

A `sym` expression.

## Examples

``` r
x <- Symbol("x"); y <- Symbol("y")
subs(x^2 + y, x, 2)
#> y + 4
subs(x^2 + y, list(x = 1, y = 3))
#> 4
subs(x * y + 1, x * y, Symbol("z"))
#> z + 1
```
