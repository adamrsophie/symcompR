# Unification

`unify(x, y, s, variables)` lists the substitutions (named lists from
the variables to expressions) that make `x` and `y` equal, treating the
symbols in `variables` as wildcards; sums, products, unions,
intersections and finite sets match associatively, and commutatively
where SymPy's do.
`rewriterule(source, target, variables, condition, assume)` returns a
function giving the rewrites of an expression. `unify_core()` is the
generic algorithm on `unify_compound()`, `unify_variable()` and
`unify_condvariable()` trees.

## Usage

``` r
unify_compound(op, args)

unify_variable(arg)

unify_condvariable(arg, valid)

unify_core(x, y, s = list(), fns = list())

unify(x, y, s = list(), variables = list(), ...)

rewriterule(
  source,
  target,
  variables = list(),
  condition = NULL,
  assume = NULL
)
```

## Arguments

- op, args, arg, valid:

  Parts of a tree.

- x, y:

  Expressions or trees.

- s:

  A starting substitution.

- fns:

  Named list with `is_associative` and `is_commutative`.

- variables:

  Symbols that act as wildcards.

- ...:

  Unused.

- source, target:

  Patterns.

- condition:

  A predicate on the matched values.

- assume:

  A boolean assumption the match must satisfy.

## Value

A list of substitutions, or a rewrite function.

## Examples

``` r
x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
unify(Symbol("a") + Symbol("b"), x + y, variables = list(x, y))
#> [[1]]
#> [[1]]$x
#> a
#> 
#> [[1]]$y
#> b
#> 
#> 
#> [[2]]
#> [[2]]$y
#> a
#> 
#> [[2]]$x
#> b
#> 
#> 
```
