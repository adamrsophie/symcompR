# Branching strategies

Combinators for branching rules after `sympy.strategies.branch`. A
branching rule maps an expression to a list of results (SymPy's
generators). `branch_identity()` returns its input, `branch_exhaust()`
applies a rule until it has no effect, `onaction()` calls a function on
each result that differs from the input and `branch_debug()` prints the
rule, input and output of each, `multiplex()` merges the distinct
results of several rules, `branch_condition()` applies a rule when a
predicate holds, `sfilter()` keeps the results that satisfy a predicate,
`notempty()` returns the input when a rule gives nothing,
`branch_do_one()` uses the first rule that gives results,
`branch_chain()` composes rules and `yieldify()` turns a rule into a
branching rule. `branch_top_down()` applies a rule down an expression
tree, `branch_sall()` applies it to the arguments, and `branch_canon()`
applies rules top down until nothing changes.

## Usage

``` r
branch_identity(x)

branch_exhaust(brule)

onaction(brule, fn)

branch_debug(brule, file = NULL)

multiplex(...)

branch_condition(cond, brule)

sfilter(pred, brule)

notempty(brule)

branch_do_one(...)

branch_chain(...)

yieldify(rl)

branch_canon(...)

branch_top_down(brule, fns = basic_fns)

branch_sall(brule, fns = basic_fns)
```

## Arguments

- x:

  An expression.

- brule, ...:

  Branching rules.

- fn:

  A function `fn(brule, expr, result)`.

- file:

  `NULL` to print to the console, or a connection to write to.

- cond, pred:

  A predicate.

- rl:

  A rule.

- fns:

  Traversal functions (`op`, `new`, `children`, `leaf`), as
  [basic_fns](https://adamrsophie.github.io/symcompR/reference/strategies.md).

## Value

A branching rule; `branch_identity()` returns a list.

## Examples

``` r
posdec <- function(x) if (x > 0) list(x - 1L) else list(x)
branch5 <- function(x) {
  if (0 < x && x < 5) list(x - 1L) else if (5 < x && x < 10) list(x + 1L)
  else if (x == 5) list(x + 1L, x - 1L) else list(x)
}
branch_exhaust(branch5)(5L)
#> [[1]]
#> [1] 10
#> 
#> [[2]]
#> [1] 0
#> 
multiplex(posdec, branch5)(5L)
#> [[1]]
#> [1] 4
#> 
#> [[2]]
#> [1] 6
#> 
branch_chain(branch5, yieldify(function(x) x * 2L))(5L)
#> [[1]]
#> [1] 12
#> 
#> [[2]]
#> [1] 8
#> 
x <- Symbol("x")
inc <- function(e) if (is_sym(e) && e$type == "Integer") list(e + 1L) else list()
branch_top_down(inc)(x + 1)
#> [[1]]
#> x + 2
#> 
```
