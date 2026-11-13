# Strategies

Combinators for rewriting rules after `sympy.strategies`. A rule maps an
expression to an expression. `chain()` applies rules in order,
`do_one()` applies the first rule that changes the expression,
`exhaust()` repeats a rule until nothing changes, `condition()` applies
a rule when a predicate holds, `null_safe()` keeps the input when a rule
returns `NULL`, `tryit()` keeps it when a rule fails, `minimize()` keeps
the result with the smallest objective, `memoize()` caches results,
`strat_switch()` picks a rule by a key, and `canon()` applies rules top
down to a fixed point. The lists `rl` (expression rules `rm_id`, `glom`,
`sort`, `distribute`, `subs`, `unpack`, `flatten`, `rebuild`),
`traverse` (`top_down`, `bottom_up`, `top_down_once`, `bottom_up_once`,
`sall`) and `branch` (the same combinators for branching rules) mirror
SymPy's modules. `greedy(tree, objective)` evaluates a tree of rules
(lists are alternatives, `strat_seq()` are sequences) keeping the best
alternative, `allresults()` returns every result and `brute()` the best
of all.

## Usage

``` r
strat_identity(x)

strat_seq(...)

exhaust(rule)

memoize(rule)

condition(cond, rule)

chain(...)

strat_debug(rule, file = NULL)

null_safe(rule)

tryit(rule, exception = NULL)

do_one(...)

strat_switch(key, ruledict)

minimize(..., objective = strat_identity)

basic_fns

expr_fns

rl

rm_id(isid, new = function(op, ...) .strat_new(op, list(...)))

unpack(expr)

glom(key, count, combine)

distribute(A, B)

rebuild(expr)

strat_new(op, ...)

traverse

canon(..., fns = basic_fns)

typed(ruletypes)

branch

treeapply(tree, join, leaf = strat_identity)

greedy(tree, objective = strat_identity, leaf = strat_identity)

allresults(tree, leaf = branch$yieldify)

brute(tree, objective = strat_identity, leaf = branch$yieldify)
```

## Arguments

- x:

  An expression (`strat_identity()` returns it unchanged).

- ...:

  Rules.

- rule:

  A rule.

- cond:

  A predicate.

- file:

  `NULL` to print to the console, or a file name or connection to append
  the log to.

- exception:

  Ignored (all errors are caught).

- key:

  A function giving a name.

- ruledict:

  A named list of rules.

- objective:

  A function giving a number or a vector compared lexicographically.

- isid:

  For `rm_id()`: a predicate for identity elements.

- new:

  A constructor `new(op, ...)` for rebuilt expressions.

- expr:

  An expression.

- count, combine:

  For `glom()`: a function giving a count and one combining the counts.

- A, B:

  For `distribute()`: the outer and inner operations.

- op:

  An operation name.

- fns:

  Traversal functions (`op`, `new`, `children`, `leaf`).

- ruletypes:

  A named list of rules by expression type.

- tree:

  A tree of rules.

- join:

  A named list with functions for `list` and `seq` nodes.

- leaf:

  A function applied to the rules at the leaves.

## Value

A rule, or for `allresults()` a branching rule.

## Examples

``` r
x <- Symbol("x")
inc <- function(e) if (isTRUE(is_integer(e)) && as.numeric(e) < 5) e + 1 else e
exhaust(inc)(Integer(1))
#> 5
chain(inc, inc)(Integer(0))
#> 2
greedy(list(function(e) expand(e), function(e) sym_factor(e)), objective = count_ops)((x + 1)^2)
#> (x + 1)^2
```
