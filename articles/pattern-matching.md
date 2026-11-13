# Pattern matching, unification, and strategies

``` r

library(symcompR)
x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
a <- Symbol("a"); b <- Symbol("b"); p <- Symbol("p"); q <- Symbol("q")
f <- Function("f"); g <- Function("g")
```

This vignette covers pattern matching with wildcards, `sympy.unify`,
`sympy.strategies`, and the expression paths of
`sympy.simplify.epathtools`.

## Pattern matching

[`Wild()`](https://adamrsophie.github.io/symcompR/reference/Dummy.md)
creates a pattern variable.
[`sym_match()`](https://adamrsophie.github.io/symcompR/reference/sym_match.md)
returns the values of the wildcards for which the pattern equals the
expression, or `NULL` when there are none. The `exclude` argument keeps
a wildcard from matching expressions containing given symbols.

``` r

wa <- Wild("a", exclude = list(x)); wb <- Wild("b", exclude = list(x))
sym_match(3 * x + 2, wa * x + wb)
#> $a
#> 3
#> 
#> $b
#> 2
sym_match(x^2 + 5, wa * x + wb)
#> NULL
```

[`sym_replace()`](https://adamrsophie.github.io/symcompR/reference/sym_match.md)
replaces every match of a pattern:

``` r

u <- Wild("u")
sym_replace(sin(x)^2 + 1, sin(u), cos(u))
#> cos(x)^2 + 1
sym_replace(f(x) + f(y), f(u), g(u))
#> g(x) + g(y)
```

[`WildFunction()`](https://adamrsophie.github.io/symcompR/reference/WildFunction.md)
matches any function application:

``` r

F <- WildFunction("F")
sym_match(f(x), F)
#> $F
#> f(x)
sym_match(x + cos(x), x + F)
#> $F
#> cos(x)
```

## Unification

[`unify()`](https://adamrsophie.github.io/symcompR/reference/unify.md)
finds all ways to make two expressions equal by assigning the given
variables, taking commutativity into account.

``` r

unify(a + b, x + y, variables = list(x, y))
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
```

[`rewriterule()`](https://adamrsophie.github.io/symcompR/reference/unify.md)
turns a pattern and a replacement into a rule, which returns all the
rewritten forms. A condition can restrict the matches.

``` r

rewriterule(p^2, p^3, variables = list(p))(x^2)
#> [[1]]
#> x^3
rewriterule(p^2 + q^3, (p * q)^4, list(p, q))(x^2 + y^3)
#> [[1]]
#> x^4*y^4
rewriterule(x, x + 1, list(x), function(x) isTRUE(x < 10))(Integer(5))
#> [[1]]
#> 1 + 5
rewriterule(x, x + 1, list(x), function(x) isTRUE(x < 10))(Integer(15))
#> list()
```

The unification engine works on generic trees.
[`unify_compound()`](https://adamrsophie.github.io/symcompR/reference/unify.md)
builds a node from an operator and arguments,
[`unify_variable()`](https://adamrsophie.github.io/symcompR/reference/unify.md)
a variable, and
[`unify_condvariable()`](https://adamrsophie.github.io/symcompR/reference/unify.md)
a variable with a condition.
[`unify_core()`](https://adamrsophie.github.io/symcompR/reference/unify.md)
returns all the substitutions that unify two trees.

``` r

v <- unify_variable("v")
t1 <- unify_compound("Add", list(1, 2))
t2 <- unify_compound("Add", list(v, 2))
unify_core(t1, t2)
#> [[1]]
#> [[1]]$`V(Lv)`
#> [[1]]$`V(Lv)`$k
#> $arg
#> [1] "v"
#> 
#> attr(,"class")
#> [1] "unify_variable"
#> 
#> [[1]]$`V(Lv)`$v
#> [1] 1
cv <- unify_condvariable("w", function(e) e > 1)
unify_core(t1, unify_compound("Add", list(1, cv)))
#> [[1]]
#> [[1]]$`C(Lw)`
#> [[1]]$`C(Lw)`$k
#> $arg
#> [1] "w"
#> 
#> $valid
#> function (e) 
#> e > 1
#> 
#> attr(,"class")
#> [1] "unify_condvariable"
#> 
#> [[1]]$`C(Lw)`$v
#> [1] 2
unify_core(t1, unify_compound("Add", list(cv, 2)))
#> list()
```

## Rewriting with strategies

`sympy.strategies` provides combinators for building rewrite systems out
of simple rules. A rule is a function from an expression to an
expression;
[`exhaust()`](https://adamrsophie.github.io/symcompR/reference/strategies.md)
applies a rule until nothing changes,
[`chain()`](https://adamrsophie.github.io/symcompR/reference/strategies.md)
applies rules in sequence, and
[`greedy()`](https://adamrsophie.github.io/symcompR/reference/strategies.md)
explores a tree of alternatives and keeps the best result by an
objective function.

``` r

inc <- function(e) if (isTRUE(is_integer(e)) && as.numeric(e) < 5) e + 1 else e
exhaust(inc)(Integer(1))
#> 5
chain(inc, inc)(Integer(0))
#> 2
best <- greedy(list(function(e) expand(e), function(e) sym_factor(e)),
               objective = count_ops)
best(x^2 + 2 * x + 1)
#> (x + 1)^2
```

### Rules

A rule takes an expression and returns a new one. The rules of
`sympy.strategies.rl` work on any node;
[`strat_new()`](https://adamrsophie.github.io/symcompR/reference/strategies.md)
builds one with a given operator, here a generic `Basic` node.

``` r

B <- function(...) strat_new(list(type = "Basic"), ...)
e <- B(S(1), S(0), S(2))
e
#> Basic(1, 0, 2)
rm_id(function(a) a == 0)(e)
#> Basic(1, 2)
unpack(B(S(2)))
#> 2
rl$flatten(B(S(1), B(S(2), S(3))))
#> Basic(1, 2, 3)
rl$sort(function(a) as.numeric(a))(B(S(3), S(1), S(2)))
#> Basic(1, 2, 3)
glom(function(a) a, function(a) 1, function(cnt, a) cnt * a)(B(S(2), S(2), S(3)))
#> Basic(4, 3)
distribute("Mul", "Add")(strat_new(list(type = "Mul"), S(2), strat_new(list(type = "Add"), x, y)))
#> 2*x + 2*y
rebuild(strat_new(list(type = "Add"), x, x))
#> 2*x
rl$subs(S(1), S(5))(S(1))
#> 5
```

`basic_fns` and `expr_fns` tell the traversals how to take a node apart
and put it back together; `expr_fns` evaluates the rebuilt expression.

``` r

names(basic_fns)
#> [1] "op"       "new"      "leaf"     "children"
traverse$top_down(function(e) if (isTRUE(e == 1)) S(9) else e)(B(S(1), B(S(1), S(2))))
#> Basic(9, Basic(9, 2))
traverse$bottom_up(function(e) if (isTRUE(is_integer(e))) e + 1 else e)(B(S(1), S(2)))
#> Basic(2, 3)
traverse$top_down(function(e) if (isTRUE(e == 1)) S(9) else e, fns = expr_fns)(x + 1)
#> x + 9
```

### Combinators

[`condition()`](https://adamrsophie.github.io/symcompR/reference/strategies.md)
applies a rule only when a predicate holds,
[`do_one()`](https://adamrsophie.github.io/symcompR/reference/strategies.md)
applies the first rule that changes the expression,
[`minimize()`](https://adamrsophie.github.io/symcompR/reference/strategies.md)
keeps the best of several results,
[`strat_switch()`](https://adamrsophie.github.io/symcompR/reference/strategies.md)
dispatches on a key,
[`typed()`](https://adamrsophie.github.io/symcompR/reference/strategies.md)
dispatches on the type of the expression,
[`memoize()`](https://adamrsophie.github.io/symcompR/reference/strategies.md)
caches results,
[`null_safe()`](https://adamrsophie.github.io/symcompR/reference/strategies.md)
and
[`tryit()`](https://adamrsophie.github.io/symcompR/reference/strategies.md)
return the input when a rule fails, and
[`canon()`](https://adamrsophie.github.io/symcompR/reference/strategies.md)
applies rules everywhere until nothing changes.
[`strat_identity()`](https://adamrsophie.github.io/symcompR/reference/strategies.md)
is the rule that does nothing, and
[`strat_debug()`](https://adamrsophie.github.io/symcompR/reference/strategies.md)
prints each step.

``` r

inc <- function(e) e + 1; dec <- function(e) e - 1; dbl <- function(e) 2 * e
condition(function(e) e > 2, inc)(S(1))
#> 1
do_one(strat_identity, inc)(S(1))
#> 2
minimize(inc, dec, objective = as.numeric)(S(3))
#> 2
strat_switch(function(e) e$type, list(Integer = inc))(S(1))
#> 2
by_type <- typed(list(Integer = inc, Symbol = dbl))
vapply(lapply(list(S(3), x, x + 1), by_type), sstr, "")
#> [1] "4"     "2*x"   "x + 1"
memoize(inc)(S(1))
#> 2
null_safe(function(e) NULL)(S(1))
#> 1
tryit(function(e) stop("no"))(S(1))
#> 1
canon(function(e) if (identical(e$type, "Basic") && length(e$args) == 1) e$args[[1]] else e)(B(B(S(1))))
#> 1
strat_identity(x)
#> x
strat_debug(inc)(S(1))
#> Rule: inc
#> In:   1
#> Out:  2
#> 2
```

### Trees of strategies

A list of rules is a choice and
[`strat_seq()`](https://adamrsophie.github.io/symcompR/reference/strategies.md)
a sequence.
[`greedy()`](https://adamrsophie.github.io/symcompR/reference/strategies.md)
takes the best choice at each step,
[`brute()`](https://adamrsophie.github.io/symcompR/reference/strategies.md)
tries all of them, and
[`allresults()`](https://adamrsophie.github.io/symcompR/reference/strategies.md)
returns every result.
[`treeapply()`](https://adamrsophie.github.io/symcompR/reference/strategies.md)
evaluates a tree with given functions for choices and sequences, and
`branch` holds the branching versions of the combinators, which return
lists of results.

``` r

tree <- list(inc, strat_seq(dec, dbl))
greedy(tree, objective = as.numeric)(S(4))
#> 5
greedy(tree, objective = as.numeric)(S(1))
#> 0
brute(tree, objective = as.numeric)(S(4))
#> 5
allresults(tree)(S(4))
#> [[1]]
#> 5
#> 
#> [[2]]
#> 6
treeapply(list(strat_seq(3, 2), strat_seq(4, 1)), list(list = max, seq = min))
#> [1] 2
branch$multiplex(function(e) list(e + 1), function(e) list(e - 1))(S(3))
#> [[1]]
#> 4
#> 
#> [[2]]
#> 2
```

## Expression paths

[`epath()`](https://adamrsophie.github.io/symcompR/reference/EPath.md)
selects or modifies parts of an expression with a path language similar
to XPath: `/*` descends into all arguments, `/Name` filters by class,
and `/[i]` selects by position.
[`EPath()`](https://adamrsophie.github.io/symcompR/reference/EPath.md)
stores a compiled path.

``` r

EPath("/*/[0]")
#> EPath('/*/[0]')
e <- Symbol("t") + sin(x + 1) + cos(x + y + sym_E)
epath("/*/*/Symbol", e)
#> [[1]]
#> x
#> 
#> [[2]]
#> x
#> 
#> [[3]]
#> y
epath("/*/*/Symbol", e, function(s) 2 * s)
#> t + sin(2*x + 1) + cos(2*x + 2*y + E)
nested <- list(list(list(x, 1), 2), list(list(3, y), z))
epath("/*/[0]/Symbol", nested)
#> [[1]]
#> x
#> 
#> [[2]]
#> y
epath("/*/[-1]", nested)
#> [[1]]
#> [1] 2
#> 
#> [[2]]
#> z
```
