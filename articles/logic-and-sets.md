# Logic and sets

This vignette covers Boolean algebra (normal forms, minimization, and
satisfiability) and sets: intervals, finite sets, the number sets, image
sets, condition sets, ranges, set operations, and transfinite ordinals.
The predicate based assumptions system
([`ask()`](https://adamrsophie.github.io/symcompR/reference/ask.md),
[`refine()`](https://adamrsophie.github.io/symcompR/reference/refine.md),
and
[`assuming()`](https://adamrsophie.github.io/symcompR/reference/assumptions-context.md))
has its own vignette.

``` r

library(symcompR)
p <- Symbol("p"); q <- Symbol("q"); r <- Symbol("r")
w <- Symbol("w"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
n <- Symbol("n", integer = TRUE)
```

## Boolean expressions

Symbols double as Boolean variables.
[`And()`](https://adamrsophie.github.io/symcompR/reference/And.md),
[`Or()`](https://adamrsophie.github.io/symcompR/reference/And.md),
[`Not()`](https://adamrsophie.github.io/symcompR/reference/And.md),
[`Xor()`](https://adamrsophie.github.io/symcompR/reference/And.md),
[`Implies()`](https://adamrsophie.github.io/symcompR/reference/And.md),
[`Equivalent()`](https://adamrsophie.github.io/symcompR/reference/And.md),
and [`ITE()`](https://adamrsophie.github.io/symcompR/reference/And.md)
(if then else) build Boolean functions, and `&`, `|`, and `!` work as
operators.

``` r

And(p, Or(q, Not(p)))
#> p & (q | ~p)
(p & q) | !r
#> ~r | (p & q)
Implies(p, q)
#> Implies(p, q)
Equivalent(p, q)
#> Equivalent(p, q)
Xor(p, q)
#> p ^ q
And(p, TRUE)
#> p
```

Truth values are substituted with
[`subs()`](https://adamrsophie.github.io/symcompR/reference/subs.md),
and
[`truth_table()`](https://adamrsophie.github.io/symcompR/reference/truth_table.md)
lists all assignments.

``` r

subs(And(p, q), list(list(p, S_true), list(q, S_false)))
#> False
truth_table(Implies(p, q))
#>       p     q value
#> 1 FALSE FALSE  TRUE
#> 2 FALSE  TRUE  TRUE
#> 3  TRUE FALSE FALSE
#> 4  TRUE  TRUE  TRUE
```

### Normal forms and simplification

``` r

to_cnf(Or(And(p, q), r))
#> (p | r) & (q | r)
to_dnf(And(Or(p, q), r))
#> (p & r) | (q & r)
to_nnf(Not(Implies(p, q)))
#> p & ~q
eliminate_implications(Implies(p, q))
#> q | ~p
is_cnf(And(Or(p, q), r))
#> [1] TRUE
to_cnf(Equivalent(p, Or(q, r)), simplify = TRUE)
#> (p | ~q) & (p | ~r) & (q | r | ~p)
```

[`simplify_logic()`](https://adamrsophie.github.io/symcompR/reference/simplify_logic.md)
finds a minimal equivalent expression with the Quine-McCluskey
algorithm.

``` r

simplify_logic(Or(And(p, q), And(p, Not(q))))
#> p
simplify_logic((!p & !q & !r) | (!p & !q & r))
#> ~p & ~q
simplify_logic(Implies(p, q) & Implies(q, p), form = "cnf")
#> (p | ~q) & (q | ~p)
```

[`SOPform()`](https://adamrsophie.github.io/symcompR/reference/SOPform.md)
and
[`POSform()`](https://adamrsophie.github.io/symcompR/reference/SOPform.md)
build a minimal sum of products or product of sums from the minterms of
a truth table (the rows where the function is true), optionally with
“don’t care” rows. Minterms are given as row numbers, so `c(3, 7)` over
`x, y, z` means the rows `011` and `111`.

``` r

SOPform(list(x, y, z), c(3, 7))
#> y & z
SOPform(list(x, y, z), c(1, 2, 5, 6))
#> (y & ~z) | (z & ~y)
POSform(list(w, x, y, z), c(1, 3, 7, 11, 15), dontcares = c(0, 2, 5))
#> z & (y | ~w)
```

### Satisfiability

[`satisfiable()`](https://adamrsophie.github.io/symcompR/reference/satisfiable.md)
runs the DPLL algorithm. It returns a satisfying assignment, or `FALSE`
when there is none.
[`valid()`](https://adamrsophie.github.io/symcompR/reference/satisfiable.md)
tests whether an expression is a tautology.

``` r

satisfiable(And(p, Not(q)))
#> $p
#> [1] TRUE
#> 
#> $q
#> [1] FALSE
satisfiable(And(p, Not(p)))
#> [1] FALSE
length(satisfiable(Or(p, q), all_models = TRUE))
#> [1] 3
valid(Or(p, Not(p)))
#> [1] TRUE
```

The remaining connectives are
[`Nand()`](https://adamrsophie.github.io/symcompR/reference/And.md),
[`Nor()`](https://adamrsophie.github.io/symcompR/reference/And.md), and
the if then else function
[`ITE()`](https://adamrsophie.github.io/symcompR/reference/And.md).
[`is_dnf()`](https://adamrsophie.github.io/symcompR/reference/eliminate_implications.md)
and
[`is_nnf()`](https://adamrsophie.github.io/symcompR/reference/eliminate_implications.md)
test for disjunctive and negation normal form.

``` r

Nand(TRUE, FALSE)
#> True
Nor(TRUE, FALSE)
#> False
Nand(x > 0, y > 0)
#> ~((x > 0) & (y > 0))
Nor(x > 0, y > 0)
#> ~((x > 0) | (y > 0))
ITE(TRUE, p, q)
#> p
c(is_dnf(Or(And(p, q), r)), is_dnf(And(p, Or(q, r))))
#> [1]  TRUE FALSE
c(is_nnf(And(Not(p), q)), is_nnf(Not(And(p, q))))
#> [1]  TRUE FALSE
distribute_or_over_and(And(Or(Not(p), q), r))
#> (q & r) | (r & ~p)
```

[`bool_map()`](https://adamrsophie.github.io/symcompR/reference/bool_map.md)
decides whether two Boolean functions are the same up to a renaming of
variables, and returns the renaming.

``` r

bool_map(Or(And(p, q), Not(p)), Or(And(x, y), Not(x)))
#> $expr
#> q | ~p
#> 
#> $mapping
#> $mapping$p
#> x
#> 
#> $mapping$q
#> y
```

## Assumptions

The assumptions system, with
[`ask()`](https://adamrsophie.github.io/symcompR/reference/ask.md) and
[`refine()`](https://adamrsophie.github.io/symcompR/reference/refine.md),
has its own vignette,
[`vignette("assumptions", package = "symcompR")`](https://adamrsophie.github.io/symcompR/articles/assumptions.md).

## Sets

### Intervals and finite sets

`Interval(a, b)` is closed by default; `left_open` and `right_open`
change that, and
[`Interval_open()`](https://adamrsophie.github.io/symcompR/reference/sets.md),
[`Interval_Lopen()`](https://adamrsophie.github.io/symcompR/reference/sets.md),
and
[`Interval_Ropen()`](https://adamrsophie.github.io/symcompR/reference/sets.md)
are shortcuts.

``` r

Interval(0, 1)
#> Interval(0, 1)
Interval(0, 1, TRUE, FALSE)
#> Interval.Lopen(0, 1)
Interval_open(0, 1)
#> Interval.open(0, 1)
Interval_Lopen(0, 1)
#> Interval.Lopen(0, 1)
Interval_Ropen(0, 1)
#> Interval.Ropen(0, 1)
FiniteSet(1, 2, 3)
#> {1, 2, 3}
```

Set operations simplify their result whenever they can.

``` r

Union(Interval(0, 1), Interval(2, 3), Interval(1, 2))
#> Interval(0, 3)
Union(Interval(0, 1), EmptySet)
#> Interval(0, 1)
Intersection(Interval(0, 1), UniversalSet)
#> Interval(0, 1)
PowerSet(FiniteSet(1, 2))
#> FiniteSet(EmptySet, {1}, {2}, {1, 2})
DisjointUnion(FiniteSet(1, 2), FiniteSet(3))
#> DisjointUnion({1, 2}, {3})
Intersection(Interval(0, 5), Interval(3, 10))
#> Interval(3, 5)
Complement(Interval(0, 5), Interval(1, 2))
#> Union(Interval.Ropen(0, 1), Interval.Lopen(2, 5))
Complement(Reals, FiniteSet(0))
#> Union(Interval.open(-oo, 0), Interval.open(0, oo))
SymmetricDifference(Interval(0, 3), Interval(2, 5))
#> Union(Interval.Ropen(0, 2), Interval.Lopen(3, 5))
Union(FiniteSet(1, 2, 3), FiniteSet(3, 4))
#> {1, 2, 3, 4}
Intersection(FiniteSet(1, 2, 3), Interval(2, 5))
#> {2, 3}
ProductSet(FiniteSet(1, 2), FiniteSet(3, 4))
#> ProductSet({1, 2}, {3, 4})
```

Membership, inclusion, and measure:

``` r

contains(Interval(0, 1), Rational(1, 2))
#> [1] TRUE
contains(Interval(0, 1), x)
#> (0 <= x) & (x <= 1)
is_subset(Interval(0, 1), Interval(-1, 2))
#> [1] TRUE
measure(Union(Interval(0, 1), Interval(2, 4)))
#> 3
boundary(Interval(0, 1, TRUE))
#> {0, 1}
set_inf(Union(Interval(0, 1), Interval(2, 4)))
#> 0
set_sup(Union(Interval(0, 1), Interval(2, 4)))
#> 4
```

[`contains()`](https://adamrsophie.github.io/symcompR/reference/sets.md)
returns a condition when membership depends on a symbol, and
[`as_relational()`](https://adamrsophie.github.io/symcompR/reference/as_relational.md)
turns a set into the equivalent condition.

``` r

as_relational(Union(Interval(0, 1), Interval(2, 3)), x)
#> ((0 <= x) & (x <= 1)) | ((2 <= x) & (x <= 3))
```

### Number sets and ranges

`Naturals`, `Naturals0`, `Integers`, `Rationals`, `Reals`, and
`Complexes` are the standard number sets.
[`Range()`](https://adamrsophie.github.io/symcompR/reference/Range.md)
is an arithmetic progression of integers, possibly infinite.

``` r

is_subset(Naturals, Integers)
#> [1] TRUE
Naturals0
#> Naturals0
Rationals
#> Rationals
Complexes
#> Complexes
contains(Rationals, Rational(1, 2))
#> [1] TRUE
Intersection(Integers, Interval(-2.5, 2))
#> Range(-2, 3, 1)
Range(1, 10, 3)
#> Range(1, 10, 3)
contains(Range(1, 10, 3), 7)
#> [1] TRUE
Intersection(Range(1, 10, 3), Interval(3, 8))
#> Range(4, 10, 3)
Range(0, oo, 2)
#> Range(0, oo, 2)
```

### Image sets and condition sets

`ImageSet(f, S)` is the image \\f(s) : s \in S\\;
[`imageset()`](https://adamrsophie.github.io/symcompR/reference/imageset.md)
simplifies it when possible. `ConditionSet(x, cond, S)` is \\x \in S :
\text{cond}\\. These are the sets that
[`solveset()`](https://adamrsophie.github.io/symcompR/reference/solveset.md)
returns for periodic and unsolved equations.

``` r

squares <- ImageSet(Lambda(x, x^2), Naturals)
squares
#> ImageSet(Lambda(x, x^2), Naturals)
contains(squares, 4)
#> [1] TRUE
imageset(Lambda(x, 2 * x), Integers)
#> ImageSet(Lambda(x, 2*x), Integers)
imageset(Lambda(x, x^2), Interval(-2, 1))
#> Interval(0, 4)
ConditionSet(x, Eq(sin(x), x), Interval(0, 1))
#> ConditionSet(x, Eq(sin(x), x), Interval(0, 1))
```

### Complex regions

[`ComplexRegion()`](https://adamrsophie.github.io/symcompR/reference/ComplexRegion.md)
describes a region of the complex plane as a product of intervals, in
rectangular or polar form.

``` r

R1 <- ComplexRegion(ProductSet(Interval(2, 3), Interval(4, 6)))
contains(R1, Rational(5, 2) + Rational(9, 2) * sym_I)
#> [1] TRUE
Intersection(ComplexRegion(ProductSet(Interval(1, 3), Interval(-1, 1))), Reals)
#> Interval(1, 3)
ComplexRegion(ProductSet(Interval(0, 1), Interval(0, 2 * sym_pi)), polar = TRUE)
#> PolarComplexRegion(ProductSet(Interval(0, 1), Interval.Ropen(0, 2*pi)))
```

[`a_interval()`](https://adamrsophie.github.io/symcompR/reference/ComplexRegion.md),
[`b_interval()`](https://adamrsophie.github.io/symcompR/reference/ComplexRegion.md),
and
[`psets()`](https://adamrsophie.github.io/symcompR/reference/ComplexRegion.md)
return the parts of a region,
[`complex_region_from_real()`](https://adamrsophie.github.io/symcompR/reference/ComplexRegion.md)
embeds a set of reals, and
[`normalize_theta_set()`](https://adamrsophie.github.io/symcompR/reference/ComplexRegion.md)
reduces a set of angles to \[0, 2\pi).

``` r

R2 <- ComplexRegion(ProductSet(Interval(1, 2), Interval(3, 4)))
a_interval(R2)
#> Interval(1, 2)
b_interval(R2)
#> Interval(3, 4)
psets(R2)[[1]]
#> ProductSet(Interval(1, 2), Interval(3, 4))
complex_region_from_real(Interval(1, 2))
#> CartesianComplexRegion(ProductSet(Interval(1, 2), {0}))
normalize_theta_set(Interval(9 * sym_pi / 2, 5 * sym_pi))
#> Interval(pi/2, pi)
normalize_theta_set(Interval(-3 * sym_pi / 2, sym_pi / 2))
#> Interval.Ropen(0, 2*pi)
```

### Set arithmetic

[`SetExpr()`](https://adamrsophie.github.io/symcompR/reference/SetExpr.md)
wraps a set so that arithmetic acts on all of its elements, which is
interval arithmetic for intervals.

``` r

SetExpr(Interval(0, 1)) + 1
#> SetExpr(Interval(1, 2))
SetExpr(Interval(1, 2)) * SetExpr(Interval(10, 20))
#> SetExpr(Interval(10, 40))
```

## Ordinals

symcompR implements ordinal numbers below \varepsilon_0 in Cantor normal
form. `ord_omega` is \omega, and ordinal arithmetic is not commutative:

``` r

w <- ord_omega
w + 1
#> w + 1
1 + w
#> w
1 + w == w
#> [1] TRUE
(w + 3) * (w + 2)
#> w^2 + w*2 + 3
Ordinal(OmegaPower(3, 2), OmegaPower(1, 4), OmegaPower(0, 7))
#> w^3*2 + w*4 + 7
is_successor_ordinal(w^2 + 1)
#> [1] TRUE
is_limit_ordinal(w^2)
#> [1] TRUE
ord0
#> ord0
ord0 + 1
#> 1
ordinal_degree(w^2 + w * 3 + 1)
#> 2
```
