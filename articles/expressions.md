# Working with expression trees

Every symcompR expression is a tree: an operation (its type) applied to
arguments, down to atoms such as symbols and numbers. This vignette
shows how to inspect that tree, substitute into it, match patterns
against it, rewrite it with rules, and keep parts of it unevaluated.

``` r

library(symcompR)
x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
a <- Symbol("a"); b <- Symbol("b")
f <- Function("f"); g <- Function("g")
```

Many functions here return lists of expressions. A small helper prints
them on one line.

``` r

show <- function(l) vapply(l, sstr, "")
```

## Inspecting structure

``` r

e <- x^2 * sin(y) + 3
sym_type(e)
#> [1] "Add"
show(sym_args(e))
#> [1] "3"          "x^2*sin(y)"
has(e, y)
#> [1] TRUE
has(e, z)
#> [1] FALSE
show(free_symbols(e))
#> [1] "x" "y"
show(atoms(e))
#> [1] "3" "y" "x" "2"
show(atoms(e, "Symbol"))
#> [1] "y" "x"
show(atoms(e, "Integer"))
#> [1] "3" "2"
show(atoms(e, "Function"))
#> [1] "sin(y)"
```

Bound variables, such as an integration variable, are not free:

``` r

show(free_symbols(Integral(x * y, list(x, 0, 1))))
#> [1] "y"
```

[`count_ops()`](https://adamrsophie.github.io/symcompR/reference/count_ops.md)
counts the operations in an expression; with `visual = TRUE` it shows
which ones. The simplification functions use it as a measure of
complexity.

``` r

count_ops(x^2 + sin(x) / y)
#> [1] 4
count_ops(x^2 + sin(x) / y, visual = TRUE)
#> ADD + DIV + POW + SIN
```

[`postorder_traversal()`](https://adamrsophie.github.io/symcompR/reference/preorder_traversal.md)
and
[`preorder_traversal()`](https://adamrsophie.github.io/symcompR/reference/preorder_traversal.md)
visit every subexpression.
[`bottom_up()`](https://adamrsophie.github.io/symcompR/reference/core-helpers.md)
rebuilds an expression by applying a function to every subexpression,
starting from the leaves.

``` r

length(postorder_traversal(e))
#> [1] 8
vapply(preorder_traversal(sin(cos(x))), sstr, "")
#> [1] "sin(cos(x))" "cos(x)"      "x"
bottom_up(sin(2 * x) + cos(2 * x)^2, expand_trig)
#> (2*cos(x)^2 - 1)^2 + 2*sin(x)*cos(x)
```

## Substitution

[`subs()`](https://adamrsophie.github.io/symcompR/reference/subs.md)
replaces subexpressions; replacements given in a list are applied one
after the other, unless `simultaneous = TRUE`.
[`xreplace()`](https://adamrsophie.github.io/symcompR/reference/subs.md)
replaces exact subtrees only, without any mathematical matching.

``` r

subs(x^2 + y, x, 2)
#> y + 4
subs(x + y * z, y * z, a)
#> a + x
subs(x + y, list(list(x, y), list(y, x)))
#> 2*x
subs(x + y, list(list(x, y), list(y, x)), simultaneous = TRUE)
#> x + y
xreplace(x * y + 1, list(list(x * y, z)))
#> z + 1
```

[`Subs()`](https://adamrsophie.github.io/symcompR/reference/core-helpers.md)
is an unevaluated substitution, used for example for derivatives
evaluated at a point.

``` r

Subs(diff(f(x), x), x, 0)
#> Subs(Derivative(f(x), x), x, 0)
doit(Subs(sin(x) * x, x, 0))
#> 0
```

## Pattern matching and rewriting

Wildcards, unification, rule based rewriting, and expression paths are
covered in
[`vignette("pattern-matching", package = "symcompR")`](https://adamrsophie.github.io/symcompR/articles/pattern-matching.md).

## Rewriting in terms of other functions

``` r

rewrite(exp(x), "cos")
#> sinh(x) + cosh(x)
```

## Extracting parts

``` r

extract_multiplicatively(2 * x, 2)
#> x
extract_multiplicatively((x * y)^3, x^2 * y)
#> x*y^2
extract_additively(x + 3, 2)
#> x + 1
as_coefficient(2 * sym_I * x, sym_I)
#> 2*x
```

## Relations

[`Eq()`](https://adamrsophie.github.io/symcompR/reference/Eq.md),
[`Ne()`](https://adamrsophie.github.io/symcompR/reference/Eq.md),
[`Lt()`](https://adamrsophie.github.io/symcompR/reference/Eq.md),
[`Le()`](https://adamrsophie.github.io/symcompR/reference/Eq.md),
[`Gt()`](https://adamrsophie.github.io/symcompR/reference/Eq.md), and
[`Ge()`](https://adamrsophie.github.io/symcompR/reference/Eq.md) create
relations. They evaluate to `TRUE` or `FALSE` when the answer is
determined, and stay symbolic otherwise.

``` r

r <- Eq(x^2, 4)
lhs(r)
#> x^2
rhs(r)
#> 4
Ne(x, 1)
#> Ne(x, 1)
Gt(x, 1) & Lt(x, 3)
#> (x > 1) & (x < 3)
subs(Eq(x + 1, 3), x, 2)
#> True
Lt(2, 3)
#> True
```

Note that `==` compares structure, not mathematical equality:

``` r

expand((x + 1)^2) == x^2 + 2 * x + 1
#> [1] TRUE
(x + 1)^2 == x^2 + 2 * x + 1
#> [1] FALSE
```

## Unevaluated expressions

[`UnevaluatedExpr()`](https://adamrsophie.github.io/symcompR/reference/UnevaluatedExpr.md)
protects a subexpression from automatic simplification until
[`doit()`](https://adamrsophie.github.io/symcompR/reference/doit.md) is
called.

``` r

x * UnevaluatedExpr(1 / x)
#> x*1/x
doit(x * UnevaluatedExpr(1 / x))
#> 1
```

## Anonymous functions

[`Lambda()`](https://adamrsophie.github.io/symcompR/reference/Lambda.md)
creates a symbolic function;
[`rcall()`](https://adamrsophie.github.io/symcompR/reference/diffgeom.md)
applies it.

``` r

sq <- Lambda(x, x^2)
sq
#> Lambda(x, x^2)
rcall(sq, 3)
#> 9
apply_lambda(sq, 3)
#> 9
```

## Singleton constants

The `S_*` objects are the singletons of SymPy’s `S` registry. They are
ordinary symcompR expressions, and `TribonacciConstant` is one more
named constant.

``` r

show(list(S_Zero, S_One, S_NegativeOne, S_Half, S_Pi, S_E, S_I,
          S_Infinity, S_NegativeInfinity, S_ComplexInfinity, S_NaN,
          S_EmptySet))
#>  [1] "0"        "1"        "-1"       "1/2"      "pi"       "E"       
#>  [7] "I"        "oo"       "-oo"      "zoo"      "nan"      "EmptySet"
N(TribonacciConstant)
#> 1.83928675521416
nan
#> nan
zoo
#> zoo
```

## Small constructors and helpers

[`Add()`](https://adamrsophie.github.io/symcompR/reference/core-helpers.md),
[`Mul()`](https://adamrsophie.github.io/symcompR/reference/core-helpers.md),
and
[`Pow()`](https://adamrsophie.github.io/symcompR/reference/core-helpers.md)
build a sum, a product, and a power,
[`ln()`](https://adamrsophie.github.io/symcompR/reference/function-aliases.md)
is the natural logarithm,
[`Id()`](https://adamrsophie.github.io/symcompR/reference/function-aliases.md)
is the identity, and
[`Tuple()`](https://adamrsophie.github.io/symcompR/reference/Tuple.md)
holds a fixed sequence of expressions.
[`sym_name()`](https://adamrsophie.github.io/symcompR/reference/sym_name.md)
returns the name of a symbol,
[`is_sym()`](https://adamrsophie.github.io/symcompR/reference/is_sym.md)
tests for a symcompR object,
[`as_sym()`](https://adamrsophie.github.io/symcompR/reference/as_sym.md)
converts R values and strings,
[`numbered_symbols()`](https://adamrsophie.github.io/symcompR/reference/iterables.md)
makes a run of fresh symbols, and
[`sym_prod()`](https://adamrsophie.github.io/symcompR/reference/misc-helpers.md)
multiplies a list.
[`sym_poly()`](https://adamrsophie.github.io/symcompR/reference/poly-aliases.md)
is the [`poly()`](https://rdrr.io/r/stats/poly.html) function of SymPy.

``` r

Add(x, y, 2)
#> x + y + 2
Mul(2, x, y)
#> 2*x*y
Pow(x, 2)
#> x^2
ln(x)
#> log(x)
Id(x)
#> x
Tuple(1, x)
#> (1, x)
sym_name(x)
#> [1] "x"
c(is_sym(x), is_sym(2))
#> [1]  TRUE FALSE
as_sym("x + 1")
#> x + 1
show(numbered_symbols("y", 3))
#> [1] "y0" "y1" "y2"
sym_prod(list(1, 2, 3, x))
#> 6*x
sym_poly(x^2 + 2 * x + 1)
#> Poly(x^2 + 2*x + 1, x, domain='ZZ')
```

[`sym_ordered()`](https://adamrsophie.github.io/symcompR/reference/core-helpers.md)
sorts expressions into SymPy’s canonical order,
[`sym_reshape()`](https://adamrsophie.github.io/symcompR/reference/misc-helpers.md)
regroups a flat list, and
[`nfloat()`](https://adamrsophie.github.io/symcompR/reference/core-helpers.md)
turns the rational coefficients of an expression into floats while
leaving its powers exact.
[`comp()`](https://adamrsophie.github.io/symcompR/reference/core-helpers.md)
compares two numbers to a tolerance.

``` r

show(sym_ordered(list(x^2, x, S(1), y)))
#> [1] "1"   "x"   "y"   "x^2"
vapply(sym_reshape(as.list(0:5), 2), function(r) paste(unlist(r), collapse = " "), "")
#> [1] "0 1" "2 3" "4 5"
nfloat(x^2 + 2 * x / 3)
#> x^2 + 0.666666666666667*x
c(comp(1, 1.001, 0.01), comp(1, 1.1, 0.01))
#> [1]  TRUE FALSE
```

[`sym_rf()`](https://adamrsophie.github.io/symcompR/reference/function-aliases.md)
and
[`ff()`](https://adamrsophie.github.io/symcompR/reference/function-aliases.md)
are the rising and falling factorials; with an integer second argument
they expand into products.

``` r

sym_rf(x, 3)
#> x*(x + 1)*(x + 2)
ff(x, 3)
#> x*(x - 2)*(x - 1)
expand(expand_func(ff(x, 3)))
#> x^3 - 3*x^2 + 2*x
```

## Integer helpers

`integer_log(y, x)` returns the largest `e` with `x^e <= y` and whether
the power is exact,
[`num_digits()`](https://adamrsophie.github.io/symcompR/reference/integer-helpers.md)
counts digits in a base,
[`trailing()`](https://adamrsophie.github.io/symcompR/reference/integer-helpers.md)
counts the trailing zero bits, and
[`arity()`](https://adamrsophie.github.io/symcompR/reference/integer-helpers.md)
counts the arguments of an R function.

``` r

str(integer_log(125, 5))
#> List of 2
#>  $ e    : num 3
#>  $ exact: logi TRUE
str(integer_log(126, 5))
#> List of 2
#>  $ e    : num 3
#>  $ exact: logi FALSE
c(num_digits(123456), num_digits(255, 16), trailing(48))
#> [1] 6 2 4
arity(function(a, b) a)
#> [1] 2
```

## Content, primitive parts, and signs

[`add_primitive()`](https://adamrsophie.github.io/symcompR/reference/extract-methods.md)
and
[`as_content_primitive()`](https://adamrsophie.github.io/symcompR/reference/as_content_primitive.md)
pull the rational content out of a sum,
[`gcd_terms()`](https://adamrsophie.github.io/symcompR/reference/factor_terms.md)
factors out common terms without full factorization, and
[`extract_branch_factor()`](https://adamrsophie.github.io/symcompR/reference/extract-methods.md)
splits a polar exponential into a principal part and a number of turns.
[`could_extract_minus_sign()`](https://adamrsophie.github.io/symcompR/reference/could_extract_minus_sign.md)
tells whether an expression reads more naturally with its sign flipped.

``` r

show(add_primitive(2 * x + 4))
#> [1] "2"     "x + 2"
show(as_content_primitive(2 * x + 6 * y))
#> [1] "2"       "x + 3*y"
gcd_terms((x + 1)^2 * y + (x + 1) * y^2)
#> y*(x + 1)*(x + y + 1)
show(extract_branch_factor(exp_polar(5 * sym_I * sym_pi)))
#> [1] "exp_polar(I*pi)" "2"
c(could_extract_minus_sign(-x - y), could_extract_minus_sign(x - y))
#> [1]  TRUE FALSE
```

[`factor_nc()`](https://adamrsophie.github.io/symcompR/reference/factor_nc.md)
factors expressions with noncommutative symbols, keeping the order of
the factors.

``` r

P <- Symbol("P", commutative = FALSE); Qn <- Symbol("Q", commutative = FALSE)
factor_nc(P * x + P * y)
#> (x + y)*P
factor_nc(x * P * Qn + x * Qn * P)
#> x*(P*Q + Q*P)
```

## Relation constructors

The relational classes have constructors under their SymPy class names,
and
[`Rel()`](https://adamrsophie.github.io/symcompR/reference/core-helpers.md)
picks the class from an operator string.
[`Contains()`](https://adamrsophie.github.io/symcompR/reference/imageset.md)
states that an element lies in a set; it evaluates when the answer is
known.

``` r

Equality(x, 1)
#> Eq(x, 1)
GreaterThan(x, 1)
#> x >= 1
LessThan(x, 1)
#> x <= 1
StrictGreaterThan(x, 1)
#> x > 1
StrictLessThan(x, 1)
#> x < 1
Unequality(x, 1)
#> Ne(x, 1)
Le(x, 1)
#> x <= 1
Rel(x, 1, "<")
#> x < 1
Contains(x, Interval(0, 1))
#> Contains(x, Interval(0, 1))
Contains(Rational(1, 2), Interval(0, 1))
#> True
```

## Dictionaries

[`Dict()`](https://adamrsophie.github.io/symcompR/reference/Dict.md) is
SymPy’s immutable dictionary of expressions, with
[`dict_get()`](https://adamrsophie.github.io/symcompR/reference/Dict.md),
[`dict_keys()`](https://adamrsophie.github.io/symcompR/reference/Dict.md),
[`dict_values()`](https://adamrsophie.github.io/symcompR/reference/Dict.md),
and
[`dict_items()`](https://adamrsophie.github.io/symcompR/reference/Dict.md).

``` r

d <- Dict(list(x, y), list(1, 2))
d
#> {x: 1, y: 2}
dict_get(d, x)
#> 1
show(dict_keys(d))
#> [1] "x" "y"
show(dict_values(d))
#> [1] "1" "2"
show(dict_items(d))
#> [1] "(x, 1)" "(y, 2)"
```

## Exact numbers as gmp objects

[`as_bigq()`](https://adamrsophie.github.io/symcompR/reference/as_bigq.md)
and
[`as_bigz()`](https://adamrsophie.github.io/symcompR/reference/as_bigq.md)
return the exact value of a rational or integer as a `gmp` number.

``` r

as_bigq(Rational(1, 3))
#> Big Rational ('bigq') :
#> [1] 1/3
as_bigz(Integer(10)^20)
#> Big Integer ('bigz') :
#> [1] 100000000000000000000
```

## Dummy symbols

[`Dummy()`](https://adamrsophie.github.io/symcompR/reference/Dummy.md)
creates a symbol that is never equal to any other symbol, even one with
the same name. It is used internally for temporary variables.

``` r

Symbol("d") == Symbol("d")
#> [1] TRUE
Dummy("d") == Dummy("d")
#> [1] FALSE
```

## Indexed objects

[`IndexedBase()`](https://adamrsophie.github.io/symcompR/reference/Indexed.md)
and
[`Idx()`](https://adamrsophie.github.io/symcompR/reference/Indexed.md)
represent subscripted quantities such as matrix entries, and
[`get_indices()`](https://adamrsophie.github.io/symcompR/reference/Indexed.md)
finds the free indices of a product, following the summation convention
for repeated indices.

``` r

i <- Idx("i"); j <- Idx("j")
A <- IndexedBase("A"); v <- IndexedBase("v")
A(i, j) * v(j)
#> A[i, j]*v[j]
show(get_indices(A(i, j) * v(j))$indices)
#> [1] "i"
```
