# Core constructors and helpers

`Add()`, `Mul()` and `Pow()` build canonical sums, products and powers.
`Equality()`, `Unequality()`, `StrictLessThan()`, `LessThan()`,
`StrictGreaterThan()` and `GreaterThan()` are the long names of
[`Eq()`](https://adamrsophie.github.io/symcompR/reference/Eq.md) and
friends, and `Rel(a, b, op)` builds a relation from its operator string.
`Subs()` is an unevaluated substitution (evaluate with
[`doit()`](https://adamrsophie.github.io/symcompR/reference/doit.md)).
`nfloat()` evaluates the numbers in an expression and keeps symbols.
`bottom_up()` applies a function to every node, leaves first.
`default_sort_key()` returns the canonical ordering key and
`sym_ordered()` sorts expressions by it. `comp()` compares two numbers
within a tolerance. `check_assumptions()` and `failing_assumptions()`
test facts. `sstrrepr()` is
[`sstr()`](https://adamrsophie.github.io/symcompR/reference/sstr.md),
`true` and `false` are the boolean atoms.

## Usage

``` r
Add(...)

Mul(...)

Pow(a, b)

Equality(a, b)

Unequality(a, b)

StrictLessThan(a, b)

LessThan(a, b)

StrictGreaterThan(a, b)

GreaterThan(a, b)

Rel(a, b, op = "==")

Subs(e, old, new)

nfloat(e, n = 15, exponent = FALSE, dkeys = FALSE)

bottom_up(e, f)

default_sort_key(e)

comp(a, b, tol = 1e-12)

check_assumptions(e, ...)

failing_assumptions(e, ...)

sstrrepr(e)

true

false

sym_ordered(x, ...)
```

## Arguments

- ...:

  Expressions, or assumption facts as `name = value`.

- a, b, e:

  Expressions.

- op:

  One of `"=="`, `"!="`, `"<"`, `"<="`, `">"`, `">="`.

- old, new:

  Old and new values for `Subs()`.

- n:

  Significant digits for `nfloat()`.

- exponent:

  If `TRUE`, `nfloat()` also converts exponents.

- dkeys:

  If `TRUE`, `nfloat()` also converts the names of a list.

- f:

  A function.

- tol:

  Relative tolerance.

- x:

  A list of expressions.

## Value

An expression, list, logical or string.

## Examples

``` r
x <- Symbol("x")
Add(x, 1, x)
#> 2*x + 1
doit(Subs(x^2, x, 3))
#> 9
check_assumptions(x^2 + 1, positive = TRUE)
#> [1] NA
```
