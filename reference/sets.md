# Sets

Intervals, finite sets, the standard number sets and set operations,
like SymPy's sets module. `Interval(a, b)` is closed by default; use
`left_open`/`right_open` or `Interval_open()`, `Interval_Lopen()` and
`Interval_Ropen()`. Number sets are `Reals`, `Integers`, `Naturals`
(starting at 1), `Naturals0`, `Rationals`, `Complexes` and
`UniversalSet`. `Union()`, `Intersection()`, `Complement()`,
`SymmetricDifference()` and `ProductSet()` combine sets and simplify
when possible. `contains()` tests membership; `measure()`, `boundary()`,
`set_inf()`, `set_sup()` and `is_subset()` query sets.

## Usage

``` r
Interval(a, b, left_open = FALSE, right_open = FALSE)

contains(s, x)

is_subset(s, t)

Union(...)

Intersection(...)

Complement(s, t)

SymmetricDifference(s, t)

measure(s)

boundary(s)

set_inf(s)

set_sup(s)

ImageSet(lam, base, ...)

Reals

Integers

Naturals

Naturals0

Rationals

Complexes

UniversalSet

Interval_open(a, b)

Interval_Lopen(a, b)

Interval_Ropen(a, b)

FiniteSet(...)

EmptySet

ProductSet(...)

ConditionSet(sym, condition, base = UniversalSet)
```

## Arguments

- a, b:

  Endpoints.

- left_open, right_open:

  Openness of the endpoints.

- s, t:

  Sets.

- x:

  An element.

- ...:

  Elements or sets; for `ImageSet()`, further sets, one for each
  argument of a `Lambda` with several arguments.

- lam:

  A
  [`Lambda()`](https://adamrsophie.github.io/symcompR/reference/Lambda.md).

- base:

  A set.

- sym:

  A symbol.

- condition:

  A boolean condition.

## Value

A set expression, or a logical/boolean value.

## Examples

``` r
Interval(0, 1)
#> Interval(0, 1)
Union(Interval(0, 1), Interval(2, 3), Interval(1, 2))
#> Interval(0, 3)
Intersection(Interval(0, 5), Interval(3, 10))
#> Interval(3, 5)
contains(Interval(0, 1), Rational(1, 2))
#> [1] TRUE
Complement(Reals, FiniteSet(0))
#> Union(Interval.open(-oo, 0), Interval.open(0, oo))
```
