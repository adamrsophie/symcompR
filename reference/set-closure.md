# Closure, interior and power set of a set

The set methods `closure`, `interior`, `powerset` and
`symmetric_difference` of SymPy's `Set`: the closure is the union of the
set with its
[`boundary()`](https://adamrsophie.github.io/symcompR/reference/sets.md),
the interior is the set without its boundary, the power set of a finite
set is the finite set of its subsets (any other set gives an unevaluated
`PowerSet`), and the symmetric difference is
[`SymmetricDifference()`](https://adamrsophie.github.io/symcompR/reference/sets.md).

## Usage

``` r
closure(s)

interior(s)

powerset(s)

symmetric_difference(s, t)
```

## Arguments

- s, t:

  Sets.

## Value

A set.

## Examples

``` r
closure(Interval(0, 1, TRUE, TRUE))
#> Interval(0, 1)
interior(Interval(0, 1))
#> Interval.open(0, 1)
powerset(FiniteSet(1, 2))
#> FiniteSet(EmptySet, {1}, {2}, {1, 2})
symmetric_difference(Interval(0, 2), Interval(1, 3))
#> Union(Interval.Ropen(0, 1), Interval.Lopen(2, 3))
```
