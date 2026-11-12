# Where an intersection with a finite set is not empty

`not_empty_in(s, x)` takes the intersection of a finite set of
expressions in `x` with an interval or a union of intervals, and returns
the set of real `x` for which the intersection is not empty, as SymPy's
`not_empty_in()`. Points where an element is undefined are excluded.

## Usage

``` r
not_empty_in(finset_intersection, ...)
```

## Arguments

- finset_intersection:

  A `FiniteSet`, or its intersection with an interval or a union of
  intervals, or a union of such sets.

- ...:

  One symbol.

## Value

A set.

## Examples

``` r
x <- Symbol("x")
not_empty_in(Intersection(FiniteSet(x/2), Interval(0, 1)), x)
#> Interval(0, 2)
not_empty_in(Intersection(FiniteSet(x^2/(x + 2)), Interval(1, Inf)), x)
#> Union(Interval.Lopen(-2, -1), Interval(2, oo))
```
