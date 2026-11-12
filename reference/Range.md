# Ranges of integers and disjoint unions

`Range(stop)` or `Range(start, stop, step)` is the set of integers
`start, start + step, ...` before `stop`, in SymPy's canonical form (the
stop moves to the first value past the last element, and an empty range
is `Range(0, 0, 1)`). It supports membership with
[`contains()`](https://adamrsophie.github.io/symcompR/reference/sets.md),
intersections with intervals, ranges and the integers, and
`sequence_terms(r, n)` for its first elements. `DisjointUnion(...)` tags
the elements of each set with the set's index, like SymPy.

## Usage

``` r
Range(start, stop = NULL, step = 1)

DisjointUnion(...)
```

## Arguments

- start, stop, step:

  Integers or infinities.

- ...:

  Sets.

## Value

A set.

## Examples

``` r
Range(1, 10, 3)
#> Range(1, 10, 3)
contains(Range(1, 10, 3), 7)
#> [1] TRUE
Intersection(Range(1, 10, 3), Interval(3, 8))
#> Range(4, 10, 3)
```
