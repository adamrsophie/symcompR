# Set aliases

`imageset()` is the simplified image of sets under a map, as SymPy's
`imageset()`: the map is a `Lambda`, a function name such as `"cos"`, or
a symbol (or list of symbols) followed by an expression, and the sets
come after it. `Contains(x, s)` is the membership condition, and
`PowerSet(s)` the set of all subsets of a finite set.

## Usage

``` r
imageset(...)

Contains(x, s)

PowerSet(s)
```

## Arguments

- ...:

  The map, then the sets.

- x, s:

  An element and a set.

## Value

A set or boolean.

## Examples

``` r
x <- Symbol("x")
imageset(x, 2 * x, Interval(0, 1))
#> Interval(0, 2)
imageset(Lambda(x, x^2), Interval(-2, 1))
#> Interval(0, 4)
```
