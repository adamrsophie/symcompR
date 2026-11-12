# Accumulation bounds

`AccumBounds(min, max)` (also `AccumulationBounds()`) is SymPy's closed
real interval `<min, max>` of the values an expression can take, for
example `sin(oo)` or `limit(sin(x), x, oo)`. Arithmetic with numbers and
other bounds is interval arithmetic: `AccumBounds(1, 2) + 3` is
`AccumBounds(4, 5)` and `AccumBounds(-1, 2)^2` is `AccumBounds(0, 4)`.
Comparisons decide when the bounds do,
[`contains()`](https://adamrsophie.github.io/symcompR/reference/sets.md)
tests membership and
[`intersection()`](https://adamrsophie.github.io/symcompR/reference/Point.md)
intersects with bounds or a finite set.

## Usage

``` r
AccumBounds(min, max)

AccumulationBounds(min, max)

ab_union(a, b)
```

## Arguments

- min, max:

  The bounds; equal bounds give the value itself.

- a, b:

  Bounds; `ab_union()` is SymPy's `AccumBounds.union()` (base R has
  [`union()`](https://rdrr.io/r/base/sets.html)), and
  [`intersection()`](https://adamrsophie.github.io/symcompR/reference/Point.md)
  takes bounds or a finite set.

## Value

An `AccumulationBounds` expression, or a number.

## Examples

``` r
AccumBounds(1, 2) + 3
#> AccumBounds(4, 5)
AccumBounds(-1, 2) * AccumBounds(3, 4)
#> AccumBounds(-4, 8)
1/AccumBounds(0, 2)
#> AccumBounds(1/2, oo)
sin(AccumBounds(0, pi/2))
#> AccumBounds(0, 1.0)
x <- Symbol("x")
limit(sin(x), x, Inf)
#> AccumBounds(-1, 1)
```
