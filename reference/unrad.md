# Remove radicals from an equation

`unrad(eq, ...)` returns an equation free of the radicals in the symbols
`...` (all free symbols when none are given) whose roots include those
of `eq`, as SymPy's `unrad()`. The result is the pair `(eq, cov)`, where
`cov` is empty or holds a change of variables `p` and the equation
relating `p` to the original symbol; it is `NULL` when there are no
radicals to remove. An error of class `NotImplementedError` is raised
when the radicals cannot all be removed.

## Usage

``` r
unrad(eq, ...)
```

## Arguments

- eq:

  An expression equal to zero, or an
  [`Eq()`](https://adamrsophie.github.io/symcompR/reference/Eq.md).

- ...:

  Symbols of interest.

## Value

A two element `pytuple` or `NULL`.

## Examples

``` r
x <- Symbol("x")
unrad(sqrt(x) + sqrt(1 - x))
#> (2*x - 1, []) 
unrad(x^Rational(5, 6) + 2)
#> (x^5 - 64, []) 
```
