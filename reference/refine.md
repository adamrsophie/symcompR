# Simplify an expression under assumptions

`refine(expr, assumptions)` rewrites `expr` using what
[`ask()`](https://adamrsophie.github.io/symcompR/reference/ask.md) can
decide under the assumptions: `Abs(x)` becomes `x` when `x` is known to
be nonnegative, `(-1)^(x + 1)` becomes `-(-1)^x` for integer `x`, a
predicate or relation that the assumptions decide becomes `TRUE` or
`FALSE`, and so on for powers,
[`exp()`](https://rdrr.io/r/base/Log.html),
[`re()`](https://adamrsophie.github.io/symcompR/reference/elementary-functions.md),
[`im()`](https://adamrsophie.github.io/symcompR/reference/elementary-functions.md),
[`arg()`](https://adamrsophie.github.io/symcompR/reference/elementary-functions.md),
[`sign()`](https://rdrr.io/r/base/sign.html),
[`atan2()`](https://rdrr.io/r/base/Trig.html),
[`sin()`](https://rdrr.io/r/base/Trig.html),
[`cos()`](https://rdrr.io/r/base/Trig.html),
[`floor()`](https://rdrr.io/r/base/Round.html),
[`ceiling()`](https://rdrr.io/r/base/Round.html),
[`Heaviside()`](https://adamrsophie.github.io/symcompR/reference/elementary-functions.md)
and elements of symmetric matrices. The arguments are refined first.

## Usage

``` r
refine(expr, assumptions = TRUE, ...)
```

## Arguments

- expr:

  An expression.

- assumptions:

  Assumptions such as `Q$positive(x)`, or `TRUE` for none.

- ...:

  Further assumptions, joined to `assumptions` with `&`.

## Value

The refined expression.

## Examples

``` r
x <- Symbol("x")
refine(sqrt(x^2), Q$real(x))
#> Abs(x)
refine(sqrt(x^2), Q$positive(x))
#> x
refine(Q$real(x), Q$positive(x))
#> True
```
