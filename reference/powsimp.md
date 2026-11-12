# Combine powers

`powsimp()` combines powers with the same base (`x^a*x^b = x^(a+b)`)
and, when valid or with `force = TRUE`, powers with the same exponent
(`x^a*y^a = (x*y)^a`), as SymPy's `powsimp()`. `powdenest()` rewrites
`(x^a)^b` as `x^(a*b)` when valid. `logcombine()` merges sums of
logarithms.

## Usage

``` r
powsimp(e, deep = FALSE, combine = "all", force = FALSE, measure = count_ops)

powdenest(e, force = FALSE, polar = FALSE)

logcombine(e, force = FALSE)
```

## Arguments

- e:

  An expression.

- deep:

  Also simplify inside function arguments, bases and exponents.

- combine:

  `"all"`, `"exp"` or `"base"`.

- force:

  Ignore assumptions.

- measure:

  The complexity measure (unused by the combination itself, it orders
  the radical bases).

- polar:

  Treat the expression as one on the Riemann surface of the logarithm
  (`powdenest()`).

## Value

An expression.

## Examples

``` r
x <- Symbol("x"); a <- Symbol("a"); b <- Symbol("b")
powsimp(x^a * x^b)
#> x^(a + b)
logcombine(log(x) + log(Symbol("y")), force = TRUE)
#> log(x*y)
```
