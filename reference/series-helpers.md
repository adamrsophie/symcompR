# Series and summation helpers

`limit_seq(e, n)` the limit as `n` tends to infinity,
`difference_delta(e, n, step)` the forward difference, and
`residue(e, x, x0)` the residue at a pole (the coefficient of
`1/(x - x0)` in the Laurent series). `hypersimp(f, k)` returns the
rational term ratio `f(k + 1)/f(k)` of a hypergeometric term (or
`NULL`), and `hypersimilar(f, g, k)` tests whether `f/g` is rational.

## Usage

``` r
difference_delta(e, n = NULL, step = 1)

residue(e, x, x0 = 0)

hypersimp(f, k)

hypersimilar(f, g, k)
```

## Arguments

- e, f, g:

  Expressions.

- step:

  Step size.

- x, n, k:

  Variables.

- x0:

  Point.

## Value

An expression, `NULL` or a logical.

## Examples

``` r
x <- Symbol("x")
residue(1 / (x^2 + 1), x, sym_I)
#> -I/2
hypersimp(sym_factorial(Symbol("k")), Symbol("k"))
#> k + 1
```
