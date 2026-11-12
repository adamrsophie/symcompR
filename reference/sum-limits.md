# Limits of sums and products

`sum_index(s, x)` is the position of the index `x` among the limits of a
`Sum` or `Product` (from 0, as SymPy's `index`).
`reorder_limit(s, i, j)` swaps two limits and `sym_reorder(s, ...)`
(SymPy's `reorder`) applies swaps given as pairs of positions or
indices. `change_index(s, var, trafo, newvar)` substitutes a linear
transformation of an index. `reverse_order(s, ...)` reverses the limits
at the given positions or indices (Karr's convention).
`euler_maclaurin(s, m, n, eps, eval_integral)` returns the
Euler-Maclaurin approximation of a sum with `m` explicit terms and `n`
correction terms, and an error bound. `is_convergent()` and
`is_absolutely_convergent()` apply SymPy's convergence tests; they stop
when no test decides.

## Usage

``` r
sum_index(s, x)

reorder_limit(s, i, j)

sym_reorder(s, ...)

change_index(s, var, trafo, newvar = NULL)

reverse_order(s, ...)

euler_maclaurin(s, m = 0, n = 0, eps = 0, eval_integral = TRUE)

is_convergent(s)

is_absolutely_convergent(s)

continuous_domain(f, symbol, domain)
```

## Arguments

- s:

  A `Sum` or `Product`.

- x, var:

  An index symbol.

- i, j:

  Positions of limits, from 0.

- ...:

  Pairs (lists or vectors of two) of positions or indices, or positions
  or indices to reverse.

- trafo:

  A linear expression in `var`.

- newvar:

  The new index.

- m, n:

  Numbers of explicit terms and of correction terms.

- eps:

  Tolerance.

- eval_integral:

  Evaluate the integral.

- f:

  An expression.

- symbol:

  The variable.

- domain:

  A real set.

## Value

An expression, a list or `TRUE`/`FALSE`.

## Examples

``` r
n <- Symbol("n", integer = TRUE); a <- Symbol("a"); b <- Symbol("b")
x <- Symbol("x")
reverse_order(Sum(x * n, list(n, a, b)), 0)
#> Sum(-n*x, (n, b + 1, a - 1))
change_index(Sum(x^2, list(x, a, b)), x, x + 1)
#> Sum((x - 1)^2, (x, a + 1, b + 1))
is_convergent(Sum(1 / n^2, list(n, 1, oo)))
#> [1] TRUE
```
