# More polynomial functions

`sym_decompose(f)` (SymPy's
[`decompose()`](https://rdrr.io/r/stats/decompose.html); the name
clashes with
[`stats::decompose()`](https://rdrr.io/r/stats/decompose.html)) writes
`f` as a composition `g1(g2(...))` and returns the list of components.
`gff_list(f)` is the greatest factorial factorization, a list of
`list(factor, k)`. `subresultants(f, g)` returns the subresultant
polynomial remainder sequence. `sym_trunc(f, p)` reduces coefficients to
symmetric residues modulo `p`. `symmetrize(f, ...)` writes a symmetric
polynomial in elementary symmetric polynomials and returns
`list(symmetric, remainder)` (with `formal = TRUE` also the symbols
`s1`, `s2`, ... and their definitions).
`rational_interpolate(data, degnum, x)` fits a rational function with
numerator degree `degnum` through the points. `minimal_polynomial(e, x)`
(alias `minpoly()`) is the minimal polynomial of an algebraic number
built from rationals, `I`, radicals, sums, products and `cos`/`sin` of
rational multiples of `pi`. `swinnerton_dyer_poly(n, x)` has the roots
`+-sqrt(2) +- ... +- sqrt(p_n)`. `nth_power_roots_poly(f, n)` has the
`n`-th powers of the roots of `f`. `intervals(f)` isolates the real
roots in rational intervals (refined to width `eps`, within `inf` and
`sup`) by continued fractions as SymPy,
`refine_root(f, s, t, eps, steps)` narrows one, and `ground_roots(f)`
gives the rational roots with multiplicities. `hurwitz_conditions(f, x)`
and `schur_conditions(f, x)` are the stability conditions: expressions
that must all be positive. `poly_from_expr(e)` returns
`list(Poly, gens)`.

## Usage

``` r
sym_decompose(f, x = NULL)

gff_list(f, x = NULL)

subresultants(f, g, x = NULL)

sym_trunc(f, p, x = NULL)

symmetrize(f, ..., formal = FALSE)

rational_interpolate(data, degnum, x = Symbol("x"))

minimal_polynomial(e, x = Symbol("x"))

minpoly(e, x = Symbol("x"))

swinnerton_dyer_poly(n, x = Symbol("x"))

nth_power_roots_poly(f, n, x = NULL)

ground_roots(f, x = NULL)

intervals(
  f,
  x = NULL,
  eps = NULL,
  inf = NULL,
  sup = NULL,
  fast = FALSE,
  sqf = FALSE
)

refine_root(f, s, t, eps = NULL, steps = NULL, fast = FALSE, x = NULL)

hurwitz_conditions(f, x = NULL)

schur_conditions(f, x = NULL)

poly_from_expr(e, ...)
```

## Arguments

- f, g:

  Polynomials.

- x:

  The generator.

- p:

  Modulus.

- ...:

  Generators.

- formal:

  Return symbols for the elementary polynomials.

- data:

  A list of `c(x, y)` points.

- degnum:

  Numerator degree.

- e:

  An algebraic number.

- n:

  Integer.

- eps:

  Width of the refined interval.

- inf, sup:

  Bounds of the search interval.

- fast:

  Use fast continued fraction steps.

- sqf:

  Treat `f` as square-free and return bare intervals.

- s, t:

  Interval endpoints.

- steps:

  Number of refinement steps.

## Value

Expressions, lists or vectors.

## Examples

``` r
x <- Symbol("x"); y <- Symbol("y")
sym_decompose(x^4 + 2*x^3 - x - 1)
#> [[1]]
#> x^2 - x - 1
#> 
#> [[2]]
#> x^2 + x
#> 
subresultants(x^2 + 2*x + 1, x^2 - 1)
#> [[1]]
#> x^2 + 2*x + 1
#> 
#> [[2]]
#> x^2 - 1
#> 
#> [[3]]
#> -2*x - 2
#> 
symmetrize(x^2 + y^2)
#> [[1]]
#> -2*x*y + (x + y)^2
#> 
#> [[2]]
#> 0
#> 
minimal_polynomial(sqrt(S(2)) + sqrt(S(3)), x)
#> x^4 - 10*x^2 + 1
intervals(x^2 - 3)
#> [[1]]
#> [[1]][[1]]
#> [[1]][[1]][[1]]
#> -2
#> 
#> [[1]][[1]][[2]]
#> -1
#> 
#> 
#> [[1]][[2]]
#> [1] 1
#> 
#> 
#> [[2]]
#> [[2]][[1]]
#> [[2]][[1]][[1]]
#> 1
#> 
#> [[2]][[1]][[2]]
#> 2
#> 
#> 
#> [[2]][[2]]
#> [1] 1
#> 
#> 
```
