# Polynomial aliases and constructions

`sym_poly()` is
[`Poly()`](https://adamrsophie.github.io/symcompR/reference/Poly.md),
`PurePoly()` is
[`Poly()`](https://adamrsophie.github.io/symcompR/reference/Poly.md),
`GroebnerBasis()` is
[`groebner()`](https://adamrsophie.github.io/symcompR/reference/groebner.md),
and `rootof()`, `RootOf()` and `ComplexRootOf()` are
[`CRootOf()`](https://adamrsophie.github.io/symcompR/reference/CRootOf.md).
`monic()` divides by the leading coefficient, `cofactors()` returns the
gcd and both cofactors, `half_gcdex()` the Bezout coefficient and gcd,
`pexquo()` exact pseudo division, `all_roots()` all roots with
multiplicity including `CRootOf` ones. `cyclotomic_poly()`,
`symmetric_poly()`, `interpolating_poly()`, `random_poly()`,
`itermonomials()` and `viete()` build polynomials and formulas;
`legendre_poly()`, `chebyshevt_poly()`, `chebyshevu_poly()`,
`hermite_poly()`, `hermite_prob_poly()`, `laguerre_poly()` and
`jacobi_poly()` return the explicit orthogonal polynomials. `lex`,
`grlex` and `grevlex` name the monomial orders.

## Usage

``` r
PurePoly(f, ...)

GroebnerBasis(f, ...)

rootof(f, x, k)

RootOf(f, x, k)

ComplexRootOf(f, x, k)

lex

grlex

grevlex

monic(f, ...)

cofactors(f, g, ...)

half_gcdex(f, g, ...)

pexquo(f, g, ...)

all_roots(f, x = NULL)

cyclotomic_poly(n, x = Symbol("x"))

symmetric_poly(k, vars)

interpolating_poly(n, x, X = NULL, Y = NULL)

random_poly(x, n, inf = -10, sup = 10)

itermonomials(vars, max_degree, min_degree = 0)

viete(f, x = NULL, roots = NULL)

legendre_poly(n, x = Symbol("x"))

chebyshevt_poly(n, x = Symbol("x"))

chebyshevu_poly(n, x = Symbol("x"))

hermite_poly(n, x = Symbol("x"))

hermite_prob_poly(n, x = Symbol("x"))

laguerre_poly(n, x = Symbol("x"))

jacobi_poly(n, a, b, x = Symbol("x"))

sym_poly(x, ...)
```

## Arguments

- ...:

  Generators.

- x, f, g:

  Expressions.

- k:

  Order of the symmetric polynomial.

- n:

  Degree, or number of variables.

- vars:

  A list of symbols.

- X, Y:

  Interpolation nodes and values.

- inf, sup:

  Coefficient range for `random_poly()`.

- max_degree, min_degree:

  Degree bounds.

- roots:

  A list of root symbols for `viete()`.

- a, b:

  Jacobi parameters.

## Value

An expression, list or `Poly`.

## Examples

``` r
x <- Symbol("x")
cyclotomic_poly(6, x)
#> x^2 - x + 1
cofactors(x^2 - 1, x - 1)
#> $h
#> x - 1
#> 
#> $cff
#> x + 1
#> 
#> $cfg
#> 1
#> 
viete(x^2 + Symbol("p") * x + Symbol("q"), x)
#> [[1]]
#> [[1]][[1]]
#> r1 + r2
#> 
#> [[1]][[2]]
#> -p
#> 
#> 
#> [[2]]
#> [[2]][[1]]
#> r1*r2
#> 
#> [[2]][[2]]
#> q
#> 
#> 
```
