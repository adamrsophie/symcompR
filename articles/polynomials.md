# Polynomials and algebraic numbers

The polynomial module is the engine under much of symcompR: factoring,
cancellation, partial fractions, root finding, and integration of
rational functions all reduce to polynomial arithmetic. This vignette
shows the `Poly` class, division and GCDs, factorization over different
domains, Groebner bases, resultants, roots, Galois groups, and algebraic
number theory.

``` r

library(symcompR)
x <- Symbol("x"); y <- Symbol("y")
a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c")
```

## The `Poly` class

[`Poly()`](https://adamrsophie.github.io/symcompR/reference/Poly.md)
converts an expression into a polynomial in the given generators with
coefficients in a domain. The domain is inferred: `ZZ` for integers,
`QQ` for rationals, and polynomial rings such as `ZZ[y]` when other
symbols appear in the coefficients.

``` r

p <- Poly(x^3 - 2 * x + 1, x)
p
#> Poly(x^3 - 2*x + 1, x, domain='ZZ')
Poly(x^2 + 2 * x * y + 1, x, y)
#> Poly(x^2 + 2*x*y + 1, x, y, domain='ZZ')
Poly(x^2 + y * x + 1, x)
#> Poly(x^2 + y*x + 1, x, domain='ZZ[y]')
Poly(x^2 + 1, x, modulus = 5)
#> Poly(x^2 + 1, x, modulus=5)
```

Polynomials support arithmetic, and the usual accessors return degree,
coefficients, and monomials.
[`as_expr()`](https://adamrsophie.github.io/symcompR/reference/Poly.md)
converts back to an expression.

``` r

q <- Poly(x^2 - 1, x)
p + q
#> Poly(x^3 + x^2 - 2*x, x, domain='ZZ')
p * q
#> Poly(x^5 - 3*x^3 + x^2 + 2*x - 1, x, domain='ZZ')
p^2
#> Poly(x^6 - 4*x^4 + 2*x^3 + 4*x^2 - 4*x + 1, x, domain='ZZ')
degree(p)
#> [1] 3
LC(p)
#> 1
all_coeffs(p)
#> [[1]]
#> 1
#> 
#> [[2]]
#> 0
#> 
#> [[3]]
#> -2
#> 
#> [[4]]
#> 1
coeffs(p)
#> [[1]]
#> 1
#> 
#> [[2]]
#> -2
#> 
#> [[3]]
#> 1
monoms(p)
#> [[1]]
#> [1] 3
#> 
#> [[2]]
#> [1] 1
#> 
#> [[3]]
#> [1] 0
as_expr(p)
#> x^3 - 2*x + 1
```

Most functions accept plain expressions as well and convert internally.

``` r

degree(x^3 * y + x * y^4, y)
#> [1] 4
total_degree(x^3 * y + x * y^4)
#> [1] 5
coeff(3 * x^2 * y + x, x, 2)
#> 3*y
```

## Division and greatest common divisors

[`div()`](https://adamrsophie.github.io/symcompR/reference/div.md)
returns the quotient and remainder,
[`quo()`](https://adamrsophie.github.io/symcompR/reference/div.md) and
[`rem()`](https://adamrsophie.github.io/symcompR/reference/div.md)
return one of them.
[`reduced()`](https://adamrsophie.github.io/symcompR/reference/div.md)
divides by several polynomials at once (multivariate division).

``` r

div(x^3 + 2 * x + 1, x - 1)
#> $q
#> x^2 + x + 3
#> 
#> $r
#> 4
quo(x^3 + 2 * x + 1, x - 1)
#> x^2 + x + 3
rem(x^3 + 2 * x + 1, x - 1)
#> 4
div(Poly(x^3 + 2 * x + 1, x), Poly(x^2 - 1, x))
#> $q
#> Poly(x, x, domain='ZZ')
#> 
#> $r
#> Poly(3*x + 1, x, domain='ZZ')
reduced(2 * x^4 + y^2 - x^2 + y^3, list(x^3 - x, y^3 - y))
#> $q
#> $q[[1]]
#> 2*x
#> 
#> $q[[2]]
#> 1
#> 
#> 
#> $r
#> x^2 + y^2 + y
```

[`gcd()`](https://adamrsophie.github.io/symcompR/reference/gcd.md) and
[`sym_lcm()`](https://adamrsophie.github.io/symcompR/reference/gcd.md)
work for univariate and multivariate polynomials.
[`gcdex()`](https://adamrsophie.github.io/symcompR/reference/gcd.md)
solves the Bezout identity s f + t g = h,
[`invert()`](https://adamrsophie.github.io/symcompR/reference/gcd.md)
computes an inverse modulo a polynomial, and
[`cofactors()`](https://adamrsophie.github.io/symcompR/reference/poly-aliases.md)
returns the GCD together with the cofactors.

``` r

gcd(x^2 - 1, x^2 + 2 * x + 1)
#> x + 1
gcd(x^2 * y - y, x * y^2 + y^2)
#> x*y + y
sym_lcm(x^2 - 1, x^2 + 2 * x + 1)
#> x^3 + x^2 - x - 1
gcdex(x^2 - 1, x^2 + 1)
#> $s
#> -1/2
#> 
#> $t
#> 1/2
#> 
#> $h
#> 1
invert(x^2 + 1, x^3 - 1)
#> -x^2/2 + x/2 + 1/2
cofactors(x^2 - 1, x - 1)
#> $h
#> x - 1
#> 
#> $cff
#> x + 1
#> 
#> $cfg
#> 1
```

## Factorization

[`factor_sym()`](https://adamrsophie.github.io/symcompR/reference/sym_factor.md)
factors over the rationals by default, using Zassenhaus’ algorithm for
univariate and Wang’s EEZ algorithm for multivariate polynomials.

``` r

factor_sym(x^6 - 1)
#> (x - 1)*(x + 1)*(x^2 - x + 1)*(x^2 + x + 1)
factor_sym(x^4 + 4)
#> (x^2 - 2*x + 2)*(x^2 + 2*x + 2)
factor_sym(x^2 * y^2 - x^2 - y^2 + 1)
#> (x - 1)*(x + 1)*(y - 1)*(y + 1)
```

Over an algebraic extension, the Gaussian integers, or a finite field
\mathrm{GF}(p):

``` r

factor_sym(x^2 - 2, extension = sqrt(S(2)))
#> (x - sqrt(2))*(x + sqrt(2))
factor_sym(x^4 + 1, extension = sym_I)
#> (x^2 - I)*(x^2 + I)
factor_sym(x^3 + 2 * x + 2, modulus = 3)
#> x^3 - x - 1
factor_sym(x^5 + x + 1, modulus = 7)
#> (x - 2)*(x + 3)^2*(x^2 + 3*x - 2)
```

The first finite field example is irreducible modulo 3. Coefficients
modulo p are printed in the symmetric range (-p/2, p/2\].

[`factor_list()`](https://adamrsophie.github.io/symcompR/reference/sym_factor.md)
and
[`sqf_list()`](https://adamrsophie.github.io/symcompR/reference/sym_factor.md)
return the factors with multiplicities.

``` r

str(lapply(sqf_list(x^5 - 2 * x^4 + x^3)$factors,
           function(f) list(sstr(f[[1]]), f[[2]])))
#> List of 2
#>  $ :List of 2
#>   ..$ : chr "x - 1"
#>   ..$ : num 2
#>  $ :List of 2
#>   ..$ : chr "x"
#>   ..$ : num 3
```

## Resultants and discriminants

The resultant of two polynomials vanishes exactly when they have a
common root; it eliminates a variable from a system. The discriminant
vanishes when a polynomial has a repeated root.

``` r

resultant(x^2 - 1, x - 2, x)
#> 3
resultant(x^2 + y^2 - 1, x - y, x)
#> 2*y^2 - 1
discriminant(a * x^2 + b * x + c, x)
#> -4*a*c + b^2
discriminant(x^3 + a * x + b, x)
#> -4*a^3 - 27*b^2
subresultants(x^2 + 2 * x + 1, x^2 - 1)
#> [[1]]
#> x^2 + 2*x + 1
#> 
#> [[2]]
#> x^2 - 1
#> 
#> [[3]]
#> -2*x - 2
```

## Groebner bases

[`groebner()`](https://adamrsophie.github.io/symcompR/reference/groebner.md)
computes a reduced Groebner basis with Buchberger’s algorithm. With the
lexicographic order, the last element involves only the last variable,
which is how polynomial systems are solved.

``` r

G <- groebner(list(x^2 + y^2 - 1, x - y), x, y, order = "lex")
vapply(G, sstr, "")
#> [1] "x - y"     "2*y^2 - 1"
G <- groebner(list(x * y - 1, x^2 + y^2 - 4), x, y, order = "lex")
vapply(G, sstr, "")
#> [1] "x + y^3 - 4*y"   "y^4 - 4*y^2 + 1"
G <- groebner(list(x^2 + y, x * y - 1), x, y, order = "grevlex")
vapply(G, sstr, "")
#> [1] "x^2 + y" "x*y - 1" "x + y^2"
is_zero_dimensional(list(x^2 + y^2 - 1, x - y), x, y)
#> [1] TRUE
```

## Roots

[`roots()`](https://adamrsophie.github.io/symcompR/reference/roots.md)
finds roots in closed form using formulas up to degree four and special
methods for binomials, cyclotomic polynomials, and decomposable
polynomials. It returns a list of roots with their multiplicities stored
in the `"multiplicity"` attribute (SymPy returns a dictionary).

``` r

roots(x^2 - 3 * x + 2, x)
#> [[1]]
#> 2
#> 
#> [[2]]
#> 1
#> 
#> attr(,"multiplicity")
#> [1] 1 1
r <- roots((x - 1)^2 * (x + 3), x)
attr(r, "multiplicity")
#> [1] 1 2
roots(x^3 - 1, x)
#> [[1]]
#> 1
#> 
#> [[2]]
#> -1/2 - sqrt(3)*I/2
#> 
#> [[3]]
#> -1/2 + sqrt(3)*I/2
#> 
#> attr(,"multiplicity")
#> [1] 1 1 1
roots(x^4 + 1, x)
#> [[1]]
#> -sqrt(2)/2 - sqrt(2)*I/2
#> 
#> [[2]]
#> -sqrt(2)/2 + sqrt(2)*I/2
#> 
#> [[3]]
#> sqrt(2)/2 - sqrt(2)*I/2
#> 
#> [[4]]
#> sqrt(2)/2 + sqrt(2)*I/2
#> 
#> attr(,"multiplicity")
#> [1] 1 1 1 1
```

[`all_roots()`](https://adamrsophie.github.io/symcompR/reference/poly-aliases.md)
lists every root, repeated by multiplicity, as a plain list.
[`real_roots()`](https://adamrsophie.github.io/symcompR/reference/nroots.md)
isolates the real roots exactly,
[`nroots()`](https://adamrsophie.github.io/symcompR/reference/nroots.md)
approximates all roots numerically, and
[`count_roots()`](https://adamrsophie.github.io/symcompR/reference/nroots.md)
counts the roots in an interval with Sturm sequences.

``` r

all_roots(Poly(x^4 - 1, x))
#> [[1]]
#> -1
#> 
#> [[2]]
#> 1
#> 
#> [[3]]
#> -I
#> 
#> [[4]]
#> I
real_roots(x^3 - 2 * x)
#> [[1]]
#> -sqrt(2)
#> 
#> [[2]]
#> 0
#> 
#> [[3]]
#> sqrt(2)
real_roots(x^5 - x + 1)
#> [[1]]
#> CRootOf(x^5 - x + 1, 0)
nroots(x^3 - x + 1, x)
#> [[1]]
#> -1.32471795724475
#> 
#> [[2]]
#> 0.662358978622373 - 0.562279512062301*I
#> 
#> [[3]]
#> 0.662358978622373 + 0.562279512062301*I
count_roots(x^3 - 2 * x, x)
#> [1] 3
count_roots(x^4 - 1, x, -2, 0)
#> [1] 1
sturm(x^3 - 2 * x + 1, x)
#> [[1]]
#> x^3 - 2*x + 1
#> 
#> [[2]]
#> 3*x^2 - 2
#> 
#> [[3]]
#> 4*x/3 - 1
#> 
#> [[4]]
#> 5/16
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
```

### `CRootOf` and `RootSum`

`CRootOf(f, x, k)` is the k-th root of f (real roots first, in
increasing order), an exact object that can be evaluated numerically.
[`RootSum()`](https://adamrsophie.github.io/symcompR/reference/RootSum.md)
sums a function over all roots of a polynomial without computing them.

``` r

r0 <- CRootOf(x^5 - x + 1, x, 0)
r0
#> CRootOf(x^5 - x + 1, 0)
N(r0)
#> -1.16730397826142
CRootOf(x^2 - 2, x, 1)
#> CRootOf(x^2 - 2, 1)
RootSum(x^3 + x + 1, Lambda(x, x^2))
#> -2
```

The sum of the squares of the roots of x^3 + x + 1 is e_1^2 - 2e_2 = 0 -
2 = -2, found without solving the cubic.

## Special polynomials and constructions

``` r

cyclotomic_poly(12, x)
#> x^4 - x^2 + 1
legendre_poly(4, x)
#> 35*x^4/8 - 15*x^2/4 + 3/8
chebyshevt_poly(5, x)
#> 16*x^5 - 20*x^3 + 5*x
hermite_poly(4, x)
#> 16*x^4 - 48*x^2 + 12
swinnerton_dyer_poly(2, x)
#> x^4 - 10*x^2 + 1
```

[`horner()`](https://adamrsophie.github.io/symcompR/reference/compose.md)
rewrites a polynomial for efficient evaluation,
[`compose()`](https://adamrsophie.github.io/symcompR/reference/compose.md)
and
[`sym_decompose()`](https://adamrsophie.github.io/symcompR/reference/polys-extra.md)
compose and decompose polynomials,
[`interpolate()`](https://adamrsophie.github.io/symcompR/reference/compose.md)
finds the interpolating polynomial through points, and
[`symmetrize()`](https://adamrsophie.github.io/symcompR/reference/polys-extra.md)
writes a symmetric polynomial in terms of elementary symmetric
polynomials.

``` r

horner(x^3 + 2 * x^2 + 3 * x + 4, x)
#> x*(x*(x + 2) + 3) + 4
compose(x^2 + x, x - 1, x)
#> x^2 - x
sym_decompose(x^4 + 2 * x^3 - x - 1)
#> [[1]]
#> x^2 - x - 1
#> 
#> [[2]]
#> x^2 + x
interpolate(list(c(1, 1), c(2, 4), c(3, 9)), x)
#> x^2
symmetrize(x^2 + y^2)
#> [[1]]
#> -2*x*y + (x + y)^2
#> 
#> [[2]]
#> 0
viete(x^2 + a * x + b, x)
#> [[1]]
#> [[1]][[1]]
#> r1 + r2
#> 
#> [[1]][[2]]
#> -a
#> 
#> 
#> [[2]]
#> [[2]][[1]]
#> r1*r2
#> 
#> [[2]][[2]]
#> b
monic(3 * x^2 + 6 * x + 9)
#> x^2 + 2*x + 3
```

[`symmetrize()`](https://adamrsophie.github.io/symcompR/reference/polys-extra.md)
shows x^2 + y^2 = (x + y)^2 - 2xy with a zero remainder.
[`viete()`](https://adamrsophie.github.io/symcompR/reference/poly-aliases.md)
gives Vieta’s formulas: r_1 + r_2 = -a and r_1 r_2 = b.

More classical polynomials, given by degree and variable, and the roots
of Chebyshev polynomials:

``` r

chebyshevu_poly(3, x)
#> 8*x^3 - 4*x
hermite_prob_poly(3, x)
#> x^3 - 3*x
laguerre_poly(3, x)
#> -x^3/6 + 3*x^2/2 - 3*x + 1
jacobi_poly(2, 1, 2, x)
#> 21*x^2/4 - 3*x/2 - 3/4
chebyshevt_root(3, 1)
#> 0
chebyshevu_root(3, 1)
#> 0
symmetric_poly(2, list(x, y, a))
#> a*x + a*y + x*y
```

[`interpolating_poly()`](https://adamrsophie.github.io/symcompR/reference/poly-aliases.md)
gives the Lagrange interpolating polynomial through given points,
[`itermonomials()`](https://adamrsophie.github.io/symcompR/reference/poly-aliases.md)
lists all monomials up to a total degree, and
[`random_poly()`](https://adamrsophie.github.io/symcompR/reference/poly-aliases.md)
draws integer coefficients at random (so its output depends on the R
random seed).

``` r

expand(interpolating_poly(3, x, list(1, 2, 3), list(1, 4, 9)))
#> x^2
show <- function(l) vapply(l, sstr, "")
show(sym_ordered(itermonomials(list(x, y), 2)))
#> [1] "1"   "x"   "y"   "x^2" "y^2" "x*y"
set.seed(1)
random_poly(x, 3)
#> -7*x^3 - 4*x^2 - 10*x - 9
```

## More polynomial functions

The functional interface mirrors the methods of `Poly`.
[`LM()`](https://adamrsophie.github.io/symcompR/reference/degree.md) and
[`LT()`](https://adamrsophie.github.io/symcompR/reference/degree.md) are
the leading monomial and term,
[`degree_list()`](https://adamrsophie.github.io/symcompR/reference/degree.md)
gives the degree in each generator,
[`poly_terms()`](https://adamrsophie.github.io/symcompR/reference/degree.md)
lists the terms, and
[`terms_gcd()`](https://adamrsophie.github.io/symcompR/reference/primitive.md)
pulls out the common monomial.
[`content()`](https://adamrsophie.github.io/symcompR/reference/primitive.md)
and
[`primitive()`](https://adamrsophie.github.io/symcompR/reference/primitive.md)
split off the integer content.

``` r

LM(4 * x^2 + 2 * x * y^2 + x * y + 3 * y, x, y)
#> x^2
LT(4 * x^2 + 2 * x * y^2 + x * y + 3 * y, x, y)
#> 4*x^2
degree_list(x^2 + y * x + 1)
#> [1] 2 1
vapply(poly_terms(x^2 + 2 * x * y^2 - y, x, y), function(t)
  paste0("(", paste(t$monom, collapse = ", "), "): ", sstr(t$coeff)), "")
#> [1] "(2, 0): 1"  "(1, 2): 2"  "(0, 1): -1"
terms_gcd(x^6 * y^2 + x^3 * y, x, y)
#> x^3*y*(x^3*y + 1)
content(2 * x + 4)
#> 2
show(primitive(2 * x + 4))
#>   content primitive 
#>       "2"   "x + 2"
is_polynomial(x^2 + 1, x)
#> [1] TRUE
is_polynomial(1 / x, x)
#> [1] FALSE
```

Exact and pseudo division:
[`exquo()`](https://adamrsophie.github.io/symcompR/reference/div.md)
divides exactly,
[`pdiv()`](https://adamrsophie.github.io/symcompR/reference/div.md),
[`pquo()`](https://adamrsophie.github.io/symcompR/reference/div.md), and
[`prem()`](https://adamrsophie.github.io/symcompR/reference/div.md)
divide after multiplying by a power of the leading coefficient so that
no fractions appear, and
[`pexquo()`](https://adamrsophie.github.io/symcompR/reference/poly-aliases.md)
is the exact pseudo quotient.
[`half_gcdex()`](https://adamrsophie.github.io/symcompR/reference/poly-aliases.md)
returns half of the extended Euclidean algorithm, and
[`gcd_list()`](https://adamrsophie.github.io/symcompR/reference/gcd.md)
and
[`lcm_list()`](https://adamrsophie.github.io/symcompR/reference/gcd.md)
work on lists.

``` r

exquo(x^2 - 1, x - 1)
#> x + 1
show(pdiv(x^2 + 1, 2 * x - 4))
#>         q         r 
#> "2*x + 4"      "20"
pquo(x^2 + 1, 2 * x - 4)
#> 2*x + 4
prem(x^2 + 1, 2 * x - 4)
#> 20
pexquo(x^2 - 1, 2 * x - 2)
#> 2*x + 2
show(half_gcdex(x^4 - 2 * x^3 - 6 * x^2 + 12 * x + 15, x^3 + x^2 - 4 * x - 4))
#>           s           h 
#> "3/5 - x/5"     "x + 1"
gcd_list(list(x^3 - 1, x^2 - 1, x^2 - 3 * x + 2))
#> x - 1
lcm_list(list(x^3 - 1, x^2 - 1, x^2 - 3 * x + 2))
#> x^5 - x^4 - 2*x^3 - x^2 + x + 2
```

[`gff_list()`](https://adamrsophie.github.io/symcompR/reference/polys-extra.md)
is the greatest factorial factorization,
[`sqf_norm()`](https://adamrsophie.github.io/symcompR/reference/misc-helpers.md)
the square free norm over an algebraic extension,
[`nth_power_roots_poly()`](https://adamrsophie.github.io/symcompR/reference/polys-extra.md)
a polynomial whose roots are the n-th powers of the roots of the input,
[`ground_roots()`](https://adamrsophie.github.io/symcompR/reference/polys-extra.md)
the roots in the coefficient domain with multiplicities, and
[`sym_trunc()`](https://adamrsophie.github.io/symcompR/reference/polys-extra.md)
reduces the coefficients modulo an integer into a symmetric range.

``` r

vapply(gff_list(x^5 + 2 * x^4 - x^3 - 2 * x^2), function(p)
  paste0("(", sstr(p[[1]]), ", ", p[[2]], ")"), "")
#> [1] "(x, 1)"     "(x + 2, 4)"
show(sqf_norm(x^2 - 2, sqrt(S(3)))[2:3])
#> [1] "x^2 - 2*sqrt(3)*x + 1" "x^4 - 10*x^2 + 1"
nth_power_roots_poly(x^4 - x^2 + 1, 2)
#> x^4 - 2*x^3 + 3*x^2 - 2*x + 1
minpoly(sqrt(S(2)) + sqrt(S(3)), x)
#> x^4 - 10*x^2 + 1
ground_roots(x^6 - 4 * x^4 + 4 * x^3 - x^2)
#> 0 1 
#> 2 2
sym_trunc(2 * x^3 + 3 * x^2 + 5 * x + 7, 3)
#> -x^3 - x + 1
```

[`poly_from_expr()`](https://adamrsophie.github.io/symcompR/reference/polys-extra.md)
and
[`parallel_poly_from_expr()`](https://adamrsophie.github.io/symcompR/reference/misc-helpers.md)
build `Poly` objects together with their options,
[`hurwitz_conditions()`](https://adamrsophie.github.io/symcompR/reference/polys-extra.md)
and
[`schur_conditions()`](https://adamrsophie.github.io/symcompR/reference/polys-extra.md)
give the stability conditions of a polynomial, and
[`rootof()`](https://adamrsophie.github.io/symcompR/reference/poly-aliases.md)
indexes the roots of a polynomial.

``` r

poly_from_expr(x^2 + 2 * x)[[1]]
#> Poly(x^2 + 2*x, x, domain='ZZ')
parallel_poly_from_expr(list(x^2 + 1, x * y))[[1]]
#> [[1]]
#> Poly(x^2 + 1, x, y, domain='ZZ')
#> 
#> [[2]]
#> Poly(x*y, x, y, domain='ZZ')
show(hurwitz_conditions(x^3 + 2 * x^2 + 3 * x + 1))
#> [1] "2" "5" "1"
show(schur_conditions(x^2 + x / 2 + Rational(1, 4)))
#> [1] "18" "21"
rootof(x^3 - x - 1, x, 0)
#> CRootOf(x^3 - x - 1, 0)
show(refine_root(x^2 - 3, 1, 2, eps = 1e-2))
#> [1] "19/11" "26/15"
```

[`refine_root()`](https://adamrsophie.github.io/symcompR/reference/polys-extra.md)
shrinks an isolating interval of a real root, here around \sqrt 3
\approx 1.732.
[`gff()`](https://adamrsophie.github.io/symcompR/reference/misc-helpers.md),
the greatest factorial factorization as an expression, is not
implemented in SymPy either and signals an error:

``` r

tryCatch(gff(x^5 + 2 * x^4 - x^3 - 2 * x^2), error = conditionMessage)
#> [1] "symbolic falling factorial"
```

[`apart_list()`](https://adamrsophie.github.io/symcompR/reference/apart.md)
returns the structure of a partial fraction decomposition,
[`assemble_partfrac_list()`](https://adamrsophie.github.io/symcompR/reference/apart.md)
turns it back into an expression, and
[`numer()`](https://adamrsophie.github.io/symcompR/reference/fraction.md)
and
[`denom()`](https://adamrsophie.github.io/symcompR/reference/fraction.md)
take a fraction apart.

``` r

pf <- apart_list((2 * x^3 - 2 * x) / (x^2 - 2 * x + 1))
assemble_partfrac_list(pf)
#> 2*x + 4 + 4/(x - 1)
numer((x + 1) / (x - 1))
#> x + 1
denom((x + 1) / (x - 1))
#> x - 1
```

## Square free parts, Appell sequences and interpolation

[`sqf_part()`](https://adamrsophie.github.io/symcompR/reference/sym_factor.md)
is the square free part, the product of the distinct irreducible
factors.
[`gcdex_steps()`](https://adamrsophie.github.io/symcompR/reference/gcdex_steps.md)
lists the steps of the extended Euclidean algorithm as triples
`(s, t, r)` with `s f + t g = r`.

``` r

sqf_part(x^3 - 3 * x - 2)
#> x^2 - x - 2
as_expr(sqf_part(Poly(x^3 - 3 * x - 2)))
#> x^2 - x - 2
gcdex_steps(x^4 - 2 * x^3 - 6 * x^2 + 12 * x + 15, x^3 + x^2 - 4 * x - 4)
#> [[1]]
#> [[1]][[1]]
#> 0
#> 
#> [[1]][[2]]
#> 1
#> 
#> [[1]][[3]]
#> x^3 + x^2 - 4*x - 4
#> 
#> 
#> [[2]]
#> [[2]][[1]]
#> 1
#> 
#> [[2]][[2]]
#> 3 - x
#> 
#> [[2]][[3]]
#> x^2 + 4*x + 3
#> 
#> 
#> [[3]]
#> [[3]][[1]]
#> 3 - x
#> 
#> [[3]][[2]]
#> x^2 - 6*x + 10
#> 
#> [[3]][[3]]
#> 5*x + 5
```

The Appell sequences of `sympy.polys.appellseqs`:

``` r

bernoulli_poly(4, x)
#> x^4 - 2*x^3 + x^2 - 1/30
bernoulli_c_poly(4, x)
#> x^4 - 2*x^2 + 7/15
genocchi_poly(5, x)
#> -5*x^4 + 10*x^3 - 5*x
euler_poly(3, x)
#> x^3 - 3*x^2/2 + 1/4
andre_poly(4, x)
#> x^4 - 6*x^2 + 5
```

[`rational_interpolate()`](https://adamrsophie.github.io/symcompR/reference/polys-extra.md)
fits a rational function with a numerator of the given degree to data
points:

``` r

rational_interpolate(list(c(1, 2), c(2, 3), c(3, 5), c(4, 7)), 2)
#> (2*x^2 - 3*x + 1)/(x - 1)
rational_interpolate(list(c(1, 2), c(2, 3), c(3, 5), c(4, 7)), 2, x = y)
#> (2*y^2 - 3*y + 1)/(y - 1)
```

[`PurePoly()`](https://adamrsophie.github.io/symcompR/reference/poly-aliases.md)
is a polynomial that does not depend on the names of its generators, and
[`solve_poly_system()`](https://adamrsophie.github.io/symcompR/reference/solver-helpers.md)
solves a system of polynomial equations.

``` r

PurePoly(x^2 + 1) == PurePoly(y^2 + 1)
#> [1] TRUE
as_expr(PurePoly(x^2 + 1))
#> _0^2 + 1
degree(PurePoly(x^2 + 1, x))
#> [1] 2
solve_poly_system(list(x * y - 2 * y, 2 * y^2 - x^2), x, y)
#> [[1]]
#> (0, 0) 
#> 
#> [[2]]
#> (2, -sqrt(2)) 
#> 
#> [[3]]
#> (2, sqrt(2))
```

[`PolyQuintic()`](https://adamrsophie.github.io/symcompR/reference/PolyQuintic.md)
holds the resolvents used to solve a quintic in radicals, and
[`roots_quintic()`](https://adamrsophie.github.io/symcompR/reference/PolyQuintic.md)
returns the roots, or an empty list when the quintic is not solvable in
radicals.

``` r

as_expr(PolyQuintic(x^5 - 5 * x + 12)$f20)
#> x^6 - 40*x^5 + 1000*x^4 - 20000*x^3 + 250000*x^2 - 66400000*x + 976000000
roots_quintic(x^5 - x + 1)
#> list()
```

Coefficient domains and rings are covered in
[`vignette("polynomial-domains", package = "symcompR")`](https://adamrsophie.github.io/symcompR/articles/polynomial-domains.md),
and Galois groups and algebraic number fields in
[`vignette("galois-number-fields", package = "symcompR")`](https://adamrsophie.github.io/symcompR/articles/galois-number-fields.md).
