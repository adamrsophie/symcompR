# Polynomial domains and rings

``` r

library(symcompR)
x <- Symbol("x"); y <- Symbol("y")
```

This vignette covers `sympy.polys.domains`, `sympy.polys.rings`,
`sympy.polys.fields`, and `sympy.polys.orderings`: the coefficient
domains behind `Poly`, and the sparse polynomial rings that the
algorithms use.

## Domains

Coefficient domains are objects: `ZZ`, `QQ`, `GF(p)`, polynomial rings
`ZZ[x]`, fraction fields, and algebraic fields.

``` r

ZZ[x]
#> ZZ[x]
GF(7)(3) * GF(7)(5)
#> 1 mod 7
frac_field(QQ, x)
#> QQ(x)
algebraic_field(QQ, sqrt(S(2)))
#> QQ<sqrt(2)>
construct_domain(list(1, Rational(1, 2)))[[1]]
#> QQ
```

[`ring()`](https://adamrsophie.github.io/symcompR/reference/rings.md)
creates a sparse polynomial ring and its generators, the low level
representation that the algorithms use.

``` r

R <- ring("x,y", ZZ)
R[[1]]
#> Polynomial ring in x, y over ZZ with lex order
X <- R[[2]]; Y <- R[[3]]
(X + Y)^2
#> x^2 + 2*x*y + y^2
```

### Domain constructors

Each SymPy domain class has a constructor. The integer and rational
domains exist in “python” and “gmpy” flavors in SymPy; here they all
give the same domain.

``` r

IntegerRing()
#> ZZ
RationalField()
#> QQ
RealField()
#> RR
ComplexField()
#> CC
RR
#> RR
CC
#> CC
ZZ_I
#> ZZ_I
QQ_I
#> QQ_I
EX
#> EX
EXRAW
#> EXRAW
ExpressionDomain()
#> EX
FiniteField(5)
#> GF(5)
FF(5)
#> GF(5)
FractionField(ZZ, x)
#> ZZ(x)
PolynomialRing(QQ, x, y)
#> QQ[x,y]
poly_ring(ZZ, x)
#> ZZ[x]
AlgebraicField(QQ, sqrt(S(3)))
#> QQ<sqrt(3)>
CyclotomicField(5)
#> QQ<zeta5>
Domain(ZZ)
#> [1] TRUE
Domain(3)
#> [1] FALSE
list(ZZ_python(), ZZ_gmpy(), GMPYIntegerRing(), PythonIntegerRing(),
     QQ_python(), QQ_gmpy(), GMPYRationalField())
#> [[1]]
#> ZZ
#> 
#> [[2]]
#> ZZ
#> 
#> [[3]]
#> ZZ
#> 
#> [[4]]
#> ZZ
#> 
#> [[5]]
#> QQ
#> 
#> [[6]]
#> QQ
#> 
#> [[7]]
#> QQ
list(FF_python(3), FF_gmpy(3), GMPYFiniteField(3), PythonFiniteField(3))
#> [[1]]
#> GF(3)
#> 
#> [[2]]
#> GF(3)
#> 
#> [[3]]
#> GF(3)
#> 
#> [[4]]
#> GF(3)
PythonRational(3, 6)
#> Big Rational ('bigq') :
#> [1] 1/2
```

[`from_sympy()`](https://adamrsophie.github.io/symcompR/reference/domains.md)
converts an expression into a domain element and
[`to_sympy()`](https://adamrsophie.github.io/symcompR/reference/domains.md)
converts it back.

``` r

a <- from_sympy(QQ, Rational(3, 4))
a
#> Big Rational ('bigq') :
#> [1] 3/4
to_sympy(QQ, a)
#> 3/4
to_sympy(GF(7), GF(7)(10))
#> 3
```

## Polynomial rings and fields

[`ring()`](https://adamrsophie.github.io/symcompR/reference/rings.md)
and
[`xring()`](https://adamrsophie.github.io/symcompR/reference/rings.md)
create a ring and its generators,
[`vring()`](https://adamrsophie.github.io/symcompR/reference/rings.md)
also assigns the generators to variables, and
[`sring()`](https://adamrsophie.github.io/symcompR/reference/rings.md)
builds the ring from expressions.
[`field()`](https://adamrsophie.github.io/symcompR/reference/rings.md),
[`xfield()`](https://adamrsophie.github.io/symcompR/reference/rings.md),
[`vfield()`](https://adamrsophie.github.io/symcompR/reference/rings.md),
and
[`sfield()`](https://adamrsophie.github.io/symcompR/reference/rings.md)
are the counterparts for rational function fields.

``` r

R <- xring("x,y", ZZ)
R[[1]]
#> Polynomial ring in x, y over ZZ with lex order
X <- R[[2]][[1]]; Y <- R[[2]][[2]]
(X + Y)^2
#> x^2 + 2*x*y + y^2
vring("u,v", QQ)
#> Polynomial ring in u, v over QQ with lex order
(u - v)^3
#> u^3 - 3*u^2*v + 3*u*v^2 - v^3
sring(list(x^2 + 2 * y, x / 3))
#> [[1]]
#> Polynomial ring in x, y over QQ with lex order
#> 
#> [[2]]
#> [[2]][[1]]
#> x^2 + 2*y
#> 
#> [[2]][[2]]
#> x/3
K <- xfield("s,t", ZZ)
K[[2]][[1]] / K[[2]][[2]]
#> s/t
vfield("p", QQ)
#> Rational function field in p over QQ with lex order
1 / (p + 1)
#> 1/(p + 1)
sfield(list(x / y, 1 / x))
#> [[1]]
#> Rational function field in x, y over ZZ with lex order
#> 
#> [[2]]
#> [[2]][[1]]
#> x/y
#> 
#> [[2]][[2]]
#> 1/x
Kf <- field("u,v", QQ)
Kf[[1]]
#> Rational function field in u, v over QQ with lex order
Kf[[2]] / Kf[[3]] + 1
#> (u + v)/v
```

## Monomials and orderings

[`Monomial()`](https://adamrsophie.github.io/symcompR/reference/Monomial.md)
is an exponent vector with optional generators. The ordering keys
[`ilex()`](https://adamrsophie.github.io/symcompR/reference/misc-helpers.md),
[`igrlex()`](https://adamrsophie.github.io/symcompR/reference/misc-helpers.md),
and
[`igrevlex()`](https://adamrsophie.github.io/symcompR/reference/misc-helpers.md)
sort monomials in the inverse of the lexicographic, graded
lexicographic, and graded reverse lexicographic orders.

``` r

m <- Monomial(c(1, 2))
m
#> Monomial((1, 2))
Monomial(c(2, 0, 1), list(x, y, Symbol("z")))
#> Monomial((2, 0, 1))
ilex(c(1, 2))
#> [1] -1 -2
igrlex(c(1, 2))
#> [1] -3 -1 -2
igrevlex(c(1, 2))
#> [1] -3  2  1
```

[`Options()`](https://adamrsophie.github.io/symcompR/reference/Options.md)
is the option manager that `Poly` uses to process generators and flags.

``` r

Options(list(x, y), list(domain = "QQ"))
#> $gens
#> $gens[[1]]
#> x
#> 
#> $gens[[2]]
#> y
#> 
#> 
#> $expand
#> [1] TRUE
#> 
#> $greedy
#> [1] TRUE
#> 
#> $auto
#> [1] TRUE
#> 
#> $strict
#> [1] FALSE
#> 
#> $frac
#> [1] FALSE
#> 
#> $formal
#> [1] FALSE
#> 
#> $polys
#> NULL
#> 
#> $include
#> [1] FALSE
#> 
#> $all
#> [1] FALSE
#> 
#> $order
#> [1] "lex"
#> 
#> $domain
#> QQ
#> 
#> attr(,"class")
#> [1] "sym_poly_options"
```
