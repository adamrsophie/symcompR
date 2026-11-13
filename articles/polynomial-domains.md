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
#> 1/3*x
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

## Domain matrices

A `DomainMatrix` is a matrix whose entries are elements of a domain, the
fast matrix type that SymPy uses internally.
[`DM()`](https://adamrsophie.github.io/symcompR/reference/DomainMatrix.md)
builds one from a list of rows; the other constructors build it from
flat lists, dictionaries, SymPy matrices or expressions, or as special
matrices.

``` r

A <- DM(list(list(1, 2), list(3, 4)), ZZ)
A
#> DomainMatrix([[1, 2], [3, 4]], (2, 2), ZZ)
A$det()
#> Big Integer ('bigz') :
#> [1] -2
A$rank()
#> [1] 2
A$rref()
#> (DomainMatrix([[1, 0], [0, 1]], (2, 2), QQ), (0, 1))
A$charpoly()
#> [[1]]
#> Big Integer ('bigz') :
#> [1] 1
#> 
#> [[2]]
#> Big Integer ('bigz') :
#> [1] -5
#> 
#> [[3]]
#> Big Integer ('bigz') :
#> [1] -2
A$transpose()
#> DomainMatrix([[1, 3], [2, 4]], (2, 2), ZZ)
A$convert_to(QQ)$inv()
#> DomainMatrix([[-2, 1], [3/2, -1/2]], (2, 2), QQ)
A * A
#> DomainMatrix([[7, 10], [15, 22]], (2, 2), ZZ)
A$to_Matrix()
#> Matrix(2, 2)
#> [1, 2]
#> [3, 4]
DomainMatrix(list(list(ZZ(1), ZZ(2)), list(ZZ(3), ZZ(4))), c(2, 2), ZZ)
#> DomainMatrix([[1, 2], [3, 4]], (2, 2), ZZ)
DomainMatrix_from_list(list(list(1, 2), list(3, 4)), ZZ)
#> DomainMatrix([[1, 2], [3, 4]], (2, 2), ZZ)
DomainMatrix_from_list_flat(list(ZZ(1), ZZ(2), ZZ(3), ZZ(4)), c(2, 2), ZZ)
#> DomainMatrix([[1, 2], [3, 4]], (2, 2), ZZ)
DomainMatrix_from_dod(list(`0` = list(`0` = ZZ(5)), `1` = list(`1` = ZZ(7))), c(2, 2), ZZ)
#> DomainMatrix({0: {0: 5}, 1: {1: 7}}, (2, 2), ZZ)
DomainMatrix_from_dok(list(`0,1` = ZZ(3)), c(2, 2), ZZ)
#> DomainMatrix({0: {1: 3}}, (2, 2), ZZ)
DomainMatrix_from_Matrix(Matrix(list(list(1, x), list(0, 2))))
#> DomainMatrix({0: {0: 1, 1: x}, 1: {1: 2}}, (2, 2), ZZ[x])
DomainMatrix_from_list_sympy(1, 2, list(list(Rational(1, 2), x)))
#> DomainMatrix([[1/2, x]], (1, 2), QQ[x])
DomainMatrix_from_dict_sympy(2, 2, list(`0` = list(`0` = x), `1` = list(`1` = Rational(1, 2))))
#> DomainMatrix({0: {0: x}, 1: {1: 1/2}}, (2, 2), QQ[x])
DomainMatrix_diag(list(ZZ(1), ZZ(2)), ZZ)
#> DomainMatrix({0: {0: 1}, 1: {1: 2}}, (2, 2), ZZ)
DomainMatrix_eye(2, QQ)
#> DomainMatrix({0: {0: 1}, 1: {1: 1}}, (2, 2), QQ)
DomainMatrix_ones(c(2, 3), ZZ)
#> DomainMatrix([[1, 1, 1], [1, 1, 1]], (2, 3), ZZ)
DomainMatrix_zeros(c(2, 2), ZZ)
#> DomainMatrix({}, (2, 2), ZZ)
```

The dense (`DDM`) and sparse (`SDM`) representations can be used
directly, and
[`DomainMatrix_from_rep()`](https://adamrsophie.github.io/symcompR/reference/DomainMatrix.md)
wraps one:

``` r

D <- DDM(list(list(ZZ(1), ZZ(2)), list(ZZ(3), ZZ(4))), c(2, 2), ZZ)
D
#> DDM([[1, 2], [3, 4]], (2, 2), ZZ)
D$det()
#> Big Integer ('bigz') :
#> [1] -2
D$charpoly()
#> [[1]]
#> Big Integer ('bigz') :
#> [1] 1
#> 
#> [[2]]
#> Big Integer ('bigz') :
#> [1] -5
#> 
#> [[3]]
#> Big Integer ('bigz') :
#> [1] -2
Sp <- SDM(list(`0` = list(`1` = ZZ(2))), c(2, 2), ZZ)
Sp
#> SDM({0: {1: 2}}, (2, 2), ZZ)
Sp$to_ddm()
#> DDM([[0, 2], [0, 0]], (2, 2), ZZ)
DomainMatrix_from_rep(Sp)
#> DomainMatrix({0: {1: 2}}, (2, 2), ZZ)
```

The Smith normal form of an integer matrix, with the unimodular matrices
that produce it:

``` r

smith_normal_decomp(DM(list(list(12, 6, 4), list(3, 9, 6), list(2, 16, 14)), ZZ))
#> (DomainMatrix([[1, 0, 0], [0, 10, 0], [0, 0, 30]], (3, 3), ZZ), DomainMatrix([[0, 1, -1], [1, -4, 0], [0, -2, 3]], (3, 3), ZZ), DomainMatrix([[1, 1, 10], [0, -1, -2], [0, 1, 3]], (3, 3), ZZ))
is_smith_normal_form(DM(list(list(1, 0), list(0, 2)), ZZ))
#> [1] TRUE
is_smith_normal_form(DM(list(list(2, 0), list(0, 3)), ZZ))
#> [1] FALSE
```

A `DomainScalar` is a domain element that knows its domain, so that
arithmetic unifies the domains:

``` r

DomainScalar(ZZ(3), ZZ) + DomainScalar(QQ(1, 2), QQ)
#> MPQ(7,2)
DomainScalar(ZZ(3), ZZ)$to_sympy()
#> 3
DomainScalar_from_sympy(Rational(3, 4))
#> MPQ(3,4)
DomainScalar_from_sympy(Rational(3, 4))$domain
#> QQ
```

## Power series rings

[`PythonPowerSeriesRingZZ()`](https://adamrsophie.github.io/symcompR/reference/power_series_rings.md)
and
[`PythonPowerSeriesRingQQ()`](https://adamrsophie.github.io/symcompR/reference/power_series_rings.md)
are rings of truncated univariate power series with integer and rational
coefficients:

``` r

Rz <- PythonPowerSeriesRingZZ()
Rz$pretty(Rz$pow_int(Rz(c(1, 2, 3)), 2))
#> [1] "1 + 4*x + 10*x**2 + 12*x**3 + 9*x**4"
Rz$pretty(Rz$inverse(Rz(c(1, 1))))
#> [1] "1 - x + x**2 - x**3 + x**4 - x**5 + O(x**6)"
Rq <- PythonPowerSeriesRingQQ(4)
Rq$pretty(Rq$exp(Rq$gen))
#> [1] "1 + x + 1/2*x**2 + 1/6*x**3 + O(x**4)"
Rq$pretty(Rq$log(Rq(c(1, 1))))
#> [1] "x - 1/2*x**2 + 1/3*x**3 + O(x**4)"
```
