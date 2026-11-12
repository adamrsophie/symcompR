# Polynomial domains

SymPy's coefficient domains. `ZZ`, `QQ`, `RR`, `CC`, `EX`, `EXRAW`,
`ZZ_I` (Gaussian integers) and `QQ_I` (Gaussian rationals) are domain
objects; `GF(p)` (alias `FF`) is the finite field of `p` elements,
`poly_ring(K, ...)` (also `K[x, y]` and `PolynomialRing()`) a polynomial
ring, `frac_field(K, ...)` (`FractionField()`) a field of rational
functions, and `algebraic_field(K, ...)` (`AlgebraicField()`) an
algebraic number field such as `QQ<sqrt(2)>`; `CyclotomicField(n)` is
`QQ<zeta_n>`. Calling a domain converts its arguments to an element:
`ZZ(3)` is a `bigz`, `QQ(1, 2)` a `bigq`, `GF(7)(10)` the residue
`3 mod 7`, and `ZZ[x](x^2 + 1)` a
[`Poly()`](https://adamrsophie.github.io/symcompR/reference/Poly.md).
`to_sympy()` and `from_sympy()` convert between elements and
expressions, and `Domain(K)` tests for a domain. The names
`IntegerRing()`, `RationalField()`, `RealField()`, `ComplexField()`,
`ExpressionDomain()`, `FiniteField()` and the Python and GMPY variants
of SymPy (`ZZ_python`, `QQ_gmpy`, `FF_python`, ...) return the same
domains.

## Usage

``` r
ZZ(...)

QQ(...)

RR(...)

CC(...)

EX(...)

EXRAW(...)

ZZ_I(...)

QQ_I(...)

GF(p)

FF(p)

poly_ring(K, ...)

frac_field(K, ...)

algebraic_field(K, ...)

CyclotomicField(n)

Domain(x)

IntegerRing()

RationalField()

RealField(prec = 53)

ComplexField(prec = 53)

ExpressionDomain()

FiniteField(p)

PolynomialRing(K, ...)

FractionField(K, ...)

AlgebraicField(K, ...)

PythonIntegerRing()

GMPYIntegerRing()

ZZ_python()

ZZ_gmpy()

GMPYRationalField()

QQ_python()

QQ_gmpy()

PythonFiniteField(p)

GMPYFiniteField(p)

FF_python(p)

FF_gmpy(p)

PythonRational(p, q = 1)

to_sympy(K, a)

from_sympy(K, a)
```

## Arguments

- ...:

  Generators, extensions or values.

- p:

  A prime modulus, or the numerator for `PythonRational()`.

- K:

  A domain.

- n:

  The order of the root of unity.

- x:

  A domain.

- prec:

  Precision in bits (kept for SymPy compatibility).

- q:

  The denominator for `PythonRational()`.

- a:

  An element or expression.

## Value

A domain, an element, or an expression.

## Examples

``` r
ZZ(3); QQ(1, 2); GF(7)(10)
#> Big Integer ('bigz') :
#> [1] 3
#> Big Rational ('bigq') :
#> [1] 1/2
#> 3 mod 7
x <- Symbol("x")
ZZ[x]
#> ZZ[x]
frac_field(QQ, x)
#> QQ(x)
algebraic_field(QQ, sqrt(Integer(2)))
#> QQ<sqrt(2)>
GF(7)(3) * GF(7)(5)
#> 1 mod 7
```
