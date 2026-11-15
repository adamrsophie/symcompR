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
`to_sympy()` and `from_sympy()` convert between elements and expressions
(with one argument, `to_sympy(m)` is the quantum matrix helper of
[`matrix_tensor_product()`](https://adamrsophie.github.io/symcompR/reference/quantum-matrixutils.md)),
and `Domain(K)` tests for a domain. The names `IntegerRing()`,
`RationalField()`, `RealField()`, `ComplexField()`,
`ExpressionDomain()`, `ExpressionRawDomain()`, `FiniteField()` and the
Python and GMPY variants of SymPy (`ZZ_python`, `QQ_gmpy`, `FF_python`,
...) return the same domains. The ground type constructors make the
element types SymPy uses: `PythonInteger()` and `GMPYInteger()` a
`bigz`, `PythonRational()` and `GMPYRational()` a `bigq`,
`SymPyInteger()`, `SymPyRational()` and `SymPyReal()` the SymPy numbers,
`PythonReal()` a double and `PythonComplex()` a complex number.

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

RealField(prec = 53, dps = NULL)

ComplexField(prec = 53, dps = NULL)

ExpressionDomain()

ExpressionRawDomain()

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

PythonRationalField()

GMPYRational(p, q = 1)

PythonInteger(x)

GMPYInteger(x)

SymPyInteger(x)

SymPyRational(p, q = 1)

SymPyReal(x, ...)

PythonReal(x)

PythonComplex(real = 0, imag = 0)

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

- dps:

  Precision in decimal digits; when given it sets `prec`.

- q:

  The denominator for `PythonRational()`.

- real, imag:

  Real and imaginary parts.

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
