# Galois groups and algebraic number fields

``` r

library(symcompR)
x <- Symbol("x")
```

This vignette covers `sympy.polys.numberfields`, including Galois groups
of polynomials, algebraic numbers, minimal polynomials, and rings of
integers.

## Galois groups

[`galois_group()`](https://adamrsophie.github.io/symcompR/reference/galois_group.md)
computes the Galois group of an irreducible polynomial over the
rationals, up to degree six. With `by_name = TRUE` it returns the name
of the transitive group, and a flag that is `TRUE` when the group is
contained in the alternating group.

``` r

galois_group(x^4 - 2, by_name = TRUE)
#> (<S4TransitiveSubgroups.D4: 'D4'>, False)
galois_group(x^3 - 3 * x + 1, by_name = TRUE)
#> (<S3TransitiveSubgroups.A3: 'A3'>, True)
galois_group(x^5 - x + 1, by_name = TRUE)
#> (<S5TransitiveSubgroups.S5: 'S5'>, False)
group_order(galois_group(x^3 + 2)[[1]])
#> [1] 6
```

x^4 - 2 has the dihedral group of order 8, x^3 - 3x + 1 is cyclic, and
x^5 - x + 1 has the full symmetric group, so its roots cannot be written
with radicals.

### The transitive groups

The Galois group is identified among the transitive subgroups of S_n.
`S1TransitiveSubgroups` through `S6TransitiveSubgroups` list their
names, and
[`get_perm_group()`](https://adamrsophie.github.io/symcompR/reference/get_perm_group.md)
returns a concrete permutation group for a name.

``` r

S4TransitiveSubgroups
#> $C4
#> S4TransitiveSubgroups.C4
#> 
#> $V
#> S4TransitiveSubgroups.V
#> 
#> $D4
#> S4TransitiveSubgroups.D4
#> 
#> $A4
#> S4TransitiveSubgroups.A4
#> 
#> $S4
#> S4TransitiveSubgroups.S4
S5TransitiveSubgroups
#> $C5
#> S5TransitiveSubgroups.C5
#> 
#> $D5
#> S5TransitiveSubgroups.D5
#> 
#> $M20
#> S5TransitiveSubgroups.M20
#> 
#> $A5
#> S5TransitiveSubgroups.A5
#> 
#> $S5
#> S5TransitiveSubgroups.S5
length(S6TransitiveSubgroups)
#> [1] 16
S1TransitiveSubgroups$S1
#> S1TransitiveSubgroups.S1
S2TransitiveSubgroups$S2
#> S2TransitiveSubgroups.S2
S3TransitiveSubgroups$A3
#> S3TransitiveSubgroups.A3
group_order(get_perm_group(S4TransitiveSubgroups$D4))
#> [1] 8
is_cyclic(get_perm_group(S4TransitiveSubgroups$V))
#> [1] FALSE
four_group()
#> PermutationGroup([
#>     (0 1)(2 3),
#>     (0 2)(1 3)])
```

The groups used for degree five and six have their own constructors:

``` r

M20()
#> PermutationGroup([
#>     (0 1 2 3 4),
#>     (1 2 4 3)])
vapply(list(M20(), PSL2F5(), PGL2F5(), S3_in_S6(), A4_in_S6(), S4m(),
            S4p(), A4xC2(), S4xC2(), G18(), G36m(), G36p(), G72()),
       group_order, 0)
#>  [1]  20  60 120   6  12  24  24  24  48  18  36  36  72
```

[`tschirnhausen_transformation()`](https://adamrsophie.github.io/symcompR/reference/tschirnhausen_transformation.md)
replaces a polynomial by another one with the same splitting field,
which the algorithm uses when a resolvent has repeated roots.

``` r

tschirnhausen_transformation(Poly(x^2 + x + 1, x))
#> [[1]]
#> Poly(x^2 + 1, x, domain='ZZ')
#> 
#> [[2]]
#> Poly(_X^2 - _X + 1, _X, domain='ZZ')
```

## Algebraic numbers

[`minimal_polynomial()`](https://adamrsophie.github.io/symcompR/reference/polys-extra.md)
finds the minimal polynomial of an algebraic number, and
[`primitive_element()`](https://adamrsophie.github.io/symcompR/reference/primitive_element.md)
finds a single generator for a field extension such as
\mathbb{Q}(\sqrt2, \sqrt3) = \mathbb{Q}(\sqrt2 + \sqrt3).

``` r

minimal_polynomial(sqrt(S(2)) + sqrt(S(3)), x)
#> x^4 - 10*x^2 + 1
minimal_polynomial(S(2)^Rational(1, 3) + 1, x)
#> x^3 - 3*x^2 + 3*x - 3
minimal_polynomial(cos(sym_pi / 7), x)
#> 8*x^3 - 4*x^2 - 4*x + 1
primitive_element(list(sqrt(S(2)), sqrt(S(3))), x)
#> [[1]]
#> x^4 - 10*x^2 + 1
#> 
#> [[2]]
#> [[2]][[1]]
#> 1
#> 
#> [[2]][[2]]
#> 1
```

[`field_isomorphism()`](https://adamrsophie.github.io/symcompR/reference/primitive_element.md)
expresses one generator in terms of another: below, \sqrt2 =
\tfrac12\theta^3 - \tfrac92\theta for \theta = \sqrt2 + \sqrt3.

``` r

field_isomorphism(sqrt(S(2)), sqrt(S(2)) + sqrt(S(3)))
#> [[1]]
#> 1/2
#> 
#> [[2]]
#> 0
#> 
#> [[3]]
#> -9/2
#> 
#> [[4]]
#> 0
to_number_field(sqrt(S(2)), sqrt(S(2)) + sqrt(S(3)))
#> sqrt(2)
```

[`AlgebraicNumber()`](https://adamrsophie.github.io/symcompR/reference/AlgebraicNumber.md)
represents an element of a number field by its coordinates in a power
basis.

``` r

an <- AlgebraicNumber(sqrt(S(2)), c(1, 1))
an
#> 1 + sqrt(2)
minimal_polynomial(an, x)
#> x^2 - 2*x - 1
```

### Rings of integers and prime decomposition

[`round_two()`](https://adamrsophie.github.io/symcompR/reference/round_two.md)
computes an integral basis of the ring of integers of a number field,
together with the field discriminant.
[`prime_decomp()`](https://adamrsophie.github.io/symcompR/reference/round_two.md)
factors a rational prime into prime ideals.

``` r

round_two(x^3 + x^2 - 2 * x + 8)
#> (Submodule[[2, 0, 0], [0, 2, 0], [0, 1, 1]]/2, -503)
prime_decomp(2, x^3 + x^2 - 2 * x + 8)
#> [[ (2, x^2 + 1) e=1, f=1 ], [ (2, (x^2 + 3*x + 2)/2) e=1, f=1 ], [ (2, (3*x^2 + 3*x)/2) e=1, f=1 ]]
```

The field defined by x^3 + x^2 - 2x + 8 is a classic example of
Dedekind: the prime 2 splits into three distinct prime ideals, but this
cannot be seen from the factorization of the polynomial modulo 2,
because 2 divides the index of \mathbb{Z}\[\theta\] in the ring of
integers.

### Working with field elements

An `AlgebraicNumber` is an element \sum c_i \theta^i of the field
\mathbb{Q}(\theta).
[`field_element()`](https://adamrsophie.github.io/symcompR/reference/AlgebraicNumber.md)
makes another element of the same field,
[`to_root()`](https://adamrsophie.github.io/symcompR/reference/AlgebraicNumber.md)
gives the number as an expression, and
[`minpoly_of_element()`](https://adamrsophie.github.io/symcompR/reference/AlgebraicNumber.md)
its minimal polynomial.
[`isolate()`](https://adamrsophie.github.io/symcompR/reference/AlgebraicNumber.md)
finds a rational isolating interval for a real algebraic number.

``` r

b <- AlgebraicNumber(sqrt(S(2)), c(1, 1))
field_element(b, c(1, 2))
#> sqrt(2) + 2
to_root(b)
#> 1 + sqrt(2)
minpoly_of_element(b, x)
#> x^2 - 2*x - 1
is_primitive_element(b)
#> [1] FALSE
to_primitive_element(b)
#> 1 + sqrt(2)
isolate(sqrt(S(2)))
#> [[1]]
#> 1
#> 
#> [[2]]
#> 2
```

[`extract_fundamental_discriminant()`](https://adamrsophie.github.io/symcompR/reference/extract_fundamental_discriminant.md)
writes a discriminant as a fundamental discriminant times a square.
[`prime_decomp()`](https://adamrsophie.github.io/symcompR/reference/round_two.md)
splits a rational prime into prime ideals, and
[`prime_valuation()`](https://adamrsophie.github.io/symcompR/reference/round_two.md)
gives the valuation of an ideal at a prime ideal; the whole ring of
integers has valuation 0.

``` r

extract_fundamental_discriminant(-432)
#> $D
#> $D$`3`
#> [1] 1
#> 
#> $D$`-1`
#> [1] 1
#> 
#> 
#> $F
#> $F$`2`
#> [1] 2
#> 
#> $F$`3`
#> [1] 1
T <- Poly(x^3 + x^2 - 2 * x + 8, x)
P <- prime_decomp(2, T)
prime_valuation(P[[1]]$ZK, P[[1]])
#> [1] 0
```
