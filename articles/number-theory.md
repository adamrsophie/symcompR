# Number theory

symcompR ports `sympy.ntheory`: primality testing and prime generation,
integer factorization, arithmetic functions, modular arithmetic and
residues, discrete logarithms, the Chinese remainder theorem, continued
fractions, and partitions. Integers of any size are handled with the
[gmp](https://CRAN.R-project.org/package=gmp) package; pass a `bigz` (or
a string) when a number does not fit in a double.

``` r

library(symcompR)
big <- function(s) gmp::as.bigz(s)
```

## Primes

[`isprime()`](https://adamrsophie.github.io/symcompR/reference/isprime.md)
is a strong probable prime test that is deterministic for numbers below
2^{64} and has no known counterexample above it.

``` r

isprime(561)
#> [1] FALSE
isprime(big(2)^61 - 1)
#> [1] TRUE
isprime(big(2)^127 - 1)
#> [1] TRUE
is_mersenne_prime(127)
#> [1] TRUE
is_gaussian_prime(3 + 2 * sym_I)
#> [1] TRUE
```

561 = 3 \cdot 11 \cdot 17 is the smallest Carmichael number: it fools
the Fermat test for every base coprime to it, but not
[`isprime()`](https://adamrsophie.github.io/symcompR/reference/isprime.md).

``` r

is_carmichael(561)
#> [1] TRUE
is_fermat_pseudoprime(561, 2)
#> [1] TRUE
```

Generating primes:

``` r

nextprime(100)
#> [1] 101
prevprime(100)
#> [1] 97
prime(1000)
#> [1] 7919
primerange(50, 100)
#>  [1] 53 59 61 67 71 73 79 83 89 97
primepi(1000)
#> 168
primorial(5)
#> [1] 2310
mersenne_prime_exponent(10)
#> [1] 89
```

## Factorization

[`factorint()`](https://adamrsophie.github.io/symcompR/reference/factorint.md)
returns the prime factorization as a named vector of exponents. It
combines trial division, Pollard’s rho and p - 1 methods, and the
elliptic curve method.

``` r

factorint(360)
#> 2 3 5 
#> 3 2 1
factorint(big(2)^64 - 1)
#>       3       5      17     257     641   65537 6700417 
#>       1       1       1       1       1       1       1
factorint(big(2)^67 - 1)
#>    193707721 761838257287 
#>            1            1
primefactors(360)
#> [1] 2 3 5
factorrat(Rational(-24, 35))
#> -1  2  3  5  7 
#>  1  3  1 -1 -1
```

2^{67} - 1 is the Mersenne number that Cole factored by hand in 1903.
The individual methods are available too:
[`pollard_rho()`](https://adamrsophie.github.io/symcompR/reference/pollard_rho.md),
[`pollard_pm1()`](https://adamrsophie.github.io/symcompR/reference/pollard_rho.md),
[`ecm()`](https://adamrsophie.github.io/symcompR/reference/ecm.md)
(elliptic curves), and
[`qs()`](https://adamrsophie.github.io/symcompR/reference/qs_factor.md)
(the quadratic sieve).

``` r

pollard_rho(8051)
#> [1] 97
pollard_pm1(21477639576571, B = 2000)
#> [1] 4410317
ecm(big("398883434337287"))
#> [1]  4009823 99476569
```

Perfect powers and multiplicities:

``` r

perfect_power(6^5)
#> [1] 6 5
integer_nthroot(big(10)^40 + 1, 2)
#> $root
#> Big Integer ('bigz') :
#> [1] 100000000000000000000
#> 
#> $exact
#> [1] FALSE
multiplicity(2, 96)
#> [1] 5
multiplicity_in_factorial(2, 10)
#> [1] 8
smoothness(2^7 * 3^2)
#> [1]   3 128
```

## Divisors and arithmetic functions

``` r

divisors(60)
#>  [1]  1  2  3  4  5  6 10 12 15 20 30 60
divisor_count(60)
#> [1] 12
proper_divisors(28)
#> [1]  1  2  4  7 14
divisor_sigma(60)
#> 168
divisor_sigma(60, 2)
#> 5460
totient(60)
#> 16
reduced_totient(60)
#> 4
mobius(30)
#> -1
mobius(12)
#> 0
primenu(60)
#> 3
primeomega(60)
#> 4
```

`totient(60)` counts the integers below 60 coprime to it;
[`reduced_totient()`](https://adamrsophie.github.io/symcompR/reference/totient.md)
is the Carmichael function, the exponent of the multiplicative group
modulo 60. These functions also accept symbols, and stay unevaluated:

``` r

n <- Symbol("n", integer = TRUE, positive = TRUE)
totient(n)
#> totient(n)
```

Perfect, abundant, and amicable numbers:

``` r

is_perfect(28)
#> [1] TRUE
abundance(12)
#> [1] 4
is_amicable(220, 284)
#> [1] TRUE
```

## Modular arithmetic

``` r

igcd(36, 48, 60)
#> [1] 12
ilcm(4, 6, 10)
#> [1] 60
mod_inverse(3, 11)
#> [1] 4
Integer(17) %% 5
#> 2
```

### Chinese remainder theorem

`crt(m, v)` finds the x with x \equiv v_i \pmod{m_i}; it returns the
solution and the modulus.
[`solve_congruence()`](https://adamrsophie.github.io/symcompR/reference/crt.md)
takes the pairs `list(v, m)` directly.

``` r

crt(c(3, 5, 7), c(2, 3, 2))
#> [1]  23 105
solve_congruence(list(2, 3), list(3, 5), list(2, 7))
#> [1]  23 105
```

This is the classic puzzle of Sunzi: a number that leaves remainders 2,
3, and 2 when divided by 3, 5, and 7 is 23 + 105k.

### Residues and roots modulo n

``` r

legendre_symbol(3, 7)
#> -1
jacobi_symbol(3, 35)
#> 1
quadratic_residues(11)
#> [1] 0 1 3 4 5 9
is_quad_residue(5, 11)
#> [1] TRUE
sqrt_mod(10, 13)
#> [1] 6
sqrt_mod(10, 13, all_roots = TRUE)
#> [1] 6 7
nthroot_mod(11, 3, 13, all_roots = TRUE)
#> numeric(0)
```

11 is not a cube modulo 13, so the last list is empty.

### Orders, primitive roots, and discrete logarithms

``` r

n_order(3, 7)
#> [1] 6
primitive_root(23)
#> [1] 5
is_primitive_root(5, 23)
#> [1] TRUE
discrete_log(41, 15, 7)
#> [1] 3
```

`discrete_log(41, 15, 7)` solves 7^k \equiv 15 \pmod{41}, giving k = 3.

## Continued fractions

[`continued_fraction()`](https://adamrsophie.github.io/symcompR/reference/continued_fraction.md)
expands a rational number or a quadratic irrational. The expansion of
\sqrt{7} is periodic: the second element is the repeating block.

``` r

continued_fraction(Rational(415, 93))
#> [1] 4 2 6 7
continued_fraction(sqrt(S(7)))
#> [[1]]
#> [1] 2
#> 
#> [[2]]
#> [1] 1 1 1 4
continued_fraction_periodic(0, 1, 7)
#> [[1]]
#> [1] 2
#> 
#> [[2]]
#> [1] 1 1 1 4
continued_fraction_reduce(c(4, 2, 6, 7))
#> 415/93
continued_fraction_iterator(sym_pi, 6)
#> [1]   3   7  15   1 292   1
continued_fraction_convergents(c(3, 7, 15, 1, 292))
#> [[1]]
#> 3
#> 
#> [[2]]
#> 22/7
#> 
#> [[3]]
#> 333/106
#> 
#> [[4]]
#> 355/113
#> 
#> [[5]]
#> 103993/33102
```

The convergents of \pi include the famous approximations 22/7 and
355/113.

[`egyptian_fraction()`](https://adamrsophie.github.io/symcompR/reference/egyptian_fraction.md)
writes a rational number as a sum of distinct unit fractions:

``` r

egyptian_fraction(Rational(4, 13))
#> [1]   4  18 468
```

so 4/13 = 1/4 + 1/18 + 1/468.

[`egypt_harmonic()`](https://adamrsophie.github.io/symcompR/reference/egypt_harmonic.md)
writes a rational number as a sum of the first terms of the harmonic
series and a remainder:

``` r

egypt_harmonic(Rational(5, 2))
#> [[1]]
#> [[1]][[1]]
#> 1
#> 
#> [[1]][[2]]
#> 2
#> 
#> [[1]][[3]]
#> 3
#> 
#> [[1]][[4]]
#> 4
#> 
#> [[1]][[5]]
#> 5
#> 
#> [[1]][[6]]
#> 6
#> 
#> 
#> [[2]]
#> 1/20
egypt_harmonic(Rational(1, 3))
#> [[1]]
#> list()
#> 
#> [[2]]
#> 1/3
```

## Partitions and digits

[`npartitions()`](https://adamrsophie.github.io/symcompR/reference/npartitions.md)
counts the partitions of an integer with the Hardy-Ramanujan-Rademacher
formula, so large arguments are fast.

``` r

npartitions(100)
#> [1] 190569292
npartitions(1000)
#> Big Integer ('bigz') :
#> [1] 24061467864032622473692149727991
```

``` r

digits(123, 10)
#> [1] 10  1  2  3
digits(255, 16)
#> [1] 16 15 15
count_digits(1122334)
#> 1 2 3 4 
#> 2 2 2 1
is_palindromic(12321)
#> [1] TRUE
pi_hex_digits(0, 20)
#> [1] "3243f6a8885a308d3131"
```

[`digits()`](https://adamrsophie.github.io/symcompR/reference/digits.md)
returns the base followed by the digits.
[`pi_hex_digits()`](https://adamrsophie.github.io/symcompR/reference/pi_hex_digits.md)
uses the Bailey-Borwein-Plouffe formula to compute hexadecimal digits of
\pi starting at any position, without computing the earlier ones.

## Integer sequences

The classical combinatorial numbers are functions that evaluate for
integer arguments and stay symbolic otherwise.

``` r

fibonacci(50)
#> 12586269025
lucas(20)
#> 15127
catalan(10)
#> 16796
bell(10)
#> 115975
bernoulli(10)
#> 5/66
euler(10)
#> -50521
harmonic(10)
#> 7381/2520
partition(100)
#> 190569292
fibonacci(n)
#> fibonacci(n)
```

## Sieves

`sieve` is a shared, growing table of primes, totients, and Möbius
values.
[`Sieve()`](https://adamrsophie.github.io/symcompR/reference/Sieve.md)
creates a fresh one.

``` r

s <- Sieve()
s
#> prime sieve (6): 2, 3, 5, ... 11, 13
#> totient sieve (6): 0, 1, 1, ... 2, 4
#> mobius sieve (6): 0, 1, -1, ... 0, -1
class(sieve)
#> [1] "Sieve"
```

## More primality tests

The probable prime tests behind
[`isprime()`](https://adamrsophie.github.io/symcompR/reference/isprime.md)
are available individually: Miller-Rabin with given bases
([`mr()`](https://adamrsophie.github.io/symcompR/reference/isprime.md)),
the Lucas tests, the Euler and Euler-Jacobi tests, and Proth’s test for
numbers k 2^n + 1.

``` r

mr(1373651, c(2, 3))
#> [1] FALSE
mr(1373653, c(2, 3))
#> [1] TRUE
is_lucas_prp(989)
#> [1] FALSE
is_strong_lucas_prp(5459)
#> [1] TRUE
is_extra_strong_lucas_prp(989)
#> [1] TRUE
is_euler_pseudoprime(341, 2)
#> [1] TRUE
is_euler_jacobi_pseudoprime(561, 2)
#> [1] TRUE
proth_test(41)
#> [1] TRUE
find_first_n_carmichaels(3)
#> [1]  561 1105 1729
find_carmichael_numbers_in_range(1, 3000)
#> [1]  561 1105 1729 2465 2821
carmichael$is_carmichael(561)
#> [1] TRUE
composite(1)
#> [1] 4
composite(1000)
#> [1] 1197
compositepi(100)
#> [1] 74
p <- randprime(1000, 2000)
isprime(p) && p >= 1000 && p < 2000
#> [1] TRUE
```

## More divisor functions

Unitary divisors d of n have \gcd(d, n/d) = 1, and the antidivisors of n
are the numbers that do not divide n by the largest possible margin.

``` r

udivisors(60)
#> [1]  1  3  4  5 12 15 20 60
udivisor_count(120)
#> [1] 8
udivisor_sigma(18, 3)
#> 6570
antidivisors(24)
#> [1]  7 16
antidivisor_count(24)
#> [1] 2
proper_divisor_count(12)
#> [1] 5
is_abundant(12)
#> [1] TRUE
is_deficient(13)
#> [1] TRUE
dra(19, 12)
#> [1] 8
drm(277777788888899, 10)
#> [1] 0
core(24)
#> [1] 6
core(64, 3)
#> [1] 1
```

`core(n, t)` is the t-th power free part of n: 24 is 2^2 \cdot 6, so its
square free part is 6.

[`smoothness_p()`](https://adamrsophie.github.io/symcompR/reference/smoothness.md)
reports, for each prime factor p of n, how smooth p - 1 (or p + 1 with
`m = 1`) is; this decides how well the Pollard p - 1 and Williams p + 1
methods work.

``` r

smoothness_p(10431, m = 1)
#> [[1]]
#> [1] 1
#> 
#> [[2]]
#> [[2]][[1]]
#> [1] 3 2 2 4
#> 
#> [[2]][[2]]
#> [1] 19  1  5  5
#> 
#> [[2]][[3]]
#> [1] 61  1 31 31
```

## More modular arithmetic

[`crt1()`](https://adamrsophie.github.io/symcompR/reference/crt.md) and
[`crt2()`](https://adamrsophie.github.io/symcompR/reference/crt.md)
split the Chinese remainder computation into a precomputation that
depends only on the moduli and a cheap second step.
[`symmetric_residue()`](https://adamrsophie.github.io/symcompR/reference/crt.md)
is the residue closest to zero,
[`sqrt_mod_iter()`](https://adamrsophie.github.io/symcompR/reference/sqrt_mod.md)
iterates over square roots, and
[`kronecker_symbol()`](https://adamrsophie.github.io/symcompR/reference/totient.md)
extends the Jacobi symbol to all integers.

``` r

pre <- crt1(c(18, 42, 6))
crt2(c(18, 42, 6), c(0, 0, 0), pre[[1]], pre[[2]], pre[[3]])
#> [1]    0 4536
symmetric_residue(5, 7)
#> [1] -2
unlist(sqrt_mod_iter(11, 43))
#> [1] 21 22
kronecker_symbol(3, 8)
#> -1
is_nthpow_residue(2, 3, 7)
#> [1] FALSE
quadratic_congruence(2, 5, 3, 7)
#> [1] 2 6
polynomial_congruence("x^6 - 2*x^5 - 35", 6125)
#> [1] 3257
binomial_mod(10, 2, 6)
#> [1] 3
```

[`cycle_length()`](https://adamrsophie.github.io/symcompR/reference/prime.md)
finds the tail length \mu and cycle length \lambda of the sequence
x\_{k+1} = f(x_k) with Brent’s algorithm:

``` r

cycle_length(function(x) (x^2 + 12) %% 17, 4)
#> $lambda
#> [1] 3
#> 
#> $mu
#> [1] 0
```

## Binomial and multinomial coefficients

``` r

binomial_coefficients_list(5)
#> [1]  1  5 10 10  5  1
binomial_coefficients(3)
#> 0,3 1,2 2,1 3,0 
#>   1   3   3   1
multinomial_coefficients(3, 2)
#> 0,0,2 0,1,1 0,2,0 1,0,1 1,1,0 2,0,0 
#>     1     2     1     2     2     1
multinomial_coefficients_iterator(2, 3)
#> 0,3 1,2 2,1 3,0 
#>   1   3   3   1
```

## Factoring with the quadratic sieve

[`qs()`](https://adamrsophie.github.io/symcompR/reference/qs_factor.md)
is the self-initializing quadratic sieve, with a bound on the factor
base and a sieve interval; it returns the factors it finds, and
[`qs_factor()`](https://adamrsophie.github.io/symcompR/reference/qs_factor.md)
returns them with their multiplicities.
[`FactorCache()`](https://adamrsophie.github.io/symcompR/reference/FactorCache.md)
creates a cache of factorizations, and `factor_cache` is the cache that
[`factorint()`](https://adamrsophie.github.io/symcompR/reference/factorint.md)
uses.

``` r

qs(10009202107, 100, 10000)
#> [1] 100043 100049
qs_factor(10009202107, 100, 10000)
#> 100043 100049 
#>      1      1
fc <- FactorCache(maxsize = 10)
class(fc)
#> [1] "FactorCache"
class(factor_cache)
#> [1] "FactorCache"
```

## Points on elliptic curves

[`EllipticCurvePoint()`](https://adamrsophie.github.io/symcompR/reference/EllipticCurve.md)
creates a point in projective coordinates.
[`point_order()`](https://adamrsophie.github.io/symcompR/reference/EllipticCurve.md)
finds the order of a point; on y^2 = x^3 + 1 the rational torsion group
is cyclic of order 6.

``` r

E <- EllipticCurve(0, 1)
point_order(E$point(2, 3))
#> [1] 6
point_order(E$point(0, 1))
#> [1] 3
point_order(E$point(-1, 0))
#> [1] 2
EllipticCurvePoint(2, 3, 1, E)
#> (2, 3)
```
