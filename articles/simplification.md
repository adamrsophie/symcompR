# Simplification

Symbolic results are only useful when they can be put in a convenient
form. symcompR ports SymPy’s simplification module: a general purpose
[`simplify()`](https://adamrsophie.github.io/symcompR/reference/simplify.md),
and a set of targeted functions that each perform one well defined
transformation. The targeted functions are predictable and fast;
[`simplify()`](https://adamrsophie.github.io/symcompR/reference/simplify.md)
tries many of them and keeps the shortest result.

``` r

library(symcompR)
x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c")
n <- Symbol("n", integer = TRUE); k <- Symbol("k", integer = TRUE)
```

## `simplify()`

[`simplify()`](https://adamrsophie.github.io/symcompR/reference/simplify.md)
applies a sequence of heuristics (cancellation, trigonometric, power and
gamma function simplification, and others) and returns the simplest form
it finds, measured by
[`count_ops()`](https://adamrsophie.github.io/symcompR/reference/count_ops.md).

``` r

simplify(sin(x)^2 + cos(x)^2)
#> 1
simplify((x^3 + x^2 - x - 1) / (x^2 + 2 * x + 1))
#> x - 1
simplify(gamma(x) / gamma(x - 2))
#> (x - 2)*(x - 1)
simplify(exp(x) * exp(y) / exp(x + y - 1))
#> E
simplify(2 * sin(x) * cos(x))
#> sin(2*x)
count_ops(sin(x)^2 + cos(x)^2)
#> [1] 5
```

[`simplify()`](https://adamrsophie.github.io/symcompR/reference/simplify.md)
never applies an identity that is false for some values of the symbols.
\sqrt{x^2} = x only holds for nonnegative x, so it is left alone until
the symbol carries an assumption:

``` r

simplify(sqrt(x^2))
#> sqrt(x^2)
p <- Symbol("p", positive = TRUE)
simplify(sqrt(p^2))
#> p
```

Because the result of
[`simplify()`](https://adamrsophie.github.io/symcompR/reference/simplify.md)
depends on heuristics, code that needs a particular form should call the
targeted function directly. The rest of this vignette describes those
functions.

## Expansion

[`expand()`](https://adamrsophie.github.io/symcompR/reference/expand.md)
multiplies out products and powers of sums.

``` r

expand((x + 1)^2)
#> x^2 + 2*x + 1
expand((x + 2) * (x - 3))
#> x^2 - x - 6
expand((x + y)^2 * (x - y))
#> x^3 + x^2*y - x*y^2 - y^3
expand((x + 1)^2 / (x - 1))
#> x^2/(x - 1) + 2*x/(x - 1) + 1/(x - 1)
```

It accepts hints that switch on other kinds of expansion. Each hint also
exists as a standalone function:
[`expand_log()`](https://adamrsophie.github.io/symcompR/reference/expand.md),
[`expand_trig()`](https://adamrsophie.github.io/symcompR/reference/expand.md),
[`expand_power_exp()`](https://adamrsophie.github.io/symcompR/reference/expand.md),
[`expand_power_base()`](https://adamrsophie.github.io/symcompR/reference/expand.md),
[`expand_complex()`](https://adamrsophie.github.io/symcompR/reference/expand.md),
[`expand_func()`](https://adamrsophie.github.io/symcompR/reference/expand.md),
[`expand_mul()`](https://adamrsophie.github.io/symcompR/reference/expand.md),
and
[`expand_multinomial()`](https://adamrsophie.github.io/symcompR/reference/expand.md).

Logarithms and powers of products only split when the symbols allow it,
or when `force = TRUE` tells symcompR to ignore assumptions:

``` r

expand(log(x * y), log = TRUE)
#> log(x*y)
expand_log(log(x * y), force = TRUE)
#> log(x) + log(y)
expand_log(log(a^2 * b), force = TRUE)
#> 2*log(a) + log(b)
expand_power_base((x * y)^z, force = TRUE)
#> x^z*y^z
expand_power_exp(exp(x + y))
#> exp(x)*exp(y)
```

Trigonometric, complex, and special function expansion:

``` r

expand_trig(sin(x + y))
#> sin(x)*cos(y) + sin(y)*cos(x)
expand_trig(cos(2 * x))
#> 2*cos(x)^2 - 1
expand_complex(exp(x + sym_I * y))
#> I*exp(re(x) - im(y))*sin(re(y) + im(x)) + exp(re(x) - im(y))*cos(re(y) + im(x))
expand_func(gamma(x + 3))
#> x*(x + 1)*(x + 2)*gamma(x)
expand_mul(x * (y + 1))
#> x*y + x
expand_multinomial((x + y + 1)^2)
#> x^2 + 2*x*y + 2*x + y^2 + 2*y + 1
```

## Factoring

[`factor_sym()`](https://adamrsophie.github.io/symcompR/reference/sym_factor.md)
(also
[`sym_factor()`](https://adamrsophie.github.io/symcompR/reference/sym_factor.md))
factors a polynomial into irreducible factors over the rationals. It
works for multivariate polynomials too.

``` r

factor_sym(x^3 - x^2 + x - 1)
#> (x - 1)*(x^2 + 1)
factor_sym(x^2 * z + 4 * x * y * z + 4 * y^2 * z)
#> z*(x + 2*y)^2
factor_sym(x^4 - 4)
#> (x^2 - 2)*(x^2 + 2)
```

Other domains are selected with options: an algebraic extension, the
Gaussian integers, or a finite field.

``` r

factor_sym(x^4 - 4, extension = sqrt(S(2)))
#> (x - sqrt(2))*(x + sqrt(2))*(x^2 + 2)
factor_sym(x^2 + 1, gaussian = TRUE)
#> (x - I)*(x + I)
factor_sym(x^4 + 1, modulus = 2)
#> (x + 1)^4
```

[`factor_list()`](https://adamrsophie.github.io/symcompR/reference/sym_factor.md)
returns the content and a list of factors, each paired with its
multiplicity.
[`sqf()`](https://adamrsophie.github.io/symcompR/reference/sym_factor.md)
computes only the square free decomposition, which is much cheaper.

``` r

fl <- factor_list(x^5 - x^3)
fl$coeff
#> 1
vapply(fl$factors, function(f) paste0("(", sstr(f[[1]]), ")^", f[[2]]), "")
#> [1] "(x - 1)^1" "(x + 1)^1" "(x)^3"
sqf(x^5 - x^3)
#> x^3*(x^2 - 1)
```

[`factor_terms()`](https://adamrsophie.github.io/symcompR/reference/factor_terms.md)
pulls out common factors without factoring the remaining polynomial, and
[`collect_const()`](https://adamrsophie.github.io/symcompR/reference/collect_const.md)
groups numerical coefficients.

``` r

factor_terms(3 * x + 3 * y)
#> 3*(x + y)
factor_terms(x * (x + x * y) + x^2)
#> x^2*(y + 2)
collect_const(2 * x + 2 * y + 3 * z)
#> 3*z + 2*(x + y)
```

## Collecting terms

[`collect()`](https://adamrsophie.github.io/symcompR/reference/collect.md)
groups the terms of a sum by powers of an expression.
[`coeff()`](https://adamrsophie.github.io/symcompR/reference/degree.md)
extracts one coefficient.

``` r

e <- x * y + x - 3 + 2 * x^2 - z * x^2 + x^3
collect(e, x)
#> x^3 + x^2*(2 - z) + x*(y + 1) - 3
coeff(e, x, 2)
#> 2 - z
collect(a * x^2 + b * x^2 + a * x - b * x + c, x)
#> c + x^2*(a + b) + x*(a - b)
collect(a * sin(2 * x) + b * sin(2 * x), sin(2 * x))
#> (a + b)*sin(2*x)
rcollect((x + y * x) * (1 + x) + x, x)
#> x*((x + 1)*(y + 1) + 1)
```

## Rational functions

[`cancel()`](https://adamrsophie.github.io/symcompR/reference/cancel.md)
puts a rational function in the canonical form p/q, with p and q
expanded polynomials without common factors.
[`together()`](https://adamrsophie.github.io/symcompR/reference/together.md)
combines terms over a common denominator without expanding, and
[`ratsimp()`](https://adamrsophie.github.io/symcompR/reference/simplify-helpers.md)
does both.

``` r

cancel((x^2 + 2 * x + 1) / (x^2 + x))
#> (x + 1)/x
cancel(1 / x + (3 * x / 2 - 2) / (x - 4))
#> (3*x^2 - 2*x - 8)/(2*x^2 - 8*x)
cancel((x * y^2 - 2 * x * y * z + x * z^2 + y^2 - 2 * y * z + z^2) / (x^2 - 1))
#> (y^2 - 2*y*z + z^2)/(x - 1)
together(1 / x + 1 / y + 1 / z)
#> (x*y + x*z + y*z)/(x*y*z)
together(1 / (x + 1) + 1 / (x - 1))
#> 2*x/((x - 1)*(x + 1))
ratsimp(1 / x + 1 / y)
#> (x + y)/(x*y)
fraction((x + 1) / x^2)
#> [[1]]
#> x + 1
#> 
#> [[2]]
#> x^2
```

[`apart()`](https://adamrsophie.github.io/symcompR/reference/apart.md)
performs the partial fraction decomposition, the inverse of
[`together()`](https://adamrsophie.github.io/symcompR/reference/together.md).
With `full = TRUE` it uses the Bronstein algorithm and splits over the
roots of the denominator;
[`doit()`](https://adamrsophie.github.io/symcompR/reference/doit.md)
writes the result out explicitly.

``` r

apart((4 * x^3 + 21 * x^2 + 10 * x + 12) / (x^4 + 5 * x^3 + 5 * x^2 + 4 * x))
#> (2*x - 1)/(x^2 + x + 1) - 1/(x + 4) + 3/x
apart(1 / (x^2 + 1) / (x - 1))
#> -(x + 1)/(2*(x^2 + 1)) + 1/(2*(x - 1))
doit(apart(1 / (x^3 - 2), full = TRUE))
#> (-2^(1/3)/2 - 2^(1/3)*sqrt(3)*I/2)/(6*(x + 2^(1/3)/2 + 2^(1/3)*sqrt(3)*I/2)) + (-2^(1/3)/2 + 2^(1/3)*sqrt(3)*I/2)/(6*(x + 2^(1/3)/2 - 2^(1/3)*sqrt(3)*I/2)) + 2^(1/3)/(6*(x - 2^(1/3)))
```

## Trigonometric simplification

[`trigsimp()`](https://adamrsophie.github.io/symcompR/reference/trigsimp.md)
simplifies expressions with trigonometric and hyperbolic functions using
the Fu et al. algorithm, and
[`exptrigsimp()`](https://adamrsophie.github.io/symcompR/reference/trigsimp.md)
converts between exponential and hyperbolic forms.

``` r

trigsimp(sin(x)^4 - 2 * cos(x)^2 * sin(x)^2 + cos(x)^4)
#> cos(4*x)/2 + 1/2
trigsimp(sin(x) * tan(x) / sec(x))
#> sin(x)^2
trigsimp(cosh(x)^2 + sinh(x)^2)
#> cosh(2*x)
exptrigsimp(exp(z) + exp(-z))
#> 2*cosh(z)
```

### `fu()` and the TR rules

[`fu()`](https://adamrsophie.github.io/symcompR/reference/fu.md) is the
full Fu, Zhong and Zeng algorithm. It searches combinations of
transformation rules for the result with the fewest trigonometric
functions.

``` r

fu(sin(S(50))^2 + cos(S(50))^2 + sin(sym_pi / 6))
#> 3/2
fu(sqrt(S(3)) * cos(x) / 2 + sin(x) / 2)
#> sin(x + pi/3)
fu(1 - sin(2 * x)^2 / 4 - sin(y)^2 - cos(x)^4)
#> -cos(x)^2 + cos(y)^2
fu(cos(4 * sym_pi / 9))
#> sin(pi/18)
fu(cos(x) * cos(2 * x) * cos(4 * x))
#> sin(8*x)/(8*sin(x))
```

Each individual rule is exported under its name in the paper, from
[`TR0()`](https://adamrsophie.github.io/symcompR/reference/fu.md) to
[`TR22()`](https://adamrsophie.github.io/symcompR/reference/fu.md),
together with
[`TRmorrie()`](https://adamrsophie.github.io/symcompR/reference/fu.md),
[`TR10i()`](https://adamrsophie.github.io/symcompR/reference/fu.md), and
others, so a transformation can be applied by hand:

``` r

TR2(tan(x))              # tan and cot to sin and cos
#> sin(x)/cos(x)
TR5(sin(x)^4)            # sin^2 to 1 - cos^2
#> (1 - cos(x)^2)^2
TR8(cos(2 * x) * cos(3 * x))  # products to sums
#> cos(x)/2 + cos(5*x)/2
TR10(sin(a + b))         # expand sums of angles
#> sin(a)*cos(b) + sin(b)*cos(a)
TR10i(cos(S(1)) * cos(S(3)) + sin(S(1)) * sin(S(3)))  # and back
#> cos(2)
TR11(sin(2 * x))         # double angles
#> 2*sin(x)*cos(x)
TRmorrie(cos(x) * cos(2 * x))
#> sin(4*x)/(4*sin(x))
TR22(1 + tan(x)^2)
#> sec(x)^2
```

The remaining rules cover reciprocal functions, products of tangents,
squares, and powers:

``` r

TR0(x * (x + 1))                    # plain simplification step
#> x^2 + x
TR1(2 * csc(x) + sec(x))            # sec and csc to 1/cos and 1/sin
#> 1/cos(x) + 2/sin(x)
TR111(1 - 1 / tan(x)^2)             # negative powers to reciprocal functions
#> 1 - cot(x)^2
TR2i(sin(x) / cos(x))               # ratios of sin and cos to tan
#> tan(x)
TR3(cos(sym_pi / 2 + x))            # induced angles
#> -sin(x)
lapply(list(sin, cos, tan), function(f) TR4(f(sym_pi / 6)))  # special angles
#> [[1]]
#> 1/2
#> 
#> [[2]]
#> sqrt(3)/2
#> 
#> [[3]]
#> sqrt(3)/3
TR6(cos(x)^4)                       # cos^2 to 1 - sin^2
#> (1 - sin(x)^2)^2
TR7(cos(x)^2)                       # lower the degree of cos^2
#> cos(2*x)/2 + 1/2
TR9(cos(S(1)) + cos(S(2)))          # sums of cos or sin to products
#> 2*cos(1/2)*cos(3/2)
TR12(tan(x + y))                    # expand tan of a sum
#> (tan(x) + tan(y))/(-tan(x)*tan(y) + 1)
TR12i((tan(x) + tan(y)) / (1 - tan(x) * tan(y)))  # and back
#> tan(x + y)
TR13(tan(S(3)) * tan(S(2)))         # products of tan or cot
#> -tan(2)/tan(5) - tan(3)/tan(5) + 1
TR14((cos(x) - 1) * (cos(x) + 1))   # factored powers of sin and cos
#> -sin(x)^2
TR15(1 - 1 / sin(x)^2)              # negative powers of sin
#> -cot(x)^2
TR16(1 - 1 / cos(x)^2)              # negative powers of cos
#> -tan(x)^2
TRpower(sin(x)^3)                   # powers to sums of multiple angles
#> 3*sin(x)/4 - sin(3*x)/4
```

[`futrig()`](https://adamrsophie.github.io/symcompR/reference/trigsimp.md)
applies a fixed sequence of these rules,
[`trigsimp_old()`](https://adamrsophie.github.io/symcompR/reference/trigsimp.md)
is the older pattern based
[`trigsimp()`](https://adamrsophie.github.io/symcompR/reference/trigsimp.md),
and
[`trigsimp_groebner()`](https://adamrsophie.github.io/symcompR/reference/trigsimp.md)
simplifies with a Groebner basis of the trigonometric identities.
[`hyper_as_trig()`](https://adamrsophie.github.io/symcompR/reference/fu.md)
rewrites hyperbolic functions as trigonometric ones and returns a
function that undoes the change.

``` r

futrig(sin(x)^4 - cos(x)^4)
#> -cos(2*x)
futrig(sinh(x) / cosh(x))
#> tanh(x)
trigsimp_old(sin(x)^2 + cos(x)^2)
#> 1
trigsimp_groebner(sin(x) / cos(x), hints = list(tan))
#> tan(x)
h <- hyper_as_trig(sinh(x)^2 + cosh(x)^2)
h[[2]](TR8(h[[1]]))
#> cosh(2*x)
```

## Powers and logarithms

[`powsimp()`](https://adamrsophie.github.io/symcompR/reference/powsimp.md)
combines powers with the same base or the same exponent. Combining bases
is only valid for suitable assumptions, so it needs `force = TRUE` for
general symbols. `combine` selects which of the two transformations to
apply.

``` r

powsimp(x^a * x^b)
#> x^(a + b)
powsimp(x^a * y^a)
#> x^a*y^a
powsimp(x^a * y^a, force = TRUE)
#> (x*y)^a
powsimp(exp(x) * exp(y))
#> exp(x + y)
powsimp(x^a * x^b * y^a * y^b, combine = "exp")
#> x^(a + b)*y^(a + b)
powsimp(x^a * x^b * y^a * y^b, combine = "base", force = TRUE)
#> (x*y)^a*(x*y)^b
```

[`powdenest()`](https://adamrsophie.github.io/symcompR/reference/powsimp.md)
simplifies powers of powers, and
[`logcombine()`](https://adamrsophie.github.io/symcompR/reference/powsimp.md)
merges sums of logarithms.

``` r

powdenest((x^a)^b, force = TRUE)
#> x^(a*b)
powdenest(sqrt(x^2), force = TRUE)
#> x
logcombine(log(x) + log(y), force = TRUE)
#> log(x*y)
logcombine(a * log(x) + log(y) - log(z), force = TRUE)
#> log(x^a*y/z)
```

## Radicals

[`radsimp()`](https://adamrsophie.github.io/symcompR/reference/radsimp.md)
rationalizes denominators, and
[`sqrtdenest()`](https://adamrsophie.github.io/symcompR/reference/radsimp.md)
removes nested square roots when possible.

``` r

radsimp(1 / (2 + sqrt(S(2))))
#> (2 - sqrt(2))/2
radsimp(1 / (1 + sqrt(S(2)) + sqrt(S(3))))
#> (-sqrt(6) + sqrt(2) + 2)/4
radsimp(x / (sqrt(x) + 1))
#> (x^(3/2) - x)/(x - 1)
sqrtdenest(sqrt(5 + 2 * sqrt(S(6))))
#> sqrt(2) + sqrt(3)
sqrtdenest(sqrt(3 + 2 * sqrt(S(2))))
#> 1 + sqrt(2)
```

## Recognizing numbers

[`nsimplify()`](https://adamrsophie.github.io/symcompR/reference/nsimplify.md)
finds a simple exact expression that is close to a floating point
number. A list of constants can be offered as building blocks.

``` r

nsimplify(0.333333333333333)
#> 1/3
nsimplify(1.4142135623730951)
#> sqrt(2)
nsimplify(2.0943951023931957, list(sym_pi))
#> 2*pi/3
nsimplify(0.5 + 1.5 * sym_I)
#> 1/2 + 3*I/2
```

## Gamma functions and combinatorial expressions

[`gammasimp()`](https://adamrsophie.github.io/symcompR/reference/gammasimp.md)
simplifies products and ratios of gamma functions, and
[`combsimp()`](https://adamrsophie.github.io/symcompR/reference/gammasimp.md)
does the same for factorials and binomial coefficients.

``` r

gammasimp(gamma(x) / gamma(x - 3))
#> (x - 3)*(x - 2)*(x - 1)
gammasimp(gamma(x) * gamma(1 - x))
#> pi/sin(pi*x)
combsimp(sym_factorial(n) / sym_factorial(n - 3))
#> n*(n - 2)*(n - 1)
combsimp(sym_binomial(n + 1, k + 1) / sym_binomial(n, k))
#> (n + 1)/(k + 1)
```

## Special functions

[`hyperexpand()`](https://adamrsophie.github.io/symcompR/reference/hyperexpand.md)
rewrites generalized hypergeometric functions in terms of elementary or
simpler special functions, and
[`besselsimp()`](https://adamrsophie.github.io/symcompR/reference/besselsimp.md)
simplifies Bessel functions.

``` r

hyperexpand(hyper(list(1, 1), list(2), z))
#> -log(1 - z)/z
besselsimp(besselj(Rational(1, 2), z))
#> sqrt(2)*sin(z)/(sqrt(pi)*sqrt(z))
```

## Other helpers

[`signsimp()`](https://adamrsophie.github.io/symcompR/reference/signsimp.md)
moves minus signs into a canonical position,
[`separatevars()`](https://adamrsophie.github.io/symcompR/reference/separatevars.md)
writes an expression as a product of factors that each depend on fewer
symbols, and
[`cse()`](https://adamrsophie.github.io/symcompR/reference/cse.md) finds
common subexpressions, which is useful before generating code.

``` r

signsimp(exp(-(x - y)))
#> exp(-x + y)
separatevars(2 * x^2 * z * sin(y) + 2 * z * x^2)
#> 2*x^2*z*(sin(y) + 1)
res <- cse(list((x + y)^2 + sin(x + y), cos(x + y)))
res$replacements[[1]]
#> [[1]]
#> x0
#> 
#> [[2]]
#> x + y
res$reduced
#> [[1]]
#> x0^2 + sin(x0)
#> 
#> [[2]]
#> cos(x0)
```

[`cse()`](https://adamrsophie.github.io/symcompR/reference/cse.md)
returns the replacements as pairs (here x_0 = x + y) and the reduced
expressions written in terms of the new symbols.

## More specialized simplifiers

[`collect_sqrt()`](https://adamrsophie.github.io/symcompR/reference/radsimp.md)
collects terms with the same square root,
[`collect_abs()`](https://adamrsophie.github.io/symcompR/reference/simplify-helpers.md)
merges products of absolute values, and
[`split_surds()`](https://adamrsophie.github.io/symcompR/reference/radsimp.md)
and
[`rad_rationalize()`](https://adamrsophie.github.io/symcompR/reference/radsimp.md)
are the helpers behind
[`radsimp()`](https://adamrsophie.github.io/symcompR/reference/radsimp.md).

``` r

collect_sqrt(sqrt(S(2)) * x + sqrt(S(2)) * y + sqrt(S(3)))
#> sqrt(2)*(x + y) + sqrt(3)
collect_abs(Abs(x) * Abs(y) + 1 / Abs(x))
#> Abs(1/x) + Abs(x*y)
show <- function(l) vapply(l, sstr, "")
show(split_surds(3 * sqrt(S(3)) + sqrt(S(5)) / 7 + sqrt(S(6)) +
                   sqrt(S(10)) + sqrt(S(15))))
#> [1] "3"                     "sqrt(2) + sqrt(5) + 3" "sqrt(5)/7 + sqrt(10)"
show(rad_rationalize(1, sqrt(S(2)) + 1))
#> [1] "-1 + sqrt(2)" "1"
```

[`inversecombine()`](https://adamrsophie.github.io/symcompR/reference/simplify-helpers.md)
cancels a function against its inverse,
[`nc_simplify()`](https://adamrsophie.github.io/symcompR/reference/nc_simplify.md)
shortens products of noncommutative symbols,
[`sum_simplify()`](https://adamrsophie.github.io/symcompR/reference/simplify-helpers.md)
and
[`product_simplify()`](https://adamrsophie.github.io/symcompR/reference/simplify-helpers.md)
join adjacent sums and products, and
[`kroneckersimp()`](https://adamrsophie.github.io/symcompR/reference/kroneckersimp.md)
simplifies products of Kronecker deltas.

``` r

inversecombine(asin(sin(x)))
#> x
inversecombine(exp(log(x)) + log(exp(y)))
#> x + y
A <- Symbol("A", commutative = FALSE); B <- Symbol("B", commutative = FALSE)
nc_simplify(A * B * A * B * A * B)
#> (A*B)^3
sum_simplify(Sum(x, list(x, 1, n)) + Sum(x, list(x, n + 1, 2 * n)))
#> Sum(x, (x, 1, 2*n))
product_simplify(Product(x, list(x, 1, n)) * Product(x, list(x, n + 1, 2 * n)))
#> Product(x, (x, 1, 2*n))
kroneckersimp(1 + KroneckerDelta(0, k) * KroneckerDelta(1, k))
#> 1
```

[`ratsimpmodprime()`](https://adamrsophie.github.io/symcompR/reference/ratsimpmodprime.md)
simplifies a rational function modulo the ideal generated by a list of
polynomials.

``` r

ratsimpmodprime((x + y^5 + y) / (x - y), list(x * y^5 - x - y), x, y,
                order = "lex")
#> (-x^2 - x*y - x - y)/(-x^2 + x*y)
```

[`posify()`](https://adamrsophie.github.io/symcompR/reference/signsimp.md)
replaces symbols without a sign by positive dummies, which allows
simplifications that are only valid for positive values, and returns the
substitution that undoes it.

``` r

p <- posify(1 + x)
p$expr
#> _x + 1
show(p$reps[[1]])
#> [1] "_x" "x"
```

[`piecewise_simplify()`](https://adamrsophie.github.io/symcompR/reference/piecewise_simplify.md)
merges and drops branches of a `Piecewise` expression, and
[`simplify_univariate()`](https://adamrsophie.github.io/symcompR/reference/piecewise_simplify.md)
is the helper that rewrites a condition in one variable as a union of
intervals.
[`denoms()`](https://adamrsophie.github.io/symcompR/reference/denoms.md)
lists the denominators of an expression, and
[`distribute_or_over_and()`](https://adamrsophie.github.io/symcompR/reference/piecewise_simplify.md)
converts a logical expression to disjunctive normal form.

``` r

piecewise_simplify(Piecewise(list(x, Eq(x, 0)), list(0, TRUE)))
#> 0
simplify_univariate(Or(And(x > 0, x < 2), And(x > 1, x < 3)))
#> ((1 <= x) & (x <= 2)) | ((1 <= x) & (x < 3)) | ((x <= 2) & (0 < x))
show(denoms(1 / x + 1 / y + 1 / (x * y)))
#> [1] "x" "y"
distribute_or_over_and(And(Or(Not(a), b), c))
#> (b & c) | (c & ~a)
```
