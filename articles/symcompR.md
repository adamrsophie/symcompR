# Getting started with symcompR

symcompR is a computer algebra system written in R. It is a native port
of the Python library [SymPy](https://www.sympy.org): expressions are
ordinary R objects, exact integers and rationals are stored with the
[gmp](https://CRAN.R-project.org/package=gmp) package, and every
algorithm runs in R. Results follow SymPy closely, down to the order in
which terms print, so code and answers move between the two systems with
little translation.

This vignette introduces the building blocks: symbols, exact numbers,
constants, printing, substitution, numerical evaluation, and turning
expressions into R functions. The other vignettes cover each area of the
package in depth; they are listed at the end.

``` r

library(symcompR)
```

## Symbols

A symbol is a named mathematical variable. Create one with
[`Symbol()`](https://adamrsophie.github.io/symcompR/reference/Symbol.md),
or several at once with
[`sym_symbols()`](https://adamrsophie.github.io/symcompR/reference/Symbol.md),
which returns a named list.

``` r

x <- Symbol("x")
y <- Symbol("y")
x + y + x
#> 2*x + y
xyz <- sym_symbols("a b c")
xyz$a * xyz$b^2
#> a*b^2
sym_symbols("t0:3")
#> $t0
#> t0
#> 
#> $t1
#> t1
#> 
#> $t2
#> t2
```

[`sym_var()`](https://adamrsophie.github.io/symcompR/reference/Symbol.md)
creates symbols and assigns them in the calling environment, the way
[`var()`](https://rdrr.io/r/stats/cor.html) does in SymPy.

``` r

sym_var("u v")
u^2 - v^2
#> u^2 - v^2
```

### Assumptions

Symbols can carry assumptions such as `positive`, `integer`, or `real`.
Assumptions change what simplifications are valid. In general \sqrt{x^2}
\neq x, but it is true for a positive x:

``` r

sqrt(x^2)
#> sqrt(x^2)
p <- Symbol("p", positive = TRUE)
sqrt(p^2)
#> p
n <- Symbol("n", integer = TRUE)
cos(2 * sym_pi * n)
#> 1
(-1)^(2 * n)
#> 1
```

The queries
[`is_positive()`](https://adamrsophie.github.io/symcompR/reference/assumption-queries.md),
[`is_integer()`](https://adamrsophie.github.io/symcompR/reference/assumption-queries.md),
[`is_real()`](https://adamrsophie.github.io/symcompR/reference/assumption-queries.md),
and friends follow SymPy’s three valued logic: `TRUE`, `FALSE`, or `NA`
when the answer is not known.

``` r

is_integer(n)
#> [1] TRUE
is_positive(p + 1)
#> [1] TRUE
is_real(x)
#> [1] NA
is_even(2 * n)
#> [1] TRUE
```

The vignette on logic and sets shows the second assumptions system,
built on
[`ask()`](https://adamrsophie.github.io/symcompR/reference/ask.md) and
`Q`.

## Building expressions

R’s arithmetic operators work on symbolic objects, and the usual
mathematical functions (`sin`, `exp`, `log`, `sqrt`, …) accept them. The
result is a `sym` object: a tree whose nodes are typed (`Add`, `Mul`,
`Pow`, function applications, …). Expressions are put in a canonical
form as they are built, so equal terms are combined automatically.

``` r

e <- (x + 1)^2 * y + 3 * y * (x + 1)^2
e
#> 4*y*(x + 1)^2
sym_type(e)
#> [1] "Mul"
sym_args(e)
#> [[1]]
#> 4
#> 
#> [[2]]
#> y
#> 
#> [[3]]
#> (x + 1)^2
free_symbols(sin(x) * exp(y))
#> [[1]]
#> x
#> 
#> [[2]]
#> y
```

Expressions are not expanded unless you ask; see the simplification
vignette for
[`expand()`](https://adamrsophie.github.io/symcompR/reference/expand.md),
[`factor_sym()`](https://adamrsophie.github.io/symcompR/reference/sym_factor.md),
and friends.

``` r

expand((x + y)^3)
#> x^3 + 3*x^2*y + 3*x*y^2 + y^3
```

`==` tests structural equality, as in SymPy. To build an equation use
[`Eq()`](https://adamrsophie.github.io/symcompR/reference/Eq.md);
[`Ne()`](https://adamrsophie.github.io/symcompR/reference/Eq.md),
[`Lt()`](https://adamrsophie.github.io/symcompR/reference/Eq.md),
[`Le()`](https://adamrsophie.github.io/symcompR/reference/Eq.md),
[`Gt()`](https://adamrsophie.github.io/symcompR/reference/Eq.md), and
[`Ge()`](https://adamrsophie.github.io/symcompR/reference/Eq.md) build
the other relations.

``` r

(x + 1)^2 == x^2 + 2 * x + 1
#> [1] FALSE
expand((x + 1)^2) == x^2 + 2 * x + 1
#> [1] TRUE
Eq(x^2, 4)
#> Eq(x^2, 4)
Lt(x, 3)
#> x < 3
```

Undefined functions are created with
[`Function()`](https://adamrsophie.github.io/symcompR/reference/Function.md).
They stay unevaluated and are useful in differential equations.

``` r

f <- Function("f")
f(x) + f(y)^2
#> f(x) + f(y)^2
diff(f(x^2), x)
#> 2*x*Subs(Derivative(f(_xi_1), _xi_1), _xi_1, x^2)
```

## Exact numbers

symcompR keeps numbers exact. Integers have arbitrary size and rationals
never lose precision.

``` r

Integer(2)^100
#> 1267650600228229401496703205376
Rational(1, 3) + Rational(1, 6)
#> 1/2
Integer("123456789012345678901234567890") * 3
#> 370370367037037036703703703670
sym_factorial(Integer(25))
#> 15511210043330985984000000
```

R evaluates `1/3` as a double before symcompR sees it, so the result is
inexact. Write `Rational(1, 3)`, or make one operand symbolic first with
[`S()`](https://adamrsophie.github.io/symcompR/reference/as_sym.md) (an
alias of
[`as_sym()`](https://adamrsophie.github.io/symcompR/reference/as_sym.md)),
and the division stays exact.

``` r

x + 1/3
#> x + 0.333333333333333
x + Rational(1, 3)
#> x + 1/3
S(1) / 3
#> 1/3
sqrt(S(8))
#> 2*sqrt(2)
```

The same applies to functions of plain R numbers: `sqrt(8)` and
`sin(pi/6)` are computed by base R, while `sqrt(S(8))` and
`sin(sym_pi/6)` are symbolic.

``` r

sqrt(8)
#> [1] 2.828427
sin(sym_pi / 6)
#> 1/2
cos(sym_pi / 4)
#> sqrt(2)/2
```

Floating point numbers are arbitrary precision `Float`s, as in SymPy.

``` r

Float("1.5", 30)
#> 1.50000000000000000000000000000
x / 2 + 0.5
#> x/2 + 0.5
```

## Constants

Constants whose names would clash with R carry a `sym_` prefix:
`sym_pi`, `sym_E`, and `sym_I` (the imaginary unit). Infinity is `oo`,
complex infinity is `zoo`, and “not a number” is `nan`. `EulerGamma`,
`Catalan`, and `GoldenRatio` are also available.

``` r

sym_I^2
#> -1
exp(sym_I * sym_pi)
#> -1
log(sym_E^3)
#> 3
1 / oo
#> 0
oo - oo
#> nan
GoldenRatio^2 - GoldenRatio
#> -GoldenRatio + GoldenRatio^2
N(GoldenRatio, 20)
#> 1.6180339887498948482
```

## The `sym_` prefix convention

Attaching symcompR masks no base R function. SymPy functions whose names
clash with base R carry a `sym_` prefix:
[`sym_factor()`](https://adamrsophie.github.io/symcompR/reference/sym_factor.md),
[`sym_integrate()`](https://adamrsophie.github.io/symcompR/reference/sym_integrate.md),
[`sym_symbols()`](https://adamrsophie.github.io/symcompR/reference/Symbol.md),
[`sym_var()`](https://adamrsophie.github.io/symcompR/reference/Symbol.md),
[`sym_det()`](https://adamrsophie.github.io/symcompR/reference/sym_det.md),
[`sym_factorial()`](https://adamrsophie.github.io/symcompR/reference/special-functions.md),
[`sym_binomial()`](https://adamrsophie.github.io/symcompR/reference/special-functions.md),
[`sym_beta()`](https://adamrsophie.github.io/symcompR/reference/special-functions.md),
[`sym_Mod()`](https://adamrsophie.github.io/symcompR/reference/elementary-functions.md),
[`sym_trace()`](https://adamrsophie.github.io/symcompR/reference/matrix-builders.md),
and a few others. Where base R has a generic (`diff`, `sin`, `exp`, `t`,
`solve`, …), symcompR adds an S3 method instead, so the base name works
on symbolic input. A few functions also have a second name with the
suffix `_sym`, such as
[`factor_sym()`](https://adamrsophie.github.io/symcompR/reference/sym_factor.md).

``` r

sym_factor(x^3 - x)
#> x*(x - 1)*(x + 1)
factor_sym(x^3 - x)
#> x*(x - 1)*(x + 1)
sym_binomial(x, 2)
#> binomial(x, 2)
expand_func(sym_binomial(x, 2))
#> x*(x - 1)/2
sym_Mod(x + 5, 2)
#> Mod(x + 1, 2)
Integer(17) %% 5
#> 2
```

Base R behavior is untouched for ordinary R values:

``` r

diff(c(1, 4, 9, 16))
#> [1] 3 5 7
sym_factorial(5)
#> [1] 120
```

## Printing

A `sym` object prints with SymPy’s
[`str()`](https://rdrr.io/r/utils/str.html) form, with powers shown as
`^`. The printers return strings:

- [`sstr()`](https://adamrsophie.github.io/symcompR/reference/sstr.md)
  gives the same one line form as
  [`print()`](https://rdrr.io/r/base/print.html).
- [`srepr()`](https://adamrsophie.github.io/symcompR/reference/sstr.md)
  shows the full expression tree.
- [`latex()`](https://adamrsophie.github.io/symcompR/reference/latex.md)
  gives LaTeX.
- [`pprint()`](https://adamrsophie.github.io/symcompR/reference/pretty.md)
  draws a two dimensional picture in the console.

``` r

e <- x^2 / (y + 1) + sqrt(x)
sstr(e)
#> [1] "sqrt(x) + x^2/(y + 1)"
srepr(x^2 + 1)
#> [1] "Add(Integer(1), Pow(Symbol('x'), Integer(2)))"
latex(e)
#> [1] "\\sqrt{x} + \\frac{x^{2}}{y + 1}"
pprint(e, use_unicode = FALSE)
#>           2  
#>   ___    x   
#> \/ x  + -----
#>         y + 1
pprint(Rational(1, 3) * x^3 + exp(-x))
#>  3      
#> x     -x
#> ── + ℯ  
#> 3
```

The printing vignette covers MathML, tables, and code generation for C,
Fortran, Python, Julia, and other languages.

## Substitution

[`subs()`](https://adamrsophie.github.io/symcompR/reference/subs.md)
replaces a subexpression and simplifies the result. It takes either an
old and a new value, or a list of pairs.

``` r

e <- x^2 + y
subs(e, x, 2)
#> y + 4
subs(e, list(list(x, 2), list(y, 3)))
#> 7
subs(sin(x), x, y^2)
#> sin(y^2)
subs(x^4 + x^2, x^2, y)
#> y^2 + y
```

[`xreplace()`](https://adamrsophie.github.io/symcompR/reference/subs.md)
does exact structural replacement without any mathematical matching. It
takes a named list from symbol names to new values.

``` r

xreplace(x + 2 * y, list(x = 3))
#> 2*y + 3
```

## Numerical evaluation

[`N()`](https://adamrsophie.github.io/symcompR/reference/N.md) (also
available as
[`evalf()`](https://adamrsophie.github.io/symcompR/reference/N.md))
evaluates an expression to a floating point number with any number of
significant digits. The default is 15.

``` r

N(sym_pi)
#> 3.14159265358979
N(sym_pi, 50)
#> 3.1415926535897932384626433832795028841971693993751
N(sqrt(S(2)), 30)
#> 1.41421356237309504880168872421
N(exp(sym_pi * sqrt(S(163))), 40)
#> 262537412640768743.9999999999992500725972
```

The last result is famous: e^{\pi\sqrt{163}} is within 10^{-12} of an
integer, which only shows at high precision.
[`N()`](https://adamrsophie.github.io/symcompR/reference/N.md) can
substitute values before evaluating, and leaves remaining symbols alone.

``` r

N(sin(x), subs = list(x = 1))
#> 0.841470984807897
N(x^2 + sym_pi, 10)
#> x^2 + 3.141592654
evalf(sqrt(S(2)) * x, 5)
#> 1.4142*x
```

To get an ordinary R number, use
[`as.numeric()`](https://rdrr.io/r/base/numeric.html).

``` r

as.numeric(sqrt(S(2)) + 1)
#> [1] 2.414214
```

## From expressions to R functions

[`lambdify()`](https://adamrsophie.github.io/symcompR/reference/lambdify.md)
compiles an expression into a vectorized R function, which is the fast
way to evaluate an expression many times, to plot it, or to pass it to
an optimizer.

``` r

g <- lambdify(x, sin(x)^2 + cos(x))
g(c(0, 1, 2))
#> [1] 1.000000 1.248376 0.410675
h <- lambdify(list(x, y), x * y + 1)
h(1:3, 2)
#> [1] 3 5 7
body(h)
#> 1 + x * y
```

For example, the derivative of an expression can be handed to
[`uniroot()`](https://rdrr.io/r/stats/uniroot.html):

``` r

fx <- x^3 - 2 * x - 5
dfx <- lambdify(x, diff(fx, x))
dfx(2)
#> [1] 10
uniroot(lambdify(x, fx), c(2, 3))$root
#> [1] 2.094541
```

## Strings and parsing

[`sympify()`](https://adamrsophie.github.io/symcompR/reference/sympify.md)
turns a string into an expression. It accepts both Python’s `**` and R’s
`^` for powers.

``` r

sympify("x**2 + sin(y)")
#> x^2 + sin(y)
sympify("x^2 + 2*x + 1")
#> x^2 + 2*x + 1
```

[`parse_expr()`](https://adamrsophie.github.io/symcompR/reference/parse_expr.md),
[`parse_latex()`](https://adamrsophie.github.io/symcompR/reference/parse_latex.md),
and
[`parse_mathematica()`](https://adamrsophie.github.io/symcompR/reference/parse_mathematica.md)
are described in the printing, code generation, and parsing vignette.

## Where to go next

Core algebra:

- [`vignette("expressions")`](https://adamrsophie.github.io/symcompR/articles/expressions.md):
  traversing, substituting, and constructing expression trees.
- [`vignette("simplification")`](https://adamrsophie.github.io/symcompR/articles/simplification.md):
  expansion, factoring, and simplification.
- [`vignette("assumptions")`](https://adamrsophie.github.io/symcompR/articles/assumptions.md):
  assumptions on symbols,
  [`ask()`](https://adamrsophie.github.io/symcompR/reference/ask.md),
  and
  [`refine()`](https://adamrsophie.github.io/symcompR/reference/refine.md).
- [`vignette("pattern-matching")`](https://adamrsophie.github.io/symcompR/articles/pattern-matching.md):
  pattern matching, unification, and rewriting strategies.
- [`vignette("polynomials")`](https://adamrsophie.github.io/symcompR/articles/polynomials.md):
  polynomial algebra.
- [`vignette("polynomial-domains")`](https://adamrsophie.github.io/symcompR/articles/polynomial-domains.md):
  coefficient domains and polynomial rings.
- [`vignette("matrices")`](https://adamrsophie.github.io/symcompR/articles/matrices.md):
  linear algebra and matrix expressions.
- [`vignette("logic-and-sets")`](https://adamrsophie.github.io/symcompR/articles/logic-and-sets.md):
  Boolean algebra, sets, and ordinals.

Calculus and equations:

- [`vignette("calculus")`](https://adamrsophie.github.io/symcompR/articles/calculus.md):
  derivatives, limits, series, sums, and products.
- [`vignette("integration")`](https://adamrsophie.github.io/symcompR/articles/integration.md):
  integrals and integral transforms.
- [`vignette("series-and-sequences")`](https://adamrsophie.github.io/symcompR/articles/series-and-sequences.md):
  sequences, Fourier series, and formal power series.
- [`vignette("holonomic")`](https://adamrsophie.github.io/symcompR/articles/holonomic.md):
  holonomic functions.
- [`vignette("special-functions")`](https://adamrsophie.github.io/symcompR/articles/special-functions.md):
  special functions and `Piecewise`.
- [`vignette("solvers")`](https://adamrsophie.github.io/symcompR/articles/solvers.md):
  algebraic, differential, and Diophantine equations.

Discrete mathematics and algebra:

- [`vignette("number-theory")`](https://adamrsophie.github.io/symcompR/articles/number-theory.md):
  primes, factorization, and modular arithmetic.
- [`vignette("combinatorics")`](https://adamrsophie.github.io/symcompR/articles/combinatorics.md):
  partitions, subsets, Gray codes, and other combinatorial objects.
- [`vignette("permutation-groups")`](https://adamrsophie.github.io/symcompR/articles/permutation-groups.md):
  permutations, permutation groups, and polyhedra.
- [`vignette("fp-groups")`](https://adamrsophie.github.io/symcompR/articles/fp-groups.md):
  free groups and finitely presented groups.
- [`vignette("galois-number-fields")`](https://adamrsophie.github.io/symcompR/articles/galois-number-fields.md):
  Galois groups and algebraic number fields.
- [`vignette("lie-algebras")`](https://adamrsophie.github.io/symcompR/articles/lie-algebras.md):
  root systems and Weyl groups.
- [`vignette("categories")`](https://adamrsophie.github.io/symcompR/articles/categories.md):
  categories and commutative diagrams.
- [`vignette("cryptography")`](https://adamrsophie.github.io/symcompR/articles/cryptography.md):
  classical ciphers and key exchange.

Geometry and physics:

- [`vignette("geometry")`](https://adamrsophie.github.io/symcompR/articles/geometry.md):
  points, lines, conics, polygons, and quaternions.
- [`vignette("vector-calculus")`](https://adamrsophie.github.io/symcompR/articles/vector-calculus.md):
  vector calculus with `CoordSys3D`.
- [`vignette("tensors")`](https://adamrsophie.github.io/symcompR/articles/tensors.md):
  arrays, index notation, and tensor canonicalization.
- [`vignette("differential-geometry")`](https://adamrsophie.github.io/symcompR/articles/differential-geometry.md):
  manifolds, fields, and forms.
- [`vignette("physics-units")`](https://adamrsophie.github.io/symcompR/articles/physics-units.md):
  physical units and dimensions.

Input, output, and utilities:

- [`vignette("printing")`](https://adamrsophie.github.io/symcompR/articles/printing.md):
  string, pretty, LaTeX, and MathML output.
- [`vignette("code-printers")`](https://adamrsophie.github.io/symcompR/articles/code-printers.md):
  code generation for other languages.
- [`vignette("parsing")`](https://adamrsophie.github.io/symcompR/articles/parsing.md):
  parsing strings, LaTeX, and Mathematica input.
- [`vignette("other-modules")`](https://adamrsophie.github.io/symcompR/articles/other-modules.md):
  statistics, discrete transforms, plotting, and utilities.
- [`vignette("aliases-and-conditions")`](https://adamrsophie.github.io/symcompR/articles/aliases-and-conditions.md):
  alternative names and error conditions.
