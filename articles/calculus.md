# Calculus: derivatives, limits, series, and sums

This vignette covers differentiation, limits, series expansions,
symbolic sums and products, finite differences, and the tools for
studying a function on an interval: singularities, extrema,
monotonicity, periodicity, and interval arithmetic with `AccumBounds`.
Integration has its own vignette.

``` r

library(symcompR)
x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); t <- Symbol("t")
h <- Symbol("h")
n <- Symbol("n", integer = TRUE); k <- Symbol("k", integer = TRUE)
f <- Function("f"); g <- Function("g")
```

## Derivatives

symcompR provides a [`diff()`](https://rdrr.io/r/base/diff.html) method
for symbolic expressions. The variables follow the expression; a number
after a variable gives the order of differentiation with respect to it.

``` r

diff(sin(x) * exp(x), x)
#> exp(x)*sin(x) + exp(x)*cos(x)
diff(atan(x), x)
#> 1/(x^2 + 1)
diff(x^4, x, 3)
#> 24*x
diff(x^2 * y^3, x, y)
#> 6*x*y^2
diff(x^2 * y^3, x, 2, y)
#> 6*y^2
```

Derivatives of undefined functions stay unevaluated, and the product and
chain rules are applied as usual.

``` r

diff(f(x) * g(x), x)
#> f(x)*Derivative(g(x), x) + g(x)*Derivative(f(x), x)
diff(f(x), x, 2)
#> Derivative(f(x), (x, 2))
```

[`Derivative()`](https://adamrsophie.github.io/symcompR/reference/Derivative.md)
builds an unevaluated derivative, which
[`doit()`](https://adamrsophie.github.io/symcompR/reference/doit.md)
evaluates. This is useful to display a derivative or to build
differential equations.

``` r

d <- Derivative(cos(x)^2, x)
d
#> Derivative(cos(x)^2, x)
doit(d)
#> -2*sin(x)*cos(x)
```

[`idiff()`](https://adamrsophie.github.io/symcompR/reference/idiff.md)
differentiates implicitly: for the circle x^2 + y^2 = 4 it returns dy/dx
with y treated as a function of x.

``` r

idiff(x^2 + y^2 - 4, y, x)
#> -x/y
```

## Limits

`limit(e, x, x0)` computes the limit of `e` as `x` approaches `x0`. It
uses the Gruntz algorithm, so it handles essential singularities and
limits at infinity. The default direction is from the right (`"+"`);
`"-"` takes the limit from the left.

``` r

limit(sin(x) / x, x, 0)
#> 1
limit((1 + 1 / x)^x, x, oo)
#> E
limit((x - sin(x)) / x^3, x, 0)
#> 1/6
limit(x^x, x, 0)
#> 1
limit(sqrt(x^2 + x) - x, x, oo)
#> 1/2
limit(log(x) / x, x, oo)
#> 0
limit(sym_factorial(n) / n^n, n, oo)
#> 0
```

One sided limits differ at a pole or a jump:

``` r

limit(1 / x, x, 0, "+")
#> oo
limit(1 / x, x, 0, "-")
#> -oo
limit(exp(-1 / x), x, 0, "-")
#> oo
```

With `"+-"` the two sided limit is required, and an error is raised when
the sides disagree.

``` r

tryCatch(limit(abs(x) / x, x, 0, "+-"), error = function(e) conditionMessage(e))
#> [1] "The limit does not exist since left hand limit = -1 and right hand limit = 1"
```

When an expression oscillates, the limit is the set of accumulation
points, returned as an `AccumBounds` interval (see below).

``` r

limit(sin(1 / x), x, 0)
#> AccumBounds(-1, 1)
```

[`Limit()`](https://adamrsophie.github.io/symcompR/reference/limit.md)
is the unevaluated form, and
[`gruntz()`](https://adamrsophie.github.io/symcompR/reference/gruntz.md)
calls the core algorithm directly.

``` r

L <- Limit(sin(x) / x, x, 0)
L
#> Limit(sin(x)/x, x, 0, dir='+')
doit(L)
#> 1
gruntz(exp(x) / x^100, x, oo)
#> oo
```

## Series expansions

`series(e, x, x0, n)` computes the expansion of `e` around `x0` up to,
but not including, order `n`. The default is `x0 = 0` and `n = 6`. The
result carries an order term
[`O()`](https://adamrsophie.github.io/symcompR/reference/O.md) that
records the truncation.

``` r

series(cos(x), x)
#> 1 - x^2/2 + x^4/24 + O(x^6)
series(exp(x), x, 0, 4)
#> 1 + x + x^2/2 + x^3/6 + O(x^4)
series(tan(x), x, 0, 8)
#> x + x^3/3 + 2*x^5/15 + 17*x^7/315 + O(x^8)
series(sin(x) / x, x, 0, 6)
#> 1 - x^2/6 + x^4/120 + O(x^6)
series(exp(x), x, 1, 3)
#> E + E*(x - 1) + E*(x - 1)^2/2 + O((x - 1)^3, (x, 1))
series(cos(x), x, sym_pi, 4)
#> -1 + (x - pi)^2/2 + O((x - pi)^4, (x, pi))
```

Laurent series, fractional powers, logarithmic terms, and expansions at
infinity are supported.

``` r

series(1 / sin(x), x, 0, 4)
#> 1/x + x/6 + 7*x^3/360 + O(x^4)
series(sqrt(x) * exp(x), x, 0, 3)
#> sqrt(x) + x^(3/2) + x^(5/2)/2 + O(x^3)
series(x^x, x, 0, 3)
#> 1 + x*log(x) + x^2*log(x)^2/2 + O(x^3*log(x)^3)
series(gamma(x), x, 0, 2)
#> 1/x - EulerGamma + x*(EulerGamma^2/2 + pi^2/12) + O(x^2)
series(exp(1 / x), x, oo, 3)
#> 1/(2*x^2) + 1/x + 1 + O(x^(-3), (x, oo))
series(atan(x), x, oo, 6)
#> -1/(5*x^5) + 1/(3*x^3) - 1/x + pi/2 + O(x^(-6), (x, oo))
```

[`removeO()`](https://adamrsophie.github.io/symcompR/reference/O.md)
drops the order term,
[`getO()`](https://adamrsophie.github.io/symcompR/reference/O.md)
extracts it, and order terms can be manipulated directly: they absorb
all higher order terms.

``` r

removeO(series(log(1 + x), x, 0, 5))
#> -x^4/4 + x^3/3 - x^2/2 + x
getO(x + x^2 + O(x^4))
#> O(x^4)
O(x^2) + x^3
#> O(x^2)
```

The leading term and a single Taylor coefficient are available directly.

``` r

as_leading_term(sin(x) / x, x)
#> 1
leadterm(1 / x + x^2, x)
#> $coeff
#> 1
#> 
#> $exponent
#> -1
taylor_term(sin(x), x, 5)
#> x^5/120
```

[`residue()`](https://adamrsophie.github.io/symcompR/reference/series-helpers.md)
computes the residue of a function at a pole.

``` r

residue(1 / (x^2 + 1), x, sym_I)
#> -I/2
```

## Sums

`summation(f, list(k, a, b))` evaluates a sum symbolically. The bounds
may be symbolic or infinite.

``` r

summation(k, list(k, 1, n))
#> n^2/2 + n/2
summation(k^2, list(k, 1, n))
#> n^3/3 + n^2/2 + n/6
factor_sym(summation(k^3, list(k, 1, n)))
#> n^2*(n + 1)^2/4
summation(2^k, list(k, 0, n - 1))
#> 2^n - 1
summation(1 / (k * (k + 1)), list(k, 1, n))
#> 1 - 1/(n + 1)
summation(1 / k^2, list(k, 1, oo))
#> pi^2/6
summation((-1)^k / (2 * k + 1), list(k, 0, oo))
#> pi/4
summation(x^k / sym_factorial(k), list(k, 0, oo))
#> exp(x)
```

[`Sum()`](https://adamrsophie.github.io/symcompR/reference/summation.md)
is the unevaluated form.
[`doit()`](https://adamrsophie.github.io/symcompR/reference/doit.md)
evaluates it, and
[`is_convergent()`](https://adamrsophie.github.io/symcompR/reference/sum-limits.md)
and
[`is_absolutely_convergent()`](https://adamrsophie.github.io/symcompR/reference/sum-limits.md)
test infinite series.

``` r

s <- Sum(1 / k^3, list(k, 1, oo))
s
#> Sum(k^(-3), (k, 1, oo))
doit(s)
#> zeta(3)
is_convergent(Sum(1 / k^2, list(k, 1, oo)))
#> [1] TRUE
is_convergent(Sum(1 / k, list(k, 1, oo)))
#> [1] FALSE
is_convergent(Sum((-1)^k / k, list(k, 1, oo)))
#> [1] TRUE
is_absolutely_convergent(Sum((-1)^k / k, list(k, 1, oo)))
#> [1] FALSE
```

[`gosper_sum()`](https://adamrsophie.github.io/symcompR/reference/summation.md)
applies Gosper’s algorithm for hypergeometric terms directly.

``` r

gosper_sum(k * sym_factorial(k), list(k, 0, n))
#> n*factorial(n) + factorial(n) - 1
```

A `Sum` can be transformed before it is evaluated.
[`change_index()`](https://adamrsophie.github.io/symcompR/reference/sum-limits.md)
substitutes the summation variable,
[`reverse_order()`](https://adamrsophie.github.io/symcompR/reference/sum-limits.md)
and
[`reorder_limit()`](https://adamrsophie.github.io/symcompR/reference/sum-limits.md)
reverse or swap limits,
[`sym_reorder()`](https://adamrsophie.github.io/symcompR/reference/sum-limits.md)
reorders the limits by variable, and
[`sum_index()`](https://adamrsophie.github.io/symcompR/reference/sum-limits.md)
gives the position of a variable among the limits.
[`euler_maclaurin()`](https://adamrsophie.github.io/symcompR/reference/sum-limits.md)
returns the Euler-Maclaurin approximation of a sum and an error
estimate.

``` r

change_index(Sum(x, list(x, 1, n)), x, x + 1, y)
#> Sum(y - 1, (y, 2, n + 1))
reverse_order(Sum(x, list(x, 1, 3)), 0)
#> Sum(-x, (x, 4, 0))
reorder_limit(Sum(x * y, list(x, 1, 3), list(y, 1, 2)), 0, 1)
#> Sum(x*y, (y, 1, 2), (x, 1, 3))
sym_reorder(Sum(x * y, list(x, 1, 3), list(y, 1, 2)), list(x, y))
#> Sum(x*y, (y, 1, 2), (x, 1, 3))
sum_index(Sum(x * y, list(x, 1, 3), list(y, 1, 2)), y)
#> [1] 1
euler_maclaurin(Sum(k, list(k, 1, 10)))[[1]]
#> 55
euler_maclaurin(Sum(1 / k^2, list(k, 1, oo)), 0, 2)[[1]]
#> 49/30
```

[`deltasummation()`](https://adamrsophie.github.io/symcompR/reference/summation.md)
and
[`deltaproduct()`](https://adamrsophie.github.io/symcompR/reference/summation.md)
handle sums and products of Kronecker deltas without evaluating the
general case.

``` r

i <- Symbol("i", integer = TRUE); j <- Symbol("j", integer = TRUE)
deltasummation(KroneckerDelta(i, k), list(k, 0, oo))
#> Piecewise((1, i >= 0), (0, True))
deltasummation(KroneckerDelta(i, k), list(k, 1, 3))
#> Piecewise((1, (i >= 1) & (i <= 3)), (0, True))
deltaproduct(KroneckerDelta(i, j), list(j, 1, 3))
#> 0
```

## Products

[`product()`](https://adamrsophie.github.io/symcompR/reference/summation.md)
and
[`Product()`](https://adamrsophie.github.io/symcompR/reference/summation.md)
work the same way for products.

``` r

product(k, list(k, 1, n))
#> factorial(n)
product(2, list(k, 1, n))
#> 2^n
product(1 + 1 / k, list(k, 1, n))
#> RisingFactorial(2, n)/factorial(n)
doit(Product(k, list(k, 1, 5)))
#> 120
```

## Finite differences

[`differentiate_finite()`](https://adamrsophie.github.io/symcompR/reference/differentiate_finite.md)
replaces derivatives by central differences.
[`finite_diff_weights()`](https://adamrsophie.github.io/symcompR/reference/calculus-utils.md)
computes the weights of a finite difference formula on any grid, and
[`apply_finite_diff()`](https://adamrsophie.github.io/symcompR/reference/calculus-utils.md)
applies them to function values.

``` r

differentiate_finite(f(x), x)
#> -f(x - 1/2) + f(x + 1/2)
w <- finite_diff_weights(1, list(-1, 0, 1), 0)
vapply(w[[2]][[3]], sstr, "")
#> [1] "-1/2" "0"    "1/2"
apply_finite_diff(1, list(x - h, x, x + h), list(f(x - h), f(x), f(x + h)), x)
#> -f(-h + x)/(2*h) + f(h + x)/(2*h)
difference_delta(n^2, n)
#> -n^2 + (n + 1)^2
```

The weights -1/2, 0, 1/2 are the familiar central difference f'(x)
\approx (f(x+h) - f(x-h)) / (2h).

## Studying a function

The calculus utilities answer questions about a function of one variable
on a domain, which defaults to the real line.

``` r

singularities(1 / (x^2 - 1), x)
#> {-1, 1}
singularities(log(x), x)
#> {0}
stationary_points(x^3 - 3 * x, x)
#> {-1, 1}
stationary_points(sin(x), x, Interval(0, 2 * sym_pi))
#> {pi/2, 3*pi/2}
minimum(x^2 + 2 * x + 3, x)
#> 2
minimum(x^3 - 3 * x, x, Interval(-2, 3))
#> -2
maximum(-x^2 + 4 * x, x)
#> 4
maximum(sin(x), x, Interval(0, sym_pi))
#> 1
```

Monotonicity and convexity:

``` r

is_increasing(x^3, Interval(-oo, oo))
#> [1] TRUE
is_strictly_increasing(exp(x))
#> [1] TRUE
is_decreasing(1 / x, Interval(0, oo, TRUE))
#> [1] TRUE
is_convex(x^2, x)
#> [1] TRUE
is_strictly_decreasing(1 / (x^2 - 3 * x), Interval_open(Rational(3, 2), 3))
#> [1] TRUE
is_monotonic(exp(x))
#> [1] TRUE
```

[`periodicity()`](https://adamrsophie.github.io/symcompR/reference/calculus-utils.md)
returns the fundamental period, or `NULL` for a function that is not
periodic, and
[`continuous_domain()`](https://adamrsophie.github.io/symcompR/reference/sum-limits.md)
finds where a function is continuous.

``` r

periodicity(sin(2 * x), x)
#> pi
periodicity(sin(x) + cos(3 * x / 2), x)
#> 4*pi
periodicity(exp(x), x)
#> NULL
continuous_domain(1 / x, x, Reals)
#> Union(Interval.open(-oo, 0), Interval.open(0, oo))
continuous_domain(log(x) + sqrt(1 - x), x, Reals)
#> Interval.Lopen(0, 1)
```

[`function_range()`](https://adamrsophie.github.io/symcompR/reference/function_range.md)
gives the range of a function on a domain,
[`monotonicity_helper()`](https://adamrsophie.github.io/symcompR/reference/monotonicity_helper.md)
tests the sign of the derivative on an interval, and
[`lcim()`](https://adamrsophie.github.io/symcompR/reference/lcim.md) is
the least common multiple of numbers whose ratios are rational (`NULL`
when there is none):

``` r

function_range(sin(x), x, Interval(0, 2 * sym_pi))
#> Interval(-1, 1)
function_range(1 / x, x, Reals)
#> Union(Interval.open(-oo, 0), Interval.open(0, oo))
function_range(exp(x), x, Interval(-1, 1))
#> Interval(exp(-1), E)
monotonicity_helper(x^3, function(d) Ge(d, 0))
#> [1] TRUE
monotonicity_helper(x^2, function(d) Ge(d, 0))
#> [1] FALSE
monotonicity_helper(x^2, function(d) Ge(d, 0), Interval(0, oo))
#> [1] TRUE
lcim(list(Rational(1, 2), 2, 3))
#> 6
lcim(list(sym_pi, sym_pi / 2, sym_pi / 4))
#> pi
lcim(list(1, sym_pi))
#> NULL
```

[`euler_equations()`](https://adamrsophie.github.io/symcompR/reference/calculus-utils.md)
derives the Euler-Lagrange equations of a functional. For the harmonic
oscillator Lagrangian L = \dot f^2/2 - f^2/2:

``` r

euler_equations(diff(f(t), t)^2 / 2 - f(t)^2 / 2, f(t), t)
#> [[1]]
#> Eq(-f(t) - Derivative(f(t), (t, 2)), 0)
```

## Interval arithmetic with `AccumBounds`

`AccumBounds(a, b)` represents an unknown value somewhere in \[a, b\].
Arithmetic on these bounds gives guaranteed enclosures, and functions
such as [`sin()`](https://rdrr.io/r/base/Trig.html) return them at
infinity.

``` r

AccumBounds(-1, 1) + AccumBounds(2, 3)
#> AccumBounds(1, 4)
AccumBounds(-1, 2)^2
#> AccumBounds(0, 4)
1 / AccumBounds(1, 2)
#> AccumBounds(1/2, 1)
ab_union(AccumBounds(1, 3), AccumBounds(2, 4))
#> AccumBounds(1, 4)
AccumulationBounds(-1, 1)
#> AccumBounds(-1, 1)
AccumBounds(0, 1) * AccumBounds(-1, 3)
#> AccumBounds(-1, 3)
sin(oo)
#> AccumBounds(-1, 1)
limit(sin(x), x, oo)
#> AccumBounds(-1, 1)
```

[`not_empty_in()`](https://adamrsophie.github.io/symcompR/reference/not_empty_in.md)
finds the values of a parameter for which a set expression is nonempty:
here, where at least one of x and 2x lies in \[1, 2\].

``` r

not_empty_in(Intersection(FiniteSet(x, 2 * x), Interval(1, 2)), x)
#> Interval(1/2, 2)
```
