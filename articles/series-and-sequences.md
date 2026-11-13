# Sequences, Fourier series, and formal power series

The calculus vignette covers Taylor and Laurent series with
[`series()`](https://adamrsophie.github.io/symcompR/reference/series.md).
This vignette covers the other series objects of SymPy: sequences,
limits of sequences, Fourier series, formal power series with an
explicit formula for the coefficients, holonomic functions, and rational
approximation.

``` r

library(symcompR)
x <- Symbol("x"); z <- Symbol("z"); k <- Symbol("k")
n <- Symbol("n", integer = TRUE)
```

## Sequences

[`SeqFormula()`](https://adamrsophie.github.io/symcompR/reference/sequences.md)
defines a sequence by a formula,
[`SeqPer()`](https://adamrsophie.github.io/symcompR/reference/sequences.md)
a periodic sequence.
[`sequence_terms()`](https://adamrsophie.github.io/symcompR/reference/sequences.md)
lists the first terms.

``` r

s <- SeqFormula(n^2, list(n, 0, oo))
s
#> SeqFormula(n^2, (n, 0, oo))
sequence_terms(s, 5)
#> [[1]]
#> 0
#> 
#> [[2]]
#> 1
#> 
#> [[3]]
#> 4
#> 
#> [[4]]
#> 9
#> 
#> [[5]]
#> 16
sequence_terms(SeqPer(list(1, 2, 3)), 7)
#> [[1]]
#> 1
#> 
#> [[2]]
#> 2
#> 
#> [[3]]
#> 3
#> 
#> [[4]]
#> 1
#> 
#> [[5]]
#> 2
#> 
#> [[6]]
#> 3
#> 
#> [[7]]
#> 1
```

Sequences add and multiply term by term:

``` r

SeqFormula(n^2, list(n, 0, oo)) + SeqFormula(n^3, list(n, 0, oo))
#> SeqFormula(n^3 + n^2, (n, 0, oo))
sequence_terms(SeqFormula(n^2, list(n, 0, oo)) + SeqFormula(n, list(n, 0, oo)), 5)
#> [[1]]
#> 0
#> 
#> [[2]]
#> 2
#> 
#> [[3]]
#> 6
#> 
#> [[4]]
#> 12
#> 
#> [[5]]
#> 20
alt <- SeqFormula(n, list(n, 0, oo)) * SeqPer(list(1, -1), list(n, 0, oo))
sequence_terms(alt, 6)
#> [[1]]
#> 0
#> 
#> [[2]]
#> -1
#> 
#> [[3]]
#> 2
#> 
#> [[4]]
#> -3
#> 
#> [[5]]
#> 4
#> 
#> [[6]]
#> -5
```

[`SeqAdd()`](https://adamrsophie.github.io/symcompR/reference/sequences.md)
and
[`SeqMul()`](https://adamrsophie.github.io/symcompR/reference/sequences.md)
are the classes behind `+` and `*`, and `EmptySequence` is the neutral
sequence with no terms.

``` r

sq_seq <- SeqFormula(n^2, list(n, 0, oo)); pp <- SeqPer(list(1, 2), list(n, 0, oo))
unlist(lapply(sequence_terms(SeqAdd(sq_seq, pp), 5), sstr))
#> [1] "1"  "3"  "5"  "11" "17"
unlist(lapply(sequence_terms(SeqMul(sq_seq, pp), 5), sstr))
#> [1] "0"  "2"  "4"  "18" "16"
EmptySequence
#> EmptySequence
SeqAdd(sq_seq, EmptySequence)
#> SeqFormula(n^2, (n, 0, oo))
```

[`find_linear_recurrence()`](https://adamrsophie.github.io/symcompR/reference/sequences.md)
guesses a linear recurrence from the first terms. For the Fibonacci
numbers it finds a_n = a\_{n-1} + a\_{n-2}, and with `gfvar` it also
returns the generating function.

``` r

fib <- sym_sequence(fibonacci(n))
find_linear_recurrence(fib, 10)
#> [[1]]
#> 1
#> 
#> [[2]]
#> 1
find_linear_recurrence(fib, 10, gfvar = x)
#> [[1]]
#> [[1]][[1]]
#> 1
#> 
#> [[1]][[2]]
#> 1
#> 
#> 
#> [[2]]
#> -x/(x^2 + x - 1)
```

## Limits of sequences

[`limit_seq()`](https://adamrsophie.github.io/symcompR/reference/limit_seq.md)
computes the limit of a sequence as n \to \infty. It handles oscillating
sequences (returning the range of accumulation points) and sequences
defined by sums.

``` r

limit_seq(n / (n + 1), n)
#> 1
limit_seq((-1)^n, n)
#> AccumBounds(-1, 1)
limit_seq(harmonic(n) / log(n), n)
#> 1
limit_seq((5^n + 3^n) / (3^n + 5^n * n), n)
#> 0
limit_seq(Sum(1 / k, list(k, 1, n)) / n, n)
#> 0
```

[`dominant()`](https://adamrsophie.github.io/symcompR/reference/limit_seq.md)
finds the dominant term of a sum as n \to \infty, and
[`difference_delta()`](https://adamrsophie.github.io/symcompR/reference/series-helpers.md)
is the forward difference operator.

``` r

dominant(5 * n^3 + 4 * n^2 + n + 1, n)
#> 5*n^3
dominant(2^n + n^5, n)
#> 2^n
difference_delta(n^2, n)
#> -n^2 + (n + 1)^2
difference_delta(n * (n + 1), n, 2)
#> 4*n + 6
```

[`hypersimp()`](https://adamrsophie.github.io/symcompR/reference/series-helpers.md)
returns the ratio of consecutive terms when it is a rational function,
which shows that a term is hypergeometric.

``` r

hypersimp(sym_factorial(k), k)
#> k + 1
hypersimp(sym_binomial(n, k), k)
#> (-k + n)/(k + 1)
hypersimilar(sym_factorial(k), sym_factorial(k + 3), k)
#> [1] TRUE
```

## Fourier series

`fourier_series(f, list(x, a, b))` computes the Fourier series of `f` on
the interval \[a, b\]; `truncate(s, n)` returns its first `n` nonzero
terms.

``` r

s <- fourier_series(x^2, list(x, -sym_pi, sym_pi))
truncate(s, 4)
#> -4*cos(x) + cos(2*x) - 4*cos(3*x)/9 + pi^2/3
truncate(fourier_series(x, list(x, -sym_pi, sym_pi)), 3)
#> 2*sin(x) - sin(2*x) + 2*sin(3*x)/3
truncate(fourier_series(sin(x)^2, list(x, -sym_pi, sym_pi)), 3)
#> 1/2 - cos(2*x)/2
truncate(fourier_series(x, list(x, -1, 1)), 3)
#> 2*sin(pi*x)/pi - sin(2*pi*x)/pi + 2*sin(3*pi*x)/(3*pi)
```

The square wave:

``` r

sq <- Piecewise(list(0, Lt(x, 0)), list(1, TRUE))
truncate(fourier_series(sq, list(x, -sym_pi, sym_pi)), 4)
#> 2*sin(x)/pi + 2*sin(3*x)/(3*pi) + 2*sin(5*x)/(5*pi) + 1/2
```

Series can be shifted and scaled, and
[`sigma_approximation()`](https://adamrsophie.github.io/symcompR/reference/fourier_series.md)
applies Lanczos sigma factors that reduce the Gibbs phenomenon.

``` r

truncate(shift(s, 1), 3)
#> -4*cos(x) + cos(2*x) + 1 + pi^2/3
truncate(scale(s, 2), 3)
#> -8*cos(x) + 2*cos(2*x) + 2*pi^2/3
sigma_approximation(fourier_series(x, list(x, -sym_pi, sym_pi)), 4)
#> 2*sin(x)*sinc(pi/4) - 2*sin(2*x)/pi + 2*sin(3*x)*sinc(3*pi/4)/3
truncate(shiftx(fourier_series(x, list(x, -sym_pi, sym_pi)), 1), 3)
#> 2*sin(x + 1) - sin(2*x + 2) + 2*sin(3*x + 3)/3
truncate(scalex(fourier_series(x, list(x, -sym_pi, sym_pi)), 2), 3)
#> 2*sin(2*x) - sin(4*x) + 2*sin(6*x)/3
```

## Formal power series

[`fps()`](https://adamrsophie.github.io/symcompR/reference/fps.md)
computes a formal power series: unlike
[`series()`](https://adamrsophie.github.io/symcompR/reference/series.md),
it finds a formula for the general coefficient, so the series can be
expanded to any order.

``` r

fps_truncate(fps(exp(x), x))
#> 1 + x + x^2/2 + x^3/6 + x^4/24 + x^5/120 + O(x^6)
fps_truncate(fps(log(1 + x), x), 4)
#> x - x^2/2 + x^3/3 + O(x^4)
fps_truncate(fps(sin(x), x), 8)
#> x - x^3/6 + x^5/120 - x^7/5040 + O(x^8)
fps_truncate(fps(atan(x), x), 8)
#> x - x^3/3 + x^5/5 - x^7/7 + O(x^8)
fps_truncate(fps(1 / (1 - x), x), 5)
#> 1 + x + x^2 + x^3 + x^4 + O(x^5)
fps_term(fps(exp(x), x), 3)
#> x^3/6
fps_polynomial(fps(exp(x), x), 4)
#> x^3/6 + x^2/2 + x + 1
```

Formal power series can be integrated, multiplied, composed, and
inverted:

``` r

fps_truncate(fps_integrate(fps(exp(x), x), x), 4)
#> 1 + x + x^2/2 + x^3/6 + O(x^4)
fps_truncate(fps_product(fps(sin(x), x), fps(cos(x), x), 6), 6)
#> x - 2*x^3/3 + 2*x^5/15 + O(x^6)
fps_truncate(fps_compose(fps(exp(x), x), fps(sin(x), x), 5), 5)
#> 1 + x + x^2/2 - x^4/8 + O(x^5)
fps_truncate(fps_inverse(fps(exp(x), x), 5), 5)
#> 1 - x + x^2/2 - x^3/6 + x^4/24 + O(x^5)
```

[`compute_fps()`](https://adamrsophie.github.io/symcompR/reference/fps.md)
is the engine behind
[`fps()`](https://adamrsophie.github.io/symcompR/reference/fps.md). It
returns the coefficient sequence, the sequence of powers, and the
independent terms.

``` r

cf <- compute_fps(exp(x), x)
cf[[1]]
#> SeqFormula(Piecewise((1/factorial(_k), Eq(Mod(_k, 1), 0)), (0, True)), (_k, 1, oo))
cf[[2]]
#> SeqFormula(x^_k, (_k, 0, oo))
```

## Holonomic functions

Holonomic functions have their own vignette,
[`vignette("holonomic", package = "symcompR")`](https://adamrsophie.github.io/symcompR/articles/holonomic.md).

## Rational approximation

[`approximants()`](https://adamrsophie.github.io/symcompR/reference/approximants.md)
returns the successive Padé approximants of a sequence of coefficients.
It returns a list;
[`sstr()`](https://adamrsophie.github.io/symcompR/reference/sstr.md)
prints each element on one line. The Catalan numbers have generating
function (1 - \sqrt{1 - 4x}) / (2x), and the approximants converge to
it:

``` r

vapply(approximants(list(1, 1, 2, 5, 14, 42, 132)), sstr, "")
#> [1] "1"                                                 
#> [2] "-1/(x - 1)"                                        
#> [3] "(x - 1)/(2*x - 1)"                                 
#> [4] "(1 - 2*x)/(x^2 - 3*x + 1)"                         
#> [5] "(x^2 - 3*x + 1)/(3*x^2 - 4*x + 1)"                 
#> [6] "(-3*x^2 + 4*x - 1)/(x^3 - 6*x^2 + 5*x - 1)"        
#> [7] "(x^3 - 6*x^2 + 5*x - 1)/(4*x^3 - 10*x^2 + 6*x - 1)"
```

For the Fibonacci numbers the approximants stop changing once they reach
the exact rational generating function:

``` r

vapply(approximants(lapply(1:9, fibonacci), x), sstr, "")
#> [1] "1"                 "-1/(x - 1)"        "(x - 1)/(2*x - 1)"
#> [4] "-1/(x^2 + x - 1)"
```

[`pade()`](https://adamrsophie.github.io/symcompR/reference/pade.md)
computes the \[m/n\] Padé approximant of a function at zero from its
Taylor series. `pade_approximant(f, x, x0, m, n)` is SymPy’s version
from `sympy.series.approximants`; it expands at `x0`:

``` r

pade(cos(x), x, 2, 2)
#> (1 - 5*x^2/12)/(x^2/12 + 1)
pade_approximant(cos(x), x, 0, 2, 2)
#> (1 - 5*x^2/12)/(x^2/12 + 1)
```

## Order terms and other expansions

[`Order()`](https://adamrsophie.github.io/symcompR/reference/O.md) (also
[`O()`](https://adamrsophie.github.io/symcompR/reference/O.md)) is the
order symbol, which absorbs higher order terms.
[`nseries()`](https://adamrsophie.github.io/symcompR/reference/series.md)
is the expansion engine behind
[`series()`](https://adamrsophie.github.io/symcompR/reference/series.md),
and
[`aseries()`](https://adamrsophie.github.io/symcompR/reference/series.md)
computes asymptotic expansions at infinity.

``` r

Order(x^2)
#> O(x^2)
Order(x^2) + x^3
#> O(x^2)
nseries(sin(x), x, 0, 6)
#> x - x^3/6 + x^5/120 + O(x^6)
nseries(exp(x), x, 0, 4)
#> 1 + x + x^2/2 + x^3/6 + O(x^4)
aseries(exp(1 / x + exp(-x)) - exp(1 / x), x, 3)
#> (1/(2*x^2) + 1/x + 1 + O(x^(-3), (x, oo)))*exp(-x)
```

## Residues

``` r

residue(1 / (x^2 + 1), x, sym_I)
#> -I/2
residue(exp(z) / z^3, z, 0)
#> 1/2
residue(cot(sym_pi * z), z, 0)
#> 1/pi
```
