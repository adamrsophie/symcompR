# Assumptions

``` r

library(symcompR)
x <- Symbol("x"); y <- Symbol("y")
n <- Symbol("n", integer = TRUE)
```

symcompR has both assumptions systems of SymPy: facts attached to
symbols (the “old” assumptions of `sympy.core.assumptions`) and
predicates queried with
[`ask()`](https://adamrsophie.github.io/symcompR/reference/ask.md)
(`sympy.assumptions`).

## Assumptions on symbols

[`Symbol()`](https://adamrsophie.github.io/symcompR/reference/Symbol.md)
accepts assumptions as arguments. The `is_*()` queries return `TRUE`,
`FALSE`, or `NA` when the answer is not known, and
[`ask_fact()`](https://adamrsophie.github.io/symcompR/reference/ask_fact.md)
queries any fact by name.

``` r

p <- Symbol("p", positive = TRUE)
r <- Symbol("r", real = TRUE)
is_negative(1 - sqrt(S(2)))
#> [1] TRUE
is_nonnegative(p)
#> [1] TRUE
is_nonpositive(-p)
#> [1] TRUE
is_nonzero(p)
#> [1] TRUE
is_zero(p - p)
#> [1] TRUE
is_finite(r)
#> [1] TRUE
is_infinite(oo * sym_I)
#> [1] TRUE
is_complex(r)
#> [1] TRUE
is_imaginary(sym_I * r)
#> [1] NA
is_imaginary(2 * sym_I)
#> [1] TRUE
is_rational(Rational(1, 3))
#> [1] TRUE
is_irrational(sqrt(S(2)))
#> [1] TRUE
is_algebraic(sqrt(S(2)))
#> [1] TRUE
is_transcendental(sym_pi)
#> [1] TRUE
is_odd(2 * n + 1)
#> [1] TRUE
is_prime_expr(S(7))
#> [1] TRUE
is_commutative(x)
#> [1] TRUE
ask_fact(p + 1, "positive")
#> [1] TRUE
ask_fact(x, "positive")
#> [1] NA
```

[`check_assumptions()`](https://adamrsophie.github.io/symcompR/reference/core-helpers.md)
tests an expression against a set of assumptions,
[`failing_assumptions()`](https://adamrsophie.github.io/symcompR/reference/core-helpers.md)
lists the ones that fail, and
[`assumptions()`](https://adamrsophie.github.io/symcompR/reference/assumptions.md)
lists the known assumptions of an expression.

``` r

check_assumptions(S(-2), integer = TRUE, positive = TRUE)
#> [1] FALSE
check_assumptions(p + 1, positive = TRUE)
#> [1] TRUE
failing_assumptions(S(-2), integer = TRUE, positive = TRUE)
#> $positive
#> [1] FALSE
assumptions(x)
#> commutative 
#>        TRUE
```

## Predicates with `ask()` and `refine()`

The second assumptions system keeps facts separate from symbols:
predicates such as `Q$positive(x)` are queried with
[`ask()`](https://adamrsophie.github.io/symcompR/reference/ask.md),
which returns `TRUE`, `FALSE`, or `NA` (unknown). Further arguments are
assumptions to use.

``` r

ask(Q$positive(x))
#> [1] NA
ask(Q$positive(x^2 + 1), Q$real(x))
#> [1] TRUE
ask(Q$integer(x + y), Q$integer(x), Q$integer(y))
#> [1] TRUE
ask(Q$positive(x * y), Q$negative(x), Q$negative(y))
#> [1] TRUE
ask(Q$even(2 * n))
#> [1] TRUE
ask(Q$prime(7))
#> [1] TRUE
ask(Q$rational(sqrt(S(2))))
#> [1] FALSE
```

[`refine()`](https://adamrsophie.github.io/symcompR/reference/refine.md)
simplifies an expression using assumptions.

``` r

refine(sqrt(x^2), Q$positive(x))
#> x
refine(Abs(x), Q$negative(x))
#> -x
refine((-1)^n, Q$even(n))
#> 1
```

[`assuming()`](https://adamrsophie.github.io/symcompR/reference/assumptions-context.md)
evaluates an expression with assumptions added to the global context,
and removes them afterwards.

``` r

assuming(Q$positive(x), refine(sqrt(x^2)))
#> x
assuming(Q$positive(x), ask(Q$real(x)))
#> [1] TRUE
ask(Q$real(x))
#> [1] NA
```

An `AssumptionsContext` is a set of facts. `global_assumptions` is the
context that
[`ask()`](https://adamrsophie.github.io/symcompR/reference/ask.md)
always consults;
[`assuming()`](https://adamrsophie.github.io/symcompR/reference/assumptions-context.md)
adds to it temporarily.

``` r

ctx <- AssumptionsContext(Q$positive(x))
ctx
#> AssumptionsContext({Q.positive(x)})
global_assumptions
#> AssumptionsContext()
```
