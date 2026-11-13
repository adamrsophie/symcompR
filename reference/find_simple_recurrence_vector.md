# Guessing recurrences and generating functions

`find_simple_recurrence_vector(l)` returns the coefficients of a linear
recurrence with constant coefficients satisfied by the rational sequence
`l`, or `0` when none is found; `find_simple_recurrence(v, A, N)`
returns the recurrence as an expression in `A(N), A(N + 1), ...`.
`guess_generating_function_rational(v, X)` guesses a rational ordinary
generating function, and
`guess_generating_function(v, X, types, maxsqrtn)` tries the generating
function types `"ogf"`, `"egf"`, `"lgf"`, `"hlgf"`, `"lgdogf"` and
`"lgdegf"` and returns a named list.
`guess(l, all, evaluate, niter, variables)` guesses a closed form built
from rational functions and nested products, as SymPy's `guess()`.
`rationalize(x, maxcoeff)` recovers a rational number from a float by
its continued fraction, stopping at the first partial quotient of at
least `maxcoeff`; an expression that is not a number is passed to
[`radsimp()`](https://adamrsophie.github.io/symcompR/reference/radsimp.md).

## Usage

``` r
find_simple_recurrence_vector(l)

find_simple_recurrence(v, A = Function("a"), N = Symbol("n"))

rationalize(x, maxcoeff = 10000)

guess_generating_function_rational(v, X = Symbol("x"))

guess_generating_function(v, X = Symbol("x"), types = "all", maxsqrtn = 2)

guess(l, all = FALSE, evaluate = TRUE, niter = 2, variables = NULL)
```

## Arguments

- l, v:

  A vector or list of rational numbers.

- A:

  A function made by
  [`Function()`](https://adamrsophie.github.io/symcompR/reference/Function.md).

- N, X:

  Symbols.

- x:

  A number.

- maxcoeff:

  Bound on the partial quotients.

- types:

  Generating function types, or `"all"`.

- maxsqrtn:

  The largest root of a rational function tried.

- all:

  Return every formula found instead of the first ones.

- evaluate:

  Evaluate the products.

- niter:

  Number of iterations.

- variables:

  Symbols to use instead of `i0, i1, ...`.

## Value

A list, an expression, a named list or `NULL`.

## Examples

``` r
find_simple_recurrence_vector(c(0, 1, 1, 2, 3, 5, 8, 13, 21, 34))
#> [[1]]
#> 1
#> 
#> [[2]]
#> -1
#> 
#> [[3]]
#> -1
#> 
find_simple_recurrence(c(0, 1, 1, 2, 3, 5, 8, 13, 21, 34))
#> -a(n) - a(n + 1) + a(n + 2)
guess_generating_function_rational(c(5, 8, 13, 21, 34, 55, 89, 144))
#> (3*x + 5)/(-x^2 - x + 1)
guess(c(1, 2, 6, 24, 120))
#> [[1]]
#> RisingFactorial(2, i0 - 1)
#> 
rationalize(pi, maxcoeff = 250)
#> 355/113
```
