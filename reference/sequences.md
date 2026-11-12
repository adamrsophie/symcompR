# Sequences

`SeqFormula(formula, limits)` is the sequence of values of `formula` for
the variable running over `limits = list(n, start, stop)` (default
`list(n, 0, oo)`). `SeqPer(periodical, limits)` repeats a list of
values. `sym_sequence(seq, limits)` picks one of them, like SymPy's
[`sequence()`](https://rdrr.io/r/base/sequence.html) (renamed so it does
not mask [`base::sequence()`](https://rdrr.io/r/base/sequence.html)).
`+`, `-` and `*` combine sequences termwise into `SeqAdd` and `SeqMul`,
reducing like terms; `-s` negates. `coeff(s, pt)` is the term at `pt`
and `sequence_terms(s, n)` the first `n` terms (SymPy's `s[:n]`).
`find_linear_recurrence(s, n, d, gfvar)` finds a linear recurrence among
the first `n` terms and, with `gfvar`, the generating function.

## Usage

``` r
SeqFormula(formula, limits = NULL)

SeqPer(periodical, limits = NULL)

sym_sequence(seq, limits = NULL)

EmptySequence

sequence_terms(s, n = 10)

SeqAdd(...)

SeqMul(...)

find_linear_recurrence(s, n, d = NULL, gfvar = NULL)
```

## Arguments

- formula, seq:

  An expression, or a list of values for `SeqPer()`.

- limits:

  `list(n, start, stop)` or `list(start, stop)`.

- periodical:

  A list of values repeated by the sequence.

- s:

  A sequence.

- n:

  Number of terms.

- ...:

  Sequences.

- d:

  Largest recurrence order to try.

- gfvar:

  A symbol for the generating function.

## Value

A sequence, a list of terms, or for `find_linear_recurrence()` the
coefficients (with `gfvar`, a list of the coefficients and the
generating function).

## Examples

``` r
n <- Symbol("n")
s <- SeqFormula(n^2, list(n, 0, Inf))
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
#> 
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
#> 
SeqFormula(n^2) + SeqFormula(n^3)
#> SeqFormula(n^3 + n^2, (n, 0, oo))
find_linear_recurrence(sym_sequence(fibonacci(n)), 10)
#> [[1]]
#> 1
#> 
#> [[2]]
#> 1
#> 
```
