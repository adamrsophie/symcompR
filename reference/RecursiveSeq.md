# Recursively defined sequences

`RecursiveSeq(recurrence, yn, n, initial, start)` is the sequence `yn`
(an undefined function applied to the symbol `n`) defined by
`yn = recurrence` for `n >= start`, where `recurrence` uses earlier
terms `y(n - k)`, with the first terms `initial` (default dummy
symbols), as SymPy's `RecursiveSeq`. `coeff(s, i)` and
`sequence_terms(s, k)` give its terms; `recursive_seq_recurrence(s)` is
the equation `yn = recurrence`.

## Usage

``` r
RecursiveSeq(recurrence, yn, n, initial = NULL, start = 0)

recursive_seq_recurrence(s)
```

## Arguments

- recurrence:

  The expression for `yn` in terms of earlier terms.

- yn:

  The term, an undefined function applied to `n`.

- n:

  The index symbol.

- initial:

  A list of the first terms, as many as the largest shift.

- start:

  The first index.

- s:

  A `RecursiveSeq`.

## Value

A `RecursiveSeq` sequence, or an equation.

## Examples

``` r
n <- Symbol("n"); y <- Function("y")
fib <- RecursiveSeq(y(n - 1) + y(n - 2), y(n), n, list(0, 1))
sequence_terms(fib, 10)
#> [[1]]
#> 0
#> 
#> [[2]]
#> 1
#> 
#> [[3]]
#> 1
#> 
#> [[4]]
#> 2
#> 
#> [[5]]
#> 3
#> 
#> [[6]]
#> 5
#> 
#> [[7]]
#> 8
#> 
#> [[8]]
#> 13
#> 
#> [[9]]
#> 21
#> 
#> [[10]]
#> 34
#> 
```
