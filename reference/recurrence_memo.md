# Memoized recurrences

`recurrence_memo(initial)` returns a decorator: given `f(n, prev)`,
which computes term `n` from the list `prev` of the terms before it, it
returns a function of `n` that caches every term computed, starting from
the terms in `initial`. `assoc_recurrence_memo(base_seq)` does the same
for a two index sequence `F(n, m)` with `F(i, 0) = base_seq(i)` and
`f(i, j, prev)` giving `F(i, j)` for `1 <= j <= i`, where
`prev[[i + 1]]` holds the terms `F(i, 0), ..., F(i, j - 1)` (indices
start at 0 as in SymPy, so the list positions are shifted by one). As in
SymPy, only `m <= n` is cached. The function made by `recurrence_memo()`
carries the attributes `cache_length` and `fetch_item`.

## Usage

``` r
recurrence_memo(initial)

assoc_recurrence_memo(base_seq)
```

## Arguments

- initial:

  A list of the first terms.

- base_seq:

  A function of `i` giving `F(i, 0)`.

## Value

A decorator: a function taking `f` and returning the memoized function.

## Examples

``` r
fib <- recurrence_memo(list(0, 1))(function(n, prev) prev[[n]] + prev[[n - 1]])
fib(30)
#> [1] 832040
```
