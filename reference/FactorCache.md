# Cache of prime factors

`factor_cache` remembers one prime factor (the largest seen) for
integers that
[`factorint()`](https://adamrsophie.github.io/symcompR/reference/factorint.md)
has split, and
[`factorint()`](https://adamrsophie.github.io/symcompR/reference/factorint.md)
consults it first. `FactorCache(maxsize)` creates another cache; the
oldest entries are dropped beyond `maxsize`. Methods: `get(n)` returns a
stored factor or `NULL` (primes within the sieve return themselves),
`set(n, factor)` stores a prime factor, `add(n, factors)` stores every
prime factor of `n` given in `factors`, `remove(n)` drops an entry,
`cache_clear()` empties the cache, `set_maxsize(m)` changes the bound
and
[`size()`](https://adamrsophie.github.io/symcompR/reference/codegen-fnodes.md)
counts entries.

## Usage

``` r
FactorCache(maxsize = NULL)

factor_cache
```

## Format

`factor_cache` is the shared cache, bounded at 1000 entries.

## Arguments

- maxsize:

  Maximum number of entries, or `NULL` for no bound.

## Value

A `FactorCache` environment.

## Examples

``` r
fc <- FactorCache(10)
fc$set(35, 7)
fc$get(35)
#> [1] 7
```
