# Smoothness of an integer

`smoothness(n)` returns `c(B, Bpow)`: the largest prime factor of `n`
and the largest prime power dividing it.
`smoothness_p(n, m, power, visual)` reports the smoothness of `p + m`
for every prime power `p^i` in `n`, sorted by `B` (`power = 0`) or
`Bpow` (`power = 1`). It returns `list(m, entries)` where each entry is
`c(p, i, B, Bpow)`; with `visual = TRUE` it returns the report as a
string, which it also accepts back as input.

## Usage

``` r
smoothness(n)

smoothness_p(n, m = -1, power = 0, visual = NULL)
```

## Arguments

- n:

  An integer, a factor vector, a report string or a returned list.

- m:

  Offset added to each prime.

- power:

  0 or 1, the sort key.

- visual:

  `TRUE` for the string report, `FALSE` for the list.

## Value

As described above.

## Examples

``` r
smoothness(2^7 * 3^2)
#> [1]   3 128
smoothness_p(10431, m = 1)
#> [[1]]
#> [1] 1
#> 
#> [[2]]
#> [[2]][[1]]
#> [1] 3 2 2 4
#> 
#> [[2]][[2]]
#> [1] 19  1  5  5
#> 
#> [[2]][[3]]
#> [1] 61  1 31 31
#> 
#> 
cat(smoothness_p(10431, m = 1, visual = TRUE))
#> p**i=3**2 has p+1 B=2, B-pow=4
#> p**i=19**1 has p+1 B=5, B-pow=5
#> p**i=61**1 has p+1 B=31, B-pow=31
```
