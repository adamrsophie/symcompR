# Certificate of a graph

Port of `sympy.combinatorics.testutil.graph_certificate()`. Each vertex
of an undirected graph becomes a symmetric tensor with one index per
edge, contracted along the edges; the canonical form of the product
(from
[`canonicalize()`](https://adamrsophie.github.io/symcompR/reference/canonicalize.md))
is the certificate. Isomorphic graphs have equal certificates.

## Usage

``` r
graph_certificate(gr)
```

## Arguments

- gr:

  An adjacency list: a list named by the vertices `0, 1, ...` whose
  entries are the neighbors of each vertex.

## Value

An integer vector.

## Examples

``` r
gr1 <- list(`0` = c(1, 2, 3, 5), `1` = c(0, 2, 4), `2` = c(0, 1, 3, 4),
            `3` = c(0, 2, 4), `4` = c(1, 2, 3, 5), `5` = c(0, 4))
gr2 <- list(`0` = c(1, 5), `1` = c(0, 2, 3, 4), `2` = c(1, 3, 5),
            `3` = c(1, 2, 4, 5), `4` = c(1, 3, 5), `5` = c(0, 2, 3, 4))
graph_certificate(gr1)
#>  [1]  0  2  4  6  1  8 10 12  3 14 16 18  5  9 15  7 11 17 13 19 20 21
identical(graph_certificate(gr1), graph_certificate(gr2))
#> [1] TRUE
```
