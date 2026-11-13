# Prufer sequences

`Prufer(x, n)` is SymPy's `Prufer`: a labeled tree on the nodes
`0, ..., n - 1`, given either as a list of edges (pairs) or as its
Prufer sequence. The fields `$prufer_repr`, `$tree_repr`, `$nodes`,
`$rank` and `$size` mirror SymPy's. `prufer_next()` and `prufer_prev()`
move through the trees in rank order, `prufer_unrank(rank, n)` builds
the tree of a rank, `prufer_edges(...)` turns runs of connected nodes
into a list of edges and a node count (relabeled from 0), and
`prufer_to_prufer()` and `prufer_to_tree()` convert between the two
representations.

## Usage

``` r
Prufer(x, n = NULL)

prufer_next(p, delta = 1)

prufer_prev(p, delta = 1)

prufer_unrank(rank, n)

prufer_to_prufer(tree, n)

prufer_to_tree(prufer)

prufer_edges(...)
```

## Arguments

- x:

  A list of edges or a Prufer sequence.

- n:

  Number of nodes.

- p:

  A `Prufer`.

- delta:

  Steps to move.

- rank:

  A rank.

- tree:

  A list of edges.

- prufer:

  A Prufer sequence.

- ...:

  Runs of nodes (vectors).

## Value

A `Prufer`, a sequence, or a list of edges.

## Examples

``` r
a <- Prufer(list(c(0, 1), c(0, 2), c(0, 3), c(0, 4)))
a$prufer_repr
#> [1] 0 0 0
a$rank
#> [1] 0
prufer_next(a)$tree_repr
#> [[1]]
#> [1] 0 2
#> 
#> [[2]]
#> [1] 0 3
#> 
#> [[3]]
#> [1] 0 1
#> 
#> [[4]]
#> [1] 1 4
#> 
Prufer(c(1, 1))$tree_repr
#> [[1]]
#> [1] 0 1
#> 
#> [[2]]
#> [1] 1 2
#> 
#> [[3]]
#> [1] 1 3
#> 
prufer_edges(c(0, 1, 2, 3), c(1, 4, 5), c(1, 4, 6))
#> [[1]]
#> [[1]][[1]]
#> [1] 0 1
#> 
#> [[1]][[2]]
#> [1] 1 2
#> 
#> [[1]][[3]]
#> [1] 1 4
#> 
#> [[1]][[4]]
#> [1] 2 3
#> 
#> [[1]][[5]]
#> [1] 4 5
#> 
#> [[1]][[6]]
#> [1] 4 6
#> 
#> 
#> [[2]]
#> [1] 7
#> 
```
