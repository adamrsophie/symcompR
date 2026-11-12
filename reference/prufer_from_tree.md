# Prufer sequences

`prufer_from_tree(edges, n)` encodes a labeled tree (0 based edges) and
`tree_from_prufer(seq)` decodes it.

## Usage

``` r
prufer_from_tree(edges, n)

tree_from_prufer(seq)
```

## Arguments

- edges:

  A list of length 2 integer vectors.

- n:

  Number of nodes.

- seq:

  A Prufer sequence.

## Value

An integer vector or a list of edges.
