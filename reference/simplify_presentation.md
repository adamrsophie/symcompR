# Simplifying presentations

`simplify_presentation(G)` returns a simplified isomorphic copy of the
`FpGroup` `G`; `simplify_presentation(gens, rels)` simplifies the lists
of generators and relators and returns both. With `change_gens = TRUE`
redundant generators are removed. `elimination_technique_1()` and
`elimination_technique_2()` are the Tietze transformations SymPy uses.

## Usage

``` r
simplify_presentation(..., change_gens = FALSE)

elimination_technique_1(gens, rels, identity)

elimination_technique_2(C)
```

## Arguments

- ...:

  An `FpGroup`, or the generators and the relators.

- change_gens:

  Remove redundant generators.

- gens, rels:

  Lists of words.

- identity:

  The identity of the free group.

- C:

  A coset table with Reidemeister data.

## Value

An `FpGroup` or a list of generators and relators.

## Examples

``` r
g <- free_group("x, y")
x <- g[[2]]; y <- g[[3]]
fp_str(simplify_presentation(FpGroup(g[[1]], list(x^2, x^5, y^3)))$relators)
#> [1] "[x, y^3]"
```
