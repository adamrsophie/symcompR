# Coset representatives from transversals

Port of SymPy's `transversal2coset()`: for each point `i < size`, the
sorted transversal elements when `i` is a base point, else the identity;
trailing identity entries are dropped.

## Usage

``` r
transversal2coset(size, base, transversal)
```

## Arguments

- size:

  Permutation size.

- base:

  Base points (0 based).

- transversal:

  Transversals as returned by
  [`get_transversals()`](https://adamrsophie.github.io/symcompR/reference/tensor-bsgs.md).

## Value

A list of lists of array forms.

## Examples

``` r
b <- riemann_bsgs[[1L]]
transversal2coset(6, b, get_transversals(b, riemann_bsgs[[2L]]))
#> [[1]]
#> [[1]][[1]]
#> [1] 0 1 2 3 4 5
#> 
#> [[1]][[2]]
#> [1] 1 0 2 3 5 4
#> 
#> [[1]][[3]]
#> [1] 2 3 0 1 4 5
#> 
#> [[1]][[4]]
#> [1] 3 2 0 1 5 4
#> 
#> 
#> [[2]]
#> [[2]][[1]]
#> [1] 0 1 2 3 4 5
#> 
#> 
#> [[3]]
#> [[3]][[1]]
#> [1] 0 1 2 3 4 5
#> 
#> [[3]][[2]]
#> [1] 0 1 3 2 5 4
#> 
#> 
```
