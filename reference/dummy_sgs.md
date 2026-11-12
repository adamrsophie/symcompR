# Strong generators of the dummy index symmetry

Port of SymPy's `dummy_sgs()`. `dummies[2k]` and `dummies[2k+1]` are
paired (contravariant, covariant) dummy indices in consecutive
positions. The generators exchange the two indices of a pair (with a
sign change for an antisymmetric metric) and rename consecutive pairs.

## Usage

``` r
dummy_sgs(dummies, sym, n)
```

## Arguments

- dummies:

  Integer vector of dummy indices (0 based).

- sym:

  Metric symmetry: `NULL` (none), `0` (symmetric) or `1`
  (antisymmetric).

- n:

  Number of indices.

## Value

A list of array forms (integer vectors of length `n + 2`).

## Examples

``` r
dummy_sgs(2:7, 0, 8)
#> [[1]]
#>  [1] 0 1 3 2 4 5 6 7 8 9
#> 
#> [[2]]
#>  [1] 0 1 2 3 5 4 6 7 8 9
#> 
#> [[3]]
#>  [1] 0 1 2 3 4 5 7 6 8 9
#> 
#> [[4]]
#>  [1] 0 1 4 5 2 3 6 7 8 9
#> 
#> [[5]]
#>  [1] 0 1 2 3 6 7 4 5 8 9
#> 
```
