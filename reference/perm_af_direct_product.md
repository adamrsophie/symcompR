# Direct products of permutation generators

Port of SymPy's `perm_af_direct_product()`. With `signed = TRUE` the
last two points of each array form carry the sign.

## Usage

``` r
perm_af_direct_product(gens1, gens2, signed = TRUE)
```

## Arguments

- gens1, gens2:

  Lists of array forms.

- signed:

  Whether the permutations are signed.

## Value

A list of array forms.

## Examples

``` r
perm_af_direct_product(list(c(1, 0, 2, 3), c(0, 1, 3, 2)), list(c(1, 0)), FALSE)
#> [[1]]
#> [1] 1 0 2 3 4 5
#> 
#> [[2]]
#> [1] 0 1 3 2 4 5
#> 
#> [[3]]
#> [1] 0 1 2 3 5 4
#> 
```
