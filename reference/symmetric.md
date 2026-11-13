# Permutation group elements

`symmetric(n)`, `cyclic(n)`, `alternating(n)` and `dihedral(n)` list the
elements of the symmetric, cyclic, alternating and dihedral groups of
degree `n` as permutations, in SymPy's order (SymPy's generators of
`sympy.combinatorics.generators`).

## Usage

``` r
symmetric(n)

cyclic(n)

alternating(n)

dihedral(n)
```

## Arguments

- n:

  The degree.

## Value

A list of `Permutation` objects.

## Examples

``` r
cyclic(4)
#> [[1]]
#> (3)
#> 
#> [[2]]
#> (0 1 2 3)
#> 
#> [[3]]
#> (0 2)(1 3)
#> 
#> [[4]]
#> (0 3 2 1)
#> 
length(alternating(4))
#> [1] 12
dihedral(3)
#> [[1]]
#> (2)
#> 
#> [[2]]
#> (0 2)
#> 
#> [[3]]
#> (0 1 2)
#> 
#> [[4]]
#> (1 2)
#> 
#> [[5]]
#> (0 2 1)
#> 
#> [[6]]
#> (2)(0 1)
#> 
```
