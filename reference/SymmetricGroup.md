# Named permutation groups

SymPy's named groups, with the same generators and the properties SymPy
records on them: `SymmetricGroup(n)`, `AlternatingGroup(n)`,
`CyclicGroup(n)`, `DihedralGroup(n)`, `AbelianGroup(...)` (the direct
product of cyclic groups of the given orders) and `RubikGroup(n)`.
`DirectProduct(...)` is the direct product of permutation groups;
`rubik(n)` and `rubik_cube_generators()` give the generators of the
Rubik's cube groups.

## Usage

``` r
SymmetricGroup(n)

AlternatingGroup(n)

CyclicGroup(n)

DihedralGroup(n)

AbelianGroup(...)

DirectProduct(...)

RubikGroup(n)

rubik_cube_generators()

rubik(n)
```

## Arguments

- n:

  Degree, or the cube dimension.

- ...:

  Cyclic orders for `AbelianGroup()`, groups for `DirectProduct()`.

## Value

A `PermutationGroup`, or a list of permutations.

## Examples

``` r
group_order(AlternatingGroup(5))
#> [1] 60
group_generators(SymmetricGroup(4))
#> [[1]]
#> (0 1 2 3)
#> 
#> [[2]]
#> (3)(0 1)
#> 
group_order(AbelianGroup(2, 3, 4))
#> [1] 24
```
