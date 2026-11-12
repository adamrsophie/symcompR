# Names of the transitive subgroups of S1 to S6

Each of these lists holds the names of the transitive subgroups of the
symmetric group of the given degree, up to conjugacy, as in Cohen's *A
Course in Computational Algebraic Number Theory* and SymPy's
`sympy.combinatorics.galois`. A name prints like SymPy's enum member
(`S4TransitiveSubgroups.V`);
[`get_perm_group()`](https://adamrsophie.github.io/symcompR/reference/get_perm_group.md)
returns the group itself.
[`galois_group()`](https://adamrsophie.github.io/symcompR/reference/galois_group.md)
returns these names when `by_name = TRUE`.

## Usage

``` r
S1TransitiveSubgroups

S2TransitiveSubgroups

S3TransitiveSubgroups

S4TransitiveSubgroups

S5TransitiveSubgroups

S6TransitiveSubgroups
```

## Format

Named lists of `TransitiveSubgroupName` objects.

## Examples

``` r
S4TransitiveSubgroups$V
#> S4TransitiveSubgroups.V
names(S6TransitiveSubgroups)
#>  [1] "C6"     "S3"     "D6"     "A4"     "G18"    "A4xC2"  "S4m"    "S4p"   
#>  [9] "G36m"   "G36p"   "S4xC2"  "PSL2F5" "G72"    "PGL2F5" "A6"     "S6"    
```
