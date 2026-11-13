# Lie algebras

``` r

library(symcompR)
```

This vignette covers `sympy.liealgebras`: Cartan types, Cartan matrices,
Dynkin diagrams, root systems, and Weyl groups of the simple Lie
algebras.

## Cartan types and root systems

The root systems of the simple Lie algebras are available through
[`CartanType()`](https://adamrsophie.github.io/symcompR/reference/liealgebras.md);
Weyl groups act by reflections, and
[`generators()`](https://adamrsophie.github.io/symcompR/reference/liealgebras.md)
names the simple reflections. Roots are returned as lists of
coordinates, shown here with
[`sstr()`](https://adamrsophie.github.io/symcompR/reference/sstr.md).

``` r

ct <- CartanType("B3")
vapply(simple_root(ct, 3), sstr, "")
#> [1] "0" "0" "1"
cartan_matrix(ct)
#> Matrix(3, 3)
#> [2 , -1, 0 ]
#> [-1, 2 , -2]
#> [0 , -1, 2 ]
cartan_matrix(CartanType("G2"))
#> Matrix(2, 2)
#> [2 , -1]
#> [-3, 2 ]
dimension(CartanType("A3"))
#> [1] 4
basis(CartanType("A2"))
#> [1] 3
length(positive_roots(CartanType("E8")))
#> [1] 120
vapply(highest_root(CartanType("A2")), sstr, "")
#> [1] "1"  "0"  "-1"
W <- WeylGroup("A3")
group_order(W)
#> [1] 24
generators(W)
#> [1] "r1" "r2" "r3"
element_order(W, "r1*r2*r3")
#> [1] 4
matrix_form(WeylGroup("A2"), "r1*r2")
#> Matrix(3, 3)
#> [0, 0, 1]
#> [1, 0, 0]
#> [0, 1, 0]
group_name(WeylGroup("B2"))
#> [1] "The hyperoctahedral group acting on 4 elements."
```

[`CartanMatrix()`](https://adamrsophie.github.io/symcompR/reference/liealgebras.md)
and
[`DynkinDiagram()`](https://adamrsophie.github.io/symcompR/reference/liealgebras.md)
accept a type given as a string or as a list of a series letter and a
rank.
[`cartan_rank()`](https://adamrsophie.github.io/symcompR/reference/liealgebras.md)
is the rank and
[`lie_algebra()`](https://adamrsophie.github.io/symcompR/reference/liealgebras.md)
the name of the classical Lie algebra.

``` r

CartanMatrix(list("G", 2))
#> Matrix(2, 2)
#> [2 , -1]
#> [-3, 2 ]
cat(DynkinDiagram(list("B", 3)), "\n")
#> 0---0=>=0
#> 1   2   3
cat(dynkin_diagram(CartanType("F4")), "\n")
#> 0---0=>=0---0
#> 1   2   3   4
cartan_rank(CartanType("E7"))
#> [1] 7
lie_algebra(CartanType("C4"))
#> [1] "sp(8)"
lie_algebra(CartanType("A3"))
#> [1] "su(4)"
basic_root(CartanType("A3"), 1, 3)
#> [[1]]
#> 0
#> 
#> [[2]]
#> 1
#> 
#> [[3]]
#> 0
#> 
#> [[4]]
#> -1
```

[`RootSystem()`](https://adamrsophie.github.io/symcompR/reference/liealgebras.md)
gives access to the simple roots, to all roots, and to sums of roots.
[`add_simple_roots()`](https://adamrsophie.github.io/symcompR/reference/liealgebras.md)
adds two simple roots and
[`add_as_roots()`](https://adamrsophie.github.io/symcompR/reference/liealgebras.md)
adds two roots, returning the sum only when it is again a root.

``` r

rs <- RootSystem("B2")
simple_roots(rs)
#> $`1`
#> $`1`[[1]]
#> 1
#> 
#> $`1`[[2]]
#> -1
#> 
#> 
#> $`2`
#> $`2`[[1]]
#> 0
#> 
#> $`2`[[2]]
#> 1
root_space(rs)
#> [1] "alpha[1] + alpha[2]"
vapply(all_roots(rs), function(v) paste(vapply(v, sstr, ""), collapse = ", "), "")
#>        1        2        3        4        5        6        7        8 
#>  "1, -1"   "1, 1"   "1, 0"   "0, 1"  "-1, 1" "-1, -1"  "-1, 0"  "0, -1"
add_simple_roots(rs, 1, 2)
#> [[1]]
#> 1
#> 
#> [[2]]
#> 0
add_as_roots(RootSystem("A3"), c(1, 0, -1, 0), c(0, 0, 1, -1))
#> [[1]]
#> 1
#> 
#> [[2]]
#> 0
#> 
#> [[3]]
#> 0
#> 
#> [[4]]
#> -1
```

## Weyl groups

[`coxeter_diagram()`](https://adamrsophie.github.io/symcompR/reference/liealgebras.md)
draws the Coxeter diagram of a Weyl group, and
[`delete_doubles()`](https://adamrsophie.github.io/symcompR/reference/liealgebras.md)
cancels repeated adjacent reflections in a word.

``` r

cat(coxeter_diagram(WeylGroup("B2")), "\n")
#> 0===0
#> 1   2
delete_doubles(WeylGroup("G2"), c("1", "1", "2", "1", "2", "2"))
#> [1] "2" "1"
```
