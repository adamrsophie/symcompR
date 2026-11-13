# Vectors and dyadics of classical mechanics

`PhysVector()` builds a vector of `sympy.physics.vector` from a list of
`list(matrix, frame)` pairs, and `PhysVector(0)` is the zero vector;
vectors are usually made from the unit vectors of a
[`ReferenceFrame()`](https://adamrsophie.github.io/symcompR/reference/ReferenceFrame.md),
as `x * N$x + y * N$y`. `PhysDyadic()` builds a dyadic from a list of
`list(coefficient, vector1, vector2)` triples; `N$x | N$y` is the outer
product. The names differ from SymPy's `Vector` and `Dyadic` because
`sympy.vector` has classes of the same names.

## Usage

``` r
PhysVector(inlist)

PhysDyadic(inlist)
```

## Arguments

- inlist:

  `0`, or a list of pairs (vectors) or triples (dyadics).

## Value

A `PhysVector` or `PhysDyadic`.

## Details

Operators follow SymPy: `+`, `-`, scalar `*` and `/`, `==`, `&` (dot
product), `^` (cross product) and `|` (outer product). Methods such as
`dot`, `cross`, `outer`, `express`, `dt`, `diff`, `magnitude`,
`normalize`, `subs`, `simplify`, `to_matrix`, `doit`, `applyfunc`,
`xreplace`, `separate` and `free_symbols` are reached with `$`.

## Examples

``` r
N <- ReferenceFrame("N")
x <- Symbol("x")
v <- x * N$x + 2 * N$y
v
#> x*N.x + 2*N.y
v$magnitude()
#> sqrt(x^2 + 4)
(N$x | N$y) & N$y
#> N.x
PhysVector(0)
#> 0
```
