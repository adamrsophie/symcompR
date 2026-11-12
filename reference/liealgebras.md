# Lie algebras

`CartanType("A3")` (or `CartanType(list("A", 3))`) is a simple Lie
algebra of series A to G. `cartan_rank()`, `dimension()`,
`basic_root()`, `simple_root()`, `positive_roots()`, `highest_root()`
(series A),
[`roots()`](https://adamrsophie.github.io/symcompR/reference/roots.md)
(the number of roots), `cartan_matrix()`, `basis()`, `lie_algebra()` and
`dynkin_diagram()` are SymPy's methods. `CartanMatrix(t)` and
`DynkinDiagram(t)` take a type. `RootSystem(t)` has `simple_roots()`,
[`all_roots()`](https://adamrsophie.github.io/symcompR/reference/poly-aliases.md),
`root_space()`, `add_simple_roots()` and `add_as_roots()`.
`WeylGroup(t)` has `generators()`,
[`group_order()`](https://adamrsophie.github.io/symcompR/reference/PermutationGroup.md),
`group_name()`, `element_order()`, `delete_doubles()`, `matrix_form()`
and `coxeter_diagram()`; Weyl group elements are strings such as
`"r1*r2"`.

## Usage

``` r
CartanType(t)

cartan_rank(x)

dimension(ct)

basic_root(ct, i, j)

simple_root(ct, i)

positive_roots(ct)

highest_root(ct)

cartan_matrix(ct)

basis(ct)

lie_algebra(ct)

dynkin_diagram(ct)

CartanMatrix(t)

DynkinDiagram(t)

RootSystem(t)

simple_roots(x)

root_space(x)

add_simple_roots(x, root1, root2)

add_as_roots(x, root1, root2)

WeylGroup(t)

generators(x)

group_name(x)

delete_doubles(x, reflections)

matrix_form(x, weylelt)

element_order(x, weylelt)

coxeter_diagram(x)
```

## Arguments

- x:

  A Cartan type, root system or Weyl group.

- ct, t:

  A Cartan type, or a string such as `"B4"`.

- i, j:

  Indices from 1 (`simple_root`) or 0 (`basic_root`).

- root1, root2:

  Roots or simple root indices.

- reflections:

  A vector of reflection numbers.

- weylelt:

  A Weyl group element as a string.

## Value

A Cartan type, root system or Weyl group, a list of roots, a matrix, a
string or a number.

## Examples

``` r
ct <- CartanType("B3")
simple_root(ct, 3)
#> [[1]]
#> 0
#> 
#> [[2]]
#> 0
#> 
#> [[3]]
#> 1
#> 
cartan_matrix(ct)
#> Matrix(3, 3)
#> [2 , -1, 0 ]
#> [-1, 2 , -2]
#> [0 , -1, 2 ]
cat(dynkin_diagram(CartanType("D5")))
#>         4
#>         0
#>         |
#>         |
#> 0---0---0---0
#> 1   2   3   5
element_order(WeylGroup("A3"), "r1*r2*r3")
#> [1] 4
```
