# Dense and sparse domain matrices

`DDM(rows, shape, domain)` is a dense matrix over a polynomial domain
(SymPy's `sympy.polys.matrices.ddm.DDM`),
`SDM(elemsdict, shape, domain)` a sparse one
(`sympy.polys.matrices.sdm.SDM`) given as a nested list
`list(i = list(j = element))` with 0-based names, or as a list of
`list(i, list(list(j, element), ...))` pairs. They are the internal
representations of
[`DomainMatrix()`](https://adamrsophie.github.io/symcompR/reference/DomainMatrix.md)
(`A$rep`), with the same methods as in SymPy, called with `$`:
`M$rref()`, `M$det()`, `M$matmul(N)`, ... Indices passed to methods are
0-based as in SymPy.

## Usage

``` r
DDM(rows, shape, domain)

# S3 method for class 'DDM'
x$name

# S3 method for class 'DDM'
format(x, ...)

# S3 method for class 'DDM'
print(x, ...)

# S3 method for class 'DDM'
Ops(e1, e2)

SDM(elemsdict, shape, domain)

# S3 method for class 'SDM'
x$name

# S3 method for class 'SDM'
format(x, ...)

# S3 method for class 'SDM'
print(x, ...)

# S3 method for class 'SDM'
Ops(e1, e2)
```

## Arguments

- rows:

  A list of rows, each a list of domain elements.

- shape:

  The number of rows and columns.

- domain:

  The domain of the elements.

- x:

  A `DDM` or `SDM`.

- name:

  A field or method name.

- ...:

  Unused.

- e1, e2:

  Matrices or elements.

- elemsdict:

  The nonzero elements by row and column.

## Value

A `DDM` or `SDM` object.

## Examples

``` r
A <- DDM(list(list(ZZ(1), ZZ(2)), list(ZZ(3), ZZ(4))), c(2, 2), ZZ)
A$det()
#> Big Integer ('bigz') :
#> [1] -2
A$charpoly()
#> [[1]]
#> Big Integer ('bigz') :
#> [1] 1
#> 
#> [[2]]
#> Big Integer ('bigz') :
#> [1] -5
#> 
#> [[3]]
#> Big Integer ('bigz') :
#> [1] -2
#> 
B <- SDM(list(`0` = list(`1` = ZZ(2))), c(2, 2), ZZ)
B$to_ddm()
#> DDM([[0, 2], [0, 0]], (2, 2), ZZ)
```
