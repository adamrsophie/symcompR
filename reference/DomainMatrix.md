# Matrices over polynomial domains

`DomainMatrix(rows, shape, domain)` is SymPy's
`sympy.polys.matrices.DomainMatrix`: a matrix whose entries are elements
of a domain such as `ZZ`, `QQ`, `GF(p)`, `ZZ_I`, `QQ[x]`, `QQ(x)`,
`QQ<sqrt(2)>`, `RR` or `EX`. `rows` is a list of rows (each a list of
elements) for a dense matrix, or a nested named list
`list(i = list(j = element))` with 0-based indices for a sparse one.
`DM(rows, domain)` converts plain values with the domain (`c(1, 2)` is
the pair `QQ(1, 2)` or `ZZ_I(1, 2)`). The matrix methods of SymPy are
called with `$`: `A$rref()`, `A$rref_den()`, `A$det()`, `A$inv()`,
`A$charpoly()`, `A$nullspace()`, `A$lu()`, `A$fflu()`,
`A$convert_to(QQ)`, `A$to_Matrix()` and so on, with 0-based indices and
Python tuples and lists as results. `A + B`, `A - B`, `A * B` (matrix
product or scalar multiple), `A^n` and `A / x` work as in SymPy;
`A[i, j]` indexes from 1 and returns a
[`DomainScalar()`](https://adamrsophie.github.io/symcompR/reference/DomainScalar.md)
or a submatrix. `DomainMatrix_eye()`, `DomainMatrix_zeros()`,
`DomainMatrix_ones()`, `DomainMatrix_diag()`,
`DomainMatrix_from_list()`, `DomainMatrix_from_list_sympy()`,
`DomainMatrix_from_dict_sympy()`, `DomainMatrix_from_Matrix()`,
`DomainMatrix_from_rep()`, `DomainMatrix_from_list_flat()`,
`DomainMatrix_from_dod()` and `DomainMatrix_from_dok()` are the class
methods. Matrices over flint (`DFM`) are not available.

## Usage

``` r
DomainMatrix(rows, shape, domain, fmt = NULL)

DM(rows, domain)

DomainMatrix_from_rep(rep)

DomainMatrix_from_list(rows, domain)

DomainMatrix_from_list_sympy(nrows, ncols, rows, ...)

DomainMatrix_from_dict_sympy(nrows, ncols, elemsdict, ...)

DomainMatrix_from_Matrix(M, fmt = "sparse", ...)

DomainMatrix_eye(shape, domain)

DomainMatrix_diag(diagonal, domain, shape = NULL)

DomainMatrix_zeros(shape, domain, fmt = "sparse")

DomainMatrix_ones(shape, domain)

DomainMatrix_from_list_flat(elements, shape, domain)

DomainMatrix_from_dod(dod, shape, domain)

DomainMatrix_from_dok(dok, shape, domain)

# S3 method for class 'DomainMatrix'
x$name

# S3 method for class 'DomainMatrix'
format(x, ...)

# S3 method for class 'DomainMatrix'
print(x, ...)

# S3 method for class 'DomainMatrix'
x[i, j, drop = FALSE]

# S3 method for class 'DomainMatrix'
Ops(e1, e2)

# S3 method for class 'DomainMatrix'
chooseOpsMethod(x, y, mx, my, cl, reverse)
```

## Arguments

- rows:

  A list of rows, or a nested named list of nonzero entries.

- shape:

  The number of rows and columns.

- domain:

  The domain of the entries.

- fmt:

  `"dense"` or `"sparse"` to force a format.

- rep:

  A `DDM` or `SDM`.

- nrows, ncols:

  The shape.

- ...:

  Options passed to
  [`construct_domain()`](https://adamrsophie.github.io/symcompR/reference/construct_domain.md),
  or unused.

- elemsdict:

  A nested named list of entries.

- M:

  A
  [`Matrix()`](https://adamrsophie.github.io/symcompR/reference/Matrix.md).

- diagonal:

  A list of diagonal elements.

- elements:

  A flat list of elements.

- dod, dok:

  A dict of dicts or a dict of keys.

- x:

  A `DomainMatrix`.

- name:

  A method name.

- i, j:

  Row and column indices (from 1).

- drop:

  Unused.

- e1, e2:

  Matrices or scalars.

- y, mx, my, cl, reverse:

  Arguments of
  [`base::chooseOpsMethod()`](https://rdrr.io/r/base/chooseOpsMethod.html).

## Value

A `DomainMatrix`, or the value of the method called.

## Examples

``` r
A <- DM(list(list(1, 2), list(3, 4)), ZZ)
A
#> DomainMatrix([[1, 2], [3, 4]], (2, 2), ZZ)
A$det()
#> Big Integer ('bigz') :
#> [1] -2
A$rref()
#> (DomainMatrix([[1, 0], [0, 1]], (2, 2), QQ), (0, 1)) 
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
A$convert_to(QQ)$inv()
#> DomainMatrix([[-2, 1], [3/2, -1/2]], (2, 2), QQ)
A * A
#> DomainMatrix([[7, 10], [15, 22]], (2, 2), ZZ)
A[1, 2]
#> 2
```
