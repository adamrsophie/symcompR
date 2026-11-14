# Characteristic polynomial, eigenvalues and eigenvectors

`charpoly(M, lambda)` returns `sym_det(lambda*I - M)`. `eigenvals()`
returns the roots of the characteristic polynomial with multiplicities;
`eigenvects()` returns a list of `list(value, multiplicity, vectors)`.
`diagonalize()` returns `list(P, D)` with `M = P D P^-1`.
`is_diagonalizable()` tests whether that is possible.

## Usage

``` r
charpoly(M, lambda = Symbol("lambda"))

eigenvals(
  M,
  error_when_incomplete = TRUE,
  simplify = FALSE,
  multiple = FALSE,
  rational = FALSE,
  ...
)

eigenvects(
  M,
  error_when_incomplete = TRUE,
  iszerofunc = NULL,
  simplify = FALSE,
  chop = FALSE,
  ...
)

is_diagonalizable(M, reals_only = FALSE)

diagonalize(M, reals_only = FALSE, sort = FALSE, normalize = FALSE)
```

## Arguments

- M:

  A square `symmatrix`.

- lambda:

  Symbol for the characteristic polynomial.

- error_when_incomplete:

  Signal a `MatrixError` when some eigenvalues cannot be found.

- simplify:

  `FALSE`, `TRUE` or a function applied to each result.

- multiple:

  Return the eigenvalues as a list repeated by multiplicity instead of a
  list with attribute `"multiplicity"`.

- rational:

  Replace floating point entries by rationals first.

- ...:

  Flags passed to
  [`roots()`](https://adamrsophie.github.io/symcompR/reference/roots.md).

- iszerofunc:

  A zero test returning `TRUE`, `FALSE` or `NA`; `NULL` uses the default
  test.

- chop:

  Accepted for SymPy compatibility; currently not used.

- reals_only:

  Treat a matrix with complex eigenvalues as not diagonalizable.

- sort:

  Sort the eigenvalues on the diagonal.

- normalize:

  Scale each column of `P` to unit norm.

## Value

An expression or list.

## Examples

``` r
M <- Matrix(list(list(2, 1), list(1, 2)))
charpoly(M)
#> lambda^2 - 4*lambda + 3
eigenvals(M)
#> [[1]]
#> 3
#> 
#> [[2]]
#> 1
#> 
#> attr(,"multiplicity")
#> [1] 1 1
```
