# Matrices and linear algebra

symcompR has dense matrices of symbolic expressions with exact linear
algebra (determinants, inverses, row reduction, eigenvalues, and
canonical forms and decompositions), and symbolic matrix expressions in
which a whole matrix is a single symbol.

``` r

library(symcompR)
x <- Symbol("x"); y <- Symbol("y")
a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c")
n <- Symbol("n", integer = TRUE)
```

## Building matrices

[`Matrix()`](https://adamrsophie.github.io/symcompR/reference/Matrix.md)
takes a list of rows, an R matrix, a vector (which becomes a column), or
dimensions and a function of the row and column index. Following R,
indices start at 1 (SymPy starts at 0).

``` r

M <- Matrix(list(list(1, 2), list(3, 4)))
M
#> Matrix(2, 2)
#> [1, 2]
#> [3, 4]
Matrix(matrix(1:6, nrow = 2))
#> Matrix(2, 3)
#> [1, 3, 5]
#> [2, 4, 6]
Matrix(c(1, 2, 3))
#> Matrix(3, 1)
#> [1]
#> [2]
#> [3]
Matrix(2, 3, function(i, j) 10 * i + j)
#> Matrix(2, 3)
#> [11, 12, 13]
#> [21, 22, 23]
```

Special matrices:

``` r

eye(3)
#> Matrix(3, 3)
#> [1, 0, 0]
#> [0, 1, 0]
#> [0, 0, 1]
zeros(2, 3)
#> Matrix(2, 3)
#> [0, 0, 0]
#> [0, 0, 0]
ones(2, 2)
#> Matrix(2, 2)
#> [1, 1]
#> [1, 1]
sym_diag(1, 2, 3)
#> Matrix(3, 3)
#> [1, 0, 0]
#> [0, 2, 0]
#> [0, 0, 3]
jordan_cell(a, 3)
#> Matrix(3, 3)
#> [a, 1, 0]
#> [0, a, 1]
#> [0, 0, a]
rot_axis3(sym_pi / 2)
#> Matrix(3, 3)
#> [0 , 1, 0]
#> [-1, 0, 0]
#> [0 , 0, 1]
```

[`sym_diag()`](https://adamrsophie.github.io/symcompR/reference/matrix-builders.md)
is SymPy’s [`diag()`](https://rdrr.io/r/base/diag.html); it carries a
prefix so that it does not mask
[`base::diag()`](https://rdrr.io/r/base/diag.html).

## Indexing and shape

Matrices are indexed like R matrices, and assignment works the same way.

``` r

M[1, 2]
#> 2
M[2, ]
#> Matrix(1, 2)
#> [3, 4]
M[, 1]
#> Matrix(2, 1)
#> [1]
#> [3]
dim(M)
#> [1] 2 2
M2 <- M
M2[1, 1] <- x
M2
#> Matrix(2, 2)
#> [x, 2]
#> [3, 4]
```

[`as.matrix()`](https://rdrr.io/r/base/matrix.html) converts a numeric
symbolic matrix into an ordinary R matrix.

``` r

as.matrix(M)
#>      [,1] [,2]
#> [1,]    1    2
#> [2,]    3    4
```

## Arithmetic

`+`, `-`, and multiplication by a scalar work as usual. Because `*`
between two R matrices is elementwise, the same holds here; the matrix
product is `%*%` and the matrix power `%^%`.

``` r

M + M
#> Matrix(2, 2)
#> [2, 4]
#> [6, 8]
2 * M
#> Matrix(2, 2)
#> [2, 4]
#> [6, 8]
M * M
#> Matrix(2, 2)
#> [1, 4 ]
#> [9, 16]
M %*% M
#> Matrix(2, 2)
#> [7 , 10]
#> [15, 22]
M %^% 2
#> Matrix(2, 2)
#> [7 , 10]
#> [15, 22]
M %^% -1
#> Matrix(2, 2)
#> [-2 , 1   ]
#> [3/2, -1/2]
t(M)
#> Matrix(2, 2)
#> [1, 3]
#> [2, 4]
```

Functions such as
[`subs()`](https://adamrsophie.github.io/symcompR/reference/subs.md),
[`diff()`](https://rdrr.io/r/base/diff.html),
[`N()`](https://adamrsophie.github.io/symcompR/reference/N.md), and
[`simplify()`](https://adamrsophie.github.io/symcompR/reference/simplify.md)
work entry by entry, and
[`applyfunc()`](https://adamrsophie.github.io/symcompR/reference/dot.md)
applies any function to every entry.

``` r

P <- Matrix(list(list(1, x), list(x^2, sin(x))))
subs(P, x, 2)
#> Matrix(2, 2)
#> [1, 2     ]
#> [4, sin(2)]
diff(P, x)
#> Matrix(2, 2)
#> [0  , 1     ]
#> [2*x, cos(x)]
applyfunc(Matrix(list(list(sin(x)^2 + cos(x)^2, 1))), simplify)
#> Matrix(1, 2)
#> [1, 1]
N(Matrix(list(list(sqrt(S(2)), sym_pi))), 5)
#> Matrix(1, 2)
#> [1.4142, 3.1416]
```

## Determinant, inverse, and row reduction

``` r

sym_det(M)
#> -2
inverse(M)
#> Matrix(2, 2)
#> [-2 , 1   ]
#> [3/2, -1/2]
S <- Matrix(list(list(a, b), list(c, 1)))
sym_det(S)
#> a - b*c
inverse(S)
#> Matrix(2, 2)
#> [1/(a - b*c) , -b/(a - b*c)]
#> [-c/(a - b*c), a/(a - b*c) ]
adjugate(M)
#> Matrix(2, 2)
#> [4 , -2]
#> [-3, 1 ]
cofactor(M, 1, 2)
#> -3
permanent(M)
#> 10
```

[`sym_det()`](https://adamrsophie.github.io/symcompR/reference/sym_det.md)
uses Bareiss’ fraction free elimination by default;
`method = "berkowitz"` and `method = "laplace"` are also available. A
linear system Mx = v is solved with
[`solve()`](https://rdrr.io/r/base/solve.html):

``` r

solve(M, Matrix(c(5, 6)))
#> Matrix(2, 1)
#> [-4 ]
#> [9/2]
```

[`rref()`](https://adamrsophie.github.io/symcompR/reference/sym_det.md)
returns the reduced row echelon form and the pivot columns.
[`rank_matrix()`](https://adamrsophie.github.io/symcompR/reference/sym_det.md),
[`nullspace()`](https://adamrsophie.github.io/symcompR/reference/sym_det.md),
[`columnspace()`](https://adamrsophie.github.io/symcompR/reference/sym_det.md),
and
[`rowspace()`](https://adamrsophie.github.io/symcompR/reference/sym_det.md)
describe the fundamental subspaces.

``` r

S3 <- Matrix(list(list(1, 2, 3), list(4, 5, 6), list(7, 8, 9)))
rref(S3)
#> $matrix
#> Matrix(3, 3)
#> [1, 0, -1]
#> [0, 1, 2 ]
#> [0, 0, 0 ]
#> 
#> $pivots
#> [1] 1 2
rank_matrix(S3)
#> [1] 2
nullspace(S3)
#> [[1]]
#> Matrix(3, 1)
#> [1 ]
#> [-2]
#> [1 ]
columnspace(S3)
#> [[1]]
#> Matrix(3, 1)
#> [1]
#> [4]
#> [7]
#> 
#> [[2]]
#> Matrix(3, 1)
#> [2]
#> [5]
#> [8]
sym_trace(S3)
#> 15
```

## Eigenvalues and diagonalization

[`charpoly()`](https://adamrsophie.github.io/symcompR/reference/charpoly.md)
gives the characteristic polynomial and
[`eigenvals()`](https://adamrsophie.github.io/symcompR/reference/charpoly.md)
its roots. Exact eigenvalues can involve radicals.

``` r

A <- Matrix(list(list(2, 1), list(1, 2)))
charpoly(A, x)
#> x^2 - 4*x + 3
eigenvals(A)
#> [[1]]
#> 3
#> 
#> [[2]]
#> 1
#> 
#> attr(,"multiplicity")
#> [1] 1 1
eigenvals(M)
#> [[1]]
#> 5/2 - sqrt(33)/2
#> 
#> [[2]]
#> 5/2 + sqrt(33)/2
#> 
#> attr(,"multiplicity")
#> [1] 1 1
```

[`eigenvects()`](https://adamrsophie.github.io/symcompR/reference/charpoly.md)
returns, for each eigenvalue, its algebraic multiplicity and a basis of
the eigenspace.

``` r

ev <- eigenvects(A)
ev[[1]]$value
#> 1
ev[[1]]$vectors[[1]]
#> Matrix(2, 1)
#> [-1]
#> [1 ]
is_diagonalizable(A)
#> [1] TRUE
```

[`diagonalize()`](https://adamrsophie.github.io/symcompR/reference/charpoly.md)
returns `list(P, D)` with A = P D P^{-1} and D diagonal. The order of
the eigenvalues on the diagonal is not canonical, so the check below
multiplies the factors back together.

``` r

dg <- diagonalize(A)
dg$P %*% dg$D %*% inverse(dg$P) == A
#> [1] TRUE
sort(vapply(1:2, function(i) sstr(dg$D[i, i]), ""))
#> [1] "1" "3"
```

### Jordan normal form

Not every matrix is diagonalizable.
[`jordan_form()`](https://adamrsophie.github.io/symcompR/reference/jordan_form.md)
returns `list(P, J)` with M = P J P^{-1} and J in Jordan normal form.

``` r

J3 <- Matrix(list(list(1, 1, 0), list(0, 1, 0), list(0, 0, 2)))
is_diagonalizable(J3)
#> [1] FALSE
jordan_form(J3)
#> $P
#> Matrix(3, 3)
#> [0, 1, 0]
#> [0, 0, 1]
#> [1, 0, 0]
#> 
#> $J
#> Matrix(3, 3)
#> [2, 0, 0]
#> [0, 1, 1]
#> [0, 0, 1]
K <- Matrix(list(list(5, 4, 2, 1), list(0, 1, -1, -1),
                 list(-1, -1, 3, 0), list(1, 1, -1, 2)))
jordan_form(K, calc_transform = FALSE)
#> Matrix(4, 4)
#> [2, 0, 0, 0]
#> [0, 1, 0, 0]
#> [0, 0, 4, 1]
#> [0, 0, 0, 4]
```

The last matrix has eigenvalues 1, 2, and 4, with a 2 by 2 Jordan block
for the repeated eigenvalue 4.

[`matrix_exp()`](https://adamrsophie.github.io/symcompR/reference/matrix_exp.md)
computes the matrix exponential of a diagonalizable matrix:

``` r

matrix_exp(Matrix(list(list(0, 1), list(-1, 0))))
#> Matrix(2, 2)
#> [cos(1) , sin(1)]
#> [-sin(1), cos(1)]
```

## Decompositions

``` r

lu <- LUdecomposition(M)
lu$L
#> Matrix(2, 2)
#> [1, 0]
#> [3, 1]
lu$U
#> Matrix(2, 2)
#> [1, 2 ]
#> [0, -2]
qr <- QRdecomposition(M)
qr$Q
#> Matrix(2, 2)
#> [sqrt(10)/10  , 3*sqrt(10)/10]
#> [3*sqrt(10)/10, -sqrt(10)/10 ]
qr$R
#> Matrix(2, 2)
#> [sqrt(10), 7*sqrt(10)/5]
#> [0       , sqrt(10)/5  ]
SPD <- Matrix(list(list(4, 2), list(2, 3)))
cholesky(SPD)
#> Matrix(2, 2)
#> [2, 0      ]
#> [1, sqrt(2)]
LDLdecomposition(SPD)
#> $L
#> Matrix(2, 2)
#> [1  , 0]
#> [1/2, 1]
#> 
#> $D
#> Matrix(2, 2)
#> [4, 0]
#> [0, 2]
```

Check one of them: QR should give back M.

``` r

simplify(qr$Q %*% qr$R)
#> Matrix(2, 2)
#> [1, 2]
#> [3, 4]
```

## Integer normal forms

The Hermite normal form is the integer analogue of the reduced row
echelon form, and the Smith normal form diagonalizes an integer matrix
by unimodular row and column operations. Its diagonal, the invariant
factors, classifies finitely generated abelian groups.

``` r

hermite_normal_form(Matrix(list(list(2, 7, 17, 29, 41),
                                list(3, 11, 19, 31, 43),
                                list(5, 13, 23, 37, 47))))
#> Matrix(3, 3)
#> [1, 0, 0]
#> [0, 2, 1]
#> [0, 0, 1]
SN <- Matrix(list(list(12, 6, 4), list(3, 9, 6), list(2, 16, 14)))
smith_normal_form(SN)
#> Matrix(3, 3)
#> [1, 0 , 0 ]
#> [0, 10, 0 ]
#> [0, 0 , 30]
invariant_factors(SN)
#> (1, 10, 30)
```

## Calculus with matrices

``` r

jacobian(list(x^2 * y, 5 * x + sin(y)), list(x, y))
#> Matrix(2, 2)
#> [2*x*y, x^2   ]
#> [5    , cos(y)]
hessian(x^3 + x * y^2, list(x, y))
#> Matrix(2, 2)
#> [6*x, 2*y]
#> [2*y, 2*x]
wronskian(list(sin(x), cos(x)), x)
#> -sin(x)^2 - cos(x)^2
gradient(x^2 * y, list(x, y))
#> Matrix(2, 1)
#> [2*x*y]
#> [x^2  ]
```

## Vectors

Column vectors are one column matrices.

``` r

u <- Matrix(c(1, 2, 3)); v <- Matrix(c(4, 5, 6))
dot(u, v)
#> 32
cross(Matrix(c(1, 0, 0)), Matrix(c(0, 1, 0)))
#> Matrix(3, 1)
#> [0]
#> [0]
#> [1]
norm_matrix(Matrix(c(3, 4)))
#> 5
kronecker_product(M, eye(2))
#> Matrix(4, 4)
#> [1, 0, 2, 0]
#> [0, 1, 0, 2]
#> [3, 0, 4, 0]
#> [0, 3, 0, 4]
```

## Matrix expressions

`MatrixSymbol(name, n, m)` is a matrix whose entries are not specified.
Products, sums, powers, transposes, and inverses of matrix symbols form
matrix expressions, which respect noncommutativity and check shapes.

``` r

A <- MatrixSymbol("A", 2, 2); B <- MatrixSymbol("B", 2, 2)
A * B + A * B
#> 2*A*B
A * B - B * A
#> A*B - B*A
t(A * B)
#> B.T*A.T
A^-1 * A
#> I
mat_shape(A * B)
#> [[1]]
#> 2
#> 
#> [[2]]
#> 2
```

Here `*` is the matrix product, because matrix symbols have no
elementwise interpretation;
[`HadamardProduct()`](https://adamrsophie.github.io/symcompR/reference/matrix-expressions.md)
is the elementwise product.

``` r

HadamardProduct(A, B)
#> HadamardProduct(A, B)
MatrixElement(A, 0, 1)
#> A[0, 1]
as_explicit(A * B)
#> Matrix(2, 2)
#> [A[0, 0]*B[0, 0] + A[0, 1]*B[1, 0], A[0, 0]*B[0, 1] + A[0, 1]*B[1, 1]]
#> [A[1, 0]*B[0, 0] + A[1, 1]*B[1, 0], A[1, 0]*B[0, 1] + A[1, 1]*B[1, 1]]
```

[`MatrixElement()`](https://adamrsophie.github.io/symcompR/reference/matrix-expressions.md)
uses SymPy’s zero based indices, which is also how the entries of an
explicit matrix symbol print.
[`Trace()`](https://adamrsophie.github.io/symcompR/reference/matrix-expressions.md),
[`Determinant()`](https://adamrsophie.github.io/symcompR/reference/matrix-expressions.md),
and
[`Inverse()`](https://adamrsophie.github.io/symcompR/reference/matrix-expressions.md)
are unevaluated until
[`doit()`](https://adamrsophie.github.io/symcompR/reference/doit.md).

``` r

doit(Trace(A + B))
#> Trace(A) + Trace(B)
X <- MatrixSymbol("X", n, n)
X * Identity(n)
#> X
X + ZeroMatrix(n, n)
#> X
t(Inverse(X))
#> X.T^(-1)
MatrixSlice(X, list(0, 2), list(0, 2))
#> X[:2, :2]
```

Block matrices are built from matrix expressions.

``` r

BlockMatrix(list(list(A, ZeroMatrix(2, 2)), list(ZeroMatrix(2, 2), B)))
#> Matrix([
#> [A, 0],
#> [0, B]])
```

More matrix expression classes.
[`MatAdd()`](https://adamrsophie.github.io/symcompR/reference/matrix-expressions.md),
[`MatMul()`](https://adamrsophie.github.io/symcompR/reference/matrix-expressions.md),
and
[`MatPow()`](https://adamrsophie.github.io/symcompR/reference/matrix-expressions.md)
build unevaluated sums, products, and powers,
[`Transpose()`](https://adamrsophie.github.io/symcompR/reference/matrix-expressions.md)
and
[`Adjoint()`](https://adamrsophie.github.io/symcompR/reference/matrix-expressions.md)
the (conjugate) transposes,
[`HadamardPower()`](https://adamrsophie.github.io/symcompR/reference/matrix-expressions.md)
the elementwise power, and
[`KroneckerProduct()`](https://adamrsophie.github.io/symcompR/reference/matrix-expressions.md)
the Kronecker product.
[`matrix_symbols()`](https://adamrsophie.github.io/symcompR/reference/matrix-expressions.md)
lists the matrix symbols in an expression.

``` r

C <- MatrixSymbol("C", 2, 2)
MatAdd(A, B)
#> A + B
MatMul(2, A, B)
#> 2*A*B
MatPow(A, 3)
#> A^3
as_explicit(MatPow(A, 2))
#> Matrix(2, 2)
#> [A[0, 0]^2 + A[0, 1]*A[1, 0]      , A[0, 0]*A[0, 1] + A[0, 1]*A[1, 1]]
#> [A[0, 0]*A[1, 0] + A[1, 0]*A[1, 1], A[0, 1]*A[1, 0] + A[1, 1]^2      ]
Transpose(A * B)
#> (A*B).T
doit(Transpose(A * B))
#> B.T*A.T
Adjoint(A)
#> Adjoint(A)
HadamardPower(A, 2)
#> HadamardPower(A, 2)
mat_shape(KroneckerProduct(A, B))
#> [[1]]
#> 4
#> 
#> [[2]]
#> 4
matrix_symbols(A * B)
#> [[1]]
#> A
#> 
#> [[2]]
#> B
doit(Determinant(A))
#> Determinant(A)
```

Structured matrices:
[`DotProduct()`](https://adamrsophie.github.io/symcompR/reference/matrix-expressions.md)
of two column vectors, the diagonal of a matrix
([`DiagonalOf()`](https://adamrsophie.github.io/symcompR/reference/matrix-expressions.md)),
a diagonal matrix from a matrix or vector
([`DiagonalMatrix()`](https://adamrsophie.github.io/symcompR/reference/matrix-expressions.md),
[`DiagMatrix()`](https://adamrsophie.github.io/symcompR/reference/matrix-expressions.md),
[`diagonalize_vector()`](https://adamrsophie.github.io/symcompR/reference/matrix-expressions.md)),
the matrix of ones, matrices defined by a function of the indices,
permutation matrices, and permuted rows
([`MatrixPermute()`](https://adamrsophie.github.io/symcompR/reference/matrix-expressions.md)).

``` r

u <- MatrixSymbol("u", 3, 1); w <- MatrixSymbol("w", 3, 1)
DotProduct(u, w)
#> DotProduct(u, w)
doit(DotProduct(u, w))
#> u[0, 0]*w[0, 0] + u[1, 0]*w[1, 0] + u[2, 0]*w[2, 0]
as_explicit(DiagonalOf(A))
#> Matrix(2, 1)
#> [A[0, 0]]
#> [A[1, 1]]
as_explicit(DiagonalMatrix(A))
#> Matrix(2, 2)
#> [A[0, 0], 0      ]
#> [0      , A[1, 1]]
as_explicit(DiagMatrix(u))
#> Matrix(3, 3)
#> [u[0, 0], 0      , 0      ]
#> [0      , u[1, 0], 0      ]
#> [0      , 0      , u[2, 0]]
as_explicit(diagonalize_vector(u))
#> Matrix(3, 3)
#> [u[0, 0], 0      , 0      ]
#> [0      , u[1, 0], 0      ]
#> [0      , 0      , u[2, 0]]
as_explicit(OneMatrix(2, 3))
#> Matrix(2, 3)
#> [1, 1, 1]
#> [1, 1, 1]
as_explicit(FunctionMatrix(2, 2, Lambda(list(x, y), x + y)))
#> Matrix(2, 2)
#> [0, 1]
#> [1, 2]
as_explicit(PermutationMatrix(Permutation(c(1, 2, 0))))
#> Matrix(3, 3)
#> [0, 1, 0]
#> [0, 0, 1]
#> [1, 0, 0]
as_explicit(MatrixPermute(A, Permutation(c(1, 0))))
#> Matrix(2, 2)
#> [A[1, 0], A[1, 1]]
#> [A[0, 0], A[0, 1]]
doit(Permanent(Matrix(list(list(1, 2), list(3, 4)))))
#> 10
```

Block matrices:
[`BlockDiagMatrix()`](https://adamrsophie.github.io/symcompR/reference/matrix-expressions.md)
is block diagonal,
[`block_collapse()`](https://adamrsophie.github.io/symcompR/reference/matrix-expressions.md)
evaluates block arithmetic, and
[`blockcut()`](https://adamrsophie.github.io/symcompR/reference/matrix-expressions.md)
cuts a matrix into blocks.

``` r

BlockDiagMatrix(A, B)
#> Matrix([
#> [A, 0],
#> [0, B]])
block_collapse(BlockMatrix(list(list(A, B))) * BlockDiagMatrix(B, C))
#> Matrix([[A*B, B*C]])
block_collapse(BlockDiagMatrix(A, B) + BlockDiagMatrix(B, C))
#> Matrix([
#> [A + B,     0],
#> [    0, B + C]])
as_explicit(blockcut(Matrix(4, 4, f = 0:15), c(2, 2), c(2, 2)))
#> Matrix(4, 4)
#> [0 , 1 , 2 , 3 ]
#> [4 , 5 , 6 , 7 ]
#> [8 , 9 , 10, 11]
#> [12, 13, 14, 15]
```

## More matrix functions

``` r

M <- Matrix(list(list(1, 2), list(3, 4)))
transpose(M)
#> Matrix(2, 2)
#> [1, 3]
#> [2, 4]
det_quick(M)
#> -2
inv_quick(M)
#> Matrix(2, 2)
#> [-2 , 1   ]
#> [3/2, -1/2]
per(M)
#> 10
rowspace(Matrix(list(list(1, 2), list(2, 4))))[[1]]
#> Matrix(1, 2)
#> [1, 2]
set.seed(1)
randMatrix(2, 3)
#> Matrix(2, 3)
#> [67, 38, 0 ]
#> [33, 86, 42]
cofactor_matrix(M)
#> Matrix(2, 2)
#> [4 , -3]
#> [-2, 1 ]
minor(M, 1, 2)
#> 3
trace_matrix(M)
#> 5
is_square(M)
#> [1] TRUE
is_symmetric(M)
#> [1] FALSE
is_symmetric(Matrix(list(list(1, 2), list(2, 1))))
#> [1] TRUE
is_zero_matrix(zeros(2, 2))
#> [1] TRUE
vec(M)
#> Matrix(4, 1)
#> [1]
#> [3]
#> [2]
#> [4]
hadamard_product(M, M)
#> Matrix(2, 2)
#> [1, 4 ]
#> [9, 16]
matrix_multiply_elementwise(M, M)
#> Matrix(2, 2)
#> [1, 4 ]
#> [9, 16]
diag_matrix(1, 2, 3)
#> Matrix(3, 3)
#> [1, 0, 0]
#> [0, 2, 0]
#> [0, 0, 3]
adjoint(M)
#> Matrix(2, 2)
#> [1, 3]
#> [2, 4]
pinv(Matrix(list(list(1, 2, 3))))
#> Matrix(3, 1)
#> [1/14]
#> [1/7 ]
#> [3/14]
pinv(Matrix(list(list(1, 0), list(0, 2), list(0, 0))))
#> Matrix(2, 3)
#> [1, 0  , 0]
#> [0, 1/2, 0]
GramSchmidt(list(Matrix(list(1, 1, 0)), Matrix(list(1, 0, 1))))
#> [[1]]
#> Matrix(3, 1)
#> [1]
#> [1]
#> [0]
#> 
#> [[2]]
#> Matrix(3, 1)
#> [1/2 ]
#> [-1/2]
#> [1   ]
GramSchmidt(list(Matrix(list(1, 1, 0)), Matrix(list(1, 0, 1))), orthonormal = TRUE)
#> [[1]]
#> Matrix(3, 1)
#> [sqrt(2)/2]
#> [sqrt(2)/2]
#> [0        ]
#> 
#> [[2]]
#> Matrix(3, 1)
#> [sqrt(6)/6 ]
#> [-sqrt(6)/6]
#> [sqrt(6)/3 ]
```

`minor(M, i, j)` uses R’s indices, starting at 1.

Sparse and banded matrices:

``` r

SparseMatrix(2, 2, list(list(0, 0, 1), list(1, 1, 5)))
#> Matrix(2, 2)
#> [1, 0]
#> [0, 5]
ImmutableSparseMatrix(2, 2, list(list(0, 1, 3)))
#> Matrix(2, 2)
#> [0, 3]
#> [0, 0]
MutableSparseMatrix(2, 2, list(list(1, 0, 7)))
#> Matrix(2, 2)
#> [0, 0]
#> [7, 0]
banded(3, list("0" = 1, "-1" = c(2, 3)))
#> Matrix(3, 3)
#> [1, 0, 0]
#> [2, 1, 0]
#> [0, 3, 1]
```

Rotation matrices:
[`rot_axis1()`](https://adamrsophie.github.io/symcompR/reference/matrix-constructions.md)
to
[`rot_axis3()`](https://adamrsophie.github.io/symcompR/reference/matrix-constructions.md)
rotate the coordinate frame,
[`rot_ccw_axis1()`](https://adamrsophie.github.io/symcompR/reference/matrix-constructions.md)
to
[`rot_ccw_axis3()`](https://adamrsophie.github.io/symcompR/reference/matrix-constructions.md)
rotate vectors counterclockwise, and
[`rot_givens()`](https://adamrsophie.github.io/symcompR/reference/matrix-builders.md)
is a Givens rotation in a coordinate plane.

``` r

rot_axis1(sym_pi / 2)
#> Matrix(3, 3)
#> [1, 0 , 0]
#> [0, 0 , 1]
#> [0, -1, 0]
rot_axis2(sym_pi / 2)
#> Matrix(3, 3)
#> [0, 0, -1]
#> [0, 1, 0 ]
#> [1, 0, 0 ]
rot_ccw_axis1(sym_pi / 2)
#> Matrix(3, 3)
#> [1, 0, 0 ]
#> [0, 0, -1]
#> [0, 1, 0 ]
rot_ccw_axis2(sym_pi / 2)
#> Matrix(3, 3)
#> [0 , 0, 1]
#> [0 , 1, 0]
#> [-1, 0, 0]
rot_ccw_axis3(sym_pi / 2)
#> Matrix(3, 3)
#> [0, -1, 0]
#> [1, 0 , 0]
#> [0, 0 , 1]
rot_givens(0, 1, sym_pi / 3)
#> Matrix(3, 3)
#> [1/2       , sqrt(3)/2, 0]
#> [-sqrt(3)/2, 1/2      , 0]
#> [0         , 0        , 1]
```

Linear systems and other constructions:
[`linear_eq_to_matrix()`](https://adamrsophie.github.io/symcompR/reference/solver-helpers.md)
extracts the coefficient matrix and right hand side,
[`casoratian()`](https://adamrsophie.github.io/symcompR/reference/matrix-constructions.md)
is the Casoratian of sequences (the discrete Wronskian), and
[`symarray()`](https://adamrsophie.github.io/symcompR/reference/matrix-constructions.md)
is a matrix of indexed symbols.

``` r

k <- Symbol("k", integer = TRUE)
linear_eq_to_matrix(list(x + 2 * y - 3, 3 * x - y + 1), list(x, y))
#> $A
#> Matrix(2, 2)
#> [1, 2 ]
#> [3, -1]
#> 
#> $b
#> Matrix(2, 1)
#> [3 ]
#> [-1]
casoratian(list(1, 2^k), k)
#> 1
symarray("a", c(2, 2))
#> Matrix(2, 2)
#> [a_0_0, a_0_1]
#> [a_1_0, a_1_1]
```
