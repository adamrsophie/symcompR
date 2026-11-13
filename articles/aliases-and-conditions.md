# Aliases and error conditions

This reference collects the exports that do not need their own example:
alternative names for functions documented elsewhere, and the condition
classes that mirror SymPy’s exception hierarchy.

``` r

library(symcompR)
x <- Symbol("x"); y <- Symbol("y")
```

## Alternative names

SymPy often offers several names for the same object, and symcompR keeps
them so that code translated from Python works unchanged.

| Alias | Same as |
|----|----|
| [`RootOf()`](https://adamrsophie.github.io/symcompR/reference/poly-aliases.md), [`ComplexRootOf()`](https://adamrsophie.github.io/symcompR/reference/poly-aliases.md) | [`rootof()`](https://adamrsophie.github.io/symcompR/reference/poly-aliases.md), which calls [`CRootOf()`](https://adamrsophie.github.io/symcompR/reference/CRootOf.md) |
| [`Eijk()`](https://adamrsophie.github.io/symcompR/reference/function-aliases.md) | [`LeviCivita()`](https://adamrsophie.github.io/symcompR/reference/function-aliases.md) |
| [`ImmutableDenseMatrix()`](https://adamrsophie.github.io/symcompR/reference/matrix-constructions.md), [`MutableMatrix()`](https://adamrsophie.github.io/symcompR/reference/matrix-constructions.md), [`MutableDenseMatrix()`](https://adamrsophie.github.io/symcompR/reference/matrix-constructions.md) | [`ImmutableMatrix()`](https://adamrsophie.github.io/symcompR/reference/matrix-constructions.md), which calls [`Matrix()`](https://adamrsophie.github.io/symcompR/reference/Matrix.md) |
| [`DenseNDimArray()`](https://adamrsophie.github.io/symcompR/reference/array-aliases.md), [`ImmutableDenseNDimArray()`](https://adamrsophie.github.io/symcompR/reference/array-aliases.md), [`ImmutableSparseNDimArray()`](https://adamrsophie.github.io/symcompR/reference/array-aliases.md), [`MutableDenseNDimArray()`](https://adamrsophie.github.io/symcompR/reference/array-aliases.md), [`MutableSparseNDimArray()`](https://adamrsophie.github.io/symcompR/reference/array-aliases.md), [`SparseNDimArray()`](https://adamrsophie.github.io/symcompR/reference/array-aliases.md) | [`NDimArray()`](https://adamrsophie.github.io/symcompR/reference/array-aliases.md) |
| [`GroebnerBasis()`](https://adamrsophie.github.io/symcompR/reference/poly-aliases.md) | [`groebner()`](https://adamrsophie.github.io/symcompR/reference/groebner.md) |
| `true`, `false` | `S_true`, `S_false` |
| `lex`, `grlex`, `grevlex` | the strings `"lex"`, `"grlex"`, and `"grevlex"`, monomial order names |
| [`deg()`](https://adamrsophie.github.io/symcompR/reference/geometry-aliases.md), [`rad()`](https://adamrsophie.github.io/symcompR/reference/geometry-aliases.md) | conversion between radians and degrees |

``` r

c(identical(RootOf, rootof), identical(ComplexRootOf, rootof),
  identical(Eijk, LeviCivita))
#> [1] TRUE TRUE TRUE
RootOf(x^3 - x - 1, x, 0)
#> CRootOf(x^3 - x - 1, 0)
LeviCivita(1, 2, 3)
#> 1
Eijk(1, 3, 2)
#> -1
c(identical(ImmutableDenseMatrix, ImmutableMatrix),
  identical(MutableMatrix, ImmutableMatrix),
  identical(MutableDenseMatrix, ImmutableMatrix))
#> [1] TRUE TRUE TRUE
vapply(list(DenseNDimArray, ImmutableDenseNDimArray, ImmutableSparseNDimArray,
            MutableDenseNDimArray, MutableSparseNDimArray, SparseNDimArray),
       identical, NA, NDimArray)
#> [1] TRUE TRUE TRUE TRUE TRUE TRUE
vapply(GroebnerBasis(list(x * y - 2 * y, 2 * y^2 - x^2), x, y, order = "lex"),
       sstr, "")
#> [1] "x^2 - 2*y^2" "x*y - 2*y"   "y^3 - 2*y"
true
#> True
false
#> False
c(lex, grlex, grevlex)
#> [1] "lex"     "grlex"   "grevlex"
rad(180)
#> pi
deg(sym_pi)
#> 180
```

## Error conditions

Each SymPy exception class has a constructor that returns an R error
condition. Its class vector follows SymPy’s inheritance, so a handler in
[`tryCatch()`](https://rdrr.io/r/base/conditions.html) may name a base
class: a `NonSquareMatrixError` is caught by `ShapeError` and by
`MatrixError`, and a `MultivariatePolynomialError` by `PolynomialError`
and by `BasePolynomialError`. Every condition also has the class
`"sympy_error"`.

``` r

tryCatch(sympify("x +* 1"), SympifyError = function(e) "not an expression")
#> [1] "not an expression"
tryCatch(stop(MultivariatePolynomialError("x*y")),
         PolynomialError = function(e) conditionMessage(e))
#> [1] "x*y"
tryCatch(sym_det(Matrix(list(list(1, 2)))),
         MatrixError = function(e) class(e)[1])
#> [1] "NonSquareMatrixError"
```

The table lists every condition class with its parent classes.

``` r

conds <- c("ArgumentIndexError", "BasePolynomialError", "CoercionFailed",
  "ComputationFailed", "DomainError", "EvaluationFailed",
  "ExactQuotientFailed", "ExtraneousFactors", "FlagError", "GeneratorsError",
  "GeneratorsNeeded", "GeometryError", "HeuristicGCDFailed",
  "HomomorphismFailed", "IsomorphismFailed", "MatrixError",
  "ModularGCDFailed", "MultivariatePolynomialError", "NonInvertibleMatrixError",
  "NonPositiveDefiniteMatrixError", "NonSquareMatrixError", "NotAlgebraic",
  "NotInvertible", "NotReversible", "OperationNotSupported", "OptionError",
  "PoleError", "PolificationFailed", "PolynomialDivisionFailed",
  "PolynomialError", "PrecisionExhausted", "RefinementFailed", "ShapeError",
  "SympifyError", "UnificationFailed", "UnivariatePolynomialError",
  "UnsolvableFactorError")
parents <- vapply(conds, function(nm) {
  cl <- class(do.call(nm, list("")))
  paste(setdiff(cl, c(nm, "sympy_error", "error", "condition")), collapse = ", ")
}, "")
data.frame(condition = conds, parents = parents, row.names = NULL)
#>                         condition                              parents
#> 1              ArgumentIndexError                           ValueError
#> 2             BasePolynomialError                                     
#> 3                  CoercionFailed                  BasePolynomialError
#> 4               ComputationFailed                  BasePolynomialError
#> 5                     DomainError                  BasePolynomialError
#> 6                EvaluationFailed                  BasePolynomialError
#> 7             ExactQuotientFailed                  BasePolynomialError
#> 8               ExtraneousFactors                  BasePolynomialError
#> 9                       FlagError     OptionError, BasePolynomialError
#> 10                GeneratorsError                  BasePolynomialError
#> 11               GeneratorsNeeded GeneratorsError, BasePolynomialError
#> 12                  GeometryError                           ValueError
#> 13             HeuristicGCDFailed                  BasePolynomialError
#> 14             HomomorphismFailed                  BasePolynomialError
#> 15              IsomorphismFailed                  BasePolynomialError
#> 16                    MatrixError                                     
#> 17               ModularGCDFailed                  BasePolynomialError
#> 18    MultivariatePolynomialError PolynomialError, BasePolynomialError
#> 19       NonInvertibleMatrixError              ValueError, MatrixError
#> 20 NonPositiveDefiniteMatrixError              ValueError, MatrixError
#> 21           NonSquareMatrixError  ShapeError, ValueError, MatrixError
#> 22                   NotAlgebraic                  BasePolynomialError
#> 23                  NotInvertible                  BasePolynomialError
#> 24                  NotReversible                  BasePolynomialError
#> 25          OperationNotSupported                  BasePolynomialError
#> 26                    OptionError                  BasePolynomialError
#> 27                      PoleError                                     
#> 28             PolificationFailed PolynomialError, BasePolynomialError
#> 29       PolynomialDivisionFailed                  BasePolynomialError
#> 30                PolynomialError                  BasePolynomialError
#> 31             PrecisionExhausted                      ArithmeticError
#> 32               RefinementFailed                  BasePolynomialError
#> 33                     ShapeError              ValueError, MatrixError
#> 34                   SympifyError                           ValueError
#> 35              UnificationFailed                  BasePolynomialError
#> 36      UnivariatePolynomialError PolynomialError, BasePolynomialError
#> 37          UnsolvableFactorError                  BasePolynomialError
```
