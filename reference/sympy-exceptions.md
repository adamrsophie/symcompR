# SymPy exceptions

Constructors for SymPy's exception classes. Each returns an R error
condition for [`stop()`](https://rdrr.io/r/base/stop.html), classed with
the exception name and its SymPy base classes, so handlers in
[`tryCatch()`](https://rdrr.io/r/base/conditions.html) may name a base
class: a `MultivariatePolynomialError` is caught by `PolynomialError`
and by `BasePolynomialError`, and a `NonSquareMatrixError` by
`ShapeError` and `MatrixError`. All carry the class `"sympy_error"`. The
package raises `SympifyError` from
[`sympify()`](https://adamrsophie.github.io/symcompR/reference/sympify.md),
`PoleError` from series expansions at singular points, and
`NonSquareMatrixError` and `NonInvertibleMatrixError` from matrix
inversion and determinants.

## Usage

``` r
SympifyError(message = "")

PoleError(message = "")

PrecisionExhausted(message = "")

ArgumentIndexError(message = "")

GeometryError(message = "")

MatrixError(message = "")

ShapeError(message = "")

NonSquareMatrixError(message = "")

NonInvertibleMatrixError(message = "")

NonPositiveDefiniteMatrixError(message = "")

BasePolynomialError(message = "")

ExactQuotientFailed(message = "")

PolynomialDivisionFailed(message = "")

OperationNotSupported(message = "")

HeuristicGCDFailed(message = "")

ModularGCDFailed(message = "")

HomomorphismFailed(message = "")

IsomorphismFailed(message = "")

ExtraneousFactors(message = "")

EvaluationFailed(message = "")

RefinementFailed(message = "")

CoercionFailed(message = "")

NotInvertible(message = "")

NotReversible(message = "")

NotAlgebraic(message = "")

DomainError(message = "")

PolynomialError(message = "")

UnificationFailed(message = "")

UnsolvableFactorError(message = "")

GeneratorsError(message = "")

GeneratorsNeeded(message = "")

ComputationFailed(message = "")

UnivariatePolynomialError(message = "")

MultivariatePolynomialError(message = "")

PolificationFailed(message = "")

OptionError(message = "")

FlagError(message = "")
```

## Arguments

- message:

  The error message.

## Value

An error condition.

## Examples

``` r
r <- tryCatch(sympify("x +* 1"), SympifyError = function(e) "not an expression")
r
#> [1] "not an expression"
tryCatch(stop(MultivariatePolynomialError("x*y")),
         PolynomialError = function(e) conditionMessage(e))
#> [1] "x*y"
```
