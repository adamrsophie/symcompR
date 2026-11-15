# SymPy's exception classes as R conditions, after sympy/polys/polyerrors.py,
# sympy/core/sympify.py, sympy/core/function.py, sympy/core/evalf.py,
# sympy/geometry/exceptions.py, sympy/matrices/exceptions.py and
# sympy/holonomic/holonomicerrors.py. Each
# constructor returns an error condition whose classes follow SymPy's
# inheritance, so tryCatch(..., PolynomialError = h) also catches
# MultivariatePolynomialError.

.SYMPY_EXCEPTIONS <- list(
  SympifyError = "ValueError", PoleError = character(), PrecisionExhausted = "ArithmeticError",
  ArgumentIndexError = "ValueError", GeometryError = "ValueError",
  MatrixError = character(), ShapeError = c("ValueError", "MatrixError"),
  NonSquareMatrixError = c("ShapeError", "ValueError", "MatrixError"),
  NonInvertibleMatrixError = c("ValueError", "MatrixError"),
  NonPositiveDefiniteMatrixError = c("ValueError", "MatrixError"),
  BasePolynomialError = character(),
  ExactQuotientFailed = "BasePolynomialError", PolynomialDivisionFailed = "BasePolynomialError",
  OperationNotSupported = "BasePolynomialError", HeuristicGCDFailed = "BasePolynomialError",
  ModularGCDFailed = "BasePolynomialError", HomomorphismFailed = "BasePolynomialError",
  IsomorphismFailed = "BasePolynomialError", ExtraneousFactors = "BasePolynomialError",
  EvaluationFailed = "BasePolynomialError", RefinementFailed = "BasePolynomialError",
  CoercionFailed = "BasePolynomialError", NotInvertible = "BasePolynomialError",
  NotReversible = "BasePolynomialError", NotAlgebraic = "BasePolynomialError",
  DomainError = "BasePolynomialError", PolynomialError = "BasePolynomialError",
  UnificationFailed = "BasePolynomialError", UnsolvableFactorError = "BasePolynomialError",
  GeneratorsError = "BasePolynomialError", GeneratorsNeeded = c("GeneratorsError", "BasePolynomialError"),
  ComputationFailed = "BasePolynomialError",
  UnivariatePolynomialError = c("PolynomialError", "BasePolynomialError"),
  MultivariatePolynomialError = c("PolynomialError", "BasePolynomialError"),
  PolificationFailed = c("PolynomialError", "BasePolynomialError"),
  OptionError = "BasePolynomialError", FlagError = c("OptionError", "BasePolynomialError"),
  BaseHolonomicError = character(), NotHyperSeriesError = "BaseHolonomicError",
  # sympy/polys/matrices/exceptions.py, sympy/physics/quantum,
  # sympy/polys/numberfields/exceptions.py and sympy/parsing/smtlib.py
  DMError = character(),
  DMBadInputError = "DMError",
  DMDomainError = "DMError",
  DMFormatError = "DMError",
  DMRankError = "DMError",
  DMShapeError = "DMError",
  DMNotAField = c("DMDomainError", "DMError"),
  DMNonInvertibleMatrixError = "DMError",
  DMNonSquareMatrixError = c("DMShapeError", "DMError"),
  DMValueError = "DMError",
  QuantumError = character(),
  HilbertSpaceError = "QuantumError",
  ClosureFailure = character(),
  StructureError = character(),
  MissingUnityError = "StructureError",
  SMTLibSyntaxError = character(),
  UnknownSMTLibCommandError = character(),
  UnknownSMTLibOperatorError = character())

.sympy_condition <- function(name, message) {
  structure(list(message = paste(as.character(message), collapse = ""), call = NULL),
            class = c(name, .SYMPY_EXCEPTIONS[[name]], "sympy_error", "error", "condition"))
}

#' SymPy exceptions
#'
#' Constructors for SymPy's exception classes. Each returns an R error
#' condition for [stop()], classed with the exception name and its SymPy
#' base classes, so handlers in [tryCatch()] may name a base class: a
#' `MultivariatePolynomialError` is caught by `PolynomialError` and by
#' `BasePolynomialError`, and a `NonSquareMatrixError` by `ShapeError` and
#' `MatrixError`. All carry the class `"sympy_error"`. The package raises
#' `SympifyError` from [sympify()], `PoleError` from series expansions at
#' singular points, and `NonSquareMatrixError` and `NonInvertibleMatrixError`
#' from matrix inversion and determinants.
#'
#' @param message The error message.
#' @return An error condition.
#' @examples
#' r <- tryCatch(sympify("x +* 1"), SympifyError = function(e) "not an expression")
#' r
#' tryCatch(stop(MultivariatePolynomialError("x*y")),
#'          PolynomialError = function(e) conditionMessage(e))
#' @name sympy-exceptions
NULL

#' @rdname sympy-exceptions
#' @export
SympifyError <- function(message = "") .sympy_condition("SympifyError", message)

#' @rdname sympy-exceptions
#' @export
PoleError <- function(message = "") .sympy_condition("PoleError", message)

#' @rdname sympy-exceptions
#' @export
PrecisionExhausted <- function(message = "") .sympy_condition("PrecisionExhausted", message)

#' @rdname sympy-exceptions
#' @export
ArgumentIndexError <- function(message = "") .sympy_condition("ArgumentIndexError", message)

#' @rdname sympy-exceptions
#' @export
GeometryError <- function(message = "") .sympy_condition("GeometryError", message)

#' @rdname sympy-exceptions
#' @export
MatrixError <- function(message = "") .sympy_condition("MatrixError", message)

#' @rdname sympy-exceptions
#' @export
ShapeError <- function(message = "") .sympy_condition("ShapeError", message)

#' @rdname sympy-exceptions
#' @export
NonSquareMatrixError <- function(message = "") .sympy_condition("NonSquareMatrixError", message)

#' @rdname sympy-exceptions
#' @export
NonInvertibleMatrixError <- function(message = "") .sympy_condition("NonInvertibleMatrixError", message)

#' @rdname sympy-exceptions
#' @export
NonPositiveDefiniteMatrixError <- function(message = "") .sympy_condition("NonPositiveDefiniteMatrixError", message)

#' @rdname sympy-exceptions
#' @export
BasePolynomialError <- function(message = "") .sympy_condition("BasePolynomialError", message)

#' @rdname sympy-exceptions
#' @export
ExactQuotientFailed <- function(message = "") .sympy_condition("ExactQuotientFailed", message)

#' @rdname sympy-exceptions
#' @export
PolynomialDivisionFailed <- function(message = "") .sympy_condition("PolynomialDivisionFailed", message)

#' @rdname sympy-exceptions
#' @export
OperationNotSupported <- function(message = "") .sympy_condition("OperationNotSupported", message)

#' @rdname sympy-exceptions
#' @export
HeuristicGCDFailed <- function(message = "") .sympy_condition("HeuristicGCDFailed", message)

#' @rdname sympy-exceptions
#' @export
ModularGCDFailed <- function(message = "") .sympy_condition("ModularGCDFailed", message)

#' @rdname sympy-exceptions
#' @export
HomomorphismFailed <- function(message = "") .sympy_condition("HomomorphismFailed", message)

#' @rdname sympy-exceptions
#' @export
IsomorphismFailed <- function(message = "") .sympy_condition("IsomorphismFailed", message)

#' @rdname sympy-exceptions
#' @export
ExtraneousFactors <- function(message = "") .sympy_condition("ExtraneousFactors", message)

#' @rdname sympy-exceptions
#' @export
EvaluationFailed <- function(message = "") .sympy_condition("EvaluationFailed", message)

#' @rdname sympy-exceptions
#' @export
RefinementFailed <- function(message = "") .sympy_condition("RefinementFailed", message)

#' @rdname sympy-exceptions
#' @export
CoercionFailed <- function(message = "") .sympy_condition("CoercionFailed", message)

#' @rdname sympy-exceptions
#' @export
NotInvertible <- function(message = "") .sympy_condition("NotInvertible", message)

#' @rdname sympy-exceptions
#' @export
NotReversible <- function(message = "") .sympy_condition("NotReversible", message)

#' @rdname sympy-exceptions
#' @export
NotAlgebraic <- function(message = "") .sympy_condition("NotAlgebraic", message)

#' @rdname sympy-exceptions
#' @export
DomainError <- function(message = "") .sympy_condition("DomainError", message)

#' @rdname sympy-exceptions
#' @export
PolynomialError <- function(message = "") .sympy_condition("PolynomialError", message)

#' @rdname sympy-exceptions
#' @export
UnificationFailed <- function(message = "") .sympy_condition("UnificationFailed", message)

#' @rdname sympy-exceptions
#' @export
UnsolvableFactorError <- function(message = "") .sympy_condition("UnsolvableFactorError", message)

#' @rdname sympy-exceptions
#' @export
GeneratorsError <- function(message = "") .sympy_condition("GeneratorsError", message)

#' @rdname sympy-exceptions
#' @export
GeneratorsNeeded <- function(message = "") .sympy_condition("GeneratorsNeeded", message)

#' @rdname sympy-exceptions
#' @export
ComputationFailed <- function(message = "") .sympy_condition("ComputationFailed", message)

#' @rdname sympy-exceptions
#' @export
UnivariatePolynomialError <- function(message = "") .sympy_condition("UnivariatePolynomialError", message)

#' @rdname sympy-exceptions
#' @export
MultivariatePolynomialError <- function(message = "") .sympy_condition("MultivariatePolynomialError", message)

#' @rdname sympy-exceptions
#' @export
PolificationFailed <- function(message = "") .sympy_condition("PolificationFailed", message)

#' @rdname sympy-exceptions
#' @export
OptionError <- function(message = "") .sympy_condition("OptionError", message)

#' @rdname sympy-exceptions
#' @export
FlagError <- function(message = "") .sympy_condition("FlagError", message)

#' @rdname sympy-exceptions
#' @export
DMError <- function(message = "") .sympy_condition("DMError", message)

#' @rdname sympy-exceptions
#' @export
DMBadInputError <- function(message = "") .sympy_condition("DMBadInputError", message)

#' @rdname sympy-exceptions
#' @export
DMDomainError <- function(message = "") .sympy_condition("DMDomainError", message)

#' @rdname sympy-exceptions
#' @export
DMFormatError <- function(message = "") .sympy_condition("DMFormatError", message)

#' @rdname sympy-exceptions
#' @export
DMRankError <- function(message = "") .sympy_condition("DMRankError", message)

#' @rdname sympy-exceptions
#' @export
DMShapeError <- function(message = "") .sympy_condition("DMShapeError", message)

#' @rdname sympy-exceptions
#' @export
DMNotAField <- function(message = "") .sympy_condition("DMNotAField", message)

#' @rdname sympy-exceptions
#' @export
DMNonInvertibleMatrixError <- function(message = "") .sympy_condition("DMNonInvertibleMatrixError", message)

#' @rdname sympy-exceptions
#' @export
DMNonSquareMatrixError <- function(message = "") .sympy_condition("DMNonSquareMatrixError", message)

#' @rdname sympy-exceptions
#' @export
DMValueError <- function(message = "") .sympy_condition("DMValueError", message)

#' @rdname sympy-exceptions
#' @export
QuantumError <- function(message = "") .sympy_condition("QuantumError", message)

#' @rdname sympy-exceptions
#' @export
HilbertSpaceError <- function(message = "") .sympy_condition("HilbertSpaceError", message)

#' @rdname sympy-exceptions
#' @export
ClosureFailure <- function(message = "") .sympy_condition("ClosureFailure", message)

#' @rdname sympy-exceptions
#' @export
StructureError <- function(message = "") .sympy_condition("StructureError", message)

#' @rdname sympy-exceptions
#' @export
MissingUnityError <- function(message = "") .sympy_condition("MissingUnityError", message)

#' @rdname sympy-exceptions
#' @export
SMTLibSyntaxError <- function(message = "") .sympy_condition("SMTLibSyntaxError", message)

#' @rdname sympy-exceptions
#' @export
UnknownSMTLibCommandError <- function(message = "") .sympy_condition("UnknownSMTLibCommandError", message)

#' @rdname sympy-exceptions
#' @export
UnknownSMTLibOperatorError <- function(message = "") .sympy_condition("UnknownSMTLibOperatorError", message)
