# Kinds of expressions, after sympy/core/kind.py, sympy/matrices/kind.py,
# sympy/tensor/array/array_derivatives.py (ArrayKind), sympy/sets/sets.py
# (SetKind), sympy/core/containers.py (TupleKind), sympy/vector/kind.py and
# sympy/physics/quantum/kind.py. A kind is a small S3 object with a name
# and, for the parametrized kinds, the kinds of the elements.

.new_kind <- function(name, elems = list()) structure(list(name = name, elems = elems), class = "sym_kind")

#' Kinds of expressions
#'
#' `expr_kind(e)` is SymPy's `e.kind`: the mathematical kind of an object,
#' independent of its class. Numbers, commutative symbols and arithmetic of
#' them are `NumberKind`; Booleans and relations `BooleanKind`; matrices
#' `MatrixKind(element kind)`; arrays `ArrayKind(...)`; sets
#' `SetKind(...)`; tuples `TupleKind(...)`; vectors `VectorKind(...)`;
#' quantum states and operators `KetKind`, `BraKind` and `OperatorKind`;
#' anything else, including applied functions other than `exp`,
#' `UndefinedKind`. (`kind()` is the Fortran intrinsic of the code
#' generation nodes.)
#'
#' @param e An expression, matrix, array, set or vector.
#' @param element_kind The kind of the elements.
#' @param ... Kinds of the tuple entries.
#' @return A kind.
#' @examples
#' x <- Symbol("x")
#' expr_kind(x + 1)
#' expr_kind(sin(x))
#' expr_kind(x > 1)
#' expr_kind(Matrix(list(1, x)))
#' expr_kind(FiniteSet(1, x))
#' MatrixKind() == MatrixKind(NumberKind)
#' @name kinds
#' @export
expr_kind <- function(e) {
  if (inherits(e, "sym_kind")) return(e)
  if (inherits(e, "symmatrix")) return(MatrixKind(.common_kind(lapply(e$data, expr_kind))))
  if (inherits(e, "SymArray")) return(ArrayKind(.common_kind(lapply(unclass(e)$data, expr_kind))))
  if (inherits(e, "sym_vec")) return(VectorKind())
  if (is.logical(e)) return(BooleanKind)
  if (is.numeric(e)) return(NumberKind)
  if (!is_sym(e)) return(UndefinedKind)
  t <- e$type
  if (t %in% c("Integer", "Rational", "Float", "Constant", "NumberSymbol", "ImaginaryUnit", "Infinity",
               "NegativeInfinity", "ComplexInfinity", "NaN", "AlgebraicNumber", "BaseScalar"))
    return(NumberKind)
  if (t %in% c("Symbol", "Dummy", "Wild")) return(if (isFALSE(ask_fact(e, "commutative"))) UndefinedKind else NumberKind)
  if (t %in% c(.REL_TYPES, .BOOL_TYPES, "BooleanTrue", "BooleanFalse", "Contains"))
    return(BooleanKind)
  if (.is_matexpr(e)) return(MatrixKind(NumberKind))
  if (t %in% .SET_TYPES) {
    if (t == "EmptySet") return(SetKind())
    if (t == "FiniteSet") return(SetKind(.common_kind(lapply(e$args, expr_kind))))
    return(SetKind(NumberKind))
  }
  if (t == "Tuple") return(do.call(TupleKind, lapply(e$args, expr_kind)))
  if (t %in% c("Ket", "Bra", "Operator") || length(.qmro(e))) {
    if (.qisa(e, "KetBase")) return(KetKind)
    if (.qisa(e, "BraBase")) return(BraKind)
    if (.qisa(e, "Operator")) return(OperatorKind)
  }
  if (t == "Add") {
    ks <- lapply(e$args, expr_kind)
    return(if (all(vapply(ks, function(k) k == ks[[1L]], NA))) ks[[1L]] else UndefinedKind)
  }
  if (t == "Mul") return(Reduce(.mul_kind, lapply(e$args, expr_kind)))
  if (t == "Pow") {
    kb <- expr_kind(e$args[[1L]]); ke <- expr_kind(e$args[[2L]])
    if (ke != NumberKind) return(UndefinedKind)
    return(if (kb == NumberKind || kb$name == "MatrixKind") kb else UndefinedKind)
  }
  if (t == "Function" && e$name == "exp") return(expr_kind(e$args[[1L]]))
  if (t == "Function" && e$name == "Mod") return(.common_kind(lapply(e$args, expr_kind)))
  if (t %in% c("Integral", "Sum", "Product", "Derivative")) return(expr_kind(e$args[[1L]]))
  UndefinedKind
}

# Mul._kind_dispatch: numbers scale every kind; any other pair is undefined
.mul_kind <- function(a, b) {
  if (a == NumberKind) return(b)
  if (b == NumberKind) return(a)
  UndefinedKind
}
.common_kind <- function(ks) {
  if (!length(ks)) return(UndefinedKind)
  if (all(vapply(ks, function(k) k == ks[[1L]], NA))) ks[[1L]] else UndefinedKind
}

#' @rdname kinds
#' @format The kinds without parameters are objects.
#' @export
NumberKind <- .new_kind("NumberKind")
#' @rdname kinds
#' @export
UndefinedKind <- .new_kind("UndefinedKind")
#' @rdname kinds
#' @export
BooleanKind <- .new_kind("BooleanKind")
#' @rdname kinds
#' @export
KetKind <- .new_kind("KetKind")
#' @rdname kinds
#' @export
BraKind <- .new_kind("BraKind")
#' @rdname kinds
#' @export
OperatorKind <- .new_kind("OperatorKind")
#' @rdname kinds
#' @export
MatrixKind <- function(element_kind = NumberKind) .new_kind("MatrixKind", list(element_kind))
#' @rdname kinds
#' @export
ArrayKind <- function(element_kind = NumberKind) .new_kind("ArrayKind", list(element_kind))
#' @rdname kinds
#' @export
VectorKind <- function(element_kind = NumberKind) .new_kind("VectorKind", list(element_kind))
#' @rdname kinds
#' @export
SetKind <- function(element_kind = NULL) .new_kind("SetKind", if (is.null(element_kind)) list() else list(element_kind))
#' @rdname kinds
#' @export
TupleKind <- function(...) .new_kind("TupleKind", list(...))

#' @export
format.sym_kind <- function(x, ...) {
  if (!length(x$elems) && !x$name %in% c("SetKind", "TupleKind")) return(x$name)
  paste0(x$name, "(", paste(vapply(x$elems, format, ""), collapse = ", "), ")")
}
#' @export
print.sym_kind <- function(x, ...) { cat(format(x), "\n"); invisible(x) }
#' @export
Ops.sym_kind <- function(e1, e2) {
  eq <- inherits(e1, "sym_kind") && inherits(e2, "sym_kind") && identical(format(e1), format(e2))
  switch(.Generic, `==` = eq, `!=` = !eq, stop("operator ", .Generic, " is not defined for kinds", call. = FALSE))
}
