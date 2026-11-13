# TensorFlow and PyTorch code printers, after sympy/printing/tensorflow.py
# and sympy/printing/pytorch.py. Both derive from ArrayPrinter (in
# R/codegen-printer-pythonarray.R) and AbstractPythonCodePrinter, and print
# functions through a table from SymPy class names to library functions.

.prt_tf_mapping <- c(
  Abs = "tensorflow.math.abs", sign = "tensorflow.math.sign", ceiling = "tensorflow.math.ceil",
  floor = "tensorflow.math.floor", log = "tensorflow.math.log", exp = "tensorflow.math.exp",
  Sqrt = "tensorflow.math.sqrt", cos = "tensorflow.math.cos", acos = "tensorflow.math.acos",
  sin = "tensorflow.math.sin", asin = "tensorflow.math.asin", tan = "tensorflow.math.tan",
  atan = "tensorflow.math.atan", atan2 = "tensorflow.math.atan2", cosh = "tensorflow.math.cosh",
  acosh = "tensorflow.math.acosh", sinh = "tensorflow.math.sinh", asinh = "tensorflow.math.asinh",
  tanh = "tensorflow.math.tanh", atanh = "tensorflow.math.atanh", re = "tensorflow.math.real",
  im = "tensorflow.math.imag", arg = "tensorflow.math.angle", erf = "tensorflow.math.erf",
  loggamma = "tensorflow.math.lgamma", Equality = "tensorflow.math.equal",
  Unequality = "tensorflow.math.not_equal", StrictGreaterThan = "tensorflow.math.greater",
  StrictLessThan = "tensorflow.math.less", LessThan = "tensorflow.math.less_equal",
  GreaterThan = "tensorflow.math.greater_equal", And = "tensorflow.math.logical_and",
  Or = "tensorflow.math.logical_or", Not = "tensorflow.math.logical_not", Max = "tensorflow.math.maximum",
  Min = "tensorflow.math.minimum", MatAdd = "tensorflow.math.add", HadamardProduct = "tensorflow.math.multiply",
  Trace = "tensorflow.linalg.trace", Determinant = "tensorflow.linalg.det")

# the class name SymPy keys the mapping with: type(expr)
.prt_type_name <- function(e) {
  if (is_sym(e) && e$type == "Function" && !is.null(.cgp_fun_names[e$name]) && !is.na(.cgp_fun_names[e$name]))
    return(.cgp_fun_names[[e$name]])
  .cgp_clsname(e)
}

# TensorflowPrinter._print_Function (and the aliases set to it): the mapped
# function, folded over the arguments when there are several
.prt_print_mapped <- function(cls) function(p, e) {
  op <- .cgp_a(p, "mapping")[.prt_type_name(e)]
  if (is.na(op)) return(.cgp_super(p, cls, "_print_Basic", e))
  children <- lapply(e$args, function(a) .cgpp(p, a))
  if (length(children) == 1L) return(sprintf("%s(%s)", .cgp_module_format(p, op), children[[1L]]))
  .prt_fold(p, op, children)
}

# version_tuple(a) < version_tuple(b)
.prt_version_lt <- function(a, b) utils::compareVersion(a, b) < 0L

# Derivative.variables: the variables, each repeated by its count
.prt_deriv_vars <- function(e) {
  out <- list()
  for (vc in e$args[-1L]) {
    if (is_sym(vc) && vc$type == "Tuple") {
      n <- as.integer(as.numeric(vc$args[[2L]]$value))
      out <- c(out, rep(list(vc$args[[1L]]), n))
    } else out[[length(out) + 1L]] <- vc
  }
  out
}

# the printed entries of a matrix as nested Python lists
.prt_matrix_data <- function(p, e)
  paste0("[", paste(vapply(.prt_matrix_tolist(e), function(r) paste0("[", .prt_join(p, r, ", "), "]"), ""), collapse = ", "), "]")

# MatMul in TensorFlow and PyTorch: the scalar factors times the folded
# matrix product
.prt_matmul <- function(op) function(p, e) {
  ism <- vapply(e$args, function(a) .is_matexpr(a) || (is_sym(a) && a$type == "ImmutableMatrix"), NA)
  mat_args <- e$args[ism]
  args <- e$args[!ism]
  if (length(args))
    return(sprintf("%s*%s", .cgp_paren(p, .mul(args), .CG_PRECEDENCE[["Mul"]]), .prt_fold(p, op, mat_args)))
  .prt_fold(p, op, mat_args)
}
.prt_matpow <- function(op) function(p, e) .prt_fold(p, op, rep(list(e$args[[1L]]), as.integer(as.numeric(e$args[[2L]]$value))))

# _print_Piecewise of both printers: nested where() calls
.prt_where <- function(p, e, fn) {
  first <- e$args[[1L]]
  if (length(e$args) == 1L)
    return(sprintf("%s(%s, %s, %s)", .cgp_module_format(p, fn), .cgpp(p, first$args[[2L]]), .cgpp(p, first$args[[1L]]), 0L))
  rest <- do.call(Piecewise, lapply(e$args[-1L], function(pr) list(pr$args[[1L]], pr$args[[2L]])))
  sprintf("%s(%s, %s, %s)", .cgp_module_format(p, fn), .cgpp(p, first$args[[2L]]), .cgpp(p, first$args[[1L]]), .cgpp(p, rest))
}

.prt_tf_fun <- .prt_print_mapped("TensorflowPrinter")
.cgp_def("TensorflowPrinter", c("ArrayPrinter", "AbstractPythonCodePrinter"),
         attrs = list(printmethod = "_tensorflowcode", module = "sympy.printing.tensorflow", mapping = .prt_tf_mapping,
                      `_default_settings` = c(.cgp_classes$AbstractPythonCodePrinter$attrs[["_default_settings"]],
                                              list(tensorflow_version = NULL)),
                      `_module` = "tensorflow", `_einsum` = "linalg.einsum", `_add` = "math.add",
                      `_transpose` = "transpose", `_ones` = "ones", `_zeros` = "zeros"),
         init = function(p, settings, nxt) {
           nxt(settings)
           p$tensorflow_version <- .cgp_s(p, "tensorflow_version")
           invisible(NULL)
         },
         methods = list(
           `_print_Function` = .prt_tf_fun, `_print_Expr` = .prt_tf_fun, `_print_Application` = .prt_tf_fun,
           `_print_MatrixExpr` = .prt_tf_fun, `_print_Relational` = .prt_tf_fun, `_print_Not` = .prt_tf_fun,
           `_print_And` = .prt_tf_fun, `_print_Or` = .prt_tf_fun, `_print_HadamardProduct` = .prt_tf_fun,
           `_print_Trace` = .prt_tf_fun, `_print_Determinant` = .prt_tf_fun,
           `_print_Inverse` = function(p, e) sprintf("%s(%s)", .cgp_module_format(p, "tensorflow.linalg.inv"), .cgpp(p, e$args[[1L]])),
           `_print_Transpose` = function(p, e) {
             v <- p$tensorflow_version
             op <- if (!is.null(v) && nzchar(v) && .prt_version_lt(v, "1.14")) "tensorflow.matrix_transpose" else "tensorflow.linalg.matrix_transpose"
             sprintf("%s(%s)", .cgp_module_format(p, op), .cgpp(p, e$args[[1L]]))
           },
           `_print_Derivative` = function(p, e) {
             vars <- .prt_deriv_vars(e)
             unfold <- function(ex, args) {
               if (!length(args)) return(.cgpp(p, ex))
               n <- length(args)
               sprintf("%s(%s, %s)[0]", .cgp_module_format(p, "tensorflow.gradients"), unfold(ex, args[-n]), .cgpp(p, args[[n]]))
             }
             unfold(e$args[[1L]], vars)
           },
           `_print_Piecewise` = function(p, e) {
             v <- p$tensorflow_version
             fn <- if (!is.null(v) && nzchar(v) && .prt_version_lt(v, "1.0")) "tensorflow.select" else "tensorflow.where"
             .prt_where(p, e, fn)
           },
           `_print_Pow` = function(p, e, ...) {
             if (identical(e$args[[2L]]$key, S_Half$key))
               return(sprintf("%s(%s)", .cgp_module_format(p, "tensorflow.math.sqrt"), .cgpp(p, e$args[[1L]])))
             sprintf("%s(%s, %s)", .cgp_module_format(p, "tensorflow.math.pow"), .cgpp(p, e$args[[1L]]), .cgpp(p, e$args[[2L]]))
           },
           `_print_MatrixBase` = function(p, e) {
             f <- if (length(free_symbols(e))) "tensorflow.Variable" else "tensorflow.constant"
             sprintf("%s(%s)", .cgp_module_format(p, f), .prt_matrix_data(p, e))
           },
           `_print_MatMul` = .prt_matmul("tensorflow.linalg.matmul"),
           `_print_MatPow` = .prt_matpow("tensorflow.linalg.matmul"),
           `_print_CodeBlock` = function(p, e) paste(vapply(e$args, function(a) .cgpp(p, a), ""), collapse = "\n"),
           `_print_isnan` = function(p, e) sprintf("tensorflow.math.is_nan(%s)", .cgpp(p, e$args[[1L]])),
           `_print_isinf` = function(p, e) sprintf("tensorflow.math.is_inf(%s)", .cgpp(p, e$args[[1L]]))
         ))

#' TensorFlow and PyTorch code printing
#'
#' Code printers after SymPy's `TensorflowPrinter`
#' (`sympy.printing.tensorflow`) and `TorchPrinter`
#' (`sympy.printing.pytorch`). `tensorflow_code()` and `torch_code()` print an
#' expression as Python code calling the TensorFlow or PyTorch functions:
#' elementary functions, relationals, logic, `Piecewise` (as `where()`),
#' derivatives (as gradient calls), matrices and matrix expressions, and array
#' expressions (as einsum and transpose calls).
#'
#' @param settings A named list of printer settings: `tensorflow_version`
#'   for `TensorflowPrinter()`; `requires_grad` and `dtype` for
#'   `TorchPrinter()`.
#' @param expr An expression.
#' @param ... Printer settings, such as `tensorflow_version = "1.13"`.
#' @param requires_grad Whether printed tensors require gradients.
#' @param dtype The dtype of printed tensors.
#' @param assign_to Optional variable to assign to.
#' @return A printer, or a character scalar of code.
#' @examples
#' x <- Symbol("x")
#' y <- Symbol("y")
#' tensorflow_code(sin(x) + x^2)
#' tensorflow_code(Piecewise(list(x, x < 0), list(y, TRUE)))
#' torch_code(Max(x, y))
#' torch_code(Matrix(list(list(x, 1), list(2, y))), requires_grad = TRUE)
#' printer_doprint(TensorflowPrinter(), exp(x))
#' @name tensorflow-torch-code
NULL

#' @rdname tensorflow-torch-code
#' @export
TensorflowPrinter <- function(settings = NULL) .cgp_new("TensorflowPrinter", settings)

#' @rdname tensorflow-torch-code
#' @export
tensorflow_code <- function(expr, ...) printer_doprint(TensorflowPrinter(list(...)), expr)
