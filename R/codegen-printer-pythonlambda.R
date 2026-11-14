# The printers of sympy/printing/lambdarepr.py that the base port lacks:
# NumExprPrinter, IntervalPrinter and lambdarepr(). LambdaPrinter itself is
# in R/codegen-printer-python.R.

.prt_numexpr_functions <- c(sin = "sin", cos = "cos", tan = "tan", asin = "arcsin", acos = "arccos", atan = "arctan",
                            atan2 = "arctan2", sinh = "sinh", cosh = "cosh", tanh = "tanh", asinh = "arcsinh",
                            acosh = "arccosh", atanh = "arctanh", ln = "log", log = "log", exp = "exp", sqrt = "sqrt",
                            Abs = "abs", conjugate = "conj", im = "imag", re = "real", where = "where",
                            complex = "complex", contains = "contains")

.prt_is_codegen_ast <- function(e) is_sym(e) && "CodegenAST" %in% .cgp_mro(e, full = TRUE)
.prt_is_nee <- function(e) is_sym(e) && e$type == "NumExprEvaluate"

# NumExprPrinter._print_Function
.prt_numexpr_function <- function(p, e) {
  fname <- .cgp_fname(e)
  nstr <- .prt_numexpr_functions[fname]
  if (is.na(nstr)) .cg_err("TypeError", sprintf("numexpr does not support function '%s'", fname))
  sprintf("%s(%s)", nstr, .cgp_call(p, "_print_seq", e$args))
}
.prt_blacklisted <- function(p, e) .cg_err("TypeError", sprintf("numexpr cannot be used with %s", .cgp_clsname(e)))

.cgp_def("NumExprPrinter", "LambdaPrinter",
         attrs = list(printmethod = "_numexprcode", `_numexpr_module` = "numexpr"),
         methods = c(list(
           `_print_ImaginaryUnit` = function(p, e) "1j",
           `_print_seq` = function(p, seq, delimiter = ", ") {
             s <- vapply(seq, function(i) .cgpp(p, i), "")
             if (length(s)) paste(s, collapse = delimiter) else ""
           },
           `_print_Function` = .prt_numexpr_function,
           `_print_Piecewise` = function(p, e) {
             exprs <- vapply(e$args, function(a) .cgpp(p, a$args[[1L]]), "")
             conds <- vapply(e$args, function(a) .cgpp(p, a$args[[2L]]), "")
             ans <- character(); pc <- 0L; last_true <- FALSE
             for (k in seq_along(conds)) {
               if (conds[k] == "True") { ans <- c(ans, exprs[k]); last_true <- TRUE; break }
               ans <- c(ans, sprintf("where(%s, %s, ", conds[k], exprs[k]))
               pc <- pc + 1L
             }
             if (!last_true) ans <- c(ans, "log(-1)")
             paste0(paste(ans, collapse = ""), strrep(")", pc))
           },
           `_print_ITE` = function(p, e) .cgpp(p, .cg_ite_piecewise(e)),
           `_print_NumExprEvaluate` = function(p, e)
             sprintf("%s('%s', truediv=True)", .cgp_module_format(p, paste0(.cgp_a(p, "_numexpr_module"), ".evaluate")),
                     .cgpp(p, .cg_vals(e)$expr)),
           doprint = function(p, e, assign_to = NULL) {
             if (!.prt_is_codegen_ast(e)) e <- NumExprEvaluate(e)
             .cgp_doprint(p, e, assign_to)
           },
           `_print_Return` = function(p, e) {
             r <- e$args[[1L]]
             if (!.prt_is_nee(r)) e <- sym_rebuild(e, list(NumExprEvaluate(r)))
             .cgp_super(p, "NumExprPrinter", "_print_Return", e)
           },
           `_print_Assignment` = function(p, e) {
             if (!.prt_is_nee(e$args[[2L]])) e <- Assignment(e$args[[1L]], NumExprEvaluate(e$args[[2L]]))
             .cgp_super(p, "NumExprPrinter", "_print_Assignment", e)
           },
           # super()._print_CodeBlock(self, ...) passes the printer twice in SymPy
           `_print_CodeBlock` = function(p, e)
             .cg_err("TypeError", "CodePrinter._print_CodeBlock() takes 2 positional arguments but 3 were given")
         ),
         stats::setNames(rep(list(.prt_blacklisted), 13L),
                  paste0("_print_", c("SparseRepMatrix", "MutableSparseMatrix", "ImmutableSparseMatrix", "Matrix",
                                      "DenseMatrix", "MutableDenseMatrix", "ImmutableMatrix", "ImmutableDenseMatrix",
                                      "list", "tuple", "Tuple", "dict", "Dict"))),
         stats::setNames(rep(list(.prt_numexpr_function), length(.prt_numexpr_functions)),
                  paste0("_print_", names(.prt_numexpr_functions)))))

# IntervalPrinter(MpmathPrinter, LambdaPrinter): numbers as mpi intervals
.cgp_def("IntervalPrinter", c("MpmathPrinter", "LambdaPrinter"),
         methods = list(
           `_print_Integer` = function(p, e) sprintf("mpi('%s')", .cgp_super(p, "PythonCodePrinter", "_print_Integer", e)),
           `_print_Rational` = function(p, e) sprintf("mpi('%s')", .cgp_super(p, "PythonCodePrinter", "_print_Rational", e)),
           `_print_Half` = function(p, e) sprintf("mpi('%s')", .cgp_super(p, "PythonCodePrinter", "_print_Rational", e)),
           `_print_Pow` = function(p, e, ...) .cgp_super(p, "MpmathPrinter", "_print_Pow", e, rational = TRUE)
         ))

#' Lambda code printers
#'
#' Printers after `sympy.printing.lambdarepr`. `lambdarepr()` prints an
#' expression as the code `lambdify()` evaluates, with `LambdaPrinter()`.
#' `NumExprPrinter()` prints code for the numexpr library: the expression is
#' wrapped in `numexpr.evaluate('...', truediv=True)`, and functions numexpr
#' lacks are an error. `IntervalPrinter()` prints numbers as mpmath `mpi`
#' intervals.
#'
#' @param expr An expression.
#' @param ... Printer settings.
#' @param settings A named list of printer settings.
#' @return A character scalar of code, or a printer.
#' @examples
#' x <- Symbol("x")
#' y <- Symbol("y")
#' lambdarepr(Piecewise(list(1, x < 1), list(2, TRUE)))
#' printer_doprint(NumExprPrinter(), sin(x) + cos(y))
#' printer_doprint(IntervalPrinter(), x/2 + 1)
#' @name lambda-printers
NULL

#' @rdname lambda-printers
#' @export
lambdarepr <- function(expr, ...) printer_doprint(LambdaPrinter(list(...)), expr)
#' @rdname lambda-printers
#' @export
NumExprPrinter <- function(settings = NULL) .cgp_new("NumExprPrinter", settings)
#' @rdname lambda-printers
#' @export
IntervalPrinter <- function(settings = NULL) .cgp_new("IntervalPrinter", settings)
