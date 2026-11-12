#' symcompR: Symbolic Computation in R
#'
#' A native R computer algebra system ported from the Python library SymPy,
#' using 'gmp' for exact integers and rationals. The SymPy Development Team
#' is gratefully acknowledged.
#'
#' Attaching the package masks no function of base R. SymPy names that clash
#' with base R carry a `sym_` prefix, such as [sym_factor()],
#' [sym_integrate()], [sym_symbols()], [sym_var()] and [sym_det()]; the
#' constants are `sym_pi`, `sym_E` and `sym_I`.
#'
#' @keywords internal
#' @importFrom gmp as.bigq as.bigz
"_PACKAGE"
