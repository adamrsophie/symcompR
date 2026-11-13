# Algorithms as abstract syntax trees, after sympy/codegen/algorithms.py.

#' Newton's method as an abstract syntax tree
#'
#' `newtons_method(expr, wrt)` builds the [CodeBlock()] of a Newton
#' iteration solving `expr == 0` for `wrt`, as SymPy's `newtons_method()`.
#' `newtons_method_function(expr, wrt)` wraps it in a
#' [FunctionDefinition()] returning the root.
#'
#' @param expr The expression whose root is sought.
#' @param wrt The symbol to solve for.
#' @param atol The absolute tolerance.
#' @param delta `NULL` (a [Dummy()] inside a [Scope()]), or the symbol of
#'   the step.
#' @param rtol The relative tolerance.
#' @param debug Whether to print the iterates.
#' @param itermax `NULL`, or the largest number of iterations.
#' @param counter `NULL`, or the symbol counting the iterations.
#' @param delta_fn An R function of the expression and the variable giving
#'   the step; the default is the Newton step `-e/diff(e, x)`.
#' @param cse Whether to apply common subexpression elimination to the step.
#' @param handle_nan `NULL`, or a node executed when the step is not a
#'   number.
#' @param bounds `NULL`, or a list of the lower and upper bounds of `wrt`.
#' @return A [CodeBlock()] (in a [Scope()] when `delta` is `NULL`), or for
#'   `newtons_method_function()` a [FunctionDefinition()].
#' @examples
#' x <- Symbol("x", real = TRUE)
#' algo <- newtons_method(cos(x) - x^3, x)
#' cat(ccode(newtons_method_function(cos(x) - x^3, x)), "\n")
#' @export
newtons_method <- function(expr, wrt, atol = 1e-12, delta = NULL, rtol = 4e-16, debug = FALSE,
                           itermax = NULL, counter = NULL, delta_fn = function(e, x) -e / diff(e, x),
                           cse = FALSE, handle_nan = NULL, bounds = NULL) {
  expr <- as_sym(expr)
  if (is.null(delta)) {
    delta <- Dummy()
    wrapper <- Scope
    name_d <- "delta"
  } else {
    wrapper <- identity
    name_d <- delta$name
  }
  delta_expr <- delta_fn(expr, wrt)
  if (cse) {
    r <- cse(list(sym_factor(delta_expr)))
    whl_bdy <- lapply(r$replacements, function(pr) Assignment(pr[[1L]], pr[[2L]]))
    whl_bdy <- c(whl_bdy, list(Assignment(delta, r$reduced[[1L]])))
  } else {
    whl_bdy <- list(Assignment(delta, delta_expr))
  }
  if (!is.null(handle_nan)) whl_bdy <- c(whl_bdy, list(While(isnan(delta), CodeBlock(handle_nan, break_))))
  whl_bdy <- c(whl_bdy, list(AddAugmentedAssignment(wrt, delta)))
  if (!is.null(bounds)) whl_bdy <- c(whl_bdy, list(Assignment(wrt, Min(Max(wrt, bounds[[1L]]), bounds[[2L]]))))
  if (debug) {
    prnt <- Print(list(wrt, delta), paste0(wrt$name, "=%12.5g ", name_d, "=%12.5g\\n"))
    whl_bdy <- c(whl_bdy, list(prnt))
  }
  req <- Gt(Abs(delta), as_sym(atol) + as_sym(rtol) * Abs(wrt))
  declars <- list(Declaration(Variable(delta, type = real, value = S_Infinity)))
  if (!is.null(itermax)) {
    if (is.null(counter)) counter <- Dummy(integer = TRUE)
    v_counter <- variable_deduced(counter, 0L)
    declars <- c(declars, list(Declaration(v_counter)))
    whl_bdy <- c(whl_bdy, list(AddAugmentedAssignment(counter, 1L)))
    req <- And(req, Lt(counter, itermax))
  }
  whl <- While(req, do.call(CodeBlock, whl_bdy))
  blck <- declars
  if (debug) blck <- c(blck, list(Print(list(wrt), paste0(wrt$name, "=%12.5g\\n"))))
  blck <- c(blck, list(whl))
  wrapper(do.call(CodeBlock, blck))
}

.cg_symbol_of <- function(arg) {
  if (is_sym(arg) && arg$type == "Declaration") arg <- .cg_vals(arg)$variable
  if (is_sym(arg) && arg$type %in% c("Variable", "Pointer")) arg <- .cg_vals(arg)$symbol
  arg
}

#' @rdname newtons_method
#' @param params `NULL` (only `wrt`), or a list of the parameters (symbols,
#'   [Variable()] or [Pointer()] nodes).
#' @param func_name The name of the function.
#' @param attrs The attributes of the function.
#' @param ... Further arguments of `newtons_method()`.
#' @export
newtons_method_function <- function(expr, wrt, params = NULL, func_name = "newton", attrs = list(),
                                    delta = NULL, ...) {
  expr <- as_sym(expr)
  if (is.null(params)) params <- list(wrt)
  pointer_subs <- list()
  for (p in params) if (is_sym(p) && p$type == "Pointer") {
    s <- .cg_vals(p)$symbol
    pointer_subs[[length(pointer_subs) + 1L]] <- list(s, Symbol(sprintf("(*%s)", s$name)))
  }
  if (is.null(delta)) {
    delta <- Symbol(paste0("d_", wrt$name))
    if (has(expr, delta)) delta <- NULL
  }
  algo <- newtons_method(expr, wrt, delta = delta, ...)
  if (length(pointer_subs)) algo <- xreplace(algo, pointer_subs)
  if (algo$type == "Scope") algo <- .cg_vals(algo)$body
  known <- vapply(params, function(p) .k(.cg_symbol_of(p)), "")
  not_in_params <- Filter(function(s) !(.k(s) %in% known), free_symbols(expr))
  if (length(not_in_params))
    .cg_err("ValueError", paste0("Missing symbols in params: ",
                                 paste(vapply(not_in_params, .py_str_of, ""), collapse = ", ")))
  declars <- lapply(params, function(p) Variable(p, real))
  body <- CodeBlock(algo, Return(wrt))
  FunctionDefinition(real, func_name, declars, body, attrs = attrs)
}
