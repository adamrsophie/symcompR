# implemented_function() and lambdastr(), after sympy/utilities/lambdify.py.

.cg_imp_key <- function(imp) if (is_sym(imp)) imp$key else paste(deparse(imp), collapse = " ")

# _imp_namespace: two implementations under one function name are an error
.cg_imp_namespace <- function(e) {
  seen <- list()
  walk <- function(z) {
    if (is.list(z) && !is_sym(z)) { for (a in z) walk(a); return(invisible()) }
    if (!is_sym(z)) return(invisible())
    if (z$type == "Function" && !is.null(z$extra$imp)) {
      ik <- .cg_imp_key(z$extra$imp)
      if (!is.null(seen[[z$name]]) && !identical(seen[[z$name]], ik))
        .cg_err("ValueError", sprintf("We found more than one implementation with name \"%s\"", z$name))
      seen[[z$name]] <<- ik
    }
    for (a in z$args) walk(a)
  }
  walk(e)
  invisible()
}

# f(*args) for an undefined function `name` that carries the implementation
# `imp` (SymPy's _imp_ class attribute). The implementation is part of the
# key, as two such functions with different implementations differ.
.cg_imp_call <- function(name, imp, args) {
  r <- fcall(name, lapply(args, as_sym))
  if (!(r$type == "Function" && identical(r$name, name))) return(r)
  ik <- .cg_imp_key(imp)
  new_sym("Function", r$args, name = name, extra = list(imp = imp), key = paste0(r$key, "{imp:", ik, "}"))
}

#' Add a numerical implementation to an undefined function
#'
#' `implemented_function(symfunc, implementation)` returns an undefined
#' function named after `symfunc` whose applications carry `implementation`,
#' as SymPy's `implemented_function()`. [lambdify()] evaluates the
#' applications with the implementation, and the code printers inline an
#' implementation that is a [Lambda()].
#'
#' @param symfunc A function name, or an undefined function made by
#'   [Function()].
#' @param implementation A [Lambda()] or an R function.
#' @return An undefined function (class `sym_function`).
#' @examples
#' x <- Symbol("x")
#' f <- implemented_function("f", Lambda(x, 2 * x))
#' ccode(f(x))
#' @export
implemented_function <- function(symfunc, implementation) {
  if (is.character(symfunc) && length(symfunc) == 1L) name <- symfunc
  else if (inherits(symfunc, "sym_function") && !is.null(attr(symfunc, "sym_name"))) name <- attr(symfunc, "sym_name")
  else .cg_err("ValueError", "\nsymfunc should be either a string or an UndefinedFunction instance.")
  force(implementation)
  fn <- function(...) .cg_imp_call(name, implementation, list(...))
  structure(fn, class = c("sym_function", "function"), sym_name = name, `_imp_` = implementation)
}

.cg_iterable <- function(a) is.list(a) && !is_sym(a) || (is_sym(a) && a$type == "Tuple")
.cg_items <- function(a) if (is_sym(a)) a$args else a

# _recursive_to_string
.cg_rec_to_string <- function(doprint, arg) {
  if (is_sym(arg) || inherits(arg, "symmatrix")) return(doprint(arg))
  if (is.character(arg)) return(arg)
  if (is.list(arg)) return(paste0("[", paste(vapply(arg, function(a) .cg_rec_to_string(doprint, a), ""), collapse = ", "), "]"))
  doprint(arg)
}

#' Python source of a lambda function
#'
#' `lambdastr(args, expr)` returns the source of a Python lambda function
#' that evaluates `expr` at `args`, as SymPy's `lambdastr()`.
#'
#' @param args A symbol, or a list of symbols (possibly nested).
#' @param expr An expression, or a list of expressions.
#' @param printer A function of an expression giving its source; the default
#'   is SymPy's `LambdaPrinter`.
#' @param dummify Whether to replace the arguments by dummy symbols. The
#'   default dummifies when an argument contains a function or a derivative.
#' @return A character scalar.
#' @examples
#' x <- Symbol("x"); y <- Symbol("y")
#' lambdastr(x, x^2)
#' lambdastr(list(x, list(y)), x + y)
#' @export
lambdastr <- function(args, expr, printer = NULL, dummify = NULL) {
  lambdarepr <- if (is.null(printer)) .cg_lambdarepr
                else if (inherits(printer, "cg_printer")) function(e) .cgp_doprint(printer, e, NULL)
                else printer
  dstr <- function(d) sstr(d)
  sub_args <- function(a, dd) {
    if (is.character(a)) return(list(s = a, dd = dd))
    if (.cg_iterable(a)) {
      parts <- character()
      for (z in .cg_items(a)) { r <- sub_args(z, dd); dd <- r$dd; parts <- c(parts, r$s) }
      return(list(s = paste(parts, collapse = ","), dd = dd))
    }
    if (is_sym(a) && (a$type %in% c("Symbol", "Dummy", "Function", "Derivative") || .cgp_is_function(a))) {
      d <- Dummy()
      dd[[length(dd) + 1L]] <- list(a, d)
      return(list(s = dstr(d), dd = dd))
    }
    list(s = .py_str_of(a), dd = dd)
  }
  sub_expr <- function(e, dd) {
    if (is.list(e) && !is_sym(e)) return(lapply(e, sub_expr, dd = dd))
    if (is.character(e)) return(e)
    xreplace(as_sym(e), dd)
  }
  isiter <- function(l) .cg_iterable(l)
  flat <- function(a) if (isiter(a)) do.call(c, lapply(.cg_items(a), flat)) else list(a)
  flat_indexes <- function(it) {
    out <- list(); n <- 0L
    for (el in .cg_items(it)) {
      if (isiter(el)) for (nd in flat_indexes(el)) out[[length(out) + 1L]] <- c(n, nd)
      else out[[length(out) + 1L]] <- n
      n <- n + 1L
    }
    out
  }
  if (is.null(dummify)) {
    al <- if (isiter(args)) .cg_items(args) else list(args)
    dummify <- any(vapply(al, function(a) is_sym(a) && length(atoms(a, c("Function", "Derivative"))) > 0L, TRUE))
  }
  if (isiter(args) && any(vapply(.cg_items(args), isiter, TRUE))) {
    dum_args <- vapply(seq_along(.cg_items(args)) - 1L, function(i) dstr(Dummy(as.character(i))), "")
    indexed <- paste(vapply(flat_indexes(args), function(ind)
      paste0(dum_args[ind[1L] + 1L], paste(sprintf("[%d]", ind[-1L]), collapse = "")), ""), collapse = ",")
    lstr <- lambdastr(flat(args), expr, printer = printer, dummify = dummify)
    return(sprintf("lambda %s: (%s)(%s)", paste(dum_args, collapse = ","), lstr, indexed))
  }
  dd <- list()
  if (dummify) {
    r <- sub_args(args, dd); args <- r$s; dd <- r$dd
  } else if (!is.character(args) && isiter(args)) {
    args <- paste(vapply(.cg_items(args), .py_str_of, ""), collapse = ",")
  } else if (!is.character(args)) args <- .py_str_of(args)
  if (dummify && !is.character(expr)) expr <- sub_expr(expr, dd)
  sprintf("lambda %s: (%s)", args, .cg_rec_to_string(lambdarepr, expr))
}
.cg_lambdarepr <- function(e) .cgp_doprint(.cgp_new("LambdaPrinter"), e, NULL)
