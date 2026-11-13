# Utilities shared by the differential equation solvers, after
# sympy/solvers/deutils.py: _preprocess, ode_order and _desolve.

# Session state of the ODE solvers.
.ode_state <- new.env()

# An error condition with SymPy's class name and its base classes.
.ode_cond <- function(cls, msg, parents = character()) {
  structure(list(message = paste(as.character(msg), collapse = ""), call = NULL),
            class = c(cls, parents, "sympy_error", "error", "condition"))
}
.ode_stop <- function(cls, msg, parents = character()) stop(.ode_cond(cls, msg, parents))
.ode_nie <- function(msg) .ode_stop("NotImplementedError", msg)

# Python truthiness of an expression: only a zero number is false.
.ode_true <- function(e) {
  if (is.null(e)) return(FALSE)
  if (is.logical(e)) return(isTRUE(e))
  if (is.list(e) && !is_sym(e)) return(length(e) > 0L)
  if (!is_sym(e)) return(as.logical(e))
  if (e$type == "BooleanAtom") return(e$name == "true")
  .sp_bool(e)
}

# f(x).diff(x, n), with f(x) itself for n = 0
.ode_d <- function(fx, x, n = 1L) {
  n <- as.integer(n)
  if (n == 0L) return(fx)
  .diff_spec_eval(fx, list(list(x, n)))
}

# The undefined function applications in an expression.
.ode_undef_atoms <- function(e) {
  out <- Filter(.is_undefined_application, atoms(e, "Function"))
  if (length(out)) out[!duplicated(vapply(out, .k, ""))] else out
}

# Derivative.variables, with each variable repeated by its count
.ode_deriv_vars <- function(d) {
  out <- list()
  for (t in d$args[-1L]) {
    n <- as.integer(as.numeric(t$args[[2L]]$value))
    for (i in seq_len(n)) out[[length(out) + 1L]] <- t$args[[1L]]
  }
  out
}

# Derivative.variable_count as a list of list(variable, count)
.ode_variable_count <- function(d) lapply(d$args[-1L], function(t) list(t$args[[1L]], as.integer(as.numeric(t$args[[2L]]$value))))

# sympy.solvers.deutils._preprocess
.ode_preprocess <- function(expr, func = NULL, hint = "_Integral") {
  if (expr$type == "Pow") {
    # if f(x)**p=0 then f(x)=0 (p>0)
    if (isTRUE(ask_fact(expr$args[[2L]], "positive"))) expr <- expr$args[[1L]]
  }
  derivs <- atoms(expr, "Derivative")
  if (is.null(func)) {
    funcs <- list()
    for (d in derivs) funcs <- c(funcs, .ode_undef_atoms(d))
    if (length(funcs)) funcs <- funcs[!duplicated(vapply(funcs, .k, ""))]
    if (length(funcs) != 1L)
      .ode_stop("ValueError", paste0("The function cannot be automatically detected for ", .str(expr), "."))
    func <- funcs[[1L]]
  }
  fvars <- vapply(func$args, .k, "")
  if (is.null(hint)) return(list(expr, func))
  reps <- list()
  for (d in derivs) {
    dv <- vapply(.ode_deriv_vars(d), .k, "")
    if (!endsWith(hint, "_Integral") || .has1(d, func) || any(dv %in% fvars))
      reps[[length(reps) + 1L]] <- list(d, doit(d))
  }
  eq <- if (length(reps)) subs(expr, reps) else expr
  list(eq, func)
}

#' Order of a differential equation
#'
#' `ode_order(expr, func)` returns the order of the differential equation
#' `expr` in the function `func`, as SymPy's `ode_order()`: a derivative of
#' an expression adds its own order to the highest order among the
#' arguments of that expression.
#'
#' @param expr An expression or equation.
#' @param func The function, such as `f(x)`.
#' @return An integer.
#' @examples
#' x <- Symbol("x"); f <- Function("f")
#' ode_order(Derivative(f(x), x, 3) + f(x), f(x))
#' @export
ode_order <- function(expr, func) {
  expr <- as_sym(expr)
  func <- as_sym(func)
  .ode_order(expr, func)
}

.ode_order <- function(expr, func) {
  # a = Wild('a', exclude=[func]); expr.match(a)
  if (!.has1(expr, func)) return(0L)
  if (expr$type == "Derivative") {
    rv <- length(.ode_deriv_vars(expr))
    if (identical(expr$args[[1L]]$key, func$key)) return(rv)
    args <- expr$args[[1L]]$args
    if (length(args)) rv <- rv + max(vapply(args, .ode_order, 0L, func = func))
    return(rv)
  }
  if (!length(expr$args)) return(0L)
  max(vapply(expr$args, .ode_order, 0L, func = func))
}

# sympy.solvers.deutils._desolve for ordinary differential equations: the
# dictionary of hints that dsolve() goes on to solve.
.ode_desolve <- function(eq, func = NULL, hint = "default", ics = NULL, simplify = TRUE, prep = TRUE, kw = list()) {
  if (eq$type == "Equality") eq <- .sub(eq$args[[1L]], eq$args[[2L]])
  if (prep || is.null(func)) {
    pp <- .ode_preprocess(eq, func)
    eq <- pp[[1L]]
    func <- pp[[2L]]
    prep <- FALSE
  }
  xi <- kw$xi
  eta <- kw$eta
  x0 <- if (is.null(kw$x0)) S_Zero else kw$x0
  terms <- kw$n
  if (!isFALSE(kw$classify)) {
    hints <- .ode_classify(eq, func, dict = TRUE, ics = ics, prep = prep, xi = xi, eta = eta, n = terms,
                           kw = list(x0 = x0, hint = hint))
  } else {
    hints <- kw$hints
    if (is.null(hints)) {
      hints <- list(default = hint)
      hints[[hint]] <- kw$match
      hints$order <- kw$order
    }
  }
  if (is.null(hints$default)) {
    if (!(hint %in% .ODE_ALLHINTS) && hint != "default")
      .ode_stop("ValueError", paste0("Hint not recognized: ", hint))
    else if (!(hint %in% hints$ordered_hints) && hint != "default")
      .ode_stop("ValueError", paste0("ODE ", .str(eq), " does not match hint ", hint))
    else if (identical(hints$order, 0L))
      .ode_stop("ValueError", paste0(.str(eq), " is not a solvable differential equation in ", .str(func)))
    else .ode_nie(paste0("solve: Cannot solve ", .str(eq)))
  }
  if (hint == "default") {
    return(.ode_desolve(eq, func, ics = ics, hint = hints$default, simplify = simplify, prep = prep,
                        kw = list(x0 = x0, classify = FALSE, order = hints$order, match = hints[[hints$default]],
                                  xi = xi, eta = eta, n = terms)))
  } else if (hint %in% c("all", "all_Integral", "best")) {
    retdict <- list()
    gethints <- setdiff(names(hints), c("order", "default", "ordered_hints"))
    if (hint == "all_Integral") {
      for (i in names(hints))
        if (endsWith(i, "_Integral")) gethints <- setdiff(gethints, sub("_Integral$", "", i))
      gethints <- setdiff(gethints, c("1st_homogeneous_coeff_best", "1st_power_series", "lie_group",
                                      "2nd_power_series_ordinary", "2nd_power_series_regular"))
    }
    for (i in gethints) {
      retdict[[i]] <- .ode_desolve(eq, func, ics = ics, hint = i, simplify = simplify, prep = prep,
                                   kw = list(x0 = x0, classify = FALSE, n = terms, order = hints$order,
                                             match = hints[[i]]))
    }
    retdict$all <- TRUE
    retdict$eq <- eq
    return(retdict)
  } else if (!(hint %in% .ODE_ALLHINTS)) {
    .ode_stop("ValueError", paste0("Hint not recognized: ", hint))
  } else if (is.null(hints[[hint]]) && !(hint %in% names(hints))) {
    .ode_stop("ValueError", paste0("ODE ", .str(eq), " does not match hint ", hint))
  } else {
    hints$hint <- hint
  }
  hints$func <- func
  hints$eq <- eq
  hints
}
