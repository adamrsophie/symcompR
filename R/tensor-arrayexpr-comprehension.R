# Array comprehensions, after sympy/tensor/array/array_comprehension.py.

.ac_items <- function(l) if (is_sym(l)) l$args else as.list(l)

# ArrayComprehension._check_limits_validity
.ac_check_limits <- function(lims) {
  lapply(lims, function(l) {
    it <- .ac_items(l)
    var <- as_sym(it[[1L]]); inf <- as_sym(it[[2L]])
    sup <- if (is.list(it[[3L]]) && !is_sym(it[[3L]])) do.call(Tuple, lapply(it[[3L]], as_sym)) else as_sym(it[[3L]])
    ok_bound <- function(b) b$type != "Tuple" && all(vapply(.ac_atoms(b), function(a) a$type %in% c("Symbol", "Integer"), NA))
    if (!ok_bound(inf) || !ok_bound(sup))
      .perr("TypeError", "Bounds should be an Expression(combination of Integer and Symbol)")
    if (isTRUE(.facts(.sub(inf, sup))[["positive"]]))
      .perr("ValueError", "Lower bound should be inferior to upper bound")
    if (.has1(inf, var) || .has1(sup, var)) .perr("ValueError", "Variable should not be part of its bounds")
    Tuple(var, inf, sup)
  })
}
.ac_atoms <- function(e) if (!length(e$args)) list(e) else do.call(c, lapply(e$args, .ac_atoms))

.ac_new <- function(type, fun, lims, extra = list()) {
  if (any(vapply(lims, function(l) length(.ac_items(l)) != 3L, NA)))
    .perr("ValueError", "ArrayComprehension requires values lower and upper bound for the expression")
  lims <- .ac_check_limits(lims)
  shape <- lapply(lims, function(l) .add(list(.sub(l$args[[3L]], l$args[[2L]]), S_One)))
  new_sym(type, c(if (!is.null(fun)) list(fun), lims), extra = c(list(shape = shape), extra))
}

#' Array comprehensions
#'
#' `ArrayComprehension(fun, ...)` is the array of the values of `fun` when
#' the variables of the limits `list(var, lower, upper)` run over their
#' ranges; `ArrayComprehensionMap(f, ...)` takes an R function of one
#' argument applied to the product of the variable values. [doit()]
#' expands an array with numeric bounds into an explicit array.
#'
#' @param fun An expression.
#' @param f An R function of zero or one argument.
#' @param ... Limits, each `list(var, lower, upper)` or a [Tuple()].
#' @return A `sym` array comprehension.
#' @examples
#' i <- Symbol("i"); j <- Symbol("j")
#' a <- ArrayComprehension(i * j, list(i, 1, 3), list(j, 2, 4))
#' doit(a)
#' @export
ArrayComprehension <- function(fun, ...) .ac_new("ArrayComprehension", as_sym(fun), list(...))

#' @rdname ArrayComprehension
#' @export
ArrayComprehensionMap <- function(f, ...) {
  if (!is.function(f)) .perr("ValueError", "Data type not supported")
  e <- .ac_new("ArrayComprehensionMap", NULL, list(...), list(fun = f))
  e$key <- paste0(e$key, "|fn:", paste(deparse(f), collapse = ""))
  e
}

.ac_limits <- function(e) if (e$type == "ArrayComprehension") e$args[-1L] else e$args
register_bound_vars("ArrayComprehension", function(e) vapply(e$args[-1L], function(l) l$args[[1L]]$key, ""))
register_rebuilder("ArrayComprehension", function(e, args) do.call(ArrayComprehension, args))
register_rebuilder("ArrayComprehensionMap", function(e, args) do.call(ArrayComprehensionMap, c(list(e$extra$fun), args)))
register_str_printer("ArrayComprehensionMap", function(e)
  paste0("ArrayComprehensionMap(", paste(vapply(e$args, .str, ""), collapse = ", "), ")"))

.ac_is_numeric_shape <- function(e)
  all(vapply(.ac_limits(e), function(l) !length(.free_sym_list(l$args[[2L]])) && !length(.free_sym_list(l$args[[3L]])), NA))

.ac_element <- function(e, vals) {
  if (e$type == "ArrayComprehensionMap") {
    f <- e$extra$fun
    if (!length(formals(f))) return(as_sym(f()))
    return(as_sym(f(.mul(lapply(vals, .int)))))
  }
  vars <- lapply(.ac_limits(e), function(l) l$args[[1L]])
  t <- e$args[[1L]]
  for (k in seq_along(vars)) t <- subs(t, vars[[k]], .int(vals[k]))
  t
}

.ac_expand <- function(e) {
  lims <- .ac_limits(e)
  rng <- lapply(lims, function(l) seq(.int_val(l$args[[2L]]), .int_val(l$args[[3L]])))
  grid <- as.matrix(expand.grid(rev(rng)))[, rev(seq_along(rng)), drop = FALSE]
  data <- lapply(seq_len(nrow(grid)), function(r) .ac_element(e, grid[r, ]))
  Array(data, vapply(e$extra$shape, .int_val, 0))
}

register_doit("ArrayComprehension", function(e) if (.ac_is_numeric_shape(e)) .ac_expand(e) else e)
register_doit("ArrayComprehensionMap", function(e) if (.ac_is_numeric_shape(e)) .ac_expand(e) else e)
