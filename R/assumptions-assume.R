# Predicates and applied predicates, after sympy/assumptions/assume.py and
# the AssumptionKeys class of sympy/assumptions/ask.py.
#
# A predicate is a sym node of type "Predicate" (a predicate of Q, carrying
# its SymPy class name for sorting) or "UndefinedPredicate" (made by
# Predicate("name")). Applying one gives an "AppliedPredicate" node, or an
# "AppliedBinaryRelation" for Q$eq, Q$ne, Q$gt, Q$ge, Q$lt and Q$le, whose
# first argument is the predicate. These are boolean expressions: they go
# inside And(), Or() and Not() and print as in SymPy.

.PRED_TYPES <- c("Predicate", "UndefinedPredicate")
.APPLIED_TYPES <- c("AppliedPredicate", "AppliedBinaryRelation")
.BINREL_NAMES <- c("eq", "ne", "gt", "ge", "lt", "le")

.is_predicate <- function(e) is_sym(e) && e$type %in% .PRED_TYPES
.is_applied <- function(e) is_sym(e) && e$type %in% .APPLIED_TYPES

# The predicate nodes of Q, made on first use.
.pred_cache <- new.env(parent = emptyenv())
.builtin_pred <- function(name) {
  p <- .pred_cache[[name]]
  if (is.null(p)) {
    p <- new_sym("Predicate", name = name, extra = list(sympy_class = .ASK_PRED_CLASS[[name]]))
    assign(name, p, envir = .pred_cache)
  }
  p
}
.qp <- function(name) .builtin_pred(name)

# AppliedPredicate.function and .arguments
.ap_fun <- function(e) e$args[[1L]]
.ap_args <- function(e) e$args[-1L]
.ap_name <- function(e) e$args[[1L]]$name

# _sympify for predicate arguments: an explicit matrix becomes immutable
.ask_sympify <- function(a) {
  if (inherits(a, "symmatrix")) return(.as_block(a))
  as_sym(a)
}

.ask_type_error <- function(msg) stop(.sympy_condition("TypeError", msg))
.ask_value_error <- function(msg) stop(.sympy_condition("ValueError", msg))

.as_pred <- function(p) {
  if (inherits(p, "sym_predicate_fn")) return(attr(p, "predicate"))
  if (.is_predicate(p)) return(p)
  .ask_type_error(paste0(if (is_sym(p)) .str(p) else class(p)[1L], " is not a Predicate."))
}

# Predicate.__call__, with the overrides of IsTruePredicate, BinaryRelation
# and UndefinedPredicate
.apply_pred <- function(pred, args) {
  args <- lapply(args, .ask_sympify)
  if (pred$type == "Predicate") {
    nm <- pred$name
    if (nm == "is_true") {
      if (length(args) != 1L)
        .ask_type_error(sprintf("IsTruePredicate.__call__() takes 2 positional arguments but %d were given", length(args) + 1L))
      a <- args[[1L]]
      if (.is_applied(a)) return(a)
      if (is_relational(a)) return(.rel_to_binrel(a))
    } else if (nm %in% .BINREL_NAMES) {
      if (length(args) != 2L)
        .ask_type_error(sprintf("Q.%s takes two arguments, but got %d.", nm, length(args)))
      return(new_sym("AppliedBinaryRelation", c(list(pred), args)))
    }
  } else if (length(args) != 1L) {
    .ask_type_error(sprintf("UndefinedPredicate.__call__() takes 2 positional arguments but %d were given", length(args) + 1L))
  }
  new_sym("AppliedPredicate", c(list(pred), args))
}

.REL_TO_BINREL <- c(Equality = "eq", Unequality = "ne", StrictGreaterThan = "gt",
                    StrictLessThan = "lt", GreaterThan = "ge", LessThan = "le")
.rel_to_binrel <- function(e) .apply_pred(.qp(.REL_TO_BINREL[[e$type]]), e$args)

# A callable R function standing for a predicate: Q$positive(x) applies it,
# and as_sym(Q$positive) gives the predicate node.
.pred_fn <- function(pred) {
  f <- function(...) .apply_pred(pred, list(...))
  attr(f, "predicate") <- pred
  class(f) <- c("sym_predicate_fn", "function")
  f
}

#' @export
print.sym_predicate_fn <- function(x, ...) {
  cat(.str(attr(x, "predicate")), "\n", sep = "")
  invisible(x)
}

#' @export
format.sym_predicate_fn <- function(x, ...) .str(attr(x, "predicate"))

# Q$even | Q$odd and !Q$even combine the predicate nodes
#' @export
Ops.sym_predicate_fn <- function(e1, e2) {
  f <- get(.Generic, envir = baseenv())
  if (missing(e2)) f(as_sym(e1)) else f(as_sym(e1), as_sym(e2))
}

.applied_str <- function(e) paste0(.str(e$args[[1L]]), "(", paste(vapply(e$args[-1L], .str, ""), collapse = ", "), ")")
.pred_latex <- function(e) paste0("\\operatorname{Q}_{\\text{", gsub("_", "\\_", e$name, fixed = TRUE), "}}")
.applied_latex <- function(e) paste0(.pred_latex(e$args[[1L]]), "(", paste(vapply(e$args[-1L], .latex, ""), collapse = ", "), ")")

register_str_printer("Predicate", function(e) paste0("Q.", e$name))
register_str_printer("UndefinedPredicate", function(e) paste0("Q.", e$name))
register_str_printer("AppliedPredicate", .applied_str)
register_str_printer("AppliedBinaryRelation", .applied_str)
register_latex_printer("Predicate", .pred_latex)
register_latex_printer("UndefinedPredicate", .pred_latex)
register_latex_printer("AppliedPredicate", .applied_latex)
register_latex_printer("AppliedBinaryRelation", .applied_latex)

#' Assumption predicates
#'
#' `Q` holds the predicates of SymPy's `Q` object, such as `Q$positive`,
#' `Q$integer`, `Q$invertible` or `Q$eq`. Calling one gives an applied
#' predicate, a boolean expression that prints as in SymPy (`Q.positive(x)`)
#' and can be combined with `&`, `|` and `!` and passed to [ask()],
#' [refine()] and [assuming()]. Asking for a predicate that `Q` does not
#' have raises an `AttributeError`; `Q$name <- p` adds a predicate made with
#' [Predicate()], and `Q$name <- NULL` removes it again.
#'
#' The predicates are `hermitian`, `antihermitian`, `real`, `extended_real`,
#' `imaginary`, `complex`, `algebraic`, `transcendental`, `integer`,
#' `noninteger`, `rational`, `irrational`, `finite`, `infinite`,
#' `positive_infinite`, `negative_infinite`, `positive`, `negative`, `zero`,
#' `extended_positive`, `extended_negative`, `nonzero`, `nonpositive`,
#' `nonnegative`, `extended_nonzero`, `extended_nonpositive`,
#' `extended_nonnegative`, `even`, `odd`, `prime`, `composite`,
#' `commutative`, `is_true`, the matrix predicates `symmetric`,
#' `invertible`, `orthogonal`, `unitary`, `positive_definite`,
#' `upper_triangular`, `lower_triangular`, `diagonal`, `fullrank`, `square`,
#' `integer_elements`, `real_elements`, `complex_elements`, `singular`,
#' `normal`, `triangular`, `unit_triangular`, and the binary relations `eq`,
#' `ne`, `gt`, `ge`, `lt` and `le`.
#'
#' @format An object of class `AssumptionKeys`.
#' @examples
#' x <- Symbol("x")
#' Q$positive(x)
#' Q$positive(x) & !Q$integer(x)
#' Q$eq(x, 1)
#' srepr(Q$positive(x))
#' ask(Q$positive(x^2 + 1), Q$real(x))
#' @export
Q <- structure(new.env(parent = emptyenv()), class = "AssumptionKeys")

#' @export
`$.AssumptionKeys` <- function(x, name) {
  f <- get0(name, envir = x, inherits = FALSE)
  if (!is.null(f)) return(f)
  if (name %in% names(.ASK_PRED_CLASS)) {
    f <- .pred_fn(.builtin_pred(name))
    assign(name, f, envir = x)
    return(f)
  }
  stop(.sympy_condition("AttributeError", sprintf("'AssumptionKeys' object has no attribute '%s'", name)))
}

#' @export
`[[.AssumptionKeys` <- function(x, i, ...) `$.AssumptionKeys`(x, i)

#' @export
`$<-.AssumptionKeys` <- function(x, name, value) {
  if (is.null(value)) {
    if (exists(name, envir = x, inherits = FALSE)) rm(list = name, envir = x)
  } else {
    if (!inherits(value, "sym_predicate_fn")) value <- .pred_fn(.as_pred(value))
    assign(name, value, envir = x)
  }
  x
}

#' @export
names.AssumptionKeys <- function(x) unique(c(names(.ASK_PRED_CLASS), ls(x, sorted = TRUE)))

#' @export
print.AssumptionKeys <- function(x, ...) {
  cat("<AssumptionKeys: ", paste(names(x), collapse = ", "), ">\n", sep = "")
  invisible(x)
}

# Handlers of predicates made with Predicate(name, dispatch = TRUE), and
# handlers added to the predicates of Q, by predicate key.
.ask_custom <- new.env(parent = emptyenv())

#' Define an assumption predicate
#'
#' `Predicate(name)` makes a new predicate, as SymPy's `Predicate("name")`
#' does: it prints as `Q.name`, can be applied to one expression, and
#' [ask()] answers it only from the assumptions it is given. With
#' `dispatch = TRUE` the predicate also has its own handler, so
#' [register_predicate()] can attach evaluation rules to it, as for a
#' subclass of SymPy's `Predicate`. `AppliedPredicate(p, ...)` applies a
#' predicate to arguments.
#'
#' @param name The predicate name, a string.
#' @param dispatch Whether the predicate has a handler for
#'   [register_predicate()].
#' @param predicate A predicate such as `Q$positive` or one made by
#'   `Predicate()`.
#' @param ... The arguments of the applied predicate.
#' @return `Predicate()` returns a callable predicate; `AppliedPredicate()`
#'   returns the applied predicate as a `sym` expression.
#' @examples
#' x <- Symbol("x")
#' P <- Predicate("P")
#' P(x)
#' ask(P(x), P(x) & Q$positive(x))
#' AppliedPredicate(Q$integer, x)
#' @export
Predicate <- function(name, dispatch = FALSE) {
  name <- as.character(name)
  if (!isTRUE(dispatch))
    return(.pred_fn(new_sym("UndefinedPredicate", name = name, extra = list(sympy_class = "UndefinedPredicate"))))
  p <- new_sym("Predicate", name = name, extra = list(sympy_class = name, custom = TRUE),
               key = paste0("CustomPredicate[", nchar(name), ":", name, "]()"))
  .pred_fn(p)
}

#' @rdname Predicate
#' @export
AppliedPredicate <- function(predicate, ...) .apply_pred(.as_pred(predicate), list(...))

#' Register an assumption handler
#'
#' Adds an evaluation rule to a predicate, as SymPy's
#' `Predicate.register(*types)` decorator does. `fn` is called as
#' `fn(expr, assumptions)` for a unary predicate, or with one argument per
#' type for a polyadic one, when the arguments have the given SymPy
#' classes (such as `"Symbol"`, `"Integer"`, `"Add"` or `"MatrixSymbol"`),
#' or subclasses of them. It returns `TRUE`, `FALSE` or `NA` for unknown.
#'
#' @param predicate A predicate made by `Predicate(name, dispatch = TRUE)`,
#'   or one of `Q`.
#' @param types A character vector of SymPy class names, one per argument.
#' @param fn The handler function.
#' @return The predicate, invisibly.
#' @examples
#' x <- Symbol("x")
#' P <- Predicate("my_key", dispatch = TRUE)
#' register_predicate(P, "Symbol", function(expr, assumptions) TRUE)
#' ask(P(x))
#' ask(P(x + 1))
#' @export
register_predicate <- function(predicate, types, fn) {
  p <- .as_pred(predicate)
  if (p$type == "UndefinedPredicate")
    .ask_type_error("<class 'sympy.assumptions.assume.UndefinedPredicate'> cannot be dispatched.")
  sig <- paste(types, collapse = ",")
  h <- .ask_custom[[p$key]]
  if (is.null(h)) h <- list(order = character(), funs = list())
  if (!(sig %in% h$order)) {
    # more specific signatures come first, as multipledispatch orders them
    depth <- function(s) sum(vapply(strsplit(s, ",", fixed = TRUE)[[1L]], function(t) length(.py_mro(t)), 0))
    o <- c(h$order, sig)
    h$order <- o[order(-vapply(o, depth, 0))]
  }
  h$funs[[sig]] <- fn
  assign(p$key, h, envir = .ask_custom)
  invisible(predicate)
}
