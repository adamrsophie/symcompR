# Assumption contexts, after sympy/assumptions/assume.py: a set of
# predicates that ask() and refine() add to the ones they are given.

#' Assumption contexts
#'
#' `AssumptionsContext(...)` is a mutable set of predicates with `$add(...)`,
#' `$remove(...)` and `$clear()`. `global_assumptions` is the context that
#' [ask()] and [refine()] always consult, as in SymPy.
#' `assuming(..., code)` evaluates its last argument with the other arguments
#' added to `global_assumptions`, and restores the previous state afterwards,
#' like SymPy's `with assuming(...):` block.
#'
#' @param ... Predicates such as `Q$positive(x)`; for `assuming()`, the last
#'   argument is the code to evaluate.
#' @return `AssumptionsContext()` returns a context; `assuming()` returns the
#'   value of the code.
#' @examples
#' x <- Symbol("x")
#' assuming(Q$positive(x), ask(Q$real(x)))
#' assuming(Q$positive(x), refine(sqrt(x^2)))
#' ask(Q$real(x))
#' @name assumptions-context
NULL

.pred_key <- function(p) paste0(if (isFALSE(p$value)) "~", p$fact, ":", p$expr$key)

#' @rdname assumptions-context
#' @export
AssumptionsContext <- function(...) {
  ctx <- new.env(parent = emptyenv())
  ctx$preds <- list()
  ctx$add <- function(...) {
    for (p in list(...)) {
      if (!inherits(p, "sym_predicate")) stop("assumptions must be predicates such as Q$positive(x)", call. = FALSE)
      ctx$preds[[.pred_key(p)]] <- p
    }
    invisible(ctx)
  }
  ctx$remove <- function(...) {
    for (p in list(...)) {
      k <- .pred_key(p)
      if (is.null(ctx$preds[[k]])) stop("predicate not in the context: ", .pred_str(p), call. = FALSE)
      ctx$preds[[k]] <- NULL
    }
    invisible(ctx)
  }
  ctx$clear <- function() { ctx$preds <- list(); invisible(ctx) }
  ctx$add(...)
  class(ctx) <- "AssumptionsContext"
  ctx
}

.pred_str <- function(p) paste0(if (isFALSE(p$value)) "~", "Q.", p$fact, "(", .str(p$expr), ")")

#' @export
print.AssumptionsContext <- function(x, ...) {
  cat("AssumptionsContext(", paste(vapply(unname(x$preds), .pred_str, ""), collapse = ", "), ")\n", sep = "")
  invisible(x)
}

#' @rdname assumptions-context
#' @export
global_assumptions <- AssumptionsContext()

#' @rdname assumptions-context
#' @export
assuming <- function(...) {
  cl <- as.list(substitute(list(...)))[-1L]
  if (!length(cl)) stop("assuming() needs code to evaluate", call. = FALSE)
  env <- parent.frame()
  # the arguments are the caller's own code, evaluated in the caller's frame
  # as base::local() does; nothing comes from strings or outside input
  preds <- lapply(cl[-length(cl)], eval, env)
  old <- global_assumptions$preds
  on.exit(global_assumptions$preds <- old)
  do.call(global_assumptions$add, preds)
  eval(cl[[length(cl)]], env)
}
