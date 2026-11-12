# Assumption contexts, after sympy/assumptions/assume.py: a set of
# assumptions that ask() and refine() add to the ones they are given.

#' Assumption contexts
#'
#' `AssumptionsContext(...)` is a mutable set of assumptions with
#' `$add(...)`, `$remove(...)` and `$clear()`; its `$preds` field lists
#' them. `global_assumptions` is the context that [ask()] and [refine()]
#' consult by default, as in SymPy. `assuming(..., code)` evaluates its
#' last argument with the other arguments added to `global_assumptions`,
#' and restores the previous state afterwards, like SymPy's
#' `with assuming(...):` block. Removing an assumption that is not in the
#' context raises a `KeyError`, as for a Python set.
#'
#' @param ... Assumptions such as `Q$positive(x)` or `x > 0`; for
#'   `assuming()`, the last argument is the code to evaluate.
#' @return `AssumptionsContext()` returns a context; `assuming()` returns the
#'   value of the code.
#' @examples
#' x <- Symbol("x")
#' assuming(Q$positive(x), ask(Q$real(x)))
#' assuming(Q$positive(x), refine(sqrt(x^2)))
#' ask(Q$real(x))
#' ctx <- AssumptionsContext(Q$integer(x))
#' ctx
#' ask(Q$rational(x), context = ctx)
#' @name assumptions-context
NULL

#' @rdname assumptions-context
#' @export
AssumptionsContext <- function(...) {
  ctx <- new.env(parent = emptyenv())
  ctx$preds <- list()
  ctx$add <- function(...) {
    for (p in list(...)) {
      p <- .ask_sympify(p)
      ctx$preds[[p$key]] <- p
    }
    invisible(NULL)
  }
  ctx$remove <- function(...) {
    for (p in list(...)) {
      p <- .ask_sympify(p)
      if (is.null(ctx$preds[[p$key]])) stop(.sympy_condition("KeyError", .str(p)))
      ctx$preds[[p$key]] <- NULL
    }
    invisible(NULL)
  }
  ctx$clear <- function() {
    ctx$preds <- list()
    invisible(NULL)
  }
  ctx$add(...)
  class(ctx) <- "AssumptionsContext"
  ctx
}

.ctx_preds <- function(ctx) unname(ctx$preds)
.ctx_contains <- function(ctx, p) !is.null(ctx$preds[[.ask_sympify(p)$key]])

#' @export
format.AssumptionsContext <- function(x, ...) {
  p <- .ctx_preds(x)
  if (!length(p)) return("AssumptionsContext()")
  p <- .key_sort(lapply(p, .sort_key), p)
  paste0("AssumptionsContext({", paste(vapply(p, .str, ""), collapse = ", "), "})")
}

#' @export
print.AssumptionsContext <- function(x, ...) {
  cat(format(x), "\n", sep = "")
  invisible(x)
}

#' @export
length.AssumptionsContext <- function(x) length(x$preds)

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
