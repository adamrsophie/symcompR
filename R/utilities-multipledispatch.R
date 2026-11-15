# Multiple dispatch, after sympy/multipledispatch (dispatcher.py, core.py,
# conflict.py). A Dispatcher holds implementations keyed by signatures of
# class names. An argument matches a signature entry when the name is in
# its SymPy method resolution order (.cgp_mro) or in its R class vector;
# among matching signatures the most specific one wins, where a signature
# is more specific than another when each of its entries comes no later in
# the argument's resolution order. Ties raise SymPy's AmbiguityWarning.

.md_order <- function(a) {
  m <- tryCatch(.cgp_mro(a, full = TRUE), error = function(e) character())
  unique(c(m, class(a), if (is.numeric(a)) c("int", "float", "Number"), "object", "ANY"))
}

.md_ns <- new.env()
.md_state <- new.env()
.md_state$halted <- FALSE

#' Multiple dispatch
#'
#' A port of `sympy.multipledispatch`. `Dispatcher(name)` makes a function
#' that dispatches on the classes of all its arguments; `dispatcher_add(d,
#' signature, fun)` registers an implementation for a vector of class names,
#' which may be SymPy class names (`"Integer"`, `"Expr"`, `"Basic"`), R
#' classes, or `"object"` for any value. `dispatch(..., fun)` registers
#' `fun` under its name in a namespace and returns the dispatcher, as the
#' decorator does. An implementation that raises `MDNotImplementedError`
#' passes the call on to the next most specific one.
#' `halt_ordering()` and `restart_ordering()` suspend and resume the
#' ambiguity warnings while many implementations are added.
#'
#' @param name The function name.
#' @param d A dispatcher.
#' @param signature A character vector of class names.
#' @param fun The implementation.
#' @param ... Class names of the signature.
#' @param fname The name to register `fun` under.
#' @param namespace An environment that holds the dispatchers.
#' @param message The error message.
#' @return A dispatcher (a function); `dispatcher_add()` returns it
#'   invisibly.
#' @examples
#' f <- dispatch("Integer", fname = "f", fun = function(x) x + 1)
#' f <- dispatch("Rational", fname = "f", fun = function(x) x * 2)
#' f(Integer(3))
#' f(Rational(1, 3))
#' @name multipledispatch
#' @export
Dispatcher <- function(name) {
  env <- new.env()
  env$name <- name
  env$funcs <- list()
  env$sigs <- list()
  f <- function(...) .md_call(env, list(...))
  structure(f, class = c("Dispatcher", "function"), env = env)
}

#' @export
print.Dispatcher <- function(x, ...) { cat("<dispatched ", attr(x, "env")$name, ">\n", sep = ""); invisible(x) }

#' @rdname multipledispatch
#' @export
dispatcher_add <- function(d, signature, fun) {
  env <- attr(d, "env")
  key <- paste(signature, collapse = ",")
  i <- match(key, vapply(env$sigs, paste, "", collapse = ","))
  if (is.na(i)) i <- length(env$sigs) + 1L
  env$sigs[[i]] <- as.character(signature)
  env$funcs[[i]] <- fun
  invisible(d)
}

.md_candidates <- function(env, args) {
  ords <- lapply(args, .md_order)
  sc <- list(); idx <- integer()
  for (i in seq_along(env$sigs)) {
    s <- env$sigs[[i]]
    if (length(s) != length(args)) next
    pos <- vapply(seq_along(s), function(k) { p <- match(s[[k]], ords[[k]]); if (is.na(p)) NA_integer_ else p }, 0L)
    if (anyNA(pos)) next
    sc[[length(sc) + 1L]] <- pos; idx <- c(idx, i)
  }
  if (!length(idx)) return(integer())
  # order by dominance: a candidate comes before every candidate it beats
  dominates <- function(a, b) all(a <= b) && any(a < b)
  ord <- integer(); rest <- seq_along(idx)
  while (length(rest)) {
    top <- rest[vapply(rest, function(r) !any(vapply(rest, function(o) o != r && dominates(sc[[o]], sc[[r]]), NA)), NA)]
    if (length(top) > 1L) {
      t1 <- top[order(vapply(top, function(r) sum(sc[[r]]), 0))]
      if (!.md_state$halted && sum(sc[[t1[1L]]]) == sum(sc[[t1[2L]]]))
        warning("AmbiguityWarning: ambiguous signatures for ", env$name, ": ",
                paste(vapply(t1[1:2], function(r) paste0("(", paste(env$sigs[[idx[r]]], collapse = ", "), ")"), ""), collapse = " and "),
                call. = FALSE)
      top <- t1
    }
    ord <- c(ord, top); rest <- setdiff(rest, top)
  }
  idx[ord]
}

.md_call <- function(env, args) {
  cand <- .md_candidates(env, args)
  if (!length(cand))
    stop("NotImplementedError: Could not find signature for ", env$name, ": <",
         paste(vapply(args, function(a) .md_order(a)[1L], ""), collapse = ", "), ">", call. = FALSE)
  for (i in cand) {
    r <- tryCatch(list(v = do.call(env$funcs[[i]], args)), MDNotImplementedError = function(e) NULL)
    if (!is.null(r)) return(r$v)
  }
  stop("NotImplementedError: Matching functions for ", env$name, " found, but none completed successfully", call. = FALSE)
}

#' @rdname multipledispatch
#' @export
dispatch <- function(..., fun, fname = NULL, namespace = NULL) {
  if (is.null(namespace)) namespace <- .md_ns
  if (is.null(fname)) fname <- deparse(substitute(fun))
  d <- if (exists(fname, envir = namespace, inherits = FALSE)) get(fname, envir = namespace) else Dispatcher(fname)
  assign(fname, d, envir = namespace)
  dispatcher_add(d, c(...), fun)
  d
}

#' @rdname multipledispatch
#' @export
MDNotImplementedError <- function(message = "") .sympy_condition("MDNotImplementedError", message)

#' @rdname multipledispatch
#' @export
halt_ordering <- function() { .md_state$halted <- TRUE; invisible(NULL) }

#' @rdname multipledispatch
#' @export
restart_ordering <- function() { .md_state$halted <- FALSE; invisible(NULL) }
