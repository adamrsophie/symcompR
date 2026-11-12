# Small utilities after sympy/utilities (misc.py, decorator.py,
# iterables.py, timeutils.py) and sympy/interactive.

#' Utility helpers
#'
#' `filldedent(s, w)` strips leading and trailing empty lines, removes the
#' common indentation and refills the text to width `w`, starting with a
#' newline, as SymPy's `filldedent()`. `capture(func)` returns what
#' `func()` prints. `threaded(func)` and `xthreaded(func)` return a
#' function that maps `func` over matrices, lists and the two sides of a
#' relation; `threaded()` also maps over the terms of a sum. `public(obj)`
#' returns `obj` (R namespaces are declared in `NAMESPACE`).
#' `memoize_property(f)` caches `f(x)` per argument. `timed(func)` times
#' `func()` adaptively and returns the number of runs, the time per run in
#' seconds, and that time in the unit given last.
#'
#' @param s A string.
#' @param w The line width.
#' @param func,f A function.
#' @param obj Any object.
#' @param limit The largest number of runs per timing in `timed()`.
#' @return A string, a function, or a list.
#' @examples
#' filldedent("
#'     A long docstring line
#'     continues here.")
#' capture(function() print(Symbol("x")^2))
#' d <- threaded(function(e) diff(e, Symbol("x")))
#' d(list(Symbol("x")^2, sin(Symbol("x"))))
#' @name utility-helpers
NULL

#' @rdname utility-helpers
#' @export
filldedent <- function(s, w = 70) {
  lines <- strsplit(as.character(s), "\n", fixed = TRUE)[[1L]]
  # textwrap.dedent: drop the common leading whitespace of non-blank lines
  nb <- lines[nzchar(trimws(lines))]
  ind <- if (length(nb)) min(nchar(sub("^([ \t]*).*$", "\\1", nb))) else 0L
  lines <- substring(lines, ind + 1L)
  txt <- gsub("^\n+|\n+$", "", paste(lines, collapse = "\n"))
  paste0("\n", paste(strwrap(txt, width = w + 1L), collapse = "\n"))
}

#' @rdname utility-helpers
#' @export
capture <- function(func) {
  out <- utils::capture.output(invisible(func()))
  if (length(out)) paste0(paste(out, collapse = "\n"), "\n") else ""
}

.threaded_factory <- function(func, use_add) {
  force(func)
  function(expr, ...) {
    if (inherits(expr, "symmatrix")) return(.mat_map(expr, function(z) func(z, ...)))
    if (is.list(expr) && !is_sym(expr)) return(lapply(expr, function(z) func(z, ...)))
    expr <- as_sym(expr)
    if (use_add && expr$type == "Add") return(.add(lapply(expr$args, function(z) func(z, ...))))
    if (is_relational(expr)) return(sym_rebuild(expr, list(func(expr$args[[1L]], ...), func(expr$args[[2L]], ...))))
    func(expr, ...)
  }
}

#' @rdname utility-helpers
#' @export
threaded <- function(func) .threaded_factory(func, TRUE)

#' @rdname utility-helpers
#' @export
xthreaded <- function(func) .threaded_factory(func, FALSE)

#' @rdname utility-helpers
#' @export
public <- function(obj) {
  if (!is.function(obj)) stop("expected a function, got ", class(obj)[1L], call. = FALSE)
  obj
}

#' @rdname utility-helpers
#' @export
memoize_property <- function(f) {
  force(f)
  cache <- new.env(hash = TRUE)
  function(x) {
    k <- if (is_sym(x)) x$key else paste(deparse(x), collapse = "")
    if (exists(k, envir = cache, inherits = FALSE)) return(get(k, envir = cache))
    v <- f(x)
    assign(k, v, envir = cache)
    v
  }
}

#' @rdname utility-helpers
#' @export
timed <- function(func, limit = NULL) {
  run <- function(number) { t0 <- proc.time()[["elapsed"]]; for (i in seq_len(number)) func(); proc.time()[["elapsed"]] - t0 }
  number <- 1
  for (i in 1:9) {
    if (run(number) >= 0.2) break
    if (!is.null(limit) && number >= limit) break
    number <- number * 10
  }
  time <- min(vapply(1:3, function(i) run(number), 0)) / number
  units <- c("s", "ms", "us", "ns"); scales <- c(1, 1e3, 1e6, 1e9)
  order <- if (time > 0) min(-floor(log10(time)) %/% 3, 3) else 3
  order <- max(order, 0)
  list(number = number, time = time, scaled = time * scales[order + 1], unit = units[order + 1])
}

#' Interactive sessions
#'
#' `init_printing()` chooses how symbolic expressions print at the console,
#' as SymPy's `init_printing()`: as two dimensional pretty printing
#' (`pretty_print = TRUE`, the default), as LaTeX (`use_latex = TRUE`), or as
#' plain strings (`pretty_print = FALSE`). `init_session()` also defines the
#' symbols `x`, `y`, `z`, `t`, the integer symbols `k`, `m`, `n` and the
#' functions `f`, `g`, `h` in `envir`, as SymPy's session does.
#' `interactive_traversal(expr)` walks an expression tree, asking at each
#' node which argument to enter, and returns the chosen subexpression.
#'
#' @param pretty_print Pretty print expressions.
#' @param use_unicode Unicode characters in pretty printing (`NULL` keeps
#'   the current choice).
#' @param use_latex Print LaTeX.
#' @param quiet Do not print the banner.
#' @param envir The environment receiving the symbols.
#' @param expr An expression.
#' @return `init_printing()` and `init_session()` return `NULL` invisibly;
#'   `interactive_traversal()` returns an expression.
#' @examples
#' init_printing(pretty_print = FALSE)
#' interactive_traversal(Symbol("x") + 1)
#' @name interactive-session
NULL

#' @rdname interactive-session
#' @export
init_printing <- function(pretty_print = TRUE, use_unicode = NULL, use_latex = FALSE) {
  mode <- if (isTRUE(use_latex)) "latex" else if (isTRUE(pretty_print)) "pretty" else "str"
  options(symcompR.print = mode)
  if (!is.null(use_unicode)) pprint_use_unicode(use_unicode)
  invisible(NULL)
}

#' @rdname interactive-session
#' @export
init_session <- function(pretty_print = TRUE, use_unicode = NULL, quiet = FALSE, envir = parent.frame()) {
  init_printing(pretty_print = pretty_print, use_unicode = use_unicode)
  for (nm in c("x", "y", "z", "t")) assign(nm, Symbol(nm), envir = envir)
  for (nm in c("k", "m", "n")) assign(nm, Symbol(nm, integer = TRUE), envir = envir)
  for (nm in c("f", "g", "h")) assign(nm, Function(nm), envir = envir)
  if (!quiet) message("These commands were executed:\n",
                      "x, y, z, t <- symbols; k, m, n <- integer symbols; f, g, h <- functions\n",
                      "init_printing(pretty_print = ", pretty_print, ")")
  invisible(NULL)
}

#' @rdname interactive-session
#' @export
interactive_traversal <- function(expr) {
  expr <- as_sym(expr)
  repeat {
    if (!interactive() || !length(expr$args)) return(expr)
    cat("Expression:", .str(expr), "\n")
    for (i in seq_along(expr$args)) cat(sprintf("  %d: %s\n", i, .str(expr$args[[i]])))
    ans <- readline("Choose an argument (empty to stop): ")
    k <- suppressWarnings(as.integer(ans))
    if (is.na(k) || k < 1L || k > length(expr$args)) return(expr)
    expr <- expr$args[[k]]
  }
}
