# Registries shared across files. Defined first so that any file may
# register hooks at load time regardless of collation order.

.rebuilders <- new.env()
register_rebuilder <- function(type, fn) assign(type, fn, envir = .rebuilders)

.str_printers <- new.env()
register_str_printer <- function(type, fn) assign(type, fn, envir = .str_printers)

.latex_printers <- new.env()
register_latex_printer <- function(type, fn) assign(type, fn, envir = .latex_printers)

# mul_post(args): Mul constructor postprocessing when noncommutative factors
# are present; pow_pre(b, e): a hook run first by the Pow constructor
.arith_hooks <- new.env()

# expand hints beyond the built in ones: .expand_custom[[hint]](e) returns
# the expanded node or NULL
.expand_custom <- new.env()

.neval_hooks <- new.env()
register_neval <- function(type, fn) assign(type, fn, envir = .neval_hooks)

.bound_hooks <- new.env()
register_bound_vars <- function(type, fn) assign(type, fn, envir = .bound_hooks)

.subs_hooks <- new.env()
register_subs_hook <- function(type, fn) assign(type, fn, envir = .subs_hooks)

.diff_hooks <- new.env()
register_diff_hook <- function(type, fn) assign(type, fn, envir = .diff_hooks)

.parse_functions <- new.env(hash = TRUE)
register_parse_function <- function(name, fn) assign(name, fn, envir = .parse_functions)

.facts_hooks <- new.env()

# String keyed hash tables for the caches. Unlike environment names, keys have
# no 10000 byte limit, lookups cost the same for any key length and clearing
# releases the memory. Kept behind these wrappers so the store can be swapped.
.htab <- function() utils::hashtab()
.hget <- function(h, k) utils::gethash(h, k)
.hset <- function(h, k, v) {
  if (utils::numhash(h) >= 50000L) utils::clrhash(h)
  utils::sethash(h, k, v)
}
.hclear <- function(h) utils::clrhash(h)
.fact_cache <- .htab()

# SymPy's @cacheit for pure functions of expressions: results keyed by the
# argument keys (and the assumptions in force), dropped wholesale when the
# cache fills up. ponytail: no LRU order; a full cache is simply cleared.
.cacheit_store <- .htab()
# global_parameters.evaluate: evaluate(FALSE, expr) turns it off
.EVAL <- new.env()
.EVAL$off <- FALSE
.ev_default <- function() !isTRUE(.EVAL$off)

.cacheit <- function(name, key, expr) {
  if (isTRUE(.EVAL$off)) return(expr)
  ctx <- global_assumptions$preds
  k <- paste(name, key, if (length(ctx)) paste(vapply(ctx, function(p) if (is.null(p$key)) "?" else p$key, ""), collapse = "&"), sep = "\r")
  hit <- .hget(.cacheit_store, k)
  if (!is.null(hit)) return(hit[[1L]])
  v <- expr
  .hset(.cacheit_store, k, list(v))
  v
}
.ckey <- function(z) if (is.null(z)) "NULL" else if (is_sym(z)) z$key else paste(format(z), collapse = ",")
.fs_cache <- .htab()
.ops_hooks <- new.env()
.expand_hooks <- new.env()
.doit_hooks <- new.env()
register_doit <- function(type, fn) assign(type, fn, envir = .doit_hooks)

.REL_TYPES <- c("Equality", "Unequality", "StrictLessThan", "LessThan",
                "StrictGreaterThan", "GreaterThan")

.FUNS <- new.env(hash = TRUE)

def_function <- function(name, nargs = 1L, eval = NULL, diff = NULL,
                         numeric = NULL, facts = NULL, latex = NULL,
                         str = NULL, ...) {
  assign(name, list(name = name, nargs = nargs, eval = eval, diff = diff,
                    numeric = numeric, facts = facts, latex = latex,
                    str = str, ...), envir = .FUNS)
  invisible(name)
}

#' Clear the caches
#'
#' SymPy's `clear_cache()`: empties the caches of `as_leading_term`,
#' `powsimp`, `Order.contains` and the other cached computations.
#'
#' @return `NULL`, invisibly.
#' @examples
#' clear_cache()
#' @export
clear_cache <- function() {
  for (h in list(.cacheit_store, .fact_cache, .fs_cache, .comm_tab)) .hclear(h)
  invisible(NULL)
}

.cacheit_ids <- new.env()
.cacheit_ids$n <- 0L

#' Cache a function
#'
#' SymPy's `cacheit` decorator: `cacheit(f)` returns a function that
#' remembers its results by the arguments, in the cache that
#' [clear_cache()] empties.
#'
#' @param f A function.
#' @return The caching function.
#' @examples
#' slow_square <- function(x) { Sys.sleep(0.1); x^2 }
#' fast <- cacheit(slow_square)
#' fast(Symbol("x")); fast(Symbol("x"))
#' @export
cacheit <- function(f) {
  .cacheit_ids$n <- .cacheit_ids$n + 1L
  id <- paste0("cacheit:", .cacheit_ids$n)
  function(...) {
    a <- list(...)
    .cacheit(id, paste(vapply(a, .ckey, ""), names(a), collapse = "\t"), f(...))
  }
}
