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
.fact_cache <- new.env(hash = TRUE)

# SymPy's @cacheit for pure functions of expressions: results keyed by the
# argument keys (and the assumptions in force), dropped wholesale when the
# cache fills up. ponytail: no LRU order; a full cache is simply cleared.
.cacheit_store <- new.env(hash = TRUE, parent = emptyenv())
.cacheit <- function(name, key, expr) {
  ctx <- global_assumptions$preds
  k <- paste(name, key, if (length(ctx)) paste(vapply(ctx, function(p) if (is.null(p$key)) "?" else p$key, ""), collapse = "&"), sep = "\r")
  hit <- .cacheit_store[[k]]
  if (!is.null(hit)) return(hit[[1L]])
  v <- expr
  if (length(.cacheit_store) >= 50000L) rm(list = ls(.cacheit_store, all.names = TRUE), envir = .cacheit_store)
  .cacheit_store[[k]] <- list(v)
  v
}
.ckey <- function(z) if (is.null(z)) "NULL" else if (is_sym(z)) z$key else paste(format(z), collapse = ",")
.fs_cache <- new.env(hash = TRUE)
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
