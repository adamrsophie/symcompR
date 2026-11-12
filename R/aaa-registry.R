# Registries shared across files. Defined first so that any file may
# register hooks at load time regardless of collation order.

.rebuilders <- new.env()
register_rebuilder <- function(type, fn) assign(type, fn, envir = .rebuilders)

.str_printers <- new.env()
register_str_printer <- function(type, fn) assign(type, fn, envir = .str_printers)

.latex_printers <- new.env()
register_latex_printer <- function(type, fn) assign(type, fn, envir = .latex_printers)

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
