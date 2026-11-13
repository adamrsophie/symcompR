# Branching strategies, after sympy/strategies/branch (core.py, tools.py and
# traverse.py). A branching rule maps an expression to a list of results, the
# values SymPy's generators yield, in the same order.

.prt_bkey <- function(x) if (is_sym(x)) x$key else paste(deparse(x), collapse = "")
.prt_bstr <- function(x) if (is_sym(x)) sstr(x) else if (is.character(x)) x else paste(format(x), collapse = " ")

#' Branching strategies
#'
#' Combinators for branching rules after `sympy.strategies.branch`. A
#' branching rule maps an expression to a list of results (SymPy's
#' generators). `branch_identity()` returns its input, `branch_exhaust()`
#' applies a rule until it has no effect, `onaction()` calls a function on
#' each result that differs from the input and `branch_debug()` prints the
#' rule, input and output of each, `multiplex()` merges the distinct results
#' of several rules, `branch_condition()` applies a rule when a predicate
#' holds, `sfilter()` keeps the results that satisfy a predicate,
#' `notempty()` returns the input when a rule gives nothing,
#' `branch_do_one()` uses the first rule that gives results,
#' `branch_chain()` composes rules and `yieldify()` turns a rule into a
#' branching rule. `branch_top_down()` applies a rule down an expression tree,
#' `branch_sall()` applies it to the arguments, and `branch_canon()` applies
#' rules top down until nothing changes.
#'
#' @param x An expression.
#' @param brule,... Branching rules.
#' @param fn A function `fn(brule, expr, result)`.
#' @param file `NULL` to print to the console, or a connection to write to.
#' @param cond,pred A predicate.
#' @param rl A rule.
#' @param fns Traversal functions (`op`, `new`, `children`, `leaf`), as
#'   [basic_fns].
#' @return A branching rule; `branch_identity()` returns a list.
#' @examples
#' posdec <- function(x) if (x > 0) list(x - 1L) else list(x)
#' branch5 <- function(x) {
#'   if (0 < x && x < 5) list(x - 1L) else if (5 < x && x < 10) list(x + 1L)
#'   else if (x == 5) list(x + 1L, x - 1L) else list(x)
#' }
#' branch_exhaust(branch5)(5L)
#' multiplex(posdec, branch5)(5L)
#' branch_chain(branch5, yieldify(function(x) x * 2L))(5L)
#' x <- Symbol("x")
#' inc <- function(e) if (is_sym(e) && e$type == "Integer") list(e + 1L) else list()
#' branch_top_down(inc)(x + 1)
#' @name branch-strategies
NULL

#' @rdname branch-strategies
#' @export
branch_identity <- function(x) list(x)

#' @rdname branch-strategies
#' @export
branch_exhaust <- function(brule) {
  force(brule)
  rec <- function(expr) {
    seen <- .prt_bkey(expr)
    out <- list()
    for (ne in brule(expr)) {
      k <- .prt_bkey(ne)
      if (!(k %in% seen)) {
        seen <- c(seen, k)
        out <- c(out, rec(ne))
      }
    }
    if (length(seen) == 1L) out <- c(out, list(expr))
    out
  }
  rec
}

#' @rdname branch-strategies
#' @export
onaction <- function(brule, fn) {
  force(brule); force(fn)
  function(expr) {
    out <- brule(expr)
    for (r in out) if (!.strat_eq(r, expr)) fn(brule, expr, r)
    out
  }
}

#' @rdname branch-strategies
#' @export
branch_debug <- function(brule, file = NULL) {
  name <- paste(deparse(substitute(brule)), collapse = "")
  force(brule)
  write <- function(brl, expr, result) {
    s <- sprintf("Rule: %s\nIn: %s\nOut: %s\n\n", name, .prt_bstr(expr), .prt_bstr(result))
    if (is.null(file)) cat(s) else cat(s, file = file, sep = "")
  }
  onaction(brule, write)
}

#' @rdname branch-strategies
#' @export
multiplex <- function(...) {
  brules <- list(...)
  function(expr) {
    seen <- character(); out <- list()
    for (b in brules) for (ne in b(expr)) {
      k <- .prt_bkey(ne)
      if (!(k %in% seen)) { seen <- c(seen, k); out[[length(out) + 1L]] <- ne }
    }
    out
  }
}

#' @rdname branch-strategies
#' @export
branch_condition <- function(cond, brule) {
  force(cond); force(brule)
  function(expr) if (isTRUE(cond(expr))) brule(expr) else list()
}

#' @rdname branch-strategies
#' @export
sfilter <- function(pred, brule) {
  force(pred); force(brule)
  function(expr) Filter(function(r) isTRUE(pred(r)), brule(expr))
}

#' @rdname branch-strategies
#' @export
notempty <- function(brule) {
  force(brule)
  function(expr) {
    out <- brule(expr)
    if (length(out)) out else list(expr)
  }
}

#' @rdname branch-strategies
#' @export
branch_do_one <- function(...) {
  brules <- list(...)
  function(expr) {
    for (b in brules) {
      out <- b(expr)
      if (length(out)) return(out)
    }
    list()
  }
}

#' @rdname branch-strategies
#' @export
branch_chain <- function(...) {
  brules <- list(...)
  function(expr) {
    if (!length(brules)) return(list(expr))
    tail <- do.call(branch_chain, brules[-1L])
    out <- list()
    for (ne in brules[[1L]](expr)) out <- c(out, tail(ne))
    out
  }
}

#' @rdname branch-strategies
#' @export
yieldify <- function(rl) {
  force(rl)
  function(expr) list(rl(expr))
}

#' @rdname branch-strategies
#' @export
branch_canon <- function(...) branch_exhaust(do.call(multiplex, lapply(list(...), branch_top_down)))

#' @rdname branch-strategies
#' @export
branch_top_down <- function(brule, fns = basic_fns) {
  force(brule); force(fns)
  branch_chain(branch_do_one(brule, branch_identity), function(expr) branch_sall(branch_top_down(brule, fns), fns)(expr))
}

# itertools.product of lists, the last one varying fastest
.prt_product <- function(opts) {
  out <- list(list())
  for (o in opts) {
    nxt <- list()
    for (pre in out) for (v in o) nxt[[length(nxt) + 1L]] <- c(pre, list(v))
    out <- nxt
  }
  out
}

#' @rdname branch-strategies
#' @export
branch_sall <- function(brule, fns = basic_fns) {
  force(brule); force(fns)
  function(expr) {
    if (isTRUE(fns$leaf(expr))) return(list(expr))
    myop <- fns$op(expr)
    lapply(.prt_product(lapply(fns$children(expr), brule)), function(args) do.call(fns$new, c(list(myop), args)))
  }
}
