# Rule-based programming, after sympy/strategies: rules are functions from
# an expression to an expression; branching rules return a list of
# results (SymPy's generators). Trees for greedy() and allresults() use
# plain lists for alternatives (SymPy lists) and strat_seq() for sequences
# (SymPy tuples).

#' Strategies
#'
#' Combinators for rewriting rules after `sympy.strategies`. A rule maps an
#' expression to an expression. `chain()` applies rules in order,
#' `do_one()` applies the first rule that changes the expression,
#' `exhaust()` repeats a rule until nothing changes, `condition()` applies a
#' rule when a predicate holds, `null_safe()` keeps the input when a rule
#' returns `NULL`, `tryit()` keeps it when a rule fails, `minimize()` keeps
#' the result with the smallest objective, `memoize()` caches results,
#' `strat_switch()` picks a rule by a key, and `canon()` applies rules top
#' down to a fixed point. The lists `rl` (expression rules `rm_id`, `glom`,
#' `sort`, `distribute`, `subs`, `unpack`, `flatten`, `rebuild`),
#' `traverse` (`top_down`, `bottom_up`, `top_down_once`, `bottom_up_once`,
#' `sall`) and `branch` (the same combinators for branching rules) mirror
#' SymPy's modules. `greedy(tree, objective)` evaluates a tree of rules
#' (lists are alternatives, `strat_seq()` are sequences) keeping the best
#' alternative, `allresults()` returns every result and `brute()` the best
#' of all.
#'
#' @param ... Rules.
#' @param rule A rule.
#' @param cond A predicate.
#' @param objective A function giving a number or a vector compared
#'   lexicographically.
#' @param key A function giving a name.
#' @param ruledict A named list of rules.
#' @param exception Ignored (all errors are caught).
#' @param ruletypes A named list of rules by expression type.
#' @param tree A tree of rules.
#' @param join A named list with functions for `list` and `seq` nodes.
#' @param leaf A function applied to the rules at the leaves.
#' @param fns Traversal functions (`op`, `new`, `children`, `leaf`).
#' @param file Ignored; debugging output goes to the console.
#' @param x An expression (`strat_identity()` returns it unchanged).
#' @param isid For `rm_id()`: a predicate for identity elements.
#' @param new A constructor `new(op, ...)` for rebuilt expressions.
#' @param expr An expression.
#' @param count,combine For `glom()`: a function giving a count and one
#'   combining the counts.
#' @param A,B For `distribute()`: the outer and inner operations.
#' @param op An operation name.
#' @return A rule, or for `allresults()` a branching rule.
#' @examples
#' x <- Symbol("x")
#' inc <- function(e) if (isTRUE(is_integer(e)) && as.numeric(e) < 5) e + 1 else e
#' exhaust(inc)(Integer(1))
#' chain(inc, inc)(Integer(0))
#' greedy(list(function(e) expand(e), function(e) sym_factor(e)), objective = count_ops)((x + 1)^2)
#' @name strategies
NULL

#' @rdname strategies
#' @export
strat_identity <- function(x) x

#' @rdname strategies
#' @export
strat_seq <- function(...) structure(list(...), class = "strat_seq")

.strat_eq <- function(a, b) {
  if (is_sym(a) && is_sym(b)) return(identical(a$key, b$key))
  identical(a, b)
}

#' @rdname strategies
#' @export
exhaust <- function(rule) {
  force(rule)
  function(expr) {
    new <- rule(expr); old <- expr
    while (!.strat_eq(new, old)) { old <- new; new <- rule(new) }
    new
  }
}

#' @rdname strategies
#' @export
memoize <- function(rule) {
  force(rule)
  cache <- new.env(hash = TRUE)
  function(expr) {
    k <- if (is_sym(expr)) expr$key else paste(deparse(expr), collapse = "")
    if (!is.null(cache[[k]])) return(cache[[k]])
    r <- rule(expr); assign(k, r, envir = cache); r
  }
}

#' @rdname strategies
#' @export
condition <- function(cond, rule) { force(cond); force(rule); function(expr) if (isTRUE(cond(expr))) rule(expr) else expr }

#' @rdname strategies
#' @export
chain <- function(...) {
  rules <- list(...)
  function(expr) { for (r in rules) expr <- r(expr); expr }
}

#' @rdname strategies
#' @export
strat_debug <- function(rule, file = NULL) {
  nm <- paste(deparse(substitute(rule)), collapse = "")
  function(...) {
    expr <- list(...)[[1L]]
    result <- rule(...)
    if (!.strat_eq(result, expr))
      cat("Rule: ", nm, "\nIn:   ", format(expr), "\nOut:  ", format(result), "\n\n", sep = "")
    result
  }
}

#' @rdname strategies
#' @export
null_safe <- function(rule) { force(rule); function(expr) { r <- rule(expr); if (is.null(r)) expr else r } }

#' @rdname strategies
#' @export
tryit <- function(rule, exception = NULL) { force(rule); function(expr) tryCatch(rule(expr), error = function(err) expr) }

#' @rdname strategies
#' @export
do_one <- function(...) {
  rules <- list(...)
  function(expr) {
    for (r in rules) { res <- r(expr); if (!.strat_eq(res, expr)) return(res) }
    expr
  }
}

#' @rdname strategies
#' @export
strat_switch <- function(key, ruledict) {
  force(key); force(ruledict)
  function(expr) {
    k <- key(expr)
    r <- if (!is.null(k) && !is.null(ruledict[[k]])) ruledict[[k]] else strat_identity
    r(expr)
  }
}

# lexicographic comparison of objective values (numbers, logicals, vectors)
.obj_lt <- function(a, b) {
  a <- unlist(lapply(a, function(v) if (is_sym(v)) num_value_dbl(v) else as.numeric(v)))
  b <- unlist(lapply(b, function(v) if (is_sym(v)) num_value_dbl(v) else as.numeric(v)))
  for (i in seq_len(min(length(a), length(b)))) if (a[i] != b[i]) return(a[i] < b[i])
  length(a) < length(b)
}

#' @rdname strategies
#' @export
minimize <- function(..., objective = strat_identity) {
  rules <- list(...)
  function(expr) {
    rs <- lapply(rules, function(r) r(expr))
    best <- 1L; bv <- objective(rs[[1L]])
    for (i in seq_along(rs)[-1L]) { v <- objective(rs[[i]]); if (.obj_lt(v, bv)) { best <- i; bv <- v } }
    rs[[best]]
  }
}

# ---------------------------------------------------------------------------
# Expression rules (strategies/rl.py) and traversal (traverse.py)

.strat_op <- function(e) list(type = e$type, name = e$name, extra = e$extra, value = e$value)
.strat_new <- function(op, args) new_sym(op$type, args, name = op$name, value = op$value, extra = op$extra)
.strat_leaf <- function(e) !is_sym(e) || !length(e$args)

#' @rdname strategies
#' @export
basic_fns <- list(op = .strat_op, new = function(op, ...) .strat_new(op, list(...)),
                  leaf = .strat_leaf, children = function(e) e$args)
#' @rdname strategies
#' @export
expr_fns <- list(op = .strat_op, new = function(op, ...) sym_rebuild(.strat_new(op, list(...)), list(...)),
                 leaf = .strat_leaf, children = function(e) e$args)

#' @rdname strategies
#' @export
rl <- list(
  rm_id = function(isid, new = function(op, ...) .strat_new(op, list(...))) function(expr) {
    ids <- vapply(expr$args, function(a) isTRUE(isid(a)), TRUE)
    if (!any(ids)) return(expr)
    if (!all(ids)) return(do.call(new, c(list(.strat_op(expr)), expr$args[!ids])))
    do.call(new, list(.strat_op(expr), expr$args[[1L]]))
  },
  glom = function(key, count, combine) function(expr) {
    ks <- lapply(expr$args, key); kk <- vapply(ks, function(k) if (is_sym(k)) k$key else paste(deparse(k), collapse = ""), "")
    groups <- split(seq_along(kk), factor(kk, levels = unique(kk)))
    newargs <- lapply(groups, function(ix) combine(sum(vapply(expr$args[ix], function(a) as.numeric(count(a)), 0)), ks[[ix[1L]]]))
    newargs <- unname(newargs)
    if (setequal(vapply(newargs, .k, ""), vapply(expr$args, .k, ""))) return(expr)
    .strat_new(.strat_op(expr), newargs)
  },
  sort = function(key, new = function(op, ...) .strat_new(op, list(...))) function(expr) {
    ks <- lapply(expr$args, key)
    o <- if (all(vapply(ks, function(k) is.numeric(k) || is.character(k), TRUE)) && all(lengths(ks) == 1L)) order(unlist(ks))
         else order(vapply(ks, function(k) paste(format(k), collapse = ""), ""))
    do.call(new, c(list(.strat_op(expr)), expr$args[o]))
  },
  distribute = function(A, B) function(expr) {
    for (i in seq_along(expr$args)) {
      a <- expr$args[[i]]
      if (a$type == B) {
        first <- expr$args[seq_len(i - 1L)]; tail <- expr$args[-seq_len(i)]
        return(sym_rebuild(a, lapply(a$args, function(b) sym_rebuild(expr, c(first, list(b), tail)))))
      }
    }
    expr
  },
  subs = function(a, b) function(expr) if (.strat_eq(expr, a)) b else expr,
  unpack = function(expr) if (length(expr$args) == 1L) expr$args[[1L]] else expr,
  flatten = function(expr, new = function(op, ...) .strat_new(op, list(...))) {
    args <- list()
    for (a in expr$args) if (a$type == expr$type && identical(a$name, expr$name)) args <- c(args, a$args) else args[[length(args) + 1L]] <- a
    do.call(new, c(list(.strat_op(expr)), args))
  },
  rebuild = function(expr) if (!length(expr$args)) expr else sym_rebuild(expr, lapply(expr$args, rl$rebuild))
)

#' @rdname strategies
#' @export
rm_id <- rl$rm_id
#' @rdname strategies
#' @export
unpack <- rl$unpack
#' @rdname strategies
#' @export
glom <- rl$glom
#' @rdname strategies
#' @export
distribute <- rl$distribute
#' @rdname strategies
#' @export
rebuild <- rl$rebuild
#' @rdname strategies
#' @export
strat_new <- function(op, ...) .strat_new(op, list(...))

.sall <- function(rule, fns = basic_fns) function(expr) {
  if (fns$leaf(expr)) return(expr)
  do.call(fns$new, c(list(fns$op(expr)), lapply(fns$children(expr), rule)))
}
.top_down <- function(rule, fns = basic_fns) chain(rule, function(expr) .sall(.top_down(rule, fns), fns)(expr))
.bottom_up <- function(rule, fns = basic_fns) chain(function(expr) .sall(.bottom_up(rule, fns), fns)(expr), rule)

#' @rdname strategies
#' @export
traverse <- list(
  top_down = .top_down,
  bottom_up = .bottom_up,
  top_down_once = function(rule, fns = basic_fns) do_one(rule, function(expr) .sall(.top_down(rule, fns), fns)(expr)),
  bottom_up_once = function(rule, fns = basic_fns) do_one(function(expr) .sall(.bottom_up(rule, fns), fns)(expr), rule),
  sall = .sall)

#' @rdname strategies
#' @export
canon <- function(..., fns = basic_fns) exhaust(.top_down(exhaust(do_one(...)), fns))

#' @rdname strategies
#' @export
typed <- function(ruletypes) strat_switch(function(e) e$type, ruletypes)

# ---------------------------------------------------------------------------
# Branching rules (strategies/branch): functions returning lists

.uniq_add <- function(seen, x) {
  k <- if (is_sym(x)) x$key else paste(deparse(x), collapse = "")
  if (k %in% seen$k) return(FALSE)
  seen$k <- c(seen$k, k); TRUE
}

.b_identity <- function(x) list(x)
.b_chain <- function(...) {
  brules <- list(...)
  function(expr) {
    if (!length(brules)) return(list(expr))
    out <- list()
    for (ne in brules[[1L]](expr)) out <- c(out, do.call(.b_chain, brules[-1L])(ne))
    out
  }
}
.b_multiplex <- function(...) {
  brules <- list(...)
  function(expr) {
    seen <- new.env(); seen$k <- character(); out <- list()
    for (b in brules) for (ne in b(expr)) if (.uniq_add(seen, ne)) out[[length(out) + 1L]] <- ne
    out
  }
}
.b_do_one <- function(...) {
  brules <- list(...)
  function(expr) {
    for (b in brules) { r <- b(expr); if (length(r)) return(r) }
    list()
  }
}
.b_exhaust <- function(brule) {
  rec <- function(expr, seen) {
    out <- list(); any_new <- FALSE
    for (ne in brule(expr)) {
      k <- if (is_sym(ne)) ne$key else paste(deparse(ne), collapse = "")
      if (!(k %in% seen$k)) { seen$k <- c(seen$k, k); any_new <- TRUE; out <- c(out, rec(ne, seen)) }
    }
    if (!any_new) out <- c(out, list(expr))
    out
  }
  function(expr) { seen <- new.env(); seen$k <- if (is_sym(expr)) expr$key else paste(deparse(expr), collapse = ""); rec(expr, seen) }
}
.b_sall <- function(brule, fns = basic_fns) function(expr) {
  if (fns$leaf(expr)) return(list(expr))
  opts <- lapply(fns$children(expr), brule)
  combos <- expand.grid(lapply(opts, seq_along), KEEP.OUT.ATTRS = FALSE)
  # itertools.product order: the last argument varies fastest
  combos <- combos[do.call(order, as.list(combos)), , drop = FALSE]
  lapply(seq_len(nrow(combos)), function(r) do.call(fns$new, c(list(fns$op(expr)), lapply(seq_along(opts), function(j) opts[[j]][[combos[r, j]]]))))
}
.b_top_down <- function(brule, fns = basic_fns) .b_chain(.b_do_one(brule, .b_identity), function(expr) .b_sall(.b_top_down(brule, fns), fns)(expr))

#' @rdname strategies
#' @export
branch <- list(
  identity = .b_identity,
  exhaust = .b_exhaust,
  onaction = function(brule, fn) function(expr) { rs <- brule(expr); for (r in rs) if (!.strat_eq(r, expr)) fn(brule, expr, r); rs },
  debug = function(brule, file = NULL) function(expr) { rs <- brule(expr); for (r in rs) if (!.strat_eq(r, expr)) cat("In: ", format(expr), "\nOut: ", format(r), "\n\n", sep = ""); rs },
  multiplex = .b_multiplex,
  condition = function(cond, brule) function(expr) if (isTRUE(cond(expr))) brule(expr) else list(),
  sfilter = function(pred, brule) function(expr) Filter(function(r) isTRUE(pred(r)), brule(expr)),
  notempty = function(brule) function(expr) { rs <- brule(expr); if (length(rs)) rs else list(expr) },
  do_one = .b_do_one,
  chain = .b_chain,
  yieldify = function(rule) function(expr) list(rule(expr)),
  canon = function(...) .b_exhaust(do.call(.b_multiplex, lapply(list(...), .b_top_down))),
  top_down = .b_top_down,
  sall = .b_sall)

# ---------------------------------------------------------------------------
# Trees of rules (strategies/tree.py)

#' @rdname strategies
#' @export
treeapply <- function(tree, join, leaf = strat_identity) {
  if (inherits(tree, "strat_seq")) return(do.call(join$seq, lapply(tree, treeapply, join = join, leaf = leaf)))
  if (is.list(tree) && !is_sym(tree)) return(do.call(join$list, lapply(tree, treeapply, join = join, leaf = leaf)))
  leaf(tree)
}

#' @rdname strategies
#' @export
greedy <- function(tree, objective = strat_identity, leaf = strat_identity)
  treeapply(tree, list(list = function(...) minimize(..., objective = objective), seq = chain), leaf = leaf)

#' @rdname strategies
#' @export
allresults <- function(tree, leaf = branch$yieldify)
  treeapply(tree, list(list = .b_multiplex, seq = .b_chain), leaf = leaf)

#' @rdname strategies
#' @export
brute <- function(tree, objective = strat_identity, leaf = branch$yieldify) {
  force(tree); force(objective)
  function(expr) {
    rs <- allresults(tree, leaf = leaf)(expr)
    best <- 1L; bv <- objective(rs[[1L]])
    for (i in seq_along(rs)[-1L]) { v <- objective(rs[[i]]); if (.obj_lt(v, bv)) { best <- i; bv <- v } }
    rs[[best]]
  }
}
