# Tree traversal and substitution, after sympy/core/basic.py.


# Types whose some arguments bind variables. Each hook returns the list of
# bound variable keys for an expression of that type.

.free_sym_list <- function(e) {
  t <- e$type
  if (t == "Symbol" || t == "Dummy") return(stats::setNames(list(e), e$key))
  if (!length(e$args)) return(list())
  cacheable <- nchar(e$key) < 9000L
  if (cacheable) {
    hit <- .fs_cache[[e$key]]
    if (!is.null(hit)) return(hit)
  }
  out <- list()
  add <- function(z) {
    s <- .free_sym_list(z)
    if (length(s)) out[names(s)] <<- s
  }
  if (t %in% c("Integral", "Sum", "Product")) {
    # A limit variable binds the integrand and earlier limits only; its
    # endpoints stay free, as ExprWithLimits.free_symbols in SymPy.
    # A composite variable such as f(x) is masked by a Dummy first, so
    # the symbols inside it do not leak out of the integrand.
    f <- e$args[[1L]]
    masks <- character()
    for (tp in e$args[-1L]) {
      v <- tp$args[[1L]]
      if (v$type %in% c("Symbol", "Dummy")) next
      d <- Dummy()
      f <- xreplace(f, list(list(v, d)))
      masks <- c(masks, d$key)
    }
    add(f)
    out <- out[setdiff(names(out), masks)]
    for (tp in e$args[-1L]) {
      v <- tp$args[[1L]]
      if (length(tp$args) == 1L) {
        # an indefinite composite variable such as f(x) counts itself
        if (v$type %in% c("Symbol", "Dummy")) add(v) else out[[v$key]] <- v
        next
      }
      out[[v$key]] <- NULL
      for (z in tp$args[-1L]) add(z)
    }
  } else if (t %in% names(.TRANSFORMS)) {
    # the integration variable binds the function only
    fs <- .free_sym_list(e$args[[1L]])
    fs[[e$args[[2L]]$key]] <- NULL
    out <- fs
    add(e$args[[3L]])
  } else if (t %in% c("Limit", "Subs")) {
    add(e$args[[1L]])
    bv <- e$args[[2L]]
    for (v in if (bv$type == "Tuple") bv$args else list(bv)) out[[v$key]] <- NULL
    for (z in e$args[-(1:2)]) add(z)
  } else {
    for (a in e$args) add(a)
    bh <- .bound_hooks[[t]]
    if (!is.null(bh)) {
      bound <- bh(e)
      if (length(bound)) out <- out[setdiff(names(out), bound)]
    }
  }
  if (cacheable) {
    if (length(.fs_cache) > 50000L) rm(list = ls(.fs_cache, all.names = TRUE), envir = .fs_cache)
    assign(e$key, out, envir = .fs_cache)
  }
  out
}

#' Free symbols of an expression
#'
#' Returns the symbols that appear unbound in `e`, sorted by name. Variables
#' of integration, summation and lambda arguments are bound and not
#' included.
#'
#' @param e A `sym` expression.
#' @return A list of symbols.
#' @examples
#' x <- Symbol("x"); y <- Symbol("y")
#' free_symbols(x^2 + y)
#' @export
free_symbols <- function(e) {
  if (inherits(e, "symmatrix")) {
    out <- list()
    for (a in e$data) out <- c(out, free_symbols(a))
    out <- out[!duplicated(vapply(out, .k, ""))]
    return(.sort_syms(out))
  }
  e <- as_sym(e)
  .sort_syms(unname(.free_sym_list(e)))
}

.sort_syms <- function(lst) {
  if (length(lst) < 2L) return(lst)
  nm <- vapply(lst, function(s) s$name, "")
  lst[order(nm, vapply(lst, .k, ""), method = "radix")]
}

#' Atoms of an expression
#'
#' Collects the leaves of the expression tree. With `types`, only atoms whose
#' type is listed are returned (for example `"Symbol"`, `"Integer"`,
#' `"Rational"`, `"Constant"`, or `"Function"` to collect function
#' applications).
#'
#' @param e A `sym` expression.
#' @param types Optional character vector of types.
#' @return A list of unique `sym` objects.
#' @export
atoms <- function(e, types = NULL) {
  e <- as_sym(e)
  out <- list()
  walk <- function(x) {
    leaf <- !length(x$args)
    if ((is.null(types) && leaf) || (!is.null(types) && x$type %in% types))
      out[[x$key]] <<- x
    for (a in x$args) walk(a)
  }
  walk(e)
  unname(out)
}

#' Preorder and postorder traversal
#'
#' @param e A `sym` expression.
#' @return A list of all subexpressions.
#' @export
preorder_traversal <- function(e) {
  e <- as_sym(e)
  out <- list()
  walk <- function(x) {
    out[[length(out) + 1L]] <<- x
    for (a in x$args) walk(a)
  }
  walk(e)
  out
}

#' @rdname preorder_traversal
#' @export
postorder_traversal <- function(e) {
  e <- as_sym(e)
  out <- list()
  walk <- function(x) {
    for (a in x$args) walk(a)
    out[[length(out) + 1L]] <<- x
  }
  walk(e)
  out
}

#' Test whether an expression contains a subexpression
#'
#' @param e A `sym` expression.
#' @param ... Subexpressions to look for. A symbol matches when it is free in
#'   `e`.
#' @return `TRUE` or `FALSE`.
#' @export
has <- function(e, ...) {
  e <- as_sym(e)
  targets <- lapply(list(...), as_sym)
  for (tg in targets) if (.has1(e, tg)) return(TRUE)
  FALSE
}

.has1 <- function(e, tg) {
  if (identical(e$key, tg$key)) return(TRUE)
  if (tg$type == "Symbol" || tg$type == "Dummy") {
    return(!is.null(.free_sym_list(e)[[tg$key]]))
  }
  for (a in e$args) if (.has1(a, tg)) return(TRUE)
  FALSE
}

# Fast containment check by key without the free symbol shortcut.
.contains_key <- function(e, k) {
  if (identical(e$key, k)) return(TRUE)
  for (a in e$args) if (.contains_key(a, k)) return(TRUE)
  FALSE
}

.has_type <- function(e, types) {
  if (e$type %in% types) return(TRUE)
  for (a in e$args) if (.has_type(a, types)) return(TRUE)
  FALSE
}

.has_function <- function(e, names) {
  if (e$type == "Function" && e$name %in% names) return(TRUE)
  for (a in e$args) if (.has_function(a, names)) return(TRUE)
  FALSE
}

# ---------------------------------------------------------------------------
# Substitution

.subs_pairs <- function(...) {
  args <- list(...)
  plain <- function(z) !is.list(z) || inherits(z, "sym")
  if (length(args) == 2L && plain(args[[1L]]) && plain(args[[2L]]))
    return(list(list(as_sym(args[[1L]]), as_sym(args[[2L]]))))
  if (length(args) == 1L) args <- args[[1L]]
  if (inherits(args, "sym")) stop("subs() needs old and new values")
  out <- list()
  nms <- names(args)
  for (i in seq_along(args)) {
    a <- args[[i]]
    if (!is.null(nms) && nzchar(nms[i])) {
      out[[length(out) + 1L]] <- list(nms[i], as_sym(a))
    } else if (is.list(a) && !inherits(a, "sym") && length(a) == 2L) {
      out[[length(out) + 1L]] <- list(as_sym(a[[1L]]), as_sym(a[[2L]]))
    } else if (inherits(a, "sym") && is_relational(a) && a$type == "Equality") {
      out[[length(out) + 1L]] <- list(a$args[[1L]], a$args[[2L]])
    } else {
      stop("could not interpret substitution argument ", i)
    }
  }
  out
}

# Resolve a name given in subs(list(x = 1)) to the matching free symbol.
.resolve_name <- function(e, nm) {
  fs <- .free_sym_list(e)
  for (s in fs) if (s$name == nm) return(s)
  Symbol(nm)
}

#' Substitute values or expressions
#'
#' `subs()` replaces subexpressions and re-evaluates the result, like SymPy's
#' `expr.subs()`. Substitutions can be given as `subs(e, old, new)`, as a
#' named list `subs(e, list(x = 1, y = 2))` (names refer to free symbols),
#' as a list of pairs `list(list(old, new), ...)`, or as equations. Pairs are
#' applied in sequence unless `simultaneous = TRUE`. `xreplace()` performs
#' exact structural replacement without algebraic matching.
#'
#' @param e A `sym` expression.
#' @param ... Substitutions.
#' @param simultaneous Apply all substitutions at once.
#' @return A `sym` expression.
#' @examples
#' x <- Symbol("x"); y <- Symbol("y")
#' subs(x^2 + y, x, 2)
#' subs(x^2 + y, list(x = 1, y = 3))
#' subs(x * y + 1, x * y, Symbol("z"))
#' @export
subs <- function(e, ..., simultaneous = FALSE) {
  if (inherits(e, "symmatrix")) {
    a <- list(...)
    return(.mat_map(e, function(z) do.call(subs, c(list(z), a, list(simultaneous = simultaneous)))))
  }
  if (is.list(e) && !inherits(e, "sym")) return(lapply(e, function(z) subs(z, ..., simultaneous = simultaneous)))
  e <- as_sym(e)
  pairs <- .subs_pairs(...)
  pairs <- lapply(pairs, function(p) {
    if (is.character(p[[1L]])) list(.resolve_name(e, p[[1L]]), p[[2L]]) else p
  })
  old_mode <- .subs_state$structural
  .subs_state$structural <- .has_type(e, names(.TRANSFORMS))
  on.exit(.subs_state$structural <- old_mode)
  if (simultaneous && length(pairs) > 1L) {
    dums <- lapply(seq_along(pairs), function(i) Dummy(paste0("s", i)))
    for (i in seq_along(pairs)) e <- .subs1(e, pairs[[i]][[1L]], dums[[i]])
    return(xreplace(e, lapply(seq_along(pairs), function(i) list(dums[[i]], pairs[[i]][[2L]]))))
  }
  for (p in pairs) e <- .subs1(e, p[[1L]], p[[2L]])
  e
}


#' @rdname subs
#' @param rule A named list or list of pairs mapping old to new expressions.
#' @export
xreplace <- function(e, rule) {
  e <- as_sym(e)
  pairs <- .subs_pairs(rule)
  pairs <- lapply(pairs, function(p) {
    if (is.character(p[[1L]])) list(.resolve_name(e, p[[1L]]), p[[2L]]) else p
  })
  keys <- vapply(pairs, function(p) p[[1L]]$key, "")
  vals <- lapply(pairs, `[[`, 2L)
  walk <- function(x) {
    i <- match(x$key, keys)
    if (!is.na(i)) return(vals[[i]])
    if (!length(x$args)) return(x)
    args <- lapply(x$args, walk)
    if (identical(vapply(args, .k, ""), vapply(x$args, .k, ""))) return(x)
    sym_rebuild(x, args)
  }
  walk(e)
}

# Apply f bottom up to every node.
.map_bottom_up <- function(e, f) {
  if (length(e$args)) {
    args <- lapply(e$args, .map_bottom_up, f)
    if (!identical(vapply(args, .k, ""), vapply(e$args, .k, ""))) e <- sym_rebuild(e, args)
  }
  f(e)
}

#' Count operations in an expression
#'
#' The number of operations in an expression, as SymPy's `count_ops()`: a
#' subtraction, a division and a negation count once each, and every named
#' function or power counts once. With `visual = TRUE` the result is a sum
#' of operation names such as `2*ADD + SIN`.
#'
#' @param e A `sym` expression.
#' @param visual Return the operations as a symbolic sum.
#' @return An integer, or an expression when `visual = TRUE`.
#' @examples
#' x <- Symbol("x"); y <- Symbol("y")
#' count_ops(x^2 + sin(x)/y)
#' count_ops(x^2 + sin(x)/y, visual = TRUE)
#' @export
count_ops <- function(e, visual = FALSE) {
  e <- as_sym(e)
  ops <- .cops_tab(e)
  if (!visual) return(as.integer(sum(ops)))
  ops <- ops[ops != 0]
  if (!length(ops)) return(S_Zero)
  .add(Map(function(n, k) .mul(list(.int(k), Symbol(n))), names(ops), ops))
}

.cops <- function(e) as.integer(sum(.cops_tab(e)))

.REL_OPNAMES <- c(Equality = "EQ", Unequality = "NE", StrictLessThan = "LT", LessThan = "LE",
                  StrictGreaterThan = "GT", GreaterThan = "GE")

# SymPy's count_ops, as a table of operation counts.
.cops_tab <- function(expr) {
  ops <- integer()
  add_op <- function(nm, k = 1L) ops[nm] <<- (if (is.na(ops[nm])) 0L else ops[nm]) + as.integer(k)
  rel <- expr$type %in% names(.REL_OPNAMES) || expr$type %in% c("And", "Or", "Not", "Xor", "Implies", "Equivalent", "ITE", "Nand", "Nor")
  if (rel) {
    for (a in expr$args) { t <- .cops_tab(a); for (k in names(t)) add_op(k, t[[k]]) }
    add_op(if (!is.null(.REL_OPNAMES[expr$type]) && !is.na(.REL_OPNAMES[expr$type])) .REL_OPNAMES[[expr$type]] else toupper(expr$type))
    return(ops)
  }
  negc <- function(a) (is_number_type(a) && num_sign(a) < 0) || (a$type == "Mul" && is_number_type(a$args[[1L]]) && num_sign(a$args[[1L]]) < 0)
  args <- list(expr)
  while (length(args)) {
    a <- args[[length(args)]]; args[[length(args)]] <- NULL
    t <- a$type
    if (t %in% c("Integer", "Rational")) {
      if (!is_one_num(a)) {
        if (a$value < 0) add_op("NEG")
        if (t == "Rational") add_op("DIV")
        next
      }
    } else if (t == "Mul") {
      if (negc(a)) {
        add_op("NEG")
        a <- if (identical(a$args[[1L]]$key, S_NegativeOne$key)) .mul(a$args[-1L]) else .neg(a)
      }
      nd <- fraction(a, exact = TRUE)
      n <- nd[[1L]]; d <- nd[[2L]]
      if (n$type == "Integer") {
        add_op("DIV")
        if (n$value < 0) add_op("NEG")
        args[[length(args) + 1L]] <- d
        next
      } else if (!is_one_num(d)) {
        if (d$type != "Integer") args[[length(args) + 1L]] <- d
        add_op("DIV")
        args[[length(args) + 1L]] <- n
        next
      }
      t <- a$type
    } else if (t == "Add") {
      negs <- 0L
      for (i in seq_along(a$args)) {
        ai <- a$args[[i]]
        if (negc(ai)) { negs <- negs + 1L; args[[length(args) + 1L]] <- .neg(ai); if (i > 1L) add_op("SUB") }
        else { args[[length(args) + 1L]] <- ai; if (i > 1L) add_op("ADD") }
      }
      if (negs == length(a$args)) add_op("NEG")
      else if (negc(a$args[[1L]])) { add_op("SUB"); add_op("ADD", -1L) }
      next
    }
    if (t == "Pow" && identical(a$args[[2L]]$key, S_NegativeOne$key)) { add_op("DIV"); args[[length(args) + 1L]] <- a$args[[1L]]; next }
    if (is_const(a, "E")) { add_op("EXP"); next }
    if (t == "Mul") add_op("MUL", length(a$args) - 1L)
    else if (t == "Function" && a$name %in% c("Min", "Max")) add_op(toupper(a$name), length(a$args) - 1L)
    else if (length(a$args) && t %in% c("Pow", "Function", "Derivative", "Integral", "Sum", "Product", "Piecewise")) {
      nm <- if (t == "Function") { if (is.null(.FUNS[[a$name]])) paste0("FUNC_", toupper(a$name)) else toupper(a$name) } else toupper(t)
      add_op(nm)
    }
    if (!(t %in% c("Symbol", "Dummy"))) args <- c(args, a$args)
  }
  ops
}
