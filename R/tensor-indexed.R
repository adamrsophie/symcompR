# Indexed objects and index analysis, after sympy/tensor/indexed.py and
# sympy/tensor/index_methods.py.

#' Indexed objects
#'
#' `Idx(name, range)` is an integer index symbol, as SymPy's `Idx`; `range`
#' is an upper bound `n` (giving `0` to `n - 1`) or a pair of bounds.
#' `IndexedBase(name, shape)` returns a function building indexed objects
#' `A(i, j)` (SymPy's `A[i, j]`), and `Indexed(base, ...)` builds one
#' directly. `get_indices(e)` returns the free indices of an expression and
#' its symmetries; repeated indices in a product are summed (Einstein
#' convention), and the terms of a sum must have the same free indices.
#' `get_contraction_structure(e)` lists the summations: index tuples with
#' the terms summed over them, `NULL` for terms without summation, and
#' nested structures for subexpressions, as SymPy's dictionary.
#'
#' @param name The name, or a symbol to take it from (an integer symbol for
#'   `Idx()`).
#' @param range An upper bound or a pair of bounds.
#' @param shape Optional dimensions.
#' @param base An `IndexedBase` function or a name.
#' @param ... Indices.
#' @param e An expression.
#' @return `Idx()` and `Indexed()` return expressions; `IndexedBase()` a
#'   function; `get_indices()` a list of the free indices and the
#'   symmetries; `get_contraction_structure()` a list of `key`, `value`
#'   pairs.
#' @examples
#' i <- Idx("i"); j <- Idx("j")
#' A <- IndexedBase("A"); x <- IndexedBase("x"); y <- IndexedBase("y")
#' A(i, j)
#' get_indices(A(i, j) * x(j))$indices
#' format_contraction_structure(get_contraction_structure(x(i)*y(i) + A(j, j)))
#' @name Indexed
NULL

#' @rdname Indexed
#' @export
Idx <- function(name, range = NULL) {
  if (is_sym(name) && name$type == "Symbol") {
    if (!isTRUE(ask_fact(name, "integer"))) .cg_err("TypeError", "Idx object requires an integer label.")
    name <- name$name
  }
  s <- Symbol(name, integer = TRUE)
  lohi <- NULL
  if (!is.null(range)) {
    lohi <- if (length(range) == 2L) lapply(as.list(range), as_sym) else list(S_Zero, .sub(as_sym(range), S_One))
  }
  new_sym("Symbol", name = name, extra = c(s$extra, list(idx = TRUE, range = lohi)),
          key = paste0("Idx:", s$key, if (!is.null(lohi)) paste0("[", lohi[[1L]]$key, ",", lohi[[2L]]$key, "]") else ""))
}
.is_idx <- function(z) z$type == "Symbol" && isTRUE(z$extra$idx)

#' @rdname Indexed
#' @export
IndexedBase <- function(name, shape = NULL) {
  if (is_sym(name) && name$type == "Symbol") name <- name$name
  if (is_sym(shape) && shape$type == "Tuple") shape <- shape$args
  base <- new_sym("IndexedBase", list(Symbol(name)), extra = list(shape = if (is.null(shape)) NULL else lapply(as.list(shape), as_sym)))
  f <- function(...) new_sym("Indexed", c(list(base), lapply(list(...), as_sym)))
  attr(f, "base") <- base
  f
}

#' @rdname Indexed
#' @export
Indexed <- function(base, ...) {
  b <- if (is.function(base)) attr(base, "base") else if (is.character(base)) attr(IndexedBase(base), "base") else as_sym(base)
  new_sym("Indexed", c(list(b), lapply(list(...), as_sym)))
}

register_str_printer("IndexedBase", function(e) .str(e$args[[1L]]))
register_latex_printer("IndexedBase", function(e) .latex(e$args[[1L]]))
register_str_printer("Indexed", function(e) paste0(.str(e$args[[1L]]), "[", paste(vapply(e$args[-1L], .str, ""), collapse = ", "), "]"))
register_latex_printer("Indexed", function(e) paste0("{", .latex(e$args[[1L]]), "}_{", paste(vapply(e$args[-1L], .latex, ""), collapse = ","), "}"))

# Sets of indices are kept as lists without duplicates, sorted for printing.
.iset <- function(l) {
  if (!length(l)) return(list())
  l <- l[!duplicated(vapply(l, .k, ""))]
  .key_sort(lapply(l, .sort_key), l)
}
.iset_eq <- function(a, b) setequal(vapply(a, .k, ""), vapply(b, .k, ""))

# _unique_and_repeated, in order of first appearance
.idx_uniq_rep <- function(inds) {
  ks <- vapply(inds, .k, "")
  first <- !duplicated(ks)
  rep_k <- unique(ks[duplicated(ks)])
  list(u = inds[first & !(ks %in% rep_k)], r = inds[first & ks %in% rep_k])
}

.idx_mul <- function(e) {
  parts <- lapply(e$args, get_indices)
  inds <- do.call(c, lapply(parts, `[[`, "indices"))
  ur <- .idx_uniq_rep(if (is.null(inds)) list() else inds)
  list(indices = .iset(ur$u), symmetries = list(), dummies = ur$r)
}

#' @rdname Indexed
#' @export
get_indices <- function(e) {
  if (inherits(e, "TensExpr")) return(.tn_get_indices(e))
  e <- as_sym(e)
  if (e$type == "Indexed") return(list(indices = .iset(.idx_uniq_rep(e$args[-1L])$u), symmetries = list()))
  if (.is_idx(e)) return(list(indices = list(e), symmetries = list()))
  if (!length(e$args)) return(list(indices = list(), symmetries = list()))
  if (e$type == "Mul") { r <- .idx_mul(e); return(r[c("indices", "symmetries")]) }
  if (e$type == "Add") {
    inds <- lapply(e$args, function(a) get_indices(a)$indices)
    ns <- Filter(length, inds)
    if (!length(ns)) return(list(indices = list(), symmetries = list()))
    if (!all(vapply(ns[-1L], .iset_eq, TRUE, b = ns[[1L]])))
      stop(structure(class = c("IndexConformanceException", "error", "condition"),
                     list(message = paste0("Indices are not consistent: ", .str(e)), call = NULL)))
    return(list(indices = ns[[1L]], symmetries = list()))
  }
  if (e$type == "Pow" || (e$type == "Function" && e$name == "exp")) {
    be <- .base_exp(e)
    return(list(indices = .iset(c(get_indices(be[[1L]])$indices, get_indices(be[[2L]])$indices)), symmetries = list()))
  }
  if (e$type == "Piecewise") return(list(indices = list(), symmetries = list()))
  if (e$type == "Function") return(list(indices = .iset(do.call(c, lapply(e$args, function(a) get_indices(a)$indices))), symmetries = list()))
  if (!.has_type(e, "Indexed")) return(list(indices = list(), symmetries = list()))
  stop("no specialized handling of ", e$type, " in get_indices()", call. = FALSE)
}

# a contraction structure: list of list(key, value); key NULL, a list of
# indices, or an expression; value a list of terms or of nested structures
.cs_put <- function(cs, key, value) {
  kk <- .cs_key(key)
  for (i in seq_along(cs)) if (identical(.cs_key(cs[[i]]$key), kk)) {
    cs[[i]]$value <- .iset(c(cs[[i]]$value, value)); return(cs)
  }
  cs[[length(cs) + 1L]] <- list(key = key, value = value)
  cs
}
.cs_key <- function(k) if (is.null(k)) "None" else if (is_sym(k)) k$key else paste(vapply(k, .k, ""), collapse = ",")
.cs_trivial <- function(d) length(d) == 1L && is.null(d[[1L]]$key)

#' @rdname Indexed
#' @export
get_contraction_structure <- function(e) {
  e <- as_sym(e)
  one <- function(key, val) list(list(key = key, value = val))
  if (e$type == "Indexed") {
    r <- .idx_uniq_rep(e$args[-1L])$r
    return(one(if (length(r)) r else NULL, list(e)))
  }
  if (!length(e$args)) return(one(NULL, list(e)))
  nested_of <- function(subs) Filter(Negate(.cs_trivial), lapply(subs, get_contraction_structure))
  if (e$type == "Mul") {
    d <- .idx_mul(e)$dummies
    res <- one(if (length(d)) d else NULL, list(e))
    nested <- nested_of(e$args)
    if (length(nested)) res[[2L]] <- list(key = e, value = nested)
    return(res)
  }
  if (e$type == "Pow" || (e$type == "Function" && e$name == "exp")) {
    be <- .base_exp(e)
    res <- one(NULL, list(e))
    nested <- nested_of(be)
    if (length(nested)) res[[2L]] <- list(key = e, value = nested)
    return(res)
  }
  if (e$type == "Add") {
    res <- list()
    # terms in SymPy's argument order (class first), as its dict insertion
    terms <- .key_sort(lapply(e$args, .sort_key), e$args)
    for (t in terms) for (p in get_contraction_structure(t)) res <- .cs_put(res, p$key, p$value)
    return(res)
  }
  if (e$type == "Piecewise") return(one(NULL, e))
  if (e$type == "Function") {
    res <- one(NULL, list(e))
    nested <- nested_of(e$args)
    if (length(nested)) res[[2L]] <- list(key = e, value = nested)
    return(res)
  }
  if (!.has_type(e, "Indexed")) return(one(NULL, list(e)))
  stop("no specialized handling of ", e$type, " in get_contraction_structure()", call. = FALSE)
}

#' @rdname Indexed
#' @param cs A contraction structure.
#' @export
format_contraction_structure <- function(cs) {
  fk <- function(k) if (is.null(k)) "None" else if (is_sym(k)) .str(k)
                    else paste0("(", paste(vapply(k, .str, ""), collapse = ", "), if (length(k) == 1L) ",)" else ")")
  fv <- function(p) {
    if (!is.null(p$key) && is_sym(p$key) && !(p$key$type %in% c("Symbol"))) {
      return(paste0("[", paste(vapply(p$value, format_contraction_structure, ""), collapse = ", "), "]"))
    }
    paste0("{", paste(vapply(.iset(p$value), .str, ""), collapse = ", "), "}")
  }
  paste0("{", paste(vapply(cs, function(p) paste0(fk(p$key), ": ", fv(p)), ""), collapse = ", "), "}")
}
