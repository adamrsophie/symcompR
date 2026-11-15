# Pattern matching of tensor expressions, after sympy/tensor/tensor.py:
# WildTensorIndex, WildTensorHead, WildTensor, _WildTensExpr and the matches
# methods of Tensor, TensMul and TensAdd, with Basic.replace for tensor
# expressions. A match is the list used by .matches() in pattern.R: entries
# list(w = pattern object, v = value), named by the key of the pattern
# object.

# ---------------------------------------------------------------------------
# WildTensorIndex

.wti_new <- function(name, typ, is_up, ignore_updown) {
  structure(list(name = name, tensor_index_type = typ, is_up = is_up, ignore_updown = ignore_updown,
                 key = paste0("WildTensorIndex(", name, ",", .subset2(typ, "key"), ",", ignore_updown, ",", is_up, ")")),
            class = c("WildTensorIndex", "TensorIndex"))
}

#' Wild tensors and tensor pattern matching
#'
#' Ports of SymPy's `WildTensorIndex`, `WildTensorHead` and `WildTensor`.
#' A `WildTensorIndex` matches any index of its index type, of the same
#' variance unless `ignore_updown = TRUE`. A `WildTensorHead` is called like
#' a tensor head; the result, a `WildTensor`, matches any tensor expression
#' whose free indices match its indices (in any order when
#' `unordered_indices = TRUE`), and with no indices it matches any tensor
#' expression. [sym_match()] matches a tensor expression against a pattern
#' and [sym_replace()] replaces the matches, as `expr.match()` and
#' `expr.replace()` do in SymPy; a matched wild tensor head stands for the
#' matched expression with its free indices renamed (SymPy's
#' `_WildTensExpr`).
#'
#' @param name The name of the index or head.
#' @param tensor_index_type A `TensorIndexType`.
#' @param is_up Whether the index is contravariant.
#' @param ignore_updown Whether the index matches indices of either variance.
#' @param index_types,symmetry,comm As for [TensorHead()]; the index types
#'   of a wild tensor are those of its indices, and wild heads have no
#'   symmetry.
#' @param unordered_indices Whether the indices of the wild tensor match in
#'   any order.
#' @return A wild index or a wild tensor head.
#' @examples
#' R3 <- TensorIndexType("R3", dim = 3)
#' idx <- tensor_indices("p q", R3)
#' p <- idx$p; q <- idx$q
#' a <- WildTensorIndex("a", R3, ignore_updown = TRUE)
#' K <- TensorHead("K", list(R3)); V <- TensorHead("V", list(R3))
#' W <- WildTensorHead("W", unordered_indices = TRUE)
#' sym_match(K(p) * V(q), W(a, p))
#' sym_replace(K(p) * V(-p), W(a) * V(-a), 1)
#' sym_replace(K(p) * V(q), W(q) * K(p), W(p) * W(q))
#' @name wild-tensors
#' @export
WildTensorIndex <- function(name, tensor_index_type, is_up = TRUE, ignore_updown = FALSE) {
  if (isTRUE(name)) {
    k <- .subset2(tensor_index_type, "key")
    n <- get0(k, envir = .ti_autogen, inherits = FALSE, ifnotfound = 0L)
    assign(k, n + 1L, envir = .ti_autogen)
    name <- paste0("_i", n)
  } else if (inherits(name, "sym") && name$type == "Symbol") {
    name <- name$name
  } else if (!is.character(name)) {
    .tn_error("ValueError", "invalid name")
  }
  if (inherits(is_up, "sym")) is_up <- identical(is_up$key, S_true$key)
  if (inherits(ignore_updown, "sym")) ignore_updown <- identical(ignore_updown$key, S_true$key)
  .wti_new(name, tensor_index_type, isTRUE(as.logical(is_up)), isTRUE(as.logical(ignore_updown)))
}

# WildTensorIndex.matches
.wti_matches <- function(w, e, d) {
  if (!inherits(e, "TensorIndex")) return(NULL)
  if (!identical(w$tensor_index_type$key, e$tensor_index_type$key)) return(NULL)
  if (!isTRUE(w$ignore_updown) && !identical(w$is_up, e$is_up)) return(NULL)
  d[[w$key]] <- list(w = w, v = e)
  d
}

# ---------------------------------------------------------------------------
# WildTensorHead and WildTensor

.wth_new <- function(name, index_types, comm, unordered) {
  d <- .th_new(name, index_types, TensorSymmetry.no_symmetry(length(index_types)), comm)
  d$key <- paste0("Wild", d$key, ";", unordered)
  d$wild <- TRUE
  d$unordered_indices <- unordered
  d
}

#' @rdname wild-tensors
#' @export
WildTensorHead <- function(name, index_types = NULL, symmetry = NULL, comm = 0, unordered_indices = FALSE) {
  if (inherits(name, "sym") && name$type == "Symbol") name <- name$name
  else if (!is.character(name)) .tn_error("ValueError", "invalid name")
  if (is.null(index_types)) index_types <- list()
  if (inherits(index_types, "TensorIndexType")) index_types <- list(index_types)
  if (!is.null(symmetry)) {
    if (symmetry$rank != length(index_types)) .tn_error("AssertionError", "")
    if (!identical(symmetry$key, TensorSymmetry.no_symmetry(length(index_types))$key))
      .tn_error("NotImplementedError", "Wild matching based on symmetry is not implemented.")
  }
  .th_callable(.wth_new(name, index_types, as_sym(comm), isTRUE(unordered_indices)))
}

# WildTensor.__new__: the index types come from the indices
.wt_new <- function(d, indices) {
  for (i in indices) if (!inherits(i, "TensorIndex")) .tn_error("TypeError", paste0("wrong type for index: ", .tn_py_str(i)))
  d <- .wth_new(d$name, lapply(indices, `[[`, "tensor_index_type"), d$comm, d$unordered_indices)
  structure(list(head = d, indices = indices, is = .tis_from_indices(indices), is_canon_bp = FALSE,
                 key = paste0("Tensor(", d$key, ";", paste(.tn_keys(indices), collapse = ","), ")")),
            class = c("WildTensor", "Tensor", "TensExpr"))
}

# _WildTensExpr: a matched expression, called with new free indices
.wte_new <- function(e) structure(list(expr = e, key = paste0("_WildTensExpr(", .tn_key(e), ")")), class = "WildTensExpr")
.wte_call <- function(w, indices) {
  fi <- .tn_get_free_indices(w$expr)
  repl <- list()
  for (k in seq_along(fi)) repl[[fi[[k]]$key]] <- indices[[k]]
  .tn_replace_indices(w$expr, repl)
}
.wte_add <- function(a, b) .wte_new(a$expr + b$expr)
#' @export
format.WildTensExpr <- function(x, ...) paste0("_WildTensExpr(", .tn_str(x$expr), ")")
#' @export
print.WildTensExpr <- function(x, ...) { cat(format(x), "\n"); invisible(x) }

# _get_wilds: the wild objects in a tensor expression
.tn_wilds <- function(x) {
  if (inherits(x, "WildTensor")) return(c(list(x), Filter(function(i) inherits(i, "WildTensorIndex"), x$indices)))
  if (inherits(x, "Tensor")) return(Filter(function(i) inherits(i, "WildTensorIndex"), x$indices))
  if (inherits(x, "TensMul") || inherits(x, "TensAdd")) return(do.call("c", c(list(list()), lapply(x$args, .tn_wilds))))
  if (inherits(x, "sym") && .has_type(x, "Wild")) return(list(x))
  list()
}
.tn_has_wildtensor <- function(x) any(vapply(.tn_wilds(x), function(w) inherits(w, "WildTensor"), NA))

.tn_merge <- function(d, m) { for (k in names(m)) d[[k]] <- m[[k]]; d }

# ---------------------------------------------------------------------------
# matches

.tn_matches <- function(p, e, d = list(), old = FALSE) {
  if (is.null(d)) d <- list()
  if (inherits(p, "WildTensorIndex")) return(.wti_matches(p, e, d))
  if (inherits(p, "TensorIndex")) return(if (inherits(e, "TensorIndex") && identical(p$key, e$key)) d else NULL)
  if (inherits(p, "WildTensor")) return(.wt_matches(p, e, d))
  if (inherits(p, "Tensor")) return(.tensor_matches(p, e, d))
  if (inherits(p, "TensMul")) return(.tensmul_matches(p, e, d, old))
  if (inherits(p, "TensAdd")) return(.tensadd_matches(p, e, d, old))
  if (.tn_is(e) || inherits(e, "TensorIndex")) return(NULL)
  .matches(as_sym(p), as_sym(e), d, old)
}

# TensExpr._matches_simple: no wilds, so equal up to canonicalization
.tn_matches_simple <- function(p, e, d) {
  pf <- sort(.tn_keys(.tn_get_free_indices(p)))
  if (!.tn_is(e)) {
    if (length(pf)) return(NULL)
  } else if (!identical(unique(pf), unique(sort(.tn_keys(.tn_get_free_indices(e)))))) {
    return(NULL)
  }
  diff <- .tn_canon_bp(.tn_add(p, .tn_neg(e)))
  if (!.tn_is(diff) && .tn_is_num(diff, 0)) d else NULL
}

# WildTensor.matches
.wt_matches <- function(w, e, d) {
  if (!.tn_is(e) && !.tn_is_num(e, 1)) return(NULL)
  if (length(w$indices)) {
    if (!.tn_is(e)) return(NULL)
    ei <- .tn_get_free_indices(e)
    if (length(ei) != length(w$indices)) return(NULL)
    if (isTRUE(w$head$unordered_indices)) {
      m <- .wt_match_unordered(w, e)
      if (is.null(m)) return(NULL)
      d <- .tn_merge(d, m)
    } else {
      for (i in seq_along(ei)) {
        m <- .tn_matches(w$indices[[i]], ei[[i]])
        if (is.null(m)) return(NULL)
        d <- .tn_merge(d, m)
      }
    }
    d[[w$head$key]] <- list(w = w$head, v = .wte_new(e))
  } else {
    d[[w$key]] <- list(w = w, v = e)
  }
  d
}

# WildTensor._match_indices_ignoring_order
.wt_match_unordered <- function(w, e) {
  d <- list()
  kind <- vapply(w$indices, function(i) if (!inherits(i, "WildTensorIndex")) "nonwild" else if (isTRUE(i$ignore_updown)) "wild, updown" else "wild", "")
  matched <- list()
  remaining <- .tn_get_indices(e)
  in_matched <- function(i) i$key %in% .tn_keys(matched)
  for (ind in w$indices[kind == "nonwild"]) {
    ok <- FALSE
    for (ei in remaining) {
      if (in_matched(ei)) next
      m <- .tn_matches(ind, ei)
      if (!is.null(m)) { ok <- TRUE; d <- .tn_merge(d, m); matched[[length(matched) + 1L]] <- ei; break }
    }
    if (!ok) return(NULL)
  }
  for (grp in c("wild", "wild, updown")) {
    remaining <- Filter(function(i) !in_matched(i), remaining)
    for (ind in w$indices[kind == grp]) {
      ok <- FALSE
      for (ei in remaining) {
        m <- .tn_matches(ind, ei)
        if (is.null(m)) next
        nk <- .tn_neg_key(ind$key)
        if (!is.null(d[[nk]]) && !identical((-d[[nk]]$v)$key, m[[ind$key]]$v$key)) return(NULL)
        ok <- TRUE; d <- .tn_merge(d, m); matched[[length(matched) + 1L]] <- ei
        break
      }
      if (!ok) return(NULL)
    }
  }
  if (length(matched) < length(w$indices)) return(NULL)
  d
}

# Tensor._get_symmetrized_forms: the tensor under the generators of its
# symmetry, with their signs
.tensor_symmetrized_forms <- function(t) {
  gens <- lapply(t$head$symmetry$generators, .tc_af)
  rank <- .th_rank(t$head)
  forms <- list(t); keys <- .tn_key(t)
  k <- 1L
  while (k <= length(forms)) {
    tens <- forms[[k]]
    for (g in gens) {
      inds <- .tn_get_indices(tens)
      repl <- list()
      for (i in seq_len(rank)) repl[[inds[[i]]$key]] <- inds[[g[i] + 1L]]
      nt <- .tn_replace_indices(tens, repl)
      if (g[rank + 1L] != rank) nt <- .tn_neg(nt)
      if (!(.tn_key(nt) %in% keys)) { forms[[length(forms) + 1L]] <- nt; keys <- c(keys, .tn_key(nt)) }
    }
    k <- k + 1L
  }
  forms
}

.tensor_matches <- function(p, e, d) {
  if (.tn_is(e) && .tn_eq(p, e)) return(d)
  if (!inherits(e, "Tensor") || !identical(p$head$key, e$head$key)) return(NULL)
  for (ne in .tensor_symmetrized_forms(e)) {
    m <- .tensor_matches0(p, ne, d)
    if (!is.null(m)) return(.tn_merge(d, m))
  }
  NULL
}

# Tensor._matches: the indices in order
.tensor_matches0 <- function(p, e, d) {
  if (.tn_is(e) && .tn_eq(p, e)) return(d)
  if (!inherits(e, "Tensor") || !identical(p$head$key, e$head$key)) return(NULL)
  si <- .tn_get_indices(p); ei <- .tn_get_indices(e)
  if (length(si) != length(ei)) return(NULL)
  for (i in seq_along(si)) {
    s <- si[[i]]
    m <- .tn_matches(s, ei[[i]])
    if (is.null(m)) return(NULL)
    nk <- .tn_neg_key(s$key)
    if (!is.null(d[[nk]]) && !is.null(m[[s$key]]) && !identical((-d[[nk]]$v)$key, m[[s$key]]$v$key)) return(NULL)
    d <- .tn_merge(d, m)
  }
  d
}

.tn_scalar_part <- function(t) {
  r <- S_One
  for (a in t$args) if (!.tn_is(a)) r <- .mul(list(r, a))
  r
}
.tn_drop <- function(args, x) Filter(function(a) !(.tn_is(a) && .tn_eq(a, x)), args)

.tensmul_matches <- function(p, e, d, old) {
  if (inherits(e, "TensMul")) {
    for (a in e$args) if (inherits(a, "Tensor") && .th_comm(a$head) != 0)
      .tn_error("NotImplementedError", "Tensor matching not implemented for non-commuting tensors")
  }
  .tensmul_matches_commutative(p, e, d, old)
}

# TensMul._matches_commutative
.tensmul_matches_commutative <- function(p, e, d, old) {
  if (!inherits(e, "TensMul")) {
    if (.tn_is(e) || inherits(e, "sym")) e <- .tensmul_new(list(e)) else return(NULL)
  }
  if (length(d)) {
    vals <- Filter(function(v) inherits(v, "TensorIndex"), lapply(d, `[[`, "v"))
    if (length(vals)) {
      r <- .tm_dedupe_indices(p, vals); if (!is.null(r)) p <- r
      r <- .tm_dedupe_indices(e, vals); if (!is.null(r)) e <- r
    }
  }
  e_coeff <- .tn_scalar_part(e)
  if (.tn_eq(p, e)) return(d)
  if (!length(.tn_wilds(p))) return(.tn_matches_simple(p, e, d))
  q_wild <- Filter(function(a) inherits(a, "WildTensor"), p$args)
  q_tens <- Filter(function(a) .tn_is(a) && !inherits(a, "WildTensor"), p$args)
  e_tens <- Filter(function(a) .tn_is(a) && !inherits(a, "WildTensor"), e$args)
  heads <- function(x) paste(vapply(.tn_components(x), `[[`, "", "key"), collapse = ";")
  if (!all(vapply(q_tens, heads, "") %in% vapply(e_tens, heads, ""))) return(NULL)
  if (length(q_tens)) {
    q_tensor <- q_tens[[1L]]
    forms <- if (inherits(q_tensor, "Tensor")) .tensor_symmetrized_forms(q_tensor) else list(q_tensor)
    for (qf in forms) {
      sign <- if (inherits(qf, "TensMul")) -1L else 1L
      qs <- if (sign < 0L) .tn_neg(qf) else qf
      for (et in e_tens) {
        m <- if (inherits(qs, "Tensor")) .tensor_matches0(qs, et, list()) else .tn_matches(qs, et)
        if (is.null(m)) next
        rem_query <- .tn_doit(.tensmul_new(c(list(as_sym(sign)), .tn_drop(p$args, q_tensor))), FALSE)
        rem_expr <- .tn_doit(.tensmul_new(.tn_drop(e$args, et)), FALSE)
        rem_m <- .tn_matches(rem_query, rem_expr, .tn_merge(d, m), old)
        if (is.null(rem_m)) next
        ok <- TRUE
        for (k in names(rem_m)) {
          if (!inherits(rem_m[[k]]$w, "TensorIndex")) next
          nk <- .tn_neg_key(k)
          if (!is.null(rem_m[[nk]]) && !identical(rem_m[[nk]]$v$key, (-rem_m[[k]]$v)$key)) { ok <- FALSE; break }
        }
        if (ok) return(.tn_merge(d, rem_m))
      }
    }
    return(NULL)
  }
  matched <- list()
  remaining <- e_tens
  indexless <- Filter(function(w) !length(.tn_get_free_indices(w)), q_wild)
  wilds <- Filter(function(w) length(.tn_get_free_indices(w)) > 0L, q_wild)
  for (w in wilds) {
    fw <- .tn_get_free_indices(w)
    try_t <- Filter(function(t) {
      all(vapply(.tn_get_free_indices(t), function(i) any(vapply(fw, function(j) !is.null(.tn_matches(j, i)), NA)), NA))
    }, remaining)
    m <- .wt_matches(w, .tn_doit(.tensmul_new(try_t), FALSE), list())
    if (is.null(m)) return(NULL)
    matched <- c(matched, try_t)
    d <- .tn_merge(d, m)
  }
  remaining <- Filter(function(t) !(.tn_key(t) %in% .tn_keys(matched)), e_tens)
  if (length(indexless)) {
    m <- .wt_matches(indexless[[1L]], .tn_doit(.tensmul_new(c(list(S_One), remaining)), FALSE), list())
    if (is.null(m)) return(NULL)
    d <- .tn_merge(d, m)
  } else if (length(remaining)) {
    return(NULL)
  }
  m <- .matches(.tn_scalar_part(p), e_coeff, list(), old)
  if (is.null(m)) return(NULL)
  .tn_merge(d, m)
}

# TensAdd.matches
.tensadd_matches <- function(p, e, d, old) {
  if (!inherits(e, "TensAdd")) return(NULL)
  if (!length(.tn_wilds(p))) return(.tn_matches_simple(p, e, d))
  sift <- function(a) {
    ws <- .tn_wilds(a)
    if (!length(ws)) return("nonwild")
    wt <- Filter(function(w) inherits(w, "WildTensor"), ws)
    if (length(wt)) return(if (any(vapply(wt, function(w) !length(w$indices), NA))) "indexless_wildtensor" else "wildtensor")
    "otherwild"
  }
  qk <- vapply(p$args, sift, ""); ek <- vapply(e$args, sift, "")
  e_non <- e$args[ek == "nonwild"]
  matched <- character()
  for (qt in p$args[qk == "nonwild"]) {
    ok <- FALSE
    for (et in e_non) {
      if (.tn_key(et) %in% matched) next
      m <- .tn_matches(qt, et, d, old)
      if (is.null(m)) next
      ok <- TRUE; d <- .tn_merge(d, m); matched <- c(matched, .tn_key(et))
      break
    }
    if (!ok) return(NULL)
  }
  for (grp in c("otherwild", "wildtensor", "indexless_wildtensor")) {
    remaining <- Filter(function(t) !(.tn_key(t) %in% matched), e_non)
    for (w in p$args[qk == grp]) {
      wk <- if (inherits(w, "WildTensor") && length(w$indices)) w$head$key else .tn_key(w)
      for (et in remaining) {
        m <- .tn_matches(w, et)
        if (is.null(m)) next
        matched <- c(matched, .tn_key(et))
        if (!is.null(d[[wk]]) && !is.null(m[[wk]])) {
          a <- d[[wk]]$v; b <- m[[wk]]$v
          d[[wk]]$v <- if (inherits(a, "WildTensExpr")) .wte_add(a, b) else a + b
          m[[wk]] <- NULL
        }
        d <- .tn_merge(d, m)
      }
    }
  }
  if (any(!(vapply(e_non, .tn_key, "") %in% matched))) return(NULL)
  d
}

# ---------------------------------------------------------------------------
# replace

# value.subs(result) for a tensor or scalar value: wild indices, wild heads
# and indexless wild tensors take their matched values
.tn_wsubs <- function(x, d) {
  if (inherits(x, "TensorIndex")) return(if (!is.null(d[[x$key]])) d[[x$key]]$v else x)
  if (inherits(x, "WildTensor")) {
    if (!length(x$indices)) return(if (!is.null(d[[x$key]])) d[[x$key]]$v else x)
    inds <- lapply(x$indices, .tn_wsubs, d = d)
    hv <- d[[x$head$key]]$v
    if (inherits(hv, "WildTensExpr")) return(.wte_call(hv, inds))
    return(.wt_new(x$head, inds))
  }
  if (inherits(x, "Tensor")) {
    inds <- lapply(x$indices, .tn_wsubs, d = d)
    if (identical(.tn_keys(inds), .tn_keys(x$indices))) return(x)
    return(.tensor_doit(.tensor_new(x$head, inds)))
  }
  if (inherits(x, "TensMul")) {
    # TensMul._eval_subs: first rename dummies that clash with new indices
    new <- Filter(function(v) inherits(v, "TensorIndex"), lapply(d, `[[`, "v"))
    if (length(new)) { r <- .tm_dedupe_indices(x, new); if (!is.null(r)) x <- r }
    return(.tn_doit(.tensmul_new(lapply(x$args, .tn_wsubs, d = d)), FALSE))
  }
  if (inherits(x, "TensAdd")) return(.tn_doit(.tensadd_new(lapply(x$args, .tn_wsubs, d = d)), FALSE))
  if (inherits(x, "sym")) {
    sd <- Filter(function(z) inherits(z$w, "sym"), d)
    return(.rd_xrep(x, sd))
  }
  x
}

# Basic.replace with a pattern: bottom-up, each subexpression that matches
# the pattern is replaced by the value with the match substituted
.tn_replace <- function(e, pattern, value, exact = NULL) {
  if (!.tn_is(value) && !is.function(value)) value <- as_sym(value)
  if (is.null(exact)) exact <- sum(vapply(.tn_wilds(pattern), function(w) inherits(w, "sym") && w$type == "Wild", NA)) > 1L
  f <- function(x) {
    m <- .tn_matches(pattern, x)
    if (is.null(m)) return(x)
    if (exact && any(vapply(m, function(z) inherits(z$v, "sym") && .tn_is_num(z$v, 0), NA))) return(x)
    if (is.function(value)) value(x, m) else .tn_wsubs(value, m)
  }
  walk <- function(x) {
    if (inherits(x, "TensMul") || inherits(x, "TensAdd")) {
      na <- lapply(x$args, walk)
      if (!identical(vapply(na, .tn_key, ""), vapply(x$args, .tn_key, "")))
        x <- if (inherits(x, "TensMul")) .tn_doit(.tensmul_new(na), FALSE) else .tn_doit(.tensadd_new(na), FALSE)
    }
    f(x)
  }
  walk(e)
}

# expr.match(pattern) for tensor expressions, named by the pattern objects
.tn_match_named <- function(e, pattern) {
  d <- .tn_matches(pattern, e)
  if (is.null(d)) return(NULL)
  out <- list()
  for (z in d) {
    nm <- if (inherits(z$w, "sym")) z$w$name else if (is.function(z$w) || !is.null(z$w$index_types)) format.TensorHead(z$w) else format(z$w)
    out[[nm]] <- z$v
  }
  out
}
