# Abstract tensors with index notation, after sympy/tensor/tensor.py:
# index types, indices, monoterm symmetries, tensor heads, tensors, their
# products (TensMul) and sums (TensAdd), Butler-Portugal canonicalization
# (canon_bp, via R/tensor-can.R), metric and delta contraction, index
# substitution and replacement of tensors by arrays.
#
# Tensor expressions are S3 objects of class "TensExpr" (not sym nodes); their
# coefficients are sym expressions. A TensorHead is a function: calling it
# with indices gives a Tensor.

# ---------------------------------------------------------------------------
# small helpers

.tn_is <- function(x) inherits(x, "TensExpr")
.tn_sym <- function(x) if (.tn_is(x)) x else as_sym(x)
.tn_key <- function(x) if (is.function(x)) attr(x, "th")$key else if (is.list(x)) x$key else as_sym(x)$key
.tn_eq <- function(a, b) identical(.tn_key(a), .tn_key(b))
.tn_keys <- function(l) vapply(l, .tn_key, "")
.tn_error <- function(cls, msg) stop(.sympy_condition(cls, msg))
.tn_is_num <- function(x, v) inherits(x, "sym") && identical(x$key, as_sym(v)$key)
.tn_uniq <- function(l) l[!duplicated(.tn_keys(l))]

# stable insertion sort with a strict `lt` comparator, as Python's sort
.tn_sort <- function(items, lt) {
  o <- integer()
  for (i in seq_along(items)) {
    j <- length(o)
    while (j > 0L && lt(items[[i]], items[[o[j]]])) j <- j - 1L
    o <- append(o, i, after = j)
  }
  items[o]
}

# str() of a Permutation
.tn_perm_str <- function(p) {
  a <- .tc_af(p)
  n <- length(a)
  if (!n) return("()")
  cf <- cyclic_form(.tc_af_new(a))
  s <- paste0(vapply(cf, function(c) paste0("(", paste(c, collapse = " "), ")"), ""), collapse = "")
  if (!((n - 1L) %in% unlist(cf))) s <- paste0("(", n - 1L, ")", s)
  s
}

# ---------------------------------------------------------------------------
# TensorManager: commutation groups of tensors

.TM <- new.env()
.tm_init <- function() {
  .TM$comm <- list(list(`0` = 0, `1` = 0, `2` = 0), list(`0` = 0, `1` = 1, `2` = NA),
                   list(`0` = 0, `1` = NA))
  .TM$s2i <- list(`0` = 0L, `1` = 1L, `2` = 2L)
  .TM$i2s <- list(0, 1, 2)
  invisible(NULL)
}
.tm_init()

.tm_label_key <- function(i) {
  if (!inherits(i, "sym")) i <- as_sym(i)
  if (i$type == "Integer") return(as.character(i$value))
  i$key
}

.tm_add_group <- function(i) {
  n <- length(.TM$comm)
  .TM$comm[[n + 1L]] <- list(`0` = 0)
  .TM$comm[[1L]][[as.character(n)]] <- 0
  .TM$s2i[[.tm_label_key(i)]] <- n
  .TM$i2s[[n + 1L]] <- i
  n
}

.tm_symbols2i <- function(i) {
  k <- .tm_label_key(i)
  if (is.null(.TM$s2i[[k]])) return(.tm_add_group(i))
  .TM$s2i[[k]]
}

.tm_get_comm <- function(i, j) {
  row <- .TM$comm[[i + 1L]]
  v <- row[[as.character(j)]]
  if (is.null(v)) return(if (i == 0L || j == 0L) 0 else NULL)
  if (is.na(v)) NULL else v
}

.tm_set_comm <- function(i, j, c) {
  if (!(is.null(c) || (length(c) == 1L && c %in% c(0, 1))))
    .tn_error("ValueError", "`c` can assume only the values 0, 1 or None")
  i <- as_sym(i)
  j <- as_sym(j)
  if (is.null(.TM$s2i[[.tm_label_key(i)]])) .tm_add_group(i)
  if (is.null(.TM$s2i[[.tm_label_key(j)]])) .tm_add_group(j)
  ni <- .TM$s2i[[.tm_label_key(i)]]
  nj <- .TM$s2i[[.tm_label_key(j)]]
  v <- if (is.null(c)) NA else c
  .TM$comm[[ni + 1L]][[as.character(nj)]] <- v
  .TM$comm[[nj + 1L]][[as.character(ni)]] <- v
  invisible(NULL)
}

#' Commutation groups of tensors
#'
#' Port of SymPy's `TensorManager`. Each tensor head belongs to a
#' commutation group: `0` commutes with everything, `1` anticommutes within
#' itself, `2` commutes only with group `0`. `TensorManager$set_comm(i, j,
#' c)` sets the commutation (`0` commuting, `1` anticommuting, `NULL`
#' neither) of the groups labeled `i` and `j`; `get_comm(i, j)` reads it
#' from group numbers; `comm_symbols2i(i)` gives the group number of a label
#' (creating it); `set_comms(...)` takes `list(i, j, c)` triples; `clear()`
#' restores the defaults.
#'
#' @format An environment of functions.
#' @examples
#' TensorManager$set_comm("Gcomm", "GHcomm", 0)
#' TensorManager$get_comm(TensorManager$comm_symbols2i("Gcomm"),
#'                        TensorManager$comm_symbols2i("GHcomm"))
#' TensorManager$clear()
#' @export
TensorManager <- new.env()
TensorManager$set_comm <- function(i, j, c) .tm_set_comm(i, j, c)
TensorManager$set_comms <- function(...) {
  for (a in list(...)) .tm_set_comm(a[[1L]], a[[2L]], a[[3L]])
  invisible(NULL)
}
TensorManager$get_comm <- function(i, j) .tm_get_comm(i, j)
TensorManager$comm_symbols2i <- function(i) .tm_symbols2i(i)
TensorManager$comm_i2symbol <- function(i) .TM$i2s[[i + 1L]]
TensorManager$clear <- function() .tm_init()

# ---------------------------------------------------------------------------
# TensorIndexType

#' Tensor index types, indices and symmetries
#'
#' Ports of SymPy's `TensorIndexType`, `TensorIndex`, `tensor_indices()` and
#' `TensorSymmetry`. An index type has a name, a dummy name (the head of
#' generated dummy index names such as `L_0`), a dimension and a metric;
#' `L$metric`, `L$delta` and `L$epsilon` are its metric, Kronecker delta
#' and Levi-Civita tensor heads. A `TensorIndex` is contravariant (`is_up`)
#' or covariant; unary minus flips it. `TensorSymmetry(base, gens)` is a
#' monoterm slot symmetry given by a BSGS; `TensorSymmetry.fully_symmetric(n)`
#' (antisymmetric for `n < 0`), `TensorSymmetry.direct_product(...)`,
#' `TensorSymmetry.riemann()` and `TensorSymmetry.no_symmetry(n)` build the
#' usual ones.
#'
#' @param name Name (a string).
#' @param dummy_name Head of dummy index names (default: first letter).
#' @param dim Dimension (default: symbol `dim_<dummy_name>`).
#' @param eps_dim Dimension of the epsilon tensor (default `dim`).
#' @param metric_symmetry `1` symmetric, `-1` antisymmetric, `0` no
#'   symmetry, `NULL` no metric.
#' @param metric_name Name of the metric tensor.
#' @param tensor_index_type,typ A `TensorIndexType`.
#' @param is_up `TRUE` for a contravariant index.
#' @param s A string of index names, as in [symbols()].
#' @param ... For `TensorSymmetry()`: a BSGS `list(base, gens)` or `base`
#'   and `gens`; for `TensorSymmetry.direct_product()`: integers, positive
#'   for symmetric and negative for antisymmetric groups of indices.
#' @param rank Number of indices; negative for antisymmetric.
#' @return A `TensorIndexType`, `TensorIndex`, list of indices or
#'   `TensorSymmetry`.
#' @examples
#' L <- TensorIndexType("L", dummy_name = "L")
#' idx <- tensor_indices("i, j", L)
#' -idx$j
#' L$metric
#' TensorSymmetry.fully_symmetric(-2)
#' TensorSymmetry.riemann()
#' @name tensor-indices
NULL

#' @rdname tensor-indices
#' @export
TensorIndexType <- function(name, dummy_name = NULL, dim = NULL, eps_dim = NULL,
                            metric_symmetry = 1, metric_name = "metric") {
  if (inherits(name, "sym")) name <- name$name
  if (is.null(dummy_name)) dummy_name <- substr(name, 1L, 1L)
  if (inherits(dummy_name, "sym")) dummy_name <- dummy_name$name
  dim <- if (is.null(dim)) Symbol(paste0("dim_", dummy_name)) else as_sym(dim)
  eps_dim <- if (is.null(eps_dim)) dim else as_sym(eps_dim)
  ms <- if (is.null(metric_symmetry)) NULL else as_sym(metric_symmetry)
  if (inherits(metric_name, "sym")) metric_name <- metric_name$name
  key <- paste0("TensorIndexType(", name, ",", dummy_name, ",", dim$key, ",", eps_dim$key, ",",
                if (is.null(ms)) "None" else ms$key, ",", metric_name, ")")
  x <- structure(list(name = name, dummy_name = dummy_name, dim = dim, eps_dim = eps_dim,
                      metric_symmetry = ms, metric_name = metric_name, key = key),
                 class = "TensorIndexType")
  x[["skey"]] <- .tn_tit_key(x)
  x
}

#' @export
`$.TensorIndexType` <- function(x, name) {
  switch(name,
         metric = .tit_metric(x),
         delta = TensorHead("KD", list(x, x), TensorSymmetry.fully_symmetric(2)),
         epsilon = {
           ed <- .subset2(x, "eps_dim")
           if (ed$type != "Integer") NULL
           else {
             n <- as.integer(ed$value)
             TensorHead("Eps", rep(list(x), n), TensorSymmetry.fully_symmetric(-n))
           }
         },
         .subset2(x, name))
}

.tit_metric <- function(x) {
  ms <- .subset2(x, "metric_symmetry")
  if (is.null(ms)) return(NULL)
  v <- as.numeric(ms$value)
  symmetry <- if (v == 0) TensorSymmetry.no_symmetry(2) else if (v == 1) TensorSymmetry.fully_symmetric(2)
              else if (v == -1) TensorSymmetry.fully_symmetric(-2)
              else .tn_error("UnboundLocalError", "cannot access local variable 'symmetry' where it is not associated with a value")
  TensorHead(.subset2(x, "metric_name"), list(x, x), symmetry)
}

# structural equality (Basic.__eq__) of index types, symmetries and heads
.tn_struct_ops <- function(gen, e1, e2) {
  if (gen %in% c("==", "!=") && !missing(e2)) {
    eq <- (is.list(e2) || is.function(e2)) && (is.list(e1) || is.function(e1)) && identical(.tn_key(e1), .tn_key(e2))
    return(if (gen == "==") eq else !eq)
  }
  stop("unsupported operator: ", gen)
}

#' @export
Ops.TensorIndexType <- function(e1, e2) .tn_struct_ops(.Generic, e1, e2)
#' @export
Ops.TensorSymmetry <- function(e1, e2) .tn_struct_ops(.Generic, e1, e2)
#' @export
Ops.TensorHead <- function(e1, e2) .tn_struct_ops(.Generic, e1, e2)

#' @export
format.TensorIndexType <- function(x, ...) .subset2(x, "name")
#' @export
print.TensorIndexType <- function(x, ...) {
  cat(format(x), "\n")
  invisible(x)
}

# ---------------------------------------------------------------------------
# TensorIndex

.ti_new <- function(name, typ, is_up) {
  structure(list(name = name, tensor_index_type = typ, is_up = is_up,
                 key = paste0("TensorIndex(", name, ",", .subset2(typ, "key"), ",", is_up, ")")),
            class = "TensorIndex")
}

.ti_autogen <- new.env()

#' @rdname tensor-indices
#' @export
TensorIndex <- function(name, tensor_index_type, is_up = TRUE) {
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
  .ti_new(name, tensor_index_type, isTRUE(as.logical(is_up)))
}

#' @export
Ops.TensorIndex <- function(e1, e2) {
  if (missing(e2)) {
    if (.Generic == "-") return(.ti_new(e1$name, e1$tensor_index_type, !e1$is_up))
    if (.Generic == "+") return(e1)
  }
  if (.Generic %in% c("==", "!=")) {
    eq <- inherits(e1, "TensorIndex") && inherits(e2, "TensorIndex") && identical(e1$key, e2$key)
    return(if (.Generic == "==") eq else !eq)
  }
  stop("unsupported operator for TensorIndex: ", .Generic)
}

#' @export
format.TensorIndex <- function(x, ...) if (x$is_up) x$name else paste0("-", x$name)
#' @export
print.TensorIndex <- function(x, ...) {
  cat(format(x), "\n")
  invisible(x)
}

# TensorIndex.__lt__: (index type, name)
.ti_lt <- function(a, b) {
  ta <- a$tensor_index_type
  tb <- b$tensor_index_type
  if (identical(ta$key, tb$key)) return(.cmp_codepoints(a$name, b$name) < 0L)
  .cmp_codepoints(ta$name, tb$name) < 0L
}

#' @rdname tensor-indices
#' @export
tensor_indices <- function(s, typ) {
  if (!is.character(s)) .tn_error("ValueError", "expecting a string")
  a <- .expand_symbol_names(s)
  tilist <- lapply(a, function(n) TensorIndex(n, typ))
  if (length(tilist) == 1L) return(tilist[[1L]])
  names(tilist) <- a
  tilist
}

# ---------------------------------------------------------------------------
# TensorSymmetry

.tsym_new <- function(base, gens) {
  base <- as.integer(base)
  gens <- lapply(gens, .tc_af)
  structure(list(base = base, generators = lapply(gens, .tc_af_new), rank = length(gens[[1L]]) - 2L,
                 key = paste0("TensorSymmetry(", paste(base, collapse = ","), ";",
                              paste(vapply(gens, .tc_key, ""), collapse = ";"), ")")),
            class = "TensorSymmetry")
}

#' @rdname tensor-indices
#' @export
TensorSymmetry <- function(...) {
  args <- list(...)
  if (length(args) == 1L) return(.tsym_new(args[[1L]][[1L]], args[[1L]][[2L]]))
  if (length(args) == 2L) return(.tsym_new(args[[1L]], args[[2L]]))
  .tn_error("TypeError", "bsgs required, either two separate parameters or one tuple")
}

#' @rdname tensor-indices
#' @export
TensorSymmetry.fully_symmetric <- function(rank) {
  bsgs <- if (rank > 0) get_symmetric_group_sgs(rank, FALSE) else if (rank < 0) get_symmetric_group_sgs(-rank, TRUE)
          else list(integer(), list(.tc_af_new(0:1)))
  TensorSymmetry(bsgs)
}

#' @rdname tensor-indices
#' @export
TensorSymmetry.direct_product <- function(...) {
  base <- integer()
  sgs <- list(.tc_af_new(0:1))
  for (arg in list(...)) {
    if (arg > 0) bsgs2 <- get_symmetric_group_sgs(arg, FALSE)
    else if (arg < 0) bsgs2 <- get_symmetric_group_sgs(-arg, TRUE)
    else next
    r <- bsgs_direct_product(base, sgs, bsgs2[[1L]], bsgs2[[2L]])
    base <- r[[1L]]
    sgs <- r[[2L]]
  }
  TensorSymmetry(base, sgs)
}

#' @rdname tensor-indices
#' @export
TensorSymmetry.riemann <- function() TensorSymmetry(riemann_bsgs)

#' @rdname tensor-indices
#' @export
TensorSymmetry.no_symmetry <- function(rank) TensorSymmetry(integer(), list(.tc_af_new(seq_len(rank + 2L) - 1L)))

#' @export
format.TensorSymmetry <- function(x, ...) {
  tup <- function(el) if (length(el) == 1L) paste0("(", el, ",)") else paste0("(", paste(el, collapse = ", "), ")")
  paste0("TensorSymmetry(", tup(as.character(x$base)), ", ", tup(vapply(x$generators, .tn_perm_str, "")), ")")
}
#' @export
print.TensorSymmetry <- function(x, ...) {
  cat(format(x), "\n")
  invisible(x)
}

# ---------------------------------------------------------------------------
# TensorHead

.th <- function(x) if (is.function(x)) attr(x, "th") else x

.th_new <- function(name, index_types, symmetry, comm) {
  d <- list(name = name, index_types = index_types, symmetry = symmetry, comm = comm,
            key = paste0("TensorHead(", name, ";", paste(.tn_keys(index_types), collapse = ","), ";",
                         symmetry$key, ";", comm$key, ")"))
  d$skey <- .tn_th_key(d)
  d
}

.th_callable <- function(d) {
  f <- function(...) .th_call(d, list(...))
  structure(f, class = c("TensorHead", "function"), th = d)
}

.th_rank <- function(d) length(d$index_types)
.th_comm <- function(d) .tm_symbols2i(d$comm)
.th_commutes_with <- function(a, b) .tm_get_comm(.th_comm(a), .th_comm(b))

#' Tensor heads
#'
#' Port of SymPy's `TensorHead` and `tensor_heads()`. A tensor head has a
#' name, the index types of its slots, a monoterm symmetry (default: none)
#' and a commutation group `comm` (see [TensorManager]). The head is a
#' function: `A(i, -j)` builds the tensor with indices `i` (contravariant)
#' and `-j` (covariant); indices may also be strings such as `"-j"`.
#' `A$name`, `A$index_types`, `A$symmetry`, `A$rank` and `A$comm` give its
#' properties.
#'
#' @param name Name; for `tensor_heads()` a string of names as in
#'   [symbols()].
#' @param index_types List of `TensorIndexType` objects.
#' @param symmetry A `TensorSymmetry`, or `NULL` for no symmetry.
#' @param comm Commutation group label.
#' @param s A string of names.
#' @return A `TensorHead` function, or a list of them.
#' @examples
#' L <- TensorIndexType("L", dummy_name = "L")
#' i <- TensorIndex("i", L); j <- TensorIndex("j", L)
#' A <- TensorHead("A", list(L, L), TensorSymmetry.fully_symmetric(2))
#' A(i, -i)
#' A
#' @export
TensorHead <- function(name, index_types, symmetry = NULL, comm = 0) {
  if (inherits(name, "sym") && name$type == "Symbol") name <- name$name
  else if (!is.character(name)) .tn_error("ValueError", "invalid name")
  if (inherits(index_types, "TensorIndexType")) index_types <- list(index_types)
  if (is.null(symmetry)) symmetry <- TensorSymmetry.no_symmetry(length(index_types))
  else if (symmetry$rank != length(index_types)) .tn_error("AssertionError", "")
  if (is.null(comm)) comm <- 2L
  .th_callable(.th_new(name, index_types, symmetry, as_sym(comm)))
}

#' @export
`$.TensorHead` <- function(x, name) {
  d <- attr(x, "th")
  switch(name, comm = .th_comm(d), rank = .th_rank(d), d[[name]])
}

#' @export
format.TensorHead <- function(x, ...) {
  d <- .th(x)
  paste0(d$name, "(", paste(vapply(d$index_types, function(t) t$name, ""), collapse = ","), ")")
}
#' @export
print.TensorHead <- function(x, ...) {
  cat(format(x), "\n")
  invisible(x)
}

#' @rdname TensorHead
#' @export
tensor_heads <- function(s, index_types, symmetry = NULL, comm = 0) {
  if (!is.character(s)) .tn_error("ValueError", "expecting a string")
  names <- .expand_symbol_names(s)
  thlist <- lapply(names, function(n) TensorHead(n, index_types, symmetry, comm))
  if (length(thlist) == 1L) return(thlist[[1L]])
  names(thlist) <- names
  thlist
}

.th_call <- function(d, indices) {
  types <- d$index_types
  updated <- list()
  for (k in seq_len(min(length(indices), length(types)))) {
    idx <- indices[[k]]
    typ <- types[[k]]
    if (is.character(idx)) {
      idx <- gsub(" ", "", trimws(idx), fixed = TRUE)
      updated[[k]] <- if (startsWith(idx, "-")) TensorIndex(substring(idx, 2L), typ, FALSE) else TensorIndex(idx, typ)
    } else {
      updated[[k]] <- idx
    }
  }
  if (length(indices) > length(updated)) updated <- c(updated, indices[(length(updated) + 1L):length(indices)])
  .tensor_doit(.tensor_new(d, updated))
}

# ---------------------------------------------------------------------------
# _IndexStructure: free indices as list(index, pos), dummies as c(pos_up, pos_down)

.tis_new <- function(free, dum, index_types, indices) {
  if (length(dum)) dum <- dum[order(vapply(dum, `[`, 0L, 1L))]
  list(free = free, dum = dum, index_types = index_types, indices = indices,
       ext_rank = length(free) + 2L * length(dum))
}

.tis_free_dum_from_indices <- function(indices) {
  n <- length(indices)
  if (n == 1L) return(list(free = list(list(indices[[1L]], 0L)), dum = list()))
  free <- rep(TRUE, n)
  seen <- list()
  dum <- list()
  for (i in seq_len(n)) {
    index <- indices[[i]]
    contr <- index$is_up
    k <- paste0(index$name, "\r", index$tensor_index_type$key)
    v <- seen[[k]]
    if (!is.null(v)) {
      is_contr <- v[[1L]]
      pos <- v[[2L]]
      if (is_contr) {
        if (contr) .tn_error("ValueError", sprintf("two equal contravariant indices in slots %d and %d", pos, i - 1L))
      } else if (!contr) {
        .tn_error("ValueError", sprintf("two equal covariant indices in slots %d and %d", pos, i - 1L))
      }
      free[pos + 1L] <- FALSE
      free[i] <- FALSE
      dum[[length(dum) + 1L]] <- if (contr) c(i - 1L, pos) else c(pos, i - 1L)
    } else {
      seen[[k]] <- list(contr, i - 1L)
    }
  }
  fr <- lapply(which(free), function(i) list(indices[[i]], i - 1L))
  # free.sort(): by index (type, name), then position
  fr <- .tn_sort(fr, function(a, b) .ti_lt(a[[1L]], b[[1L]]) || (!.ti_lt(b[[1L]], a[[1L]]) && a[[2L]] < b[[2L]]))
  list(free = fr, dum = dum)
}

.tis_get_free_indices <- function(is) {
  fr <- is$free
  if (!length(fr)) return(list())
  fr <- fr[order(vapply(fr, `[[`, 0L, 2L))]
  lapply(fr, `[[`, 1L)
}

# generator of dummy names, starting above the numbers of free indices named
# like dummies
.tis_dummy_gen <- function(free) {
  cdt <- new.env()
  for (f in free) {
    indx <- f[[1L]]
    typ <- indx$tensor_index_type
    parts <- strsplit(indx$name, "_", fixed = TRUE)[[1L]]
    if (length(parts) && parts[1L] == typ$dummy_name) {
      if (length(parts) < 2L) .tn_error("IndexError", "list index out of range")
      nv <- suppressWarnings(as.integer(parts[2L]))
      if (is.na(nv)) .tn_error("ValueError", sprintf("invalid literal for int() with base 10: '%s'", parts[2L]))
      cur <- get0(typ$key, envir = cdt, inherits = FALSE, ifnotfound = 0L)
      assign(typ$key, max(cur, nv + 1L), envir = cdt)
    }
  }
  function(typ) {
    nd <- get0(typ$key, envir = cdt, inherits = FALSE, ifnotfound = 0L)
    assign(typ$key, nd + 1L, envir = cdt)
    paste0(typ$dummy_name, "_", nd)
  }
}

.tis_replace_dummy_names <- function(indices, free, dum) {
  if (length(dum)) dum <- dum[order(vapply(dum, `[`, 0L, 1L))]
  stopifnot(length(indices) == length(free) + 2L * length(dum))
  gen <- .tis_dummy_gen(free)
  for (d in dum) {
    typ1 <- indices[[d[1L] + 1L]]$tensor_index_type
    nm <- gen(typ1)
    indices[[d[1L] + 1L]] <- .ti_new(nm, typ1, TRUE)
    indices[[d[2L] + 1L]] <- .ti_new(nm, typ1, FALSE)
  }
  indices
}

.tis_from_indices <- function(indices) {
  fd <- .tis_free_dum_from_indices(indices)
  index_types <- lapply(indices, `[[`, "tensor_index_type")
  indices <- .tis_replace_dummy_names(indices, fd$free, fd$dum)
  .tis_new(fd$free, fd$dum, index_types, indices)
}

.tis_generate_indices <- function(free, dum, index_types) {
  indices <- vector("list", length(free) + 2L * length(dum))
  for (f in free) indices[[f[[2L]] + 1L]] <- f[[1L]]
  gen <- .tis_dummy_gen(free)
  for (d in dum) {
    typ1 <- index_types[[d[1L] + 1L]]
    nm <- gen(typ1)
    indices[[d[1L] + 1L]] <- .ti_new(nm, typ1, TRUE)
    indices[[d[2L] + 1L]] <- .ti_new(nm, typ1, FALSE)
  }
  .tis_replace_dummy_names(indices, free, dum)
}

.tis_from_components_free_dum <- function(components, free, dum) {
  index_types <- do.call("c", lapply(components, function(cm) cm$index_types))
  if (is.null(index_types)) index_types <- list()
  indices <- .tis_generate_indices(free, dum, index_types)
  .tis_new(free, dum, index_types, indices)
}

.tis_sorted_free_for_canon <- function(is) .tn_sort(is$free, function(a, b) .ti_lt(a[[1L]], b[[1L]]))
.tis_sorted_dum_for_canon <- function(is) if (length(is$dum)) is$dum[order(vapply(is$dum, `[`, 0L, 1L))] else list()

.metric_symmetry_to_msym <- function(metric) {
  if (is.null(metric)) return(NULL)
  sym <- .th(metric)$symmetry
  if (identical(sym$key, TensorSymmetry.fully_symmetric(2)$key)) return(0)
  if (identical(sym$key, TensorSymmetry.fully_symmetric(-2)$key)) return(1)
  NULL
}

# (g, dummies, msym), the arguments of canonicalize()
.tis_indices_canon_args <- function(is) {
  n <- is$ext_rank
  g <- c(rep(NA_integer_, n), n, n + 1L)
  sf <- .tis_sorted_free_for_canon(is)
  for (i in seq_along(sf)) g[sf[[i]][[2L]] + 1L] <- i - 1L
  pos <- length(is$free)
  j <- length(is$free)
  dummies <- list()
  prev <- NULL
  a <- integer()
  msym <- list()
  for (d in .tis_sorted_dum_for_canon(is)) {
    g[d[1L] + 1L] <- j
    g[d[2L] + 1L] <- j + 1L
    j <- j + 2L
    typ <- is$index_types[[d[1L] + 1L]]
    if (is.null(prev) || !identical(typ$key, prev$key)) {
      if (length(a)) dummies[[length(dummies) + 1L]] <- a
      a <- c(pos, pos + 1L)
      prev <- typ
      msym <- c(msym, list(.metric_symmetry_to_msym(typ$metric)))
    } else {
      a <- c(a, pos, pos + 1L)
    }
    pos <- pos + 2L
  }
  if (length(a)) dummies[[length(dummies) + 1L]] <- a
  list(g = g, dummies = dummies, msym = msym)
}

.tis_perm2tensor <- function(is, g) {
  sorted_free <- lapply(.tis_sorted_free_for_canon(is), `[[`, 1L)
  perm <- .tis_indices_canon_args(is)$g
  rank <- is$ext_rank
  lex_index_types <- vector("list", rank)
  lex_indices <- vector("list", rank)
  for (i in seq_len(rank)) {
    lex_index_types[[perm[i] + 1L]] <- is$index_types[[i]]
    lex_indices[[perm[i] + 1L]] <- is$indices[[i]]
  }
  nfree <- length(sorted_free)
  dum <- rep(list(c(NA_integer_, NA_integer_)), (rank - nfree) %/% 2L)
  free <- list()
  index_types <- vector("list", rank)
  indices <- vector("list", rank)
  for (i in seq_len(rank)) {
    gi <- g[i]
    index_types[[i]] <- lex_index_types[[gi + 1L]]
    indices[[i]] <- lex_indices[[gi + 1L]]
    if (gi < nfree) {
      free[[length(free) + 1L]] <- list(sorted_free[[gi + 1L]], i - 1L)
    } else {
      jj <- gi - nfree
      idum <- jj %/% 2L
      if (jj %% 2L) dum[[idum + 1L]][2L] <- i - 1L else dum[[idum + 1L]][1L] <- i - 1L
    }
  }
  .tis_new(free, dum, index_types, indices)
}

.components_canon_args <- function(components) {
  numtyp <- list()
  prev <- NULL
  for (t in components) {
    if (!is.null(prev) && identical(t$key, prev$key)) {
      numtyp[[length(numtyp)]][[2L]] <- numtyp[[length(numtyp)]][[2L]] + 1L
    } else {
      prev <- t
      numtyp[[length(numtyp) + 1L]] <- list(prev, 1L)
    }
  }
  lapply(numtyp, function(hn) {
    h <- hn[[1L]]
    cm <- .th_comm(h)
    comm <- if (cm %in% c(0L, 1L)) cm else .tm_get_comm(cm, cm)
    list(h$symmetry$base, h$symmetry$generators, hn[[2L]], comm)
  })
}

# ---------------------------------------------------------------------------
# generic accessors (get_indices, get_free_indices, ... in tensor.py)

.tn_get_indices <- function(t) {
  if (inherits(t, "Tensor") || inherits(t, "TensMul") || inherits(t, "PartialDerivative")) return(t$indices)
  if (inherits(t, "TensAdd")) {
    out <- list()
    for (a in t$args) {
      ii <- .tn_get_indices(a)
      out <- c(out, ii[!(.tn_keys(ii) %in% .tn_keys(out))])
    }
    return(out)
  }
  list()
}

.tn_get_free_indices <- function(t) {
  if (inherits(t, "Tensor") || inherits(t, "TensMul")) return(.tis_get_free_indices(t$is))
  if (inherits(t, "PartialDerivative")) return(.tis_get_free_indices(list(free = t$free)))
  if (inherits(t, "TensAdd")) return(if (.tn_is(t$args[[1L]])) .tn_get_free_indices(t$args[[1L]]) else list())
  list()
}

.tn_get_dummy_indices <- function(t) {
  if (!.tn_is(t)) return(list())
  inds <- .tn_get_indices(t)
  fk <- .tn_keys(.tn_get_free_indices(t))
  inds[!(.tn_keys(inds) %in% fk)]
}

.tn_coeff <- function(t) {
  if (inherits(t, "TensMul")) return(t$coeff)
  if (.tn_is(t)) return(S_One)
  t
}

.tn_components <- function(t) {
  if (inherits(t, "Tensor")) return(list(t$head))
  if (inherits(t, "TensMul")) {
    out <- list()
    for (a in t$args) if (inherits(a, "Tensor")) out <- c(out, list(a$head))
    return(out)
  }
  list()
}

.tn_free <- function(t) t$is$free
.tn_dum <- function(t) t$is$dum

# ---------------------------------------------------------------------------
# Tensor

.tn_parse_indices <- function(d, indices) {
  if (inherits(indices, "TensorIndex")) .tn_error("TypeError", "indices should be an array, got <class 'sympy.tensor.tensor.TensorIndex'>")
  for (i in seq_along(indices)) {
    index <- indices[[i]]
    if (inherits(index, "sym") && index$type == "Symbol") {
      indices[[i]] <- TensorIndex(index$name, d$index_types[[i]], TRUE)
    } else if (inherits(index, "sym") && index$type == "Mul") {
      ce <- as_coeff_mul_split(index)
      if (identical(ce[[1L]]$key, S_NegativeOne$key) && ce[[2L]]$type == "Symbol") {
        indices[[i]] <- TensorIndex(ce[[2L]]$name, d$index_types[[i]], FALSE)
      } else {
        .tn_error("ValueError", paste0("index not understood: ", sstr(index)))
      }
    } else if (!inherits(index, "TensorIndex")) {
      .tn_error("TypeError", paste0("wrong type for index: ", .tn_py_str(index), " is ", .tn_py_class(index)))
    }
  }
  indices
}

# str() and type() of a non tensor value, for messages
.tn_py_str <- function(x) if (inherits(x, "sym")) sstr(x) else if (is.list(x) && !is.null(x$key)) format(x) else paste(format(x), collapse = " ")
.tn_py_class <- function(x) {
  if (!inherits(x, "sym")) {
    cl <- if (is.logical(x)) "bool" else if (is.character(x)) "str" else if (is.numeric(x) && all(x == round(x))) "int" else if (is.numeric(x)) "float" else class(x)[1L]
    return(paste0("<class '", cl, "'>"))
  }
  t <- x$type
  path <- switch(t, Add = "sympy.core.add.Add", Mul = "sympy.core.mul.Mul", Pow = "sympy.core.power.Pow",
                 Integer = paste0("sympy.core.numbers.", .sympy_class(x)), Rational = paste0("sympy.core.numbers.", .sympy_class(x)),
                 Float = "sympy.core.numbers.Float", Symbol = "sympy.core.symbol.Symbol",
                 Function = paste0("sympy.core.function.", x$name), paste0("sympy.", t))
  paste0("<class '", path, "'>")
}

.tensor_new <- function(d, indices, is_canon_bp = FALSE) {
  d <- .th(d)
  indices <- .tn_parse_indices(d, indices)
  is <- .tis_from_indices(indices)
  if (.th_rank(d) != length(indices)) .tn_error("ValueError", "wrong number of indices")
  structure(list(head = d, indices = indices, is = is, is_canon_bp = is_canon_bp,
                 key = paste0("Tensor(", d$key, ";", paste(.tn_keys(indices), collapse = ","), ")")),
            class = c("Tensor", "TensExpr"))
}

.tensor_doit <- function(t) .tm_contract_indices(list(t))$args[[1L]]

# xreplace of indices: repl is a list named by index keys
.tn_replace_indices <- function(t, repl) {
  if (inherits(t, "Tensor")) {
    ks <- .tn_keys(t$indices)
    hit <- ks %in% names(repl)
    if (!any(hit)) return(t)
    new <- t$indices
    for (k in which(hit)) new[[k]] <- repl[[ks[k]]]
    if (identical(.tn_keys(new), ks)) return(t)
    return(.tensor_new(t$head, new))
  }
  if (inherits(t, "TensMul")) return(.tensmul_new(lapply(t$args, function(a) if (.tn_is(a)) .tn_replace_indices(a, repl) else a)))
  if (inherits(t, "TensAdd")) return(.tensadd_new(lapply(t$args, function(a) if (.tn_is(a)) .tn_replace_indices(a, repl) else a)))
  if (inherits(t, "PartialDerivative")) {
    expr <- if (.tn_is(t$args[[1L]])) .tn_replace_indices(t$args[[1L]], repl) else t$args[[1L]]
    mirrored <- list()
    for (k in names(repl)) mirrored[[.tn_neg_key(k)]] <- -repl[[k]]
    variables <- lapply(t$args[-1L], function(v) if (.tn_is(v)) .tn_replace_indices(v, mirrored) else v)
    return(.pd_new(expr, variables))
  }
  t
}

# the key of -i from the key of a TensorIndex i
.tn_neg_key <- function(k) {
  if (endsWith(k, ",TRUE)")) paste0(substr(k, 1L, nchar(k) - 6L), ",FALSE)") else paste0(substr(k, 1L, nchar(k) - 7L), ",TRUE)")
}

.tensor_set_indices <- function(t, indices, is_canon_bp = FALSE) {
  if (length(indices) != t$is$ext_rank) .tn_error("ValueError", "indices length mismatch")
  .tensor_doit(.tensor_new(t$head, indices, is_canon_bp))
}

.tn_set_new_index_structure <- function(t, im, is_canon_bp = FALSE) {
  if (inherits(t, "Tensor")) return(.tensor_set_indices(t, im$indices, is_canon_bp))
  .tensmul_set_indices(t, im$indices, is_canon_bp)
}

.tn_perm2tensor <- function(t, g, is_canon_bp = FALSE) {
  if (!.tn_is(t)) return(t)
  if (!(inherits(t, "Tensor") || inherits(t, "TensMul"))) .tn_error("NotImplementedError", "")
  nim <- .tis_perm2tensor(t$is, g)
  res <- .tn_set_new_index_structure(t, nim, is_canon_bp)
  if (g[length(g)] != length(g) - 1L) return(.tn_neg(res))
  res
}

.tensor_canon_bp <- function(t) {
  if (isTRUE(t$is_canon_bp)) return(t)
  ca <- .tis_indices_canon_args(t$is)
  v <- .components_canon_args(list(t$head))
  can <- do.call(canonicalize, c(list(ca$g, ca$dummies, ca$msym), v))
  if (identical(can, 0)) return(S_Zero)
  .tn_perm2tensor(t, can, TRUE)
}

.tensor_substitute_indices <- function(t, index_tuples) {
  indices <- list()
  for (index in t$indices) {
    rep <- NULL
    for (tp in index_tuples) {
      old <- tp[[1L]]
      if (index$name == old$name && identical(index$tensor_index_type$key, old$tensor_index_type$key)) {
        rep <- if (index$is_up == old$is_up) tp[[2L]] else -tp[[2L]]
        break
      }
    }
    indices[[length(indices) + 1L]] <- if (is.null(rep)) index else rep
  }
  .th_call(t$head, indices)
}

.tn_metric_antisym <- function(g) {
  sk <- g$symmetry$key
  if (identical(sk, TensorSymmetry.fully_symmetric(-2)$key)) return(1)
  if (identical(sk, TensorSymmetry.fully_symmetric(2)$key)) return(0)
  if (identical(sk, TensorSymmetry.no_symmetry(2)$key)) return(NULL)
  .tn_error("NotImplementedError", "")
}

.tensor_contract_metric <- function(t, g) {
  g <- .th(g)
  if (!identical(t$head$key, g$key)) return(t)
  if (length(t$is$free)) return(t)
  antisym <- .tn_metric_antisym(g)
  sign <- S_One
  typ <- g$index_types[[1L]]
  sign <- .mul(list(sign, typ$dim))
  if (isTRUE(antisym == 1)) {
    dp <- t$is$dum[[1L]]
    if (dp[1L] < dp[2L]) sign <- .neg(sign)
  }
  sign
}

# ---------------------------------------------------------------------------
# TensMul

.tm_flatten_args <- function(args) {
  flat <- list()
  for (a in args) {
    if (inherits(a, "TensMul")) flat <- c(flat, a$args)
    else if (inherits(a, "sym") && a$type == "Mul") flat <- c(flat, .sympy_args(a))
    else flat[[length(flat) + 1L]] <- a
  }
  flat
}

.tm_indices_to_free_dum <- function(args_indices) {
  fkeys <- character()
  fidx <- list()
  fpos1 <- integer()
  fpos2 <- integer()
  dummy_data <- list()
  indices <- list()
  pos2 <- 0L
  for (pos1 in seq_along(args_indices) - 1L) {
    for (index in args_indices[[pos1 + 1L]]) {
      if (!inherits(index, "TensorIndex")) .tn_error("TypeError", "expected TensorIndex")
      m <- match(.tn_neg_key(index$key), fkeys)
      if (!is.na(m)) {
        other_pos1 <- fpos1[m]
        other_pos2 <- fpos2[m]
        fkeys <- fkeys[-m]; fidx <- fidx[-m]; fpos1 <- fpos1[-m]; fpos2 <- fpos2[-m]
        dummy_data[[length(dummy_data) + 1L]] <- if (index$is_up) list(index, pos1, other_pos1, pos2, other_pos2)
                                                 else list(-index, other_pos1, pos1, other_pos2, pos2)
      } else if (index$key %in% fkeys) {
        .tn_error("ValueError", paste0("Repeated index: ", format(index)))
      } else {
        fkeys <- c(fkeys, index$key); fidx[[length(fidx) + 1L]] <- index
        fpos1 <- c(fpos1, pos1); fpos2 <- c(fpos2, pos2)
      }
      indices[[length(indices) + 1L]] <- index
      pos2 <- pos2 + 1L
    }
  }
  free <- lapply(seq_along(fidx), function(k) list(fidx[[k]], fpos2[k]))
  if (length(dummy_data)) dummy_data <- dummy_data[order(vapply(dummy_data, `[[`, 0L, 4L))]
  list(indices = indices, free = free, free_names = vapply(fidx, `[[`, "", "name"), dummy_data = dummy_data)
}

.tm_contract_indices <- function(args, replace_indices = TRUE) {
  r <- .tm_indices_to_free_dum(lapply(args, .tn_get_indices))
  if (replace_indices) {
    cdt <- new.env()
    gen <- function(typ) {
      nd <- get0(typ$key, envir = cdt, inherits = FALSE, ifnotfound = 0L)
      assign(typ$key, nd + 1L, envir = cdt)
      paste0(typ$dummy_name, "_", nd)
    }
    replacements <- rep(list(list()), length(args))
    for (dd in r$dummy_data) {
      old_index <- dd[[1L]]
      index_type <- old_index$tensor_index_type
      repeat {
        dummy_name <- gen(index_type)
        if (!(dummy_name %in% r$free_names)) break
      }
      dummy <- .ti_new(dummy_name, index_type, old_index$is_up)
      replacements[[dd[[2L]] + 1L]][[old_index$key]] <- dummy
      replacements[[dd[[3L]] + 1L]][[(-old_index)$key]] <- -dummy
    }
    args <- lapply(seq_along(args), function(k) {
      a <- args[[k]]
      if (.tn_is(a) && length(replacements[[k]])) .tn_replace_indices(a, replacements[[k]])
      else if (.tn_is(a)) .tn_replace_indices(a, list())
      else a
    })
    r <- .tm_indices_to_free_dum(lapply(args, .tn_get_indices))
  }
  dum <- lapply(r$dummy_data, function(x) c(x[[4L]], x[[5L]]))
  list(args = args, indices = r$indices, free = r$free, dum = dum)
}

.tm_dedupe_indices <- function(new, exclude) {
  exk <- .tn_keys(exclude)
  dums_new <- .tn_uniq(.tn_get_dummy_indices(new))
  free_new <- .tn_get_free_indices(new)
  conflicts <- dums_new[.tn_keys(dums_new) %in% exk]
  if (!length(conflicts)) return(NULL)
  exclude <- .tn_uniq(c(exclude, dums_new, free_new))
  gen <- .tis_dummy_gen(lapply(exclude, function(i) list(i, NA)))
  repl <- list()
  for (d in conflicts) {
    if ((-d)$key %in% names(repl)) next
    newname <- gen(d$tensor_index_type)
    new_d <- .ti_new(newname, d$tensor_index_type, d$is_up)
    repl[[d$key]] <- new_d
    repl[[(-d)$key]] <- -new_d
  }
  if (!length(repl)) return(NULL)
  .tn_replace_indices(new, repl)
}

.tensmul_new <- function(args, is_canon_bp = FALSE) {
  args <- lapply(args, .tn_sym)
  free <- .tn_uniq(do.call("c", c(list(list()), lapply(args, .tn_get_free_indices))))
  fk <- .tn_keys(free)
  newargs <- list()
  for (arg in args) {
    dum_this <- .tn_get_dummy_indices(arg)
    if (any(.tn_keys(dum_this) %in% fk)) {
      dum_other <- do.call("c", c(list(list()), lapply(newargs, .tn_get_dummy_indices)))
      exclude <- .tn_uniq(c(.tn_get_free_indices(arg), free, dum_other))
      newarg <- .tm_dedupe_indices(arg, exclude)
    } else {
      newarg <- arg
    }
    newargs[[length(newargs) + 1L]] <- newarg
  }
  args <- .tm_flatten_args(newargs)
  r <- .tm_contract_indices(args, replace_indices = FALSE)
  index_types <- lapply(r$indices, `[[`, "tensor_index_type")
  structure(list(args = args, indices = r$indices, index_types = index_types,
                 is = .tis_new(r$free, r$dum, index_types, r$indices), coeff = S_One,
                 is_canon_bp = is_canon_bp,
                 key = paste0("TensMul(", paste(.tn_keys(args), collapse = ";"), ")")),
            class = c("TensMul", "TensExpr"))
}

.tn_doit <- function(x, deep = TRUE) {
  if (inherits(x, "Tensor")) return(.tensor_doit(x))
  if (inherits(x, "PartialDerivative")) return(.pd_new(x$args[[1L]], x$args[-1L]))
  if (inherits(x, "TensMul")) return(.tensmul_doit(x, deep))
  if (inherits(x, "TensAdd")) return(.tensadd_doit(x, deep))
  doit(x, deep = deep)
}

.tensmul_doit <- function(t, deep = TRUE) {
  is_canon_bp <- t$is_canon_bp
  if (deep) {
    args <- lapply(t$args, .tn_doit)
    # rename dummies of the new args that collide with other indices
    exclude <- t$indices
    newargs <- args
    for (k in seq_along(args)) {
      if (.tn_eq(t$args[[k]], args[[k]])) next
      nr <- .tm_dedupe_indices(args[[k]], exclude)
      if (!is.null(nr)) {
        newargs[[k]] <- nr
        exclude <- c(exclude, .tn_get_indices(nr))
      }
    }
    args <- newargs
  } else {
    args <- t$args
  }
  args <- .tm_flatten_args(Filter(function(a) .tn_is(a) || !.tn_is_num(a, 1), args))
  scal <- Filter(function(a) !.tn_is(a), args)
  coeff <- S_One
  for (a in scal) coeff <- .mul(list(coeff, a))
  args <- Filter(.tn_is, args)
  if (!length(args)) return(coeff)
  if (!.tn_is_num(coeff, 1)) args <- c(list(coeff), args)
  if (.tn_is_num(coeff, 0)) return(S_Zero)
  if (length(args) == 1L) return(args[[1L]])
  r <- .tm_contract_indices(args)
  index_types <- lapply(r$indices, `[[`, "tensor_index_type")
  obj <- .tensmul_new(r$args)
  obj$index_types <- index_types
  obj$is <- .tis_new(r$free, r$dum, index_types, r$indices)
  obj$coeff <- coeff
  obj$is_canon_bp <- is_canon_bp
  obj
}

.tensmul_nocoeff <- function(t) .tensmul_doit(.tensmul_new(c(t$args, list(.div(S_One, t$coeff)))), FALSE)

.tensmul_set_indices <- function(t, indices, is_canon_bp = FALSE) {
  if (length(indices) != t$is$ext_rank) .tn_error("ValueError", "indices length mismatch")
  args <- t$args
  pos <- 0L
  for (i in seq_along(args)) {
    a <- args[[i]]
    if (!.tn_is(a)) next
    er <- a$is$ext_rank
    args[[i]] <- .tensor_set_indices(a, indices[pos + seq_len(er)])
    pos <- pos + er
  }
  .tensmul_doit(.tensmul_new(args, is_canon_bp), FALSE)
}

.tensmul_sort_args_for_sorted_components <- function(t) {
  cv <- Filter(.tn_is, t$args)
  sign <- 1L
  n <- length(cv) - 1L
  tkey <- function(x) {
    tys <- .tn_uniq(x$head$index_types)
    tys <- .tn_sort(tys, function(a, b) .cmp_codepoints(a$name, b$name) < 0L)
    list(vapply(tys, `[[`, "", "name"), x$head$name)
  }
  gt <- function(k1, k2) {
    a <- k1[[1L]]; b <- k2[[1L]]
    m <- min(length(a), length(b))
    for (q in seq_len(m)) if (a[q] != b[q]) return(.cmp_codepoints(a[q], b[q]) > 0L)
    if (length(a) != length(b)) return(length(a) > length(b))
    .cmp_codepoints(k1[[2L]], k2[[2L]]) > 0L
  }
  for (i in seq_len(max(n, 0L)) - 1L) {
    for (j in seq(n, i + 1L, by = -1L)) {
      cm <- .tn_commutes_with(cv[[j]], cv[[j + 1L]])
      if (is.null(cm) || !(cm %in% c(0, 1))) next
      if (gt(tkey(cv[[j]]), tkey(cv[[j + 1L]]))) {
        tmp <- cv[[j]]; cv[[j]] <- cv[[j + 1L]]; cv[[j + 1L]] <- tmp
        if (cm) sign <- -sign
      }
    }
  }
  coeff <- .mul(list(as_sym(sign), t$coeff))
  if (!.tn_is_num(coeff, 1)) return(c(list(coeff), cv))
  cv
}

.tn_commutes_with <- function(a, b) {
  if (!.tn_is(b)) return(0)
  if (inherits(a, "Tensor") && inherits(b, "Tensor")) return(.th_commutes_with(a$head, b$head))
  .tn_error("NotImplementedError", "commutes_with between these types")
}

.tensmul_canon_bp <- function(t) {
  if (isTRUE(t$is_canon_bp)) return(t)
  expr <- .tn_expand(t)
  if (inherits(expr, "TensAdd")) return(.tensadd_canon_bp(expr))
  if (!length(.tn_components(expr))) return(expr)
  expr <- .tn_doit(expr, FALSE)
  tt <- if (inherits(expr, "TensMul")) .tn_doit(.tensmul_new(.tensmul_sort_args_for_sorted_components(expr)), FALSE) else expr
  ca <- .tis_indices_canon_args(tt$is)
  v <- .components_canon_args(.tn_components(tt))
  can <- do.call(canonicalize, c(list(ca$g, ca$dummies, ca$msym), v))
  if (identical(can, 0)) return(S_Zero)
  .tn_perm2tensor(tt, can, TRUE)
}

.tensmul_indices_to_args_pos <- function(t) {
  pm <- integer()
  for (arg_i in seq_along(t$args)) {
    a <- t$args[[arg_i]]
    if (!.tn_is(a)) next
    pm <- c(pm, rep(arg_i - 1L, a$is$ext_rank))
  }
  pm
}

.tensmul_contract_metric <- function(t, g) {
  g <- .th(g)
  expr <- .tn_doit(.tn_expand(t), FALSE)
  if (!.tn_eq(t, expr)) return(.tn_contract_metric(.tn_canon_bp(expr), g))
  pm <- .tensmul_indices_to_args_pos(t)
  pmap <- function(p) pm[p + 1L]
  args <- t$args
  antisym <- .tn_metric_antisym(g)
  gpos <- which(vapply(t$args, function(x) inherits(x, "Tensor") && identical(x$head$key, g$key), NA)) - 1L
  if (!length(gpos)) return(t)
  sign <- S_One
  dum <- t$is$dum
  free <- t$is$free
  elim <- integer()
  in_list <- function(x, l) any(vapply(l, function(y) identical(y, x), NA))
  free_in <- function(x, l) any(vapply(l, function(y) identical(y[[2L]], x[[2L]]) && identical(y[[1L]]$key, x[[1L]]$key), NA))
  for (gposx in gpos) {
    if (gposx %in% elim) next
    free1 <- Filter(function(x) pmap(x[[2L]]) == gposx, free)
    dum1 <- Filter(function(x) pmap(x[1L]) == gposx || pmap(x[2L]) == gposx, dum)
    if (!length(dum1)) next
    elim <- c(elim, gposx)
    args[[gposx + 1L]] <- S_One
    if (length(dum1) == 2L) {
      dum10 <- dum1[[1L]]
      dum11 <- dum1[[2L]]
      if (is.null(antisym) || antisym == 0) {
        p0 <- if (pmap(dum10[2L]) == gposx) dum10[1L] else dum10[2L]
        p1 <- if (pmap(dum11[2L]) == gposx) dum11[1L] else dum11[2L]
      } else {
        if (pmap(dum10[2L]) == gposx) {
          p0 <- dum10[1L]
          if (dum10[2L] == 1L) sign <- .neg(sign)
        } else {
          p0 <- dum10[2L]
          if (dum10[1L] == 0L) sign <- .neg(sign)
        }
        if (pmap(dum11[2L]) == gposx) {
          p1 <- dum11[1L]
          sign <- .neg(sign)
        } else {
          p1 <- dum11[2L]
        }
      }
      dum[[length(dum) + 1L]] <- c(p0, p1)
    } else if (length(dum1) == 1L) {
      dp <- dum1[[1L]]
      if (pmap(dp[1L]) == pmap(dp[2L])) {
        sign <- .mul(list(sign, g$index_types[[1L]]$dim))
        if (isTRUE(antisym == 1) && dp[1L] < dp[2L]) sign <- .neg(sign)
      } else {
        if (pmap(dp[1L]) == gposx) {
          p1 <- dp[2L]
          if (isTRUE(antisym == 1) && dp[1L] == 0L) sign <- .neg(sign)
        } else {
          p1 <- dp[1L]
        }
        free[[length(free) + 1L]] <- list(free1[[1L]][[1L]], p1)
      }
    }
    dum <- Filter(function(x) !in_list(x, dum1), dum)
    free <- Filter(function(x) !free_in(x, free1), free)
  }
  # shift positions
  shift <- 0L
  shifts <- integer(length(args))
  for (i in seq_along(args)) {
    if ((i - 1L) %in% elim) {
      shift <- shift + 2L
      next
    }
    shifts[i] <- shift
  }
  free <- lapply(Filter(function(x) !(pmap(x[[2L]]) %in% elim), free),
                 function(x) list(x[[1L]], x[[2L]] - shifts[pmap(x[[2L]]) + 1L]))
  dum <- lapply(Filter(function(x) !(pmap(x[1L]) %in% elim) && !(pmap(x[2L]) %in% elim), dum),
                function(x) c(x[1L] - shifts[pmap(x[1L]) + 1L], x[2L] - shifts[pmap(x[2L]) + 1L]))
  res <- .tn_doit(.tn_mul(sign, .tensmul_new(args)), FALSE)
  if (!.tn_is(res)) return(res)
  im <- .tis_from_components_free_dum(.tn_components(res), free, dum)
  .tn_set_new_index_structure(res, im)
}

.tensmul_from_data <- function(coeff, components, free, dum) {
  is <- .tis_from_components_free_dum(components, free, dum)
  indices <- is$indices
  tensors <- list()
  pos <- 0L
  for (cm in components) {
    r <- .th_rank(cm)
    tensors[[length(tensors) + 1L]] <- .tensor_new(cm, indices[pos + seq_len(r)])
    pos <- pos + r
  }
  .tn_doit(.tensmul_new(c(list(coeff), tensors)), FALSE)
}

# ---------------------------------------------------------------------------
# TensAdd

.tensadd_new <- function(args) {
  args <- lapply(Filter(function(x) .tn_is(x) || !isTRUE(as_sym(x)$type %in% c("Integer", "Rational", "Float") && as_sym(x)$value == 0), args), .tn_sym)
  a <- list()
  for (x in args) {
    if (inherits(x, "TensAdd")) a <- c(a, x$args)
    else if (inherits(x, "sym") && x$type == "Add") a <- c(a, .sympy_args(x))
    else a[[length(a) + 1L]] <- x
  }
  a <- Filter(function(x) !(.tn_is(x) && .tn_is_num(.tn_coeff(x), 0)), a)
  if (length(a) > 1L) {
    keys <- lapply(a, .tn_sort_key)
    o <- .tn_sort(seq_along(a), function(i, j) .tn_key_cmp(keys[[i]], keys[[j]]) < 0)
    a <- a[unlist(o)]
  }
  if (!length(a)) return(S_Zero)
  if (length(a) == 1L) return(a[[1L]])
  structure(list(args = a, key = paste0("TensAdd(", paste(.tn_keys(a), collapse = ";"), ")")),
            class = c("TensAdd", "TensExpr"))
}

.tensadd_check <- function(args) {
  s0 <- sort(.tn_keys(.tn_get_free_indices(args[[1L]])))
  for (a in args[-1L]) {
    if (!identical(sort(.tn_keys(.tn_get_free_indices(a))), s0))
      .tn_error("ValueError", "all tensors must have the same indices")
  }
}

.tensadd_collect_terms <- function(args) {
  setk <- function(x) sort(unique(.tn_keys(.tn_get_free_indices(x))))
  free_indices <- if (.tn_is(args[[1L]])) setk(args[[1L]]) else character()
  tkeys <- character()
  terms <- list()
  coeffs <- list()
  scalars <- S_Zero
  for (arg in args) {
    if (!.tn_is(arg)) {
      if (length(free_indices)) .tn_error("ValueError", "wrong valence")
      scalars <- .add(list(scalars, arg))
      next
    }
    if (!identical(free_indices, setk(arg))) .tn_error("ValueError", "wrong valence")
    nc <- if (inherits(arg, "TensMul")) .tensmul_nocoeff(arg) else arg
    k <- .tn_key(nc)
    m <- match(k, tkeys)
    if (is.na(m)) {
      tkeys <- c(tkeys, k)
      terms[[length(terms) + 1L]] <- nc
      coeffs[[length(coeffs) + 1L]] <- list(.tn_coeff(arg))
    } else {
      coeffs[[m]] <- c(coeffs[[m]], list(.tn_coeff(arg)))
    }
  }
  new_args <- list()
  for (m in seq_along(terms)) {
    cs <- .add(coeffs[[m]])
    if (!.tn_is_num(cs, 0)) new_args[[length(new_args) + 1L]] <- .tn_doit(.tensmul_new(list(cs, terms[[m]])), FALSE)
  }
  if (scalars$type == "Add") new_args <- c(.sympy_args(scalars), new_args)
  else if (!.tn_is_num(scalars, 0)) new_args <- c(list(scalars), new_args)
  new_args
}

.tensadd_doit <- function(t, deep = TRUE) {
  args <- if (deep) lapply(t$args, .tn_doit) else t$args
  args <- Filter(function(a) .tn_is(a) || !.tn_is_num(a, 0), args)
  if (!length(args)) return(S_Zero)
  if (length(args) == 1L) return(args[[1L]])
  .tensadd_check(args)
  args <- .tensadd_collect_terms(args)
  if (!length(args)) return(S_Zero)
  if (length(args) == 1L) return(args[[1L]])
  if (!any(vapply(args, .tn_is, NA))) return(.add(args))
  .tensadd_new(args)
}

.tensadd_canon_bp <- function(t) {
  expr <- .tn_expand(t)
  if (inherits(expr, "TensAdd")) {
    args <- lapply(expr$args, .tn_canon_bp)
    return(.tn_doit(.tensadd_new(args), FALSE))
  }
  .tn_canon_bp(expr)
}

# ---------------------------------------------------------------------------
# expand (Expr.expand with the mul hint applied to TensMul)

.tn_expand <- function(e) {
  if (!.tn_is(e)) return(expand(e))
  if (inherits(e, "Tensor")) return(e)
  sargs <- lapply(e$args, .tn_expand)
  hit <- !identical(.tn_keys(sargs), .tn_keys(e$args))
  if (inherits(e, "PartialDerivative")) return(if (hit) .pd_new(sargs[[1L]], sargs[-1L]) else e)
  if (inherits(e, "TensAdd")) return(if (hit) .tensadd_new(sargs) else e)
  if (hit) e <- .tensmul_new(sargs)
  repeat {
    new <- .tensmul_expand_mul(e)
    if (.tn_eq(new, e) || !inherits(new, "TensMul")) return(new)
    e <- new
  }
}

.tensmul_expand_mul <- function(e) {
  args1 <- lapply(e$args, function(a) {
    if (inherits(a, "TensAdd")) a$args
    else if (inherits(a, "sym") && a$type == "Add") .sympy_args(a)
    else list(a)
  })
  if (all(lengths(args1) == 1L)) return(.tensadd_new(list(.tn_doit(.tensmul_new(lapply(args1, `[[`, 1L)), FALSE))))
  combos <- list(list())
  for (choices in args1) combos <- do.call("c", lapply(combos, function(cmb) lapply(choices, function(ch) c(cmb, list(ch)))))
  .tensadd_new(lapply(combos, function(cmb) .tn_doit(.tensmul_new(cmb), FALSE)))
}

# ---------------------------------------------------------------------------
# arithmetic

.tn_neg <- function(x) {
  if (inherits(x, "TensMul")) return(.tn_doit(.tensmul_new(list(S_NegativeOne, x), x$is_canon_bp), FALSE))
  if (.tn_is(x)) return(.tn_doit(.tensmul_new(list(x, S_NegativeOne)), FALSE))
  .neg(as_sym(x))
}
.tn_add <- function(a, b) .tn_doit(.tensadd_new(list(a, b)), FALSE)
.tn_mul <- function(a, b) .tn_doit(.tensmul_new(list(a, b)), FALSE)

.tn_binop <- function(gen, e1, e2) {
  switch(gen,
         "+" = .tn_add(e1, e2),
         "-" = .tn_add(e1, .tn_neg(e2)),
         "*" = .tn_mul(e1, e2),
         "/" = {
           if (.tn_is(e2) || !.tn_is(e1)) .tn_error("ValueError", "cannot divide by a tensor")
           .tn_mul(e1, .div(S_One, as_sym(e2)))
         },
         "==" = .tn_is(e1) && .tn_is(e2) && .tn_eq(e1, e2),
         "!=" = !(.tn_is(e1) && .tn_is(e2) && .tn_eq(e1, e2)),
         stop("unsupported operator for tensor expressions: ", gen))
}

#' @export
chooseOpsMethod.TensExpr <- function(x, y, mx, my, cl, reverse) TRUE

#' @export
Ops.TensExpr <- function(e1, e2) {
  if (missing(e2)) {
    if (.Generic == "-") return(.tn_neg(e1))
    if (.Generic == "+") return(e1)
    stop("unsupported unary operator for tensor expressions: ", .Generic)
  }
  .tn_binop(.Generic, e1, e2)
}

# a sym on the left of a tensor reaches Ops.sym, which asks the hook first
local({
  old <- .ops_hooks$dispatch
  .ops_hooks$dispatch <- function(gen, e1, e2) {
    if (inherits(e1, "TensExpr") || inherits(e2, "TensExpr")) return(.tn_binop(gen, e1, e2))
    if (is.null(old)) NULL else old(gen, e1, e2)
  }
})

# ---------------------------------------------------------------------------
# printing (StrPrinter._print_Tensor, _print_TensMul, _print_TensAdd)

.tn_str <- function(x) {
  if (inherits(x, "Tensor")) {
    if (.th_rank(x$head) > 0L) return(paste0(x$head$name, "(", paste(vapply(x$indices, format, ""), collapse = ", "), ")"))
    return(x$head$name)
  }
  if (inherits(x, "TensMul")) {
    args <- x$args
    sign <- ""
    if (could_extract_minus_sign(x$coeff)) {
      sign <- "-"
      if (inherits(args[[1L]], "sym") && identical(args[[1L]]$key, S_NegativeOne$key)) args <- args[-1L]
      else args[[1L]] <- .neg(args[[1L]])
    }
    parts <- vapply(args, function(a) {
      if (inherits(a, "TensAdd")) paste0("(", .tn_str(a), ")")
      else if (.tn_is(a)) .tn_str(a)
      else .paren(a, PREC$Mul)
    }, "")
    return(paste0(sign, paste(parts, collapse = "*")))
  }
  if (inherits(x, "TensAdd")) {
    s <- paste(vapply(x$args, .tn_str, ""), collapse = " + ")
    return(gsub("+ -", "- ", s, fixed = TRUE))
  }
  if (inherits(x, "PartialDerivative")) return(paste0("PartialDerivative(", paste(vapply(x$args, .tn_str, ""), collapse = ", "), ")"))
  sstr(x)
}

#' @export
format.TensExpr <- function(x, ...) .tn_str(x)
#' @export
print.TensExpr <- function(x, ...) {
  cat(.tn_str(x), "\n")
  invisible(x)
}
#' @export
as.character.TensExpr <- function(x, ...) .tn_str(x)

# ---------------------------------------------------------------------------
# default_sort_key of tensor objects (TensAdd orders its terms by it)

# .key_cmp, skipping equal subtrees
.tn_key_cmp <- function(a, b) {
  if (identical(a, b)) return(0)
  if (is.list(a) && is.list(b)) {
    for (i in seq_len(min(length(a), length(b)))) {
      r <- .tn_key_cmp(a[[i]], b[[i]])
      if (r != 0) return(r)
    }
    return(sign(length(a) - length(b)))
  }
  .key_cmp(a, b)
}

.tn_basic_key <- function(cls, keys) list(list(5, 0, cls), list(length(keys), keys), .ONE_SORT_KEY, 1)
.tn_symbol_key <- function(name) list(list(2, 0, "Symbol"), list(1, list(name)), .ONE_SORT_KEY, 1)
.tn_tit_key <- function(t) {
  if (!is.null(.subset2(t, "skey"))) return(.subset2(t, "skey"))
  ms <- .subset2(t, "metric_symmetry")
  .tn_basic_key("TensorIndexType", list(.tn_symbol_key(.subset2(t, "name")), .tn_symbol_key(.subset2(t, "dummy_name")),
                                        .sort_key(.subset2(t, "dim")), .sort_key(.subset2(t, "eps_dim")),
                                        if (is.null(ms)) "None" else .sort_key(ms), .tn_symbol_key(.subset2(t, "metric_name"))))
}
.tn_th_key <- function(d) {
  if (!is.null(d$skey)) return(d$skey)
  s <- d$symmetry
  skey <- .tn_basic_key("TensorSymmetry", list(
    .tn_basic_key("Tuple", lapply(s$base, function(b) .sort_key(as_sym(b)))),
    .tn_basic_key("Tuple", lapply(s$generators, function(p) list(list(2, 0, "Permutation"), list(1, list(.tn_perm_str(p))), .ONE_SORT_KEY, 1)))))
  .tn_basic_key("TensorHead", list(.tn_symbol_key(d$name), .tn_basic_key("Tuple", lapply(d$index_types, .tn_tit_key)),
                                   skey, .sort_key(d$comm)))
}
.tn_ti_key <- function(i) {
  .tn_basic_key("TensorIndex", list(.tn_symbol_key(i$name), .tn_tit_key(i$tensor_index_type),
                                    list(list(5, 0, if (i$is_up) "BooleanTrue" else "BooleanFalse"), list(0, list()), .ONE_SORT_KEY, 1)))
}
.tn_sort_key <- function(x) {
  if (inherits(x, "Tensor")) return(.tn_basic_key("Tensor", list(.tn_th_key(x$head), .tn_basic_key("Tuple", lapply(x$indices, .tn_ti_key)))))
  if (inherits(x, "TensMul")) return(.tn_basic_key("TensMul", lapply(x$args, .tn_sort_key)))
  if (inherits(x, "TensAdd")) return(.tn_basic_key("TensAdd", lapply(x$args, .tn_sort_key)))
  if (inherits(x, "PartialDerivative")) return(.tn_basic_key("PartialDerivative", lapply(x$args, .tn_sort_key)))
  .sort_key(as_sym(x))
}

# ---------------------------------------------------------------------------
# public functions

.tn_canon_bp <- function(p) {
  if (inherits(p, "Tensor")) return(.tensor_canon_bp(p))
  if (inherits(p, "TensMul")) return(.tensmul_canon_bp(p))
  if (inherits(p, "TensAdd")) return(.tensadd_canon_bp(p))
  p
}

.tn_contract_metric <- function(t, g) {
  if (inherits(t, "Tensor")) return(.tensor_contract_metric(t, g))
  if (inherits(t, "TensMul")) return(.tensmul_contract_metric(t, g))
  if (inherits(t, "TensAdd")) {
    args <- lapply(t$args, function(x) .tn_contract_metric(x, g))
    return(.tn_canon_bp(.tn_doit(.tensadd_new(args), FALSE)))
  }
  t
}

#' Tensor expressions
#'
#' Ports of SymPy's tensor expression functions. Products and sums of
#' tensors are built with `*`, `+` and `-` (or [TensMul()], [TensAdd()]);
#' repeated indices, one contravariant and one covariant, are contracted
#' and renamed as dummies (`L_0`, `L_1`, ...). `canon_bp()` brings an
#' expression to canonical form with the Butler-Portugal algorithm;
#' `contract_metric(t, g)` contracts the metric `g`; `contract_delta(t,
#' delta)` contracts a Kronecker delta; `substitute_indices(t, ...)` replaces
#' free indices given as `list(old, new)` pairs; `get_free_indices()` and
#' `get_dummy_indices()` list indices; `riemann_cyclic()` applies the cyclic
#' identity of the Riemann tensor; `tensor_mul()` multiplies its arguments.
#'
#' @param t,p A tensor expression (other values are returned unchanged).
#' @param g,delta A metric or delta tensor head.
#' @param ... For `substitute_indices()`, pairs `list(old, new)`; for
#'   `tensor_mul()`, `TensMul()` and `TensAdd()`, the factors or terms.
#' @param t_r A Riemann tensor expression.
#' @param t2 A tensor expression.
#' @param is_canon_bp Whether the product is already canonical.
#' @return A tensor expression or a `sym`; `get_free_indices()` and
#'   `get_dummy_indices()` return lists of `TensorIndex`.
#' @examples
#' L <- TensorIndexType("L", dummy_name = "L")
#' idx <- tensor_indices("i, j", L)
#' A <- TensorHead("A", list(L, L), TensorSymmetry.fully_symmetric(-2))
#' B <- TensorHead("B", list(L, L), TensorSymmetry.fully_symmetric(2))
#' canon_bp(A(idx$i, idx$j) * B(-idx$j, -idx$i))
#' p <- TensorHead("p", list(L))
#' contract_metric(L$metric(idx$i, idx$j) * p(-idx$j), L$metric)
#' @name tensor-expressions
NULL

#' @rdname tensor-expressions
#' @export
canon_bp <- function(p) .tn_canon_bp(p)

#' @rdname tensor-expressions
#' @export
contract_metric <- function(t, g) .tn_contract_metric(t, g)

#' @rdname tensor-expressions
#' @export
contract_delta <- function(t, delta) {
  if (inherits(t, "TensAdd")) {
    args <- lapply(t$args, function(x) if (.tn_is(x)) contract_delta(x, delta) else x)
    return(.tn_canon_bp(.tn_doit(.tensadd_new(args), FALSE)))
  }
  .tn_contract_metric(t, delta)
}

#' @rdname tensor-expressions
#' @export
substitute_indices <- function(t, ...) {
  tuples <- list(...)
  if (inherits(t, "Tensor")) return(.tensor_substitute_indices(t, tuples))
  if (inherits(t, "TensMul") || inherits(t, "TensAdd")) {
    new_args <- lapply(t$args, function(a) if (.tn_is(a)) do.call(substitute_indices, c(list(a), tuples)) else a)
    if (inherits(t, "TensMul")) return(.tn_doit(.tensmul_new(new_args), FALSE))
    return(.tn_doit(.tensadd_new(new_args), FALSE))
  }
  t
}

#' @rdname tensor-expressions
#' @export
get_free_indices <- function(t) .tn_get_free_indices(t)

#' @rdname tensor-expressions
#' @export
get_dummy_indices <- function(t) .tn_get_dummy_indices(t)

#' @rdname tensor-expressions
#' @export
TensMul <- function(..., is_canon_bp = FALSE) .tensmul_new(list(...), is_canon_bp)

#' @rdname tensor-expressions
#' @export
TensAdd <- function(...) .tensadd_new(list(...))

#' @rdname tensor-expressions
#' @export
tensor_mul <- function(...) {
  a <- list(...)
  if (!length(a)) return(.tensmul_from_data(S_One, list(), list(), list()))
  t <- a[[1L]]
  for (tx in a[-1L]) t <- .tn_mul(t, tx)
  t
}

.tn_split <- function(t) {
  if (!inherits(t, "TensMul")) return(list(t))
  splitp <- list()
  res <- S_One
  for (arg in t$args) {
    if (inherits(arg, "Tensor")) {
      splitp[[length(splitp) + 1L]] <- .tn_mul(res, arg)
      res <- S_One
    } else {
      res <- .mul(list(res, arg))
    }
  }
  splitp
}

#' @rdname tensor-expressions
#' @export
riemann_cyclic_replace <- function(t_r) {
  free <- t_r$is$free
  free <- free[order(vapply(free, `[[`, 0L, 2L))]
  m <- free[[1L]][[1L]]; n <- free[[2L]][[1L]]; p <- free[[3L]][[1L]]; q <- free[[4L]][[1L]]
  t0 <- .tn_mul(t_r, make_rational(gmp::as.bigq(2, 3)))
  t1 <- .tn_mul(.tn_neg(substitute_indices(t_r, list(m, m), list(n, q), list(p, n), list(q, p))),
                make_rational(gmp::as.bigq(1, 3)))
  t2 <- .tn_mul(substitute_indices(t_r, list(m, m), list(n, p), list(p, n), list(q, q)),
                make_rational(gmp::as.bigq(1, 3)))
  .tn_add(.tn_add(t0, t1), t2)
}

#' @rdname tensor-expressions
#' @export
riemann_cyclic <- function(t2) {
  t2 <- .tn_expand(t2)
  args <- if (inherits(t2, "TensMul") || inherits(t2, "Tensor")) list(t2) else t2$args
  a2 <- lapply(args, function(x) lapply(.tn_split(x), riemann_cyclic_replace))
  a3 <- lapply(a2, function(v) do.call(tensor_mul, v))
  t3 <- .tn_doit(.tensadd_new(a3), FALSE)
  if (!.tn_is(t3) && .tn_is_num(t3, 0)) return(t3)
  .tn_canon_bp(t3)
}

# ---------------------------------------------------------------------------
# replace_with_arrays

.tn_arr <- function(x) if (inherits(x, "SymArray")) x else if (inherits(x, "symmatrix")) {
  Array(lapply(seq_len(x$nrow), function(i) lapply(seq_len(x$ncol), function(j) .mget(x, i, j))))
} else Array(x)

# str() of an array (nested lists) and of a replacement dict, for messages
.tn_arr_str <- function(a) {
  if (!inherits(a, "SymArray")) return(sstr(a))
  nest <- function(data, shape) {
    if (length(shape) == 1L) return(paste0("[", paste(vapply(data, sstr, ""), collapse = ", "), "]"))
    m <- length(data) %/% shape[1L]
    paste0("[", paste(vapply(seq_len(shape[1L]), function(k) nest(data[(k - 1L) * m + seq_len(m)], shape[-1L]), ""), collapse = ", "), "]")
  }
  if (length(a$shape)) nest(a$data, a$shape) else sstr(a$data[[1L]])
}
.tn_rd_str <- function(rd) {
  paste0("{", paste(vapply(rd, function(p) paste0(if (inherits(p[[1L]], "TensorIndexType")) format(p[[1L]]) else .tn_str(p[[1L]]),
                                                    ": ", .tn_arr_str(p[[2L]])), ""), collapse = ", "), "}")
}

# permutedims with a 0 based permutation, validated as SymPy's Permutation
.tn_permutedims <- function(a, perm) {
  if (length(perm) && !identical(sort(as.integer(perm)), seq_len(max(perm) + 1L) - 1L))
    .tn_error("ValueError", sprintf("Integers 0 through %d must be present.", max(perm)))
  if (length(perm) != length(a$shape)) .tn_error("ValueError", "wrong permutation size")
  if (!length(perm)) return(a)
  permutedims(a, perm + 1L)
}

.tn_arr_zeros <- function(shape) Array(rep(list(S_Zero), prod(shape)), shape)
.tn_arr_add <- function(a, b) Array(Map(function(u, v) .add(list(u, v)), a$data, b$data), a$shape)
.tn_arr_scale <- function(c, a) {
  if (!inherits(a, "SymArray")) return(.mul(list(c, as_sym(a))))
  Array(lapply(a$data, function(u) .mul(list(c, u))), a$shape)
}

# tensorcontraction(a, *pairs) with 0 based axes, all pairs at once
.tn_contract_pairs <- function(a, pairs) {
  if (!length(pairs)) return(a)
  sh <- a$shape
  cax <- unlist(lapply(pairs, function(p) p + 1L))
  rest <- setdiff(seq_along(sh), cax)
  csh <- vapply(pairs, function(p) sh[p[1L] + 1L], 0L)
  cind <- .arr_indices(csh)
  entry <- function(ridx) {
    .add(lapply(seq_len(nrow(cind)), function(r) {
      idx <- integer(length(sh))
      idx[rest] <- ridx
      for (k in seq_along(pairs)) idx[pairs[[k]] + 1L] <- cind[r, k]
      arr_get(a, idx)
    }))
  }
  if (!length(rest)) return(entry(integer()))
  rind <- .arr_indices(sh[rest])
  Array(lapply(seq_len(nrow(rind)), function(r) entry(rind[r, ])), sh[rest])
}

.tn_data_contract_dum <- function(arrays, dum) {
  arrays <- lapply(arrays, function(a) if (inherits(a, "SymArray")) a else Array(list(as_sym(a)), integer()))
  prod <- do.call(tensorproduct, arrays)
  .tn_contract_pairs(prod, dum)
}

.tn_inverse_matrix <- function(a) {
  m <- Matrix(lapply(seq_len(a$shape[1L]), function(i) lapply(seq_len(a$shape[2L]), function(j) arr_get(a, c(i, j)))))
  .tn_arr(inverse(m))
}

.tn_contract_and_permute_with_metric <- function(metric, array, pos, dim) {
  array <- .tn_contract_pairs(tensorproduct(metric, array), list(c(1L, 2L + pos)))
  permu <- seq_len(dim) - 1L
  permu <- append(permu[-1L], permu[1L], after = pos)
  .tn_permutedims(array, permu)
}

.tn_rd_get <- function(rd, key) {
  for (p in rd) if (identical(.tn_key(p[[1L]]), key)) return(p[[2L]])
  NULL
}

.tn_match_indices_with_other_tensor <- function(array, free_ind1, free_ind2, rd) {
  index_types1 <- lapply(free_ind1, `[[`, "tensor_index_type")
  pos2up <- integer()
  pos2down <- integer()
  rem <- .tn_keys(free_ind2)
  incompatible <- function() .tn_error("ValueError", paste0("incompatible indices: [", paste(vapply(free_ind1, format, ""), collapse = ", "),
                                                            "] and [", paste(vapply(free_ind2, format, ""), collapse = ", "), "]"))
  for (pos1 in seq_along(free_ind1)) {
    index1 <- free_ind1[[pos1]]
    m <- match(index1$key, rem)
    if (!is.na(m)) {
      rem[m] <- NA
      next
    }
    m <- match((-index1)$key, rem)
    if (!is.na(m)) {
      rem[m] <- NA
      free_ind2[[m]] <- index1
      if (index1$is_up) pos2up <- c(pos2up, m - 1L) else pos2down <- c(pos2down, m - 1L)
    } else {
      if (pos1 > length(rem) || is.na(rem[pos1])) incompatible()
      index2 <- free_ind2[[pos1]]
      rem[pos1] <- NA
      free_ind2[[pos1]] <- index1
      if (xor(index1$is_up, index2$is_up)) {
        if (index1$is_up) pos2up <- c(pos2up, pos1 - 1L) else pos2down <- c(pos2down, pos1 - 1L)
      }
    }
  }
  if (sum(unique(.tn_keys(free_ind1)) %in% .tn_keys(free_ind2)) < length(free_ind1)) incompatible()
  for (pos in pos2up) {
    metric <- .tn_rd_get(rd, index_types1[[pos + 1L]]$key)
    if (is.null(metric)) .tn_error("ValueError", "No metric provided to lower index")
    array <- .tn_contract_and_permute_with_metric(.tn_inverse_matrix(metric), array, pos, length(free_ind1))
  }
  for (pos in pos2down) {
    metric <- .tn_rd_get(rd, index_types1[[pos + 1L]]$key)
    if (is.null(metric)) .tn_error("ValueError", "No metric provided to lower index")
    array <- .tn_contract_and_permute_with_metric(metric, array, pos, length(free_ind1))
  }
  if (length(free_ind1)) {
    permutation <- match(.tn_keys(free_ind1), .tn_keys(free_ind2)) - 1L
    array <- .tn_permutedims(array, permutation)
  }
  if (inherits(array, "SymArray") && !length(array$shape)) array <- array$data[[1L]]
  list(free_ind2, array)
}

.tn_extract_data <- function(t, rd) {
  if (inherits(t, "Tensor")) {
    other <- NULL
    for (p in rd) if (inherits(p[[1L]], "Tensor") && identical(p[[1L]]$head$key, t$head$key)) {
      other <- p[[1L]]
      array <- p[[2L]]
      break
    }
    if (is.null(other)) .tn_error("ValueError", paste0(.tn_str(t), " not found in ", .tn_rd_str(rd)))
    dum1 <- t$is$dum
    dum2 <- other$is$dum
    if (length(dum2)) {
      for (pair in dum2) if (!any(vapply(dum1, function(x) identical(x, pair), NA)))
        .tn_error("NotImplementedError", paste0(.tn_str(other), " with contractions is not implemented"))
      dum1 <- Filter(function(x) !any(vapply(dum2, function(y) identical(x, y), NA)), dum1)
    }
    if (length(dum1)) {
      indices1 <- t$indices
      indices2 <- other$indices
      repl <- list()
      for (d in dum1) {
        repl[[indices2[[d[2L] + 1L]]$key]] <- -indices2[[d[1L] + 1L]]
        for (pos in d) {
          if (xor(indices1[[pos + 1L]]$is_up, indices2[[pos + 1L]]$is_up)) {
            metric <- .tn_rd_get(rd, indices1[[pos + 1L]]$tensor_index_type$key)
            if (indices1[[pos + 1L]]$is_up) metric <- .tn_inverse_matrix(metric)
            array <- .tn_contract_and_permute_with_metric(metric, array, pos, length(indices2))
          }
        }
      }
      other <- .tensor_doit(.tn_replace_indices(other, repl))
      array <- .tn_data_contract_dum(list(array), dum1)
    }
    return(.tn_match_indices_with_other_tensor(array, .tn_get_free_indices(t), .tn_get_free_indices(other), rd))
  }
  if (inherits(t, "TensMul")) {
    targs <- Filter(.tn_is, t$args)
    ex <- lapply(targs, function(a) .tn_extract_data(a, rd))
    coeff <- S_One
    for (a in Filter(function(a) !.tn_is(a), t$args)) coeff <- .mul(list(coeff, a))
    r <- .tm_indices_to_free_dum(lapply(ex, `[[`, 1L))
    dum <- lapply(r$dummy_data, function(x) c(x[[4L]], x[[5L]]))
    free <- r$free
    if (length(free)) free <- free[order(vapply(free, `[[`, 0L, 2L))]
    return(list(lapply(free, `[[`, 1L), .tn_arr_scale(coeff, .tn_data_contract_dum(lapply(ex, `[[`, 2L), dum))))
  }
  if (inherits(t, "TensAdd")) {
    ex <- lapply(t$args, function(a) if (.tn_is(a)) .tn_extract_data(a, rd) else list(list(), a))
    arrays <- lapply(ex, function(e) .tn_arr(e[[2L]]))
    ref_indices <- ex[[1L]][[1L]]
    last_shape <- arrays[[length(arrays)]]$shape
    for (i in seq_along(ex)[-1L]) {
      perm <- match(.tn_keys(ref_indices), .tn_keys(ex[[i]][[1L]])) - 1L
      arrays[[i]] <- .tn_permutedims(arrays[[i]], perm)
    }
    total <- .tn_arr_zeros(last_shape)
    for (a in arrays) total <- .tn_arr_add(total, a)
    return(list(ref_indices, total))
  }
  if (inherits(t, "PartialDerivative")) return(.pd_extract_data(t, rd))
  list(list(), t)
}

#' Replace tensors by arrays
#'
#' Port of SymPy's `TensExpr.replace_with_arrays()`. `replacement_dict` is
#' a list of pairs `list(key, array)`, where `key` is a tensor such as
#' `A(i, j)` (its indices fix the meaning of the array axes) or a
#' `TensorIndexType` (its covariant metric). The result is an [Array()]
#' with axes ordered as `indices` (default: the free indices of `expr`), or
#' an expression when no free index is left.
#'
#' @param expr A tensor expression.
#' @param replacement_dict List of `list(key, array)` pairs.
#' @param indices Optional list of `TensorIndex` giving the axis order.
#' @return A `SymArray` or a `sym`.
#' @examples
#' L <- TensorIndexType("L", dim = 2)
#' i <- TensorIndex("i", L); j <- TensorIndex("j", L)
#' A <- TensorHead("A", list(L))
#' replace_with_arrays(A(i) * A(j), list(list(A(i), Array(list(1, 2)))), list(i, j))
#' @export
replace_with_arrays <- function(expr, replacement_dict, indices = NULL) {
  if (is.null(indices)) indices <- list()
  if (inherits(indices, "TensorIndex")) indices <- list(indices)
  rd <- lapply(replacement_dict, function(p) list(p[[1L]], .tn_arr(p[[2L]])))
  for (p in rd) {
    tensor <- p[[1L]]
    array <- p[[2L]]
    expected <- if (inherits(tensor, "TensorIndexType")) list(tensor$dim, tensor$dim) else lapply(tensor$head$index_types, function(t) t$dim)
    ok <- length(expected) == length(array$shape) &&
      all(vapply(seq_along(expected), function(k) !(expected[[k]]$type %in% c("Integer", "Rational", "Float")) ||
                   identical(expected[[k]]$key, as_sym(array$shape[k])$key), NA))
    if (!ok) .tn_error("ValueError", paste0("shapes for tensor ", if (inherits(tensor, "TensorIndexType")) format(tensor) else .tn_str(tensor),
                                            " expected to be [", paste(vapply(expected, sstr, ""), collapse = ", "),
                                            "], replacement array shape is {array.shape}"))
  }
  ex <- .tn_extract_data(expr, rd)
  .tn_match_indices_with_other_tensor(ex[[2L]], indices, ex[[1L]], rd)[[2L]]
}

# ---------------------------------------------------------------------------
# PartialDerivative (sympy/tensor/toperators.py)

.pd_contract_indices <- function(expr, variables) {
  vov <- list()
  for (v in variables) {
    if (inherits(v, "Tensor")) {
      fr <- .tn_get_free_indices(v)
      repl <- list()
      for (k in fr) repl[[k$key]] <- -k
      vov[[length(vov) + 1L]] <- .tn_replace_indices(v, repl)
    } else if (inherits(v, "sym") && v$type == "Symbol") {
      vov[[length(vov) + 1L]] <- v
    }
  }
  r <- .tm_contract_indices(c(list(expr), vov), replace_indices = TRUE)
  for (i in seq_along(r$args)[-1L]) {
    a <- r$args[[i]]
    if (inherits(a, "Tensor")) {
      repl <- list()
      for (k in .tn_get_free_indices(a)) repl[[k$key]] <- -k
      r$args[[i]] <- .tn_replace_indices(a, repl)
    }
  }
  r
}

.pd_new <- function(expr, variables) {
  if (inherits(expr, "PartialDerivative")) {
    variables <- c(expr$args[-1L], variables)
    expr <- expr$args[[1L]]
  }
  r <- .pd_contract_indices(.tn_sym(expr), variables)
  structure(list(args = r$args, indices = r$indices, free = r$free, dum = r$dum,
                 key = paste0("PartialDerivative(", paste(.tn_keys(r$args), collapse = ";"), ")")),
            class = c("PartialDerivative", "TensExpr"))
}

#' Partial derivative of tensor expressions
#'
#' Port of SymPy's `PartialDerivative` (sympy/tensor/toperators.py): the
#' derivative of a tensor expression `expr` with respect to tensors or
#' symbols. The indices of a tensor variable enter with opposite valence;
#' repeated indices are contracted. Nested derivatives are flattened.
#' `pd_perform_derivative()` evaluates the derivative with Kronecker deltas
#' and metrics; `pd_expand()` applies the sum and product rules.
#'
#' @param expr A tensor expression or a scalar expression.
#' @param ... Variables: tensors or symbols.
#' @param e A `PartialDerivative` or tensor expression.
#' @return A `PartialDerivative`, tensor expression or `sym`.
#' @examples
#' L <- TensorIndexType("L")
#' i <- TensorIndex("i", L); j <- TensorIndex("j", L)
#' A <- TensorHead("A", list(L))
#' PartialDerivative(A(i), A(j))
#' pd_perform_derivative(PartialDerivative(A(i), A(j)))
#' @export
PartialDerivative <- function(expr, ...) .pd_new(expr, list(...))

.tn_free_symbols_nonempty <- function(x) .tn_is(x) || length(free_symbols(as_sym(x))) > 0L

#' @rdname PartialDerivative
#' @export
pd_expand <- function(e) {
  if (inherits(e, "PartialDerivative")) {
    obj <- .pd_new(e$args[[1L]], e$args[-1L])
    expr <- obj$args[[1L]]
    variables <- obj$args[-1L]
    if (!.tn_free_symbols_nonempty(expr)) return(S_Zero)
    if (inherits(expr, "TensAdd")) {
      return(.tensadd_new(lapply(expr$args, function(a) pd_expand(.pd_new(a, variables)))))
    }
    if (inherits(expr, "TensMul")) {
      if (length(variables) == 1L) {
        terms <- list()
        mulargs <- expr$args
        for (ind in seq_along(mulargs)) {
          a <- mulargs[[ind]]
          if (!.tn_is(a) && a$type %in% c("Integer", "Rational", "Float")) next
          d <- pd_expand(.pd_new(a, variables))
          new <- mulargs
          new[[ind]] <- d
          terms[[length(terms) + 1L]] <- .tensmul_new(new)
        }
        return(.tensadd_new(terms))
      }
      result <- expr
      for (v in variables) result <- pd_expand(.pd_new(result, list(v)))
      return(result)
    }
    return(obj)
  }
  if (inherits(e, "Tensor")) return(.tensor_new(e$head, e$indices))
  if (inherits(e, "TensMul")) return(.tensmul_new(lapply(e$args, function(a) if (.tn_is(a)) pd_expand(a) else a)))
  if (inherits(e, "TensAdd")) return(.tensadd_new(lapply(e$args, function(a) if (.tn_is(a)) pd_expand(a) else a)))
  e
}

.tn_diff_wrt <- function(v) inherits(v, "sym") && v$type == "Symbol"

.tn_eval_partial_derivative <- function(t, s) {
  if (inherits(t, "Tensor")) {
    if (!inherits(s, "Tensor")) return(S_Zero)
    if (!identical(t$head$key, s$head$key)) return(S_Zero)
    kd <- list(S_One)
    fs <- .tn_get_free_indices(t)
    fo <- .tn_get_free_indices(s)
    for (count in seq_len(min(length(fs), length(fo)))) {
      iself <- fs[[count]]
      iother <- fo[[count]]
      if (!identical(iself$tensor_index_type$key, iother$tensor_index_type$key)) .tn_error("ValueError", "index types not compatible")
      tit <- iself$tensor_index_type
      dummy <- TensorIndex(paste0("d_", count - 1L), tit, is_up = iself$is_up)
      kd[[length(kd) + 1L]] <- if (iself$is_up == iother$is_up) tit$delta(iself, -iother)
                               else .tensmul_new(list(tit$metric(iself, dummy), tit$delta(-dummy, -iother)))
    }
    return(.tn_doit(.tensmul_new(kd), FALSE))
  }
  if (inherits(t, "TensMul")) {
    terms <- list()
    for (i in seq_along(t$args)) {
      arg <- t$args[[i]]
      d <- if (.tn_is(arg)) .tn_eval_partial_derivative(arg, s) else if (.tn_diff_wrt(s)) diff(arg, s) else S_Zero
      if (.tn_is(d) || !.tn_is_num(d, 0)) {
        new <- t$args
        new[[i]] <- d
        terms[[length(terms) + 1L]] <- .tn_doit(.tensmul_new(new), FALSE)
      }
    }
    return(.tn_doit(.tensadd_new(terms), FALSE))
  }
  if (inherits(t, "TensAdd")) {
    adds <- list()
    for (a in t$args) {
      if (.tn_is(a)) adds[[length(adds) + 1L]] <- .tn_eval_partial_derivative(a, s)
      else if (.tn_diff_wrt(s)) adds[[length(adds) + 1L]] <- diff(a, s)
    }
    return(.tn_doit(.tensadd_new(adds), FALSE))
  }
  .tn_error("NotImplementedError", "partial derivative of this expression")
}

#' @rdname PartialDerivative
#' @export
pd_perform_derivative <- function(e) {
  result <- e$args[[1L]]
  for (v in e$args[-1L]) {
    if (.tn_is(result)) result <- .tn_eval_partial_derivative(result, v)
    else result <- if (.tn_diff_wrt(v)) diff(result, v) else S_Zero
  }
  result
}

.pd_extract_data <- function(t, rd) {
  ex <- .tn_extract_data(t$args[[1L]], rd)
  indices <- ex[[1L]]
  array <- .tn_arr(ex[[2L]])
  for (variable in t$args[-1L]) {
    vx <- .tn_extract_data(variable, rd)
    var_indices <- lapply(vx[[1L]], function(i) -i)
    va <- .tn_arr(vx[[2L]])
    cm <- lapply(va$data, as_coeff_mul_split)
    coeff_array <- lapply(cm, `[[`, 1L)
    var_array <- lapply(cm, `[[`, 2L)
    dim_before <- length(array$shape)
    array <- derive_by_array(array, var_array)
    dim_increase <- length(array$shape) - dim_before
    array <- .tn_permutedims(array, c(seq_len(dim_before) - 1L + dim_increase, seq_len(dim_increase) - 1L))
    varindex <- var_indices[[1L]]
    ind <- .arr_indices(array$shape)
    for (i in seq_along(coeff_array)) {
      for (r in which(ind[, 1L] == i)) {
        k <- .arr_index(array$shape, ind[r, ])
        array$data[[k]] <- .div(array$data[[k]], coeff_array[[i]])
      }
    }
    pos <- match((-varindex)$key, .tn_keys(indices))
    if (!is.na(pos)) {
      array <- .tn_contract_pairs(array, list(c(0L, pos)))
      indices <- indices[-pos]
    } else {
      indices[[length(indices) + 1L]] <- varindex
    }
  }
  list(indices, array)
}
