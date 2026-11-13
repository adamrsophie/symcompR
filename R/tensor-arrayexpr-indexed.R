# Conversions between array expressions and indexed expressions, after
# sympy/tensor/array/expressions/from_indexed_to_array.py and
# from_array_to_indexed.py.
#
# An index is a sym (a symbol or an integer) or, for indices identified by
# KroneckerDelta factors, a set of symbols: a list with the attribute
# pyset, sorted by key, standing for SymPy's frozenset.

.ix_is_set <- function(i) isTRUE(attr(i, "pyset"))
.ix_key <- function(i) if (.ix_is_set(i)) paste0("{", paste(vapply(i, .k, ""), collapse = ","), "}") else i$key
.ix_set <- function(syms) {
  syms <- syms[!duplicated(vapply(syms, .k, ""))]
  structure(syms[order(vapply(syms, .k, ""), method = "radix")], pyset = TRUE)
}
.ix_keys <- function(l) vapply(l, .ix_key, "")
.ix_in <- function(i, l) .ix_key(i) %in% .ix_keys(l)
.ix_is_int <- function(i) !.ix_is_set(i) && is_sym(i) && i$type == "Integer"

# utils._get_diagonal_indices
.ix_diagonal_indices <- function(flat) {
  keys <- .ix_keys(flat)
  pos <- list()
  for (i in seq_along(flat)) {
    if (.ix_is_int(flat[[i]])) next
    pos[[keys[i]]] <- c(pos[[keys[i]]], i - 1L)
  }
  pos <- pos[vapply(pos, length, 0L) > 1L]
  ret <- flat[!(keys %in% names(pos))]
  dk <- names(pos)[order(match(names(pos), keys))]
  diag <- lapply(dk, function(k) pos[[k]])
  ret <- c(ret, lapply(dk, function(k) flat[[match(k, keys)]]))
  list(diag, ret)
}

.ix_argindex <- function(sub, ind) {
  for (i in seq_along(sub)) {
    s <- sub[[i]]
    if (identical(.ix_key(s), .ix_key(ind))) return(i - 1L)
    if (.ix_is_set(s) && .ix_in(ind, s)) return(i - 1L)
  }
  .perr("IndexError", paste0(.str(ind), " not found"))
}

.ix_is_kd <- function(e) is_sym(e) && e$type == "Function" && identical(e$name, "KroneckerDelta")

#' Convert between indexed and array expressions
#'
#' `convert_indexed_to_array(expr, first_indices)` turns an expression in
#' indexed matrix or array elements, with sums over repeated indices, into an
#' array expression whose axes follow `first_indices`.
#' `convert_array_to_indexed(expr, indices)` writes an array expression as
#' an indexed expression with the given free indices.
#'
#' @param expr An expression.
#' @param first_indices Optional list of indices giving the order of the
#'   axes.
#' @param indices A list of index symbols, one per axis.
#' @return An array expression or an indexed expression.
#' @examples
#' k <- Symbol("k"); i <- Symbol("i"); j <- Symbol("j"); d <- Symbol("d")
#' M <- MatrixSymbol("M", k, k); N <- MatrixSymbol("N", k, k)
#' convert_indexed_to_array(Sum(MatrixElement(M, i, d) * MatrixElement(N, d, j), list(d, 0, k - 1)), list(i, j))
#' convert_array_to_indexed(ArraySymbol("A", list(3, 3)), list(i, j))
#' @export
convert_indexed_to_array <- function(expr, first_indices = NULL) {
  r <- .convert_indexed_to_array(expr)
  result <- r[[1L]]; indices <- r[[2L]]
  if (any(vapply(indices, .ix_is_int, NA))) {
    result <- ArrayElement(result, indices)
    indices <- list()
  }
  if (!length(first_indices)) return(result)
  check_in <- function(el, ind) .ix_in(el, ind) || any(vapply(ind, function(z) .ix_is_set(z) && .ix_in(el, z), NA))
  repl <- list()
  for (s in Filter(.ix_is_set, indices)) for (j in s) repl[[j$key]] <- s
  first <- lapply(first_indices, function(i) { i <- as_sym(i); if (!is.null(repl[[i$key]])) repl[[i$key]] else i })
  first <- Filter(function(i) check_in(i, indices), first)
  first <- c(first, Filter(function(i) !check_in(i, first), indices))
  get_pos <- function(el, ind) {
    if (.ix_in(el, ind)) return(match(.ix_key(el), .ix_keys(ind)) - 1L)
    for (i in seq_along(ind)) if (.ix_is_set(ind[[i]]) && .ix_in(el, ind[[i]])) return(i - 1L)
    .perr("ValueError", "not found")
  }
  perm <- .ae_inv(vapply(indices, get_pos, 0L, first))
  if (is_sym(result) && result$type == "ArrayAdd")
    return(do.call(.ae_array_add, lapply(result$args, .ae_permute_dims, perm)))
  .ae_permute_dims(result, perm)
}

.ix_is_array_expression <- function(e) .is_arrexpr(e)

.convert_indexed_to_array <- function(expr) {
  t <- if (is_sym(expr)) expr$type else ""
  if (t == "Sum") return(.ix_sum(expr))
  if (t == "Mul") {
    conv <- lapply(expr$args, .convert_indexed_to_array)
    args <- lapply(conv, `[[`, 1L); inds <- lapply(conv, `[[`, 2L)
    repl <- list()
    for (a in args) {
      if (!.ix_is_kd(a)) next
      ki <- a$args[1:2]
      for (z in ki) if (!is.null(repl[[z$key]])) ki <- c(ki, repl[[z$key]])
      ks <- .ix_set(ki)
      for (z in ks) repl[[z$key]] <- ks
    }
    keep <- !vapply(args, .ix_is_kd, NA)
    flat <- lapply(do.call(c, c(list(list()), inds[keep])), function(j) if (!.ix_is_set(j) && !is.null(repl[[j$key]])) repl[[j$key]] else j)
    di <- .ix_diagonal_indices(flat)
    tp <- do.call(.ae_tensor_product, args[keep])
    if (length(di[[1L]])) return(list(do.call(.ae_diagonal, c(list(tp), di[[1L]])), di[[2L]]))
    return(list(tp, di[[2L]]))
  }
  if (t == "MatrixElement" || t == "ArrayElement" || t == "Indexed") {
    inds <- if (t == "MatrixElement") expr$args[2:3] else if (t == "ArrayElement") expr$args[[2L]]$args else expr$args[-1L]
    base <- expr$args[[1L]]
    di <- .ix_diagonal_indices(inds)
    if (length(di[[1L]])) return(list(do.call(.ae_diagonal, c(list(base), di[[1L]])), di[[2L]]))
    return(list(base, di[[2L]]))
  }
  if (t == "IndexedBase") .perr("NotImplementedError", "")
  if (.ix_is_kd(expr)) return(list(expr, expr$args[1:2]))
  if (t == "Add") {
    conv <- lapply(expr$args, .convert_indexed_to_array)
    args <- lapply(conv, `[[`, 1L); inds <- lapply(conv, `[[`, 2L)
    if (!any(vapply(inds, length, 0L) > 0L) && !any(vapply(args, .ix_is_array_expression, NA))) return(list(expr, list()))
    index0 <- list(); shape0 <- list()
    for (k in seq_along(args)) {
      ai <- inds[[k]]
      sh <- get_shape(args[[k]])
      for (i in seq_along(ai)) if (!.ix_in(ai[[i]], index0) && !(.ix_key(ai[[i]]) %in% .ix_keys(ai[seq_len(i - 1L)]))) {
        index0[[length(index0) + 1L]] <- ai[[i]]
        shape0[[length(shape0) + 1L]] <- sh[[i]]
      }
    }
    for (k in seq_along(args)) {
      ai <- inds[[k]]
      if (length(ai) < length(index0)) {
        miss <- which(!vapply(index0, function(z) .ix_in(z, ai), NA))
        ai <- c(index0[miss], ai)
        args[[k]] <- .ae_tensor_product(do.call(OneArray, shape0[miss]), args[[k]])
      }
      perm <- vapply(index0, function(j) match(.ix_key(j), .ix_keys(ai)) - 1L, 0L)
      args[[k]] <- .ae_permute_dims(args[[k]], perm)
    }
    return(list(do.call(.ae_array_add, args), index0))
  }
  if (t == "Pow") {
    r <- .convert_indexed_to_array(expr$args[[1L]])
    sub <- r[[1L]]; si <- r[[2L]]
    if (!length(si) && !.ix_is_array_expression(sub)) return(list(expr, list()))
    ex <- expr$args[[2L]]
    if (ex$type == "Integer" && .int_val(ex) >= 1) {
      n <- .int_val(ex)
      if (n == 1) return(list(sub, si))
      if (!length(si)) return(list(do.call(.ae_tensor_product, rep(list(sub), n)), list()))
      diags <- list(2L * (seq_len(n) - 1L), 2L * (seq_len(n) - 1L) + 1L)
      return(list(do.call(.ae_diagonal, c(list(do.call(.ae_tensor_product, rep(list(sub), n))), diags)), si))
    }
    syms <- unlist(lapply(si, function(z) if (.ix_is_set(z)) vapply(z, .k, "") else names(.free_sym_list(z))))
    if (length(syms) && any(vapply(.free_sym_list(ex), function(s) s$key %in% syms, NA)))
      .perr("NotImplementedError", paste0("cannot convert power with exponent depending on the array indices: ", .str(expr)))
    d <- Dummy("d")
    return(list(ArrayElementwiseApplyFunc(Lambda(d, .pow(d, ex)), sub), si))
  }
  if (t == "Function") {
    if (length(expr$args) != 1L)
      .perr("NotImplementedError", paste0("cannot convert multi-argument function to array expression: ", .str(expr)))
    r <- .convert_indexed_to_array(expr$args[[1L]])
    if (!length(r[[2L]]) && !.ix_is_array_expression(r[[1L]])) return(list(expr, list()))
    nm <- expr$name
    return(list(ArrayElementwiseApplyFunc(function(z) fcall(nm, list(z)), r[[1L]]), r[[2L]]))
  }
  list(expr, list())
}

.ix_sum <- function(expr) {
  fun <- expr$args[[1L]]
  lims <- expr$args[-1L]
  r <- .convert_indexed_to_array(fun)
  sub <- r[[1L]]; si <- r[[2L]]
  sets <- list()
  for (s in Filter(.ix_is_set, si)) for (j in s) sets[[j$key]] <- s
  sumi <- lapply(lims, function(l) { v <- l$args[[1L]]; if (!is.null(sets[[v$key]])) sets[[v$key]] else v })
  sumi <- sumi[!duplicated(.ix_keys(sumi))]
  sumi <- sumi[.ae_default_sort_order(lapply(sumi, function(z) if (.ix_is_set(z)) z[[1L]] else z))]
  kron <- list()
  if (fun$type == "Mul") for (a in fun$args) {
    if (!.ix_is_kd(a)) next
    ai <- a$args[1:2]
    ai <- ai[!duplicated(vapply(ai, .k, ""))]
    if (length(ai) == 2L) kron <- c(kron, ai)
  }
  kron_keys <- unique(vapply(kron, .k, ""))
  shape <- get_shape(sub)
  if (length(shape)) for (l in lims) {
    i <- .ix_argindex(si, l$args[[1L]])
    if (!is_zero_num(l$args[[2L]]) || !identical(.add(list(l$args[[3L]], S_One))$key, shape[[i + 1L]]$key))
      .perr("ValueError", paste0("summation index and array dimension mismatch: ", .str(l$args[[1L]])))
  }
  contr <- list()
  if (is_sym(sub) && sub$type == "ArrayDiagonal") {
    dg <- .ae_diag_ind(sub)
    nd <- length(dg)
    dind <- si[(length(si) - nd + 1L):length(si)]
    si <- si[seq_len(length(si) - nd)]
    inner <- length(si) + sum(vapply(dg, length, 0L))
    removed <- list(); remaining <- list(); remd <- list()
    for (g in seq_along(dg)) {
      if (.ix_in(dind[[g]], sumi)) removed[[length(removed) + 1L]] <- dg[[g]]
      else { remaining[[length(remaining) + 1L]] <- dg[[g]]; remd[[length(remd) + 1L]] <- dind[[g]] }
    }
    free_pos <- .ae_diag_push_down(dg, as.list(seq_along(si) - 1L), inner)
    if (length(remaining)) {
      sub <- do.call(.ae_diagonal, c(list(sub$args[[1L]]), remaining))
      up <- .ae_diag_push_up(remaining, list(free_pos, lapply(removed, as.list), lapply(remaining, function(g) g[[1L]])), inner)
      free_pos <- up[[1L]]; removed <- lapply(up[[2L]], unlist); rem_pos <- up[[3L]]
      new_nd <- inner - sum(vapply(remaining, length, 0L) - 1L)
    } else {
      sub <- sub$args[[1L]]
      rem_pos <- list()
      new_nd <- inner
    }
    contr <- c(contr, removed)
    ns <- vector("list", new_nd)
    for (k in seq_along(si)) ns[[free_pos[[k]] + 1L]] <- si[[k]]
    for (k in seq_along(rem_pos)) ns[[rem_pos[[k]] + 1L]] <- remd[[k]]
    si <- ns
  }
  ax <- list()
  for (i in seq_along(si)) {
    ind <- si[[i]]
    if (is.null(ind)) next
    include <- if (.ix_is_set(ind)) all(!(vapply(ind, .k, "") %in% kron_keys)) else !(ind$key %in% kron_keys)
    if (.ix_in(ind, sumi) && include) {
      k <- .ix_key(ind)
      if (is.null(ax[[k]])) ax[[k]] <- list(ind = ind, pos = integer())
      ax[[k]]$pos <- c(ax[[k]]$pos, i - 1L)
      si[i] <- list(NULL)
    }
  }
  for (a in ax) {
    ind <- a$ind
    if (if (.ix_is_set(ind)) any(vapply(ind, .k, "") %in% kron_keys) else ind$key %in% kron_keys) next
    contr[[length(contr) + 1L]] <- a$pos
  }
  free <- Filter(Negate(is.null), si)
  list(do.call(.ae_contraction, c(list(sub), contr)), free)
}

#' @rdname convert_indexed_to_array
#' @export
convert_array_to_indexed <- function(expr, indices) {
  st <- new.env()
  st$count <- 0L
  .ix_do_convert(st, .ae_sympify(expr), lapply(indices, function(i) if (is.null(i)) NULL else as_sym(i)))
}

.ix_do_convert <- function(st, expr, indices) {
  t <- if (is_sym(expr)) expr$type else ""
  if (t == "ArrayTensorProduct") {
    nd <- vapply(expr$args, get_ndim, 0L)
    cum <- c(0L, cumsum(nd))
    return(.mul(lapply(seq_along(expr$args), function(i)
      .ix_do_convert(st, expr$args[[i]], indices[seq_len(cum[i + 1L] - cum[i]) + cum[i]]))))
  }
  if (t == "ArrayContraction") {
    inner <- expr$args[[1L]]
    ni <- vector("list", get_ndim(inner))
    lims <- list()
    bs <- get_shape(inner)
    for (g in .ae_contr_ind(expr)) {
      d <- Dummy(paste0("d", st$count))
      st$count <- st$count + 1L
      lims[[length(lims) + 1L]] <- list(d, S_Zero, .sub(bs[[g[1L] + 1L]], S_One))
      for (i in g) ni[[i + 1L]] <- d
    }
    j <- 1L
    for (i in seq_along(ni)) if (is.null(ni[[i]])) { ni[[i]] <- indices[[j]]; j <- j + 1L }
    return(do.call(Sum, c(list(.ix_do_convert(st, inner, ni)), lims)))
  }
  if (t == "ArrayDiagonal") {
    inner <- expr$args[[1L]]
    ni <- vector("list", get_ndim(inner))
    pos <- .ae_diag_push_down(.ae_diag_ind(expr), as.list(seq_along(indices) - 1L), get_ndim(inner))
    for (k in seq_along(indices)) for (j in unlist(pos[[k]])) ni[[j + 1L]] <- indices[[k]]
    return(.ix_do_convert(st, inner, ni))
  }
  if (t == "PermuteDims") {
    p <- .ae_paf(expr$args[[2L]])
    ni <- vector("list", length(p))
    for (i in seq_along(indices)) ni[[p[i] + 1L]] <- indices[[i]]
    return(.ix_do_convert(st, expr$args[[1L]], ni))
  }
  if (t == "ArrayAdd") return(.add(lapply(expr$args, function(a) .ix_do_convert(st, a, indices))))
  if (t == "Add") return(.add(lapply(expr$args, function(a) if (.ae_is_plain_scalar(a)) a else .ix_do_convert(st, a, indices))))
  if (t == "ArrayElementwiseApplyFunc") return(apply_lambda(expr$args[[1L]], .ix_do_convert(st, expr$args[[2L]], indices)))
  if (t == "Reshape") {
    up <- expr$args[[2L]]$args
    down <- get_shape(expr$args[[1L]])
    cum <- Reduce(function(a, b) .mul(list(a, b)), rev(up), S_One, accumulate = TRUE)
    one <- .add(Map(function(i, s) .mul(list(i, s)), rev(indices), cum[seq_along(indices)]))
    dest <- vector("list", length(down))
    c <- S_One
    rd <- rev(down)
    for (i in seq_along(rd)) {
      e <- rd[[i]]
      last <- i == length(rd)
      dest[[i]] <- if (is_one_num(c)) { if (last) one else fcall("Mod", list(one, e)) }
                   else if (last) fcall("floor", list(.div(one, c)))
                   else fcall("Mod", list(fcall("floor", list(.div(one, c))), e))
      c <- .mul(list(c, e))
    }
    return(.ix_do_convert(st, expr$args[[1L]], rev(dest)))
  }
  if (t == "ZeroArray") return(S_Zero)
  if (t == "OneArray") return(S_One)
  ArrayElement(expr, indices)
}

# Internal names of sympy/tensor/array/expressions, kept for code ported
# from SymPy that calls them directly.
.array_tensor_product <- function(...) .ae_tensor_product(...)
.array_contraction <- function(expr, ...) .ae_contraction(expr, ...)
.array_diagonal <- function(expr, ...) .ae_diagonal(expr, ...)
.permute_dims <- function(expr, permutation) .ae_permute_dims(expr, permutation)
.array_add <- function(...) .ae_array_add(...)
.array2matrix <- function(expr) .a2m(expr)
.remove_trivial_dims <- function(expr) {
  r <- .ae_remove_trivial_dims(expr)
  list(r[[1L]], as.list(r[[2L]]))
}
.support_function_tp1_recognize <- function(contraction_indices, args) .ae_tp1_recognize(contraction_indices, args)
.split_multiple_contractions <- function(expr) .ae_split_multiple_contractions(expr)
.array_diag2contr_diagmatrix <- function(expr) .ae_diag2contr_diagmatrix(expr)
.array_contraction_to_diagonal_multiple_identity <- function(expr) .ae_contraction_to_diagonal_multiple_identity(expr)
.combine_removed <- function(dim, removed1, removed2)
  as.list(.ae_combine_removed(dim, as.integer(unlist(removed1)), as.integer(unlist(removed2))))
.split_scalar_coefficient <- function(expr) .ae_split_coeff(expr)
