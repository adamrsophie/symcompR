# The array contraction editor of sympy/tensor/array/expressions:
# _ArgE (an argument with the contraction index of each of its axes, NA
# for a free axis and a negative number for a diagonal group) and
# _EditArrayContraction, which turns an ArrayContraction, ArrayDiagonal or
# ArrayTensorProduct into a mutable list of such arguments and back.
# Arguments are environments, so they are compared by identity as in SymPy.

.ae_arge <- function(element, indices = NULL) {
  e <- new.env(parent = emptyenv())
  e$element <- element
  e$indices <- if (is.null(indices)) rep(NA_integer_, get_ndim(element)) else as.integer(indices)
  e
}

# _get_mapping_from_sub_ndim_list: flat axis -> (argument, relative axis), 0 based
.ae_mapping <- function(sub) {
  arg <- rep(seq_along(sub) - 1L, sub)
  rel <- unlist(lapply(sub, function(r) seq_len(r) - 1L))
  list(arg = as.integer(arg), rel = as.integer(rel))
}

.ae_editor <- function(base) {
  t <- base$type
  diagd <- list(); ci <- list()
  if (t == "ArrayContraction") {
    expr <- base$args[[1L]]
    ci <- .ae_contr_ind(base)
  } else if (t == "ArrayDiagonal") {
    inner <- base$args[[1L]]
    if (inner$type == "ArrayContraction") {
      expr <- inner$args[[1L]]
      diagd <- lapply(.ae_contr_push_down(.ae_contr_ind(inner), lapply(.ae_diag_ind(base), as.list)), function(g) as.integer(unlist(g)))
      ci <- .ae_contr_ind(inner)
    } else {
      expr <- inner
      diagd <- .ae_diag_ind(base)
    }
  } else if (t == "ArrayTensorProduct") {
    expr <- base
  } else .perr("NotImplementedError", "")
  args <- if (expr$type == "ArrayTensorProduct") expr$args else list(expr)
  ed <- new.env(parent = emptyenv())
  ed$args <- lapply(args, .ae_arge)
  mp <- .ae_mapping(.ae_sub_ndim_list(base))
  for (i in seq_along(ci)) for (j in ci[[i]]) {
    a <- ed$args[[mp$arg[j + 1L] + 1L]]
    a$indices[mp$rel[j + 1L] + 1L] <- i - 1L
  }
  ed$n_contr <- length(ci)
  ed$track <- NULL
  for (i in seq_along(diagd)) for (j in diagd[[i]]) {
    a <- ed$args[[mp$arg[j + 1L] + 1L]]
    a$indices[mp$rel[j + 1L] + 1L] <- -i
  }
  ed
}

.ae_ed_pos <- function(ed, arg) which(vapply(ed$args, identical, NA, arg))[1L]

.ae_ed_insert_after <- function(ed, arg, new_arg) {
  pos <- .ae_ed_pos(ed, arg)
  ed$args <- append(ed$args, list(new_arg), after = pos)
}
.ae_ed_remove <- function(ed, arg) ed$args[[.ae_ed_pos(ed, arg)]] <- NULL

.ae_ed_new_index <- function(ed) {
  ed$n_contr <- ed$n_contr + 1L
  ed$n_contr - 1L
}

.ae_ed_refresh <- function(ed) {
  used <- sort(unique(unlist(lapply(ed$args, function(a) a$indices[!is.na(a$indices)]))))
  ed$n_contr <- length(used)
  for (a in ed$args) a$indices <- match(a$indices, used) - 1L
}

.ae_ed_merge_scalars <- function(ed) {
  sc <- Filter(function(a) !length(a$indices), ed$args)
  ed$args <- Filter(function(a) length(a$indices) > 0L, ed$args)
  elements <- lapply(sc, function(a) {
    el <- a$element
    if (el$type == "ArrayElementwiseApplyFunc" && get_ndim(el) == 0L) {
      conv <- .a2m(el)
      if (!.ae_has_shape(conv)) el <- conv
    }
    el
  })
  scalar <- .mul(elements)
  if (!length(ed$args)) ed$args <- list(.ae_arge(scalar))
  else ed$args[[1L]]$element <- .a2m_tensor_product(scalar, ed$args[[1L]]$element)
}

# hasattr(x, "shape")
.ae_has_shape <- function(x) inherits(x, c("symmatrix", "SymArray")) || .is_arrexpr(x) || .is_matexpr(x)

.ae_ed_to_array_contraction <- function(ed) {
  counter <- 0L
  dkeys <- integer(); dvals <- list()
  allind <- unlist(lapply(ed$args, function(a) a$indices))
  freq <- function(i) if (is.na(i)) sum(is.na(allind)) else sum(allind == i, na.rm = TRUE)
  free_count <- sum(is.na(allind))
  inv1 <- integer(); inv2 <- integer(); done <- integer(); counter4 <- 0L
  for (a in ed$args) {
    counter2 <- 0L
    for (i in a$indices) {
      if (is.na(i)) {
        inv1 <- c(inv1, counter4); counter2 <- counter2 + 1L; counter4 <- counter4 + 1L
        next
      }
      if (i >= 0L) next
      k <- -1L - i
      w <- match(k, dkeys)
      if (is.na(w)) { dkeys <- c(dkeys, k); dvals[[length(dvals) + 1L]] <- counter + counter2 }
      else dvals[[w]] <- c(dvals[[w]], counter + counter2)
      if (freq(i) == 1L && !(i %in% done)) { inv1 <- c(inv1, free_count - 1L - i); done <- c(done, i) }
      else if (!(i %in% done)) { inv2 <- c(inv2, free_count - 1L - i); done <- c(done, i) }
      counter2 <- counter2 + 1L
    }
    a$indices <- ifelse(!is.na(a$indices) & a$indices >= 0L, a$indices, NA_integer_)
    counter <- counter + sum(is.na(a$indices))
  }
  perm <- .ae_inv(c(inv1, inv2))
  diag_f <- Filter(function(v) length(v) > 1L, dvals)
  .ae_ed_merge_scalars(ed)
  .ae_ed_refresh(ed)
  args <- lapply(ed$args, function(a) a$element)
  ci <- .ae_ed_contraction_indices(ed)
  expr <- do.call(.ae_contraction, c(list(do.call(.ae_tensor_product, args)), ci))
  expr2 <- do.call(.ae_diagonal, c(list(expr), diag_f))
  if (!is.null(ed$track)) expr2 <- .ae_permute_dims(expr2, .ae_inv(unlist(ed$track)))
  .ae_permute_dims(expr2, perm)
}

.ae_ed_contraction_indices <- function(ed) {
  out <- rep(list(integer()), ed$n_contr)
  pos <- 0L
  for (a in ed$args) for (j in a$indices) {
    if (!is.na(j)) out[[j + 1L]] <- c(out[[j + 1L]], pos)
    pos <- pos + 1L
  }
  out
}

# list of (argument position, relative axis), 0 based
.ae_ed_mapping_for_index <- function(ed, ind) {
  if (ind >= ed$n_contr) .perr("ValueError", "index value exceeding the index range")
  out <- list()
  for (i in seq_along(ed$args)) {
    ix <- ed$args[[i]]$indices
    for (j in seq_along(ix)) if (!is.na(ix[j]) && ix[j] == ind) out[[length(out) + 1L]] <- c(i - 1L, j - 1L)
  }
  out
}

.ae_in <- function(index, v) if (is.na(index)) anyNA(v) else any(!is.na(v) & v == index)
.ae_ed_count_args_with_index <- function(ed, index) sum(vapply(ed$args, function(a) .ae_in(index, a$indices), NA))
.ae_ed_args_with_index <- function(ed, index) Filter(function(a) .ae_in(index, a$indices), ed$args)
.ae_ed_n_diag <- function(ed) length(unique(unlist(lapply(ed$args, function(a) a$indices[!is.na(a$indices) & a$indices < 0L]))))

.ae_ed_track_start <- function(ed) {
  perm <- list(); pdiag <- integer(); counter <- 0L; counter2 <- -1L
  for (a in ed$args) {
    p <- integer()
    for (i in a$indices) {
      if (!is.na(i)) {
        if (i < 0L) { pdiag <- c(pdiag, counter2); counter2 <- counter2 - 1L }
        next
      }
      p <- c(p, counter); counter <- counter + 1L
    }
    perm[[length(perm) + 1L]] <- p
  }
  max_ind <- if (length(perm)) max(vapply(perm, function(p) if (length(p)) max(p) else -1L, 0L)) else -1L
  ed$track <- c(perm, list(max_ind - pdiag))
}

.ae_ed_track_merge <- function(ed, dest, from) {
  i <- .ae_ed_pos(ed, dest); j <- .ae_ed_pos(ed, from)
  ed$track[[i]] <- c(ed$track[[i]], ed$track[[j]])
  ed$track[[j]] <- NULL
}

.ae_ed_abs_free_range <- function(ed, arg) {
  counter <- 0L
  for (a in ed$args) {
    n <- sum(is.na(a$indices))
    if (identical(a, arg)) return(c(counter, counter + n))
    counter <- counter + n
  }
  .perr("IndexError", "argument not found")
}

.ae_ed_abs_range <- function(ed, arg) {
  counter <- 0L
  for (a in ed$args) {
    n <- length(a$indices)
    if (identical(a, arg)) return(c(counter, counter + n))
    counter <- counter + n
  }
  .perr("IndexError", "argument not found")
}

#' @noRd
.ae_split_multiple_contractions <- function(self) {
  ed <- .ae_editor(self)
  ci <- .ae_contr_ind(self)
  onearray_insert <- list()
  for (indl in seq_along(ci) - 1L) {
    links <- ci[[indl + 1L]]
    if (length(links) <= 2L) next
    positions <- .ae_ed_mapping_for_index(ed, indl)
    aig <- vapply(positions, `[`, 0L, 1L)
    if (anyDuplicated(aig)) next
    cur_dim <- .ae_shape(self$args[[1L]])[[links[1L] + 1L]]
    not_vectors <- list(); vectors <- list()
    for (p in positions) {
      arg <- ed$args[[p[1L] + 1L]]
      mat <- arg$element
      rng <- .ae_ed_abs_range(ed, arg)
      other_abs <- rng[1L] + 1L - p[2L]
      msh <- .ae_shape(mat)
      one_in <- any(vapply(msh, is_one_num, NA))
      is11 <- length(msh) == 2L && all(vapply(msh, is_one_num, NA))
      other_in <- any(vapply(seq_along(ci) - 1L, function(li) li != indl && other_abs %in% ci[[li + 1L]], NA))
      if (!one_in || (is_one_num(cur_dim) && !is11) || other_in) not_vectors[[length(not_vectors) + 1L]] <- list(arg, p[2L])
      else vectors[[length(vectors) + 1L]] <- list(arg, p[2L])
    }
    if (length(not_vectors) > 2L) next
    loop <- c(not_vectors[seq_len(min(1L, length(not_vectors)))], vectors, not_vectors[-seq_len(min(1L, length(not_vectors)))])
    inner <- loop[-c(1L, length(loop))]
    for (v in inner) v[[1L]]$element <- diagonalize_vector(v[[1L]]$element)
    first <- loop[[1L]]
    new_index <- first[[1L]]$indices[first[[2L]] + 1L]
    for (v in inner) {
      a <- v[[1L]]; rel <- v[[2L]]
      a$indices[rel + 1L] <- new_index
      new_index <- .ae_ed_new_index(ed)
      a$indices[which(is.na(a$indices))[1L]] <- new_index
      onearray_insert[[length(onearray_insert) + 1L]] <- a
    }
    last <- loop[[length(loop)]]
    last[[1L]]$indices[last[[2L]] + 1L] <- new_index
  }
  for (v in onearray_insert) .ae_ed_insert_after(ed, v, .ae_arge(OneArray(1L), NA_integer_))
  .ae_ed_to_array_contraction(ed)
}

.ae_flatten_contraction_of_diagonal <- function(self) {
  inner <- self$args[[1L]]
  if (inner$type != "ArrayDiagonal") return(self)
  nd <- get_ndim(inner$args[[1L]])
  down <- .ae_diag_push_down(.ae_diag_ind(inner), lapply(.ae_contr_ind(self), as.list), nd)
  down <- lapply(down, function(g) as.integer(unlist(g)))
  new_ci <- list()
  diag <- .ae_diag_ind(inner)
  for (g in down) {
    grp <- g
    for (j in g) {
      w <- vapply(diag, function(k) j %in% k, NA)
      grp <- c(grp, unlist(diag[w]))
      diag <- diag[!w]
    }
    new_ci[[length(new_ci) + 1L]] <- sort(unique(grp))
  }
  new_ci <- lapply(.ae_diag_push_up(diag, lapply(new_ci, as.list), nd), function(g) as.integer(unlist(g)))
  do.call(.ae_contraction, c(list(do.call(.ae_diagonal, c(list(inner$args[[1L]]), diag))), new_ci))
}
