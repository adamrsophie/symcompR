# Dirac gamma matrices as tensors, after sympy/physics/hep/gamma_matrices.py.
# The objects are built on first use since the tensor module is collated
# after this file.

#' Gamma matrices
#'
#' `GammaMatrix` is the tensor head of the Dirac gamma matrices, with one
#' index of type `LorentzIndex` (dimension 4, dummy name `L`), as in
#' `sympy.physics.hep.gamma_matrices`. `kahane_simplify()` removes the
#' contracted gamma matrices of a product by Kahane's algorithm,
#' `gamma_trace()` is the trace of a product of gamma matrices,
#' `simplify_gamma_expression()` simplifies the gamma matrices of a tensor
#' product, `simplify_gpgp()` contracts products of gamma matrices with
#' momenta and `extract_type_tens()` splits a product into the tensors of
#' one head and the rest.
#'
#' @param expression,t,ex Tensor expressions.
#' @param component A tensor head.
#' @param sort Whether to sort the components first.
#' @return A tensor expression, a matrix or, for `extract_type_tens()`, a
#'   list of two expressions.
#' @examples
#' i0 <- TensorIndex("i0", LorentzIndex); i1 <- TensorIndex("i1", LorentzIndex)
#' kahane_simplify(GammaMatrix(i0) * GammaMatrix(i1) * GammaMatrix(-i0))
#' gamma_trace(GammaMatrix(i0) * GammaMatrix(i1))
#' @name gamma-matrices
NULL

#' @rdname gamma-matrices
#' @export
LorentzIndex <- NULL
delayedAssign("LorentzIndex", TensorIndexType("LorentzIndex", dim = 4, dummy_name = "L"))

#' @rdname gamma-matrices
#' @export
GammaMatrix <- NULL
delayedAssign("GammaMatrix", TensorHead("GammaMatrix", list(LorentzIndex), TensorSymmetry.no_symmetry(1), comm = NULL))

.is_head <- function(t, comp) inherits(t, "Tensor") && identical(t$head$key, .th(comp)$key)

#' @rdname gamma-matrices
#' @export
extract_type_tens <- function(expression, component) {
  if (inherits(expression, "Tensor")) sp <- list(expression)
  else if (inherits(expression, "TensMul")) sp <- expression$args
  else .qstop("ValueError", "wrong type")
  new_expr <- S_One
  residual_expr <- S_One
  for (i in sp) {
    if (.is_head(i, component)) new_expr <- new_expr * i
    else residual_expr <- residual_expr * i
  }
  list(new_expr, residual_expr)
}

#' @rdname gamma-matrices
#' @export
simplify_gamma_expression <- function(expression) {
  er <- extract_type_tens(expression, GammaMatrix)
  .simplify_single_line(er[[1L]]) * er[[2L]]
}

.hep_sorted_components <- function(t) {
  if (!inherits(t, "TensMul")) return(t)
  .tn_doit(.tensmul_new(.tensmul_sort_args_for_sorted_components(t)), FALSE)
}

#' @rdname gamma-matrices
#' @export
simplify_gpgp <- function(ex, sort = TRUE) {
  gk <- .th(GammaMatrix)$key
  one <- function(ex) {
    components <- .tn_components(ex)
    a <- list()
    comp_map <- integer()
    for (i in seq_along(components)) comp_map <- c(comp_map, rep(i - 1L, .th_rank(components[[i]])))
    dum <- lapply(.tn_dum(ex), function(d) c(d[1L], d[2L], comp_map[d[1L] + 1L], comp_map[d[2L] + 1L]))
    for (i in seq_along(components) - 1L) {
      if (!identical(components[[i + 1L]]$key, gk)) next
      for (dx in dum) {
        if (dx[3L] == i) p_pos1 <- dx[4L]
        else if (dx[4L] == i) p_pos1 <- dx[3L]
        else next
        comp1 <- components[[p_pos1 + 1L]]
        if (identical(.th_comm(comp1), 0L) || identical(.th_comm(comp1), 0)) {
          if (.th_rank(comp1) == 1L) a[[length(a) + 1L]] <- c(i, p_pos1)
        }
      }
    }
    if (!length(a)) return(ex)
    elim <- integer()
    tv <- list()
    hit <- TRUE
    coeff <- S_One
    ta <- NULL
    while (hit) {
      hit <- FALSE
      for (i in seq_len(length(a) - 1L)) {
        ai <- a[[i]]
        if (ai[1L] %in% elim) next
        if (ai[1L] != a[[i + 1L]][1L] - 1L) next
        if (!identical(components[[ai[2L] + 1L]]$key, components[[a[[i + 1L]][2L] + 1L]]$key)) next
        elim <- c(elim, ai[1L], ai[2L], a[[i + 1L]][1L], a[[i + 1L]][2L])
        if (is.null(ta)) {
          ta <- .tn_split(ex)
          mu <- TensorIndex("mu", LorentzIndex)
        }
        hit <- TRUE
        if (i == 1L) coeff <- .tn_coeff(ex)
        comp <- components[[ai[2L] + 1L]]
        tx <- .th_call(comp, list(mu)) * .th_call(comp, list(-mu))
        if (length(a) == 2L) tx <- tx * 4
        tv[[length(tv) + 1L]] <- tx
        break
      }
    }
    if (length(tv)) {
      keep <- ta[!((seq_along(ta) - 1L) %in% elim)]
      return(do.call(tensor_mul, c(keep, tv)) * coeff)
    }
    ex
  }
  if (sort) ex <- .hep_sorted_components(ex)
  repeat {
    t <- one(ex)
    if (!.hep_same(t, ex)) ex <- t else return(t)
  }
}

.hep_same <- function(a, b) {
  if (.tn_is(a) && .tn_is(b)) return(.tn_eq(a, b))
  if (.tn_is(a) || .tn_is(b)) return(FALSE)
  .qeq(a, b)
}

#' @rdname gamma-matrices
#' @export
gamma_trace <- function(t) {
  if (inherits(t, "TensAdd")) return(.tensadd_new(lapply(t$args, gamma_trace)))
  t <- .simplify_single_line(t)
  .trace_single_line(t)
}

.simplify_single_line <- function(expression) {
  t12 <- extract_type_tens(expression, GammaMatrix)
  t1 <- t12[[1L]]
  if (!(is_sym(t1) && .qeq(t1, S_One))) t1 <- kahane_simplify(t1)
  t1 * t12[[2L]]
}

.trace_single_line <- function(t) {
  gk <- .th(GammaMatrix)$key
  one <- function(t) {
    t <- .hep_sorted_components(t)
    components <- .tn_components(t)
    ncomps <- length(components)
    g <- LorentzIndex$metric
    hit <- 0L
    i <- 0L
    for (i in seq_len(ncomps) - 1L) {
      if (identical(components[[i + 1L]]$key, gk)) { hit <- 1L; break }
    }
    j <- ncomps
    if (i + hit < ncomps) for (jj in seq.int(i + hit, ncomps - 1L)) {
      if (!identical(components[[jj + 1L]]$key, gk)) { j <- jj; break }
    }
    numG <- j - i
    if (numG == 0L) {
      tcoeff <- .tn_coeff(t)
      return(if (.qtruthy(tcoeff)) .hep_nocoeff(t) else t)
    }
    if (numG %% 2L == 1L) return(.tensmul_from_data(S_Zero, list(), list(), list()))
    a <- .tn_split(t)
    if (numG > 4L) {
      ind1 <- .tn_get_indices(a[[i + 1L]])[[1L]]
      ind2 <- .tn_get_indices(a[[i + 2L]])[[1L]]
      aa <- c(a[seq_len(i)], a[-seq_len(i + 2L)])
      t1 <- do.call(tensor_mul, aa) * g(ind1, ind2)
      t1 <- contract_metric(t1, g)
      args <- list(t1)
      sign <- 1
      if (i + 2L < j) for (k in seq.int(i + 2L, j - 1L)) {
        sign <- -sign
        ind2 <- .tn_get_indices(a[[k + 1L]])[[1L]]
        aa <- c(a[seq_len(i)], a[seq.int(i + 2L, length.out = k - i - 1L)], a[-seq_len(k + 1L)])
        t2 <- sign * do.call(tensor_mul, aa) * g(ind1, ind2)
        t2 <- contract_metric(t2, g)
        t2 <- simplify_gpgp(t2, FALSE)
        args[[length(args) + 1L]] <- t2
      }
      t3 <- .tensadd_new(args)
      return(.trace_single_line(t3))
    }
    t1 <- do.call(.gamma_trace1, a[seq.int(i + 1L, length.out = j - i)])
    a2 <- c(a[seq_len(i)], a[-seq_len(j)])
    t2 <- do.call(tensor_mul, a2)
    t3 <- t1 * t2
    if (!.tn_is(t3) && !.qtruthy(t3)) return(t3)
    contract_metric(t3, g)
  }
  t <- .tn_expand(t)
  if (inherits(t, "TensAdd")) return(.tensadd_new(lapply(t$args, function(x) one(x) * .tn_coeff(x))))
  if (inherits(t, c("Tensor", "TensMul"))) return(.tn_coeff(t) * one(t))
  sym_trace(t)
}

.hep_nocoeff <- function(t) if (inherits(t, "TensMul")) .tensmul_nocoeff(t) else t

.gamma_trace1 <- function(...) {
  a <- list(...)
  gctr <- 4L
  g <- LorentzIndex$metric
  if (!length(a)) return(as_sym(gctr))
  n <- length(a)
  if (n %% 2L == 1L) return(S_Zero)
  ind <- lapply(a, function(x) .tn_get_indices(x)[[1L]])
  if (n == 2L) return(gctr * g(ind[[1L]], ind[[2L]]))
  if (n == 4L)
    return(gctr * (g(ind[[1L]], ind[[2L]]) * g(ind[[3L]], ind[[4L]]) - g(ind[[1L]], ind[[3L]]) * g(ind[[2L]], ind[[4L]]) +
                     g(ind[[1L]], ind[[4L]]) * g(ind[[2L]], ind[[3L]])))
  NULL
}

#' @rdname gamma-matrices
#' @export
kahane_simplify <- function(expression) {
  if (is_sym(expression) && expression$type == "Mul") return(expression)
  if (inherits(expression, "TensAdd")) return(.tensadd_new(lapply(expression$args, kahane_simplify)))
  if (inherits(expression, "Tensor")) return(expression)
  if (!inherits(expression, "TensMul")) .qstop("AssertionError", "")
  for (gm in expression$args) if (!.is_head(gm, GammaMatrix)) .qstop("AssertionError", "")
  free <- .tn_free(expression)
  lk <- LorentzIndex$key
  dum <- list()
  for (dp in .tn_dum(expression)) {
    if (identical(expression$index_types[[dp[1L] + 1L]]$key, lk)) dum[[length(dum) + 1L]] <- c(dp[1L], dp[2L])
  }
  if (length(dum)) dum <- dum[order(vapply(dum, `[`, 0L, 1L), vapply(dum, `[`, 0L, 2L))]
  if (!length(dum)) return(expression)
  first_dum_pos <- min(vapply(dum, min, 0L))
  total_number <- length(free) + length(dum) * 2L
  number_of_contractions <- length(dum)
  free_pos <- vector("list", total_number)
  for (f in free) free_pos[f[[2L]] + 1L] <- list(f[[1L]])
  index_is_free <- rep(FALSE, total_number)
  for (f in free) index_is_free[f[[2L]] + 1L] <- TRUE
  # links: a dict keyed by position, possibly negative, kept in insertion order
  links <- list()
  K <- function(i) as.character(i)
  has <- function(i) K(i) %in% names(links)
  add <- function(i, v) links[[K(i)]] <<- c(links[[K(i)]], v)
  for (i in seq.int(first_dum_pos, total_number - 1L)) links[K(i)] <- list(integer())
  cum_sign <- -1L
  cum_sign_list <- rep(NA_integer_, total_number)
  block_free_count <- 0L
  resulting_coeff <- S_One
  resulting_indices <- list(list())
  connected_components <- 1L
  for (i in seq_len(total_number) - 1L) {
    if (i < first_dum_pos) next
    if (index_is_free[i + 1L]) {
      block_free_count <- block_free_count + 1L
      if (block_free_count > 1L) { add(i - 1L, i); add(i, i - 1L) }
    } else {
      cum_sign <- cum_sign * (if (block_free_count %% 2L) 1L else -1L)
      if (block_free_count == 0L && i != first_dum_pos) {
        if (cum_sign == -1L) {
          links[[K(-1L - i)]] <- -1L - i + 1L
          links[[K(-1L - i + 1L)]] <- -1L - i
        }
      }
      if (has(i - cum_sign)) {
        if (i != first_dum_pos) add(i, i - cum_sign)
        if (block_free_count != 0L) {
          if (i - cum_sign < total_number && index_is_free[i - cum_sign + 1L]) add(i - cum_sign, i)
        }
      }
      block_free_count <- 0L
    }
    cum_sign_list[i + 1L] <- cum_sign
  }
  for (d in dum) {
    pos1 <- d[1L]; pos2 <- d[2L]
    add(pos1, pos2); add(pos2, pos1)
    linkpos1 <- pos1 + cum_sign_list[pos1 + 1L]
    linkpos2 <- pos2 + cum_sign_list[pos2 + 1L]
    if (linkpos1 >= total_number || linkpos2 >= total_number) next
    if (linkpos1 < first_dum_pos || linkpos2 < first_dum_pos) next
    if (has(-1L - linkpos1)) linkpos1 <- -1L - linkpos1
    if (has(-1L - linkpos2)) linkpos2 <- -1L - linkpos2
    if (linkpos1 >= 0L && !index_is_free[linkpos1 + 1L]) linkpos1 <- pos1
    if (linkpos2 >= 0L && !index_is_free[linkpos2 + 1L]) linkpos2 <- pos2
    if (!(linkpos2 %in% links[[K(linkpos1)]])) add(linkpos1, linkpos2)
    if (!(linkpos1 %in% links[[K(linkpos2)]])) add(linkpos2, linkpos1)
  }
  pop <- function(p) { v <- links[[K(p)]]; links[[K(p)]] <<- NULL; v }
  rm1 <- function(v, x) { k <- match(x, v); if (is.na(k)) v else v[-k] }
  pointer <- first_dum_pos
  previous_pointer <- 0L
  repeat {
    if (has(pointer)) next_ones <- pop(pointer) else break
    if (previous_pointer %in% next_ones) next_ones <- rm1(next_ones, previous_pointer)
    previous_pointer <- pointer
    if (length(next_ones)) pointer <- next_ones[1L] else break
    if (pointer == previous_pointer) break
    if (pointer >= 0L && !is.null(free_pos[[pointer + 1L]]))
      resulting_indices <- lapply(resulting_indices, function(ri) c(ri, list(free_pos[[pointer + 1L]])))
  }
  while (length(links)) {
    connected_components <- connected_components + 1L
    pointer <- min(as.integer(names(links)))
    previous_pointer <- pointer
    prepend_indices <- list()
    repeat {
      if (has(pointer)) next_ones <- pop(pointer) else break
      if (previous_pointer %in% next_ones && length(next_ones) > 1L) next_ones <- rm1(next_ones, previous_pointer)
      previous_pointer <- pointer
      if (length(next_ones)) pointer <- next_ones[1L]
      if (pointer >= first_dum_pos && !is.null(free_pos[[pointer + 1L]]))
        prepend_indices <- c(list(free_pos[[pointer + 1L]]), prepend_indices)
    }
    if (!length(prepend_indices)) {
      resulting_coeff <- resulting_coeff * 2L
    } else {
      expr1 <- prepend_indices
      expr2 <- rev(prepend_indices)
      new <- list()
      for (ri in resulting_indices) for (ex in list(expr1, expr2)) new[[length(new) + 1L]] <- c(ex, ri)
      resulting_indices <- new
    }
  }
  resulting_coeff <- resulting_coeff * (if ((number_of_contractions - connected_components + 1L) %% 2L) -1L else 1L)
  resulting_coeff <- resulting_coeff * 2L^number_of_contractions
  head <- if (first_dum_pos > 0L) free_pos[seq_len(first_dum_pos)] else list()
  resulting_indices <- lapply(resulting_indices, function(ri) c(head, ri))
  resulting_expr <- S_Zero
  for (ri in resulting_indices) {
    temp_expr <- S_One
    for (j in ri) temp_expr <- temp_expr * GammaMatrix(j)
    resulting_expr <- resulting_expr + temp_expr
  }
  t <- resulting_coeff * resulting_expr
  if (!(inherits(t, "TensAdd") || inherits(t, "TensMul"))) t <- eye(4) * t
  t
}

# TensExpr.equals
.hep_equals <- function(self, other) {
  if (inherits(self, "symmatrix")) return(.qeq_mat(self, other))
  if (!.tn_is(self)) return(.qeq(self, other))
  if (is_sym(other) && .qeq(other, S_Zero)) {
    if (inherits(self, "TensAdd")) return(all(vapply(self$args, function(x) .qeq(.tn_coeff(x), S_Zero), NA)))
    return(.qeq(.tn_coeff(self), S_Zero))
  }
  if (!.tn_is(other)) return(.qeq(.tn_coeff(self), other))
  .tn_eq(.tn_canon_bp(self), .tn_canon_bp(other)) || {
    d <- .tn_canon_bp(self - other)
    !.tn_is(d) && .qeq(d, S_Zero)
  }
}

.qeq_mat <- function(a, b) inherits(b, "symmatrix") && a$nrow == b$nrow && a$ncol == b$ncol &&
  all(vapply(seq_along(a$data), function(k) .qeq(simplify(a$data[[k]] - b$data[[k]]), S_Zero), NA))
