# Conversions between matrix expressions and array expressions, after
# sympy/tensor/array/expressions/from_matrix_to_array.py and
# from_array_to_matrix.py.

.ae_is_type <- function(e, types) is_sym(e) && e$type %in% types

# ---------------------------------------------------------------------------
# matrix expression -> array expression

.ae_elementwise_apply <- function(lam, element) {
  if (!(.is_matexpr(element) || .is_arrexpr(element))) return(.ae_lambda_call(lam, element))
  ArrayElementwiseApplyFunc(lam, element)
}

#' Convert between matrix and array expressions
#'
#' `convert_matrix_to_array()` rewrites a matrix expression as an array
#' expression: products become contractions of tensor products, transposes
#' become axis permutations, traces contractions and Hadamard products
#' diagonals. `convert_array_to_matrix()` recognizes matrix expressions in
#' an array expression, as SymPy's functions of the same names.
#'
#' @param expr A matrix or array expression.
#' @return An array or matrix expression.
#' @examples
#' k <- Symbol("k")
#' A <- MatrixSymbol("A", k, k); B <- MatrixSymbol("B", k, k)
#' cg <- convert_matrix_to_array(A * t(B))
#' cg
#' convert_array_to_matrix(cg)
#' @export
convert_matrix_to_array <- function(expr) .m2a(.ae_sympify(expr))

#' @rdname convert_matrix_to_array
#' @export
convert_array_to_matrix <- function(expr) {
  rec <- .a2m_absorb_scalars(.a2m(.ae_sympify(expr)))
  rec <- .ae_remove_trivial_dims(rec)[[1L]]
  .a2m_absorb_scalars(rec)
}

.m2a <- function(expr) {
  if (!is_sym(expr)) return(expr)
  switch(expr$type,
    MatMul = {
      nonmat <- lapply(Filter(function(a) !.is_matexpr(a), expr$args), .m2a)
      mats <- Filter(.is_matexpr, expr$args)
      contr <- lapply(seq_len(length(mats) - 1L), function(i) c(2L * i - 1L, 2L * i))
      scalar <- if (length(nonmat)) do.call(.ae_tensor_product, nonmat) else S_One
      tprod <- if (is_sym(scalar) && is_one_num(scalar)) do.call(.ae_tensor_product, lapply(mats, .m2a))
               else do.call(.ae_tensor_product, c(list(scalar), lapply(mats, .m2a)))
      do.call(.ae_contraction, c(list(tprod), contr))
    },
    Add = , MatAdd = do.call(.ae_array_add, lapply(expr$args, .m2a)),
    Transpose = .ae_permute_dims(.m2a(expr$args[[1L]]), c(1L, 0L)),
    Trace = {
      inner <- .m2a(expr$args[[1L]])
      .ae_contraction(inner, c(0L, length(.ae_shape(inner)) - 1L))
    },
    Mul = do.call(.ae_tensor_product, lapply(expr$args, .m2a)),
    Pow = {
      base <- .m2a(expr$args[[1L]]); ex <- expr$args[[2L]]
      if (!length(.ae_shape(base))) {
        if (!(.is_matexpr(base) || .is_arrexpr(base))) return(.pow(base, ex))
        if (ex$type == "Integer" && .int_val(ex) > 0) return(do.call(.ae_tensor_product, rep(list(base), .int_val(ex))))
        d <- Dummy("d")
        return(.ae_elementwise_apply(Lambda(d, .pow(d, ex)), base))
      }
      if (.ae_rel_true(ex, "positive")) return(do.call(.ae_tensor_product, rep(list(base), .int_val(ex))))
      expr
    },
    MatPow = {
      base <- .m2a(expr$args[[1L]]); ex <- expr$args[[2L]]
      if (ex$type != "Integer") {
        sh <- .mshape(expr)
        if (is_one_num(sh[[1L]]) && is_one_num(sh[[2L]])) {
          b <- Dummy("b")
          return(.ae_elementwise_apply(Lambda(b, .pow(b, ex)), base))
        }
        return(expr)
      }
      if (.int_val(ex) > 0) return(.m2a(do.call(MatMul, rep(list(expr$args[[1L]]), .int_val(ex)))))
      expr
    },
    HadamardProduct = {
      n <- length(expr$args)
      tp <- do.call(.ae_tensor_product, lapply(expr$args, .m2a))
      .ae_diagonal(tp, 2L * (seq_len(n) - 1L), 2L * (seq_len(n) - 1L) + 1L)
    },
    HadamardPower = {
      base <- expr$args[[1L]]; ex <- expr$args[[2L]]
      if (ex$type == "Integer" && .int_val(ex) > 0) return(.m2a(do.call(HadamardProduct, rep(list(base), .int_val(ex)))))
      d <- Dummy("d")
      .ae_elementwise_apply(Lambda(d, .pow(d, ex)), base)
    },
    KroneckerProduct = {
      kp <- lapply(expr$args, .m2a)
      n <- length(kp)
      perm <- c(2L * (seq_len(n) - 1L), 2L * (seq_len(n) - 1L) + 1L)
      Reshape(.ae_permute_dims(do.call(.ae_tensor_product, kp), perm), .mshape(expr))
    },
    expr)
}

# ---------------------------------------------------------------------------
# array expression -> matrix expression

.a2m_absorb_scalars <- function(expr) {
  if (.ae_is_type(expr, "ArrayTensorProduct")) {
    ca <- .ae_term_coeff_arrays(expr)
    if (length(ca[[2L]]) == 1L && .is_matexpr(ca[[2L]][[1L]]) && .ae_comm(ca[[1L]])) return(.ae_mmul(ca[[1L]], ca[[2L]][[1L]]))
  }
  expr
}

.a2m_tensor_product <- function(...) {
  scalars <- list(); arrays <- list()
  for (a in list(...)) {
    if (.ae_is_type(a, "ArrayElementwiseApplyFunc") && get_ndim(a) == 0L) a <- .a2m(a)
    if (.is_matexpr(a) || .is_arrexpr(a) || inherits(a, c("symmatrix", "SymArray")) || .ae_is_type(a, "ImmutableDenseNDimArray"))
      arrays[[length(arrays) + 1L]] <- a
    else scalars[[length(scalars) + 1L]] <- a
  }
  scalar <- .mul(scalars)
  if (!length(arrays)) return(scalar)
  if (!is_one_num(scalar)) {
    if (.is_arrexpr(arrays[[1L]]) || !.ae_comm(scalar)) arrays <- c(list(scalar), arrays)
    else arrays[[1L]] <- .ae_mmul(arrays[[1L]], scalar)
  }
  if (length(arrays) == 1L) return(arrays[[1L]])
  do.call(.ae_tensor_product, arrays)
}

.a2m_add <- function(...) {
  args <- list(...)
  if (!any(vapply(args, .is_arrexpr, NA))) {
    if (!any(vapply(args, .is_matexpr, NA))) return(.add(args))
    return(doit(do.call(MatAdd, args)))
  }
  do.call(.ae_array_add, args)
}

.a2m_trace <- function(arg) if (.is_arrexpr(arg)) .ae_contraction(arg, c(0L, 1L)) else Trace(arg)
.a2m_transpose <- function(arg) if (.is_arrexpr(arg)) .ae_permute_dims(arg, c(1L, 0L)) else doit(Transpose(arg))

# Trace._normalize: rotate the cyclic product to start at its smallest factor
.ae_trace_normalize <- function(tr) {
  arg <- tr$args[[1L]]
  if (arg$type != "MatMul") return(tr)
  keyof <- function(a) .sort_key(if (a$type == "Transpose") a$args[[1L]] else a)
  imin <- .ae_default_argmin(lapply(arg$args, keyof))
  if (arg$args[[imin]]$type == "Transpose") {
    arg <- doit(Transpose(arg))
    imin <- .ae_default_argmin(lapply(arg$args, .sort_key))
  }
  n <- length(arg$args)
  Trace(do.call(MatMul, arg$args[c(imin:n, seq_len(imin - 1L))]))
}
.ae_default_argmin <- function(keys) {
  best <- 1L
  for (i in seq_along(keys)[-1L]) if (.key_cmp(keys[[i]], keys[[best]]) < 0) best <- i
  best
}

.a2m <- function(expr) {
  if (!is_sym(expr)) return(expr)
  switch(expr$type,
    ZeroArray = if (get_ndim(expr) == 2L) do.call(ZeroMatrix, .ae_shape(expr)) else expr,
    ArrayTensorProduct = do.call(.a2m_tensor_product, lapply(expr$args, .a2m)),
    ArraySum = .a2m_arraysum(expr),
    ArrayContraction = .a2m_contraction(expr),
    ArrayDiagonal = {
      pexpr <- do.call(.ae_diagonal, c(list(.a2m(expr$args[[1L]])), .ae_diag_ind(expr)))
      pexpr <- .ae_identify_hadamard_products(pexpr)
      if (.ae_is_type(pexpr, "ArrayDiagonal")) pexpr <- .ae_diag2contr_diagmatrix(pexpr)
      if (.ae_eq(expr, pexpr)) return(expr)
      .a2m(pexpr)
    },
    PermuteDims = .a2m_permutedims(expr),
    ArrayAdd = , Add = do.call(.a2m_add, lapply(expr$args, .a2m)),
    ArrayElementwiseApplyFunc = {
      sub <- .a2m(expr$args[[2L]])
      lam <- expr$args[[1L]]
      if (.is_matexpr(sub)) {
        sh <- .mshape(sub)
        if (!(is_one_num(sh[[1L]]) && is_one_num(sh[[2L]]))) {
          m <- .ae_match_power(lam)
          if (!is.null(m)) return(.ae_mmul(m[[1L]], HadamardPower(sub, m[[2L]])))
        }
        return(ElementwiseApplyFunction(lam, sub))
      }
      if (!.ae_has_shape(sub)) return(.ae_lambda_call(lam, sub))
      ArrayElementwiseApplyFunc(lam, sub)
    },
    ArrayElement = {
      ret <- .a2m(expr$args[[1L]])
      if (.is_matexpr(ret)) return(MatrixElement(ret, expr$args[[2L]]$args[[1L]], expr$args[[2L]]$args[[2L]]))
      ArrayElement(ret, expr$args[[2L]]$args)
    },
    expr)
}

# Lambda(d, w*d**p): list(w, p) or NULL
.ae_match_power <- function(lam) {
  d <- lam$args[[1L]]$args[[1L]]
  body <- lam$args[[2L]]
  split <- function(b) {
    if (identical(b$key, d$key)) return(list(S_One, S_One))
    if (b$type == "Pow" && identical(b$args[[1L]]$key, d$key) && !.has1(b$args[[2L]], d)) return(list(S_One, b$args[[2L]]))
    NULL
  }
  r <- split(body)
  if (!is.null(r)) return(r)
  if (body$type == "Mul") {
    hit <- which(vapply(body$args, function(f) .has1(f, d), NA))
    if (length(hit) == 1L) {
      r <- split(body$args[[hit]])
      if (!is.null(r)) return(list(.mul(body$args[-hit]), r[[2L]]))
    }
  }
  NULL
}

.a2m_arraysum <- function(expr) {
  nf <- .a2m_absorb_scalars(.a2m(expr$args[[1L]]))
  out <- do.call(ArraySum, c(list(nf), expr$args[-1L]))
  if (.is_matexpr(nf) && all(vapply(expr$args[-1L], function(l) length(l$args) == 3L &&
      l$args[[2L]]$type == "Integer" && l$args[[3L]]$type == "Integer", NA)))
    return(.ae_arraysum_doit(out, FALSE))
  out
}

.a2m_contraction <- function(expr) {
  expr <- .ae_flatten_contraction_of_diagonal(expr)
  expr <- .ae_identify_removable_identity_matrices(expr)
  if (!.ae_is_type(expr, "ArrayContraction")) return(.a2m(expr))
  expr <- .ae_split_multiple_contractions(expr)
  expr <- .ae_identify_hadamard_products(expr)
  if (!.ae_is_type(expr, "ArrayContraction")) return(.a2m(expr))
  sub <- expr$args[[1L]]
  ci <- .ae_contr_ind(expr)
  ck <- vapply(ci, paste, "", collapse = ",")
  if (get_ndim(sub) == 2L && (identical(ck, c("0", "1")) || identical(ck, "0") || identical(ck, "1"))) {
    shape <- .ae_shape(sub)
    sub <- .a2m(sub)
    if (.is_matexpr(sub)) {
      ret <- sub
      if ("0" %in% ck) ret <- .ae_mmul(OneMatrix(1L, shape[[1L]]), ret)
      if ("1" %in% ck) ret <- .ae_mmul(ret, OneMatrix(shape[[2L]], 1L))
      return(ret)
    }
  }
  if (sub$type == "ArrayTensorProduct") {
    newexpr <- do.call(.ae_contraction, c(list(.a2m(sub)), ci))
    if (!.ae_is_type(newexpr, c("ArrayContraction", "ArrayAdd"))) return(.a2m(newexpr))
    if (newexpr$type == "ArrayContraction") ci <- .ae_contr_ind(newexpr)
    if (newexpr$type == "ArrayContraction" && any(.ae_sub_ndim_list(newexpr) > 2L)) {
      choices <- lapply(sub$args, function(a) if (a$type == "ArrayAdd") a$args else list(a))
      combos <- .ae_product(choices)
      addends <- do.call(.ae_array_add, lapply(combos, function(j) do.call(.a2m_tensor_product, j)))
      newexpr <- do.call(.ae_contraction, c(list(addends), ci))
    }
    if (.ae_is_type(newexpr, "ArrayAdd")) return(.a2m(newexpr))
    if (!.ae_is_type(newexpr, "ArrayContraction")) return(.a2m(newexpr))
    inner <- newexpr$args[[1L]]
    args <- if (inner$type == "ArrayTensorProduct") inner$args else list(inner)
    return(.ae_tp1_recognize(ci, args))
  }
  ret <- .a2m(sub)
  if (.is_matexpr(ret) && identical(ck, "0,1")) return(.a2m_trace(ret))
  do.call(.ae_contraction, c(list(ret), ci))
}

# itertools.product over a list of lists
.ae_product <- function(choices) {
  out <- list(list())
  for (ch in choices) {
    nxt <- list()
    for (o in out) for (c in ch) nxt[[length(nxt) + 1L]] <- c(o, list(c))
    out <- nxt
  }
  out
}

.a2m_permutedims <- function(expr) {
  perm <- .ae_paf(expr$args[[2L]])
  inner <- expr$args[[1L]]
  if (identical(perm, c(1L, 0L))) return(.a2m_transpose(.a2m(inner)))
  if (inner$type == "ArrayTensorProduct") {
    ranks <- .ae_sub_ndim_list(inner)
    invp <- .ae_inv(perm)
    newrange <- invp[seq_len(sum(ranks))]
    cum <- cumsum(c(0L, ranks))
    newargs <- list(); newperm <- integer(); scalars <- list()
    for (i in seq_along(ranks)) {
      pos <- newrange[seq_len(ranks[i]) + cum[i]]
      arg <- inner$args[[i]]
      if (!length(pos)) scalars[[length(scalars) + 1L]] <- .a2m(arg)
      else if (identical(pos, sort(pos))) { newargs[[length(newargs) + 1L]] <- .a2m(arg); newperm <- c(newperm, pos) }
      else if (length(pos) == 2L) { newargs[[length(newargs) + 1L]] <- .a2m_transpose(.a2m(arg)); newperm <- c(newperm, rev(pos)) }
      else .perr("NotImplementedError", "")
    }
    return(.ae_permute_dims(do.call(.a2m_tensor_product, c(scalars, newargs)), .ae_inv(newperm)))
  }
  if (inner$type == "ArrayContraction") {
    lines <- .a2m(inner)
    if (!.ae_is_type(lines, "ArrayTensorProduct")) return(.ae_permute_dims(lines, perm))
    nd <- vapply(lines$args, get_ndim, 0L)
    if (any(!nd %in% c(0L, 2L))) return(.ae_permute_dims(lines, perm))
    scalars <- lines$args[nd == 0L]
    mats <- lines$args[nd == 2L]
    p <- .ae_pmul(.ae_id(2L * length(mats)), perm)
    out <- vector("list", length(mats))
    for (i in seq_along(mats)) {
      p1 <- p[2L * i - 1L]; p2 <- p[2L * i]
      if (p1 %/% 2L != p2 %/% 2L) return(.ae_permute_dims(lines, p))
      out[[i]] <- if (p1 > p2) .a2m_transpose(mats[[p1 %/% 2L + 1L]]) else mats[[p1 %/% 2L + 1L]]
    }
    return(do.call(.a2m_tensor_product, c(scalars, out)))
  }
  expr
}

# ---------------------------------------------------------------------------
# matrix multiplication lines inside contractions

.ae_count_index_slots <- function(ed, index) sum(vapply(ed$args, function(a) sum(!is.na(a$indices) & a$indices == index), 0L))

# isinstance(e, MatrixExpr): an explicit matrix is one too
.ae_ismat <- function(z) .is_matexpr(z) || inherits(z, "symmatrix")

.ae_candidate_for_matmul <- function(scan, remaining) {
  scan <- scan[!is.na(scan)]
  if (!length(scan)) return(list(NULL, FALSE, -1L))
  transpose <- FALSE; cand <- NULL; cidx <- -1L
  for (a2 in remaining) {
    if (!.ae_ismat(a2$element)) next
    for (index in scan) {
      if (cidx != -1L && cidx != index) next
      if (.ae_in(index, a2$indices)) {
        if (!anyNA(a2$indices) && all(a2$indices == index)) { cand <- NULL; break }
        if (is.null(cand)) {
          cand <- a2; cidx <- index
          transpose <- identical(index, a2$indices[2L])
        } else { cand <- NULL; break }
      }
    }
  }
  list(cand, transpose, cidx)
}

.ae_tp1_recognize <- function(ci, args) {
  if (!length(ci)) return(do.call(.a2m_tensor_product, args))
  ac <- do.call(.ae_contraction, c(list(do.call(.ae_tensor_product, args)), ci))
  ed <- .ae_editor(ac)
  .ae_ed_track_start(ed)
  repeat {
    stop_flag <- TRUE
    for (i in seq_along(ed$args)) {
      a <- ed$args[[i]]
      if (!.ae_ismat(a$element)) next
      fi <- a$indices[1L]; si <- a$indices[2L]
      ff <- .ae_ed_count_args_with_index(ed, fi)
      sf <- .ae_ed_count_args_with_index(ed, si)
      if (!is.na(fi) && ff == 1L && identical(fi, si)) {
        stop_flag <- FALSE
        a$element <- .ae_trace_normalize(Trace(a$element))
        a$indices <- integer()
        break
      }
      scan <- integer()
      if (ff == 2L && !is.na(fi) && .ae_count_index_slots(ed, fi) == 2L) scan <- c(scan, fi)
      if (sf == 2L && !is.na(si) && .ae_count_index_slots(ed, si) == 2L) scan <- c(scan, si)
      cr <- .ae_candidate_for_matmul(scan, ed$args[-seq_len(i)])
      cand <- cr[[1L]]
      if (!is.null(cand)) {
        stop_flag <- FALSE
        .ae_ed_track_merge(ed, a, cand)
        t1 <- identical(cr[[3L]], fi)
        other <- cand$element
        if (cr[[2L]]) { other <- Transpose(other); oi <- cand$indices[1L] } else oi <- cand$indices[2L]
        new_el <- .ae_mmul(if (t1) Transpose(a$element) else a$element, other)
        .ae_ed_remove(ed, cand)
        na <- .ae_arge(new_el)
        na$indices <- if (t1) c(si, oi) else c(fi, oi)
        u <- unique(na$indices)
        if (length(u) == 1L && !is.na(u)) {
          rem <- sum(vapply(ed$args, function(b) if (identical(b, a)) 0L else sum(!is.na(b$indices) & b$indices == u), 0L))
          if (rem == 0L) {
            na$element <- .ae_trace_normalize(Trace(na$element))
            na$indices <- integer()
          }
        }
        ed$args[[i]] <- na
        break
      }
    }
    if (stop_flag) break
  }
  .ae_ed_refresh(ed)
  .a2m_absorb_scalars(.ae_ed_to_array_contraction(ed))
}

# ---------------------------------------------------------------------------
# Hadamard products and identity matrices

.ae_check_transpose <- function(x) { x <- ifelse(x >= 0L, x, -1L - x); identical(x, sort(x)) }

.ae_identify_hadamard_products <- function(expr) {
  if (!.ae_is_type(expr, c("ArrayContraction", "ArrayDiagonal", "ArrayTensorProduct"))) return(expr)
  ed <- .ae_editor(expr)
  keys <- character(); groups <- list(); kinds <- list()
  allind <- unlist(lapply(ed$args, function(a) a$indices))
  nind <- function(ind) if (is.na(ind)) sum(is.na(allind)) else sum(allind == ind, na.rm = TRUE)
  for (a in ed$args) {
    if (anyNA(a$indices)) next
    k <- sort(unique(a$indices))
    ks <- paste(k, collapse = ",")
    w <- match(ks, keys)
    if (is.na(w)) { keys <- c(keys, ks); groups[[length(groups) + 1L]] <- list(a); kinds[[length(kinds) + 1L]] <- k }
    else groups[[w]][[length(groups[[w]]) + 1L]] <- a
  }
  for (g in seq_along(keys)) {
    k <- kinds[[g]]; v <- groups[[g]]
    make_trace <- FALSE
    if (length(k) == 1L && k < 0L) {
      ind <- k
      ids <- Filter(function(i) i$element$type == "Identity" && sum(!is.na(i$indices) & i$indices == ind) == 1L && anyNA(i$indices), ed$args)
      total <- sum(vapply(ed$args, function(i) sum(!is.na(i$indices) & i$indices == ind), 0L))
      if (length(ids) == 1L && total == 2L * length(v) + 1L) {
        ida <- ids[[1L]]
        if (all(vapply(v, function(i) .is_matexpr(i$element) && .ae_shape_eq(.ae_shape(i$element), .ae_shape(ida$element)), NA))) {
          ida$element <- do.call(hadamard_product, c(lapply(v, function(i) i$element), list(ida$element)))
          for (i in v) .ae_ed_remove(ed, i)
        }
      }
      next
    }
    if (length(k) == 1L && k >= 0L && sum(vapply(kinds, function(kk) k %in% kk, NA)) == 1L) {
      ind <- k
      if (nind(ind) != 2L * length(v)) {
        if (length(v) > 1L && all(vapply(v, function(i) .is_matexpr(i$element), NA))) {
          hp <- do.call(hadamard_product, lapply(v, function(i) i$element))
          .ae_ed_insert_after(ed, v[[1L]], .ae_arge(hp, c(ind, ind)))
          for (i in v) .ae_ed_remove(ed, i)
        }
        next
      }
      make_trace <- TRUE
      first_element <- S_One
    } else if (length(k) != 2L) next
    if (length(v) == 1L) next
    if (all(vapply(k, function(j) nind(j) == length(v) && j >= 0L, NA))) {
      make_trace <- TRUE
      first_element <- v[[1L]]$element
      if (!.ae_check_transpose(v[[1L]]$indices)) first_element <- .mtranspose(first_element)
      factors <- v[-1L]
    } else factors <- v
    elems <- lapply(factors, function(i) if (.ae_check_transpose(i$indices)) i$element else Transpose(i$element))
    sh <- .ae_shape(elems[[1L]])
    hp <- if (length(sh) == 2L && is_one_num(sh[[1L]]) && is_one_num(sh[[2L]])) do.call(.mm, elems) else do.call(hadamard_product, elems)
    hp_ind <- v[[1L]]$indices
    if (!.ae_check_transpose(factors[[1L]]$indices)) hp_ind <- rev(hp_ind)
    if (make_trace) {
      hp <- .ae_trace_normalize(Trace(.ae_mmul(first_element, .mtranspose(hp))))
      hp_ind <- integer()
    }
    .ae_ed_insert_after(ed, v[[1L]], .ae_arge(hp, hp_ind))
    for (i in v) .ae_ed_remove(ed, i)
  }
  .ae_ed_to_array_contraction(ed)
}

# ask(Q.diagonal(x)) for the cases met here
.ae_is_diagonal <- function(x) {
  if (!.is_matexpr(x)) return(FALSE)
  if (x$type %in% c("Identity", "ZeroMatrix", "DiagMatrix", "DiagonalMatrix")) return(TRUE)
  sh <- .mshape(x)
  is_one_num(sh[[1L]]) && is_one_num(sh[[2L]])
}

.ae_identify_removable_identity_matrices <- function(expr) {
  ed <- .ae_editor(expr)
  flag <- TRUE
  while (flag) {
    flag <- FALSE
    for (a in ed$args) {
      el <- a$element
      if (el$type == "Identity") {
        k <- .mshape(el)[[1L]]
        if (all(is.na(a$indices))) next
        if (anyNA(a$indices)) {
          ind <- a$indices[!is.na(a$indices)][1L]
          cnt <- .ae_ed_count_args_with_index(ed, ind)
          if (cnt == 1L) {
            .ae_ed_insert_after(ed, a, .ae_arge(OneArray(k)))
            .ae_ed_remove(ed, a)
            flag <- TRUE
            break
          }
          next
        }
        if (a$indices[1L] == a$indices[2L]) {
          cnt <- .ae_ed_count_args_with_index(ed, a$indices[1L])
          if (cnt > 1L) { .ae_ed_remove(ed, a); flag <- TRUE; break }
        }
      } else if (.ae_is_diagonal(el)) {
        if (anyNA(a$indices)) next
        if (a$indices[1L] == a$indices[2L]) {
          ind <- a$indices[1L]
          if (.ae_ed_count_args_with_index(ed, ind) == 3L) {
            ind_new <- .ae_ed_new_index(ed)
            others <- Filter(function(j) !identical(j, a), ed$args)
            o <- others[[2L]]
            o$indices <- ifelse(!is.na(o$indices) & o$indices == ind, ind_new, o$indices)
            a$indices <- c(ind, ind_new)
            flag <- TRUE
            break
          }
        }
      }
    }
  }
  .ae_ed_to_array_contraction(ed)
}

.ae_remove_identity_matrices <- function(expr) {
  ed <- .ae_editor(expr)
  removed <- integer()
  pmap_k <- integer(); pmap_v <- integer()
  fi <- cumsum(c(0L, vapply(ed$args, function(a) sum(is.na(a$indices)), 0L)))
  free_of <- function(arg) fi[.ae_ed_pos(ed, arg)]
  free_map <- lapply(ed$args, function(a) list(a, free_of(a)))
  free_lookup <- function(arg) { for (p in free_map) if (identical(p[[1L]], arg)) return(p[[2L]]); NA_integer_ }
  pairs <- list()
  for (ind in seq_len(ed$n_contr) - 1L) {
    args <- .ae_ed_args_with_index(ed, ind)
    isid <- vapply(args, function(i) i$element$type == "Identity", NA)
    nid <- sum(isid)
    if (nid != length(args) - 1L || nid == 0L) next
    nonid <- args[!isid][[1L]]
    ids <- args[isid]
    if (any(vapply(ids, function(i) !anyNA(i$indices), NA))) next
    if (sum(!is.na(nonid$indices) & nonid$indices == ind) != 1L) next
    for (i in ids) {
      i$element <- NULL
      f0 <- free_lookup(i)
      removed <- c(removed, f0 + seq_len(sum(is.na(i$indices))) - 1L)
    }
    last_removed <- removed[length(removed)]
    removed <- removed[-length(removed)]
    pairs[[length(pairs) + 1L]] <- list(last_removed, ind, nonid, nonid$indices)
    nonid$indices <- ifelse(!is.na(nonid$indices) & nonid$indices == ind, NA_integer_, nonid$indices)
  }
  removed <- sort(removed)
  nd <- get_ndim(expr)
  shifts <- cumsum(as.integer((seq_len(nd) - 1L) %in% removed))
  for (p in pairs) {
    freed <- which(!is.na(p[[4L]]) & p[[4L]] == p[[2L]])
    freed <- freed[1L] - 1L
    j <- 0L
    for (a in ed$args) {
      if (is.null(a$element)) next
      if (identical(a, p[[3L]])) { j <- j + sum(is.na(a$indices[seq_len(freed)])); break }
      j <- j + sum(is.na(a$indices))
    }
    pmap_k <- c(pmap_k, j); pmap_v <- c(pmap_v, p[[1L]])
  }
  ed$args <- Filter(function(a) !is.null(a$element), ed$args)
  ret <- .ae_ed_to_array_contraction(ed)
  perm <- integer(); counter <- 0L; counter2 <- 0L
  for (j in seq_len(nd) - 1L) {
    if (j %in% removed) next
    w <- match(counter2, pmap_k)
    if (!is.na(w)) {
      target <- pmap_v[w]
      perm <- c(perm, target - shifts[target + 1L])
      counter2 <- counter2 + 1L
    } else {
      while (counter %in% pmap_v) counter <- counter + 1L
      perm <- c(perm, counter)
      counter <- counter + 1L
      counter2 <- counter2 + 1L
    }
  }
  list(.ae_permute_dims(ret, .ae_inv(perm)), removed)
}

.ae_combine_removed <- function(dim, r1, r2) {
  r1 <- sort(r1); r2 <- sort(r2)
  i <- 0L; j <- 0L; out <- integer()
  repeat {
    if (j >= length(r2)) {
      while (i < length(r1)) { out <- c(out, r1[i + 1L]); i <- i + 1L }
      break
    } else if (i < length(r1) && r1[i + 1L] <= i + r2[j + 1L]) {
      out <- c(out, r1[i + 1L]); i <- i + 1L
    } else {
      out <- c(out, i + r2[j + 1L]); j <- j + 1L
    }
  }
  out
}

.ae_contraction_to_diagonal_multiple_identity <- function(expr) {
  ed <- .ae_editor(expr)
  .ae_ed_track_start(ed)
  dcount <- 0L
  for (i in seq_len(ed$n_contr) - 1L) {
    ids <- list(); args <- list()
    for (a in ed$args) {
      if (!.ae_in(i, a$indices)) next
      if (is_sym(a$element) && a$element$type == "Identity") ids[[length(ids) + 1L]] <- a else args[[length(args) + 1L]] <- a
    }
    if (!length(ids)) next
    if (length(args) + length(ids) < 3L) next
    new_diag <- -1L - dcount
    dcount <- dcount + 1L
    flag <- TRUE
    i1 <- 1L
    for (i1 in seq_along(ids)) {
      id1 <- ids[[i1]]
      if (!anyNA(id1$indices)) { flag <- TRUE; break }
      free_pos <- .ae_ed_abs_free_range(ed, id1)[1L]
      ed$track[[length(ed$track)]] <- c(ed$track[[length(ed$track)]], free_pos)
      id1$element <- NULL
      flag <- FALSE
      break
    }
    if (flag) next
    for (a in c(ids[-i1], args)) a$indices <- ifelse(!is.na(a$indices) & a$indices == i, new_diag, a$indices)
  }
  keep <- vapply(ed$args, function(e) !is.null(e$element), NA)
  ed$track <- c(ed$track[seq_along(ed$args)][keep], ed$track[-seq_along(ed$args)])
  all_vals <- sort(unique(unlist(ed$track)))
  ed$track <- lapply(ed$track, function(t) match(t, all_vals) - 1L)
  ed$args <- ed$args[keep]
  list(.ae_ed_to_array_contraction(ed), integer())
}

.ae_diag2contr_diagmatrix <- function(expr) {
  if (expr$args[[1L]]$type != "ArrayTensorProduct") return(expr)
  args <- expr$args[[1L]]$args
  diag <- .ae_diag_ind(expr)
  mp <- .ae_mapping(vapply(args, .ae_sub_ndim, 0L))
  links <- lapply(diag, function(g) lapply(g, function(j) c(mp$arg[j + 1L], mp$rel[j + 1L])))
  contr <- list()
  total <- get_ndim(expr)
  replaced <- rep(FALSE, length(args))
  converted <- rep(FALSE, length(diag))
  partner_used <- integer()
  for (i in seq_along(diag)) {
    if (length(diag[[i]]) != 2L) next
    p1 <- links[[i]][[1L]]; p2 <- links[[i]][[2L]]
    a1 <- args[[p1[1L] + 1L]]; a2 <- args[[p2[1L] + 1L]]
    if (get_ndim(a1) != 2L || get_ndim(a2) != 2L) next
    s1 <- .ae_shape(a1); s2 <- .ae_shape(a2)
    in1 <- 1L - p1[2L]; in2 <- 1L - p2[2L]
    if (is_one_num(s1[[in1 + 1L]]) && !(p1[1L] %in% partner_used)) {
      d1 <- if (!is_one_num(s1[[p1[2L] + 1L]])) DiagMatrix(a1) else a1
      args[[length(args) + 1L]] <- d1
      contr[[length(contr) + 1L]] <- list(p2, c(length(args) - 1L, p1[2L]))
      total <- total + 1L
      diag[i] <- list(NULL); converted[i] <- TRUE
      args[[p1[1L] + 1L]] <- OneArray(s1[[in1 + 1L]])
      replaced[p1[1L] + 1L] <- TRUE
      partner_used <- c(partner_used, p2[1L])
    } else if (is_one_num(s2[[in2 + 1L]]) && !(p2[1L] %in% partner_used)) {
      d2 <- if (!is_one_num(s2[[p2[2L] + 1L]])) DiagMatrix(a2) else a2
      args[[length(args) + 1L]] <- d2
      contr[[length(contr) + 1L]] <- list(p1, c(length(args) - 1L, p2[2L]))
      total <- total + 1L
      diag[i] <- list(NULL); converted[i] <- TRUE
      args[[p2[1L] + 1L]] <- OneArray(s2[[in2 + 1L]])
      replaced[p2[1L] + 1L] <- TRUE
      partner_used <- c(partner_used, p1[1L])
    }
  }
  if (!any(converted)) return(expr)
  cum <- cumsum(c(0L, vapply(args, get_ndim, 0L)))
  contr2 <- lapply(contr, function(c) vapply(c, function(ab) cum[ab[1L] + 1L] + ab[2L], 0L))
  if (any(replaced)) {
    caxes <- sort(unlist(contr2))
    newd <- list()
    for (i in seq_along(diag)) {
      if (is.null(diag[[i]])) next
      rel <- links[[i]]
      if (any(vapply(rel, function(r) replaced[r[1L] + 1L], NA)) && length(diag[[i]]) != 2L) return(expr)
      np <- vapply(rel, function(r) cum[r[1L] + 1L] + (if (replaced[r[1L] + 1L]) 0L else r[2L]), 0L)
      newd[[length(newd) + 1L]] <- vapply(np, function(p) p - sum(caxes < p), 0L)
    }
  } else newd <- Filter(Negate(is.null), diag)
  tc <- do.call(.ae_contraction, c(list(do.call(.ae_tensor_product, args)), contr2))
  td <- do.call(.ae_diagonal, c(list(tc), newd))
  ng <- length(diag)
  nfree <- get_ndim(expr) - ng
  nconv <- sum(converted)
  perm <- seq_len(nfree) - 1L
  cr <- 0L; sr <- 0L
  newpos <- integer(ng)
  for (i in seq_len(ng)) {
    if (converted[i]) { newpos[i] <- nfree + cr; cr <- cr + 1L }
    else { newpos[i] <- nfree + nconv + sr; sr <- sr + 1L }
  }
  perm <- c(perm, integer(ng))
  for (i in seq_len(ng)) perm[newpos[i] + 1L] <- nfree + i - 1L
  if (!identical(perm, seq_along(perm) - 1L)) td <- .ae_permute_dims(td, perm)
  td
}

# ---------------------------------------------------------------------------
# _remove_trivial_dims: drop axes of dimension one, list(expr, removed axes)

.ae_remove_trivial_dims <- function(expr) {
  if (!is_sym(expr)) return(list(expr, integer()))
  switch(expr$type,
    ArrayTensorProduct = .ae_rtd_tp(expr),
    ArrayAdd = {
      rec <- lapply(expr$args, .ae_remove_trivial_dims)
      newargs <- lapply(rec, `[[`, 1L)
      if (length(unique(vapply(newargs, function(a) .ae_skey(.ae_shape(a)), ""))) > 1L) return(list(expr, integer()))
      list(do.call(.a2m_add, newargs), rec[[1L]][[2L]])
    },
    PermuteDims = {
      r <- .ae_remove_trivial_dims(expr$args[[1L]])
      sub <- r[[1L]]; subrem <- r[[2L]]
      p <- .ae_paf(expr$args[[2L]])
      pinv <- .ae_inv(p)
      shift <- cumsum(as.integer((seq_along(p) - 1L) %in% subrem))
      premoved <- sort(pinv[subrem + 1L])
      p2 <- vapply(p[!p %in% subrem], function(e) e - shift[e + 1L], 0L)
      newexpr <- .ae_permute_dims(sub, p2)
      if (!.ae_eq(newexpr, expr)) {
        r2 <- .ae_remove_trivial_dims(.a2m(newexpr))
        newexpr <- r2[[1L]]
        premoved <- .ae_combine_removed(-1L, premoved, r2[[2L]])
      }
      list(newexpr, premoved)
    },
    ArrayContraction = .ae_rtd_contraction(expr),
    ArrayDiagonal = .ae_rtd_diagonal(expr),
    ElementwiseApplyFunction = {
      r <- .ae_remove_trivial_dims(expr$args[[2L]])
      sub <- r[[1L]]
      if (.is_matexpr(sub) && all(vapply(.mshape(sub), is_one_num, NA)))
        {
        # MatrixElement.__new__ evaluates the entry of an explicit matrix only
        el <- if (.is_explicit(sub)) MatrixElement(sub, 0L, 0L) else new_sym("MatrixElement", list(sub, S_Zero, S_Zero))
        return(list(.ae_lambda_call(expr$args[[1L]], el), c(r[[2L]], 0L, 1L)))
      }
      list(ElementwiseApplyFunction(expr$args[[1L]], sub), integer())
    },
    ArrayElementwiseApplyFunc = {
      r <- .ae_remove_trivial_dims(expr$args[[2L]])
      list(ArrayElementwiseApplyFunc(expr$args[[1L]], r[[1L]]), r[[2L]])
    },
    list(expr, integer()))
}

.ae_is11 <- function(a) { s <- .ae_shape(a); length(s) == 2L && is_one_num(s[[1L]]) && is_one_num(s[[2L]]) }
.ae_has_one <- function(a) any(vapply(.ae_shape(a), is_one_num, NA))

.ae_rtd_tp <- function(expr) {
  removed <- integer(); newargs <- list()
  cum <- cumsum(c(0L, vapply(expr$args, get_ndim, 0L)))
  pending <- NULL; prev_i <- NULL
  for (i in seq_along(expr$args)) {
    arg <- expr$args[[i]]
    rng <- seq.int(cum[i], length.out = cum[i + 1L] - cum[i])
    if (arg$type == "OneArray") { removed <- c(removed, rng); next }
    if (!(.is_matexpr(arg) || inherits(arg, "symmatrix"))) {
      r <- .ae_remove_trivial_dims(arg)
      remd <- rng[r[[2L]] + 1L]
      removed <- c(removed, remd)
      if (!.is_matexpr(r[[1L]]) || !.ae_has_one(r[[1L]])) { newargs[[length(newargs) + 1L]] <- r[[1L]]; next }
      rng <- rng[!rng %in% remd]
      arg <- r[[1L]]
    }
    sh <- .ae_shape(arg)
    if (arg$type == "Identity" && .ae_is11(arg)) { removed <- c(removed, rng); next }
    if (.ae_is11(arg)) {
      r <- .ae_remove_trivial_dims(arg)
      arg <- r[[1L]]
      if (length(r[[2L]])) removed <- c(removed, rng[r[[2L]] + 1L])
      if (!length(newargs)) {
        newargs[[1L]] <- arg; pending <- 1L; prev_i <- i
      } else if (.ae_has_one(newargs[[length(newargs)]])) {
        last <- newargs[[length(newargs)]]
        newargs[[length(newargs)]] <- if (is_one_num(.ae_shape(last)[[2L]])) .ae_mmul(last, arg) else .ae_mmul(arg, last)
        removed <- c(removed, rng)
      } else {
        newargs[[length(newargs) + 1L]] <- arg; pending <- 1L; prev_i <- i
      }
    } else if (any(vapply(sh, is_one_num, NA))) {
      k <- Filter(function(s) !is_one_num(s), sh)[[1L]]
      if (is.null(pending)) {
        pending <- k; prev_i <- i
        newargs[[length(newargs) + 1L]] <- arg
      } else if (is_sym(pending) && is_one_num(pending) || identical(pending, 1L)) {
        prev <- newargs[[length(newargs)]]
        newargs[[length(newargs)]] <- if (is_one_num(sh[[1L]])) .ae_mmul(prev, arg) else .ae_mmul(arg, prev)
        removed <- c(removed, setdiff(c(cum[prev_i], cum[prev_i] + 1L), removed))
        pending <- k; prev_i <- i
      } else if (identical(pending$key, k$key)) {
        prev <- newargs[[length(newargs)]]
        if (is_one_num(.ae_shape(prev)[[1L]])) { d1 <- cum[prev_i]; prev <- .a2m_transpose(prev) } else d1 <- cum[prev_i] + 1L
        if (is_one_num(sh[[2L]])) { d2 <- cum[i] + 1L; arg <- .a2m_transpose(arg) } else d2 <- cum[i]
        newargs[[length(newargs)]] <- .ae_mmul(prev, arg)
        pending <- NULL
        removed <- c(removed, d1, d2)
      } else {
        newargs[[length(newargs) + 1L]] <- arg
        pending <- k; prev_i <- i
      }
    } else {
      newargs[[length(newargs) + 1L]] <- arg
      pending <- NULL
    }
  }
  newexpr <- do.call(.a2m_tensor_product, newargs)
  newremoved <- sort(unique(removed))
  if (.ae_is_type(newexpr, "ArrayTensorProduct")) {
    r <- .ae_find_trivial_matrices_rewrite(newexpr)
    newexpr <- r[[1L]]; newremoved <- .ae_combine_removed(-1L, newremoved, r[[2L]])
  }
  if (.ae_is_type(newexpr, "ArrayTensorProduct")) {
    r <- .ae_find_trivial_kronecker_broadcast(newexpr)
    newexpr <- r[[1L]]; newremoved <- .ae_combine_removed(-1L, newremoved, r[[2L]])
  }
  list(newexpr, newremoved)
}

.ae_find_trivial_matrices_rewrite <- function(expr) {
  triv <- list(); pos <- NULL; first <- NULL; second <- NULL
  removed <- integer(); counter <- 0L
  args <- expr$args
  keep <- rep(TRUE, length(args))
  for (i in seq_along(expr$args)) {
    arg <- expr$args[[i]]
    if (.is_matexpr(arg)) {
      if (.ae_is11(arg)) {
        triv[[length(triv) + 1L]] <- arg
        keep[i] <- FALSE
        removed <- c(removed, counter, counter + 1L)
      } else if (is.null(pos) && arg$type == "MatMul") {
        m <- arg$args
        for (j in seq_along(m)) {
          e <- m[[j]]
          if (.is_matexpr(e) && is_one_num(.mshape(e)[[2L]])) {
            pos <- i
            first <- do.call(MatMul, m[seq_len(j)])
            second <- do.call(MatMul, m[-seq_len(j)])
            break
          }
        }
      }
    }
    counter <- counter + get_ndim(arg)
  }
  if (is.null(pos)) return(list(expr, integer()))
  args[[pos]] <- .ae_mmul(first, do.call(MatMul, triv), second)
  list(do.call(.ae_tensor_product, args[keep]), removed)
}

.ae_find_trivial_kronecker_broadcast <- function(expr) {
  newargs <- list(); removed <- integer(); count <- 0L
  for (arg in expr$args) {
    count <- count + get_ndim(arg)
    sh <- .ae_shape(arg)
    cur <- count - rev(seq_along(sh))
    last <- if (length(newargs)) newargs[[length(newargs)]] else NULL
    if (.ae_is11(arg) && !is.null(last) && !.ae_has_one(last) && .is_matexpr(last) && .is_matexpr(arg)) {
      newargs[[length(newargs)]] <- KroneckerProduct(last, arg)
      removed <- c(removed, cur)
    } else if (!any(vapply(sh, is_one_num, NA)) && !is.null(last) && .ae_is11(last)) {
      newargs[[length(newargs)]] <- KroneckerProduct(last, arg)
      prev <- setdiff(seq_len(min(cur)) - 1L, removed)
      removed <- c(removed, utils::tail(prev, 2L))
    } else newargs[[length(newargs) + 1L]] <- arg
  }
  list(do.call(.ae_tensor_product, newargs), removed)
}

.ae_rtd_contraction <- function(expr) {
  r0 <- .ae_contraction_to_diagonal_multiple_identity(expr)
  if (!.ae_eq(r0[[1L]], expr)) {
    r1 <- .ae_remove_trivial_dims(.a2m(r0[[1L]]))
    return(list(r1[[1L]], .ae_combine_removed(-1L, r0[[2L]], r1[[2L]])))
  }
  rank1 <- get_ndim(expr)
  ri <- .ae_remove_identity_matrices(expr)
  expr <- ri[[1L]]; removed1 <- ri[[2L]]
  if (!.ae_is_type(expr, "ArrayContraction")) {
    r2 <- .ae_remove_trivial_dims(expr)
    return(list(r2[[1L]], .ae_combine_removed(rank1, removed1, r2[[2L]])))
  }
  r2 <- .ae_remove_trivial_dims(expr$args[[1L]])
  newexpr <- r2[[1L]]; removed2 <- r2[[2L]]
  ci <- .ae_contr_ind(expr)
  shifts <- cumsum(as.integer((seq_len(get_ndim(expr$args[[1L]])) - 1L) %in% removed2))
  nci <- lapply(ci, function(g) g[!g %in% removed2])
  nci <- Filter(length, nci)
  removed2 <- removed2[!removed2 %in% unlist(ci)]
  nci <- lapply(nci, function(g) g - shifts[g + 1L])
  removed2 <- unlist(.ae_contr_push_up(ci, as.list(removed2)))
  removed <- .ae_combine_removed(rank1, removed1, if (is.null(removed2)) integer() else removed2)
  list(do.call(.ae_contraction, c(list(newexpr), nci)), removed)
}

.ae_remove_diagonalized_identity_matrices <- function(expr) {
  ed <- .ae_editor(expr)
  nd <- .ae_ed_n_diag(ed)
  removed <- integer()
  for (a in ed$args) {
    if (is.null(a$element) || a$element$type != "Identity") next
    if (!(anyNA(a$indices) && any(!is.na(a$indices) & a$indices < 0L))) next
    dind <- a$indices[!is.na(a$indices)][1L]
    group <- Filter(function(j) .ae_in(dind, j$indices), ed$args)
    other <- Filter(function(j) !identical(j, a), group)
    if (!length(other)) next
    other <- other[[1L]]
    if (!.is_matexpr(other$element)) next
    if (!.ae_has_one(other$element)) next
    if (!anyNA(other$indices)) next
    a$element <- NULL
    none_index <- which(is.na(other$indices))[1L] - 1L
    other$element <- DiagMatrix(other$element)
    rng <- .ae_ed_abs_range(ed, other)
    removed <- c(removed, rng[1L] + none_index)
  }
  ed$args <- Filter(function(a) !is.null(a$element), ed$args)
  removed <- unlist(.ae_diag_push_up(.ae_diag_ind(expr), as.list(removed), get_ndim(expr$args[[1L]])))
  list(.ae_ed_to_array_contraction(ed), if (is.null(removed)) integer() else removed)
}

.ae_rtd_diagonal <- function(expr) {
  r <- .ae_remove_trivial_dims(expr$args[[1L]])
  newexpr <- r[[1L]]; removed <- r[[2L]]
  nd <- get_ndim(expr$args[[1L]])
  shifts <- cumsum(c(0L, as.integer((seq_len(nd) - 1L) %in% removed)))
  diag <- .ae_diag_ind(expr)
  newmap <- lapply(diag, function(g) g[!g %in% removed])
  for (i in seq_along(diag)) if (length(newmap[[i]]) == 1L) removed <- removed[!removed %in% diag[[i]]]
  new_diag <- lapply(newmap, function(g) g - shifts[g + 1L])
  removed <- unlist(.ae_diag_push_up(diag, as.list(removed), nd))
  removed <- sort(unique(if (is.null(removed)) integer() else removed))
  new_diag <- Filter(length, new_diag)
  newexpr2 <- if (length(new_diag)) do.call(.ae_diagonal, c(list(newexpr), new_diag, list(allow_trivial_diags = TRUE))) else newexpr
  if (.ae_is_type(newexpr2, "ArrayDiagonal")) {
    r3 <- .ae_remove_diagonalized_identity_matrices(newexpr2)
    return(list(r3[[1L]], .ae_combine_removed(-1L, removed, r3[[2L]])))
  }
  list(newexpr2, removed)
}
