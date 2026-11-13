# Explicit forms of array expressions: the as_explicit() methods of
# sympy/tensor/array/expressions/array_expressions.py, and the elementwise
# arithmetic of explicit arrays (NDimArray.__add__ and friends).

# arg.as_explicit() if hasattr(arg, "as_explicit") else arg
.ae_explicit_arg <- function(a) {
  if (!is_sym(a)) return(a)
  if (a$type == "ImmutableDenseNDimArray") return(.ae_ndarray_r(a))
  if (.is_arrexpr(a)) return(.ae_as_explicit(a))
  if (.is_matexpr(a)) return(as_explicit(a))
  a
}

# _arrayfy: an explicit matrix becomes a rank 2 array
.ae_arrayfy <- function(a) {
  if (inherits(a, "symmatrix")) return(Array(a$data, c(a$nrow, a$ncol)))
  if (is_sym(a) && a$type == "ImmutableDenseNDimArray") return(.ae_ndarray_r(a))
  a
}

.ae_int_shape <- function(sh) all(vapply(sh, function(s) is_sym(s) && s$type == "Integer", NA))
.ae_int_dims <- function(sh) vapply(sh, function(s) as.integer(.int_val(s)), 0L)

# all index tuples of a shape, one per row, 1 based
.ae_index_rows <- function(dims) if (!length(dims)) matrix(integer(), 1L, 0L) else .arr_indices(dims)

# tensorproduct() on explicit operands
.ae_np_tensorproduct <- function(args) {
  args <- lapply(args, .ae_arrayfy)
  if (!any(vapply(args, inherits, NA, "SymArray"))) return(.mul(args))
  do.call(tensorproduct, args)
}

# tensorcontraction(array, *groups) with 0 based axis groups
.ae_np_contraction <- function(a, groups) {
  a <- .ae_arrayfy(a)
  sh <- a$shape
  for (g in groups) if (length(unique(sh[g + 1L])) != 1L)
    .perr("ValueError", "contracting indices of different dimensions")
  rest <- setdiff(seq_along(sh), unlist(groups) + 1L)
  outer <- .ae_index_rows(sh[rest])
  inner <- .ae_index_rows(vapply(groups, function(g) sh[g[1L] + 1L], 0L))
  vals <- lapply(seq_len(nrow(outer)), function(r) .add(lapply(seq_len(nrow(inner)), function(s) {
    idx <- integer(length(sh))
    idx[rest] <- outer[r, ]
    for (k in seq_along(groups)) idx[groups[[k]] + 1L] <- inner[s, k]
    arr_get(a, idx)
  })))
  if (!length(rest)) return(vals[[1L]])
  Array(vals, sh[rest])
}

# expr.applyfunc(f) for an explicit matrix or array
.ae_np_applyfunc <- function(x, f) {
  if (inherits(x, "SymArray")) return(Array(lapply(x$data, f), x$shape))
  if (inherits(x, "symmatrix")) return(.new_mat(lapply(x$data, f), x$nrow, x$ncol))
  f(x)
}

.ae_as_explicit <- function(e) {
  a <- e$args
  sh <- .ae_shape(e)
  switch(e$type,
    ArraySymbol = {
      if (!.ae_int_shape(sh)) .perr("ValueError", "cannot express explicit array with symbolic shape")
      dims <- .ae_int_dims(sh)
      ind <- .ae_index_rows(dims) - 1L
      Array(lapply(seq_len(nrow(ind)), function(r) ArrayElement(e, lapply(ind[r, ], .int))), dims)
    },
    ZeroArray = , OneArray = {
      if (!.ae_int_shape(sh)) .perr("ValueError", "Cannot return explicit form for symbolic shape.")
      dims <- .ae_int_dims(sh)
      Array(rep(list(if (e$type == "ZeroArray") S_Zero else S_One), prod(dims)), dims)
    },
    ArrayTensorProduct = .ae_np_tensorproduct(lapply(a, .ae_explicit_arg)),
    ArrayAdd = {
      terms <- lapply(a, .ae_explicit_arg)
      arr <- vapply(terms, inherits, NA, "SymArray")
      # mixing matrices and N-dim arrays: everything becomes an N-dim array
      if (any(arr) && !all(arr))
        terms <- lapply(terms, function(t) if (inherits(t, "SymArray")) t else if (inherits(t, "symmatrix")) .ae_arrayfy(t) else Array(list(t), integer()))
      Reduce(`+`, terms)
    },
    PermuteDims = {
      x <- .ae_arrayfy(.ae_explicit_arg(a[[1L]]))
      p <- .ae_paf(a[[2L]])
      n <- length(x$shape)
      if (length(p) < n) p <- c(p, seq.int(length(p), n - 1L))
      permutedims(x, p + 1L)
    },
    ArrayDiagonal = {
      groups <- lapply(a[-1L], .ae_v)
      if (any(lengths(groups) == 1L)) {
        # tensordiagonal() takes no trivial groups; the canonical form is equivalent
        r <- .ae_diag_canon(e)
        return(.ae_explicit_arg(r))
      }
      do.call(tensordiagonal, c(list(.ae_arrayfy(.ae_explicit_arg(a[[1L]]))), lapply(groups, function(g) g + 1L)))
    },
    ArrayContraction = .ae_np_contraction(.ae_explicit_arg(a[[1L]]), lapply(a[-1L], .ae_v)),
    ArrayElementwiseApplyFunc = {
      x <- .ae_explicit_arg(a[[2L]])
      f <- function(z) .ae_lambda_call(a[[1L]], z)
      if (!length(.ae_shape(a[[2L]]))) {
        if (inherits(x, "SymArray")) x <- x$data[[1L]]
        return(f(x))
      }
      .ae_np_applyfunc(x, f)
    },
    Reshape = {
      x <- .ae_explicit_arg(a[[1L]])
      if (inherits(x, "symmatrix")) x <- .ae_arrayfy(x)
      else if (.is_matexpr(x)) return(e)
      dims <- .ae_int_dims(a[[2L]]$args)
      if (prod(dims) != length(x$data)) .perr("ValueError", paste0("Invalid reshape parameters ", .str(a[[2L]])))
      Array(x$data, dims)
    },
    ArraySum = {
      lims <- a[-1L]
      if (!all(vapply(lims, function(l) length(l$args) == 3L && l$args[[2L]]$type == "Integer" && l$args[[3L]]$type == "Integer", NA)))
        .perr("ValueError", "cannot express explicit form of array sum with symbolic limits")
      if (!.ae_int_shape(sh)) .perr("ValueError", "cannot express explicit array with symbolic shape")
      lo <- vapply(lims, function(l) .int_val(l$args[[2L]]), 0)
      hi <- vapply(lims, function(l) .int_val(l$args[[3L]]), 0)
      if (any(hi < lo)) return(Array(rep(list(S_Zero), prod(.ae_int_dims(sh))), .ae_int_dims(sh)))
      grid <- as.matrix(expand.grid(rev(lapply(seq_along(lo), function(i) lo[i]:hi[i]))))[, rev(seq_along(lo)), drop = FALSE]
      terms <- lapply(seq_len(nrow(grid)), function(r) {
        t <- subs(a[[1L]], lapply(seq_along(lims), function(i) list(lims[[i]]$args[[1L]], .int(grid[r, i]))))
        t <- .ae_arrayfy(.ae_explicit_arg(t))
        if (inherits(t, "SymArray")) t else Array(list(t), integer())
      })
      Reduce(`+`, terms)
    },
    .perr("NotImplementedError", paste0("as_explicit() is not implemented for ", e$type)))
}

# NDimArray arithmetic: elementwise sums, products with scalars

#' @export
chooseOpsMethod.SymArray <- function(x, y, mx, my, cl, reverse) TRUE

#' @export
Ops.SymArray <- function(e1, e2) {
  isA <- function(z) inherits(z, "SymArray")
  if (missing(e2)) {
    if (.Generic == "-") return(Array(lapply(e1$data, .neg), e1$shape))
    if (.Generic == "+") return(e1)
  }
  if (.Generic %in% c("==", "!=")) {
    same <- isA(e1) && isA(e2) && identical(e1$shape, e2$shape) &&
      identical(vapply(e1$data, .k, ""), vapply(e2$data, .k, ""))
    return(if (.Generic == "==") same else !same)
  }
  if (.Generic %in% c("+", "-")) {
    if (!isA(e1) || !isA(e2)) .perr("TypeError", paste0("unsupported operand types for ", .Generic))
    if (!identical(e1$shape, e2$shape)) .perr("ValueError", "array shape mismatch")
    f <- if (.Generic == "+") function(x, y) .add(list(x, y)) else .sub
    return(Array(Map(f, e1$data, e2$data), e1$shape))
  }
  if (.Generic %in% c("*", "/")) {
    if (isA(e1) && isA(e2)) .perr("ValueError", "scalar expected, use tensorproduct(...) for tensorial product")
    if (.Generic == "/" && !isA(e1)) .perr("TypeError", "unsupported operand types for /")
    if (isA(e1)) {
      s <- as_sym(e2)
      f <- if (.Generic == "*") function(x) .mul(list(x, s)) else function(x) .div(x, s)
      return(Array(lapply(e1$data, f), e1$shape))
    }
    s <- as_sym(e1)
    return(Array(lapply(e2$data, function(x) .mul(list(s, x))), e2$shape))
  }
  .perr("TypeError", paste0("unsupported operand types for ", .Generic))
}
