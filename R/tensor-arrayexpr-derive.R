# Derivatives of matrix and array expressions, after
# sympy/tensor/array/expressions/arrayexpr_derivatives.py (array_derive),
# sympy/tensor/array/array_derivatives.py (ArrayDerivative) and the
# _eval_derivative methods of the matrix expression classes.

# free symbols in SymPy's sense: matrix and array symbols are atoms
.ae_free_keys <- function(e) {
  if (inherits(e, "symmatrix")) return(unique(unlist(lapply(e$data, .ae_free_keys))))
  if (inherits(e, "SymArray")) return(unique(unlist(lapply(e$data, .ae_free_keys))))
  if (!is_sym(e)) return(character())
  if (e$type %in% c("ArraySymbol", "MatrixSymbol", "Symbol", "Dummy")) return(e$key)
  if (e$type == "Lambda") return(setdiff(.ae_free_keys(e$args[[2L]]), vapply(e$args[[1L]]$args, .k, "")))
  if (e$type %in% c("Sum", "ArraySum")) {
    bound <- vapply(e$args[-1L], function(l) l$args[[1L]]$key, "")
    inner <- setdiff(.ae_free_keys(e$args[[1L]]), bound)
    return(unique(c(inner, unlist(lapply(e$args[-1L], function(l) unlist(lapply(l$args[-1L], .ae_free_keys)))))))
  }
  if (!length(e$args)) return(character())
  unique(unlist(lapply(e$args, .ae_free_keys)))
}

.ae_is_explicit <- function(e) inherits(e, "symmatrix") || .ae_is_type(e, "ImmutableMatrix")
.ae_is_ndarray <- function(e) inherits(e, "SymArray") || .ae_is_type(e, "ImmutableDenseNDimArray")
.ae_as_mat <- function(e) if (inherits(e, "symmatrix")) e else as_explicit(e)
.ae_as_arr <- function(e) {
  if (inherits(e, "SymArray")) return(e)
  if (.ae_is_type(e, "ImmutableDenseNDimArray")) return(.ae_ndarray_r(e))
  M <- .ae_as_mat(e)
  Array(M$data, c(M$nrow, M$ncol))
}
.ae_is_matrix_like <- function(e) .is_matexpr(e) || inherits(e, "symmatrix")
.ae_is_scalar <- function(e) !(.ae_is_matrix_like(e) || .ae_is_ndarray(e) || .ae_is_type(e, c("ImmutableMatrix")) ||
                                 (.ae_is_type(e, "Derivative") && .ae_derivative_is_array(e)))
.ae_is_zero <- function(e) is_sym(e) && (is_zero_num(e) || e$type == "ZeroMatrix")

.ae_derivative_is_array <- function(e) {
  parts <- c(list(e$args[[1L]]), lapply(e$args[-1L], function(t) t$args[[1L]]))
  any(vapply(parts, function(p) .ae_is_matrix_like(p) || .ae_is_ndarray(p) || .ae_is_type(p, "ImmutableMatrix"), NA))
}
.ae_derivative_shape <- function(e) {
  shape <- list()
  for (t in e$args[-1L]) {
    v <- t$args[[1L]]
    sh <- .ae_shape(v)
    if (length(sh)) for (i in seq_len(.int_val(t$args[[2L]]))) shape <- c(shape, sh)
  }
  c(shape, .ae_shape(e$args[[1L]]))
}

# d expr / d x for any kind of expression and variable
.ae_diff <- function(e, x) .diff_spec_eval(as_sym_any(e), list(list(x, 1L)))
as_sym_any <- function(e) if (is.numeric(e)) as_sym(e) else e

# Does diff() of f by the variables of spec go through ArrayDerivative?
.ae_derivative_applies <- function(f, spec) {
  arr <- function(z) .ae_is_matrix_like(z) || .ae_is_ndarray(z) || .ae_is_type(z, "ImmutableMatrix")
  arr(f) || any(vapply(spec, function(s) arr(s[[1L]]) || .ae_is_type(s[[1L]], "MatrixElement"), NA))
}

.ae_zero_like <- function(expr) {
  if (.ae_is_explicit(expr)) { M <- .ae_as_mat(expr); return(zeros(M$nrow, M$ncol)) }
  if (.ae_is_ndarray(expr)) { a <- .ae_as_arr(expr); return(Array(rep(list(S_Zero), length(a$data)), a$shape)) }
  if (.is_matexpr(expr)) return(do.call(ZeroMatrix, .mshape(expr)))
  .perr("RuntimeError", "Unable to determine shape of array-derivative.")
}

# Derivative.__new__ with evaluate=True, for ArrayDerivative
.ae_derivative_eval <- function(f, spec) {
  free <- .ae_free_keys(f)
  for (s in spec) {
    v <- s[[1L]]; c <- s[[2L]]
    vfree <- .ae_free_keys(v)
    if (!(is.numeric(c) && c > 0) || !length(vfree)) next
    if (.is_matexpr(v) || .ae_is_type(v, "MatrixElement")) break
    if (.ae_is_type(v, "Symbol") && !(v$key %in% free)) return(.ae_zero_like(f))
    if (!length(intersect(free, vfree))) return(.ae_zero_like(f))
  }
  expr <- f
  for (i in seq_along(spec)) {
    v <- spec[[i]][[1L]]; cnt <- spec[[i]][[2L]]
    obj <- .ae_dispatch_n(expr, v, cnt)
    if (!is.null(obj) && .ae_is_zero(obj)) return(obj)
    if (is.null(obj)) {
      if (inherits(expr, "symmatrix")) expr <- .as_block(expr)
      return(.derivative_node(expr, spec[i:length(spec)], sort = FALSE))
    }
    expr <- obj
  }
  expr
}

.ae_applyfunc <- function(M, f) {
  if (.ae_is_ndarray(M)) { a <- .ae_as_arr(M); return(Array(lapply(a$data, f), a$shape)) }
  M <- .ae_as_mat(M)
  .new_mat(lapply(M$data, f), M$nrow, M$ncol)
}

.ae_dispatch_n <- function(expr, v, count) {
  if (!is.numeric(count) || count <= 0) return(NULL)
  if (.ae_is_scalar(expr)) {
    if (.ae_is_explicit(v)) result <- .ae_applyfunc(v, function(z) .ae_diff(expr, z))
    else if (.is_matexpr(v)) result <- if (.ae_has(expr, v)) .ae_matrix_derivative(expr, v) else do.call(ZeroMatrix, .mshape(v))
    else if (.ae_is_ndarray(v)) result <- .ae_applyfunc(v, function(z) .ae_diff(expr, z))
    else if (.ae_is_type(v, "MatrixElement")) {
      result <- expr
      for (k in seq_len(count)) result <- .me_diff_by_elem(result, v)
      return(result)
    }
    else if (.ae_is_scalar(v)) return(.diff_count(expr, v, count))
    else return(NULL)
  } else if (.ae_is_scalar(v)) {
    if (.ae_is_explicit(expr)) result <- .ae_applyfunc(expr, function(z) .ae_diff(z, v))
    else if (.is_matexpr(expr)) result <- .me_eval_derivative(expr, v)
    else if (.ae_is_ndarray(expr)) result <- .ae_applyfunc(expr, function(z) .ae_diff(z, v))
    else return(NULL)
  } else {
    if (.ae_is_explicit(expr) && .ae_is_explicit(v)) result <- derive_by_array(.ae_as_arr(expr), .ae_as_arr(v))
    else if (.ae_is_explicit(expr) && .is_matexpr(v)) result <- .ae_applyfunc(expr, function(z) .ae_diff(z, v))
    else if (.ae_is_explicit(expr) || .ae_is_explicit(v)) {
      # derive_by_array leaves a symbolic matrix expression whole
      if (.is_matexpr(expr) && !.ae_is_explicit(expr)) {
        a <- .ae_as_arr(v)
        result <- Array(lapply(a$data, function(z) .ae_diff(expr, z)), a$shape)
      } else result <- derive_by_array(.ae_as_arr(expr), .ae_as_arr(v))
    }
    else if (.is_matexpr(expr) && .is_matexpr(v)) result <- if (.ae_has(expr, v)) .ae_matrix_derivative(expr, v) else NULL
    else if (.is_matexpr(expr) || .is_matexpr(v)) return(NULL)
    else result <- derive_by_array(.ae_as_arr(expr), .ae_as_arr(v))
  }
  if (is.null(result)) return(NULL)
  if (count == 1) result else .ae_dispatch_n(result, v, count - 1L)
}

.ae_has <- function(e, x) {
  if (inherits(e, c("symmatrix", "SymArray"))) return(any(vapply(e$data, .ae_has, NA, x)))
  .has1(e, x)
}

# _matrix_derivative
.ae_matrix_derivative <- function(expr, x) {
  if (.ae_is_explicit(expr) || .ae_is_explicit(x)) {
    if (.ae_is_explicit(expr) && .ae_is_explicit(x)) return(derive_by_array(.ae_as_arr(expr), .ae_as_arr(x)))
    if (.ae_is_explicit(expr)) return(.ae_applyfunc(expr, function(z) .ae_diff(z, x)))
    return(.ae_applyfunc(x, function(z) .ae_diff(expr, z)))
  }
  convert_array_to_matrix(.ad(.m2a(expr), x))
}

# ---------------------------------------------------------------------------
# _eval_derivative of matrix expressions by a scalar; NULL when not computed

.me_eval_derivative <- function(e, x) {
  # a type with its own _eval_derivative, such as codegen's MatrixSolve
  h <- .diff_hooks[[e$type]]
  if (!is.null(h)) return(h(e, x))
  has_x <- .ae_has(e, x) || (.ae_is_type(x, "MatrixElement") && .ae_has(e, x$args[[1L]]))
  if (e$type == "MatrixSymbol") {
    if (length(intersect(.ae_free_keys(e), .ae_free_keys(x)))) {
      if (.ae_is_type(x, "MatrixElement") && identical(x$args[[1L]]$key, e$key))
        return(MatrixUnit(.mshape(e)[[1L]], .mshape(e)[[2L]], x$args[[2L]], x$args[[3L]]))
      return(NULL)
    }
    return(do.call(ZeroMatrix, .mshape(e)))
  }
  if (e$type == "Transpose") {
    d <- .me_eval_derivative(e$args[[1L]], x)
    return(if (is.null(d)) NULL else .mtranspose(d))
  }
  if (e$type == "MatPow" || e$type == "Inverse") {
    base <- e$args[[1L]]
    ex <- if (e$type == "Inverse") S_NegativeOne else e$args[[2L]]
    if (!is_zero_num(.ae_diff(ex, x))) return(NULL)
    if (ex$type == "Integer") {
      if (.int_val(ex) > 0) return(.me_eval_derivative(do.call(MatMul, rep(list(base), .int_val(ex))), x))
      if (.int_val(ex) == 0) return(do.call(ZeroMatrix, .mshape(base)))
    }
    return(NULL)
  }
  if (e$type == "HadamardPower") {
    base <- e$args[[1L]]; ex <- e$args[[2L]]
    dexp <- .ae_diff(ex, x)
    logbase <- ElementwiseApplyFunction(function(z) fcall("log", list(z)), base)
    dlbase <- .ae_diff(logbase, x)
    return(hadamard_product(.ae_madd(.ae_mmul(dexp, logbase), .ae_mmul(ex, dlbase)), e))
  }
  if (e$type == "ElementwiseApplyFunction") {
    dexpr <- .ae_diff(e$args[[2L]], x)
    fd <- .ae_function_fdiff(e$args[[1L]])
    return(hadamard_product(dexpr, ElementwiseApplyFunction(fd, e$args[[2L]])))
  }
  if (e$type == "HadamardProduct") {
    args <- e$args
    terms <- lapply(seq_along(args), function(i) {
      f <- args; f[[i]] <- .ae_diff(args[[i]], x)
      do.call(hadamard_product, f)
    })
    return(do.call(.ae_madd, terms))
  }
  if (!has_x) return(do.call(ZeroMatrix, .mshape(e)))
  if (e$type == "MatAdd") return(doit(do.call(MatAdd, lapply(e$args, .ae_diff, x))))
  if (e$type == "MatMul") {
    args <- e$args
    terms <- list()
    for (i in seq_along(args)) {
      d <- .ae_diff(args[[i]], x)
      if (is_sym(d) && is_zero_num(d)) next
      f <- args; f[[i]] <- d
      terms[[length(terms) + 1L]] <- do.call(.ae_mmul, f)
    }
    return(do.call(.ae_madd, terms))
  }
  NULL
}

# Add.fromiter of matrix terms (the canonical MatAdd), scalars added up
.ae_madd <- function(...) {
  ts <- Filter(function(t) !(is_sym(t) && is_zero_num(t)), list(...))
  if (!length(ts)) return(S_Zero)
  if (length(ts) == 1L) return(ts[[1L]])
  if (any(vapply(ts, .is_matexpr, NA))) return(doit(do.call(MatAdd, ts)))
  .add(ts)
}

# ---------------------------------------------------------------------------
# array_derive

#' Derivatives of array and matrix expressions
#'
#' `array_derive(expr, x)` differentiates an array expression by an array or
#' matrix symbol, giving an array expression whose first axes are those of
#' `x`. `matrix_derive(expr, x)` differentiates a matrix expression and
#' converts the result back to a matrix expression where possible. These
#' follow `sympy.tensor.array.expressions.arrayexpr_derivatives`; [diff()]
#' uses them for matrix expressions.
#'
#' @param expr An array or matrix expression.
#' @param x A matrix or array symbol.
#' @return An array or matrix expression.
#' @examples
#' k <- Symbol("k")
#' X <- MatrixSymbol("X", k, k); A <- MatrixSymbol("A", k, k)
#' array_derive(convert_matrix_to_array(A * X), X)
#' matrix_derive(Trace(A * X), X)
#' @export
array_derive <- function(expr, x) .ad(.ae_sympify(expr), .ae_sympify(x))

#' @rdname array_derive
#' @export
matrix_derive <- function(expr, x) doit(convert_array_to_matrix(.ad(.m2a(.ae_sympify(expr)), .ae_sympify(x))))

.ae_shapes <- function(...) do.call(c, c(list(list()), lapply(list(...), .ae_shape)))

.ad <- function(expr, x) {
  if (.ae_is_explicit(expr)) return(.ad_matrixbase(expr, x))
  if (.ae_is_ndarray(expr)) return(.ad_ndarray(expr, x))
  rank_x <- length(.ae_shape(x))
  switch(expr$type,
    ArrayTensorProduct = .ad_tp(expr, x),
    ArraySymbol = {
      if (.ae_eq(expr, x)) {
        sh <- .ae_shape(expr); n <- length(sh)
        return(.ae_permute_dims(do.call(ArrayTensorProduct, lapply(sh, Identity)), c(2L * (seq_len(n) - 1L), 2L * (seq_len(n) - 1L) + 1L)))
      }
      do.call(ZeroArray, .ae_shapes(x, expr))
    },
    MatrixSymbol = {
      if (.ae_eq(expr, x)) {
        sh <- .mshape(expr)
        return(.ae_permute_dims(.ae_tensor_product(Identity(sh[[1L]]), Identity(sh[[2L]])), c(0L, 2L, 1L, 3L)))
      }
      do.call(ZeroArray, .ae_shapes(x, expr))
    },
    Determinant = {
      arg <- expr$args[[1L]]
      darg <- .ad(arg, x)
      tp <- .ae_tensor_product(expr, .ae_minv(arg), darg)
      .ae_contraction(tp, c(0L, rank_x + 3L), c(1L, rank_x + 2L))
    },
    Identity = , OneMatrix = , ZeroMatrix = , ZeroArray = do.call(ZeroArray, .ae_shapes(x, expr)),
    Transpose = .ae_permute_dims(.ad(expr$args[[1L]], x), c(seq_len(rank_x) - 1L, rank_x + 1L, rank_x)),
    Inverse = {
      mat <- expr$args[[1L]]
      dexpr <- .ad(mat, x)
      tp <- .ae_tensor_product(.ae_mmul(S_NegativeOne, expr), dexpr, expr)
      mp <- .ae_contraction(tp, c(1L, rank_x + 2L), c(rank_x + 3L, rank_x + 4L))
      .ae_permute_dims(mp, c(seq_len(rank_x), 0L, rank_x + 1L))
    },
    ElementwiseApplyFunction = {
      fd <- .ae_function_fdiff(expr$args[[1L]])
      dexpr <- .ad(expr$args[[2L]], x)
      tp <- .ae_tensor_product(ElementwiseApplyFunction(fd, expr$args[[2L]]), dexpr)
      .ae_diagonal(tp, c(0L, 4L), c(1L, 5L))
    },
    ArrayElementwiseApplyFunc = {
      fd <- .ae_function_fdiff(expr$args[[1L]])
      sub <- expr$args[[2L]]
      tp <- .ae_tensor_product(.ad(sub, x), ArrayElementwiseApplyFunc(fd, sub))
      b <- get_ndim(x); rk <- get_ndim(expr)
      do.call(.ae_diagonal, c(list(tp), lapply(seq_len(rk) - 1L, function(i) c(b + i, b + rk + i))))
    },
    MatPow = .ad_matpow(expr, x),
    DiagMatrix = {
      vec <- expr$args[[1L]]
      darg <- .ad(vec, x)
      if (.ae_is_type(darg, "ZeroArray")) return(do.call(ZeroArray, .ae_shapes(x, expr)))
      n <- .mshape(expr)[[1L]]
      vaxis <- if (!is_one_num(.ae_shape(vec)[[1L]])) rank_x else rank_x + 1L
      td <- .ae_diagonal(.ae_tensor_product(darg, Identity(n)), c(vaxis, rank_x + 2L))
      Reshape(td, .ae_shapes(x, expr))
    },
    DiagonalMatrix = {
      sh <- .mshape(expr)
      if (!identical(sh[[1L]]$key, sh[[2L]]$key)) .perr("NotImplementedError", "derivative of a non-square DiagonalMatrix")
      darg <- .ad(expr$args[[1L]], x)
      if (.ae_is_type(darg, "ZeroArray")) return(do.call(ZeroArray, .ae_shapes(x, expr)))
      .ae_diagonal(.ae_tensor_product(darg, Identity(sh[[1L]])), c(rank_x, rank_x + 1L, rank_x + 2L))
    },
    HadamardProduct = .perr("NotImplementedError", ""),
    ArrayContraction = {
      fd <- .ad(expr$args[[1L]], x)
      do.call(.ae_contraction, c(list(fd), lapply(.ae_contr_ind(expr), function(g) g + rank_x)))
    },
    ArrayDiagonal = {
      fd <- .ad(expr$args[[1L]], x)
      do.call(.ae_diagonal, c(list(fd), lapply(.ae_diag_ind(expr), function(g) g + rank_x)))
    },
    ArrayAdd = do.call(.ae_array_add, lapply(expr$args, .ad, x)),
    PermuteDims = {
      de <- .ad(expr$args[[1L]], x)
      .ae_permute_dims(de, c(seq_len(rank_x) - 1L, .ae_paf(expr$args[[2L]]) + rank_x))
    },
    Reshape = Reshape(.ad(expr$args[[1L]], x), c(.ae_shape(x), expr$args[[2L]]$args)),
    Trace = .ad(.m2a(expr), x),
    Function = .ad_function(expr, x),
    Pow = .ad_pow(expr, x),
    Mul = {
      a <- expr$args
      do.call(ArrayAdd, lapply(seq_along(a), function(i)
        .ae_tensor_product(.mul(a[seq_len(i - 1L)]), .ad(a[[i]], x), .mul(a[-seq_len(i)]))))
    },
    Add = do.call(.ae_array_add, lapply(expr$args, .ad, x)),
    {
      if (.is_matexpr(expr)) {
        cg <- .m2a(expr)
        if (.ae_eq(cg, expr)) .perr("NotImplementedError", "")
        return(.ad(cg, x))
      }
      .ad_expr(expr, x)
    })
}

# isinstance(..., Expr) fallback
.ad_expr <- function(expr, x) {
  if (length(intersect(.ae_free_keys(expr), .ae_free_keys(x)))) {
    if (expr$type == "MatrixElement") {
      par <- expr$args[[1L]]
      if (.ae_eq(par, x) && x$type == "MatrixSymbol")
        return(MatrixUnit(.mshape(x)[[1L]], .mshape(x)[[2L]], expr$args[[2L]], expr$args[[3L]]))
      dpar <- .ad(par, x)
      rank_x <- length(.ae_shape(x))
      sh <- .mshape(par)
      tp <- .ae_tensor_product(dpar, MatrixUnit(sh[[1L]], sh[[2L]], expr$args[[2L]], expr$args[[3L]]))
      return(.ae_contraction(tp, c(rank_x, rank_x + 2L), c(rank_x + 1L, rank_x + 3L)))
    }
    .perr("NotImplementedError", "algorithm not implemented for this case")
  }
  do.call(ZeroArray, .ae_shape(x))
}

# arg.inv(): the _eval_inverse methods
.ae_minv <- function(arg) .me_inverse(arg)

.ad_tp <- function(expr, x) {
  args <- expr$args
  addends <- list()
  rx <- length(.ae_shape(x))
  for (i in seq_along(args)) {
    darg <- .ad(args[[i]], x)
    if (is_sym(darg) && is_zero_num(darg)) next
    prev <- args[seq_len(i - 1L)]; succ <- args[-seq_len(i)]
    sp <- sum(vapply(prev, get_ndim, 0L)); ss <- sum(vapply(succ, get_ndim, 0L))
    addend <- do.call(.ae_tensor_product, c(prev, list(darg), succ))
    t1 <- rx; t2 <- t1 + sp; t3 <- t2 + get_ndim(args[[i]]); t4 <- t3 + ss
    perm <- c(seq.int(t1, length.out = t2 - t1), seq_len(t1) - 1L, seq.int(t2, length.out = t3 - t2), seq.int(t3, length.out = t4 - t3))
    addends[[length(addends) + 1L]] <- .ae_permute_dims(addend, .ae_inv(perm))
  }
  if (length(addends) == 1L) return(addends[[1L]])
  if (!length(addends)) return(S_Zero)
  do.call(.ae_array_add, addends)
}

.ad_matpow <- function(expr, x) {
  base <- expr$args[[1L]]; ex <- expr$args[[2L]]
  dex <- .ad(ex, x)
  if (!.ae_is_type(dex, "ZeroArray") || is_zero_num(dex)) .perr("NotImplementedError", "")
  rank_x <- length(.ae_shape(x))
  if (is_zero_num(ex)) return(do.call(ZeroArray, .ae_shapes(x, expr)))
  if (ex$type == "Integer" && .int_val(ex) < 0) return(.ad(Inverse(doit(MatPow(base, .int(-.int_val(ex))))), x))
  dbase <- .ad(base, x)
  d <- Dummy("d")
  tp <- .ae_tensor_product(MatPow(base, d), dbase, MatPow(base, .sub(.sub(ex, d), S_One)))
  tc <- .ae_contraction(tp, c(1L, rank_x + 2L), c(rank_x + 3L, rank_x + 4L))
  pd <- .ae_permute_dims(tc, c(seq_len(rank_x), 0L, rank_x + 1L))
  if (.ae_is_type(pd, "ZeroArray")) return(pd)
  ArraySum(pd, Tuple(d, S_Zero, .sub(ex, S_One)))
}

.ae_fdiff <- function(e, i) {
  d <- .FUNS[[e$name]]
  if (!is.null(d) && !is.null(d$diff)) {
    r <- d$diff(e$args, i)
    if (!is.null(r)) return(r)
  }
  D <- Dummy("xi")
  a2 <- e$args; a2[[i]] <- D
  subs(.diff1(sym_rebuild(e, a2), D), D, e$args[[i]])
}

.ad_function <- function(expr, x) {
  addends <- list()
  for (i in seq_along(expr$args)) {
    arg <- expr$args[[i]]
    if (!length(intersect(.ae_free_keys(arg), .ae_free_keys(x)))) next
    if (length(.ae_shape(arg))) .perr("NotImplementedError", "derivative of a function of a non-scalar argument is not supported")
    darg <- .ad(arg, x)
    if (.ae_is_type(darg, "ZeroArray") || (is_sym(darg) && is_zero_num(darg))) next
    addends[[length(addends) + 1L]] <- .ae_tensor_product(.ae_fdiff(expr, i), darg)
  }
  if (!length(addends)) return(do.call(ZeroArray, .ae_shape(x)))
  do.call(.ae_array_add, addends)
}

# Pow._eval_derivative with an array variable
.ad_pow <- function(expr, x) {
  b <- expr$args[[1L]]; p <- expr$args[[2L]]
  db <- .ae_diff(b, x); dp <- .ae_diff(p, x)
  smul <- function(...) {
    a <- list(...)
    if (any(vapply(a, .ae_is_zero, NA))) return(S_Zero)
    if (any(vapply(a, .is_matexpr, NA))) return(do.call(.ae_mmul, a))
    if (any(vapply(a, .is_arrexpr, NA))) return(do.call(.ae_tensor_product, a))
    .mul(a)
  }
  sadd <- function(u, v) {
    if (.ae_is_zero(u)) return(v)
    if (.ae_is_zero(v)) return(u)
    if (.is_matexpr(u) || .is_matexpr(v)) return(.ae_madd(u, v))
    if (.is_arrexpr(u) || .is_arrexpr(v)) return(.ae_array_add(u, v))
    .add(list(u, v))
  }
  smul(expr, sadd(smul(dp, fcall("log", list(b))), smul(db, p, .pow(b, S_NegativeOne))))
}

.ad_matrixbase <- function(expr, x) {
  if (!length(intersect(.ae_free_keys(expr), .ae_free_keys(x)))) return(do.call(ZeroArray, .ae_shapes(x, expr)))
  if (.is_matexpr(x) && all(vapply(.mshape(x), function(s) s$type == "Integer", NA))) x <- as_explicit(x)
  if (.ae_is_explicit(x)) return(derive_by_array(.ae_as_arr(expr), .ae_as_arr(x)))
  .perr("NotImplementedError", "could not determine derivative")
}

.ad_ndarray <- function(expr, x) {
  if (!length(intersect(.ae_free_keys(expr), .ae_free_keys(x)))) return(do.call(ZeroArray, .ae_shapes(x, expr)))
  if ((.is_matexpr(x) || .ae_is_type(x, "ArraySymbol")) && all(vapply(.ae_shape(x), function(s) s$type == "Integer", NA)))
    x <- as_explicit(x)
  derive_by_array(.ae_as_arr(expr), if (.ae_is_explicit(x)) .ae_as_arr(x) else x)
}

# diff() arguments as list(list(var, count), ...), keeping matrices as given
.ae_diff_spec <- function(vars) {
  out <- list()
  i <- 1L
  while (i <= length(vars)) {
    v <- vars[[i]]
    if (is.list(v) && !is_sym(v) && !inherits(v, c("symmatrix", "SymArray"))) {
      out[[length(out) + 1L]] <- list(.ae_var(v[[1L]]), .dcount(as_sym(v[[2L]])))
      i <- i + 1L
      next
    }
    cnt <- 1L
    if (i < length(vars)) {
      nx <- vars[[i + 1L]]
      if (is.numeric(nx) || (is_sym(nx) && nx$type == "Integer")) {
        cnt <- .dcount(as_sym(nx))
        i <- i + 1L
      }
    }
    out[[length(out) + 1L]] <- list(.ae_var(v), cnt)
    i <- i + 1L
  }
  out
}
.ae_var <- function(v) if (is.numeric(v)) as_sym(v) else v

#' @rdname array_derive
#' @param ... Variables of differentiation, each optionally followed by a
#'   count.
#' @param evaluate If `TRUE`, compute the derivative; by default
#'   `ArrayDerivative()` stays unevaluated, as in SymPy.
#' @export
ArrayDerivative <- function(expr, ..., evaluate = FALSE) {
  spec <- .ae_diff_spec(list(...))
  # Derivative.__new__: every variable must have _diff_wrt True
  for (sp in spec) {
    v <- sp[[1L]]
    if (is_sym(v) && v$type %in% c("ArraySymbol", "Add", "Mul", "Pow"))
      .perr("ValueError", paste0("Can't calculate derivative wrt ", sstr(v), "."))
  }
  if (isTRUE(evaluate)) return(.ae_derivative_eval(expr, spec))
  spec <- lapply(spec, function(s) list(.ae_sympify(s[[1L]]), s[[2L]]))
  .derivative_node(.ae_sympify(expr), spec, sort = FALSE)
}

# Trace._eval_derivative
register_diff_hook("Trace", function(e, x) {
  if (.ae_is_type(x, "MatrixElement")) return(.me_diff_by_elem(.me_trace_sum(e), x))
  ex <- doit(e)
  if (is_sym(ex) && ex$type == "Trace") return(doit(Trace(.ae_diff(ex$args[[1L]], x))))
  .diff1(ex, x)
})

# Trace._eval_rewrite_as_Sum
.me_trace_sum <- function(e) {
  A <- e$args[[1L]]
  taken <- vapply(.free_sym_list(A), function(z) z$name, "")
  nm <- "i"
  k <- 0L
  while (nm %in% taken) { nm <- paste0("i", k); k <- k + 1L }
  i <- Symbol(nm)
  doit(Sum(.me_getitem(A, i, i), list(i, 0, .sub(.mshape(A)[[1L]], S_One))))
}

# MatrixElement._eval_derivative for a scalar variable
register_diff_hook("MatrixElement", function(e, x) {
  par <- e$args[[1L]]
  if (x$type != "MatrixElement") {
    pd <- .ae_diff(par, x)
    if (inherits(pd, "symmatrix")) pd <- .as_mexpr(pd)
    if (!.is_matexpr(pd) || !.ae_shape_eq(.mshape(pd), .mshape(par))) return(.derivative_node(e, list(list(x, 1L))))
    return(MatrixElement(pd, e$args[[2L]], e$args[[3L]]))
  }
  .me_elem_by_elem(e, x)
})

# d M[i, j] / d v[k, l] for matrix elements
.me_elem_by_elem <- function(e, v) {
  M <- e$args[[1L]]
  sh <- .mshape(M)
  if (identical(M$key, v$args[[1L]]$key))
    return(.mul(list(KroneckerDelta(e$args[[2L]], v$args[[2L]], Tuple(S_Zero, .sub(sh[[1L]], S_One))),
                     KroneckerDelta(e$args[[3L]], v$args[[3L]], Tuple(S_Zero, .sub(sh[[2L]], S_One))))))
  if (M$type == "Inverse") {
    i <- e$args[[2L]]; j <- e$args[[3L]]
    i1 <- Dummy("z1"); i2 <- Dummy("z2")
    Y <- M$args[[1L]]
    ys <- .mshape(Y)
    body <- .mul(list(MatrixElement(M, i, i1), .ae_diff(MatrixElement(Y, i1, i2), v), MatrixElement(M, i2, j)))
    return(.neg(Sum(body, list(i1, 0, .sub(ys[[1L]], S_One)), list(i2, 0, .sub(ys[[2L]], S_One)))))
  }
  if (.has1(e, v$args[[1L]])) return(.derivative_node(e, list(list(v, 1L))))
  S_Zero
}

# The derivative of a scalar expression with respect to a matrix element:
# the Add, Mul, Pow, function and Sum rules over MatrixElement leaves.
.me_diff_by_elem <- function(e, v) {
  if (identical(e$key, v$key)) return(S_One)
  if (!length(e$args) || !.has1(e, v$args[[1L]])) return(S_Zero)
  switch(e$type,
    MatrixElement = .me_elem_by_elem(e, v),
    Add = .add(lapply(e$args, .me_diff_by_elem, v)),
    Mul = {
      fs <- e$args
      terms <- list()
      for (i in seq_along(fs)) {
        d <- .me_diff_by_elem(fs[[i]], v)
        if (is_zero_num(d)) next
        fi <- fs; fi[[i]] <- d
        terms[[length(terms) + 1L]] <- .mul(fi)
      }
      .add(terms)
    },
    Pow = {
      b <- e$args[[1L]]; p <- e$args[[2L]]
      db <- .me_diff_by_elem(b, v); dp <- .me_diff_by_elem(p, v)
      if (is_zero_num(dp)) .mul(list(e, .div(.mul(list(db, p)), b)))
      else .mul(list(e, .add(list(.mul(list(dp, fcall("log", list(b)))), .div(.mul(list(p, db)), b)))))
    },
    Function = .add(lapply(seq_along(e$args), function(i) {
      d <- .me_diff_by_elem(e$args[[i]], v)
      if (is_zero_num(d)) S_Zero else .mul(list(.ae_fdiff(e, i), d))
    })),
    Sum = {
      lims <- e$args[-1L]
      last <- lims[[length(lims)]]
      f <- if (length(lims) > 1L) do.call(Sum, c(list(e$args[[1L]]), lapply(lims[-length(lims)], function(l) l$args))) else e$args[[1L]]
      if (.has1(last$args[[2L]], v) || .has1(last$args[[3L]], v)) return(.derivative_node(e, list(list(v, 1L))))
      Sum(.me_diff_by_elem(f, v), last$args)
    },
    {
      # other objects differentiate themselves (_eval_derivative)
      h <- .diff_hooks[[e$type]]
      if (!is.null(h)) h(e, v) else .derivative_node(e, list(list(v, 1L)))
    })
}

# ArrayElement._eval_derivative
register_diff_hook("ArrayElement", function(e, x) {
  if (x$type != "ArrayElement") return(S_Zero)
  if (identical(e$key, x$key)) return(S_One)
  if (!identical(e$args[[1L]]$key, x$args[[1L]]$key)) return(S_Zero)
  .mul(Map(function(i, j) KroneckerDelta(i, j), e$args[[2L]]$args, x$args[[2L]]$args))
})
