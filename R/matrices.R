# Symbolic matrices, after sympy/matrices.

.new_mat <- function(data, nr, nc) {
  structure(list(data = data, nrow = as.integer(nr), ncol = as.integer(nc)),
            class = "symmatrix")
}

.mat_map <- function(M, f) .keep_mclass(.new_mat(lapply(M$data, f), M$nrow, M$ncol), M)

# SymPy's matrix classes. A symmatrix without the "pyclass" attribute is a
# MutableDenseMatrix; methods keep the class of self, and binary operations
# take the class of higher _class_priority.
.MCLASS_PRIORITY <- c(MutableSparseMatrix = 3, MutableDenseMatrix = 4,
                      ImmutableDenseMatrix = 8, ImmutableSparseMatrix = 9)
.mclass <- function(M) {
  k <- attr(M, "pyclass")
  if (is.null(k) || !k %in% names(.MCLASS_PRIORITY)) "MutableDenseMatrix" else k
}
.set_mclass <- function(M, k) {
  attr(M, "pyclass") <- if (identical(k, "MutableDenseMatrix")) NULL else k
  M
}
# R gets the class of M (a matrix argument of the call that made R)
.keep_mclass <- function(R, M) {
  if (inherits(R, "symmatrix") && inherits(M, "symmatrix") && .mclass(M) != "MutableDenseMatrix" && is.null(attr(R, "pyclass")))
    attr(R, "pyclass") <- .mclass(M)
  R
}
.join_mclass <- function(R, A, B) {
  if (!inherits(R, "symmatrix") || !is.null(attr(R, "pyclass"))) return(R)
  ks <- c(if (inherits(A, "symmatrix")) .mclass(A), if (inherits(B, "symmatrix")) .mclass(B))
  .set_mclass(R, ks[which.max(.MCLASS_PRIORITY[ks])])
}
.is_immutable_mat <- function(M) .mclass(M) %in% c("ImmutableDenseMatrix", "ImmutableSparseMatrix")
.is_sparse_mat <- function(M) .mclass(M) %in% c("MutableSparseMatrix", "ImmutableSparseMatrix")

# internal coercion to a symmatrix, keeping a matrix (and its class) as is
.as_mat <- function(x) if (inherits(x, "symmatrix")) x else Matrix(x)

.mget <- function(M, i, j) M$data[[(i - 1L) * M$ncol + j]]

.mrows <- function(M) lapply(seq_len(M$nrow), function(i) M$data[(i - 1L) * M$ncol + seq_len(M$ncol)])

.from_rows <- function(rows) {
  nr <- length(rows)
  nc <- if (nr) length(rows[[1L]]) else 0L
  .new_mat(unlist(rows, recursive = FALSE), nr, nc)
}

#' Symbolic matrices
#'
#' `Matrix()` builds a matrix of expressions from a list of rows, from an R
#' matrix, from a vector (a column), or from dimensions and a function
#' `f(i, j)` or a data list in row-major order. `zeros()`, `ones()`, `eye()`
#' and `diag_matrix()` build special matrices. `+`, `-` and scalar `*` work as
#' usual; `*` between two matrices of equal shape is elementwise, `%*%` is the
#' matrix product and `%^%` the matrix power. `[` indexes like R matrices.
#' A scalar added to or subtracted from a matrix acts on every entry, as R
#' arithmetic does; SymPy raises a `TypeError` for `Matrix + 2` instead.
#'
#' @param data A list of rows, an R matrix, a vector, or a number of rows.
#' @param ncol Number of columns, when `data` is a number of rows.
#' @param f A function of `(i, j)` or a list of entries in row-major order.
#' @param n,m Dimensions.
#' @param ... Diagonal entries for `diag_matrix()`.
#' @return A `symmatrix`.
#' @examples
#' x <- Symbol("x")
#' M <- Matrix(list(list(1, x), list(x, 1)))
#' sym_det(M)
#' M %*% M
#' Matrix(2, 2, function(i, j) i + j)
#' @export
Matrix <- function(data, ncol = NULL, f = NULL) {
  # Matrix(M) is a new MutableDenseMatrix, whatever the class of M
  if (inherits(data, "symmatrix")) return(if (.mclass(data) == "MutableDenseMatrix") data else .set_mclass(data, "MutableDenseMatrix"))
  if (is.numeric(data) && length(data) == 1L && !is.null(ncol)) {
    nr <- as.integer(data)
    nc <- as.integer(ncol)
    if (is.function(f)) {
      d <- vector("list", nr * nc)
      for (i in seq_len(nr)) for (j in seq_len(nc)) d[[(i - 1L) * nc + j]] <- as_sym(f(i, j))
      return(.new_mat(d, nr, nc))
    }
    if (is.null(f)) return(zeros(nr, nc))
    if (length(f) != nr * nc) stop("ValueError: List length should be equal to rows*columns", call. = FALSE)
    return(.new_mat(lapply(as.list(f), as_sym), nr, nc))
  }
  if (is.matrix(data)) {
    nr <- nrow(data)
    nc <- ncol(data)
    d <- vector("list", nr * nc)
    for (i in seq_len(nr)) for (j in seq_len(nc)) d[[(i - 1L) * nc + j]] <- as_sym(data[[i, j]])
    return(.new_mat(d, nr, nc))
  }
  if (is.list(data) && !inherits(data, "sym")) {
    if (!length(data)) return(.new_mat(list(), 0L, 0L))
    if (all(vapply(data, function(r) is.list(r) && !inherits(r, "sym") || (is.atomic(r) && length(r) > 1L), TRUE))) {
      rows <- lapply(data, function(r) lapply(as.list(r), as_sym))
      lens <- vapply(rows, length, 0L)
      if (length(unique(lens)) != 1L) stop("rows must have equal length")
      return(.from_rows(rows))
    }
    return(.new_mat(lapply(data, as_sym), length(data), 1L))
  }
  if (is.atomic(data)) return(.new_mat(lapply(as.list(data), as_sym), length(data), 1L))
  if (.is_matexpr(data)) return(as_explicit(data))
  if (inherits(data, "sym")) return(.new_mat(list(data), 1L, 1L))
  stop("cannot build a Matrix from ", class(data)[1])
}

#' @rdname Matrix
#' @export
zeros <- function(n, m = n) .new_mat(rep(list(S_Zero), n * m), n, m)
#' @rdname Matrix
#' @export
ones <- function(n, m = n) .new_mat(rep(list(S_One), n * m), n, m)
#' @rdname Matrix
#' @export
eye <- function(n, m = n) {
  M <- zeros(n, m)
  for (i in seq_len(min(n, m))) M$data[[(i - 1L) * m + i]] <- S_One
  M
}
#' @rdname Matrix
#' @export
diag_matrix <- function(...) {
  d <- list(...)
  if (length(d) == 1L && is.list(d[[1L]]) && !inherits(d[[1L]], c("sym", "symmatrix"))) d <- d[[1L]]
  # Matrix.diag: a matrix argument is a block on the diagonal
  d <- lapply(d, function(v) if (inherits(v, "symmatrix")) v else .new_mat(list(as_sym(v)), 1L, 1L))
  nr <- sum(vapply(d, function(b) b$nrow, 0L))
  nc <- sum(vapply(d, function(b) b$ncol, 0L))
  M <- zeros(nr, nc)
  r0 <- 0L
  c0 <- 0L
  for (b in d) {
    for (i in seq_len(b$nrow)) for (j in seq_len(b$ncol))
      M$data[[(r0 + i - 1L) * nc + c0 + j]] <- .mget(b, i, j)
    r0 <- r0 + b$nrow
    c0 <- c0 + b$ncol
  }
  M
}

#' @export
dim.symmatrix <- function(x) c(x$nrow, x$ncol)

#' @export
length.symmatrix <- function(x) x$nrow * x$ncol

#' @export
`[.symmatrix` <- function(x, i, j, drop = TRUE) {
  nr <- x$nrow
  nc <- x$ncol
  if (missing(j) && !missing(i) && nargs() == 2L) {
    idx <- if (is.logical(i)) which(i) else i
    vals <- lapply(idx, function(k) {
      r <- (k - 1L) %% nr + 1L
      cc <- (k - 1L) %/% nr + 1L
      .mget(x, r, cc)
    })
    if (length(vals) == 1L) return(vals[[1L]])
    return(vals)
  }
  ri <- if (missing(i)) seq_len(nr) else if (is.logical(i)) which(i) else i
  ci <- if (missing(j)) seq_len(nc) else if (is.logical(j)) which(j) else j
  if (any(ri < 0)) ri <- setdiff(seq_len(nr), -ri)
  if (any(ci < 0)) ci <- setdiff(seq_len(nc), -ci)
  if (length(ri) == 1L && length(ci) == 1L && drop) return(.mget(x, ri, ci))
  d <- vector("list", length(ri) * length(ci))
  k <- 0L
  for (a in ri) for (b in ci) {
    k <- k + 1L
    d[[k]] <- .mget(x, a, b)
  }
  .keep_mclass(.new_mat(d, length(ri), length(ci)), x)
}

#' @export
`[<-.symmatrix` <- function(x, i, j, value) {
  if (.is_immutable_mat(x))
    stop("TypeError: Cannot set values of <class 'sympy.matrices.immutable.", .mclass(x), "'>", call. = FALSE)
  nr <- x$nrow
  nc <- x$ncol
  ri <- if (missing(i)) seq_len(nr) else i
  ci <- if (missing(j)) seq_len(nc) else j
  vals <- if (inherits(value, "symmatrix")) value$data else if (is.list(value) && !inherits(value, "sym")) value else list(value)
  vals <- lapply(vals, as_sym)
  data <- unclass(x)$data
  k <- 0L
  for (a in ri) for (b in ci) {
    k <- k + 1L
    data[[(a - 1L) * nc + b]] <- vals[[(k - 1L) %% length(vals) + 1L]]
  }
  .keep_mclass(.new_mat(data, nr, nc), x)
}

#' @export
`$<-.symmatrix` <- function(x, name, value) {
  y <- unclass(x)
  y[[name]] <- value
  structure(y, class = "symmatrix")
}

#' @export
format.symmatrix <- function(x, ...) {
  if (!x$nrow || !x$ncol) return("Matrix(0, 0, [])")
  strs <- vapply(x$data, .str, "")
  m <- matrix(strs, x$nrow, x$ncol, byrow = TRUE)
  w <- apply(m, 2L, function(col) max(nchar(col)))
  m[] <- vapply(seq_along(m), function(k) formatC(m[k], width = w[(k - 1L) %/% nrow(m) + 1L], flag = "-"), "")
  return(apply(m, 1L, function(r) paste0("[", paste(r, collapse = ", "), "]")))
  lines <- apply(m, 1L, function(r) paste0("[", paste(formatC(r, width = w, flag = "-"), collapse = ", "), "]"))
  lines
}

#' @export
print.symmatrix <- function(x, ...) {
  cat("Matrix(", x$nrow, ", ", x$ncol, ")\n", sep = "")
  cat(format(x), sep = "\n")
  invisible(x)
}

#' @export
as.matrix.symmatrix <- function(x, ...) {
  v <- lapply(x$data, function(e) tryCatch(.neval(e, list()), error = function(err) NA_real_))
  v <- unlist(v)
  if (anyNA(v)) warning("entries that are not numbers became NA; use N() or subs() first", call. = FALSE)
  matrix(v, x$nrow, x$ncol, byrow = TRUE)
}

#' @export
as.list.symmatrix <- function(x, ...) x$data

#' @export
t.symmatrix <- function(x) {
  d <- vector("list", length(x$data))
  for (i in seq_len(x$nrow)) for (j in seq_len(x$ncol)) d[[(j - 1L) * x$nrow + i]] <- .mget(x, i, j)
  .keep_mclass(.new_mat(d, x$ncol, x$nrow), x)
}

#' @export
chooseOpsMethod.symmatrix <- function(x, y, mx, my, cl, reverse) TRUE

#' @export
Ops.symmatrix <- function(e1, e2) {
  if (missing(e2)) {
    if (.Generic == "-") return(.mat_map(e1, .neg))
    return(e1)
  }
  m1 <- inherits(e1, "symmatrix")
  m2 <- inherits(e2, "symmatrix")
  # an explicit matrix with a matrix expression builds a matrix expression
  if (.is_matexpr(e1) || .is_matexpr(e2)) return(.matexpr_ops(.Generic, e1, e2))
  if (.Generic %in% c("==", "!=")) {
    eq <- m1 && m2 && e1$nrow == e2$nrow && e1$ncol == e2$ncol &&
      all(vapply(seq_along(e1$data), function(k) identical(e1$data[[k]]$key, e2$data[[k]]$key), TRUE))
    return(if (.Generic == "==") eq else !eq)
  }
  if (m1 && m2) {
    if (e1$nrow != e2$nrow || e1$ncol != e2$ncol) stop("matrix dimensions do not match")
    op <- switch(.Generic, "+" = function(a, b) .add(list(a, b)), "-" = .sub,
                 "*" = function(a, b) .mul(list(a, b)), "/" = .div,
                 stop("unsupported matrix operator ", .Generic))
    return(.join_mclass(.new_mat(Map(op, e1$data, e2$data), e1$nrow, e1$ncol), e1, e2))
  }
  M <- if (m1) e1 else e2
  s <- as_sym(if (m1) e2 else e1)
  op <- switch(.Generic,
               "+" = function(a) .add(list(a, s)),
               "-" = if (m1) function(a) .sub(a, s) else function(a) .sub(s, a),
               "*" = function(a) .mul(list(a, s)),
               "/" = if (m1) function(a) .div(a, s) else function(a) .div(s, a),
               "^" = if (m1) function(a) .pow(a, s) else function(a) .pow(s, a),
               stop("unsupported matrix operator ", .Generic))
  .mat_map(M, op)
}

.ops_hooks$dispatch <- function(gen, e1, e2) {
  # SetExpr has the higher _op_priority: its operators apply set arithmetic
  if (.is_setexpr(e1) || .is_setexpr(e2)) return(.setexpr_ops(gen, e1, e2))
  # Morphism.__mul__ composes
  if (gen == "*" && .is_morphism(e1) && .is_morphism(e2)) return(CompositeMorphism(e2, e1))
  if (inherits(e1, "symmatrix") || inherits(e2, "symmatrix")) {
    return(.ops_matrix_dispatch(gen, e1, e2))
  }
  if (.is_seq(e1) || .is_seq(e2)) return(.seq_ops(gen, e1, e2))
  if (.is_fourier(e1) || .is_fourier(e2)) return(.fourier_ops(gen, e1, e2))
  if (.is_matexpr(e1) || .is_matexpr(e2)) return(.matexpr_ops(gen, e1, e2))
  if (.fps_is(e1) || .fps_is(e2)) { r <- .ops_fps(gen, e1, e2); if (!is.null(r)) return(r) }
  NULL
}

.ops_matrix_dispatch <- function(gen, e1, e2) {
  f <- Ops.symmatrix
  environment(f) <- list2env(list(.Generic = gen), parent = environment(Ops.symmatrix))
  f(e1, e2)
}

.matmul <- function(A, B) {
  if (A$ncol != B$nrow) stop("non-conformable matrices")
  d <- vector("list", A$nrow * B$ncol)
  for (i in seq_len(A$nrow)) for (j in seq_len(B$ncol)) {
    d[[(i - 1L) * B$ncol + j]] <- .add(lapply(seq_len(A$ncol), function(k)
      .mul(list(.mget(A, i, k), .mget(B, k, j)))))
  }
  .new_mat(d, A$nrow, B$ncol)
}

#' @export
matrixOps.symmatrix <- function(x, y) {
  if (.Generic == "%*%" && (.is_matexpr(x) || .is_matexpr(y))) return(doit(MatMul(x, y)))
  if (.Generic == "%*%") {
    x <- .as_mat(x)
    y <- .as_mat(y)
    return(.join_mclass(.matmul(x, y), x, y))
  }
  stop("unsupported operator ", .Generic)
}

#' @rdname Matrix
#' @param A A matrix.
#' @export
`%^%` <- function(A, n) matrix_pow(A, n)

#' @export
cbind.symmatrix <- function(..., deparse.level = 1) {
  ms <- lapply(list(...), .as_mat)
  nr <- ms[[1L]]$nrow
  rows <- lapply(seq_len(nr), function(i) unlist(lapply(ms, function(M) .mrows(M)[[i]]), recursive = FALSE))
  .keep_mclass(.from_rows(rows), ms[[1L]])
}

#' @export
rbind.symmatrix <- function(..., deparse.level = 1) {
  ms <- lapply(list(...), .as_mat)
  .keep_mclass(.from_rows(unlist(lapply(ms, .mrows), recursive = FALSE)), ms[[1L]])
}

# ---------------------------------------------------------------------------
# Determinant, inverse, row reduction

# Bareiss fraction-free determinant of a list of rows.
# Bareiss elimination directly on the polynomials, for a matrix whose entries
# are polynomials in commutative symbols with rational coefficients. Every
# division is exact in that ring; NULL leaves other matrices to the
# expression path below.
.det_bareiss_qp <- function(rows) {
  n <- length(rows)
  ents <- unlist(rows, recursive = FALSE)
  if (any(vapply(ents, function(e) .has_type(e, "Float"), NA))) return(NULL)
  gens <- .find_gens(ents)
  if (!length(gens) || !all(vapply(gens, function(g) g$type %in% c("Symbol", "Dummy") && !identical(g$extra$facts[["commutative"]], FALSE), NA))) return(NULL)
  P <- vector("list", n)
  for (i in seq_len(n)) {
    P[[i]] <- vector("list", n)
    for (j in seq_len(n)) {
      q <- expr_to_qpoly(rows[[i]][[j]], gens, strict = TRUE)
      if (is.null(q) || any(q$e < 0L)) return(NULL)
      P[[i]][[j]] <- q
    }
  }
  sign <- 1
  prev <- NULL
  for (k in seq_len(n - 1L)) {
    if (qp_is_zero(P[[k]][[k]])) {
      sw <- NA
      for (i in (k + 1L):n) if (!qp_is_zero(P[[i]][[k]])) {
        sw <- i
        break
      }
      if (is.na(sw)) return(S_Zero)
      tmp <- P[[k]]; P[[k]] <- P[[sw]]; P[[sw]] <- tmp
      sign <- -sign
    }
    for (i in (k + 1L):n) {
      for (j in (k + 1L):n) {
        num <- qp_sub(qp_mul(P[[i]][[j]], P[[k]][[k]]), qp_mul(P[[i]][[k]], P[[k]][[j]]))
        if (!is.null(prev)) {
          num <- qp_exquo(num, prev)
          if (is.null(num)) return(NULL)
        }
        P[[i]][[j]] <- num
      }
    }
    prev <- P[[k]][[k]]
  }
  r <- P[[n]][[n]]
  if (sign < 0) r <- qp_neg(r)
  qpoly_to_expr(r, gens)
}

.det_bareiss <- function(rows) {
  n <- length(rows)
  if (n == 0L) return(S_One)
  if (n == 1L) return(rows[[1L]][[1L]])
  if (n > 2L && !is.null(r <- .det_bareiss_qp(rows))) return(r)
  if (n == 2L) return(expand(.sub(.mul(list(rows[[1L]][[1L]], rows[[2L]][[2L]])),
                                  .mul(list(rows[[1L]][[2L]], rows[[2L]][[1L]])))))
  M <- rows
  sign <- 1
  prev <- S_One
  for (k in seq_len(n - 1L)) {
    if (is_zero_num(cancel(M[[k]][[k]]))) {
      sw <- NA
      for (i in (k + 1L):n) if (!is_zero_num(cancel(M[[i]][[k]]))) {
        sw <- i
        break
      }
      if (is.na(sw)) return(S_Zero)
      tmp <- M[[k]]; M[[k]] <- M[[sw]]; M[[sw]] <- tmp
      sign <- -sign
    }
    for (i in (k + 1L):n) {
      for (j in (k + 1L):n) {
        num <- .sub(.mul(list(M[[i]][[j]], M[[k]][[k]])), .mul(list(M[[i]][[k]], M[[k]][[j]])))
        M[[i]][[j]] <- cancel(.div(num, prev))
      }
    }
    prev <- M[[k]][[k]]
  }
  r <- M[[n]][[n]]
  if (sign < 0) r <- .neg(r)
  cancel(r)
}

.det_berkowitz <- function(M) {
  cp <- .charpoly_coeffs(M)
  n <- M$nrow
  r <- cp[[n + 1L]]
  if (n %% 2L == 1L) r <- .neg(r)
  expand(r)
}

#' @noRd
det <- function(x, ...) {
  if (inherits(x, "symmatrix")) return(.det(x, ...))
  if (.is_matexpr(x)) return(.det_doit(x))
  base::det(x, ...)
}

.det <- function(M, method = "bareiss", iszerofunc = NULL) {
  if (is.character(iszerofunc) && !iszerofunc %in% c("_is_zero_after_expand_mul"))
    stop(.sympy_condition("ValueError", sprintf("Zero testing method '%s' unrecognized", iszerofunc)))
  if (M$nrow != M$ncol) stop(NonSquareMatrixError(""))
  method <- tolower(method)
  if (method == "bareis") method <- "bareiss" else if (method == "det_lu") method <- "lu"
  else if (method %in% c("domain-ge", "domain_ge")) method <- "bareiss"
  switch(method,
         bareiss = .det_bareiss(.mrows(M)),
         berkowitz = .det_berkowitz(M),
         laplace = .det_laplace(M),
         lu = .mg_det_lu(M),
         bird = .mg_det_bird(M),
         stop(.sympy_condition("ValueError", sprintf("Determinant method '%s' unrecognized", method))))
}

.det_laplace <- function(M) {
  n <- M$nrow
  if (n == 1L) return(.mget(M, 1L, 1L))
  terms <- lapply(seq_len(n), function(j) {
    s <- if (j %% 2L == 1L) S_One else S_NegativeOne
    .mul(list(s, .mget(M, 1L, j), .det_laplace(M[-1L, -j, drop = FALSE])))
  })
  expand(.add(terms))
}

#' @rdname sym_det
#' @param iszerofunc A function returning `TRUE`, `FALSE` or `NA` for an
#'   entry that is zero, nonzero or undecided; `NULL` uses the assumptions.
#' @param try_block_diag Invert each diagonal block separately.
#' @export
inverse <- function(M, method = NULL, iszerofunc = NULL, try_block_diag = FALSE) {
  # Function.inverse() of the characteristic curves
  if (is_sym(M) && M$type == "Function" && !is.null(.FUNS[[M$name]]$inverse)) return(.FUNS[[M$name]]$inverse())
  if (.is_matexpr(M) && !.is_explicit(M)) return(.me_inverse(M))
  .keep_mclass(.matint_inv(.as_mat(M), method, if (is.null(iszerofunc)) .matint_iszero else iszerofunc, try_block_diag), M)
}

#' @rdname sym_det
#' @param a A `symmatrix`.
#' @param b Optional right hand side.
#' @export
solve.symmatrix <- function(a, b, ...) {
  if (missing(b)) return(inverse(a))
  b <- .as_mat(b)
  n <- a$ncol
  aug <- cbind.symmatrix(a, b)
  R <- .rref(aug)
  if (any(R$pivots > n)) stop("linear system is inconsistent")
  if (length(R$pivots) < n) stop("linear system has infinitely many solutions; use linsolve()")
  R$M[seq_len(n), (n + 1L):(n + b$ncol), drop = FALSE]
}

# Matrix.rref() with its defaults: list(M = , pivots = ), 1 based pivots
.rref <- function(M) .matint_rref(M)

#' @rdname sym_det
#' @param simplify `TRUE`, `FALSE` or a function used to decide pivots.
#' @param pivots Also return the pivot columns.
#' @param normalize_last Scale the pivot rows to 1 at the end, as SymPy does.
#' @export
rref <- function(M, iszerofunc = NULL, simplify = FALSE, pivots = TRUE, normalize_last = TRUE) {
  r <- .matint_rref(.as_mat(M), if (is.null(iszerofunc)) .matint_iszero else iszerofunc, simplify, normalize_last)
  if (!pivots) return(r$M)
  list(matrix = r$M, pivots = r$pivots)
}

#' @rdname sym_det
#' @export
rank_matrix <- function(M) length(.rref(.as_mat(M))$pivots)

#' @rdname sym_det
#' @export
nullspace <- function(M, simplify = FALSE, iszerofunc = NULL) {
  M <- .as_mat(M)
  r <- .matint_rref(M, if (is.null(iszerofunc)) .matint_iszero else iszerofunc, simplify, TRUE)
  free <- setdiff(seq_len(M$ncol), r$pivots)
  lapply(free, function(fv) {
    v <- rep(list(S_Zero), M$ncol)
    v[[fv]] <- S_One
    for (k in seq_along(r$pivots)) v[[r$pivots[k]]] <- .neg(.mget(r$M, k, fv))
    .new_mat(v, M$ncol, 1L)
  })
}

#' @rdname sym_det
#' @export
columnspace <- function(M) {
  M <- .as_mat(M)
  r <- .rref(M)
  lapply(r$pivots, function(c) M[, c, drop = FALSE])
}

#' @rdname sym_det
#' @export
rowspace <- function(M) {
  r <- .rref(.as_mat(M))
  lapply(seq_along(r$pivots), function(k) r$M[k, , drop = FALSE])
}

#' @rdname sym_det
#' @export
minor <- function(M, i, j, method = "berkowitz") {
  if (inherits(M, "GeoEllipse")) return(.ellipse_minor(M))
  M <- .as_mat(M)
  if (M$nrow != M$ncol || M$nrow < 1L) stop(NonSquareMatrixError(""))
  if (i < 1L || i > M$nrow || j < 1L || j > M$ncol)
    stop(.sympy_condition("ValueError", sprintf("`i` and `j` must satisfy 1 <= i <= ``M.rows`` (%d)and 1 <= j <= ``M.cols`` (%d).", M$nrow, M$ncol)))
  .det(M[-i, -j, drop = FALSE], method)
}

#' @rdname sym_det
#' @export
cofactor <- function(M, i, j, method = "berkowitz") {
  m <- minor(M, i, j, method)
  if ((i + j) %% 2L == 1L) .neg(m) else m
}

#' @rdname sym_det
#' @export
cofactor_matrix <- function(M, method = "berkowitz") {
  M <- .as_mat(M)
  if (M$nrow != M$ncol) stop(NonSquareMatrixError(""))
  Matrix(M$nrow, M$ncol, function(i, j) cofactor(M, i, j, method))
}

#' @rdname sym_det
#' @export
adjugate <- function(M, method = "berkowitz") .keep_mclass(t(cofactor_matrix(M, method)), M)

#' @rdname sym_det
#' @export
trace_matrix <- function(M) {
  M <- .as_mat(M)
  .add(lapply(seq_len(min(M$nrow, M$ncol)), function(i) .mget(M, i, i)))
}

# ---------------------------------------------------------------------------
# Characteristic polynomial and eigenvalues

# Berkowitz algorithm: coefficients of det(lambda*I - M), highest first.
.charpoly_coeffs <- function(M) {
  n <- M$nrow
  if (n == 0L) return(list(S_One))
  A <- .mrows(M)
  vec <- list(S_One, .neg(A[[1L]][[1L]]))
  if (n == 1L) return(vec)
  for (k in 2:n) {
    R <- A[[k]][seq_len(k - 1L)]
    C <- lapply(seq_len(k - 1L), function(i) A[[i]][[k]])
    a <- A[[k]][[k]]
    Am <- lapply(seq_len(k - 1L), function(i) A[[i]][seq_len(k - 1L)])
    items <- list(S_One, .neg(a))
    cur <- C
    terms <- list(.neg(.add(lapply(seq_len(k - 1L), function(i) .mul(list(R[[i]], cur[[i]]))))))
    for (p in seq_len(k - 2L)) {
      cur <- lapply(seq_len(k - 1L), function(i) .add(lapply(seq_len(k - 1L), function(j)
        .mul(list(Am[[i]][[j]], cur[[j]])))))
      terms[[length(terms) + 1L]] <- .neg(.add(lapply(seq_len(k - 1L), function(i) .mul(list(R[[i]], cur[[i]])))))
    }
    T <- c(items, terms)
    # Toeplitz multiply: new_vec = T_matrix (k+1 x k) %*% vec
    newv <- vector("list", k + 1L)
    for (i in seq_len(k + 1L)) {
      acc <- list()
      for (j in seq_len(min(i, length(vec)))) {
        ti <- i - j + 1L
        if (ti <= length(T)) acc[[length(acc) + 1L]] <- .mul(list(T[[ti]], vec[[j]]))
      }
      newv[[i]] <- expand(.add(acc))
    }
    vec <- newv
  }
  vec
}

#' Characteristic polynomial, eigenvalues and eigenvectors
#'
#' `charpoly(M, lambda)` returns `sym_det(lambda*I - M)`. `eigenvals()` returns
#' the roots of the characteristic polynomial with multiplicities;
#' `eigenvects()` returns a list of `list(value, multiplicity, vectors)`.
#' `diagonalize()` returns `list(P, D)` with `M = P D P^-1`. `is_diagonalizable()`
#' tests whether that is possible.
#'
#' @param M A square `symmatrix`.
#' @param lambda Symbol for the characteristic polynomial.
#' @return An expression or list.
#' @examples
#' M <- Matrix(list(list(2, 1), list(1, 2)))
#' charpoly(M)
#' eigenvals(M)
#' @export
charpoly <- function(M, lambda = Symbol("lambda")) {
  M <- .as_mat(M)
  cf <- .charpoly_coeffs(M)
  lambda <- as_sym(lambda)
  n <- length(cf) - 1L
  # PurePoly.as_expr(): each coefficient expanded, collected by power
  .add(lapply(seq_along(cf), function(k) .mul(list(expand(cf[[k]]), .pow(lambda, .int(n - k + 1L))))))
}

#' @rdname charpoly
#' @param error_when_incomplete Signal a `MatrixError` when some eigenvalues
#'   cannot be found.
#' @param simplify `FALSE`, `TRUE` or a function applied to each result.
#' @param multiple Return the eigenvalues as a list repeated by multiplicity
#'   instead of a list with attribute `"multiplicity"`.
#' @param rational Replace floating point entries by rationals first.
#' @param ... Flags passed to [roots()].
#' @param iszerofunc A zero test returning `TRUE`, `FALSE` or `NA`; `NULL`
#'   uses the default test.
#' @param chop Accepted for SymPy compatibility; currently not used.
#' @export
eigenvals <- function(M, error_when_incomplete = TRUE, simplify = FALSE, multiple = FALSE, rational = FALSE, ...)
  .matint_eigenvals(M, error_when_incomplete, simplify, multiple, rational, list(...))

#' @rdname charpoly
#' @export
eigenvects <- function(M, error_when_incomplete = TRUE, iszerofunc = NULL, simplify = FALSE, chop = FALSE, ...)
  .matint_eigenvects(M, error_when_incomplete, iszerofunc, simplify, chop, list(...))

#' @rdname charpoly
#' @param reals_only Treat a matrix with complex eigenvalues as not
#'   diagonalizable.
#' @export
is_diagonalizable <- function(M, reals_only = FALSE) .matint_is_diagonalizable(M, reals_only)

#' @rdname charpoly
#' @param sort Sort the eigenvalues on the diagonal.
#' @param normalize Scale each column of `P` to unit norm.
#' @export
diagonalize <- function(M, reals_only = FALSE, sort = FALSE, normalize = FALSE)
  .matint_diagonalize(M, reals_only, sort, normalize)

# Apply a scalar function to a diagonalizable matrix.
.matrix_function <- function(M, f) {
  d <- diagonalize(M)
  n <- M$nrow
  fd <- diag_matrix(lapply(seq_len(n), function(i) f(.mget(d$D, i, i))))
  .mat_map(.matmul(.matmul(d$P, fd), inverse(d$P)), function(z) simplify(z))
}

#' Matrix exponential
#'
#' `matrix_exp(M)` computes `exp(M)` from the Jordan normal form, as SymPy's
#' `Matrix.exp()`.
#'
#' @param M A square `symmatrix`.
#' @return A `symmatrix`.
#' @export
matrix_exp <- function(M) {
  M <- .as_mat(M)
  if (M$nrow != M$ncol) stop(NonSquareMatrixError("Exponentiation is valid only for square matrices"))
  jf <- tryCatch(jordan_form(M), MatrixError = function(e)
    stop(.sympy_condition("NotImplementedError", "Exponentiation is implemented only for matrices for which the Jordan normal form can be computed")))
  J <- jf$J
  n <- J$nrow
  # the Jordan cells, and exp of each: exp(l)/i! on the i-th superdiagonal
  starts <- c(1L, which(vapply(seq_len(n - 1L), function(i) is_zero_num(.mget(J, i, i + 1L)), TRUE)) + 1L)
  sizes <- diff(c(starts, n + 1L))
  blocks <- lapply(seq_along(starts), function(k) {
    sz <- sizes[k]
    exp_l <- fcall("exp", list(.mget(J, starts[k], starts[k])))
    B <- zeros(sz, sz)
    for (i in seq_len(sz)) for (j in i:sz) B$data[[(i - 1L) * sz + j]] <- .div(exp_l, factorial(.int(j - i)))
    B
  })
  ret <- .matmul(.matmul(jf$P, do.call(diag_matrix, blocks)), inverse(jf$P))
  vals <- Filter(function(z) !is_zero_num(z), M$data)
  if (all(vapply(vals, function(z) isTRUE(ask_fact(z, "real")), TRUE))) ret <- .mat_map(ret, re)
  ret
}

# ---------------------------------------------------------------------------
# Decompositions

#' Matrix decompositions
#'
#' `LUdecomposition()` returns `list(L, U, perm)` with `P*M = L*U`, where
#' `perm` lists the row swaps (1 based) that make `P`;
#' `QRdecomposition()` returns `list(Q, R)` by Gram-Schmidt;
#' `cholesky()` returns `L` with `M = L*t(L)` for symmetric positive definite
#' matrices; `LDLdecomposition()` returns `list(L, D)`. `GramSchmidt()`
#' orthogonalizes a list of vectors.
#'
#' @param M A `symmatrix`.
#' @param vlist A list of column vectors.
#' @param orthonormal Normalize the vectors.
#' @param iszerofunc,simpfunc,rankcheck For `LUdecomposition()`: the zero
#'   test, an optional simplification for the pivot search, and whether to
#'   raise for a rank deficient matrix.
#' @param hermitian For `cholesky()` and `LDLdecomposition()`: require a Hermitian
#'   matrix and use conjugates (SymPy's default); `FALSE` requires a symmetric matrix.
#' @return A list of matrices.
#' @export
LUdecomposition <- function(M, iszerofunc = NULL, simpfunc = NULL, rankcheck = FALSE) {
  M <- .as_mat(M)
  r <- .mg_lu(M, iszerofunc, simpfunc, rankcheck)
  cm <- r$lu
  L <- .mg_new(cm$nrow, cm$nrow, function(i, j) if (i < j) S_Zero else if (i == j) S_One else if (j <= cm$ncol) .mget(cm, i, j) else S_Zero)
  U <- .mg_new(cm$nrow, cm$ncol, function(i, j) if (i > j) S_Zero else .mget(cm, i, j))
  list(L = L, U = U, perm = r$swaps)
}

#' @rdname LUdecomposition
#' @export
GramSchmidt <- function(vlist, orthonormal = FALSE) .mg_orthogonalize(vlist, normalize = orthonormal, rankcheck = TRUE)

.dot <- function(a, b) .add(Map(function(x, y) .mul(list(fcall("conjugate", list(x)), y)), a$data, b$data))

#' @rdname LUdecomposition
#' @export
QRdecomposition <- function(M) .mg_qr(.as_mat(M), normalize = TRUE)

#' @rdname LUdecomposition
#' @export
cholesky <- function(M, hermitian = TRUE) {
  M <- .ldl_check(.as_mat(M), hermitian)
  n <- M$nrow
  L <- .mrows(zeros(n))
  cj <- if (hermitian) function(z) fcall("conjugate", list(z)) else identity
  for (i in seq_len(n)) {
    for (j in seq_len(i - 1L)) {
      s <- .add(lapply(seq_len(j - 1L), function(k) .mul(list(L[[i]][[k]], cj(L[[j]][[k]])))))
      L[[i]][[j]] <- .mul(list(.div(S_One, L[[j]][[j]]), .sub(.mget(M, i, j), s)))
    }
    s <- .add(lapply(seq_len(i - 1L), function(k) if (hermitian) .mul(list(L[[i]][[k]], cj(L[[i]][[k]]))) else .pow(L[[i]][[k]], S_Two)))
    d <- .sub(.mget(M, i, i), s)
    if (hermitian && isFALSE(ask_fact(d, "positive"))) stop(NonPositiveDefiniteMatrixError("Matrix must be positive-definite"))
    L[[i]][[i]] <- .sqrt(d)
  }
  .from_rows(L)
}

# _cholesky and _LDLdecomposition: square, then Hermitian (or symmetric)
.ldl_check <- function(M, hermitian) {
  if (M$nrow != M$ncol) stop(NonSquareMatrixError("Matrix must be square."))
  if (hermitian && !.is_hermitian(M)) stop(.sympy_condition("ValueError", "Matrix must be Hermitian."))
  if (!hermitian && !is_symmetric(M)) stop(.sympy_condition("ValueError", "Matrix must be symmetric."))
  M
}
.is_hermitian <- function(M) all(vapply(seq_along(M$data), function(k) {
  i <- (k - 1L) %/% M$ncol + 1L
  j <- (k - 1L) %% M$ncol + 1L
  is_zero_num(simplify(.sub(.mget(M, i, j), fcall("conjugate", list(.mget(M, j, i))))))
}, TRUE))

#' @rdname LUdecomposition
#' @export
LDLdecomposition <- function(M, hermitian = TRUE) {
  M <- .ldl_check(.as_mat(M), hermitian)
  n <- M$nrow
  L <- .mrows(eye(n))
  D <- rep(list(S_Zero), n)
  cj <- if (hermitian) function(z) fcall("conjugate", list(z)) else identity
  for (i in seq_len(n)) {
    for (j in seq_len(i - 1L)) {
      s <- .add(lapply(seq_len(j - 1L), function(k) .mul(list(L[[i]][[k]], cj(L[[j]][[k]]), D[[k]]))))
      L[[i]][[j]] <- .mul(list(.div(S_One, D[[j]]), .sub(.mget(M, i, j), s)))
    }
    s <- .add(lapply(seq_len(i - 1L), function(k)
      if (hermitian) .mul(list(L[[i]][[k]], cj(L[[i]][[k]]), D[[k]])) else .mul(list(.pow(L[[i]][[k]], S_Two), D[[k]]))))
    D[[i]] <- .sub(.mget(M, i, i), s)
    if (hermitian && isFALSE(ask_fact(D[[i]], "positive"))) stop(NonPositiveDefiniteMatrixError("Matrix must be positive-definite"))
  }
  list(L = .from_rows(L), D = diag_matrix(D))
}

# ---------------------------------------------------------------------------
# Misc

#' Matrix utilities
#'
#' `dot()` and `cross()` for vectors, `norm_matrix()` (Frobenius or vector
#' 2-norm), `applyfunc()` elementwise, `is_symmetric()`, `is_square()`,
#' `is_zero_matrix()`, `hadamard_product()`, `kronecker_product()`, `vec()`
#' and `pinv()` (Moore-Penrose pseudo inverse for full rank matrices).
#' Given an integer rather than a matrix, `is_square()` tests for a perfect
#' square, as `sympy.ntheory.primetest.is_square` does.
#'
#' @param a,b Vectors as matrices, or [CoordSys3D()] vectors, dyadics and
#'   the [Del()] operator, for which the methods of `sympy.vector` apply.
#' @param doit For the del operator: evaluate the derivatives.
#' @param hermitian Conjugate one of the vectors.
#' @param conjugate_convention Which vector to conjugate: `"maths"` or
#'   `"left"` the first, `"physics"` or `"right"` the second.
#' @param M A matrix.
#' @param f A function applied to each entry.
#' @return An expression or matrix.
#' @export
dot <- function(a, b, doit = FALSE, hermitian = NULL, conjugate_convention = NULL) {
  if (inherits(a, "GeoPoint")) return(.mg_pdot(a, if (inherits(b, "GeoPoint")) b else Point(b)))
  if (.pv_dispatch(a, b)) return(.pv_dot_dispatch(a, b))
  if (inherits(a, "Del")) return(.vec_divergence(b, doit = doit))
  if (.is_bd(a) || .is_bd(b)) return(.bd_dot_method(a, b))
  a <- .as_mat(a)
  b <- .as_mat(b)
  # MatrixBase.dot: two vectors of equal length
  if (!(1L %in% c(a$nrow, a$ncol)) || !(1L %in% c(b$nrow, b$ncol)) || length(a$data) != length(b$data))
    stop(sprintf("ShapeError: Dimensions incorrect for dot product: (%d, %d), (%d, %d)", a$nrow, a$ncol, b$nrow, b$ncol), call. = FALSE)
  if (!is.null(conjugate_convention) && is.null(hermitian)) hermitian <- TRUE
  if (isTRUE(hermitian) && is.null(conjugate_convention)) conjugate_convention <- "maths"
  if (isTRUE(hermitian)) {
    if (conjugate_convention %in% c("maths", "left", "math")) a <- conjugate(a)
    else if (conjugate_convention %in% c("physics", "right")) b <- conjugate(b)
    else stop(.sympy_condition("ValueError", "Unknown conjugate_convention was entered. conjugate_convention must be one of the following: math, maths, left, physics or right."))
  }
  .add(Map(function(x, y) .mul(list(x, y)), a$data, b$data))
}

#' @rdname dot
#' @export
cross <- function(a, b, doit = FALSE) {
  if (.pv_dispatch(a, b)) return(.pv_cross_dispatch(a, b))
  if (inherits(a, "Del")) return(.vec_curl(b, doit = doit))
  if (.is_bd(a) || .is_bd(b)) return(.vec_cross_method(as_sym(a), b))
  a <- .as_mat(a)$data
  b <- .as_mat(b)$data
  if (length(a) != 3L || length(b) != 3L) stop("cross product needs 3-vectors")
  Matrix(list(.sub(.mul(list(a[[2]], b[[3]])), .mul(list(a[[3]], b[[2]]))),
              .sub(.mul(list(a[[3]], b[[1]])), .mul(list(a[[1]], b[[3]]))),
              .sub(.mul(list(a[[1]], b[[2]])), .mul(list(a[[2]], b[[1]])))))
}

#' @rdname dot
#' @export
norm_matrix <- function(M) {
  M <- .as_mat(M)
  .sqrt(.add(lapply(M$data, function(z) .sq(fcall("Abs", list(z))))))
}

#' @rdname dot
#' @export
applyfunc <- function(M, f) {
  if (.is_matexpr(M) && !.is_explicit(M)) return(ElementwiseApplyFunction(f, M))
  M <- .as_mat(M)
  g <- function(z) as_sym(if (is.function(f)) f(z) else apply_lambda(f, z))
  # SparseRepMatrix.applyfunc maps only the stored (nonzero) entries
  if (.is_sparse_mat(M)) return(.mat_map(M, function(z) if (is_zero_num(z)) z else g(z)))
  .mat_map(M, g)
}

#' @rdname dot
#' @export
is_square <- function(M) {
  int_like <- (is_sym(M) && M$type == "Integer") || inherits(M, "bigz") ||
    (is.numeric(M) && length(M) == 1L && is.null(dim(M))) ||
    (is.character(M) && length(M) == 1L && grepl("^\\s*[-+]?[0-9]+\\s*$", M))
  if (int_like) return(.is_square_int(M))
  .as_mat(M)$nrow == .as_mat(M)$ncol
}

#' @rdname dot
#' @export
is_symmetric <- function(M) {
  M <- .as_mat(M)
  is_square(M) && all(vapply(seq_along(M$data), function(k) {
    i <- (k - 1L) %/% M$ncol + 1L
    j <- (k - 1L) %% M$ncol + 1L
    is_zero_num(simplify(.sub(.mget(M, i, j), .mget(M, j, i))))
  }, TRUE))
}

#' @rdname dot
#' @export
is_zero_matrix <- function(M) {
  # MatrixBase._eval_is_zero_matrix: FALSE if an entry is known nonzero,
  # NA (SymPy's None) if an entry is undecided
  z <- vapply(.as_mat(M)$data, function(e) as.logical(ask_fact(e, "zero")), NA)
  if (any(z %in% FALSE)) return(FALSE)
  if (anyNA(z)) return(NA)
  TRUE
}

#' @rdname dot
#' @export
vec <- function(M) {
  M <- .as_mat(M)
  .new_mat(t(M)$data, M$nrow * M$ncol, 1L)
}

#' @rdname dot
#' @export
pinv <- function(M) {
  # _pinv_rank_decomposition of sympy/matrices/inverse.py
  M <- .as_mat(M)
  H <- .mat_h(M)
  if (all(vapply(M$data, is_zero_num, TRUE))) return(H)
  if (M$nrow == 1L || M$ncol == 1L)
    return(.mat_map(H, function(z) .div(z, .add(lapply(M$data, function(a) .mul(list(fcall("conjugate", list(a)), a)))))))
  rd <- .rref(M)
  rank <- length(rd$pivots)
  B <- M[, rd$pivots, drop = FALSE]
  C <- rd$M[seq_len(rank), , drop = FALSE]
  .matmul(.pinv_full_rank(C), .pinv_full_rank(B))
}

.mat_h <- function(M) t(.mat_map(M, function(z) fcall("conjugate", list(z))))

.pinv_full_rank <- function(M) {
  H <- .mat_h(M)
  if (all(vapply(M$data, is_zero_num, TRUE))) return(H)
  if (M$nrow >= M$ncol) .matmul(inverse(.matmul(H, M)), H) else .matmul(H, inverse(.matmul(M, H)))
}

#' @export
diff.symmatrix <- function(x, ...) {
  if (.ae_derivative_applies(NULL, lapply(list(...), list))) return(.ae_derivative_eval(x, .ae_diff_spec(list(...))))
  .mat_map(x, function(z) diff(z, ...))
}
