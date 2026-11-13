# Symbolic matrices, after sympy/matrices.

.new_mat <- function(data, nr, nc) {
  structure(list(data = data, nrow = as.integer(nr), ncol = as.integer(nc)),
            class = "symmatrix")
}

.mat_map <- function(M, f) .new_mat(lapply(M$data, f), M$nrow, M$ncol)

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
  if (inherits(data, "symmatrix")) return(data)
  if (is.numeric(data) && length(data) == 1L && !is.null(ncol)) {
    nr <- as.integer(data)
    nc <- as.integer(ncol)
    if (is.function(f)) {
      d <- vector("list", nr * nc)
      for (i in seq_len(nr)) for (j in seq_len(nc)) d[[(i - 1L) * nc + j]] <- as_sym(f(i, j))
      return(.new_mat(d, nr, nc))
    }
    if (is.null(f)) return(zeros(nr, nc))
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
  .new_mat(d, length(ri), length(ci))
}

#' @export
`[<-.symmatrix` <- function(x, i, j, value) {
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
  .new_mat(data, nr, nc)
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
  v <- vapply(x$data, function(e) {
    r <- tryCatch(.neval(e, list()), error = function(err) NA_real_)
    if (is.complex(r)) NA_real_ else r
  }, 0)
  matrix(v, x$nrow, x$ncol, byrow = TRUE)
}

#' @export
as.list.symmatrix <- function(x, ...) x$data

#' @export
t.symmatrix <- function(x) {
  d <- vector("list", length(x$data))
  for (i in seq_len(x$nrow)) for (j in seq_len(x$ncol)) d[[(j - 1L) * x$nrow + i]] <- .mget(x, i, j)
  .new_mat(d, x$ncol, x$nrow)
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
    return(.new_mat(Map(op, e1$data, e2$data), e1$nrow, e1$ncol))
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
    x <- Matrix(x)
    y <- Matrix(y)
    return(.matmul(x, y))
  }
  stop("unsupported operator ", .Generic)
}

#' @rdname Matrix
#' @param A A matrix.
#' @export
`%^%` <- function(A, n) {
  A <- Matrix(A)
  if (A$nrow != A$ncol) stop(NonSquareMatrixError("matrix power needs a square matrix"))
  n <- as_sym(n)
  if (n$type != "Integer") return(.matrix_function(A, function(l, t) .pow(l, n)))
  k <- gmp::numerator(n$value)
  if (k < 0) {
    A <- inverse(A)
    k <- -k
  }
  res <- eye(A$nrow)
  base <- A
  while (k > 0) {
    if (k %% 2 == 1) res <- .mat_map(.matmul(res, base), expand)
    k <- k %/% 2
    if (k > 0) base <- .mat_map(.matmul(base, base), expand)
  }
  res
}

#' @export
cbind.symmatrix <- function(..., deparse.level = 1) {
  ms <- lapply(list(...), Matrix)
  nr <- ms[[1L]]$nrow
  rows <- lapply(seq_len(nr), function(i) unlist(lapply(ms, function(M) .mrows(M)[[i]]), recursive = FALSE))
  .from_rows(rows)
}

#' @export
rbind.symmatrix <- function(..., deparse.level = 1) {
  ms <- lapply(list(...), Matrix)
  .from_rows(unlist(lapply(ms, .mrows), recursive = FALSE))
}

# ---------------------------------------------------------------------------
# Determinant, inverse, row reduction

# Bareiss fraction-free determinant of a list of rows.
.det_bareiss <- function(rows) {
  n <- length(rows)
  if (n == 0L) return(S_One)
  if (n == 1L) return(rows[[1L]][[1L]])
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

.det <- function(M, method = "bareiss") {
  if (M$nrow != M$ncol) stop(NonSquareMatrixError("determinant needs a square matrix"))
  switch(method,
         bareiss = .det_bareiss(.mrows(M)),
         berkowitz = .det_berkowitz(M),
         laplace = .det_laplace(M),
         stop("unknown determinant method"))
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
#' @export
inverse <- function(M, method = "GE") {
  # Function.inverse() of the characteristic curves
  if (is_sym(M) && M$type == "Function" && !is.null(.FUNS[[M$name]]$inverse)) return(.FUNS[[M$name]]$inverse())
  if (.is_matexpr(M) && !.is_explicit(M)) return(.me_inverse(M))
  M <- Matrix(M)
  n <- M$nrow
  if (n != M$ncol) stop(NonSquareMatrixError("inverse needs a square matrix"))
  if (method == "ADJ") {
    d <- det(M)
    if (is_zero_num(d)) stop(NonInvertibleMatrixError("matrix is not invertible"))
    return(.mat_map(adjugate(M), function(z) cancel(.div(z, d))))
  }
  aug <- cbind.symmatrix(M, eye(n))
  R <- .rref(aug)
  if (length(R$pivots) < n || any(R$pivots[seq_len(n)] != seq_len(n))) stop(NonInvertibleMatrixError("matrix is not invertible"))
  R$M[, (n + 1L):(2L * n), drop = FALSE]
}

#' @rdname sym_det
#' @param a A `symmatrix`.
#' @param b Optional right hand side.
#' @export
solve.symmatrix <- function(a, b, ...) {
  if (missing(b)) return(inverse(a))
  b <- Matrix(b)
  n <- a$ncol
  aug <- cbind.symmatrix(a, b)
  R <- .rref(aug)
  if (any(R$pivots > n)) stop("linear system is inconsistent")
  if (length(R$pivots) < n) stop("linear system has infinitely many solutions; use linsolve()")
  R$M[seq_len(n), (n + 1L):(n + b$ncol), drop = FALSE]
}

.rref <- function(M, simplify_fn = cancel) {
  A <- .mrows(M)
  nr <- M$nrow
  nc <- M$ncol
  pivots <- integer()
  r <- 1L
  for (c in seq_len(nc)) {
    if (r > nr) break
    for (i in r:nr) A[[i]][[c]] <- simplify_fn(A[[i]][[c]])
    fp <- .coremm_find_pivot(lapply(r:nr, function(i) A[[i]][[c]]))
    for (nd in fp[[3L]]) A[[r + nd[[1L]] - 1L]][[c]] <- nd[[2L]]
    if (is.na(fp[[1L]])) next
    pr <- r + fp[[1L]] - 1L
    if (pr != r) {
      tmp <- A[[r]]; A[[r]] <- A[[pr]]; A[[pr]] <- tmp
    }
    pv <- A[[r]][[c]]
    A[[r]] <- lapply(A[[r]], function(z) simplify_fn(.div(z, pv)))
    for (i in seq_len(nr)) {
      if (i == r) next
      f <- A[[i]][[c]]
      if (is_zero_num(f)) next
      A[[i]] <- lapply(seq_len(nc), function(j) simplify_fn(.sub(A[[i]][[j]], .mul(list(f, A[[r]][[j]])))))
    }
    pivots <- c(pivots, c)
    r <- r + 1L
  }
  list(M = .from_rows(A), pivots = pivots)
}

#' @rdname sym_det
#' @export
rref <- function(M) {
  r <- .rref(Matrix(M))
  list(matrix = r$M, pivots = r$pivots)
}

#' @rdname sym_det
#' @export
rank_matrix <- function(M) length(.rref(Matrix(M))$pivots)

#' @rdname sym_det
#' @export
nullspace <- function(M) {
  M <- Matrix(M)
  r <- .rref(M)
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
  M <- Matrix(M)
  r <- .rref(M)
  lapply(r$pivots, function(c) M[, c, drop = FALSE])
}

#' @rdname sym_det
#' @export
rowspace <- function(M) {
  r <- .rref(Matrix(M))
  lapply(seq_along(r$pivots), function(k) r$M[k, , drop = FALSE])
}

#' @rdname sym_det
#' @export
minor <- function(M, i, j) {
  if (inherits(M, "GeoEllipse")) return(.ellipse_minor(M))
  det(M[-i, -j, drop = FALSE])
}

#' @rdname sym_det
#' @export
cofactor <- function(M, i, j) {
  m <- minor(M, i, j)
  if ((i + j) %% 2L == 1L) .neg(m) else m
}

#' @rdname sym_det
#' @export
cofactor_matrix <- function(M) {
  M <- Matrix(M)
  Matrix(M$nrow, M$ncol, function(i, j) cofactor(M, i, j))
}

#' @rdname sym_det
#' @export
adjugate <- function(M) t(cofactor_matrix(M))

#' @rdname sym_det
#' @export
trace_matrix <- function(M) {
  M <- Matrix(M)
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
  M <- Matrix(M)
  cf <- .charpoly_coeffs(M)
  lambda <- as_sym(lambda)
  n <- length(cf) - 1L
  expand(.add(lapply(seq_along(cf), function(k) .mul(list(cf[[k]], .pow(lambda, .int(n - k + 1L)))))))
}

#' @rdname charpoly
#' @export
eigenvals <- function(M) {
  M <- Matrix(M)
  lam <- Dummy("lambda")
  roots(charpoly(M, lam), lam)
}

#' @rdname charpoly
#' @export
eigenvects <- function(M) {
  M <- Matrix(M)
  ev <- eigenvals(M)
  # sorted by the default_sort_key of the eigenvalue
  ord <- unlist(.key_sort(lapply(ev, .sort_key), as.list(seq_along(ev))))
  lapply(ord, function(k) {
    val <- ev[[k]]
    A <- M - .mat_map(eye(M$nrow), function(z) .mul(list(z, val)))
    ns <- nullspace(.mat_map(A, function(z) simplify(z)))
    list(value = val, multiplicity = attr(ev, "multiplicity")[k], vectors = ns)
  })
}

#' @rdname charpoly
#' @export
is_diagonalizable <- function(M) {
  ev <- eigenvects(M)
  sum(vapply(ev, function(e) length(e$vectors), 0L)) == Matrix(M)$nrow
}

#' @rdname charpoly
#' @export
diagonalize <- function(M) {
  M <- Matrix(M)
  if (M$nrow != M$ncol) stop(NonSquareMatrixError(""))
  ev <- eigenvects(M)
  vecs <- list()
  vals <- list()
  for (e in ev) {
    # eigenvects(simplify=True) divides each vector by the gcd of its entries
    if (e$multiplicity != length(e$vectors)) stop(MatrixError("Matrix is not diagonalizable"))
    for (v in e$vectors) {
      g <- Reduce(gcd, v$data)
      vecs[[length(vecs) + 1L]] <- .mat_map(v, function(z) simplify(.div(z, g)))
    }
    vals <- c(vals, rep(list(e$value), e$multiplicity))
  }
  P <- do.call(cbind.symmatrix, vecs)
  list(P = P, D = diag_matrix(vals))
}

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
  M <- Matrix(M)
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
#' `LUdecomposition()` returns `list(L, U, perm)` with `P*M = L*U`;
#' `QRdecomposition()` returns `list(Q, R)` by Gram-Schmidt;
#' `cholesky()` returns `L` with `M = L*t(L)` for symmetric positive definite
#' matrices; `LDLdecomposition()` returns `list(L, D)`. `GramSchmidt()`
#' orthogonalizes a list of vectors.
#'
#' @param M A `symmatrix`.
#' @param vlist A list of column vectors.
#' @param orthonormal Normalize the vectors.
#' @return A list of matrices.
#' @export
LUdecomposition <- function(M) {
  M <- Matrix(M)
  n <- M$nrow
  m <- M$ncol
  U <- .mrows(M)
  L <- .mrows(eye(n))
  perm <- seq_len(n)
  for (k in seq_len(min(n, m))) {
    if (is_zero_num(cancel(U[[k]][[k]]))) {
      sw <- NA
      if (k < n) for (i in (k + 1L):n) if (!is_zero_num(cancel(U[[i]][[k]]))) {
        sw <- i
        break
      }
      if (is.na(sw)) next
      tmp <- U[[k]]; U[[k]] <- U[[sw]]; U[[sw]] <- tmp
      if (k > 1L) for (j in seq_len(k - 1L)) {
        tmp <- L[[k]][[j]]; L[[k]][[j]] <- L[[sw]][[j]]; L[[sw]][[j]] <- tmp
      }
      tp <- perm[k]; perm[k] <- perm[sw]; perm[sw] <- tp
    }
    if (k < n) for (i in (k + 1L):n) {
      f <- cancel(.div(U[[i]][[k]], U[[k]][[k]]))
      L[[i]][[k]] <- f
      U[[i]] <- lapply(seq_len(m), function(j) cancel(.sub(U[[i]][[j]], .mul(list(f, U[[k]][[j]])))))
    }
  }
  list(L = .from_rows(L), U = .from_rows(U), perm = perm)
}

#' @rdname LUdecomposition
#' @export
GramSchmidt <- function(vlist, orthonormal = FALSE) {
  vlist <- lapply(vlist, Matrix)
  out <- list()
  for (v in vlist) {
    w <- v
    for (u in out) {
      coef <- .div(.dot(u, v), .dot(u, u))
      w <- w - .mat_map(u, function(z) .mul(list(z, coef)))
    }
    w <- .mat_map(w, simplify)
    if (all(vapply(w$data, is_zero_num, TRUE))) next
    out[[length(out) + 1L]] <- w
  }
  if (orthonormal) out <- lapply(out, function(w) {
    nrm <- .sqrt(.dot(w, w))
    .mat_map(w, function(z) simplify(.div(z, nrm)))
  })
  out
}

.dot <- function(a, b) .add(Map(function(x, y) .mul(list(fcall("conjugate", list(x)), y)), a$data, b$data))

#' @rdname LUdecomposition
#' @export
QRdecomposition <- function(M) {
  M <- Matrix(M)
  cols <- lapply(seq_len(M$ncol), function(j) M[, j, drop = FALSE])
  Q <- GramSchmidt(cols, orthonormal = TRUE)
  if (!length(Q)) return(list(Q = zeros(M$nrow, 0L), R = zeros(0L, M$ncol)))
  Qm <- do.call(cbind.symmatrix, Q)
  R <- .mat_map(.matmul(t(.mat_map(Qm, function(z) fcall("conjugate", list(z)))), M), simplify)
  list(Q = Qm, R = R)
}

#' @rdname LUdecomposition
#' @export
cholesky <- function(M) {
  M <- Matrix(M)
  n <- M$nrow
  L <- .mrows(zeros(n))
  for (i in seq_len(n)) {
    for (j in seq_len(i)) {
      s <- .add(lapply(seq_len(j - 1L), function(k) .mul(list(L[[i]][[k]], fcall("conjugate", list(L[[j]][[k]]))))))
      if (i == j) {
        L[[i]][[j]] <- simplify(.sqrt(.sub(.mget(M, i, i), s)))
      } else {
        L[[i]][[j]] <- simplify(.div(.sub(.mget(M, i, j), s), L[[j]][[j]]))
      }
    }
  }
  .from_rows(L)
}

#' @rdname LUdecomposition
#' @export
LDLdecomposition <- function(M) {
  M <- Matrix(M)
  n <- M$nrow
  L <- .mrows(eye(n))
  D <- rep(list(S_Zero), n)
  for (i in seq_len(n)) {
    for (j in seq_len(i)) {
      # Hermitian recurrence, as SymPy's default hermitian = TRUE
      s <- .add(lapply(seq_len(j - 1L), function(k) .mul(list(L[[i]][[k]], conjugate(L[[j]][[k]]), D[[k]]))))
      if (i == j) D[[i]] <- simplify(.sub(.mget(M, i, i), s))
      else L[[i]][[j]] <- simplify(.div(.sub(.mget(M, i, j), s), D[[j]]))
    }
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
#' @param M A matrix.
#' @param f A function applied to each entry.
#' @return An expression or matrix.
#' @export
dot <- function(a, b, doit = FALSE) {
  if (.pv_dispatch(a, b)) return(.pv_dot_dispatch(a, b))
  if (inherits(a, "Del")) return(.vec_divergence(b, doit = doit))
  if (.is_bd(a) || .is_bd(b)) return(.bd_dot_method(a, b))
  a <- Matrix(a)
  b <- Matrix(b)
  .add(Map(function(x, y) .mul(list(x, y)), a$data, b$data))
}

#' @rdname dot
#' @export
cross <- function(a, b, doit = FALSE) {
  if (.pv_dispatch(a, b)) return(.pv_cross_dispatch(a, b))
  if (inherits(a, "Del")) return(.vec_curl(b, doit = doit))
  if (.is_bd(a) || .is_bd(b)) return(.vec_cross_method(as_sym(a), b))
  a <- Matrix(a)$data
  b <- Matrix(b)$data
  if (length(a) != 3L || length(b) != 3L) stop("cross product needs 3-vectors")
  Matrix(list(.sub(.mul(list(a[[2]], b[[3]])), .mul(list(a[[3]], b[[2]]))),
              .sub(.mul(list(a[[3]], b[[1]])), .mul(list(a[[1]], b[[3]]))),
              .sub(.mul(list(a[[1]], b[[2]])), .mul(list(a[[2]], b[[1]])))))
}

#' @rdname dot
#' @export
norm_matrix <- function(M) {
  M <- Matrix(M)
  .sqrt(.add(lapply(M$data, function(z) .sq(fcall("Abs", list(z))))))
}

#' @rdname dot
#' @export
applyfunc <- function(M, f) {
  if (.is_matexpr(M) && !.is_explicit(M)) return(ElementwiseApplyFunction(f, M))
  .mat_map(Matrix(M), function(z) as_sym(if (is.function(f)) f(z) else apply_lambda(f, z)))
}

#' @rdname dot
#' @export
is_square <- function(M) {
  int_like <- (is_sym(M) && M$type == "Integer") || inherits(M, "bigz") ||
    (is.numeric(M) && length(M) == 1L && is.null(dim(M))) ||
    (is.character(M) && length(M) == 1L && grepl("^\\s*[-+]?[0-9]+\\s*$", M))
  if (int_like) return(.is_square_int(M))
  Matrix(M)$nrow == Matrix(M)$ncol
}

#' @rdname dot
#' @export
is_symmetric <- function(M) {
  M <- Matrix(M)
  is_square(M) && all(vapply(seq_along(M$data), function(k) {
    i <- (k - 1L) %/% M$ncol + 1L
    j <- (k - 1L) %% M$ncol + 1L
    is_zero_num(simplify(.sub(.mget(M, i, j), .mget(M, j, i))))
  }, TRUE))
}

#' @rdname dot
#' @export
is_zero_matrix <- function(M) all(vapply(Matrix(M)$data, function(z) is_zero_num(simplify(z)), TRUE))

#' @rdname dot
#' @export
vec <- function(M) {
  M <- Matrix(M)
  .new_mat(t(M)$data, M$nrow * M$ncol, 1L)
}

#' @rdname dot
#' @export
pinv <- function(M) {
  # _pinv_rank_decomposition of sympy/matrices/inverse.py
  M <- Matrix(M)
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
