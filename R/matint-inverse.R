# Row reduction and matrix inversion after sympy/matrices/reductions.py,
# sympy/matrices/determinant.py (_find_reasonable_pivot), decompositions.py
# (_LUdecomposition_Simple), solvers.py (_LUsolve) and inverse.py.

# _iszero: x.is_zero (TRUE, FALSE or NA)
.matint_iszero <- function(x) {
  if (is_zero_num(x)) return(TRUE)
  if (x$type == "Float") return(isTRUE(as.numeric(.neval(x, list())) == 0))
  if (.is_number(x) && is_number_type(x)) return(FALSE)
  z <- ask_fact(x, "zero")
  if (is.null(z) || length(z) != 1L) NA else z
}

.matint_call_iszero <- function(f, x) {
  z <- f(x)
  if (is.null(z) || length(z) != 1L) NA else as.logical(z)
}

# Expr.equals(0): TRUE, FALSE or NA
.matint_equals0 <- function(x) {
  r <- tryCatch(.expr_equals(x, S_Zero), error = function(e) NULL)
  if (is.null(r)) NA else r
}

# _find_reasonable_pivot: list(offset (0 based, NA for none), value,
# assumed_nonzero, newly_determined = list of list(offset, value))
.matint_find_pivot <- function(col, iszerofunc = .matint_iszero, simpfunc = simplify) {
  newly <- list()
  types <- vapply(col, function(z) z$type, "")
  if (length(col) && all(types %in% c("Float", "Integer")) && any(types == "Float")) {
    col_abs <- lapply(col, function(z) fcall("Abs", list(z)))
    av <- vapply(col, function(z) abs(as.numeric(.neval(z, list()))), 0)
    im <- which.max(av)
    max_value <- col_abs[[im]]
    if (isTRUE(.matint_call_iszero(iszerofunc, max_value))) {
      if (!is_zero_num(max_value))
        newly <- lapply(which(!vapply(col, is_zero_num, TRUE)), function(i) list(i - 1L, S_Zero))
      return(list(offset = NA_integer_, value = NULL, assumed = FALSE, newly = newly))
    }
    return(list(offset = im - 1L, value = col[[im]], assumed = FALSE, newly = newly))
  }
  possible <- logical(length(col))
  for (i in seq_along(col)) {
    z <- .matint_call_iszero(iszerofunc, col[[i]])
    if (isFALSE(z)) return(list(offset = i - 1L, value = col[[i]], assumed = FALSE, newly = newly))
    possible[i] <- z
  }
  if (all(possible %in% TRUE)) return(list(offset = NA_integer_, value = NULL, assumed = FALSE, newly = newly))
  for (i in seq_along(col)) {
    if (!is.na(possible[i])) next
    simped <- simpfunc(col[[i]])
    z <- .matint_call_iszero(iszerofunc, simped)
    if (!is.na(z)) newly[[length(newly) + 1L]] <- list(i - 1L, simped)
    if (isFALSE(z)) return(list(offset = i - 1L, value = simped, assumed = FALSE, newly = newly))
    possible[i] <- z
  }
  if (all(possible %in% TRUE)) return(list(offset = NA_integer_, value = NULL, assumed = FALSE, newly = newly))
  for (i in seq_along(col)) {
    if (!is.na(possible[i])) next
    if (isTRUE(.matint_equals0(col[[i]]))) {
      possible[i] <- TRUE
      newly[[length(newly) + 1L]] <- list(i - 1L, S_Zero)
    }
  }
  if (all(possible %in% TRUE)) return(list(offset = NA_integer_, value = NULL, assumed = FALSE, newly = newly))
  i <- which(is.na(possible))[1L]
  list(offset = i - 1L, value = col[[i]], assumed = TRUE, newly = newly)
}

# _find_reasonable_pivot_naive
.matint_find_pivot_naive <- function(col, iszerofunc = .matint_iszero, simpfunc = simplify) {
  ind <- list()
  for (i in seq_along(col)) {
    z <- .matint_call_iszero(iszerofunc, col[[i]])
    if (isFALSE(z)) return(list(offset = i - 1L, value = col[[i]], assumed = FALSE, newly = list()))
    if (is.na(z)) ind[[length(ind) + 1L]] <- list(i - 1L, col[[i]])
  }
  if (!length(ind)) return(list(offset = NA_integer_, value = NULL, assumed = FALSE, newly = list()))
  if (is.null(simpfunc)) return(list(offset = ind[[1L]][[1L]], value = ind[[1L]][[2L]], assumed = TRUE, newly = list()))
  newly <- list()
  for (p in ind) {
    tmp <- simpfunc(p[[2L]])
    if (!identical(tmp$key, p[[2L]]$key)) {
      newly[[length(newly) + 1L]] <- list(p[[1L]], tmp)
      if (isFALSE(.matint_call_iszero(iszerofunc, tmp)))
        return(list(offset = p[[1L]], value = tmp, assumed = FALSE, newly = newly))
    }
  }
  list(offset = ind[[1L]][[1L]], value = ind[[1L]][[2L]], assumed = TRUE, newly = newly)
}

# _row_reduce_list on a row major list of entries; isimp is dotprodsimp
.matint_row_reduce_list <- function(mat, rows, cols, iszerofunc, simpfunc, normalize_last = TRUE,
                                    normalize = TRUE, zero_above = TRUE) {
  ix <- function(i, j) i * cols + j + 1L
  isimp <- dotprodsimp
  row_swap <- function(i, j) {
    a <- ix(i, 0L):ix(i, cols - 1L); b <- ix(j, 0L):ix(j, cols - 1L)
    tmp <- mat[a]; mat[a] <<- mat[b]; mat[b] <<- tmp
  }
  cross_cancel <- function(a, i, b, j) {
    q <- (j - i) * cols
    for (p in ix(i, 0L):ix(i, cols - 1L))
      mat[[p]] <<- isimp(.sub(.mul(list(a, mat[[p]])), .mul(list(b, mat[[p + q]]))))
  }
  piv_row <- 0L; piv_col <- 0L
  pivot_cols <- integer(); swaps <- list()
  while (piv_col < cols && piv_row < rows) {
    col <- mat[ix(piv_row:(rows - 1L), piv_col)]
    fp <- .matint_find_pivot(col, iszerofunc, simpfunc)
    for (nd in fp$newly) mat[[ix(nd[[1L]] + piv_row, piv_col)]] <- nd[[2L]]
    if (is.na(fp$offset)) { piv_col <- piv_col + 1L; next }
    pivot_cols <- c(pivot_cols, piv_col)
    if (fp$offset != 0L) {
      row_swap(piv_row, fp$offset + piv_row)
      swaps[[length(swaps) + 1L]] <- c(piv_row, fp$offset + piv_row)
    }
    pivot_val <- fp$value
    if (!normalize_last || !.is_comm(pivot_val)) {
      i <- piv_row; j <- piv_col
      mat[[ix(i, j)]] <- S_One
      if (j + 1L < cols) for (p in ix(i, j + 1L):ix(i, cols - 1L))
        mat[[p]] <- isimp(.mul(list(.pow(pivot_val, S_NegativeOne), mat[[p]])))
      pivot_val <- S_One
    }
    for (row in seq_len(rows) - 1L) {
      if (row == piv_row) next
      if (!zero_above && row < piv_row) next
      val <- mat[[ix(row, piv_col)]]
      if (isTRUE(.matint_call_iszero(iszerofunc, val))) next
      cross_cancel(pivot_val, row, val, piv_row)
    }
    piv_row <- piv_row + 1L
  }
  if (normalize_last && normalize) {
    for (k in seq_along(pivot_cols)) {
      piv_i <- k - 1L; piv_j <- pivot_cols[k]
      pivot_val <- mat[[ix(piv_i, piv_j)]]
      mat[[ix(piv_i, piv_j)]] <- S_One
      if (piv_j + 1L < cols) for (p in ix(piv_i, piv_j + 1L):ix(piv_i, cols - 1L))
        mat[[p]] <- isimp(.mul(list(.pow(pivot_val, S_NegativeOne), mat[[p]])))
    }
  }
  list(mat = mat, pivots = pivot_cols, swaps = swaps)
}

# exact rref of a rational matrix (_rref_dm): the reduced form is unique
.matint_rref_qq <- function(M) {
  A <- .mrows(M); nr <- M$nrow; nc <- M$ncol
  pivots <- integer(); r <- 1L
  for (cc in seq_len(nc)) {
    if (r > nr) break
    p <- which(!vapply(r:nr, function(i) is_zero_num(A[[i]][[cc]]), TRUE))
    if (!length(p)) next
    pr <- r + p[1L] - 1L
    if (pr != r) { tmp <- A[[r]]; A[[r]] <- A[[pr]]; A[[pr]] <- tmp }
    pv <- A[[r]][[cc]]
    A[[r]] <- lapply(A[[r]], function(z) .div(z, pv))
    for (i in seq_len(nr)) {
      if (i == r) next
      f <- A[[i]][[cc]]
      if (is_zero_num(f)) next
      A[[i]] <- lapply(seq_len(nc), function(j) .sub(A[[i]][[j]], .mul(list(f, A[[r]][[j]]))))
    }
    pivots <- c(pivots, cc)
    r <- r + 1L
  }
  list(M = .from_rows(A), pivots = pivots)
}

# _rref: list(M = , pivots = ) with 1 based pivots
.matint_rref <- function(M, iszerofunc = .matint_iszero, simplify = FALSE, normalize_last = TRUE) {
  M <- Matrix(M)
  if (.matint_is_qq(M)) return(.matint_rref_qq(M))
  simpfunc <- if (is.function(simplify)) simplify else function(e) simplify(e)
  r <- .matint_row_reduce_list(M$data, M$nrow, M$ncol, iszerofunc, simpfunc, normalize_last)
  list(M = .new_mat(r$mat, M$nrow, M$ncol), pivots = r$pivots + 1L)
}

# _LUdecomposition_Simple: list(lu = row major entries, swaps)
.matint_lu_simple <- function(M, iszerofunc = .matint_iszero, simpfunc = NULL, rankcheck = FALSE) {
  rows <- M$nrow; cols <- M$ncol
  rank_err <- function() stop(.sympy_condition("ValueError", paste0(
    "Rank of matrix is strictly less than number of rows or columns. ",
    "Pass keyword argument rankcheck=False to compute the LU decomposition of this matrix.")))
  if (rows == 0L || cols == 0L) return(list(lu = zeros(rows, cols)$data, swaps = list()))
  lu <- M$data
  ix <- function(i, j) i * cols + j + 1L
  swaps <- list()
  pivot_col <- 0L
  if (rows >= 2L) for (pivot_row in 0:(rows - 2L)) {
    iszeropivot <- TRUE
    fp <- NULL
    while (pivot_col != cols && iszeropivot) {
      fp <- .matint_find_pivot_naive(lu[ix(pivot_row:(rows - 1L), pivot_col)], iszerofunc, simpfunc)
      iszeropivot <- is.null(fp$value)
      if (iszeropivot) pivot_col <- pivot_col + 1L
    }
    if (rankcheck && pivot_col != pivot_row) rank_err()
    cand <- if (is.na(fp$offset)) NA_integer_ else pivot_row + fp$offset
    if (is.na(cand) && iszeropivot) return(list(lu = lu, swaps = swaps))
    for (nd in fp$newly) lu[[ix(pivot_row + nd[[1L]], pivot_col)]] <- nd[[2L]]
    if (pivot_row != cand) {
      swaps[[length(swaps) + 1L]] <- c(pivot_row, cand)
      if (pivot_row > 0L) {
        a <- ix(pivot_row, 0L):ix(pivot_row, pivot_row - 1L); b <- ix(cand, 0L):ix(cand, pivot_row - 1L)
        tmp <- lu[a]; lu[a] <- lu[b]; lu[b] <- tmp
      }
      a <- ix(pivot_row, pivot_col):ix(pivot_row, cols - 1L); b <- ix(cand, pivot_col):ix(cand, cols - 1L)
      tmp <- lu[a]; lu[a] <- lu[b]; lu[b] <- tmp
    }
    start_col <- pivot_col + 1L
    for (row in (pivot_row + 1L):(rows - 1L)) {
      lu[[ix(row, pivot_row)]] <- .div(lu[[ix(row, pivot_col)]], lu[[ix(pivot_row, pivot_col)]])
      if (start_col < cols) for (cc in start_col:(cols - 1L))
        lu[[ix(row, cc)]] <- .sub(lu[[ix(row, cc)]], .mul(list(lu[[ix(row, pivot_row)]], lu[[ix(pivot_row, cc)]])))
    }
    if (pivot_row != pivot_col) for (row in (pivot_row + 1L):(rows - 1L)) lu[[ix(row, pivot_col)]] <- S_Zero
    pivot_col <- pivot_col + 1L
    if (pivot_col == cols) return(list(lu = lu, swaps = swaps))
  }
  if (rankcheck) {
    k <- min(rows, cols) - 1L
    if (isTRUE(.matint_call_iszero(iszerofunc, lu[[ix(k, k)]]))) rank_err()
  }
  list(lu = lu, swaps = swaps)
}

# _LUsolve
.matint_lusolve <- function(M, rhs, iszerofunc = .matint_iszero) {
  m <- M$nrow; n <- M$ncol
  if (rhs$nrow != m) stop(.sympy_condition("ShapeError", "``M`` and ``rhs`` must have the same number of rows."))
  if (m < n) stop(.sympy_condition("NotImplementedError", "Underdetermined systems not supported."))
  dec <- tryCatch(.matint_lu_simple(M, iszerofunc, rankcheck = TRUE), ValueError = function(e)
    stop(NonInvertibleMatrixError("Matrix det == 0; not invertible.")))
  A <- .new_mat(dec$lu, m, n)
  b <- .mrows(rhs)
  for (s in dec$swaps) { tmp <- b[[s[1L] + 1L]]; b[[s[1L] + 1L]] <- b[[s[2L] + 1L]]; b[[s[2L] + 1L]] <- tmp }
  for (i in seq_len(m)) for (j in seq_len(min(i - 1L, n))) {
    scale <- .mget(A, i, j)
    b[[i]] <- Map(function(x, y) .sub(x, .mul(list(scale, y))), b[[i]], b[[j]])
  }
  if (m > n) {
    for (i in (n + 1L):m) for (z in b[[i]]) if (!isTRUE(.matint_call_iszero(iszerofunc, z)))
      stop(.sympy_condition("ValueError", "The system is inconsistent."))
    b <- b[seq_len(n)]
  }
  for (i in rev(seq_len(n))) {
    if (i < n) for (j in (i + 1L):n) {
      scale <- .mget(A, i, j)
      b[[i]] <- Map(function(x, y) .sub(x, .mul(list(scale, y))), b[[i]], b[[j]])
    }
    scale <- .mget(A, i, i)
    b[[i]] <- lapply(b[[i]], function(x) .mul(list(.pow(scale, S_NegativeOne), x)))
  }
  .from_rows(b)
}

# _eval_get_diag_blocks
.matint_get_diag_blocks <- function(M) {
  out <- list()
  rec <- function(A) {
    n <- A$nrow
    for (i in seq_len(n)) {
      right <- if (i < A$ncol) A[seq_len(i), (i + 1L):A$ncol, drop = FALSE]$data else list()
      bottom <- if (i < n) A[(i + 1L):n, seq_len(i), drop = FALSE]$data else list()
      if (!all(vapply(c(right, bottom), is_zero_num, TRUE))) next
      out[[length(out) + 1L]] <<- A[seq_len(i), seq_len(i), drop = FALSE]
      if (i < n) rec(A[(i + 1L):n, (i + 1L):A$ncol, drop = FALSE])
      return(invisible())
    }
  }
  rec(M)
  out
}

# _inv_GE
.matint_inv_ge <- function(M, iszerofunc = .matint_iszero) {
  n <- M$nrow
  big <- cbind.symmatrix(M, eye(n))
  red <- .matint_rref(big, iszerofunc, simplify = TRUE)$M
  for (j in seq_len(red$nrow)) if (isTRUE(.matint_call_iszero(iszerofunc, .mget(red, j, j))))
    stop(NonInvertibleMatrixError("Matrix det == 0; not invertible."))
  red[, (n + 1L):(2L * n), drop = FALSE]
}

# _inv_ADJ: the berkowitz adjugate times 1/det, uncancelled
.matint_inv_adj <- function(M, iszerofunc = .matint_iszero) {
  d <- .det(M)
  if (isTRUE(.matint_call_iszero(iszerofunc, d))) stop(NonInvertibleMatrixError("Matrix det == 0; not invertible."))
  n <- M$nrow
  adj <- if (n == 1L) .new_mat(list(S_One), 1L, 1L) else Matrix(n, n, function(i, j) {
    m <- .det(M[-j, -i, drop = FALSE], method = "berkowitz")
    if ((i + j) %% 2L == 1L) .neg(m) else m
  })
  inv_d <- .pow(d, S_NegativeOne)
  .mat_map(adj, function(z) .mul(list(z, inv_d)))
}

# _inv_block
.matint_inv_block <- function(M, iszerofunc = .matint_iszero) {
  i <- M$nrow
  if (i <= 20L) return(.matint_inv(M, method = "LU"))
  h <- i %/% 2L
  A <- M[seq_len(h), seq_len(h), drop = FALSE]; B <- M[seq_len(h), (h + 1L):i, drop = FALSE]
  C <- M[(h + 1L):i, seq_len(h), drop = FALSE]; D <- M[(h + 1L):i, (h + 1L):i, drop = FALSE]
  lu <- function() .matint_inv(M, method = "LU")
  D_inv <- tryCatch(.matint_inv_block(D), NonInvertibleMatrixError = function(e) NULL)
  if (is.null(D_inv)) return(lu())
  B_D_i <- B %*% D_inv
  A_n <- tryCatch(.matint_inv_block(A - B_D_i %*% C), NonInvertibleMatrixError = function(e) NULL)
  if (is.null(A_n)) return(lu())
  B_n <- -(A_n %*% B_D_i)
  dc <- D_inv %*% C
  C_n <- -(dc %*% A_n)
  D_n <- D_inv + dc %*% (-B_n)
  rbind(cbind.symmatrix(A_n, B_n), cbind.symmatrix(C_n, D_n))
}

# _try_DM and _inv_DM: the fraction free inverse over the domain of the
# entries, NULL when that domain is EX
.matint_inv_dm <- function(M) {
  dM <- tryCatch(DomainMatrix_from_Matrix(M), error = function(e) NULL)
  if (is.null(dM)) return(NULL)
  K <- .dm_dom(dM)
  kind <- .de_kind(K)
  if (kind %in% c("EX", "EXRAW")) return(NULL)
  # ponytail: RR[x] and RR(x) would be converted to QQ[x] first in SymPy;
  # they take the GE route here
  if (!(kind %in% c("RR", "CC")) && grepl("RR|CC", format(K))) return(NULL)
  r <- tryCatch(dM$inv_den(), DMNonInvertibleMatrixError = function(e)
    stop(NonInvertibleMatrixError("Matrix det == 0; not invertible.")))
  dMi <- r[[1L]]; den <- r[[2L]]
  Kd <- K
  dMi <- dMi$to_field()
  .dm_div(dMi, den, Kd)$to_Matrix()
}

# _inv
.matint_inv <- function(M, method = NULL, iszerofunc = .matint_iszero, try_block_diag = FALSE) {
  if (M$nrow != M$ncol) stop(NonSquareMatrixError("A Matrix must be square to invert."))
  if (try_block_diag) {
    blocks <- .matint_get_diag_blocks(M)
    return(do.call(diag_matrix, lapply(blocks, function(b) .matint_inv(b, method, iszerofunc))))
  }
  if (is.null(method) && identical(iszerofunc, .matint_iszero)) {
    r <- .matint_inv_dm(M)
    if (!is.null(r)) return(r)
  }
  if (is.null(method)) method <- "GE"
  switch(method,
    GE = .matint_inv_ge(M, iszerofunc),
    LU = .matint_lusolve(M, eye(M$nrow), iszerofunc),
    ADJ = .matint_inv_adj(M, iszerofunc),
    BLOCK = .matint_inv_block(M, iszerofunc),
    stop(.sympy_condition("ValueError", "Inversion method unrecognized")))
}
