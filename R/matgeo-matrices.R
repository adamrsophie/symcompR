# Dense matrix methods after sympy/matrices/matrixbase.py, decompositions.py,
# solvers.py, reductions.py, subspaces.py, inverse.py, determinant.py and
# graph.py: matrix powers by Jordan blocks, structural predicates, row and
# column operations, the linear solvers and the remaining decompositions.

# bool(expr) of a SymPy entry: FALSE only for a numeric zero
.mg_bool <- function(e) {
  if (is_exact_number(e)) return(e$value != 0)
  if (e$type == "Float") return(!num_is_zero(e))
  if (identical(e$key, S_false$key)) return(FALSE)
  TRUE
}

# the structurally nonzero entries as list(i, j, v), row major (iter_items)
.mg_items <- function(M) {
  out <- list()
  for (i in seq_len(M$nrow)) for (j in seq_len(M$ncol)) {
    v <- .mget(M, i, j)
    if (.mg_bool(v)) out[[length(out) + 1L]] <- list(i = i, j = j, v = v)
  }
  out
}

.mg_iszero <- function(v) .matint_iszero(v)

.mg_simplify <- function(e) simplify(e)

.mg_all_zero <- function(M, keep) {
  for (it in .mg_items(M)) if (keep(it$i, it$j) && !isTRUE(.mg_iszero(it$v))) return(FALSE)
  TRUE
}

.mg_mat <- function(M) {
  if (.is_matexpr(M) && !.is_explicit(M)) return(M)
  .as_mat(M)
}

.mg_stop <- function(cls, msg = "") stop(.sympy_condition(cls, msg))

# ---------------------------------------------------------------------------
# Structural predicates

#' Structural tests of a matrix
#'
#' `is_upper()`, `is_lower()`, `is_diagonal()`, `is_hermitian()`,
#' `is_anti_symmetric()`, `is_upper_hessenberg()`, `is_lower_hessenberg()`,
#' `is_weakly_diagonally_dominant()`, `is_strongly_diagonally_dominant()`,
#' `is_Identity()`, `is_symbolic()` and `is_echelon()` follow the properties
#' and methods of SymPy's `MatrixBase`. The fuzzy tests return `NA` where
#' SymPy returns `None`.
#'
#' @param M A `symmatrix`.
#' @param simplify For `is_anti_symmetric()`: simplify the entries before
#'   testing (`TRUE`), or a function to do so.
#' @param iszerofunc For `is_echelon()`: the zero test.
#' @return `TRUE`, `FALSE` or `NA`.
#' @examples
#' M <- Matrix(list(list(1, 2), list(0, 3)))
#' is_upper(M)
#' is_lower(M)
#' is_diagonal(M)
#' @export
is_upper <- function(M) {
  M <- .as_mat(M)
  .mg_all_zero(M, function(i, j) i > j)
}

#' @rdname is_upper
#' @export
is_lower <- function(M) {
  M <- .as_mat(M)
  .mg_all_zero(M, function(i, j) i < j)
}

#' @rdname is_upper
#' @export
is_upper_hessenberg <- function(M) {
  M <- .as_mat(M)
  .mg_all_zero(M, function(i, j) i > j + 1L)
}

#' @rdname is_upper
#' @export
is_lower_hessenberg <- function(M) {
  M <- .as_mat(M)
  .mg_all_zero(M, function(i, j) i + 1L < j)
}

#' @rdname is_upper
#' @export
is_diagonal <- function(M) {
  M <- .as_mat(M)
  .matint_fuzzy_and(lapply(Filter(function(it) it$i != it$j, .mg_items(M)), function(it) .mg_iszero(it$v)))
}

#' @rdname is_upper
#' @export
is_hermitian <- function(M) {
  M <- .as_mat(M)
  if (M$nrow != M$ncol) return(FALSE)
  .matint_fuzzy_and(lapply(.mg_items(M), function(it)
    .matint_fact(simplify(.sub(.mget(M, it$i, it$j), conjugate(.mget(M, it$j, it$i)))), "zero")))
}

#' @rdname is_upper
#' @export
is_anti_symmetric <- function(M, simplify = TRUE) {
  M <- .as_mat(M)
  sf <- if (is.function(simplify)) simplify else if (isTRUE(simplify)) .mg_simplify else identity
  if (M$nrow != M$ncol) return(FALSE)
  .matint_fuzzy_and(lapply(.mg_items(M), function(it)
    .matint_fact(sf(.add(list(.mget(M, it$i, it$j), .mget(M, it$j, it$i)))), "zero")))
}

#' @rdname is_upper
#' @export
is_weakly_diagonally_dominant <- function(M) .matint_diag_dominant(.as_mat(M), "nonnegative")

#' @rdname is_upper
#' @export
is_strongly_diagonally_dominant <- function(M) .matint_diag_dominant(.as_mat(M), "positive")

#' @rdname is_upper
#' @export
is_Identity <- function(M) {
  M <- .as_mat(M)
  if (M$nrow != M$ncol) return(FALSE)
  for (i in seq_len(M$nrow)) for (j in seq_len(M$ncol)) {
    v <- .mget(M, i, j)
    if (i == j && !identical(v$key, S_One$key)) return(FALSE)
    if (i != j && !identical(v$key, S_Zero$key)) return(FALSE)
  }
  TRUE
}

#' @rdname is_upper
#' @export
is_symbolic <- function(M) any(vapply(.as_mat(M)$data, function(v) length(free_symbols(v)) > 0L, TRUE))

#' @rdname is_upper
#' @export
is_echelon <- function(M, iszerofunc = NULL) {
  M <- .as_mat(M)
  iz <- if (is.null(iszerofunc)) .mg_iszero else iszerofunc
  rec <- function(A) {
    if (A$nrow <= 0L || A$ncol <= 0L) return(TRUE)
    zeros_below <- all(vapply(seq_len(A$nrow - 1L) + 1L, function(i) isTRUE(iz(.mget(A, i, 1L))), TRUE))
    if (isTRUE(iz(.mget(A, 1L, 1L)))) {
      return(zeros_below && rec(A[, -1L, drop = FALSE]))
    }
    zeros_below && (A$nrow == 1L || A$ncol == 1L || rec(A[-1L, -1L, drop = FALSE]))
  }
  rec(M)
}

# ---------------------------------------------------------------------------
# Matrix powers

#' Power of a square matrix
#'
#' `matrix_pow(M, n, method)` is SymPy's `Matrix.pow()`: integer powers by
#' repeated squaring with intermediate simplification, symbolic and
#' rational powers from the Jordan normal form (`_matrix_pow_by_jordan_blocks`),
#' so a matrix that is not diagonalizable works too. `M %^% n` calls it.
#'
#' @param M A square `symmatrix`.
#' @param n The exponent.
#' @param method `NULL` (choose), `"multiply"`, `"mulsimp"`, `"jordan"` or
#'   `"cayley"`.
#' @return A `symmatrix`, or an unevaluated `MatPow`.
#' @examples
#' n <- Symbol("n")
#' matrix_pow(Matrix(list(list(1, 1), list(0, 1))), n)
#' Matrix(list(list(2, 1), list(0, 2))) %^% n
#' @export
matrix_pow <- function(M, n, method = NULL) {
  .keep_mclass((function() {
    if (!is.null(method) && !(method %in% c("multiply", "mulsimp", "jordan", "cayley")))
      .mg_stop("TypeError", "No such method")
    a <- .as_mat(M)
    if (a$nrow != a$ncol) stop(NonSquareMatrixError())
    e <- as_sym(n)
    if (isTRUE(.mg_iszero(e))) return(eye(a$nrow))
    if (identical(e$key, S_One$key)) return(a)
    if (isTRUE(is_diagonal(a))) {
      out <- a
      for (i in seq_len(a$nrow)) for (j in seq_len(a$ncol))
        out$data[[(i - 1L) * a$ncol + j]] <- if (i == j) .pow(.mget(a, i, j), e) else S_Zero
      return(out)
    }
    is_int_num <- (e$type %in% c("Integer", "Rational") && gmp::denominator(e$value) == 1) ||
      (e$type == "Float" && .mg_float_is_int(e))
    if (is_int_num) {
      if (a$nrow == 1L) return(.new_mat(list(.pow(.mget(a, 1L, 1L), e)), 1L, 1L))
      if (.mg_num(e) < 0) {
        e <- .neg(e)
        a <- inverse(a)
      }
    }
    if (identical(method, "jordan")) return(.mg_pow_jordan(a, e))
    if (identical(method, "cayley")) {
      if (e$type != "Integer") .mg_stop("ValueError", "cayley method is only valid for integer powers")
      return(.mg_pow_cayley(a, gmp::numerator(e$value)))
    }
    if (identical(method, "mulsimp")) {
      if (e$type != "Integer") .mg_stop("ValueError", "mulsimp method is only valid for integer powers")
      return(.mg_pow_dotprodsimp(a, gmp::numerator(e$value)))
    }
    if (identical(method, "multiply")) {
      if (e$type != "Integer") .mg_stop("ValueError", "multiply method is only valid for integer powers")
      return(.mg_pow_recursion(a, gmp::numerator(e$value)))
    }
    if (is.null(method) && e$type %in% c("Integer", "Float") && is_int_num) {
      k <- if (e$type == "Float") gmp::as.bigz(round(.mg_num(e))) else gmp::numerator(e$value)
      if (a$nrow == 2L && k > 100000) return(.mg_pow_jordan(a, .int(k)))
      return(.mg_pow_dotprodsimp(a, k))
    }
    r <- tryCatch(.mg_pow_jordan(a, e), NonInvertibleMatrixError = function(err) {
      if (isFALSE(.matint_fact(e, "integer")) || isFALSE(.matint_fact(e, "nonnegative"))) stop(err)
      NULL
    })
    if (!is.null(r)) return(r)
    MatPow(a, e)
  })(), M)
}

.mg_num <- function(e) as.numeric(.neval(e, list()))

.mg_float_is_int <- function(e) {
  v <- .mg_num(e)
  is.finite(v) && v == round(v)
}

.mg_pow_recursion <- function(a, k) {
  if (k == 1) return(a)
  if (k %% 2 == 1) {
    x <- a
    y <- .mg_pow_recursion(a, k - 1)
  } else {
    x <- y <- .mg_pow_recursion(a, k %/% 2)
  }
  .matmul(x, y)
}

# _eval_pow_by_recursion_dotprodsimp: prevsimp is shared through the recursion
.mg_pow_dotprodsimp <- function(a, k) {
  st <- new.env()
  st$prev <- rep(TRUE, length(a$data))
  rec <- function(k) {
    if (k == 1) return(a)
    if (k %% 2 == 1) {
      x <- a
      y <- rec(k - 1)
    } else {
      x <- y <- rec(k %/% 2)
    }
    m <- .matmul(x, y)
    for (i in seq_along(m$data)) if (st$prev[i]) {
      r <- dotprodsimp(m$data[[i]], withsimp = TRUE)
      m$data[[i]] <- r[[1L]]
      st$prev[i] <- isTRUE(r[[2L]])
    }
    m
  }
  rec(k)
}

# _eval_pow_by_cayley
.mg_pow_cayley <- function(a, k) {
  row <- a$nrow
  p <- charpoly(a)
  coeffs <- lapply(all_coeffs(-p)[-1L], as_sym)
  coeffs <- .mg_linrec_coeffs(coeffs, k)
  new_mat <- eye(row)
  ans <- zeros(row)
  for (i in seq_len(row)) {
    ans <- ans + .mat_map(new_mat, function(z) .mul(list(as_sym(coeffs[[i]]), z)))
    new_mat <- .matmul(new_mat, a)
  }
  ans
}

# _matrix_pow_by_jordan_blocks
.mg_pow_jordan <- function(A, num) {
  cell_power <- function(jc, n) {
    N <- jc$nrow
    l <- .mget(jc, 1L, 1L)
    set <- function(i, j, v) jc$data[[(i - 1L) * N + j]] <<- v
    if (isTRUE(.mg_iszero(l))) {
      if (N == 1L && isTRUE(.matint_fact(n, "nonnegative"))) {
        set(1L, 1L, .pow(l, n))
      } else if (!(isTRUE(.matint_fact(n, "integer")) && isTRUE(.matint_fact(n, "nonnegative")))) {
        stop(NonInvertibleMatrixError("Non-invertible matrix can only be raised to a nonnegative integer"))
      } else {
        for (i in seq_len(N)) set(1L, i, KroneckerDelta(.int(i - 1L), n))
      }
    } else {
      for (i in seq_len(N)) {
        bn <- binomial(n, .int(i - 1L))
        if (is_sym(bn) && bn$type == "Function" && identical(bn$name, "binomial")) bn <- expand_func(bn)
        set(1L, i, .mul(list(.pow(l, .sub(n, .int(i - 1L))), bn)))
      }
    }
    for (i in seq_len(N)) for (j in seq_len(N - i)) set(j + 1L, i + j, .mget(jc, j, i + j - 1L))
    jc
  }
  jf <- jordan_form(A)
  cells <- lapply(.matint_get_diag_blocks(jf$J), cell_power, n = num)
  .matmul(.matmul(jf$P, do.call(diag_matrix, cells)), inverse(jf$P))
}

# linrec_coeffs of sympy/discrete/recurrences.py
.mg_linrec_coeffs <- function(cf, n) {
  k <- length(cf)
  sq_reduce <- function(u, offset) {
    w <- rep(list(S_Zero), 2L * length(u) - 1L + offset)
    for (i in seq_along(u)) for (j in seq_along(u))
      w[[offset + i + j - 1L]] <- .add(list(w[[offset + i + j - 1L]], .mul(list(u[[i]], u[[j]]))))
    if (length(w) > k) for (j in seq(length(w), k + 1L)) for (i in seq_len(k))
      w[[j - i]] <- .add(list(w[[j - i]], .mul(list(w[[j]], cf[[i]]))))
    w[seq_len(k)]
  }
  fin <- function(n) {
    if (n < k) return(c(rep(list(S_Zero), as.integer(n)), list(S_One), rep(list(S_Zero), k - as.integer(n) - 1L)))
    sq_reduce(fin(n %/% 2), as.integer(n %% 2))
  }
  fin(gmp::as.bigz(n))
}

# ---------------------------------------------------------------------------
# Rows, columns and reshaping

.mg_new <- function(nr, nc, entry) {
  d <- vector("list", nr * nc)
  for (i in seq_len(nr)) for (j in seq_len(nc)) d[[(i - 1L) * nc + j]] <- as_sym(entry(i, j))
  .new_mat(d, nr, nc)
}

#' Row and column operations
#'
#' Functional versions of SymPy's `Matrix` row and column methods. Indices
#' are 1 based, as everywhere in R. `row_insert()` and `col_insert()` insert
#' the rows or columns of `other` before position `pos`; `row_del()` and
#' `col_del()` delete one; `row_join()` and `col_join()` concatenate;
#' `extract()` picks rows and columns; `permute()`, `permute_rows()` and
#' `permute_cols()` apply a permutation given as a list of swaps or as an
#' image vector (a SymPy `Permutation` is taken 0 based); `elementary_row_op()`
#' and `elementary_col_op()` apply `"n->kn"`, `"n<->m"` or `"n->n+km"`;
#' `row_swap()` and `col_swap()` swap two rows or columns. `diagonal()`,
#' `upper_triangular()` and `lower_triangular()` take the `k`th diagonal and
#' the triangles; `vech()` stacks the lower triangle; `rot90()` rotates;
#' `mat_row()`, `mat_col()`, `tolist()` and `matrix_values()` read entries;
#' `multiply_elementwise()` is the Hadamard product.
#'
#' @param M A `symmatrix`.
#' @param pos,i,j,i1,i2 Positions (1 based).
#' @param other A `symmatrix`.
#' @param rows,cols Row and column indices, or logical vectors.
#' @param perm A list of swaps, an image vector or a `Permutation`.
#' @param orientation `"rows"` or `"cols"`.
#' @param direction `"forward"` or `"backward"`.
#' @param op The elementary operation.
#' @param row,col,row1,row2,col1,col2 The rows or columns it acts on.
#' @param k A scale factor, or for `diagonal()` and the triangles the offset
#'   of the diagonal; for `rot90()` the number of quarter turns.
#' @param diagonal For `vech()`: include the diagonal.
#' @param check_symmetry For `vech()`: require a symmetric matrix.
#' @return A `symmatrix` (`tolist()` and `matrix_values()` return lists).
#' @examples
#' M <- Matrix(list(list(1, 2), list(3, 4)))
#' row_insert(M, 2, Matrix(list(list(0, 0))))
#' col_del(M, 1)
#' permute_rows(M, list(list(1, 2)))
#' elementary_row_op(M, "n->n+km", row = 2, k = -3, row2 = 1)
#' diagonal(M)
#' @export
row_insert <- function(M, pos, other) {
  .keep_mclass((function() {
    M <- .as_mat(M); other <- .as_mat(other)
    if (!length(M$data)) return(other)
    pos <- as.integer(pos) - 1L
    if (pos < 0L) pos <- M$nrow + pos + 1L
    if (pos < 0L) pos <- 0L else if (pos > M$nrow) pos <- M$nrow
    if (M$ncol != other$ncol)
      .mg_stop("ShapeError", sprintf("The matrices have incompatible number of columns (%d and %d)", M$ncol, other$ncol))
    ins <- pos * M$ncol
    .new_mat(append(M$data, other$data, after = ins), M$nrow + other$nrow, M$ncol)
  })(), M)
}

#' @rdname row_insert
#' @export
col_insert <- function(M, pos, other) {
  .keep_mclass((function() {
    M <- .as_mat(M); other <- .as_mat(other)
    if (!length(M$data)) return(other)
    pos <- as.integer(pos) - 1L
    if (pos < 0L) pos <- M$ncol + pos + 1L
    if (pos < 0L) pos <- 0L else if (pos > M$ncol) pos <- M$ncol
    if (M$nrow != other$nrow)
      .mg_stop("ShapeError", sprintf("The matrices have incompatible number of rows (%d and %d)", M$nrow, other$nrow))
    .mg_new(M$nrow, M$ncol + other$ncol, function(i, j) {
      if (j <= pos) .mget(M, i, j) else if (j <= pos + other$ncol) .mget(other, i, j - pos) else .mget(M, i, j - other$ncol)
    })
  })(), M)
}

#' @rdname row_insert
#' @export
row_del <- function(M, i) {
  .keep_mclass((function() {
    M <- .as_mat(M)
    i <- as.integer(i)
    if (i < 0L) i <- i + M$nrow + 1L
    if (i < 1L || i > M$nrow) .mg_stop("IndexError", sprintf("Row %d is out of range.", i))
    M[-i, , drop = FALSE]
  })(), M)
}

#' @rdname row_insert
#' @export
col_del <- function(M, j) {
  .keep_mclass((function() {
    M <- .as_mat(M)
    j <- as.integer(j)
    if (j < 0L) j <- j + M$ncol + 1L
    if (j < 1L || j > M$ncol) .mg_stop("IndexError", sprintf("Column %d is out of range.", j))
    M[, -j, drop = FALSE]
  })(), M)
}

#' @rdname row_insert
#' @export
row_join <- function(M, other) {
  .keep_mclass((function() {
    M <- .as_mat(M); other <- .as_mat(other)
    if (M$ncol == 0L && M$nrow != other$nrow) return(other)
    if (M$nrow != other$nrow)
      .mg_stop("ShapeError", sprintf("The matrices have incompatible number of rows (%d and %d)", M$nrow, other$nrow))
    cbind.symmatrix(M, other)
  })(), M)
}

#' @rdname row_insert
#' @export
col_join <- function(M, other) {
  .keep_mclass((function() {
    M <- .as_mat(M); other <- .as_mat(other)
    if (M$nrow == 0L && M$ncol != other$ncol) return(other)
    if (M$ncol != other$ncol)
      .mg_stop("ShapeError", sprintf("The matrices have incompatible number of columns (%d and %d)", M$ncol, other$ncol))
    rbind.symmatrix(M, other)
  })(), M)
}

#' @rdname row_insert
#' @export
extract <- function(M, rows, cols) {
  M <- .as_mat(M)
  idx <- function(v, n) {
    v <- unlist(v)
    if (is.logical(v)) return(which(v))
    v <- as.integer(v)
    v[v < 0L] <- v[v < 0L] + n + 1L
    if (any(v < 1L | v > n)) .mg_stop("IndexError", "Index out of range")
    v
  }
  ri <- idx(rows, M$nrow); ci <- idx(cols, M$ncol)
  .mg_new(length(ri), length(ci), function(i, j) .mget(M, ri[i], ci[j]))
}

#' @rdname row_insert
#' @export
permute <- function(M, perm, orientation = "rows", direction = "forward") {
  M <- .as_mat(M)
  if (identical(direction, "forwards")) direction <- "forward"
  if (identical(direction, "backwards")) direction <- "backward"
  if (identical(orientation, "columns")) orientation <- "cols"
  if (!(direction %in% c("forward", "backward")))
    .mg_stop("TypeError", sprintf("direction='%s' is an invalid kwarg. Try 'forward' or 'backward'", direction))
  if (!(orientation %in% c("rows", "cols")))
    .mg_stop("TypeError", sprintf("orientation='%s' is an invalid kwarg. Try 'rows' or 'cols'", orientation))
  max_index <- if (orientation == "rows") M$nrow else M$ncol
  if (inherits(perm, "Permutation")) {
    p <- Permutation(perm, size = max_index + 1L)
  } else {
    flat <- as.integer(unlist(perm))
    if (!all(flat >= 1L & flat <= max_index + 1L)) .mg_stop("IndexError", "`swap` indices out of range.")
    if (length(perm) && is.list(perm) && length(perm[[1L]]) > 1L || (is.list(perm) && length(perm) && is.list(perm[[1L]]))) {
      sw <- lapply(perm, function(s) as.integer(unlist(s)) - 1L)
      if (direction == "forward") sw <- rev(sw)
      p <- Permutation(sw, size = max_index + 1L)
    } else {
      p <- Permutation(flat - 1L, size = max_index + 1L)
    }
  }
  mapping <- array_form(p) + 1L
  if (orientation == "rows") .mg_new(M$nrow, M$ncol, function(i, j) .mget(M, mapping[i], j))
  else .mg_new(M$nrow, M$ncol, function(i, j) .mget(M, i, mapping[j]))
}

#' @rdname row_insert
#' @export
permute_rows <- function(M, perm, direction = "forward") permute(M, perm, "rows", direction)

#' @rdname row_insert
#' @export
permute_cols <- function(M, perm, direction = "forward") permute(M, perm, "cols", direction)

.mg_op_args <- function(M, op, col, k, col1, col2, what) {
  if (!(op %in% c("n->kn", "n<->m", "n->n+km")))
    .mg_stop("ValueError", sprintf("Unknown %s operation '%s'. Valid col operations are 'n->kn', 'n<->m', 'n->n+km'", what, op))
  n <- if (what == "col") M$ncol else M$nrow
  bad <- function(c) .mg_stop("ValueError", sprintf("This matrix does not have a %s '%s'", what, c))
  if (op == "n->kn") {
    if (is.null(col)) col <- col1
    if (is.null(col) || is.null(k))
      .mg_stop("ValueError", sprintf("For a %s operation 'n->kn' you must provide the kwargs `%s` and `k`", what, what))
    if (!(col >= 1 && col <= n)) bad(col)
  } else if (op == "n<->m") {
    cs <- unique(unlist(list(col, if (is.numeric(k)) k, col1, col2)))
    if (length(cs) > 2L) cs <- unique(unlist(list(col, col1, col2)))
    if (length(cs) != 2L)
      .mg_stop("ValueError", sprintf("For a %s operation 'n<->m' you must provide the kwargs `%s1` and `%s2`", what, what, what))
    col1 <- cs[1L]; col2 <- cs[2L]
    if (!(col1 >= 1 && col1 <= n)) bad(col1)
    if (!(col2 >= 1 && col2 <= n)) bad(col2)
  } else {
    if (is.null(col)) col <- col1
    if (is.null(col2)) col2 <- col1
    if (is.null(col) || is.null(col2) || is.null(k))
      .mg_stop("ValueError", sprintf("For a %s operation 'n->n+km' you must provide the kwargs `%s`, `k`, and `%s2`", what, what, what))
    if (col == col2)
      .mg_stop("ValueError", sprintf("For a %s operation 'n->n+km' `%s` and `%s2` must be different.", what, what, what))
    if (!(col >= 1 && col <= n)) bad(col)
    if (!(col2 >= 1 && col2 <= n)) bad(col2)
  }
  list(col = col, k = k, col1 = col1, col2 = col2)
}

#' @rdname row_insert
#' @export
elementary_row_op <- function(M, op = "n->kn", row = NULL, k = NULL, row1 = NULL, row2 = NULL) {
  M <- .as_mat(M)
  a <- .mg_op_args(M, op, row, k, row1, row2, "row")
  switch(op,
    "n->kn" = .mg_new(M$nrow, M$ncol, function(i, j) if (i == a$col) .mul(list(as_sym(a$k), .mget(M, i, j))) else .mget(M, i, j)),
    "n<->m" = .mg_new(M$nrow, M$ncol, function(i, j)
      if (i == a$col1) .mget(M, a$col2, j) else if (i == a$col2) .mget(M, a$col1, j) else .mget(M, i, j)),
    .mg_new(M$nrow, M$ncol, function(i, j)
      if (i == a$col) .add(list(.mget(M, i, j), .mul(list(as_sym(a$k), .mget(M, a$col2, j))))) else .mget(M, i, j)))
}

#' @rdname row_insert
#' @export
elementary_col_op <- function(M, op = "n->kn", col = NULL, k = NULL, col1 = NULL, col2 = NULL) {
  M <- .as_mat(M)
  a <- .mg_op_args(M, op, col, k, col1, col2, "col")
  switch(op,
    "n->kn" = .mg_new(M$nrow, M$ncol, function(i, j) if (j == a$col) .mul(list(as_sym(a$k), .mget(M, i, j))) else .mget(M, i, j)),
    "n<->m" = .mg_new(M$nrow, M$ncol, function(i, j)
      if (j == a$col1) .mget(M, i, a$col2) else if (j == a$col2) .mget(M, i, a$col1) else .mget(M, i, j)),
    .mg_new(M$nrow, M$ncol, function(i, j)
      if (j == a$col) .add(list(.mget(M, i, j), .mul(list(as_sym(a$k), .mget(M, i, a$col2))))) else .mget(M, i, j)))
}

#' @rdname row_insert
#' @export
row_swap <- function(M, i1, i2) elementary_row_op(M, "n<->m", row1 = i1, row2 = i2)

#' @rdname row_insert
#' @export
col_swap <- function(M, i1, i2) elementary_col_op(M, "n<->m", col1 = i1, col2 = i2)

#' @rdname row_insert
#' @export
diagonal <- function(M, k = 0) {
  M <- .as_mat(M)
  k <- as.integer(k)
  r <- if (max(0L, -k) < min(M$nrow, M$ncol - k)) seq(max(0L, -k), min(M$nrow, M$ncol - k) - 1L) else integer()
  .new_mat(lapply(r, function(i) .mget(M, i + 1L, i + k + 1L)), 1L, length(r))
}

#' @rdname row_insert
#' @export
upper_triangular <- function(M, k = 0) {
  M <- .as_mat(M)
  .mg_new(M$nrow, M$ncol, function(i, j) if (i + k <= j) .mget(M, i, j) else S_Zero)
}

#' @rdname row_insert
#' @export
lower_triangular <- function(M, k = 0) {
  M <- .as_mat(M)
  .mg_new(M$nrow, M$ncol, function(i, j) if (i + k >= j) .mget(M, i, j) else S_Zero)
}

#' @rdname row_insert
#' @export
vech <- function(M, diagonal = TRUE, check_symmetry = TRUE) {
  M <- .as_mat(M)
  if (M$nrow != M$ncol) stop(NonSquareMatrixError())
  if (check_symmetry && !isTRUE(is_symmetric(M))) .mg_stop("ValueError", "The matrix is not symmetric.")
  n <- M$ncol
  v <- list()
  for (j in seq_len(n)) for (i in seq_len(n)) if (if (diagonal) i >= j else i > j) v[[length(v) + 1L]] <- .mget(M, i, j)
  .new_mat(v, length(v), 1L)
}

#' @rdname row_insert
#' @export
rot90 <- function(M, k = 1) {
  M <- .as_mat(M)
  switch(as.character(k %% 4),
    "0" = M,
    "1" = t(M[rev(seq_len(M$nrow)), , drop = FALSE]),
    "2" = M[rev(seq_len(M$nrow)), rev(seq_len(M$ncol)), drop = FALSE],
    "3" = t(M[, rev(seq_len(M$ncol)), drop = FALSE]))
}

#' @rdname row_insert
#' @export
mat_row <- function(M, i) .as_mat(M)[i, , drop = FALSE]

#' @rdname row_insert
#' @export
mat_col <- function(M, j) .as_mat(M)[, j, drop = FALSE]

#' @rdname row_insert
#' @export
tolist <- function(M) .mrows(.as_mat(M))

#' @rdname row_insert
#' @export
matrix_values <- function(M) lapply(.mg_items(.as_mat(M)), function(it) it$v)

#' @rdname row_insert
#' @export
multiply_elementwise <- function(M, other) {
  M <- .as_mat(M); other <- .as_mat(other)
  if (M$nrow != other$nrow || M$ncol != other$ncol)
    .mg_stop("ShapeError", sprintf("Matrix shapes must agree (%d, %d) != (%d, %d)", M$nrow, M$ncol, other$nrow, other$ncol))
  .new_mat(Map(function(a, b) .mul(list(a, b)), M$data, other$data), M$nrow, M$ncol)
}

#' Special matrices
#'
#' `jordan_block(size, eigenvalue, band)` is a Jordan block,
#' `companion(p)` the companion matrix of a monic univariate `Poly`,
#' `wilkinson(n)` the pair of Wilkinson matrices `list(W_minus, W_plus)`, and
#' `get_diag_blocks(M)` the square blocks along the diagonal, as SymPy's
#' `Matrix` class methods.
#'
#' @param size The size of the block.
#' @param eigenvalue The diagonal entry.
#' @param band `"upper"` or `"lower"`: where the ones go.
#' @param eigenval An alias of `eigenvalue`.
#' @param p A monic univariate `Poly`.
#' @param n The Wilkinson order: the matrices are `2n + 1` square.
#' @param M A `symmatrix`.
#' @return A `symmatrix`, or a list of them.
#' @examples
#' jordan_block(3, 2)
#' x <- Symbol("x")
#' companion(Poly(x^2 - 3*x + 2, x))
#' get_diag_blocks(diag_matrix(1, Matrix(list(list(1, 2), list(3, 4)))))
#' @export
jordan_block <- function(size = NULL, eigenvalue = NULL, band = "upper", eigenval = NULL) {
  if (!is.null(eigenvalue)) {
    if (!is.null(eigenval) && !identical(as_sym(eigenvalue)$key, as_sym(eigenval)$key))
      .mg_stop("ValueError", "Cannot supply both 'eigenvalue' and 'eigenval'")
    e <- as_sym(eigenvalue)
  } else if (!is.null(eigenval)) {
    e <- as_sym(eigenval)
  } else .mg_stop("ValueError", "Must supply an eigenvalue")
  if (is.null(size)) .mg_stop("ValueError", "Must supply a matrix size")
  size <- as.integer(size)
  if (band == "lower") .mg_new(size, size, function(i, j) if (i == j) e else if (j + 1L == i) S_One else S_Zero)
  else .mg_new(size, size, function(i, j) if (i == j) e else if (i + 1L == j) S_One else S_Zero)
}

#' @rdname jordan_block
#' @export
companion <- function(p) {
  if (!inherits(p, "Poly"))
    .mg_stop("ValueError", sprintf("%s must be a Poly instance.", sstr(p)))
  cf <- lapply(all_coeffs(p), as_sym)
  if (!identical(cf[[1L]]$key, S_One$key)) .mg_stop("ValueError", sprintf("%s must be a monic polynomial.", sstr(p)))
  if (.poly_ngens(p) != 1L) .mg_stop("ValueError", sprintf("%s must be a univariate polynomial.", sstr(p)))
  size <- length(cf) - 1L
  if (size < 1L) .mg_stop("ValueError", sprintf("%s must have degree not less than 1.", sstr(p)))
  .mg_new(size, size, function(i, j) {
    if (j == size) .neg(cf[[length(cf) - i + 1L]]) else if (i == j + 1L) S_One else S_Zero
  })
}

#' @rdname jordan_block
#' @export
wilkinson <- function(n) {
  n <- as.integer(n)
  N <- 2L * n + 1L
  D <- .mg_new(N, N, function(i, j) if (i + 1L == j) S_One else S_Zero)
  dg <- do.call(diag_matrix, lapply(seq(-n, n), .int))
  dga <- do.call(diag_matrix, lapply(abs(seq(-n, n)), .int))
  list(dg + D + t(D), dga + D + t(D))
}

#' @rdname jordan_block
#' @export
get_diag_blocks <- function(M) .matint_get_diag_blocks(.as_mat(M))

# ---------------------------------------------------------------------------
# Decompositions

# is_zero_matrix as a fuzzy property: fuzzy_and of the entries' is_zero
.mg_fuzzy_zero_matrix <- function(M) .matint_fuzzy_and(lapply(M$data, .mg_iszero))

# Matrix.norm() of a vector: sqrt of the sum of |v|^2 over the nonzero entries
.mg_vnorm <- function(v) {
  vals <- Filter(.mg_bool, v$data)
  if (!length(vals)) vals <- list(S_Zero)
  .sqrt(.add(lapply(vals, function(z) .pow(Abs(z), .int(2L)))))
}

# u.dot(v, hermitian = TRUE): conjugate the first vector
.mg_hdot <- function(u, v) .add(Map(function(a, b) .mul(list(conjugate(a), b)), u$data, v$data))

# _QRdecomposition_optional
.mg_qr <- function(M, normalize = TRUE) {
  dps <- function(e) expand_mul(e)
  A <- M
  nr <- A$nrow; nc <- A$ncol
  Q <- A
  R <- zeros(nc)
  colof <- function(X, j) X[, j, drop = FALSE]
  setcol <- function(X, j, v) { for (i in seq_len(nr)) X$data[[(i - 1L) * nc + j]] <- v$data[[i]]; X }
  setR <- function(i, j, v) R$data[[(i - 1L) * nc + j]] <<- v
  ranked <- integer()
  for (j in seq_len(nc)) {
    for (i in seq_len(j - 1L)) {
      qi <- colof(Q, i)
      z <- .mg_fuzzy_zero_matrix(qi)
      if (isTRUE(z)) next
      nrm <- dps(.mg_hdot(qi, qi))
      if (is.na(z) && isTRUE(.matint_equals0(nrm))) next
      rij <- dps(.div(.mg_hdot(qi, colof(Q, j)), nrm))
      setR(i, j, rij)
      Q <- setcol(Q, j, .mat_map(colof(Q, j) - .mat_map(qi, function(e) .mul(list(e, rij))), identity))
    }
    Q <- setcol(Q, j, .mat_map(colof(Q, j), dps))
    qj <- colof(Q, j)
    z <- .mg_fuzzy_zero_matrix(qj)
    if (isTRUE(z) || isTRUE(.matint_equals0(dps(.mg_hdot(qj, qj))))) next
    ranked <- c(ranked, j)
    setR(j, j, S_One)
  }
  Q <- Q[seq_len(nr), ranked, drop = FALSE]
  R <- R[ranked, seq_len(nc), drop = FALSE]
  if (normalize) for (i in seq_len(Q$ncol)) {
    nrm <- .mg_vnorm(Q[, i, drop = FALSE])
    for (r in seq_len(nr)) Q$data[[(r - 1L) * Q$ncol + i]] <- .div(Q$data[[(r - 1L) * Q$ncol + i]], nrm)
    for (c in seq_len(R$ncol)) R$data[[(i - 1L) * R$ncol + c]] <- .mul(list(R$data[[(i - 1L) * R$ncol + c]], nrm))
  }
  list(Q = Q, R = R)
}

# _orthogonalize
.mg_orthogonalize <- function(vecs, normalize = FALSE, rankcheck = FALSE) {
  if (!length(vecs)) return(list())
  vecs <- lapply(vecs, Matrix)
  all_row <- vecs[[1L]]$nrow == 1L
  M <- do.call(cbind.symmatrix, lapply(vecs, vec))
  qr <- .mg_qr(M, normalize)
  if (rankcheck && qr$Q$ncol < length(vecs))
    .mg_stop("ValueError", "GramSchmidt: vector set not linearly independent")
  lapply(seq_len(qr$Q$ncol), function(i) {
    col <- qr$Q[, i, drop = FALSE]
    if (all_row) t(col) else col
  })
}

.mg_lu <- function(M, iszerofunc = NULL, simpfunc = NULL, rankcheck = FALSE) {
  M <- .as_mat(M)
  iz <- if (is.null(iszerofunc)) .matint_iszero else iszerofunc
  r <- .matint_lu_simple(M, iz, simpfunc, rankcheck)
  list(lu = .new_mat(r$lu, M$nrow, M$ncol), swaps = lapply(r$swaps, function(s) as.integer(s) + 1L))
}

#' More matrix decompositions
#'
#' `LUdecomposition_Simple()` returns `list(LU, perm)`: the combined `L`
#' and `U` factors and the row swaps; `LUdecompositionFF()` the fraction
#' free `list(P, L, DD, U)` with `P*M = L*DD^-1*U`;
#' `rank_decomposition()` `list(C, F)` with `M = C*F`;
#' `singular_value_decomposition()` `list(U, S, V)` with `M = U*S*V^H`;
#' `upper_hessenberg_decomposition()` `list(H, P)` with `M = P*H*P^H`;
#' `bidiagonal_decomposition()` `list(U, B, V)` and `bidiagonalize()` the
#' bidiagonal form, all as SymPy's `Matrix` methods. Row swaps are 1 based.
#'
#' @param M A `symmatrix`.
#' @param iszerofunc The zero test (`x$is_zero` by default).
#' @param simpfunc An optional simplification function for pivot search.
#' @param rankcheck Raise when the matrix is rank deficient.
#' @param simplify For `rank_decomposition()`: simplify during row reduction.
#' @param upper For the bidiagonal forms: upper (`TRUE`) or lower.
#' @return A list of matrices, or a matrix.
#' @examples
#' M <- Matrix(list(list(1, 2), list(3, 4)))
#' LUdecomposition_Simple(M)
#' rank_decomposition(Matrix(list(list(1, 2), list(2, 4))))
#' upper_hessenberg_decomposition(Matrix(list(list(1, 2, 3), list(4, 5, 6), list(7, 8, 9))))$H
#' @export
LUdecomposition_Simple <- function(M, iszerofunc = NULL, simpfunc = NULL, rankcheck = FALSE) {
  r <- .mg_lu(M, iszerofunc, simpfunc, rankcheck)
  list(r$lu, r$swaps)
}

#' @rdname LUdecomposition_Simple
#' @export
LUdecompositionFF <- function(M) {
  M <- .as_mat(M)
  n <- M$nrow; m <- M$ncol
  U <- .mrows(M); L <- .mrows(eye(n)); P <- .mrows(eye(n))
  DD <- .mrows(zeros(n, n))
  oldpivot <- S_One
  for (k in seq_len(n - 1L)) {
    if (!.mg_bool(U[[k]][[k]])) {
      kp <- NA
      for (kk in seq(k + 1L, n)) if (.mg_bool(U[[kk]][[k]])) { kp <- kk; break }
      if (is.na(kp)) .mg_stop("ValueError", "Matrix is not full rank")
      for (j in k:m) { tmp <- U[[k]][[j]]; U[[k]][[j]] <- U[[kp]][[j]]; U[[kp]][[j]] <- tmp }
      for (j in seq_len(k - 1L)) { tmp <- L[[k]][[j]]; L[[k]][[j]] <- L[[kp]][[j]]; L[[kp]][[j]] <- tmp }
      tmp <- P[[k]]; P[[k]] <- P[[kp]]; P[[kp]] <- tmp
    }
    Ukk <- U[[k]][[k]]
    L[[k]][[k]] <- Ukk
    DD[[k]][[k]] <- .mul(list(oldpivot, Ukk))
    for (i in seq(k + 1L, n)) {
      Uik <- U[[i]][[k]]
      L[[i]][[k]] <- Uik
      if (k + 1L <= m) for (j in seq(k + 1L, m))
        U[[i]][[j]] <- .div(.sub(.mul(list(Ukk, U[[i]][[j]])), .mul(list(U[[k]][[j]], Uik))), oldpivot)
      U[[i]][[k]] <- S_Zero
    }
    oldpivot <- Ukk
  }
  DD[[n]][[n]] <- oldpivot
  list(P = .from_rows(P), L = .from_rows(L), DD = .from_rows(DD), U = .from_rows(U))
}

#' @rdname LUdecomposition_Simple
#' @export
rank_decomposition <- function(M, iszerofunc = NULL, simplify = FALSE) {
  M <- .as_mat(M)
  r <- .matint_rref(M, if (is.null(iszerofunc)) .matint_iszero else iszerofunc, simplify)
  rank <- length(r$pivots)
  C <- M[seq_len(M$nrow), r$pivots, drop = FALSE]
  F <- r$M[seq_len(rank), , drop = FALSE]
  list(C = C, F = F)
}

#' @rdname LUdecomposition_Simple
#' @export
singular_value_decomposition <- function(M) {
  A <- .as_mat(M)
  AH <- .matint_H(A)
  half <- function(X) {
    d <- diagonalize(X)
    ranked <- which(vapply(seq_len(d$D$nrow), function(i) !isTRUE(.mg_iszero(.mget(d$D, i, i))), TRUE))
    V <- d$P[, ranked, drop = FALSE]
    sv <- lapply(ranked, function(i) .sqrt(.mget(d$D, i, i)))
    S <- if (length(sv)) do.call(diag_matrix, sv) else zeros(0L, 0L)
    list(V = .mg_qr(V)$Q, S = S)
  }
  if (A$nrow >= A$ncol) {
    h <- half(.matmul(AH, A))
    U <- .matmul(.matmul(A, h$V), inverse(h$S))
    list(U = U, S = h$S, V = h$V)
  } else {
    h <- half(.matmul(A, AH))
    V <- .matmul(.matmul(AH, h$V), inverse(h$S))
    list(U = h$V, S = h$S, V = V)
  }
}

#' @rdname LUdecomposition_Simple
#' @export
upper_hessenberg_decomposition <- function(M) {
  H <- .as_mat(M)
  if (H$nrow != H$ncol) stop(NonSquareMatrixError("Matrix must be square."))
  n <- H$ncol
  P <- eye(n)
  for (j in seq_len(max(n - 2L, 0L))) {
    rr <- (j + 1L):n
    u <- H[rr, j, drop = FALSE]
    if (u$nrow > 1L && isTRUE(.mg_fuzzy_zero_matrix(u[-1L, , drop = FALSE]))) next
    if (u$nrow == 1L) next
    s <- fcall("sign", list(u$data[[1L]]))
    nu <- .mg_vnorm(u)
    u$data[[1L]] <- if (!is_zero_num(s)) .add(list(u$data[[1L]], .mul(list(s, nu)))) else .add(list(u$data[[1L]], nu))
    nu <- .mg_vnorm(u)
    v <- .mat_map(u, function(e) .div(e, nu))
    vH <- .matint_H(v)
    H[rr, ] <- H[rr, , drop = FALSE] - .mat_map(.matmul(v, .matmul(vH, H[rr, , drop = FALSE])), function(e) .mul(list(.int(2L), e)))
    H[, rr] <- H[, rr, drop = FALSE] - .matmul(.matmul(H[, rr, drop = FALSE], .mat_map(v, function(e) .mul(list(.int(2L), e)))), vH)
    P[, rr] <- P[, rr, drop = FALSE] - .matmul(.matmul(P[, rr, drop = FALSE], .mat_map(v, function(e) .mul(list(.int(2L), e)))), vH)
  }
  list(H = H, P = P)
}

# _householder_vector
.mg_householder <- function(x) {
  if (x$ncol != 1L) .mg_stop("ValueError", "Input must be a column matrix")
  v <- x; vp <- x; vm <- x
  x0 <- .mget(x, 1L, 1L)
  q <- .div(x0, Abs(x0))
  nx <- .mg_vnorm(x)
  vp$data[[1L]] <- .add(list(x0, .mul(list(q, nx))))
  vm$data[[1L]] <- .sub(x0, .mul(list(q, nx)))
  rest <- if (x$nrow > 1L) .mg_vnorm(x[-1L, , drop = FALSE]) else S_Zero
  if (is_zero_num(rest)) {
    bet <- S_Zero
    v$data[[1L]] <- S_One
  } else {
    v <- if (isTRUE(.is_true(Le(.mg_vnorm(vp), .mg_vnorm(vm))))) vp else vm
    v0 <- v$data[[1L]]
    v <- .mat_map(v, function(e) .div(e, v0))
    bet <- .div(.int(2L), .pow(.mg_vnorm(v), .int(2L)))
  }
  list(v = v, bet = bet)
}

.mg_bidiag <- function(M, with_uv) {
  m <- M$nrow; n <- M$ncol
  A <- M
  U <- eye(m); V <- eye(n)
  hh <- function(v, bet, k) eye(k) - .mat_map(.matmul(v, .matint_H(v)), function(e) .mul(list(bet, e)))
  for (i in seq_len(min(m, n))) {
    h <- .mg_householder(A[i:m, i, drop = FALSE])
    H1 <- hh(h$v, h$bet, m - i + 1L)
    A[i:m, i:n] <- .matmul(H1, A[i:m, i:n, drop = FALSE])
    if (with_uv) { tmp <- eye(m); tmp[i:m, i:m] <- H1; U <- .matmul(U, tmp) }
    if (i + 1L <= n - 1L) {
      h <- .mg_householder(t(A[i, (i + 1L):n, drop = FALSE]))
      H2 <- hh(h$v, h$bet, n - i)
      A[i:m, (i + 1L):n] <- .matmul(A[i:m, (i + 1L):n, drop = FALSE], H2)
      if (with_uv) { tmp <- eye(n); tmp[(i + 1L):n, (i + 1L):n] <- H2; V <- .matmul(tmp, V) }
    }
  }
  if (with_uv) list(U, A, V) else A
}

#' @rdname LUdecomposition_Simple
#' @export
bidiagonal_decomposition <- function(M, upper = TRUE) {
  if (!is.logical(upper)) .mg_stop("ValueError", "upper must be a boolean")
  M <- .as_mat(M)
  if (upper) return(.mg_bidiag(M, TRUE))
  X <- .mg_bidiag(.matint_H(M), TRUE)
  list(.matint_H(X[[3L]]), .matint_H(X[[2L]]), .matint_H(X[[1L]]))
}

#' @rdname LUdecomposition_Simple
#' @export
bidiagonalize <- function(M, upper = TRUE) {
  if (!is.logical(upper)) .mg_stop("ValueError", "upper must be a boolean")
  M <- .as_mat(M)
  if (upper) return(.mg_bidiag(M, FALSE))
  .matint_H(.mg_bidiag(.matint_H(M), FALSE))
}

#' Condition number of a matrix
#'
#' `condition_number(M)` is the ratio of the largest to the smallest
#' singular value, as SymPy's `Matrix.condition_number()`.
#'
#' @param M A `symmatrix`.
#' @return An expression.
#' @examples
#' condition_number(Matrix(list(list(1, 0), list(0, 10))))
#' @export
condition_number <- function(M) {
  M <- .as_mat(M)
  if (!length(M$data)) return(S_Zero)
  sv <- singular_values(M)
  .div(do.call(Max, sv), do.call(Min, sv))
}

# ---------------------------------------------------------------------------
# Linear solvers

#' Linear solvers of a matrix
#'
#' The solvers of SymPy's `Matrix`: `diagonal_solve()`,
#' `lower_triangular_solve()`, `upper_triangular_solve()`, `cholesky_solve()`,
#' `LDLsolve()`, `LUsolve()`, `QRsolve()`, `gauss_jordan_solve()` (returns
#' `list(sol, params)`, plus the free variable indices with
#' `freevar = TRUE`), `pinv_solve()`, `cramer_solve()`, `matrix_solve()`
#' (SymPy's `Matrix.solve(rhs, method)`) and `solve_least_squares()`.
#'
#' @param M The coefficient `symmatrix`.
#' @param rhs,B The right hand side.
#' @param iszerofunc The zero test.
#' @param freevar Also return the indices of the free variables (1 based).
#' @param arbitrary_matrix For `pinv_solve()`: the arbitrary part; dummies by default.
#' @param det_method For `cramer_solve()`: `"laplace"`, `"bird"`, another
#'   `det()` method or a function.
#' @param method The method: `"GJ"`, `"LU"`, `"CH"`, `"QR"`, `"LDL"`,
#'   `"PINV"`, `"CRAMER"`, or an `inverse()` method.
#' @return A `symmatrix` (a list for `gauss_jordan_solve()`).
#' @examples
#' A <- Matrix(list(list(2, 1), list(1, 3)))
#' b <- Matrix(list(1, 2))
#' LUsolve(A, b)
#' cholesky_solve(A, b)
#' gauss_jordan_solve(Matrix(list(list(1, 2), list(2, 4))), Matrix(list(3, 6)))
#' @export
diagonal_solve <- function(M, rhs) {
  M <- .as_mat(M); rhs <- .as_mat(rhs)
  if (!isTRUE(is_diagonal(M))) .mg_stop("TypeError", "Matrix should be diagonal")
  if (rhs$nrow != M$nrow) .mg_stop("TypeError", "Size mismatch")
  .mg_new(rhs$nrow, rhs$ncol, function(i, j) .div(.mget(rhs, i, j), .mget(M, i, i)))
}

#' @rdname diagonal_solve
#' @export
lower_triangular_solve <- function(M, rhs) {
  M <- .as_mat(M); rhs <- .as_mat(rhs)
  if (M$nrow != M$ncol) stop(NonSquareMatrixError("Matrix must be square."))
  if (rhs$nrow != M$nrow) .mg_stop("ShapeError", "Matrices size mismatch.")
  if (!is_lower(M)) .mg_stop("ValueError", "Matrix must be lower triangular.")
  X <- zeros(M$nrow, rhs$ncol)
  for (j in seq_len(rhs$ncol)) for (i in seq_len(M$nrow)) {
    if (!.mg_bool(.mget(M, i, i))) .mg_stop("TypeError", "Matrix must be non-singular.")
    s <- .add(c(list(S_Zero), lapply(seq_len(i - 1L), function(k) .mul(list(.mget(M, i, k), .mget(X, k, j))))))
    X$data[[(i - 1L) * X$ncol + j]] <- .div(.sub(.mget(rhs, i, j), s), .mget(M, i, i))
  }
  X
}

#' @rdname diagonal_solve
#' @export
upper_triangular_solve <- function(M, rhs) {
  M <- .as_mat(M); rhs <- .as_mat(rhs)
  if (M$nrow != M$ncol) stop(NonSquareMatrixError("Matrix must be square."))
  if (rhs$nrow != M$nrow) .mg_stop("ShapeError", "Matrix size mismatch.")
  if (!is_upper(M)) .mg_stop("TypeError", "Matrix is not upper triangular.")
  n <- M$nrow
  X <- zeros(n, rhs$ncol)
  for (j in seq_len(rhs$ncol)) for (i in rev(seq_len(n))) {
    if (!.mg_bool(.mget(M, i, i))) .mg_stop("ValueError", "Matrix must be non-singular.")
    ks <- if (i < n) (i + 1L):n else integer()
    s <- .add(c(list(S_Zero), lapply(ks, function(k) .mul(list(.mget(M, i, k), .mget(X, k, j))))))
    X$data[[(i - 1L) * X$ncol + j]] <- .div(.sub(.mget(rhs, i, j), s), .mget(M, i, i))
  }
  X
}

# the reform step shared by _cholesky_solve and _LDLsolve
.mg_ch_prepare <- function(M, rhs) {
  if (M$nrow < M$ncol) .mg_stop("NotImplementedError", "Under-determined System. Try M.gauss_jordan_solve(rhs)")
  hermitian <- TRUE
  reform <- FALSE
  if (isTRUE(is_symmetric(M))) hermitian <- FALSE
  else if (!isTRUE(is_hermitian(M))) reform <- TRUE
  if (reform || isFALSE(.matint_fuzzy_pd(M, FALSE))) {
    H <- .matint_H(M)
    M <- .matmul(H, M)
    rhs <- .matmul(H, rhs)
    hermitian <- !isTRUE(is_symmetric(M))
  }
  list(M = M, rhs = rhs, hermitian = hermitian)
}

#' @rdname diagonal_solve
#' @export
cholesky_solve <- function(M, rhs) {
  p <- .mg_ch_prepare(.as_mat(M), .as_mat(rhs))
  L <- cholesky(p$M, hermitian = p$hermitian)
  Y <- lower_triangular_solve(L, p$rhs)
  upper_triangular_solve(if (p$hermitian) .matint_H(L) else t(L), Y)
}

#' @rdname diagonal_solve
#' @export
LDLsolve <- function(M, rhs) {
  p <- .mg_ch_prepare(.as_mat(M), .as_mat(rhs))
  ld <- LDLdecomposition(p$M, hermitian = p$hermitian)
  Y <- lower_triangular_solve(ld$L, p$rhs)
  Z <- diagonal_solve(ld$D, Y)
  upper_triangular_solve(if (p$hermitian) .matint_H(ld$L) else t(ld$L), Z)
}

#' @rdname diagonal_solve
#' @export
LUsolve <- function(M, rhs, iszerofunc = NULL) {
  # SymPy builds the solution with the class of rhs
  rc <- if (inherits(rhs, "symmatrix")) .mclass(rhs) else "MutableDenseMatrix"
  .set_mclass((function() {
    .matint_lusolve(.as_mat(M), .as_mat(rhs), if (is.null(iszerofunc)) .matint_iszero else iszerofunc)
  })(), rc)
}

#' @rdname diagonal_solve
#' @export
QRsolve <- function(M, B) {
  M <- .as_mat(M); B <- .as_mat(B)
  qr <- .mg_qr(M)
  y <- .matmul(t(qr$Q), B)
  n <- qr$R$nrow
  x <- vector("list", n)
  for (j in rev(seq_len(n))) {
    tmp <- y[j, , drop = FALSE]
    ks <- if (j < n) (j + 1L):n else integer()
    for (k in ks) tmp <- tmp - .mat_map(x[[k]], function(e) .mul(list(.mget(qr$R, j, k), e)))
    tmp <- .mat_map(tmp, expand_mul)
    x[[j]] <- .mat_map(tmp, function(e) .div(e, .mget(qr$R, j, j)))
  }
  do.call(rbind.symmatrix, x)
}

#' @rdname diagonal_solve
#' @export
gauss_jordan_solve <- function(M, B, freevar = FALSE) {
  M <- .as_mat(M); B <- .as_mat(B)
  aug <- cbind.symmatrix(M, B)
  bc <- B$ncol
  row <- M$nrow; col <- M$ncol
  r <- .matint_rref(aug, .matint_iszero, TRUE)
  A <- r$M[, seq_len(col), drop = FALSE]
  v <- r$M[, col + seq_len(bc), drop = FALSE]
  pivots <- r$pivots[r$pivots <= col]
  rank <- length(pivots)
  free <- setdiff(seq_len(col), pivots)
  perm <- c(pivots, free)
  if (rank < row && !isTRUE(.mg_fuzzy_zero_matrix(v[(rank + 1L):row, , drop = FALSE])))
    .mg_stop("ValueError", "Linear system has no solution")
  name <- .mg_unique_tau(aug)
  ntau <- (col - rank) * bc
  taus <- lapply(seq_len(ntau) - 1L, function(k) Symbol(paste0(name, k)))
  tau <- if (ntau) .new_mat(taus, col - rank, bc) else zeros(col - rank, bc)
  vt <- if (rank) v[seq_len(rank), , drop = FALSE] else zeros(0L, bc)
  V <- if (rank && length(free)) A[seq_len(rank), free, drop = FALSE] else zeros(rank, length(free))
  top <- if (rank) (if (length(free)) vt - .matmul(V, tau) else vt) else NULL
  free_sol <- if (is.null(top)) tau else if (tau$nrow) rbind.symmatrix(top, tau) else top
  sol <- zeros(col, bc)
  for (k in seq_len(col)) for (j in seq_len(bc)) sol$data[[(perm[k] - 1L) * bc + j]] <- .mget(free_sol, k, j)
  if (freevar) return(list(sol, tau, free))
  list(sol, tau)
}

# uniquely_named_symbol('tau', [aug], compare = strip trailing digits, modify = '_' + s)
.mg_unique_tau <- function(aug) {
  used <- unique(unlist(lapply(aug$data, function(e) vapply(free_symbols(e), function(s) sub("[0-9]+$", "", s$name), ""))))
  name <- "tau"
  while (name %in% used) name <- paste0("_", name)
  name
}

#' @rdname diagonal_solve
#' @export
pinv_solve <- function(M, B, arbitrary_matrix = NULL) {
  M <- .as_mat(M); B <- .as_mat(B)
  Ap <- pinv(M)
  if (is.null(arbitrary_matrix)) {
    rows <- M$ncol; cols <- B$ncol
    w <- list()
    for (i in seq_len(rows) - 1L) for (j in seq_len(cols) - 1L) w[[length(w) + 1L]] <- Dummy(paste0("w", i, "_", j))
    arbitrary_matrix <- t(.new_mat(w, cols, rows))
  }
  .matmul(Ap, B) + .matmul(eye(M$ncol) - .matmul(Ap, M), Matrix(arbitrary_matrix))
}

#' @rdname diagonal_solve
#' @export
cramer_solve <- function(M, rhs, det_method = "laplace") {
  M <- .as_mat(M); rhs <- .as_mat(rhs)
  dfun <- if (is.function(det_method)) det_method else function(X) .det(X, method = det_method)
  dM <- dfun(M)
  x <- zeros(rhs$nrow, rhs$ncol)
  for (s in seq_len(rhs$ncol)) for (cc in seq_len(rhs$nrow)) {
    X <- .mg_new(M$nrow, M$ncol, function(i, j) if (j == cc) .mget(rhs, i, s) else .mget(M, i, j))
    x$data[[(cc - 1L) * x$ncol + s]] <- .div(dfun(X), dM)
  }
  x
}

#' @rdname diagonal_solve
#' @export
matrix_solve <- function(M, rhs, method = "GJ") {
  M <- .as_mat(M); rhs <- .as_mat(rhs)
  if (method %in% c("GJ", "GE")) {
    r <- tryCatch(gauss_jordan_solve(M, rhs), ValueError = function(e)
      stop(NonInvertibleMatrixError("Matrix det == 0; not invertible.")))
    if (length(r[[2L]]$data))
      stop(NonInvertibleMatrixError("Matrix det == 0; not invertible. Try ``M.gauss_jordan_solve(rhs)`` to obtain a parametric solution."))
    return(r[[1L]])
  }
  switch(method,
    LU = LUsolve(M, rhs), CH = cholesky_solve(M, rhs), QR = QRsolve(M, rhs), LDL = LDLsolve(M, rhs),
    PINV = pinv_solve(M, rhs), CRAMER = cramer_solve(M, rhs),
    .matmul(inverse(M, method = method), rhs))
}

#' @rdname diagonal_solve
#' @export
solve_least_squares <- function(M, rhs, method = "CH") {
  M <- .as_mat(M); rhs <- .as_mat(rhs)
  switch(method,
    CH = cholesky_solve(M, rhs), QR = QRsolve(M, rhs), LDL = LDLsolve(M, rhs), PINV = pinv_solve(M, rhs), {
      tt <- .matint_H(M)
      matrix_solve(.matmul(tt, M), .matmul(tt, rhs), method = method)
    })
}

# _verify_invertible
.mg_verify_invertible <- function(M, iszerofunc = .matint_iszero) {
  if (M$nrow != M$ncol) stop(NonSquareMatrixError("A Matrix must be square to invert."))
  d <- .det(M, method = "berkowitz")
  zero <- .matint_equals0(d)
  if (is.na(zero)) {
    ok <- .matint_rref(M, iszerofunc, TRUE)$M
    zero <- any(vapply(seq_len(ok$nrow), function(j) isTRUE(.matint_call_iszero(iszerofunc, .mget(ok, j, j))), TRUE))
  }
  if (isTRUE(zero)) stop(NonInvertibleMatrixError("Matrix det == 0; not invertible."))
  d
}

.mg_inv_ch <- function(M, iszerofunc) { .mg_verify_invertible(M, iszerofunc); cholesky_solve(M, eye(M$nrow)) }
.mg_inv_ldl <- function(M, iszerofunc) { .mg_verify_invertible(M, iszerofunc); LDLsolve(M, eye(M$nrow)) }
.mg_inv_qr <- function(M, iszerofunc) { .mg_verify_invertible(M, iszerofunc); QRsolve(M, eye(M$nrow)) }

#' Inverse by a named method
#'
#' `inverse_GE()`, `inverse_LU()`, `inverse_ADJ()`, `inverse_CH()`,
#' `inverse_LDL()`, `inverse_QR()` and `inverse_BLOCK()` are
#' `inverse(M, method = )` for each of SymPy's inversion methods.
#'
#' @param M A square `symmatrix`.
#' @param iszerofunc The zero test.
#' @return A `symmatrix`.
#' @examples
#' A <- Matrix(list(list(2, 1), list(1, 3)))
#' inverse_CH(A)
#' inverse_ADJ(A)
#' @export
inverse_GE <- function(M, iszerofunc = NULL) inverse(M, method = "GE", iszerofunc = iszerofunc)

#' @rdname inverse_GE
#' @export
inverse_LU <- function(M, iszerofunc = NULL) inverse(M, method = "LU", iszerofunc = iszerofunc)

#' @rdname inverse_GE
#' @export
inverse_ADJ <- function(M, iszerofunc = NULL) inverse(M, method = "ADJ", iszerofunc = iszerofunc)

#' @rdname inverse_GE
#' @export
inverse_CH <- function(M, iszerofunc = NULL) inverse(M, method = "CH", iszerofunc = iszerofunc)

#' @rdname inverse_GE
#' @export
inverse_LDL <- function(M, iszerofunc = NULL) inverse(M, method = "LDL", iszerofunc = iszerofunc)

#' @rdname inverse_GE
#' @export
inverse_QR <- function(M, iszerofunc = NULL) inverse(M, method = "QR", iszerofunc = iszerofunc)

#' @rdname inverse_GE
#' @export
inverse_BLOCK <- function(M, iszerofunc = NULL) inverse(M, method = "BLOCK", iszerofunc = iszerofunc)

# ---------------------------------------------------------------------------
# Reductions

#' Echelon forms
#'
#' `echelon_form(M)` is the row echelon form of SymPy's
#' `Matrix.echelon_form()` (with `with_pivots = TRUE` also the 1 based pivot
#' columns); `rref_rhs(M, rhs)` reduces `M` augmented with the identity and
#' `rhs` and returns `list(reduced M, reduced rhs)`.
#'
#' @param M A `symmatrix`.
#' @param iszerofunc The zero test.
#' @param simplify A simplification function, or `FALSE` for `simplify()`.
#' @param with_pivots Also return the pivot columns.
#' @param rhs The right hand side.
#' @return A `symmatrix` or a list.
#' @examples
#' echelon_form(Matrix(list(list(1, 2), list(3, 4))))
#' rref_rhs(Matrix(list(list(1, 2), list(3, 4))), Matrix(list(5, 6)))
#' @export
echelon_form <- function(M, iszerofunc = NULL, simplify = FALSE, with_pivots = FALSE) {
  M <- .as_mat(M)
  iz <- if (is.null(iszerofunc)) .matint_iszero else iszerofunc
  sf <- if (is.function(simplify)) simplify else .mg_simplify
  r <- .matint_row_reduce_list(M$data, M$nrow, M$ncol, iz, sf, normalize_last = TRUE, normalize = FALSE, zero_above = FALSE)
  out <- .new_mat(r$mat, M$nrow, M$ncol)
  if (with_pivots) return(list(out, r$pivots + 1L))
  out
}

#' @rdname echelon_form
#' @export
rref_rhs <- function(M, rhs) {
  M <- .as_mat(M); rhs <- .as_mat(rhs)
  r <- .matint_rref(do.call(cbind.symmatrix, list(M, eye(M$nrow), rhs)))$M
  list(r[, seq_len(M$ncol), drop = FALSE], r[, (r$ncol - rhs$ncol + 1L):r$ncol, drop = FALSE])
}

# ---------------------------------------------------------------------------
# Determinant methods

# _det_LU
.mg_det_lu <- function(M, iszerofunc = .matint_iszero, simpfunc = NULL) {
  if (M$nrow != M$ncol) stop(NonSquareMatrixError())
  if (M$nrow == 0L) return(S_One)
  r <- .mg_lu(M, iszerofunc, simpfunc)
  n <- r$lu$nrow
  if (isTRUE(.matint_call_iszero(iszerofunc, .mget(r$lu, n, n)))) return(S_Zero)
  d <- if (length(r$swaps) %% 2L) S_NegativeOne else S_One
  for (k in seq_len(n)) d <- .mul(list(d, .mget(r$lu, k, k)))
  d
}

# _det_bird
.mg_det_bird <- function(M) {
  n <- M$nrow
  if (n == 0L) return(S_One)
  A <- .mrows(M)
  mu <- function(X) {
    total <- S_Zero
    ds <- list(S_Zero)
    for (i in rev(seq_len(n - 1L) + 1L)) {
      total <- .sub(total, X[[i]][[i]])
      ds[[length(ds) + 1L]] <- total
    }
    ds <- rev(ds)
    lapply(seq_len(n), function(i) c(rep(list(S_Zero), i - 1L), list(ds[[i]]), if (i < n) X[[i]][(i + 1L):n] else list()))
  }
  mm <- function(X, Y) lapply(seq_len(n), function(i) lapply(seq_len(n), function(j)
    .add(lapply(seq_len(n), function(k) .mul(list(X[[i]][[k]], Y[[k]][[j]]))))))
  F <- A
  for (s in seq_len(n - 1L)) F <- mm(mu(F), A)
  d <- F[[1L]][[1L]]
  if (n %% 2L == 0L) d <- .neg(d)
  d
}

#' Minor submatrix
#'
#' `minor_submatrix(M, i, j)` deletes row `i` and column `j` (1 based), as
#' SymPy's `Matrix.minor_submatrix()`.
#'
#' @param M A `symmatrix`.
#' @param i,j The row and column to delete.
#' @return A `symmatrix`.
#' @examples
#' minor_submatrix(Matrix(list(list(1, 2, 3), list(4, 5, 6), list(7, 8, 9))), 2, 2)
#' @export
minor_submatrix <- function(M, i, j) {
  M <- .as_mat(M)
  i <- as.integer(i); j <- as.integer(j)
  if (i < 0L) i <- i + M$nrow + 1L
  if (j < 0L) j <- j + M$ncol + 1L
  if (i < 1L || i > M$nrow || j < 1L || j > M$ncol)
    .mg_stop("ValueError", sprintf("`i` and `j` must satisfy 0 <= i < ``M.rows`` (%d)and 0 <= j < ``M.cols`` (%d).", M$nrow, M$ncol))
  M[-i, -j, drop = FALSE]
}

# ---------------------------------------------------------------------------
# Matrix functions by Jordan blocks

#' Matrix logarithm and analytic functions of a matrix
#'
#' `matrix_log(M)` is SymPy's `Matrix.log()`: the logarithm of each Jordan
#' block, transformed back (by default the result is `cancel`led).
#' `analytic_func(M, f, x)` evaluates the analytic function `f(x)` at the
#' matrix `M` by interpolating on the eigenvalues, as `Matrix.analytic_func()`.
#'
#' @param M A square `symmatrix`.
#' @param simplify A function applied to `M` and to the result, or `FALSE`.
#' @param f An expression in `x`.
#' @param x The symbol of `f`.
#' @return A `symmatrix`.
#' @examples
#' matrix_log(Matrix(list(list(1, 1), list(0, 1))))
#' x <- Symbol("x")
#' analytic_func(Matrix(list(list(2, 1), list(0, 2))), exp(x), x)
#' @export
matrix_log <- function(M, simplify = cancel) {
  M <- .as_mat(M)
  if (M$nrow != M$ncol) stop(NonSquareMatrixError("Logarithm is valid only for square matrices"))
  sf <- if (is.function(simplify)) function(X) .mat_map(X, simplify) else NULL
  jf <- tryCatch(jordan_form(if (is.null(sf)) M else sf(M)), MatrixError = function(e)
    .mg_stop("NotImplementedError", "Logarithm is implemented only for matrices for which the Jordan normal form can be computed"))
  cells <- .matint_get_diag_blocks(jf$J)
  blocks <- lapply(cells, function(cell) {
    size <- cell$nrow
    l <- .mget(cell, 1L, 1L)
    if (isTRUE(.mg_iszero(l)))
      .mg_stop("MatrixError", sprintf("Could not take logarithm or reciprocal for the given eigenvalue %s", sstr(l)))
    .mg_new(size, size, function(i, j) {
      d <- j - i
      if (d < 0L) S_Zero else if (d == 0L) log(l) else .neg(.div(.pow(.neg(l), .int(-d)), .int(d)))
    })
  })
  eJ <- do.call(diag_matrix, blocks)
  if (!is.null(sf)) return(sf(.matmul(.matmul(jf$P, eJ), sf(inverse(jf$P)))))
  .matmul(.matmul(jf$P, eJ), inverse(jf$P))
}

#' @rdname matrix_log
#' @export
analytic_func <- function(M, f, x) {
  M <- .as_mat(M)
  f <- as_sym(f); x <- as_sym(x)
  if (M$nrow != M$ncol) stop(NonSquareMatrixError())
  if (!(x$type %in% c("Symbol", "Dummy"))) .mg_stop("ValueError", sprintf("%s must be a symbol.", sstr(x)))
  if (!has(f, x)) .mg_stop("ValueError", sprintf("%s must be a parameter of %s.", sstr(x), sstr(f)))
  if (any(vapply(M$data, function(e) has(e, x), TRUE)))
    .mg_stop("ValueError", sprintf("%s must not be a parameter of %s.", sstr(x), sstr(M)))
  ev <- .matint_eigenvals_dict(M)
  vals <- unclass(ev)
  attr(vals, "multiplicity") <- NULL
  muls <- as.list(attr(ev, "multiplicity"))
  max_mul <- max(unlist(muls))
  derivative <- list()
  dd <- f
  for (ii in seq_len(max_mul - 1L)) {
    dd <- diff(dd, x)
    derivative[[ii]] <- dd
  }
  n <- M$nrow
  r <- rep(list(rep(list(S_Zero), n)), n)
  f_val <- rep(list(S_Zero), n)
  row <- 1L
  for (k in seq_along(vals)) {
    i <- vals[[k]]
    mul <- muls[[k]]
    f_val[[row]] <- subs(f, x, i)
    if (.mg_not_complex(f_val[[row]]))
      .mg_stop("ValueError", sprintf("Cannot evaluate the function because the function %s is not analytic at the given eigenvalue %s", sstr(f), sstr(f_val[[row]])))
    val <- S_One
    for (a in seq_len(n)) { r[[row]][[a]] <- val; val <- .mul(list(val, i)) }
    if (mul > 1L) {
      coe <- rep(1L, n)
      deri <- 1L
      while (mul > 1L) {
        row <- row + 1L
        mul <- mul - 1L
        d_i <- subs(derivative[[deri]], x, i)
        if (.mg_not_complex(d_i))
          .mg_stop("ValueError", sprintf("Cannot evaluate the function because the derivative %s is not analytic at the given eigenvalue %s", sstr(derivative[[deri]]), sstr(d_i)))
        f_val[[row]] <- d_i
        for (a in seq_len(n)) {
          ap <- a - 1L
          if (ap - deri + 1L <= 0L) { r[[row]][[a]] <- S_Zero; coe[a] <- 0L; next }
          coe[a] <- coe[a] * (ap - deri + 1L)
          r[[row]][[a]] <- .mul(list(.int(coe[a]), .pow(i, .int(ap - deri))))
        }
        deri <- deri + 1L
      }
    }
    row <- row + 1L
  }
  cc <- matrix_solve(.from_rows(r), .new_mat(f_val, n, 1L))
  ans <- zeros(n)
  pre <- eye(n)
  for (i2 in seq_len(n)) {
    ci <- cc$data[[i2]]
    ans <- ans + .mat_map(pre, function(e) .mul(list(ci, e)))
    pre <- .matmul(pre, M)
  }
  ans
}

.mg_not_complex <- function(e) isTRUE(.is_number(e)) && isFALSE(.matint_fact(e, "complex"))

# ---------------------------------------------------------------------------
# Graph decompositions

.mg_graph <- function(M) {
  E <- lapply(.mg_items(M), function(it) list(it$i, it$j))
  list(as.list(seq_len(M$nrow)), E)
}

.mg_components <- function(M, strong) {
  M <- .as_mat(M)
  g <- .mg_graph(M)
  cc <- if (strong) strongly_connected_components(g) else connected_components(g)
  lapply(cc, function(b) as.integer(unlist(b)))
}

#' Block decompositions from the graph of a matrix
#'
#' `connected_components_decomposition(M)` returns `list(P, B)` with a
#' `PermutationMatrix` `P` and a `BlockDiagMatrix` `B` such that
#' `M = P^T B P`; `strongly_connected_components_decomposition(M, lower)`
#' does the same with a block triangular `BlockMatrix`, as SymPy.
#'
#' @param M A square `symmatrix`.
#' @param lower Make the blocks lower triangular (`TRUE`) or upper.
#' @return A list of two matrix expressions.
#' @examples
#' A <- Matrix(list(list(1, 0, 2), list(0, 3, 0), list(4, 0, 5)))
#' connected_components_decomposition(A)
#' @export
connected_components_decomposition <- function(M) {
  M <- .as_mat(M)
  ib <- .mg_components(M, FALSE)
  P <- PermutationMatrix(Permutation(unlist(ib) - 1L))
  B <- do.call(BlockDiagMatrix, lapply(ib, function(b) M[b, b, drop = FALSE]))
  list(P, B)
}

#' @rdname connected_components_decomposition
#' @export
strongly_connected_components_decomposition <- function(M, lower = TRUE) {
  M <- .as_mat(M)
  ib <- .mg_components(M, TRUE)
  if (!lower) ib <- rev(ib)
  P <- PermutationMatrix(Permutation(unlist(ib) - 1L))
  B <- BlockMatrix(lapply(ib, function(a) lapply(ib, function(b) M[a, b, drop = FALSE])))
  list(P, B)
}

#' Normalized vector and nonzero pattern
#'
#' `normalized_vector(v)` divides a row or column vector by its norm, as
#' SymPy's `Matrix.normalized()`; `print_nonzero(M, symb)` prints the
#' pattern of nonzero entries.
#'
#' @param v A vector `symmatrix`.
#' @param iszerofunc The zero test for the norm.
#' @param M A `symmatrix`.
#' @param symb The character marking a nonzero entry.
#' @return A `symmatrix`; `print_nonzero()` returns `M` invisibly.
#' @examples
#' normalized_vector(Matrix(list(3, 4)))
#' print_nonzero(eye(3))
#' @export
normalized_vector <- function(v, iszerofunc = NULL) {
  v <- Matrix(v)
  if (v$nrow != 1L && v$ncol != 1L) .mg_stop("ShapeError", "A Matrix must be a vector to normalize.")
  nrm <- .mg_vnorm(v)
  iz <- if (is.null(iszerofunc)) .matint_iszero else iszerofunc
  if (isTRUE(iz(nrm))) return(zeros(v$nrow, v$ncol))
  .mat_map(v, function(e) .div(e, nrm))
}

#' @rdname normalized_vector
#' @export
print_nonzero <- function(M, symb = "X") {
  M <- .as_mat(M)
  for (i in seq_len(M$nrow))
    cat("[", paste(vapply(seq_len(M$ncol), function(j) if (is_zero_num(.mget(M, i, j))) " " else symb, ""), collapse = ""), "]\n", sep = "")
  invisible(M)
}
