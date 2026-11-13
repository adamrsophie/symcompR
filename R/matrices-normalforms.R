# Smith and Hermite normal forms over ZZ, after
# sympy/polys/matrices/normalforms.py and sympy/matrices/normalforms.py.
# Integer matrices are lists of rows, each a bigz vector.

# igcdex (sympy/core/intfunc.py): x*a + y*b = g with g = gcd(a, b) >= 0
.igcdex <- function(a, b) {
  a <- gmp::as.bigz(a); b <- gmp::as.bigz(b)
  if (a == 0 && b == 0) return(list(gmp::as.bigz(0), gmp::as.bigz(1), gmp::as.bigz(0)))
  if (a == 0) return(list(gmp::as.bigz(0), gmp::divq.bigz(b, abs(b)), abs(b)))
  if (b == 0) return(list(gmp::divq.bigz(a, abs(a)), gmp::as.bigz(0), abs(a)))
  xs <- if (a < 0) { a <- -a; -1L } else 1L
  ys <- if (b < 0) { b <- -b; -1L } else 1L
  x <- gmp::as.bigz(1); y <- gmp::as.bigz(0); r <- gmp::as.bigz(0); s <- gmp::as.bigz(1)
  while (b != 0) {
    c <- gmp::mod.bigz(a, b); q <- gmp::divq.bigz(a, b)
    nr <- x - q * r; ns <- y - q * s
    a <- b; b <- c; x <- r; y <- s; r <- nr; s <- ns
  }
  list(x * xs, y * ys, a)
}
# _gcdex of normalforms.py
.nf_gcdex <- function(a, b) {
  r <- .igcdex(a, b)
  if (a != 0 && gmp::mod.bigz(b, a) == 0) { r[[2L]] <- gmp::as.bigz(0); r[[1L]] <- gmp::as.bigz(if (a < 0) -1L else 1L) }
  r
}
# Python // on integers (gmp's divq is floor division)
.zm_floordiv <- function(a, b) gmp::divq.bigz(a, b)
.zm_add_columns <- function(m, i, j, a, b, c, d) {
  for (k in seq_along(m)) {
    e <- m[[k]][i]
    m[[k]][i] <- a * e + b * m[[k]][j]
    m[[k]][j] <- c * e + d * m[[k]][j]
  }
  m
}
.zm_from_sym <- function(M) {
  nr <- M$nrow; nc <- M$ncol
  lapply(seq_len(nr), function(i) gmp::as.bigz(vapply(seq_len(nc), function(j) {
    v <- .mget(M, i, j)
    if (!(v$type == "Integer")) stop("Matrix must be over domain ZZ.", call. = FALSE)
    as.character(v$value)
  }, "")))
}
.zm_to_sym <- function(m, nc = NULL) {
  nr <- length(m)
  nc <- nc %||% (if (nr) length(m[[1L]]) else 0L)
  Matrix(lapply(seq_len(nr), function(i) lapply(seq_len(nc), function(j) .int(as.character(m[[i]][j])))))
}

# _hermite_normal_form
.hnf <- function(A, ncols) {
  m <- length(A); n <- ncols
  k <- n
  for (i in rev(seq_len(m))) {
    if (k == 0L) break
    k <- k - 1L
    if (k >= 1L) for (j in rev(seq_len(k))) {
      if (A[[i]][j] != 0) {
        g <- .nf_gcdex(A[[i]][k + 1L], A[[i]][j])
        r <- .zm_floordiv(A[[i]][k + 1L], g[[3L]]); s <- .zm_floordiv(A[[i]][j], g[[3L]])
        A <- .zm_add_columns(A, k + 1L, j, g[[1L]], g[[2L]], -s, r)
      }
    }
    b <- A[[i]][k + 1L]
    if (b < 0) { A <- .zm_add_columns(A, k + 1L, k + 1L, -1, 0, -1, 0); b <- -b }
    if (b == 0) k <- k + 1L
    else if (k + 2L <= n) for (j in (k + 2L):n) {
      q <- .zm_floordiv(A[[i]][j], b)
      A <- .zm_add_columns(A, j, k + 1L, 1, -q, 0, 1)
    }
  }
  lapply(A, function(row) if (k < n) row[(k + 1L):n] else gmp::as.bigz(integer(0)))
}
.symmetric_residue <- function(a, m) { if (a <= gmp::divq.bigz(m, 2)) a else a - m }
# _hermite_normal_form_modulo_D
.hnf_mod_D <- function(A, ncols, D) {
  m <- length(A); n <- ncols
  if (n < m) stop("Matrix must have at least as many columns as rows.", call. = FALSE)
  addc <- function(M, R, i, j, a, b, c, d) {
    for (kk in seq_along(M)) {
      e <- M[[kk]][i]
      M[[kk]][i] <- .symmetric_residue(gmp::mod.bigz(a * e + b * M[[kk]][j], R), R)
      M[[kk]][j] <- .symmetric_residue(gmp::mod.bigz(c * e + d * M[[kk]][j], R), R)
    }
    M
  }
  W <- lapply(seq_len(m), function(i) gmp::as.bigz(rep(0L, m)))
  k <- n; R <- gmp::as.bigz(D)
  for (i in rev(seq_len(m))) {
    k <- k - 1L
    if (k >= 1L) for (j in rev(seq_len(k))) {
      if (A[[i]][j] != 0) {
        g <- .nf_gcdex(A[[i]][k + 1L], A[[i]][j])
        r <- .zm_floordiv(A[[i]][k + 1L], g[[3L]]); s <- .zm_floordiv(A[[i]][j], g[[3L]])
        A <- addc(A, R, k + 1L, j, g[[1L]], g[[2L]], -s, r)
      }
    }
    b <- A[[i]][k + 1L]
    if (b == 0) { A[[i]][k + 1L] <- R; b <- R }
    g <- .nf_gcdex(b, R)
    for (ii in seq_len(m)) W[[ii]][i] <- gmp::mod.bigz(g[[1L]] * A[[ii]][k + 1L], R)
    if (W[[i]][i] == 0) W[[i]][i] <- R
    if (i + 1L <= m) for (j in (i + 1L):m) {
      q <- .zm_floordiv(W[[i]][j], W[[i]][i])
      W <- .zm_add_columns(W, j, i, 1, -q, 0, 1)
    }
    R <- gmp::divq.bigz(R, g[[3L]])
  }
  W
}

#' Hermite and Smith normal forms
#'
#' `hermite_normal_form(A, D)` is the Hermite normal form of an integer
#' matrix (SymPy's `hermite_normal_form`, optionally modulo `D`);
#' `smith_normal_form(A)` and `invariant_factors(A)` give the Smith normal
#' form and its diagonal.
#'
#' @param A An integer matrix.
#' @param D Optional positive integer modulus (a multiple of the determinant).
#' @param check_rank With `D`, use the modulo method only for full row rank.
#' @return A matrix, or for `invariant_factors()` a Tuple of integers.
#' @examples
#' A <- Matrix(list(list(2, 7, 17, 29, 41), list(3, 11, 19, 31, 43), list(5, 13, 23, 37, 47)))
#' hermite_normal_form(A)
#' smith_normal_form(Matrix(list(list(12, 6, 4), list(3, 9, 6), list(2, 16, 14))))
#' @export
hermite_normal_form <- function(A, D = NULL, check_rank = FALSE) {
  if (inherits(A, "DomainMatrix")) return(.dm_hnf(A, D, check_rank))
  m <- .zm_from_sym(A)
  if (!is.null(D) && (!check_rank || rank(A) == A$nrow)) return(.zm_to_sym(.hnf_mod_D(m, A$ncol, D)))
  h <- .hnf(m, A$ncol)
  .zm_to_sym(h, if (length(h)) length(h[[1L]]) else 0L)
}

# _smith_normal_decomp(full = FALSE)
.snf_invariants <- function(m, rows, cols) {
  if (rows == 0L || cols == 0L) return(list())
  zero <- gmp::as.bigz(0)
  add_rows <- function(M, i, j, a, b, c, d) {
    e <- M[[i]]; M[[i]] <- a * e + b * M[[j]]; M[[j]] <- c * e + d * M[[j]]; M
  }
  clear_column <- function(M) {
    pivot <- M[[1L]][1L]
    if (rows >= 2L) for (j in 2:rows) {
      if (M[[j]][1L] == 0) next
      d <- .zm_floordiv(M[[j]][1L], pivot); r <- M[[j]][1L] - d * pivot
      if (r == 0) M <- add_rows(M, 1L, j, 1, 0, -d, 1)
      else {
        g <- .igcdex(pivot, M[[j]][1L])
        d0 <- gmp::divq.bigz(M[[j]][1L], g[[3L]]); dj <- gmp::divq.bigz(pivot, g[[3L]])
        M <- add_rows(M, 1L, j, g[[1L]], g[[2L]], d0, -dj)
        pivot <- g[[3L]]
      }
    }
    M
  }
  clear_row <- function(M) {
    pivot <- M[[1L]][1L]
    if (cols >= 2L) for (j in 2:cols) {
      if (M[[1L]][j] == 0) next
      d <- .zm_floordiv(M[[1L]][j], pivot); r <- M[[1L]][j] - d * pivot
      if (r == 0) M <- .zm_add_columns(M, 1L, j, 1, 0, -d, 1)
      else {
        g <- .igcdex(pivot, M[[1L]][j])
        d0 <- gmp::divq.bigz(M[[1L]][j], g[[3L]]); dj <- gmp::divq.bigz(pivot, g[[3L]])
        M <- .zm_add_columns(M, 1L, j, g[[1L]], g[[2L]], d0, -dj)
        pivot <- g[[3L]]
      }
    }
    M
  }
  ind <- which(vapply(m, function(r) r[1L] != 0, NA))
  if (length(ind) && ind[1L] != 1L) { t <- m[[1L]]; m[[1L]] <- m[[ind[1L]]]; m[[ind[1L]]] <- t }
  else {
    ind <- which(as.logical(m[[1L]] != 0))
    if (length(ind) && ind[1L] != 1L) m <- lapply(m, function(r) { t <- r[1L]; r[1L] <- r[ind[1L]]; r[ind[1L]] <- t; r })
  }
  more <- function(M) (cols >= 2L && any(as.logical(M[[1L]][2:cols] != 0))) || (rows >= 2L && any(vapply(M[-1L], function(r) r[1L] != 0, NA)))
  while (more(m)) { m <- clear_column(m); m <- clear_row(m) }
  if (m[[1L]][1L] != 0 && m[[1L]][1L] < 0) m[[1L]][1L] <- -m[[1L]][1L]
  invs <- if (rows == 1L || cols == 1L) list() else .snf_invariants(lapply(m[-1L], function(r) r[-1L]), rows - 1L, cols - 1L)
  if (m[[1L]][1L] != 0) {
    result <- c(list(m[[1L]][1L]), invs)
    for (i in seq_len(length(result) - 1L)) {
      a <- result[[i]]; b <- result[[i + 1L]]
      if (b != 0 && gmp::mod.bigz(b, a) != 0) {
        d <- gmp::gcd.bigz(a, b)
        alpha <- .zm_floordiv(a, d)
        result[[i + 1L]] <- b * alpha
        result[[i]] <- d
      } else break
    }
  } else result <- c(invs, list(m[[1L]][1L]))
  result
}

#' @rdname hermite_normal_form
#' @export
invariant_factors <- function(A) {
  if (inherits(A, "DomainMatrix")) return(.dm_invariant_factors(A))
  r <- .snf_invariants(.zm_from_sym(A), A$nrow, A$ncol)
  do.call(Tuple, lapply(r, function(v) .int(as.character(v))))
}

#' @rdname hermite_normal_form
#' @export
smith_normal_form <- function(A) {
  if (inherits(A, "DomainMatrix")) return(.dm_smith_normal_form(A))
  inv <- invariant_factors(A)$args
  nr <- A$nrow; nc <- A$ncol
  Matrix(lapply(seq_len(nr), function(i) lapply(seq_len(nc), function(j) if (i == j && i <= length(inv)) inv[[i]] else S_Zero)))
}
