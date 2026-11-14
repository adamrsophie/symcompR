# Definiteness tests after sympy/matrices/eigen.py (_is_positive_definite,
# _is_positive_semidefinite, _is_negative_definite, _is_negative_semidefinite,
# _is_indefinite and their helpers) and the diagonal dominance properties of
# sympy/matrices/matrixbase.py.

# fuzzy_and over TRUE / FALSE / NA
.matint_fuzzy_and <- function(v) {
  v <- unlist(v)
  if (any(v %in% FALSE)) FALSE else if (anyNA(v)) NA else TRUE
}

# fuzzy_or over TRUE / FALSE / NA
.matint_fuzzy_or <- function(v) {
  v <- unlist(v)
  if (any(v %in% TRUE)) TRUE else if (anyNA(v)) NA else FALSE
}

.matint_fact <- function(e, fact) {
  r <- ask_fact(e, fact)
  if (is.null(r) || length(r) != 1L) NA else as.logical(r)
}

# M.H
.matint_H <- function(M) transpose(.mat_map(M, conjugate))

# _eval_is_matrix_hermitian(simplify): each entry minus the adjoint of its
# mirror simplifies to zero
.matint_hermitian <- function(M) {
  if (M$nrow != M$ncol) return(FALSE)
  .matint_fuzzy_and(lapply(seq_len(M$nrow), function(i) lapply(seq_len(M$ncol), function(j)
    .matint_fact(simplify(.sub(.mget(M, i, j), conjugate(.mget(M, j, i)))), "zero"))))
}

.matint_diag_fact <- function(M, fact) .matint_fuzzy_and(lapply(seq_len(M$nrow), function(i) .matint_fact(.mget(M, i, i), fact)))

# is_strongly_diagonally_dominant (fact "positive") and
# is_weakly_diagonally_dominant (fact "nonnegative")
.matint_diag_dominant <- function(M, fact) {
  if (M$nrow != M$ncol) return(FALSE)
  .matint_fuzzy_and(lapply(seq_len(M$nrow), function(i) {
    s <- .add(lapply(setdiff(seq_len(M$ncol), i), function(j) Abs(.mget(M, i, j))))
    .matint_fact(.sub(Abs(.mget(M, i, i)), s), fact)
  }))
}

# _fuzzy_positive_definite and _fuzzy_positive_semidefinite
.matint_fuzzy_pd <- function(M, semi) {
  d <- .matint_diag_fact(M, if (semi) "nonnegative" else "positive")
  if (isFALSE(d)) return(FALSE)
  if (isTRUE(d) && isTRUE(.matint_diag_dominant(M, if (semi) "nonnegative" else "positive"))) return(TRUE)
  NA
}

# _is_positive_definite_GE
.matint_pd_ge <- function(M) {
  n <- M$nrow
  rows <- lapply(seq_len(n), function(i) lapply(seq_len(n), function(j) .mget(M, i, j)))
  for (i in seq_len(n)) {
    p <- .matint_fact(rows[[i]][[i]], "positive")
    if (!isTRUE(p)) return(p)
    if (i < n) for (j in (i + 1L):n) for (k in (i + 1L):n)
      rows[[j]][[k]] <- .sub(.mul(list(rows[[i]][[i]], rows[[j]][[k]])), .mul(list(rows[[j]][[i]], rows[[i]][[k]])))
  }
  TRUE
}

# _is_positive_semidefinite_cholesky
.matint_psd_cholesky <- function(M) {
  n <- M$nrow
  A <- lapply(seq_len(n), function(i) lapply(seq_len(n), function(j) .mget(M, i, j)))
  for (k in seq_len(n)) {
    diags <- lapply(k:n, function(i) A[[i]][[i]])
    fp <- .matint_find_pivot(diags, .matint_iszero, .simplify_default)
    if (isTRUE(fp$assumed)) return(NA)
    if (is.na(fp$offset)) {
      if (k < n) for (i in (k + 1L):n) for (j in k:n) {
        iz <- .matint_fact(A[[i]][[j]], "zero")
        if (is.na(iz)) return(NA)
        if (!iz) return(FALSE)
      }
      return(TRUE)
    }
    if (isTRUE(.matint_fact(A[[k]][[k]], "negative")) || isTRUE(.matint_fact(fp$value, "negative"))) return(FALSE)
    if (!(isTRUE(.matint_fact(A[[k]][[k]], "nonnegative")) && isTRUE(.matint_fact(fp$value, "nonnegative")))) return(NA)
    if (fp$offset > 0L) {
      q <- k + fp$offset
      for (i in seq_len(n)) { tmp <- A[[i]][[k]]; A[[i]][[k]] <- A[[i]][[q]]; A[[i]][[q]] <- tmp }
      tmp <- A[[k]]; A[[k]] <- A[[q]]; A[[q]] <- tmp
    }
    A[[k]][[k]] <- sqrt(A[[k]][[k]])
    if (k < n) {
      for (j in (k + 1L):n) A[[k]][[j]] <- .div(A[[k]][[j]], A[[k]][[k]])
      for (i in (k + 1L):n) for (j in (k + 1L):n)
        A[[i]][[j]] <- .sub(A[[i]][[j]], .mul(list(conjugate(A[[k]][[i]]), A[[k]][[j]])))
    }
  }
  .matint_fact(A[[n]][[n]], "nonnegative")
}

.matint_is_pd <- function(M, semi) {
  M <- Matrix(M)
  if (!isTRUE(.matint_hermitian(M))) {
    if (M$nrow != M$ncol) return(FALSE)
    M <- M + .matint_H(M)
  }
  fz <- .matint_fuzzy_pd(M, semi)
  if (!is.na(fz)) return(fz)
  if (semi) .matint_psd_cholesky(M) else .matint_pd_ge(M)
}

.matint_neg <- function(M) .mat_map(Matrix(M), .neg)

#' Definiteness of a matrix
#'
#' Fuzzy tests of whether a matrix is positive or negative (semi)definite, or
#' indefinite. A matrix that is not Hermitian is tested through `M + M^H`.
#' Each returns `TRUE`, `FALSE` or `NA` when the answer cannot be decided.
#'
#' @param M A matrix.
#' @return `TRUE`, `FALSE` or `NA`.
#' @examples
#' is_positive_definite(Matrix(list(list(2, -1), list(-1, 2))))
#' is_indefinite(Matrix(list(list(1, 2), list(2, 1))))
#' @export
is_positive_definite <- function(M) .matint_is_pd(M, FALSE)

#' @rdname is_positive_definite
#' @export
is_positive_semidefinite <- function(M) .matint_is_pd(M, TRUE)

#' @rdname is_positive_definite
#' @export
is_negative_definite <- function(M) .matint_is_pd(.matint_neg(M), FALSE)

#' @rdname is_positive_definite
#' @export
is_negative_semidefinite <- function(M) .matint_is_pd(.matint_neg(M), TRUE)

#' @rdname is_positive_definite
#' @export
is_indefinite <- function(M) {
  M <- Matrix(M)
  if (isTRUE(.matint_hermitian(M))) {
    ev <- eigenvals(M)
    pos <- .matint_fuzzy_or(lapply(ev, .matint_fact, fact = "positive"))
    neg <- .matint_fuzzy_or(lapply(ev, .matint_fact, fact = "negative"))
    return(.matint_fuzzy_and(list(pos, neg)))
  }
  if (M$nrow == M$ncol) return(is_indefinite(M + .matint_H(M)))
  FALSE
}
