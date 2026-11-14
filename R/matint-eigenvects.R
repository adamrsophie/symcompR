# Eigenvectors after sympy/matrices/eigen.py (_eigenvects, _eigenvects_DOM,
# _eigenvects_sympy, _eigenspace, _left_eigenvects, _singular_values,
# _is_diagonalizable and _diagonalize) and polys/matrices/eigen.py
# (dom_eigenvects and dom_eigenvects_to_sympy).

# Gauss-Jordan over a field given by its zero test, inverse and reduction;
# ddm_irref takes the first nonzero entry of each column as the pivot
.matint_field_rref <- function(rows, ncol, red, inv) {
  nr <- length(rows)
  piv <- integer()
  i <- 1L
  for (j in seq_len(ncol)) {
    if (i > nr) break
    k <- NA_integer_
    for (r in i:nr) if (!is_zero_num(rows[[r]][[j]])) { k <- r; break }
    if (is.na(k)) next
    if (k != i) { tmp <- rows[[i]]; rows[[i]] <- rows[[k]]; rows[[k]] <- tmp }
    pinv <- inv(rows[[i]][[j]])
    rows[[i]] <- lapply(rows[[i]], function(z) red(.mul(list(z, pinv))))
    for (r in seq_len(nr)) {
      if (r == i || is_zero_num(rows[[r]][[j]])) next
      fct <- rows[[r]][[j]]
      rows[[r]] <- lapply(seq_len(ncol), function(c) red(.sub(rows[[r]][[c]], .mul(list(fct, rows[[i]][[c]])))))
    }
    piv <- c(piv, j)
    i <- i + 1L
  }
  list(rows = rows, pivots = piv)
}

# DomainMatrix.nullspace(divide_last=True) from a field rref
.matint_nullspace_rows <- function(rr, ncol) {
  lapply(setdiff(seq_len(ncol), rr$pivots), function(fv) {
    v <- rep(list(S_Zero), ncol)
    v[[fv]] <- S_One
    for (k in seq_along(rr$pivots)) v[[rr$pivots[k]]] <- .neg(rr$rows[[k]][[fv]])
    v
  })
}

# The domain of DomainMatrix.from_Matrix(M, field=True) when it is QQ or a
# fraction field over it; NULL for EX and for algebraic number fields
.matint_dom_field <- function(M) {
  if (.matint_is_qq(M)) return("QQ")
  K <- tryCatch(construct_domain(M$data, field = TRUE)[[1L]], error = function(e) NULL)
  if (is.null(K) || !identical(.de_kind(K), "FracField")) return(NULL)
  "frac"
}

# FiniteExtension(minpoly) element in canonical form: the remainder modulo
# the minimal polynomial with each coefficient a reduced fraction
.matint_fe_reduce <- function(z, base, l) {
  p <- rem(expand(z), base, l)
  if (is_zero_num(p)) return(S_Zero)
  cf <- .coeff_list(p, l)
  k <- length(cf)
  .add(lapply(seq_len(k), function(i) .mul(list(cancel(cf[[i]]), .pow(l, .int(k - i))))))
}

# dom_eigenvects and dom_eigenvects_to_sympy for a matrix over QQ or a
# fraction field over QQ; an irreducible factor of degree two or more is
# handled in K[l]/(factor)
.matint_eigenvects_dom <- function(M, flags = list(), frac = FALSE) {
  n <- M$nrow
  l <- Dummy("lambda")
  rows0 <- lapply(seq_len(n), function(i) lapply(seq_len(n), function(j) .mget(M, i, j)))
  cp <- charpoly(M, l)
  res <- list()
  canon <- if (frac) cancel else identity
  for (fm in factor_list(cp, l)$factors) {
    base <- fm[[1L]]; mult <- as.integer(fm[[2L]])
    deg <- as.integer(degree(base, l))
    if (deg == 0L) next
    if (deg == 1L) {
      cf <- .coeff_list(base, l)
      ev <- canon(.neg(.div(cf[[2L]], cf[[1L]])))
      rows <- lapply(seq_len(n), function(i) {
        r <- rows0[[i]]; r[[i]] <- canon(.sub(r[[i]], ev)); r
      })
      rr <- .matint_field_rref(rows, n, canon, function(z) canon(.pow(z, S_NegativeOne)))
      vecs <- lapply(.matint_nullspace_rows(rr, n), function(v) .new_mat(v, n, 1L))
      res[[length(res) + 1L]] <- list(value = ev, multiplicity = mult, vectors = vecs)
    } else {
      red <- function(z) .matint_fe_reduce(z, base, l)
      inv <- function(z) red(invert(z, base, l))
      rows <- lapply(seq_len(n), function(i) {
        r <- rows0[[i]]; r[[i]] <- .sub(r[[i]], l); r
      })
      rr <- .matint_field_rref(rows, n, red, inv)
      ns <- .matint_nullspace_rows(rr, n)
      evs <- tryCatch(do.call(roots, c(list(base, l), flags)), error = function(e) list())
      if (length(evs) != deg) evs <- .all_roots_rootof(base, l)
      for (ev in evs) {
        vecs <- lapply(ns, function(v) .new_mat(lapply(v, function(z) subs(z, l, ev)), n, 1L))
        res[[length(res) + 1L]] <- list(value = ev, multiplicity = mult, vectors = vecs)
      }
    }
  }
  ord <- unlist(.key_sort(lapply(res, function(e) .sort_key(e$value)), as.list(seq_along(res))))
  res[ord]
}

# _eigenspace
.matint_eigenspace <- function(M, val, iszerofunc, simplify) {
  m <- M - .mat_map(eye(M$nrow), function(z) .mul(list(z, val)))
  ret <- nullspace(m, iszerofunc = iszerofunc)
  if (length(ret) == 0L && !isFALSE(simplify)) ret <- nullspace(m, simplify = TRUE, iszerofunc = iszerofunc)
  if (length(ret) == 0L) stop(.sympy_condition("NotImplementedError", paste0("Can't evaluate eigenvector for eigenvalue ", sstr(val))))
  ret
}

# _eigenvects_sympy
.matint_eigenvects_sympy <- function(M, iszerofunc, simplify, flags) {
  ev <- do.call(.matint_eigenvals, list(M, TRUE, FALSE, FALSE, FALSE, flags))
  for (v in ev) if (.has_type(v, "CRootOf"))
    stop(MatrixError("Eigenvector computation is not implemented if the matrix have eigenvalues in CRootOf form"))
  ord <- unlist(.key_sort(lapply(ev, .sort_key), as.list(seq_along(ev))))
  lapply(ord, function(k) list(value = ev[[k]], multiplicity = attr(ev, "multiplicity")[k],
                               vectors = .matint_eigenspace(M, ev[[k]], iszerofunc, simplify)))
}

# _eigenvects
.matint_eigenvects <- function(M, error_when_incomplete = TRUE, iszerofunc = NULL, simplify = FALSE,
                               chop = FALSE, flags = list()) {
  M <- Matrix(M)
  if (is.null(iszerofunc)) iszerofunc <- .matint_iszero
  primitive <- !isFALSE(simplify)
  simpfunc <- if (is.function(simplify)) simplify else if (isTRUE(simplify)) .simplify_default else identity
  has_floats <- any(vapply(M$data, function(z) .has_type(z, "Float"), TRUE))
  # ponytail: matrices of Float numbers go through nsimplify here; SymPy uses
  # _eigenvects_mpmath for them, port it when numeric eigenvectors matter
  if (has_floats) M <- .mat_map(M, function(z) nsimplify(z, rational = TRUE))
  dom <- if (M$nrow == M$ncol) .matint_dom_field(M)
  # ponytail: algebraic number fields (QQ<sqrt(3), I>) take the eigenvals
  # path; port DomainMatrix over AlgebraicField to follow SymPy there too
  ret <- if (!is.null(dom)) .matint_eigenvects_dom(M, flags, frac = dom == "frac")
         else .matint_eigenvects_sympy(M, iszerofunc, simplify, flags)
  if (primitive) ret <- lapply(ret, function(e) {
    e$vectors <- lapply(e$vectors, function(v) {
      g <- Reduce(gcd, v$data)
      .mat_map(v, function(z) simpfunc(.div(z, g)))
    })
    e
  })
  if (has_floats) ret <- lapply(ret, function(e) {
    e$value <- evalf(e$value)
    e$vectors <- lapply(e$vectors, function(v) .mat_map(v, evalf))
    e
  })
  ret
}

# _is_diagonalizable_with_eigen
.matint_diag_with_eigen <- function(M, reals_only = FALSE) {
  if (M$nrow != M$ncol) return(list(FALSE, list()))
  ev <- .matint_eigenvects(M, simplify = TRUE)
  for (e in ev) {
    if (reals_only && !isTRUE(ask_fact(e$value, "real"))) return(list(FALSE, ev))
    if (e$multiplicity != length(e$vectors)) return(list(FALSE, ev))
  }
  list(TRUE, ev)
}

.matint_all_fact <- function(M, fact) all(vapply(M$data, function(z) isTRUE(ask_fact(z, fact)), TRUE))

# _is_diagonalizable
.matint_is_diagonalizable <- function(M, reals_only = FALSE) {
  M <- Matrix(M)
  if (M$nrow != M$ncol) return(FALSE)
  if (.matint_all_fact(M, "real") && isTRUE(is_symmetric(M))) return(TRUE)
  if (.matint_all_fact(M, "complex") && isTRUE(.matint_hermitian(M))) return(TRUE)
  .matint_diag_with_eigen(M, reals_only)[[1L]]
}

# _diagonalize
.matint_diagonalize <- function(M, reals_only = FALSE, sort = FALSE, normalize = FALSE) {
  M <- Matrix(M)
  if (M$nrow != M$ncol) stop(NonSquareMatrixError(""))
  de <- .matint_diag_with_eigen(M, reals_only)
  if (!de[[1L]]) stop(MatrixError("Matrix is not diagonalizable"))
  ev <- de[[2L]]
  if (sort) {
    ord <- unlist(.key_sort(lapply(ev, function(e) .sort_key(e$value)), as.list(seq_along(ev))))
    ev <- ev[ord]
  }
  cols <- list(); dg <- list()
  for (e in ev) {
    dg <- c(dg, rep(list(e$value), e$multiplicity))
    cols <- c(cols, e$vectors)
  }
  if (normalize) cols <- lapply(cols, function(v) .mat_map(v, function(z) .div(z, norm_matrix(v))))
  list(P = do.call(cbind.symmatrix, cols), D = diag_matrix(dg))
}

#' Left eigenvectors and singular values
#'
#' `left_eigenvects()` returns the eigenvectors of the transpose as row
#' vectors, and `singular_values()` returns the square roots of the
#' eigenvalues of `M^H M` (or `M M^H` for a wide matrix) in decreasing order.
#'
#' @param M A matrix.
#' @param ... Flags passed to [eigenvects()].
#' @return For `left_eigenvects()`, a list of `list(value, multiplicity,
#'   vectors)`; for `singular_values()`, a list of expressions.
#' @examples
#' left_eigenvects(Matrix(list(list(0, 1), list(1, 0))))
#' singular_values(Matrix(list(list(1, 0), list(0, 2))))
#' @export
left_eigenvects <- function(M, ...) {
  ev <- eigenvects(transpose(Matrix(M)), ...)
  lapply(ev, function(e) { e$vectors <- lapply(e$vectors, transpose); e })
}

#' @rdname left_eigenvects
#' @export
singular_values <- function(M) {
  M <- Matrix(M)
  H <- transpose(.mat_map(M, conjugate))
  vp <- if (M$nrow >= M$ncol) eigenvals(.matmul(H, M)) else eigenvals(.matmul(M, H))
  vals <- list()
  for (k in seq_along(vp)) vals <- c(vals, rep(list(sqrt(vp[[k]])), attr(vp, "multiplicity")[k]))
  if (length(vals) < M$ncol) vals <- c(vals, rep(list(S_Zero), M$ncol - length(vals)))
  ord <- unlist(.key_sort(lapply(vals, .sort_key), as.list(seq_along(vals))))
  rev(vals[ord])
}
