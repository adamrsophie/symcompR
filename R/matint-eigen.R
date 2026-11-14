# Eigenvalues after sympy/matrices/eigen.py: _eigenvals, _eigenvals_dict and
# _eigenvals_list. The characteristic polynomial of each strongly connected
# block is factored with factor_list(), each factor is solved with roots()
# and Poly.all_roots() takes over when roots() does not find every root.

# M._rep.domain in (ZZ, QQ): a Matrix whose entries are all rational
.matint_is_qq <- function(M) all(vapply(M$data, function(z) z$type %in% c("Integer", "Rational"), TRUE))

# Matrix.strongly_connected_components(): DomainMatrix.scc() on the sparse
# rep, whose rows list their nonzero columns in increasing order
.matint_scc <- function(M) {
  n <- M$nrow
  V <- as.character(seq_len(n) - 1L)
  Gmap <- stats::setNames(lapply(seq_len(n), function(i) {
    nz <- which(!vapply(seq_len(n), function(j) is_zero_num(.mget(M, i, j)), TRUE))
    as.character(nz - 1L)
  }), V)
  lapply(.py_scc(V, Gmap), function(b) as.integer(b) + 1L)
}

# Matrix.charpoly() as an expression in lam: for the EXRAW domain (entries
# that are not polynomials over a ground domain) each coefficient goes
# through the simplify function, as _charpoly does
.matint_charpoly <- function(M, lam, simpfunc = NULL) {
  cf <- .charpoly_coeffs(M)
  dom <- tryCatch(construct_domain(M$data)[[1L]], error = function(e) NULL)
  exraw <- is.null(dom) || identical(.de_kind(dom), "EX") || identical(.de_kind(dom), "EXRAW")
  if (!is.null(simpfunc)) cf <- lapply(cf, simpfunc)
  else if (exraw) cf <- lapply(cf, simplify)
  n <- length(cf) - 1L
  expand(.add(lapply(seq_along(cf), function(k) .mul(list(cf[[k]], .pow(lam, .int(n - k + 1L)))))))
}

.matint_eigenvals_error <- paste0(
  "It is not always possible to express the eigenvalues of a matrix ",
  "of size 5x5 or higher in radicals. ",
  "We have CRootOf, but domains other than the rationals are not ",
  "currently supported. ",
  "If there are no symbols in the matrix, ",
  "it should still be possible to compute numeric approximations ",
  "of the eigenvalues using ",
  "M.evalf().eigenvals() or M.charpoly().nroots().")

# _eigenvals_dict: a list of the eigenvalues with attribute "multiplicity"
.matint_eigenvals_dict <- function(M, error_when_incomplete = TRUE, simplify = FALSE, flags = list()) {
  keys <- list(); mult <- integer(); kkeys <- character()
  exp_keys <- character(); exp_vals <- list()
  add_eig <- function(k, v) {
    i <- match(k$key, kkeys)
    if (is.na(i)) {
      keys[[length(keys) + 1L]] <<- k; mult <<- c(mult, as.integer(v)); kkeys <<- c(kkeys, k$key)
    } else mult[i] <<- mult[i] + as.integer(v)
  }
  is_dom <- .matint_is_qq(M)
  simpfunc <- if (is.function(simplify)) simplify else NULL
  for (b in .matint_scc(M)) {
    if (is_dom && length(b) == 1L) {
      add_eig(.mget(M, b, b), 1L)
      next
    }
    block <- M[b, b, drop = FALSE]
    lam <- Dummy("lambda")
    cp <- .matint_charpoly(block, lam, simpfunc)
    factors <- factor_list(cp, lam)$factors
    for (fm in factors) {
      fac <- fm[[1L]]; multiplicity <- as.integer(fm[[2L]])
      eigs <- do.call(roots, c(list(fac, lam), flags))
      em <- attr(eigs, "multiplicity")
      if (is.null(em)) em <- rep(1L, length(eigs))
      deg <- as.integer(degree(fac, lam))
      if (sum(em) != deg) {
        ar <- tryCatch(.all_roots_rootof(fac, lam), error = function(e) NULL)
        if (is.null(ar)) {
          if (error_when_incomplete) stop(MatrixError(.matint_eigenvals_error))
          eigs <- list(); em <- integer()
        } else {
          # dict(factor.all_roots(multiple=False))
          eigs <- list(); em <- integer(); ek <- character()
          for (r in ar) {
            i <- match(r$key, ek)
            if (is.na(i)) { eigs[[length(eigs) + 1L]] <- r; em <- c(em, 1L); ek <- c(ek, r$key) }
            else em[i] <- em[i] + 1L
          }
        }
      }
      for (j in seq_along(eigs)) {
        k <- eigs[[j]]
        # canonicalize equal eigenvalues through their expanded form
        ke <- expand(k)$key
        i <- match(ke, exp_keys)
        if (!is.na(i)) k <- exp_vals[[i]]
        else { exp_keys <- c(exp_keys, ke); exp_vals[[length(exp_vals) + 1L]] <- k }
        add_eig(k, em[j] * multiplicity)
      }
    }
  }
  if (!isFALSE(simplify)) {
    sf <- if (is.function(simplify)) simplify else .simplify_default
    keys <- lapply(keys, sf)
  }
  structure(keys, multiplicity = mult)
}

.simplify_default <- function(e) simplify(e)

# _eigenvals
.matint_eigenvals <- function(M, error_when_incomplete = TRUE, simplify = FALSE, multiple = FALSE,
                              rational = FALSE, flags = list()) {
  M <- Matrix(M)
  if (M$nrow == 0L || M$ncol == 0L) return(if (multiple) list() else structure(list(), multiplicity = integer()))
  if (M$nrow != M$ncol) stop(NonSquareMatrixError(paste0(.matint_mat_str(M), " must be a square matrix.")))
  if (rational) M <- .mat_map(M, function(z) if (.has_type(z, "Float")) nsimplify(z, rational = TRUE) else z)
  d <- .matint_eigenvals_dict(M, error_when_incomplete, simplify, flags)
  if (!multiple) return(d)
  unlist(lapply(seq_along(d), function(i) rep(list(d[[i]]), attr(d, "multiplicity")[i])), recursive = FALSE)
}

# str(Matrix) on one line, as SymPy prints it inside error messages
.matint_mat_str <- function(M) {
  rows <- vapply(seq_len(M$nrow), function(i)
    paste0("[", paste(vapply(seq_len(M$ncol), function(j) sstr(.mget(M, i, j)), ""), collapse = ", "), "]"), "")
  paste0("Matrix([", paste(rows, collapse = ", "), "])")
}
