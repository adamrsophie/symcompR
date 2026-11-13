# Matrix expressions, after sympy/matrices/expressions: symbolic matrices
# with shapes, and the canonical MatAdd, MatMul and MatPow products.

.MATEXPR <- c("MatrixSymbol", "Identity", "ZeroMatrix", "OneMatrix", "MatAdd", "MatMul", "MatPow",
              "Transpose", "Inverse", "Adjoint", "HadamardProduct", "HadamardPower", "KroneckerProduct",
              "BlockMatrix", "BlockDiagMatrix", "FunctionMatrix", "DiagMatrix", "DiagonalOf",
              "DiagonalMatrix", "PermutationMatrix", "MatrixPermute", "ImmutableMatrix", "MatrixSlice", "MatrixSolve",
              "ElementwiseApplyFunction", "MatrixUnit", "GenericIdentity", "GenericZeroMatrix", "CompanionMatrix",
              "DFT", "IDFT", "LofLU", "UofLU", "LofCholesky", "UofCholesky", "QofQR", "RofQR", "EigenVectors",
              "EigenValues", "UofSVD", "SofSVD", "VofSVD", "ExpectationMatrix", "VarianceMatrix",
              "CrossCovarianceMatrix")
.is_matexpr <- function(e) is_sym(e) && e$type %in% .MATEXPR

.mshape <- function(e) {
  if (inherits(e, "symmatrix")) return(list(.int(e$nrow), .int(e$ncol)))
  a <- e$args
  switch(e$type,
    MatrixSymbol = , ZeroMatrix = , OneMatrix = , FunctionMatrix = , ImmutableMatrix = list(a[[1L]], a[[2L]]),
    MatrixSlice = lapply(2:3, function(k) {
      s <- a[[k]]$args
      d <- .sub(s[[2L]], s[[1L]])
      if (is_one_num(s[[3L]])) d else fcall("ceiling", list(.div(d, s[[3L]])))
    }),
    Identity = , PermutationMatrix = list(a[[1L]], a[[1L]]),
    MatAdd = , HadamardProduct = , MatPow = , DiagonalMatrix = , Inverse = ,
    MatrixPermute = .mshape(a[[1L]]),
    MatrixSolve = .mshape(a[[2L]]),
    MatMul = {
      ms <- Filter(.is_matexpr, a)
      list(.mshape(ms[[1L]])[[1L]], .mshape(ms[[length(ms)]])[[2L]])
    },
    Transpose = , Adjoint = rev(.mshape(a[[1L]])),
    KroneckerProduct = {
      sh <- lapply(a, .mshape)
      list(.mul(lapply(sh, `[[`, 1L)), .mul(lapply(sh, `[[`, 2L)))
    },
    BlockMatrix = {
      rows <- a[[1L]]$args
      list(.add(lapply(rows, function(r) .mshape(r$args[[1L]])[[1L]])),
           .add(lapply(rows[[1L]]$args, function(b) .mshape(b)[[2L]])))
    },
    BlockDiagMatrix = list(.add(lapply(a, function(b) .mshape(b)[[1L]])), .add(lapply(a, function(b) .mshape(b)[[2L]]))),
    DiagMatrix = {
      s <- .mshape(a[[1L]])
      k <- if (identical(s[[1L]]$key, S_One$key)) s[[2L]] else s[[1L]]
      list(k, k)
    },
    DiagonalOf = {
      s <- .mshape(a[[1L]])
      list(if (identical(s[[1L]]$key, s[[2L]]$key)) s[[1L]] else Min(s[[1L]], s[[2L]]), S_One)
    },
    { h <- .me_shape_hooks[[e$type]]; if (!is.null(h)) h(e) })
}

.dim_eq <- function(a, b) identical(a$key, b$key)
.is_square_expr <- function(e) { s <- .mshape(e); .dim_eq(s[[1L]], s[[2L]]) }
.dim_known_ne <- function(a, b) is_exact_number(a) && is_exact_number(b) && !.dim_eq(a, b)
.shape_error <- function(...) stop("ShapeError: ", ..., call. = FALSE)
.mdims <- function(x) if (is.numeric(x)) .int(x) else as_sym(x)

#' Matrix expressions
#'
#' `MatrixSymbol(name, n, m)` is a symbolic `n` by `m` matrix. `Identity(n)`,
#' `ZeroMatrix(m, n)` and `OneMatrix(m, n)` are the special matrices. `+`,
#' `-`, `*` and `^` build the canonical `MatAdd`, `MatMul` and `MatPow`
#' (like terms combine, identities and zeros simplify, shapes are checked
#' and mismatches raise a `ShapeError`). [t()] transposes and distributes;
#' `Transpose()`, `Inverse()` and `Adjoint()` build unevaluated nodes
#' (`A^-1` is `Inverse(A)`). `Trace()`, `Determinant()` and `Permanent()`
#' are scalar functions of a matrix; [doit()] applies their rules. `MatrixElement(A, i, j)`
#' is the entry with SymPy's zero based indices. `as_explicit()` turns an
#' expression with integer shape into an explicit [Matrix()].
#' `mat_shape()` gives the shape and `matrix_symbols()` the matrix symbols.
#'
#' @param name The name.
#' @param n,m Dimensions.
#' @param ... Matrices, or for `MatMul()` matrices and scalars.
#' @param e,A,B A matrix expression.
#' @param k An exponent.
#' @param i,j Zero based indices.
#' @param f A [Lambda()] of two variables giving the entries.
#' @return A matrix expression, a scalar expression or an explicit matrix.
#' @examples
#' A <- MatrixSymbol("A", 2, 2); B <- MatrixSymbol("B", 2, 2)
#' A * B + A * B
#' t(A * B)
#' A^-1 * A
#' as_explicit(A * B)
#' doit(Trace(A + B))
#' @name matrix-expressions
NULL

#' @rdname matrix-expressions
#' @export
MatrixSymbol <- function(name, n, m) {
  n <- .mdims(n); m <- .mdims(m)
  .me_check_dim(m); .me_check_dim(n)
  new_sym("MatrixSymbol", list(n, m), name = name,
          key = paste0("MS", nchar(name), ":", name, ":", n$key, ":", m$key))
}

#' @rdname matrix-expressions
#' @export
Identity <- function(n) new_sym("Identity", list(.me_check_dim(.mdims(n))))
#' @rdname matrix-expressions
#' @export
ZeroMatrix <- function(n, m) {
  n <- .mdims(n); m <- .mdims(m)
  .me_check_dim(n); .me_check_dim(m)
  new_sym("ZeroMatrix", list(n, m))
}
#' @rdname matrix-expressions
#' @export
OneMatrix <- function(n, m) {
  n <- .mdims(n); m <- .mdims(m)
  .me_check_dim(n); .me_check_dim(m)
  new_sym("OneMatrix", list(n, m))
}

#' @rdname matrix-expressions
#' @export
mat_shape <- function(e) .mshape(as_matexpr(e))

as_matexpr <- function(e) if (inherits(e, "symmatrix") || .is_matexpr(e)) e else as_sym(e)

# base and exponent of a matrix factor
.mbase_exp <- function(m) {
  if (m$type == "MatPow") return(list(m$args[[1L]], m$args[[2L]]))
  if (m$type == "Inverse") return(list(m$args[[1L]], S_NegativeOne))
  list(m, S_One)
}

# An explicit matrix inside a matrix expression is an ImmutableMatrix node.
.is_explicit <- function(a) inherits(a, "symmatrix") || (is_sym(a) && a$type == "ImmutableMatrix")
.as_mexpr <- function(a) if (inherits(a, "symmatrix")) .as_block(a) else a
# a lone explicit result is an explicit matrix
.unwrap_explicit <- function(a) if (is_sym(a) && a$type == "ImmutableMatrix") as_explicit(a) else a

#' @rdname matrix-expressions
#' @param evaluate If `TRUE`, `MatAdd()` and `MatMul()` put their arguments
#'   in canonical form, as [doit()] does; by default they keep them as
#'   given, as in SymPy. A `MatMul()` of scalars only is their product.
#' @export
MatMul <- function(..., evaluate = FALSE) {
  if (!...length()) return(GenericIdentity())
  if (!isTRUE(evaluate)) {
    args <- Filter(function(a) !(is_sym(a) && a$type == "GenericIdentity"), lapply(lapply(list(...), as_matexpr), .as_mexpr))
    if (!any(vapply(args, .is_matexpr, NA))) return(.mul(args))
    return(new_sym("MatMul", args))
  }
  flat <- list()
  for (a in list(...)) {
    a <- as_matexpr(a)
    if (is_sym(a) && a$type == "MatMul") flat <- c(flat, a$args) else flat[[length(flat) + 1L]] <- a
  }
  ismat <- vapply(flat, function(a) .is_matexpr(a) || inherits(a, "symmatrix"), TRUE)
  ex <- vapply(flat, .is_explicit, TRUE)
  if (any(ex) && !all(ex[ismat])) {
    # MatMul's merge_explicit: neighbouring explicit matrices and numbers
    # multiply out
    en <- function(z) .is_explicit(z) || (is_sym(z) && is_number_type(z))
    times <- function(u, v) {
      if (!.is_explicit(u)) return(.mat_map(as_explicit(v), function(z) .mul(list(u, z))))
      if (!.is_explicit(v)) return(.mat_map(as_explicit(u), function(z) .mul(list(z, v))))
      .matmul(as_explicit(u), as_explicit(v))
    }
    out <- list(); last <- flat[[1L]]
    for (a in flat[-1L]) {
      if (en(a) && en(last)) last <- if (!.is_explicit(a) && !.is_explicit(last)) .mul(list(last, a)) else times(last, a)
      else { out[[length(out) + 1L]] <- last; last <- a }
    }
    flat <- lapply(c(out, list(last)), .as_mexpr)
  }
  # MatMul's distribute_monom: a rational times a MatAdd distributes
  if (length(flat) == 2L) {
    ra <- vapply(flat, function(a) is_sym(a) && a$type %in% c("Integer", "Rational"), TRUE)
    ad <- vapply(flat, function(a) is_sym(a) && a$type == "MatAdd", TRUE)
    if (any(ra) && any(ad) && which(ra) != which(ad))
      return(do.call(.madd, lapply(flat[ad][[1L]]$args, function(t) if (ra[1L]) .mm(flat[[1L]], t) else .mm(t, flat[[2L]]))))
  }
  coef <- S_One
  mats <- list()
  for (a in flat) if (.is_matexpr(a) || inherits(a, "symmatrix")) mats[[length(mats) + 1L]] <- a else coef <- .mul(list(coef, a))
  if (!length(mats)) return(coef)
  if (all(vapply(mats, .is_explicit, TRUE))) {
    P <- Reduce(.matmul, lapply(mats, as_explicit))
    return(.mat_map(P, function(z) .mul(list(coef, z))))
  }
  if (length(mats) > 1L) for (i in 2:length(mats)) {
    c1 <- .mshape(mats[[i - 1L]])[[2L]]; r2 <- .mshape(mats[[i]])[[1L]]
    if (.dim_known_ne(c1, r2)) .shape_error("Matrices are not aligned (", .str(c1), ", ", .str(r2), ")")
  }
  shp <- list(.mshape(mats[[1L]])[[1L]], .mshape(mats[[length(mats)]])[[2L]])
  if (is_zero_num(coef) || any(vapply(mats, function(m) m$type == "ZeroMatrix", TRUE))) return(ZeroMatrix(shp[[1L]], shp[[2L]]))
  # merge equal neighbours into powers, X^-1 * X cancels
  out <- list()
  # remove_ids: identities, and ones matrices known to be 1x1 (x.is_Identity,
  # which an explicit identity matrix has too)
  is_id <- function(m) m$type == "Identity" || (m$type == "OneMatrix" && is_one_num(m$args[[1L]]) && is_one_num(m$args[[2L]])) ||
    (.is_explicit(m) && .me_is_identity(m))
  for (m in Filter(Negate(is_id), mats)) {
    if (length(out)) {
      last <- out[[length(out)]]
      # combine_one_matrices: OneMatrix(a, b)*OneMatrix(b, c) is b*OneMatrix(a, c)
      if (last$type == "OneMatrix" && m$type == "OneMatrix") {
        coef <- .mul(list(coef, last$args[[2L]]))
        out[[length(out)]] <- OneMatrix(last$args[[1L]], m$args[[2L]])
        next
      }
      be1 <- .mbase_exp(last); be2 <- .mbase_exp(m)
      if (identical(be1[[1L]]$key, be2[[1L]]$key)) {
        p <- .add(list(be1[[2L]], be2[[2L]]))
        out[[length(out)]] <- NULL
        if (!is_zero_num(p)) out[[length(out) + 1L]] <- MatPow(be1[[1L]], p)
        next
      }
      # combine_powers: B_base.inverse() == A_base cancels too
      if (!.is_explicit(be2[[1L]])) {
        inv <- tryCatch(.me_inverse(be2[[1L]]), error = function(err) NULL)
        if (is_sym(inv) && identical(inv$key, be1[[1L]]$key)) {
          p <- .sub(be1[[2L]], be2[[2L]])
          out[[length(out)]] <- NULL
          if (!is_zero_num(p)) out[[length(out) + 1L]] <- MatPow(be1[[1L]], p)
          next
        }
      }
    }
    out[[length(out) + 1L]] <- m
  }
  if (!length(out)) out <- list(Identity(shp[[1L]]))
  if (identical(coef$key, S_One$key) && length(out) == 1L) return(.unwrap_explicit(out[[1L]]))
  new_sym("MatMul", c(if (!identical(coef$key, S_One$key)) list(coef), out))
}

# MatMul and MatAdd in canonical form, for internal use
.mm <- function(...) MatMul(..., evaluate = TRUE)
.madd <- function(...) MatAdd(..., evaluate = TRUE)

.mterm_split <- function(t) {
  if (t$type == "MatMul" && !.is_matexpr(t$args[[1L]])) {
    rest <- t$args[-1L]
    return(list(t$args[[1L]], if (length(rest) == 1L) rest[[1L]] else new_sym("MatMul", rest)))
  }
  list(S_One, t)
}

#' @rdname matrix-expressions
#' @export
MatAdd <- function(..., evaluate = FALSE) {
  if (!...length()) return(GenericZeroMatrix())
  if (!isTRUE(evaluate)) {
    args <- Filter(function(a) !(is_sym(a) && a$type == "GenericZeroMatrix"), lapply(lapply(list(...), as_matexpr), .as_mexpr))
    if (!all(vapply(args, .is_matexpr, TRUE))) stop(.sympy_condition("TypeError", "Mix of Matrix and Scalar symbols"))
    return(new_sym("MatAdd", args))
  }
  flat <- list()
  for (a in list(...)) {
    a <- as_matexpr(a)
    if (is_sym(a) && a$type == "MatAdd") flat <- c(flat, a$args) else flat[[length(flat) + 1L]] <- a
  }
  ex <- vapply(flat, .is_explicit, TRUE)
  if (all(ex)) return(Reduce(`+`, lapply(flat, as_explicit)))
  # MatAdd's merge_explicit: the explicit terms add up to one
  if (any(ex)) flat <- c(flat[!ex], list(.as_block(Reduce(`+`, lapply(flat[ex], as_explicit)))))
  if (!all(vapply(flat, .is_matexpr, TRUE))) stop("cannot add a matrix and a scalar", call. = FALSE)
  s1 <- .mshape(flat[[1L]])
  for (f in flat[-1L]) {
    s <- .mshape(f)
    if (.dim_known_ne(s[[1L]], s1[[1L]]) || .dim_known_ne(s[[2L]], s1[[2L]]))
      .shape_error("Matrices ", .str(flat[[1L]]), " and ", .str(f), " are not aligned")
  }
  coefs <- list(); parts <- list()
  for (f in Filter(function(f) f$type != "ZeroMatrix", flat)) {
    sp <- .mterm_split(f)
    k <- sp[[2L]]$key
    if (is.null(parts[[k]])) { parts[[k]] <- sp[[2L]]; coefs[[k]] <- sp[[1L]] }
    else coefs[[k]] <- .add(list(coefs[[k]], sp[[1L]]))
  }
  terms <- list()
  for (k in names(parts)) if (!is_zero_num(coefs[[k]])) terms[[length(terms) + 1L]] <-
    if (parts[[k]]$type == "ImmutableMatrix") .as_block(.mat_map(as_explicit(parts[[k]]), function(z) .mul(list(coefs[[k]], z))))
    else .mm(coefs[[k]], parts[[k]])
  if (!length(terms)) return(ZeroMatrix(s1[[1L]], s1[[2L]]))
  if (length(terms) == 1L) return(.unwrap_explicit(terms[[1L]]))
  # canonicalize's last rule: sort(default_sort_key)
  new_sym("MatAdd", .key_sort(lapply(terms, .sort_key), terms))
}

#' @rdname matrix-expressions
#' @export
MatPow <- function(A, k) {
  A <- as_matexpr(A); k <- as_sym(k)
  if (inherits(A, "symmatrix")) {
    if (k$type != "Integer") stop("MatPow of an explicit matrix needs an integer power", call. = FALSE)
    M <- if (.int_val(k) < 0) inverse(A) else A
    return(if (.int_val(k) == 0) eye(A$nrow) else Reduce(.matmul, rep(list(M), abs(.int_val(k)))))
  }
  if (A$type == "GenericIdentity") return(A)
  s <- .mshape(A)
  if (.dim_known_ne(s[[1L]], s[[2L]])) .perr("NonSquareMatrixError", paste0("Power of non-square matrix ", .str(A)))
  if (is_zero_num(k)) return(Identity(s[[1L]]))
  if (identical(k$key, S_One$key)) return(A)
  if (A$type == "Identity") return(A)
  if (A$type == "ZeroMatrix" && isTRUE(.facts(k)[["positive"]])) return(A)
  # ZeroMatrix._eval_inverse ends its message with a period, ZeroMatrix._eval_power does not
  if (A$type == "ZeroMatrix" && isTRUE(.facts(k)[["negative"]]))
    .perr("NonInvertibleMatrixError", paste0("Matrix det == 0; not invertible", if (identical(k$key, S_NegativeOne$key)) "."))
  if (A$type %in% c("MatPow", "Inverse") && k$type == "Integer") {
    be <- .mbase_exp(A)
    return(MatPow(be[[1L]], .mul(list(be[[2L]], k))))
  }
  if (identical(k$key, S_NegativeOne$key)) return(new_sym("Inverse", list(A)))
  new_sym("MatPow", list(A, k))
}

#' @rdname matrix-expressions
#' @export
Transpose <- function(A) new_sym("Transpose", list(.as_mexpr(as_matexpr(A))))
#' @rdname matrix-expressions
#' @export
Inverse <- function(A) {
  A <- as_matexpr(A)
  if (inherits(A, "symmatrix")) return(inverse(A))
  s <- .mshape(A)
  if (.dim_known_ne(s[[1L]], s[[2L]])) .perr("NonSquareMatrixError", paste0("Inverse of non-square matrix ", .str(A)))
  new_sym("Inverse", list(A))
}
#' @rdname matrix-expressions
#' @export
Adjoint <- function(A) {
  A <- as_matexpr(A)
  if (inherits(A, "symmatrix")) return(.mat_map(t(A), function(z) fcall("conjugate", list(z))))
  new_sym("Adjoint", list(A))
}

# The adjoint() function: MatAdd._eval_adjoint and MatMul._eval_adjoint
# adjoin the arguments (a MatMul in reverse order) and doit(); other
# matrix expressions give an Adjoint node.
.madjoint <- function(e) {
  if (inherits(e, "symmatrix")) return(.mat_map(t(e), function(z) fcall("conjugate", list(z))))
  if (!.is_matexpr(e)) return(fcall("conjugate", list(e)))
  switch(e$type,
    ImmutableMatrix = .madjoint(as_explicit(e)),
    Adjoint = .unwrap_explicit(e$args[[1L]]),
    MatAdd = do.call(.madd, lapply(e$args, function(a) .as_mexpr(.madjoint(a)))),
    MatMul = do.call(.mm, lapply(rev(e$args), function(a) .as_mexpr(.madjoint(a)))),
    { h <- .me_adjoint_hooks[[e$type]]; if (is.null(h)) Adjoint(e) else h(e) })
}

# Evaluating transpose, as SymPy's .T.
.mtranspose <- function(e) {
  if (inherits(e, "symmatrix")) return(t(e))
  switch(e$type,
    Transpose = .unwrap_explicit(e$args[[1L]]),
    ImmutableMatrix = t(as_explicit(e)),
    Identity = , DiagMatrix = e,
    ZeroMatrix = ZeroMatrix(e$args[[2L]], e$args[[1L]]),
    OneMatrix = OneMatrix(e$args[[2L]], e$args[[1L]]),
    MatAdd = do.call(.madd, lapply(e$args, .mtranspose)),
    MatMul = {
      sp <- .mterm_split(e)
      ms <- if (sp[[2L]]$type == "MatMul") sp[[2L]]$args else list(sp[[2L]])
      do.call(.mm, c(list(sp[[1L]]), rev(lapply(ms, .mtranspose))))
    },
    MatPow = MatPow(.mtranspose(e$args[[1L]]), e$args[[2L]]),
    Inverse = Inverse(.mtranspose(e$args[[1L]])),
    HadamardProduct = do.call(HadamardProduct, lapply(e$args, .mtranspose)),
    BlockMatrix = {
      rows <- lapply(e$args[[1L]]$args, function(r) r$args)
      nr <- length(rows); nc <- length(rows[[1L]])
      BlockMatrix(lapply(seq_len(nc), function(j) lapply(seq_len(nr), function(i) .mtranspose(rows[[i]][[j]]))))
    },
    { h <- .me_transpose_hooks[[e$type]]; if (is.null(h)) Transpose(e) else h(e) })
}

#' @export
t.sym <- function(x) {
  if (!.is_matexpr(x)) stop("t() needs a matrix", call. = FALSE)
  .mtranspose(x)
}

#' @rdname matrix-expressions
#' @export
Trace <- function(A) {
  A <- .as_mexpr(as_matexpr(A))
  s <- .mshape(A)
  if (.dim_known_ne(s[[1L]], s[[2L]])) .perr("NonSquareMatrixError", "Trace of a non-square matrix")
  new_sym("Trace", list(A))
}
.trace_explicit <- function(M) .add(lapply(seq_len(min(M$nrow, M$ncol)), function(i) .mget(M, i, i)))

#' @rdname matrix-expressions
#' @export
Determinant <- function(A) {
  A <- as_matexpr(A)
  if (inherits(A, "symmatrix")) return(det(A))
  s <- .mshape(A)
  if (.dim_known_ne(s[[1L]], s[[2L]])) .perr("NonSquareMatrixError", "Det of a non-square matrix")
  new_sym("Determinant", list(A))
}

#' @rdname matrix-expressions
#' @export
Permanent <- function(A) {
  A <- as_matexpr(A)
  if (inherits(A, "symmatrix")) return(permanent(A))
  new_sym("Permanent", list(A))
}

.trace_doit <- function(A) {
  switch(A$type,
    MatAdd = .add(lapply(A$args, .trace_doit)),
    MatMul = {
      sp <- .mterm_split(A)
      if (!identical(sp[[1L]]$key, S_One$key)) .mul(list(sp[[1L]], .trace_doit(sp[[2L]]))) else Trace(A)
    },
    Identity = A$args[[1L]],
    ZeroMatrix = S_Zero,
    Transpose = .trace_doit(A$args[[1L]]),
    KroneckerProduct = .mul(lapply(A$args, .trace_doit)),
    MatPow = if (.is_explicit(A$args[[1L]])) .trace_explicit(as_explicit(A)) else Trace(A),
    {
      s <- .mshape(A)
      if (is_exact_number(s[[1L]]) && A$type != "MatrixSymbol") .trace_explicit(as_explicit(A)) else Trace(A)
    })
}

.det_doit <- function(A) {
  switch(A$type,
    Identity = S_One,
    ZeroMatrix = S_Zero,
    Transpose = .det_doit(A$args[[1L]]),
    Inverse = .div(S_One, .det_doit(A$args[[1L]])),
    MatPow = Determinant(A),
    MatMul = {
      sp <- .mterm_split(A)
      ms <- if (sp[[2L]]$type == "MatMul") sp[[2L]]$args else list(sp[[2L]])
      n <- .mshape(A)[[1L]]
      # factor determinants need square factors
      if (!all(vapply(ms, .is_square_expr, TRUE))) return(.mul(list(.pow(sp[[1L]], n), Determinant(sp[[2L]]))))
      .mul(c(list(.pow(sp[[1L]], n)), lapply(ms, .det_doit)))
    },
    MatrixSymbol = , MatAdd = Determinant(A),
    KroneckerProduct = .me_kron_det(A),
    {
      s <- .mshape(A)
      if (is_exact_number(s[[1L]])) det(as_explicit(A)) else Determinant(A)
    })
}

register_doit("Trace", function(e) .trace_doit(doit(e$args[[1L]])))
register_doit("Determinant", function(e) .det_doit(doit(e$args[[1L]])))
register_doit("Permanent", function(e) {
  s <- .mshape(e$args[[1L]])
  if (!is_exact_number(s[[1L]]) || !is_exact_number(s[[2L]])) return(e)
  permanent(as_explicit(e$args[[1L]]))
})
register_doit("Transpose", function(e) .mtranspose(doit(e$args[[1L]])))
# Inverse.doit: arg.inverse()
register_doit("Inverse", function(e) .me_inverse(doit(e$args[[1L]])))

#' @rdname matrix-expressions
#' @export
HadamardProduct <- function(..., evaluate = FALSE) .me_hadamard_new(list(...), evaluate)
#' @rdname matrix-expressions
#' @export
HadamardPower <- function(A, k) .me_hadamard_power_new(A, k)
#' @rdname matrix-expressions
#' @export
KroneckerProduct <- function(...) .me_kronecker_new(list(...))

#' @rdname matrix-expressions
#' @export
MatrixElement <- function(A, i, j) {
  A <- as_matexpr(A); i <- as_sym(i); j <- as_sym(j)
  if (.is_explicit(A)) {
    # an explicit matrix gives its entry for integer indices (MatrixElement.__new__)
    if (is_exact_number(i) && is_exact_number(j)) {
      M <- as_explicit(A)
      ix <- .me_a2idx(.int_val(i), M$nrow); jx <- .me_a2idx(.int_val(j), M$ncol)
      return(.mget(M, ix + 1L, jx + 1L))
    }
    A <- .as_mexpr(A)
  }
  if ((is_exact_number(i) || is_exact_number(j)) && !.me_valid_index(A, i, j)) .perr("IndexError", "indices out of range")
  h <- .me_entry_hooks[[A$type]]
  if (!is.null(h)) return(h(A, i, j))
  # structured matrices give their entries without being materialized (_entry)
  if (A$type == "ZeroMatrix") return(S_Zero)
  if (A$type == "OneMatrix") return(S_One)
  if (A$type == "Identity") {
    if (is_exact_number(i) && is_exact_number(j)) return(if (i$value == j$value) S_One else S_Zero)
    return(fcall("KroneckerDelta", list(i, j)))
  }
  if (A$type == "MatAdd") return(.add(lapply(A$args, function(B) MatrixElement(B, i, j))))
  if (A$type == "MatrixSymbol" || !is_exact_number(i) || !is_exact_number(j) || !is_exact_number(.mshape(A)[[1L]]))
    return(new_sym("MatrixElement", list(A, i, j)))
  M <- tryCatch(as_explicit(A), ValueError = function(err) NULL)
  if (is.null(M)) return(new_sym("MatrixElement", list(A, i, j)))
  .mget(M, .int_val(i) + 1L, .int_val(j) + 1L)
}

#' @rdname matrix-expressions
#' @export
FunctionMatrix <- function(n, m, f) new_sym("FunctionMatrix", list(.mdims(n), .mdims(m), as_sym(f)))
#' @rdname matrix-expressions
#' @export
DiagMatrix <- function(A) new_sym("DiagMatrix", list(as_matexpr(A)))
#' @rdname matrix-expressions
#' @export
diagonalize_vector <- function(A) doit(DiagMatrix(A))
#' @rdname matrix-expressions
#' @export
DiagonalOf <- function(A) new_sym("DiagonalOf", list(as_matexpr(A)))
#' @rdname matrix-expressions
#' @export
DiagonalMatrix <- function(A) new_sym("DiagonalMatrix", list(as_matexpr(A)))
#' @rdname matrix-expressions
#' @param perm A [Permutation()] or a vector of zero based images.
#' @export
PermutationMatrix <- function(perm) {
  arr <- if (inherits(perm, "Permutation")) unclass(perm)$array else as.integer(perm)
  new_sym("PermutationMatrix", list(.int(length(arr)), do.call(Tuple, lapply(arr, .int))))
}
#' @rdname matrix-expressions
#' @param axis 0 to permute rows, 1 for columns.
#' @export
MatrixPermute <- function(A, perm, axis = 0) {
  arr <- if (inherits(perm, "Permutation")) unclass(perm)$array else as.integer(perm)
  new_sym("MatrixPermute", list(as_matexpr(A), do.call(Tuple, lapply(arr, .int)), .int(axis)))
}
#' @rdname matrix-expressions
#' @export
DotProduct <- function(A, B) {
  args <- list(as_matexpr(A), as_matexpr(B))
  for (k in 1:2) {
    if (!.me_is_matrix(args[[k]])) .perr("TypeError", paste0("Argument ", k, " of DotProduct is not a matrix"))
    if (!any(vapply(.mshape(args[[k]]), is_one_num, NA))) .perr("TypeError", paste0("Argument ", k, " of DotProduct is not a vector"))
  }
  keys <- lapply(args, function(a) sort(unique(vapply(.mshape(a), .k, ""))))
  if (!identical(keys[[1L]], keys[[2L]])) .perr("TypeError", "DotProduct arguments are not the same length")
  new_sym("DotProduct", lapply(args, .as_block))
}
register_doit("DotProduct", function(e) {
  # DotProduct.doit: the 1x1 product of the vectors, entry 0
  a <- e$args[[1L]]; b <- e$args[[2L]]
  sa <- .mshape(a); sb <- .mshape(b)
  rowa <- is_one_num(sa[[1L]])
  mul <- if (identical(vapply(sa, .k, ""), vapply(sb, .k, ""))) {
    if (rowa) .mm(a, .mtranspose(b)) else .mm(.mtranspose(a), b)
  } else {
    if (rowa) .mm(a, b) else .mm(.mtranspose(a), .mtranspose(b))
  }
  if (inherits(mul, "symmatrix")) return(mul$data[[1L]])
  if (.is_explicit(mul)) return(as_explicit(mul)$data[[1L]])
  .me_getitem(mul, 0L, 0L)
})

# SymPy's slice normalize(): an index, c(start, stop) or c(start, stop,
# step), with NULL or NA for a default end and negative values from the end.
.mslice_norm <- function(i, size) {
  if (is_sym(i) || length(i) == 1L) {
    i <- as_sym(i)
    if (.is_true(Lt(i, S_Zero))) i <- .add(list(i, size))
    return(list(i, .add(list(i, S_One)), S_One))
  }
  i <- as.list(i)
  if (length(i) == 2L) i[[3L]] <- 1
  dflt <- function(v, d) if (is.null(v) || (length(v) == 1L && !is_sym(v) && is.na(v))) d else as_sym(v)
  start <- dflt(i[[1L]], S_Zero); stop <- dflt(i[[2L]], size); step <- dflt(i[[3L]], S_One)
  if (.is_true(Lt(start, S_Zero))) start <- .add(list(start, size))
  if (.is_true(Lt(stop, S_Zero))) stop <- .add(list(stop, size))
  if (.is_true(Lt(.mul(list(.sub(stop, start), step)), S_One))) stop("IndexError: empty slice", call. = FALSE)
  list(start, stop, step)
}

#' @rdname matrix-expressions
#' @param rowslice,colslice A zero based index, `c(start, stop)` or
#'   `c(start, stop, step)`; `NA` (or `NULL` in a list) takes the default
#'   end, and negative values count from the end, as Python slices.
#' @export
MatrixSlice <- function(A, rowslice, colslice) {
  A <- .as_block(A)
  sh <- .mshape(A)
  rs <- .mslice_norm(rowslice, sh[[1L]]); cs <- .mslice_norm(colslice, sh[[2L]])
  if (.is_true(Gt(S_Zero, rs[[1L]])) || .is_true(Lt(sh[[1L]], rs[[2L]])) ||
      .is_true(Gt(S_Zero, cs[[1L]])) || .is_true(Lt(sh[[2L]], cs[[2L]]))) stop("IndexError: slice out of range", call. = FALSE)
  if (A$type == "MatrixSlice") {
    # a slice of a slice is a slice of the parent
    sos <- function(s, t) {
      start <- .add(list(s[[1L]], .mul(list(t[[1L]], s[[3L]]))))
      stp <- .add(list(s[[1L]], .mul(list(s[[3L]], t[[2L]]))))
      if (.is_true(Gt(stp, s[[2L]]))) stop("IndexError: slice out of range", call. = FALSE)
      list(start, stp, .mul(list(s[[3L]], t[[3L]])))
    }
    return(MatrixSlice(A$args[[1L]], sos(A$args[[2L]]$args, rs), sos(A$args[[3L]]$args, cs)))
  }
  new_sym("MatrixSlice", list(A, do.call(Tuple, rs), do.call(Tuple, cs)))
}

#' @rdname matrix-expressions
#' @param rowsizes,colsizes Block sizes along rows and columns.
#' @export
blockcut <- function(A, rowsizes, colsizes) {
  # exact (possibly symbolic) block boundaries
  bnds <- function(sz) {
    sz <- lapply(as.list(sz), as_sym)
    out <- list(); lo <- S_Zero
    for (z in sz) { hi <- .add(list(lo, z)); out[[length(out) + 1L]] <- list(lo, hi); lo <- hi }
    out
  }
  BlockMatrix(lapply(bnds(rowsizes), function(rb) lapply(bnds(colsizes), function(cb) MatrixSlice(A, rb, cb))))
}

#' @rdname matrix-expressions
#' @param blocks A list of rows, each a list of matrices.
#' @export
BlockMatrix <- function(blocks) {
  rows <- lapply(blocks, function(r) do.call(Tuple, lapply(if (is.list(r) && !is_sym(r) && !inherits(r, "symmatrix")) r else list(r), .as_block)))
  new_sym("BlockMatrix", list(do.call(Tuple, rows)))
}
#' @rdname matrix-expressions
#' @export
BlockDiagMatrix <- function(...) new_sym("BlockDiagMatrix", lapply(list(...), .as_block))

# An explicit matrix inside a block wrapper becomes a keyed ImmutableMatrix node.
.as_block <- function(b) {
  if (!inherits(b, "symmatrix")) return(as_matexpr(b))
  new_sym("ImmutableMatrix", list(.int(b$nrow), .int(b$ncol), do.call(Tuple, b$data)))
}

# SymPy's is_Identity for block matrices: identity diagonal blocks, zero elsewhere.
.is_identity_block <- function(e) {
  if (!(.is_matexpr(e) && e$type %in% c("BlockMatrix", "BlockDiagMatrix"))) return(FALSE)
  R <- .block_rows(e)
  if (length(R) != length(R[[1L]])) return(FALSE)
  all(vapply(seq_along(R), function(i) all(vapply(seq_along(R[[i]]), function(j) {
    b <- R[[i]][[j]]
    if (i == j) b$type == "Identity" else b$type == "ZeroMatrix"
  }, TRUE)), TRUE))
}

.block_rows <- function(e) {
  if (e$type == "BlockMatrix") return(lapply(e$args[[1L]]$args, function(r) r$args))
  if (e$type == "BlockDiagMatrix") {
    bs <- e$args
    return(lapply(seq_along(bs), function(i) lapply(seq_along(bs), function(j) {
      if (i == j) bs[[i]] else ZeroMatrix(.mshape(bs[[i]])[[1L]], .mshape(bs[[j]])[[2L]])
    })))
  }
  NULL
}

#' @rdname matrix-expressions
#' @export
block_collapse <- function(e) {
  r <- .block_collapse(e)
  # SymPy's bc_unpack: a 1x1 block matrix is its block
  if (.is_matexpr(r) && r$type == "BlockMatrix") {
    R <- .block_rows(r)
    if (length(R) == 1L && length(R[[1L]]) == 1L) return(R[[1L]][[1L]])
  }
  r
}

.block_collapse <- function(e) {
  e <- as_matexpr(e)
  if (!.is_matexpr(e) || !length(e$args)) return(e)
  args <- lapply(e$args, function(a) if (.is_matexpr(a)) .block_collapse(a) else a)
  isb <- vapply(args, function(a) .is_matexpr(a) && a$type %in% c("BlockMatrix", "BlockDiagMatrix"), TRUE)
  if (e$type == "MatMul" && any(isb)) {
    coef <- S_One; mats <- list()
    for (a in args) if (.is_matexpr(a)) mats[[length(mats) + 1L]] <- a else coef <- .mul(list(coef, a))
    # identity block matrices drop out of the product, as SymPy's remove_ids
    ids <- vapply(mats, .is_identity_block, TRUE)
    if (any(ids)) {
      shp <- .mshape(e)
      mats <- mats[!ids]
      if (!length(mats)) return(.mm(coef, Identity(shp[[1L]])))
      if (length(mats) == 1L) return(.block_collapse(.mm(coef, mats[[1L]])))
    }
    if (all(vapply(mats, function(a) a$type %in% c("BlockMatrix", "BlockDiagMatrix"), TRUE))) {
      sizes <- function(R, rows) vapply(if (rows) R else R[[1L]], function(b) .k((if (rows) .mshape(b[[1L]]) else .mshape(b))[[if (rows) 1L else 2L]]), "")
      # the column partition of each factor must match the next one's row partition
      for (i in seq_len(length(mats) - 1L))
        if (!identical(sizes(.block_rows(mats[[i]]), FALSE), sizes(.block_rows(mats[[i + 1L]]), TRUE))) return(sym_rebuild(e, args))
      R <- .block_rows(mats[[1L]])
      for (B in mats[-1L]) {
        C <- .block_rows(B)
        R <- lapply(seq_along(R), function(i) lapply(seq_along(C[[1L]]), function(j)
          do.call(.madd, lapply(seq_along(C), function(k) .mm(R[[i]][[k]], C[[k]][[j]])))))
      }
      return(BlockMatrix(lapply(R, function(r) lapply(r, function(b) .mm(coef, b)))))
    }
  }
  if (e$type == "MatAdd" && all(isb)) {
    Rs <- lapply(args, .block_rows)
    R <- Rs[[1L]]
    grid <- function(B) paste(vapply(B, function(r) paste(vapply(r, function(b) paste(vapply(.mshape(b), .k, ""), collapse = "x"), ""), collapse = ","), ""), collapse = ";")
    if (length(unique(vapply(Rs, grid, ""))) > 1L) return(sym_rebuild(e, args))
    return(BlockMatrix(lapply(seq_along(R), function(i) lapply(seq_along(R[[1L]]), function(j)
      do.call(.madd, lapply(Rs, function(B) B[[i]][[j]]))))))
  }
  if (e$type == "Transpose" && isb[[1L]]) return(.mtranspose(BlockMatrix(.block_rows(args[[1L]]))))
  # bc_matmul: a positive integer power is the product of its copies
  if (e$type == "MatPow" && isb[[1L]] && args[[2L]]$type == "Integer" && .int_val(args[[2L]]) > 0)
    return(.block_collapse(do.call(MatMul, rep(args[1L], .int_val(args[[2L]])))))
  # bc_inverse: expr.inverse() for a BlockDiagMatrix argument, blockinverse_1x1
  if (e$type == "Inverse" && isb[[1L]]) {
    if (args[[1L]]$type == "BlockDiagMatrix") return(.me_inverse(new_sym("Inverse", args)))
    R <- .block_rows(args[[1L]])
    if (length(R) == 1L && length(R[[1L]]) == 1L) return(BlockMatrix(list(list(.me_inverse(R[[1L]][[1L]])))))
  }
  sym_rebuild(e, args)
}

#' @rdname matrix-expressions
#' @export
matrix_symbols <- function(e) {
  out <- list()
  walk <- function(z) {
    if (z$type == "MatrixSymbol") out[[z$key]] <<- z
    for (a in z$args) walk(a)
  }
  walk(as_matexpr(e))
  out <- unname(out)
  out[order(vapply(out, function(z) z$name, ""), method = "radix")]
}

#' @rdname matrix-expressions
#' @export
as_explicit <- function(e) {
  if (inherits(e, "SymArray") || .is_arrexpr(e) || (is_sym(e) && e$type == "ImmutableDenseNDimArray"))
    return(.ae_explicit_arg(e))
  if (is_sym(e) && e$type == "Trace") return(doit(Trace(as_explicit(e$args[[1L]]))))
  e <- as_matexpr(e)
  if (inherits(e, "symmatrix")) return(e)
  s <- .mshape(e)
  if (!is_exact_number(s[[1L]]) || !is_exact_number(s[[2L]]))
    .perr("ValueError", "Matrix with symbolic shape cannot be represented explicitly.")
  nr <- .int_val(s[[1L]]); nc <- .int_val(s[[2L]])
  build <- function(f) .new_mat(unlist(lapply(seq_len(nr), function(i) lapply(seq_len(nc), function(j) f(i, j))),
                                       recursive = FALSE), nr, nc)
  a <- e$args
  switch(e$type,
    MatrixSymbol = build(function(i, j) new_sym("MatrixElement", list(e, .int(i - 1L), .int(j - 1L)))),
    ImmutableMatrix = .new_mat(a[[3L]]$args, nr, nc),
    MatrixSlice = {
      rs <- a[[2L]]$args; cs <- a[[3L]]$args
      P <- a[[1L]]
      at <- function(s, i) .add(list(.mul(list(.int(i - 1L), s[[3L]])), s[[1L]]))
      if (inherits(P, "symmatrix") || P$type == "ImmutableMatrix") {
        Pe <- as_explicit(P)
        build(function(i, j) .mget(Pe, .int_val(at(rs, i)) + 1L, .int_val(at(cs, j)) + 1L))
      } else build(function(i, j) MatrixElement(P, at(rs, i), at(cs, j)))
    },
    Identity = eye(nr),
    ZeroMatrix = zeros(nr, nc),
    OneMatrix = ones(nr, nc),
    MatAdd = Reduce(`+`, lapply(a, as_explicit)),
    MatMul = {
      sp <- .mterm_split(e)
      ms <- if (sp[[2L]]$type == "MatMul") sp[[2L]]$args else list(sp[[2L]])
      P <- Reduce(function(X, Y) .matmul(X, Y), lapply(ms, as_explicit))
      .mat_map(P, function(z) expand(.mul(list(sp[[1L]], z))))
    },
    MatPow = {
      k <- a[[2L]]
      if (k$type != "Integer") stop("as_explicit() needs an integer power", call. = FALSE)
      M <- as_explicit(a[[1L]])
      if (.int_val(k) < 0) M <- inverse(M)
      Reduce(.matmul, rep(list(M), abs(.int_val(k))))
    },
    Transpose = t(as_explicit(a[[1L]])),
    Inverse = inverse(as_explicit(a[[1L]])),
    Adjoint = .mat_map(t(as_explicit(a[[1L]])), function(z) fcall("conjugate", list(z))),
    HadamardProduct = Reduce(function(X, Y) .new_mat(Map(function(u, v) .mul(list(u, v)), X$data, Y$data), X$nrow, X$ncol),
                             lapply(a, as_explicit)),
    KroneckerProduct = Reduce(kronecker_product, lapply(a, as_explicit)),
    FunctionMatrix = build(function(i, j) apply_lambda(a[[3L]], .int(i - 1L), .int(j - 1L))),
    DiagMatrix = {
      v <- as_explicit(a[[1L]])
      build(function(i, j) if (i == j) v$data[[i]] else S_Zero)
    },
    DiagonalOf = {
      M <- as_explicit(a[[1L]])
      build(function(i, j) .mget(M, i, i))
    },
    DiagonalMatrix = {
      M <- as_explicit(a[[1L]])
      build(function(i, j) if (i == j) .mget(M, i, j) else S_Zero)
    },
    PermutationMatrix = {
      arr <- vapply(a[[2L]]$args, .int_val, 0)
      build(function(i, j) if (arr[i] == j - 1L) S_One else S_Zero)
    },
    MatrixPermute = {
      M <- as_explicit(a[[1L]])
      arr <- vapply(a[[2L]]$args, .int_val, 0) + 1L
      if (.int_val(a[[3L]]) == 0L) build(function(i, j) .mget(M, arr[i], j)) else build(function(i, j) .mget(M, i, arr[j]))
    },
    BlockMatrix = , BlockDiagMatrix = {
      R <- .block_rows(e)
      do.call(rbind, lapply(R, function(r) do.call(cbind, lapply(r, as_explicit))))
    },
    {
      h <- .me_explicit_hooks[[e$type]]
      if (!is.null(h)) return(h(e))
      h <- .me_entry_hooks[[e$type]]
      if (is.null(h)) stop("as_explicit() is not implemented for ", e$type, call. = FALSE)
      build(function(i, j) h(e, .int(i - 1L), .int(j - 1L)))
    })
}

# ---------------------------------------------------------------------------
# Arithmetic and printing

.matexpr_ops <- function(gen, e1, e2) {
  me <- function(z) if (is.numeric(z)) as_sym(z) else z
  e1 <- me(e1); e2 <- me(e2)
  # MatrixExpr.__add__ and friends: the canonical node, then doit()
  switch(gen,
    "+" = doit(MatAdd(e1, e2)),
    "-" = doit(MatAdd(e1, MatMul(S_NegativeOne, e2))),
    "*" = doit(MatMul(e1, e2)),
    "/" = if (!.is_matexpr(e2)) doit(MatMul(.div(S_One, as_sym(e2)), e1))
          else if (.is_matexpr(e1)) doit(MatMul(e1, doit(Inverse(e2))))
          else .perr("NotImplementedError", ""),
    "^" = MatPow(e1, e2),
    "==" = identical(e1$key, e2$key),
    "!=" = !identical(e1$key, e2$key),
    stop("unsupported matrix expression operator ", gen, call. = FALSE))
}

.mparen <- function(a) {
  s <- .str(a)
  if (a$type %in% c("MatAdd")) paste0("(", s, ")") else s
}
.mpost <- function(a) {
  s <- .str(a)
  if (a$type %in% c("MatrixSymbol", "Identity", "ZeroMatrix", "OneMatrix", "Transpose", "FunctionMatrix",
                    "DiagMatrix", "DiagonalOf", "Adjoint", "ImmutableMatrix", "GenericIdentity", "GenericZeroMatrix",
                    "ElementwiseApplyFunction", "MatrixUnit", "CompanionMatrix", "DFT", "IDFT", .me_factor_types)) s else paste0("(", s, ")")
}

register_str_printer("MatrixSymbol", function(e) e$name)
register_str_printer("Identity", function(e) "I")
# MatrixBase._format_str, how StrPrinter prints a matrix inside an
# expression: columns right justified, one row per line.
.mat_format_str <- function(nr, nc, cells) {
  if (nr == 0L || nc == 0L) return(sprintf("Matrix(%d, %d, [])", nr, nc))
  s <- matrix(vapply(cells, .str, ""), nr, nc, byrow = TRUE)
  w <- apply(nchar(s), 2L, max)
  rows <- vapply(seq_len(nr), function(i) paste0("[", paste(sprintf("%*s", w, s[i, ]), collapse = ", "), "]"), "")
  if (nr == 1L) paste0("Matrix([", rows, "])") else paste0("Matrix([\n", paste(rows, collapse = ",\n"), "])")
}
register_str_printer("ImmutableMatrix", function(e)
  .mat_format_str(.int_val(e$args[[1L]]), .int_val(e$args[[2L]]), e$args[[3L]]$args))
register_str_printer("ZeroMatrix", function(e) "0")
register_str_printer("OneMatrix", function(e) "1")
# StrPrinter._print_MatMul: a negative number coefficient prints as a sign,
# and each factor is parenthesized at the precedence of Mul
register_str_printer("MatMul", function(e) {
  args <- e$args
  isc <- !vapply(args, .is_matexpr, NA)
  c <- .mul(args[isc])
  sign <- ""
  if (!length(.free_sym_list(c)) && could_extract_minus_sign(c)) {
    sign <- "-"
    c <- .neg(c)
    args <- c(if (!is_one_num(c)) list(c), args[!isc])
  }
  paste0(sign, paste(vapply(args, .paren, "", PREC$Mul), collapse = "*"))
})
# Add.as_ordered_terms for matrix terms: by the sort keys of the matrix
# factors (class name, then the factor), so explicit matrices come first.
# the printed order of the terms: Expr.as_ordered_terms
.matadd_print_order <- function(ts) .me_ordered_terms(ts)
register_str_printer("MatAdd", function(e) {
  parts <- vapply(.matadd_print_order(e$args), .str, "")
  out <- parts[1L]
  for (p in parts[-1L]) out <- if (startsWith(p, "-")) paste0(out, " - ", substring(p, 2L)) else paste0(out, " + ", p)
  out
})
register_str_printer("MatPow", function(e) paste0(.mpost(e$args[[1L]]), "^", .paren(e$args[[2L]], PREC$Pow)))
register_str_printer("Inverse", function(e) paste0(.mpost(e$args[[1L]]), "^(-1)"))
register_str_printer("Transpose", function(e) paste0(.mpost(e$args[[1L]]), ".T"))
for (.t in c("Adjoint", "Trace", "Determinant", "Permanent", "HadamardProduct", "HadamardPower", "KroneckerProduct",
             "DiagMatrix", "DiagonalOf", "DiagonalMatrix", "DotProduct", "FunctionMatrix"))
  register_str_printer(.t, local({
    tt <- .t
    function(e) paste0(tt, "(", paste(vapply(e$args, .str, ""), collapse = ", "), ")")
  }))
rm(.t)
register_str_printer("MatrixElement", function(e) paste0(.mpost(e$args[[1L]]), "[", .str(e$args[[2L]]), ", ", .str(e$args[[3L]]), "]"))
# SymPy's StrPrinter means to print a 1x1 block matrix as its block but
# lacks the return, so it always prints the grid; this matches the output.
# SymPy's str of a matrix: one row inline, several rows one per line.
.block_str_all <- function(e) {
  R <- .block_rows(e)
  .mat_format_str(length(R), length(R[[1L]]), do.call(c, R))
}

.mslice_str <- function(s, dim, pr) {
  x <- s$args
  if (is_one_num(x[[3L]])) x <- x[1:2]
  out <- vapply(x, pr, "")
  if (is_zero_num(x[[1L]])) out[1L] <- ""
  if (identical(x[[2L]]$key, dim$key)) out[2L] <- ""
  paste(out, collapse = ":")
}
register_str_printer("MatrixSlice", function(e) {
  sh <- .mshape(e$args[[1L]])
  paste0(.mpost(e$args[[1L]]), "[", .mslice_str(e$args[[2L]], sh[[1L]], .str), ", ", .mslice_str(e$args[[3L]], sh[[2L]], .str), "]")
})
register_str_printer("BlockMatrix", .block_str_all)
register_str_printer("BlockDiagMatrix", .block_str_all)
register_str_printer("PermutationMatrix", function(e) paste0("PermutationMatrix(", .str(e$args[[2L]]), ")"))
register_str_printer("MatrixPermute", function(e) paste0("MatrixPermute(", .str(e$args[[1L]]), ", ", .str(e$args[[2L]]), ", ", .str(e$args[[3L]]), ")"))

.lmpost <- function(a) {
  s <- .latex(a)
  if (a$type %in% c("MatrixSymbol", "Identity", "ZeroMatrix", "OneMatrix")) s else paste0("\\left(", s, "\\right)")
}
register_latex_printer("MatrixSymbol", function(e) e$name)
register_latex_printer("Identity", function(e) "\\mathbb{I}")
register_latex_printer("ZeroMatrix", function(e) "0")
register_latex_printer("OneMatrix", function(e) "1")
register_latex_printer("MatMul", function(e) {
  sp <- .mterm_split(e)
  ms <- if (sp[[2L]]$type == "MatMul") sp[[2L]]$args else list(sp[[2L]])
  body <- paste(vapply(ms, function(a) if (a$type == "MatAdd") paste0("\\left(", .latex(a), "\\right)") else .latex(a), ""), collapse = " ")
  c <- sp[[1L]]
  if (identical(c$key, S_One$key)) body else if (identical(c$key, S_NegativeOne$key)) paste0("- ", body) else paste(.latex(c), body)
})
register_latex_printer("MatAdd", function(e) {
  parts <- vapply(e$args, .latex, "")
  out <- parts[1L]
  for (p in parts[-1L]) out <- if (startsWith(p, "- ")) paste0(out, " - ", substring(p, 3L)) else paste0(out, " + ", p)
  out
})
register_latex_printer("MatPow", function(e) paste0(.lmpost(e$args[[1L]]), "^{", .latex(e$args[[2L]]), "}"))
register_latex_printer("Inverse", function(e) paste0(.lmpost(e$args[[1L]]), "^{-1}"))
register_latex_printer("Transpose", function(e) paste0(.lmpost(e$args[[1L]]), "^{T}"))
register_latex_printer("Adjoint", function(e) paste0(.lmpost(e$args[[1L]]), "^{\\dagger}"))
register_latex_printer("Trace", function(e) paste0("\\operatorname{tr}\\left(", .latex(e$args[[1L]]), " \\right)"))
register_latex_printer("Determinant", function(e) paste0("\\left|{", .latex(e$args[[1L]]), "}\\right|"))
register_latex_printer("HadamardProduct", function(e) paste(vapply(e$args, .latex, ""), collapse = " \\circ "))
register_latex_printer("KroneckerProduct", function(e) paste(vapply(e$args, .latex, ""), collapse = " \\otimes "))
register_latex_printer("MatrixSlice", function(e) {
  sh <- .mshape(e$args[[1L]])
  paste0(.latex(e$args[[1L]]), "\\left[", .mslice_str(e$args[[2L]], sh[[1L]], .latex), ", ", .mslice_str(e$args[[3L]], sh[[2L]], .latex), "\\right]")
})
register_latex_printer("MatrixElement", function(e) paste0("{", .latex(e$args[[1L]]), "}_{", .latex(e$args[[2L]]), ",", .latex(e$args[[3L]]), "}"))

register_rebuilder("MatAdd", function(e, args) do.call(MatAdd, args))
register_rebuilder("MatMul", function(e, args) do.call(MatMul, args))
# MatAdd.doit and MatMul.doit: canonicalize after the arguments' doit()
register_doit("MatAdd", function(e) do.call(.madd, e$args))
register_doit("MatMul", function(e) do.call(.mm, e$args))
register_rebuilder("MatPow", function(e, args) MatPow(args[[1L]], args[[2L]]))
register_rebuilder("MatrixElement", function(e, args) MatrixElement(args[[1L]], args[[2L]], args[[3L]]))

# ---------------------------------------------------------------------------
# Explicit matrix constructors

#' Matrix constructors
#'
#' `sym_diag(...)` is SymPy's `diag()` (block diagonal from scalars and
#' matrices; renamed so it does not mask [base::diag()]). `banded(...)`
#' builds a banded matrix from a named list of diagonals, `list("0" = 1,
#' "-1" = c(2, 3))`, with the size given as `banded(n, d)` or `banded(r, c,
#' d)` or inferred. `rot_givens(i, j, theta, dim)` is the Givens rotation in
#' the `(i, j)` plane (zero based). `sym_trace(M)` is the trace (renamed so
#' it does not mask [base::trace()]). `SparseMatrix(r, c, entries)` takes
#' `list(list(i, j, value), ...)` with zero based indices.
#'
#' @param ... Diagonal entries or blocks; for `banded()`, the size and the
#'   diagonals.
#' @param i,j Plane indices.
#' @param theta The angle.
#' @param dim The dimension.
#' @param M A matrix or matrix expression.
#' @param r,c Dimensions.
#' @param entries Nonzero entries.
#' @return A matrix.
#' @examples
#' sym_diag(1, 2, 3)
#' banded(3, list("0" = 1, "1" = 2))
#' rot_givens(0, 1, pi/3)
#' @name matrix-builders
NULL

#' @rdname matrix-builders
#' @export
sym_diag <- function(...) diag_matrix(...)

#' @rdname matrix-builders
#' @export
sym_trace <- function(M) {
  if (inherits(M, "symmatrix")) return(.trace_explicit(M))
  doit(Trace(M))
}

#' @rdname matrix-builders
#' @export
banded <- function(...) {
  args <- list(...)
  d <- args[[length(args)]]
  dims <- unlist(args[-length(args)])
  ks <- as.integer(names(d))
  isseq <- function(v) !is_sym(v) && (length(v) > 1L || is.list(v))
  vals <- lapply(d, function(v) if (isseq(v)) lapply(v, as_sym) else as_sym(v))
  islist <- vapply(d, isseq, TRUE)
  if (length(dims) == 2L) { nr <- dims[1L]; nc <- dims[2L] }
  else if (length(dims) == 1L) { nr <- nc <- dims[1L] }
  else {
    if (!any(islist)) stop("unable to infer the size of the banded matrix", call. = FALSE)
    nr <- max(vapply(which(islist), function(i) length(vals[[i]]) + max(0L, -ks[i]), 0))
    nc <- max(vapply(which(islist), function(i) length(vals[[i]]) + max(0L, ks[i]), 0))
  }
  data <- rep(list(S_Zero), nr * nc)
  for (q in seq_along(ks)) {
    k <- ks[q]
    i0 <- max(0L, -k); j0 <- max(0L, k)
    len <- min(nr - i0, nc - j0)
    if (len <= 0L) next
    for (t in seq_len(len)) {
      v <- if (islist[q]) { if (t > length(vals[[q]])) next; vals[[q]][[t]] } else vals[[q]]
      data[[(i0 + t - 1L) * nc + j0 + t]] <- v
    }
  }
  .new_mat(data, nr, nc)
}

#' @rdname matrix-builders
#' @export
rot_givens <- function(i, j, theta, dim = 3) {
  theta <- as_sym(theta)
  M <- .mrows(eye(dim))
  c <- fcall("cos", list(theta)); s <- fcall("sin", list(theta))
  a <- i + 1L; b <- j + 1L
  M[[a]][[a]] <- c; M[[b]][[b]] <- c
  M[[a]][[b]] <- s; M[[b]][[a]] <- .neg(s)
  .from_rows(M)
}

#' @rdname matrix-builders
#' @export
SparseMatrix <- function(r, c = NULL, entries = NULL) {
  if (is.null(entries)) return(Matrix(r, c))
  data <- rep(list(S_Zero), r * c)
  for (en in entries) data[[en[[1L]] * c + en[[2L]] + 1L]] <- as_sym(en[[3L]])
  .new_mat(data, r, c)
}
#' @rdname matrix-builders
#' @export
MutableSparseMatrix <- SparseMatrix
#' @rdname matrix-builders
#' @export
ImmutableSparseMatrix <- SparseMatrix
