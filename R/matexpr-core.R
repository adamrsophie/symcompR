# Further matrix expression classes, after sympy/matrices/expressions:
# applyfunc.py (ElementwiseApplyFunction), hadamard.py, kronecker.py,
# special.py (GenericIdentity, GenericZeroMatrix, MatrixUnit), companion.py,
# factorizations.py, fourier.py and sets.py.

# Hooks that matrices-expr.R reads for the classes defined in these files:
# the shape, MatrixExpr._entry, _eval_transpose and as_explicit.
.me_shape_hooks <- new.env(parent = emptyenv())
.me_entry_hooks <- new.env(parent = emptyenv())
.me_transpose_hooks <- new.env(parent = emptyenv())
.me_explicit_hooks <- new.env(parent = emptyenv())
.me_adjoint_hooks <- new.env(parent = emptyenv())
.me_conjugate_hooks <- new.env(parent = emptyenv())

register_rebuilder("MatrixSymbol", function(e, args) MatrixSymbol(e$name, args[[1L]], args[[2L]]))

# Expr.as_ordered_terms with the default order: the terms sorted by their
# commutative monomials (higher degree first), then by the sort keys of their
# noncommutative factors, then by their numeric coefficients
.me_ordered_terms <- function(ts) {
  if (length(ts) < 2L) return(ts)
  info <- lapply(ts, function(t) {
    coeff <- S_One
    rest <- t
    if (t$type == "MatMul" && is_number_type(t$args[[1L]])) {
      coeff <- t$args[[1L]]
      r <- t$args[-1L]
      rest <- if (length(r) == 1L) r[[1L]] else new_sym("MatMul", r)
    }
    facs <- if (rest$type %in% c("MatMul", "Mul")) rest$args else list(rest)
    cp <- list(); nc <- list()
    for (f in facs) {
      if (.is_matexpr(f) || .is_arrexpr(f)) nc[[length(nc) + 1L]] <- f
      else if (.expr_is_number(f)) coeff <- .mul(list(coeff, f))
      else {
        be <- .base_exp(f)
        k <- be[[1L]]$key
        cp[[k]] <- list(be[[1L]], if (is.null(cp[[k]])) be[[2L]] else .add(list(cp[[k]][[2L]], be[[2L]])))
      }
    }
    list(coeff = coeff, cp = cp, nc = nc)
  })
  gens <- list()
  for (z in info) for (g in z$cp) gens[[g[[1L]]$key]] <- g[[1L]]
  gens <- .key_sort(lapply(gens, .sort_key), unname(gens))
  num <- function(v) if (is_exact_number(v) || v$type == "Float") num_value_dbl(v) else 0
  keys <- lapply(info, function(z) {
    monom <- vapply(gens, function(g) { h <- z$cp[[g$key]]; if (is.null(h)) 0 else -num(h[[2L]]) }, 0)
    list(as.list(monom), lapply(z$nc, .sort_key), list(num(z$coeff)))
  })
  .key_sort(keys, ts)
}

# DiagMatrix.doit: a diagonal argument (1x1, an identity) is the result
.me_known_diagonal <- function(v) {
  if (!.is_matexpr(v) && !inherits(v, "symmatrix")) return(FALSE)
  sh <- .mshape(.as_mexpr(v))
  if (is_one_num(sh[[1L]]) && is_one_num(sh[[2L]])) return(TRUE)
  if (inherits(v, "symmatrix")) v <- .as_mexpr(v)
  if (v$type %in% c("Identity", "DiagMatrix", "DiagonalMatrix")) return(TRUE)
  if (v$type == "ImmutableMatrix") {
    M <- as_explicit(v)
    if (M$nrow != M$ncol) return(FALSE)
    for (i in seq_len(M$nrow)) for (j in seq_len(M$ncol)) if (i != j && !is_zero_num(.mget(M, i, j))) return(FALSE)
    return(TRUE)
  }
  FALSE
}
register_doit("DiagMatrix", function(e) {
  v <- e$args[[1L]]
  if (.me_known_diagonal(v)) return(.unwrap_explicit(v))
  if (.is_explicit(v)) {
    M <- as_explicit(v)
    n <- max(M$nrow, M$ncol)
    return(.new_mat(lapply(seq_len(n * n) - 1L, function(k) if (k %/% n == k %% n) M$data[[k %/% n + 1L]] else S_Zero), n, n))
  }
  if (v$type == "MatMul") {
    mats <- Filter(.is_matexpr, v$args)
    sc <- Filter(Negate(.is_matexpr), v$args)
    if (length(sc)) return(.mul(sc) * doit(DiagMatrix(doit(do.call(MatMul, mats)))))
  }
  if (v$type == "Transpose") v <- v$args[[1L]]
  DiagMatrix(v)
})

# MatrixExpr._check_dim
.me_check_dim <- function(d) {
  ok <- if (d$type == "Float") FALSE else {
    f <- .facts(d)
    if (identical(f[["integer"]], FALSE) || identical(f[["nonnegative"]], FALSE)) FALSE else TRUE
  }
  if (!ok) .perr("ValueError", paste0("The dimension specification ", sstr(d), " should be a nonnegative integer."))
  invisible(d)
}

.me_is_matrix <- function(e) .is_matexpr(e) || inherits(e, "symmatrix")
# The constructor postprocessors of Mul and Add with matrix arguments: a
# MatMul or MatAdd, with the scalars absorbed into an explicit matrix
.me_postprocess <- function(cls, args) {
  args <- lapply(args, function(a) if (.me_is_matrix(a)) a else as_sym(a))
  ident <- if (cls == "Mul") S_One else S_Zero
  args <- Filter(function(a) !(is_sym(a) && identical(a$key, ident$key)), args)
  ism <- vapply(args, .me_is_matrix, NA)
  mats <- args[ism]
  non <- args[!ism]
  if (cls == "Add") {
    m <- .unwrap_explicit(do.call(.madd, mats))
    if (!length(non)) return(m)
    # Add(scalar, matrix) stays an Add, so that algorithms can put matrices
    # in place of noncommutative symbols
    s <- .add(non)
    if (is_zero_num(s)) return(m)
    return(new_sym("Add", c(if (s$type == "Add") s$args else list(s), list(.as_mexpr(m)))))
  }
  if (length(non)) for (i in seq_along(mats)) if (.is_explicit(mats[[i]])) {
    mats[[i]] <- .mat_map(as_explicit(mats[[i]]), function(z) .mul(c(list(z), non)))
    non <- list()
    break
  }
  .unwrap_explicit(do.call(.mm, c(if (length(non)) list(.mul(non)), mats)))
}
.me_shape_str <- function(v) paste0("(", paste(vapply(v, .str, ""), collapse = ", "), if (length(v) == 1L) ",", ")")

# _shape.validate_matadd_integer
.me_validate_matadd <- function(args) {
  sh <- lapply(args, .mshape)
  for (k in 1:2) {
    d <- lapply(sh, `[[`, k)
    ints <- unique(vapply(Filter(function(z) z$type == "Integer", d), .k, ""))
    if (length(ints) > 1L) .perr("ShapeError", paste0("Matrices have mismatching shape: ", .me_shape_str(d)))
  }
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# GenericIdentity and GenericZeroMatrix, the identities of MatMul and MatAdd

#' Generic identity and zero matrices
#'
#' `GenericIdentity()` and `GenericZeroMatrix()` are the matrix identity and
#' zero of unspecified shape that `MatMul()` and `MatAdd()` of no arguments
#' give. They print as `I` and `0`; asking for their shape is an error.
#'
#' @return A `sym` matrix expression.
#' @examples
#' GenericIdentity()
#' MatMul()
#' MatAdd()
#' @export
GenericIdentity <- function() new_sym("GenericIdentity", list())

#' @rdname GenericIdentity
#' @export
GenericZeroMatrix <- function() new_sym("GenericZeroMatrix", list())

register_str_printer("GenericIdentity", function(e) "I")
register_str_printer("GenericZeroMatrix", function(e) "0")
.me_shape_hooks$GenericIdentity <- function(e) .perr("TypeError", "GenericIdentity does not have a specified shape")
.me_shape_hooks$GenericZeroMatrix <- function(e) .perr("TypeError", "GenericZeroMatrix does not have a specified shape")

# ---------------------------------------------------------------------------
# MatrixUnit

#' Matrix unit
#'
#' `MatrixUnit(rows, cols, i, j)` is the `rows` by `cols` matrix with a one
#' at position `(i, j)` (zero based) and zeros elsewhere: the derivative of
#' the element `A[i, j]` with respect to `A`.
#'
#' @param rows,cols The shape.
#' @param i,j The position of the unit entry.
#' @return A `sym` matrix expression.
#' @examples
#' k <- Symbol("k")
#' MatrixUnit(k, k, 0, 1)
#' as_explicit(MatrixUnit(2, 2, 0, 1))
#' @export
MatrixUnit <- function(rows, cols, i, j) new_sym("MatrixUnit", lapply(list(rows, cols, i, j), as_sym))

.me_shape_hooks$MatrixUnit <- function(e) e$args[1:2]
.me_entry_hooks$MatrixUnit <- function(e, i, j, ...) {
  a <- e$args
  .mul(list(KroneckerDelta(i, a[[3L]], Tuple(S_Zero, .sub(a[[1L]], S_One))),
            KroneckerDelta(j, a[[4L]], Tuple(S_Zero, .sub(a[[2L]], S_One)))))
}
.me_transpose_hooks$MatrixUnit <- function(e) MatrixUnit(e$args[[2L]], e$args[[1L]], e$args[[4L]], e$args[[3L]])

# ---------------------------------------------------------------------------
# ElementwiseApplyFunction

#' Elementwise application of a function to a matrix expression
#'
#' `ElementwiseApplyFunction(fun, expr)` applies a scalar function to every
#' entry of the matrix `expr`, as SymPy's `expr.applyfunc(fun)` does for a
#' matrix expression. `fun` is an R function of one argument or a
#' [Lambda()]; it is kept as a `Lambda`. The result prints as
#' `Lambda(_d, f(_d)).(expr)`.
#'
#' @param fun A function of one argument or a `Lambda`.
#' @param expr A matrix expression or an explicit matrix.
#' @return A `sym` matrix expression.
#' @examples
#' k <- Symbol("k")
#' X <- MatrixSymbol("X", k, k)
#' ElementwiseApplyFunction(sin, X)
#' @export
ElementwiseApplyFunction <- function(fun, expr) {
  if (is.numeric(expr)) expr <- as_sym(expr)
  if (!.me_is_matrix(expr)) .perr("ValueError", paste0(.str(expr), " must be a matrix instance."))
  expr <- .as_mexpr(expr)
  sh <- .mshape(expr)
  if (is_one_num(sh[[1L]]) && is_one_num(sh[[2L]])) {
    ret <- tryCatch(.me_lift(if (is.function(fun)) fun(expr) else apply_lambda(fun, expr)), error = function(e) NULL)
    if (.is_matexpr(ret)) return(ret)
  }
  lam <- .ae_as_lambda(fun)
  if (length(lam$args[[1L]]$args) != 1L)
    .perr("ValueError", paste0(.str(lam), " should be able to accept 1 arguments."))
  new_sym("ElementwiseApplyFunction", list(lam, expr))
}

register_str_printer("ElementwiseApplyFunction", function(e) paste0(.str(e$args[[1L]]), ".(", .str(e$args[[2L]]), ")"))
register_rebuilder("ElementwiseApplyFunction", function(e, args) ElementwiseApplyFunction(args[[1L]], args[[2L]]))
.me_shape_hooks$ElementwiseApplyFunction <- function(e) .mshape(e$args[[2L]])
.me_entry_hooks$ElementwiseApplyFunction <- function(e, i, j, ...) apply_lambda(e$args[[1L]], .me_entry(e$args[[2L]], i, j, ...))
.me_transpose_hooks$ElementwiseApplyFunction <- function(e) ElementwiseApplyFunction(e$args[[1L]], doit(Transpose(e$args[[2L]])))

.me_lambda_is_identity <- function(lam) {
  v <- lam$args[[1L]]$args
  length(v) == 1L && identical(v[[1L]]$key, lam$args[[2L]]$key)
}

register_doit("ElementwiseApplyFunction", function(e, deep = TRUE) {
  expr <- e$args[[2L]]
  lam <- e$args[[1L]]
  if (.me_lambda_is_identity(lam)) return(.unwrap_explicit(expr))
  if (.is_explicit(expr)) return(.mat_map(as_explicit(expr), function(z) apply_lambda(lam, z)))
  if (expr$type == "ElementwiseApplyFunction") {
    inner <- expr$args[[1L]]
    return(doit(ElementwiseApplyFunction(function(x) apply_lambda(lam, apply_lambda(inner, x)), expr$args[[2L]]), deep = deep))
  }
  e
})

# ---------------------------------------------------------------------------
# Hadamard products and powers

.me_hadamard_new <- function(args, evaluate = FALSE) {
  args <- lapply(args, function(a) if (is.numeric(a)) as_sym(a) else .as_mexpr(a))
  if (!length(args)) .perr("ValueError", "HadamardProduct needs at least one argument")
  if (!all(vapply(args, .is_matexpr, NA))) .perr("TypeError", "Mix of Matrix and Scalar symbols")
  .me_validate_matadd(args)
  e <- new_sym("HadamardProduct", args)
  if (isTRUE(evaluate)) e <- doit(e, deep = FALSE)
  e
}

# hadamard.canonicalize
.me_hadamard_canon <- function(x) {
  if (x$type != "HadamardProduct") return(x)
  repeat {
    fl <- do.call(c, lapply(x$args, function(a) if (a$type == "HadamardProduct") a$args else list(a)))
    if (length(fl) == length(x$args) && all(vapply(x$args, function(a) a$type != "HadamardProduct", NA))) break
    x <- new_sym("HadamardProduct", fl)
  }
  ids <- vapply(x$args, function(a) a$type == "OneMatrix", NA)
  if (any(ids)) x <- new_sym("HadamardProduct", if (all(ids)) x$args[1L] else x$args[!ids])
  if (any(vapply(x$args, function(a) a$type == "ZeroMatrix", NA))) return(do.call(ZeroMatrix, .mshape(x)))
  keys <- vapply(x$args, .k, "")
  u <- unique(keys)
  newargs <- lapply(u, function(k) {
    n <- sum(keys == k)
    b <- x$args[[match(k, keys)]]
    if (n == 1L) b else HadamardPower(b, n)
  })
  newargs <- newargs[.ae_default_sort_order(newargs)]
  if (length(newargs) == 1L) return(newargs[[1L]])
  new_sym("HadamardProduct", newargs)
}

register_doit("HadamardProduct", function(e, deep = TRUE) {
  args <- e$args
  ex <- vapply(args, .is_explicit, NA)
  if (any(ex)) {
    em <- lapply(args[ex], as_explicit)
    M <- .new_mat(lapply(seq_along(em[[1L]]$data), function(k) .mul(lapply(em, function(m) m$data[[k]]))),
                  em[[1L]]$nrow, em[[1L]]$ncol)
    e <- new_sym("HadamardProduct", c(list(.as_block(M)), args[!ex]))
  }
  .unwrap_explicit(.me_hadamard_canon(e))
})
register_rebuilder("HadamardProduct", function(e, args) .me_hadamard_new(args))
.me_entry_hooks$HadamardProduct <- function(e, i, j, ...) .mul(lapply(e$args, function(a) .me_entry(a, i, j, ...)))

#' Hadamard (elementwise) products and powers
#'
#' `hadamard_product(...)` multiplies matrices elementwise and puts the
#' result in canonical form; `hadamard_power(base, exp)` raises a matrix
#' elementwise to a power. With explicit matrices only, the result is an
#' explicit matrix.
#'
#' @param ... Matrices or matrix expressions.
#' @param base A matrix expression or a scalar.
#' @param exp The exponent, a scalar.
#' @return A matrix expression or an explicit matrix.
#' @examples
#' k <- Symbol("k")
#' A <- MatrixSymbol("A", k, k); B <- MatrixSymbol("B", k, k)
#' hadamard_product(B, A, A)
#' hadamard_power(A, 2)
#' @export
hadamard_product <- function(...) {
  ms <- list(...)
  if (!length(ms)) .perr("TypeError", "Empty Hadamard product is undefined")
  if (length(ms) == 1L) return(ms[[1L]])
  doit(.me_hadamard_new(ms))
}

#' @rdname hadamard_product
#' @export
hadamard_power <- function(base, exp) {
  base <- if (is.numeric(base)) as_sym(base) else base
  exp <- if (is.numeric(exp)) as_sym(exp) else exp
  if (is_sym(exp) && is_one_num(exp)) return(base)
  if (!.me_is_matrix(base)) return(.pow(base, exp))
  if (.me_is_matrix(exp)) .perr("ValueError", "cannot raise expression to a matrix")
  HadamardPower(base, exp)
}

.me_hadamard_power_new <- function(base, exp) {
  base <- if (is.numeric(base)) as_sym(base) else .as_mexpr(base)
  exp <- if (is.numeric(exp)) as_sym(exp) else .as_mexpr(exp)
  if (!.is_matexpr(base) && !.is_matexpr(exp)) return(.pow(base, exp))
  if (.is_matexpr(base) && .is_matexpr(exp)) .me_validate_matadd(list(base, exp))
  new_sym("HadamardPower", list(base, exp))
}
.me_shape_hooks$HadamardPower <- function(e) .mshape(if (.is_matexpr(e$args[[1L]])) e$args[[1L]] else e$args[[2L]])
.me_entry_hooks$HadamardPower <- function(e, i, j, ...) {
  ent <- function(z) if (.is_matexpr(z)) .me_entry(z, i, j, ...) else z
  .pow(ent(e$args[[1L]]), ent(e$args[[2L]]))
}
.me_transpose_hooks$HadamardPower <- function(e) HadamardPower(.mtranspose(e$args[[1L]]), e$args[[2L]])
register_rebuilder("HadamardPower", function(e, args) .me_hadamard_power_new(args[[1L]], args[[2L]]))

# ---------------------------------------------------------------------------
# Kronecker products

.me_is_identity <- function(a) {
  if (inherits(a, "symmatrix")) a <- .as_block(a)
  if (a$type %in% c("Identity", "GenericIdentity")) return(TRUE)
  if (a$type == "ImmutableMatrix") {
    M <- as_explicit(a)
    if (M$nrow != M$ncol) return(FALSE)
    for (i in seq_len(M$nrow)) for (j in seq_len(M$ncol))
      if (!identical(.mget(M, i, j)$key, (if (i == j) S_One else S_Zero)$key)) return(FALSE)
    return(TRUE)
  }
  FALSE
}

.me_kronecker_new <- function(args) {
  args <- lapply(args, function(a) if (is.numeric(a)) as_sym(a) else a)
  if (length(args) && all(vapply(args, function(a) .me_is_matrix(a) && .me_is_identity(a), NA))) {
    n <- .mul(lapply(args, function(a) .mshape(.as_mexpr(a))[[1L]]))
    if (all(vapply(args, .is_explicit, NA))) return(eye(.int_val(n)))
    return(Identity(n))
  }
  if (!all(vapply(args, .me_is_matrix, NA))) .perr("TypeError", "Mix of Matrix and Scalar symbols")
  new_sym("KroneckerProduct", lapply(args, .as_mexpr))
}

# Expr.args_cnc of a factor of a Kronecker product
.me_cnc <- function(a) {
  if (a$type == "MatMul") {
    sc <- Filter(function(z) !.is_matexpr(z), a$args)
    mt <- Filter(.is_matexpr, a$args)
    return(list(sc, if (length(mt) == 1L) mt[[1L]] else new_sym("MatMul", mt)))
  }
  list(list(), a)
}

# one step of kronecker.canonicalize: the first rule that changes x
.me_kron_step <- function(x) {
  a <- x$args
  if (length(a) == 1L) return(a[[1L]])
  if (all(vapply(a, .is_explicit, NA))) return(Reduce(.kron2, lapply(a, as_explicit)))
  if (any(vapply(a, function(z) z$type == "KroneckerProduct", NA)))
    return(new_sym("KroneckerProduct", do.call(c, lapply(a, function(z) if (z$type == "KroneckerProduct") z$args else list(z)))))
  cs <- lapply(a, .me_cnc)
  c <- .mul(do.call(c, lapply(cs, `[[`, 1L)))
  if (!is_one_num(c)) return(.ae_mmul(c, new_sym("KroneckerProduct", lapply(cs, `[[`, 2L))))
  single <- vapply(a, function(z) {
    sh <- .mshape(z)
    is_one_num(sh[[1L]]) && is_one_num(sh[[2L]]) && .me_commutative_entries(z)
  }, NA)
  if (sum(single) > 1L) return(new_sym("KroneckerProduct", c(a[!single], list(doit(do.call(MatMul, a[single]))))))
  x
}
.me_commutative_entries <- function(z) {
  if (z$type == "ImmutableMatrix") return(all(vapply(as_explicit(z)$data, function(v) !identical(.facts(v)[["commutative"]], FALSE), NA)))
  TRUE
}
.me_kron_canon <- function(x) {
  while (is_sym(x) && x$type == "KroneckerProduct") {
    y <- .me_kron_step(x)
    if (is_sym(y) && identical(y$key, x$key)) break
    x <- y
  }
  .unwrap_explicit(x)
}

register_doit("KroneckerProduct", function(e, deep = TRUE) .me_kron_canon(.me_kronecker_new(e$args)))
register_rebuilder("KroneckerProduct", function(e, args) .me_kronecker_new(args))
.me_entry_hooks$KroneckerProduct <- function(e, i, j, ...) {
  res <- S_One
  for (m in rev(e$args)) {
    sh <- .mshape(m)
    mi <- fcall("Mod", list(i, sh[[1L]])); mj <- fcall("Mod", list(j, sh[[2L]]))
    i <- fcall("floor", list(.div(i, sh[[1L]]))); j <- fcall("floor", list(.div(j, sh[[2L]])))
    res <- .mul(list(res, .me_getitem(m, mi, mj)))
  }
  res
}
.me_transpose_hooks$KroneckerProduct <- function(e) doit(.me_kronecker_new(lapply(e$args, .mtranspose)))
.me_adjoint_hooks$KroneckerProduct <- function(e) doit(.me_kronecker_new(lapply(e$args, .madjoint)))
.me_conjugate_hooks$KroneckerProduct <- function(e) doit(.me_kronecker_new(lapply(e$args, .me_conjugate)))
# _eval_adjoint and _eval_conjugate of the special matrices
.me_adjoint_hooks$ZeroMatrix <- function(e) ZeroMatrix(e$args[[2L]], e$args[[1L]])
.me_adjoint_hooks$OneMatrix <- function(e) OneMatrix(e$args[[2L]], e$args[[1L]])
.me_adjoint_hooks$Identity <- function(e) e
.me_adjoint_hooks$GenericIdentity <- function(e) e
.me_adjoint_hooks$MatrixUnit <- function(e) .mtranspose(e)
.me_conjugate_hooks$MatrixUnit <- function(e) e

# MatrixExpr.as_real_imag, for re() and im() of a matrix
.me_real_imag <- function(M) {
  if (inherits(M, "symmatrix"))
    return(list(.mat_map(M, function(z) fcall("re", list(z))), .mat_map(M, function(z) fcall("im", list(z)))))
  if (M$type == "ImmutableMatrix") return(.me_real_imag(as_explicit(M)))
  cj <- .me_conjugate(M)
  list(S_Half * (M + cj), (M - cj) / (2 * S_I))
}

# KroneckerProduct._eval_determinant
.me_kron_det <- function(A) {
  sq <- vapply(A$args, function(a) { s <- .mshape(a); identical(s[[1L]]$key, s[[2L]]$key) }, NA)
  if (!all(sq)) return(Determinant(A))
  m <- .mshape(A)[[1L]]
  .mul(lapply(A$args, function(a) .pow(.det_doit(a), .div(m, .mshape(a)[[1L]]))))
}

# conjugate() of a matrix: the _eval_conjugate methods
.me_conjugate <- function(e) {
  if (inherits(e, "symmatrix")) return(.mat_map(e, function(z) fcall("conjugate", list(z))))
  if (!.is_matexpr(e)) return(fcall("conjugate", list(as_sym(e))))
  a <- e$args
  switch(e$type,
    ImmutableMatrix = .me_conjugate(as_explicit(e)),
    ZeroMatrix = , Identity = , OneMatrix = , GenericIdentity = , GenericZeroMatrix = e,
    Adjoint = .mtranspose(a[[1L]]),
    Transpose = .madjoint(a[[1L]]),
    Inverse = Inverse(.as_mexpr(.me_conjugate(a[[1L]]))),
    MatPow = if (isTRUE(.facts(a[[2L]])[["integer"]]) && isTRUE(.facts(a[[2L]])[["nonnegative"]]))
      MatPow(.as_mexpr(.me_conjugate(a[[1L]])), a[[2L]]) else fnode("conjugate", list(e)),
    MatMul = do.call(MatMul, lapply(a, function(z) .as_mexpr(.me_conjugate(z)))),
    MatAdd = do.call(MatAdd, lapply(a, function(z) .as_mexpr(.me_conjugate(z)))),
    { h <- .me_conjugate_hooks[[e$type]]; if (is.null(h)) Adjoint(Transpose(e)) else h(e) })
}

#' Kronecker products of matrices
#'
#' `kronecker_product(...)` is the Kronecker (tensor) product of matrices or
#' matrix expressions, in canonical form: scalar factors come out in front,
#' nested products flatten and explicit matrices multiply out.
#' `combine_kronecker(expr)` merges sums, products and powers of Kronecker
#' products inside a matrix expression.
#'
#' @param ... Matrices or matrix expressions.
#' @param expr A matrix expression.
#' @return A matrix expression or an explicit matrix.
#' @examples
#' n <- Symbol("n"); m <- Symbol("m")
#' A <- MatrixSymbol("A", n, m); B <- MatrixSymbol("B", n, m)
#' kronecker_product(2*A, B)
#' combine_kronecker(kronecker_product(A, B) + kronecker_product(B, B))
#' @export
kronecker_product <- function(...) {
  ms <- list(...)
  if (!length(ms)) .perr("TypeError", "Empty Kronecker product is undefined")
  if (length(ms) == 1L) return(ms[[1L]])
  if (all(vapply(ms, function(m) inherits(m, "symmatrix") || is.matrix(m) || is.list(m) && !is_sym(m), NA)))
    return(Reduce(.kron2, lapply(ms, Matrix)))
  doit(.me_kronecker_new(ms))
}

# the Kronecker product of two explicit matrices
.kron2 <- function(a, b) {
  Matrix(a$nrow * b$nrow, a$ncol * b$ncol, function(i, j) {
    ai <- (i - 1L) %/% b$nrow + 1L
    bi <- (i - 1L) %% b$nrow + 1L
    aj <- (j - 1L) %/% b$ncol + 1L
    bj <- (j - 1L) %% b$ncol + 1L
    .mul(list(.mget(a, ai, aj), .mget(b, bi, bj)))
  })
}

.me_coeff_kron <- function(e) {
  if (e$type == "KroneckerProduct") return(list(S_One, e))
  if (e$type == "MatMul") {
    ms <- Filter(.is_matexpr, e$args)
    if (length(ms) == 1L && ms[[1L]]$type == "KroneckerProduct") return(list(.mul(Filter(Negate(.is_matexpr), e$args)), ms[[1L]]))
  }
  NULL
}
.me_kron_dims_key <- function(e) {
  ck <- .me_coeff_kron(e)
  if (is.null(ck)) return("0")
  paste(vapply(ck[[2L]]$args, function(a) paste(vapply(.mshape(a), .k, ""), collapse = ","), ""), collapse = ";")
}
.me_structurally_equal <- function(a, b) {
  if (b$type != "KroneckerProduct" || length(a$args) != length(b$args)) return(FALSE)
  sk <- function(e) paste(vapply(.mshape(e), .k, ""), collapse = ",")
  sk(a) == sk(b) && all(mapply(function(u, v) sk(u) == sk(v), a$args, b$args))
}
.me_kron_merge_add <- function(p1, p2) {
  k1 <- p1[[2L]]; k2 <- p2[[2L]]
  if (!.me_structurally_equal(k1, k2)) return(NULL)
  d <- which(!mapply(function(u, v) identical(u$key, v$key), k1$args, k2$args))
  if (!length(d)) return(list(.add(list(p1[[1L]], p2[[1L]])), k1))
  if (length(d) == 1L) {
    na <- k1$args
    na[[d]] <- .ae_madd(.ae_mmul(p1[[1L]], k1$args[[d]]), .ae_mmul(p2[[1L]], k2$args[[d]]))
    return(list(S_One, .me_kronecker_new(na)))
  }
  NULL
}
.me_kron_mat_add <- function(e) {
  keys <- vapply(e$args, .me_kron_dims_key, "")
  nonkrons <- e$args[keys == "0"]
  groups <- unique(keys[keys != "0"])
  if (!length(groups)) return(e)
  krons <- list()
  for (g in groups) {
    merged <- list()
    for (t in e$args[keys == g]) {
      pair <- .me_coeff_kron(t)
      done <- FALSE
      for (i in seq_along(merged)) {
        comb <- .me_kron_merge_add(merged[[i]], pair)
        if (!is.null(comb)) { merged[[i]] <- comb; done <- TRUE; break }
      }
      if (!done) merged[[length(merged) + 1L]] <- pair
    }
    krons <- c(krons, lapply(merged, function(p) .ae_mmul(p[[1L]], p[[2L]])))
  }
  # MatAdd(*krons) + nonkrons adds a list to a matrix expression in SymPy
  if (length(nonkrons)) .perr("TypeError", "unsupported operand type(s) for +: 'MatAdd' and 'list'")
  do.call(MatAdd, krons)
}
.me_kron_mat_mul <- function(e) {
  factor <- .mul(Filter(Negate(.is_matexpr), e$args))
  ms <- Filter(.is_matexpr, e$args)
  i <- 1L
  while (i < length(ms)) {
    A <- ms[[i]]; B <- ms[[i + 1L]]
    if (A$type == "KroneckerProduct" && B$type == "KroneckerProduct") {
      ok <- length(A$args) == length(B$args) && identical(.mshape(A)[[2L]]$key, .mshape(B)[[1L]]$key) &&
        all(mapply(function(a, b) identical(.mshape(a)[[2L]]$key, .mshape(b)[[1L]]$key), A$args, B$args))
      ms[[i]] <- if (ok) .me_kronecker_new(Map(.ae_mmul, A$args, B$args)) else .ae_mmul(A, B)
      ms[[i + 1L]] <- NULL
    } else i <- i + 1L
  }
  .ae_mmul(factor, do.call(MatMul, ms))
}
.me_kron_mat_pow <- function(e) {
  b <- e$args[[1L]]
  if (b$type == "KroneckerProduct" && all(vapply(b$args, .is_square_expr, NA)))
    return(.me_kronecker_new(lapply(b$args, function(a) new_sym("MatPow", list(a, e$args[[2L]])))))
  e
}

#' @rdname kronecker_product
#' @export
combine_kronecker <- function(expr) {
  haskron <- function(z) .is_matexpr(z) && .has_type(z, "KroneckerProduct")
  rule1 <- function(z) {
    repeat {
      if (!haskron(z)) return(z)
      y <- switch(z$type, MatAdd = .me_kron_mat_add(z), MatMul = .me_kron_mat_mul(z), MatPow = .me_kron_mat_pow(z), z)
      if (identical(y$key, z$key)) return(z)
      z <- y
    }
  }
  bottom_up <- function(z) {
    if (!is_sym(z) || !length(z$args)) return(rule1(z))
    na <- lapply(z$args, bottom_up)
    if (!identical(vapply(na, .k, ""), vapply(z$args, .k, ""))) z <- new_sym(z$type, na, name = z$name, extra = z$extra)
    rule1(z)
  }
  res <- expr
  repeat {
    y <- bottom_up(res)
    if (identical(y$key, res$key)) break
    res <- y
  }
  doit(res)
}

# ---------------------------------------------------------------------------
# CompanionMatrix, factorizations, DFT and IDFT

#' Structured matrix expressions
#'
#' `CompanionMatrix(poly)` is the companion matrix of a monic univariate
#' polynomial. `DFT(n)` and `IDFT(n)` are the unitary discrete Fourier
#' transform matrix and its inverse. `LofLU()`, `UofLU()`, `LofCholesky()`,
#' `UofCholesky()`, `QofQR()`, `RofQR()`, `EigenVectors()`, `EigenValues()`,
#' `UofSVD()`, `SofSVD()` and `VofSVD()` are the factors of a matrix
#' factorization, with the shape of their argument; `mat_lu()`, `mat_qr()`,
#' `mat_eig()` and `mat_svd()` give the factors as a list.
#'
#' @param poly A monic univariate [Poly()].
#' @param n The size, a nonnegative integer.
#' @param arg,expr A matrix expression.
#' @return A `sym` matrix expression, or a list of them.
#' @examples
#' x <- Symbol("x")
#' as_explicit(DFT(2))
#' QofQR(MatrixSymbol("X", 3, 3))
#' @name structured-matrices
NULL

#' @rdname structured-matrices
#' @export
CompanionMatrix <- function(poly) {
  if (!inherits(poly, "Poly")) .perr("ValueError", paste0(.fmt_any(poly), " must be a Poly instance."))
  cs <- .me_poly_coeffs(poly)
  if (!is_one_num(cs[[1L]])) .perr("ValueError", paste0(.fmt_any(poly), " must be a monic polynomial."))
  if (.poly_ngens(poly) != 1L) .perr("ValueError", paste0(.fmt_any(poly), " must be a univariate polynomial."))
  if (length(cs) < 2L) .perr("ValueError", paste0(.fmt_any(poly), " must have degree not less than 1."))
  new_sym("CompanionMatrix", list(as_expr(poly)), extra = list(coeffs = cs, poly = paste(format(poly), collapse = " "), gen = poly$gens[[1L]]))
}
register_rebuilder("CompanionMatrix", function(e, args) CompanionMatrix(Poly(args[[1L]], e$extra$gen)))
# str() of a Python tuple
.fmt_tuple <- function(x) {
  if (is_sym(x)) return(.str(x))
  s <- vapply(as.list(x), function(z) if (is_sym(z)) .str(z) else format(z), "")
  paste0("(", paste(s, collapse = ", "), if (length(s) == 1L) ",", ")")
}
.fmt_any <- function(x) if (is_sym(x)) .str(x) else if (inherits(x, "Poly")) paste(format(x), collapse = " ") else paste(format(x), collapse = " ")
.me_poly_coeffs <- function(p) all_coeffs(p)
register_str_printer("CompanionMatrix", function(e) paste0("CompanionMatrix(", e$extra$poly, ")"))
.me_shape_hooks$CompanionMatrix <- function(e) { d <- .int(length(e$extra$coeffs) - 1L); list(d, d) }
.me_entry_hooks$CompanionMatrix <- function(e, i, j, ...) {
  cs <- e$extra$coeffs
  n <- length(cs) - 1L
  if (is_exact_number(j) && .int_val(j) == n - 1L) return(.neg(cs[[length(cs) - .int_val(i)]]))
  if (identical(i$key, .add(list(j, S_One))$key)) return(S_One)
  S_Zero
}

.me_factor_types <- c("LofLU", "UofLU", "LofCholesky", "UofCholesky", "QofQR", "RofQR", "EigenVectors",
                      "EigenValues", "UofSVD", "SofSVD", "VofSVD")
for (.ft in .me_factor_types) .me_shape_hooks[[.ft]] <- function(e) .mshape(e$args[[1L]])
rm(.ft)
.me_factor <- function(type, arg) new_sym(type, list(.as_mexpr(as_matexpr(arg))))

#' @rdname structured-matrices
#' @export
LofLU <- function(arg) .me_factor("LofLU", arg)
#' @rdname structured-matrices
#' @export
UofLU <- function(arg) .me_factor("UofLU", arg)
#' @rdname structured-matrices
#' @export
LofCholesky <- function(arg) .me_factor("LofCholesky", arg)
#' @rdname structured-matrices
#' @export
UofCholesky <- function(arg) .me_factor("UofCholesky", arg)
#' @rdname structured-matrices
#' @export
QofQR <- function(arg) .me_factor("QofQR", arg)
#' @rdname structured-matrices
#' @export
RofQR <- function(arg) .me_factor("RofQR", arg)
#' @rdname structured-matrices
#' @export
EigenVectors <- function(arg) .me_factor("EigenVectors", arg)
#' @rdname structured-matrices
#' @export
EigenValues <- function(arg) .me_factor("EigenValues", arg)
#' @rdname structured-matrices
#' @export
UofSVD <- function(arg) .me_factor("UofSVD", arg)
#' @rdname structured-matrices
#' @export
SofSVD <- function(arg) .me_factor("SofSVD", arg)
#' @rdname structured-matrices
#' @export
VofSVD <- function(arg) .me_factor("VofSVD", arg)
#' @rdname structured-matrices
#' @export
mat_lu <- function(expr) list(LofLU(expr), UofLU(expr))
#' @rdname structured-matrices
#' @export
mat_qr <- function(expr) list(QofQR(expr), RofQR(expr))
#' @rdname structured-matrices
#' @export
mat_eig <- function(expr) list(EigenValues(expr), EigenVectors(expr))
#' @rdname structured-matrices
#' @export
mat_svd <- function(expr) list(UofSVD(expr), SofSVD(expr), VofSVD(expr))

#' @rdname structured-matrices
#' @export
DFT <- function(n) new_sym("DFT", list(.me_check_dim(as_sym(n))))
#' @rdname structured-matrices
#' @export
IDFT <- function(n) new_sym("IDFT", list(.me_check_dim(as_sym(n))))
.me_shape_hooks$DFT <- .me_shape_hooks$IDFT <- function(e) list(e$args[[1L]], e$args[[1L]])
.me_dft_entry <- function(sgn) function(e, i, j, ...) {
  n <- e$args[[1L]]
  w <- fcall("exp", list(.div(.mul(list(.int(-2L), S_Pi, S_I)), n)))
  .div(.pow(w, .mul(list(.int(sgn), i, j))), .sqrt(n))
}
.me_entry_hooks$DFT <- .me_dft_entry(1L)
.me_entry_hooks$IDFT <- .me_dft_entry(-1L)

# ---------------------------------------------------------------------------
# MatrixSet

#' The set of matrices of a shape with entries in a set
#'
#' `MatrixSet(n, m, set)` is the set of `n` by `m` matrices whose entries
#' lie in `set`. `matrix_set_contains(M, X)` tests membership, as SymPy's
#' `X in M`.
#'
#' @param n,m The shape.
#' @param set A set.
#' @param M A `MatrixSet`.
#' @param X A matrix or matrix expression.
#' @return A `sym` set, or `TRUE`, `FALSE` or `NA`.
#' @examples
#' M <- MatrixSet(2, 2, Reals)
#' matrix_set_contains(M, Matrix(list(list(1, 2), list(3, 4))))
#' @export
MatrixSet <- function(n, m, set) {
  n <- as_sym(n); m <- as_sym(m)
  .me_check_dim(n); .me_check_dim(m)
  if (!is_sym(set) || !is_set(set)) .perr("TypeError", paste0(.fmt_tuple(set), " should be an instance of Set."))
  new_sym("MatrixSet", list(n, m, set))
}
#' @rdname MatrixSet
#' @export
matrix_set_contains <- function(M, X) {
  if (!.me_is_matrix(X)) .perr("TypeError", paste0(.fmt_any(X), " should be an instance of MatrixExpr."))
  X <- .as_mexpr(X)
  sh <- .mshape(X); ms <- M$args[1:2]
  if (!identical(sh[[1L]]$key, ms[[1L]]$key) || !identical(sh[[2L]]$key, ms[[2L]]$key)) {
    if (any(vapply(c(sh, ms), function(z) z$type == "Symbol", NA))) return(NA)
    return(FALSE)
  }
  E <- as_explicit(X)
  vals <- vapply(E$data, function(z) {
    r <- contains(M$args[[3L]], z)
    if (is.logical(r)) r else if (identical(r$key, S_true$key)) TRUE else if (identical(r$key, S_false$key)) FALSE else NA
  }, NA)
  if (any(vals %in% FALSE)) FALSE else if (all(vals %in% TRUE)) TRUE else NA
}

#' Simplify a sum of products
#'
#' Simplification for the expression blowup that occurs while summing
#' products, as in matrix multiplication. Only algebraic structure is
#' handled; non-algebraic subexpressions are left alone.
#'
#' @param expr An expression.
#' @param withsimp If `TRUE`, return `list(expr, flag)`, where the flag tells
#'   whether further simplification may help.
#' @return The simplified expression, or a list as described for `withsimp`.
#' @examples
#' x <- Symbol("x")
#' dotprodsimp(x/(x - 1) - 1/(x - 1))
#' @export
dotprodsimp <- function(expr, withsimp = FALSE) {
  expr <- as_sym(expr)
  simplified <- FALSE
  if (expr$type %in% c("Add", "Mul", "Pow")) {
    expr2 <- expand(expr, deep = TRUE, power_base = FALSE, power_exp = FALSE, mul = TRUE,
                    log = FALSE, multinomial = TRUE, basic = FALSE)
    if (!identical(expr2$key, expr$key)) { expr <- expr2; simplified <- TRUE }
    co <- .dps_count_ops(expr)
    exprops <- co[[1L]]
    if (exprops >= 6L) {
      if (co[[2L]]) {
        dummies <- new.env()
        dummies$keys <- character(); dummies$exprs <- list(); dummies$syms <- list()
        expr2 <- .dps_subs_dummies(expr, dummies)
        if (identical(expr2$key, expr$key) || .dps_count_ops(expr2)[[1L]] >= 6L) {
          expr3 <- cancel(expr2)
          if (!identical(expr3$key, expr2$key)) {
            expr <- if (length(dummies$syms)) subs(expr3, Map(list, dummies$syms, dummies$exprs)) else expr3
            simplified <- TRUE
          }
        }
      }
    } else if (exprops == 5L && expr$type == "Add" && all(vapply(.dps_pyargs(expr)[1:2], function(a)
                 a$type == "Mul" && (l <- .dps_pyargs(a)[[length(a$args)]])$type == "Pow" && is_negone_num(l$args[[2L]]), NA))) {
      expr2 <- together(expr)
      if (.dps_count_ops(expr2)[[1L]] < exprops) { expr <- expr2; simplified <- TRUE }
    } else simplified <- TRUE
  }
  if (withsimp) list(expr, simplified) else expr
}

# the arguments in SymPy's order
.dps_pyargs <- function(e) .key_sort(lapply(e$args, .sort_key), e$args)

# count_ops_alg of simplify.dotprodsimp: list(ops, ratfunc)
.dps_count_ops <- function(expr) {
  ops <- 0L
  args <- list(expr)
  ratfunc <- FALSE
  while (length(args)) {
    a <- args[[length(args)]]
    args[[length(args)]] <- NULL
    if (!is_sym(a)) next
    t <- a$type
    if (t %in% c("Integer", "Rational")) {
      if (!is_one_num(a)) ops <- ops + could_extract_minus_sign(a) + (t == "Rational")
    } else if (t == "Mul") {
      if (could_extract_minus_sign(a)) {
        ops <- ops + 1L
        a <- if (is_negone_num(a$args[[1L]])) .mul(a$args[-1L]) else .neg(a)
      }
      nd <- fraction(a)
      n <- nd[[1L]]; d <- nd[[2L]]
      if (n$type == "Integer") {
        ops <- ops + 1L + could_extract_minus_sign(n)
        args[[length(args) + 1L]] <- d
      } else if (!is_one_num(d)) {
        if (d$type != "Integer") { args[[length(args) + 1L]] <- d; ratfunc <- TRUE }
        ops <- ops + 1L
        args[[length(args) + 1L]] <- n
      } else {
        ops <- ops + length(a$args) - 1L
        args <- c(args, a$args)
      }
    } else if (t == "Add") {
      negs <- 0L
      for (ai in a$args) {
        if (could_extract_minus_sign(ai)) { negs <- negs + 1L; ai <- .neg(ai) }
        args[[length(args) + 1L]] <- ai
      }
      ops <- ops + length(a$args) - (negs != length(a$args))
    } else if (t == "Pow") {
      ops <- ops + 1L
      args[[length(args) + 1L]] <- a$args[[1L]]
      if (!ratfunc) ratfunc <- !isFALSE(ask_fact(a$args[[2L]], "negative"))
    }
  }
  list(ops, ratfunc)
}

.dps_subs_dummies <- function(expr, dummies) {
  if (!length(expr$args)) return(expr)
  if (expr$type %in% c("Add", "Mul", "Pow")) {
    args <- lapply(expr$args, .dps_subs_dummies, dummies)
    if (all(mapply(function(a, b) identical(a$key, b$key), args, expr$args))) return(expr)
    return(sym_rebuild(expr, args))
  }
  i <- match(expr$key, dummies$keys)
  if (!is.na(i)) return(dummies$syms[[i]])
  d <- Dummy()
  dummies$keys <- c(dummies$keys, expr$key)
  dummies$exprs[[length(dummies$exprs) + 1L]] <- expr
  dummies$syms[[length(dummies$syms) + 1L]] <- d
  d
}

# ---------------------------------------------------------------------------
# Element access: MatrixExpr.__getitem__, valid_index and the _entry methods.

.me_is_false <- function(r) is_sym(r) && identical(r$key, S_false$key)

# MatrixExpr.valid_index
.me_valid_index <- function(A, i, j) {
  s <- .mshape(A)
  ok <- function(idx, n) !.me_is_false(Ge(idx, .neg(n))) && !.me_is_false(Lt(idx, n))
  ok(i, s[[1L]]) && ok(j, s[[2L]])
}

# MatrixExpr.__getitem__ for A[i, j] and A[k]
.me_getitem <- function(A, ...) {
  key <- list(...)
  if (.is_explicit(A) || !.is_matexpr(A)) return(do.call(MatrixElement, c(list(A), key)))
  if (length(key) == 2L) {
    i <- as_sym(key[[1L]]); j <- as_sym(key[[2L]])
    if (.me_valid_index(A, i, j)) return(.me_entry(A, i, j))
    .perr("IndexError", paste0("Invalid indices (", .str(i), ", ", .str(j), ")"))
  }
  k <- key[[1L]]
  if (is.numeric(k) || (is_sym(k) && k$type == "Integer")) {
    cols <- .mshape(A)[[2L]]
    if (cols$type != "Integer")
      .perr("IndexError", "\nSingle indexing is only supported when the number of columns is known.")
    kv <- as.numeric(.int_val(as_sym(k))); cv <- as.numeric(.int_val(cols))
    i <- .int(kv %/% cv); j <- .int(kv %% cv)
    if (.me_valid_index(A, i, j)) return(.me_entry(A, i, j))
    .perr("IndexError", paste0("Invalid index ", .str(as_sym(k))))
  }
  if (is_sym(k)) .perr("IndexError", "\nOnly integers may be used when addressing the matrix with a single\nindex.")
  .perr("IndexError", paste0("Invalid index, wanted ", .str(A), "[i,j]"))
}

# The _entry methods; expand = NULL stands for an argument not given.
.me_entry <- function(A, i, j, expand = NULL, gen = NULL) {
  if (inherits(A, "symmatrix") || .is_explicit(A)) return(MatrixElement(A, i, j))
  switch(A$type,
    MatrixSymbol = new_sym("MatrixElement", list(A, i, j)),
    MatAdd = .add(lapply(A$args, .me_entry, i, j, expand, gen)),
    MatMul = .me_matmul_entry(A, i, j, if (is.null(expand)) TRUE else expand, gen),
    Transpose = .me_entry(A$args[[1L]], j, i, if (is.null(expand)) FALSE else expand, gen),
    Adjoint = fcall("conjugate", list(.me_entry(A$args[[1L]], j, i, expand, gen))),
    MatPow = , Inverse = .me_matpow_entry(A, i, j),
    Identity = {
      eq <- Eq(i, j)
      if (identical(eq$key, S_true$key)) S_One else if (.me_is_false(eq)) S_Zero
      else KroneckerDelta(i, j, Tuple(S_Zero, .sub(.mshape(A)[[2L]], S_One)))
    },
    ZeroMatrix = S_Zero,
    OneMatrix = S_One,
    DiagonalOf = .me_entry(A$args[[1L]], i, i, expand, gen),
    DiagMatrix = {
      v <- A$args[[1L]]
      col <- !is_one_num(.mshape(v)[[1L]])
      r <- if (col) .me_entry(v, i, S_Zero, expand, gen) else .me_entry(v, S_Zero, j, expand, gen)
      if (!identical(i$key, j$key)) r <- .mul(list(r, KroneckerDelta(i, j)))
      r
    },
    {
      h <- .me_entry_hooks[[A$type]]
      if (!is.null(h)) h(A, i, j, expand = expand, gen = gen) else MatrixElement(A, i, j)
    })
}

# MatMul._entry
.me_matmul_entry <- function(A, i, j, expand, gen) {
  cm <- .ae_as_coeff_matrices(A)
  coeff <- cm[[1L]]; mats <- cm[[2L]]
  if (length(mats) == 1L) return(.mul(list(coeff, .me_getitem(mats[[1L]], i, j))))
  if (is.null(gen)) { gen <- new.env(); gen$n <- 0L }
  idx <- vector("list", length(mats) + 1L)
  idx[[1L]] <- i; idx[[length(idx)]] <- j
  for (k in seq_len(length(mats) - 1L)) {
    gen$n <- gen$n + 1L
    idx[[k + 1L]] <- Dummy(paste0("i_", gen$n))
  }
  rng <- lapply(mats[-length(mats)], function(m) .sub(.mshape(m)[[2L]], S_One))
  ents <- lapply(seq_along(mats), function(k) .me_entry(mats[[k]], idx[[k]], idx[[k + 1L]], gen = gen))
  if (any(vapply(ents, function(v) .has_type(v, "ImmutableMatrix"), NA))) expand <- TRUE
  lims <- lapply(seq_along(rng), function(k) list(idx[[k + 1L]], S_Zero, rng[[k]]))
  res <- do.call(Sum, c(list(.mul(ents)), lims))
  if (!any(vapply(rng, function(r) r$type == "Integer", NA))) expand <- FALSE
  if (expand) doit(res) else res
}

# MatPow._entry, which Inverse inherits
.me_matpow_entry <- function(A, i, j) {
  B <- doit(A)
  if (is_sym(B) && B$type %in% c("MatPow", "Inverse")) {
    be <- .mbase_exp(B)
    if (be[[2L]]$type == "Integer" && .int_val(be[[2L]]) > 0) {
      B <- do.call(MatMul, rep(list(be[[1L]]), .int_val(be[[2L]])))
    } else if (all(vapply(.mshape(A), function(d) d$type == "Integer", NA))) {
      return(.me_getitem(as_explicit(B), i, j))
    } else return(new_sym("MatrixElement", list(A, i, j)))
  }
  .me_getitem(B, i, j)
}

register_doit("MatrixElement", function(e) {
  a <- lapply(e$args, function(z) if (is_sym(z)) doit(z) else z)
  .me_getitem(a[[1L]], a[[2L]], a[[3L]])
})

# MatrixExpr.inverse and the _eval_inverse methods
.me_inverse <- function(M) {
  if (.is_explicit(M)) return(inverse(M))
  s <- .mshape(M)
  if (.dim_known_ne(s[[1L]], s[[2L]])) .perr("NonSquareMatrixError", "Inverse of non-square matrix")
  switch(M$type,
    Inverse = M$args[[1L]],
    Identity = M,
    ZeroMatrix = .perr("NonInvertibleMatrixError", "Matrix det == 0; not invertible."),
    MatPow = MatPow(M$args[[1L]], .neg(M$args[[2L]])),
    DFT = IDFT(M$args[[1L]]),
    IDFT = DFT(M$args[[1L]]),
    MatMul = {
      sq <- vapply(Filter(.is_matexpr, M$args), function(a) { d <- .mshape(a); identical(d[[1L]]$key, d[[2L]]$key) }, NA)
      if (all(sq)) doit(do.call(MatMul, lapply(rev(M$args), function(a) if (.is_matexpr(a)) .me_inverse(a) else .pow(a, S_NegativeOne))))
      else Inverse(M)
    },
    KroneckerProduct = tryCatch(do.call(KroneckerProduct, lapply(M$args, .me_inverse)),
                                ShapeError = function(e) Inverse(M)),
    # BlockDiagMatrix._eval_inverse
    BlockDiagMatrix = {
      if (!all(vapply(M$args, .is_square_expr, NA))) .perr("NonInvertibleMatrixError", "Matrix det == 0; not invertible.")
      do.call(BlockDiagMatrix, lapply(M$args, .me_inverse))
    },
    Inverse(M))
}

# a2idx: a nonnegative index after validating it against n
.me_a2idx <- function(j, n) {
  if (j < 0) j <- j + n
  if (!(j >= 0 && j < n)) .perr("IndexError", paste0("Index out of range: a[", format(j, scientific = FALSE), "]"))
  as.integer(j)
}

# A scalar expression built over matrix leaves (a lambda applied to a 1x1
# matrix expression) rebuilt with the matrix operations: Pow of a matrix is
# MatPow, Mul and Add go through their matrix postprocessors.
.me_lift <- function(e) {
  if (!is_sym(e) || .is_matexpr(e) || !length(e$args)) return(e)
  args <- lapply(e$args, .me_lift)
  if (!any(vapply(args, .me_is_matrix, NA))) return(e)
  switch(e$type,
    Pow = if (.me_is_matrix(args[[1L]])) MatPow(args[[1L]], args[[2L]]) else e,
    Mul = do.call(Mul, args),
    Add = do.call(Add, args),
    e)
}

# _eval_is_eq(MatrixExpr, MatrixExpr): FALSE for different shapes, TRUE when
# lhs - rhs is a ZeroMatrix, NULL otherwise
.me_is_eq <- function(a, b) {
  if (!identical(vapply(.mshape(a), .k, ""), vapply(.mshape(b), .k, ""))) return(FALSE)
  d <- tryCatch(doit(MatAdd(a, MatMul(S_NegativeOne, b, evaluate = TRUE))), error = function(err) NULL)
  if (is_sym(d) && d$type == "ZeroMatrix") return(TRUE)
  NULL
}

# a term of an Add with a matrix or array expression factor
.me_has_nc_factor <- function(t) {
  if (.is_matexpr(t) || .is_arrexpr(t)) return(TRUE)
  t$type == "Mul" && any(vapply(t$args, function(f) .is_matexpr(f) || .is_arrexpr(f), NA))
}

# OneMatrix.doit: OneMatrix(1, 1) evaluates to Identity(1)
register_doit("OneMatrix", function(e) {
  if (is_one_num(e$args[[1L]]) && is_one_num(e$args[[2L]])) return(Identity(1L))
  e
})
