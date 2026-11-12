# Handlers of sympy/assumptions/handlers/matrices.py: square, symmetric,
# invertible, orthogonal, unitary, fullrank, positive_definite,
# upper_triangular, lower_triangular, diagonal and the element predicates.
# The handlers for SymPy's Factorization and DFT classes have no
# counterpart here, since symcompR has no such matrix expressions.

# MatrixExpr.is_square: TRUE, FALSE or NA
.mis_square <- function(e) {
  s <- .mshape(e)
  if (s[[1L]]$type == "Integer" && s[[2L]]$type == "Integer") return(s[[1L]]$value == s[[2L]]$value)
  if (.keq(s[[1L]], s[[2L]])) return(TRUE)
  NA
}
.mshape_is <- function(e, r, c) {
  s <- .mshape(e)
  .keq(s[[1L]], .int(r)) && .keq(s[[2L]], .int(c))
}
.is_1x1 <- function(e) .mshape_is(e, 1L, 1L)
.is_empty_or_1x1 <- function(e) .mshape_is(e, 0L, 0L) || .mshape_is(e, 1L, 1L)
.mparent <- function(e) e$args[[1L]]
.on_diag <- function(e) .keq(e$args[[2L]], e$args[[3L]])
# The exponent of a MatPow; Inverse, a MatPow in SymPy, keeps it implicit
.matpow_exp <- function(e) if (e$type == "Inverse") S_NegativeOne else e$args[[2L]]
.is_matrix_arg <- function(x) .is_matexpr(x)

# as_coeff_matrices: the scalar coefficient and the matrix factors
.coeff_matrices <- function(e) {
  isM <- vapply(e$args, .is_matrix_arg, TRUE)
  list(.mul(e$args[!isM]), e$args[isM])
}

.all_ask <- function(name, xs, assumptions) .py_all(length(xs), function(i) .askr(name, xs[[i]], assumptions))
.any_false <- function(name, xs, assumptions) .py_any(length(xs), function(i) identical(.askr(name, xs[[i]], assumptions), FALSE))
.closed_group_list <- function(xs, assumptions, key)
  .fuzzy_group_lazy(length(xs), function(i) .askr(key, xs[[i]], assumptions), quick_exit = TRUE)

# The MatPow handler shared by symmetric, the triangular predicates,
# diagonal and the element predicates.
.matpow_nonneg_or_invertible <- function(key) {
  force(key)
  function(expr, assumptions) {
    base <- expr$args[[1L]]; x <- .matpow_exp(expr)
    if (!isTRUE(.askr("integer", x, assumptions))) return(NA)
    non_negative <- .ask_recursive(.QN("negative", x), assumptions)
    if (isTRUE(non_negative) || (identical(non_negative, FALSE) && isTRUE(.askr("invertible", base, assumptions))))
      return(.askr(key, base, assumptions))
    NA
  }
}

# Handlers that ask the same predicate of the argument or of the parent.
.of_arg <- function(key) { force(key); function(expr, assumptions) .askr(key, expr$args[[1L]], assumptions) }
.of_parent_on_diag <- function(key) {
  force(key)
  function(expr, assumptions) if (!.on_diag(expr)) NA else .askr(key, .mparent(expr), assumptions)
}
.in_assumptions <- function(key) {
  force(key)
  function(expr, assumptions) if (.in_conjuncts(.QA(key, expr), assumptions)) TRUE else NA
}
.is_1x1_handler <- function(expr, assumptions) .is_1x1(expr)
.true_handler <- function(expr, assumptions) TRUE
.false_handler <- function(expr, assumptions) FALSE

# ---------------------------------------------------------------------------
# square

.ask_register("square", "MatrixExpr", function(expr, assumptions) {
  s <- .mshape(expr)
  .keq(s[[1L]], s[[2L]])
})

# ---------------------------------------------------------------------------
# symmetric

.ask_register("symmetric", "MatMul", function(expr, assumptions) {
  cm <- .coeff_matrices(expr)
  mm <- cm[[2L]]
  if (.all_ask("symmetric", mm, assumptions)) return(TRUE)
  if (isTRUE(.askr("diagonal", expr, assumptions))) return(TRUE)
  n <- length(mm)
  if (n >= 2L && .keq(mm[[1L]], .mtranspose(mm[[n]]))) {
    if (n == 2L) return(TRUE)
    return(.askr("symmetric", do.call(MatMul, mm[2:(n - 1L)]), assumptions))
  }
  NA
})

.ask_register("symmetric", "MatPow", .matpow_nonneg_or_invertible("symmetric"))

.ask_register("symmetric", "MatAdd", function(expr, assumptions) .all_ask("symmetric", expr$args, assumptions))

.ask_register("symmetric", "MatrixSymbol", function(expr, assumptions) {
  if (!isTRUE(.mis_square(expr))) return(FALSE)
  if (isTRUE(.askr("diagonal", expr, assumptions))) return(TRUE)
  if (.in_conjuncts(.QA("symmetric", expr), assumptions)) return(TRUE)
  NA
})

.ask_register("symmetric", c("OneMatrix", "ZeroMatrix"), function(expr, assumptions) .askr("square", expr, assumptions))
.ask_register("symmetric", c("Inverse", "Transpose"), .of_arg("symmetric"))

.ask_register("symmetric", "MatrixSlice", function(expr, assumptions) {
  if (isTRUE(.askr("diagonal", expr, assumptions))) return(TRUE)
  if (!.on_diag(expr)) return(NA)
  .askr("symmetric", .mparent(expr), assumptions)
})

.ask_register("symmetric", "Identity", .true_handler)

# ---------------------------------------------------------------------------
# invertible

.ask_register("invertible", "MatMul", function(expr, assumptions) {
  mm <- .coeff_matrices(expr)[[2L]]
  if (.all_ask("invertible", mm, assumptions)) return(TRUE)
  if (.any_false("invertible", mm, assumptions)) return(FALSE)
  NA
})

.ask_register("invertible", "MatPow", function(expr, assumptions) {
  base <- expr$args[[1L]]; x <- .matpow_exp(expr)
  if (!isTRUE(.askr("integer", x, assumptions))) return(NA)
  if (identical(.fact(x, "negative"), FALSE)) return(.askr("invertible", base, assumptions))
  NA
})

.ask_register("invertible", "MatAdd", function(expr, assumptions) NA)

.ask_register("invertible", "MatrixSymbol", function(expr, assumptions) {
  if (!isTRUE(.mis_square(expr))) return(FALSE)
  if (.in_conjuncts(.QA("invertible", expr), assumptions)) return(TRUE)
  NA
})

.ask_register("invertible", c("Identity", "Inverse"), .true_handler)
.ask_register("invertible", "ZeroMatrix", .false_handler)
.ask_register("invertible", "OneMatrix", .is_1x1_handler)
.ask_register("invertible", "Transpose", .of_arg("invertible"))
.ask_register("invertible", "MatrixSlice", .of_parent_on_diag("invertible"))

.ask_register("invertible", "MatrixBase", function(expr, assumptions) {
  if (!isTRUE(.mis_square(expr))) return(FALSE)
  rank_matrix(as_explicit(expr)) == .int_val(expr$args[[1L]])
})

.ask_register("invertible", "MatrixExpr", function(expr, assumptions) {
  if (!isTRUE(.mis_square(expr))) return(FALSE)
  NA
})

# X - Y * Z.I * W for blocks: explicit arithmetic when all four are
# explicit, a matrix expression otherwise.
.block_schur <- function(X, Y, Z, W) {
  if (all(vapply(list(X, Y, Z, W), .is_matrixbase, TRUE))) {
    ex <- lapply(list(X, Y, Z, W), as_explicit)
    r <- ex[[1L]] - .matmul(.matmul(ex[[2L]], inverse(ex[[3L]])), ex[[4L]])
    return(.as_block(r))
  }
  zi <- if (.is_matrixbase(Z)) .as_block(inverse(as_explicit(Z))) else Inverse(Z)
  X - Y * zi * W
}

.ask_register("invertible", "BlockMatrix", function(expr, assumptions) {
  if (!isTRUE(.mis_square(expr))) return(FALSE)
  R <- .block_rows(expr)
  if (length(R) == 1L && length(R[[1L]]) == 1L) return(.askr("invertible", R[[1L]][[1L]], assumptions))
  # reblock_2x2 changes only block matrices with more than two blocks
  # along both axes, which are left unanswered here
  if (length(R) == 2L && length(R[[1L]]) == 2L) {
    A <- R[[1L]][[1L]]; B <- R[[1L]][[2L]]; C <- R[[2L]][[1L]]; D <- R[[2L]][[2L]]
    if (isTRUE(.askr("invertible", A, assumptions))) {
      inv <- .askr("invertible", .block_schur(D, C, A, B), assumptions)
      if (!is.na(inv)) return(inv)
    }
    if (isTRUE(.askr("invertible", B, assumptions))) {
      inv <- .askr("invertible", .block_schur(C, D, B, A), assumptions)
      if (!is.na(inv)) return(inv)
    }
    if (isTRUE(.askr("invertible", C, assumptions))) {
      inv <- .askr("invertible", .block_schur(B, A, C, D), assumptions)
      if (!is.na(inv)) return(inv)
    }
    if (isTRUE(.askr("invertible", D, assumptions))) {
      inv <- .askr("invertible", .block_schur(A, B, D, C), assumptions)
      if (!is.na(inv)) return(inv)
    }
  }
  NA
})

.ask_register("invertible", "BlockDiagMatrix", function(expr, assumptions) {
  rs <- lapply(expr$args, function(b) .mshape(b)[[1L]]$key)
  cs <- lapply(expr$args, function(b) .mshape(b)[[2L]]$key)
  if (!identical(rs, cs)) return(NA)
  .fuzzy_and_vals(lapply(expr$args, function(a) .askr("invertible", a, assumptions)))
})

# ---------------------------------------------------------------------------
# orthogonal and unitary

.ask_register("orthogonal", "MatMul", function(expr, assumptions) {
  cm <- .coeff_matrices(expr)
  if (.all_ask("orthogonal", cm[[2L]], assumptions) && .keq(cm[[1L]], S_One)) return(TRUE)
  if (.any_false("invertible", cm[[2L]], assumptions)) return(FALSE)
  NA
})

.ask_register("orthogonal", "MatPow", function(expr, assumptions) {
  if (isTRUE(.askr("integer", .matpow_exp(expr), assumptions))) return(.askr("orthogonal", expr$args[[1L]], assumptions))
  NA
})

.ask_register("orthogonal", "MatAdd", function(expr, assumptions) {
  if (length(expr$args) == 1L && isTRUE(.askr("orthogonal", expr$args[[1L]], assumptions))) return(TRUE)
  NA
})

.unitary_symbol <- function(key) {
  force(key)
  function(expr, assumptions) {
    if (!isTRUE(.mis_square(expr)) || identical(.askr("invertible", expr, assumptions), FALSE)) return(FALSE)
    if (.in_conjuncts(.QA(key, expr), assumptions)) return(TRUE)
    NA
  }
}

.ask_register("orthogonal", "MatrixSymbol", .unitary_symbol("orthogonal"))
.ask_register("orthogonal", "Identity", .true_handler)
.ask_register("orthogonal", "ZeroMatrix", .false_handler)
.ask_register("orthogonal", c("Inverse", "Transpose"), .of_arg("orthogonal"))
.ask_register("orthogonal", "MatrixSlice", .of_parent_on_diag("orthogonal"))

.ask_register("unitary", "MatMul", function(expr, assumptions) {
  cm <- .coeff_matrices(expr)
  if (.all_ask("unitary", cm[[2L]], assumptions) && .keq(Abs(cm[[1L]]), S_One)) return(TRUE)
  if (.any_false("invertible", cm[[2L]], assumptions)) return(FALSE)
  NA
})

.ask_register("unitary", "MatPow", function(expr, assumptions) {
  if (isTRUE(.askr("integer", .matpow_exp(expr), assumptions))) return(.askr("unitary", expr$args[[1L]], assumptions))
  NA
})

.ask_register("unitary", "MatrixSymbol", .unitary_symbol("unitary"))
.ask_register("unitary", c("Inverse", "Transpose"), .of_arg("unitary"))
.ask_register("unitary", "MatrixSlice", .of_parent_on_diag("unitary"))
.ask_register("unitary", "Identity", .true_handler)
.ask_register("unitary", "ZeroMatrix", .false_handler)

# ---------------------------------------------------------------------------
# fullrank

.ask_register("fullrank", "MatMul", function(expr, assumptions) {
  if (.all_ask("fullrank", expr$args, assumptions)) return(TRUE)
  NA
})

.ask_register("fullrank", "MatPow", function(expr, assumptions) {
  base <- expr$args[[1L]]; x <- .matpow_exp(expr)
  if (isTRUE(.askr("integer", x, assumptions)) && isTRUE(.ask_recursive(.QN("negative", x), assumptions)))
    return(.askr("fullrank", base, assumptions))
  NA
})

.ask_register("fullrank", "Identity", .true_handler)
.ask_register("fullrank", "ZeroMatrix", .false_handler)
.ask_register("fullrank", "OneMatrix", .is_1x1_handler)
.ask_register("fullrank", c("Inverse", "Transpose"), .of_arg("fullrank"))
.ask_register("fullrank", "MatrixSlice", function(expr, assumptions) {
  if (isTRUE(.askr("orthogonal", .mparent(expr), assumptions))) return(TRUE)
  NA
})

# ---------------------------------------------------------------------------
# positive_definite

.ask_register("positive_definite", "MatMul", function(expr, assumptions) {
  cm <- .coeff_matrices(expr)
  mm <- cm[[2L]]
  if (.all_ask("positive_definite", mm, assumptions) && .rel_truth(Gt(cm[[1L]], S_Zero))) return(TRUE)
  n <- length(mm)
  if (n >= 2L && .keq(mm[[1L]], .mtranspose(mm[[n]])) && isTRUE(.askr("fullrank", mm[[1L]], assumptions))) {
    # MatMul() of nothing is the GenericIdentity, which is positive definite
    if (n == 2L) return(TRUE)
    return(.askr("positive_definite", do.call(MatMul, mm[seq_len(n - 2L) + 1L]), assumptions))
  }
  NA
})

.ask_register("positive_definite", "MatPow", function(expr, assumptions) {
  if (isTRUE(.askr("positive_definite", expr$args[[1L]], assumptions))) return(TRUE)
  NA
})

.ask_register("positive_definite", "MatAdd", function(expr, assumptions) {
  if (.all_ask("positive_definite", expr$args, assumptions)) return(TRUE)
  NA
})

.ask_register("positive_definite", "MatrixSymbol", function(expr, assumptions) {
  if (!isTRUE(.mis_square(expr))) return(FALSE)
  if (.in_conjuncts(.QA("positive_definite", expr), assumptions)) return(TRUE)
  NA
})

.ask_register("positive_definite", "Identity", .true_handler)
.ask_register("positive_definite", "ZeroMatrix", .false_handler)
.ask_register("positive_definite", "OneMatrix", .is_1x1_handler)
.ask_register("positive_definite", c("Inverse", "Transpose"), .of_arg("positive_definite"))
.ask_register("positive_definite", "MatrixSlice", .of_parent_on_diag("positive_definite"))

# ---------------------------------------------------------------------------
# upper_triangular and lower_triangular

for (.tri in c("upper_triangular", "lower_triangular")) local({
  key <- .tri
  other <- if (key == "upper_triangular") "lower_triangular" else "upper_triangular"
  .ask_register(key, "MatMul", function(expr, assumptions) {
    if (.all_ask(key, .coeff_matrices(expr)[[2L]], assumptions)) return(TRUE)
    NA
  })
  .ask_register(key, "MatAdd", function(expr, assumptions) {
    if (.all_ask(key, expr$args, assumptions)) return(TRUE)
    NA
  })
  .ask_register(key, "MatPow", .matpow_nonneg_or_invertible(key))
  .ask_register(key, "MatrixSymbol", .in_assumptions(key))
  .ask_register(key, c("Identity", "ZeroMatrix"), .true_handler)
  .ask_register(key, "OneMatrix", .is_1x1_handler)
  .ask_register(key, "Transpose", .of_arg(other))
  .ask_register(key, "Inverse", .of_arg(key))
  .ask_register(key, "MatrixSlice", .of_parent_on_diag(key))
})
rm(.tri)

# ---------------------------------------------------------------------------
# diagonal

.ask_register("diagonal", "MatMul", function(expr, assumptions) {
  if (.is_empty_or_1x1(expr)) return(TRUE)
  if (.all_ask("diagonal", .coeff_matrices(expr)[[2L]], assumptions)) return(TRUE)
  NA
})

.ask_register("diagonal", "MatPow", .matpow_nonneg_or_invertible("diagonal"))

.ask_register("diagonal", "MatAdd", function(expr, assumptions) {
  if (.all_ask("diagonal", expr$args, assumptions)) return(TRUE)
  NA
})

.ask_register("diagonal", "MatrixSymbol", function(expr, assumptions) {
  if (.is_empty_or_1x1(expr)) return(TRUE)
  if (.in_conjuncts(.QA("diagonal", expr), assumptions)) return(TRUE)
  NA
})

.ask_register("diagonal", "OneMatrix", .is_1x1_handler)
.ask_register("diagonal", c("Inverse", "Transpose"), .of_arg("diagonal"))

.ask_register("diagonal", "MatrixSlice", function(expr, assumptions) {
  if (.is_empty_or_1x1(expr)) return(TRUE)
  if (!.on_diag(expr)) return(NA)
  .askr("diagonal", .mparent(expr), assumptions)
})

.ask_register("diagonal", c("DiagonalMatrix", "DiagMatrix", "Identity", "ZeroMatrix"), .true_handler)

# ---------------------------------------------------------------------------
# integer_elements, real_elements and complex_elements

.bm_elements <- function(key) {
  force(key)
  function(expr, assumptions) {
    blocks <- unlist(.block_rows(expr), recursive = FALSE)
    .py_all(length(blocks), function(i) .askr(key, blocks[[i]], assumptions))
  }
}

.matmul_elements <- function(mkey, skey) {
  force(mkey); force(skey)
  function(expr, assumptions) {
    isM <- vapply(expr$args, .is_matrix_arg, TRUE)
    .fuzzy_and_vals(list(.closed_group_list(expr$args[!isM], assumptions, skey),
                         .closed_group_list(expr$args[isM], assumptions, mkey)))
  }
}

.closed_elements <- function(key) { force(key); function(expr, assumptions) .test_closed_group(expr, assumptions, key) }

.ask_register("integer_elements", c("Determinant", "HadamardProduct", "MatAdd", "Trace", "Transpose"),
              .closed_elements("integer_elements"))
.ask_register("integer_elements", "MatPow", function(expr, assumptions) {
  base <- expr$args[[1L]]; x <- .matpow_exp(expr)
  if (!isTRUE(.askr("integer", x, assumptions))) return(NA)
  if (identical(.fact(x, "negative"), FALSE)) return(.askr("integer_elements", base, assumptions))
  NA
})
.ask_register("integer_elements", c("Identity", "OneMatrix", "ZeroMatrix"), .true_handler)
.ask_register("integer_elements", "MatMul", .matmul_elements("integer_elements", "integer"))
.ask_register("integer_elements", "MatrixSlice", .of_arg("integer_elements"))
.ask_register("integer_elements", "BlockMatrix", .bm_elements("integer_elements"))

.ask_register("real_elements", c("Determinant", "HadamardProduct", "MatAdd", "Trace", "Transpose"),
              .closed_elements("real_elements"))
.ask_register("real_elements", "MatPow", .matpow_nonneg_or_invertible("real_elements"))
.ask_register("real_elements", "MatMul", .matmul_elements("real_elements", "real"))
.ask_register("real_elements", "MatrixSlice", .of_arg("real_elements"))
.ask_register("real_elements", "BlockMatrix", .bm_elements("real_elements"))

.ask_register("complex_elements", c("Determinant", "HadamardProduct", "Inverse", "MatAdd", "Trace", "Transpose"),
              .closed_elements("complex_elements"))
.ask_register("complex_elements", "MatPow", .matpow_nonneg_or_invertible("complex_elements"))
.ask_register("complex_elements", "MatMul", .matmul_elements("complex_elements", "complex"))
.ask_register("complex_elements", "MatrixSlice", .of_arg("complex_elements"))
.ask_register("complex_elements", "BlockMatrix", .bm_elements("complex_elements"))
