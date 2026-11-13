skip_unless_full_tests()

# Expected values are the str() outputs of SymPy 1.15 for the same calls, taken
# from the SymPy test files named in each test title. Python ** is printed as ^.

.bat_fmt0 <- function(v) {
  if (is.null(v)) return("None")
  if (inherits(v, "pytuple")) {
    s <- vapply(v, .bat_fmt0, "")
    return(paste0("(", paste(s, collapse = ", "), if (length(s) == 1L) ",", ")"))
  }
  if (inherits(v, "pydict")) {
    if (!length(v$keys)) return("{}")
    return(paste0("{", paste(vapply(seq_along(v$keys), function(i)
      paste0(.bat_fmt0(v$keys[[i]]), ": ", .bat_fmt0(v$values[[i]])), ""), collapse = ", "), "}"))
  }
  if (isTRUE(attr(v, "pyset"))) return(paste0("{", paste(sort(vapply(v, .bat_fmt0, ""), method = "radix"), collapse = ", "), "}"))
  if (inherits(v, "sym") || inherits(v, "symmatrix") || inherits(v, "SymArray")) return(sstr(v))
  if (inherits(v, "Poly")) return(paste(format(v), collapse = " "))
  # a Python str prints as its text
  if (is.character(v) && length(v) == 1L) return(v)
  if (is.numeric(v) && !length(v)) return("[]")
  if (is.character(v)) {
    s <- paste0("'", v, "'")
    return(paste0("(", paste(s, collapse = ", "), if (length(s) == 1L) ",", ")"))
  }
  if (is.logical(v) && length(v) == 1L) return(if (is.na(v)) "None" else if (v) "True" else "False")
  if (is.numeric(v) && length(v) == 1L) return(format(v))
  if (is.list(v) && !is.object(v)) {
    nm <- names(v)
    if (length(v) && !is.null(nm) && all(nzchar(nm)))
      return(paste0("{", paste(vapply(seq_along(v), function(i) paste0(nm[i], ": ", .bat_fmt0(v[[i]])), ""),
                               collapse = ", "), "}"))
    return(paste0("[", paste(vapply(v, .bat_fmt0, ""), collapse = ", "), "]"))
  }
  paste(format(v), collapse = " ")
}

.err_class <- function(e) {
  cl <- setdiff(class(e), c("simpleError", "error", "condition", "sympy_error"))
  if (length(cl)) cl[1L] else "Error"
}

# Python's obj[...] and set results, for the transpiled cases
pyitem <- function(x, ...) {
  idx <- list(...)
  iv <- function(k) { k <- if (is.numeric(k)) k else as.numeric(.int_val(k)); k }
  if (inherits(x, "symmatrix")) {
    isint <- function(k) is.numeric(k) || (inherits(k, "sym") && k$type == "Integer")
    if (length(idx) == 2L && !all(vapply(idx, isint, NA))) {
      # Matrix.__getitem__ with a symbolic index
      tr <- function(r) inherits(r, "sym") && identical(r$key, S_true$key)
      i <- as_sym(idx[[1L]]); j <- as_sym(idx[[2L]])
      if (tr(Lt(j, 0)) || tr(Ge(j, x$ncol)) || tr(Lt(i, 0)) || tr(Ge(i, x$nrow)))
        stop(structure(class = c("ValueError", "error", "condition"), list(message = "index out of boundary", call = NULL)))
      return(MatrixElement(x, i, j))
    }
    if (length(idx) == 1L) { k <- iv(idx[[1L]]); if (k < 0) k <- length(x$data) + k; return(x$data[[k + 1L]]) }
    i <- iv(idx[[1L]]); j <- iv(idx[[2L]])
    if (i < 0) i <- x$nrow + i
    if (j < 0) j <- x$ncol + j
    return(x$data[[i * x$ncol + j + 1L]])
  }
  if (inherits(x, "SymArray")) {
    sh <- x$shape
    if (length(idx) == 1L && length(sh) > 1L) return(x$data[[iv(idx[[1L]]) + 1L]])
    k <- 0
    for (d in seq_along(sh)) k <- k * sh[d] + iv(idx[[d]])
    return(x$data[[k + 1L]])
  }
  if (inherits(x, "sym")) {
    if (x$type == "Tuple") { k <- iv(idx[[1L]]); if (k < 0) k <- length(x$args) + k; return(x$args[[k + 1L]]) }
    if (.is_matexpr(x)) {
      if (exists(".me_getitem", envir = asNamespace("symcompR"))) return(do.call(.me_getitem, c(list(x), idx)))
      if (length(idx) == 1L) {
        sh <- get_shape(x); k <- iv(idx[[1L]]); nc <- iv(sh[[2L]])
        return(MatrixElement(x, k %/% nc, k %% nc))
      }
      return(MatrixElement(x, idx[[1L]], idx[[2L]]))
    }
    if (exists(".ae_getitem", envir = asNamespace("symcompR"))) return(.ae_getitem(x, idx))
    return(ArrayElement(x, idx))
  }
  if (is.list(x) || is.atomic(x)) {
    k <- iv(idx[[1L]]); if (k < 0) k <- length(x) + k
    return(x[[k + 1L]])
  }
  stop("cannot index ", class(x)[1L])
}
pyset <- function(x) structure(x, pyset = TRUE)
.bat_fmt <- function(v) gsub("\n", "\\n", gsub("**", "^", .bat_fmt0(v), fixed = TRUE), fixed = TRUE)

test_that("test_derivatives: test_matrix_derivative_by_scalar matches SymPy", {
  expect_equal(.bat_fmt(local({v_i <- Symbol("i"); v_A <- MatrixSymbol("A", Symbol("k"), Symbol("k")); diff(v_A, v_i)})), "0")
  expect_equal(.bat_fmt(local({v_k <- Symbol("k"); ZeroMatrix(v_k, v_k)})), "0")
  expect_equal(.bat_fmt(local({v_i <- Symbol("i"); v_A <- MatrixSymbol("A", Symbol("k"), Symbol("k")); v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_B <- MatrixSymbol("B", Symbol("k"), Symbol("k")); v_c <- MatrixSymbol("c", Symbol("k"), Integer(1)); diff(((v_A * (v_X + v_B)) * v_c), v_i)})), "0")
  expect_equal(.bat_fmt(local({v_k <- Symbol("k"); ZeroMatrix(v_k, 1)})), "0")
  expect_equal(.bat_fmt(local({v_i <- Symbol("i"); v_x <- MatrixSymbol("x", Symbol("k"), Integer(1)); diff(v_x, v_i)})), "0")
  expect_equal(.bat_fmt(local({v_i <- Symbol("i"); v_x <- MatrixSymbol("x", Symbol("k"), Integer(1)); v_y <- MatrixSymbol("y", Symbol("k"), Integer(1)); diff((t(v_x) * v_y), v_i)})), "0")
  expect_equal(.bat_fmt(local({ZeroMatrix(1, 1)})), "0")
  expect_equal(.bat_fmt(local({v_i <- Symbol("i"); v_x <- MatrixSymbol("x", Symbol("k"), Integer(1)); diff((v_x * t(v_x)), v_i)})), "0")
  expect_equal(.bat_fmt(local({v_i <- Symbol("i"); v_x <- MatrixSymbol("x", Symbol("k"), Integer(1)); v_y <- MatrixSymbol("y", Symbol("k"), Integer(1)); diff((v_x + v_y), v_i)})), "0")
  expect_equal(.bat_fmt(local({v_i <- Symbol("i"); v_x <- MatrixSymbol("x", Symbol("k"), Integer(1)); diff(hadamard_power(v_x, 2), v_i)})), "0")
  expect_equal(.bat_fmt(local({v_i <- Symbol("i"); v_x <- MatrixSymbol("x", Symbol("k"), Integer(1)); diff(hadamard_power(v_x, v_i), v_i)})), "HadamardProduct(Lambda(_d, log(_d)).(x), HadamardPower(x, i))")
  expect_equal(.bat_fmt(local({v_x <- MatrixSymbol("x", Symbol("k"), Integer(1)); v_i <- Symbol("i"); HadamardProduct(applyfunc(v_x, log), HadamardPower(v_x, v_i))})), "HadamardProduct(Lambda(_d, log(_d)).(x), HadamardPower(x, i))")
  expect_equal(.bat_fmt(local({v_i <- Symbol("i"); v_x <- MatrixSymbol("x", Symbol("k"), Integer(1)); v_y <- MatrixSymbol("y", Symbol("k"), Integer(1)); diff(hadamard_product(v_x, v_y), v_i)})), "0")
  expect_equal(.bat_fmt(local({v_i <- Symbol("i"); v_k <- Symbol("k"); v_x <- MatrixSymbol("x", Symbol("k"), Integer(1)); v_y <- MatrixSymbol("y", Symbol("k"), Integer(1)); diff(hadamard_product((v_i * OneMatrix(v_k, 1)), v_x, v_y), v_i)})), "HadamardProduct(x, y)")
  expect_equal(.bat_fmt(local({v_x <- MatrixSymbol("x", Symbol("k"), Integer(1)); v_y <- MatrixSymbol("y", Symbol("k"), Integer(1)); hadamard_product(v_x, v_y)})), "HadamardProduct(x, y)")
  expect_equal(.bat_fmt(local({v_i <- Symbol("i"); v_x <- MatrixSymbol("x", Symbol("k"), Integer(1)); diff((v_i * v_x), v_i)})), "x")
  expect_equal(.bat_fmt(local({v_x <- MatrixSymbol("x", Symbol("k"), Integer(1)); v_x})), "x")
  expect_equal(.bat_fmt(local({v_i <- Symbol("i"); v_A <- MatrixSymbol("A", Symbol("k"), Symbol("k")); v_B <- MatrixSymbol("B", Symbol("k"), Symbol("k")); v_x <- MatrixSymbol("x", Symbol("k"), Integer(1)); diff((((sin(v_i) * v_A) * v_B) * v_x), v_i)})), "cos(i)*A*B*x")
  expect_equal(.bat_fmt(local({v_i <- Symbol("i"); v_A <- MatrixSymbol("A", Symbol("k"), Symbol("k")); v_B <- MatrixSymbol("B", Symbol("k"), Symbol("k")); v_x <- MatrixSymbol("x", Symbol("k"), Integer(1)); (((cos(v_i) * v_A) * v_B) * v_x)})), "cos(i)*A*B*x")
  expect_equal(.bat_fmt(local({v_i <- Symbol("i"); v_x <- MatrixSymbol("x", Symbol("k"), Integer(1)); diff(applyfunc(v_x, sin), v_i)})), "0")
  expect_equal(.bat_fmt(local({v_i <- Symbol("i"); v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); diff(Trace(((v_i ^ 2) * v_X)), v_i)})), "2*i*Trace(X)")
  expect_equal(.bat_fmt(local({v_i <- Symbol("i"); v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); ((2 * v_i) * Trace(v_X))})), "2*i*Trace(X)")
  expect_equal(.bat_fmt(local({v_x <- MatrixSymbol("x", Symbol("k"), Integer(1)); v_expr <- MatMul((Integer(2)*Symbol("mu")), MatrixSymbol("x", Symbol("k"), Integer(1))); diff(v_expr, v_x)})), "(2*mu)*I")
  expect_equal(.bat_fmt(local({v_mu <- Symbol("mu"); v_k <- Symbol("k"); ((2 * v_mu) * Identity(v_k))})), "(2*mu)*I")
})

test_that("test_derivatives: test_one_matrix matches SymPy", {
  expect_equal(.bat_fmt(local({v_x <- MatrixSymbol("x", Symbol("k"), Integer(1)); v_k <- Symbol("k"); diff(MatMul(t(v_x), OneMatrix(v_k, 1)), v_x)})), "1")
  expect_equal(.bat_fmt(local({v_k <- Symbol("k"); OneMatrix(v_k, 1)})), "1")
})

test_that("test_derivatives: test_matrix_derivative_non_matrix_result matches SymPy", {
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("k"), Symbol("k")); diff(v_A, v_A)})), "PermuteDims(ArrayTensorProduct(I, I), (3)(1 2))")
  expect_equal(.bat_fmt(local({v_AdA <- PermuteDims(ArrayTensorProduct(Identity(Symbol("k")), Identity(Symbol("k"))), Permutation(c(0, 2, 1, 3))); v_AdA})), "PermuteDims(ArrayTensorProduct(I, I), (3)(1 2))")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("k"), Symbol("k")); diff(t(v_A), v_A)})), "PermuteDims(ArrayTensorProduct(I, I), (1 2 3))")
  expect_equal(.bat_fmt(local({v_I <- Identity(Symbol("k")); v__perm1 <- Permutation(c(0, 2, 3, 1)); PermuteDims(ArrayTensorProduct(v_I, v_I), v__perm1)})), "PermuteDims(ArrayTensorProduct(I, I), (1 2 3))")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("k"), Symbol("k")); diff((2 * v_A), v_A)})), "PermuteDims(ArrayTensorProduct(2, I, I), (3)(1 2))")
  expect_equal(.bat_fmt(local({v_I <- Identity(Symbol("k")); v__perm1 <- Permutation(c(0, 2, 1, 3)); PermuteDims(ArrayTensorProduct(2, v_I, v_I), v__perm1)})), "PermuteDims(ArrayTensorProduct(2, I, I), (3)(1 2))")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("k"), Symbol("k")); diff(MatAdd(v_A, v_A), v_A)})), "PermuteDims(ArrayTensorProduct(2, I, I), (3)(1 2))")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("k"), Symbol("k")); v_B <- MatrixSymbol("B", Symbol("k"), Symbol("k")); diff((v_A + v_B), v_A)})), "PermuteDims(ArrayTensorProduct(I, I), (3)(1 2))")
})

test_that("test_derivatives: test_matrix_derivative_trivial_cases matches SymPy", {
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("k"), Symbol("k")); v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); diff(v_X, v_A)})), "Derivative(X, A)")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_A <- MatrixSymbol("A", Symbol("k"), Symbol("k")); ArrayDerivative(v_X, v_A)})), "Derivative(X, A)")
})

test_that("test_derivatives: test_matrix_derivative_of_determinant matches SymPy", {
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_expr <- Determinant(MatrixSymbol("X", Symbol("k"), Symbol("k"))); diff(v_expr, v_X)})), "Determinant(X)*X.T^(-1)")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); (Determinant(v_X) * inverse(t(v_X)))})), "Determinant(X)*X.T^(-1)")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_expr <- Determinant(MatMul(MatrixSymbol("A", Symbol("k"), Symbol("k")), MatrixSymbol("X", Symbol("k"), Symbol("k")), MatrixSymbol("B", Symbol("k"), Symbol("k")))); diff(v_expr, v_X)})), "(Determinant(A)*Determinant(B)*Determinant(X))*X.T^(-1)")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("k"), Symbol("k")); v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_B <- MatrixSymbol("B", Symbol("k"), Symbol("k")); (Determinant(((v_A * v_X) * v_B)) * inverse(t(v_X)))})), "(Determinant(A)*Determinant(B)*Determinant(X))*X.T^(-1)")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_expr <- Determinant(MatMul(Transpose(MatrixSymbol("X", Symbol("k"), Symbol("k"))), MatrixSymbol("A", Symbol("k"), Symbol("k")), MatrixSymbol("X", Symbol("k"), Symbol("k")))); diff(v_expr, v_X)})), "(2*Determinant(A)*Determinant(X)^2)*X.T^(-1)")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_A <- MatrixSymbol("A", Symbol("k"), Symbol("k")); (((2 * (Determinant(v_X) ^ 2)) * Determinant(v_A)) * inverse(t(v_X)))})), "(2*Determinant(A)*Determinant(X)^2)*X.T^(-1)")
  expect_equal(.bat_fmt(local({v_Z <- MatrixSymbol("Z", Symbol("k"), Symbol("j")); v_expr <- Determinant(MatMul(Transpose(MatrixSymbol("Z", Symbol("k"), Symbol("j"))), MatrixSymbol("A", Symbol("k"), Symbol("k")), MatrixSymbol("Z", Symbol("k"), Symbol("j")))); diff(v_expr, v_Z)})), "Determinant(Z.T*A*Z)*A*Z*(Z.T*A*Z)^(-1) + Determinant(Z.T*A*Z)*A.T*Z*(Z.T*A.T*Z)^(-1)")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_expr <- log(Determinant(MatMul(Transpose(MatrixSymbol("X", Symbol("k"), Symbol("k"))), MatrixSymbol("X", Symbol("k"), Symbol("k"))))); diff(v_expr, v_X)})), "2*X.T^(-1)")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); (2 * inverse(t(v_X)))})), "2*X.T^(-1)")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_expr <- log(Determinant(MatrixSymbol("X", Symbol("k"), Symbol("k")))); diff(v_expr, v_X)})), "X.T^(-1)")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); inverse(t(v_X))})), "X.T^(-1)")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_expr <- Determinant(MatPow(MatrixSymbol("X", Symbol("k"), Symbol("k")), Symbol("j"))); doit(diff(v_expr, v_X))})), "(j*Determinant(X^j))*X.T^(-1)")
  expect_equal(.bat_fmt(local({v_j <- Symbol("j"); v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); ((v_j * Determinant((v_X ^ v_j))) * t(inverse(v_X)))})), "(j*Determinant(X^j))*X.T^(-1)")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_expr <- Determinant(MatPow(MatrixSymbol("X", Symbol("k"), Symbol("k")), Integer(5))); diff(v_expr, v_X)})), "(5*Determinant(X^5))*X.T^(-1)")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); ((5 * Determinant((v_X ^ 5))) * t(inverse(v_X)))})), "(5*Determinant(X^5))*X.T^(-1)")
})

test_that("test_derivatives: test_issue_28708 matches SymPy", {
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_expr <- Determinant(MatPow(MatrixSymbol("X", Symbol("k"), Symbol("k")), Integer(0))); diff(v_expr, v_X)})), "0")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_expr <- Determinant(MatPow(MatrixSymbol("X", Symbol("k"), Symbol("k")), Integer(-1))); diff(v_expr, v_X)})), "(-1/Determinant(X))*X.T^(-1)")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); ((-inverse(t(v_X))) / Determinant(v_X))})), "(-1/Determinant(X))*X.T^(-1)")
})

test_that("test_derivatives: test_issue_28838_det_in_non_matrix_expr matches SymPy", {
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_expr <- (Symbol("k")*Determinant(MatrixSymbol("X", Symbol("k"), Symbol("k")))); diff(v_expr, v_X)})), "(k*Determinant(X))*X.T^(-1)")
  expect_equal(.bat_fmt(local({v_k <- Symbol("k"); v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); ((v_k * Determinant(v_X)) * (t(v_X) ^ (-1)))})), "(k*Determinant(X))*X.T^(-1)")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_expr <- (Determinant(MatrixSymbol("X", Symbol("k"), Symbol("k"))))^(Integer(-1)); diff(v_expr, v_X)})), "(-1/Determinant(X))*X.T^(-1)")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); (((-1) / Determinant(v_X)) * (t(v_X) ^ (-1)))})), "(-1/Determinant(X))*X.T^(-1)")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_expr <- (Symbol("k")*(Determinant(MatrixSymbol("X", Symbol("k"), Symbol("k"))))^(Integer(-1))); diff(v_expr, v_X)})), "(-k/Determinant(X))*X.T^(-1)")
  expect_equal(.bat_fmt(local({v_k <- Symbol("k"); v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); (((-v_k) / Determinant(v_X)) * (t(v_X) ^ (-1)))})), "(-k/Determinant(X))*X.T^(-1)")
})

test_that("test_derivatives: test_matrix_derivative_with_inverse matches SymPy", {
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_expr <- MatMul(Transpose(MatrixSymbol("a", Symbol("k"), Integer(1))), Inverse(MatrixSymbol("X", Symbol("k"), Symbol("k"))), MatrixSymbol("b", Symbol("k"), Integer(1))); diff(v_expr, v_X)})), "-X.T^(-1)*a*b.T*X.T^(-1)")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_a <- MatrixSymbol("a", Symbol("k"), Integer(1)); v_b <- MatrixSymbol("b", Symbol("k"), Integer(1)); ((((-t(Inverse(v_X))) * v_a) * t(v_b)) * t(Inverse(v_X)))})), "-X.T^(-1)*a*b.T*X.T^(-1)")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_expr <- Determinant(Inverse(MatrixSymbol("X", Symbol("k"), Symbol("k")))); diff(v_expr, v_X)})), "(-1/Determinant(X))*X.T^(-1)")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); ((-Determinant(inverse(v_X))) * t(inverse(v_X)))})), "(-1/Determinant(X))*X.T^(-1)")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_expr <- Trace(MatMul(MatrixSymbol("A", Symbol("k"), Symbol("k")), Inverse(MatrixSymbol("X", Symbol("k"), Symbol("k"))), MatrixSymbol("B", Symbol("k"), Symbol("k")))); diff(v_expr, v_X)})), "-X.T^(-1)*A.T*B.T*X.T^(-1)")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_B <- MatrixSymbol("B", Symbol("k"), Symbol("k")); v_A <- MatrixSymbol("A", Symbol("k"), Symbol("k")); (-t(((((v_X ^ (-1)) * v_B) * v_A) * (v_X ^ (-1)))))})), "-X.T^(-1)*A.T*B.T*X.T^(-1)")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_expr <- Trace(Inverse(MatAdd(MatrixSymbol("A", Symbol("k"), Symbol("k")), MatrixSymbol("X", Symbol("k"), Symbol("k"))))); diff(v_expr, v_X)})), "-(A.T + X.T)^(-2)")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_A <- MatrixSymbol("A", Symbol("k"), Symbol("k")); (-(t(Inverse((v_X + v_A))) ^ 2))})), "-(A.T + X.T)^(-2)")
})

test_that("test_derivatives: test_matrix_derivative_vectors_and_scalars matches SymPy", {
  expect_equal(.bat_fmt(local({v_x <- MatrixSymbol("x", Symbol("k"), Integer(1)); diff(v_x, v_x)})), "I")
  expect_equal(.bat_fmt(local({v_k <- Symbol("k"); Identity(v_k)})), "I")
  expect_equal(.bat_fmt(local({v_x <- MatrixSymbol("x", Symbol("k"), Integer(1)); v_m <- Symbol("m"); v_i <- Symbol("i"); doit(diff(pyitem(v_x, v_i, 0), pyitem(v_x, v_m, 0)))})), "KroneckerDelta(i, m, (0, k - 1))")
  expect_equal(.bat_fmt(local({v_m <- Symbol("m"); v_i <- Symbol("i"); v_k <- Symbol("k"); KroneckerDelta(v_m, v_i, list(0, (v_k - 1)))})), "KroneckerDelta(i, m, (0, k - 1))")
  expect_equal(.bat_fmt(local({v_x <- MatrixSymbol("x", Symbol("k"), Integer(1)); diff(t(v_x), v_x)})), "I")
  expect_equal(.bat_fmt(local({v_x <- MatrixSymbol("x", Symbol("k"), Integer(1)); v_expr <- MatMul(Transpose(MatrixSymbol("x", Symbol("k"), Integer(1))), MatrixSymbol("a", Symbol("k"), Integer(1))); diff(v_expr, v_x)})), "a")
  expect_equal(.bat_fmt(local({v_a <- MatrixSymbol("a", Symbol("k"), Integer(1)); v_a})), "a")
  expect_equal(.bat_fmt(local({v_x <- MatrixSymbol("x", Symbol("k"), Integer(1)); v_m <- Symbol("m"); v_expr <- MatMul(Transpose(MatrixSymbol("x", Symbol("k"), Integer(1))), MatrixSymbol("a", Symbol("k"), Integer(1))); doit(diff(pyitem(v_expr, 0, 0), pyitem(v_x, v_m, 0)))})), "a[m, 0]")
  expect_equal(.bat_fmt(local({v_a <- MatrixSymbol("a", Symbol("k"), Integer(1)); v_m <- Symbol("m"); pyitem(v_a, v_m, 0)})), "a[m, 0]")
  expect_equal(.bat_fmt(local({v_x <- MatrixSymbol("x", Symbol("k"), Integer(1)); v_expr <- MatMul(Transpose(MatrixSymbol("a", Symbol("k"), Integer(1))), MatrixSymbol("x", Symbol("k"), Integer(1))); diff(v_expr, v_x)})), "a")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_expr <- MatMul(Transpose(MatrixSymbol("a", Symbol("k"), Integer(1))), MatrixSymbol("X", Symbol("k"), Symbol("k")), MatrixSymbol("b", Symbol("k"), Integer(1))); diff(v_expr, v_X)})), "a*b.T")
  expect_equal(.bat_fmt(local({v_a <- MatrixSymbol("a", Symbol("k"), Integer(1)); v_b <- MatrixSymbol("b", Symbol("k"), Integer(1)); (v_a * t(v_b))})), "a*b.T")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_expr <- MatMul(Transpose(MatrixSymbol("a", Symbol("k"), Integer(1))), Transpose(MatrixSymbol("X", Symbol("k"), Symbol("k"))), MatrixSymbol("b", Symbol("k"), Integer(1))); diff(v_expr, v_X)})), "b*a.T")
  expect_equal(.bat_fmt(local({v_b <- MatrixSymbol("b", Symbol("k"), Integer(1)); v_a <- MatrixSymbol("a", Symbol("k"), Integer(1)); (v_b * t(v_a))})), "b*a.T")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_expr <- MatMul(Transpose(MatrixSymbol("a", Symbol("k"), Integer(1))), MatrixSymbol("X", Symbol("k"), Symbol("k")), MatrixSymbol("a", Symbol("k"), Integer(1))); diff(v_expr, v_X)})), "a*a.T")
  expect_equal(.bat_fmt(local({v_a <- MatrixSymbol("a", Symbol("k"), Integer(1)); (v_a * t(v_a))})), "a*a.T")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_expr <- MatMul(Transpose(MatrixSymbol("a", Symbol("k"), Integer(1))), Transpose(MatrixSymbol("X", Symbol("k"), Symbol("k"))), MatrixSymbol("a", Symbol("k"), Integer(1))); diff(v_expr, v_X)})), "a*a.T")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_expr <- MatMul(Transpose(MatrixSymbol("b", Symbol("k"), Integer(1))), Transpose(MatrixSymbol("X", Symbol("k"), Symbol("k"))), MatrixSymbol("X", Symbol("k"), Symbol("k")), MatrixSymbol("c", Symbol("k"), Integer(1))); diff(v_expr, v_X)})), "X*b*c.T + X*c*b.T")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_b <- MatrixSymbol("b", Symbol("k"), Integer(1)); v_c <- MatrixSymbol("c", Symbol("k"), Integer(1)); (((v_X * v_b) * t(v_c)) + ((v_X * v_c) * t(v_b)))})), "X*b*c.T + X*c*b.T")
  expect_equal(.bat_fmt(local({v_x <- MatrixSymbol("x", Symbol("k"), Integer(1)); v_expr <- MatMul(MatAdd(MatMul(Transpose(MatrixSymbol("x", Symbol("k"), Integer(1))), Transpose(MatrixSymbol("B", Symbol("k"), Symbol("k")))), Transpose(MatrixSymbol("b", Symbol("k"), Integer(1)))), MatrixSymbol("C", Symbol("k"), Symbol("k")), MatAdd(MatMul(MatrixSymbol("D", Symbol("k"), Symbol("k")), MatrixSymbol("x", Symbol("k"), Integer(1))), MatrixSymbol("d", Symbol("k"), Integer(1)))); diff(v_expr, v_x)})), "B.T*C*(D*x + d) + D.T*C.T*(B*x + b)")
  expect_equal(.bat_fmt(local({v_B <- MatrixSymbol("B", Symbol("k"), Symbol("k")); v_C <- MatrixSymbol("C", Symbol("k"), Symbol("k")); v_D <- MatrixSymbol("D", Symbol("k"), Symbol("k")); v_x <- MatrixSymbol("x", Symbol("k"), Integer(1)); v_d <- MatrixSymbol("d", Symbol("k"), Integer(1)); v_b <- MatrixSymbol("b", Symbol("k"), Integer(1)); (((t(v_B) * v_C) * ((v_D * v_x) + v_d)) + ((t(v_D) * t(v_C)) * ((v_B * v_x) + v_b)))})), "B.T*C*(D*x + d) + D.T*C.T*(B*x + b)")
  expect_equal(.bat_fmt(local({v_x <- MatrixSymbol("x", Symbol("k"), Integer(1)); v_expr <- MatMul(Transpose(MatrixSymbol("x", Symbol("k"), Integer(1))), MatrixSymbol("B", Symbol("k"), Symbol("k")), MatrixSymbol("x", Symbol("k"), Integer(1))); diff(v_expr, v_x)})), "B*x + B.T*x")
  expect_equal(.bat_fmt(local({v_B <- MatrixSymbol("B", Symbol("k"), Symbol("k")); v_x <- MatrixSymbol("x", Symbol("k"), Integer(1)); ((v_B * v_x) + (t(v_B) * v_x))})), "B*x + B.T*x")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_expr <- MatMul(Transpose(MatrixSymbol("b", Symbol("k"), Integer(1))), Transpose(MatrixSymbol("X", Symbol("k"), Symbol("k"))), MatrixSymbol("D", Symbol("k"), Symbol("k")), MatrixSymbol("X", Symbol("k"), Symbol("k")), MatrixSymbol("c", Symbol("k"), Integer(1))); diff(v_expr, v_X)})), "D*X*c*b.T + D.T*X*b*c.T")
  expect_equal(.bat_fmt(local({v_D <- MatrixSymbol("D", Symbol("k"), Symbol("k")); v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_b <- MatrixSymbol("b", Symbol("k"), Integer(1)); v_c <- MatrixSymbol("c", Symbol("k"), Integer(1)); ((((t(v_D) * v_X) * v_b) * t(v_c)) + (((v_D * v_X) * v_c) * t(v_b)))})), "D*X*c*b.T + D.T*X*b*c.T")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_expr <- MatMul(MatAdd(MatMul(Transpose(MatrixSymbol("b", Symbol("k"), Integer(1))), Transpose(MatrixSymbol("X", Symbol("k"), Symbol("k")))), Transpose(MatrixSymbol("c", Symbol("k"), Integer(1)))), MatrixSymbol("D", Symbol("k"), Symbol("k")), MatAdd(MatMul(MatrixSymbol("X", Symbol("k"), Symbol("k")), MatrixSymbol("b", Symbol("k"), Integer(1))), MatrixSymbol("c", Symbol("k"), Integer(1)))); diff(v_expr, v_X)})), "D*(X*b + c)*b.T + D.T*(X*b + c)*b.T")
  expect_equal(.bat_fmt(local({v_D <- MatrixSymbol("D", Symbol("k"), Symbol("k")); v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_b <- MatrixSymbol("b", Symbol("k"), Integer(1)); v_c <- MatrixSymbol("c", Symbol("k"), Integer(1)); (((v_D * ((v_X * v_b) + v_c)) * t(v_b)) + ((t(v_D) * ((v_X * v_b) + v_c)) * t(v_b)))})), "D*(X*b + c)*b.T + D.T*(X*b + c)*b.T")
  expect_equal(.bat_fmt(local({v_x <- MatrixSymbol("x", Symbol("k"), Integer(1)); v_expr <- MatMul(MatrixSymbol("x", Symbol("k"), Integer(1)), Transpose(MatrixSymbol("x", Symbol("k"), Integer(1))), MatrixSymbol("x", Symbol("k"), Integer(1))); diff(v_expr, v_x)})), "KroneckerProduct(I, x.T*x) + 2*x*x.T")
  expect_equal(.bat_fmt(local({v_I <- Identity(Symbol("k")); v_x <- MatrixSymbol("x", Symbol("k"), Integer(1)); (KroneckerProduct(v_I, (t(v_x) * v_x)) + ((2 * v_x) * t(v_x)))})), "KroneckerProduct(I, x.T*x) + 2*x*x.T")
})

test_that("test_derivatives: test_matrix_derivatives_of_traces matches SymPy", {
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("k"), Symbol("k")); v_expr <- MatMul(Trace(MatrixSymbol("A", Symbol("k"), Symbol("k"))), MatrixSymbol("A", Symbol("k"), Symbol("k"))); diff(v_expr, v_A)})), "ArrayAdd(ArrayTensorProduct(I, A), PermuteDims(ArrayTensorProduct(Trace(A), I, I), (3)(1 2)))")
  expect_equal(.bat_fmt(local({v_I <- Identity(Symbol("k")); v_A <- MatrixSymbol("A", Symbol("k"), Symbol("k")); v__perm1 <- Permutation(c(0, 2, 1, 3)); ArrayAdd(ArrayTensorProduct(v_I, v_A), PermuteDims(ArrayTensorProduct(Trace(v_A), v_I, v_I), v__perm1))})), "ArrayAdd(ArrayTensorProduct(I, A), PermuteDims(ArrayTensorProduct(Trace(A), I, I), (3)(1 2)))")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("k"), Symbol("k")); v_m <- Symbol("m"); v_n <- Symbol("n"); v_expr <- MatMul(Trace(MatrixSymbol("A", Symbol("k"), Symbol("k"))), MatrixSymbol("A", Symbol("k"), Symbol("k"))); v_i <- Symbol("i"); v_j <- Symbol("j"); doit(diff(pyitem(v_expr, v_i, v_j), pyitem(v_A, v_m, v_n)))})), "KroneckerDelta(i, m, (0, k - 1))*KroneckerDelta(j, n, (0, k - 1))*Trace(A) + KroneckerDelta(m, n, (0, k - 1))*A[i, j]")
  expect_equal(.bat_fmt(local({v_i <- Symbol("i"); v_m <- Symbol("m"); v_k <- Symbol("k"); v_j <- Symbol("j"); v_n <- Symbol("n"); v_A <- MatrixSymbol("A", Symbol("k"), Symbol("k")); (((KroneckerDelta(v_i, v_m, list(0, (v_k - 1))) * KroneckerDelta(v_j, v_n, list(0, (v_k - 1)))) * Trace(v_A)) + (KroneckerDelta(v_m, v_n, list(0, (v_k - 1))) * pyitem(v_A, v_i, v_j)))})), "KroneckerDelta(i, m, (0, k - 1))*KroneckerDelta(j, n, (0, k - 1))*Trace(A) + KroneckerDelta(m, n, (0, k - 1))*A[i, j]")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_expr <- Trace(MatrixSymbol("X", Symbol("k"), Symbol("k"))); diff(v_expr, v_X)})), "I")
  expect_equal(.bat_fmt(local({v_m <- Symbol("m"); v_n <- Symbol("n"); v_k <- Symbol("k"); KroneckerDelta(v_m, v_n, list(0, (v_k - 1)))})), "KroneckerDelta(m, n, (0, k - 1))")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_expr <- Trace(MatMul(MatrixSymbol("X", Symbol("k"), Symbol("k")), MatrixSymbol("A", Symbol("k"), Symbol("k")))); diff(v_expr, v_X)})), "A.T")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("k"), Symbol("k")); t(v_A)})), "A.T")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("k"), Symbol("k")); v_n <- Symbol("n"); v_m <- Symbol("m"); pyitem(v_A, v_n, v_m)})), "A[n, m]")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_expr <- Trace(MatMul(MatrixSymbol("A", Symbol("k"), Symbol("k")), MatrixSymbol("X", Symbol("k"), Symbol("k")), MatrixSymbol("B", Symbol("k"), Symbol("k")))); diff(v_expr, v_X)})), "A.T*B.T")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("k"), Symbol("k")); v_B <- MatrixSymbol("B", Symbol("k"), Symbol("k")); (t(v_A) * t(v_B))})), "A.T*B.T")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("k"), Symbol("k")); v_B <- MatrixSymbol("B", Symbol("k"), Symbol("k")); v_m <- Symbol("m"); v_n <- Symbol("n"); pyitem((t(v_A) * t(v_B)), v_m, v_n)})), "Sum(A[_i_1, m]*B[n, _i_1], (_i_1, 0, k - 1))")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_expr <- Trace(MatMul(MatrixSymbol("A", Symbol("k"), Symbol("k")), Transpose(MatrixSymbol("X", Symbol("k"), Symbol("k"))), MatrixSymbol("B", Symbol("k"), Symbol("k")))); diff(v_expr, v_X)})), "B*A")
  expect_equal(.bat_fmt(local({v_B <- MatrixSymbol("B", Symbol("k"), Symbol("k")); v_A <- MatrixSymbol("A", Symbol("k"), Symbol("k")); (v_B * v_A)})), "B*A")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_expr <- Trace(MatMul(Transpose(MatrixSymbol("X", Symbol("k"), Symbol("k"))), MatrixSymbol("A", Symbol("k"), Symbol("k")))); diff(v_expr, v_X)})), "A")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("k"), Symbol("k")); v_A})), "A")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_expr <- Trace(MatMul(MatrixSymbol("A", Symbol("k"), Symbol("k")), Transpose(MatrixSymbol("X", Symbol("k"), Symbol("k"))))); diff(v_expr, v_X)})), "A")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_expr <- Trace(MatPow(MatrixSymbol("X", Symbol("k"), Symbol("k")), Integer(2))); diff(v_expr, v_X)})), "2*X.T")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); (2 * t(v_X))})), "2*X.T")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_expr <- Trace(MatMul(MatPow(MatrixSymbol("X", Symbol("k"), Symbol("k")), Integer(2)), MatrixSymbol("B", Symbol("k"), Symbol("k")))); diff(v_expr, v_X)})), "B.T*X.T + X.T*B.T")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_B <- MatrixSymbol("B", Symbol("k"), Symbol("k")); t(((v_X * v_B) + (v_B * v_X)))})), "B.T*X.T + X.T*B.T")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_expr <- Trace(MatMul(MatrixSymbol("X", Symbol("k"), Symbol("k")), MatrixSymbol("X", Symbol("k"), Symbol("k")), MatrixSymbol("B", Symbol("k"), Symbol("k")))); diff(v_expr, v_X)})), "B.T*X.T + X.T*B.T")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_expr <- Trace(MatMul(Transpose(MatrixSymbol("X", Symbol("k"), Symbol("k"))), MatrixSymbol("B", Symbol("k"), Symbol("k")), MatrixSymbol("X", Symbol("k"), Symbol("k")))); diff(v_expr, v_X)})), "B*X + B.T*X")
  expect_equal(.bat_fmt(local({v_B <- MatrixSymbol("B", Symbol("k"), Symbol("k")); v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); ((v_B * v_X) + (t(v_B) * v_X))})), "B*X + B.T*X")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_expr <- Trace(MatMul(MatrixSymbol("B", Symbol("k"), Symbol("k")), MatrixSymbol("X", Symbol("k"), Symbol("k")), Transpose(MatrixSymbol("X", Symbol("k"), Symbol("k"))))); diff(v_expr, v_X)})), "B*X + B.T*X")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_expr <- Trace(MatMul(MatrixSymbol("X", Symbol("k"), Symbol("k")), Transpose(MatrixSymbol("X", Symbol("k"), Symbol("k"))), MatrixSymbol("B", Symbol("k"), Symbol("k")))); diff(v_expr, v_X)})), "B*X + B.T*X")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_expr <- Trace(MatMul(MatrixSymbol("X", Symbol("k"), Symbol("k")), MatrixSymbol("B", Symbol("k"), Symbol("k")), Transpose(MatrixSymbol("X", Symbol("k"), Symbol("k"))))); diff(v_expr, v_X)})), "X*B + X*B.T")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_B <- MatrixSymbol("B", Symbol("k"), Symbol("k")); ((v_X * t(v_B)) + (v_X * v_B))})), "X*B + X*B.T")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_expr <- Trace(MatMul(MatrixSymbol("B", Symbol("k"), Symbol("k")), Transpose(MatrixSymbol("X", Symbol("k"), Symbol("k"))), MatrixSymbol("X", Symbol("k"), Symbol("k")))); diff(v_expr, v_X)})), "X*B + X*B.T")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_expr <- Trace(MatMul(Transpose(MatrixSymbol("X", Symbol("k"), Symbol("k"))), MatrixSymbol("X", Symbol("k"), Symbol("k")), MatrixSymbol("B", Symbol("k"), Symbol("k")))); diff(v_expr, v_X)})), "X*B + X*B.T")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_expr <- Trace(MatMul(MatrixSymbol("A", Symbol("k"), Symbol("k")), MatrixSymbol("X", Symbol("k"), Symbol("k")), MatrixSymbol("B", Symbol("k"), Symbol("k")), MatrixSymbol("X", Symbol("k"), Symbol("k")))); diff(v_expr, v_X)})), "A.T*X.T*B.T + B.T*X.T*A.T")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("k"), Symbol("k")); v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_B <- MatrixSymbol("B", Symbol("k"), Symbol("k")); (((t(v_A) * t(v_X)) * t(v_B)) + ((t(v_B) * t(v_X)) * t(v_A)))})), "A.T*X.T*B.T + B.T*X.T*A.T")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_expr <- Trace(MatMul(Transpose(MatrixSymbol("X", Symbol("k"), Symbol("k"))), MatrixSymbol("X", Symbol("k"), Symbol("k")))); diff(v_expr, v_X)})), "2*X")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); (2 * v_X)})), "2*X")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_expr <- Trace(MatMul(MatrixSymbol("X", Symbol("k"), Symbol("k")), Transpose(MatrixSymbol("X", Symbol("k"), Symbol("k"))))); diff(v_expr, v_X)})), "2*X")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_expr <- Trace(MatMul(Transpose(MatrixSymbol("B", Symbol("k"), Symbol("k"))), Transpose(MatrixSymbol("X", Symbol("k"), Symbol("k"))), MatrixSymbol("C", Symbol("k"), Symbol("k")), MatrixSymbol("X", Symbol("k"), Symbol("k")), MatrixSymbol("B", Symbol("k"), Symbol("k")))); diff(v_expr, v_X)})), "C*X*B*B.T + C.T*X*B*B.T")
  expect_equal(.bat_fmt(local({v_C <- MatrixSymbol("C", Symbol("k"), Symbol("k")); v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_B <- MatrixSymbol("B", Symbol("k"), Symbol("k")); ((((t(v_C) * v_X) * v_B) * t(v_B)) + (((v_C * v_X) * v_B) * t(v_B)))})), "C*X*B*B.T + C.T*X*B*B.T")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_expr <- Trace(MatMul(Transpose(MatrixSymbol("X", Symbol("k"), Symbol("k"))), MatrixSymbol("B", Symbol("k"), Symbol("k")), MatrixSymbol("X", Symbol("k"), Symbol("k")), MatrixSymbol("C", Symbol("k"), Symbol("k")))); diff(v_expr, v_X)})), "B*X*C + B.T*X*C.T")
  expect_equal(.bat_fmt(local({v_B <- MatrixSymbol("B", Symbol("k"), Symbol("k")); v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_C <- MatrixSymbol("C", Symbol("k"), Symbol("k")); (((v_B * v_X) * v_C) + ((t(v_B) * v_X) * t(v_C)))})), "B*X*C + B.T*X*C.T")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_expr <- Trace(MatMul(MatrixSymbol("A", Symbol("k"), Symbol("k")), MatrixSymbol("X", Symbol("k"), Symbol("k")), MatrixSymbol("B", Symbol("k"), Symbol("k")), Transpose(MatrixSymbol("X", Symbol("k"), Symbol("k"))), MatrixSymbol("C", Symbol("k"), Symbol("k")))); diff(v_expr, v_X)})), "C*A*X*B + A.T*C.T*X*B.T")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("k"), Symbol("k")); v_C <- MatrixSymbol("C", Symbol("k"), Symbol("k")); v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_B <- MatrixSymbol("B", Symbol("k"), Symbol("k")); ((((t(v_A) * t(v_C)) * v_X) * t(v_B)) + (((v_C * v_A) * v_X) * v_B))})), "C*A*X*B + A.T*C.T*X*B.T")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_expr <- Trace(MatMul(MatAdd(MatMul(MatrixSymbol("A", Symbol("k"), Symbol("k")), MatrixSymbol("X", Symbol("k"), Symbol("k")), MatrixSymbol("B", Symbol("k"), Symbol("k"))), MatrixSymbol("C", Symbol("k"), Symbol("k"))), MatAdd(MatMul(Transpose(MatrixSymbol("B", Symbol("k"), Symbol("k"))), Transpose(MatrixSymbol("X", Symbol("k"), Symbol("k"))), Transpose(MatrixSymbol("A", Symbol("k"), Symbol("k")))), Transpose(MatrixSymbol("C", Symbol("k"), Symbol("k")))))); diff(v_expr, v_X)})), "2*A.T*(A*X*B + C)*B.T")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("k"), Symbol("k")); v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_B <- MatrixSymbol("B", Symbol("k"), Symbol("k")); v_C <- MatrixSymbol("C", Symbol("k"), Symbol("k")); (((2 * t(v_A)) * (((v_A * v_X) * v_B) + v_C)) * t(v_B))})), "2*A.T*(A*X*B + C)*B.T")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_expr <- Trace(MatMul(Transpose(MatrixSymbol("B", Symbol("k"), Symbol("k"))), Transpose(MatrixSymbol("X", Symbol("k"), Symbol("k"))), MatrixSymbol("C", Symbol("k"), Symbol("k")), MatrixSymbol("X", Symbol("k"), Symbol("k")), Transpose(MatrixSymbol("X", Symbol("k"), Symbol("k"))), MatrixSymbol("C", Symbol("k"), Symbol("k")), MatrixSymbol("X", Symbol("k"), Symbol("k")), MatrixSymbol("B", Symbol("k"), Symbol("k")))); diff(v_expr, v_X)})), "C*X*B*B.T*X.T*C*X + C*X*X.T*C*X*B*B.T + C.T*X*B*B.T*X.T*C.T*X + C.T*X*X.T*C.T*X*B*B.T")
  expect_equal(.bat_fmt(local({v_C <- MatrixSymbol("C", Symbol("k"), Symbol("k")); v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_B <- MatrixSymbol("B", Symbol("k"), Symbol("k")); (((((((((v_C * v_X) * t(v_X)) * v_C) * v_X) * v_B) * t(v_B)) + ((((((t(v_C) * v_X) * v_B) * t(v_B)) * t(v_X)) * t(v_C)) * v_X)) + ((((((v_C * v_X) * v_B) * t(v_B)) * t(v_X)) * v_C) * v_X)) + ((((((t(v_C) * v_X) * t(v_X)) * t(v_C)) * v_X) * v_B) * t(v_B)))})), "C*X*B*B.T*X.T*C*X + C*X*X.T*C*X*B*B.T + C.T*X*B*B.T*X.T*C.T*X + C.T*X*X.T*C.T*X*B*B.T")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_A <- MatrixSymbol("A", Symbol("k"), Symbol("k")); v_B <- MatrixSymbol("B", Symbol("k"), Symbol("k")); ((((-t(Inverse(v_X))) * t(v_A)) * t(v_B)) * t(Inverse(v_X)))})), "-X.T^(-1)*A.T*B.T*X.T^(-1)")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_expr <- Trace(MatMul(Inverse(MatrixSymbol("X", Symbol("k"), Symbol("k"))), Inverse(MatrixSymbol("C", Symbol("k"), Symbol("k"))), Inverse(Transpose(MatrixSymbol("X", Symbol("k"), Symbol("k")))), MatrixSymbol("A", Symbol("k"), Symbol("k")))); diff(v_expr, v_X)})), "-X.T^(-1)*A*X^(-1)*C^(-1)*X.T^(-1) - X.T^(-1)*A.T*X^(-1)*C.T^(-1)*X.T^(-1)")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_A <- MatrixSymbol("A", Symbol("k"), Symbol("k")); v_C <- MatrixSymbol("C", Symbol("k"), Symbol("k")); ((((((-t(inverse(v_X))) * t(v_A)) * inverse(v_X)) * t(inverse(v_C))) * t(inverse(v_X))) - ((((t(inverse(v_X)) * v_A) * inverse(v_X)) * inverse(v_C)) * t(inverse(v_X))))})), "-X.T^(-1)*A*X^(-1)*C^(-1)*X.T^(-1) - X.T^(-1)*A.T*X^(-1)*C.T^(-1)*X.T^(-1)")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_expr <- Trace(MatMul(Inverse(MatrixSymbol("X", Symbol("k"), Symbol("k"))), Inverse(MatrixSymbol("C", Symbol("k"), Symbol("k"))), MatrixSymbol("B", Symbol("k"), Symbol("k")), MatrixSymbol("X", Symbol("k"), Symbol("k")))); diff(v_expr, v_X)})), "0")
  expect_equal(.bat_fmt(local({v_C <- MatrixSymbol("C", Symbol("k"), Symbol("k")); v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_B <- MatrixSymbol("B", Symbol("k"), Symbol("k")); (((((((((-2) * v_C) * v_X) * inverse(((t(v_X) * v_C) * v_X))) * t(v_X)) * v_B) * v_X) * inverse(((t(v_X) * v_C) * v_X))) + (((2 * v_B) * v_X) * inverse(((t(v_X) * v_C) * v_X))))})), "0")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_expr <- Trace(MatMul(Inverse(MatAdd(MatMul(Transpose(MatrixSymbol("X", Symbol("k"), Symbol("k"))), MatrixSymbol("C", Symbol("k"), Symbol("k")), MatrixSymbol("X", Symbol("k"), Symbol("k"))), MatrixSymbol("A", Symbol("k"), Symbol("k")))), Transpose(MatrixSymbol("X", Symbol("k"), Symbol("k"))), MatrixSymbol("B", Symbol("k"), Symbol("k")), MatrixSymbol("X", Symbol("k"), Symbol("k")))); diff(v_expr, v_X)})), "B*X*(A + X.T*C*X)^(-1) - C*X*(A + X.T*C*X)^(-1)*X.T*B*X*(A + X.T*C*X)^(-1) + B.T*X*(A.T + X.T*C.T*X)^(-1) - C.T*X*(A.T + X.T*C.T*X)^(-1)*X.T*B.T*X*(A.T + X.T*C.T*X)^(-1)")
  expect_equal(.bat_fmt(local({v_B <- MatrixSymbol("B", Symbol("k"), Symbol("k")); v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_A <- MatrixSymbol("A", Symbol("k"), Symbol("k")); v_C <- MatrixSymbol("C", Symbol("k"), Symbol("k")); (((((v_B * v_X) * Inverse((v_A + ((t(v_X) * v_C) * v_X)))) - ((((((v_C * v_X) * Inverse((v_A + ((t(v_X) * v_C) * v_X)))) * t(v_X)) * v_B) * v_X) * Inverse((v_A + ((t(v_X) * v_C) * v_X))))) - ((((((t(v_C) * v_X) * Inverse((t(v_A) + (t((v_C * v_X)) * v_X)))) * t(v_X)) * t(v_B)) * v_X) * Inverse((t(v_A) + (t((v_C * v_X)) * v_X))))) + ((t(v_B) * v_X) * Inverse((t(v_A) + (t((v_C * v_X)) * v_X)))))})), "B*X*(A + X.T*C*X)^(-1) - C*X*(A + X.T*C*X)^(-1)*X.T*B*X*(A + X.T*C*X)^(-1) + B.T*X*(A.T + X.T*C.T*X)^(-1) - C.T*X*(A.T + X.T*C.T*X)^(-1)*X.T*B.T*X*(A.T + X.T*C.T*X)^(-1)")
})

test_that("test_derivatives: test_derivatives_of_complicated_matrix_expr matches SymPy", {
  expect_equal(.bat_fmt(local({v_result <- MatAdd(MatMul(Transpose(MatrixSymbol("A", Symbol("k"), Symbol("k"))), MatrixSymbol("a", Symbol("k"), Integer(1)), Transpose(MatrixSymbol("b", Symbol("k"), Integer(1))), MatAdd(MatMul(Transpose(MatrixSymbol("A", Symbol("k"), Symbol("k"))), Transpose(MatrixSymbol("X", Symbol("k"), Symbol("k")))), MatMul(Transpose(MatrixSymbol("B", Symbol("k"), Symbol("k"))), MatrixSymbol("X", Symbol("k"), Symbol("k"))))), MatMul(MatAdd(MatMul(MatrixSymbol("B", Symbol("k"), Symbol("k")), MatAdd(MatMul(MatrixSymbol("X", Symbol("k"), Symbol("k")), Transpose(MatrixSymbol("D", Symbol("k"), Symbol("k")))), MatMul(Integer(-3), MatrixSymbol("X", Symbol("k"), Symbol("k")), MatrixSymbol("D", Symbol("k"), Symbol("k"))), MatMul(MatrixSymbol("B", Symbol("k"), Symbol("k")), MatrixSymbol("A", Symbol("k"), Symbol("k")), MatrixSymbol("X", Symbol("k"), Symbol("k")), Transpose(MatrixSymbol("A", Symbol("k"), Symbol("k")))))), MatMul(MatrixSymbol("a", Symbol("k"), Integer(1)), Transpose(MatrixSymbol("b", Symbol("k"), Integer(1))), MatAdd(MatMul(MatrixSymbol("X", Symbol("k"), Symbol("k")), MatrixSymbol("D", Symbol("k"), Symbol("k")), Transpose(MatrixSymbol("X", Symbol("k"), Symbol("k")))), MatMul(Integer(-1), Transpose(MatrixSymbol("X", Symbol("k"), Symbol("k"))), Transpose(MatrixSymbol("C", Symbol("k"), Symbol("k"))), MatrixSymbol("A", Symbol("k"), Symbol("k"))), MatMul(MatrixSymbol("X", Symbol("k"), Symbol("k")), MatAdd(MatMul(MatrixSymbol("A", Symbol("k"), Symbol("k")), MatrixSymbol("X", Symbol("k"), Symbol("k"))), MatMul(Transpose(MatrixSymbol("X", Symbol("k"), Symbol("k"))), MatrixSymbol("B", Symbol("k"), Symbol("k")))), MatrixSymbol("D", Symbol("k"), Symbol("k")), MatrixSymbol("B", Symbol("k"), Symbol("k")))), MatrixSymbol("B", Symbol("k"), Symbol("k")))), MatrixSymbol("B", Symbol("k"), Symbol("k")), MatrixSymbol("b", Symbol("k"), Integer(1)), Transpose(MatrixSymbol("a", Symbol("k"), Integer(1))), Transpose(MatrixSymbol("B", Symbol("k"), Symbol("k")))), MatMul(MatrixSymbol("B", Symbol("k"), Symbol("k")), MatrixSymbol("b", Symbol("k"), Integer(1)), Transpose(MatrixSymbol("a", Symbol("k"), Integer(1))), MatrixSymbol("A", Symbol("k"), Symbol("k")), MatrixSymbol("X", Symbol("k"), Symbol("k"))), MatMul(Transpose(MatrixSymbol("X", Symbol("k"), Symbol("k"))), Transpose(MatrixSymbol("A", Symbol("k"), Symbol("k"))), MatrixSymbol("a", Symbol("k"), Integer(1)), Transpose(MatrixSymbol("b", Symbol("k"), Integer(1))), Transpose(MatrixSymbol("A", Symbol("k"), Symbol("k")))), MatMul(MatrixSymbol("a", Symbol("k"), Integer(1)), Transpose(MatrixSymbol("b", Symbol("k"), Integer(1))), MatAdd(MatMul(Integer(42), MatrixSymbol("X", Symbol("k"), Symbol("k"))), MatMul(Integer(42), Transpose(MatrixSymbol("X", Symbol("k"), Symbol("k"))))), MatrixSymbol("A", Symbol("k"), Symbol("k")), MatrixSymbol("X", Symbol("k"), Symbol("k")), Transpose(MatrixSymbol("B", Symbol("k"), Symbol("k")))), MatMul(Transpose(MatrixSymbol("A", Symbol("k"), Symbol("k"))), MatAdd(MatMul(Integer(42), MatrixSymbol("X", Symbol("k"), Symbol("k"))), MatMul(Integer(42), Transpose(MatrixSymbol("X", Symbol("k"), Symbol("k"))))), MatrixSymbol("b", Symbol("k"), Integer(1)), Transpose(MatrixSymbol("a", Symbol("k"), Integer(1))), MatrixSymbol("X", Symbol("k"), Symbol("k")), MatrixSymbol("B", Symbol("k"), Symbol("k"))), MatMul(Integer(42), MatrixSymbol("A", Symbol("k"), Symbol("k")), MatrixSymbol("X", Symbol("k"), Symbol("k")), Transpose(MatrixSymbol("B", Symbol("k"), Symbol("k"))), Transpose(MatrixSymbol("X", Symbol("k"), Symbol("k"))), MatrixSymbol("a", Symbol("k"), Integer(1)), Transpose(MatrixSymbol("b", Symbol("k"), Integer(1)))), MatMul(Integer(42), MatrixSymbol("b", Symbol("k"), Integer(1)), Transpose(MatrixSymbol("a", Symbol("k"), Integer(1))), MatrixSymbol("X", Symbol("k"), Symbol("k")), MatrixSymbol("B", Symbol("k"), Symbol("k")), Transpose(MatrixSymbol("X", Symbol("k"), Symbol("k"))), Transpose(MatrixSymbol("A", Symbol("k"), Symbol("k")))), MatMul(Transpose(MatrixSymbol("B", Symbol("k"), Symbol("k"))), MatrixSymbol("X", Symbol("k"), Symbol("k")), MatrixSymbol("B", Symbol("k"), Symbol("k")), MatrixSymbol("a", Symbol("k"), Integer(1)), Transpose(MatrixSymbol("b", Symbol("k"), Integer(1))), Transpose(MatrixSymbol("B", Symbol("k"), Symbol("k"))), MatrixSymbol("D", Symbol("k"), Symbol("k"))), MatMul(Integer(-3), Transpose(MatrixSymbol("B", Symbol("k"), Symbol("k"))), MatrixSymbol("X", Symbol("k"), Symbol("k")), MatrixSymbol("B", Symbol("k"), Symbol("k")), MatrixSymbol("a", Symbol("k"), Integer(1)), Transpose(MatrixSymbol("b", Symbol("k"), Integer(1))), Transpose(MatrixSymbol("B", Symbol("k"), Symbol("k"))), Transpose(MatrixSymbol("D", Symbol("k"), Symbol("k")))), MatMul(Transpose(MatrixSymbol("A", Symbol("k"), Symbol("k"))), MatPow(Transpose(MatrixSymbol("B", Symbol("k"), Symbol("k"))), Integer(2)), MatrixSymbol("X", Symbol("k"), Symbol("k")), MatrixSymbol("B", Symbol("k"), Symbol("k")), MatrixSymbol("a", Symbol("k"), Integer(1)), Transpose(MatrixSymbol("b", Symbol("k"), Integer(1))), Transpose(MatrixSymbol("B", Symbol("k"), Symbol("k"))), MatrixSymbol("A", Symbol("k"), Symbol("k"))), MatMul(MatPow(MatrixSymbol("B", Symbol("k"), Symbol("k")), Integer(2)), MatrixSymbol("b", Symbol("k"), Integer(1)), Transpose(MatrixSymbol("a", Symbol("k"), Integer(1))), Transpose(MatrixSymbol("B", Symbol("k"), Symbol("k"))), Transpose(MatrixSymbol("X", Symbol("k"), Symbol("k"))), MatrixSymbol("a", Symbol("k"), Integer(1)), Transpose(MatrixSymbol("b", Symbol("k"), Integer(1))), MatrixSymbol("X", Symbol("k"), Symbol("k")), MatrixSymbol("D", Symbol("k"), Symbol("k"))), MatMul(MatrixSymbol("b", Symbol("k"), Integer(1)), Transpose(MatrixSymbol("a", Symbol("k"), Integer(1))), MatrixSymbol("X", Symbol("k"), Symbol("k")), MatrixSymbol("B", Symbol("k"), Symbol("k")), MatrixSymbol("a", Symbol("k"), Integer(1)), Transpose(MatrixSymbol("b", Symbol("k"), Integer(1))), MatPow(Transpose(MatrixSymbol("B", Symbol("k"), Symbol("k"))), Integer(2)), MatrixSymbol("X", Symbol("k"), Symbol("k")), Transpose(MatrixSymbol("D", Symbol("k"), Symbol("k")))), MatMul(MatrixSymbol("b", Symbol("k"), Integer(1)), Transpose(MatrixSymbol("a", Symbol("k"), Integer(1))), MatrixSymbol("X", Symbol("k"), Symbol("k")), MatrixSymbol("B", Symbol("k"), Symbol("k")), MatrixSymbol("a", Symbol("k"), Integer(1)), Transpose(MatrixSymbol("b", Symbol("k"), Integer(1))), MatPow(Transpose(MatrixSymbol("B", Symbol("k"), Symbol("k"))), Integer(3)), Transpose(MatrixSymbol("D", Symbol("k"), Symbol("k"))), MatAdd(MatMul(Transpose(MatrixSymbol("B", Symbol("k"), Symbol("k"))), MatrixSymbol("X", Symbol("k"), Symbol("k"))), MatMul(Transpose(MatrixSymbol("X", Symbol("k"), Symbol("k"))), Transpose(MatrixSymbol("A", Symbol("k"), Symbol("k")))))), MatMul(Integer(-1), Transpose(MatrixSymbol("C", Symbol("k"), Symbol("k"))), MatrixSymbol("A", Symbol("k"), Symbol("k")), MatPow(MatrixSymbol("B", Symbol("k"), Symbol("k")), Integer(2)), MatrixSymbol("b", Symbol("k"), Integer(1)), Transpose(MatrixSymbol("a", Symbol("k"), Integer(1))), Transpose(MatrixSymbol("B", Symbol("k"), Symbol("k"))), Transpose(MatrixSymbol("X", Symbol("k"), Symbol("k"))), MatrixSymbol("a", Symbol("k"), Integer(1)), Transpose(MatrixSymbol("b", Symbol("k"), Integer(1)))), MatMul(MatrixSymbol("B", Symbol("k"), Symbol("k")), MatrixSymbol("D", Symbol("k"), Symbol("k")), MatPow(MatrixSymbol("B", Symbol("k"), Symbol("k")), Integer(3)), MatrixSymbol("b", Symbol("k"), Integer(1)), Transpose(MatrixSymbol("a", Symbol("k"), Integer(1))), Transpose(MatrixSymbol("B", Symbol("k"), Symbol("k"))), Transpose(MatrixSymbol("X", Symbol("k"), Symbol("k"))), MatrixSymbol("a", Symbol("k"), Integer(1)), Transpose(MatrixSymbol("b", Symbol("k"), Integer(1))), MatrixSymbol("X", Symbol("k"), Symbol("k"))), MatMul(Transpose(MatrixSymbol("A", Symbol("k"), Symbol("k"))), Transpose(MatrixSymbol("X", Symbol("k"), Symbol("k"))), MatrixSymbol("b", Symbol("k"), Integer(1)), Transpose(MatrixSymbol("a", Symbol("k"), Integer(1))), MatrixSymbol("X", Symbol("k"), Symbol("k")), MatrixSymbol("B", Symbol("k"), Symbol("k")), MatrixSymbol("a", Symbol("k"), Integer(1)), Transpose(MatrixSymbol("b", Symbol("k"), Integer(1))), MatPow(Transpose(MatrixSymbol("B", Symbol("k"), Symbol("k"))), Integer(3)), Transpose(MatrixSymbol("D", Symbol("k"), Symbol("k"))))); v_result})), "(B*(B*A*X*A.T - 3*X*D + X*D.T) + a*b.T*(X*(A*X + X.T*B)*D*B + X*D*X.T - X.T*C.T*A)*B)*B*b*a.T*B.T + B^2*b*a.T*B.T*X.T*a*b.T*X*D + 42*A*X*B.T*X.T*a*b.T + B*D*B^3*b*a.T*B.T*X.T*a*b.T*X + B*b*a.T*A*X + a*b.T*(42*X + 42*X.T)*A*X*B.T + b*a.T*X*B*a*b.T*B.T^2*X*D.T + b*a.T*X*B*a*b.T*B.T^3*D.T*(B.T*X + X.T*A.T) + 42*b*a.T*X*B*X.T*A.T + A.T*(42*X + 42*X.T)*b*a.T*X*B + A.T*B.T^2*X*B*a*b.T*B.T*A + A.T*a*b.T*(A.T*X.T + B.T*X) + A.T*X.T*b*a.T*X*B*a*b.T*B.T^3*D.T + B.T*X*B*a*b.T*B.T*D - 3*B.T*X*B*a*b.T*B.T*D.T - C.T*A*B^2*b*a.T*B.T*X.T*a*b.T + X.T*A.T*a*b.T*A.T")
})

test_that("test_derivatives: test_mixed_deriv_mixed_expressions matches SymPy", {
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("k"), Symbol("k")); v_expr <- (Integer(3)*Trace(MatrixSymbol("A", Symbol("k"), Symbol("k")))); diff(v_expr, v_A)})), "3*I")
  expect_equal(.bat_fmt(local({v_k <- Symbol("k"); (3 * Identity(v_k))})), "3*I")
  expect_equal(.bat_fmt(local({v_deriv <- ZeroMatrix(Symbol("k"), Symbol("k")); v_deriv})), "0")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("k"), Symbol("k")); v_expr <- (Trace(MatrixSymbol("A", Symbol("k"), Symbol("k"))))^(Integer(2)); diff(v_expr, v_A)})), "(2*Trace(A))*I")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("k"), Symbol("k")); v_k <- Symbol("k"); ((2 * Trace(v_A)) * Identity(v_k))})), "(2*Trace(A))*I")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("k"), Symbol("k")); v_expr <- Trace(MatMul(Trace(MatrixSymbol("A", Symbol("k"), Symbol("k"))), MatrixSymbol("A", Symbol("k"), Symbol("k")))); diff(v_expr, v_A)})), "(2*Trace(A))*I")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("k"), Symbol("k")); v_expr <- Trace(MatMul((Trace(MatrixSymbol("A", Symbol("k"), Symbol("k"))))^(Integer(2)), MatrixSymbol("A", Symbol("k"), Symbol("k")))); diff(v_expr, v_A)})), "(3*Trace(A)^2)*I")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("k"), Symbol("k")); v_k <- Symbol("k"); ((3 * (Trace(v_A) ^ 2)) * Identity(v_k))})), "(3*Trace(A)^2)*I")
})

test_that("test_derivatives: test_derivatives_matrix_norms matches SymPy", {
  expect_equal(.bat_fmt(local({v_x <- MatrixSymbol("x", Symbol("k"), Integer(1)); v_expr <- MatMul(Transpose(MatrixSymbol("x", Symbol("k"), Integer(1))), MatrixSymbol("y", Symbol("k"), Integer(1))); diff(v_expr, v_x)})), "y")
  expect_equal(.bat_fmt(local({v_y <- MatrixSymbol("y", Symbol("k"), Integer(1)); v_y})), "y")
  expect_equal(.bat_fmt(local({v_x <- MatrixSymbol("x", Symbol("k"), Integer(1)); v_m <- Symbol("m"); v_expr <- MatMul(Transpose(MatrixSymbol("x", Symbol("k"), Integer(1))), MatrixSymbol("y", Symbol("k"), Integer(1))); doit(diff(pyitem(v_expr, 0, 0), pyitem(v_x, v_m, 0)))})), "y[m, 0]")
  expect_equal(.bat_fmt(local({v_y <- MatrixSymbol("y", Symbol("k"), Integer(1)); v_m <- Symbol("m"); pyitem(v_y, v_m, 0)})), "y[m, 0]")
  expect_equal(.bat_fmt(local({v_x <- MatrixSymbol("x", Symbol("k"), Integer(1)); v_expr <- MatPow(MatMul(Transpose(MatrixSymbol("x", Symbol("k"), Integer(1))), MatrixSymbol("y", Symbol("k"), Integer(1))), Rational(1, 2)); diff(v_expr, v_x)})), "(1/2)*y*(x.T*y)^(-1/2)")
  expect_equal(.bat_fmt(local({v_y <- MatrixSymbol("y", Symbol("k"), Integer(1)); v_x <- MatrixSymbol("x", Symbol("k"), Integer(1)); (v_y / (2 * sqrt((t(v_x) * v_y))))})), "(1/2)*y*(x.T*y)^(-1/2)")
  expect_equal(.bat_fmt(local({v_x <- MatrixSymbol("x", Symbol("k"), Integer(1)); v_expr <- MatPow(MatMul(Transpose(MatrixSymbol("x", Symbol("k"), Integer(1))), MatrixSymbol("x", Symbol("k"), Integer(1))), Rational(1, 2)); diff(v_expr, v_x)})), "x*(x.T*x)^(-1/2)")
  expect_equal(.bat_fmt(local({v_x <- MatrixSymbol("x", Symbol("k"), Integer(1)); (v_x * ((t(v_x) * v_x) ^ Rational((-1), 2)))})), "x*(x.T*x)^(-1/2)")
  expect_equal(.bat_fmt(local({v_x <- MatrixSymbol("x", Symbol("k"), Integer(1)); v_expr <- MatPow(MatMul(Transpose(MatrixSymbol("c", Symbol("k"), Integer(1))), MatrixSymbol("a", Symbol("k"), Integer(1)), Transpose(MatrixSymbol("x", Symbol("k"), Integer(1))), MatrixSymbol("b", Symbol("k"), Integer(1))), Rational(1, 2)); diff(v_expr, v_x)})), "(1/2)*b*a.T*c*(c.T*a*x.T*b)^(-1/2)")
  expect_equal(.bat_fmt(local({v_b <- MatrixSymbol("b", Symbol("k"), Integer(1)); v_a <- MatrixSymbol("a", Symbol("k"), Integer(1)); v_c <- MatrixSymbol("c", Symbol("k"), Integer(1)); v_x <- MatrixSymbol("x", Symbol("k"), Integer(1)); ((((v_b * t(v_a)) * v_c) / sqrt((((t(v_c) * v_a) * t(v_x)) * v_b))) / 2)})), "(1/2)*b*a.T*c*(c.T*a*x.T*b)^(-1/2)")
  expect_equal(.bat_fmt(local({v_x <- MatrixSymbol("x", Symbol("k"), Integer(1)); v_expr <- MatPow(MatMul(Transpose(MatrixSymbol("c", Symbol("k"), Integer(1))), MatrixSymbol("a", Symbol("k"), Integer(1)), Transpose(MatrixSymbol("x", Symbol("k"), Integer(1))), MatrixSymbol("b", Symbol("k"), Integer(1))), Rational(1, 3)); diff(v_expr, v_x)})), "(1/3)*b*a.T*c*(c.T*a*x.T*b)^(-2/3)")
  expect_equal(.bat_fmt(local({v_b <- MatrixSymbol("b", Symbol("k"), Integer(1)); v_a <- MatrixSymbol("a", Symbol("k"), Integer(1)); v_c <- MatrixSymbol("c", Symbol("k"), Integer(1)); v_x <- MatrixSymbol("x", Symbol("k"), Integer(1)); ((((v_b * t(v_a)) * v_c) * ((((t(v_c) * v_a) * t(v_x)) * v_b) ^ Rational((-2), 3))) / 3)})), "(1/3)*b*a.T*c*(c.T*a*x.T*b)^(-2/3)")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_expr <- MatPow(MatMul(Transpose(MatrixSymbol("a", Symbol("k"), Integer(1))), MatrixSymbol("X", Symbol("k"), Symbol("k")), MatrixSymbol("b", Symbol("k"), Integer(1))), Rational(1, 2)); diff(v_expr, v_X)})), "(1/2)*a*(a.T*X*b)^(-1/2)*b.T")
  expect_equal(.bat_fmt(local({v_a <- MatrixSymbol("a", Symbol("k"), Integer(1)); v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_b <- MatrixSymbol("b", Symbol("k"), Integer(1)); ((v_a / (2 * sqrt(((t(v_a) * v_X) * v_b)))) * t(v_b))})), "(1/2)*a*(a.T*X*b)^(-1/2)*b.T")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_expr <- MatMul(Transpose(MatrixSymbol("d", Symbol("k"), Integer(1))), MatrixSymbol("x", Symbol("k"), Integer(1)), MatPow(MatMul(Transpose(MatrixSymbol("a", Symbol("k"), Integer(1))), MatrixSymbol("X", Symbol("k"), Symbol("k")), MatrixSymbol("b", Symbol("k"), Integer(1))), Rational(1, 2)), Transpose(MatrixSymbol("y", Symbol("k"), Integer(1))), MatrixSymbol("c", Symbol("k"), Integer(1))); diff(v_expr, v_X)})), "(1/2)*a*x.T*d*y.T*c*(b.T*X.T*a)^(-1/2)*b.T")
  expect_equal(.bat_fmt(local({v_a <- MatrixSymbol("a", Symbol("k"), Integer(1)); v_x <- MatrixSymbol("x", Symbol("k"), Integer(1)); v_d <- MatrixSymbol("d", Symbol("k"), Integer(1)); v_y <- MatrixSymbol("y", Symbol("k"), Integer(1)); v_c <- MatrixSymbol("c", Symbol("k"), Integer(1)); v_b <- MatrixSymbol("b", Symbol("k"), Integer(1)); v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); ((((((v_a * t(v_x)) * v_d) * t(v_y)) * v_c) / (2 * sqrt(((t(v_b) * t(v_X)) * v_a)))) * t(v_b))})), "(1/2)*a*x.T*d*y.T*c*(b.T*X.T*a)^(-1/2)*b.T")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_expr <- (Trace(MatMul(Transpose(MatrixSymbol("X", Symbol("k"), Symbol("k"))), MatrixSymbol("X", Symbol("k"), Symbol("k")))))^(Rational(1, 2)); diff(v_expr, v_X)})), "1/sqrt(Trace(X*X.T))*X")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); (v_X / sqrt(Trace((v_X * t(v_X)))))})), "1/sqrt(Trace(X*X.T))*X")
})

test_that("test_derivatives: test_derivatives_elementwise_applyfunc matches SymPy", {
  expect_equal(.bat_fmt(local({.d1 <- Dummy("d"); v_x <- MatrixSymbol("x", Symbol("k"), Integer(1)); v_expr <- ElementwiseApplyFunction(Lambda(.d1, tan(.d1)), MatrixSymbol("x", Symbol("k"), Integer(1))); diff(v_expr, v_x)})), "DiagMatrix(Lambda(_d, tan(_d)^2 + 1).(x))")
  expect_equal(.bat_fmt(local({v_x <- MatrixSymbol("x", Symbol("k"), Integer(1)); DiagMatrix(applyfunc(v_x, function(x) ((tan(x) ^ 2) + 1)))})), "DiagMatrix(Lambda(_d, tan(_d)^2 + 1).(x))")
  expect_equal(.bat_fmt(local({.d1 <- Dummy("d"); v_x <- MatrixSymbol("x", Symbol("k"), Integer(1)); v_m <- Symbol("m"); v_expr <- ElementwiseApplyFunction(Lambda(.d1, tan(.d1)), MatrixSymbol("x", Symbol("k"), Integer(1))); v_i <- Symbol("i"); doit(diff(pyitem(v_expr, v_i, 0), pyitem(v_x, v_m, 0)))})), "(tan(x[i, 0])^2 + 1)*KroneckerDelta(i, m, (0, k - 1))")
  expect_equal(.bat_fmt(local({v_x <- MatrixSymbol("x", Symbol("k"), Integer(1)); v_i <- Symbol("i"); v_m <- Symbol("m"); v_k <- Symbol("k"); (((tan(pyitem(v_x, v_i, 0)) ^ 2) + 1) * KroneckerDelta(v_i, v_m, list(0, (v_k - 1))))})), "(tan(x[i, 0])^2 + 1)*KroneckerDelta(i, m, (0, k - 1))")
  expect_equal(.bat_fmt(local({.d1 <- Dummy("d"); v_i <- Symbol("i"); v_expr <- ElementwiseApplyFunction(Lambda(.d1, sin(.d1)), MatMul((Symbol("i"))^(Integer(2)), MatrixSymbol("x", Symbol("k"), Integer(1)))); diff(v_expr, v_i)})), "HadamardProduct((2*i)*x, Lambda(_d, cos(_d)).(i^2*x))")
  expect_equal(.bat_fmt(local({v_i <- Symbol("i"); v_x <- MatrixSymbol("x", Symbol("k"), Integer(1)); HadamardProduct(((2 * v_i) * v_x), applyfunc(((v_i ^ 2) * v_x), cos))})), "HadamardProduct((2*i)*x, Lambda(_d, cos(_d)).(i^2*x))")
  expect_equal(.bat_fmt(local({.d1 <- Dummy("d"); v_i <- Symbol("i"); v_expr <- ElementwiseApplyFunction(Lambda(.d1, sin(.d1)), MatMul((Symbol("i"))^(Integer(2)), MatrixSymbol("x", Symbol("k"), Integer(1)))); doit(diff(pyitem(v_expr, v_i, 0), v_i))})), "2*i*cos(i^2*x[i, 0])*x[i, 0]")
  expect_equal(.bat_fmt(local({v_i <- Symbol("i"); v_x <- MatrixSymbol("x", Symbol("k"), Integer(1)); (((2 * v_i) * pyitem(v_x, v_i, 0)) * cos(((v_i ^ 2) * pyitem(v_x, v_i, 0))))})), "2*i*cos(i^2*x[i, 0])*x[i, 0]")
  expect_equal(.bat_fmt(local({.d1 <- Dummy("d"); v_i <- Symbol("i"); v_expr <- ElementwiseApplyFunction(Lambda(.d1, sin(.d1)), MatMul(log(Symbol("i")), MatrixSymbol("A", Symbol("k"), Symbol("k")), MatrixSymbol("B", Symbol("k"), Symbol("k")))); diff(v_expr, v_i)})), "HadamardProduct(1/i*A*B, Lambda(_d, cos(_d)).(log(i)*A*B))")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("k"), Symbol("k")); v_B <- MatrixSymbol("B", Symbol("k"), Symbol("k")); v_i <- Symbol("i"); HadamardProduct(((v_A * v_B) / v_i), applyfunc(((log(v_i) * v_A) * v_B), cos))})), "HadamardProduct(1/i*A*B, Lambda(_d, cos(_d)).(log(i)*A*B))")
  expect_equal(.bat_fmt(local({.d1 <- Dummy("d"); v_x <- MatrixSymbol("x", Symbol("k"), Integer(1)); v_expr <- MatAdd(MatMul(Symbol("k"), ElementwiseApplyFunction(Lambda(.d1, sin(.d1)), Transpose(MatrixSymbol("y", Symbol("k"), Integer(1)))), MatrixSymbol("x", Symbol("k"), Integer(1))), MatMul(Transpose(MatrixSymbol("x", Symbol("k"), Integer(1))), MatrixSymbol("A", Symbol("k"), Symbol("k")), MatrixSymbol("x", Symbol("k"), Integer(1)))); diff(v_expr, v_x)})), "k*Lambda(_d, sin(_d)).(y) + A*x + A.T*x")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("k"), Symbol("k")); v_x <- MatrixSymbol("x", Symbol("k"), Integer(1)); v_k <- Symbol("k"); v_y <- MatrixSymbol("y", Symbol("k"), Integer(1)); (((t(v_A) * v_x) + (v_A * v_x)) + (v_k * applyfunc(v_y, sin)))})), "k*Lambda(_d, sin(_d)).(y) + A*x + A.T*x")
  expect_equal(.bat_fmt(local({.d1 <- Dummy("d"); v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_expr <- ElementwiseApplyFunction(Lambda(.d1, sin(.d1)), MatMul(Transpose(MatrixSymbol("a", Symbol("k"), Integer(1))), MatrixSymbol("X", Symbol("k"), Symbol("k")), MatrixSymbol("b", Symbol("k"), Integer(1)))); diff(v_expr, v_X)})), "a*Lambda(_d, cos(_d)).(a.T*X*b)*b.T")
  expect_equal(.bat_fmt(local({v_a <- MatrixSymbol("a", Symbol("k"), Integer(1)); v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_b <- MatrixSymbol("b", Symbol("k"), Integer(1)); ((v_a * applyfunc(((t(v_a) * v_X) * v_b), cos)) * t(v_b))})), "a*Lambda(_d, cos(_d)).(a.T*X*b)*b.T")
  expect_equal(.bat_fmt(local({.d1 <- Dummy("d"); v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_expr <- MatMul(Transpose(MatrixSymbol("a", Symbol("k"), Integer(1))), ElementwiseApplyFunction(Lambda(.d1, sin(.d1)), MatrixSymbol("X", Symbol("k"), Symbol("k"))), MatrixSymbol("b", Symbol("k"), Integer(1))); diff(v_expr, v_X)})), "DiagMatrix(a)*Lambda(_d, cos(_d)).(X)*DiagMatrix(b)")
  expect_equal(.bat_fmt(local({v_a <- MatrixSymbol("a", Symbol("k"), Integer(1)); v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_b <- MatrixSymbol("b", Symbol("k"), Integer(1)); ((DiagMatrix(v_a) * applyfunc(v_X, cos)) * DiagMatrix(v_b))})), "DiagMatrix(a)*Lambda(_d, cos(_d)).(X)*DiagMatrix(b)")
  expect_equal(.bat_fmt(local({.d1 <- Dummy("d"); v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_expr <- MatMul(Transpose(MatrixSymbol("a", Symbol("k"), Integer(1))), ElementwiseApplyFunction(Lambda(.d1, sin(.d1)), MatMul(MatrixSymbol("A", Symbol("k"), Symbol("k")), MatrixSymbol("X", Symbol("k"), Symbol("k")), MatrixSymbol("B", Symbol("k"), Symbol("k")))), MatrixSymbol("b", Symbol("k"), Integer(1))); diff(v_expr, v_X)})), "A.T*DiagMatrix(a)*Lambda(_d, cos(_d)).(A*X*B)*DiagMatrix(b)*B.T")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("k"), Symbol("k")); v_a <- MatrixSymbol("a", Symbol("k"), Integer(1)); v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_B <- MatrixSymbol("B", Symbol("k"), Symbol("k")); v_b <- MatrixSymbol("b", Symbol("k"), Integer(1)); ((((t(v_A) * DiagMatrix(v_a)) * applyfunc(((v_A * v_X) * v_B), cos)) * DiagMatrix(v_b)) * t(v_B))})), "A.T*DiagMatrix(a)*Lambda(_d, cos(_d)).(A*X*B)*DiagMatrix(b)*B.T")
  expect_equal(.bat_fmt(local({.d1 <- Dummy("d"); v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_expr <- MatMul(Transpose(MatrixSymbol("a", Symbol("k"), Integer(1))), MatrixSymbol("A", Symbol("k"), Symbol("k")), ElementwiseApplyFunction(Lambda(.d1, sin(.d1)), MatrixSymbol("X", Symbol("k"), Symbol("k"))), MatrixSymbol("B", Symbol("k"), Symbol("k")), MatrixSymbol("b", Symbol("k"), Integer(1))); diff(v_expr, v_X)})), "HadamardProduct(A.T*a*b.T*B.T, Lambda(_d, cos(_d)).(X))")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("k"), Symbol("k")); v_a <- MatrixSymbol("a", Symbol("k"), Integer(1)); v_b <- MatrixSymbol("b", Symbol("k"), Integer(1)); v_B <- MatrixSymbol("B", Symbol("k"), Symbol("k")); v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); HadamardProduct((((t(v_A) * v_a) * t(v_b)) * t(v_B)), applyfunc(v_X, cos))})), "HadamardProduct(A.T*a*b.T*B.T, Lambda(_d, cos(_d)).(X))")
})

test_that("test_derivatives: test_derivatives_of_hadamard_expressions matches SymPy", {
  expect_equal(.bat_fmt(local({v_x <- MatrixSymbol("x", Symbol("k"), Integer(1)); v_expr <- HadamardProduct(MatrixSymbol("a", Symbol("k"), Integer(1)), MatrixSymbol("b", Symbol("k"), Integer(1)), MatrixSymbol("x", Symbol("k"), Integer(1))); diff(v_expr, v_x)})), "DiagMatrix(HadamardProduct(a, b))")
  expect_equal(.bat_fmt(local({v_b <- MatrixSymbol("b", Symbol("k"), Integer(1)); v_a <- MatrixSymbol("a", Symbol("k"), Integer(1)); DiagMatrix(hadamard_product(v_b, v_a))})), "DiagMatrix(HadamardProduct(a, b))")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_expr <- MatMul(Transpose(MatrixSymbol("a", Symbol("k"), Integer(1))), HadamardProduct(MatrixSymbol("A", Symbol("k"), Symbol("k")), MatrixSymbol("B", Symbol("k"), Symbol("k")), MatrixSymbol("X", Symbol("k"), Symbol("k"))), MatrixSymbol("b", Symbol("k"), Integer(1))); diff(v_expr, v_X)})), "HadamardProduct(a*b.T, A, B)")
  expect_equal(.bat_fmt(local({v_a <- MatrixSymbol("a", Symbol("k"), Integer(1)); v_b <- MatrixSymbol("b", Symbol("k"), Integer(1)); v_A <- MatrixSymbol("A", Symbol("k"), Symbol("k")); v_B <- MatrixSymbol("B", Symbol("k"), Symbol("k")); HadamardProduct((v_a * t(v_b)), v_A, v_B)})), "HadamardProduct(a*b.T, A, B)")
  expect_equal(.bat_fmt(local({v_x <- MatrixSymbol("x", Symbol("k"), Integer(1)); v_expr <- HadamardPower(MatrixSymbol("x", Symbol("k"), Integer(1)), Integer(2)); doit(diff(v_expr, v_x))})), "2*DiagMatrix(x)")
  expect_equal(.bat_fmt(local({v_x <- MatrixSymbol("x", Symbol("k"), Integer(1)); (2 * DiagMatrix(v_x))})), "2*DiagMatrix(x)")
  expect_equal(.bat_fmt(local({v_x <- MatrixSymbol("x", Symbol("k"), Integer(1)); v_expr <- HadamardPower(Transpose(MatrixSymbol("x", Symbol("k"), Integer(1))), Integer(2)); doit(diff(v_expr, v_x))})), "2*DiagMatrix(x)")
  expect_equal(.bat_fmt(local({v_x <- MatrixSymbol("x", Symbol("k"), Integer(1)); v_expr <- HadamardPower(MatrixSymbol("x", Symbol("k"), Integer(1)), Rational(1, 2)); diff(v_expr, v_x)})), "(1/2)*DiagMatrix(HadamardPower(x, -1/2))")
  expect_equal(.bat_fmt(local({v_x <- MatrixSymbol("x", Symbol("k"), Integer(1)); (S_Half * DiagMatrix(hadamard_power(v_x, Rational((-1), 2))))})), "(1/2)*DiagMatrix(HadamardPower(x, -1/2))")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_expr <- HadamardPower(MatMul(Transpose(MatrixSymbol("a", Symbol("k"), Integer(1))), MatrixSymbol("X", Symbol("k"), Symbol("k")), MatrixSymbol("b", Symbol("k"), Integer(1))), Integer(2)); diff(v_expr, v_X)})), "2*a*a.T*X*b*b.T")
  expect_equal(.bat_fmt(local({v_a <- MatrixSymbol("a", Symbol("k"), Integer(1)); v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_b <- MatrixSymbol("b", Symbol("k"), Integer(1)); (((((2 * v_a) * t(v_a)) * v_X) * v_b) * t(v_b))})), "2*a*a.T*X*b*b.T")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("k"), Symbol("k")); v_expr <- HadamardPower(MatMul(Transpose(MatrixSymbol("a", Symbol("k"), Integer(1))), MatrixSymbol("X", Symbol("k"), Symbol("k")), MatrixSymbol("b", Symbol("k"), Integer(1))), Rational(1, 2)); diff(v_expr, v_X)})), "(1/2)*a*(a.T*X*b)^(-1/2)*b.T")
})

test_that("test_derivatives: test_fix_issue_24021 matches SymPy", {
  expect_equal(.bat_fmt(local({v_E <- MatrixSymbol("E", Integer(3), Integer(3)); v_expr <- ((Trace(MatrixSymbol("E", Integer(3), Integer(3))))^(Integer(2)) + Trace(MatMul(MatrixSymbol("E", Integer(3), Integer(3)), Transpose(MatrixSymbol("E", Integer(3), Integer(3)))))); diff(v_expr, v_E)})), "(2*Trace(E))*I + 2*E")
  expect_equal(.bat_fmt(local({v_E <- MatrixSymbol("E", Integer(3), Integer(3)); (((2 * Trace(v_E)) * Identity(3)) + (2 * v_E))})), "(2*Trace(E))*I + 2*E")
})

test_that("test_derivatives: test_fix_issue_on_stackoverflow matches SymPy", {
  expect_equal(.bat_fmt(local({v_x <- MatrixSymbol("x", Integer(2), Integer(1)); v_fx <- MatAdd(MatMul(Integer(6), Identity(Integer(1))), MatMul(Transpose(MatrixSymbol("x", Integer(2), Integer(1))), Matrix(list(list(Integer(3)), list(Integer(5))))), MatMul(Transpose(MatrixSymbol("x", Integer(2), Integer(1))), Matrix(list(list(Integer(1), Integer(2)), list(Integer(4), Integer(7)))), MatrixSymbol("x", Integer(2), Integer(1)))); diff(v_fx, v_x)})), "Matrix([\\n[3],\\n[5]]) + Matrix([\\n[1, 2],\\n[4, 7]])*x + Matrix([\\n[1, 4],\\n[2, 7]])*x")
  expect_equal(.bat_fmt(local({v_m1 <- Matrix(list(list(Integer(1), Integer(2)), list(Integer(4), Integer(7)))); v_x <- MatrixSymbol("x", Integer(2), Integer(1)); v_m3 <- Matrix(list(list(Integer(3)), list(Integer(5)))); (((v_m1 * v_x) + (t(v_m1) * v_x)) + v_m3)})), "Matrix([\\n[3],\\n[5]]) + Matrix([\\n[1, 2],\\n[4, 7]])*x + Matrix([\\n[1, 4],\\n[2, 7]])*x")
})

test_that("test_derivatives: test_issue_25118 matches SymPy", {
  expect_equal(.bat_fmt(local({v_F <- MatrixSymbol("F", Integer(3), Integer(3)); diff((sym_trace(v_F) - 1), v_F)})), "I")
  expect_equal(.bat_fmt(local({Identity(3)})), "I")
  expect_equal(.bat_fmt(local({v_P <- MatMul(Integer(2), MatrixSymbol("F", Integer(3), Integer(3))); v_P})), "2*F")
  expect_equal(.bat_fmt(local({v_F <- MatrixSymbol("F", Integer(3), Integer(3)); (2 * v_F)})), "2*F")
})

test_that("test_derivatives: test_issue_24725 matches SymPy", {
  expect_equal(.bat_fmt(local({v_d <- MatAdd(MatMul((Integer(-1)*Symbol("m")), OneMatrix(Symbol("N"), Integer(1))), MatrixSymbol("x", Symbol("N"), Integer(1))); v_m <- Symbol("m", commutative = TRUE); diff(Trace((t(v_d) * v_d)), v_m)})), "-Trace(((-m)*1 + x.T)*1) - Trace(1*((-m)*1 + x))")
  expect_equal(.bat_fmt(local({v_m <- Symbol("m", commutative = TRUE); v_one <- OneMatrix(Symbol("N"), Integer(1)); v_x <- MatrixSymbol("x", Symbol("N"), Integer(1)); ((-Trace(((((-v_m) * t(v_one)) + t(v_x)) * v_one))) - Trace((t(v_one) * (((-v_m) * v_one) + v_x))))})), "-Trace(((-m)*1 + x.T)*1) - Trace(1*((-m)*1 + x))")
  expect_equal(.bat_fmt(local({v_x <- MatrixSymbol("x", Symbol("N"), Integer(1)); v_w <- MatrixSymbol("w", Symbol("N"), Integer(1)); v_b <- Symbol("b"); diff((Trace((t(v_x) * v_w)) + v_b), v_w)})), "x")
  expect_equal(.bat_fmt(local({v_x <- MatrixSymbol("x", Symbol("N"), Integer(1)); v_x})), "x")
  expect_equal(.bat_fmt(local({v_w <- MatrixSymbol("w", Symbol("N"), Integer(1)); diff(pyitem(v_w, 0, 0), v_w)})), "MatrixUnit(N, 1, 0, 0)")
  expect_equal(.bat_fmt(local({v_w <- MatrixSymbol("w", Symbol("N"), Integer(1)); MatrixUnit(pyitem(get_shape(v_w), 0), pyitem(get_shape(v_w), 1), 0, 0)})), "MatrixUnit(N, 1, 0, 0)")
})

test_that("test_derivatives: test_issue_23492 matches SymPy", {
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Integer(3), Integer(3)); diff(sym_trace(v_A), v_A)})), "I")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Integer(3), Integer(3)); diff(as_explicit(sym_trace(v_A)), v_A)})), "MatrixUnit(3, 3, 0, 0) + MatrixUnit(3, 3, 1, 1) + MatrixUnit(3, 3, 2, 2)")
  expect_equal(.bat_fmt(local({((MatrixUnit(3, 3, 0, 0) + MatrixUnit(3, 3, 1, 1)) + MatrixUnit(3, 3, 2, 2))})), "MatrixUnit(3, 3, 0, 0) + MatrixUnit(3, 3, 1, 1) + MatrixUnit(3, 3, 2, 2)")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Integer(3), Integer(3)); diff(sym_trace(v_A), as_explicit(v_A))})), "Matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])")
  expect_equal(.bat_fmt(local({eye(3)})), "Matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Integer(3), Integer(3)); diff((sym_trace(v_A) + sym_trace((v_A ^ 2))), v_A)})), "I + 2*A.T")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Integer(3), Integer(3)); (Identity(3) + (2 * t(v_A)))})), "I + 2*A.T")
})

test_that("test_derivatives: test_issue_23364 matches SymPy", {
  expect_equal(.bat_fmt(local({v_x <- Symbol("x"); v_F <- MatrixSymbol("F", Integer(3), Integer(3)); diff(((v_x * v_F) ^ 2), v_x)})), "(2*x)*F^2")
  expect_equal(.bat_fmt(local({v_x <- Symbol("x"); v_F <- MatrixSymbol("F", Integer(3), Integer(3)); ((2 * v_x) * (v_F ^ 2))})), "(2*x)*F^2")
})

test_that("test_derivatives: test_issue_21597 matches SymPy", {
  expect_equal(.bat_fmt(local({v_b <- MatrixSymbol("b", Integer(3), Integer(1)); v_F <- MatrixSymbol("F", Integer(3), Integer(3)); v_c <- MatrixSymbol("c", Integer(3), Integer(1)); diff((v_b - (v_F * v_c)), v_b)})), "I")
})

test_that("test_derivatives: test_issue_23931 matches SymPy", {
  expect_equal(.bat_fmt(local({v_F <- MatrixSymbol("F", Integer(3), Integer(3)); v_expr <- (Trace(MatrixSymbol("F", Integer(3), Integer(3))))^(Integer(-1)); diff(v_expr, v_F)})), "(-1/Trace(F)^2)*I")
  expect_equal(.bat_fmt(local({v_result <- MatMul((Integer(-1)*(Trace(MatrixSymbol("F", Integer(3), Integer(3))))^(Integer(-2))), Identity(Integer(3))); v_result})), "(-1/Trace(F)^2)*I")
  expect_equal(.bat_fmt(local({v_x <- MatrixSymbol("x", Symbol("d", integer = TRUE, positive = TRUE), Integer(1)); v_expr <- MatMul((Trace(MatMul(Transpose(MatrixSymbol("x", Symbol("d", integer = TRUE, positive = TRUE), Integer(1))), MatrixSymbol("x", Symbol("d", integer = TRUE, positive = TRUE), Integer(1)))))^(Integer(-1)), MatrixSymbol("x", Symbol("d", integer = TRUE, positive = TRUE), Integer(1))); diff(v_expr, v_x)})), "(-2/Trace(x*x.T)^2)*x*x.T + 1/Trace(x*x.T)*I")
  expect_equal(.bat_fmt(local({v_result <- MatAdd(MatMul((Trace(MatMul(MatrixSymbol("x", Symbol("d", integer = TRUE, positive = TRUE), Integer(1)), Transpose(MatrixSymbol("x", Symbol("d", integer = TRUE, positive = TRUE), Integer(1))))))^(Integer(-1)), Identity(Symbol("d", integer = TRUE, positive = TRUE))), MatMul((Integer(-2)*(Trace(MatMul(MatrixSymbol("x", Symbol("d", integer = TRUE, positive = TRUE), Integer(1)), Transpose(MatrixSymbol("x", Symbol("d", integer = TRUE, positive = TRUE), Integer(1))))))^(Integer(-2))), MatrixSymbol("x", Symbol("d", integer = TRUE, positive = TRUE), Integer(1)), Transpose(MatrixSymbol("x", Symbol("d", integer = TRUE, positive = TRUE), Integer(1))))); v_result})), "(-2/Trace(x*x.T)^2)*x*x.T + 1/Trace(x*x.T)*I")
})

test_that("test_derivatives: test_issue_17229 matches SymPy", {
  expect_equal(.bat_fmt(local({v_M <- MatrixSymbol("M", Integer(3), Integer(3)); diff(pyitem(v_M, 1, 2), v_M)})), "MatrixUnit(3, 3, 1, 2)")
  expect_equal(.bat_fmt(local({MatrixUnit(3, 3, 1, 2)})), "MatrixUnit(3, 3, 1, 2)")
  expect_equal(.bat_fmt(local({v_M <- MatrixSymbol("M", Integer(3), Integer(3)); diff(pyitem(v_M, 1, 2), as_explicit(v_M))})), "Matrix([[0, 0, 0], [0, 0, 1], [0, 0, 0]])")
  expect_equal(.bat_fmt(local({Matrix(list(list(0, 0, 0), list(0, 0, 1), list(0, 0, 0)))})), "Matrix([[0, 0, 0], [0, 0, 1], [0, 0, 0]])")
  expect_equal(.bat_fmt(local({v_M <- MatrixSymbol("M", Integer(3), Integer(3)); diff(v_M, pyitem(v_M, 1, 1))})), "MatrixUnit(3, 3, 1, 1)")
  expect_equal(.bat_fmt(local({MatrixUnit(3, 3, 1, 1)})), "MatrixUnit(3, 3, 1, 1)")
  expect_equal(.bat_fmt(local({v_M <- MatrixSymbol("M", Integer(3), Integer(3)); diff(as_explicit(v_M), pyitem(v_M, 1, 2))})), "Matrix([[0, 0, 0], [0, 0, 1], [0, 0, 0]])")
  expect_equal(.bat_fmt(local({v_M <- MatrixSymbol("M", Integer(3), Integer(3)); diff(v_M, as_explicit(v_M))})), "[[MatrixUnit(3, 3, 0, 0), MatrixUnit(3, 3, 0, 1), MatrixUnit(3, 3, 0, 2)], [MatrixUnit(3, 3, 1, 0), MatrixUnit(3, 3, 1, 1), MatrixUnit(3, 3, 1, 2)], [MatrixUnit(3, 3, 2, 0), MatrixUnit(3, 3, 2, 1), MatrixUnit(3, 3, 2, 2)]]")
  expect_equal(.bat_fmt(local({v_M <- MatrixSymbol("M", Integer(3), Integer(3)); diff(as_explicit(v_M), v_M)})), "Matrix([[MatrixUnit(3, 3, 0, 0), MatrixUnit(3, 3, 0, 1), MatrixUnit(3, 3, 0, 2)], [MatrixUnit(3, 3, 1, 0), MatrixUnit(3, 3, 1, 1), MatrixUnit(3, 3, 1, 2)], [MatrixUnit(3, 3, 2, 0), MatrixUnit(3, 3, 2, 1), MatrixUnit(3, 3, 2, 2)]])")
})

test_that("test_derivatives: test_issue_15651 matches SymPy", {
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Integer(3), Integer(3)); v_expr <- MatMul(Matrix(list(list(MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(0), Integer(0)), MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(0), Integer(1)), MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(0), Integer(2))), list(MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(1), Integer(0)), MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(1), Integer(1)), MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(1), Integer(2))), list(MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(2), Integer(0)), MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(2), Integer(1)), MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(2), Integer(2))))), MatrixSymbol("X", Integer(3), Integer(3))); diff(v_expr, v_X)})), "ArrayAdd(ArrayContraction(ArrayTensorProduct([[[[1, 0, 0], [0, 0, 0], [0, 0, 0]], [[0, 1, 0], [0, 0, 0], [0, 0, 0]], [[0, 0, 1], [0, 0, 0], [0, 0, 0]]], [[[0, 0, 0], [1, 0, 0], [0, 0, 0]], [[0, 0, 0], [0, 1, 0], [0, 0, 0]], [[0, 0, 0], [0, 0, 1], [0, 0, 0]]], [[[0, 0, 0], [0, 0, 0], [1, 0, 0]], [[0, 0, 0], [0, 0, 0], [0, 1, 0]], [[0, 0, 0], [0, 0, 0], [0, 0, 1]]]], X), (3, 4)), PermuteDims(ArrayTensorProduct(Matrix([\\n[X[0, 0], X[1, 0], X[2, 0]],\\n[X[0, 1], X[1, 1], X[2, 1]],\\n[X[0, 2], X[1, 2], X[2, 2]]]), I), (3)(1 2)))")
  expect_equal(.bat_fmt(local({v_expected <- ArrayAdd(ArrayContraction(ArrayTensorProduct(Array(list(Integer(1), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(1), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(1), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(1), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(1), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(1), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(1), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(1), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(1)), c(3, 3, 3, 3)), MatrixSymbol("X", Integer(3), Integer(3))), Tuple(Integer(3), Integer(4))), PermuteDims(ArrayTensorProduct(Matrix(list(list(MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(0), Integer(0)), MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(1), Integer(0)), MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(2), Integer(0))), list(MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(0), Integer(1)), MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(1), Integer(1)), MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(2), Integer(1))), list(MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(0), Integer(2)), MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(1), Integer(2)), MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(2), Integer(2))))), Identity(Integer(3))), Permutation(c(0, 2, 1, 3)))); v_expected})), "ArrayAdd(ArrayContraction(ArrayTensorProduct([[[[1, 0, 0], [0, 0, 0], [0, 0, 0]], [[0, 1, 0], [0, 0, 0], [0, 0, 0]], [[0, 0, 1], [0, 0, 0], [0, 0, 0]]], [[[0, 0, 0], [1, 0, 0], [0, 0, 0]], [[0, 0, 0], [0, 1, 0], [0, 0, 0]], [[0, 0, 0], [0, 0, 1], [0, 0, 0]]], [[[0, 0, 0], [0, 0, 0], [1, 0, 0]], [[0, 0, 0], [0, 0, 0], [0, 1, 0]], [[0, 0, 0], [0, 0, 0], [0, 0, 1]]]], X), (3, 4)), PermuteDims(ArrayTensorProduct(Matrix([\\n[X[0, 0], X[1, 0], X[2, 0]],\\n[X[0, 1], X[1, 1], X[2, 1]],\\n[X[0, 2], X[1, 2], X[2, 2]]]), I), (3)(1 2)))")
})

test_that("test_derivatives: test_matpow_derivative matches SymPy", {
  expect_equal(.bat_fmt(local({v_m <- Symbol("m"); v_M <- MatPow(Matrix(list(list(Symbol("m"), Integer(0)), list(Integer(0), Symbol("m")))), Integer(2)); diff(v_M, v_m)})), "Matrix([[2*m, 0], [0, 2*m]])")
  expect_equal(.bat_fmt(local({v_m <- Symbol("m"); Matrix(list(list((2 * v_m), 0), list(0, (2 * v_m))))})), "Matrix([[2*m, 0], [0, 2*m]])")
  expect_equal(.bat_fmt(local({v_m <- Symbol("m"); v_M <- MatPow(Matrix(list(list(Integer(0), Symbol("m")), list(Symbol("n"), Integer(0)))), Integer(3)); diff(v_M, v_m)})), "Matrix([[0, 2*m*n], [n^2, 0]])")
  expect_equal(.bat_fmt(local({v_m <- Symbol("m"); v_n <- Symbol("n"); Matrix(list(list(0, ((2 * v_m) * v_n)), list((v_n ^ 2), 0)))})), "Matrix([[0, 2*m*n], [n^2, 0]])")
  expect_equal(.bat_fmt(local({ZeroMatrix(2, 2)})), "0")
})

test_that("test_derivatives: test_matexpr_derivative_by_matrix_element matches SymPy", {
  expect_equal(.bat_fmt(local({v_Y <- MatrixSymbol("Y", Integer(2), Integer(2)); v_X <- MatrixSymbol("X", Integer(2), Integer(2)); doit(diff((v_X + v_Y), pyitem(v_Y, 0, 0)))})), "MatrixUnit(2, 2, 0, 0)")
  expect_equal(.bat_fmt(local({MatrixUnit(2, 2, 0, 0)})), "MatrixUnit(2, 2, 0, 0)")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Integer(2), Integer(2)); v_a <- Symbol("a"); v_Y <- MatrixSymbol("Y", Integer(2), Integer(2)); diff((v_X + (v_a * v_Y)), pyitem(v_X, v_a, 0))})), "MatrixUnit(2, 2, a, 0)")
  expect_equal(.bat_fmt(local({v_a <- Symbol("a"); MatrixUnit(2, 2, v_a, 0)})), "MatrixUnit(2, 2, a, 0)")
  expect_equal(.bat_fmt(local({v_Y <- MatrixSymbol("Y", Integer(2), Integer(2)); v_a <- Symbol("a"); v_X <- MatrixSymbol("X", Integer(2), Integer(2)); diff((v_X + (v_a * v_Y)), pyitem(v_Y, v_a, 0))})), "a*MatrixUnit(2, 2, a, 0)")
  expect_equal(.bat_fmt(local({v_a <- Symbol("a"); (v_a * MatrixUnit(2, 2, v_a, 0))})), "a*MatrixUnit(2, 2, a, 0)")
  expect_equal(.bat_fmt(local({v_Z <- MatrixSymbol("Z", Integer(2), Integer(2)); v_a <- Symbol("a"); v_X <- MatrixSymbol("X", Integer(2), Integer(2)); v_Y <- MatrixSymbol("Y", Integer(2), Integer(2)); diff((v_X + (v_a * v_Y)), pyitem(v_Z, v_a, 0))})), "0")
  expect_equal(.bat_fmt(local({v_a <- Symbol("a"); v_X <- MatrixSymbol("X", Integer(2), Integer(2)); v_Y <- MatrixSymbol("Y", Integer(2), Integer(2)); diff((v_X + (v_a * v_Y)), v_a)})), "Y")
  expect_equal(.bat_fmt(local({v_Y <- MatrixSymbol("Y", Integer(2), Integer(2)); (v_Y + ZeroMatrix(2, 2))})), "Y")
  expect_equal(.bat_fmt(local({v_Y <- MatrixSymbol("Y", Integer(2), Integer(2)); v_a <- Symbol("a"); v_A <- MatrixSymbol("a", Symbol("a"), Symbol("a")); diff(v_A, pyitem(v_Y, v_a, v_a))})), "0")
  expect_equal(.bat_fmt(local({v_a <- Symbol("a"); ZeroMatrix(v_a, v_a)})), "0")
})

test_that("test_derivatives: test_matrix_element_derivative_by_matrix_symbol matches SymPy", {
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Integer(3), Integer(3)); diff(pyitem(v_X, 0, 1), v_X)})), "MatrixUnit(3, 3, 0, 1)")
  expect_equal(.bat_fmt(local({MatrixUnit(3, 3, 0, 1)})), "MatrixUnit(3, 3, 0, 1)")
  expect_equal(.bat_fmt(local({v_t <- Symbol("t"); v_X <- MatrixSymbol("X", Integer(3), Integer(3)); diff(pyitem((v_t * v_X), 0, 1), v_t)})), "X[0, 1]")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Integer(3), Integer(3)); pyitem(v_X, 0, 1)})), "X[0, 1]")
  expect_equal(.bat_fmt(local({v_d <- MatAdd(MatMul(MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(0), Integer(1)), Identity(Integer(3))), MatMul(Trace(MatrixSymbol("X", Integer(3), Integer(3))), MatrixUnit(Integer(3), Integer(3), Integer(0), Integer(1)))); v_d})), "X[0, 1]*I + Trace(X)*MatrixUnit(3, 3, 0, 1)")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Integer(3), Integer(3)); ((pyitem(v_X, 0, 1) * Identity(3)) + (Trace(v_X) * MatrixUnit(3, 3, 0, 1)))})), "X[0, 1]*I + Trace(X)*MatrixUnit(3, 3, 0, 1)")
  expect_equal(.bat_fmt(local({v_d <- MatAdd(MatMul(MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(0), Integer(1)), Identity(Integer(3))), MatMul(Trace(MatrixSymbol("X", Integer(3), Integer(3))), MatrixUnit(Integer(3), Integer(3), Integer(0), Integer(1)))); as_explicit(v_d)})), "Matrix([[X[0, 1], Trace(X), 0], [0, X[0, 1], 0], [0, 0, X[0, 1]]])")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Integer(3), Integer(3)); Matrix(list(list(pyitem(v_X, 0, 1), Trace(v_X), 0), list(0, pyitem(v_X, 0, 1), 0), list(0, 0, pyitem(v_X, 0, 1))))})), "Matrix([[X[0, 1], Trace(X), 0], [0, X[0, 1], 0], [0, 0, X[0, 1]]])")
})

test_that("test_derivatives: test_derivatives_of_scalar_functions_of_traces matches SymPy", {
  expect_equal(.bat_fmt(local({v_A2 <- MatrixSymbol("A", Integer(2), Integer(2)); diff(sin(Trace(v_A2)), v_A2)})), "cos(Trace(A))*I")
  expect_equal(.bat_fmt(local({v_A2 <- MatrixSymbol("A", Integer(2), Integer(2)); (cos(Trace(v_A2)) * Identity(2))})), "cos(Trace(A))*I")
  expect_equal(.bat_fmt(local({v_d <- MatMul(exp(Trace(MatMul(MatrixSymbol("A", Integer(2), Integer(2)), MatrixSymbol("B", Integer(2), Integer(2))))), Transpose(MatrixSymbol("B", Integer(2), Integer(2)))); get_shape(v_d)})), "[2, 2]")
  expect_equal(.bat_fmt(local({list(2, 2)})), "[2, 2]")
})

test_that("test_derivatives: test_derivative_negative_matrix_power matches SymPy", {
  expect_equal(.bat_fmt(local({v_A2 <- MatrixSymbol("A", Integer(2), Integer(2)); as_explicit(diff(MatPow(v_A2, (-1)), v_A2))})), "[[[[-A[1, 1]^2/(A[0, 0]*A[1, 1] - A[0, 1]*A[1, 0])^2, A[0, 1]*A[1, 1]/(A[0, 0]*A[1, 1] - A[0, 1]*A[1, 0])^2], [A[1, 0]*A[1, 1]/(A[0, 0]*A[1, 1] - A[0, 1]*A[1, 0])^2, -A[0, 1]*A[1, 0]/(A[0, 0]*A[1, 1] - A[0, 1]*A[1, 0])^2]], [[A[1, 0]*A[1, 1]/(A[0, 0]*A[1, 1] - A[0, 1]*A[1, 0])^2, -A[0, 0]*A[1, 1]/(A[0, 0]*A[1, 1] - A[0, 1]*A[1, 0])^2], [-A[1, 0]^2/(A[0, 0]*A[1, 1] - A[0, 1]*A[1, 0])^2, A[0, 0]*A[1, 0]/(A[0, 0]*A[1, 1] - A[0, 1]*A[1, 0])^2]]], [[[A[0, 1]*A[1, 1]/(A[0, 0]*A[1, 1] - A[0, 1]*A[1, 0])^2, -A[0, 1]^2/(A[0, 0]*A[1, 1] - A[0, 1]*A[1, 0])^2], [-A[0, 0]*A[1, 1]/(A[0, 0]*A[1, 1] - A[0, 1]*A[1, 0])^2, A[0, 0]*A[0, 1]/(A[0, 0]*A[1, 1] - A[0, 1]*A[1, 0])^2]], [[-A[0, 1]*A[1, 0]/(A[0, 0]*A[1, 1] - A[0, 1]*A[1, 0])^2, A[0, 0]*A[0, 1]/(A[0, 0]*A[1, 1] - A[0, 1]*A[1, 0])^2], [A[0, 0]*A[1, 0]/(A[0, 0]*A[1, 1] - A[0, 1]*A[1, 0])^2, -A[0, 0]^2/(A[0, 0]*A[1, 1] - A[0, 1]*A[1, 0])^2]]]]")
  expect_equal(.bat_fmt(local({v_A2 <- MatrixSymbol("A", Integer(2), Integer(2)); as_explicit(diff(Inverse(v_A2), v_A2))})), "[[[[-A[1, 1]^2/(A[0, 0]*A[1, 1] - A[0, 1]*A[1, 0])^2, A[0, 1]*A[1, 1]/(A[0, 0]*A[1, 1] - A[0, 1]*A[1, 0])^2], [A[1, 0]*A[1, 1]/(A[0, 0]*A[1, 1] - A[0, 1]*A[1, 0])^2, -A[0, 1]*A[1, 0]/(A[0, 0]*A[1, 1] - A[0, 1]*A[1, 0])^2]], [[A[1, 0]*A[1, 1]/(A[0, 0]*A[1, 1] - A[0, 1]*A[1, 0])^2, -A[0, 0]*A[1, 1]/(A[0, 0]*A[1, 1] - A[0, 1]*A[1, 0])^2], [-A[1, 0]^2/(A[0, 0]*A[1, 1] - A[0, 1]*A[1, 0])^2, A[0, 0]*A[1, 0]/(A[0, 0]*A[1, 1] - A[0, 1]*A[1, 0])^2]]], [[[A[0, 1]*A[1, 1]/(A[0, 0]*A[1, 1] - A[0, 1]*A[1, 0])^2, -A[0, 1]^2/(A[0, 0]*A[1, 1] - A[0, 1]*A[1, 0])^2], [-A[0, 0]*A[1, 1]/(A[0, 0]*A[1, 1] - A[0, 1]*A[1, 0])^2, A[0, 0]*A[0, 1]/(A[0, 0]*A[1, 1] - A[0, 1]*A[1, 0])^2]], [[-A[0, 1]*A[1, 0]/(A[0, 0]*A[1, 1] - A[0, 1]*A[1, 0])^2, A[0, 0]*A[0, 1]/(A[0, 0]*A[1, 1] - A[0, 1]*A[1, 0])^2], [A[0, 0]*A[1, 0]/(A[0, 0]*A[1, 1] - A[0, 1]*A[1, 0])^2, -A[0, 0]^2/(A[0, 0]*A[1, 1] - A[0, 1]*A[1, 0])^2]]]]")
  expect_equal(.bat_fmt(local({v_got <- Array(list(Integer(310), Integer(-122), Integer(-183), Integer(72), Integer(-183), Integer(66), Integer(108), Integer(-39), Integer(-122), Integer(48), Integer(66), Integer(-26), Integer(72), Integer(-26), Integer(-39), Integer(14)), c(2, 2, 2, 2)); v_got})), "[[[[310, -122], [-183, 72]], [[-183, 66], [108, -39]]], [[[-122, 48], [66, -26]], [[72, -26], [-39, 14]]]]")
  expect_equal(.bat_fmt(local({v_td <- MatAdd(MatMul(Integer(-1), Inverse(Transpose(MatrixSymbol("A", Integer(2), Integer(2)))), Transpose(MatrixSymbol("B", Integer(2), Integer(2))), MatPow(Transpose(MatrixSymbol("A", Integer(2), Integer(2))), Integer(-2))), MatMul(Integer(-1), MatPow(Transpose(MatrixSymbol("A", Integer(2), Integer(2))), Integer(-2)), Transpose(MatrixSymbol("B", Integer(2), Integer(2))), Inverse(Transpose(MatrixSymbol("A", Integer(2), Integer(2)))))); get_shape(v_td)})), "[2, 2]")
})

test_that("test_derivatives: test_derivative_trace_of_hadamard_shape matches SymPy", {
  expect_equal(.bat_fmt(local({v_d <- HadamardProduct(Identity(Integer(2)), Transpose(MatrixSymbol("A", Integer(2), Integer(2)))); get_shape(v_d)})), "[2, 2]")
  expect_equal(.bat_fmt(local({v_d <- HadamardProduct(Identity(Integer(2)), Transpose(MatrixSymbol("A", Integer(2), Integer(2)))); as_explicit(v_d)})), "Matrix([[A[0, 0], 0], [0, A[1, 1]]])")
  expect_equal(.bat_fmt(local({v_A2 <- MatrixSymbol("A", Integer(2), Integer(2)); Matrix(list(list(pyitem(v_A2, 0, 0), 0), list(0, pyitem(v_A2, 1, 1))))})), "Matrix([[A[0, 0], 0], [0, A[1, 1]]])")
  expect_equal(.bat_fmt(local({v_d2 <- MatMul(Integer(2), HadamardProduct(Identity(Integer(2)), Transpose(MatrixSymbol("X", Integer(2), Integer(2))))); get_shape(v_d2)})), "[2, 2]")
  expect_equal(.bat_fmt(local({v_d2 <- MatMul(Integer(2), HadamardProduct(Identity(Integer(2)), Transpose(MatrixSymbol("X", Integer(2), Integer(2))))); as_explicit(v_d2)})), "Matrix([[2*X[0, 0], 0], [0, 2*X[1, 1]]])")
  expect_equal(.bat_fmt(local({v_X2 <- MatrixSymbol("X", Integer(2), Integer(2)); Matrix(list(list((2 * pyitem(v_X2, 0, 0)), 0), list(0, (2 * pyitem(v_X2, 1, 1)))))})), "Matrix([[2*X[0, 0], 0], [0, 2*X[1, 1]]])")
  expect_equal(.bat_fmt(local({v_d3 <- MatMul(Integer(3), HadamardProduct(HadamardPower(Transpose(MatrixSymbol("X", Integer(2), Integer(2))), Integer(2)), Identity(Integer(2)))); get_shape(v_d3)})), "[2, 2]")
  expect_equal(.bat_fmt(local({v_d3 <- MatMul(Integer(3), HadamardProduct(HadamardPower(Transpose(MatrixSymbol("X", Integer(2), Integer(2))), Integer(2)), Identity(Integer(2)))); as_explicit(v_d3)})), "Matrix([[3*X[0, 0]^2, 0], [0, 3*X[1, 1]^2]])")
  expect_equal(.bat_fmt(local({v_X2 <- MatrixSymbol("X", Integer(2), Integer(2)); Matrix(list(list((3 * (pyitem(v_X2, 0, 0) ^ 2)), 0), list(0, (3 * (pyitem(v_X2, 1, 1) ^ 2)))))})), "Matrix([[3*X[0, 0]^2, 0], [0, 3*X[1, 1]^2]])")
})

test_that("test_derivatives: test_derivative_trace_of_vector_hadamard matches SymPy", {
  expect_equal(.bat_fmt(local({list(((2 * 2) * 5), ((2 * 3) * 7))})), "[20, 42]")
})

test_that("test_derivatives: test_derivative_diag_matrix matches SymPy", {
  expect_equal(.bat_fmt(local({v_d <- Reshape(ArrayDiagonal(ArrayTensorProduct(Identity(Integer(2)), Identity(Integer(1)), Identity(Integer(2))), Tuple(Integer(1), Integer(4))), Tuple(Integer(2), Integer(1), Integer(2), Integer(2))); get_shape(v_d)})), "[2, 1, 2, 2]")
  expect_equal(.bat_fmt(local({list(2, 1, 2, 2)})), "[2, 1, 2, 2]")
  expect_equal(.bat_fmt(local({v_d <- Reshape(ArrayDiagonal(ArrayTensorProduct(Identity(Integer(2)), Identity(Integer(1)), Identity(Integer(2))), Tuple(Integer(1), Integer(4))), Tuple(Integer(2), Integer(1), Integer(2), Integer(2))); as_explicit(v_d)})), "[[[[1, 0], [0, 0]]], [[[0, 0], [0, 1]]]]")
  expect_equal(.bat_fmt(local({.d1 <- Dummy("d"); v_dd <- Reshape(PermuteDims(ArrayDiagonal(ArrayTensorProduct(ElementwiseApplyFunction(Lambda(.d1, (Integer(-1)*sin(.d1))), MatrixSymbol("x", Integer(2), Integer(1))), Identity(Integer(2)), Identity(Integer(1)), Identity(Integer(2))), Tuple(Integer(1), Integer(5)), Tuple(Integer(0), Integer(3), Integer(6))), Permutation(c(0, 1, 3, 2, 4))), Tuple(Integer(2), Integer(1), Integer(2), Integer(2))); get_shape(v_dd)})), "[2, 1, 2, 2]")
  expect_equal(.bat_fmt(local({.d1 <- Dummy("d"); v_dd <- Reshape(PermuteDims(ArrayDiagonal(ArrayTensorProduct(ElementwiseApplyFunction(Lambda(.d1, (Integer(-1)*sin(.d1))), MatrixSymbol("x", Integer(2), Integer(1))), Identity(Integer(2)), Identity(Integer(1)), Identity(Integer(2))), Tuple(Integer(1), Integer(5)), Tuple(Integer(0), Integer(3), Integer(6))), Permutation(c(0, 1, 3, 2, 4))), Tuple(Integer(2), Integer(1), Integer(2), Integer(2))); as_explicit(v_dd)})), "[[[[-sin(x[0, 0]), 0], [0, 0]]], [[[0, 0], [0, -sin(x[1, 0])]]]]")
})

test_that("test_derivatives: test_derivative_hadamard_of_outer_product matches SymPy", {
  expect_equal(.bat_fmt(local({v_got <- Array(list(Integer(10), Integer(14), Integer(21), Integer(0), Integer(0), Integer(10), Integer(15), Integer(70)), c(2, 1, 2, 2)); v_got})), "[[[[10, 14], [21, 0]]], [[[0, 10], [15, 70]]]]")
  expect_equal(.bat_fmt(local({v_truth_n <- Array(list(Integer(10), Integer(14), Integer(21), Integer(0), Integer(0), Integer(10), Integer(15), Integer(70)), c(2, 1, 2, 2)); v_truth_n})), "[[[[10, 14], [21, 0]]], [[[0, 10], [15, 70]]]]")
  expect_equal(.bat_fmt(local({v_x2 <- MatrixSymbol("x", Integer(2), Integer(1)); v_a2 <- MatrixSymbol("a", Integer(2), Integer(1)); v_B2 <- MatrixSymbol("B", Integer(2), Integer(2)); get_shape(diff((HadamardProduct(v_a2, v_x2) + (v_B2 * v_x2)), v_x2))})), "[2, 2]")
})

test_that("test_derivatives: test_derivative_object_routing matches SymPy", {
  expect_equal(.bat_fmt(local({v_A2 <- MatrixSymbol("A", Integer(2), Integer(2)); v_X2 <- MatrixSymbol("X", Integer(2), Integer(2)); doit(Derivative(Trace((v_A2 * v_X2)), v_X2))})), "A.T")
  expect_equal(.bat_fmt(local({v_A2 <- MatrixSymbol("A", Integer(2), Integer(2)); t(v_A2)})), "A.T")
  expect_equal(.bat_fmt(local({v_X2 <- MatrixSymbol("X", Integer(2), Integer(2)); doit(Derivative(Trace((v_X2 * v_X2)), v_X2, 2))})), "PermuteDims(ArrayTensorProduct(2, I, I), (1 2 3))")
  expect_equal(.bat_fmt(local({v_X2 <- MatrixSymbol("X", Integer(2), Integer(2)); diff(diff(Trace((v_X2 * v_X2)), v_X2), v_X2)})), "PermuteDims(ArrayTensorProduct(2, I, I), (1 2 3))")
})

test_that("test_matexpr: test_matrix_symbol_creation matches SymPy", {
  expect_equal(.bat_fmt(local({MatrixSymbol("A", 2, 2)})), "A")
  expect_equal(.bat_fmt(local({MatrixSymbol("A", 0, 0)})), "A")
  expect_error(local({MatrixSymbol("A", (-1), 2)}), "The dimension specification -1 should be a nonnegative integer.", fixed = TRUE, class = "ValueError")
  expect_error(local({MatrixSymbol("A", 2, (-1))}), "The dimension specification -1 should be a nonnegative integer.", fixed = TRUE, class = "ValueError")
  expect_equal(.bat_fmt(local({v_n <- Symbol("n"); MatrixSymbol("A", v_n, v_n)})), "A")
  expect_error(local({v_n <- Symbol("n", integer = FALSE); MatrixSymbol("A", v_n, v_n)}), "The dimension specification n should be a nonnegative integer.", fixed = TRUE, class = "ValueError")
  expect_error(local({v_n <- Symbol("n", negative = TRUE); MatrixSymbol("A", v_n, v_n)}), "The dimension specification n should be a nonnegative integer.", fixed = TRUE, class = "ValueError")
})

test_that("test_matexpr: test_matexpr_properties matches SymPy", {
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n", integer = TRUE), Symbol("m", integer = TRUE)); get_shape(v_A)})), "[n, m]")
  expect_equal(.bat_fmt(local({v_n <- Symbol("n", integer = TRUE); v_m <- Symbol("m", integer = TRUE); list(v_n, v_m)})), "[n, m]")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n", integer = TRUE), Symbol("m", integer = TRUE)); v_B <- MatrixSymbol("B", Symbol("m", integer = TRUE), Symbol("l", integer = TRUE)); get_shape((v_A * v_B))})), "[n, l]")
  expect_equal(.bat_fmt(local({v_n <- Symbol("n", integer = TRUE); v_l <- Symbol("l", integer = TRUE); list(v_n, v_l)})), "[n, l]")
  expect_equal(.bat_fmt(local({list(0, 1)})), "[0, 1]")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n", integer = TRUE), Symbol("m", integer = TRUE)); v_A})), "A")
})

test_that("test_matexpr: test_matexpr matches SymPy", {
  expect_equal(.bat_fmt(local({v_x <- Symbol("x"); v_A <- MatrixSymbol("A", Symbol("n", integer = TRUE), Symbol("m", integer = TRUE)); get_shape((v_x * v_A))})), "[n, m]")
  expect_equal(.bat_fmt(local({v_x <- Symbol("x"); v_A <- MatrixSymbol("A", Symbol("n", integer = TRUE), Symbol("m", integer = TRUE)); sym_type((v_x * v_A))})), "MatMul")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n", integer = TRUE), Symbol("m", integer = TRUE)); (((2 * v_A) - v_A) - v_A)})), "0")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n", integer = TRUE), Symbol("m", integer = TRUE)); do.call(ZeroMatrix, c(get_shape(v_A)))})), "0")
})

test_that("test_matexpr: test_matexpr_subs matches SymPy", {
  expect_equal(.bat_fmt(local({v_n <- Symbol("n", integer = TRUE); v_m <- Symbol("m", integer = TRUE); v_A <- MatrixSymbol("A", Symbol("n", integer = TRUE), Symbol("m", integer = TRUE)); get_shape(subs(v_A, v_n, v_m))})), "[m, m]")
  expect_equal(.bat_fmt(local({v_m <- Symbol("m", integer = TRUE); list(v_m, v_m)})), "[m, m]")
  expect_equal(.bat_fmt(local({v_B <- MatrixSymbol("B", Symbol("m", integer = TRUE), Symbol("l", integer = TRUE)); v_C <- MatrixSymbol("C", Symbol("m", integer = TRUE), Symbol("l", integer = TRUE)); v_A <- MatrixSymbol("A", Symbol("n", integer = TRUE), Symbol("m", integer = TRUE)); subs((v_A * v_B), v_B, v_C)})), "A*C")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n", integer = TRUE), Symbol("m", integer = TRUE)); v_C <- MatrixSymbol("C", Symbol("m", integer = TRUE), Symbol("l", integer = TRUE)); (v_A * v_C)})), "A*C")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Integer(2), Integer(2)); v_Y <- MatrixSymbol("Y", Integer(1), Integer(2)); subs(v_X, v_X, v_Y)})), "Y")
  expect_equal(.bat_fmt(local({v_Y <- MatrixSymbol("Y", Integer(1), Integer(2)); v_Y})), "Y")
  expect_equal(.bat_fmt(local({v_y <- "y"; list(v_y, 2, 2)})), "[y, 2, 2]")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Integer(2), Integer(2)); v_W <- MatrixSymbol("W", Integer(3), Integer(3)); subs(pyitem(v_X, 1, 1), v_X, v_W)})), "W[1, 1]")
  expect_equal(.bat_fmt(local({v_W <- MatrixSymbol("W", Integer(3), Integer(3)); pyitem(v_W, 1, 1)})), "W[1, 1]")
  expect_error(local({v_X <- MatrixSymbol("X", Integer(2), Integer(2)); v_Y <- MatrixSymbol("Y", Integer(1), Integer(2)); subs(pyitem(v_X, 1, 1), v_X, v_Y)}), "indices out of range", fixed = TRUE, class = "IndexError")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Integer(2), Integer(2)); v_Y <- MatrixSymbol("Y", Integer(1), Integer(2)); subs(pyitem(v_X, 0, 1), v_X, v_Y)})), "Y[0, 1]")
  expect_equal(.bat_fmt(local({v_Y <- MatrixSymbol("Y", Integer(1), Integer(2)); pyitem(v_Y, 0, 1)})), "Y[0, 1]")
  expect_equal(.bat_fmt(local({v_W <- MatrixSymbol("W", Integer(3), Integer(3)); v_Z <- MatrixSymbol("Z", Symbol("n", integer = TRUE), Integer(2)); subs(pyitem(v_W, 2, 1), v_W, v_Z)})), "Z[2, 1]")
  expect_equal(.bat_fmt(local({v_Z <- MatrixSymbol("Z", Symbol("n", integer = TRUE), Integer(2)); pyitem(v_Z, 2, 1)})), "Z[2, 1]")
  expect_error(local({v_W <- MatrixSymbol("W", Integer(3), Integer(3)); v_Z <- MatrixSymbol("Z", Symbol("n", integer = TRUE), Integer(2)); subs(pyitem(v_W, 2, 2), v_W, v_Z)}), "indices out of range", fixed = TRUE, class = "IndexError")
  expect_error(local({v_W <- MatrixSymbol("W", Integer(3), Integer(3)); subs(pyitem(v_W, 2, 2), v_W, zeros(2))}), "Index out of range: a[2]", fixed = TRUE, class = "IndexError")
  expect_equal(.bat_fmt(local({v_C <- MatrixSymbol("C", Integer(2), Integer(2)); v_A <- Matrix(list(list(Integer(1), Integer(2)), list(Integer(3), Integer(4)))); v_D <- MatrixSymbol("D", Integer(2), Integer(2)); v_B <- Matrix(list(list(Integer(1), Integer(2)), list(Integer(3), Integer(4)))); subs((v_C * v_D), list(list(v_C, v_A), list(v_D, v_B)))})), "Matrix([\\n[1, 2],\\n[3, 4]])*Matrix([\\n[1, 2],\\n[3, 4]])")
  expect_equal(.bat_fmt(local({v_A <- Matrix(list(list(Integer(1), Integer(2)), list(Integer(3), Integer(4)))); v_B <- Matrix(list(list(Integer(1), Integer(2)), list(Integer(3), Integer(4)))); MatMul(v_A, v_B)})), "Matrix([\\n[1, 2],\\n[3, 4]])*Matrix([\\n[1, 2],\\n[3, 4]])")
})

test_that("test_matexpr: test_addition matches SymPy", {
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n", integer = TRUE), Symbol("m", integer = TRUE)); v_B <- MatrixSymbol("B", Symbol("n", integer = TRUE), Symbol("m", integer = TRUE)); (v_A + v_B)})), "A + B")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n", integer = TRUE), Symbol("m", integer = TRUE)); v_B <- MatrixSymbol("B", Symbol("n", integer = TRUE), Symbol("m", integer = TRUE)); get_shape((v_A + v_B))})), "[n, m]")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n", integer = TRUE), Symbol("m", integer = TRUE)); v_B <- MatrixSymbol("B", Symbol("n", integer = TRUE), Symbol("m", integer = TRUE)); ((v_A - v_A) + (2 * v_B))})), "2*B")
  expect_error(local({v_A <- MatrixSymbol("A", Symbol("n", integer = TRUE), Symbol("m", integer = TRUE)); (v_A + 1)}), "Mix of Matrix and Scalar symbols", fixed = TRUE, class = "TypeError")
  expect_error(local({v_A <- MatrixSymbol("A", Symbol("n", integer = TRUE), Symbol("m", integer = TRUE)); (5 + v_A)}), "Mix of Matrix and Scalar symbols", fixed = TRUE, class = "TypeError")
  expect_error(local({v_A <- MatrixSymbol("A", Symbol("n", integer = TRUE), Symbol("m", integer = TRUE)); (5 - v_A)}), "Mix of Matrix and Scalar symbols", fixed = TRUE, class = "TypeError")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n", integer = TRUE), Symbol("m", integer = TRUE)); v_n <- Symbol("n", integer = TRUE); v_m <- Symbol("m", integer = TRUE); ((v_A + ZeroMatrix(v_n, v_m)) - v_A)})), "0")
  expect_equal(.bat_fmt(local({v_n <- Symbol("n", integer = TRUE); v_m <- Symbol("m", integer = TRUE); ZeroMatrix(v_n, v_m)})), "0")
  expect_error(local({v_n <- Symbol("n", integer = TRUE); v_m <- Symbol("m", integer = TRUE); (ZeroMatrix(v_n, v_m) + S_Zero)}), "Mix of Matrix and Scalar symbols", fixed = TRUE, class = "TypeError")
})

test_that("test_matexpr: test_multiplication matches SymPy", {
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n", integer = TRUE), Symbol("m", integer = TRUE)); v_B <- MatrixSymbol("B", Symbol("m", integer = TRUE), Symbol("l", integer = TRUE)); get_shape(((2 * v_A) * v_B))})), "[n, l]")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n", integer = TRUE), Symbol("m", integer = TRUE)); v_B <- MatrixSymbol("B", Symbol("m", integer = TRUE), Symbol("l", integer = TRUE)); ((v_A * 0) * v_B)})), "0")
  expect_equal(.bat_fmt(local({v_n <- Symbol("n", integer = TRUE); v_l <- Symbol("l", integer = TRUE); ZeroMatrix(v_n, v_l)})), "0")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n", integer = TRUE), Symbol("m", integer = TRUE)); get_shape((2 * v_A))})), "[n, m]")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n", integer = TRUE), Symbol("m", integer = TRUE)); v_m <- Symbol("m", integer = TRUE); v_B <- MatrixSymbol("B", Symbol("m", integer = TRUE), Symbol("l", integer = TRUE)); ((v_A * ZeroMatrix(v_m, v_m)) * v_B)})), "0")
  expect_equal(.bat_fmt(local({v_C <- MatrixSymbol("C", Symbol("n", integer = TRUE), Symbol("n", integer = TRUE)); v_n <- Symbol("n", integer = TRUE); ((v_C * Identity(v_n)) * inverse(v_C))})), "I")
  expect_equal(.bat_fmt(local({v_n <- Symbol("n", integer = TRUE); Identity(v_n)})), "I")
  expect_equal(.bat_fmt(local({v_B <- MatrixSymbol("B", Symbol("m", integer = TRUE), Symbol("l", integer = TRUE)); (v_B / 2)})), "(1/2)*B")
  expect_equal(.bat_fmt(local({v_B <- MatrixSymbol("B", Symbol("m", integer = TRUE), Symbol("l", integer = TRUE)); (S_Half * v_B)})), "(1/2)*B")
  expect_error(local({v_B <- MatrixSymbol("B", Symbol("m", integer = TRUE), Symbol("l", integer = TRUE)); (2 / v_B)}), "", fixed = TRUE, class = "NotImplementedError")
  expect_equal(.bat_fmt(local({v_n <- Symbol("n", integer = TRUE); v_A <- MatrixSymbol("A", Symbol("n", integer = TRUE), Symbol("n", integer = TRUE)); v_B <- MatrixSymbol("B", Symbol("n", integer = TRUE), Symbol("n", integer = TRUE)); (Identity(v_n) * (v_A + v_B))})), "A + B")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n", integer = TRUE), Symbol("n", integer = TRUE)); v_B <- MatrixSymbol("B", Symbol("n", integer = TRUE), Symbol("n", integer = TRUE)); (v_A + v_B)})), "A + B")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n", integer = TRUE), Symbol("n", integer = TRUE)); ((v_A ^ 2) * v_A)})), "A^3")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n", integer = TRUE), Symbol("n", integer = TRUE)); (v_A ^ 3)})), "A^3")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n", integer = TRUE), Symbol("n", integer = TRUE)); ((v_A ^ 2) * (inverse(v_A) ^ 3))})), "A^(-1)")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n", integer = TRUE), Symbol("n", integer = TRUE)); inverse(v_A)})), "A^(-1)")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n", integer = TRUE), Symbol("n", integer = TRUE)); ((v_A ^ 3) * (inverse(v_A) ^ 2))})), "A")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n", integer = TRUE), Symbol("n", integer = TRUE)); v_A})), "A")
})

test_that("test_matexpr: test_MatPow matches SymPy", {
  expect_equal(.bat_fmt(local({v_n <- Symbol("n", integer = TRUE); v_n})), "n")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n", integer = TRUE), Symbol("n", integer = TRUE)); (v_A ^ 0)})), "I")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n", integer = TRUE), Symbol("n", integer = TRUE)); (v_A ^ 1)})), "A")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n", integer = TRUE), Symbol("n", integer = TRUE)); (v_A ^ 2)})), "A^2")
  expect_equal(.bat_fmt(local({v_AA <- MatPow(MatrixSymbol("A", Symbol("n", integer = TRUE), Symbol("n", integer = TRUE)), Integer(2)); v_AA})), "A^2")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n", integer = TRUE), Symbol("n", integer = TRUE)); (v_A ^ (-1))})), "A^(-1)")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n", integer = TRUE), Symbol("n", integer = TRUE)); Inverse(v_A)})), "A^(-1)")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n", integer = TRUE), Symbol("n", integer = TRUE)); ((v_A ^ (-1)) ^ (-1))})), "A")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n", integer = TRUE), Symbol("n", integer = TRUE)); ((v_A ^ 2) ^ 3)})), "A^6")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n", integer = TRUE), Symbol("n", integer = TRUE)); (v_A ^ 6)})), "A^6")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n", integer = TRUE), Symbol("n", integer = TRUE)); (v_A ^ S_Half)})), "A^(1/2)")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n", integer = TRUE), Symbol("n", integer = TRUE)); sqrt(v_A)})), "A^(1/2)")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n", integer = TRUE), Symbol("n", integer = TRUE)); (v_A ^ Rational(1, 3))})), "A^(1/3)")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n", integer = TRUE), Symbol("n", integer = TRUE)); cbrt(v_A)})), "A^(1/3)")
  expect_error(local({(MatrixSymbol("B", 3, 2) ^ 2)}), "Power of non-square matrix B", fixed = TRUE, class = "NonSquareMatrixError")
})

test_that("test_matexpr: test_MatrixSymbol matches SymPy", {
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("n"), Symbol("m")); get_shape(v_X)})), "[n, m]")
  expect_equal(.bat_fmt(local({v_n <- Symbol("n"); v_m <- Symbol("m"); list(v_n, v_m)})), "[n, m]")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("n"), Symbol("m")); doit(v_X)})), "X")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("n"), Symbol("m")); v_X})), "X")
})

test_that("test_matexpr: test_dense_conversion matches SymPy", {
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Integer(2), Integer(2)); Matrix(v_X)})), "Matrix([[X[0, 0], X[0, 1]], [X[1, 0], X[1, 1]]])")
})

test_that("test_matexpr: test_free_symbols matches SymPy", {
  expect_equal(.bat_fmt(local({v_C <- MatrixSymbol("C", Symbol("n", integer = TRUE), Symbol("n", integer = TRUE)); v_D <- MatrixSymbol("D", Symbol("n", integer = TRUE), Symbol("n", integer = TRUE)); pyset(free_symbols((v_C * v_D)))})), "{C, D}")
})

test_that("test_matexpr: test_zero_matmul matches SymPy", {
  expect_equal(.bat_fmt(local({(S_Zero * MatrixSymbol("X", 2, 2))})), "0")
})

test_that("test_matexpr: test_matadd_simplify matches SymPy", {
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Integer(1), Integer(1)); v_x <- Symbol("x"); simplify(MatAdd(v_A, Matrix(list(list(((sin(v_x) ^ 2) + (cos(v_x) ^ 2)))))))})), "Matrix([[1]]) + A")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Integer(1), Integer(1)); MatAdd(v_A, Matrix(list(list(1))))})), "Matrix([[1]]) + A")
})

test_that("test_matexpr: test_matmul_simplify matches SymPy", {
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Integer(1), Integer(1)); v_x <- Symbol("x"); simplify(MatMul(v_A, Matrix(list(list(((sin(v_x) ^ 2) + (cos(v_x) ^ 2)))))))})), "A*Matrix([[1]])")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Integer(1), Integer(1)); MatMul(v_A, Matrix(list(list(1))))})), "A*Matrix([[1]])")
})

test_that("test_matexpr: test_single_indexing matches SymPy", {
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Integer(2), Integer(3)); pyitem(v_A, 1)})), "A[0, 1]")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Integer(2), Integer(3)); pyitem(v_A, 0, 1)})), "A[0, 1]")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Integer(2), Integer(3)); pyitem(v_A, 3)})), "A[1, 0]")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Integer(2), Integer(3)); pyitem(v_A, 1, 0)})), "A[1, 0]")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Integer(2), Integer(3)); list(pyitem(v_A, 0, 0), pyitem(v_A, 0, 1), pyitem(v_A, 1, 0), pyitem(v_A, 1, 1))})), "[A[0, 0], A[0, 1], A[1, 0], A[1, 1]]")
  expect_error(local({v_A <- MatrixSymbol("A", Integer(2), Integer(3)); pyitem(v_A, 6)}), "Invalid index 6", fixed = TRUE, class = "IndexError")
  expect_error(local({v_A <- MatrixSymbol("A", Integer(2), Integer(3)); v_n <- Symbol("n", integer = TRUE); pyitem(v_A, v_n)}), "\nOnly integers may be used when addressing the matrix with a single\nindex.", fixed = TRUE, class = "IndexError")
  expect_error(local({v_B <- MatrixSymbol("B", Symbol("n", integer = TRUE), Symbol("m", integer = TRUE)); pyitem(v_B, 1)}), "\nSingle indexing is only supported when the number of columns is known.", fixed = TRUE, class = "IndexError")
  expect_equal(.bat_fmt(local({v_B <- MatrixSymbol("B", Symbol("n", integer = TRUE), Integer(3)); pyitem(v_B, 3)})), "B[1, 0]")
  expect_equal(.bat_fmt(local({v_B <- MatrixSymbol("B", Symbol("n", integer = TRUE), Integer(3)); pyitem(v_B, 1, 0)})), "B[1, 0]")
})

test_that("test_matexpr: test_MatrixElement_commutative matches SymPy", {
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n", integer = TRUE), Symbol("m", integer = TRUE)); (pyitem(v_A, 0, 1) * pyitem(v_A, 1, 0))})), "A[0, 1]*A[1, 0]")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n", integer = TRUE), Symbol("m", integer = TRUE)); (pyitem(v_A, 1, 0) * pyitem(v_A, 0, 1))})), "A[0, 1]*A[1, 0]")
})

test_that("test_matexpr: test_MatrixSymbol_determinant matches SymPy", {
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Integer(4), Integer(4)); sym_det(as_explicit(v_A))})), "A[0, 0]*A[1, 1]*A[2, 2]*A[3, 3] - A[0, 0]*A[1, 1]*A[2, 3]*A[3, 2] - A[0, 0]*A[1, 2]*A[2, 1]*A[3, 3] + A[0, 0]*A[1, 2]*A[2, 3]*A[3, 1] + A[0, 0]*A[1, 3]*A[2, 1]*A[3, 2] - A[0, 0]*A[1, 3]*A[2, 2]*A[3, 1] - A[0, 1]*A[1, 0]*A[2, 2]*A[3, 3] + A[0, 1]*A[1, 0]*A[2, 3]*A[3, 2] + A[0, 1]*A[1, 2]*A[2, 0]*A[3, 3] - A[0, 1]*A[1, 2]*A[2, 3]*A[3, 0] - A[0, 1]*A[1, 3]*A[2, 0]*A[3, 2] + A[0, 1]*A[1, 3]*A[2, 2]*A[3, 0] + A[0, 2]*A[1, 0]*A[2, 1]*A[3, 3] - A[0, 2]*A[1, 0]*A[2, 3]*A[3, 1] - A[0, 2]*A[1, 1]*A[2, 0]*A[3, 3] + A[0, 2]*A[1, 1]*A[2, 3]*A[3, 0] + A[0, 2]*A[1, 3]*A[2, 0]*A[3, 1] - A[0, 2]*A[1, 3]*A[2, 1]*A[3, 0] - A[0, 3]*A[1, 0]*A[2, 1]*A[3, 2] + A[0, 3]*A[1, 0]*A[2, 2]*A[3, 1] + A[0, 3]*A[1, 1]*A[2, 0]*A[3, 2] - A[0, 3]*A[1, 1]*A[2, 2]*A[3, 0] - A[0, 3]*A[1, 2]*A[2, 0]*A[3, 1] + A[0, 3]*A[1, 2]*A[2, 1]*A[3, 0]")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Integer(4), Integer(4)); ((((((((((((((((((((((((((pyitem(v_A, 0, 0) * pyitem(v_A, 1, 1)) * pyitem(v_A, 2, 2)) * pyitem(v_A, 3, 3)) - (((pyitem(v_A, 0, 0) * pyitem(v_A, 1, 1)) * pyitem(v_A, 2, 3)) * pyitem(v_A, 3, 2))) - (((pyitem(v_A, 0, 0) * pyitem(v_A, 1, 2)) * pyitem(v_A, 2, 1)) * pyitem(v_A, 3, 3))) + (((pyitem(v_A, 0, 0) * pyitem(v_A, 1, 2)) * pyitem(v_A, 2, 3)) * pyitem(v_A, 3, 1))) + (((pyitem(v_A, 0, 0) * pyitem(v_A, 1, 3)) * pyitem(v_A, 2, 1)) * pyitem(v_A, 3, 2))) - (((pyitem(v_A, 0, 0) * pyitem(v_A, 1, 3)) * pyitem(v_A, 2, 2)) * pyitem(v_A, 3, 1))) - (((pyitem(v_A, 0, 1) * pyitem(v_A, 1, 0)) * pyitem(v_A, 2, 2)) * pyitem(v_A, 3, 3))) + (((pyitem(v_A, 0, 1) * pyitem(v_A, 1, 0)) * pyitem(v_A, 2, 3)) * pyitem(v_A, 3, 2))) + (((pyitem(v_A, 0, 1) * pyitem(v_A, 1, 2)) * pyitem(v_A, 2, 0)) * pyitem(v_A, 3, 3))) - (((pyitem(v_A, 0, 1) * pyitem(v_A, 1, 2)) * pyitem(v_A, 2, 3)) * pyitem(v_A, 3, 0))) - (((pyitem(v_A, 0, 1) * pyitem(v_A, 1, 3)) * pyitem(v_A, 2, 0)) * pyitem(v_A, 3, 2))) + (((pyitem(v_A, 0, 1) * pyitem(v_A, 1, 3)) * pyitem(v_A, 2, 2)) * pyitem(v_A, 3, 0))) + (((pyitem(v_A, 0, 2) * pyitem(v_A, 1, 0)) * pyitem(v_A, 2, 1)) * pyitem(v_A, 3, 3))) - (((pyitem(v_A, 0, 2) * pyitem(v_A, 1, 0)) * pyitem(v_A, 2, 3)) * pyitem(v_A, 3, 1))) - (((pyitem(v_A, 0, 2) * pyitem(v_A, 1, 1)) * pyitem(v_A, 2, 0)) * pyitem(v_A, 3, 3))) + (((pyitem(v_A, 0, 2) * pyitem(v_A, 1, 1)) * pyitem(v_A, 2, 3)) * pyitem(v_A, 3, 0))) + (((pyitem(v_A, 0, 2) * pyitem(v_A, 1, 3)) * pyitem(v_A, 2, 0)) * pyitem(v_A, 3, 1))) - (((pyitem(v_A, 0, 2) * pyitem(v_A, 1, 3)) * pyitem(v_A, 2, 1)) * pyitem(v_A, 3, 0))) - (((pyitem(v_A, 0, 3) * pyitem(v_A, 1, 0)) * pyitem(v_A, 2, 1)) * pyitem(v_A, 3, 2))) + (((pyitem(v_A, 0, 3) * pyitem(v_A, 1, 0)) * pyitem(v_A, 2, 2)) * pyitem(v_A, 3, 1))) + (((pyitem(v_A, 0, 3) * pyitem(v_A, 1, 1)) * pyitem(v_A, 2, 0)) * pyitem(v_A, 3, 2))) - (((pyitem(v_A, 0, 3) * pyitem(v_A, 1, 1)) * pyitem(v_A, 2, 2)) * pyitem(v_A, 3, 0))) - (((pyitem(v_A, 0, 3) * pyitem(v_A, 1, 2)) * pyitem(v_A, 2, 0)) * pyitem(v_A, 3, 1))) + (((pyitem(v_A, 0, 3) * pyitem(v_A, 1, 2)) * pyitem(v_A, 2, 1)) * pyitem(v_A, 3, 0)))})), "A[0, 0]*A[1, 1]*A[2, 2]*A[3, 3] - A[0, 0]*A[1, 1]*A[2, 3]*A[3, 2] - A[0, 0]*A[1, 2]*A[2, 1]*A[3, 3] + A[0, 0]*A[1, 2]*A[2, 3]*A[3, 1] + A[0, 0]*A[1, 3]*A[2, 1]*A[3, 2] - A[0, 0]*A[1, 3]*A[2, 2]*A[3, 1] - A[0, 1]*A[1, 0]*A[2, 2]*A[3, 3] + A[0, 1]*A[1, 0]*A[2, 3]*A[3, 2] + A[0, 1]*A[1, 2]*A[2, 0]*A[3, 3] - A[0, 1]*A[1, 2]*A[2, 3]*A[3, 0] - A[0, 1]*A[1, 3]*A[2, 0]*A[3, 2] + A[0, 1]*A[1, 3]*A[2, 2]*A[3, 0] + A[0, 2]*A[1, 0]*A[2, 1]*A[3, 3] - A[0, 2]*A[1, 0]*A[2, 3]*A[3, 1] - A[0, 2]*A[1, 1]*A[2, 0]*A[3, 3] + A[0, 2]*A[1, 1]*A[2, 3]*A[3, 0] + A[0, 2]*A[1, 3]*A[2, 0]*A[3, 1] - A[0, 2]*A[1, 3]*A[2, 1]*A[3, 0] - A[0, 3]*A[1, 0]*A[2, 1]*A[3, 2] + A[0, 3]*A[1, 0]*A[2, 2]*A[3, 1] + A[0, 3]*A[1, 1]*A[2, 0]*A[3, 2] - A[0, 3]*A[1, 1]*A[2, 2]*A[3, 0] - A[0, 3]*A[1, 2]*A[2, 0]*A[3, 1] + A[0, 3]*A[1, 2]*A[2, 1]*A[3, 0]")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Integer(4), Integer(4)); v_B <- MatrixSymbol("B", Integer(4), Integer(4)); doit(Determinant((v_A + v_B)))})), "Determinant(A + B)")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Integer(4), Integer(4)); v_B <- MatrixSymbol("B", Integer(4), Integer(4)); sym_det((v_A + v_B))})), "Determinant(A + B)")
})

test_that("test_matexpr: test_MatrixElement_diff matches SymPy", {
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n", integer = TRUE), Symbol("m", integer = TRUE)); diff((pyitem(v_A, 3, 0) * pyitem(v_A, 0, 0)), pyitem(v_A, 0, 0))})), "A[3, 0]")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n", integer = TRUE), Symbol("m", integer = TRUE)); pyitem(v_A, 3, 0)})), "A[3, 0]")
})

test_that("test_matexpr: test_MatrixElement_doit matches SymPy", {
  expect_equal(.bat_fmt(local({v_u <- MatrixSymbol("u", Integer(2), Integer(1)); v_v <- Matrix(list(list(Integer(3)), list(Integer(5)))); doit(subs(pyitem(v_u, 0, 0), v_u, v_v))})), "3")
  expect_equal(.bat_fmt(local({v_v <- Matrix(list(list(Integer(3)), list(Integer(5)))); pyitem(v_v, 0, 0)})), "3")
})

test_that("test_matexpr: test_identity_powers matches SymPy", {
  expect_equal(.bat_fmt(local({v_M <- Identity(Symbol("n", integer = TRUE)); doit(MatPow(v_M, 3))})), "I")
  expect_equal(.bat_fmt(local({v_M <- Identity(Symbol("n", integer = TRUE)); (v_M ^ 3)})), "I")
  expect_equal(.bat_fmt(local({v_M <- Identity(Symbol("n", integer = TRUE)); v_n <- Symbol("n", integer = TRUE); (v_M ^ v_n)})), "I")
  expect_equal(.bat_fmt(local({v_M <- Identity(Symbol("n", integer = TRUE)); v_M})), "I")
  expect_equal(.bat_fmt(local({v_M <- Identity(Symbol("n", integer = TRUE)); doit(MatPow(v_M, 0))})), "I")
  expect_equal(.bat_fmt(local({v_M <- Identity(Symbol("n", integer = TRUE)); (v_M ^ 2)})), "I")
  expect_equal(.bat_fmt(local({v_M <- Identity(Symbol("n", integer = TRUE)); (v_M ^ (-2))})), "I")
  expect_equal(.bat_fmt(local({v_M <- Identity(Symbol("n", integer = TRUE)); doit(MatPow(v_M, (-2)))})), "I")
  expect_equal(.bat_fmt(local({v_M <- Identity(Symbol("n", integer = TRUE)); (v_M ^ 0)})), "I")
  expect_equal(.bat_fmt(local({v_N <- Identity(Integer(3)); doit(MatPow(v_N, 2))})), "I")
  expect_equal(.bat_fmt(local({v_N <- Identity(Integer(3)); v_n <- Symbol("n", integer = TRUE); (v_N ^ v_n)})), "I")
  expect_equal(.bat_fmt(local({v_N <- Identity(Integer(3)); doit(MatPow(v_N, 3))})), "I")
  expect_equal(.bat_fmt(local({v_N <- Identity(Integer(3)); v_N})), "I")
  expect_equal(.bat_fmt(local({v_N <- Identity(Integer(3)); doit(MatPow(v_N, (-2)))})), "I")
  expect_equal(.bat_fmt(local({v_N <- Identity(Integer(3)); (v_N ^ 4)})), "I")
  expect_equal(.bat_fmt(local({v_N <- Identity(Integer(3)); (v_N ^ 0)})), "I")
})

test_that("test_matexpr: test_Zero_power matches SymPy", {
  expect_equal(.bat_fmt(local({v_z1 <- ZeroMatrix(Symbol("n", integer = TRUE), Symbol("n", integer = TRUE)); (v_z1 ^ 4)})), "0")
  expect_equal(.bat_fmt(local({v_z1 <- ZeroMatrix(Symbol("n", integer = TRUE), Symbol("n", integer = TRUE)); v_z1})), "0")
  expect_error(local({v_z1 <- ZeroMatrix(Symbol("n", integer = TRUE), Symbol("n", integer = TRUE)); (v_z1 ^ (-2))}), "Matrix det == 0; not invertible", fixed = TRUE, class = "NonInvertibleMatrixError")
  expect_equal(.bat_fmt(local({v_z1 <- ZeroMatrix(Symbol("n", integer = TRUE), Symbol("n", integer = TRUE)); (v_z1 ^ 0)})), "I")
  expect_equal(.bat_fmt(local({v_z1 <- ZeroMatrix(Symbol("n", integer = TRUE), Symbol("n", integer = TRUE)); doit(MatPow(v_z1, 2))})), "0")
  expect_equal(.bat_fmt(local({v_z1 <- ZeroMatrix(Symbol("n", integer = TRUE), Symbol("n", integer = TRUE)); (v_z1 ^ 2)})), "0")
  expect_error(local({v_z1 <- ZeroMatrix(Symbol("n", integer = TRUE), Symbol("n", integer = TRUE)); doit(MatPow(v_z1, (-2)))}), "Matrix det == 0; not invertible", fixed = TRUE, class = "NonInvertibleMatrixError")
  expect_equal(.bat_fmt(local({v_z2 <- ZeroMatrix(Integer(3), Integer(3)); doit(MatPow(v_z2, 4))})), "0")
  expect_equal(.bat_fmt(local({v_z2 <- ZeroMatrix(Integer(3), Integer(3)); (v_z2 ^ 4)})), "0")
  expect_error(local({v_z2 <- ZeroMatrix(Integer(3), Integer(3)); (v_z2 ^ (-3))}), "Matrix det == 0; not invertible", fixed = TRUE, class = "NonInvertibleMatrixError")
  expect_equal(.bat_fmt(local({v_z2 <- ZeroMatrix(Integer(3), Integer(3)); (v_z2 ^ 3)})), "0")
  expect_equal(.bat_fmt(local({v_z2 <- ZeroMatrix(Integer(3), Integer(3)); doit(MatPow(v_z2, 3))})), "0")
  expect_equal(.bat_fmt(local({v_z2 <- ZeroMatrix(Integer(3), Integer(3)); (v_z2 ^ 0)})), "I")
  expect_error(local({v_z2 <- ZeroMatrix(Integer(3), Integer(3)); doit(MatPow(v_z2, (-1)))}), "Matrix det == 0; not invertible.", fixed = TRUE, class = "NonInvertibleMatrixError")
})

test_that("test_matexpr: test_matrixelement_diff matches SymPy", {
  expect_equal(.bat_fmt(local({v_w <- MatrixSymbol("w", Symbol("n", integer = TRUE), Integer(1)); v_k <- Symbol("k", integer = TRUE); v_p <- Symbol("p", integer = TRUE); diff(pyitem(v_w, v_k, v_p), pyitem(v_w, v_k, v_p))})), "1")
  expect_equal(.bat_fmt(local({v_w <- MatrixSymbol("w", Symbol("n", integer = TRUE), Integer(1)); v_k <- Symbol("k", integer = TRUE); v_p <- Symbol("p", integer = TRUE); diff(pyitem(v_w, v_k, v_p), pyitem(v_w, 0, 0))})), "KroneckerDelta(0, k, (0, n - 1))*KroneckerDelta(0, p, (0, 0))")
  expect_equal(.bat_fmt(local({v_k <- Symbol("k", integer = TRUE); v_n <- Symbol("n", integer = TRUE); v_p <- Symbol("p", integer = TRUE); (KroneckerDelta(0, v_k, list(0, (v_n - 1))) * KroneckerDelta(0, v_p, list(0, 0)))})), "KroneckerDelta(0, k, (0, n - 1))*KroneckerDelta(0, p, (0, 0))")
  expect_equal(.bat_fmt(local({.d1 <- Dummy("i_1"); v_dexpr <- Sum((KroneckerDelta(.d1, Symbol("p", integer = TRUE), Tuple(Integer(0), (Integer(-1) + Symbol("n", integer = TRUE))))*MatrixElement(MatrixSymbol("D", Symbol("n", integer = TRUE), Symbol("n", integer = TRUE)), Symbol("k", integer = TRUE), .d1)), list(.d1, Integer(0), (Integer(-1) + Symbol("n", integer = TRUE)))); v_dexpr})), "Sum(KroneckerDelta(_i_1, p, (0, n - 1))*D[k, _i_1], (_i_1, 0, n - 1))")
  expect_equal(.bat_fmt(local({.d1 <- Dummy("_i_1"); v__i_1 <- .d1; v_p <- Symbol("p", integer = TRUE); v_n <- Symbol("n", integer = TRUE); v_D <- MatrixSymbol("D", Symbol("n", integer = TRUE), Symbol("n", integer = TRUE)); v_k <- Symbol("k", integer = TRUE); Sum((KroneckerDelta(v__i_1, v_p, list(0, (v_n - 1))) * pyitem(v_D, v_k, v__i_1)), list(v__i_1, 0, (v_n - 1)))})), "Sum(KroneckerDelta(__i_1, p, (0, n - 1))*D[k, __i_1], (__i_1, 0, n - 1))")
  expect_equal(.bat_fmt(local({.d1 <- Dummy("i_1"); v_dexpr <- Sum((KroneckerDelta(.d1, Symbol("p", integer = TRUE), Tuple(Integer(0), (Integer(-1) + Symbol("n", integer = TRUE))))*MatrixElement(MatrixSymbol("D", Symbol("n", integer = TRUE), Symbol("n", integer = TRUE)), Symbol("k", integer = TRUE), .d1)), list(.d1, Integer(0), (Integer(-1) + Symbol("n", integer = TRUE)))); doit(v_dexpr)})), "D[k, p]")
  expect_equal(.bat_fmt(local({v_D <- MatrixSymbol("D", Symbol("n", integer = TRUE), Symbol("n", integer = TRUE)); v_k <- Symbol("k", integer = TRUE); v_p <- Symbol("p", integer = TRUE); pyitem(v_D, v_k, v_p)})), "D[k, p]")
})

test_that("test_matexpr: test_MatrixElement_with_values matches SymPy", {
  expect_equal(.bat_fmt(local({v_Mij <- MatrixElement(Matrix(list(list(Symbol("x"), Symbol("y")), list(Symbol("z"), Symbol("w")))), Symbol("i"), Symbol("j")); v_Mij})), "Matrix([\\n[x, y],\\n[z, w]])[i, j]")
  expect_equal(.bat_fmt(local({v_msij <- MatrixElement(Matrix(list(list(Integer(2), Integer(3)), list(Integer(4), Integer(5)))), Symbol("i"), Symbol("j")); v_msij})), "Matrix([\\n[2, 3],\\n[4, 5]])[i, j]")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Integer(2), Integer(2)); v_M <- Matrix(list(list(Symbol("x"), Symbol("y")), list(Symbol("z"), Symbol("w")))); subs(pyitem(v_A, 0, 0), v_A, v_M)})), "x")
  expect_equal(.bat_fmt(local({v_x <- Symbol("x"); v_x})), "x")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Integer(2), Integer(2)); v_M <- Matrix(list(list(Symbol("x"), Symbol("y")), list(Symbol("z"), Symbol("w")))); v_i <- Symbol("i"); v_j <- Symbol("j"); subs(pyitem(v_A, v_i, v_j), v_A, v_M)})), "Matrix([\\n[x, y],\\n[z, w]])[i, j]")
  expect_equal(.bat_fmt(local({v_M <- Matrix(list(list(Symbol("x"), Symbol("y")), list(Symbol("z"), Symbol("w")))); v_i <- Symbol("i"); v_j <- Symbol("j"); pyitem(v_M, v_i, v_j)})), "Matrix([\\n[x, y],\\n[z, w]])[i, j]")
  expect_equal(.bat_fmt(local({v_M <- Matrix(list(list(Symbol("x"), Symbol("y")), list(Symbol("z"), Symbol("w")))); v_A <- MatrixSymbol("A", Integer(2), Integer(2)); v_i <- Symbol("i"); v_j <- Symbol("j"); subs(pyitem(v_M, v_i, v_j), v_M, v_A)})), "A[i, j]")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Integer(2), Integer(2)); v_i <- Symbol("i"); v_j <- Symbol("j"); pyitem(v_A, v_i, v_j)})), "A[i, j]")
  expect_equal(.bat_fmt(local({v_M <- Matrix(list(list(Symbol("x"), Symbol("y")), list(Symbol("z"), Symbol("w")))); v_i <- Symbol("i"); v_j <- Symbol("j"); pyitem(v_M, ((3 * v_i) - 2), v_j)})), "Matrix([\\n[x, y],\\n[z, w]])[3*i - 2, j]")
  expect_equal(.bat_fmt(local({v_i <- Symbol("i"); v_j <- Symbol("j"); v_M <- Matrix(list(list(Symbol("x"), Symbol("y")), list(Symbol("z"), Symbol("w")))); subs(pyitem(v_M, ((3 * v_i) - 2), v_j), list(list(v_i, 1), list(v_j, 0)))})), "z")
  expect_equal(.bat_fmt(local({v_M <- Matrix(list(list(Symbol("x"), Symbol("y")), list(Symbol("z"), Symbol("w")))); pyitem(v_M, 1, 0)})), "z")
  expect_equal(.bat_fmt(local({v_M <- Matrix(list(list(Symbol("x"), Symbol("y")), list(Symbol("z"), Symbol("w")))); v_i <- Symbol("i"); pyitem(v_M, v_i, 0)})), "Matrix([\\n[x, y],\\n[z, w]])[i, 0]")
  expect_equal(.bat_fmt(local({v_i <- Symbol("i"); v_M <- Matrix(list(list(Symbol("x"), Symbol("y")), list(Symbol("z"), Symbol("w")))); subs(pyitem(v_M, v_i, 0), v_i, 0)})), "x")
  expect_equal(.bat_fmt(local({v_M <- Matrix(list(list(Symbol("x"), Symbol("y")), list(Symbol("z"), Symbol("w")))); pyitem(v_M, 0, 0)})), "x")
  expect_equal(.bat_fmt(local({v_i <- Symbol("i"); v_M <- Matrix(list(list(Symbol("x"), Symbol("y")), list(Symbol("z"), Symbol("w")))); subs(pyitem(v_M, 0, v_i), v_i, 1)})), "y")
  expect_equal(.bat_fmt(local({v_M <- Matrix(list(list(Symbol("x"), Symbol("y")), list(Symbol("z"), Symbol("w")))); pyitem(v_M, 0, 1)})), "y")
  expect_equal(.bat_fmt(local({v_x <- Symbol("x"); v_M <- Matrix(list(list(Symbol("x"), Symbol("y")), list(Symbol("z"), Symbol("w")))); v_i <- Symbol("i"); v_j <- Symbol("j"); diff(pyitem(v_M, v_i, v_j), v_x)})), "Matrix([\\n[1, 0],\\n[0, 0]])[i, j]")
  expect_equal(.bat_fmt(local({v_i <- Symbol("i"); v_j <- Symbol("j"); pyitem(Matrix(list(list(1, 0), list(0, 0))), v_i, v_j)})), "Matrix([\\n[1, 0],\\n[0, 0]])[i, j]")
  expect_error(local({v_M <- Matrix(list(list(Symbol("x"), Symbol("y")), list(Symbol("z"), Symbol("w")))); v_i <- Symbol("i"); pyitem(v_M, v_i, 2)}), "index out of boundary", fixed = TRUE, class = "ValueError")
  expect_error(local({v_M <- Matrix(list(list(Symbol("x"), Symbol("y")), list(Symbol("z"), Symbol("w")))); v_i <- Symbol("i"); pyitem(v_M, v_i, (-1))}), "index out of boundary", fixed = TRUE, class = "ValueError")
  expect_error(local({v_M <- Matrix(list(list(Symbol("x"), Symbol("y")), list(Symbol("z"), Symbol("w")))); v_i <- Symbol("i"); pyitem(v_M, 2, v_i)}), "index out of boundary", fixed = TRUE, class = "ValueError")
  expect_error(local({v_M <- Matrix(list(list(Symbol("x"), Symbol("y")), list(Symbol("z"), Symbol("w")))); v_i <- Symbol("i"); pyitem(v_M, (-1), v_i)}), "index out of boundary", fixed = TRUE, class = "ValueError")
})

test_that("test_matexpr: test_inv matches SymPy", {
  expect_equal(.bat_fmt(local({v_B <- MatrixSymbol("B", Integer(3), Integer(3)); inverse(v_B)})), "B^(-1)")
  expect_equal(.bat_fmt(local({v_B <- MatrixSymbol("B", Integer(3), Integer(3)); (v_B ^ (-1))})), "B^(-1)")
  expect_equal(.bat_fmt(local({v_X <- Matrix(list(list(MatrixElement(MatrixSymbol("X", Integer(1), Integer(1)), Integer(0), Integer(0))))); inverse(v_X)})), "Matrix([[1/X[0, 0]]])")
  expect_equal(.bat_fmt(local({v_X <- Matrix(list(list(MatrixElement(MatrixSymbol("X", Integer(1), Integer(1)), Integer(0), Integer(0))))); Matrix(list(list((1 / pyitem(v_X, 0, 0)))))})), "Matrix([[1/X[0, 0]]])")
  expect_equal(.bat_fmt(local({v_X <- Matrix(list(list(MatrixElement(MatrixSymbol("X", Integer(2), Integer(2)), Integer(0), Integer(0)), MatrixElement(MatrixSymbol("X", Integer(2), Integer(2)), Integer(0), Integer(1))), list(MatrixElement(MatrixSymbol("X", Integer(2), Integer(2)), Integer(1), Integer(0)), MatrixElement(MatrixSymbol("X", Integer(2), Integer(2)), Integer(1), Integer(1))))); inverse(v_X)})), "Matrix([[X[1, 1]/(X[0, 0]*X[1, 1] - X[0, 1]*X[1, 0]), -X[0, 1]/(X[0, 0]*X[1, 1] - X[0, 1]*X[1, 0])], [-X[1, 0]/(X[0, 0]*X[1, 1] - X[0, 1]*X[1, 0]), X[0, 0]/(X[0, 0]*X[1, 1] - X[0, 1]*X[1, 0])]])")
  expect_equal(.bat_fmt(local({v_invX <- Matrix(list(list(((((MatrixElement(MatrixSymbol("X", Integer(2), Integer(2)), Integer(0), Integer(0))*MatrixElement(MatrixSymbol("X", Integer(2), Integer(2)), Integer(1), Integer(1))) + (Integer(-1)*MatrixElement(MatrixSymbol("X", Integer(2), Integer(2)), Integer(0), Integer(1))*MatrixElement(MatrixSymbol("X", Integer(2), Integer(2)), Integer(1), Integer(0)))))^(Integer(-1))*MatrixElement(MatrixSymbol("X", Integer(2), Integer(2)), Integer(1), Integer(1))), (Integer(-1)*(((MatrixElement(MatrixSymbol("X", Integer(2), Integer(2)), Integer(0), Integer(0))*MatrixElement(MatrixSymbol("X", Integer(2), Integer(2)), Integer(1), Integer(1))) + (Integer(-1)*MatrixElement(MatrixSymbol("X", Integer(2), Integer(2)), Integer(0), Integer(1))*MatrixElement(MatrixSymbol("X", Integer(2), Integer(2)), Integer(1), Integer(0)))))^(Integer(-1))*MatrixElement(MatrixSymbol("X", Integer(2), Integer(2)), Integer(0), Integer(1)))), list((Integer(-1)*(((MatrixElement(MatrixSymbol("X", Integer(2), Integer(2)), Integer(0), Integer(0))*MatrixElement(MatrixSymbol("X", Integer(2), Integer(2)), Integer(1), Integer(1))) + (Integer(-1)*MatrixElement(MatrixSymbol("X", Integer(2), Integer(2)), Integer(0), Integer(1))*MatrixElement(MatrixSymbol("X", Integer(2), Integer(2)), Integer(1), Integer(0)))))^(Integer(-1))*MatrixElement(MatrixSymbol("X", Integer(2), Integer(2)), Integer(1), Integer(0))), ((((MatrixElement(MatrixSymbol("X", Integer(2), Integer(2)), Integer(0), Integer(0))*MatrixElement(MatrixSymbol("X", Integer(2), Integer(2)), Integer(1), Integer(1))) + (Integer(-1)*MatrixElement(MatrixSymbol("X", Integer(2), Integer(2)), Integer(0), Integer(1))*MatrixElement(MatrixSymbol("X", Integer(2), Integer(2)), Integer(1), Integer(0)))))^(Integer(-1))*MatrixElement(MatrixSymbol("X", Integer(2), Integer(2)), Integer(0), Integer(0)))))); v_invX})), "Matrix([[X[1, 1]/(X[0, 0]*X[1, 1] - X[0, 1]*X[1, 0]), -X[0, 1]/(X[0, 0]*X[1, 1] - X[0, 1]*X[1, 0])], [-X[1, 0]/(X[0, 0]*X[1, 1] - X[0, 1]*X[1, 0]), X[0, 0]/(X[0, 0]*X[1, 1] - X[0, 1]*X[1, 0])]])")
})

test_that("test_matexpr: test_factor_expand matches SymPy", {
  expect_equal(.bat_fmt(local({v_expr1 <- MatMul(MatAdd(MatrixSymbol("A", Symbol("n", integer = TRUE), Symbol("n", integer = TRUE)), MatrixSymbol("B", Symbol("n", integer = TRUE), Symbol("n", integer = TRUE))), MatAdd(MatrixSymbol("C", Symbol("n", integer = TRUE), Symbol("n", integer = TRUE)), MatrixSymbol("D", Symbol("n", integer = TRUE), Symbol("n", integer = TRUE)))); v_expr1})), "(A + B)*(C + D)")
  expect_equal(.bat_fmt(local({v_expr2 <- MatAdd(MatMul(MatrixSymbol("A", Symbol("n", integer = TRUE), Symbol("n", integer = TRUE)), MatrixSymbol("C", Symbol("n", integer = TRUE), Symbol("n", integer = TRUE))), MatMul(MatrixSymbol("A", Symbol("n", integer = TRUE), Symbol("n", integer = TRUE)), MatrixSymbol("D", Symbol("n", integer = TRUE), Symbol("n", integer = TRUE))), MatMul(MatrixSymbol("B", Symbol("n", integer = TRUE), Symbol("n", integer = TRUE)), MatrixSymbol("C", Symbol("n", integer = TRUE), Symbol("n", integer = TRUE))), MatMul(MatrixSymbol("B", Symbol("n", integer = TRUE), Symbol("n", integer = TRUE)), MatrixSymbol("D", Symbol("n", integer = TRUE), Symbol("n", integer = TRUE)))); v_expr2})), "A*C + A*D + B*C + B*D")
  expect_equal(.bat_fmt(local({v_I <- Identity(Symbol("n", integer = TRUE)); v_C <- MatrixSymbol("C", Symbol("n", integer = TRUE), Symbol("n", integer = TRUE)); v_B <- MatrixSymbol("B", Symbol("n", integer = TRUE), Symbol("n", integer = TRUE)); v_A <- MatrixSymbol("A", Symbol("n", integer = TRUE), Symbol("n", integer = TRUE)); list((v_I - v_C), (((v_B ^ (-1)) * ((((v_A ^ (-1)) * (v_I - v_C)) * (v_B ^ (-1))) ^ (-1))) * (v_A ^ (-1))))})), "[I - C, (I - C)^(-1)]")
})

test_that("test_matexpr: test_issue_2749 matches SymPy", {
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Integer(5), Integer(2)); as_explicit(inverse((t(v_A) * v_A)))})), "Matrix([[(A[0, 1]^2 + A[1, 1]^2 + A[2, 1]^2 + A[3, 1]^2 + A[4, 1]^2)/(A[0, 0]^2*A[1, 1]^2 + A[0, 0]^2*A[2, 1]^2 + A[0, 0]^2*A[3, 1]^2 + A[0, 0]^2*A[4, 1]^2 - 2*A[0, 0]*A[0, 1]*A[1, 0]*A[1, 1] - 2*A[0, 0]*A[0, 1]*A[2, 0]*A[2, 1] - 2*A[0, 0]*A[0, 1]*A[3, 0]*A[3, 1] - 2*A[0, 0]*A[0, 1]*A[4, 0]*A[4, 1] + A[0, 1]^2*A[1, 0]^2 + A[0, 1]^2*A[2, 0]^2 + A[0, 1]^2*A[3, 0]^2 + A[0, 1]^2*A[4, 0]^2 + A[1, 0]^2*A[2, 1]^2 + A[1, 0]^2*A[3, 1]^2 + A[1, 0]^2*A[4, 1]^2 - 2*A[1, 0]*A[1, 1]*A[2, 0]*A[2, 1] - 2*A[1, 0]*A[1, 1]*A[3, 0]*A[3, 1] - 2*A[1, 0]*A[1, 1]*A[4, 0]*A[4, 1] + A[1, 1]^2*A[2, 0]^2 + A[1, 1]^2*A[3, 0]^2 + A[1, 1]^2*A[4, 0]^2 + A[2, 0]^2*A[3, 1]^2 + A[2, 0]^2*A[4, 1]^2 - 2*A[2, 0]*A[2, 1]*A[3, 0]*A[3, 1] - 2*A[2, 0]*A[2, 1]*A[4, 0]*A[4, 1] + A[2, 1]^2*A[3, 0]^2 + A[2, 1]^2*A[4, 0]^2 + A[3, 0]^2*A[4, 1]^2 - 2*A[3, 0]*A[3, 1]*A[4, 0]*A[4, 1] + A[3, 1]^2*A[4, 0]^2), (-A[0, 0]*A[0, 1] - A[1, 0]*A[1, 1] - A[2, 0]*A[2, 1] - A[3, 0]*A[3, 1] - A[4, 0]*A[4, 1])/(A[0, 0]^2*A[1, 1]^2 + A[0, 0]^2*A[2, 1]^2 + A[0, 0]^2*A[3, 1]^2 + A[0, 0]^2*A[4, 1]^2 - 2*A[0, 0]*A[0, 1]*A[1, 0]*A[1, 1] - 2*A[0, 0]*A[0, 1]*A[2, 0]*A[2, 1] - 2*A[0, 0]*A[0, 1]*A[3, 0]*A[3, 1] - 2*A[0, 0]*A[0, 1]*A[4, 0]*A[4, 1] + A[0, 1]^2*A[1, 0]^2 + A[0, 1]^2*A[2, 0]^2 + A[0, 1]^2*A[3, 0]^2 + A[0, 1]^2*A[4, 0]^2 + A[1, 0]^2*A[2, 1]^2 + A[1, 0]^2*A[3, 1]^2 + A[1, 0]^2*A[4, 1]^2 - 2*A[1, 0]*A[1, 1]*A[2, 0]*A[2, 1] - 2*A[1, 0]*A[1, 1]*A[3, 0]*A[3, 1] - 2*A[1, 0]*A[1, 1]*A[4, 0]*A[4, 1] + A[1, 1]^2*A[2, 0]^2 + A[1, 1]^2*A[3, 0]^2 + A[1, 1]^2*A[4, 0]^2 + A[2, 0]^2*A[3, 1]^2 + A[2, 0]^2*A[4, 1]^2 - 2*A[2, 0]*A[2, 1]*A[3, 0]*A[3, 1] - 2*A[2, 0]*A[2, 1]*A[4, 0]*A[4, 1] + A[2, 1]^2*A[3, 0]^2 + A[2, 1]^2*A[4, 0]^2 + A[3, 0]^2*A[4, 1]^2 - 2*A[3, 0]*A[3, 1]*A[4, 0]*A[4, 1] + A[3, 1]^2*A[4, 0]^2)], [(-A[0, 0]*A[0, 1] - A[1, 0]*A[1, 1] - A[2, 0]*A[2, 1] - A[3, 0]*A[3, 1] - A[4, 0]*A[4, 1])/(A[0, 0]^2*A[1, 1]^2 + A[0, 0]^2*A[2, 1]^2 + A[0, 0]^2*A[3, 1]^2 + A[0, 0]^2*A[4, 1]^2 - 2*A[0, 0]*A[0, 1]*A[1, 0]*A[1, 1] - 2*A[0, 0]*A[0, 1]*A[2, 0]*A[2, 1] - 2*A[0, 0]*A[0, 1]*A[3, 0]*A[3, 1] - 2*A[0, 0]*A[0, 1]*A[4, 0]*A[4, 1] + A[0, 1]^2*A[1, 0]^2 + A[0, 1]^2*A[2, 0]^2 + A[0, 1]^2*A[3, 0]^2 + A[0, 1]^2*A[4, 0]^2 + A[1, 0]^2*A[2, 1]^2 + A[1, 0]^2*A[3, 1]^2 + A[1, 0]^2*A[4, 1]^2 - 2*A[1, 0]*A[1, 1]*A[2, 0]*A[2, 1] - 2*A[1, 0]*A[1, 1]*A[3, 0]*A[3, 1] - 2*A[1, 0]*A[1, 1]*A[4, 0]*A[4, 1] + A[1, 1]^2*A[2, 0]^2 + A[1, 1]^2*A[3, 0]^2 + A[1, 1]^2*A[4, 0]^2 + A[2, 0]^2*A[3, 1]^2 + A[2, 0]^2*A[4, 1]^2 - 2*A[2, 0]*A[2, 1]*A[3, 0]*A[3, 1] - 2*A[2, 0]*A[2, 1]*A[4, 0]*A[4, 1] + A[2, 1]^2*A[3, 0]^2 + A[2, 1]^2*A[4, 0]^2 + A[3, 0]^2*A[4, 1]^2 - 2*A[3, 0]*A[3, 1]*A[4, 0]*A[4, 1] + A[3, 1]^2*A[4, 0]^2), (A[0, 0]^2 + A[1, 0]^2 + A[2, 0]^2 + A[3, 0]^2 + A[4, 0]^2)/(A[0, 0]^2*A[1, 1]^2 + A[0, 0]^2*A[2, 1]^2 + A[0, 0]^2*A[3, 1]^2 + A[0, 0]^2*A[4, 1]^2 - 2*A[0, 0]*A[0, 1]*A[1, 0]*A[1, 1] - 2*A[0, 0]*A[0, 1]*A[2, 0]*A[2, 1] - 2*A[0, 0]*A[0, 1]*A[3, 0]*A[3, 1] - 2*A[0, 0]*A[0, 1]*A[4, 0]*A[4, 1] + A[0, 1]^2*A[1, 0]^2 + A[0, 1]^2*A[2, 0]^2 + A[0, 1]^2*A[3, 0]^2 + A[0, 1]^2*A[4, 0]^2 + A[1, 0]^2*A[2, 1]^2 + A[1, 0]^2*A[3, 1]^2 + A[1, 0]^2*A[4, 1]^2 - 2*A[1, 0]*A[1, 1]*A[2, 0]*A[2, 1] - 2*A[1, 0]*A[1, 1]*A[3, 0]*A[3, 1] - 2*A[1, 0]*A[1, 1]*A[4, 0]*A[4, 1] + A[1, 1]^2*A[2, 0]^2 + A[1, 1]^2*A[3, 0]^2 + A[1, 1]^2*A[4, 0]^2 + A[2, 0]^2*A[3, 1]^2 + A[2, 0]^2*A[4, 1]^2 - 2*A[2, 0]*A[2, 1]*A[3, 0]*A[3, 1] - 2*A[2, 0]*A[2, 1]*A[4, 0]*A[4, 1] + A[2, 1]^2*A[3, 0]^2 + A[2, 1]^2*A[4, 0]^2 + A[3, 0]^2*A[4, 1]^2 - 2*A[3, 0]*A[3, 1]*A[4, 0]*A[4, 1] + A[3, 1]^2*A[4, 0]^2)]])")
})

test_that("test_matexpr: test_issue_2750 matches SymPy", {
  expect_equal(.bat_fmt(local({v_x <- MatrixSymbol("x", Integer(1), Integer(1)); (as_explicit((t(v_x) * v_x)) ^ (-1))})), "Matrix([[x[0, 0]^(-2)]])")
  expect_equal(.bat_fmt(local({v_x <- MatrixSymbol("x", Integer(1), Integer(1)); Matrix(list(list((pyitem(v_x, 0, 0) ^ (-2)))))})), "Matrix([[x[0, 0]^(-2)]])")
})

test_that("test_matexpr: test_issue_7842 matches SymPy", {
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Integer(3), Integer(1)); v_B <- MatrixSymbol("B", Integer(2), Integer(1)); Eq(v_A, v_B)})), "False")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Integer(3), Integer(1)); v_B <- MatrixSymbol("B", Integer(2), Integer(1)); sym_type(Eq(pyitem(v_A, 1, 0), pyitem(v_B, 1, 0)))})), "Equality")
  expect_equal(.bat_fmt(local({v_A <- ZeroMatrix(Integer(2), Integer(3)); v_B <- ZeroMatrix(Integer(2), Integer(3)); Eq(v_A, v_B)})), "True")
})

test_that("test_matexpr: test_issue_21195 matches SymPy", {
  expect_equal(.bat_fmt(local({v_x <- Function("x")(Symbol("t")); v_A <- Matrix(list(list(((Derivative(Function("x")(Symbol("t")), Symbol("t"))*cos(Function("x")(Symbol("t")))) + cos(Function("x")(Symbol("t"))))), list(((Derivative(Function("x")(Symbol("t")), list(Symbol("t"), Integer(2)))*tan(Function("x")(Symbol("t")))) + sin(Function("x")(Symbol("t"))))), list((Derivative(Function("x")(Symbol("t")), list(Symbol("t"), Integer(3)))*sin(Symbol("t"))*sin(Function("x")(Symbol("t"))))))); diff(v_A, v_x)})), "Matrix([[-sin(x(t))*Derivative(x(t), t) - sin(x(t))], [(tan(x(t))^2 + 1)*Derivative(x(t), (t, 2)) + cos(x(t))], [sin(t)*cos(x(t))*Derivative(x(t), (t, 3))]])")
  expect_equal(.bat_fmt(local({v_B <- Matrix(list(list(((Integer(-1)*sin(Function("x")(Symbol("t")))) + (Integer(-1)*Derivative(Function("x")(Symbol("t")), Symbol("t"))*sin(Function("x")(Symbol("t")))))), list((((Integer(1) + (tan(Function("x")(Symbol("t"))))^(Integer(2)))*Derivative(Function("x")(Symbol("t")), list(Symbol("t"), Integer(2)))) + cos(Function("x")(Symbol("t"))))), list((Derivative(Function("x")(Symbol("t")), list(Symbol("t"), Integer(3)))*cos(Function("x")(Symbol("t")))*sin(Symbol("t")))))); v_B})), "Matrix([[-sin(x(t))*Derivative(x(t), t) - sin(x(t))], [(tan(x(t))^2 + 1)*Derivative(x(t), (t, 2)) + cos(x(t))], [sin(t)*cos(x(t))*Derivative(x(t), (t, 3))]])")
})

test_that("test_matexpr: test_issue_24859 matches SymPy", {
  expect_equal(.bat_fmt(local({v_x <- Symbol("x"); v_Jk <- Matrix(list(list(((MatrixElement(MatAdd(MatMul(Symbol("x"), MatrixSymbol("u", Integer(2), Integer(3))), MatrixSymbol("A", Integer(2), Integer(3))), Integer(1), Integer(0))*MatrixElement(MatrixSymbol("B", Integer(3), Integer(2)), Integer(0), Integer(1))) + (MatrixElement(MatAdd(MatMul(Symbol("x"), MatrixSymbol("u", Integer(2), Integer(3))), MatrixSymbol("A", Integer(2), Integer(3))), Integer(1), Integer(1))*MatrixElement(MatrixSymbol("B", Integer(3), Integer(2)), Integer(1), Integer(1))) + (MatrixElement(MatAdd(MatMul(Symbol("x"), MatrixSymbol("u", Integer(2), Integer(3))), MatrixSymbol("A", Integer(2), Integer(3))), Integer(1), Integer(2))*MatrixElement(MatrixSymbol("B", Integer(3), Integer(2)), Integer(2), Integer(1)))), ((Integer(-1)*MatrixElement(MatAdd(MatMul(Symbol("x"), MatrixSymbol("u", Integer(2), Integer(3))), MatrixSymbol("A", Integer(2), Integer(3))), Integer(0), Integer(0))*MatrixElement(MatrixSymbol("B", Integer(3), Integer(2)), Integer(0), Integer(1))) + (Integer(-1)*MatrixElement(MatAdd(MatMul(Symbol("x"), MatrixSymbol("u", Integer(2), Integer(3))), MatrixSymbol("A", Integer(2), Integer(3))), Integer(0), Integer(1))*MatrixElement(MatrixSymbol("B", Integer(3), Integer(2)), Integer(1), Integer(1))) + (Integer(-1)*MatrixElement(MatAdd(MatMul(Symbol("x"), MatrixSymbol("u", Integer(2), Integer(3))), MatrixSymbol("A", Integer(2), Integer(3))), Integer(0), Integer(2))*MatrixElement(MatrixSymbol("B", Integer(3), Integer(2)), Integer(2), Integer(1))))), list(((Integer(-1)*MatrixElement(MatAdd(MatMul(Symbol("x"), MatrixSymbol("u", Integer(2), Integer(3))), MatrixSymbol("A", Integer(2), Integer(3))), Integer(1), Integer(0))*MatrixElement(MatrixSymbol("B", Integer(3), Integer(2)), Integer(0), Integer(0))) + (Integer(-1)*MatrixElement(MatAdd(MatMul(Symbol("x"), MatrixSymbol("u", Integer(2), Integer(3))), MatrixSymbol("A", Integer(2), Integer(3))), Integer(1), Integer(1))*MatrixElement(MatrixSymbol("B", Integer(3), Integer(2)), Integer(1), Integer(0))) + (Integer(-1)*MatrixElement(MatAdd(MatMul(Symbol("x"), MatrixSymbol("u", Integer(2), Integer(3))), MatrixSymbol("A", Integer(2), Integer(3))), Integer(1), Integer(2))*MatrixElement(MatrixSymbol("B", Integer(3), Integer(2)), Integer(2), Integer(0)))), ((MatrixElement(MatAdd(MatMul(Symbol("x"), MatrixSymbol("u", Integer(2), Integer(3))), MatrixSymbol("A", Integer(2), Integer(3))), Integer(0), Integer(0))*MatrixElement(MatrixSymbol("B", Integer(3), Integer(2)), Integer(0), Integer(0))) + (MatrixElement(MatAdd(MatMul(Symbol("x"), MatrixSymbol("u", Integer(2), Integer(3))), MatrixSymbol("A", Integer(2), Integer(3))), Integer(0), Integer(1))*MatrixElement(MatrixSymbol("B", Integer(3), Integer(2)), Integer(1), Integer(0))) + (MatrixElement(MatAdd(MatMul(Symbol("x"), MatrixSymbol("u", Integer(2), Integer(3))), MatrixSymbol("A", Integer(2), Integer(3))), Integer(0), Integer(2))*MatrixElement(MatrixSymbol("B", Integer(3), Integer(2)), Integer(2), Integer(0))))))); diff(pyitem(v_Jk, 0, 0), v_x)})), "B[0, 1]*u[1, 0] + B[1, 1]*u[1, 1] + B[2, 1]*u[1, 2]")
  expect_equal(.bat_fmt(local({v_expected <- ((MatrixElement(MatrixSymbol("B", Integer(3), Integer(2)), Integer(0), Integer(1))*MatrixElement(MatrixSymbol("u", Integer(2), Integer(3)), Integer(1), Integer(0))) + (MatrixElement(MatrixSymbol("B", Integer(3), Integer(2)), Integer(1), Integer(1))*MatrixElement(MatrixSymbol("u", Integer(2), Integer(3)), Integer(1), Integer(1))) + (MatrixElement(MatrixSymbol("B", Integer(3), Integer(2)), Integer(2), Integer(1))*MatrixElement(MatrixSymbol("u", Integer(2), Integer(3)), Integer(1), Integer(2)))); v_expected})), "B[0, 1]*u[1, 0] + B[1, 1]*u[1, 1] + B[2, 1]*u[1, 2]")
  expect_equal(.bat_fmt(local({v_Jk <- Matrix(list(list(((MatrixElement(MatAdd(MatMul(Symbol("x"), MatrixSymbol("u", Integer(2), Integer(3))), MatrixSymbol("A", Integer(2), Integer(3))), Integer(1), Integer(0))*MatrixElement(MatrixSymbol("B", Integer(3), Integer(2)), Integer(0), Integer(1))) + (MatrixElement(MatAdd(MatMul(Symbol("x"), MatrixSymbol("u", Integer(2), Integer(3))), MatrixSymbol("A", Integer(2), Integer(3))), Integer(1), Integer(1))*MatrixElement(MatrixSymbol("B", Integer(3), Integer(2)), Integer(1), Integer(1))) + (MatrixElement(MatAdd(MatMul(Symbol("x"), MatrixSymbol("u", Integer(2), Integer(3))), MatrixSymbol("A", Integer(2), Integer(3))), Integer(1), Integer(2))*MatrixElement(MatrixSymbol("B", Integer(3), Integer(2)), Integer(2), Integer(1)))), ((Integer(-1)*MatrixElement(MatAdd(MatMul(Symbol("x"), MatrixSymbol("u", Integer(2), Integer(3))), MatrixSymbol("A", Integer(2), Integer(3))), Integer(0), Integer(0))*MatrixElement(MatrixSymbol("B", Integer(3), Integer(2)), Integer(0), Integer(1))) + (Integer(-1)*MatrixElement(MatAdd(MatMul(Symbol("x"), MatrixSymbol("u", Integer(2), Integer(3))), MatrixSymbol("A", Integer(2), Integer(3))), Integer(0), Integer(1))*MatrixElement(MatrixSymbol("B", Integer(3), Integer(2)), Integer(1), Integer(1))) + (Integer(-1)*MatrixElement(MatAdd(MatMul(Symbol("x"), MatrixSymbol("u", Integer(2), Integer(3))), MatrixSymbol("A", Integer(2), Integer(3))), Integer(0), Integer(2))*MatrixElement(MatrixSymbol("B", Integer(3), Integer(2)), Integer(2), Integer(1))))), list(((Integer(-1)*MatrixElement(MatAdd(MatMul(Symbol("x"), MatrixSymbol("u", Integer(2), Integer(3))), MatrixSymbol("A", Integer(2), Integer(3))), Integer(1), Integer(0))*MatrixElement(MatrixSymbol("B", Integer(3), Integer(2)), Integer(0), Integer(0))) + (Integer(-1)*MatrixElement(MatAdd(MatMul(Symbol("x"), MatrixSymbol("u", Integer(2), Integer(3))), MatrixSymbol("A", Integer(2), Integer(3))), Integer(1), Integer(1))*MatrixElement(MatrixSymbol("B", Integer(3), Integer(2)), Integer(1), Integer(0))) + (Integer(-1)*MatrixElement(MatAdd(MatMul(Symbol("x"), MatrixSymbol("u", Integer(2), Integer(3))), MatrixSymbol("A", Integer(2), Integer(3))), Integer(1), Integer(2))*MatrixElement(MatrixSymbol("B", Integer(3), Integer(2)), Integer(2), Integer(0)))), ((MatrixElement(MatAdd(MatMul(Symbol("x"), MatrixSymbol("u", Integer(2), Integer(3))), MatrixSymbol("A", Integer(2), Integer(3))), Integer(0), Integer(0))*MatrixElement(MatrixSymbol("B", Integer(3), Integer(2)), Integer(0), Integer(0))) + (MatrixElement(MatAdd(MatMul(Symbol("x"), MatrixSymbol("u", Integer(2), Integer(3))), MatrixSymbol("A", Integer(2), Integer(3))), Integer(0), Integer(1))*MatrixElement(MatrixSymbol("B", Integer(3), Integer(2)), Integer(1), Integer(0))) + (MatrixElement(MatAdd(MatMul(Symbol("x"), MatrixSymbol("u", Integer(2), Integer(3))), MatrixSymbol("A", Integer(2), Integer(3))), Integer(0), Integer(2))*MatrixElement(MatrixSymbol("B", Integer(3), Integer(2)), Integer(2), Integer(0))))))); v_x <- Symbol("x"); doit(diff(pyitem(v_Jk, 0, 0), v_x))})), "B[0, 1]*u[1, 0] + B[1, 1]*u[1, 1] + B[2, 1]*u[1, 2]")
})

test_that("test_matexpr: test_MatMul_postprocessor matches SymPy", {
  expect_equal(.bat_fmt(local({v_z <- Matrix(list(list(Integer(0), Integer(0)), list(Integer(0), Integer(0)))); Mul(0, v_z)})), "Matrix([[0, 0], [0, 0]])")
  expect_equal(.bat_fmt(local({v_z <- Matrix(list(list(Integer(0), Integer(0)), list(Integer(0), Integer(0)))); Mul(v_z, 0)})), "Matrix([[0, 0], [0, 0]])")
  expect_equal(.bat_fmt(local({v_z <- Matrix(list(list(Integer(0), Integer(0)), list(Integer(0), Integer(0)))); v_z1 <- ZeroMatrix(Integer(2), Integer(2)); list(v_z, v_z1)})), "[Matrix([[0, 0], [0, 0]]), 0]")
  expect_equal(.bat_fmt(local({v_x <- Symbol("x"); v_M <- Matrix(list(list(Integer(1), Integer(2)), list(Integer(3), Integer(4)))); Mul(v_x, v_M)})), "Matrix([[x, 2*x], [3*x, 4*x]])")
  expect_equal(.bat_fmt(local({v_M <- Matrix(list(list(Integer(1), Integer(2)), list(Integer(3), Integer(4)))); v_x <- Symbol("x"); Mul(v_M, v_x)})), "Matrix([[x, 2*x], [3*x, 4*x]])")
  expect_equal(.bat_fmt(local({v_Mx <- Matrix(list(list(Symbol("x"), (Integer(2)*Symbol("x"))), list((Integer(3)*Symbol("x")), (Integer(4)*Symbol("x"))))); v_Mx})), "Matrix([[x, 2*x], [3*x, 4*x]])")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Integer(2), Integer(2)); v_M <- Matrix(list(list(Integer(1), Integer(2)), list(Integer(3), Integer(4)))); Mul(v_A, v_M)})), "A*Matrix([\\n[1, 2],\\n[3, 4]])")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Integer(2), Integer(2)); v_M <- Matrix(list(list(Integer(1), Integer(2)), list(Integer(3), Integer(4)))); MatMul(v_A, v_M)})), "A*Matrix([\\n[1, 2],\\n[3, 4]])")
  expect_equal(.bat_fmt(local({v_M <- Matrix(list(list(Integer(1), Integer(2)), list(Integer(3), Integer(4)))); v_A <- MatrixSymbol("A", Integer(2), Integer(2)); Mul(v_M, v_A)})), "Matrix([\\n[1, 2],\\n[3, 4]])*A")
  expect_equal(.bat_fmt(local({v_M <- Matrix(list(list(Integer(1), Integer(2)), list(Integer(3), Integer(4)))); v_A <- MatrixSymbol("A", Integer(2), Integer(2)); MatMul(v_M, v_A)})), "Matrix([\\n[1, 2],\\n[3, 4]])*A")
  expect_equal(.bat_fmt(local({v_a <- MatMul(Matrix(list(list(Symbol("x"), (Integer(2)*Symbol("x"))), list((Integer(3)*Symbol("x")), (Integer(4)*Symbol("x"))))), MatrixSymbol("A", Integer(2), Integer(2))); v_a})), "Matrix([\\n[  x, 2*x],\\n[3*x, 4*x]])*A")
  expect_equal(.bat_fmt(local({v_b <- MatMul(Matrix(list(list(Symbol("x"), (Integer(2)*Symbol("x"))), list((Integer(3)*Symbol("x")), (Integer(4)*Symbol("x"))))), MatrixSymbol("A", Integer(2), Integer(2))); v_b})), "Matrix([\\n[  x, 2*x],\\n[3*x, 4*x]])*A")
  expect_equal(.bat_fmt(local({v_c <- MatMul(Matrix(list(list(Symbol("x"), (Integer(2)*Symbol("x"))), list((Integer(3)*Symbol("x")), (Integer(4)*Symbol("x"))))), MatrixSymbol("A", Integer(2), Integer(2))); v_c})), "Matrix([\\n[  x, 2*x],\\n[3*x, 4*x]])*A")
  expect_equal(.bat_fmt(local({v_Mx <- Matrix(list(list(Symbol("x"), (Integer(2)*Symbol("x"))), list((Integer(3)*Symbol("x")), (Integer(4)*Symbol("x"))))); v_A <- MatrixSymbol("A", Integer(2), Integer(2)); MatMul(v_Mx, v_A)})), "Matrix([\\n[  x, 2*x],\\n[3*x, 4*x]])*A")
  expect_equal(.bat_fmt(local({v_a <- MatMul(MatrixSymbol("A", Integer(2), Integer(2)), Matrix(list(list(Symbol("x"), (Integer(2)*Symbol("x"))), list((Integer(3)*Symbol("x")), (Integer(4)*Symbol("x")))))); v_a})), "A*Matrix([\\n[  x, 2*x],\\n[3*x, 4*x]])")
  expect_equal(.bat_fmt(local({v_b <- MatMul(MatrixSymbol("A", Integer(2), Integer(2)), Matrix(list(list(Symbol("x"), (Integer(2)*Symbol("x"))), list((Integer(3)*Symbol("x")), (Integer(4)*Symbol("x")))))); v_b})), "A*Matrix([\\n[  x, 2*x],\\n[3*x, 4*x]])")
  expect_equal(.bat_fmt(local({v_c <- MatMul(MatrixSymbol("A", Integer(2), Integer(2)), Matrix(list(list(Symbol("x"), (Integer(2)*Symbol("x"))), list((Integer(3)*Symbol("x")), (Integer(4)*Symbol("x")))))); v_c})), "A*Matrix([\\n[  x, 2*x],\\n[3*x, 4*x]])")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Integer(2), Integer(2)); v_Mx <- Matrix(list(list(Symbol("x"), (Integer(2)*Symbol("x"))), list((Integer(3)*Symbol("x")), (Integer(4)*Symbol("x"))))); MatMul(v_A, v_Mx)})), "A*Matrix([\\n[  x, 2*x],\\n[3*x, 4*x]])")
  expect_equal(.bat_fmt(local({v_M <- Matrix(list(list(Integer(1), Integer(2)), list(Integer(3), Integer(4)))); Mul(v_M, v_M)})), "Matrix([[7, 10], [15, 22]])")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Integer(2), Integer(2)); v_M <- Matrix(list(list(Integer(1), Integer(2)), list(Integer(3), Integer(4)))); Mul(v_A, v_M, v_M)})), "A*Matrix([\\n[ 7, 10],\\n[15, 22]])")
  expect_equal(.bat_fmt(local({v_M <- Matrix(list(list(Integer(1), Integer(2)), list(Integer(3), Integer(4)))); v_A <- MatrixSymbol("A", Integer(2), Integer(2)); Mul(v_M, v_M, v_A)})), "Matrix([\\n[ 7, 10],\\n[15, 22]])*A")
  expect_equal(.bat_fmt(local({v_M <- Matrix(list(list(Integer(1), Integer(2)), list(Integer(3), Integer(4)))); v_A <- MatrixSymbol("A", Integer(2), Integer(2)); Mul(v_M, v_A, v_M)})), "Matrix([\\n[1, 2],\\n[3, 4]])*A*Matrix([\\n[1, 2],\\n[3, 4]])")
  expect_equal(.bat_fmt(local({v_M <- Matrix(list(list(Integer(1), Integer(2)), list(Integer(3), Integer(4)))); v_A <- MatrixSymbol("A", Integer(2), Integer(2)); MatMul(v_M, v_A, v_M)})), "Matrix([\\n[1, 2],\\n[3, 4]])*A*Matrix([\\n[1, 2],\\n[3, 4]])")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Integer(2), Integer(2)); v_x <- Symbol("x"); v_M <- Matrix(list(list(Integer(1), Integer(2)), list(Integer(3), Integer(4)))); Mul(v_A, v_x, v_M, v_M, v_x)})), "A*Matrix([\\n[ 7*x^2, 10*x^2],\\n[15*x^2, 22*x^2]])")
})

test_that("test_matexpr: test_MatAdd_postprocessor_xfail matches SymPy", {
  expect_equal(.bat_fmt(local({v_z <- Matrix(list(list(Integer(0), Integer(0)), list(Integer(0), Integer(0)))); Add(v_z, nan)})), "nan + Matrix([\\n[0, 0],\\n[0, 0]])")
})

test_that("test_matexpr: test_MatAdd_postprocessor matches SymPy", {
  expect_equal(.bat_fmt(local({v_z <- Matrix(list(list(Integer(0), Integer(0)), list(Integer(0), Integer(0)))); Add(0, v_z)})), "Matrix([[0, 0], [0, 0]])")
  expect_equal(.bat_fmt(local({v_z <- Matrix(list(list(Integer(0), Integer(0)), list(Integer(0), Integer(0)))); Add(v_z, 0)})), "Matrix([[0, 0], [0, 0]])")
  expect_equal(.bat_fmt(local({v_z <- Matrix(list(list(Integer(0), Integer(0)), list(Integer(0), Integer(0)))); v_z})), "Matrix([[0, 0], [0, 0]])")
  expect_equal(.bat_fmt(local({v_z <- Matrix(list(list(Integer(0), Integer(0)), list(Integer(0), Integer(0)))); Add(v_z, oo)})), "oo + Matrix([\\n[0, 0],\\n[0, 0]])")
  expect_equal(.bat_fmt(local({v_z <- Matrix(list(list(Integer(0), Integer(0)), list(Integer(0), Integer(0)))); list(oo, v_z)})), "[oo, Matrix([[0, 0], [0, 0]])]")
  expect_equal(.bat_fmt(local({v_z <- Matrix(list(list(Integer(0), Integer(0)), list(Integer(0), Integer(0)))); Add(v_z, zoo)})), "zoo + Matrix([\\n[0, 0],\\n[0, 0]])")
  expect_equal(.bat_fmt(local({v_z <- Matrix(list(list(Integer(0), Integer(0)), list(Integer(0), Integer(0)))); list(zoo, v_z)})), "[zoo, Matrix([[0, 0], [0, 0]])]")
  expect_equal(.bat_fmt(local({v_z <- Matrix(list(list(Integer(0), Integer(0)), list(Integer(0), Integer(0)))); list(nan, v_z)})), "[nan, Matrix([[0, 0], [0, 0]])]")
  expect_equal(.bat_fmt(local({v_M <- Matrix(list(list(Integer(1), Integer(2)), list(Integer(3), Integer(4)))); v_x <- Symbol("x"); Add(v_M, v_x)})), "x + Matrix([\\n[1, 2],\\n[3, 4]])")
  expect_equal(.bat_fmt(local({v_x <- Symbol("x"); v_M <- Matrix(list(list(Integer(1), Integer(2)), list(Integer(3), Integer(4)))); list(v_x, v_M)})), "[x, Matrix([[1, 2], [3, 4]])]")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Integer(2), Integer(2)); v_M <- Matrix(list(list(Integer(1), Integer(2)), list(Integer(3), Integer(4)))); Add(v_A, v_M)})), "Matrix([\\n[1, 2],\\n[3, 4]]) + A")
  expect_equal(.bat_fmt(local({v_M <- Matrix(list(list(Integer(1), Integer(2)), list(Integer(3), Integer(4)))); v_A <- MatrixSymbol("A", Integer(2), Integer(2)); Add(v_M, v_A)})), "Matrix([\\n[1, 2],\\n[3, 4]]) + A")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Integer(2), Integer(2)); v_M <- Matrix(list(list(Integer(1), Integer(2)), list(Integer(3), Integer(4)))); (v_A + v_M)})), "Matrix([\\n[1, 2],\\n[3, 4]]) + A")
  expect_equal(.bat_fmt(local({v_M <- Matrix(list(list(Integer(1), Integer(2)), list(Integer(3), Integer(4)))); v_x <- Symbol("x"); v_A <- MatrixSymbol("A", Integer(2), Integer(2)); Add(v_M, v_x, v_A)})), "x + (Matrix([\\n[1, 2],\\n[3, 4]]) + A)")
  expect_equal(.bat_fmt(local({v_M <- Matrix(list(list(Integer(1), Integer(2)), list(Integer(3), Integer(4)))); v_A <- MatrixSymbol("A", Integer(2), Integer(2)); v_x <- Symbol("x"); Add(v_M, v_A, v_x)})), "x + (Matrix([\\n[1, 2],\\n[3, 4]]) + A)")
  expect_equal(.bat_fmt(local({v_x <- Symbol("x"); v_A <- MatrixSymbol("A", Integer(2), Integer(2)); v_M <- Matrix(list(list(Integer(1), Integer(2)), list(Integer(3), Integer(4)))); Add(v_x, v_A, v_M)})), "x + (Matrix([\\n[1, 2],\\n[3, 4]]) + A)")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Integer(2), Integer(2)); v_x <- Symbol("x"); v_M <- Matrix(list(list(Integer(1), Integer(2)), list(Integer(3), Integer(4)))); Add(v_A, v_x, v_M)})), "x + (Matrix([\\n[1, 2],\\n[3, 4]]) + A)")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Integer(2), Integer(2)); v_M <- Matrix(list(list(Integer(1), Integer(2)), list(Integer(3), Integer(4)))); v_x <- Symbol("x"); Add(v_A, v_M, v_x)})), "x + (Matrix([\\n[1, 2],\\n[3, 4]]) + A)")
  expect_equal(.bat_fmt(local({v_x <- Symbol("x"); v_A <- MatrixSymbol("A", Integer(2), Integer(2)); v_M <- Matrix(list(list(Integer(1), Integer(2)), list(Integer(3), Integer(4)))); list(v_x, (v_A + v_M))})), "[x, Matrix([\\n[1, 2],\\n[3, 4]]) + A]")
  expect_equal(.bat_fmt(local({v_M <- Matrix(list(list(Integer(1), Integer(2)), list(Integer(3), Integer(4)))); Add(v_M, v_M)})), "Matrix([[2, 4], [6, 8]])")
  expect_equal(.bat_fmt(local({v_M <- Matrix(list(list(Integer(1), Integer(2)), list(Integer(3), Integer(4)))); (2 * v_M)})), "Matrix([[2, 4], [6, 8]])")
  expect_equal(.bat_fmt(local({v_M <- Matrix(list(list(Integer(1), Integer(2)), list(Integer(3), Integer(4)))); v_A <- MatrixSymbol("A", Integer(2), Integer(2)); Add(v_M, v_A, v_M)})), "Matrix([\\n[2, 4],\\n[6, 8]]) + A")
  expect_equal(.bat_fmt(local({v_M <- Matrix(list(list(Integer(1), Integer(2)), list(Integer(3), Integer(4)))); v_A <- MatrixSymbol("A", Integer(2), Integer(2)); Add(v_M, v_M, v_A)})), "Matrix([\\n[2, 4],\\n[6, 8]]) + A")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Integer(2), Integer(2)); v_M <- Matrix(list(list(Integer(1), Integer(2)), list(Integer(3), Integer(4)))); Add(v_A, v_M, v_M)})), "Matrix([\\n[2, 4],\\n[6, 8]]) + A")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Integer(2), Integer(2)); v_M <- Matrix(list(list(Integer(1), Integer(2)), list(Integer(3), Integer(4)))); (v_A + (2 * v_M))})), "Matrix([\\n[2, 4],\\n[6, 8]]) + A")
  expect_equal(.bat_fmt(local({v_x <- Symbol("x"); v_A <- MatrixSymbol("A", Integer(2), Integer(2)); v_M <- Matrix(list(list(Integer(1), Integer(2)), list(Integer(3), Integer(4)))); list((2 * v_x), (v_A + (2 * v_M)))})), "[2*x, Matrix([\\n[2, 4],\\n[6, 8]]) + A]")
})

test_that("test_matexpr: test_simplify_matrix_expressions matches SymPy", {
  expect_equal(.bat_fmt(local({v_C <- MatrixSymbol("C", Symbol("n", integer = TRUE), Symbol("n", integer = TRUE)); v_D <- MatrixSymbol("D", Symbol("n", integer = TRUE), Symbol("n", integer = TRUE)); sym_type(gcd_terms(((v_C * v_D) + (v_D * v_C))))})), "MatAdd")
  expect_equal(.bat_fmt(local({v_a <- MatAdd(MatMul(Integer(2), MatrixSymbol("C", Symbol("n", integer = TRUE), Symbol("n", integer = TRUE)), MatrixSymbol("D", Symbol("n", integer = TRUE), Symbol("n", integer = TRUE))), MatMul(Integer(4), MatrixSymbol("D", Symbol("n", integer = TRUE), Symbol("n", integer = TRUE)), MatrixSymbol("C", Symbol("n", integer = TRUE), Symbol("n", integer = TRUE)))); sym_type(v_a)})), "MatAdd")
  expect_equal(.bat_fmt(local({v_a <- MatAdd(MatMul(Integer(2), MatrixSymbol("C", Symbol("n", integer = TRUE), Symbol("n", integer = TRUE)), MatrixSymbol("D", Symbol("n", integer = TRUE), Symbol("n", integer = TRUE))), MatMul(Integer(4), MatrixSymbol("D", Symbol("n", integer = TRUE), Symbol("n", integer = TRUE)), MatrixSymbol("C", Symbol("n", integer = TRUE), Symbol("n", integer = TRUE)))); sym_args(v_a)})), "[2*C*D, 4*D*C]")
  expect_equal(.bat_fmt(local({v_C <- MatrixSymbol("C", Symbol("n", integer = TRUE), Symbol("n", integer = TRUE)); v_D <- MatrixSymbol("D", Symbol("n", integer = TRUE), Symbol("n", integer = TRUE)); list(((2 * v_C) * v_D), ((4 * v_D) * v_C))})), "[2*C*D, 4*D*C]")
})

test_that("test_matexpr: test_as_explicit matches SymPy", {
  expect_equal(.bat_fmt(local({v_Z <- MatrixSymbol("Z", Integer(2), Integer(3)); as_explicit(v_Z)})), "Matrix([[Z[0, 0], Z[0, 1], Z[0, 2]], [Z[1, 0], Z[1, 1], Z[1, 2]]])")
  expect_equal(.bat_fmt(local({v_Z <- MatrixSymbol("Z", Integer(2), Integer(3)); Matrix(list(list(pyitem(v_Z, 0, 0), pyitem(v_Z, 0, 1), pyitem(v_Z, 0, 2)), list(pyitem(v_Z, 1, 0), pyitem(v_Z, 1, 1), pyitem(v_Z, 1, 2))))})), "Matrix([[Z[0, 0], Z[0, 1], Z[0, 2]], [Z[1, 0], Z[1, 1], Z[1, 2]]])")
  expect_error(local({v_A <- MatrixSymbol("A", Symbol("n", integer = TRUE), Symbol("m", integer = TRUE)); as_explicit(v_A)}), "Matrix with symbolic shape cannot be represented explicitly.", fixed = TRUE, class = "ValueError")
})

test_that("test_matexpr: test_MatrixSet matches SymPy", {
  expect_equal(.bat_fmt(local({v_M <- MatrixSet(Integer(2), Integer(2), Reals); get_shape(v_M)})), "[2, 2]")
  expect_equal(.bat_fmt(local({Reals})), "Reals")
  expect_equal(.bat_fmt(local({v_X <- Matrix(list(list(Integer(1), Integer(2)), list(Integer(3), Integer(4)))); v_X})), "Matrix([[1, 2], [3, 4]])")
  expect_equal(.bat_fmt(local({v_M <- MatrixSet(Integer(2), Integer(2), Reals); v_M})), "MatrixSet(2, 2, Reals)")
  expect_equal(.bat_fmt(local({v_X <- ZeroMatrix(Integer(2), Integer(2)); v_X})), "0")
  expect_equal(.bat_fmt(local({v_M <- MatrixSet(Symbol("n", integer = TRUE), Symbol("m", integer = TRUE), Reals); v_M})), "MatrixSet(n, m, Reals)")
  expect_equal(.bat_fmt(local({S_false})), "False")
  expect_error(local({MatrixSet(2, (-2), Reals)}), "The dimension specification -2 should be a nonnegative integer.", fixed = TRUE, class = "ValueError")
  expect_error(local({MatrixSet(2.4, (-1), Reals)}), "The dimension specification 2.40000000000000 should be a nonnegative integer.", fixed = TRUE, class = "ValueError")
  expect_error(local({MatrixSet(2, 2, list(1, 2, 3))}), "(1, 2, 3) should be an instance of Set.", fixed = TRUE, class = "TypeError")
})

test_that("test_matexpr: test_matrixsymbol_solving matches SymPy", {
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Integer(2), Integer(2)); v_B <- MatrixSymbol("B", Integer(2), Integer(2)); (((-((-v_A) + v_B)) - v_A) + v_B)})), "0")
  expect_equal(.bat_fmt(local({v_Z <- ZeroMatrix(Integer(2), Integer(2)); v_Z})), "0")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Integer(2), Integer(2)); v_B <- MatrixSymbol("B", Integer(2), Integer(2)); simplify((((-((-v_A) + v_B)) - v_A) + v_B))})), "0")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Integer(2), Integer(2)); v_B <- MatrixSymbol("B", Integer(2), Integer(2)); expand((((-((-v_A) + v_B)) - v_A) + v_B))})), "0")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Integer(2), Integer(2)); v_B <- MatrixSymbol("B", Integer(2), Integer(2)); v_Z <- ZeroMatrix(Integer(2), Integer(2)); simplify(((((-((-v_A) + v_B)) - v_A) + v_B) - v_Z))})), "0")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Integer(2), Integer(2)); v_B <- MatrixSymbol("B", Integer(2), Integer(2)); v_Z <- ZeroMatrix(Integer(2), Integer(2)); expand(((((-((-v_A) + v_B)) - v_A) + v_B) - v_Z))})), "0")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Integer(2), Integer(2)); v_B <- MatrixSymbol("B", Integer(2), Integer(2)); ((((v_A ^ 2) + (v_A * v_B)) + (v_B * t(v_A))) + (v_B * t(v_B)))})), "A^2 + A*B + B*A.T + B*B.T")
})

test_that("test_companion: test_creation matches SymPy", {
  expect_error(local({CompanionMatrix(1)}), "1 must be a Poly instance.", fixed = TRUE, class = "ValueError")
})

test_that("test_companion: test_shape matches SymPy", {
  expect_equal(.bat_fmt(local({list(1, 1)})), "[1, 1]")
  expect_equal(.bat_fmt(local({list(3, 3)})), "[3, 3]")
})

test_that("test_companion: test_entry matches SymPy", {
  expect_equal(.bat_fmt(local({v_A <- CompanionMatrix(Poly((Symbol("c0") + (Symbol("x"))^(Integer(3)) + (Symbol("c1")*Symbol("x")) + (Symbol("c2")*(Symbol("x"))^(Integer(2)))), Symbol("x"))); pyitem(v_A, 0, 0)})), "0")
  expect_equal(.bat_fmt(local({v_A <- CompanionMatrix(Poly((Symbol("c0") + (Symbol("x"))^(Integer(3)) + (Symbol("c1")*Symbol("x")) + (Symbol("c2")*(Symbol("x"))^(Integer(2)))), Symbol("x"))); pyitem(v_A, 1, 0)})), "1")
  expect_equal(.bat_fmt(local({v_A <- CompanionMatrix(Poly((Symbol("c0") + (Symbol("x"))^(Integer(3)) + (Symbol("c1")*Symbol("x")) + (Symbol("c2")*(Symbol("x"))^(Integer(2)))), Symbol("x"))); pyitem(v_A, 1, 1)})), "0")
  expect_equal(.bat_fmt(local({v_A <- CompanionMatrix(Poly((Symbol("c0") + (Symbol("x"))^(Integer(3)) + (Symbol("c1")*Symbol("x")) + (Symbol("c2")*(Symbol("x"))^(Integer(2)))), Symbol("x"))); pyitem(v_A, 2, 1)})), "1")
  expect_equal(.bat_fmt(local({v_A <- CompanionMatrix(Poly((Symbol("c0") + (Symbol("x"))^(Integer(3)) + (Symbol("c1")*Symbol("x")) + (Symbol("c2")*(Symbol("x"))^(Integer(2)))), Symbol("x"))); pyitem(v_A, 0, 2)})), "-c0")
  expect_equal(.bat_fmt(local({v_c0 <- Symbol("c0"); (-v_c0)})), "-c0")
  expect_equal(.bat_fmt(local({v_A <- CompanionMatrix(Poly((Symbol("c0") + (Symbol("x"))^(Integer(3)) + (Symbol("c1")*Symbol("x")) + (Symbol("c2")*(Symbol("x"))^(Integer(2)))), Symbol("x"))); pyitem(v_A, 1, 2)})), "-c1")
  expect_equal(.bat_fmt(local({v_c1 <- Symbol("c1"); (-v_c1)})), "-c1")
  expect_equal(.bat_fmt(local({v_A <- CompanionMatrix(Poly((Symbol("c0") + (Symbol("x"))^(Integer(3)) + (Symbol("c1")*Symbol("x")) + (Symbol("c2")*(Symbol("x"))^(Integer(2)))), Symbol("x"))); pyitem(v_A, 2, 2)})), "-c2")
  expect_equal(.bat_fmt(local({v_c2 <- Symbol("c2"); (-v_c2)})), "-c2")
})

test_that("test_companion: test_as_explicit matches SymPy", {
  expect_equal(.bat_fmt(local({v_c0 <- Symbol("c0"); Matrix(list((-v_c0)))})), "Matrix([[-c0]])")
  expect_equal(.bat_fmt(local({v_c0 <- Symbol("c0"); v_c1 <- Symbol("c1"); Matrix(list(list(0, (-v_c0)), list(1, (-v_c1))))})), "Matrix([[0, -c0], [1, -c1]])")
  expect_equal(.bat_fmt(local({v_c0 <- Symbol("c0"); v_c1 <- Symbol("c1"); v_c2 <- Symbol("c2"); Matrix(list(list(0, 0, (-v_c0)), list(1, 0, (-v_c1)), list(0, 1, (-v_c2))))})), "Matrix([[0, 0, -c0], [1, 0, -c1], [0, 1, -c2]])")
})

test_that("test_factorizations: test_LU matches SymPy", {
  expect_equal(.bat_fmt(local({v_L <- LofLU(MatrixSymbol("X", Symbol("n"), Symbol("n"))); get_shape(v_L)})), "[n, n]")
  expect_equal(.bat_fmt(local({v_U <- UofLU(MatrixSymbol("X", Symbol("n"), Symbol("n"))); get_shape(v_U)})), "[n, n]")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("n"), Symbol("n")); get_shape(v_X)})), "[n, n]")
})

test_that("test_factorizations: test_QR matches SymPy", {
  expect_equal(.bat_fmt(local({v_Q_ <- QofQR(MatrixSymbol("X", Symbol("n"), Symbol("n"))); get_shape(v_Q_)})), "[n, n]")
  expect_equal(.bat_fmt(local({v_R <- RofQR(MatrixSymbol("X", Symbol("n"), Symbol("n"))); get_shape(v_R)})), "[n, n]")
})

test_that("test_factorizations: test_svd matches SymPy", {
  expect_equal(.bat_fmt(local({v_U <- UofSVD(MatrixSymbol("X", Symbol("n"), Symbol("n"))); get_shape(v_U)})), "[n, n]")
  expect_equal(.bat_fmt(local({v_S <- SofSVD(MatrixSymbol("X", Symbol("n"), Symbol("n"))); get_shape(v_S)})), "[n, n]")
  expect_equal(.bat_fmt(local({v_V <- VofSVD(MatrixSymbol("X", Symbol("n"), Symbol("n"))); get_shape(v_V)})), "[n, n]")
})

test_that("test_fourier: test_dft_creation matches SymPy", {
  expect_equal(.bat_fmt(local({DFT(2)})), "DFT(2)")
  expect_equal(.bat_fmt(local({DFT(0)})), "DFT(0)")
  expect_error(local({DFT((-1))}), "The dimension specification -1 should be a nonnegative integer.", fixed = TRUE, class = "ValueError")
  expect_equal(.bat_fmt(local({v_n <- Symbol("n"); DFT(v_n)})), "DFT(n)")
  expect_error(local({v_n <- Symbol("n", integer = FALSE); DFT(v_n)}), "The dimension specification n should be a nonnegative integer.", fixed = TRUE, class = "ValueError")
  expect_error(local({v_n <- Symbol("n", negative = TRUE); DFT(v_n)}), "The dimension specification n should be a nonnegative integer.", fixed = TRUE, class = "ValueError")
})

test_that("test_fourier: test_dft matches SymPy", {
  expect_equal(.bat_fmt(local({get_shape(DFT(4))})), "[4, 4]")
  expect_equal(.bat_fmt(local({list(4, 4)})), "[4, 4]")
  expect_equal(.bat_fmt(local({Abs(simplify(sym_det(Matrix(DFT(4)))))})), "1")
  expect_equal(.bat_fmt(local({v_n <- Symbol("n"); (DFT(v_n) * IDFT(v_n))})), "I")
  expect_equal(.bat_fmt(local({v_n <- Symbol("n"); Identity(v_n)})), "I")
  expect_equal(.bat_fmt(local({v_n <- Symbol("n"); v_i <- Symbol("i"); v_j <- Symbol("j"); pyitem(DFT(v_n), v_i, v_j)})), "exp(-2*I*pi/n)^(i*j)/sqrt(n)")
  expect_equal(.bat_fmt(local({v_I <- sym_I; v_n <- Symbol("n"); v_i <- Symbol("i"); v_j <- Symbol("j"); ((exp(((((-2) * S_Pi) * v_I) / v_n)) ^ (v_i * v_j)) / sqrt(v_n))})), "exp(-2*I*pi/n)^(i*j)/sqrt(n)")
})

test_that("test_fourier: test_dft2 matches SymPy", {
  expect_equal(.bat_fmt(local({as_explicit(DFT(1))})), "Matrix([[1]])")
  expect_equal(.bat_fmt(local({Matrix(list(list(1)))})), "Matrix([[1]])")
  expect_equal(.bat_fmt(local({as_explicit(DFT(2))})), "Matrix([[sqrt(2)/2, sqrt(2)/2], [sqrt(2)/2, -sqrt(2)/2]])")
  expect_equal(.bat_fmt(local({as_explicit(DFT(4))})), "Matrix([[1/2, 1/2, 1/2, 1/2], [1/2, -I/2, -1/2, I/2], [1/2, -1/2, 1/2, -1/2], [1/2, I/2, -1/2, -I/2]])")
  expect_equal(.bat_fmt(local({v_I <- sym_I; Matrix(list(list(S_Half, S_Half, S_Half, S_Half), list(S_Half, ((-v_I) / 2), Rational((-1), 2), (v_I / 2)), list(S_Half, Rational((-1), 2), S_Half, Rational((-1), 2)), list(S_Half, (v_I / 2), Rational((-1), 2), ((-v_I) / 2))))})), "Matrix([[1/2, 1/2, 1/2, 1/2], [1/2, -I/2, -1/2, I/2], [1/2, -1/2, 1/2, -1/2], [1/2, I/2, -1/2, -I/2]])")
})

test_that("test_hadamard: test_HadamardProduct matches SymPy", {
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n"), Symbol("m")); v_B <- MatrixSymbol("B", Symbol("n"), Symbol("m")); get_shape(HadamardProduct(v_A, v_B, v_A))})), "[n, m]")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n"), Symbol("m")); get_shape(v_A)})), "[n, m]")
  expect_error(local({v_A <- MatrixSymbol("A", Symbol("n"), Symbol("m")); v_n <- Symbol("n"); HadamardProduct(v_A, v_n)}), "Mix of Matrix and Scalar symbols", fixed = TRUE, class = "TypeError")
  expect_error(local({v_A <- MatrixSymbol("A", Symbol("n"), Symbol("m")); HadamardProduct(v_A, 1)}), "Mix of Matrix and Scalar symbols", fixed = TRUE, class = "TypeError")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n"), Symbol("m")); v_B <- MatrixSymbol("B", Symbol("n"), Symbol("m")); pyitem(HadamardProduct(v_A, (2 * v_B), (-v_A)), 1, 1)})), "-2*A[1, 1]^2*B[1, 1]")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n"), Symbol("m")); v_B <- MatrixSymbol("B", Symbol("n"), Symbol("m")); ((((-2) * pyitem(v_A, 1, 1)) * pyitem(v_B, 1, 1)) * pyitem(v_A, 1, 1))})), "-2*A[1, 1]^2*B[1, 1]")
  expect_equal(.bat_fmt(local({v_mix <- MatMul(HadamardProduct(MatMul(MatrixSymbol("Z", Symbol("n"), Symbol("n")), MatrixSymbol("A", Symbol("n"), Symbol("m"))), MatrixSymbol("B", Symbol("n"), Symbol("m"))), MatrixSymbol("C", Symbol("m"), Symbol("k"))); get_shape(v_mix)})), "[n, k]")
  expect_equal(.bat_fmt(local({v_n <- Symbol("n"); v_k <- Symbol("k"); list(v_n, v_k)})), "[n, k]")
})

test_that("test_hadamard: test_HadamardProduct_isnt_commutative matches SymPy", {
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n"), Symbol("m")); v_B <- MatrixSymbol("B", Symbol("n"), Symbol("m")); HadamardProduct(v_A, v_B)})), "HadamardProduct(A, B)")
  expect_equal(.bat_fmt(local({v_B <- MatrixSymbol("B", Symbol("n"), Symbol("m")); v_A <- MatrixSymbol("A", Symbol("n"), Symbol("m")); HadamardProduct(v_B, v_A)})), "HadamardProduct(B, A)")
})

test_that("test_hadamard: test_mixed_indexing matches SymPy", {
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Integer(2), Integer(2)); v_Y <- MatrixSymbol("Y", Integer(2), Integer(2)); v_Z <- MatrixSymbol("Z", Integer(2), Integer(2)); pyitem((v_X * HadamardProduct(v_Y, v_Z)), 0, 0)})), "X[0, 0]*Y[0, 0]*Z[0, 0] + X[0, 1]*Y[1, 0]*Z[1, 0]")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Integer(2), Integer(2)); v_Y <- MatrixSymbol("Y", Integer(2), Integer(2)); v_Z <- MatrixSymbol("Z", Integer(2), Integer(2)); (((pyitem(v_X, 0, 0) * pyitem(v_Y, 0, 0)) * pyitem(v_Z, 0, 0)) + ((pyitem(v_X, 0, 1) * pyitem(v_Y, 1, 0)) * pyitem(v_Z, 1, 0)))})), "X[0, 0]*Y[0, 0]*Z[0, 0] + X[0, 1]*Y[1, 0]*Z[1, 0]")
})

test_that("test_hadamard: test_canonicalize matches SymPy", {
  expect_equal(.bat_fmt(local({v_expr <- HadamardProduct(MatrixSymbol("X", Integer(2), Integer(2))); v_expr})), "HadamardProduct(X)")
  expect_equal(.bat_fmt(local({v_expr2 <- MatrixSymbol("X", Integer(2), Integer(2)); v_expr2})), "X")
  expect_equal(.bat_fmt(local({v_Z <- ZeroMatrix(Integer(2), Integer(2)); v_X <- MatrixSymbol("X", Integer(2), Integer(2)); doit(HadamardProduct(v_Z, v_X))})), "0")
  expect_equal(.bat_fmt(local({v_U <- OneMatrix(Integer(2), Integer(2)); v_X <- MatrixSymbol("X", Integer(2), Integer(2)); doit(HadamardProduct(v_U, v_X, v_X, v_U))})), "HadamardPower(X, 2)")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Integer(2), Integer(2)); HadamardPower(v_X, 2)})), "HadamardPower(X, 2)")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Integer(2), Integer(2)); v_U <- OneMatrix(Integer(2), Integer(2)); v_Y <- MatrixSymbol("Y", Integer(2), Integer(2)); doit(HadamardProduct(v_X, v_U, v_Y))})), "HadamardProduct(X, Y)")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Integer(2), Integer(2)); v_Y <- MatrixSymbol("Y", Integer(2), Integer(2)); HadamardProduct(v_X, v_Y)})), "HadamardProduct(X, Y)")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Integer(2), Integer(2)); v_Z <- ZeroMatrix(Integer(2), Integer(2)); v_U <- OneMatrix(Integer(2), Integer(2)); v_Y <- MatrixSymbol("Y", Integer(2), Integer(2)); doit(HadamardProduct(v_X, v_Z, v_U, v_Y))})), "0")
})

test_that("test_hadamard: test_hadamard matches SymPy", {
  expect_error(local({hadamard_product()}), "Empty Hadamard product is undefined", fixed = TRUE, class = "TypeError")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("m", integer = TRUE), Symbol("n", integer = TRUE)); hadamard_product(v_A)})), "A")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("m", integer = TRUE), Symbol("n", integer = TRUE)); v_A})), "A")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("m", integer = TRUE), Symbol("n", integer = TRUE)); v_B <- MatrixSymbol("B", Symbol("m", integer = TRUE), Symbol("n", integer = TRUE)); hadamard_product(v_A, v_B)})), "HadamardProduct(A, B)")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("m", integer = TRUE), Symbol("n", integer = TRUE)); v_B <- MatrixSymbol("B", Symbol("m", integer = TRUE), Symbol("n", integer = TRUE)); doit(hadamard_product(v_A, v_B))})), "HadamardProduct(A, B)")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("m", integer = TRUE), Symbol("m", integer = TRUE)); v_I <- Identity(Symbol("m", integer = TRUE)); hadamard_product(v_X, v_I)})), "HadamardProduct(I, X)")
  expect_equal(.bat_fmt(local({v_I <- Identity(Symbol("m", integer = TRUE)); v_X <- MatrixSymbol("X", Symbol("m", integer = TRUE), Symbol("m", integer = TRUE)); HadamardProduct(v_I, v_X)})), "HadamardProduct(I, X)")
  expect_equal(.bat_fmt(local({v_expr <- HadamardProduct(MatAdd(ZeroMatrix(Symbol("k"), Integer(1)), OneMatrix(Symbol("k"), Integer(1))), MatrixSymbol("a", Symbol("k"), Integer(1))); doit(v_expr)})), "a")
  expect_error(local({HadamardProduct()}), "HadamardProduct needs at least one argument", fixed = TRUE, class = "ValueError")
})

test_that("test_hadamard: test_hadamard_product_with_explicit_mat matches SymPy", {
  expect_equal(.bat_fmt(local({v_expr <- Matrix(list(list((MatrixElement(MatrixSymbol("A", Integer(3), Integer(3)), Integer(0), Integer(0))*MatrixElement(MatrixSymbol("B", Integer(3), Integer(3)), Integer(0), Integer(0))), (MatrixElement(MatrixSymbol("A", Integer(3), Integer(3)), Integer(0), Integer(1))*MatrixElement(MatrixSymbol("B", Integer(3), Integer(3)), Integer(0), Integer(1))), (MatrixElement(MatrixSymbol("A", Integer(3), Integer(3)), Integer(0), Integer(2))*MatrixElement(MatrixSymbol("B", Integer(3), Integer(3)), Integer(0), Integer(2)))), list((MatrixElement(MatrixSymbol("A", Integer(3), Integer(3)), Integer(1), Integer(0))*MatrixElement(MatrixSymbol("B", Integer(3), Integer(3)), Integer(1), Integer(0))), (MatrixElement(MatrixSymbol("A", Integer(3), Integer(3)), Integer(1), Integer(1))*MatrixElement(MatrixSymbol("B", Integer(3), Integer(3)), Integer(1), Integer(1))), (MatrixElement(MatrixSymbol("A", Integer(3), Integer(3)), Integer(1), Integer(2))*MatrixElement(MatrixSymbol("B", Integer(3), Integer(3)), Integer(1), Integer(2)))), list((MatrixElement(MatrixSymbol("A", Integer(3), Integer(3)), Integer(2), Integer(0))*MatrixElement(MatrixSymbol("B", Integer(3), Integer(3)), Integer(2), Integer(0))), (MatrixElement(MatrixSymbol("A", Integer(3), Integer(3)), Integer(2), Integer(1))*MatrixElement(MatrixSymbol("B", Integer(3), Integer(3)), Integer(2), Integer(1))), (MatrixElement(MatrixSymbol("A", Integer(3), Integer(3)), Integer(2), Integer(2))*MatrixElement(MatrixSymbol("B", Integer(3), Integer(3)), Integer(2), Integer(2)))))); v_expr})), "Matrix([[A[0, 0]*B[0, 0], A[0, 1]*B[0, 1], A[0, 2]*B[0, 2]], [A[1, 0]*B[1, 0], A[1, 1]*B[1, 1], A[1, 2]*B[1, 2]], [A[2, 0]*B[2, 0], A[2, 1]*B[2, 1], A[2, 2]*B[2, 2]]])")
  expect_equal(.bat_fmt(local({v_ret <- Matrix(list(list((MatrixElement(MatrixSymbol("A", Integer(3), Integer(3)), Integer(0), Integer(0))*MatrixElement(MatrixSymbol("B", Integer(3), Integer(3)), Integer(0), Integer(0))), (MatrixElement(MatrixSymbol("A", Integer(3), Integer(3)), Integer(0), Integer(1))*MatrixElement(MatrixSymbol("B", Integer(3), Integer(3)), Integer(0), Integer(1))), (MatrixElement(MatrixSymbol("A", Integer(3), Integer(3)), Integer(0), Integer(2))*MatrixElement(MatrixSymbol("B", Integer(3), Integer(3)), Integer(0), Integer(2)))), list((MatrixElement(MatrixSymbol("A", Integer(3), Integer(3)), Integer(1), Integer(0))*MatrixElement(MatrixSymbol("B", Integer(3), Integer(3)), Integer(1), Integer(0))), (MatrixElement(MatrixSymbol("A", Integer(3), Integer(3)), Integer(1), Integer(1))*MatrixElement(MatrixSymbol("B", Integer(3), Integer(3)), Integer(1), Integer(1))), (MatrixElement(MatrixSymbol("A", Integer(3), Integer(3)), Integer(1), Integer(2))*MatrixElement(MatrixSymbol("B", Integer(3), Integer(3)), Integer(1), Integer(2)))), list((MatrixElement(MatrixSymbol("A", Integer(3), Integer(3)), Integer(2), Integer(0))*MatrixElement(MatrixSymbol("B", Integer(3), Integer(3)), Integer(2), Integer(0))), (MatrixElement(MatrixSymbol("A", Integer(3), Integer(3)), Integer(2), Integer(1))*MatrixElement(MatrixSymbol("B", Integer(3), Integer(3)), Integer(2), Integer(1))), (MatrixElement(MatrixSymbol("A", Integer(3), Integer(3)), Integer(2), Integer(2))*MatrixElement(MatrixSymbol("B", Integer(3), Integer(3)), Integer(2), Integer(2)))))); v_ret})), "Matrix([[A[0, 0]*B[0, 0], A[0, 1]*B[0, 1], A[0, 2]*B[0, 2]], [A[1, 0]*B[1, 0], A[1, 1]*B[1, 1], A[1, 2]*B[1, 2]], [A[2, 0]*B[2, 0], A[2, 1]*B[2, 1], A[2, 2]*B[2, 2]]])")
  expect_equal(.bat_fmt(local({v_expr <- HadamardProduct(Matrix(list(list((MatrixElement(MatrixSymbol("A", Integer(3), Integer(3)), Integer(0), Integer(0))*MatrixElement(MatrixSymbol("B", Integer(3), Integer(3)), Integer(0), Integer(0))), (MatrixElement(MatrixSymbol("A", Integer(3), Integer(3)), Integer(0), Integer(1))*MatrixElement(MatrixSymbol("B", Integer(3), Integer(3)), Integer(0), Integer(1))), (MatrixElement(MatrixSymbol("A", Integer(3), Integer(3)), Integer(0), Integer(2))*MatrixElement(MatrixSymbol("B", Integer(3), Integer(3)), Integer(0), Integer(2)))), list((MatrixElement(MatrixSymbol("A", Integer(3), Integer(3)), Integer(1), Integer(0))*MatrixElement(MatrixSymbol("B", Integer(3), Integer(3)), Integer(1), Integer(0))), (MatrixElement(MatrixSymbol("A", Integer(3), Integer(3)), Integer(1), Integer(1))*MatrixElement(MatrixSymbol("B", Integer(3), Integer(3)), Integer(1), Integer(1))), (MatrixElement(MatrixSymbol("A", Integer(3), Integer(3)), Integer(1), Integer(2))*MatrixElement(MatrixSymbol("B", Integer(3), Integer(3)), Integer(1), Integer(2)))), list((MatrixElement(MatrixSymbol("A", Integer(3), Integer(3)), Integer(2), Integer(0))*MatrixElement(MatrixSymbol("B", Integer(3), Integer(3)), Integer(2), Integer(0))), (MatrixElement(MatrixSymbol("A", Integer(3), Integer(3)), Integer(2), Integer(1))*MatrixElement(MatrixSymbol("B", Integer(3), Integer(3)), Integer(2), Integer(1))), (MatrixElement(MatrixSymbol("A", Integer(3), Integer(3)), Integer(2), Integer(2))*MatrixElement(MatrixSymbol("B", Integer(3), Integer(3)), Integer(2), Integer(2)))))), MatrixSymbol("X", Integer(3), Integer(3))); v_expr})), "HadamardProduct(Matrix([\\n[A[0, 0]*B[0, 0], A[0, 1]*B[0, 1], A[0, 2]*B[0, 2]],\\n[A[1, 0]*B[1, 0], A[1, 1]*B[1, 1], A[1, 2]*B[1, 2]],\\n[A[2, 0]*B[2, 0], A[2, 1]*B[2, 1], A[2, 2]*B[2, 2]]]), X)")
  expect_equal(.bat_fmt(local({v_ret <- Matrix(list(list((MatrixElement(MatrixSymbol("A", Integer(3), Integer(3)), Integer(0), Integer(0))*MatrixElement(MatrixSymbol("B", Integer(3), Integer(3)), Integer(0), Integer(0))), (MatrixElement(MatrixSymbol("A", Integer(3), Integer(3)), Integer(0), Integer(1))*MatrixElement(MatrixSymbol("B", Integer(3), Integer(3)), Integer(0), Integer(1))), (MatrixElement(MatrixSymbol("A", Integer(3), Integer(3)), Integer(0), Integer(2))*MatrixElement(MatrixSymbol("B", Integer(3), Integer(3)), Integer(0), Integer(2)))), list((MatrixElement(MatrixSymbol("A", Integer(3), Integer(3)), Integer(1), Integer(0))*MatrixElement(MatrixSymbol("B", Integer(3), Integer(3)), Integer(1), Integer(0))), (MatrixElement(MatrixSymbol("A", Integer(3), Integer(3)), Integer(1), Integer(1))*MatrixElement(MatrixSymbol("B", Integer(3), Integer(3)), Integer(1), Integer(1))), (MatrixElement(MatrixSymbol("A", Integer(3), Integer(3)), Integer(1), Integer(2))*MatrixElement(MatrixSymbol("B", Integer(3), Integer(3)), Integer(1), Integer(2)))), list((MatrixElement(MatrixSymbol("A", Integer(3), Integer(3)), Integer(2), Integer(0))*MatrixElement(MatrixSymbol("B", Integer(3), Integer(3)), Integer(2), Integer(0))), (MatrixElement(MatrixSymbol("A", Integer(3), Integer(3)), Integer(2), Integer(1))*MatrixElement(MatrixSymbol("B", Integer(3), Integer(3)), Integer(2), Integer(1))), (MatrixElement(MatrixSymbol("A", Integer(3), Integer(3)), Integer(2), Integer(2))*MatrixElement(MatrixSymbol("B", Integer(3), Integer(3)), Integer(2), Integer(2)))))); v_X <- MatrixSymbol("X", Integer(3), Integer(3)); HadamardProduct(v_ret, v_X)})), "HadamardProduct(Matrix([\\n[A[0, 0]*B[0, 0], A[0, 1]*B[0, 1], A[0, 2]*B[0, 2]],\\n[A[1, 0]*B[1, 0], A[1, 1]*B[1, 1], A[1, 2]*B[1, 2]],\\n[A[2, 0]*B[2, 0], A[2, 1]*B[2, 1], A[2, 2]*B[2, 2]]]), X)")
  expect_equal(.bat_fmt(local({v_expr <- Matrix(list(list(MatrixElement(MatrixSymbol("A", Integer(3), Integer(3)), Integer(0), Integer(0)), Integer(0), Integer(0)), list(Integer(0), MatrixElement(MatrixSymbol("A", Integer(3), Integer(3)), Integer(1), Integer(1)), Integer(0)), list(Integer(0), Integer(0), MatrixElement(MatrixSymbol("A", Integer(3), Integer(3)), Integer(2), Integer(2))))); v_expr})), "Matrix([[A[0, 0], 0, 0], [0, A[1, 1], 0], [0, 0, A[2, 2]]])")
  expect_equal(.bat_fmt(local({v_A <- Matrix(list(list(MatrixElement(MatrixSymbol("A", Integer(3), Integer(3)), Integer(0), Integer(0)), MatrixElement(MatrixSymbol("A", Integer(3), Integer(3)), Integer(0), Integer(1)), MatrixElement(MatrixSymbol("A", Integer(3), Integer(3)), Integer(0), Integer(2))), list(MatrixElement(MatrixSymbol("A", Integer(3), Integer(3)), Integer(1), Integer(0)), MatrixElement(MatrixSymbol("A", Integer(3), Integer(3)), Integer(1), Integer(1)), MatrixElement(MatrixSymbol("A", Integer(3), Integer(3)), Integer(1), Integer(2))), list(MatrixElement(MatrixSymbol("A", Integer(3), Integer(3)), Integer(2), Integer(0)), MatrixElement(MatrixSymbol("A", Integer(3), Integer(3)), Integer(2), Integer(1)), MatrixElement(MatrixSymbol("A", Integer(3), Integer(3)), Integer(2), Integer(2))))); Matrix(list(list(pyitem(v_A, 0, 0), 0, 0), list(0, pyitem(v_A, 1, 1), 0), list(0, 0, pyitem(v_A, 2, 2))))})), "Matrix([[A[0, 0], 0, 0], [0, A[1, 1], 0], [0, 0, A[2, 2]]])")
  expect_equal(.bat_fmt(local({v_expr <- Matrix(list(list(Integer(1), Integer(0), Integer(0)), list(Integer(0), Integer(1), Integer(0)), list(Integer(0), Integer(0), Integer(1)))); v_expr})), "Matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])")
})

test_that("test_hadamard: test_hadamard_power matches SymPy", {
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("m", integer = TRUE), Symbol("n", integer = TRUE)); hadamard_power(v_A, 1)})), "A")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("m", integer = TRUE), Symbol("n", integer = TRUE)); hadamard_power(v_A, 2)})), "HadamardPower(A, 2)")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("m", integer = TRUE), Symbol("n", integer = TRUE)); v_n <- Symbol("n", integer = TRUE); t(hadamard_power(v_A, v_n))})), "HadamardPower(A.T, n)")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("m", integer = TRUE), Symbol("n", integer = TRUE)); v_n <- Symbol("n", integer = TRUE); hadamard_power(t(v_A), v_n)})), "HadamardPower(A.T, n)")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("m", integer = TRUE), Symbol("n", integer = TRUE)); v_n <- Symbol("n", integer = TRUE); pyitem(hadamard_power(v_A, v_n), 0, 0)})), "A[0, 0]^n")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("m", integer = TRUE), Symbol("n", integer = TRUE)); v_n <- Symbol("n", integer = TRUE); (pyitem(v_A, 0, 0) ^ v_n)})), "A[0, 0]^n")
  expect_equal(.bat_fmt(local({v_m <- Symbol("m", integer = TRUE); v_n <- Symbol("n", integer = TRUE); hadamard_power(v_m, v_n)})), "m^n")
  expect_equal(.bat_fmt(local({v_m <- Symbol("m", integer = TRUE); v_n <- Symbol("n", integer = TRUE); (v_m ^ v_n)})), "m^n")
  expect_error(local({v_A <- MatrixSymbol("A", Symbol("m", integer = TRUE), Symbol("n", integer = TRUE)); hadamard_power(v_A, v_A)}), "cannot raise expression to a matrix", fixed = TRUE, class = "ValueError")
})

test_that("test_hadamard: test_hadamard_power_explicit matches SymPy", {
  expect_equal(.bat_fmt(local({v_a <- Symbol("a"); v_b <- Symbol("b"); HadamardPower(v_a, v_b)})), "a^b")
  expect_equal(.bat_fmt(local({v_a <- Symbol("a"); v_b <- Symbol("b"); (v_a ^ v_b)})), "a^b")
  expect_equal(.bat_fmt(local({v_a <- Symbol("a"); v_B <- MatrixSymbol("B", Integer(2), Integer(2)); as_explicit(HadamardPower(v_a, v_B))})), "Matrix([[a^B[0, 0], a^B[0, 1]], [a^B[1, 0], a^B[1, 1]]])")
  expect_equal(.bat_fmt(local({v_a <- Symbol("a"); v_B <- MatrixSymbol("B", Integer(2), Integer(2)); Matrix(list(list((v_a ^ pyitem(v_B, 0, 0)), (v_a ^ pyitem(v_B, 0, 1))), list((v_a ^ pyitem(v_B, 1, 0)), (v_a ^ pyitem(v_B, 1, 1)))))})), "Matrix([[a^B[0, 0], a^B[0, 1]], [a^B[1, 0], a^B[1, 1]]])")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Integer(2), Integer(2)); v_b <- Symbol("b"); as_explicit(HadamardPower(v_A, v_b))})), "Matrix([[A[0, 0]^b, A[0, 1]^b], [A[1, 0]^b, A[1, 1]^b]])")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Integer(2), Integer(2)); v_b <- Symbol("b"); Matrix(list(list((pyitem(v_A, 0, 0) ^ v_b), (pyitem(v_A, 0, 1) ^ v_b)), list((pyitem(v_A, 1, 0) ^ v_b), (pyitem(v_A, 1, 1) ^ v_b))))})), "Matrix([[A[0, 0]^b, A[0, 1]^b], [A[1, 0]^b, A[1, 1]^b]])")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Integer(2), Integer(2)); v_B <- MatrixSymbol("B", Integer(2), Integer(2)); as_explicit(HadamardPower(v_A, v_B))})), "Matrix([[A[0, 0]^B[0, 0], A[0, 1]^B[0, 1]], [A[1, 0]^B[1, 0], A[1, 1]^B[1, 1]]])")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Integer(2), Integer(2)); v_B <- MatrixSymbol("B", Integer(2), Integer(2)); Matrix(list(list((pyitem(v_A, 0, 0) ^ pyitem(v_B, 0, 0)), (pyitem(v_A, 0, 1) ^ pyitem(v_B, 0, 1))), list((pyitem(v_A, 1, 0) ^ pyitem(v_B, 1, 0)), (pyitem(v_A, 1, 1) ^ pyitem(v_B, 1, 1)))))})), "Matrix([[A[0, 0]^B[0, 0], A[0, 1]^B[0, 1]], [A[1, 0]^B[1, 0], A[1, 1]^B[1, 1]]])")
})

test_that("test_hadamard: test_shape_error matches SymPy", {
  expect_error(local({v_A <- MatrixSymbol("A", Integer(2), Integer(3)); v_B <- MatrixSymbol("B", Integer(3), Integer(3)); HadamardProduct(v_A, v_B)}), "Matrices have mismatching shape: (2, 3)", fixed = TRUE, class = "ShapeError")
  expect_error(local({v_A <- MatrixSymbol("A", Integer(2), Integer(3)); v_B <- MatrixSymbol("B", Integer(3), Integer(3)); HadamardPower(v_A, v_B)}), "Matrices have mismatching shape: (2, 3)", fixed = TRUE, class = "ShapeError")
  expect_error(local({v_A <- MatrixSymbol("A", Integer(3), Integer(2)); v_B <- MatrixSymbol("B", Integer(3), Integer(3)); HadamardProduct(v_A, v_B)}), "Matrices have mismatching shape: (2, 3)", fixed = TRUE, class = "ShapeError")
  expect_error(local({v_A <- MatrixSymbol("A", Integer(3), Integer(2)); v_B <- MatrixSymbol("B", Integer(3), Integer(3)); HadamardPower(v_A, v_B)}), "Matrices have mismatching shape: (2, 3)", fixed = TRUE, class = "ShapeError")
})

test_that("test_kronecker: test_KroneckerProduct matches SymPy", {
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n"), Symbol("m")); v_B <- MatrixSymbol("B", Symbol("n"), Symbol("m")); KroneckerProduct(v_A, v_B)})), "KroneckerProduct(A, B)")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n"), Symbol("m")); v_C <- MatrixSymbol("C", Symbol("m"), Symbol("k")); v_B <- MatrixSymbol("B", Symbol("n"), Symbol("m")); subs(KroneckerProduct(v_A, v_B), v_A, v_C)})), "KroneckerProduct(C, B)")
  expect_equal(.bat_fmt(local({v_C <- MatrixSymbol("C", Symbol("m"), Symbol("k")); v_B <- MatrixSymbol("B", Symbol("n"), Symbol("m")); KroneckerProduct(v_C, v_B)})), "KroneckerProduct(C, B)")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n"), Symbol("m")); v_C <- MatrixSymbol("C", Symbol("m"), Symbol("k")); get_shape(KroneckerProduct(v_A, v_C))})), "[m*n, k*m]")
  expect_equal(.bat_fmt(local({v_n <- Symbol("n"); v_m <- Symbol("m"); v_k <- Symbol("k"); list((v_n * v_m), (v_m * v_k))})), "[m*n, k*m]")
})

test_that("test_kronecker: test_KroneckerProduct_identity matches SymPy", {
  expect_equal(.bat_fmt(local({v_m <- Symbol("m"); v_n <- Symbol("n"); KroneckerProduct(Identity(v_m), Identity(v_n))})), "I")
  expect_equal(.bat_fmt(local({v_m <- Symbol("m"); v_n <- Symbol("n"); Identity((v_m * v_n))})), "I")
  expect_equal(.bat_fmt(local({KroneckerProduct(eye(2), eye(3))})), "Matrix([[1, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0], [0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 1]])")
  expect_equal(.bat_fmt(local({eye(6)})), "Matrix([[1, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0], [0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 1]])")
})

test_that("test_kronecker: test_KroneckerProduct_explicit matches SymPy", {
  expect_equal(.bat_fmt(local({v_kp <- KroneckerProduct(MatrixSymbol("X", Integer(2), Integer(2)), MatrixSymbol("Y", Integer(2), Integer(2))); get_shape(v_kp)})), "[4, 4]")
  expect_equal(.bat_fmt(local({v_kp <- KroneckerProduct(MatrixSymbol("X", Integer(2), Integer(2)), MatrixSymbol("Y", Integer(2), Integer(2))); as_explicit(v_kp)})), "Matrix([[X[0, 0]*Y[0, 0], X[0, 0]*Y[0, 1], X[0, 1]*Y[0, 0], X[0, 1]*Y[0, 1]], [X[0, 0]*Y[1, 0], X[0, 0]*Y[1, 1], X[0, 1]*Y[1, 0], X[0, 1]*Y[1, 1]], [X[1, 0]*Y[0, 0], X[1, 0]*Y[0, 1], X[1, 1]*Y[0, 0], X[1, 1]*Y[0, 1]], [X[1, 0]*Y[1, 0], X[1, 0]*Y[1, 1], X[1, 1]*Y[1, 0], X[1, 1]*Y[1, 1]]])")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Integer(2), Integer(2)); v_Y <- MatrixSymbol("Y", Integer(2), Integer(2)); Matrix(list(list((pyitem(v_X, 0, 0) * pyitem(v_Y, 0, 0)), (pyitem(v_X, 0, 0) * pyitem(v_Y, 0, 1)), (pyitem(v_X, 0, 1) * pyitem(v_Y, 0, 0)), (pyitem(v_X, 0, 1) * pyitem(v_Y, 0, 1))), list((pyitem(v_X, 0, 0) * pyitem(v_Y, 1, 0)), (pyitem(v_X, 0, 0) * pyitem(v_Y, 1, 1)), (pyitem(v_X, 0, 1) * pyitem(v_Y, 1, 0)), (pyitem(v_X, 0, 1) * pyitem(v_Y, 1, 1))), list((pyitem(v_X, 1, 0) * pyitem(v_Y, 0, 0)), (pyitem(v_X, 1, 0) * pyitem(v_Y, 0, 1)), (pyitem(v_X, 1, 1) * pyitem(v_Y, 0, 0)), (pyitem(v_X, 1, 1) * pyitem(v_Y, 0, 1))), list((pyitem(v_X, 1, 0) * pyitem(v_Y, 1, 0)), (pyitem(v_X, 1, 0) * pyitem(v_Y, 1, 1)), (pyitem(v_X, 1, 1) * pyitem(v_Y, 1, 0)), (pyitem(v_X, 1, 1) * pyitem(v_Y, 1, 1)))))})), "Matrix([[X[0, 0]*Y[0, 0], X[0, 0]*Y[0, 1], X[0, 1]*Y[0, 0], X[0, 1]*Y[0, 1]], [X[0, 0]*Y[1, 0], X[0, 0]*Y[1, 1], X[0, 1]*Y[1, 0], X[0, 1]*Y[1, 1]], [X[1, 0]*Y[0, 0], X[1, 0]*Y[0, 1], X[1, 1]*Y[0, 0], X[1, 1]*Y[0, 1]], [X[1, 0]*Y[1, 0], X[1, 0]*Y[1, 1], X[1, 1]*Y[1, 0], X[1, 1]*Y[1, 1]]])")
})

test_that("test_kronecker: test_tensor_product_adjoint matches SymPy", {
  expect_equal(.bat_fmt(local({v_I <- sym_I; v_A <- MatrixSymbol("A", Symbol("n"), Symbol("m")); v_B <- MatrixSymbol("B", Symbol("n"), Symbol("m")); adjoint(KroneckerProduct((v_I * v_A), v_B))})), "-I*(KroneckerProduct(Adjoint(A), Adjoint(B)))")
  expect_equal(.bat_fmt(local({v_I <- sym_I; v_A <- MatrixSymbol("A", Symbol("n"), Symbol("m")); v_B <- MatrixSymbol("B", Symbol("n"), Symbol("m")); ((-v_I) * KroneckerProduct(adjoint(v_A), adjoint(v_B)))})), "-I*(KroneckerProduct(Adjoint(A), Adjoint(B)))")
  expect_equal(.bat_fmt(local({v_mat1 <- Matrix(list(list(Integer(1), (Integer(2)*sym_I)), list((Integer(1) + sym_I), Integer(3)))); v_mat2 <- Matrix(list(list((Integer(2)*sym_I), Integer(3)), list((Integer(4)*sym_I), Integer(2)))); adjoint(KroneckerProduct(v_mat1, v_mat2))})), "Matrix([[-2*I, -4*I, -2*I*(1 - I), -4*I*(1 - I)], [3, 2, 3 - 3*I, 2 - 2*I], [-4, -8, -6*I, -12*I], [-6*I, -4*I, 9, 6]])")
  expect_equal(.bat_fmt(local({v_mat1 <- Matrix(list(list(Integer(1), (Integer(2)*sym_I)), list((Integer(1) + sym_I), Integer(3)))); v_mat2 <- Matrix(list(list((Integer(2)*sym_I), Integer(3)), list((Integer(4)*sym_I), Integer(2)))); kronecker_product(adjoint(v_mat1), adjoint(v_mat2))})), "Matrix([[-2*I, -4*I, -2*I*(1 - I), -4*I*(1 - I)], [3, 2, 3 - 3*I, 2 - 2*I], [-4, -8, -6*I, -12*I], [-6*I, -4*I, 9, 6]])")
})

test_that("test_kronecker: test_tensor_product_conjugate matches SymPy", {
  expect_equal(.bat_fmt(local({v_I <- sym_I; v_A <- MatrixSymbol("A", Symbol("n"), Symbol("m")); v_B <- MatrixSymbol("B", Symbol("n"), Symbol("m")); conjugate(KroneckerProduct((v_I * v_A), v_B))})), "-I*(KroneckerProduct(Adjoint(A.T), Adjoint(B.T)))")
  expect_equal(.bat_fmt(local({v_I <- sym_I; v_A <- MatrixSymbol("A", Symbol("n"), Symbol("m")); v_B <- MatrixSymbol("B", Symbol("n"), Symbol("m")); ((-v_I) * KroneckerProduct(conjugate(v_A), conjugate(v_B)))})), "-I*(KroneckerProduct(Adjoint(A.T), Adjoint(B.T)))")
  expect_equal(.bat_fmt(local({v_mat1 <- Matrix(list(list(Integer(1), (Integer(2)*sym_I)), list((Integer(1) + sym_I), Integer(3)))); v_mat2 <- Matrix(list(list((Integer(2)*sym_I), Integer(3)), list((Integer(4)*sym_I), Integer(2)))); conjugate(KroneckerProduct(v_mat1, v_mat2))})), "Matrix([[-2*I, 3, -4, -6*I], [-4*I, 2, -8, -4*I], [-2*I*(1 - I), 3 - 3*I, -6*I, 9], [-4*I*(1 - I), 2 - 2*I, -12*I, 6]])")
  expect_equal(.bat_fmt(local({v_mat1 <- Matrix(list(list(Integer(1), (Integer(2)*sym_I)), list((Integer(1) + sym_I), Integer(3)))); v_mat2 <- Matrix(list(list((Integer(2)*sym_I), Integer(3)), list((Integer(4)*sym_I), Integer(2)))); kronecker_product(conjugate(v_mat1), conjugate(v_mat2))})), "Matrix([[-2*I, 3, -4, -6*I], [-4*I, 2, -8, -4*I], [-2*I*(1 - I), 3 - 3*I, -6*I, 9], [-4*I*(1 - I), 2 - 2*I, -12*I, 6]])")
})

test_that("test_kronecker: test_tensor_product_transpose matches SymPy", {
  expect_equal(.bat_fmt(local({v_I <- sym_I; v_A <- MatrixSymbol("A", Symbol("n"), Symbol("m")); v_B <- MatrixSymbol("B", Symbol("n"), Symbol("m")); t(KroneckerProduct((v_I * v_A), v_B))})), "I*(KroneckerProduct(A.T, B.T))")
  expect_equal(.bat_fmt(local({v_I <- sym_I; v_A <- MatrixSymbol("A", Symbol("n"), Symbol("m")); v_B <- MatrixSymbol("B", Symbol("n"), Symbol("m")); (v_I * KroneckerProduct(t(v_A), t(v_B)))})), "I*(KroneckerProduct(A.T, B.T))")
  expect_equal(.bat_fmt(local({v_mat1 <- Matrix(list(list(Integer(1), (Integer(2)*sym_I)), list((Integer(1) + sym_I), Integer(3)))); v_mat2 <- Matrix(list(list((Integer(2)*sym_I), Integer(3)), list((Integer(4)*sym_I), Integer(2)))); t(KroneckerProduct(v_mat1, v_mat2))})), "Matrix([[2*I, 4*I, 2*I*(1 + I), 4*I*(1 + I)], [3, 2, 3 + 3*I, 2 + 2*I], [-4, -8, 6*I, 12*I], [6*I, 4*I, 9, 6]])")
  expect_equal(.bat_fmt(local({v_mat1 <- Matrix(list(list(Integer(1), (Integer(2)*sym_I)), list((Integer(1) + sym_I), Integer(3)))); v_mat2 <- Matrix(list(list((Integer(2)*sym_I), Integer(3)), list((Integer(4)*sym_I), Integer(2)))); kronecker_product(t(v_mat1), t(v_mat2))})), "Matrix([[2*I, 4*I, 2*I*(1 + I), 4*I*(1 + I)], [3, 2, 3 + 3*I, 2 + 2*I], [-4, -8, 6*I, 12*I], [6*I, 4*I, 9, 6]])")
})

test_that("test_kronecker: test_KroneckerProduct_is_associative matches SymPy", {
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n"), Symbol("m")); v_B <- MatrixSymbol("B", Symbol("n"), Symbol("m")); v_C <- MatrixSymbol("C", Symbol("m"), Symbol("k")); kronecker_product(v_A, kronecker_product(v_B, v_C))})), "KroneckerProduct(A, B, C)")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n"), Symbol("m")); v_B <- MatrixSymbol("B", Symbol("n"), Symbol("m")); v_C <- MatrixSymbol("C", Symbol("m"), Symbol("k")); kronecker_product(kronecker_product(v_A, v_B), v_C)})), "KroneckerProduct(A, B, C)")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n"), Symbol("m")); v_B <- MatrixSymbol("B", Symbol("n"), Symbol("m")); v_C <- MatrixSymbol("C", Symbol("m"), Symbol("k")); KroneckerProduct(v_A, v_B, v_C)})), "KroneckerProduct(A, B, C)")
})

test_that("test_kronecker: test_KroneckerProduct_is_bilinear matches SymPy", {
  expect_equal(.bat_fmt(local({v_x <- Symbol("x"); v_A <- MatrixSymbol("A", Symbol("n"), Symbol("m")); v_B <- MatrixSymbol("B", Symbol("n"), Symbol("m")); kronecker_product((v_x * v_A), v_B)})), "x*(KroneckerProduct(A, B))")
  expect_equal(.bat_fmt(local({v_x <- Symbol("x"); v_A <- MatrixSymbol("A", Symbol("n"), Symbol("m")); v_B <- MatrixSymbol("B", Symbol("n"), Symbol("m")); (v_x * kronecker_product(v_A, v_B))})), "x*(KroneckerProduct(A, B))")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n"), Symbol("m")); v_x <- Symbol("x"); v_B <- MatrixSymbol("B", Symbol("n"), Symbol("m")); kronecker_product(v_A, (v_x * v_B))})), "x*(KroneckerProduct(A, B))")
})

test_that("test_kronecker: test_KroneckerProduct_determinant matches SymPy", {
  expect_equal(.bat_fmt(local({v_kp <- KroneckerProduct(MatrixSymbol("W", Symbol("m"), Symbol("m")), MatrixSymbol("Z", Symbol("n"), Symbol("n"))); sym_det(v_kp)})), "Determinant(W)^n*Determinant(Z)^m")
  expect_equal(.bat_fmt(local({v_W <- MatrixSymbol("W", Symbol("m"), Symbol("m")); v_n <- Symbol("n"); v_Z <- MatrixSymbol("Z", Symbol("n"), Symbol("n")); v_m <- Symbol("m"); ((sym_det(v_W) ^ v_n) * (sym_det(v_Z) ^ v_m))})), "Determinant(W)^n*Determinant(Z)^m")
})

test_that("test_kronecker: test_KroneckerProduct_trace matches SymPy", {
  expect_equal(.bat_fmt(local({v_kp <- KroneckerProduct(MatrixSymbol("W", Symbol("m"), Symbol("m")), MatrixSymbol("Z", Symbol("n"), Symbol("n"))); sym_trace(v_kp)})), "Trace(W)*Trace(Z)")
  expect_equal(.bat_fmt(local({v_W <- MatrixSymbol("W", Symbol("m"), Symbol("m")); v_Z <- MatrixSymbol("Z", Symbol("n"), Symbol("n")); (sym_trace(v_W) * sym_trace(v_Z))})), "Trace(W)*Trace(Z)")
})

test_that("test_kronecker: test_KroneckerProduct_isnt_commutative matches SymPy", {
  expect_equal(.bat_fmt(local({v_B <- MatrixSymbol("B", Symbol("n"), Symbol("m")); v_A <- MatrixSymbol("A", Symbol("n"), Symbol("m")); KroneckerProduct(v_B, v_A)})), "KroneckerProduct(B, A)")
})

test_that("test_kronecker: test_KroneckerProduct_extracts_commutative_part matches SymPy", {
  expect_equal(.bat_fmt(local({v_x <- Symbol("x"); v_A <- MatrixSymbol("A", Symbol("n"), Symbol("m")); v_B <- MatrixSymbol("B", Symbol("n"), Symbol("m")); kronecker_product((v_x * v_A), (2 * v_B))})), "(2*x)*(KroneckerProduct(A, B))")
  expect_equal(.bat_fmt(local({v_x <- Symbol("x"); v_A <- MatrixSymbol("A", Symbol("n"), Symbol("m")); v_B <- MatrixSymbol("B", Symbol("n"), Symbol("m")); ((v_x * 2) * KroneckerProduct(v_A, v_B))})), "(2*x)*(KroneckerProduct(A, B))")
})

test_that("test_kronecker: test_KroneckerProduct_inverse matches SymPy", {
  expect_equal(.bat_fmt(local({v_kp <- KroneckerProduct(MatrixSymbol("W", Symbol("m"), Symbol("m")), MatrixSymbol("Z", Symbol("n"), Symbol("n"))); inverse(v_kp)})), "KroneckerProduct(W^(-1), Z^(-1))")
  expect_equal(.bat_fmt(local({v_W <- MatrixSymbol("W", Symbol("m"), Symbol("m")); v_Z <- MatrixSymbol("Z", Symbol("n"), Symbol("n")); kronecker_product(inverse(v_W), inverse(v_Z))})), "KroneckerProduct(W^(-1), Z^(-1))")
})

test_that("test_kronecker: test_KroneckerProduct_combine_add matches SymPy", {
  expect_equal(.bat_fmt(local({v_kp1 <- KroneckerProduct(MatrixSymbol("A", Symbol("n"), Symbol("m")), MatrixSymbol("B", Symbol("n"), Symbol("m"))); v_kp2 <- KroneckerProduct(MatrixSymbol("C", Symbol("m"), Symbol("k")), MatrixSymbol("W", Symbol("m"), Symbol("m"))); combine_kronecker((v_kp1 * v_kp2))})), "KroneckerProduct(A*C, B*W)")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n"), Symbol("m")); v_C <- MatrixSymbol("C", Symbol("m"), Symbol("k")); v_B <- MatrixSymbol("B", Symbol("n"), Symbol("m")); v_W <- MatrixSymbol("W", Symbol("m"), Symbol("m")); kronecker_product((v_A * v_C), (v_B * v_W))})), "KroneckerProduct(A*C, B*W)")
})

test_that("test_kronecker: test_KroneckerProduct_combine_mul matches SymPy", {
  expect_equal(.bat_fmt(local({v_kp1 <- KroneckerProduct(MatrixSymbol("A", Symbol("n"), Symbol("m")), MatrixSymbol("X", Symbol("m"), Symbol("n"))); v_kp2 <- KroneckerProduct(MatrixSymbol("B", Symbol("n"), Symbol("m")), MatrixSymbol("Y", Symbol("m"), Symbol("n"))); combine_kronecker((v_kp1 + v_kp2))})), "KroneckerProduct(A, X) + KroneckerProduct(B, Y)")
  expect_equal(.bat_fmt(local({v_kp1 <- KroneckerProduct(MatrixSymbol("A", Symbol("n"), Symbol("m")), MatrixSymbol("X", Symbol("m"), Symbol("n"))); v_kp2 <- KroneckerProduct(MatrixSymbol("B", Symbol("n"), Symbol("m")), MatrixSymbol("Y", Symbol("m"), Symbol("n"))); (v_kp1 + v_kp2)})), "KroneckerProduct(A, X) + KroneckerProduct(B, Y)")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n"), Symbol("m")); v_X <- MatrixSymbol("X", Symbol("m"), Symbol("n")); v_B <- MatrixSymbol("B", Symbol("n"), Symbol("m")); combine_kronecker((kronecker_product(v_A, v_X) + kronecker_product(v_B, v_X)))})), "KroneckerProduct(A + B, X)")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n"), Symbol("m")); v_B <- MatrixSymbol("B", Symbol("n"), Symbol("m")); v_X <- MatrixSymbol("X", Symbol("m"), Symbol("n")); kronecker_product((v_A + v_B), v_X)})), "KroneckerProduct(A + B, X)")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n"), Symbol("m")); v_X <- MatrixSymbol("X", Symbol("m"), Symbol("n")); v_Y <- MatrixSymbol("Y", Symbol("m"), Symbol("n")); combine_kronecker((kronecker_product(v_A, v_X) + kronecker_product(v_A, v_Y)))})), "KroneckerProduct(A, X + Y)")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n"), Symbol("m")); v_X <- MatrixSymbol("X", Symbol("m"), Symbol("n")); v_Y <- MatrixSymbol("Y", Symbol("m"), Symbol("n")); kronecker_product(v_A, (v_X + v_Y))})), "KroneckerProduct(A, X + Y)")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n"), Symbol("m")); v_X <- MatrixSymbol("X", Symbol("m"), Symbol("n")); v_B <- MatrixSymbol("B", Symbol("n"), Symbol("m")); combine_kronecker(((2 * kronecker_product(v_A, v_X)) + (3 * kronecker_product(v_B, v_X))))})), "KroneckerProduct(2*A + 3*B, X)")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n"), Symbol("m")); v_B <- MatrixSymbol("B", Symbol("n"), Symbol("m")); v_X <- MatrixSymbol("X", Symbol("m"), Symbol("n")); kronecker_product(((2 * v_A) + (3 * v_B)), v_X)})), "KroneckerProduct(2*A + 3*B, X)")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n"), Symbol("m")); v_X <- MatrixSymbol("X", Symbol("m"), Symbol("n")); combine_kronecker(((2 * kronecker_product(v_A, v_X)) + (3 * kronecker_product(v_A, v_X))))})), "5*(KroneckerProduct(A, X))")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n"), Symbol("m")); v_X <- MatrixSymbol("X", Symbol("m"), Symbol("n")); (5 * kronecker_product(v_A, v_X))})), "5*(KroneckerProduct(A, X))")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n"), Symbol("m")); v_X <- MatrixSymbol("X", Symbol("m"), Symbol("n")); v_B <- MatrixSymbol("B", Symbol("n"), Symbol("m")); v_Y <- MatrixSymbol("Y", Symbol("m"), Symbol("n")); combine_kronecker(((kronecker_product(v_A, v_X) + kronecker_product(v_B, v_X)) + kronecker_product(v_A, v_Y)))})), "KroneckerProduct(A, X + Y) + KroneckerProduct(B, X)")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n"), Symbol("m")); v_B <- MatrixSymbol("B", Symbol("n"), Symbol("m")); v_X <- MatrixSymbol("X", Symbol("m"), Symbol("n")); v_Y <- MatrixSymbol("Y", Symbol("m"), Symbol("n")); list((kronecker_product((v_A + v_B), v_X) + kronecker_product(v_A, v_Y)), (kronecker_product(v_A, (v_X + v_Y)) + kronecker_product(v_B, v_X)))})), "[KroneckerProduct(A + B, X) + KroneckerProduct(A, Y), KroneckerProduct(A, X + Y) + KroneckerProduct(B, X)]")
})

test_that("test_kronecker: test_KroneckerProduct_combine_pow matches SymPy", {
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("n"), Symbol("n")); v_Y <- MatrixSymbol("Y", Symbol("n"), Symbol("n")); v_x <- Symbol("x"); combine_kronecker((KroneckerProduct(v_X, v_Y) ^ v_x))})), "KroneckerProduct(X^x, Y^x)")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("n"), Symbol("n")); v_x <- Symbol("x"); v_Y <- MatrixSymbol("Y", Symbol("n"), Symbol("n")); KroneckerProduct((v_X ^ v_x), (v_Y ^ v_x))})), "KroneckerProduct(X^x, Y^x)")
  expect_equal(.bat_fmt(local({v_x <- Symbol("x"); v_X <- MatrixSymbol("X", Symbol("n"), Symbol("n")); v_Y <- MatrixSymbol("Y", Symbol("n"), Symbol("n")); combine_kronecker((v_x * (KroneckerProduct(v_X, v_Y) ^ 2)))})), "x*(KroneckerProduct(X^2, Y^2))")
  expect_equal(.bat_fmt(local({v_x <- Symbol("x"); v_X <- MatrixSymbol("X", Symbol("n"), Symbol("n")); v_Y <- MatrixSymbol("Y", Symbol("n"), Symbol("n")); (v_x * KroneckerProduct((v_X ^ 2), (v_Y ^ 2)))})), "x*(KroneckerProduct(X^2, Y^2))")
  expect_equal(.bat_fmt(local({v_x <- Symbol("x"); v_X <- MatrixSymbol("X", Symbol("n"), Symbol("n")); v_Y <- MatrixSymbol("Y", Symbol("n"), Symbol("n")); v_A <- MatrixSymbol("A", Symbol("n"), Symbol("m")); v_B <- MatrixSymbol("B", Symbol("n"), Symbol("m")); combine_kronecker(((v_x * (KroneckerProduct(v_X, v_Y) ^ 2)) * KroneckerProduct(v_A, v_B)))})), "x*(KroneckerProduct(X^2*A, Y^2*B))")
  expect_equal(.bat_fmt(local({v_x <- Symbol("x"); v_X <- MatrixSymbol("X", Symbol("n"), Symbol("n")); v_A <- MatrixSymbol("A", Symbol("n"), Symbol("m")); v_Y <- MatrixSymbol("Y", Symbol("n"), Symbol("n")); v_B <- MatrixSymbol("B", Symbol("n"), Symbol("m")); (v_x * KroneckerProduct(((v_X ^ 2) * v_A), ((v_Y ^ 2) * v_B)))})), "x*(KroneckerProduct(X^2*A, Y^2*B))")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n"), Symbol("m")); v_B <- MatrixSymbol("B", Symbol("n"), Symbol("m")); v_m <- Symbol("m"); combine_kronecker((KroneckerProduct(v_A, t(v_B)) ^ v_m))})), "(KroneckerProduct(A, B.T))^m")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n"), Symbol("m")); v_B <- MatrixSymbol("B", Symbol("n"), Symbol("m")); v_m <- Symbol("m"); (KroneckerProduct(v_A, t(v_B)) ^ v_m)})), "(KroneckerProduct(A, B.T))^m")
})

test_that("test_kronecker: test_KroneckerProduct_expand matches SymPy", {
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Symbol("n"), Symbol("n")); v_Y <- MatrixSymbol("Y", Symbol("n"), Symbol("n")); v_Z <- MatrixSymbol("Z", Symbol("n"), Symbol("n")); (((KroneckerProduct(v_X, v_Y) + KroneckerProduct(v_X, v_Z)) + KroneckerProduct(v_Y, v_Y)) + KroneckerProduct(v_Y, v_Z))})), "KroneckerProduct(X, Y) + KroneckerProduct(X, Z) + KroneckerProduct(Y, Y) + KroneckerProduct(Y, Z)")
})

test_that("test_kronecker: test_KroneckerProduct_singletons matches SymPy", {
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Integer(3), Integer(3)); v_A1 <- MatrixSymbol("A1", Integer(1), Integer(1)); v_B <- MatrixSymbol("B", Integer(3), Integer(3)); v_B1 <- MatrixSymbol("B1", Integer(1), Integer(1)); doit(KroneckerProduct(v_A, v_A1, v_B, v_B1))})), "KroneckerProduct(A, B, A1*B1)")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Integer(3), Integer(3)); v_B <- MatrixSymbol("B", Integer(3), Integer(3)); v_A1 <- MatrixSymbol("A1", Integer(1), Integer(1)); v_B1 <- MatrixSymbol("B1", Integer(1), Integer(1)); KroneckerProduct(v_A, v_B, (v_A1 * v_B1))})), "KroneckerProduct(A, B, A1*B1)")
  expect_equal(.bat_fmt(local({v_A1 <- MatrixSymbol("A1", Integer(1), Integer(1)); v_A <- MatrixSymbol("A", Integer(3), Integer(3)); v_B1 <- MatrixSymbol("B1", Integer(1), Integer(1)); doit(KroneckerProduct(v_A1, v_A, v_B1))})), "KroneckerProduct(A, A1*B1)")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Integer(3), Integer(3)); v_A1 <- MatrixSymbol("A1", Integer(1), Integer(1)); v_B1 <- MatrixSymbol("B1", Integer(1), Integer(1)); KroneckerProduct(v_A, (v_A1 * v_B1))})), "KroneckerProduct(A, A1*B1)")
  expect_equal(.bat_fmt(local({v_A1 <- MatrixSymbol("A1", Integer(1), Integer(1)); v_B1 <- MatrixSymbol("B1", Integer(1), Integer(1)); doit(KroneckerProduct(v_A1, v_B1))})), "A1*B1")
  expect_equal(.bat_fmt(local({v_A1 <- MatrixSymbol("A1", Integer(1), Integer(1)); v_B1 <- MatrixSymbol("B1", Integer(1), Integer(1)); (v_A1 * v_B1)})), "A1*B1")
  expect_equal(.bat_fmt(local({v_A1 <- MatrixSymbol("A1", Integer(1), Integer(1)); v_A <- MatrixSymbol("A", Integer(3), Integer(3)); doit(KroneckerProduct(v_A1, v_A))})), "KroneckerProduct(A1, A)")
  expect_equal(.bat_fmt(local({v_A1 <- MatrixSymbol("A1", Integer(1), Integer(1)); v_A <- MatrixSymbol("A", Integer(3), Integer(3)); KroneckerProduct(v_A1, v_A)})), "KroneckerProduct(A1, A)")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Integer(3), Integer(3)); v_A1 <- MatrixSymbol("A1", Integer(1), Integer(1)); doit(KroneckerProduct(Identity(1), v_A, v_A1))})), "KroneckerProduct(A, A1)")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Integer(3), Integer(3)); v_A1 <- MatrixSymbol("A1", Integer(1), Integer(1)); KroneckerProduct(v_A, v_A1)})), "KroneckerProduct(A, A1)")
  expect_equal(.bat_fmt(local({v_A1 <- MatrixSymbol("A1", Integer(1), Integer(1)); v_A <- MatrixSymbol("A", Integer(3), Integer(3)); v_B1 <- MatrixSymbol("B1", Integer(1), Integer(1)); doit(KroneckerProduct((2 * v_A1), v_A, (3 * v_B1)))})), "6*(KroneckerProduct(A, A1*B1))")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Integer(3), Integer(3)); v_A1 <- MatrixSymbol("A1", Integer(1), Integer(1)); v_B1 <- MatrixSymbol("B1", Integer(1), Integer(1)); (6 * KroneckerProduct(v_A, (v_A1 * v_B1)))})), "6*(KroneckerProduct(A, A1*B1))")
  expect_equal(.bat_fmt(local({v_expr <- KroneckerProduct(MatrixSymbol("A", Integer(3), Integer(3)), MatrixSymbol("A1", Integer(1), Integer(1)), MatrixSymbol("B", Integer(3), Integer(3)), MatrixSymbol("B1", Integer(1), Integer(1))); as_explicit(doit(v_expr))})), "Matrix([[A[0, 0]*A1[0, 0]*B[0, 0]*B1[0, 0], A[0, 0]*A1[0, 0]*B[0, 1]*B1[0, 0], A[0, 0]*A1[0, 0]*B[0, 2]*B1[0, 0], A[0, 1]*A1[0, 0]*B[0, 0]*B1[0, 0], A[0, 1]*A1[0, 0]*B[0, 1]*B1[0, 0], A[0, 1]*A1[0, 0]*B[0, 2]*B1[0, 0], A[0, 2]*A1[0, 0]*B[0, 0]*B1[0, 0], A[0, 2]*A1[0, 0]*B[0, 1]*B1[0, 0], A[0, 2]*A1[0, 0]*B[0, 2]*B1[0, 0]], [A[0, 0]*A1[0, 0]*B[1, 0]*B1[0, 0], A[0, 0]*A1[0, 0]*B[1, 1]*B1[0, 0], A[0, 0]*A1[0, 0]*B[1, 2]*B1[0, 0], A[0, 1]*A1[0, 0]*B[1, 0]*B1[0, 0], A[0, 1]*A1[0, 0]*B[1, 1]*B1[0, 0], A[0, 1]*A1[0, 0]*B[1, 2]*B1[0, 0], A[0, 2]*A1[0, 0]*B[1, 0]*B1[0, 0], A[0, 2]*A1[0, 0]*B[1, 1]*B1[0, 0], A[0, 2]*A1[0, 0]*B[1, 2]*B1[0, 0]], [A[0, 0]*A1[0, 0]*B[2, 0]*B1[0, 0], A[0, 0]*A1[0, 0]*B[2, 1]*B1[0, 0], A[0, 0]*A1[0, 0]*B[2, 2]*B1[0, 0], A[0, 1]*A1[0, 0]*B[2, 0]*B1[0, 0], A[0, 1]*A1[0, 0]*B[2, 1]*B1[0, 0], A[0, 1]*A1[0, 0]*B[2, 2]*B1[0, 0], A[0, 2]*A1[0, 0]*B[2, 0]*B1[0, 0], A[0, 2]*A1[0, 0]*B[2, 1]*B1[0, 0], A[0, 2]*A1[0, 0]*B[2, 2]*B1[0, 0]], [A[1, 0]*A1[0, 0]*B[0, 0]*B1[0, 0], A[1, 0]*A1[0, 0]*B[0, 1]*B1[0, 0], A[1, 0]*A1[0, 0]*B[0, 2]*B1[0, 0], A[1, 1]*A1[0, 0]*B[0, 0]*B1[0, 0], A[1, 1]*A1[0, 0]*B[0, 1]*B1[0, 0], A[1, 1]*A1[0, 0]*B[0, 2]*B1[0, 0], A[1, 2]*A1[0, 0]*B[0, 0]*B1[0, 0], A[1, 2]*A1[0, 0]*B[0, 1]*B1[0, 0], A[1, 2]*A1[0, 0]*B[0, 2]*B1[0, 0]], [A[1, 0]*A1[0, 0]*B[1, 0]*B1[0, 0], A[1, 0]*A1[0, 0]*B[1, 1]*B1[0, 0], A[1, 0]*A1[0, 0]*B[1, 2]*B1[0, 0], A[1, 1]*A1[0, 0]*B[1, 0]*B1[0, 0], A[1, 1]*A1[0, 0]*B[1, 1]*B1[0, 0], A[1, 1]*A1[0, 0]*B[1, 2]*B1[0, 0], A[1, 2]*A1[0, 0]*B[1, 0]*B1[0, 0], A[1, 2]*A1[0, 0]*B[1, 1]*B1[0, 0], A[1, 2]*A1[0, 0]*B[1, 2]*B1[0, 0]], [A[1, 0]*A1[0, 0]*B[2, 0]*B1[0, 0], A[1, 0]*A1[0, 0]*B[2, 1]*B1[0, 0], A[1, 0]*A1[0, 0]*B[2, 2]*B1[0, 0], A[1, 1]*A1[0, 0]*B[2, 0]*B1[0, 0], A[1, 1]*A1[0, 0]*B[2, 1]*B1[0, 0], A[1, 1]*A1[0, 0]*B[2, 2]*B1[0, 0], A[1, 2]*A1[0, 0]*B[2, 0]*B1[0, 0], A[1, 2]*A1[0, 0]*B[2, 1]*B1[0, 0], A[1, 2]*A1[0, 0]*B[2, 2]*B1[0, 0]], [A[2, 0]*A1[0, 0]*B[0, 0]*B1[0, 0], A[2, 0]*A1[0, 0]*B[0, 1]*B1[0, 0], A[2, 0]*A1[0, 0]*B[0, 2]*B1[0, 0], A[2, 1]*A1[0, 0]*B[0, 0]*B1[0, 0], A[2, 1]*A1[0, 0]*B[0, 1]*B1[0, 0], A[2, 1]*A1[0, 0]*B[0, 2]*B1[0, 0], A[2, 2]*A1[0, 0]*B[0, 0]*B1[0, 0], A[2, 2]*A1[0, 0]*B[0, 1]*B1[0, 0], A[2, 2]*A1[0, 0]*B[0, 2]*B1[0, 0]], [A[2, 0]*A1[0, 0]*B[1, 0]*B1[0, 0], A[2, 0]*A1[0, 0]*B[1, 1]*B1[0, 0], A[2, 0]*A1[0, 0]*B[1, 2]*B1[0, 0], A[2, 1]*A1[0, 0]*B[1, 0]*B1[0, 0], A[2, 1]*A1[0, 0]*B[1, 1]*B1[0, 0], A[2, 1]*A1[0, 0]*B[1, 2]*B1[0, 0], A[2, 2]*A1[0, 0]*B[1, 0]*B1[0, 0], A[2, 2]*A1[0, 0]*B[1, 1]*B1[0, 0], A[2, 2]*A1[0, 0]*B[1, 2]*B1[0, 0]], [A[2, 0]*A1[0, 0]*B[2, 0]*B1[0, 0], A[2, 0]*A1[0, 0]*B[2, 1]*B1[0, 0], A[2, 0]*A1[0, 0]*B[2, 2]*B1[0, 0], A[2, 1]*A1[0, 0]*B[2, 0]*B1[0, 0], A[2, 1]*A1[0, 0]*B[2, 1]*B1[0, 0], A[2, 1]*A1[0, 0]*B[2, 2]*B1[0, 0], A[2, 2]*A1[0, 0]*B[2, 0]*B1[0, 0], A[2, 2]*A1[0, 0]*B[2, 1]*B1[0, 0], A[2, 2]*A1[0, 0]*B[2, 2]*B1[0, 0]]])")
  expect_equal(.bat_fmt(local({v_expr <- KroneckerProduct(MatrixSymbol("A", Integer(3), Integer(3)), MatrixSymbol("A1", Integer(1), Integer(1)), MatrixSymbol("B", Integer(3), Integer(3)), MatrixSymbol("B1", Integer(1), Integer(1))); as_explicit(v_expr)})), "Matrix([[A[0, 0]*A1[0, 0]*B[0, 0]*B1[0, 0], A[0, 0]*A1[0, 0]*B[0, 1]*B1[0, 0], A[0, 0]*A1[0, 0]*B[0, 2]*B1[0, 0], A[0, 1]*A1[0, 0]*B[0, 0]*B1[0, 0], A[0, 1]*A1[0, 0]*B[0, 1]*B1[0, 0], A[0, 1]*A1[0, 0]*B[0, 2]*B1[0, 0], A[0, 2]*A1[0, 0]*B[0, 0]*B1[0, 0], A[0, 2]*A1[0, 0]*B[0, 1]*B1[0, 0], A[0, 2]*A1[0, 0]*B[0, 2]*B1[0, 0]], [A[0, 0]*A1[0, 0]*B[1, 0]*B1[0, 0], A[0, 0]*A1[0, 0]*B[1, 1]*B1[0, 0], A[0, 0]*A1[0, 0]*B[1, 2]*B1[0, 0], A[0, 1]*A1[0, 0]*B[1, 0]*B1[0, 0], A[0, 1]*A1[0, 0]*B[1, 1]*B1[0, 0], A[0, 1]*A1[0, 0]*B[1, 2]*B1[0, 0], A[0, 2]*A1[0, 0]*B[1, 0]*B1[0, 0], A[0, 2]*A1[0, 0]*B[1, 1]*B1[0, 0], A[0, 2]*A1[0, 0]*B[1, 2]*B1[0, 0]], [A[0, 0]*A1[0, 0]*B[2, 0]*B1[0, 0], A[0, 0]*A1[0, 0]*B[2, 1]*B1[0, 0], A[0, 0]*A1[0, 0]*B[2, 2]*B1[0, 0], A[0, 1]*A1[0, 0]*B[2, 0]*B1[0, 0], A[0, 1]*A1[0, 0]*B[2, 1]*B1[0, 0], A[0, 1]*A1[0, 0]*B[2, 2]*B1[0, 0], A[0, 2]*A1[0, 0]*B[2, 0]*B1[0, 0], A[0, 2]*A1[0, 0]*B[2, 1]*B1[0, 0], A[0, 2]*A1[0, 0]*B[2, 2]*B1[0, 0]], [A[1, 0]*A1[0, 0]*B[0, 0]*B1[0, 0], A[1, 0]*A1[0, 0]*B[0, 1]*B1[0, 0], A[1, 0]*A1[0, 0]*B[0, 2]*B1[0, 0], A[1, 1]*A1[0, 0]*B[0, 0]*B1[0, 0], A[1, 1]*A1[0, 0]*B[0, 1]*B1[0, 0], A[1, 1]*A1[0, 0]*B[0, 2]*B1[0, 0], A[1, 2]*A1[0, 0]*B[0, 0]*B1[0, 0], A[1, 2]*A1[0, 0]*B[0, 1]*B1[0, 0], A[1, 2]*A1[0, 0]*B[0, 2]*B1[0, 0]], [A[1, 0]*A1[0, 0]*B[1, 0]*B1[0, 0], A[1, 0]*A1[0, 0]*B[1, 1]*B1[0, 0], A[1, 0]*A1[0, 0]*B[1, 2]*B1[0, 0], A[1, 1]*A1[0, 0]*B[1, 0]*B1[0, 0], A[1, 1]*A1[0, 0]*B[1, 1]*B1[0, 0], A[1, 1]*A1[0, 0]*B[1, 2]*B1[0, 0], A[1, 2]*A1[0, 0]*B[1, 0]*B1[0, 0], A[1, 2]*A1[0, 0]*B[1, 1]*B1[0, 0], A[1, 2]*A1[0, 0]*B[1, 2]*B1[0, 0]], [A[1, 0]*A1[0, 0]*B[2, 0]*B1[0, 0], A[1, 0]*A1[0, 0]*B[2, 1]*B1[0, 0], A[1, 0]*A1[0, 0]*B[2, 2]*B1[0, 0], A[1, 1]*A1[0, 0]*B[2, 0]*B1[0, 0], A[1, 1]*A1[0, 0]*B[2, 1]*B1[0, 0], A[1, 1]*A1[0, 0]*B[2, 2]*B1[0, 0], A[1, 2]*A1[0, 0]*B[2, 0]*B1[0, 0], A[1, 2]*A1[0, 0]*B[2, 1]*B1[0, 0], A[1, 2]*A1[0, 0]*B[2, 2]*B1[0, 0]], [A[2, 0]*A1[0, 0]*B[0, 0]*B1[0, 0], A[2, 0]*A1[0, 0]*B[0, 1]*B1[0, 0], A[2, 0]*A1[0, 0]*B[0, 2]*B1[0, 0], A[2, 1]*A1[0, 0]*B[0, 0]*B1[0, 0], A[2, 1]*A1[0, 0]*B[0, 1]*B1[0, 0], A[2, 1]*A1[0, 0]*B[0, 2]*B1[0, 0], A[2, 2]*A1[0, 0]*B[0, 0]*B1[0, 0], A[2, 2]*A1[0, 0]*B[0, 1]*B1[0, 0], A[2, 2]*A1[0, 0]*B[0, 2]*B1[0, 0]], [A[2, 0]*A1[0, 0]*B[1, 0]*B1[0, 0], A[2, 0]*A1[0, 0]*B[1, 1]*B1[0, 0], A[2, 0]*A1[0, 0]*B[1, 2]*B1[0, 0], A[2, 1]*A1[0, 0]*B[1, 0]*B1[0, 0], A[2, 1]*A1[0, 0]*B[1, 1]*B1[0, 0], A[2, 1]*A1[0, 0]*B[1, 2]*B1[0, 0], A[2, 2]*A1[0, 0]*B[1, 0]*B1[0, 0], A[2, 2]*A1[0, 0]*B[1, 1]*B1[0, 0], A[2, 2]*A1[0, 0]*B[1, 2]*B1[0, 0]], [A[2, 0]*A1[0, 0]*B[2, 0]*B1[0, 0], A[2, 0]*A1[0, 0]*B[2, 1]*B1[0, 0], A[2, 0]*A1[0, 0]*B[2, 2]*B1[0, 0], A[2, 1]*A1[0, 0]*B[2, 0]*B1[0, 0], A[2, 1]*A1[0, 0]*B[2, 1]*B1[0, 0], A[2, 1]*A1[0, 0]*B[2, 2]*B1[0, 0], A[2, 2]*A1[0, 0]*B[2, 0]*B1[0, 0], A[2, 2]*A1[0, 0]*B[2, 1]*B1[0, 0], A[2, 2]*A1[0, 0]*B[2, 2]*B1[0, 0]]])")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Integer(3), Integer(3)); v_P <- Matrix(list(list(Symbol("p", commutative = FALSE)))); v_B <- MatrixSymbol("B", Integer(3), Integer(3)); v_Q <- Matrix(list(list(Symbol("q", commutative = FALSE)))); doit(KroneckerProduct(v_A, v_P, v_B, v_Q))})), "KroneckerProduct(A, Matrix([[p]]), B, Matrix([[q]]))")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Integer(3), Integer(3)); v_P <- Matrix(list(list(Symbol("p", commutative = FALSE)))); v_B <- MatrixSymbol("B", Integer(3), Integer(3)); v_Q <- Matrix(list(list(Symbol("q", commutative = FALSE)))); KroneckerProduct(v_A, v_P, v_B, v_Q)})), "KroneckerProduct(A, Matrix([[p]]), B, Matrix([[q]]))")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Integer(3), Integer(3)); v_P <- Matrix(list(list(Symbol("p", commutative = FALSE)))); v_B <- MatrixSymbol("B", Integer(3), Integer(3)); v_A1 <- MatrixSymbol("A1", Integer(1), Integer(1)); doit(KroneckerProduct(v_A, v_P, v_B, v_A1))})), "KroneckerProduct(A, Matrix([[p]]), B, A1)")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Integer(3), Integer(3)); v_P <- Matrix(list(list(Symbol("p", commutative = FALSE)))); v_B <- MatrixSymbol("B", Integer(3), Integer(3)); v_A1 <- MatrixSymbol("A1", Integer(1), Integer(1)); KroneckerProduct(v_A, v_P, v_B, v_A1)})), "KroneckerProduct(A, Matrix([[p]]), B, A1)")
  expect_equal(.bat_fmt(local({v_P <- Matrix(list(list(Symbol("p", commutative = FALSE)))); v_A1 <- MatrixSymbol("A1", Integer(1), Integer(1)); v_A <- MatrixSymbol("A", Integer(3), Integer(3)); v_B1 <- MatrixSymbol("B1", Integer(1), Integer(1)); doit(KroneckerProduct(v_P, v_A1, v_A, v_B1))})), "KroneckerProduct(Matrix([[p]]), A, A1*B1)")
  expect_equal(.bat_fmt(local({v_P <- Matrix(list(list(Symbol("p", commutative = FALSE)))); v_A <- MatrixSymbol("A", Integer(3), Integer(3)); v_A1 <- MatrixSymbol("A1", Integer(1), Integer(1)); v_B1 <- MatrixSymbol("B1", Integer(1), Integer(1)); KroneckerProduct(v_P, v_A, (v_A1 * v_B1))})), "KroneckerProduct(Matrix([[p]]), A, A1*B1)")
})

test_that("test_special: test_zero_matrix_creation matches SymPy", {
  expect_error(local({ZeroMatrix((-1), 2)}), "The dimension specification -1 should be a nonnegative integer.", fixed = TRUE, class = "ValueError")
  expect_error(local({ZeroMatrix(2, (-1))}), "The dimension specification -1 should be a nonnegative integer.", fixed = TRUE, class = "ValueError")
  expect_error(local({v_n <- Symbol("n", integer = FALSE); ZeroMatrix(v_n, v_n)}), "The dimension specification n should be a nonnegative integer.", fixed = TRUE, class = "ValueError")
  expect_error(local({v_n <- Symbol("n", negative = TRUE); ZeroMatrix(v_n, v_n)}), "The dimension specification n should be a nonnegative integer.", fixed = TRUE, class = "ValueError")
})

test_that("test_special: test_generic_zero_matrix matches SymPy", {
  expect_equal(.bat_fmt(local({v_z <- GenericZeroMatrix(); v_z})), "0")
  expect_error(local({v_z <- GenericZeroMatrix(); get_shape(v_z)}), "GenericZeroMatrix does not have a specified shape", fixed = TRUE, class = "TypeError")
  expect_error(local({v_z <- GenericZeroMatrix(); get_shape(v_z)[[1]]}), "GenericZeroMatrix does not have a specified shape", fixed = TRUE, class = "TypeError")
  expect_error(local({v_z <- GenericZeroMatrix(); get_shape(v_z)[[2]]}), "GenericZeroMatrix does not have a specified shape", fixed = TRUE, class = "TypeError")
  expect_equal(.bat_fmt(local({MatAdd()})), "0")
  expect_equal(.bat_fmt(local({v_z <- GenericZeroMatrix(); v_A <- MatrixSymbol("A", Symbol("n", integer = TRUE), Symbol("n", integer = TRUE)); MatAdd(v_z, v_A)})), "A")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n", integer = TRUE), Symbol("n", integer = TRUE)); MatAdd(v_A)})), "A")
})

test_that("test_special: test_identity_matrix_creation matches SymPy", {
  expect_equal(.bat_fmt(local({Identity(2)})), "I")
  expect_equal(.bat_fmt(local({Identity(0)})), "I")
  expect_error(local({Identity((-1))}), "The dimension specification -1 should be a nonnegative integer.", fixed = TRUE, class = "ValueError")
  expect_error(local({v_n <- Symbol("n", integer = FALSE); Identity(v_n)}), "The dimension specification n should be a nonnegative integer.", fixed = TRUE, class = "ValueError")
  expect_error(local({v_n <- Symbol("n", negative = TRUE); Identity(v_n)}), "The dimension specification n should be a nonnegative integer.", fixed = TRUE, class = "ValueError")
})

test_that("test_special: test_generic_identity matches SymPy", {
  expect_equal(.bat_fmt(local({v_I <- GenericIdentity(); v_I})), "I")
  expect_equal(.bat_fmt(local({v_I <- GenericIdentity(); (v_I ^ (-1))})), "I")
  expect_error(local({v_I <- GenericIdentity(); get_shape(v_I)}), "GenericIdentity does not have a specified shape", fixed = TRUE, class = "TypeError")
  expect_error(local({v_I <- GenericIdentity(); get_shape(v_I)[[1]]}), "GenericIdentity does not have a specified shape", fixed = TRUE, class = "TypeError")
  expect_error(local({v_I <- GenericIdentity(); get_shape(v_I)[[2]]}), "GenericIdentity does not have a specified shape", fixed = TRUE, class = "TypeError")
  expect_equal(.bat_fmt(local({MatMul()})), "I")
  expect_equal(.bat_fmt(local({v_I <- GenericIdentity(); v_A <- MatrixSymbol("A", Symbol("n", integer = TRUE), Symbol("n", integer = TRUE)); MatMul(v_I, v_A)})), "A")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n", integer = TRUE), Symbol("n", integer = TRUE)); MatMul(v_A)})), "A")
})

test_that("test_special: test_one_matrix_creation matches SymPy", {
  expect_equal(.bat_fmt(local({OneMatrix(2, 2)})), "1")
  expect_equal(.bat_fmt(local({OneMatrix(0, 0)})), "1")
  expect_equal(.bat_fmt(local({Eq(OneMatrix(1, 1), Identity(1))})), "True")
  expect_error(local({OneMatrix((-1), 2)}), "The dimension specification -1 should be a nonnegative integer.", fixed = TRUE, class = "ValueError")
  expect_error(local({OneMatrix(2, (-1))}), "The dimension specification -1 should be a nonnegative integer.", fixed = TRUE, class = "ValueError")
  expect_equal(.bat_fmt(local({v_n <- Symbol("n"); OneMatrix(v_n, v_n)})), "1")
  expect_error(local({v_n <- Symbol("n", integer = FALSE); OneMatrix(v_n, v_n)}), "The dimension specification n should be a nonnegative integer.", fixed = TRUE, class = "ValueError")
  expect_error(local({v_n <- Symbol("n", negative = TRUE); OneMatrix(v_n, v_n)}), "The dimension specification n should be a nonnegative integer.", fixed = TRUE, class = "ValueError")
})

test_that("test_special: test_ZeroMatrix matches SymPy", {
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n", integer = TRUE), Symbol("m", integer = TRUE)); v_Z <- ZeroMatrix(Symbol("n", integer = TRUE), Symbol("m", integer = TRUE)); (v_A + v_Z)})), "A")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n", integer = TRUE), Symbol("m", integer = TRUE)); v_Z <- ZeroMatrix(Symbol("n", integer = TRUE), Symbol("m", integer = TRUE)); (v_A * t(v_Z))})), "0")
  expect_equal(.bat_fmt(local({v_n <- Symbol("n", integer = TRUE); ZeroMatrix(v_n, v_n)})), "0")
  expect_equal(.bat_fmt(local({v_Z <- ZeroMatrix(Symbol("n", integer = TRUE), Symbol("m", integer = TRUE)); v_A <- MatrixSymbol("A", Symbol("n", integer = TRUE), Symbol("m", integer = TRUE)); (v_Z * t(v_A))})), "0")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n", integer = TRUE), Symbol("m", integer = TRUE)); (v_A - v_A)})), "0")
  expect_equal(.bat_fmt(local({v_Z <- ZeroMatrix(Symbol("n", integer = TRUE), Symbol("m", integer = TRUE)); v_Z})), "0")
  expect_equal(.bat_fmt(local({v_Z <- ZeroMatrix(Symbol("n", integer = TRUE), Symbol("m", integer = TRUE)); t(v_Z)})), "0")
  expect_equal(.bat_fmt(local({v_m <- Symbol("m", integer = TRUE); v_n <- Symbol("n", integer = TRUE); ZeroMatrix(v_m, v_n)})), "0")
  expect_equal(.bat_fmt(local({v_Z <- ZeroMatrix(Symbol("n", integer = TRUE), Symbol("m", integer = TRUE)); conjugate(v_Z)})), "0")
  expect_equal(.bat_fmt(local({v_Z <- ZeroMatrix(Symbol("n", integer = TRUE), Symbol("m", integer = TRUE)); adjoint(v_Z)})), "0")
  expect_equal(.bat_fmt(local({v_Z <- ZeroMatrix(Symbol("n", integer = TRUE), Symbol("m", integer = TRUE)); re(v_Z)})), "0")
  expect_equal(.bat_fmt(local({v_Z <- ZeroMatrix(Symbol("n", integer = TRUE), Symbol("m", integer = TRUE)); im(v_Z)})), "0")
  expect_equal(.bat_fmt(local({v_n <- Symbol("n", integer = TRUE); (ZeroMatrix(v_n, v_n) ^ 0)})), "I")
  expect_equal(.bat_fmt(local({as_explicit(ZeroMatrix(3, 3))})), "Matrix([[0, 0, 0], [0, 0, 0], [0, 0, 0]])")
})

test_that("test_special: test_ZeroMatrix_doit matches SymPy", {
  expect_equal(.bat_fmt(local({v_Znn <- ZeroMatrix((Symbol("n", integer = TRUE) + Symbol("n", integer = TRUE)), Symbol("n", integer = TRUE)); doit(v_Znn)})), "0")
  expect_equal(.bat_fmt(local({v_n <- Symbol("n", integer = TRUE); ZeroMatrix((2 * v_n), v_n)})), "0")
  expect_equal(.bat_fmt(local({v_Znn <- ZeroMatrix((Symbol("n", integer = TRUE) + Symbol("n", integer = TRUE)), Symbol("n", integer = TRUE)); get_shape(doit(v_Znn))[[1]]})), "2*n")
})

test_that("test_special: test_OneMatrix matches SymPy", {
  expect_equal(.bat_fmt(local({v_U <- OneMatrix(Symbol("n", integer = TRUE), Symbol("m", integer = TRUE)); get_shape(v_U)})), "[n, m]")
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n", integer = TRUE), Symbol("m", integer = TRUE)); v_U <- OneMatrix(Symbol("n", integer = TRUE), Symbol("m", integer = TRUE)); (v_A + v_U)})), "A + 1")
  expect_equal(.bat_fmt(local({v_U <- OneMatrix(Symbol("n", integer = TRUE), Symbol("m", integer = TRUE)); t(v_U)})), "1")
  expect_equal(.bat_fmt(local({v_m <- Symbol("m", integer = TRUE); v_n <- Symbol("n", integer = TRUE); OneMatrix(v_m, v_n)})), "1")
  expect_equal(.bat_fmt(local({v_U <- OneMatrix(Symbol("n", integer = TRUE), Symbol("m", integer = TRUE)); conjugate(v_U)})), "1")
  expect_equal(.bat_fmt(local({v_U <- OneMatrix(Symbol("n", integer = TRUE), Symbol("m", integer = TRUE)); v_U})), "1")
  expect_equal(.bat_fmt(local({v_U <- OneMatrix(Symbol("n", integer = TRUE), Symbol("m", integer = TRUE)); adjoint(v_U)})), "1")
  expect_equal(.bat_fmt(local({v_U <- OneMatrix(Symbol("n", integer = TRUE), Symbol("m", integer = TRUE)); re(v_U)})), "1")
  expect_equal(.bat_fmt(local({v_U <- OneMatrix(Symbol("n", integer = TRUE), Symbol("m", integer = TRUE)); im(v_U)})), "0")
  expect_equal(.bat_fmt(local({v_n <- Symbol("n", integer = TRUE); (OneMatrix(v_n, v_n) ^ 0)})), "I")
  expect_equal(.bat_fmt(local({v_U <- OneMatrix(Symbol("n", integer = TRUE), Symbol("n", integer = TRUE)); pyitem(v_U, 1, 2)})), "1")
  expect_equal(.bat_fmt(local({v_U <- OneMatrix(Integer(2), Integer(3)); as_explicit(v_U)})), "Matrix([[1, 1, 1], [1, 1, 1]])")
})

test_that("test_special: test_OneMatrix_doit matches SymPy", {
  expect_equal(.bat_fmt(local({v_Unn <- OneMatrix((Symbol("n", integer = TRUE) + Symbol("n", integer = TRUE)), Symbol("n", integer = TRUE)); doit(v_Unn)})), "1")
  expect_equal(.bat_fmt(local({v_n <- Symbol("n", integer = TRUE); OneMatrix((2 * v_n), v_n)})), "1")
  expect_equal(.bat_fmt(local({v_Unn <- OneMatrix((Symbol("n", integer = TRUE) + Symbol("n", integer = TRUE)), Symbol("n", integer = TRUE)); get_shape(doit(v_Unn))[[1]]})), "2*n")
})

test_that("test_special: test_OneMatrix_mul matches SymPy", {
  expect_equal(.bat_fmt(local({v_n <- Symbol("n", integer = TRUE); v_m <- Symbol("m", integer = TRUE); v_k <- Symbol("k", integer = TRUE); (OneMatrix(v_n, v_m) * OneMatrix(v_m, v_k))})), "m*1")
  expect_equal(.bat_fmt(local({v_n <- Symbol("n", integer = TRUE); v_k <- Symbol("k", integer = TRUE); v_m <- Symbol("m", integer = TRUE); (OneMatrix(v_n, v_k) * v_m)})), "m*1")
  expect_equal(.bat_fmt(local({v_w <- MatrixSymbol("w", Symbol("n", integer = TRUE), Integer(1)); (v_w * OneMatrix(1, 1))})), "w")
  expect_equal(.bat_fmt(local({v_w <- MatrixSymbol("w", Symbol("n", integer = TRUE), Integer(1)); v_w})), "w")
  expect_equal(.bat_fmt(local({v_w <- MatrixSymbol("w", Symbol("n", integer = TRUE), Integer(1)); (OneMatrix(1, 1) * t(v_w))})), "w.T")
  expect_equal(.bat_fmt(local({v_w <- MatrixSymbol("w", Symbol("n", integer = TRUE), Integer(1)); t(v_w)})), "w.T")
})

test_that("test_special: test_MatrixUnit matches SymPy", {
  expect_equal(.bat_fmt(local({v_mu <- MatrixUnit(Integer(3), Integer(3), Integer(1), Integer(2)); as_explicit(v_mu)})), "Matrix([[0, 0, 0], [0, 0, 1], [0, 0, 0]])")
  expect_equal(.bat_fmt(local({v_mu <- MatrixUnit(Integer(3), Integer(3), Integer(1), Integer(2)); t(v_mu)})), "MatrixUnit(3, 3, 2, 1)")
  expect_equal(.bat_fmt(local({MatrixUnit(3, 3, 2, 1)})), "MatrixUnit(3, 3, 2, 1)")
  expect_equal(.bat_fmt(local({v_mu <- MatrixUnit(Integer(3), Integer(3), Integer(1), Integer(2)); pyset(free_symbols(v_mu))})), "{}")
  expect_equal(.bat_fmt(local({v_mu <- MatrixUnit(Integer(3), Integer(3), Integer(1), Integer(2)); get_shape(v_mu)})), "[3, 3]")
  expect_equal(.bat_fmt(local({v_mu <- MatrixUnit(Integer(3), Integer(3), Integer(1), Integer(2)); get_shape(v_mu)[[1]]})), "3")
  expect_equal(.bat_fmt(local({v_mu <- MatrixUnit(Integer(3), Integer(3), Integer(1), Integer(2)); get_shape(v_mu)[[2]]})), "3")
  expect_equal(.bat_fmt(local({v_ma <- MatrixUnit(Integer(3), Integer(4), Integer(1), Integer(2)); t(v_ma)})), "MatrixUnit(4, 3, 2, 1)")
  expect_equal(.bat_fmt(local({MatrixUnit(4, 3, 2, 1)})), "MatrixUnit(4, 3, 2, 1)")
  expect_equal(.bat_fmt(local({v_ma <- MatrixUnit(Integer(3), Integer(4), Integer(1), Integer(2)); adjoint(v_ma)})), "MatrixUnit(4, 3, 2, 1)")
  expect_equal(.bat_fmt(local({v_ma <- MatrixUnit(Integer(3), Integer(4), Integer(1), Integer(2)); conjugate(v_ma)})), "MatrixUnit(3, 4, 1, 2)")
  expect_equal(.bat_fmt(local({v_ma <- MatrixUnit(Integer(3), Integer(4), Integer(1), Integer(2)); v_ma})), "MatrixUnit(3, 4, 1, 2)")
  expect_equal(.bat_fmt(local({v_ma <- MatrixUnit(Integer(3), Integer(4), Integer(1), Integer(2)); get_shape(v_ma)})), "[3, 4]")
  expect_equal(.bat_fmt(local({list(3, 4)})), "[3, 4]")
  expect_equal(.bat_fmt(local({v_ma <- MatrixUnit(Integer(3), Integer(4), Integer(1), Integer(2)); get_shape(v_ma)[[1]]})), "3")
  expect_equal(.bat_fmt(local({v_ma <- MatrixUnit(Integer(3), Integer(4), Integer(1), Integer(2)); get_shape(v_ma)[[2]]})), "4")
})

test_that("test_special: test_Identity matches SymPy", {
  expect_equal(.bat_fmt(local({v_A <- MatrixSymbol("A", Symbol("n", integer = TRUE), Symbol("m", integer = TRUE)); v_Im <- Identity(Symbol("m", integer = TRUE)); (v_A * v_Im)})), "A")
  expect_equal(.bat_fmt(local({v_In <- Identity(Symbol("n", integer = TRUE)); v_A <- MatrixSymbol("A", Symbol("n", integer = TRUE), Symbol("m", integer = TRUE)); (v_In * v_A)})), "A")
  expect_equal(.bat_fmt(local({v_In <- Identity(Symbol("n", integer = TRUE)); t(v_In)})), "I")
  expect_equal(.bat_fmt(local({v_In <- Identity(Symbol("n", integer = TRUE)); v_In})), "I")
  expect_equal(.bat_fmt(local({v_In <- Identity(Symbol("n", integer = TRUE)); inverse(v_In)})), "I")
  expect_equal(.bat_fmt(local({v_In <- Identity(Symbol("n", integer = TRUE)); conjugate(v_In)})), "I")
  expect_equal(.bat_fmt(local({v_In <- Identity(Symbol("n", integer = TRUE)); adjoint(v_In)})), "I")
  expect_equal(.bat_fmt(local({v_In <- Identity(Symbol("n", integer = TRUE)); re(v_In)})), "I")
  expect_equal(.bat_fmt(local({v_In <- Identity(Symbol("n", integer = TRUE)); im(v_In)})), "0")
  expect_equal(.bat_fmt(local({v_In <- Identity(Symbol("n", integer = TRUE)); v_i <- Symbol("i"); v_j <- Symbol("j"); pyitem(v_In, v_i, v_j)})), "KroneckerDelta(i, j, (0, n - 1))")
  expect_equal(.bat_fmt(local({v_n <- Symbol("n", integer = TRUE); v_In <- Identity(Symbol("n", integer = TRUE)); v_i <- Symbol("i"); v_j <- Symbol("j"); doit(subs(Sum(pyitem(v_In, v_i, v_j), list(v_i, 0, (v_n - 1)), list(v_j, 0, (v_n - 1))), v_n, 3))})), "3")
  expect_equal(.bat_fmt(local({v_n <- Symbol("n", integer = TRUE); v_In <- Identity(Symbol("n", integer = TRUE)); v_i <- Symbol("i"); v_j <- Symbol("j"); doit(subs(Sum(Sum(pyitem(v_In, v_i, v_j), list(v_i, 0, (v_n - 1))), list(v_j, 0, (v_n - 1))), v_n, 3))})), "3")
  expect_equal(.bat_fmt(local({v_expr <- Sum(KroneckerDelta(Symbol("i"), Symbol("j"), Tuple(Integer(0), (Integer(-1) + Symbol("n", integer = TRUE)))), list(Symbol("i"), Integer(0), (Integer(-1) + Symbol("n", integer = TRUE)))); doit(v_expr)})), "1")
  expect_equal(.bat_fmt(local({v_expr <- Sum(KroneckerDelta(Symbol("i"), Symbol("j"), Tuple(Integer(0), (Integer(-1) + Symbol("n", integer = TRUE)))), list(Symbol("i"), Integer(0), (Integer(-2) + Symbol("n", integer = TRUE)))); doit(v_expr)})), "Piecewise((1, (j >= 0) & (j <= n - 2)), (0, True))")
  expect_equal(.bat_fmt(local({v_expr <- Sum(KroneckerDelta(Symbol("i"), Symbol("j"), Tuple(Integer(0), (Integer(-1) + Symbol("n", integer = TRUE)))), list(Symbol("i"), Integer(1), (Integer(-1) + Symbol("n", integer = TRUE)))); doit(v_expr)})), "Piecewise((1, (j >= 1) & (j <= n - 1)), (0, True))")
  expect_equal(.bat_fmt(local({as_explicit(Identity(3))})), "Matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])")
})

test_that("test_special: test_Identity_doit matches SymPy", {
  expect_equal(.bat_fmt(local({v_Inn <- Identity((Symbol("n", integer = TRUE) + Symbol("n", integer = TRUE))); doit(v_Inn)})), "I")
  expect_equal(.bat_fmt(local({v_n <- Symbol("n", integer = TRUE); Identity((2 * v_n))})), "I")
  expect_equal(.bat_fmt(local({v_Inn <- Identity((Symbol("n", integer = TRUE) + Symbol("n", integer = TRUE))); get_shape(doit(v_Inn))[[1]]})), "2*n")
})

test_that("test_applyfunc: test_applyfunc_matrix matches SymPy", {
  expect_equal(.bat_fmt(local({.d1 <- Dummy("x"); v_expr <- ElementwiseApplyFunction(Lambda(.d1, (.d1)^(Integer(2))), Matrix(list(list(MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(0), Integer(0)), MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(0), Integer(1)), MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(0), Integer(2))), list(MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(1), Integer(0)), MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(1), Integer(1)), MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(1), Integer(2))), list(MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(2), Integer(0)), MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(2), Integer(1)), MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(2), Integer(2)))))); v_expr})), "Lambda(_x, _x^2).(Matrix([\\n[X[0, 0], X[0, 1], X[0, 2]],\\n[X[1, 0], X[1, 1], X[1, 2]],\\n[X[2, 0], X[2, 1], X[2, 2]]]))")
  expect_equal(.bat_fmt(local({.d1 <- Dummy("x"); v_expr <- ElementwiseApplyFunction(Lambda(.d1, (.d1)^(Integer(2))), Matrix(list(list(MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(0), Integer(0)), MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(0), Integer(1)), MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(0), Integer(2))), list(MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(1), Integer(0)), MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(1), Integer(1)), MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(1), Integer(2))), list(MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(2), Integer(0)), MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(2), Integer(1)), MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(2), Integer(2)))))); doit(v_expr)})), "Matrix([[X[0, 0]^2, X[0, 1]^2, X[0, 2]^2], [X[1, 0]^2, X[1, 1]^2, X[1, 2]^2], [X[2, 0]^2, X[2, 1]^2, X[2, 2]^2]])")
  expect_equal(.bat_fmt(local({v_Xd <- Matrix(list(list(MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(0), Integer(0)), MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(0), Integer(1)), MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(0), Integer(2))), list(MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(1), Integer(0)), MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(1), Integer(1)), MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(1), Integer(2))), list(MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(2), Integer(0)), MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(2), Integer(1)), MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(2), Integer(2))))); applyfunc(v_Xd, function(x) (x ^ 2))})), "Matrix([[X[0, 0]^2, X[0, 1]^2, X[0, 2]^2], [X[1, 0]^2, X[1, 1]^2, X[1, 2]^2], [X[2, 0]^2, X[2, 1]^2, X[2, 2]^2]])")
  expect_equal(.bat_fmt(local({.d1 <- Dummy("x"); v_expr <- ElementwiseApplyFunction(Lambda(.d1, (.d1)^(Integer(2))), Matrix(list(list(MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(0), Integer(0)), MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(0), Integer(1)), MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(0), Integer(2))), list(MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(1), Integer(0)), MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(1), Integer(1)), MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(1), Integer(2))), list(MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(2), Integer(0)), MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(2), Integer(1)), MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(2), Integer(2)))))); get_shape(v_expr)})), "[3, 3]")
  expect_equal(.bat_fmt(local({.d1 <- Dummy("x"); v_expr <- ElementwiseApplyFunction(Lambda(.d1, (.d1)^(Integer(2))), Matrix(list(list(MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(0), Integer(0)), MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(0), Integer(1)), MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(0), Integer(2))), list(MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(1), Integer(0)), MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(1), Integer(1)), MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(1), Integer(2))), list(MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(2), Integer(0)), MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(2), Integer(1)), MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(2), Integer(2)))))); simplify(v_expr)})), "Lambda(_x, _x^2).(Matrix([\\n[X[0, 0], X[0, 1], X[0, 2]],\\n[X[1, 0], X[1, 1], X[1, 2]],\\n[X[2, 0], X[2, 1], X[2, 2]]]))")
  expect_equal(.bat_fmt(local({.d1 <- Dummy("x"); v_expr <- ElementwiseApplyFunction(Lambda(.d1, (.d1)^(Integer(2))), Matrix(list(list(MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(0), Integer(0)), MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(0), Integer(1)), MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(0), Integer(2))), list(MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(1), Integer(0)), MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(1), Integer(1)), MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(1), Integer(2))), list(MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(2), Integer(0)), MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(2), Integer(1)), MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(2), Integer(2)))))); pyitem(v_expr, 0, 0)})), "X[0, 0]^2")
  expect_equal(.bat_fmt(local({.d1 <- Dummy("x"); v_expr <- ElementwiseApplyFunction(Lambda(.d1, (.d1)^(Integer(2))), MatrixSymbol("X", Integer(3), Integer(3))); v_expr})), "Lambda(_x, _x^2).(X)")
  expect_equal(.bat_fmt(local({.d1 <- Dummy("x"); v_expr <- ElementwiseApplyFunction(Lambda(.d1, (.d1)^(Integer(2))), MatrixSymbol("X", Integer(3), Integer(3))); doit(v_expr)})), "Lambda(_x, _x^2).(X)")
  expect_equal(.bat_fmt(local({.d1 <- Dummy("x"); v_double <- Lambda(.d1, (.d1)^(Integer(2))); v_X <- MatrixSymbol("X", Integer(3), Integer(3)); applyfunc(v_X, v_double)})), "Lambda(_x, _x^2).(X)")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Integer(3), Integer(3)); v_Y <- MatrixSymbol("Y", Integer(3), Integer(3)); (v_X * v_Y)})), "X*Y")
  expect_equal(.bat_fmt(local({.d1 <- Dummy("x"); v_x <- .d1; Lambda(v_x, exp(v_x))})), "Lambda(_x, exp(_x))")
  expect_equal(.bat_fmt(local({.d1 <- Dummy("d"); v_expr <- ElementwiseApplyFunction(Lambda(.d1, exp(.d1)), MatMul(MatrixSymbol("X", Integer(3), Integer(3)), MatrixSymbol("Y", Integer(3), Integer(3)))); v_expr})), "Lambda(_d, exp(_d)).(X*Y)")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Integer(3), Integer(3)); v_Y <- MatrixSymbol("Y", Integer(3), Integer(3)); applyfunc((v_X * v_Y), exp)})), "Lambda(_d, exp(_d)).(X*Y)")
  expect_equal(.bat_fmt(local({.d1 <- Dummy("d"); v_X <- MatrixSymbol("X", Integer(3), Integer(3)); v_expr <- ElementwiseApplyFunction(Lambda(.d1, exp(.d1)), MatMul(MatrixSymbol("X", Integer(3), Integer(3)), MatrixSymbol("Y", Integer(3), Integer(3)))); (v_X * v_expr)})), "X*Lambda(_d, exp(_d)).(X*Y)")
  expect_equal(.bat_fmt(local({.d1 <- Dummy("d"); v_X <- MatrixSymbol("X", Integer(3), Integer(3)); v_expr <- ElementwiseApplyFunction(Lambda(.d1, exp(.d1)), MatMul(MatrixSymbol("X", Integer(3), Integer(3)), MatrixSymbol("Y", Integer(3), Integer(3)))); get_shape((v_X * v_expr))})), "[3, 3]")
  expect_equal(.bat_fmt(local({.d1 <- Dummy("d"); v_Z <- MatrixSymbol("Z", Integer(2), Integer(3)); v_expr <- ElementwiseApplyFunction(Lambda(.d1, exp(.d1)), MatMul(MatrixSymbol("X", Integer(3), Integer(3)), MatrixSymbol("Y", Integer(3), Integer(3)))); get_shape((v_Z * v_expr))})), "[2, 3]")
  expect_equal(.bat_fmt(local({list(2, 3)})), "[2, 3]")
  expect_equal(.bat_fmt(local({.d1 <- Dummy("d"); .d2 <- Dummy("d"); v_expr <- MatMul(ElementwiseApplyFunction(Lambda(.d1, exp(.d1)), Transpose(MatrixSymbol("Z", Integer(2), Integer(3)))), ElementwiseApplyFunction(Lambda(.d2, exp(.d2)), MatrixSymbol("Z", Integer(2), Integer(3)))); get_shape(v_expr)})), "[3, 3]")
  expect_equal(.bat_fmt(local({.d1 <- Dummy("d"); .d2 <- Dummy("d"); v_expr <- MatMul(ElementwiseApplyFunction(Lambda(.d1, exp(.d1)), MatrixSymbol("Z", Integer(2), Integer(3))), ElementwiseApplyFunction(Lambda(.d2, exp(.d2)), Transpose(MatrixSymbol("Z", Integer(2), Integer(3))))); get_shape(v_expr)})), "[2, 2]")
  expect_equal(.bat_fmt(local({.d1 <- Dummy("d"); .d2 <- Dummy("x"); v_expr <- ElementwiseApplyFunction(Lambda(.d1, sin(.d1)), Matrix(list(list(.d2, Symbol("y")), list(Symbol("z"), Symbol("t"))))); v_expr})), "Lambda(_d, sin(_d)).(Matrix([\\n[_x, y],\\n[ z, t]]))")
  expect_equal(.bat_fmt(local({.d1 <- Dummy("x"); v_x <- .d1; Lambda(v_x, sin(v_x))})), "Lambda(_x, sin(_x))")
  expect_equal(.bat_fmt(local({.d1 <- Dummy("x"); v_M <- Matrix(list(list(.d1, Symbol("y")), list(Symbol("z"), Symbol("t")))); v_M})), "Matrix([[_x, y], [z, t]])")
  expect_equal(.bat_fmt(local({.d1 <- Dummy("d"); .d2 <- Dummy("x"); v_expr <- ElementwiseApplyFunction(Lambda(.d1, sin(.d1)), Matrix(list(list(.d2, Symbol("y")), list(Symbol("z"), Symbol("t"))))); doit(v_expr)})), "Matrix([[sin(_x), sin(y)], [sin(z), sin(t)]])")
  expect_equal(.bat_fmt(local({.d1 <- Dummy("x"); v_M <- Matrix(list(list(.d1, Symbol("y")), list(Symbol("z"), Symbol("t")))); applyfunc(v_M, sin)})), "Matrix([[sin(_x), sin(y)], [sin(z), sin(t)]])")
  expect_equal(.bat_fmt(local({.d1 <- Dummy("x"); v_x <- .d1; v_y <- Symbol("y"); v_z <- Symbol("z"); v_t <- Symbol("t"); Matrix(list(list(sin(v_x), sin(v_y)), list(sin(v_z), sin(v_t))))})), "Matrix([[sin(_x), sin(y)], [sin(z), sin(t)]])")
  expect_equal(.bat_fmt(local({.d1 <- Dummy("x"); v_expr <- ElementwiseApplyFunction(Lambda(.d1, (.d1)^(Integer(2))), MatrixSymbol("X", Symbol("k"), Symbol("k"))); doit(v_expr)})), "Lambda(_x, _x^2).(X)")
  expect_equal(.bat_fmt(local({.d1 <- Dummy("x"); v_expr <- ElementwiseApplyFunction(Lambda(.d1, (.d1)^(Integer(2))), MatrixSymbol("X", Symbol("k"), Symbol("k"))); v_expr})), "Lambda(_x, _x^2).(X)")
  expect_equal(.bat_fmt(local({.d1 <- Dummy("x"); v_k <- Symbol("k"); v_expr <- ElementwiseApplyFunction(Lambda(.d1, (.d1)^(Integer(2))), MatrixSymbol("X", Symbol("k"), Symbol("k"))); get_shape(subs(v_expr, v_k, 2))})), "[2, 2]")
  expect_equal(.bat_fmt(local({.d1 <- Dummy("x"); v_expr <- ElementwiseApplyFunction(Lambda(.d1, (.d1)^(Integer(2))), MatrixSymbol("X", Symbol("k"), Symbol("k"))); get_shape((v_expr * v_expr))})), "[k, k]")
  expect_equal(.bat_fmt(local({v_k <- Symbol("k"); list(v_k, v_k)})), "[k, k]")
  expect_equal(.bat_fmt(local({.d1 <- Dummy("x"); v_expr2 <- MatMul(Transpose(MatrixSymbol("M", Symbol("k"), Symbol("t"))), ElementwiseApplyFunction(Lambda(.d1, (.d1)^(Integer(2))), MatrixSymbol("X", Symbol("k"), Symbol("k"))), MatrixSymbol("M", Symbol("k"), Symbol("t"))); v_expr2})), "M.T*Lambda(_x, _x^2).(X)*M")
  expect_equal(.bat_fmt(local({.d1 <- Dummy("x"); v_expr2 <- MatMul(Transpose(MatrixSymbol("M", Symbol("k"), Symbol("t"))), ElementwiseApplyFunction(Lambda(.d1, (.d1)^(Integer(2))), MatrixSymbol("X", Symbol("k"), Symbol("k"))), MatrixSymbol("M", Symbol("k"), Symbol("t"))); pyitem(sym_args(v_expr2), 1)})), "Lambda(_x, _x^2).(X)")
  expect_equal(.bat_fmt(local({.d1 <- Dummy("x"); v_expr2 <- MatMul(Transpose(MatrixSymbol("M", Symbol("k"), Symbol("t"))), ElementwiseApplyFunction(Lambda(.d1, (.d1)^(Integer(2))), MatrixSymbol("X", Symbol("k"), Symbol("k"))), MatrixSymbol("M", Symbol("k"), Symbol("t"))); get_shape(v_expr2)})), "[t, t]")
  expect_equal(.bat_fmt(local({v_t <- Symbol("t"); list(v_t, v_t)})), "[t, t]")
  expect_equal(.bat_fmt(local({.d1 <- Dummy("x"); v_expr3 <- MatMul(ElementwiseApplyFunction(Lambda(.d1, (.d1)^(Integer(2))), MatrixSymbol("X", Symbol("k"), Symbol("k"))), MatrixSymbol("M", Symbol("k"), Symbol("t"))); get_shape(v_expr3)})), "[k, t]")
  expect_equal(.bat_fmt(local({v_k <- Symbol("k"); v_t <- Symbol("t"); list(v_k, v_t)})), "[k, t]")
  expect_equal(.bat_fmt(local({.d1 <- Dummy("d"); v_expr1 <- ElementwiseApplyFunction(Lambda(.d1, (Integer(1) + .d1)), MatrixSymbol("X", Symbol("k"), Symbol("k"))); v_expr1})), "Lambda(_d, _d + 1).(X)")
  expect_equal(.bat_fmt(local({v_expr2 <- ElementwiseApplyFunction(Lambda(Symbol("x"), Symbol("x")), MatrixSymbol("X", Symbol("k"), Symbol("k"))); v_expr2})), "Lambda(x, x).(X)")
})

test_that("test_applyfunc: test_applyfunc_entry matches SymPy", {
  expect_equal(.bat_fmt(local({.d1 <- Dummy("d"); v_af <- ElementwiseApplyFunction(Lambda(.d1, sin(.d1)), MatrixSymbol("X", Integer(3), Integer(3))); pyitem(v_af, 0, 0)})), "sin(X[0, 0])")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Integer(3), Integer(3)); sin(pyitem(v_X, 0, 0))})), "sin(X[0, 0])")
  expect_equal(.bat_fmt(local({v_af <- Matrix(list(list(sin(MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(0), Integer(0))), sin(MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(0), Integer(1))), sin(MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(0), Integer(2)))), list(sin(MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(1), Integer(0))), sin(MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(1), Integer(1))), sin(MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(1), Integer(2)))), list(sin(MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(2), Integer(0))), sin(MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(2), Integer(1))), sin(MatrixElement(MatrixSymbol("X", Integer(3), Integer(3)), Integer(2), Integer(2)))))); pyitem(v_af, 0, 0)})), "sin(X[0, 0])")
})

test_that("test_applyfunc: test_applyfunc_as_explicit matches SymPy", {
  expect_equal(.bat_fmt(local({.d1 <- Dummy("d"); v_af <- ElementwiseApplyFunction(Lambda(.d1, sin(.d1)), MatrixSymbol("X", Integer(3), Integer(3))); as_explicit(v_af)})), "Matrix([[sin(X[0, 0]), sin(X[0, 1]), sin(X[0, 2])], [sin(X[1, 0]), sin(X[1, 1]), sin(X[1, 2])], [sin(X[2, 0]), sin(X[2, 1]), sin(X[2, 2])]])")
  expect_equal(.bat_fmt(local({v_X <- MatrixSymbol("X", Integer(3), Integer(3)); Matrix(list(list(sin(pyitem(v_X, 0, 0)), sin(pyitem(v_X, 0, 1)), sin(pyitem(v_X, 0, 2))), list(sin(pyitem(v_X, 1, 0)), sin(pyitem(v_X, 1, 1)), sin(pyitem(v_X, 1, 2))), list(sin(pyitem(v_X, 2, 0)), sin(pyitem(v_X, 2, 1)), sin(pyitem(v_X, 2, 2)))))})), "Matrix([[sin(X[0, 0]), sin(X[0, 1]), sin(X[0, 2])], [sin(X[1, 0]), sin(X[1, 1]), sin(X[1, 2])], [sin(X[2, 0]), sin(X[2, 1]), sin(X[2, 2])]])")
})

test_that("test_applyfunc: test_applyfunc_transpose matches SymPy", {
  expect_equal(.bat_fmt(local({.d1 <- Dummy("d"); v_af <- ElementwiseApplyFunction(Lambda(.d1, sin(.d1)), MatrixSymbol("X", Symbol("k"), Symbol("k"))); t(v_af)})), "Lambda(_d, sin(_d)).(X.T)")
  expect_equal(.bat_fmt(local({v_Xk <- MatrixSymbol("X", Symbol("k"), Symbol("k")); applyfunc(t(v_Xk), sin)})), "Lambda(_d, sin(_d)).(X.T)")
})

test_that("test_applyfunc: test_applyfunc_shape_11_matrices matches SymPy", {
  expect_equal(.bat_fmt(local({.d1 <- Dummy("d"); v_expr <- ElementwiseApplyFunction(Lambda(.d1, sin(.d1)), MatrixSymbol("M", Integer(1), Integer(1))); v_expr})), "Lambda(_d, sin(_d)).(M)")
  expect_equal(.bat_fmt(local({v_expr <- MatMul(Integer(2), MatrixSymbol("M", Integer(1), Integer(1))); v_expr})), "2*M")
  expect_equal(.bat_fmt(local({v_M <- MatrixSymbol("M", Integer(1), Integer(1)); (2 * v_M)})), "2*M")
})

