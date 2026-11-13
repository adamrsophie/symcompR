# Matrix distributions (sympy/stats/matrix_distributions.py): random
# matrices with the matrix gamma, Wishart, matrix normal and matrix
# Student t distributions. A random matrix is a MatrixSymbol node that
# carries its probability space.

.defclass("MatrixPSpace", "PSpace")
.defclass("MatrixDistribution", "Distribution")
.defclass("MatrixDomain", "RandomDomain")

# ImmutableMatrix(arg) of SymPy: nested lists are rows, a flat list is a column
.mat_list <- function(a) {
  if (inherits(a, "symmatrix")) return(.as_block(a))
  if (is.list(a) && !is_sym(a)) {
    nested <- vapply(a, function(r) is.list(r) && !is_sym(r), NA)
    if (any(nested)) {
      lens <- vapply(a, function(r) if (is.list(r) && !is_sym(r)) length(r) else 1L, 0L)
      if (!all(nested) || length(unique(lens)) != 1L) .serr("ValueError", "mismatched dimensions")
    }
    return(.as_block(Matrix(a)))
  }
  as_sym(a)
}

.RandomMatrixSymbol <- function(sym, n, m, pspace) {
  sym <- .symbol_converter(sym)
  n <- as_sym(n); m <- as_sym(m)
  new_sym("MatrixSymbol", list(n, m), name = sym$name,
          extra = list(rv_matrix = TRUE, pspace = pspace, symbol = sym, sympy_class = "RandomMatrixSymbol"),
          key = paste0("RMS[", sym$key, ":", n$key, ":", m$key, ";", pspace$key, "]"))
}

# a positive definite test (Sylvester's criterion on M or on M + M^T)
.mat_is_pd <- function(a) {
  if (!(is_sym(a) && a$type == "ImmutableMatrix")) return(NA)
  M <- as_explicit(a)
  if (M$nrow != M$ncol) return(FALSE)
  sym <- all(vapply(seq_len(M$nrow), function(i) all(vapply(seq_len(M$ncol), function(j)
    identical(.mget(M, i, j)$key, .mget(M, j, i)$key), NA)), NA))
  if (!sym) M <- M + t(M)
  .is_psd(M, strict = TRUE)
}
.mat_is_square <- function(a) { s <- .mshape(a); if (identical(s[[1L]]$key, s[[2L]]$key)) TRUE else if (is_number_type(s[[1L]]) && is_number_type(s[[2L]])) FALSE else NA }

.MatrixPSpace <- function(sym, distribution, dim_n, dim_m) {
  sym <- .symbol_converter(sym)
  dim_n <- as_sym(dim_n); dim_m <- as_sym(dim_m)
  if (!(isTRUE(ask_fact(dim_n, "integer")) && isTRUE(ask_fact(dim_m, "integer"))))
    .serr("ValueError", "Dimensions should be integers")
  new_sym("MatrixPSpace", list(sym, distribution, dim_n, dim_m))
}
.defm("MatrixPSpace", "distribution", function(ps) ps$args[[2L]])
.defm("MatrixPSpace", "symbol", function(ps) ps$args[[1L]])
.defm("MatrixPSpace", "value", function(ps) .RandomMatrixSymbol(ps$args[[1L]], ps$args[[3L]], ps$args[[4L]], ps))
.defm("MatrixPSpace", "values", function(ps) list(.callm(ps, "value")))
.defm("MatrixPSpace", "domain", function(ps) new_sym("MatrixDomain", list(ps$args[[1L]], .dist_set(ps$args[[2L]]))))
.defm("MatrixPSpace", "compute_density", function(ps, expr, ...) {
  rms <- Filter(function(z) !is.null(z$extra$rv_matrix), .uniq(.all_subexprs(expr)))
  if (length(rms) > 1L || is.null(expr$extra$rv_matrix))
    .serr("NotImplementedError", "Currently, no algorithm has been implemented to handle general expressions containing multiple matrix distributions.")
  .callm(ps$args[[2L]], "pdf", expr)
})

.matrix_rv <- function(symbol, cls, args) {
  args <- lapply(args, .mat_list)
  dist <- new_sym(cls, args)
  do.call(.findm(dist, "check"), args)
  dim <- .callm(dist, "dimension")
  .callm(.MatrixPSpace(symbol, dist, dim[[1L]], dim[[2L]]), "value")
}

.defm("MatrixDistribution", "check", function(...) invisible(NULL))
.defm("MatrixDistribution", "call", function(d, expr, ...) .callm(d, "pdf", expr))
.mat_pdf_arg <- function(x) {
  if (inherits(x, "symmatrix") || (is.list(x) && !is_sym(x))) return(.mat_list(x))
  if (!(is_sym(x) && .is_matexpr(x))) .serr("ValueError", paste0(.rv_pystr(x), " should be an isinstance of Matrix or MatrixSymbol"))
  x
}
# SymPy evaluates a product of matrix expressions (MatMul(...).doit()), so
# the inverse of an explicit matrix and a scalar times an explicit matrix
# come out explicit
.is_imm <- function(A) is_sym(A) && A$type == "ImmutableMatrix"
.minv <- function(A) if (.is_imm(A)) .as_block(inverse(as_explicit(A))) else Inverse(A)
.mscal <- function(A, c) if (.is_imm(A)) .as_block(.mat_map(as_explicit(A), function(z) z * c)) else c * A
.mprod <- function(A, B) if (.is_imm(A) && .is_imm(B)) .as_block(.matmul(as_explicit(A), as_explicit(B))) else A * B
.msub <- function(x, M) if (.is_imm(x) && .is_imm(M)) .as_block(as_explicit(x) - as_explicit(M)) else if (.is_imm(M)) .mscal(M, -1) + x else x - M
.mt <- function(A) if (.is_imm(A)) .as_block(t(as_explicit(A))) else .mtranspose(A)
.mset <- function(n, m) new_sym("MatrixSet", list(as_sym(n), as_sym(m), Reals))
.sq_set <- function(d) { k <- .mshape(.da(d, "scale_matrix"))[[1L]]; .mset(k, k) }
.loc_set <- function(d) { s <- .mshape(.da(d, "location_matrix")); .mset(s[[1L]], s[[2L]]) }

.defdist("MatrixGammaDistribution", "MatrixDistribution", c("alpha", "beta", "scale_matrix"),
  check = function(alpha, beta, scale_matrix) {
    if (scale_matrix$type != "MatrixSymbol") .value_check(.mat_is_pd(scale_matrix), "The shape matrix must be positive definite.")
    .value_check(.mat_is_square(scale_matrix), "Should be square matrix")
    .value_check(ask_fact(alpha, "positive"), "Shape parameter should be positive.")
    .value_check(ask_fact(beta, "positive"), "Scale parameter should be positive.")
  },
  set = .sq_set,
  dimension = function(d) .mshape(.da(d, "scale_matrix")),
  pdf = function(d, x) {
    alpha <- d$args[[1L]]; beta <- d$args[[2L]]; sm <- d$args[[3L]]
    p <- .mshape(sm)[[1L]]
    x <- .mat_pdf_arg(x)
    sigma_inv_x <- .mprod(.mscal(.minv(sm), -1 / beta), x)
    term1 <- exp(Trace(sigma_inv_x)) / (beta^(p * alpha) * multigamma(alpha, p))
    term1 * Determinant(sm)^(-alpha) * Determinant(x)^(alpha - (p + 1) / 2)
  })

.defdist("WishartDistribution", "MatrixDistribution", c("n", "scale_matrix"),
  check = function(n, scale_matrix) {
    if (scale_matrix$type != "MatrixSymbol") .value_check(.mat_is_pd(scale_matrix), "The shape matrix must be positive definite.")
    .value_check(.mat_is_square(scale_matrix), "Should be square matrix")
    .value_check(ask_fact(n, "positive"), "Shape parameter should be positive.")
  },
  set = .sq_set,
  dimension = function(d) .mshape(.da(d, "scale_matrix")),
  pdf = function(d, x) {
    n <- d$args[[1L]]; sm <- d$args[[2L]]
    p <- .mshape(sm)[[1L]]
    x <- .mat_pdf_arg(x)
    sigma_inv_x <- .mprod(.mscal(.minv(sm), Rational(-1, 2)), x)
    term1 <- exp(Trace(sigma_inv_x)) / (as_sym(2)^(p * n / as_sym(2)) * multigamma(n / as_sym(2), p))
    term1 * Determinant(sm)^(-n / as_sym(2)) * Determinant(x)^((n - p - 1) / 2)
  })

.defdist("MatrixNormalDistribution", "MatrixDistribution", c("location_matrix", "scale_matrix_1", "scale_matrix_2"),
  check = function(location_matrix, scale_matrix_1, scale_matrix_2) {
    if (scale_matrix_1$type != "MatrixSymbol") .value_check(.mat_is_pd(scale_matrix_1), "The shape matrix must be positive definite.")
    if (scale_matrix_2$type != "MatrixSymbol") .value_check(.mat_is_pd(scale_matrix_2), "The shape matrix must be positive definite.")
    .value_check(.mat_is_square(scale_matrix_1), "Scale matrix 1 should be be square matrix")
    .value_check(.mat_is_square(scale_matrix_2), "Scale matrix 2 should be be square matrix")
    s <- .mshape(location_matrix)
    .value_check(Eq(.mshape(scale_matrix_1)[[1L]], s[[1L]]), paste0("Scale matrix 1 should be of shape ", .str(s[[1L]]), " x ", .str(s[[1L]])))
    .value_check(Eq(.mshape(scale_matrix_2)[[1L]], s[[2L]]), paste0("Scale matrix 2 should be of shape ", .str(s[[2L]]), " x ", .str(s[[2L]])))
  },
  set = .loc_set,
  dimension = function(d) .mshape(.da(d, "location_matrix")),
  pdf = function(d, x) {
    M <- d$args[[1L]]; U <- d$args[[2L]]; V <- d$args[[3L]]
    s <- .mshape(M); n <- s[[1L]]; p <- s[[2L]]
    x <- .mat_pdf_arg(x)
    d <- .msub(x, M)
    term1 <- .mprod(.mprod(.mprod(.minv(V), .mt(d)), .minv(U)), d)
    num <- exp(-Trace(term1) / as_sym(2))
    den <- (2 * S_Pi)^(n * p / 2) * Determinant(U)^(p / 2) * Determinant(V)^(n / 2)
    num / den
  })

.defdist("MatrixStudentTDistribution", "MatrixDistribution", c("nu", "location_matrix", "scale_matrix_1", "scale_matrix_2"),
  check = function(nu, location_matrix, scale_matrix_1, scale_matrix_2) {
    ne_false <- function(v) !identical(v, FALSE)
    if (scale_matrix_1$type != "MatrixSymbol") .value_check(ne_false(.mat_is_pd(scale_matrix_1)), "The shape matrix must be positive definite.")
    if (scale_matrix_2$type != "MatrixSymbol") .value_check(ne_false(.mat_is_pd(scale_matrix_2)), "The shape matrix must be positive definite.")
    .value_check(ne_false(.mat_is_square(scale_matrix_1)), "Scale matrix 1 should be be square matrix")
    .value_check(ne_false(.mat_is_square(scale_matrix_2)), "Scale matrix 2 should be be square matrix")
    s <- .mshape(location_matrix)
    .value_check(Eq(.mshape(scale_matrix_1)[[1L]], s[[2L]]), paste0("Scale matrix 1 should be of shape ", .str(s[[2L]]), " x ", .str(s[[2L]])))
    .value_check(Eq(.mshape(scale_matrix_2)[[1L]], s[[1L]]), paste0("Scale matrix 2 should be of shape ", .str(s[[1L]]), " x ", .str(s[[1L]])))
    .value_check(ne_false(ask_fact(nu, "positive")), "Degrees of freedom must be positive")
  },
  set = .loc_set,
  dimension = function(d) .mshape(.da(d, "location_matrix")),
  pdf = function(d, x) {
    x <- .mat_pdf_arg(x)
    nu <- d$args[[1L]]; M <- d$args[[2L]]; Omega <- d$args[[3L]]; Sigma <- d$args[[4L]]
    s <- .mshape(M); n <- s[[1L]]; p <- s[[2L]]
    # n, p = M.shape: Python ints for an explicit matrix, so n*p/2, -n/2 and
    # -p/2 are Python floats there
    half <- if (.is_imm(M)) function(v) Float(.int_val(v) / 2) else function(v) v / 2
    K <- multigamma((nu + n + p - 1) / 2, p) * Determinant(Omega)^half(-n) * Determinant(Sigma)^half(-p) /
      (S_Pi^half(n * p) * multigamma((nu + p - 1) / 2, p))
    I <- .as_block(eye(as.integer(n$value)))
    d <- .msub(x, M)
    K * Determinant(I + .mprod(.mprod(.mprod(.minv(Sigma), d), .minv(Omega)), .mt(d)))^(-(nu + n + p - 1) / 2)
  })

#' Matrix distributions
#'
#' Random matrices of `sympy.stats`: `MatrixGamma()`, `Wishart()`,
#' `MatrixNormal()` and `MatrixStudentT()`. Matrix arguments are explicit
#' matrices, nested lists (rows) or matrix symbols. `rv_density(M)` is the
#' distribution, which is called on a matrix with [apply_lambda()].
#'
#' @param symbol The name of the random matrix.
#' @param alpha,beta Shape and scale parameters.
#' @param n Degrees of freedom of the Wishart distribution.
#' @param nu Degrees of freedom of the matrix Student t distribution.
#' @param scale_matrix,scale_matrix_1,scale_matrix_2 Scale matrices.
#' @param location_matrix The location (mean) matrix.
#' @return A random matrix symbol.
#' @examples
#' M <- MatrixGamma("M", 1, 2, list(list(1, 0), list(0, 1)))
#' rv_attr(rv_attr(rv_attr(M, "pspace"), "distribution"), "set")
#' @name matrix_distributions
NULL

#' @rdname matrix_distributions
#' @export
MatrixGamma <- function(symbol, alpha, beta, scale_matrix) .matrix_rv(symbol, "MatrixGammaDistribution", list(alpha, beta, scale_matrix))
#' @rdname matrix_distributions
#' @export
Wishart <- function(symbol, n, scale_matrix) .matrix_rv(symbol, "WishartDistribution", list(n, scale_matrix))
#' @rdname matrix_distributions
#' @export
MatrixNormal <- function(symbol, location_matrix, scale_matrix_1, scale_matrix_2)
  .matrix_rv(symbol, "MatrixNormalDistribution", list(location_matrix, scale_matrix_1, scale_matrix_2))
#' @rdname matrix_distributions
#' @export
MatrixStudentT <- function(symbol, nu, location_matrix, scale_matrix_1, scale_matrix_2)
  .matrix_rv(symbol, "MatrixStudentTDistribution", list(nu, location_matrix, scale_matrix_1, scale_matrix_2))
