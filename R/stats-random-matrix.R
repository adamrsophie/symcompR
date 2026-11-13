# Random matrix models (sympy/stats/random_matrix.py and
# random_matrix_models.py): the Gaussian and circular ensembles.

.defclass("RandomMatrixPSpace", "PSpace")
.defclass("RandomMatrixEnsembleModel")
.defclass("GaussianEnsembleModel", "RandomMatrixEnsembleModel")
.defclass("GaussianUnitaryEnsembleModel", "GaussianEnsembleModel")
.defclass("GaussianOrthogonalEnsembleModel", "GaussianEnsembleModel")
.defclass("GaussianSymplecticEnsembleModel", "GaussianEnsembleModel")
.defclass("CircularEnsembleModel", "RandomMatrixEnsembleModel")
.defclass("CircularUnitaryEnsembleModel", "CircularEnsembleModel")
.defclass("CircularOrthogonalEnsembleModel", "CircularEnsembleModel")
.defclass("CircularSymplecticEnsembleModel", "CircularEnsembleModel")

.RandomMatrixPSpace <- function(sym, model = NULL) {
  sym <- .symbol_converter(sym)
  new_sym("RandomMatrixPSpace", if (is.null(model)) list(sym) else list(sym, model))
}
.defm("RandomMatrixPSpace", "model", function(ps) if (length(ps$args) > 1L) ps$args[[2L]] else NULL)
.defm("RandomMatrixPSpace", "compute_density", function(ps, expr, ...) {
  if (is.null(expr$extra$rv_matrix))
    .serr("NotImplementedError", "Currently, no algorithm has been implemented to handle general expressions containing multiple random matrices.")
  .callm(ps$args[[2L]], "density", expr)
})

.rm_model <- function(cls, sym, dim) {
  sym <- .symbol_converter(sym)
  dim <- as_sym(dim)
  isint <- .st_fact(dim, "integer")
  if (identical(isint, FALSE))
    .serr("ValueError", paste0("Dimension of the random matrices must be integers, received ", .rv_pystr(dim), " instead."))
  new_sym(cls, list(sym, dim))
}
.ensemble <- function(cls, sym, dim) {
  model <- .rm_model(cls, sym, dim)
  sym <- model$args[[1L]]; dim <- model$args[[2L]]
  .RandomMatrixSymbol(sym, dim, dim, .RandomMatrixPSpace(sym, model))
}

.defm("RandomMatrixEnsembleModel", "density", function(m, expr) Density(expr))
.defm("RandomMatrixEnsembleModel", "call", function(m, expr) .callm(m, "density", expr))

.gaussian_norm <- function(beta, n) {
  j <- Dummy("j", integer = TRUE, positive = TRUE)
  term1 <- doit(Product(gamma(1 + beta * j / 2) / gamma(S_One + beta / as_sym(2)), list(j, 1, n)))
  term1 * (2 / (beta * n))^(beta * n * (n - 1) / 4 + n / 2) * (2 * S_Pi)^(n / 2)
}
.defm("GaussianEnsembleModel", "joint_eigen", function(m, beta) {
  n <- m$args[[2L]]
  Zbn <- .gaussian_norm(beta, n)
  i <- Dummy("i", integer = TRUE, positive = TRUE)
  j <- Dummy("j", integer = TRUE, positive = TRUE)
  k <- Dummy("k", integer = TRUE, positive = TRUE)
  term1 <- exp(-n / 2 * doit(Sum(Indexed("l", k)^2, list(k, 1, n))))
  inner <- doit(Product(Abs(Indexed("l", j) - Indexed("l", i))^beta, list(j, i + 1, n)))
  term2 <- doit(Product(inner, list(i, 1, n - 1)))
  syms <- lapply(seq_len(as.integer(n$value)), function(q) Indexed("l", q))
  Lambda(syms, term1 * term2 / Zbn)
})
.gauss_density <- function(m, expr, coef, Z) {
  # Lambda(H, f)(expr) of SymPy is f with expr in place of H
  n <- m$args[[2L]]
  exp(-n * coef * Trace(expr^2)) / Z
}
.defm("GaussianUnitaryEnsembleModel", "normalization_constant", function(m) { n <- m$args[[2L]]; 2^(n / 2) * S_Pi^(n^2 / 2) })
.defm("GaussianUnitaryEnsembleModel", "density", function(m, expr) .gauss_density(m, expr, Rational(1, 2), .callm(m, "normalization_constant")))
.defm("GaussianUnitaryEnsembleModel", "joint_eigen_distribution", function(m) .callm(m, "joint_eigen", as_sym(2)))
.defm("GaussianUnitaryEnsembleModel", "level_spacing_distribution", function(m) {
  s <- Dummy("s"); Lambda(s, 32 / S_Pi^2 * s^2 * exp(-4 / S_Pi * s^2))
})
.gauss_mat_norm <- function(m, coef) {
  n <- m$args[[2L]]
  H <- MatrixSymbol("_H", n, n)
  Integral(exp(-n * coef * Trace(H^2)), H)
}
.defm("GaussianOrthogonalEnsembleModel", "normalization_constant", function(m) .gauss_mat_norm(m, Rational(1, 4)))
.defm("GaussianOrthogonalEnsembleModel", "density", function(m, expr) .gauss_density(m, expr, Rational(1, 4), .callm(m, "normalization_constant")))
.defm("GaussianOrthogonalEnsembleModel", "joint_eigen_distribution", function(m) .callm(m, "joint_eigen", S_One))
.defm("GaussianOrthogonalEnsembleModel", "level_spacing_distribution", function(m) {
  s <- Dummy("s"); Lambda(s, S_Pi / 2 * s * exp(-S_Pi / 4 * s^2))
})
.defm("GaussianSymplecticEnsembleModel", "normalization_constant", function(m) .gauss_mat_norm(m, S_One))
.defm("GaussianSymplecticEnsembleModel", "density", function(m, expr) .gauss_density(m, expr, S_One, .callm(m, "normalization_constant")))
.defm("GaussianSymplecticEnsembleModel", "joint_eigen_distribution", function(m) .callm(m, "joint_eigen", as_sym(4)))
.defm("GaussianSymplecticEnsembleModel", "level_spacing_distribution", function(m) {
  s <- Dummy("s"); Lambda(s, as_sym(2)^18 / (as_sym(3)^6 * S_Pi^3) * s^4 * exp(-64 / (9 * S_Pi) * s^2))
})

.defm("CircularEnsembleModel", "density", function(m, expr)
  .serr("NotImplementedError", paste0("Support for Haar measure hasn't been implemented yet, therefore the density of ", .str(m), " cannot be computed.")))
.defm("CircularEnsembleModel", "joint_eigen", function(m, beta) {
  n <- m$args[[2L]]
  Zbn <- (2 * S_Pi)^n * (gamma(beta * n / 2 + 1) / gamma(beta / 2 + 1)^n)
  j <- Dummy("j", integer = TRUE); k <- Dummy("k", integer = TRUE)
  syms <- lapply(seq_len(as.integer(n$value)), function(q) Indexed("t", q))
  inner <- doit(Product(Abs(exp(S_I * Indexed("t", k)) - exp(S_I * Indexed("t", j)))^beta, list(j, k + 1, n)))
  f <- doit(Product(inner, list(k, 1, n - 1)))
  Lambda(syms, f / Zbn)
})
.defm("CircularUnitaryEnsembleModel", "joint_eigen_distribution", function(m) .callm(m, "joint_eigen", as_sym(2)))
.defm("CircularOrthogonalEnsembleModel", "joint_eigen_distribution", function(m) .callm(m, "joint_eigen", S_One))
.defm("CircularSymplecticEnsembleModel", "joint_eigen_distribution", function(m) .callm(m, "joint_eigen", as_sym(4)))

#' Random matrix ensembles
#'
#' The Gaussian and circular ensembles of `sympy.stats`.
#' `joint_eigen_distribution()` and `level_spacing_distribution()` give the
#' joint distribution of the eigenvalues and the level spacing density;
#' `JointEigenDistribution()` is the joint distribution of the eigenvalues
#' of a matrix with random entries.
#'
#' @param sym The name of the random matrix.
#' @param dim Its dimension.
#' @param mat A random matrix symbol, or a matrix of random variables.
#' @param n,m Dimensions of a random matrix symbol.
#' @param pspace A probability space, or `NULL`.
#' @return A random matrix symbol, a `Lambda`, or a distribution.
#' @examples
#' G <- GaussianUnitaryEnsemble("U", 2)
#' level_spacing_distribution(G)
#' @name random_matrix
NULL

#' @rdname random_matrix
#' @export
GaussianEnsemble <- function(sym, dim) .ensemble("GaussianEnsembleModel", sym, dim)
#' @rdname random_matrix
#' @export
GaussianUnitaryEnsemble <- function(sym, dim) .ensemble("GaussianUnitaryEnsembleModel", sym, dim)
#' @rdname random_matrix
#' @export
GaussianOrthogonalEnsemble <- function(sym, dim) .ensemble("GaussianOrthogonalEnsembleModel", sym, dim)
#' @rdname random_matrix
#' @export
GaussianSymplecticEnsemble <- function(sym, dim) .ensemble("GaussianSymplecticEnsembleModel", sym, dim)
#' @rdname random_matrix
#' @export
CircularEnsemble <- function(sym, dim) .ensemble("CircularEnsembleModel", sym, dim)
#' @rdname random_matrix
#' @export
CircularUnitaryEnsemble <- function(sym, dim) .ensemble("CircularUnitaryEnsembleModel", sym, dim)
#' @rdname random_matrix
#' @export
CircularOrthogonalEnsemble <- function(sym, dim) .ensemble("CircularOrthogonalEnsembleModel", sym, dim)
#' @rdname random_matrix
#' @export
CircularSymplecticEnsemble <- function(sym, dim) .ensemble("CircularSymplecticEnsembleModel", sym, dim)

.rm_model_of <- function(mat) {
  if (!(is_sym(mat) && !is.null(mat$extra$rv_matrix)))
    .serr("ValueError", paste0(.rv_pystr(mat), " is not of type, RandomMatrixSymbol."))
  mat$extra$pspace$args[[2L]]
}
#' @rdname random_matrix
#' @export
joint_eigen_distribution <- function(mat) .callm(.rm_model_of(mat), "joint_eigen_distribution")
#' @rdname random_matrix
#' @export
level_spacing_distribution <- function(mat) .callm(mat$extra$pspace$args[[2L]], "level_spacing_distribution")

#' @rdname random_matrix
#' @export
JointEigenDistribution <- function(mat) {
  lam <- Dummy("lambda")
  vals <- roots(charpoly(Matrix(mat), lam), lam, multiple = TRUE)
  if (!all(vapply(vals, is_random, NA)))
    .serr("ValueError", "Eigen values do not have any random expression, joint distribution cannot be generated.")
  new_sym("JointDistributionHandmade", vals)
}

#' @rdname random_matrix
#' @export
RandomMatrixSymbol <- function(sym, n, m, pspace = NULL) {
  if (is.null(pspace)) pspace <- .PSpace0()
  .RandomMatrixSymbol(sym, n, m, pspace)
}
