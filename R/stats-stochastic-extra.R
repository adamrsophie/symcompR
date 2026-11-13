# Exported constructors of sympy.stats that the stochastic process module
# relies on: the handmade distributions, sampling_density(), and the matrix
# forms of Expectation, Variance and Covariance
# (sympy/stats/symbolic_multivariate_probability.py).

#' Handmade distributions
#'
#' `ContinuousDistributionHandmade()`, `DiscreteDistributionHandmade()` and
#' `FiniteDistributionHandmade()` of `sympy.stats`: distributions given
#' directly by a density. The density of a continuous or discrete one is a
#' `Lambda`; a finite one takes a `Dict` (or a named list) from values to
#' probabilities. Read `"pdf"`, `"set"`, `"dict"` or `"pmf"` with
#' [rv_attr()].
#'
#' @param pdf A `Lambda` giving the density.
#' @param set The support; by default the real line, or the integers for a
#'   discrete distribution.
#' @param density A `Dict` or a named list of probabilities.
#' @return A distribution.
#' @examples
#' x <- Symbol("x")
#' CH <- ContinuousDistributionHandmade(Lambda(x, 2 * x), Interval(0, 1))
#' rv_attr(CH, "pdf", x)
#' DH <- DiscreteDistributionHandmade(Lambda(x, 2^(-x)), Range(1, oo))
#' rv_attr(DH, "pdf", 3)
#' FH <- FiniteDistributionHandmade(Dict(list(1, 2), list(Rational(1, 2), Rational(1, 2))))
#' rv_attr(FH, "set")
#' @name handmade_distributions
NULL

#' @rdname handmade_distributions
#' @export
ContinuousDistributionHandmade <- function(pdf, set = NULL) {
  if (!is.null(set)) set <- as_sym(set)
  .ContinuousDistributionHandmade(as_sym(pdf), set)
}

#' @rdname handmade_distributions
#' @export
DiscreteDistributionHandmade <- function(pdf, set = Integers) .DiscreteDistributionHandmade(as_sym(pdf), as_sym(set))

#' @rdname handmade_distributions
#' @export
FiniteDistributionHandmade <- function(density) new_sym("FiniteDistributionHandmade", list(.as_fd(density)))

#' Sampling density
#'
#' `sampling_density(expr)` of `sympy.stats`: the frequencies of sampled
#' values. SymPy samples through scipy, numpy or pymc; this function follows
#' SymPy without those libraries and raises SymPy's error. Use [rv_sample()]
#' for samples drawn with R's random number generator.
#'
#' @param expr A random expression.
#' @param given_condition An optional condition.
#' @param library The sampling library SymPy would use.
#' @param numsamples Number of samples.
#' @param seed Ignored.
#' @param ... Ignored.
#' @return Does not return; raises a `ValueError` condition.
#' @examples
#' r <- tryCatch(sampling_density(Normal("X", 0, 1)), error = conditionMessage)
#' r
#' @export
sampling_density <- function(expr, given_condition = NULL, library = "scipy", numsamples = 1, seed = NULL, ...)
  .serr("ValueError", paste0("Failed to import ", library))

#' Matrix expectation, variance and cross covariance
#'
#' `ExpectationMatrix()`, `VarianceMatrix()` and `CrossCovarianceMatrix()` of
#' `sympy.stats`: the unevaluated expectation of a random matrix, and the
#' variance and cross covariance matrices of random vectors. [Expectation()],
#' [Variance()] and [Covariance()] build them for matrix arguments;
#' [expand()] applies linearity.
#'
#' @param expr,arg,arg1,arg2 Matrix expressions (vectors for the variance and
#'   cross covariance).
#' @param condition An optional condition.
#' @return A matrix expression.
#' @examples
#' k <- Symbol("k", integer = TRUE)
#' A <- MatrixSymbol("A", k, k)
#' X <- RandomMatrixSymbol("X", k, 1)
#' expand(ExpectationMatrix(A * X))
#' expand(VarianceMatrix(A * X))
#' @name matrix_expectation
NULL

#' @rdname matrix_expectation
#' @export
ExpectationMatrix <- function(expr, condition = NULL) .ExpectationMatrix(expr, condition)

#' @rdname matrix_expectation
#' @export
VarianceMatrix <- function(arg, condition = NULL) .VarianceMatrix(arg, condition)

#' @rdname matrix_expectation
#' @export
CrossCovarianceMatrix <- function(arg1, arg2, condition = NULL) .CrossCovarianceMatrix(arg1, arg2, condition)

# ---------------------------------------------------------------------------
# Shapes, attributes, doit and expand of the matrix forms

.stp_vm_shape <- function(sh) if (identical(sh[[2L]]$key, S_One$key)) list(sh[[1L]], sh[[1L]]) else list(sh[[2L]], sh[[2L]])
.stp_ccm_shape <- function(s1, s2) {
  if (identical(s1[[2L]]$key, S_One$key) && identical(s2[[2L]]$key, S_One$key)) list(s1[[1L]], s2[[1L]]) else list(S_One, S_One)
}
.me_shape_hooks$ExpectationMatrix <- function(e) .mshape(e$args[[1L]])
.me_shape_hooks$VarianceMatrix <- function(e) .stp_vm_shape(.mshape(e$args[[1L]]))
.me_shape_hooks$CrossCovarianceMatrix <- function(e) .stp_ccm_shape(.mshape(e$args[[1L]]), .mshape(e$args[[2L]]))
for (.stp_cls in c("ExpectationMatrix", "VarianceMatrix", "CrossCovarianceMatrix")) {
  .defm(.stp_cls, "shape", function(e) .pytuple(.mshape(e)))
  .defm(.stp_cls, "rows", function(e) .mshape(e)[[1L]])
  .defm(.stp_cls, "cols", function(e) .mshape(e)[[2L]])
}
rm(.stp_cls)

.stp_is_rms <- function(e) is_sym(e) && (!is.null(e$extra$rv_matrix) || .is_rs(e))
.stp_cond <- function(e, i) if (length(e$args) >= i) e$args[[i]] else NULL
# Mul.fromiter of scalars and matrices
.stp_mprod <- function(l) if (!length(l)) S_One else Reduce(function(a, b) a * b, l)
.stp_msum <- function(l) if (!length(l)) S_Zero else Reduce(function(a, b) a + b, l)
.stp_transpose <- function(e) if (.is_matexpr(e) || inherits(e, "symmatrix")) .mtranspose(e) else e
.stp_args <- function(e, types) if (is_sym(e) && e$type %in% types) e$args else NULL
.stp_zero <- function(sh) ZeroMatrix(sh[[1L]], sh[[2L]])

# Expectation applied to a matrix expression: the doit of Expectation for
# an explicit matrix, entry by entry
.stp_em_doit <- function(e, deep = TRUE) {
  expr <- e$args[[1L]]
  condition <- .stp_cond(e, 2L)
  if (deep) expr <- doit(expr)
  if (!is_random(expr) || (is_sym(expr) && expr$type %in% c("Expectation", "ExpectationMatrix")))
    return(if (deep && is_sym(expr) && expr$type == "ImmutableMatrix") as_explicit(expr) else expr)
  if (!is.null(condition)) return(doit(Expectation(given(expr, condition))))
  if (.is_explicit(expr)) {
    M <- if (inherits(expr, "symmatrix")) expr else as_explicit(expr)
    return(.mat_map(M, function(z) doit(Expectation(z))))
  }
  e
}
register_doit("ExpectationMatrix", function(e, deep = TRUE) .stp_em_doit(e, deep))

# expand() of a matrix product: MatMul distributes over MatAdd
.stp_mexpand <- function(e) {
  if (!(is_sym(e) && e$type == "MatMul")) return(e)
  terms <- list(list())
  for (a in e$args) {
    a <- .stp_mexpand(a)
    parts <- if (is_sym(a) && a$type == "MatAdd") a$args else list(a)
    terms <- unlist(lapply(terms, function(t) lapply(parts, function(q) c(t, list(q)))), recursive = FALSE)
  }
  if (length(terms) == 1L) return(e)
  .stp_msum(lapply(terms, .stp_mprod))
}

.stp_em_expand <- function(e) {
  expr <- e$args[[1L]]
  condition <- .stp_cond(e, 2L)
  if (!is_random(expr)) return(expr)
  each <- function(a) .stp_matrix_expand(Expectation(a, condition))
  if (expr$type %in% c("Add", "MatAdd")) return(.stp_msum(lapply(expr$args, each)))
  ee <- .stp_mexpand(expr)
  if (is_sym(ee) && ee$type %in% c("Add", "MatAdd")) return(.stp_msum(lapply(ee$args, each)))
  if (expr$type %in% c("Mul", "MatMul")) {
    rv <- list(); nonrv <- list(); postnon <- list()
    for (a in expr$args) {
      if (is_random(a)) {
        if (length(rv)) rv <- c(rv, postnon) else nonrv <- c(nonrv, postnon)
        postnon <- list()
        rv[[length(rv) + 1L]] <- a
      } else if (.is_matexpr(a) || inherits(a, "symmatrix")) {
        postnon[[length(postnon) + 1L]] <- a
      } else {
        nonrv[[length(nonrv) + 1L]] <- a
      }
    }
    if (!length(nonrv)) return(e)
    return(.stp_mprod(nonrv) * Expectation(.stp_mprod(rv), condition) * .stp_mprod(postnon))
  }
  e
}

.stp_vm_expand <- function(e) {
  arg <- e$args[[1L]]
  condition <- .stp_cond(e, 2L)
  if (!is_random(arg)) return(.stp_zero(.mshape(e)))
  if (.stp_is_rms(arg)) return(e)
  if (arg$type %in% c("Add", "MatAdd")) {
    rv <- Filter(is_random, arg$args)
    variances <- .stp_msum(lapply(rv, function(xv) .stp_matrix_expand(Variance(xv, condition))))
    covs <- list()
    if (length(rv) >= 2L) for (i in 1:(length(rv) - 1L)) for (j in (i + 1L):length(rv))
      covs[[length(covs) + 1L]] <- 2 * .stp_matrix_expand(Covariance(rv[[i]], rv[[j]], condition))
    return(variances + .stp_msum(covs))
  }
  if (arg$type %in% c("Mul", "MatMul")) {
    rv <- Filter(is_random, arg$args)
    nonrv <- Filter(function(a) !is_random(a), arg$args)
    if (!length(rv)) return(.stp_zero(.mshape(e)))
    if (!length(nonrv) || length(rv) > 1L) return(e)
    nr <- .stp_mprod(nonrv)
    return(nr * Variance(.stp_mprod(rv), condition) * .stp_transpose(nr))
  }
  e
}

.stp_mul_nonrv_rv <- function(m) {
  list(.stp_mprod(Filter(function(a) !is_random(a), m$args)), .stp_mprod(Filter(is_random, m$args)))
}
.stp_expand_single <- function(expr) {
  if (.stp_is_rms(expr)) return(list(list(S_One, expr)))
  if (expr$type %in% c("Add", "MatAdd")) {
    out <- list()
    for (a in expr$args) {
      if (a$type %in% c("Mul", "MatMul")) out[[length(out) + 1L]] <- .stp_mul_nonrv_rv(a)
      else if (is_random(a)) out[[length(out) + 1L]] <- list(S_One, a)
    }
    return(out)
  }
  if (expr$type %in% c("Mul", "MatMul")) return(list(.stp_mul_nonrv_rv(expr)))
  if (is_random(expr)) return(list(list(S_One, expr)))
  .serr("TypeError", "'NoneType' object is not iterable")
}

.stp_ccm_expand <- function(e) {
  arg1 <- e$args[[1L]]; arg2 <- e$args[[2L]]
  condition <- .stp_cond(e, 3L)
  if (identical(arg1$key, arg2$key)) return(.stp_vm_expand(VarianceMatrix(arg1, condition)))
  if (!is_random(arg1) || !is_random(arg2)) return(.stp_zero(.mshape(e)))
  if (.stp_is_rms(arg1) && .stp_is_rms(arg2)) return(CrossCovarianceMatrix(arg1, arg2, condition))
  l1 <- .stp_expand_single(expand(arg1))
  l2 <- .stp_expand_single(expand(arg2))
  addends <- list()
  for (p1 in l1) for (p2 in l2)
    addends[[length(addends) + 1L]] <- p1[[1L]] * CrossCovarianceMatrix(p1[[2L]], p2[[2L]], condition) * .stp_transpose(p2[[1L]])
  .stp_msum(addends)
}

# expand() of the matrix forms (their Python expand methods)
.stp_matrix_expand <- function(e) {
  if (!is_sym(e)) return(e)
  switch(e$type,
         ExpectationMatrix = .stp_em_expand(e),
         VarianceMatrix = .stp_vm_expand(e),
         CrossCovarianceMatrix = .stp_ccm_expand(e),
         e)
}
