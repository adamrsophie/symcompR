# Unevaluated probability objects (sympy/stats/symbolic_probability.py) and
# the moment functions of sympy/stats/rv_interface.py.

.is_default_pspace <- function(X) {
  ps <- tryCatch(pspace(X), error = function(e) NULL)
  !is.null(ps) && is_sym(ps) && ps$type == "PSpace"
}

.cond_args <- function(condition) if (is.null(condition)) list() else list(.as_bool_sym(as_sym(condition)))

#' Unevaluated probability expressions
#'
#' `Probability()`, `Expectation()`, `Variance()`, `Covariance()`,
#' `Moment()` and `CentralMoment()` build SymPy's unevaluated objects.
#' Evaluate them with [doit()], or rewrite them with [rewrite()] to
#' `"Expectation"`, `"Probability"`, `"Integral"` or `"Sum"`; [expand()]
#' applies linearity.
#'
#' @param prob,condition A condition and an optional given condition.
#' @param expr,arg,arg1,arg2,X Random expressions.
#' @param n The order of a moment.
#' @param c The center of a moment.
#' @param evaluate For `Covariance()`, `FALSE` keeps the argument order.
#' @param ... Ignored.
#' @return A `sym` expression.
#' @examples
#' X <- Normal("X", 0, 1)
#' e <- Expectation(X^2)
#' doit(e)
#' rewrite(Variance(X), "Expectation")
#' @name symbolic_probability
NULL

#' @rdname symbolic_probability
#' @export
Probability <- function(prob, condition = NULL, ...) new_sym("Probability", c(list(.as_bool_sym(as_sym(prob))), .cond_args(condition)))

#' @rdname symbolic_probability
#' @export
Expectation <- function(expr, condition = NULL, ...) {
  if (inherits(expr, "symmatrix") || (is_sym(expr) && .is_matrix_expr(expr))) return(.ExpectationMatrix(expr, condition))
  expr <- as_sym(expr)
  if (is.null(condition) && !is_random(expr)) return(expr)
  new_sym("Expectation", c(list(expr), .cond_args(condition)))
}

#' @rdname symbolic_probability
#' @export
Variance <- function(arg, condition = NULL, ...) {
  if (inherits(arg, "symmatrix") || (is_sym(arg) && .is_matrix_expr(arg))) return(.VarianceMatrix(arg, condition))
  new_sym("Variance", c(list(as_sym(arg)), .cond_args(condition)))
}

#' @rdname symbolic_probability
#' @export
Covariance <- function(arg1, arg2, condition = NULL, evaluate = TRUE, ...) {
  if (inherits(arg1, "symmatrix") || inherits(arg2, "symmatrix") || (is_sym(arg1) && .is_matrix_expr(arg1)) || (is_sym(arg2) && .is_matrix_expr(arg2)))
    return(.CrossCovarianceMatrix(arg1, arg2, condition))
  arg1 <- as_sym(arg1); arg2 <- as_sym(arg2)
  if (evaluate) { s <- .sorted_basic(list(arg1, arg2)); arg1 <- s[[1L]]; arg2 <- s[[2L]] }
  new_sym("Covariance", c(list(arg1, arg2), .cond_args(condition)))
}

#' @rdname symbolic_probability
#' @export
Moment <- function(X, n, c = 0, condition = NULL, ...) new_sym("Moment", c(list(as_sym(X), as_sym(n), as_sym(c)), .cond_args(condition)))

#' @rdname symbolic_probability
#' @export
CentralMoment <- function(X, n, condition = NULL, ...) new_sym("CentralMoment", c(list(as_sym(X), as_sym(n)), .cond_args(condition)))

.is_matrix_expr <- function(e) isTRUE(tryCatch(.mexpr_is_matrix(e), error = function(err) FALSE))
.mexpr_is_matrix <- function(e) e$type %in% c("MatrixSymbol", "MatAdd", "MatMul", "MatPow", "Transpose", "Inverse", "ImmutableMatrix", "Identity", "ZeroMatrix", "ExpectationMatrix", "VarianceMatrix", "CrossCovarianceMatrix") || !is.null(e$extra$rv_matrix)

# The matrix forms of Expectation, Variance and Covariance
# (sympy/stats/symbolic_multivariate_probability.py), kept unevaluated
.is_vec_shape <- function(sh) any(vapply(sh, function(d) identical(d$key, S_One$key), NA))
.ExpectationMatrix <- function(expr, condition = NULL) {
  expr <- .as_mexpr(expr)
  if (is.null(condition) && !is_random(expr)) return(.unwrap_explicit(expr))
  new_sym("ExpectationMatrix", c(list(expr), .cond_args(condition)))
}
.VarianceMatrix <- function(arg, condition = NULL) {
  arg <- .as_mexpr(arg)
  if (!.is_vec_shape(.mshape(arg))) .serr("ShapeError", "Expression is not a vector")
  new_sym("VarianceMatrix", c(list(arg), .cond_args(condition)))
}
.CrossCovarianceMatrix <- function(arg1, arg2, condition = NULL) {
  arg1 <- .as_mexpr(arg1); arg2 <- .as_mexpr(arg2)
  s1 <- .mshape(arg1); s2 <- .mshape(arg2)
  if (!.is_vec_shape(s1) || !.is_vec_shape(s2) || !identical(s1[[2L]]$key, s2[[2L]]$key))
    .serr("ShapeError", "Expression is not a vector")
  new_sym("CrossCovarianceMatrix", c(list(arg1, arg2), .cond_args(condition)))
}

.cond_of <- function(e, i) if (length(e$args) >= i) e$args[[i]] else NULL

# ---------------------------------------------------------------------------
# Probability.doit

.probability_doit <- function(p, evaluate = TRUE, ...) {
  condition <- p$args[[1L]]
  given_condition <- .cond_of(p, 2L)
  if (condition$type == "Not")
    return(.sub(S_One, .probability_doit(Probability(condition$args[[1L]], given_condition), evaluate = evaluate)))
  if (.has_rvindexed(condition))
    return(.callm(pspace(condition), "probability", condition, given_condition, evaluate = evaluate))
  if (!is.null(given_condition) && .is_rs(given_condition)) {
    condrv <- random_symbols(condition)
    if (length(condrv) == 1L && identical(condrv[[1L]]$key, given_condition$key))
      return(new_sym("BernoulliDistribution", list(.probability_doit(Probability(condition)), S_Zero, S_One)))
    if (any(vapply(condrv, function(rv) isTRUE(dependent(rv, given_condition)), NA)))
      return(Probability(condition, given_condition))
    return(.probability_doit(Probability(condition)))
  }
  if (!is.null(given_condition) && !.is_cond(given_condition))
    .serr("ValueError", paste0(.str(given_condition), " is not a relational or combination of relationals"))
  if ((!is.null(given_condition) && .is_false(given_condition)) || .is_false(condition)) return(S_Zero)
  if (!.is_cond(condition) || .is_rs(condition))
    .serr("ValueError", paste0(.str(condition), " is not a relational or combination of relationals"))
  if (.is_true(condition)) return(S_One)
  if (!is.null(given_condition)) return(.probability_doit(Probability(given(condition, given_condition))))
  ps <- pspace(condition)
  if (ps$type == "PSpace") return(Probability(condition, given_condition))
  result <- .callm(ps, "probability", condition)
  if (evaluate && is_sym(result)) return(doit(result))
  result
}
register_doit("Probability", function(e, deep = TRUE) .probability_doit(e))

.has_rvindexed <- function(e) is_sym(e) && any(vapply(.rv_atoms(e), function(r) identical(r$extra$rv_class, "RandomIndexedSymbol"), NA))

# ---------------------------------------------------------------------------
# Expectation.doit

.expectation_doit <- function(e, deep = TRUE, evaluate = TRUE) {
  if (e$type != "Expectation") return(if (deep) doit(e) else e)
  condition <- .cond_of(e, 2L)
  expr <- e$args[[1L]]
  if (deep) expr <- doit(expr)
  if (!is_random(expr) || (is_sym(expr) && expr$type == "Expectation")) return(expr)
  if (.has_rvindexed(expr)) return(.callm(pspace(expr), "compute_expectation", expr, condition))
  if (!is.null(condition)) return(.expectation_doit(Expectation(given(expr, condition)), deep = deep, evaluate = evaluate))
  if (expr$type == "Add")
    return(.add(lapply(expr$args, function(a) {
      if (a$type == "Expectation") return(Expectation(a, condition))
      ea <- Expectation(a, condition)
      if (is_sym(ea) && ea$type == "Expectation") .expectation_doit(ea, deep = deep, evaluate = evaluate) else doit(ea)
    })))
  if (expr$type == "Mul" && any(vapply(.all_subexprs(expr), function(z) z$type == "Expectation", NA))) return(expr)
  ps <- pspace(expr)
  if (ps$type == "PSpace") return(Expectation(expr))
  result <- .callm(ps, "compute_expectation", expr, evaluate = evaluate)
  if (evaluate && is_sym(result)) return(doit(result))
  result
}
register_doit("Expectation", function(e, deep = TRUE) .expectation_doit(e, deep = deep))

# Moment.doit and CentralMoment.doit
register_doit("Moment", function(e, deep = TRUE) doit(.rw_moment_E(e$args[[1L]], e$args[[2L]], e$args[[3L]], .cond_of(e, 4L))))
register_doit("CentralMoment", function(e, deep = TRUE) doit(.rw_cmoment_E(e$args[[1L]], e$args[[2L]], .cond_of(e, 3L))))

# ---------------------------------------------------------------------------
# expand()

.stats_expand <- function(e) {
  switch(e$type,
         Expectation = .expectation_expand(e),
         Variance = .variance_expand(e),
         Covariance = .covariance_expand(e),
         e)
}

.expectation_expand <- function(e) {
  expr <- e$args[[1L]]
  condition <- .cond_of(e, 2L)
  if (!is_random(expr)) return(expr)
  ex_each <- function(a) { r <- Expectation(a, condition); if (is_sym(r) && r$type == "Expectation") .expectation_expand(r) else r }
  if (expr$type == "Add") return(.add(lapply(expr$args, ex_each)))
  ee <- expand(expr)
  if (ee$type == "Add") return(.add(lapply(ee$args, ex_each)))
  if (expr$type == "Mul") {
    rv <- Filter(is_random, expr$args)
    nonrv <- Filter(function(a) !is_random(a), expr$args)
    return(.mul(c(nonrv, list(Expectation(.mul(rv), condition)))))
  }
  e
}

.variance_expand <- function(e) {
  arg <- e$args[[1L]]
  condition <- .cond_of(e, 2L)
  if (!is_random(arg)) return(S_Zero)
  if (.is_rs(arg)) return(e)
  if (arg$type == "Add") {
    rv <- Filter(is_random, arg$args)
    variances <- .add(lapply(rv, function(xv) .variance_expand(Variance(xv, condition))))
    covs <- list()
    if (length(rv) >= 2L) for (i in 1:(length(rv) - 1L)) for (j in (i + 1L):length(rv))
      covs[[length(covs) + 1L]] <- .mul(list(as_sym(2), .covariance_expand(Covariance(rv[[i]], rv[[j]], condition))))
    return(.add(list(variances, .add(covs))))
  }
  if (arg$type == "Mul") {
    rv <- Filter(is_random, arg$args)
    nonrv <- lapply(Filter(function(a) !is_random(a), arg$args), function(a) a^2)
    if (!length(rv)) return(S_Zero)
    return(.mul(c(nonrv, list(Variance(.mul(rv), condition)))))
  }
  e
}

.cov_single <- function(expr) {
  mul_split <- function(m) {
    rv <- Filter(is_random, m$args); nonrv <- Filter(function(a) !is_random(a), m$args)
    list(.mul(nonrv), .mul(rv))
  }
  if (.is_rs(expr)) return(list(list(S_One, expr)))
  if (expr$type == "Add") {
    out <- list()
    for (a in expr$args) {
      if (a$type == "Mul") out[[length(out) + 1L]] <- mul_split(a)
      else if (is_random(a)) out[[length(out) + 1L]] <- list(S_One, a)
    }
    return(out)
  }
  if (expr$type == "Mul") return(list(mul_split(expr)))
  if (is_random(expr)) return(list(list(S_One, expr)))
  list()
}

.covariance_expand <- function(e) {
  a1 <- e$args[[1L]]; a2 <- e$args[[2L]]
  condition <- .cond_of(e, 3L)
  if (identical(a1$key, a2$key)) return(.variance_expand(Variance(a1, condition)))
  if (!is_random(a1) || !is_random(a2)) return(S_Zero)
  s <- .sorted_basic(list(a1, a2)); a1 <- s[[1L]]; a2 <- s[[2L]]
  if (.is_rs(a1) && .is_rs(a2)) return(Covariance(a1, a2, condition))
  l1 <- .cov_single(expand(a1)); l2 <- .cov_single(expand(a2))
  adds <- list()
  for (p in l1) for (q in l2) {
    rr <- .sorted_basic(list(p[[2L]], q[[2L]]))
    adds[[length(adds) + 1L]] <- .mul(list(p[[1L]], q[[1L]], Covariance(rr[[1L]], rr[[2L]], condition)))
  }
  .add(adds)
}

# ---------------------------------------------------------------------------
# rewrite() rules

.rw_moment_E <- function(X, n, c = S_Zero, condition = NULL) Expectation((X - c)^n, condition)
.rw_cmoment_E <- function(X, n, condition = NULL) {
  mu <- Expectation(X, condition)
  .rw_moment_E(X, n, mu, condition)
}

.rw_E_P <- function(arg, condition = NULL) {
  rvs <- .rv_atoms(arg)
  if (length(rvs) > 1L) .serr("NotImplementedError", "")
  if (!length(rvs)) return(arg)
  rv <- rvs[[1L]]
  ps <- .rs_pspace(rv)
  nm <- .rs_symbol(rv)$name
  symbol <- if (grepl("^[[:upper:]]", nm)) Symbol(tolower(nm)) else Symbol(paste0(nm, "_1"))
  body <- .mul(list(sym_replace(arg, rv, symbol), Probability(Eq(rv, symbol), condition)))
  dset <- .dom_set(.ps_domain(ps))
  if (.ps_cont(ps)) return(Integral(body, list(symbol, set_inf(dset), set_sup(dset))))
  if (.ps_fin(ps)) .serr("NotImplementedError", "")
  Sum(body, list(symbol, set_inf(dset), set_sup(.callm(ps, "set"))))
}

.rw_P_int <- function(arg, condition = NULL) .probability_doit(Probability(arg, condition), evaluate = FALSE)
.rw_E_int <- function(arg, condition = NULL) .expectation_doit(Expectation(arg, condition), deep = FALSE, evaluate = FALSE)
.rw_V_E <- function(arg, condition = NULL) Expectation(arg^2, condition) - Expectation(arg, condition)^2
.rw_C_E <- function(arg1, arg2, condition = NULL) Expectation(arg1 * arg2, condition) - Expectation(arg1, condition) * Expectation(arg2, condition)

local({
  add_rule <- function(target, type, fn) {
    cur <- .rewrite_rules[[target]]
    if (is.null(cur)) cur <- list()
    cur[[type]] <- fn
    .rewrite_rules[[target]] <<- cur
  }
  for (tg in c("Integral", "Sum")) {
    add_rule(tg, "Probability", .rw_P_int)
    add_rule(tg, "Expectation", .rw_E_int)
    add_rule(tg, "Variance", function(arg, condition = NULL) variance(arg, condition, evaluate = FALSE))
    add_rule(tg, "Covariance", function(arg1, arg2, condition = NULL) covariance(arg1, arg2, condition, evaluate = FALSE))
    add_rule(tg, "Moment", function(X, n, c = S_Zero, condition = NULL) rewrite(.rw_moment_E(X, n, c, condition), "Integral"))
    add_rule(tg, "CentralMoment", function(X, n, condition = NULL) rewrite(.rw_cmoment_E(X, n, condition), "Integral"))
  }
  add_rule("Expectation", "Variance", .rw_V_E)
  add_rule("Expectation", "Covariance", .rw_C_E)
  add_rule("Expectation", "Moment", .rw_moment_E)
  add_rule("Expectation", "CentralMoment", .rw_cmoment_E)
  add_rule("Probability", "Expectation", .rw_E_P)
  add_rule("Probability", "Variance", function(arg, condition = NULL) rewrite(.rw_V_E(arg, condition), "Probability"))
  add_rule("Probability", "Covariance", function(arg1, arg2, condition = NULL) rewrite(.rw_C_E(arg1, arg2, condition), "Probability"))
  add_rule("Probability", "Moment", function(X, n, c = S_Zero, condition = NULL) rewrite(.rw_moment_E(X, n, c, condition), "Probability"))
  add_rule("Probability", "CentralMoment", function(X, n, condition = NULL) rewrite(.rw_cmoment_E(X, n, condition), "Probability"))
})

# ---------------------------------------------------------------------------
# rv_interface

#' Moments and summary statistics of random expressions
#'
#' `variance()`, `standard_deviation()` (alias `std()`), `covariance()`,
#' `correlation()`, `moment()`, `cmoment()` (central), `smoment()`
#' (standardized), `skewness()`, `kurtosis()`, `factorial_moment()`,
#' `coskewness()`, `entropy()` (alias `H()`) and `rv_median()` (SymPy's
#' `median()`), following `sympy.stats`.
#'
#' @param X,Y,Z Random expressions.
#' @param n The order of a moment.
#' @param c The center of a raw moment.
#' @param condition An optional condition.
#' @param evaluate `FALSE` returns unevaluated integrals or sums.
#' @param ... Passed on to [E()].
#' @return An expression; `rv_median()` returns a set.
#' @examples
#' X <- Normal("X", 0, 1)
#' variance(X)
#' moment(X, 4)
#' skewness(Exponential("T", 1))
#' rv_median(Die("D", 6))
#' @name rv_moments
NULL

#' @rdname rv_moments
#' @export
moment <- function(X, n, c = 0, condition = NULL, evaluate = TRUE, ...) {
  m <- Moment(X, n, c, condition)
  if (evaluate) return(doit(m))
  rewrite(m, "Integral")
}

#' @rdname rv_moments
#' @export
variance <- function(X, condition = NULL, ...) {
  if (is_random(X) && .is_default_pspace(X)) return(Variance(X, condition))
  cmoment(X, 2, condition, ...)
}

#' @rdname rv_moments
#' @export
standard_deviation <- function(X, condition = NULL, ...) sqrt(as_sym(variance(X, condition, ...)))

#' @rdname rv_moments
#' @export
std <- standard_deviation

#' @rdname rv_moments
#' @param b Logarithm base for the entropy.
#' @export
entropy <- function(X, condition = NULL, b = NULL, ...) {
  # sympy.physics.quantum.entropy of a density operator or a matrix
  if (inherits(X, "symmatrix") || (is_sym(X) && .qisa(X, "QDensity"))) return(.qentropy(X))
  pdf <- rv_density(X, condition, ...)
  base <- if (is.null(b)) exp(S_One) else as_sym(b)
  lg <- function(v) if (identical(base$key, S_E$key)) log(v) else log(v) / log(base)
  if (is_sym(pdf) && pdf$type == "FiniteDensity")
    return(.add(c(list(S_Zero), lapply(.fd_vals(pdf), function(p) -p * lg(p)))))
  E(-lg(.rv_lcall(pdf, X)))
}

#' @rdname rv_moments
#' @export
H <- entropy

#' @rdname rv_moments
#' @export
covariance <- function(X, Y, condition = NULL, ...) {
  if ((is_random(X) && .is_default_pspace(X)) || (is_random(Y) && .is_default_pspace(Y))) return(Covariance(X, Y, condition))
  E((X - E(X, condition, ...)) * (Y - E(Y, condition, ...)), condition, ...)
}

#' @rdname rv_moments
#' @export
correlation <- function(X, Y, condition = NULL, ...)
  as_sym(covariance(X, Y, condition, ...)) / (std(X, condition, ...) * std(Y, condition, ...))

#' @rdname rv_moments
#' @export
cmoment <- function(X, n, condition = NULL, evaluate = TRUE, ...) {
  m <- CentralMoment(X, n, condition)
  if (evaluate) return(doit(m))
  rewrite(m, "Integral")
}

#' @rdname rv_moments
#' @export
smoment <- function(X, n, condition = NULL, ...) {
  sigma <- std(X, condition, ...)
  (1 / sigma)^n * cmoment(X, n, condition, ...)
}

#' @rdname rv_moments
#' @export
skewness <- function(X, condition = NULL, ...) smoment(X, 3, condition = condition, ...)

#' @rdname rv_moments
#' @export
kurtosis <- function(X, condition = NULL, ...) smoment(X, 4, condition = condition, ...)

#' @rdname rv_moments
#' @export
factorial_moment <- function(X, n, condition = NULL, ...) E(FallingFactorial(X, n), condition = condition, ...)

#' @rdname rv_moments
#' @export
rv_median <- function(X, evaluate = TRUE, ...) {
  if (!is_random(X)) return(X)
  ps <- pspace(X)
  if (.isa(ps, "FinitePSpace")) {
    cdf <- .callm(ps, "compute_cdf", X)
    res <- list()
    for (kv in cdf$args) {
      key <- kv$args[[1L]]; value <- kv$args[[2L]]
      c1 <- .truth(value >= S_Half)
      if (c1 && .truth(1 - value + .callm(ps, "probability", Eq(X, key)) >= S_Half)) res[[length(res) + 1L]] <- key
    }
    return(if (length(res)) do.call(FiniteSet, res) else S_EmptySet)
  }
  if (.isa(ps, "ContinuousPSpace") || .isa(ps, "DiscretePSpace")) {
    cdf <- .callm(ps, "compute_cdf", X)
    x <- Dummy("x")
    return(solveset(piecewise_fold(.rv_lcall(cdf, x) - S_Half), x, .callm(ps, "set")))
  }
  .serr("NotImplementedError", paste0("The median of ", .str(ps), " is not implemented."))
}

#' @rdname rv_moments
#' @export
coskewness <- function(X, Y, Z, condition = NULL, ...) {
  num <- E((X - E(X, condition, ...)) * (Y - E(Y, condition, ...)) * (Z - E(Z, condition, ...)), condition, ...)
  den <- std(X, condition, ...) * std(Y, condition, ...) * std(Z, condition, ...)
  as_sym(num) / den
}
