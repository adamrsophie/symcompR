# Joint (multivariate) random variables: sympy/stats/joint_rv.py and
# sympy/stats/joint_rv_types.py.

.defclass("JointDistribution", "Distribution")
.defclass("JointPSpace", "ProductPSpace")
.defclass("MarginalDistribution", "Distribution")

# list arguments become (column) ImmutableMatrix nodes, as in
# JointDistribution.__new__
.jarg <- function(a) {
  if (is.function(a) && !is.null(attr(a, "base"))) return(attr(a, "base"))
  if (inherits(a, "symmatrix")) return(.as_block(a))
  if (is.list(a) && !is_sym(a)) return(.as_block(Matrix(a)))
  as_sym(a)
}
# explicit matrix of an argument (or the matrix expression itself)
.mx <- function(a) if (is_sym(a) && a$type == "ImmutableMatrix") as_explicit(a) else a
.is_mat <- function(a) inherits(a, "symmatrix") || (is_sym(a) && .is_matexpr(a))
.mentries <- function(a) {
  m <- .mx(a)
  if (inherits(m, "symmatrix")) return(m$data)
  if (is.list(m) && !is_sym(m)) return(lapply(m, as_sym))
  .serr("TypeError", paste0("'", .sympy_class(m), "' object is not iterable"))
}
.mrows_n <- function(a) { m <- .mx(a); if (inherits(m, "symmatrix")) .int(m$nrow) else .mshape(m)[[1L]] }

.rv_set_pow <- function(s, k) {
  k <- as_sym(k)
  if (k$type != "Integer" || k$value < 0) .serr("ValueError", paste0(.str(k), ": Exponent must be a positive Integer"))
  do.call(ProductSet, rep(list(s), as.integer(k$value)))
}

.JointPSpace <- function(sym, dist) {
  if (.isa(dist, "SingleContinuousDistribution")) return(.SingleContinuousPSpace(sym, dist))
  if (.isa(dist, "SingleDiscreteDistribution")) return(.SingleDiscretePSpace(sym, dist))
  new_sym("JointPSpace", list(.symbol_converter(sym), dist))
}

.jps_dist <- function(ps) ps$args[[2L]]
.defm("JointPSpace", "symbol", function(ps) ps$args[[1L]])
.defm("JointPSpace", "distribution", function(ps) ps$args[[2L]])
.defm("JointPSpace", "set", function(ps) .dom_set(.ps_domain(ps)))
.defm("JointPSpace", "value", function(ps) .RandomSymbol(ps$args[[1L]], ps, cls = "JointRandomSymbol"))
.defm("JointPSpace", "component_count", function(ps) {
  s <- .dist_set(.jps_dist(ps))
  if (is_sym(s) && s$type == "ProductSet") return(.int(length(s$args)))
  if (is_sym(s) && s$type == "Product") return(s$args[[2L]]$args[[3L]])
  S_One
})
.jidx <- function(sym, i) Indexed(sym$name, i)
.defm("JointPSpace", "pdf", function(ps) {
  n <- .int_val(.callm(ps, "component_count"))
  syms <- lapply(seq_len(n) - 1L, function(i) .jidx(ps$args[[1L]], i))
  do.call(.callm, c(list(.jps_dist(ps), "pdf"), syms))
})
.defm("JointPSpace", "domain", function(ps) {
  d <- .jps_dist(ps)
  rvs <- random_symbols(d)
  if (!length(rvs)) return(new_sym("SingleDomain", list(ps$args[[1L]], .dist_set(d))))
  do.call(.ProductDomain, lapply(rvs, function(rv) .ps_domain(.rs_pspace(rv))))
})
.defm("JointPSpace", "component_domain", function(ps, index) .callm(ps, "set")$args[[.int_val(as_sym(index)) + 1L]])
.defm("JointPSpace", "where", function(ps, condition) .serr("NotImplementedError", ""))
.defm("JointPSpace", "compute_density", function(ps, expr, ...) .serr("NotImplementedError", ""))
.defm("JointPSpace", "probability", function(ps, condition, ...) .serr("NotImplementedError", ""))

# set limits: an Interval-like set gives its ends
.lim_from_set <- function(v, s) list(v, set_inf(s), set_sup(s))

.defm("JointPSpace", "marginal_distribution", function(ps, ...) {
  indices <- vapply(list(...), function(i) as.integer(.int_val(as_sym(i))), 0L)
  count <- .callm(ps, "component_count")
  if (length(free_symbols(count))) .serr("ValueError", "Marginal distributions cannot be computed for symbolic dimensions. It is a work under progress.")
  count <- .int_val(count)
  sym0 <- ps$args[[1L]]
  orig <- lapply(seq_len(count) - 1L, function(i) .jidx(sym0, i))
  all_syms <- lapply(orig, function(o) Symbol(.str(o)))
  sym <- lapply(indices, function(i) Symbol(.str(.jidx(sym0, i))))
  dist <- .jps_dist(ps)
  dset <- .dist_set(dist)
  limits <- list()
  for (i in seq_len(count) - 1L) if (!(i %in% indices)) limits[[length(limits) + 1L]] <- .lim_from_set(all_syms[[i + 1L]], dset$args[[i + 1L]])
  body <- do.call(.callm, c(list(dist, "pdf"), all_syms))
  if (isTRUE(.callm(dist, "is_Continuous"))) f <- Lambda(sym, do.call(integrate, c(list(body), limits)))
  else f <- Lambda(sym, do.call(summation, c(list(body), limits)))
  xreplace(f, lapply(seq_along(all_syms), function(k) list(all_syms[[k]], orig[[k]])))
})

.defm("JointPSpace", "compute_expectation", function(ps, expr, rvs = NULL, evaluate = FALSE, ...) {
  n <- .int_val(.callm(ps, "component_count"))
  value <- .callm(ps, "value")
  syms <- lapply(seq_len(n) - 1L, function(i) rv_getitem(value, i))
  if (is.null(rvs) || !length(rvs)) rvs <- syms
  if (!any(vapply(syms, function(s) .in_list(s, rvs), NA))) return(expr)
  expr <- expr * .callm(ps, "pdf")
  for (rv in rvs) {
    if (rv$type == "Indexed") expr <- xreplace(expr, list(list(rv, Indexed(.str(rv$args[[1L]]), rv$args[[2L]]))))
    else if (.is_rs(rv)) expr <- xreplace(expr, list(list(rv, .rs_symbol(rv))))
  }
  if (.in_list(value, random_symbols(expr)))
    .serr("NotImplementedError", "\nExpectations of expression with unindexed joint random symbols cannot\nbe calculated yet.")
  dset <- .dist_set(.jps_dist(ps))
  limits <- lapply(syms, function(rv) .lim_from_set(Indexed(.str(rv$args[[1L]]), rv$args[[2L]]), dset$args[[.int_val(rv$args[[2L]]) + 1L]]))
  do.call(Integral, c(list(expr), limits))
})

.joint_getitem <- function(x, key) {
  ps <- .rs_pspace(x)
  key <- as_sym(key)
  cc <- .callm(ps, "component_count")
  if (isTRUE(.truth_or_na(Le(cc, key))))
    .serr("ValueError", paste0("Index keys for ", .str(.rs_symbol(x)), " can only up to ", .str(cc - 1), "."))
  new_sym("Indexed", list(x, key))
}
.truth_or_na <- function(c) if (is.logical(c)) c else if (is_sym(c) && c$type == "BooleanAtom") c$name == "true" else NA

# ---------------------------------------------------------------------------
# JointDistribution

.defm("JointDistribution", "domain", function(d) .ProductDomain(.callm(d, "symbols")))
.defm("JointDistribution", "is_Continuous", function(d) FALSE)
.defm("JointDistribution", "is_Discrete", function(d) FALSE)
.defjd <- function(cls, argnames, continuous = TRUE, ...) {
  .defdist(cls, "JointDistribution", argnames, ...)
  .defm(cls, "is_Continuous", function(d) continuous)
  .defm(cls, "is_Discrete", function(d) !continuous)
}
.jd_new <- function(cls, args) new_sym(cls, lapply(args, .jarg))

# multivariate_rv(cls, sym, *args)
.mvrv <- function(cls, sym, args) {
  dist <- .jd_new(cls, args)
  ck <- .findm(dist, "check")
  if (!is.null(ck)) do.call(ck, dist$args)
  .ps_value(.JointPSpace(sym, dist))
}

#' @rdname joint_rv
#' @export
marginal_distribution <- function(rv, ...) {
  indices <- list(...)
  indices <- lapply(indices, function(i) if (is_sym(i) && i$type == "Indexed") i$args[[2L]] else i)
  ps <- .rs_pspace(rv)
  if (!length(indices)) .serr("ValueError", "At least one component for marginal density is needed.")
  dist <- .callm(ps, "distribution")
  md <- .findm(dist, "_marginal_distribution")
  if (!is.null(md)) return(md(dist, indices, .rs_symbol(rv)))
  do.call(.callm, c(list(ps, "marginal_distribution"), indices))
}

.py_is0 <- function(i) !is.list(i) && isTRUE(identical(as_sym(i)$key, S_Zero$key))

#' Joint random variables
#'
#' Multivariate distributions of `sympy.stats`: `MultivariateNormal()`,
#' `MultivariateLaplace()`, `MultivariateT()`, `NormalGamma()`,
#' `MultivariateBeta()` (alias `Dirichlet()`), `MultivariateEwens()`,
#' `GeneralizedMultivariateLogGamma()` and its `Omega` form,
#' `Multinomial()`, `NegativeMultinomial()` and the hand-made `JointRV()`.
#' `marginal_distribution(rv, ...)` gives the marginal density of the
#' listed components as a `Lambda`. Index a joint random symbol with
#' [rv_getitem()].
#'
#' @param name,syms,symbol,sym The name of the random symbol.
#' @param mu,sigma Mean (location) vector and covariance (shape)
#'   matrix, as lists or matrices.
#' @param v Degrees of freedom or shape parameters.
#' @param lamda,alpha,beta,delta,omega,theta,n,k0,pdf,set Distribution
#'   parameters, following SymPy.
#' @param rv A joint random symbol.
#' @param ... Components (indices or indexed symbols) for
#'   `marginal_distribution()`, or further parameters.
#' @return A random symbol, or a `Lambda` for `marginal_distribution()`.
#' @examples
#' A <- MultivariateNormal("A", list(1, 2), list(list(1, 0), list(0, 1)))
#' apply_lambda(rv_density(A), 1, 2)
#' M <- Multinomial("M", 3, list(Rational(1, 2), Rational(1, 2)))
#' rv_getitem(M, 0)
#' @name joint_rv
NULL

#' @rdname joint_rv
#' @export
JointRV <- function(symbol, pdf, set = NULL) {
  symbol <- .symbol_converter(symbol)
  pdf <- as_sym(pdf)
  syms <- Filter(function(i) i$type == "Indexed" && identical(i$args[[1L]]$key, attr(IndexedBase(symbol$name), "base")$key),
                 .all_subexprs(pdf))
  syms <- .uniq(syms)
  syms <- syms[order(vapply(syms, function(i) as.numeric(N(i$args[[2L]])), 0))]
  set <- .rv_set_pow(Reals, length(syms))
  dist <- new_sym("JointDistributionHandmade", list(Lambda(syms, pdf), set))
  jrv <- .ps_value(.JointPSpace(symbol, dist))
  if (length(random_symbols(pdf))) {
    dist <- new_sym("MarginalDistribution", list(dist, Tuple(jrv)))
    return(.ps_value(.JointPSpace(symbol, dist)))
  }
  jrv
}
.defjd("JointDistributionHandmade", c("pdf"),
  set = function(d) d$args[[2L]],
  pdf = function(d, ...) apply_lambda(d$args[[1L]], ...))

# MultivariateNormal

.is_psd <- function(M, strict = FALSE) {
  M <- .mx(M)
  if (!inherits(M, "symmatrix")) return(NA)
  n <- M$nrow
  minors <- list()
  if (strict) {
    for (k in seq_len(n)) minors[[k]] <- det(Matrix(lapply(seq_len(k), function(i) lapply(seq_len(k), function(j) .mget(M, i, j)))))
    res <- lapply(minors, function(m) .fuzzy_bool(m > 0))
  } else {
    subsets <- unlist(lapply(seq_len(n), function(k) utils::combn(n, k, simplify = FALSE)), recursive = FALSE)
    res <- lapply(subsets, function(s) .fuzzy_bool(det(Matrix(lapply(s, function(i) lapply(s, function(j) .mget(M, i, j))))) >= 0))
  }
  if (any(vapply(res, identical, NA, FALSE))) return(FALSE)
  if (all(vapply(res, isTRUE, NA))) return(TRUE)
  NA
}

.defjd("MultivariateNormalDistribution", c("mu", "sigma"),
  set = function(d) .rv_set_pow(Reals, .mrows_n(d$args[[1L]])),
  check = function(mu, sigma) {
    .value_check(Eq(.mrows_n(mu), .mrows_n(sigma)), "Size of the mean vector and covariance matrix are incorrect.")
    if (!(is_sym(sigma) && sigma$type == "MatrixSymbol"))
      .value_check(.is_psd(sigma), "The covariance matrix must be positive semi definite. ")
  },
  pdf = function(d, ...) {
    args <- list(...)
    mu <- .mx(d$args[[1L]]); sigma <- .mx(d$args[[2L]])
    k <- .mrows_n(d$args[[1L]])
    if (length(args) == 1L && .is_mat(args[[1L]])) xv <- args[[1L]] else xv <- Matrix(args)
    if (inherits(mu, "symmatrix") && inherits(sigma, "symmatrix") && inherits(xv, "symmatrix")) {
      x <- xv - mu
      q <- (Rational(-1, 2) * t(x)) %*% (inverse(sigma) %*% x)
      return(S_One / sqrt((2 * S_Pi)^k * det(sigma)) * exp(.mget(q, 1L, 1L)))
    }
    x <- MatAdd(xv, -mu)
    density <- S_One / sqrt((2 * S_Pi)^k * Determinant(sigma)) *
      exp(MatMul(Rational(-1, 2), Transpose(x), MatMul(Inverse(sigma), x)))
    # the density is a scalar times exp of a 1 x 1 matrix expression
    new_sym("MatrixElement", list(density, S_Zero, S_Zero))
  },
  `_marginal_distribution` = function(d, indices, sym) {
    idx <- vapply(indices, function(i) as.integer(.int_val(as_sym(i))), 0L)
    symv <- Matrix(lapply(idx, function(i) Indexed(sym$name, i)))
    mu <- .mx(d$args[[1L]]); sigma <- .mx(d$args[[2L]])
    k <- mu$nrow
    for (i in seq_len(k) - 1L) if (!(i %in% idx)) {
      mu <- .row_del(mu, i)
      sigma <- .row_del(.col_del(sigma, i), i)
    }
    dv <- mu - symv
    q <- (Rational(-1, 2) * t(dv)) %*% (inverse(sigma) %*% dv)
    Lambda(symv$data, S_One / sqrt((2 * S_Pi)^length(mu$data) * det(sigma)) * exp(.mget(q, 1L, 1L)))
  })
.row_del <- function(M, i) Matrix(lapply(setdiff(seq_len(M$nrow), i + 1L), function(r) lapply(seq_len(M$ncol), function(c) .mget(M, r, c))))
.col_del <- function(M, i) Matrix(lapply(seq_len(M$nrow), function(r) lapply(setdiff(seq_len(M$ncol), i + 1L), function(c) .mget(M, r, c))))

#' @rdname joint_rv
#' @export
MultivariateNormal <- function(name, mu, sigma) .mvrv("MultivariateNormalDistribution", name, list(mu, sigma))

.defjd("MultivariateLaplaceDistribution", c("mu", "sigma"),
  set = function(d) .rv_set_pow(Reals, .mrows_n(d$args[[1L]])),
  check = function(mu, sigma) {
    .value_check(Eq(.mrows_n(mu), .mrows_n(sigma)), "Size of the mean vector and covariance matrix are incorrect.")
    if (!(is_sym(sigma) && sigma$type == "MatrixSymbol"))
      .value_check(.is_psd(sigma, strict = TRUE), "The covariance matrix must be positive definite. ")
  },
  pdf = function(d, ...) {
    mu <- .mx(d$args[[1L]]); sigma <- .mx(d$args[[2L]])
    mu_T <- t(mu)
    k <- as_sym(mu$nrow)
    sigma_inv <- inverse(sigma)
    args <- Matrix(list(...))
    args_T <- t(args)
    x <- .mget(mu_T %*% sigma_inv %*% mu, 1L, 1L)
    y <- .mget(args_T %*% sigma_inv %*% args, 1L, 1L)
    v <- 1 - k / 2
    2 * (y / (2 + x))^(v / 2) * besselk(v, sqrt((2 + x) * y)) * exp(.mget(args_T %*% sigma_inv %*% mu, 1L, 1L)) /
      ((2 * S_Pi)^(k / 2) * sqrt(det(sigma)))
  })

#' @rdname joint_rv
#' @export
MultivariateLaplace <- function(name, mu, sigma) .mvrv("MultivariateLaplaceDistribution", name, list(mu, sigma))

.defjd("MultivariateTDistribution", c("mu", "shape_mat", "dof"),
  set = function(d) .rv_set_pow(Reals, .mrows_n(d$args[[1L]])),
  check = function(mu, sigma, v) {
    .value_check(Eq(.mrows_n(mu), .mrows_n(sigma)), "Size of the location vector and shape matrix are incorrect.")
    if (!(is_sym(sigma) && sigma$type == "MatrixSymbol"))
      .value_check(.is_psd(sigma, strict = TRUE), "The shape matrix must be positive definite. ")
  },
  pdf = function(d, ...) {
    mu <- .mx(d$args[[1L]]); sigma <- .mx(d$args[[2L]])
    v <- d$args[[3L]]
    k <- as_sym(mu$nrow)
    sigma_inv <- inverse(sigma)
    x <- Matrix(list(...)) - mu
    gamma((k + v) / 2) / (gamma(v / 2) * (v * S_Pi)^(k / 2) * sqrt(det(sigma))) *
      (1 + 1 / v * .mget(t(x) %*% sigma_inv %*% x, 1L, 1L))^((-v - k) / 2)
  })

#' @rdname joint_rv
#' @export
MultivariateT <- function(syms, mu, sigma, v) .mvrv("MultivariateTDistribution", syms, list(mu, sigma, v))

.defjd("NormalGammaDistribution", c("mu", "lamda", "alpha", "beta"),
  check = function(mu, lamda, alpha, beta) {
    .value_check(ask_fact(mu, "real"), "Location must be real.")
    .value_check(lamda > 0, "Lambda must be positive")
    .value_check(alpha > 0, "alpha must be positive")
    .value_check(beta > 0, "beta must be positive")
  },
  set = function(d) ProductSet(Reals, Interval(0, S_Infinity)),
  pdf = function(d, x, tau) {
    p <- .dl(d)
    p$beta^p$alpha * sqrt(p$lamda) / (gamma(p$alpha) * sqrt(2 * S_Pi)) * tau^(p$alpha - S_Half) *
      exp(-1 * p$beta * tau) * exp(-1 * (p$lamda * tau * (x - p$mu)^2) / as_sym(2))
  },
  `_marginal_distribution` = function(d, indices, ...) {
    sym <- list(...)
    if (length(indices) == 2L) return(do.call(.callm, c(list(d, "pdf"), sym)))
    p <- .dl(d)
    if (.py_is0(indices[[1L]])) {
      x <- sym[[1L]]
      v <- p$alpha - S_Half; mu <- p$mu; sigma <- p$beta / (p$lamda * p$alpha)
      return(Lambda(sym, gamma((v + 1) / 2) / (gamma(v / 2) * sqrt(S_Pi * v) * sigma) *
                      (1 + 1 / v * ((x - mu) / sigma)^2)^((-v - 1) / 2)))
    }
    Lambda(sym, .dist_pdf(new_sym("GammaDistribution", list(p$alpha, p$beta)), sym[[1L]]))
  })

#' @rdname joint_rv
#' @export
NormalGamma <- function(sym, mu, lamda, alpha, beta) .mvrv("NormalGammaDistribution", sym, list(mu, lamda, alpha, beta))

.defjd("MultivariateBetaDistribution", c("alpha"),
  check = function(alpha) {
    a <- .mentries(alpha)
    .value_check(length(a) >= 2L, "At least two categories should be passed.")
    for (a_k in a) .value_check(!identical(.fuzzy_bool(a_k > 0), FALSE), "Each concentration parameter should be positive.")
  },
  set = function(d) .rv_set_pow(Interval(0, 1), length(.mentries(d$args[[1L]]))),
  pdf = function(d, ...) {
    alpha <- .mentries(d$args[[1L]])
    syms <- list(...)
    B <- .mul(lapply(alpha, gamma)) / gamma(.add(alpha))
    .mul(lapply(seq_along(alpha), function(i) as_sym(syms[[i]])^(alpha[[i]] - 1))) / B
  })

#' @rdname joint_rv
#' @export
MultivariateBeta <- function(syms, ...) {
  alpha <- list(...)
  if (!(is.list(alpha[[1L]]) && !is_sym(alpha[[1L]])) && !inherits(alpha[[1L]], "symmatrix")) alpha <- list(alpha)
  .mvrv("MultivariateBetaDistribution", syms, list(alpha[[1L]]))
}

#' @rdname joint_rv
#' @export
Dirichlet <- MultivariateBeta

.defjd("MultivariateEwensDistribution", c("n", "theta"), continuous = FALSE,
  check = function(n, theta) {
    .value_check(n > 0, "sample size should be positive integer.")
    .value_check(ask_fact(theta, "positive"), "mutation rate should be positive.")
  },
  set = function(d) {
    n <- d$args[[1L]]
    if (n$type != "Integer") {
      i <- Symbol("i", integer = TRUE, positive = TRUE)
      return(Product(Intersection(Naturals0, Interval(0, floor(n / i))), list(i, 1, n)))
    }
    nn <- as.integer(n$value)
    sets <- c(list(Range(0, nn + 1L)), if (nn >= 2L) lapply(2:nn, function(i) Range(0, nn %/% i + 1L)) else list())
    do.call(ProductSet, sets)
  },
  pdf = function(d, ...) {
    syms <- list(...)
    n <- d$args[[1L]]; theta <- d$args[[2L]]
    condi <- n$type == "Integer"
    s1 <- syms[[1L]]
    if (is.function(s1) && !is.null(attr(s1, "base"))) s1 <- attr(s1, "base")
    isbase <- is_sym(s1) && s1$type == "IndexedBase"
    if (!(isbase || condi)) .serr("ValueError", "Please use IndexedBase object for syms as the dimension is symbolic")
    term_1 <- sym_factorial(n) / RisingFactorial(theta, n)
    if (condi) {
      nn <- as.integer(n$value)
      term_2 <- .mul(lapply(seq_len(nn) - 1L, function(j) { s <- as_sym(syms[[j + 1L]]); theta^s / ((j + 1L)^s * sym_factorial(s)) }))
      cond <- Eq(.add(c(list(S_Zero), lapply(seq_len(nn) - 1L, function(k) (k + 1L) * as_sym(syms[[k + 1L]])))), n)
      return(Piecewise(list(term_1 * term_2, cond), list(S_Zero, S_true)))
    }
    j <- Symbol("j", positive = TRUE, integer = TRUE); k <- Symbol("k", positive = TRUE, integer = TRUE)
    ix <- function(v) new_sym("Indexed", list(s1, v))
    term_2 <- Product(theta^ix(j) / ((j + 1)^ix(j) * sym_factorial(ix(j))), list(j, 0, n - 1))
    cond <- Eq(Sum((k + 1) * ix(k), list(k, 0, n - 1)), n)
    Piecewise(list(term_1 * term_2, cond), list(S_Zero, S_true))
  })

#' @rdname joint_rv
#' @export
MultivariateEwens <- function(syms, n, theta) .mvrv("MultivariateEwensDistribution", syms, list(n, theta))

.defjd("GeneralizedMultivariateLogGammaDistribution", c("delta", "v", "lamda", "mu"),
  check = function(delta, v, l, mu) {
    .value_check(list(delta >= 0, delta <= 1), "delta must be in range [0, 1].")
    .value_check(v > 0, "v must be positive")
    for (lk in .mentries(l)) .value_check(lk > 0, "lamda must be a positive vector.")
    for (muk in .mentries(mu)) .value_check(muk > 0, "mu must be a positive vector.")
    .value_check(length(.mentries(l)) > 1L, "the distribution should have at least two random variables.")
  },
  set = function(d) .rv_set_pow(Reals, length(.mentries(d$args[[3L]]))),
  pdf = function(d, ...) {
    y <- lapply(list(...), as_sym)
    dd <- d$args[[1L]]; v <- d$args[[2L]]; l <- .mentries(d$args[[3L]]); mu <- .mentries(d$args[[4L]])
    n <- Symbol("n", negative = FALSE, integer = TRUE)
    k <- length(l)
    sterm1 <- (1 - dd)^n / (gamma(v + n)^(k - 1) * gamma(v) * gamma(n + 1))
    sterm2 <- .mul(lapply(seq_len(k), function(i) mu[[i]] * l[[i]]^(-v - n)))
    term1 <- sterm1 * sterm2
    sterm3 <- (v + n) * .add(c(list(S_Zero), lapply(seq_len(k), function(i) mu[[i]] * y[[i]])))
    sterm4 <- .add(c(list(S_Zero), lapply(seq_len(k), function(i) exp(mu[[i]] * y[[i]]) / l[[i]])))
    term2 <- exp(sterm3 - sterm4)
    dd^v * Sum(term1 * term2, list(n, 0, S_Infinity))
  })

#' @rdname joint_rv
#' @export
GeneralizedMultivariateLogGamma <- function(syms, delta, v, lamda, mu)
  .mvrv("GeneralizedMultivariateLogGammaDistribution", syms, list(delta, v, lamda, mu))

#' @rdname joint_rv
#' @export
GeneralizedMultivariateLogGammaOmega <- function(syms, omega, v, lamda, mu) {
  om <- .mx(.jarg(omega))
  ismat <- inherits(om, "symmatrix")
  .value_check(list(ismat && om$nrow == om$ncol, ismat), "omega must be a square matrix")
  for (val in om$data) .value_check(list(val >= 0, val <= 1), "all values in matrix must be between 0 and 1(both inclusive).")
  dg <- lapply(seq_len(om$nrow), function(i) .mget(om, i, i))
  .value_check(all(vapply(dg, function(z) identical(z$key, S_One$key), NA)), "all the elements of diagonal should be 1.")
  nl <- length(.mentries(.jarg(lamda))); nm <- length(.mentries(.jarg(mu)))
  .value_check(list(om$nrow == nl, nl == nm), "lamda, mu should be of same length and omega should  be of shape (length of lamda, length of mu)")
  .value_check(nl > 1L, "the distribution should have at least two random variables.")
  delta <- as_sym(det(om))^Rational(1, nl - 1L)
  GeneralizedMultivariateLogGamma(syms, delta, v, lamda, mu)
}

.defjd("MultinomialDistribution", c("n", "p"), continuous = FALSE,
  check = function(n, p) {
    .value_check(n > 0, "number of trials must be a positive integer")
    pe <- .mentries(p)
    for (p_k in pe) .value_check(list(p_k >= 0, p_k <= 1), "probability must be in range [0, 1]")
    .value_check(Eq(.add(c(list(S_Zero), pe)), 1), "probabilities must sum to 1")
  },
  set = function(d) .rv_set_pow(Intersection(Naturals0, Interval(0, d$args[[1L]])), length(.mentries(d$args[[2L]]))),
  pdf = function(d, ...) {
    x <- lapply(list(...), as_sym)
    n <- d$args[[1L]]; p <- .mentries(d$args[[2L]])
    term_1 <- sym_factorial(n) / .mul(lapply(x, sym_factorial))
    term_2 <- .mul(lapply(seq_along(p), function(i) p[[i]]^x[[i]]))
    Piecewise(list(term_1 * term_2, Eq(.add(c(list(S_Zero), x)), n)), list(S_Zero, S_true))
  })

#' @rdname joint_rv
#' @export
Multinomial <- function(syms, n, ...) {
  p <- list(...)
  if (!(is.list(p[[1L]]) && !is_sym(p[[1L]])) && !inherits(p[[1L]], "symmatrix")) p <- list(p)
  .mvrv("MultinomialDistribution", syms, list(n, p[[1L]]))
}

.defjd("NegativeMultinomialDistribution", c("k0", "p"), continuous = FALSE,
  check = function(k0, p) {
    .value_check(k0 > 0, "number of failures must be a positive integer")
    pe <- .mentries(p)
    for (p_k in pe) .value_check(list(p_k >= 0, p_k <= 1), "probability must be in range [0, 1].")
    .value_check(.add(c(list(S_Zero), pe)) <= 1, "success probabilities must not be greater than 1.")
  },
  set = function(d) .rv_set_pow(Range(0, S_Infinity), length(.mentries(d$args[[2L]]))),
  pdf = function(d, ...) {
    k <- lapply(list(...), as_sym)
    k0 <- d$args[[1L]]; p <- .mentries(d$args[[2L]])
    term_1 <- gamma(k0 + .add(c(list(S_Zero), k))) * (1 - .add(c(list(S_Zero), p)))^k0 / gamma(k0)
    term_2 <- .mul(lapply(seq_along(p), function(i) p[[i]]^k[[i]] / sym_factorial(k[[i]])))
    term_1 * term_2
  })

#' @rdname joint_rv
#' @export
NegativeMultinomial <- function(syms, k0, ...) {
  p <- list(...)
  if (!(is.list(p[[1L]]) && !is_sym(p[[1L]])) && !inherits(p[[1L]], "symmatrix")) p <- list(p)
  .mvrv("NegativeMultinomialDistribution", syms, list(k0, p[[1L]]))
}
