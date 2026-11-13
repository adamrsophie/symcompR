# Discrete (countably infinite) probability spaces, after
# sympy/stats/drv.py and drv_types.py.

.defclass("DiscreteDistribution", "Distribution")
.defclass("SingleDiscreteDistribution", "DiscreteDistribution")
.defclass("DiscreteDomain", "RandomDomain")
.defclass("SingleDiscreteDomain", c("DiscreteDomain", "SingleDomain"))
.defclass("ConditionalDiscreteDomain", c("DiscreteDomain", "ConditionalDomain"))
.defclass("ProductDiscreteDomain", c("ProductDomain", "DiscreteDomain"))
.defclass("DiscretePSpace", "PSpace")
.defclass("SingleDiscretePSpace", c("DiscretePSpace", "SinglePSpace"))

.defm("DiscreteDomain", "is_Discrete", function(d) TRUE)
.defm("DiscretePSpace", "is_Discrete", function(ps) TRUE)
.defm("DiscretePSpace", "is_real", function(ps) TRUE)

.defm("SingleDiscreteDistribution", "set", function(d) Integers)

.defm("SingleDiscreteDistribution", "compute_cdf", function(d, ...) {
  x <- Dummy("x", integer = TRUE); z <- Dummy("z", real = TRUE)
  left <- set_inf(.dist_set(d))
  pdf <- .dist_pdf(d, x)
  cdf <- summation(pdf, list(x, left, floor(z)))
  cdf <- Piecewise(list(cdf, Ge(z, left)), list(S_Zero, S_true))
  Lambda(z, cdf)
})
.defm("SingleDiscreteDistribution", "cdf", function(d, x, ...) .dist_cdf(d, x, ...))

.defm("SingleDiscreteDistribution", "compute_characteristic_function", function(d, ...) {
  x <- Dummy("x", real = TRUE); t <- Dummy("t", real = TRUE)
  pdf <- .dist_pdf(d, x)
  s <- .dist_set(d)
  Lambda(t, summation(.mul(list(exp(.mul(list(S_I, t, x))), pdf)), list(x, set_inf(s), set_sup(s))))
})
.defm("SingleDiscreteDistribution", "characteristic_function", function(d, t, ...) {
  if (!...length()) {
    c <- .opt_method(d, "_characteristic_function", t)
    if (!is.null(c)) return(c)
  }
  apply_lambda(.callm(d, "compute_characteristic_function", ...), t)
})

.defm("SingleDiscreteDistribution", "compute_moment_generating_function", function(d, ...) {
  t <- Dummy("t", real = TRUE); x <- Dummy("x", integer = TRUE)
  pdf <- .dist_pdf(d, x)
  s <- .dist_set(d)
  Lambda(t, summation(.mul(list(exp(.mul(list(t, x))), pdf)), list(x, set_inf(s), set_sup(s))))
})
.defm("SingleDiscreteDistribution", "moment_generating_function", function(d, t, ...) {
  if (!...length()) {
    c <- .opt_method(d, "_moment_generating_function", t)
    if (!is.null(c)) return(c)
  }
  apply_lambda(.callm(d, "compute_moment_generating_function", ...), t)
})

.defm("SingleDiscreteDistribution", "compute_quantile", function(d, ...) {
  x <- Dummy("x", integer = TRUE); p <- Dummy("p", real = TRUE)
  left <- set_inf(.dist_set(d))
  pdf <- .dist_pdf(d, x)
  cdf <- summation(pdf, list(x, left, x))
  Lambda(p, Piecewise(list(x, Le(p, cdf))))
})
.defm("SingleDiscreteDistribution", "quantile", function(d, x, ...) {
  if (!...length()) {
    q <- .opt_method(d, "_quantile", x)
    if (!is.null(q)) return(q)
  }
  apply_lambda(.callm(d, "compute_quantile", ...), x)
})

.defm("SingleDiscreteDistribution", "expectation", function(d, expr, var, evaluate = TRUE, ...) {
  s <- .dist_set(d)
  if (evaluate) {
    p <- .poly_or_null(expr, var)
    if (is.null(p)) return(summation(.mul(list(expr, .dist_pdf(d, var))), list(var, set_inf(s), set_sup(s))))
    t <- Dummy("t", real = TRUE)
    mgf <- .callm(d, "moment_generating_function", t)
    cs <- .poly_coeffs_by_degree(p)
    deg <- length(cs) - 1L
    taylor <- removeO(series(mgf, t, 0, deg + 1L))
    tp <- .poly_or_null(taylor, t)
    if (is.null(tp)) return(summation(.mul(list(expr, .dist_pdf(d, var))), list(var, set_inf(s), set_sup(s))))
    tc <- .poly_coeffs_by_degree(tp)
    result <- S_Zero
    for (k in 0:deg) {
      tk <- if (k + 1L <= length(tc)) tc[[k + 1L]] else S_Zero
      result <- .add(list(result, .mul(list(cs[[k + 1L]], tk, sym_factorial(as_sym(k))))))
    }
    return(result)
  }
  Sum(.mul(list(expr, .dist_pdf(d, var))), list(var, set_inf(s), set_sup(s)))
})
.defm("DiscreteDistribution", "call", function(d, ...) .callm(d, "pdf", ...))

# ---------------------------------------------------------------------------
# Domains

.defm("SingleDiscreteDomain", "as_boolean", function(d) Contains(d$args[[1L]], d$args[[2L]]))
.ConditionalDiscreteDomain <- function(fulldomain, condition) .ConditionalDomain_new("ConditionalDiscreteDomain", fulldomain, condition)
.defm("ConditionalDiscreteDomain", "set", function(d) {
  syms <- .dom_symbols(d)
  if (length(syms) > 1L) .serr("NotImplementedError", "\nMultivariate conditional domains are not yet implemented.")
  Intersection(.reduce_wrap(d$args[[2L]], syms[[1L]]), .dom_set(d$args[[1L]]))
})
.defm("ProductDiscreteDomain", "as_boolean", function(d) .serr("TypeError", "can't convert method to Boolean"))

# ---------------------------------------------------------------------------
# DiscretePSpace

.defm("DiscretePSpace", "pdf", function(ps) do.call(.rv_lcall, c(list(.callm(ps, "density")), .ps_symbols(ps))))

.defm("DiscretePSpace", "where", function(ps, condition) {
  rvs <- random_symbols(condition)
  if (length(rvs) > 1L) .serr("NotImplementedError", "Multivariate discrete\nrandom variables are not yet supported.")
  cd <- .reduce_wrap(condition, rvs[[1L]])
  cd <- Intersection(cd, .dom_set(.ps_domain(ps)))
  .domain_new("SingleDiscreteDomain", .rs_symbol(rvs[[1L]]), cd)
})

.defm("DiscretePSpace", "probability", function(ps, condition, ...) {
  complement <- condition$type == "Unequality"
  if (complement) condition <- Eq(condition$args[[1L]], condition$args[[2L]])
  prob <- tryCatch({
    dom <- .dom_set(.callm(ps, "where", condition))
    if (.is_false(condition) || dom$type == "EmptySet") return(S_Zero)
    if (.is_true(condition) || identical(dom$key, .dom_set(.ps_domain(ps))$key)) return(S_One)
    .dps_eval_prob(ps, dom)
  }, NotImplementedError = function(e) {
    expr <- .sub(condition$args[[1L]], condition$args[[2L]])
    dens <- rv_density(expr)
    if (!.isa(dens, "DiscreteDistribution")) dens <- .DiscreteDistributionHandmade(dens)
    z <- Dummy("z", real = TRUE)
    space <- .SingleDiscretePSpace(z, dens)
    .callm(space, "probability", .rel_same(condition, .ps_value(space), S_Zero))
  })
  if (is.null(prob)) prob <- Probability(condition)
  if (complement) .sub(S_One, prob) else prob
})

.dps_eval_prob <- function(ps, dom) {
  sym <- .ps_symbols(ps)[[1L]]
  pdf <- .ps_pdf(ps)
  if (dom$type == "Range") {
    n <- Symbol("n", integer = TRUE)
    a <- dom$args
    inf <- a[[1L]]; sup <- a[[2L]]; step <- a[[3L]]
    summand <- subs(pdf, sym, n * step)
    return(doit(summation(summand, list(n, inf / step, sup / step - 1))))
  }
  if (dom$type == "FiniteSet") {
    f <- Lambda(sym, pdf)
    return(.add(c(list(S_Zero), lapply(dom$args, function(x) apply_lambda(f, x)))))
  }
  if (dom$type == "Union") return(.add(c(list(S_Zero), lapply(dom$args, function(x) .dps_eval_prob(ps, x)))))
  NULL
}

.defm("DiscretePSpace", "conditional_space", function(ps, condition, ...) {
  density <- Lambda(.ps_symbols(ps), .div(.ps_pdf(ps), .callm(ps, "probability", condition)))
  condition <- xreplace(.as_bool_sym(condition), lapply(.ps_values(ps), function(rv) list(rv, .rs_symbol(rv))))
  domain <- .ConditionalDiscreteDomain(.ps_domain(ps), condition)
  new_sym("DiscretePSpace", list(domain, density))
})

.SingleDiscretePSpace <- function(s, distribution) new_sym("SingleDiscretePSpace", list(.symbol_converter(s), distribution))

.defm("SingleDiscretePSpace", "set", function(ps) .dist_set(ps$args[[2L]]))
.defm("SingleDiscretePSpace", "domain", function(ps) .domain_new("SingleDiscreteDomain", ps$args[[1L]], .dist_set(ps$args[[2L]])))
.defm("SingleDiscretePSpace", "density", function(ps) ps$args[[2L]])
.defm("SingleDiscretePSpace", "pdf", function(ps) .dist_pdf(ps$args[[2L]], ps$args[[1L]]))

.defm("SingleDiscretePSpace", "compute_expectation", function(ps, expr, rvs = NULL, evaluate = TRUE, ...) {
  val <- .ps_value(ps)
  if (is.null(rvs) || !length(rvs)) rvs <- list(val)
  if (!.in_list(val, rvs)) return(expr)
  expr <- xreplace(as_sym(expr), lapply(rvs, function(rv) list(rv, .rs_symbol(rv))))
  x <- .rs_symbol(val)
  tryCatch(.callm(ps$args[[2L]], "expectation", expr, x, evaluate = evaluate, ...),
           NotImplementedError = function(e) {
             s <- .callm(ps, "set")
             Sum(.mul(list(expr, .ps_pdf(ps))), list(x, set_inf(s), set_sup(s)))
           })
})
.defm("SingleDiscretePSpace", "compute_cdf", function(ps, expr, ...) {
  if (.eq_struct(expr, .ps_value(ps))) {
    x <- Dummy("x", real = TRUE)
    return(Lambda(x, .dist_cdf(ps$args[[2L]], x, ...)))
  }
  .serr("NotImplementedError", "")
})
.defm("SingleDiscretePSpace", "compute_density", function(ps, expr, ...) {
  if (.eq_struct(expr, .ps_value(ps))) return(ps$args[[2L]])
  .serr("NotImplementedError", "")
})
.defm("SingleDiscretePSpace", "compute_characteristic_function", function(ps, expr, ...) {
  if (.eq_struct(expr, .ps_value(ps))) {
    t <- Dummy("t", real = TRUE)
    return(Lambda(t, .callm(ps$args[[2L]], "characteristic_function", t, ...)))
  }
  .serr("NotImplementedError", "")
})
.defm("SingleDiscretePSpace", "compute_moment_generating_function", function(ps, expr, ...) {
  if (.eq_struct(expr, .ps_value(ps))) {
    t <- Dummy("t", real = TRUE)
    return(Lambda(t, .callm(ps$args[[2L]], "moment_generating_function", t, ...)))
  }
  .serr("NotImplementedError", "")
})
.defm("SingleDiscretePSpace", "compute_quantile", function(ps, expr, ...) {
  if (.eq_struct(expr, .ps_value(ps))) {
    p <- Dummy("p", real = TRUE)
    return(Lambda(p, .callm(ps$args[[2L]], "quantile", p, ...)))
  }
  .serr("NotImplementedError", "")
})

# ---------------------------------------------------------------------------
# Discrete distribution types

.defdd <- function(.cls, .argnames, ...) .defdist(.cls, "SingleDiscreteDistribution", .argnames, ...)
.rvd <- function(name, cls, args, check = TRUE) .rv_make(name, cls, args, check = check, single = "SingleDiscretePSpace")

#' Discrete random variables
#'
#' Constructors for the countably infinite discrete distributions of
#' `sympy.stats`: `DiscreteRV()` from a probability mass function,
#' `FlorySchulz()`, `Geometric()`, `Hermite()`, `Logarithmic()`,
#' `NegativeBinomial()`, `Poisson()`, `Skellam()`, `YuleSimon()` and
#' `Zeta()`.
#'
#' @param name The name of the random variable.
#' @param symbol,density,set For `DiscreteRV()`: the symbol, the mass
#'   function in it, and its support (the integers by default).
#' @param check Whether to check the arguments.
#' @param a,p,a1,a2,r,lamda,mu1,mu2,rho,s Distribution parameters, named as
#'   in SymPy.
#' @return A random symbol.
#' @examples
#' X <- Poisson("X", 3)
#' E(X)
#' variance(Geometric("G", Rational(1, 4)))
#' @name discrete_rv
NULL

.defdd("DiscreteDistributionHandmade", c("pdf"),
  set = function(d) d$args[[2L]],
  pdf = function(d, x) .rv_lcall(d$args[[1L]], x),
  check = function(pdf, set) {
    x <- Dummy("x")
    val <- doit(Sum(.rv_lcall(pdf, x), list(x, set_inf(set), set_sup(set))))
    .value_check(!identical(.fuzzy_bool(Eq(val, 1)), FALSE), "The pdf is incorrect on the given set.")
  })

.DiscreteDistributionHandmade <- function(pdf, set = Integers) new_sym("DiscreteDistributionHandmade", list(pdf, set))

#' @rdname discrete_rv
#' @export
DiscreteRV <- function(symbol, density, set = Integers, check = FALSE) {
  set <- as_sym(set)
  pdf <- Piecewise(list(density, as_relational(set, symbol)), list(0, S_true))
  pdf <- Lambda(symbol, pdf)
  .rvd(symbol$name, "DiscreteDistributionHandmade", list(pdf, set), check = check)
}

.defdd("FlorySchulzDistribution", c("a"),
  set = function(d) Naturals,
  check = function(a) .value_check(list(0 < a, a < 1), "a must be between 0 and 1"),
  pdf = function(d, k) { a <- .da(d, "a"); a^2 * k * (1 - a)^(k - 1) },
  `_characteristic_function` = function(d, t) { a <- .da(d, "a"); a^2 * exp(S_I * t) / (1 + (a - 1) * exp(S_I * t))^2 },
  `_moment_generating_function` = function(d, t) { a <- .da(d, "a"); a^2 * exp(t) / (1 + (a - 1) * exp(t))^2 })

#' @rdname discrete_rv
#' @export
FlorySchulz <- function(name, a) .rvd(name, "FlorySchulzDistribution", list(a))

.defdd("GeometricDistribution", c("p"),
  set = function(d) Naturals,
  check = function(p) .value_check(list(0 < p, p <= 1), "p must be between 0 and 1"),
  pdf = function(d, k) { p <- .da(d, "p"); (1 - p)^(k - 1) * p },
  `_characteristic_function` = function(d, t) { p <- .da(d, "p"); p * exp(S_I * t) / (1 - (1 - p) * exp(S_I * t)) },
  `_moment_generating_function` = function(d, t) { p <- .da(d, "p"); p * exp(t) / (1 - (1 - p) * exp(t)) })

#' @rdname discrete_rv
#' @export
Geometric <- function(name, p) .rvd(name, "GeometricDistribution", list(p))

.defdd("HermiteDistribution", c("a1", "a2"),
  set = function(d) Naturals0,
  check = function(a1, a2) {
    .value_check(.st_fact(a1, "nonnegative"), "Parameter a1 must be >= 0.")
    .value_check(.st_fact(a2, "nonnegative"), "Parameter a2 must be >= 0.")
  },
  pdf = function(d, k) {
    p <- .dl(d)
    term1 <- exp(-(p$a1 + p$a2))
    j <- Dummy("j", integer = TRUE)
    num <- p$a1^(k - 2 * j) * p$a2^j
    den <- sym_factorial(k - 2 * j) * sym_factorial(j)
    term1 * doit(Sum(num / den, list(j, 0, floor(as_sym(k) / 2))))
  },
  `_moment_generating_function` = function(d, t) { p <- .dl(d); exp(p$a1 * (exp(t) - 1) + p$a2 * (exp(2 * t) - 1)) },
  `_characteristic_function` = function(d, t) { p <- .dl(d); exp(p$a1 * (exp(S_I * t) - 1) + p$a2 * (exp(2 * S_I * t) - 1)) })

#' @rdname discrete_rv
#' @export
Hermite <- function(name, a1, a2) .rvd(name, "HermiteDistribution", list(a1, a2))

.defdd("LogarithmicDistribution", c("p"),
  set = function(d) Naturals,
  check = function(p) .value_check(list(p > 0, p < 1), "p should be between 0 and 1"),
  pdf = function(d, k) { p <- .da(d, "p"); -1 * p^k / (k * log(1 - p)) },
  `_characteristic_function` = function(d, t) { p <- .da(d, "p"); log(1 - p * exp(S_I * t)) / log(1 - p) },
  `_moment_generating_function` = function(d, t) { p <- .da(d, "p"); log(1 - p * exp(t)) / log(1 - p) })

#' @rdname discrete_rv
#' @export
Logarithmic <- function(name, p) .rvd(name, "LogarithmicDistribution", list(p))

.defdd("NegativeBinomialDistribution", c("r", "p"),
  set = function(d) Naturals0,
  check = function(r, p) {
    .value_check(r > 0, "r should be positive")
    .value_check(list(p > 0, p < 1), "p should be between 0 and 1")
  },
  pdf = function(d, k) { q <- .dl(d); sym_binomial(k + q$r - 1, k) * (1 - q$p)^k * q$p^q$r },
  `_characteristic_function` = function(d, t) { q <- .dl(d); (q$p / (1 - (1 - q$p) * exp(S_I * t)))^q$r },
  `_moment_generating_function` = function(d, t) { q <- .dl(d); (q$p / (1 - (1 - q$p) * exp(t)))^q$r })

#' @rdname discrete_rv
#' @export
NegativeBinomial <- function(name, r, p) .rvd(name, "NegativeBinomialDistribution", list(r, p))

.defdd("PoissonDistribution", c("lamda"),
  set = function(d) Naturals0,
  check = function(lamda) .value_check(lamda > 0, "Lambda must be positive"),
  pdf = function(d, k) { l <- .da(d, "lamda"); l^k / sym_factorial(k) * exp(-l) },
  `_characteristic_function` = function(d, t) exp(.da(d, "lamda") * (exp(S_I * t) - 1)),
  `_moment_generating_function` = function(d, t) exp(.da(d, "lamda") * (exp(t) - 1)),
  expectation = function(d, expr, var, evaluate = TRUE, ...) {
    if (evaluate) {
      if (.eq_struct(expr, var)) return(.da(d, "lamda"))
      if (expr$type == "Function" && expr$name == "FallingFactorial" && isTRUE(.st_fact(expr$args[[2L]], "integer")) &&
          isTRUE(.st_fact(expr$args[[2L]], "positive")) && .eq_struct(expr$args[[1L]], var))
        return(.da(d, "lamda")^expr$args[[2L]])
    }
    .STM[["SingleDiscreteDistribution"]]$expectation(d, expr, var, evaluate, ...)
  })

#' @rdname discrete_rv
#' @export
Poisson <- function(name, lamda) .rvd(name, "PoissonDistribution", list(lamda))

.defdd("SkellamDistribution", c("mu1", "mu2"),
  set = function(d) Integers,
  check = function(mu1, mu2) {
    .value_check(mu1 >= 0, "Parameter mu1 must be >= 0")
    .value_check(mu2 >= 0, "Parameter mu2 must be >= 0")
  },
  pdf = function(d, k) {
    p <- .dl(d)
    exp(-(p$mu1 + p$mu2)) * (p$mu1 / p$mu2)^(k / 2) * besseli(k, 2 * sqrt(p$mu1 * p$mu2))
  },
  `_cdf` = function(d, x) .serr("NotImplementedError", "Skellam doesn't have closed form for the CDF."),
  `_characteristic_function` = function(d, t) { p <- .dl(d); exp(-(p$mu1 + p$mu2) + p$mu1 * exp(S_I * t) + p$mu2 * exp(-S_I * t)) },
  `_moment_generating_function` = function(d, t) { p <- .dl(d); exp(-(p$mu1 + p$mu2) + p$mu1 * exp(t) + p$mu2 * exp(-t)) })

#' @rdname discrete_rv
#' @export
Skellam <- function(name, mu1, mu2) .rvd(name, "SkellamDistribution", list(mu1, mu2))

.defdd("YuleSimonDistribution", c("rho"),
  set = function(d) Naturals,
  check = function(rho) .value_check(rho > 0, "rho should be positive"),
  pdf = function(d, k) { rho <- .da(d, "rho"); rho * sym_beta(k, rho + 1) },
  `_cdf` = function(d, x) { rho <- .da(d, "rho"); Piecewise(list(1 - floor(x) * sym_beta(floor(x), rho + 1), x >= 1), list(0, S_true)) },
  `_characteristic_function` = function(d, t) { rho <- .da(d, "rho"); rho * hyper(list(1, 1), list(rho + 2), exp(S_I * t)) * exp(S_I * t) / (rho + 1) },
  `_moment_generating_function` = function(d, t) { rho <- .da(d, "rho"); rho * hyper(list(1, 1), list(rho + 2), exp(t)) * exp(t) / (rho + 1) })

#' @rdname discrete_rv
#' @export
YuleSimon <- function(name, rho) .rvd(name, "YuleSimonDistribution", list(rho))

.defdd("ZetaDistribution", c("s"),
  set = function(d) Naturals,
  check = function(s) .value_check(s > 1, "s should be greater than 1"),
  pdf = function(d, k) { s <- .da(d, "s"); 1 / (k^s * zeta(s)) },
  `_characteristic_function` = function(d, t) { s <- .da(d, "s"); polylog(s, exp(S_I * t)) / zeta(s) },
  `_moment_generating_function` = function(d, t) { s <- .da(d, "s"); polylog(s, exp(t)) / zeta(s) })

#' @rdname discrete_rv
#' @export
Zeta <- function(name, s) .rvd(name, "ZetaDistribution", list(s))
