# Compound distributions (sympy/stats/compound_rv.py): a distribution whose
# parameters are themselves random symbols.

.defclass("CompoundDistribution", "Distribution")
.defclass("CompoundPSpace", "PSpace")

.CompoundDistribution <- function(dist) {
  if (!(.isa(dist, "ContinuousDistribution") || .isa(dist, "SingleFiniteDistribution") || .isa(dist, "DiscreteDistribution")))
    .serr("NotImplementedError", paste0("Compound Distribution for ", .str(dist), " is not implemented yet"))
  randoms <- unlist(lapply(dist$args, function(a) if (is_sym(a)) random_symbols(a) else list()), recursive = FALSE)
  if (!length(randoms)) return(dist)
  new_sym("CompoundDistribution", list(dist))
}

.CompoundPSpace <- function(s, distribution) {
  s <- .symbol_converter(s)
  if (.isa(distribution, "ContinuousDistribution")) return(.SingleContinuousPSpace(s, distribution))
  if (.isa(distribution, "DiscreteDistribution")) return(.SingleDiscretePSpace(s, distribution))
  if (.isa(distribution, "SingleFiniteDistribution")) return(new_sym("SingleFinitePSpace", list(s, distribution)))
  if (!.isa(distribution, "CompoundDistribution"))
    .serr("ValueError", paste0(.str(distribution), " should be an isinstance of CompoundDistribution"))
  new_sym("CompoundPSpace", list(s, distribution))
}

.cd_parent <- function(d) d$args[[1L]]
.defm("CompoundDistribution", "set", function(d) {
  s <- .dist_set(.cd_parent(d))
  if (is.list(s) && !is_sym(s)) do.call(FiniteSet, s) else s
})
.defm("CompoundDistribution", "is_Continuous", function(d) .isa(.cd_parent(d), "ContinuousDistribution"))
.defm("CompoundDistribution", "is_Finite", function(d) .isa(.cd_parent(d), "SingleFiniteDistribution"))
.defm("CompoundDistribution", "is_Discrete", function(d) .isa(.cd_parent(d), "DiscreteDistribution"))
.defm("CompoundDistribution", "pdf", function(d, x, evaluate = FALSE) {
  dist <- .cd_parent(d)
  randoms <- Filter(is_random, dist$args)
  if (.isa(dist, "SingleFiniteDistribution")) {
    y <- Dummy("y", integer = TRUE, negative = FALSE)
    expr <- .callm(dist, "pmf", y)
  } else {
    y <- Dummy("y")
    expr <- .dist_pdf(dist, y)
  }
  for (rv in randoms) expr <- .cd_marginalise(expr, rv, evaluate)
  apply_lambda(Lambda(y, expr), x)
})

.cd_marginalise <- function(expr, rv, evaluate) {
  ps <- .rs_pspace(rv)
  rd <- .callm(ps, "distribution")
  rv_dens <- if (.isa(rd, "SingleFiniteDistribution")) .callm(rd, "pmf", rv) else .dist_pdf(rd, rv)
  rv_dom <- .dom_set(.ps_domain(ps))
  lims <- list(rv, set_inf(rv_dom), set_sup(rv_dom))
  expr <- if (.ps_disc(ps) || .ps_fin(ps)) Sum(expr * rv_dens, lims) else Integral(expr * rv_dens, lims)
  if (evaluate) doit(expr) else expr
}

.defm("CompoundPSpace", "value", function(ps) .RandomSymbol(ps$args[[1L]], ps))
.defm("CompoundPSpace", "symbol", function(ps) ps$args[[1L]])
.defm("CompoundPSpace", "distribution", function(ps) ps$args[[2L]])
.defm("CompoundPSpace", "is_Continuous", function(ps) .callm(ps$args[[2L]], "is_Continuous"))
.defm("CompoundPSpace", "is_Finite", function(ps) .callm(ps$args[[2L]], "is_Finite"))
.defm("CompoundPSpace", "is_Discrete", function(ps) .callm(ps$args[[2L]], "is_Discrete"))
.defm("CompoundPSpace", "pdf", function(ps) .callm(ps$args[[2L]], "pdf", ps$args[[1L]]))
.defm("CompoundPSpace", "set", function(ps) .dist_set(ps$args[[2L]]))
.defm("CompoundPSpace", "domain", function(ps) .ps_domain(.cps_new(ps)))
.defm("CompoundPSpace", "symbols", function(ps) list(ps$args[[1L]]))

.cps_new <- function(ps, evaluate = FALSE) {
  x <- Dummy("x")
  parent <- .cd_parent(ps$args[[2L]])
  func <- Lambda(x, .callm(ps$args[[2L]], "pdf", x, evaluate = evaluate))
  sym <- ps$args[[1L]]
  pdf <- Lambda(sym, apply_lambda(func, sym))
  set <- .dist_set(parent)
  if (.isa(parent, "ContinuousDistribution")) return(.SingleContinuousPSpace(sym, .ContinuousDistributionHandmade(pdf, set)))
  if (.isa(parent, "DiscreteDistribution")) return(.SingleDiscretePSpace(sym, .DiscreteDistributionHandmade(pdf, set)))
  if (.isa(parent, "SingleFiniteDistribution")) {
    ks <- if (is.list(set) && !is_sym(set)) set else .finite_elems(set)
    dens <- .FD(ks, lapply(ks, function(k) apply_lambda(pdf, k)))
    return(new_sym("SingleFinitePSpace", list(sym, new_sym("FiniteDistributionHandmade", list(dens)))))
  }
  .serr("NotImplementedError", paste0("Compound Distribution for ", .str(parent), " is not implemented yet"))
}

.cps_swap <- function(ps, newps, expr) subs(expr, .ps_value(ps), .ps_value(newps))

.defm("CompoundPSpace", "compute_density", function(ps, expr, compound_evaluate = TRUE, ...) {
  np <- .cps_new(ps, compound_evaluate)
  .callm(np, "compute_density", .cps_swap(ps, np, expr), ...)
})
.defm("CompoundPSpace", "compute_cdf", function(ps, expr, compound_evaluate = TRUE, ...) {
  np <- .cps_new(ps, compound_evaluate)
  .callm(np, "compute_cdf", .cps_swap(ps, np, expr), ...)
})
.defm("CompoundPSpace", "compute_expectation", function(ps, expr, rvs = NULL, evaluate = FALSE, ...) {
  np <- .cps_new(ps, evaluate)
  expr <- .cps_swap(ps, np, expr)
  if (!is.null(rvs) && is_sym(rvs)) rvs <- .cps_swap(ps, np, rvs)
  if (.isa(np, "SingleFinitePSpace")) return(.callm(np, "compute_expectation", expr, rvs, ...))
  .callm(np, "compute_expectation", expr, rvs, evaluate = evaluate, ...)
})
.defm("CompoundPSpace", "probability", function(ps, condition, compound_evaluate = TRUE, ...) {
  np <- .cps_new(ps, compound_evaluate)
  .callm(np, "probability", .cps_swap(ps, np, condition))
})
.defm("CompoundPSpace", "conditional_space", function(ps, condition, compound_evaluate = TRUE, ...) {
  np <- .cps_new(ps, compound_evaluate)
  .callm(np, "conditional_space", .cps_swap(ps, np, condition))
})

# A distribution object from its class name and parameters (SymPy's
# XxxDistribution(...) constructors), used for the compound distributions.
.mkdist <- function(cls, ...) {
  args <- list(...)
  if (cls %in% c("MultivariateNormalDistribution", "MultivariateLaplaceDistribution", "MultivariateTDistribution",
                 "NormalGammaDistribution", "MultivariateBetaDistribution", "MultivariateEwensDistribution",
                 "GeneralizedMultivariateLogGammaDistribution", "MultinomialDistribution", "NegativeMultinomialDistribution"))
    return(.jd_new(cls, args))
  if (cls %in% c("MatrixGammaDistribution", "WishartDistribution", "MatrixNormalDistribution", "MatrixStudentTDistribution"))
    return(new_sym(cls, lapply(args, .mat_list)))
  new_sym(cls, lapply(args, as_sym))
}
