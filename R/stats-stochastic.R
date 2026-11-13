# Stochastic processes (sympy/stats/stochastic_process.py and
# sympy/stats/stochastic_process_types.py): the Bernoulli process and the
# counting processes (Poisson, Wiener, gamma). A process is an R function,
# so that X(t) indexes a continuous time process; X[t] is rv_getitem(X, t).

.defclass("StochasticProcess")
.defclass("DiscreteTimeStochasticProcess", "StochasticProcess")
.defclass("ContinuousTimeStochasticProcess", "StochasticProcess")
.defclass("BernoulliProcess", "DiscreteTimeStochasticProcess")
.defclass("CountingProcess", "ContinuousTimeStochasticProcess")
.defclass("PoissonProcess", "CountingProcess")
.defclass("WienerProcess", "CountingProcess")
.defclass("GammaProcess", "CountingProcess")
.defclass("StochasticPSpace", "ProductPSpace")
.defclass("RandomIndexedSymbol")

.mkproc <- function(node) {
  f <- function(time) .proc_call(node, time)
  structure(f, class = c("rv_process", "function"), node = node)
}
.pnode <- function(x) if (inherits(x, "rv_process")) attr(x, "node") else x

#' @export
print.rv_process <- function(x, ...) {
  cat(.str(attr(x, "node")), "\n")
  invisible(x)
}

.process_pspace <- function(x) .serr("ValueError", "Expression containing Random Variable expected, not a process")
.process_attr <- function(x, name, ...) {
  node <- attr(x, "node")
  f <- .findm(node, name)
  if (!is.null(f)) return(f(node, ...))
  .serr("AttributeError", paste0("'", node$type, "' object has no attribute '", name, "'"))
}
.process_getitem <- function(x, i) .callm(attr(x, "node"), "getitem", i)
.proc_call <- function(node, time) .callm(node, "call", time)

# RandomIndexedSymbol(idx_obj, pspace): a node whose arguments are the
# indexed object and the space, so substitution reaches the index
.RandomIndexedSymbol <- function(idx_obj, pspace = NULL) {
  if (is.null(pspace)) pspace <- .PSpace0()
  if (!(is_sym(idx_obj) && idx_obj$type %in% c("Indexed", "Function")))
    .serr("TypeError", paste0("An Function or Indexed object is expected not ", .str(idx_obj)))
  new_sym("RandomIndexedSymbol", list(idx_obj, pspace),
          extra = list(rv_class = "RandomIndexedSymbol", sympy_class = "RandomIndexedSymbol"))
}
register_str_printer("RandomIndexedSymbol", function(e) .str(e$args[[1L]]))
.facts_hooks$RandomIndexedSymbol <- function(e) c(commutative = TRUE, finite = TRUE)
# the space is not part of the free symbols; the symbol itself is
.bound_hooks$RandomIndexedSymbol <- function(e) setdiff(names(.free_sym_list(e$args[[2L]])), names(.free_sym_list(e$args[[1L]])))
.ris_key <- function(x) {
  s <- x$args[[1L]]
  if (s$type == "Indexed") s$args[[2L]] else s$args[[1L]]
}

#' @rdname stochastic_process
#' @param idx_obj An indexed object or an applied function.
#' @export
RandomIndexedSymbol <- function(idx_obj, pspace = NULL) .RandomIndexedSymbol(as_sym(idx_obj), pspace)

.StochasticPSpace <- function(sym, process, distribution = NULL) {
  if (is.null(distribution)) distribution <- new_sym("Distribution", list())
  new_sym("StochasticPSpace", list(.symbol_converter(sym), process, distribution))
}
.defm("StochasticPSpace", "process", function(ps) ps$args[[2L]])
.defm("StochasticPSpace", "symbol", function(ps) ps$args[[1L]])
.defm("StochasticPSpace", "distribution", function(ps) ps$args[[3L]])
.defm("StochasticPSpace", "domain", function(ps) .ProductDomain(.callm(ps$args[[2L]], "index_set"), .callm(ps$args[[2L]], "state_space")))
.defm("StochasticPSpace", "probability", function(ps, condition, given_condition = NULL, evaluate = TRUE, ...)
  .callm(ps$args[[2L]], "probability", condition, given_condition, evaluate = evaluate))
.defm("StochasticPSpace", "compute_expectation", function(ps, expr, condition = NULL, evaluate = TRUE, ...)
  .callm(ps$args[[2L]], "expectation", expr, condition, evaluate = evaluate))

# ---------------------------------------------------------------------------
# StochasticProcess

.defm("StochasticProcess", "symbol", function(p) p$args[[1L]])
.defm("StochasticProcess", "index_set", function(p) Reals)
.defm("StochasticProcess", "distribution", function(p, key = NULL) new_sym("Distribution", list()))
.defm("StochasticProcess", "call", function(p, time) .serr("NotImplementedError", "Use [] for indexing discrete time stochastic process."))
.defm("StochasticProcess", "getitem", function(p, time) .serr("NotImplementedError", "Use () for indexing continuous time stochastic process."))
.defm("StochasticProcess", "probability", function(p, condition, ...) .serr("NotImplementedError", ""))
.defm("StochasticProcess", "expectation", function(p, ...) .serr("NotImplementedError", "Abstract method for expectation queries."))
.defm("StochasticProcess", "joint_distribution", function(p, ...) {
  args <- list(...)
  for (i in seq_along(args)) {
    a <- args[[i]]
    if (is.numeric(a) || (is_sym(a) && is_number_type(a))) {
      args[[i]] <- if (isTRUE(.is_subset_int(.callm(p, "index_set")))) .callm(p, "getitem", a) else .callm(p, "call", a)
    } else if (!(is_sym(a) && a$type == "RandomIndexedSymbol")) {
      .serr("ValueError", paste0("Expected a RandomIndexedSymbol or key not  ", class(a)[1]))
    }
  }
  if (identical(.rs_pspace(args[[1L]])$args[[3L]]$type, "Distribution")) return(new_sym("JointDistribution", args))
  dens <- Lambda(args, .mul(lapply(args, function(a) .callm(.rs_pspace(a)$args[[2L]], "density", a))))
  new_sym("JointDistributionHandmade", list(dens))
})
.is_subset_int <- function(s) identical(s$key, Naturals0$key) || identical(s$key, Integers$key) || identical(s$key, Naturals$key)

.defm("DiscreteTimeStochasticProcess", "getitem", function(p, time) {
  time <- as_sym(time)
  if (!(time$type %in% c("Symbol", "Dummy"))) {
    inside <- .fuzzy_bool(Contains(time, .callm(p, "index_set")))
    if (!isTRUE(inside)) .serr("IndexError", paste0(.str(time), " is not in the index set of ", .str(p$args[[1L]])))
  }
  idx <- Indexed(p$args[[1L]]$name, time)
  .RandomIndexedSymbol(idx, .StochasticPSpace(p$args[[1L]], p, .callm(p, "distribution", time)))
})
.defm("ContinuousTimeStochasticProcess", "call", function(p, time) {
  time <- as_sym(time)
  if (!(time$type %in% c("Symbol", "Dummy"))) {
    inside <- .fuzzy_bool(Contains(time, .callm(p, "index_set")))
    if (!isTRUE(inside)) .serr("IndexError", paste0(.str(time), " is not in the index set of ", .str(p$args[[1L]])))
  }
  fo <- fcall(p$args[[1L]]$name, list(time))
  .RandomIndexedSymbol(fo, .StochasticPSpace(p$args[[1L]], p, .callm(p, "distribution", time)))
})

# _SubstituteRV
.rvindexed_subs <- function(expr, condition = NULL) {
  swap <- function(e) {
    if (is.null(e) || !is_sym(e)) return(e)
    for (rv in random_symbols(e))
      if (rv$type == "RandomIndexedSymbol") e <- subs(e, rv, .callm(.rs_pspace(rv)$args[[2L]], "simple_rv", rv))
    e
  }
  list(swap(expr), swap(condition))
}

.subrv_expectation <- function(expr, condition = NULL, evaluate = TRUE) {
  nc <- .rvindexed_subs(expr, condition)
  new_expr <- nc[[1L]]; new_cond <- nc[[2L]]
  if (!is_random(new_expr)) return(new_expr)
  new_ps <- pspace(new_expr)
  if (!is.null(new_cond)) new_expr <- given(new_expr, new_cond)
  if (new_expr$type == "Add")
    return(.add(lapply(new_expr$args, function(a) .callm(new_ps, "compute_expectation", a, evaluate = evaluate))))
  .callm(new_ps, "compute_expectation", new_expr, evaluate = evaluate)
}

.subrv_probability <- function(condition, given_condition = NULL, evaluate = TRUE) {
  nc <- .rvindexed_subs(as_sym(condition), if (is.null(given_condition)) NULL else as_sym(given_condition))
  new_cond <- nc[[1L]]; new_given <- nc[[2L]]
  if (!is.null(new_given) && .is_rs(new_given)) {
    condrv <- random_symbols(new_cond)
    if (length(condrv) == 1L && identical(condrv[[1L]]$key, new_given$key))
      return(new_sym("BernoulliDistribution", list(.subrv_probability(new_cond), S_Zero, S_One)))
    if (any(vapply(condrv, function(rv) isTRUE(dependent(rv, new_given)), NA))) return(Probability(new_cond, new_given))
    return(.subrv_probability(new_cond))
  }
  if (!is.null(new_given) && !.is_cond(new_given))
    .serr("ValueError", paste0(.str(new_given), " is not a relational or combination of relationals"))
  if ((!is.null(new_given) && .is_false(new_given)) || .is_false(new_cond)) return(S_Zero)
  if (.is_true(new_cond)) return(S_One)
  if (!.is_cond(new_cond))
    .serr("ValueError", paste0(.str(new_cond), " is not a relational or combination of relationals"))
  if (!is.null(new_given)) return(.subrv_probability(given(new_cond, new_given)))
  result <- .callm(pspace(new_cond), "probability", new_cond)
  if (evaluate && is_sym(result)) doit(result) else result
}

# ---------------------------------------------------------------------------
# BernoulliProcess

.bp_arg <- function(p, i) p$args[[i]]
.defm("BernoulliProcess", "index_set", function(p) Naturals0)
.defm("BernoulliProcess", "p", function(p) p$args[[2L]])
.defm("BernoulliProcess", "success", function(p) p$args[[3L]])
.defm("BernoulliProcess", "failure", function(p) p$args[[4L]])
.defm("BernoulliProcess", "state_space", function(p) FiniteSet(p$args[[3L]], p$args[[4L]]))
.defm("BernoulliProcess", "distribution", function(p, key = NULL) {
  if (is.null(key)) return(new_sym("BernoulliDistribution", list(p$args[[2L]], S_One, S_Zero)))
  new_sym("BernoulliDistribution", list(p$args[[2L]], p$args[[3L]], p$args[[4L]]))
})
.defm("BernoulliProcess", "simple_rv", function(p, rv) Bernoulli(.str(rv$args[[1L]]), p$args[[2L]], p$args[[3L]], p$args[[4L]]))
.defm("BernoulliProcess", "expectation", function(p, expr, condition = NULL, evaluate = TRUE, ...) .subrv_expectation(expr, condition, evaluate))
.defm("BernoulliProcess", "probability", function(p, condition, given_condition = NULL, evaluate = TRUE, ...) .subrv_probability(condition, given_condition, evaluate))
.defm("BernoulliProcess", "density", function(p, x)
  Piecewise(list(p$args[[2L]], Eq(x, p$args[[3L]])), list(1 - p$args[[2L]], Eq(x, p$args[[4L]])), list(S_Zero, S_true)))

.sym_sympify <- function(a) if (is.character(a)) Symbol(a) else as_sym(a)

#' Stochastic processes
#'
#' `BernoulliProcess()`, `PoissonProcess()`, `WienerProcess()` and
#' `GammaProcess()` of `sympy.stats`. A discrete time process is indexed
#' with [rv_getitem()] (SymPy's `X[t]`); a continuous time process is
#' called, `X(t)`. Process attributes such as `"state_space"`,
#' `"index_set"` or `"joint_distribution"` are read with [rv_attr()].
#'
#' @param sym The name of the process.
#' @param p,success,failure Probability and outcomes of a Bernoulli trial.
#' @param lamda,gamma Rate and shape parameters.
#' @param pspace A probability space.
#' @return A process (an R function), or a random indexed symbol.
#' @examples
#' B <- BernoulliProcess("B", p = Rational(1, 2))
#' E(rv_getitem(B, 1) + rv_getitem(B, 2))
#' X <- PoissonProcess("X", 3)
#' E(X(2))
#' @name stochastic_process
NULL

#' @rdname stochastic_process
#' @export
BernoulliProcess <- function(sym, p, success = 1, failure = 0) {
  p <- as_sym(p)
  .value_check(list(p >= 0, p <= 1), "Value of p must be between 0 and 1.")
  .mkproc(new_sym("BernoulliProcess", list(.symbol_converter(sym), p, .sym_sympify(success), .sym_sympify(failure))))
}

# ---------------------------------------------------------------------------
# Counting processes

.defm("CountingProcess", "index_set", function(p) Interval(0, S_Infinity))
.defm("CountingProcess", "expectation", function(p, expr, condition = NULL, evaluate = TRUE, ...) {
  if (!is.null(condition)) {
    sw <- .timerv_swaps(expr, condition)
    intervals <- sw[[1L]]
    disjoint <- length(intervals) == 1L || all(unlist(lapply(seq_along(intervals), function(i) lapply(seq_along(intervals), function(j)
      if (i < j) identical(Intersection(intervals[[i]], intervals[[j]])$key, S_EmptySet$key) else TRUE))))
    if (disjoint) {
      if (expr$type == "Add") return(.add(lapply(expr$args, function(a) .callm(p, "expectation", a, condition))))
      for (pr in sw[[2L]]) expr <- subs(expr, pr[[1L]], pr[[2L]])
    } else {
      return(Expectation(expr, condition))
    }
  }
  .subrv_expectation(expr, evaluate = evaluate)
})

.timerv_swaps <- function(expr, condition) {
  # Contains is a Boolean in SymPy
  if (!(.is_cond(condition) || (is_sym(condition) && condition$type == "Contains")))
    .serr("ValueError", paste0(.str(condition), " is not a relational or combination of relationals"))
  expr_syms <- Filter(function(z) z$type == "RandomIndexedSymbol", .uniq(.all_subexprs(expr)))
  args <- if (condition$type %in% c("And", "Or")) condition$args else list(condition)
  swaps <- list(); intervals <- list()
  for (es in expr_syms) for (a in args) {
    key <- .ris_key(es)
    if (.has1(a, key) && key$type %in% c("Symbol", "Dummy")) {
      intv <- a$args[[2L]]
      dk <- set_sup(intv) - set_inf(intv)
      if (identical(dk$key, S_Infinity$key)) .serr("ValueError", paste0(.str(es), " should have finite bounds"))
      if (identical(dk$key, S_Zero$key)) dk <- set_sup(intv)
      swaps[[length(swaps) + 1L]] <- list(es, subs(es, key, dk))
      intervals[[length(intervals) + 1L]] <- intv
    }
  }
  list(intervals, swaps)
}

.defm("CountingProcess", "probability", function(p, condition, given_condition = NULL, evaluate = TRUE, ...) {
  cond_args <- if (condition$type %in% c("And", "Or")) condition$args else list(condition)
  check_numeric <- all(vapply(cond_args, function(a) .rv_is_number(.ris_key(a$args[[1L]])), NA))
  if (!is.null(given_condition)) {
    gargs <- if (given_condition$type %in% c("And", "Or")) given_condition$args else list(given_condition)
    check_given <- !.has_type(given_condition, "Contains")
    if (check_numeric && check_given) .serr("NotImplementedError", "numeric given conditions")
    if (!all(vapply(gargs, function(a) .has_type(a, "Contains"), NA)))
      .serr("ValueError", "If given condition is passed with `Contains`, then please pass the evaluated condition with its corresponding information in terms of intervals of each time stamp to be passed in given condition.")
    sw <- .timerv_swaps(condition, given_condition)
    intervals <- sw[[1L]]
    disjoint <- length(intervals) == 1L || all(unlist(lapply(seq_along(intervals), function(i) lapply(seq_along(intervals), function(j)
      if (i < j) identical(Intersection(intervals[[i]], intervals[[j]])$key, S_EmptySet$key) else TRUE))))
    if (!disjoint) return(Probability(condition, given_condition))
    if (condition$type == "And") return(.mul(lapply(condition$args, function(a) .callm(p, "probability", a, given_condition))))
    if (condition$type == "Or") return(.add(lapply(condition$args, function(a) .callm(p, "probability", a, given_condition))))
    for (pr in sw[[2L]]) condition <- subs(condition, pr[[1L]], pr[[2L]])
  }
  .subrv_probability(condition, evaluate = evaluate)
})

.defm("PoissonProcess", "lamda", function(p) p$args[[2L]])
.defm("PoissonProcess", "state_space", function(p) Naturals0)
.defm("PoissonProcess", "distribution", function(p, key) new_sym("PoissonDistribution", list(p$args[[2L]] * as_sym(key))))
.defm("PoissonProcess", "density", function(p, x) { k <- .ris_key(x); (p$args[[2L]] * k)^x / sym_factorial(x) * exp(-(p$args[[2L]] * k)) })
.defm("PoissonProcess", "simple_rv", function(p, rv) Poisson(.str(rv), p$args[[2L]] * .ris_key(rv)))

#' @rdname stochastic_process
#' @export
PoissonProcess <- function(sym, lamda) {
  lamda <- as_sym(lamda)
  .value_check(lamda > 0, "lamda should be a positive number.")
  .mkproc(new_sym("PoissonProcess", list(.symbol_converter(sym), lamda)))
}

.defm("WienerProcess", "state_space", function(p) Reals)
.defm("WienerProcess", "distribution", function(p, key) new_sym("NormalDistribution", list(S_Zero, sqrt(as_sym(key)))))
.defm("WienerProcess", "density", function(p, x) { k <- .ris_key(x); exp(-x^2 / (2 * k)) / (sqrt(2 * S_Pi) * sqrt(k)) })
.defm("WienerProcess", "simple_rv", function(p, rv) Normal(.str(rv), 0, sqrt(.ris_key(rv))))

#' @rdname stochastic_process
#' @export
WienerProcess <- function(sym) .mkproc(new_sym("WienerProcess", list(.symbol_converter(sym))))

.defm("GammaProcess", "lamda", function(p) p$args[[2L]])
.defm("GammaProcess", "gamma", function(p) p$args[[3L]])
.defm("GammaProcess", "state_space", function(p) Interval(0, S_Infinity))
.defm("GammaProcess", "distribution", function(p, key) new_sym("GammaDistribution", list(p$args[[3L]] * as_sym(key), 1 / p$args[[2L]])))
.defm("GammaProcess", "density", function(p, x) {
  k <- p$args[[3L]] * .ris_key(x); theta <- 1 / p$args[[2L]]
  x^(k - 1) * exp(-x / theta) / (gamma(k) * theta^k)
})
.defm("GammaProcess", "simple_rv", function(p, rv) GammaRV(.str(rv), p$args[[3L]] * .ris_key(rv), 1 / p$args[[2L]]))

#' @rdname stochastic_process
#' @export
GammaProcess <- function(sym, lamda, gamma) {
  lamda <- as_sym(lamda); gamma <- as_sym(gamma)
  .value_check(lamda > 0, "lamda should be a positive number")
  .value_check(gamma > 0, "gamma should be a positive number")
  .mkproc(new_sym("GammaProcess", list(.symbol_converter(sym), lamda, gamma)))
}
