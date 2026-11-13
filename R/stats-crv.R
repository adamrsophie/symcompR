# Continuous probability spaces, after sympy/stats/crv.py.

.defclass("ContinuousDomain", "RandomDomain")
.defclass("SingleContinuousDomain", c("ContinuousDomain", "SingleDomain"))
.defclass("ProductContinuousDomain", c("ProductDomain", "ContinuousDomain"))
.defclass("ConditionalContinuousDomain", c("ContinuousDomain", "ConditionalDomain"))
.defclass("ContinuousDistribution", "Distribution")
.defclass("SingleContinuousDistribution", "ContinuousDistribution")
.defclass("ContinuousPSpace", "PSpace")
.defclass("SingleContinuousPSpace", c("ContinuousPSpace", "SinglePSpace"))

.defm("ContinuousDomain", "is_Continuous", function(d) TRUE)
.defm("ContinuousPSpace", "is_Continuous", function(ps) TRUE)
.defm("ContinuousPSpace", "is_real", function(ps) TRUE)
.defm("PSpace", "is_Continuous", function(ps) NA)
.defm("PSpace", "is_Discrete", function(ps) NA)
.defm("PSpace", "is_Finite", function(ps) NA)
.defm("PSpace", "is_real", function(ps) NA)

# Integral limits over a set: SymPy turns (x, Interval) into (x, inf, sup)
.set_lims <- function(x, s) {
  if (s$type %in% c("Interval", "SpecialSet")) return(list(x, set_inf(s), set_sup(s)))
  list(x, s)
}

.domain_new <- function(cls, symbol, set) new_sym(cls, list(symbol, set))

.defm("ContinuousDomain", "as_boolean", function(d) .serr("NotImplementedError", "Not Implemented for generic Domains"))

.defm("SingleContinuousDomain", "compute_expectation", function(d, expr, variables = NULL, evaluate = TRUE, ...) {
  if (is.null(variables)) variables <- .dom_symbols(d)
  # a set of symbols may come as a FiniteSet
  if (is_sym(variables) && variables$type == "FiniteSet") variables <- variables$args
  if (!length(variables)) return(expr)
  vk <- sort(unique(vapply(variables, .k, ""))); sk <- sort(unique(vapply(.dom_symbols(d), .k, "")))
  if (!identical(vk, sk)) .serr("ValueError", "Values should be equal")
  do.call(Integral, c(list(expr), list(.set_lims(d$args[[1L]], d$args[[2L]]))))
})
.defm("SingleContinuousDomain", "as_boolean", function(d) as_relational(d$args[[2L]], d$args[[1L]]))

.defm("ProductContinuousDomain", "compute_expectation", function(d, expr, variables = NULL, ...) {
  if (is.null(variables)) variables <- .dom_symbols(d)
  for (domain in d$args) {
    dv <- Filter(function(v) .in_list(v, .dom_symbols(domain)), variables)
    if (length(dv)) expr <- .callm(domain, "compute_expectation", expr, dv, ...)
  }
  expr
})
.defm("ProductContinuousDomain", "as_boolean", function(d) do.call(And, lapply(d$args, .dom_as_boolean)))

.ConditionalContinuousDomain <- function(fulldomain, condition) .ConditionalDomain_new("ConditionalContinuousDomain", fulldomain, condition)

.integral_parts <- function(I) {
  if (is_sym(I) && I$type == "Integral") return(list(I$args[[1L]], lapply(I$args[-1L], function(t) t$args)))
  list(I, list())
}

.defm("ConditionalContinuousDomain", "compute_expectation", function(d, expr, variables = NULL, ...) {
  if (is.null(variables)) variables <- .dom_symbols(d)
  if (!length(variables)) return(expr)
  full <- .callm(d$args[[1L]], "compute_expectation", expr, variables)
  pr <- .integral_parts(full)
  integrand <- pr[[1L]]; limits <- pr[[2L]]
  conditions <- list(d$args[[2L]])
  syms <- .dom_symbols(d)
  while (length(conditions)) {
    cond <- conditions[[length(conditions)]]
    conditions[[length(conditions)]] <- NULL
    if (is_sym(cond) && is_relational(cond)) {
      if (cond$type == "Equality") {
        integrand <- .mul(list(integrand, DiracDelta(.sub(cond$args[[1L]], cond$args[[2L]]))))
      } else {
        fs <- Filter(function(s) .in_list(s, syms), free_symbols(cond))
        if (length(fs) != 1L) .serr("NotImplementedError", "Multivariate Inequalities not yet implemented")
        symbol <- fs[[1L]]
        for (i in seq_along(limits)) {
          lim <- limits[[i]]
          if (identical(lim[[1L]]$key, symbol$key)) {
            cintvl <- .reduce_wrap(cond, symbol)
            lintvl <- Interval(lim[[2L]], lim[[3L]])
            intvl <- Intersection(cintvl, lintvl)
            limits[[i]] <- list(symbol, .set_left(intvl), .set_right(intvl))
          }
        }
      }
    } else if (is_sym(cond) && cond$type %in% c("And", "Or", "Not", "BooleanAtom")) {
      if (cond$type == "And") conditions <- c(conditions, .sympy_args(cond))
      else if (cond$type == "Or") .serr("NotImplementedError", "Or not implemented here")
    } else {
      .serr("TypeError", paste0("Condition ", .str(cond), " is not a relational or Boolean"))
    }
  }
  do.call(Integral, c(list(integrand), limits))
})
.defm("ConditionalContinuousDomain", "as_boolean", function(d) And(.dom_as_boolean(d$args[[1L]]), d$args[[2L]]))
.defm("ConditionalContinuousDomain", "set", function(d) {
  syms <- .dom_symbols(d)
  if (length(syms) == 1L) return(Intersection(.dom_set(d$args[[1L]]), .reduce_wrap(d$args[[2L]], syms[[1L]])))
  .serr("NotImplementedError", "Set of Conditional Domain not Implemented")
})

# Interval.left and .right
.set_left <- function(s) if (s$type == "Interval") .iv(s)$a else set_inf(s)
.set_right <- function(s) if (s$type == "Interval") .iv(s)$b else set_sup(s)

# reduce_rational_inequalities_wrap
.reduce_ineq <- function(cond, var) {
  tryCatch(reduce_rational_inequalities(list(list(cond)), var, relational = FALSE),
           error = function(e) {
             if (inherits(e, c("ValueError", "NotImplementedError", "TypeError"))) stop(e)
             .serr("ValueError", paste0("Reduction of condition failed ", .rv_pystr(list(cond)), "\n"))
           })
}

.reduce_wrap <- function(condition, var) {
  condition <- .as_bool_sym(condition)
  if (is_relational(condition)) return(.reduce_ineq(condition, var))
  if (condition$type == "Or") return(do.call(Union, lapply(.sympy_args(condition), .reduce_ineq, var = var)))
  if (condition$type == "And") {
    ints <- lapply(.sympy_args(condition), .reduce_ineq, var = var)
    I <- ints[[1L]]
    for (i in ints) I <- Intersection(I, i)
    return(I)
  }
  NULL
}

# ---------------------------------------------------------------------------
# SingleContinuousDistribution

.defm("SingleContinuousDistribution", "set", function(d) Interval(S_NegativeInfinity, S_Infinity))
.defm("ContinuousDistribution", "call", function(d, ...) .callm(d, "pdf", ...))

.has_m <- function(d, name) !is.null(.findm(d, name))

.dsym <- function(name, ...) Dummy(name, ...)

.defm("SingleContinuousDistribution", "compute_cdf", function(d, ...) {
  x <- Dummy("x", real = TRUE); z <- Dummy("z", real = TRUE)
  left <- set_inf(.dist_set(d))
  pdf <- .dist_pdf(d, x)
  cdf <- integrate(doit(pdf), list(x, left, z))
  cdf <- Piecewise(list(cdf, Ge(z, left)), list(S_Zero, S_true))
  Lambda(z, cdf)
})

.opt_method <- function(d, name, ...) {
  f <- .findm(d, name)
  if (is.null(f)) return(NULL)
  f(d, ...)
}

.dist_cdf <- function(d, x, ...) {
  if (!...length()) {
    c <- .opt_method(d, "_cdf", x)
    if (!is.null(c)) return(c)
  }
  apply_lambda(.callm(d, "compute_cdf", ...), x)
}
.defm("SingleContinuousDistribution", "cdf", function(d, x, ...) .dist_cdf(d, x, ...))

.defm("SingleContinuousDistribution", "compute_characteristic_function", function(d, ...) {
  x <- Dummy("x", real = TRUE); t <- Dummy("t", real = TRUE)
  pdf <- .dist_pdf(d, x)
  cf <- do.call(integrate, list(.mul(list(exp(.mul(list(S_I, t, x))), pdf)), .set_lims(x, .dist_set(d))))
  Lambda(t, cf)
})
.defm("SingleContinuousDistribution", "characteristic_function", function(d, t, ...) {
  if (!...length()) {
    c <- .opt_method(d, "_characteristic_function", t)
    if (!is.null(c)) return(c)
  }
  apply_lambda(.callm(d, "compute_characteristic_function", ...), t)
})

.defm("SingleContinuousDistribution", "compute_moment_generating_function", function(d, ...) {
  x <- Dummy("x", real = TRUE); t <- Dummy("t", real = TRUE)
  pdf <- .dist_pdf(d, x)
  mgf <- do.call(integrate, list(.mul(list(exp(.mul(list(t, x))), pdf)), .set_lims(x, .dist_set(d))))
  Lambda(t, mgf)
})
.defm("SingleContinuousDistribution", "moment_generating_function", function(d, t, ...) {
  if (!...length()) {
    c <- .opt_method(d, "_moment_generating_function", t)
    if (!is.null(c)) return(c)
  }
  apply_lambda(.callm(d, "compute_moment_generating_function", ...), t)
})

# poly(expr, var): NULL when expr is not a polynomial in var
.poly_or_null <- function(expr, var) {
  if (!isTRUE(tryCatch(is_polynomial(expr, var), error = function(e) FALSE))) return(NULL)
  tryCatch(Poly(expr, var), error = function(e) NULL)
}

.poly_coeffs_by_degree <- function(p) {
  cs <- all_coeffs(p)
  rev(cs)
}

.defm("SingleContinuousDistribution", "expectation", function(d, expr, var, evaluate = TRUE, ...) {
  if (evaluate) {
    p <- .poly_or_null(expr, var)
    if (is.null(p)) return(do.call(integrate, list(.mul(list(expr, .dist_pdf(d, var))), .set_lims(var, .dist_set(d)))))
    cs <- .poly_coeffs_by_degree(p)
    if (all(vapply(cs, is_zero_num, NA))) return(S_Zero)
    t <- Dummy("t", real = TRUE)
    mgf <- .opt_method(d, "_moment_generating_function", t)
    if (is.null(mgf)) return(do.call(integrate, list(.mul(list(expr, .dist_pdf(d, var))), .set_lims(var, .dist_set(d)))))
    deg <- length(cs) - 1L
    taylor <- removeO(series(mgf, t, 0, deg + 1L))
    tp <- .poly_or_null(taylor, t)
    if (is.null(tp)) return(do.call(integrate, list(.mul(list(expr, .dist_pdf(d, var))), .set_lims(var, .dist_set(d)))))
    tc <- .poly_coeffs_by_degree(tp)
    result <- S_Zero
    for (k in 0:deg) {
      ck <- cs[[k + 1L]]
      tk <- if (k + 1L <= length(tc)) tc[[k + 1L]] else S_Zero
      result <- .add(list(result, .mul(list(ck, tk, sym_factorial(as_sym(k))))))
    }
    return(result)
  }
  do.call(Integral, list(.mul(list(expr, .dist_pdf(d, var))), .set_lims(var, .dist_set(d))))
})

.defm("SingleContinuousDistribution", "compute_quantile", function(d, ...) {
  x <- Dummy("x", real = TRUE); p <- Dummy("p", real = TRUE)
  left <- set_inf(.dist_set(d))
  pdf <- .dist_pdf(d, x)
  cdf <- integrate(pdf, list(x, left, x))
  q <- solveset(.sub(cdf, p), x, .dist_set(d))
  Lambda(p, Piecewise(list(q, And(Ge(p, 0), Le(p, 1))), list(S_NaN, S_true)))
})
.defm("SingleContinuousDistribution", "quantile", function(d, x, ...) {
  if (!...length()) {
    q <- .opt_method(d, "_quantile", x)
    if (!is.null(q)) return(q)
  }
  apply_lambda(.callm(d, "compute_quantile", ...), x)
})

# ---------------------------------------------------------------------------
# ContinuousPSpace

.rv_lcall <- function(f, ...) {
  if (is_sym(f) && f$type == "Lambda") return(apply_lambda(f, ...))
  .rv_call(f, list(...))
}

.defm("ContinuousPSpace", "pdf", function(ps) do.call(.rv_lcall, c(list(.callm(ps, "density")), .dom_symbols(.ps_domain(ps)))))

.defm("ContinuousPSpace", "compute_expectation", function(ps, expr, rvs = NULL, evaluate = FALSE, ...) {
  if (is.null(rvs)) rvs <- .ps_values(ps)
  expr <- xreplace(as_sym(expr), lapply(rvs, function(rv) list(rv, .rs_symbol(rv))))
  dsyms <- lapply(rvs, .rs_symbol)
  .callm(.ps_domain(ps), "compute_expectation", .mul(list(.ps_pdf(ps), expr)), dsyms, ...)
})

.defm("ContinuousPSpace", "compute_density", function(ps, expr, ...) {
  vals <- .ps_values(ps)
  if (.in_list(expr, vals)) {
    others <- Filter(function(r) !identical(r$key, expr$key), vals)
    syms <- lapply(others, .rs_symbol)
    pdf <- .callm(.ps_domain(ps), "compute_expectation", .ps_pdf(ps), syms, ...)
    return(Lambda(.rs_symbol(expr), pdf))
  }
  z <- Dummy("z", real = TRUE)
  Lambda(z, .callm(ps, "compute_expectation", DiracDelta(.sub(expr, z)), ...))
})

.ps_set_is_interval <- function(ps) {
  s <- tryCatch(.dom_set(.ps_domain(ps)), error = function(e) NULL)
  !is.null(s) && (s$type == "Interval" || (s$type == "SpecialSet" && s$name == "Reals"))
}

.defm("ContinuousPSpace", "compute_cdf", function(ps, expr, ...) {
  if (!.ps_set_is_interval(ps)) .serr("ValueError", "CDF not well defined on multivariate expressions")
  dd <- .callm(ps, "compute_density", expr, ...)
  x <- Dummy("x", real = TRUE); z <- Dummy("z", real = TRUE)
  left <- set_inf(.dom_set(.ps_domain(ps)))
  cdf <- integrate(.rv_lcall(dd, x), list(x, left, z))
  cdf <- Piecewise(list(cdf, Ge(z, left)), list(S_Zero, S_true))
  Lambda(z, cdf)
})

.defm("ContinuousPSpace", "compute_characteristic_function", function(ps, expr, ...) {
  if (!.ps_set_is_interval(ps)) .serr("NotImplementedError", "Characteristic function of multivariate expressions not implemented")
  dd <- .callm(ps, "compute_density", expr, ...)
  x <- Dummy("x", real = TRUE); t <- Dummy("t", real = TRUE)
  cf <- integrate(.mul(list(exp(.mul(list(S_I, t, x))), .rv_lcall(dd, x))), list(x, S_NegativeInfinity, S_Infinity))
  Lambda(t, cf)
})

.defm("ContinuousPSpace", "compute_moment_generating_function", function(ps, expr, ...) {
  if (!.ps_set_is_interval(ps)) .serr("NotImplementedError", "Moment generating function of multivariate expressions not implemented")
  dd <- .callm(ps, "compute_density", expr, ...)
  x <- Dummy("x", real = TRUE); t <- Dummy("t", real = TRUE)
  mgf <- integrate(.mul(list(exp(.mul(list(t, x))), .rv_lcall(dd, x))), list(x, S_NegativeInfinity, S_Infinity))
  Lambda(t, mgf)
})

.defm("ContinuousPSpace", "compute_quantile", function(ps, expr, ...) {
  if (!.ps_set_is_interval(ps)) .serr("ValueError", "Quantile not well defined on multivariate expressions")
  dd <- .callm(ps, "compute_cdf", expr, ...)
  x <- Dummy("x", real = TRUE); p <- Dummy("p", positive = TRUE)
  q <- solveset(.sub(apply_lambda(dd, x), p), x, .callm(ps, "set"))
  Lambda(p, q)
})

.defm("ContinuousPSpace", "probability", function(ps, condition, ...) {
  z <- Dummy("z", real = TRUE)
  cond_inv <- FALSE
  if (condition$type == "Unequality") {
    condition <- Eq(condition$args[[1L]], condition$args[[2L]])
    cond_inv <- TRUE
  }
  res <- tryCatch({
    domain <- .callm(ps, "where", condition)
    rv <- Filter(function(r) identical(.rs_symbol(r)$key, domain$args[[1L]]$key), .ps_values(ps))[[1L]]
    pdf <- .callm(ps, "compute_density", rv, ...)
    dset <- .dom_set(domain)
    if (dset$type %in% c("EmptySet", "FiniteSet")) return(if (!cond_inv) S_Zero else S_One)
    if (dset$type == "Union") {
      parts <- lapply(Filter(function(s) s$type == "Interval", dset$args),
                      function(s) Integral(.rv_lcall(pdf, z), .set_lims(z, s)))
      return(.add(c(list(S_Zero), parts)))
    }
    Integral(.rv_lcall(pdf, z), .set_lims(z, dset))
  }, NotImplementedError = function(e) NULL)
  if (!is.null(res)) return(res)
  expr <- .sub(condition$args[[1L]], condition$args[[2L]])
  if (!is_random(expr)) {
    dens <- .callm(ps, "density")
    comp <- condition$args[[2L]]
  } else {
    dens <- rv_density(expr, ...)
    comp <- S_Zero
  }
  if (!.isa(dens, "ContinuousDistribution")) dens <- .ContinuousDistributionHandmade(dens, .dom_set(.ps_domain(ps)))
  space <- .SingleContinuousPSpace(z, dens)
  result <- .callm(space, "probability", .rel_same(condition, .ps_value(space), comp))
  if (!cond_inv) result else .sub(S_One, result)
})

.defm("ContinuousPSpace", "where", function(ps, condition) {
  rvs <- random_symbols(condition)
  if (!(length(rvs) == 1L && .in_list(rvs[[1L]], .ps_values(ps))))
    .serr("NotImplementedError", "Multiple continuous random variables not supported")
  rv <- rvs[[1L]]
  interval <- .reduce_wrap(condition, rv)
  interval <- Intersection(interval, .dom_set(.ps_domain(ps)))
  .domain_new("SingleContinuousDomain", .rs_symbol(rv), interval)
})

.defm("ContinuousPSpace", "conditional_space", function(ps, condition, normalize = TRUE, ...) {
  condition <- xreplace(.as_bool_sym(condition), lapply(.ps_values(ps), function(rv) list(rv, .rs_symbol(rv))))
  domain <- .ConditionalContinuousDomain(.ps_domain(ps), condition)
  density <- NULL
  if (normalize) {
    replacement <- lapply(.ps_symbols(ps), function(s) list(s, Dummy(.str(s))))
    norm <- .callm(domain, "compute_expectation", .ps_pdf(ps), ...)
    pdf <- .div(.ps_pdf(ps), xreplace(norm, replacement))
    density <- Lambda(.dom_symbols(domain), pdf)
  }
  new_sym("ContinuousPSpace", list(domain, density))
})

# ---------------------------------------------------------------------------
# SingleContinuousPSpace

.SingleContinuousPSpace <- function(s, distribution) new_sym("SingleContinuousPSpace", list(.symbol_converter(s), distribution))

.defm("SingleContinuousPSpace", "set", function(ps) .dist_set(ps$args[[2L]]))
.defm("SingleContinuousPSpace", "domain", function(ps) .domain_new("SingleContinuousDomain", ps$args[[1L]], .dist_set(ps$args[[2L]])))
.defm("SingleContinuousPSpace", "density", function(ps) ps$args[[2L]])
.defm("SingleContinuousPSpace", "pdf", function(ps) .dist_pdf(ps$args[[2L]], ps$args[[1L]]))

.defm("SingleContinuousPSpace", "compute_expectation", function(ps, expr, rvs = NULL, evaluate = FALSE, ...) {
  val <- .ps_value(ps)
  if (is.null(rvs) || !length(rvs)) rvs <- list(val)
  if (!.in_list(val, rvs)) return(expr)
  expr <- xreplace(as_sym(expr), lapply(rvs, function(rv) list(rv, .rs_symbol(rv))))
  x <- .rs_symbol(val)
  tryCatch(.callm(ps$args[[2L]], "expectation", expr, x, evaluate = evaluate, ...),
           PoleError = function(e) do.call(Integral, list(.mul(list(expr, .ps_pdf(ps))), .set_lims(x, .callm(ps, "set")))))
})

.defm("SingleContinuousPSpace", "compute_cdf", function(ps, expr, ...) {
  if (.eq_struct(expr, .ps_value(ps))) {
    z <- Dummy("z", real = TRUE)
    return(Lambda(z, .dist_cdf(ps$args[[2L]], z, ...)))
  }
  .STM[["ContinuousPSpace"]]$compute_cdf(ps, expr, ...)
})
.defm("SingleContinuousPSpace", "compute_characteristic_function", function(ps, expr, ...) {
  if (.eq_struct(expr, .ps_value(ps))) {
    t <- Dummy("t", real = TRUE)
    return(Lambda(t, .callm(ps$args[[2L]], "characteristic_function", t, ...)))
  }
  .STM[["ContinuousPSpace"]]$compute_characteristic_function(ps, expr, ...)
})
.defm("SingleContinuousPSpace", "compute_moment_generating_function", function(ps, expr, ...) {
  if (.eq_struct(expr, .ps_value(ps))) {
    t <- Dummy("t", real = TRUE)
    return(Lambda(t, .callm(ps$args[[2L]], "moment_generating_function", t, ...)))
  }
  .STM[["ContinuousPSpace"]]$compute_moment_generating_function(ps, expr, ...)
})
.defm("SingleContinuousPSpace", "compute_density", function(ps, expr, ...) {
  val <- .ps_value(ps)
  if (.eq_struct(expr, val)) return(ps$args[[2L]])
  y <- Dummy("y", real = TRUE)
  gs <- solveset(.sub(expr, y), val, Reals)
  if (gs$type == "Intersection" && length(gs$args) == 2L && identical(gs$args[[1L]]$key, Reals$key)) gs <- gs$args[[2L]]
  if (gs$type != "FiniteSet") .serr("ValueError", paste0("Can not solve ", .str(expr), " for ", .str(val)))
  fx <- .callm(ps, "compute_density", val)
  fy <- .add(c(list(S_Zero), lapply(gs$args, function(g) .mul(list(.rv_lcall(fx, g), Abs(diff(g, y)))))))
  Lambda(y, fy)
})
.defm("SingleContinuousPSpace", "compute_quantile", function(ps, expr, ...) {
  if (.eq_struct(expr, .ps_value(ps))) {
    p <- Dummy("p", real = TRUE)
    return(Lambda(p, .callm(ps$args[[2L]], "quantile", p, ...)))
  }
  .STM[["ContinuousPSpace"]]$compute_quantile(ps, expr, ...)
})
