# manualintegrate, after sympy/integrals/manualintegrate.py (part 1): the
# integration rules as objects with eval() and contains_dont_know().
# A rule is list(type, integrand, variable, ...) of class "mrule". The rule
# functions of SymPy are generators that yield subintegrals to the
# IntegrationSolver; here they receive the solver S and call S$solve().

.mr <- function(type, integrand, variable, ...) structure(c(list(type = type, integrand = if (is.null(integrand)) NULL else as_sym(integrand),
                                                                 variable = if (is.null(variable)) NULL else as_sym(variable)), list(...)), class = "mrule")
.is_mr <- function(r, type = NULL) inherits(r, "mrule") && (is.null(type) || r$type %in% type)

.MR_ATOMIC <- c("ConstantRule", "PowerRule", "NestedPowRule", "SinRule", "CosRule", "SinhRule", "CoshRule", "ExpRule",
                "ReciprocalRule", "ArcsinRule", "ArcsinhRule", "ReciprocalSqrtQuadraticRule", "SqrtQuadraticDenomRule",
                "RatintRule", "DerivativeRule", "DiracDeltaRule", "ArctanRule", "JacobiRule", "GegenbauerRule",
                "ChebyshevTRule", "ChebyshevURule", "LegendreRule", "HermiteRule", "LaguerreRule", "AssocLaguerreRule",
                "CiRule", "ChiRule", "EiRule", "SiRule", "ShiRule", "LiRule", "ErfRule", "OwensTRule", "FresnelCRule",
                "FresnelSRule", "PolylogRule", "UpperGammaRule", "EllipticFRule", "EllipticERule")

mr_contains_dont_know <- function(r) {
  if (r$type == "DontKnowRule") return(TRUE)
  if (r$type %in% .MR_ATOMIC) return(FALSE)
  switch(r$type,
    ConstantTimesRule = , URule = , ReparameterizationRule = , SqrtQuadraticRule = , RewriteRule = , CompleteSquareRule = ,
    HeavisideRule = mr_contains_dont_know(r$substep),
    AddRule = any(vapply(r$substeps, mr_contains_dont_know, NA)),
    PartsRule = mr_contains_dont_know(r$v_step) || (!is.null(r$second_step) && mr_contains_dont_know(r$second_step)),
    CyclicPartsRule = any(vapply(r$parts_rules, mr_contains_dont_know, NA)),
    AlternativeRule = any(vapply(r$alternatives, mr_contains_dont_know, NA)),
    PiecewiseRule = any(vapply(r$subfunctions, function(sf) mr_contains_dont_know(sf[[1L]]), NA)),
    FALSE)
}

mr_eval <- function(r) {
  x <- r$variable; f <- r$integrand
  switch(r$type,
    ConstantRule = .mul(list(f, x)),
    ConstantTimesRule = .mul(list(r$constant, mr_eval(r$substep))),
    PowerRule = Piecewise(list(.div(.pow(r$base, .add(list(r$exp, S_One))), .add(list(r$exp, S_One))), Ne(r$exp, -1)),
                          list(fcall("log", list(r$base)), S_true)),
    NestedPowRule = {
      m <- .mul(list(r$base, f))
      Piecewise(list(.div(m, .add(list(r$exp, S_One))), Ne(r$exp, -1)), list(.mul(list(m, fcall("log", list(r$base)))), S_true))
    },
    AddRule = .add(lapply(r$substeps, mr_eval)),
    URule = {
      result <- mr_eval(r$substep)
      if (r$u_func$type == "Pow") {
        be <- .as_base_exp(r$u_func)
        if (identical(be[[2L]]$key, S_NegativeOne$key)) result <- subs(result, fcall("log", list(r$u_var)), .neg(fcall("log", list(be[[1L]]))))
      }
      subs(result, r$u_var, r$u_func)
    },
    ReparameterizationRule = xreplace(mr_eval(r$substep), r$replacements),
    PartsRule = .sub(.mul(list(r$u, mr_eval(r$v_step))), mr_eval(r$second_step)),
    CyclicPartsRule = {
      res <- list(); sg <- 1L
      for (pr in r$parts_rules) { res[[length(res) + 1L]] <- .mul(list(.int(sg), pr$u, mr_eval(pr$v_step))); sg <- -sg }
      .div(.add(res), .sub(S_One, r$coefficient))
    },
    SinRule = .neg(cos(x)), CosRule = sin(x), SinhRule = cosh(x), CoshRule = sinh(x),
    ExpRule = .div(f, fcall("log", list(r$base))),
    ReciprocalRule = fcall("log", list(r$base)),
    ArcsinRule = asin(x), ArcsinhRule = asinh(x),
    ReciprocalSqrtQuadraticRule = {
      a <- r$a; b <- r$b; cc <- r$c
      q <- .add(list(a, .mul(list(b, x)), .mul(list(cc, .sq(x)))))
      inner <- .add(list(.mul(list(S_Two, .sqrt(cc), .sqrt(q))), b, .mul(list(S_Two, cc, x))))
      .div(fcall("log", list(inner)), .sqrt(cc))
    },
    SqrtQuadraticDenomRule = {
      a <- r$a; b <- r$b; cc <- r$c
      red <- .sqrt_quadratic_denom_reduce(r$coeffs, a, b, cc)
      s <- .sqrt(.add(list(a, .mul(list(b, x)), .mul(list(cc, .sq(x))))))
      I0 <- if (!is.null(r$i0_step)) .mul(list(red[[3L]], mr_eval(r$i0_step))) else S_Zero
      n <- length(r$coeffs); rc <- red[[1L]]
      terms <- c(lapply(seq_along(rc), function(i) .mul(list(rc[[i]], .pow(x, .int(n - 1L - i))))), list(.div(red[[2L]], cc)))
      .add(list(.mul(list(.add(terms), s)), I0))
    },
    SqrtQuadraticRule = , RewriteRule = , CompleteSquareRule = mr_eval(r$substep),
    RatintRule = ratint(f, x),
    AlternativeRule = mr_eval(r$alternatives[[1L]]),
    DontKnowRule = Integral(f, x),
    DerivativeRule = .mr_eval_derivative(r),
    PiecewiseRule = piecewise_fold(do.call(Piecewise, lapply(r$subfunctions, function(sf) list(mr_eval(sf[[1L]]), sf[[2L]])))),
    HeavisideRule = {
      result <- mr_eval(r$substep)
      .mul(list(Heaviside(r$harg), .sub(result, subs(result, x, r$ibnd))))
    },
    DiracDeltaRule = if (is_zero_num(r$n)) .div(Heaviside(.add(list(r$a, .mul(list(r$b, x))))), r$b)
                     else .div(fcall("DiracDelta", list(.add(list(r$a, .mul(list(r$b, x)))), .sub(r$n, S_One))), r$b),
    ArctanRule = .mul(list(.div(.div(r$a, r$b), .sqrt(.div(r$c, r$b))), atan(.div(x, .sqrt(.div(r$c, r$b)))))),
    JacobiRule = {
      n <- r$n; a <- r$a; b <- r$b
      Piecewise(list(.div(.mul(list(S_Two, fcall("jacobi", list(.add(list(n, S_One)), .sub(a, S_One), .sub(b, S_One), x)))), .add(list(n, a, b))), Ne(.add(list(n, a, b)), 0)),
                list(x, Eq(n, 0)),
                list(.add(list(.div(.mul(list(.add(list(a, b, S_Two)), .sq(x))), .int(4)), .div(.mul(list(.sub(a, b), x)), S_Two))), Eq(n, 1)))
    },
    GegenbauerRule = {
      n <- r$n; a <- r$a
      Piecewise(list(.div(fcall("gegenbauer", list(.add(list(n, S_One)), .sub(a, S_One), x)), .mul(list(S_Two, .sub(a, S_One)))), Ne(a, 1)),
                list(.div(fcall("chebyshevt", list(.add(list(n, S_One)), x)), .add(list(n, S_One))), Ne(n, -1)), list(S_Zero, S_true))
    },
    ChebyshevTRule = {
      n <- r$n
      Piecewise(list(.div(.sub(.div(fcall("chebyshevt", list(.add(list(n, S_One)), x)), .add(list(n, S_One))),
                               .div(fcall("chebyshevt", list(.sub(n, S_One), x)), .sub(n, S_One))), S_Two), Ne(Abs(n), 1)),
                list(.div(.sq(x), S_Two), S_true))
    },
    ChebyshevURule = Piecewise(list(.div(fcall("chebyshevt", list(.add(list(r$n, S_One)), x)), .add(list(r$n, S_One))), Ne(r$n, -1)), list(S_Zero, S_true)),
    LegendreRule = .div(.sub(fcall("legendre", list(.add(list(r$n, S_One)), x)), fcall("legendre", list(.sub(r$n, S_One), x))), .add(list(.mul(list(S_Two, r$n)), S_One))),
    HermiteRule = .div(fcall("hermite", list(.add(list(r$n, S_One)), x)), .mul(list(S_Two, .add(list(r$n, S_One))))),
    LaguerreRule = .sub(fcall("laguerre", list(r$n, x)), fcall("laguerre", list(.add(list(r$n, S_One)), x))),
    AssocLaguerreRule = .neg(fcall("assoc_laguerre", list(.add(list(r$n, S_One)), .sub(r$a, S_One), x))),
    CiRule = .sub(.mul(list(cos(r$b), Ci(.mul(list(r$a, x))))), .mul(list(sin(r$b), Si(.mul(list(r$a, x)))))),
    ChiRule = .add(list(.mul(list(cosh(r$b), Chi(.mul(list(r$a, x))))), .mul(list(sinh(r$b), Shi(.mul(list(r$a, x))))))),
    EiRule = .mul(list(exp(r$b), Ei(.mul(list(r$a, x))))),
    SiRule = .add(list(.mul(list(sin(r$b), Ci(.mul(list(r$a, x))))), .mul(list(cos(r$b), Si(.mul(list(r$a, x))))))),
    ShiRule = .add(list(.mul(list(sinh(r$b), Chi(.mul(list(r$a, x))))), .mul(list(cosh(r$b), Shi(.mul(list(r$a, x))))))),
    LiRule = .div(fcall("li", list(.add(list(.mul(list(r$a, x)), r$b)))), r$a),
    ErfRule = .mr_eval_erf(r),
    OwensTRule = .mul(list(.div(.mul(list(.int(-2), .sqrt(S_Pi))), r$a), fcall("owens_t", list(.mul(list(.sqrt(S_Two), .add(list(.mul(list(r$a, x)), r$b)))), r$y)))),
    FresnelCRule = .mr_eval_fresnel(r, TRUE),
    FresnelSRule = .mr_eval_fresnel(r, FALSE),
    PolylogRule = fcall("polylog", list(.add(list(r$b, S_One)), .mul(list(r$a, x)))),
    UpperGammaRule = .div(.mul(list(.pow(x, r$e), .pow(.mul(list(.neg(r$a), x)), .neg(r$e)), fcall("uppergamma", list(.add(list(r$e, S_One)), .mul(list(.neg(r$a), x)))))), r$a),
    EllipticFRule = .div(fcall("elliptic_f", list(x, .div(r$d, r$a))), .sqrt(r$a)),
    EllipticERule = .mul(list(fcall("elliptic_e", list(x, .div(r$d, r$a))), .sqrt(r$a))),
    stop("unknown manualintegrate rule ", r$type))
}

.mr_eval_derivative <- function(r) {
  d <- r$integrand
  vc <- d$args[-1L]
  for (k in seq_along(vc)) {
    t <- vc[[k]]
    v <- if (t$type == "Tuple") t$args[[1L]] else t
    cnt <- if (t$type == "Tuple") t$args[[2L]] else S_One
    if (identical(v$key, r$variable$key)) { vc[[k]] <- Tuple(v, .sub(cnt, S_One)); break }
  }
  # Derivative() drops zero counts and returns the expression when none remain
  vc <- Filter(function(t) !(t$type == "Tuple" && is_zero_num(t$args[[2L]])), vc)
  if (!length(vc)) return(d$args[[1L]])
  new_sym("Derivative", c(list(d$args[[1L]]), vc))
}

.mr_eval_erf <- function(r) {
  a <- r$a; b <- r$b; c <- r$c; x <- r$variable
  common <- function(sq) .mul(list(.div(.div(.sqrt(S_Pi), sq), S_Two), exp(.sub(c, .div(.sq(b), .mul(list(.int(4), a)))))))
  ei <- .mul(list(common(.sqrt(a)), erfi(.div(.add(list(.mul(list(S_Two, a, x)), b)), .mul(list(S_Two, .sqrt(a)))))))
  if (isTRUE(ask_fact(a, "extended_real"))) {
    ef <- .mul(list(common(.sqrt(.neg(a))), erf(.div(.sub(.mul(list(.int(-2), a, x)), b), .mul(list(S_Two, .sqrt(.neg(a))))))))
    return(Piecewise(list(ef, Lt(a, 0)), list(ei, S_true)))
  }
  ei
}

.mr_eval_fresnel <- function(r, is_c) {
  a <- r$a; b <- r$b; c <- r$c; x <- r$variable
  arg <- .div(.add(list(.mul(list(S_Two, a, x)), b)), .sqrt(.mul(list(S_Two, a, S_Pi))))
  ang <- .sub(.div(.sq(b), .mul(list(.int(4), a))), c)
  pre <- .div(.sqrt(S_Pi), .sqrt(.mul(list(S_Two, a))))
  if (is_c) return(.mul(list(pre, .add(list(.mul(list(cos(ang), fresnelc(arg))), .mul(list(sin(ang), fresnels(arg))))))))
  .mul(list(pre, .sub(.mul(list(cos(ang), fresnels(arg))), .mul(list(sin(ang), fresnelc(arg))))))
}

.sqrt_quadratic_denom_reduce <- function(coeffs, a, b, c) {
  result <- list()
  if (length(coeffs) > 2L) for (i in seq_len(length(coeffs) - 2L)) {
    n <- length(coeffs) - i
    coeff <- .div(coeffs[[i]], .mul(list(c, .int(n))))
    result[[length(result) + 1L]] <- coeff
    coeffs[[i + 1L]] <- .sub(coeffs[[i + 1L]], .mul(list(.div(.mul(list(.int(2L * n - 1L), b)), S_Two), coeff)))
    coeffs[[i + 2L]] <- .sub(coeffs[[i + 2L]], .mul(list(.int(n - 1L), a, coeff)))
  }
  d <- coeffs[[length(coeffs)]]; e <- coeffs[[length(coeffs) - 1L]]
  list(result, e, .sub(d, .div(.mul(list(b, e)), .mul(list(S_Two, c)))))
}

# ---------------------------------------------------------------------------
# helpers

.if_zero_implies_zero <- function(P, Q) {
  P <- expand(P); Q <- expand(Q)
  pz <- ask_fact(P, "zero")
  if (identical(pz, FALSE)) return(TRUE)
  if (isTRUE(pz)) return(isTRUE(ask_fact(Q, "zero")))
  num_p <- .as_numer_denom(P)[[1L]]; num_q <- .as_numer_denom(Q)[[1L]]
  fq <- sym_factor(num_q)
  for (fk in factor_list(num_p)$factors) if (!is_zero_num(subs(fq, fk[[1L]], S_Zero))) return(FALSE)
  TRUE
}

manual_diff <- function(f, symbol) {
  f <- as_sym(f)
  if (length(f$args)) {
    a <- f$args[[1L]]
    if (f$type == "Function") {
      if (f$name == "tan") return(.mul(list(.diff1(a, symbol), .sq(fcall("sec", list(a))))))
      if (f$name == "cot") return(.neg(.mul(list(.diff1(a, symbol), .sq(fcall("csc", list(a)))))))
      if (f$name == "sec") return(.mul(list(.diff1(a, symbol), fcall("sec", list(a)), fcall("tan", list(a)))))
      if (f$name == "csc") return(.neg(.mul(list(.diff1(a, symbol), fcall("csc", list(a)), fcall("cot", list(a))))))
    }
    if (f$type == "Add") return(.add(lapply(f$args, manual_diff, symbol = symbol)))
    if (f$type == "Mul" && length(f$args) == 2L && is_number_type(f$args[[1L]])) return(.mul(list(f$args[[1L]], manual_diff(f$args[[2L]], symbol))))
  }
  .diff1(f, symbol)
}

manual_subs <- function(expr, old, new) {
  seq <- list(list(old, new))
  new_subs <- list()
  for (p in seq) {
    o <- p[[1L]]
    if (o$type == "Function" && o$name == "log") {
      x0 <- o$args[[1L]]
      expr <- .map_bottom_up(expr, function(z) if (z$type == "Pow" && identical(z$args[[1L]]$key, x0$key)) fcall("exp", list(.mul(list(z$args[[2L]], p[[2L]])))) else z)
      new_subs[[length(new_subs) + 1L]] <- list(x0, fcall("exp", list(p[[2L]])))
    }
  }
  subs(expr, c(seq, new_subs))
}

.INV_TRIG <- c("atan", "asin", "acos", "acot", "acsc", "asec")

find_substitutions <- function(integrand, symbol, u_var) {
  results <- list()
  test_subterm <- function(u, u_diff) {
    if (is_zero_num(u_diff)) return(NULL)
    substituted <- .div(integrand, u_diff)
    substituted <- cancel(manual_subs(substituted, u, u_var))
    if (.has1(substituted, symbol) && !.is_ratfun(integrand, symbol) && isTRUE(is_polynomial(u, symbol)) && degree(u, symbol) == 1) {
      a <- coeff(u, symbol); b <- subs(u, symbol, S_Zero)
      symbol_of_u <- .div(.sub(u_var, b), a)
      if (!is_zero_num(b) && !is_zero_num(a)) substituted <- cancel(manual_subs(substituted, symbol, symbol_of_u))
    }
    if (.has1(substituted, symbol)) return(NULL)
    if (.is_ratfun(integrand, symbol) && .is_ratfun(substituted, u_var)) {
      db <- max(vapply(.as_numer_denom(integrand), function(t) .ri_deg_in(t, symbol), 0))
      da <- max(vapply(.as_numer_denom(substituted), function(t) .ri_deg_in(t, u_var), 0))
      if (da > db) return(NULL)
    }
    # as_independent(u_var, as_Add = False)
    sp <- .as_coeff_mul_dep(substituted, u_var)
    list(sp[[1L]], .mul(sp[[2L]]))
  }
  exp_subterms <- function(term) {
    lin <- list(); terms <- list()
    n <- Wild("n", properties = list(function(v) v$type == "Integer"))
    for (ex in .uniq_syms(Filter(function(z) z$type == "Function" && z$name == "exp", .all_subexprs(term)))) {
      arg <- ex$args[[1L]]
      if (!.has1(arg, symbol)) next
      m <- tryCatch(sym_match(arg, n * symbol), error = function(e) NULL)
      if (!is.null(m) && length(m)) lin[[length(lin) + 1L]] <- m$n
      else terms[[length(terms) + 1L]] <- ex
    }
    if (length(lin)) terms[[length(terms) + 1L]] <- fcall("exp", list(.mul(list(gcd_list(lin), symbol))))
    terms
  }
  possible_subterms <- function(term) {
    if (term$type == "Function") {
      nm <- term$name
      if (nm %in% c(.TRIG_FUNS, .HYP_FUNS, .INV_TRIG, "exp", "log", "Heaviside")) return(list(term$args[[1L]]))
      if (nm %in% c("chebyshevt", "chebyshevu", "legendre", "hermite", "laguerre")) return(list(term$args[[2L]]))
      if (nm %in% c("gegenbauer", "assoc_laguerre")) return(list(term$args[[3L]]))
      if (nm == "jacobi") return(list(term$args[[4L]]))
      return(list())
    }
    if (term$type == "Mul") {
      r <- list()
      for (u in .sympy_args(term)) { r[[length(r) + 1L]] <- u; r <- c(r, possible_subterms(u)) }
      return(r)
    }
    if (term$type == "Pow") {
      r <- Filter(function(a) .has1(a, symbol), term$args)
      ex <- term$args[[2L]]
      if (ex$type == "Integer") {
        ev <- abs(.int_val(ex))
        pf <- as.integer(as.character(primefactors(.int(.int_val(ex)))))
        for (d in pf) if (1 < d && d < ev) r[[length(r) + 1L]] <- .pow(term$args[[1L]], .int(d))
        if (term$args[[1L]]$type == "Add") r <- c(r, Filter(function(t) t$type == "Pow", possible_subterms(term$args[[1L]])))
      }
      return(r)
    }
    if (term$type == "Add") {
      r <- list()
      for (a in .sympy_args(term)) { r[[length(r) + 1L]] <- a; r <- c(r, possible_subterms(a)) }
      return(r)
    }
    list()
  }
  cands <- .uniq_syms(c(possible_subterms(integrand), exp_subterms(integrand)))
  sub_x <- subs(integrand, symbol, u_var)
  for (u in cands) {
    if (identical(u$key, symbol$key)) next
    u_diff <- manual_diff(u, symbol)
    ni <- test_subterm(u, u_diff)
    if (!is.null(ni)) {
      if (identical(ni[[2L]]$key, sub_x$key)) next
      sb <- list(u, ni[[1L]], ni[[2L]])
      if (!any(vapply(results, function(r) identical(r[[1L]]$key, sb[[1L]]$key) && identical(r[[2L]]$key, sb[[2L]]$key) && identical(r[[3L]]$key, sb[[3L]]$key), NA)))
        results[[length(results) + 1L]] <- sb
    }
  }
  if (length(results) > 1L) results <- results[order(vapply(results, function(s) count_ops(s[[3L]]), 0))]
  results
}

# rewriter(condition, rewrite) and proxy_rewriter
.mi_rewriter <- function(condition, rewrite_fn) function(integral, S) {
  integrand <- integral[[1L]]; symbol <- integral[[2L]]
  if (condition(integrand, symbol)) {
    rewritten <- rewrite_fn(integrand, symbol)
    if (!identical(rewritten$key, integrand$key)) {
      substep <- S$solve(rewritten, symbol)
      if (!.is_mr(substep, "DontKnowRule") && !is.null(substep)) return(.mr("RewriteRule", integrand, symbol, rewritten = rewritten, substep = substep))
    }
  }
  NULL
}
.mi_proxy_rewriter <- function(condition, rewrite_fn) function(criteria, integral, S) {
  integrand <- integral[[1L]]; symbol <- integral[[2L]]
  args <- c(criteria, list(integrand, symbol))
  if (do.call(condition, args)) {
    rewritten <- do.call(rewrite_fn, args)
    if (!identical(rewritten$key, integrand$key)) {
      substep <- S$solve(rewritten, symbol)
      return(.mr("RewriteRule", integrand, symbol, rewritten = rewritten, substep = substep))
    }
  }
  NULL
}

.mr_equal <- function(a, b) identical(.mr_key(a), .mr_key(b))
.mr_key <- function(r) {
  if (inherits(r, "mrule")) return(paste0(r$type, "(", paste(vapply(r[-1L], .mr_key, ""), collapse = ","), ")"))
  if (is_sym(r)) return(r$key)
  if (is.list(r)) return(paste0("[", paste(vapply(r, .mr_key, ""), collapse = ","), "]"))
  if (is.null(r)) return("NULL")
  paste(format(r), collapse = ",")
}

.mi_alternatives <- function(rules, branch = FALSE) function(integral, S) {
  alts <- list(); first_partial <- NULL
  for (rule in rules) {
    result <- rule(integral, S)
    if (is.null(result) || .is_mr(result, "DontKnowRule")) next
    if (!branch) {
      if (!mr_contains_dont_know(result)) return(result)
      if (is.null(first_partial)) first_partial <- result
      next
    }
    if (!any(vapply(alts, function(a) .mr_equal(a, result), NA))) alts[[length(alts) + 1L]] <- result
  }
  if (!branch) return(first_partial)
  if (length(alts) == 1L) return(alts[[1L]])
  if (length(alts)) {
    doable <- Filter(function(r) !mr_contains_dont_know(r), alts)
    return(.mr("AlternativeRule", integral[[1L]], integral[[2L]], alternatives = if (length(doable)) doable else alts))
  }
  NULL
}

constant_rule <- function(integral, S) .mr("ConstantRule", integral[[1L]], integral[[2L]])

power_rule <- function(integral, S) {
  integrand <- integral[[1L]]; symbol <- integral[[2L]]
  be <- .as_base_exp(integrand); base <- be[[1L]]; expt <- be[[2L]]
  if (!.has1(expt, symbol) && base$type %in% c("Symbol", "Dummy")) {
    if (is_zero_num(simplify(.add(list(expt, S_One))))) return(.mr("ReciprocalRule", integrand, symbol, base = base))
    return(.mr("PowerRule", integrand, symbol, base = base, exp = expt))
  }
  if (!.has1(base, symbol) && expt$type %in% c("Symbol", "Dummy")) {
    rule <- .mr("ExpRule", integrand, symbol, base = base, exp = expt)
    lz <- ask_fact(fcall("log", list(base)), "zero")
    if (identical(lz, FALSE)) return(rule)
    if (isTRUE(lz)) return(.mr("ConstantRule", S_One, symbol))
    return(.mr("PiecewiseRule", integrand, symbol, subfunctions = list(list(rule, Ne(fcall("log", list(base)), 0)), list(.mr("ConstantRule", S_One, symbol), S_true))))
  }
  NULL
}

exp_rule <- function(integral, S) {
  integrand <- integral[[1L]]
  if (integrand$args[[1L]]$type %in% c("Symbol", "Dummy")) return(.mr("ExpRule", integrand, integral[[2L]], base = S_E, exp = integrand$args[[1L]]))
  NULL
}

combine_power_rule <- function(integral, S) {
  integrand <- integral[[1L]]; symbol <- integral[[2L]]
  a <- Wild("a", exclude = list(symbol, S_Zero)); b <- Wild("b", exclude = list(symbol)); k <- Wild("k", exclude = list(symbol)); rest <- Wild("rest")
  m <- tryCatch(sym_match(integrand, rest * k^(a * symbol^2) * k^(b * symbol)), error = function(e) NULL)
  if (is.null(m) || !length(m)) return(NULL)
  simplified <- powsimp(integrand, combine = "exp")
  if (!identical(simplified$key, integrand$key)) return(.mr("RewriteRule", integrand, symbol, rewritten = simplified, substep = S$solve(simplified, symbol)))
  NULL
}

orthogonal_poly_rule <- function(integral, S) {
  integrand <- integral[[1L]]; symbol <- integral[[2L]]
  cls <- c(jacobi = "JacobiRule", gegenbauer = "GegenbauerRule", chebyshevt = "ChebyshevTRule", chebyshevu = "ChebyshevURule",
           legendre = "LegendreRule", hermite = "HermiteRule", laguerre = "LaguerreRule", assoc_laguerre = "AssocLaguerreRule")
  vi <- c(jacobi = 4L, gegenbauer = 3L, assoc_laguerre = 3L)
  if (integrand$type != "Function" || !(integrand$name %in% names(cls))) return(NULL)
  k <- integrand$name
  idx <- if (is.na(vi[k])) 2L else vi[[k]]
  if (identical(integrand$args[[idx]]$key, symbol$key) && !any(vapply(integrand$args[seq_len(idx - 1L)], function(v) .has1(v, symbol), NA))) {
    pars <- integrand$args[seq_len(idx - 1L)]
    nm <- c("n", "a", "b")[seq_along(pars)]
    return(do.call(.mr, c(list(cls[[k]], integrand, symbol), stats::setNames(pars, nm))))
  }
  NULL
}

.add_degenerate_step <- function(generic_cond, generic_step, degenerate_step) {
  if (is.null(degenerate_step)) return(generic_step)
  if (.is_mr(generic_step, "PiecewiseRule")) subfunctions <- lapply(generic_step$subfunctions, function(sf) list(sf[[1L]], simplify(And(sf[[2L]], generic_cond))))
  else subfunctions <- list(list(generic_step, generic_cond))
  if (.is_mr(degenerate_step, "PiecewiseRule")) subfunctions <- c(subfunctions, degenerate_step$subfunctions)
  else subfunctions[[length(subfunctions) + 1L]] <- list(degenerate_step, S_true)
  .mr("PiecewiseRule", generic_step$integrand, generic_step$variable, subfunctions = subfunctions)
}

.mi_sf_state <- new.env()
special_function_rule <- function(integral, S) {
  integrand <- integral[[1L]]; symbol <- integral[[2L]]
  if (is.null(.mi_sf_state$patterns)) {
    xs <- Dummy("x"); .mi_sf_state$x <- xs
    a <- Wild("a", exclude = list(xs), properties = list(function(v) !isTRUE(ask_fact(v, "zero"))))
    b <- Wild("b", exclude = list(xs)); c <- Wild("c", exclude = list(xs))
    d <- Wild("d", exclude = list(xs), properties = list(function(v) !isTRUE(ask_fact(v, "zero"))))
    e <- Wild("e", exclude = list(xs), properties = list(function(v) !(isTRUE(ask_fact(v, "nonnegative")) && isTRUE(ask_fact(v, "integer")))))
    y <- Wild("y", exclude = list(xs))
    .mi_sf_state$wilds <- list(a = a, b = b, c = c, d = d, e = e, y = y)
    lin <- a * xs + b; quad <- a * xs^2 + b * xs + c
    fn <- function(n, ...) fnode(n, lapply(list(...), as_sym))
    .mi_sf_state$patterns <- list(
      list("Mul", fn("exp", lin) / xs, NULL, "EiRule"),
      list("Mul", fn("cos", lin) / xs, NULL, "CiRule"),
      list("Mul", fn("cosh", lin) / xs, NULL, "ChiRule"),
      list("Mul", fn("sin", lin) / xs, NULL, "SiRule"),
      list("Mul", fn("sinh", lin) / xs, NULL, "ShiRule"),
      list("Pow", 1 / fn("log", lin), NULL, "LiRule"),
      list("exp", fn("exp", quad), NULL, "ErfRule"),
      list("Mul", fn("exp", -lin^2) * erf(y * lin), function(w) !is_one_num(w$y) && !identical(w$y$key, S_NegativeOne$key), "OwensTRule"),
      list("sin", fn("sin", quad), NULL, "FresnelSRule"),
      list("cos", fn("cos", quad), NULL, "FresnelCRule"),
      list("Mul", xs^e * fn("exp", a * xs), NULL, "UpperGammaRule"),
      list("Mul", fn("polylog", b, a * xs) / xs, NULL, "PolylogRule"),
      list("Pow", 1 / sqrt(a - d * fn("sin", xs)^2), function(w) !identical(w$a$key, w$d$key), "EllipticFRule"),
      list("Pow", sqrt(a - d * fn("sin", xs)^2), function(w) !identical(w$a$key, w$d$key), "EllipticERule"))
  }
  xs <- .mi_sf_state$x; W <- .mi_sf_state$wilds
  ti <- subs(integrand, symbol, xs)
  tkey <- if (ti$type == "Function") ti$name else ti$type
  slots <- list(EiRule = c("a", "b"), CiRule = c("a", "b"), ChiRule = c("a", "b"), SiRule = c("a", "b"), ShiRule = c("a", "b"), LiRule = c("a", "b"),
                ErfRule = c("a", "b", "c"), OwensTRule = c("a", "b", "y"), FresnelSRule = c("a", "b", "c"), FresnelCRule = c("a", "b", "c"),
                UpperGammaRule = c("a", "e"), PolylogRule = c("a", "b"), EllipticFRule = c("a", "d"), EllipticERule = c("a", "d"))
  for (pt in .mi_sf_state$patterns) {
    if (tkey != pt[[1L]]) next
    m <- tryCatch(sym_match(ti, pt[[2L]]), error = function(e) NULL)
    if (is.null(m) || !length(m)) next
    if (!is.null(pt[[3L]]) && !pt[[3L]](m)) next
    rname <- pt[[4L]]
    vals <- lapply(slots[[rname]], function(k) m[[k]])
    names(vals) <- slots[[rname]]
    if (any(vapply(vals, is.null, NA))) next
    step <- do.call(.mr, c(list(rname, integrand, symbol), vals))
    for (w in c("a", "d")) {
      val <- m[[w]]
      if (!is.null(val) && !identical(ask_fact(val, "zero"), FALSE)) {
        dstep <- S$solve(subs(integrand, val, S_Zero), symbol)
        step <- .add_degenerate_step(Ne(val, 0), step, dstep)
      }
    }
    if (rname %in% c("EllipticFRule", "EllipticERule")) {
      av <- m$a; dv <- m$d
      if (!is.null(av) && !is.null(dv) && !identical(ask_fact(.sub(av, dv), "zero"), FALSE)) {
        dstep <- S$solve(subs(integrand, dv, av), symbol)
        step <- .add_degenerate_step(Ne(av, dv), step, dstep)
      }
    }
    return(step)
  }
  NULL
}

nested_pow_rule <- function(integral, S) {
  integrand <- integral[[1L]]; x <- integral[[2L]]
  a_ <- Wild("a", exclude = list(x)); b_ <- Wild("b", exclude = list(x, S_Zero))
  generic_cond <- S_true
  nomatch <- function() stop(structure(class = c("NoMatch", "error", "condition"), list(message = "nomatch", call = NULL)))
  gbe <- function(expr) {
    if (!.has1(expr, x)) return(list(S_One, S_Zero))
    if (expr$type == "Mul") {
      terms <- as_coeff_mul_split(expr)
      tl <- if (is_one_num(terms[[1L]]) && FALSE) list() else (if (terms[[2L]]$type == "Mul") .sympy_args(terms[[2L]]) else list(terms[[2L]]))
      if (!length(tl)) return(list(S_One, S_Zero))
      res <- lapply(tl, gbe)
      bases <- .uniq_syms(Filter(function(b) !is_one_num(b), lapply(res, `[[`, 1L)))
      if (length(bases) == 1L) return(list(bases[[1L]], .add(lapply(res, `[[`, 2L))))
      nomatch()
    }
    if (expr$type == "Pow") {
      if (.has1(expr$args[[2L]], x)) nomatch()
      r <- gbe(expr$args[[1L]])
      return(list(r[[1L]], .mul(list(r[[2L]], expr$args[[2L]]))))
    }
    m <- tryCatch(sym_match(expr, a_ + b_ * x), error = function(e) NULL)
    if (!is.null(m) && length(m) && !is.null(m$b)) {
      av <- if (is.null(m$a)) S_Zero else m$a
      generic_cond <<- Ne(m$b, 0)
      return(list(.add(list(x, .div(av, m$b))), S_One))
    }
    nomatch()
  }
  be <- tryCatch(gbe(integrand), NoMatch = function(e) NULL)
  if (is.null(be)) return(NULL)
  degenerate_step <- if (.is_true(generic_cond)) NULL else .mr("ConstantRule", subs(integrand, x, S_Zero), x)
  .add_degenerate_step(generic_cond, .mr("NestedPowRule", integrand, x, base = be[[1L]], exp = be[[2L]]), degenerate_step)
}

add_rule <- function(integral, S) {
  integrand <- integral[[1L]]; symbol <- integral[[2L]]
  results <- lapply(.as_ordered_terms(integrand), function(g) S$solve(g, symbol))
  if (any(vapply(results, is.null, NA))) return(NULL)
  .mr("AddRule", integrand, symbol, substeps = results)
}
# Add.as_ordered_terms: the printing order of the terms
.as_ordered_terms <- function(e) if (e$type == "Add") .add_print_order(e$args) else list(e)

mul_rule <- function(integral, S) {
  integrand <- integral[[1L]]; symbol <- integral[[2L]]
  sp <- .as_coeff_mul_dep(integrand, symbol)
  coeff <- sp[[1L]]; f <- .mul(sp[[2L]])
  if (!is_one_num(coeff)) {
    nxt <- S$solve(f, symbol)
    if (!is.null(nxt)) return(.mr("ConstantTimesRule", integrand, symbol, constant = coeff, other = f, substep = nxt))
  }
  NULL
}
