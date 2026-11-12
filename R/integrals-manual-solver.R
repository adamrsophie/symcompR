# manualintegrate, after sympy/integrals/manualintegrate.py (part 2): the
# remaining rule functions, the IntegrationSolver with its rule strategy,
# integral_steps and manualintegrate.

.mi_fn_is <- function(e, names) e$type == "Function" && e$name %in% names
.mi_args <- function(e) if (length(e$args)) .sympy_args(e) else list()
.mi_make_args_mul <- function(e) if (e$type == "Mul") .sympy_args(e) else list(e)
.mi_make_args_add <- function(e) if (e$type == "Add") .sympy_args(e) else list(e)
.mi_has_fn <- function(e, names) .has_function(e, names)
.mi_eq <- function(a, b) identical(as_sym(a)$key, as_sym(b)$key)
# Expr.equals: structural or symbolic equality
.mi_equals <- function(a, b) .mi_eq(a, b) || isTRUE(tryCatch(is_zero_num(simplify(.sub(a, b))), error = function(e) FALSE))
.mi_has_infinite <- function(e) .has1(e, S_ComplexInfinity) || .has1(e, S_Infinity) || .has1(e, S_NegativeInfinity) || .has1(e, S_NaN)
.mi_is_integer_num <- function(e) e$type == "Integer"
.mi_is_rational_num <- function(e) e$type %in% c("Integer", "Rational")
.mi_q <- function(e) if (e$type == "Integer") 1L else as.integer(as.character(gmp::denominator(e$value)))
.mi_p <- function(e) if (e$type == "Integer") .int_val(e) else as.integer(as.character(gmp::numerator(e$value)))
.mi_as_coeff_Mul <- function(e) as_coeff_mul_split(e)
.mi_rewrite_only <- function(e, names, target) {
  rules <- .rewrite_rules[[target]]
  .map_bottom_up(e, function(z) {
    if (z$type == "Function" && z$name %in% names && !is.null(rules[[z$name]])) {
      r <- do.call(rules[[z$name]], z$args)
      if (!is.null(r)) return(r)
    }
    z
  })
}
# number of subexpressions that are the function `name` (Basic.count)
.mi_count_fn <- function(e, name) {
  n <- 0L
  rec <- function(z) { if (z$type == "Function" && z$name == name) n <<- n + 1L; for (a in z$args) rec(a) }
  rec(e); n
}
# Expr.is_polynomial(*gens) and is_rational_function(*gens), structural as in
# SymPy (no expansion); gens may be function atoms
.mi_struct_check <- function(e, gens, rational) {
  keys <- vapply(gens, `[[`, "", "key")
  hasg <- function(z) any(vapply(gens, function(g) .has1(z, g), NA))
  rec <- function(z) {
    if (z$key %in% keys) return(TRUE)
    if (!hasg(z)) return(TRUE)
    if (z$type %in% c("Add", "Mul")) return(all(vapply(z$args, rec, NA)))
    if (z$type == "Pow") {
      if (hasg(z$args[[2L]])) return(FALSE)
      ex <- z$args[[2L]]
      ok <- if (rational) ex$type == "Integer" else ex$type == "Integer" && .int_val(ex) >= 0
      return(ok && rec(z$args[[1L]]))
    }
    FALSE
  }
  rec(e)
}
.mi_is_poly <- function(e, gens) .mi_struct_check(e, if (is_sym(gens)) list(gens) else gens, FALSE)
.mi_is_ratfun <- function(e, gens) .mi_struct_check(e, if (is_sym(gens)) list(gens) else gens, TRUE)
# Expr.is_algebraic_expr(symbol)
.mi_is_algebraic <- function(e, x) {
  if (!.has1(e, x) || identical(e$key, x$key)) return(TRUE)
  if (e$type %in% c("Add", "Mul")) return(all(vapply(e$args, .mi_is_algebraic, NA, x = x)))
  if (e$type == "Pow") {
    if (.has1(e$args[[2L]], x)) return(FALSE)
    return(.mi_is_algebraic(e$args[[1L]], x) && .mi_is_rational_num(e$args[[2L]]))
  }
  FALSE
}
# Expr.as_powers_dict()[x]
.mi_power_of <- function(e, x) {
  tot <- NULL
  for (t in .mi_make_args_mul(e)) {
    be <- .as_base_exp(t)
    if (identical(be[[1L]]$key, x$key)) tot <- if (is.null(tot)) be[[2L]] else .add(list(tot, be[[2L]]))
  }
  tot
}
.mi_poly_deg <- function(e, x) { p <- P_or_null(e, x); if (is.null(p)) NA_integer_ else P_deg(p) }
.mi_sort_ordered <- function(l) if (length(l) > 1L) ordered(l) else l
.mi_find_pows <- function(e) .uniq_syms(Filter(function(z) z$type == "Pow", .all_subexprs(e)))
.mi_find_fns <- function(e, names) .uniq_syms(Filter(function(z) .mi_fn_is(z, names), .all_subexprs(e)))

# ---------------------------------------------------------------------------
# parts

.MI_SPECIAL_ERR <- c("erf", "erfc", "erfi", "fresnelc", "fresnels", "Ci", "Chi", "Si", "Shi", "Ei", "li", "owens_t")
.MI_ORTHO <- c("jacobi", "gegenbauer", "chebyshevt", "chebyshevu", "legendre", "assoc_legendre", "hermite", "hermite_prob", "laguerre", "assoc_laguerre")

.parts_rule_gen <- function(integrand, symbol, S) {
  pull_out_algebraic <- function(ig) {
    ig <- together(cancel(ig))
    algebraic <- if (ig$type == "Piecewise" || ig$type != "Mul") list() else Filter(function(a) .mi_is_algebraic(a, symbol), .sympy_args(ig))
    if (length(algebraic)) { u <- .mul(algebraic); return(list(u, cancel(.div(ig, u)))) }
    NULL
  }
  pull_out_u <- function(fns) function(ig) {
    if (.mi_has_fn(ig, fns)) {
      args <- Filter(function(a) .mi_fn_is(a, fns), .mi_args(ig))
      if (length(args)) { u <- .mul(args); return(list(u, .div(ig, u))) }
    }
    NULL
  }
  pull_out_dv <- function(fns) function(ig) {
    power <- .mi_power_of(ig, symbol)
    if (!is.null(power) && power$type == "Integer" && .int_val(power) >= 2 && .mi_has_fn(ig, fns)) {
      for (target in .mi_args(ig)) {
        if (!.mi_fn_is(target, fns)) next
        inner <- target$args[[1L]]
        if (isTRUE(is_polynomial(inner, symbol)) && identical(.mi_poly_deg(inner, symbol), 2L)) {
          dv <- .mul(list(target, symbol)); u <- .div(ig, dv)
          if (isTRUE(is_polynomial(u, symbol))) return(list(u, dv))
        }
      }
    }
    NULL
  }
  liate <- list(pull_out_dv("exp"), pull_out_dv(c("sin", "cos")), pull_out_u(.MI_SPECIAL_ERR), pull_out_u("log"),
                pull_out_u(.INV_TRIG), pull_out_algebraic, pull_out_u(c("sin", "cos")), pull_out_u(c("sinh", "cosh")), pull_out_u("exp"))
  ALG <- 6L
  dummy <- Dummy("temporary")
  if (.mi_fn_is(integrand, c("log", .INV_TRIG, .MI_SPECIAL_ERR))) integrand <- .mul(list(dummy, integrand))
  for (index in seq_along(liate)) {
    result <- liate[[index]](integrand)
    if (is.null(result)) next
    u <- result[[1L]]; dv <- result[[2L]]
    if (!.has1(u, symbol) && !.has1(u, dummy)) return(NULL)
    u <- subs(u, dummy, S_One); dv <- subs(dv, dummy, S_One)
    if (index == ALG && !isTRUE(is_polynomial(u, symbol))) return(NULL)
    if (.mi_fn_is(u, c("log", .MI_SPECIAL_ERR))) {
      rec_dv <- .inv(dv)
      if (isTRUE(is_polynomial(rec_dv, symbol)) && identical(.mi_poly_deg(rec_dv, symbol), 1L)) return(NULL)
    }
    if (index == ALG) {
      if (dv$type == "Derivative" || .mi_has_fn(dv, c(.TRIG_FUNS, .HYP_FUNS)) || .mi_fn_is(dv, .MI_ORTHO)) {
        v_step <- S$solve(dv, symbol)
        if (mr_contains_dont_know(v_step)) return(NULL)
        return(list(u, dv, mr_eval(v_step), .diff1(u, symbol), v_step))
      }
    }
    accept <- FALSE
    if (index < ALG) accept <- TRUE
    else if (index == ALG && length(dv$args) && all(vapply(dv$args, function(a) .mi_fn_is(a, c("sin", "cos", "exp")), NA))) accept <- TRUE
    else if (index < length(liate)) {
      for (lr in liate[(index + 1L):length(liate)]) {
        r <- lr(integrand)
        if (!is.null(r) && .mi_equals(subs(r[[1L]], dummy, S_One), dv)) { accept <- TRUE; break }
      }
    }
    if (accept) {
      du <- .diff1(u, symbol)
      v_step <- S$solve(simplify(dv), symbol)
      if (!mr_contains_dont_know(v_step)) return(list(u, dv, mr_eval(v_step), du, v_step))
    }
  }
  NULL
}

parts_rule <- function(integral, S) {
  symbol <- integral[[2L]]
  cm <- .mi_as_coeff_Mul(integral[[1L]])
  constant <- cm[[1L]]; integrand <- cm[[2L]]
  result <- .parts_rule_gen(integrand, symbol, S)
  steps <- list()
  if (!is.null(result)) {
    u <- result[[1L]]; dv <- result[[2L]]; v <- result[[3L]]; du <- result[[4L]]
    steps[[1L]] <- result
    if (v$type == "Integral") return(NULL)
    v_du <- .mul(list(v, du))
    if (.mi_has_fn(integrand, .TRIG_FUNS) && .mi_has_fn(integrand, "log") && .mi_has_fn(v_du, "log")) {
      if (.mi_count_fn(v_du, "log") > .mi_count_fn(integrand, "log")) return(NULL)
    }
    if (.mi_fn_is(u, c("sin", "cos", "exp", "sinh", "cosh"))) {
      ck <- xreplace(u, list(list(symbol, .mi_cache_dummy())))$key
      if (!S$check_parts_u(ck)) return(NULL)
    }
    for (it in 1:4) {
      if (is_one_num(dv)) break
      coefficient <- cancel(.div(.mul(list(v, du)), integrand))
      if (is_one_num(coefficient)) break
      if (!.has1(coefficient, symbol)) {
        prs <- lapply(steps, function(st) .mr("PartsRule", NULL, NULL, u = st[[1L]], dv = st[[2L]], v_step = st[[5L]], second_step = NULL))
        rule <- .mr("CyclicPartsRule", integrand, symbol, parts_rules = prs, coefficient = .mul(list(.int((-1L)^length(steps)), coefficient)))
        if (!is_one_num(constant)) rule <- .mr("ConstantTimesRule", .mul(list(constant, integrand)), symbol, constant = constant, other = integrand, substep = rule)
        return(rule)
      }
      nc <- .mi_as_coeff_Mul(.mul(list(v, du)))
      result <- .parts_rule_gen(nc[[2L]], symbol, S)
      if (!is.null(result)) {
        u <- .mul(list(result[[1L]], nc[[1L]])); dv <- result[[2L]]; v <- result[[3L]]; du <- .mul(list(result[[4L]], nc[[1L]]))
        steps[[length(steps) + 1L]] <- list(u, dv, v, du, result[[5L]])
      } else break
    }
  }
  make_second_step <- function(steps, ig) {
    if (length(steps)) {
      st <- steps[[1L]]
      substep <- make_second_step(steps[-1L], .mul(list(st[[3L]], st[[4L]])))
      return(.mr("PartsRule", ig, symbol, u = st[[1L]], dv = st[[2L]], v_step = st[[5L]], second_step = substep))
    }
    S$solve(ig, symbol)
  }
  if (length(steps)) {
    st <- steps[[1L]]
    second <- make_second_step(steps[-1L], .mul(list(st[[3L]], st[[4L]])))
    rule <- .mr("PartsRule", integrand, symbol, u = st[[1L]], dv = st[[2L]], v_step = st[[5L]], second_step = second)
    if (!is_one_num(constant)) rule <- .mr("ConstantTimesRule", .mul(list(constant, integrand)), symbol, constant = constant, other = integrand, substep = rule)
    return(rule)
  }
  NULL
}

# ---------------------------------------------------------------------------
# trig

.mi_trig_base_rule <- function(integral, S) {
  original <- integral[[1L]]; symbol <- integral[[2L]]
  integrand <- original
  be <- .as_base_exp(integrand)
  recip <- c(sin = "csc", cos = "sec", tan = "cot", cot = "tan", sec = "cos", csc = "sin")
  if (identical(be[[2L]]$key, S_NegativeOne$key) && .mi_fn_is(be[[1L]], names(recip))) integrand <- fcall(recip[[be[[1L]]$name]], list(symbol))
  f1 <- function(n) fcall(n, list(symbol))
  if (.mi_eq(integrand, sin(symbol))) step <- .mr("SinRule", integrand, symbol)
  else if (.mi_eq(integrand, cos(symbol))) step <- .mr("CosRule", integrand, symbol)
  else {
    u_var <- Dummy("u")
    if (.mi_eq(integrand, f1("tan"))) { rewritten <- .div(sin(symbol), cos(symbol)); u_func <- cos(symbol); substituted <- .neg(.inv(u_var)) }
    else if (.mi_eq(integrand, f1("cot"))) { rewritten <- .div(cos(symbol), sin(symbol)); u_func <- sin(symbol); substituted <- .inv(u_var) }
    else if (.mi_eq(integrand, f1("sec"))) {
      se <- f1("sec"); ta <- f1("tan")
      rewritten <- .div(.add(list(.sq(se), .mul(list(ta, se)))), .add(list(se, ta))); u_func <- .add(list(se, ta)); substituted <- .inv(u_var)
    } else if (.mi_eq(integrand, f1("csc"))) {
      cs <- f1("csc"); co <- f1("cot")
      rewritten <- .div(.add(list(.sq(cs), .mul(list(co, cs)))), .add(list(cs, co))); u_func <- .add(list(cs, co)); substituted <- .neg(.inv(u_var))
    } else return(NULL)
    step <- .mi_trig_substitution_rule(list(integrand, symbol), rewritten, S_One, u_var, u_func, substituted, S)
  }
  if (!.mi_eq(integrand, original)) return(.mr("RewriteRule", original, symbol, rewritten = integrand, substep = step))
  step
}

trig_rule <- function(integral, S) {
  integrand <- integral[[1L]]; symbol <- integral[[2L]]
  be <- .as_base_exp(integrand)
  if (.mi_fn_is(integrand, .TRIG_FUNS)) argument <- integrand$args[[1L]]
  else if (identical(be[[2L]]$key, S_NegativeOne$key) && .mi_fn_is(be[[1L]], .TRIG_FUNS)) argument <- be[[1L]]$args[[1L]]
  else return(NULL)
  p <- P_or_null(argument, symbol)
  if (is.null(p) || P_deg(p) != 1L) return(NULL)
  if (.mi_eq(argument, symbol)) return(.mi_trig_base_rule(integral, S))
  coefficient <- P_nth(p, 1L)
  u_var <- Dummy("u")
  standard <- xreplace(integrand, list(list(argument, u_var)))
  substep <- .mi_trig_base_rule(list(standard, u_var), S)
  if (!is_one_num(coefficient)) {
    constant <- .inv(coefficient)
    substep <- .mr("ConstantTimesRule", .mul(list(constant, standard)), u_var, constant = constant, other = standard, substep = substep)
  }
  generic <- .mr("URule", integrand, symbol, u_var = u_var, u_func = argument, substep = substep)
  .mi_add_trig_degenerate_step(integral, coefficient, generic)
}

.mi_add_trig_degenerate_step <- function(integral, coefficient, generic_step) {
  integrand <- integral[[1L]]; symbol <- integral[[2L]]
  cond <- Ne(coefficient, 0)
  if (.is_true(cond)) return(generic_step)
  di <- subs(subs(integrand, symbol, S_Zero), coefficient, S_Zero)
  if (.mi_has_infinite(di)) return(generic_step)
  dstep <- .mr("RewriteRule", integrand, symbol, rewritten = di, substep = .mr("ConstantRule", di, symbol))
  .add_degenerate_step(cond, generic_step, dstep)
}

.mi_trig_substitution_rule <- function(integral, rewritten, coefficient, u_var, u_func, substituted, S) {
  integrand <- integral[[1L]]; symbol <- integral[[2L]]
  substep <- S$solve(substituted, u_var)
  generic <- .mr("URule", rewritten, symbol, u_var = u_var, u_func = u_func, substep = substep)
  if (!.mi_eq(rewritten, integrand)) generic <- .mr("RewriteRule", integrand, symbol, rewritten = rewritten, substep = generic)
  .mi_add_trig_degenerate_step(integral, coefficient, generic)
}

.mi_collect_trig_powers <- function(integrand, symbol, fns) {
  powers <- stats::setNames(lapply(fns, function(f) list()), fns)
  coefficients <- list()
  for (factor in .mi_make_args_mul(integrand)) {
    be <- .as_base_exp(factor)
    base <- be[[1L]]; power <- be[[2L]]
    if (!.mi_fn_is(base, fns)) return(NULL)
    if (.has1(power, symbol)) return(NULL)
    argument <- base$args[[1L]]
    p <- P_or_null(argument, symbol)
    if (is.null(p) || P_deg(p) != 1L) return(NULL)
    coefficients[[argument$key]] <- list(argument, P_nth(p, 1L))
    fp <- powers[[base$name]]
    old <- fp[[argument$key]]
    fp[[argument$key]] <- if (is.null(old)) power else .add(list(old, power))
    powers[[base$name]] <- fp
  }
  list(powers, coefficients)
}
.mi_pow_get <- function(powers, fn, key) { v <- powers[[fn]][[key]]; if (is.null(v)) S_Zero else v }

.mi_rewrite_to_reciprocal <- function(integrand, fn, recip) {
  .mul(lapply(.mi_make_args_mul(integrand), function(factor) {
    be <- .as_base_exp(factor)
    if (.mi_fn_is(be[[1L]], fn) && isTRUE(ask_fact(be[[2L]], "integer")) && isTRUE(ask_fact(be[[2L]], "negative")))
      return(.pow(fcall(recip, list(be[[1L]]$args[[1L]])), .neg(be[[2L]])))
    factor
  }))
}

.mi_odd <- function(e) isTRUE(ask_fact(e, "odd"))
.mi_even <- function(e) isTRUE(ask_fact(e, "even"))
.mi_neg <- function(e) isTRUE(ask_fact(e, "negative"))

.mi_sincos_sinodd <- function(integral, argument, coefficient, m, n, S) {
  half <- .div(.sub(m, S_One), S_Two)
  rewritten <- .mul(list(.pow(.sub(S_One, .sq(cos(argument))), half), sin(argument), .pow(cos(argument), n)))
  u <- Dummy("u")
  substituted <- .neg(.div(.mul(list(.pow(.sub(S_One, .sq(u)), half), .pow(u, n))), coefficient))
  .mi_trig_substitution_rule(integral, rewritten, coefficient, u, cos(argument), substituted, S)
}
.mi_sincos_cosodd <- function(integral, argument, coefficient, m, n, S) {
  half <- .div(.sub(n, S_One), S_Two)
  rewritten <- .mul(list(.pow(.sub(S_One, .sq(sin(argument))), half), cos(argument), .pow(sin(argument), m)))
  u <- Dummy("u")
  substituted <- .div(.mul(list(.pow(.sub(S_One, .sq(u)), half), .pow(u, m))), coefficient)
  .mi_trig_substitution_rule(integral, rewritten, coefficient, u, sin(argument), substituted, S)
}
.mi_sincos_product_to_sum <- function(integral, S) {
  integrand <- integral[[1L]]; symbol <- integral[[2L]]
  rewritten <- sincos_to_sum(integrand)
  if (!.mi_eq(rewritten, integrand)) return(.mr("RewriteRule", integrand, symbol, rewritten = rewritten, substep = S$solve(rewritten, symbol)))
  NULL
}
.mi_recip_pair <- function(t, s) list(f = function(a) fcall(t, list(a)), g = function(a) fcall(s, list(a)))
.mi_tansec_tanodd <- function(integral, argument, coefficient, m, n, S) {
  ta <- fcall("tan", list(argument)); se <- fcall("sec", list(argument))
  half <- .div(.sub(m, S_One), S_Two)
  rewritten <- .mul(list(.pow(.sub(.sq(se), S_One), half), ta, .pow(se, n)))
  u <- Dummy("u")
  substituted <- .div(.mul(list(.pow(.sub(.sq(u), S_One), half), .pow(u, .sub(n, S_One)))), coefficient)
  .mi_trig_substitution_rule(integral, rewritten, coefficient, u, se, substituted, S)
}
.mi_tansec_seceven <- function(integral, argument, coefficient, m, n, S) {
  ta <- fcall("tan", list(argument)); se <- fcall("sec", list(argument))
  half <- .div(.sub(n, S_Two), S_Two)
  rewritten <- .mul(list(.pow(.add(list(S_One, .sq(ta))), half), .sq(se), .pow(ta, m)))
  u <- Dummy("u")
  substituted <- .div(.mul(list(.pow(.add(list(S_One, .sq(u))), half), .pow(u, m))), coefficient)
  .mi_trig_substitution_rule(integral, rewritten, coefficient, u, ta, substituted, S)
}
.mi_cotcsc_cotodd <- function(integral, argument, coefficient, m, n, S) {
  co <- fcall("cot", list(argument)); cs <- fcall("csc", list(argument))
  half <- .div(.sub(m, S_One), S_Two)
  rewritten <- .mul(list(.pow(.sub(.sq(cs), S_One), half), co, .pow(cs, n)))
  u <- Dummy("u")
  substituted <- .neg(.div(.mul(list(.pow(.sub(.sq(u), S_One), half), .pow(u, .sub(n, S_One)))), coefficient))
  .mi_trig_substitution_rule(integral, rewritten, coefficient, u, cs, substituted, S)
}
.mi_cotcsc_csceven <- function(integral, argument, coefficient, m, n, S) {
  co <- fcall("cot", list(argument)); cs <- fcall("csc", list(argument))
  half <- .div(.sub(n, S_Two), S_Two)
  rewritten <- .mul(list(.pow(.add(list(S_One, .sq(co))), half), .sq(cs), .pow(co, m)))
  u <- Dummy("u")
  substituted <- .neg(.div(.mul(list(.pow(.add(list(S_One, .sq(u))), half), .pow(u, m))), coefficient))
  .mi_trig_substitution_rule(integral, rewritten, coefficient, u, co, substituted, S)
}
# tansec_tan_reduction and cotcsc_cot_reduction
.mi_tan_reduction <- function(integral, argument, coefficient, m, degenerate, S, cot) {
  integrand <- integral[[1L]]; symbol <- integral[[2L]]
  t1 <- fcall(if (cot) "cot" else "tan", list(argument)); s1 <- fcall(if (cot) "csc" else "sec", list(argument))
  reduced <- .pow(t1, .sub(m, S_Two))
  sub_ig <- .mul(list(reduced, .sq(s1)))
  rewritten <- .sub(sub_ig, reduced)
  u <- Dummy("u")
  substituted <- .div(.pow(u, .sub(m, S_Two)), coefficient)
  if (cot) substituted <- .neg(substituted)
  substep <- S$solve(substituted, u)
  sstep <- .mr("URule", sub_ig, symbol, u_var = u, u_func = t1, substep = substep)
  rstep <- if (.int_val(m) == 2L) .mr("ConstantRule", S_One, symbol)
           else (if (cot) trig_cotcsc_rule else trig_tansec_rule)(list(reduced, symbol), S, degenerate = FALSE)
  rstep <- .mr("ConstantTimesRule", .neg(reduced), symbol, constant = S_NegativeOne, other = reduced, substep = rstep)
  generic <- .mr("RewriteRule", integrand, symbol, rewritten = rewritten, substep = .mr("AddRule", rewritten, symbol, substeps = list(sstep, rstep)))
  if (!degenerate) return(generic)
  .mi_add_trig_degenerate_step(integral, coefficient, generic)
}

.mi_single_arg <- function(powers, fa, fb) {
  ks <- unique(c(names(powers[[fa]]), names(powers[[fb]])))
  ks
}

trig_sincos_rule <- function(integral, S) {
  integrand <- integral[[1L]]; symbol <- integral[[2L]]
  col <- .mi_collect_trig_powers(integrand, symbol, c("sin", "cos"))
  if (is.null(col)) return(NULL)
  powers <- col[[1L]]; coefficients <- col[[2L]]
  ks <- .mi_single_arg(powers, "sin", "cos")
  if (length(ks) == 1L) {
    k <- ks[[1L]]
    argument <- coefficients[[k]][[1L]]; coefficient <- coefficients[[k]][[2L]]
    sp <- .mi_pow_get(powers, "sin", k); cp <- .mi_pow_get(powers, "cos", k)
    if (.mi_odd(sp)) return(.mi_sincos_sinodd(integral, argument, coefficient, sp, cp, S))
    if (.mi_odd(cp)) return(.mi_sincos_cosodd(integral, argument, coefficient, sp, cp, S))
    if (.mi_even(sp) && .mi_even(cp) && (.mi_neg(sp) || .mi_neg(cp))) {
      rp <- .neg(.add(list(sp, cp)))
      if (.mi_neg(sp)) { rewritten <- .mul(list(.pow(fcall("cot", list(argument)), cp), .pow(fcall("csc", list(argument)), rp))); rule <- trig_cotcsc_rule }
      else { rewritten <- .mul(list(.pow(fcall("tan", list(argument)), sp), .pow(fcall("sec", list(argument)), rp))); rule <- trig_tansec_rule }
      substep <- rule(list(rewritten, symbol), S)
      if (!is.null(substep)) return(.mr("RewriteRule", integrand, symbol, rewritten = rewritten, substep = substep))
    }
  }
  allp <- unlist(lapply(powers, function(fp) vapply(fp, function(p) p$type == "Integer" && .int_val(p) >= 0, NA)))
  if (all(allp)) return(.mi_sincos_product_to_sum(integral, S))
  NULL
}

trig_tansec_rule <- function(integral, S, degenerate = TRUE) {
  original <- integral[[1L]]; symbol <- integral[[2L]]
  integrand <- .mi_rewrite_to_reciprocal(original, "cos", "sec")
  integrand <- .mi_rewrite_to_reciprocal(integrand, "cot", "tan")
  col <- .mi_collect_trig_powers(integrand, symbol, c("tan", "sec"))
  if (is.null(col)) return(NULL)
  powers <- col[[1L]]; coefficients <- col[[2L]]
  ks <- .mi_single_arg(powers, "tan", "sec")
  if (length(ks) != 1L) return(NULL)
  k <- ks[[1L]]
  argument <- coefficients[[k]][[1L]]; coefficient <- coefficients[[k]][[2L]]
  tp <- .mi_pow_get(powers, "tan", k); sp <- .mi_pow_get(powers, "sec", k)
  if (.mi_neg(tp) && is_zero_num(sp)) return(NULL)
  integral <- list(integrand, symbol)
  if (isTRUE(ask_fact(sp, "integer")) && .mi_neg(sp) && is_zero_num(tp)) {
    rewritten <- .mi_rewrite_to_reciprocal(integrand, "sec", "cos")
    substep <- trig_sincos_rule(list(rewritten, symbol), S)
    if (is.null(substep)) return(NULL)
    step <- .mr("RewriteRule", integrand, symbol, rewritten = rewritten, substep = substep)
  } else if (.mi_odd(tp)) step <- .mi_tansec_tanodd(integral, argument, coefficient, tp, sp, S)
  else if (.mi_even(sp) && !is_zero_num(sp)) step <- .mi_tansec_seceven(integral, argument, coefficient, tp, sp, S)
  else if (is_zero_num(sp) && tp$type == "Integer" && .int_val(tp) >= 2L) step <- .mi_tan_reduction(integral, argument, coefficient, tp, degenerate, S, FALSE)
  else return(NULL)
  if (!.mi_eq(integrand, original)) return(.mr("RewriteRule", original, symbol, rewritten = integrand, substep = step))
  step
}

trig_cotcsc_rule <- function(integral, S, degenerate = TRUE) {
  original <- integral[[1L]]; symbol <- integral[[2L]]
  integrand <- .mi_rewrite_to_reciprocal(original, "sin", "csc")
  integrand <- .mi_rewrite_to_reciprocal(integrand, "tan", "cot")
  col <- .mi_collect_trig_powers(integrand, symbol, c("cot", "csc"))
  if (is.null(col)) return(NULL)
  powers <- col[[1L]]; coefficients <- col[[2L]]
  ks <- .mi_single_arg(powers, "cot", "csc")
  if (length(ks) != 1L) return(NULL)
  k <- ks[[1L]]
  argument <- coefficients[[k]][[1L]]; coefficient <- coefficients[[k]][[2L]]
  tp <- .mi_pow_get(powers, "cot", k); sp <- .mi_pow_get(powers, "csc", k)
  integral <- list(integrand, symbol)
  if (isTRUE(ask_fact(sp, "integer")) && .mi_neg(sp) && is_zero_num(tp)) {
    rewritten <- .mi_rewrite_to_reciprocal(integrand, "csc", "sin")
    substep <- trig_sincos_rule(list(rewritten, symbol), S)
    if (is.null(substep)) return(NULL)
    step <- .mr("RewriteRule", integrand, symbol, rewritten = rewritten, substep = substep)
  } else if (.mi_odd(tp)) step <- .mi_cotcsc_cotodd(integral, argument, coefficient, tp, sp, S)
  else if (.mi_even(sp) && !is_zero_num(sp)) step <- .mi_cotcsc_csceven(integral, argument, coefficient, tp, sp, S)
  else if (is_zero_num(sp) && tp$type == "Integer" && .int_val(tp) >= 2L) step <- .mi_tan_reduction(integral, argument, coefficient, tp, degenerate, S, TRUE)
  else return(NULL)
  if (!.mi_eq(integrand, original)) return(.mr("RewriteRule", original, symbol, rewritten = integrand, substep = step))
  step
}

trig_powers_products_rule <- function(integral, S) {
  for (rule in list(trig_sincos_rule, trig_tansec_rule, trig_cotcsc_rule)) {
    r <- rule(integral, S)
    if (!is.null(r)) return(r)
  }
  NULL
}

trig_cmplx_exp_rule <- function(integral, S) {
  integrand <- integral[[1L]]; symbol <- integral[[2L]]
  factors <- .mi_make_args_mul(integrand)
  if (!any(vapply(factors, function(t) .mi_fn_is(t, "exp") && isTRUE(is_polynomial(t$args[[1L]], symbol)), NA))) return(NULL)
  prim <- c("sin", "cos", "sinh", "cosh")
  is_quad <- function(arg) { p <- P_or_null(arg, symbol); !is.null(p) && P_deg(p) == 2L }
  trig_arg <- function(t) {
    if (.mi_fn_is(t, prim)) return(t$args[[1L]])
    if (t$type == "Pow" && .mi_fn_is(t$args[[1L]], prim)) return(t$args[[1L]]$args[[1L]])
    NULL
  }
  quadratic_phase <- any(vapply(factors, function(t) { a <- trig_arg(t); !is.null(a) && is_quad(a) }, NA))
  a <- Wild("a", exclude = list(symbol, S_Zero)); b <- Wild("b", exclude = list(symbol)); f <- WildFunction("f")
  m <- tryCatch(sym_match(integrand, f * fnode("exp", list(a * symbol + b)) / symbol), error = function(e) NULL)
  trigexp <- !is.null(m) && length(m) > 0L
  gauss <- any(vapply(factors, function(t) .mi_fn_is(t, "exp") && is_quad(t$args[[1L]]), NA))
  if (!(gauss || trigexp || quadratic_phase)) return(NULL)
  rewritten <- .mi_rewrite_only(integrand, prim, "exp")
  if (!.mi_eq(rewritten, integrand)) return(.mr("RewriteRule", integrand, symbol, rewritten = rewritten, substep = S$solve(rewritten, symbol)))
  NULL
}

trig_poly_mul_rule <- function(integral, S) {
  integrand <- integral[[1L]]; symbol <- integral[[2L]]
  if (integrand$type != "Mul") return(NULL)
  a_ <- Wild("a", exclude = list(symbol, S_Zero)); b_ <- Wild("b", exclude = list(symbol, S_Zero)); c_ <- Wild("c", exclude = list(symbol))
  is_sin <- NULL; m <- NULL; rest <- list()
  for (factor in .sympy_args(integrand)) {
    if (is.null(m) && .mi_fn_is(factor, c("sin", "cos"))) {
      mm <- tryCatch(sym_match(factor$args[[1L]], a_ * symbol^2 + b_ * symbol + c_), error = function(e) NULL)
      if (!is.null(mm) && length(mm) && !is.null(mm$a) && !is.null(mm$b)) { m <- mm; is_sin <- factor$name == "sin"; next }
    }
    rest[[length(rest) + 1L]] <- factor
  }
  if (is.null(m)) return(NULL)
  poly <- .mul(rest)
  if (!.has1(poly, symbol) || !isTRUE(is_polynomial(poly, symbol))) return(NULL)
  a <- m$a; b <- m$b; cc <- if (is.null(m$c)) S_Zero else m$c
  cond <- Ne(a, 0)
  dstep <- NULL
  if (!.is_true(cond)) {
    lin <- .add(list(.mul(list(b, symbol)), cc))
    dstep <- S$solve(.mul(list(poly, if (is_sin) sin(lin) else cos(lin))), symbol)
  }
  h <- .neg(.div(b, .mul(list(S_Two, a))))
  k <- .sub(cc, .div(.sq(b), .mul(list(.int(4), a))))
  u <- Dummy("u")
  u_func <- .sub(symbol, h)
  poly_u <- expand(subs(poly, symbol, .add(list(u, h))))
  au2 <- .mul(list(a, .sq(u)))
  et <- if (is_sin) .add(list(.mul(list(sin(au2), cos(k))), .mul(list(cos(au2), sin(k)))))
        else .sub(.mul(list(cos(au2), cos(k))), .mul(list(sin(au2), sin(k))))
  ni <- expand(.mul(list(poly_u, et)))
  rewritten <- subs(ni, u, u_func)
  substep <- S$solve(ni, u)
  if (mr_contains_dont_know(substep)) return(NULL)
  substep <- .mr("URule", rewritten, symbol, u_var = u, u_func = u_func, substep = substep)
  generic <- .mr("CompleteSquareRule", integrand, symbol, rewritten = rewritten, substep = substep)
  .add_degenerate_step(cond, generic, dstep)
}

trig_product_to_sum_rule <- function(integral, S) {
  integrand <- integral[[1L]]; symbol <- integral[[2L]]
  if (integrand$type != "Mul") return(NULL)
  tf <- list(); rest <- list()
  for (f in .sympy_args(integrand)) if (.mi_fn_is(f, c("sin", "cos"))) tf[[length(tf) + 1L]] <- f else rest[[length(rest) + 1L]] <- f
  if (length(tf) != 2L) return(NULL)
  A <- tf[[1L]]$args[[1L]]; B <- tf[[2L]]$args[[1L]]
  if (!(.has1(A, symbol) && .has1(B, symbol))) return(NULL)
  lin <- function(e) { d <- .mi_poly_deg(e, symbol); !is.na(d) && d <= 1L }
  if (lin(A) && lin(B)) return(NULL)
  rewritten <- expand(.mul(list(.mul(rest), sincos_to_sum(.mul(tf)))))
  substep <- S$solve(rewritten, symbol)
  if (mr_contains_dont_know(substep)) return(NULL)
  .mr("RewriteRule", integrand, symbol, rewritten = rewritten, substep = substep)
}

hyperbolic_rule <- function(integral, S) {
  integrand <- integral[[1L]]; symbol <- integral[[2L]]
  if (!(.mi_fn_is(integrand, .HYP_FUNS) && .mi_eq(integrand$args[[1L]], symbol))) return(NULL)
  nm <- integrand$name
  if (nm == "sinh") return(.mr("SinhRule", integrand, symbol))
  if (nm == "cosh") return(.mr("CoshRule", integrand, symbol))
  u <- Dummy("u")
  if (nm == "tanh") {
    rw <- .div(sinh(symbol), cosh(symbol))
    return(.mr("RewriteRule", integrand, symbol, rewritten = rw, substep = .mr("URule", rw, symbol, u_var = u, u_func = cosh(symbol), substep = .mr("ReciprocalRule", .inv(u), u, base = u))))
  }
  if (nm == "coth") {
    rw <- .div(cosh(symbol), sinh(symbol))
    return(.mr("RewriteRule", integrand, symbol, rewritten = rw, substep = .mr("URule", rw, symbol, u_var = u, u_func = sinh(symbol), substep = .mr("ReciprocalRule", .inv(u), u, base = u))))
  }
  rw <- rewrite(integrand, "tanh")
  th <- fcall("tanh", list(.div(symbol, S_Two)))
  if (nm == "sech") return(.mr("RewriteRule", integrand, symbol, rewritten = rw, substep = .mr("URule", rw, symbol, u_var = u, u_func = th,
                         substep = .mr("ArctanRule", .div(S_Two, .add(list(.sq(u), S_One))), u, a = S_Two, b = S_One, c = S_One))))
  if (nm == "csch") return(.mr("RewriteRule", integrand, symbol, rewritten = rw, substep = .mr("URule", rw, symbol, u_var = u, u_func = th,
                         substep = .mr("ReciprocalRule", .inv(u), u, base = u))))
  NULL
}

heaviside_rule <- function(integral, S) {
  integrand <- integral[[1L]]; symbol <- integral[[2L]]
  m <- Wild("m", exclude = list(symbol)); b <- Wild("b", exclude = list(symbol)); g <- Wild("g")
  mm <- tryCatch(sym_match(integrand, fnode("Heaviside", list(m * symbol + b)) * g), error = function(e) NULL)
  if (is.null(mm) || !length(mm) || is.null(mm$m)) return(NULL)
  gv <- if (is.null(mm$g)) S_One else mm$g
  if (is_zero_num(gv)) return(NULL)
  bv <- if (is.null(mm$b)) S_Zero else mm$b
  substep <- S$solve(gv, symbol)
  .mr("HeavisideRule", integrand, symbol, harg = .add(list(.mul(list(mm$m, symbol)), bv)), ibnd = .neg(.div(bv, mm$m)), substep = substep)
}

dirac_delta_rule <- function(integral, S) {
  integrand <- integral[[1L]]; x <- integral[[2L]]
  n <- if (length(integrand$args) == 1L) S_Zero else integrand$args[[2L]]
  if (n$type != "Integer" || .int_val(n) < 0) return(NULL)
  a <- Wild("a", exclude = list(x)); b <- Wild("b", exclude = list(x, S_Zero))
  m <- tryCatch(sym_match(integrand$args[[1L]], a + b * x), error = function(e) NULL)
  if (is.null(m) || !length(m) || is.null(m$b)) return(NULL)
  av <- if (is.null(m$a)) S_Zero else m$a; bv <- m$b
  cond <- Ne(bv, 0)
  dstep <- if (.is_true(cond)) NULL else .mr("ConstantRule", fcall("DiracDelta", list(av, n)), x)
  .add_degenerate_step(cond, .mr("DiracDeltaRule", integrand, x, n = n, a = av, b = bv), dstep)
}

derivative_rule <- function(integral, S) {
  integrand <- integral[[1L]]; symbol <- integral[[2L]]
  dvars <- lapply(integrand$args[-1L], function(t) if (t$type == "Tuple") t$args[[1L]] else t)
  if (.has1(integrand$args[[1L]], symbol)) {
    if (any(vapply(dvars, function(v) .mi_eq(v, symbol), NA))) return(.mr("DerivativeRule", integrand, symbol))
    return(.mr("DontKnowRule", integrand, symbol))
  }
  .mr("ConstantRule", integrand, symbol)
}

rewrites_rule <- function(integral, S) {
  integrand <- integral[[1L]]; symbol <- integral[[2L]]
  # Expr.match without wilds returns an empty (falsy) dict on success, so
  # SymPy never takes this branch
  m <- tryCatch(sym_match(integrand, 1 / cos(symbol)), error = function(e) NULL)
  if (!is.null(m) && length(m)) {
    rewritten <- subs(integrand, .inv(cos(symbol)), fcall("sec", list(symbol)))
    return(.mr("RewriteRule", integrand, symbol, rewritten = rewritten, substep = S$solve(rewritten, symbol)))
  }
  NULL
}

partial_fractions_rule <- function(integral, S) {
  integrand <- integral[[1L]]; symbol <- integral[[2L]]
  if (!.mi_is_ratfun(integrand, symbol)) return(NULL)
  rewritten <- apart(integrand, symbol)
  if (.mi_eq(rewritten, integrand)) return(.mr("RatintRule", integrand, symbol))
  substep <- S$solve(rewritten, symbol)
  if (!.is_mr(substep, "DontKnowRule")) return(.mr("RewriteRule", integrand, symbol, rewritten = rewritten, substep = substep))
  NULL
}
cancel_rule <- .mi_rewriter(function(ig, s) !.mi_is_ratfun(ig, s), function(ig, s) cancel(ig))
distribute_expand_rule <- .mi_rewriter(
  function(ig, s) (ig$type %in% c("Pow", "Mul") || all(vapply(.mi_args(ig), function(a) a$type == "Pow" || .mi_is_poly(a, s), NA))) && !.mi_is_ratfun(ig, s),
  function(ig, s) expand(ig))
trig_expand_rule <- .mi_rewriter(
  function(ig, s) length(unique(vapply(.mi_find_fns(ig, .TRIG_FUNS), function(a) a$args[[1L]]$key, ""))) > 1L,
  function(ig, s) expand(ig, trig = TRUE))

# ---------------------------------------------------------------------------
# substitution

substitution_rule <- function(integral, S) {
  integrand <- integral[[1L]]; symbol <- integral[[2L]]
  u_var <- Dummy("u")
  subst <- find_substitutions(integrand, symbol, u_var)
  if (!length(subst)) return(NULL)
  branch <- S$branch
  ways <- list()
  factored <- sym_factor(integrand)
  denom_integrand <- .as_numer_denom(factored)[[2L]]
  for (sb in subst) {
    u_func <- sb[[1L]]; cc <- sb[[2L]]; substituted <- sb[[3L]]
    subrule <- S$solve(substituted, u_var)
    if (mr_contains_dont_know(subrule)) next
    if (!is_zero_num(simplify(.sub(cc, S_One)))) {
      denom_c <- .as_numer_denom(cc)[[2L]]
      subrule <- .mr("ConstantTimesRule", .mul(list(cc, substituted)), u_var, constant = cc, other = substituted, substep = subrule)
      if (length(free_symbols(denom_c))) {
        pieces <- list()
        for (fk in factor_list(denom_c)$factors) {
          pole <- fk[[1L]]
          if (!.if_zero_implies_zero(pole, denom_integrand)) {
            rewritten <- manual_subs(factored, pole, S_Zero)
            if (.mi_has_infinite(rewritten)) next
            substep <- S$solve(rewritten, symbol)
            substep <- .mr("RewriteRule", integrand, symbol, rewritten = rewritten, substep = substep)
            pieces[[length(pieces) + 1L]] <- list(substep, Eq(pole, 0))
          }
        }
        if (length(pieces)) {
          pieces[[length(pieces) + 1L]] <- list(subrule, S_true)
          subrule <- .mr("PiecewiseRule", substituted, symbol, subfunctions = pieces)
        }
      }
    }
    rule <- .mr("URule", integrand, symbol, u_var = u_var, u_func = u_func, substep = subrule)
    if (!branch) return(rule)
    ways[[length(ways) + 1L]] <- rule
  }
  if (length(ways) > 1L) return(.mr("AlternativeRule", integrand, symbol, alternatives = ways))
  if (length(ways)) return(ways[[1L]])
  NULL
}

# ---------------------------------------------------------------------------
# inverse trig and quadratic denominators

inverse_trig_rule <- function(integral, S, degenerate = TRUE) {
  integrand <- integral[[1L]]; symbol <- integral[[2L]]
  be <- .as_base_exp(integrand); base <- be[[1L]]; ex <- be[[2L]]
  a <- Wild("a", exclude = list(symbol)); b <- Wild("b", exclude = list(symbol)); c <- Wild("c", exclude = list(symbol, S_Zero))
  m <- tryCatch(sym_match(base, a + b * symbol + c * symbol^2), error = function(e) NULL)
  if (is.null(m) || !length(m) || is.null(m$c)) return(NULL)
  make_inverse_trig <- function(rtype, a, sign_a, c, sign_c, h) {
    u_var <- Dummy("u")
    rewritten <- .inv(.sqrt(.add(list(.mul(list(.int(sign_a), a)), .mul(list(.int(sign_c), c, .sq(.sub(symbol, h))))))))
    qb <- .mul(list(.sqrt(.div(c, a)), .sub(symbol, h)))
    constant <- .inv(.sqrt(c))
    u_func <- NULL
    if (!.mi_eq(qb, symbol)) { u_func <- qb; qb <- u_var }
    std <- .inv(.sqrt(.add(list(.int(sign_a), .mul(list(.int(sign_c), .sq(qb)))))))
    substep <- .mr(rtype, std, qb)
    if (!is_one_num(constant)) substep <- .mr("ConstantTimesRule", .mul(list(constant, std)), symbol, constant = constant, other = std, substep = substep)
    if (!is.null(u_func)) substep <- .mr("URule", rewritten, symbol, u_var = u_var, u_func = u_func, substep = substep)
    if (!is_zero_num(h)) substep <- .mr("CompleteSquareRule", integrand, symbol, rewritten = rewritten, substep = substep)
    substep
  }
  av <- if (is.null(m$a)) S_Zero else m$a; bv <- if (is.null(m$b)) S_Zero else m$b; cv <- m$c
  cond <- Ne(cv, 0)
  if (!degenerate || .is_true(cond)) dstep <- NULL
  else if (is_zero_num(bv)) dstep <- .mr("ConstantRule", .pow(av, ex), symbol)
  else dstep <- sqrt_fractional_linear_rule(list(.pow(.add(list(av, .mul(list(bv, symbol)))), ex), symbol), S)
  if (is_zero_num(simplify(.add(list(.mul(list(S_Two, ex)), S_One))))) {
    h <- .neg(.div(bv, .mul(list(S_Two, cv)))); k <- .sub(av, .div(.sq(bv), .mul(list(.int(4), cv))))
    nsc <- Ne(k, 0)
    square_step <- NULL
    if (!.is_true(nsc)) square_step <- .mr("NestedPowRule", .inv(.sqrt(.mul(list(cv, .sq(.sub(symbol, h)))))), symbol, base = .sub(symbol, h), exp = S_NegativeOne)
    if (identical(nsc$key, S_false$key)) return(square_step)
    generic <- .mr("ReciprocalSqrtQuadraticRule", integrand, symbol, a = av, b = bv, c = cv)
    step <- .add_degenerate_step(nsc, generic, square_step)
    if (isTRUE(ask_fact(k, "real")) && isTRUE(ask_fact(cv, "real"))) {
      rules <- list()
      cands <- list(list(list("ArcsinRule", k, 1L, .neg(cv), -1L, h), And(Gt(k, 0), Lt(cv, 0))),
                    list(list("ArcsinhRule", k, 1L, cv, 1L, h), And(Gt(k, 0), Gt(cv, 0))))
      for (cd in cands) {
        if (.is_true(cd[[2L]])) return(do.call(make_inverse_trig, cd[[1L]]))
        if (!identical(cd[[2L]]$key, S_false$key)) rules[[length(rules) + 1L]] <- list(do.call(make_inverse_trig, cd[[1L]]), cd[[2L]])
      }
      if (length(rules)) {
        if (!isTRUE(ask_fact(k, "positive"))) rules[[length(rules) + 1L]] <- list(generic, S_true)
        step <- .mr("PiecewiseRule", integrand, symbol, subfunctions = rules)
      } else step <- generic
    }
    return(.add_degenerate_step(cond, step, dstep))
  }
  if (.mi_eq(ex, S_Half)) {
    substep <- sqrt_quadratic_rule(list(integrand, symbol), S, degenerate = FALSE)
    step <- .mr("SqrtQuadraticRule", integrand, symbol, a = av, b = bv, c = cv, substep = substep)
    return(.add_degenerate_step(cond, step, dstep))
  }
  NULL
}

quadratic_denom_rule <- function(integral, S) {
  integrand <- integral[[1L]]; symbol <- integral[[2L]]
  if (!.mi_is_ratfun(integrand, symbol)) return(NULL)
  nd <- .as_numer_denom(integrand); num <- nd[[1L]]; den <- nd[[2L]]
  if (is_one_num(den)) return(NULL)
  ci <- .as_coeff_mul_dep(den, symbol)
  if (!is_one_num(ci[[1L]])) num <- .div(num, ci[[1L]])
  den <- .mul(ci[[2L]])
  if (den$type == "Pow") { q <- den$args[[1L]]; n <- den$args[[2L]] } else { q <- den; n <- S_One }
  den_poly <- .P(q, symbol); num_poly <- .P(num, symbol)
  if (P_deg(den_poly) != 2L) return(NULL)
  if (!P_is_zero(num_poly) && P_deg(num_poly) >= 2L) return(NULL)
  deg_num <- if (P_is_zero(num_poly)) -Inf else P_deg(num_poly)

  arctan_match <- function(B, a, c, sym, degenerate = TRUE) {
    ig <- .div(B, .add(list(.mul(list(a, .sq(sym))), c)))
    pieces <- list()
    if (degenerate && !.if_zero_implies_zero(a, c)) {
      st <- subs(ig, a, S_Zero)
      pieces[[length(pieces) + 1L]] <- list(.mr("RewriteRule", ig, sym, rewritten = st, substep = S$solve(st, sym)), Eq(a, 0))
    }
    if (degenerate && !.if_zero_implies_zero(c, a)) {
      st <- subs(ig, c, S_Zero)
      pieces[[length(pieces) + 1L]] <- list(.mr("RewriteRule", ig, sym, rewritten = st, substep = S$solve(st, sym)), Eq(c, 0))
    }
    if (isTRUE(ask_fact(a, "extended_real")) && isTRUE(ask_fact(c, "extended_real"))) {
      pc <- Gt(.div(c, a), 0)
      if (!.is_true(pc)) {
        coeff <- .div(B, .mul(list(S_Two, .sqrt(.neg(c)), .sqrt(a))))
        constant <- .sqrt(.neg(.div(c, a)))
        r1 <- .inv(.sub(sym, constant)); r2 <- .inv(.add(list(sym, constant)))
        log_steps <- list(.mr("ReciprocalRule", r1, sym, base = .sub(sym, constant)),
                          .mr("ConstantTimesRule", .neg(r2), sym, constant = S_NegativeOne, other = r2, substep = .mr("ReciprocalRule", r2, sym, base = .add(list(sym, constant)))))
        sub <- .sub(r1, r2); rewritten <- sub
        neg_step <- .mr("AddRule", sub, sym, substeps = log_steps)
        if (!is_one_num(coeff)) {
          rewritten <- .mul(list(coeff, sub))
          neg_step <- .mr("ConstantTimesRule", rewritten, sym, constant = coeff, other = sub, substep = neg_step)
        }
        neg_step <- .mr("RewriteRule", ig, sym, rewritten = rewritten, substep = neg_step)
        if (identical(pc$key, S_false$key)) {
          pieces[[length(pieces) + 1L]] <- list(neg_step, S_true)
          return(.mr("PiecewiseRule", ig, sym, subfunctions = pieces))
        }
        pieces[[length(pieces) + 1L]] <- list(neg_step, Lt(.div(c, a), 0))
      }
    }
    general <- .mr("ArctanRule", ig, sym, a = B, b = a, c = c)
    if (length(pieces)) {
      pieces[[length(pieces) + 1L]] <- list(general, S_true)
      return(.mr("PiecewiseRule", ig, sym, subfunctions = pieces))
    }
    general
  }

  complete_square <- function(B, a, b, c, n, sym, degenerate_a = TRUE, degenerate_discriminant = TRUE) {
    pieces <- list()
    disc <- expand(.sub(.mul(list(.int(4), a, c)), .sq(b)))
    denominator <- .add(list(.mul(list(a, .sq(sym))), .mul(list(b, sym)), c))
    ig <- .div(B, .pow(denominator, n))
    if (degenerate_a && !.if_zero_implies_zero(a, denominator)) {
      st <- if (isTRUE(ask_fact(disc, "zero"))) subs(ig, list(list(a, S_Zero), list(b, S_Zero)), simultaneous = TRUE) else subs(ig, a, S_Zero)
      pieces[[length(pieces) + 1L]] <- list(.mr("RewriteRule", ig, sym, rewritten = st, substep = S$solve(st, sym)), Eq(a, 0))
    }
    if (degenerate_discriminant && !.if_zero_implies_zero(disc, denominator)) {
      u <- Dummy("u")
      u_func <- .add(list(sym, .div(b, .mul(list(S_Two, a)))))
      rewritten <- .mul(list(.div(B, .pow(a, n)), .pow(u_func, .mul(list(.int(-2), n)))))
      subexpr <- .mul(list(.div(B, .pow(a, n)), .pow(u, .mul(list(.int(-2), n)))))
      rule <- .mr("RewriteRule", ig, sym, rewritten = rewritten, substep = .mr("URule", rewritten, sym, u_var = u, u_func = u_func, substep = S$solve(subexpr, u)))
      if (isTRUE(ask_fact(disc, "zero"))) {
        if (length(pieces)) {
          pieces[[length(pieces) + 1L]] <- list(rule, S_true)
          return(.mr("PiecewiseRule", ig, sym, subfunctions = pieces))
        }
        return(rule)
      }
      pieces[[length(pieces) + 1L]] <- list(rule, Eq(disc, 0))
    }
    if (is_one_num(n)) {
      u <- Dummy("u")
      u_func <- .add(list(sym, .div(b, .mul(list(S_Two, a)))))
      substep <- arctan_match(B, a, .div(disc, .mul(list(.int(4), a))), u, degenerate = FALSE)
      general <- .mr("URule", ig, sym, u_var = u, u_func = u_func, substep = substep)
    } else {
      n1 <- .sub(n, S_One)
      T <- .add(list(.mul(list(S_Two, a, sym)), b))
      F <- .div(.mul(list(B, T)), .mul(list(n1, disc, .pow(denominator, n1))))
      derivative <- new_sym("Derivative", list(F, Tuple(sym, S_One)))
      coeff <- .div(.mul(list(S_Two, a, .sub(.mul(list(S_Two, n)), .int(3)))), .mul(list(n1, disc)))
      remainder <- .div(B, .pow(denominator, n1))
      scaled <- .mul(list(coeff, remainder))
      rewritten <- .add(list(derivative, scaled))
      dstep <- .mr("DerivativeRule", derivative, sym)
      rstep <- complete_square(B, a, b, c, n1, sym, degenerate_a = FALSE, degenerate_discriminant = FALSE)
      sstep <- .mr("ConstantTimesRule", scaled, sym, constant = coeff, other = remainder, substep = rstep)
      general <- .mr("RewriteRule", ig, sym, rewritten = rewritten, substep = .mr("AddRule", rewritten, sym, substeps = list(dstep, sstep)))
    }
    if (length(pieces)) {
      pieces[[length(pieces) + 1L]] <- list(general, S_true)
      return(.mr("PiecewiseRule", ig, sym, subfunctions = pieces))
    }
    general
  }

  split_sum <- function(A, B, a, b, c, n, sym) {
    pieces <- list()
    denominator <- .add(list(.mul(list(a, .sq(sym))), .mul(list(b, sym)), c))
    ig <- .div(.add(list(.mul(list(A, sym)), B)), .pow(denominator, n))
    disc <- expand(.sub(.mul(list(.int(4), a, c)), .sq(b)))
    if (!.if_zero_implies_zero(a, denominator)) {
      st <- if (isTRUE(ask_fact(disc, "zero"))) subs(ig, list(list(a, S_Zero), list(b, S_Zero)), simultaneous = TRUE) else subs(ig, a, S_Zero)
      pieces[[length(pieces) + 1L]] <- list(.mr("RewriteRule", ig, sym, rewritten = st, substep = S$solve(st, sym)), Eq(a, 0))
    }
    const <- .div(A, .mul(list(S_Two, a)))
    numer1 <- .add(list(.mul(list(S_Two, a, sym)), b))
    numer2 <- .add(list(.neg(.mul(list(const, b))), B))
    qprime <- .div(numer1, .pow(denominator, n))
    u <- Dummy("u")
    step1 <- .mr("URule", qprime, sym, u_var = u, u_func = denominator, substep = S$solve(.pow(u, .neg(n)), u))
    if (!is_one_num(const)) step1 <- .mr("ConstantTimesRule", .mul(list(const, qprime)), sym, constant = const, other = qprime, substep = step1)
    if (isTRUE(ask_fact(numer2, "zero"))) {
      rewritten <- .mul(list(const, qprime))
      general <- .mr("RewriteRule", ig, sym, rewritten = rewritten, substep = step1)
    } else {
      step2 <- complete_square(numer2, a, b, c, n, sym, degenerate_a = FALSE)
      rewritten <- .add(list(.mul(list(const, qprime)), .div(numer2, .pow(denominator, n))))
      general <- .mr("RewriteRule", ig, sym, rewritten = rewritten, substep = .mr("AddRule", rewritten, sym, substeps = list(step1, step2)))
    }
    if (length(pieces)) {
      pieces[[length(pieces) + 1L]] <- list(general, S_true)
      return(.mr("PiecewiseRule", ig, sym, subfunctions = pieces))
    }
    general
  }

  B <- P_nth(num_poly, 0L)
  a <- P_nth(den_poly, 2L); b <- P_nth(den_poly, 1L); c <- P_nth(den_poly, 0L)
  normalized <- .div(P_expr(num_poly), .pow(P_expr(den_poly), n))
  if (is_zero_num(b) && deg_num == 0 && is_one_num(n)) step <- arctan_match(B, a, c, symbol)
  else if (deg_num == 1) step <- split_sum(P_nth(num_poly, 1L), B, a, b, c, n, symbol)
  else step <- complete_square(B, a, b, c, n, symbol)
  if (!.mi_eq(normalized, integrand)) step <- .mr("RewriteRule", integrand, symbol, rewritten = normalized, substep = step)
  step
}
