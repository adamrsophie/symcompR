# manualintegrate, after sympy/integrals/manualintegrate.py (part 3): the
# algebraic substitutions (fractional linear roots, Euler, square roots of
# quadratics, Chebyshev binomials), Bioche's rules, the IntegrationSolver,
# integral_steps and manualintegrate.

# sympy.simplify.fu.sincos_to_sum
sincos_to_sum <- function(expr) {
  expr <- as_sym(expr)
  if (!.has_function(expr, c("sin", "cos"))) return(expr)
  TR8(expand_mul(TRpower(expr)))
}

.mi_pow_ratio <- function(pw, base0, ex) cancel(powsimp(.div(pw, .pow(base0, ex)), force = FALSE))

sqrt_fractional_linear_rule <- function(integral, S) {
  integrand <- integral[[1L]]; x <- integral[[2L]]
  a <- Wild("a", exclude = list(x)); b <- Wild("b", exclude = list(x)); c <- Wild("c", exclude = list(x)); d <- Wild("d", exclude = list(x))
  base0 <- NULL
  powers <- list(); exps <- list(); ratios <- list()
  cbs <- list()
  for (pw in .mi_sort_ordered(.mi_find_pows(integrand))) {
    base <- pw$args[[1L]]; ex <- pw$args[[2L]]
    if (ex$type == "Integer" || !.has1(base, x)) next
    if (!.mi_is_rational_num(ex)) return(NULL)
    nd <- .as_numer_denom(base)
    mn <- tryCatch(sym_match(nd[[1L]], a * x + b), error = function(e) NULL)
    md <- tryCatch(sym_match(nd[[2L]], c * x + d), error = function(e) NULL)
    if (is.null(mn) || is.null(md) || !length(mn) || !length(md)) next
    g <- function(m, k) if (is.null(m[[k]])) S_Zero else m[[k]]
    aa <- g(mn, "a"); bb <- g(mn, "b"); cc <- g(md, "c"); dd <- g(md, "d")
    if (isTRUE(ask_fact(cc, "zero")) && isTRUE(ask_fact(dd, "zero"))) return(NULL)
    det <- .sub(.mul(list(aa, dd)), .mul(list(bb, cc)))
    if (isTRUE(ask_fact(det, "zero"))) {
      cbs[[length(cbs) + 1L]] <- list(base, if (!isTRUE(ask_fact(cc, "zero"))) .div(aa, cc) else .div(bb, dd))
      next
    }
    if (is.null(base0)) {
      base0 <- base; a0 <- aa; b0 <- bb; c0 <- cc; d0 <- dd
      powers[[1L]] <- pw; exps[[1L]] <- ex; ratios[[1L]] <- S_One
    } else {
      pr <- .mi_pow_ratio(pw, base0, ex)
      if (.has1(pr, x)) return(NULL)
      powers[[length(powers) + 1L]] <- pw; exps[[length(exps) + 1L]] <- ex; ratios[[length(ratios) + 1L]] <- pr
    }
  }
  if (is.null(base0) && !length(cbs)) return(NULL)
  if (length(cbs)) integrand <- xreplace(integrand, cbs)
  if (is.null(base0)) {
    substep <- S$solve(integrand, x)
    if (!mr_contains_dont_know(substep)) return(.mr("RewriteRule", integral[[1L]], x, rewritten = integrand, substep = substep))
    return(NULL)
  }
  q0 <- Reduce(function(p, q) { g <- function(a, b) if (b == 0) a else g(b, a %% b); p * q / g(p, q) }, vapply(exps, .mi_q, 0L))
  u <- Dummy("u")
  u_x <- .pow(base0, Rational(1, q0))
  u_pow <- .pow(u, .int(q0))
  x_u <- .div(.sub(b0, .mul(list(d0, u_pow))), .sub(.mul(list(c0, u_pow)), a0))
  dx_u <- .div(.mul(list(.int(q0), .sub(.mul(list(a0, d0)), .mul(list(b0, c0))), .pow(u, .int(q0 - 1L)))), .sq(.sub(.mul(list(c0, u_pow)), a0)))
  sd <- lapply(seq_along(powers), function(i) list(powers[[i]], .mul(list(ratios[[i]], .pow(u, .mul(list(.int(q0), exps[[i]])))))))
  substituted <- .mul(list(xreplace(xreplace(integrand, sd), list(list(x, x_u))), dx_u))
  substep <- S$solve(substituted, u)
  if (mr_contains_dont_know(substep)) return(NULL)
  pieces <- list()
  det <- .sub(.mul(list(a0, d0)), .mul(list(b0, c0)))
  b0den <- .as_numer_denom(base0)[[2L]]
  if (!.if_zero_implies_zero(det, b0den)) {
    d0c0 <- .if_zero_implies_zero(d0, c0); c0d0 <- .if_zero_implies_zero(c0, d0)
    if (!d0c0 || (c0d0 && d0c0)) {
      cv <- .div(a0, c0)
      sa <- lapply(seq_along(powers), function(i) list(powers[[i]], .mul(list(ratios[[i]], .pow(cv, exps[[i]])))))
      pieces[[length(pieces) + 1L]] <- list(S$solve(xreplace(integrand, sa), x), And(Eq(det, 0), Ne(c0, 0)))
    }
    if (!c0d0) {
      cv <- .div(b0, d0)
      sb <- lapply(seq_along(powers), function(i) list(powers[[i]], .mul(list(ratios[[i]], .pow(cv, exps[[i]])))))
      simp <- subs(xreplace(integrand, sb), list(list(a0, S_Zero), list(c0, S_Zero)))
      pieces[[length(pieces) + 1L]] <- list(S$solve(simp, x), And(Eq(det, 0), Eq(c0, 0)))
    }
  }
  step <- .mr("URule", integrand, x, u_var = u, u_func = u_x, substep = substep)
  if (length(pieces)) {
    pieces[[length(pieces) + 1L]] <- list(step, S_true)
    step <- .mr("PiecewiseRule", integrand, x, subfunctions = pieces)
  }
  if (length(cbs)) return(.mr("RewriteRule", integral[[1L]], x, rewritten = integrand, substep = step))
  step
}

euler_substitution_rule <- function(integral, S) {
  integrand <- integral[[1L]]; x <- integral[[2L]]
  base0 <- NULL
  powers <- list(); exps <- list(); ratios <- list()
  for (pw in .mi_sort_ordered(.mi_find_pows(integrand))) {
    base <- pw$args[[1L]]; ex <- pw$args[[2L]]
    if (ex$type == "Integer" || !.has1(base, x)) next
    if (!.mi_is_rational_num(ex)) return(NULL)
    if (.mi_q(ex) != 2L) return(NULL)
    bp <- P_or_null(base, x)
    if (is.null(bp) || P_deg(bp) != 2L) return(NULL)
    R <- P_expr(bp)
    if (is.null(base0)) {
      base0 <- R; a0 <- P_nth(bp, 0L); b0 <- P_nth(bp, 1L); c0 <- P_nth(bp, 2L)
      powers[[1L]] <- pw; exps[[1L]] <- ex; ratios[[1L]] <- S_One
    } else {
      pr <- cancel(powsimp(.div(.pow(R, ex), .pow(base0, ex)), force = FALSE))
      if (.has1(pr, x)) return(NULL)
      powers[[length(powers) + 1L]] <- pw; exps[[length(exps) + 1L]] <- ex; ratios[[length(ratios) + 1L]] <- pr
    }
  }
  if (is.null(base0)) return(NULL)
  pieces <- list()
  delta <- .sub(.mul(list(.int(4), a0, c0)), .sq(b0))
  delta_zero_step <- function() {
    shift <- .add(list(x, .div(b0, .mul(list(S_Two, c0)))))
    rb <- .mul(list(c0, .sq(shift)))
    sd <- lapply(seq_along(powers), function(i) list(powers[[i]], .mul(list(ratios[[i]], .pow(rb, exps[[i]])))))
    rewritten <- xreplace(integrand, sd)
    .mr("RewriteRule", integrand, x, rewritten = rewritten, substep = S$solve(rewritten, x))
  }
  c_zero_step <- function() {
    di <- subs(integrand, c0, S_Zero)
    if (isTRUE(ask_fact(b0, "zero"))) return(S$solve(di, x))
    st <- sqrt_fractional_linear_rule(list(di, x), S)
    if (is.null(st)) st <- .mr("DontKnowRule", di, x)
    st
  }
  general_euler_step <- function() {
    s <- Dummy("s")
    sd <- lapply(seq_along(powers), function(i) list(powers[[i]], .mul(list(ratios[[i]], .pow(s, .mul(list(S_Two, exps[[i]])))))))
    rewritten <- xreplace(integrand, sd)
    nd <- .as_numer_denom(rewritten)
    if (!.mi_is_poly(nd[[1L]], list(x, s)) || !.mi_is_poly(nd[[2L]], list(x, s))) return(NULL)
    u <- Dummy("u")
    sc <- .sqrt(c0)
    den <- .add(list(b0, .mul(list(S_Two, sc, u))))
    x_u <- .div(.sub(.sq(u), a0), den)
    s_u <- .sub(u, .mul(list(sc, x_u)))
    dx_u <- .div(.mul(list(S_Two, .add(list(.mul(list(b0, u)), .mul(list(sc, .add(list(.sq(u), a0)))))))), .sq(den))
    substituted <- .mul(list(xreplace(rewritten, list(list(x, x_u), list(s, s_u))), dx_u))
    substep <- S$solve(substituted, u)
    .mr("URule", integrand, x, u_var = u, u_func = .add(list(.sqrt(base0), .mul(list(sc, x)))), substep = substep)
  }
  dz <- Eq(delta, 0)
  if (.is_true(dz)) {
    general <- delta_zero_step()
    if (mr_contains_dont_know(general)) return(NULL)
  } else {
    general <- general_euler_step()
    if (is.null(general) || mr_contains_dont_know(general)) return(NULL)
  }
  if (is.na(ask_fact(c0, "zero"))) pieces[[length(pieces) + 1L]] <- list(c_zero_step(), Eq(c0, 0))
  if (is.na(ask_fact(delta, "zero"))) pieces[[length(pieces) + 1L]] <- list(delta_zero_step(), dz)
  if (length(pieces)) {
    pieces[[length(pieces) + 1L]] <- list(general, S_true)
    general <- .mr("PiecewiseRule", integrand, x, subfunctions = pieces)
  }
  general
}

sqrt_quadratic_rule <- function(integral, S, degenerate = TRUE) {
  starting <- integral[[1L]]; x <- integral[[2L]]
  f <- S_One; root_base <- NULL; root_exp <- S_Zero
  for (factor in .mi_make_args_mul(starting)) {
    if (!.has1(factor, x)) { f <- .mul(list(f, factor)); next }
    be <- .as_base_exp(factor); base <- be[[1L]]; ex <- be[[2L]]
    if (ex$type == "Integer") { f <- .mul(list(f, factor)); next }
    if (!.mi_is_rational_num(ex)) return(NULL)
    bp <- P_or_null(base, x)
    if (is.null(bp) || P_deg(bp) != 2L) return(NULL)
    base <- P_expr(bp)
    if (is.null(root_base)) { root_base <- base; root_exp <- ex; next }
    ratio <- cancel(powsimp(.div(factor, .pow(root_base, ex)), force = FALSE))
    if (.has1(ratio, x)) return(NULL)
    f <- .mul(list(f, ratio)); root_exp <- .add(list(root_exp, ex))
  }
  if (is.null(root_base)) return(NULL)
  f_poly <- P_or_null(f, x)
  if (is.null(f_poly)) return(NULL)
  n <- .mul(list(S_Two, root_exp))
  if (n$type != "Integer" || .int_val(n) %% 2L == 0L) return(NULL)
  nv <- .int_val(n)
  rp <- .P(root_base, x)
  a <- P_nth(rp, 0L); b <- P_nth(rp, 1L); c <- P_nth(rp, 2L)
  root_base <- .add(list(a, .mul(list(b, x)), .mul(list(c, .sq(x)))))
  integrand <- .mul(list(f, .pow(root_base, Rational(nv, 2L))))
  cond <- Ne(c, 0)
  if (!degenerate || .is_true(cond)) dstep <- NULL
  else if (isTRUE(ask_fact(b, "zero"))) dstep <- S$solve(.mul(list(f, .pow(.sqrt(a), n))), x)
  else {
    di <- .mul(list(f, .pow(.sqrt(.add(list(a, .mul(list(b, x))))), n)))
    dstep <- sqrt_fractional_linear_rule(list(di, x), S)
    if (is.null(dstep)) dstep <- .mr("DontKnowRule", di, x)
  }
  denom_rule <- function(numer_poly, ig) {
    denom <- .sqrt(.add(list(a, .mul(list(b, x)), .mul(list(c, .sq(x))))))
    deg <- P_deg(numer_poly)
    if (deg <= 1L) {
      if (deg == 1L) { e <- P_nth(numer_poly, 1L); d <- P_nth(numer_poly, 0L) } else { e <- S_Zero; d <- P_expr(numer_poly) }
      A <- .div(e, .mul(list(S_Two, c))); B <- .sub(d, .mul(list(A, b)))
      pre <- .div(.add(list(.mul(list(S_Two, c, x)), b)), denom)
      cstep <- NULL; lstep <- NULL
      if (!is_zero_num(A)) {
        u <- Dummy("u")
        pr <- .mr("PowerRule", .inv(.sqrt(u)), u, base = u, exp = Rational(-1, 2))
        lstep <- .mr("URule", pre, x, u_var = u, u_func = .add(list(a, .mul(list(b, x)), .mul(list(c, .sq(x))))), substep = pr)
        if (!is_one_num(A)) lstep <- .mr("ConstantTimesRule", .mul(list(A, pre)), x, constant = A, other = pre, substep = lstep)
      }
      if (!is_zero_num(B)) {
        cstep <- inverse_trig_rule(list(.inv(denom), x), S, degenerate = FALSE)
        if (!is_one_num(B)) cstep <- .mr("ConstantTimesRule", .div(B, denom), x, constant = B, other = .inv(denom), substep = cstep)
      }
      if (!is.null(lstep) && !is.null(cstep)) {
        add <- .add(list(.mul(list(A, pre)), .div(B, denom)))
        return(.mr("RewriteRule", ig, x, rewritten = add, substep = .mr("AddRule", add, x, substeps = list(lstep, cstep))))
      }
      return(if (!is.null(lstep)) lstep else cstep)
    }
    coeffs <- P_all_coeffs(numer_poly)
    constant <- .sqrt_quadratic_denom_reduce(coeffs, a, b, c)[[3L]]
    i0 <- NULL
    if (!is_zero_num(constant)) i0 <- inverse_trig_rule(list(.inv(denom), x), S, degenerate = FALSE)
    .mr("SqrtQuadraticDenomRule", ig, x, a = a, b = b, c = c, coeffs = coeffs, i0_step = i0)
  }
  reduction_rule <- function(ig, nv, const) {
    k <- (-nv - 1L) %/% 2L
    delta <- .sub(.mul(list(.int(4), a, c)), .sq(b))
    R <- .add(list(.mul(list(c, .sq(x))), .mul(list(b, x)), a))
    double_root_step <- function() {
      sb <- .add(list(.mul(list(.sqrt(c), x)), .div(b, .mul(list(S_Two, .sqrt(c))))))
      nested <- new_sym("Pow", list(new_sym("Pow", list(sb, S_Two)), Rational(nv, 2L)))
      rewritten <- .mul(list(const, nested))
      .mr("RewriteRule", ig, x, rewritten = rewritten, substep = nested_pow_rule(list(rewritten, x), S))
    }
    if (isTRUE(ask_fact(delta, "zero"))) return(double_root_step())
    term_denom <- .mul(list(.int(2L * k - 1L), delta, .pow(R, Rational(2L * k - 1L, 2L))))
    constant_term <- .div(.mul(list(const, S_Two, .add(list(.mul(list(S_Two, c, x)), b)))), term_denom)
    coeff <- .div(.mul(list(.int(8), c, .int(k - 1L))), .mul(list(.int(2L * k - 1L), delta)))
    expr <- .mul(list(const, .pow(R, .sub(S_Half, .int(k)))))
    derive_expr <- new_sym("Derivative", list(constant_term, Tuple(x, S_One)))
    rewrite_expr <- .add(list(derive_expr, .mul(list(coeff, expr))))
    derive_step <- S$solve(derive_expr, x)
    if (is_zero_num(coeff)) substep <- derive_step
    else {
      nxt <- S$solve(expr, x)
      substep <- .mr("AddRule", rewrite_expr, x, substeps = list(derive_step, .mr("ConstantTimesRule", .mul(list(coeff, expr)), x, constant = coeff, other = expr, substep = nxt)))
    }
    nondeg <- .mr("RewriteRule", ig, x, rewritten = rewrite_expr, substep = substep)
    if (is.na(ask_fact(delta, "zero"))) return(.add_degenerate_step(Ne(delta, 0), nondeg, double_root_step()))
    nondeg
  }
  poly_reduction_rule <- function() {
    terms <- list(); steps <- list()
    rb <- .add(list(.mul(list(c, .sq(x))), .mul(list(b, x)), a))
    qr <- P_div(f_poly, .P(rb, x))
    quotient <- qr[[1L]]; rest <- qr[[2L]]
    if (!P_is_zero(quotient)) {
      qi <- .mul(list(P_expr(quotient), .pow(.sqrt(rb), .int(nv + 2L))))
      terms[[length(terms) + 1L]] <- qi
      steps[[length(steps) + 1L]] <- sqrt_quadratic_rule(list(qi, x), S, degenerate = FALSE)
    }
    if (!P_is_zero(rest)) {
      e <- P_nth(rest, 1L); d <- P_nth(rest, 0L)
      A <- .div(e, .mul(list(S_Two, c))); B <- .sub(d, .mul(list(A, b)))
      if (!is_zero_num(A)) {
        u <- Dummy("u")
        base <- .mul(list(.add(list(.mul(list(S_Two, c, x)), b)), .pow(.sqrt(rb), n)))
        term <- .mul(list(A, base))
        pr <- .mr("PowerRule", .pow(u, Rational(nv, 2L)), u, base = u, exp = Rational(nv, 2L))
        ustep <- .mr("URule", base, x, u_var = u, u_func = rb, substep = pr)
        if (!is_one_num(A)) ustep <- .mr("ConstantTimesRule", term, x, constant = A, other = base, substep = ustep)
        terms[[length(terms) + 1L]] <- term; steps[[length(steps) + 1L]] <- ustep
      }
      if (!is_zero_num(B)) {
        term <- .mul(list(B, .pow(.sqrt(rb), n)))
        terms[[length(terms) + 1L]] <- term; steps[[length(steps) + 1L]] <- reduction_rule(term, nv, B)
      }
    }
    rewritten <- .add(terms)
    substep <- if (length(steps) == 1L) steps[[1L]] else .mr("AddRule", rewritten, x, substeps = steps)
    .mr("RewriteRule", integrand, x, rewritten = rewritten, substep = substep)
  }
  if (nv > 0L) {
    numer_poly <- P_mul(f_poly, P_pow(.P(.add(list(a, .mul(list(b, x)), .mul(list(c, .sq(x))))), x), (nv + 1L) %/% 2L))
    rewritten <- .div(P_expr(numer_poly), .sqrt(.add(list(a, .mul(list(b, x)), .mul(list(c, .sq(x)))))))
    generic <- .mr("RewriteRule", integrand, x, rewritten = rewritten, substep = denom_rule(numer_poly, rewritten))
  } else if (nv == -1L) generic <- denom_rule(f_poly, integrand)
  else if (P_deg(f_poly) <= 0L) generic <- reduction_rule(integrand, nv, f)
  else generic <- poly_reduction_rule()
  step <- .add_degenerate_step(cond, generic, dstep)
  if (!.mi_eq(integrand, starting)) return(.mr("RewriteRule", starting, x, rewritten = integrand, substep = step))
  step
}

chebyshev_substitution_rule <- function(integral, S) {
  integrand <- integral[[1L]]; x <- integral[[2L]]
  if (.mi_is_ratfun(integrand, x)) return(NULL)
  c <- S_One; m <- S_Zero; ref <- NULL; p <- S_Zero
  for (factor in .mi_make_args_mul(integrand)) {
    if (!.has1(factor, x)) { c <- .mul(list(c, factor)); next }
    be <- .as_base_exp(factor); base <- be[[1L]]; ex <- be[[2L]]
    if (!.mi_is_rational_num(ex)) return(NULL)
    if (.mi_eq(base, x)) { m <- .add(list(m, ex)); next }
    if (base$type != "Add") return(NULL)
    ia <- .mi_indep_add(base, x)
    ma <- ia[[1L]]; xterm <- ia[[2L]]
    if (isTRUE(ask_fact(ma, "zero")) || xterm$type == "Add") return(NULL)
    im <- .as_coeff_mul_dep(xterm, x)
    mb <- im[[1L]]; xpower <- .mul(im[[2L]])
    xb <- .as_base_exp(xpower)
    mn <- xb[[2L]]
    if (!.mi_eq(xb[[1L]], x) || !.mi_is_rational_num(mn) || is_zero_num(mn)) return(NULL)
    if (is.null(ref)) { ref <- base; aa <- ma; bb <- mb; nn <- mn; p <- ex; next }
    if (!.mi_eq(mn, nn)) return(NULL)
    ratio <- cancel(powsimp(.div(factor, .pow(ref, ex)), force = FALSE))
    if (.has1(ratio, x)) return(NULL)
    c <- .mul(list(c, ratio)); p <- .add(list(p, ex))
  }
  if (is.null(ref)) return(NULL)
  canonical <- .mul(list(c, .pow(x, m), .pow(ref, p)))
  denom_c <- .as_numer_denom(c)[[2L]]
  if (p$type == "Integer" && .int_val(p) >= 0L) {
    rewritten <- expand(canonical)
    if (.mi_eq(rewritten, integrand)) return(NULL)
    return(.mr("RewriteRule", integrand, x, rewritten = rewritten, substep = S$solve(rewritten, x)))
  }
  if (p$type == "Integer") {
    L <- .mi_lcm(.mi_q(m), .mi_q(nn))
    if (L == 1L) return(.mr("RewriteRule", integrand, x, rewritten = canonical, substep = S$solve(canonical, x)))
    u <- Dummy("u")
    transformed <- .mul(list(c, .int(L), .pow(u, .sub(.mul(list(.int(L), .add(list(m, S_One)))), S_One)),
                             .pow(.add(list(aa, .mul(list(bb, .pow(u, .mul(list(.int(L), nn))))))), p)))
    return(.mr("URule", integrand, x, u_var = u, u_func = .pow(x, Rational(1L, L)), substep = S$solve(transformed, u)))
  }
  P <- .mi_p(p); q <- .mi_q(p)
  r <- .div(.add(list(m, S_One)), nn)
  if (r$type == "Integer") {
    u <- Dummy("u")
    transformed <- .mul(list(c, .div(.int(q), .mul(list(bb, nn))), .pow(u, .int(P + q - 1L)), .pow(.div(.sub(.pow(u, .int(q)), aa), bb), .sub(r, S_One))))
    general <- .mr("URule", integrand, x, u_var = u, u_func = .pow(ref, Rational(1L, q)), substep = S$solve(transformed, u))
    if (.if_zero_implies_zero(bb, denom_c) || (.mi_neg(p) && .if_zero_implies_zero(bb, ref))) return(general)
    di <- subs(canonical, bb, S_Zero)
    dstep <- .mr("RewriteRule", integrand, x, rewritten = di, substep = S$solve(di, x))
    return(.mr("PiecewiseRule", integrand, x, subfunctions = list(list(dstep, Eq(bb, 0)), list(general, S_true))))
  }
  s <- .add(list(r, p))
  if (s$type == "Integer") {
    u <- Dummy("u")
    transformed <- .mul(list(.neg(c), .div(.int(q), .mul(list(aa, nn))), .pow(u, .int(P + q - 1L)), .pow(.div(aa, .sub(.pow(u, .int(q)), bb)), .add(list(s, S_One)))))
    u_func <- .mul(list(.pow(ref, Rational(1L, q)), .pow(x, .neg(.div(nn, .int(q))))))
    general <- .mr("URule", integrand, x, u_var = u, u_func = u_func, substep = S$solve(transformed, u))
    if (.if_zero_implies_zero(aa, denom_c) || (.mi_neg(p) && .if_zero_implies_zero(aa, ref))) return(general)
    di <- subs(canonical, aa, S_Zero)
    dstep <- .mr("RewriteRule", integrand, x, rewritten = di, substep = S$solve(di, x))
    return(.mr("PiecewiseRule", integrand, x, subfunctions = list(list(dstep, Eq(aa, 0)), list(general, S_true))))
  }
  NULL
}
.mi_lcm <- function(a, b) { g <- function(a, b) if (b == 0) a else g(b, a %% b); as.integer(a * b / g(a, b)) }
# Expr.as_independent(x, as_Add = True)
.mi_indep_add <- function(e, x) {
  ts <- .mi_make_args_add(e)
  dep <- vapply(ts, function(t) .has1(t, x), NA)
  list(.add(ts[!dep]), .add(ts[dep]))
}

# ---------------------------------------------------------------------------
# Bioche's rules

bioche_substitution <- function(integral, S) {
  integrand <- integral[[1L]]; x <- integral[[2L]]
  TRIG <- .TRIG_FUNS
  atoms <- .mi_sort_ordered(.mi_find_fns(integrand, TRIG))
  if (!length(atoms)) return(NULL)
  if (all(vapply(atoms, function(a) a$name %in% c("sin", "cos"), NA)) && .mi_is_poly(integrand, atoms)) return(NULL)
  if (!.mi_is_ratfun(integrand, atoms)) return(NULL)
  masked <- xreplace(integrand, lapply(atoms, function(a) list(a, Dummy())))
  if (.has1(masked, x)) return(NULL)
  td <- list()
  for (atom in atoms) {
    ap <- P_or_null(atom$args[[1L]], x)
    if (is.null(ap) || P_deg(ap) > 1L) return(NULL)
    td[[length(td) + 1L]] <- list(atom, P_nth(ap, 1L), P_nth(ap, 0L))
  }
  ref <- NULL
  for (t in td) if (!is_zero_num(t[[2L]])) { ref <- t; break }
  if (is.null(ref)) return(NULL)
  coeff0 <- ref[[2L]]; phase0 <- ref[[3L]]
  ratios <- list(); residuals <- list()
  for (t in td) {
    ratio <- cancel(.div(t[[2L]], coeff0))
    if (!.mi_is_rational_num(ratio)) return(NULL)
    ratios[[length(ratios) + 1L]] <- ratio
    residuals[[length(residuals) + 1L]] <- cancel(.sub(t[[3L]], .mul(list(ratio, phase0))))
  }
  dl <- Reduce(.mi_lcm, vapply(ratios, .mi_q, 0L))
  omega <- cancel(.div(coeff0, .int(dl)))
  harmonics <- lapply(ratios, function(r) .mul(list(r, .int(dl))))
  phase_shift <- cancel(.div(phase0, .int(dl)))
  phase_base <- NULL
  for (r in residuals) if (!isTRUE(ask_fact(r, "zero"))) { phase_base <- r; break }
  phase <- NULL
  if (is.null(phase_base)) pm <- lapply(td, function(t) S_Zero)
  else {
    pr <- lapply(residuals, function(r) cancel(.div(r, phase_base)))
    if (!all(vapply(pr, .mi_is_rational_num, NA))) return(NULL)
    pl <- Reduce(.mi_lcm, vapply(pr, .mi_q, 0L))
    phase <- cancel(.div(phase_base, .int(pl)))
    pm <- lapply(pr, function(r) .mul(list(r, .int(pl))))
    cands <- list(S_Zero)
    for (i in seq_along(pm)) for (j in seq_len(i - 1L)) {
      for (nd in list(list(.sub(pm[[i]], pm[[j]]), .sub(harmonics[[i]], harmonics[[j]])), list(.add(list(pm[[i]], pm[[j]])), .add(list(harmonics[[i]], harmonics[[j]]))))) {
        if (!is_zero_num(nd[[2L]])) {
          v <- .div(nd[[1L]], nd[[2L]])
          cands <- c(cands, list(floor(v), ceiling(v)))
        }
      }
    }
    cands <- .mi_sort_ordered(.uniq_syms(cands))
    score <- function(k) {
      devs <- sort(vapply(seq_along(pm), function(i) as.numeric(abs(.sub(pm[[i]], .mul(list(harmonics[[i]], k))))$value), 0), decreasing = TRUE)
      c(devs, abs(as.numeric(k$value)))
    }
    best <- cands[[1L]]; bs <- score(best)
    for (k in cands[-1L]) {
      sc <- score(k)
      lt <- FALSE
      for (z in seq_along(sc)) { if (sc[z] < bs[z]) { lt <- TRUE; break }; if (sc[z] > bs[z]) break }
      if (lt) { best <- k; bs <- sc }
    }
    phase_shift <- cancel(.add(list(phase_shift, .mul(list(best, phase)))))
    pm <- lapply(seq_along(pm), function(i) .sub(pm[[i]], .mul(list(harmonics[[i]], best))))
  }
  u_func <- .add(list(.mul(list(omega, x)), phase_shift))
  u <- Dummy("u"); v <- Dummy("v"); s <- Dummy("s"); c <- Dummy("c")
  reps <- lapply(seq_along(td), function(i) list(td[[i]][[1L]], fcall(td[[i]][[1L]]$name, list(.add(list(.mul(list(harmonics[[i]], u)), .mul(list(pm[[i]], v))))))))
  expr_u <- cancel(expand_trig(xreplace(integrand, reps)))
  trig_reps <- function(angle, sv, cv) list(list(sin(angle), sv), list(cos(angle), cv), list(fcall("tan", list(angle)), .div(sv, cv)),
                                            list(fcall("cot", list(angle)), .div(cv, sv)), list(fcall("sec", list(angle)), .inv(cv)), list(fcall("csc", list(angle)), .inv(sv)))
  rewrite_even <- function(expr, squares, reps) {
    nd <- .as_numer_denom(cancel(expr))
    gens <- lapply(squares, `[[`, 1L)
    pr <- reps
    parts <- list()
    for (part in nd) {
      dict <- .mi_poly_dict(part, gens)
      if (is.null(dict)) return(NULL)
      for (mono in dict$monoms) {
        if (any(mono %% 2L != 0L)) return(NULL)
        for (z in seq_along(gens)) if (mono[z] > 0L)
          pr[[length(pr) + 1L]] <- list(.pow(gens[[z]], .int(mono[z])), .pow(squares[[z]][[2L]], .int(mono[z] %/% 2L)))
      }
      parts[[length(parts) + 1L]] <- dict$expr
    }
    pr <- .mi_dedupe_pairs(pr)
    cancel(.div(xreplace(parts[[1L]], pr), xreplace(parts[[2L]], pr)))
  }
  expanded <- expand_mul(expr_u)
  eatoms <- .mi_sort_ordered(.mi_find_fns(expanded, TRIG))
  polynomial <- all(vapply(eatoms, function(a) a$name %in% c("sin", "cos"), NA)) && (!length(eatoms) || .mi_is_poly(expanded, eatoms))
  separable <- !is.null(phase_base) && all(vapply(.mi_make_args_add(expanded), function(t) !.has1(.mul(.as_coeff_mul_dep(t, u)[[2L]]), v), NA))
  if (!is.null(phase_base)) expanded <- xreplace(expanded, list(list(v, phase)))
  termwise <- polynomial || separable
  parametrize_phase <- function(expr_u) {
    if (is.null(phase_base)) return(list(expr_u, list(), S_false, NULL))
    z <- Dummy("z")
    expr_sc <- xreplace(expr_u, trig_reps(v, s, c))
    methods <- list(
      list(expr_sc, list(list(s, .div(.sub(S_One, z), S_Two)), list(c, .div(.add(list(S_One, z)), S_Two))), list(), cos(.mul(list(S_Two, phase))), S_false),
      list(expr_sc, list(list(c, .sub(S_One, .sq(z)))), list(list(s, z)), sin(phase), S_false),
      list(expr_sc, list(list(s, .sub(S_One, .sq(z)))), list(list(c, z)), cos(phase), S_false),
      list(xreplace(expr_sc, list(list(s, .mul(list(z, c))))), list(list(c, .inv(.add(list(S_One, .sq(z)))))), list(), fcall("tan", list(phase)), Eq(cos(phase), 0)))
    found <- FALSE
    for (md in methods) {
      tr <- rewrite_even(md[[1L]], md[[2L]], md[[3L]])
      if (!is.null(tr)) { expr_u <- tr; psub <- list(list(z, md[[4L]])); sing <- md[[5L]]; found <- TRUE; break }
    }
    if (!found) {
      expr_u <- cancel(xreplace(expr_sc, list(list(s, .div(.mul(list(S_Two, z)), .add(list(S_One, .sq(z))))), list(c, .div(.sub(S_One, .sq(z)), .add(list(S_One, .sq(z))))))))
      psub <- list(list(z, fcall("tan", list(.div(phase, S_Two)))))
      sing <- Eq(cos(.div(phase, S_Two)), 0)
    }
    sing_u <- NULL
    if (!identical(sing$key, S_false$key)) {
      w <- Dummy("w")
      sing_u <- xreplace(cancel(xreplace(expr_u, list(list(z, .inv(w))))), list(list(w, S_Zero)))
      if (.has1(sing_u, S_ComplexInfinity) || .has1(sing_u, S_NaN)) return(NULL)
    }
    if (.is_true(sing)) { expr_u <- sing_u; psub <- list(); sing_u <- NULL }
    list(expr_u, psub, sing, sing_u)
  }
  psub <- list(); sing <- S_false; sing_u <- NULL
  if (!termwise) {
    pz <- parametrize_phase(expr_u)
    if (is.null(pz)) return(NULL)
    expr_u <- pz[[1L]]; psub <- pz[[2L]]; sing <- pz[[3L]]; sing_u <- pz[[4L]]
  }
  if (termwise) {
    if (isTRUE(ask_fact(cancel(.sub(u_func, x)), "zero"))) {
      ex <- xreplace(expanded, list(list(u, x)))
      generic <- .mr("RewriteRule", integrand, x, rewritten = ex, substep = S$solve(ex, x))
    } else generic <- .mr("URule", integrand, x, u_var = u, u_func = u_func, substep = S$solve(.div(expanded, omega), u))
  } else {
    t <- Dummy("t")
    expr_sc <- xreplace(expr_u, trig_reps(u, s, c))
    methods <- list(
      list(.neg(.div(expr_sc, .mul(list(.int(4), omega, s, c)))), list(list(s, .div(.sub(S_One, t), S_Two)), list(c, .div(.add(list(S_One, t)), S_Two))), list(), cos(.mul(list(S_Two, u_func)))),
      list(.div(expr_sc, .mul(list(omega, c))), list(list(c, .sub(S_One, .sq(t)))), list(list(s, t)), sin(u_func)),
      list(.neg(.div(expr_sc, .mul(list(omega, s)))), list(list(s, .sub(S_One, .sq(t)))), list(list(c, t)), cos(u_func)),
      list(xreplace(.div(.mul(list(expr_sc, .sq(c))), omega), list(list(s, .mul(list(t, c))))), list(list(c, .inv(.add(list(S_One, .sq(t)))))), list(), fcall("tan", list(u_func))))
    generic <- NULL
    for (md in methods) {
      tr <- rewrite_even(md[[1L]], md[[2L]], md[[3L]])
      if (is.null(tr)) next
      st <- .mr("URule", integrand, x, u_var = t, u_func = md[[4L]], substep = S$solve(tr, t))
      if (!mr_contains_dont_know(st)) { generic <- st; break }
    }
    if (is.null(generic)) {
      tr <- xreplace(expr_sc, list(list(s, .div(.mul(list(S_Two, t)), .add(list(S_One, .sq(t))))), list(c, .div(.sub(S_One, .sq(t)), .add(list(S_One, .sq(t)))))))
      tr <- cancel(.div(.mul(list(tr, S_Two)), .mul(list(omega, .add(list(S_One, .sq(t)))))))
      generic <- .mr("URule", integrand, x, u_var = t, u_func = fcall("tan", list(.div(u_func, S_Two))), substep = S$solve(tr, t))
    }
    if (!is.null(sing_u)) {
      ss <- .mr("URule", integrand, x, u_var = u, u_func = u_func, substep = S$solve(cancel(.div(sing_u, omega)), u))
      generic <- .mr("PiecewiseRule", integrand, x, subfunctions = list(list(ss, sing), list(generic, S_true)))
    }
    if (length(psub)) generic <- .mr("ReparameterizationRule", integrand, x, replacements = psub, substep = generic)
  }
  if (identical(ask_fact(omega, "zero"), FALSE)) return(generic)
  zr <- lapply(td, function(t) list(t[[1L]], fcall(t[[1L]]$name, list(t[[3L]]))))
  zi <- xreplace(integrand, zr)
  zstep <- .mr("RewriteRule", integrand, x, rewritten = zi, substep = S$solve(zi, x))
  .mr("PiecewiseRule", integrand, x, subfunctions = list(list(zstep, Eq(omega, 0)), list(generic, S_true)))
}
.mi_dedupe_pairs <- function(pr) pr[!duplicated(vapply(pr, function(p) p[[1L]]$key, ""))]
# monomials of a polynomial in gens (Poly(e, *gens).monoms()); NULL if not a
# polynomial
.mi_poly_dict <- function(e, gens) {
  e <- expand(e)
  if (!.mi_is_poly(e, gens)) return(NULL)
  monoms <- lapply(.mi_make_args_add(e), function(t) {
    vapply(gens, function(g) {
      tot <- 0L
      for (f in .mi_make_args_mul(t)) {
        be <- .as_base_exp(f)
        if (.mi_eq(be[[1L]], g)) tot <- tot + as.integer(.int_val(be[[2L]]))
        else if (.mi_eq(f, g)) tot <- tot + 1L
      }
      tot
    }, 0L)
  })
  list(monoms = monoms, expr = e)
}

# ---------------------------------------------------------------------------
# the solver

.mi_state_md <- new.env()
.mi_cache_dummy <- function() { if (is.null(.mi_state_md$z)) .mi_state_md$z <- Dummy("z"); .mi_state_md$z }

.mr_xreplace <- function(r, map) {
  if (inherits(r, "mrule")) {
    for (nm in setdiff(names(r), "type")) if (!is.null(r[[nm]])) r[[nm]] <- .mr_xreplace(r[[nm]], map)
    return(r)
  }
  if (is_sym(r)) return(xreplace(r, map))
  if (is.list(r)) { nms <- names(r); out <- lapply(r, .mr_xreplace, map = map); names(out) <- nms; return(out) }
  r
}

.mi_integral_key <- function(integrand, symbol) {
  if (!.has1(integrand, symbol)) return("Number")
  t <- integrand$type
  if (t %in% c("Symbol", "Dummy")) return("Symbol")
  if (t == "Function") {
    if (integrand$name %in% .TRIG_FUNS) return("TrigonometricFunction")
    if (integrand$name %in% .MI_ORTHO) return("OrthogonalPolynomial")
    return(integrand$name)
  }
  t
}

.mi_cond <- function(keys, rule) function(integral, S) {
  if (.mi_integral_key(integral[[1L]], integral[[2L]]) %in% keys) rule(integral, S) else NULL
}

.mi_switch <- function() list(
  Pow = list(power_rule, trig_rule, trig_powers_products_rule, inverse_trig_rule, quadratic_denom_rule, sqrt_quadratic_rule,
             sqrt_fractional_linear_rule, chebyshev_substitution_rule, euler_substitution_rule, bioche_substitution),
  Symbol = list(power_rule), exp = list(exp_rule), Add = list(add_rule),
  Mul = list(mul_rule, trig_powers_products_rule, heaviside_rule, quadratic_denom_rule, sqrt_quadratic_rule, sqrt_fractional_linear_rule,
             chebyshev_substitution_rule, euler_substitution_rule, trig_cmplx_exp_rule, trig_product_to_sum_rule, trig_poly_mul_rule, bioche_substitution),
  Derivative = list(derivative_rule), TrigonometricFunction = list(trig_rule), Heaviside = list(heaviside_rule),
  DiracDelta = list(dirac_delta_rule), OrthogonalPolynomial = list(orthogonal_poly_rule), Number = list(constant_rule))

.mi_strategy <- function(integral, S) {
  r <- special_function_rule(integral, S)
  if (!is.null(r)) return(r)
  key <- .mi_integral_key(integral[[1L]], integral[[2L]])
  for (rl in S$switch[[key]]) { r <- rl(integral, S); if (!is.null(r)) return(r) }
  r <- hyperbolic_rule(integral, S)
  if (!is.null(r)) return(r)
  r <- S$alternatives(integral, S)
  if (!is.null(r)) return(r)
  if (key %in% c("Mul", "Pow")) { r <- nested_pow_rule(integral, S); if (!is.null(r)) return(r) }
  .mr("DontKnowRule", integral[[1L]], integral[[2L]])
}

.mi_solver <- function(max_depth = 100L, branch = FALSE) {
  S <- new.env()
  S$max_depth <- max_depth; S$branch <- branch
  S$active <- character(); S$parts_u_count <- list(); S$solved <- new.env(hash = TRUE); S$starved <- FALSE
  S$switch <- .mi_switch()
  S$alternatives <- .mi_alternatives(list(
    rewrites_rule, substitution_rule,
    .mi_cond(c("Mul", "Pow"), partial_fractions_rule), .mi_cond(c("Mul", "Pow"), cancel_rule), .mi_cond("Mul", combine_power_rule),
    .mi_cond(c("Mul", "log", .INV_TRIG, .MI_SPECIAL_ERR), parts_rule), .mi_cond(c("Mul", "Pow"), distribute_expand_rule),
    trig_expand_rule), branch = branch)
  S$solve <- function(integrand, symbol) {
    cd <- .mi_cache_dummy()
    key <- xreplace(integrand, list(list(symbol, cd)))$key
    hit <- S$solved[[key]]
    if (!is.null(hit)) return(.mr_xreplace(hit, list(list(cd, symbol))))
    if (!is.null(S$max_depth) && length(S$active) >= S$max_depth) { S$starved <- TRUE; return(.mr("DontKnowRule", integrand, symbol)) }
    if (key %in% S$active) return(.mr("DontKnowRule", integrand, symbol))
    had_ancestor <- length(S$active) > 0L
    S$active <- c(S$active, key)
    outer_starved <- if (had_ancestor) S$starved else FALSE
    S$starved <- FALSE
    on.exit({ S$active <- S$active[S$active != key]; S$starved <- outer_starved || S$starved })
    result <- .mi_strategy(list(integrand, symbol), S)
    if (!S$starved) assign(key, .mr_xreplace(result, list(list(symbol, cd))), envir = S$solved)
    result
  }
  S$check_parts_u <- function(k) {
    n <- S$parts_u_count[[k]]
    if (is.null(n)) n <- 0L
    if (n > 2L) return(FALSE)
    S$parts_u_count[[k]] <- n + 1L
    TRUE
  }
  S
}

#' Integration by the manual (step-by-step) method
#'
#' `integral_steps(f, x)` returns the tree of integration rules that
#' SymPy's `manualintegrate` would apply; `manualintegrate(f, x)` evaluates
#' it to an antiderivative.
#'
#' @param f Expression.
#' @param var Integration variable.
#' @param ... Options of the IntegrationSolver (`max_depth`, `branch`).
#' @return A rule object (`integral_steps`) or a `sym` expression.
#' @export
integral_steps <- function(f, var, ...) {
  S <- .mi_solver(...)
  old <- options(expressions = max(getOption("expressions"), 50000L))
  on.exit(options(old))
  S$solve(as_sym(f), as_sym(var))
}

#' @rdname integral_steps
#' @export
manualintegrate <- function(f, var) {
  f <- as_sym(f); var <- as_sym(var)
  result <- mr_eval(integral_steps(f, var))
  if (result$type == "Piecewise" && length(result$args) == 2L) {
    cond <- result$args[[1L]]$args[[2L]]
    if (cond$type == "Equality" && .is_true(result$args[[2L]]$args[[2L]]))
      result <- Piecewise(list(result$args[[2L]]$args[[1L]], Ne(cond$args[[1L]], cond$args[[2L]])), list(result$args[[1L]]$args[[1L]], S_true))
  }
  has_erf_trig_mul <- function(expr) {
    for (sub in Filter(function(z) z$type == "Mul", .all_subexprs(expr)))
      if (.has_function(sub, c("erf", "erfc", "erfi")) && .has_function(sub, c("sin", "cos", "sinh", "cosh"))) return(TRUE)
    FALSE
  }
  if (has_erf_trig_mul(f) && has_erf_trig_mul(result)) result <- factor_terms(result)
  rac <- function(expr) {
    if (!.has1(expr, var)) return(S_Zero)
    if (expr$type == "Add") return(.add(lapply(expr$args, rac)))
    if (expr$type == "Mul") {
      sp <- .as_coeff_mul_dep(expr, var)
      dep <- .mul(sp[[2L]])
      if (dep$type == "Add") return(.mul(list(sp[[1L]], rac(dep))))
      return(expr)
    }
    expr
  }
  if (!.has_type(f, "Derivative")) result <- rac(result)
  result
}

#' @export
format.mrule <- function(x, ...) {
  fl <- setdiff(names(x), "type")
  parts <- vapply(fl, function(nm) {
    v <- x[[nm]]
    s <- if (is.null(v)) "None" else if (inherits(v, "mrule")) format(v)
         else if (is_sym(v)) .str(v)
         else if (is.list(v)) paste0("[", paste(vapply(v, function(z) if (inherits(z, "mrule")) format(z) else if (is_sym(z)) .str(z) else if (is.list(z)) paste0("(", paste(vapply(z, function(w) if (inherits(w, "mrule")) format(w) else .str(w), ""), collapse = ", "), ")") else paste(z), ""), collapse = ", "), "]")
         else paste(v)
    paste0(nm, "=", s)
  }, "")
  paste0(x$type, "(", paste(parts, collapse = ", "), ")")
}
#' @export
print.mrule <- function(x, ...) { cat(format(x), "\n"); invisible(x) }
