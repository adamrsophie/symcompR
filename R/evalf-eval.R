# sympy/core/evalf.py, continued: sums and products (hypsum, evalf_sum,
# evalf_prod), the _eval_evalf methods that the evalf fallback calls, and
# EvalfMixin.evalf, the top level of N().

# check_convergence(numer, denom, n): list(h, g, p)
.ev_check_convergence <- function(numer, denom, n) {
  npol <- Poly(numer, n); dpol <- Poly(denom, n)
  p <- degree(npol); q <- degree(dpol)
  rate <- q - p
  if (rate != 0) return(list(rate, NULL, NULL))
  constant <- .div(LC(dpol), LC(npol))
  one <- if (is_exact_number(constant)) abs(gmp::as.bigq(constant$value)) == 1
         else constant$type == "Float" && abs(constant$value) == 1
  if (!isTRUE(one)) return(list(rate, constant, NULL))
  if (p == 0 && q == 0) return(list(rate, constant, S_Zero))
  pc <- all_coeffs(npol)[[2L]]; qc <- all_coeffs(dpol)[[2L]]
  list(rate, constant, .div(.sub(qc, pc), LC(dpol)))
}

# lambdify(n, f) for a polynomial f in n: its value at the integer k as an
# integer (MPZ of the value, truncated)
.ev_int_poly <- function(f, n) {
  cf <- lapply(all_coeffs(Poly(f, n)), function(c) gmp::as.bigq(c$value))
  function(k) {
    v <- gmp::as.bigq(0)
    for (c in cf) v <- v * k + c
    q <- gmp::divq.bigz(abs(gmp::numerator(v)), gmp::denominator(v))
    if (v < 0) -q else q
  }
}

.ev_hypsum <- function(expr, n, start, prec) {
  if (!is.finite(prec)) .ev_nie("does not support inf prec")
  if (start != 0) expr <- subs(expr, n, .add(list(n, .int(start))))
  hs <- tryCatch(hypersimp(expr, n), error = function(err) NULL)
  if (is.null(hs)) .ev_nie("a hypergeometric series is required")
  nd <- .as_numer_denom(together(hs))
  num <- nd[[1L]]; den <- nd[[2L]]
  func1 <- .ev_int_poly(num, n); func2 <- .ev_int_poly(den, n)
  cc <- .ev_check_convergence(num, den, n)
  h <- cc[[1L]]; g <- cc[[2L]]; p <- cc[[3L]]
  if (h < 0) .ev_valueerror(sprintf("Sum diverges like (n!)^%i", as.integer(-h)))
  term <- subs(expr, n, S_Zero)
  if (!(term$type %in% c("Integer", "Rational")))
    .ev_nie("Non rational term functionality is not implemented.")
  tq <- gmp::as.bigq(term$value)
  gq <- if (!is.null(g)) gmp::as.bigq(g$value) else NULL
  # direct summation if geometric or faster
  if (h > 0 || (h == 0 && abs(gq) > 1)) {
    term <- .ev_floordiv(gmp::numerator(tq) * .pow2(prec), gmp::denominator(tq))
    s <- term
    k <- 1
    while (abs(term) > 5) {
      term <- .ev_floordiv(term * func1(k - 1), func2(k - 1))
      s <- s + term
      k <- k + 1
    }
    return(.mp_norm(list(m = s, e = -prec)))
  }
  alt <- gq < 0
  if (abs(gq) < 1)
    .ev_valueerror(sprintf("Sum diverges like (%i)^n", as.integer(as.numeric(abs(1 / gq)))))
  pv <- gmp::as.bigq(p$value)
  if (pv < 1 || (pv == 1 && !alt)) .ev_valueerror(sprintf("Sum diverges like n^%i", as.integer(as.numeric(-pv))))
  # polynomial convergence: Richardson extrapolation at quadruple precision,
  # repeated at doubled precision until two results agree to ndig digits
  vold <- NULL
  ndig <- .prec_to_dps(prec)
  old <- .MPX$prec
  on.exit(.MPX$prec <- old, add = TRUE)
  repeat {
    prec2 <- 4 * prec
    cur <- .ev_floordiv(gmp::numerator(tq) * .pow2(prec2), gmp::denominator(tq))
    summand <- function(k) {
      k <- .mv_dbl(k)
      if (k != 0) cur <<- .ev_floordiv(cur * func1(k - 1), func2(k - 1))
      .mv_mpf(.mp_norm(list(m = cur, e = -prec2)))
    }
    .MPX$prec <- prec
    v <- .mv_convert(.mp_nsum(summand, list(.mv_int(0), .mv_inf()), method = "richardson"))$re
    vf <- .mp_round(v$m, v$e, .dps_to_prec(ndig))
    if (!is.null(vold) && .mp_cmp(vold, vf) == 0) break
    prec <- prec + prec
    vold <- vf
  }
  v
}
# Python's // on integers (floor)
.ev_floordiv <- function(a, b) {
  a <- gmp::as.bigz(a); b <- gmp::as.bigz(b)
  q <- gmp::divq.bigz(abs(a), abs(b))
  if ((a < 0) != (b < 0)) { if (q * abs(b) != abs(a)) q <- q + 1; -q } else q
}

.ev_evalf_sum <- function(expr, prec, options) {
  if (!is.null(options$subs)) expr <- .ev_subs_all(expr, options$subs)
  func <- expr$args[[1L]]
  lims <- expr$args[-1L]
  if (length(lims) != 1L || length(lims[[1L]]$args) != 3L) .ev_nie()
  if (func$type == "Integer" && func$value == 0) return(.ev_res(NULL, NULL, prec, NULL))
  prec2 <- prec + 10
  n <- lims[[1L]]$args[[1L]]; a <- lims[[1L]]$args[[2L]]; b <- lims[[1L]]$args[[3L]]
  r <- tryCatch({
    if (!(b$type == "Constant" && b$name == "oo") || (a$type == "Constant" && a$name == "-oo") || a$type != "Integer")
      .ev_nie()
    a0 <- as.integer(as.character(a$value))
    v <- .ev_hypsum(func, n, a0, prec2)
    delta <- prec - .ev_fastlog(v)
    if (.ev_fastlog(v) < -10) v <- .ev_hypsum(func, n, a0, delta)
    .ev_res(if (v$m == 0) NULL else v, NULL, min(prec, delta), NULL)
  }, NotImplementedError = function(err) NULL)
  if (!is.null(r)) return(r)
  eps <- 2^(-prec)
  for (i in 1:4) {
    m <- 2^i * prec
    em <- euler_maclaurin(expr, m = m, n = m, eps = eps, eval_integral = FALSE)
    s <- em[[1L]]; err <- N(em[[2L]])
    if (is_nan(err)) .ev_nie()
    if (isTRUE(num_value_dbl(err) <= eps)) break
  }
  e2 <- .ev_evalf(fcall("Abs", list(err)), 20, options)
  err <- .ev_fastlog(e2$re)
  r <- .ev_evalf(s, prec2, options)
  if (.ev_is_zoo(r)) return(r)
  .ev_res(r$re, r$im, if (is.null(r$re_acc)) -err else r$re_acc, if (is.null(r$im_acc)) -err else r$im_acc)
}

.ev_evalf_prod <- function(expr, prec, options) {
  lims <- expr$args[-1L]
  if (all(vapply(lims, function(l) .sub(l$args[[2L]], l$args[[3L]])$type == "Integer", NA)))
    return(.ev_evalf(doit(expr), prec, options))
  # Product._eval_rewrite_as_Sum: exp(Sum(log(f), *limits))
  .ev_evalf(fcall("exp", list(.sum_new(fcall("log", list(expr$args[[1L]])), .lims(expr), "Sum"))),
            prec, options)
}

# ---------------------------------------------------------------------------
# _eval_evalf: the value of an expression at binary precision prec with the
# numbers it holds evaluated (NULL when the class has none)

.ev_is_number <- function(e) .expr_is_number(e)

# x._evalf(prec): _eval_evalf, or x itself
.ev_evalf_or_self <- function(x, prec) {
  r <- .ev_eval_evalf(x, prec)
  if (is.null(r)) x else r
}

.ev_eval_evalf <- function(e, prec) {
  switch(e$type,
    Integer = , Rational = .ev_float(.mp_from_q(e$value, prec), prec),
    Float = .ev_float(.float_mp(e), prec),
    Constant = .ev_const_evalf(e, prec),
    Symbol = , Dummy = NULL,
    Add = .ev_assocop_evalf(e, prec),
    Mul = .ev_mul_eval_evalf(e, prec),
    Pow = .ev_pow_eval_evalf(e, prec),
    Function = .ev_function_evalf(e, prec),
    CRootOf = .ev_crootof_evalf(e, prec),
    RootSum = .ev_rootsum_evalf(e, prec),
    AlgebraicNumber = .ev_eval_evalf(.an_value(e), prec),
    Piecewise = sym_rebuild(e, lapply(e$args, function(pc)
      sym_rebuild(pc, list(.ev_evalf_or_self(pc$args[[1L]], prec), pc$args[[2L]])))),
    {
      if (e$type %in% .REL_TYPES) return(sym_rebuild(e, lapply(e$args, .ev_evalf_or_self, prec)))
      NULL
    })
}

.ev_const_evalf <- function(e, prec) {
  switch(EXPR = e$name,
    pi = .ev_float(.ev_pi(prec), prec),
    E = .ev_float(.ev_e(prec), prec),
    oo = S_Infinity, `-oo` = S_NegativeInfinity, nan = S_NaN,
    I = NULL, zoo = NULL,
    TribonacciConstant = {
      t <- .div(.add(list(S_One, .pow(.add(list(.int(19), .neg(.mul(list(.int(3), .pow(.int(33), S_Half)))))), .rat(1, 3)),
                          .pow(.add(list(.int(19), .mul(list(.int(3), .pow(.int(33), S_Half))))), .rat(1, 3)))), .int(3))
      r <- .ev_evalf(t, prec + 4, .ev_options())
      .ev_float(r$re, prec)
    },
    {
      v <- .mp_const(e$name, prec + 10)
      if (is.null(v)) NULL else .ev_float(v, prec)
    })
}

# AssocOp._eval_evalf
.ev_assocop_evalf <- function(e, prec) {
  isdep <- function(a) length(free_symbols(a)) > 0L || .has_undef(a)
  args <- .sympy_args(e)
  if (e$type %in% c("Add", "Mul")) {
    dep <- vapply(args, isdep, NA)
    identity <- if (e$type == "Add") S_Zero else S_One
    build <- function(as) if (e$type == "Add") .add(as) else .mul(as)
    if (any(dep) && !all(dep)) {
      x <- build(args[!dep])
      x <- if (identical(x$key, identity$key)) identity else .ev_evalf_or_self(x, prec)
      tail_args <- args[dep]
      nargs <- lapply(tail_args, .ev_evalf_or_self, prec)
      return(build(c(list(x), nargs)))
    }
  }
  sym_rebuild(e, lapply(args, .ev_evalf_or_self, prec))
}

# Mul._eval_evalf
.ev_mul_eval_evalf <- function(e, prec) {
  cm <- .ev_as_coeff_Mul(e)
  c <- cm[[1L]]; m <- cm[[2L]]
  if (c$type == "Integer" && c$value == -1) {
    if (m$type == "Mul") rv <- .neg(.ev_assocop_evalf(m, prec))
    else rv <- .neg(.ev_evalf_or_self(m, prec))
  } else rv <- .ev_assocop_evalf(e, prec)
  if (.ev_is_number(rv)) return(expand(rv))
  rv
}

# Pow._eval_evalf
.ev_pow_eval_evalf <- function(e, prec) {
  base <- e$args[[1L]]; ex <- e$args[[2L]]
  if (base$type == "Constant" && base$name == "E") return(.ev_function_evalf(fnode("exp", list(ex)), prec))
  base <- .ev_evalf_or_self(base, prec)
  if (ex$type != "Integer") ex <- .ev_evalf_or_self(ex, prec)
  if (isTRUE(.fact(ex, "negative")) && .ev_is_number(base) && isFALSE(.fact(base, "extended_real"))) {
    cb <- conjugate(base)
    base <- .div(cb, .ev_evalf_or_self(.mul(list(base, cb)), prec))
    return(expand(.pow(base, .neg(ex))))
  }
  .pow(base, ex)
}

# Function._eval_evalf: the mpmath function of the arguments converted at
# prec + 5, computed at prec; NULL when an argument is not a number
.ev_function_evalf <- function(e, prec) {
  if (.is_undefined_application(e)) {
    imp <- e$extra$imp
    if (is.null(imp)) return(NULL)
    v <- tryCatch(imp(lapply(e$args, function(a) N(a, prec))), error = function(err) NULL)
    if (is.null(v)) return(NULL)
    v <- as_sym(v)
    if (!(v$type %in% .NUMBER_TYPES)) return(NULL)
    return(Float(v, dps = prec))
  }
  # MinMaxBase.evalf maps over the arguments, numbers or not
  if (e$name %in% c("Max", "Min")) return(sym_rebuild(e, lapply(e$args, .ev_evalf_or_self, prec)))
  if (e$name == "polar_lift") return(.ev_evalf_or_self(e$args[[1L]], prec))
  if (!all(vapply(e$args, .ev_is_number, NA))) return(NULL)
  if (e$name == "exp_polar") return(.ev_exp_polar_evalf(e, prec))
  # catalan._eval_evalf and mpmath's rf and ff: through gamma
  if (e$name %in% c("catalan", "RisingFactorial", "FallingFactorial")) {
    x <- e$args[[1L]]
    g <- switch(e$name,
      catalan = .mul(list(.pow(.int(4), x), fcall("gamma", list(.add(list(x, S_Half)))),
                          .pow(.mul(list(.sqrt(S_Pi), fcall("gamma", list(.add(list(x, .int(2))))))), S_NegativeOne))),
      RisingFactorial = .div(fcall("gamma", list(.add(list(x, e$args[[2L]])))), fcall("gamma", list(x))),
      FallingFactorial = .div(fcall("gamma", list(.add(list(x, S_One)))),
                              fcall("gamma", list(.add(list(x, .neg(e$args[[2L]]), S_One))))))
    if (e$name == "catalan") return(.ev_evalf_or_self(g, prec))
    r <- tryCatch(.ev_evalf(g, prec, .ev_options()), NotImplementedError = function(err) NULL)
    if (is.null(r) || .ev_is_zoo(r)) return(NULL)
    return(.ev_parts_sym(r, prec))
  }
  fe <- .ev_mp_funs(e$name)
  if (!is.null(fe)) {
    v <- tryCatch(fe(e, prec), NotImplementedError = function(err) NULL)
    if (is.null(v) || is_sym(v)) return(v)
    return(.mp_value_sym(v, prec))
  }
  if (!is.null(.MP_SYMPY[[e$name]])) {
    v <- tryCatch(.MP_SYMPY[[e$name]](e, prec, list()), NotImplementedError = function(err) NULL)
    if (is.null(v)) return(NULL)
    return(.mp_value_sym(v, prec))
  }
  # an integer sequence at moderate integer arguments (fibonacci, euler,
  # andre, ... held unevaluated): its exact value
  if (all(vapply(e$args, function(a) a$type == "Integer" && abs(gmp::as.bigz(a$value)) <= 10000, NA))) {
    d <- tryCatch(fcall(e$name, e$args), error = function(err) NULL)
    if (!is.null(d) && is_exact_number(d)) return(.ev_eval_evalf(d, prec))
  }
  h <- .FUNS[[e$name]]$evalf
  if (!is.null(h)) {
    r <- tryCatch(.ev_evalf(h(e), prec, .ev_options()), NotImplementedError = function(err) NULL)
    if (is.null(r)) return(NULL)
    if (.ev_is_zoo(r)) return(S_ComplexInfinity)
    return(.ev_parts_sym(r, prec))
  }
  args <- lapply(e$args, function(a) tryCatch(.ev_evalf(a, prec + 5, .ev_options()), NotImplementedError = function(err) NULL))
  if (any(vapply(args, is.null, NA)) || any(vapply(args, .ev_is_zoo, NA))) return(NULL)
  z0 <- function(a) if (is.null(a)) .mp_zero else a
  if (length(args) != 1L) return(NULL)
  a <- args[[1L]]
  if (.ev_sp(a$re) || .ev_sp(a$im)) return(NULL)
  z <- list(re = z0(a$re), im = z0(a$im))
  nm <- e$name
  # mpmath's asech(z) = acosh(1/z) and acsch(z) = asinh(1/z)
  if (nm %in% c("asech", "acsch")) {
    if (z$re$m == 0 && z$im$m == 0) return(NULL)
    z <- .mpc_div(.mpc(.mpf(1, 0)), z, prec + 20)
    nm <- if (nm == "asech") "acosh" else "asinh"
  }
  v <- tryCatch(.mpc_fun(nm, z, prec), error = function(err) NULL)
  if (is.null(v)) return(NULL)
  .mpc_sym(v, prec, zero_int = TRUE)
}

# exp_polar._eval_evalf: exp of the argument when its imaginary part is in
# (-pi, pi], with a negative imaginary part from rounding near pi dropped
.ev_exp_polar_evalf <- function(e, prec) {
  i <- fcall("im", list(e$args[[1L]]))
  bad <- tryCatch(isTRUE(.ev_lt(i, .neg(S_Pi), TRUE)) || isTRUE(.ev_lt(S_Pi, i, FALSE)), error = function(err) TRUE)
  if (bad) return(e)
  res <- .ev_function_evalf(fnode("exp", e$args), prec)
  if (is.null(res)) return(NULL)
  if (isTRUE(.ev_lt(S_Zero, i, FALSE))) {
    ri <- .ev_float_parts(res)
    if (!is.null(ri) && !is.null(ri$im) && ri$im$m < 0) return(.ev_float(ri$re, prec))
  }
  res
}
# a < b (or a <= b with eq) for real numbers, by their values
.ev_lt <- function(a, b, eq) {
  d <- .ev_evalf(.sub(b, a), 53, .ev_options())
  v <- if (is.null(d$re)) 0 else .mp_sign(d$re)
  if (eq) v >= 0 else v > 0
}

# _eval_evalf of the number sequences of
# sympy/functions/combinatorial/numbers.py and of polygamma: the mpmath
# value at the context precision, a sym (zoo, nan), or NULL
.ev_mp_funs <- function(name) {
  switch(name,
    bernoulli = function(e, prec) {
      a <- .ev_mp_args(e$args, prec, 2L, .mv_int(1))
      if (is.null(a)) return(NULL)
      n <- a[[1L]]; x <- a[[2L]]
      .ev_wp(prec, {
        if (.hx_iszero(n)) .hx_one()
        else if (.mv_eq(n, 1)) .hx_sub(x, 0.5)
        else if (.hx_isint(n) && .mv_cmp(n, 0) >= 0) {
          if (.mv_eq(x, 1)) .ev_bernoulli_num(.mv_dbl(n)) else .zt_bernpoly(.mv_dbl(n), x)
        } else .hx_mul(.hx_neg(n), .zt_zeta(.hx_sub(1, n), x))
      })
    },
    harmonic = function(e, prec) {
      a <- .ev_mp_args(e$args, prec, 2L, .mv_int(1))
      if (is.null(a)) return(NULL)
      n <- a[[1L]]; m <- a[[2L]]
      if (.hx_isint(n) && .mv_cmp(n, 0) < 0) return(S_NaN)
      .ev_wp(prec, {
        if (.mv_eq(m, 1)) {
          if (.mv_convert(n)$k == "mpf") .zt_harmonic(n)
          else .hx_add(.mp_psi(0, .hx_add(n, 1)), .hx_euler())
        } else .hx_sub(.zt_zeta(m), .zt_zeta(m, .hx_add(n, 1)))
      })
    },
    euler = function(e, prec) {
      a <- .ev_mp_args(e$args, prec, 2L, NULL)
      if (is.null(a)) return(NULL)
      m <- a[[1L]]; x <- if (length(a) > 1L) a[[2L]] else NULL
      if (.hx_isint(m) && .mv_cmp(m, 0) >= 0 && is.null(x)) {
        d <- tryCatch(fcall("euler", list(.int(.mv_dbl(m)))), error = function(err) NULL)
        if (!is.null(d) && is_exact_number(d)) return(.ev_eval_evalf(d, prec))
        return(NULL)
      }
      .ev_wp(prec, {
        if (.hx_isint(m) && .mv_cmp(m, 0) >= 0) .ev_eulerpoly(.mv_dbl(m), x)
        else if (.mv_eq(m, -1)) {
          if (is.null(x)) .hx_mul(.hx_pi(), .hx_pow(2, m))
          else .hx_sub(.mp_psi(0, .hx_div(.hx_add(x, 1), 2)), .mp_psi(0, .hx_div(x, 2)))
        } else {
          y <- if (is.null(x)) .mv_mpf(.mpf(1, -1)) else x
          r <- .hx_mul(2, .hx_sub(.zt_zeta(.hx_neg(m), y),
                                  .hx_mul(.hx_pow(2, .hx_add(m, 1)), .zt_zeta(.hx_neg(m), .hx_div(.hx_add(y, 1), 2)))))
          if (is.null(x)) .hx_mul(r, .hx_pow(2, m)) else r
        }
      })
    },
    genocchi = function(e, prec) {
      n <- e$args[[1L]]
      x <- if (length(e$args) > 1L) e$args[[2L]] else S_One
      B <- function(a, b) fnode("bernoulli", if (is.null(b)) list(a) else list(a, b))
      g <- if (identical(x$key, S_One$key) && isTRUE(.fact(n, "integer")) && isTRUE(.fact(n, "nonnegative")))
        .mul(list(S_Two, .sub(S_One, .pow(S_Two, n)), B(n, NULL)))
      else .mul(list(S_Two, .sub(B(n, x), .mul(list(.pow(S_Two, n), B(n, .div(.add(list(x, S_One)), S_Two)))))))
      .ev_evalf_or_self(g, prec)
    },
    andre = function(e, prec) {
      s <- .ev_to_mpmath(e$args[[1L]], prec + 12)
      if (is.null(s)) return(NULL)
      .ev_wp(prec + 12, {
        s2 <- .hx_div(s, 2)
        sp <- .hx_sinpi(s2); cp <- .hx_cospi(s2)
        .hx_mul(2, .ev_dirichlet(.hx_neg(s), list(.hx_neg(sp), cp, sp, .hx_neg(cp))))
      })
    },
    polygamma = function(e, prec) {
      s <- .ev_to_mpmath(e$args[[1L]], prec + 12)
      z <- .ev_to_mpmath(e$args[[2L]], prec + 12)
      if (is.null(s) || is.null(z)) return(NULL)
      if (.hx_isint(z) && .mv_cmp(z, 0) <= 0) return(S_ComplexInfinity)
      .ev_wp(prec + 12, {
        if (.hx_isint(s) && .mv_cmp(s, 0) >= 0) .mp_psi(.mv_dbl(s), z)
        else {
          zt <- .zt_zeta(.hx_add(s, 1), z)
          dzt <- .zt_zeta(.hx_add(s, 1), z, 1)
          .hx_mul(.hx_add(dzt, .hx_mul(.hx_add(.hx_euler(), .mp_psi(0, .hx_neg(s))), zt)), .hx_rgamma(.hx_neg(s)))
        }
      })
    },
    NULL)
}
# mpmath's _hurwitz_em for a derivative order d >= 1: the truncated
# L-series of (-log(a + k))^d*(a + k)^(-s) and the Euler-Maclaurin tail
.ev_hurwitz_em_d <- function(s, a, d, prec) {
  a <- .mv_convert(a)
  tol <- -prec
  M1 <- 0; M2 <- prec %/% 3; N <- M2
  lsum <- .mv_int(0)
  if (.hx_isint(s)) s <- .mv_int(.hx_int_value(.hx_re(s)))
  s1 <- .hx_sub(s, 1)
  negs <- .hx_neg(s)
  repeat {
    terms <- lapply(seq_len(M2 - M1) - 1, function(k) {
      w <- .hx_add(.hx_add(M1, a), k)
      .hx_mul(.hx_pow(.hx_log(w), d), .hx_pow(w, negs))
    })
    l <- .hx_fsum(terms)
    if (d %% 2 == 1) l <- .hx_neg(l)
    lsum <- .hx_add(lsum, l)
    M2a <- .hx_add(M2, a)
    logM2a <- .hx_log(M2a)
    logs <- list(.hx_pow(logM2a, d))
    logr <- .hx_div(1, logM2a)
    rM2a <- .hx_div(1, M2a)
    M2as <- .hx_pow(M2a, negs)
    # gammainc(d + 1, s1*log(M2a))/s1^(d + 1); for d = 1 the upper incomplete
    # gamma is (x + 1)*exp(-x)
    x <- .hx_mul(s1, logM2a)
    G <- if (d == 1) .hx_mul(.hx_add(x, 1), .hx_pow(M2a, .hx_neg(s1))) else .ev_nie("higher derivatives")
    tailsum <- .hx_div(G, .hx_pow(s1, d + 1))
    tailsum <- .hx_add(tailsum, .hx_mul(.hx_mul(0.5, logs[[1L]]), M2as))
    U <- list(.mv_int(1))
    r <- M2as
    fact <- gmp::as.bigz(2)
    for (j in seq_len(N)) {
      j2 <- 2 * j
      upds <- if (j == 1) 1 else c(j2 - 2, j2 - 1)
      for (m in upds) {
        D <- min(m, d + 1)
        if (m <= d) logs[[length(logs) + 1L]] <- .hx_mul(logs[[length(logs)]], logr)
        Un <- rep(list(.mv_int(0)), D + 1)
        for (i in seq_len(D)) Un[[i]] <- .hx_mul(.hx_sub(1 - m, s), U[[i]])
        for (i in seq_len(D)) if (i + 1 <= length(Un) && i <= length(U))
          Un[[i + 1]] <- .hx_add(Un[[i + 1]], .hx_mul(d - (i - 1), U[[i]]))
        U <- Un
        r <- .hx_mul(r, rM2a)
      }
      dot <- .hx_fsum(lapply(seq_len(min(length(U), length(logs))), function(i) .hx_mul(U[[i]], logs[[i]])))
      t <- .hx_div(.hx_mul(.hx_mul(dot, r), .zt_bernoulli(j2)), .mv_int(-fact))
      tailsum <- .hx_add(tailsum, t)
      if (.hx_mag(t) < tol) return(list(lsum, if (d %% 2 == 1) .hx_neg(tailsum) else tailsum))
      fact <- fact * (j2 + 1) * (j2 + 2)
    }
    M1 <- M2; M2 <- M2 * 2
    if (.mv_cmp(.hx_re(s), 0) < 0) N <- N + N %/% 2
  }
}

# the arguments as mpmath values at prec (_to_mpmath), padded to n with a
# default; NULL when one is not a number
.ev_mp_args <- function(args, prec, n, default) {
  if (!all(vapply(args, .ev_is_number, NA))) return(NULL)
  a <- lapply(args, .ev_to_mpmath, prec)
  if (any(vapply(a, is.null, NA))) return(NULL)
  if (length(a) < n && !is.null(default)) a[[n]] <- default
  a
}
.ev_wp <- function(p, expr) {
  old <- .MPX$prec
  on.exit(.MPX$prec <- old)
  .MPX$prec <- p
  force(expr)
}
# mpmath's bernoulli(n): exact below 5000, else 2*n!*zeta(n)/(2*pi)^n with
# the sign (-1)^(n/2 + 1)
.ev_bernoulli_num <- function(n) {
  if (n < 5000) return(.zt_bernoulli(n))
  if (n %% 2 == 1) return(.hx_zero())
  wp <- .MPX$prec
  .MPX$prec <- wp + 20
  v <- .hx_div(.hx_mul(.hx_mul(2, .hx_gamma(n + 1)), .zt_zeta(n)), .hx_pow(.hx_mul(2, .hx_pi()), n))
  .MPX$prec <- wp
  v <- .mv_pos(v)
  if ((n / 2) %% 2 == 0) .hx_neg(v) else v
}
# mpmath's eulerpoly(n, z) through E_n(z) = 2*(B_(n+1)(z) - 2^(n+1)*B_(n+1)(z/2))/(n+1)
.ev_eulerpoly <- function(n, z) {
  if (n == 0) return(.hx_pow(z, 0))
  if (n == 1) return(.hx_sub(z, 0.5))
  if (n == 2) return(.hx_mul(z, .hx_sub(z, 1)))
  m <- n + 1
  wp <- .MPX$prec
  .MPX$prec <- wp + 20
  v <- .hx_div(.hx_mul(2, .hx_sub(.zt_bernpoly(m, z), .hx_mul(.hx_pow(2, m), .zt_bernpoly(m, .hx_div(z, 2))))), m)
  .MPX$prec <- wp
  .mv_pos(v)
}
# mpmath's dirichlet(s, chi): sum of chi[p mod q]*zeta(s, p/q), p = 1..q,
# over q^s
.ev_dirichlet <- function(s, chi) {
  q <- length(chi)
  wp <- .MPX$prec
  .MPX$prec <- wp + 10
  # the pole of zeta at s = 1 cancels unless every character value is 0 or 1
  if (isTRUE(.mv_eq(s, 1))) {
    have_pole <- TRUE
    for (x in chi) if (!.hx_iszero(x) && !isTRUE(.mv_eq(x, 1))) {
      have_pole <- FALSE
      h <- .mv_pos(.hx_eps())
      .MPX$prec <- .MPX$prec * 2
      s <- .hx_add(s, h)
    }
    if (have_pole) { .MPX$prec <- wp; return(.mv_inf()) }
  }
  z <- .hx_zero()
  for (p in seq_len(q)) {
    cp <- chi[[p %% q + 1L]]
    if (!.hx_iszero(cp)) z <- .hx_add(z, .hx_mul(cp, .zt_zeta(s, .mv_q(p, q))))
  }
  z <- .hx_div(z, .hx_pow(q, s))
  .MPX$prec <- wp
  .mv_pos(z)
}

# Float(re, prec) + Float(im, prec)*I from an evalf result
.ev_parts_sym <- function(r, prec) {
  re <- if (is.null(r$re)) S_Zero else .ev_float(r$re, prec)
  if (is.null(r$im)) return(re)
  .add(list(re, .mul(list(.ev_float(r$im, prec), S_I))))
}

# _to_mpmath(prec): an mpmath value of a number expression, NULL if none
.ev_to_mpmath <- function(x, prec) {
  r <- tryCatch(.ev_evalf(x, prec, .ev_options()), NotImplementedError = function(err) NULL,
                KeyError = function(err) NULL)
  if (is.null(r)) {
    v <- .ev_eval_evalf(x, prec)
    if (is.null(v)) return(NULL)
    parts <- .ev_float_parts(v)
    if (is.null(parts)) return(NULL)
    r <- .ev_res(parts$re, parts$im)
  }
  if (.ev_is_zoo(r)) return(NULL)
  if (.ev_sp(r$re) && is.null(r$im)) return(.mv_special(r$re$sp))
  if (.ev_sp(r$re) || .ev_sp(r$im)) return(NULL)
  if (!is.null(r$im)) return(.mv_mpc(if (is.null(r$re)) .mp_zero else r$re, r$im))
  .mv_mpf(if (is.null(r$re)) .mp_zero else r$re)
}

# ---------------------------------------------------------------------------
# CRootOf and RootSum

# CRootOf._eval_evalf: the root to prec bits. A real root is bisected on
# an exact sign change, a complex one refined by Newton's method; a root on
# the imaginary axis keeps an exact zero real part.
.ev_crootof_evalf <- function(e, prec) {
  f <- e$args[[1L]]
  x <- Symbol("x")
  k <- as.integer(as.character(e$args[[2L]]$value))
  rs <- .numeric_roots(f, x)
  z0 <- rs[[k + 1L]]
  cf <- lapply(.coeff_list(f, x), function(c) gmp::as.bigq(c$value))
  P <- prec + 10
  if (!is.complex(z0)) return(.ev_float(.ev_real_root(cf, z0, P), prec))
  if (abs(Re(z0)) < 1e-8 * Mod(z0)) {
    y <- .ev_imag_root(cf, Im(z0), P)
    if (!is.null(y)) return(.mul(list(.ev_float(y, prec), S_I)))
  }
  z <- .ev_newton_complex(cf, z0, P)
  .add(list(.ev_float(z$re, prec), .mul(list(.ev_float(z$im, prec), S_I))))
}
# the sign of the polynomial with rational coefficients cf (highest first)
# at the dyadic m * 2^e, exactly
.ev_poly_sign <- function(cf, m, e) {
  d <- gmp::as.bigz(1)
  for (c in cf) d <- d * gmp::denominator(c) / gmp::gcd(d, gmp::denominator(c))
  zc <- lapply(cf, function(c) gmp::numerator(c * d))
  n <- length(zc) - 1L
  if (e >= 0) {
    t <- m * .pow2(e); s <- gmp::as.bigz(0)
    for (c in zc) s <- s * t + c
  } else {
    # sum c_i m^(n-i) 2^(-e*i): the denominators 2^(-e*n) cleared
    sc <- .pow2(-e)
    s <- zc[[1L]]
    for (i in seq_len(n)) s <- s * m + zc[[i + 1L]] * sc^i
  }
  if (s > 0) 1 else if (s < 0) -1 else 0
}
.ev_real_root <- function(cf, r0, P) {
  # a bracket around the double approximation, widened until the sign changes
  mp <- .mp_from_dbl(r0)
  w <- max(abs(r0), 1e-300) * 1e-12
  sgn <- function(a) .ev_poly_sign(cf, a$m, a$e)
  repeat {
    d <- .mp_from_dbl(w)
    lo <- .mp_sub(mp, d, 200); hi <- .mp_add(mp, d, 200)
    slo <- sgn(lo); shi <- sgn(hi)
    if (slo == 0) return(lo)
    if (shi == 0) return(hi)
    if (slo != shi) break
    w <- w * 16
  }
  repeat {
    mid <- .mp_ldexp(.mp_add(lo, hi, Inf), -1)
    s <- sgn(mid)
    if (s == 0) return(mid)
    if (s == slo) lo <- mid else hi <- mid
    width <- .mp_sub(hi, lo, Inf)
    mag <- max(.mp_mag(lo), .mp_mag(hi))
    if (.mp_mag(width) < mag - P - 4) return(.mp_round(mid$m, mid$e, P + 4))
  }
}
# a root i*y of f with y real: f(i*y) = A(y) + i*B(y), y a root of gcd(A, B)
.ev_imag_root <- function(cf, y0, P) {
  n <- length(cf) - 1L
  y <- Symbol("y")
  A <- S_Zero; B <- S_Zero
  for (i in seq_along(cf)) {
    d <- n - i + 1L
    c <- make_rational(cf[[i]])
    t <- .mul(list(c, .pow(y, .int(d))))
    switch(d %% 4 + 1, A <- .add(list(A, t)), B <- .add(list(B, t)), A <- .sub(A, t), B <- .sub(B, t))
  }
  g <- tryCatch(gcd(A, B), error = function(err) NULL)
  if (is.null(g) || length(free_symbols(g)) == 0L) return(NULL)
  gc <- lapply(.coeff_list(expand(g), y), function(c) gmp::as.bigq(c$value))
  r <- tryCatch(.ev_real_root(gc, y0, P), error = function(err) NULL)
  r
}
.ev_newton_complex <- function(cf, z0, P) {
  cfm <- lapply(cf, function(c) .mpc(.mp_from_q(c, P + 40)))
  dcf <- lapply(seq_len(length(cf) - 1L), function(i) .mpc(.mp_from_q(cf[[i]] * (length(cf) - i), P + 40)))
  hor <- function(cs, z, w) { s <- .mpc(.mp_zero); for (c in cs) s <- .mpc_add(.mpc_mul(s, z, w), c, w); s }
  z <- .mpc(.mp_from_dbl(Re(z0)), .mp_from_dbl(Im(z0)))
  w <- 60
  repeat {
    w <- min(2 * w, P + 40)
    for (it in 1:3) {
      st <- .mpc_div(hor(cfm, z, w), hor(dcf, z, w), w)
      z <- .mpc_sub(z, st, w)
    }
    if (w >= P + 40) {
      mz <- max(.mp_mag(z$re), .mp_mag(z$im)); ms <- max(.mp_mag(st$re), .mp_mag(st$im))
      if (!is.finite(ms) || ms < mz - P - 8) return(z)
      if (w > 20 * P + 400) return(z)
    }
  }
}

# RootSum._eval_evalf: the function summed over the numerical roots
.ev_rootsum_evalf <- function(e, prec) {
  f <- e$args[[1L]]
  fs <- free_symbols(f)
  if (length(fs) != 1L) return(e)
  x <- fs[[1L]]
  d <- degree(f, x)
  lam <- e$args[[2L]]
  cf <- tryCatch(lapply(all_coeffs(Poly(f, x)), as_sym), error = function(err) NULL)
  rts <- if (!is.null(cf) && all(vapply(cf, function(c) c$type %in% c("Integer", "Rational"), NA)))
    .poly_nroots(cf, n = .prec_to_dps(prec))
  else {
    f0 <- subs(f, x, Symbol("x"))
    lapply(seq_len(d) - 1L, function(k) .ev_crootof_evalf(new_sym("CRootOf", list(expand(f0), .int(k))), prec))
  }
  .add(lapply(rts, function(r) subs(lam$args[[2L]], .rootsum_w, r)))
}

# ---------------------------------------------------------------------------
# EvalfMixin.evalf

.ev_N <- function(e, n = 15, subs = NULL, maxn = 100, chop = FALSE, strict = FALSE, quad = NULL, verbose = FALSE) {
  prec <- .dps_to_prec(n)
  options <- .ev_options(maxprec = max(prec, as.integer(maxn * .ev_LG10)), chop = chop,
                         strict = strict, verbose = verbose)
  if (!is.null(subs)) assign("subs", subs, envir = options)
  if (!is.null(quad)) assign("quad", quad, envir = options)
  result <- tryCatch(.ev_evalf(e, prec + 4, options), NotImplementedError = function(err) NULL)
  if (is.null(result)) {
    v <- if (length(subs)) .ev_eval_evalf(.ev_subs_all(e, subs), prec) else .ev_eval_evalf(e, prec)
    if (is.null(v)) return(structure(list(e), class = "ev_unevaluated"))
    if (!.ev_is_number(v)) return(v)
    result <- tryCatch(.ev_evalf(v, prec, options), NotImplementedError = function(err) NULL)
    if (is.null(result)) return(structure(list(v), class = "ev_unevaluated"))
  }
  if (.ev_is_zoo(result)) return(S_ComplexInfinity)
  re <- result$re; im <- result$im
  re <- if (!is.null(re)) .ev_float(re, max(min(prec, result$re_acc), 1)) else S_Zero
  if (!is.null(im)) {
    im <- .ev_float(im, max(min(prec, result$im_acc), 1))
    return(.add(list(re, .mul(list(im, S_I)))))
  }
  re
}
