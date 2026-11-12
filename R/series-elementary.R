# taylor_term, _eval_as_leading_term and _eval_nseries of the trigonometric
# and hyperbolic functions and their inverses, after
# sympy/functions/elementary/trigonometric.py and hyperbolic.py.

.se_F <- function(name, ...) fcall(name, lapply(list(...), as_sym))
.se_fact <- function(k) fcall("factorial", list(.int(k)))
.se_x0 <- function(arg, x) cancel(subs(arg, x, S_Zero))
.se_isint <- function(v) isTRUE(ask_fact(v, "integer"))
.se_finite <- function(v) isTRUE(ask_fact(v, "finite"))
.se_im <- function(v) fcall("im", list(v))
.se_re <- function(v) fcall("re", list(v))
.se_neg <- function(v) isTRUE(ask_fact(v, "negative"))
.se_pos <- function(v) isTRUE(ask_fact(v, "positive"))
.se_cdir <- function(cdir) { cdir <- as_sym(cdir); if (is_zero_num(cdir)) S_One else cdir }
.se_limdir <- function(cdir) if (isTRUE(ask_fact(.se_re(as_sym(cdir)), "negative"))) "-" else "+"
.se_in <- function(v, l) any(vapply(l, function(z) identical(v$key, z$key), NA))
.se_rewrite_log_lt <- function(e, x, logx, cdir) .so_eval_alt(rewrite(e, "log"), x, logx, cdir)
.se_rewrite_log_ns <- function(e, x, n, logx, cdir) .ns_eval(rewrite(e, "log"), x, n, logx, cdir)

# ---------------------------------------------------------------------------
# taylor terms

.ns_taylor_hooks$sin <- function(k, a, prev, e) {
  if (k < 0L || k %% 2L == 0L) return(S_Zero)
  .div(.mul(list(.pow(S_NegativeOne, .int(k %/% 2L)), .pow(a, .int(k)))), .se_fact(k))
}
.ns_taylor_hooks$cos <- function(k, a, prev, e) {
  if (k < 0L || k %% 2L == 1L) return(S_Zero)
  .div(.mul(list(.pow(S_NegativeOne, .int(k %/% 2L)), .pow(a, .int(k)))), .se_fact(k))
}
.ns_taylor_hooks$tan <- function(k, a, prev, e) {
  if (k < 0L || k %% 2L == 0L) return(S_Zero)
  b <- .pow(S_Two, .int(k + 1L))
  .mul(list(.pow(S_NegativeOne, .int((k - 1L) %/% 2L)), b, .sub(b, S_One), .div(.se_F("bernoulli", k + 1L), .se_fact(k + 1L)), .pow(a, .int(k))))
}
.ns_taylor_hooks$cot <- function(k, a, prev, e) {
  if (k == 0L) return(.inv(a))
  if (k < 0L || k %% 2L == 0L) return(S_Zero)
  .mul(list(.pow(S_NegativeOne, .int((k + 1L) %/% 2L)), .pow(S_Two, .int(k + 1L)), .div(.se_F("bernoulli", k + 1L), .se_fact(k + 1L)), .pow(a, .int(k))))
}
.ns_taylor_hooks$sec <- function(k, a, prev, e) {
  if (k < 0L || k %% 2L == 1L) return(S_Zero)
  kk <- k %/% 2L
  .mul(list(.pow(S_NegativeOne, .int(kk)), .div(.se_F("euler", 2L * kk), .se_fact(2L * kk)), .pow(a, .int(2L * kk))))
}
.ns_taylor_hooks$csc <- function(k, a, prev, e) {
  if (k == 0L) return(.inv(a))
  if (k < 0L || k %% 2L == 0L) return(S_Zero)
  kk <- k %/% 2L + 1L
  .div(.mul(list(.pow(S_NegativeOne, .int(kk - 1L)), S_Two, .sub(.pow(S_Two, .int(2L * kk - 1L)), S_One), .se_F("bernoulli", 2L * kk), .pow(a, .int(2L * kk - 1L)))), .se_fact(2L * kk))
}
.se_asin_term <- function(k, a, sgn_alt) {
  kk <- (k - 1L) %/% 2L
  R <- .se_F("RisingFactorial", S_Half, kk)
  s <- if (sgn_alt) .pow(S_NegativeOne, .int(kk)) else S_One
  .div(.mul(list(s, R, .pow(a, .int(k)))), .mul(list(.se_fact(kk), .int(k))))
}
.ns_taylor_hooks$asin <- function(k, a, prev, e) if (k < 0L || k %% 2L == 0L) S_Zero else .se_asin_term(k, a, FALSE)
.ns_taylor_hooks$acos <- function(k, a, prev, e) {
  if (k == 0L) return(.div(S_Pi, S_Two))
  if (k < 0L || k %% 2L == 0L) return(S_Zero)
  .neg(.se_asin_term(k, a, FALSE))
}
.ns_taylor_hooks$atan <- function(k, a, prev, e) {
  if (k < 0L || k %% 2L == 0L) return(S_Zero)
  .div(.mul(list(.pow(S_NegativeOne, .int((k - 1L) %/% 2L)), .pow(a, .int(k)))), .int(k))
}
.ns_taylor_hooks$acot <- function(k, a, prev, e) {
  if (k == 0L) return(.div(S_Pi, S_Two))
  if (k < 0L || k %% 2L == 0L) return(S_Zero)
  .div(.mul(list(.pow(S_NegativeOne, .int((k + 1L) %/% 2L)), .pow(a, .int(k)))), .int(k))
}
.ns_taylor_hooks$sinh <- function(k, a, prev, e) if (k < 0L || k %% 2L == 0L) S_Zero else .div(.pow(a, .int(k)), .se_fact(k))
.ns_taylor_hooks$cosh <- function(k, a, prev, e) if (k < 0L || k %% 2L == 1L) S_Zero else .div(.pow(a, .int(k)), .se_fact(k))
.ns_taylor_hooks$tanh <- function(k, a, prev, e) {
  if (k < 0L || k %% 2L == 0L) return(S_Zero)
  b <- .pow(S_Two, .int(k + 1L))
  .mul(list(b, .sub(b, S_One), .div(.se_F("bernoulli", k + 1L), .se_fact(k + 1L)), .pow(a, .int(k))))
}
.ns_taylor_hooks$coth <- function(k, a, prev, e) {
  if (k == 0L) return(.inv(a))
  if (k < 0L || k %% 2L == 0L) return(S_Zero)
  .mul(list(.pow(S_Two, .int(k + 1L)), .div(.se_F("bernoulli", k + 1L), .se_fact(k + 1L)), .pow(a, .int(k))))
}
.ns_taylor_hooks$csch <- function(k, a, prev, e) {
  if (k == 0L) return(.inv(a))
  if (k < 0L || k %% 2L == 0L) return(S_Zero)
  .mul(list(S_Two, .sub(S_One, .pow(S_Two, .int(k))), .div(.se_F("bernoulli", k + 1L), .se_fact(k + 1L)), .pow(a, .int(k))))
}
.ns_taylor_hooks$sech <- function(k, a, prev, e) if (k < 0L || k %% 2L == 1L) S_Zero else .mul(list(.div(.se_F("euler", k), .se_fact(k)), .pow(a, .int(k))))
.ns_taylor_hooks$asinh <- function(k, a, prev, e) if (k < 0L || k %% 2L == 0L) S_Zero else .se_asin_term(k, a, TRUE)
.ns_taylor_hooks$acosh <- function(k, a, prev, e) {
  if (k == 0L) return(.div(.mul(list(S_I, S_Pi)), S_Two))
  if (k < 0L || k %% 2L == 0L) return(S_Zero)
  .neg(.mul(list(S_I, .se_asin_term(k, a, FALSE))))
}
.ns_taylor_hooks$atanh <- function(k, a, prev, e) if (k < 0L || k %% 2L == 0L) S_Zero else .div(.pow(a, .int(k)), .int(k))

# ---------------------------------------------------------------------------
# leading terms

# sin, cos, tan, cot, sec, csc: the zero (or pole) nearest the limit point
.se_trig_alt <- function(nm) { force(nm); function(e, x, logx, cdir) {
  arg <- e$args[[1L]]
  x0 <- .se_x0(arg, x)
  n <- switch(nm, sin = , csc = .div(x0, S_Pi), cos = , sec = .div(.add(list(x0, .div(S_Pi, S_Two))), S_Pi), tan = , cot = .div(.mul(list(S_Two, x0)), S_Pi))
  if (.se_isint(n)) {
    lt <- switch(nm,
      sin = , csc = as_leading_term(.sub(arg, .mul(list(n, S_Pi))), x),
      cos = , sec = as_leading_term(.add(list(.sub(arg, .mul(list(n, S_Pi))), .div(S_Pi, S_Two))), x),
      tan = , cot = as_leading_term(.sub(arg, .div(.mul(list(n, S_Pi)), S_Two)), x))
    return(switch(nm,
      sin = , cos = .mul(list(.pow(S_NegativeOne, n), lt)),
      csc = , sec = .div(.pow(S_NegativeOne, n), lt),
      tan = if (isTRUE(ask_fact(n, "even"))) lt else .neg(.inv(lt)),
      cot = if (isTRUE(ask_fact(n, "even"))) .inv(lt) else .neg(lt)))
  }
  if (identical(x0$key, S_ComplexInfinity$key)) x0 <- limit(arg, x, S_Zero, .se_limdir(cdir))
  if (is_const(x0, "oo") || is_const(x0, "-oo"))
    return(if (nm %in% c("sin", "cos")) AccumBounds(-1, 1) else AccumBounds(S_NegativeInfinity, S_Infinity))
  if (.se_finite(x0)) fcall(nm, list(x0)) else e
} }
for (.nm in c("sin", "cos", "tan", "cot", "sec", "csc")) .so_alt_hooks[[.nm]] <- .se_trig_alt(.nm)
rm(.nm)

.so_alt_hooks$sinh <- function(e, x, logx, cdir) {
  arg <- as_leading_term(e$args[[1L]], x, logx, cdir)
  arg0 <- subs(arg, x, S_Zero)
  if (is_nan(arg0)) arg0 <- limit(arg, x, S_Zero, if (.se_neg(as_sym(cdir))) "-" else "+")
  if (isTRUE(ask_fact(arg0, "zero"))) return(arg)
  if (.se_finite(arg0)) return(fcall("sinh", list(arg0)))
  e
}
.so_alt_hooks$cosh <- function(e, x, logx, cdir) {
  arg <- as_leading_term(e$args[[1L]], x, logx, cdir)
  arg0 <- subs(arg, x, S_Zero)
  if (is_nan(arg0)) arg0 <- limit(arg, x, S_Zero, if (.se_neg(as_sym(cdir))) "-" else "+")
  if (isTRUE(ask_fact(arg0, "zero"))) return(S_One)
  if (.se_finite(arg0)) return(fcall("cosh", list(arg0)))
  e
}
.so_alt_hooks$tanh <- function(e, x, logx, cdir) {
  arg <- as_leading_term(e$args[[1L]], x)
  if (.has1(arg, x) && isTRUE(.so_contains(.so_new(S_One, list(x)), arg))) return(arg)
  fcall("tanh", list(arg))
}
.so_alt_hooks$coth <- function(e, x, logx, cdir) {
  arg <- as_leading_term(e$args[[1L]], x)
  if (.has1(arg, x) && isTRUE(.so_contains(.so_new(S_One, list(x)), arg))) return(.inv(arg))
  fcall("coth", list(arg))
}
.so_alt_hooks$sech <- function(e, x, logx, cdir) .so_eval_alt(.inv(fcall("cosh", list(e$args[[1L]]))), x, logx, cdir)
.so_alt_hooks$csch <- function(e, x, logx, cdir) .so_eval_alt(.inv(fcall("sinh", list(e$args[[1L]]))), x, logx, cdir)

.so_alt_hooks$asin <- function(e, x, logx, cdir) {
  arg <- e$args[[1L]]; x0 <- .se_x0(arg, x)
  if (is_nan(x0)) return(fcall("asin", list(as_leading_term(arg, x))))
  if (isTRUE(ask_fact(x0, "zero"))) return(as_leading_term(arg, x))
  if (.se_in(x0, list(S_NegativeOne, S_One, S_ComplexInfinity))) return(expand(.se_rewrite_log_lt(e, x, logx, cdir)))
  if (.se_neg(.sub(S_One, .sq(x0)))) {
    ndir <- .so_dir(arg, x, .se_cdir(cdir))
    if (.se_neg(.se_im(ndir))) { if (.se_neg(x0)) return(.sub(.neg(S_Pi), fcall("asin", list(x0)))) }
    else if (.se_pos(.se_im(ndir))) { if (.se_pos(x0)) return(.sub(S_Pi, fcall("asin", list(x0)))) }
    else return(expand(.se_rewrite_log_lt(e, x, logx, cdir)))
  }
  fcall("asin", list(x0))
}
.so_alt_hooks$acos <- function(e, x, logx, cdir) {
  arg <- e$args[[1L]]; x0 <- .se_x0(arg, x)
  if (is_nan(x0)) return(fcall("acos", list(as_leading_term(arg, x))))
  if (is_one_num(x0)) return(.mul(list(.sqrt(S_Two), .sqrt(as_leading_term(.sub(S_One, arg), x)))))
  if (.se_in(x0, list(S_NegativeOne, S_ComplexInfinity))) return(.se_rewrite_log_lt(e, x, logx, cdir))
  if (.se_neg(.sub(S_One, .sq(x0)))) {
    ndir <- .so_dir(arg, x, .se_cdir(cdir))
    if (.se_neg(.se_im(ndir))) { if (.se_neg(x0)) return(.sub(.mul(list(S_Two, S_Pi)), fcall("acos", list(x0)))) }
    else if (.se_pos(.se_im(ndir))) { if (.se_pos(x0)) return(.neg(fcall("acos", list(x0)))) }
    else return(expand(.se_rewrite_log_lt(e, x, logx, cdir)))
  }
  fcall("acos", list(x0))
}
.so_alt_hooks$atan <- function(e, x, logx, cdir) {
  arg <- e$args[[1L]]; x0 <- .se_x0(arg, x)
  if (is_nan(x0)) return(fcall("atan", list(as_leading_term(arg, x))))
  if (isTRUE(ask_fact(x0, "zero"))) return(as_leading_term(arg, x))
  if (.se_in(x0, list(.neg(S_I), S_I, S_ComplexInfinity))) return(expand(.se_rewrite_log_lt(e, x, logx, cdir)))
  if (.se_neg(.add(list(S_One, .sq(x0))))) {
    ndir <- .so_dir(arg, x, .se_cdir(cdir))
    if (.se_neg(.se_re(ndir))) { if (.se_pos(.se_im(x0))) return(.sub(fcall("atan", list(x0)), S_Pi)) }
    else if (.se_pos(.se_re(ndir))) { if (.se_neg(.se_im(x0))) return(.add(list(fcall("atan", list(x0)), S_Pi))) }
    else return(expand(.se_rewrite_log_lt(e, x, logx, cdir)))
  }
  fcall("atan", list(x0))
}
.so_alt_hooks$acot <- function(e, x, logx, cdir) {
  arg <- e$args[[1L]]; x0 <- .se_x0(arg, x)
  if (is_nan(x0)) return(fcall("acot", list(as_leading_term(arg, x))))
  if (identical(x0$key, S_ComplexInfinity$key)) return(as_leading_term(.inv(arg), x))
  if (.se_in(x0, list(.neg(S_I), S_I, S_Zero))) return(expand(.se_rewrite_log_lt(e, x, logx, cdir)))
  if (isTRUE(ask_fact(x0, "imaginary")) && .se_pos(.add(list(S_One, .sq(x0))))) {
    ndir <- .so_dir(arg, x, .se_cdir(cdir))
    if (.se_pos(.se_re(ndir))) { if (.se_pos(.se_im(x0))) return(.add(list(fcall("acot", list(x0)), S_Pi))) }
    else if (.se_neg(.se_re(ndir))) { if (.se_neg(.se_im(x0))) return(.sub(fcall("acot", list(x0)), S_Pi)) }
    else return(expand(.se_rewrite_log_lt(e, x, logx, cdir)))
  }
  fcall("acot", list(x0))
}
.so_alt_hooks$asinh <- function(e, x, logx, cdir) {
  arg <- e$args[[1L]]; x0 <- .se_x0(arg, x)
  if (isTRUE(ask_fact(x0, "zero"))) return(as_leading_term(arg, x))
  if (is_nan(x0)) { ex <- fcall("asinh", list(as_leading_term(arg, x))); return(if (.se_finite(ex)) ex else e) }
  if (.se_in(x0, list(.neg(S_I), S_I, S_ComplexInfinity))) return(.se_rewrite_log_lt(e, x, logx, cdir))
  if (.se_neg(.add(list(S_One, .sq(x0))))) {
    ndir <- .so_dir(arg, x, .se_cdir(cdir))
    if (.se_pos(.se_re(ndir))) { if (.se_neg(.se_im(x0))) return(.sub(.neg(fcall("asinh", list(x0))), .mul(list(S_I, S_Pi)))) }
    else if (.se_neg(.se_re(ndir))) { if (.se_pos(.se_im(x0))) return(.add(list(.neg(fcall("asinh", list(x0))), .mul(list(S_I, S_Pi))))) }
    else return(.se_rewrite_log_lt(e, x, logx, cdir))
  }
  fcall("asinh", list(x0))
}
.so_alt_hooks$acosh <- function(e, x, logx, cdir) {
  arg <- e$args[[1L]]; x0 <- .se_x0(arg, x)
  if (.se_in(x0, list(S_NegativeOne, S_Zero, S_One, S_ComplexInfinity))) return(.se_rewrite_log_lt(e, x, logx, cdir))
  if (is_nan(x0)) { ex <- fcall("acosh", list(as_leading_term(arg, x))); return(if (.se_finite(ex)) ex else e) }
  if (.se_neg(.sub(x0, S_One))) {
    ndir <- .so_dir(arg, x, .se_cdir(cdir))
    if (.se_neg(.se_im(ndir))) {
      if (.se_neg(.add(list(x0, S_One)))) return(.sub(fcall("acosh", list(x0)), .mul(list(S_Two, S_I, S_Pi))))
      return(.neg(fcall("acosh", list(x0))))
    } else if (!.se_pos(.se_im(ndir))) return(.se_rewrite_log_lt(e, x, logx, cdir))
  }
  fcall("acosh", list(x0))
}
.so_alt_hooks$atanh <- function(e, x, logx, cdir) {
  arg <- e$args[[1L]]; x0 <- .se_x0(arg, x)
  if (isTRUE(ask_fact(x0, "zero"))) return(as_leading_term(arg, x))
  if (is_nan(x0)) { ex <- fcall("atanh", list(as_leading_term(arg, x))); return(if (.se_finite(ex)) ex else e) }
  if (.se_in(x0, list(S_NegativeOne, S_One, S_ComplexInfinity))) return(.se_rewrite_log_lt(e, x, logx, cdir))
  if (.se_neg(.sub(S_One, .sq(x0)))) {
    ndir <- .so_dir(arg, x, .se_cdir(cdir))
    if (.se_neg(.se_im(ndir))) { if (.se_neg(x0)) return(.sub(fcall("atanh", list(x0)), .mul(list(S_I, S_Pi)))) }
    else if (.se_pos(.se_im(ndir))) { if (.se_pos(x0)) return(.add(list(fcall("atanh", list(x0)), .mul(list(S_I, S_Pi))))) }
    else return(.se_rewrite_log_lt(e, x, logx, cdir))
  }
  fcall("atanh", list(x0))
}

# ---------------------------------------------------------------------------
# nseries

.se_logx_subs <- function(arg, x, logx) if (is.null(logx)) arg else subs(arg, fcall("log", list(x)), logx)
for (.nm in c("sin", "cos")) .ns_hooks[[.nm]] <- local({
  nm <- .nm
  function(e, x, n, logx, cdir) {
    arg <- .se_logx_subs(e$args[[1L]], x, logx)
    a0 <- subs(arg, x, S_Zero)
    if (.has1(a0, S_NaN) || .has1(a0, S_ComplexInfinity)) .pole_error(paste("Cannot expand", .str(e), "around 0"))
    .ns_function(e, x, n, logx, cdir)
  }
})
rm(.nm)
.ns_hooks$tan <- function(e, x, n, logx, cdir) {
  i <- .div(.mul(list(limit(e$args[[1L]], x, S_Zero), S_Two)), S_Pi)
  if (!is_zero_num(i) && i$type == "Integer") return(.ns_eval(rewrite(e, "cos"), x, n, logx))
  .ns_function(e, x, n, logx, 0)
}
.ns_hooks$cot <- function(e, x, n, logx, cdir) {
  i <- .div(limit(e$args[[1L]], x, S_Zero), S_Pi)
  if (!is_zero_num(i) && i$type == "Integer") return(.ns_eval(rewrite(e, "cos"), x, n, logx))
  .ns_eval(rewrite(e, "tan"), x, n, logx)
}
.ns_hooks$sec <- function(e, x, n, logx, cdir) .ns_eval(.inv(fcall("cos", list(e$args[[1L]]))), x, n, logx)
.ns_hooks$csc <- function(e, x, n, logx, cdir) .ns_eval(.inv(fcall("sin", list(e$args[[1L]]))), x, n, logx)
.ns_hooks$sinc <- function(e, x, n, logx, cdir) .ns_eval(.div(sin(e$args[[1L]]), e$args[[1L]]), x, n, logx)

# asin and acos at the branch points 1 and -1
.se_asinacos_ns <- function(nm) { force(nm); function(e, x, n, logx, cdir) {
  arg <- e$args[[1L]]
  arg0 <- subs(arg, x, S_Zero)
  for (pt in list(S_One, S_NegativeOne)) if (identical(arg0$key, pt$key)) {
    t <- Dummy("t", positive = TRUE)
    sgn <- if (is_one_num(pt)) S_One else S_NegativeOne
    ser <- .ns_series_driver(rewrite(fcall(nm, list(.sub(pt, .mul(list(sgn, .sq(t)))))), "log"), t, S_Zero, .mul(list(S_Two, n)), "+", NULL, 0)
    arg1 <- if (is_one_num(pt)) .sub(S_One, arg) else .add(list(S_One, arg))
    f <- as_leading_term(arg1, x)
    g <- .div(.sub(arg1, f), f)
    if (!isTRUE(.so_is_meromorphic(g, x, S_Zero))) {
      if (is_zero_num(n)) return(.so_new(S_One))
      const <- if (nm == "asin") (if (is_one_num(pt)) .div(S_Pi, S_Two) else .div(.neg(S_Pi), S_Two)) else (if (is_one_num(pt)) S_Zero else S_Pi)
      return(.add(list(const, .so_new(.sqrt(x), list(x)))))
    }
    res1 <- .ns_eval(.sqrt(.add(list(S_One, g))), x, n, logx)
    res <- expand(.mul(list(.so_removeO(res1), .sqrt(f))))
    return(.add(list(powsimp(expand(subs(.so_removeO(ser), t, res))), .so_new(.pow(x, n), list(x)))))
  }
  res <- .ns_function(e, x, n, logx, 0)
  if (identical(arg0$key, S_ComplexInfinity$key)) return(res)
  if (.se_neg(.sub(S_One, .sq(arg0)))) {
    ndir <- .so_dir(arg, x, .se_cdir(cdir))
    if (.se_neg(.se_im(ndir))) {
      if (.se_neg(arg0)) return(if (nm == "asin") .sub(.neg(S_Pi), res) else .sub(.mul(list(S_Two, S_Pi)), res))
    } else if (.se_pos(.se_im(ndir))) {
      if (.se_pos(arg0)) return(if (nm == "asin") .sub(S_Pi, res) else .neg(res))
    } else return(.se_rewrite_log_ns(e, x, n, logx, cdir))
  }
  res
} }
.ns_hooks$asin <- .se_asinacos_ns("asin")
.ns_hooks$acos <- .se_asinacos_ns("acos")
.ns_hooks$atan <- function(e, x, n, logx, cdir) {
  arg0 <- subs(e$args[[1L]], x, S_Zero)
  if (.se_in(arg0, list(S_I, .neg(S_I)))) return(.se_rewrite_log_ns(e, x, n, logx, cdir))
  res <- .ns_function(e, x, n, logx, 0)
  ndir <- .so_dir(e$args[[1L]], x, .se_cdir(cdir))
  if (identical(arg0$key, S_ComplexInfinity$key)) return(if (.is_true(Gt(.se_re(ndir), 0))) .sub(res, S_Pi) else res)
  if (.se_neg(.add(list(S_One, .sq(arg0))))) {
    if (.se_neg(.se_re(ndir))) { if (.se_pos(.se_im(arg0))) return(.sub(res, S_Pi)) }
    else if (.se_pos(.se_re(ndir))) { if (.se_neg(.se_im(arg0))) return(.add(list(res, S_Pi))) }
    else return(.se_rewrite_log_ns(e, x, n, logx, cdir))
  }
  res
}
.ns_hooks$acot <- function(e, x, n, logx, cdir) {
  arg0 <- subs(e$args[[1L]], x, S_Zero)
  if (.se_in(arg0, list(S_I, .neg(S_I)))) return(.se_rewrite_log_ns(e, x, n, logx, cdir))
  res <- .ns_function(e, x, n, logx, 0)
  if (identical(arg0$key, S_ComplexInfinity$key)) return(res)
  ndir <- .so_dir(e$args[[1L]], x, .se_cdir(cdir))
  if (isTRUE(ask_fact(arg0, "zero"))) return(if (.is_true(Lt(.se_re(ndir), 0))) .sub(res, S_Pi) else res)
  if (isTRUE(ask_fact(arg0, "imaginary")) && .se_pos(.add(list(S_One, .sq(arg0))))) {
    if (.se_pos(.se_re(ndir))) { if (.se_pos(.se_im(arg0))) return(.add(list(res, S_Pi))) }
    else if (.se_neg(.se_re(ndir))) { if (.se_neg(.se_im(arg0))) return(.sub(res, S_Pi)) }
    else return(.se_rewrite_log_ns(e, x, n, logx, cdir))
  }
  res
}
.ns_hooks$asinh <- function(e, x, n, logx, cdir) {
  arg <- e$args[[1L]]; arg0 <- subs(arg, x, S_Zero)
  if (.se_in(arg0, list(S_I, .neg(S_I)))) return(.se_rewrite_log_ns(e, x, n, logx, cdir))
  res <- .ns_function(e, x, n, logx, 0)
  if (identical(arg0$key, S_ComplexInfinity$key)) return(res)
  if (.se_neg(.add(list(S_One, .sq(arg0))))) {
    ndir <- .so_dir(arg, x, .se_cdir(cdir))
    if (.se_pos(.se_re(ndir))) { if (.se_neg(.se_im(arg0))) return(.sub(.neg(res), .mul(list(S_I, S_Pi)))) }
    else if (.se_neg(.se_re(ndir))) { if (.se_pos(.se_im(arg0))) return(.add(list(.neg(res), .mul(list(S_I, S_Pi))))) }
    else return(.se_rewrite_log_ns(e, x, n, logx, cdir))
  }
  res
}
.ns_hooks$acosh <- function(e, x, n, logx, cdir) {
  arg <- e$args[[1L]]; arg0 <- subs(arg, x, S_Zero)
  if (.se_in(arg0, list(S_One, S_NegativeOne))) return(.se_rewrite_log_ns(e, x, n, logx, cdir))
  res <- .ns_function(e, x, n, logx, 0)
  if (identical(arg0$key, S_ComplexInfinity$key)) return(res)
  if (.se_neg(.sub(arg0, S_One))) {
    ndir <- .so_dir(arg, x, .se_cdir(cdir))
    if (.se_neg(.se_im(ndir))) {
      if (.se_neg(.add(list(arg0, S_One)))) return(.sub(res, .mul(list(S_Two, S_I, S_Pi))))
      return(.neg(res))
    } else if (!.se_pos(.se_im(ndir))) return(.se_rewrite_log_ns(e, x, n, logx, cdir))
  }
  res
}
.ns_hooks$atanh <- function(e, x, n, logx, cdir) {
  arg <- e$args[[1L]]; arg0 <- subs(arg, x, S_Zero)
  if (.se_in(arg0, list(S_One, S_NegativeOne))) return(.se_rewrite_log_ns(e, x, n, logx, cdir))
  res <- .ns_function(e, x, n, logx, 0)
  if (identical(arg0$key, S_ComplexInfinity$key)) return(res)
  if (.se_neg(.sub(S_One, .sq(arg0)))) {
    ndir <- .so_dir(arg, x, .se_cdir(cdir))
    if (.se_neg(.se_im(ndir))) { if (.se_neg(arg0)) return(.sub(res, .mul(list(S_I, S_Pi)))) }
    else if (.se_pos(.se_im(ndir))) { if (.se_pos(arg0)) return(.add(list(res, .mul(list(S_I, S_Pi))))) }
    else return(.se_rewrite_log_ns(e, x, n, logx, cdir))
  }
  res
}
.ns_hooks$sech <- function(e, x, n, logx, cdir) .ns_eval(.inv(fcall("cosh", list(e$args[[1L]]))), x, n, logx)
.ns_hooks$csch <- function(e, x, n, logx, cdir) .ns_eval(.inv(fcall("sinh", list(e$args[[1L]]))), x, n, logx)

# Expr.is_meromorphic: the _eval_is_meromorphic methods of AtomicExpr, Add,
# Mul, Pow, Function and Piecewise; None is NA
.so_fuzzy_group <- function(vals) {
  saw_other <- FALSE
  for (v in vals) {
    if (isTRUE(v)) next
    if (is.na(v)) return(NA)
    if (saw_other) return(NA)
    saw_other <- TRUE
  }
  !saw_other
}
.so_is_meromorphic <- function(e, x, a) {
  if (!length(e$args) && !(e$type %in% c("Function"))) {
    if (e$type == "AccumulationBounds") return(FALSE)
    return(!is_number_type(e) && !(e$type == "Constant" && e$name %in% c("oo", "-oo", "zoo", "nan")) || isTRUE(ask_fact(e, "finite")))
  }
  if (e$type %in% c("Add", "Mul")) return(.so_fuzzy_group(lapply(e$args, .so_is_meromorphic, x = x, a = a)))
  if (e$type == "Pow") {
    bm <- .so_is_meromorphic(e$args[[1L]], x, a)
    ei <- e$args[[2L]]$type == "Integer"
    if (ei) return(bm)
    em <- .so_is_meromorphic(e$args[[2L]], x, a)
    if (isFALSE(bm)) return(if (isTRUE(em)) FALSE else NA)
    if (is.na(bm)) return(NA)
    b <- subs(e$args[[1L]], x, a)
    bz <- ask_fact(b, "zero")
    log_defined <- if (isTRUE(bz)) FALSE else {
      f <- ask_fact(b, "finite"); nz <- if (is.na(bz)) NA else !bz
      if (isFALSE(f) || isFALSE(nz)) FALSE else if (isTRUE(f) && isTRUE(nz)) TRUE else NA
    }
    if (isFALSE(log_defined)) return(ei)
    if (is.na(log_defined)) return(NA)
    if (!isTRUE(em)) return(em)
    return(ask_fact(subs(e$args[[2L]], x, a), "finite"))
  }
  if (e$type == "Function") {
    if (!length(e$args)) return(TRUE)
    if (any(vapply(e$args[-1L], function(z) .has1(z, x), NA))) return(FALSE)
    if (!isTRUE(.so_is_meromorphic(e$args[[1L]], x, a))) return(NA)
    sing <- .fn_is_singular(e$name, subs(e$args[[1L]], x, a))
    return(if (is.na(sing)) NA else !sing)
  }
  if (e$type == "Piecewise") {
    if (!isTRUE(ask_fact(a, "real"))) return(NA)
    for (p in e$args) {
      c <- p$args[[2L]]
      cond <- subs(c, x, a)
      if (cond$type %in% .REL_TYPES) return(NA)
      bd <- tryCatch(boundary(.conv_cond_set(c, x)), error = function(err) NULL)
      if (is.null(bd) || isTRUE(.is_true(contains(bd, a)))) return(NA)
      if (.is_true(cond)) return(.so_is_meromorphic(p$args[[1L]], x, a))
    }
    return(NA)
  }
  if (identical(e$key, x$key)) return(TRUE)
  if (!.has1(e, x)) return(TRUE)
  NA
}
# FunctionClass.is_singular with the _singularities of each class
.fn_is_singular <- function(name, a) {
  ss <- switch(name,
    gamma = , sin = , cos = , tan = , cot = , sec = , csc = , sinc = , exp = list(S_ComplexInfinity),
    asin = , acos = , asec = , acsc = list(S_One, S_NegativeOne, S_Zero, S_ComplexInfinity),
    atan = , acot = list(S_I, .neg(S_I)),
    log = list(S_Zero, S_ComplexInfinity),
    LambertW = list(.neg(.inv(S_E)), S_ComplexInfinity),
    Abs = , re = , im = , sign = , arg = , conjugate = , polar_lift = , periodic_argument = , principal_branch = TRUE,
    NULL)
  if (is.null(ss)) return(NA)
  if (isTRUE(ss)) return(TRUE)
  out <- FALSE
  for (s in ss) {
    v <- if (identical(s$key, S_ComplexInfinity$key)) ask_fact(a, "infinite") else ask_fact(.sub(a, s), "zero")
    if (isTRUE(v)) return(TRUE)
    if (is.na(v)) out <- NA
  }
  out
}

# asec, acsc, acoth, asech and acsch
.se_recip_RF <- function(k, n) {
  # R = RisingFactorial(1/2, k)*n and F = factorial(k)*n//2*n//2 (integer division)
  Fz <- ((gmp::factorialZ(k) * n) %/% 2L * n) %/% 2L
  list(.mul(list(.se_F("RisingFactorial", S_Half, k), .int(n))), .int(as.character(Fz)))
}
.ns_taylor_hooks$asec <- function(k, a, prev, e) {
  if (k == 0L) return(.mul(list(S_I, log(.div(S_Two, a)))))
  if (k < 0L || k %% 2L == 1L) return(S_Zero)
  rf <- .se_recip_RF(k %/% 2L, k)
  .div(.mul(list(.neg(S_I), rf[[1L]], .pow(a, .int(k)))), .mul(list(rf[[2L]], .int(4L))))
}
.ns_taylor_hooks$acsc <- function(k, a, prev, e) {
  if (k == 0L) return(.add(list(.div(S_Pi, S_Two), .neg(.mul(list(S_I, log(S_Two)))), .mul(list(S_I, log(a))))))
  if (k < 0L || k %% 2L == 1L) return(S_Zero)
  rf <- .se_recip_RF(k %/% 2L, k)
  .div(.mul(list(S_I, rf[[1L]], .pow(a, .int(k)))), .mul(list(rf[[2L]], .int(4L))))
}
.ns_taylor_hooks$acoth <- function(k, a, prev, e) {
  if (k == 0L) return(.div(.mul(list(.neg(S_I), S_Pi)), S_Two))
  if (k < 0L || k %% 2L == 0L) return(S_Zero)
  .div(.pow(a, .int(k)), .int(k))
}
.ns_taylor_hooks$asech <- function(k, a, prev, e) {
  if (k == 0L) return(log(.div(S_Two, a)))
  if (k < 0L || k %% 2L == 1L) return(S_Zero)
  rf <- .se_recip_RF(k %/% 2L, k)
  .div(.mul(list(S_NegativeOne, rf[[1L]], .pow(a, .int(k)))), .mul(list(rf[[2L]], .int(4L))))
}
.ns_taylor_hooks$acsch <- function(k, a, prev, e) {
  if (k == 0L) return(log(.div(S_Two, a)))
  if (k < 0L || k %% 2L == 1L) return(S_Zero)
  kk <- k %/% 2L
  rf <- .se_recip_RF(kk, k)
  .div(.mul(list(.pow(S_NegativeOne, .int(kk + 1L)), rf[[1L]], .pow(a, .int(k)))), .mul(list(rf[[2L]], .int(4L))))
}

.se_ndir <- function(arg, x, cdir) .so_dir(arg, x, .se_cdir(cdir))
.so_alt_hooks$asec <- function(e, x, logx, cdir) {
  arg <- e$args[[1L]]; x0 <- .se_x0(arg, x)
  if (is_nan(x0)) return(fcall("asec", list(as_leading_term(arg, x))))
  if (is_one_num(x0)) return(.mul(list(.sqrt(S_Two), .sqrt(as_leading_term(.sub(arg, S_One), x)))))
  if (.se_in(x0, list(S_NegativeOne, S_Zero))) return(.se_rewrite_log_lt(e, x, logx, cdir))
  if (isTRUE(ask_fact(x0, "real")) && .se_pos(.sub(S_One, .sq(x0)))) {
    ndir <- .se_ndir(arg, x, cdir)
    if (.se_neg(.se_im(ndir))) { if (.se_pos(x0)) return(.neg(fcall("asec", list(x0)))) }
    else if (.se_pos(.se_im(ndir))) { if (.se_neg(x0)) return(.sub(.mul(list(S_Two, S_Pi)), fcall("asec", list(x0)))) }
    else return(expand(.se_rewrite_log_lt(e, x, logx, cdir)))
  }
  fcall("asec", list(x0))
}
.so_alt_hooks$acsc <- function(e, x, logx, cdir) {
  arg <- e$args[[1L]]; x0 <- .se_x0(arg, x)
  if (is_nan(x0)) return(fcall("acsc", list(as_leading_term(arg, x))))
  if (.se_in(x0, list(S_NegativeOne, S_One, S_Zero))) return(expand(.se_rewrite_log_lt(e, x, logx, cdir)))
  if (identical(x0$key, S_ComplexInfinity$key)) return(as_leading_term(.inv(arg), x))
  if (isTRUE(ask_fact(x0, "real")) && .se_pos(.sub(S_One, .sq(x0)))) {
    ndir <- .se_ndir(arg, x, cdir)
    if (.se_neg(.se_im(ndir))) { if (.se_pos(x0)) return(.sub(S_Pi, fcall("acsc", list(x0)))) }
    else if (.se_pos(.se_im(ndir))) { if (.se_neg(x0)) return(.sub(.neg(S_Pi), fcall("acsc", list(x0)))) }
    else return(expand(.se_rewrite_log_lt(e, x, logx, cdir)))
  }
  fcall("acsc", list(x0))
}
.so_alt_hooks$acoth <- function(e, x, logx, cdir) {
  arg <- e$args[[1L]]; x0 <- .se_x0(arg, x)
  if (identical(x0$key, S_ComplexInfinity$key)) return(as_leading_term(.inv(arg), x))
  if (is_nan(x0)) { ex <- fcall("acoth", list(as_leading_term(arg, x))); return(if (.se_finite(ex)) ex else e) }
  if (.se_in(x0, list(S_NegativeOne, S_One, S_Zero))) return(.se_rewrite_log_lt(e, x, logx, cdir))
  if (isTRUE(ask_fact(x0, "real")) && .se_pos(.sub(S_One, .sq(x0)))) {
    ndir <- .se_ndir(arg, x, cdir)
    if (.se_neg(.se_im(ndir))) { if (.se_pos(x0)) return(.add(list(fcall("acoth", list(x0)), .mul(list(S_I, S_Pi))))) }
    else if (.se_pos(.se_im(ndir))) { if (.se_neg(x0)) return(.sub(fcall("acoth", list(x0)), .mul(list(S_I, S_Pi)))) }
    else return(.se_rewrite_log_lt(e, x, logx, cdir))
  }
  fcall("acoth", list(x0))
}
.so_alt_hooks$asech <- function(e, x, logx, cdir) {
  arg <- e$args[[1L]]; x0 <- .se_x0(arg, x)
  if (.se_in(x0, list(S_NegativeOne, S_Zero, S_One, S_ComplexInfinity))) return(.se_rewrite_log_lt(e, x, logx, cdir))
  if (is_nan(x0)) { ex <- fcall("asech", list(as_leading_term(arg, x))); return(if (.se_finite(ex)) ex else e) }
  if (.se_neg(x0) || .se_neg(.sub(S_One, x0))) {
    ndir <- .se_ndir(arg, x, cdir)
    if (.se_pos(.se_im(ndir))) {
      if (.se_pos(x0) || .se_neg(.add(list(x0, S_One)))) return(.neg(fcall("asech", list(x0))))
      return(.sub(fcall("asech", list(x0)), .mul(list(S_Two, S_I, S_Pi))))
    } else if (!.se_neg(.se_im(ndir))) return(.se_rewrite_log_lt(e, x, logx, cdir))
  }
  fcall("asech", list(x0))
}
.so_alt_hooks$acsch <- function(e, x, logx, cdir) {
  arg <- e$args[[1L]]; x0 <- .se_x0(arg, x)
  if (.se_in(x0, list(.neg(S_I), S_I, S_Zero))) return(.se_rewrite_log_lt(e, x, logx, cdir))
  if (is_nan(x0)) { ex <- fcall("acsch", list(as_leading_term(arg, x))); return(if (.se_finite(ex)) ex else e) }
  if (identical(x0$key, S_ComplexInfinity$key)) return(as_leading_term(.inv(arg), x))
  if (isTRUE(ask_fact(x0, "imaginary")) && .se_pos(.add(list(S_One, .sq(x0))))) {
    ndir <- .se_ndir(arg, x, cdir)
    if (.se_pos(.se_re(ndir))) { if (.se_pos(.se_im(x0))) return(.sub(.neg(fcall("acsch", list(x0))), .mul(list(S_I, S_Pi)))) }
    else if (.se_neg(.se_re(ndir))) { if (.se_neg(.se_im(x0))) return(.add(list(.neg(fcall("acsch", list(x0))), .mul(list(S_I, S_Pi))))) }
    else return(.se_rewrite_log_lt(e, x, logx, cdir))
  }
  fcall("acsch", list(x0))
}

# expansion at a branch point: f(base(t)) rewritten as log, t = sqrt of the
# leading term of arg1 (the _eval_nseries branches of asin, asec, asech, ...)
.se_branch_ns <- function(nm, e, x, n, logx, base, arg1, fallback) {
  t <- Dummy("t", positive = TRUE)
  ser <- .ns_series_driver(rewrite(fcall(nm, list(base(t))), "log"), t, S_Zero, .mul(list(S_Two, n)), "+", NULL, 0)
  f <- as_leading_term(arg1, x)
  g <- .div(.sub(arg1, f), f)
  if (!is.null(fallback) && !isTRUE(.so_is_meromorphic(g, x, S_Zero)))
    return(if (is_zero_num(as_sym(n))) .so_new(S_One) else .add(list(fallback, .so_new(.sqrt(x), list(x)))))
  res1 <- .ns_eval(.sqrt(.add(list(S_One, g))), x, n, logx)
  res <- expand(.mul(list(.so_removeO(res1), .sqrt(f))))
  .add(list(powsimp(expand(subs(.so_removeO(ser), t, res))), .so_new(.pow(x, n), list(x))))
}
.se_recip_ns <- function(nm) { force(nm); function(e, x, n, logx, cdir) {
  arg <- e$args[[1L]]; arg0 <- subs(arg, x, S_Zero)
  if (is_one_num(arg0)) return(.se_branch_ns(nm, e, x, n, logx, function(t) .add(list(S_One, .sq(t))), .add(list(S_NegativeOne, arg)), NULL))
  if (identical(arg0$key, S_NegativeOne$key)) return(.se_branch_ns(nm, e, x, n, logx, function(t) .sub(S_NegativeOne, .sq(t)), .sub(S_NegativeOne, arg), NULL))
  res <- .ns_function(e, x, n, logx, 0)
  if (identical(arg0$key, S_ComplexInfinity$key)) return(res)
  if (isTRUE(ask_fact(arg0, "real")) && .se_pos(.sub(S_One, .sq(arg0)))) {
    ndir <- .se_ndir(arg, x, cdir)
    if (.se_neg(.se_im(ndir))) { if (.se_pos(arg0)) return(if (nm == "asec") .neg(res) else .sub(S_Pi, res)) }
    else if (.se_pos(.se_im(ndir))) { if (.se_neg(arg0)) return(if (nm == "asec") .sub(.mul(list(S_Two, S_Pi)), res) else .sub(.neg(S_Pi), res)) }
    else return(.se_rewrite_log_ns(e, x, n, logx, cdir))
  }
  res
} }
.ns_hooks$asec <- .se_recip_ns("asec")
.ns_hooks$acsc <- .se_recip_ns("acsc")
.ns_hooks$acoth <- function(e, x, n, logx, cdir) {
  arg <- e$args[[1L]]; arg0 <- subs(arg, x, S_Zero)
  if (.se_in(arg0, list(S_One, S_NegativeOne))) return(.se_rewrite_log_ns(e, x, n, logx, cdir))
  res <- .ns_function(e, x, n, logx, 0)
  if (identical(arg0$key, S_ComplexInfinity$key)) return(res)
  if (isTRUE(ask_fact(arg0, "real")) && .se_pos(.sub(S_One, .sq(arg0)))) {
    ndir <- .se_ndir(arg, x, cdir)
    if (.se_neg(.se_im(ndir))) { if (.se_pos(arg0)) return(.add(list(res, .mul(list(S_I, S_Pi))))) }
    else if (.se_pos(.se_im(ndir))) { if (.se_neg(arg0)) return(.sub(res, .mul(list(S_I, S_Pi)))) }
    else return(.se_rewrite_log_ns(e, x, n, logx, cdir))
  }
  res
}
.ns_hooks$asech <- function(e, x, n, logx, cdir) {
  arg <- e$args[[1L]]; arg0 <- subs(arg, x, S_Zero)
  if (is_one_num(arg0)) return(.se_branch_ns("asech", e, x, n, logx, function(t) .sub(S_One, .sq(t)), .sub(S_One, arg), S_Zero))
  if (identical(arg0$key, S_NegativeOne$key)) return(.se_branch_ns("asech", e, x, n, logx, function(t) .add(list(S_NegativeOne, .sq(t))), .add(list(S_One, arg)), .mul(list(S_I, S_Pi))))
  res <- .ns_function(e, x, n, logx, 0)
  if (identical(arg0$key, S_ComplexInfinity$key)) return(res)
  if (.se_neg(arg0) || .se_neg(.sub(S_One, arg0))) {
    ndir <- .se_ndir(arg, x, cdir)
    if (.se_pos(.se_im(ndir))) {
      if (.se_pos(arg0) || .se_neg(.add(list(arg0, S_One)))) return(.neg(res))
      return(.sub(res, .mul(list(S_Two, S_I, S_Pi))))
    } else if (!.se_neg(.se_im(ndir))) return(.se_rewrite_log_ns(e, x, n, logx, cdir))
  }
  res
}
.ns_hooks$acsch <- function(e, x, n, logx, cdir) {
  arg <- e$args[[1L]]; arg0 <- subs(arg, x, S_Zero)
  if (identical(arg0$key, S_I$key)) return(.se_branch_ns("acsch", e, x, n, logx, function(t) .add(list(S_I, .sq(t))), .sub(arg, S_I), .div(.mul(list(.neg(S_I), S_Pi)), S_Two)))
  if (identical(arg0$key, .neg(S_I)$key)) return(.se_branch_ns("acsch", e, x, n, logx, function(t) .add(list(.neg(S_I), .sq(t))), .add(list(S_I, arg)), .div(.mul(list(S_I, S_Pi)), S_Two)))
  res <- .ns_function(e, x, n, logx, 0)
  if (identical(arg0$key, S_ComplexInfinity$key)) return(res)
  if (isTRUE(ask_fact(arg0, "imaginary")) && .se_pos(.add(list(S_One, .sq(arg0))))) {
    ndir <- .se_ndir(arg, x, cdir)
    if (.se_pos(.se_re(ndir))) { if (.se_pos(.se_im(arg0))) return(.sub(.neg(res), .mul(list(S_I, S_Pi)))) }
    else if (.se_neg(.se_re(ndir))) { if (.se_neg(.se_im(arg0))) return(.add(list(.neg(res), .mul(list(S_I, S_Pi))))) }
    else return(.se_rewrite_log_ns(e, x, n, logx, cdir))
  }
  res
}
