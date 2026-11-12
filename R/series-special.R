# Series methods of the non-elementary functions and of Piecewise,
# Derivative, Subs and Integral, after the _eval_nseries,
# _eval_as_leading_term, _eval_aseries and taylor_term methods in
# sympy/functions/elementary/{complexes,integers,piecewise}.py,
# sympy/functions/special/{gamma_functions,error_functions}.py,
# sympy/functions/combinatorial/factorials.py, sympy/core/{mod,function}.py
# and sympy/integrals/integrals.py.

.sp_lim0 <- function(e, x, dir = "+") limit(e, x, S_Zero, dir)
.sp_is <- function(v, k) isTRUE(identical(v$key, k$key))
.sp_inf_pts <- function(p) .sp_is(p, S_Infinity) || .sp_is(p, S_NegativeInfinity)
.sp_finite <- function(v) isTRUE(ask_fact(v, "finite"))
.sp_infinite <- function(v) isTRUE(ask_fact(v, "infinite"))
.sp_zero <- function(v) isTRUE(ask_fact(v, "zero"))
.sp_redir <- function(cdir) if (isTRUE(ask_fact(.se_re(as_sym(cdir)), "negative"))) "-" else "+"
.sp_logx <- function(x, logx) if (is.null(logx)) fcall("log", list(x)) else logx
.sp_sum <- function(l) .add(l)
.sp_o <- function(e, x) .so_new(e, list(x))
.sp_nn <- function(n) as.integer(num_value_dbl(as_sym(n)))

# Function._eval_aseries dispatch
.ns_aseries <- function(e, n, args0, x, logx) {
  h <- .ns_aseries_hooks[[e$name]]
  if (!is.null(h)) return(h(e, n, args0, x, logx))
  .pole_error(paste0("Asymptotic expansion of ", e$name, " around [", paste(vapply(args0, .str, ""), collapse = ", "), "] is not implemented."))
}

# ---------------------------------------------------------------------------
# complexes: sign, Abs, arg

.ns_hooks$sign <- function(e, x, n, logx, cdir) {
  arg0 <- e$args[[1L]]
  x0 <- subs(arg0, x, S_Zero)
  if (!is_zero_num(x0)) return(fcall("sign", list(x0)))
  cdir <- as_sym(cdir)
  if (!is_zero_num(cdir)) cdir <- .so_dir(arg0, x, cdir)
  if (.is_true(Lt(.se_re(cdir), 0))) S_NegativeOne else S_One
}
.ns_hooks$Abs <- function(e, x, n, logx, cdir) {
  direction <- .so_leadterm(e$args[[1L]], x)[[1L]]
  lx <- fcall("log", list(x))
  if (.has1(direction, lx)) direction <- subs(direction, lx, logx)
  s <- .ns_eval(e$args[[1L]], x, n, logx)
  expand(.mul(list(fcall("sign", list(direction)), s)))
}
.so_alt_hooks$arg <- function(e, x, logx, cdir) {
  t <- Dummy("t", positive = TRUE)
  cdir <- as_sym(cdir)
  if (is_zero_num(cdir)) cdir <- S_One
  z <- subs(e$args[[1L]], x, .mul(list(cdir, t)))
  if (isTRUE(ask_fact(z, "positive"))) return(S_Zero)
  if (isTRUE(ask_fact(z, "negative"))) return(S_Pi)
  .pole_error(paste("Cannot expand", .str(e), "around 0"))
}
.ns_hooks$arg <- function(e, x, n, logx, cdir) {
  if (num_value_dbl(as_sym(n)) <= 0) return(.so_new(S_One))
  .so_alt_hooks$arg(e, x, logx, cdir)
}

# ---------------------------------------------------------------------------
# floor, ceiling, frac and Mod

.sp_round_lt <- function(nm) { force(nm); function(e, x, logx, cdir) {
  arg <- e$args[[1L]]
  arg0 <- subs(arg, x, S_Zero)
  r <- subs(e, x, S_Zero)
  if (is_nan(arg0) || arg0$type == "AccumBounds") {
    arg0 <- limit(arg, x, S_Zero, .sp_redir(cdir))
    r <- fcall(nm, list(arg0))
  }
  if (.sp_finite(arg0)) {
    if (.sp_is(arg0, r)) {
      ndir <- .so_dir(arg, x, if (is_zero_num(as_sym(cdir))) S_One else as_sym(cdir))
      if (isTRUE(ask_fact(ndir, "negative"))) return(if (nm == "floor") .sub(r, S_One) else r)
      if (isTRUE(ask_fact(ndir, "positive"))) return(if (nm == "floor") r else .add(list(r, S_One)))
      .notimpl_error(paste("Not sure of sign of", .str(ndir)))
    }
    return(r)
  }
  as_leading_term(arg, x, logx, cdir)
} }
.so_alt_hooks$floor <- .sp_round_lt("floor")
.so_alt_hooks$ceiling <- .sp_round_lt("ceiling")
.sp_round_ns <- function(nm) { force(nm); function(e, x, n, logx, cdir) {
  arg <- e$args[[1L]]
  arg0 <- subs(arg, x, S_Zero)
  r <- subs(e, x, S_Zero)
  if (is_nan(arg0)) {
    arg0 <- limit(arg, x, S_Zero, .sp_redir(cdir))
    r <- fcall(nm, list(arg0))
  }
  if (.sp_infinite(arg0)) {
    s <- .ns_eval(arg, x, n, logx, cdir)
    o <- if (num_value_dbl(as_sym(n)) <= 0) .so_new(S_One, list(x)) else if (nm == "floor") AccumBounds(-1, 0) else AccumBounds(0, 1)
    return(.add(list(s, o)))
  }
  if (.sp_is(arg0, r)) {
    ndir <- .so_dir(arg, x, if (is_zero_num(as_sym(cdir))) S_One else as_sym(cdir))
    if (isTRUE(ask_fact(ndir, "negative"))) return(if (nm == "floor") .sub(r, S_One) else r)
    if (isTRUE(ask_fact(ndir, "positive"))) return(if (nm == "floor") r else .add(list(r, S_One)))
    .notimpl_error(paste("Not sure of sign of", .str(ndir)))
  }
  r
} }
.ns_hooks$floor <- .sp_round_ns("floor")
.ns_hooks$ceiling <- .sp_round_ns("ceiling")
.so_alt_hooks$frac <- function(e, x, logx, cdir) {
  arg <- e$args[[1L]]
  arg0 <- subs(arg, x, S_Zero)
  r <- subs(e, x, S_Zero)
  if (.sp_finite(arg0)) {
    if (.sp_zero(r)) {
      ndir <- .so_dir(arg, x, cdir)
      if (isTRUE(ask_fact(ndir, "negative"))) return(S_One)
      return(as_leading_term(.sub(arg, arg0), x, logx, cdir))
    }
    return(r)
  }
  if (.se_in(arg0, list(S_ComplexInfinity, S_Infinity, S_NegativeInfinity))) return(AccumBounds(0, 1))
  as_leading_term(arg, x, logx, cdir)
}
.ns_hooks$frac <- function(e, x, n, logx, cdir) {
  arg <- e$args[[1L]]
  arg0 <- subs(arg, x, S_Zero)
  r <- subs(e, x, S_Zero)
  if (.sp_infinite(arg0)) {
    if (num_value_dbl(as_sym(n)) <= 0) return(.so_new(S_One, list(x)))
    return(.add(list(AccumBounds(0, 1), .so_new(.pow(x, n), list(x)))))
  }
  res <- .ns_eval(.sub(arg, arg0), x, n, logx, cdir)
  if (.sp_zero(r)) {
    ndir <- .so_dir(arg, x, cdir)
    .add(list(res, if (isTRUE(ask_fact(ndir, "negative"))) S_One else S_Zero))
  } else .add(list(res, r))
}
.sp_mod_floor <- function(e) .sub(e$args[[1L]], .mul(list(e$args[[2L]], fcall("floor", list(.div(e$args[[1L]], e$args[[2L]]))))))
.so_alt_hooks$Mod <- function(e, x, logx, cdir) .so_eval_alt(.sp_mod_floor(e), x, logx, cdir)
.ns_hooks$Mod <- function(e, x, n, logx, cdir) .ns_eval(.sp_mod_floor(e), x, n, logx, cdir)

# ---------------------------------------------------------------------------
# Piecewise, Derivative, Subs, Integral, SingularityFunction

.so_alt_hooks$Piecewise <- function(e, x, logx, cdir) {
  for (p in e$args) {
    c <- p$args[[2L]]
    if (isTRUE(c$type == "BooleanTrue") || isTRUE(subs(c, x, S_Zero)$type == "BooleanTrue")) return(as_leading_term(p$args[[1L]], x))
  }
  NULL
}
.ns_hooks$Piecewise <- function(e, x, n, logx, cdir) {
  .piecewise(lapply(e$args, function(p) list(.ns_eval(p$args[[1L]], x, n, logx), p$args[[2L]])))
}
.sp_lseries_first <- function(e, x) {
  # the first nonzero term of Expr.lseries
  for (t in .so_lseries(e, x, 8L)) if (!is_zero_num(t)) return(t)
  S_Zero
}
.so_alt_hooks$Derivative <- function(e, x, logx, cdir) {
  lt <- .sp_lseries_first(e$args[[1L]], x)
  do.call(diff, c(list(lt), lapply(e$args[-1L], function(v) if (v$type == "Tuple") v$args[[1L]] else v)))
}
.ns_hooks$Derivative <- function(e, x, n, logx, cdir) {
  arg <- .ns_nseries(e$args[[1L]], x, n, logx)
  o <- .so_getO(arg)
  rest <- .so_removeO(arg)
  terms <- if (rest$type == "Add") rest$args else list(rest)
  rv <- lapply(terms, function(a) sym_rebuild(e, c(list(a), e$args[-1L])))
  if (!is.null(o)) rv <- c(rv, list(.div(o, x)))
  .add(rv)
}
.so_alt_hooks$Subs <- function(e, x, logx, cdir) {
  vars <- e$args[[2L]]$args; pts <- e$args[[3L]]$args
  i <- match(x$key, vapply(pts, .k, ""))
  if (!is.na(i)) return(as_leading_term(e$args[[1L]], vars[[i]]))
  if (x$key %in% vapply(vars, .k, "")) return(e)
  as_leading_term(e$args[[1L]], x)
}
.ns_hooks$Subs <- function(e, x, n, logx, cdir) {
  vars <- e$args[[2L]]$args; pts <- e$args[[3L]]$args
  i <- match(x$key, vapply(pts, .k, ""))
  other <- if (!is.na(i)) vars[[i]] else x
  arg <- .ns_nseries(e$args[[1L]], other, n, logx)
  o <- .so_getO(arg)
  rest <- .so_removeO(arg)
  terms <- if (rest$type == "Add") rest$args else list(rest)
  rv <- .add(lapply(terms, function(a) sym_rebuild(e, c(list(a), e$args[-1L]))))
  if (!is.null(o)) rv <- .add(list(rv, subs(o, other, x)))
  rv
}
.ns_hooks$Integral <- function(e, x, n, logx, cdir) {
  symb <- x
  lims <- e$args[-1L]
  for (l in lims) if (length(l$args) > 1L && x$key %in% vapply(l$args[-1L], .k, "")) { symb <- l$args[[1L]]; break }
  s <- .ns_nseries(e$args[[1L]], symb, n, logx)
  terms <- if (s$type == "Add") s$args else list(s)
  is_o <- vapply(terms, function(t) t$type == "Order", NA)
  order <- lapply(terms[is_o], function(o) subs(o, symb, x))
  .add(list(do.call(integrate, c(list(.add(terms[!is_o])), lapply(lims, .int_limit_arg))), .mul(list(.add(order), x))))
}
.int_limit_arg <- function(l) if (length(l$args) == 1L) l$args[[1L]] else lapply(l$args, identity)
.so_alt_hooks$Integral <- function(e, x, logx, cdir) {
  lt <- .sp_lseries_first(e$args[[1L]], x)
  do.call(integrate, c(list(lt), lapply(e$args[-1L], .int_limit_arg)))
}
.so_alt_hooks$SingularityFunction <- function(e, x, logx, cdir) {
  z <- e$args[[1L]]; a <- e$args[[2L]]; nn <- e$args[[3L]]
  shift <- subs(.sub(z, a), x, S_Zero)
  if (.is_true(Lt(nn, 0))) return(S_Zero)
  if (.sp_zero(nn) && .sp_zero(shift)) return(if (.sp_is(as_sym(cdir), S_NegativeOne)) S_Zero else S_One)
  if (isTRUE(ask_fact(shift, "positive"))) return(.pow(shift, nn))
  S_Zero
}
.ns_hooks$SingularityFunction <- function(e, x, n, logx, cdir) {
  z <- e$args[[1L]]; a <- e$args[[2L]]; nn <- e$args[[3L]]
  shift <- subs(.sub(z, a), x, S_Zero)
  if (.is_true(Lt(nn, 0))) return(S_Zero)
  if (.sp_zero(nn) && .sp_zero(shift)) return(if (.sp_is(as_sym(cdir), S_NegativeOne)) S_Zero else S_One)
  if (isTRUE(ask_fact(shift, "positive"))) return(.ns_eval(.pow(.sub(z, a), nn), x, nn, logx, cdir))
  S_Zero
}

# ---------------------------------------------------------------------------
# gamma family and factorials

.ns_hooks$gamma <- function(e, x, n, logx, cdir) {
  x0 <- .sp_lim0(e$args[[1L]], x)
  if (!(x0$type == "Integer" && num_value_dbl(x0) <= 0)) return(.ns_function(e, x, n, logx, 0))
  t <- .sub(e$args[[1L]], x0)
  .ns_eval(.div(fcall("gamma", list(.add(list(t, S_One)))), fcall("RisingFactorial", list(e$args[[1L]], .add(list(.neg(x0), S_One))))), x, n, logx)
}
.so_alt_hooks$gamma <- function(e, x, logx, cdir) {
  arg <- e$args[[1L]]
  x0 <- subs(arg, x, S_Zero)
  if (isTRUE(ask_fact(x0, "integer")) && isTRUE(ask_fact(x0, "nonpositive"))) {
    nn <- .neg(x0)
    res <- .div(.pow(S_NegativeOne, nn), fcall("gamma", list(.add(list(nn, S_One)))))
    return(.div(res, as_leading_term(.add(list(arg, nn)), x)))
  }
  if (!.sp_infinite(x0)) return(fcall("gamma", list(x0)))
  .pole_error("")
}
.so_alt_hooks$polygamma <- function(e, x, logx, cdir) {
  nz <- lapply(e$args, function(a) as_leading_term(a, x))
  o <- .so_new(nz[[2L]], list(x))
  if (is_zero_num(nz[[1L]]) && isTRUE(.so_contains(o, .inv(x)))) return(.mul(list(.so_getn(o), .sp_logx(x, logx))))
  fcall("polygamma", nz)
}
.so_alt_hooks$digamma <- function(e, x, logx, cdir) as_leading_term(fcall("polygamma", list(S_Zero, e$args[[1L]])), x)
.so_alt_hooks$trigamma <- function(e, x, logx, cdir) as_leading_term(fcall("polygamma", list(S_One, e$args[[1L]])), x)
.ns_hooks$loggamma <- function(e, x, n, logx, cdir) {
  x0 <- .sp_lim0(e$args[[1L]], x)
  if (.sp_zero(x0)) return(.ns_eval(fcall("log", list(fcall("gamma", list(e$args[[1L]])))), x, n, logx))
  .ns_function(e, x, n, logx, 0)
}
.so_alt_hooks$factorial <- function(e, x, logx, cdir) {
  arg <- as_leading_term(e$args[[1L]], x)
  arg0 <- subs(arg, x, S_Zero)
  if (.sp_zero(arg0)) return(S_One)
  if (!.sp_infinite(arg0)) return(fcall("factorial", list(arg)))
  .pole_error(paste("Cannot expand", .str(e), "around 0"))
}
.so_alt_hooks$binomial <- function(e, x, logx, cdir) .so_eval_alt(rewrite(e, "gamma"), x, logx, cdir)
.so_alt_hooks$LambertW <- function(e, x, logx, cdir) {
  if (length(e$args) != 1L) return(.so_alt_function(e, x, logx, cdir))
  arg <- e$args[[1L]]
  arg0 <- cancel(subs(arg, x, S_Zero))
  if (!.sp_zero(arg0)) return(fcall("LambertW", list(arg0)))
  as_leading_term(arg, x)
}
.ns_hooks$LambertW <- function(e, x, n, logx, cdir) {
  if (length(e$args) != 1L) return(.ns_function(e, x, n, logx, 0))
  arg <- .ns_nseries(e$args[[1L]], x, n, logx)
  lt <- as_leading_term(arg, x, logx)
  lte <- if (lt$type == "Pow") lt$args[[2L]] else S_One
  top <- .ns_ceiling(.div(as_sym(n), lte))
  s <- S_Zero
  if (.is_true(Ge(top, 1))) {
    kk <- seq_len(max(as.integer(num_value_dbl(top)) - 1L, 0L))
    s <- .add(lapply(kk, function(k) .div(.mul(list(.pow(S_NegativeOne, .int(k - 1L)), .pow(.int(k), .int(k - 2L)), .pow(arg, .int(k)))), .se_fact(k - 1L))))
    s <- expand_multinomial(s)
  }
  .add(list(s, .so_new(.pow(x, n), list(x))))
}

# ---------------------------------------------------------------------------
# asymptotic expansions of the gamma family

.ns_aseries_hooks$lowergamma <- function(e, n, args0, x, logx) {
  s <- e$args[[1L]]; z <- e$args[[2L]]
  if (.sp_is(args0[[1L]], S_Infinity) && !.has1(z, x)) {
    nn <- .sp_nn(n)
    coeff <- .mul(list(.pow(z, s), exp(.neg(z))))
    sm <- .add(lapply(seq_len(max(nn - 1L, 0L)) - 1L, function(k) .div(.pow(z, .int(k)), fcall("RisingFactorial", list(s, .int(k + 1L))))))
    return(.add(list(.mul(list(coeff, sm)), .so_new(.mul(list(.pow(z, s), .pow(s, .neg(as_sym(n))))))) ))
  }
  .pole_error(paste0("Asymptotic expansion of lowergamma around [", paste(vapply(args0, .str, ""), collapse = ", "), "] is not implemented."))
}
.ns_aseries_hooks$polygamma <- function(e, n, args0, x, logx) {
  N <- e$args[[1L]]; z <- e$args[[2L]]
  if (!.sp_is(args0[[2L]], S_Infinity) || !(N$type == "Integer" && num_value_dbl(N) >= 0))
    .pole_error(paste0("Asymptotic expansion of polygamma around [", paste(vapply(args0, .str, ""), collapse = ", "), "] is not implemented."))
  nn <- .sp_nn(n)
  if (is_zero_num(N)) {
    r <- .sub(fcall("log", list(z)), .div(S_One, .mul(list(S_Two, z))))
    if (nn < 2L) o <- .sp_o(.inv(z), x)
    else {
      m <- (nn + 1L) %/% 2L
      l <- lapply(seq_len(max(m - 1L, 0L)), function(k) .div(.se_F("bernoulli", 2L * k), .mul(list(.int(2L * k), .pow(z, .int(2L * k))))))
      r <- .sub(r, .add(l))
      o <- .sp_o(.pow(z, .int(-nn)), x)
    }
    return(.add(list(.ns_eval(r, x, n, logx), o)))
  }
  Ni <- as.integer(num_value_dbl(N))
  fac <- fcall("gamma", list(N))
  e0 <- .add(list(fac, .div(.mul(list(N, fac)), .mul(list(S_Two, z)))))
  m <- (nn + 1L) %/% 2L
  for (k in seq_len(max(m - 1L, 0L))) {
    fac <- .div(.mul(list(fac, .int(2L * k + Ni - 1L), .int(2L * k + Ni - 2L))), .int(2L * k * (2L * k - 1L)))
    e0 <- .add(list(e0, .div(.mul(list(.se_F("bernoulli", 2L * k), fac)), .pow(z, .int(2L * k)))))
  }
  o <- .sp_o(.pow(z, .int(-2L * m)), x)
  if (nn == 0L) o <- .sp_o(.inv(z), x) else if (nn == 1L) o <- .sp_o(.pow(z, .int(-2L)), x)
  r <- .add(list(.ns_eval(e0, z, n, logx), o))
  .ns_eval(.mul(list(S_NegativeOne, .pow(.neg(.inv(z)), N), r)), x, n, logx)
}
.ns_aseries_hooks$loggamma <- function(e, n, args0, x, logx) {
  if (!.sp_is(args0[[1L]], S_Infinity))
    .pole_error(paste0("Asymptotic expansion of loggamma around [", .str(args0[[1L]]), "] is not implemented."))
  z <- e$args[[1L]]
  nn <- .sp_nn(n)
  r <- .add(list(.mul(list(fcall("log", list(z)), .sub(z, S_Half))), .neg(z), .div(fcall("log", list(.mul(list(S_Two, S_Pi)))), S_Two)))
  l <- lapply(seq_len(max(nn - 1L, 0L)), function(k) .div(.se_F("bernoulli", 2L * k), .mul(list(.int(2L * k * (2L * k - 1L)), .pow(z, .int(2L * k - 1L))))))
  o <- if (nn == 0L) .sp_o(S_One, x) else .sp_o(.pow(z, .int(-nn)), x)
  .add(list(.ns_eval(.add(c(list(r), l)), x, n, logx), o))
}
.ns_aseries_hooks$digamma <- function(e, n, args0, x, logx) .ns_aseries_hooks$polygamma(fcall("polygamma", list(S_Zero, e$args[[1L]])), n, c(list(S_Zero), args0), x, logx)
.ns_aseries_hooks$trigamma <- function(e, n, args0, x, logx) .ns_aseries_hooks$polygamma(fcall("polygamma", list(S_One, e$args[[1L]])), n, c(list(S_One), args0), x, logx)

# ---------------------------------------------------------------------------
# error functions and exponential integrals

.ns_taylor_hooks$erf <- function(k, a, prev, e) {
  if (k < 0L || k %% 2L == 0L) return(S_Zero)
  kk <- (k - 1L) %/% 2L
  .div(.mul(list(S_Two, .pow(S_NegativeOne, .int(kk)), .pow(a, .int(k)))), .mul(list(.int(k), .se_fact(kk), .sqrt(S_Pi))))
}
.ns_taylor_hooks$erfc <- function(k, a, prev, e) {
  if (k == 0L) return(S_One)
  if (k < 0L || k %% 2L == 0L) return(S_Zero)
  .neg(.ns_taylor_hooks$erf(k, a, prev, e))
}
.ns_taylor_hooks$erfi <- function(k, a, prev, e) {
  if (k < 0L || k %% 2L == 0L) return(S_Zero)
  kk <- (k - 1L) %/% 2L
  .div(.mul(list(S_Two, .pow(a, .int(k)))), .mul(list(.int(k), .se_fact(kk), .sqrt(S_Pi))))
}
.so_alt_hooks$erf <- function(e, x, logx, cdir) {
  arg <- as_leading_term(e$args[[1L]], x, logx, cdir)
  arg0 <- subs(arg, x, S_Zero)
  if (.sp_is(arg0, S_ComplexInfinity)) arg0 <- limit(arg, x, S_Zero, if (.sp_is(as_sym(cdir), S_NegativeOne)) "-" else "+")
  if (.has1(arg, x) && .sp_zero(arg0)) return(.div(.mul(list(S_Two, arg)), .sqrt(S_Pi)))
  fcall("erf", list(arg0))
}
.so_alt_hooks$erfc <- function(e, x, logx, cdir) {
  arg <- as_leading_term(e$args[[1L]], x, logx, cdir)
  arg0 <- subs(arg, x, S_Zero)
  if (.sp_is(arg0, S_ComplexInfinity)) arg0 <- limit(arg, x, S_Zero, if (.sp_is(as_sym(cdir), S_NegativeOne)) "-" else "+")
  if (.sp_zero(arg0)) return(S_One)
  fcall("erfc", list(arg0))
}
.so_alt_hooks$erfi <- function(e, x, logx, cdir) {
  arg <- as_leading_term(e$args[[1L]], x, logx, cdir)
  arg0 <- subs(arg, x, S_Zero)
  if (.has1(arg, x) && .sp_zero(arg0)) return(.div(.mul(list(S_Two, arg)), .sqrt(S_Pi)))
  if (.sp_finite(arg0)) return(fcall("erfi", list(arg0)))
  fcall("erfi", list(arg))
}
.ns_aseries_hooks$erf <- function(e, n, args0, x, logx) {
  if (.sp_inf_pts(args0[[1L]])) {
    z <- e$args[[1L]]
    ex <- tryCatch(.so_leadterm(z, x)[[2L]], error = function(err) NULL)
    if (is.null(ex)) return(e)
    ex <- .neg(ex)
    if (isTRUE(ask_fact(ex, "positive"))) {
      newn <- as.integer(num_value_dbl(.ns_ceiling(.div(as_sym(n), ex))))
      s <- c(lapply(seq_len(newn) - 1L, function(k) .div(.mul(list(.pow(S_NegativeOne, .int(k)), .se_F("factorial2", 2L * k - 1L))), .mul(list(.pow(z, .int(2L * k + 1L)), .pow(S_Two, .int(k)))))),
             list(.sp_o(.inv(.pow(z, .int(newn))), x)))
      return(.sub(S_One, .mul(list(exp(.neg(.sq(z))), .inv(.sqrt(S_Pi)), .add(s)))))
    }
  }
  .pole_error(paste0("Asymptotic expansion of erf around [", .str(args0[[1L]]), "] is not implemented."))
}
.ns_aseries_hooks$erfc <- function(e, n, args0, x, logx) .sub(S_One, .ns_aseries_hooks$erf(fcall("erf", e$args), n, args0, x, logx))
.ns_aseries_hooks$erfi <- function(e, n, args0, x, logx) {
  if (.sp_is(args0[[1L]], S_Infinity)) {
    z <- e$args[[1L]]; nn <- .sp_nn(n)
    s <- c(lapply(seq_len(nn) - 1L, function(k) .div(.se_F("factorial2", 2L * k - 1L), .mul(list(.pow(S_Two, .int(k)), .pow(z, .int(2L * k + 1L)))))),
           list(.sp_o(.pow(z, .int(-nn)), x)))
    return(.add(list(.neg(S_I), .mul(list(exp(.sq(z)), .inv(.sqrt(S_Pi)), .add(s))))))
  }
  .pole_error(paste0("Asymptotic expansion of erfi around [", .str(args0[[1L]]), "] is not implemented."))
}
.sp_log_lt <- function(arg, x, logx, extra = S_Zero) {
  ce <- .so_as_coeff_exponent(arg, x)
  .add(list(fcall("log", list(ce[[1L]])), .mul(list(ce[[2L]], .sp_logx(x, logx))), S_EulerGamma, extra))
}
.so_alt_hooks$Ei <- function(e, x, logx, cdir) {
  x0 <- .sp_lim0(e$args[[1L]], x)
  arg <- as_leading_term(e$args[[1L]], x, NULL, cdir)
  cd <- .so_dir(arg, x, if (is_zero_num(as_sym(cdir))) S_One else as_sym(cdir))
  if (.sp_zero(x0)) return(.sp_log_lt(arg, x, logx, if (isTRUE(ask_fact(.se_re(cd), "negative"))) .neg(.mul(list(S_I, S_Pi))) else S_Zero))
  .so_alt_function(e, x, logx, cdir)
}
.ns_hooks$Ei <- function(e, x, n, logx, cdir) {
  x0 <- .sp_lim0(e$args[[1L]], x)
  if (.sp_zero(x0)) return(.ns_eval(.sp_ei_as_si(e$args[[1L]]), x, n, logx))
  .ns_function(e, x, n, logx, 0)
}
.ns_aseries_hooks$Ei <- function(e, n, args0, x, logx) {
  if (.sp_inf_pts(args0[[1L]])) {
    z <- e$args[[1L]]; nn <- .sp_nn(n)
    s <- c(lapply(seq_len(nn) - 1L, function(k) .div(.se_fact(k), .pow(z, .int(k)))), list(.sp_o(.pow(z, .int(-nn)), x)))
    return(.mul(list(exp(z), .inv(z), .add(s))))
  }
  .pole_error(paste0("Asymptotic expansion of Ei around [", .str(args0[[1L]]), "] is not implemented."))
}
.ns_hooks$expint <- function(e, x, n, logx, cdir) {
  nu <- e$args[[1L]]
  if (!.has1(nu, x)) {
    if (is_one_num(nu)) return(.ns_eval(.sp_expint_as_si(nu, e$args[[2L]]), x, n, logx))
    if (nu$type == "Integer" && num_value_dbl(nu) > 1) return(.ns_eval(.sp_expint_as_ei(nu, e$args[[2L]]), x, n, logx))
  }
  .ns_function(e, x, n, logx, 0)
}
.ns_aseries_hooks$expint <- function(e, n, args0, x, logx) {
  nu <- e$args[[1L]]
  if (.sp_is(args0[[2L]], S_Infinity)) {
    z <- e$args[[2L]]; nn <- .sp_nn(n)
    s <- c(lapply(seq_len(nn) - 1L, function(k) .div(.mul(list(.pow(S_NegativeOne, .int(k)), fcall("RisingFactorial", list(nu, .int(k))))), .pow(z, .int(k)))), list(.sp_o(.pow(z, .int(-nn)), x)))
    return(.mul(list(exp(.neg(z)), .inv(z), .add(s))))
  }
  .pole_error(paste0("Asymptotic expansion of expint around [", paste(vapply(args0, .str, ""), collapse = ", "), "] is not implemented."))
}
.ns_hooks$li <- function(e, x, n, logx, cdir) {
  z <- e$args[[1L]]; nn <- .sp_nn(n)
  s <- lapply(seq_len(max(nn - 1L, 0L)), function(k) .div(.pow(fcall("log", list(z)), .int(k)), .mul(list(.se_fact(k), .int(k)))))
  .add(c(list(S_EulerGamma, fcall("log", list(fcall("log", list(z))))), s))
}
.ns_hooks$Li <- function(e, x, n, logx, cdir) .ns_eval(.sp_li_as_li(e$args[[1L]]), x, n, logx)

# TrigonometricIntegral._eval_nseries (Si, Ci, Shi, Chi)
.sp_trigint_ns <- function(trig) { force(trig); function(e, x, n, logx, cdir) {
  if (!is_zero_num(subs(e$args[[1L]], x, S_Zero))) return(.ns_function(e, x, n, logx, 0))
  f0 <- !is_zero_num(fcall(trig, list(S_Zero)))
  base <- .ns_eval(fcall(trig, list(x)), x, n, logx)
  if (f0) base <- .sub(base, S_One)
  base <- .map_bottom_up(base, function(t) if (t$type == "Pow") .div(t, t$args[[2L]]) else t)
  if (f0) base <- .add(list(base, S_EulerGamma, fcall("log", list(x))))
  .ns_eval(subs(base, x, e$args[[1L]]), x, n, logx)
} }
.ns_hooks$Si <- .sp_trigint_ns("sin")
.ns_hooks$Ci <- .sp_trigint_ns("cos")
.ns_hooks$Shi <- .sp_trigint_ns("sinh")
.ns_hooks$Chi <- .sp_trigint_ns("cosh")
.sp_si_lt <- function(nm, use_logx) { force(nm); force(use_logx); function(e, x, logx, cdir) {
  arg <- if (nm == "Shi") as_leading_term(e$args[[1L]], x) else as_leading_term(e$args[[1L]], x, logx, cdir)
  arg0 <- subs(arg, x, S_Zero)
  if (is_nan(arg0)) arg0 <- limit(arg, x, S_Zero, .sp_redir(cdir))
  if (.sp_zero(arg0)) return(if (use_logx) .sp_log_lt(arg, x, logx) else arg)
  if (if (use_logx) .sp_finite(arg0) else !.sp_infinite(arg0)) return(fcall(nm, list(arg0)))
  e
} }
.so_alt_hooks$Si <- .sp_si_lt("Si", FALSE)
.so_alt_hooks$Shi <- .sp_si_lt("Shi", FALSE)
.so_alt_hooks$Ci <- .sp_si_lt("Ci", TRUE)
.so_alt_hooks$Chi <- .sp_si_lt("Chi", TRUE)
.sp_si_pq <- function(z, nn, x) {
  p <- c(lapply(seq_len(nn %/% 2L + 1L) - 1L, function(k) .div(.mul(list(.pow(S_NegativeOne, .int(k)), .se_fact(2L * k))), .pow(z, .int(2L * k + 1L)))), list(.sp_o(.pow(z, .int(-nn)), x)))
  q <- c(lapply(seq_len(nn %/% 2L) - 1L, function(k) .div(.mul(list(.pow(S_NegativeOne, .int(k)), .se_fact(2L * k + 1L))), .pow(z, .int(2L * (k + 1L))))), list(.sp_o(.pow(z, .int(-nn)), x)))
  list(.add(p), .add(q))
}
.ns_aseries_hooks$Si <- function(e, n, args0, x, logx) {
  if (.sp_is(args0[[1L]], S_Infinity)) {
    z <- e$args[[1L]]; pq <- .sp_si_pq(z, .sp_nn(n), x)
    return(.add(list(.div(S_Pi, S_Two), .neg(.mul(list(cos(z), pq[[1L]]))), .neg(.mul(list(sin(z), pq[[2L]]))))))
  }
  .pole_error(paste0("Asymptotic expansion of Si around [", .str(args0[[1L]]), "] is not implemented."))
}
.ns_aseries_hooks$Ci <- function(e, n, args0, x, logx) {
  if (.sp_inf_pts(args0[[1L]])) {
    z <- e$args[[1L]]; pq <- .sp_si_pq(z, .sp_nn(n), x)
    r <- .sub(.mul(list(sin(z), pq[[1L]])), .mul(list(cos(z), pq[[2L]])))
    if (.sp_is(args0[[1L]], S_NegativeInfinity)) r <- .add(list(r, .mul(list(S_I, S_Pi))))
    return(r)
  }
  .pole_error(paste0("Asymptotic expansion of Ci around [", .str(args0[[1L]]), "] is not implemented."))
}

# Fresnel integrals
.ns_taylor_hooks$fresnels <- function(k, a, prev, e) {
  if (k < 0L) return(S_Zero)
  .div(.mul(list(.pow(a, .int(3L)), .pow(.neg(.pow(a, .int(4L))), .int(k)), .pow(S_Two, .int(-2L * k - 1L)), .pow(S_Pi, .int(2L * k + 1L)))), .mul(list(.int(4L * k + 3L), .se_fact(2L * k + 1L))))
}
.ns_taylor_hooks$fresnelc <- function(k, a, prev, e) {
  if (k < 0L) return(S_Zero)
  .div(.mul(list(a, .pow(.neg(.pow(a, .int(4L))), .int(k)), .pow(S_Two, .int(-2L * k)), .pow(S_Pi, .int(2L * k)))), .mul(list(.int(4L * k + 1L), .se_fact(2L * k))))
}
.sp_fresnel_lt <- function(nm) { force(nm); function(e, x, logx, cdir) {
  arg <- as_leading_term(e$args[[1L]], x, logx, cdir)
  arg0 <- subs(arg, x, S_Zero)
  if (.sp_is(arg0, S_ComplexInfinity)) arg0 <- limit(arg, x, S_Zero, .sp_redir(cdir))
  if (.sp_zero(arg0)) return(if (nm == "fresnels") .div(.mul(list(S_Pi, .pow(arg, .int(3L)))), .int(6L)) else arg)
  if (.sp_inf_pts(arg0)) return(.add(list(if (.sp_is(arg0, S_Infinity)) S_Half else .neg(S_Half), .sp_o(x, x))))
  fcall(nm, list(arg0))
} }
.so_alt_hooks$fresnels <- .sp_fresnel_lt("fresnels")
.so_alt_hooks$fresnelc <- .sp_fresnel_lt("fresnelc")
.sp_fresnel_as <- function(nm) { force(nm); function(e, n, args0, x, logx) {
  pt <- args0[[1L]]
  if (!.sp_inf_pts(pt)) .pole_error(paste0("Asymptotic expansion of ", nm, " around [", .str(pt), "] is not implemented."))
  z <- e$args[[1L]]; nn <- .sp_nn(n)
  kp <- Filter(function(k) 4L * k + 3L < nn, seq_len(nn) - 1L)
  p <- lapply(kp, function(k) .div(.mul(list(.pow(S_NegativeOne, .int(k)), .se_fact(4L * k + 1L))), .mul(list(.pow(S_Two, .int(2L * k + 2L)), .pow(z, .int(4L * k + 3L)), .pow(S_Two, .int(2L * k)), .se_fact(2L * k)))))
  kq <- Filter(function(k) 4L * k + 1L < nn, seq_len(max(nn - 1L, 0L)))
  q <- c(list(.inv(.mul(list(S_Two, z)))), lapply(kq, function(k) .div(.mul(list(.pow(S_NegativeOne, .int(k)), .se_fact(4L * k - 1L))), .mul(list(.pow(S_Two, .int(2L * k + 1L)), .pow(z, .int(4L * k + 1L)), .pow(S_Two, .int(2L * k - 1L)), .se_fact(2L * k - 1L))))))
  c0 <- .sqrt(.div(S_Two, S_Pi))
  p <- lapply(p, function(t) .neg(.mul(list(c0, t))))
  q <- lapply(q, function(t) if (nm == "fresnels") .neg(.mul(list(c0, t))) else .mul(list(c0, t)))
  s <- if (.sp_is(pt, S_Infinity)) S_Half else .neg(S_Half)
  z2 <- .sq(z)
  body <- if (nm == "fresnels") .add(list(.mul(list(sin(z2), .add(p))), .mul(list(cos(z2), .add(q)))))
          else .add(list(.mul(list(cos(z2), .add(p))), .mul(list(sin(z2), .add(q)))))
  .add(list(s, subs(body, x, .mul(list(c0, x))), .sp_o(.pow(z, .int(-nn)), x)))
} }
.ns_aseries_hooks$fresnels <- .sp_fresnel_as("fresnels")
.ns_aseries_hooks$fresnelc <- .sp_fresnel_as("fresnelc")

# rewrites used above: Mod._eval_rewrite_as_floor, Ei._eval_rewrite_as_Si,
# expint._eval_rewrite_as_Si and _eval_rewrite_as_Ei, Li._eval_rewrite_as_li
.sp_ei_as_si <- function(z) {
  r <- .add(list(fcall("Shi", list(z)), fcall("Chi", list(z))))
  if (isTRUE(ask_fact(z, "negative"))) .sub(r, .mul(list(S_I, S_Pi))) else r
}
.sp_expint_as_si <- function(nu, z) if (!is_one_num(nu)) NULL else .sub(fcall("Shi", list(z)), fcall("Chi", list(z)))
.sp_expint_as_ei <- function(nu, z) {
  if (is_one_num(nu)) return(.sub(.neg(fcall("Ei", list(.mul(list(z, fcall("exp_polar", list(.neg(.mul(list(S_I, S_Pi)))))))))), .mul(list(S_I, S_Pi))))
  if (nu$type == "Integer" && num_value_dbl(nu) > 1) {
    v <- as.integer(num_value_dbl(nu))
    xx <- .neg(unpolarify(z))
    e1 <- .sp_expint_as_ei(S_One, z)
    return(.add(list(.mul(list(.div(.pow(xx, .int(v - 1L)), .se_fact(v - 1L)), e1)),
                     .mul(list(.div(exp(xx), .se_fact(v - 1L)), .add(lapply(seq_len(v - 1L) - 1L, function(k) .mul(list(.se_fact(v - k - 2L), .pow(xx, .int(k)))))))))))
  }
  NULL
}
.sp_li_as_li <- function(z) .sub(fcall("li", list(z)), fcall("li", list(S_Two)))
.rewrite_rules$floor$Mod <- function(a, b) .sub(a, .mul(list(b, fcall("floor", list(.div(a, b))))))
for (.t in c("Si", "Ci", "Chi", "Shi")) {
  .rewrite_rules[[.t]]$Ei <- .sp_ei_as_si
  .rewrite_rules[[.t]]$expint <- .sp_expint_as_si
}
rm(.t)
.rewrite_rules$Ei$expint <- .sp_expint_as_ei
.rewrite_rules$li$Li <- .sp_li_as_li

# _erfs and _eis: helpers that make erf, Ei and li tractable for the
# Gruntz algorithm (sympy/functions/special/error_functions.py)
def_function("_erfs",
  eval = function(args) if (isTRUE(ask_fact(args[[1L]], "zero"))) S_One else NULL,
  diff = function(args, i) .add(list(.div(.neg(S_Two), .sqrt(S_Pi)), .mul(list(S_Two, args[[1L]], fcall("_erfs", args))))))
def_function("_eis",
  diff = function(args, i) .sub(.inv(args[[1L]]), fcall("_eis", args)))
.ns_aseries_hooks$`_erfs` <- function(e, n, args0, x, logx) {
  pt <- args0[[1L]]
  t <- if (.sp_is(pt, S_Infinity)) S_Infinity else tryCatch(extract_multiplicatively(pt, S_I), error = function(err) NULL)
  if (!is.null(t) && .sp_is(t, S_Infinity)) {
    z <- e$args[[1L]]; nn <- .sp_nn(n)
    l <- lapply(seq_len(nn) - 1L, function(k) .mul(list(.inv(.sqrt(S_Pi)), .se_fact(2L * k), .pow(.int(-4L), .int(-k)), .inv(.se_fact(k)), .pow(.inv(z), .int(2L * k + 1L)))))
    o <- .sp_o(.pow(z, .int(-(2L * nn + 1L))), x)
    return(.add(list(.ns_eval(.add(l), x, n, logx), o)))
  }
  .pole_error(paste0("Asymptotic expansion of _erfs around [", .str(pt), "] is not implemented."))
}
.ns_aseries_hooks$`_eis` <- function(e, n, args0, x, logx) {
  if (!.sp_inf_pts(args0[[1L]])) .pole_error(paste0("Asymptotic expansion of _eis around [", .str(args0[[1L]]), "] is not implemented."))
  z <- e$args[[1L]]; nn <- .sp_nn(n)
  l <- lapply(seq_len(nn) - 1L, function(k) .mul(list(.se_fact(k), .pow(.inv(z), .int(k + 1L)))))
  o <- .sp_o(.pow(z, .int(-(nn + 1L))), x)
  .add(list(.ns_eval(.add(l), x, n, logx), o))
}
.so_alt_hooks$`_eis` <- function(e, x, logx, cdir) {
  x0 <- .sp_lim0(e$args[[1L]], x)
  if (.sp_zero(x0)) return(.so_eval_alt(.mul(list(exp(.neg(e$args[[1L]])), fcall("Ei", e$args))), x, logx, cdir))
  .so_alt_function(e, x, logx, cdir)
}
.ns_hooks$`_eis` <- function(e, x, n, logx, cdir) {
  x0 <- .sp_lim0(e$args[[1L]], x)
  if (.sp_zero(x0)) return(.ns_eval(.mul(list(exp(.neg(e$args[[1L]])), fcall("Ei", e$args))), x, n, logx))
  .ns_function(e, x, n, logx, 0)
}
# atan and acot at infinity
.ns_aseries_hooks$atan <- function(e, n, args0, x, logx) {
  if (.sp_inf_pts(args0[[1L]])) return(.ns_eval(.sub(.div(S_Pi, S_Two), fcall("atan", list(.inv(e$args[[1L]])))), x, n, logx))
  .pole_error(paste0("Asymptotic expansion of atan around [", .str(args0[[1L]]), "] is not implemented."))
}
.ns_aseries_hooks$acot <- function(e, n, args0, x, logx) {
  if (.sp_inf_pts(args0[[1L]])) return(.ns_eval(fcall("atan", list(.inv(e$args[[1L]]))), x, n, logx))
  .pole_error(paste0("Asymptotic expansion of acot around [", .str(args0[[1L]]), "] is not implemented."))
}

# hyper, meijerg, elliptic integrals, polylog, zeta and jtheta
.so_alt_hooks$hyper <- function(e, x, logx, cdir) {
  arg <- e$args[[3L]]
  x0 <- subs(arg, x, S_Zero)
  if (is_nan(x0)) x0 <- limit(arg, x, S_Zero, .sp_redir(cdir))
  if (is_zero_num(x0)) return(S_One)
  .so_alt_function(e, x, logx, cdir)
}
.ns_hooks$hyper <- function(e, x, n, logx, cdir) {
  arg <- e$args[[3L]]
  x0 <- .sp_lim0(arg, x)
  ap <- e$args[[1L]]$args; bq <- e$args[[2L]]$args
  if (!(identical(arg$key, x$key) && is_zero_num(x0))) return(hyperexpand(.ns_function(e, x, n, logx, 0)))
  terms <- lapply(seq_len(.sp_nn(n)) - 1L, function(i) {
    num <- .mul(lapply(ap, function(a) fcall("RisingFactorial", list(a, .int(i)))))
    den <- .mul(lapply(bq, function(b) fcall("RisingFactorial", list(b, .int(i)))))
    .div(.mul(list(num, .pow(arg, .int(i)))), .mul(list(den, .se_fact(i))))
  })
  .add(c(terms, list(.so_new(.pow(x, n), list(x)))))
}
.so_alt_hooks$meijerg <- function(e, x, logx, cdir) as_leading_term(hyperexpand(e), x, logx, cdir)
.ns_hooks$elliptic_k <- function(e, x, n, logx, cdir)
  hyperexpand(.ns_eval(.mul(list(S_Pi, S_Half, hyper(list(S_Half, S_Half), list(S_One), e$args[[1L]]))), x, n, logx))
.ns_hooks$elliptic_e <- function(e, x, n, logx, cdir) {
  if (length(e$args) == 1L)
    return(hyperexpand(.ns_eval(.mul(list(.div(S_Pi, S_Two), hyper(list(.neg(S_Half), S_Half), list(S_One), e$args[[1L]]))), x, n, logx)))
  .ns_function(e, x, n, logx, 0)
}
.ns_hooks$polylog <- function(e, x, n, logx, cdir) {
  nu <- e$args[[1L]]; z <- e$args[[2L]]
  z0 <- subs(z, x, S_Zero)
  if (is_nan(z0)) z0 <- limit(z, x, S_Zero, .sp_redir(cdir))
  if (.sp_zero(z0)) {
    ex <- tryCatch(.so_leadterm(z, x)[[2L]], ValueError = function(err) NULL, NotImplementedError = function(err) NULL)
    if (is.null(ex)) return(e)
    if (isTRUE(ask_fact(ex, "positive"))) {
      newn <- .sp_nn(.ns_ceiling(.div(as_sym(n), ex)))
      o <- .so_new(.pow(x, n), list(x))
      r <- .so_removeO(.ns_eval(z, x, n, logx, cdir))
      if (is_zero_num(r)) return(o)
      term <- r; s <- list(term)
      for (k in seq_len(max(newn - 2L, 0L)) + 1L) {
        term <- .mul(list(term, r))
        s[[length(s) + 1L]] <- .div(term, .pow(.int(k), nu))
      }
      return(.add(c(s, list(o))))
    }
  }
  .ns_function(e, x, n, logx, 0)
}
.so_alt_hooks$zeta <- function(e, x, logx, cdir) {
  s <- e$args[[1L]]
  a <- if (length(e$args) == 2L) e$args[[2L]] else S_One
  ce <- tryCatch(.so_leadterm(a, x), NotImplementedError = function(err) NULL)
  if (is.null(ce)) return(e)
  if (isTRUE(ask_fact(ce[[2L]], "negative")) && !isTRUE(ask_fact(s, "positive"))) .notimpl_error("")
  .so_alt_function(e, x, logx, cdir)
}
.ns_hooks$jtheta <- function(e, x, n, logx, cdir) {
  index <- e$args[[1L]]; z <- e$args[[2L]]; q <- e$args[[3L]]
  d <- if (length(e$args) > 3L) e$args[[4L]] else S_Zero
  if (!.has1(q, x) || !isTRUE(ask_fact(subs(q, x, S_Zero), "zero"))) return(.ns_function(e, x, n, logx, 0))
  if (!(index$type == "Integer" && num_value_dbl(index) %in% 1:4) || d$type != "Integer" || num_value_dbl(d) < 0)
    .notimpl_error("jtheta series requires an explicit index and derivative order")
  qo <- .so_as_coeff_exponent(as_leading_term(q, x, logx, cdir), x)[[2L]]
  if (!.expr_is_number(qo) || !isTRUE(ask_fact(qo, "positive"))) .notimpl_error("jtheta series requires a nome with explicit positive order")
  idx <- as.integer(num_value_dbl(index)); nn <- num_value_dbl(as_sym(n)); qd <- num_value_dbl(qo)
  terms <- list()
  if (idx %in% 3:4) {
    if (is_zero_num(d)) terms[[1L]] <- S_One
    k <- 1L
    while (k^2 * qd < nn) {
      cf <- .mul(list(S_Two, .pow(.int(2L * k), d), cos(.add(list(.mul(list(.int(2L * k), z)), .div(.mul(list(d, S_Pi)), S_Two))))))
      if (idx == 4L) cf <- .mul(list(cf, .pow(S_NegativeOne, .int(k))))
      terms[[length(terms) + 1L]] <- .mul(list(cf, .pow(q, .int(k^2))))
      k <- k + 1L
    }
  } else {
    k <- 0L
    while ((k + 0.5)^2 * qd < nn) {
      h <- .int(2L * k + 1L)
      ph <- .add(list(.mul(list(h, z)), .div(.mul(list(d, S_Pi)), S_Two)))
      cf <- if (idx == 1L) .mul(list(S_Two, .pow(S_NegativeOne, .int(k)), .pow(h, d), sin(ph))) else .mul(list(S_Two, .pow(h, d), cos(ph)))
      terms[[length(terms) + 1L]] <- .mul(list(cf, .pow(q, .pow(.add(list(.int(k), S_Half)), S_Two))))
      k <- k + 1L
    }
  }
  .add(list(.ns_eval(.add(terms), x, n, logx, cdir), .so_new(.pow(x, n), list(x))))
}
