# Series expansion at x = 0, after the _eval_nseries methods of
# sympy/core/add.py, mul.py, power.py, function.py and
# functions/elementary/exponential.py, with Expr.nseries, Expr.leadterm and
# Expr.as_coeff_exponent.


# Expr.nseries for x0 = 0 from the right
.ns_nseries <- function(e, x, n, logx = NULL, cdir = 0) {
  if (!(x$key %in% names(.free_sym_list(e)))) return(e)
  .ns_eval(e, x, n, logx, cdir)
}

# _eval_nseries dispatch
.ns_eval <- function(e, x, n, logx = NULL, cdir = 0) {
  n <- as_sym(n)
  switch(e$type,
    Integer = , Rational = , Float = , Constant = , Symbol = , Dummy = , Wild = e,
    Order = e,
    Add = .add(lapply(.sympy_args(e), function(t) .ns_nseries(t, x, n, logx, cdir))),
    Mul = .ns_mul(e, x, n, logx, cdir),
    Pow = .ns_pow(e, x, n, logx, cdir),
    Function = {
      h <- .ns_hooks[[e$name]]
      if (!is.null(h)) h(e, x, n, logx, cdir) else .ns_function(e, x, n, logx, cdir)
    },
    {
      h <- .ns_hooks[[e$type]]
      if (!is.null(h)) h(e, x, n, logx, cdir) else .ns_no_nseries(e)
    })
}

# Expr._eval_nseries (sympy/core/expr.py:3495), whose message names self.func
.ns_no_nseries <- function(e) {
  cl <- switch(e$type, AccumulationBounds = "sympy.calculus.accumulationbounds.AccumulationBounds", paste0("sympy.", e$type))
  .notimpl_error(filldedent(sprintf("
                     The _eval_nseries method should be added to
                     <class '%s'> to give terms up to O(x**n) at x=0
                     from the positive direction so it is available when
                     nseries calls it.", cl)))
}

.ns_ceiling <- function(v) { v <- as_sym(v); if (is_number_type(v)) .int(ceiling(num_value_dbl(v))) else fcall("ceiling", list(v)) }
.ns_lt <- function(a, b) .is_true(Lt(a, b))

# Expr.as_coeff_exponent
.so_as_coeff_exponent <- function(e, x) {
  s <- collect(e, x)
  cm <- .as_coeff_mul_dep(s, x)
  if (length(cm[[2L]]) == 1L) {
    be <- .as_base_exp(cm[[2L]][[1L]])
    if (identical(be[[1L]]$key, x$key)) return(list(cm[[1L]], be[[2L]]))
  }
  list(s, S_Zero)
}

# Expr.leadterm
.so_leadterm <- function(e, x, logx = NULL, cdir = 0) {
  l <- as_leading_term(e, x, logx, cdir)
  d <- Dummy("logx")
  lx <- fcall("log", list(x))
  if (.has1(l, lx)) l <- subs(l, lx, d)
  ce <- .so_as_coeff_exponent(l, x)
  if (x$key %in% names(.free_sym_list(ce[[1L]])))
    stop(structure(class = c("ValueError", "error", "condition"), list(message = paste0("cannot compute leadterm(", .str(e), ", ", .str(x), "). The coefficient should have been free of ", .str(x), " but got ", .str(ce[[1L]])), call = NULL)))
  list(subs(ce[[1L]], d, lx), ce[[2L]])
}

# Mul._eval_nseries
.ns_mul <- function(e, x, n, logx, cdir) {
  margs <- .sympy_args(e)
  coeff_exp <- function(term) {
    lt <- .so_as_coeff_exponent(term, x)
    if (.has1(lt[[1L]], x)) lt <- tryCatch(.so_leadterm(term, x), ValueError = function(err) list(term, S_Zero))
    lt
  }
  ords <- list()
  facs <- tryCatch({
    for (t in margs) {
      ce <- .so_leadterm(t, x)
      if (!.has1(ce[[1L]], x)) ords[[length(ords) + 1L]] <- list(t, ce[[2L]]) else stop(structure(class = c("ValueError", "error", "condition"), list(message = "", call = NULL)))
    }
    n0 <- .add(lapply(Filter(function(o) is_number_type(o[[2L]]), ords), `[[`, 2L))
    fl <- list()
    for (o in ords) {
      m <- o[[2L]]
      n1 <- .ns_ceiling(.add(list(n, .neg(n0), if (is_number_type(m)) m else S_Zero)))
      s <- .ns_nseries(o[[1L]], x, n1, logx, cdir)
      ns <- .so_getn(s)
      if (!is.null(ns) && .ns_lt(ns, n1)) n <- .sub(n, .sub(n1, ns))
      fl[[length(fl) + 1L]] <- s
    }
    fl
  }, ValueError = function(err) NULL, NotImplementedError = function(err) NULL, PoleError = function(err) NULL)
  if (is.null(facs)) {
    n0 <- .add(lapply(Filter(function(o) is_number_type(o[[2L]]), ords), `[[`, 2L))
    if (isTRUE(ask_fact(n0, "nonnegative"))) n0 <- S_Zero
    facs <- lapply(margs, function(t) .ns_nseries(t, x, .ns_ceiling(.sub(n, n0)), logx, cdir))
    res <- powsimp(expand(.mul(facs)), combine = "exp", deep = TRUE)
    if (.has_type(res, "Order")) res <- .add(list(res, .so_new(.pow(x, n), list(x))))
    return(res)
  }
  res <- S_Zero
  ords2 <- lapply(facs, .mi_make_args_add)
  idx <- rep(1L, length(ords2))
  repeat {
    fac <- lapply(seq_along(ords2), function(i) ords2[[i]][[idx[i]]])
    ce <- lapply(fac, coeff_exp)
    power <- .add(lapply(ce, `[[`, 2L))
    if (isTRUE(ask_fact(.sub(power, n), "negative"))) res <- .add(list(res, .mul(c(lapply(ce, `[[`, 1L), list(.pow(x, power))))))
    k <- length(idx)
    while (k >= 1L) { idx[k] <- idx[k] + 1L; if (idx[k] <= length(ords2[[k]])) break; idx[k] <- 1L; k <- k - 1L }
    if (k < 1L) break
  }
  max_degree <- function(z) {
    if (identical(z$key, x$key)) return(S_One)
    if (!length(z$args)) return(S_Zero)
    if (z$type == "Add") return(Reduce(function(a, b) if (.is_true(Ge(a, b))) a else b, lapply(z$args, max_degree)))
    if (z$type == "Mul") return(.add(lapply(z$args, max_degree)))
    if (z$type == "Pow") return(.mul(list(max_degree(z$args[[1L]]), z$args[[2L]])))
    S_Zero
  }
  if (.mi_is_poly(e, x)) {
    ok <- tryCatch({
      if (.is_true(Ge(max_degree(e), n)) || degree(e, x) != degree(res, x)) res <- .add(list(res, .so_new(.pow(x, n), list(x))))
      TRUE
    }, error = function(err) FALSE)
    if (ok) return(res)
  }
  if (!identical(res$key, e$key)) {
    if (is_zero_num(subs(.sub(e, res), x, S_Zero)) && .is_true(Gt(n, 0))) {
      lt <- .so_eval_alt(e, x, logx, cdir)
      if (is_zero_num(lt)) return(res)
    }
    res <- .add(list(res, .so_new(.pow(x, n), list(x))))
  }
  res
}

# Pow._eval_nseries
.ns_pow <- function(e, x, n, logx, cdir) {
  e <- trigsimp(powdenest(e, force = TRUE))
  if (e$type == "Function" && e$name == "exp") return(.ns_eval(e, x, n, logx, cdir))
  if (e$type != "Pow") return(.ns_nseries(e, x, n, logx, cdir))
  b <- e$args[[1L]]; ex <- e$args[[2L]]
  if (.it_illegal(ex)) .pole_error("")
  if (.has1(ex, x)) return(.ns_eval(fcall("exp", list(.mul(list(ex, fcall("log", list(b)))))), x, n, logx, cdir))
  if (!is.null(logx) && .has_function(b, "log")) {
    c <- Wild("c", exclude = list(x)); exw <- Wild("ex", exclude = list(x))
    b <- .map_bottom_up(b, function(z) {
      if (z$type == "Function" && z$name == "log") {
        m <- tryCatch(sym_match(z$args[[1L]], c * x^exw), error = function(err) NULL)
        if (!is.null(m) && !is.null(m$c)) return(.add(list(fcall("log", list(m$c)), .mul(list(m$ex %||% S_One, logx)))))
      }
      z
    })
    e <- .pow(b, ex)
  }
  b <- .so_removeO(b)
  m <- tryCatch({
    if ((.has_function(b, "polygamma") || .has1(b, S_EulerGamma)) && !is.null(logx)) stop(structure(class = c("ValueError", "error", "condition"), list(message = "", call = NULL)))
    .so_leadterm(b, x)[[2L]]
  }, ValueError = function(err) NULL, NotImplementedError = function(err) NULL, PoleError = function(err) NULL)
  if (is.null(m)) {
    b <- .so_removeO(.ns_eval(b, x, .int(max(2, num_value_dbl(n))), logx, cdir))
    if (.has1(b, S_NaN) || .has1(b, S_ComplexInfinity)) .notimpl_error("")
    m <- .so_leadterm(b, x)[[2L]]
  }
  if (.has_function(ex, "log")) ex <- cancel(logcombine(ex))
  if (!(is_zero_num(m) || (.expr_is_number(ex) && isTRUE(ask_fact(ex, "real"))))) {
    if (identical(e$key, .so_eval_alt(e, x, logx, cdir)$key)) {
      el <- fcall("exp", list(.mul(list(ex, fcall("log", list(b))))))
      res <- .ns_eval(el, x, n, logx, cdir)
      if (identical(res$key, el$key)) return(e)
      return(res)
    }
  }
  f <- as_leading_term(b, x, logx)
  g <- .div(cancel(.sub(expand_mul(expand_multinomial(b)), f)), f)
  if (!is_number_type(m)) .notimpl_error("")
  maxpow <- .sub(n, .mul(list(m, ex)))
  if (.has_type(maxpow, "Symbol")) maxpow <- n
  if (isTRUE(ask_fact(maxpow, "negative"))) return(.so_new(.pow(x, .mul(list(m, ex))), list(x)))
  if (isTRUE(ask_fact(g, "zero"))) {
    r <- .pow(f, ex)
    if (!identical(r$key, e$key)) r <- .add(list(r, .so_new(.pow(x, n), list(x))))
    return(r)
  }
  coeff_exp <- function(term) {
    coeff <- S_One; expo <- S_Zero
    for (factor in .mi_make_args_mul(term)) {
      if (.has1(factor, x)) {
        be <- .as_base_exp(factor)
        if (!identical(be[[1L]]$key, x$key)) return(tryCatch(.so_leadterm(term, x), ValueError = function(err) list(term, S_Zero)))
        expo <- be[[2L]]
      } else coeff <- .mul(list(coeff, factor))
    }
    list(coeff, expo)
  }
  dmul <- function(d1, d2) {
    res <- list()
    for (k1 in names(d1)) for (k2 in names(d2)) {
      exk <- .add(list(d1[[k1]][[1L]], d2[[k2]][[1L]]))
      if (.ns_lt(exk, maxpow)) {
        kk <- exk$key
        res[[kk]] <- list(exk, .add(list(if (is.null(res[[kk]])) S_Zero else res[[kk]][[2L]], .mul(list(d1[[k1]][[2L]], d2[[k2]][[2L]])))))
      }
    }
    res
  }
  cd <- tryCatch(.so_leadterm(g, x, logx), ValueError = function(err) NULL, NotImplementedError = function(err) NULL)
  if (is.null(cd)) {
    if (is_zero_num(limit(.div(g, .pow(x, maxpow)), x, S_Zero))) return(.add(list(.pow(f, ex), .mul(list(ex, .pow(f, ex), g)))))
    .notimpl_error("")
  }
  c <- cd[[1L]]; d <- cd[[2L]]
  if (c$type == "Float" && is_zero_num(d)) {
    g <- .map_bottom_up(g, function(z) if (z$type == "Float") nsimplify(z, rational = TRUE) else z)
    d <- .so_leadterm(g, x, logx)[[2L]]
  }
  if (!isTRUE(ask_fact(d, "positive"))) {
    g <- simplify(g)
    if (isTRUE(ask_fact(g, "zero"))) return(.pow(f, ex))
    d <- .so_leadterm(g, x, logx)[[2L]]
    if (!isTRUE(ask_fact(d, "positive"))) {
      g <- expand(.div(.sub(b, f), f))
      d <- .so_leadterm(g, x, logx)[[2L]]
      if (!isTRUE(ask_fact(d, "positive"))) .notimpl_error("")
    }
  }
  gpoly <- .so_removeO(.ns_eval(g, x, .ns_ceiling(maxpow), logx, cdir))
  gterms <- list()
  for (term in .mi_make_args_add(gpoly)) {
    ce <- coeff_exp(term)
    kk <- ce[[2L]]$key
    gterms[[kk]] <- list(ce[[2L]], .add(list(if (is.null(gterms[[kk]])) S_Zero else gterms[[kk]][[2L]], ce[[1L]])))
  }
  k <- S_One
  terms <- list()
  terms[[S_Zero$key]] <- list(S_Zero, S_One)
  tk <- gterms
  while (isTRUE(ask_fact(.sub(.mul(list(k, d)), maxpow), "negative"))) {
    cf <- .div(fcall("FallingFactorial", list(ex, k)), fcall("factorial", list(k)))
    for (kk in names(tk)) terms[[kk]] <- list(tk[[kk]][[1L]], .add(list(if (is.null(terms[[kk]])) S_Zero else terms[[kk]][[2L]], .mul(list(cf, tk[[kk]][[2L]])))))
    tk <- dmul(tk, gterms)
    k <- .add(list(k, S_One))
  }
  if (!isTRUE(ask_fact(ex, "integer")) && is_zero_num(m) && isTRUE(ask_fact(f, "negative"))) {
    ndir <- .so_dir(.sub(b, f), x, cdir)
    if (isTRUE(ask_fact(fcall("im", list(ndir)), "negative"))) ie <- coeff_exp(.mul(list(.pow(f, ex), .pow(S_NegativeOne, .mul(list(.int(-2), ex))))))
    else if (isTRUE(ask_fact(fcall("im", list(ndir)), "zero"))) ie <- coeff_exp(as_leading_term(fcall("exp", list(.mul(list(ex, fcall("log", list(b)))))), x, logx, cdir))
    else ie <- coeff_exp(.pow(f, ex))
  } else ie <- coeff_exp(.pow(f, ex))
  inco <- ie[[1L]]; inex <- ie[[2L]]
  res <- .add(lapply(terms, function(tt) .mul(list(tt[[2L]], inco, .pow(x, .add(list(tt[[1L]], inex)))))))
  if (!(isTRUE(ask_fact(ex, "integer")) && isTRUE(ask_fact(ex, "positive")) && isTRUE(ask_fact(.sub(.mul(list(ex, d)), n), "nonpositive")) &&
        identical(res$key, expand_mul(expand_multinomial(e))$key))) {
    r2 <- tryCatch(.add(list(res, .so_new(.pow(x, n), list(x)))), NotImplementedError = function(err) NULL)
    if (is.null(r2)) return(.ns_eval(fcall("exp", list(.mul(list(ex, fcall("log", list(b)))))), x, n, logx, cdir))
    res <- r2
  }
  res
}

# Function._eval_nseries
.ns_function <- function(e, x, n, logx, cdir) {
  args <- e$args
  args0 <- lapply(args, function(t) limit(t, x, S_Zero))
  if (any(vapply(args0, function(t) identical(ask_fact(t, "finite"), FALSE), NA))) {
    a <- lapply(args, function(t) as_leading_term(t, x, logx))
    a0 <- lapply(a, function(t) limit(t, x, S_Zero))
    if (any(vapply(a0, .it_illegal, NA))) return(.ns_aseries(e, n, args0, x, logx))
    a <- lapply(args, function(t) .ns_eval(t, x, n, logx))
    z <- lapply(seq_along(a), function(i) .sub(a[[i]], a0[[i]]))
    p <- lapply(z, function(zz) Dummy())
    q <- list(); v <- NULL; zi <- NULL
    for (i in seq_along(a0)) {
      if (.has1(z[[i]], x)) {
        if (!is.null(v)) .notimpl_error("")
        q[[i]] <- .add(list(a0[[i]], p[[i]])); v <- p[[i]]; zi <- z[[i]]
      } else q[[i]] <- a0[[i]]
    }
    e1 <- fcall(e$name, q)
    if (is.null(v)) return(e1)
    s <- .ns_eval(e1, v, n, logx)
    o <- .so_getO(s); s <- .so_removeO(s)
    return(.add(list(expand(subs(s, v, zi)), .so_new(subs(o$args[[1L]], v, zi), list(x)))))
  }
  nargs <- .FUNS[[e$name]]$nargs
  if (is.null(nargs) || (identical(as.integer(nargs), 1L) && !is_zero_num(args0[[1L]])) || any(nargs > 1L)) {
    e1 <- expand(e)
    if (identical(e1$key, e$key)) {
      ee <- e
      if (length(ee$args) == 1L) ee <- fcall(ee$name, list(cancel(ee$args[[1L]])))
      term <- subs(ee, x, S_Zero)
      if (identical(ask_fact(term, "finite"), FALSE) || is_nan(term)) .pole_error(paste("Cannot expand", .str(e), "around 0"))
      series <- term
      fact <- S_One
      .x <- Dummy("xi")
      ee <- subs(ee, x, .x)
      nn <- as.integer(num_value_dbl(n))
      if (nn > 1L) for (i in seq_len(nn - 1L)) {
        fact <- .mul(list(fact, .int(i)))
        ee <- .diff1(ee, .x)
        sb <- subs(ee, .x, S_Zero)
        if (is_nan(sb)) sb <- limit(ee, .x, S_Zero)
        if (identical(ask_fact(sb, "finite"), FALSE)) .pole_error(paste("Cannot expand", .str(e), "around 0"))
        series <- .add(list(series, expand(.div(.mul(list(sb, .pow(x, .int(i)))), fact))))
      }
      return(.add(list(series, .so_new(.pow(x, n), list(x)))))
    }
    return(.ns_nseries(e1, x, n, logx))
  }
  arg <- args[[1L]]
  l <- list(); g <- NULL
  nterms <- as.integer(num_value_dbl(n)) + 2L
  cf <- .so_getn(.so_new(as_leading_term(arg, x), list(x)))
  if (!is_zero_num(cf)) nterms <- as.integer(ceiling(num_value_dbl(n) / num_value_dbl(cf)))
  tt <- .ns_taylor_term_fn(e$name)
  for (i in seq_len(max(nterms, 0L)) - 1L) {
    g <- tt(i, arg, g, e)
    g <- .ns_nseries(g, x, n, logx)
    l[[length(l) + 1L]] <- g
  }
  .add(c(l, list(.so_new(.pow(x, n), list(x)))))
}

# taylor_term of a function class (previous term passed as `prev`)
.ns_taylor_term_fn <- function(name) {
  h <- .ns_taylor_hooks[[name]]
  if (!is.null(h)) return(h)
  # Expr.taylor_term: the n-th derivative at 0
  function(k, a, prev, e) {
    .x <- Dummy("x")
    fe <- fcall(name, list(.x))
    d <- if (k == 0L) fe else diff(fe, .x, k)
    .div(.mul(list(subs(d, .x, S_Zero), .pow(a, .int(k)))), fcall("factorial", list(.int(k))))
  }
}

# Function._taylor
.ns_taylor <- function(e, x, n) {
  l <- list(); g <- NULL
  tt <- .ns_taylor_term_fn(e$name)
  nn <- as.integer(num_value_dbl(n))
  for (i in seq_len(max(nn, 0L)) - 1L) {
    g <- tt(i, e$args[[1L]], g, e)
    g <- .ns_nseries(g, x, n)
    l[[length(l) + 1L]] <- .so_removeO(g)
  }
  .add(l)
}

# exp
.ns_taylor_hooks$exp <- function(k, a, prev, e) {
  if (k < 0L) return(S_Zero)
  if (k == 0L) return(S_One)
  if (!is.null(prev)) return(.div(.mul(list(prev, a)), .int(k)))
  .div(.pow(a, .int(k)), fcall("factorial", list(.int(k))))
}
.ns_hooks$exp <- function(e, x, n, logx, cdir) {
  arg <- e$args[[1L]]
  arg_series <- .ns_eval(arg, x, n, logx)
  if (.so_is_order(arg_series)) return(.add(list(S_One, arg_series)))
  arg0 <- limit(.so_removeO(arg_series), x, S_Zero)
  if (is_const(arg0, "-oo")) return(.so_new(.pow(x, n), list(x)))
  if (is_const(arg0, "oo")) return(e)
  if (isTRUE(ask_fact(arg0, "infinite"))) .pole_error(paste("Cannot expand", .str(e), "around 0"))
  if (any(vapply(arg0$args, function(a) a$type == "Function" && a$name == "sign", NA))) return(e)
  t <- Dummy("t")
  nterms <- n
  cf <- tryCatch(.so_getn(.so_new(as_leading_term(arg, x, logx), list(x))), NotImplementedError = function(err) S_Zero, PoleError = function(err) S_Zero)
  if (!is.null(cf) && !is_zero_num(cf) && isTRUE(ask_fact(cf, "positive"))) nterms <- .ns_ceiling(.div(n, cf))
  exp_series <- .ns_taylor(fcall("exp", list(t)), t, nterms)
  r <- .mul(list(fcall("exp", list(arg0)), subs(exp_series, t, .sub(arg_series, arg0))))
  rr <- if (!is.null(logx)) subs(r, logx, fcall("log", list(x))) else r
  if (identical(rr$key, e$key)) return(r)
  if (!is.null(cf) && !is_zero_num(cf) && isTRUE(ask_fact(.sub(cf, S_One), "positive")))
    r <- .add(list(r, .div(.so_new(.pow(.sub(arg_series, arg0), n), list(x)), .pow(x, .mul(list(.sub(cf, S_One), n))))))
  else r <- .add(list(r, .so_new(.pow(.sub(arg_series, arg0), n), list(x))))
  r <- expand(r)
  r <- powsimp(r, deep = TRUE, combine = "exp")
  .map_bottom_up(r, function(z) {
    if (z$type == "Pow" && identical(z$args[[1L]]$key, S_NegativeOne$key) && z$args[[2L]]$type == "Rational" &&
        .mi_q(z$args[[2L]]) %in% c(3L, 4L, 6L)) return(expand_complex(z))
    z
  })
}
.so_alt_hooks$exp <- function(e, x, logx, cdir) {
  arg <- as_leading_term(cancel(e$args[[1L]]), x, logx)
  arg0 <- subs(arg, x, S_Zero)
  if (is_nan(arg)) return(S_NaN)
  if (arg0$type == "AccumulationBounds") {
    if (isTRUE(ask_fact(fcall("re", list(as_sym(cdir))), "negative"))) return(fcall("exp", list(.neg(arg0))))
    return(fcall("exp", list(arg0)))
  }
  if (is_nan(arg0)) arg0 <- limit(arg, x, S_Zero)
  if (identical(ask_fact(arg0, "infinite"), FALSE)) return(fcall("exp", list(arg0)))
  .pole_error(paste("Cannot expand", .str(e), "around 0"))
}

# log
.ns_taylor_hooks$log <- function(k, a, prev, e) {
  if (k < 0L) return(S_Zero)
  if (k == 0L) return(a)
  if (!is.null(prev)) return(powsimp(.div(.mul(list(.int(-k), prev, a)), .int(k + 1L)), deep = TRUE, combine = "exp"))
  .div(.mul(list(.int(1L - 2L * (k %% 2L)), .pow(a, .int(k + 1L)))), .int(k + 1L))
}
.ns_hooks$log <- function(e, x, n, logx, cdir) {
  if (length(e$args) != 1L) return(.ns_function(e, x, n, logx, cdir))
  arg <- e$args[[1L]]
  if (identical(arg$key, x$key)) return(if (is.null(logx)) fcall("log", list(x)) else logx)
  t <- Dummy("t", positive = TRUE)
  cdir <- as_sym(cdir)
  if (is_zero_num(cdir)) cdir <- S_One
  z <- subs(arg, x, .mul(list(cdir, t)))
  kw <- Wild("k"); lw <- Wild("l")
  r <- tryCatch(sym_match(z, kw * t^lw), error = function(err) NULL)
  if (!is.null(r) && length(r)) {
    k <- r$k %||% S_One; l <- r$l %||% S_One
    if (!is_zero_num(l) && !.has1(l, t) && !.has1(k, t)) {
      rr <- if (is.null(logx)) .mul(list(l, fcall("log", list(x)))) else .mul(list(l, logx))
      return(.add(list(rr, fcall("log", list(k)), .neg(.mul(list(l, fcall("log", list(cdir))))))))
    }
  }
  coeff_exp <- function(term, v) {
    coeff <- S_One; expo <- S_Zero
    for (factor in .mi_make_args_mul(term)) {
      if (.has1(factor, v)) {
        be <- .as_base_exp(factor)
        if (!identical(be[[1L]]$key, v$key)) return(tryCatch(.so_leadterm(term, v), ValueError = function(err) list(term, S_Zero)))
        expo <- be[[2L]]
      } else coeff <- .mul(list(coeff, factor))
    }
    list(coeff, expo)
  }
  ab <- tryCatch(.so_leadterm(z, t, logx, 1), ValueError = function(err) NULL, NotImplementedError = function(err) NULL, PoleError = function(err) NULL)
  if (is.null(ab)) {
    s <- .ns_eval(z, t, n, logx, 1)
    while (.so_is_order(s)) { n <- .add(list(n, S_One)); s <- .ns_eval(z, t, n, logx, 1) }
    ab <- tryCatch(.so_leadterm(.so_removeO(s), t, NULL, 1), ValueError = function(err) list(as_leading_term(.so_removeO(s), t, NULL, 1), S_Zero))
  }
  a <- ab[[1L]]; b <- ab[[2L]]
  p <- .ns_eval(cancel(.sub(.div(z, .mul(list(a, .pow(t, b)))), S_One)), t, n, logx, 1)
  if (.has_function(p, "exp")) p <- logcombine(p)
  if (.so_is_order(p)) n <- .so_getn(p)
  d <- coeff_exp(p, t)[[2L]]
  lx <- if (is.null(logx)) fcall("log", list(x)) else logx
  if (!isTRUE(ask_fact(d, "positive"))) {
    res <- .add(list(fcall("log", list(a)), .neg(.mul(list(b, fcall("log", list(cdir))))), .mul(list(b, lx))))
    lf <- function(z) expand(z, deep = TRUE, log = TRUE, mul = FALSE, power_exp = FALSE, power_base = FALSE, multinomial = FALSE, basic = FALSE, force = TRUE, factor = FALSE)
    ex <- lf(e)
    .res <- if (!could_extract_minus_sign(a) && could_extract_minus_sign(lx)) lf(subs(res, .neg(lx), .neg(fcall("log", list(x))))) else lf(subs(res, lx, fcall("log", list(x))))
    if (identical(.res$key, ex$key)) return(res)
    return(.add(list(res, .so_new(.pow(x, n), list(x)))))
  }
  pterms <- list()
  for (term in .mi_make_args_add(.so_removeO(p))) {
    ce <- coeff_exp(term, t)
    kk <- ce[[2L]]$key
    pterms[[kk]] <- list(ce[[2L]], .add(list(if (is.null(pterms[[kk]])) S_Zero else pterms[[kk]][[2L]], ce[[1L]])))
  }
  dmul <- function(d1, d2) {
    res <- list()
    for (k1 in names(d1)) for (k2 in names(d2)) {
      exk <- .add(list(d1[[k1]][[1L]], d2[[k2]][[1L]]))
      if (.ns_lt(exk, n)) {
        kk <- exk$key
        res[[kk]] <- list(exk, .add(list(if (is.null(res[[kk]])) S_Zero else res[[kk]][[2L]], .mul(list(d1[[k1]][[2L]], d2[[k2]][[2L]])))))
      }
    }
    res
  }
  k <- S_One; terms <- list(); pk <- pterms
  while (.ns_lt(.mul(list(k, d)), n)) {
    cf <- .neg(.div(.pow(S_NegativeOne, k), k))
    for (kk in names(pk)) terms[[kk]] <- list(pk[[kk]][[1L]], .add(list(if (is.null(terms[[kk]])) S_Zero else terms[[kk]][[2L]], .mul(list(cf, pk[[kk]][[2L]])))))
    pk <- dmul(pk, pterms)
    k <- .add(list(k, S_One))
  }
  res <- .add(list(fcall("log", list(a)), .neg(.mul(list(b, fcall("log", list(cdir))))), .mul(list(b, lx))))
  for (tt in terms) res <- .add(list(res, .mul(list(cancel(tt[[2L]]), .pow(t, tt[[1L]])))))
  # a negative leading coefficient with a nonreal approach crosses the cut
  if (isTRUE(ask_fact(a, "negative")) && !is_zero_num(fcall("im", list(z)))) {
    ls <- .so_lseries(z, t, 6L)
    i <- 0L; term <- NULL
    for (term in ls) { if (!isTRUE(ask_fact(term, "real")) || i == 5L) break; i <- i + 1L }
    if (i < 5L && !is.null(term)) {
      cf <- .so_as_coeff_exponent(term, t)[[1L]]
      res <- .add(list(res, .mul(list(.int(-2L), S_I, S_Pi, fcall("Heaviside", list(.neg(fcall("im", list(cf))), S_Zero))))))
    }
  }
  res <- subs(res, t, .div(x, cdir))
  .add(list(res, .so_new(.pow(x, n), list(x))))
}
.so_alt_hooks$log <- function(e, x, logx, cdir) {
  if (length(e$args) != 1L) return(.so_alt_function(e, x, logx, cdir))
  arg0 <- together(e$args[[1L]])
  t <- Dummy("t", positive = TRUE)
  cdir <- as_sym(cdir)
  if (is_zero_num(cdir)) cdir <- S_One
  z <- subs(arg0, x, .mul(list(cdir, t)))
  ce <- tryCatch(.so_leadterm(z, t, logx, 1), ValueError = function(err) NULL)
  if (is.null(ce)) return(fcall("log", list(as_leading_term(arg0, x, logx, cdir))))
  c <- ce[[1L]]; ex <- ce[[2L]]
  if (.has1(c, t)) {
    c <- subs(c, t, .div(x, cdir))
    if (!is_zero_num(ex)) .pole_error(paste("Cannot expand", .str(e), "around 0"))
    return(fcall("log", list(c)))
  }
  if (is_one_num(c) && is_zero_num(ex)) return(as_leading_term(.sub(arg0, S_One), x, logx))
  res <- .sub(fcall("log", list(c)), .mul(list(ex, fcall("log", list(cdir)))))
  lx <- if (is.null(logx)) fcall("log", list(x)) else logx
  res <- .add(list(res, .mul(list(ex, lx))))
  # STEP 5: a negative leading coefficient with a nonreal approach
  if (isTRUE(ask_fact(c, "negative")) && !is_zero_num(fcall("im", list(z)))) {
    terms <- .so_lseries(z, t, 6L)
    i <- 0L; term <- NULL
    for (term in terms) { if (!isTRUE(ask_fact(term, "real")) || i == 5L) break; i <- i + 1L }
    if (i < 5L && !is.null(term)) {
      cf <- .so_as_coeff_exponent(term, t)[[1L]]
      res <- .add(list(res, .mul(list(.int(-2L), S_I, S_Pi, fcall("Heaviside", list(.neg(fcall("im", list(cf))), S_Zero))))))
    }
  }
  res
}

# Expr.lseries with the default _eval_lseries: the first k terms, in order.
# SymPy raises n without bound; internal callers stop at n = nmax.
.so_lseries <- function(e, x, k, logx = NULL, cdir = 0, nmax = 30L) {
  out <- list()
  push <- function(si) {
    if (si$type != "Add") { out[[length(out) + 1L]] <<- si; return(invisible()) }
    # yield_lseries: the terms of a sum one order at a time
    yielded <- S_Zero
    o <- .mul(list(.so_new(si, list(x)), x))
    ndid <- 0L; ndo <- length(si$args)
    for (guard in 1:40) {
      do <- .so_removeO(.add(list(si, .neg(yielded), o)))
      o <- .mul(list(o, x))
      if (is_zero_num(do) || do$type == "Order") next
      ndid <- ndid + (if (do$type == "Add") length(do$args) else 1L)
      out[[length(out) + 1L]] <<- do
      if (ndid == ndo || length(out) >= k) break
      yielded <- .add(list(yielded, do))
    }
  }
  n <- 0L
  sr <- .ns_eval(e, x, .int(n), logx, cdir)
  while (sr$type == "Order" && n < nmax) { n <- n + 1L; sr <- .ns_eval(e, x, .int(n), logx, cdir) }
  cur <- .so_removeO(sr)
  push(cur)
  if (is_zero_num(cur)) return(out)
  while (length(out) < k && n < nmax) {
    repeat {
      n <- n + 1L
      sr <- .so_removeO(.ns_eval(e, x, .int(n), logx, cdir))
      if (!identical(sr$key, cur$key) || n >= nmax) break
      if (is_zero_num(cancel(.sub(sr, e)))) return(out)
    }
    if (identical(sr$key, cur$key)) break
    push(.sub(sr, cur))
    cur <- sr
  }
  out[seq_len(min(k, length(out)))]
}
