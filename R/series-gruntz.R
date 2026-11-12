# The Gruntz algorithm for limits, after sympy/series/gruntz.py. Limits are
# moved to x -> oo; the most rapidly varying subexpressions (the mrv set)
# are found, rewritten in terms of a single w -> 0+, and the leading term of
# the series in w decides the limit.
#' Limits by the Gruntz algorithm
#'
#' `gruntz(e, z, z0, dir)` computes the limit of `e` as `z` approaches
#' `z0` with the Gruntz algorithm, as SymPy's `gruntz()`. `z0` may be `oo`
#' or `-oo`; `dir` (`"+"` or `"-"`) gives the side for finite `z0`.
#' [limit()] tries this algorithm first.
#'
#' @param e An expression.
#' @param z The variable, a symbol.
#' @param z0 The limit point.
#' @param dir Direction, `"+"` or `"-"`.
#' @return An expression.
#' @examples
#' x <- Symbol("x")
#' gruntz(exp(x) * (exp(1/x - exp(-x)) - exp(1/x)), x, oo)
#' gruntz((3^x + 5^x)^(1/x), x, oo)
#' gruntz(x^x, x, 0)
#' @export
gruntz <- function(e, z, z0, dir = "+") {
  e <- as_sym(e); z <- as_sym(z); z0 <- as_sym(z0)
  if (!(z$type %in% c("Symbol", "Dummy"))) stop("Second argument must be a Symbol", call. = FALSE)
  if (is.null(.gr_state$cache)) {
    .gr_state$cache <- new.env(hash = TRUE); .gr_state$depth <- 0L
    on.exit(.gr_state$cache <- NULL)
  }
  if (.has_type(e, "Function") && any(c("Min", "Max") %in% .fun_names(e))) e <- rewrite(e, "Abs")
  e0 <- if (is_const(z0, "oo")) e
        else if (is_const(z0, "-oo")) subs(e, z, .neg(z))
        else if (dir == "-") subs(e, z, .sub(z0, .inv(z)))
        else if (dir == "+") subs(e, z, .add(list(z0, .inv(z))))
        else stop("dir must be '+' or '-'", call. = FALSE)
  r <- .gr_limitinf(e0, z)
  .gr_intractable(r)
}

.gr_state <- new.env()
.fun_names <- function(e) {
  out <- character()
  walk <- function(z) { if (z$type == "Function") out <<- c(out, z$name); for (a in z$args) walk(a) }
  walk(e)
  unique(out)
}
.gr_memo <- function(tag, e, x, f) {
  cache <- .gr_state$cache
  if (is.null(cache)) return(f())
  k <- paste(tag, e$key, x$key)
  if (!is.null(r <- cache[[k]])) return(r)
  r <- f()
  assign(k, r, envir = cache)
  r
}

.gr_is_exp <- function(e) (e$type == "Function" && e$name == "exp") || (e$type == "Pow" && is_const(e$args[[1L]], "E"))
.gr_exp_arg <- function(e) if (e$type == "Function") e$args[[1L]] else e$args[[2L]]
.gr_exp <- function(a) fcall("exp", list(a))
.gr_log <- function(a) fcall("log", list(a))
.gr_rep <- function(e, old, new) xreplace(e, list(list(old, new)))
.gr_infinite <- function(v) {
  fs <- if (v$type == "Mul") v$args else list(v)
  any(vapply(fs, is_inf_const, TRUE))
}

# SubsSet: (expr, dummy) pairs in insertion order, plus rewrites keyed by dummy.
.ss_new <- function() list(exprs = list(), vars = list(), rewrites = list())
.ss_keys <- function(s) vapply(s$exprs, .k, "")
.ss_get <- function(s, e) {
  i <- match(e$key, .ss_keys(s))
  if (!is.na(i)) return(list(s = s, var = s$vars[[i]]))
  v <- Dummy("d")
  s$exprs <- c(s$exprs, list(e)); s$vars <- c(s$vars, list(v))
  list(s = s, var = v)
}
.ss_do_subs <- function(s, e) {
  for (i in seq_along(s$exprs)) e <- .gr_rep(e, s$vars[[i]], s$exprs[[i]])
  e
}
.ss_union <- function(a, b, exps = NULL) {
  res <- a
  tr <- list()
  ka <- .ss_keys(a)
  for (i in seq_along(b$exprs)) {
    j <- match(b$exprs[[i]]$key, ka)
    if (!is.na(j)) {
      if (!is.null(exps)) exps <- .gr_rep(exps, b$vars[[i]], res$vars[[j]])
      tr[[length(tr) + 1L]] <- list(b$vars[[i]], res$vars[[j]])
    } else {
      res$exprs <- c(res$exprs, b$exprs[i]); res$vars <- c(res$vars, b$vars[i])
    }
  }
  for (k in names(b$rewrites)) res$rewrites[[k]] <- if (length(tr)) xreplace(b$rewrites[[k]], tr) else b$rewrites[[k]]
  list(s = res, exps = exps)
}

.gr_compare <- function(a, b, x) {
  la <- if (.gr_is_exp(a)) .gr_exp_arg(a) else .gr_log(a)
  lb <- if (.gr_is_exp(b)) .gr_exp_arg(b) else .gr_log(b)
  c <- .gr_limitinf(.div(la, lb), x)
  if (is_zero_num(c)) "<" else if (isTRUE(ask_fact(c, "infinite"))) ">" else "="
}

.mrv_max3 <- function(f, expsf, g, expsg, union, expsboth, x) {
  if (!length(f$exprs)) return(list(s = g, e = expsg))
  if (!length(g$exprs)) return(list(s = f, e = expsf))
  if (any(.ss_keys(f) %in% .ss_keys(g))) return(list(s = union, e = expsboth))
  switch(.gr_compare(f$exprs[[1L]], g$exprs[[1L]], x),
         ">" = list(s = f, e = expsf), "<" = list(s = g, e = expsg), list(s = union, e = expsboth))
}
.mrv_max1 <- function(f, g, exps, x) {
  u <- .ss_union(f, g, exps)
  .mrv_max3(f, .ss_do_subs(g, exps), g, .ss_do_subs(f, exps), u$s, u$exps, x)
}

# The mrv set of e and e written in terms of its dummies.
.mrv <- function(e, x) {
  e <- powsimp(e, deep = TRUE, combine = "exp")
  if (!.has1(e, x)) return(list(s = .ss_new(), e = e))
  if (identical(e$key, x$key)) { g <- .ss_get(.ss_new(), x); return(list(s = g$s, e = g$var)) }
  if (e$type %in% c("Mul", "Add")) {
    op <- if (e$type == "Mul") .mul else .add
    sargs <- .sympy_args(e)
    dep <- Filter(function(a) .has1(a, x), sargs)
    ind <- op(Filter(function(a) !.has1(a, x), sargs))
    if (length(dep) == 1L) { r <- .mrv(dep[[1L]], x); return(list(s = r$s, e = op(list(ind, r$e)))) }
    r1 <- .mrv(dep[[1L]], x); r2 <- .mrv(op(dep[-1L]), x)
    return(.mrv_max1(r1$s, r2$s, op(list(ind, r1$e, r2$e)), x))
  }
  if (e$type == "Pow" && !is_const(e$args[[1L]], "E")) {
    e1 <- S_One; b1 <- e
    while (b1$type == "Pow") { e1 <- .mul(list(e1, b1$args[[2L]])); b1 <- b1$args[[1L]] }
    if (is_one_num(b1)) return(list(s = .ss_new(), e = b1))
    if (.has1(e1, x)) return(.mrv(.gr_exp(.mul(list(e1, .gr_log(b1)))), x))
    r <- .mrv(b1, x)
    return(list(s = r$s, e = .pow(r$e, e1)))
  }
  if (e$type == "Function" && e$name == "log") {
    r <- .mrv(e$args[[1L]], x)
    return(list(s = r$s, e = .gr_log(r$e)))
  }
  if (.gr_is_exp(e)) {
    a <- .gr_exp_arg(e)
    if (a$type == "Function" && a$name == "log") return(.mrv(a$args[[1L]], x))
    if (.gr_infinite(.gr_limitinf(a, x))) {
      g <- .ss_get(.ss_new(), e)
      r2 <- .mrv(a, x)
      su <- .ss_union(g$s, r2$s)$s
      su$rewrites[[g$var$key]] <- .gr_exp(r2$e)
      return(.mrv_max3(g$s, g$var, r2$s, .gr_exp(r2$e), su, g$var, x))
    }
    r <- .mrv(a, x)
    return(list(s = r$s, e = .gr_exp(r$e)))
  }
  if (e$type == "Function") {
    if (is.null(.FUNS[[e$name]])) .gr_value_error("MRV set computation for UndefinedFunction is not allowed")
    if (e$name %in% c("Min", "Max")) return(.mrv(rewrite(e, "Abs"), x))
    l <- lapply(e$args, .mrv, x = x)
    l2 <- Filter(function(r) length(r$s$exprs) > 0L, l)
    if (length(l2) != 1L) .notimpl_error("MRV set computation for functions in several variables not implemented.")
    return(list(s = l2[[1L]]$s, e = sym_rebuild(e, lapply(l, `[[`, "e"))))
  }
  if (e$type == "Derivative") .notimpl_error("MRV set computation for derivatives not implemented yet.")
  .notimpl_error(paste0("Don't know how to calculate the mrv of '", .str(e), "'"))
}

# Sign of e(x) for large x, as an Integer (or an unevaluated sign).
.gr_sign <- function(e, x) .gr_memo("sign", e, x, function() {
  if (isTRUE(ask_fact(e, "positive"))) return(S_One)
  if (isTRUE(ask_fact(e, "negative"))) return(S_NegativeOne)
  if (isTRUE(ask_fact(e, "zero"))) return(S_Zero)
  if (!.has1(e, x)) {
    # _sign of logcombine(e), with certified evaluation and an exact zero test
    e2 <- logcombine(e)
    s <- .gr_const_sign(e2)
    return(if (is.na(s)) fcall("sign", list(e2)) else .int(s))
  }
  if (identical(e$key, x$key)) return(S_One)
  if (e$type == "Mul") {
    sa <- .gr_sign(e$args[[1L]], x)
    if (is_zero_num(sa)) return(S_Zero)
    return(.mul(list(sa, .gr_sign(.mul(e$args[-1L]), x))))
  }
  if (.gr_is_exp(e)) return(S_One)
  if (e$type == "Pow") {
    s <- .gr_sign(e$args[[1L]], x)
    if (is_one_num(s)) return(S_One)
    if (e$args[[2L]]$type == "Integer") return(.pow(s, e$args[[2L]]))
  }
  if (e$type == "Function" && e$name == "log" && isTRUE(ask_fact(e$args[[1L]], "positive")))
    return(.gr_sign(.sub(e$args[[1L]], S_One), x))
  .gr_sign(.mrv_leadterm(e, x)[[1L]], x)
})

.gr_value_error <- function(msg) stop(structure(class = c("ValueError", "error", "condition"), list(message = msg, call = NULL)))

# limitinf
.gr_limitinf <- function(e, x) .gr_memo("lim", e, x, function() {
  old <- e
  if (!.has1(e, x)) return(e)
  if (.has_type(e, "Order")) e <- removeO(expand(e))
  if (!isTRUE(ask_fact(x, "positive")) || isTRUE(ask_fact(x, "integer"))) {
    p <- Dummy("p", positive = TRUE)
    e <- subs(e, x, p); x <- p
  }
  e <- powdenest(.gr_tractable(e, x))
  if (e$type == "AccumulationBounds") {
    lmin <- .mrv_leadterm(e$args[[1L]], x); lmax <- .mrv_leadterm(e$args[[2L]], x)
    if (!identical(lmin[[1L]]$key, lmax[[1L]]$key) || !identical(lmin[[2L]]$key, lmax[[2L]]$key)) .notimpl_error("")
    lt <- lmin
  } else lt <- .mrv_leadterm(e, x)
  c0 <- lt[[1L]]
  sig <- .gr_sign(lt[[2L]], x)
  if (is_one_num(sig)) return(S_Zero)
  if (identical(sig$key, S_NegativeOne$key)) {
    if (.gr_is_imag_multiple(c0)) return(.mul(list(c0, S_Infinity)))
    s <- .gr_sign(c0, x)
    if (is_zero_num(s)) .gr_value_error("Leading term should not be 0")
    return(.mul(list(s, S_Infinity)))
  }
  if (is_zero_num(sig)) {
    if (identical(c0$key, old$key)) c0 <- cancel(c0)
    return(.gr_limitinf(c0, x))
  }
  .gr_value_error(paste(.str(sig), "could not be evaluated"))
})

# rewrite("tractable", deep=True, limitvar=x): the _eval_rewrite_as_tractable
# methods of the gamma family, error functions, exponential integrals,
# factorials, harmonic and fibonacci numbers, hyperbolic functions and the
# inverse trigonometric and hyperbolic functions (as log)
.GR_LOGREW <- c("asin", "acos", "atan", "acot", "asec", "acsc", "asinh", "acosh", "atanh", "acoth")
.gr_tractable <- function(e, x = NULL) .map_bottom_up(e, function(z) {
  if (z$type != "Function") return(z)
  r <- .gr_tractable1(z, x)
  if (is.null(r)) z else r
})
.gr_tractable1 <- function(z, x) {
  a <- z$args[[1L]]
  rec <- function(v) .gr_tractable(v, x)
  switch(z$name,
    gamma = .gr_exp(fcall("loggamma", list(a))),
    uppergamma = .sub(.gr_exp(fcall("loggamma", list(a))), fcall("lowergamma", z$args)),
    erf = {
      if (!is.null(x) && is_const(limit(a, x, S_Infinity), "-oo"))
        .add(list(S_NegativeOne, .mul(list(fcall("_erfs", list(.neg(a))), .gr_exp(.neg(.sq(a)))))))
      else .sub(S_One, .mul(list(fcall("_erfs", list(a)), .gr_exp(.neg(.sq(a))))))
    },
    erfc = rec(.sub(S_One, fcall("erf", list(a)))),
    erfi = rec(.mul(list(.neg(S_I), fcall("erf", list(.mul(list(S_I, a))))))),
    Ei = .mul(list(.gr_exp(a), fcall("_eis", list(a)))),
    li = .mul(list(a, fcall("_eis", list(.gr_log(a))))),
    Li = rec(.sp_li_as_li(a)),
    fibonacci = if (length(z$args) == 1L) .div(.sub(.pow(S_GoldenRatio, a), .div(cos(.mul(list(S_Pi, a))), .pow(S_GoldenRatio, a))), .sqrt(.int(5L))) else NULL,
    harmonic = {
      pg <- .harmonic_as_polygamma(z)
      if (!(pg$type == "Function" && pg$name == "harmonic")) rec(pg)
      else {
        m <- if (length(z$args) > 1L) z$args[[2L]] else S_One
        if (isTRUE(ask_fact(.sub(m, S_One), "nonzero"))) rec(.sub(fcall("zeta", list(m)), fcall("zeta", list(m, .add(list(a, S_One)))))) else NULL
      }
    },
    RisingFactorial = , FallingFactorial = {
      k <- z$args[[2L]]
      kl <- if (is.null(x)) NULL else subs(k, x, S_Infinity)
      if (!is.null(kl) && is_const(kl, "oo")) {
        if (z$name == "RisingFactorial") .div(rec(fcall("gamma", list(.add(list(a, k))))), fcall("gamma", list(a)))
        else .div(.mul(list(.pow(S_NegativeOne, k), rec(fcall("gamma", list(.sub(k, a)))))), fcall("gamma", list(.neg(a))))
      } else if (!is.null(kl) && is_const(kl, "-oo")) {
        if (z$name == "RisingFactorial") .div(.mul(list(.pow(S_NegativeOne, k), fcall("gamma", list(.sub(S_One, a))))), rec(fcall("gamma", list(.add(list(.neg(k), .neg(a), S_One))))))
        else .div(fcall("gamma", list(.add(list(a, S_One)))), rec(fcall("gamma", list(.add(list(a, .neg(k), S_One))))))
      } else rec(rewrite(z, "gamma"))
    },
    binomial = rec(rewrite(z, "gamma")),
    besseli = if (isTRUE(ask_fact(z$args[[2L]], "extended_real"))) .mul(list(.gr_exp(z$args[[2L]]), fcall("_besseli", z$args))) else NULL,
    besselk = if (isTRUE(ask_fact(z$args[[2L]], "extended_real"))) .mul(list(.gr_exp(.neg(z$args[[2L]])), fcall("_besselk", z$args))) else NULL,
    sinh = .div(.sub(.gr_exp(a), .gr_exp(.neg(a))), S_Two),
    cosh = .div(.add(list(.gr_exp(a), .gr_exp(.neg(a)))), S_Two),
    tanh = .div(.sub(.gr_exp(a), .gr_exp(.neg(a))), .add(list(.gr_exp(a), .gr_exp(.neg(a))))),
    coth = .div(.add(list(.gr_exp(a), .gr_exp(.neg(a)))), .sub(.gr_exp(a), .gr_exp(.neg(a)))),
    csch = .inv(.div(.sub(.gr_exp(a), .gr_exp(.neg(a))), S_Two)),
    sech = .inv(.div(.add(list(.gr_exp(a), .gr_exp(.neg(a)))), S_Two)),
    if (z$name %in% .GR_LOGREW) rewrite(z, "log") else NULL)
}
# harmonic._eval_rewrite_as_polygamma
.harmonic_as_polygamma <- function(z) {
  n <- z$args[[1L]]; m <- if (length(z$args) > 1L) z$args[[2L]] else S_One
  if (isTRUE(ask_fact(m, "integer")) && isTRUE(ask_fact(m, "positive")))
    return(.mul(list(.div(.pow(S_NegativeOne, m), fcall("factorial", list(.sub(m, S_One)))),
                     .sub(fcall("polygamma", list(.sub(m, S_One), S_One)), fcall("polygamma", list(.sub(m, S_One), .add(list(n, S_One))))))))
  z
}
.rewrite_rules$polygamma$harmonic <- function(n, m = S_One) { r <- .harmonic_as_polygamma(fcall("harmonic", list(n, m))); if (r$type == "Function" && r$name == "harmonic") NULL else r }

# rewrite("intractable", deep=True)
.gr_intractable <- function(e) .map_bottom_up(e, function(z) {
  if (z$type != "Function") return(z)
  a <- z$args[[1L]]
  switch(z$name,
    loggamma = .gr_log(fcall("gamma", list(a))),
    `_erfs` = .mul(list(.sub(S_One, fcall("erf", list(a))), .gr_exp(.sq(a)))),
    `_eis` = .mul(list(.gr_exp(.neg(a)), fcall("Ei", list(a)))),
    `_besseli` = .mul(list(.gr_exp(.neg(z$args[[2L]])), fcall("besseli", z$args))),
    `_besselk` = .mul(list(.gr_exp(z$args[[2L]]), fcall("besselk", z$args))),
    z)
})

# c0.match(I*Wild("a", exclude=[I]))
.gr_is_imag_multiple <- function(c) {
  if (identical(c$key, S_I$key)) return(TRUE)
  c$type == "Mul" && any(vapply(c$args, function(a) identical(a$key, S_I$key), TRUE)) &&
    !.contains_key(.div(c, S_I), S_I$key)
}

# mrv_leadterm
.mrv_leadterm <- function(e, x) .gr_memo("lt", e, x, function() {
  if (!.has1(e, x)) return(list(e, S_Zero))
  r <- .mrv(e, x)
  Omega <- r$s; exps <- r$e
  if (!length(Omega$exprs)) return(list(exps, S_Zero))
  if (x$key %in% .ss_keys(Omega)) {
    ex <- .gr_exp(x)
    Omega$exprs <- lapply(Omega$exprs, .gr_rep, old = x, new = ex)
    Omega$rewrites <- lapply(Omega$rewrites, .gr_rep, old = x, new = ex)
    exps <- .gr_rep(exps, x, ex)
  }
  w <- Dummy("w", positive = TRUE)
  rw <- .gr_rewrite(exps, Omega, x, w)
  f <- .map_bottom_up(rw$f, function(z) {
    if (z$type == "Pow" && .has1(z, x)) .gr_exp(.mul(list(.gr_log(z$args[[1L]]), z$args[[2L]]))) else z
  })
  logw <- rw$logw
  lt <- tryCatch(.so_leadterm(f, w, logw), NotImplementedError = function(err) NULL, PoleError = function(err) NULL,
                 ValueError = function(err) NULL)
  if (is.null(lt)) {
    n0 <- 1L; incr <- 1L
    sr <- .so_new(S_One)
    while (.so_is_order(sr)) {
      sr <- .ns_eval(f, w, .int(n0 + incr), logw, 0)
      incr <- incr * 2L
    }
    sr <- .so_removeO(expand(sr))
    lt <- tryCatch(.so_leadterm(sr, w, logw), NotImplementedError = function(err) NULL, PoleError = function(err) NULL,
                   ValueError = function(err) NULL)
    if (is.null(lt)) {
      lt <- .so_as_coeff_exponent(f, w)
      if (.has1(lt[[1L]], w)) {
        be <- .as_base_exp(f)
        b <- .so_as_coeff_exponent(be[[1L]], w)
        lt <- list(.pow(b[[1L]], be[[2L]]), .mul(list(b[[2L]], be[[2L]])))
      }
    }
  }
  list(subs(lt[[1L]], .gr_log(w), logw), lt[[2L]])
})

# Sign of a real constant: -1, 0, 1, or NA when unknown.
.gr_const_sign <- function(d) {
  if (is_zero_num(d)) return(0)
  if (isTRUE(ask_fact(d, "positive"))) return(1)
  if (isTRUE(ask_fact(d, "negative"))) return(-1)
  if (length(free_symbols(d))) return(NA)
  s <- .mp_num_sign(d)
  if (!is.na(s)) return(s)
  # a constant that evaluates to no significant digits may be an exact
  # zero: logs of powers combined, as logcombine before _sign
  z <- tryCatch(simplify(expand_log(d, force = TRUE)), error = function(err) d)
  if (is_zero_num(z)) return(0)
  NA
}

.gr_rewrite <- function(e, Omega, x, wsym) {
  if (!length(Omega$exprs)) .gr_value_error("Length cannot be 0")
  if (!all(vapply(Omega$exprs, .gr_is_exp, TRUE))) .gr_value_error("Value should be exp")
  rewrites <- Omega$rewrites
  n <- length(Omega$exprs)
  vkeys <- vapply(Omega$vars, .k, "")
  before <- lapply(seq_len(n), function(i) {
    r <- rewrites[[vkeys[i]]]
    if (is.null(r)) integer() else which(vapply(Omega$vars, function(v) .has1(r, v), TRUE))
  })
  ht <- function(i) 1 + sum(vapply(before[[i]], ht, 0))
  ord <- order(-vapply(seq_len(n), ht, 0))
  exprs <- Omega$exprs[ord]; vars <- Omega$vars[ord]
  for (g in exprs) {
    sig <- .gr_sign(.gr_exp_arg(g), x)
    if (!is_one_num(sig) && !identical(sig$key, S_NegativeOne$key) && !.has_type(sig, "AccumulationBounds"))
      .notimpl_error(paste("Result depends on the sign of", .str(sig)))
  }
  w <- wsym
  if (is_one_num(sig)) wsym <- .inv(wsym)
  garg <- .gr_exp_arg(g)
  O2 <- list(); dens <- gmp::as.bigz(1)
  for (i in seq_along(exprs)) {
    c <- .gr_limitinf(.div(.gr_exp_arg(exprs[[i]]), garg), x)
    if (c$type %in% c("Integer", "Rational")) dens <- gmp::lcm.bigz(dens, gmp::denominator(gmp::as.bigq(c$value)))
    arg <- .gr_exp_arg(exprs[[i]])
    rv <- rewrites[[vars[[i]]$key]]
    if (!is.null(rv)) {
      if (!.gr_is_exp(rv)) .gr_value_error("Value should be exp")
      arg <- .gr_exp_arg(rv)
    }
    O2[[i]] <- list(vars[[i]], .mul(list(.gr_exp(.sub(arg, .mul(list(c, garg)))), .pow(wsym, c))))
  }
  f <- powsimp(e, deep = TRUE, combine = "exp")
  for (p in O2) f <- .gr_rep(f, p[[1L]], p[[2L]])
  logw <- if (is_one_num(sig)) .neg(garg) else garg
  if (dens != 1) {
    f <- subs(f, w, .pow(w, .int(as.character(dens))))
    logw <- .div(logw, .int(as.character(dens)))
  }
  # bottom_up(f, Expr.normal)
  f <- .map_bottom_up(f, function(z) {
    if (!length(z$args) || !(z$type %in% c("Add", "Mul", "Pow", "Function"))) return(z)
    .expr_normal(z)
  })
  list(f = f, logw = logw)
}
