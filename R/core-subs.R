# Basic._subs and the _eval_subs methods of Number, Add, Mul, Pow and exp,
# after sympy/core/basic.py, add.py, mul.py, power.py, numbers.py and
# functions/elementary/exponential.py.

.subs_number_like <- function(e) is_number_type(e) || (e$type == "Constant" && e$name %in% c("oo", "-oo", "zoo", "nan"))

.subs1 <- function(e, old, new) {
  if (identical(e$key, old$key)) return(new)
  if (!length(e$args)) {
    # Number._eval_subs: (-2).subs(2, y) = -y
    if (.subs_number_like(e) && .subs_number_like(old) && identical(.neg(e)$key, old$key)) return(.neg(new))
    # Symbol._eval_subs: x as Pow(x, 1) against a power, x.subs(1/x, y) = 1/y
    if (e$type %in% c("Symbol", "Dummy") && old$type == "Pow") {
      r <- .pow_eval_subs(e, S_One, old, new, old$args, identical(old$args[[1L]]$key, S_E$key))
      if (!is.null(r)) return(r)
    }
    return(e)
  }
  # no subtree can match when a free symbol of old is not free in e (with
  # integral transforms, whose variables are substituted, by containment)
  if (isTRUE(.subs_state$structural)) {
    if (old$type %in% c("Symbol", "Dummy") && !.contains_key(e, old$key)) return(e)
  } else {
    ofs <- names(.free_sym_list(old))
    if (length(ofs) && !all(ofs %in% names(.free_sym_list(e)))) return(e)
  }
  rv <- .eval_subs(e, old, new)
  if (is.null(rv)) rv <- .subs_fallback(e, old, new)
  rv
}

.subs_fallback <- function(e, old, new) {
  args <- lapply(e$args, .subs1, old, new)
  for (i in seq_along(args)) if (!identical(args[[i]]$key, e$args[[i]]$key)) return(sym_rebuild(e, args))
  e
}

.eval_subs <- function(e, old, new) {
  h <- .subs_hooks[[e$type]]
  if (!is.null(h)) {
    r <- h(e, old, new)
    if (!is.null(r)) return(r)
  }
  switch(e$type,
    Add = .add_eval_subs(e, old, new),
    Mul = .mul_eval_subs(e, old, new),
    Pow = .pow_eval_subs(e$args[[1L]], e$args[[2L]], old, new,
                         if (old$type == "Pow") old$args else NULL,
                         (old$type == "Function" && old$name == "exp") || (old$type == "Pow" && identical(old$args[[1L]]$key, S_E$key))),
    Function = if (e$name == "exp") .exp_eval_subs(e, old, new) else NULL,
    NULL)
}

# ---------------------------------------------------------------------------
# Add

.subs_is_rat <- function(e) e$type %in% c("Integer", "Rational")
.subs_as_coeff_Add <- function(e) {
  if (e$type == "Add") {
    a <- .sympy_args(e)
    if (is_number_type(a[[1L]])) return(list(a[[1L]], .add(a[-1L])))
    return(list(S_Zero, e))
  }
  if (is_number_type(e)) return(list(e, S_Zero))
  list(S_Zero, e)
}
.subs_add_args <- function(e) if (e$type == "Add") e$args else list(e)

.add_eval_subs <- function(e, old, new) {
  if (old$type != "Add") {
    if (identical(old$key, S_Infinity$key) && any(vapply(e$args, function(a) identical(a$key, S_NegativeInfinity$key), NA)))
      return(xreplace(e, list(list(S_NegativeInfinity, .neg(new)))))
    return(NULL)
  }
  cs <- .subs_as_coeff_Add(e); co <- .subs_as_coeff_Add(old)
  coeff_self <- cs[[1L]]; terms_self <- cs[[2L]]; coeff_old <- co[[1L]]; terms_old <- co[[2L]]
  if (.subs_is_rat(coeff_self) && .subs_is_rat(coeff_old)) {
    if (identical(terms_self$key, terms_old$key)) return(.add(list(new, coeff_self, .neg(coeff_old))))
    if (identical(terms_self$key, .neg(terms_old)$key)) return(.add(list(.neg(new), coeff_self, coeff_old)))
  }
  if ((.subs_is_rat(coeff_self) && .subs_is_rat(coeff_old)) || identical(coeff_self$key, coeff_old$key)) {
    args_old <- .subs_add_args(terms_old); args_self <- .subs_add_args(terms_self)
    if (length(args_old) < length(args_self)) {
      sk <- vapply(args_self, .k, "")
      proper <- function(ok) all(ok %in% sk) && length(setdiff(unique(sk), ok)) > 0L
      ok <- vapply(args_old, .k, "")
      if (proper(ok)) {
        rest <- args_self[!(sk %in% ok)]
        return(.add(c(list(new, coeff_self, .neg(coeff_old)), lapply(rest, .subs1, old, new))))
      }
      ok <- vapply(.subs_add_args(.neg(terms_old)), .k, "")
      if (proper(ok)) {
        rest <- args_self[!(sk %in% ok)]
        return(.add(c(list(.neg(new), coeff_self, coeff_old), lapply(rest, .subs1, old, new))))
      }
    }
  }
  NULL
}

# ---------------------------------------------------------------------------
# Mul

.subs_base_exp <- function(a) {
  if (a$type == "Pow") return(list(a$args[[1L]], a$args[[2L]]))
  if (a$type == "Function" && a$name == "exp") return(list(S_E, a$args[[1L]]))
  list(a, S_One)
}
# Expr.as_coeff_mul(): the numeric coefficient and the other factors
.subs_coeff_mul <- function(e) {
  if (e$type == "Mul") {
    a <- .sympy_args(e)
    if (is_number_type(a[[1L]])) return(list(a[[1L]], a[-1L]))
    return(list(S_One, a))
  }
  if (is_number_type(e)) return(if (is_zero_num(e)) list(S_One, list(e)) else list(e, list()))
  list(S_One, list(e))
}
# Mul._eval_subs.breakup for commutative factors: base key -> (base, exponent)
.subs_breakup <- function(eq) {
  out <- list()
  for (a in .margs(eq)) {
    a <- powdenest(a)
    be <- .subs_base_exp(a)
    b <- be[[1L]]; ex <- be[[2L]]
    if (!identical(ex$key, S_One$key)) {
      co <- .subs_coeff_mul(ex)[[1L]]
      b <- .pow(b, .div(ex, co)); ex <- co
    }
    k <- b$key
    out[[k]] <- if (is.null(out[[k]])) list(b, ex) else list(b, .add(list(out[[k]][[2L]], ex)))
  }
  out
}
.subs_rejoin <- function(b, co) { be <- .subs_base_exp(b); .pow(be[[1L]], .mul(list(be[[2L]], co))) }
.subs_q <- function(a) if (a$type == "Integer") gmp::as.bigz(1) else gmp::denominator(a$value)
.subs_ndiv <- function(a, b) {
  if (!.subs_is_rat(a) || !.subs_is_rat(b)) return(0L)
  qa <- .subs_q(a); qb <- .subs_q(b)
  if (qa %% qb == 0 || qb %% qa == 0) {
    v <- a$value / b$value
    return(as.integer(as.character(sign(v) * floor(abs(v)))))
  }
  0L
}

.mul_eval_subs <- function(e, old, new) {
  if (old$type != "Mul") return(NULL)
  oa <- .sympy_args(old); ea <- .sympy_args(e)
  if (is_number_type(oa[[1L]]) && isTRUE(ask_fact(oa[[1L]], "negative"))) {
    if (is_number_type(ea[[1L]])) {
      if (isTRUE(ask_fact(ea[[1L]], "negative"))) return(.subs1(e, .neg(old), .neg(new)))
      return(NULL)
    }
  }
  # ponytail: the noncommutative branch of Mul._eval_subs is not ported;
  # products with noncommutative factors fall back to argument-wise subs
  if (!.is_comm(e) || !.is_comm(old)) return(NULL)
  rv <- NULL
  fr <- fraction(e)
  self2 <- e
  if (!is_one_num(fr[[2L]])) {
    self2 <- .div(.subs1(fr[[1L]], old, new), .subs1(fr[[2L]], old, new))
    if (self2$type != "Mul") return(.subs1(self2, old, new))
    if (!identical(self2$key, e$key)) rv <- self2
  }
  co_self <- .sympy_args(self2)[[1L]]
  co_old <- oa[[1L]]
  co_xmul <- NULL
  if (.subs_is_rat(co_old) && .subs_is_rat(co_self)) {
    if (!identical(co_old$key, co_self$key)) co_xmul <- extract_multiplicatively(co_self, co_old)
  } else if (.subs_is_rat(co_old)) return(rv)
  c_ <- .subs_breakup(self2); old_c <- .subs_breakup(old)
  co_residual <- S_One
  if (!is.null(co_xmul) && .subs_is_rat(co_xmul) && !is_one_num(abs(co_old))) {
    mult <- as_sym(multiplicity(abs(co_old), co_self))
    c_[[co_self$key]] <- NULL
    k <- co_old$key
    c_[[k]] <- if (is.null(c_[[k]])) list(co_old, mult) else list(co_old, .add(list(c_[[k]][[2L]], mult)))
    co_residual <- .div(co_self, .pow(co_old, mult))
  }
  if (length(old_c) > length(c_)) return(rv)
  if (length(setdiff(names(old_c), names(c_)))) return(rv)
  for (k in names(old_c)) if (!identical(sign(c_[[k]][[2L]])$key, sign(old_c[[k]][[2L]])$key)) return(rv)
  rat <- integer()
  for (k in names(old_c)) {
    r <- .subs_ndiv(c_[[k]][[2L]], old_c[[k]][[2L]])
    if (!r) return(rv)
    rat <- c(rat, r)
  }
  cdid <- min(rat)
  margs <- list()
  for (k in names(c_)) {
    b <- c_[[k]][[1L]]; ce <- c_[[k]][[2L]]
    margs[[length(margs) + 1L]] <- if (k %in% names(old_c)) .subs_rejoin(b, .sub(ce, .mul(list(old_c[[k]][[2L]], .int(cdid)))))
                                   else .subs_rejoin(.subs1(b, old, new), ce)
  }
  margs <- c(list(.pow(new, .int(cdid))), margs)
  .mul(list(co_residual, .mul(margs)))
}

# ---------------------------------------------------------------------------
# Pow and exp

.subs_has_symbol <- function(e) {
  if (e$type %in% c("Symbol", "Dummy", "Wild")) return(TRUE)
  for (a in e$args) if (.subs_has_symbol(a)) return(TRUE)
  FALSE
}
# Expr.as_independent(Symbol, as_Add = False)
.subs_indep_symbol <- function(e) {
  if (e$type == "Mul") {
    a <- .sympy_args(e)
    dep <- vapply(a, .subs_has_symbol, NA)
    return(list(.mul(a[!dep]), .mul(a[dep])))
  }
  if (.subs_has_symbol(e)) list(S_One, e) else list(e, S_One)
}
.subs_terms_key <- function(t) if (is_sym(t)) t$key else paste(vapply(t, .k, ""), collapse = "\u001f")
.subs_int_valued <- function(v) {
  if (v$type == "Integer") return(TRUE)
  if (v$type == "Float") { d <- num_value_dbl(v); return(is.finite(d) && d == round(d)) }
  FALSE
}
# Pow._eval_subs._check for commutative old
.subs_pow_check <- function(ct1, ct2, old_parts) {
  if (identical(.subs_terms_key(ct1[[2L]]), .subs_terms_key(ct2[[2L]]))) {
    pow <- .div(ct1[[1L]], ct2[[1L]])
    combines <- .subs_int_valued(pow)
    if (!combines) {
      b <- old_parts[[1L]]; e <- old_parts[[2L]]
      combines <- (isTRUE(ask_fact(b, "positive")) && isTRUE(ask_fact(e, "real"))) ||
                  (isTRUE(ask_fact(b, "nonnegative")) && isTRUE(ask_fact(e, "nonnegative")))
    }
    return(list(ok = combines, pow = pow))
  }
  list(ok = FALSE)
}

.pow_eval_subs <- function(base, ex, old, new, old_parts, old_is_exp) {
  if (identical(old$key, base$key)) return(.pow(new, .subs1(ex, old, new)))
  if (!is.null(old_parts) && identical(ex$key, old_parts[[2L]]$key)) {
    l <- fcall("log", list(base, old_parts[[1L]]))
    if (is_number_type(l)) return(.pow(new, l))
  }
  if (!is.null(old_parts) && identical(base$key, old_parts[[1L]]$key)) {
    if (ex$type != "Add") {
      chk <- .subs_pow_check(.subs_indep_symbol(ex), .subs_indep_symbol(old_parts[[2L]]), old_parts)
      if (chk$ok) return(.pow(new, chk$pow))
    } else {
      new_l <- list(); o_al <- list()
      ct2 <- .subs_coeff_mul(old_parts[[2L]])
      for (a in .sympy_args(ex)) {
        newa <- .subs1(a, old, new)
        chk <- .subs_pow_check(.subs_coeff_mul(newa), ct2, old_parts)
        if (chk$ok) { new_l[[length(new_l) + 1L]] <- .pow(new, chk$pow); next }
        o_al[[length(o_al) + 1L]] <- newa
      }
      if (length(new_l)) {
        expo <- .add(o_al)
        new_l[[length(new_l) + 1L]] <- if (!is_one_num(expo)) .pow(base, expo) else base
        return(.mul(new_l))
      }
    }
  }
  if (old_is_exp && isTRUE(ask_fact(ex, "extended_real")) && isTRUE(ask_fact(base, "positive"))) {
    op <- .subs_base_exp(old)
    ct1 <- .subs_indep_symbol(op[[2L]])
    ct2 <- .subs_indep_symbol(.mul(list(ex, fcall("log", list(base)))))
    chk <- .subs_pow_check(ct1, ct2, op)
    if (chk$ok) return(.pow(new, chk$pow))
  }
  NULL
}

# exp._eval_subs: power-like substitutions go through Pow._eval_subs with
# exp(a) seen as E**a
.exp_eval_subs <- function(e, old, new) {
  if (old$type == "Pow") old <- fcall("exp", list(.mul(list(old$args[[2L]], fcall("log", list(old$args[[1L]]))))))
  is_exp <- old$type == "Function" && old$name == "exp"
  if (is_exp || identical(old$key, S_E$key)) {
    parts <- if (is_exp) list(S_E, old$args[[1L]]) else NULL
    oldp <- if (is_exp) new_sym("Pow", list(S_E, old$args[[1L]])) else old
    return(.pow_eval_subs(S_E, e$args[[1L]], oldp, new, parts, is_exp))
  }
  NULL
}

# ---------------------------------------------------------------------------
# Derivative._eval_subs

.diff_wrt_ok <- function(v) v$type %in% c("Symbol", "Dummy", "Quantity") || .is_undefined_application(v) || v$type == "Derivative" ||
  (v$type == "Function" && v$name == "Indexed")
.derivative_eval_subs <- function(e, old, new) {
  expr <- e$args[[1L]]; vc <- e$args[-1L]
  wrt <- lapply(vc, function(t) t$args[[1L]])
  if (any(vapply(wrt, function(v) identical(v$key, old$key), NA))) {
    nvc <- lapply(vc, function(t) Tuple(t$args[[1L]], .subs1(t$args[[2L]], old, new)))
    if (!identical(vapply(nvc, .k, ""), vapply(vc, .k, ""))) {
      d2 <- new_sym("Derivative", c(list(expr), nvc))
      r <- .derivative_eval_subs(d2, old, new)
      return(if (is.null(r)) d2 else r)
    }
    if (!.diff_wrt_ok(new)) {
      if (old$type %in% c("Symbol", "Dummy")) return(Subs(e, old, new))
      xi <- Dummy("xi")
      return(Subs(xreplace(e, list(list(old, xi))), xi, new))
    }
  }
  if (old$type == "Derivative" && identical(old$args[[1L]]$key, expr$key)) {
    if (identical(e$key, old$key)) return(new)
    cnt <- function(d) { out <- list(); for (t in d$args[-1L]) { k <- t$args[[1L]]$key; out[[k]] <- list(t$args[[1L]], .add(list(if (is.null(out[[k]])) S_Zero else out[[k]][[2L]], t$args[[2L]]))) }; out }
    ov <- cnt(old); sv <- cnt(e)
    if (all(names(ov) %in% names(sv)) && all(vapply(names(ov), function(k) .is_true(Le(ov[[k]][[2L]], sv[[k]][[2L]])), NA))) {
      rest <- lapply(names(sv), function(k) { d <- if (is.null(ov[[k]])) sv[[k]][[2L]] else .sub(sv[[k]][[2L]], ov[[k]][[2L]]); if (is_zero_num(d)) NULL else list(sv[[k]][[1L]], d) })
      rest <- Filter(Negate(is.null), rest)
      if (!length(rest)) return(new)
      return(do.call(Derivative, c(list(new), lapply(rest, function(r) list(r[[1L]], r[[2L]])))))
    }
  }
  nexpr <- .subs1(expr, old, new)
  nvc <- lapply(vc, function(t) Tuple(.subs1(t$args[[1L]], old, new), .subs1(t$args[[2L]], old, new)))
  if (identical(expr$key, old$key)) return(.derivative_rebuild(nexpr, nvc))
  if (!identical(nexpr$key, expr$key)) {
    wk <- vapply(wrt, .k, "")
    forbidden <- intersect(names(.free_sym_list(expr)), wk)
    nfree <- setdiff(names(.free_sym_list(new)), names(.free_sym_list(old)))
    if (length(intersect(nfree, forbidden))) return(Subs(e, old, new))
  }
  changed <- !identical(vapply(nvc, function(t) t$args[[1L]]$key, ""), vapply(vc, function(t) t$args[[1L]]$key, ""))
  if (changed) {
    at <- .uniq_syms(Filter(function(z) !length(z$args) || .is_undefined_application(z), .all_subexprs(expr)))
    for (i in seq_along(nvc)) {
      vi <- nvc[[i]]$args[[1L]]
      if (!identical(vi$key, vc[[i]]$args[[1L]]$key) && any(vapply(at, function(a) identical(a$key, vi$key), NA))) return(Subs(e, old, new))
    }
    newe <- expr; subsl <- list()
    for (i in seq_along(nvc)) {
      vi <- nvc[[i]]$args[[1L]]
      if (!.diff_wrt_ok(vi)) {
        xi <- Dummy(paste0("xi_", i - 1L))
        newe <- xreplace(newe, list(list(wrt[[i]], xi)))
        nvc[[i]] <- Tuple(xi, nvc[[i]]$args[[2L]])
        subsl[[length(subsl) + 1L]] <- list(xi, vi)
      }
    }
    if (length(subsl)) {
      newe <- .subs1(newe, old, new)
      d <- .derivative_rebuild(newe, nvc)
      if (length(subsl) == 1L) return(Subs(d, subsl[[1L]][[1L]], subsl[[1L]][[2L]]))
      return(Subs(d, do.call(Tuple, lapply(subsl, `[[`, 1L)), do.call(Tuple, lapply(subsl, `[[`, 2L))))
    }
  }
  .derivative_rebuild(nexpr, nvc)
}
# _derivative_dispatch: Derivative(expr, *variable_count), unevaluated, with
# zero counts dropped
.derivative_rebuild <- function(expr, vc) {
  vc <- Filter(function(t) !is_zero_num(t$args[[2L]]), vc)
  if (!length(vc)) return(expr)
  new_sym("Derivative", c(list(expr), vc))
}
register_subs_hook("Derivative", .derivative_eval_subs)

.subs_state <- new.env()
