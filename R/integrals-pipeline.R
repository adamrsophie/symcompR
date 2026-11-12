# SymPy's indefinite integration pipeline: Integral._eval_integral of
# sympy/integrals/integrals.py with trigintegrate (trigonometry.py) and
# deltaintegrate (deltafunctions.py).

.int_is_integer_inst <- function(v) v$type == "Integer"

# trigintegrate
trigintegrate <- function(f, x, conds = "piecewise") {
  f <- as_sym(f); x <- as_sym(x)
  a <- Wild("a", exclude = list(x))
  n <- Wild("n", exclude = list(x), properties = list(.int_is_integer_inst))
  m <- Wild("m", exclude = list(x), properties = list(.int_is_integer_inst))
  f <- .rewrite_sincos(f)
  M <- tryCatch(sym_match(f, sin(a * x)^n * cos(a * x)^m), error = function(e) NULL)
  if (is.null(M) || !length(M)) return(NULL)
  nn <- M$n; mm <- M$m; aa <- M$a
  if (is.null(nn)) nn <- S_Zero
  if (is.null(mm)) mm <- S_Zero
  if (is_zero_num(nn) && is_zero_num(mm)) return(x)
  zz <- if (is_zero_num(nn)) x else S_Zero
  ni <- .int_val(nn); mi <- .int_val(mm)
  odd <- function(k) k %% 2 != 0
  if (odd(ni) || odd(mi)) {
    u <- .int_u()
    n_ <- odd(ni); m_ <- odd(mi)
    if (n_ && m_) {
      if (ni < 0 && mi > 0) { m_ <- TRUE; n_ <- FALSE }
      else if (mi < 0 && ni > 0) { n_ <- TRUE; m_ <- FALSE }
      else if (ni < 0 && mi < 0) { n_ <- ni > mi; m_ <- !(ni > mi) }
      else { n_ <- ni < mi; m_ <- !(ni < mi) }
    }
    if (n_) { ff <- .neg(.mul(list(.pow(.sub(S_One, .sq(u)), Rational(ni - 1, 2)), .pow(u, mm)))); uu <- cos(.mul(list(aa, x))) }
    else { ff <- .mul(list(.pow(u, nn), .pow(.sub(S_One, .sq(u)), Rational(mi - 1, 2)))); uu <- sin(.mul(list(aa, x))) }
    fi <- integrate(ff, u)
    fx <- subs(fi, u, uu)
    if (conds == "piecewise") return(Piecewise(list(.div(fx, aa), Ne(aa, 0)), list(zz, S_true)))
    return(.div(fx, aa))
  }
  n_ <- abs(ni) > abs(mi); m_ <- abs(mi) > abs(ni)
  res <- S_Zero
  sx <- sin(x); cx <- cos(x)
  if (n_) {
    if (mi > 0) {
      for (i in 0:(mi %/% 2)) res <- .add(list(res, .mul(list(.int((-1)^i), binomial(.int(mi %/% 2), .int(i)), .sin_pow_integrate(ni + 2 * i, x)))))
    } else if (mi == 0) res <- .sin_pow_integrate(ni, x)
    else res <- .add(list(.mul(list(Rational(-1, mi + 1), .pow(cx, .int(mi + 1)), .pow(sx, .int(ni - 1)))),
                          .mul(list(Rational(ni - 1, mi + 1), trigintegrate(.mul(list(.pow(cx, .int(mi + 2)), .pow(sx, .int(ni - 2)))), x)))))
  } else if (m_) {
    if (ni > 0) {
      for (i in 0:(ni %/% 2)) res <- .add(list(res, .mul(list(.int((-1)^i), binomial(.int(ni %/% 2), .int(i)), .cos_pow_integrate(mi + 2 * i, x)))))
    } else if (ni == 0) res <- .cos_pow_integrate(mi, x)
    else res <- .add(list(.mul(list(Rational(1, ni + 1), .pow(cx, .int(mi - 1)), .pow(sx, .int(ni + 1)))),
                          .mul(list(Rational(mi - 1, ni + 1), trigintegrate(.mul(list(.pow(cx, .int(mi - 2)), .pow(sx, .int(ni + 2)))), x)))))
  } else if (mi == ni) {
    res <- integrate(.pow(.mul(list(sin(.mul(list(S_Two, x))), S_Half)), mm), x)
  } else if (mi == -ni) {
    if (ni < 0) res <- .add(list(.mul(list(Rational(1, ni + 1), .pow(cx, .int(mi - 1)), .pow(sx, .int(ni + 1)))),
                                 .mul(list(Rational(mi - 1, ni + 1), integrate(.mul(list(.pow(cx, .int(mi - 2)), .pow(sx, .int(ni + 2)))), x)))))
    else res <- .add(list(.mul(list(Rational(-1, mi + 1), .pow(cx, .int(mi + 1)), .pow(sx, .int(ni - 1)))),
                          .mul(list(Rational(ni - 1, mi + 1), integrate(.mul(list(.pow(cx, .int(mi + 2)), .pow(sx, .int(ni - 2)))), x)))))
  }
  if (conds == "piecewise") return(Piecewise(list(.div(subs(res, x, .mul(list(aa, x))), aa), Ne(aa, 0)), list(zz, S_true)))
  .div(subs(res, x, .mul(list(aa, x))), aa)
}
.int_u <- function() { if (is.null(.mi_state$int_u)) .mi_state$int_u <- Dummy("u"); .mi_state$int_u }
.rewrite_sincos <- function(f) .map_bottom_up(f, function(z) {
  if (z$type != "Function") return(z)
  a <- z$args[[1L]]
  switch(z$name,
    tan = .div(sin(a), cos(a)), cot = .div(cos(a), sin(a)), sec = .inv(cos(a)), csc = .inv(sin(a)), z)
})

.sin_pow_integrate <- function(n, x) {
  if (n > 0) {
    if (n == 1) return(.neg(cos(x)))
    return(.add(list(.mul(list(Rational(-1, n), cos(x), .pow(sin(x), .int(n - 1)))), .mul(list(Rational(n - 1, n), .sin_pow_integrate(n - 2, x))))))
  }
  if (n < 0) {
    if (n == -1) return(trigintegrate(.inv(sin(x)), x))
    return(.add(list(.mul(list(Rational(1, n + 1), cos(x), .pow(sin(x), .int(n + 1)))), .mul(list(Rational(n + 2, n + 1), .sin_pow_integrate(n + 2, x))))))
  }
  x
}
.cos_pow_integrate <- function(n, x) {
  if (n > 0) {
    if (n == 1) return(sin(x))
    return(.add(list(.mul(list(Rational(1, n), sin(x), .pow(cos(x), .int(n - 1)))), .mul(list(Rational(n - 1, n), .cos_pow_integrate(n - 2, x))))))
  }
  if (n < 0) {
    if (n == -1) return(trigintegrate(.inv(cos(x)), x))
    return(.add(list(.mul(list(Rational(-1, n + 1), sin(x), .pow(cos(x), .int(n + 1)))), .mul(list(Rational(n + 2, n + 1), .cos_pow_integrate(n + 2, x))))))
  }
  x
}

# ---------------------------------------------------------------------------
# deltaintegrate

.is_dirac <- function(e) e$type == "Function" && e$name == "DiracDelta"
.dirac_is_simple <- function(d, x) { p <- P_or_null(d$args[[1L]], x); !is.null(p) && P_deg(p) == 1L }
# DiracDelta._eval_expand_diracdelta(wrt = x)
.dirac_expand <- function(d, x) {
  if (!.has1(d$args[[1L]], x) || (length(d$args) > 1L && !is_zero_num(d$args[[2L]]))) return(d)
  rs <- tryCatch(roots(d$args[[1L]], x), error = function(e) NULL)
  if (is.null(rs)) return(d)
  darg <- Abs(.diff1(d$args[[1L]], x))
  mult <- attr(rs, "multiplicity")
  res <- list()
  for (k in seq_along(rs)) {
    r <- rs[[k]]
    if (!identical(ask_fact(r, "real"), FALSE) && mult[k] == 1L) res[[length(res) + 1L]] <- .div(fcall("DiracDelta", list(.sub(x, r))), subs(darg, x, r))
    else return(d)
  }
  .add(res)
}
.expand_dirac_in <- function(e, x) .map_bottom_up(e, function(z) if (.is_dirac(z)) .dirac_expand(z, x) else z)

.change_mul <- function(node, x) {
  args <- if (node$type == "Mul") node$args else list(node)
  comm <- vapply(args, .is_comm, NA)
  c <- args[comm]; nc <- args[!comm]
  sorted_args <- c(.key_sort(lapply(c, .sort_key), c), nc)
  new_args <- list(); dirac <- NULL
  for (arg in sorted_args) {
    if (arg$type == "Pow" && .is_dirac(arg$args[[1L]])) {
      new_args[[length(new_args) + 1L]] <- .pow(arg$args[[1L]], .sub(arg$args[[2L]], S_One))
      arg <- arg$args[[1L]]
    }
    if (is.null(dirac) && .is_dirac(arg) && .dirac_is_simple(arg, x)) dirac <- arg
    else new_args[[length(new_args) + 1L]] <- arg
  }
  if (is.null(dirac)) {
    na <- lapply(sorted_args, function(arg) {
      if (.is_dirac(arg)) return(.dirac_expand(arg, x))
      if (arg$type == "Pow" && .is_dirac(arg$args[[1L]])) return(.pow(.dirac_expand(arg$args[[1L]], x), arg$args[[2L]]))
      arg
    })
    same <- length(na) == length(sorted_args) && all(vapply(seq_along(na), function(i) identical(na[[i]]$key, sorted_args[[i]]$key), NA))
    return(list(NULL, if (same) NULL else expand(.mul(na))))
  }
  list(dirac, .mul(new_args))
}

deltaintegrate <- function(f, x) {
  f <- as_sym(f); x <- as_sym(x)
  if (!.has_function(f, "DiracDelta")) return(NULL)
  if (.is_dirac(f)) {
    h <- .dirac_expand(f, x)
    if (identical(h$key, f$key)) {
      if (.dirac_is_simple(f, x)) {
        if (length(f$args) <= 1L || is_zero_num(f$args[[2L]])) return(Heaviside(f$args[[1L]]))
        return(.div(fcall("DiracDelta", list(f$args[[1L]], .sub(f$args[[2L]], S_One))), P_LC(.P(f$args[[1L]], x))))
      }
    } else return(integrate(h, x))
  } else if (f$type %in% c("Mul", "Pow")) {
    g <- expand(f)
    if (!identical(f$key, g$key)) {
      fh <- integrate(g, x)
      if (!is.null(fh) && fh$type != "Integral") return(fh)
    } else {
      cm <- .change_mul(f, x)
      deltaterm <- cm[[1L]]; rest_mult <- cm[[2L]]
      if (is.null(deltaterm)) {
        if (!is.null(rest_mult)) return(integrate(rest_mult, x))
      } else {
        deltaterm <- .dirac_expand(deltaterm, x)
        if (deltaterm$type == "Mul") {
          cm2 <- .change_mul(deltaterm, x)
          deltaterm <- cm2[[1L]]; rest_mult <- .mul(list(rest_mult, cm2[[2L]]))
        }
        point <- solve(deltaterm$args[[1L]], x)[[1L]]
        n <- if (length(deltaterm$args) == 1L) 0L else .int_val(deltaterm$args[[2L]])
        m <- 0L
        while (n >= 0) {
          r <- .mul(list(.pow(S_NegativeOne, .int(n)), subs(if (n) diff(rest_mult, x, n) else rest_mult, x, point)))
          if (isTRUE(ask_fact(r, "zero"))) { n <- n - 1L; m <- m + 1L }
          else if (m == 0L) return(.mul(list(r, Heaviside(.sub(x, point)))))
          else return(.mul(list(r, fcall("DiracDelta", list(x, .int(m - 1L))))))
        }
        return(S_Zero)
      }
    }
  }
  NULL
}

# ---------------------------------------------------------------------------
# Integral._eval_integral

.sympy_eval_integral <- function(f, x, meijerg = NULL, risch = NULL, manual = NULL, heurisch = NULL, conds = "piecewise") {
  if (isTRUE(risch)) return(tryCatch(risch_integrate(f, x, conds = conds), NotImplementedError = function(e) NULL))
  if (isTRUE(manual)) {
    r <- tryCatch(manualintegrate(f, x), error = function(e) NULL)
    if (!is.null(r) && r$type != "Integral") return(r)
  }
  kw <- list(meijerg = meijerg, risch = risch, manual = manual, heurisch = heurisch, conds = conds)
  if (f$type == "Piecewise")
    return(do.call(Piecewise, lapply(f$args, function(p) list(do.call(integrate, c(list(p$args[[1L]], x), kw)), p$args[[2L]]))))
  if (!.has1(f, x)) return(.mul(list(f, x)))
  plain <- !isTRUE(manual) && !isTRUE(meijerg) && !isTRUE(risch)
  if (plain) {
    p <- P_or_null(f, x)
    if (!is.null(p)) return(P_expr(.dp_integrate(p)))
  }
  if (!identical(risch, FALSE)) {
    r <- tryCatch(risch_integrate(f, x, separate_integral = TRUE, conds = conds), NotImplementedError = function(e) NULL)
    if (!is.null(r)) {
      result <- r[[1L]]; i <- r[[2L]]
      if (!is_zero_num(i)) {
        if (is_zero_num(result)) return(.int_doit_norisch(NonElementaryIntegral(f, x), meijerg, manual, heurisch, conds))
        return(.add(list(result, .int_doit_norisch(i, meijerg, manual, heurisch, conds))))
      }
      return(result)
    }
  }
  parts <- list()
  args <- if (f$type == "Add") .sympy_args(f) else list(f)
  for (g in args) {
    cm <- .as_coeff_mul_dep(g, x)
    coeff <- cm[[1L]]; g <- .mul(cm[[2L]])
    if (is_one_num(g) && !isTRUE(meijerg)) { parts[[length(parts) + 1L]] <- .mul(list(coeff, x)); next }
    # g(x) = expr + O(x**n)
    order_term <- .so_getO(g)
    if (!is.null(order_term)) {
      h <- .sympy_eval_integral(.so_removeO(g), x, meijerg, risch, manual, heurisch, conds)
      if (!is.null(h)) {
        ho <- .sympy_eval_integral(order_term$args[[1L]], x, meijerg, risch, manual, heurisch, conds)
        if (!is.null(ho)) {
          parts[[length(parts) + 1L]] <- .mul(list(coeff, .add(list(h, .so_new(ho, .so_vars(order_term))))))
          next
        }
      }
      return(NULL)
    }
    if (g$type == "Pow" && !.has1(g$args[[2L]], x) && !isTRUE(meijerg)) {
      a <- Wild("a", exclude = list(x)); b <- Wild("b", exclude = list(x))
      M <- tryCatch(sym_match(g$args[[1L]], a * x + b), error = function(e) NULL)
      if (!is.null(M) && length(M) && !is.null(M$a)) {
        base <- g$args[[1L]]; ex <- g$args[[2L]]
        h <- if (identical(ex$key, S_NegativeOne$key)) fcall("log", list(base))
             else if (conds != "piecewise") .div(.pow(base, .add(list(ex, S_One))), .add(list(ex, S_One)))
             else Piecewise(list(.div(.pow(base, .add(list(ex, S_One))), .add(list(ex, S_One))), Ne(ex, -1)), list(fcall("log", list(base)), S_true))
        parts[[length(parts) + 1L]] <- .div(.mul(list(coeff, h)), M$a)
        next
      }
    }
    if (plain && .is_ratfun(g, x)) { parts[[length(parts) + 1L]] <- .mul(list(coeff, ratint(g, x))); next }
    h <- NULL
    if (plain) {
      h <- trigintegrate(g, x, conds = conds)
      if (!is.null(h)) { parts[[length(parts) + 1L]] <- .mul(list(coeff, h)); next }
      h <- deltaintegrate(g, x)
      if (!is.null(h)) { parts[[length(parts) + 1L]] <- .mul(list(coeff, h)); next }
      h <- singularityintegrate(g, x)
      if (!is.null(h)) { parts[[length(parts) + 1L]] <- .mul(list(coeff, h)); next }
      if (!identical(risch, FALSE)) {
        r <- tryCatch(risch_integrate(g, x, separate_integral = TRUE, conds = conds), NotImplementedError = function(e) NULL)
        if (!is.null(r)) {
          h <- r[[1L]]
          if (!is_zero_num(r[[2L]])) h <- .add(list(h, .int_doit_norisch(r[[2L]], meijerg, manual, heurisch, conds)))
          parts[[length(parts) + 1L]] <- .mul(list(coeff, h))
          next
        }
      }
      if (!identical(heurisch, FALSE))
        h <- tryCatch(if (conds == "piecewise") heurisch_wrapper(g, x, hints = list()) else heurisch(g, x, hints = list()), error = function(e) NULL)
    }
    if (!identical(meijerg, FALSE) && is.null(h)) {
      h <- tryCatch(meijerint_indefinite(g, x), NotImplementedError = function(e) NULL)
      if (!is.null(h)) { parts[[length(parts) + 1L]] <- .mul(list(coeff, h)); next }
    }
    if (is.null(h) && !identical(manual, FALSE)) {
      result <- tryCatch(manualintegrate(g, x), error = function(e) NULL)
      if (!is.null(result) && result$type != "Integral") {
        if (.has_type(result, "Integral") && !isTRUE(manual)) {
          kw2 <- kw; kw2$manual <- FALSE
          result <- expand(sym_rebuild(result, lapply(result$args, function(a) if (.has_type(a, "Integral")) .int_doit_all(a, kw2) else a)),
                           multinomial = FALSE, log = FALSE, power_exp = FALSE, power_base = FALSE)
        }
        if (!.has_type(result, "Integral")) { parts[[length(parts) + 1L]] <- .mul(list(coeff, result)); next }
      }
    }
    if ((is.null(h) || is_zero_num(h)) && length(args) == 1L) {
      ff <- expand(sincos_to_sum(f), mul = TRUE, deep = FALSE)
      if (ff$type == "Add") return(do.call(.sympy_eval_integral, c(list(ff, x), kw)))
    }
    if (!is.null(h)) parts[[length(parts) + 1L]] <- .mul(list(coeff, h))
    else return(NULL)
  }
  .add(parts)
}

# expr.doit(**kwargs) for the Integral nodes of an expression
.int_doit_all <- function(e, kw) .map_bottom_up(e, function(z) {
  if (z$type != "Integral" || length(z$args) != 2L) return(z)
  lim <- z$args[[2L]]
  if (lim$type == "Tuple" && length(lim$args) != 1L) return(z)
  v <- if (lim$type == "Tuple") lim$args[[1L]] else lim
  r <- tryCatch(do.call(.sympy_eval_integral, c(list(z$args[[1L]], v), kw)), error = function(err) NULL)
  if (is.null(r)) z else r
})

# NonElementaryIntegral(...).doit(risch=False)
.int_doit_norisch <- function(i, meijerg, manual, heurisch, conds) {
  f <- i$args[[1L]]; lim <- i$args[[2L]]
  x <- if (lim$type == "Tuple") lim$args[[1L]] else lim
  r <- .sympy_eval_integral(f, x, meijerg = meijerg, risch = FALSE, manual = manual, heurisch = heurisch, conds = conds)
  if (is.null(r)) i else r
}

