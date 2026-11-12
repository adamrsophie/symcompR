# Integral.doit after sympy/integrals/integrals.py, with Expr._eval_interval
# (sympy/core/expr.py), Piecewise._eval_integral, Piecewise._eval_interval,
# Piecewise._handle_irel and _clip (sympy/functions/elementary/piecewise.py),
# _add_atan_floor_terms and _atan_floor_correction.

.it_illegal <- function(C) .has1(C, S_NaN) || .has1(C, S_Infinity) || .has1(C, S_NegativeInfinity) ||
  .has1(C, S_ComplexInfinity) || .has_type(C, "AccumulationBounds")

# Expr._eval_interval; NULL bounds stand for None
.expr_eval_interval <- function(F, x, a, b) {
  if (is.null(a) && is.null(b)) stop("Both interval ends cannot be None.", call. = FALSE)
  endpoint <- function(left) {
    c <- if (left) a else b
    if (is.null(c)) return(S_Zero)
    C <- subs(F, x, c)
    if (.it_illegal(C)) {
      fwd <- is.null(a) || is.null(b) || !identical(Lt(a, b)$key, S_false$key)
      dir <- if (fwd) (if (left) "+" else "-") else (if (left) "-" else "+")
      C <- limit(F, x, c, dir)
      if (C$type == "Limit") .notimpl_error("Could not compute limit")
    }
    C
  }
  if (!is.null(a) && !is.null(b) && identical(a$key, b$key)) return(S_Zero)
  A <- endpoint(TRUE)
  if (is_nan(A)) return(A)
  B <- endpoint(FALSE)
  if (is.null(a) || (!is_zero_num(a) && is.null(b))) return(.sub(B, A))
  value <- .sub(B, A)
  if (.is_comparable(a) && .is_comparable(b)) {
    domain <- if (.is_true(Lt(a, b))) Interval(a, b) else Interval(b, a)
    sing <- solveset(.as_numer_denom(cancel(F))[[2L]], x, domain = domain)
    for (lt in .uniq_syms(Filter(function(z) z$type == "Function" && z$name == "log", .all_subexprs(F))))
      sing <- Union(sing, solveset(lt$args[[1L]], x, domain = domain))
    if (sing$type == "FiniteSet") {
      for (s in sing$args) {
        if (is_nan(value)) break
        if (!.is_comparable(s)) next
        if (.is_true(Lt(a, s)) && .is_true(Lt(s, b)))
          value <- .add(list(value, .neg(limit(F, x, s, "+")), limit(F, x, s, "-")))
        else if (.is_true(Lt(b, s)) && .is_true(Lt(s, a)))
          value <- .add(list(value, limit(F, x, s, "+"), .neg(limit(F, x, s, "-"))))
      }
    }
  }
  value
}

.eval_interval <- function(F, x, a, b) if (F$type == "Piecewise") .pw_eval_interval(F, x, a, b) else .expr_eval_interval(F, x, a, b)

# _clip(A, B, k)
.pw_clip <- function(A, B, k) {
  a <- B[[1L]]; b <- B[[2L]]; c <- A[[1L]]; d <- A[[2L]]
  c <- Min(Max(c, a), b); d <- Min(Max(d, a), b)
  a <- Min(a, b)
  p <- list()
  if (!.mi_eq(a, c)) p[[length(p) + 1L]] <- list(a, c, -1L)
  if (!.mi_eq(c, d)) p[[length(p) + 1L]] <- list(c, d, k)
  if (!.mi_eq(b, d)) {
    if (.mi_eq(d, c) && length(p) && p[[length(p)]][[3L]] == -1L) p[[length(p)]] <- list(p[[length(p)]][[1L]], b, -1L)
    else p[[length(p) + 1L]] <- list(d, b, -1L)
  }
  p
}
# the replacement of done[j] by the clipped pieces of each interval
.pw_cover <- function(abei, lo, hi) {
  done <- list(list(lo, hi, -1L))
  for (k in seq_along(abei)) {
    p <- abei[[k]]
    if (identical(p[[1L]]$key, S_NegativeInfinity$key) && identical(p[[2L]]$key, S_Infinity$key)) {
      for (j in seq_along(done)) if (done[[j]][[3L]] == -1L) done[[j]][[3L]] <- k
      break
    }
    for (j in rev(seq_along(done))) {
      if (done[[j]][[3L]] == -1L) {
        cl <- .pw_clip(list(p[[1L]], p[[2L]]), list(done[[j]][[1L]], done[[j]][[2L]]), k)
        done <- c(done[seq_len(j - 1L)], cl, done[seq_along(done) > j])
      }
    }
  }
  Filter(function(d) !.mi_eq(d[[1L]], d[[2L]]), done)
}

# Piecewise._handle_irel: the relationals free of x, decided in turn
.pw_handle_irel <- function(pw, x, handler) {
  rel <- Filter(function(r) r$type %in% .REL_TYPES, .uniq_syms(.all_subexprs(pw)))
  irel <- .mi_sort_ordered(Filter(function(r) !.has1(r, x), rel))
  if (!length(irel)) return(NULL)
  n <- length(irel)
  args <- list(); order_keys <- character(); exprs <- list()
  expr <- NULL
  for (code in seq_len(2L^n) - 1L) {
    truth <- rev(as.integer(intToBits(2L^n - 1L - code))[seq_len(n)])
    reps <- lapply(seq_len(n), function(i) list(irel[[i]], if (truth[i]) S_true else S_false))
    if (!any(truth == 1L)) cond <- NULL
    else {
      andargs <- irel[truth == 1L]
      free <- .uniq_syms(unlist(lapply(andargs, free_symbols), recursive = FALSE))
      if (length(free) == 1L) {
        t <- tryCatch(reduce_inequalities(andargs, free[[1L]]), error = function(e) NULL)
        if (is.null(t)) t <- do.call(And, lapply(andargs, function(a) .solve_inequality(a, free[[1L]], linear = TRUE)))
      } else t <- if (length(andargs) == 1L) andargs[[1L]] else do.call(And, andargs)
      if (identical(t$key, S_false$key)) next
      cond <- t
    }
    expr <- handler(xreplace(pw, reps))
    if (expr$type == "Piecewise" && length(expr$args) == 1L) {
      ec <- expr$args[[1L]]$args[[2L]]; expr <- expr$args[[1L]]$args[[1L]]
      cond <- And(ec, if (is.null(cond)) S_true else cond)
    }
    if (!is.null(cond)) {
      k <- expr$key
      if (is.null(args[[k]])) { args[[k]] <- list(); exprs[[k]] <- expr }
      args[[k]][[length(args[[k]]) + 1L]] <- cond
      order_keys <- c(order_keys, k)
    }
  }
  pairs <- lapply(unique(order_keys), function(k) list(exprs[[k]], if (length(args[[k]]) == 1L) args[[k]][[1L]] else do.call(Or, args[[k]])))
  pairs[[length(pairs) + 1L]] <- list(expr, S_true)
  do.call(Piecewise, pairs)
}

# Piecewise._eval_integral (indefinite)
.pw_eval_integral <- function(pw, x, hints = list(), first = TRUE) {
  if (first) {
    irv <- .pw_handle_irel(pw, x, function(ipw) {
      if (ipw$type == "Piecewise") .pw_eval_integral(ipw, x, hints, first = FALSE)
      else do.call(integrate, c(list(ipw, x), hints))
    })
    if (!is.null(irv)) return(irv)
  }
  iv <- .pw_intervals(pw, x)
  if (!isTRUE(iv[[1L]])) return(Integral(pw, x))
  abei <- iv[[2L]]
  done <- .pw_cover(abei, S_NegativeInfinity, S_Infinity)
  if (any(vapply(done, function(d) d[[3L]] == -1L, NA))) {
    abei[[length(abei) + 1L]] <- list(S_NegativeInfinity, S_Infinity, S_NaN, -1L)
    done <- lapply(done, function(d) { if (d[[3L]] == -1L) d[[3L]] <- length(abei); d })
  }
  out <- list(); sum <- NULL
  for (d in done) {
    a <- d[[1L]]; b <- d[[2L]]; i <- d[[3L]]
    anti <- do.call(integrate, c(list(abei[[i]][[3L]], x), hints))
    if (is.null(sum)) sum <- anti
    else {
      sum <- subs(sum, x, a)
      e <- .eval_interval(anti, x, a, x)
      if (.it_illegal(sum) || .it_illegal(e)) sum <- anti else sum <- .add(list(sum, e))
    }
    idx <- abei[[i]][[4L]]
    cond <- if (identical(b$key, S_Infinity$key)) S_true
            else if (idx >= 1L && identical(subs(pw$args[[idx]]$args[[2L]], x, b)$key, S_false$key)) Lt(x, b)
            else Le(x, b)
    out[[length(out) + 1L]] <- list(sum, cond)
  }
  do.call(Piecewise, out)
}

# Piecewise._eval_interval
.pw_eval_interval <- function(pw, x, lo, hi, first = TRUE) {
  if (is.null(lo) || is.null(hi)) return(.expr_eval_interval(pw, x, lo, hi))
  if (isTRUE(first)) {
    irv <- .pw_handle_irel(pw, x, function(ipw) {
      if (ipw$type == "Piecewise") .pw_eval_interval(ipw, x, lo, hi, first = NA) else .expr_eval_interval(ipw, x, lo, hi)
    })
    if (!is.null(irv)) return(irv)
    lt <- Lt(lo, hi)
    if (identical(lt$key, S_false$key) || identical(lo$key, S_Infinity$key) || identical(hi$key, S_NegativeInfinity$key)) {
      rv <- .pw_eval_interval(pw, x, hi, lo, first = FALSE)
      if (rv$type == "Piecewise") return(do.call(Piecewise, lapply(rv$args, function(p) list(.neg(p$args[[1L]]), p$args[[2L]]))))
      return(.neg(rv))
    }
    if (!(.is_true(lt) || identical(hi$key, S_Infinity$key) || identical(lo$key, S_NegativeInfinity$key))) {
      .a <- Dummy("lo"); .b <- Dummy("hi")
      a <- if (.is_comparable(lo)) lo else .a
      b <- if (.is_comparable(hi)) hi else .b
      pos <- .pw_eval_interval(pw, x, a, b, first = FALSE)
      if (.mi_eq(a, .a) && .mi_eq(b, .b)) {
        neg <- .neg(xreplace(pos, list(list(.a, hi), list(.b, lo))))
        pos <- xreplace(pos, list(list(.a, lo), list(.b, hi)))
      } else {
        neg <- .neg(.pw_eval_interval(pw, x, hi, lo, first = FALSE))
        pos <- xreplace(pos, list(list(.a, lo), list(.b, hi)))
      }
      p <- Dummy("", positive = TRUE)
      if (lo$type %in% c("Symbol", "Dummy")) {
        pos <- xreplace(xreplace(pos, list(list(lo, .sub(hi, p)))), list(list(p, .sub(hi, lo))))
        neg <- xreplace(xreplace(neg, list(list(lo, .add(list(hi, p))))), list(list(p, .sub(lo, hi))))
      } else if (hi$type %in% c("Symbol", "Dummy")) {
        pos <- xreplace(xreplace(pos, list(list(hi, .add(list(lo, p))))), list(list(p, .sub(hi, lo))))
        neg <- xreplace(xreplace(neg, list(list(hi, .sub(lo, p)))), list(list(p, .sub(lo, hi))))
      }
      touch <- function(z) .map_bottom_up(z, function(w) if (w$type == "Function" && w$name %in% c("Min", "Max")) do.call(get(w$name), w$args) else w)
      neg <- touch(neg); pos <- touch(pos)
      rv <- if (.mi_eq(a, .a)) Piecewise(list(pos, Lt(lo, hi)), list(neg, S_true)) else Piecewise(list(neg, Lt(hi, lo)), list(pos, S_true))
      if (is_nan(rv)) stop("Can't integrate across undefined region.", call. = FALSE)
      if (pos$type == "Piecewise" || neg$type == "Piecewise") rv <- piecewise_fold(rv)
      return(rv)
    }
  }
  iv <- .pw_intervals(pw, x)
  if (!isTRUE(iv[[1L]])) return(Integral(.diff1(pw, x), list(x, lo, hi)))
  abei <- iv[[2L]]
  done <- .pw_cover(abei, lo, hi)
  sum <- S_Zero; upto <- NULL
  for (d in done) {
    a <- d[[1L]]; b <- d[[2L]]; i <- d[[3L]]
    if (i == -1L) {
      if (is.null(upto)) return(S_NaN)
      return(Piecewise(list(sum, Le(hi, upto)), list(S_NaN, S_true)))
    }
    k <- .eval_interval(abei[[i]][[3L]], x, a, b)
    if (is_nan(k)) return(Piecewise(list(sum, Le(hi, upto)), list(S_NaN, S_true)))
    sum <- .add(list(sum, k)); upto <- b
  }
  sum
}

# ---------------------------------------------------------------------------
# atan floor terms

.atan_floor_correction <- function(atan_arg, x) {
  parts <- .mi_sort_ordered(Filter(function(z) z$type == "Function" && z$name %in% c("tan", "cot") && .has1(z, x), .uniq_syms(.all_subexprs(atan_arg))))
  if (!length(parts)) return(NULL)
  dums <- lapply(parts, function(p) Dummy("p"))
  masked <- xreplace(atan_arg, lapply(seq_along(parts), function(i) list(parts[[i]], dums[[i]])))
  me <- expand(masked)
  if (!.mi_is_poly(me, dums)) return(NULL)
  # total degree at most one: the coefficient of each part and the rest
  coefs <- lapply(dums, function(dv) coeff(me, dv))
  d <- expand(xreplace(me, lapply(dums, function(dv) list(dv, S_Zero))))
  if (any(vapply(coefs, function(cf) any(vapply(dums, function(dv) .has1(cf, dv), NA)), NA))) return(NULL)
  if (!is_zero_num(expand(.sub(me, .add(c(list(d), lapply(seq_along(dums), function(i) .mul(list(coefs[[i]], dums[[i]]))))))))) return(NULL)
  n <- Dummy("n", integer = TRUE)
  terms <- list()
  for (i in seq_along(parts)) {
    part <- parts[[i]]; cf <- coefs[[i]]; a <- part$args[[1L]]
    if (!isTRUE(ask_fact(cf, "extended_real")) || !isTRUE(ask_fact(a, "extended_real"))) return(NULL)
    if (.has1(d, x) && !isTRUE(ask_fact(d, "finite"))) {
      m <- .diff1(a, x)
      if (.has1(m, x) || !identical(ask_fact(m, "zero"), FALSE)) return(NULL)
      poles <- .div(.sub(.mul(list(n, S_Pi)), subs(a, x, S_Zero)), m)
      if (part$name == "tan") poles <- .add(list(poles, .div(S_Pi, .mul(list(S_Two, m)))))
      if (!isTRUE(ask_fact(subs(d, x, poles), "finite"))) return(NULL)
    }
    terms[[length(terms) + 1L]] <- if (part$name == "tan") list(cf, .div(.add(list(a, .div(S_Pi, S_Two))), S_Pi)) else list(.neg(cf), .div(a, S_Pi))
  }
  sgn <- function(z) fcall("sign", list(z))
  if (length(terms) == 1L) {
    cf <- terms[[1L]][[1L]]; phi <- terms[[1L]][[2L]]
    if (.has1(sgn(cf), x)) return(NULL)
    return(.mul(list(sgn(cf), S_Pi, fcall("floor", list(phi)))))
  }
  groups <- list(); overlaps <- list()
  for (tm in terms) {
    cf <- tm[[1L]]; phi <- tm[[2L]]
    m <- .diff1(phi, x)
    if (.has1(cf, x) || .has1(m, x) || !identical(ask_fact(m, "zero"), FALSE)) return(NULL)
    x0 <- .div(.neg(subs(phi, x, S_Zero)), m); period <- .inv(m)
    new <- list(coeff = cf, phi = phi, x0 = x0, period = period, id = length(groups) + 1L)
    new_ov <- list(); merged <- FALSE
    for (gi in seq_along(groups)) {
      g <- groups[[gi]]
      ratio <- .div(g$period, period)
      if (!(ratio$type %in% c("Integer", "Rational"))) return(NULL)
      rq <- .mi_q(ratio); rp <- .mi_p(ratio)
      common <- .div(period, .int(rq))
      shift <- expand(.div(.sub(x0, g$x0), common))
      si <- ask_fact(shift, "integer")
      if (is.na(si)) return(NULL)
      if (!si) next
      if (is_one_num(abs(ratio))) { groups[[gi]]$coeff <- .add(list(g$coeff, .mul(list(ratio, cf)))); merged <- TRUE; break }
      if (rq == 1L) k <- 0L
      else if (shift$type == "Integer") {
        inv <- as.integer(as.character(gmp::inv.bigz(gmp::as.bigz(rp %% rq), gmp::as.bigz(rq))))
        k <- (.int_val(shift) * inv) %% rq
      } else return(NULL)
      new_ov[[length(new_ov) + 1L]] <- list(gi, length(groups) + 1L, ratio, k)
    }
    if (!merged) { groups[[length(groups) + 1L]] <- new; overlaps <- c(overlaps, new_ov) }
  }
  result <- lapply(Filter(function(g) !isTRUE(ask_fact(g$coeff, "zero")), groups),
                   function(g) .mul(list(sgn(g$coeff), S_Pi, fcall("floor", list(g$phi)))))
  seen <- integer()
  for (ov in overlaps) {
    if (ov[[1L]] %in% seen || ov[[2L]] %in% seen) return(NULL)
    seen <- c(seen, ov[[1L]], ov[[2L]])
    g <- groups[[ov[[1L]]]]; h <- groups[[ov[[2L]]]]; ratio <- ov[[3L]]
    cg <- g$coeff; ch <- .div(h$coeff, ratio)
    result[[length(result) + 1L]] <- .mul(list(.sub(.sub(sgn(.add(list(cg, ch))), sgn(cg)), sgn(ch)), S_Pi,
                                               fcall("floor", list(.div(.sub(g$phi, .int(ov[[4L]])), .int(.mi_q(ratio)))))))
  }
  .add(result)
}

.add_atan_floor_terms <- function(antideriv, x) {
  if (antideriv$type == "Piecewise")
    return(do.call(Piecewise, lapply(antideriv$args, function(p) list(.add_atan_floor_terms(p$args[[1L]], x), p$args[[2L]]))))
  if (!isTRUE(ask_fact(x, "extended_real"))) return(antideriv)
  A <- Dummy("A")
  for (at in .uniq_syms(Filter(function(z) z$type == "Function" && z$name == "atan", .all_subexprs(antideriv)))) {
    corr <- .atan_floor_correction(at$args[[1L]], x)
    if (is.null(corr)) next
    p <- P_or_null(xreplace(antideriv, list(list(at, A))), A)
    if (is.null(p)) next
    cs <- P_all_coeffs(p)
    coeffs <- lapply(cs[-length(cs)], cancel)
    if (any(vapply(coeffs, function(c) .has1(c, x), NA))) next
    rc <- rev(coeffs)
    jump <- .div(.add(lapply(seq_along(rc), function(k) .mul(list(rc[[k]], .sub(.pow(.div(S_Pi, S_Two), .int(k)), .pow(.div(.neg(S_Pi), S_Two), .int(k))))))), S_Pi)
    if (!is_zero_num(jump)) antideriv <- .add(list(antideriv, .mul(list(jump, corr))))
  }
  antideriv
}

# ---------------------------------------------------------------------------
# Integral.doit

# AddWithLimits._eval_factor
.integral_factor <- function(I) {
  lims <- I$args[-1L]
  vars <- lapply(lims, function(l) l$args[[1L]])
  if (length(lims) == 1L) {
    summand <- sym_factor(I$args[[1L]])
    if (summand$type == "Mul") {
      ind <- vapply(summand$args, function(w) .is_comm(w) && !any(vapply(vars, function(v) .has1(w, v), NA)), NA)
      return(.mul(list(.mul(summand$args[ind]), new_sym("Integral", c(list(.mul(summand$args[!ind])), lims)))))
    }
    return(I)
  }
  summand <- .integral_factor(new_sym("Integral", c(list(I$args[[1L]]), lims[-length(lims)])))
  last <- lims[[length(lims)]]
  if (!.has1(summand, last$args[[1L]])) return(.mul(list(.integral_doit(new_sym("Integral", list(S_One, last))), summand)))
  if (summand$type == "Mul") return(.integral_factor(new_sym("Integral", list(summand, last))))
  I
}

# Integral._eval_is_zero
.integral_is_zero <- function(I) {
  f <- I$args[[1L]]; lims <- I$args[-1L]
  fz <- ask_fact(f, "zero")
  if (isTRUE(fz)) return(TRUE)
  for (l in lims) if (length(l$args) == 3L) {
    if (.mi_eq(l$args[[2L]], l$args[[3L]]) || isTRUE(ask_fact(.sub(l$args[[2L]], l$args[[3L]]), "zero"))) return(TRUE)
  }
  free <- names(.free_sym_list(f))
  for (l in lims) {
    v <- l$args[[1L]]
    if (length(l$args) == 1L) { free <- union(free, v$key); next }
    if (length(l$args) == 2L && !(v$key %in% free) && isTRUE(ask_fact(l$args[[2L]], "zero"))) return(TRUE)
    free <- setdiff(free, v$key)
    for (z in l$args[-1L]) free <- union(free, names(.free_sym_list(z)))
  }
  NA
}

.integral_doit <- function(I, hints = list(), deep = TRUE, final = TRUE) {
  meijerg <- hints$meijerg; risch <- hints$risch; heurisch <- hints$heurisch; manual <- hints$manual
  conds <- if (is.null(hints$conds)) "piecewise" else hints$conds
  if (sum(vapply(list(manual, meijerg, risch, heurisch), isTRUE, NA)) > 1L)
    stop("At most one of manual, meijerg, risch, heurisch can be True", call. = FALSE)
  if (isTRUE(manual)) { meijerg <- FALSE; risch <- FALSE; heurisch <- FALSE }
  else if (isTRUE(meijerg)) { manual <- FALSE; risch <- FALSE; heurisch <- FALSE }
  else if (isTRUE(risch)) { manual <- FALSE; meijerg <- FALSE; heurisch <- FALSE }
  else if (isTRUE(heurisch)) { manual <- FALSE; meijerg <- FALSE; risch <- FALSE }
  kw <- list(meijerg = meijerg, risch = risch, manual = manual, heurisch = heurisch, conds = conds)
  hints <- kw
  if (!(conds %in% c("separate", "piecewise", "none"))) stop("conds must be one of \"separate\", \"piecewise\", \"none\", got: ", conds, call. = FALSE)
  lims <- I$args[-1L]
  if (isTRUE(risch) && any(vapply(lims, function(l) length(l$args) > 1L, NA))) stop("risch=True is only allowed for indefinite integrals.", call. = FALSE)
  if (isTRUE(.integral_is_zero(I))) return(S_Zero)
  function_ <- I$args[[1L]]
  if (function_$type == "Sum") {
    sl <- function_$args[-1L]
    if (any(vapply(lims, function(l) any(vapply(sl[[1L]]$args, function(z) .mi_eq(z, l$args[[1L]]), NA)), NA)))
      stop("Limit of the sum cannot be an integration variable.", call. = FALSE)
    if (any(vapply(sl[[1L]]$args[-1L], function(z) isTRUE(ask_fact(z, "infinite")), NA))) return(I)
    inner <- .integral_doit(new_sym("Integral", c(list(function_$args[[1L]]), lims)), hints)
    return(doit(sym_rebuild(function_, c(list(inner), sl))))
  }
  function_ <- .map_bottom_up(function_, function(z) {
    if (z$type == "Function" && z$name == "Heaviside" && length(z$args) == 2L && !is_one_num(.mul(list(z$args[[2L]], S_Two))))
      return(fcall("Heaviside", list(z$args[[1L]])))
    z
  })
  if (deep) function_ <- doit(function_)
  if (isTRUE(ask_fact(function_, "zero"))) return(S_Zero)
  if (inherits(function_, "symmatrix")) return(.mat_map(function_, function(f) .integral_doit(new_sym("Integral", c(list(f), lims)), hints)))
  # real, positive or negative integration variables for real limits
  reps <- list()
  for (l in lims) {
    if (length(l$args) != 3L) next
    x <- l$args[[1L]]; ab <- l$args[2:3]
    d <- NULL
    if (all(vapply(ab, function(i) isTRUE(ask_fact(i, "nonnegative")), NA)) && !isTRUE(ask_fact(x, "nonnegative"))) d <- Dummy(positive = TRUE)
    else if (all(vapply(ab, function(i) isTRUE(ask_fact(i, "nonpositive")), NA)) && !isTRUE(ask_fact(x, "nonpositive"))) d <- Dummy(negative = TRUE)
    else if (all(vapply(ab, function(i) isTRUE(ask_fact(i, "real")), NA)) && !isTRUE(ask_fact(x, "real"))) d <- Dummy(real = TRUE)
    if (!is.null(d)) reps[[length(reps) + 1L]] <- list(x, d)
  }
  if (length(reps)) {
    undo <- lapply(reps, function(p) list(p[[2L]], p[[1L]]))
    did <- .integral_doit(xreplace(new_sym("Integral", c(list(function_), lims)), reps), hints, deep = deep, final = final)
    if (is_sym(did)) return(xreplace(did, undo))
    return(lapply(did, function(z) xreplace(z, undo)))
  }
  undone <- list(); ulj <- character()
  for (xab in lims) {
    la <- length(xab$args)
    x <- xab$args[[1L]]
    uli <- if (la == 1L) x$key else if (la == 2L) names(.free_sym_list(xab$args[[2L]])) else union(names(.free_sym_list(xab$args[[2L]])), names(.free_sym_list(xab$args[[3L]])))
    if (x$key %in% ulj || any(vapply(undone, function(v) v$args[[1L]]$key %in% uli, NA))) {
      undone[[length(undone) + 1L]] <- xab
      ulj <- union(ulj, uli)
      function_ <- new_sym("Integral", list(function_, xab))
      ff <- .integral_factor(function_)
      if (ff$type != "Integral") function_ <- ff
      next
    }
    real_ends <- if (la < 3L) all(vapply(xab$args, function(z) isTRUE(ask_fact(z, "extended_real")), NA))
                 else all(vapply(xab$args[2:3], function(z) isTRUE(ask_fact(z, "extended_real")) && !isTRUE(ask_fact(z, "infinite")), NA))
    if (.has_function(function_, c("Abs", "sign")) && real_ends) {
      xr <- Dummy("xr", real = TRUE)
      function_ <- xreplace(rewrite(xreplace(function_, list(list(x, xr))), "Piecewise"), list(list(xr, x)))
    } else if (.has_function(function_, c("Min", "Max"))) function_ <- rewrite(function_, "Piecewise")
    if (.has_type(function_, "Piecewise") && function_$type != "Piecewise") function_ <- piecewise_fold(function_)
    antideriv <- NULL
    if (function_$type == "Piecewise") {
      antideriv <- if (la == 1L) .pw_eval_integral(function_, x, kw) else do.call(.sympy_eval_integral, c(list(function_, x), kw))
    } else {
      try_meijerg <- function(f, xab) {
        if (length(xab$args) != 3L || identical(meijerg, FALSE)) return(NULL)
        a <- xab$args[[2L]]; b <- xab$args[[3L]]
        res <- tryCatch(meijerint_definite(f, x, a, b), NotImplementedError = function(e) NULL)
        if (is.null(res)) return(NULL)
        if (conds == "piecewise") return(.piecewise_noeval(list(list(res[[1L]], res[[2L]]), list(new_sym("Integral", list(f, xab)), S_true))))
        if (conds == "separate") {
          if (length(lims) != 1L) stop("conds=separate not supported in multiple integrals", call. = FALSE)
          return(res)
        }
        res[[1L]]
      }
      meijerg1 <- meijerg
      infl <- function(z) .contains_key(z, S_Infinity$key) || .contains_key(z, S_NegativeInfinity$key)
      if (!identical(meijerg, FALSE) && la == 3L && isTRUE(ask_fact(xab$args[[2L]], "extended_real")) && isTRUE(ask_fact(xab$args[[3L]], "extended_real")) &&
          (infl(xab$args[[2L]]) || infl(xab$args[[3L]]))) {
        ret <- try_meijerg(function_, xab)
        if (!is.null(ret)) {
          if (!is_sym(ret)) return(ret)
          function_ <- ret
          next
        }
        meijerg1 <- FALSE
      }
      if (identical(meijerg1, FALSE) && isTRUE(meijerg)) antideriv <- NULL
      else {
        antideriv <- .antideriv_eval(function_, x, kw)
        if (is.null(antideriv) && isTRUE(meijerg)) {
          ret <- try_meijerg(function_, xab)
          if (!is.null(ret)) {
            if (!is_sym(ret)) return(ret)
            function_ <- ret
            next
          }
        }
      }
    }
    if (la == 3L && final && !is.null(antideriv) && antideriv$type != "Integral") antideriv <- .add_atan_floor_terms(antideriv, x)
    if (is.null(antideriv)) {
      undone[[length(undone) + 1L]] <- xab
      function_ <- .integral_factor(new_sym("Integral", list(function_, xab)))
      ff <- .integral_factor(function_)
      if (ff$type != "Integral") function_ <- ff
      next
    }
    if (la == 1L) { function_ <- antideriv; next }
    if (la == 3L) { a <- xab$args[[2L]]; b <- xab$args[[3L]] } else { a <- NULL; b <- xab$args[[2L]] }
    if (deep) { if (!is.null(a)) a <- doit(a); b <- doit(b) }
    is_indef_int <- function(g) g$type == "Integral" && any(vapply(g$args[-1L], function(l) length(l$args) == 1L && .mi_eq(l$args[[1L]], x), NA))
    eval_factored <- function(f) .mul(lapply(.margs(f), function(g) if (is_indef_int(g)) .expr_eval_interval(g, x, a, b) else g))
    integrals <- list(); others <- list(); pws <- list()
    for (f in .mi_make_args_add(antideriv)) {
      ms <- .margs(f)
      if (any(vapply(ms, is_indef_int, NA))) integrals[[length(integrals) + 1L]] <- f
      else if (any(vapply(ms, function(g) g$type == "Piecewise", NA))) pws[[length(pws) + 1L]] <- piecewise_fold(f)
      else others[[length(others) + 1L]] <- f
    }
    uneval <- .add(lapply(integrals, eval_factored))
    r <- tryCatch({
      evalued <- .expr_eval_interval(.add(others), x, a, b)
      evalued_pw <- .eval_interval(piecewise_fold(.add(pws)), x, a, b)
      .add(list(uneval, evalued, evalued_pw))
    }, NotImplementedError = function(e) NULL)
    if (is.null(r)) {
      undone[[length(undone) + 1L]] <- xab
      function_ <- new_sym("Integral", list(function_, xab))
      ff <- .integral_factor(function_)
      if (ff$type != "Integral") function_ <- ff
    } else function_ <- r
  }
  function_
}

# Integral._eval_integral with a Piecewise integrand: piecewise_integrate
.antideriv_eval <- function(f, x, kw) {
  if (f$type == "Piecewise")
    return(do.call(Piecewise, lapply(f$args, function(p) list(do.call(integrate, c(list(p$args[[1L]], x), kw)), p$args[[2L]]))))
  .antideriv_sympy(f, x, kw)
}
