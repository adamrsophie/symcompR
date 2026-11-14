# Adaptive numerical evaluation after sympy/core/evalf.py. evalf(x, prec,
# options) returns list(re, im, re_acc, im_acc): re and im are binary
# floats list(m, e) (evalf-mp.R) or NULL for an exact zero, and the
# accuracies are log2 of the relative error of each part. Complex infinity
# is the string "zoo". The options are an environment, since SymPy mutates
# options["maxprec"] and the symbol cache during the evaluation.

.ev_MINUS_INF <- -Inf
.ev_DEFAULT_MAXPREC <- 333
.ev_LG10 <- log2(10)

.ev_fone <- list(m = gmp::as.bigz(1), e = 0)
.ev_fhalf <- list(m = gmp::as.bigz(1), e = -1)
.ev_fnone <- list(m = gmp::as.bigz(-1), e = 0)
# mpmath's special values: zero mantissa with a marker exponent and bit count
.ev_finf <- list(m = gmp::as.bigz(0), e = -456, bc = -2, sp = "inf")
.ev_fninf <- list(m = gmp::as.bigz(0), e = -789, bc = -3, sp = "-inf")
.ev_fnan <- list(m = gmp::as.bigz(0), e = -123, bc = -1, sp = "nan")

.ev_res <- function(re = NULL, im = NULL, re_acc = NULL, im_acc = NULL)
  list(re = re, im = im, re_acc = re_acc, im_acc = im_acc)
.ev_none <- function() .ev_res()
.ev_is_zoo <- function(r) is.character(r)

.ev_condition <- function(cls, msg) structure(class = c(cls, "sympy_error", "error", "condition"),
                                               list(message = msg, call = NULL))
.ev_nie <- function(msg = "") stop(.ev_condition("NotImplementedError", msg))
.ev_valueerror <- function(msg) stop(.ev_condition("ValueError", msg))

.ev_sp <- function(x) !is.null(x) && !is.null(x$sp)
.ev_is_fzero <- function(x) !is.null(x) && is.null(x$sp) && !isTRUE(x$sz) && x$m == 0
.ev_bc <- function(x) if (!is.null(x$bc)) x$bc else .bitlen(abs(x$m))
.ev_sign <- function(x) if (isTRUE(x$sz)) x$s else if (x$m < 0) 1L else 0L

# fastlog: exponent plus mantissa width
.ev_fastlog <- function(x) {
  if (is.null(x) || .ev_is_fzero(x)) return(.ev_MINUS_INF)
  x$e + .ev_bc(x)
}

# iszero: an exact zero, or with scaled = TRUE a scaled zero
.ev_iszero <- function(x, scaled = FALSE) {
  if (!scaled) return(is.null(x) || .ev_is_fzero(x))
  !is.null(x) && isTRUE(x$sz)
}

# scaled_zero(mag, sign): the power of two 2^mag marked as a zero of that
# magnitude with precision -1; scaled_zero(z) unmarks it
.ev_scaled_zero <- function(mag, sign = 1) {
  if (is.list(mag)) {
    if (!isTRUE(mag$sz)) .ev_valueerror("scaled zero expects int or scaled_zero tuple.")
    return(list(m = gmp::as.bigz(if (mag$s == 1L) -1 else 1), e = mag$e))
  }
  if (!(sign %in% c(-1, 1))) .ev_valueerror("sign must be +/-1")
  list(list(m = gmp::as.bigz(1), e = mag, sz = TRUE, s = if (sign == 1) 0L else 1L), -1)
}

# rounding to p bits toward zero (mpmath's round_down, the default
# round_fast of libmpf) or to nearest
.ev_round <- function(m, e, p, rnd = "n") {
  if (rnd == "n") return(.mp_round(m, e, p))
  m <- gmp::as.bigz(m)
  if (m == 0) return(.mp_zero)
  b <- .bitlen(abs(m))
  if (b > p) {
    s <- b - p
    q <- gmp::divq.bigz(abs(m), .pow2(s))
    m <- if (m < 0) -q else q
    e <- e + s
  }
  .mp_round(m, e, .bitlen(abs(m)))
}
.ev_norm <- function(a, p, rnd = "n") if (.ev_sp(a)) a else .ev_round(a$m, a$e, p, rnd)
# a value computed with 20 guard bits, then rounded down
.ev_down <- function(a, p) .ev_round(a$m, a$e, p, "d")
.ev_mul <- function(a, b, p, rnd = "d") .ev_round(a$m * b$m, a$e + b$e, p, rnd)
.ev_add <- function(a, b, p, rnd = "d") {
  if (a$m == 0) return(.ev_round(b$m, b$e, p, rnd))
  if (b$m == 0) return(.ev_round(a$m, a$e, p, rnd))
  if (a$e < b$e) { t <- a; a <- b; b <- t }
  if (.mp_mag(a) - .mp_mag(b) > p + 4) {
    k <- p + 4 + max(0, p - .bitlen(abs(a$m)))
    return(.ev_round(a$m * .pow2(k) + sign(as.numeric(b$m)), a$e - k, p, rnd))
  }
  .ev_round(a$m * .pow2(a$e - b$e) + b$m, b$e, p, rnd)
}
.ev_neg <- function(a) {
  if (.ev_sp(a)) return(switch(a$sp, inf = .ev_fninf, `-inf` = .ev_finf, a))
  list(m = -a$m, e = a$e)
}
.ev_abs <- function(a) {
  if (.ev_sp(a)) return(if (a$sp == "-inf") .ev_finf else a)
  list(m = abs(a$m), e = a$e)
}
.ev_cmp <- function(a, b) {
  if (.ev_sp(a) || .ev_sp(b)) {
    v <- function(x) if (.ev_sp(x)) switch(x$sp, inf = Inf, `-inf` = -Inf, NaN) else .mp_to_dbl(x)
    va <- v(a); vb <- v(b)
    return(if (is.nan(va) || is.nan(vb)) 0 else sign(va - vb))
  }
  .mp_cmp(a, b)
}
.ev_from_int <- function(n, p) .ev_round(gmp::as.bigz(n), 0, p, "d")
.ev_from_rational <- function(q, p) {
  q <- gmp::as.bigq(q)
  if (gmp::denominator(q) == 1) return(.ev_from_int(gmp::numerator(q), p))
  .mp_from_q_down(q, p)
}
# to_int(x, round_nearest)
.ev_to_int <- function(x, rnd = "n") {
  if (.ev_sp(x)) .ev_valueerror("cannot convert inf or nan to int")
  if (is.null(x) || x$m == 0) return(gmp::as.bigz(0))
  if (x$e >= 0) return(x$m * .pow2(x$e))
  if (rnd == "d") { q <- gmp::divq.bigz(abs(x$m), .pow2(-x$e)); return(if (x$m < 0) -q else q) }
  .shr_round(x$m, -x$e)
}
.ev_pi <- function(p) .ev_down(.mp_pi(p + 20), p)
.ev_e <- function(p) .ev_down(.mp_exp(.mpf(1, 0), p + 20), p)
.ev_sqrt <- function(a, p) .ev_down(.mp_sqrt(a, p + 20), p)
.ev_exp <- function(a, p) .ev_down(.mp_exp(a, p + 20), p)
.ev_ln <- function(a, p) {
  if (a$m == 0) return(.ev_fninf)
  .mp_log(a, p)
}
.ev_pow_int <- function(a, n, p) {
  r <- .mp_pow_int(a, n, p + 20)
  .ev_down(r, p)
}
.ev_cos <- function(a, p) .mp_cos_sin(a, p)$cos
.ev_sin <- function(a, p) .mp_cos_sin(a, p)$sin
.ev_tan <- function(a, p) {
  cs <- .mp_cos_sin(a, p + 20)
  .mp_div(cs$sin, cs$cos, p)
}
.ev_atan2 <- function(y, x, p) .ev_down(.mp_atan2(y, x, p + 20), p)
.ev_mpc_down <- function(z, p) list(re = .ev_down(z$re, p), im = .ev_down(z$im, p))

# The sym for a special value or a Float of the given precision
.ev_float <- function(x, p) {
  if (.ev_sp(x)) return(switch(x$sp, inf = S_Infinity, `-inf` = S_NegativeInfinity, S_NaN))
  if (isTRUE(x$sz)) x <- .ev_scaled_zero(x)
  if (x$m == 0) return(S_Zero)
  if (abs(x$e) > 2^52) stop("exponent beyond the exact range of multiprecision numbers", call. = FALSE)
  .mp_float(.mp_round(x$m, x$e, p), p)
}

# pure_complex(v): list(h, c) for v = h + c*I with Numbers h and c
.ev_is_Number <- function(e) e$type %in% .NUMBER_TYPES ||
  (e$type == "Constant" && e$name %in% c("oo", "-oo", "nan"))
.ev_as_coeff_Mul <- function(t) {
  if (t$type == "Mul") {
    a <- .sympy_args(t)
    if (.ev_is_Number(a[[1L]])) {
      rest <- a[-1L]
      return(list(a[[1L]], if (length(rest) == 1L) rest[[1L]] else new_sym("Mul", rest)))
    }
  }
  if (.ev_is_Number(t)) return(list(t, S_One))
  list(S_One, t)
}
.ev_pure_complex <- function(v, or_real = FALSE) {
  if (v$type == "Add") {
    a <- .sympy_args(v)
    if (.ev_is_Number(a[[1L]])) {
      h <- a[[1L]]; rest <- a[-1L]
      t <- if (length(rest) == 1L) rest[[1L]] else new_sym("Add", rest)
    } else { h <- S_Zero; t <- v }
  } else if (.ev_is_Number(v)) {
    h <- v; t <- S_Zero
  } else { h <- S_Zero; t <- v }
  if (!(t$type == "Integer" && t$value == 0)) {
    ci <- .ev_as_coeff_Mul(t)
    if (ci[[2L]]$type == "Constant" && ci[[2L]]$name == "I") return(list(h, ci[[1L]]))
  } else if (or_real) return(list(h, S_Zero))
  NULL
}

# ---------------------------------------------------------------------------
# Accuracy bookkeeping

.ev_complex_accuracy <- function(result) {
  if (.ev_is_zoo(result)) return(Inf)
  re <- result$re; im <- result$im
  if (is.null(im)) {
    if (is.null(re)) return(Inf)
    return(result$re_acc)
  }
  if (is.null(re)) return(result$im_acc)
  re_size <- .ev_fastlog(re); im_size <- .ev_fastlog(im)
  absolute_error <- max(re_size - result$re_acc, im_size - result$im_acc)
  relative_error <- absolute_error - max(re_size, im_size)
  -relative_error
}

.ev_get_abs <- function(expr, prec, options) {
  result <- .ev_evalf(expr, prec + 2, options)
  if (.ev_is_zoo(result)) return(.ev_res(.ev_finf, NULL, prec, NULL))
  re <- result$re; im <- result$im; re_acc <- result$re_acc; im_acc <- result$im_acc
  if (is.null(re)) { t <- re; re <- im; im <- t; t <- re_acc; re_acc <- im_acc; im_acc <- t }
  if (!is.null(im)) {
    if (.ev_is_number(expr)) {
      a <- .mpc_abs(list(re = .ev_abs(re), im = .ev_abs(im)), prec)
      acc <- min(prec, .ev_complex_accuracy(result))
      return(.ev_res(a, NULL, acc, NULL))
    }
    if (!is.null(options$subs)) return(.ev_res(.ev_down(.mpc_abs(list(re = re, im = im), prec + 20), prec), NULL, re_acc, NULL))
    .ev_nie()
  }
  if (!is.null(re)) return(.ev_res(.ev_abs(re), NULL, re_acc, NULL))
  .ev_none()
}

.ev_get_complex_part <- function(expr, no, prec, options) {
  workprec <- prec
  i <- 0
  repeat {
    res <- .ev_evalf(expr, workprec, options)
    if (.ev_is_zoo(res)) return(.ev_res(.ev_fnan, NULL, prec, NULL))
    value <- if (no == 0) res$re else res$im
    accuracy <- if (no == 0) res$re_acc else res$im_acc
    if (is.null(value) || accuracy >= prec || -value$e > prec)
      return(.ev_res(value, NULL, accuracy, NULL))
    workprec <- workprec + max(30, 2^i)
    i <- i + 1
  }
}

.ev_finalize_complex <- function(re, im, prec) {
  if (.ev_is_fzero(re) && .ev_is_fzero(im)) .ev_valueerror("got complex zero with unknown accuracy")
  if (.ev_is_fzero(re)) return(.ev_res(NULL, im, NULL, prec))
  if (.ev_is_fzero(im)) return(.ev_res(re, NULL, prec, NULL))
  size_re <- .ev_fastlog(re); size_im <- .ev_fastlog(im)
  if (size_re > size_im) {
    re_acc <- prec; im_acc <- prec + min(-(size_re - size_im), 0)
  } else {
    im_acc <- prec; re_acc <- prec + min(-(size_im - size_re), 0)
  }
  .ev_res(re, im, re_acc, im_acc)
}

.ev_chop_parts <- function(value, prec) {
  if (.ev_is_zoo(value)) return(value)
  re <- value$re; im <- value$im; re_acc <- value$re_acc; im_acc <- value$im_acc
  if (!is.null(re) && !.ev_sp(re) && .ev_fastlog(re) < -prec + 4) { re <- NULL; re_acc <- NULL }
  if (!is.null(im) && !.ev_sp(im) && .ev_fastlog(im) < -prec + 4) { im <- NULL; im_acc <- NULL }
  if (!is.null(re) && !is.null(im)) {
    delta <- .ev_fastlog(re) - .ev_fastlog(im)
    if (re_acc < 2 && (delta - re_acc <= -prec + 4)) { re <- NULL; re_acc <- NULL }
    if (im_acc < 2 && (delta - im_acc >= prec - 4)) { im <- NULL; im_acc <- NULL }
  }
  .ev_res(re, im, re_acc, im_acc)
}

.ev_check_target <- function(expr, result, prec) {
  a <- .ev_complex_accuracy(result)
  if (a < prec)
    stop(.ev_condition("PrecisionExhausted", paste0(
      "Failed to distinguish the expression: \n\n", sstr(expr),
      "\n\nfrom zero. Try simplifying the input, using chop=True, or providing a higher maxn for evalf")))
}

# ---------------------------------------------------------------------------
# Integer parts

.ev_get_integer_part <- function(expr, no, options, return_ints = FALSE) {
  assumed_size <- 30
  result <- .ev_evalf(expr, assumed_size, options)
  if (.ev_is_zoo(result)) .ev_valueerror("Cannot get integer part of Complex Infinity")
  ire <- result$re; iim <- result$im; ire_acc <- result$re_acc; iim_acc <- result$im_acc
  if (!is.null(ire) && !is.null(iim)) {
    gap <- max(.ev_fastlog(ire) - ire_acc, .ev_fastlog(iim) - iim_acc)
  } else if (!is.null(ire)) {
    gap <- .ev_fastlog(ire) - ire_acc
  } else if (!is.null(iim)) {
    gap <- .ev_fastlog(iim) - iim_acc
  } else {
    if (return_ints) return(list(gmp::as.bigz(0), gmp::as.bigz(0)))
    return(.ev_none())
  }
  margin <- 10
  if (gap >= -margin) {
    prec <- margin + assumed_size + gap
    r <- .ev_evalf(expr, prec, options)
    ire <- r$re; iim <- r$im; ire_acc <- r$re_acc; iim_acc <- r$im_acc
  } else {
    prec <- assumed_size
  }
  calc_part <- function(re_im, nexpr) {
    is_int <- nexpr$e == 0
    nint <- .ev_to_int(nexpr)
    if (is_int) {
      r <- .ev_evalf(.add(list(re_im, make_rational(gmp::as.bigq(-nint)))), 10, options)
      ire <- r$re
      size <- -.ev_fastlog(ire) + 2
      if (size > prec && is.finite(size)) {
        r <- .ev_evalf(re_im, size, options)
        ire <- r$re
        nexpr <- ire
      }
      nint <- .ev_to_int(nexpr)
      is_int <- !is.null(ire) && ire$e == 0
    }
    if (!is_int) {
      s <- options$subs
      if (length(s)) {
        is_int_reim <- function(v) {
          v <- as_sym(v)
          ok <- function(z) z$type == "Integer" || (z$type == "Float" && .mp_is_int(.float_mp(z)))
          if (ok(v)) return(TRUE)
          ri <- tryCatch(.as_real_imag(v), error = function(err) NULL)
          !is.null(ri) && ok(ri[[1L]]) && ok(ri[[2L]])
        }
        if (all(vapply(s, is_int_reim, NA))) re_im <- .ev_subs_all(re_im, s)
      }
      re_im <- new_sym("Add", list(re_im, make_rational(gmp::as.bigq(-nint))))
      r <- .ev_evalf(re_im, 10, options)
      x <- r$re
      chk <- tryCatch({ .ev_check_target(re_im, .ev_res(x, NULL, r$re_acc, NULL), 3); TRUE },
                      PrecisionExhausted = function(err) FALSE)
      if (!chk) {
        # re_im.equals(0); SymPy's equals simplifies re() and im() of
        # a real expression, so the evaluated part is tried as well
        z <- function(v) isTRUE(tryCatch(.expr_equals(v, S_Zero), error = function(err) NA))
        if (!z(re_im) && !z(tryCatch(.add(lapply(re_im$args, function(a)
              if (a$type == "Function" && a$name %in% c("re", "im")) fcall(a$name, a$args) else a)),
              error = function(err) re_im)))
          stop(.ev_condition("PrecisionExhausted", ""))
        x <- .mp_zero
      }
      xx <- if (is.null(x)) .mp_zero else x
      nint <- nint + as.integer(no * (.ev_cmp(xx, .mp_zero) == no))
    }
    list(.ev_from_int(nint, Inf), Inf)
  }
  re_ <- im_ <- re_acc <- im_acc <- NULL
  if (!is.null(ire) && !.ev_is_fzero(ire)) {
    cp <- calc_part(fnode("re", list(expr)), ire); re_ <- cp[[1L]]; re_acc <- cp[[2L]]
  }
  if (!is.null(iim) && !.ev_is_fzero(iim)) {
    cp <- calc_part(fnode("im", list(expr)), iim); im_ <- cp[[1L]]; im_acc <- cp[[2L]]
  }
  if (return_ints) return(list(.ev_to_int(re_), .ev_to_int(im_)))
  .ev_res(re_, im_, re_acc, im_acc)
}
.ev_from_int_exact <- function(n) .mp_round(gmp::as.bigz(n), 0, Inf)

# ---------------------------------------------------------------------------
# Arithmetic

.ev_add_terms <- function(terms, prec, target_prec) {
  terms <- Filter(function(t) !.ev_iszero(t[[1L]]), terms)
  if (!length(terms)) return(list(NULL, NULL))
  if (length(terms) == 1L) return(terms[[1L]])
  special <- Filter(.ev_sp, lapply(terms, `[[`, 1L))
  if (length(special)) {
    # evalf(Add(*special)): oo + -oo or a nan is nan
    sps <- vapply(special, `[[`, "", "sp")
    v <- if ("nan" %in% sps || ("inf" %in% sps && "-inf" %in% sps)) .ev_fnan else if ("inf" %in% sps) .ev_finf else .ev_fninf
    return(list(v, prec + 4))
  }
  working_prec <- 2 * prec
  sum_man <- gmp::as.bigz(0); sum_exp <- 0
  absolute_err <- numeric(0)
  for (t in terms) {
    x <- t[[1L]]; accuracy <- t[[2L]]
    if (isTRUE(x$sz)) {
      # a scaled zero's sign is a list, which is always true in Python
      man <- gmp::as.bigz(-1); ex <- x$e; bc <- 1
    } else {
      man <- x$m; ex <- x$e; bc <- .bitlen(abs(man))
    }
    absolute_err <- c(absolute_err, bc + ex - accuracy)
    delta <- ex - sum_exp
    if (ex >= sum_exp) {
      if ((delta > working_prec) && (sum_man == 0 || delta - .bitlen(abs(sum_man)) > working_prec)) {
        sum_man <- man; sum_exp <- ex
      } else {
        sum_man <- sum_man + man * .pow2(delta)
      }
    } else {
      delta <- -delta
      if (delta - bc > working_prec) {
        if (sum_man == 0) { sum_man <- man; sum_exp <- ex }
      } else {
        sum_man <- sum_man * .pow2(delta) + man
        sum_exp <- ex
      }
    }
  }
  absolute_error <- max(absolute_err)
  if (sum_man == 0) return(.ev_scaled_zero(absolute_error))
  sum_bc <- .bitlen(abs(sum_man))
  sum_accuracy <- sum_exp + sum_bc - absolute_error
  list(.mp_round(sum_man, sum_exp, target_prec), sum_accuracy)
}

.ev_evalf_add <- function(v, prec, options) {
  res <- .ev_pure_complex(v)
  if (!is.null(res)) {
    a <- .ev_evalf(res[[1L]], prec, options)
    b <- .ev_evalf(res[[2L]], prec, options)
    return(.ev_res(a$re, b$re, a$re_acc, b$re_acc))
  }
  oldmaxprec <- if (is.null(options$maxprec)) .ev_DEFAULT_MAXPREC else options$maxprec
  i <- 0
  target_prec <- prec
  args <- .sympy_args(v)
  repeat {
    options$maxprec <- min(oldmaxprec, 2 * prec)
    terms <- lapply(args, .ev_evalf, prec + 10, options)
    n <- sum(vapply(terms, .ev_is_zoo, NA))
    if (n >= 2) return(.ev_res(.ev_fnan, NULL, prec, NULL))
    tl <- Filter(function(a) !.ev_is_zoo(a), terms)
    ra <- .ev_add_terms(lapply(Filter(function(a) !is.null(a$re), tl), function(a) list(a$re, a$re_acc)), prec, target_prec)
    ia <- .ev_add_terms(lapply(Filter(function(a) !is.null(a$im), tl), function(a) list(a$im, a$im_acc)), prec, target_prec)
    re <- ra[[1L]]; re_acc <- ra[[2L]]; im <- ia[[1L]]; im_acc <- ia[[2L]]
    if (n == 1) {
      if (.ev_sp(re) || .ev_sp(im)) return(.ev_res(.ev_fnan, NULL, prec, NULL))
      return("zoo")
    }
    acc <- .ev_complex_accuracy(.ev_res(re, im, re_acc, im_acc))
    if (acc >= target_prec) {
      if (isTRUE(options$verbose)) cat("ADD: wanted", target_prec, "accurate bits, got", re_acc, im_acc, "\n")
      break
    } else {
      if ((prec - target_prec) > options$maxprec) break
      prec <- prec + max(10 + 2^i, target_prec - acc)
      i <- i + 1
      if (isTRUE(options$verbose)) cat("ADD: restarting with prec", prec, "\n")
    }
  }
  options$maxprec <- oldmaxprec
  if (.ev_iszero(re, scaled = TRUE)) re <- .ev_scaled_zero(re)
  if (.ev_iszero(im, scaled = TRUE)) im <- .ev_scaled_zero(im)
  .ev_res(re, im, re_acc, im_acc)
}

.ev_evalf_mul <- function(v, prec, options) {
  res <- .ev_pure_complex(v)
  if (!is.null(res)) {
    b <- .ev_evalf(res[[2L]], prec, options)
    return(.ev_res(NULL, b$re, NULL, b$re_acc))
  }
  args <- .sympy_args(v)
  flat_factors <- list()
  # the exact coefficient (a + b*I)/q
  ca <- gmp::as.bigz(1); cb <- gmp::as.bigz(0); cq <- gmp::as.bigz(1)
  have_exact <- FALSE
  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]; i <- i + 1L
    if (arg$type == "Mul") { args <- c(args, .sympy_args(arg)); next }
    pair <- .ev_pure_complex(arg)
    if (!is.null(pair) && !(pair[[1L]]$type == "Integer" && pair[[1L]]$value == 0) &&
        all(vapply(pair, is_exact_number, NA))) {
      r <- gmp::as.bigq(pair[[1L]]$value); im <- gmp::as.bigq(pair[[2L]]$value)
      c <- gmp::numerator(r) * gmp::denominator(im); d <- gmp::numerator(im) * gmp::denominator(r)
      s <- gmp::denominator(r) * gmp::denominator(im)
      real <- ca * c - cb * d; imag <- ca * d + cb * c; denom <- cq * s
      g <- gmp::gcd(gmp::gcd(real, imag), denom)
      ca <- gmp::divq.bigz(real, g); cb <- gmp::divq.bigz(imag, g); cq <- gmp::divq.bigz(denom, g)
      have_exact <- TRUE
      next
    }
    flat_factors[[length(flat_factors) + 1L]] <- arg
  }
  has_zero <- FALSE
  special <- list()
  for (arg in flat_factors) {
    result <- .ev_evalf(arg, prec, options)
    if (.ev_is_zoo(result)) { special[[length(special) + 1L]] <- S_ComplexInfinity; next }
    if (is.null(result$re)) {
      if (is.null(result$im)) has_zero <- TRUE
      next
    }
    num <- .ev_float(result$re, 1)
    if (is_nan(num)) return(.ev_res(.ev_fnan, NULL, prec, NULL))
    if (is_inf_const(num)) special[[length(special) + 1L]] <- num
  }
  if (length(special)) {
    if (has_zero) return(.ev_res(.ev_fnan, NULL, prec, NULL))
    return(.ev_evalf(.mul(special), prec + 4, .ev_options()))
  }
  if (has_zero) return(.ev_none())
  acc <- prec
  working_prec <- prec + length(flat_factors) + as.integer(have_exact) + 5
  man <- gmp::as.bigz(1); ex <- 0; bc <- 1
  direction <- 0
  complex_factors <- list()
  LAST <- length(flat_factors)
  for (k in seq_len(LAST + 1L)) {
    if (k <= LAST) {
      r <- .ev_evalf(flat_factors[[k]], working_prec, options)
      re <- r$re; im <- r$im; re_acc <- r$re_acc; im_acc <- r$im_acc
    } else if (have_exact) {
      a <- .ev_evalf(make_rational(gmp::as.bigq(ca, cq)), working_prec, options)
      b <- .ev_evalf(make_rational(gmp::as.bigq(cb, cq)), working_prec, options)
      re <- a$re; re_acc <- a$re_acc; im <- b$re; im_acc <- b$re_acc
    } else break
    if (!is.null(re) && !is.null(im)) {
      complex_factors[[length(complex_factors) + 1L]] <- list(re, im, re_acc, im_acc)
      next
    } else if (!is.null(re)) {
      w <- re; w_acc <- re_acc
    } else if (!is.null(im)) {
      w <- im; w_acc <- im_acc
      direction <- direction + 1
    } else return(.ev_none())
    direction <- direction + 2 * .ev_sign(w)
    man <- man * abs(w$m); ex <- ex + w$e; bc <- bc + .ev_bc(w)
    while (bc > 3 * working_prec) {
      man <- gmp::divq.bigz(man, .pow2(working_prec))
      ex <- ex + working_prec; bc <- bc - working_prec
    }
    acc <- min(acc, w_acc)
  }
  sgn <- (direction %% 4) %/% 2
  if (!length(complex_factors)) {
    v <- .mp_round(if (sgn == 1) -man else man, ex, prec)
    if (direction %% 2 == 1) return(.ev_res(NULL, v, NULL, acc))
    return(.ev_res(v, NULL, acc, NULL))
  }
  if (!(man == 1 && ex == 0 && bc == 1)) {
    re <- list(m = if (sgn == 1) -man else man, e = ex); im <- .mp_zero
    i0 <- 1L
  } else {
    w <- complex_factors[[1L]]
    acc <- min(acc, .ev_complex_accuracy(.ev_res(w[[1L]], w[[2L]], w[[3L]], w[[4L]])))
    re <- w[[1L]]; im <- w[[2L]]
    i0 <- 2L
  }
  for (w in complex_factors[seq_len(length(complex_factors) - i0 + 1L) + i0 - 1L]) {
    acc <- min(acc, .ev_complex_accuracy(.ev_res(w[[1L]], w[[2L]], w[[3L]], w[[4L]])))
    use_prec <- working_prec
    A <- .ev_mul(re, w[[1L]], use_prec)
    B <- .ev_mul(.ev_neg(im), w[[2L]], use_prec)
    C <- .ev_mul(re, w[[2L]], use_prec)
    D <- .ev_mul(im, w[[1L]], use_prec)
    re <- .ev_add(A, B, use_prec)
    im <- .ev_add(C, D, use_prec)
  }
  if (isTRUE(options$verbose)) cat("MUL: wanted", prec, "accurate bits, got", acc, "\n")
  if (direction %% 2 == 1) { t <- .ev_neg(im); im <- re; re <- t }
  .ev_res(re, im, acc, acc)
}

.ev_evalf_pow <- function(v, prec, options) {
  target_prec <- prec
  base <- v$args[[1L]]; ex <- v$args[[2L]]
  if (ex$type == "Integer") {
    p <- gmp::numerator(gmp::as.bigq(ex$value))
    if (p == 0) return(.ev_res(.ev_fone, NULL, prec, NULL))
    prec <- prec + floor(log2(abs(as.numeric(p))))
    result <- .ev_evalf(base, prec + 5, options)
    if (.ev_is_zoo(result)) {
      if (p < 0) return(.ev_none())
      return(result)
    }
    re <- result$re; im <- result$im
    if (!is.null(re) && is.null(im)) return(.ev_res(.ev_pow_int_sp(re, p, target_prec), NULL, target_prec, NULL))
    if (!is.null(im) && is.null(re)) {
      z <- .ev_pow_int_sp(im, p, target_prec)
      case <- as.integer(gmp::mod.bigz(p, 4))
      return(switch(case + 1L, .ev_res(z, NULL, target_prec, NULL), .ev_res(NULL, z, NULL, target_prec),
                    .ev_res(.ev_neg(z), NULL, target_prec, NULL), .ev_res(NULL, .ev_neg(z), NULL, target_prec)))
    }
    if (is.null(re)) {
      if (p < 0) return("zoo")
      return(.ev_none())
    }
    z <- .mpc_pow_int(list(re = re, im = im), p, prec)
    return(.ev_finalize_complex(z$re, z$im, target_prec))
  }
  result <- .ev_evalf(base, prec + 5, options)
  if (.ev_is_zoo(result)) {
    if (is_exact_number(ex)) {
      if (gmp::as.bigq(ex$value) < 0) return(.ev_none())
      return(result)
    }
    .ev_nie()
  }
  if (ex$type == "Rational" && gmp::as.bigq(ex$value) == gmp::as.bigq(1, 2)) {
    xre <- result$re; xim <- result$im
    if (!is.null(xim)) {
      z <- .mpc_sqrt(list(re = if (is.null(xre)) .mp_zero else xre, im = xim), prec)
      return(.ev_finalize_complex(z$re, z$im, prec))
    }
    if (is.null(xre)) return(.ev_none())
    if (.ev_cmp(xre, .mp_zero) < 0) return(.ev_res(NULL, .ev_sqrt(.ev_neg(xre), prec), NULL, prec))
    return(.ev_res(.ev_sqrt_sp(xre, prec), NULL, prec, NULL))
  }
  prec <- prec + 10
  result <- .ev_evalf(ex, prec, options)
  if (.ev_is_zoo(result)) return(.ev_res(.ev_fnan, NULL, prec, NULL))
  yre <- result$re; yim <- result$im
  if (is.null(yre) && is.null(yim)) return(.ev_res(.ev_fone, NULL, prec, NULL))
  ysize <- .ev_fastlog(yre)
  if (ysize > 5) {
    prec <- prec + ysize
    r <- .ev_evalf(ex, prec, options)
    yre <- r$re; yim <- r$im
  }
  if (base$type == "Constant" && base$name == "E") {
    if (!is.null(yim)) {
      z <- .ev_mpc_down(.mpc_exp(list(re = if (is.null(yre)) .mp_zero else yre, im = yim), prec + 20), prec)
      return(.ev_finalize_complex(z$re, z$im, target_prec))
    }
    return(.ev_res(.ev_exp_sp(yre, target_prec), NULL, target_prec, NULL))
  }
  r <- .ev_evalf(base, prec + 5, options)
  xre <- r$re; xim <- r$im
  if (is.null(xre) && is.null(xim)) {
    if (!is.null(yim)) return(.ev_res(.ev_fnan, NULL, prec, NULL))
    if (.ev_sign(yre) == 1L) return("zoo")
    return(.ev_none())
  }
  z0 <- function(a) if (is.null(a)) .mp_zero else a
  if (!is.null(yim)) {
    z <- .ev_mpc_down(.mpc_pow(list(re = z0(xre), im = z0(xim)), list(re = z0(yre), im = yim), target_prec + 20), target_prec)
    return(.ev_finalize_complex(z$re, z$im, target_prec))
  }
  if (!is.null(xim) || .ev_cmp(xre, .mp_zero) < 0) {
    z <- .ev_mpc_down(.mpc_pow_mpf(list(re = z0(xre), im = z0(xim)), yre, target_prec + 20), target_prec)
    return(.ev_finalize_complex(z$re, z$im, target_prec))
  }
  .ev_res(.ev_pow_real(xre, yre, target_prec), NULL, target_prec, NULL)
}
# the libmpf functions on special values
.ev_pow_int_sp <- function(a, n, p) {
  if (.ev_sp(a)) {
    if (a$sp == "nan") return(.ev_fnan)
    if (n < 0) return(.mp_zero)
    if (a$sp == "-inf" && gmp::mod.bigz(n, 2) == 0) return(.ev_finf)
    return(a)
  }
  if (a$m == 0 && n < 0) return(.ev_finf)
  .ev_pow_int(a, n, p)
}
.ev_sqrt_sp <- function(a, p) if (.ev_sp(a)) a else .ev_sqrt(a, p)
.ev_exp_sp <- function(a, p) {
  if (.ev_sp(a)) return(switch(a$sp, inf = .ev_finf, `-inf` = .mp_zero, .ev_fnan))
  .ev_exp(a, p)
}
.ev_pow_real <- function(x, y, p) {
  if (.ev_sp(x) || .ev_sp(y)) {
    xs <- if (.ev_sp(x)) x$sp else "num"; ys <- if (.ev_sp(y)) y$sp else "num"
    if (xs == "nan" || ys == "nan") return(.ev_fnan)
    if (ys == "inf") return(if (.ev_cmp(x, .ev_fone) > 0) .ev_finf else .mp_zero)
    if (ys == "-inf") return(if (.ev_cmp(x, .ev_fone) > 0) .mp_zero else .ev_finf)
    return(if (y$m > 0) .ev_finf else .mp_zero)
  }
  .ev_down(.mp_pow_real(x, y, p + 20), p)
}

# ---------------------------------------------------------------------------
# Functions

.ev_evalf_exp <- function(expr, prec, options)
  .ev_evalf_pow(new_sym("Pow", list(S_E, expr$args[[1L]])), prec, options)

.ev_evalf_trig <- function(v, prec, options) {
  func <- switch(v$name, cos = .ev_cos, sin = .ev_sin, tan = .ev_tan, .ev_nie())
  arg <- v$args[[1L]]
  xprec <- prec + 20
  r <- .ev_evalf(arg, xprec, options)
  if (.ev_is_zoo(r)) .ev_nie()
  re <- r$re; im <- r$im
  if (!is.null(im)) {
    if (!is.null(options$subs)) v <- .ev_subs_all(v, options$subs)
    xe <- .ev_eval_evalf(v, prec)
    if (is.null(xe)) .ev_nie()
    return(.ev_evalf(xe, prec, options))
  }
  if (is.null(re)) {
    if (v$name == "cos") return(.ev_res(.ev_fone, NULL, prec, NULL))
    return(.ev_none())
  }
  if (.ev_sp(re)) .ev_nie()
  xsize <- .ev_fastlog(re)
  if (xsize < 1) return(.ev_res(func(re, prec), NULL, prec, NULL))
  if (xsize >= 10) {
    xprec <- prec + xsize
    r <- .ev_evalf(arg, xprec, options); re <- r$re
  }
  maxprec <- if (is.null(options$maxprec)) .ev_DEFAULT_MAXPREC else options$maxprec
  repeat {
    y <- func(re, prec)
    ysize <- .ev_fastlog(y)
    gap <- -ysize
    accuracy <- (xprec - xsize) - gap
    if (accuracy < prec) {
      if (isTRUE(options$verbose)) cat("SIN/COS/TAN", accuracy, "wanted", prec, "gap", gap, "\n")
      if (xprec > maxprec) return(.ev_res(y, NULL, accuracy, NULL))
      xprec <- xprec + gap
      r <- .ev_evalf(arg, xprec, options); re <- r$re
      next
    }
    return(.ev_res(y, NULL, prec, NULL))
  }
}

.ev_evalf_log <- function(expr, prec, options) {
  if (length(expr$args) > 1L) {
    # log(x, b).doit(): the evaluated log, or log(x)/log(b)
    v <- fcall("log", expr$args)
    if (v$type == "Function" && v$name == "log" && length(v$args) > 1L)
      v <- new_sym("Mul", list(fnode("log", list(v$args[[1L]])), new_sym("Pow", list(fnode("log", list(v$args[[2L]])), S_NegativeOne))))
    return(.ev_evalf(v, prec, options))
  }
  arg <- expr$args[[1L]]
  workprec <- prec + 10
  result <- .ev_evalf(arg, workprec, options)
  if (.ev_is_zoo(result)) return(result)
  xre <- result$re; xim <- result$im
  if (is.null(xre) && is.null(xim)) xre <- .mp_zero
  if (!is.null(xim)) {
    re <- .ev_evalf_log(fnode("log", list(fnode("Abs", list(arg)))), prec, options)
    im <- .ev_atan2(xim, if (is.null(xre)) .mp_zero else xre, prec)
    return(.ev_res(re$re, im, re$re_acc, prec))
  }
  if (.ev_sp(xre)) {
    re <- switch(xre$sp, nan = .ev_fnan, .ev_finf)
    if (xre$sp == "-inf") return(.ev_res(.ev_finf, .ev_pi(prec), prec, prec))
    return(.ev_res(re, NULL, prec, NULL))
  }
  imaginary_term <- .ev_cmp(xre, .mp_zero) < 0
  re <- .ev_ln(.ev_abs(xre), prec)
  size <- .ev_fastlog(re)
  if (prec - size > workprec) {
    add <- new_sym("Add", list(S_NegativeOne, arg))
    r <- .ev_evalf_add(add, prec, options)
    xre <- r$re; xre_acc <- r$re_acc
    if (!.ev_is_fzero(xre) && !is.null(xre) && (is.null(xre_acc) || xre_acc > 1)) {
      prec2 <- workprec - .ev_fastlog(xre)
      re <- .ev_ln(.ev_abs(.ev_add(xre, .ev_fone, prec2)), prec)
    }
  }
  re_acc <- prec
  if (imaginary_term) return(.ev_res(re, .ev_pi(prec), re_acc, prec))
  .ev_res(re, NULL, re_acc, NULL)
}

.ev_evalf_atan <- function(v, prec, options) {
  arg <- v$args[[1L]]
  result <- .ev_evalf(arg, prec + 5, options)
  if (.ev_is_zoo(result)) return(result)
  if (is.null(result$re) && is.null(result$im)) return(.ev_none())
  if (!is.null(result$im)) .ev_nie()
  xre <- result$re
  if (.ev_sp(xre)) {
    if (xre$sp == "nan") return(.ev_res(.ev_fnan, NULL, prec, NULL))
    h <- .mp_ldexp(.mp_pi(prec), -1)
    return(.ev_res(if (xre$sp == "inf") h else .ev_neg(h), NULL, prec, NULL))
  }
  .ev_res(.mp_atan(xre, prec), NULL, prec, NULL)
}

# evalf_subs: the Floats of subs at precision prec
.ev_evalf_subs <- function(prec, subs) {
  lapply(subs, function(b) {
    b <- as_sym(b)
    if (b$type == "Float") .ev_float(.float_mp(b), prec) else b
  })
}
.ev_subs_all <- function(e, subs) {
  if (!length(subs)) return(e)
  pairs <- lapply(names(subs), function(nm) list(Symbol(nm), as_sym(subs[[nm]])))
  subs(e, pairs)
}

.ev_evalf_piecewise <- function(expr, prec, options) {
  if (!is.null(options$subs)) {
    expr <- .ev_subs_all(expr, .ev_evalf_subs(prec, options$subs))
    newopts <- .ev_copy_options(options)
    rm("subs", envir = newopts)
    return(.ev_evalf(expr, prec, newopts))
  }
  .ev_nie()
}

.ev_evalf_alg_num <- function(a, prec, options) .ev_evalf(.an_value(a), prec, options)

.ev_evalf_symbol <- function(x, prec, options) {
  val <- options$subs[[x$name]]
  if (is.null(val)) stop(.ev_condition("KeyError", x$name))
  if (is.null(options$`_cache`)) options$`_cache` <- list()
  cached <- options$`_cache`[[x$name]]
  if (!is.null(cached) && cached[[2L]] >= prec) return(cached[[1L]])
  v <- .ev_evalf(as_sym(val), prec, options)
  cache <- options$`_cache`
  cache[[x$name]] <- list(v, prec)
  options$`_cache` <- cache
  v
}

# ---------------------------------------------------------------------------
# Dispatch

.ev_options <- function(...) {
  o <- new.env(parent = emptyenv())
  a <- list(...)
  for (nm in names(a)) assign(nm, a[[nm]], envir = o)
  o
}
.ev_copy_options <- function(options) {
  o <- new.env(parent = emptyenv())
  for (nm in ls(options, all.names = TRUE)) assign(nm, get(nm, envir = options), envir = o)
  o
}

.ev_table <- function(x) {
  t <- x$type
  switch(EXPR = t,
    Symbol = , Dummy = .ev_evalf_symbol,
    Float = function(x, prec, options) {
      f <- .float_mp(x)
      if (f$m == 0) return(.ev_res(.mp_zero, NULL, prec, NULL))
      .ev_res(f, NULL, prec, NULL)
    },
    Rational = function(x, prec, options) {
      q <- gmp::as.bigq(x$value)
      if (q == gmp::as.bigq(1, 2)) return(.ev_res(.ev_fhalf, NULL, prec, NULL))
      .ev_res(.ev_from_rational(q, prec), NULL, prec, NULL)
    },
    Integer = function(x, prec, options) {
      n <- gmp::numerator(gmp::as.bigq(x$value))
      if (n == 0) return(.ev_res(NULL, NULL, prec, NULL))
      if (n == 1) return(.ev_res(.ev_fone, NULL, prec, NULL))
      if (n == -1) return(.ev_res(.ev_fnone, NULL, prec, NULL))
      .ev_res(.ev_from_int(n, prec), NULL, prec, NULL)
    },
    Constant = switch(EXPR = x$name,
      pi = function(x, prec, options) .ev_res(.ev_pi(prec), NULL, prec, NULL),
      E = function(x, prec, options) .ev_res(.ev_e(prec), NULL, prec, NULL),
      I = function(x, prec, options) .ev_res(NULL, .ev_fone, NULL, prec),
      zoo = function(x, prec, options) "zoo",
      nan = function(x, prec, options) .ev_res(.ev_fnan, NULL, prec, NULL),
      NULL),
    Add = .ev_evalf_add,
    Mul = .ev_evalf_mul,
    Pow = .ev_evalf_pow,
    Function = switch(EXPR = x$name,
      exp = .ev_evalf_exp,
      cos = , sin = , tan = .ev_evalf_trig,
      log = .ev_evalf_log,
      atan = .ev_evalf_atan,
      Abs = function(expr, prec, options) .ev_get_abs(expr$args[[1L]], prec, options),
      re = function(expr, prec, options) .ev_get_complex_part(expr$args[[1L]], 0, prec, options),
      im = function(expr, prec, options) .ev_get_complex_part(expr$args[[1L]], 1, prec, options),
      floor = function(expr, prec, options) .ev_get_integer_part(expr$args[[1L]], -1, options),
      ceiling = function(expr, prec, options) .ev_get_integer_part(expr$args[[1L]], 1, options),
      NULL),
    Integral = .ev_evalf_integral,
    Sum = .ev_evalf_sum,
    Product = .ev_evalf_prod,
    Piecewise = .ev_evalf_piecewise,
    AlgebraicNumber = .ev_evalf_alg_num,
    NULL)
}

# evalf(x, prec, options)
.ev_evalf <- function(x, prec, options) {
  rf <- .ev_table(x)
  r <- if (is.null(rf)) .ev_fallback(x, prec, options) else
    tryCatch(rf(x, prec, options), KeyError = function(err) .ev_fallback(x, prec, options))
  if (isTRUE(options$verbose)) {
    cat("### input", sstr(x), "\n")
    cat("### output", if (.ev_is_zoo(r)) "zoo" else if (is.null(r$re)) "0" else .mp_to_str(.ev_norm(r$re, 200), 50), "\n\n")
  }
  chop <- options$chop
  if (!is.null(chop) && !isFALSE(chop)) {
    chop_prec <- if (isTRUE(chop)) prec else {
      cp <- as.integer(round(-3.321 * log10(as.numeric(chop)) + 2.5))
      if (cp == 3L) cp - 1L else cp
    }
    r <- .ev_chop_parts(r, chop_prec)
  }
  if (isTRUE(options$strict)) .ev_check_target(x, r, prec)
  r
}

# The fallback of evalf for classes without a table entry: _eval_evalf at
# prec, then the real and imaginary parts converted at prec
.ev_fallback <- function(x, prec, options) {
  if (!is.null(options$subs)) x <- .ev_subs_all(x, .ev_evalf_subs(prec, options$subs))
  xe <- .ev_eval_evalf(x, prec)
  if (is.null(xe)) .ev_nie()
  if (!is_sym(xe)) .ev_nie()
  if (xe$type %in% .ev_set_types()) .ev_nie()
  ri <- .ev_float_parts(xe, prec)
  if (is.null(ri)) .ev_nie()
  .ev_res(ri$re, ri$im, if (is.null(ri$re)) NULL else prec, if (is.null(ri$im)) NULL else prec)
}
.ev_set_types <- function() c("FiniteSet", "Interval", "Union", "Intersection", "Complement", "EmptySet",
                               "Tuple", .REL_TYPES, .BOOL_TYPES)

# The parts of a number Float, Float*I or Float + Float*I (as_real_imag of
# the value of _eval_evalf); NULL when it is not such a number. With prec,
# any other number part is converted by _to_mpmath(prec), as evalf's
# fallback does (pi in log(-1) + 0.1, say)
.ev_float_parts <- function(v, prec = NULL) {
  part <- function(z) {
    if (z$type == "Integer" && z$value == 0) return(list(ok = TRUE, v = NULL))
    if (z$type == "Float") { f <- .float_mp(z); return(list(ok = TRUE, v = if (f$m == 0) NULL else f)) }
    if (is_exact_number(z)) return(list(ok = TRUE, v = .mp_from_q(z$value, 53)))
    if (z$type == "Constant" && z$name %in% c("oo", "-oo", "nan"))
      return(list(ok = TRUE, v = switch(z$name, oo = .ev_finf, `-oo` = .ev_fninf, .ev_fnan)))
    if (!is.null(prec) && .expr_is_number(z)) {
      q <- tryCatch(.ev_evalf(z, prec, .ev_options()), error = function(e) NULL)
      if (!is.null(q) && is.null(q$im) && !.ev_is_zoo(q)) return(list(ok = TRUE, v = q$re))
    }
    list(ok = FALSE)
  }
  re <- S_Zero; im <- S_Zero
  if (!is.null(prec)) {
    # xe.as_real_imag(), refused when re() or im() stays unevaluated
    re <- tryCatch(re(v), error = function(e) NULL); im <- tryCatch(im(v), error = function(e) NULL)
    if (is.null(re) || is.null(im) || .has_function(re, "re") || .has_function(im, "im")) return(NULL)
  } else if (v$type == "Constant" && v$name == "I") im <- S_One
  else if (v$type == "Add") {
    if (length(v$args) > 2L) return(NULL)
    for (a in v$args) {
      pc <- .ev_pure_complex(a)
      if (!is.null(pc) && pc[[1L]]$type == "Integer" && pc[[1L]]$value == 0) im <- pc[[2L]] else re <- a
    }
  } else if (v$type == "Mul") {
    pc <- .ev_pure_complex(v)
    if (is.null(pc)) return(NULL)
    im <- pc[[2L]]
  } else re <- v
  pr <- part(re); pim <- part(im)
  if (!pr$ok || !pim$ok) return(NULL)
  list(re = pr$v, im = pim$v)
}

# ---------------------------------------------------------------------------
# Integrals, sums and products

.ev_as_mpmath <- function(x, prec, options) {
  x <- as_sym(x)
  if ((x$type == "Integer" && x$value == 0) || (x$type == "Float" && x$value == 0)) return(.mv_mpf(.mp_zero))
  if (x$type == "Constant" && x$name == "oo") return(.mv_inf())
  if (x$type == "Constant" && x$name == "-oo") return(.mv_ninf())
  .ev_quad_to_mpmath(.ev_evalf(x, prec, options))
}
.ev_quad_to_mpmath <- function(q) {
  if (.ev_is_zoo(q)) .ev_nie()
  sp <- function(a) if (.ev_sp(a)) .mv_special(a$sp) else NULL
  if (!is.null(q$im)) return(.mv_mpc(if (is.null(q$re)) .mp_zero else q$re, q$im))
  if (!is.null(q$re)) return(if (.ev_sp(q$re)) sp(q$re) else .mv_mpf(q$re))
  .mv_mpf(.mp_zero)
}

.ev_do_integral <- function(expr, prec, options) {
  func <- expr$args[[1L]]
  lim <- expr$args[[2L]]$args
  x <- lim[[1L]]; xlow <- lim[[2L]]; xhigh <- lim[[3L]]
  if (identical(xlow$key, xhigh$key)) {
    xlow <- xhigh <- S_Zero
  } else if (!(x$name %in% vapply(free_symbols(func), `[[`, "", "name"))) {
    fl <- vapply(free_symbols(xlow), `[[`, "", "name"); fh <- vapply(free_symbols(xhigh), `[[`, "", "name")
    if (length(intersect(fl, fh))) {
      d <- .sub(xhigh, xlow)
      if (!length(free_symbols(d))) { xlow <- S_Zero; xhigh <- d }
    }
  }
  oldmaxprec <- if (is.null(options$maxprec)) .ev_DEFAULT_MAXPREC else options$maxprec
  options$maxprec <- min(oldmaxprec, 2 * prec)
  old <- .MPX$prec
  on.exit(.MPX$prec <- old, add = TRUE)
  .MPX$prec <- prec + 5
  lo <- .ev_as_mpmath(xlow, prec + 15, options)
  hi <- .ev_as_mpmath(xhigh, prec + 15, options)
  have_part <- c(FALSE, FALSE)
  max_real_term <- -Inf; max_imag_term <- -Inf
  xname <- x$name
  f <- function(t) {
    sub <- list(.ev_mv_sym(t)); names(sub) <- xname
    r <- .ev_evalf(func, .MPX$prec, .ev_options(subs = sub))
    if (.ev_is_zoo(r)) .ev_nie()
    have_part[1L] <<- have_part[1L] || !is.null(r$re)
    have_part[2L] <<- have_part[2L] || !is.null(r$im)
    max_real_term <<- max(max_real_term, .ev_fastlog(r$re))
    max_imag_term <<- max(max_imag_term, .ev_fastlog(r$im))
    z0 <- function(a) if (is.null(a)) .mp_zero else a
    if (!is.null(r$im)) return(.mv_mpc(z0(r$re), r$im))
    .mv_mpf(z0(r$re))
  }
  if (identical(options$quad, "osc")) {
    A <- .ev_osc_frequency(func, x)
    if (is.null(A))
      .ev_valueerror("An integrand of the form sin(A*x+B)*f(x) or cos(A*x+B)*f(x) is required for oscillatory quadrature")
    period <- .ev_as_mpmath(.div(.mul(list(S_Two, S_Pi)), A), prec + 15, options)
    result <- .mv_as_mpc(.ev_quadosc(f, lo, hi, period))
    # quadosc does not estimate its error
    quadrature_error <- -Inf
  } else {
    qr <- .mp_quad(f, list(lo, hi), error = TRUE)
    result <- .mv_as_mpc(qr$value)
    quadrature_error <- .ev_fastlog(.mv_convert(qr$err)$re)
  }
  .MPX$prec <- old
  options$maxprec <- oldmaxprec
  re <- im <- re_acc <- im_acc <- NULL
  if (have_part[1L]) {
    re <- result$re
    if (re$m == 0) {
      sz <- .ev_scaled_zero(as.integer(-max(prec, max_real_term, quadrature_error)))
      re_acc <- sz[[2L]]; re <- .ev_scaled_zero(sz[[1L]])
    } else re_acc <- as.integer(-max(max_real_term - .ev_fastlog(re) - prec, quadrature_error))
  }
  if (have_part[2L]) {
    im <- result$im
    if (im$m == 0) {
      sz <- .ev_scaled_zero(as.integer(-max(prec, max_imag_term, quadrature_error)))
      im_acc <- sz[[2L]]; im <- .ev_scaled_zero(sz[[1L]])
    } else im_acc <- as.integer(-max(max_imag_term - .ev_fastlog(im) - prec, quadrature_error))
  }
  .ev_res(re, im, re_acc, im_acc)
}
# A of a factor sin(A*x + B) or cos(A*x + B) of the integrand
# (func.match(cos(A*x + B)*D), then sin), NULL when there is none
.ev_osc_frequency <- function(func, x) {
  fs <- if (func$type == "Mul") func$args else list(func)
  for (nm in c("cos", "sin")) for (t in fs) {
    if (t$type == "Function" && t$name == nm) {
      A <- diff(t$args[[1L]], x)
      if (!.has1(A, x) && !is_zero_num(A)) return(A)
    }
  }
  NULL
}
# mpmath's quadosc with a period: the integral to the first zero n*period/2
# and the sum of the integrals between consecutive zeros by nsum
.ev_quadosc <- function(f, a, b, period) {
  isinf <- function(v, s) { v <- .mv_convert(v); identical(v$sp, s) }
  if (isinf(a, "-inf") && isinf(b, "inf"))
    return(.hx_add(.ev_quadosc(f, a, .mv_int(0), period), .ev_quadosc(f, .mv_int(0), b, period)))
  if (isinf(a, "-inf"))
    return(.ev_quadosc(function(t) f(.hx_neg(t)), .hx_neg(b), .mv_inf(), period))
  if (!isinf(b, "inf")) .ev_valueerror("quadosc requires an infinite integration interval")
  zeros <- function(n) .hx_div(.hx_mul(n, period), .mv_int(2))
  s <- .mp_quad(f, list(a, zeros(.mv_int(1))))
  term <- function(k) .mp_quad(f, list(zeros(k), zeros(.hx_add(k, .mv_int(1)))))
  .hx_add(s, .mp_nsum(term, list(.mv_int(1), .mv_inf())))
}

# an mpmath value as a sym (a quadrature node substituted for the variable)
.ev_mv_sym <- function(t) {
  t <- .mv_convert(t)
  if (t$k == "mpc") return(.mpc_sym(.mpc(t$re, t$im), .MPX$prec))
  .ev_float(t$re, .MPX$prec)
}

.ev_evalf_integral <- function(expr, prec, options) {
  lims <- expr$args[-1L]
  if (length(lims) != 1L || length(lims[[1L]]$args) != 3L) .ev_nie()
  workprec <- prec
  i <- 0
  maxprec <- if (is.null(options$maxprec)) Inf else options$maxprec
  repeat {
    result <- .ev_do_integral(expr, workprec, options)
    accuracy <- .ev_complex_accuracy(result)
    if (accuracy >= prec) break
    if (workprec >= maxprec) break
    if (accuracy == -1) workprec <- workprec * 2 else workprec <- workprec + max(prec, 2^i)
    workprec <- min(workprec, maxprec)
    i <- i + 1
  }
  result
}
