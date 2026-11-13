# Elementary functions, after sympy/functions/elementary.

.real_arg <- function(afs) isTRUE(afs[[1L]][["real"]])

# ---------------------------------------------------------------------------
# exp and log

def_function(
  "exp",
  eval = function(args) {
    a <- args[[1L]]
    if (is_zero_num(a)) return(S_One)
    if (is_one_num(a)) return(S_E)
    if (is_const(a, "oo")) return(S_Infinity)
    if (is_const(a, "-oo")) return(S_Zero)
    if (is_const(a, "zoo")) return(S_NaN)
    if (a$type == "Float") return(.float_eval("exp", args))
    if (a$type == "Function" && a$name == "log") return(a$args[[1L]])
    if (a$type == "Mul") {
      # exp.eval: rational multiples of pi*I, reduced to (-pi, pi]
      coeff <- as_coefficient(a, .mul(list(S_Pi, S_I)))
      if (!is.null(coeff) && !is_zero_num(coeff)) {
        if (isTRUE(ask_fact(.mul(list(S_Two, coeff)), "integer"))) {
          if (isTRUE(ask_fact(coeff, "even"))) return(S_One)
          if (isTRUE(ask_fact(coeff, "odd"))) return(S_NegativeOne)
          ch <- .add(list(coeff, S_Half))
          if (isTRUE(ask_fact(ch, "even"))) return(.neg(S_I))
          if (isTRUE(ask_fact(ch, "odd"))) return(S_I)
        } else if (coeff$type %in% c("Integer", "Rational")) {
          q <- coeff$value
          nq <- q - 2 * gmp::as.bigq(floor(as.double(q / 2)))
          while (nq < 0) nq <- nq + 2
          while (nq >= 2) nq <- nq - 2
          if (nq > 1) nq <- nq - 2
          if (nq != q) return(fcall("exp", list(.mul(list(make_rational(nq), S_Pi, S_I)))))
        }
      }
      cm <- as_coeff_mul_split(a); co <- cm[[1L]]; terms <- cm[[2L]]
      # oo and -oo are the Number coefficient of a SymPy Mul
      io <- which(vapply(a$args, function(t) is_const(t, "oo") || is_const(t, "-oo"), NA))
      if (length(io)) {
        co <- a$args[[io[1L]]]; rest <- a$args[-io[1L]]
        terms <- if (length(rest) == 1L) rest[[1L]] else .mul(rest)
      }
      if (is_const(co, "oo") || is_const(co, "-oo")) {
        if (!length(free_symbols(terms))) {
          if (is_const(co, "-oo")) terms <- .neg(terms)
          ri <- .re_im(terms)
          if (is.null(ri)) return(NULL)
          if (is_zero_num(ri[[1L]]) && !is_zero_num(terms)) return(S_NaN)
          sg <- if (isTRUE(ask_fact(ri[[1L]], "positive"))) 1 else if (isTRUE(ask_fact(ri[[1L]], "negative"))) -1
                else .mp_num_sign(ri[[1L]])
          if (isTRUE(sg > 0) && !is_zero_num(ri[[2L]])) return(S_ComplexInfinity)
          if (isTRUE(sg < 0)) return(S_Zero)
        }
        return(NULL)
      }
      coeffs <- list(co); log_term <- NULL
      for (term in .margs(terms)) {
        tl <- logcombine(term)
        if (tl$type == "Function" && tl$name == "log") {
          if (is.null(log_term)) log_term <- tl$args[[1L]] else return(NULL)
        } else if (.is_comparable(term)) coeffs[[length(coeffs) + 1L]] <- term
        else return(NULL)
      }
      return(if (is.null(log_term)) NULL else .pow(log_term, .mul(coeffs)))
    }
    if (a$type == "Add") {
      out <- list(); add <- list(); changed <- FALSE
      for (t in a$args) {
        if (is_one_num(t)) { add[[length(add) + 1L]] <- t; next }
        nt <- fcall("exp", list(t))
        if (nt$type == "Function" && nt$name == "exp") {
          if (!identical(nt$args[[1L]]$key, t$key)) { add[[length(add) + 1L]] <- nt$args[[1L]]; changed <- TRUE }
          else add[[length(add) + 1L]] <- t
        } else out[[length(out) + 1L]] <- nt
      }
      # Mul.flatten rebuilds exp factors, so only a lone exp stays unevaluated
      if (length(out)) return(.mul(c(out, list(fcall("exp", list(.add(add)))))))
      if (changed) return(fnode("exp", list(.add(add))))
    }
    if (isTRUE(ask_fact(a, "zero"))) return(S_One)
    NULL
  },
  diff = function(args, i) fcall("exp", args),
  numeric = exp,
  facts = function(args, afs) {
    f <- afs[[1L]]
    out <- c(commutative = !identical(f[["commutative"]], FALSE))
    if (isTRUE(f[["real"]])) out["positive"] <- TRUE
    if (isTRUE(f[["finite"]])) {
      out["finite"] <- TRUE
      out["zero"] <- FALSE
    }
    if (isTRUE(f[["rational"]]) && isTRUE(f[["nonzero"]])) out["transcendental"] <- TRUE
    out
  },
  latex = function(a) paste0("e^{", a[1], "}")
)

.log_exact <- function(x, b) {
  if (!is_exact_number(x) || !is_exact_number(b)) return(NULL)
  if (x$value <= 0 || b$value <= 0 || b$value == 1) return(NULL)
  lx <- log(abs(as.double(x$value)))
  lb <- log(abs(as.double(b$value)))
  if (!is.finite(lx) || !is.finite(lb)) return(NULL)
  k <- round(lx / lb)
  if (abs(k) > 10000) return(NULL)
  if (.bigq_pow_int(b$value, k) == x$value) return(.int(k))
  NULL
}

def_function(
  "log", nargs = 1:2,
  eval = function(args) {
    if (length(args) == 2L) {
      x <- args[[1L]]
      b <- args[[2L]]
      ex <- .log_exact(x, b)
      if (!is.null(ex)) return(ex)
      if (is_const(b, "E")) return(fcall("log", list(x)))
      return(.div(fcall("log", list(x)), fcall("log", list(b))))
    }
    a <- args[[1L]]
    t <- a$type
    if (t == "Function" && a$name == "exp_polar") return(unpolarify(a$args[[1L]]))
    if (t == "Integer" || t == "Rational") {
      if (a$value == 1) return(S_Zero)
      if (a$value == 0) return(S_ComplexInfinity)
      if (a$value < 0) return(.add(list(fcall("log", list(num_neg(a))),
                                        .mul(list(S_I, S_Pi)))))
      if (gmp::numerator(a$value) == 1)
        return(.neg(fcall("log", list(make_rational(gmp::as.bigq(gmp::denominator(a$value)))))))
      return(NULL)
    }
    if (t == "Float") {
      # log.eval: zero gives zoo, a negative number log(-a) + I*pi
      if (a$value == 0) return(S_ComplexInfinity)
      if (a$value < 0) return(.add(list(.mul(list(S_Pi, S_I)), fcall("log", list(num_neg(a))))))
      if (.float_prec(a) != 53 || !is.null(a$extra$mp)) return(.float_eval("log", args))
      if (a$value > 0) return(Float(log(a$value)))
      return(.complex_to_sym(log(as.complex(a$value))))
    }
    if (is_const(a, "E")) return(S_One)
    if (is_const(a, "zoo")) return(S_ComplexInfinity)
    if (is_const(a, "oo") || is_const(a, "-oo")) return(S_Infinity)
    if (is_const(a, "I")) return(.mul(list(S_Half, S_I, S_Pi)))
    # SymPy: a negative number gives log(-a) + I*pi
    if (t %in% c("Mul", "Add", "Pow", "Constant", "Function") && !length(free_symbols(a)) && isTRUE(ask_fact(a, "negative")))
      return(.add(list(.mul(list(S_Pi, S_I)), fcall("log", list(.neg(a))))))
    if (t == "Function" && a$name == "exp" && isTRUE(ask_fact(a$args[[1L]], "real")))
      return(a$args[[1L]])
    ic <- .imag_coeff(a)
    if (!is.null(ic) && is_exact_number(ic)) {
      if (ic$value > 0) return(.add(list(fcall("log", list(ic)), .mul(list(S_Half, S_I, S_Pi)))))
      return(.add(list(fcall("log", list(num_neg(ic))), .mul(list(.rat(-1, 2), S_I, S_Pi)))))
    }
    if (.has1(a, S_I) && !length(free_symbols(a)) && isTRUE(ask_fact(a, "algebraic"))) {
      r <- .log_algebraic(a)
      if (!is.null(r)) return(r)
    }
    NULL
  },
  diff = function(args, i) {
    if (length(args) == 2L) {
      x <- args[[1L]]
      b <- args[[2L]]
      if (i == 1L) return(.inv(.mul(list(x, fcall("log", list(b))))))
      return(.neg(.div(fcall("log", list(x)),
                       .mul(list(b, .sq(fcall("log", list(b))))))))
    }
    .inv(args[[1L]])
  },
  numeric = function(x, b) {
    if (missing(b)) {
      if (is.complex(x) || x < 0) return(log(as.complex(x)))
      return(log(x))
    }
    log(as.complex(x)) / log(as.complex(b))
  },
  facts = function(args, afs) {
    f <- afs[[1L]]
    out <- c(commutative = TRUE)
    if (isTRUE(f[["positive"]])) {
      out["real"] <- TRUE
      if (isTRUE(f[["finite"]])) out["finite"] <- TRUE
    }
    if (isTRUE(f[["negative"]])) out["real"] <- FALSE
    # SymPy: the sign of log(a) is the sign of a - 1
    if (length(args) == 1L && isTRUE(f[["positive"]])) {
      d <- .sub(args[[1L]], S_One)
      if (isTRUE(ask_fact(d, "positive"))) out["positive"] <- TRUE
      else if (isTRUE(ask_fact(d, "negative"))) out["negative"] <- TRUE
      else if (isTRUE(ask_fact(d, "zero"))) out["zero"] <- TRUE
      else if (isTRUE(ask_fact(d, "nonnegative"))) out["nonnegative"] <- TRUE
    }
    if (is_exact_number(args[[1L]]) && args[[1L]]$value > 0 &&
        !is_one_num(args[[1L]])) out["transcendental"] <- TRUE
    out
  },
  latex = function(a) if (length(a) == 2L) paste0("\\log_{", a[2], "}{\\left(", a[1], " \\right)}")
                      else paste0("\\log{\\left(", a[1], " \\right)}")
)

# log.eval for an algebraic number r + I*i: log of the modulus plus I times
# the angle from _log_atan_table
.log_atan_table <- function() {
  if (!is.null(.mi_state$log_atan)) return(.mi_state$log_atan)
  s <- function(n) sqrt(as_sym(n)); P <- S_Pi; R <- function(p, q) Rational(p, q)
  tb <- list(
    list(s(3), P / 3), list(S_One, P / 4), list(sqrt(5 - 2 * s(5)), P / 5),
    list(s(2) * sqrt(5 - s(5)) / (1 + s(5)), P / 5),
    list(sqrt(5 + 2 * s(5)), P * R(2, 5)),
    list(s(2) * sqrt(s(5) + 5) / (-1 + s(5)), P * R(2, 5)),
    list(s(3) / 3, P / 6), list(s(2) - 1, P / 8),
    list(sqrt(2 - s(2)) / sqrt(s(2) + 2), P / 8),
    list(s(2) + 1, P * R(3, 8)),
    list(sqrt(s(2) + 2) / sqrt(2 - s(2)), P * R(3, 8)),
    list(sqrt(1 - 2 * s(5) / 5), P / 10),
    list((-s(2) + s(10)) / (2 * sqrt(s(5) + 5)), P / 10),
    list(sqrt(1 + 2 * s(5) / 5), P * R(3, 10)),
    list((s(2) + s(10)) / (2 * sqrt(5 - s(5))), P * R(3, 10)),
    list(2 - s(3), P / 12), list((-1 + s(3)) / (1 + s(3)), P / 12),
    list(2 + s(3), P * R(5, 12)), list((1 + s(3)) / (-1 + s(3)), P * R(5, 12)))
  out <- lapply(tb, `[[`, 2L); names(out) <- vapply(tb, function(p) as_sym(p[[1L]])$key, "")
  .mi_state$log_atan <- out
  out
}
.log_algebraic <- function(a) {
  sp <- .as_coeff_mul_dep(a, S_I)
  coeff <- sp[[1L]]; arg_ <- .mul(sp[[2L]])
  if (isTRUE(ask_fact(coeff, "negative"))) { coeff <- .neg(coeff); arg_ <- .neg(arg_) }
  arg_ <- expand_mul(arg_, deep = FALSE)
  ts <- if (arg_$type == "Add") arg_$args else list(arg_)
  dep <- vapply(ts, function(z) .has1(z, S_I), NA)
  r_ <- .add(ts[!dep]); i_ <- extract_multiplicatively(.add(ts[dep]), S_I)
  if (is.null(i_) || .has1(i_, S_I)) return(NULL)
  if (!(isTRUE(ask_fact(coeff, "real")) && !is_zero_num(i_) && isTRUE(ask_fact(i_, "real")) && isTRUE(ask_fact(r_, "real")))) return(NULL)
  if (is_zero_num(r_)) {
    if (isTRUE(ask_fact(i_, "positive"))) return(.add(list(.mul(list(S_Pi, S_I, S_Half)), fcall("log", list(.mul(list(coeff, i_)))))))
    if (isTRUE(ask_fact(i_, "negative"))) return(.add(list(.mul(list(.neg(S_Pi), S_I, S_Half)), fcall("log", list(.mul(list(coeff, .neg(i_))))))))
    return(NULL)
  }
  t <- cancel(.div(i_, r_)); t1 <- cancel(.neg(t))
  tab <- .log_atan_table()
  modulus <- function() ratsimp(.mul(list(coeff, fcall("Abs", list(arg_)))))
  if (!is.null(tab[[t$key]])) {
    ang <- tab[[t$key]]
    if (isTRUE(ask_fact(r_, "positive"))) return(.add(list(fcall("log", list(modulus())), .mul(list(S_I, ang)))))
    return(.add(list(fcall("log", list(modulus())), .mul(list(S_I, .sub(ang, S_Pi))))))
  }
  if (!is.null(tab[[t1$key]])) {
    ang <- tab[[t1$key]]
    if (isTRUE(ask_fact(r_, "positive"))) return(.add(list(fcall("log", list(modulus())), .mul(list(S_I, .neg(ang))))))
    return(.add(list(fcall("log", list(modulus())), .mul(list(S_I, .sub(S_Pi, ang))))))
  }
  NULL
}

# ---------------------------------------------------------------------------
# Trigonometric functions

.cos_table <- function(r) {
  # cos(r*pi) for rational r in [0, 1/2]
  key <- as.character(r)
  s <- function(n) .sqrt(.int(n))
  switch(key,
         "0" = S_One,
         "1/2" = S_Zero,
         "1/3" = S_Half,
         "1/4" = .mul(list(S_Half, s(2))),
         "1/6" = .mul(list(S_Half, s(3))),
         "1/5" = .mul(list(.rat(1, 4), .add(list(S_One, s(5))))),
         "2/5" = .mul(list(.rat(1, 4), .add(list(S_NegativeOne, s(5))))),
         "1/8" = .sqrt(.add(list(S_Half, .mul(list(.rat(1, 4), s(2)))))),
         "3/8" = .sqrt(.sub(S_Half, .mul(list(.rat(1, 4), s(2))))),
         "1/10" = .sqrt(.add(list(.rat(5, 8), .mul(list(.rat(1, 8), s(5)))))),
         "3/10" = .sqrt(.sub(.rat(5, 8), .mul(list(.rat(1, 8), s(5))))),
         "1/12" = .mul(list(.rat(1, 4), .add(list(s(2), s(6))))),
         "5/12" = .mul(list(.rat(1, 4), .sub(s(6), s(2)))),
         NULL)
}

.mod2 <- function(q) q - 2 * gmp::as.bigq(.bigq_floor(q / 2))

.cos_rational <- function(q) {
  r <- .mod2(q)
  sgn <- 1
  if (r > 1) r <- 2 - r
  if (r > gmp::as.bigq(1, 2)) {
    r <- 1 - r
    sgn <- -1
  }
  v <- .cos_table(r)
  if (is.null(v)) return(NULL)
  if (sgn < 0) .neg(v) else v
}

.tan_table <- function(r) {
  s <- function(n) .sqrt(.int(n))
  switch(as.character(r),
         "0" = S_Zero,
         "1/6" = .mul(list(.rat(1, 3), s(3))),
         "1/4" = S_One,
         "1/3" = s(3),
         "1/2" = S_ComplexInfinity,
         "1/12" = .sub(S_Two, s(3)),
         "5/12" = .add(list(S_Two, s(3))),
         "1/8" = .sub(s(2), S_One),
         "3/8" = .add(list(S_One, s(2))),
         "1/5" = .sqrt(.sub(.int(5), .mul(list(S_Two, s(5))))),
         "2/5" = .sqrt(.add(list(.int(5), .mul(list(S_Two, s(5)))))),
         "1/10" = .sqrt(.sub(S_One, .mul(list(.rat(2, 5), s(5))))),
         "3/10" = .sqrt(.add(list(S_One, .mul(list(.rat(2, 5), s(5)))))),
         NULL)
}

.tan_rational <- function(q) {
  r <- q - gmp::as.bigq(.bigq_floor(q))
  sgn <- 1
  if (r > gmp::as.bigq(1, 2)) {
    r <- 1 - r
    sgn <- -1
  }
  v <- .tan_table(r)
  if (is.null(v)) return(NULL)
  if (sgn < 0 && !is_const(v, "zoo")) .neg(v) else v
}

.trig_inverse_simplify <- function(name, a) {
  if (a$type != "Function") return(NULL)
  inner <- a$name
  if (identical(inner, "atan2")) {
    y <- a$args[[1L]]; x <- a$args[[2L]]
    r <- function() .sqrt(.add(list(.sq(x), .sq(y))))
    return(switch(name, sin = .div(y, r()), cos = .div(x, r()), tan = .div(y, x), cot = .div(x, y), NULL))
  }
  z <- a$args[[1L]]
  # the trig.eval tables: sqrt(1 - z^2), sqrt(1 + z^2), sqrt(1 - 1/z^2) and
  # sqrt(1 + 1/z^2)
  om <- function() .sqrt(.sub(S_One, .sq(z)))
  op <- function() .sqrt(.add(list(S_One, .sq(z))))
  omi <- function() .sqrt(.sub(S_One, .pow(z, Integer(-2))))
  opi <- function() .sqrt(.add(list(S_One, .pow(z, Integer(-2)))))
  switch(paste(name, inner),
         "sin asin" = z, "cos acos" = z, "tan atan" = z, "cot acot" = z,
         "sec asec" = z, "csc acsc" = z,
         "sin atan" = .div(z, op()), "sin acos" = om(), "sin acot" = .inv(.mul(list(opi(), z))),
         "sin acsc" = .inv(z), "sin asec" = omi(),
         "cos atan" = .inv(op()), "cos asin" = om(), "cos acot" = .inv(opi()),
         "cos acsc" = omi(), "cos asec" = .inv(z),
         "tan asin" = .div(z, om()), "tan acos" = .div(om(), z), "tan acot" = .inv(z),
         "tan acsc" = .inv(.mul(list(omi(), z))), "tan asec" = .mul(list(omi(), z)),
         "cot atan" = .inv(z), "cot asin" = .div(om(), z), "cot acos" = .div(z, om()),
         "cot acsc" = .mul(list(omi(), z)), "cot asec" = .inv(.mul(list(omi(), z))),
         "sinh asinh" = z, "cosh acosh" = z, "tanh atanh" = z,
         "coth acoth" = z,
         NULL)
}

.trig_eval <- function(name) {
  force(name)
  even <- name %in% c("cos", "sec")
  recip <- name %in% c("sec", "csc")
  base_eval <- if (recip) .trig_eval(if (name == "sec") "cos" else "sin")
  function(args) {
    a <- args[[1L]]
    if (a$type == "Float") return(.float_eval(name, args))
    if (is_nan(a) || is_const(a, "zoo")) return(S_NaN)
    if (!recip && (is_const(a, "oo") || is_const(a, "-oo")))
      return(if (name %in% c("sin", "cos")) AccumBounds(S_NegativeOne, S_One)
             else AccumBounds(S_NegativeInfinity, S_Infinity))
    # SymPy's order: the sign, then (sec, csc) the reflection of a rational
    # multiple of pi and the reciprocal of cos.eval or sin.eval, then the
    # rational multiples of pi, I*x, a symbolic multiple of pi, and a peeled
    # multiple of pi/2 in a sum.
    if (could_extract_minus_sign(a))
      return(if (even) fcall(name, list(.neg(a))) else .neg(fcall(name, list(.neg(a)))))
    inv <- .trig_inverse_simplify(name, a)
    if (!is.null(inv)) return(inv)
    one <- gmp::as.bigq(1); half <- gmp::as.bigq(1, 2)
    qp <- function(r) .mul(list(make_rational(r), S_Pi))
    if (recip) {
      pc <- .pi_coeff_sym(a)
      if (!is.null(pc) && is_exact_number(pc) && gmp::denominator(2 * pc$value) != 1) {
        c0 <- gmp::as.bigq(pc$value)
        if (c0 > one) return(.neg(fcall(name, list(qp(c0 - one)))))
        if (2 * c0 > one) {
          v <- fcall(name, list(qp(one - c0)))
          return(if (even) .neg(v) else v)
        }
      }
      t <- base_eval(list(a))
      if (is.null(t)) return(NULL)
      neg <- t$type == "Mul" && length(t$args) == 2L && identical(t$args[[1L]]$key, S_NegativeOne$key)
      u <- if (neg) t$args[[2L]] else t
      if (u$type == "Function" && u$name %in% c("cos", "sin")) {
        v <- fcall(if (u$name == "cos") "sec" else "csc", u$args)
        return(if (neg) .neg(v) else v)
      }
      return(.inv(t))
    }
    r <- .pi_coeff(a)
    if (is.null(r)) {
      # _pi_coeff recasts exact binary Float coefficients (0.5*pi)
      r2 <- .pi_coeff_sym(a)
      if (!is.null(r2) && is_exact_number(r2) && .has_float(a)) r <- r2
    }
    if (!is.null(r)) {
      q <- r$value
      v <- switch(name,
                  cos = .cos_rational(q),
                  sin = .cos_rational(half - q),
                  tan = .tan_rational(q),
                  cot = .tan_rational(half - q))
      if (!is.null(v)) return(v)
      # SymPy's reduction of a rational multiple of pi into [0, pi/2]; cos
      # does not rewrite cos(15*pi/7) as cos(pi/7)
      md <- function(v, m) v - m * .bigq_floor(v / m)
      if (name == "sin") {
        rr <- md(q, gmp::as.bigq(2))
        if (rr > one) return(.neg(fcall(name, list(qp(rr - one)))))
        if (rr > half) return(fcall(name, list(qp(one - rr))))
        if (rr != q) return(fcall(name, list(qp(rr))))
      } else if (name == "cos") {
        rr <- md(q, gmp::as.bigq(2))
        if (rr > one) return(.neg(fcall(name, list(qp(rr - one)))))
        if (rr > half) return(.neg(fcall(name, list(qp(one - rr)))))
      } else {
        rr <- md(q, one)
        if (rr > half) return(.neg(fcall(name, list(qp(one - rr)))))
        if (rr != q) return(fcall(name, list(qp(rr))))
      }
      return(NULL)
    }
    ic <- .imag_coeff(a)
    if (!is.null(ic)) {
      return(switch(name,
                    sin = .mul(list(S_I, fcall("sinh", list(ic)))),
                    cos = fcall("cosh", list(ic)),
                    tan = .mul(list(S_I, fcall("tanh", list(ic)))),
                    cot = .mul(list(S_NegativeOne, S_I, fcall("coth", list(ic))))))
    }
    pc <- .pi_coeff_sym(a)
    if (!is.null(pc)) {
      if (isTRUE(ask_fact(pc, "integer")))
        return(switch(name, sin = S_Zero, tan = S_Zero, cos = .pow(S_NegativeOne, pc), cot = S_ComplexInfinity))
      if (name %in% c("sin", "cos") && isTRUE(ask_fact(.mul(list(S_Two, pc)), "integer")) &&
          identical(ask_fact(pc, "even"), FALSE))
        return(if (name == "sin") .pow(S_NegativeOne, .sub(pc, S_Half)) else S_Zero)
      if (!is_exact_number(pc)) {
        narg <- .mul(list(pc, S_Pi))
        return(if (identical(narg$key, a$key)) NULL else fcall(name, list(narg)))
      }
    }
    if (a$type == "Add") {
      pp <- .peeloff_pi(a)
      if (!is.null(pp)) {
        x <- pp[[1L]]
        mpi <- .mul(list(pp[[2L]], S_Pi))
        return(switch(name,
                      sin = .add(list(.mul(list(fcall("sin", list(mpi)), fcall("cos", list(x)))),
                                      .mul(list(fcall("cos", list(mpi)), fcall("sin", list(x)))))),
                      cos = .sub(.mul(list(fcall("cos", list(mpi)), fcall("cos", list(x)))),
                                 .mul(list(fcall("sin", list(mpi)), fcall("sin", list(x))))),
                      tan = if (is_const(fcall("tan", list(mpi)), "zoo")) .neg(fcall("cot", list(x)))
                            else fcall("tan", list(x)),
                      cot = if (is_const(fcall("cot", list(mpi)), "zoo")) fcall("cot", list(x))
                            else .neg(fcall("tan", list(x)))))
      }
    }
    NULL
  }
}

.trig_facts <- function(name) {
  force(name)
  function(args, afs) {
    f <- afs[[1L]]
    out <- c(commutative = TRUE)
    if (isTRUE(f[["real"]])) {
      out["real"] <- TRUE
      if (name %in% c("sin", "cos")) out["finite"] <- TRUE
    }
    if (isTRUE(f[["rational"]]) && isTRUE(f[["nonzero"]]) && name %in% c("sin", "cos", "tan"))
      out["transcendental"] <- TRUE
    out
  }
}

def_function("sin", eval = .trig_eval("sin"),
             diff = function(args, i) fcall("cos", args),
             numeric = sin, facts = .trig_facts("sin"))
def_function("cos", eval = .trig_eval("cos"),
             diff = function(args, i) .neg(fcall("sin", args)),
             numeric = cos, facts = .trig_facts("cos"))
def_function("tan", eval = .trig_eval("tan"),
             diff = function(args, i) .add(list(S_One, .sq(fcall("tan", args)))),
             numeric = tan, facts = .trig_facts("tan"))
def_function("cot", eval = .trig_eval("cot"),
             diff = function(args, i) .sub(S_NegativeOne, .sq(fcall("cot", args))),
             numeric = function(x) 1 / tan(x), facts = .trig_facts("cot"))
def_function("sec", eval = .trig_eval("sec"),
             diff = function(args, i) .mul(list(fcall("sec", args), fcall("tan", args))),
             numeric = function(x) 1 / cos(x), facts = .trig_facts("sec"))
def_function("csc", eval = .trig_eval("csc"),
             diff = function(args, i) .neg(.mul(list(fcall("cot", args), fcall("csc", args)))),
             numeric = function(x) 1 / sin(x), facts = .trig_facts("csc"))

def_function("sinc",
             eval = function(args) {
               a <- args[[1L]]
               if (is_zero_num(a)) return(S_One)
               if (a$type == "Float") return(.float_eval("sinc", args))
               r <- .pi_coeff(a)
               if (!is.null(r) && gmp::denominator(r$value) == 1) return(S_Zero)
               # half-integer multiples of pi: (-1)^(k - 1/2)/a
               if (!is.null(r) && gmp::denominator(r$value) == 2)
                 return(.div(.pow(S_NegativeOne, .sub(r, S_Half)), a))
               .even_eval("sinc", a)
             },
             diff = function(args, i) {
               x <- args[[1L]]
               .sub(.div(fcall("cos", list(x)), x), .div(fcall("sin", list(x)), .sq(x)))
             },
             numeric = function(x) if (x == 0) 1 else sin(x) / x,
             facts = .trig_facts("sin"))

# ---------------------------------------------------------------------------
# Inverse trigonometric functions

.asin_table <- function(a) {
  s <- function(n) .sqrt(.int(n))
  vals <- list(
    list(S_Zero, S_Zero), list(S_One, .mul(list(S_Half, S_Pi))),
    list(S_Half, .mul(list(.rat(1, 6), S_Pi))),
    list(.mul(list(S_Half, s(2))), .mul(list(.rat(1, 4), S_Pi))),
    list(.mul(list(S_Half, s(3))), .mul(list(.rat(1, 3), S_Pi))),
    list(.mul(list(.rat(1, 4), .sub(s(6), s(2)))), .mul(list(.rat(1, 12), S_Pi))),
    list(.mul(list(.rat(1, 4), .add(list(s(6), s(2))))), .mul(list(.rat(5, 12), S_Pi))),
    list(.mul(list(.rat(1, 4), .sub(s(5), S_One))), .mul(list(.rat(1, 10), S_Pi))),
    list(.mul(list(.rat(1, 4), .add(list(s(5), S_One)))), .mul(list(.rat(3, 10), S_Pi)))
  )
  for (v in vals) if (identical(v[[1L]]$key, a$key)) return(v[[2L]])
  NULL
}

.atan_table <- function(a) {
  s <- function(n) .sqrt(.int(n))
  vals <- list(
    list(S_Zero, S_Zero), list(S_One, .mul(list(.rat(1, 4), S_Pi))),
    list(s(3), .mul(list(.rat(1, 3), S_Pi))),
    list(.mul(list(.rat(1, 3), s(3))), .mul(list(.rat(1, 6), S_Pi))),
    list(.sub(S_Two, s(3)), .mul(list(.rat(1, 12), S_Pi))),
    list(.add(list(S_Two, s(3))), .mul(list(.rat(5, 12), S_Pi))),
    list(.sub(s(2), S_One), .mul(list(.rat(1, 8), S_Pi))),
    list(.add(list(S_One, s(2))), .mul(list(.rat(3, 8), S_Pi))),
    list(S_Infinity, .mul(list(S_Half, S_Pi)))
  )
  for (v in vals) if (identical(v[[1L]]$key, a$key)) return(v[[2L]])
  NULL
}

def_function("asin",
             eval = function(args) {
               a <- args[[1L]]
               if (a$type == "Float") return(.float_eval("asin", args))
               v <- .asin_table(a)
               if (!is.null(v)) return(v)
               if (could_extract_minus_sign(a)) {
                 v <- .asin_table(.neg(a))
                 if (!is.null(v)) return(.neg(v))
               }
               ic <- .imag_coeff(a)
               if (!is.null(ic)) return(.mul(list(S_I, fcall("asinh", list(ic)))))
               .odd_eval("asin", a)
             },
             diff = function(args, i) .inv(.sqrt(.sub(S_One, .sq(args[[1L]])))),
             numeric = function(x) if (is.complex(x) || abs(x) > 1) asin(as.complex(x)) else asin(x))

def_function("acos",
             eval = function(args) {
               a <- args[[1L]]
               if (a$type == "Float") return(.float_eval("acos", args))
               # SymPy's acos.eval: no reflection of negative arguments
               if (is_const(a, "oo")) return(.mul(list(S_Infinity, S_I)))
               if (is_const(a, "-oo")) return(.mul(list(S_NegativeInfinity, S_I)))
               if (is_const(a, "zoo")) return(S_ComplexInfinity)
               v <- .asin_table(a)
               if (!is.null(v)) return(.sub(.mul(list(S_Half, S_Pi)), v))
               if (could_extract_minus_sign(a)) {
                 v <- .asin_table(.neg(a))
                 if (!is.null(v)) return(.add(list(.mul(list(S_Half, S_Pi)), v)))
               }
               if (!is.null(.imag_coeff(a))) return(.sub(.mul(list(S_Half, S_Pi)), fcall("asin", list(a))))
               NULL
             },
             diff = function(args, i) .neg(.inv(.sqrt(.sub(S_One, .sq(args[[1L]]))))),
             numeric = function(x) if (is.complex(x) || abs(x) > 1) acos(as.complex(x)) else acos(x))

def_function("atan",
             eval = function(args) {
               a <- args[[1L]]
               if (a$type == "Float") return(.float_eval("atan", args))
               v <- .atan_table(a)
               if (!is.null(v)) return(v)
               if (could_extract_minus_sign(a)) {
                 v <- .atan_table(.neg(a))
                 if (!is.null(v)) return(.neg(v))
                 return(.neg(fcall("atan", list(.neg(a)))))
               }
               ic <- .imag_coeff(a)
               if (!is.null(ic)) return(.mul(list(S_I, fcall("atanh", list(ic)))))
               NULL
             },
             diff = function(args, i) .inv(.add(list(S_One, .sq(args[[1L]])))),
             numeric = atan,
             facts = function(args, afs) {
               out <- c(commutative = TRUE)
               if (isTRUE(afs[[1L]][["real"]])) {
                 out["real"] <- TRUE
                 out["finite"] <- TRUE
               }
               if (isTRUE(afs[[1L]][["positive"]])) out["positive"] <- TRUE
               if (isTRUE(afs[[1L]][["negative"]])) out["negative"] <- TRUE
               out
             })

def_function("acot",
             eval = function(args) {
               a <- args[[1L]]
               if (a$type == "Float") return(.float_eval("acot", args))
               if (is_zero_num(a)) return(.mul(list(S_Half, S_Pi)))
               if (is_const(a, "oo")) return(S_Zero)
               if (is_const(a, "-oo") || is_const(a, "zoo")) return(S_Zero)
               v <- .atan_table(a)
               if (!is.null(v)) return(.sub(.mul(list(S_Half, S_Pi)), v))
               # acot(-x) = -acot(x) only away from zero, where acot(0) = pi/2
               if (isTRUE(ask_fact(a, "nonzero"))) {
                 r <- .odd_eval("acot", a)
                 if (!is.null(r)) return(r)
               }
               ic <- .imag_coeff(a)
               if (!is.null(ic)) return(.mul(list(.neg(S_I), fcall("acoth", list(ic)))))
               NULL
             },
             diff = function(args, i) .neg(.inv(.add(list(S_One, .sq(args[[1L]]))))),
             numeric = function(x) atan(1 / x))

def_function("asec",
             eval = function(args) {
               a <- args[[1L]]
               if (a$type == "Float") return(.float_eval("asec", args))
               if (is_one_num(a)) return(S_Zero)
               if (is_negone_num(a)) return(S_Pi)
               if (is_zero_num(a)) return(S_ComplexInfinity)
               if (is_inf_const(a)) return(.mul(list(S_Half, S_Pi)))
               # SymPy's acsc table, read through asin(1/a)
               if (!length(free_symbols(a))) {
                 v <- .asin_table(.inv(a))
                 if (!is.null(v)) return(.sub(.mul(list(S_Half, S_Pi)), v))
                 v <- .asin_table(.inv(.neg(a)))
                 if (!is.null(v)) return(.add(list(.mul(list(S_Half, S_Pi)), v)))
               }
               NULL
             },
             diff = function(args, i) {
               x <- args[[1L]]
               .inv(.mul(list(.sq(x), .sqrt(.sub(S_One, .inv(.sq(x)))))))
             },
             numeric = function(x) acos(as.complex(1 / x)))

def_function("acsc",
             eval = function(args) {
               a <- args[[1L]]
               if (a$type == "Float") return(.float_eval("acsc", args))
               if (is_one_num(a)) return(.mul(list(S_Half, S_Pi)))
               if (is_negone_num(a)) return(.mul(list(.rat(-1, 2), S_Pi)))
               if (is_zero_num(a)) return(S_ComplexInfinity)
               if (is_inf_const(a)) return(S_Zero)
               r <- .odd_eval("acsc", a)
               if (!is.null(r)) return(r)
               if (!length(free_symbols(a))) return(.asin_table(.inv(a)))
               NULL
             },
             diff = function(args, i) {
               x <- args[[1L]]
               .neg(.inv(.mul(list(.sq(x), .sqrt(.sub(S_One, .inv(.sq(x))))))))
             },
             numeric = function(x) asin(as.complex(1 / x)))

def_function("atan2", nargs = 2L,
             eval = function(args) {
               y <- args[[1L]]
               x <- args[[2L]]
               if (is_zero_num(y) && is_zero_num(x)) return(S_NaN)
               if (isTRUE(ask_fact(x, "positive")) && isTRUE(ask_fact(y, "real")))
                 return(fcall("atan", list(.div(y, x))))
               if (is_zero_num(x)) {
                 if (isTRUE(ask_fact(y, "positive"))) return(.mul(list(S_Half, S_Pi)))
                 if (isTRUE(ask_fact(y, "negative"))) return(.mul(list(.rat(-1, 2), S_Pi)))
               }
               if (isTRUE(ask_fact(x, "negative"))) {
                 if (isTRUE(ask_fact(y, "nonnegative")))
                   return(.add(list(fcall("atan", list(.div(y, x))), S_Pi)))
                 if (isTRUE(ask_fact(y, "negative")))
                   return(.sub(fcall("atan", list(.div(y, x))), S_Pi))
               }
               NULL
             },
             diff = function(args, i) {
               y <- args[[1L]]
               x <- args[[2L]]
               den <- .add(list(.sq(x), .sq(y)))
               if (i == 1L) .div(x, den) else .neg(.div(y, den))
             },
             numeric = function(y, x) atan2(Re(y), Re(x)))

# ---------------------------------------------------------------------------
# Hyperbolic functions

.hyp_eval <- function(name) {
  force(name)
  function(args) {
    a <- args[[1L]]
    if (a$type == "Float") return(.float_eval(name, args))
    if (is_nan(a)) return(S_NaN)
    inv <- .trig_inverse_simplify(name, a)
    if (!is.null(inv)) return(inv)
    if (is_zero_num(a)) {
      return(switch(name, sinh = S_Zero, cosh = S_One, tanh = S_Zero,
                    coth = S_ComplexInfinity, sech = S_One, csch = S_ComplexInfinity))
    }
    if (is_const(a, "oo")) {
      return(switch(name, sinh = S_Infinity, cosh = S_Infinity, tanh = S_One,
                    coth = S_One, sech = S_Zero, csch = S_Zero))
    }
    if (is_const(a, "-oo")) {
      return(switch(name, sinh = S_NegativeInfinity, cosh = S_Infinity,
                    tanh = S_NegativeOne, coth = S_NegativeOne, sech = S_Zero,
                    csch = S_Zero))
    }
    ic <- .imag_coeff(a)
    if (!is.null(ic)) {
      return(switch(name,
                    sinh = .mul(list(S_I, fcall("sin", list(ic)))),
                    cosh = fcall("cos", list(ic)),
                    tanh = .mul(list(S_I, fcall("tan", list(ic)))),
                    coth = .mul(list(S_NegativeOne, S_I, fcall("cot", list(ic)))),
                    sech = fcall("sec", list(ic)),
                    csch = .mul(list(S_NegativeOne, S_I, fcall("csc", list(ic))))))
    }
    if (name %in% c("cosh", "sech")) return(.even_eval(name, a))
    .odd_eval(name, a)
  }
}

.hyp_facts <- function(name) {
  force(name)
  function(args, afs) {
    f <- afs[[1L]]
    out <- c(commutative = TRUE)
    if (isTRUE(f[["real"]])) {
      out["real"] <- TRUE
      if (name %in% c("cosh", "sech")) out["positive"] <- TRUE
      if (name %in% c("sinh", "tanh", "cosh", "sech")) out["finite"] <- TRUE
    }
    if (isTRUE(f[["positive"]]) && name %in% c("sinh", "tanh", "coth", "csch"))
      out["positive"] <- TRUE
    if (isTRUE(f[["negative"]]) && name %in% c("sinh", "tanh", "coth", "csch"))
      out["negative"] <- TRUE
    out
  }
}

def_function("sinh", eval = .hyp_eval("sinh"),
             diff = function(args, i) fcall("cosh", args), numeric = sinh,
             facts = .hyp_facts("sinh"))
def_function("cosh", eval = .hyp_eval("cosh"),
             diff = function(args, i) fcall("sinh", args), numeric = cosh,
             facts = .hyp_facts("cosh"))
def_function("tanh", eval = .hyp_eval("tanh"),
             diff = function(args, i) .sub(S_One, .sq(fcall("tanh", args))),
             numeric = tanh, facts = .hyp_facts("tanh"))
def_function("coth", eval = .hyp_eval("coth"),
             diff = function(args, i) .neg(.inv(.sq(fcall("sinh", args)))),
             numeric = function(x) 1 / tanh(x), facts = .hyp_facts("coth"))
def_function("sech", eval = .hyp_eval("sech"),
             diff = function(args, i) .neg(.mul(list(fcall("tanh", args), fcall("sech", args)))),
             numeric = function(x) 1 / cosh(x), facts = .hyp_facts("sech"))
def_function("csch", eval = .hyp_eval("csch"),
             diff = function(args, i) .neg(.mul(list(fcall("coth", args), fcall("csch", args)))),
             numeric = function(x) 1 / sinh(x), facts = .hyp_facts("csch"))

.ahyp <- function(name, zero, diff, numeric, odd = TRUE, extra = NULL) {
  def_function(name,
               eval = function(args) {
                 a <- args[[1L]]
                 if (a$type == "Float") return(.float_eval(name, args))
                 if (is_zero_num(a)) return(zero)
                 if (!is.null(extra)) {
                   r <- extra(a)
                   if (!is.null(r)) return(r)
                 }
                 if (odd) .odd_eval(name, a) else NULL
               },
               diff = diff, numeric = numeric)
}

.ahyp("asinh", S_Zero,
      function(args, i) .inv(.sqrt(.add(list(.sq(args[[1L]]), S_One)))),
      asinh,
      extra = function(a) {
        if (is_one_num(a)) return(fcall("log", list(.add(list(S_One, .sqrt(S_Two))))))
        if (is_const(a, "oo")) return(S_Infinity)
        ic <- .imag_coeff(a)
        if (!is.null(ic)) return(.mul(list(S_I, fcall("asin", list(ic)))))
        NULL
      })
.ahyp("acosh", .mul(list(S_Half, S_I, S_Pi)),
      function(args, i) .inv(.mul(list(.sqrt(.sub(args[[1L]], S_One)), .sqrt(.add(list(args[[1L]], S_One)))))),
      function(x) if (is.complex(x) || x < 1) acosh(as.complex(x)) else acosh(x), odd = FALSE,
      extra = function(a) {
        if (is_one_num(a)) return(S_Zero)
        if (is_negone_num(a)) return(.mul(list(S_I, S_Pi)))
        if (is_const(a, "oo")) return(S_Infinity)
        NULL
      })
.ahyp("atanh", S_Zero,
      function(args, i) .inv(.sub(S_One, .sq(args[[1L]]))),
      function(x) if (is.complex(x) || abs(x) > 1) atanh(as.complex(x)) else atanh(x),
      extra = function(a) {
        if (is_one_num(a)) return(S_Infinity)
        if (is_negone_num(a)) return(S_NegativeInfinity)
        ic <- .imag_coeff(a)
        if (!is.null(ic)) return(.mul(list(S_I, fcall("atan", list(ic)))))
        NULL
      })
.ahyp("acoth", .mul(list(S_Half, S_I, S_Pi)),
      function(args, i) .inv(.sub(S_One, .sq(args[[1L]]))),
      function(x) atanh(as.complex(1 / x)),
      extra = function(a) {
        if (is_one_num(a)) return(S_Infinity)
        if (is_const(a, "oo")) return(S_Zero)
        NULL
      })
.ahyp("asech", S_Infinity,
      function(args, i) {
        # asech.fdiff: -sqrt(1/(z + 1))/(z*sqrt(1 - z))
        z <- args[[1L]]
        .neg(.div(.sqrt(.inv(.add(list(z, S_One)))), .mul(list(z, .sqrt(.sub(S_One, z))))))
      },
      function(x) acosh(as.complex(1 / x)), odd = FALSE,
      extra = function(a) if (is_one_num(a)) S_Zero else NULL)
.ahyp("acsch", S_ComplexInfinity,
      function(args, i) {
        x <- args[[1L]]
        .neg(.inv(.mul(list(.sq(x), .sqrt(.add(list(S_One, .inv(.sq(x)))))))))
      },
      function(x) asinh(1 / x),
      extra = function(a) {
        if (is_one_num(a)) return(fcall("log", list(.add(list(S_One, .sqrt(S_Two))))))
        NULL
      })

# ---------------------------------------------------------------------------
# Complex functions: Abs, sign, re, im, arg, conjugate

def_function("Abs",
             eval = function(args) .abs_eval(args[[1L]]),
             diff = function(args, i) {
               x <- args[[1L]]
               if (isTRUE(ask_fact(x, "real"))) return(fcall("sign", list(x)))
               NULL
             },
             numeric = function(x) abs(x),
             # Abs: extended real and nonnegative; real only for a finite
             # argument, positive when the argument is nonzero
             facts = function(args, afs) {
               f <- afs[[1L]]
               out <- c(commutative = TRUE, extended_real = TRUE, extended_negative = FALSE,
                        extended_nonnegative = TRUE)
               if (isTRUE(f[["finite"]])) out["real"] <- TRUE
               if (isTRUE(f[["imaginary"]])) out["zero"] <- FALSE
               z <- f[["zero"]]
               if (!is.na(z)) { out["zero"] <- z; out["extended_positive"] <- !z }
               else if (isTRUE(f[["imaginary"]])) out["extended_positive"] <- TRUE
               if (isTRUE(f[["extended_real"]]))
                 for (nm in c("integer", "rational", "even", "odd")) if (!is.na(f[[nm]])) out[nm] <- f[[nm]]
               if (!is.na(f[["algebraic"]])) out["algebraic"] <- f[["algebraic"]]
               out
             },
             latex = function(a) paste0("\\left|{", a[1], "}\\right|"))

def_function("sign",
             # sign.eval
             eval = function(args) {
               a <- args[[1L]]
               if (a$type == "Mul") {
                 c <- if (is_number_type(a$args[[1L]])) a$args[[1L]] else S_One
                 fs <- if (is_number_type(a$args[[1L]])) a$args[-1L] else a$args
                 s <- fcall("sign", list(c))
                 unk <- list()
                 for (f in fs) {
                   ff <- .facts(f)
                   if (isTRUE(ff[["extended_negative"]])) s <- .neg(s)
                   else if (isTRUE(ff[["extended_positive"]])) next
                   else if (isTRUE(ff[["imaginary"]])) {
                     ai <- fcall("im", list(f))
                     if (.expr_is_number(ai) && isTRUE(ask_fact(ai, "real"))) {
                       s <- .mul(list(s, S_I))
                       if (isTRUE(ask_fact(ai, "extended_negative"))) s <- .neg(s)
                     } else unk[[length(unk) + 1L]] <- f
                   } else unk[[length(unk) + 1L]] <- f
                 }
                 if (is_one_num(c) && length(unk) == length(fs)) return(NULL)
                 return(.mul(list(s, fcall("sign", list(if (length(unk)) new_sym("Mul", unk) else S_One)))))
               }
               if (is_nan(a)) return(S_NaN)
               fa <- .facts(a)
               if (isTRUE(fa[["zero"]])) return(S_Zero)
               if (isTRUE(fa[["extended_positive"]])) return(S_One)
               if (isTRUE(fa[["extended_negative"]])) return(S_NegativeOne)
               if (a$type == "Function" && a$name == "sign") return(a)
               if (isTRUE(fa[["imaginary"]])) {
                 if (a$type == "Pow" && identical(a$args[[2L]]$key, S_Half$key)) return(S_I)
                 a2 <- .mul(list(.neg(S_I), a))
                 if (isTRUE(ask_fact(a2, "extended_positive"))) return(S_I)
                 if (isTRUE(ask_fact(a2, "extended_negative"))) return(.neg(S_I))
               }
               NULL
             },
             diff = function(args, i) {
               x <- args[[1L]]
               if (isTRUE(ask_fact(x, "real"))) return(.mul(list(S_Two, fcall("DiracDelta", list(x)))))
               S_Zero
             },
             numeric = function(x) if (is.complex(x)) {
               if (x == 0) 0 else x / Mod(x)
             } else sign(x),
             facts = function(args, afs) {
               f <- afs[[1L]]
               out <- c(commutative = TRUE, finite = TRUE)
               if (isTRUE(f[["real"]])) out["integer"] <- TRUE
               if (isTRUE(f[["nonzero"]])) out["zero"] <- FALSE
               if (isTRUE(f[["zero"]])) out["zero"] <- TRUE
               out
             })

# Abs.eval
.abs_eval <- function(a) {
  t <- a$type
  if (t == "Integer" || t == "Rational") return(make_rational(abs(a$value)))
  if (t == "Float") return(Float(abs(a$value)))
  # the argument's _eval_Abs
  if (t == "Function") {
    if (a$name == "sign" && identical(ask_fact(a$args[[1L]], "zero"), FALSE)) return(S_One)
    if (a$name %in% c("conjugate", "polar_lift")) return(fcall("Abs", a$args))
    if (a$name == "exp_polar") return(fcall("exp", list(fcall("re", a$args))))
  }
  # signsimp(arg, evaluate = False)
  if (t == "Add" && could_extract_minus_sign(a)) return(fcall("Abs", list(.neg(a))))
  nd <- .as_numer_denom(a)
  if (length(free_symbols(nd[[2L]])) && !length(free_symbols(nd[[1L]])))
    return(.div(fcall("Abs", list(nd[[1L]])), fcall("Abs", list(nd[[2L]]))))
  is_abs <- function(z) z$type == "Function" && z$name == "Abs"
  if (t == "Mul") {
    known <- list(); unk <- list()
    for (f in a$args) {
      if (f$type == "Pow" && isTRUE(ask_fact(f$args[[2L]], "integer")) && isTRUE(ask_fact(f$args[[2L]], "negative"))) {
        bnew <- fcall("Abs", list(f$args[[1L]]))
        if (is_abs(bnew)) unk[[length(unk) + 1L]] <- f else known[[length(known) + 1L]] <- .pow(bnew, f$args[[2L]])
      } else {
        tnew <- fcall("Abs", list(f))
        if (is_abs(tnew)) unk[[length(unk) + 1L]] <- f else known[[length(known) + 1L]] <- tnew
      }
    }
    u <- if (length(unk)) fnode("Abs", list(.mul(unk))) else S_One
    return(.mul(c(known, list(u))))
  }
  if (is_nan(a)) return(S_NaN)
  if (is_const(a, "zoo")) return(S_Infinity)
  if (t == "Pow") {
    base <- a$args[[1L]]; ex <- a$args[[2L]]
    if (isTRUE(ask_fact(base, "extended_real"))) {
      if (isTRUE(ask_fact(ex, "integer"))) {
        if (isTRUE(ask_fact(ex, "even"))) return(a)
        if (is_negone_num(base)) return(S_One)
        return(.pow(fcall("Abs", list(base)), ex))
      }
      if (isTRUE(ask_fact(base, "extended_nonnegative"))) return(.pow(base, fcall("re", list(ex))))
      if (isTRUE(ask_fact(base, "extended_negative")))
        return(.mul(list(.pow(.neg(base), fcall("re", list(ex))), fcall("exp", list(.neg(.mul(list(S_Pi, fcall("im", list(ex))))))))))
      return(NULL)
    }
    if (!.has_type(base, c("Symbol", "Dummy"))) {
      ri <- .re_im(fcall("log", list(base)))
      if (!is.null(ri)) {
        lz <- .add(list(ri[[1L]], .mul(list(S_I, ri[[2L]]))))
        return(fcall("exp", list(fcall("re", list(.mul(list(ex, lz)))))))
      }
    }
  }
  if (t == "Function" && a$name == "exp") return(fcall("exp", list(fcall("re", a$args))))
  if (t == "Function" && is.null(.FUNS[[a$name]])) {
    if (isTRUE(ask_fact(a, "positive"))) return(a)
    if (isTRUE(ask_fact(a, "negative"))) return(.neg(a))
    return(NULL)
  }
  if (t == "Add" && (.has1(a, S_Infinity) || .has1(a, S_NegativeInfinity))) {
    ri <- .re_im(a)
    if (!is.null(ri) && any(vapply(ri, function(p) isTRUE(ask_fact(p, "infinite")), TRUE))) return(S_Infinity)
  }
  if (isTRUE(ask_fact(a, "zero"))) return(S_Zero)
  if (isTRUE(ask_fact(a, "extended_nonnegative"))) return(a)
  if (isTRUE(ask_fact(a, "extended_nonpositive"))) return(.neg(a))
  if (isTRUE(ask_fact(a, "imaginary"))) {
    a2 <- .mul(list(.neg(S_I), a))
    if (isTRUE(ask_fact(a2, "extended_nonnegative"))) return(a2)
  }
  if (isTRUE(ask_fact(a, "extended_real"))) return(NULL)
  conj <- fcall("conjugate", list(a))
  if (conj$type == "Add" && could_extract_minus_sign(conj)) conj <- .neg(.neg(conj))
  ca <- function(z) { out <- list(); walk <- function(w) { if (w$type == "Function" && w$name == "conjugate") out[[w$key]] <<- w; for (b in w$args) walk(b) }; walk(z); out }
  newc <- ca(conj); oldc <- ca(a)
  newc <- newc[setdiff(names(newc), names(oldc))]
  if (length(newc) && all(vapply(newc, function(w) .has1(a, w$args[[1L]]), TRUE))) return(NULL)
  if (!identical(a$key, conj$key) && !identical(a$key, .neg(conj)$key)) {
    absa <- atoms(a, "Function"); absa <- Filter(function(w) w$name == "Abs", absa)
    afree <- if (length(absa)) xreplace(a, lapply(absa, function(w) list(w, Dummy("r", real = TRUE)))) else a
    unk <- Filter(function(s) is.na(ask_fact(s, "extended_real")), free_symbols(afree))
    if (!length(unk) || !all(vapply(unk, function(u) .has1(conj, fcall("conjugate", list(u))), TRUE)))
      return(.sqrt(expand_mul(.mul(list(a, conj)))))
  }
  NULL
}

# Split e into real and imaginary parts when every piece is decidable.
.re_im <- function(e) {
  t <- e$type
  fe <- .facts(e)
  if (isTRUE(fe[["real"]])) return(list(e, S_Zero))
  if (is_const(e, "I")) return(list(S_Zero, S_One))
  if (isTRUE(fe[["imaginary"]])) {
    ic <- .imag_coeff(e)
    if (!is.null(ic) && isTRUE(ask_fact(ic, "real"))) return(list(S_Zero, ic))
  }
  if (t == "Add") {
    re <- list()
    im <- list()
    for (a in e$args) {
      r <- .re_im(a)
      re[[length(re) + 1L]] <- r[[1L]]
      im[[length(im) + 1L]] <- r[[2L]]
    }
    return(list(.add(re), .add(im)))
  }
  if (t == "Mul") {
    ri <- list(S_One, S_Zero)
    for (a in e$args) {
      r <- .re_im(a)
      ri <- list(.sub(.mul(list(ri[[1L]], r[[1L]])), .mul(list(ri[[2L]], r[[2L]]))),
                 .add(list(.mul(list(ri[[1L]], r[[2L]])), .mul(list(ri[[2L]], r[[1L]])))))
    }
    return(ri)
  }
  if (t == "Pow" && e$args[[2L]]$type == "Integer" && e$args[[2L]]$value > 0 &&
      e$args[[2L]]$value < 50) {
    r <- .re_im(e$args[[1L]])
    out <- list(S_One, S_Zero)
    for (k in seq_len(as.integer(as.character(e$args[[2L]]$value)))) {
      out <- list(.sub(.mul(list(out[[1L]], r[[1L]])), .mul(list(out[[2L]], r[[2L]]))),
                  .add(list(.mul(list(out[[1L]], r[[2L]])), .mul(list(out[[2L]], r[[1L]])))))
    }
    return(out)
  }
  list(fnode("re", list(e)), fnode("im", list(e)))
}

.unevaluated_re_im <- function(ri) {
  has_re_im <- function(x) {
    if (x$type == "Function" && x$name %in% c("re", "im")) return(TRUE)
    if (!length(x$args)) return(FALSE)
    any(vapply(x$args, has_re_im, TRUE))
  }
  if (has_re_im(ri[[1L]]) || has_re_im(ri[[2L]])) TRUE else NULL
}

.re_im_eval <- function(which) {
  force(which)
  function(args) {
    a <- args[[1L]]
    if (is_nan(a)) return(S_NaN)
    if (a$type == "Function" && a$name %in% c("re", "im")) {
      return(if (which == 1L) a else S_Zero)
    }
    ri <- .re_im(a)
    out <- ri[[which]]
    if (out$type == "Function" && out$name %in% c("re", "im") &&
        identical(out$args[[1L]]$key, a$key)) return(NULL)
    out
  }
}

# re and im: extended real; finite (so complex) for a finite argument; re is
# zero when the argument is zero or imaginary, im when it is extended real
.re_im_facts <- function(f, which) {
  out <- c(commutative = TRUE, extended_real = TRUE)
  if (isTRUE(f[["finite"]])) { out["finite"] <- TRUE; out["complex"] <- TRUE }
  z <- if (which == "re") {
    if (isTRUE(f[["imaginary"]]) || isTRUE(f[["zero"]])) TRUE
    else if (isFALSE(f[["imaginary"]]) && isFALSE(f[["zero"]])) FALSE else NA
  } else f[["extended_real"]]
  if (!is.na(z)) out["zero"] <- z
  if (!is.na(f[["algebraic"]])) out["algebraic"] <- f[["algebraic"]]
  out
}

def_function("re", eval = function(args) .re_eval_s(args), numeric = Re,
             diff = function(args, i) {
               if (isTRUE(ask_fact(args[[1L]], "real"))) return(S_One)
               NULL
             },
             facts = function(args, afs) .re_im_facts(afs[[1L]], "re"),
             latex = function(a) paste0("\\operatorname{re}{\\left(", a[1], " \\right)}"))
def_function("im", eval = function(args) .im_eval_s(args), numeric = Im,
             diff = function(args, i) {
               if (isTRUE(ask_fact(args[[1L]], "real"))) return(S_Zero)
               NULL
             },
             facts = function(args, afs) .re_im_facts(afs[[1L]], "im"),
             latex = function(a) paste0("\\operatorname{im}{\\left(", a[1], " \\right)}"))

def_function("conjugate",
             eval = function(args) {
               a <- args[[1L]]
               t <- a$type
               if (isTRUE(ask_fact(a, "real"))) return(a)
               if (is_const(a, "I")) return(.neg(S_I))
               if (t == "Function" && a$name == "conjugate") return(a$args[[1L]])
               if (!is.null(.qnode_def(t))) return(.eval_conjugate(a))
               if (isTRUE(ask_fact(a, "imaginary"))) return(.neg(a))
               if (t == "Add") return(.add(lapply(a$args, function(z) fcall("conjugate", list(z)))))
               if (t == "Mul" && all(vapply(a$args, .is_comm, TRUE)))
                 return(.mul(lapply(a$args, function(z) fcall("conjugate", list(z)))))
               if (t == "Pow" && isTRUE(ask_fact(a$args[[2L]], "integer")))
                 return(.pow(fcall("conjugate", list(a$args[[1L]])), a$args[[2L]]))
               if (t == "Function" && a$name %in% c("exp", "sin", "cos", "tan", "sinh", "cosh", "tanh"))
                 return(fcall(a$name, list(fcall("conjugate", a$args))))
               # a function's own _eval_conjugate
               if (t == "Function" && !is.null(.FUNS[[a$name]]$conjugate)) {
                 r <- .FUNS[[a$name]]$conjugate(a$args)
                 if (!is.null(r)) return(r)
               }
               NULL
             },
             numeric = Conj,
             facts = function(args, afs) {
               f <- afs[[1L]]
               out <- c(commutative = !identical(f[["commutative"]], FALSE))
               for (nm in c("real", "imaginary", "finite", "zero", "positive", "negative"))
                 if (!is.na(f[[nm]])) out[nm] <- f[[nm]]
               out
             },
             latex = function(a) paste0("\\overline{", a[1], "}"))

def_function("arg",
             eval = function(args) {
               a <- args[[1L]]
               fa <- .facts(a)
               if (isTRUE(fa[["positive"]])) return(S_Zero)
               if (isTRUE(fa[["negative"]])) return(S_Pi)
               if (is_const(a, "I")) return(.mul(list(S_Half, S_Pi)))
               if (length(free_symbols(a)) == 0L) {
                 ri <- .re_im(a)
                 if (is.null(.unevaluated_re_im(ri))) return(fcall("atan2", list(ri[[2L]], ri[[1L]])))
               }
               NULL
             },
             numeric = Arg,
             facts = function(args, afs) c(commutative = TRUE, real = TRUE, finite = TRUE))

# ---------------------------------------------------------------------------
# Integer parts

.round_fn <- function(name, rfun) {
  force(name)
  force(rfun)
  function(args) {
    a <- args[[1L]]
    t <- a$type
    if (t == "Integer") return(a)
    if (t == "Rational") {
      n <- gmp::numerator(a$value)
      d <- gmp::denominator(a$value)
      q <- gmp::divq.bigz(n, d)
      if (name == "ceiling" && n %% d != 0) q <- q + 1
      return(make_rational(gmp::as.bigq(q)))
    }
    if (t == "Float") return(.int(rfun(a$value)))
    if (is_inf_const(a) || is_nan(a)) return(a)
    if (isTRUE(ask_fact(a, "integer"))) return(a)
    if (t == "Function" && a$name %in% c("floor", "ceiling") &&
        isTRUE(ask_fact(a$args[[1L]], "real"))) return(a)
    if (is_const(a, "I")) return(a)
    if (t == "Add") {
      ints <- list()
      rest <- list()
      for (z in a$args) {
        if (isTRUE(ask_fact(z, "integer"))) ints[[length(ints) + 1L]] <- z
        else rest[[length(rest) + 1L]] <- z
      }
      if (length(ints)) return(.add(c(ints, list(fcall(name, list(.add(rest)))))))
    }
    if (length(free_symbols(a)) == 0L) {
      # get_integer_part
      if (.contains_key(a, S_ComplexInfinity$key) && !identical(ask_fact(a, "finite"), FALSE))
        stop(structure(class = c("ValueError", "error", "condition"), list(message = "Cannot get integer part of Complex Infinity", call = NULL)))
      r <- .round_number_mp(name, a)
      if (!is.null(r)) return(r)
    }
    NULL
  }
}

# get_integer_part of evalf.py: the integer parts of the real and imaginary
# parts of a constant, from an evaluation at increasing precision; an argument
# within rounding error of an integer is settled by the sign of arg - k at
# higher precision. NULL when the precision runs out (PrecisionExhausted).
.round_number_mp <- function(name, a) {
  dv <- tryCatch(.neval(a, list()), error = function(e) NA)
  mag <- if (is.na(dv) || !is.finite(Mod(dv))) 64 else max(0, ceiling(log2(Mod(dv) + 1)))
  prec <- mag + 40
  v <- tryCatch(.mpeval_adaptive(a, prec, list()), error = function(e) NULL)
  if (is.null(v) || !isTRUE(attr(v, "converged"))) return(NULL)
  part <- function(x, ex) {
    r <- .mp_near_integer(name, ex, x, prec, list())
    if (!is.null(r)) return(.mp_floor(r$re))
    if (name == "floor") .mp_floor(x) else -.mp_floor(list(m = -x$m, e = x$e))
  }
  re <- .int(part(v$re, a))
  if (v$im$m == 0) return(re)
  .add(list(re, .mul(list(.int(part(v$im, fcall("im", list(a)))), S_I))))
}

def_function("floor", eval = .round_fn("floor", floor),
             numeric = floor, diff = function(args, i) S_Zero,
             facts = function(args, afs) {
               out <- c(commutative = TRUE)
               if (isTRUE(afs[[1L]][["real"]])) out["integer"] <- TRUE
               out
             },
             latex = function(a) paste0("\\left\\lfloor{", a[1], "}\\right\\rfloor"))
def_function("ceiling", eval = .round_fn("ceiling", ceiling),
             numeric = ceiling, diff = function(args, i) S_Zero,
             facts = function(args, afs) {
               out <- c(commutative = TRUE)
               if (isTRUE(afs[[1L]][["real"]])) out["integer"] <- TRUE
               out
             },
             latex = function(a) paste0("\\left\\lceil{", a[1], "}\\right\\rceil"))
def_function("frac",
             eval = function(args) {
               a <- args[[1L]]
               if (is_number_type(a) || isTRUE(ask_fact(a, "integer")))
                 return(.sub(a, fcall("floor", list(a))))
               NULL
             },
             numeric = function(x) x - floor(x),
             facts = function(args, afs) {
               out <- c(commutative = TRUE)
               if (isTRUE(afs[[1L]][["real"]])) out["nonnegative"] <- TRUE
               out
             })

# ---------------------------------------------------------------------------
# Min and Max

# StrPrinter._print_LatticeOp: arguments by default_sort_key
.print_lattice_op <- function(e) {
  args <- .key_sort(lapply(e$args, .sort_key), e$args)
  paste0(e$name, "(", paste(vapply(args, .str, ""), collapse = ", "), ")")
}

.minmax_eval <- function(is_max) {
  force(is_max)
  name <- if (is_max) "Max" else "Min"
  function(args) {
    flat <- list()
    for (a in args) {
      if (a$type == "Function" && a$name == name) flat <- c(flat, a$args)
      else flat[[length(flat) + 1L]] <- a
    }
    flat <- flat[!duplicated(vapply(flat, .k, ""))]
    # SymPy MinMaxBase: oo absorbs Max and is dropped from Min; -oo the reverse
    zero <- if (is_max) "oo" else "-oo"
    ident <- if (is_max) "-oo" else "oo"
    if (any(vapply(flat, is_const, TRUE, nm = zero))) return(new_constant(zero))
    flat <- Filter(function(a) !is_const(a, ident), flat)
    if (!length(flat)) return(new_constant(ident))
    keep <- rep(TRUE, length(flat))
    for (i in seq_along(flat)) {
      if (!keep[i]) next
      for (j in seq_along(flat)) {
        if (i == j || !keep[j]) next
        d <- .sub(flat[[i]], flat[[j]])
        fd <- .facts(d)
        # MinMaxBase._is_connected: Ge and Le decide by is_extended_nonnegative
        dom <- if (is_max) isTRUE(fd[["extended_nonnegative"]]) else isTRUE(fd[["extended_nonpositive"]])
        # comparable constants are compared numerically, as SymPy's MinMaxBase
        if (!dom && !length(free_symbols(d)) && isTRUE(fd[["real"]]))
          dom <- tryCatch(.is_true(if (is_max) Ge(flat[[i]], flat[[j]]) else Le(flat[[i]], flat[[j]])), error = function(e) FALSE)
        if (dom) keep[j] <- FALSE
      }
    }
    flat <- flat[keep]
    if (length(flat) == 1L) return(flat[[1L]])
    # Expr.__new__(cls, *ordered(args))
    new_sym("Function", .bool_print_order(flat), name = name)
  }
}

def_function("Max", nargs = NULL, eval = .minmax_eval(TRUE),
             str = function(e) .print_lattice_op(e),
             numeric = function(...) max(Re(c(...))),
             diff = function(args, i) {
               others <- args[-i]
               fcall("Heaviside", list(.sub(args[[i]], fcall("Max", others))))
             },
             facts = function(args, afs) {
               out <- c(commutative = TRUE)
               if (.all_true(afs, "real")) out["real"] <- TRUE
               if (.any_true(afs, "positive")) out["positive"] <- TRUE
               if (.any_true(afs, "nonnegative")) out["nonnegative"] <- TRUE
               if (.all_true(afs, "integer")) out["integer"] <- TRUE
               out
             })
def_function("Min", nargs = NULL, eval = .minmax_eval(FALSE),
             str = function(e) .print_lattice_op(e),
             numeric = function(...) min(Re(c(...))),
             diff = function(args, i) {
               others <- args[-i]
               fcall("Heaviside", list(.sub(fcall("Min", others), args[[i]])))
             },
             facts = function(args, afs) {
               out <- c(commutative = TRUE)
               if (.all_true(afs, "real")) out["real"] <- TRUE
               if (.any_true(afs, "negative")) out["negative"] <- TRUE
               if (.all_true(afs, "positive")) out["positive"] <- TRUE
               if (.all_true(afs, "integer")) out["integer"] <- TRUE
               out
             })

# ---------------------------------------------------------------------------
# Heaviside, DiracDelta, KroneckerDelta, LambertW

def_function("Heaviside", nargs = 1:2,
             eval = function(args) {
               a <- args[[1L]]
               h0 <- if (length(args) == 2L) args[[2L]] else S_Half
               fa <- .facts(a)
               if (isTRUE(fa[["extended_negative"]])) return(S_Zero)
               if (isTRUE(fa[["extended_positive"]])) return(S_One)
               if (isTRUE(fa[["zero"]])) return(h0)
               if (is_nan(a)) return(S_NaN)
               NULL
             },
             diff = function(args, i) fcall("DiracDelta", list(args[[1L]])),
             numeric = function(x, h0 = 0.5) ifelse(Re(x) > 0, 1, ifelse(Re(x) < 0, 0, Re(h0))),
             facts = function(args, afs) c(commutative = TRUE, real = TRUE,
                                           nonnegative = TRUE, finite = TRUE),
             latex = function(a) paste0("\\theta\\left(", a[1], "\\right)"))

def_function("DiracDelta", nargs = 1:2,
             eval = function(args) {
               a <- args[[1L]]
               fa <- .facts(a)
               if (isTRUE(fa[["nonzero"]])) return(S_Zero)
               if (is_nan(a)) return(S_NaN)
               if (length(args) == 2L && is_zero_num(args[[2L]])) return(fcall("DiracDelta", list(a)))
               # a product led by a negative number (args_cnc()[0] is -1)
               if (length(args) == 1L && a$type == "Mul" && is_number_type(.sympy_args(a)[[1L]]) &&
                   num_sign(.sympy_args(a)[[1L]]) < 0)
                 return(fcall("DiracDelta", list(.neg(a))))
               NULL
             },
             diff = function(args, i) {
               k <- if (length(args) == 2L) args[[2L]] else S_Zero
               fcall("DiracDelta", list(args[[1L]], .add(list(k, S_One))))
             },
             numeric = function(x, k = 0) ifelse(x == 0, Inf, 0),
             facts = function(args, afs) c(commutative = TRUE, real = TRUE),
             latex = function(a) paste0("\\delta\\left(", a[1], "\\right)"))

def_function("KroneckerDelta", nargs = 2:3,
             eval = function(args) {
               # an index outside delta_range gives zero
               if (length(args) == 3L) {
                 r <- args[[3L]]$args
                 pos <- function(z) isTRUE(.facts(z)[["positive"]])
                 if (pos(.sub(r[[1L]], args[[1L]])) || pos(.sub(r[[1L]], args[[2L]])) ||
                     pos(.sub(args[[1L]], r[[2L]])) || pos(.sub(args[[2L]], r[[2L]]))) return(S_Zero)
               }
               d <- .sub(args[[1L]], args[[2L]])
               if (is_zero_num(d)) return(S_One)
               fd <- .facts(d)
               if (isTRUE(fd[["nonzero"]])) return(S_Zero)
               if (isTRUE(fd[["zero"]])) return(S_One)
               if (isTRUE(ask_fact(args[[1L]], "integer")) && isTRUE(ask_fact(args[[2L]], "integer")) &&
                   identical(fd[["integer"]], TRUE) && is_number_type(d)) return(S_Zero)
               # indices on opposite sides of the Fermi level
               fl <- function(e, nm) isTRUE(e$extra$other[nm][[1L]])
               if ((fl(args[[1L]], "below_fermi") && fl(args[[2L]], "above_fermi")) ||
                   (fl(args[[2L]], "below_fermi") && fl(args[[1L]], "above_fermi"))) return(S_Zero)
               # canonical order: default_sort_key(j) < default_sort_key(i) swaps them
               if (.key_cmp(.sort_key(args[[2L]]), .sort_key(args[[1L]])) < 0)
                 return(fnode("KroneckerDelta", c(args[2:1], args[-(1:2)])))
               NULL
             },
             numeric = function(i, j) as.numeric(i == j),
             facts = function(args, afs) c(commutative = TRUE, integer = TRUE,
                                           nonnegative = TRUE),
             latex = function(a) paste0("\\delta_{", a[1], " ", a[2], "}"))

def_function("LambertW", nargs = 1:2,
             eval = function(args) {
               a <- args[[1L]]
               k <- if (length(args) == 2L) args[[2L]] else S_Zero
               if (length(args) == 2L && is_zero_num(k)) return(fcall("LambertW", list(a)))
               if (!is_zero_num(k)) return(NULL)
               if (is_zero_num(a)) return(S_Zero)
               if (is_const(a, "E")) return(S_One)
               if (identical(a$key, .neg(.inv(S_E))$key)) return(S_NegativeOne)
               if (identical(a$key, .mul(list(.rat(-1, 2), fcall("log", list(S_Two))))$key))
                 return(.neg(fcall("log", list(S_Two))))
               if (is_const(a, "oo")) return(S_Infinity)
               if (a$type == "Float") return(.float_eval("LambertW", args))
               NULL
             },
             diff = function(args, i) {
               x <- args[[1L]]
               w <- fcall("LambertW", args)
               .div(w, .mul(list(x, .add(list(S_One, w)))))
             },
             numeric = function(x, k = 0) .lambertw_num(x, k),
             facts = function(args, afs) {
               out <- c(commutative = TRUE)
               if (length(args) == 1L && isTRUE(afs[[1L]][["positive"]])) {
                 out["positive"] <- TRUE
               }
               out
             })

.lambertw_num <- function(x, k = 0) {
  k <- Re(k)
  if (!is.finite(k) || k != round(k)) return(NaN)
  xr <- Re(x)
  xi <- if (is.complex(x)) Im(x) else 0
  if (xi == 0) {
    if (k == 0 && xr >= -exp(-1)) return(.lambertw_real(xr, 0))
    if (k == -1 && xr < 0 && xr >= -exp(-1)) return(.lambertw_real(xr, -1))
  }
  # the complex branches, as mpmath.lambertw
  z <- complex(real = xr, imaginary = xi)
  if (z == 0) return(if (k == 0) 0 else -Inf)
  if ((k == 0 || k == -1) && Mod(z + exp(-1)) < 0.3) {
    p <- sqrt(2 * (exp(1) * z + 1))
    if (k == -1) p <- -p
    w <- -1 + p - p^2 / 3
  } else if (k == 0 && Mod(z) < 3) {
    w <- log(1 + z)
    if (!is.finite(Re(w)) || !is.finite(Im(w))) w <- complex(real = -1, imaginary = 1)
  } else {
    L1 <- log(z) + 2i * pi * k
    w <- L1 - log(L1)
  }
  for (i in 1:100) {
    ew <- exp(w)
    f <- w * ew - z
    wn <- w - f / (ew * (w + 1) - (w + 2) * f / (2 * w + 2))
    if (!is.finite(Re(wn)) || !is.finite(Im(wn))) return(NaN)
    if (Mod(wn - w) < 1e-15 * (1 + Mod(wn))) return(wn)
    w <- wn
  }
  w
}

# the real branches k = 0 (x >= -1/e) and k = -1 (-1/e <= x < 0)
.lambertw_real <- function(x, k) {
  if (k == -1) {
    w <- if (x < -0.25) -1 - sqrt(2 * (1 + exp(1) * x)) else log(-x) - log(-log(-x))
  } else {
    w <- if (x < 1) x else log(x)
  }
  for (i in 1:100) {
    ew <- exp(w)
    f <- w * ew - x
    wn <- w - f / (ew * (w + 1) - (w + 2) * f / (2 * w + 2))
    if (!is.finite(wn)) return(w)
    if (abs(wn - w) < 1e-15 * (1 + abs(wn))) return(wn)
    w <- wn
  }
  w
}

# Basic.is_comparable: a number with a real value that can be computed
.is_comparable <- function(e) {
  if (is_nan(e)) return(FALSE)
  if (is_number_type(e) || is_inf_const(e)) return(TRUE)
  if (identical(ask_fact(e, "extended_real"), FALSE) || length(free_symbols(e))) return(FALSE)
  v <- tryCatch(.neval(e, list()), error = function(err) NA)
  length(v) == 1L && !is.na(v) && is.finite(Re(v)) && abs(Im(v)) == 0
}
