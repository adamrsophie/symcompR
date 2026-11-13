# An emulation of mpmath's number model, the layer below its special
# functions (mpmath/ctx_mp_python.py, ctx_mp.py, ctx_base.py). A value
# ("mpv") has one of the kinds mpmath code meets:
#   int    a Python int, exact (field q, a bigq with denominator 1)
#   mpq    an MPQ rational, exact (field q)
#   float  a Python float (field d, an R double)
#   mpf    an mpf (field re, a list(m, e) as in evalf-mp.R)
#   mpc    an mpc (fields re and im)
# An mpf or mpc may instead be infinite or NaN (field sp: "inf", "-inf",
# "nan"). Arithmetic follows Python's coercions: int with int stays int,
# int / int is a float, an MPQ meeting an mpf is first rounded to the
# context precision, ints and floats convert to mpf exactly. Every mpf or
# mpc operation rounds to nearest at the context precision .MPX$prec, as
# mpmath's context does. Values print through format() for debugging.

.MPX <- new.env(parent = emptyenv())
.MPX$prec <- 53

.mv_new <- function(k, ...) structure(list(k = k, ...), class = "mpv")
.hx <- function(re, im = .mp_zero, cx = FALSE) {
  if (cx || im$m != 0) .mv_new("mpc", re = re, im = im) else .mv_new("mpf", re = re, im = .mp_zero)
}
.mv_mpf <- function(re) .mv_new("mpf", re = re, im = .mp_zero)
.mv_mpc <- function(re, im) .mv_new("mpc", re = re, im = im)
.mv_int <- function(n) .mv_new("int", q = gmp::as.bigq(gmp::as.bigz(n)))
.mv_q <- function(n, d = 1) {
  q <- gmp::as.bigq(n, d)
  .mv_new("mpq", q = q)
}
.mv_float <- function(x) .mv_new("float", d = as.numeric(x))
.mv_special <- function(sp, cx = FALSE) .mv_new(if (cx) "mpc" else "mpf", re = .mp_zero, im = .mp_zero, sp = sp)
.mv_inf <- function() .mv_special("inf")
.mv_ninf <- function() .mv_special("-inf")
.mv_nan <- function(cx = FALSE) .mv_special("nan", cx)
.hx_int <- function(n) .mv_int(n)
.hx_q <- function(n, d) .mv_q(n, d)
.hx_one <- function() .mv_mpf(.mpf(1, 0))
.hx_zero <- function() .mv_mpf(.mp_zero)
.hx_j <- function() .mv_mpc(.mp_zero, .mpf(1, 0))
.hx_pi <- function() .mv_mpf(.mp_pi(.MPX$prec))
.hx_eps <- function() .mv_mpf(.mpf(1, 1 - .MPX$prec))
.hx_euler <- function() .mv_mpf(.mp_euler(.MPX$prec))
.hx_from_mpc <- function(z, cx = NULL) .hx(z$re, z$im, if (is.null(cx)) z$im$m != 0 else cx)
.hx_from_dbl <- function(x) .mv_float(x)

# Conditions with mpmath's exception names; .sympy_condition gives them the
# classes the battery compares (ValueError, ZeroDivisionError, ...).
.mp_error <- function(cls, msg) stop(.mp_condition(cls, msg))
.mp_condition <- function(cls, msg) {
  cond <- .sympy_condition(cls, msg)
  class(cond) <- c(class(cond)[1L], "mpmath_error", class(cond)[-1L])
  cond
}
.hx_noconv <- function(msg = "") .mp_error("NoConvergence", msg)
.hx_zerodiv <- function(msg = "pole in hypergeometric series") .mp_error("ZeroDivisionError", msg)
.hx_valueerror <- function(msg) .mp_error("ValueError", msg)
.hx_is_noconv <- function(err) inherits(err, "NoConvergence")

# ---------------------------------------------------------------------------
# Kinds and conversion

.mv_is <- function(x) inherits(x, "mpv")
# R numbers in ported code stand for Python literals: whole numbers are ints,
# others floats
.mv <- function(x) {
  if (inherits(x, "mpv")) return(x)
  if (inherits(x, "bigz")) return(.mv_int(x))
  if (inherits(x, "bigq")) return(.mv_new("mpq", q = x))
  if (is.numeric(x)) {
    if (is.finite(x) && x == round(x) && abs(x) < 2^53) return(.mv_int(x))
    return(.mv_float(x))
  }
  if (is.complex(x)) return(.mv_mpc(.mp_from_dbl(Re(x)), .mp_from_dbl(Im(x))))
  stop("not an mpmath value")
}
.mv_exact <- function(x) x$k %in% c("int", "mpq")
.mv_isnum <- function(x) x$k %in% c("mpf", "mpc")
.mv_cx <- function(x) x$k == "mpc"
.mv_sp <- function(x) !is.null(x$sp)

# ctx.convert: an mpf or mpc; ints and floats exactly, MPQ rounded to the
# context precision
.mv_convert <- function(x) {
  x <- .mv(x)
  switch(x$k,
    int = .mv_mpf(.mp_round(gmp::numerator(x$q), 0, .bitlen(gmp::numerator(x$q)) + 1)),
    mpq = .mv_mpf(.mp_from_q(x$q, .MPX$prec)),
    float = if (!is.finite(x$d)) .mv_special(if (is.nan(x$d)) "nan" else if (x$d > 0) "inf" else "-inf") else .mv_mpf(.mp_from_dbl(x$d)),
    x)
}
.hx_convert <- .mv_convert
.mv_as_mpc <- function(x) { x <- .mv_convert(x); if (x$k == "mpc") x else .mv_new("mpc", re = x$re, im = .mp_zero, sp = x$sp) }

# Python float value of a number (int / int division, float literals)
.mv_dbl <- function(x) {
  x <- .mv(x)
  switch(x$k, float = x$d, int = , mpq = .bigq_dbl(x$q),
         if (.mv_sp(x)) switch(x$sp, inf = Inf, `-inf` = -Inf, NaN) else .mp_to_dbl(x$re))
}

# ---------------------------------------------------------------------------
# Rounding helpers

.mp_pos <- function(a, p) .mp_round(a$m, a$e, p)
.mv_pos <- function(x, p = .MPX$prec) {
  x <- .mv(x)
  if (!.mv_isnum(x) || .mv_sp(x)) return(x)
  if (x$k == "mpf") .mv_mpf(.mp_pos(x$re, p)) else .mv_mpc(.mp_pos(x$re, p), .mp_pos(x$im, p))
}
.hx_round <- function(a, p = .MPX$prec) .mv_pos(a, p)

# ---------------------------------------------------------------------------
# Arithmetic

.mv_special_arith <- function(op, a, b) {
  # minimal IEEE-like rules for infinities and NaN
  val <- function(x) if (.mv_sp(x)) switch(x$sp, inf = Inf, `-inf` = -Inf, NaN) else {
    if (x$k == "mpc") complex(real = .mp_to_dbl(x$re), imaginary = .mp_to_dbl(x$im)) else .mv_dbl(x)
  }
  r <- switch(op, "+" = val(a) + val(b), "-" = val(a) - val(b), "*" = val(a) * val(b), "/" = val(a) / val(b))
  cx <- a$k == "mpc" || b$k == "mpc"
  if (is.complex(r)) {
    if (is.nan(Re(r)) || is.nan(Im(r))) return(.mv_nan(TRUE))
    if (is.infinite(Re(r)) && Im(r) == 0) return(.mv_special(if (Re(r) > 0) "inf" else "-inf", TRUE))
    if (is.infinite(Re(r)) || is.infinite(Im(r))) return(.mv_special("inf", TRUE))
    return(.mv_mpc(.mp_from_dbl(Re(r)), .mp_from_dbl(Im(r))))
  }
  if (is.nan(r)) return(.mv_nan(cx))
  if (is.infinite(r)) return(.mv_special(if (r > 0) "inf" else "-inf", cx))
  if (cx) .mv_mpc(.mp_from_dbl(r), .mp_zero) else .mv_mpf(.mp_from_dbl(r))
}

.hx_add <- function(a, b) {
  a <- .mv(a); b <- .mv(b)
  if (.mv_exact(a) && .mv_exact(b)) {
    r <- a$q + b$q
    return(if (a$k == "int" && b$k == "int") .mv_new("int", q = r) else .mv_new("mpq", q = r))
  }
  if (a$k == "float" && b$k != "mpf" && b$k != "mpc") return(.mv_float(a$d + .mv_dbl(b)))
  if (b$k == "float" && a$k != "mpf" && a$k != "mpc") return(.mv_float(.mv_dbl(a) + b$d))
  a <- .mv_convert(a); b <- .mv_convert(b)
  if (.mv_sp(a) || .mv_sp(b)) return(.mv_special_arith("+", a, b))
  p <- .MPX$prec
  if (a$k == "mpf" && b$k == "mpf") return(.mv_mpf(.mp_add(a$re, b$re, p)))
  # mpc_add_mpf leaves the imaginary part of the complex operand unrounded
  if (b$k == "mpf") return(.mv_mpc(.mp_add(a$re, b$re, p), a$im))
  if (a$k == "mpf") return(.mv_mpc(.mp_add(a$re, b$re, p), b$im))
  .mv_mpc(.mp_add(a$re, b$re, p), .mp_add(a$im, b$im, p))
}
.hx_neg <- function(a) {
  a <- .mv(a)
  switch(a$k, int = , mpq = .mv_new(a$k, q = -a$q), float = .mv_float(-a$d),
    if (.mv_sp(a)) .mv_special_arith("*", a, .mv_int(-1))
    else if (a$k == "mpf") .mv_mpf(.mp_pos(.mp_neg(a$re), .MPX$prec))
    else .mv_mpc(.mp_pos(.mp_neg(a$re), .MPX$prec), .mp_pos(.mp_neg(a$im), .MPX$prec)))
}
# fneg(x, exact=True)
.hx_neg_exact <- function(a) {
  a <- .mv(a)
  if (a$k %in% c("mpf", "mpc") && !.mv_sp(a)) return(.mv_new(a$k, re = .mp_neg(a$re), im = .mp_neg(a$im)))
  .hx_neg(a)
}
.hx_sub <- function(a, b) {
  a <- .mv(a); b <- .mv(b)
  if (.mv_exact(a) && .mv_exact(b)) {
    r <- a$q - b$q
    return(if (a$k == "int" && b$k == "int") .mv_new("int", q = r) else .mv_new("mpq", q = r))
  }
  if (a$k == "float" && b$k != "mpf" && b$k != "mpc") return(.mv_float(a$d - .mv_dbl(b)))
  if (b$k == "float" && a$k != "mpf" && a$k != "mpc") return(.mv_float(.mv_dbl(a) - b$d))
  a <- .mv_convert(a); b <- .mv_convert(b)
  if (.mv_sp(a) || .mv_sp(b)) return(.mv_special_arith("-", a, b))
  p <- .MPX$prec
  if (a$k == "mpf" && b$k == "mpf") return(.mv_mpf(.mp_sub(a$re, b$re, p)))
  if (b$k == "mpf") return(.mv_mpc(.mp_sub(a$re, b$re, p), a$im))           # mpc_sub_mpf
  if (a$k == "mpf") return(.mv_mpc(.mp_sub(a$re, b$re, p), .mp_pos(.mp_neg(b$im), p)))  # mpc_mpf_sub
  .mv_mpc(.mp_sub(a$re, b$re, p), .mp_sub(a$im, b$im, p))
}
.hx_mul <- function(a, b) {
  a <- .mv(a); b <- .mv(b)
  if (.mv_exact(a) && .mv_exact(b)) {
    r <- a$q * b$q
    return(if (a$k == "int" && b$k == "int") .mv_new("int", q = r) else .mv_new("mpq", q = r))
  }
  if (a$k == "float" && b$k != "mpf" && b$k != "mpc") return(.mv_float(a$d * .mv_dbl(b)))
  if (b$k == "float" && a$k != "mpf" && a$k != "mpc") return(.mv_float(.mv_dbl(a) * b$d))
  a <- .mv_convert(a); b <- .mv_convert(b)
  if (.mv_sp(a) || .mv_sp(b)) return(.mv_special_arith("*", a, b))
  p <- .MPX$prec
  if (a$k == "mpf" && b$k == "mpf") return(.mv_mpf(.mp_mul(a$re, b$re, p)))
  if (b$k == "mpf") return(.mv_mpc(.mp_mul(a$re, b$re, p), .mp_mul(a$im, b$re, p)))
  if (a$k == "mpf") return(.mv_mpc(.mp_mul(a$re, b$re, p), .mp_mul(a$re, b$im, p)))
  .mv_mpc(.mp_sub(.mp_xmul(a$re, b$re), .mp_xmul(a$im, b$im), p),
          .mp_add(.mp_xmul(a$re, b$im), .mp_xmul(a$im, b$re), p))
}
.hx_div <- function(a, b) {
  a <- .mv(a); b <- .mv(b)
  if (.mv_exact(a) && .mv_exact(b)) {
    if (b$q == 0) .mp_error("ZeroDivisionError", "division by zero")
    if (a$k == "int" && b$k == "int") return(.mv_float(.bigq_dbl(a$q / b$q)))
    return(.mv_new("mpq", q = a$q / b$q))
  }
  if (a$k == "float" && b$k != "mpf" && b$k != "mpc") {
    if (.mv_dbl(b) == 0) .mp_error("ZeroDivisionError", "float division by zero")
    return(.mv_float(a$d / .mv_dbl(b)))
  }
  if (b$k == "float" && a$k != "mpf" && a$k != "mpc") {
    if (b$d == 0) .mp_error("ZeroDivisionError", "float division by zero")
    return(.mv_float(.mv_dbl(a) / b$d))
  }
  a <- .mv_convert(a); b <- .mv_convert(b)
  if (.mv_sp(a) || .mv_sp(b)) return(.mv_special_arith("/", a, b))
  p <- .MPX$prec
  if (b$k == "mpf") {
    if (b$re$m == 0) {
      # mpf division by zero raises ZeroDivisionError
      .mp_error("ZeroDivisionError", "")
    }
    if (a$k == "mpf") return(.mv_mpf(.mp_div(a$re, b$re, p)))
    return(.mv_mpc(.mp_div(a$re, b$re, p), .mp_div(a$im, b$re, p)))
  }
  if (b$re$m == 0 && b$im$m == 0) .mp_error("ZeroDivisionError", "")
  if (a$k == "mpf") {
    # mpc_mpf_div
    m <- .mp_add(.mp_xmul(b$re, b$re), .mp_xmul(b$im, b$im), p + 10)
    return(.mv_mpc(.mp_div(.mp_xmul(b$re, a$re), m, p), .mp_div(.mp_neg(.mp_xmul(b$im, a$re)), m, p)))
  }
  wp <- p + 10
  mag <- .mp_add(.mp_xmul(b$re, b$re), .mp_xmul(b$im, b$im), wp)
  t <- .mp_add(.mp_xmul(a$re, b$re), .mp_xmul(a$im, b$im), wp)
  u <- .mp_sub(.mp_xmul(a$im, b$re), .mp_xmul(a$re, b$im), wp)
  .mv_mpc(.mp_div(t, mag, p), .mp_div(u, mag, p))
}
.hx_ldexp <- function(a, k) {
  a <- .mv_convert(a)
  if (.mv_sp(a)) return(a)
  .mv_new(a$k, re = .mp_ldexp(a$re, k), im = if (a$k == "mpc") .mp_ldexp(a$im, k) else .mp_zero)
}
.hx_sum <- function(xs) { s <- .mv_int(0); for (x in xs) s <- .hx_add(s, x); s }
.hx_prod <- function(xs) { s <- .mv_int(1); for (x in xs) s <- .hx_mul(s, x); s }

# x ** y with Python's and mpmath's typing
.hx_pow <- function(a, b) {
  a <- .mv(a); b <- .mv(b)
  if (.mv_exact(a) && b$k == "int") {
    n <- gmp::numerator(b$q)
    if (n >= 0) return(.mv_new(a$k, q = a$q^n))
    if (a$k == "mpq") { if (a$q == 0) .mp_error("ZeroDivisionError", "Fraction(0, 1) ** -1"); return(.mv_new("mpq", q = gmp::as.bigq(1) / a$q^(-n))) }
    if (a$q == 0) .mp_error("ZeroDivisionError", "0.0 cannot be raised to a negative power")
    return(.mv_float(.bigq_dbl(gmp::as.bigq(1) / a$q^(-n))))
  }
  if ((a$k == "float" || .mv_exact(a)) && (b$k == "float" || .mv_exact(b)) && !(a$k == "mpq" && b$k == "mpq")) {
    ad <- .mv_dbl(a); bd <- .mv_dbl(b)
    if (ad < 0 && bd != round(bd)) {
      r <- as.complex(ad)^bd
      return(.mv_mpc(.mp_from_dbl(Re(r)), .mp_from_dbl(Im(r))))
    }
    return(.mv_float(ad^bd))
  }
  p <- .MPX$prec
  if (a$k %in% c("mpf", "mpc") && b$k == "int") {
    if (.mv_sp(a)) return(.mv_special_arith("*", a, .mv_int(1)))
    n <- gmp::numerator(b$q)
    if (a$k == "mpf") {
      if (a$re$m == 0 && n < 0) .mp_error("ZeroDivisionError", "")
      return(.mv_mpf(.mp_pow_int(a$re, n, p)))
    }
    r <- .mpc_pow_int(.mpc(a$re, a$im), n, p)
    return(.mv_mpc(r$re, r$im))
  }
  a <- .mv_convert(a); b <- .mv_convert(b)
  if (.mv_sp(a) || .mv_sp(b)) return(.mv_nan(a$k == "mpc" || b$k == "mpc"))
  if (a$k == "mpf" && b$k == "mpf") {
    s <- a$re; t <- .mp_round(b$re$m, b$re$e, .bitlen(b$re$m) + 1)
    if (t$m == 0) return(.mv_mpf(.mpf(1, 0)))
    if (!(s$m < 0 && t$e < 0)) {
      if (s$m == 0) {
        if (t$m < 0) .mp_error("ZeroDivisionError", "")
        return(.mv_mpf(.mp_zero))
      }
      return(.mv_mpf(.mp_pow_real(s, t, p)))
    }
    r <- .mpc_pow_mpf(.mpc(s), t, p)
    return(.mv_mpc(r$re, r$im))
  }
  z <- .mv_as_mpc(a)
  if (z$re$m == 0 && z$im$m == 0) {
    # 0 ** w with complex w is NaN unless Re(w) > 0 and Im(w) = 0
    if (b$k == "mpf" && b$re$m > 0) return(.mv_mpc(.mp_zero, .mp_zero))
    if (b$k == "mpf" && b$re$m == 0) return(.mv_mpc(.mpf(1, 0), .mp_zero))
    return(.mv_nan(TRUE))
  }
  if (b$k == "mpf") {
    r <- .mpc_pow_mpf(.mpc(z$re, z$im), b$re, p)
  } else {
    l <- .mpc_log(.mpc(z$re, z$im), p + 10)
    r <- .mpc_exp(.mpc_mul(l, .mpc(b$re, b$im), p + 10), p)
  }
  .mv_mpc(r$re, r$im)
}
# mpf_pow for real s, t (s >= 0 or t an integer or half integer handled)
.mp_pow_real <- function(s, t, p) {
  if (t$e >= 0) return(.mp_pow_int(s, t$m * .pow2(t$e), p))
  if (t$e == -1) {
    if (abs(t$m) == 1) {
      if (t$m < 0) return(.mp_div(.mpf(1, 0), .mp_sqrt(s, p + 10), p))
      return(.mp_sqrt(s, p))
    }
    return(.mp_pow_int(.mp_sqrt(s, p + 10), t$m, p))
  }
  if (.mp_cmp(s, .mpf(1, 0)) == 0) return(.mpf(1, 0))
  .mp_exp(.mp_xmul(t, .mp_log(s, p + 10)), p)
}
.mpc_pow_mpf <- function(z, t, p) {
  t <- .mp_round(t$m, t$e, .bitlen(t$m) + 1)
  if (t$m == 0) return(.mpc_one)
  if (t$e >= 0) return(.mpc_pow_int(z, t$m * .pow2(t$e), p))
  if (t$e == -1) return(.mpc_pow_int(.mpc_sqrt(z, p + 10), t$m, p))
  l <- .mpc_log(z, p + 10)
  .mpc_exp(.mpc(.mp_mul(l$re, t, p + 10), .mp_mul(l$im, t, p + 10)), p)
}

# ---------------------------------------------------------------------------
# Comparisons (exact, as mpf_cmp; an MPQ against an mpf is rounded first)

.mv_cmp <- function(a, b) {
  a <- .mv(a); b <- .mv(b)
  if (.mv_exact(a) && .mv_exact(b)) return(sign(as.numeric(sign(a$q - b$q))))
  if (a$k == "float" && b$k == "float") return(sign(a$d - b$d))
  if (a$k == "mpc" || b$k == "mpc") stop(.mp_condition("TypeError", "no ordering relation is defined for complex numbers"))
  if (.mv_sp(a) || .mv_sp(b)) {
    d <- .mv_dbl(a) - .mv_dbl(b)
    if (is.nan(d)) return(NA_real_)
    return(sign(d))
  }
  if (a$k == "mpf" && b$k %in% c("int", "float")) return(.mp_cmp(a$re, .mv_convert(b)$re))
  if (b$k == "mpf" && a$k %in% c("int", "float")) return(.mp_cmp(.mv_convert(a)$re, b$re))
  if (.mv_exact(a) && b$k == "float") return(sign(as.numeric(sign(a$q - gmp::as.bigq(.mp_to_q(.mp_from_dbl(b$d)))))))
  if (a$k == "float" && .mv_exact(b)) return(-.mv_cmp(b, a))
  .mp_cmp(.mv_convert(a)$re, .mv_convert(b)$re)
}
.mv_eq <- function(a, b) {
  a <- .mv(a); b <- .mv(b)
  if (a$k == "mpc" || b$k == "mpc") {
    a <- .mv_as_mpc(a); b <- .mv_as_mpc(b)
    if (.mv_sp(a) || .mv_sp(b)) return(identical(a$sp, b$sp) && !identical(a$sp, "nan"))
    return(.mp_cmp(a$re, b$re) == 0 && .mp_cmp(a$im, b$im) == 0)
  }
  r <- .mv_cmp(a, b)
  !is.na(r) && r == 0
}
.hx_numeq <- .mv_eq

# ---------------------------------------------------------------------------
# Predicates and parts

.hx_iszero <- function(x) {
  x <- .mv(x)
  switch(x$k, int = , mpq = x$q == 0, float = x$d == 0,
         !.mv_sp(x) && x$re$m == 0 && x$im$m == 0)
}
.mv_bool <- function(x) !.hx_iszero(x)
.hx_isnan <- function(x) { x <- .mv(x); (x$k == "float" && is.nan(x$d)) || identical(x$sp, "nan") }
.hx_isinf <- function(x) { x <- .mv(x); (x$k == "float" && is.infinite(x$d)) || (.mv_sp(x) && x$sp != "nan") }
.hx_isfinite <- function(x) !.hx_isinf(x) && !.hx_isnan(x)
.mp_norm <- function(a) if (a$m == 0) .mp_zero else .mp_round(a$m, a$e, .bitlen(a$m) + 1)
.mp_isint <- function(a) a$m == 0 || .mp_norm(a)$e >= 0
# ctx.isint: integer valued (an mpc only when its imaginary part is zero)
.hx_isint <- function(x) {
  x <- .mv(x)
  switch(x$k, int = TRUE, mpq = gmp::denominator(x$q) == 1, float = is.finite(x$d) && x$d == round(x$d),
    if (.mv_sp(x)) FALSE else .mp_isint(x$re) && x$im$m == 0)
}
# ctx.isnpint: a nonpositive integer
.hx_isnpint <- function(x) {
  x <- .mv(x)
  if (.hx_iszero(x)) return(TRUE)
  switch(x$k, int = x$q <= 0, mpq = gmp::denominator(x$q) == 1 && x$q <= 0, float = is.finite(x$d) && x$d < 0 && x$d == round(x$d),
    if (.mv_sp(x)) FALSE else x$im$m == 0 && x$re$m < 0 && .mp_isint(x$re))
}
.hx_int_value <- function(x) {
  x <- .mv(x)
  if (.mv_exact(x)) return(gmp::numerator(x$q))
  if (x$k == "float") return(gmp::as.bigz(x$d))
  .mp_floor(x$re)
}
# _is_real_type
.hx_is_real_type <- function(x) .mv(x)$k != "mpc"
.hx_is_complex_type <- function(x) .mv(x)$k == "mpc"
.hx_re <- function(x) {
  x <- .mv(x)
  if (x$k != "mpc") return(x)
  if (.mv_sp(x)) return(.mv_special(x$sp))
  .mv_mpf(x$re)
}
.hx_im <- function(x) {
  x <- .mv(x)
  if (x$k != "mpc") return(if (x$k %in% c("int", "mpq", "float")) .mv_int(0) else .mv_mpf(.mp_zero))
  if (.mv_sp(x)) return(.mv_special(x$sp))
  .mv_mpf(x$im)
}
# ctx.re / ctx.im return mpf zero for real input as well
.hx_im0 <- function(x) { x <- .mv_convert(x); if (x$k == "mpc") .mv_mpf(x$im) else .mv_mpf(.mp_zero) }
.hx_conj <- function(x) { x <- .mv(x); if (x$k != "mpc") x else .mv_mpc(x$re, .mp_neg(x$im)) }

# abs(x) rounds, as mpf_abs and mpc_abs
.hx_abs <- function(a) {
  a <- .mv(a)
  switch(a$k, int = , mpq = .mv_new(a$k, q = abs(a$q)), float = .mv_float(abs(a$d)),
    if (.mv_sp(a)) (if (a$sp == "nan") .mv_nan() else .mv_inf())
    else if (a$k == "mpf") .mv_mpf(.mp_pos(list(m = abs(a$re$m), e = a$re$e), .MPX$prec))
    else .mv_mpf(.mpc_abs(.mpc(a$re, a$im), .MPX$prec)))
}
.hx_sign <- function(x) {
  x <- .mv(x)
  if (x$k == "mpc") stop("sign of a complex number")
  if (.mv_exact(x)) return(as.numeric(sign(x$q)))
  if (x$k == "float") return(sign(x$d))
  if (.mv_sp(x)) return(switch(x$sp, inf = 1, `-inf` = -1, NaN))
  .mp_sign(x$re)
}

# ctx.mag: |x| <= 2^mag
.hx_mag <- function(a) {
  a <- .mv(a)
  if (a$k == "int") { n <- abs(gmp::numerator(a$q)); return(if (n == 0) -Inf else .bitlen(n)) }
  if (a$k == "mpq") {
    n <- gmp::numerator(a$q); d <- gmp::denominator(a$q)
    return(if (n == 0) -Inf else 1 + .bitlen(abs(n)) - .bitlen(d))
  }
  a <- .mv_convert(a)
  if (.mv_sp(a)) return(if (a$sp == "nan") NaN else Inf)
  if (a$re$m == 0 && a$im$m == 0) return(-Inf)
  if (a$re$m == 0) return(.mp_mag(a$im))
  if (a$im$m == 0) return(.mp_mag(a$re))
  1 + max(.mp_mag(a$re), .mp_mag(a$im))
}

# ctx.nint_distance: the nearest integer n to x and an estimate d of
# log2 |x - n| (MPQ exactly)
.hx_nint_distance <- function(x) {
  x <- .mv(x)
  if (x$k == "int") return(list(n = gmp::numerator(x$q), d = -Inf))
  if (x$k == "mpq") {
    p <- gmp::numerator(x$q); q <- gmp::denominator(x$q)
    n <- gmp::divq.bigz(p, q); r <- p - n * q
    if (2 * r >= q) n <- n + 1
    else if (r == 0) return(list(n = n, d = -Inf))
    return(list(n = n, d = .bitlen(abs(p - n * q)) - .bitlen(q)))
  }
  x <- .mv_convert(x)
  im_dist <- if (x$k == "mpc" && x$im$m != 0) .mp_mag(.mp_norm(x$im)) else -Inf
  re <- .mp_norm(x$re)
  if (re$m == 0) return(list(n = gmp::as.bigz(0), d = im_dist))
  man <- re$m; ex <- re$e
  mag <- ex + .bitlen(man)
  if (mag < 0) return(list(n = gmp::as.bigz(0), d = max(mag, im_dist)))
  neg <- man < 0; man <- abs(man)
  if (ex >= 0) {
    n <- man * .pow2(ex); re_dist <- -Inf
  } else if (ex == -1) {
    n <- gmp::divq.bigz(man, 2) + 1; re_dist <- 0
  } else {
    d <- -ex - 1
    t <- gmp::divq.bigz(man, .pow2(d))
    if (gmp::mod.bigz(t, 2) == 1) {
      t <- t + 1
      man <- t * .pow2(d) - man
    } else man <- man - t * .pow2(d)
    n <- gmp::divq.bigz(t, 2)
    re_dist <- ex + .bitlen(man)
  }
  if (neg) n <- -n
  list(n = n, d = max(re_dist, im_dist))
}

# ---------------------------------------------------------------------------
# Sums and products

# fsum: exact accumulation of the parts, one rounding (mpf_sum)
.hx_fsum <- function(xs) {
  p <- .MPX$prec
  re <- list(); im <- list(); cx <- FALSE
  for (x in xs) {
    x <- .mv_convert(x)
    if (.mv_sp(x)) return(Reduce(.hx_add, xs))
    re[[length(re) + 1L]] <- x$re
    if (x$k == "mpc") { cx <- TRUE; if (x$im$m != 0) im[[length(im) + 1L]] <- x$im }
  }
  s <- .mp_sum(re, p)
  if (length(im) || cx) {
    if (!length(im)) return(.mv_mpc(s, .mp_zero))
    return(.mv_mpc(s, .mp_sum(im, p)))
  }
  .mv_mpf(s)
}
# mpf_sum: exact sum of the terms whose magnitude lies above the rounding
# position of the largest (smaller ones decide only the rounding), then
# rounded once
# mpf_sum: exact accumulation, except that a term more than 2*prec bits
# away from the running sum is dropped (or replaces it), one final rounding
.mp_sum <- function(xs, p) {
  max_extra_prec <- if (p) 2 * p else 1000000
  man <- gmp::as.bigz(0); ex <- 0
  for (x in xs) {
    if (x$m == 0) next
    xman <- x$m; xexp <- x$e
    delta <- xexp - ex
    if (xexp >= ex) {
      if (delta > max_extra_prec && (man == 0 || delta - .bitlen(abs(man)) > max_extra_prec)) {
        man <- xman; ex <- xexp
      } else {
        man <- man + xman * .pow2(delta)
      }
    } else {
      delta <- -delta
      if (delta - .bitlen(abs(xman)) > max_extra_prec) {
        if (man == 0) { man <- xman; ex <- xexp }
      } else {
        man <- man * .pow2(delta) + xman
        ex <- xexp
      }
    }
  }
  .mp_round(man, ex, p)
}
# fprod: the product at the context precision, rounded at the end
.hx_fprod <- function(xs) { v <- .hx_one(); for (x in xs) v <- .hx_mul(v, x); .mv_pos(v) }

# sum_accurately(terms): adds terms with increasing precision until the
# cancellation is covered; terms is a function returning a list
.hx_sum_accurately <- function(terms, check_step = 1) {
  prec <- .MPX$prec
  on.exit(.MPX$prec <- prec)
  extraprec <- 10
  repeat {
    .MPX$prec <- prec + extraprec + 5
    max_mag <- -Inf; sum_mag <- 0
    s <- .hx_zero()
    k <- 0
    for (term in terms()) {
      s <- .hx_add(s, term)
      if (k %% check_step == 0 && !.hx_iszero(term)) {
        term_mag <- .hx_mag(term)
        max_mag <- max(max_mag, term_mag)
        sum_mag <- .hx_mag(s)
        if (sum_mag - term_mag > .MPX$prec) break
      }
      k <- k + 1
    }
    cancellation <- max_mag - sum_mag
    if (is.nan(cancellation)) break
    if (cancellation < extraprec) break
    extraprec <- extraprec + min(.MPX$prec, cancellation)
  }
  s
}

# ---------------------------------------------------------------------------
# Elementary functions at the context precision

.hx_exp <- function(a) {
  a <- .mv_convert(a); p <- .MPX$prec
  if (.mv_sp(a)) return(switch(a$sp, inf = .mv_inf(), `-inf` = .hx_zero(), .mv_nan()))
  if (a$k == "mpf") return(.mv_mpf(.mp_exp(a$re, p)))
  r <- .mpc_exp(.mpc(a$re, a$im), p); .mv_mpc(r$re, r$im)
}
.hx_log <- function(a) {
  a <- .mv_convert(a); p <- .MPX$prec
  if (.mv_sp(a)) return(if (a$sp == "inf") .mv_inf() else .mv_nan())
  if (a$k == "mpf") {
    if (a$re$m == 0) return(.mv_ninf())
    if (a$re$m > 0) return(.mv_mpf(.mp_log(a$re, p)))
  }
  if (a$re$m == 0 && a$im$m == 0) return(.mv_ninf())
  r <- .mpc_log(.mpc(a$re, a$im), p); .mv_mpc(r$re, r$im)
}
.hx_sqrt <- function(a) {
  a <- .mv_convert(a); p <- .MPX$prec
  if (.mv_sp(a)) return(a)
  if (a$k == "mpf") {
    if (a$re$m >= 0) return(.mv_mpf(.mp_sqrt(a$re, p)))
    return(.mv_mpc(.mp_zero, .mp_sqrt(.mp_neg(a$re), p)))
  }
  r <- .mpc_sqrt(.mpc(a$re, a$im), p); .mv_mpc(r$re, r$im)
}
.hx_cos_sin <- function(a) {
  a <- .mv_convert(a); p <- .MPX$prec
  if (a$k == "mpf") { cs <- .mp_cos_sin(a$re, p); return(list(cos = .mv_mpf(cs$cos), sin = .mv_mpf(cs$sin))) }
  z <- .mpc(a$re, a$im)
  list(cos = .hx_from_mpc(.mpc_cos(z, p), TRUE), sin = .hx_from_mpc(.mpc_sin(z, p), TRUE))
}
.hx_sin <- function(a) .hx_cos_sin(a)$sin
.hx_cos <- function(a) .hx_cos_sin(a)$cos
.hx_cospi_sinpi <- function(a) {
  a <- .mv_convert(a); p <- .MPX$prec
  if (a$k == "mpf") { cs <- .mp_cos_sin(a$re, p, pi = TRUE); return(list(cos = .mv_mpf(cs$cos), sin = .mv_mpf(cs$sin))) }
  cs <- .mpc_cos_sin_pi(.mpc(a$re, a$im), p)
  list(cos = .hx_from_mpc(cs$cos, TRUE), sin = .hx_from_mpc(cs$sin, TRUE))
}
.hx_sinpi <- function(a) .hx_cospi_sinpi(a)$sin
.hx_cospi <- function(a) .hx_cospi_sinpi(a)$cos
# exp(i x) and exp(i pi x), as mpc_expj and mpc_expjpi (complex results)
.hx_expj <- function(a) {
  a <- .mv_convert(a); p <- .MPX$prec
  re <- a$re; im <- if (a$k == "mpc") a$im else .mp_zero
  if (im$m == 0) { cs <- .mp_cos_sin(re, p); return(.mv_mpc(cs$cos, cs$sin)) }
  if (re$m == 0) return(.mv_mpc(.mp_exp(.mp_neg(im), p), .mp_zero))
  ey <- .mp_exp(.mp_neg(im), p + 10)
  cs <- .mp_cos_sin(re, p + 10)
  .mv_mpc(.mp_mul(ey, cs$cos, p), .mp_mul(ey, cs$sin, p))
}
.hx_expjpi <- function(a) {
  a <- .mv_convert(a); p <- .MPX$prec
  re <- a$re; im <- if (a$k == "mpc") a$im else .mp_zero
  if (im$m == 0) { cs <- .mp_cos_sin(re, p, pi = TRUE); return(.mv_mpc(cs$cos, cs$sin)) }
  wp <- p + 10 + max(0, .mp_mag(im))
  im <- .mp_neg(.mp_mul(.mp_pi(wp), im, wp))
  if (re$m == 0) return(.mv_mpc(.mp_exp(im, p), .mp_zero))
  ey <- .mp_exp(im, p + 10)
  cs <- .mp_cos_sin(re, p + 10, pi = TRUE)
  .mv_mpc(.mp_mul(ey, cs$cos, p), .mp_mul(ey, cs$sin, p))
}
.hx_atan2 <- function(y, x) .mv_mpf(.mp_atan2(.mv_convert(y)$re, .mv_convert(x)$re, .MPX$prec))
.hx_arg <- function(z) { z <- .mv_as_mpc(z); .mv_mpf(.mp_atan2(z$im, z$re, .MPX$prec)) }

# gamma functions: a pole raises ValueError, as mpmath's gamma
.hx_gamma <- function(a) {
  a <- .mv_convert(a); p <- .MPX$prec
  if (.hx_isnpint(a)) .hx_valueerror("gamma function pole")
  if (a$k == "mpf") return(.mv_mpf(.mp_gamma_r(a$re, p)))
  r <- .mpc_gamma_r(.mpc(a$re, a$im), p); .mv_mpc(r$re, r$im)
}
.hx_rgamma <- function(a) {
  a <- .mv_convert(a); p <- .MPX$prec
  if (.hx_isnpint(a)) return(if (a$k == "mpc") .mv_mpc(.mp_zero, .mp_zero) else .hx_zero())
  if (a$k == "mpf") return(.mv_mpf(.mp_gamma_r(a$re, p, TRUE)))
  r <- .mpc_gamma_r(.mpc(a$re, a$im), p, TRUE); .mv_mpc(r$re, r$im)
}

# gammaprod(a, b): prod gamma(a) / prod gamma(b) with poles cancelled; an
# excess of poles in the numerator gives an infinity (signed with infsign)
.hx_gammaprod <- function(a, b, infsign = FALSE) {
  a <- lapply(a, .mv_convert); b <- lapply(b, .mv_convert)
  pn <- Filter(.hx_isnpint, a); rn <- Filter(Negate(.hx_isnpint), a)
  pd <- Filter(.hx_isnpint, b); rd <- Filter(Negate(.hx_isnpint), b)
  if (length(pn) < length(pd)) return(.hx_zero())
  if (length(pn) > length(pd)) {
    if (!infsign) return(.mv_inf())
    s <- .hx_gammaprod(c(pn[seq_len(length(pd))], rn), c(pd, rd))
    return(if (.hx_sign(.hx_re(s)) < 0) .mv_ninf() else .mv_inf())
  }
  orig <- .MPX$prec
  on.exit(.MPX$prec <- orig)
  .MPX$prec <- orig + 15
  p <- .hx_one()
  while (length(pn)) {
    i <- pn[[length(pn)]]; pn[[length(pn)]] <- NULL
    j <- pd[[length(pd)]]; pd[[length(pd)]] <- NULL
    sg <- if (gmp::mod.bigz(.hx_int_value(i) + .hx_int_value(j), 2) == 0) 1 else -1
    p <- .hx_mul(p, .hx_div(.hx_mul(sg, .hx_gamma(.hx_sub(1, j))), .hx_gamma(.hx_sub(1, i))))
  }
  for (x in rn) p <- .hx_mul(p, .hx_gamma(x))
  for (x in rd) p <- .hx_div(p, .hx_gamma(x))
  .MPX$prec <- orig
  .mv_pos(p)
}

# ---------------------------------------------------------------------------
# Operators, so ported formulas read like the Python

Ops.mpv <- function(e1, e2) {
  if (missing(e2)) {
    return(switch(.Generic, "-" = .hx_neg(e1), "+" = .mv_pos(e1), "!" = !.mv_bool(e1),
                  stop("unsupported unary operator for mpmath values")))
  }
  switch(.Generic,
    "+" = .hx_add(e1, e2), "-" = .hx_sub(e1, e2), "*" = .hx_mul(e1, e2), "/" = .hx_div(e1, e2),
    "^" = .hx_pow(e1, e2),
    "==" = .mv_eq(e1, e2), "!=" = !.mv_eq(e1, e2),
    "<" = isTRUE(.mv_cmp(e1, e2) < 0), ">" = isTRUE(.mv_cmp(e1, e2) > 0),
    "<=" = isTRUE(.mv_cmp(e1, e2) <= 0), ">=" = isTRUE(.mv_cmp(e1, e2) >= 0),
    stop("unsupported operator for mpmath values: ", .Generic))
}
Math.mpv <- function(x, ...) {
  switch(.Generic, abs = .hx_abs(x), sqrt = .hx_sqrt(x), exp = .hx_exp(x), log = .hx_log(x),
         sin = .hx_sin(x), cos = .hx_cos(x), gamma = .hx_gamma(x),
         stop("unsupported function for mpmath values: ", .Generic))
}
format.mpv <- function(x, ...) {
  switch(x$k, int = paste0("int(", as.character(x$q), ")"), mpq = paste0("mpq(", as.character(x$q), ")"),
         float = paste0("float(", format(x$d, digits = 17), ")"),
         if (.mv_sp(x)) paste0(x$k, "(", x$sp, ")")
         else if (x$k == "mpf") paste0("mpf(", .mp_to_str(x$re, 20), ")")
         else paste0("mpc(", .mp_to_str(x$re, 20), ", ", .mp_to_str(x$im, 20), ")"))
}
print.mpv <- function(x, ...) { cat(format(x), "\n"); invisible(x) }
