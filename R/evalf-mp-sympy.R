# SymPy's evalf for functions it evaluates through mpmath. As
# Function._eval_evalf does, the arguments are converted at prec + 5
# (Integers stay exact Python ints), the mpmath function runs with the
# context precision prec, and Expr._from_mpmath rounds the result to prec.
# Classes with their own _eval_evalf (meijerg, lowergamma, ...) follow it.
# At the top level EvalfMixin.evalf evaluates once at prec + 4 and rounds the
# parts to prec.

# arg._to_mpmath(p): an mpmath value, or NULL when the argument has no value
.mp_sympy_arg <- function(x, p, env, allow_ints = TRUE) {
  x <- as_sym(x)
  t <- x$type
  if ((t == "Symbol" || t == "Dummy") && !is.null(env[[x$name]])) return(.mp_sympy_arg(env[[x$name]], p, env, allow_ints))
  if (t == "Integer" && allow_ints) return(.mv_int(gmp::numerator(x$value)))
  if (t == "Integer" || t == "Rational") return(.mv_mpf(.mp_from_q(gmp::as.bigq(x$value), p)))
  if (t == "Float") { f <- .float_mp(x); return(.mv_mpf(.mp_round(f$m, f$e, p))) }
  if (t == "Constant" && x$name == "I") return(.mv_mpc(.mp_zero, .mpf(1, 0)))
  if (is_inf_const(x) || is_nan(x)) {
    if (identical(x, S_Infinity)) return(.mv_inf())
    if (identical(x, S_NegativeInfinity)) return(.mv_ninf())
    return(NULL)
  }
  # h + c*I with numbers h, c (pure_complex): evalf_add and evalf_mul take
  # each part through evalf, where a Rational is rounded down (round_fast)
  # and a Float keeps all its bits
  pc <- .mp_pure_complex(x)
  if (!is.null(pc)) {
    part <- function(n) if (n$type == "Float") .float_mp(n) else .mp_from_q_down(n$value, p)
    return(.mv_mpc(part(pc[[1L]]), part(pc[[2L]])))
  }
  # evalf(x, p) as in Expr._to_mpmath (evalf-core.R)
  if (!length(env)) {
    v <- .ev_to_mpmath(x, p)
    if (!is.null(v)) return(v)
  }
  r <- .mpeval_adaptive(x, p, env)
  if (is.null(r)) return(NULL)
  if (r$im$m == 0) .mv_mpf(r$re) else .mv_mpc(r$re, r$im)
}
# SymPy's pure_complex: list(h, c) for h + c*I with numbers h and c (h may
# be zero), else NULL
.mp_pure_complex <- function(x) {
  isnum <- function(n) n$type %in% c("Integer", "Rational", "Float")
  imag <- function(t) {
    if (t$type == "Constant" && t$name == "I") return(S_One)
    if (t$type == "Mul" && length(t$args) == 2L) {
      for (k in 1:2) if (isnum(t$args[[k]]) && t$args[[3L - k]]$type == "Constant" && t$args[[3L - k]]$name == "I")
        return(t$args[[k]])
    }
    NULL
  }
  if (x$type == "Add" && length(x$args) == 2L) {
    for (k in 1:2) if (isnum(x$args[[k]]) && !is.null(c <- imag(x$args[[3L - k]]))) return(list(x$args[[k]], c))
    return(NULL)
  }
  if (x$type == "Mul") {
    c <- imag(x)
    return(if (is.null(c)) NULL else list(S_Zero, c))
  }
  NULL
}
# mpmath's from_rational(p, q, prec, round_down)
.mp_from_q_down <- function(q, p) {
  q <- gmp::as.bigq(q)
  n <- gmp::numerator(q); d <- gmp::denominator(q)
  if (n == 0) return(.mp_zero)
  neg <- n < 0; n <- abs(n)
  k <- p + .bitlen(d) - .bitlen(n) + 4
  m <- if (k >= 0) gmp::divq.bigz(n * .pow2(k), d) else gmp::divq.bigz(n, d * .pow2(-k))
  b <- .bitlen(m)
  if (b > p) { m <- gmp::divq.bigz(m, .pow2(b - p)); k <- k - (b - p) }
  .mp_round(if (neg) -m else m, -k, Inf)
}
.mp_sympy_args <- function(xs, p, env) {
  out <- lapply(xs, .mp_sympy_arg, p, env)
  if (any(vapply(out, is.null, NA))) NULL else out
}

# The registry: SymPy function name to an evaluator taking the node, the
# binary precision and the environment, returning an mpmath value (or NULL)
.MP_SYMPY <- new.env(parent = emptyenv())

.mp_reg <- function(names, f) for (nm in names) assign(nm, f, envir = .MP_SYMPY)

# Function._eval_evalf with mpmath function f of the converted arguments
.mp_generic <- function(f, nargs = NULL) {
  force(f)
  function(e, p, env) {
    if (!is.null(nargs) && !(length(e$args) %in% nargs)) return(NULL)
    args <- .mp_sympy_args(e$args, p + 5, env)
    if (is.null(args)) return(NULL)
    .mp_workprec(p, do.call(f, args))
  }
}
.mp_workprec <- function(p, expr) {
  old <- .MPX$prec
  .MPX$prec <- p
  on.exit(.MPX$prec <- old)
  force(expr)
}

# hyper: Function._eval_evalf with Tuple arguments
.mp_reg("hyper", function(e, p, env) {
  ap <- .mp_sympy_args(e$args[[1L]]$args, p + 5, env)
  bq <- .mp_sympy_args(e$args[[2L]]$args, p + 5, env)
  z <- .mp_sympy_arg(e$args[[3L]], p + 5, env)
  if (is.null(ap) || is.null(bq) || is.null(z)) return(NULL)
  .mp_workprec(p, .hx_hyper(ap, bq, z))
})

# meijerg._eval_evalf: a polar argument z*exp_polar(I*b) with |b| > pi is
# evaluated as G(w**n) with w = z**(1/n)*exp(I*b/n) and n = ceiling(|b/pi|) + 1
.mp_reg("meijerg", function(e, p, env) {
  zarg <- e$args[[3L]]
  for (s in free_symbols(zarg)) if (!is.null(env[[s$name]])) zarg <- subs(zarg, s, env[[s$name]])
  fs <- if (zarg$type == "Mul") zarg$args else list(zarg)
  coeff <- list(); branch <- list()
  for (f in fs) {
    if (f$type == "Function" && f$name == "exp_polar") {
      # exp_polar._eval_evalf: an argument with -pi < Im <= pi becomes exp
      iv <- .mpeval_adaptive(f$args[[1L]], p + 20, env)
      if (is.null(iv)) return(NULL)
      i <- iv$im
      pi_ <- .mp_pi(p + 20)
      bad <- .mp_cmp(i, .mp_neg(pi_)) <= 0 || .mp_cmp(i, pi_) > 0
      if (bad) branch[[length(branch) + 1L]] <- f else coeff[[length(coeff) + 1L]] <- fcall("exp", f$args)
    } else coeff[[length(coeff) + 1L]] <- f
  }
  if (length(branch) > 1L) return(NULL)
  n <- 1
  if (length(branch)) {
    bs <- .mul(list(.neg(S_I), branch[[1L]]$args[[1L]]))
    b <- .mpeval_adaptive(bs, p + 20, env)
    if (is.null(b)) return(NULL)
    # n = ceiling(abs(branch/pi)) + 1, exactly when branch/pi is rational;
    # a complex branch (exp_polar of a non imaginary argument) uses |branch/pi|
    ratio <- .mul(list(bs, .pow(S_Pi, S_NegativeOne)))
    if (ratio$type %in% c("Integer", "Rational")) {
      q <- abs(gmp::as.bigq(ratio$value))
      n <- as.numeric(-gmp::divq.bigz(-gmp::numerator(q), gmp::denominator(q))) + 1
    } else {
      rv <- .mpeval_adaptive(ratio, p + 20, env)
      n <- ceiling(Mod(complex(real = .mp_to_dbl(rv$re), imaginary = .mp_to_dbl(rv$im)))) + 1
    }
  }
  c <- .mp_sympy_arg(if (length(coeff)) .mul(coeff) else S_One, p, env)
  if (is.null(c)) return(NULL)
  # znum = argument._eval_evalf(prec) holds prec bits
  c <- .mv_pos(c, p)
  if (n > 1 && b$im$m != 0) {
    # a complex branch: the expression znum**(1/n)*exp(I*branch/n) itself
    # goes through _to_mpmath(prec), that is evalf at prec
    zs <- .mul(list(.pow(if (length(coeff)) .mul(coeff) else S_One, .rat(1, n)),
                    fcall("exp", list(.mul(list(S_I, bs, .rat(1, n)))))))
    z <- .mp_sympy_arg(zs, p, env)
    if (is.null(z)) return(NULL)
    r <- .mv_mpf(.mp_from_q(gmp::as.bigq(1, n), p))
  } else if (n > 1) {
    # znum**(1/n)*exp(I*branch/n), then _to_mpmath(prec)
    old <- .MPX$prec; .MPX$prec <- p; on.exit(.MPX$prec <- old)
    w <- .hx_pow(c, .mv_mpf(.mp_from_q(gmp::as.bigq(1, n), p)))
    .MPX$prec <- p + 20
    w <- .hx_mul(w, .hx_expj(.hx_div(.mv_mpf(b$re), n)))
    .MPX$prec <- old
    z <- .mv_pos(w, p)
    r <- .mv_mpf(.mp_from_q(gmp::as.bigq(1, n), p))
  } else {
    z <- c; r <- .mv_int(1)
  }
  A <- e$args[[1L]]$args; B <- e$args[[2L]]$args
  an <- .mp_sympy_args(A[[1L]]$args, p, env); ap <- .mp_sympy_args(A[[2L]]$args, p, env)
  bm <- .mp_sympy_args(B[[1L]]$args, p, env); bq <- .mp_sympy_args(B[[2L]]$args, p, env)
  if (is.null(an) || is.null(ap) || is.null(bm) || is.null(bq)) return(NULL)
  .mp_workprec(p, .hx_meijerg(an, ap, bm, bq, z, r))
})

# Called from .mpeval: FALSE when SymPy does not evaluate the function
# through mpmath, else an mpc at precision p or NULL
.mp_sympy_fun <- function(e, p, env) {
  f <- .MP_SYMPY[[e$name]]
  if (is.null(f)) return(FALSE)
  v <- f(e, p, env)
  if (is.null(v)) return(NULL)
  v <- .mv(v)
  if (.mv_sp(v) || !.mv_isnum(v)) {
    if (!.mv_isnum(v) && !.mv_sp(v)) v <- .mv_convert(v) else return(NULL)
  }
  if (v$k == "mpf") .mpc(.mp_round(v$re$m, v$re$e, p)) else .mpc(.mp_round(v$re$m, v$re$e, p), .mp_round(v$im$m, v$im$e, p))
}

# A value from Expr._from_mpmath as a sym: infinities and NaN included
.mp_value_sym <- function(v, prec) {
  v <- .mv(v)
  if (.mv_sp(v)) return(switch(v$sp, inf = S_Infinity, `-inf` = S_NegativeInfinity, S_NaN))
  if (!.mv_isnum(v)) v <- .mv_convert(v)
  re <- .mp_round(v$re$m, v$re$e, prec)
  im <- if (v$k == "mpc") .mp_round(v$im$m, v$im$e, prec) else .mp_zero
  .mpc_sym(.mpc(re, im), prec, zero_int = TRUE)
}

# evalf_table: Pi and Exp1 evaluate to mpf_pi(prec) and mpf_e(prec), whose
# default rounding (round_fast) is down: the constant at prec + 20 bits cut
# to prec bits. NULL for any other expression.
.mp_sympy_const <- function(e, prec) {
  if (e$type != "Constant" || !(e$name %in% c("pi", "E"))) return(NULL)
  v <- if (e$name == "pi") .mp_pi(prec + 20) else .mp_exp(.mpf(1, 0), prec + 20)
  s <- .bitlen(v$m) - prec
  if (s > 0) v <- .mp_round(gmp::divq.bigz(v$m, .pow2(s)), v$e + s, prec)
  r <- .mpc(v)
  attr(r, "converged") <- TRUE
  r
}

# errors raised by the mpmath port propagate out of N(), as in SymPy
.mp_error_handler <- function(err) if (inherits(err, "mpmath_error")) stop(err) else NULL

# ---------------------------------------------------------------------------
# Special functions (evalf-mp-functions.R): Function._eval_evalf with the
# mpmath function named in MPMATH_TRANSLATIONS, or the class's own _eval_evalf

# _eval_evalf through a rewrite: the rewritten expression evaluated at p
.mp_rewrite_eval <- function(e, target, p, env) {
  r <- .mpeval(rewrite(e, target), p, env)
  if (is.null(r)) NULL else .mv_mpc(r$re, r$im)
}
# local_workprec(prec): the arguments at prec and a context of the pool
.mp_local <- function(f) {
  force(f)
  function(e, p, env) {
    args <- .mp_sympy_args(e$args, p, env)
    if (is.null(args)) return(NULL)
    .mp_with_ctx("pool", .mp_workprec(p, do.call(f, args)))
  }
}
.mp_int_order <- function(e) e$args[[1L]]$type == "Integer"

.mp_reg("besselj", .mp_generic(function(n, z) .sf_besselj(n, z)))
.mp_reg("bessely", .mp_generic(function(n, z) .sf_bessely(n, z)))
.mp_reg("besseli", .mp_generic(function(n, z) .sf_besseli(n, z)))
.mp_reg("besselk", .mp_generic(function(n, z) .sf_besselk(n, z)))
.mp_reg("hankel1", .mp_generic(function(n, z) .sf_hankel(n, z, 1)))
.mp_reg("hankel2", .mp_generic(function(n, z) .sf_hankel(n, z, -1)))
.mp_reg(c("jn", "hn1", "hn2"), function(e, p, env) if (.mp_int_order(e)) .mp_sph_eval(e, "besselj", p, env) else NULL)
.mp_reg("yn", function(e, p, env) if (.mp_int_order(e)) .mp_sph_eval(e, "bessely", p, env) else NULL)

# jn, yn, hn1, hn2 of integer order: SymPy's rewrite(...)._eval_evalf(prec)
# computes in Float arithmetic. Every number of the rewritten product becomes
# a Float at prec, the Floats are multiplied in argument order (Mul.flatten,
# the coefficient of a Float times I last) and the product is distributed
# over a complex Float sum (Mul._eval_evalf expands), each product rounded.
# A complex z leaves sqrt(base) unevaluated: as_real_imag then fails at
# prec + 4, and evalf retries _eval_evalf at prec and evaluates the Float
# terms times sqrt(base) accurately.
.mp_sph_eval <- function(e, target, p, env) {
  r <- rewrite(e, target)
  v <- .mp_sph_float(r, p, env)
  if (is.null(v)) return(.mp_rewrite_eval(e, target, p, env))
  if (is.null(v$pw)) return(.mv_mpc(v$re, v$im))
  P <- p - 4L
  v <- .mp_sph_float(r, P, env)
  w <- .mpc_mul(.mpc(v$re, v$im), .mpc_sqrt(v$pw, P + 100L), P + 100L)
  .mv_mpc(.mp_round(w$re$m, w$re$e, P), .mp_round(w$im$m, w$im$e, P))
}

# A value in SymPy's Float arithmetic at precision P: real and imaginary
# Floats and the complex base of an unevaluated sqrt factor (or NULL); NULL
# for a shape this emulation does not cover
.mp_sph_float <- function(e, P, env) {
  val <- function(re, im = .mp_zero, pw = NULL) list(re = re, im = im, pw = pw)
  t <- e$type
  if (t %in% c("Integer", "Rational")) return(val(.mp_from_q(e$value, P)))
  if (t == "Constant" && e$name == "pi") return(val(.mp_pi(P)))
  if (t == "Constant" && e$name == "I") return(val(.mp_zero, .mpf(1, 0)))
  if (t == "Function") {
    f <- .MP_SYMPY[[e$name]]
    v <- if (is.null(f)) NULL else f(e, P, env)
    if (is.null(v)) return(NULL)
    v <- .mv(v)
    if (.mv_sp(v)) return(NULL)
    if (!.mv_isnum(v)) v <- .mv_convert(v)
    im <- if (v$k == "mpc") .mp_round(v$im$m, v$im$e, P) else .mp_zero
    return(val(.mp_round(v$re$m, v$re$e, P), im))
  }
  if (t == "Pow") {
    ex <- e$args[[2L]]
    if (!(ex$type == "Rational" && ex$value == gmp::as.bigq(1, 2))) return(NULL)
    b <- .mp_sph_float(e$args[[1L]], P, env)
    if (is.null(b) || !is.null(b$pw)) return(NULL)
    if (b$im$m == 0 && b$re$m > 0) return(val(.mp_sqrt(b$re, P)))
    if (b$im$m != 0) return(val(.mpf(1, 0), pw = .mpc(b$re, b$im)))
    return(NULL)
  }
  if (t == "Add") {
    re <- NULL; im <- NULL
    for (a in e$args) {
      v <- .mp_sph_float(a, P, env)
      if (is.null(v) || !is.null(v$pw)) return(NULL)
      if (v$re$m != 0) re <- if (is.null(re)) v$re else .mp_add(re, v$re, P)
      if (v$im$m != 0) im <- if (is.null(im)) v$im else .mp_add(im, v$im, P)
    }
    return(val(if (is.null(re)) .mp_zero else re, if (is.null(im)) .mp_zero else im))
  }
  if (t == "Mul") {
    coefs <- list(); later <- list(); ipow <- 0L; pair <- NULL; pw <- NULL
    for (a in e$args) {
      v <- .mp_sph_float(a, P, env)
      if (is.null(v)) return(NULL)
      if (!is.null(v$pw)) {
        if (!is.null(pw)) return(NULL)
        pw <- v$pw
      } else if (v$im$m == 0) {
        coefs[[length(coefs) + 1L]] <- v$re
      } else if (v$re$m == 0) {
        # Float*I: Mul.flatten appends its factors to the end of the sequence
        later[[length(later) + 1L]] <- v$im; ipow <- ipow + 1L
      } else {
        if (!is.null(pair)) return(NULL)
        pair <- v
      }
    }
    coef <- NULL
    for (c in c(coefs, later)) coef <- if (is.null(coef)) c else .mp_mul(coef, c, P)
    if (is.null(pair)) {
      re <- if (is.null(coef)) .mpf(1, 0) else coef; im <- .mp_zero
    } else if (is.null(coef)) {
      re <- pair$re; im <- pair$im
    } else {
      re <- .mp_mul(coef, pair$re, P); im <- .mp_mul(coef, pair$im, P)
    }
    for (k in seq_len(ipow %% 4L)) { tmp <- re; re <- .mp_neg(im); im <- tmp }
    return(val(re, im, pw))
  }
  NULL
}
.mp_reg("airyai", .mp_generic(function(z) .sf_airy(z, "ai"), 1L))
.mp_reg("airybi", .mp_generic(function(z) .sf_airy(z, "bi"), 1L))
.mp_reg("airyaiprime", .mp_local(function(z) .sf_airy(z, "ai", 1)))
.mp_reg("airybiprime", .mp_local(function(z) .sf_airy(z, "bi", 1)))
.mp_reg("erf", .mp_generic(.sf_erf, 1L))
.mp_reg("erfc", .mp_generic(.sf_erfc, 1L))
.mp_reg("erfi", .mp_generic(.sf_erfi, 1L))
.mp_reg("expint", .mp_generic(.sf_expint, 2L))
.mp_reg("Ei", .mp_generic(function(z) .sf_ei(z), 1L))
.mp_reg("Si", .mp_generic(.sf_si, 1L))
.mp_reg("Ci", .mp_generic(.sf_ci, 1L))
.mp_reg("Shi", .mp_generic(.sf_shi, 1L))
.mp_reg("Chi", .mp_generic(.sf_chi, 1L))
.mp_reg("li", .mp_generic(.sf_li, 1L))
# Li._eval_evalf: rewrite(li).evalf(n = prec_to_dps(prec))
.mp_reg("Li", function(e, p, env) {
  r <- .mpeval_adaptive(rewrite(e, "li"), .dps_to_prec(.prec_to_dps(p)), env)
  if (is.null(r)) NULL else .mv_mpc(r$re, r$im)
})
.mp_reg("fresnels", .mp_generic(function(z) .sf_fresnel(z, "s"), 1L))
.mp_reg("fresnelc", .mp_generic(function(z) .sf_fresnel(z, "c"), 1L))
.mp_reg("uppergamma", .mp_local(function(a, z) .sf_gammainc(a, z, .mv_inf())))
.mp_reg("lowergamma", .mp_local(function(a, z) .sf_gammainc(a, 0, z)))
.mp_reg("legendre", .mp_generic(.sf_legendre, 2L))
.mp_reg("assoc_legendre", .mp_generic(.sf_legenp, 3L))
.mp_reg("chebyshevt", .mp_generic(.sf_chebyt, 2L))
.mp_reg("chebyshevu", .mp_generic(.sf_chebyu, 2L))
.mp_reg("hermite", .mp_generic(.sf_hermite, 2L))
.mp_reg("gegenbauer", .mp_generic(.sf_gegenbauer, 3L))
.mp_reg("jacobi", .mp_generic(.sf_jacobi, 4L))
.mp_reg("elliptic_k", .mp_generic(.sf_ellipk, 1L))
.mp_reg("elliptic_f", .mp_generic(.sf_ellipf, 2L))
.mp_reg("elliptic_e", .mp_generic(.sf_ellipe, 1:2))
.mp_reg("elliptic_pi", .mp_generic(.sf_ellippi, 2:3))
.mp_reg("LambertW", .mp_generic(function(z, k = 0) .sf_lambertw(z, k), 1:2))
