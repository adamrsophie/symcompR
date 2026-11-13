# arg.eval, after sympy/functions/elementary/complexes.py (class arg).

.coremm_undefs <- function(e) {
  out <- list()
  walk <- function(x) {
    if (.is_undefined_application(x)) out[[x$key]] <<- x
    for (a in x$args) walk(a)
  }
  walk(e)
  unname(out)
}

.coremm_arg_eval <- function(arg) {
  is_arg <- function(z) z$type == "Function" && identical(z$name, "arg")
  a <- arg
  broke <- FALSE
  for (i in 0:2) {
    if (is_arg(a)) {
      a <- a$args[[1L]]
    } else {
      if (i == 2L && isTRUE(ask_fact(a, "extended_real"))) return(S_NaN)
      broke <- TRUE
      break
    }
  }
  if (!broke) return(S_NaN)
  if (arg$type == "Function" && identical(arg$name, "exp_polar"))
    return(fcall("periodic_argument", list(arg, S_Infinity)))
  if (arg$type == "Function" && identical(arg$name, "exp")) {
    i_ <- fcall("im", list(arg$args[[1L]]))
    if (.is_comparable(i_)) {
      i_ <- Mod(i_, .mul(list(S_Two, S_Pi)))
      if (isTRUE(ask_fact(.sub(i_, S_Pi), "positive"))) i_ <- .sub(i_, .mul(list(S_Two, S_Pi)))
      return(i_)
    }
  }
  if (length(arg$args)) {
    cm <- .as_coeff_mul1(factor_terms(arg))
    c <- cm[[1L]]; arg_ <- cm[[2L]]
    if (arg_$type == "Mul") {
      arg_ <- .mul(lapply(.sympy_args(arg_), function(f) {
        s <- fcall("sign", list(f))
        if (is_one_num(s) || is_negone_num(s)) s else f
      }))
    }
    arg_ <- .mul(list(fcall("sign", list(c)), arg_))
  } else {
    arg_ <- arg
  }
  if (any(vapply(.coremm_undefs(arg_), function(u) is.na(ask_fact(u, "extended_positive")), TRUE)))
    return(NULL)
  ri <- .as_real_imag(arg_)
  rv <- fcall("atan2", list(ri[[2L]], ri[[1L]]))
  if (.is_number(rv)) return(rv)
  if (!identical(arg_$key, arg$key)) return(fnode("arg", list(arg_)))
  NULL
}

# atan2.eval, after sympy/functions/elementary/trigonometric.py (class atan2).
.coremm_atan2_eval <- function(y, x) {
  A <- function(e, f) ask_fact(e, f)
  if (is_const(x, "-oo")) {
    if (isTRUE(A(y, "zero"))) return(S_Pi)
    return(.sub(.mul(list(S_Two, S_Pi, fcall("Heaviside", list(fcall("re", list(y)))))), S_Pi))
  } else if (is_const(x, "oo")) {
    return(S_Zero)
  } else if (isTRUE(A(x, "imaginary")) && isTRUE(A(y, "imaginary")) && .is_number(x) && .is_number(y)) {
    x <- fcall("im", list(x))
    y <- fcall("im", list(y))
  }
  if (isTRUE(A(x, "extended_real")) && isTRUE(A(y, "extended_real"))) {
    if (isTRUE(A(x, "positive"))) {
      return(fcall("atan", list(.div(y, x))))
    } else if (isTRUE(A(x, "negative"))) {
      if (isTRUE(A(y, "negative"))) return(.sub(fcall("atan", list(.div(y, x))), S_Pi))
      else if (isTRUE(A(y, "nonnegative"))) return(.add(list(fcall("atan", list(.div(y, x))), S_Pi)))
    } else if (isTRUE(A(x, "zero"))) {
      if (isTRUE(A(y, "positive"))) return(.mul(list(S_Half, S_Pi)))
      else if (isTRUE(A(y, "negative"))) return(.mul(list(.rat(-1, 2), S_Pi)))
      else if (isTRUE(A(y, "zero"))) return(S_NaN)
    }
  }
  if (isTRUE(A(y, "zero"))) {
    if (isTRUE(A(x, "extended_nonzero")))
      return(.mul(list(S_Pi, .sub(S_One, fcall("Heaviside", list(x))))))
    if (.is_number(x))
      return(Piecewise(list(S_Pi, Lt(fcall("re", list(x)), S_Zero)), list(S_Zero, Ne(x, S_Zero)), list(S_NaN, S_true)))
  }
  if (.is_number(x) && .is_number(y))
    return(.mul(list(.neg(S_I), fcall("log", list(.div(.add(list(x, .mul(list(S_I, y)))),
                                                        .sqrt(.add(list(.sq(x), .sq(y))))))))))
  NULL
}

# The infinite branch of is_eq (sympy/core/relational.py): NULL when neither
# side is known to be infinite, else TRUE, FALSE or NA (undecided).
.coremm_is_eq_inf <- function(lhs, rhs, assumptions = NULL) {
  W <- function(e, f) .awf(e, f, assumptions)
  linf <- W(lhs, "infinite"); rinf <- W(rhs, "infinite")
  if (!(isTRUE(linf) || isTRUE(rinf))) return(NULL)
  if (isTRUE(.fuzzy_xor(list(linf, rinf)))) return(FALSE)
  lr <- W(lhs, "extended_real"); rr <- W(rhs, "extended_real")
  if (isTRUE(.fuzzy_xor(list(lr, rr)))) return(FALSE)
  if (isTRUE(.fuzzy_and_vals(list(lr, rr)))) {
    if (isTRUE(.fuzzy_and_vals(list(linf, rinf))))
      return(.fuzzy_xor(list(W(lhs, "extended_positive"), .fuzzy_not(W(rhs, "extended_positive")))))
    if ((isTRUE(W(lhs, "extended_positive")) && isTRUE(W(rhs, "extended_nonpositive"))) ||
        (isTRUE(W(lhs, "extended_negative")) && isTRUE(W(rhs, "extended_nonnegative"))) ||
        (isTRUE(W(lhs, "extended_nonpositive")) && isTRUE(W(rhs, "extended_positive"))) ||
        (isTRUE(W(lhs, "extended_nonnegative")) && isTRUE(W(rhs, "extended_negative")))) return(FALSE)
    return(NA)
  }
  split_ri <- function(e) {
    out <- list(real = list(), imag = list(), none = list())
    for (t in .add_terms(e)) {
      g <- if (isTRUE(W(t, "extended_real"))) "real"
           else if (isTRUE(W(.mul(list(S_I, t)), "extended_real"))) "imag" else "none"
      out[[g]][[length(out[[g]]) + 1L]] <- t
    }
    out
  }
  lri <- split_ri(lhs)
  if (!length(lri$none)) {
    rri <- split_ri(rhs)
    if (!length(rri$none)) {
      eq_real <- .is_eq_a(.add(lri$real), .add(rri$real), assumptions)
      eq_imag <- .is_eq_a(.mul(list(S_I, .add(lri$imag))), .mul(list(S_I, .add(rri$imag))), assumptions)
      return(.fuzzy_and_vals(list(eq_real, eq_imag)))
    }
  }
  al <- fcall("arg", list(lhs)); ar <- fcall("arg", list(rhs))
  if (!(.is_nan_sym(al) && .is_nan_sym(ar))) return(.is_eq_a(al, ar, assumptions))
  NA
}

# Poly(f, n).all_coeffs() for a univariate f: the coefficients are converted
# to the domain construct_domain() finds, so over RR or RR[...] every numeric
# coefficient becomes a Float (1.0*pi, 3.00000000000000).
# ponytail: RR and RR[...] grounds only; CC and fraction fields keep the
# expression coefficients.
.coremm_all_coeffs <- function(f, n) {
  cl <- .coeff_list(f, n)
  K <- construct_domain(cl)[[1L]]
  g <- attr(K, "dom")
  ground <- if (identical(g$kind, "PolyRing")) attr(g$dom, "dom")$name else g$name
  if (!identical(ground, "RR") || !(g$kind %in% c("RR", "PolyRing"))) return(cl)
  lapply(cl, function(c) {
    if (is_zero_num(c)) return(S_Zero)
    .add(lapply(.aargs(expand(c)), function(t) {
      sp <- as_coeff_mul_split(t)
      cf <- if (sp[[1L]]$type == "Float") sp[[1L]] else N(sp[[1L]], 15)
      .coeff_times(cf, sp[[2L]])
    }))
  })
}

# Pow._eval_is_imaginary, after sympy/core/power.py (exp(x) is a function
# here, so the base E branch is not needed).
.coremm_pow_imaginary <- function(b, x, fb, fx) {
  g <- function(f, nm) { v <- f[nm]; if (length(v) != 1L) NA else unname(v) }
  if (identical(g(fb, "commutative"), FALSE)) return(FALSE)
  if (isTRUE(g(fb, "imaginary")) && isTRUE(g(fx, "integer"))) return(g(fx, "odd"))
  if (isTRUE(g(fx, "imaginary"))) {
    imlog <- ask_fact(fcall("log", list(b)), "imaginary")
    if (!is.na(imlog)) return(FALSE)
  }
  if (isTRUE(g(fb, "extended_real")) && isTRUE(g(fx, "extended_real"))) {
    if (isTRUE(g(fb, "positive"))) return(FALSE)
    rat <- g(fx, "rational")
    if (!isTRUE(rat)) return(rat)
    if (isTRUE(g(fx, "integer"))) return(FALSE)
    half <- ask_fact(.mul(list(S_Two, x)), "integer")
    if (isTRUE(half)) return(g(fb, "negative"))
    return(half)
  }
  if (identical(g(fb, "extended_real"), FALSE)) {
    i <- .div(.mul(list(fcall("arg", list(b)), x)), S_Pi)
    isodd <- ask_fact(.mul(list(S_Two, i)), "odd")
    if (!is.na(isodd)) return(isodd)
  }
  NA
}

# terms_gcd over RR (sympy/polys/polytools.py): with Float coefficients the
# content is 1, only the common monomial comes out, and every coefficient of
# the cofactor is a Float (y*(1.0*C1 - 51.05*a)). Other inexact domains keep
# the expression.
.coremm_terms_gcd_rr <- function(p, gens, orig) {
  real <- vapply(p$terms, function(t) is_exact_number(t$coeff) || t$coeff$type == "Float", TRUE)
  if (!all(real)) return(orig)
  J <- Reduce(pmin, lapply(p$terms, `[[`, "monom"))
  if (all(J == 0L)) return(orig)
  rest <- .add(lapply(p$terms, function(t) {
    m <- t$monom - J
    cf <- if (t$coeff$type == "Float") t$coeff else Float(as.double(gmp::as.bigq(t$coeff$value)))
    .mul(c(list(cf), Map(function(g, k) .pow(g, .int(k)), gens[m > 0L], m[m > 0L])))
  }))
  term <- .mul(Map(function(g, k) .pow(g, .int(k)), gens[J > 0L], J[J > 0L]))
  .keep_coeff(S_One, .mul(list(term, rest)), clear = FALSE)
}

# terms_gcd over ZZ_I: I is part of the coefficients, not a generator, and
# the content is the Gaussian gcd of the coefficients (dup_content), so
# -I*C1 + I*x^4 keeps its form.
# ponytail: QQ_I (Gaussian rationals) keeps the expression; clear_denoms
# over QQ_I when a case needs it.
.coremm_terms_gcd_zzi <- function(f, gens, clear, orig) {
  gx <- Filter(function(g) !is_const(g, "I"), gens)
  if (!length(gx)) return(orig)
  d <- Dummy("i")
  p <- .sp_as_poly_noexpand(xreplace(f, list(list(S_I, d))), c(gx, list(d)))
  if (is.null(p) || !length(p$terms)) return(orig)
  n <- length(gx)
  cfs <- list(); mons <- list()
  for (t in p$terms) {
    if (!is_exact_number(t$coeff)) return(orig)
    k <- paste(t$monom[seq_len(n)], collapse = ",")
    v <- .mul(list(t$coeff, .pow(S_I, .int(t$monom[n + 1L]))))
    if (is.null(cfs[[k]])) { cfs[[k]] <- v; mons[[k]] <- t$monom[seq_len(n)] }
    else cfs[[k]] <- .add(list(cfs[[k]], v))
  }
  parts <- lapply(cfs, .gi_parts)
  if (!all(vapply(parts, function(v) all(gmp::denominator(v) == 1), NA))) return(orig)
  cont <- Reduce(.gi_gcd, parts, gmp::as.bigq(c(0, 0)))
  J <- Reduce(pmin, mons)
  rest <- .add(lapply(names(cfs), function(k) {
    m <- mons[[k]] - J
    q <- .gc_div(parts[[k]], cont)
    .mul(c(list(.gi_make(q[1], q[2])), Map(function(g, e) .pow(g, .int(e)), gx[m > 0L], m[m > 0L])))
  }))
  term <- .mul(Map(function(g, e) .pow(g, .int(e)), gx[J > 0L], J[J > 0L]))
  coeff <- .gi_make(cont[1], cont[2])
  if (is_one_num(coeff)) {
    coeff <- S_One
    if (is_one_num(term)) return(orig)
  }
  if (clear) return(.keep_coeff(coeff, .mul(list(term, rest))))
  kc <- as_coeff_mul_split(.keep_coeff(coeff, rest, clear = FALSE))
  .keep_coeff(kc[[1L]], .mul(list(term, kc[[2L]])), clear = FALSE)
}

# _find_reasonable_pivot, after sympy/matrices/determinant.py: a definite
# nonzero entry first, then the undecided entries after simplify(), and an
# entry still undecided is assumed nonzero. Returns list(index, value,
# newly) with index NA when the column has no pivot; newly holds the entries
# simplify() decided, as list(index, value).
.coremm_find_pivot <- function(col) {
  newly <- list()
  z <- vapply(col, function(x) {
    if (is_zero_num(x)) return(TRUE)
    if (.is_number(x) && is_number_type(x)) return(FALSE)
    ask_fact(x, "zero")
  }, NA)
  hit <- which(z %in% FALSE)
  if (length(hit)) return(list(hit[1L], col[[hit[1L]]], newly))
  if (all(z %in% TRUE)) return(list(NA_integer_, NULL, newly))
  for (i in which(is.na(z))) {
    s <- simplify(col[[i]])
    zi <- if (is_zero_num(s)) TRUE else ask_fact(s, "zero")
    if (!is.na(zi)) newly[[length(newly) + 1L]] <- list(i, if (zi) S_Zero else s)
    if (isFALSE(zi)) return(list(i, s, newly))
    z[i] <- zi
  }
  if (all(z %in% TRUE)) return(list(NA_integer_, NULL, newly))
  i <- which(is.na(z))[1L]
  list(i, col[[i]], newly)
}

# Pow._eval_conjugate, after sympy/core/power.py: NULL leaves conjugate()
# unevaluated.
.coremm_pow_conjugate <- function(a) {
  b <- a$args[[1L]]; x <- a$args[[2L]]
  i <- ask_fact(x, "integer"); p <- ask_fact(b, "positive")
  if (isTRUE(i)) return(.pow(fcall("conjugate", list(b)), x))
  if (isTRUE(p)) return(.pow(b, fcall("conjugate", list(x))))
  if (isFALSE(i) && isFALSE(p)) {
    ex <- expand_complex(a)
    if (!identical(ex$key, a$key)) return(fcall("conjugate", list(ex)))
  }
  if (isTRUE(ask_fact(a, "extended_real"))) return(a)
  NULL
}
