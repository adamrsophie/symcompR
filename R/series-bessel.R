# Series methods of the Bessel functions and the Gruntz helpers _besseli
# and _besselk, after sympy/functions/special/bessel.py.

.sb_mexpand <- function(e) expand_mul(expand_multinomial(e))
.sb_trunc <- function(e, o) .so_removeO(.add(list(.sb_mexpand(e), o)))
.sb_lt_exp <- function(z, x) tryCatch(.so_leadterm(z, x)[[2L]], ValueError = function(err) NULL, NotImplementedError = function(err) NULL)
.sb_ratio <- function(a, b) .div(a, .int(b))
.sb_nn <- function(n) as.integer(num_value_dbl(as_sym(n)))

# besselj and besseli at a zero of the argument
.sb_ji_ns <- function(nm) { force(nm); function(e, x, n, logx, cdir) {
  nu <- e$args[[1L]]; z <- e$args[[2L]]
  ex <- .sb_lt_exp(z, x)
  if (is.null(ex)) return(e)
  if (isTRUE(ask_fact(ex, "positive"))) {
    newn <- .sb_nn(.ns_ceiling(.div(as_sym(n), ex)))
    o <- .so_new(.pow(x, n), list(x))
    r <- .so_removeO(.ns_eval(.div(z, S_Two), x, n, logx, cdir))
    if (is_zero_num(r)) return(o)
    t <- .sb_trunc(.sq(r), o)
    term <- .div(.pow(r, nu), fcall("gamma", list(.add(list(nu, S_One)))))
    s <- list(term)
    for (k in seq_len(max((newn + 1L) %/% 2L - 1L, 0L))) {
      fac <- .div(t, .mul(list(.int(k), .add(list(nu, .int(k))))))
      term <- .mul(list(term, if (nm == "besselj") .neg(fac) else fac))
      term <- .sb_trunc(term, o)
      s[[length(s) + 1L]] <- term
    }
    return(.add(c(s, list(o))))
  }
  .ns_function(e, x, n, logx, 0)
} }
.ns_hooks$besselj <- .sb_ji_ns("besselj")
.ns_hooks$besseli <- .sb_ji_ns("besseli")

.sb_lt_arg <- function(z, x) tryCatch(as_leading_term(z, x), NotImplementedError = function(err) NULL)
.so_alt_hooks$besselj <- function(e, x, logx, cdir) {
  nu <- e$args[[1L]]; z <- e$args[[2L]]
  arg <- .sb_lt_arg(z, x)
  if (is.null(arg)) return(e)
  ce <- .so_as_coeff_exponent(arg, x)
  if (isTRUE(ask_fact(ce[[2L]], "positive"))) return(.div(.pow(arg, nu), .mul(list(.pow(S_Two, nu), fcall("gamma", list(.add(list(nu, S_One))))))))
  if (isTRUE(ask_fact(ce[[2L]], "negative"))) {
    cd <- if (is_zero_num(as_sym(cdir))) S_One else as_sym(cdir)
    if (!isTRUE(ask_fact(.mul(list(ce[[1L]], .pow(cd, ce[[2L]]))), "negative")))
      return(.div(.mul(list(.sqrt(S_Two), cos(.sub(z, .div(.mul(list(S_Pi, .add(list(.mul(list(S_Two, nu)), S_One)))), .int(4L)))))), .sqrt(.mul(list(S_Pi, z)))))
    return(e)
  }
  .so_alt_function(e, x, logx, cdir)
}
.so_alt_hooks$besseli <- function(e, x, logx, cdir) {
  nu <- e$args[[1L]]; z <- e$args[[2L]]
  arg <- .sb_lt_arg(z, x)
  if (is.null(arg)) return(e)
  ce <- .so_as_coeff_exponent(arg, x)
  if (isTRUE(ask_fact(ce[[2L]], "positive"))) return(.div(.pow(arg, nu), .mul(list(.pow(S_Two, nu), fcall("gamma", list(.add(list(nu, S_One))))))))
  if (isTRUE(ask_fact(ce[[2L]], "negative"))) {
    cd <- if (is_zero_num(as_sym(cdir))) S_One else as_sym(cdir)
    if (!isTRUE(ask_fact(.mul(list(ce[[1L]], .pow(cd, ce[[2L]]))), "negative")))
      return(.div(exp(z), .sqrt(.mul(list(S_Two, S_Pi, z)))))
    return(e)
  }
  .so_alt_function(e, x, logx, cdir)
}
.so_alt_hooks$bessely <- function(e, x, logx, cdir) {
  nu <- e$args[[1L]]; z <- e$args[[2L]]
  arg <- .sb_lt_arg(z, x)
  if (is.null(arg)) return(e)
  ce <- .so_as_coeff_exponent(arg, x)
  if (isTRUE(ask_fact(ce[[2L]], "positive"))) {
    t1 <- .mul(list(.div(S_Two, S_Pi), log(.div(z, S_Two)), fcall("besselj", list(nu, z))))
    t2 <- if (isTRUE(ask_fact(nu, "positive"))) .div(.neg(.mul(list(.pow(.div(z, S_Two), .neg(nu)), fcall("factorial", list(.sub(nu, S_One)))))), S_Pi) else S_Zero
    t3 <- .mul(list(.neg(.div(.pow(.div(z, S_Two), nu), .mul(list(S_Pi, fcall("factorial", list(nu)))))), .sub(fcall("polygamma", list(S_Zero, .add(list(nu, S_One)))), S_EulerGamma)))
    return(as_leading_term(.add(list(t1, t2, t3)), x, logx))
  }
  if (isTRUE(ask_fact(ce[[2L]], "negative"))) {
    cd <- if (is_zero_num(as_sym(cdir))) S_One else as_sym(cdir)
    if (!isTRUE(ask_fact(.mul(list(ce[[1L]], .pow(cd, ce[[2L]]))), "negative"))) {
      ph <- .add(list(.div(.mul(list(S_Pi, nu)), S_Two), .neg(z), .div(S_Pi, .int(4L))))
      return(.div(.mul(list(.sqrt(S_Two), .add(list(.neg(sin(ph)), .div(.mul(list(.int(3L), cos(ph))), .mul(list(.int(8L), z))))), .sqrt(.inv(z)))), .sqrt(S_Pi)))
    }
    return(e)
  }
  .so_alt_function(e, x, logx, cdir)
}
.ns_hooks$bessely <- function(e, x, n, logx, cdir) {
  nu <- e$args[[1L]]; z <- e$args[[2L]]
  ex <- .sb_lt_exp(z, x)
  if (is.null(ex)) return(e)
  if (isTRUE(ask_fact(ex, "positive")) && isTRUE(ask_fact(nu, "integer"))) {
    newn <- .sb_nn(.ns_ceiling(.div(as_sym(n), ex)))
    bn <- fcall("besselj", list(nu, z))
    a <- .ns_eval(.mul(list(.div(S_Two, S_Pi), log(.div(z, S_Two)), bn)), x, n, logx, cdir)
    b <- list(); cc <- list()
    o <- .so_new(.pow(x, n), list(x))
    r <- .so_removeO(.ns_eval(.div(z, S_Two), x, n, logx, cdir))
    if (is_zero_num(r)) return(o)
    t <- .sb_trunc(.sq(r), o)
    nui <- as.integer(num_value_dbl(nu))
    if (nui > 0L) {
      term <- .div(.mul(list(.pow(r, .neg(nu)), fcall("factorial", list(.sub(nu, S_One))))), S_Pi)
      b[[1L]] <- term
      for (k in seq_len(nui - 1L)) {
        den <- (nui - k) * k
        term <- .mul(list(term, if (den == 0L) .div(t, .int(k)) else .div(t, .int(den))))
        term <- .sb_trunc(term, o)
        b[[length(b) + 1L]] <- term
      }
    }
    p <- .div(.pow(r, nu), .mul(list(S_Pi, fcall("factorial", list(nu)))))
    cc[[1L]] <- .mul(list(p, .sub(fcall("polygamma", list(S_Zero, .add(list(nu, S_One)))), S_EulerGamma)))
    for (k in seq_len(max((newn + 1L) %/% 2L - 1L, 0L))) {
      p <- .mul(list(p, .neg(.div(t, .mul(list(.int(k), .add(list(.int(k), nu))))))))
      p <- .sb_trunc(p, o)
      cc[[length(cc) + 1L]] <- .mul(list(p, .add(list(fcall("polygamma", list(S_Zero, .add(list(.int(k), nu, S_One)))), fcall("polygamma", list(S_Zero, .int(k + 1L)))))))
    }
    return(.sub(.sub(a, .add(b)), .add(cc)))
  }
  .ns_function(e, x, n, logx, 0)
}
.so_alt_hooks$besselk <- function(e, x, logx, cdir) {
  nu <- e$args[[1L]]; z <- e$args[[2L]]
  arg <- .sb_lt_arg(z, x)
  if (is.null(arg)) return(e)
  ex <- .so_as_coeff_exponent(arg, x)[[2L]]
  if (isTRUE(ask_fact(ex, "positive"))) {
    if (isTRUE(ask_fact(nu, "zero"))) term <- .add(list(.neg(log(z)), .neg(S_EulerGamma), log(S_Two)))
    else if (isTRUE(ask_fact(nu, "nonzero"))) term <- .div(.mul(list(fcall("gamma", list(Abs(nu))), .pow(.div(z, S_Two), .neg(Abs(nu))))), S_Two)
    else .notimpl_error(paste0("Cannot proceed without knowing if ", .str(nu), " is zero or not."))
    return(as_leading_term(term, x, logx))
  }
  if (isTRUE(ask_fact(ex, "negative"))) return(.div(.mul(list(.sqrt(S_Pi), exp(.neg(arg)))), .sqrt(.mul(list(S_Two, arg)))))
  fcall("besselk", list(nu, arg))
}
.ns_hooks$besselk <- function(e, x, n, logx, cdir) {
  nu <- e$args[[1L]]; z <- e$args[[2L]]
  ex <- .sb_lt_exp(z, x)
  if (is.null(ex)) return(e)
  if (isTRUE(ask_fact(ex, "positive"))) {
    r <- .so_removeO(.ns_eval(.div(z, S_Two), x, n, logx, cdir))
    if (is_zero_num(r)) return(.so_new(.add(list(.pow(z, .neg(nu)), .pow(z, nu))), list(x)))
    o <- .so_new(.pow(x, n), list(x))
    if (isTRUE(ask_fact(nu, "integer"))) {
      newn <- .sb_nn(.ns_ceiling(.div(as_sym(n), ex)))
      bn <- fcall("besseli", list(nu, z))
      a <- .ns_eval(.mul(list(.pow(S_NegativeOne, .sub(nu, S_One)), log(.div(z, S_Two)), bn)), x, n, logx, cdir)
      b <- list(); cc <- list()
      t <- .sb_mexpand(.sq(r))
      nui <- as.integer(num_value_dbl(nu))
      if (nui > 0L) {
        term <- .div(.mul(list(.pow(r, .neg(nu)), fcall("factorial", list(.sub(nu, S_One))))), S_Two)
        b[[1L]] <- term
        for (k in seq_len(nui - 1L)) {
          term <- .mul(list(term, .div(t, .int((k - nui) * k))))
          term <- .sb_trunc(term, o)
          b[[length(b) + 1L]] <- term
        }
      }
      p <- .div(.mul(list(.pow(r, nu), .pow(S_NegativeOne, nu))), .mul(list(S_Two, fcall("factorial", list(nu)))))
      cc[[1L]] <- .mul(list(p, .sub(fcall("polygamma", list(S_Zero, .add(list(nu, S_One)))), S_EulerGamma)))
      for (k in seq_len(max((newn + 1L) %/% 2L - 1L, 0L))) {
        p <- .mul(list(p, .div(t, .mul(list(.int(k), .add(list(.int(k), nu)))))))
        p <- .sb_trunc(p, o)
        cc[[length(cc) + 1L]] <- .mul(list(p, .add(list(fcall("polygamma", list(S_Zero, .add(list(.int(k), nu, S_One)))), fcall("polygamma", list(S_Zero, .int(k + 1L)))))))
      }
      return(.add(c(list(a), b, cc, list(o))))
    }
    if (isTRUE(ask_fact(nu, "noninteger"))) {
      na <- .sb_nn(.ns_ceiling(.div(.add(list(as_sym(n), nu)), ex)))
      nb <- .sb_nn(.ns_ceiling(.div(.sub(as_sym(n), nu), ex)))
      a <- lapply(seq_len(max((na + 1L) %/% 2L, 0L)) - 1L, function(k)
        .sb_mexpand(.div(.mul(list(fcall("gamma", list(nu)), .pow(r, .sub(.int(2L * k), nu)))), .mul(list(S_Two, fcall("RisingFactorial", list(.sub(S_One, nu), .int(k))), .se_fact(k))))))
      b <- lapply(seq_len(max((nb + 1L) %/% 2L, 0L)) - 1L, function(k)
        .sb_mexpand(.div(.mul(list(fcall("gamma", list(.neg(nu))), .pow(r, .add(list(.int(2L * k), nu))))), .mul(list(S_Two, fcall("RisingFactorial", list(.add(list(nu, S_One)), .int(k))), .se_fact(k))))))
      return(.add(c(a, b, list(o))))
    }
    .notimpl_error("besselk expansion is only implemented for real order")
  }
  .ns_function(e, x, n, logx, 0)
}
.sb_aseries <- function(nm, sgn, pref) { force(nm); force(sgn); force(pref); function(e, n, args0, x, logx) {
  pt <- args0[[2L]]
  if (!.sp_inf_pts(pt)) .pole_error(paste0("Asymptotic expansion of ", nm, " around [", paste(vapply(args0, .str, ""), collapse = ", "), "] is not implemented."))
  nu <- e$args[[1L]]; z <- e$args[[2L]]; nn <- .sb_nn(n)
  l <- lapply(seq_len(nn) - 1L, function(k)
    .div(.mul(list(fcall("RisingFactorial", list(.div(.sub(.mul(list(S_Two, nu)), S_One), S_Two), .int(k))),
                   fcall("RisingFactorial", list(.div(.add(list(.mul(list(S_Two, nu)), S_One)), S_Two), .int(k))))),
         .mul(list(.pow(.int(sgn * 2L), .int(k)), .pow(z, .div(.int(2L * k + 1L), S_Two)), .se_fact(k)))))
  o <- .so_new(.inv(.pow(z, .div(.int(2L * nn + 1L), S_Two))), list(x))
  switch(pref,
    besseli = .mul(list(exp(z), .inv(.sqrt(.mul(list(S_Two, S_Pi)))), .add(c(l, list(o))))),
    besselk = .mul(list(exp(.neg(z)), .sqrt(.div(S_Pi, S_Two)), .add(c(l, list(o))))),
    .add(list(.mul(list(.sqrt(.div(S_Pi, S_Two)), .add(l))), o)))
} }
.ns_aseries_hooks$besseli <- .sb_aseries("besseli", 1L, "besseli")
.ns_aseries_hooks$besselk <- .sb_aseries("besselk", -1L, "besselk")
.ns_aseries_hooks$`_besseli` <- .sb_aseries("_besseli", 1L, "helper")
.ns_aseries_hooks$`_besselk` <- .sb_aseries("_besselk", -1L, "helper")

def_function("_besseli", nargs = 2L)
def_function("_besselk", nargs = 2L)
.sb_helper_ns <- function(nm) { force(nm); function(e, x, n, logx, cdir) {
  # SymPy tests the order argument here
  if (isTRUE(ask_fact(.sp_lim0(e$args[[1L]], x), "zero"))) {
    f <- if (nm == "_besseli") .mul(list(exp(.neg(e$args[[2L]])), fcall("besseli", e$args))) else .mul(list(exp(e$args[[2L]]), fcall("besselk", e$args)))
    return(.ns_eval(f, x, n, logx))
  }
  .ns_function(e, x, n, logx, 0)
} }
.ns_hooks$`_besseli` <- .sb_helper_ns("_besseli")
.ns_hooks$`_besselk` <- .sb_helper_ns("_besselk")

# airyai and airybi taylor terms
.ns_taylor_hooks$airyai <- function(k, a, prev, e) {
  if (k < 0L) return(S_Zero)
  c3 <- .pow(.int(3L), .div(S_One, .int(3L)))
  .mul(list(.inv(.mul(list(.pow(.int(3L), .div(S_Two, .int(3L))), S_Pi))), fcall("gamma", list(.div(.int(k + 1L), .int(3L)))),
            sin(.mul(list(.div(S_Two, .int(3L)), S_Pi, .int(k + 1L)))), .inv(.se_fact(k)), .pow(.mul(list(c3, a)), .int(k))))
}
.ns_taylor_hooks$airybi <- function(k, a, prev, e) {
  if (k < 0L) return(S_Zero)
  c3 <- .pow(.int(3L), .div(S_One, .int(3L)))
  .mul(list(.inv(.mul(list(.pow(.int(3L), .div(S_One, .int(6L))), S_Pi))), fcall("gamma", list(.div(.int(k + 1L), .int(3L)))),
            Abs(sin(.mul(list(.div(S_Two, .int(3L)), S_Pi, .int(k + 1L))))), .inv(.se_fact(k)), .pow(.mul(list(c3, a)), .int(k))))
}
