# Expr.as_real_imag and the evaluation of re() and im(), after
# sympy/functions/elementary/complexes.py and the as_real_imag methods of
# Add, Mul, Pow, exp, log, the trigonometric and hyperbolic functions.

.ari_zero <- function(e) isTRUE(ask_fact(e, "zero"))
.ari_struct_zero <- function(e) is_zero_num(e)
.ari_expand <- function(e) expand(e)

# Expr.as_real_imag(deep, ignore): list(re, im), or NULL when the result
# would be `ignore` itself
.as_real_imag <- function(e, deep = TRUE, ignore = NULL) {
  ign <- function(z) !is.null(ignore) && identical(z$key, ignore$key)
  default <- function() if (ign(e)) NULL else list(fcall("re", list(e)), fcall("im", list(e)))
  t <- e$type
  if (t == "Add") {
    parts <- lapply(.sympy_args(e), .as_real_imag, deep = deep)
    return(list(.add(lapply(parts, `[[`, 1L)), .add(lapply(parts, `[[`, 2L))))
  }
  if (t == "Mul") return(.ari_mul(e, ignore))
  if (t == "Pow") return(.ari_pow(e, deep, ignore))
  if (t == "Function") {
    a <- e$args[[1L]]
    real_arg <- function() isTRUE(ask_fact(a, "extended_real"))
    arg_ri <- function() if (deep) .as_real_imag(.ari_expand(a)) else .as_real_imag(a)
    nm <- e$name
    switch(nm,
      re = , im = return(list(e, S_Zero)),
      exp = {
        ri <- .as_real_imag(a)
        re_ <- ri[[1L]]; im_ <- ri[[2L]]
        if (deep) { re_ <- .ari_expand(re_); im_ <- .ari_expand(im_) }
        return(list(.mul(list(fcall("exp", list(re_)), fcall("cos", list(im_)))),
                    .mul(list(fcall("exp", list(re_)), fcall("sin", list(im_))))))
      },
      log = {
        sarg <- if (deep) .ari_expand(a) else a
        sabs <- fcall("Abs", list(sarg))
        if (identical(sabs$key, sarg$key)) return(list(e, S_Zero))
        return(list(fcall("log", list(sabs)), fcall("arg", list(sarg))))
      },
      sin = , cos = , tan = , cot = {
        ri <- if (real_arg()) list(if (deep) .ari_expand(a) else a, S_Zero) else arg_ri()
        re_ <- ri[[1L]]; im_ <- ri[[2L]]
        f <- function(n, z) fcall(n, list(z))
        if (nm == "sin") return(list(.mul(list(f("sin", re_), f("cosh", im_))), .mul(list(f("cos", re_), f("sinh", im_)))))
        if (nm == "cos") return(list(.mul(list(f("cos", re_), f("cosh", im_))), .neg(.mul(list(f("sin", re_), f("sinh", im_))))))
        if (is_zero_num(im_)) return(list(f(nm, re_), S_Zero))
        tr <- .mul(list(S_Two, re_)); ti <- .mul(list(S_Two, im_))
        if (nm == "tan") {
          den <- .add(list(f("cos", tr), f("cosh", ti)))
          return(list(.div(f("sin", tr), den), .div(f("sinh", ti), den)))
        }
        den <- .sub(f("cos", tr), f("cosh", ti))
        return(list(.neg(.div(f("sin", tr), den)), .div(f("sinh", ti), den)))
      },
      sec = return(.as_real_imag(.inv(fcall("cos", list(a))), deep)),
      csc = return(.as_real_imag(.inv(fcall("sin", list(a))), deep)),
      sinh = , cosh = , tanh = , coth = {
        if (real_arg()) return(list(if (deep) .ari_expand(e) else e, S_Zero))
        ri <- arg_ri(); re_ <- ri[[1L]]; im_ <- ri[[2L]]
        f <- function(n, z) fcall(n, list(z))
        if (nm == "sinh") return(list(.mul(list(f("sinh", re_), f("cos", im_))), .mul(list(f("cosh", re_), f("sin", im_)))))
        if (nm == "cosh") return(list(.mul(list(f("cosh", re_), f("cos", im_))), .mul(list(f("sinh", re_), f("sin", im_)))))
        if (nm == "tanh") {
          den <- .add(list(.sq(f("sinh", re_)), .sq(f("cos", im_))))
          return(list(.div(.mul(list(f("sinh", re_), f("cosh", re_))), den), .div(.mul(list(f("sin", im_), f("cos", im_))), den)))
        }
        den <- .add(list(.sq(f("sinh", re_)), .sq(f("sin", im_))))
        return(list(.div(.mul(list(f("sinh", re_), f("cosh", re_))), den), .neg(.div(.mul(list(f("sin", im_), f("cos", im_))), den))))
      },
      sech = return(.as_real_imag(.inv(fcall("cosh", list(a))), deep)),
      csch = return(.as_real_imag(.inv(fcall("sinh", list(a))), deep)),
      airyai = , airybi = , airyaiprime = , airybiprime = {
        zc <- fcall("conjugate", list(a))
        fz <- fcall(nm, list(a)); fzc <- fcall(nm, list(zc))
        return(list(.div(.add(list(fz, fzc)), S_Two), .div(.mul(list(S_I, .sub(fzc, fz))), S_Two)))
      })
    return(default())
  }
  default()
}

# Mul.as_real_imag
.ari_mul <- function(e, ignore) {
  other <- list(); coeffr <- list(); coeffi <- list(); addterms <- S_One
  for (a in .sympy_args(e)) {
    ri <- .as_real_imag(a)
    r <- ri[[1L]]; i <- ri[[2L]]
    if (.ari_zero(i)) coeffr[[length(coeffr) + 1L]] <- r
    else if (.ari_zero(r)) coeffi[[length(coeffi) + 1L]] <- .mul(list(i, S_I))
    else if (.is_comm(a)) {
      aconj <- if (length(other)) fcall("conjugate", list(a)) else NULL
      hit <- 0L
      if (!is.null(aconj)) for (k in seq_along(other)) if (identical(other[[k]]$key, aconj$key)) { hit <- k; break }
      if (hit) { coeffr[[length(coeffr) + 1L]] <- .sq(fcall("Abs", list(other[[hit]]))); other[[hit]] <- NULL }
      else if (a$type == "Add") addterms <- .mul(list(addterms, a))
      else other[[length(other) + 1L]] <- a
    } else other[[length(other) + 1L]] <- a
  }
  m <- .mul(other)
  if (!is.null(ignore) && identical(m$key, ignore$key)) return(NULL)
  if (length(coeffi) %% 2L) { imco <- fcall("im", list(coeffi[[1L]])); coeffi <- coeffi[-1L] }
  else imco <- S_Zero
  reco <- .mul(c(coeffr, coeffi))
  r <- .mul(list(reco, fcall("re", list(m)))); i <- .mul(list(reco, fcall("im", list(m))))
  if (is_one_num(addterms)) {
    if (is_one_num(m)) {
      if (.ari_zero(imco)) return(list(reco, S_Zero))
      return(list(S_Zero, .mul(list(reco, imco))))
    }
    if (is_zero_num(imco)) return(list(r, i))
    return(list(.mul(list(.neg(imco), i)), .mul(list(imco, r))))
  }
  ai <- .as_real_imag(expand_mul(addterms, deep = FALSE))
  addre <- ai[[1L]]; addim <- ai[[2L]]
  if (is_zero_num(imco))
    return(list(.sub(.mul(list(r, addre)), .mul(list(i, addim))), .add(list(.mul(list(i, addre)), .mul(list(r, addim))))))
  r2 <- .mul(list(.neg(imco), i)); i2 <- .mul(list(imco, r))
  list(.sub(.mul(list(r2, addre)), .mul(list(i2, addim))), .add(list(.mul(list(r2, addim)), .mul(list(i2, addre)))))
}

# Pow.as_real_imag
.ari_pow <- function(e, deep, ignore) {
  b <- e$args[[1L]]; ex <- e$args[[2L]]
  if (ex$type == "Integer") {
    n <- .int_val(ex)
    if (n < 0 && .ari_zero(b)) return(list(S_NaN, S_NaN))
    ri <- .as_real_imag(b, deep = deep)
    re_e <- ri[[1L]]; im_e <- ri[[2L]]
    if (is_zero_num(im_e)) return(list(e, S_Zero))
    if (n < 0) {
      mag <- .add(list(.sq(re_e), .sq(im_e)))
      re_e <- .div(re_e, mag); im_e <- .div(.neg(im_e), mag)
    }
    if (is_number_type(re_e) && is_number_type(im_e)) {
      base <- if (n >= 0) b else .add(list(re_e, .mul(list(im_e, S_I))))
      expr <- expand_multinomial(.pow(base, .int(abs(n))))
      if (!identical(expr$key, e$key)) return(.as_real_imag(expr))
    }
    N <- abs(n)
    re_t <- list(); im1 <- list(); im3 <- list()
    for (k in 0:N) {
      cc <- .int(choose(N, k))
      if (k %% 2L == 0L) re_t[[length(re_t) + 1L]] <- .mul(list(cc, .pow(re_e, .int(N - k)), .pow(.mul(list(S_I, im_e)), .int(k))))
      else if (k %% 4L == 1L) im1[[length(im1) + 1L]] <- .mul(list(cc, .pow(re_e, .int(N - k)), .pow(im_e, .int(k))))
      else im3[[length(im3) + 1L]] <- .mul(list(cc, .pow(re_e, .int(N - k)), .pow(.neg(im_e), .int(k))))
    }
    return(list(.add(re_t), .add(list(.add(im1), .add(im3)))))
  }
  if (ex$type == "Rational") {
    ri <- .as_real_imag(b, deep = deep)
    re_e <- ri[[1L]]; im_e <- ri[[2L]]
    if (.ari_zero(im_e) && identical(ex$key, S_Half$key)) {
      if (isTRUE(ask_fact(re_e, "extended_nonnegative"))) return(list(e, S_Zero))
      if (isTRUE(ask_fact(re_e, "extended_nonpositive"))) return(list(S_Zero, .pow(.neg(b), ex)))
    }
    r <- .pow(.add(list(.pow(re_e, S_Two), .pow(im_e, S_Two))), S_Half)
    t <- fcall("atan2", list(im_e, re_e))
    rp <- .pow(r, ex); tp <- .mul(list(t, ex))
    return(list(.mul(list(rp, fcall("cos", list(tp)))), .mul(list(rp, fcall("sin", list(tp))))))
  }
  if (deep) {
    expanded <- expand(e)
    if (!is.null(ignore) && identical(expanded$key, ignore$key)) return(NULL)
    return(list(fcall("re", list(expanded)), fcall("im", list(expanded))))
  }
  list(fcall("re", list(e)), fcall("im", list(e)))
}

# re.eval and im.eval
.re_im_eval_sympy <- function(which) {
  force(which)
  function(args) {
    a <- args[[1L]]
    if (is_nan(a) || is_const(a, "zoo")) return(S_NaN)
    if (isTRUE(ask_fact(a, "extended_real"))) return(if (which == 1L) a else S_Zero)
    if (isTRUE(ask_fact(a, "imaginary")) || isTRUE(ask_fact(.mul(list(S_I, a)), "extended_real")))
      return(if (which == 1L) S_Zero else .mul(list(.neg(S_I), a)))
    if (a$type == "Function" && a$name == "conjugate")
      return(if (which == 1L) fcall("re", list(a$args[[1L]])) else .neg(fcall("im", list(a$args[[1L]]))))
    included <- list(); reverted <- list(); excluded <- list()
    terms <- if (a$type == "Add") .sympy_args(a) else list(a)
    for (term in terms) {
      coeff <- as_coefficient(term, S_I)
      if (!is.null(coeff)) {
        if (!isTRUE(ask_fact(coeff, "extended_real"))) reverted[[length(reverted) + 1L]] <- coeff
        else if (which == 2L) excluded[[length(excluded) + 1L]] <- coeff
        next
      }
      has_i <- .contains_key(term, S_I$key)
      if (which == 1L && !has_i && isTRUE(ask_fact(term, "extended_real"))) {
        excluded[[length(excluded) + 1L]] <- term
        next
      }
      if (which == 2L && !(has_i || !isTRUE(ask_fact(term, "extended_real")))) next
      ri <- .as_real_imag(term, ignore = a)
      if (!is.null(ri)) excluded[[length(excluded) + 1L]] <- ri[[which]]
      else included[[length(included) + 1L]] <- term
    }
    if (length(terms) == length(included)) return(NULL)
    A <- .add(included); B <- .add(reverted); C <- .add(excluded)
    if (which == 1L) .add(list(fcall("re", list(A)), .neg(fcall("im", list(B))), C))
    else .add(list(fcall("im", list(A)), fcall("re", list(B)), C))
  }
}

.re_eval_s <- .re_im_eval_sympy(1L)
.im_eval_s <- .re_im_eval_sympy(2L)
