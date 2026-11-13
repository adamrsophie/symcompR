# Special and combinatorial functions, after sympy/functions/special and
# sympy/functions/combinatorial.

.is_nonpos_int <- function(a) a$type == "Integer" && a$value <= 0
.is_pos_int <- function(a) a$type == "Integer" && a$value > 0
.int_val <- function(a) as.numeric(gmp::numerator(a$value))
.bigz_val <- function(a) gmp::numerator(a$value)

.half_int <- function(a) a$type == "Rational" && gmp::denominator(a$value) == 2

# ---------------------------------------------------------------------------
# Gamma family

def_function("gamma",
  eval = function(args) {
    a <- args[[1L]]
    if (a$type == "Integer") {
      if (a$value <= 0) return(S_ComplexInfinity)
      return(make_rational(gmp::as.bigq(gmp::factorialZ(.int_val(a) - 1))))
    }
    if (.half_int(a)) {
      # gamma(n + 1/2) = (2n)! sqrt(pi) / (4^n n!)
      n <- as.numeric(.bigq_floor(a$value))
      if (n >= 0) {
        cf <- gmp::as.bigq(gmp::factorialZ(2 * n), gmp::as.bigz(4)^n * gmp::factorialZ(n))
      } else {
        m <- -n
        cf <- gmp::as.bigq(gmp::as.bigz((-1)^m) * gmp::as.bigz(4)^m * gmp::factorialZ(m), gmp::factorialZ(2 * m))
      }
      return(.mul(list(make_rational(cf), .sqrt(S_Pi))))
    }
    if (a$type == "Float") return(.float_eval("gamma", args))
    if (is_const(a, "oo")) return(S_Infinity)
    NULL
  },
  diff = function(args, i) .mul(list(fcall("gamma", args), fcall("polygamma", list(S_Zero, args[[1L]])))),
  numeric = function(x) if (is.complex(x)) .cgamma(x) else gamma(x),
  facts = function(args, afs) {
    f <- afs[[1L]]
    out <- c(commutative = TRUE)
    if (isTRUE(f[["positive"]])) out["positive"] <- TRUE
    if (isTRUE(f[["real"]]) && isTRUE(f[["noninteger"]])) out["real"] <- TRUE
    if (isTRUE(f[["integer"]]) && isTRUE(f[["positive"]])) out["integer"] <- TRUE
    # gamma._eval_is_positive: floor(x).is_even for a noninteger x
    if (!isTRUE(f[["positive"]]) && isTRUE(f[["noninteger"]])) {
      fe <- ask_fact(fcall("floor", list(args[[1L]])), "even")
      if (!is.na(fe)) out["positive"] <- fe
    }
    out
  },
  latex = function(a) paste0("\\Gamma\\left(", a[1], "\\right)"))

# Lanczos approximation for complex gamma.
.cgamma <- function(z) {
  g <- 7
  p <- c(0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313,
         -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6,
         1.5056327351493116e-7)
  if (Re(z) < 0.5) return(pi / (sin(pi * z) * .cgamma(1 - z)))
  z <- z - 1
  x <- p[1]
  for (i in 1:8) x <- x + p[i + 1] / (z + i)
  t <- z + g + 0.5
  sqrt(2 * pi) * t^(z + 0.5) * exp(-t) * x
}

def_function("loggamma",
  eval = function(args) {
    a <- args[[1L]]
    if (a$type == "Float") return(.float_eval("loggamma", args))
    if (isTRUE(ask_fact(a, "integer"))) {
      if (isTRUE(ask_fact(a, "nonpositive"))) return(S_Infinity)
      if (isTRUE(ask_fact(a, "positive"))) return(fcall("log", list(fcall("gamma", args))))
    } else if (a$type == "Rational" && gmp::denominator(a$value) == 2 && gmp::numerator(a$value) > 0) {
      p <- .int(gmp::numerator(a$value))
      return(fcall("log", list(.mul(list(.sqrt(S_Pi), .pow(S_Two, .sub(S_One, p)), fcall("gamma", list(p)),
                                         .inv(fcall("gamma", list(.mul(list(.add(list(p, S_One)), S_Half))))))))))
    }
    if (is_const(a, "oo")) return(S_Infinity)
    # abs(z) is oo
    if (is_inf_const(a) || isTRUE(is_infinite(a))) return(S_ComplexInfinity)
    if (is_nan(a)) return(S_NaN)
    NULL
  },
  diff = function(args, i) fcall("polygamma", list(S_Zero, args[[1L]])),
  # loggamma._eval_is_real
  facts = function(args, afs) {
    f <- afs[[1L]]
    out <- c(commutative = TRUE)
    if (isTRUE(f[["positive"]])) out["real"] <- TRUE
    else if (isTRUE(f[["nonpositive"]])) out["real"] <- FALSE
    out
  },
  numeric = function(x) {
    if (Im(x) != 0) stop("complex loggamma is not implemented")
    x <- Re(x)
    # principal branch on the negative axis, as SymPy/mpmath
    if (x < 0 && x != trunc(x)) return(complex(real = lgamma(x), imaginary = -pi * ceiling(-x)))
    lgamma(x)
  })

def_function("polygamma", nargs = 2L,
  eval = function(args) .polygamma_eval(args[[1L]], args[[2L]]),
  diff = function(args, i) {
    if (i == 2L) return(fcall("polygamma", list(.add(list(args[[1L]], S_One)), args[[2L]])))
    NULL
  },
  numeric = function(n, z) psigamma(Re(z), Re(n)),
  str = function(e) paste0("polygamma(", .str(e$args[[1L]]), ", ", .str(e$args[[2L]]), ")"))

def_function("beta", nargs = 2L,
  eval = function(args) {
    a <- args[[1L]]; b <- args[[2L]]
    if (is_number_type(a) && is_number_type(b) && (a$type != "Float" && b$type != "Float") &&
        (.is_pos_int(a) || .half_int(a)) && (.is_pos_int(b) || .half_int(b)))
      return(.div(.mul(list(fcall("gamma", list(a)), fcall("gamma", list(b)))), fcall("gamma", list(.add(list(a, b))))))
    if (a$type == "Float" || b$type == "Float") return(.float_eval("beta", args))
    # SymPy keeps the arguments in the order given
    NULL
  },
  diff = function(args, i) {
    a <- args[[1L]]; b <- args[[2L]]
    .mul(list(fcall("beta", args), .sub(fcall("polygamma", list(S_Zero, args[[i]])),
                                         fcall("polygamma", list(S_Zero, .add(list(a, b)))))))
  },
  numeric = function(a, b) beta(Re(a), Re(b)))

def_function("lowergamma", nargs = 2L,
  eval = function(args) {
    # lowergamma.eval
    a <- args[[1L]]; x <- args[[2L]]
    if (is_zero_num(x)) return(S_Zero)
    bf <- extract_branch_factor(x); nx <- bf[[1L]]; n <- bf[[2L]]
    ai <- isTRUE(ask_fact(a, "integer"))
    if (ai && isTRUE(ask_fact(a, "positive"))) {
      nx <- unpolarify(x)
      if (!identical(nx$key, x$key)) return(fcall("lowergamma", list(a, nx)))
    } else if (ai && isTRUE(ask_fact(a, "nonpositive"))) {
      if (!is_zero_num(n))
        return(.add(list(.div(.mul(list(S_Two, S_Pi, S_I, n, .pow(S_NegativeOne, .neg(a)))), fcall("factorial", list(.neg(a)))),
                         fcall("lowergamma", list(a, nx)))))
    } else if (!is_zero_num(n)) {
      return(.mul(list(fcall("exp", list(.mul(list(S_Two, S_Pi, S_I, n, a)))), fcall("lowergamma", list(a, nx)))))
    }
    if (is_number_type(a)) {
      if (is_one_num(a)) return(.sub(S_One, fcall("exp", list(.neg(x)))))
      if (identical(a$key, S_Half$key)) return(.mul(list(.sqrt(S_Pi), fcall("erf", list(.sqrt(x))))))
      if (a$type == "Integer" || .mul(list(S_Two, a))$type == "Integer") {
        b <- .sub(a, S_One)
        if (isTRUE(ask_fact(b, "positive"))) {
          if (a$type == "Integer") {
            fb <- fcall("factorial", list(b))
            return(.sub(fb, .mul(list(fcall("exp", list(.neg(x))), fb,
              .add(lapply(seq_len(.int_val(a)) - 1L, function(k) .div(.pow(x, .int(k)), fcall("factorial", list(.int(k))))))))))
          }
          kmax <- as.integer(.int_val(.add(list(a, S_Half))))
          terms <- lapply(seq_len(kmax - 1L), function(k)
            .div(.pow(x, .sub(.int(k), S_Half)), fcall("gamma", list(.add(list(S_Half, .int(k)))))))
          inner <- .sub(.div(fcall("lowergamma", list(S_Half, x)), .sqrt(S_Pi)), .mul(list(fcall("exp", list(.neg(x))), .add(terms))))
          return(.mul(list(fcall("gamma", list(a)), inner)))
        }
        if (a$type != "Integer") {
          kmax <- as.integer(.int_val(.sub(.rat(3L, 2L), a)))
          terms <- lapply(seq_len(kmax - 1L), function(k)
            .div(.mul(list(.pow(x, .add(list(.int(k), a, S_NegativeOne))), fcall("gamma", list(a)))), fcall("gamma", list(.add(list(a, .int(k)))))))
          first <- .div(.mul(list(.pow(S_NegativeOne, .sub(S_Half, a)), S_Pi, fcall("erf", list(.sqrt(x))))), fcall("gamma", list(.sub(S_One, a))))
          return(.add(list(first, .mul(list(fcall("exp", list(.neg(x))), .add(terms))))))
        }
      }
    }
    if (isTRUE(ask_fact(x, "zero"))) return(S_Zero)
    NULL
  },
  diff = function(args, i) {
    if (i == 2L) return(.mul(list(.pow(args[[2L]], .sub(args[[1L]], S_One)), fcall("exp", list(.neg(args[[2L]]))))))
    NULL
  },
  numeric = function(s, x) stats::pgamma(Re(x), Re(s)) * gamma(Re(s)))

def_function("uppergamma", nargs = 2L,
  eval = function(args) {
    # uppergamma.eval
    a <- args[[1L]]; z <- args[[2L]]
    if (is_number_type(z) || is_nan(z) || is_const(z, "oo")) {
      if (is_nan(z)) return(S_NaN)
      if (is_const(z, "oo")) return(S_Zero)
      if (is_zero_num(z) && isTRUE(ask_fact(fcall("re", list(a)), "positive"))) return(fcall("gamma", list(a)))
    }
    bf <- extract_branch_factor(z); nx <- bf[[1L]]; n <- bf[[2L]]
    ai <- isTRUE(ask_fact(a, "integer"))
    if (ai && isTRUE(ask_fact(a, "positive"))) {
      nx <- unpolarify(z)
      if (!identical(nx$key, z$key)) return(fcall("uppergamma", list(a, nx)))
    } else if (ai && isTRUE(ask_fact(a, "nonpositive"))) {
      if (!is_zero_num(n))
        return(.add(list(.neg(.div(.mul(list(S_Two, S_Pi, S_I, n, .pow(S_NegativeOne, .neg(a)))), fcall("factorial", list(.neg(a))))),
                         fcall("uppergamma", list(a, nx)))))
    } else if (!is_zero_num(n)) {
      e2 <- fcall("exp", list(.mul(list(S_Two, S_Pi, S_I, n, a))))
      return(.add(list(.mul(list(fcall("gamma", list(a)), .sub(S_One, e2))), .mul(list(e2, fcall("uppergamma", list(a, nx)))))))
    }
    if (is_number_type(a)) {
      if (is_zero_num(a) && isTRUE(ask_fact(z, "positive"))) return(.neg(fcall("Ei", list(.neg(z)))))
      if (is_one_num(a)) return(fcall("exp", list(.neg(z))))
      if (identical(a$key, S_Half$key)) return(.mul(list(.sqrt(S_Pi), fcall("erfc", list(.sqrt(z))))))
      if (a$type == "Integer" || .mul(list(S_Two, a))$type == "Integer") {
        b <- .sub(a, S_One)
        if (isTRUE(ask_fact(b, "positive"))) {
          if (a$type == "Integer")
            return(.mul(list(fcall("exp", list(.neg(z))), fcall("factorial", list(b)),
              .add(lapply(seq_len(.int_val(a)) - 1L, function(k) .div(.pow(z, .int(k)), fcall("factorial", list(.int(k)))))))))
          kmax <- as.integer(.int_val(.sub(a, S_Half)))
          return(.add(list(.mul(list(fcall("gamma", list(a)), fcall("erfc", list(.sqrt(z))))),
            .mul(list(.pow(S_NegativeOne, .sub(a, .rat(3L, 2L))), fcall("exp", list(.neg(z))), .sqrt(z),
              .add(lapply(seq_len(kmax) - 1L, function(k)
                .div(.mul(list(fcall("gamma", list(.sub(.neg(S_Half), .int(k)))), .pow(.neg(z), .int(k)))), fcall("gamma", list(.sub(S_One, a)))))))))))
        }
        if (b$type == "Integer") return(.mul(list(fcall("expint", list(.neg(b), z)), .pow(unpolarify(z), .add(list(b, S_One))))))
        if (a$type != "Integer") {
          kmax <- as.integer(.int_val(.sub(S_Half, a)))
          return(.sub(.div(.mul(list(.pow(S_NegativeOne, .sub(S_Half, a)), S_Pi, fcall("erfc", list(.sqrt(z))))), fcall("gamma", list(.sub(S_One, a)))),
            .mul(list(.pow(z, a), fcall("exp", list(.neg(z))), .add(lapply(seq_len(kmax) - 1L, function(k)
              .div(.mul(list(.pow(z, .int(k)), fcall("gamma", list(a)))), fcall("gamma", list(.add(list(a, .int(k), S_One)))))))))))
        }
      }
    }
    if (isTRUE(ask_fact(a, "zero")) && isTRUE(ask_fact(z, "positive"))) return(.neg(fcall("Ei", list(.neg(z)))))
    if (isTRUE(ask_fact(z, "zero")) && isTRUE(ask_fact(fcall("re", list(a)), "positive"))) return(fcall("gamma", list(a)))
    NULL
  },
  diff = function(args, i) {
    if (i == 2L) return(.neg(.mul(list(.pow(args[[2L]], .sub(args[[1L]], S_One)), fcall("exp", list(.neg(args[[2L]])))))))
    NULL
  },
  numeric = function(s, x) stats::pgamma(Re(x), Re(s), lower.tail = FALSE) * gamma(Re(s)))

# ---------------------------------------------------------------------------
# Factorials and binomials

def_function("factorial",
  eval = function(args) {
    n <- args[[1L]]
    if (n$type == "Integer") {
      if (n$value < 0) return(S_ComplexInfinity)
      if (n$value > 100000) return(NULL)
      return(make_rational(gmp::as.bigq(gmp::factorialZ(.int_val(n)))))
    }
    if (n$type == "Float") return(.float_eval("factorial", args))
    if (is_const(n, "oo")) return(S_Infinity)
    NULL
  },
  diff = function(args, i) .mul(list(fcall("gamma", list(.add(list(args[[1L]], S_One)))),
                                     fcall("polygamma", list(S_Zero, .add(list(args[[1L]], S_One)))))),
  numeric = function(n) if (is.complex(n)) .cgamma(n + 1) else gamma(n + 1),
  facts = function(args, afs) {
    f <- afs[[1L]]
    out <- c(commutative = TRUE)
    if (isTRUE(f[["integer"]]) && isTRUE(f[["nonnegative"]])) {
      out["integer"] <- TRUE
      out["positive"] <- TRUE
    }
    out
  },
  str = function(e) paste0("factorial(", .str(e$args[[1L]]), ")"),
  latex = function(a) paste0(a[1], "!"))

def_function("factorial2",
  eval = function(args) {
    n <- args[[1L]]
    if (n$type == "Integer") {
      v <- .int_val(n)
      if (v < -1) return(NULL)
      if (v <= 0) return(S_One)
      r <- gmp::as.bigz(1)
      for (k in seq(v, 1, by = -2)) r <- r * k
      return(make_rational(gmp::as.bigq(r)))
    }
    NULL
  },
  numeric = function(n) prod(seq(Re(n), 1, by = -2)))

def_function("subfactorial",
  eval = function(args) {
    n <- args[[1L]]
    # subfactorial._eval: z1, z2 = 1, 0; z1, z2 = z2, (i - 1)*(z2 + z1)
    if (n$type == "Integer" && n$value >= 0) {
      v <- .int_val(n)
      if (v == 0) return(S_One)
      a <- gmp::as.bigz(1); b <- gmp::as.bigz(0)
      for (k in seq_len(v)[-1L]) {
        c <- (k - 1) * (a + b)
        a <- b; b <- c
      }
      return(make_rational(gmp::as.bigq(b)))
    }
    if (is_const(n, "nan")) return(S_NaN)
    if (is_const(n, "oo")) return(S_Infinity)
    NULL
  },
  numeric = function(n) round(gamma(Re(n) + 1) / exp(1)))

def_function("binomial", nargs = 2L,
  # binomial.eval and binomial._eval
  eval = function(args) {
    n <- args[[1L]]; k <- args[[2L]]
    d <- .sub(n, k)
    nn <- ask_fact(n, "nonnegative"); ni <- ask_fact(n, "integer")
    nn_or <- isTRUE(nn) || isFALSE(ni)
    if (isTRUE(ask_fact(k, "zero")) || (nn_or && isTRUE(ask_fact(d, "zero")))) return(S_One)
    if (isTRUE(ask_fact(.sub(k, S_One), "zero")) || (nn_or && isTRUE(ask_fact(.sub(d, S_One), "zero")))) return(n)
    if (isTRUE(ask_fact(k, "integer"))) {
      if (isTRUE(ask_fact(k, "negative")) || (isTRUE(nn) && isTRUE(ni) && isTRUE(ask_fact(d, "negative")))) return(S_Zero)
      if (!length(free_symbols(n)) && k$type == "Integer") {
        kv <- .bigz_val(k)
        if (n$type == "Integer" && n$value >= 0) {
          nv <- .bigz_val(n)
          if (kv > nv) return(S_Zero)
          return(make_rational(gmp::as.bigq(gmp::chooseZ(nv, kv))))
        }
        if (kv > 5000) return(NULL)
        res <- S_One; dd <- d
        for (i in seq_len(as.integer(as.character(kv)))) { dd <- .add(list(dd, S_One)); res <- .mul(list(res, dd)) }
        return(expand(.div(res, make_rational(gmp::as.bigq(gmp::factorialZ(as.integer(as.character(kv))))))))
      }
    } else if (isFALSE(nn) && isTRUE(ni)) {
      return(S_ComplexInfinity)
    } else if (!length(free_symbols(k))) {
      g <- function(v) fcall("gamma", list(v))
      return(.div(g(.add(list(n, S_One))), .mul(list(g(.add(list(k, S_One))), g(.add(list(.sub(n, k), S_One)))))))
    }
    NULL
  },
  diff = function(args, i) {
    n <- args[[1L]]; k <- args[[2L]]
    if (i == 1L) return(.mul(list(fcall("binomial", args), .sub(fcall("polygamma", list(S_Zero, .add(list(n, S_One)))),
                                                                fcall("polygamma", list(S_Zero, .add(list(.sub(n, k), S_One))))))))
    .mul(list(fcall("binomial", args), .sub(fcall("polygamma", list(S_Zero, .add(list(.sub(n, k), S_One)))),
                                             fcall("polygamma", list(S_Zero, .add(list(k, S_One)))))))
  },
  numeric = function(n, k) choose(Re(n), Re(k)),
  facts = function(args, afs) {
    out <- c(commutative = TRUE)
    if (.all_true(afs, "integer")) out["integer"] <- TRUE
    if (isTRUE(afs[[1L]][["integer"]]) && isTRUE(afs[[2L]][["integer"]]) && isTRUE(afs[[1L]][["nonnegative"]]))
      out["nonnegative"] <- TRUE
    out
  },
  latex = function(a) paste0("{\\binom{", a[1], "}{", a[2], "}}"))

def_function("RisingFactorial", nargs = 2L,
  # RisingFactorial.eval
  eval = function(args) {
    x <- args[[1L]]; k <- args[[2L]]
    if (is_nan(x) || is_nan(k)) return(S_NaN)
    if (is_one_num(x)) return(fcall("factorial", list(k)))
    if (k$type == "Integer") {
      kv <- .int_val(k)
      if (kv == 0) return(S_One)
      # ponytail: products stop at |k| = 5000, larger k stay unevaluated
      if (abs(kv) > 5000) return(NULL)
      if (kv > 0) {
        if (is_const(x, "oo")) return(S_Infinity)
        if (is_const(x, "-oo")) return(if (kv %% 2 == 1) S_NegativeInfinity else S_Infinity)
        return(.mul(lapply(0:(kv - 1), function(j) .add(list(x, .int(j))))))
      }
      if (is_inf_const(x)) return(S_Infinity)
      return(.inv(.mul(lapply(1:(-kv), function(j) .sub(x, .int(j))))))
    }
    if (isFALSE(ask_fact(k, "integer")) && isTRUE(ask_fact(x, "integer")) && isTRUE(ask_fact(x, "negative"))) return(S_Zero)
    NULL
  },
  numeric = function(x, k) gamma(Re(x) + Re(k)) / gamma(Re(x)),
  str = function(e) paste0("RisingFactorial(", .str(e$args[[1L]]), ", ", .str(e$args[[2L]]), ")"))

def_function("FallingFactorial", nargs = 2L,
  eval = function(args) {
    x <- args[[1L]]; k <- args[[2L]]
    if (k$type == "Integer") {
      kv <- .int_val(k)
      if (kv == 0) return(S_One)
      if (kv > 0 && kv <= 100) return(.mul(lapply(0:(kv - 1), function(j) .sub(x, .int(j)))))
      if (kv < 0 && kv >= -100) return(.inv(.mul(lapply(1:(-kv), function(j) .add(list(x, .int(j)))))))
    }
    NULL
  },
  numeric = function(x, k) gamma(Re(x) + 1) / gamma(Re(x) - Re(k) + 1),
  str = function(e) paste0("FallingFactorial(", .str(e$args[[1L]]), ", ", .str(e$args[[2L]]), ")"))

# ---------------------------------------------------------------------------
# Combinatorial numbers

.bernoulli_cache <- new.env()
.bernoulli_num <- function(n) {
  key <- as.character(n)
  if (!is.null(.bernoulli_cache[[key]])) return(.bernoulli_cache[[key]])
  # Akiyama-Tanigawa
  A <- vector("list", n + 1L)
  for (m in 0:n) {
    A[[m + 1L]] <- gmp::as.bigq(1, m + 1)
    if (m >= 1) for (j in m:1) A[[j]] <- j * (A[[j]] - A[[j + 1L]])
  }
  r <- A[[1L]]
  if (n == 1) r <- -r
  assign(key, r, envir = .bernoulli_cache)
  r
}

.seq_numbers <- function(name, fn, numeric = NULL, nargs = 1L) {
  def_function(name, nargs = nargs,
               eval = function(args) {
                 n <- args[[1L]]
                 if (n$type == "Integer" && n$value >= 0 && n$value < 5000 &&
                     (all(vapply(args[-1L], is_exact_number, TRUE)) || name == "bernoulli"))
                   return(fn(.int_val(n), args))
                 NULL
               },
               numeric = numeric)
}

.seq_numbers("bernoulli", function(n, args) {
  # SymPy's numbers take B1 = +1/2; the polynomials keep B1(x) = x - 1/2
  if (length(args) == 1L) return(if (n == 1) S_Half else make_rational(.bernoulli_num(n)))
  # bernoulli_poly(n, x): the polynomial in the generator x, as Poly.as_expr
  x <- args[[2L]]
  .add(lapply(0:n, function(k) {
    c <- gmp::as.bigq(gmp::chooseZ(n, k)) * .bernoulli_num(k)
    if (c == 0) S_Zero else .mul(list(make_rational(c), .pow(x, .int(n - k))))
  }))
}, nargs = 1:2)

.seq_numbers("fibonacci", function(n, args) {
  a <- gmp::as.bigz(0); b <- gmp::as.bigz(1)
  if (n == 0) return(S_Zero)
  for (i in seq_len(n - 1)) { t <- a + b; a <- b; b <- t }
  make_rational(gmp::as.bigq(b))
}, numeric = function(n) round(((1 + sqrt(5)) / 2)^Re(n) / sqrt(5)))

.seq_numbers("lucas", function(n, args) {
  a <- gmp::as.bigz(2); b <- gmp::as.bigz(1)
  if (n == 0) return(S_Two)
  for (i in seq_len(n - 1)) { t <- a + b; a <- b; b <- t }
  make_rational(gmp::as.bigq(b))
})

.seq_numbers("tribonacci", function(n, args) {
  v <- gmp::as.bigz(c(0, 1, 1))
  if (n < 3) return(.int(c(0, 1, 1)[n + 1]))
  for (i in 3:n) v <- c(v[2], v[3], v[1] + v[2] + v[3])
  make_rational(gmp::as.bigq(v[3]))
})

.seq_numbers("catalan", function(n, args) {
  make_rational(gmp::as.bigq(gmp::chooseZ(2 * n, n), n + 1))
}, numeric = function(n) gamma(2 * Re(n) + 1) / (gamma(Re(n) + 2) * gamma(Re(n) + 1)))

# bell.eval: Bell numbers, Bell polynomials B_n(x) and the incomplete Bell
# polynomials B_{n,k}(x_1, ..., x_{n-k+1})
def_function("bell", nargs = 1:3, eval = function(args) {
  n <- args[[1L]]
  if (is_const(n, "oo")) {
    if (length(args) == 1L) return(S_Infinity)
    stop(.sympy_condition("ValueError", "Bell polynomial is not defined"))
  }
  if (isTRUE(ask_fact(n, "negative")) || isFALSE(ask_fact(n, "integer")))
    stop(.sympy_condition("ValueError", "a non-negative integer expected"))
  if (n$type != "Integer") return(NULL)
  m <- .int_val(n)
  if (length(args) == 1L) {
    # _bell with recurrence_memo([1, 1])
    prev <- list(gmp::as.bigz(1), gmp::as.bigz(1))
    if (m >= 2L) for (j in 2:m) {
      s <- gmp::as.bigz(1); a <- gmp::as.bigz(1)
      for (k in seq_len(j - 1L)) { a <- gmp::divq.bigz(a * (j - k), k); s <- s + a * prev[[k + 1L]] }
      prev[[j + 1L]] <- s
    }
    return(make_rational(gmp::as.bigq(prev[[m + 1L]])))
  }
  if (length(args) == 2L) {
    # _bell_poly_coeffs: each row from the previous one, leading coefficient first
    row <- list(gmp::as.bigz(1))
    if (m >= 1L) for (j in seq_len(m)) {
      i <- j; nr <- list(gmp::as.bigz(1))
      if (length(row) >= 2L) for (t in seq_len(length(row) - 1L)) { i <- i - 1L; nr[[length(nr) + 1L]] <- row[[t]] * i + row[[t + 1L]] }
      row <- c(nr, list(gmp::as.bigz(0)))
    }
    d <- length(row) - 1L
    return(.add(lapply(seq_along(row), function(t) .mul(list(make_rational(gmp::as.bigq(row[[t]])), .pow(args[[2L]], .int(d - t + 1L)))))))
  }
  if (args[[2L]]$type != "Integer") return(NULL)
  syms <- if (args[[3L]]$type == "Tuple") args[[3L]]$args else list(args[[3L]])
  .bell_incomplete(m, .int_val(args[[2L]]), syms)
})

# harmonic.eval
def_function("harmonic", nargs = 1:2, eval = function(args) {
  n <- args[[1L]]
  if (length(args) == 2L && is_one_num(args[[2L]])) return(fcall("harmonic", list(n)))
  m <- if (length(args) == 2L) args[[2L]] else S_One
  if (isTRUE(ask_fact(n, "zero"))) return(S_Zero)
  if (isTRUE(ask_fact(m, "zero"))) return(n)
  if (is_const(n, "oo")) {
    if (isTRUE(ask_fact(m, "negative"))) return(S_NaN)
    if (isTRUE(ask_fact(.sub(m, S_One), "nonpositive"))) return(S_Infinity)
    if (isTRUE(ask_fact(.sub(m, S_One), "positive"))) return(fcall("zeta", list(m)))
    return(NULL)
  }
  if (m$type == "Integer" && m$value <= 0) {
    k <- .sub(S_One, m)
    return(.div(.sub(fcall("bernoulli", list(k, .add(list(n, S_One)))), fcall("bernoulli", list(k))), k))
  }
  if (n$type == "Integer") {
    if (n$value < 0 && (isFALSE(ask_fact(m, "integer")) || isFALSE(ask_fact(m, "nonpositive"))))
      return(if (is_one_num(m)) S_ComplexInfinity else S_NaN)
    # ponytail: exact sums stop at n = 5000, larger n stay unevaluated
    if (n$value >= 0 && n$value < 5000) {
      nn <- .int_val(n)
      if (m$type == "Integer") {
        s <- gmp::as.bigq(0); mm <- as.numeric(m$value)
        for (k in seq_len(nn)) s <- s + .bigq_pow_int(gmp::as.bigq(k), -mm)
        return(make_rational(s))
      }
      return(.add(lapply(seq_len(nn), function(k) .pow(.int(k), .neg(m)))))
    }
  }
  NULL
})

.seq_numbers("euler", function(n, args) {
  # Euler numbers via the recurrence E_n = -sum_{k<n, even} C(n,k) E_k
  if (n %% 2 == 1) return(S_Zero)
  E <- list(gmp::as.bigz(1))
  if (n == 0) return(S_One)
  for (m in seq(2, n, by = 2)) {
    s <- gmp::as.bigz(0)
    for (k in seq(0, m - 2, by = 2)) s <- s + gmp::chooseZ(m, k) * E[[k / 2 + 1]]
    E[[m / 2 + 1]] <- -s
  }
  make_rational(gmp::as.bigq(E[[n / 2 + 1]]))
})

.seq_numbers("genocchi", function(n, args) {
  if (n == 0) return(S_Zero)
  make_rational(2 * (1 - gmp::as.bigq(2)^n) * .bernoulli_num(n))
})

.seq_numbers("partition", function(n, args) make_rational(gmp::as.bigq(.npartitions(n))))

# motzkin.eval: Integer(cls._motzkin(n - 1)) for an integer n >= 0, where
# _motzkin(-1) reads the last entry of a fresh recurrence_memo cache [1, 1]
def_function("motzkin", eval = function(args) {
  n <- args[[1L]]
  if (n$type != "Integer" || n$value < 0)
    .op_value_error("The provided number must be a positive integer")
  k <- .int_val(n) - 1L
  M <- list(gmp::as.bigz(1), gmp::as.bigz(1))
  if (k < 2L) return(S_One)
  for (m in 2:k) M[[m + 1L]] <- ((2 * m + 1) * M[[m]] + (3 * m - 3) * M[[m - 1L]]) %/% (m + 2)
  make_rational(gmp::as.bigq(M[[k + 1L]]))
})

# ---------------------------------------------------------------------------
# Error functions and exponential integrals

def_function("erf",
  eval = function(args) {
    a <- args[[1L]]
    if (is_zero_num(a)) return(S_Zero)
    if (is_const(a, "oo")) return(S_One)
    if (is_const(a, "-oo")) return(S_NegativeOne)
    if (is_nan(a)) return(S_NaN)
    if (a$type == "Function" && a$name == "erfinv") return(a$args[[1L]])
    if (a$type == "Function" && a$name == "erfcinv") return(.sub(S_One, a$args[[1L]]))
    if (a$type == "Function" && a$name == "erf2inv" && is_zero_num(a$args[[1L]])) return(a$args[[2L]])
    if (a$type == "Float" && (.float_prec(a) != 53 || !is.null(a$extra$mp))) return(.float_eval("erf", args))
    if (a$type == "Float") return(Float(2 * stats::pnorm(a$value * sqrt(2)) - 1))
    # erf.eval: erf(I*oo) is I*oo
    t <- extract_multiplicatively(a, S_I)
    if (!is.null(t) && (is_const(t, "oo") || is_const(t, "-oo"))) return(a)
    .odd_eval("erf", a)
  },
  diff = function(args, i) .mul(list(S_Two, .inv(.sqrt(S_Pi)), fcall("exp", list(.neg(.sq(args[[1L]])))))),
  numeric = function(x) if (Im(x) == 0) 2 * stats::pnorm(Re(x) * sqrt(2)) - 1 else 1 - .erfc_complex(x),
  facts = function(args, afs) {
    out <- c(commutative = TRUE)
    if (isTRUE(afs[[1L]][["real"]])) { out["real"] <- TRUE; out["finite"] <- TRUE }
    out
  })

def_function("erfc",
  eval = function(args) {
    a <- args[[1L]]
    if (is_zero_num(a)) return(S_One)
    if (is_const(a, "oo")) return(S_Zero)
    if (is_nan(a)) return(S_NaN)
    if (a$type == "Function" && a$name == "erfinv") return(.sub(S_One, a$args[[1L]]))
    if (a$type == "Function" && a$name == "erfcinv") return(a$args[[1L]])
    if (a$type == "Float" && (.float_prec(a) != 53 || !is.null(a$extra$mp))) return(.float_eval("erfc", args))
    if (a$type == "Float") return(Float(2 * stats::pnorm(-a$value * sqrt(2))))
    # erfc.eval: erfc(I*oo) is -I*oo, erfc(-x) is 2 - erfc(x)
    t <- extract_multiplicatively(a, S_I)
    if (!is.null(t) && (is_const(t, "oo") || is_const(t, "-oo"))) return(.neg(a))
    if (could_extract_minus_sign(a)) return(.sub(S_Two, fcall("erfc", list(.neg(a)))))
    NULL
  },
  diff = function(args, i) .mul(list(.int(-2), .inv(.sqrt(S_Pi)), fcall("exp", list(.neg(.sq(args[[1L]])))))),
  numeric = function(x) if (Im(x) == 0) 2 * stats::pnorm(-Re(x) * sqrt(2)) else .erfc_complex(x))

# erfc of a complex argument: Maclaurin series of erf near zero, the
# Laplace continued fraction for large |z|, and erfc(-z) = 2 - erfc(z).
.erfc_complex <- function(z) {
  if (Re(z) < 0) return(2 - .erfc_complex(-z))
  if (Mod(z) < 2) {
    s <- 0; term <- z; n <- 0
    repeat {
      add <- term / (2 * n + 1)
      s <- s + add
      if (Mod(add) < 1e-17 * Mod(s)) break
      n <- n + 1
      term <- -term * z^2 / n
    }
    return(1 - 2 / sqrt(pi) * s)
  }
  f <- z
  for (k in 60:1) f <- z + (k / 2) / f
  exp(-z^2) / sqrt(pi) / f
}

def_function("erfi",
  eval = function(args) {
    a <- args[[1L]]
    if (is_zero_num(a)) return(S_Zero)
    if (is_nan(a)) return(S_NaN)
    if (is_const(a, "oo")) return(S_Infinity)
    r <- .odd_eval("erfi", a)
    if (!is.null(r)) return(r)
    nz <- extract_multiplicatively(a, S_I)
    if (!is.null(nz)) {
      if (is_const(nz, "oo")) return(S_I)
      if (nz$type == "Function" && nz$name == "erfinv") return(.mul(list(S_I, nz$args[[1L]])))
      if (nz$type == "Function" && nz$name == "erfcinv") return(.mul(list(S_I, .sub(S_One, nz$args[[1L]]))))
      if (nz$type == "Function" && nz$name == "erf2inv" && is_zero_num(nz$args[[1L]])) return(.mul(list(S_I, nz$args[[2L]])))
    }
    NULL
  },
  diff = function(args, i) .mul(list(S_Two, .inv(.sqrt(S_Pi)), fcall("exp", list(.sq(args[[1L]]))))),
  numeric = function(x) {
    x <- Re(x)
    stats::integrate(function(t) exp(t^2), 0, x)$value * 2 / sqrt(pi)
  })

def_function("erfinv",
  eval = function(args) {
    a <- args[[1L]]
    if (is_zero_num(a)) return(S_Zero)
    if (is_one_num(a)) return(S_Infinity)
    if (is_negone_num(a)) return(S_NegativeInfinity)
    if (is_nan(a)) return(S_NaN)
    if (a$type == "Function" && a$name == "erf" && isTRUE(ask_fact(a$args[[1L]], "extended_real"))) return(a$args[[1L]])
    nz <- extract_multiplicatively(a, S_NegativeOne)
    if (!is.null(nz) && nz$type == "Function" && nz$name == "erf" && isTRUE(ask_fact(nz$args[[1L]], "extended_real"))) return(.neg(nz$args[[1L]]))
    NULL
  },
  diff = function(args, i) .mul(list(S_Half, .sqrt(S_Pi), fcall("exp", list(.sq(fcall("erfinv", args)))))),
  numeric = function(x) stats::qnorm((Re(x) + 1) / 2) / sqrt(2))

def_function("erfcinv",
  eval = function(args) {
    a <- args[[1L]]
    if (is_one_num(a)) return(S_Zero)
    if (is_zero_num(a)) return(S_Infinity)
    NULL
  },
  diff = function(args, i) .mul(list(.rat(-1, 2), .sqrt(S_Pi), fcall("exp", list(.sq(fcall("erfcinv", args)))))),
  numeric = function(x) -stats::qnorm(Re(x) / 2) / sqrt(2))

def_function("Ei",
  eval = function(args) .Ei_eval(args[[1L]]),
  diff = function(args, i) .div(fcall("exp", args), args[[1L]]),
  numeric = function(x) {
    x <- Re(x)
    if (x < -1) {
      # Ei(x) = -E1(-x), continued fraction of E1
      t <- -x; f <- t + 201
      for (k in 100:1) f <- t + 2 * k - 1 - k^2 / f
      return(-exp(-t) / f)
    }
    if (x > 40) {
      # asymptotic series e^x/x * sum k!/x^k, truncated at its smallest term
      s <- 1; term <- 1
      for (k in 1:40) {
        nt <- term * k / x
        if (nt > term) break
        term <- nt
        s <- s + term
      }
      return(exp(x) / x * s)
    }
    s <- 0.57721566490153286 + log(abs(x))
    term <- 1
    for (k in 1:500) {
      term <- term * x / k
      s <- s + term / k
      if (abs(term / k) < 1e-17 * abs(s)) break
    }
    s
  })

def_function("li",
  eval = function(args) {
    a <- args[[1L]]
    if (is_zero_num(a)) return(S_Zero)
    if (is_one_num(a)) return(S_NegativeInfinity)
    if (is_const(a, "oo")) return(S_Infinity)
    NULL
  },
  diff = function(args, i) .inv(fcall("log", args)),
  numeric = function(x) .neval(fcall("Ei", list(Float(log(Re(x))))), list()))

def_function("Si",
  eval = function(args) .ti_eval("Si", args[[1L]]),
  diff = function(args, i) .div(fcall("sin", args), args[[1L]]),
  numeric = function(x) stats::integrate(function(t) ifelse(t == 0, 1, sin(t) / t), 0, Re(x))$value)

def_function("Ci",
  eval = function(args) .ti_eval("Ci", args[[1L]]),
  diff = function(args, i) .div(fcall("cos", args), args[[1L]]),
  numeric = function(x) 0.57721566490153286 + log(Re(x)) +
    stats::integrate(function(t) (cos(t) - 1) / t, 0, Re(x))$value)

def_function("Shi",
  eval = function(args) .ti_eval("Shi", args[[1L]]),
  diff = function(args, i) .div(fcall("sinh", args), args[[1L]]),
  numeric = function(x) stats::integrate(function(t) ifelse(t == 0, 1, sinh(t) / t), 0, Re(x))$value)

def_function("Chi",
  eval = function(args) .ti_eval("Chi", args[[1L]]),
  diff = function(args, i) .div(fcall("cosh", args), args[[1L]]),
  numeric = function(x) 0.57721566490153286 + log(Re(x)) +
    stats::integrate(function(t) (cosh(t) - 1) / t, 0, Re(x))$value)

def_function("expint", nargs = 2L,
  eval = function(args) .expint_eval(args[[1L]], args[[2L]]),
  diff = function(args, i) if (i == 2L) .neg(fcall("expint", list(.sub(args[[1L]], S_One), args[[2L]]))) else NULL,
  numeric = function(nu, z) {
    nu <- Re(nu); z <- Re(z)
    if (z > 0) return(stats::integrate(function(t) exp(-z * t) / t^nu, 1, Inf, rel.tol = 1e-12)$value)
    if (z < 0 && nu >= 1 && nu == trunc(nu)) {
      # principal branch: E1(-x) = -Ei(x) - i*pi, then E(n+1) = (exp(-z) - z*E(n))/n
      e <- complex(real = -.neval(fcall("Ei", list(Float(-z))), list()), imaginary = -pi)
      if (nu > 1) for (k in seq_len(nu - 1)) e <- (exp(-z) - z * e) / k
      return(e)
    }
    stop("expint() is not implemented for this argument")
  })

def_function("fresnels",
  eval = function(args) .fresnel_eval("fresnels", args[[1L]]),
  diff = function(args, i) fcall("sin", list(.mul(list(S_Half, S_Pi, .sq(args[[1L]]))))),
  numeric = function(x) stats::integrate(function(t) sin(pi * t^2 / 2), 0, Re(x))$value)

def_function("fresnelc",
  eval = function(args) .fresnel_eval("fresnelc", args[[1L]]),
  diff = function(args, i) fcall("cos", list(.mul(list(S_Half, S_Pi, .sq(args[[1L]]))))),
  numeric = function(x) stats::integrate(function(t) cos(pi * t^2 / 2), 0, Re(x))$value)

# ---------------------------------------------------------------------------
# Zeta functions

def_function("zeta", nargs = 1:2,
  eval = function(args) .zeta_eval(args),
  numeric = function(s, a = 1) .zeta_em(s, a),
  str = function(e) paste0("zeta(", paste(vapply(e$args, .str, ""), collapse = ", "), ")"))

# Hurwitz zeta by Euler-Maclaurin summation, valid for every s != 1.
.ZETA_B2J <- c(1/6, -1/30, 1/42, -1/30, 5/66, -691/2730, 7/6, -3617/510, 43867/798, -174611/330)
.zeta_em <- function(s, a = 1) {
  if (Im(s) == 0) s <- Re(s)
  if (Im(a) == 0) a <- Re(a)
  N <- 30
  tot <- sum((0:(N - 1) + a)^(-s)) + (N + a)^(1 - s) / (s - 1) + 0.5 * (N + a)^(-s)
  poch <- s
  for (j in seq_along(.ZETA_B2J)) {
    tot <- tot + .ZETA_B2J[j] / factorial(2 * j) * poch * (N + a)^(-s - 2 * j + 1)
    poch <- poch * (s + 2 * j - 1) * (s + 2 * j)
  }
  tot
}

def_function("dirichlet_eta",
  eval = function(args) {
    s <- args[[1L]]
    if (is_one_num(s)) return(fcall("log", list(S_Two)))
    z <- fcall("zeta", list(s))
    if (z$type != "Function") return(.mul(list(.sub(S_One, .pow(S_Two, .sub(S_One, s))), z)))
    NULL
  },
  numeric = function(s) (1 - 2^(1 - Re(s))) * .neval(fcall("zeta", list(Float(Re(s)))), list()))

def_function("polylog", nargs = 2L,
  eval = function(args) {
    # polylog.eval
    s <- args[[1L]]; z <- args[[2L]]
    if (!length(free_symbols(z))) {
      if (is_one_num(z)) return(fcall("zeta", list(s)))
      if (is_negone_num(z)) return(.neg(fcall("dirichlet_eta", list(s))))
      if (is_zero_num(z)) return(S_Zero)
      if (identical(s$key, S_Two$key)) {
        v <- .dilogtable()[[z$key]]
        if (!is.null(v)) return(v)
      }
    }
    if (isTRUE(ask_fact(z, "zero"))) return(S_Zero)
    zone <- .equals_one(z)
    if (isTRUE(zone)) return(fcall("zeta", list(s)))
    if (identical(zone, FALSE)) {
      if (is_zero_num(s)) return(.div(z, .sub(S_One, z)))
      if (is_negone_num(s)) return(.div(z, .pow(.sub(S_One, z), S_Two)))
      if (isTRUE(ask_fact(s, "zero"))) return(.div(z, .sub(S_One, z)))
    }
    if (.has_function(z, c("exp_polar", "polar_lift")) && (isTRUE(zone) || .is_true(Le(fcall("Abs", list(z)), S_One))))
      return(fcall("polylog", list(s, unpolarify(z))))
    NULL
  },
  diff = function(args, i) if (i == 2L) .div(fcall("polylog", list(.sub(args[[1L]], S_One), args[[2L]])), args[[2L]]) else NULL,
  numeric = function(s, z) { k <- 1:2000; sum(Re(z)^k / k^Re(s)) })

def_function("lerchphi", nargs = 3L,
  numeric = function(z, s, a) { k <- 0:2000; sum(Re(z)^k / (k + Re(a))^Re(s)) })

# Expr.equals(1) for polylog.eval: FALSE for a nonconstant expression or a
# number away from 1, NA when undecided
.equals_one <- function(z) {
  if (is_one_num(z)) return(TRUE)
  if (length(free_symbols(z))) return(FALSE)
  v <- tryCatch(.neval(z, list()), error = function(err) NA)
  if (is.na(v)) return(NA)
  if (Mod(v - 1) > 1e-10) FALSE else NA
}

.dilogtable <- local({
  tab <- NULL
  function() {
    if (is.null(tab)) {
      s5 <- .sqrt(.int(5)); l2 <- fcall("log", list(S_Two)); pi2 <- .sq(S_Pi); C <- S_Catalan
      ent <- list(
        list(S_Half, .sub(.div(pi2, .int(12)), .div(.sq(l2), S_Two))),
        list(S_Two, .sub(.div(pi2, .int(4)), .mul(list(S_I, S_Pi, l2)))),
        list(.neg(.div(.sub(s5, S_One), S_Two)), .add(list(.neg(.div(pi2, .int(15))), .div(.sq(fcall("log", list(.div(.sub(s5, S_One), S_Two)))), S_Two)))),
        list(.neg(.div(.add(list(s5, S_One)), S_Two)), .sub(.neg(.div(pi2, .int(10))), .sq(fcall("log", list(.div(.add(list(s5, S_One)), S_Two)))))),
        list(.div(.sub(.int(3), s5), S_Two), .sub(.div(pi2, .int(15)), .sq(fcall("log", list(.div(.sub(s5, S_One), S_Two)))))),
        list(.div(.sub(s5, S_One), S_Two), .sub(.div(pi2, .int(10)), .sq(fcall("log", list(.div(.sub(s5, S_One), S_Two)))))),
        list(S_I, .sub(.mul(list(S_I, C)), .div(pi2, .int(48)))),
        list(.neg(S_I), .sub(.neg(.mul(list(S_I, C))), .div(pi2, .int(48)))),
        list(.sub(S_One, S_I), .add(list(.div(pi2, .int(16)), .neg(.mul(list(S_I, C))), .neg(.mul(list(S_Pi, S_I, l2, .rat(1, 4))))))),
        list(.add(list(S_One, S_I)), .add(list(.div(pi2, .int(16)), .mul(list(S_I, C)), .mul(list(S_Pi, S_I, l2, .rat(1, 4)))))),
        list(.div(.sub(S_One, S_I), S_Two), .add(list(.neg(.div(.sq(l2), .int(8))), .mul(list(S_Pi, S_I, l2, .rat(1, 8))),
                                                      .mul(list(.rat(5, 96), pi2)), .neg(.mul(list(S_I, C)))))))
      tab <<- stats::setNames(lapply(ent, `[[`, 2L), vapply(ent, function(e) e[[1L]]$key, ""))
    }
    tab
  }
})

# polylog._eval_expand_func
.expand_hooks$polylog <- function(e) {
  s <- e$args[[1L]]; z <- e$args[[2L]]
  if (is_one_num(s)) return(.neg(fcall("log", list(.sub(S_One, z)))))
  if (s$type == "Integer" && s$value <= 0) {
    u <- Dummy("u")
    start <- .div(u, .sub(S_One, u))
    for (k in seq_len(-.int_val(s))) start <- .mul(list(u, .diff1(start, u)))
    return(subs(expand_mul(start), u, z))
  }
  NULL
}

# lerchphi._eval_expand_func
.expand_hooks$lerchphi <- function(e) {
  z <- e$args[[1L]]; s <- e$args[[2L]]; a <- e$args[[3L]]
  if (is_one_num(z)) return(fcall("zeta", list(s, a)))
  if (s$type == "Integer" && s$value <= 0) {
    t <- Dummy("t")
    cs <- rev(.hx_coeffs(.pow(.add(list(t, a)), .neg(s)), t))
    start <- .inv(.sub(S_One, t))
    res <- S_Zero
    for (c in rev(cs)) {
      res <- .add(list(res, .mul(list(c, start))))
      start <- .mul(list(t, .diff1(start, t)))
    }
    return(subs(res, t, z))
  }
  if (a$type %in% c("Integer", "Rational")) {
    add <- S_Zero; mul <- S_One
    if (a$value > 1) {
      n <- floor(as.double(a$value))
      if (n == as.double(a$value)) n <- n - 1
      a <- .sub(a, .int(n))
      mul <- .pow(z, .int(-n))
      add <- .add(lapply(seq_len(n) - 1L, function(k) .neg(.div(.pow(z, .int(k - n)), .pow(.add(list(a, .int(k))), s)))))
    } else if (a$value <= 0) {
      n <- floor(-as.double(a$value)) + 1
      a <- .add(list(a, .int(n)))
      mul <- .pow(z, .int(n))
      add <- .add(lapply(seq_len(n) - 1L, function(k) .div(.pow(z, .int(n - 1L - k)), .pow(.sub(a, .int(k + 1L)), s))))
    }
    m <- .int(as.character(gmp::numerator(a$value))); n <- as.integer(gmp::denominator(a$value))
    zet <- fcall("exp_polar", list(.div(.mul(list(S_Two, S_Pi, S_I)), .int(n))))
    root <- .pow(z, .rat(1L, n))
    up_zet <- unpolarify(zet)
    addargs <- lapply(seq_len(n) - 1L, function(k) {
      p <- fcall("polylog", list(s, .mul(list(.pow(zet, .int(k)), root))))
      if (p$type == "Function" && p$name == "polylog") { q <- .expand_hooks$polylog(p); if (!is.null(q)) p <- q }
      .div(p, .pow(.mul(list(.pow(up_zet, .int(k)), root)), m))
    })
    return(.add(list(add, .mul(list(mul, .pow(.int(n), .sub(s, S_One)), .add(addargs))))))
  }
  zr <- NULL
  if (z$type == "Function" && z$name == "exp") {
    r <- .div(z$args[[1L]], .mul(list(S_Pi, S_I)))
    if (r$type %in% c("Integer", "Rational")) zr <- .div(z$args[[1L]], .mul(list(S_Two, S_Pi, S_I)))
  }
  if (!is.null(zr) || is_negone_num(z) || identical(z$key, S_I$key) || identical(z$key, .neg(S_I)$key)) {
    pq <- if (is_negone_num(z)) c(1L, 2L) else if (identical(z$key, S_I$key)) c(1L, 4L) else if (is.null(zr)) c(-1L, 4L)
          else c(as.integer(gmp::numerator(zr$value)), as.integer(gmp::denominator(zr$value)))
    p <- pq[1L]; q <- pq[2L]
    return(.add(lapply(seq_len(q) - 1L, function(k)
      .mul(list(fcall("exp", list(.mul(list(S_Two, S_Pi, S_I, .rat(k * p, q))))), .pow(.int(q), .neg(s)),
                fcall("zeta", list(s, .div(.add(list(.int(k), a)), .int(q)))))))))
  }
  NULL
}

# ---------------------------------------------------------------------------
# Bessel and Airy functions

# After sympy/functions/special/bessel.py: BesselBase with the constants
# _a and _b of fdiff and _eval_expand_func, the eval rules of besselj,
# bessely, besseli and besselk, and the spherical functions jn, yn, hn1, hn2.
.BES_TEX <- c(besselj = "J", bessely = "Y", besseli = "I", besselk = "K", hankel1 = "H^{(1)}", hankel2 = "H^{(2)}")
.BES_A <- c(besselj = 1, bessely = 1, besseli = -1, besselk = 1, hankel1 = 1, hankel2 = 1)
.BES_B <- c(besselj = 1, bessely = 1, besseli = 1, besselk = -1, hankel1 = 1, hankel2 = 1)
.fq <- function(z, nm) ask_fact(z, nm)
.re_f <- function(z, nm) ask_fact(fcall("re", list(z)), nm)

.bessel_eval <- function(name) function(args) {
  nu <- args[[1L]]; z <- args[[2L]]
  if (z$type == "Float" && is_number_type(nu)) return(.float_eval(name, args))
  if (name %in% c("hankel1", "hankel2")) return(NULL)
  bf <- function(n2, z2) fcall(name, list(n2, z2))
  if (isTRUE(.fq(z, "zero"))) {
    if (name %in% c("besselj", "besseli")) {
      if (isTRUE(.fq(nu, "zero"))) return(S_One)
      if ((isTRUE(.fq(nu, "integer")) && isFALSE(.fq(nu, "zero"))) || isTRUE(.re_f(nu, "positive"))) return(S_Zero)
      if (isTRUE(.re_f(nu, "negative")) && !isTRUE(.fq(nu, "integer"))) return(S_ComplexInfinity)
      if (isTRUE(.fq(nu, "imaginary"))) return(S_NaN)
    } else {
      if (isTRUE(.fq(nu, "zero"))) return(if (name == "bessely") S_NegativeInfinity else S_Infinity)
      rz <- .re_f(nu, "zero")
      if (isFALSE(rz)) return(S_ComplexInfinity)
      if (isTRUE(rz)) return(S_NaN)
    }
  }
  ioo <- .mul(list(S_I, S_Infinity)); mioo <- .mul(list(S_I, S_NegativeInfinity))
  is_k <- function(a, b) identical(a$key, b$key)
  if (name == "besselj") {
    if (is_const(z, "oo") || is_const(z, "-oo")) return(S_Zero)
    if (could_extract_minus_sign(z)) return(.mul(list(.pow(z, nu), .pow(.neg(z), .neg(nu)), bf(nu, .neg(z)))))
    if (isTRUE(.fq(nu, "integer"))) {
      if (could_extract_minus_sign(nu)) return(.mul(list(.pow(S_NegativeOne, .neg(nu)), bf(.neg(nu), z))))
      nz <- extract_multiplicatively(z, S_I)
      if (!is.null(nz) && !is_zero_num(nz)) return(.mul(list(.pow(S_I, nu), fcall("besseli", list(nu, nz)))))
    }
  } else if (name == "besseli") {
    imz <- fcall("im", list(z))
    if (is_const(imz, "oo") || is_const(imz, "-oo")) return(S_Zero)
    if (is_const(z, "oo")) return(S_Infinity)
    if (is_const(z, "-oo")) return(.mul(list(.pow(S_NegativeOne, nu), S_Infinity)))
    if (could_extract_minus_sign(z)) return(.mul(list(.pow(z, nu), .pow(.neg(z), .neg(nu)), bf(nu, .neg(z)))))
    if (isTRUE(.fq(nu, "integer"))) {
      if (could_extract_minus_sign(nu)) return(bf(.neg(nu), z))
      nz <- extract_multiplicatively(z, S_I)
      if (!is.null(nz) && !is_zero_num(nz)) return(.mul(list(.pow(S_I, .neg(nu)), fcall("besselj", list(nu, .neg(nz))))))
    }
  } else if (name == "bessely") {
    if (is_const(z, "oo") || is_const(z, "-oo")) return(S_Zero)
    if (is_k(z, ioo)) return(.mul(list(fcall("exp", list(.mul(list(S_I, S_Pi, .add(list(nu, S_One)), S_Half)))), S_Infinity)))
    if (is_k(z, mioo)) return(.mul(list(fcall("exp", list(.mul(list(.neg(S_I), S_Pi, .add(list(nu, S_One)), S_Half)))), S_Infinity)))
    if (isTRUE(.fq(nu, "integer")) && could_extract_minus_sign(nu)) return(.mul(list(.pow(S_NegativeOne, .neg(nu)), bf(.neg(nu), z))))
    return(NULL)
  } else if (name == "besselk") {
    if (is_const(z, "oo") || is_k(z, ioo) || is_k(z, mioo)) return(S_Zero)
    if (isTRUE(.fq(nu, "integer")) && could_extract_minus_sign(nu)) return(bf(.neg(nu), z))
    return(NULL)
  }
  # besselj and besseli: polar arguments and orders
  if (isTRUE(.fq(nu, "integer"))) {
    nz <- unpolarify(z)
    if (!is_k(nz, z)) return(bf(nu, nz))
  } else {
    eb <- extract_branch_factor(z)
    if (!is_zero_num(eb[[2L]])) return(.mul(list(fcall("exp", list(.mul(list(S_Two, eb[[2L]], S_Pi, nu, S_I)))), bf(nu, eb[[1L]]))))
  }
  nnu <- unpolarify(nu)
  if (!is_k(nnu, nu)) return(bf(nnu, z))
  NULL
}

.bessel_facts <- function(name) function(args, afs) {
  nu <- args[[1L]]; z <- args[[2L]]
  out <- c(commutative = TRUE)
  # _eval_is_extended_real (finiteness stays unknown)
  if (name %in% c("besselj", "besseli") && isTRUE(.fq(nu, "integer")) && isTRUE(.fq(z, "real"))) out["extended_real"] <- TRUE
  if (name %in% c("bessely", "besselk") && isTRUE(.fq(nu, "integer")) && isTRUE(.fq(z, "positive"))) out["extended_real"] <- TRUE
  out
}

# BesselBase._eval_expand_func: shift the order towards (-1, 1]
.bessel_expand <- function(name, nu, z) {
  a <- .BES_A[[name]]; b <- .BES_B[[name]]
  f <- function(n2) { v <- fcall(name, list(n2, z)); if (v$type == "Function" && v$name == name) .bessel_expand(name, v$args[[1L]], v$args[[2L]]) else v }
  if (isTRUE(.fq(nu, "real"))) {
    if (isTRUE(.fq(.sub(nu, S_One), "positive")))
      return(.add(list(.mul(list(.int(-a * b), f(.sub(nu, S_Two)))),
                       .div(.mul(list(.int(2 * a), .sub(nu, S_One), f(.sub(nu, S_One)))), z))))
    if (isTRUE(.fq(.add(list(nu, S_One)), "negative")))
      return(.sub(.div(.mul(list(.int(2 * b), .add(list(nu, S_One)), f(.add(list(nu, S_One))))), z),
                  .mul(list(.int(a * b), f(.add(list(nu, S_Two)))))))
  }
  fcall(name, list(nu, z))
}

.bessel_def <- function(name, rfun) {
  def_function(name, nargs = 2L,
    eval = .bessel_eval(name),
    diff = function(args, i) {
      if (i != 2L) return(NULL)
      nu <- args[[1L]]; z <- args[[2L]]
      .sub(.mul(list(.rat(.BES_B[[name]], 2), fcall(name, list(.sub(nu, S_One), z)))),
           .mul(list(.rat(.BES_A[[name]], 2), fcall(name, list(.add(list(nu, S_One)), z)))))
    },
    facts = .bessel_facts(name),
    conjugate = function(args) {
      # z.is_extended_negative is False
      z <- args[[2L]]
      if (!isFALSE(ask_fact(z, "negative"))) return(NULL)
      to <- switch(name, hankel1 = "hankel2", hankel2 = "hankel1", name)
      fcall(to, list(fcall("conjugate", list(args[[1L]])), fcall("conjugate", list(z))))
    },
    # LatexPrinter._hprint_BesselBase
    latex = function(a) paste0(.BES_TEX[[name]], "_{", a[1L], "}\\left(", a[2L], "\\right)"),
    numeric = rfun)
  .expand_hooks[[name]] <- function(e) {
    r <- .bessel_expand(name, e$args[[1L]], e$args[[2L]])
    if (identical(r$key, e$key)) NULL else r
  }
}
.bessel_def("besselj", function(nu, z) besselJ(Re(z), Re(nu)))
.bessel_def("bessely", function(nu, z) besselY(Re(z), Re(nu)))
.bessel_def("besseli", function(nu, z) besselI(Re(z), Re(nu)))
.bessel_def("besselk", function(nu, z) besselK(Re(z), Re(nu)))
.bessel_def("hankel1", function(nu, z) complex(real = besselJ(Re(z), Re(nu)), imaginary = besselY(Re(z), Re(nu))))
.bessel_def("hankel2", function(nu, z) complex(real = besselJ(Re(z), Re(nu)), imaginary = -besselY(Re(z), Re(nu))))

# spherical_bessel_fn(n, z): fn(0) = 1/z, fn(1) = 1/z^2,
# fn(n - 1) + fn(n + 1) = (2n + 1)/z fn(n), as a polynomial in 1/z
.sph_fn <- function(n, z) {
  # dup coefficient lists in t = 1/z, highest power first
  lshift <- function(p) c(p, gmp::as.bigz(0))
  sub <- function(p, q) { m <- max(length(p), length(q)); p <- c(gmp::as.bigz(rep(0, m - length(p))), p); q <- c(gmp::as.bigz(rep(0, m - length(q))), q); p - q }
  if (n >= 0) {
    if (n < 1) p <- gmp::as.bigz(c(1, 0)) else {
      m2 <- gmp::as.bigz(1); m1 <- gmp::as.bigz(c(1, 0))
      if (n >= 2) for (i in 2:n) { tmp <- m1; m1 <- sub(lshift(m1) * (2 * i - 1), m2); m2 <- tmp }
      p <- lshift(m1)
    }
  } else {
    k <- -n
    m2 <- gmp::as.bigz(c(1, 0)); m1 <- gmp::as.bigz(0)
    if (k >= 2) for (i in 2:k) { tmp <- m1; m1 <- sub(lshift(m1) * (3 - 2 * i), m2); m2 <- tmp }
    p <- m1
  }
  t <- .inv(z); d <- length(p) - 1L
  .add(lapply(seq_along(p), function(i) .mul(list(.int(p[i]), .pow(t, .int(d - i + 1L))))))
}
.sph_jn <- function(n, z) .add(list(.mul(list(.sph_fn(n, z), fcall("sin", list(z)))),
                                    .mul(list(.pow(S_NegativeOne, .int(n + 1)), .sph_fn(-n - 1, z), fcall("cos", list(z))))))
.sph_yn <- function(n, z) .sub(.mul(list(.pow(S_NegativeOne, .int(n + 1)), .sph_fn(-n - 1, z), fcall("sin", list(z)))),
                               .mul(list(.sph_fn(n, z), fcall("cos", list(z)))))

.sph_def <- function(name, numeric, zero_eval = FALSE) {
  def_function(name, nargs = 2L,
    eval = function(args) {
      nu <- args[[1L]]; z <- args[[2L]]
      if (z$type == "Float" && is_number_type(nu)) return(.float_eval(name, args))
      if (name == "jn") {
        if (isTRUE(.fq(z, "zero"))) {
          if (isTRUE(.fq(nu, "zero"))) return(S_One)
          if (isTRUE(.fq(nu, "integer"))) return(if (isTRUE(.fq(nu, "positive"))) S_Zero else S_ComplexInfinity)
        }
        if (is_const(z, "oo") || is_const(z, "-oo")) return(S_Zero)
      }
      NULL
    },
    diff = function(args, i) {
      if (i != 2L) return(NULL)
      nu <- args[[1L]]; z <- args[[2L]]
      .sub(fcall(name, list(.sub(nu, S_One), z)), .div(.mul(list(fcall(name, list(nu, z)), .add(list(nu, S_One)))), z))
    },
    numeric = numeric)
  .expand_hooks[[name]] <- function(e) {
    nu <- e$args[[1L]]; z <- e$args[[2L]]
    if (name %in% c("hn1", "hn2")) {
      hks <- if (name == "hn1") S_I else .neg(S_I)
      if (nu$type == "Integer") { n <- as.integer(nu$value); return(expand(.add(list(.sph_jn(n, z), .mul(list(hks, .sph_yn(n, z))))))) }
      return(.add(list(fcall("jn", list(nu, z)), .mul(list(hks, fcall("yn", list(nu, z)))))))
    }
    if (nu$type != "Integer") return(NULL)
    if (name == "jn") .sph_jn(as.integer(nu$value), z) else .sph_yn(as.integer(nu$value), z)
  }
}
.sph_def("jn", function(n, z) sqrt(pi / (2 * Re(z))) * besselJ(Re(z), Re(n) + 0.5))
.sph_def("yn", function(n, z) sqrt(pi / (2 * Re(z))) * besselY(Re(z), Re(n) + 0.5))
.sph_def("hn1", .sph_hankel(1))
.sph_def("hn2", .sph_hankel(-1))

.airy_num <- function(z, kind) {
  z <- Re(z)
  if (kind == "ai") {
    if (z > 0) return(sqrt(z / 3) / pi * besselK(2 / 3 * z^1.5, 1 / 3))
    if (z == 0) return(1 / (3^(2 / 3) * gamma(2 / 3)))
    w <- 2 / 3 * (-z)^1.5
    return(sqrt(-z) / 3 * (besselJ(w, 1 / 3) + besselJ(w, -1 / 3)))
  }
  if (z > 0) return(sqrt(z / 3) * (besselI(2 / 3 * z^1.5, -1 / 3) + besselI(2 / 3 * z^1.5, 1 / 3)))
  if (z == 0) return(1 / (3^(1 / 6) * gamma(2 / 3)))
  w <- 2 / 3 * (-z)^1.5
  sqrt(-z / 3) * (besselJ(w, -1 / 3) - besselJ(w, 1 / 3))
}

def_function("airyai",
  eval = function(args) if (is_zero_num(args[[1L]])) .inv(.mul(list(.pow(.int(3), .rat(2, 3)), fcall("gamma", list(.rat(2, 3)))))) else NULL,
  diff = function(args, i) fcall("airyaiprime", args), numeric = function(z) .airy_num(z, "ai"))
def_function("airybi",
  eval = function(args) if (is_zero_num(args[[1L]])) .inv(.mul(list(.pow(.int(3), .rat(1, 6)), fcall("gamma", list(.rat(2, 3)))))) else NULL,
  diff = function(args, i) fcall("airybiprime", args), numeric = function(z) .airy_num(z, "bi"))
def_function("airyaiprime",
  # airyaiprime.eval
  eval = function(args) {
    a <- args[[1L]]
    if (is_nan(a)) return(S_NaN)
    if (is_const(a, "oo")) return(S_Zero)
    if (isTRUE(ask_fact(a, "zero"))) return(.div(S_NegativeOne, .mul(list(.pow(.int(3), .rat(1, 3)), fcall("gamma", list(.rat(1, 3)))))))
    NULL
  },
  diff = function(args, i) .mul(list(args[[1L]], fcall("airyai", args))),
  numeric = function(z) (.airy_num(Re(z) + 1e-6, "ai") - .airy_num(Re(z) - 1e-6, "ai")) / 2e-6)
def_function("airybiprime",
  # airybiprime.eval
  eval = function(args) {
    a <- args[[1L]]
    if (is_nan(a)) return(S_NaN)
    if (is_const(a, "oo")) return(S_Infinity)
    if (is_const(a, "-oo")) return(S_Zero)
    if (isTRUE(ask_fact(a, "zero"))) return(.div(.pow(.int(3), .rat(1, 6)), fcall("gamma", list(.rat(1, 3)))))
    NULL
  },
  diff = function(args, i) .mul(list(args[[1L]], fcall("airybi", args))),
  numeric = function(z) (.airy_num(Re(z) + 1e-6, "bi") - .airy_num(Re(z) - 1e-6, "bi")) / 2e-6)

# ---------------------------------------------------------------------------
# Elliptic integrals

# m is one of the given constants, as SymPy's `m in (S.Infinity, ...)`
.ell_in <- function(m, vals) m$key %in% vapply(vals, .k, "")
.ell_infs <- function() list(S_Infinity, S_NegativeInfinity)

def_function("elliptic_k",
  # elliptic_k.eval
  eval = function(args) {
    m <- args[[1L]]
    if (isTRUE(ask_fact(m, "zero"))) return(.mul(list(S_Pi, S_Half)))
    if (identical(m$key, S_Half$key))
      return(.div(.mul(list(.int(8), .pow(S_Pi, .rat(3, 2)))), .sq(fcall("gamma", list(.rat(-1, 4))))))
    if (is_one_num(m)) return(S_ComplexInfinity)
    if (identical(m$key, S_NegativeOne$key))
      return(.div(.sq(fcall("gamma", list(.rat(1, 4)))), .mul(list(.int(4), .sqrt(.mul(list(S_Two, S_Pi)))))))
    if (.ell_in(m, list(S_Infinity, S_NegativeInfinity, .mul(list(S_I, S_Infinity)),
                        .mul(list(S_I, S_NegativeInfinity)), S_ComplexInfinity))) return(S_Zero)
    NULL
  },
  diff = function(args, i) {
    m <- args[[1L]]
    .div(.sub(fcall("elliptic_e", args), .mul(list(.sub(S_One, m), fcall("elliptic_k", args)))),
         .mul(list(S_Two, m, .sub(S_One, m))))
  },
  numeric = function(m) { a <- 1; b <- sqrt(1 - Re(m)); for (i in 1:40) { t <- (a + b) / 2; b <- sqrt(a * b); a <- t }; pi / (2 * a) })

def_function("elliptic_e", nargs = 1:2,
  # elliptic_e.eval
  eval = function(args) {
    if (length(args) == 2L) {
      z <- args[[1L]]; m <- args[[2L]]
      k <- .div(.mul(list(S_Two, z)), S_Pi)
      if (isTRUE(ask_fact(m, "zero"))) return(z)
      if (isTRUE(ask_fact(z, "zero"))) return(S_Zero)
      if (isTRUE(ask_fact(k, "integer"))) return(.mul(list(k, fcall("elliptic_e", list(m)))))
      if (.ell_in(m, .ell_infs())) return(S_ComplexInfinity)
      if (could_extract_minus_sign(z)) return(.neg(fcall("elliptic_e", list(.neg(z), m))))
      return(NULL)
    }
    m <- args[[1L]]
    if (isTRUE(ask_fact(m, "zero"))) return(.div(S_Pi, S_Two))
    if (is_one_num(m)) return(S_One)
    if (is_const(m, "oo")) return(.mul(list(S_I, S_Infinity)))
    if (is_const(m, "-oo")) return(S_Infinity)
    if (is_const(m, "zoo")) return(S_ComplexInfinity)
    NULL
  },
  diff = function(args, i) {
    if (length(args) == 1L) {
      m <- args[[1L]]
      return(.div(.sub(fcall("elliptic_e", args), fcall("elliptic_k", args)), .mul(list(S_Two, m))))
    }
    NULL
  },
  numeric = function(a, m) {
    # elliptic_e(m) is complete; elliptic_e(phi, m) runs to the amplitude phi
    if (missing(m)) return(stats::integrate(function(t) sqrt(1 - Re(a) * sin(t)^2), 0, pi / 2)$value)
    stats::integrate(function(t) sqrt(1 - Re(m) * sin(t)^2), 0, Re(a), rel.tol = 1e-12)$value
  })

def_function("elliptic_f", nargs = 2L,
  # elliptic_f.eval
  eval = function(args) {
    z <- args[[1L]]; m <- args[[2L]]
    if (isTRUE(ask_fact(z, "zero"))) return(S_Zero)
    if (isTRUE(ask_fact(m, "zero"))) return(z)
    k <- .div(.mul(list(S_Two, z)), S_Pi)
    if (isTRUE(ask_fact(k, "integer"))) return(.mul(list(k, fcall("elliptic_k", list(m)))))
    if (.ell_in(m, .ell_infs())) return(S_Zero)
    if (could_extract_minus_sign(z)) return(.neg(fcall("elliptic_f", list(.neg(z), m))))
    NULL
  },
  diff = function(args, i) if (i == 1L) .inv(.sqrt(.sub(S_One, .mul(list(args[[2L]], .sq(fcall("sin", list(args[[1L]])))))))) else NULL,
  numeric = function(phi, m) stats::integrate(function(t) 1 / sqrt(1 - Re(m) * sin(t)^2), 0, Re(phi))$value)

def_function("elliptic_pi", nargs = 2:3,
  # elliptic_pi.eval
  eval = function(args) {
    n <- args[[1L]]
    real_inf <- function(v) isTRUE(ask_fact(v, "extended_real")) && isTRUE(ask_fact(v, "infinite"))
    if (length(args) == 3L) {
      z <- args[[2L]]; m <- args[[3L]]
      if (isTRUE(ask_fact(n, "zero"))) return(fcall("elliptic_f", list(z, m)))
      if (is_one_num(n))
        return(.add(list(fcall("elliptic_f", list(z, m)),
                         .div(.sub(.mul(list(.sqrt(.sub(S_One, .mul(list(m, .sq(fcall("sin", list(z))))))), fcall("tan", list(z)))),
                                   fcall("elliptic_e", list(z, m))), .sub(S_One, m)))))
      k <- .div(.mul(list(S_Two, z)), S_Pi)
      if (isTRUE(ask_fact(k, "integer"))) return(.mul(list(k, fcall("elliptic_pi", list(n, m)))))
      if (isTRUE(ask_fact(m, "zero")))
        return(.div(fcall("atanh", list(.mul(list(.sqrt(.sub(n, S_One)), fcall("tan", list(z)))))), .sqrt(.sub(n, S_One))))
      if (identical(n$key, m$key))
        return(.add(list(fcall("elliptic_f", list(z, n)), .neg(fcall("elliptic_pi", list(S_One, z, n))),
                         .div(fcall("tan", list(z)), .sqrt(.sub(S_One, .mul(list(n, .sq(fcall("sin", list(z)))))))))))
      if (.ell_in(n, .ell_infs()) || .ell_in(m, .ell_infs())) return(S_Zero)
      if (could_extract_minus_sign(z)) return(.neg(fcall("elliptic_pi", list(n, .neg(z), m))))
      if (real_inf(m) || real_inf(n)) return(S_Zero)
      return(NULL)
    }
    m <- args[[2L]]
    if (isTRUE(ask_fact(n, "zero"))) return(fcall("elliptic_k", list(m)))
    if (is_one_num(n)) return(S_ComplexInfinity)
    if (isTRUE(ask_fact(m, "zero"))) return(.div(S_Pi, .mul(list(S_Two, .sqrt(.sub(S_One, n))))))
    if (is_one_num(m)) return(.div(S_NegativeInfinity, fcall("sign", list(.sub(n, S_One)))))
    if (identical(n$key, m$key)) return(.div(fcall("elliptic_e", list(n)), .sub(S_One, n)))
    if (.ell_in(n, .ell_infs()) || .ell_in(m, .ell_infs())) return(S_Zero)
    if (real_inf(m) || real_inf(n)) return(S_Zero)
    NULL
  },
  numeric = function(n, a, m) {
    # elliptic_pi(n, m) is complete; elliptic_pi(n, phi, m) runs to phi
    phi <- if (missing(m)) pi / 2 else Re(a)
    if (missing(m)) m <- a
    stats::integrate(function(t) 1 / ((1 - Re(n) * sin(t)^2) * sqrt(1 - Re(m) * sin(t)^2)), 0, phi, rel.tol = 1e-12)$value
  })

# ---------------------------------------------------------------------------
# Orthogonal polynomials, after sympy/functions/special/polynomials.py and
# sympy/polys/orthopolys.py. Integer degrees build the dense coefficient list
# of SymPy's dup_* recurrence over the smallest field holding the parameters
# (elements kept canonical by cancel and expand), then Poly.as_expr().

.op_K <- function(e) if (is_exact_number(e)) e else expand(cancel(e))
.op_add <- function(f, g) {
  n <- max(length(f), length(g))
  f <- c(rep(list(S_Zero), n - length(f)), f); g <- c(rep(list(S_Zero), n - length(g)), g)
  Map(function(u, v) .op_K(.add(list(u, v))), f, g)
}
.op_sub <- function(f, g) .op_add(f, lapply(g, .neg))
.op_mulc <- function(f, c) lapply(f, function(u) .op_K(.mul(list(u, c))))
.op_lshift <- function(f) c(f, list(S_Zero))
.op_mul <- function(f, g) {
  out <- rep(list(S_Zero), length(f) + length(g) - 1L)
  for (i in seq_along(f)) for (j in seq_along(g))
    out[[i + j - 1L]] <- .add(list(out[[i + j - 1L]], .mul(list(f[[i]], g[[j]]))))
  lapply(out, .op_K)
}
# Poly.new(dup, x).as_expr()
.op_expr <- function(f, x) {
  d <- length(f)
  terms <- list()
  for (i in seq_len(d)) if (!is_zero_num(f[[i]]))
    terms[[length(terms) + 1L]] <- .mul(list(f[[i]], .pow(x, .int(d - i))))
  .add(terms)
}
.op_rec2 <- function(n, m2, m1, step) {
  if (n < 1L) return(list(S_One))
  for (i in seq_len(n)[-1L]) { nx <- step(i, m2, m1); m2 <- m1; m1 <- nx }
  m1
}
.dup_jacobi <- function(n, a, b) {
  ab <- .add(list(a, b))
  .op_rec2(n, list(S_One), list(.op_K(.add(list(.div(ab, S_Two), S_One))), .op_K(.div(.sub(a, b), S_Two))),
    function(i, m2, m1) {
      I <- .int(i)
      den <- .op_K(.mul(list(I, .add(list(ab, I)), .add(list(ab, .int(2 * i - 2))))))
      f0 <- .op_K(.div(.mul(list(.add(list(ab, .int(2 * i - 1))), .sub(.mul(list(a, a)), .mul(list(b, b))))), .mul(list(S_Two, den))))
      f1 <- .op_K(.div(.mul(list(.add(list(ab, .int(2 * i - 1))), .add(list(ab, .int(2 * i - 2))), .add(list(ab, .int(2 * i))))), .mul(list(S_Two, den))))
      f2 <- .op_K(.div(.mul(list(.add(list(a, .int(i - 1))), .add(list(b, .int(i - 1))), .add(list(ab, .int(2 * i))))), den))
      .op_sub(.op_add(.op_mulc(m1, f0), .op_mulc(.op_lshift(m1), f1)), .op_mulc(m2, f2))
    })
}
.dup_gegenbauer <- function(n, a) {
  .op_rec2(n, list(S_One), list(.op_K(.mul(list(S_Two, a))), S_Zero), function(i, m2, m1) {
    q <- .op_K(.div(.mul(list(S_Two, .sub(a, S_One))), .int(i)))
    .op_sub(.op_mulc(.op_lshift(m1), .op_K(.add(list(q, S_Two)))), .op_mulc(m2, .op_K(.add(list(q, S_One)))))
  })
}
.dup_laguerre <- function(n, alpha) {
  m2 <- list(S_Zero); m1 <- list(S_One)
  for (i in seq_len(n)) {
    q <- .op_K(.div(.sub(alpha, S_One), .int(i)))
    a <- .op_mul(m1, list(.rat(-1, i), .op_K(.add(list(q, S_Two)))))
    b <- .op_mulc(m2, .op_K(.add(list(q, S_One))))
    m2 <- m1; m1 <- .op_sub(a, b)
  }
  m1
}
.dup_chebyshevt <- function(n) .op_rec2(n, list(S_One), list(S_One, S_Zero),
  function(i, m2, m1) .op_sub(.op_mulc(.op_lshift(m1), S_Two), m2))
.dup_chebyshevu <- function(n) .op_rec2(n, list(S_One), list(S_Two, S_Zero),
  function(i, m2, m1) .op_sub(.op_mulc(.op_lshift(m1), S_Two), m2))
.dup_hermite <- function(n) .op_rec2(n, list(S_One), list(S_Two, S_Zero),
  function(i, m2, m1) .op_mulc(.op_sub(.op_lshift(m1), .op_mulc(m2, .int(i - 1))), S_Two))
.dup_hermite_prob <- function(n) .op_rec2(n, list(S_One), list(S_One, S_Zero),
  function(i, m2, m1) .op_sub(.op_lshift(m1), .op_mulc(m2, .int(i - 1))))
.dup_legendre <- function(n) .op_rec2(n, list(S_One), list(S_One, S_Zero),
  function(i, m2, m1) .op_sub(.op_mulc(.op_lshift(m1), .rat(2 * i - 1, i)), .op_mulc(m2, .rat(i - 1, i))))

# expr.is_Number, which holds for the infinities and nan too
.is_Number <- function(e) is_number_type(e) || (e$type == "Constant" && e$name %in% c("oo", "-oo", "zoo", "nan"))
# OrthogonalPolynomial._eval_at_order: the polynomial for a nonnegative integer n
.op_at_order <- function(n, x, gen) {
  if (n$type != "Integer" || n$value < 0) return(NULL)
  .op_expr(gen(.int_val(n)), x)
}
.op_neg1n <- function(n) .pow(S_NegativeOne, n)
.op_value_error <- function(msg) stop(structure(class = c("ValueError", "error", "condition"),
                                               list(message = msg, call = NULL)))

def_function("jacobi", nargs = 4L,
  eval = function(args) {
    n <- args[[1L]]; a <- args[[2L]]; b <- args[[3L]]; x <- args[[4L]]
    if (identical(a$key, b$key)) {
      if (identical(a$key, .rat(-1, 2)$key))
        return(.mul(list(.div(fcall("RisingFactorial", list(S_Half, n)), fcall("factorial", list(n))), fcall("chebyshevt", list(n, x)))))
      if (isTRUE(ask_fact(a, "zero"))) return(fcall("legendre", list(n, x)))
      if (identical(a$key, S_Half$key))
        return(.mul(list(.div(fcall("RisingFactorial", list(.rat(3, 2), n)), fcall("factorial", list(.add(list(n, S_One))))), fcall("chebyshevu", list(n, x)))))
      return(.mul(list(.div(fcall("RisingFactorial", list(.add(list(a, S_One)), n)), fcall("RisingFactorial", list(.add(list(.mul(list(S_Two, a)), S_One)), n))),
                       fcall("gegenbauer", list(n, .add(list(a, S_Half)), x)))))
    }
    if (identical(b$key, .neg(a)$key))
      return(.mul(list(.div(fcall("gamma", list(.add(list(n, a, S_One)))), fcall("gamma", list(.add(list(n, S_One))))),
                       .div(.pow(.add(list(S_One, x)), .div(a, S_Two)), .pow(.sub(S_One, x), .div(a, S_Two))),
                       fcall("assoc_legendre", list(n, .neg(a), x)))))
    if (!.is_Number(n)) {
      if (could_extract_minus_sign(x)) return(.mul(list(.op_neg1n(n), fcall("jacobi", list(n, b, a, .neg(x))))))
      if (isTRUE(ask_fact(x, "zero")))
        return(.mul(list(.pow(S_Two, .neg(n)), fcall("gamma", list(.add(list(a, n, S_One)))),
                         .inv(.mul(list(fcall("gamma", list(.add(list(a, S_One)))), fcall("factorial", list(n))))),
                         fcall("hyper", list(fcall("Tuple", list(.sub(.neg(b), n), .neg(n))), fcall("Tuple", list(.add(list(a, S_One)))), S_NegativeOne)))))
      if (is_one_num(x)) return(.div(fcall("RisingFactorial", list(.add(list(a, S_One)), n)), fcall("factorial", list(n))))
      if (is_const(x, "oo") && isTRUE(ask_fact(n, "positive"))) {
        if (isTRUE(ask_fact(.add(list(a, b, .mul(list(S_Two, n)))), "integer")))
          .op_value_error("Error. a + b + 2*n should not be an integer.")
        return(.mul(list(fcall("RisingFactorial", list(.add(list(a, b, n, S_One)), n)), S_Infinity)))
      }
      return(NULL)
    }
    # named_poly builds the polynomial of index int(n), for a Rational n too
    if (!(n$type %in% c("Integer", "Rational"))) return(NULL)
    if (n$value < 0) .op_value_error(paste0("Cannot generate Jacobi polynomial of index ", .str(n)))
    .op_expr(.dup_jacobi(as.integer(trunc(as.numeric(n$value))), a, b), x)
  })

def_function("gegenbauer", nargs = 3L,
  eval = function(args) {
    n <- args[[1L]]; a <- args[[2L]]; x <- args[[3L]]
    if (isTRUE(ask_fact(n, "negative"))) return(S_Zero)
    if (identical(a$key, S_Half$key)) return(fcall("legendre", list(n, x)))
    if (is_one_num(a)) return(fcall("chebyshevu", list(n, x)))
    if (identical(a$key, S_NegativeOne$key)) return(S_Zero)
    if (!.is_Number(n)) {
      if (identical(x$key, S_NegativeOne$key)) {
        if (identical(Gt(fcall("re", list(a)), S_Half)$key, S_true$key)) return(S_ComplexInfinity)
        return(.div(.mul(list(fcall("cos", list(.mul(list(S_Pi, .add(list(a, n)))))), fcall("sec", list(.mul(list(S_Pi, a)))),
                              fcall("gamma", list(.add(list(.mul(list(S_Two, a)), n)))))),
                    .mul(list(fcall("gamma", list(.mul(list(S_Two, a)))), fcall("gamma", list(.add(list(n, S_One))))))))
      }
      if (could_extract_minus_sign(x)) return(.mul(list(.op_neg1n(n), fcall("gegenbauer", list(n, a, .neg(x))))))
      if (isTRUE(ask_fact(x, "zero"))) {
        num <- .mul(list(.pow(S_Two, n), .sqrt(S_Pi), fcall("gamma", list(.add(list(a, .mul(list(S_Half, n))))))))
        den <- .mul(list(fcall("gamma", list(.div(.sub(S_One, n), S_Two))), fcall("gamma", list(.add(list(n, S_One)))), fcall("gamma", list(a))))
        return(.div(num, den))
      }
      if (is_one_num(x))
        return(.div(fcall("gamma", list(.add(list(.mul(list(S_Two, a)), n)))),
                    .mul(list(fcall("gamma", list(.mul(list(S_Two, a)))), fcall("gamma", list(.add(list(n, S_One))))))))
      if (is_const(x, "oo") && isTRUE(ask_fact(n, "positive")))
        return(.mul(list(fcall("RisingFactorial", list(a, n)), S_Infinity)))
      return(NULL)
    }
    if (!(n$type %in% c("Integer", "Rational"))) return(NULL)
    .op_expr(.dup_gegenbauer(as.integer(trunc(as.numeric(n$value))), a), x)
  })

def_function("chebyshevt", nargs = 2L,
  eval = function(args) {
    n <- args[[1L]]; x <- args[[2L]]
    if (!.is_Number(n)) {
      if (could_extract_minus_sign(x)) return(.mul(list(.op_neg1n(n), fcall("chebyshevt", list(n, .neg(x))))))
      if (could_extract_minus_sign(n)) return(fcall("chebyshevt", list(.neg(n), x)))
      if (isTRUE(ask_fact(x, "zero"))) return(fcall("cos", list(.mul(list(S_Half, S_Pi, n)))))
      if (is_one_num(x)) return(S_One)
      if (is_const(x, "oo")) return(S_Infinity)
      return(NULL)
    }
    if (isTRUE(ask_fact(n, "negative"))) n <- .neg(n)
    .op_at_order(n, x, .dup_chebyshevt)
  },
  diff = function(args, i) if (i == 2L) .mul(list(args[[1L]], fcall("chebyshevu", list(.sub(args[[1L]], S_One), args[[2L]])))) else NULL)

def_function("chebyshevu", nargs = 2L,
  eval = function(args) {
    n <- args[[1L]]; x <- args[[2L]]
    if (!.is_Number(n)) {
      if (could_extract_minus_sign(x)) return(.mul(list(.op_neg1n(n), fcall("chebyshevu", list(n, .neg(x))))))
      if (could_extract_minus_sign(n)) {
        if (identical(n$key, S_NegativeOne$key)) return(S_Zero)
        if (!could_extract_minus_sign(.sub(.neg(n), S_Two))) return(.neg(fcall("chebyshevu", list(.sub(.neg(n), S_Two), x))))
      }
      if (isTRUE(ask_fact(x, "zero"))) return(fcall("cos", list(.mul(list(S_Half, S_Pi, n)))))
      if (is_one_num(x)) return(.add(list(S_One, n)))
      if (is_const(x, "oo")) return(S_Infinity)
      return(NULL)
    }
    if (isTRUE(ask_fact(n, "negative"))) {
      if (identical(n$key, S_NegativeOne$key)) return(S_Zero)
      r <- .op_at_order(.sub(.neg(n), S_Two), x, .dup_chebyshevu)
      return(if (is.null(r)) NULL else .neg(r))
    }
    .op_at_order(n, x, .dup_chebyshevu)
  })

def_function("legendre", nargs = 2L,
  eval = function(args) {
    n <- args[[1L]]; x <- args[[2L]]
    if (!.is_Number(n)) {
      if (could_extract_minus_sign(x)) return(.mul(list(.op_neg1n(n), fcall("legendre", list(n, .neg(x))))))
      if (could_extract_minus_sign(n) && !could_extract_minus_sign(.sub(.neg(n), S_One)))
        return(fcall("legendre", list(.sub(.neg(n), S_One), x)))
      if (isTRUE(ask_fact(x, "zero")))
        return(.div(.sqrt(S_Pi), .mul(list(fcall("gamma", list(.sub(S_Half, .div(n, S_Two)))), fcall("gamma", list(.add(list(S_One, .div(n, S_Two)))))))))
      if (is_one_num(x)) return(S_One)
      if (is_const(x, "oo")) return(S_Infinity)
      return(NULL)
    }
    if (isTRUE(ask_fact(n, "negative"))) n <- .sub(.neg(n), S_One)
    .op_at_order(n, x, .dup_legendre)
  },
  diff = function(args, i) {
    if (i != 2L) return(NULL)
    n <- args[[1L]]; x <- args[[2L]]
    .div(.mul(list(n, .sub(.mul(list(x, fcall("legendre", args))), fcall("legendre", list(.sub(n, S_One), x))))), .sub(.sq(x), S_One))
  })

def_function("hermite", nargs = 2L,
  eval = function(args) {
    n <- args[[1L]]; x <- args[[2L]]
    if (!.is_Number(n)) {
      if (could_extract_minus_sign(x)) return(.mul(list(.op_neg1n(n), fcall("hermite", list(n, .neg(x))))))
      if (isTRUE(ask_fact(x, "zero"))) return(.div(.mul(list(.pow(S_Two, n), .sqrt(S_Pi))), fcall("gamma", list(.div(.sub(S_One, n), S_Two)))))
      if (is_const(x, "oo")) return(S_Infinity)
      return(NULL)
    }
    if (isTRUE(ask_fact(n, "negative"))) .op_value_error(paste0("The index n must be nonnegative integer (got ", .str(n), ")"))
    .op_at_order(n, x, .dup_hermite)
  },
  diff = function(args, i) if (i == 2L) .mul(list(S_Two, args[[1L]], fcall("hermite", list(.sub(args[[1L]], S_One), args[[2L]])))) else NULL)

def_function("hermite_prob", nargs = 2L,
  eval = function(args) {
    n <- args[[1L]]; x <- args[[2L]]
    if (!.is_Number(n)) {
      if (could_extract_minus_sign(x)) return(.mul(list(.op_neg1n(n), fcall("hermite_prob", list(n, .neg(x))))))
      if (isTRUE(ask_fact(x, "zero"))) return(.div(.sqrt(S_Pi), fcall("gamma", list(.div(.sub(S_One, n), S_Two)))))
      if (is_const(x, "oo")) return(S_Infinity)
      return(NULL)
    }
    if (isTRUE(ask_fact(n, "negative"))) return(NULL)
    .op_at_order(n, x, .dup_hermite_prob)
  })

def_function("laguerre", nargs = 2L,
  eval = function(args) {
    n <- args[[1L]]; x <- args[[2L]]
    if (isFALSE(ask_fact(n, "integer"))) .op_value_error("Error: n should be an integer.")
    if (!.is_Number(n)) {
      if (could_extract_minus_sign(n) && !could_extract_minus_sign(.sub(.neg(n), S_One)))
        return(.mul(list(fcall("exp", list(x)), fcall("laguerre", list(.sub(.neg(n), S_One), .neg(x))))))
      if (isTRUE(ask_fact(x, "zero"))) return(S_One)
      if (is_const(x, "-oo")) return(S_Infinity)
      if (is_const(x, "oo")) return(.mul(list(.op_neg1n(n), S_Infinity)))
      return(NULL)
    }
    if (isTRUE(ask_fact(n, "negative")))
      return(.mul(list(fcall("exp", list(x)), fcall("laguerre", list(.sub(.neg(n), S_One), .neg(x))))))
    .op_at_order(n, x, function(k) .dup_laguerre(k, S_Zero))
  },
  diff = function(args, i) if (i == 2L) .neg(fcall("assoc_laguerre", list(.sub(args[[1L]], S_One), S_One, args[[2L]]))) else NULL)

def_function("assoc_laguerre", nargs = 3L,
  eval = function(args) {
    n <- args[[1L]]; alpha <- args[[2L]]; x <- args[[3L]]
    if (isTRUE(ask_fact(alpha, "zero"))) return(fcall("laguerre", list(n, x)))
    if (!.is_Number(n)) {
      if (isTRUE(ask_fact(x, "zero"))) return(fcall("binomial", list(.add(list(n, alpha)), alpha)))
      if (is_const(x, "oo") && isTRUE(ask_fact(n, "positive"))) return(.mul(list(.op_neg1n(n), S_Infinity)))
      if (is_const(x, "-oo") && isTRUE(ask_fact(n, "positive"))) return(S_Infinity)
      return(NULL)
    }
    if (isTRUE(ask_fact(n, "negative"))) .op_value_error(paste0("The index n must be nonnegative integer (got ", .str(n), ")"))
    if (!(n$type %in% c("Integer", "Rational"))) return(NULL)
    .op_expr(.dup_laguerre(as.integer(trunc(as.numeric(n$value))), alpha), x)
  })

.op_diff <- function(f) {
  d <- length(f) - 1L
  if (d < 1L) return(list(S_Zero))
  lapply(seq_len(d), function(i) .op_K(.mul(list(f[[i]], .int(d - i + 1L)))))
}

def_function("assoc_legendre", nargs = 3L,
  eval = function(args) {
    n <- args[[1L]]; m <- args[[2L]]; x <- args[[3L]]
    if (could_extract_minus_sign(m))
      return(.mul(list(.pow(S_NegativeOne, .neg(m)),
                       .div(fcall("factorial", list(.add(list(m, n)))), fcall("factorial", list(.sub(n, m)))),
                       fcall("assoc_legendre", list(n, .neg(m), x)))))
    if (is_zero_num(m)) return(fcall("legendre", list(n, x)))
    if (is_zero_num(x))
      return(.div(.mul(list(.pow(S_Two, m), .sqrt(S_Pi))),
                  .mul(list(fcall("gamma", list(.div(.sub(.sub(S_One, m), n), S_Two))),
                            fcall("gamma", list(.sub(S_One, .div(.sub(m, n), S_Two))))))))
    if (.is_Number(n) && .is_Number(m) && isTRUE(ask_fact(n, "integer")) && isTRUE(ask_fact(m, "integer"))) {
      if (isTRUE(ask_fact(n, "negative")))
        .op_value_error(paste0("assoc_legendre : 1st index must be nonnegative integer (got ", .str(n), ")"))
      nv <- .int_val(n); mv <- abs(.int_val(m))
      if (mv > nv)
        .op_value_error(paste0("assoc_legendre : abs('2nd index') must be <= '1st index' (got ", .str(n), ", ", .str(m), ")"))
      P <- .dup_legendre(nv)
      for (k in seq_len(mv)) P <- .op_diff(P)
      return(.mul(list(.pow(S_NegativeOne, .int(mv)), .pow(.sub(S_One, .sq(x)), .rat(mv, 2)), .op_expr(P, x))))
    }
    NULL
  })

def_function("Ynm", nargs = 4L,
  # Ynm.eval: only the signs of m, theta and phi are normalized
  eval = function(args) {
    n <- args[[1L]]; m <- args[[2L]]; theta <- args[[3L]]; phi <- args[[4L]]
    e2 <- function(m, phi) fcall("exp", list(.mul(list(.int(-2), S_I, m, phi))))
    if (could_extract_minus_sign(m)) {
      m <- .neg(m)
      return(.mul(list(.pow(S_NegativeOne, m), e2(m, phi), fcall("Ynm", list(n, m, theta, phi)))))
    }
    if (could_extract_minus_sign(theta)) return(fcall("Ynm", list(n, m, .neg(theta), phi)))
    if (could_extract_minus_sign(phi)) {
      phi <- .neg(phi)
      return(.mul(list(e2(m, phi), fcall("Ynm", list(n, m, theta, phi)))))
    }
    NULL
  },
  diff = function(args, i) {
    n <- args[[1L]]; m <- args[[2L]]; theta <- args[[3L]]; phi <- args[[4L]]
    if (i == 3L)
      return(.add(list(.mul(list(m, fcall("cot", list(theta)), fcall("Ynm", args))),
                       .mul(list(.sqrt(.mul(list(.sub(n, m), .add(list(n, m, S_One))))), fcall("exp", list(.mul(list(S_NegativeOne, S_I, phi)))),
                                 fcall("Ynm", list(n, .add(list(m, S_One)), theta, phi)))))))
    if (i == 4L) return(.mul(list(S_I, m, fcall("Ynm", args))))
    NULL
  },
  conjugate = function(args) .mul(list(.pow(S_NegativeOne, args[[2L]]), fcall("Ynm", list(args[[1L]], .neg(args[[2L]]), args[[3L]], args[[4L]])))))

# Ynm._eval_expand_func
.expand_hooks$Ynm <- function(e) {
  n <- e$args[[1L]]; m <- e$args[[2L]]; theta <- e$args[[3L]]; phi <- e$args[[4L]]
  rv <- .mul(list(.sqrt(.mul(list(.div(.add(list(.mul(list(S_Two, n)), S_One)), .mul(list(.int(4), S_Pi))),
                                  .div(fcall("factorial", list(.sub(n, m))), fcall("factorial", list(.add(list(n, m)))))))),
                  fcall("exp", list(.mul(list(S_I, m, phi)))), fcall("assoc_legendre", list(n, m, fcall("cos", list(theta))))))
  subs(rv, .sqrt(.add(list(.neg(.sq(fcall("cos", list(theta)))), S_One))), fcall("sin", list(theta)))
}

# Znm.eval
def_function("Znm", nargs = 4L, eval = function(args) {
  m <- args[[2L]]
  y <- fcall("Ynm", args)
  yc <- fcall("conjugate", list(fcall("Ynm", args)))
  if (isTRUE(ask_fact(m, "positive"))) return(.div(.add(list(y, yc)), .sqrt(S_Two)))
  if (isTRUE(ask_fact(m, "zero"))) return(y)
  if (isTRUE(ask_fact(m, "negative"))) return(.div(.sub(y, yc), .mul(list(.sqrt(S_Two), S_I))))
  NULL
})

# ---------------------------------------------------------------------------
# Hypergeometric and Meijer G (representation and numerics).

# hyper() and meijerg() are defined in functions-extra.R.

#' Special functions
#'
#' Symbolic special functions with special values, derivatives and numeric
#' evaluation: the gamma family (`gamma`, `loggamma`, `polygamma`,
#' `digamma`, `trigamma`, `beta`, `lowergamma`, `uppergamma`), factorials
#' (`factorial`, `factorial2`, `subfactorial`, `binomial`,
#' `RisingFactorial`, `FallingFactorial`), combinatorial numbers
#' (`bernoulli`, `fibonacci`, `lucas`, `tribonacci`, `catalan`, `bell`,
#' `harmonic`, `euler`, `genocchi`, `partition`, `motzkin`), error and
#' exponential integrals (`erf`, `erfc`, `erfi`, `erfinv`, `erfcinv`, `Ei`,
#' `li`, `Si`, `Ci`, `Shi`, `Chi`, `expint`, `fresnels`, `fresnelc`), zeta
#' functions (`zeta`, `dirichlet_eta`, `polylog`, `lerchphi`), Bessel and Airy
#' functions, elliptic integrals and orthogonal polynomials (`legendre`,
#' `assoc_legendre`, `chebyshevt`, `chebyshevu`, `hermite`, `hermite_prob`,
#' `laguerre`, `assoc_laguerre`, `gegenbauer`, `jacobi`, `Ynm`).
#'
#' `sym_factorial()` and `sym_beta()` dispatch to the base R functions for numeric
#' input.
#'
#' @param x,a,b,n,k Arguments.
#' @param k_sym,symbols For `bell()`: the second index of the incomplete Bell
#'   polynomial and the list of its symbols.
#' @param ... Arguments.
#' @return A `sym` expression.
#' @examples
#' sym_factorial(S(20))
#' sym_binomial(Symbol("n"), 2)
#' zeta(S(2))
#' legendre(3, Symbol("x"))
#' @name special-functions
NULL

.sf <- function(name) {
  force(name)
  function(...) fcall(name, lapply(list(...), as_sym))
}

#' @noRd
factorial <- function(x) if (is_sym(x)) fcall("factorial", list(x)) else base::factorial(x)
#' @noRd
beta <- function(a, b) if (is_sym(a) || is_sym(b)) fcall("beta", list(as_sym(a), as_sym(b))) else base::beta(a, b)
#' @rdname special-functions
#' @export
loggamma <- .sf("loggamma")
#' @rdname special-functions
#' @export
polygamma <- .sf("polygamma")
#' @rdname special-functions
#' @export
digamma_sym <- function(x) fcall("polygamma", list(S_Zero, as_sym(x)))
#' @rdname special-functions
#' @export
lowergamma <- .sf("lowergamma")
#' @rdname special-functions
#' @export
uppergamma <- .sf("uppergamma")
#' @rdname special-functions
#' @export
factorial2 <- .sf("factorial2")
#' @rdname special-functions
#' @export
subfactorial <- .sf("subfactorial")
#' @noRd
binomial <- function(n, k, ...) {
  if (missing(n) || missing(k) || ...length() || is.character(n)) return(.base_call(quote(stats::binomial)))
  fcall("binomial", list(as_sym(n), as_sym(k)))
}
#' @rdname special-functions
#' @export
RisingFactorial <- .sf("RisingFactorial")
#' @rdname special-functions
#' @export
FallingFactorial <- .sf("FallingFactorial")
#' @rdname special-functions
#' @export
bernoulli <- .sf("bernoulli")
#' @rdname special-functions
#' @export
fibonacci <- .sf("fibonacci")
#' @rdname special-functions
#' @export
lucas <- .sf("lucas")
#' @rdname special-functions
#' @export
tribonacci <- .sf("tribonacci")
#' @rdname special-functions
#' @export
catalan <- .sf("catalan")
#' @rdname special-functions
#' @export
bell <- function(n, k_sym = NULL, symbols = NULL) {
  args <- list(as_sym(n))
  if (!is.null(k_sym)) args[[2L]] <- as_sym(k_sym)
  if (!is.null(symbols)) args[[3L]] <- if (is_sym(symbols)) symbols else do.call(Tuple, as.list(symbols))
  fcall("bell", args)
}
#' @rdname special-functions
#' @export
harmonic <- .sf("harmonic")
#' @rdname special-functions
#' @export
euler <- .sf("euler")
#' @rdname special-functions
#' @export
genocchi <- .sf("genocchi")
#' @rdname special-functions
#' @export
partition <- .sf("partition")
#' @rdname special-functions
#' @export
motzkin <- .sf("motzkin")
#' @rdname special-functions
#' @export
erf <- .sf("erf")
#' @rdname special-functions
#' @export
erfc <- .sf("erfc")
#' @rdname special-functions
#' @export
erfi <- .sf("erfi")
#' @rdname special-functions
#' @export
erfinv <- .sf("erfinv")
#' @rdname special-functions
#' @export
erfcinv <- .sf("erfcinv")
#' @rdname special-functions
#' @export
Ei <- .sf("Ei")
#' @rdname special-functions
#' @export
li <- .sf("li")
#' @rdname special-functions
#' @export
Si <- .sf("Si")
#' @rdname special-functions
#' @export
Ci <- .sf("Ci")
#' @rdname special-functions
#' @export
Shi <- .sf("Shi")
#' @rdname special-functions
#' @export
Chi <- .sf("Chi")
#' @rdname special-functions
#' @export
expint <- .sf("expint")
#' @rdname special-functions
#' @export
fresnels <- .sf("fresnels")
#' @rdname special-functions
#' @export
fresnelc <- .sf("fresnelc")
#' @rdname special-functions
#' @export
zeta <- .sf("zeta")
#' @rdname special-functions
#' @export
dirichlet_eta <- .sf("dirichlet_eta")
#' @rdname special-functions
#' @export
polylog <- .sf("polylog")
#' @rdname special-functions
#' @export
lerchphi <- .sf("lerchphi")
#' @rdname special-functions
#' @export
besselj <- .sf("besselj")
#' @rdname special-functions
#' @export
bessely <- .sf("bessely")
#' @rdname special-functions
#' @export
besseli <- .sf("besseli")
#' @rdname special-functions
#' @export
besselk <- .sf("besselk")
#' @rdname special-functions
#' @export
hankel1 <- .sf("hankel1")
#' @rdname special-functions
#' @export
hankel2 <- .sf("hankel2")
#' @rdname special-functions
#' @export
jn <- .sf("jn")
#' @rdname special-functions
#' @export
yn <- .sf("yn")
#' @rdname special-functions
#' @export
airyai <- .sf("airyai")
#' @rdname special-functions
#' @export
airybi <- .sf("airybi")
#' @rdname special-functions
#' @export
airyaiprime <- .sf("airyaiprime")
#' @rdname special-functions
#' @export
airybiprime <- .sf("airybiprime")
#' @rdname special-functions
#' @export
elliptic_k <- .sf("elliptic_k")
#' @rdname special-functions
#' @export
elliptic_e <- .sf("elliptic_e")
#' @rdname special-functions
#' @export
elliptic_f <- .sf("elliptic_f")
#' @rdname special-functions
#' @export
elliptic_pi <- .sf("elliptic_pi")
#' @rdname special-functions
#' @export
legendre <- .sf("legendre")
#' @rdname special-functions
#' @export
assoc_legendre <- .sf("assoc_legendre")
#' @rdname special-functions
#' @export
chebyshevt <- .sf("chebyshevt")
#' @rdname special-functions
#' @export
chebyshevu <- .sf("chebyshevu")
#' @rdname special-functions
#' @export
hermite <- .sf("hermite")
#' @rdname special-functions
#' @export
hermite_prob <- .sf("hermite_prob")
#' @rdname special-functions
#' @export
laguerre <- .sf("laguerre")
#' @rdname special-functions
#' @export
assoc_laguerre <- .sf("assoc_laguerre")
#' @rdname special-functions
#' @export
gegenbauer <- .sf("gegenbauer")
#' @rdname special-functions
#' @export
jacobi <- .sf("jacobi")
#' @rdname special-functions
#' @export
Ynm <- .sf("Ynm")

# expand_func rules
.expand_hooks$gamma <- function(e) {
  a <- e$args[[1L]]
  if (a$type == "Add") {
    nums <- Filter(is_int_num, a$args)
    if (length(nums)) {
      n <- .int_val(nums[[1L]])
      rest <- .sub(a, nums[[1L]])
      if (n > 0 && n <= 50) return(.mul(c(lapply(0:(n - 1), function(k) .add(list(rest, .int(k)))), list(fcall("gamma", list(rest))))))
      if (n < 0 && n >= -50) return(.div(fcall("gamma", list(rest)), .mul(lapply(1:(-n), function(k) .sub(rest, .int(k))))))
    }
  }
  NULL
}
# factorial has no _eval_expand_func in SymPy; binomial._eval_expand_func
.expand_hooks$binomial <- function(e) {
  n <- e$args[[1L]]; k <- e$args[[2L]]
  if (.is_Number(n)) return(NULL)
  if (.sub(n, k)$type == "Integer") k <- .sub(n, k)
  if (k$type != "Integer") return(NULL)
  if (k$value == 0) return(S_One)
  if (k$value < 0) return(S_Zero)
  kv <- .int_val(k)
  res <- S_One
  for (i in seq_len(kv)) res <- .mul(list(res, .add(list(.sub(n, k), .int(i)))))
  .div(res, make_rational(gmp::as.bigq(gmp::factorialZ(kv))))
}
# harmonic._eval_expand_func
.expand_hooks$harmonic <- function(e) {
  n <- e$args[[1L]]
  if (length(e$args) == 2L && !is_one_num(e$args[[2L]])) return(NULL)
  if (n$type == "Add") {
    off <- n$args[[1L]]
    nnew <- .sub(n, off)
    if (off$type == "Integer" && off$value > 0)
      return(.add(c(lapply(seq(.int_val(off), 1L), function(i) .inv(.add(list(nnew, .int(i))))), list(fcall("harmonic", list(nnew))))))
    if (off$type == "Integer" && off$value < 0)
      return(.add(c(lapply(seq(0L, .int_val(off) + 1L), function(i) .neg(.inv(.add(list(nnew, .int(i)))))), list(fcall("harmonic", list(nnew))))))
  }
  if (n$type %in% c("Integer", "Rational")) {
    p <- gmp::numerator(gmp::as.bigq(n$value)); q <- gmp::denominator(gmp::as.bigq(n$value))
    u <- p %/% q; p <- p - u * q
    if (u >= 0 && p > 0 && q > 0 && p < q) {
      k <- Dummy("k"); P <- .zs(p); Q <- .zs(q)
      t1 <- .mul(list(Q, Sum(.inv(.add(list(.mul(list(Q, k)), P))), list(k, S_Zero, .zs(u)))))
      t2 <- .mul(list(S_Two, Sum(.mul(list(fcall("cos", list(.div(.mul(list(S_Two, S_Pi, P, k)), Q))),
                                           fcall("log", list(fcall("sin", list(.div(.mul(list(S_Pi, k)), Q))))))),
                                 list(k, S_One, fcall("floor", list(.div(.sub(Q, S_One), S_Two)))))))
      t3 <- .add(list(.mul(list(.div(S_Pi, S_Two), fcall("cot", list(.div(.mul(list(S_Pi, P)), Q))))), fcall("log", list(.mul(list(S_Two, Q))))))
      return(.sub(.add(list(t1, t2)), t3))
    }
  }
  NULL
}

# TrigonometricIntegral.eval (Si, Ci, Shi, Chi)
.TI <- list(
  Si = list(trig0 = 0, atzero = function() S_Zero, atinf = function() .mul(list(S_Pi, S_Half)), atneginf = function() .neg(.mul(list(S_Pi, S_Half))),
            minus = function(z) .neg(fcall("Si", list(z))), If = function(z, sg) .mul(list(S_I, fcall("Shi", list(z)), .int(sg)))),
  Ci = list(trig0 = 1, atzero = function() S_ComplexInfinity, atinf = function() S_Zero, atneginf = function() .mul(list(S_I, S_Pi)),
            minus = function(z) .add(list(fcall("Ci", list(z)), .mul(list(S_I, S_Pi)))),
            If = function(z, sg) .add(list(fcall("Chi", list(z)), .mul(list(S_I, S_Pi, S_Half, .int(sg)))))),
  Shi = list(trig0 = 0, atzero = function() S_Zero, atinf = function() S_Infinity, atneginf = function() S_NegativeInfinity,
             minus = function(z) .neg(fcall("Shi", list(z))), If = function(z, sg) .mul(list(S_I, fcall("Si", list(z)), .int(sg)))),
  Chi = list(trig0 = 1, atzero = function() S_ComplexInfinity, atinf = function() S_Infinity, atneginf = function() S_Infinity,
             minus = function(z) .add(list(fcall("Chi", list(z)), .mul(list(S_I, S_Pi)))),
             If = function(z, sg) .add(list(fcall("Ci", list(z)), .mul(list(S_I, S_Pi, S_Half, .int(sg)))))))

.ti_eval <- function(nm, z) {
  cl <- .TI[[nm]]
  if (is_zero_num(z)) return(cl$atzero())
  if (is_const(z, "oo")) return(cl$atinf())
  if (is_const(z, "-oo")) return(cl$atneginf())
  if (isTRUE(ask_fact(z, "zero"))) return(cl$atzero())
  nz <- extract_multiplicatively(z, fcall("polar_lift", list(S_I)))
  if (is.null(nz) && cl$trig0 == 0) nz <- extract_multiplicatively(z, S_I)
  if (!is.null(nz)) return(cl$If(nz, 1L))
  nz <- extract_multiplicatively(z, fcall("polar_lift", list(.neg(S_I))))
  if (!is.null(nz)) return(cl$If(nz, -1L))
  nz <- extract_multiplicatively(z, fcall("polar_lift", list(S_NegativeOne)))
  if (is.null(nz) && cl$trig0 == 0) nz <- extract_multiplicatively(z, S_NegativeOne)
  if (!is.null(nz)) return(cl$minus(nz))
  bf <- extract_branch_factor(z)
  nz <- bf[[1L]]; n <- bf[[2L]]
  if (is_zero_num(n) && identical(nz$key, z$key)) return(NULL)
  .add(list(.mul(list(S_Two, S_Pi, S_I, n, .int(cl$trig0))), fcall(nm, list(nz))))
}

# FresnelIntegral.eval
.fresnel_eval <- function(nm, z) {
  if (is_const(z, "oo")) return(S_Half)
  if (isTRUE(ask_fact(z, "zero"))) return(S_Zero)
  sgn <- if (nm == "fresnels") S_NegativeOne else S_One
  prefact <- S_One; newarg <- z; changed <- FALSE
  nz <- extract_multiplicatively(newarg, S_NegativeOne)
  if (!is.null(nz)) { prefact <- .neg(prefact); newarg <- nz; changed <- TRUE }
  nz <- extract_multiplicatively(newarg, S_I)
  if (!is.null(nz)) { prefact <- .mul(list(sgn, S_I, prefact)); newarg <- nz; changed <- TRUE }
  if (changed) return(.mul(list(prefact, fcall(nm, list(newarg)))))
  NULL
}

# Ei.eval
.Ei_eval <- function(z) {
  if (isTRUE(ask_fact(z, "zero"))) return(S_NegativeInfinity)
  if (is_const(z, "oo")) return(S_Infinity)
  if (is_const(z, "-oo")) return(S_Zero)
  bf <- extract_branch_factor(z)
  if (!is_zero_num(bf[[2L]])) return(.add(list(fcall("Ei", list(bf[[1L]])), .mul(list(S_Two, S_I, S_Pi, bf[[2L]])))))
  NULL
}

# expint.eval
.expint_eval <- function(nu, z) {
  nu2 <- unpolarify(nu)
  if (!identical(nu$key, nu2$key)) return(fcall("expint", list(nu2, z)))
  if ((nu$type == "Integer" && .int_val(nu) <= 0) || (nu$type != "Integer" && .mul(list(S_Two, nu))$type == "Integer"))
    return(unpolarify(expand_mul(.mul(list(.pow(z, .sub(nu, S_One)), fcall("uppergamma", list(.sub(S_One, nu), z)))))))
  bf <- extract_branch_factor(z)
  z <- bf[[1L]]; n <- bf[[2L]]
  if (is_zero_num(n)) return(NULL)
  if (isTRUE(ask_fact(nu, "integer"))) {
    if (!.is_true(Gt(nu, 0))) return(NULL)
    return(.sub(fcall("expint", list(nu, z)), .mul(list(S_Two, S_Pi, S_I, n, .pow(S_NegativeOne, .sub(nu, S_One)),
                                                     .inv(fcall("factorial", list(.sub(nu, S_One)))), .pow(unpolarify(z), .sub(nu, S_One))))))
  }
  .add(list(.mul(list(.sub(fcall("exp", list(.mul(list(S_Two, S_I, S_Pi, nu, n)))), S_One), .pow(z, .sub(nu, S_One)), fcall("gamma", list(.sub(S_One, nu))))),
            fcall("expint", list(nu, z))))
}

# polygamma.eval
.polygamma_eval <- function(n, z) {
  if (is_nan(n) || is_nan(z)) return(S_NaN)
  if (is_const(z, "oo")) return(if (isTRUE(ask_fact(n, "zero"))) S_Infinity else S_Zero)
  if (z$type == "Integer" && z$value <= 0) return(S_ComplexInfinity)
  if (identical(n$key, S_NegativeOne$key))
    return(.sub(fcall("loggamma", list(z)), .div(fcall("log", list(.mul(list(S_Two, S_Pi)))), S_Two)))
  if (isTRUE(ask_fact(n, "zero"))) {
    zi <- extract_multiplicatively(z, S_I)
    if (is_const(z, "-oo") || (!is.null(zi) && (is_const(zi, "oo") || is_const(zi, "-oo")))) return(S_Infinity)
    if (z$type == "Integer") return(.sub(fcall("harmonic", list(.sub(z, S_One))), S_EulerGamma))
    if (z$type == "Rational") {
      if (gmp::denominator(gmp::as.bigq(z$value)) <= 6) return(.polygamma_expand(S_Zero, z))
    }
    if (z$type == "Float") return(.float_eval("polygamma", list(n, z)))
    return(NULL)
  }
  if (isTRUE(ask_fact(n, "integer")) && isTRUE(ask_fact(n, "nonnegative"))) {
    nz <- unpolarify(z)
    if (!identical(nz$key, z$key)) return(fcall("polygamma", list(n, nz)))
    if (z$type == "Integer")
      return(.mul(list(.pow(S_NegativeOne, .add(list(n, S_One))), fcall("factorial", list(n)), fcall("zeta", list(.add(list(n, S_One)), z)))))
    if (identical(z$key, S_Half$key))
      return(.mul(list(.pow(S_NegativeOne, .add(list(n, S_One))), fcall("factorial", list(n)),
                       .sub(.pow(S_Two, .add(list(n, S_One))), S_One), fcall("zeta", list(.add(list(n, S_One)))))))
  }
  if (z$type == "Float" && n$type == "Integer") return(.float_eval("polygamma", list(n, z)))
  NULL
}

# polygamma._eval_expand_func
.polygamma_expand <- function(n, z) {
  if (n$type == "Integer" && n$value >= 0) {
    if (z$type == "Add") {
      coeff <- .sympy_args(z)[[1L]]
      if (coeff$type == "Integer") {
        e <- .neg(.add(list(n, S_One)))
        cv <- .int_val(coeff)
        tail <- if (cv > 0) .add(lapply(seq_len(cv), function(i) .pow(.sub(z, .int(i)), e)))
                else .neg(.add(lapply(seq_len(-cv) - 1L, function(i) .pow(.add(list(z, .int(i))), e))))
        return(.add(list(fcall("polygamma", list(n, .sub(z, coeff))),
                         .mul(list(.pow(S_NegativeOne, n), fcall("factorial", list(n)), tail)))))
      }
    } else if (z$type == "Mul") {
      ag <- .sympy_args(z)
      coeff <- ag[[1L]]; zz <- .mul(ag[-1L])
      if (coeff$type == "Integer" && coeff$value > 0) {
        cv <- .int_val(coeff)
        tl <- lapply(seq_len(cv) - 1L, function(i) fcall("polygamma", list(n, .add(list(zz, Rational(i, cv))))))
        if (is_zero_num(n)) return(.add(list(.div(.add(tl), coeff), fcall("log", list(coeff)))))
        return(.div(.add(tl), .pow(coeff, .add(list(n, S_One)))))
      }
    }
  }
  if (is_zero_num(n) && z$type %in% c("Integer", "Rational")) {
    q <- gmp::as.bigq(z$value); P <- .zs(gmp::numerator(q)); Q <- as.integer(as.character(gmp::denominator(q)))
    part1 <- .add(c(list(.neg(S_EulerGamma), .neg(.div(.mul(list(S_Pi, fcall("cot", list(.div(.mul(list(P, S_Pi)), .int(Q)))))), S_Two)),
                         .neg(fcall("log", list(.int(Q))))),
                    lapply(seq_len(Q - 1L), function(k) .mul(list(fcall("cos", list(.div(.mul(list(.int(2L * k), S_Pi, P)), .int(Q)))),
                                                                  fcall("log", list(.mul(list(S_Two, fcall("sin", list(.div(.mul(list(.int(k), S_Pi)), .int(Q)))))))))))))
    if (num_sign(z) > 0) {
      nn <- floor(as.double(q)); z0 <- .sub(z, .int(nn))
      return(.add(c(list(part1), lapply(seq_len(nn) - 1L, function(k) .inv(.add(list(z0, .int(k))))))))
    }
    if (num_sign(z) < 0) {
      nn <- floor(1 - as.double(q)); z0 <- .add(list(z, .int(nn)))
      return(.sub(part1, .add(lapply(seq_len(nn) - 1L, function(k) .inv(.sub(.sub(z0, S_One), .int(k)))))))
    }
  }
  if (identical(n$key, S_NegativeOne$key))
    return(.sub(fcall("loggamma", list(z)), .div(fcall("log", list(.mul(list(S_Two, S_Pi)))), S_Two)))
  NULL
}
.expand_hooks$polygamma <- function(e) .polygamma_expand(e$args[[1L]], e$args[[2L]])

# zeta.eval
.zeta_eval <- function(args) {
  s <- args[[1L]]; a <- if (length(args) == 2L) args[[2L]] else NULL
  if (!is.null(a) && is_one_num(a)) return(fcall("zeta", list(s)))
  if (is_nan(s) || (!is.null(a) && is_nan(a))) return(S_NaN)
  if (is_one_num(s)) return(S_ComplexInfinity)
  if (is_const(s, "oo")) return(S_One)
  if (!is.null(a) && is_const(a, "oo")) return(S_Zero)
  sint <- s$type == "Integer"
  aa <- if (is.null(a)) S_One else a
  if (sint && s$value <= 0) return(.div(fcall("bernoulli", list(.sub(S_One, s), aa)), .sub(s, S_One)))
  if (is_one_num(aa)) {
    if (sint && .int_val(s) %% 2 == 0)
      return(.neg(.div(.mul(list(.pow(.mul(list(S_Two, S_Pi, S_I)), s), fcall("bernoulli", list(s)))),
                       .mul(list(S_Two, fcall("factorial", list(s)))))))
  } else if (sint && aa$type == "Integer" && aa$value > 0) {
    return(.sub(fcall("zeta", list(s)), fcall("harmonic", list(.sub(aa, S_One), s))))
  } else if (aa$type == "Integer" && aa$value <= 0 &&
             (identical(ask_fact(s, "integer"), FALSE) || identical(ask_fact(s, "nonpositive"), FALSE))) {
    return(S_NaN)
  }
  if (s$type == "Float" && (is.null(a) || is_number_type(a))) return(.float_eval("zeta", args))
  NULL
}
