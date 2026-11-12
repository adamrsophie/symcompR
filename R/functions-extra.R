# Further special functions of sympy/functions: Rem, andre, multigamma,
# stieltjes, riemann_xi, owens_t, erf2inv, marcumq, hyper, meijerg,
# appellf1, spherical Hankel functions, SingularityFunction, incomplete beta,
# Jacobi theta, normalized Jacobi polynomials, conjugate and real spherical
# harmonics, Mathieu functions and the polar-number helpers.

.fx <- function(name) {
  force(name)
  function(...) fcall(name, lapply(list(...), as_sym))
}

#' More special functions
#'
#' `Rem(p, q)` is the remainder with the sign of `p`. `andre(n)` gives the
#' Euler zigzag numbers. `multigamma(x, p)` is the multivariate gamma
#' function, `stieltjes(n, a)` the Stieltjes constants, `riemann_xi(s)`
#' Riemann's xi function, `owens_t(h, a)` Owen's T function, `erf2inv(x, y)`
#' the inverse of `erf2()` in its second argument, `marcumq(m, a, b)` the
#' Marcum Q function. `hyper(ap, bq, z)` is the generalized hypergeometric
#' function and `meijerg(a, b, z)` the Meijer G function (with `a` and `b`
#' lists of two lists), `appellf1()` Appell's F1. `hn1()` and `hn2()` are the
#' spherical Hankel functions and `jn_zeros(n, k)` the first `k` zeros of
#' [jn()]. `SingularityFunction(x, a, n)` is the Macaulay bracket. `betainc()`
#' and `betainc_regularized()` are the generalized incomplete beta function
#' between two limits. `jtheta(n, z, q)` is the Jacobi theta function,
#' `jacobi_normalized()` the orthonormal Jacobi polynomial, `Ynm_c()` and
#' `Znm()` the conjugate and real spherical harmonics. `mathieus()`,
#' `mathieuc()`, `mathieusprime()` and `mathieucprime()` are symbolic only.
#' `exp_polar()`, `polar_lift()`, `unpolarify()`, `periodic_argument()`,
#' `unbranched_argument()` and `principal_branch()` work with numbers on the
#' Riemann surface of the logarithm. `piecewise_exclusive()` makes the
#' conditions of a [Piecewise()] mutually exclusive.
#'
#' @param ... Arguments of the symbolic function, in SymPy's order.
#' @param n,x,a,b,m,z,k,theta,phi,period,e Arguments.
#' @param ap,bq Lists of upper and lower parameters.
#' @param subs For `unpolarify()`: a named list of values substituted for
#'   symbols.
#' @param exponents_only For `unpolarify()`: only unpolarify exponents.
#' @return An expression, or a list of Floats for `jn_zeros()`.
#' @examples
#' x <- Symbol("x")
#' Rem(-7, 3)
#' andre(8)
#' multigamma(x, 2)
#' SingularityFunction(5, 2, 3)
#' N(hyper(list(1, 2), list(3), Rational(1, 2)))
#' jn_zeros(1, 2)
#' @name special-functions-extra
NULL

def_function("Rem", nargs = 2L, eval = function(args) {
  p <- args[[1L]]; q <- args[[2L]]
  if (is_zero_num(q)) return(S_NaN)
  if (is_exact_number(p) && is_exact_number(q)) {
    r <- p$value / q$value
    tq <- gmp::as.bigq(gmp::numerator(r) %/% gmp::denominator(r))
    if (r < 0 && tq != r) tq <- tq + 1
    return(make_rational(p$value - q$value * tq))
  }
  NULL
}, numeric = function(p, q) p - q * trunc(p / q))
#' @rdname special-functions-extra
#' @export
Rem <- .fx("Rem")

# Euler zigzag numbers by the Seidel-Entringer triangle.
.andre <- function(n) {
  row <- list(gmp::as.bigz(1))
  if (n == 0) return(row[[1L]])
  for (k in seq_len(n)) {
    nr <- list(gmp::as.bigz(0))
    for (j in seq_len(k)) nr[[j + 1L]] <- nr[[j]] + row[[k - j + 1L]]
    row <- nr
  }
  row[[length(row)]]
}
def_function("andre", eval = function(args) {
  n <- args[[1L]]
  if (n$type == "Integer" && n$value >= 0) return(make_rational(gmp::as.bigq(.andre(as.numeric(n$value)))))
  NULL
})
#' @rdname special-functions-extra
#' @export
andre <- .fx("andre")

def_function("multigamma", nargs = 2L, eval = function(args) {
  x <- args[[1L]]; p <- args[[2L]]
  if (p$type == "Integer" && p$value >= 1) {
    pv <- as.numeric(p$value)
    g <- lapply(seq_len(pv), function(j) fcall("gamma", list(.add(list(x, .rat(1 - j, 2))))))
    return(.mul(c(list(.pow(S_Pi, .rat(pv * (pv - 1), 4))), g)))
  }
  NULL
})
#' @rdname special-functions-extra
#' @export
multigamma <- .fx("multigamma")

def_function("stieltjes", nargs = 1:2, eval = function(args) {
  if (is_zero_num(args[[1L]]) && (length(args) == 1L || is_one_num(args[[2L]]))) return(S_EulerGamma)
  NULL
})
#' @rdname special-functions-extra
#' @export
stieltjes <- .fx("stieltjes")

.riemann_xi_formula <- function(s) .mul(list(s, .sub(s, S_One), S_Half, .pow(S_Pi, .neg(.div(s, S_Two))),
                                             fcall("gamma", list(.div(s, S_Two))), fcall("zeta", list(s))))
def_function("riemann_xi", eval = function(args) {
  s <- args[[1L]]
  if (is_zero_num(s) || is_one_num(s)) return(S_Half)
  if (s$type == "Integer") {
    z <- fcall("zeta", list(s))
    if (z$type != "Function") return(simplify(.riemann_xi_formula(s)))
  }
  NULL
}, numeric = function(s) .neval(.riemann_xi_formula(as_sym(s)), list()))
#' @rdname special-functions-extra
#' @export
riemann_xi <- .fx("riemann_xi")

def_function("owens_t", nargs = 2L, eval = function(args) {
  h <- args[[1L]]; a <- args[[2L]]
  if (is_zero_num(h)) return(.div(fcall("atan", list(a)), .mul(list(S_Two, S_Pi))))
  if (is_zero_num(a)) return(S_Zero)
  NULL
}, numeric = function(h, a) stats::integrate(function(t) exp(-Re(h)^2 * (1 + t^2) / 2) / (1 + t^2),
                                             0, Re(a), rel.tol = 1e-12)$value / (2 * pi))
#' @rdname special-functions-extra
#' @export
owens_t <- .fx("owens_t")

def_function("erf2inv", nargs = 2L, eval = function(args) {
  x <- args[[1L]]; y <- args[[2L]]
  if (is_zero_num(x)) return(fcall("erfinv", list(y)))
  if (is_zero_num(y)) return(x)
  NULL
}, numeric = function(x, y) stats::qnorm((Re(y) + 2 * stats::pnorm(Re(x) * sqrt(2)) - 1 + 1) / 2) / sqrt(2))
#' @rdname special-functions-extra
#' @export
erf2inv <- .fx("erf2inv")

def_function("marcumq", nargs = 3L, eval = function(args) {
  m <- args[[1L]]; a <- args[[2L]]; b <- args[[3L]]
  if (is_one_num(m) && is_zero_num(a)) return(fcall("exp", list(.neg(.div(.sq(b), S_Two)))))
  if (is_zero_num(a) && isTRUE(ask_fact(m, "integer")) && isTRUE(ask_fact(m, "positive")))
    return(.div(fcall("uppergamma", list(m, .div(.sq(b), S_Two))), fcall("gamma", list(m))))
  NULL
}, numeric = function(m, a, b) {
  m <- Re(m); a <- Re(a); b <- Re(b)
  if (a == 0) return(stats::pgamma(b^2 / 2, m, lower.tail = FALSE))
  stats::integrate(function(t) t * (t / a)^(m - 1) * exp(-(t - a)^2 / 2) * besselI(a * t, m - 1, expon.scaled = TRUE),
                   b, Inf, rel.tol = 1e-12)$value
})
#' @rdname special-functions-extra
#' @export
marcumq <- .fx("marcumq")

# Generalized hypergeometric series, |z| < 1 or terminating.
.pfq <- function(ap, bq, z) {
  s <- 1; term <- 1
  for (k in 0:20000) {
    term <- term * prod(ap + k) / prod(bq + k) * z / (k + 1)
    s <- s + term
    if (Mod(term) < 1e-17 * Mod(s)) return(s)
  }
  stop("hyper() series did not converge; |z| is too close to 1", call. = FALSE)
}
# Tuples evaluate to numeric vectors, for hyper() parameters.
register_neval("Tuple", function(e, env) unlist(lapply(e$args, function(a) .neval(a, env))))


# appellf1.eval: b1 <= b2 in default_sort_key order (x and y swapped with
# them), and appellf1(..., 0, 0) = 1
.appellf1_eval <- function(args) {
  a <- args[[1L]]; b1 <- args[[2L]]; b2 <- args[[3L]]; c <- args[[4L]]; x <- args[[5L]]; y <- args[[6L]]
  if (.key_cmp(.sort_key(b1), .sort_key(b2)) > 0) return(fcall("appellf1", list(a, b2, b1, c, y, x)))
  if (identical(b1$key, b2$key) && .key_cmp(.sort_key(x), .sort_key(y)) > 0) return(fcall("appellf1", list(a, b1, b2, c, y, x)))
  if (is_zero_num(x) && is_zero_num(y)) return(S_One)
  NULL
}
.appellf1_diff <- function(args, i) {
  a <- args[[1L]]; b1 <- args[[2L]]; b2 <- args[[3L]]; c <- args[[4L]]
  one <- function(v) .add(list(v, S_One))
  if (i == 5L) return(.mul(list(a, b1, .inv(c), fcall("appellf1", list(one(a), one(b1), b2, one(c), args[[5L]], args[[6L]])))))
  if (i == 6L) return(.mul(list(a, b2, .inv(c), fcall("appellf1", list(one(a), b1, one(b2), one(c), args[[5L]], args[[6L]])))))
  NULL
}
def_function("appellf1", nargs = 6L, eval = .appellf1_eval, diff = .appellf1_diff, numeric = function(a, b1, b2, c, x, y) {
  K <- 400
  # signed Pochhammer ratios by recurrence
  r <- cumprod(c(1, (a + 0:(K - 1)) / (c + 0:(K - 1))))
  p1 <- cumprod(c(1, (b1 + 0:(K - 1)) / (1:K))) * x^(0:K)
  p2 <- cumprod(c(1, (b2 + 0:(K - 1)) / (1:K))) * y^(0:K)
  s <- 0
  for (i in 0:K) {
    js <- 0:(K - i)
    s <- s + sum(r[i + js + 1] * p1[i + 1] * p2[js + 1])
  }
  s
})
#' @rdname special-functions-extra
#' @export
appellf1 <- .fx("appellf1")

.sph_hankel <- function(sgn) function(nu, z) {
  nu <- Re(nu); z <- Re(z)
  sqrt(pi / (2 * z)) * complex(real = besselJ(z, nu + 0.5), imaginary = sgn * besselY(z, nu + 0.5))
}
def_function("hn1", nargs = 2L, numeric = .sph_hankel(1))
def_function("hn2", nargs = 2L, numeric = .sph_hankel(-1))
#' @rdname special-functions-extra
#' @export
hn1 <- .fx("hn1")
#' @rdname special-functions-extra
#' @export
hn2 <- .fx("hn2")

#' @rdname special-functions-extra
#' @export
jn_zeros <- function(n, k) {
  n <- as.numeric(n); k <- as.integer(k)
  f <- function(x) sqrt(pi / (2 * x)) * besselJ(x, n + 0.5)
  out <- numeric()
  # zeros of j_n lie above n + 1/2; start there so besselJ does not underflow
  lo <- max(1e-3, n)
  while (length(out) < k) {
    hi <- lo + 0.25
    fl <- f(lo); fh <- f(hi)
    if (fl != 0 && fh != 0 && sign(fl) != sign(fh)) out <- c(out, stats::uniroot(f, c(lo, hi), tol = 1e-15)$root)
    lo <- hi
  }
  lapply(out, Float)
}

def_function("SingularityFunction", nargs = 3L,
  eval = function(args) {
    x <- args[[1L]]; a <- args[[2L]]; n <- args[[3L]]
    d <- .sub(x, a)
    if (!is_number_type(n)) return(NULL)
    neg <- num_sign(n) < 0
    if (is_zero_num(d)) return(if (neg) S_Infinity else if (is_zero_num(n)) S_One else S_Zero)
    pos <- ask_fact(d, "positive")
    if (isTRUE(pos)) return(if (neg) S_Zero else .pow(d, n))
    if (identical(pos, FALSE) && isTRUE(ask_fact(d, "nonzero"))) return(S_Zero)
    NULL
  },
  diff = function(args, i) {
    if (i != 1L) return(NULL)
    n <- args[[3L]]
    if (isTRUE(ask_fact(n, "positive"))) return(.mul(list(n, fcall("SingularityFunction", list(args[[1L]], args[[2L]], .sub(n, S_One))))))
    fcall("SingularityFunction", list(args[[1L]], args[[2L]], .sub(n, S_One)))
  })
#' @rdname special-functions-extra
#' @export
SingularityFunction <- .fx("SingularityFunction")

.betainc_num <- function(a, b, x1, x2) {
  if (Im(x1) == 0 && Im(x2) == 0)
    return(stats::integrate(function(t) t^(Re(a) - 1) * (1 - t)^(Re(b) - 1), Re(x1), Re(x2), rel.tol = 1e-12)$value)
  # complex limits: integrate along the straight segment from x1 to x2
  x1 <- as.complex(x1); d <- as.complex(x2) - x1
  g <- function(u) { t <- x1 + u * d; t^(a - 1) * (1 - t)^(b - 1) * d }
  complex(real = stats::integrate(function(u) Re(g(u)), 0, 1, rel.tol = 1e-12)$value,
          imaginary = stats::integrate(function(u) Im(g(u)), 0, 1, rel.tol = 1e-12)$value)
}
def_function("betainc", nargs = 4L, numeric = .betainc_num)
def_function("betainc_regularized", nargs = 4L, numeric = function(a, b, x1, x2)
  .betainc_num(a, b, x1, x2) / base::beta(Re(a), Re(b)))
#' @rdname special-functions-extra
#' @export
betainc <- .fx("betainc")
#' @rdname special-functions-extra
#' @export
betainc_regularized <- .fx("betainc_regularized")

def_function("jtheta", nargs = 3L, numeric = function(n, z, q) {
  n <- Re(n); k <- 0:200
  switch(as.character(n),
         "1" = 2 * sum((-1)^k * q^((k + 0.5)^2) * sin((2 * k + 1) * z)),
         "2" = 2 * sum(q^((k + 0.5)^2) * cos((2 * k + 1) * z)),
         "3" = 1 + 2 * sum(q^((k[-1])^2) * cos(2 * k[-1] * z)),
         "4" = 1 + 2 * sum((-1)^k[-1] * q^((k[-1])^2) * cos(2 * k[-1] * z)),
         stop("jtheta() order must be 1, 2, 3 or 4"))
})
#' @rdname special-functions-extra
#' @export
jtheta <- .fx("jtheta")

#' @rdname special-functions-extra
#' @export
jacobi_normalized <- function(n, a, b, x) {
  n <- as_sym(n); a <- as_sym(a); b <- as_sym(b)
  nf <- .div(.mul(list(.pow(S_Two, .add(list(a, b, S_One))), fcall("gamma", list(.add(list(n, a, S_One)))),
                       fcall("gamma", list(.add(list(n, b, S_One)))))),
             .mul(list(.add(list(.mul(list(S_Two, n)), a, b, S_One)), fcall("gamma", list(.add(list(n, a, b, S_One)))),
                       fcall("factorial", list(n)))))
  .div(jacobi(n, a, b, x), .sqrt(nf))
}

#' @rdname special-functions-extra
#' @export
Ynm_c <- function(n, m, theta, phi) {
  m <- as_sym(m)
  .mul(list(.pow(S_NegativeOne, m), fcall("Ynm", list(as_sym(n), .neg(m), as_sym(theta), as_sym(phi)))))
}

#' @rdname special-functions-extra
#' @export
Znm <- function(n, m, theta, phi) {
  m <- as_sym(m)
  y <- fcall("Ynm", list(as_sym(n), m, as_sym(theta), as_sym(phi)))
  if (is_zero_num(m)) return(y)
  yc <- Ynm_c(n, m, theta, phi)
  if (isTRUE(ask_fact(m, "positive"))) return(.div(.add(list(y, yc)), .sqrt(S_Two)))
  if (isTRUE(ask_fact(m, "negative"))) return(.div(.sub(y, yc), .mul(list(.sqrt(S_Two), S_I))))
  fcall("Znm", list(as_sym(n), m, as_sym(theta), as_sym(phi)))
}

for (.mf in c("mathieus", "mathieuc", "mathieusprime", "mathieucprime")) def_function(.mf, nargs = 3L)
#' @rdname special-functions-extra
#' @export
mathieus <- .fx("mathieus")
#' @rdname special-functions-extra
#' @export
mathieuc <- .fx("mathieuc")
#' @rdname special-functions-extra
#' @export
mathieusprime <- .fx("mathieusprime")
#' @rdname special-functions-extra
#' @export
mathieucprime <- .fx("mathieucprime")

# ---------------------------------------------------------------------------
# Polar numbers

# exp_polar._eval_evalf: only arguments with imaginary part in (-pi, pi]
.exp_polar_num <- function(x) {
  i <- Im(x)
  if (i <= -pi || i > pi) stop("exp_polar() cannot be evaluated off the principal branch", call. = FALSE)
  r <- exp(x)
  if (i > 0 && Im(r) < 0) Re(r) else r
}
def_function("exp_polar", numeric = .exp_polar_num, facts = function(args, afs) {
  out <- c(commutative = TRUE)
  f <- afs[[1L]]
  # exp_polar._eval_is_extended_real and ExpBase._eval_is_finite
  if (isTRUE(f[["real"]])) { out["real"] <- TRUE; out["positive"] <- TRUE }
  if (isTRUE(f[["finite"]])) out["finite"] <- TRUE
  # Expr._eval_is_extended_positive_negative: the sign of the value
  if (!length(free_symbols(args[[1L]]))) {
    v <- tryCatch(.exp_polar_num(.neval(args[[1L]], list())), error = function(err) NULL)
    if (!is.null(v) && is.finite(Re(v))) {
      if (abs(Im(v)) > 1e-12 * max(1, abs(Re(v)))) { out["extended_positive"] <- FALSE; out["extended_negative"] <- FALSE }
      else if (Re(v) > 0) out["positive"] <- TRUE
      else if (Re(v) < 0) out["negative"] <- TRUE
    }
  }
  out
})
#' @rdname special-functions-extra
#' @export
exp_polar <- .fx("exp_polar")

# polar_lift.eval: numbers of argument 0, pi/2, -pi/2 or pi, and the polar
# and positive factors of a product
.is_polar <- function(z) {
  if (z$type == "Function" && z$name %in% c("exp_polar", "polar_lift")) return(TRUE)
  if (z$type == "Mul") return(all(vapply(z$args, function(f) .is_polar(f) || isTRUE(ask_fact(f, "positive")), TRUE)) &&
                              any(vapply(z$args, .is_polar, TRUE)))
  if (z$type == "Pow") return(.is_polar(z$args[[1L]]))
  is_symbol_like(z) && isTRUE(z$extra$facts["polar"])
}
def_function("polar_lift", eval = function(args) {
  z <- args[[1L]]
  if (!length(free_symbols(z)) && !.has_type(z, "AccumulationBounds")) {
    ar <- fcall("arg", list(z))
    for (v in list(S_Zero, .div(S_Pi, S_Two), .neg(.div(S_Pi, S_Two)), S_Pi))
      if (identical(ar$key, v$key)) return(.mul(list(fcall("exp_polar", list(.mul(list(S_I, ar)))), fcall("Abs", list(z)))))
  }
  fs <- if (z$type == "Mul") z$args else list(z)
  inc <- list(); exc <- list(); pos <- list()
  for (f in fs) {
    if (.is_polar(f)) inc[[length(inc) + 1L]] <- f
    else if (isTRUE(ask_fact(f, "positive"))) pos[[length(pos) + 1L]] <- f
    else exc[[length(exc) + 1L]] <- f
  }
  if (length(exc) < length(fs)) {
    if (length(exc)) return(.mul(c(inc, pos, list(fcall("polar_lift", list(.mul(exc)))))))
    if (length(inc)) return(.mul(c(inc, pos)))
    return(.mul(c(pos, list(fcall("exp_polar", list(S_Zero))))))
  }
  NULL
}, numeric = function(z) z)
#' @rdname special-functions-extra
#' @export
polar_lift <- .fx("polar_lift")

#' @rdname special-functions-extra
#' @export
unpolarify <- function(e, subs = NULL, exponents_only = FALSE) {
  if (is.logical(e)) return(e)
  e <- as_sym(e)
  if (!is.null(subs)) return(unpolarify(subs(e, subs), exponents_only = exponents_only))
  repeat {
    res <- .unpolarify(e, exponents_only, exponents_only)
    if (identical(res$key, e$key)) break
    e <- res
  }
  # exp_polar(0) is 1 and polar_lift(0) is 0
  .map_bottom_up(res, function(z) {
    if (z$type == "Function" && z$name %in% c("exp_polar", "polar_lift") && is_zero_num(z$args[[1L]]))
      return(if (z$name == "exp_polar") S_One else S_Zero)
    z
  })
}

# functions with unbranched = True in SymPy
.UNBRANCHED <- c("mathieus", "mathieuc", "mathieusprime", "mathieucprime", "beta", "betainc",
                 "betainc_regularized", "gamma", "multigamma", "erf", "erfc", "erfi", "fresnels",
                 "fresnelc", "airyai", "airybi", "airyaiprime", "airybiprime", "re", "im", "Abs",
                 "sin", "cos", "tan", "cot", "sec", "csc", "exp", "exp_polar",
                 "sinh", "cosh", "tanh", "coth", "sech", "csch")

# SymPy _unpolarify
.unpolarify <- function(eq, exponents_only, pause = FALSE) {
  if (!length(eq$args)) return(eq)
  rec <- function(x, p = FALSE) .unpolarify(x, exponents_only, p)
  t <- eq$type
  if (!pause) {
    if (t == "Function" && eq$name == "exp_polar") return(fcall("exp", list(rec(eq$args[[1L]]))))
    if (t == "Function" && eq$name == "principal_branch" && identical(eq$args[[2L]]$key, .mul(list(S_Two, S_Pi))$key))
      return(rec(eq$args[[1L]]))
    if (t %in% c("Add", "Mul") || t %in% .BOOL_TYPES ||
        (is_relational(eq) && (!(t %in% c("Equality", "Unequality")) || any(vapply(eq$args, is_zero_num, TRUE)))))
      return(sym_rebuild(eq, lapply(eq$args, rec)))
    if (t == "Function" && eq$name == "polar_lift") return(rec(eq$args[[1L]]))
  }
  if (t == "Pow") {
    expo <- rec(eq$args[[2L]])
    base <- rec(eq$args[[1L]], !(isTRUE(ask_fact(expo, "integer")) && !pause))
    return(.pow(base, expo))
  }
  if (t == "Function" && eq$name %in% .UNBRANCHED)
    return(sym_rebuild(eq, lapply(eq$args, function(x) rec(x, exponents_only))))
  sym_rebuild(eq, lapply(eq$args, function(x) rec(x, TRUE)))
}

# Argument of z on the Riemann surface: exp_polar(I*t) contributes t.
# periodic_argument._getunbranched
.getunbranched <- function(ar) {
  fs <- if (ar$type == "Mul") ar$args else list(ar)
  u <- S_Zero
  for (a in fs) {
    if (!.is_polar(a)) u <- .add(list(u, fcall("arg", list(a))))
    else if (a$type == "Function" && a$name == "exp_polar") u <- .add(list(u, .re_im(a$args[[1L]])[[2L]]))
    else if (a$type == "Pow") {
      ri <- .re_im(a$args[[2L]])
      u <- .add(list(u, .mul(list(ri[[1L]], unbranched_argument(a$args[[1L]]))), .mul(list(ri[[2L]], fcall("log", list(Abs(a$args[[1L]])))))))
    } else if (a$type == "Function" && a$name == "polar_lift") u <- .add(list(u, fcall("arg", list(a$args[[1L]]))))
    else return(NULL)
  }
  u
}
# periodic_argument.eval
.periodic_argument_eval <- function(args) {
  ar <- args[[1L]]; period <- args[[2L]]
  if (!isTRUE(ask_fact(period, "extended_positive"))) return(NULL)
  if (is_const(period, "oo") && ar$type == "Function" && ar$name == "principal_branch") return(fcall("periodic_argument", ar$args))
  if (ar$type == "Function" && ar$name == "polar_lift" && .is_true(Ge(period, .mul(list(S_Two, S_Pi)))))
    return(fcall("periodic_argument", list(ar$args[[1L]], period)))
  if (ar$type == "Mul") {
    newargs <- Filter(function(z) !isTRUE(ask_fact(z, "positive")), ar$args)
    if (length(newargs) != length(ar$args)) return(fcall("periodic_argument", list(.mul(newargs), period)))
  }
  u <- .getunbranched(ar)
  if (is.null(u)) return(NULL)
  if (.has_function(u, "periodic_argument") || .has_function(u, "atan2") || .has_function(u, "atan")) return(NULL)
  if (is_const(period, "oo")) return(u)
  nn <- .mul(list(fcall("ceiling", list(.sub(.div(u, period), S_Half))), period))
  if (!.has_function(nn, "ceiling")) return(.sub(u, nn))
  NULL
}
def_function("periodic_argument", nargs = 2L, eval = .periodic_argument_eval)
.unbranched_arg <- function(z) fcall("periodic_argument", list(z, S_Infinity))

#' @rdname special-functions-extra
#' @export
unbranched_argument <- function(z) .unbranched_arg(as_sym(z))

#' @rdname special-functions-extra
#' @export
periodic_argument <- function(z, period) fcall("periodic_argument", list(as_sym(z), as_sym(period)))

def_function("principal_branch", nargs = 2L)

#' @rdname special-functions-extra
#' @export
principal_branch <- function(z, period) {
  z <- as_sym(z); period <- as_sym(period)
  # principal_branch.eval: polar_lift is dropped, an infinite period is no-op
  if (z$type == "Function" && z$name == "polar_lift") return(principal_branch(z$args[[1L]], period))
  if (is_const(period, "oo")) return(z)
  fs <- if (z$type == "Mul") z$args else list(z)
  mag <- list()
  for (f in fs) {
    if (f$type == "Function" && f$name == "exp_polar") next
    if (isTRUE(ask_fact(f, "positive"))) { mag[[length(mag) + 1L]] <- f; next }
    if (is_number_type(f) && num_sign(f) < 0) { mag[[length(mag) + 1L]] <- .neg(f); next }
    # unknown phase: leave unevaluated, as SymPy does
    return(fcall("principal_branch", list(z, period)))
  }
  .mul(c(mag, list(fcall("exp_polar", list(.mul(list(S_I, periodic_argument(z, period))))))))
}

# Drop bounds made redundant by tighter ones on the same symbol:
# x >= 1 & x >= 2 is x >= 2.
.prune_bounds <- function(a) {
  args <- a$args
  info <- lapply(args, function(r) {
    if (!(r$type %in% c("GreaterThan", "StrictGreaterThan", "LessThan", "StrictLessThan"))) return(NULL)
    if (!is_symbol_like(r$args[[1L]]) || !is_number_type(r$args[[2L]])) return(NULL)
    list(x = r$args[[1L]]$key, lower = r$type %in% c("GreaterThan", "StrictGreaterThan"),
         v = r$args[[2L]], strict = startsWith(r$type, "Strict"))
  })
  keep <- rep(TRUE, length(args))
  for (i in seq_along(args)) for (j in seq_along(args)) {
    a1 <- info[[i]]; a2 <- info[[j]]
    if (i == j || is.null(a1) || is.null(a2) || !keep[j] || a1$x != a2$x || a1$lower != a2$lower) next
    c21 <- .sgn_cmp(a2$v, a1$v)
    tighter <- if (a1$lower) c21 > 0 || (c21 == 0 && a2$strict && !a1$strict) else
      c21 < 0 || (c21 == 0 && a2$strict && !a1$strict)
    if (tighter) keep[i] <- FALSE
  }
  if (sum(keep) == 1L) args[keep][[1L]] else do.call(And, args[keep])
}

#' @rdname special-functions-extra
#' @export
piecewise_exclusive <- function(e) {
  e <- as_sym(e)
  if (e$type != "Piecewise") return(e)
  prev <- list()
  pieces <- lapply(e$args, function(pc) {
    ex <- pc$args[[1L]]; cond <- pc$args[[2L]]
    nc <- if (length(prev)) And(cond, Not(do.call(Or, prev))) else cond
    nc <- tryCatch(simplify_logic(nc), error = function(err) nc)
    if (nc$type == "And") nc <- .prune_bounds(nc)
    prev[[length(prev) + 1L]] <<- cond
    list(ex, nc)
  })
  # Piecewise(*newargs, evaluate=False)
  new_sym("Piecewise", lapply(pieces, function(p) new_sym("ExprCondPair", list(as_sym(p[[1L]]), .as_bool(p[[2L]])))))
}
