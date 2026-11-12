# Symbolic statistics and plotting, after sympy/stats and sympy/plotting.

.rv_counter <- new.env()
.rv_counter$n <- 0L

.new_rv <- function(name, dist) {
  .rv_counter$n <- .rv_counter$n + 1L
  facts <- .symbol_facts(c(real = TRUE))
  new_sym("Symbol", name = name, extra = list(facts = facts, dist = dist),
          key = paste0("S", nchar(name), ":", name, "#rv", .rv_counter$n))
}

is_rv <- function(e) e$type == "Symbol" && !is.null(e$extra$dist)

.rvs_in <- function(e) Filter(is_rv, free_symbols(e))

.cdist <- function(pdf, lo, hi, moment = NULL, cdf = NULL) {
  list(kind = "continuous", pdf = pdf, lo = lo, hi = hi, moment = moment, cdf = cdf)
}

.ddist <- function(pmf, values = NULL, lo = NULL, hi = NULL, moment = NULL) {
  list(kind = "discrete", pmf = pmf, values = values, lo = lo, hi = hi, moment = moment)
}

#' Random variables
#'
#' Constructors for random symbols with exact symbolic distributions, after
#' `sympy.stats`: continuous `Normal`, `Uniform`, `Exponential`, `GammaRV`,
#' `BetaRV`, `ChiSquared`, `StudentT`, `Cauchy`, `Laplace`, `LogNormal`,
#' `Weibull`, `Pareto`; discrete `Bernoulli`, `Binomial`, `Poisson`,
#' `Geometric`, `Die` and `Coin`. Parameters may be symbolic; give them
#' assumptions such as `positive = TRUE` so that integrals converge.
#'
#' `E(expr)` (alias `expectation()`) is the expected value, using closed
#' form moments when known and integration or summation otherwise.
#' `variance()`, `std()`, `moment()`, `cmoment()`, `skewness()`,
#' `kurtosis()`, `covariance()` and `correlation()` build on it.
#' `P(condition)` is the probability of an event in one random variable.
#' `rv_density(X)` gives the density or mass function and `rv_cdf(X)` the
#' distribution function, both as expressions in the variable `z`.
#' Random variables in a product or sum are treated as independent.
#' `rv_sample(X, n)` draws numeric samples.
#'
#' @param name Display name.
#' @param mu,sigma,a,b,rate,k,theta,alpha,beta,nu,x0,gamma,p,n,sides,lamda
#'   Distribution parameters.
#' @param expr An expression in random variables.
#' @param X,Y Random variables or expressions in them.
#' @param condition A relational in one random variable.
#' @param z Symbol for density and distribution functions.
#' @param size Number of draws.
#' @return A random symbol, expression or numeric vector.
#' @examples
#' mu <- Symbol("mu", real = TRUE); s <- Symbol("s", positive = TRUE)
#' X <- Normal("X", mu, s)
#' E(X^2)
#' variance(2 * X + 1)
#' Z <- Normal("Z", 0, 1)
#' P(Z > 0)
#' D <- Die("D", 6)
#' E(D)
#' @export
Normal <- function(name, mu, sigma) {
  mu <- as_sym(mu); sigma <- as_sym(sigma)
  pdf <- function(z) .mul(list(.inv(.mul(list(sigma, .sqrt(.mul(list(S_Two, S_Pi)))))),
                               fcall("exp", list(.neg(.div(.sq(.sub(z, mu)), .mul(list(S_Two, .sq(sigma)))))))))
  moment <- function(k) {
    # E[X^k] = sum C(k, 2j) mu^(k-2j) sigma^(2j) (2j-1)!!
    .add(lapply(0:(k %/% 2), function(j) {
      df <- if (j == 0) S_One else fcall("factorial2", list(.int(2 * j - 1)))
      .mul(list(fcall("binomial", list(.int(k), .int(2 * j))), .pow(mu, .int(k - 2 * j)), .pow(sigma, .int(2 * j)), df))
    }))
  }
  cdf <- function(z) .mul(list(S_Half, .add(list(S_One, fcall("erf", list(.div(.sub(z, mu), .mul(list(sigma, .sqrt(S_Two))))))))))
  .new_rv(name, .cdist(pdf, S_NegativeInfinity, S_Infinity, moment, cdf))
}

#' @rdname Normal
#' @export
Uniform <- function(name, a, b) {
  a <- as_sym(a); b <- as_sym(b)
  pdf <- function(z) .inv(.sub(b, a))
  moment <- function(k) .div(.sub(.pow(b, .int(k + 1)), .pow(a, .int(k + 1))), .mul(list(.int(k + 1), .sub(b, a))))
  cdf <- function(z) Piecewise(list(S_Zero, Lt(z, a)), list(.div(.sub(z, a), .sub(b, a)), Le(z, b)), list(S_One, S_true))
  d <- .cdist(pdf, a, b, moment, cdf)
  # SymPy's density carries the support; integration keeps the plain pdf
  d$density <- function(z) Piecewise(list(.inv(.sub(b, a)), And(Ge(z, a), Le(z, b))), list(S_Zero, S_true))
  .new_rv(name, d)
}

#' @rdname Normal
#' @export
Exponential <- function(name, rate) {
  r <- as_sym(rate)
  pdf <- function(z) .mul(list(r, fcall("exp", list(.neg(.mul(list(r, z)))))))
  moment <- function(k) .div(fcall("factorial", list(.int(k))), .pow(r, .int(k)))
  cdf <- function(z) Piecewise(list(.sub(S_One, fcall("exp", list(.neg(.mul(list(r, z)))))), Ge(z, S_Zero)), list(S_Zero, S_true))
  .new_rv(name, .cdist(pdf, S_Zero, S_Infinity, moment, cdf))
}

#' @rdname Normal
#' @export
GammaRV <- function(name, k, theta) {
  k <- as_sym(k); theta <- as_sym(theta)
  pdf <- function(z) .div(.mul(list(.pow(z, .sub(k, S_One)), fcall("exp", list(.neg(.div(z, theta)))))),
                          .mul(list(fcall("gamma", list(k)), .pow(theta, k))))
  moment <- function(m) .mul(list(.pow(theta, .int(m)), fcall("RisingFactorial", list(k, .int(m)))))
  .new_rv(name, .cdist(pdf, S_Zero, S_Infinity, moment))
}

#' @rdname Normal
#' @export
BetaRV <- function(name, alpha, beta) {
  a <- as_sym(alpha); b <- as_sym(beta)
  pdf <- function(z) .div(.mul(list(.pow(z, .sub(a, S_One)), .pow(.sub(S_One, z), .sub(b, S_One)))), fcall("beta", list(a, b)))
  moment <- function(m) .div(fcall("RisingFactorial", list(a, .int(m))), fcall("RisingFactorial", list(.add(list(a, b)), .int(m))))
  .new_rv(name, .cdist(pdf, S_Zero, S_One, moment))
}

#' @rdname Normal
#' @export
ChiSquared <- function(name, k) {
  k <- as_sym(k)
  X <- GammaRV(name, .div(k, S_Two), S_Two)
  X
}

#' @rdname Normal
#' @export
StudentT <- function(name, nu) {
  nu <- as_sym(nu)
  ex <- .div(.add(list(nu, S_One)), S_Two)
  pdf <- function(z) .mul(list(.inv(.mul(list(.sqrt(nu), fcall("beta", list(S_Half, .div(nu, S_Two)))))),
                               .pow(.add(list(S_One, .div(.sq(z), nu))), .neg(ex))))
  .new_rv(name, .cdist(pdf, S_NegativeInfinity, S_Infinity))
}

#' @rdname Normal
#' @export
Cauchy <- function(name, x0, gamma) {
  x0 <- as_sym(x0); g <- as_sym(gamma)
  pdf <- function(z) .inv(.mul(list(S_Pi, g, .add(list(S_One, .sq(.div(.sub(z, x0), g)))))))
  cdf <- function(z) .add(list(.mul(list(.inv(S_Pi), fcall("atan", list(.div(.sub(z, x0), g))))), S_Half))
  .new_rv(name, .cdist(pdf, S_NegativeInfinity, S_Infinity, NULL, cdf))
}

#' @rdname Normal
#' @export
Laplace <- function(name, mu, b) {
  mu <- as_sym(mu); b <- as_sym(b)
  pdf <- function(z) .mul(list(.inv(.mul(list(S_Two, b))), fcall("exp", list(.neg(.div(fcall("Abs", list(.sub(z, mu))), b))))))
  moment <- function(k) .add(lapply(seq(0, k, by = 2), function(j)
    .mul(list(fcall("binomial", list(.int(k), .int(j))), .pow(mu, .int(k - j)), fcall("factorial", list(.int(j))), .pow(b, .int(j))))))
  .new_rv(name, .cdist(pdf, S_NegativeInfinity, S_Infinity, moment))
}

#' @rdname Normal
#' @export
LogNormal <- function(name, mu, sigma) {
  mu <- as_sym(mu); s <- as_sym(sigma)
  pdf <- function(z) .mul(list(.inv(.mul(list(z, s, .sqrt(.mul(list(S_Two, S_Pi)))))),
                               fcall("exp", list(.neg(.div(.sq(.sub(fcall("log", list(z)), mu)), .mul(list(S_Two, .sq(s)))))))))
  moment <- function(k) fcall("exp", list(.add(list(.mul(list(.int(k), mu)), .div(.mul(list(.int(k * k), .sq(s))), S_Two)))))
  .new_rv(name, .cdist(pdf, S_Zero, S_Infinity, moment))
}

#' @rdname Normal
#' @export
Weibull <- function(name, alpha, beta) {
  a <- as_sym(alpha); b <- as_sym(beta)
  pdf <- function(z) .mul(list(.div(b, a), .pow(.div(z, a), .sub(b, S_One)), fcall("exp", list(.neg(.pow(.div(z, a), b))))))
  moment <- function(k) .mul(list(.pow(a, .int(k)), fcall("gamma", list(.add(list(S_One, .div(.int(k), b)))))))
  .new_rv(name, .cdist(pdf, S_Zero, S_Infinity, moment))
}

#' @rdname Normal
#' @export
Pareto <- function(name, x0, alpha) {
  xm <- as_sym(x0); a <- as_sym(alpha)
  pdf <- function(z) .div(.mul(list(a, .pow(xm, a))), .pow(z, .add(list(a, S_One))))
  moment <- function(k) {
    m <- .div(.mul(list(a, .pow(xm, .int(k)))), .sub(a, .int(k)))
    d <- .sub(a, .int(k))
    if (isTRUE(ask_fact(d, "positive"))) return(m)
    if (isTRUE(ask_fact(d, "nonpositive"))) return(S_Infinity)
    Piecewise(list(m, Gt(a, .int(k))), list(S_Infinity, TRUE))
  }
  .new_rv(name, .cdist(pdf, xm, S_Infinity, moment))
}

#' @rdname Normal
#' @export
Bernoulli <- function(name, p) {
  p <- as_sym(p)
  .new_rv(name, .ddist(NULL, values = list(list(S_One, p), list(S_Zero, .sub(S_One, p)))))
}

#' @rdname Normal
#' @export
Coin <- function(name, p = Rational(1, 2)) Bernoulli(name, p)

#' @rdname Normal
#' @export
Die <- function(name, sides = 6) {
  s <- as.integer(sides)
  .new_rv(name, .ddist(NULL, values = lapply(seq_len(s), function(k) list(.int(k), .rat(1L, s)))))
}

#' @rdname Normal
#' @export
Binomial <- function(name, n, p) {
  p <- as_sym(p)
  n <- as.integer(n)
  .new_rv(name, .ddist(NULL, values = lapply(0:n, function(k)
    list(.int(k), .mul(list(fcall("binomial", list(.int(n), .int(k))), .pow(p, .int(k)), .pow(.sub(S_One, p), .int(n - k))))))))
}

#' @rdname Normal
#' @export
Poisson <- function(name, lamda) {
  l <- as_sym(lamda)
  pmf <- function(k) .mul(list(.pow(l, k), fcall("exp", list(.neg(l))), .inv(fcall("factorial", list(k)))))
  moment <- function(m) .add(lapply(0:m, function(j) .mul(list(stirling(m, j), .pow(l, .int(j))))))
  .new_rv(name, .ddist(pmf, lo = S_Zero, hi = S_Infinity, moment = moment))
}

#' @rdname Normal
#' @export
Geometric <- function(name, p) {
  p <- as_sym(p)
  pmf <- function(k) .mul(list(p, .pow(.sub(S_One, p), .sub(k, S_One))))
  moment <- function(m) {
    if (m == 1) return(.inv(p))
    if (m == 2) return(.div(.sub(S_Two, p), .sq(p)))
    NULL
  }
  .new_rv(name, .ddist(pmf, lo = S_One, hi = S_Infinity, moment = moment))
}

# E[g(X)] for one random variable.
.E1 <- function(g, X) {
  d <- X$extra$dist
  g <- expand(g)
  if (is_polynomial(g, X)) {
    cl <- .coeff_list(g, X)
    n <- length(cl) - 1L
    if (!is.null(d$moment)) {
      used <- which(!vapply(cl, is_zero_num, TRUE))
      pows <- n - used + 1L
      ms <- lapply(pows, function(p) if (p == 0L) S_One else d$moment(p))
      if (!any(vapply(ms, is.null, TRUE))) {
        return(.add(Map(function(k, m) .mul(list(cl[[k]], m)), used, ms)))
      }
    }
  }
  if (d$kind == "discrete") {
    if (!is.null(d$values)) return(.add(lapply(d$values, function(v) .mul(list(subs(g, X, v[[1L]]), v[[2L]])))))
    k <- Dummy("k", integer = TRUE, nonnegative = TRUE)
    return(summation(.mul(list(subs(g, X, k), d$pmf(k))), list(k, d$lo, d$hi)))
  }
  z <- Dummy("z", real = TRUE)
  integrate(.mul(list(subs(g, X, z), d$pdf(z))), list(z, d$lo, d$hi))
}

#' @rdname Normal
#' @export
E <- function(expr) {
  expr <- expand(as_sym(expr))
  terms <- if (expr$type == "Add") expr$args else list(expr)
  out <- lapply(terms, function(t) {
    rvs <- .rvs_in(t)
    if (!length(rvs)) return(t)
    facs <- if (t$type == "Mul") t$args else list(t)
    res <- list()
    for (r in rvs) {
      dep <- Filter(function(f) .has1(f, r), facs)
      if (length(Filter(function(f) length(.rvs_in(f)) > 1L, dep))) stop("E() needs terms that factor over independent variables")
      res[[length(res) + 1L]] <- .E1(.mul(dep), r)
    }
    constant <- Filter(function(f) !length(.rvs_in(f)), facs)
    .mul(c(constant, res))
  })
  simplify_if_cheap(.add(out))
}

#' @rdname Normal
#' @export
expectation <- E

#' @rdname Normal
#' @export
moment <- function(X, k) E(.pow(as_sym(X), .int(k)))

#' @rdname Normal
#' @export
cmoment <- function(X, k) {
  X <- as_sym(X)
  E(.pow(.sub(X, E(X)), .int(k)))
}

#' @rdname Normal
#' @export
variance <- function(X) {
  m2 <- E(.sq(as_sym(X)))
  if (is_const(m2, "oo")) return(S_Infinity)
  simplify_if_cheap(.sub(m2, .sq(E(as_sym(X)))))
}

#' @rdname Normal
#' @export
std <- function(X) .sqrt(variance(X))

#' @rdname Normal
#' @export
covariance <- function(X, Y) simplify_if_cheap(.sub(E(.mul(list(as_sym(X), as_sym(Y)))), .mul(list(E(X), E(Y)))))

#' @rdname Normal
#' @export
correlation <- function(X, Y) simplify_if_cheap(.div(covariance(X, Y), .mul(list(std(X), std(Y)))))

#' @rdname Normal
#' @export
skewness <- function(X) simplify_if_cheap(.div(cmoment(X, 3), .pow(variance(X), .rat(3, 2))))

#' @rdname Normal
#' @export
kurtosis <- function(X) simplify_if_cheap(.div(cmoment(X, 4), .sq(variance(X))))

#' @rdname Normal
#' @export
rv_density <- function(X, z = Symbol("z")) {
  d <- X$extra$dist
  z <- as_sym(z)
  if (d$kind == "continuous") return(if (is.null(d$density)) d$pdf(z) else d$density(z))
  if (!is.null(d$pmf)) return(d$pmf(z))
  stop("rv_density(): finite distributions are given by their value list")
}

#' @rdname Normal
#' @export
rv_cdf <- function(X, z = Symbol("z")) {
  d <- X$extra$dist
  z <- as_sym(z)
  if (d$kind == "continuous" && !is.null(d$cdf)) return(d$cdf(z))
  stop("rv_cdf(): no closed form available for this distribution")
}

#' @rdname Normal
#' @export
P <- function(condition) {
  condition <- as_sym(condition)
  if (condition$type == "BooleanAtom") return(if (.is_true(condition)) S_One else S_Zero)
  rvs <- .rvs_in(condition)
  if (length(rvs) != 1L) stop("P() supports events in one random variable")
  X <- rvs[[1L]]
  d <- X$extra$dist
  if (d$kind == "discrete") {
    if (is.null(d$values)) stop("P() needs a finite discrete distribution")
    terms <- list()
    for (v in d$values) {
      c <- subs(condition, X, v[[1L]])
      if (.is_true(c)) terms[[length(terms) + 1L]] <- v[[2L]]
      else if (!.is_false(c)) terms[[length(terms) + 1L]] <- Piecewise(list(v[[2L]], c), list(S_Zero, TRUE))
    }
    return(.add(terms))
  }
  if (condition$type == "Equality") return(S_Zero)
  z <- Dummy("z", real = TRUE)
  cz <- subs(condition, X, z)
  support <- Interval(d$lo, d$hi)
  set <- solve_univariate_inequality(cz, z, relational = FALSE, domain = support)
  .set_prob(set, d, z)
}

.set_prob <- function(s, d, z) {
  pdf <- d$pdf(z)
  switch(s$type,
         EmptySet = , FiniteSet = S_Zero,
         SpecialSet = S_One,
         Interval = {
           a <- s$args[[1L]]; b <- s$args[[2L]]
           if (!is.null(d$cdf)) {
             Fb <- if (is_const(b, "oo")) S_One else if (identical(b$key, d$hi$key)) S_One else d$cdf(b)
             Fa <- if (is_const(a, "-oo")) S_Zero else if (identical(a$key, d$lo$key)) S_Zero else d$cdf(a)
             return(simplify_if_cheap(.sub(Fb, Fa)))
           }
           integrate(pdf, list(z, a, b))
         },
         Union = .add(lapply(s$args, .set_prob, d = d, z = z)),
         stop("cannot integrate over ", s$type))
}

#' @rdname Normal
#' @export
rv_sample <- function(X, size = 1) {
  d <- X$extra$dist
  if (d$kind == "discrete" && !is.null(d$values)) {
    vals <- vapply(d$values, function(v) as.numeric(v[[1L]]), 0)
    probs <- vapply(d$values, function(v) as.numeric(N(v[[2L]])), 0)
    return(sample(vals, size, replace = TRUE, prob = probs))
  }
  if (d$kind == "continuous") {
    lo <- as.numeric(N(d$lo)); hi <- as.numeric(N(d$hi))
    z <- Symbol("z")
    f <- lambdify(z, d$pdf(z))
    # grid around the mean, 40 standard deviations wide, inside the support
    m <- tryCatch(as.numeric(N(d$moment(1))), error = function(e) NA)
    sdv <- tryCatch(sqrt(as.numeric(N(d$moment(2))) - m^2), error = function(e) NA)
    if (!is.finite(m) || !is.finite(sdv) || sdv <= 0) { m <- 0; sdv <- 1.25 }
    a <- if (is.finite(lo)) lo else m - 40 * sdv
    b <- if (is.finite(hi)) hi else m + 40 * sdv
    grid <- seq(a, b, length.out = 20001L)
    w <- rep_len(f(grid), length(grid))
    w[!is.finite(w)] <- 0
    return(sample(grid, size, replace = TRUE, prob = w))
  }
  stop("rv_sample(): not available for this distribution")
}

# ---------------------------------------------------------------------------
# Plotting

#' Plot symbolic expressions with base graphics
#'
#' `plot(expr, range)` draws a function of one variable, with `range` given
#' as `list(x, a, b)` (default `-10` to `10`). `plot_parametric()`,
#' `plot3d()` (surface) and `plot_implicit()` (zero contour) are the
#' analogues of SymPy's plotting functions. Extra arguments go to the base
#' graphics functions.
#'
#' @param x A `sym` expression, for `plot()`.
#' @param range `list(var, a, b)`.
#' @param ... Graphical parameters.
#' @param fx,fy Coordinate expressions.
#' @param e An expression.
#' @param range1,range2 Ranges for two variables.
#' @param n Number of sample points.
#' @return The sampled data, invisibly.
#' @examples
#' x <- Symbol("x")
#' if (interactive()) plot(sin(x) / x, list(x, -10, 10))
#' @export
plot.sym <- function(x, range = NULL, ..., n = 400) {
  if (is.null(range)) {
    fs <- free_symbols(x)
    if (length(fs) != 1L) stop("plot() needs a range for multivariate expressions")
    range <- list(fs[[1L]], -10, 10)
  }
  v <- as_sym(range[[1L]])
  xs <- seq(as.numeric(N(as_sym(range[[2L]]))), as.numeric(N(as_sym(range[[3L]]))), length.out = n)
  ys <- .eval_grid(x, list(v), list(xs))
  graphics::plot(xs, ys, type = "l", xlab = .str(v), ylab = .str(x), ...)
  invisible(data.frame(x = xs, y = ys))
}

.eval_grid <- function(e, vars, vals) {
  out <- rep(NA_real_, length(vals[[1L]]))
  f <- tryCatch(do.call(lambdify, list(vars, e)), error = function(err) NULL)
  if (!is.null(f)) {
    r <- tryCatch(suppressWarnings(do.call(f, vals)), error = function(err) NULL)
    if (!is.null(r)) {
      r <- rep_len(r, length(out))
      return(if (is.complex(r)) ifelse(abs(Im(r)) < 1e-12, Re(r), NA_real_) else as.numeric(r))
    }
  }
  nms <- vapply(vars, function(s) s$name, "")
  for (i in seq_along(out)) {
    env <- stats::setNames(lapply(vals, `[`, i), nms)
    v <- tryCatch(.neval(e, env), error = function(err) NA_real_)
    out[i] <- if (is.complex(v)) { if (abs(Im(v)) < 1e-12) Re(v) else NA_real_ } else as.numeric(v)
  }
  out
}

#' @rdname plot.sym
#' @export
plot_parametric <- function(fx, fy, range, ..., n = 400) {
  v <- as_sym(range[[1L]])
  ts <- seq(as.numeric(N(as_sym(range[[2L]]))), as.numeric(N(as_sym(range[[3L]]))), length.out = n)
  xs <- .eval_grid(as_sym(fx), list(v), list(ts))
  ys <- .eval_grid(as_sym(fy), list(v), list(ts))
  graphics::plot(xs, ys, type = "l", xlab = .str(as_sym(fx)), ylab = .str(as_sym(fy)), ...)
  invisible(data.frame(t = ts, x = xs, y = ys))
}

.grid2 <- function(e, r1, r2, n) {
  v1 <- as_sym(r1[[1L]]); v2 <- as_sym(r2[[1L]])
  xs <- seq(as.numeric(N(as_sym(r1[[2L]]))), as.numeric(N(as_sym(r1[[3L]]))), length.out = n)
  ys <- seq(as.numeric(N(as_sym(r2[[2L]]))), as.numeric(N(as_sym(r2[[3L]]))), length.out = n)
  g <- expand.grid(xs, ys)
  z <- matrix(.eval_grid(as_sym(e), list(v1, v2), list(g[[1L]], g[[2L]])), n, n)
  list(x = xs, y = ys, z = z)
}

#' @rdname plot.sym
#' @export
plot3d <- function(e, range1, range2, ..., n = 50) {
  g <- .grid2(e, range1, range2, n)
  graphics::persp(g$x, g$y, g$z, xlab = .str(as_sym(range1[[1L]])), ylab = .str(as_sym(range2[[1L]])),
                  zlab = "z", theta = 30, phi = 25, ...)
  invisible(g)
}

#' @rdname plot.sym
#' @export
plot_implicit <- function(e, range1, range2, ..., n = 200) {
  e <- .eq_to_expr(e)
  g <- .grid2(e, range1, range2, n)
  graphics::contour(g$x, g$y, g$z, levels = 0, drawlabels = FALSE,
                    xlab = .str(as_sym(range1[[1L]])), ylab = .str(as_sym(range2[[1L]])), ...)
  invisible(g)
}
