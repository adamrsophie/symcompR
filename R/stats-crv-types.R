# Continuous distributions, after sympy/stats/crv_types.py.

.Rinf <- function() Interval(S_NegativeInfinity, S_Infinity)
.R0inf <- function() Interval(S_Zero, S_Infinity)

# A continuous distribution class with its named arguments and methods
.defcd <- function(.cls, .argnames, ...) .defdist(.cls, "SingleContinuousDistribution", .argnames, ...)

.rvc <- function(name, cls, args, check = TRUE) .rv_make(name, cls, args, check = check, single = "SingleContinuousPSpace")

.nie <- function(msg) .serr("NotImplementedError", msg)

# ---------------------------------------------------------------------------
# ContinuousRV and handmade distributions

.defcd("ContinuousDistributionHandmade", c("pdf"),
  set = function(d) d$args[[2L]],
  pdf = function(d, x) .rv_lcall(d$args[[1L]], x),
  check = function(pdf, set) {
    x <- Dummy("x")
    val <- do.call(integrate, list(.rv_lcall(pdf, x), .set_lims(x, set)))
    .value_check(!identical(.fuzzy_bool(Eq(val, 1)), FALSE), "The pdf on the given set is incorrect.")
  })

.ContinuousDistributionHandmade <- function(pdf, set = NULL) {
  if (is.null(set)) set <- .Rinf()
  new_sym("ContinuousDistributionHandmade", list(pdf, set))
}

#' Continuous random variables
#'
#' Constructors for the continuous distributions of `sympy.stats`. Each
#' returns a random symbol named `name` on a `SingleContinuousPSpace`.
#' `ContinuousRV(symbol, density, set)` builds one from a density
#' expression in `symbol`. `BetaRV()` and `GammaRV()` are SymPy's `Beta()`
#' and `Gamma()`, renamed so they do not mask base R; `Wald()` is
#' `GaussianInverse()` and `Erlang(k, l)` is `GammaRV(k, 1/l)`.
#'
#' @param name The name of the random variable (a string or symbol).
#' @param symbol A symbol for `ContinuousRV()`.
#' @param density A density expression in `symbol`.
#' @param set The support, an interval (the real line by default).
#' @param check Whether to check the arguments.
#' @param a,b,c,d,alpha,beta,sigma,lamda,left,right,x0,gamma,k,l,p,n,mu,mean,std,rate,d1,d2,s,m,theta,minimum,eta,omega,shape,xm,xi,nu,R Distribution parameters, named as in SymPy.
#' @return A random symbol.
#' @examples
#' X <- Normal("X", 0, 1)
#' rv_density(X)
#' E(Exponential("T", 2))
#' variance(Uniform("U", 0, 1))
#' @name continuous_rv
NULL

#' @rdname continuous_rv
#' @export
ContinuousRV <- function(symbol, density, set = NULL, check = FALSE) {
  if (is.null(set)) set <- .Rinf()
  pdf <- Piecewise(list(density, as_relational(set, symbol)), list(S_Zero, S_true))
  pdf <- Lambda(symbol, pdf)
  .rvc(symbol$name, "ContinuousDistributionHandmade", list(pdf, set), check = check)
}

# ---------------------------------------------------------------------------

.defcd("ArcsinDistribution", c("a", "b"),
  set = function(d) Interval(.da(d, "a"), .da(d, "b")),
  pdf = function(d, x) { a <- .da(d, "a"); b <- .da(d, "b"); 1 / (S_Pi * sqrt((x - a) * (b - x))) },
  `_cdf` = function(d, x) {
    a <- .da(d, "a"); b <- .da(d, "b")
    Piecewise(list(S_Zero, x < a), list(2 * asin(sqrt((x - a) / (b - a))) / S_Pi, x <= b), list(S_One, S_true))
  })

#' @rdname continuous_rv
#' @export
Arcsin <- function(name, a = 0, b = 1) .rvc(name, "ArcsinDistribution", list(a, b))

.defcd("BeniniDistribution", c("alpha", "beta", "sigma"),
  check = function(alpha, beta, sigma) {
    .value_check(alpha > 0, "Shape parameter Alpha must be positive.")
    .value_check(beta > 0, "Shape parameter Beta must be positive.")
    .value_check(sigma > 0, "Scale parameter Sigma must be positive.")
  },
  set = function(d) Interval(.da(d, "sigma"), S_Infinity),
  pdf = function(d, x) {
    p <- .dl(d)
    exp(-p$alpha * log(x / p$sigma) - p$beta * log(x / p$sigma)^2) * (p$alpha / x + 2 * p$beta * log(x / p$sigma) / x)
  },
  `_moment_generating_function` = function(d, t) .nie("The moment generating function of the Benini distribution does not exist."))

#' @rdname continuous_rv
#' @export
Benini <- function(name, alpha, beta, sigma) .rvc(name, "BeniniDistribution", list(alpha, beta, sigma))

.defcd("BetaDistribution", c("alpha", "beta"),
  set = function(d) Interval(0, 1),
  check = function(alpha, beta) {
    .value_check(alpha > 0, "Shape parameter Alpha must be positive.")
    .value_check(beta > 0, "Shape parameter Beta must be positive.")
  },
  pdf = function(d, x) { p <- .dl(d); x^(p$alpha - 1) * (1 - x)^(p$beta - 1) / sym_beta(p$alpha, p$beta) },
  `_characteristic_function` = function(d, t) { p <- .dl(d); hyper(list(p$alpha), list(p$alpha + p$beta), S_I * t) },
  `_moment_generating_function` = function(d, t) { p <- .dl(d); hyper(list(p$alpha), list(p$alpha + p$beta), t) })

#' @rdname continuous_rv
#' @export
BetaRV <- function(name, alpha, beta) .rvc(name, "BetaDistribution", list(alpha, beta))

.defcd("BetaNoncentralDistribution", c("alpha", "beta", "lamda"),
  set = function(d) Interval(0, 1),
  check = function(alpha, beta, lamda) {
    .value_check(alpha > 0, "Shape parameter Alpha must be positive.")
    .value_check(beta > 0, "Shape parameter Beta must be positive.")
    .value_check(lamda >= 0, "Noncentrality parameter Lambda must be positive")
  },
  pdf = function(d, x) {
    p <- .dl(d); k <- Dummy("k")
    Sum(exp(-p$lamda / 2) * (p$lamda / 2)^k * x^(p$alpha + k - 1) * (1 - x)^(p$beta - 1) /
          (sym_factorial(k) * sym_beta(p$alpha + k, p$beta)), list(k, 0, S_Infinity))
  })

#' @rdname continuous_rv
#' @export
BetaNoncentral <- function(name, alpha, beta, lamda) .rvc(name, "BetaNoncentralDistribution", list(alpha, beta, lamda))

.defcd("BetaPrimeDistribution", c("alpha", "beta"),
  check = function(alpha, beta) {
    .value_check(alpha > 0, "Shape parameter Alpha must be positive.")
    .value_check(beta > 0, "Shape parameter Beta must be positive.")
  },
  set = function(d) .R0inf(),
  pdf = function(d, x) { p <- .dl(d); x^(p$alpha - 1) * (1 + x)^(-p$alpha - p$beta) / sym_beta(p$alpha, p$beta) })

#' @rdname continuous_rv
#' @export
BetaPrime <- function(name, alpha, beta) .rvc(name, "BetaPrimeDistribution", list(alpha, beta))

.defcd("BoundedParetoDistribution", c("alpha", "left", "right"),
  set = function(d) Interval(.da(d, "left"), .da(d, "right")),
  check = function(alpha, left, right) {
    .value_check(.st_fact(alpha, "positive"), "Shape must be positive.")
    .value_check(.st_fact(left, "positive"), "Left value should be positive.")
    .value_check(right > left, "Right should be greater than left.")
  },
  pdf = function(d, x) {
    p <- .dl(d)
    num <- p$alpha * p$left^p$alpha * x^(-p$alpha - 1)
    den <- 1 - (p$left / p$right)^p$alpha
    num / den
  })

#' @rdname continuous_rv
#' @export
BoundedPareto <- function(name, alpha, left, right) .rvc(name, "BoundedParetoDistribution", list(alpha, left, right))

.defcd("CauchyDistribution", c("x0", "gamma"),
  check = function(x0, gamma) {
    .value_check(gamma > 0, "Scale parameter Gamma must be positive.")
    .value_check(.st_fact(x0, "real"), "Location parameter must be real.")
  },
  pdf = function(d, x) { p <- .dl(d); 1 / (S_Pi * p$gamma * (1 + ((x - p$x0) / p$gamma)^2)) },
  `_cdf` = function(d, x) { p <- .dl(d); 1 / S_Pi * atan((x - p$x0) / p$gamma) + S_Half },
  `_characteristic_function` = function(d, t) { p <- .dl(d); exp(p$x0 * S_I * t - p$gamma * Abs(t)) },
  `_moment_generating_function` = function(d, t) .nie("The moment generating function for the Cauchy distribution does not exist."),
  `_quantile` = function(d, q) { p <- .dl(d); p$x0 + p$gamma * tan(S_Pi * (q - S_Half)) })

#' @rdname continuous_rv
#' @export
Cauchy <- function(name, x0, gamma) .rvc(name, "CauchyDistribution", list(x0, gamma))

.defcd("ChiDistribution", c("k"),
  check = function(k) {
    .value_check(k > 0, "Number of degrees of freedom (k) must be positive.")
    .value_check(.st_fact(k, "integer"), "Number of degrees of freedom (k) must be an integer.")
  },
  set = function(d) .R0inf(),
  pdf = function(d, x) { k <- .da(d, "k"); 2^(1 - k / 2) * x^(k - 1) * exp(-x^2 / 2) / gamma(k / 2) },
  `_characteristic_function` = function(d, t) {
    k <- .da(d, "k")
    part_1 <- hyper(list(k / 2), list(S_Half), -t^2 / 2)
    part_2 <- S_I * t * sqrt(as_sym(2)) * gamma((k + 1) / 2) / gamma(k / 2)
    part_3 <- hyper(list((k + 1) / 2), list(Rational(3, 2)), -t^2 / 2)
    part_1 + part_2 * part_3
  },
  `_moment_generating_function` = function(d, t) {
    k <- .da(d, "k")
    part_1 <- hyper(list(k / 2), list(S_Half), t^2 / 2)
    part_2 <- t * sqrt(as_sym(2)) * gamma((k + 1) / 2) / gamma(k / 2)
    part_3 <- hyper(list((k + 1) / 2), list(Rational(3, 2)), t^2 / 2)
    part_1 + part_2 * part_3
  })

#' @rdname continuous_rv
#' @export
ChiRV <- function(name, k) .rvc(name, "ChiDistribution", list(k))

.defcd("ChiNoncentralDistribution", c("k", "l"),
  check = function(k, l) {
    .value_check(k > 0, "Number of degrees of freedom (k) must be positive.")
    .value_check(.st_fact(k, "integer"), "Number of degrees of freedom (k) must be an integer.")
    .value_check(l > 0, "Shift parameter Lambda must be positive.")
  },
  set = function(d) .R0inf(),
  pdf = function(d, x) {
    k <- .da(d, "k"); l <- .da(d, "l")
    exp(-(x^2 + l^2) / 2) * x^k * l / (l * x)^(k / 2) * besseli(k / 2 - 1, l * x)
  })

#' @rdname continuous_rv
#' @export
ChiNoncentral <- function(name, k, l) .rvc(name, "ChiNoncentralDistribution", list(k, l))

.defcd("ChiSquaredDistribution", c("k"),
  check = function(k) {
    .value_check(k > 0, "Number of degrees of freedom (k) must be positive.")
    .value_check(.st_fact(k, "integer"), "Number of degrees of freedom (k) must be an integer.")
  },
  set = function(d) .R0inf(),
  pdf = function(d, x) { k <- .da(d, "k"); 1 / (2^(k / 2) * gamma(k / 2)) * x^(k / 2 - 1) * exp(-x / 2) },
  `_cdf` = function(d, x) { k <- .da(d, "k"); Piecewise(list(S_One / gamma(k / 2) * lowergamma(k / 2, x / 2), x >= 0), list(0, S_true)) },
  `_characteristic_function` = function(d, t) (1 - 2 * S_I * t)^(-.da(d, "k") / 2),
  `_moment_generating_function` = function(d, t) (1 - 2 * t)^(-.da(d, "k") / 2))

#' @rdname continuous_rv
#' @export
ChiSquared <- function(name, k) .rvc(name, "ChiSquaredDistribution", list(k))

.defcd("DagumDistribution", c("p", "a", "b"),
  set = function(d) .R0inf(),
  check = function(p, a, b) {
    .value_check(p > 0, "Shape parameter p must be positive.")
    .value_check(a > 0, "Shape parameter a must be positive.")
    .value_check(b > 0, "Scale parameter b must be positive.")
  },
  pdf = function(d, x) { q <- .dl(d); q$a * q$p / x * ((x / q$b)^(q$a * q$p) / ((x / q$b)^q$a + 1)^(q$p + 1)) },
  `_cdf` = function(d, x) { q <- .dl(d); Piecewise(list((S_One + (as_sym(x) / q$b)^(-q$a))^(-q$p), x >= 0), list(S_Zero, S_true)) })

#' @rdname continuous_rv
#' @export
Dagum <- function(name, p, a, b) .rvc(name, "DagumDistribution", list(p, a, b))

.defcd("DavisDistribution", c("b", "n", "mu"),
  set = function(d) .R0inf(),
  check = function(b, n, mu) {
    .value_check(b > 0, "Scale parameter b must be positive.")
    .value_check(n > 1, "Shape parameter n must be above 1.")
    .value_check(mu > 0, "Location parameter mu must be positive.")
  },
  pdf = function(d, x) {
    p <- .dl(d)
    dividend <- p$b^p$n * (x - p$mu)^(-1 - p$n)
    divisor <- (exp(p$b / (x - p$mu)) - 1) * (gamma(p$n) * zeta(p$n))
    dividend / divisor
  })

#' @rdname continuous_rv
#' @export
Davis <- function(name, b, n, mu) .rvc(name, "DavisDistribution", list(b, n, mu))

#' @rdname continuous_rv
#' @export
Erlang <- function(name, k, l) .rvc(name, "GammaDistribution", list(k, S_One / as_sym(l)))

.defcd("ExGaussianDistribution", c("mean", "std", "rate"),
  set = function(d) .Rinf(),
  check = function(mean, std, rate) {
    .value_check(std > 0, "Standard deviation of ExGaussian must be positive.")
    .value_check(rate > 0, "Rate of ExGaussian must be positive.")
  },
  pdf = function(d, x) {
    p <- .dl(d)
    term1 <- p$rate / 2
    term2 <- exp(p$rate * (2 * p$mean + p$rate * p$std^2 - 2 * x) / 2)
    term3 <- erfc((p$mean + p$rate * p$std^2 - x) / (sqrt(as_sym(2)) * p$std))
    term1 * term2 * term3
  },
  `_cdf` = function(d, x) {
    p <- .dl(d)
    u <- p$rate * (x - p$mean)
    v <- p$rate * p$std
    G1 <- apply_lambda(rv_cdf(Normal("x", 0, v)), u)
    G2 <- apply_lambda(rv_cdf(Normal("x", v^2, v)), u)
    G1 - exp(-u + v^2 / 2 + log(G2))
  },
  `_characteristic_function` = function(d, t) {
    p <- .dl(d)
    (1 - S_I * t / p$rate)^(-1) * exp(S_I * p$mean * t - p$std^2 * t^2 / 2)
  },
  `_moment_generating_function` = function(d, t) {
    p <- .dl(d)
    (1 - t / p$rate)^(-1) * exp(p$mean * t + p$std^2 * t^2 / 2)
  })

#' @rdname continuous_rv
#' @export
ExGaussian <- function(name, mean, std, rate) .rvc(name, "ExGaussianDistribution", list(mean, std, rate))

.defcd("ExponentialDistribution", c("rate"),
  set = function(d) .R0inf(),
  check = function(rate) .value_check(rate > 0, "Rate must be positive."),
  pdf = function(d, x) { r <- .da(d, "rate"); r * exp(-r * x) },
  `_cdf` = function(d, x) { r <- .da(d, "rate"); Piecewise(list(S_One - exp(-r * x), x >= 0), list(0, S_true)) },
  `_characteristic_function` = function(d, t) { r <- .da(d, "rate"); r / (r - S_I * t) },
  `_moment_generating_function` = function(d, t) { r <- .da(d, "rate"); r / (r - t) },
  `_quantile` = function(d, p) -log(1 - p) / .da(d, "rate"))

#' @rdname continuous_rv
#' @export
Exponential <- function(name, rate) .rvc(name, "ExponentialDistribution", list(rate))

.defcd("ExponentialPowerDistribution", c("mu", "alpha", "beta"),
  set = function(d) .Rinf(),
  check = function(mu, alpha, beta) {
    .value_check(alpha > 0, "Scale parameter alpha must be positive.")
    .value_check(beta > 0, "Shape parameter beta must be positive.")
  },
  pdf = function(d, x) {
    p <- .dl(d)
    num <- p$beta * exp(-(Abs(x - p$mu) / p$alpha)^p$beta)
    den <- 2 * p$alpha * gamma(1 / p$beta)
    num / den
  },
  `_cdf` = function(d, x) {
    p <- .dl(d)
    num <- lowergamma(1 / p$beta, (Abs(x - p$mu) / p$alpha)^p$beta)
    den <- 2 * gamma(1 / p$beta)
    sign(x - p$mu) * num / den + S_Half
  })

#' @rdname continuous_rv
#' @export
ExponentialPower <- function(name, mu, alpha, beta) .rvc(name, "ExponentialPowerDistribution", list(mu, alpha, beta))

.defcd("FDistributionDistribution", c("d1", "d2"),
  set = function(d) .R0inf(),
  check = function(d1, d2) {
    .value_check(list(d1 > 0, .st_fact(d1, "integer")), "Degrees of freedom d1 must be positive integer.")
    .value_check(list(d2 > 0, .st_fact(d2, "integer")), "Degrees of freedom d2 must be positive integer.")
  },
  pdf = function(d, x) {
    p <- .dl(d)
    sqrt((p$d1 * x)^p$d1 * p$d2^p$d2 / (p$d1 * x + p$d2)^(p$d1 + p$d2)) / (x * sym_beta(p$d1 / 2, p$d2 / 2))
  },
  `_moment_generating_function` = function(d, t) .nie("The moment generating function for the F-distribution does not exist."))

#' @rdname continuous_rv
#' @export
FDistribution <- function(name, d1, d2) .rvc(name, "FDistributionDistribution", list(d1, d2))

.defcd("FisherZDistribution", c("d1", "d2"),
  set = function(d) .Rinf(),
  check = function(d1, d2) {
    .value_check(d1 > 0, "Degree of freedom d1 must be positive.")
    .value_check(d2 > 0, "Degree of freedom d2 must be positive.")
  },
  pdf = function(d, x) {
    p <- .dl(d)
    2 * p$d1^(p$d1 / 2) * p$d2^(p$d2 / 2) / sym_beta(p$d1 / 2, p$d2 / 2) * exp(p$d1 * x) / (p$d1 * exp(2 * x) + p$d2)^((p$d1 + p$d2) / 2)
  })

#' @rdname continuous_rv
#' @export
FisherZ <- function(name, d1, d2) .rvc(name, "FisherZDistribution", list(d1, d2))

.defcd("FrechetDistribution", c("a", "s", "m"),
  set = function(d) .R0inf(),
  check = function(a, s, m) {
    .value_check(a > 0, "Shape parameter alpha must be positive.")
    .value_check(s > 0, "Scale parameter s must be positive.")
  },
  pdf = function(d, x) { p <- .dl(d); p$a / p$s * ((x - p$m) / p$s)^(-1 - p$a) * exp(-((x - p$m) / p$s)^(-p$a)) },
  `_cdf` = function(d, x) { p <- .dl(d); Piecewise(list(exp(-((x - p$m) / p$s)^(-p$a)), x >= p$m), list(S_Zero, S_true)) })

#' @rdname continuous_rv
#' @export
Frechet <- function(name, a, s = 1, m = 0) .rvc(name, "FrechetDistribution", list(a, s, m))

.defcd("GammaDistribution", c("k", "theta"),
  set = function(d) .R0inf(),
  check = function(k, theta) {
    .value_check(k > 0, "k must be positive")
    .value_check(theta > 0, "Theta must be positive")
  },
  pdf = function(d, x) { p <- .dl(d); x^(p$k - 1) * exp(-x / p$theta) / (gamma(p$k) * p$theta^p$k) },
  `_cdf` = function(d, x) { p <- .dl(d); Piecewise(list(lowergamma(p$k, as_sym(x) / p$theta) / gamma(p$k), x > 0), list(S_Zero, S_true)) },
  `_characteristic_function` = function(d, t) { p <- .dl(d); (1 - p$theta * S_I * t)^(-p$k) },
  `_moment_generating_function` = function(d, t) { p <- .dl(d); (1 - p$theta * t)^(-p$k) })

#' @rdname continuous_rv
#' @export
GammaRV <- function(name, k, theta) .rvc(name, "GammaDistribution", list(k, theta))

.defcd("GammaInverseDistribution", c("a", "b"),
  set = function(d) .R0inf(),
  check = function(a, b) {
    .value_check(a > 0, "alpha must be positive")
    .value_check(b > 0, "beta must be positive")
  },
  pdf = function(d, x) { p <- .dl(d); p$b^p$a / gamma(p$a) * x^(-p$a - 1) * exp(-p$b / x) },
  `_cdf` = function(d, x) { p <- .dl(d); Piecewise(list(uppergamma(p$a, p$b / x) / gamma(p$a), x > 0), list(S_Zero, S_true)) },
  `_characteristic_function` = function(d, t) {
    p <- .dl(d)
    2 * (-S_I * p$b * t)^(p$a / 2) * besselk(p$a, sqrt(-4 * S_I * p$b * t)) / gamma(p$a)
  },
  `_moment_generating_function` = function(d, t) .nie("The moment generating function for the gamma inverse distribution does not exist."))

#' @rdname continuous_rv
#' @export
GammaInverse <- function(name, a, b) .rvc(name, "GammaInverseDistribution", list(a, b))

.gumbel_pw <- function(d, fmin, fmax) {
  m <- .da(d, "minimum")
  Piecewise(list(fmin, m), list(fmax, as_sym(!.truth(m))))
}

.defcd("GumbelDistribution", c("beta", "mu", "minimum"),
  set = function(d) .Rinf(),
  check = function(beta, mu, minimum) .value_check(beta > 0, "Scale parameter beta must be positive."),
  pdf = function(d, x) {
    p <- .dl(d); z <- (x - p$mu) / p$beta
    .gumbel_pw(d, 1 / p$beta * exp(z - exp(z)), 1 / p$beta * exp(-z - exp(-z)))
  },
  `_cdf` = function(d, x) {
    p <- .dl(d); z <- (x - p$mu) / p$beta
    .gumbel_pw(d, 1 - exp(-exp(z)), exp(-exp(-z)))
  },
  `_characteristic_function` = function(d, t) {
    p <- .dl(d)
    .gumbel_pw(d, gamma(1 + S_I * p$beta * t) * exp(S_I * p$mu * t), gamma(1 - S_I * p$beta * t) * exp(S_I * p$mu * t))
  },
  `_moment_generating_function` = function(d, t) {
    p <- .dl(d)
    .gumbel_pw(d, gamma(1 + p$beta * t) * exp(p$mu * t), gamma(1 - p$beta * t) * exp(p$mu * t))
  })

#' @rdname continuous_rv
#' @export
Gumbel <- function(name, beta, mu, minimum = FALSE) .rvc(name, "GumbelDistribution", list(beta, mu, minimum))

.defcd("GompertzDistribution", c("b", "eta"),
  set = function(d) .R0inf(),
  check = function(b, eta) {
    .value_check(b > 0, "b must be positive")
    .value_check(eta > 0, "eta must be positive")
  },
  pdf = function(d, x) { p <- .dl(d); p$b * p$eta * exp(p$b * x) * exp(p$eta) * exp(-p$eta * exp(p$b * x)) },
  `_cdf` = function(d, x) { p <- .dl(d); 1 - exp(p$eta) * exp(-p$eta * exp(p$b * x)) },
  `_moment_generating_function` = function(d, t) { p <- .dl(d); p$eta * exp(p$eta) * expint(t / p$b, p$eta) })

#' @rdname continuous_rv
#' @export
Gompertz <- function(name, b, eta) .rvc(name, "GompertzDistribution", list(b, eta))

.defcd("KumaraswamyDistribution", c("a", "b"),
  set = function(d) .R0inf(),
  check = function(a, b) {
    .value_check(a > 0, "a must be positive")
    .value_check(b > 0, "b must be positive")
  },
  pdf = function(d, x) { p <- .dl(d); p$a * p$b * x^(p$a - 1) * (1 - x^p$a)^(p$b - 1) },
  `_cdf` = function(d, x) {
    p <- .dl(d)
    Piecewise(list(S_Zero, x < S_Zero), list(1 - (1 - x^p$a)^p$b, x <= S_One), list(S_One, S_true))
  })

#' @rdname continuous_rv
#' @export
Kumaraswamy <- function(name, a, b) .rvc(name, "KumaraswamyDistribution", list(a, b))

.defcd("LaplaceDistribution", c("mu", "b"),
  set = function(d) .Rinf(),
  check = function(mu, b) {
    .value_check(b > 0, "Scale parameter b must be positive.")
    .value_check(.st_fact(mu, "real"), "Location parameter mu should be real")
  },
  pdf = function(d, x) { p <- .dl(d); 1 / (2 * p$b) * exp(-Abs(x - p$mu) / p$b) },
  `_cdf` = function(d, x) {
    p <- .dl(d)
    Piecewise(list(S_Half * exp((x - p$mu) / p$b), x < p$mu), list(S_One - S_Half * exp(-(x - p$mu) / p$b), x >= p$mu))
  },
  `_characteristic_function` = function(d, t) { p <- .dl(d); exp(p$mu * S_I * t) / (1 + p$b^2 * t^2) },
  `_moment_generating_function` = function(d, t) { p <- .dl(d); exp(p$mu * t) / (1 - p$b^2 * t^2) })

.is_vecarg <- function(a) .is_mat(a) || (is.list(a) && !is_sym(a))

#' @rdname continuous_rv
#' @export
Laplace <- function(name, mu, b) {
  if (.is_vecarg(mu) && .is_vecarg(b)) return(MultivariateLaplace(name, mu, b))
  .rvc(name, "LaplaceDistribution", list(mu, b))
}

.defcd("LevyDistribution", c("mu", "c"),
  set = function(d) Interval(.da(d, "mu"), S_Infinity),
  check = function(mu, c) {
    .value_check(c > 0, "c (scale parameter) must be positive")
    .value_check(.st_fact(mu, "real"), "mu (location parameter) must be real")
  },
  pdf = function(d, x) {
    p <- .dl(d)
    sqrt(p$c / (2 * S_Pi)) * exp(-p$c / (2 * (x - p$mu))) / (x - p$mu)^(S_One + S_Half)
  },
  `_cdf` = function(d, x) { p <- .dl(d); erfc(sqrt(p$c / (2 * (x - p$mu)))) },
  `_characteristic_function` = function(d, t) { p <- .dl(d); exp(S_I * p$mu * t - sqrt(-2 * S_I * p$c * t)) },
  `_moment_generating_function` = function(d, t) .nie("The moment generating function of Levy distribution does not exist."))

#' @rdname continuous_rv
#' @export
Levy <- function(name, mu, c) .rvc(name, "LevyDistribution", list(mu, c))

.defcd("LogCauchyDistribution", c("mu", "sigma"),
  set = function(d) Interval(0, S_Infinity, TRUE, TRUE),
  check = function(mu, sigma) {
    .value_check(!identical(.fuzzy_bool(sigma > 0), FALSE), "Scale parameter Gamma must be positive.")
    .value_check(!identical(.st_fact(mu, "real"), FALSE), "Location parameter must be real.")
  },
  pdf = function(d, x) { p <- .dl(d); 1 / (x * S_Pi) * (p$sigma / ((log(x) - p$mu)^2 + p$sigma^2)) },
  `_cdf` = function(d, x) { p <- .dl(d); 1 / S_Pi * atan((log(x) - p$mu) / p$sigma) + S_Half },
  `_characteristic_function` = function(d, t) .nie("The characteristic function for the Log-Cauchy distribution does not exist."),
  `_moment_generating_function` = function(d, t) .nie("The moment generating function for the Log-Cauchy distribution does not exist."))

#' @rdname continuous_rv
#' @export
LogCauchy <- function(name, mu, sigma) .rvc(name, "LogCauchyDistribution", list(mu, sigma))

.defcd("LogisticDistribution", c("mu", "s"),
  set = function(d) .Rinf(),
  check = function(mu, s) .value_check(s > 0, "Scale parameter s must be positive."),
  pdf = function(d, x) { p <- .dl(d); exp(-(x - p$mu) / p$s) / (p$s * (1 + exp(-(x - p$mu) / p$s))^2) },
  `_cdf` = function(d, x) { p <- .dl(d); S_One / (1 + exp(-(x - p$mu) / p$s)) },
  `_characteristic_function` = function(d, t) {
    p <- .dl(d)
    Piecewise(list(exp(S_I * t * p$mu) * S_Pi * p$s * t / sinh(S_Pi * p$s * t), Ne(t, 0)), list(S_One, S_true))
  },
  `_moment_generating_function` = function(d, t) { p <- .dl(d); exp(p$mu * t) * sym_beta(1 - p$s * t, 1 + p$s * t) },
  `_quantile` = function(d, q) { p <- .dl(d); p$mu - p$s * log(-S_One + S_One / q) })

#' @rdname continuous_rv
#' @export
Logistic <- function(name, mu, s) .rvc(name, "LogisticDistribution", list(mu, s))

.defcd("LogLogisticDistribution", c("alpha", "beta"),
  set = function(d) .R0inf(),
  check = function(alpha, beta) {
    .value_check(alpha > 0, "Scale parameter Alpha must be positive.")
    .value_check(beta > 0, "Shape parameter Beta must be positive.")
  },
  pdf = function(d, x) { p <- .dl(d); p$beta / p$alpha * (x / p$alpha)^(p$beta - 1) / (1 + (x / p$alpha)^p$beta)^2 },
  `_cdf` = function(d, x) { p <- .dl(d); 1 / (1 + (x / p$alpha)^(-p$beta)) },
  `_quantile` = function(d, q) { p <- .dl(d); p$alpha * (q / (1 - q))^(1 / p$beta) },
  expectation = function(d, expr, var, ...) {
    a <- d$args[[1L]]; b <- d$args[[2L]]
    Piecewise(list(S_NaN, b <= 1), list(S_Pi * a / (b * sin(S_Pi / b)), S_true))
  })

#' @rdname continuous_rv
#' @export
LogLogistic <- function(name, alpha, beta) .rvc(name, "LogLogisticDistribution", list(alpha, beta))

.defcd("LogitNormalDistribution", c("mu", "s"),
  set = function(d) Interval(0, 1, TRUE, TRUE),
  check = function(mu, s) {
    .value_check(.pyand(!identical(.st_fact(s^2, "real"), FALSE), s^2 > 0), "Squared scale parameter s must be positive.")
    .value_check(!identical(.st_fact(mu, "real"), FALSE), "Location parameter must be real")
  },
  pdf = function(d, x) {
    p <- .dl(d); lg <- log(x / (1 - x))
    exp(-(lg - p$mu)^2 / (2 * p$s^2)) * (S_One / sqrt(2 * S_Pi * p$s^2)) * (1 / (x * (1 - x)))
  },
  `_cdf` = function(d, x) {
    p <- .dl(d); lg <- log(x / (1 - x))
    S_One / 2 * (1 + erf((lg - p$mu) / sqrt(2 * p$s^2)))
  })

#' @rdname continuous_rv
#' @export
LogitNormal <- function(name, mu, s) .rvc(name, "LogitNormalDistribution", list(mu, s))

.defcd("LogNormalDistribution", c("mean", "std"),
  set = function(d) .R0inf(),
  check = function(mean, std) .value_check(std > 0, "Parameter std must be positive."),
  pdf = function(d, x) { p <- .dl(d); exp(-(log(x) - p$mean)^2 / (2 * p$std^2)) / (x * sqrt(2 * S_Pi) * p$std) },
  `_cdf` = function(d, x) {
    p <- .dl(d)
    Piecewise(list(S_Half + S_Half * erf((log(x) - p$mean) / sqrt(as_sym(2)) / p$std), x > 0), list(S_Zero, S_true))
  },
  `_moment_generating_function` = function(d, t) .nie("Moment generating function of the log-normal distribution is not defined."))

#' @rdname continuous_rv
#' @export
LogNormal <- function(name, mean, std) .rvc(name, "LogNormalDistribution", list(mean, std))

.defcd("LomaxDistribution", c("alpha", "lamda"),
  set = function(d) .R0inf(),
  check = function(alpha, lamda) {
    .value_check(.st_fact(alpha, "real"), "Shape parameter should be real.")
    .value_check(.st_fact(lamda, "real"), "Scale parameter should be real.")
    .value_check(.st_fact(alpha, "positive"), "Shape parameter should be positive.")
    .value_check(.st_fact(lamda, "positive"), "Scale parameter should be positive.")
  },
  pdf = function(d, x) { p <- .dl(d); p$alpha / p$lamda * (S_One + x / p$lamda)^(-p$alpha - 1) })

#' @rdname continuous_rv
#' @export
Lomax <- function(name, alpha, lamda) .rvc(name, "LomaxDistribution", list(alpha, lamda))

.defcd("MaxwellDistribution", c("a"),
  set = function(d) .R0inf(),
  check = function(a) .value_check(a > 0, "Parameter a must be positive."),
  pdf = function(d, x) { a <- .da(d, "a"); sqrt(2 / S_Pi) * x^2 * exp(-x^2 / (2 * a^2)) / a^3 },
  `_cdf` = function(d, x) {
    a <- .da(d, "a")
    erf(sqrt(as_sym(2)) * x / (2 * a)) - sqrt(as_sym(2)) * x * exp(-x^2 / (2 * a^2)) / (sqrt(S_Pi) * a)
  })

#' @rdname continuous_rv
#' @export
Maxwell <- function(name, a) .rvc(name, "MaxwellDistribution", list(a))

.defcd("MoyalDistribution", c("mu", "sigma"),
  check = function(mu, sigma) {
    .value_check(.st_fact(mu, "real"), "Location parameter must be real.")
    .value_check(.pyand(.st_fact(sigma, "real"), sigma > 0), "Scale parameter must be real        and positive.")
  },
  pdf = function(d, x) {
    p <- .dl(d)
    num <- exp(-(exp(-(x - p$mu) / p$sigma) + (x - p$mu) / p$sigma) / 2)
    den <- sqrt(2 * S_Pi) * p$sigma
    num / den
  },
  `_characteristic_function` = function(d, t) {
    p <- .dl(d)
    exp(S_I * t * p$mu) * (2^(-S_I * p$sigma * t) * gamma(Rational(1, 2) - S_I * t * p$sigma)) / sqrt(S_Pi)
  },
  `_moment_generating_function` = function(d, t) {
    p <- .dl(d)
    exp(t * p$mu) * (2^(-1 * p$sigma * t) * gamma(Rational(1, 2) - t * p$sigma)) / sqrt(S_Pi)
  })

#' @rdname continuous_rv
#' @export
Moyal <- function(name, mu, sigma) .rvc(name, "MoyalDistribution", list(mu, sigma))

.defcd("NakagamiDistribution", c("mu", "omega"),
  set = function(d) .R0inf(),
  check = function(mu, omega) {
    .value_check(mu >= S_Half, "Shape parameter mu must be greater than equal to 1/2.")
    .value_check(omega > 0, "Spread parameter omega must be positive.")
  },
  pdf = function(d, x) {
    p <- .dl(d)
    2 * p$mu^p$mu / (gamma(p$mu) * p$omega^p$mu) * x^(2 * p$mu - 1) * exp(-p$mu / p$omega * x^2)
  },
  `_cdf` = function(d, x) {
    p <- .dl(d)
    Piecewise(list(lowergamma(p$mu, p$mu / p$omega * x^2) / gamma(p$mu), x > 0), list(S_Zero, S_true))
  })

#' @rdname continuous_rv
#' @export
Nakagami <- function(name, mu, omega) .rvc(name, "NakagamiDistribution", list(mu, omega))

.defcd("NormalDistribution", c("mean", "std"),
  check = function(mean, std) .value_check(std > 0, "Standard deviation must be positive"),
  pdf = function(d, x) { p <- .dl(d); exp(-(x - p$mean)^2 / (2 * p$std^2)) / (sqrt(2 * S_Pi) * p$std) },
  `_cdf` = function(d, x) { p <- .dl(d); erf(sqrt(as_sym(2)) * (-p$mean + x) / (2 * p$std)) / 2 + S_Half },
  `_characteristic_function` = function(d, t) { p <- .dl(d); exp(S_I * p$mean * t - p$std^2 * t^2 / 2) },
  `_moment_generating_function` = function(d, t) { p <- .dl(d); exp(p$mean * t + p$std^2 * t^2 / 2) },
  `_quantile` = function(d, q) { p <- .dl(d); p$mean + p$std * sqrt(as_sym(2)) * erfinv(2 * q - 1) })

#' @rdname continuous_rv
#' @export
Normal <- function(name, mean, std) {
  if ((is.list(mean) && !is_sym(mean)) || .is_mat(mean) && (is.list(std) && !is_sym(std)) || .is_mat(std))
    return(MultivariateNormal(name, mean, std))
  .rvc(name, "NormalDistribution", list(mean, std))
}

.defcd("GaussianInverseDistribution", c("mean", "shape"),
  set = function(d) .R0inf(),
  check = function(mean, shape) {
    .value_check(shape > 0, "Shape parameter must be positive")
    .value_check(mean > 0, "Mean must be positive")
  },
  pdf = function(d, x) {
    p <- .dl(d)
    exp(-p$shape * (x - p$mean)^2 / (2 * x * p$mean^2)) * sqrt(p$shape / (2 * S_Pi * x^3))
  },
  `_cdf` = function(d, x) {
    p <- .dl(d)
    sn <- rv_cdf(Normal("x", 0, 1))
    first <- apply_lambda(sn, sqrt(p$shape / x) * (x / p$mean - S_One))
    second <- exp(2 * p$shape / p$mean) * apply_lambda(sn, -sqrt(p$shape / x) * (x / p$mean + S_One))
    first + second
  },
  `_characteristic_function` = function(d, t) {
    p <- .dl(d)
    exp(p$shape / p$mean * (1 - sqrt(1 - 2 * p$mean^2 * S_I * t / p$shape)))
  },
  `_moment_generating_function` = function(d, t) {
    p <- .dl(d)
    exp(p$shape / p$mean * (1 - sqrt(1 - 2 * p$mean^2 * t / p$shape)))
  })

#' @rdname continuous_rv
#' @export
GaussianInverse <- function(name, mean, shape) .rvc(name, "GaussianInverseDistribution", list(mean, shape))

#' @rdname continuous_rv
#' @export
Wald <- GaussianInverse

.defcd("ParetoDistribution", c("xm", "alpha"),
  set = function(d) Interval(.da(d, "xm"), S_Infinity),
  check = function(xm, alpha) {
    .value_check(xm > 0, "Xm must be positive")
    .value_check(alpha > 0, "Alpha must be positive")
  },
  pdf = function(d, x) { p <- .dl(d); p$alpha * p$xm^p$alpha / x^(p$alpha + 1) },
  `_cdf` = function(d, x) { p <- .dl(d); Piecewise(list(S_One - p$xm^p$alpha / x^p$alpha, x >= p$xm), list(0, S_true)) },
  `_moment_generating_function` = function(d, t) { p <- .dl(d); p$alpha * (-p$xm * t)^p$alpha * uppergamma(-p$alpha, -p$xm * t) },
  `_characteristic_function` = function(d, t) { p <- .dl(d); p$alpha * (-S_I * p$xm * t)^p$alpha * uppergamma(-p$alpha, -S_I * p$xm * t) })

#' @rdname continuous_rv
#' @export
Pareto <- function(name, xm, alpha) .rvc(name, "ParetoDistribution", list(xm, alpha))

.defcd("PowerFunctionDistribution", c("alpha", "a", "b"),
  set = function(d) Interval(.da(d, "a"), .da(d, "b")),
  check = function(alpha, a, b) {
    .value_check(.st_fact(a, "real"), "Continuous Boundary parameter should be real.")
    .value_check(.st_fact(b, "real"), "Continuous Boundary parameter should be real.")
    .value_check(a < b, " 'a' the left Boundary must be smaller than 'b' the right Boundary.")
    .value_check(.st_fact(alpha, "positive"), "Continuous Shape parameter should be positive.")
  },
  pdf = function(d, x) { p <- .dl(d); p$alpha * (x - p$a)^(p$alpha - 1) / (p$b - p$a)^p$alpha })

#' @rdname continuous_rv
#' @export
PowerFunction <- function(name, alpha, a, b) .rvc(name, "PowerFunctionDistribution", list(alpha, a, b))

.defcd("QuadraticUDistribution", c("a", "b"),
  set = function(d) Interval(.da(d, "a"), .da(d, "b")),
  check = function(a, b) .value_check(b > a, paste0("Parameter b must be in range (", .rv_pystr(a), ", oo).")),
  pdf = function(d, x) {
    p <- .dl(d)
    alpha <- 12 / (p$b - p$a)^3
    beta <- (p$a + p$b) / 2
    Piecewise(list(alpha * (x - beta)^2, And(p$a <= x, x <= p$b)), list(S_Zero, S_true))
  },
  `_moment_generating_function` = function(d, t) {
    a <- .da(d, "a"); b <- .da(d, "b")
    -3 * (exp(a * t) * (4 + (a^2 + 2 * a * (-2 + b) + b^2) * t) - exp(b * t) * (4 + (-4 * b + (a + b)^2) * t)) / ((a - b)^3 * t^2)
  },
  `_characteristic_function` = function(d, t) {
    a <- .da(d, "a"); b <- .da(d, "b")
    -3 * S_I * (exp(S_I * a * t * exp(S_I * b * t)) * (4 * S_I - (-4 * b + (a + b)^2) * t)) / ((a - b)^3 * t^2)
  })

#' @rdname continuous_rv
#' @export
QuadraticU <- function(name, a, b) .rvc(name, "QuadraticUDistribution", list(a, b))

.defcd("RaisedCosineDistribution", c("mu", "s"),
  set = function(d) { p <- .dl(d); Interval(p$mu - p$s, p$mu + p$s) },
  check = function(mu, s) .value_check(s > 0, "s must be positive"),
  pdf = function(d, x) {
    p <- .dl(d)
    Piecewise(list((1 + cos(S_Pi * (x - p$mu) / p$s)) / (2 * p$s), And(p$mu - p$s <= x, x <= p$mu + p$s)), list(S_Zero, S_true))
  },
  `_characteristic_function` = function(d, t) {
    p <- .dl(d)
    Piecewise(list(exp(-S_I * S_Pi * p$mu / p$s) / 2, Eq(t, -S_Pi / p$s)),
              list(exp(S_I * S_Pi * p$mu / p$s) / 2, Eq(t, S_Pi / p$s)),
              list(S_Pi^2 * sin(p$s * t) * exp(S_I * p$mu * t) / (p$s * t * (S_Pi^2 - p$s^2 * t^2)), S_true))
  },
  `_moment_generating_function` = function(d, t) {
    p <- .dl(d)
    S_Pi^2 * sinh(p$s * t) * exp(p$mu * t) / (p$s * t * (S_Pi^2 + p$s^2 * t^2))
  })

#' @rdname continuous_rv
#' @export
RaisedCosine <- function(name, mu, s) .rvc(name, "RaisedCosineDistribution", list(mu, s))

.defcd("RayleighDistribution", c("sigma"),
  set = function(d) .R0inf(),
  check = function(sigma) .value_check(sigma > 0, "Scale parameter sigma must be positive."),
  pdf = function(d, x) { s <- .da(d, "sigma"); x / s^2 * exp(-x^2 / (2 * s^2)) },
  `_cdf` = function(d, x) { s <- .da(d, "sigma"); 1 - exp(-(x^2 / (2 * s^2))) },
  `_characteristic_function` = function(d, t) {
    s <- .da(d, "sigma")
    1 - s * t * exp(-s^2 * t^2 / 2) * sqrt(S_Pi / 2) * (erfi(s * t / sqrt(as_sym(2))) - S_I)
  },
  `_moment_generating_function` = function(d, t) {
    s <- .da(d, "sigma")
    1 + s * t * exp(s^2 * t^2 / 2) * sqrt(S_Pi / 2) * (erf(s * t / sqrt(as_sym(2))) + 1)
  })

#' @rdname continuous_rv
#' @export
Rayleigh <- function(name, sigma) .rvc(name, "RayleighDistribution", list(sigma))

.defcd("ReciprocalDistribution", c("a", "b"),
  set = function(d) Interval(.da(d, "a"), .da(d, "b")),
  check = function(a, b) {
    .value_check(a > 0, paste0("Parameter > 0. a = ", .rv_pystr(a)))
    .value_check(a < b, paste0("Parameter b must be in range (", .rv_pystr(a), ", +oo]. b = ", .rv_pystr(b)))
  },
  pdf = function(d, x) { p <- .dl(d); 1 / (x * (log(p$b) - log(p$a))) })

#' @rdname continuous_rv
#' @export
Reciprocal <- function(name, a, b) .rvc(name, "ReciprocalDistribution", list(a, b))

.defcd("ShiftedGompertzDistribution", c("b", "eta"),
  set = function(d) .R0inf(),
  check = function(b, eta) {
    .value_check(b > 0, "b must be positive")
    .value_check(eta > 0, "eta must be positive")
  },
  pdf = function(d, x) {
    p <- .dl(d)
    p$b * exp(-p$b * x) * exp(-p$eta * exp(-p$b * x)) * (1 + p$eta * (1 - exp(-p$b * x)))
  })

#' @rdname continuous_rv
#' @export
ShiftedGompertz <- function(name, b, eta) .rvc(name, "ShiftedGompertzDistribution", list(b, eta))

.defcd("SkewNormalDistribution", c("xi", "omega", "alpha"),
  check = function(xi, omega, alpha) {
    .value_check(omega > 0, "Scale parameter must be positive")
    .value_check(!identical(.st_fact(alpha, "real"), FALSE), "Shape parameter alpha must be real.")
  },
  pdf = function(d, x) {
    p <- .dl(d); z <- (x - p$xi) / p$omega
    exp(-z^2 / 2) * (1 + erf(p$alpha * z / sqrt(as_sym(2)))) / (sqrt(2 * S_Pi) * p$omega)
  },
  `_cdf` = function(d, x) {
    p <- .dl(d); z <- (x - p$xi) / p$omega
    erf(sqrt(as_sym(2)) * z / 2) / 2 + S_Half - 2 * owens_t(z, p$alpha)
  })

#' @rdname continuous_rv
#' @export
SkewNormal <- function(name, xi, omega, alpha) .rvc(name, "SkewNormalDistribution", list(xi, omega, alpha))

.defcd("StudentTDistribution", c("nu"),
  set = function(d) .Rinf(),
  check = function(nu) .value_check(nu > 0, "Degrees of freedom nu must be positive."),
  pdf = function(d, x) { nu <- .da(d, "nu"); 1 / (sqrt(nu) * sym_beta(S_Half, nu / 2)) * (1 + x^2 / nu)^(-(nu + 1) / 2) },
  `_cdf` = function(d, x) {
    nu <- .da(d, "nu")
    S_Half + x * gamma((nu + 1) / 2) * hyper(list(S_Half, (nu + 1) / 2), list(Rational(3, 2)), -x^2 / nu) / (sqrt(S_Pi * nu) * gamma(nu / 2))
  },
  `_moment_generating_function` = function(d, t) .nie("The moment generating function for the Student-T distribution is undefined."))

#' @rdname continuous_rv
#' @export
StudentT <- function(name, nu) .rvc(name, "StudentTDistribution", list(nu))

.defcd("TrapezoidalDistribution", c("a", "b", "c", "d"),
  set = function(d) Interval(.da(d, "a"), .da(d, "d")),
  check = function(a, b, c, d) {
    .value_check(a < d, paste0("Lower bound parameter a < ", .rv_pystr(d), ". a = ", .rv_pystr(a)))
    .value_check(list(a <= b, b < c), paste0("Level start parameter b must be in range [", .rv_pystr(a), ", ", .rv_pystr(c), "). b = ", .rv_pystr(b)))
    .value_check(list(b < c, c <= d), paste0("Level end parameter c must be in range (", .rv_pystr(b), ", ", .rv_pystr(d), "]. c = ", .rv_pystr(c)))
    .value_check(d >= c, paste0("Upper bound parameter d > ", .rv_pystr(c), ". d = ", .rv_pystr(d)))
  },
  pdf = function(dd, x) {
    p <- .dl(dd); a <- p$a; b <- p$b; c <- p$c; d <- p$d
    Piecewise(list(2 * (x - a) / ((b - a) * (d + c - a - b)), And(a <= x, x < b)),
              list(2 / (d + c - a - b), And(b <= x, x < c)),
              list(2 * (d - x) / ((d - c) * (d + c - a - b)), And(c <= x, x <= d)),
              list(S_Zero, S_true))
  })

#' @rdname continuous_rv
#' @export
Trapezoidal <- function(name, a, b, c, d) .rvc(name, "TrapezoidalDistribution", list(a, b, c, d))

.defcd("TriangularDistribution", c("a", "b", "c"),
  set = function(d) Interval(.da(d, "a"), .da(d, "b")),
  check = function(a, b, c) {
    .value_check(b > a, paste0("Parameter b > ", .rv_pystr(a), ". b = ", .rv_pystr(b)))
    .value_check(list(a <= c, c <= b), paste0("Parameter c must be in range [", .rv_pystr(a), ", ", .rv_pystr(b), "]. c = ", .rv_pystr(c)))
  },
  pdf = function(d, x) {
    p <- .dl(d); a <- p$a; b <- p$b; c <- p$c
    Piecewise(list(2 * (x - a) / ((b - a) * (c - a)), And(a <= x, x < c)),
              list(2 / (b - a), Eq(x, c)),
              list(2 * (b - x) / ((b - a) * (b - c)), And(c < x, x <= b)),
              list(S_Zero, S_true))
  },
  `_characteristic_function` = function(d, t) {
    p <- .dl(d); a <- p$a; b <- p$b; c <- p$c
    -2 * ((b - c) * exp(S_I * a * t) - (b - a) * exp(S_I * c * t) + (c - a) * exp(S_I * b * t)) / ((b - a) * (c - a) * (b - c) * t^2)
  },
  `_moment_generating_function` = function(d, t) {
    p <- .dl(d); a <- p$a; b <- p$b; c <- p$c
    2 * ((b - c) * exp(a * t) - (b - a) * exp(c * t) + (c - a) * exp(b * t)) / ((b - a) * (c - a) * (b - c) * t^2)
  })

#' @rdname continuous_rv
#' @export
Triangular <- function(name, a, b, c) .rvc(name, "TriangularDistribution", list(a, b, c))

.defcd("UniformDistribution", c("left", "right"),
  set = function(d) Interval(.da(d, "left"), .da(d, "right")),
  check = function(left, right) .value_check(left < right, "Lower limit should be less than Upper limit."),
  pdf = function(d, x) {
    p <- .dl(d)
    Piecewise(list(S_One / (p$right - p$left), And(p$left <= x, x <= p$right)), list(S_Zero, S_true))
  },
  `_cdf` = function(d, x) {
    p <- .dl(d)
    Piecewise(list(S_Zero, x < p$left), list((x - p$left) / (p$right - p$left), x <= p$right), list(S_One, S_true))
  },
  `_characteristic_function` = function(d, t) {
    p <- .dl(d)
    Piecewise(list((exp(S_I * t * p$right) - exp(S_I * t * p$left)) / (S_I * t * (p$right - p$left)), Ne(t, 0)), list(S_One, S_true))
  },
  `_moment_generating_function` = function(d, t) {
    p <- .dl(d)
    Piecewise(list((exp(t * p$right) - exp(t * p$left)) / (t * (p$right - p$left)), Ne(t, 0)), list(S_One, S_true))
  },
  expectation = function(d, expr, var, ...) {
    a <- list(...); a$evaluate <- TRUE
    result <- do.call(.STM[["SingleContinuousDistribution"]]$expectation, c(list(d, expr, var), a))
    p <- .dl(d)
    subs(result, list(list(Max(p$left, p$right), p$right), list(Min(p$left, p$right), p$left)))
  })

#' @rdname continuous_rv
#' @export
Uniform <- function(name, left, right) .rvc(name, "UniformDistribution", list(left, right))

.defcd("UniformSumDistribution", c("n"),
  set = function(d) Interval(0, .da(d, "n")),
  check = function(n) .value_check(list(n > 0, .st_fact(n, "integer")), "Parameter n must be positive integer."),
  pdf = function(d, x) {
    n <- .da(d, "n"); k <- Dummy("k")
    1 / sym_factorial(n - 1) * Sum((-1)^k * sym_binomial(n, k) * (x - k)^(n - 1), list(k, 0, floor(x)))
  },
  `_cdf` = function(d, x) {
    n <- .da(d, "n"); k <- Dummy("k")
    Piecewise(list(S_Zero, x < 0),
              list(1 / sym_factorial(n) * Sum((-1)^k * sym_binomial(n, k) * (x - k)^n, list(k, 0, floor(x))), x <= n),
              list(S_One, S_true))
  },
  `_characteristic_function` = function(d, t) ((exp(S_I * t) - 1) / (S_I * t))^.da(d, "n"),
  `_moment_generating_function` = function(d, t) ((exp(t) - 1) / t)^.da(d, "n"))

#' @rdname continuous_rv
#' @export
UniformSum <- function(name, n) .rvc(name, "UniformSumDistribution", list(n))

.defcd("VonMisesDistribution", c("mu", "k"),
  set = function(d) Interval(0, 2 * S_Pi),
  check = function(mu, k) .value_check(k > 0, "k must be positive"),
  pdf = function(d, x) { p <- .dl(d); exp(p$k * cos(x - p$mu)) / (2 * S_Pi * besseli(0, p$k)) })

#' @rdname continuous_rv
#' @export
VonMises <- function(name, mu, k) .rvc(name, "VonMisesDistribution", list(mu, k))

.defcd("WeibullDistribution", c("alpha", "beta"),
  set = function(d) .R0inf(),
  check = function(alpha, beta) {
    .value_check(alpha > 0, "Alpha must be positive")
    .value_check(beta > 0, "Beta must be positive")
  },
  pdf = function(d, x) { p <- .dl(d); p$beta * (x / p$alpha)^(p$beta - 1) * exp(-(x / p$alpha)^p$beta) / p$alpha })

#' @rdname continuous_rv
#' @export
Weibull <- function(name, alpha, beta) .rvc(name, "WeibullDistribution", list(alpha, beta))

.defcd("WignerSemicircleDistribution", c("R"),
  set = function(d) { R <- .da(d, "R"); Interval(-R, R) },
  check = function(R) .value_check(R > 0, "Radius R must be positive."),
  pdf = function(d, x) { R <- .da(d, "R"); 2 / (S_Pi * R^2) * sqrt(R^2 - x^2) },
  `_characteristic_function` = function(d, t) {
    R <- .da(d, "R")
    Piecewise(list(2 * besselj(1, R * t) / (R * t), Ne(t, 0)), list(S_One, S_true))
  },
  `_moment_generating_function` = function(d, t) {
    R <- .da(d, "R")
    Piecewise(list(2 * besseli(1, R * t) / (R * t), Ne(t, 0)), list(S_One, S_true))
  })

#' @rdname continuous_rv
#' @export
WignerSemicircle <- function(name, R) .rvc(name, "WignerSemicircleDistribution", list(R))
