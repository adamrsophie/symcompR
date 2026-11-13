# Statistics and random variables

`sympy.stats` computes with random variables symbolically: densities,
probabilities, expectations, moments and conditional distributions, for
continuous, discrete, finite, joint and matrix valued random variables,
random matrices and stochastic processes. The vignette on other modules
shows the most common distributions
([`Normal()`](https://adamrsophie.github.io/symcompR/reference/continuous_rv.md),
[`Exponential()`](https://adamrsophie.github.io/symcompR/reference/continuous_rv.md),
[`Poisson()`](https://adamrsophie.github.io/symcompR/reference/discrete_rv.md),
[`Die()`](https://adamrsophie.github.io/symcompR/reference/finite_rv.md)
and others); this one covers the rest.

SymPy’s [`density()`](https://rdrr.io/r/stats/density.html) returns a
function. In R, `rv_density(X, z)` evaluates it at `z`, and
[`apply_lambda()`](https://adamrsophie.github.io/symcompR/reference/Lambda.md)
evaluates the function objects returned by `rv_density(X)`,
[`rv_quantile()`](https://adamrsophie.github.io/symcompR/reference/E.md),
[`characteristic_function()`](https://adamrsophie.github.io/symcompR/reference/E.md)
and the other functions of several arguments.

``` r

library(symcompR)
z <- Symbol("z"); x1 <- Symbol("x1"); x2 <- Symbol("x2")
X <- Normal("X", 0, 1); Y <- Normal("Y", 0, 1)
D6 <- Die("D", 6)
```

## Random variables and their probability spaces

[`random_symbols()`](https://adamrsophie.github.io/symcompR/reference/random_symbols.md)
lists the random variables of an expression,
[`is_random()`](https://adamrsophie.github.io/symcompR/reference/random_symbols.md)
tests for them, and
[`pspace()`](https://adamrsophie.github.io/symcompR/reference/random_symbols.md)
returns the probability space.
[`given()`](https://adamrsophie.github.io/symcompR/reference/given.md)
conditions on an event, and
[`where()`](https://adamrsophie.github.io/symcompR/reference/given.md)
gives the domain on which a condition holds;
[`rv_attr()`](https://adamrsophie.github.io/symcompR/reference/rv_attr.md)
reads the attributes of these objects.

``` r

random_symbols(2 * X + Y)
#> [[1]]
#> X
#> 
#> [[2]]
#> Y
is_random(X^2)
#> [1] TRUE
pspace(X)
#> SingleContinuousPSpace(X, NormalDistribution(0, 1))
E(given(X, X > 0))
#> sqrt(2)/sqrt(pi)
rv_density(given(D6, D6 > 3))
#> {4: 1/3, 5: 1/3, 6: 1/3}
rv_attr(where(X^2 <= 1), "set")
#> Interval(-1, 1)
dependent(X, 2 * X)
#> [1] TRUE
independent(X, Y)
#> [1] TRUE
probability(X > 0)
#> 1/2
rs_swap(list(X, Y), list(Y, X))
#> [[1]]
#> [[1]][[1]]
#> X
#> 
#> [[1]][[2]]
#> X
#> 
#> 
#> [[2]]
#> [[2]][[1]]
#> Y
#> 
#> [[2]][[2]]
#> Y
```

[`Density()`](https://adamrsophie.github.io/symcompR/reference/E.md) is
the unevaluated density, and
[`characteristic_function()`](https://adamrsophie.github.io/symcompR/reference/E.md),
[`moment_generating_function()`](https://adamrsophie.github.io/symcompR/reference/E.md)
and
[`rv_quantile()`](https://adamrsophie.github.io/symcompR/reference/E.md)
are the other functions that describe a distribution.

``` r

Density(X)
#> Density(X)
doit(Density(X))
#> NormalDistribution(0, 1)
apply_lambda(characteristic_function(X), z)
#> exp(-z^2/2)
apply_lambda(moment_generating_function(X), z)
#> exp(z^2/2)
apply_lambda(characteristic_function(D6), z)
#> exp(6*I*z)/6 + exp(5*I*z)/6 + exp(4*I*z)/6 + exp(3*I*z)/6 + exp(2*I*z)/6 + exp(I*z)/6
apply_lambda(rv_quantile(Exponential("T", 2)), z)
#> -log(1 - z)/2
apply_lambda(rv_quantile(D6), S(1) / 2)
#> 3
rv_median(D6)
#> {3, 4}
rv_median(Exponential("T", 2))
#> {log(2)/2}
```

## Moments and other summaries

``` r

standard_deviation(Exponential("T", 2))
#> 1/2
smoment(X, 4)
#> 3
factorial_moment(Poisson("Po", 2), 2)
#> 4
coskewness(X, X, X)
#> 0
entropy(Normal("N", 0, 2))
#> 1/2 + log(pi)/2 + 3*log(2)/2
```

The capitalized functions are unevaluated, and evaluate with
[`doit()`](https://adamrsophie.github.io/symcompR/reference/doit.md):

``` r

Expectation(X^2)
#> Expectation(X^2)
doit(Expectation(X^2))
#> 1
doit(Variance(2 * X + 1))
#> Variance(2*X + 1)
expand(Variance(X + Y))
#> 2*Covariance(X, Y) + Variance(X) + Variance(Y)
Covariance(X, Y)
#> Covariance(X, Y)
doit(Covariance(X, 2 * X))
#> Covariance(X, 2*X)
Probability(X > 0)
#> Probability(X > 0)
doit(Probability(X > 0))
#> 1/2
doit(Moment(X, 3, 1))
#> -4
doit(CentralMoment(X, 4))
#> 3
rewrite(Variance(X), "Expectation")
#> -Expectation(X)^2 + Expectation(X^2)
```

[`RandomSymbol()`](https://adamrsophie.github.io/symcompR/reference/RandomSymbol.md)
wraps a symbol as a random variable without a distribution, and
[`variance_prop()`](https://adamrsophie.github.io/symcompR/reference/variance_prop.md)
propagates the variances of independent variables through an expression.

``` r

RandomSymbol(Symbol("x"))
#> x
variance_prop(Symbol("x") * Symbol("y"))
#> x^2*Variance(y) + y^2*Variance(x)
variance_prop(Symbol("x") + Symbol("y"), include_covar = TRUE)
#> 2*Covariance(x, y) + Variance(x) + Variance(y)
```

## Continuous distributions

Each distribution function takes a name and its parameters.

``` r

rv_density(Arcsin("X", 0, 1), z)
#> 1/(pi*sqrt(z*(1 - z)))
rv_density(Benini("X", 1, 2, 3), z)
#> 3*(4*log(z/3)/z + 1/z)*exp(-2*log(z/3)^2)/z
rv_density(BetaNoncentral("X", 1, 2, 3), z)
#> Sum((3/2)^_k*z^_k*(1 - z)*exp(-3/2)/(beta(_k + 1, 2)*factorial(_k)), (_k, 0, oo))
rv_density(BetaPrime("X", 2, 3), z)
#> 12*z/(z + 1)^5
rv_density(BoundedPareto("X", 2, 1, 4), z)
#> 32/(15*z^3)
rv_density(ChiRV("X", 3), z)
#> sqrt(2)*z^2*exp(-z^2/2)/sqrt(pi)
rv_density(ChiNoncentral("X", 2, 3), z)
#> z*exp(-z^2/2 - 9/2)*besseli(0, 3*z)
rv_density(Dagum("X", 1, 2, 3), z)
#> 2*z/(9*(z^2/9 + 1)^2)
rv_density(Davis("X", 2, 3, 1), z)
#> 4/((z - 1)^4*(exp(2/(z - 1)) - 1)*zeta(3))
rv_density(Erlang("X", 2, 3), z)
#> 9*z*exp(-3*z)
E(Erlang("X", 2, 3))
#> 2/3
rv_density(ExGaussian("X", 0, 1, 2), z)
#> exp(2 - 2*z)*erfc(sqrt(2)*(2 - z)/2)
rv_density(ExponentialPower("X", 0, 1, 2), z)
#> exp(-Abs(z)^2)/sqrt(pi)
rv_density(FDistribution("X", 2, 4), z)
#> 64*sqrt(z^2/(2*z + 4)^6)/z
rv_density(FisherZ("X", 2, 4), z)
#> 128*exp(2*z)/(2*exp(2*z) + 4)^3
rv_density(Frechet("X", 2, 1, 0), z)
#> 2*exp(-1/z^2)/z^3
rv_density(GammaInverse("X", 2, 3), z)
#> 9*exp(-3/z)/z^3
rv_density(GaussianInverse("X", 1, 2), z)
#> sqrt(z^(-3))*exp(-(z - 1)^2/z)/sqrt(pi)
rv_density(Wald("X", 1, 2), z)
#> sqrt(z^(-3))*exp(-(z - 1)^2/z)/sqrt(pi)
rv_density(Gompertz("X", 1, 2), z)
#> 2*exp(2)*exp(z)*exp(-2*exp(z))
rv_density(Gumbel("X", 1, 0), z)
#> exp(-z - exp(-z))
rv_density(Kumaraswamy("X", 2, 3), z)
#> 6*z*(1 - z^2)^2
rv_density(Levy("X", 0, 1), z)
#> sqrt(2)*exp(-1/(2*z))/(2*sqrt(pi)*z^(3/2))
rv_density(LogCauchy("X", 0, 1), z)
#> 1/(pi*z*(log(z)^2 + 1))
rv_density(LogLogistic("X", 1, 2), z)
#> 2*z/(z^2 + 1)^2
rv_density(Logistic("X", 0, 1), z)
#> exp(-z)/(1 + exp(-z))^2
rv_density(LogitNormal("X", 0, 1), z)
#> sqrt(2)*exp(-log(z/(1 - z))^2/2)/(2*sqrt(pi)*z*(1 - z))
rv_density(Lomax("X", 2, 3), z)
#> 2/(3*(z/3 + 1)^3)
rv_density(Maxwell("X", 2), z)
#> sqrt(2)*z^2*exp(-z^2/8)/(8*sqrt(pi))
E(Maxwell("X", 2))
#> 4*sqrt(2)/sqrt(pi)
rv_density(Moyal("X", 0, 1), z)
#> sqrt(2)*exp(-z/2 - exp(-z)/2)/(2*sqrt(pi))
rv_density(Nakagami("X", 1, 2), z)
#> z*exp(-z^2/2)
rv_density(PowerFunction("X", 2, 0, 1), z)
#> 2*z
rv_density(QuadraticU("X", 0, 2), z)
#> Piecewise((3*(z - 1)^2/2, (z >= 0) & (z <= 2)), (0, True))
rv_density(RaisedCosine("X", 0, 1), z)
#> Piecewise((cos(pi*z)/2 + 1/2, (z >= -1) & (z <= 1)), (0, True))
rv_density(Rayleigh("X", 2), z)
#> z*exp(-z^2/8)/4
E(Rayleigh("X", 2))
#> sqrt(2)*sqrt(pi)
rv_density(Reciprocal("X", 1, 4), z)
#> 1/(z*log(4))
rv_density(ShiftedGompertz("X", 1, 2), z)
#> (3 - 2*exp(-z))*exp(-z)*exp(-2*exp(-z))
rv_density(SkewNormal("X", 0, 1, 2), z)
#> sqrt(2)*(erf(sqrt(2)*z) + 1)*exp(-z^2/2)/(2*sqrt(pi))
rv_density(Trapezoidal("X", 0, 1, 2, 3), z)
#> Piecewise((z/2, (z >= 0) & (z < 1)), (1/2, (z >= 1) & (z < 2)), (3/2 - z/2, (z >= 2) & (z <= 3)), (0, True))
rv_density(Triangular("X", 0, 2, 1), z)
#> Piecewise((z, (z >= 0) & (z < 1)), (1, Eq(z, 1)), (2 - z, (z <= 2) & (z > 1)), (0, True))
E(Triangular("X", 0, 2, 1))
#> 1
rv_density(UniformSum("X", 2), z)
#> Sum((-1)^_k*(-_k + z)*binomial(2, _k), (_k, 0, floor(z)))
rv_density(VonMises("X", 0, 1), z)
#> exp(cos(z))/(2*pi*besseli(0, 1))
rv_density(WignerSemicircle("X", 1), z)
#> 2*sqrt(1 - z^2)/pi
E(WignerSemicircle("X", 1))
#> 0
```

[`ContinuousRV()`](https://adamrsophie.github.io/symcompR/reference/continuous_rv.md)
builds a random variable from any density:

``` r

s <- Symbol("x")
Xc <- ContinuousRV(s, exp(-s), set = Interval(0, oo))
rv_density(Xc, z)
#> Piecewise((exp(-z), (z >= 0) & (z < oo)), (0, True))
E(Xc)
#> 1
```

## Discrete and finite distributions

``` r

rv_density(FlorySchulz("X", S(1) / 2), z)
#> (1/2)^(z - 1)*z/4
rv_density(Hermite("X", 1, 2), z)
#> exp(-3)*Sum(2^_j/(factorial(_j)*factorial(-2*_j + z)), (_j, 0, floor(z/2)))
E(Hermite("X", 1, 2))
#> 5
rv_density(Logarithmic("X", S(1) / 2), z)
#> (1/2)^z/(z*log(2))
rv_density(NegativeBinomial("X", 5, S(1) / 3), z)
#> (2/3)^z*binomial(z + 4, z)/243
E(NegativeBinomial("X", 5, S(1) / 3))
#> 10
rv_density(Skellam("X", 1, 2), z)
#> (1/2)^(z/2)*exp(-3)*besseli(z, 2*sqrt(2))
E(Skellam("X", 1, 2))
#> -1
rv_density(YuleSimon("X", 3), z)
#> 3*beta(z, 4)
rv_density(Zeta("X", 3), z)
#> 1/(z^3*zeta(3))
```

[`DiscreteRV()`](https://adamrsophie.github.io/symcompR/reference/discrete_rv.md)
builds a discrete random variable from a probability mass function on a
set of integers:

``` r

Xd <- DiscreteRV(s, S(1) / 2^s, set = Naturals)
Xd
#> x
```

Finite distributions list their values and probabilities:

``` r

rv_density(BetaBinomial("X", 2, 1, 1))
#> BetaBinomialDistribution(2, 1, 1)
rv_density(Coin("C"))
#> BernoulliDistribution(1/2, H, T)
rv_density(Coin("C", S(1) / 3))
#> BernoulliDistribution(1/3, H, T)
rv_density(DiscreteUniform("X", list(1, 2, 3)))
#> DiscreteUniformDistribution(1, 2, 3)
Xf <- FiniteRV("X", list("1" = S(1) / 2, "2" = S(1) / 4, "3" = S(1) / 4))
rv_density(Xf)
#> FiniteDistributionHandmade({1: 1/2, 2: 1/4, 3: 1/4})
E(Xf)
#> 7/4
rv_density(Hypergeometric("X", 10, 5, 3))
#> HypergeometricDistribution(10, 5, 3)
rv_density(IdealSoliton("X", 3))
#> IdealSolitonDistribution(3)
rv_density(Rademacher("X"))
#> RademacherDistribution()
rv_density(RobustSoliton("X", 3, S(1) / 2, S(1) / 10))
#> RobustSolitonDistribution(3, 1/2, 1/10)
```

## Joint distributions

Joint random variables have densities in several variables;
[`rv_getitem()`](https://adamrsophie.github.io/symcompR/reference/rv_attr.md)
extracts a component (with SymPy’s zero based index).

``` r

M <- Multinomial("M", 3, S(1) / 2, S(1) / 2)
apply_lambda(rv_density(M), 1, 2)
#> 3/8
rv_getitem(M, 0)
#> M[0]
apply_lambda(rv_density(Dirichlet("D", list(1, 2))), x1, x2)
#> 2*x2
apply_lambda(rv_density(MultivariateBeta("B", 1, 2)), x1, x2)
#> 2*x2
apply_lambda(rv_density(MultivariateLaplace("L", list(0, 0), list(list(1, 0), list(0, 1)))), x1, x2)
#> besselk(0, sqrt(2*x1^2 + 2*x2^2))/pi
apply_lambda(rv_density(MultivariateT("T", list(0, 0), list(list(1, 0), list(0, 1)), 2)), x1, x2)
#> 1/(2*pi*(x1^2/2 + x2^2/2 + 1)^2)
apply_lambda(rv_density(NegativeMultinomial("N", 2, S(1) / 4, S(1) / 4)), x1, x2)
#> (1/4)^x1*(1/4)^x2*gamma(x1 + x2 + 2)/(4*factorial(x1)*factorial(x2))
apply_lambda(rv_density(MultivariateEwens("E", 2, 1)), 1, 0)
#> 0
apply_lambda(rv_density(NormalGamma("G", 0, 1, 2, 3)), x1, x2)
#> 9*sqrt(2)*x2^(3/2)*exp(-3*x2)*exp(-x1^2*x2/2)/(2*sqrt(pi))
MN2 <- MultivariateNormal("A", list(1, 2), list(list(1, 0), list(0, 1)))
apply_lambda(rv_density(MN2), 1, 2)
#> 1/(2*pi)
expand(apply_lambda(rv_density(MN2), x1, x2))
#> exp(-5/2)*exp(x1)*exp(-x1^2/2)*exp(2*x2)*exp(-x2^2/2)/(2*pi)
apply_lambda(rv_density(GeneralizedMultivariateLogGamma("G", S(1) / 2, 1, list(1, 1), list(1, 1))), x1, x2)
#> Sum((1/2)^n*exp((n + 1)*(x1 + x2) - exp(x1) - exp(x2))/gamma(n + 1)^2, (n, 0, oo))/2
omega <- Matrix(list(list(1, S(1) / 2), list(S(1) / 2, 1)))
apply_lambda(rv_density(GeneralizedMultivariateLogGammaOmega("G", omega, 1, list(1, 1), list(1, 1))), x1, x2)
#> 3*Sum((1/4)^n*exp((n + 1)*(x1 + x2) - exp(x1) - exp(x2))/gamma(n + 1)^2, (n, 0, oo))/4
```

[`JointRV()`](https://adamrsophie.github.io/symcompR/reference/joint_rv.md)
builds a joint distribution from a density written with indexed symbols,
and
[`marginal_distribution()`](https://adamrsophie.github.io/symcompR/reference/joint_rv.md)
integrates out components.

``` r

J <- JointRV("J", exp(-Indexed("J", 0) - Indexed("J", 1)),
             ProductSet(Interval(0, oo), Interval(0, oo)))
apply_lambda(rv_density(J), x1, x2)
#> exp(-x1 - x2)
apply_lambda(marginal_distribution(M, 0), 1)
#> 3/8
```

## Matrix distributions and random matrices

Matrix valued random variables take values in a set of matrices:

``` r

MG <- MatrixGamma("M", 1, 2, list(list(1, 0), list(0, 1)))
apply_lambda(rv_density(MG), MatrixSymbol("Xm", 2, 2))
#> exp(Trace(Matrix([
#> [-1/2,    0],
#> [   0, -1/2]])*Xm))/(4*pi*Determinant(Matrix([
#> [1, 0],
#> [0, 1]]))*sqrt(Determinant(Xm)))
rv_attr(rv_attr(rv_attr(MG, "pspace"), "distribution"), "set")
#> MatrixSet(2, 2, Reals)
W <- Wishart("W", 3, list(list(1, 0), list(0, 1)))
apply_lambda(rv_density(W), MatrixSymbol("Xm", 2, 2))
#> exp(Trace(Matrix([
#> [-1/2,    0],
#> [   0, -1/2]])*Xm))/(4*pi*Determinant(Matrix([
#> [1, 0],
#> [0, 1]]))^(3/2))
MN <- MatrixNormal("M", list(list(1, 2)), list(1), list(list(1, 0), list(0, 1)))
rv_attr(rv_attr(rv_attr(MN, "pspace"), "distribution"), "set")
#> MatrixSet(1, 2, Reals)
MT <- MatrixStudentT("M", 2, list(list(1, 2)), list(list(1, 0), list(0, 1)), list(1))
rv_attr(rv_attr(rv_attr(MT, "pspace"), "distribution"), "set")
#> MatrixSet(1, 2, Reals)
```

The Gaussian and circular ensembles of random matrix theory, with their
joint eigenvalue and level spacing distributions:

``` r

GaussianEnsemble("U", 2)
#> U
CircularEnsemble("C", 2)
#> C
GaussianUnitaryEnsemble("U", 2)
#> U
apply_lambda(level_spacing_distribution(GaussianUnitaryEnsemble("U", 2)), z)
#> 32*z^2*exp(-4*z^2/pi)/pi^2
apply_lambda(level_spacing_distribution(GaussianOrthogonalEnsemble("U", 2)), z)
#> pi*z*exp(-pi*z^2/4)/2
apply_lambda(level_spacing_distribution(GaussianSymplecticEnsemble("U", 2)), z)
#> 262144*z^4*exp(-64*z^2/(9*pi))/(729*pi^3)
joint_eigen_distribution(GaussianUnitaryEnsemble("U", 2))
#> Lambda((l[1], l[2]), exp(-l[1]^2 - l[2]^2)*Product(Abs(l[_i] - l[_j])^2, (_j, _i + 1, 2), (_i, 1, 1))/pi)
joint_eigen_distribution(CircularUnitaryEnsemble("C", 2))
#> Lambda((t[1], t[2]), Product(Abs(exp(I*t[_j]) - exp(I*t[_k]))^2, (_j, _k + 1, 2), (_k, 1, 1))/(8*pi^2))
joint_eigen_distribution(CircularOrthogonalEnsemble("C", 2))
#> Lambda((t[1], t[2]), Product(Abs(exp(I*t[_j]) - exp(I*t[_k])), (_j, _k + 1, 2), (_k, 1, 1))/(16*pi))
joint_eigen_distribution(CircularSymplecticEnsemble("C", 2))
#> Lambda((t[1], t[2]), Product(Abs(exp(I*t[_j]) - exp(I*t[_k]))^4, (_j, _k + 1, 2), (_k, 1, 1))/(24*pi^2))
JointEigenDistribution(Matrix(list(list(Normal("a", 0, 1), 0), list(0, Normal("b", 0, 1)))))
#> JointDistributionHandmade(a, b)
RandomMatrixSymbol("R", 2, 2)
#> R
```

## Stochastic processes

``` r

B <- BernoulliProcess("B", p = S(1) / 2)
E(rv_getitem(B, 1) + rv_getitem(B, 2))
#> 1
P(rv_getitem(BernoulliProcess("B", p = S(1) / 3), 1) > 0)
#> 1/3
Xp <- PoissonProcess("X", 3)
E(Xp(2))
#> 6
P(Eq(Xp(2), 1))
#> 6*exp(-6)
E(WienerProcess("W")(3))
#> 0
variance(WienerProcess("W")(3))
#> 3
E(GammaProcess("G", 2, 3)(4))
#> 6
variance(GammaProcess("G", 2, 3)(4))
#> 3
RandomIndexedSymbol(Indexed("x", 1))
#> x[1]
```
