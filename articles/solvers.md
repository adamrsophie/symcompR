# Solving equations

symcompR solves algebraic equations and systems, inequalities,
recurrence relations, Diophantine equations, and ordinary and partial
differential equations. As in SymPy, there are two interfaces for
algebraic equations: [`solve()`](https://rdrr.io/r/base/solve.html)
returns a list of solutions, and
[`solveset()`](https://adamrsophie.github.io/symcompR/reference/solveset.md)
returns a mathematical set, which can also describe infinite families of
solutions and solutions it cannot find explicitly.

``` r

library(symcompR)
x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c")
n <- Symbol("n", integer = TRUE)
f <- Function("f")
```

An expression passed to a solver is understood as an equation with right
hand side zero; use
[`Eq()`](https://adamrsophie.github.io/symcompR/reference/Eq.md) for
other equations.

## `solve()`

[`solve()`](https://rdrr.io/r/base/solve.html) is an S3 method of base
R’s [`solve()`](https://rdrr.io/r/base/solve.html) generic for symbolic
input. It returns a list of solutions.

``` r

solve(x^2 - 4, x)
#> [[1]]
#> -2
#> 
#> [[2]]
#> 2
solve(Eq(2 * x + 1, 7), x)
#> [[1]]
#> 3
solve(x^3 - 6 * x^2 + 11 * x - 6, x)
#> [[1]]
#> 1
#> 
#> [[2]]
#> 2
#> 
#> [[3]]
#> 3
solve(x^4 - 1, x)
#> [[1]]
#> -1
#> 
#> [[2]]
#> 1
#> 
#> [[3]]
#> -I
#> 
#> [[4]]
#> I
solve(exp(x) - 3, x)
#> [[1]]
#> log(3)
solve(sqrt(x) - 2, x)
#> [[1]]
#> 4
```

`dict = TRUE` labels each solution with its variable. When the roots of
a polynomial cannot be written with radicals, they are returned as
indexed `CRootOf` objects, which can be evaluated to any precision.

``` r

solve(x^2 - 4, x, dict = TRUE)
#> [[1]]
#> [[1]]$x
#> -2
#> 
#> 
#> [[2]]
#> [[2]]$x
#> 2
r <- solve(x^5 - x + 1, x)
r[[1]]
#> CRootOf(x^5 - x + 1, 0)
N(r[[1]])
#> -1.16730397826142
```

For a system, give a list of equations and a list of unknowns. Each
solution is a named list.

``` r

solve(list(x + y - 3, x - y - 1), list(x, y))
#> $x
#> 2
#> 
#> $y
#> 1
```

## `solveset()`

`solveset(f, x, domain)` returns the solution set. The domain defaults
to the complex numbers; pass `Reals` for real solutions.

``` r

solveset(x^2 - 1, x)
#> {-1, 1}
solveset(x^3 - 1, x)
#> {1, -1/2 - sqrt(3)*I/2, -1/2 + sqrt(3)*I/2}
solveset(x^2 + 1, x)
#> {-I, I}
solveset(x^2 + 1, x, Reals)
#> EmptySet
solveset((x - 1)^2 * (x + 2), x)
#> {-2, 1}
solveset(x^2 - 2, x, Interval(0, 5))
#> {sqrt(2)}
```

The answer is always a set, so the edge cases have a precise meaning:

``` r

solveset(x - x, x)   # every complex number is a solution
#> Complexes
solveset(1 / x, x)   # no solution
#> EmptySet
```

### Periodic solutions

Trigonometric and exponential equations have infinitely many solutions.
[`solveset()`](https://adamrsophie.github.io/symcompR/reference/solveset.md)
returns them as an `ImageSet`, the image of the integers under a
`Lambda` function:

``` r

solveset(sin(x) - 1, x, Reals)
#> ImageSet(Lambda(_n, 2*_n*pi + pi/2), Integers)
solveset(sin(x) - Rational(1, 2), x, Reals)
#> Union(ImageSet(Lambda(_n, 2*_n*pi + pi/6), Integers), ImageSet(Lambda(_n, 2*_n*pi + 5*pi/6), Integers))
solveset(exp(x) - 1, x)
#> ImageSet(Lambda(_n, 2*_n*I*pi), Integers)
solveset(exp(x) - 1, x, Reals)
#> {0}
```

The first result reads \\2n\pi + \pi/2 : n \in \mathbb{Z}\\.

### `ConditionSet`

When
[`solveset()`](https://adamrsophie.github.io/symcompR/reference/solveset.md)
cannot solve an equation, it does not fail: it returns a `ConditionSet`,
the set of points of the domain where the equation holds.

``` r

solveset(x * exp(x) - 1, x, Reals)
#> ConditionSet(x, Eq(x*exp(x) - 1, 0), Reals)
solveset(sin(x) + x, x, Reals)
#> ConditionSet(x, Eq(x + sin(x), 0), Reals)
```

### Inequalities

[`solveset()`](https://adamrsophie.github.io/symcompR/reference/solveset.md)
also accepts inequalities, and returns intervals.

``` r

solveset(Gt(x, 3), x, Reals)
#> Interval.open(3, oo)
solveset(Lt(x^2, 4), x, Reals)
#> Interval.open(-2, 2)
solveset(abs(x) - 2, x, Reals)
#> {-2, 2}
```

## Linear and polynomial systems

[`linsolve()`](https://adamrsophie.github.io/symcompR/reference/linsolve.md)
solves systems of linear equations. As in SymPy, it returns a
`FiniteSet` holding one tuple of values; when the system is
underdetermined, free unknowns stay as themselves and the others are
written in terms of them.

``` r

linsolve(list(x + 2 * y - 3, 3 * x - y - 2), list(x, y))
#> {(1, 1)}
linsolve(list(x + y + z - 1, x + y + 2 * z - 3), list(x, y, z))
#> {(-y - 1, y, 2)}
```

Here z = 2, x = -y - 1, and y is free. Nonlinear polynomial systems are
solved by
[`nonlinsolve()`](https://adamrsophie.github.io/symcompR/reference/nonlinsolve.md),
which returns a set of solution tuples, and by
[`solve_poly_system()`](https://adamrsophie.github.io/symcompR/reference/solver-helpers.md)
through Groebner bases.

``` r

nonlinsolve(list(x * y - 1, x - 2), list(x, y))
#> {(2, 1/2)}
nonlinsolve(list(x^2 + y^2 - 4, x - y), list(x, y))
#> {(-sqrt(2), -sqrt(2)), (sqrt(2), sqrt(2))}
```

Several lower level solvers are available on their own.
[`solve_linear()`](https://adamrsophie.github.io/symcompR/reference/solver-helpers.md)
isolates one symbol in a linear expression,
[`solve_linear_system()`](https://adamrsophie.github.io/symcompR/reference/solver-helpers.md)
and
[`solve_linear_system_LU()`](https://adamrsophie.github.io/symcompR/reference/more-solvers.md)
take an augmented matrix, and
[`solve_undetermined_coeffs()`](https://adamrsophie.github.io/symcompR/reference/solver-helpers.md)
matches the coefficients of a polynomial identity.
[`checksol()`](https://adamrsophie.github.io/symcompR/reference/solver-helpers.md)
verifies a candidate solution, and
[`substitution()`](https://adamrsophie.github.io/symcompR/reference/substitution.md)
solves a system by successive substitution and returns a set of tuples.

``` r

show <- function(l) vapply(l, sstr, "")
rows <- function(l) vapply(l, function(r) paste(show(r), collapse = ", "), "")
show(solve_linear(x - y^2, 0, x))
#> [1] "x"   "y^2"
show(solve_linear(3 * x - 2 * y))
#> [1] "x"     "2*y/3"
show(solve_linear_system(Matrix(list(list(1, 1, 3), list(1, -1, 1))), list(x, y)))
#>   x   y 
#> "2" "1"
show(solve_linear_system_LU(Matrix(list(list(1, 2, 0, 1), list(3, 2, 2, 1),
                                        list(2, 0, 0, 1))), list(x, y, z)))
#>      x      y      z 
#>  "1/2"  "1/4" "-1/2"
show(solve_undetermined_coeffs(Eq(2 * a * x + a + b, x), list(a, b), x))
#>      a      b 
#>  "1/2" "-1/2"
checksol(x^4 - 1, list(x = 1))
#> [1] TRUE
substitution(list(x + y), list(x))
#> {(-y,)}
```

[`solve_triangulated()`](https://adamrsophie.github.io/symcompR/reference/more-solvers.md)
solves a polynomial system through a triangular Groebner basis,
[`factor_system()`](https://adamrsophie.github.io/symcompR/reference/more-solvers.md)
splits a system into simpler systems by factoring,
[`decompogen()`](https://adamrsophie.github.io/symcompR/reference/more-solvers.md)
decomposes a function into a chain of compositions, and
[`homogeneous_order()`](https://adamrsophie.github.io/symcompR/reference/more-solvers.md)
gives the order of a homogeneous function (or `NULL` if it is not
homogeneous).

``` r

rows(solve_triangulated(list(x^2 + y + z - 1, x + y^2 + z - 1,
                             x + y + z^2 - 1), list(x, y, z)))
#> [1] "0, 0, 1" "0, 1, 0" "1, 0, 0"
rows(factor_system(list(x^2 - 1, (x - 1) * y)))
#> [1] "x + 1, y" "x - 1"
show(decompogen(sin(cos(x)), x))
#> [1] "sin(x)" "cos(x)"
show(decompogen(sqrt(6 * x^2 - 5), x))
#> [1] "sqrt(x)"   "6*x^2 - 5"
homogeneous_order(x^2 + x * y, x, y)
#> 2
homogeneous_order(x^2 + 1, x)
#> NULL
```

## Inequalities

Besides
[`solveset()`](https://adamrsophie.github.io/symcompR/reference/solveset.md),
several functions reduce inequalities to their simplest form.
[`reduce_inequalities()`](https://adamrsophie.github.io/symcompR/reference/inequality-solvers.md)
handles systems and returns a logical combination of relations;
[`solve_univariate_inequality()`](https://adamrsophie.github.io/symcompR/reference/solve_univariate_inequality.md)
can return either relations or a set.

``` r

reduce_inequalities(Gt(x^2 - 4, 0), x)
#> ((-oo < x) & (x < -2)) | ((2 < x) & (x < oo))
reduce_inequalities(list(Gt(x, 0), Lt(x^2, 9)), x)
#> (0 < x) & (x < 3)
solve_univariate_inequality(Ge(x^2 - 1, 0), x)
#> ((1 <= x) & (x < oo)) | ((-oo < x) & (x <= -1))
solve_univariate_inequality(Ge(x^2 - 1, 0), x, relational = FALSE)
#> Union(Interval(-oo, -1), Interval(1, oo))
solve_poly_inequality(x^2 - 3 * x + 2, "<")
#> [[1]]
#> Interval.open(1, 2)
reduce_rational_inequalities(list(list(Gt((x - 1) / (x + 2), 0))), x)
#> ((-oo < x) & (x < -2)) | ((1 < x) & (x < oo))
reduce_abs_inequality(abs(x - 2) - 3, "<", x)
#> (-1 < x) & (x < 5)
reduce_abs_inequalities(list(list(Abs(3 * x - 5) - 7, "<"),
                             list(Abs(x + 25) - 13, ">")), x)
#> (-2/3 < x) & (x < 4) & (((-oo < x) & (x < -38)) | ((-12 < x) & (x < oo)))
solve_rational_inequalities(list(list(list(list(Poly(-x + 1, x), Poly(1, x)), ">="))))
#> Interval(-oo, 1)
```

## Numerical solutions

[`nsolve()`](https://adamrsophie.github.io/symcompR/reference/nsolve.md)
finds a root numerically with Newton’s method, starting from an initial
guess. For systems, give lists of equations, unknowns, and starting
values.

``` r

nsolve(cos(x) - x, x, 1)
#> 0.739085133215161
nsolve(x^3 - 2 * x - 5, x, 2)
#> 2.09455148154233
nsolve(list(x^2 + y^2 - 4, x * y - 1), list(x, y), list(2, 0.5))
#> $x
#> 1.93185165257814
#> 
#> $y
#> 0.517638090205041
```

## Recurrence relations

`rsolve(eq, y(n), n, init)` solves linear recurrences with constant or
polynomial coefficients. Initial conditions are pairs `list(n0, value)`.

``` r

rsolve(f(n + 1) - 2 * f(n), f(n), n, init = list(list(0, 3)))
#> 3*2^n
rsolve(f(n + 2) - 5 * f(n + 1) + 6 * f(n), f(n), n)
#> 2^n*C0 + 3^n*C1
rsolve(f(n + 1) - (n + 1) * f(n), f(n), n, init = list(list(0, 1)))
#> factorial(n)
rsolve(f(n + 1) - f(n) - 2^n, f(n), n, init = list(list(0, 1)))
#> 2^n
```

Without initial conditions the result contains arbitrary constants `C0`,
`C1`, and so on.

[`rsolve_poly()`](https://adamrsophie.github.io/symcompR/reference/rsolve-parts.md),
[`rsolve_ratio()`](https://adamrsophie.github.io/symcompR/reference/rsolve-parts.md),
and
[`rsolve_hyper()`](https://adamrsophie.github.io/symcompR/reference/rsolve-parts.md)
find the polynomial, rational, and hypergeometric solutions of a
recurrence given by its list of coefficients and its right hand side.

``` r

rsolve_poly(list(-1, 1), 4 * n^3, n)
#> C0 + n^4 - 2*n^3 + n^2
rsolve_ratio(list(-2 * n^3 + n^2 + 2 * n - 1, 2 * n^3 + n^2 - 6 * n,
                  -2 * n^3 - 11 * n^2 - 18 * n - 9,
                  2 * n^3 + 13 * n^2 + 22 * n + 8), 0, n)
#> C0*(2*n - 3)/(2*(n^2 - 1))
rsolve_hyper(list(-1, -1, 1), 0, n)
#> C0*(1/2 - sqrt(5)/2)^n + C1*(1/2 + sqrt(5)/2)^n
rsolve_hyper(list(-2, 1), 2^n, n)
#> 2^n*C0 + 2^(n - 1)*n
```

## Diophantine equations

[`diophantine()`](https://adamrsophie.github.io/symcompR/reference/diophantine.md)
finds all integer solutions of a polynomial equation. Linear equations
give a parametric family, Pythagorean type equations a parametrization,
and equations with finitely many solutions an explicit set.

``` r

diophantine(3 * x + 5 * y - 1)
#> {(5*t_0 + 2, -3*t_0 - 1)}
diophantine(x^2 - y^2 - 5)
#> {(-3, -2), (-3, 2), (3, -2), (3, 2)}
diophantine(x^2 + y^2 - 25)
#> {(-5, 0), (-4, -3), (-4, 3), (-3, -4), (-3, 4), (0, -5), (0, 5), (3, -4), (3, 4), (4, -3), (4, 3), (5, 0)}
diophantine(x^2 + y^2 - z^2)
#> {(2*p*q, p^2 - q^2, p^2 + q^2)}
```

The parametric solution of 3x + 5y = 1 is x = 5t + 2, y = -3t - 1 for
any integer t.
[`classify_diop()`](https://adamrsophie.github.io/symcompR/reference/diophantine.md)
reports the type of equation that was recognized.

``` r

classify_diop(x^2 + y^2 - z^2)[[3]]
#> [1] "homogeneous_ternary_quadratic_normal"
diop_solve(2 * x + 3 * y - 5)
#> (3*t_0 - 5, 5 - 2*t_0)
diop_solve(x^2 + y^2 - 5)
#> {(-2, -1), (-2, 1), (-1, -2), (-1, 2), (1, -2), (1, 2), (2, -1), (2, 1)}
```

## Ordinary differential equations

`dsolve(eq, f(x))` solves an ODE for the unknown function. Arbitrary
constants are named `C1`, `C2`, and so on.

``` r

dsolve(Derivative(f(x), x) - f(x), f(x))
#> Eq(f(x), C1*exp(x))
dsolve(diff(f(x), x) - 2 * f(x), f(x))
#> Eq(f(x), C1*exp(2*x))
dsolve(diff(f(x), x, 2) + f(x), f(x))
#> Eq(f(x), C1*sin(x) + C2*cos(x))
dsolve(diff(f(x), x, 2) + 4 * f(x), f(x))
#> Eq(f(x), C1*sin(2*x) + C2*cos(2*x))
dsolve(diff(f(x), x, 2) - 4 * f(x), f(x))
#> Eq(f(x), C1*exp(-2*x) + C2*exp(2*x))
dsolve(diff(f(x), x, 4) - f(x), f(x))
#> Eq(f(x), C1*exp(-x) + C2*exp(x) + C3*sin(x) + C4*cos(x))
```

First order separable, linear, Bernoulli, and exact equations:

``` r

dsolve(diff(f(x), x) + 2 * x * f(x), f(x))
#> Eq(f(x), C1*exp(-x^2))
dsolve(diff(f(x), x) - x * exp(-f(x)), f(x))
#> Eq(f(x), log(C1 + x^2/2))
dsolve(diff(f(x), x) + f(x) - x, f(x))
#> Eq(f(x), C1*exp(-x) + x - 1)
dsolve(diff(f(x), x) - f(x)^2, f(x))
#> Eq(f(x), -1/(C1 + x))
dsolve(diff(f(x), x) - f(x) / x - x, f(x))
#> Eq(f(x), x*(C1 + x))
dsolve(x * diff(f(x), x) + f(x) - x^3, f(x))
#> Eq(f(x), (C1 + x^4/4)/x)
dsolve((2 * x * f(x) + 1) + x^2 * diff(f(x), x), f(x))
#> Eq(f(x), (C1 - x)/x^2)
```

Nonhomogeneous linear equations with constant coefficients are solved
with undetermined coefficients:

``` r

dsolve(diff(f(x), x, 2) + f(x) - x, f(x))
#> Eq(f(x), C1*sin(x) + C2*cos(x) + x)
```

[`classify_ode()`](https://adamrsophie.github.io/symcompR/reference/dsolve.md)
names the solution methods (SymPy’s hints) that apply to an equation. A
homogeneous linear equation with constant coefficients has one:

``` r

classify_ode(diff(f(x), x, 3) - f(x), f(x))
#> [1] "nth_linear_constant_coeff_homogeneous"
```

### Initial conditions

`ics` is a list of conditions `list(x0, k, value)`, meaning f^{(k)}(x_0)
= \text{value}.

``` r

dsolve(Eq(diff(f(x), x), f(x)), f(x), ics = list(list(0, 0, 1)))
#> Eq(f(x), exp(x))
dsolve(diff(f(x), x) + f(x) - x, f(x), ics = list(list(0, 0, 2)))
#> Eq(f(x), x - 1 + 3*exp(-x))
dsolve(diff(f(x), x, 2) + f(x), f(x),
       ics = list(list(0, 0, 0), list(0, 1, 1)))
#> Eq(f(x), sin(x))
dsolve(diff(f(x), x, 2) + 4 * f(x), f(x),
       ics = list(list(0, 0, 1), list(0, 1, 0)))
#> Eq(f(x), cos(2*x))
```

[`ode_order()`](https://adamrsophie.github.io/symcompR/reference/dsolve.md)
gives the order of an equation, and
[`checkodesol()`](https://adamrsophie.github.io/symcompR/reference/solver-helpers.md)
verifies a solution by substitution.

``` r

ode_order(diff(f(x), x, 3) + f(x), f(x))
#> [1] 3
checkodesol(diff(f(x), x) - f(x), Eq(f(x), exp(x)), f(x))
#> (True, 0)
```

## Partial differential equations

[`pdsolve()`](https://adamrsophie.github.io/symcompR/reference/pdsolve.md)
solves first order linear PDEs with constant coefficients and some with
variable coefficients. The solution contains an arbitrary function `F`
of the characteristic variable.

``` r

pdsolve(2 * diff(f(x, y), x) + 3 * diff(f(x, y), y))
#> Eq(f(x, y), F(3*x - 2*y))
pdsolve(diff(f(x, y), x) + diff(f(x, y), y) + f(x, y))
#> Eq(f(x, y), F(x - y)*exp(-x/2 - y/2))
pdsolve(diff(f(x, y), x) - 3 * diff(f(x, y), y) + 2 * f(x, y) - x)
#> Eq(f(x, y), x/2 + F(-3*x - y)*exp(-x/5 + 3*y/5) - 1/4)
classify_pde(x * diff(f(x, y), x) + diff(f(x, y), y))
#> [1] "1st_linear_variable_coeff"
```

The first solution says that f is constant along the lines 3x - 2y =
\text{const}, the characteristics of the equation.

[`pde_separate()`](https://adamrsophie.github.io/symcompR/reference/pdsolve.md),
[`pde_separate_add()`](https://adamrsophie.github.io/symcompR/reference/pdsolve.md),
and
[`pde_separate_mul()`](https://adamrsophie.github.io/symcompR/reference/pdsolve.md)
separate the variables of a partial differential equation for an
additive or multiplicative ansatz, and
[`checkpdesol()`](https://adamrsophie.github.io/symcompR/reference/pdsolve.md)
checks a solution.

``` r

u <- Function("u"); X <- Function("X"); Tf <- Function("T"); Y <- Function("Y")
tt <- Symbol("t")
show(pde_separate(Eq(diff(u(x, tt), x), diff(u(x, tt), tt)), u(x, tt),
                  list(X(x), Tf(tt))))
#> [1] "Derivative(X(x), x)/X(x)" "Derivative(T(t), t)/T(t)"
show(pde_separate_add(Eq(diff(u(x, tt), x), exp(diff(u(x, tt), tt))),
                      u(x, tt), list(X(x), Tf(tt))))
#> [1] "Derivative(X(x), x)"      "exp(Derivative(T(t), t))"
show(pde_separate_mul(Eq(diff(u(x, y), x, 2) + diff(u(x, y), y, 2), 0),
                      u(x, y), list(X(x), Y(y))))
#> [1] "Derivative(X(x), (x, 2))/X(x)"  "-Derivative(Y(y), (y, 2))/Y(y)"
Fn <- Function("F")
checkpdesol(diff(f(x, y), x) + diff(f(x, y), y) + f(x, y),
            Eq(f(x, y), Fn(x - y) * exp(-x / 2 - y / 2)))[[1]]
#> [1] TRUE
```
