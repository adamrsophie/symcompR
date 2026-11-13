# Holonomic functions

``` r

library(symcompR)
x <- Symbol("x"); n <- Symbol("n")
```

This vignette covers `sympy.holonomic`: holonomic functions, the
differential and recurrence operators that define them, and conversions
to and from hypergeometric and Meijer G functions.

## Holonomic functions and operators

A holonomic function is a solution of a linear ordinary differential
equation with polynomial coefficients, together with initial conditions.
Sums and products of holonomic functions are holonomic, so many
manipulations reduce to operations on differential operators.

[`DifferentialOperators()`](https://adamrsophie.github.io/symcompR/reference/holonomic.md)
creates the algebra of operators in D_x. Below, p is e^x (the solution
of y' = y, y(0) = 1) and q is \sin x (the solution of y'' = -y, y(0) =
0, y'(0) = 1).

``` r

R <- DifferentialOperators(holonomic_ring(x), "Dx")
Dx <- R$Dx
p <- HolonomicFunction(Dx - 1, x, 0, list(1))
q <- HolonomicFunction(Dx^2 + 1, x, 0, list(0, 1))
to_expr_holonomic(p)
#> exp(x)
to_expr_holonomic(q * q)
#> sin(x)^2
```

[`expr_to_holonomic()`](https://adamrsophie.github.io/symcompR/reference/holonomic.md)
converts an expression, after which series, sequences of coefficients,
derivatives, and numerical values are available.

``` r

h <- expr_to_holonomic(exp(x) * sin(x), x)
series_holonomic(h, n = 6)
#> x + x^2 + x^3/3 - x^5/30 + O(x^6)
to_expr_holonomic(diff_holonomic(expr_to_holonomic(sin(x), x)))
#> cos(x)
to_expr_holonomic(expr_to_holonomic(log(1 + x), x))
#> log(x + 1)
to_hyper(expr_to_holonomic(exp(x), x))
#> hyper((), (), x)
```

[`evalf_holonomic()`](https://adamrsophie.github.io/symcompR/reference/holonomic.md)
integrates the differential equation numerically along a path of points,
with the fourth order Runge-Kutta method by default, and returns the
values at those points. Integrating y' = y from 0 to 1 in steps of 0.1
approximates e \approx 2.718282:

``` r

path <- as.list(seq(0.1, 1, by = 0.1))
vals <- evalf_holonomic(p, path)
vapply(vals, sstr, "")
#>  [1] "1.10517083333333" "1.22140257085069" "1.34985849706254" "1.49182424008069"
#>  [5] "1.64872063859684" "1.82211796209193" "2.01375162659678" "2.22553956329232"
#>  [9] "2.45960141378007" "2.71827974413517"
sstr(evalf_holonomic(expr_to_holonomic(sin(x), x), path)[[10]])
#> [1] "0.841470477800274"
```

The last value approximates \sin(1) \approx 0.841471.

## Transformations

[`change_x()`](https://adamrsophie.github.io/symcompR/reference/holonomic.md)
substitutes a polynomial for the variable,
[`shift_x()`](https://adamrsophie.github.io/symcompR/reference/holonomic.md)
shifts it, and
[`change_ics()`](https://adamrsophie.github.io/symcompR/reference/holonomic.md)
moves the point where the initial conditions are given.
[`holonomic_degree()`](https://adamrsophie.github.io/symcompR/reference/holonomic.md)
is the largest degree of a polynomial coefficient of the annihilator,
and
[`holonomic_indicial()`](https://adamrsophie.github.io/symcompR/reference/holonomic.md)
the indicial equation at the origin.

``` r

shift_x(p, 2)
#> HolonomicFunction((-1) + (1)*Dx, x, -2, [1])
change_ics(expr_to_holonomic(sin(x), x), 1)
#> HolonomicFunction((1) + (1)*Dx^2, x, 1, [sin(1), cos(1)])
change_x(p, 2 * x)
#> HolonomicFunction((-1) + (1)*Dx, 2*x, 0, [1])
holonomic_degree(expr_to_holonomic(exp(x^2), x))
#> [1] 1
holonomic_indicial(expr_to_holonomic(sqrt(x) * exp(x), x))
#> [[1]]
#> 1/2
#> 
#> attr(,"multiplicity")
#> [1] 1
```

[`composition()`](https://adamrsophie.github.io/symcompR/reference/holonomic.md)
substitutes an algebraic expression into a holonomic function. Without
initial conditions only the annihilator is computed; with them, the
result can be converted back to an expression.

``` r

composition(HolonomicFunction(Dx - 1, x), x^2 + x)
#> HolonomicFunction((-2*x - 1) + (1)*Dx, x)
to_expr_holonomic(composition(p, x^2, 0, list(1)))
#> exp(x^2)
```

[`integrate_holonomic()`](https://adamrsophie.github.io/symcompR/reference/holonomic.md)
integrates, either as a new holonomic function or over limits:

``` r

integrate_holonomic(q, list(x, 0, 1))
#> 1 - cos(1)
to_expr_holonomic(integrate_holonomic(p, list(x, 0, x)))
#> exp(x) - 1
```

## Recurrences

The coefficients of the power series of a holonomic function satisfy a
linear recurrence.
[`RecurrenceOperators()`](https://adamrsophie.github.io/symcompR/reference/holonomic.md)
builds the algebra of shift operators,
[`HolonomicSequence()`](https://adamrsophie.github.io/symcompR/reference/holonomic.md)
a sequence defined by a recurrence, and
[`to_sequence()`](https://adamrsophie.github.io/symcompR/reference/holonomic.md)
the recurrence for a holonomic function.

``` r

RS <- RecurrenceOperators(holonomic_ring(n), "Sn")
Sn <- RS$Sn
HolonomicSequence(Sn - n - 1, list(1))
#> HolonomicSequence((-n - 1) + (1)Sn, n), u(0) = 1
to_sequence(expr_to_holonomic(x^3, x))
#> [[1]]
#> [[1]][[1]]
#> HolonomicSequence((n), n), u(0) = 1 
#> 
#> [[1]][[2]]
#> 3
#> 
#> [[1]][[3]]
#> [1] 1
to_sequence(p)
#> [[1]]
#> [[1]][[1]]
#> HolonomicSequence((-1) + (n + 1)Sn, n), u(0) = 1 
#> 
#> [[1]][[2]]
#> [1] 0
```

## Hypergeometric and Meijer G functions

[`from_hyper()`](https://adamrsophie.github.io/symcompR/reference/holonomic.md)
and
[`from_meijerg()`](https://adamrsophie.github.io/symcompR/reference/holonomic.md)
turn special functions into holonomic functions, and
[`to_meijerg()`](https://adamrsophie.github.io/symcompR/reference/holonomic.md)
goes back.

``` r

from_hyper(hyper(list(1, 1), list(Rational(3, 2)), x^2 / 4))
#> HolonomicFunction((4*x) + (5*x^2 - 8)*Dx + (x^3 - 4*x)*Dx^2, x, 1, [2*sqrt(3)*pi/9, -4*sqrt(3)*pi/27 + 4/3])
from_meijerg(meijerg(list(), list(), list(0), list(), x))
#> HolonomicFunction((1) + (1)*Dx, x, 0, [1])
to_meijerg(expr_to_holonomic(sin(x), x))
#> sqrt(pi)*x*meijerg(((), ()), ((0,), (-1/2,)), x^2/4)/2
```
