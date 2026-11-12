# Diophantine equations

`diophantine(eq)` finds the integer solutions of a polynomial equation
with integer coefficients, as a `FiniteSet` of `Tuple`s ordered like the
sorted free symbols (or `syms`). Parametric families use integer
parameters such as `t_0`, `t_1`. Linear, binary quadratic (including
Pell equations), homogeneous ternary quadratic, general Pythagorean,
sums of squares and sums of even powers are supported, and products of
such factors. `permute = TRUE` adds the sign and order permutations
SymPy adds. `classify_diop(eq)` returns the variables, coefficients and
type; `diop_solve(eq)` solves one equation of a known type.

## Usage

``` r
diophantine(
  eq,
  param = Symbol("t", integer = TRUE),
  syms = NULL,
  permute = FALSE
)

classify_diop(eq)

diop_solve(eq, param = Symbol("t", integer = TRUE))
```

## Arguments

- eq:

  An equation or an expression equal to zero.

- param:

  The base name of the parameters.

- syms:

  Optional order of the variables.

- permute:

  Add sign and order permutations of the solutions.

## Value

A set of tuples; `classify_diop()` returns a list.

## Examples

``` r
x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
diophantine(2*x + 3*y - 5)
#> {(3*t_0 - 5, 5 - 2*t_0)}
diophantine(x^2 - 2*y^2 - 1)
#> {(-sqrt(2)*(3 - 2*sqrt(2))^t + 3*(3 - 2*sqrt(2))^t/2 + sqrt(2)*(2*sqrt(2) + 3)^t + 3*(2*sqrt(2) + 3)^t/2, -(3 - 2*sqrt(2))^t + 3*sqrt(2)*(3 - 2*sqrt(2))^t/4 - 3*sqrt(2)*(2*sqrt(2) + 3)^t/4 - (2*sqrt(2) + 3)^t), (-3*(3 - 2*sqrt(2))^t/2 + sqrt(2)*(3 - 2*sqrt(2))^t - 3*(2*sqrt(2) + 3)^t/2 - sqrt(2)*(2*sqrt(2) + 3)^t, -3*sqrt(2)*(3 - 2*sqrt(2))^t/4 + (3 - 2*sqrt(2))^t + (2*sqrt(2) + 3)^t + 3*sqrt(2)*(2*sqrt(2) + 3)^t/4)}
diophantine(x^2 + y^2 - z^2)
#> {(2*p*q, p^2 - q^2, p^2 + q^2)}
classify_diop(x^2 + y^2 - x*y + x + 5)
#> [[1]]
#> [[1]][[1]]
#> x
#> 
#> [[1]][[2]]
#> y
#> 
#> 
#> [[2]]
#> [[2]]$`x^2`
#> 1
#> 
#> [[2]]$`x*y`
#> -1
#> 
#> [[2]]$x
#> 1
#> 
#> [[2]]$`y^2`
#> 1
#> 
#> [[2]]$`1`
#> 5
#> 
#> 
#> [[3]]
#> [1] "binary_quadratic"
#> 
```
