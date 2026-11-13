# Classify an ordinary differential equation

`classify_ode(eq, func)` returns the hints, in SymPy's order, of the
methods
[`dsolve()`](https://adamrsophie.github.io/symcompR/reference/dsolve.md)
can use to solve the ODE `eq`, as SymPy's `classify_ode()`. With
`dict = TRUE` it returns a dictionary of the matching hints with the
keys `order`, `default` and `ordered_hints`.

## Usage

``` r
classify_ode(
  eq,
  func = NULL,
  dict = FALSE,
  ics = NULL,
  prep = TRUE,
  xi = NULL,
  eta = NULL,
  n = NULL,
  ...
)
```

## Arguments

- eq:

  The equation.

- func:

  The unknown function; detected when `NULL`.

- dict:

  Whether to return the dictionary of matches.

- ics:

  Initial conditions, as for
  [`dsolve()`](https://adamrsophie.github.io/symcompR/reference/dsolve.md).

- prep:

  Whether to preprocess the equation.

- xi, eta:

  Infinitesimals for the `lie_group` hint.

- n:

  The number of power series terms.

- ...:

  Further options: `hint` and `x0`.

## Value

A character vector of hints, or a dictionary.

## Examples

``` r
x <- Symbol("x"); f <- Function("f")
classify_ode(Derivative(f(x), x) - f(x), f(x))
#>  [1] "separable"                            
#>  [2] "1st_exact"                            
#>  [3] "1st_linear"                           
#>  [4] "Bernoulli"                            
#>  [5] "almost_linear"                        
#>  [6] "1st_power_series"                     
#>  [7] "lie_group"                            
#>  [8] "nth_linear_constant_coeff_homogeneous"
#>  [9] "separable_Integral"                   
#> [10] "1st_exact_Integral"                   
#> [11] "1st_linear_Integral"                  
#> [12] "Bernoulli_Integral"                   
#> [13] "almost_linear_Integral"               
```
