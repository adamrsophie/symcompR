# Rewrite an expression in terms of other functions

`rewrite(e, "exp")` expresses trigonometric and hyperbolic functions
with exponentials; `"cos"`, `"sin"` and `"tan"` convert between
trigonometric functions; `"log"` writes inverse trigonometric functions
as logarithms; `"gamma"` writes factorials and binomials with gamma
functions; `"factorial"` does the reverse; `"Piecewise"` rewrites `Abs`,
`sign`, `Heaviside`, `Max` and `Min`; `"erf"` rewrites `erfc` and
`erfi`. The `"erfc"`, `"sign"`, `"hyper"`, `"SingularityFunction"`,
`"asin"` and `"tractable"` targets follow the corresponding SymPy
methods. A target that no function knows returns the expression
unchanged, as in SymPy.

## Usage

``` r
rewrite(e, target, ...)
```

## Arguments

- e:

  An expression.

- target:

  The target family as a string, or a function such as `sign` standing
  for its class.

- ...:

  Options for the rewrite methods of quantum objects.

## Value

An expression.

## Examples

``` r
x <- Symbol("x")
rewrite(sin(x), "exp")
#> -I*(exp(I*x) - exp(-I*x))/2
rewrite(sym_factorial(Symbol("n")), "gamma")
#> gamma(n + 1)
```
