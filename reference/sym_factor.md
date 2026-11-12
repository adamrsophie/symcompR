# Factor polynomials

`sym_factor()` factors a polynomial or rational function into
irreducibles over the rationals, in any number of variables, like
SymPy's `sym_factor()`. `factor_list()` returns `list(coeff, factors)`
where `factors` is a list of `list(factor, multiplicity)`. `sqf()`,
`sqf_list()` and `sqf_part()` do the same for square-free decomposition.

## Usage

``` r
factor_sym(f, ..., deep = FALSE)

factor_list(f, ...)

sqf_list(f, ...)

sqf(f, ...)

sqf_part(f, ...)

sym_factor(x, ...)
```

## Arguments

- f:

  An expression.

- ...:

  Generators.

- deep:

  Factor inside function arguments.

- x:

  An expression.

## Value

An expression or a list.

## Examples

``` r
x <- Symbol("x"); y <- Symbol("y")
sym_factor(x^4 - 1)
#> (x - 1)*(x + 1)*(x^2 + 1)
sym_factor(x^2 * y - y^3)
#> y*(x - y)*(x + y)
factor_list(2 * x^3 - 2 * x)
#> $coeff
#> 2
#> 
#> $factors
#> $factors[[1]]
#> $factors[[1]][[1]]
#> x - 1
#> 
#> $factors[[1]][[2]]
#> [1] 1
#> 
#> 
#> $factors[[2]]
#> $factors[[2]][[1]]
#> x
#> 
#> $factors[[2]][[2]]
#> [1] 1
#> 
#> 
#> $factors[[3]]
#> $factors[[3]][[1]]
#> x + 1
#> 
#> $factors[[3]][[2]]
#> [1] 1
#> 
#> 
#> 
sqf((x - 1)^2 * (x + 2))
#> (x - 1)^2*(x + 2)
```
