# Powers of an algebraic integer, and a coefficient search

`AlgIntPowers(T, modulus)` is SymPy's class of that name: given the
monic irreducible integer polynomial `T` of an algebraic integer `theta`
of degree `n`, `p[[e]]` (or `alg_int_power(p, e)`) gives the
coefficients of `theta^e` over `1, theta, ..., theta^(n-1)`, lowest
first, from the linear recurrence that `T` defines, reduced modulo
`modulus` when one is given.

## Usage

``` r
AlgIntPowers(T, modulus = NULL)

alg_int_power(p, e)

# S3 method for class 'AlgIntPowers'
x[[i]]

coeff_search(m, R)
```

## Arguments

- T:

  A monic `Poly` over the integers.

- modulus:

  An optional integer modulus.

- p:

  An `AlgIntPowers` object.

- e:

  A non-negative exponent.

- x, i:

  An `AlgIntPowers` object and an exponent.

- m:

  The length of the coefficient lists.

- R:

  The initial bound.

## Value

`AlgIntPowers()` an object; the power a big integer vector;
`coeff_search()` a function of no arguments returning integer vectors.

## Details

`coeff_search(m, R)` is SymPy's generator of coefficient lists of length
`m`, all with entries bounded by `R` in absolute value before the bound
grows, the leading entry never negative, skipping repeats. R has no
generators, so it returns a function; each call gives the next list.

## Examples

``` r
x <- Symbol("x")
zp <- AlgIntPowers(Poly(cyclotomic_poly(5, x)))
zp[[4]]
#> Big Integer ('bigz') object of length 4:
#> [1] -1 -1 -1 -1
zp[[24]]
#> Big Integer ('bigz') object of length 4:
#> [1] -1 -1 -1 -1
nxt <- coeff_search(2, 1)
t(sapply(1:13, function(i) nxt()))
#>       [,1] [,2]
#>  [1,]    1    1
#>  [2,]    1    0
#>  [3,]    1   -1
#>  [4,]    0    1
#>  [5,]    2    2
#>  [6,]    2    1
#>  [7,]    2    0
#>  [8,]    2   -1
#>  [9,]    2   -2
#> [10,]    1    2
#> [11,]    1   -2
#> [12,]    0    2
#> [13,]    3    3
```
