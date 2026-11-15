# Ground type integer functions

The integer functions of `sympy.polys.domains.groundtypes`, which SymPy
takes from Python or from gmpy2: `python_gcd()`, `python_lcm()` and
`python_gcdex()` (Python's), and `gmpy_gcd()`, `gmpy_lcm()`,
`gmpy_gcdex()`, `gmpy_numer()`, `gmpy_denom()` and `gmpy_qdiv()`
(gmpy2's, here computed with gmp). `python_gcdex(a, b)` returns
`list(x, y, g)` with `x*a + y*b = g`, while `gmpy_gcdex(a, b)` returns
gmpy2's order `list(g, x, y)`. `gmpy_qdiv(x, y)` is the exact quotient,
an integer when it is one. `sqrtrem(n)` returns the integer square root
`s` and the remainder `n - s^2`.

## Usage

``` r
python_gcd(a, b)

python_lcm(a, b)

python_gcdex(a, b)

gmpy_gcd(a, b)

gmpy_lcm(a, b)

gmpy_gcdex(a, b)

gmpy_numer(x)

gmpy_denom(x)

gmpy_qdiv(x, y = 1)

sqrtrem(n)
```

## Arguments

- a, b, x, y, n:

  Integers (rationals for `gmpy_numer()`, `gmpy_denom()` and
  `gmpy_qdiv()`).

## Value

An integer, a rational, or a list of integers.

## Examples

``` r
python_gcdex(12, 18)
#> [[1]]
#> [1] -1
#> 
#> [[2]]
#> [1] 1
#> 
#> [[3]]
#> [1] 6
#> 
gmpy_gcdex(12, 18)
#> [[1]]
#> [1] 6
#> 
#> [[2]]
#> [1] -1
#> 
#> [[3]]
#> [1] 1
#> 
sqrtrem(10)
#> [[1]]
#> [1] 3
#> 
#> [[2]]
#> [1] 1
#> 
gmpy_qdiv(6, 4)
#> 3/2
```
