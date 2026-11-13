# Power series rings over the integers and the rationals

`PythonPowerSeriesRingZZ(prec)` and `PythonPowerSeriesRingQQ(prec)` are
univariate power series rings with coefficients in ZZ or QQ, as SymPy's
classes of the same names. A series is the pair `(coeffs, prec)`: the
coefficients from the highest degree down, and the order of the error
term `O(x**prec)` or `NULL` for an exact polynomial. Call the ring on a
vector of coefficients (lowest degree first) to make a series, and use
its methods with `$`: `add`, `subtract`, `multiply`, `divide`,
`pow_int`, `square`, `compose`, `inverse`, `reversion`, `truncate`,
`differentiate`, `pretty`, `to_list`, and over QQ `integrate`, `sqrt`,
`log`, `exp`, `sin`, `cos`, `tan`, `atan`, `hypot` and the hyperbolic
functions.

## Usage

``` r
PythonPowerSeriesRingZZ(prec = 6L)

PythonPowerSeriesRingQQ(prec = 6L)

# S3 method for class 'PythonPowerSeriesRing'
x$name

# S3 method for class 'PythonPowerSeriesRing'
format(x, ...)

# S3 method for class 'PythonPowerSeriesRing'
print(x, ...)
```

## Arguments

- prec:

  The ring precision, a nonnegative integer (default 6).

- x:

  A power series ring.

- name:

  The method name.

- ...:

  Unused.

## Value

A power series ring. Called on coefficients it returns a series.

## Examples

``` r
R <- PythonPowerSeriesRingZZ()
s <- R(c(1, 2, 3))
R$pretty(R$pow_int(s, 2))
#> [1] "1 + 4*x + 10*x**2 + 12*x**3 + 9*x**4"
R$pretty(R$inverse(R(c(1, 1))))
#> [1] "1 - x + x**2 - x**3 + x**4 - x**5 + O(x**6)"
Q <- PythonPowerSeriesRingQQ(4)
Q$pretty(Q$exp(Q$gen))
#> [1] "1 + x + 1/2*x**2 + 1/6*x**3 + O(x**4)"
```
