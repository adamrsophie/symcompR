# Solvable quintics

`PolyQuintic(f)` holds the special functions of SymPy's
`sympy.polys.polyquinticconst.PolyQuintic` for a quintic
`x^5 + p*x^3 + q*x^2 + r*x + s` with rational coefficients: `$f20` (the
resolvent sextic as a polynomial in `x`), `$F`, `$zeta` (the primitive
fifth roots of unity), the coefficient lists `$b`, `$o`, `$a` and `$c`
of SymPy's properties of the same names, and the functions `$l0(theta)`,
`$T(theta, d)`, `$order(theta, d)` and `$uv(theta, d)`.
`roots_quintic(f)` returns the five roots of a solvable irreducible
quintic in radicals, or an empty list when the quintic is reducible or
not solvable; `roots(f, quintics = TRUE)` uses it.

## Usage

``` r
PolyQuintic(f)

# S3 method for class 'PolyQuintic'
x$name

# S3 method for class 'PolyQuintic'
print(x, ...)

roots_quintic(f)
```

## Arguments

- f:

  A quintic polynomial expression (or `Poly`) in one variable.

- x:

  A `PolyQuintic` object.

- name:

  The name of a field or method.

- ...:

  Unused.

## Value

`PolyQuintic()` returns an object of class `PolyQuintic`;
`roots_quintic()` a list of expressions.

## Examples

``` r
x <- Symbol("x")
q <- PolyQuintic(x^5 - 5 * x + 12)
q$f20
#> Poly(x^6 - 40*x^5 + 1000*x^4 - 20000*x^3 + 250000*x^2 - 66400000*x + 976000000, x, domain='ZZ')
# \donttest{
# solving the quintic in radicals takes about 30 seconds
length(roots_quintic(x^5 - 5 * x + 12))
#> [1] 5
# }
```
