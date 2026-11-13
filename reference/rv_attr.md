# Attributes and items of stats objects

`rv_attr(x, name, ...)` reads a SymPy attribute of a random symbol,
probability space, domain or distribution (such as `"pspace"`,
`"domain"`, `"set"`, `"symbol"`, `"distribution"`, `"as_boolean"`), or
calls a method with further arguments. `rv_getitem(x, i)` indexes a
finite density, a joint random symbol or a stochastic process.

## Usage

``` r
rv_attr(x, name, ...)

rv_getitem(x, i, ...)
```

## Arguments

- x:

  A stats object.

- name:

  The attribute or method name.

- ...:

  Method arguments.

- i:

  An index or key.

## Value

The attribute value.

## Examples

``` r
X <- Normal("X", 0, 1)
rv_attr(rv_attr(X, "pspace"), "domain")
#> SingleContinuousDomain(X, Interval(-oo, oo))
rv_attr(where(X^2 <= 1), "set")
#> Interval(-1, 1)
```
