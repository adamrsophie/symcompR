# Solve an equation as a list, as solve() would

`solvify(f, symbol, domain)` solves `f = 0` with
[`solveset()`](https://adamrsophie.github.io/symcompR/reference/solveset.md)
and returns the solutions as a list; for periodic equations it lists the
solutions in one period, as SymPy's `solvify()`. It returns `NULL` when
the solutions cannot be listed.

## Usage

``` r
solvify(f, symbol, domain)
```

## Arguments

- f:

  An expression.

- symbol:

  The variable.

- domain:

  The domain, a set.

## Value

A list of expressions, or `NULL`.

## Examples

``` r
x <- Symbol("x")
solvify(x^2 - 9, x, Reals)
#> [[1]]
#> -3
#> 
#> [[2]]
#> 3
#> 
solvify(sin(x) + 1, x, Reals)
#> [[1]]
#> 3*pi/2
#> 
```
