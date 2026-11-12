# Simplify signs

`signsimp()` takes the sign out of sums inside products, as SymPy's
`signsimp()` (`-x*(-x + y) + y` becomes `x*(x - y) + y`), and `posify()`
replaces symbols with positive ones.

## Usage

``` r
signsimp(e)

posify(e)
```

## Arguments

- e:

  An expression.

## Value

An expression; `posify()` returns `list(expr, reps)`.
