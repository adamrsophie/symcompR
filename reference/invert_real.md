# Invert real and complex valued functions

`invert_real(f_x, y, x)` and `invert_complex(f_x, y, x)` reduce the
equation `f_x = y` to `g(x) = h` for the simplest `g(x)` they can reach,
and return the pair `(g(x), set of h)`, as SymPy's functions of the same
names. When `g(x)` is `x` the set holds the solutions.

## Usage

``` r
invert_real(f_x, y, x)

invert_complex(f_x, y, x, domain = Complexes)
```

## Arguments

- f_x:

  An expression in `x`.

- y:

  The right hand side, free of `x`.

- x:

  The symbol to invert for.

- domain:

  The domain, `Complexes` by default for `invert_complex()`.

## Value

A tuple `(g(x), set)`.

## Examples

``` r
x <- Symbol("x"); y <- Symbol("y")
invert_real(exp(x), y, x)
#> (x, Intersection({log(y)}, Reals)) 
invert_complex(exp(x), y, x)
#> (x, ImageSet(Lambda(_n, I*(2*_n*pi + arg(y)) + log(Abs(y))), Integers)) 
```
