# Unevaluated differential operators

`Gradient(f)`, `Divergence(v)`, `Curl(v)` and `Laplacian(f)` represent
the operators without computing them;
[`doit()`](https://adamrsophie.github.io/symcompR/reference/doit.md)
evaluates.

## Usage

``` r
Gradient(expr)

Divergence(expr)

Curl(expr)

Laplacian(expr)
```

## Arguments

- expr:

  A scalar or vector field.

## Value

A `sym` expression.

## Examples

``` r
R <- CoordSys3D("R")
g <- Gradient(R$x * R$y * R$z)
g
#> Gradient(R.x*R.y*R.z)
doit(g)
#> R.y*R.z*R.i + R.x*R.z*R.j + R.x*R.y*R.k
```
