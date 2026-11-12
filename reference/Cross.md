# Unevaluated dot and cross products

`Dot(v1, v2)` and `Cross(v1, v2)` represent products without computing
them;
[`doit()`](https://adamrsophie.github.io/symcompR/reference/doit.md)
evaluates. As in SymPy, `Dot` orders its arguments and `Cross(v2, v1)`
becomes `-Cross(v1, v2)` when `v1` sorts first.

## Usage

``` r
Cross(expr1, expr2)

Dot(expr1, expr2)
```

## Arguments

- expr1, expr2:

  Vectors.

## Value

A `sym` expression.

## Examples

``` r
C <- CoordSys3D("C")
Dot(C$i + C$j, C$k)
#> Dot(C.k, C.i + C.j)
doit(Cross(C$i, C$j))
#> C.k
```
