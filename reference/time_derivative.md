# Time derivative in a reference frame

The derivative with respect to time of a vector, dyadic or scalar as
seen from `frame`, as SymPy's `time_derivative()`. The time symbol is
`t`.

## Usage

``` r
time_derivative(expr, frame, order = 1)
```

## Arguments

- expr:

  A `PhysVector`, `PhysDyadic` or scalar expression.

- frame:

  A
  [`ReferenceFrame()`](https://adamrsophie.github.io/symcompR/reference/ReferenceFrame.md).

- order:

  Order of the derivative.

## Value

An object of the same kind as `expr`.

## Examples

``` r
N <- ReferenceFrame("N")
q <- dynamicsymbols("q")
A <- N$orientnew("A", "Axis", list(q, N$x))
time_derivative(A$y, N)
#> Derivative(q(t), t)*A.z
```
