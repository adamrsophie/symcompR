# NumPy and SciPy codegen functions

Functions after SymPy's `sympy.codegen.numpy_nodes` and `scipy_nodes`:
`logaddexp(a, b)` is `log(exp(a) + exp(b))`, `logaddexp2(a, b)` is
`log2(2^a + 2^b)`, `amin()` and `amax()` are array reductions,
`np_maximum()` and `np_minimum()` are the elementwise maximum and
minimum (SymPy's `maximum` and `minimum`), `cosm1(x)` is `cos(x) - 1`
and `powm1(x, y)` is `x^y - 1`. The Python printers print them as the
NumPy and SciPy calls;
[`rewrite()`](https://adamrsophie.github.io/symcompR/reference/rewrite.md)
turns them into elementary functions.

## Usage

``` r
logaddexp(a, b)

logaddexp2(a, b)

amin(...)

amax(...)

np_maximum(...)

np_minimum(...)

cosm1(x)

powm1(x, y)
```

## Arguments

- a, b, x, y:

  Expressions.

- ...:

  Arguments.

## Value

An expression or node.

## Examples

``` r
x <- Symbol("x"); y <- Symbol("y")
logaddexp(x, y)
#> logaddexp(x, y)
diff(logaddexp(x, y), x)
#> 1/(exp(-x + y) + 1)
rewrite(cosm1(x), "cos")
#> cos(x) - 1
```
