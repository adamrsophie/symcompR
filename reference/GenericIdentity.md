# Generic identity and zero matrices

`GenericIdentity()` and `GenericZeroMatrix()` are the matrix identity
and zero of unspecified shape that
[`MatMul()`](https://adamrsophie.github.io/symcompR/reference/matrix-expressions.md)
and
[`MatAdd()`](https://adamrsophie.github.io/symcompR/reference/matrix-expressions.md)
of no arguments give. They print as `I` and `0`; asking for their shape
is an error.

## Usage

``` r
GenericIdentity()

GenericZeroMatrix()
```

## Value

A `sym` matrix expression.

## Examples

``` r
GenericIdentity()
#> I
MatMul()
#> I
MatAdd()
#> 0
```
