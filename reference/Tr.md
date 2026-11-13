# Trace of a quantum expression

`Tr(expr, indices)` is the trace of `sympy.physics.quantum.trace`:
scalars come out of it, a sum gives a sum of traces, a matrix gives its
trace, and
[`doit()`](https://adamrsophie.github.io/symcompR/reference/doit.md)
evaluates traces of density operators, outer products and tensor
products (a partial trace over `indices`).

## Usage

``` r
Tr(...)
```

## Arguments

- ...:

  The expression and the optional indices to trace over.

## Value

An expression.

## Examples

``` r
A <- Symbol("A", commutative = FALSE)
a <- Symbol("a")
Tr(a * A)
#> a*Tr(A)
Tr(Matrix(list(list(1, 1), list(2, 2))))
#> 3
```
