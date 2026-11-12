# Symbolic sums and products

`summation(f, list(i, a, b))` computes `sum(f, i = a..b)` in closed form
with SymPy's `Sum.doit()`: direct sums of short ranges, telescoping,
Faulhaber polynomials, geometric sums (a Piecewise for ratio 1),
Gosper's algorithm, hypergeometric series through
[`hyperexpand()`](https://adamrsophie.github.io/symcompR/reference/hyperexpand.md)
(a Piecewise with the convergence condition), sums of rational functions
over infinite ranges by residues, `zeta` sums and Kronecker deltas
(`deltasummation()`). `gosper_sum(f, k)` is the indefinite Gosper sum.
`product(f, list(i, a, b))` is SymPy's `Product.doit()`: direct products
of short ranges, rising factorials over the roots of polynomial factors,
sums in exponents, and Kronecker deltas (`deltaproduct()`). `Sum()` and
`Product()` build unevaluated objects; use
[`doit()`](https://adamrsophie.github.io/symcompR/reference/doit.md).
[`is_convergent()`](https://adamrsophie.github.io/symcompR/reference/sum-limits.md)
applies standard convergence tests.

## Usage

``` r
deltaproduct(f, limit)

deltasummation(f, limit, no_piecewise = FALSE)

summation(f, ...)

Sum(f, ...)

product(f, ...)

Product(f, ...)

gosper_sum(f, ...)
```

## Arguments

- f:

  The summand or factor.

- limit:

  A limit `list(i, a, b)`.

- no_piecewise:

  Do not wrap the result in a Piecewise.

- ...:

  Limits `list(i, a, b)`.

## Value

An expression.

## Examples

``` r
i <- Symbol("i", integer = TRUE); n <- Symbol("n", integer = TRUE)
summation(i^2, list(i, 1, n))
#> n^3/3 + n^2/2 + n/6
summation(1 / i^2, list(i, 1, oo))
#> pi^2/6
product(i, list(i, 1, n))
#> factorial(n)
```
