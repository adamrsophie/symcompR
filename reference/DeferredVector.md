# Deferred vectors for lambdify

SymPy's `DeferredVector(name)`: a vector whose entries `X[i]` (from 0)
are symbols named `"X[i]"`; given to
[`lambdify()`](https://adamrsophie.github.io/symcompR/reference/lambdify.md)
as an argument, the function takes one vector `X` and reads its entries.

## Usage

``` r
DeferredVector(name)

# S3 method for class 'DeferredVector'
x[i]
```

## Arguments

- name:

  The vector name.

- x:

  A deferred vector.

- i:

  An index, from 0.

## Value

A deferred vector; indexing gives a symbol.

## Examples

``` r
X <- DeferredVector("X")
f <- lambdify(X, X[0] + 2 * X[1])
f(c(3, 4))
#> [1] 11
```
