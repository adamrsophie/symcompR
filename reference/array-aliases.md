# Array aliases and indexed objects

`NDimArray()`, `MutableDenseNDimArray()`, `ImmutableDenseNDimArray()`
and `DenseNDimArray()` are
[`Array()`](https://adamrsophie.github.io/symcompR/reference/Array.md),
and `shape()` its dimensions (the shape tuple of a matrix expression).
Indexed objects are documented in
[`Indexed()`](https://adamrsophie.github.io/symcompR/reference/Indexed.md).

## Usage

``` r
NDimArray(...)

MutableDenseNDimArray(...)

ImmutableDenseNDimArray(...)

DenseNDimArray(...)

SparseNDimArray(...)

MutableSparseNDimArray(...)

ImmutableSparseNDimArray(...)

shape(a)
```

## Arguments

- ...:

  Arguments for
  [`Array()`](https://adamrsophie.github.io/symcompR/reference/Array.md).

- a:

  An array or a matrix expression.

## Value

An array, function or `sym` expression.

## Examples

``` r
shape(Array(list(1, 2, 3)))
#> [1] 3
```
