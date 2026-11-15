# Kinds of expressions

`expr_kind(e)` is SymPy's `e.kind`: the mathematical kind of an object,
independent of its class. Numbers, commutative symbols and arithmetic of
them are `NumberKind`; Booleans and relations `BooleanKind`; matrices
`MatrixKind(element kind)`; arrays `ArrayKind(...)`; sets
`SetKind(...)`; tuples `TupleKind(...)`; vectors `VectorKind(...)`;
quantum states and operators `KetKind`, `BraKind` and `OperatorKind`;
anything else, including applied functions other than `exp`,
`UndefinedKind`.
([`kind()`](https://adamrsophie.github.io/symcompR/reference/codegen-fnodes.md)
is the Fortran intrinsic of the code generation nodes.)

## Usage

``` r
expr_kind(e)

NumberKind

UndefinedKind

BooleanKind

KetKind

BraKind

OperatorKind

MatrixKind(element_kind = NumberKind)

ArrayKind(element_kind = NumberKind)

VectorKind(element_kind = NumberKind)

SetKind(element_kind = NULL)

TupleKind(...)
```

## Format

The kinds without parameters are objects.

## Arguments

- e:

  An expression, matrix, array, set or vector.

- element_kind:

  The kind of the elements.

- ...:

  Kinds of the tuple entries.

## Value

A kind.

## Examples

``` r
x <- Symbol("x")
expr_kind(x + 1)
#> NumberKind 
expr_kind(sin(x))
#> UndefinedKind 
expr_kind(x > 1)
#> BooleanKind 
expr_kind(Matrix(list(1, x)))
#> MatrixKind(NumberKind) 
expr_kind(FiniteSet(1, x))
#> SetKind(NumberKind) 
MatrixKind() == MatrixKind(NumberKind)
#> [1] TRUE
```
