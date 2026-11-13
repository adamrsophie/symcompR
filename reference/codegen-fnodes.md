# Fortran codegen nodes

Nodes after SymPy's `sympy.codegen.fnodes`: `Program()`, `Module()`,
`Subroutine()` and `SubroutineCall()` build program units; `Do()` is a
do loop, `ArrayConstructor()` an array constructor and `ImpliedDoLoop()`
an implied do loop; `GoTo()` is a computed goto and `FortranReturn()` a
return statement; `fortran_use()` and `use_rename()` are use statements
(SymPy's `use`); `Extent()` is an array extent, `fortran_dimension()`
the dimension attribute (SymPy's `dimension`) and `fortran_array()`
declares an array (SymPy's `array`); `allocated()`, `lbound()`,
`ubound()`, `fortran_shape()`, `size()` and `fortran_reshape()` call the
intrinsic functions (SymPy's `shape`); `bind_C()` is the bind(C)
attribute; `isign()`, `dsign()`, `cmplx()`, `kind()` and
`fortran_merge()` (SymPy's `merge`) are intrinsic functions;
`literal_sp()` and `literal_dp()` are single and double precision
literals; `sum_()` and `product_()` are the array reductions.

## Usage

``` r
pure

elemental

intent_in

intent_out

intent_inout

allocatable

Program(...)

use_rename(...)

fortran_use(...)

Module(...)

Subroutine(...)

SubroutineCall(...)

Do(...)

ArrayConstructor(...)

ImpliedDoLoop(...)

Extent(...)

assumed_extent

fortran_dimension(...)

assumed_size

fortran_array(
  symbol,
  dim,
  intent = NULL,
  attrs = list(),
  value = NULL,
  type = NULL
)

allocated(array)

lbound(array, dim = NULL, kind = NULL)

ubound(array, dim = NULL, kind = NULL)

fortran_shape(source, kind = NULL)

size(array, dim = NULL, kind = NULL)

fortran_reshape(source, shape, pad = NULL, order = NULL)

bind_C(name = NULL)

GoTo(...)

FortranReturn(...)

isign(x, y)

dsign(x, y)

cmplx(x, y)

kind(x)

fortran_merge(x, y, z)

literal_sp(x)

literal_dp(x)

sum_(...)

product_(...)
```

## Arguments

- ...:

  Arguments or fields of the node.

- symbol, source, array:

  A symbol, an array or an expression.

- dim, kind, shape, pad, order, intent, attrs, value, type, name:

  Optional arguments of the Fortran functions and attributes.

- x, y, z:

  Expressions.

## Value

A node or expression.

## Examples

``` r
x <- Symbol("x"); i <- Symbol("i", integer = TRUE); n <- Symbol("n", integer = TRUE)
fcode(Do(list(AddAugmentedAssignment(x, i)), i, 1, n), source_format = "free")
#> [1] "do i = 1, n\n   x = x + i\nend do"
fcode(isign(1, x), source_format = "free")
#> [1] "isign(1, x)"
```
