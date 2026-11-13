# Attributes and methods of quantum objects

Quantum objects keep the Python attributes and methods of their SymPy
classes. `qattr(obj, name)` is SymPy's `obj.name` (a property is
evaluated), `qcall(obj, name, ...)` calls the method `obj.name(...)`,
`qisinstance(obj, cls)` is `isinstance(obj, cls)` with the class given
by its name (or a vector of names), and `qclass(obj)` is the class name.
For other expressions `qattr()` knows `args`, `is_commutative` and
`free_symbols`.

## Usage

``` r
qattr(obj, name)

qcall(obj, name, ...)

qisinstance(obj, cls)

qclass(obj)
```

## Arguments

- obj:

  A quantum object, a class name or an expression.

- name:

  An attribute or method name.

- ...:

  Arguments of the method.

- cls:

  A class name or a character vector of names.

## Value

The attribute value, the method result, or a logical.

## Examples

``` r
A <- Operator("A")
qattr(A, "label")
#> [[1]]
#> A
#> 
qisinstance(A, "QExpr")
#> [1] TRUE
qcall(HermitianOperator("H"), "inverse")
#> H^(-1)
```
