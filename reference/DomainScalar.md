# Elements of a domain with their domain

`DomainScalar(element, domain)` is SymPy's
`sympy.polys.matrices.domainscalar.DomainScalar`, as returned by
indexing a
[`DomainMatrix()`](https://adamrsophie.github.io/symcompR/reference/DomainMatrix.md).
Arithmetic unifies the domains; `$to_sympy()`, `$to_domain(K)`,
`$convert_to(K)`, `$unify(other)`, `$is_zero()` and `$is_one()` are its
methods, and `DomainScalar_from_sympy(expr)` builds one from an
expression.

## Usage

``` r
DomainScalar(element, domain)

DomainScalar_from_sympy(expr)

# S3 method for class 'DomainScalar'
x$name

# S3 method for class 'DomainScalar'
format(x, ...)

# S3 method for class 'DomainScalar'
print(x, ...)

# S3 method for class 'DomainScalar'
Ops(e1, e2)

# S3 method for class 'DomainScalar'
chooseOpsMethod(x, y, mx, my, cl, reverse)
```

## Arguments

- element:

  An element of `domain`.

- domain:

  A domain.

- expr:

  An expression.

- x:

  A `DomainScalar`.

- name:

  A method name.

- ...:

  Unused.

- e1, e2:

  Scalars or integers.

- y, mx, my, cl, reverse:

  Arguments of
  [`base::chooseOpsMethod()`](https://rdrr.io/r/base/chooseOpsMethod.html).

## Value

A `DomainScalar`, or the value of the method.

## Examples

``` r
a <- DomainScalar(ZZ(3), ZZ)
b <- DomainScalar(QQ(1, 2), QQ)
a + b
#> MPQ(7,2)
a$to_sympy()
#> 3
```
