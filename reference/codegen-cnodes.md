# Codegen nodes for C, C++, Python and lists

Nodes after SymPy's `sympy.codegen.cnodes`, `cxxnodes`, `pynodes` and
`abstract_nodes`. `alignof()` and `sizeof()` are calls to the C
operators; `CommaOperator()` is the C comma operator; `Label()` and
`goto()` are a labeled statement and a jump; `PreIncrement()`,
`PostIncrement()`, `PreDecrement()` and `PostDecrement()` are the C
increment and decrement operators; `struct()` and `c_union()` declare a
struct and a union (SymPy's `union`); `void` is the void type and
`restrict`, `volatile` and `static` are attributes. `using()` is a C++
using declaration and `constexpr` an attribute. `List()` is a list,
which prints with square brackets, and `NumExprEvaluate()` a call to
`numexpr.evaluate`.

## Usage

``` r
List(...)

NumExprEvaluate(...)

void

restrict

volatile

static

alignof(arg)

sizeof(arg)

CommaOperator(...)

Label(...)

goto(...)

PreDecrement(...)

PostDecrement(...)

PreIncrement(...)

PostIncrement(...)

struct(...)

c_union(...)

using(...)

constexpr
```

## Arguments

- ...:

  Arguments or fields of the node.

- arg:

  A type name or an expression.

## Value

A node.

## Examples

``` r
x <- Symbol("x")
ccode(sizeof("unsigned int"))
#> [1] "sizeof(unsigned int)"
ccode(CommaOperator(PreIncrement(x), 2 * x))
#> [1] "(++(x), 2*x)"
List(1, 2, 3)
#> [1, 2, 3]
```
