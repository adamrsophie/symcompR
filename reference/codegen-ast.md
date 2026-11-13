# Code generation abstract syntax tree

Nodes of an abstract syntax tree for code generation, after SymPy's
`sympy.codegen.ast`. Each constructor takes its fields positionally or
by name, as the SymPy classes do: for example
`Variable(x, type = real)`, `Assignment(x, y + 1)`,
`For(i, Range(10), list(aug_assign(s, "+", i)))` or
`FunctionDefinition(real, "f", list(x), list(Return(x^2)))`. The nodes
print through the code printers such as
[`ccode()`](https://adamrsophie.github.io/symcompR/reference/code-printers.md)
and
[`fcode()`](https://adamrsophie.github.io/symcompR/reference/code-printers.md).
`aug_assign(lhs, op, rhs)` builds the augmented assignment for the
operator `op`. The objects `real`, `sym_integer`, `complex_`, `bool_`,
`intc`, `float32`, `float64` and the rest are the predefined types;
`value_const` and `pointer_const` are attributes; `none` is the `None`
token.

## Usage

``` r
Token(...)

BreakToken(...)

ContinueToken(...)

NoneToken(...)

none

break_

continue_

String(...)

QuotedString(...)

Comment(...)

Node(...)

Type(...)

IntBaseType(...)

SignedIntType(...)

UnsignedIntType(...)

FloatBaseType(...)

FloatType(...)

ComplexBaseType(...)

ComplexType(...)

intc

intp

int8

int16

int32

int64

uint8

uint16

uint32

uint64

float16

float32

float64

float80

float128

float256

complex64

complex128

untyped

real

sym_integer

complex_

bool_

Attribute(...)

value_const

pointer_const

Variable(...)

Pointer(...)

Element(...)

Declaration(...)

Assignment(lhs, rhs)

AddAugmentedAssignment(lhs, rhs)

SubAugmentedAssignment(lhs, rhs)

MulAugmentedAssignment(lhs, rhs)

DivAugmentedAssignment(lhs, rhs)

ModAugmentedAssignment(lhs, rhs)

aug_assign(lhs, op, rhs)

CodeBlock(...)

For(...)

While(...)

Scope(...)

Stream(...)

stdout_

stderr_

Print(...)

FunctionPrototype(...)

FunctionDefinition(...)

Return(...)

FunctionCall(...)

KeywordFunctionCall(...)

Raise(...)

RuntimeError_(...)
```

## Arguments

- ...:

  The fields of the node, positionally or by name.

- lhs, rhs:

  Left and right hand sides of an assignment.

- op:

  An operator: `"+"`, `"-"`, `"*"`, `"/"` or `"%"`.

## Value

A `sym` node.

## Examples

``` r
x <- Symbol("x"); y <- Symbol("y")
Assignment(x, y + 1)
#> Assignment(x, y + 1)
Variable(x, type = real)
#> Variable(x, type=real)
CodeBlock(Assignment(x, 1), aug_assign(x, "+", y))
#> CodeBlock(
#> Assignment(x, 1),
#> AddAugmentedAssignment(x, y)
#> )
```
