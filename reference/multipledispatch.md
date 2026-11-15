# Multiple dispatch

A port of `sympy.multipledispatch`. `Dispatcher(name)` makes a function
that dispatches on the classes of all its arguments;
`dispatcher_add(d, signature, fun)` registers an implementation for a
vector of class names, which may be SymPy class names (`"Integer"`,
`"Expr"`, `"Basic"`), R classes, or `"object"` for any value.
`dispatch(..., fun)` registers `fun` under its name in a namespace and
returns the dispatcher, as the decorator does. An implementation that
raises `MDNotImplementedError` passes the call on to the next most
specific one. `halt_ordering()` and `restart_ordering()` suspend and
resume the ambiguity warnings while many implementations are added.

## Usage

``` r
Dispatcher(name)

dispatcher_add(d, signature, fun)

dispatch(..., fun, fname = NULL, namespace = NULL)

MDNotImplementedError(message = "")

halt_ordering()

restart_ordering()
```

## Arguments

- name:

  The function name.

- d:

  A dispatcher.

- signature:

  A character vector of class names.

- fun:

  The implementation.

- ...:

  Class names of the signature.

- fname:

  The name to register `fun` under.

- namespace:

  An environment that holds the dispatchers.

- message:

  The error message.

## Value

A dispatcher (a function); `dispatcher_add()` returns it invisibly.

## Examples

``` r
f <- dispatch("Integer", fname = "f", fun = function(x) x + 1)
f <- dispatch("Rational", fname = "f", fun = function(x) x * 2)
f(Integer(3))
#> 4
f(Rational(1, 3))
#> 2/3
```
