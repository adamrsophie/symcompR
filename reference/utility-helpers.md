# Utility helpers

`filldedent(s, w)` strips leading and trailing empty lines, removes the
common indentation and refills the text to width `w`, starting with a
newline, as SymPy's `filldedent()`. `capture(func)` returns what
`func()` prints. `threaded(func)` and `xthreaded(func)` return a
function that maps `func` over matrices, lists and the two sides of a
relation; `threaded()` also maps over the terms of a sum. `public(obj)`
returns `obj` (R namespaces are declared in `NAMESPACE`).
`memoize_property(f)` caches `f(x)` per argument. `timed(func)` times
`func()` adaptively and returns the number of runs, the time per run in
seconds, and that time in the unit given last.

## Usage

``` r
filldedent(s, w = 70)

capture(func)

threaded(func)

xthreaded(func)

public(obj)

memoize_property(f)

timed(func, limit = NULL)
```

## Arguments

- s:

  A string.

- w:

  The line width.

- func, f:

  A function.

- obj:

  Any object.

- limit:

  The largest number of runs per timing in `timed()`.

## Value

A string, a function, or a list.

## Examples

``` r
filldedent("
    A long docstring line
    continues here.")
#> [1] "\nA long docstring line continues here."
capture(function() print(Symbol("x")^2))
#> [1] "x^2\n"
d <- threaded(function(e) diff(e, Symbol("x")))
d(list(Symbol("x")^2, sin(Symbol("x"))))
#> [[1]]
#> 2*x
#> 
#> [[2]]
#> cos(x)
#> 
```
