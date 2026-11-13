# Manipulate expressions using paths

`EPath(path)` compiles a path and `epath(path, expr, func)` selects or
modifies the parts of an expression (or a nested list of expressions)
that the path names, in the manner of XPath. A path is a sequence of
selectors, each starting with `/`:

## Usage

``` r
EPath(path)

epath(path, expr = NULL, func = NULL, args = NULL, kwargs = NULL)
```

## Arguments

- path:

  A path string or an `EPath`.

- expr:

  An expression or a list of expressions.

- func:

  A function applied to the selected parts.

- args, kwargs:

  Extra positional and named arguments for `func`, as lists.

## Value

An `EPath`, a list of selected parts, or the modified expression.

## Details

- `/*` selects all arguments.

- `/[0]`, `/[1:5]`, `/[1:5:2]` select by Python index or slice (0-based,
  negative indices count from the end).

- `/Symbol` or `/list|tuple` select by SymPy class name, including
  parent classes such as `TrigonometricFunction`, `Function` or `Atom`.

- `/__iter__?` selects by attribute.

The arguments of a sum and a product are visited in printing order.

With no `expr`, `epath()` returns the compiled path. With no `func`, it
returns the list of selected parts; otherwise it returns the expression
with `func` applied to each selected part. An `EPath` object has the
methods `$select(expr)` and `$apply(expr, func, args, kwargs)`.

## Examples

``` r
x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); t <- Symbol("t")
expr <- list(list(list(x, 1), 2), list(list(3, y), z))
epath("/*/[0]/Symbol", expr)
#> [[1]]
#> x
#> 
#> [[2]]
#> y
#> 
epath("/*/[0]/Symbol", expr, function(e) e^2)
#> [[1]]
#> [[1]][[1]]
#> [[1]][[1]][[1]]
#> x^2
#> 
#> [[1]][[1]][[2]]
#> [1] 1
#> 
#> 
#> [[1]][[2]]
#> [1] 2
#> 
#> 
#> [[2]]
#> [[2]][[1]]
#> [[2]][[1]][[1]]
#> [1] 3
#> 
#> [[2]][[1]][[2]]
#> y^2
#> 
#> 
#> [[2]][[2]]
#> z
#> 
#> 
e2 <- t + sin(x + 1) + cos(x + y + sym_E)
epath("/*/*/Symbol", e2)
#> [[1]]
#> x
#> 
#> [[2]]
#> x
#> 
#> [[3]]
#> y
#> 
epath("/*/*/Symbol", e2, function(e) 2 * e)
#> t + sin(2*x + 1) + cos(2*x + 2*y + E)
```
