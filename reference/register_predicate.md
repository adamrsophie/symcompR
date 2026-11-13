# Register an assumption handler

Adds an evaluation rule to a predicate, as SymPy's
`Predicate.register(*types)` decorator does. `fn` is called as
`fn(expr, assumptions)` for a unary predicate, or with one argument per
type for a polyadic one, when the arguments have the given SymPy classes
(such as `"Symbol"`, `"Integer"`, `"Add"` or `"MatrixSymbol"`), or
subclasses of them. It returns `TRUE`, `FALSE` or `NA` for unknown.

## Usage

``` r
register_predicate(predicate, types, fn)
```

## Arguments

- predicate:

  A predicate made by `Predicate(name, dispatch = TRUE)`, or one of `Q`.

- types:

  A character vector of SymPy class names, one per argument.

- fn:

  The handler function.

## Value

The predicate, invisibly.

## Examples

``` r
x <- Symbol("x")
P <- Predicate("my_key", dispatch = TRUE)
register_predicate(P, "Symbol", function(expr, assumptions) TRUE)
ask(P(x))
#> [1] TRUE
ask(P(x + 1))
#> [1] NA
```
