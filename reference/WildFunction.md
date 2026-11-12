# Wild functions

`WildFunction(name, nargs)` is a wild symbol that matches any function
application, as SymPy's `WildFunction`; with `nargs` only applications
with that many arguments match.

## Usage

``` r
WildFunction(name, nargs = NULL)
```

## Arguments

- name:

  The name.

- nargs:

  `NULL` (any number of arguments) or allowed argument counts.

## Value

A wild symbol.

## Examples

``` r
x <- Symbol("x"); f <- Function("f")
F <- WildFunction("F")
sym_match(f(x), F)
#> $F
#> f(x)
#> 
sym_match(x + cos(x), x + F)
#> $F
#> cos(x)
#> 
sym_match(f(x), WildFunction("G", nargs = 2))
#> NULL
```
