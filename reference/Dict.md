# Symbolic dictionaries

`Dict(keys, values)` is SymPy's `Dict`, an immutable mapping between
expressions stored as `(key, value)` tuples. `keys` may also be a list
of `list(key, value)` pairs. A repeated key keeps its last value.
`dict_get()` looks a key up, and `dict_keys()`, `dict_values()` and
`dict_items()` list the contents in printing order.

## Usage

``` r
Dict(keys, values = NULL)

dict_get(d, key, default = NULL)

dict_keys(d)

dict_values(d)

dict_items(d)
```

## Arguments

- keys:

  Keys, or a list of key and value pairs.

- values:

  Values, one per key.

- d:

  A `Dict`.

- key:

  A key.

- default:

  Returned when the key is absent.

## Value

A `Dict`, an expression, or a list.

## Examples

``` r
x <- Symbol("x"); y <- Symbol("y")
d <- Dict(list(x, y), list(1, 2))
d
#> {x: 1, y: 2}
dict_get(d, x)
#> 1
subs(Dict(list(x), list(y^2)), y, 3)
#> {x: 9}
```
