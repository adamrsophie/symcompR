# Random symbols

`RandomSymbol(symbol, pspace)` wraps a symbol as a random variable of
the probability space `pspace`; without a space it belongs to the
abstract space `PSpace()`, as in SymPy.

## Usage

``` r
RandomSymbol(symbol, pspace = NULL)
```

## Arguments

- symbol:

  A symbol or a name.

- pspace:

  A probability space, or `NULL`.

## Value

A random symbol.

## Examples

``` r
a <- RandomSymbol(Symbol("x"))
Variance(a)
#> Variance(x)
```
