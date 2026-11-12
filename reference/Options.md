# Polynomial options

`Options(gens, args)` processes the options of polynomial functions, as
SymPy's `Options`: generators, `domain` (a domain or its name), `field`,
`greedy`, `expand`, `extension`, `modulus` and the others, with SymPy's
defaults.

## Usage

``` r
Options(gens = list(), args = list())
```

## Arguments

- gens:

  Generators.

- args:

  A named list of options.

## Value

A list of options.

## Examples

``` r
Options(list(Symbol("x")), list(domain = "ZZ"))$domain
#> ZZ
```
