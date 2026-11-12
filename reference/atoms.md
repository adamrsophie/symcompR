# Atoms of an expression

Collects the leaves of the expression tree. With `types`, only atoms
whose type is listed are returned (for example `"Symbol"`, `"Integer"`,
`"Rational"`, `"Constant"`, or `"Function"` to collect function
applications).

## Usage

``` r
atoms(e, types = NULL)
```

## Arguments

- e:

  A `sym` expression.

- types:

  Optional character vector of types.

## Value

A list of unique `sym` objects.
