# Coset table options

Sets SymPy's class wide coset table limits: `coset_table_max_limit`, the
default largest number of cosets of an enumeration, and
`max_stack_size`, the deduction stack size at which the Felsch strategy
runs a look ahead.

## Usage

``` r
coset_table_options(coset_table_max_limit = NULL, max_stack_size = NULL)
```

## Arguments

- coset_table_max_limit, max_stack_size:

  New values, or `NULL` to keep the current one.

## Value

The previous values, invisibly.

## Examples

``` r
old <- coset_table_options(max_stack_size = 10)
coset_table_options(max_stack_size = old$max_stack_size)
```
