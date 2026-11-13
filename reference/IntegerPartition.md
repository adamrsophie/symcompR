# Integer partitions

`IntegerPartition(partition, integer)` stores a partition of an integer,
as SymPy's `IntegerPartition`. The partition is a vector of positive
integers or a named vector of multiplicities (names are the parts); when
`integer` is given first, as in SymPy, the parts must add up to it. The
fields `$partition`, `$integer` and `$conjugate` mirror SymPy's.
`next_lex()` and `prev_lex()` step through the partitions of the integer
in lexical order, `as_dict()` gives the multiplicities and
`as_ferrers()` draws the Ferrers diagram. Comparisons follow SymPy's
(`<` and `<=` compare the parts from smallest to largest).

## Usage

``` r
IntegerPartition(partition, integer = NULL)

as_dict(p)

as_ferrers(p, char = "#")

prev_lex(p)

next_lex(p)

random_integer_partition(n, seed = NULL)
```

## Arguments

- partition:

  Parts, or the integer when `integer` gives the parts.

- integer:

  The parts, or `NULL`.

- p:

  An `IntegerPartition`.

- char:

  The character for the Ferrers diagram.

- n:

  A positive integer.

- seed:

  A seed for SymPy's `_randint`.

## Value

An `IntegerPartition`; `as_dict()` a named integer vector,
`as_ferrers()` a string and `random_integer_partition()` an integer
vector.

## Details

`random_integer_partition(n, seed)` draws a random partition; `seed`
follows SymPy's `_randint` (an integer seed or a vector of values to
use).

## Examples

``` r
a <- IntegerPartition(c(5, 4, 3, 1, 1))
a
#> [5, 4, 3, 1, 1]
a$conjugate
#> [1] 5 3 3 2 1
next_lex(IntegerPartition(c(3, 1)))
#> [4]
as_dict(a)
#> 5 4 3 1 
#> 1 1 1 2 
cat(as_ferrers(IntegerPartition(c(1, 2, 3))))
#> ###
#> ##
#> #
random_integer_partition(10, seed = c(1, 3, 2, 1, 5, 1))
#> [1] 5 2 1 1 1
```
