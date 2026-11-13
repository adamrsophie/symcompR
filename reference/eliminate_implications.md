# Eliminate implications and convert to normal forms

`eliminate_implications()` rewrites `Implies`, `Equivalent`, `Xor` and
`ITE` with `And`, `Or` and `Not`. `to_nnf()` gives negation normal form,
`to_cnf()` conjunctive normal form and `to_dnf()` disjunctive normal
form; with `simplify = TRUE` the result is minimized with
[`simplify_logic()`](https://adamrsophie.github.io/symcompR/reference/simplify_logic.md).
`is_cnf()`, `is_dnf()` and `is_nnf()` test the forms.

## Usage

``` r
eliminate_implications(e, form = NULL)

to_nnf(e, simplify = TRUE, form = NULL)

to_cnf(e, simplify = FALSE)

to_dnf(e, simplify = FALSE)

is_nnf(e)

is_cnf(e)

is_dnf(e)
```

## Arguments

- e:

  A boolean expression.

- form:

  For `eliminate_implications()` and `to_nnf()`: `"dnf"` writes a `Xor`
  as a disjunction, otherwise as a conjunction.

- simplify:

  For `to_cnf()` and `to_dnf()`, minimize the result; for `to_nnf()`,
  drop complementary arguments (the default).

## Value

A boolean expression or a logical value.

## Examples

``` r
a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c")
to_cnf(Or(And(a, b), c))
#> (a | c) & (b | c)
to_dnf(And(Or(a, b), c))
#> (a & c) | (b & c)
```
