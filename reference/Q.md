# Assumption predicates

`Q` holds the predicates of SymPy's `Q` object, such as `Q$positive`,
`Q$integer`, `Q$invertible` or `Q$eq`. Calling one gives an applied
predicate, a boolean expression that prints as in SymPy
(`Q.positive(x)`) and can be combined with `&`, `|` and `!` and passed
to [`ask()`](https://adamrsophie.github.io/symcompR/reference/ask.md),
[`refine()`](https://adamrsophie.github.io/symcompR/reference/refine.md)
and
[`assuming()`](https://adamrsophie.github.io/symcompR/reference/assumptions-context.md).
Asking for a predicate that `Q` does not have raises an
`AttributeError`; `Q$name <- p` adds a predicate made with
[`Predicate()`](https://adamrsophie.github.io/symcompR/reference/Predicate.md),
and `Q$name <- NULL` removes it again.

## Usage

``` r
Q
```

## Format

An object of class `AssumptionKeys`.

## Details

The predicates are `hermitian`, `antihermitian`, `real`,
`extended_real`, `imaginary`, `complex`, `algebraic`, `transcendental`,
`integer`, `noninteger`, `rational`, `irrational`, `finite`, `infinite`,
`positive_infinite`, `negative_infinite`, `positive`, `negative`,
`zero`, `extended_positive`, `extended_negative`, `nonzero`,
`nonpositive`, `nonnegative`, `extended_nonzero`,
`extended_nonpositive`, `extended_nonnegative`, `even`, `odd`, `prime`,
`composite`, `commutative`, `is_true`, the matrix predicates
`symmetric`, `invertible`, `orthogonal`, `unitary`, `positive_definite`,
`upper_triangular`, `lower_triangular`, `diagonal`, `fullrank`,
`square`, `integer_elements`, `real_elements`, `complex_elements`,
`singular`, `normal`, `triangular`, `unit_triangular`, and the binary
relations `eq`, `ne`, `gt`, `ge`, `lt` and `le`.

## Examples

``` r
x <- Symbol("x")
Q$positive(x)
#> Q.positive(x)
Q$positive(x) & !Q$integer(x)
#> Q.positive(x) & ~Q.integer(x)
Q$eq(x, 1)
#> Q.eq(x, 1)
srepr(Q$positive(x))
#> [1] "AppliedPredicate(Q.positive, Symbol('x'))"
ask(Q$positive(x^2 + 1), Q$real(x))
#> [1] TRUE
```
