# Assumption shortcuts

Each returns `TRUE`, `FALSE` or `NA` (unknown), like the SymPy
attributes `is_positive`, `is_real` and so on.

## Usage

``` r
is_positive(e)

is_negative(e)

is_zero(e)

is_nonzero(e)

is_nonnegative(e)

is_nonpositive(e)

is_real(e)

is_complex(e)

is_integer(e)

is_rational(e)

is_irrational(e)

is_even(e)

is_odd(e)

is_prime_expr(e)

is_finite(e)

is_infinite(e)

is_imaginary(e)

is_commutative(e)

is_algebraic(e)

is_transcendental(e)
```

## Arguments

- e:

  A `sym` expression.

## Value

`TRUE`, `FALSE` or `NA`.
