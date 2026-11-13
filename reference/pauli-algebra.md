# Pauli algebra

`Pauli(i, label)` is the Pauli matrix `sigma_i` (`i` in 1, 2, 3) of
`sympy.physics.paulialgebra` as a noncommutative Hermitian symbol named
`label` followed by `i`. The product of two of them with the same label
is evaluated, and so are their positive integer powers.
`evaluate_pauli_product()` evaluates the products of Pauli symbols in an
expression.

## Usage

``` r
Pauli(i, label = "sigma")

evaluate_pauli_product(arg)
```

## Arguments

- i:

  1, 2 or 3.

- label:

  The name of the symbol.

- arg:

  An expression.

## Value

A symbol or an expression.

## Examples

``` r
s1 <- Pauli(1); s2 <- Pauli(2); s3 <- Pauli(3)
s1 * s2
#> I*sigma3
evaluate_pauli_product(sym_I * s2 * s3)
#> -sigma1
```
