# Fundamental discriminant

`extract_fundamental_discriminant(a)` writes an integer `a` that is 0 or
1 modulo 4 as `a = D*F^2` with `D` a fundamental discriminant, returning
the factorizations of `D` and `F` as named lists of exponents, as
SymPy's function of the same name.

## Usage

``` r
extract_fundamental_discriminant(a)
```

## Arguments

- a:

  An integer.

## Value

A list with elements `D` and `F`.

## Examples

``` r
extract_fundamental_discriminant(-4 * 3^5 * 5^3)
#> $D
#> $D$`3`
#> [1] 1
#> 
#> $D$`5`
#> [1] 1
#> 
#> $D$`-1`
#> [1] 1
#> 
#> 
#> $F
#> $F$`2`
#> [1] 1
#> 
#> $F$`3`
#> [1] 2
#> 
#> $F$`5`
#> [1] 1
#> 
#> 
```
