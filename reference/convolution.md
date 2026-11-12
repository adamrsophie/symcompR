# Discrete transforms and convolutions

Exact transforms over symbolic sequences: `sym_fft()` and `sym_ifft()`
(discrete Fourier transform with exact roots of unity), `convolution()`
(linear or `cycle`-length cyclic), `mobius_transform()` and
`inverse_mobius_transform()` over subset lattices, and
`walsh_hadamard_transform()`.

## Usage

``` r
convolution(a, b, cycle = 0)

sym_fft(seq)

sym_ifft(seq)

mobius_transform(seq, subset = TRUE)

inverse_mobius_transform(seq, subset = TRUE)

walsh_hadamard_transform(seq, inverse = FALSE)
```

## Arguments

- a, b:

  Sequences.

- cycle:

  Length of a cyclic convolution (0 for linear).

- seq:

  A list or vector of values.

- subset:

  Sum over subsets (`TRUE`) or supersets.

- inverse:

  Compute the inverse transform.

## Value

A list of expressions.

## Examples

``` r
convolution(1:3, 4:5)
#> [[1]]
#> 4
#> 
#> [[2]]
#> 13
#> 
#> [[3]]
#> 22
#> 
#> [[4]]
#> 15
#> 
sym_fft(list(1, 2, 3, 4))
#> [[1]]
#> 10
#> 
#> [[2]]
#> -2 - 2*I
#> 
#> [[3]]
#> -2
#> 
#> [[4]]
#> -2 + 2*I
#> 
```
