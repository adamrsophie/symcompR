# Discrete transforms and convolutions

Exact transforms over symbolic sequences: `sym_fft()` and `sym_ifft()`
(discrete Fourier transform with exact roots of unity), `convolution()`
(linear or `cycle`-length cyclic), `mobius_transform()` and
`inverse_mobius_transform()` over subset lattices, and
`walsh_hadamard_transform()`.

`convolution(a, b)` is SymPy's `convolution`: the linear convolution of
two sequences, cyclic of length `cycle` when `cycle > 0`. Integer and
rational sequences use `convolution_int()` (Kronecker substitution);
other sequences use `convolution_fft()`, with `dps` decimal digits when
given. `prime` asks for `convolution_ntt()` modulo a prime of the form
`m 2^k + 1`, `dyadic = TRUE` for the dyadic (XOR) convolution
`convolution_fwht()` and `subset = TRUE` for the subset convolution
`convolution_subset()`. At most one of these may be given.

## Usage

``` r
sym_fft(seq)

sym_ifft(seq)

mobius_transform(seq, subset = TRUE)

inverse_mobius_transform(seq, subset = TRUE)

walsh_hadamard_transform(seq, inverse = FALSE)

convolution(
  a,
  b,
  cycle = 0,
  dps = NULL,
  prime = NULL,
  dyadic = NULL,
  subset = NULL
)

convolution_fft(a, b, dps = NULL)

convolution_ntt(a, b, prime)

convolution_fwht(a, b)

convolution_subset(a, b)

convolution_int(a, b)
```

## Arguments

- seq:

  A list or vector of values.

- subset:

  Use the subset convolution.

- inverse:

  Compute the inverse transform.

- a, b:

  Sequences (vectors or lists).

- cycle:

  Length of a cyclic convolution (0 for linear).

- dps:

  Decimal digits for a floating point transform.

- prime:

  A prime modulus for the number theoretic transform.

- dyadic:

  Use the dyadic convolution.

## Value

A list of expressions.

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
convolution(c(1, 2, 3), c(4, 5, 6))
#> [[1]]
#> 4
#> 
#> [[2]]
#> 13
#> 
#> [[3]]
#> 28
#> 
#> [[4]]
#> 27
#> 
#> [[5]]
#> 18
#> 
convolution(c(1, 2, 3), c(4, 5, 6), cycle = 3)
#> [[1]]
#> 31
#> 
#> [[2]]
#> 31
#> 
#> [[3]]
#> 28
#> 
convolution(c(1, 2, 3), c(4, 5, 6), prime = 19 * 2^10 + 1)
#> [[1]]
#> 4
#> 
#> [[2]]
#> 13
#> 
#> [[3]]
#> 28
#> 
#> [[4]]
#> 27
#> 
#> [[5]]
#> 18
#> 
convolution_fwht(c(1, 2, 3), c(4, 5, 6))
#> [[1]]
#> 32
#> 
#> [[2]]
#> 13
#> 
#> [[3]]
#> 18
#> 
#> [[4]]
#> 27
#> 
convolution_int(c(10, -5, 1, 3), c(-5, 6, 7))
#> [[1]]
#> -50
#> 
#> [[2]]
#> 85
#> 
#> [[3]]
#> 35
#> 
#> [[4]]
#> -44
#> 
#> [[5]]
#> 25
#> 
#> [[6]]
#> 21
#> 
```
