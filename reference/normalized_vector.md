# Normalized vector and nonzero pattern

`normalized_vector(v)` divides a row or column vector by its norm, as
SymPy's `Matrix.normalized()`; `print_nonzero(M, symb)` prints the
pattern of nonzero entries.

## Usage

``` r
normalized_vector(v, iszerofunc = NULL)

print_nonzero(M, symb = "X")
```

## Arguments

- v:

  A vector `symmatrix`.

- iszerofunc:

  The zero test for the norm.

- M:

  A `symmatrix`.

- symb:

  The character marking a nonzero entry.

## Value

A `symmatrix`; `print_nonzero()` returns `M` invisibly.

## Examples

``` r
normalized_vector(Matrix(list(3, 4)))
#> Matrix(2, 1)
#> [3/5]
#> [4/5]
print_nonzero(eye(3))
#> [X  ]
#> [ X ]
#> [  X]
```
