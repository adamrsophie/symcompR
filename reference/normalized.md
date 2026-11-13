# Use normalized Hadamard gates

`normalized(FALSE)` makes the Hadamard gate's matrix `[[1, 1], [1, -1]]`
instead of the normalized one, as SymPy's `gate.normalized`.

## Usage

``` r
normalized(normalize)
```

## Arguments

- normalize:

  `TRUE` or `FALSE`.

## Value

`NULL`, invisibly.

## Examples

``` r
normalized(TRUE)
```
