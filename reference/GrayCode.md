# Gray codes

`GrayCode(n, start, rank)` is SymPy's `GrayCode`: the reflected binary
code on `n` bits, positioned at the bit string `start` or at the code of
a given `rank`. The fields `$n`, `$selections` (`2^n`), `$current` and
`$rank` mirror SymPy's. `graycode_next(g, delta)` moves `delta` codes
on; `generate_gray(g, start, rank)` returns the codes from the current
one to the end, and with `lazy = TRUE` an iterator (a function returning
the next code, then `NULL`) during which `graycode_skip(g)` skips the
next code, as SymPy's generator and `skip()` do.
`graycode_unrank(n, rank)` gives the code of a rank.

## Usage

``` r
GrayCode(n, start = NULL, rank = NULL)

graycode_next(g, delta = 1)

graycode_skip(g)

generate_gray(g, start = NULL, rank = NULL, lazy = FALSE)

graycode_unrank(n, rank)

random_bitstring(n)

gray_to_bin(bin_list)

bin_to_gray(bin_list)

get_subset_from_bitstring(super_set, bitstring)

graycode_subsets(gray_code_set)
```

## Arguments

- n:

  Number of bits.

- start:

  A starting bit string.

- rank:

  A rank.

- g:

  A `GrayCode`.

- delta:

  Steps to move.

- lazy:

  Return an iterator instead of the codes.

- bin_list:

  A bit string.

- super_set:

  A vector, list or string of elements.

- bitstring:

  A bit string.

- gray_code_set:

  A vector or list of elements.

## Value

A `GrayCode`, a bit string, or a vector or list.

## Details

`bin_to_gray()` and `gray_to_bin()` convert bit strings,
`random_bitstring(n)` draws one,
`get_subset_from_bitstring(super_set, bitstring)` selects elements and
`graycode_subsets(set)` lists the subsets in Gray code order.

## Examples

``` r
g <- GrayCode(3)
generate_gray(g)
#> [1] "000" "001" "011" "010" "110" "111" "101" "100"
GrayCode(5, start = "10010")$rank
#> [1] 28
graycode_next(GrayCode(6), -1)$current
#> [1] "100000"
bin_to_gray("111")
#> [1] "100"
graycode_subsets(c("a", "b", "c"))
#> [[1]]
#> character(0)
#> 
#> [[2]]
#> [1] "c"
#> 
#> [[3]]
#> [1] "b" "c"
#> 
#> [[4]]
#> [1] "b"
#> 
#> [[5]]
#> [1] "a" "b"
#> 
#> [[6]]
#> [1] "a" "b" "c"
#> 
#> [[7]]
#> [1] "a" "c"
#> 
#> [[8]]
#> [1] "a"
#> 
```
