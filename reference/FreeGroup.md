# Free groups

`FreeGroup(symbols)` is the free group on the given symbols; groups on
the same symbols are the same object, as in SymPy. `free_group(symbols)`
returns a list whose first element is the group and whose remaining
elements are its generators. `symbols` is a string such as `"x, y"` or
`"x:3"`, a character vector or a list of `Symbol`s. Words multiply with
`*`, are raised to integer powers with `^` and divided with `/`; `<`
compares them in the length lexicographic order of SymPy.
[`fp_order()`](https://adamrsophie.github.io/symcompR/reference/FpGroup.md),
[`fp_is_abelian()`](https://adamrsophie.github.io/symcompR/reference/fp_to_perm_group.md),
[`fp_center()`](https://adamrsophie.github.io/symcompR/reference/fp_to_perm_group.md),
[`fp_elements()`](https://adamrsophie.github.io/symcompR/reference/fp_to_perm_group.md)
and `fp_contains()` work on free groups too, and `fg_slice()` is SymPy's
`F[i:j]`.

## Usage

``` r
FreeGroup(symbols)

free_group(symbols)

xfree_group(symbols)

fg_slice(G, from = 0L, to = G$rank)

fg_index(G, gen)

fp_contains(G, gen)
```

## Arguments

- symbols:

  The generator names.

- G:

  A free group.

- from, to:

  Zero based slice limits, as in Python.

- gen:

  A generator of `G`.

## Value

A `FreeGroup` (an environment), a list, or an integer.

## Examples

``` r
G <- free_group("x, y, z")
F <- G[[1]]; x <- G[[2]]; y <- G[[3]]
fp_str(F)
#> [1] "<free group on the generators (x, y, z)>"
fp_str(x^2 * y * x^-1 * x)
#> [1] "x^2*y"
fp_str(fg_slice(F, 1))
#> [1] "<free group on the generators (y, z)>"
```
