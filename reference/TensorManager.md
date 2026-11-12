# Commutation groups of tensors

Port of SymPy's `TensorManager`. Each tensor head belongs to a
commutation group: `0` commutes with everything, `1` anticommutes within
itself, `2` commutes only with group `0`.
`TensorManager$set_comm(i, j, c)` sets the commutation (`0` commuting,
`1` anticommuting, `NULL` neither) of the groups labeled `i` and `j`;
`get_comm(i, j)` reads it from group numbers; `comm_symbols2i(i)` gives
the group number of a label (creating it); `set_comms(...)` takes
`list(i, j, c)` triples; `clear()` restores the defaults.

## Usage

``` r
TensorManager
```

## Format

An environment of functions.

## Examples

``` r
TensorManager$set_comm("Gcomm", "GHcomm", 0)
TensorManager$get_comm(TensorManager$comm_symbols2i("Gcomm"),
                       TensorManager$comm_symbols2i("GHcomm"))
#> [1] 0
TensorManager$clear()
```
