# Trusses

`Truss()` is a plane truss of nodes joined by members. `$add_node()`,
`$remove_node()`, `$add_member()`, `$remove_member()`,
`$change_node_label()` and `$change_member_label()` build it; each takes
any number of tuples, such as `list(label, x, y)` for a node and
`list(label, start, end)` for a member.
`$apply_load(list(node, magnitude, direction))` adds a load with the
direction in degrees, `$apply_support(list(node, type))` a `"pinned"` or
`"roller"` support, and `$solve()` finds the reaction loads and the
internal forces. The properties are `$nodes`, `$node_labels`,
`$node_positions`, `$members`, `$member_lengths`, `$supports`, `$loads`,
`$reaction_loads` and `$internal_forces`. `$draw()` draws the truss with
base R graphics. After `sympy.physics.continuum_mechanics.truss`.

## Usage

``` r
Truss()
```

## Value

A `Truss` object.

## Examples

``` r
t <- Truss()
t$add_node(list("A", 0, 0), list("B", 2, 2), list("C", 3, 0))
t$add_member(list("AB", "A", "B"), list("BC", "B", "C"), list("AC", "A", "C"))
t$apply_load(list("B", Symbol("P"), -90))
t$apply_support(list("A", "pinned"), list("C", "roller"))
t$solve()
t$reaction_loads
#> {R_A_x: 0, R_A_y: P/3, R_C_y: 2*P/3} 
```
