# Rewriting systems for finitely presented groups

`RewritingSystem(G)` builds the rewriting system SymPy attaches to an
`FpGroup` (also reachable with `fp_rewriting_system(G)`), with its rules
(`R$rules`, an ordered dictionary of left and right hand sides) and its
reduction automaton. `rws_reduce()` reduces a word with the rules,
`rws_reduce_using_automaton()` with the automaton,
`rws_make_confluent()` runs the Knuth Bendix completion,
`rws_is_confluent()` tests confluence, `rws_add_rule()` adds a rule (and
the rules its overlaps imply) and `rws_set_max()` sets the largest
number of rules. `rws_rules()` lists the rules as pairs.

## Usage

``` r
RewritingSystem(group)

fp_rewriting_system(G)

rws_rules(R)

rws_set_max(R, n)

rws_is_confluent(R)

rws_add_rule(R, w1, w2, check = FALSE)

rws_make_confluent(R, check = FALSE)

rws_reduce(R, word, exclude = NULL)

rws_reduce_using_automaton(R, word)
```

## Arguments

- group, G:

  An `FpGroup`.

- R:

  A rewriting system.

- n:

  The largest number of rules.

- check:

  Only check, without adding rules.

- word, w1, w2:

  Words.

- exclude:

  A rule left hand side to skip.

## Value

A `RewritingSystem` (an environment), a word, a logical or a list.

## Examples

``` r
g <- free_group("a, b")
a <- g[[2]]; b <- g[[3]]
G <- FpGroup(g[[1]], list(a * b * a^-1 * b^-1))
R <- fp_rewriting_system(G)
rws_is_confluent(R)
#> [1] TRUE
fp_str(fp_reduce(G, b^-1 * a))
#> [1] "a*b^-1"
```
