# Free groups and finitely presented groups

``` r

library(symcompR)
```

This vignette covers the port of `sympy.combinatorics.free_groups`,
`fp_groups`, `coset_table`, `rewritingsystem`, `homomorphisms`, and
`pc_groups`. Permutation groups have their own vignette,
[`vignette("permutation-groups", package = "symcompR")`](https://adamrsophie.github.io/symcompR/articles/permutation-groups.md).

## Free groups and words

[`free_group()`](https://adamrsophie.github.io/symcompR/reference/FreeGroup.md)
returns the free group followed by its generators.
[`FreeGroup()`](https://adamrsophie.github.io/symcompR/reference/FreeGroup.md)
takes a character vector of generator names, and
[`xfree_group()`](https://adamrsophie.github.io/symcompR/reference/FreeGroup.md)
returns the group and a list of generators. Words simplify automatically
when they are multiplied.

``` r

fg <- free_group("x, y, z")
F <- fg[[1]]; x <- fg[[2]]; y <- fg[[3]]; z <- fg[[4]]
F
#> <free group on the generators (x, y, z)>
FreeGroup(c("a", "b"))$rank
#> [1] 2
xfree_group("u, v")[[2]]
#> [[1]]
#> u 
#> 
#> [[2]]
#> v 
#> 
#> attr(,"class")
#> [1] "py_tuple"
fg_index(F, y)
#> [1] 1
fg_slice(F, 1, 3)
#> <free group on the generators (y, z)>
fp_contains(F, x * y)
#> [1] TRUE
word_prod(F, list(x, y, x^-1))
#> x*y*x^-1
fp_str(word_inverse(x * y * x^-1))
#> [1] "x*y^-1*x^-1"
fp_str(word_commutator(x * y, y))
#> [1] "y^-1*x^-1*y^-1*x*y^2"
```

The representation of a word as syllables (maximal powers of a single
generator):

``` r

w <- x^5 * y * x^2 * y^-4 * x
w
#> x^5*y*x^2*y^-4*x
word_len(w)
#> [1] 13
number_syllables(w)
#> [1] 5
word_array_form(w)
#> [[1]]
#> [[1]]
#> x
#> 
#> [[2]]
#> [1] 5
#> 
#> attr(,"class")
#> [1] "py_tuple"
#> 
#> [[2]]
#> [[1]]
#> y
#> 
#> [[2]]
#> [1] 1
#> 
#> attr(,"class")
#> [1] "py_tuple"
#> 
#> [[3]]
#> [[1]]
#> x
#> 
#> [[2]]
#> [1] 2
#> 
#> attr(,"class")
#> [1] "py_tuple"
#> 
#> [[4]]
#> [[1]]
#> y
#> 
#> [[2]]
#> [1] -4
#> 
#> attr(,"class")
#> [1] "py_tuple"
#> 
#> [[5]]
#> [[1]]
#> x
#> 
#> [[2]]
#> [1] 1
#> 
#> attr(,"class")
#> [1] "py_tuple"
#> 
#> attr(,"class")
#> [1] "py_tuple"
word_ext_rep(x * z)
#> [[1]]
#> x
#> 
#> [[2]]
#> [1] 1
#> 
#> [[3]]
#> z
#> 
#> [[4]]
#> [1] 1
#> 
#> attr(,"class")
#> [1] "py_tuple"
word_letter_form(x^3 * y^-1)
#> [[1]]
#> x
#> 
#> [[2]]
#> x
#> 
#> [[3]]
#> x
#> 
#> [[4]]
#> -y
#> 
#> attr(,"class")
#> [1] "py_tuple"
word_letter_form_elm(x^2 * z^-1)
#> [[1]]
#> x 
#> 
#> [[2]]
#> x 
#> 
#> [[3]]
#> z^-1
exponent_syllable(w, 2)
#> [1] 2
generator_syllable(w, 3)
#> y
exponent_sum(w, x)
#> [1] 8
generator_count(w, y)
#> [1] 5
```

Subwords use the letter positions, counted from 0:

``` r

subword(w, 2, 6)
#> x^3*y
sub_syllables(w, 1, 2)
#> y
substituted_word(w, 3, 9, z)
#> x^3*z*y^-3*x
subword_index(x^2 * y^-4 * x, y^-2)
#> [1] 2
cyclic_subword(x * y * z, 4, 6)
#> y*z
word_contains(w, y)
#> [1] TRUE
word_index(x^2 * y * x, y)
#> [1] 2
contains_generators(x^-1 * y^3)
#> [[1]]
#> x 
#> 
#> [[2]]
#> y 
#> 
#> attr(,"class")
#> [1] "py_set"
is_dependent(x^2 * y^-4 * x, y^2)
#> [1] TRUE
is_independent(x^2 * y^-4 * x, z)
#> [1] TRUE
power_of(x^4, x^-1)
#> [1] TRUE
is_identity_word(x^2 * x^-2)
#> [1] TRUE
word_order(x^0)
#> 1
word_order(x)
#> oo
```

Cyclic reduction removes a generator and its inverse from the two ends,
and the cyclic conjugates of a word are the cyclic shifts of its
letters.

``` r

cyclic_reduction(x^2 * y * x^-1)
#> x*y
identity_cyclic_reduction(x^2 * y * x^-1)
#> x*y
is_cyclically_reduced(x * y * x)
#> [1] TRUE
is_cyclic_conjugate(x * y * z, z * x * y)
#> [1] TRUE
cyclic_conjugates(x * y * x^-1)
#> [[1]]
#> x*y*x^-1 
#> 
#> [[2]]
#> y 
#> 
#> attr(,"class")
#> [1] "py_set"
```

Eliminating a generator replaces each of its occurrences by another
word:

``` r

eliminate_word(x^2 * y^3 * x^-1 * y, x, y)
#> y^5
eliminate_words(x^2 * y^-3 * x^-1 * y, list(y, x^2))
#> x^-1
eliminate_words(x^2 * y^-3 * x^-1 * y, word_dict(list(y, x^2), list(z, y)))
#> z^-2*x^-1*z
```

## Finitely presented groups

[`FpGroup()`](https://adamrsophie.github.io/symcompR/reference/FpGroup.md)
builds a group from a free group and a list of relators.
[`fp_group()`](https://adamrsophie.github.io/symcompR/reference/FpGroup.md)
and
[`xfp_group()`](https://adamrsophie.github.io/symcompR/reference/FpGroup.md)
are the functional forms. The Todd-Coxeter coset enumeration finds the
order when it is finite. The presentation below is the dihedral group of
order 8.

``` r

fg2 <- free_group("x, y")
F2 <- fg2[[1]]; x <- fg2[[2]]; y <- fg2[[3]]
D4 <- FpGroup(F2, list(x^4, y^2, (x * y)^2))
fp_order(D4)
#> 8
fp_order(fp_group(F2, list(x^2, y^3, (x * y)^5))[[1]])
#> 60
fp_order(xfp_group(F2, list(x^2, y^2, (x * y)^3))[[1]])
#> 6
fp_is_abelian(D4)
#> [1] FALSE
fp_is_cyclic(D4)
#> [1] FALSE
fp_is_nilpotent(D4)
#> [1] TRUE
fp_is_solvable(FpGroup(F2, list(x^3, y^2, (x * y)^5)))
#> [1] FALSE
fp_abelian_invariants(D4)
#> [[1]]
#> [1] 2
#> 
#> [[2]]
#> [1] 2
fp_center(D4)
#> [[1]]
#> x^2
length(fp_elements(D4))
#> [1] 8
fp_index(D4, list(x))
#> [1] 2
fp_most_frequent_generator(FpGroup(F2, list(x^3, y^5, (x * y)^2)))
#> y
fps_contains(FpSubgroup(D4, list(x, y)), fp_random(D4))
#> [1] TRUE
```

Series and subgroups are computed through an isomorphic permutation
group
([`fp_to_perm_group()`](https://adamrsophie.github.io/symcompR/reference/fp_to_perm_group.md)):

``` r

group_order(fp_to_perm_group(D4)[[1]])
#> [1] 8
vapply(fp_composition_series(D4), group_order, 0)
#> [1] 8 4 2 1
length(fp_derived_series(D4))
#> [1] 3
length(fp_lower_central_series(D4))
#> [1] 3
fp_derived_subgroup(D4)
#> [[1]]
#> x^2
length(fp_centralizer(D4, list(x)))
#> [1] 3
fps_order(FpSubgroup(D4, fp_normal_closure(D4, list(y))))
#> 4
fp_quotient(FpGroup(F2, list(x^4, y^2)), list(x^2))$relators
#> [[1]]
#> x^4 
#> 
#> [[2]]
#> y^2 
#> 
#> [[3]]
#> x^2
```

[`FpSubgroup()`](https://adamrsophie.github.io/symcompR/reference/FpSubgroup.md)
is a subgroup given by generating words, and
[`fp_subgroup()`](https://adamrsophie.github.io/symcompR/reference/FpGroup.md)
returns a presentation of a subgroup as a new `FpGroup`.

``` r

H <- FpSubgroup(D4, list(x^2, y))
fps_order(H)
#> 4
fps_contains(H, x^2 * y)
#> [1] TRUE
fps_contains(H, x)
#> [1] FALSE
fp_order(fps_to_FpGroup(H))
#> 4
K <- fp_subgroup(D4, list(x))
fp_order(K)
#> 4
```

Word problems are solved with a rewriting system:

``` r

A4 <- FpGroup(F2, list(x^3, y^3, (x * y)^2))
fp_equals(A4, x * y, y * x)
#> NULL
fp_reduce(A4, x * y * x * y * x)
#> x
fp_make_confluent(A4)
```

## Coset enumeration

[`fp_coset_enumeration()`](https://adamrsophie.github.io/symcompR/reference/FpGroup.md)
runs the Todd-Coxeter procedure for the cosets of a subgroup. The
relator based (HLT) strategy is
[`coset_enumeration_r()`](https://adamrsophie.github.io/symcompR/reference/CosetTable.md)
and the coset table based (Felsch) strategy is
[`coset_enumeration_c()`](https://adamrsophie.github.io/symcompR/reference/CosetTable.md).
[`modified_coset_enumeration_r()`](https://adamrsophie.github.io/symcompR/reference/CosetTable.md)
also records how each coset was reached. Each row of the table lists the
images of a coset under x, x^{-1}, y, y^{-1}.

``` r

G <- FpGroup(F2, list(x^4, y^2, x * y * x^-1 * y))
fp_coset_table(G, list(x))
#> [[1]]
#> [1] 0 0 1 1
#> 
#> [[2]]
#> [1] 1 1 0 0
C <- fp_coset_enumeration(G, list(x * y))
C$table
#> [[1]]
#> [1] 1 1 1 1
#> 
#> [[2]]
#> [1] 0 0 0 0
coset_enumeration_r(G, list(x))$table
#> [[1]]
#> [1] 0 0 1 1
#> 
#> [[2]]
#> [1] 1 1 0 0
coset_enumeration_c(G, list(x))$table
#> [[1]]
#> [1] 0 0 1 1
#> 
#> [[2]]
#> [1] 1 1 0 0
modified_coset_enumeration_r(G, list(x))$table
#> [[1]]
#> [1] 0 0 1 1
#> 
#> [[2]]
#> [1] 1 1 0 0
coset_table_options(max_stack_size = 500)
```

[`CosetTable()`](https://adamrsophie.github.io/symcompR/reference/CosetTable.md)
creates an empty table, and the `ct_` functions are the building blocks
of the enumeration. The example enumerates the cosets of \langle x
\rangle in the group with relators x^2, y^3, (xy)^3 step by step.

``` r

f <- FpGroup(F2, list(x^2, y^3, (x * y)^3))
C <- CosetTable(f, list(x))
ct_A(C)
#> [[1]]
#> x 
#> 
#> [[2]]
#> x^-1 
#> 
#> [[3]]
#> y 
#> 
#> [[4]]
#> y^-1
ct_n(C)
#> [1] 1
ct_scan_and_fill(C, 0, x)
C$table
#> [[1]]
#> [1]  0  0 NA NA
ct_define(C, 0, y)
C$table
#> [[1]]
#> [1]  0  0  1 NA
#> 
#> [[2]]
#> [1] NA NA NA  0
ct_scan_and_fill(C, 0, y^3)
ct_is_complete(C)
#> [1] FALSE
for (alpha in 0:3) for (r in f$relators)
  if (alpha < length(C$table) && C$p[alpha + 1] == alpha) ct_scan_and_fill(C, alpha, r)
ct_is_complete(C)
#> [1] TRUE
ct_compress(C)
ct_standardize(C)
C$table
#> [[1]]
#> [1] 0 0 1 2
#> 
#> [[2]]
#> [1] 3 3 2 0
#> 
#> [[3]]
#> [1] 4 4 0 1
#> 
#> [[4]]
#> [1] 1 1 4 5
#> 
#> [[5]]
#> [1] 2 2 5 3
#> 
#> [[6]]
#> [1] 5 5 3 4
ct_omega(C)
#> [1] 0 1 2 3 4 5
lapply(0:3, function(i) ct_coset_representative(C, i))
#> [[1]]
#> <identity> 
#> 
#> [[2]]
#> y 
#> 
#> [[3]]
#> y^-1 
#> 
#> [[4]]
#> y*x^-1
ct_copy(C)$table
#> [[1]]
#> [1] 0 0 1 2
#> 
#> [[2]]
#> [1] 3 3 2 0
#> 
#> [[3]]
#> [1] 4 4 0 1
#> 
#> [[4]]
#> [1] 1 1 4 5
#> 
#> [[5]]
#> [1] 2 2 5 3
#> 
#> [[6]]
#> [1] 5 5 3 4
```

The other low level steps are
[`ct_define_c()`](https://adamrsophie.github.io/symcompR/reference/CosetTable.md),
[`ct_scan()`](https://adamrsophie.github.io/symcompR/reference/CosetTable.md),
[`ct_scan_c()`](https://adamrsophie.github.io/symcompR/reference/CosetTable.md),
[`ct_scan_check()`](https://adamrsophie.github.io/symcompR/reference/CosetTable.md),
[`ct_scan_and_fill_c()`](https://adamrsophie.github.io/symcompR/reference/CosetTable.md),
[`ct_look_ahead()`](https://adamrsophie.github.io/symcompR/reference/CosetTable.md),
[`ct_process_deductions()`](https://adamrsophie.github.io/symcompR/reference/CosetTable.md),
[`ct_process_deductions_check()`](https://adamrsophie.github.io/symcompR/reference/CosetTable.md),
[`ct_switch()`](https://adamrsophie.github.io/symcompR/reference/CosetTable.md),
[`ct_merge()`](https://adamrsophie.github.io/symcompR/reference/CosetTable.md),
[`ct_rep()`](https://adamrsophie.github.io/symcompR/reference/CosetTable.md),
[`ct_coincidence()`](https://adamrsophie.github.io/symcompR/reference/CosetTable.md),
[`ct_coincidence_c()`](https://adamrsophie.github.io/symcompR/reference/CosetTable.md),
[`ct_conjugates()`](https://adamrsophie.github.io/symcompR/reference/CosetTable.md),
and the modified forms
[`ct_modified_define()`](https://adamrsophie.github.io/symcompR/reference/CosetTable.md),
[`ct_modified_scan()`](https://adamrsophie.github.io/symcompR/reference/CosetTable.md),
[`ct_modified_scan_and_fill()`](https://adamrsophie.github.io/symcompR/reference/CosetTable.md),
[`ct_modified_merge()`](https://adamrsophie.github.io/symcompR/reference/CosetTable.md),
[`ct_modified_rep()`](https://adamrsophie.github.io/symcompR/reference/CosetTable.md),
and
[`ct_modified_coincidence()`](https://adamrsophie.github.io/symcompR/reference/CosetTable.md),
which keep track of the words that connect cosets.

``` r

C2 <- CosetTable(f, list())
ct_define_c(C2, 0, x)
ct_scan_c(C2, 0, x^2)
ct_scan_check(C2, 0, x^2)
#> [1] TRUE
ct_scan(C2, 0, y)
ct_scan_and_fill_c(C2, 0, y^3)
C2$table
#> [[1]]
#> [1] 1 1 0 0
#> 
#> [[2]]
#> [1]  0  0 NA NA
ct_look_ahead(C2)
ct_rep(C2, 1)
#> [1] 1
length(ct_conjugates(C2, f$relators))
#> [1] 4
Cm <- CosetTable(f, list(x))
ct_modified_define(Cm, 0, y)
ct_modified_scan(Cm, 0, y^3, Cm$grp$identity)
ct_modified_scan_and_fill(Cm, 0, x * y, Cm$grp$identity)
Cm$table
#> [[1]]
#> [1]  2 NA  1  2
#> 
#> [[2]]
#> [1] NA NA NA  0
#> 
#> [[3]]
#> [1] NA  0  0 NA
ct_modified_rep(Cm, 1)
```

When two cosets turn out to be equal,
[`ct_coincidence()`](https://adamrsophie.github.io/symcompR/reference/CosetTable.md)
(or
[`ct_coincidence_c()`](https://adamrsophie.github.io/symcompR/reference/CosetTable.md)
for the coset table based strategy) identifies them,
[`ct_merge()`](https://adamrsophie.github.io/symcompR/reference/CosetTable.md)
records the merge of two classes, and
[`ct_switch()`](https://adamrsophie.github.io/symcompR/reference/CosetTable.md)
swaps two cosets. The modified versions also carry the connecting words.
A compact printer shows each table row and the class representatives
`C$p`.

``` r

rows <- function(C) paste0(vapply(C$table, function(r)
  paste(ifelse(is.na(r), "-", r), collapse = " "), ""), collapse = " | ")
C3 <- CosetTable(f, list())
ct_define(C3, 0, x); ct_define(C3, 0, y); ct_define(C3, 1, y)
rows(C3)
#> [1] "1 - 2 - | - 0 3 - | - - - 0 | - - - 1"
ct_coincidence(C3, 1, 2)
rows(C3)
#> [1] "1 - 1 - | - 0 3 0 | - - - 0 | - - - 1"
C3$p
#> [1] 0 1 1 3
C4 <- CosetTable(f, list())
ct_define_c(C4, 0, x); ct_define_c(C4, 0, y)
ct_coincidence_c(C4, 1, 2)
rows(C4)
#> [1] "1 - 1 - | - 0 - 0 | - - - 0"
C5 <- CosetTable(f, list())
ct_define(C5, 0, x); ct_define(C5, 0, y)
unlist(ct_merge(C5, 1, 2, list()))
#> [1] 2
C5$p
#> [1] 0 1 1
C7 <- CosetTable(f, list())
ct_define(C7, 0, x); ct_define(C7, 0, y)
ct_switch(C7, 1, 2)
rows(C7)
#> [1] "2 - 1 - | - - - 0 | - 0 - -"
C6 <- CosetTable(f, list())
ct_modified_define(C6, 0, x); ct_modified_define(C6, 0, y)
unlist(ct_modified_merge(C6, 1, 2, C6$grp$identity, list()))
#> [1] 2
C8 <- CosetTable(f, list())
ct_modified_define(C8, 0, x); ct_modified_define(C8, 0, y)
ct_modified_coincidence(C8, 1, 2, C8$grp$identity)
rows(C8)
#> [1] "1 - 1 - | - 0 - 0 | - - - 0"
```

[`low_index_subgroups()`](https://adamrsophie.github.io/symcompR/reference/low_index_subgroups.md)
finds all subgroups up to a given index, and
[`first_in_class()`](https://adamrsophie.github.io/symcompR/reference/low_index_subgroups.md)
is the test it uses to discard duplicate tables.

``` r

L <- low_index_subgroups(FpGroup(F2, list(x^2, y^3, (x * y)^4)), 4)
length(L)
#> [1] 4
first_in_class(L[[2]])
#> [1] TRUE
```

## Presentations of subgroups

[`reidemeister_presentation()`](https://adamrsophie.github.io/symcompR/reference/reidemeister_presentation.md)
computes a presentation of a subgroup, and
[`simplify_presentation()`](https://adamrsophie.github.io/symcompR/reference/simplify_presentation.md)
removes redundant generators and relators with the elimination
techniques
[`elimination_technique_1()`](https://adamrsophie.github.io/symcompR/reference/simplify_presentation.md)
and
[`elimination_technique_2()`](https://adamrsophie.github.io/symcompR/reference/simplify_presentation.md).

``` r

f5 <- FpGroup(F2, list(x^3, y^5, (x * y)^2))
reidemeister_presentation(f5, list(y))
#> [[1]]
#> [[1]]
#> y_0 
#> 
#> attr(,"class")
#> [1] "py_tuple"
#> 
#> [[2]]
#> [[1]]
#> y_0^5 
#> 
#> attr(,"class")
#> [1] "py_tuple"
#> 
#> attr(,"class")
#> [1] "py_tuple"
reidemeister_presentation(FpGroup(F2, list(x^2, y^2)), list(x * y))
#> [[1]]
#> [[1]]
#> y_0 
#> 
#> attr(,"class")
#> [1] "py_tuple"
#> 
#> [[2]]
#> list()
#> attr(,"class")
#> [1] "py_tuple"
#> 
#> attr(,"class")
#> [1] "py_tuple"
simplify_presentation(FpGroup(F2, list(x^2 * y^-1, y^3)))
#> <fp group on the generators (x, y)>
elimination_technique_1(list(x, y), list(x * y^2, y^3), F2$identity)
#> [[1]]
#> [[1]][[1]]
#> y 
#> 
#> 
#> [[2]]
#> [[2]][[1]]
#> y^3 
#> 
#> 
#> attr(,"class")
#> [1] "py_tuple"
elimination_technique_2(FpGroup(F2, list(x^2 * y^-1, y^3)))
#> [[1]]
#> NULL
#> 
#> [[2]]
#> list()
#> 
#> attr(,"class")
#> [1] "py_tuple"
```

[`perm_group_presentation()`](https://adamrsophie.github.io/symcompR/reference/perm_group_presentation.md)
goes the other way and finds a presentation of a permutation group.

``` r

P <- perm_group_presentation(DihedralGroup(3))
fp_order(P)
#> 6
```

## Rewriting systems

[`fp_rewriting_system()`](https://adamrsophie.github.io/symcompR/reference/RewritingSystem.md)
and
[`RewritingSystem()`](https://adamrsophie.github.io/symcompR/reference/RewritingSystem.md)
return the rewriting system of a group. The Knuth-Bendix completion
makes it confluent.

``` r

fg3 <- free_group("a, b")
a <- fg3[[2]]; b <- fg3[[3]]
Ab <- FpGroup(fg3[[1]], list(a * b * a^-1 * b^-1))
R <- fp_rewriting_system(Ab)
rws_is_confluent(R)
#> [1] TRUE
rws_rules(R)
#> [[1]]
#> [[1]]
#> b*a 
#> 
#> [[2]]
#> a*b 
#> 
#> attr(,"class")
#> [1] "py_tuple"
#> 
#> [[2]]
#> [[1]]
#> b^-1*a^-1 
#> 
#> [[2]]
#> a^-1*b^-1 
#> 
#> attr(,"class")
#> [1] "py_tuple"
#> 
#> [[3]]
#> [[1]]
#> b*a^-1 
#> 
#> [[2]]
#> a^-1*b 
#> 
#> attr(,"class")
#> [1] "py_tuple"
#> 
#> [[4]]
#> [[1]]
#> b^-1*a 
#> 
#> [[2]]
#> a*b^-1 
#> 
#> attr(,"class")
#> [1] "py_tuple"
rws_reduce(R, b * a)
#> a*b
rws_reduce_using_automaton(R, b^-1 * a)
#> a*b^-1
rws_make_confluent(R)
#> [1] TRUE
R2 <- RewritingSystem(FpGroup(fg3[[1]], list(a^3, b^3, (a * b)^2)))
rws_set_max(R2, 20)
rws_add_rule(R2, a^3, fg3[[1]]$identity)
#> list()
#> attr(,"class")
#> [1] "py_set"
rws_is_confluent(R2)
#> [1] FALSE
```

## Homomorphisms

[`homomorphism()`](https://adamrsophie.github.io/symcompR/reference/GroupHomomorphism.md)
defines a homomorphism by the images of the generators and returns a
`GroupHomomorphism` object. The `hom_` functions query it.

``` r

A4p <- AlternatingGroup(4)
T <- homomorphism(A4, A4p, list(x, y),
                  list(Permutation(0, 1, 2, size = 4), Permutation(1, 2, 3)))
hom_is_injective(T)
#> [1] TRUE
hom_is_surjective(T)
#> [1] TRUE
hom_is_isomorphism(T)
#> [1] TRUE
hom_is_trivial(T)
#> [1] FALSE
group_order(hom_image(T))
#> [1] 12
fps_order(hom_kernel(T))
#> 1
hom_apply(T, x * y)
#> (0 2)(1 3)
hom_invert(T, Permutation(0, 1, 2, size = 4))
#> x
is_isomorphic(A4, A4p)
#> [1] TRUE
group_isomorphism(A4, A4p)[[1]]
#> [1] TRUE
```

For permutation groups, a homomorphism can be restricted, composed, and
inverted on subgroups:

``` r

D <- DihedralGroup(4)
Tr <- homomorphism(D, D, D$gens, D$gens)
V <- PermutationGroup(list(Permutation(list(c(0, 2), c(1, 3)))))
hom_images(hom_restrict_to(Tr, V))
#> $k
#> [1] "p2,3,0,1"
#> 
#> $key
#> $key[[1]]
#> (0 2)(1 3)
#> 
#> 
#> $v
#> $v[[1]]
#> (0 2)(1 3)
#> 
#> 
#> attr(,"class")
#> [1] "py_od"
group_order(hom_invert_subgroup(Tr, V))
#> [1] 2
T2 <- homomorphism(D, D, D$gens, list(D$gens[[1]]^2, D$gens[[2]]))
group_order(hom_image(hom_compose(Tr, T2)))
#> [1] 4
Tf <- homomorphism(Ab, FpGroup(fg3[[1]], list(a^2, b^2, (a * b)^2)),
                   list(a, b), list(a, a))
length(hom_factor(Tf))
#> [1] 2
```

[`orbit_homomorphism()`](https://adamrsophie.github.io/symcompR/reference/orbit_homomorphism.md)
maps a group onto its action on an orbit, and
[`block_homomorphism()`](https://adamrsophie.github.io/symcompR/reference/orbit_homomorphism.md)
onto its action on a block system.

``` r

S4 <- SymmetricGroup(4)
To <- orbit_homomorphism(S4, list(0, 1, 2, 3))
group_order(hom_image(To))
#> [1] 24
Tb <- block_homomorphism(DihedralGroup(4), c(0, 1, 0, 1))
group_order(hom_image(Tb))
#> [1] 2
group_order(hom_kernel(Tb))
#> [1] 4
```

## Polycyclic groups

A solvable permutation group has a polycyclic presentation.
[`perm_polycyclic_group()`](https://adamrsophie.github.io/symcompR/reference/PolycyclicGroup.md)
computes it and returns a `PolycyclicGroup` with a collector
(`Collector`) that rewrites words into collected normal form.

``` r

P3 <- perm_polycyclic_group(CyclicGroup(9))
pc_length(P3)
#> [1] 1
pc_is_prime_order(P3)
#> [1] FALSE
K3 <- P3$collector
x0 <- K3$free_group$generators[[1]]
collected_word(K3, x0^13)
#> x0^4
minimal_uncollected_subword(K3, x0^-2)
#> [[1]]
#> [1]  1 -2
pc_presentation(K3)
#> [[1]]
#> [[1]]
#> x0^9 
#> 
#> [[2]]
#> list()
#> attr(,"class")
#> [1] "py_tuple"
#> 
#> attr(,"class")
#> [1] "py_tuple"
pc_relations(K3)
#> [[1]]
#> $k
#> [1] "1^9"
#> 
#> $key
#> $key[[1]]
#> x0^9 
#> 
#> 
#> $v
#> $v[[1]]
#> list()
#> attr(,"class")
#> [1] "py_tuple"
#> 
#> 
#> attr(,"class")
#> [1] "py_od"
#> 
#> [[2]]
#> $k
#> character(0)
#> 
#> $key
#> list()
#> 
#> $v
#> list()
#> 
#> attr(,"class")
#> [1] "py_od"
```

The collector also gives the exponent vector, the depth, and the leading
exponent of an element, and decides membership in a subgroup through an
induced polycyclic generating sequence.

``` r

S3 <- SymmetricGroup(3)
PS <- perm_polycyclic_group(S3)
CH <- PS$collector
g <- Permutation(c(1, 0, 2))
exponent_vector(CH, g)
#> [1] 1 0
pc_depth(CH, g)
#> [1] 1
leading_exponent(CH, g)
#> [1] 1
ip <- induced_pcgs(CH, list(g))
ip
#> [[1]]
#> (2)(0 1)
constructive_membership_test(CH, ip, g)
#> [1] 1
PC <- PolycyclicGroup(PS$pcgs, PS$pc_series, PS$relative_order)
pc_length(PC)
#> [1] 2
K <- Collector(PS$pcgs, PS$pc_series, PS$relative_order)
exponent_vector(K, Permutation(c(1, 2, 0)))
#> [1] 0 2
```
