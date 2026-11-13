# Permutations and permutation groups

``` r

library(symcompR)
```

This vignette covers `sympy.combinatorics.permutations`, `perm_groups`,
`named_groups`, `group_constructs`, `polyhedron`, and the Schreier-Sims
utilities. Free and finitely presented groups, homomorphisms, and
polycyclic groups are in
[`vignette("fp-groups", package = "symcompR")`](https://adamrsophie.github.io/symcompR/articles/fp-groups.md).

## Permutations

Permutations act on the points 0, 1, \ldots, n - 1, as in SymPy. Create
one from its array form (the images of 0, 1, \ldots), from a list of
cycles, or from the points of a single cycle. Permutations print in
cycle notation.

``` r

p <- Permutation(c(1, 2, 0))
p
#> (0 1 2)
array_form(p)
#> [1] 1 2 0
cyclic_form(p)
#> [[1]]
#> [1] 0 1 2
q <- Permutation(list(c(0, 1)), size = 3)
q
#> (2)(0 1)
Permutation(1, 2, size = 4)
#> (3)(1 2)
```

The product `p * q` applies `p` first and then `q`, SymPy’s convention.

``` r

p * q
#> (1 2)
p^-1
#> (0 2 1)
p^3
#> (2)
```

Properties of a permutation:

``` r

perm_order(p)
#> [1] 3
perm_parity(p)
#> [1] 0
perm_is_even(p)
#> [1] TRUE
perm_signature(p)
#> [1] 1
perm_cycle_structure(p)
#> 3 
#> 1
perm_inversions(Permutation(c(2, 0, 1, 3)))
#> [1] 2
perm_transpositions(Permutation(c(0, 2, 1)))
#> [[1]]
#> [1] 1 2
```

Permutations of a given size are ranked in lexicographic order:

``` r

perm_rank(p)
#> [1] 3
unrank_lex(4, 10)
#> (0 1 3 2)
```

[`Cycle()`](https://adamrsophie.github.io/symcompR/reference/Cycle.md)
is a lightweight cycle object, convenient for building permutations by
composing cycles.

``` r

Cycle(1, 2)
#> (1 2)
cycle_mul(Cycle(1, 2), 2, 3)
#> (1 3 2)
```

More operations on single permutations.
[`perm_call()`](https://adamrsophie.github.io/symcompR/reference/Permutation.md)
applies a permutation to points,
[`apply_perm()`](https://adamrsophie.github.io/symcompR/reference/Permutation.md)
to one point, and
[`perm_list()`](https://adamrsophie.github.io/symcompR/reference/Permutation.md)
gives the array form, optionally resized.

``` r

r <- Permutation(c(2, 5, 1, 6, 3, 0, 4))
perm_size(r)
#> [1] 7
perm_call(r, 0:6)
#> [1] 2 5 1 6 3 0 4
apply_perm(r, 0)
#> [1] 2
perm_list(Permutation(1, 2), 5)
#> [1] 0 2 1 3 4
perm_support(Permutation(1, 3, size = 5))
#> [1] 1 3
perm_length(Permutation(1, 3, size = 5))
#> [1] 2
perm_cardinality(r)
#> [1] 5040
perm_atoms(r)
#> [1] 0 1 2 3 4 5 6
perm_max(Permutation(list(6, 5, 0:4)))
#> [1] 4
perm_min(Permutation(list(6, 5, 0:4)))
#> [1] 0
perm_cycles(r)
#> [1] 2
full_cyclic_form(Permutation(list(1, c(0, 3, 5, 6, 2, 4))))
#> [[1]]
#> [1] 0 3 5 6 2 4
#> 
#> [[2]]
#> [1] 1
perm_sign(r)
#> [1] -1
perm_is_odd(Permutation(c(0, 4, 1, 3, 2)))
#> [1] FALSE
perm_is_identity(Permutation(c(3, 2, 1, 0))^2)
#> [1] TRUE
perm_is_empty(Permutation(integer()))
#> [1] TRUE
perm_is_singleton(Permutation(0))
#> [1] TRUE
perm_resize(Permutation(0, 1, 2), 3)
#> (0 1 2)
```

Products, inverses, and commutators.
[`perm_rmul()`](https://adamrsophie.github.io/symcompR/reference/Permutation.md)
multiplies from right to left, and `perm_mul_inv(p, q)` is p^{-1} q.

``` r

perm_inverse(r)
#> (0 5 1 2)(3 4 6)
perm_commutator(r, Permutation(list(1, c(0, 3, 5, 6, 2, 4))))
#> (0 1 3 4 6 5 2)
perm_commutes_with(Permutation(0, 1, size = 4), Permutation(2, 3))
#> [1] TRUE
perm_prod(list(Permutation(0, 1), Permutation(1, 2)))
#> (0 2 1)
perm_rmul(Permutation(c(1, 2, 3, 0)), Permutation(c(3, 2, 0, 1)),
          Permutation(c(1, 3, 0, 2)))
#> (0 3 1 2)
perm_mul_inv(Permutation(c(4, 1, 3, 0, 2)), Permutation(c(1, 0, 3, 4, 2)))
#> (4)(0 1 3)
```

Statistics and distances used in sorting and ranking:

``` r

s <- Permutation(c(1, 5, 2, 0, 3, 6, 4))
perm_ascents(s)
#> [1] 0 3 4
perm_descents(s)
#> [1] 1 2 5
perm_runs(s)
#> [[1]]
#> [1] 1 5
#> 
#> [[2]]
#> [1] 2
#> 
#> [[3]]
#> [1] 0 3 6
#> 
#> [[4]]
#> [1] 4
perm_index(Permutation(c(3, 2, 1, 0)))
#> [1] 3
perm_inversion_vector(r)
#> [1] 2 4 1 3 1 0
perm_from_inversion_vector(c(2, 4, 1, 3, 1, 0))
#> (0 2 1 5)(3 6 4)
perm_adjacency_matrix(Permutation(c(1, 2, 0)))
#> Matrix(3, 3)
#> [0, 0, 0]
#> [0, 0, 1]
#> [1, 0, 0]
perm_precedence_matrix(Permutation(c(1, 2, 0)))
#> Matrix(3, 3)
#> [0, 0, 0]
#> [1, 0, 1]
#> [1, 0, 0]
p0 <- Permutation(0:3); p1 <- Permutation(3:0)
perm_adjacency_distance(p0, p1)
#> [1] 3
perm_positional_distance(p0, p1)
#> [1] 8
perm_precedence_distance(p0, p1)
#> [1] 6
```

Ranking and unranking in the lexicographic, Trotter-Johnson, and
Myrvold-Ruskey (nonlexicographic) orders:

``` r

perm_next_lex(Permutation(c(0, 2, 1, 3)))
#> (1 2 3)
perm_rank_trotterjohnson(Permutation(c(0, 2, 1, 3)))
#> [1] 7
perm_unrank_trotterjohnson(4, 5)
#> (1 3)
perm_next_trotterjohnson(Permutation(0:3))
#> (2 3)
perm_rank_nonlex(Permutation(c(0, 2, 1, 3)))
#> [1] 19
perm_unrank_nonlex(7, 41)
#> (6)(0 1 2 3 4)
perm_next_nonlex(Permutation(c(1, 0)))
#> (1)
```

Permutations from other data.
[`perm_from_sequence()`](https://adamrsophie.github.io/symcompR/reference/Permutation.md)
gives the permutation that sorts a sequence, and
[`perm_josephus()`](https://adamrsophie.github.io/symcompR/reference/Permutation.md)
gives the order in which the Josephus problem eliminates people.

``` r

perm_from_sequence("SymPy")
#> (4)(0 1 3)
perm_josephus(3, 6, 1)
#> (0 2 3 1 5)
```

The random permutation generators use a seeded generator that follows
SymPy’s:

``` r

sympy_seed(3)
perm_random(6)
#> (1 2 3 5)
```

[`cycle_as_list()`](https://adamrsophie.github.io/symcompR/reference/Cycle.md)
and
[`cycle_size()`](https://adamrsophie.github.io/symcompR/reference/Cycle.md)
work on `Cycle` objects.

``` r

cycle_as_list(Cycle(1, 2), 4)
#> [1] 0 2 1 3
cycle_size(Cycle(1, 2))
#> [1] 3
```

## Permutation groups

[`PermutationGroup()`](https://adamrsophie.github.io/symcompR/reference/PermutationGroup.md)
is generated by a set of permutations. Group algorithms use a base and
strong generating set computed with the Schreier-Sims algorithm, so they
scale to large groups such as the Rubik’s cube group.

``` r

G <- PermutationGroup(Permutation(0, 1), Permutation(2, 3))
orbits(G)
#> [[1]]
#> [1] 0 1
#> 
#> [[2]]
#> [1] 2 3
group_order(G)
#> [1] 4
group_order(RubikGroup(2))
#> [1] 3674160
```

The 2 \times 2 \times 2 cube has 3,674,160 positions.

### Named groups

``` r

D4 <- DihedralGroup(4)
D4
#> PermutationGroup([
#>     (0 1 2 3),
#>     (0 3)(1 2)])
group_order(D4)
#> [1] 8
S4 <- SymmetricGroup(4)
group_order(S4)
#> [1] 24
group_order(AlternatingGroup(5))
#> [1] 60
is_cyclic(CyclicGroup(6))
#> [1] TRUE
group_order(AbelianGroup(2, 3, 4))
#> [1] 24
group_elements(CyclicGroup(3))
#> [[1]]
#> (2)
#> 
#> [[2]]
#> (0 1 2)
#> 
#> [[3]]
#> (0 2 1)
```

### Structure

``` r

is_abelian(D4)
#> [1] FALSE
is_solvable(S4)
#> [1] TRUE
is_solvable(AlternatingGroup(5))
#> [1] FALSE
group_base(S4)
#> [1] 0 1 2
group_order(group_center(D4))
#> [1] 2
group_order(derived_subgroup(D4))
#> [1] 2
vapply(derived_series(S4), group_order, 0)
#> [1] 24 12  4  1
vapply(composition_series(S4), group_order, 0)
#> [1] 24 12  4  2  1
abelian_invariants(S4)
#> [1] 2
abelian_invariants(AbelianGroup(2, 3, 4))
#> [1] 2 3 4
```

The derived series S_4 \supset A_4 \supset V_4 \supset 1 shows that S_4
is solvable, while A_5 is the smallest group that is not.

Subgroups, orbits, stabilizers, and Sylow subgroups:

``` r

orbit(D4, 0)
#> [1] 0 1 2 3
group_order(stabilizer(D4, 0))
#> [1] 2
group_contains(S4, Permutation(0, 1, 2, 3))
#> [1] TRUE
is_subgroup(D4, S4)
#> [1] TRUE
is_transitive(D4)
#> [1] TRUE
is_normal(AlternatingGroup(4), S4)
#> [1] TRUE
group_order(sylow_subgroup(S4, 2))
#> [1] 8
group_order(sylow_subgroup(S4, 3))
#> [1] 3
length(conjugacy_classes(D4))
#> [1] 5
```

### Products and more named groups

[`DirectProduct()`](https://adamrsophie.github.io/symcompR/reference/SymmetricGroup.md)
builds the direct product of permutation groups.
[`SymmetricPermutationGroup()`](https://adamrsophie.github.io/symcompR/reference/Coset.md)
is the symmetric group as a lazy object that stores only its degree.
[`rubik()`](https://adamrsophie.github.io/symcompR/reference/SymmetricGroup.md)
and
[`rubik_cube_generators()`](https://adamrsophie.github.io/symcompR/reference/SymmetricGroup.md)
give the generators of the Rubik’s cube groups.

``` r

group_order(DirectProduct(CyclicGroup(2), CyclicGroup(3)))
#> [1] 6
orbits(DirectProduct(CyclicGroup(2), CyclicGroup(2)))
#> [[1]]
#> [1] 0 1
#> 
#> [[2]]
#> [1] 2 3
SymmetricPermutationGroup(5)$degree
#> [1] 5
group_order(PermutationGroup(rubik_cube_generators()))
#> Big Integer ('bigz') :
#> [1] 43252003274489856000
length(rubik(2))
#> [1] 3
group_degree(SymmetricGroup(5))
#> [1] 5
group_generators(CyclicGroup(4))
#> [[1]]
#> (0 1 2 3)
group_identity(CyclicGroup(4))
#> (3)
group_equals(PermutationGroup(Permutation(0, 1, 3), Permutation(0, 2, 3)),
             PermutationGroup(Permutation(0, 2, 3), Permutation(0, 1, 3)))
#> [1] TRUE
```

### Recognizing groups

``` r

is_dihedral(DihedralGroup(6))
#> [1] TRUE
is_symmetric_group(SymmetricGroup(5))
#> [1] TRUE
is_alternating_group(AlternatingGroup(5))
#> [1] TRUE
is_alt_sym(DihedralGroup(9))
#> [1] FALSE
is_nilpotent(SymmetricGroup(3))
#> [1] FALSE
is_nilpotent(DihedralGroup(4))
#> [1] TRUE
is_perfect_group(AlternatingGroup(5))
#> [1] TRUE
is_polycyclic(SymmetricGroup(4))
#> [1] TRUE
is_elementary(sylow_subgroup(AlternatingGroup(4), 2), 2)
#> [1] TRUE
is_trivial(PermutationGroup(Permutation(3)))
#> [1] TRUE
is_primitive(SymmetricGroup(5))
#> [1] TRUE
is_primitive(DihedralGroup(6))
#> [1] FALSE
transitivity_degree(SymmetricGroup(4))
#> [1] 4
transitivity_degree(DihedralGroup(6))
#> [1] 1
max_div(SymmetricGroup(10))
#> [1] 5
```

### Subgroups and quotients

``` r

S4 <- SymmetricGroup(4); A4 <- AlternatingGroup(4); D4 <- DihedralGroup(4)
group_index(S4, A4)
#> [1] 2
group_order(subgroup(S4, list(Permutation(0, 1, 2, 3))))
#> [1] 4
group_order(centralizer(S4, Permutation(0, 1, 2, 3)))
#> [1] 4
group_order(normal_closure(S4, Permutation(0, 1, 2, size = 4)))
#> [1] 12
group_order(group_commutator(S4, A4, S4))
#> [1] 12
vapply(lower_central_series(S4), group_order, 0)
#> [1] 24 12
group_order(quotient_group(S4, PermutationGroup(
  Permutation(list(c(0, 1), c(2, 3))), Permutation(list(c(0, 2), c(1, 3))))))
#> [1] 6
group_order(pointwise_stabilizer(S4, c(0, 1)))
#> [1] 2
conjugacy_class(D4, Permutation(0, 1, 2, 3))
#> [[1]]
#> (0 1 2 3)
#> 
#> [[2]]
#> (0 3 2 1)
group_order(subgroup_search(SymmetricGroup(5), perm_is_even))
#> [1] 60
```

[`minimal_block()`](https://adamrsophie.github.io/symcompR/reference/PermutationGroup.md)
finds the smallest block containing given points, and
[`minimal_blocks()`](https://adamrsophie.github.io/symcompR/reference/PermutationGroup.md)
the minimal block systems of a transitive group.

``` r

minimal_block(DihedralGroup(6), c(0, 3))
#> [1] 0 1 2 0 1 2
minimal_blocks(DihedralGroup(6))
#> [[1]]
#> [1] 0 1 0 1 0 1
#> 
#> [[2]]
#> [1] 0 1 2 0 1 2
```

### Cosets

[`coset_factor()`](https://adamrsophie.github.io/symcompR/reference/PermutationGroup.md)
writes an element as a product of coset representatives, and
[`coset_rank()`](https://adamrsophie.github.io/symcompR/reference/PermutationGroup.md)
and
[`coset_unrank()`](https://adamrsophie.github.io/symcompR/reference/PermutationGroup.md)
number the elements of a group.
[`Coset()`](https://adamrsophie.github.io/symcompR/reference/Coset.md)
and the `*` operator build cosets.

``` r

coset_rank(SymmetricGroup(4), Permutation(c(1, 0, 3, 2)))
#> [1] 21
coset_unrank(SymmetricGroup(4), 7)
#> (0 3 1 2)
coset_factor(SymmetricGroup(4), Permutation(c(1, 0, 3, 2)))
#> [[1]]
#> (0 1 2 3)
#> 
#> [[2]]
#> (1 3 2)
#> 
#> [[3]]
#> (2 3)
generator_product(SymmetricGroup(4), Permutation(0, 3))
#> [[1]]
#> (1 3 2)
#> 
#> [[2]]
#> (0 1 2 3)
coset_transversal(S4, A4)
#> [[1]]
#> (3)
#> 
#> [[2]]
#> (2 3)
coset_table(S4, A4)
#> [[1]]
#> [1] 1 1 1 1
#> 
#> [[2]]
#> [1] 0 0 0 0
cc <- Coset(Permutation(1, 2), PermutationGroup(Permutation(0, 1, size = 3)))
is_right_coset(cc)
#> [1] TRUE
is_left_coset(cc)
#> [1] FALSE
coset_as_list(cc)
#> [[1]]
#> (1 2)
#> 
#> [[2]]
#> (0 2 1)
```

### Bases, strong generating sets, and Schreier-Sims

The group algorithms rely on a base and strong generating set (BSGS)
computed by the Schreier-Sims algorithm. These functions expose its
parts.

``` r

S <- SymmetricGroup(4)
schreier_sims(S)
group_base(S)
#> [1] 0 1 2
strong_gens(S)
#> [[1]]
#> (0 1 2 3)
#> 
#> [[2]]
#> (3)(0 1)
#> 
#> [[3]]
#> (1 3 2)
#> 
#> [[4]]
#> (2 3)
basic_orbits(S)
#> [[1]]
#> [1] 0 1 2 3
#> 
#> [[2]]
#> [1] 1 2 3
#> 
#> [[3]]
#> [1] 2 3
lapply(basic_stabilizers(S), group_order)
#> [[1]]
#> [1] 24
#> 
#> [[2]]
#> [1] 6
#> 
#> [[3]]
#> [1] 2
basic_transversals(S)[[2]]
#> $`1`
#> (3)
#> 
#> $`3`
#> (1 3 2)
#> 
#> $`2`
#> (1 2 3)
schreier_vector(S, 1)
#> [1]  1 -1  0  0
orbit_rep(SymmetricGroup(5), 0, 4)
#> (0 4 3 2 1)
orbit_transversal(S, 2)
#> [[1]]
#> (3)
#> 
#> [[2]]
#> (0 1 2 3)
#> 
#> [[3]]
#> (0 2)(1 3)
#> 
#> [[4]]
#> (0 3 2 1)
bs <- schreier_sims_incremental(SymmetricGroup(5), base = c(0, 1, 2))
bs[[1]]
#> [1] 0 1 2 3
length(bs[[2]])
#> [1] 5
sw <- baseswap(S, group_base(S), strong_gens(S), 1, randomized = FALSE)
sw[[1]]
#> [1] 0 2 1
```

Randomized versions use the product replacement algorithm:

``` r

sympy_seed(1)
length(schreier_sims_random(SymmetricGroup(5), consec_succ = 5)[[1]])
#> [1] 4
perm_is_even(random_pr(AlternatingGroup(5)))
#> [1] TRUE
apply_perm(random_stab(SymmetricGroup(5), 2), 2)
#> [1] 2
group_contains(D4, group_random(D4))
#> [1] TRUE
```

### Listing the elements

[`group_generate()`](https://adamrsophie.github.io/symcompR/reference/PermutationGroup.md)
lists the elements by coset enumeration or with the Dimino algorithm
([`generate_dimino()`](https://adamrsophie.github.io/symcompR/reference/PermutationGroup.md)),
and
[`generate_schreier_sims()`](https://adamrsophie.github.io/symcompR/reference/PermutationGroup.md)
with the Schreier-Sims transversals.

``` r

length(group_generate(DihedralGroup(5)))
#> [1] 10
generate_dimino(CyclicGroup(3))
#> [[1]]
#> (2)
#> 
#> [[2]]
#> (0 1 2)
#> 
#> [[3]]
#> (0 2 1)
generate_schreier_sims(CyclicGroup(3))
#> [[1]]
#> (2)
#> 
#> [[2]]
#> (0 1 2)
#> 
#> [[3]]
#> (0 2 1)
```

### Invariant theory

[`molien()`](https://adamrsophie.github.io/symcompR/reference/PermutationGroup.md)
is the Molien series of the action of a permutation group on a
polynomial ring.

``` r

molien(SymmetricGroup(3))
#> 1/(3*(1 - t^3)) + 1/(2*(1 - t)*(1 - t^2)) + 1/(6*(1 - t)^3)
```

## Polyhedra

The symmetry groups of the Platonic solids are available as `Polyhedron`
objects, whose `pgroup` permutes the corners.

``` r

C <- cube()
poly_size(C)
#> [1] 8
group_order(poly_pgroup(C))
#> [1] 24
group_order(poly_pgroup(tetrahedron()))
#> [1] 12
group_order(poly_pgroup(icosahedron()))
#> [1] 60
```

[`tetrahedron()`](https://adamrsophie.github.io/symcompR/reference/Polyhedron.md),
[`cube()`](https://adamrsophie.github.io/symcompR/reference/Polyhedron.md),
[`octahedron()`](https://adamrsophie.github.io/symcompR/reference/Polyhedron.md),
[`dodecahedron()`](https://adamrsophie.github.io/symcompR/reference/Polyhedron.md),
and
[`icosahedron()`](https://adamrsophie.github.io/symcompR/reference/Polyhedron.md)
are the five Platonic solids. Their corners, faces, and edges are
available, and
[`poly_rotate()`](https://adamrsophie.github.io/symcompR/reference/Polyhedron.md)
applies a permutation from the group to the corners in place.

``` r

vapply(list(tetrahedron(), cube(), octahedron(), dodecahedron(),
            icosahedron()), poly_size, 0L)
#> [1]  4  8  6 20 12
polyhedron_faces("tetrahedron")
#> [[1]]
#> [[1]]
#> [1] 0
#> 
#> [[2]]
#> [1] 1
#> 
#> [[3]]
#> [1] 2
#> 
#> attr(,"class")
#> [1] "py_tuple"
#> 
#> [[2]]
#> [[1]]
#> [1] 0
#> 
#> [[2]]
#> [1] 1
#> 
#> [[3]]
#> [1] 3
#> 
#> attr(,"class")
#> [1] "py_tuple"
#> 
#> [[3]]
#> [[1]]
#> [1] 0
#> 
#> [[2]]
#> [1] 2
#> 
#> [[3]]
#> [1] 3
#> 
#> attr(,"class")
#> [1] "py_tuple"
#> 
#> [[4]]
#> [[1]]
#> [1] 1
#> 
#> [[2]]
#> [1] 2
#> 
#> [[3]]
#> [1] 3
#> 
#> attr(,"class")
#> [1] "py_tuple"
C <- cube()
poly_corners(C)
#> [[1]]
#> [1] 0
#> 
#> [[2]]
#> [1] 1
#> 
#> [[3]]
#> [1] 2
#> 
#> [[4]]
#> [1] 3
#> 
#> [[5]]
#> [1] 4
#> 
#> [[6]]
#> [1] 5
#> 
#> [[7]]
#> [1] 6
#> 
#> [[8]]
#> [1] 7
#> 
#> attr(,"class")
#> [1] "py_tuple"
poly_faces(tetrahedron())
#> [[1]]
#> [[1]]
#> [1] 0
#> 
#> [[2]]
#> [1] 1
#> 
#> [[3]]
#> [1] 2
#> 
#> attr(,"class")
#> [1] "py_tuple"
#> 
#> [[2]]
#> [[1]]
#> [1] 0
#> 
#> [[2]]
#> [1] 1
#> 
#> [[3]]
#> [1] 3
#> 
#> attr(,"class")
#> [1] "py_tuple"
#> 
#> [[3]]
#> [[1]]
#> [1] 0
#> 
#> [[2]]
#> [1] 2
#> 
#> [[3]]
#> [1] 3
#> 
#> attr(,"class")
#> [1] "py_tuple"
#> 
#> [[4]]
#> [[1]]
#> [1] 1
#> 
#> [[2]]
#> [1] 2
#> 
#> [[3]]
#> [1] 3
#> 
#> attr(,"class")
#> [1] "py_tuple"
#> 
#> attr(,"class")
#> [1] "py_set"
poly_edges(tetrahedron())
#> [[1]]
#> [[1]]
#> [1] 0
#> 
#> [[2]]
#> [1] 1
#> 
#> attr(,"class")
#> [1] "py_tuple"
#> 
#> [[2]]
#> [[1]]
#> [1] 0
#> 
#> [[2]]
#> [1] 2
#> 
#> attr(,"class")
#> [1] "py_tuple"
#> 
#> [[3]]
#> [[1]]
#> [1] 0
#> 
#> [[2]]
#> [1] 3
#> 
#> attr(,"class")
#> [1] "py_tuple"
#> 
#> [[4]]
#> [[1]]
#> [1] 1
#> 
#> [[2]]
#> [1] 2
#> 
#> attr(,"class")
#> [1] "py_tuple"
#> 
#> [[5]]
#> [[1]]
#> [1] 1
#> 
#> [[2]]
#> [1] 3
#> 
#> attr(,"class")
#> [1] "py_tuple"
#> 
#> [[6]]
#> [[1]]
#> [1] 2
#> 
#> [[2]]
#> [1] 3
#> 
#> attr(,"class")
#> [1] "py_tuple"
#> 
#> attr(,"class")
#> [1] "py_set"
poly_rotate(C, 1)
poly_corners(C)
#> [[1]]
#> [1] 4
#> 
#> [[2]]
#> [1] 0
#> 
#> [[3]]
#> [1] 3
#> 
#> [[4]]
#> [1] 7
#> 
#> [[5]]
#> [1] 5
#> 
#> [[6]]
#> [1] 1
#> 
#> [[7]]
#> [1] 2
#> 
#> [[8]]
#> [1] 6
#> 
#> attr(,"class")
#> [1] "py_tuple"
poly_array_form(C)
#> [1] 4 0 3 7 5 1 2 6
poly_cyclic_form(C)
#> [[1]]
#> [1] 0 4 5 1
#> 
#> [[2]]
#> [1] 2 3 7 6
poly_reset(C)
poly_corners(C)
#> [[1]]
#> [1] 0
#> 
#> [[2]]
#> [1] 1
#> 
#> [[3]]
#> [1] 2
#> 
#> [[4]]
#> [1] 3
#> 
#> [[5]]
#> [1] 4
#> 
#> [[6]]
#> [1] 5
#> 
#> [[7]]
#> [1] 6
#> 
#> [[8]]
#> [1] 7
#> 
#> attr(,"class")
#> [1] "py_tuple"
D <- poly_copy(C)
poly_rotate(D, poly_pgroup(D)$gens[[1]])
poly_corners(D)
#> [[1]]
#> [1] 1
#> 
#> [[2]]
#> [1] 2
#> 
#> [[3]]
#> [1] 3
#> 
#> [[4]]
#> [1] 0
#> 
#> [[5]]
#> [1] 5
#> 
#> [[6]]
#> [1] 6
#> 
#> [[7]]
#> [1] 7
#> 
#> [[8]]
#> [1] 4
#> 
#> attr(,"class")
#> [1] "py_tuple"
make_perm(poly_pgroup(C), 2, seed = 7)
#> (0 5)(1 4)(2 7)(3 6)
```

[`Polyhedron()`](https://adamrsophie.github.io/symcompR/reference/Polyhedron.md)
builds a polyhedron from its corners, faces, and a list of permutations
generating its symmetry group. Here a triangle is rotated by its only
generator.

``` r

fmt <- function(s) vapply(s, function(t) paste(unlist(t), collapse = " "), "")
tri <- Polyhedron(list("a", "b", "c"), list(c(1, 2, 0)),
                  list(Permutation(c(1, 2, 0))))
fmt(poly_faces(tri))
#> [1] "0 1 2"
fmt(poly_edges(tri))
#> [1] "0 1" "0 2" "1 2"
poly_rotate(tri, 0)
unlist(poly_corners(tri))
#> [1] "b" "c" "a"
```
