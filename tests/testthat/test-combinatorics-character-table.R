# Expected values are SymPy's outputs for the same calls (sympy/combinatorics/tests/test_character_table.py);
# the rows come from the battery that compares this port with SymPy.

test_that("character_table: generic", {
  expect_identical(.pyrun(local({ ct_tolist(character_table(get_perm_group(S2TransitiveSubgroups$S2))) })), "[[1, 1], [1, -1]]")
  expect_identical(.pyrun(local({ t <- character_table(get_perm_group(S2TransitiveSubgroups$S2)); .pytuple(ct_zeta_order(t), .pyT(ct_shape(t))) })), "(1, (2, 2))")
  expect_identical(.pyrun(local({ ct_tolist(character_table(get_perm_group(S3TransitiveSubgroups$S3))) })), "[[1, 1, 1], [2, 0, -1], [1, -1, 1]]")
  expect_identical(.pyrun(local({ t <- character_table(get_perm_group(S3TransitiveSubgroups$S3)); .pytuple(ct_zeta_order(t), .pyT(ct_shape(t))) })), "(1, (3, 3))")
  expect_identical(.pyrun(local({ ct_tolist(character_table(get_perm_group(S4TransitiveSubgroups$S4))) })), "[[1, 1, 1, 1, 1], [3, 1, -1, 0, -1], [2, 0, 2, -1, 0], [3, -1, -1, 0, 1], [1, -1, 1, 1, -1]]")
  expect_identical(.pyrun(local({ t <- character_table(get_perm_group(S4TransitiveSubgroups$S4)); .pytuple(ct_zeta_order(t), .pyT(ct_shape(t))) })), "(1, (5, 5))")
  expect_identical(.pyrun(local({ ct_tolist(character_table(get_perm_group(S5TransitiveSubgroups$S5))) })), "[[1, 1, 1, 1, 1, 1, 1], [4, 2, 0, 1, -1, 0, -1], [5, 1, 1, -1, 1, -1, 0], [6, 0, -2, 0, 0, 0, 1], [5, -1, 1, -1, -1, 1, 0], [4, -2, 0, 1, 1, 0, -1], [1, -1, 1, 1, -1, -1, 1]]")
  expect_identical(.pyrun(local({ t <- character_table(get_perm_group(S5TransitiveSubgroups$S5)); .pytuple(ct_zeta_order(t), .pyT(ct_shape(t))) })), "(1, (7, 7))")
  expect_identical(.pyrun(local({ ct_tolist(character_table(get_perm_group(S6TransitiveSubgroups$S3))) })), "[[1, 1, 1], [1, 1, -1], [2, -1, 0]]")
})

test_that("character_table: orth", {
  expect_identical(.pyrun(local({ G <- get_perm_group(S3TransitiveSubgroups$S3); t <- ct_tolist(character_table(G)); .pytuple(sstr(.add(lapply(t, function(r) .pow(r[[1L]], S_Two)))) == as.character(group_order(G)), all(vapply(t[[1L]], function(v) identical(v$key, S_One$key), TRUE))) })), "(True, True)")
  expect_identical(.pyrun(local({ G <- get_perm_group(S4TransitiveSubgroups$S4); t <- ct_tolist(character_table(G)); .pytuple(sstr(.add(lapply(t, function(r) .pow(r[[1L]], S_Two)))) == as.character(group_order(G)), all(vapply(t[[1L]], function(v) identical(v$key, S_One$key), TRUE))) })), "(True, True)")
  expect_identical(.pyrun(local({ G <- get_perm_group(S5TransitiveSubgroups$S5); t <- ct_tolist(character_table(G)); .pytuple(sstr(.add(lapply(t, function(r) .pow(r[[1L]], S_Two)))) == as.character(group_order(G)), all(vapply(t[[1L]], function(v) identical(v$key, S_One$key), TRUE))) })), "(True, True)")
})

test_that("character_table: sn", {
  expect_identical(.pyrun(local({ ct_tolist(character_table_sn(1)) })), "[[1]]")
  expect_identical(.pyrun(local({ ct_tolist(character_table_sn(2)) })), "[[1, 1], [1, -1]]")
  expect_identical(.pyrun(local({ ct_tolist(character_table_sn(3)) })), "[[1, 1, 1], [2, 0, -1], [1, -1, 1]]")
  expect_identical(.pyrun(local({ ct_tolist(character_table_sn(4)) })), "[[1, 1, 1, 1, 1], [3, 1, -1, 0, -1], [2, 0, 2, -1, 0], [3, -1, -1, 0, 1], [1, -1, 1, 1, -1]]")
  expect_identical(.pyrun(local({ ct_tolist(character_table_sn(5)) })), "[[1, 1, 1, 1, 1, 1, 1], [4, 2, 0, 1, -1, 0, -1], [5, 1, 1, -1, 1, -1, 0], [6, 0, -2, 0, 0, 0, 1], [5, -1, 1, -1, -1, 1, 0], [4, -2, 0, 1, 1, 0, -1], [1, -1, 1, 1, -1, -1, 1]]")
  expect_identical(.pyrun(local({ ct_conjugacy_class_reps(character_table_sn(1)) })), "[Permutation(0)]")
  expect_identical(.pyrun(local({ ct_as_matrix(character_table(SymmetricGroup(1))) })), "Matrix([[1]])")
  expect_identical(.pyrun(local({ ct_conjugacy_class_reps(character_table_sn(3)) })), "[Permutation(2), Permutation(2)(0, 1), Permutation(0, 1, 2)]")
  expect_identical(.pyrun(local({ ct_as_matrix(character_table(SymmetricGroup(3))) })), "Matrix([[1, 1, 1], [2, 0, -1], [1, -1, 1]])")
  expect_identical(.pyrun(local({ ct_conjugacy_class_reps(character_table_sn(4)) })), "[Permutation(3), Permutation(3)(0, 1), Permutation(0, 1)(2, 3), Permutation(3)(0, 1, 2), Permutation(0, 1, 2, 3)]")
  expect_identical(.pyrun(local({ ct_as_matrix(character_table(SymmetricGroup(4))) })), "Matrix([[1, 1, 1, 1, 1], [3, 1, -1, 0, -1], [2, 0, 2, -1, 0], [3, -1, -1, 0, 1], [1, -1, 1, 1, -1]])")
  expect_identical(.pyrun(local({ ct_conjugacy_class_reps(character_table_sn(5)) })), "[Permutation(4), Permutation(4)(0, 1), Permutation(4)(0, 1)(2, 3), Permutation(4)(0, 1, 2), Permutation(0, 1, 2)(3, 4), Permutation(4)(0, 1, 2, 3), Permutation(0, 1, 2, 3, 4)]")
  expect_identical(.pyrun(local({ ct_as_matrix(character_table(SymmetricGroup(5))) })), "Matrix([[1, 1, 1, 1, 1, 1, 1], [4, 2, 0, 1, -1, 0, -1], [5, 1, 1, -1, 1, -1, 0], [6, 0, -2, 0, 0, 0, 1], [5, -1, 1, -1, -1, 1, 0], [4, -2, 0, 1, 1, 0, -1], [1, -1, 1, 1, -1, -1, 1]])")
  expect_identical(.pyrun(local({ ct_conjugacy_class_reps(character_table_sn(6)) })), "[Permutation(5), Permutation(5)(0, 1), Permutation(5)(0, 1)(2, 3), Permutation(0, 1)(2, 3)(4, 5), Permutation(5)(0, 1, 2), Permutation(5)(0, 1, 2)(3, 4), Permutation(0, 1, 2)(3, 4, 5), Permutation(5)(0, 1, 2, 3), Permutation(0, 1, 2, 3)(4, 5), Permutation(5)(0, 1, 2, 3, 4), Permutation(0, 1, 2, 3, 4, 5)]")
  expect_identical(.pyrun(local({ character_table_sn(0) })), "ERR ValueError: n must be a positive integer")
  expect_identical(.pyrun(local({ character_table_sn(-3) })), "ERR ValueError: n must be a positive integer")
})

test_that("character_table: cyc", {
  expect_identical(.pyrun(local({ ct_tolist(character_table(CyclicGroup(2))) })), "[[1, 1], [1, -1]]")
})

test_that("character_table: dih", {
  expect_identical(.pyrun(local({ ct_tolist(character_table(DihedralGroup(3))) })), "[[1, 1, 1], [2, 0, -1], [1, -1, 1]]")
})

test_that("character_table: dixon", {
  expect_identical(.pyrun(local({ G <- PermutationGroup(Permutation(size = 2)); ct_tolist(dixon_character_table(conjugacy_classes(G))) })), "[[1]]")
  expect_identical(.pyrun(local({ G <- PermutationGroup(Permutation(size = 2)); ct_zeta_order(dixon_character_table(conjugacy_classes(G))) })), "1")
  expect_identical(.pyrun(local({ G <- SymmetricGroup(4); cc <- lapply(list(Permutation(list(3)), Permutation(list(c(0, 1), 3)), Permutation(list(c(0, 1), c(2, 3))), Permutation(list(c(0, 1, 2), 3)), Permutation(list(c(0, 1, 2, 3)))), function(c) conjugacy_class(G, c)); lapply(cc, length) })), "[1, 6, 3, 8, 6]")
  expect_identical(.pyrun(local({ G <- AlternatingGroup(4); cc <- lapply(list(Permutation(list(3)), Permutation(list(c(1, 2, 3))), Permutation(list(c(1, 3, 2))), Permutation(list(c(0, 1), c(2, 3)))), function(c) conjugacy_class(G, c)); lapply(cc, length) })), "[1, 4, 4, 3]")
  expect_identical(.pyrun(local({ G <- PermutationGroup(Permutation(list(c(0, 1, 3, 4, 2))), Permutation(list(c(1, 4, 2, 3)))); cc <- lapply(list(Permutation(list(4)), Permutation(list(c(1, 2), c(3, 4))), Permutation(list(c(1, 3, 2, 4))), Permutation(list(c(1, 4, 2, 3))), Permutation(list(c(0, 1, 3, 4, 2)))), function(c) conjugacy_class(G, c)); lapply(cc, length) })), "[1, 5, 5, 5, 4]")
  expect_identical(.pyrun(local({ cc <- list(list(Permutation(0, 1, 2), Permutation(0, 2, 1)), list(Permutation(2)), list(Permutation(0, 2), Permutation(1, 2), Permutation(list(c(0, 1), 2)))); ct_tolist(dixon_character_table(cc)) })), "[[1, 1, 1], [1, 1, -1], [-1, 2, 0]]")
  expect_identical(.pyrun(local({ cc <- list(list(Permutation(0, 1, 2), Permutation(0, 2, 1)), list(Permutation(2)), list(Permutation(0, 2), Permutation(1, 2), Permutation(list(c(0, 1), 2)))); ct_conjugacy_class_reps(dixon_character_table(cc))[[2L]] })), "(2)")
  expect_identical(.pyrun(local({ g <- Permutation(0, 1, 2, 3); cc <- list(list(g^2), list(g^3), list(g^4), list(g)); ct_conjugacy_class_reps(dixon_character_table(cc)) })), "[Permutation(0, 2)(1, 3), Permutation(0, 3, 2, 1), Permutation(3), Permutation(0, 1, 2, 3)]")
  expect_identical(.pyrun(local({ dixon_character_table(list(list(Permutation(0, 1)))) })), "ERR ValueError: No identity class found")
})

test_that("character_table: from", {
  expect_identical(.pyrun(local({ ct_tolist(character_table(PermutationGroup(Permutation(0, 1, 2), Permutation(0, 1)))) })), "[[1, 1, 1], [2, 0, -1], [1, -1, 1]]")
})

