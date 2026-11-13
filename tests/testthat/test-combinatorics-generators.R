# Expected values are SymPy's outputs for the same calls (sympy/combinatorics/tests/test_generators.py);
# the rows come from the battery that compares this port with SymPy.

test_that("generators: cyclic", {
  expect_identical(.pyrun(local({ cyclic(6) })), "[Permutation(5), Permutation(0, 1, 2, 3, 4, 5), Permutation(0, 2, 4)(1, 3, 5), Permutation(0, 3)(1, 4)(2, 5), Permutation(0, 4, 2)(1, 5, 3), Permutation(0, 5, 4, 3, 2, 1)]")
  expect_identical(.pyrun(local({ cyclic(10) })), "[Permutation(9), Permutation(0, 1, 2, 3, 4, 5, 6, 7, 8, 9), Permutation(0, 2, 4, 6, 8)(1, 3, 5, 7, 9), Permutation(0, 3, 6, 9, 2, 5, 8, 1, 4, 7), Permutation(0, 4, 8, 2, 6)(1, 5, 9, 3, 7), Permutation(0, 5)(1, 6)(2, 7)(3, 8)(4, 9), Permutation(0, 6, 2, 8, 4)(1, 7, 3, 9, 5), Permutation(0, 7, 4, 1, 8, 5, 2, 9, 6, 3), Permutation(0, 8, 6, 4, 2)(1, 9, 7, 5, 3), Permutation(0, 9, 8, 7, 6, 5, 4, 3, 2, 1)]")
  expect_identical(.pyrun(local({ cyclic(1) })), "[Permutation(0)]")
})

test_that("generators: alt", {
  expect_identical(.pyrun(local({ alternating(4) })), "[Permutation(3), Permutation(1, 2, 3), Permutation(1, 3, 2), Permutation(0, 1)(2, 3), Permutation(3)(0, 1, 2), Permutation(0, 1, 3), Permutation(3)(0, 2, 1), Permutation(0, 2, 3), Permutation(0, 2)(1, 3), Permutation(0, 3, 1), Permutation(0, 3, 2), Permutation(0, 3)(1, 2)]")
  expect_identical(.pyrun(local({ alternating(3) })), "[Permutation(2), Permutation(0, 1, 2), Permutation(0, 2, 1)]")
  expect_identical(.pyrun(local({ length(alternating(5)) })), "60")
})

test_that("generators: sym", {
  expect_identical(.pyrun(local({ symmetric(3) })), "[Permutation(2), Permutation(1, 2), Permutation(2)(0, 1), Permutation(0, 1, 2), Permutation(0, 2, 1), Permutation(0, 2)]")
  expect_identical(.pyrun(local({ symmetric(4) })), "[Permutation(3), Permutation(2, 3), Permutation(3)(1, 2), Permutation(1, 2, 3), Permutation(1, 3, 2), Permutation(1, 3), Permutation(3)(0, 1), Permutation(0, 1)(2, 3), Permutation(3)(0, 1, 2), Permutation(0, 1, 2, 3), Permutation(0, 1, 3, 2), Permutation(0, 1, 3), Permutation(3)(0, 2, 1), Permutation(0, 2, 3, 1), Permutation(3)(0, 2), Permutation(0, 2, 3), Permutation(0, 2)(1, 3), Permutation(0, 2, 1, 3), Permutation(0, 3, 2, 1), Permutation(0, 3, 1), Permutation(0, 3, 2), Permutation(0, 3), Permutation(0, 3, 1, 2), Permutation(0, 3)(1, 2)]")
  expect_identical(.pyrun(local({ length(symmetric(5)) })), "120")
  expect_identical(.pyrun(local({ lapply(symmetric(3), function(p) .pyL(array_form(p))) })), "[[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]]")
})

test_that("generators: dih", {
  expect_identical(.pyrun(local({ dihedral(1) })), "[Permutation(1), Permutation(0, 1)]")
  expect_identical(.pyrun(local({ dihedral(2) })), "[Permutation(3), Permutation(0, 1)(2, 3), Permutation(0, 2)(1, 3), Permutation(0, 3)(1, 2)]")
  expect_identical(.pyrun(local({ dihedral(3) })), "[Permutation(2), Permutation(0, 2), Permutation(0, 1, 2), Permutation(1, 2), Permutation(0, 2, 1), Permutation(2)(0, 1)]")
  expect_identical(.pyrun(local({ dihedral(5) })), "[Permutation(4), Permutation(0, 4)(1, 3), Permutation(0, 1, 2, 3, 4), Permutation(1, 4)(2, 3), Permutation(0, 2, 4, 1, 3), Permutation(0, 1)(2, 4), Permutation(0, 3, 1, 4, 2), Permutation(0, 2)(3, 4), Permutation(0, 4, 3, 2, 1), Permutation(4)(0, 3)(1, 2)]")
  expect_identical(.pyrun(local({ lapply(dihedral(4), format) })), "['(3)', '(0 3)(1 2)', '(0 1 2 3)', '(1 3)', '(0 2)(1 3)', '(0 1)(2 3)', '(0 3 2 1)', '(3)(0 2)']")
})

test_that("generators: grp", {
  expect_identical(.pyrun(local({ as.integer(group_order(PermutationGroup(cyclic(5)))) })), "5")
  expect_identical(.pyrun(local({ as.integer(group_order(PermutationGroup(dihedral(4)))) })), "8")
})

