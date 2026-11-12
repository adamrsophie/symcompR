# Expected values are SymPy's outputs for the same calls (Finitely presented groups (pc)).

skip_on_cran()

.bt <- new.env(parent = asNamespace("symcompR"))
eval(parse(text = r"---(NULL
)---"), envir = .bt)
.run <- function(code) {
  e <- .bt
  NULL
  gsub("\n", " ", get("fp_str", envir = .bt)(eval(parse(text = code), envir = e)), fixed = TRUE)
}
.setup <- function(code) invisible(eval(parse(text = code), envir = .bt))
.err <- function(code) tryCatch({ .run(code); "no error" }, error = function(e) gsub("\n", " ", conditionMessage(e), fixed = TRUE))

test_that("Finitely presented groups (pc) match SymPy", {
  .setup("pc_ok <- function(G) { PC <- perm_polycyclic_group(G); C <- PC$collector; P <- C$pc_presentation; for (i in seq_along(P$k)) { k <- P$key[[i]]; v <- P$v[[i]]; lhs <- .perm_af(integer()); for (t in seq_along(k$s)) lhs <- .fpp_mul(lhs, .fpp_pow(PC$pcgs[[k$s[t]]], k$e[t])); if (!inherits(v, \"FreeGroupElement\")) { if (!.is_identity(lhs)) return(FALSE); next }; rhs <- .perm_af(integer()); for (t in seq_along(v$s)) rhs <- .fpp_mul(rhs, .fpp_pow(PC$pcgs[[v$s[t]]], v$e[t])); if (!.fpp_eq(lhs, rhs)) return(FALSE) }; TRUE }")
  .setup("exp_ok <- function(G) { PC <- perm_polycyclic_group(G); C <- PC$collector; for (gen in G$gens) { ev <- exponent_vector(C, gen); g <- .perm_af(integer()); for (i in seq_along(ev)) if (ev[i] != 0) g <- .fpp_mul(g, .fpp_pow(PC$pcgs[[i]], ev[i])); if (!.fpp_eq(Permutation(g, size = length(gen$array)), gen)) return(FALSE) }; TRUE }")
  .setup("ipcgs_ok <- function(G) { PC <- perm_polycyclic_group(G); C <- PC$collector; ip <- induced_pcgs(C, G$gens); m <- lapply(ip, function(i) exponent_vector(C, i)); for (i in seq_along(m)) for (j in seq_len(min(i - 1L, length(m[[i]])))) if (m[[i]][j] != 0) return(FALSE); TRUE }")
  .setup("info <- function(G) { P <- perm_polycyclic_group(G); .py_tuple(list(length(P$pcgs), sort(P$relative_order), pc_is_prime_order(P), pc_length(P), .ppg_order(G))) }")
  .setup("S93 <- PermutationGroup(list(Permutation(c(1, 2, 0, 3, 4, 5, 6, 7, 8)), Permutation(c(3, 4, 5, 6, 7, 8, 0, 1, 2)))); S92 <- PermutationGroup(list(Permutation(c(1, 0, 2, 3, 4, 5, 6, 7, 8)), Permutation(c(2, 3, 0, 1, 4, 5, 6, 7, 8)), Permutation(c(4, 5, 6, 7, 0, 1, 2, 3, 8)))); S82 <- PermutationGroup(list(Permutation(c(1, 0, 2, 3, 4, 5, 6, 7)), Permutation(c(2, 3, 0, 1, 4, 5, 6, 7)), Permutation(c(4, 5, 6, 7, 0, 1, 2, 3))))")
  .setup("S255 <- PermutationGroup(list(Permutation(c(1, 2, 3, 4, 0, 5:24)), Permutation(c(5:24, 0:4))))")
  .setup("S202 <- PermutationGroup(list(Permutation(c(1, 0, 2:19)), Permutation(c(2, 3, 0, 1, 4:19)), Permutation(c(4, 5, 6, 7, 0, 1, 2, 3, 8:19)), Permutation(c(8:15, 0:7, 16:19)), Permutation(c(0:15, 17, 16, 18, 19)), Permutation(c(0:15, 18, 19, 16, 17))))")
  expect_equal(.run("vapply(list(SymmetricGroup(3), SymmetricGroup(4), S93, S92, S82, DihedralGroup(10)), pc_ok, NA)"), "[True, True, True, True, True, True]")
  expect_equal(.run("pc_ok(S255)"), "True")
  expect_equal(.run("vapply(list(SymmetricGroup(3), SymmetricGroup(4), S93, S92, S82), exp_ok, NA)"), "[True, True, True, True, True]")
  expect_equal(.run("vapply(list(S93, AlternatingGroup(4), DihedralGroup(4), DihedralGroup(10), DihedralGroup(9), SymmetricGroup(3), SymmetricGroup(4)), ipcgs_ok, NA)"), "[True, True, True, True, True, True, True]")
  expect_equal(.run("ipcgs_ok(S202)"), "True")
  expect_equal(.run("lapply(list(SymmetricGroup(3), SymmetricGroup(4), S93, S92, S82, DihedralGroup(10), AlternatingGroup(4), DihedralGroup(9)), info)"), "[(2, [2, 3], True, 2, 6), (4, [2, 2, 2, 3], True, 4, 24), (4, [3, 3, 3, 3], True, 4, 81), (7, [2, 2, 2, 2, 2, 2, 2], True, 7, 128), (7, [2, 2, 2, 2, 2, 2, 2], True, 7, 128), (3, [2, 2, 5], True, 3, 20), (3, [2, 2, 3], True, 3, 12), (2, [2, 9], False, 2, 18)]")
  expect_equal(.err("perm_polycyclic_group(SymmetricGroup(5))"), "The group must be solvable")
  .setup("H <- PermutationGroup(list(Permutation(c(1, 0, 2)))); PH <- perm_polycyclic_group(H); CH <- PH$collector; ipH <- induced_pcgs(CH, H$gens)")
  expect_equal(.run("length(ipH) > 0"), "True")
  expect_equal(.run("ipH"), "[Permutation(2)(0, 1)]")
  expect_equal(.run("I(constructive_membership_test(CH, ipH, H$gens[[1]]))"), "[1]")
  expect_equal(.run("I(constructive_membership_test(CH, ipH, group_identity(H)))"), "[0]")
  expect_equal(.run("constructive_membership_test(CH, ipH, Permutation(c(0, 2, 1)))"), "False")
  expect_equal(.run("structure(lapply(pc_presentation(CH), function(p) paste0(fp_str(p[[1]]), \": \", fp_str(p[[2]]))), class = \"py_dictlike\")"), "{x0^2: ()}")
  expect_equal(.run("PH$pcgs"), "[Permutation(2)(0, 1)]")
  expect_equal(.run("I(PH$relative_order)"), "[2]")
  expect_equal(.run("I(exponent_vector(CH, Permutation(c(1, 0, 2))))"), "[1]")
  expect_equal(.run("pc_depth(CH, Permutation(c(1, 0, 2)))"), "1")
  expect_equal(.run("pc_depth(CH, Permutation(c(0, 1, 2)))"), "2")
  expect_equal(.run("leading_exponent(CH, Permutation(c(1, 0, 2)))"), "1")
  expect_equal(.run("leading_exponent(CH, Permutation(c(0, 1, 2)))"), "None")
  .setup("C3 <- CyclicGroup(9); P3 <- perm_polycyclic_group(C3); K3 <- P3$collector")
  expect_equal(.run(".py_tuple(list(P3$pcgs, I(P3$relative_order)))"), "([Permutation(0, 1, 2, 3, 4, 5, 6, 7, 8)], [9])")
  expect_equal(.run("structure(lapply(pc_presentation(K3), function(p) paste0(fp_str(p[[1]]), \": \", fp_str(p[[2]]))), class = \"py_dictlike\")"), "{x0^9: ()}")
  expect_equal(.run("lapply(0:8, function(k) I(exponent_vector(K3, C3$gens[[1]]^k)))"), "[[0], [1], [2], [3], [4], [5], [6], [7], [8]]")
  expect_equal(.run("{.r <- minimal_uncollected_subword(K3, K3$free_group$generators[[1]]^-2); .py_tuple(lapply(.r, function(p) .py_tuple(list(Symbol(K3$free_group$symbols[p[1]]), p[2]))))}"), "((x0, -2),)")
  expect_equal(.run("collected_word(K3, K3$free_group$generators[[1]]^13)"), "x0^4")
  expect_equal(.run("{.r <- minimal_uncollected_subword(K3, K3$free_group$generators[[1]]^3); if (is.null(.r)) NULL else .py_tuple(lapply(.r, function(p) .py_tuple(list(Symbol(K3$free_group$symbols[p[1]]), p[2]))))}"), "None")
  expect_equal(.run("vapply(0:8, function(k) pc_depth(K3, C3$gens[[1]]^k), 0L)"), "[2, 1, 1, 1, 1, 1, 1, 1, 1]")
  expect_equal(.run("induced_pcgs(K3, list(C3$gens[[1]]^3))"), "[Permutation(0, 3, 6)(1, 4, 7)(2, 5, 8)]")
})
