skip_unless_full_tests()

# Expected values are SymPy's outputs for the same calls (sympy/combinatorics/tests/test_partitions.py);
# the rows come from the battery that compares this port with SymPy.

test_that("partitions: ctor matches SymPy", {
  expect_identical(.pyrun(local({ Partition(c(1, 1, 2)) })), "ERR ValueError: Partition contained duplicate elements.")
  expect_identical(.pyrun(local({ Partition(c(1, 2, 3), c(2, 3, 4)) })), "ERR ValueError: Partition contained duplicate elements.")
  expect_identical(.pyrun(local({ Partition(Integer(1), Integer(2), Integer(3)) })), "ERR ValueError: Each argument to Partition should be a list, set, or a FiniteSet")
  expect_identical(.pyrun(local({ Partition(c(1, 2, 3), c(4, 5)) == Partition(c(4, 5), c(1, 2, 3)) })), "True")
  expect_identical(.pyrun(local({ Partition(c(1, 2, 3), c(4, 5)) == Partition(c(1, 2, 3), c(4, 5)) })), "True")
  expect_identical(.pyrun(local({ Partition(FiniteSet(1, 2, 3), FiniteSet(4, 5)) == Partition(c(1, 2, 3), c(4, 5)) })), "True")
  expect_identical(.pyrun(local({ a <- FiniteSet(1, 2, 3); b <- FiniteSet(4, 5); Partition(FiniteSet(a, b)) == Partition(FiniteSet(a, b)) })), "True")
  expect_identical(.pyrun(local({ a <- FiniteSet(1, 2, 3); b <- FiniteSet(4, 5); Partition(FiniteSet(a, b)) != Partition(a, b) })), "True")
})

test_that("partitions: str matches SymPy", {
  expect_identical(.pyrun(local({ Partition(c(1, 2, 3), c(4, 5)) })), "Partition({4, 5}, {1, 2, 3})")
  expect_identical(.pyrun(local({ Partition(c(1, 2, 3), 4) })), "Partition({4}, {1, 2, 3})")
  expect_identical(.pyrun(local({ Partition(3, c(1, 2)) })), "Partition({3}, {1, 2})")
  expect_identical(.pyrun(local({ Partition(c(4, 5), c(1, 2, 3)) })), "Partition({4, 5}, {1, 2, 3})")
  expect_identical(.pyrun(local({ Partition(FiniteSet(FiniteSet(1, 2, 3), FiniteSet(4, 5))) })), "Partition(FiniteSet({4, 5}, {1, 2, 3}))")
  expect_identical(.pyrun(local({ x <- Symbol("x"); Partition(list(x), c(1, 2)) })), "Partition({x}, {1, 2})")
  expect_identical(.pyrun(local({ a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); Partition(list(a, b), list(c), 1) })), "Partition({1}, {c}, {a, b})")
})

test_that("partitions: members matches SymPy", {
  expect_identical(.pyrun(local({ .pyT(Partition(c(1, 2, 3), 4)$members) })), "(1, 2, 3, 4)")
  expect_identical(.pyrun(local({ x <- Symbol("x"); .pyT(Partition(list(x), c(1, 2))$members) })), "(1, 2, x)")
  expect_identical(.pyrun(local({ .pyT(Partition(c(3, 1), 2)$members) })), "(1, 2, 3)")
})

test_that("partitions: size matches SymPy", {
  expect_identical(.pyrun(local({ Partition(c(1, 2, 3), 4)$size })), "4")
  expect_identical(.pyrun(local({ Partition(c(1, 2), 3, c(4, 5, 6))$size })), "6")
})

test_that("partitions: partition matches SymPy", {
  expect_identical(.pyrun(local({ Partition(c(1, 2, 3), 4)$partition })), "[[1, 2, 3], [4]]")
  expect_identical(.pyrun(local({ Partition(1, c(2, 3))$partition })), "[[1], [2, 3]]")
  expect_identical(.pyrun(local({ Partition(c(3, 1), 2, c(5, 4))$partition })), "[[1, 3], [2], [4, 5]]")
  expect_identical(.pyrun(local({ x <- Symbol("x"); Partition(list(x), c(1, 2))$partition })), "ERR TypeError: cannot determine truth value of Relational: x < 1")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); Partition(list(x, y))$partition })), "[[x, y]]")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); Partition(list(x), list(y))$partition })), "ERR TypeError: cannot determine truth value of Relational: y < x")
  expect_identical(.pyrun(local({ Partition(c(1, 2), c(3, 4))$partition })), "[[1, 2], [3, 4]]")
})

test_that("partitions: sort matches SymPy", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); a <- Partition(c(1, 2, 3), 4); b <- Partition(c(1, 2), c(3, 4)); cc <- Partition(list(x)); symcompR:::.key_sort(lapply(list(a, b, cc), partition_sort_key), list(a, b, cc)) })), "[Partition({x}), Partition({4}, {1, 2, 3}), Partition({1, 2}, {3, 4})]")
})

test_that("partitions: sortkey matches SymPy", {
  expect_identical(.pyrun(local({ Partition(c(1, 2), 3) < Partition(1, c(2, 3)) })), "True")
})

test_that("partitions: cmp matches SymPy", {
  expect_identical(.pyrun(local({ Partition(c(1, 2, 3), 4) == Partition(c(1, 2), c(3, 4)) })), "False")
  expect_identical(.pyrun(local({ Partition(c(1, 2, 3), 4) <= Partition(c(1, 2), c(3, 4)) })), "True")
  expect_identical(.pyrun(local({ Partition(c(1, 2, 3), 4) > Partition(c(1, 2), c(3, 4)) })), "False")
  expect_identical(.pyrun(local({ Partition(c(1, 2, 3), 4) != Partition(c(1, 2), c(3, 4)) })), "True")
  expect_identical(.pyrun(local({ Partition(c(1, 2, 3), 4) < Partition(c(1, 2), c(3, 4)) })), "True")
  expect_identical(.pyrun(local({ Partition(c(1, 2), c(3, 4)) < Partition(c(1, 2, 3), 4) })), "False")
  expect_identical(.pyrun(local({ Partition(c(1, 2), c(3, 4)) <= Partition(c(1, 2), c(3, 4)) })), "True")
  expect_identical(.pyrun(local({ Partition(c(1, 2), 3) < Partition(c(1, 2), 3, 4) })), "True")
  expect_identical(.pyrun(local({ Partition(c(1, 2), 3) >= Partition(c(1, 2), 3) })), "True")
})

test_that("partitions: add matches SymPy", {
  expect_identical(.pyrun(local({ (Partition(c(1, 2, 3), 4) + 2)$partition })), "[[1, 2], [3, 4]]")
  expect_identical(.pyrun(local({ (Partition(c(1, 2), c(3, 4)) - 1)$partition })), "[[1, 2, 4], [3]]")
  expect_identical(.pyrun(local({ (Partition(c(1, 2, 3), 4) - 1)$partition })), "[[1, 2, 3, 4]]")
  expect_identical(.pyrun(local({ (Partition(c(1, 2, 3), 4) + 1)$partition })), "[[1, 2, 4], [3]]")
  expect_identical(.pyrun(local({ (Partition(c(1, 2), c(3, 4)) + 1)$partition })), "[[1, 2], [3], [4]]")
  expect_identical(.pyrun(local({ Partition(c(1, 2), c(3, 4)) + 7 })), "Partition({1}, {4}, {2, 3})")
  expect_identical(.pyrun(local({ Partition(c(1, 2), c(3, 4)) - 20 })), "Partition({1}, {2}, {3, 4})")
  expect_identical(.pyrun(local({ Partition(c(1, 2), 3) + 1.5 })), "ERR ValueError: 1.5 is not an integer")
  expect_identical(.pyrun(local({ a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); Partition(list(a, b), list(c)) + 1 })), "ERR TypeError: cannot determine truth value of Relational: a < c")
})

test_that("partitions: rank matches SymPy", {
  expect_identical(.pyrun(local({ Partition(c(1, 2, 3), 4)$rank })), "1")
  expect_identical(.pyrun(local({ Partition(c(1, 2), c(3, 4))$rank })), "3")
  expect_identical(.pyrun(local({ Partition(1, 2, 3, 4, 5)$rank })), "51")
  expect_identical(.pyrun(local({ Partition(c(1, 5), c(2, 4), 3)$rank })), "39")
})

test_that("partitions: RGS matches SymPy", {
  expect_identical(.pyrun(local({ .pyT(Partition(c(1, 2, 3), 4)$RGS) })), "(0, 0, 0, 1)")
  expect_identical(.pyrun(local({ .pyT(Partition(c(1, 2), c(3, 4))$RGS) })), "(0, 0, 1, 1)")
  expect_identical(.pyrun(local({ .pyT(Partition(c(1, 2), 3, 4, 5)$RGS) })), "(0, 0, 1, 2, 3)")
  expect_identical(.pyrun(local({ .pyT(Partition(c(5, 1), c(2, 4), 3)$RGS) })), "(0, 1, 2, 1, 0)")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); .pyT(Partition(list(x, y))$RGS) })), "(0, 0)")
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); Partition(list(x, y))$rank })), "0")
})

test_that("partitions: from_rgs matches SymPy", {
  expect_identical(.pyrun(local({ a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); partition_from_rgs(c(0, 1, 2, 0, 1), list(a, b, c, d, e)) })), "Partition({c}, {a, d}, {b, e})")
  expect_identical(.pyrun(local({ a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d"); e <- Symbol("e"); partition_from_rgs(c(0, 1, 2, 0, 1), list(e, d, c, b, a)) })), "Partition({c}, {a, d}, {b, e})")
  expect_identical(.pyrun(local({ partition_from_rgs(c(0, 1, 2, 0, 1), c(1, 3, 2, 4, 5)) })), "Partition({2}, {1, 4}, {3, 5})")
  expect_identical(.pyrun(local({ partition_from_rgs(0:2, 0:1) })), "ERR ValueError: mismatch in rgs and element lengths")
  expect_identical(.pyrun(local({ partition_from_rgs(1:2, 0:1) })), "ERR ValueError: some blocks of the partition were empty.")
  expect_identical(.pyrun(local({ partition_from_rgs(c(0, 0, 1), c(1, 2, 3))$rank })), "1")
})

test_that("partitions: ip matches SymPy", {
  expect_identical(.pyrun(local({ IntegerPartition(0:2) })), "ERR ValueError: All integer summands must be greater than one")
  expect_identical(.pyrun(local({ IntegerPartition(100, 1:2) })), "ERR ValueError: Partition did not add to 100")
  expect_identical(.pyrun(local({ next_lex(IntegerPartition(8, c(1, 3, 4))) })), "[4, 4]")
  expect_identical(.pyrun(local({ IntegerPartition(8, c(1, 3, 4)) == IntegerPartition(c(1, 3, 4)) })), "True")
  expect_identical(.pyrun(local({ IntegerPartition(8, c(1, 3, 4))$integer == IntegerPartition(8, c("1" = 3, "3" = 1, "2" = 1))$integer })), "True")
  expect_identical(.pyrun(local({ IntegerPartition(8, c("1" = 3, "3" = 1, "2" = 1)) })), "[3, 2, 1, 1, 1]")
  expect_identical(.pyrun(local({ .pyL(IntegerPartition(c("1" = 3, "3" = 1, "2" = 1, "5" = 0))$partition) })), "[3, 2, 1, 1, 1]")
  expect_identical(.pyrun(local({ .pyL(IntegerPartition(8, c(1, 3, 4))$conjugate) })), "[3, 2, 2, 1]")
  expect_identical(.pyrun(local({ .pyL(IntegerPartition(c(5, 4, 3, 1, 1))$conjugate) })), "[5, 3, 3, 2, 1]")
  expect_identical(.pyrun(local({ a <- IntegerPartition(8, c(1, 3, 4)); a == next_lex(a) })), "False")
  expect_identical(.pyrun(local({ a <- IntegerPartition(8, c(1, 3, 4)); a <= next_lex(a) })), "True")
  expect_identical(.pyrun(local({ a <- IntegerPartition(8, c(1, 3, 4)); a > next_lex(a) })), "False")
  expect_identical(.pyrun(local({ a <- IntegerPartition(8, c(1, 3, 4)); a != next_lex(a) })), "True")
  expect_identical(.pyrun(local({ a <- IntegerPartition(8, c(1, 3, 4)); a < next_lex(a) })), "True")
  expect_identical(.pyrun(local({ IntegerPartition(c(3, 1)) < IntegerPartition(c(2, 1, 1)) })), "False")
  expect_identical(.pyrun(local({ IntegerPartition(c(3, 1)) >= IntegerPartition(c(2, 1, 1)) })), "True")
  expect_identical(.pyrun(local({ as_ferrers(IntegerPartition(c(1, 2, 3))) })), "###\\n##\\n#")
  expect_identical(.pyrun(local({ as_ferrers(IntegerPartition(c(1, 1, 3)), "o") })), "ooo\\no\\no")
  expect_identical(.pyrun(local({ format(IntegerPartition(c(1, 1, 3))) })), "[3, 1, 1]")
  expect_identical(.pyrun(local({ IntegerPartition(c(1, 1, 3)) })), "[3, 1, 1]")
  expect_identical(.pyrun(local({ list(IntegerPartition(c(1, 1, 3))) })), "[IntegerPartition(5, (3, 1, 1))]")
  expect_identical(.pyrun(local({ .pyL(IntegerPartition(c(1, 1, 3))$partition) })), "[3, 1, 1]")
  expect_identical(.pyrun(local({ as_dict(IntegerPartition(c(1, 1, 2, 3))) })), "{3: 1, 2: 1, 1: 2}")
  expect_identical(.pyrun(local({ as_dict(IntegerPartition(c(3, 3, 3, 3, 2, 1, 1, 1))) })), "{3: 4, 2: 1, 1: 3}")
  expect_identical(.pyrun(local({ prev_lex(IntegerPartition(4)) })), "[3, 1]")
  expect_identical(.pyrun(local({ next_lex(IntegerPartition(c(1, 1, 1, 1))) })), "[2, 1, 1]")
  expect_identical(.pyrun(local({ next_lex(IntegerPartition(4)) })), "[1, 1, 1, 1]")
  expect_identical(.pyrun(local({ prev_lex(IntegerPartition(c(1, 1, 1, 1))) })), "[4]")
})

test_that("partitions: rip matches SymPy", {
  expect_identical(.pyrun(local({ .pyL(random_integer_partition(-1)) })), "ERR ValueError: n must be a positive integer")
  expect_identical(.pyrun(local({ .pyL(random_integer_partition(1)) })), "[1]")
  expect_identical(.pyrun(local({ .pyL(random_integer_partition(10, seed = c(1, 3, 2, 1, 5, 1))) })), "[5, 2, 1, 1, 1]")
  expect_identical(.pyrun(local({ .pyL(random_integer_partition(100, seed = c(1, 1, 12, 1, 2, 1, 85, 1))) })), "[85, 12, 2, 1]")
  expect_identical(.pyrun(local({ .pyL(random_integer_partition(10, seed = c(2, 1, 1, 1, 1, 2))) })), "ERR ValueError: _randint sequence was too short")
})

test_that("partitions: rgs matches SymPy", {
  expect_identical(.pyrun(local({ .pyL(RGS_unrank(-1, 3)) })), "ERR ValueError: Invalid arguments")
  expect_identical(.pyrun(local({ .pyL(RGS_unrank(3, 0)) })), "ERR ValueError: The superset size must be >= 1")
  expect_identical(.pyrun(local({ .pyL(RGS_unrank(10, 1)) })), "ERR ValueError: Invalid arguments")
  expect_identical(.pyrun(local({ RGS_enum(-1) })), "0")
  expect_identical(.pyrun(local({ RGS_enum(1) })), "1")
  expect_identical(.pyrun(local({ RGS_enum(4) })), "15")
  expect_identical(.pyrun(local({ RGS_enum(5) })), "52")
  expect_identical(.pyrun(local({ RGS_enum(6) })), "203")
  expect_identical(.pyrun(local({ .pyL(RGS_unrank(7, 5)) })), "[0, 0, 1, 0, 2]")
  expect_identical(.pyrun(local({ .pyL(RGS_unrank(23, 14)) })), "[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 2, 2]")
  expect_identical(.pyrun(local({ RGS_rank(RGS_unrank(40, 100)) })), "40")
  expect_identical(.pyrun(local({ .pyL(RGS_unrank(14, 4)) })), "[0, 1, 2, 3]")
  expect_identical(.pyrun(local({ .pyL(RGS_unrank(0, 4)) })), "[0, 0, 0, 0]")
  expect_identical(.pyrun(local({ RGS_rank(c(0, 1, 2, 1, 2, 3)) })), "154")
  expect_identical(.pyrun(local({ RGS_rank(c(0, 0, 1, 2, 3)) })), "14")
  expect_identical(.pyrun(local({ RGS_rank(c(0, 1, 2)) })), "4")
  expect_identical(.pyrun(local({ RGS_generalized(6) })), "Matrix([[1, 1, 1, 1, 1, 1, 1], [1, 2, 3, 4, 5, 6, 0], [2, 5, 10, 17, 26, 0, 0], [5, 15, 37, 77, 0, 0, 0], [15, 52, 151, 0, 0, 0, 0], [52, 203, 0, 0, 0, 0, 0], [203, 0, 0, 0, 0, 0, 0]])")
  expect_identical(.pyrun(local({ RGS_generalized(3) })), "Matrix([[1, 1, 1, 1], [1, 2, 3, 0], [2, 5, 0, 0], [5, 0, 0, 0]])")
  expect_identical(.pyrun(local({ RGS_generalized(0) })), "Matrix([[1]])")
  expect_identical(.pyrun(local({ identical(vapply(0:51, function(i) as.integer(as.character(RGS_rank(RGS_unrank(i, 5)))), 0L), 0:51) })), "True")
})

test_that("partitions: next_lex1 matches SymPy", {
  expect_identical(.pyrun(local({ a <- IntegerPartition(1); out <- list(); for (j in 1:1) { out[[j]] <- .pyL(a$partition); a <- next_lex(a) }; out })), "[[1]]")
})

test_that("partitions: prev_lex1 matches SymPy", {
  expect_identical(.pyrun(local({ a <- IntegerPartition(1); out <- list(); for (j in 1:1) { out[[j]] <- .pyL(a$partition); a <- prev_lex(a) }; out })), "[[1]]")
})

test_that("partitions: next_lex2 matches SymPy", {
  expect_identical(.pyrun(local({ a <- IntegerPartition(2); out <- list(); for (j in 1:2) { out[[j]] <- .pyL(a$partition); a <- next_lex(a) }; out })), "[[2], [1, 1]]")
})

test_that("partitions: prev_lex2 matches SymPy", {
  expect_identical(.pyrun(local({ a <- IntegerPartition(2); out <- list(); for (j in 1:2) { out[[j]] <- .pyL(a$partition); a <- prev_lex(a) }; out })), "[[2], [1, 1]]")
})

test_that("partitions: next_lex3 matches SymPy", {
  expect_identical(.pyrun(local({ a <- IntegerPartition(3); out <- list(); for (j in 1:3) { out[[j]] <- .pyL(a$partition); a <- next_lex(a) }; out })), "[[3], [1, 1, 1], [2, 1]]")
})

test_that("partitions: prev_lex3 matches SymPy", {
  expect_identical(.pyrun(local({ a <- IntegerPartition(3); out <- list(); for (j in 1:3) { out[[j]] <- .pyL(a$partition); a <- prev_lex(a) }; out })), "[[3], [2, 1], [1, 1, 1]]")
})

test_that("partitions: next_lex4 matches SymPy", {
  expect_identical(.pyrun(local({ a <- IntegerPartition(4); out <- list(); for (j in 1:5) { out[[j]] <- .pyL(a$partition); a <- next_lex(a) }; out })), "[[4], [1, 1, 1, 1], [2, 1, 1], [2, 2], [3, 1]]")
})

test_that("partitions: prev_lex4 matches SymPy", {
  expect_identical(.pyrun(local({ a <- IntegerPartition(4); out <- list(); for (j in 1:5) { out[[j]] <- .pyL(a$partition); a <- prev_lex(a) }; out })), "[[4], [3, 1], [2, 2], [2, 1, 1], [1, 1, 1, 1]]")
})

test_that("partitions: next_lex5 matches SymPy", {
  expect_identical(.pyrun(local({ a <- IntegerPartition(5); out <- list(); for (j in 1:7) { out[[j]] <- .pyL(a$partition); a <- next_lex(a) }; out })), "[[5], [1, 1, 1, 1, 1], [2, 1, 1, 1], [2, 2, 1], [3, 1, 1], [3, 2], [4, 1]]")
})

test_that("partitions: prev_lex5 matches SymPy", {
  expect_identical(.pyrun(local({ a <- IntegerPartition(5); out <- list(); for (j in 1:7) { out[[j]] <- .pyL(a$partition); a <- prev_lex(a) }; out })), "[[5], [4, 1], [3, 2], [3, 1, 1], [2, 2, 1], [2, 1, 1, 1], [1, 1, 1, 1, 1]]")
})

test_that("partitions: next_lex6 matches SymPy", {
  expect_identical(.pyrun(local({ a <- IntegerPartition(6); out <- list(); for (j in 1:11) { out[[j]] <- .pyL(a$partition); a <- next_lex(a) }; out })), "[[6], [1, 1, 1, 1, 1, 1], [2, 1, 1, 1, 1], [2, 2, 1, 1], [2, 2, 2], [3, 1, 1, 1], [3, 2, 1], [3, 3], [4, 1, 1], [4, 2], [5, 1]]")
})

test_that("partitions: prev_lex6 matches SymPy", {
  expect_identical(.pyrun(local({ a <- IntegerPartition(6); out <- list(); for (j in 1:11) { out[[j]] <- .pyL(a$partition); a <- prev_lex(a) }; out })), "[[6], [5, 1], [4, 2], [4, 1, 1], [3, 3], [3, 2, 1], [3, 1, 1, 1], [2, 2, 2], [2, 2, 1, 1], [2, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1]]")
})

test_that("partitions: next_lex7 matches SymPy", {
  expect_identical(.pyrun(local({ a <- IntegerPartition(7); out <- list(); for (j in 1:15) { out[[j]] <- .pyL(a$partition); a <- next_lex(a) }; out })), "[[7], [1, 1, 1, 1, 1, 1, 1], [2, 1, 1, 1, 1, 1], [2, 2, 1, 1, 1], [2, 2, 2, 1], [3, 1, 1, 1, 1], [3, 2, 1, 1], [3, 2, 2], [3, 3, 1], [4, 1, 1, 1], [4, 2, 1], [4, 3], [5, 1, 1], [5, 2], [6, 1]]")
})

test_that("partitions: prev_lex7 matches SymPy", {
  expect_identical(.pyrun(local({ a <- IntegerPartition(7); out <- list(); for (j in 1:15) { out[[j]] <- .pyL(a$partition); a <- prev_lex(a) }; out })), "[[7], [6, 1], [5, 2], [5, 1, 1], [4, 3], [4, 2, 1], [4, 1, 1, 1], [3, 3, 1], [3, 2, 2], [3, 2, 1, 1], [3, 1, 1, 1, 1], [2, 2, 2, 1], [2, 2, 1, 1, 1], [2, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1, 1]]")
})

test_that("partitions: next_lex8 matches SymPy", {
  expect_identical(.pyrun(local({ a <- IntegerPartition(8); out <- list(); for (j in 1:22) { out[[j]] <- .pyL(a$partition); a <- next_lex(a) }; out })), "[[8], [1, 1, 1, 1, 1, 1, 1, 1], [2, 1, 1, 1, 1, 1, 1], [2, 2, 1, 1, 1, 1], [2, 2, 2, 1, 1], [2, 2, 2, 2], [3, 1, 1, 1, 1, 1], [3, 2, 1, 1, 1], [3, 2, 2, 1], [3, 3, 1, 1], [3, 3, 2], [4, 1, 1, 1, 1], [4, 2, 1, 1], [4, 2, 2], [4, 3, 1], [4, 4], [5, 1, 1, 1], [5, 2, 1], [5, 3], [6, 1, 1], [6, 2], [7, 1]]")
})

test_that("partitions: prev_lex8 matches SymPy", {
  expect_identical(.pyrun(local({ a <- IntegerPartition(8); out <- list(); for (j in 1:22) { out[[j]] <- .pyL(a$partition); a <- prev_lex(a) }; out })), "[[8], [7, 1], [6, 2], [6, 1, 1], [5, 3], [5, 2, 1], [5, 1, 1, 1], [4, 4], [4, 3, 1], [4, 2, 2], [4, 2, 1, 1], [4, 1, 1, 1, 1], [3, 3, 2], [3, 3, 1, 1], [3, 2, 2, 1], [3, 2, 1, 1, 1], [3, 1, 1, 1, 1, 1], [2, 2, 2, 2], [2, 2, 2, 1, 1], [2, 2, 1, 1, 1, 1], [2, 1, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1, 1, 1]]")
})

test_that("partitions: enum1 matches SymPy", {
  expect_identical(.pyrun(local({ p <- Partition(seq_len(1)); lapply(seq_len(as.integer(as.character(RGS_enum(1)))) - 1L, function(k) format(p + k)) })), "['Partition({1})']")
})

test_that("partitions: enum2 matches SymPy", {
  expect_identical(.pyrun(local({ p <- Partition(seq_len(2)); lapply(seq_len(as.integer(as.character(RGS_enum(2)))) - 1L, function(k) format(p + k)) })), "['Partition({1, 2})', 'Partition({1}, {2})']")
})

test_that("partitions: enum3 matches SymPy", {
  expect_identical(.pyrun(local({ p <- Partition(seq_len(3)); lapply(seq_len(as.integer(as.character(RGS_enum(3)))) - 1L, function(k) format(p + k)) })), "['Partition({1, 2, 3})', 'Partition({3}, {1, 2})', 'Partition({2}, {1, 3})', 'Partition({1}, {2, 3})', 'Partition({1}, {2}, {3})']")
})

test_that("partitions: enum4 matches SymPy", {
  expect_identical(.pyrun(local({ p <- Partition(seq_len(4)); lapply(seq_len(as.integer(as.character(RGS_enum(4)))) - 1L, function(k) format(p + k)) })), "['Partition({1, 2, 3, 4})', 'Partition({4}, {1, 2, 3})', 'Partition({3}, {1, 2, 4})', 'Partition({1, 2}, {3, 4})', 'Partition({3}, {4}, {1, 2})', 'Partition({2}, {1, 3, 4})', 'Partition({1, 3}, {2, 4})', 'Partition({2}, {4}, {1, 3})', 'Partition({1, 4}, {2, 3})', 'Partition({1}, {2, 3, 4})', 'Partition({1}, {4}, {2, 3})', 'Partition({2}, {3}, {1, 4})', 'Partition({1}, {3}, {2, 4})', 'Partition({1}, {2}, {3, 4})', 'Partition({1}, {2}, {3}, {4})']")
})

test_that("partitions: enum5 matches SymPy", {
  expect_identical(.pyrun(local({ p <- Partition(seq_len(5)); lapply(seq_len(as.integer(as.character(RGS_enum(5)))) - 1L, function(k) format(p + k)) })), "['Partition({1, 2, 3, 4, 5})', 'Partition({5}, {1, 2, 3, 4})', 'Partition({4}, {1, 2, 3, 5})', 'Partition({4, 5}, {1, 2, 3})', 'Partition({4}, {5}, {1, 2, 3})', 'Partition({3}, {1, 2, 4, 5})', 'Partition({3, 5}, {1, 2, 4})', 'Partition({3}, {5}, {1, 2, 4})', 'Partition({3, 4}, {1, 2, 5})', 'Partition({1, 2}, {3, 4, 5})', 'Partition({5}, {1, 2}, {3, 4})', 'Partition({3}, {4}, {1, 2, 5})', 'Partition({4}, {1, 2}, {3, 5})', 'Partition({3}, {1, 2}, {4, 5})', 'Partition({3}, {4}, {5}, {1, 2})', 'Partition({2}, {1, 3, 4, 5})', 'Partition({2, 5}, {1, 3, 4})', 'Partition({2}, {5}, {1, 3, 4})', 'Partition({2, 4}, {1, 3, 5})', 'Partition({1, 3}, {2, 4, 5})', 'Partition({5}, {1, 3}, {2, 4})', 'Partition({2}, {4}, {1, 3, 5})', 'Partition({4}, {1, 3}, {2, 5})', 'Partition({2}, {1, 3}, {4, 5})', 'Partition({2}, {4}, {5}, {1, 3})', 'Partition({2, 3}, {1, 4, 5})', 'Partition({1, 4}, {2, 3, 5})', 'Partition({5}, {1, 4}, {2, 3})', 'Partition({1, 5}, {2, 3, 4})', 'Partition({1}, {2, 3, 4, 5})', 'Partition({1}, {5}, {2, 3, 4})', 'Partition({4}, {1, 5}, {2, 3})', 'Partition({1}, {4}, {2, 3, 5})', 'Partition({1}, {2, 3}, {4, 5})', 'Partition({1}, {4}, {5}, {2, 3})', 'Partition({2}, {3}, {1, 4, 5})', 'Partition({3}, {1, 4}, {2, 5})', 'Partition({2}, {1, 4}, {3, 5})', 'Partition({2}, {3}, {5}, {1, 4})', 'Partition({3}, {1, 5}, {2, 4})', 'Partition({1}, {3}, {2, 4, 5})', 'Partition({1}, {2, 4}, {3, 5})', 'Partition({1}, {3}, {5}, {2, 4})', 'Partition({2}, {1, 5}, {3, 4})', 'Partition({1}, {2, 5}, {3, 4})', 'Partition({1}, {2}, {3, 4, 5})', 'Partition({1}, {2}, {5}, {3, 4})', 'Partition({2}, {3}, {4}, {1, 5})', 'Partition({1}, {3}, {4}, {2, 5})', 'Partition({1}, {2}, {4}, {3, 5})', 'Partition({1}, {2}, {3}, {4, 5})', 'Partition({1}, {2}, {3}, {4}, {5})']")
})

