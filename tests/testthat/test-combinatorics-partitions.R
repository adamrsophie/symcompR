# Expected values are SymPy's outputs for the same calls (sympy/combinatorics/tests/test_partitions.py);
# the rows come from the battery that compares this port with SymPy.

test_that("partitions: ctor", {
  expect_identical(.pyrun(local({ Partition(Integer(1), Integer(2), Integer(3)) })), "ERR ValueError: Each argument to Partition should be a list, set, or a FiniteSet")
})

test_that("partitions: partition", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); Partition(list(x, y))$partition })), "[[x, y]]")
})

test_that("partitions: RGS", {
  expect_identical(.pyrun(local({ x <- Symbol("x"); y <- Symbol("y"); Partition(list(x, y))$rank })), "0")
})

test_that("partitions: from_rgs", {
  expect_identical(.pyrun(local({ partition_from_rgs(0:2, 0:1) })), "ERR ValueError: mismatch in rgs and element lengths")
  expect_identical(.pyrun(local({ partition_from_rgs(1:2, 0:1) })), "ERR ValueError: some blocks of the partition were empty.")
})

test_that("partitions: ip", {
  expect_identical(.pyrun(local({ IntegerPartition(0:2) })), "ERR ValueError: All integer summands must be greater than one")
  expect_identical(.pyrun(local({ next_lex(IntegerPartition(8, c(1, 3, 4))) })), "[4, 4]")
  expect_identical(.pyrun(local({ IntegerPartition(8, c(1, 3, 4)) == IntegerPartition(c(1, 3, 4)) })), "True")
  expect_identical(.pyrun(local({ IntegerPartition(8, c(1, 3, 4))$integer == IntegerPartition(8, c("1" = 3, "3" = 1, "2" = 1))$integer })), "True")
  expect_identical(.pyrun(local({ IntegerPartition(8, c("1" = 3, "3" = 1, "2" = 1)) })), "[3, 2, 1, 1, 1]")
  expect_identical(.pyrun(local({ .pyL(IntegerPartition(c("1" = 3, "3" = 1, "2" = 1, "5" = 0))$partition) })), "[3, 2, 1, 1, 1]")
  expect_identical(.pyrun(local({ a <- IntegerPartition(8, c(1, 3, 4)); a <= next_lex(a) })), "True")
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
})

test_that("partitions: rip", {
  expect_identical(.pyrun(local({ .pyL(random_integer_partition(-1)) })), "ERR ValueError: n must be a positive integer")
  expect_identical(.pyrun(local({ .pyL(random_integer_partition(100, seed = c(1, 1, 12, 1, 2, 1, 85, 1))) })), "[85, 12, 2, 1]")
  expect_identical(.pyrun(local({ .pyL(random_integer_partition(10, seed = c(2, 1, 1, 1, 1, 2))) })), "ERR ValueError: _randint sequence was too short")
})

test_that("partitions: rgs", {
  expect_identical(.pyrun(local({ .pyL(RGS_unrank(-1, 3)) })), "ERR ValueError: Invalid arguments")
  expect_identical(.pyrun(local({ .pyL(RGS_unrank(3, 0)) })), "ERR ValueError: The superset size must be >= 1")
  expect_identical(.pyrun(local({ .pyL(RGS_unrank(10, 1)) })), "ERR ValueError: Invalid arguments")
  expect_identical(.pyrun(local({ RGS_enum(1) })), "1")
  expect_identical(.pyrun(local({ RGS_enum(4) })), "15")
  expect_identical(.pyrun(local({ RGS_rank(c(0, 1, 2, 1, 2, 3)) })), "154")
  expect_identical(.pyrun(local({ RGS_rank(c(0, 1, 2)) })), "4")
  expect_identical(.pyrun(local({ RGS_generalized(0) })), "Matrix([[1]])")
})

test_that("partitions: prev_lex1", {
  expect_identical(.pyrun(local({ a <- IntegerPartition(1); out <- list(); for (j in 1:1) { out[[j]] <- .pyL(a$partition); a <- prev_lex(a) }; out })), "[[1]]")
})

test_that("partitions: next_lex3", {
  expect_identical(.pyrun(local({ a <- IntegerPartition(3); out <- list(); for (j in 1:3) { out[[j]] <- .pyL(a$partition); a <- next_lex(a) }; out })), "[[3], [1, 1, 1], [2, 1]]")
})

