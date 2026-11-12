# Reference values from SymPy.
L <- function(v) paste0("[", paste(vapply(v, sstr, ""), collapse = ", "), "]")

test_that("Range canonical form, elements and membership", {
  expect_equal(sstr(Range(10)), "Range(0, 10, 1)")
  expect_equal(sstr(Range(10, 1, -2)), "Range(10, 0, -2)")
  expect_equal(sstr(Range(5, 1)), "Range(0, 0, 1)")
  expect_equal(sstr(Range(2, 9, 3)), "Range(2, 11, 3)")
  expect_equal(sstr(Range(-Inf, 0, 2)), "Range(-oo, 0, 2)")
  expect_equal(L(sequence_terms(Range(1, 10, 3))), "[1, 4, 7]")
  expect_equal(L(sequence_terms(Range(10, 1, -3))), "[10, 7, 4]")
  expect_true(contains(Range(1, 10, 3), 7))
  expect_false(contains(Range(1, 10, 3), 8))
  expect_true(contains(Range(0, Inf, 2), 10))
  expect_equal(latex(Range(1, 10, 3)), "\\left\\{1, 4, 7\\right\\}")
  expect_equal(latex(Range(0, Inf)), "\\left\\{0, 1, \\ldots\\right\\}")
  expect_equal(sstr(as_relational(Range(1, 4), Symbol("n"))), "(n >= 1) & (n <= 3) & Eq(Mod(n, 1), 0)")
})

test_that("Range intersections", {
  expect_equal(sstr(Intersection(Range(1, 10, 3), Interval(3, 8))), "Range(4, 10, 3)")
  expect_equal(sstr(Intersection(Range(1, 10, 3), Range(2, 20, 2))), "Range(4, 10, 6)")
  expect_equal(sstr(Intersection(Range(0, 10, 2), Range(0, 10, 3))), "Range(0, 12, 6)")
  expect_equal(sstr(Intersection(Interval(0, 5), Range(-3, 3))), "Range(0, 3, 1)")
  expect_equal(sstr(Intersection(Range(0, Inf, 2), Range(1, Inf, 3))), "Range(4, oo, 6)")
  expect_identical(Intersection(Range(1, 10, 2), Range(0, 10, 2)), S_EmptySet)
  expect_equal(sstr(Intersection(FiniteSet(1, 2, 4, 5), Range(1, 10, 3))), "{1, 4}")
  expect_equal(sstr(DisjointUnion(FiniteSet(1, 2), FiniteSet(2, 3))), "DisjointUnion({1, 2}, {2, 3})")
  expect_true(contains(DisjointUnion(FiniteSet(1, 2), FiniteSet(2, 3)), Tuple(1, 0)))
})

test_that("p13 h2 fixes: empty unbounded ranges, symbolic membership, unbounded intersections", {
  expect_equal(sstr(Range(Inf, -Inf, 1)), "Range(0, 0, 1)")
  expect_equal(sstr(Range(Inf, Inf)), "Range(0, 0, 1)")
  expect_equal(sstr(Range(-Inf, Inf, -1)), "Range(0, 0, 1)")
  expect_equal(sstr(Range(Inf, -Inf, -1)), "Range(oo, -oo, -1)")
  expect_false(contains(Range(-Inf, Inf, -1), 0))
  n <- Symbol("n", integer = TRUE)
  expect_equal(sstr(contains(Range(0, n), 0)), "Contains(0, Range(0, n, 1))")
  expect_false(contains(Range(0, n), -1))
  expect_equal(sstr(Intersection(Range(-Inf, Inf), Interval(0, 3))), "Range(0, 4, 1)")
  expect_equal(sstr(Intersection(Range(-Inf, Inf), Range(1, 10, 3))), "Range(1, 10, 3)")
  expect_equal(sstr(Intersection(Range(0, Inf), Interval(0, Integer(2)^Rational(1, "2147483648")))), "Range(0, 2, 1)")
  # a cancelling bound is settled at higher precision, as SymPy's get_integer_part
  H <- Integer(2)^70
  v <- exp(H*(sqrt(H^2 + 1) - H)) + Rational(3, 5)
  expect_equal(sstr(Intersection(Range(0, Inf), Interval(0, v))), "Range(0, 3, 1)")
})

test_that("p14 h2 fixes: unbounded range intersections and certified bounds", {
  expect_equal(sstr(Intersection(Range(-Inf, Inf), Range(-Inf, 5))), "Range(-oo, 5, 1)")
  expect_equal(sstr(Intersection(Range(-Inf, Inf), Range(Inf, -Inf, -1))), "Range(oo, -oo, -1)")
  expect_equal(sstr(Intersection(Range(-Inf, Inf), Range(-Inf, 6, 2))), "Range(-oo, 6, 2)")
  expect_equal(sstr(Intersection(Range(-Inf, 10, 3), Range(-Inf, 20, 2))), "Range(-oo, 10, 6)")
  expect_equal(sstr(Range(Inf, Inf, 2)), "Range(0, 0, 1)")
  # a huge sine argument is reduced at the precision it needs
  expect_equal(sstr(Intersection(Range(0, Inf), Interval(0, 2 + sin(Integer("9007199254741001"))))), "Range(0, 3, 1)")
  expect_equal(sstr(Intersection(Range(0, Inf), Interval(0, sqrt(Integer(7)) + sin(Integer(1))))), "Range(0, 4, 1)")
})

test_that("p15 review fix: membership in a range starting at infinity", {
  expect_true(contains(Range(Inf, 0, -2), 2))
  expect_false(contains(Range(Inf, 0, -2), 3))
})
