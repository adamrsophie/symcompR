skip_unless_full_tests()

# Expected values are the outputs of the same calls in SymPy (physics.control);
# the cases come from its test files and run through the helpers in
# helper-physics-misc-battery.R.

test_that("lti test_state_space_rewrite_from_pr_27651 matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("a, x, b, c, s, g, d,    p, k, tau, zeta, wn, T, z"); a <- .py_getitem(.tmpu_, 0L); x <- .py_getitem(.tmpu_, 1L); b <- .py_getitem(.tmpu_, 2L); c <- .py_getitem(.tmpu_, 3L); s <- .py_getitem(.tmpu_, 4L); g <- .py_getitem(.tmpu_, 5L); d <- .py_getitem(.tmpu_, 6L); p <- .py_getitem(.tmpu_, 7L); k <- .py_getitem(.tmpu_, 8L); tau <- .py_getitem(.tmpu_, 9L); zeta <- .py_getitem(.tmpu_, 10L); wn <- .py_getitem(.tmpu_, 11L); T <- .py_getitem(.tmpu_, 12L); z <- .py_getitem(.tmpu_, 13L)
})
.py_try({
.tmpu_ <- .py_symbols("a0:4, b0:5, c0:4, d0:4"); a0 <- .py_getitem(.tmpu_, 0L); a1 <- .py_getitem(.tmpu_, 1L); a2 <- .py_getitem(.tmpu_, 2L); a3 <- .py_getitem(.tmpu_, 3L); b0 <- .py_getitem(.tmpu_, 4L); b1 <- .py_getitem(.tmpu_, 5L); b2 <- .py_getitem(.tmpu_, 6L); b3 <- .py_getitem(.tmpu_, 7L); b4 <- .py_getitem(.tmpu_, 8L); c0 <- .py_getitem(.tmpu_, 9L); c1 <- .py_getitem(.tmpu_, 10L); c2 <- .py_getitem(.tmpu_, 11L); c3 <- .py_getitem(.tmpu_, 12L); d0 <- .py_getitem(.tmpu_, 13L); d1 <- .py_getitem(.tmpu_, 14L); d2 <- .py_getitem(.tmpu_, 15L); d3 <- .py_getitem(.tmpu_, 16L)
})
.py_try({
TF1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
TF2 <- TransferFunction(k, 1L, s)
})
.py_try({
TF3 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
TF4 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(tau, s), 1L), s)
})
.py_try({
TF5 <- TransferFunction(k, p, s)
})
.py_try({
TF6 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(k, s), p), s)
})
p_are_floats_equals <- function(n1, n2) {
  "\n    Helper function to compare two floats.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return(.py_method(Float(n1), "is_same", Float(n2), .py_isclose))
}
p_are_complex_equals <- function(n1, n2) {
  "\n    Helper function to compare two complex numbers.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return((.py_truthy(p_are_floats_equals(re(n1), re(n2))) && .py_truthy(p_are_floats_equals(im(n1), im(n2)))))
}
.py_try({
.tmpu_ <- .py_symbols("c0, c1, c2", real = TRUE); c0 <- .py_getitem(.tmpu_, 0L); c1 <- .py_getitem(.tmpu_, 1L); c2 <- .py_getitem(.tmpu_, 2L)
})
.py_try({
.tmpu_ <- .py_symbols("a0, a1, a2, a3", real = TRUE); a0 <- .py_getitem(.tmpu_, 0L); a1 <- .py_getitem(.tmpu_, 1L); a2 <- .py_getitem(.tmpu_, 2L); a3 <- .py_getitem(.tmpu_, 3L)
})
.py_try({
d <- .py_symbols("d", real = TRUE)
})
.py_try({
s <- .py_symbols("s")
})
.py_try({
t <- .py_add(TransferFunction(d, 1L, s), TransferFunction(.py_add(.py_add(.py_mul(c2, .py_pow(s, 2L)), .py_mul(c1, .py_pow(s, 1L))), c0), .py_add(.py_add(.py_add(.py_add(.py_pow(s, 4L), .py_mul(a3, .py_pow(s, 3L))), .py_mul(a2, .py_pow(s, 2L))), .py_mul(a1, .py_pow(s, 1L))), a0), s))
})
.py_try({
expected <- StateSpace(.py_Matrix(list(list(0L, 1L, 0L, 0L), list(0L, 0L, 1L, 0L), list(0L, 0L, 0L, 1L), list(.py_neg(a0), .py_neg(a1), .py_neg(a2), .py_neg(a3)))), .py_Matrix(list(list(0L), list(0L), list(0L), list(1L))), .py_Matrix(list(list(c0, c1, c2, 0L))), .py_Matrix(list(list(d))))
})
.py_try({
.py_rec("lti:test_state_space_rewrite_from_pr_27651:68", function() .py_method(.py_method(t, "rewrite", StateSpace), "doit"))
})
.py_try({
.py_rec("lti:test_state_space_rewrite_from_pr_27651:69", function() .py_method(.py_method(t, "rewrite", StateSpace), "doit"))
}))---")
  expect_identical(got[["lti:test_state_space_rewrite_from_pr_27651:68"]], "StateSpace(Matrix([\\n[0,   0,   0,   0,   0],\\n[0,   0,   1,   0,   0],\\n[0,   0,   0,   1,   0],\\n[0,   0,   0,   0,   1],\\n[0, -a0, -a1, -a2, -a3]]), Matrix([\\n[0],\\n[0],\\n[0],\\n[0],\\n[1]]), Matrix([[0, c0, c1, c2, 0]]), Matrix([[d]]))")
  expect_identical(got[["lti:test_state_space_rewrite_from_pr_27651:69"]], "StateSpace(Matrix([\\n[0,   0,   0,   0,   0],\\n[0,   0,   1,   0,   0],\\n[0,   0,   0,   1,   0],\\n[0,   0,   0,   0,   1],\\n[0, -a0, -a1, -a2, -a3]]), Matrix([\\n[0],\\n[0],\\n[0],\\n[0],\\n[1]]), Matrix([[0, c0, c1, c2, 0]]), Matrix([[d]]))")
})

test_that("lti test_create_transfer_function matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("a, x, b, c, s, g, d,    p, k, tau, zeta, wn, T, z"); a <- .py_getitem(.tmpu_, 0L); x <- .py_getitem(.tmpu_, 1L); b <- .py_getitem(.tmpu_, 2L); c <- .py_getitem(.tmpu_, 3L); s <- .py_getitem(.tmpu_, 4L); g <- .py_getitem(.tmpu_, 5L); d <- .py_getitem(.tmpu_, 6L); p <- .py_getitem(.tmpu_, 7L); k <- .py_getitem(.tmpu_, 8L); tau <- .py_getitem(.tmpu_, 9L); zeta <- .py_getitem(.tmpu_, 10L); wn <- .py_getitem(.tmpu_, 11L); T <- .py_getitem(.tmpu_, 12L); z <- .py_getitem(.tmpu_, 13L)
})
.py_try({
.tmpu_ <- .py_symbols("a0:4, b0:5, c0:4, d0:4"); a0 <- .py_getitem(.tmpu_, 0L); a1 <- .py_getitem(.tmpu_, 1L); a2 <- .py_getitem(.tmpu_, 2L); a3 <- .py_getitem(.tmpu_, 3L); b0 <- .py_getitem(.tmpu_, 4L); b1 <- .py_getitem(.tmpu_, 5L); b2 <- .py_getitem(.tmpu_, 6L); b3 <- .py_getitem(.tmpu_, 7L); b4 <- .py_getitem(.tmpu_, 8L); c0 <- .py_getitem(.tmpu_, 9L); c1 <- .py_getitem(.tmpu_, 10L); c2 <- .py_getitem(.tmpu_, 11L); c3 <- .py_getitem(.tmpu_, 12L); d0 <- .py_getitem(.tmpu_, 13L); d1 <- .py_getitem(.tmpu_, 14L); d2 <- .py_getitem(.tmpu_, 15L); d3 <- .py_getitem(.tmpu_, 16L)
})
.py_try({
TF1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
TF2 <- TransferFunction(k, 1L, s)
})
.py_try({
TF3 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
TF4 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(tau, s), 1L), s)
})
.py_try({
TF5 <- TransferFunction(k, p, s)
})
.py_try({
TF6 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(k, s), p), s)
})
p_are_floats_equals <- function(n1, n2) {
  "\n    Helper function to compare two floats.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return(.py_method(Float(n1), "is_same", Float(n2), .py_isclose))
}
p_are_complex_equals <- function(n1, n2) {
  "\n    Helper function to compare two complex numbers.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return((.py_truthy(p_are_floats_equals(re(n1), re(n2))) && .py_truthy(p_are_floats_equals(im(n1), im(n2)))))
}
.py_try({
cont_tf1 <- create_transfer_function(.py_add(s, 1L), .py_add(.py_pow(s, 2L), 2L), s)
})
.py_try({
.py_rec("lti:test_create_transfer_function:74", function() .py_isinstance(cont_tf1, TransferFunction))
})
.py_try({
cont_tf2 <- create_transfer_function(.py_add(s, 1L), .py_add(.py_pow(s, 2L), 2L), s, sampling_time = 0L)
})
.py_try({
.py_rec("lti:test_create_transfer_function:77", function() .py_isinstance(cont_tf2, TransferFunction))
})
.py_try({
disc_tf1 <- create_transfer_function(z, .py_add(z, 1L), z, .py_flt(0.1))
})
.py_try({
.py_rec("lti:test_create_transfer_function:80", function() .py_isinstance(disc_tf1, DiscreteTransferFunction))
})
.py_try({
disc_tf2 <- create_transfer_function(z, .py_add(z, 1L), z, T)
})
.py_try({
.py_rec("lti:test_create_transfer_function:83", function() .py_isinstance(disc_tf2, DiscreteTransferFunction))
}))---")
  expect_identical(got[["lti:test_create_transfer_function:74"]], "True")
  expect_identical(got[["lti:test_create_transfer_function:77"]], "True")
  expect_identical(got[["lti:test_create_transfer_function:80"]], "True")
  expect_identical(got[["lti:test_create_transfer_function:83"]], "True")
})

test_that("lti test_TransferFunctionBase matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("a, x, b, c, s, g, d,    p, k, tau, zeta, wn, T, z"); a <- .py_getitem(.tmpu_, 0L); x <- .py_getitem(.tmpu_, 1L); b <- .py_getitem(.tmpu_, 2L); c <- .py_getitem(.tmpu_, 3L); s <- .py_getitem(.tmpu_, 4L); g <- .py_getitem(.tmpu_, 5L); d <- .py_getitem(.tmpu_, 6L); p <- .py_getitem(.tmpu_, 7L); k <- .py_getitem(.tmpu_, 8L); tau <- .py_getitem(.tmpu_, 9L); zeta <- .py_getitem(.tmpu_, 10L); wn <- .py_getitem(.tmpu_, 11L); T <- .py_getitem(.tmpu_, 12L); z <- .py_getitem(.tmpu_, 13L)
})
.py_try({
.tmpu_ <- .py_symbols("a0:4, b0:5, c0:4, d0:4"); a0 <- .py_getitem(.tmpu_, 0L); a1 <- .py_getitem(.tmpu_, 1L); a2 <- .py_getitem(.tmpu_, 2L); a3 <- .py_getitem(.tmpu_, 3L); b0 <- .py_getitem(.tmpu_, 4L); b1 <- .py_getitem(.tmpu_, 5L); b2 <- .py_getitem(.tmpu_, 6L); b3 <- .py_getitem(.tmpu_, 7L); b4 <- .py_getitem(.tmpu_, 8L); c0 <- .py_getitem(.tmpu_, 9L); c1 <- .py_getitem(.tmpu_, 10L); c2 <- .py_getitem(.tmpu_, 11L); c3 <- .py_getitem(.tmpu_, 12L); d0 <- .py_getitem(.tmpu_, 13L); d1 <- .py_getitem(.tmpu_, 14L); d2 <- .py_getitem(.tmpu_, 15L); d3 <- .py_getitem(.tmpu_, 16L)
})
.py_try({
TF1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
TF2 <- TransferFunction(k, 1L, s)
})
.py_try({
TF3 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
TF4 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(tau, s), 1L), s)
})
.py_try({
TF5 <- TransferFunction(k, p, s)
})
.py_try({
TF6 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(k, s), p), s)
})
p_are_floats_equals <- function(n1, n2) {
  "\n    Helper function to compare two floats.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return(.py_method(Float(n1), "is_same", Float(n2), .py_isclose))
}
p_are_complex_equals <- function(n1, n2) {
  "\n    Helper function to compare two complex numbers.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return((.py_truthy(p_are_floats_equals(re(n1), re(n2))) && .py_truthy(p_are_floats_equals(im(n1), im(n2)))))
}
.py_try({
.py_rec("lti:test_TransferFunctionBase:86", function() .py_method(TransferFunctionBase, "from_rational_expression", .py_div(2L, s)))
})
.py_try({
.py_rec("lti:test_TransferFunctionBase:88", function() .py_method(TransferFunctionBase, "from_coeff_lists", list(1L), list(1L), s))
})
.py_try({
.py_rec("lti:test_TransferFunctionBase:90", function() .py_method(TransferFunctionBase, "from_zpk", list(1L), list(1L), 1L, s))
}))---")
  expect_identical(got[["lti:test_TransferFunctionBase:86"]], "ERR NotImplementedError: \\n                The TransferFunctionBase class is not meant to be used directly.\\n                ")
  expect_identical(got[["lti:test_TransferFunctionBase:88"]], "ERR NotImplementedError: \\n                The TransferFunctionBase class is not meant to be used directly.\\n                ")
  expect_identical(got[["lti:test_TransferFunctionBase:90"]], "ERR NotImplementedError: \\n                The TransferFunctionBase class is not meant to be used directly.\\n                ")
})

test_that("lti test_TransferFunction_construction matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("a, x, b, c, s, g, d,    p, k, tau, zeta, wn, T, z"); a <- .py_getitem(.tmpu_, 0L); x <- .py_getitem(.tmpu_, 1L); b <- .py_getitem(.tmpu_, 2L); c <- .py_getitem(.tmpu_, 3L); s <- .py_getitem(.tmpu_, 4L); g <- .py_getitem(.tmpu_, 5L); d <- .py_getitem(.tmpu_, 6L); p <- .py_getitem(.tmpu_, 7L); k <- .py_getitem(.tmpu_, 8L); tau <- .py_getitem(.tmpu_, 9L); zeta <- .py_getitem(.tmpu_, 10L); wn <- .py_getitem(.tmpu_, 11L); T <- .py_getitem(.tmpu_, 12L); z <- .py_getitem(.tmpu_, 13L)
})
.py_try({
.tmpu_ <- .py_symbols("a0:4, b0:5, c0:4, d0:4"); a0 <- .py_getitem(.tmpu_, 0L); a1 <- .py_getitem(.tmpu_, 1L); a2 <- .py_getitem(.tmpu_, 2L); a3 <- .py_getitem(.tmpu_, 3L); b0 <- .py_getitem(.tmpu_, 4L); b1 <- .py_getitem(.tmpu_, 5L); b2 <- .py_getitem(.tmpu_, 6L); b3 <- .py_getitem(.tmpu_, 7L); b4 <- .py_getitem(.tmpu_, 8L); c0 <- .py_getitem(.tmpu_, 9L); c1 <- .py_getitem(.tmpu_, 10L); c2 <- .py_getitem(.tmpu_, 11L); c3 <- .py_getitem(.tmpu_, 12L); d0 <- .py_getitem(.tmpu_, 13L); d1 <- .py_getitem(.tmpu_, 14L); d2 <- .py_getitem(.tmpu_, 15L); d3 <- .py_getitem(.tmpu_, 16L)
})
.py_try({
TF1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
TF2 <- TransferFunction(k, 1L, s)
})
.py_try({
TF3 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
TF4 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(tau, s), 1L), s)
})
.py_try({
TF5 <- TransferFunction(k, p, s)
})
.py_try({
TF6 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(k, s), p), s)
})
p_are_floats_equals <- function(n1, n2) {
  "\n    Helper function to compare two floats.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return(.py_method(Float(n1), "is_same", Float(n2), .py_isclose))
}
p_are_complex_equals <- function(n1, n2) {
  "\n    Helper function to compare two complex numbers.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return((.py_truthy(p_are_floats_equals(re(n1), re(n2))) && .py_truthy(p_are_floats_equals(im(n1), im(n2)))))
}
.py_try({
tf <- TransferFunction(.py_add(s, 1L), .py_add(.py_add(.py_pow(s, 2L), s), 1L), s)
})
.py_try({
.py_rec("lti:test_TransferFunction_construction:95", function() .py_attr(tf, "num"))
})
.py_try({
.py_rec("lti:test_TransferFunction_construction:96", function() .py_attr(tf, "den"))
})
.py_try({
.py_rec("lti:test_TransferFunction_construction:97", function() .py_attr(tf, "args"))
})
.py_try({
tf1 <- TransferFunction(.py_add(s, 4L), .py_sub(s, 5L), s)
})
.py_try({
.py_rec("lti:test_TransferFunction_construction:100", function() .py_attr(tf1, "num"))
})
.py_try({
.py_rec("lti:test_TransferFunction_construction:101", function() .py_attr(tf1, "den"))
})
.py_try({
.py_rec("lti:test_TransferFunction_construction:102", function() .py_attr(tf1, "args"))
})
.py_try({
tf2 <- TransferFunction(.py_add(p, 3L), .py_sub(.py_pow(p, 2L), 9L), p)
})
.py_try({
.py_rec("lti:test_TransferFunction_construction:106", function() .py_attr(tf2, "num"))
})
.py_try({
.py_rec("lti:test_TransferFunction_construction:107", function() .py_attr(tf2, "den"))
})
.py_try({
.py_rec("lti:test_TransferFunction_construction:108", function() .py_attr(tf2, "args"))
})
.py_try({
tf3 <- TransferFunction(.py_add(.py_add(.py_pow(p, 3L), .py_mul(5L, .py_pow(p, 2L))), 4L), .py_add(.py_add(.py_pow(p, 4L), .py_mul(3L, p)), 1L), p)
})
.py_try({
.py_rec("lti:test_TransferFunction_construction:111", function() .py_attr(tf3, "args"))
})
.py_try({
tf4 <- TransferFunction(.py_mul(.py_add(s, 3L), .py_sub(s, 1L)), .py_mul(.py_sub(s, 1L), .py_add(s, 5L)), s)
})
.py_try({
.py_rec("lti:test_TransferFunction_construction:115", function() .py_attr(tf4, "den"))
})
.py_try({
.py_rec("lti:test_TransferFunction_construction:116", function() .py_attr(tf4, "args"))
})
.py_try({
tf4_ <- TransferFunction(.py_add(p, 2L), .py_add(p, 2L), p)
})
.py_try({
.py_rec("lti:test_TransferFunction_construction:119", function() .py_attr(tf4_, "args"))
})
.py_try({
tf5 <- TransferFunction(.py_sub(s, 1L), .py_sub(4L, p), s)
})
.py_try({
.py_rec("lti:test_TransferFunction_construction:122", function() .py_attr(tf5, "args"))
})
.py_try({
tf5_ <- TransferFunction(.py_sub(s, 1L), .py_sub(s, 1L), s)
})
.py_try({
.py_rec("lti:test_TransferFunction_construction:125", function() .py_attr(tf5_, "args"))
})
.py_try({
tf6 <- TransferFunction(5L, 6L, s)
})
.py_try({
.py_rec("lti:test_TransferFunction_construction:128", function() .py_attr(tf6, "num"))
})
.py_try({
.py_rec("lti:test_TransferFunction_construction:129", function() .py_attr(tf6, "den"))
})
.py_try({
.py_rec("lti:test_TransferFunction_construction:130", function() .py_attr(tf6, "args"))
})
.py_try({
tf6_ <- TransferFunction(.py_div(1L, 2L), 4L, s)
})
.py_try({
.py_rec("lti:test_TransferFunction_construction:133", function() .py_attr(tf6_, "num"))
})
.py_try({
.py_rec("lti:test_TransferFunction_construction:134", function() .py_attr(tf6_, "den"))
})
.py_try({
.py_rec("lti:test_TransferFunction_construction:135", function() .py_attr(tf6_, "args"))
})
.py_try({
tf7 <- TransferFunction(.py_add(.py_add(.py_mul(3L, .py_pow(s, 2L)), .py_mul(2L, p)), .py_mul(4L, s)), .py_add(.py_mul(8L, .py_pow(p, 2L)), .py_mul(7L, s)), s)
})
.py_try({
tf8 <- TransferFunction(.py_add(.py_add(.py_mul(3L, .py_pow(s, 2L)), .py_mul(2L, p)), .py_mul(4L, s)), .py_add(.py_mul(8L, .py_pow(p, 2L)), .py_mul(7L, s)), p)
})
.py_try({
.py_rec("lti:test_TransferFunction_construction:139", function() tf7)
})
.py_try({
tf7_ <- TransferFunction(.py_add(.py_add(.py_mul(a0, s), .py_mul(a1, .py_pow(s, 2L))), .py_mul(a2, .py_pow(s, 3L))), .py_sub(.py_mul(b0, p), .py_mul(b1, s)), s)
})
.py_try({
tf8_ <- TransferFunction(.py_add(.py_add(.py_mul(a0, s), .py_mul(a1, .py_pow(s, 2L))), .py_mul(a2, .py_pow(s, 3L))), .py_sub(.py_mul(b0, p), .py_mul(b1, s)), s)
})
.py_try({
.py_rec("lti:test_TransferFunction_construction:143", function() tf7_)
})
.py_try({
.py_rec("lti:test_TransferFunction_construction:144", function() (.py_truthy(.py_eq(.py_neg(.py_neg(tf7_)), tf7_)) && .py_truthy(.py_eq(tf7_, .py_neg(.py_neg(.py_neg(.py_neg(tf7_))))))))
})
.py_try({
tf9 <- TransferFunction(.py_add(.py_add(.py_add(.py_mul(a, .py_pow(s, 3L)), .py_mul(b, .py_pow(s, 2L))), .py_mul(g, s)), d), .py_add(.py_add(.py_mul(d, p), .py_mul(g, .py_pow(p, 2L))), .py_mul(g, s)), s)
})
.py_try({
.py_rec("lti:test_TransferFunction_construction:147", function() .py_attr(tf9, "args"))
})
.py_try({
tf10 <- TransferFunction(.py_add(.py_pow(p, 3L), d), .py_add(.py_add(.py_mul(g, .py_pow(s, 2L)), .py_mul(d, s)), a), p)
})
.py_try({
tf10_ <- TransferFunction(.py_add(.py_pow(p, 3L), d), .py_add(.py_add(.py_mul(g, .py_pow(s, 2L)), .py_mul(d, s)), a), p)
})
.py_try({
.py_rec("lti:test_TransferFunction_construction:151", function() .py_attr(tf10, "args"))
})
.py_try({
.py_rec("lti:test_TransferFunction_construction:152", function() tf10_)
})
.py_try({
tf11 <- TransferFunction(.py_add(.py_mul(a1, s), a0), .py_add(.py_add(.py_mul(b2, .py_pow(s, 2L)), .py_mul(b1, s)), b0), s)
})
.py_try({
.py_rec("lti:test_TransferFunction_construction:155", function() .py_attr(tf11, "num"))
})
.py_try({
.py_rec("lti:test_TransferFunction_construction:156", function() .py_attr(tf11, "den"))
})
.py_try({
.py_rec("lti:test_TransferFunction_construction:157", function() .py_attr(tf11, "args"))
})
.py_try({
tf12 <- TransferFunction(0L, .py_add(.py_sub(.py_pow(p, 2L), p), 1L), p)
})
.py_try({
.py_rec("lti:test_TransferFunction_construction:161", function() .py_attr(tf12, "args"))
})
.py_try({
tf13 <- TransferFunction(0L, 1L, s)
})
.py_try({
.py_rec("lti:test_TransferFunction_construction:164", function() .py_attr(tf13, "args"))
})
.py_try({
tf14 <- TransferFunction(.py_sub(.py_add(.py_mul(a0, .py_pow(s, .py_flt(0.5))), .py_mul(a2, .py_pow(s, .py_flt(0.6)))), a1), .py_mul(a1, .py_pow(p, .py_flt(-8.7))), s)
})
.py_try({
.py_rec("lti:test_TransferFunction_construction:168", function() .py_attr(tf14, "args"))
})
.py_try({
tf15 <- TransferFunction(.py_add(.py_mul(.py_pow(a2, 2L), .py_pow(p, .py_div(1L, 4L))), .py_mul(a1, .py_pow(s, .py_div(-4L, 5L)))), .py_sub(.py_mul(a0, s), p), p)
})
.py_try({
.py_rec("lti:test_TransferFunction_construction:171", function() .py_attr(tf15, "args"))
})
.py_try({
.tmpu_ <- .py_symbols("omega_o, k_p, k_o, k_i"); omega_o <- .py_getitem(.tmpu_, 0L); k_p <- .py_getitem(.tmpu_, 1L); k_o <- .py_getitem(.tmpu_, 2L); k_i <- .py_getitem(.tmpu_, 3L)
})
.py_try({
tf18 <- TransferFunction(.py_add(.py_add(k_p, .py_mul(k_o, s)), .py_div(k_i, s)), .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(2L, omega_o), s)), .py_pow(omega_o, 2L)), s)
})
.py_try({
.py_rec("lti:test_TransferFunction_construction:176", function() .py_attr(tf18, "num"))
})
.py_try({
.py_rec("lti:test_TransferFunction_construction:177", function() .py_attr(tf18, "args"))
})
.py_try({
.py_rec("lti:test_TransferFunction_construction:181", function() TransferFunction(4L, 0L, s))
})
.py_try({
.py_rec("lti:test_TransferFunction_construction:182", function() TransferFunction(s, 0L, s))
})
.py_try({
.py_rec("lti:test_TransferFunction_construction:183", function() TransferFunction(0L, 0L, s))
})
.py_try({
.py_rec("lti:test_TransferFunction_construction:185", function() TransferFunction(.py_Matrix(list(1L, 2L, 3L)), s, s))
})
.py_try({
.py_rec("lti:test_TransferFunction_construction:187", function() TransferFunction(.py_sub(.py_add(.py_pow(s, 2L), .py_mul(2L, s)), 1L), .py_add(s, 3L), 3L))
})
.py_try({
.py_rec("lti:test_TransferFunction_construction:188", function() TransferFunction(.py_add(p, 1L), .py_sub(5L, p), 4L))
})
.py_try({
.py_rec("lti:test_TransferFunction_construction:189", function() TransferFunction(3L, 4L, 8L))
}))---")
  expect_identical(got[["lti:test_TransferFunction_construction:95"]], "s + 1")
  expect_identical(got[["lti:test_TransferFunction_construction:96"]], "s^2 + s + 1")
  expect_identical(got[["lti:test_TransferFunction_construction:97"]], "(s + 1, s^2 + s + 1, s)")
  expect_identical(got[["lti:test_TransferFunction_construction:100"]], "s + 4")
  expect_identical(got[["lti:test_TransferFunction_construction:101"]], "s - 5")
  expect_identical(got[["lti:test_TransferFunction_construction:102"]], "(s + 4, s - 5, s)")
  expect_identical(got[["lti:test_TransferFunction_construction:106"]], "p + 3")
  expect_identical(got[["lti:test_TransferFunction_construction:107"]], "p^2 - 9")
  expect_identical(got[["lti:test_TransferFunction_construction:108"]], "(p + 3, p^2 - 9, p)")
  expect_identical(got[["lti:test_TransferFunction_construction:111"]], "(p^3 + 5*p^2 + 4, p^4 + 3*p + 1, p)")
  expect_identical(got[["lti:test_TransferFunction_construction:115"]], "(s - 1)*(s + 5)")
  expect_identical(got[["lti:test_TransferFunction_construction:116"]], "((s - 1)*(s + 3), (s - 1)*(s + 5), s)")
  expect_identical(got[["lti:test_TransferFunction_construction:119"]], "(p + 2, p + 2, p)")
  expect_identical(got[["lti:test_TransferFunction_construction:122"]], "(s - 1, 4 - p, s)")
  expect_identical(got[["lti:test_TransferFunction_construction:125"]], "(s - 1, s - 1, s)")
  expect_identical(got[["lti:test_TransferFunction_construction:128"]], "5")
  expect_identical(got[["lti:test_TransferFunction_construction:129"]], "6")
  expect_identical(got[["lti:test_TransferFunction_construction:130"]], "(5, 6, s)")
  expect_identical(got[["lti:test_TransferFunction_construction:133"]], "0.500000000000000")
  expect_identical(got[["lti:test_TransferFunction_construction:134"]], "4")
  expect_identical(got[["lti:test_TransferFunction_construction:135"]], "(0.500000000000000, 4, s)")
  expect_identical(got[["lti:test_TransferFunction_construction:139"]], "TransferFunction(2*p + 3*s^2 + 4*s, 8*p^2 + 7*s, s)")
  expect_identical(got[["lti:test_TransferFunction_construction:143"]], "TransferFunction(a0*s + a1*s^2 + a2*s^3, b0*p - b1*s, s)")
  expect_identical(got[["lti:test_TransferFunction_construction:144"]], "True")
  expect_identical(got[["lti:test_TransferFunction_construction:147"]], "(a*s^3 + b*s^2 + d + g*s, d*p + g*p^2 + g*s, s)")
  expect_identical(got[["lti:test_TransferFunction_construction:151"]], "(d + p^3, a + d*s + g*s^2, p)")
  expect_identical(got[["lti:test_TransferFunction_construction:152"]], "TransferFunction(d + p^3, a + d*s + g*s^2, p)")
  expect_identical(got[["lti:test_TransferFunction_construction:155"]], "a0 + a1*s")
  expect_identical(got[["lti:test_TransferFunction_construction:156"]], "b0 + b1*s + b2*s^2")
  expect_identical(got[["lti:test_TransferFunction_construction:157"]], "(a0 + a1*s, b0 + b1*s + b2*s^2, s)")
  expect_identical(got[["lti:test_TransferFunction_construction:161"]], "(0, p^2 - p + 1, p)")
  expect_identical(got[["lti:test_TransferFunction_construction:164"]], "(0, 1, s)")
  expect_identical(got[["lti:test_TransferFunction_construction:168"]], "(a0*s^0.5 - a1 + a2*s^0.6, a1/p^8.7, s)")
  expect_identical(got[["lti:test_TransferFunction_construction:171"]], "(a1/s^0.8 + a2^2*p^0.25, a0*s - p, p)")
  expect_identical(got[["lti:test_TransferFunction_construction:176"]], "k_i/s + k_o*s + k_p")
  expect_identical(got[["lti:test_TransferFunction_construction:177"]], "(k_i/s + k_o*s + k_p, omega_o^2 + 2*omega_o*s + s^2, s)")
  expect_identical(got[["lti:test_TransferFunction_construction:181"]], "ERR ValueError: TransferFunction cannot have a zero denominator.")
  expect_identical(got[["lti:test_TransferFunction_construction:182"]], "ERR ValueError: TransferFunction cannot have a zero denominator.")
  expect_identical(got[["lti:test_TransferFunction_construction:183"]], "ERR ValueError: TransferFunction cannot have a zero denominator.")
  expect_identical(got[["lti:test_TransferFunction_construction:185"]], "ERR TypeError: Unsupported type for numerator or denominator of\\n                            TransferFunction.")
  expect_identical(got[["lti:test_TransferFunction_construction:187"]], "ERR TypeError: Variable input must be a Symbol.")
  expect_identical(got[["lti:test_TransferFunction_construction:188"]], "ERR TypeError: Variable input must be a Symbol.")
  expect_identical(got[["lti:test_TransferFunction_construction:189"]], "ERR TypeError: Variable input must be a Symbol.")
})

test_that("lti test_TransferFunction_functions matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("a, x, b, c, s, g, d,    p, k, tau, zeta, wn, T, z"); a <- .py_getitem(.tmpu_, 0L); x <- .py_getitem(.tmpu_, 1L); b <- .py_getitem(.tmpu_, 2L); c <- .py_getitem(.tmpu_, 3L); s <- .py_getitem(.tmpu_, 4L); g <- .py_getitem(.tmpu_, 5L); d <- .py_getitem(.tmpu_, 6L); p <- .py_getitem(.tmpu_, 7L); k <- .py_getitem(.tmpu_, 8L); tau <- .py_getitem(.tmpu_, 9L); zeta <- .py_getitem(.tmpu_, 10L); wn <- .py_getitem(.tmpu_, 11L); T <- .py_getitem(.tmpu_, 12L); z <- .py_getitem(.tmpu_, 13L)
})
.py_try({
.tmpu_ <- .py_symbols("a0:4, b0:5, c0:4, d0:4"); a0 <- .py_getitem(.tmpu_, 0L); a1 <- .py_getitem(.tmpu_, 1L); a2 <- .py_getitem(.tmpu_, 2L); a3 <- .py_getitem(.tmpu_, 3L); b0 <- .py_getitem(.tmpu_, 4L); b1 <- .py_getitem(.tmpu_, 5L); b2 <- .py_getitem(.tmpu_, 6L); b3 <- .py_getitem(.tmpu_, 7L); b4 <- .py_getitem(.tmpu_, 8L); c0 <- .py_getitem(.tmpu_, 9L); c1 <- .py_getitem(.tmpu_, 10L); c2 <- .py_getitem(.tmpu_, 11L); c3 <- .py_getitem(.tmpu_, 12L); d0 <- .py_getitem(.tmpu_, 13L); d1 <- .py_getitem(.tmpu_, 14L); d2 <- .py_getitem(.tmpu_, 15L); d3 <- .py_getitem(.tmpu_, 16L)
})
.py_try({
TF1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
TF2 <- TransferFunction(k, 1L, s)
})
.py_try({
TF3 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
TF4 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(tau, s), 1L), s)
})
.py_try({
TF5 <- TransferFunction(k, p, s)
})
.py_try({
TF6 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(k, s), p), s)
})
p_are_floats_equals <- function(n1, n2) {
  "\n    Helper function to compare two floats.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return(.py_method(Float(n1), "is_same", Float(n2), .py_isclose))
}
p_are_complex_equals <- function(n1, n2) {
  "\n    Helper function to compare two complex numbers.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return((.py_truthy(p_are_floats_equals(re(n1), re(n2))) && .py_truthy(p_are_floats_equals(im(n1), im(n2)))))
}
.py_try({
expr_1 <- .py_Mul(0L, .py_Pow(s, -1L, evaluate = FALSE), evaluate = FALSE)
})
.py_try({
expr_2 <- .py_div(s, 0L)
})
.py_try({
expr_3 <- .py_div(.py_add(.py_mul(p, .py_pow(s, 2L)), .py_mul(5L, s)), .py_pow(.py_add(s, 1L), 3L))
})
.py_try({
expr_4 <- 6L
})
.py_try({
expr_5 <- .py_div(.py_mul(.py_add(2L, .py_mul(3L, s)), .py_add(5L, .py_mul(2L, s))), .py_mul(.py_add(9L, .py_mul(3L, s)), .py_add(5L, .py_mul(2L, .py_pow(s, 2L)))))
})
.py_try({
expr_6 <- .py_div(.py_add(.py_add(.py_mul(9L, .py_pow(s, 4L)), .py_mul(4L, .py_pow(s, 2L))), 8L), .py_mul(.py_add(s, 1L), .py_add(s, 9L)))
})
.py_try({
tf <- TransferFunction(.py_add(s, 1L), .py_add(.py_pow(s, 2L), 2L), s)
})
.py_try({
delay <- .py_exp(.py_div(.py_neg(s), tau))
})
.py_try({
expr_7 <- .py_mul(delay, .py_method(tf, "to_expr"))
})
.py_try({
H1 <- .py_method(TransferFunction, "from_rational_expression", expr_7, s)
})
.py_try({
H2 <- TransferFunction(.py_add(s, 1L), .py_mul(.py_add(.py_pow(s, 2L), 2L), .py_exp(.py_div(s, tau))), s)
})
.py_try({
expr_8 <- .py_Add(2L, .py_div(.py_mul(3L, s), .py_add(.py_pow(s, 2L), 1L)), evaluate = FALSE)
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:207", function() .py_method(TransferFunction, "from_rational_expression", expr_1))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:208", function() .py_method(TransferFunction, "from_rational_expression", expr_2))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:209", function() .py_method(TransferFunction, "from_rational_expression", expr_3))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:210", function() .py_method(TransferFunction, "from_rational_expression", expr_3, s))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:211", function() .py_method(TransferFunction, "from_rational_expression", expr_3, p))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:212", function() .py_method(TransferFunction, "from_rational_expression", expr_4))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:213", function() .py_method(TransferFunction, "from_rational_expression", expr_4, s))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:214", function() .py_method(TransferFunction, "from_rational_expression", expr_5, s))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:216", function() .py_method(TransferFunction, "from_rational_expression", expr_6, s))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:218", function() H1)
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:219", function() .py_method(TransferFunction, "from_rational_expression", expr_8, s))
})
.py_try({
tf1 <- .py_method(TransferFunction, "from_coeff_lists", list(1L, 2L), list(3L, 4L, 5L), s)
})
.py_try({
num2 <- list(.py_pow(p, 2L), .py_mul(2L, p))
})
.py_try({
den2 <- list(.py_pow(p, 3L), .py_add(p, 1L), 4L)
})
.py_try({
tf2 <- .py_method(TransferFunction, "from_coeff_lists", num2, den2, s)
})
.py_try({
num3 <- list(1L, 2L, 3L)
})
.py_try({
den3 <- list(0L, 0L)
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:230", function() tf1)
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:231", function() tf2)
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:232", function() .py_method(TransferFunction, "from_coeff_lists", num3, den3, s))
})
.py_try({
zeros <- list(4L)
})
.py_try({
poles <- list(.py_add(-1L, complex(real = 0.0, imaginary = 2.0)), .py_sub(-1L, complex(real = 0.0, imaginary = 2.0)))
})
.py_try({
gain <- 3L
})
.py_try({
tf1 <- .py_method(TransferFunction, "from_zpk", zeros, poles, gain, s)
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:240", function() tf1)
})
.py_try({
tf0 <- TransferFunction(.py_add(.py_add(.py_pow(s, 5L), .py_pow(s, 3L)), s), .py_sub(s, .py_pow(s, 2L)), s)
})
.py_try({
a <- TransferFunction(.py_neg(.py_add(.py_add(.py_pow(s, 4L), .py_pow(s, 2L)), 1L)), .py_sub(s, 1L), s)
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:245", function() (.py_truthy(.py_eq(.py_method(tf0, "simplify"), simplify(tf0))) && .py_truthy(.py_eq(simplify(tf0), a))))
})
.py_try({
tf1 <- TransferFunction(.py_mul(.py_add(p, 3L), .py_sub(p, 1L)), .py_mul(.py_sub(p, 1L), .py_add(p, 5L)), p)
})
.py_try({
b <- TransferFunction(.py_add(p, 3L), .py_add(p, 5L), p)
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:249", function() (.py_truthy(.py_eq(.py_method(tf1, "simplify"), simplify(tf1))) && .py_truthy(.py_eq(simplify(tf1), b))))
})
.py_try({
G1 <- TransferFunction(.py_pow(.py_sub(1L, s), 2L), .py_pow(.py_add(.py_pow(s, 2L), 1L), 2L), s)
})
.py_try({
G2 <- TransferFunction(1L, -3L, p)
})
.py_try({
c <- .py_mul(.py_add(.py_add(.py_mul(a2, .py_pow(s, p)), .py_mul(a1, .py_pow(s, s))), .py_mul(a0, .py_pow(p, p))), .py_add(.py_pow(p, s), .py_pow(s, p)))
})
.py_try({
d <- .py_mul(.py_add(.py_mul(b0, .py_pow(s, s)), .py_mul(b1, .py_pow(p, s))), .py_add(.py_mul(.py_mul(b2, s), p), .py_pow(p, p)))
})
.py_try({
e <- .py_add(.py_add(.py_add(.py_add(.py_add(.py_mul(.py_mul(a0, .py_pow(p, p)), .py_pow(p, s)), .py_mul(.py_mul(a0, .py_pow(p, p)), .py_pow(s, p))), .py_mul(.py_mul(a1, .py_pow(p, s)), .py_pow(s, s))), .py_mul(.py_mul(a1, .py_pow(s, p)), .py_pow(s, s))), .py_mul(.py_mul(a2, .py_pow(p, s)), .py_pow(s, p))), .py_mul(a2, .py_pow(s, .py_mul(2L, p))))
})
.py_try({
f <- .py_add(.py_add(.py_add(.py_mul(.py_mul(.py_mul(.py_mul(b0, b2), p), s), .py_pow(s, s)), .py_mul(.py_mul(b0, .py_pow(p, p)), .py_pow(s, s))), .py_mul(.py_mul(.py_mul(.py_mul(b1, b2), p), .py_pow(p, s)), s)), .py_mul(.py_mul(b1, .py_pow(p, p)), .py_pow(p, s)))
})
.py_try({
g <- .py_add(.py_add(.py_add(.py_mul(.py_mul(.py_mul(a1, a2), s), .py_pow(s, p)), .py_mul(.py_mul(a1, p), s)), .py_mul(.py_mul(.py_mul(.py_mul(a2, b1), p), s), .py_pow(s, p))), .py_mul(.py_mul(b1, .py_pow(p, 2L)), s))
})
.py_try({
G3 <- TransferFunction(c, d, s)
})
.py_try({
G4 <- TransferFunction(.py_sub(.py_mul(a0, .py_pow(s, s)), .py_mul(b0, .py_pow(p, p))), .py_mul(.py_add(.py_mul(a1, s), .py_mul(.py_mul(b1, s), p)), .py_add(.py_mul(a2, .py_pow(s, p)), p)), p)
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:262", function() .py_method(G1, "expand"))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:263", function() .py_method(tf1, "expand"))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:264", function() .py_method(G2, "expand"))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:265", function() .py_method(G3, "expand"))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:266", function() .py_method(G4, "expand"))
})
.py_try({
p1 <- .py_add(.py_mul(a1, s), a0)
})
.py_try({
p2 <- .py_add(.py_add(.py_mul(b2, .py_pow(s, 2L)), .py_mul(b1, s)), b0)
})
.py_try({
SP1 <- TransferFunction(p1, p2, s)
})
.py_try({
expect1 <- TransferFunction(.py_add(.py_mul(.py_flt(2.0), s), .py_flt(1.0)), .py_add(.py_add(.py_mul(.py_flt(5.0), .py_pow(s, 2L)), .py_mul(.py_flt(4.0), s)), .py_flt(3.0)), s)
})
.py_try({
expect1_ <- TransferFunction(.py_add(.py_mul(2L, s), 1L), .py_add(.py_add(.py_mul(5L, .py_pow(s, 2L)), .py_mul(4L, s)), 3L), s)
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:274", function() .py_method(SP1, "subs", .py_dict(list(a0, a1, b0, b1, b2), list(1L, 2L, 3L, 4L, 5L))))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:275", function() .py_method(.py_method(SP1, "subs", .py_dict(list(a0, a1, b0, b1, b2), list(1L, 2L, 3L, 4L, 5L))), "evalf"))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:276", function() .py_method(expect1_, "evalf"))
})
.py_try({
.tmpu_ <- .py_symbols("c1, d0:3"); c1 <- .py_getitem(.tmpu_, 0L); d0 <- .py_getitem(.tmpu_, 1L); d1 <- .py_getitem(.tmpu_, 2L); d2 <- .py_getitem(.tmpu_, 3L)
})
.py_try({
.tmpu_ <- .py_tup(.py_mul(c1, p), .py_sub(.py_add(.py_mul(d2, .py_pow(p, 3L)), .py_mul(d1, .py_pow(p, 2L))), d0)); p3 <- .py_getitem(.tmpu_, 0L); p4 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
SP2 <- TransferFunction(p3, p4, p)
})
.py_try({
expect2 <- TransferFunction(.py_mul(.py_flt(2.0), p), .py_sub(.py_add(.py_mul(.py_flt(5.0), .py_pow(p, 3L)), .py_mul(.py_flt(2.0), .py_pow(p, 2L))), .py_flt(3.0)), p)
})
.py_try({
expect2_ <- TransferFunction(.py_mul(2L, p), .py_sub(.py_add(.py_mul(5L, .py_pow(p, 3L)), .py_mul(2L, .py_pow(p, 2L))), 3L), p)
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:283", function() .py_method(SP2, "subs", .py_dict(list(c1, d0, d1, d2), list(2L, 3L, 2L, 5L))))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:284", function() .py_method(.py_method(SP2, "subs", .py_dict(list(c1, d0, d1, d2), list(2L, 3L, 2L, 5L))), "evalf"))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:285", function() .py_method(expect2_, "evalf"))
})
.py_try({
SP3 <- TransferFunction(.py_add(.py_sub(.py_add(.py_mul(a0, .py_pow(p, 3L)), .py_mul(a1, .py_pow(s, 2L))), .py_mul(b0, s)), b1), .py_add(.py_mul(a1, s), p), s)
})
.py_try({
expect3 <- TransferFunction(.py_add(.py_sub(.py_add(.py_mul(.py_flt(2.0), .py_pow(p, 3L)), .py_mul(.py_flt(4.0), .py_pow(s, 2L))), s), .py_flt(5.0)), .py_add(p, .py_mul(.py_flt(4.0), s)), s)
})
.py_try({
expect3_ <- TransferFunction(.py_add(.py_sub(.py_add(.py_mul(2L, .py_pow(p, 3L)), .py_mul(4L, .py_pow(s, 2L))), s), 5L), .py_add(p, .py_mul(4L, s)), s)
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:290", function() .py_method(SP3, "subs", .py_dict(list(a0, a1, b0, b1), list(2L, 4L, 1L, 5L))))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:291", function() .py_method(.py_method(SP3, "subs", .py_dict(list(a0, a1, b0, b1), list(2L, 4L, 1L, 5L))), "evalf"))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:292", function() .py_method(expect3_, "evalf"))
})
.py_try({
SP4 <- TransferFunction(.py_sub(s, .py_mul(a1, .py_pow(p, 3L))), .py_add(.py_mul(a0, s), p), p)
})
.py_try({
expect4 <- TransferFunction(.py_add(.py_mul(.py_flt(7.0), .py_pow(p, 3L)), s), .py_sub(p, s), p)
})
.py_try({
expect4_ <- TransferFunction(.py_add(.py_mul(7L, .py_pow(p, 3L)), s), .py_sub(p, s), p)
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:297", function() .py_method(SP4, "subs", .py_dict(list(a0, a1), list(-1L, -7L))))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:298", function() .py_method(.py_method(SP4, "subs", .py_dict(list(a0, a1), list(-1L, -7L))), "evalf"))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:299", function() .py_method(expect4_, "evalf"))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:302", function() simplify(.py_sub(.py_method(tf1, "eval_frequency", wn), .py_sub(.py_add(.py_div(.py_pow(wn, 2L), .py_sub(.py_add(.py_pow(wn, 2L), .py_mul(4L, wn)), 5L)), .py_div(.py_mul(2L, wn), .py_sub(.py_add(.py_pow(wn, 2L), .py_mul(4L, wn)), 5L))), .py_div(3L, .py_sub(.py_add(.py_pow(wn, 2L), .py_mul(4L, wn)), 5L))))))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:305", function() .py_method(G1, "eval_frequency", .py_add(1L, sym_I)))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:306", function() .py_method(G4, "eval_frequency", .py_div(.py_S(5L), 3L)))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:310", function() .py_method(tf0, "dc_gain"))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:311", function() .py_method(tf1, "dc_gain"))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:312", function() .py_method(SP2, "dc_gain"))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:313", function() .py_method(expect4, "dc_gain"))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:314", function() .py_method(expect2_, "dc_gain"))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:315", function() .py_method(TransferFunction(1L, s, s), "dc_gain"))
})
.py_try({
tf_ <- TransferFunction(.py_sub(.py_pow(x, 3L), k), k, x)
})
.py_try({
p_tf <- TransferFunction(k, .py_sub(.py_pow(x, 4L), k), x)
})
.py_try({
TF_ <- TransferFunction(.py_pow(x, 2L), .py_add(.py_add(.py_pow(x, 10L), x), .py_pow(x, 2L)), x)
})
.py_try({
p_TF <- TransferFunction(.py_add(.py_add(.py_pow(x, 10L), x), .py_pow(x, 2L)), .py_pow(x, 2L), x)
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:322", function() .py_method(G1, "poles"))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:323", function() .py_method(G2, "poles"))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:324", function() .py_method(tf1, "poles"))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:325", function() .py_method(expect4_, "poles"))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:326", function() .py_method(SP4, "poles"))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:327", function() .py_method(expect3, "poles"))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:329", function() .py_str(.py_method(expect1, "poles")))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:330", function() .py_method(p_tf, "poles"))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:331", function() .py_method(TF_, "poles"))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:334", function() .py_method(TransferFunction(.py_pow(x, 2L), .py_add(.py_add(.py_mul(a0, .py_pow(x, 10L)), x), .py_pow(x, 2L)), x), "poles"))
})
.py_try({
.tmpu_ <- .py_symbols("q, r", negative = TRUE); q <- .py_getitem(.tmpu_, 0L); r <- .py_getitem(.tmpu_, 1L)
})
.py_try({
t <- .py_symbols("t", positive = TRUE)
})
.py_try({
TF_ <- TransferFunction(.py_sub(.py_add(.py_pow(s, 2L), a0), .py_mul(a1, p)), .py_sub(.py_mul(q, s), r), s)
})
.py_try({
stable_tf <- TransferFunction(.py_sub(.py_add(.py_pow(s, 2L), a0), .py_mul(a1, p)), .py_sub(.py_mul(q, s), 1L), s)
})
.py_try({
stable_tf_ <- TransferFunction(.py_sub(.py_add(.py_pow(s, 2L), a0), .py_mul(a1, p)), .py_sub(.py_mul(q, s), t), s)
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:343", function() .py_method(G1, "is_stable"))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:344", function() .py_method(G2, "is_stable"))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:345", function() .py_method(tf1, "is_stable"))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:346", function() .py_method(tf1, "is_stable", TRUE))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:347", function() .py_method(expect2, "is_stable"))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:348", function() .py_method(expect1, "is_stable"))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:349", function() .py_method(stable_tf, "is_stable"))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:350", function() .py_method(stable_tf_, "is_stable"))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:351", function() .py_method(TF_, "is_stable"))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:352", function() .py_method(expect4_, "is_stable"))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:353", function() .py_method(SP4, "is_stable"))
})
.py_try({
generic_den <- .py_add(.py_add(.py_add(.py_add(.py_mul(b4, .py_pow(s, 4L)), .py_mul(b3, .py_pow(s, 3L))), .py_mul(b2, .py_pow(s, 2L))), .py_mul(b1, s)), b0)
})
.py_try({
stab_cond <- .py_method(TransferFunction(1L, generic_den, s), "get_asymptotic_stability_conditions")
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:358", function() stab_cond)
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:362", function() .py_method(TransferFunction(1L, .py_mul(.py_mul(.py_add(s, 1L), .py_add(s, .py_mul(2L, sym_I))), .py_sub(s, .py_mul(2L, sym_I))), s), "get_asymptotic_stability_conditions"))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:363", function() .py_method(TransferFunction(1L, .py_mul(.py_mul(.py_add(s, 1L), .py_add(s, 2L)), .py_add(s, .py_div(1L, 2L))), s), "get_asymptotic_stability_conditions"))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:364", function() .py_method(stable_tf, "get_asymptotic_stability_conditions"))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:367", function() .py_method(G1, "zeros"))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:368", function() .py_method(G2, "zeros"))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:369", function() .py_method(tf1, "zeros"))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:370", function() .py_method(expect4_, "zeros"))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:375", function() .py_method(SP4, "zeros"))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:382", function() .py_method(tf_, "zeros"))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:384", function() .py_method(p_TF, "zeros"))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:387", function() .py_method(TransferFunction(.py_add(.py_add(.py_mul(a0, .py_pow(x, 10L)), x), .py_pow(x, 2L)), .py_pow(x, 2L), x), "zeros"))
})
.py_try({
tf2 <- TransferFunction(.py_add(s, 3L), .py_add(.py_sub(.py_pow(s, 2L), .py_pow(s, 3L)), 9L), s)
})
.py_try({
tf3 <- TransferFunction(.py_add(.py_mul(-3L, p), 3L), .py_sub(1L, p), p)
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:392", function() .py_neg(tf2))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:393", function() .py_neg(tf3))
})
.py_try({
tf4 <- TransferFunction(.py_add(p, 4L), .py_sub(p, 3L), p)
})
.py_try({
tf5 <- TransferFunction(.py_add(.py_pow(s, 2L), 1L), .py_sub(1L, s), s)
})
.py_try({
expect2 <- TransferFunction(.py_pow(.py_add(.py_pow(s, 2L), 1L), 3L), .py_pow(.py_sub(1L, s), 3L), s)
})
.py_try({
expect1 <- TransferFunction(.py_pow(.py_add(p, 4L), 2L), .py_pow(.py_sub(p, 3L), 2L), p)
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:400", function() (.py_truthy(.py_eq(.py_method(.py_mul(tf4, tf4), "doit"), .py_pow(tf4, 2L))) && .py_truthy(.py_eq(.py_pow(tf4, 2L), pow(tf4, 2L))) && .py_truthy(.py_eq(pow(tf4, 2L), expect1))))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:401", function() (.py_truthy(.py_eq(.py_method(.py_mul(.py_mul(tf5, tf5), tf5), "doit"), .py_pow(tf5, 3L))) && .py_truthy(.py_eq(.py_pow(tf5, 3L), pow(tf5, 3L))) && .py_truthy(.py_eq(pow(tf5, 3L), expect2))))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:402", function() (.py_truthy(.py_eq(.py_pow(tf5, 0L), pow(tf5, 0L))) && .py_truthy(.py_eq(pow(tf5, 0L), TransferFunction(1L, 1L, s)))))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:403", function() (.py_truthy(.py_eq(.py_pow(.py_method(Series(tf4), "doit"), -1L), .py_pow(tf4, -1L))) && .py_truthy(.py_eq(.py_pow(tf4, -1L), pow(tf4, -1L))) && .py_truthy(.py_eq(pow(tf4, -1L), TransferFunction(.py_sub(p, 3L), .py_add(p, 4L), p)))))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:404", function() (.py_truthy(.py_eq(.py_pow(.py_method(.py_mul(tf5, tf5), "doit"), -1L), .py_pow(tf5, -2L))) && .py_truthy(.py_eq(.py_pow(tf5, -2L), pow(tf5, -2L))) && .py_truthy(.py_eq(pow(tf5, -2L), TransferFunction(.py_pow(.py_sub(1L, s), 2L), .py_pow(.py_add(.py_pow(s, 2L), 1L), 2L), s)))))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:406", function() .py_pow(tf4, .py_sub(.py_add(.py_pow(s, 2L), s), 1L)))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:407", function() .py_pow(tf5, s))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:408", function() .py_pow(tf4, tf5))
})
.py_try({
tf <- TransferFunction(.py_sub(s, 1L), .py_add(.py_sub(.py_pow(s, 2L), .py_mul(2L, s)), 1L), s)
})
.py_try({
tf6 <- TransferFunction(.py_add(s, p), .py_sub(.py_pow(p, 2L), 5L), s)
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:413", function() sym_factor(tf))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:414", function() (.py_truthy(.py_eq(.py_method(.py_attr(tf, "num"), "subs", s, 2L), .py_method(.py_attr(tf, "den"), "subs", s, 2L))) && .py_truthy(.py_eq(.py_method(.py_attr(tf, "den"), "subs", s, 2L), 1L))))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:416", function() .py_method(tf, "subs", s, 2L))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:417", function() .py_method(tf6, "subs", p, 3L))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:418", function() .py_method(tf3, "xreplace", .py_dict(list(p), list(s))))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:419", function() .py_method(tf3, "xreplace", .py_dict(list(p), list(.py_exp(2L)))))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:420", function() .py_method(tf3, "subs", p, .py_exp(2L)))
})
.py_try({
tf7 <- TransferFunction(.py_add(.py_mul(a0, .py_pow(s, p)), .py_mul(a1, .py_pow(p, s))), .py_sub(.py_mul(a2, p), s), s)
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:423", function() .py_method(tf7, "xreplace", .py_dict(list(s), list(k))))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:424", function() .py_method(tf7, "subs", s, k))
})
.py_try({
tf8 <- TransferFunction(.py_add(.py_add(.py_mul(a0, .py_pow(s, 5L)), .py_mul(5L, .py_pow(s, 2L))), 3L), .py_sub(.py_pow(s, 6L), 3L), s)
})
.py_try({
tf9 <- TransferFunction(.py_add(5L, s), .py_mul(.py_add(5L, s), .py_add(6L, s)), s)
})
.py_try({
tf10 <- TransferFunction(0L, 1L, s)
})
.py_try({
tf11 <- TransferFunction(1L, 1L, s)
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:431", function() .py_method(tf8, "to_expr"))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:432", function() .py_method(tf9, "to_expr"))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:433", function() .py_method(tf10, "to_expr"))
})
.py_try({
.py_rec("lti:test_TransferFunction_functions:434", function() .py_method(tf11, "to_expr"))
}))---")
  expect_identical(got[["lti:test_TransferFunction_functions:207"]], "TransferFunction(0, s, s)")
  expect_identical(got[["lti:test_TransferFunction_functions:208"]], "ERR ZeroDivisionError: TransferFunction cannot have a zero denominator.")
  expect_identical(got[["lti:test_TransferFunction_functions:209"]], "ERR ValueError: \\nConflicting values found for positional argument `var` ({p, s}).\\nSpecify it manually.")
  expect_identical(got[["lti:test_TransferFunction_functions:210"]], "TransferFunction(p*s^2 + 5*s, (s + 1)^3, s)")
  expect_identical(got[["lti:test_TransferFunction_functions:211"]], "TransferFunction(p*s^2 + 5*s, (s + 1)^3, p)")
  expect_identical(got[["lti:test_TransferFunction_functions:212"]], "ERR ValueError: \\nPositional argument `var` not found in the TransferFunction defined.\\nSpecify it manually.")
  expect_identical(got[["lti:test_TransferFunction_functions:213"]], "TransferFunction(6, 1, s)")
  expect_identical(got[["lti:test_TransferFunction_functions:214"]], "TransferFunction((2*s + 5)*(3*s + 2), (3*s + 9)*(2*s^2 + 5), s)")
  expect_identical(got[["lti:test_TransferFunction_functions:216"]], "TransferFunction(9*s^4 + 4*s^2 + 8, (s + 1)*(s + 9), s)")
  expect_identical(got[["lti:test_TransferFunction_functions:218"]], "TransferFunction(s + 1, (s^2 + 2)*exp(s/tau), s)")
  expect_identical(got[["lti:test_TransferFunction_functions:219"]], "TransferFunction(2*s^2 + 3*s + 2, s^2 + 1, s)")
  expect_identical(got[["lti:test_TransferFunction_functions:230"]], "TransferFunction(s + 2, 3*s^2 + 4*s + 5, s)")
  expect_identical(got[["lti:test_TransferFunction_functions:231"]], "TransferFunction(p^2*s + 2*p, p^3*s^2 + s*(p + 1) + 4, s)")
  expect_identical(got[["lti:test_TransferFunction_functions:232"]], "ERR ZeroDivisionError: TransferFunction cannot have a zero denominator.")
  expect_identical(got[["lti:test_TransferFunction_functions:240"]], "TransferFunction(3*s - 12, (s + 1.0 - 2.0*I)*(s + 1.0 + 2.0*I), s)")
  expect_identical(got[["lti:test_TransferFunction_functions:245"]], "True")
  expect_identical(got[["lti:test_TransferFunction_functions:249"]], "True")
  expect_identical(got[["lti:test_TransferFunction_functions:262"]], "TransferFunction(s^2 - 2*s + 1, s^4 + 2*s^2 + 1, s)")
  expect_identical(got[["lti:test_TransferFunction_functions:263"]], "TransferFunction(p^2 + 2*p - 3, p^2 + 4*p - 5, p)")
  expect_identical(got[["lti:test_TransferFunction_functions:264"]], "TransferFunction(1, -3, p)")
  expect_identical(got[["lti:test_TransferFunction_functions:265"]], "TransferFunction(a0*p^p*p^s + a0*p^p*s^p + a1*p^s*s^s + a1*s^p*s^s + a2*p^s*s^p + a2*s^(2*p), b0*b2*p*s*s^s + b0*p^p*s^s + b1*b2*p*p^s*s + b1*p^p*p^s, s)")
  expect_identical(got[["lti:test_TransferFunction_functions:266"]], "TransferFunction(a0*s^s - b0*p^p, a1*a2*s*s^p + a1*p*s + a2*b1*p*s*s^p + b1*p^2*s, p)")
  expect_identical(got[["lti:test_TransferFunction_functions:274"]], "TransferFunction(2*s + 1, 5*s^2 + 4*s + 3, s)")
  expect_identical(got[["lti:test_TransferFunction_functions:275"]], "TransferFunction(2.0*s + 1.0, 5.0*s^2 + 4.0*s + 3.0, s)")
  expect_identical(got[["lti:test_TransferFunction_functions:276"]], "TransferFunction(2.0*s + 1.0, 5.0*s^2 + 4.0*s + 3.0, s)")
  expect_identical(got[["lti:test_TransferFunction_functions:283"]], "TransferFunction(2*p, 5*p^3 + 2*p^2 - 3, p)")
  expect_identical(got[["lti:test_TransferFunction_functions:284"]], "TransferFunction(2.0*p, 5.0*p^3 + 2.0*p^2 - 3.0, p)")
  expect_identical(got[["lti:test_TransferFunction_functions:285"]], "TransferFunction(2.0*p, 5.0*p^3 + 2.0*p^2 - 3.0, p)")
  expect_identical(got[["lti:test_TransferFunction_functions:290"]], "TransferFunction(2*p^3 + 4*s^2 - s + 5, p + 4*s, s)")
  expect_identical(got[["lti:test_TransferFunction_functions:291"]], "TransferFunction(2.0*p^3 + 4.0*s^2 - s + 5.0, p + 4.0*s, s)")
  expect_identical(got[["lti:test_TransferFunction_functions:292"]], "TransferFunction(2.0*p^3 + 4.0*s^2 - s + 5.0, p + 4.0*s, s)")
  expect_identical(got[["lti:test_TransferFunction_functions:297"]], "TransferFunction(7*p^3 + s, p - s, p)")
  expect_identical(got[["lti:test_TransferFunction_functions:298"]], "TransferFunction(7.0*p^3 + s, p - s, p)")
  expect_identical(got[["lti:test_TransferFunction_functions:299"]], "TransferFunction(7.0*p^3 + s, p - s, p)")
  expect_identical(got[["lti:test_TransferFunction_functions:302"]], "0")
  expect_identical(got[["lti:test_TransferFunction_functions:305"]], "3/25 + 4*I/25")
  expect_identical(got[["lti:test_TransferFunction_functions:306"]], "a0*s^s/(a1*a2*s^(8/3) + 5*a1*s/3 + 5*a2*b1*s^(8/3)/3 + 25*b1*s/9) - 5*3^(1/3)*5^(2/3)*b0/(9*a1*a2*s^(8/3) + 15*a1*s + 15*a2*b1*s^(8/3) + 25*b1*s)")
  expect_identical(got[["lti:test_TransferFunction_functions:310"]], "1")
  expect_identical(got[["lti:test_TransferFunction_functions:311"]], "3/5")
  expect_identical(got[["lti:test_TransferFunction_functions:312"]], "0")
  expect_identical(got[["lti:test_TransferFunction_functions:313"]], "-1")
  expect_identical(got[["lti:test_TransferFunction_functions:314"]], "0")
  expect_identical(got[["lti:test_TransferFunction_functions:315"]], "oo")
  expect_identical(got[["lti:test_TransferFunction_functions:322"]], "[I, I, -I, -I]")
  expect_identical(got[["lti:test_TransferFunction_functions:323"]], "[]")
  expect_identical(got[["lti:test_TransferFunction_functions:324"]], "[-5, 1]")
  expect_identical(got[["lti:test_TransferFunction_functions:325"]], "[s]")
  expect_identical(got[["lti:test_TransferFunction_functions:326"]], "[-a0*s]")
  expect_identical(got[["lti:test_TransferFunction_functions:327"]], "[-0.25*p]")
  expect_identical(got[["lti:test_TransferFunction_functions:329"]], "'[-0.4 - 0.66332495807108*I, -0.4 + 0.66332495807108*I]'")
  expect_identical(got[["lti:test_TransferFunction_functions:330"]], "[k^(1/4), -k^(1/4), I*k^(1/4), -I*k^(1/4)]")
  expect_identical(got[["lti:test_TransferFunction_functions:331"]], "[CRootOf(x^9 + x + 1, 0), 0, CRootOf(x^9 + x + 1, 1), CRootOf(x^9 + x + 1, 2), CRootOf(x^9 + x + 1, 3), CRootOf(x^9 + x + 1, 4), CRootOf(x^9 + x + 1, 5), CRootOf(x^9 + x + 1, 6), CRootOf(x^9 + x + 1, 7), CRootOf(x^9 + x + 1, 8)]")
  expect_identical(got[["lti:test_TransferFunction_functions:334"]], "ERR NotImplementedError: CRootOf is not supported over ZZ[a0]")
  expect_identical(got[["lti:test_TransferFunction_functions:343"]], "False")
  expect_identical(got[["lti:test_TransferFunction_functions:344"]], "True")
  expect_identical(got[["lti:test_TransferFunction_functions:345"]], "False")
  expect_identical(got[["lti:test_TransferFunction_functions:346"]], "True")
  expect_identical(got[["lti:test_TransferFunction_functions:347"]], "False")
  expect_identical(got[["lti:test_TransferFunction_functions:348"]], "True")
  expect_identical(got[["lti:test_TransferFunction_functions:349"]], "True")
  expect_identical(got[["lti:test_TransferFunction_functions:350"]], "True")
  expect_identical(got[["lti:test_TransferFunction_functions:351"]], "False")
  expect_identical(got[["lti:test_TransferFunction_functions:352"]], "None")
  expect_identical(got[["lti:test_TransferFunction_functions:353"]], "None")
  expect_identical(got[["lti:test_TransferFunction_functions:358"]], "[b3*b4 > 0, -b1*b4 + b2*b3 > 0, -b0*b3^2*b4 - b1^2*b4^2 + b1*b2*b3*b4 > 0, b0*b4 > 0]")
  expect_identical(got[["lti:test_TransferFunction_functions:362"]], "[False]")
  expect_identical(got[["lti:test_TransferFunction_functions:363"]], "[True, True, True]")
  expect_identical(got[["lti:test_TransferFunction_functions:364"]], "[True]")
  expect_identical(got[["lti:test_TransferFunction_functions:367"]], "[1, 1]")
  expect_identical(got[["lti:test_TransferFunction_functions:368"]], "[]")
  expect_identical(got[["lti:test_TransferFunction_functions:369"]], "[-3, 1]")
  expect_identical(got[["lti:test_TransferFunction_functions:370"]], "[-7^(2/3)*s^(1/3)/7, 7^(2/3)*s^(1/3)/14 - sqrt(3)*7^(2/3)*I*s^(1/3)/14, 7^(2/3)*s^(1/3)/14 + sqrt(3)*7^(2/3)*I*s^(1/3)/14]")
  expect_identical(got[["lti:test_TransferFunction_functions:375"]], "[s^(1/3)/a1^(1/3), -s^(1/3)/(2*a1^(1/3)) - sqrt(3)*I*s^(1/3)/(2*a1^(1/3)), -s^(1/3)/(2*a1^(1/3)) + sqrt(3)*I*s^(1/3)/(2*a1^(1/3))]")
  expect_identical(got[["lti:test_TransferFunction_functions:382"]], "[k^(1/3), -k^(1/3)/2 - sqrt(3)*I*k^(1/3)/2, -k^(1/3)/2 + sqrt(3)*I*k^(1/3)/2]")
  expect_identical(got[["lti:test_TransferFunction_functions:384"]], "[CRootOf(x^9 + x + 1, 0), 0, CRootOf(x^9 + x + 1, 1), CRootOf(x^9 + x + 1, 2), CRootOf(x^9 + x + 1, 3), CRootOf(x^9 + x + 1, 4), CRootOf(x^9 + x + 1, 5), CRootOf(x^9 + x + 1, 6), CRootOf(x^9 + x + 1, 7), CRootOf(x^9 + x + 1, 8)]")
  expect_identical(got[["lti:test_TransferFunction_functions:387"]], "ERR NotImplementedError: CRootOf is not supported over ZZ[a0]")
  expect_identical(got[["lti:test_TransferFunction_functions:392"]], "TransferFunction(-s - 3, -s^3 + s^2 + 9, s)")
  expect_identical(got[["lti:test_TransferFunction_functions:393"]], "TransferFunction(3*p - 3, 1 - p, p)")
  expect_identical(got[["lti:test_TransferFunction_functions:400"]], "True")
  expect_identical(got[["lti:test_TransferFunction_functions:401"]], "True")
  expect_identical(got[["lti:test_TransferFunction_functions:402"]], "True")
  expect_identical(got[["lti:test_TransferFunction_functions:403"]], "True")
  expect_identical(got[["lti:test_TransferFunction_functions:404"]], "True")
  expect_identical(got[["lti:test_TransferFunction_functions:406"]], "ERR ValueError: Exponent must be an integer.")
  expect_identical(got[["lti:test_TransferFunction_functions:407"]], "ERR ValueError: Exponent must be an integer.")
  expect_identical(got[["lti:test_TransferFunction_functions:408"]], "ERR ValueError: Exponent must be an integer.")
  expect_identical(got[["lti:test_TransferFunction_functions:413"]], "TransferFunction(s - 1, (s - 1)^2, s)")
  expect_identical(got[["lti:test_TransferFunction_functions:414"]], "True")
  expect_identical(got[["lti:test_TransferFunction_functions:416"]], "TransferFunction(s - 1, s^2 - 2*s + 1, s)")
  expect_identical(got[["lti:test_TransferFunction_functions:417"]], "TransferFunction(s + 3, 4, s)")
  expect_identical(got[["lti:test_TransferFunction_functions:418"]], "TransferFunction(3 - 3*s, 1 - s, s)")
  expect_identical(got[["lti:test_TransferFunction_functions:419"]], "ERR TypeError: Variable input must be a Symbol.")
  expect_identical(got[["lti:test_TransferFunction_functions:420"]], "TransferFunction(3 - 3*p, 1 - p, p)")
  expect_identical(got[["lti:test_TransferFunction_functions:423"]], "TransferFunction(a0*k^p + a1*p^k, a2*p - k, k)")
  expect_identical(got[["lti:test_TransferFunction_functions:424"]], "TransferFunction(a0*s^p + a1*p^s, a2*p - s, s)")
  expect_identical(got[["lti:test_TransferFunction_functions:431"]], "(a0*s^5 + 5*s^2 + 3)/(s^6 - 3)")
  expect_identical(got[["lti:test_TransferFunction_functions:432"]], "(s + 5)/(((s + 5)*(s + 6)))")
  expect_identical(got[["lti:test_TransferFunction_functions:433"]], "0/1")
  expect_identical(got[["lti:test_TransferFunction_functions:434"]], "1/1")
})

test_that("lti test_TransferFunction_addition_and_subtraction matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("a, x, b, c, s, g, d,    p, k, tau, zeta, wn, T, z"); a <- .py_getitem(.tmpu_, 0L); x <- .py_getitem(.tmpu_, 1L); b <- .py_getitem(.tmpu_, 2L); c <- .py_getitem(.tmpu_, 3L); s <- .py_getitem(.tmpu_, 4L); g <- .py_getitem(.tmpu_, 5L); d <- .py_getitem(.tmpu_, 6L); p <- .py_getitem(.tmpu_, 7L); k <- .py_getitem(.tmpu_, 8L); tau <- .py_getitem(.tmpu_, 9L); zeta <- .py_getitem(.tmpu_, 10L); wn <- .py_getitem(.tmpu_, 11L); T <- .py_getitem(.tmpu_, 12L); z <- .py_getitem(.tmpu_, 13L)
})
.py_try({
.tmpu_ <- .py_symbols("a0:4, b0:5, c0:4, d0:4"); a0 <- .py_getitem(.tmpu_, 0L); a1 <- .py_getitem(.tmpu_, 1L); a2 <- .py_getitem(.tmpu_, 2L); a3 <- .py_getitem(.tmpu_, 3L); b0 <- .py_getitem(.tmpu_, 4L); b1 <- .py_getitem(.tmpu_, 5L); b2 <- .py_getitem(.tmpu_, 6L); b3 <- .py_getitem(.tmpu_, 7L); b4 <- .py_getitem(.tmpu_, 8L); c0 <- .py_getitem(.tmpu_, 9L); c1 <- .py_getitem(.tmpu_, 10L); c2 <- .py_getitem(.tmpu_, 11L); c3 <- .py_getitem(.tmpu_, 12L); d0 <- .py_getitem(.tmpu_, 13L); d1 <- .py_getitem(.tmpu_, 14L); d2 <- .py_getitem(.tmpu_, 15L); d3 <- .py_getitem(.tmpu_, 16L)
})
.py_try({
TF1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
TF2 <- TransferFunction(k, 1L, s)
})
.py_try({
TF3 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
TF4 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(tau, s), 1L), s)
})
.py_try({
TF5 <- TransferFunction(k, p, s)
})
.py_try({
TF6 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(k, s), p), s)
})
p_are_floats_equals <- function(n1, n2) {
  "\n    Helper function to compare two floats.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return(.py_method(Float(n1), "is_same", Float(n2), .py_isclose))
}
p_are_complex_equals <- function(n1, n2) {
  "\n    Helper function to compare two complex numbers.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return((.py_truthy(p_are_floats_equals(re(n1), re(n2))) && .py_truthy(p_are_floats_equals(im(n1), im(n2)))))
}
.py_try({
tf1 <- TransferFunction(.py_add(s, 6L), .py_sub(s, 5L), s)
})
.py_try({
tf2 <- TransferFunction(.py_add(s, 3L), .py_add(s, 1L), s)
})
.py_try({
tf3 <- TransferFunction(.py_add(s, 1L), .py_add(.py_add(.py_pow(s, 2L), s), 1L), s)
})
.py_try({
tf4 <- TransferFunction(p, .py_sub(2L, p), p)
})
.py_try({
.py_rec("lti:test_TransferFunction_addition_and_subtraction:444", function() .py_add(tf1, tf2))
})
.py_try({
.py_rec("lti:test_TransferFunction_addition_and_subtraction:445", function() .py_add(tf3, tf1))
})
.py_try({
.py_rec("lti:test_TransferFunction_addition_and_subtraction:446", function() .py_add(.py_add(.py_neg(tf1), tf2), tf3))
})
.py_try({
.py_rec("lti:test_TransferFunction_addition_and_subtraction:447", function() .py_add(tf1, .py_add(tf2, tf3)))
})
.py_try({
c <- .py_symbols("c", commutative = FALSE)
})
.py_try({
.py_rec("lti:test_TransferFunction_addition_and_subtraction:450", function() .py_add(tf1, .py_Matrix(list(1L, 2L, 3L))))
})
.py_try({
.py_rec("lti:test_TransferFunction_addition_and_subtraction:451", function() .py_add(tf2, c))
})
.py_try({
.py_rec("lti:test_TransferFunction_addition_and_subtraction:452", function() .py_add(tf3, tf4))
})
.py_try({
.py_rec("lti:test_TransferFunction_addition_and_subtraction:453", function() .py_add(tf1, .py_sub(s, 1L)))
})
.py_try({
.py_rec("lti:test_TransferFunction_addition_and_subtraction:454", function() .py_add(tf1, 8L))
})
.py_try({
.py_rec("lti:test_TransferFunction_addition_and_subtraction:455", function() .py_add(.py_sub(1L, .py_pow(p, 3L)), tf1))
})
.py_try({
.py_rec("lti:test_TransferFunction_addition_and_subtraction:458", function() .py_sub(tf1, tf2))
})
.py_try({
.py_rec("lti:test_TransferFunction_addition_and_subtraction:459", function() .py_sub(tf3, tf2))
})
.py_try({
.py_rec("lti:test_TransferFunction_addition_and_subtraction:460", function() .py_sub(.py_neg(tf1), tf3))
})
.py_try({
.py_rec("lti:test_TransferFunction_addition_and_subtraction:461", function() .py_add(.py_sub(tf1, tf2), tf3))
})
.py_try({
.py_rec("lti:test_TransferFunction_addition_and_subtraction:463", function() .py_sub(tf1, .py_Matrix(list(1L, 2L, 3L))))
})
.py_try({
.py_rec("lti:test_TransferFunction_addition_and_subtraction:464", function() .py_sub(tf3, tf4))
})
.py_try({
.py_rec("lti:test_TransferFunction_addition_and_subtraction:465", function() .py_sub(tf1, .py_sub(s, 1L)))
})
.py_try({
.py_rec("lti:test_TransferFunction_addition_and_subtraction:466", function() .py_sub(tf1, 8L))
})
.py_try({
.py_rec("lti:test_TransferFunction_addition_and_subtraction:467", function() .py_sub(.py_add(s, 5L), tf2))
})
.py_try({
.py_rec("lti:test_TransferFunction_addition_and_subtraction:468", function() .py_sub(.py_add(1L, .py_pow(p, 4L)), tf1))
})
.py_try({
dtf1 <- DiscreteTransferFunction(s, .py_add(s, 1L), s, .py_flt(0.01))
})
.py_try({
dtf2 <- DiscreteTransferFunction(.py_add(s, 1L), .py_add(s, 2L), s, .py_flt(0.01))
})
.py_try({
.py_rec("lti:test_TransferFunction_addition_and_subtraction:473", function() .py_add(dtf1, tf1))
})
.py_try({
.py_rec("lti:test_TransferFunction_addition_and_subtraction:474", function() .py_sub(dtf1, tf1))
})
.py_try({
.py_rec("lti:test_TransferFunction_addition_and_subtraction:475", function() .py_sub(.py_add(dtf2, dtf1), tf1))
})
.py_try({
.py_rec("lti:test_TransferFunction_addition_and_subtraction:476", function() .py_add(.py_sub(dtf1, dtf2), tf1))
}))---")
  expect_identical(got[["lti:test_TransferFunction_addition_and_subtraction:444"]], "Parallel(TransferFunction(s + 6, s - 5, s), TransferFunction(s + 3, s + 1, s))")
  expect_identical(got[["lti:test_TransferFunction_addition_and_subtraction:445"]], "Parallel(TransferFunction(s + 1, s^2 + s + 1, s), TransferFunction(s + 6, s - 5, s))")
  expect_identical(got[["lti:test_TransferFunction_addition_and_subtraction:446"]], "Parallel(TransferFunction(-s - 6, s - 5, s), TransferFunction(s + 3, s + 1, s), TransferFunction(s + 1, s^2 + s + 1, s))")
  expect_identical(got[["lti:test_TransferFunction_addition_and_subtraction:447"]], "Parallel(TransferFunction(s + 6, s - 5, s), TransferFunction(s + 3, s + 1, s), TransferFunction(s + 1, s^2 + s + 1, s))")
  expect_identical(got[["lti:test_TransferFunction_addition_and_subtraction:450"]], "ERR ValueError: <class 'sympy.physics.control.lti.TransferFunction'> cannot be added with <class 'sympy.matrices.dense.MutableDenseMatrix'>.")
  expect_identical(got[["lti:test_TransferFunction_addition_and_subtraction:451"]], "ERR ValueError: <class 'sympy.physics.control.lti.TransferFunction'> cannot be added with <class 'sympy.core.symbol.Symbol'>.")
  expect_identical(got[["lti:test_TransferFunction_addition_and_subtraction:452"]], "ERR ValueError: \\nAll the transfer functions should use the same complex variable of the\\nLaplace domain or z-domain.")
  expect_identical(got[["lti:test_TransferFunction_addition_and_subtraction:453"]], "ERR ValueError: <class 'sympy.physics.control.lti.TransferFunction'> cannot be added with <class 'sympy.core.add.Add'>.")
  expect_identical(got[["lti:test_TransferFunction_addition_and_subtraction:454"]], "ERR ValueError: <class 'sympy.physics.control.lti.TransferFunction'> cannot be added with <class 'int'>.")
  expect_identical(got[["lti:test_TransferFunction_addition_and_subtraction:455"]], "ERR ValueError: <class 'sympy.physics.control.lti.TransferFunction'> cannot be added with <class 'sympy.core.add.Add'>.")
  expect_identical(got[["lti:test_TransferFunction_addition_and_subtraction:458"]], "Parallel(TransferFunction(s + 6, s - 5, s), TransferFunction(-s - 3, s + 1, s))")
  expect_identical(got[["lti:test_TransferFunction_addition_and_subtraction:459"]], "Parallel(TransferFunction(s + 1, s^2 + s + 1, s), TransferFunction(-s - 3, s + 1, s))")
  expect_identical(got[["lti:test_TransferFunction_addition_and_subtraction:460"]], "Parallel(TransferFunction(-s - 6, s - 5, s), TransferFunction(-s - 1, s^2 + s + 1, s))")
  expect_identical(got[["lti:test_TransferFunction_addition_and_subtraction:461"]], "Parallel(TransferFunction(s + 6, s - 5, s), TransferFunction(-s - 3, s + 1, s), TransferFunction(s + 1, s^2 + s + 1, s))")
  expect_identical(got[["lti:test_TransferFunction_addition_and_subtraction:463"]], "ERR ValueError: <class 'sympy.physics.control.lti.TransferFunction'> cannot be subtracted from a <class 'sympy.matrices.dense.MutableDenseMatrix'>.")
  expect_identical(got[["lti:test_TransferFunction_addition_and_subtraction:464"]], "ERR ValueError: \\nAll the transfer functions should use the same complex variable of the\\nLaplace domain or z-domain.")
  expect_identical(got[["lti:test_TransferFunction_addition_and_subtraction:465"]], "ERR ValueError: <class 'sympy.physics.control.lti.TransferFunction'> cannot be subtracted from a <class 'sympy.core.add.Add'>.")
  expect_identical(got[["lti:test_TransferFunction_addition_and_subtraction:466"]], "ERR ValueError: <class 'sympy.physics.control.lti.TransferFunction'> cannot be subtracted from a <class 'int'>.")
  expect_identical(got[["lti:test_TransferFunction_addition_and_subtraction:467"]], "ERR ValueError: <class 'sympy.physics.control.lti.TransferFunction'> cannot be added with <class 'sympy.core.add.Add'>.")
  expect_identical(got[["lti:test_TransferFunction_addition_and_subtraction:468"]], "ERR ValueError: <class 'sympy.physics.control.lti.TransferFunction'> cannot be added with <class 'sympy.core.add.Add'>.")
  expect_identical(got[["lti:test_TransferFunction_addition_and_subtraction:473"]], "ERR TypeError: \\n                Incompatible systems detected.\\n                All systems should be either continuous-time or discrete-time.\\n                Found: DiscreteTransferFunction(s, s + 1, s, 0.01) and TransferFunction(s + 6, s - 5, s)")
  expect_identical(got[["lti:test_TransferFunction_addition_and_subtraction:474"]], "ERR TypeError: \\n                Incompatible systems detected.\\n                All systems should be either continuous-time or discrete-time.\\n                Found: DiscreteTransferFunction(s, s + 1, s, 0.01) and TransferFunction(s + 6, s - 5, s)")
  expect_identical(got[["lti:test_TransferFunction_addition_and_subtraction:475"]], "ERR TypeError: \\n                Incompatible systems detected.\\n                All systems should be either continuous-time or discrete-time.\\n                Found: Parallel(DiscreteTransferFunction(s + 1, s + 2, s, 0.01), DiscreteTransferFunction(s, s + 1, s, 0.01)) and TransferFunction(s + 6, s - 5, s)")
  expect_identical(got[["lti:test_TransferFunction_addition_and_subtraction:476"]], "ERR TypeError: \\n                Incompatible systems detected.\\n                All systems should be either continuous-time or discrete-time.\\n                Found: Parallel(DiscreteTransferFunction(s, s + 1, s, 0.01), DiscreteTransferFunction(-s - 1, s + 2, s, 0.01)) and TransferFunction(s + 6, s - 5, s)")
})

test_that("lti test_TransferFunction_multiplication_and_division matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("a, x, b, c, s, g, d,    p, k, tau, zeta, wn, T, z"); a <- .py_getitem(.tmpu_, 0L); x <- .py_getitem(.tmpu_, 1L); b <- .py_getitem(.tmpu_, 2L); c <- .py_getitem(.tmpu_, 3L); s <- .py_getitem(.tmpu_, 4L); g <- .py_getitem(.tmpu_, 5L); d <- .py_getitem(.tmpu_, 6L); p <- .py_getitem(.tmpu_, 7L); k <- .py_getitem(.tmpu_, 8L); tau <- .py_getitem(.tmpu_, 9L); zeta <- .py_getitem(.tmpu_, 10L); wn <- .py_getitem(.tmpu_, 11L); T <- .py_getitem(.tmpu_, 12L); z <- .py_getitem(.tmpu_, 13L)
})
.py_try({
.tmpu_ <- .py_symbols("a0:4, b0:5, c0:4, d0:4"); a0 <- .py_getitem(.tmpu_, 0L); a1 <- .py_getitem(.tmpu_, 1L); a2 <- .py_getitem(.tmpu_, 2L); a3 <- .py_getitem(.tmpu_, 3L); b0 <- .py_getitem(.tmpu_, 4L); b1 <- .py_getitem(.tmpu_, 5L); b2 <- .py_getitem(.tmpu_, 6L); b3 <- .py_getitem(.tmpu_, 7L); b4 <- .py_getitem(.tmpu_, 8L); c0 <- .py_getitem(.tmpu_, 9L); c1 <- .py_getitem(.tmpu_, 10L); c2 <- .py_getitem(.tmpu_, 11L); c3 <- .py_getitem(.tmpu_, 12L); d0 <- .py_getitem(.tmpu_, 13L); d1 <- .py_getitem(.tmpu_, 14L); d2 <- .py_getitem(.tmpu_, 15L); d3 <- .py_getitem(.tmpu_, 16L)
})
.py_try({
TF1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
TF2 <- TransferFunction(k, 1L, s)
})
.py_try({
TF3 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
TF4 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(tau, s), 1L), s)
})
.py_try({
TF5 <- TransferFunction(k, p, s)
})
.py_try({
TF6 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(k, s), p), s)
})
p_are_floats_equals <- function(n1, n2) {
  "\n    Helper function to compare two floats.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return(.py_method(Float(n1), "is_same", Float(n2), .py_isclose))
}
p_are_complex_equals <- function(n1, n2) {
  "\n    Helper function to compare two complex numbers.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return((.py_truthy(p_are_floats_equals(re(n1), re(n2))) && .py_truthy(p_are_floats_equals(im(n1), im(n2)))))
}
.py_try({
G1 <- TransferFunction(.py_add(s, 3L), .py_add(.py_neg(.py_pow(s, 3L)), 9L), s)
})
.py_try({
G2 <- TransferFunction(.py_add(s, 1L), .py_sub(s, 5L), s)
})
.py_try({
G3 <- TransferFunction(p, .py_sub(.py_pow(p, 4L), 6L), p)
})
.py_try({
G4 <- TransferFunction(.py_add(p, 4L), .py_sub(p, 5L), p)
})
.py_try({
G5 <- TransferFunction(.py_add(s, 6L), .py_sub(s, 5L), s)
})
.py_try({
G6 <- TransferFunction(.py_add(s, 3L), .py_add(s, 1L), s)
})
.py_try({
G7 <- TransferFunction(1L, 1L, s)
})
.py_try({
.py_rec("lti:test_TransferFunction_multiplication_and_division:488", function() .py_mul(G1, G2))
})
.py_try({
.py_rec("lti:test_TransferFunction_multiplication_and_division:489", function() .py_mul(.py_neg(G1), G5))
})
.py_try({
.py_rec("lti:test_TransferFunction_multiplication_and_division:490", function() .py_mul(.py_mul(.py_neg(G2), G5), .py_neg(G6)))
})
.py_try({
.py_rec("lti:test_TransferFunction_multiplication_and_division:491", function() .py_mul(.py_mul(.py_mul(.py_neg(G1), .py_neg(G2)), .py_neg(G5)), .py_neg(G6)))
})
.py_try({
.py_rec("lti:test_TransferFunction_multiplication_and_division:492", function() .py_mul(G3, G4))
})
.py_try({
.py_rec("lti:test_TransferFunction_multiplication_and_division:493", function() .py_mul(.py_mul(G1, G2), .py_neg(.py_mul(G5, G6))))
})
.py_try({
.py_rec("lti:test_TransferFunction_multiplication_and_division:495", function() .py_mul(.py_mul(G1, G2), .py_add(G5, G6)))
})
.py_try({
.py_rec("lti:test_TransferFunction_multiplication_and_division:498", function() .py_div(G5, G6))
})
.py_try({
.py_rec("lti:test_TransferFunction_multiplication_and_division:499", function() .py_div(.py_neg(G3), G4))
})
.py_try({
.py_rec("lti:test_TransferFunction_multiplication_and_division:500", function() .py_div(.py_mul(G5, G6), G7))
})
.py_try({
c <- .py_symbols("c", commutative = FALSE)
})
.py_try({
.py_rec("lti:test_TransferFunction_multiplication_and_division:503", function() .py_mul(G3, .py_Matrix(list(1L, 2L, 3L))))
})
.py_try({
.py_rec("lti:test_TransferFunction_multiplication_and_division:504", function() .py_mul(G1, c))
})
.py_try({
.py_rec("lti:test_TransferFunction_multiplication_and_division:505", function() .py_mul(G3, G5))
})
.py_try({
.py_rec("lti:test_TransferFunction_multiplication_and_division:506", function() .py_mul(G5, .py_sub(s, 1L)))
})
.py_try({
.py_rec("lti:test_TransferFunction_multiplication_and_division:507", function() .py_mul(9L, G5))
})
.py_try({
.py_rec("lti:test_TransferFunction_multiplication_and_division:509", function() .py_div(G3, .py_Matrix(list(1L, 2L, 3L))))
})
.py_try({
.py_rec("lti:test_TransferFunction_multiplication_and_division:510", function() .py_div(G6, 0L))
})
.py_try({
.py_rec("lti:test_TransferFunction_multiplication_and_division:511", function() .py_div(G3, G5))
})
.py_try({
.py_rec("lti:test_TransferFunction_multiplication_and_division:512", function() .py_div(G5, 2L))
})
.py_try({
.py_rec("lti:test_TransferFunction_multiplication_and_division:513", function() .py_div(G5, .py_pow(s, 2L)))
})
.py_try({
.py_rec("lti:test_TransferFunction_multiplication_and_division:514", function() .py_div(.py_sub(s, .py_mul(4L, .py_pow(s, 2L))), G2))
})
.py_try({
.py_rec("lti:test_TransferFunction_multiplication_and_division:515", function() .py_div(0L, G4))
})
.py_try({
.py_rec("lti:test_TransferFunction_multiplication_and_division:516", function() .py_div(G7, .py_add(1L, G6)))
})
.py_try({
.py_rec("lti:test_TransferFunction_multiplication_and_division:517", function() .py_div(G7, .py_mul(G5, G6)))
})
.py_try({
.py_rec("lti:test_TransferFunction_multiplication_and_division:518", function() .py_div(G7, .py_add(G7, .py_add(G5, G6))))
})
.py_try({
dtf1 <- DiscreteTransferFunction(s, .py_add(s, 1L), s, .py_flt(0.01))
})
.py_try({
dtf2 <- DiscreteTransferFunction(.py_add(s, 1L), .py_add(s, 2L), s, .py_flt(0.01))
})
.py_try({
.py_rec("lti:test_TransferFunction_multiplication_and_division:523", function() .py_mul(dtf1, G1))
})
.py_try({
.py_rec("lti:test_TransferFunction_multiplication_and_division:524", function() .py_div(dtf1, G1))
})
.py_try({
.py_rec("lti:test_TransferFunction_multiplication_and_division:525", function() .py_div(.py_mul(dtf1, dtf2), G2))
})
.py_try({
.py_rec("lti:test_TransferFunction_multiplication_and_division:526", function() .py_mul(.py_div(dtf1, dtf2), G2))
}))---")
  expect_identical(got[["lti:test_TransferFunction_multiplication_and_division:488"]], "Series(TransferFunction(s + 3, 9 - s^3, s), TransferFunction(s + 1, s - 5, s))")
  expect_identical(got[["lti:test_TransferFunction_multiplication_and_division:489"]], "Series(TransferFunction(-s - 3, 9 - s^3, s), TransferFunction(s + 6, s - 5, s))")
  expect_identical(got[["lti:test_TransferFunction_multiplication_and_division:490"]], "Series(TransferFunction(-s - 1, s - 5, s), TransferFunction(s + 6, s - 5, s), TransferFunction(-s - 3, s + 1, s))")
  expect_identical(got[["lti:test_TransferFunction_multiplication_and_division:491"]], "Series(TransferFunction(-s - 3, 9 - s^3, s), TransferFunction(-s - 1, s - 5, s), TransferFunction(-s - 6, s - 5, s), TransferFunction(-s - 3, s + 1, s))")
  expect_identical(got[["lti:test_TransferFunction_multiplication_and_division:492"]], "Series(TransferFunction(p, p^4 - 6, p), TransferFunction(p + 4, p - 5, p))")
  expect_identical(got[["lti:test_TransferFunction_multiplication_and_division:493"]], "Series(TransferFunction(s + 3, 9 - s^3, s), TransferFunction(s + 1, s - 5, s), TransferFunction(-1, 1, s), TransferFunction(s + 6, s - 5, s), TransferFunction(s + 3, s + 1, s))")
  expect_identical(got[["lti:test_TransferFunction_multiplication_and_division:495"]], "Series(TransferFunction(s + 3, 9 - s^3, s), TransferFunction(s + 1, s - 5, s), Parallel(TransferFunction(s + 6, s - 5, s), TransferFunction(s + 3, s + 1, s)))")
  expect_identical(got[["lti:test_TransferFunction_multiplication_and_division:498"]], "Series(TransferFunction(s + 6, s - 5, s), TransferFunction(s + 1, s + 3, s))")
  expect_identical(got[["lti:test_TransferFunction_multiplication_and_division:499"]], "Series(TransferFunction(-p, p^4 - 6, p), TransferFunction(p - 5, p + 4, p))")
  expect_identical(got[["lti:test_TransferFunction_multiplication_and_division:500"]], "Series(TransferFunction(s + 6, s - 5, s), TransferFunction(s + 3, s + 1, s), TransferFunction(1, 1, s))")
  expect_identical(got[["lti:test_TransferFunction_multiplication_and_division:503"]], "ERR ValueError: <class 'sympy.physics.control.lti.TransferFunction'> cannot be multiplied with <class 'sympy.matrices.dense.MutableDenseMatrix'>.")
  expect_identical(got[["lti:test_TransferFunction_multiplication_and_division:504"]], "ERR ValueError: <class 'sympy.physics.control.lti.TransferFunction'> cannot be multiplied with <class 'sympy.core.symbol.Symbol'>.")
  expect_identical(got[["lti:test_TransferFunction_multiplication_and_division:505"]], "ERR ValueError: \\nAll the transfer functions should use the same complex variable of the\\nLaplace domain or z-domain.")
  expect_identical(got[["lti:test_TransferFunction_multiplication_and_division:506"]], "ERR ValueError: <class 'sympy.physics.control.lti.TransferFunction'> cannot be multiplied with <class 'sympy.core.add.Add'>.")
  expect_identical(got[["lti:test_TransferFunction_multiplication_and_division:507"]], "ERR ValueError: <class 'sympy.physics.control.lti.TransferFunction'> cannot be multiplied with <class 'int'>.")
  expect_identical(got[["lti:test_TransferFunction_multiplication_and_division:509"]], "ERR ValueError: <class 'sympy.physics.control.lti.TransferFunction'> cannot be divided by <class 'sympy.matrices.dense.MutableDenseMatrix'>.")
  expect_identical(got[["lti:test_TransferFunction_multiplication_and_division:510"]], "ERR ValueError: <class 'sympy.physics.control.lti.TransferFunction'> cannot be divided by <class 'int'>.")
  expect_identical(got[["lti:test_TransferFunction_multiplication_and_division:511"]], "ERR ValueError: \\nAll the transfer functions should use the same complex variable of the\\nLaplace domain or z-domain.")
  expect_identical(got[["lti:test_TransferFunction_multiplication_and_division:512"]], "ERR ValueError: <class 'sympy.physics.control.lti.TransferFunction'> cannot be divided by <class 'int'>.")
  expect_identical(got[["lti:test_TransferFunction_multiplication_and_division:513"]], "ERR ValueError: <class 'sympy.physics.control.lti.TransferFunction'> cannot be divided by <class 'sympy.core.power.Pow'>.")
  expect_identical(got[["lti:test_TransferFunction_multiplication_and_division:514"]], "ERR ValueError: <class 'sympy.physics.control.lti.TransferFunction'> cannot be divided by <class 'sympy.core.add.Add'>.")
  expect_identical(got[["lti:test_TransferFunction_multiplication_and_division:515"]], "ERR ValueError: <class 'sympy.physics.control.lti.TransferFunction'> cannot be divided by <class 'int'>.")
  expect_identical(got[["lti:test_TransferFunction_multiplication_and_division:516"]], "ERR ValueError: <class 'sympy.physics.control.lti.TransferFunction'> cannot be added with <class 'int'>.")
  expect_identical(got[["lti:test_TransferFunction_multiplication_and_division:517"]], "ERR ValueError: <class 'sympy.physics.control.lti.TransferFunction'> cannot be divided by <class 'sympy.physics.control.lti.Series'>.")
  expect_identical(got[["lti:test_TransferFunction_multiplication_and_division:518"]], "ERR ValueError: <class 'sympy.physics.control.lti.TransferFunction'> cannot be divided by <class 'sympy.physics.control.lti.Parallel'>.")
  expect_identical(got[["lti:test_TransferFunction_multiplication_and_division:523"]], "ERR TypeError: \\n                Incompatible systems detected.\\n                All systems should be either continuous-time or discrete-time.\\n                Found: DiscreteTransferFunction(s, s + 1, s, 0.01) and TransferFunction(s + 3, 9 - s^3, s)")
  expect_identical(got[["lti:test_TransferFunction_multiplication_and_division:524"]], "ERR TypeError: \\n                Incompatible systems detected.\\n                All systems should be either continuous-time or discrete-time.\\n                Found: DiscreteTransferFunction(s, s + 1, s, 0.01) and TransferFunction(s + 3, 9 - s^3, s)")
  expect_identical(got[["lti:test_TransferFunction_multiplication_and_division:525"]], "ERR TypeError: \\n                Incompatible systems detected.\\n                All systems should be either continuous-time or discrete-time.\\n                Found: Series(DiscreteTransferFunction(s, s + 1, s, 0.01), DiscreteTransferFunction(s + 1, s + 2, s, 0.01)) and TransferFunction(s + 1, s - 5, s)")
  expect_identical(got[["lti:test_TransferFunction_multiplication_and_division:526"]], "ERR TypeError: \\n                Incompatible systems detected.\\n                All systems should be either continuous-time or discrete-time.\\n                Found: Series(DiscreteTransferFunction(s, s + 1, s, 0.01), DiscreteTransferFunction(s + 2, s + 1, s, 0.01)) and TransferFunction(s + 1, s - 5, s)")
})

test_that("lti test_TransferFunction_is_proper matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("a, x, b, c, s, g, d,    p, k, tau, zeta, wn, T, z"); a <- .py_getitem(.tmpu_, 0L); x <- .py_getitem(.tmpu_, 1L); b <- .py_getitem(.tmpu_, 2L); c <- .py_getitem(.tmpu_, 3L); s <- .py_getitem(.tmpu_, 4L); g <- .py_getitem(.tmpu_, 5L); d <- .py_getitem(.tmpu_, 6L); p <- .py_getitem(.tmpu_, 7L); k <- .py_getitem(.tmpu_, 8L); tau <- .py_getitem(.tmpu_, 9L); zeta <- .py_getitem(.tmpu_, 10L); wn <- .py_getitem(.tmpu_, 11L); T <- .py_getitem(.tmpu_, 12L); z <- .py_getitem(.tmpu_, 13L)
})
.py_try({
.tmpu_ <- .py_symbols("a0:4, b0:5, c0:4, d0:4"); a0 <- .py_getitem(.tmpu_, 0L); a1 <- .py_getitem(.tmpu_, 1L); a2 <- .py_getitem(.tmpu_, 2L); a3 <- .py_getitem(.tmpu_, 3L); b0 <- .py_getitem(.tmpu_, 4L); b1 <- .py_getitem(.tmpu_, 5L); b2 <- .py_getitem(.tmpu_, 6L); b3 <- .py_getitem(.tmpu_, 7L); b4 <- .py_getitem(.tmpu_, 8L); c0 <- .py_getitem(.tmpu_, 9L); c1 <- .py_getitem(.tmpu_, 10L); c2 <- .py_getitem(.tmpu_, 11L); c3 <- .py_getitem(.tmpu_, 12L); d0 <- .py_getitem(.tmpu_, 13L); d1 <- .py_getitem(.tmpu_, 14L); d2 <- .py_getitem(.tmpu_, 15L); d3 <- .py_getitem(.tmpu_, 16L)
})
.py_try({
TF1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
TF2 <- TransferFunction(k, 1L, s)
})
.py_try({
TF3 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
TF4 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(tau, s), 1L), s)
})
.py_try({
TF5 <- TransferFunction(k, p, s)
})
.py_try({
TF6 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(k, s), p), s)
})
p_are_floats_equals <- function(n1, n2) {
  "\n    Helper function to compare two floats.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return(.py_method(Float(n1), "is_same", Float(n2), .py_isclose))
}
p_are_complex_equals <- function(n1, n2) {
  "\n    Helper function to compare two complex numbers.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return((.py_truthy(p_are_floats_equals(re(n1), re(n2))) && .py_truthy(p_are_floats_equals(im(n1), im(n2)))))
}
.py_try({
.tmpu_ <- .py_symbols("omega_o, zeta, tau"); omega_o <- .py_getitem(.tmpu_, 0L); zeta <- .py_getitem(.tmpu_, 1L); tau <- .py_getitem(.tmpu_, 2L)
})
.py_try({
G1 <- TransferFunction(.py_pow(omega_o, 2L), .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(p, omega_o), zeta), s)), .py_pow(omega_o, 2L)), omega_o)
})
.py_try({
G2 <- TransferFunction(.py_sub(tau, .py_pow(s, 3L)), .py_add(tau, .py_pow(p, 4L)), tau)
})
.py_try({
G3 <- TransferFunction(.py_add(.py_sub(.py_add(.py_mul(.py_mul(a, b), .py_pow(s, 3L)), .py_pow(s, 2L)), .py_mul(a, p)), s), .py_sub(b, .py_mul(s, .py_pow(p, 2L))), p)
})
.py_try({
G4 <- TransferFunction(.py_add(.py_sub(.py_add(.py_mul(b, .py_pow(s, 2L)), .py_pow(p, 2L)), .py_mul(a, p)), s), .py_sub(b, .py_pow(p, 2L)), s)
})
.py_try({
.py_rec("lti:test_TransferFunction_is_proper:535", function() .py_attr(G1, "is_proper"))
})
.py_try({
.py_rec("lti:test_TransferFunction_is_proper:536", function() .py_attr(G2, "is_proper"))
})
.py_try({
.py_rec("lti:test_TransferFunction_is_proper:537", function() .py_attr(G3, "is_proper"))
})
.py_try({
.py_rec("lti:test_TransferFunction_is_proper:538", function() .py_attr(G4, "is_proper"))
}))---")
  expect_identical(got[["lti:test_TransferFunction_is_proper:535"]], "True")
  expect_identical(got[["lti:test_TransferFunction_is_proper:536"]], "True")
  expect_identical(got[["lti:test_TransferFunction_is_proper:537"]], "True")
  expect_identical(got[["lti:test_TransferFunction_is_proper:538"]], "False")
})

test_that("lti test_TransferFunction_is_strictly_proper matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("a, x, b, c, s, g, d,    p, k, tau, zeta, wn, T, z"); a <- .py_getitem(.tmpu_, 0L); x <- .py_getitem(.tmpu_, 1L); b <- .py_getitem(.tmpu_, 2L); c <- .py_getitem(.tmpu_, 3L); s <- .py_getitem(.tmpu_, 4L); g <- .py_getitem(.tmpu_, 5L); d <- .py_getitem(.tmpu_, 6L); p <- .py_getitem(.tmpu_, 7L); k <- .py_getitem(.tmpu_, 8L); tau <- .py_getitem(.tmpu_, 9L); zeta <- .py_getitem(.tmpu_, 10L); wn <- .py_getitem(.tmpu_, 11L); T <- .py_getitem(.tmpu_, 12L); z <- .py_getitem(.tmpu_, 13L)
})
.py_try({
.tmpu_ <- .py_symbols("a0:4, b0:5, c0:4, d0:4"); a0 <- .py_getitem(.tmpu_, 0L); a1 <- .py_getitem(.tmpu_, 1L); a2 <- .py_getitem(.tmpu_, 2L); a3 <- .py_getitem(.tmpu_, 3L); b0 <- .py_getitem(.tmpu_, 4L); b1 <- .py_getitem(.tmpu_, 5L); b2 <- .py_getitem(.tmpu_, 6L); b3 <- .py_getitem(.tmpu_, 7L); b4 <- .py_getitem(.tmpu_, 8L); c0 <- .py_getitem(.tmpu_, 9L); c1 <- .py_getitem(.tmpu_, 10L); c2 <- .py_getitem(.tmpu_, 11L); c3 <- .py_getitem(.tmpu_, 12L); d0 <- .py_getitem(.tmpu_, 13L); d1 <- .py_getitem(.tmpu_, 14L); d2 <- .py_getitem(.tmpu_, 15L); d3 <- .py_getitem(.tmpu_, 16L)
})
.py_try({
TF1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
TF2 <- TransferFunction(k, 1L, s)
})
.py_try({
TF3 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
TF4 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(tau, s), 1L), s)
})
.py_try({
TF5 <- TransferFunction(k, p, s)
})
.py_try({
TF6 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(k, s), p), s)
})
p_are_floats_equals <- function(n1, n2) {
  "\n    Helper function to compare two floats.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return(.py_method(Float(n1), "is_same", Float(n2), .py_isclose))
}
p_are_complex_equals <- function(n1, n2) {
  "\n    Helper function to compare two complex numbers.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return((.py_truthy(p_are_floats_equals(re(n1), re(n2))) && .py_truthy(p_are_floats_equals(im(n1), im(n2)))))
}
.py_try({
.tmpu_ <- .py_symbols("omega_o, zeta, tau"); omega_o <- .py_getitem(.tmpu_, 0L); zeta <- .py_getitem(.tmpu_, 1L); tau <- .py_getitem(.tmpu_, 2L)
})
.py_try({
tf1 <- TransferFunction(.py_pow(omega_o, 2L), .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(p, omega_o), zeta), s)), .py_pow(omega_o, 2L)), omega_o)
})
.py_try({
tf2 <- TransferFunction(.py_sub(tau, .py_pow(s, 3L)), .py_add(tau, .py_pow(p, 4L)), tau)
})
.py_try({
tf3 <- TransferFunction(.py_add(.py_sub(.py_add(.py_mul(.py_mul(a, b), .py_pow(s, 3L)), .py_pow(s, 2L)), .py_mul(a, p)), s), .py_sub(b, .py_mul(s, .py_pow(p, 2L))), p)
})
.py_try({
tf4 <- TransferFunction(.py_add(.py_sub(.py_add(.py_mul(b, .py_pow(s, 2L)), .py_pow(p, 2L)), .py_mul(a, p)), s), .py_sub(b, .py_pow(p, 2L)), s)
})
.py_try({
.py_rec("lti:test_TransferFunction_is_strictly_proper:547", function() .py_attr(tf1, "is_strictly_proper"))
})
.py_try({
.py_rec("lti:test_TransferFunction_is_strictly_proper:548", function() .py_attr(tf2, "is_strictly_proper"))
})
.py_try({
.py_rec("lti:test_TransferFunction_is_strictly_proper:549", function() .py_attr(tf3, "is_strictly_proper"))
})
.py_try({
.py_rec("lti:test_TransferFunction_is_strictly_proper:550", function() .py_attr(tf4, "is_strictly_proper"))
}))---")
  expect_identical(got[["lti:test_TransferFunction_is_strictly_proper:547"]], "False")
  expect_identical(got[["lti:test_TransferFunction_is_strictly_proper:548"]], "False")
  expect_identical(got[["lti:test_TransferFunction_is_strictly_proper:549"]], "True")
  expect_identical(got[["lti:test_TransferFunction_is_strictly_proper:550"]], "False")
})

test_that("lti test_TransferFunction_is_biproper matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("a, x, b, c, s, g, d,    p, k, tau, zeta, wn, T, z"); a <- .py_getitem(.tmpu_, 0L); x <- .py_getitem(.tmpu_, 1L); b <- .py_getitem(.tmpu_, 2L); c <- .py_getitem(.tmpu_, 3L); s <- .py_getitem(.tmpu_, 4L); g <- .py_getitem(.tmpu_, 5L); d <- .py_getitem(.tmpu_, 6L); p <- .py_getitem(.tmpu_, 7L); k <- .py_getitem(.tmpu_, 8L); tau <- .py_getitem(.tmpu_, 9L); zeta <- .py_getitem(.tmpu_, 10L); wn <- .py_getitem(.tmpu_, 11L); T <- .py_getitem(.tmpu_, 12L); z <- .py_getitem(.tmpu_, 13L)
})
.py_try({
.tmpu_ <- .py_symbols("a0:4, b0:5, c0:4, d0:4"); a0 <- .py_getitem(.tmpu_, 0L); a1 <- .py_getitem(.tmpu_, 1L); a2 <- .py_getitem(.tmpu_, 2L); a3 <- .py_getitem(.tmpu_, 3L); b0 <- .py_getitem(.tmpu_, 4L); b1 <- .py_getitem(.tmpu_, 5L); b2 <- .py_getitem(.tmpu_, 6L); b3 <- .py_getitem(.tmpu_, 7L); b4 <- .py_getitem(.tmpu_, 8L); c0 <- .py_getitem(.tmpu_, 9L); c1 <- .py_getitem(.tmpu_, 10L); c2 <- .py_getitem(.tmpu_, 11L); c3 <- .py_getitem(.tmpu_, 12L); d0 <- .py_getitem(.tmpu_, 13L); d1 <- .py_getitem(.tmpu_, 14L); d2 <- .py_getitem(.tmpu_, 15L); d3 <- .py_getitem(.tmpu_, 16L)
})
.py_try({
TF1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
TF2 <- TransferFunction(k, 1L, s)
})
.py_try({
TF3 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
TF4 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(tau, s), 1L), s)
})
.py_try({
TF5 <- TransferFunction(k, p, s)
})
.py_try({
TF6 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(k, s), p), s)
})
p_are_floats_equals <- function(n1, n2) {
  "\n    Helper function to compare two floats.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return(.py_method(Float(n1), "is_same", Float(n2), .py_isclose))
}
p_are_complex_equals <- function(n1, n2) {
  "\n    Helper function to compare two complex numbers.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return((.py_truthy(p_are_floats_equals(re(n1), re(n2))) && .py_truthy(p_are_floats_equals(im(n1), im(n2)))))
}
.py_try({
.tmpu_ <- .py_symbols("tau, omega_o, zeta"); tau <- .py_getitem(.tmpu_, 0L); omega_o <- .py_getitem(.tmpu_, 1L); zeta <- .py_getitem(.tmpu_, 2L)
})
.py_try({
tf1 <- TransferFunction(.py_pow(omega_o, 2L), .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(p, omega_o), zeta), s)), .py_pow(omega_o, 2L)), omega_o)
})
.py_try({
tf2 <- TransferFunction(.py_sub(tau, .py_pow(s, 3L)), .py_add(tau, .py_pow(p, 4L)), tau)
})
.py_try({
tf3 <- TransferFunction(.py_add(.py_sub(.py_add(.py_mul(.py_mul(a, b), .py_pow(s, 3L)), .py_pow(s, 2L)), .py_mul(a, p)), s), .py_sub(b, .py_mul(s, .py_pow(p, 2L))), p)
})
.py_try({
tf4 <- TransferFunction(.py_add(.py_sub(.py_add(.py_mul(b, .py_pow(s, 2L)), .py_pow(p, 2L)), .py_mul(a, p)), s), .py_sub(b, .py_pow(p, 2L)), s)
})
.py_try({
.py_rec("lti:test_TransferFunction_is_biproper:559", function() .py_attr(tf1, "is_biproper"))
})
.py_try({
.py_rec("lti:test_TransferFunction_is_biproper:560", function() .py_attr(tf2, "is_biproper"))
})
.py_try({
.py_rec("lti:test_TransferFunction_is_biproper:561", function() .py_attr(tf3, "is_biproper"))
})
.py_try({
.py_rec("lti:test_TransferFunction_is_biproper:562", function() .py_attr(tf4, "is_biproper"))
}))---")
  expect_identical(got[["lti:test_TransferFunction_is_biproper:559"]], "True")
  expect_identical(got[["lti:test_TransferFunction_is_biproper:560"]], "True")
  expect_identical(got[["lti:test_TransferFunction_is_biproper:561"]], "False")
  expect_identical(got[["lti:test_TransferFunction_is_biproper:562"]], "False")
})

test_that("lti test_DiscreteTransferFunction_construction matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("a, x, b, c, s, g, d,    p, k, tau, zeta, wn, T, z"); a <- .py_getitem(.tmpu_, 0L); x <- .py_getitem(.tmpu_, 1L); b <- .py_getitem(.tmpu_, 2L); c <- .py_getitem(.tmpu_, 3L); s <- .py_getitem(.tmpu_, 4L); g <- .py_getitem(.tmpu_, 5L); d <- .py_getitem(.tmpu_, 6L); p <- .py_getitem(.tmpu_, 7L); k <- .py_getitem(.tmpu_, 8L); tau <- .py_getitem(.tmpu_, 9L); zeta <- .py_getitem(.tmpu_, 10L); wn <- .py_getitem(.tmpu_, 11L); T <- .py_getitem(.tmpu_, 12L); z <- .py_getitem(.tmpu_, 13L)
})
.py_try({
.tmpu_ <- .py_symbols("a0:4, b0:5, c0:4, d0:4"); a0 <- .py_getitem(.tmpu_, 0L); a1 <- .py_getitem(.tmpu_, 1L); a2 <- .py_getitem(.tmpu_, 2L); a3 <- .py_getitem(.tmpu_, 3L); b0 <- .py_getitem(.tmpu_, 4L); b1 <- .py_getitem(.tmpu_, 5L); b2 <- .py_getitem(.tmpu_, 6L); b3 <- .py_getitem(.tmpu_, 7L); b4 <- .py_getitem(.tmpu_, 8L); c0 <- .py_getitem(.tmpu_, 9L); c1 <- .py_getitem(.tmpu_, 10L); c2 <- .py_getitem(.tmpu_, 11L); c3 <- .py_getitem(.tmpu_, 12L); d0 <- .py_getitem(.tmpu_, 13L); d1 <- .py_getitem(.tmpu_, 14L); d2 <- .py_getitem(.tmpu_, 15L); d3 <- .py_getitem(.tmpu_, 16L)
})
.py_try({
TF1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
TF2 <- TransferFunction(k, 1L, s)
})
.py_try({
TF3 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
TF4 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(tau, s), 1L), s)
})
.py_try({
TF5 <- TransferFunction(k, p, s)
})
.py_try({
TF6 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(k, s), p), s)
})
p_are_floats_equals <- function(n1, n2) {
  "\n    Helper function to compare two floats.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return(.py_method(Float(n1), "is_same", Float(n2), .py_isclose))
}
p_are_complex_equals <- function(n1, n2) {
  "\n    Helper function to compare two complex numbers.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return((.py_truthy(p_are_floats_equals(re(n1), re(n2))) && .py_truthy(p_are_floats_equals(im(n1), im(n2)))))
}
.py_try({
tf <- DiscreteTransferFunction(.py_add(z, 1L), .py_add(.py_add(.py_pow(z, 2L), z), 1L), z, 10L)
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_construction:567", function() .py_attr(tf, "num"))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_construction:568", function() .py_attr(tf, "den"))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_construction:569", function() .py_attr(tf, "args"))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_construction:570", function() .py_attr(tf, "sampling_time"))
})
.py_try({
tf1 <- DiscreteTransferFunction(.py_add(z, 4L), .py_sub(z, 5L), z, .py_flt(0.1))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_construction:573", function() .py_attr(tf1, "num"))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_construction:574", function() .py_attr(tf1, "den"))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_construction:575", function() .py_attr(tf1, "sampling_time"))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_construction:576", function() .py_isinstance(.py_attr(tf1, "sampling_time"), "Number"))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_construction:577", function() .py_attr(tf1, "args"))
})
.py_try({
tf2 <- DiscreteTransferFunction(.py_add(p, 3L), .py_sub(.py_pow(p, 2L), 9L), p, 12L)
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_construction:581", function() .py_attr(tf2, "num"))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_construction:582", function() .py_attr(tf2, "den"))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_construction:583", function() .py_attr(tf2, "sampling_time"))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_construction:584", function() .py_attr(tf2, "args"))
})
.py_try({
tf3 <- DiscreteTransferFunction(.py_add(.py_add(.py_pow(p, 3L), .py_mul(5L, .py_pow(p, 2L))), 4L), .py_add(.py_add(.py_pow(p, 4L), .py_mul(3L, p)), 1L), p)
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_construction:587", function() .py_attr(tf3, "args"))
})
.py_try({
tf4 <- DiscreteTransferFunction(.py_mul(.py_add(s, 3L), .py_sub(s, 1L)), .py_mul(.py_sub(s, 1L), .py_add(s, 5L)), s, .py_flt(6.4))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_construction:591", function() .py_attr(tf4, "den"))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_construction:592", function() .py_attr(tf4, "args"))
})
.py_try({
tf4_ <- DiscreteTransferFunction(.py_add(p, 2L), .py_add(p, 2L), p, .py_flt(1.2))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_construction:595", function() .py_attr(tf4_, "args"))
})
.py_try({
tf5 <- DiscreteTransferFunction(.py_sub(s, 1L), .py_sub(4L, p), s)
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_construction:598", function() .py_attr(tf5, "args"))
})
.py_try({
tf5_ <- DiscreteTransferFunction(.py_sub(s, 1L), .py_sub(s, 1L), s)
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_construction:601", function() .py_attr(tf5_, "args"))
})
.py_try({
tf6 <- DiscreteTransferFunction(5L, 6L, s)
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_construction:604", function() .py_attr(tf6, "num"))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_construction:605", function() .py_attr(tf6, "den"))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_construction:606", function() .py_attr(tf6, "args"))
})
.py_try({
tf6_ <- DiscreteTransferFunction(.py_flt(0.5), 4L, s)
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_construction:609", function() .py_attr(tf6_, "num"))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_construction:610", function() .py_attr(tf6_, "den"))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_construction:611", function() .py_attr(tf6_, "args"))
})
.py_try({
tf7 <- DiscreteTransferFunction(.py_add(.py_add(.py_mul(3L, .py_pow(s, 2L)), .py_mul(2L, p)), .py_mul(4L, s)), .py_add(.py_mul(8L, .py_pow(p, 2L)), .py_mul(7L, s)), s)
})
.py_try({
tf8 <- DiscreteTransferFunction(.py_add(.py_add(.py_mul(3L, .py_pow(s, 2L)), .py_mul(2L, p)), .py_mul(4L, s)), .py_add(.py_mul(8L, .py_pow(p, 2L)), .py_mul(7L, s)), p)
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_construction:615", function() tf7)
})
.py_try({
tf7_ <- DiscreteTransferFunction(.py_add(.py_add(.py_mul(a0, s), .py_mul(a1, .py_pow(s, 2L))), .py_mul(a2, .py_pow(s, 3L))), .py_sub(.py_mul(b0, p), .py_mul(b1, s)), s)
})
.py_try({
tf8_ <- DiscreteTransferFunction(.py_add(.py_add(.py_mul(a0, s), .py_mul(a1, .py_pow(s, 2L))), .py_mul(a2, .py_pow(s, 3L))), .py_sub(.py_mul(b0, p), .py_mul(b1, s)), s)
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_construction:619", function() tf7_)
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_construction:620", function() (.py_truthy(.py_eq(.py_neg(.py_neg(tf7_)), tf7_)) && .py_truthy(.py_eq(tf7_, .py_neg(.py_neg(.py_neg(.py_neg(tf7_))))))))
})
.py_try({
tf7_2 <- DiscreteTransferFunction(.py_add(.py_add(.py_mul(a0, s), .py_mul(a1, .py_pow(s, 2L))), .py_mul(a2, .py_pow(s, 3L))), .py_sub(.py_mul(b0, p), .py_mul(b1, s)), s, 2L)
})
.py_try({
tf8_2 <- DiscreteTransferFunction(.py_add(.py_add(.py_mul(a0, s), .py_mul(a1, .py_pow(s, 2L))), .py_mul(a2, .py_pow(s, 3L))), .py_sub(.py_mul(b0, p), .py_mul(b1, s)), s)
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_construction:624", function() tf7_2)
})
.py_try({
tf9 <- DiscreteTransferFunction(.py_add(.py_add(.py_add(.py_mul(a, .py_pow(s, 3L)), .py_mul(b, .py_pow(s, 2L))), .py_mul(g, s)), d), .py_add(.py_add(.py_mul(d, p), .py_mul(g, .py_pow(p, 2L))), .py_mul(g, s)), s)
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_construction:627", function() .py_attr(tf9, "args"))
})
.py_try({
tf10 <- DiscreteTransferFunction(.py_add(.py_pow(p, 3L), d), .py_add(.py_add(.py_mul(g, .py_pow(s, 2L)), .py_mul(d, s)), a), p, 2L)
})
.py_try({
tf10_ <- DiscreteTransferFunction(.py_add(.py_pow(p, 3L), d), .py_add(.py_add(.py_mul(g, .py_pow(s, 2L)), .py_mul(d, s)), a), p, 2L)
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_construction:631", function() .py_attr(tf10, "args"))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_construction:632", function() tf10_)
})
.py_try({
tf11 <- DiscreteTransferFunction(.py_add(.py_mul(a1, s), a0), .py_add(.py_add(.py_mul(b2, .py_pow(s, 2L)), .py_mul(b1, s)), b0), s)
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_construction:635", function() .py_attr(tf11, "num"))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_construction:636", function() .py_attr(tf11, "den"))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_construction:637", function() .py_attr(tf11, "args"))
})
.py_try({
tf12 <- DiscreteTransferFunction(0L, .py_add(.py_sub(.py_pow(p, 2L), p), 1L), p)
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_construction:641", function() .py_attr(tf12, "args"))
})
.py_try({
tf13 <- DiscreteTransferFunction(0L, 1L, s)
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_construction:644", function() .py_attr(tf13, "args"))
})
.py_try({
tf14 <- DiscreteTransferFunction(.py_sub(.py_add(.py_mul(a0, .py_pow(s, .py_flt(0.5))), .py_mul(a2, .py_pow(s, .py_flt(0.6)))), a1), .py_mul(a1, .py_pow(p, .py_flt(-8.7))), s, 6L)
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_construction:648", function() .py_attr(tf14, "args"))
})
.py_try({
tf15 <- DiscreteTransferFunction(.py_add(.py_mul(.py_pow(a2, 2L), .py_pow(p, .py_flt(0.25))), .py_mul(a1, .py_pow(s, .py_flt(-0.8)))), .py_sub(.py_mul(a0, s), p), p)
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_construction:651", function() .py_attr(tf15, "args"))
})
.py_try({
.tmpu_ <- .py_symbols("omega_o, k_p, k_o, k_i"); omega_o <- .py_getitem(.tmpu_, 0L); k_p <- .py_getitem(.tmpu_, 1L); k_o <- .py_getitem(.tmpu_, 2L); k_i <- .py_getitem(.tmpu_, 3L)
})
.py_try({
tf18 <- DiscreteTransferFunction(.py_add(.py_add(k_p, .py_mul(k_o, s)), .py_div(k_i, s)), .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(2L, omega_o), s)), .py_pow(omega_o, 2L)), s)
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_construction:656", function() .py_attr(tf18, "num"))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_construction:657", function() .py_attr(tf18, "args"))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_construction:661", function() DiscreteTransferFunction(4L, 0L, s))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_construction:662", function() DiscreteTransferFunction(s, 0L, s))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_construction:663", function() DiscreteTransferFunction(0L, 0L, s))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_construction:665", function() DiscreteTransferFunction(.py_Matrix(list(1L, 2L, 3L)), s, s))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_construction:667", function() DiscreteTransferFunction(.py_sub(.py_add(.py_pow(s, 2L), .py_mul(2L, s)), 1L), .py_add(s, 3L), 3L))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_construction:668", function() DiscreteTransferFunction(.py_add(p, 1L), .py_sub(5L, p), 4L))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_construction:669", function() DiscreteTransferFunction(3L, 4L, 8L))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_construction:671", function() DiscreteTransferFunction(.py_add(s, 1L), .py_add(.py_pow(s, 2L), 2L), s, 0L))
}))---")
  expect_identical(got[["lti:test_DiscreteTransferFunction_construction:567"]], "z + 1")
  expect_identical(got[["lti:test_DiscreteTransferFunction_construction:568"]], "z^2 + z + 1")
  expect_identical(got[["lti:test_DiscreteTransferFunction_construction:569"]], "(z + 1, z^2 + z + 1, z, 10)")
  expect_identical(got[["lti:test_DiscreteTransferFunction_construction:570"]], "10")
  expect_identical(got[["lti:test_DiscreteTransferFunction_construction:573"]], "z + 4")
  expect_identical(got[["lti:test_DiscreteTransferFunction_construction:574"]], "z - 5")
  expect_identical(got[["lti:test_DiscreteTransferFunction_construction:575"]], "0.100000000000000")
  expect_identical(got[["lti:test_DiscreteTransferFunction_construction:576"]], "True")
  expect_identical(got[["lti:test_DiscreteTransferFunction_construction:577"]], "(z + 4, z - 5, z, 0.100000000000000)")
  expect_identical(got[["lti:test_DiscreteTransferFunction_construction:581"]], "p + 3")
  expect_identical(got[["lti:test_DiscreteTransferFunction_construction:582"]], "p^2 - 9")
  expect_identical(got[["lti:test_DiscreteTransferFunction_construction:583"]], "12")
  expect_identical(got[["lti:test_DiscreteTransferFunction_construction:584"]], "(p + 3, p^2 - 9, p, 12)")
  expect_identical(got[["lti:test_DiscreteTransferFunction_construction:587"]], "(p^3 + 5*p^2 + 4, p^4 + 3*p + 1, p, 1)")
  expect_identical(got[["lti:test_DiscreteTransferFunction_construction:591"]], "(s - 1)*(s + 5)")
  expect_identical(got[["lti:test_DiscreteTransferFunction_construction:592"]], "((s - 1)*(s + 3), (s - 1)*(s + 5), s, 6.40000000000000)")
  expect_identical(got[["lti:test_DiscreteTransferFunction_construction:595"]], "(p + 2, p + 2, p, 1.20000000000000)")
  expect_identical(got[["lti:test_DiscreteTransferFunction_construction:598"]], "(s - 1, 4 - p, s, 1)")
  expect_identical(got[["lti:test_DiscreteTransferFunction_construction:601"]], "(s - 1, s - 1, s, 1)")
  expect_identical(got[["lti:test_DiscreteTransferFunction_construction:604"]], "5")
  expect_identical(got[["lti:test_DiscreteTransferFunction_construction:605"]], "6")
  expect_identical(got[["lti:test_DiscreteTransferFunction_construction:606"]], "(5, 6, s, 1)")
  expect_identical(got[["lti:test_DiscreteTransferFunction_construction:609"]], "0.500000000000000")
  expect_identical(got[["lti:test_DiscreteTransferFunction_construction:610"]], "4")
  expect_identical(got[["lti:test_DiscreteTransferFunction_construction:611"]], "(0.500000000000000, 4, s, 1)")
  expect_identical(got[["lti:test_DiscreteTransferFunction_construction:615"]], "DiscreteTransferFunction(2*p + 3*s^2 + 4*s, 8*p^2 + 7*s, s, 1)")
  expect_identical(got[["lti:test_DiscreteTransferFunction_construction:619"]], "DiscreteTransferFunction(a0*s + a1*s^2 + a2*s^3, b0*p - b1*s, s, 1)")
  expect_identical(got[["lti:test_DiscreteTransferFunction_construction:620"]], "True")
  expect_identical(got[["lti:test_DiscreteTransferFunction_construction:624"]], "DiscreteTransferFunction(a0*s + a1*s^2 + a2*s^3, b0*p - b1*s, s, 2)")
  expect_identical(got[["lti:test_DiscreteTransferFunction_construction:627"]], "(a*s^3 + b*s^2 + d + g*s, d*p + g*p^2 + g*s, s, 1)")
  expect_identical(got[["lti:test_DiscreteTransferFunction_construction:631"]], "(d + p^3, a + d*s + g*s^2, p, 2)")
  expect_identical(got[["lti:test_DiscreteTransferFunction_construction:632"]], "DiscreteTransferFunction(d + p^3, a + d*s + g*s^2, p, 2)")
  expect_identical(got[["lti:test_DiscreteTransferFunction_construction:635"]], "a0 + a1*s")
  expect_identical(got[["lti:test_DiscreteTransferFunction_construction:636"]], "b0 + b1*s + b2*s^2")
  expect_identical(got[["lti:test_DiscreteTransferFunction_construction:637"]], "(a0 + a1*s, b0 + b1*s + b2*s^2, s, 1)")
  expect_identical(got[["lti:test_DiscreteTransferFunction_construction:641"]], "(0, p^2 - p + 1, p, 1)")
  expect_identical(got[["lti:test_DiscreteTransferFunction_construction:644"]], "(0, 1, s, 1)")
  expect_identical(got[["lti:test_DiscreteTransferFunction_construction:648"]], "(a0*s^0.5 - a1 + a2*s^0.6, a1/p^8.7, s, 6)")
  expect_identical(got[["lti:test_DiscreteTransferFunction_construction:651"]], "(a1/s^0.8 + a2^2*p^0.25, a0*s - p, p, 1)")
  expect_identical(got[["lti:test_DiscreteTransferFunction_construction:656"]], "k_i/s + k_o*s + k_p")
  expect_identical(got[["lti:test_DiscreteTransferFunction_construction:657"]], "(k_i/s + k_o*s + k_p, omega_o^2 + 2*omega_o*s + s^2, s, 1)")
  expect_identical(got[["lti:test_DiscreteTransferFunction_construction:661"]], "ERR ValueError: TransferFunction cannot have a zero denominator.")
  expect_identical(got[["lti:test_DiscreteTransferFunction_construction:662"]], "ERR ValueError: TransferFunction cannot have a zero denominator.")
  expect_identical(got[["lti:test_DiscreteTransferFunction_construction:663"]], "ERR ValueError: TransferFunction cannot have a zero denominator.")
  expect_identical(got[["lti:test_DiscreteTransferFunction_construction:665"]], "ERR TypeError: Unsupported type for numerator or denominator of\\n                            TransferFunction.")
  expect_identical(got[["lti:test_DiscreteTransferFunction_construction:667"]], "ERR TypeError: Variable input must be a Symbol.")
  expect_identical(got[["lti:test_DiscreteTransferFunction_construction:668"]], "ERR TypeError: Variable input must be a Symbol.")
  expect_identical(got[["lti:test_DiscreteTransferFunction_construction:669"]], "ERR TypeError: Variable input must be a Symbol.")
  expect_identical(got[["lti:test_DiscreteTransferFunction_construction:671"]], "ERR ValueError: \\nThe sampling time cannot be zero. If you want to create a continuous\\ntransfer function, use the TransferFunction class instead.")
})

test_that("lti test_DiscreteTransferFunction_functions matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("a, x, b, c, s, g, d,    p, k, tau, zeta, wn, T, z"); a <- .py_getitem(.tmpu_, 0L); x <- .py_getitem(.tmpu_, 1L); b <- .py_getitem(.tmpu_, 2L); c <- .py_getitem(.tmpu_, 3L); s <- .py_getitem(.tmpu_, 4L); g <- .py_getitem(.tmpu_, 5L); d <- .py_getitem(.tmpu_, 6L); p <- .py_getitem(.tmpu_, 7L); k <- .py_getitem(.tmpu_, 8L); tau <- .py_getitem(.tmpu_, 9L); zeta <- .py_getitem(.tmpu_, 10L); wn <- .py_getitem(.tmpu_, 11L); T <- .py_getitem(.tmpu_, 12L); z <- .py_getitem(.tmpu_, 13L)
})
.py_try({
.tmpu_ <- .py_symbols("a0:4, b0:5, c0:4, d0:4"); a0 <- .py_getitem(.tmpu_, 0L); a1 <- .py_getitem(.tmpu_, 1L); a2 <- .py_getitem(.tmpu_, 2L); a3 <- .py_getitem(.tmpu_, 3L); b0 <- .py_getitem(.tmpu_, 4L); b1 <- .py_getitem(.tmpu_, 5L); b2 <- .py_getitem(.tmpu_, 6L); b3 <- .py_getitem(.tmpu_, 7L); b4 <- .py_getitem(.tmpu_, 8L); c0 <- .py_getitem(.tmpu_, 9L); c1 <- .py_getitem(.tmpu_, 10L); c2 <- .py_getitem(.tmpu_, 11L); c3 <- .py_getitem(.tmpu_, 12L); d0 <- .py_getitem(.tmpu_, 13L); d1 <- .py_getitem(.tmpu_, 14L); d2 <- .py_getitem(.tmpu_, 15L); d3 <- .py_getitem(.tmpu_, 16L)
})
.py_try({
TF1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
TF2 <- TransferFunction(k, 1L, s)
})
.py_try({
TF3 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
TF4 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(tau, s), 1L), s)
})
.py_try({
TF5 <- TransferFunction(k, p, s)
})
.py_try({
TF6 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(k, s), p), s)
})
p_are_floats_equals <- function(n1, n2) {
  "\n    Helper function to compare two floats.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return(.py_method(Float(n1), "is_same", Float(n2), .py_isclose))
}
p_are_complex_equals <- function(n1, n2) {
  "\n    Helper function to compare two complex numbers.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return((.py_truthy(p_are_floats_equals(re(n1), re(n2))) && .py_truthy(p_are_floats_equals(im(n1), im(n2)))))
}
.py_try({
expr_1 <- .py_Mul(0L, .py_Pow(s, -1L, evaluate = FALSE), evaluate = FALSE)
})
.py_try({
expr_2 <- .py_div(s, 0L)
})
.py_try({
expr_3 <- .py_div(.py_add(.py_mul(p, .py_pow(s, 2L)), .py_mul(5L, s)), .py_pow(.py_add(s, 1L), 3L))
})
.py_try({
expr_4 <- 6L
})
.py_try({
expr_5 <- .py_div(.py_mul(.py_add(2L, .py_mul(3L, s)), .py_add(5L, .py_mul(2L, s))), .py_mul(.py_add(9L, .py_mul(3L, s)), .py_add(5L, .py_mul(2L, .py_pow(s, 2L)))))
})
.py_try({
expr_6 <- .py_div(.py_add(.py_add(.py_mul(9L, .py_pow(s, 4L)), .py_mul(4L, .py_pow(s, 2L))), 8L), .py_mul(.py_add(s, 1L), .py_add(s, 9L)))
})
.py_try({
tf <- DiscreteTransferFunction(.py_add(s, 1L), .py_add(.py_pow(s, 2L), 2L), s)
})
.py_try({
delay <- .py_pow(s, .py_div(-1L, tau))
})
.py_try({
expr_7 <- .py_mul(delay, .py_method(tf, "to_expr"))
})
.py_try({
H1 <- .py_method(DiscreteTransferFunction, "from_rational_expression", expr_7, s)
})
.py_try({
H2 <- DiscreteTransferFunction(.py_add(s, 1L), .py_mul(.py_pow(s, .py_div(1L, tau)), .py_add(.py_pow(s, 2L), 2L)), s)
})
.py_try({
expr_8 <- .py_Add(2L, .py_div(.py_mul(3L, s), .py_add(.py_pow(s, 2L), 1L)), evaluate = FALSE)
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:705", function() .py_method(DiscreteTransferFunction, "from_rational_expression", expr_1, sampling_time = 12L))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:708", function() .py_method(DiscreteTransferFunction, "from_rational_expression", expr_2))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:710", function() .py_method(DiscreteTransferFunction, "from_rational_expression", expr_3))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:712", function() .py_method(DiscreteTransferFunction, "from_rational_expression", expr_3, s))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:714", function() .py_method(DiscreteTransferFunction, "from_rational_expression", expr_3, p, sampling_time = .py_flt(1.4)))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:717", function() .py_method(DiscreteTransferFunction, "from_rational_expression", expr_4))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:719", function() .py_method(DiscreteTransferFunction, "from_rational_expression", expr_4, s))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:721", function() .py_method(DiscreteTransferFunction, "from_rational_expression", expr_5, s))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:723", function() .py_method(DiscreteTransferFunction, "from_rational_expression", expr_6, s))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:725", function() H1)
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:726", function() .py_method(DiscreteTransferFunction, "from_rational_expression", expr_8, s, sampling_time = .py_flt(0.3)))
})
.py_try({
tf1 <- .py_method(DiscreteTransferFunction, "from_coeff_lists", list(1L, 2L), list(3L, 4L, 5L), s)
})
.py_try({
num2 <- list(.py_pow(p, 2L), .py_mul(2L, p))
})
.py_try({
den2 <- list(.py_pow(p, 3L), .py_add(p, 1L), 4L)
})
.py_try({
tf2 <- .py_method(DiscreteTransferFunction, "from_coeff_lists", num2, den2, s)
})
.py_try({
num3 <- list(1L, 2L, 3L)
})
.py_try({
den3 <- list(0L, 0L)
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:738", function() tf1)
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:739", function() tf2)
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:740", function() .py_method(DiscreteTransferFunction, "from_coeff_lists", num3, den3, s))
})
.py_try({
zeros <- list(4L)
})
.py_try({
poles <- list(.py_add(-1L, complex(real = 0.0, imaginary = 2.0)), .py_sub(-1L, complex(real = 0.0, imaginary = 2.0)))
})
.py_try({
gain <- 3L
})
.py_try({
tf1 <- .py_method(DiscreteTransferFunction, "from_zpk", zeros, poles, gain, s, 3L)
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:749", function() tf1)
})
.py_try({
tf0 <- DiscreteTransferFunction(.py_add(.py_add(.py_pow(s, 5L), .py_pow(s, 3L)), s), .py_sub(s, .py_pow(s, 2L)), s)
})
.py_try({
a <- DiscreteTransferFunction(.py_neg(.py_add(.py_add(.py_pow(s, 4L), .py_pow(s, 2L)), 1L)), .py_sub(s, 1L), s)
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:755", function() (.py_truthy(.py_eq(.py_method(tf0, "simplify"), simplify(tf0))) && .py_truthy(.py_eq(simplify(tf0), a))))
})
.py_try({
tf1 <- DiscreteTransferFunction(.py_mul(.py_add(p, 3L), .py_sub(p, 1L)), .py_mul(.py_sub(p, 1L), .py_add(p, 5L)), p)
})
.py_try({
b <- DiscreteTransferFunction(.py_add(p, 3L), .py_add(p, 5L), p)
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:759", function() (.py_truthy(.py_eq(.py_method(tf1, "simplify"), simplify(tf1))) && .py_truthy(.py_eq(simplify(tf1), b))))
})
.py_try({
G1 <- DiscreteTransferFunction(.py_pow(.py_sub(1L, s), 2L), .py_pow(.py_add(.py_pow(s, 2L), 1L), 2L), s)
})
.py_try({
G2 <- DiscreteTransferFunction(1L, -3L, p)
})
.py_try({
c <- .py_mul(.py_add(.py_add(.py_mul(a2, .py_pow(s, p)), .py_mul(a1, .py_pow(s, s))), .py_mul(a0, .py_pow(p, p))), .py_add(.py_pow(p, s), .py_pow(s, p)))
})
.py_try({
d <- .py_mul(.py_add(.py_mul(b0, .py_pow(s, s)), .py_mul(b1, .py_pow(p, s))), .py_add(.py_mul(.py_mul(b2, s), p), .py_pow(p, p)))
})
.py_try({
e <- .py_add(.py_add(.py_add(.py_add(.py_add(.py_mul(.py_mul(a0, .py_pow(p, p)), .py_pow(p, s)), .py_mul(.py_mul(a0, .py_pow(p, p)), .py_pow(s, p))), .py_mul(.py_mul(a1, .py_pow(p, s)), .py_pow(s, s))), .py_mul(.py_mul(a1, .py_pow(s, p)), .py_pow(s, s))), .py_mul(.py_mul(a2, .py_pow(p, s)), .py_pow(s, p))), .py_mul(a2, .py_pow(s, .py_mul(2L, p))))
})
.py_try({
f <- .py_add(.py_add(.py_add(.py_mul(.py_mul(.py_mul(.py_mul(b0, b2), p), s), .py_pow(s, s)), .py_mul(.py_mul(b0, .py_pow(p, p)), .py_pow(s, s))), .py_mul(.py_mul(.py_mul(.py_mul(b1, b2), p), .py_pow(p, s)), s)), .py_mul(.py_mul(b1, .py_pow(p, p)), .py_pow(p, s)))
})
.py_try({
g <- .py_add(.py_add(.py_add(.py_mul(.py_mul(.py_mul(a1, a2), s), .py_pow(s, p)), .py_mul(.py_mul(a1, p), s)), .py_mul(.py_mul(.py_mul(.py_mul(a2, b1), p), s), .py_pow(s, p))), .py_mul(.py_mul(b1, .py_pow(p, 2L)), s))
})
.py_try({
G3 <- DiscreteTransferFunction(c, d, s)
})
.py_try({
G4 <- DiscreteTransferFunction(.py_sub(.py_mul(a0, .py_pow(s, s)), .py_mul(b0, .py_pow(p, p))), .py_mul(.py_add(.py_mul(a1, s), .py_mul(.py_mul(b1, s), p)), .py_add(.py_mul(a2, .py_pow(s, p)), p)), p, .py_flt(2.4))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:774", function() .py_method(G1, "expand"))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:776", function() .py_method(tf1, "expand"))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:777", function() .py_method(G2, "expand"))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:778", function() .py_method(G3, "expand"))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:779", function() .py_method(G4, "expand"))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:780", function() .py_method(G4, "expand"))
})
.py_try({
p1 <- .py_add(.py_mul(a1, s), a0)
})
.py_try({
p2 <- .py_add(.py_add(.py_mul(b2, .py_pow(s, 2L)), .py_mul(b1, s)), b0)
})
.py_try({
SP1 <- DiscreteTransferFunction(p1, p2, s)
})
.py_try({
expect1 <- DiscreteTransferFunction(.py_add(.py_mul(.py_flt(2.0), s), .py_flt(1.0)), .py_add(.py_add(.py_mul(.py_flt(5.0), .py_pow(s, 2L)), .py_mul(.py_flt(4.0), s)), .py_flt(3.0)), s)
})
.py_try({
expect1_ <- DiscreteTransferFunction(.py_add(.py_mul(2L, s), 1L), .py_add(.py_add(.py_mul(5L, .py_pow(s, 2L)), .py_mul(4L, s)), 3L), s)
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:788", function() .py_method(SP1, "subs", .py_dict(list(a0, a1, b0, b1, b2), list(1L, 2L, 3L, 4L, 5L))))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:789", function() .py_method(.py_method(SP1, "subs", .py_dict(list(a0, a1, b0, b1, b2), list(1L, 2L, 3L, 4L, 5L))), "evalf"))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:790", function() .py_method(expect1_, "evalf"))
})
.py_try({
.tmpu_ <- .py_symbols("c1, d0:3"); c1 <- .py_getitem(.tmpu_, 0L); d0 <- .py_getitem(.tmpu_, 1L); d1 <- .py_getitem(.tmpu_, 2L); d2 <- .py_getitem(.tmpu_, 3L)
})
.py_try({
.tmpu_ <- .py_tup(.py_mul(c1, p), .py_sub(.py_add(.py_mul(d2, .py_pow(p, 3L)), .py_mul(d1, .py_pow(p, 2L))), d0)); p3 <- .py_getitem(.tmpu_, 0L); p4 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
SP2 <- DiscreteTransferFunction(p3, p4, p)
})
.py_try({
expect2 <- DiscreteTransferFunction(.py_mul(.py_flt(2.0), p), .py_sub(.py_add(.py_mul(.py_flt(5.0), .py_pow(p, 3L)), .py_mul(.py_flt(2.0), .py_pow(p, 2L))), .py_flt(3.0)), p)
})
.py_try({
expect2_ <- DiscreteTransferFunction(.py_mul(2L, p), .py_sub(.py_add(.py_mul(5L, .py_pow(p, 3L)), .py_mul(2L, .py_pow(p, 2L))), 3L), p)
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:797", function() .py_method(SP2, "subs", .py_dict(list(c1, d0, d1, d2), list(2L, 3L, 2L, 5L))))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:798", function() .py_method(.py_method(SP2, "subs", .py_dict(list(c1, d0, d1, d2), list(2L, 3L, 2L, 5L))), "evalf"))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:799", function() .py_method(expect2_, "evalf"))
})
.py_try({
SP3 <- DiscreteTransferFunction(.py_add(.py_sub(.py_add(.py_mul(a0, .py_pow(p, 3L)), .py_mul(a1, .py_pow(s, 2L))), .py_mul(b0, s)), b1), .py_add(.py_mul(a1, s), p), s)
})
.py_try({
expect3 <- DiscreteTransferFunction(.py_add(.py_sub(.py_add(.py_mul(.py_flt(2.0), .py_pow(p, 3L)), .py_mul(.py_flt(4.0), .py_pow(s, 2L))), s), .py_flt(5.0)), .py_add(p, .py_mul(.py_flt(4.0), s)), s)
})
.py_try({
expect3_ <- DiscreteTransferFunction(.py_add(.py_sub(.py_add(.py_mul(2L, .py_pow(p, 3L)), .py_mul(4L, .py_pow(s, 2L))), s), 5L), .py_add(p, .py_mul(4L, s)), s)
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:804", function() .py_method(SP3, "subs", .py_dict(list(a0, a1, b0, b1), list(2L, 4L, 1L, 5L))))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:805", function() .py_method(.py_method(SP3, "subs", .py_dict(list(a0, a1, b0, b1), list(2L, 4L, 1L, 5L))), "evalf"))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:806", function() .py_method(expect3_, "evalf"))
})
.py_try({
SP4 <- DiscreteTransferFunction(.py_sub(s, .py_mul(a1, .py_pow(p, 3L))), .py_add(.py_mul(a0, s), p), p)
})
.py_try({
expect4 <- DiscreteTransferFunction(.py_add(.py_mul(.py_flt(7.0), .py_pow(p, 3L)), s), .py_sub(p, s), p)
})
.py_try({
expect4_ <- DiscreteTransferFunction(.py_add(.py_mul(7L, .py_pow(p, 3L)), s), .py_sub(p, s), p)
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:811", function() .py_method(SP4, "subs", .py_dict(list(a0, a1), list(-1L, -7L))))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:812", function() .py_method(.py_method(SP4, "subs", .py_dict(list(a0, a1), list(-1L, -7L))), "evalf"))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:813", function() .py_method(expect4_, "evalf"))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:816", function() .py_method(tf1, "eval_frequency", wn))
})
.py_try({
dtf_fr <- DiscreteTransferFunction(.py_sub(z, 1L), .py_add(.py_add(.py_pow(z, 2L), z), 1L), z)
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:818", function() .py_method(dtf_fr, "eval_frequency", .py_add(1L, sym_I)))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:819", function() .py_method(G4, "eval_frequency", .py_div(.py_S(5L), 3L)))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:827", function() .py_method(DiscreteTransferFunction(z, .py_sub(z, 1L), z), "dc_gain"))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:828", function() .py_method(DiscreteTransferFunction(.py_sub(z, 1L), .py_add(.py_add(.py_pow(z, 2L), .py_mul(3L, z)), 1L), z), "dc_gain"))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:829", function() .py_method(DiscreteTransferFunction(z, .py_add(z, 2L), z, .py_flt(0.2)), "dc_gain"))
})
.py_try({
tf_ <- DiscreteTransferFunction(.py_sub(.py_pow(x, 3L), k), k, x, 3L)
})
.py_try({
p_tf <- DiscreteTransferFunction(k, .py_sub(.py_pow(x, 4L), k), x, .py_flt(1.6))
})
.py_try({
TF_ <- DiscreteTransferFunction(.py_pow(x, 2L), .py_add(.py_add(.py_pow(x, 10L), x), .py_pow(x, 2L)), x)
})
.py_try({
p_TF <- DiscreteTransferFunction(.py_add(.py_add(.py_pow(x, 10L), x), .py_pow(x, 2L)), .py_pow(x, 2L), x)
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:836", function() .py_method(G1, "poles"))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:837", function() .py_method(G2, "poles"))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:838", function() .py_method(tf1, "poles"))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:839", function() .py_method(expect4_, "poles"))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:840", function() .py_method(SP4, "poles"))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:841", function() .py_method(expect3, "poles"))
})
.py_try({
exp_poles1 <- list(.py_sub(.py_flt(-0.4), .py_mul(.py_flt(0.66332495807108), sym_I)), .py_add(.py_flt(-0.4), .py_mul(.py_flt(0.66332495807108), sym_I)))
})
.py_try({
for (.it_ in .py_iter(.py_cap(.py_enumerate(.py_method(expect1, "poles")), 60L))) {
  i <- .py_getitem(.it_, 0L); pol <- .py_getitem(.it_, 1L)
  .py_try({
  .py_rec("lti:test_DiscreteTransferFunction_functions:845", function() p_are_complex_equals(pol, .py_getitem(exp_poles1, i)))
  })
}
})
.py_try({
exp_poles2 <- list(.py_flt(0.729001428685125), .py_sub(.py_flt(-0.564500714342563), .py_mul(.py_flt(0.710198984796332), sym_I)), .py_add(.py_flt(-0.564500714342563), .py_mul(.py_flt(0.710198984796332), sym_I)))
})
.py_try({
for (.it_ in .py_iter(.py_cap(.py_enumerate(.py_method(expect2, "poles")), 60L))) {
  i <- .py_getitem(.it_, 0L); pol <- .py_getitem(.it_, 1L)
  .py_try({
  .py_rec("lti:test_DiscreteTransferFunction_functions:850", function() p_are_complex_equals(pol, .py_getitem(exp_poles2, i)))
  })
}
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:852", function() .py_method(p_tf, "poles"))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:854", function() .py_method(TF_, "poles"))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:860", function() .py_method(DiscreteTransferFunction(1L, .py_add(.py_add(.py_mul(b0, .py_pow(x, 2L)), .py_mul(b1, x)), b2), x), "poles"))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:863", function() .py_method(DiscreteTransferFunction(.py_pow(x, 2L), .py_add(.py_add(.py_mul(a0, .py_pow(x, 10L)), x), .py_pow(x, 2L)), x), "poles"))
})
.py_try({
stable_tf <- DiscreteTransferFunction(z, .py_mul(.py_sub(z, Rational(1L, 2L)), .py_add(z, Rational(9L, 10L))), z, .py_flt(0.1))
})
.py_try({
marginally_stable_tf <- DiscreteTransferFunction(z, .py_mul(.py_sub(z, Rational(1L, 4L)), .py_add(z, 1L)), z)
})
.py_try({
unstable_tf <- DiscreteTransferFunction(z, .py_mul(.py_sub(z, Rational(1L, 4L)), .py_add(z, 2L)), z, T)
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:875", function() .py_method(stable_tf, "is_stable"))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:876", function() .py_method(marginally_stable_tf, "is_stable"))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:877", function() .py_method(unstable_tf, "is_stable"))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:878", function() .py_method(DiscreteTransferFunction(z, .py_mul(.py_add(z, Rational(1L, 8L)), .py_add(z, k)), z), "is_stable"))
})
.py_try({
generic_den <- .py_add(.py_add(.py_add(.py_add(.py_mul(b4, .py_pow(s, 4L)), .py_mul(b3, .py_pow(s, 3L))), .py_mul(b2, .py_pow(s, 2L))), .py_mul(b1, s)), b0)
})
.py_try({
stab_cond <- .py_method(DiscreteTransferFunction(1L, generic_den, s), "get_asymptotic_stability_conditions")
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:885", function() stab_cond)
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:903", function() .py_method(G1, "zeros"))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:904", function() .py_method(G2, "zeros"))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:905", function() .py_method(tf1, "zeros"))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:906", function() .py_method(expect4_, "zeros"))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:911", function() .py_method(SP4, "zeros"))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:921", function() .py_method(tf_, "zeros"))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:925", function() .py_method(p_TF, "zeros"))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:931", function() .py_method(DiscreteTransferFunction(.py_add(.py_add(.py_mul(a0, .py_pow(x, 10L)), x), .py_pow(x, 2L)), .py_pow(x, 2L), x), "zeros"))
})
.py_try({
tf2 <- DiscreteTransferFunction(.py_add(s, 3L), .py_add(.py_sub(.py_pow(s, 2L), .py_pow(s, 3L)), 9L), s)
})
.py_try({
tf3 <- DiscreteTransferFunction(.py_add(.py_mul(-3L, p), 3L), .py_sub(1L, p), p)
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:937", function() .py_neg(tf2))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:938", function() .py_neg(tf3))
})
.py_try({
tf4 <- DiscreteTransferFunction(.py_add(p, 4L), .py_sub(p, 3L), p)
})
.py_try({
tf5 <- DiscreteTransferFunction(.py_add(.py_pow(s, 2L), 1L), .py_sub(1L, s), s)
})
.py_try({
expect2 <- DiscreteTransferFunction(.py_pow(.py_add(.py_pow(s, 2L), 1L), 3L), .py_pow(.py_sub(1L, s), 3L), s)
})
.py_try({
expect1 <- DiscreteTransferFunction(.py_pow(.py_add(p, 4L), 2L), .py_pow(.py_sub(p, 3L), 2L), p)
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:945", function() (.py_truthy(.py_eq(.py_method(.py_mul(tf4, tf4), "doit"), .py_pow(tf4, 2L))) && .py_truthy(.py_eq(.py_pow(tf4, 2L), pow(tf4, 2L))) && .py_truthy(.py_eq(pow(tf4, 2L), expect1))))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:946", function() (.py_truthy(.py_eq(.py_method(.py_mul(.py_mul(tf5, tf5), tf5), "doit"), .py_pow(tf5, 3L))) && .py_truthy(.py_eq(.py_pow(tf5, 3L), pow(tf5, 3L))) && .py_truthy(.py_eq(pow(tf5, 3L), expect2))))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:947", function() (.py_truthy(.py_eq(.py_pow(tf5, 0L), pow(tf5, 0L))) && .py_truthy(.py_eq(pow(tf5, 0L), DiscreteTransferFunction(1L, 1L, s)))))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:948", function() (.py_truthy(.py_eq(.py_pow(.py_method(Series(tf4), "doit"), -1L), .py_pow(tf4, -1L))) && .py_truthy(.py_eq(.py_pow(tf4, -1L), pow(tf4, -1L))) && .py_truthy(.py_eq(pow(tf4, -1L), DiscreteTransferFunction(.py_sub(p, 3L), .py_add(p, 4L), p)))))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:950", function() (.py_truthy(.py_eq(.py_pow(.py_method(.py_mul(tf5, tf5), "doit"), -1L), .py_pow(tf5, -2L))) && .py_truthy(.py_eq(.py_pow(tf5, -2L), pow(tf5, -2L))) && .py_truthy(.py_eq(pow(tf5, -2L), DiscreteTransferFunction(.py_pow(.py_sub(1L, s), 2L), .py_pow(.py_add(.py_pow(s, 2L), 1L), 2L), s)))))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:953", function() .py_pow(tf4, .py_sub(.py_add(.py_pow(s, 2L), s), 1L)))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:954", function() .py_pow(tf5, s))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:955", function() .py_pow(tf4, tf5))
})
.py_try({
tf <- DiscreteTransferFunction(.py_sub(s, 1L), .py_add(.py_sub(.py_pow(s, 2L), .py_mul(2L, s)), 1L), s)
})
.py_try({
tf6 <- DiscreteTransferFunction(.py_add(s, p), .py_sub(.py_pow(p, 2L), 5L), s)
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:960", function() sym_factor(tf))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:961", function() (.py_truthy(.py_eq(.py_method(.py_attr(tf, "num"), "subs", s, 2L), .py_method(.py_attr(tf, "den"), "subs", s, 2L))) && .py_truthy(.py_eq(.py_method(.py_attr(tf, "den"), "subs", s, 2L), 1L))))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:963", function() .py_method(tf, "subs", s, 2L))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:964", function() .py_method(tf6, "subs", p, 3L))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:965", function() .py_method(tf3, "xreplace", .py_dict(list(p), list(s))))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:966", function() .py_method(tf3, "xreplace", .py_dict(list(p), list(.py_exp(2L)))))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:967", function() .py_method(tf3, "subs", p, .py_exp(2L)))
})
.py_try({
tf7 <- DiscreteTransferFunction(.py_add(.py_mul(a0, .py_pow(s, p)), .py_mul(a1, .py_pow(p, s))), .py_sub(.py_mul(a2, p), s), s)
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:970", function() .py_method(tf7, "xreplace", .py_dict(list(s), list(k))))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:972", function() .py_method(tf7, "subs", s, k))
})
.py_try({
tf8 <- DiscreteTransferFunction(.py_add(.py_add(.py_mul(a0, .py_pow(s, 5L)), .py_mul(5L, .py_pow(s, 2L))), 3L), .py_sub(.py_pow(s, 6L), 3L), s)
})
.py_try({
tf9 <- DiscreteTransferFunction(.py_add(5L, s), .py_mul(.py_add(5L, s), .py_add(6L, s)), s)
})
.py_try({
tf10 <- DiscreteTransferFunction(0L, 1L, s)
})
.py_try({
tf11 <- DiscreteTransferFunction(1L, 1L, s)
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:979", function() .py_method(tf8, "to_expr"))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:982", function() .py_method(tf9, "to_expr"))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:985", function() .py_method(tf10, "to_expr"))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_functions:987", function() .py_method(tf11, "to_expr"))
}))---")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:705"]], "DiscreteTransferFunction(0, s, s, 12)")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:708"]], "ERR ZeroDivisionError: TransferFunction cannot have a zero denominator.")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:710"]], "ERR ValueError: \\nConflicting values found for positional argument `var` ({p, s}).\\nSpecify it manually.")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:712"]], "DiscreteTransferFunction(p*s^2 + 5*s, (s + 1)^3, s, 1)")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:714"]], "DiscreteTransferFunction(p*s^2 + 5*s, (s + 1)^3, p, 1.4)")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:717"]], "ERR ValueError: \\nPositional argument `var` not found in the TransferFunction defined.\\nSpecify it manually.")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:719"]], "DiscreteTransferFunction(6, 1, s, 1)")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:721"]], "DiscreteTransferFunction((2*s + 5)*(3*s + 2), (3*s + 9)*(2*s^2 + 5), s, 1)")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:723"]], "DiscreteTransferFunction(9*s^4 + 4*s^2 + 8, (s + 1)*(s + 9), s, 1)")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:725"]], "DiscreteTransferFunction(s + 1, s^(1/tau)*(s^2 + 2), s, 1)")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:726"]], "DiscreteTransferFunction(2*s^2 + 3*s + 2, s^2 + 1, s, 0.3)")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:738"]], "DiscreteTransferFunction(s + 2, 3*s^2 + 4*s + 5, s, 1)")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:739"]], "DiscreteTransferFunction(p^2*s + 2*p, p^3*s^2 + s*(p + 1) + 4, s, 1)")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:740"]], "ERR ZeroDivisionError: TransferFunction cannot have a zero denominator.")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:749"]], "DiscreteTransferFunction(3*s - 12, (s + 1.0 - 2.0*I)*(s + 1.0 + 2.0*I), s, 3)")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:755"]], "True")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:759"]], "True")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:774"]], "DiscreteTransferFunction(s^2 - 2*s + 1, s^4 + 2*s^2 + 1, s, 1)")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:776"]], "DiscreteTransferFunction(p^2 + 2*p - 3, p^2 + 4*p - 5, p, 1)")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:777"]], "DiscreteTransferFunction(1, -3, p, 1)")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:778"]], "DiscreteTransferFunction(a0*p^p*p^s + a0*p^p*s^p + a1*p^s*s^s + a1*s^p*s^s + a2*p^s*s^p + a2*s^(2*p), b0*b2*p*s*s^s + b0*p^p*s^s + b1*b2*p*p^s*s + b1*p^p*p^s, s, 1)")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:779"]], "DiscreteTransferFunction(a0*s^s - b0*p^p, a1*a2*s*s^p + a1*p*s + a2*b1*p*s*s^p + b1*p^2*s, p, 2.4)")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:780"]], "DiscreteTransferFunction(a0*s^s - b0*p^p, a1*a2*s*s^p + a1*p*s + a2*b1*p*s*s^p + b1*p^2*s, p, 2.4)")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:788"]], "DiscreteTransferFunction(2*s + 1, 5*s^2 + 4*s + 3, s, 1)")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:789"]], "DiscreteTransferFunction(2.0*s + 1.0, 5.0*s^2 + 4.0*s + 3.0, s, 1)")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:790"]], "DiscreteTransferFunction(2.0*s + 1.0, 5.0*s^2 + 4.0*s + 3.0, s, 1)")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:797"]], "DiscreteTransferFunction(2*p, 5*p^3 + 2*p^2 - 3, p, 1)")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:798"]], "DiscreteTransferFunction(2.0*p, 5.0*p^3 + 2.0*p^2 - 3.0, p, 1)")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:799"]], "DiscreteTransferFunction(2.0*p, 5.0*p^3 + 2.0*p^2 - 3.0, p, 1)")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:804"]], "DiscreteTransferFunction(2*p^3 + 4*s^2 - s + 5, p + 4*s, s, 1)")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:805"]], "DiscreteTransferFunction(2.0*p^3 + 4.0*s^2 - s + 5.0, p + 4.0*s, s, 1)")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:806"]], "DiscreteTransferFunction(2.0*p^3 + 4.0*s^2 - s + 5.0, p + 4.0*s, s, 1)")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:811"]], "DiscreteTransferFunction(7*p^3 + s, p - s, p, 1)")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:812"]], "DiscreteTransferFunction(7.0*p^3 + s, p - s, p, 1)")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:813"]], "DiscreteTransferFunction(7.0*p^3 + s, p - s, p, 1)")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:816"]], "wn/(wn + 5) + 3/(wn + 5)")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:818"]], "3/13 + 2*I/13")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:819"]], "a0*s^s/(a1*a2*s^(8/3) + 5*a1*s/3 + 5*a2*b1*s^(8/3)/3 + 25*b1*s/9) - 5*3^(1/3)*5^(2/3)*b0/(9*a1*a2*s^(8/3) + 15*a1*s + 15*a2*b1*s^(8/3) + 25*b1*s)")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:827"]], "oo")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:828"]], "0")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:829"]], "1/3")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:836"]], "[I, I, -I, -I]")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:837"]], "[]")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:838"]], "[-5, 1]")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:839"]], "[s]")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:840"]], "[-a0*s]")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:841"]], "[-0.25*p]")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:845"]], "True")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:845#2"]], "True")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:850"]], "True")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:850#2"]], "True")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:850#3"]], "True")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:852"]], "[k^(1/4), -k^(1/4), I*k^(1/4), -I*k^(1/4)]")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:854"]], "[CRootOf(x^9 + x + 1, 0), 0, CRootOf(x^9 + x + 1, 1), CRootOf(x^9 + x + 1, 2), CRootOf(x^9 + x + 1, 3), CRootOf(x^9 + x + 1, 4), CRootOf(x^9 + x + 1, 5), CRootOf(x^9 + x + 1, 6), CRootOf(x^9 + x + 1, 7), CRootOf(x^9 + x + 1, 8)]")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:860"]], "[-b1/(2*b0) - sqrt(-4*b0*b2 + b1^2)/(2*b0), -b1/(2*b0) + sqrt(-4*b0*b2 + b1^2)/(2*b0)]")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:863"]], "ERR NotImplementedError: CRootOf is not supported over ZZ[a0]")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:875"]], "True")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:876"]], "False")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:877"]], "False")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:878"]], "None")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:885"]], "[-4*b0^2 + 6*b0*b1 - 4*b0*b2 + 2*b0*b3 - 2*b1^2 + 2*b1*b2 - 2*b1*b4 - 2*b2*b3 + 4*b2*b4 + 2*b3^2 - 6*b3*b4 + 4*b4^2 > 0, -20*b0^2 + 10*b0*b1 + 12*b0*b2 - 18*b0*b3 - 2*b1^2 - 2*b1*b2 + 18*b1*b4 + 2*b2*b3 - 12*b2*b4 + 2*b3^2 - 10*b3*b4 + 20*b4^2 > 0, 64*b0^4 - 64*b0^3*b1 - 64*b0^3*b3 + 64*b0^2*b1*b2 + 64*b0^2*b1*b3 + 64*b0^2*b1*b4 - 64*b0^2*b2^2 + 64*b0^2*b2*b3 - 64*b0^2*b3^2 + 64*b0^2*b3*b4 - 128*b0^2*b4^2 - 64*b0*b1^2*b3 - 64*b0*b1^2*b4 + 64*b0*b1*b2*b3 - 128*b0*b1*b2*b4 + 128*b0*b1*b3*b4 + 64*b0*b1*b4^2 + 128*b0*b2^2*b4 - 64*b0*b2*b3^2 - 128*b0*b2*b3*b4 + 64*b0*b3^3 - 64*b0*b3^2*b4 + 64*b0*b3*b4^2 + 64*b1^3*b4 - 64*b1^2*b2*b4 - 64*b1^2*b4^2 + 64*b1*b2*b3*b4 + 64*b1*b2*b4^2 - 64*b1*b3^2*b4 + 64*b1*b3*b4^2 - 64*b1*b4^3 - 64*b2^2*b4^2 + 64*b2*b3*b4^2 - 64*b3*b4^3 + 64*b4^4 > 0, b0^2 + 2*b0*b2 + 2*b0*b4 - b1^2 - 2*b1*b3 + b2^2 + 2*b2*b4 - b3^2 + b4^2 > 0]")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:903"]], "[1, 1]")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:904"]], "[]")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:905"]], "[-3, 1]")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:906"]], "[-7^(2/3)*s^(1/3)/7, 7^(2/3)*s^(1/3)/14 - sqrt(3)*7^(2/3)*I*s^(1/3)/14, 7^(2/3)*s^(1/3)/14 + sqrt(3)*7^(2/3)*I*s^(1/3)/14]")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:911"]], "[s^(1/3)/a1^(1/3), -s^(1/3)/(2*a1^(1/3)) - sqrt(3)*I*s^(1/3)/(2*a1^(1/3)), -s^(1/3)/(2*a1^(1/3)) + sqrt(3)*I*s^(1/3)/(2*a1^(1/3))]")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:921"]], "[k^(1/3), -k^(1/3)/2 - sqrt(3)*I*k^(1/3)/2, -k^(1/3)/2 + sqrt(3)*I*k^(1/3)/2]")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:925"]], "[CRootOf(x^9 + x + 1, 0), 0, CRootOf(x^9 + x + 1, 1), CRootOf(x^9 + x + 1, 2), CRootOf(x^9 + x + 1, 3), CRootOf(x^9 + x + 1, 4), CRootOf(x^9 + x + 1, 5), CRootOf(x^9 + x + 1, 6), CRootOf(x^9 + x + 1, 7), CRootOf(x^9 + x + 1, 8)]")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:931"]], "ERR NotImplementedError: CRootOf is not supported over ZZ[a0]")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:937"]], "DiscreteTransferFunction(-s - 3, -s^3 + s^2 + 9, s, 1)")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:938"]], "DiscreteTransferFunction(3*p - 3, 1 - p, p, 1)")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:945"]], "True")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:946"]], "True")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:947"]], "True")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:948"]], "True")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:950"]], "True")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:953"]], "ERR ValueError: Exponent must be an integer.")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:954"]], "ERR ValueError: Exponent must be an integer.")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:955"]], "ERR ValueError: Exponent must be an integer.")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:960"]], "DiscreteTransferFunction(s - 1, (s - 1)^2, s, 1)")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:961"]], "True")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:963"]], "DiscreteTransferFunction(s - 1, s^2 - 2*s + 1, s, 1)")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:964"]], "DiscreteTransferFunction(s + 3, 4, s, 1)")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:965"]], "DiscreteTransferFunction(3 - 3*s, 1 - s, s, 1)")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:966"]], "ERR TypeError: Variable input must be a Symbol.")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:967"]], "DiscreteTransferFunction(3 - 3*p, 1 - p, p, 1)")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:970"]], "DiscreteTransferFunction(a0*k^p + a1*p^k, a2*p - k, k, 1)")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:972"]], "DiscreteTransferFunction(a0*s^p + a1*p^s, a2*p - s, s, 1)")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:979"]], "(a0*s^5 + 5*s^2 + 3)/(s^6 - 3)")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:982"]], "(s + 5)/(((s + 5)*(s + 6)))")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:985"]], "0/1")
  expect_identical(got[["lti:test_DiscreteTransferFunction_functions:987"]], "1/1")
})

test_that("lti test_DiscreteTransferFunction_addition_and_subtraction matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("a, x, b, c, s, g, d,    p, k, tau, zeta, wn, T, z"); a <- .py_getitem(.tmpu_, 0L); x <- .py_getitem(.tmpu_, 1L); b <- .py_getitem(.tmpu_, 2L); c <- .py_getitem(.tmpu_, 3L); s <- .py_getitem(.tmpu_, 4L); g <- .py_getitem(.tmpu_, 5L); d <- .py_getitem(.tmpu_, 6L); p <- .py_getitem(.tmpu_, 7L); k <- .py_getitem(.tmpu_, 8L); tau <- .py_getitem(.tmpu_, 9L); zeta <- .py_getitem(.tmpu_, 10L); wn <- .py_getitem(.tmpu_, 11L); T <- .py_getitem(.tmpu_, 12L); z <- .py_getitem(.tmpu_, 13L)
})
.py_try({
.tmpu_ <- .py_symbols("a0:4, b0:5, c0:4, d0:4"); a0 <- .py_getitem(.tmpu_, 0L); a1 <- .py_getitem(.tmpu_, 1L); a2 <- .py_getitem(.tmpu_, 2L); a3 <- .py_getitem(.tmpu_, 3L); b0 <- .py_getitem(.tmpu_, 4L); b1 <- .py_getitem(.tmpu_, 5L); b2 <- .py_getitem(.tmpu_, 6L); b3 <- .py_getitem(.tmpu_, 7L); b4 <- .py_getitem(.tmpu_, 8L); c0 <- .py_getitem(.tmpu_, 9L); c1 <- .py_getitem(.tmpu_, 10L); c2 <- .py_getitem(.tmpu_, 11L); c3 <- .py_getitem(.tmpu_, 12L); d0 <- .py_getitem(.tmpu_, 13L); d1 <- .py_getitem(.tmpu_, 14L); d2 <- .py_getitem(.tmpu_, 15L); d3 <- .py_getitem(.tmpu_, 16L)
})
.py_try({
TF1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
TF2 <- TransferFunction(k, 1L, s)
})
.py_try({
TF3 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
TF4 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(tau, s), 1L), s)
})
.py_try({
TF5 <- TransferFunction(k, p, s)
})
.py_try({
TF6 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(k, s), p), s)
})
p_are_floats_equals <- function(n1, n2) {
  "\n    Helper function to compare two floats.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return(.py_method(Float(n1), "is_same", Float(n2), .py_isclose))
}
p_are_complex_equals <- function(n1, n2) {
  "\n    Helper function to compare two complex numbers.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return((.py_truthy(p_are_floats_equals(re(n1), re(n2))) && .py_truthy(p_are_floats_equals(im(n1), im(n2)))))
}
.py_try({
tf1 <- DiscreteTransferFunction(.py_add(s, 6L), .py_sub(s, 5L), s)
})
.py_try({
tf2 <- DiscreteTransferFunction(.py_add(s, 3L), .py_add(s, 1L), s)
})
.py_try({
tf3 <- DiscreteTransferFunction(.py_add(s, 1L), .py_add(.py_add(.py_pow(s, 2L), s), 1L), s)
})
.py_try({
tf4 <- DiscreteTransferFunction(p, .py_sub(2L, p), p)
})
.py_try({
tf5 <- DiscreteTransferFunction(p, .py_add(p, 5L), p, 20L)
})
.py_try({
tf6 <- DiscreteTransferFunction(1L, .py_sub(p, 1L), p, 20L)
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_addition_and_subtraction:999", function() .py_add(tf1, tf2))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_addition_and_subtraction:1000", function() .py_add(tf3, tf1))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_addition_and_subtraction:1001", function() .py_add(.py_add(.py_neg(tf1), tf2), tf3))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_addition_and_subtraction:1002", function() .py_add(tf1, .py_add(tf2, tf3)))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_addition_and_subtraction:1003", function() .py_add(tf5, tf6))
})
.py_try({
c <- .py_symbols("c", commutative = FALSE)
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_addition_and_subtraction:1006", function() .py_add(tf1, .py_Matrix(list(1L, 2L, 3L))))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_addition_and_subtraction:1007", function() .py_add(tf2, c))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_addition_and_subtraction:1008", function() .py_add(tf3, tf4))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_addition_and_subtraction:1009", function() .py_add(tf1, .py_sub(s, 1L)))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_addition_and_subtraction:1010", function() .py_add(tf1, 8L))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_addition_and_subtraction:1011", function() .py_add(.py_sub(1L, .py_pow(p, 3L)), tf1))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_addition_and_subtraction:1014", function() .py_sub(tf1, tf2))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_addition_and_subtraction:1015", function() .py_sub(tf3, tf2))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_addition_and_subtraction:1016", function() .py_sub(.py_neg(tf1), tf3))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_addition_and_subtraction:1017", function() .py_add(.py_sub(tf1, tf2), tf3))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_addition_and_subtraction:1018", function() .py_sub(tf5, tf6))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_addition_and_subtraction:1020", function() .py_sub(tf1, .py_Matrix(list(1L, 2L, 3L))))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_addition_and_subtraction:1021", function() .py_sub(tf3, tf4))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_addition_and_subtraction:1022", function() .py_sub(tf1, .py_sub(s, 1L)))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_addition_and_subtraction:1023", function() .py_sub(tf1, 8L))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_addition_and_subtraction:1024", function() .py_sub(.py_add(s, 5L), tf2))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_addition_and_subtraction:1025", function() .py_sub(.py_add(1L, .py_pow(p, 4L)), tf1))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_addition_and_subtraction:1027", function() .py_add(tf5, tf1))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_addition_and_subtraction:1028", function() .py_sub(tf5, tf1))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_addition_and_subtraction:1029", function() .py_add(.py_add(tf5, tf6), tf1))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_addition_and_subtraction:1030", function() .py_sub(.py_sub(tf5, tf6), tf1))
}))---")
  expect_identical(got[["lti:test_DiscreteTransferFunction_addition_and_subtraction:999"]], "Parallel(DiscreteTransferFunction(s + 6, s - 5, s, 1), DiscreteTransferFunction(s + 3, s + 1, s, 1))")
  expect_identical(got[["lti:test_DiscreteTransferFunction_addition_and_subtraction:1000"]], "Parallel(DiscreteTransferFunction(s + 1, s^2 + s + 1, s, 1), DiscreteTransferFunction(s + 6, s - 5, s, 1))")
  expect_identical(got[["lti:test_DiscreteTransferFunction_addition_and_subtraction:1001"]], "Parallel(DiscreteTransferFunction(-s - 6, s - 5, s, 1), DiscreteTransferFunction(s + 3, s + 1, s, 1), DiscreteTransferFunction(s + 1, s^2 + s + 1, s, 1))")
  expect_identical(got[["lti:test_DiscreteTransferFunction_addition_and_subtraction:1002"]], "Parallel(DiscreteTransferFunction(s + 6, s - 5, s, 1), DiscreteTransferFunction(s + 3, s + 1, s, 1), DiscreteTransferFunction(s + 1, s^2 + s + 1, s, 1))")
  expect_identical(got[["lti:test_DiscreteTransferFunction_addition_and_subtraction:1003"]], "Parallel(DiscreteTransferFunction(p, p + 5, p, 20), DiscreteTransferFunction(1, p - 1, p, 20))")
  expect_identical(got[["lti:test_DiscreteTransferFunction_addition_and_subtraction:1006"]], "ERR ValueError: <class 'sympy.physics.control.lti.DiscreteTransferFunction'> cannot be added with <class 'sympy.matrices.dense.MutableDenseMatrix'>.")
  expect_identical(got[["lti:test_DiscreteTransferFunction_addition_and_subtraction:1007"]], "ERR ValueError: <class 'sympy.physics.control.lti.DiscreteTransferFunction'> cannot be added with <class 'sympy.core.symbol.Symbol'>.")
  expect_identical(got[["lti:test_DiscreteTransferFunction_addition_and_subtraction:1008"]], "ERR ValueError: \\nAll the transfer functions should use the same complex variable of the\\nLaplace domain or z-domain.")
  expect_identical(got[["lti:test_DiscreteTransferFunction_addition_and_subtraction:1009"]], "ERR ValueError: <class 'sympy.physics.control.lti.DiscreteTransferFunction'> cannot be added with <class 'sympy.core.add.Add'>.")
  expect_identical(got[["lti:test_DiscreteTransferFunction_addition_and_subtraction:1010"]], "ERR ValueError: <class 'sympy.physics.control.lti.DiscreteTransferFunction'> cannot be added with <class 'int'>.")
  expect_identical(got[["lti:test_DiscreteTransferFunction_addition_and_subtraction:1011"]], "ERR ValueError: <class 'sympy.physics.control.lti.DiscreteTransferFunction'> cannot be added with <class 'sympy.core.add.Add'>.")
  expect_identical(got[["lti:test_DiscreteTransferFunction_addition_and_subtraction:1014"]], "Parallel(DiscreteTransferFunction(s + 6, s - 5, s, 1), DiscreteTransferFunction(-s - 3, s + 1, s, 1))")
  expect_identical(got[["lti:test_DiscreteTransferFunction_addition_and_subtraction:1015"]], "Parallel(DiscreteTransferFunction(s + 1, s^2 + s + 1, s, 1), DiscreteTransferFunction(-s - 3, s + 1, s, 1))")
  expect_identical(got[["lti:test_DiscreteTransferFunction_addition_and_subtraction:1016"]], "Parallel(DiscreteTransferFunction(-s - 6, s - 5, s, 1), DiscreteTransferFunction(-s - 1, s^2 + s + 1, s, 1))")
  expect_identical(got[["lti:test_DiscreteTransferFunction_addition_and_subtraction:1017"]], "Parallel(DiscreteTransferFunction(s + 6, s - 5, s, 1), DiscreteTransferFunction(-s - 3, s + 1, s, 1), DiscreteTransferFunction(s + 1, s^2 + s + 1, s, 1))")
  expect_identical(got[["lti:test_DiscreteTransferFunction_addition_and_subtraction:1018"]], "Parallel(DiscreteTransferFunction(p, p + 5, p, 20), DiscreteTransferFunction(-1, p - 1, p, 20))")
  expect_identical(got[["lti:test_DiscreteTransferFunction_addition_and_subtraction:1020"]], "ERR ValueError: <class 'sympy.physics.control.lti.DiscreteTransferFunction'> cannot be subtracted from a <class 'sympy.matrices.dense.MutableDenseMatrix'>.")
  expect_identical(got[["lti:test_DiscreteTransferFunction_addition_and_subtraction:1021"]], "ERR ValueError: \\nAll the transfer functions should use the same complex variable of the\\nLaplace domain or z-domain.")
  expect_identical(got[["lti:test_DiscreteTransferFunction_addition_and_subtraction:1022"]], "ERR ValueError: <class 'sympy.physics.control.lti.DiscreteTransferFunction'> cannot be subtracted from a <class 'sympy.core.add.Add'>.")
  expect_identical(got[["lti:test_DiscreteTransferFunction_addition_and_subtraction:1023"]], "ERR ValueError: <class 'sympy.physics.control.lti.DiscreteTransferFunction'> cannot be subtracted from a <class 'int'>.")
  expect_identical(got[["lti:test_DiscreteTransferFunction_addition_and_subtraction:1024"]], "ERR ValueError: <class 'sympy.physics.control.lti.DiscreteTransferFunction'> cannot be added with <class 'sympy.core.add.Add'>.")
  expect_identical(got[["lti:test_DiscreteTransferFunction_addition_and_subtraction:1025"]], "ERR ValueError: <class 'sympy.physics.control.lti.DiscreteTransferFunction'> cannot be added with <class 'sympy.core.add.Add'>.")
  expect_identical(got[["lti:test_DiscreteTransferFunction_addition_and_subtraction:1027"]], "ERR TypeError: \\n                Incompatible sampling times detected.\\n                All systems should have the same sampling time.\\n                Found sampling time: 20 and 1")
  expect_identical(got[["lti:test_DiscreteTransferFunction_addition_and_subtraction:1028"]], "ERR TypeError: \\n                Incompatible sampling times detected.\\n                All systems should have the same sampling time.\\n                Found sampling time: 20 and 1")
  expect_identical(got[["lti:test_DiscreteTransferFunction_addition_and_subtraction:1029"]], "ERR TypeError: \\n                Incompatible sampling times detected.\\n                All systems should have the same sampling time.\\n                Found sampling time: 20 and 1")
  expect_identical(got[["lti:test_DiscreteTransferFunction_addition_and_subtraction:1030"]], "ERR TypeError: \\n                Incompatible sampling times detected.\\n                All systems should have the same sampling time.\\n                Found sampling time: 20 and 1")
})

test_that("lti test_DiscreteTransferFunction_multiplication_and_division matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("a, x, b, c, s, g, d,    p, k, tau, zeta, wn, T, z"); a <- .py_getitem(.tmpu_, 0L); x <- .py_getitem(.tmpu_, 1L); b <- .py_getitem(.tmpu_, 2L); c <- .py_getitem(.tmpu_, 3L); s <- .py_getitem(.tmpu_, 4L); g <- .py_getitem(.tmpu_, 5L); d <- .py_getitem(.tmpu_, 6L); p <- .py_getitem(.tmpu_, 7L); k <- .py_getitem(.tmpu_, 8L); tau <- .py_getitem(.tmpu_, 9L); zeta <- .py_getitem(.tmpu_, 10L); wn <- .py_getitem(.tmpu_, 11L); T <- .py_getitem(.tmpu_, 12L); z <- .py_getitem(.tmpu_, 13L)
})
.py_try({
.tmpu_ <- .py_symbols("a0:4, b0:5, c0:4, d0:4"); a0 <- .py_getitem(.tmpu_, 0L); a1 <- .py_getitem(.tmpu_, 1L); a2 <- .py_getitem(.tmpu_, 2L); a3 <- .py_getitem(.tmpu_, 3L); b0 <- .py_getitem(.tmpu_, 4L); b1 <- .py_getitem(.tmpu_, 5L); b2 <- .py_getitem(.tmpu_, 6L); b3 <- .py_getitem(.tmpu_, 7L); b4 <- .py_getitem(.tmpu_, 8L); c0 <- .py_getitem(.tmpu_, 9L); c1 <- .py_getitem(.tmpu_, 10L); c2 <- .py_getitem(.tmpu_, 11L); c3 <- .py_getitem(.tmpu_, 12L); d0 <- .py_getitem(.tmpu_, 13L); d1 <- .py_getitem(.tmpu_, 14L); d2 <- .py_getitem(.tmpu_, 15L); d3 <- .py_getitem(.tmpu_, 16L)
})
.py_try({
TF1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
TF2 <- TransferFunction(k, 1L, s)
})
.py_try({
TF3 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
TF4 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(tau, s), 1L), s)
})
.py_try({
TF5 <- TransferFunction(k, p, s)
})
.py_try({
TF6 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(k, s), p), s)
})
p_are_floats_equals <- function(n1, n2) {
  "\n    Helper function to compare two floats.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return(.py_method(Float(n1), "is_same", Float(n2), .py_isclose))
}
p_are_complex_equals <- function(n1, n2) {
  "\n    Helper function to compare two complex numbers.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return((.py_truthy(p_are_floats_equals(re(n1), re(n2))) && .py_truthy(p_are_floats_equals(im(n1), im(n2)))))
}
.py_try({
G1 <- DiscreteTransferFunction(.py_add(s, 3L), .py_add(.py_neg(.py_pow(s, 3L)), 9L), s)
})
.py_try({
G2 <- DiscreteTransferFunction(.py_add(s, 1L), .py_sub(s, 5L), s)
})
.py_try({
G3 <- DiscreteTransferFunction(p, .py_sub(.py_pow(p, 4L), 6L), p)
})
.py_try({
G4 <- DiscreteTransferFunction(.py_add(p, 4L), .py_sub(p, 5L), p)
})
.py_try({
G5 <- DiscreteTransferFunction(.py_add(s, 6L), .py_sub(s, 5L), s)
})
.py_try({
G6 <- DiscreteTransferFunction(.py_add(s, 3L), .py_add(s, 1L), s)
})
.py_try({
G7 <- DiscreteTransferFunction(1L, 1L, s)
})
.py_try({
G8 <- DiscreteTransferFunction(p, .py_add(p, 1L), p, 12L)
})
.py_try({
G9 <- DiscreteTransferFunction(1L, p, p, 12L)
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_multiplication_and_division:1044", function() .py_mul(G1, G2))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_multiplication_and_division:1045", function() .py_mul(.py_neg(G1), G5))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_multiplication_and_division:1046", function() .py_mul(.py_mul(.py_neg(G2), G5), .py_neg(G6)))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_multiplication_and_division:1047", function() .py_mul(.py_mul(.py_mul(.py_neg(G1), .py_neg(G2)), .py_neg(G5)), .py_neg(G6)))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_multiplication_and_division:1048", function() .py_mul(G3, G4))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_multiplication_and_division:1049", function() .py_mul(.py_mul(G1, G2), .py_neg(.py_mul(G5, G6))))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_multiplication_and_division:1051", function() .py_mul(.py_mul(G1, G2), .py_add(G5, G6)))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_multiplication_and_division:1052", function() .py_mul(G8, G9))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_multiplication_and_division:1055", function() .py_div(G5, G6))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_multiplication_and_division:1056", function() .py_div(.py_neg(G3), G4))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_multiplication_and_division:1057", function() .py_div(.py_mul(G5, G6), G7))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_multiplication_and_division:1058", function() .py_div(G8, G9))
})
.py_try({
c <- .py_symbols("c", commutative = FALSE)
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_multiplication_and_division:1061", function() .py_mul(G3, .py_Matrix(list(1L, 2L, 3L))))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_multiplication_and_division:1062", function() .py_mul(G1, c))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_multiplication_and_division:1063", function() .py_mul(G3, G5))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_multiplication_and_division:1064", function() .py_mul(G5, .py_sub(s, 1L)))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_multiplication_and_division:1065", function() .py_mul(9L, G5))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_multiplication_and_division:1067", function() .py_div(G3, .py_Matrix(list(1L, 2L, 3L))))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_multiplication_and_division:1068", function() .py_div(G6, 0L))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_multiplication_and_division:1069", function() .py_div(G3, G5))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_multiplication_and_division:1070", function() .py_div(G5, 2L))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_multiplication_and_division:1071", function() .py_div(G5, .py_pow(s, 2L)))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_multiplication_and_division:1072", function() .py_div(.py_sub(s, .py_mul(4L, .py_pow(s, 2L))), G2))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_multiplication_and_division:1073", function() .py_div(0L, G4))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_multiplication_and_division:1074", function() .py_div(G7, .py_add(1L, G6)))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_multiplication_and_division:1075", function() .py_div(G7, .py_mul(G5, G6)))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_multiplication_and_division:1076", function() .py_div(G7, .py_add(G7, .py_add(G5, G6))))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_multiplication_and_division:1079", function() .py_mul(G8, G2))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_multiplication_and_division:1080", function() .py_mul(.py_mul(G8, G9), G2))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_multiplication_and_division:1081", function() .py_mul(.py_div(G8, G9), G2))
}))---")
  expect_identical(got[["lti:test_DiscreteTransferFunction_multiplication_and_division:1044"]], "Series(DiscreteTransferFunction(s + 3, 9 - s^3, s, 1), DiscreteTransferFunction(s + 1, s - 5, s, 1))")
  expect_identical(got[["lti:test_DiscreteTransferFunction_multiplication_and_division:1045"]], "Series(DiscreteTransferFunction(-s - 3, 9 - s^3, s, 1), DiscreteTransferFunction(s + 6, s - 5, s, 1))")
  expect_identical(got[["lti:test_DiscreteTransferFunction_multiplication_and_division:1046"]], "Series(DiscreteTransferFunction(-s - 1, s - 5, s, 1), DiscreteTransferFunction(s + 6, s - 5, s, 1), DiscreteTransferFunction(-s - 3, s + 1, s, 1))")
  expect_identical(got[["lti:test_DiscreteTransferFunction_multiplication_and_division:1047"]], "Series(DiscreteTransferFunction(-s - 3, 9 - s^3, s, 1), DiscreteTransferFunction(-s - 1, s - 5, s, 1), DiscreteTransferFunction(-s - 6, s - 5, s, 1), DiscreteTransferFunction(-s - 3, s + 1, s, 1))")
  expect_identical(got[["lti:test_DiscreteTransferFunction_multiplication_and_division:1048"]], "Series(DiscreteTransferFunction(p, p^4 - 6, p, 1), DiscreteTransferFunction(p + 4, p - 5, p, 1))")
  expect_identical(got[["lti:test_DiscreteTransferFunction_multiplication_and_division:1049"]], "Series(DiscreteTransferFunction(s + 3, 9 - s^3, s, 1), DiscreteTransferFunction(s + 1, s - 5, s, 1), DiscreteTransferFunction(-1, 1, s, 1), DiscreteTransferFunction(s + 6, s - 5, s, 1), DiscreteTransferFunction(s + 3, s + 1, s, 1))")
  expect_identical(got[["lti:test_DiscreteTransferFunction_multiplication_and_division:1051"]], "Series(DiscreteTransferFunction(s + 3, 9 - s^3, s, 1), DiscreteTransferFunction(s + 1, s - 5, s, 1), Parallel(DiscreteTransferFunction(s + 6, s - 5, s, 1), DiscreteTransferFunction(s + 3, s + 1, s, 1)))")
  expect_identical(got[["lti:test_DiscreteTransferFunction_multiplication_and_division:1052"]], "Series(DiscreteTransferFunction(p, p + 1, p, 12), DiscreteTransferFunction(1, p, p, 12))")
  expect_identical(got[["lti:test_DiscreteTransferFunction_multiplication_and_division:1055"]], "Series(DiscreteTransferFunction(s + 6, s - 5, s, 1), DiscreteTransferFunction(s + 1, s + 3, s, 1))")
  expect_identical(got[["lti:test_DiscreteTransferFunction_multiplication_and_division:1056"]], "Series(DiscreteTransferFunction(-p, p^4 - 6, p, 1), DiscreteTransferFunction(p - 5, p + 4, p, 1))")
  expect_identical(got[["lti:test_DiscreteTransferFunction_multiplication_and_division:1057"]], "Series(DiscreteTransferFunction(s + 6, s - 5, s, 1), DiscreteTransferFunction(s + 3, s + 1, s, 1), DiscreteTransferFunction(1, 1, s, 1))")
  expect_identical(got[["lti:test_DiscreteTransferFunction_multiplication_and_division:1058"]], "Series(DiscreteTransferFunction(p, p + 1, p, 12), DiscreteTransferFunction(p, 1, p, 12))")
  expect_identical(got[["lti:test_DiscreteTransferFunction_multiplication_and_division:1061"]], "ERR ValueError: <class 'sympy.physics.control.lti.DiscreteTransferFunction'> cannot be multiplied with <class 'sympy.matrices.dense.MutableDenseMatrix'>.")
  expect_identical(got[["lti:test_DiscreteTransferFunction_multiplication_and_division:1062"]], "ERR ValueError: <class 'sympy.physics.control.lti.DiscreteTransferFunction'> cannot be multiplied with <class 'sympy.core.symbol.Symbol'>.")
  expect_identical(got[["lti:test_DiscreteTransferFunction_multiplication_and_division:1063"]], "ERR ValueError: \\nAll the transfer functions should use the same complex variable of the\\nLaplace domain or z-domain.")
  expect_identical(got[["lti:test_DiscreteTransferFunction_multiplication_and_division:1064"]], "ERR ValueError: <class 'sympy.physics.control.lti.DiscreteTransferFunction'> cannot be multiplied with <class 'sympy.core.add.Add'>.")
  expect_identical(got[["lti:test_DiscreteTransferFunction_multiplication_and_division:1065"]], "ERR ValueError: <class 'sympy.physics.control.lti.DiscreteTransferFunction'> cannot be multiplied with <class 'int'>.")
  expect_identical(got[["lti:test_DiscreteTransferFunction_multiplication_and_division:1067"]], "ERR ValueError: <class 'sympy.physics.control.lti.DiscreteTransferFunction'> cannot be divided by <class 'sympy.matrices.dense.MutableDenseMatrix'>.")
  expect_identical(got[["lti:test_DiscreteTransferFunction_multiplication_and_division:1068"]], "ERR ValueError: <class 'sympy.physics.control.lti.DiscreteTransferFunction'> cannot be divided by <class 'int'>.")
  expect_identical(got[["lti:test_DiscreteTransferFunction_multiplication_and_division:1069"]], "ERR ValueError: \\nAll the transfer functions should use the same complex variable of the\\nLaplace domain or z-domain.")
  expect_identical(got[["lti:test_DiscreteTransferFunction_multiplication_and_division:1070"]], "ERR ValueError: <class 'sympy.physics.control.lti.DiscreteTransferFunction'> cannot be divided by <class 'int'>.")
  expect_identical(got[["lti:test_DiscreteTransferFunction_multiplication_and_division:1071"]], "ERR ValueError: <class 'sympy.physics.control.lti.DiscreteTransferFunction'> cannot be divided by <class 'sympy.core.power.Pow'>.")
  expect_identical(got[["lti:test_DiscreteTransferFunction_multiplication_and_division:1072"]], "ERR ValueError: <class 'sympy.physics.control.lti.DiscreteTransferFunction'> cannot be divided by <class 'sympy.core.add.Add'>.")
  expect_identical(got[["lti:test_DiscreteTransferFunction_multiplication_and_division:1073"]], "ERR ValueError: <class 'sympy.physics.control.lti.DiscreteTransferFunction'> cannot be divided by <class 'int'>.")
  expect_identical(got[["lti:test_DiscreteTransferFunction_multiplication_and_division:1074"]], "ERR ValueError: <class 'sympy.physics.control.lti.DiscreteTransferFunction'> cannot be added with <class 'int'>.")
  expect_identical(got[["lti:test_DiscreteTransferFunction_multiplication_and_division:1075"]], "ERR ValueError: <class 'sympy.physics.control.lti.DiscreteTransferFunction'> cannot be divided by <class 'sympy.physics.control.lti.Series'>.")
  expect_identical(got[["lti:test_DiscreteTransferFunction_multiplication_and_division:1076"]], "ERR ValueError: <class 'sympy.physics.control.lti.DiscreteTransferFunction'> cannot be divided by <class 'sympy.physics.control.lti.Parallel'>.")
  expect_identical(got[["lti:test_DiscreteTransferFunction_multiplication_and_division:1079"]], "ERR TypeError: \\n                Incompatible sampling times detected.\\n                All systems should have the same sampling time.\\n                Found sampling time: 12 and 1")
  expect_identical(got[["lti:test_DiscreteTransferFunction_multiplication_and_division:1080"]], "ERR TypeError: \\n                Incompatible sampling times detected.\\n                All systems should have the same sampling time.\\n                Found sampling time: 12 and 1")
  expect_identical(got[["lti:test_DiscreteTransferFunction_multiplication_and_division:1081"]], "ERR TypeError: \\n                Incompatible sampling times detected.\\n                All systems should have the same sampling time.\\n                Found sampling time: 12 and 1")
})

test_that("lti test_DiscreteTransferFunction_is_proper matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("a, x, b, c, s, g, d,    p, k, tau, zeta, wn, T, z"); a <- .py_getitem(.tmpu_, 0L); x <- .py_getitem(.tmpu_, 1L); b <- .py_getitem(.tmpu_, 2L); c <- .py_getitem(.tmpu_, 3L); s <- .py_getitem(.tmpu_, 4L); g <- .py_getitem(.tmpu_, 5L); d <- .py_getitem(.tmpu_, 6L); p <- .py_getitem(.tmpu_, 7L); k <- .py_getitem(.tmpu_, 8L); tau <- .py_getitem(.tmpu_, 9L); zeta <- .py_getitem(.tmpu_, 10L); wn <- .py_getitem(.tmpu_, 11L); T <- .py_getitem(.tmpu_, 12L); z <- .py_getitem(.tmpu_, 13L)
})
.py_try({
.tmpu_ <- .py_symbols("a0:4, b0:5, c0:4, d0:4"); a0 <- .py_getitem(.tmpu_, 0L); a1 <- .py_getitem(.tmpu_, 1L); a2 <- .py_getitem(.tmpu_, 2L); a3 <- .py_getitem(.tmpu_, 3L); b0 <- .py_getitem(.tmpu_, 4L); b1 <- .py_getitem(.tmpu_, 5L); b2 <- .py_getitem(.tmpu_, 6L); b3 <- .py_getitem(.tmpu_, 7L); b4 <- .py_getitem(.tmpu_, 8L); c0 <- .py_getitem(.tmpu_, 9L); c1 <- .py_getitem(.tmpu_, 10L); c2 <- .py_getitem(.tmpu_, 11L); c3 <- .py_getitem(.tmpu_, 12L); d0 <- .py_getitem(.tmpu_, 13L); d1 <- .py_getitem(.tmpu_, 14L); d2 <- .py_getitem(.tmpu_, 15L); d3 <- .py_getitem(.tmpu_, 16L)
})
.py_try({
TF1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
TF2 <- TransferFunction(k, 1L, s)
})
.py_try({
TF3 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
TF4 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(tau, s), 1L), s)
})
.py_try({
TF5 <- TransferFunction(k, p, s)
})
.py_try({
TF6 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(k, s), p), s)
})
p_are_floats_equals <- function(n1, n2) {
  "\n    Helper function to compare two floats.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return(.py_method(Float(n1), "is_same", Float(n2), .py_isclose))
}
p_are_complex_equals <- function(n1, n2) {
  "\n    Helper function to compare two complex numbers.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return((.py_truthy(p_are_floats_equals(re(n1), re(n2))) && .py_truthy(p_are_floats_equals(im(n1), im(n2)))))
}
.py_try({
.tmpu_ <- .py_symbols("omega_o, zeta, tau"); omega_o <- .py_getitem(.tmpu_, 0L); zeta <- .py_getitem(.tmpu_, 1L); tau <- .py_getitem(.tmpu_, 2L)
})
.py_try({
G1 <- DiscreteTransferFunction(.py_pow(omega_o, 2L), .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(p, omega_o), zeta), s)), .py_pow(omega_o, 2L)), omega_o)
})
.py_try({
G2 <- DiscreteTransferFunction(.py_sub(tau, .py_pow(s, 3L)), .py_add(tau, .py_pow(p, 4L)), tau)
})
.py_try({
G3 <- DiscreteTransferFunction(.py_add(.py_sub(.py_add(.py_mul(.py_mul(a, b), .py_pow(s, 3L)), .py_pow(s, 2L)), .py_mul(a, p)), s), .py_sub(b, .py_mul(s, .py_pow(p, 2L))), p)
})
.py_try({
G4 <- DiscreteTransferFunction(.py_add(.py_sub(.py_add(.py_mul(b, .py_pow(s, 2L)), .py_pow(p, 2L)), .py_mul(a, p)), s), .py_sub(b, .py_pow(p, 2L)), s)
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_is_proper:1091", function() .py_attr(G1, "is_proper"))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_is_proper:1092", function() .py_attr(G2, "is_proper"))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_is_proper:1093", function() .py_attr(G3, "is_proper"))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_is_proper:1094", function() .py_attr(G4, "is_proper"))
}))---")
  expect_identical(got[["lti:test_DiscreteTransferFunction_is_proper:1091"]], "True")
  expect_identical(got[["lti:test_DiscreteTransferFunction_is_proper:1092"]], "True")
  expect_identical(got[["lti:test_DiscreteTransferFunction_is_proper:1093"]], "True")
  expect_identical(got[["lti:test_DiscreteTransferFunction_is_proper:1094"]], "False")
})

test_that("lti test_DiscreteTransferFunction_is_strictly_proper matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("a, x, b, c, s, g, d,    p, k, tau, zeta, wn, T, z"); a <- .py_getitem(.tmpu_, 0L); x <- .py_getitem(.tmpu_, 1L); b <- .py_getitem(.tmpu_, 2L); c <- .py_getitem(.tmpu_, 3L); s <- .py_getitem(.tmpu_, 4L); g <- .py_getitem(.tmpu_, 5L); d <- .py_getitem(.tmpu_, 6L); p <- .py_getitem(.tmpu_, 7L); k <- .py_getitem(.tmpu_, 8L); tau <- .py_getitem(.tmpu_, 9L); zeta <- .py_getitem(.tmpu_, 10L); wn <- .py_getitem(.tmpu_, 11L); T <- .py_getitem(.tmpu_, 12L); z <- .py_getitem(.tmpu_, 13L)
})
.py_try({
.tmpu_ <- .py_symbols("a0:4, b0:5, c0:4, d0:4"); a0 <- .py_getitem(.tmpu_, 0L); a1 <- .py_getitem(.tmpu_, 1L); a2 <- .py_getitem(.tmpu_, 2L); a3 <- .py_getitem(.tmpu_, 3L); b0 <- .py_getitem(.tmpu_, 4L); b1 <- .py_getitem(.tmpu_, 5L); b2 <- .py_getitem(.tmpu_, 6L); b3 <- .py_getitem(.tmpu_, 7L); b4 <- .py_getitem(.tmpu_, 8L); c0 <- .py_getitem(.tmpu_, 9L); c1 <- .py_getitem(.tmpu_, 10L); c2 <- .py_getitem(.tmpu_, 11L); c3 <- .py_getitem(.tmpu_, 12L); d0 <- .py_getitem(.tmpu_, 13L); d1 <- .py_getitem(.tmpu_, 14L); d2 <- .py_getitem(.tmpu_, 15L); d3 <- .py_getitem(.tmpu_, 16L)
})
.py_try({
TF1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
TF2 <- TransferFunction(k, 1L, s)
})
.py_try({
TF3 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
TF4 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(tau, s), 1L), s)
})
.py_try({
TF5 <- TransferFunction(k, p, s)
})
.py_try({
TF6 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(k, s), p), s)
})
p_are_floats_equals <- function(n1, n2) {
  "\n    Helper function to compare two floats.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return(.py_method(Float(n1), "is_same", Float(n2), .py_isclose))
}
p_are_complex_equals <- function(n1, n2) {
  "\n    Helper function to compare two complex numbers.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return((.py_truthy(p_are_floats_equals(re(n1), re(n2))) && .py_truthy(p_are_floats_equals(im(n1), im(n2)))))
}
.py_try({
.tmpu_ <- .py_symbols("omega_o, zeta, tau"); omega_o <- .py_getitem(.tmpu_, 0L); zeta <- .py_getitem(.tmpu_, 1L); tau <- .py_getitem(.tmpu_, 2L)
})
.py_try({
tf1 <- DiscreteTransferFunction(.py_pow(omega_o, 2L), .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(p, omega_o), zeta), s)), .py_pow(omega_o, 2L)), omega_o)
})
.py_try({
tf2 <- DiscreteTransferFunction(.py_sub(tau, .py_pow(s, 3L)), .py_add(tau, .py_pow(p, 4L)), tau)
})
.py_try({
tf3 <- DiscreteTransferFunction(.py_add(.py_sub(.py_add(.py_mul(.py_mul(a, b), .py_pow(s, 3L)), .py_pow(s, 2L)), .py_mul(a, p)), s), .py_sub(b, .py_mul(s, .py_pow(p, 2L))), p)
})
.py_try({
tf4 <- DiscreteTransferFunction(.py_add(.py_sub(.py_add(.py_mul(b, .py_pow(s, 2L)), .py_pow(p, 2L)), .py_mul(a, p)), s), .py_sub(b, .py_pow(p, 2L)), s)
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_is_strictly_proper:1104", function() .py_attr(tf1, "is_strictly_proper"))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_is_strictly_proper:1105", function() .py_attr(tf2, "is_strictly_proper"))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_is_strictly_proper:1106", function() .py_attr(tf3, "is_strictly_proper"))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_is_strictly_proper:1107", function() .py_attr(tf4, "is_strictly_proper"))
}))---")
  expect_identical(got[["lti:test_DiscreteTransferFunction_is_strictly_proper:1104"]], "False")
  expect_identical(got[["lti:test_DiscreteTransferFunction_is_strictly_proper:1105"]], "False")
  expect_identical(got[["lti:test_DiscreteTransferFunction_is_strictly_proper:1106"]], "True")
  expect_identical(got[["lti:test_DiscreteTransferFunction_is_strictly_proper:1107"]], "False")
})

test_that("lti test_DiscreteTransferFunction_is_biproper matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("a, x, b, c, s, g, d,    p, k, tau, zeta, wn, T, z"); a <- .py_getitem(.tmpu_, 0L); x <- .py_getitem(.tmpu_, 1L); b <- .py_getitem(.tmpu_, 2L); c <- .py_getitem(.tmpu_, 3L); s <- .py_getitem(.tmpu_, 4L); g <- .py_getitem(.tmpu_, 5L); d <- .py_getitem(.tmpu_, 6L); p <- .py_getitem(.tmpu_, 7L); k <- .py_getitem(.tmpu_, 8L); tau <- .py_getitem(.tmpu_, 9L); zeta <- .py_getitem(.tmpu_, 10L); wn <- .py_getitem(.tmpu_, 11L); T <- .py_getitem(.tmpu_, 12L); z <- .py_getitem(.tmpu_, 13L)
})
.py_try({
.tmpu_ <- .py_symbols("a0:4, b0:5, c0:4, d0:4"); a0 <- .py_getitem(.tmpu_, 0L); a1 <- .py_getitem(.tmpu_, 1L); a2 <- .py_getitem(.tmpu_, 2L); a3 <- .py_getitem(.tmpu_, 3L); b0 <- .py_getitem(.tmpu_, 4L); b1 <- .py_getitem(.tmpu_, 5L); b2 <- .py_getitem(.tmpu_, 6L); b3 <- .py_getitem(.tmpu_, 7L); b4 <- .py_getitem(.tmpu_, 8L); c0 <- .py_getitem(.tmpu_, 9L); c1 <- .py_getitem(.tmpu_, 10L); c2 <- .py_getitem(.tmpu_, 11L); c3 <- .py_getitem(.tmpu_, 12L); d0 <- .py_getitem(.tmpu_, 13L); d1 <- .py_getitem(.tmpu_, 14L); d2 <- .py_getitem(.tmpu_, 15L); d3 <- .py_getitem(.tmpu_, 16L)
})
.py_try({
TF1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
TF2 <- TransferFunction(k, 1L, s)
})
.py_try({
TF3 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
TF4 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(tau, s), 1L), s)
})
.py_try({
TF5 <- TransferFunction(k, p, s)
})
.py_try({
TF6 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(k, s), p), s)
})
p_are_floats_equals <- function(n1, n2) {
  "\n    Helper function to compare two floats.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return(.py_method(Float(n1), "is_same", Float(n2), .py_isclose))
}
p_are_complex_equals <- function(n1, n2) {
  "\n    Helper function to compare two complex numbers.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return((.py_truthy(p_are_floats_equals(re(n1), re(n2))) && .py_truthy(p_are_floats_equals(im(n1), im(n2)))))
}
.py_try({
.tmpu_ <- .py_symbols("tau, omega_o, zeta"); tau <- .py_getitem(.tmpu_, 0L); omega_o <- .py_getitem(.tmpu_, 1L); zeta <- .py_getitem(.tmpu_, 2L)
})
.py_try({
tf1 <- DiscreteTransferFunction(.py_pow(omega_o, 2L), .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(p, omega_o), zeta), s)), .py_pow(omega_o, 2L)), omega_o)
})
.py_try({
tf2 <- DiscreteTransferFunction(.py_sub(tau, .py_pow(s, 3L)), .py_add(tau, .py_pow(p, 4L)), tau)
})
.py_try({
tf3 <- DiscreteTransferFunction(.py_add(.py_sub(.py_add(.py_mul(.py_mul(a, b), .py_pow(s, 3L)), .py_pow(s, 2L)), .py_mul(a, p)), s), .py_sub(b, .py_mul(s, .py_pow(p, 2L))), p)
})
.py_try({
tf4 <- DiscreteTransferFunction(.py_add(.py_sub(.py_add(.py_mul(b, .py_pow(s, 2L)), .py_pow(p, 2L)), .py_mul(a, p)), s), .py_sub(b, .py_pow(p, 2L)), s)
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_is_biproper:1116", function() .py_attr(tf1, "is_biproper"))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_is_biproper:1117", function() .py_attr(tf2, "is_biproper"))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_is_biproper:1118", function() .py_attr(tf3, "is_biproper"))
})
.py_try({
.py_rec("lti:test_DiscreteTransferFunction_is_biproper:1119", function() .py_attr(tf4, "is_biproper"))
}))---")
  expect_identical(got[["lti:test_DiscreteTransferFunction_is_biproper:1116"]], "True")
  expect_identical(got[["lti:test_DiscreteTransferFunction_is_biproper:1117"]], "True")
  expect_identical(got[["lti:test_DiscreteTransferFunction_is_biproper:1118"]], "False")
  expect_identical(got[["lti:test_DiscreteTransferFunction_is_biproper:1119"]], "False")
})

test_that("lti test_PIDController matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("a, x, b, c, s, g, d,    p, k, tau, zeta, wn, T, z"); a <- .py_getitem(.tmpu_, 0L); x <- .py_getitem(.tmpu_, 1L); b <- .py_getitem(.tmpu_, 2L); c <- .py_getitem(.tmpu_, 3L); s <- .py_getitem(.tmpu_, 4L); g <- .py_getitem(.tmpu_, 5L); d <- .py_getitem(.tmpu_, 6L); p <- .py_getitem(.tmpu_, 7L); k <- .py_getitem(.tmpu_, 8L); tau <- .py_getitem(.tmpu_, 9L); zeta <- .py_getitem(.tmpu_, 10L); wn <- .py_getitem(.tmpu_, 11L); T <- .py_getitem(.tmpu_, 12L); z <- .py_getitem(.tmpu_, 13L)
})
.py_try({
.tmpu_ <- .py_symbols("a0:4, b0:5, c0:4, d0:4"); a0 <- .py_getitem(.tmpu_, 0L); a1 <- .py_getitem(.tmpu_, 1L); a2 <- .py_getitem(.tmpu_, 2L); a3 <- .py_getitem(.tmpu_, 3L); b0 <- .py_getitem(.tmpu_, 4L); b1 <- .py_getitem(.tmpu_, 5L); b2 <- .py_getitem(.tmpu_, 6L); b3 <- .py_getitem(.tmpu_, 7L); b4 <- .py_getitem(.tmpu_, 8L); c0 <- .py_getitem(.tmpu_, 9L); c1 <- .py_getitem(.tmpu_, 10L); c2 <- .py_getitem(.tmpu_, 11L); c3 <- .py_getitem(.tmpu_, 12L); d0 <- .py_getitem(.tmpu_, 13L); d1 <- .py_getitem(.tmpu_, 14L); d2 <- .py_getitem(.tmpu_, 15L); d3 <- .py_getitem(.tmpu_, 16L)
})
.py_try({
TF1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
TF2 <- TransferFunction(k, 1L, s)
})
.py_try({
TF3 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
TF4 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(tau, s), 1L), s)
})
.py_try({
TF5 <- TransferFunction(k, p, s)
})
.py_try({
TF6 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(k, s), p), s)
})
p_are_floats_equals <- function(n1, n2) {
  "\n    Helper function to compare two floats.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return(.py_method(Float(n1), "is_same", Float(n2), .py_isclose))
}
p_are_complex_equals <- function(n1, n2) {
  "\n    Helper function to compare two complex numbers.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return((.py_truthy(p_are_floats_equals(re(n1), re(n2))) && .py_truthy(p_are_floats_equals(im(n1), im(n2)))))
}
.py_try({
.tmpu_ <- .py_symbols("kp ki kd tf"); kp <- .py_getitem(.tmpu_, 0L); ki <- .py_getitem(.tmpu_, 1L); kd <- .py_getitem(.tmpu_, 2L); tf <- .py_getitem(.tmpu_, 3L)
})
.py_try({
p1 <- PIDController(kp, ki, kd, tf)
})
.py_try({
p2 <- PIDController()
})
.py_try({
.py_rec("lti:test_PIDController:1128", function() .py_isinstance(p1, PIDController))
})
.py_try({
.py_rec("lti:test_PIDController:1129", function() .py_isinstance(p1, TransferFunction))
})
.py_try({
.py_rec("lti:test_PIDController:1132", function() p1)
})
.py_try({
.py_rec("lti:test_PIDController:1133", function() p2)
})
.py_try({
.py_rec("lti:test_PIDController:1134", function() .py_attr(p1, "num"))
})
.py_try({
.py_rec("lti:test_PIDController:1135", function() .py_attr(p1, "den"))
})
.py_try({
.py_rec("lti:test_PIDController:1136", function() .py_attr(p1, "var"))
})
.py_try({
.py_rec("lti:test_PIDController:1137", function() .py_attr(p1, "kp"))
})
.py_try({
.py_rec("lti:test_PIDController:1138", function() .py_attr(p1, "ki"))
})
.py_try({
.py_rec("lti:test_PIDController:1139", function() .py_attr(p1, "kd"))
})
.py_try({
.py_rec("lti:test_PIDController:1140", function() .py_attr(p1, "tf"))
})
.py_try({
.py_rec("lti:test_PIDController:1143", function() .py_method(p1, "doit"))
})
.py_try({
.py_rec("lti:test_PIDController:1144", function() .py_attr(p1, "is_proper"))
})
.py_try({
.py_rec("lti:test_PIDController:1145", function() .py_attr(p1, "is_biproper"))
})
.py_try({
.py_rec("lti:test_PIDController:1146", function() .py_attr(p1, "is_strictly_proper"))
})
.py_try({
.py_rec("lti:test_PIDController:1147", function() .py_method(p2, "doit"))
})
.py_try({
tf1 <- TransferFunction(s, .py_add(s, 1L), s)
})
.py_try({
par1 <- Parallel(p1, tf1)
})
.py_try({
ser1 <- Series(p1, tf1)
})
.py_try({
fed1 <- Feedback(p1, tf1)
})
.py_try({
.py_rec("lti:test_PIDController:1154", function() par1)
})
.py_try({
.py_rec("lti:test_PIDController:1155", function() ser1)
})
.py_try({
.py_rec("lti:test_PIDController:1156", function() fed1)
})
.py_try({
.py_rec("lti:test_PIDController:1157", function() .py_method(par1, "doit"))
})
.py_try({
.py_rec("lti:test_PIDController:1159", function() .py_method(ser1, "doit"))
})
.py_try({
.py_rec("lti:test_PIDController:1161", function() .py_method(fed1, "doit"))
}))---")
  expect_identical(got[["lti:test_PIDController:1128"]], "True")
  expect_identical(got[["lti:test_PIDController:1129"]], "True")
  expect_identical(got[["lti:test_PIDController:1132"]], "PIDController(kd*s^2 + ki*s*tf + ki + kp*s^2*tf + kp*s, s^2*tf + s, s)")
  expect_identical(got[["lti:test_PIDController:1133"]], "PIDController(kd*s^2 + ki + kp*s, s, s)")
  expect_identical(got[["lti:test_PIDController:1134"]], "kd*s^2 + ki*s*tf + ki + kp*s^2*tf + kp*s")
  expect_identical(got[["lti:test_PIDController:1135"]], "s^2*tf + s")
  expect_identical(got[["lti:test_PIDController:1136"]], "s")
  expect_identical(got[["lti:test_PIDController:1137"]], "kp")
  expect_identical(got[["lti:test_PIDController:1138"]], "ki")
  expect_identical(got[["lti:test_PIDController:1139"]], "kd")
  expect_identical(got[["lti:test_PIDController:1140"]], "tf")
  expect_identical(got[["lti:test_PIDController:1143"]], "TransferFunction(kd*s^2 + ki*s*tf + ki + kp*s^2*tf + kp*s, s^2*tf + s, s)")
  expect_identical(got[["lti:test_PIDController:1144"]], "True")
  expect_identical(got[["lti:test_PIDController:1145"]], "True")
  expect_identical(got[["lti:test_PIDController:1146"]], "False")
  expect_identical(got[["lti:test_PIDController:1147"]], "TransferFunction(kd*s^2 + ki + kp*s, s, s)")
  expect_identical(got[["lti:test_PIDController:1154"]], "Parallel(PIDController(kd*s^2 + ki*s*tf + ki + kp*s^2*tf + kp*s, s^2*tf + s, s), TransferFunction(s, s + 1, s))")
  expect_identical(got[["lti:test_PIDController:1155"]], "Series(PIDController(kd*s^2 + ki*s*tf + ki + kp*s^2*tf + kp*s, s^2*tf + s, s), TransferFunction(s, s + 1, s))")
  expect_identical(got[["lti:test_PIDController:1156"]], "Feedback(PIDController(kd*s^2 + ki*s*tf + ki + kp*s^2*tf + kp*s, s^2*tf + s, s), TransferFunction(s, s + 1, s), -1)")
  expect_identical(got[["lti:test_PIDController:1157"]], "TransferFunction(s*(s^2*tf + s) + (s + 1)*(kd*s^2 + ki*s*tf + ki + kp*s^2*tf + kp*s), (s + 1)*(s^2*tf + s), s)")
  expect_identical(got[["lti:test_PIDController:1159"]], "TransferFunction(s*(kd*s^2 + ki*s*tf + ki + kp*s^2*tf + kp*s), (s + 1)*(s^2*tf + s), s)")
  expect_identical(got[["lti:test_PIDController:1161"]], "TransferFunction((s + 1)*(s^2*tf + s)*(kd*s^2 + ki*s*tf + ki + kp*s^2*tf + kp*s), (s*(kd*s^2 + ki*s*tf + ki + kp*s^2*tf + kp*s) + (s + 1)*(s^2*tf + s))*(s^2*tf + s), s)")
})

test_that("lti test_Series_construction matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("a, x, b, c, s, g, d,    p, k, tau, zeta, wn, T, z"); a <- .py_getitem(.tmpu_, 0L); x <- .py_getitem(.tmpu_, 1L); b <- .py_getitem(.tmpu_, 2L); c <- .py_getitem(.tmpu_, 3L); s <- .py_getitem(.tmpu_, 4L); g <- .py_getitem(.tmpu_, 5L); d <- .py_getitem(.tmpu_, 6L); p <- .py_getitem(.tmpu_, 7L); k <- .py_getitem(.tmpu_, 8L); tau <- .py_getitem(.tmpu_, 9L); zeta <- .py_getitem(.tmpu_, 10L); wn <- .py_getitem(.tmpu_, 11L); T <- .py_getitem(.tmpu_, 12L); z <- .py_getitem(.tmpu_, 13L)
})
.py_try({
.tmpu_ <- .py_symbols("a0:4, b0:5, c0:4, d0:4"); a0 <- .py_getitem(.tmpu_, 0L); a1 <- .py_getitem(.tmpu_, 1L); a2 <- .py_getitem(.tmpu_, 2L); a3 <- .py_getitem(.tmpu_, 3L); b0 <- .py_getitem(.tmpu_, 4L); b1 <- .py_getitem(.tmpu_, 5L); b2 <- .py_getitem(.tmpu_, 6L); b3 <- .py_getitem(.tmpu_, 7L); b4 <- .py_getitem(.tmpu_, 8L); c0 <- .py_getitem(.tmpu_, 9L); c1 <- .py_getitem(.tmpu_, 10L); c2 <- .py_getitem(.tmpu_, 11L); c3 <- .py_getitem(.tmpu_, 12L); d0 <- .py_getitem(.tmpu_, 13L); d1 <- .py_getitem(.tmpu_, 14L); d2 <- .py_getitem(.tmpu_, 15L); d3 <- .py_getitem(.tmpu_, 16L)
})
.py_try({
TF1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
TF2 <- TransferFunction(k, 1L, s)
})
.py_try({
TF3 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
TF4 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(tau, s), 1L), s)
})
.py_try({
TF5 <- TransferFunction(k, p, s)
})
.py_try({
TF6 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(k, s), p), s)
})
p_are_floats_equals <- function(n1, n2) {
  "\n    Helper function to compare two floats.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return(.py_method(Float(n1), "is_same", Float(n2), .py_isclose))
}
p_are_complex_equals <- function(n1, n2) {
  "\n    Helper function to compare two complex numbers.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return((.py_truthy(p_are_floats_equals(re(n1), re(n2))) && .py_truthy(p_are_floats_equals(im(n1), im(n2)))))
}
.py_try({
tf <- TransferFunction(.py_sub(.py_add(.py_mul(a0, .py_pow(s, 3L)), .py_mul(a1, .py_pow(s, 2L))), .py_mul(a2, s)), .py_sub(.py_add(.py_mul(b0, .py_pow(p, 4L)), .py_mul(b1, .py_pow(p, 3L))), .py_mul(.py_mul(b2, s), p)), s)
})
.py_try({
tf2 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
tf3 <- TransferFunction(.py_sub(.py_add(.py_mul(a0, p), .py_pow(p, a1)), s), p, p)
})
.py_try({
tf4 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
inp <- Function("X_d")(s)
})
.py_try({
out <- Function("X")(s)
})
.py_try({
dtf <- DiscreteTransferFunction(.py_sub(.py_add(.py_mul(a0, .py_pow(s, 3L)), .py_mul(a1, .py_pow(s, 2L))), .py_mul(a2, s)), .py_sub(.py_add(.py_mul(b0, .py_pow(p, 4L)), .py_mul(b1, .py_pow(p, 3L))), .py_mul(.py_mul(b2, s), p)), s, .py_flt(0.1))
})
.py_try({
dtf2 <- DiscreteTransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s, .py_flt(0.1))
})
.py_try({
dtf3 <- DiscreteTransferFunction(.py_sub(.py_add(.py_mul(a0, p), .py_pow(p, a1)), s), p, p, .py_flt(0.1))
})
.py_try({
dtf4 <- DiscreteTransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s, .py_flt(0.1))
})
.py_try({
.py_rec("lti:test_Series_construction:1179", function() .py_attr(tf, "is_continuous"))
})
.py_try({
.py_rec("lti:test_Series_construction:1180", function() .py_attr(dtf, "is_continuous"))
})
.py_try({
s0 <- Series(tf, tf2)
})
.py_try({
.py_rec("lti:test_Series_construction:1184", function() .py_attr(s0, "args"))
})
.py_try({
.py_rec("lti:test_Series_construction:1185", function() .py_attr(s0, "var"))
})
.py_try({
s1 <- Series(Parallel(tf, .py_neg(tf2)), tf2)
})
.py_try({
.py_rec("lti:test_Series_construction:1188", function() .py_attr(s1, "args"))
})
.py_try({
.py_rec("lti:test_Series_construction:1189", function() .py_attr(s1, "var"))
})
.py_try({
tf3_ <- TransferFunction(inp, 1L, s)
})
.py_try({
tf4_ <- TransferFunction(.py_neg(out), 1L, s)
})
.py_try({
s2 <- Series(tf, Parallel(tf3_, tf4_), tf2)
})
.py_try({
.py_rec("lti:test_Series_construction:1194", function() .py_attr(s2, "args"))
})
.py_try({
s3 <- Series(tf, tf2, tf4)
})
.py_try({
.py_rec("lti:test_Series_construction:1197", function() .py_attr(s3, "args"))
})
.py_try({
s4 <- Series(tf3_, tf4_)
})
.py_try({
.py_rec("lti:test_Series_construction:1200", function() .py_attr(s4, "args"))
})
.py_try({
.py_rec("lti:test_Series_construction:1201", function() .py_attr(s4, "var"))
})
.py_try({
s6 <- Series(tf2, tf4, Parallel(tf2, .py_neg(tf)), tf4)
})
.py_try({
.py_rec("lti:test_Series_construction:1204", function() .py_attr(s6, "args"))
})
.py_try({
s7 <- Series(tf, tf2)
})
.py_try({
.py_rec("lti:test_Series_construction:1207", function() s0)
})
.py_try({
.py_rec("lti:test_Series_construction:1208", function() s0)
})
.py_try({
.py_rec("lti:test_Series_construction:1210", function() Series(tf, tf3))
})
.py_try({
.py_rec("lti:test_Series_construction:1211", function() Series(tf, tf2, tf3, tf4))
})
.py_try({
.py_rec("lti:test_Series_construction:1212", function() Series(.py_neg(tf3), tf2))
})
.py_try({
.py_rec("lti:test_Series_construction:1213", function() Series(2L, tf, tf4))
})
.py_try({
.py_rec("lti:test_Series_construction:1214", function() Series(.py_add(.py_pow(s, 2L), .py_mul(p, s)), tf3, tf2))
})
.py_try({
.py_rec("lti:test_Series_construction:1215", function() Series(tf3, .py_Matrix(list(1L, 2L, 3L, 4L))))
})
.py_try({
ds0 <- Series(dtf, dtf2)
})
.py_try({
.py_rec("lti:test_Series_construction:1219", function() .py_attr(ds0, "args"))
})
.py_try({
.py_rec("lti:test_Series_construction:1220", function() .py_attr(ds0, "var"))
})
.py_try({
.py_rec("lti:test_Series_construction:1221", function() .py_attr(ds0, "sampling_time"))
})
.py_try({
ds1 <- Series(Parallel(dtf, .py_neg(dtf2)), dtf2)
})
.py_try({
.py_rec("lti:test_Series_construction:1224", function() .py_attr(ds1, "args"))
})
.py_try({
.py_rec("lti:test_Series_construction:1225", function() .py_attr(ds1, "var"))
})
.py_try({
.py_rec("lti:test_Series_construction:1226", function() .py_attr(ds1, "sampling_time"))
})
.py_try({
dtf3_ <- DiscreteTransferFunction(inp, 1L, s, .py_flt(0.1))
})
.py_try({
dtf4_ <- DiscreteTransferFunction(.py_neg(out), 1L, s, .py_flt(0.1))
})
.py_try({
ds2 <- Series(dtf, Parallel(dtf3_, dtf4_), dtf2)
})
.py_try({
.py_rec("lti:test_Series_construction:1231", function() .py_attr(ds2, "args"))
})
.py_try({
.py_rec("lti:test_Series_construction:1232", function() .py_attr(ds2, "sampling_time"))
})
.py_try({
ds3 <- Series(dtf, dtf2, dtf4)
})
.py_try({
.py_rec("lti:test_Series_construction:1235", function() .py_attr(ds3, "args"))
})
.py_try({
.py_rec("lti:test_Series_construction:1236", function() .py_attr(ds3, "sampling_time"))
})
.py_try({
ds4 <- Series(dtf3_, dtf4_)
})
.py_try({
.py_rec("lti:test_Series_construction:1239", function() .py_attr(ds4, "args"))
})
.py_try({
.py_rec("lti:test_Series_construction:1240", function() .py_attr(ds4, "var"))
})
.py_try({
.py_rec("lti:test_Series_construction:1241", function() .py_attr(ds4, "sampling_time"))
})
.py_try({
ds6 <- Series(dtf2, dtf4, Parallel(dtf2, .py_neg(dtf)), dtf4)
})
.py_try({
.py_rec("lti:test_Series_construction:1244", function() .py_attr(ds6, "args"))
})
.py_try({
.py_rec("lti:test_Series_construction:1245", function() .py_attr(ds6, "sampling_time"))
})
.py_try({
ds7 <- Series(dtf, dtf2)
})
.py_try({
.py_rec("lti:test_Series_construction:1248", function() ds0)
})
.py_try({
.py_rec("lti:test_Series_construction:1249", function() ds0)
})
.py_try({
.py_rec("lti:test_Series_construction:1250", function() .py_attr(ds7, "sampling_time"))
})
.py_try({
dtf5 <- DiscreteTransferFunction(s, .py_sub(s, 1L), s, 1L)
})
.py_try({
.py_rec("lti:test_Series_construction:1253", function() Series(dtf, dtf3))
})
.py_try({
.py_rec("lti:test_Series_construction:1254", function() Series(dtf, dtf2, dtf3, dtf4))
})
.py_try({
.py_rec("lti:test_Series_construction:1255", function() Series(.py_neg(dtf3), dtf2))
})
.py_try({
.py_rec("lti:test_Series_construction:1256", function() Series(dtf5, dtf))
})
.py_try({
.py_rec("lti:test_Series_construction:1257", function() Series(dtf5, tf))
})
.py_try({
.py_rec("lti:test_Series_construction:1258", function() Series(2L, dtf, dtf4))
})
.py_try({
.py_rec("lti:test_Series_construction:1259", function() Series(.py_add(.py_pow(s, 2L), .py_mul(p, s)), dtf3, dtf2))
})
.py_try({
.py_rec("lti:test_Series_construction:1260", function() Series(dtf3, .py_Matrix(list(1L, 2L, 3L, 4L))))
}))---")
  expect_identical(got[["lti:test_Series_construction:1179"]], "True")
  expect_identical(got[["lti:test_Series_construction:1180"]], "False")
  expect_identical(got[["lti:test_Series_construction:1184"]], "(TransferFunction(a0*s^3 + a1*s^2 - a2*s, b0*p^4 + b1*p^3 - b2*p*s, s), TransferFunction(a2*p - s, a2*s + p, s))")
  expect_identical(got[["lti:test_Series_construction:1185"]], "s")
  expect_identical(got[["lti:test_Series_construction:1188"]], "(Parallel(TransferFunction(a0*s^3 + a1*s^2 - a2*s, b0*p^4 + b1*p^3 - b2*p*s, s), TransferFunction(-a2*p + s, a2*s + p, s)), TransferFunction(a2*p - s, a2*s + p, s))")
  expect_identical(got[["lti:test_Series_construction:1189"]], "s")
  expect_identical(got[["lti:test_Series_construction:1194"]], "(TransferFunction(a0*s^3 + a1*s^2 - a2*s, b0*p^4 + b1*p^3 - b2*p*s, s), Parallel(TransferFunction(X_d(s), 1, s), TransferFunction(-X(s), 1, s)), TransferFunction(a2*p - s, a2*s + p, s))")
  expect_identical(got[["lti:test_Series_construction:1197"]], "(TransferFunction(a0*s^3 + a1*s^2 - a2*s, b0*p^4 + b1*p^3 - b2*p*s, s), TransferFunction(a2*p - s, a2*s + p, s), TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s))")
  expect_identical(got[["lti:test_Series_construction:1200"]], "(TransferFunction(X_d(s), 1, s), TransferFunction(-X(s), 1, s))")
  expect_identical(got[["lti:test_Series_construction:1201"]], "s")
  expect_identical(got[["lti:test_Series_construction:1204"]], "(TransferFunction(a2*p - s, a2*s + p, s), TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s), Parallel(TransferFunction(a2*p - s, a2*s + p, s), TransferFunction(-a0*s^3 - a1*s^2 + a2*s, b0*p^4 + b1*p^3 - b2*p*s, s)), TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s))")
  expect_identical(got[["lti:test_Series_construction:1207"]], "Series(TransferFunction(a0*s^3 + a1*s^2 - a2*s, b0*p^4 + b1*p^3 - b2*p*s, s), TransferFunction(a2*p - s, a2*s + p, s))")
  expect_identical(got[["lti:test_Series_construction:1208"]], "Series(TransferFunction(a0*s^3 + a1*s^2 - a2*s, b0*p^4 + b1*p^3 - b2*p*s, s), TransferFunction(a2*p - s, a2*s + p, s))")
  expect_identical(got[["lti:test_Series_construction:1210"]], "ERR ValueError: \\nAll transfer functions should use the same complex variable of the\\nLaplace transform or z-transform. 2 different values found.")
  expect_identical(got[["lti:test_Series_construction:1211"]], "ERR ValueError: \\nAll transfer functions should use the same complex variable of the\\nLaplace transform or z-transform. 2 different values found.")
  expect_identical(got[["lti:test_Series_construction:1212"]], "ERR ValueError: \\nAll transfer functions should use the same complex variable of the\\nLaplace transform or z-transform. 2 different values found.")
  expect_identical(got[["lti:test_Series_construction:1213"]], "ERR TypeError: All arguments must be of type <class 'sympy.physics.control.lti.SISOLinearTimeInvariant'>.")
  expect_identical(got[["lti:test_Series_construction:1214"]], "ERR TypeError: All arguments must be of type <class 'sympy.physics.control.lti.SISOLinearTimeInvariant'>.")
  expect_identical(got[["lti:test_Series_construction:1215"]], "ERR TypeError: All arguments must be of type <class 'sympy.physics.control.lti.SISOLinearTimeInvariant'>.")
  expect_identical(got[["lti:test_Series_construction:1219"]], "(DiscreteTransferFunction(a0*s^3 + a1*s^2 - a2*s, b0*p^4 + b1*p^3 - b2*p*s, s, 0.1), DiscreteTransferFunction(a2*p - s, a2*s + p, s, 0.1))")
  expect_identical(got[["lti:test_Series_construction:1220"]], "s")
  expect_identical(got[["lti:test_Series_construction:1221"]], "0.100000000000000")
  expect_identical(got[["lti:test_Series_construction:1224"]], "(Parallel(DiscreteTransferFunction(a0*s^3 + a1*s^2 - a2*s, b0*p^4 + b1*p^3 - b2*p*s, s, 0.1), DiscreteTransferFunction(-a2*p + s, a2*s + p, s, 0.1)), DiscreteTransferFunction(a2*p - s, a2*s + p, s, 0.1))")
  expect_identical(got[["lti:test_Series_construction:1225"]], "s")
  expect_identical(got[["lti:test_Series_construction:1226"]], "0.100000000000000")
  expect_identical(got[["lti:test_Series_construction:1231"]], "(DiscreteTransferFunction(a0*s^3 + a1*s^2 - a2*s, b0*p^4 + b1*p^3 - b2*p*s, s, 0.1), Parallel(DiscreteTransferFunction(X_d(s), 1, s, 0.1), DiscreteTransferFunction(-X(s), 1, s, 0.1)), DiscreteTransferFunction(a2*p - s, a2*s + p, s, 0.1))")
  expect_identical(got[["lti:test_Series_construction:1232"]], "0.100000000000000")
  expect_identical(got[["lti:test_Series_construction:1235"]], "(DiscreteTransferFunction(a0*s^3 + a1*s^2 - a2*s, b0*p^4 + b1*p^3 - b2*p*s, s, 0.1), DiscreteTransferFunction(a2*p - s, a2*s + p, s, 0.1), DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.1))")
  expect_identical(got[["lti:test_Series_construction:1236"]], "0.100000000000000")
  expect_identical(got[["lti:test_Series_construction:1239"]], "(DiscreteTransferFunction(X_d(s), 1, s, 0.1), DiscreteTransferFunction(-X(s), 1, s, 0.1))")
  expect_identical(got[["lti:test_Series_construction:1240"]], "s")
  expect_identical(got[["lti:test_Series_construction:1241"]], "0.100000000000000")
  expect_identical(got[["lti:test_Series_construction:1244"]], "(DiscreteTransferFunction(a2*p - s, a2*s + p, s, 0.1), DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.1), Parallel(DiscreteTransferFunction(a2*p - s, a2*s + p, s, 0.1), DiscreteTransferFunction(-a0*s^3 - a1*s^2 + a2*s, b0*p^4 + b1*p^3 - b2*p*s, s, 0.1)), DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.1))")
  expect_identical(got[["lti:test_Series_construction:1245"]], "0.100000000000000")
  expect_identical(got[["lti:test_Series_construction:1248"]], "Series(DiscreteTransferFunction(a0*s^3 + a1*s^2 - a2*s, b0*p^4 + b1*p^3 - b2*p*s, s, 0.1), DiscreteTransferFunction(a2*p - s, a2*s + p, s, 0.1))")
  expect_identical(got[["lti:test_Series_construction:1249"]], "Series(DiscreteTransferFunction(a0*s^3 + a1*s^2 - a2*s, b0*p^4 + b1*p^3 - b2*p*s, s, 0.1), DiscreteTransferFunction(a2*p - s, a2*s + p, s, 0.1))")
  expect_identical(got[["lti:test_Series_construction:1250"]], "0.100000000000000")
  expect_identical(got[["lti:test_Series_construction:1253"]], "ERR ValueError: \\nAll transfer functions should use the same complex variable of the\\nLaplace transform or z-transform. 2 different values found.")
  expect_identical(got[["lti:test_Series_construction:1254"]], "ERR ValueError: \\nAll transfer functions should use the same complex variable of the\\nLaplace transform or z-transform. 2 different values found.")
  expect_identical(got[["lti:test_Series_construction:1255"]], "ERR ValueError: \\nAll transfer functions should use the same complex variable of the\\nLaplace transform or z-transform. 2 different values found.")
  expect_identical(got[["lti:test_Series_construction:1256"]], "ERR TypeError: \\n                Incompatible sampling times detected.\\n                All systems should have the same sampling time.\\n                Found sampling time: 1 and 0.100000000000000")
  expect_identical(got[["lti:test_Series_construction:1257"]], "ERR TypeError: \\n                Incompatible systems detected.\\n                All systems should be either continuous-time or discrete-time.\\n                Found: DiscreteTransferFunction(s, s - 1, s, 1) and TransferFunction(a0*s^3 + a1*s^2 - a2*s, b0*p^4 + b1*p^3 - b2*p*s, s)")
  expect_identical(got[["lti:test_Series_construction:1258"]], "ERR TypeError: All arguments must be of type <class 'sympy.physics.control.lti.SISOLinearTimeInvariant'>.")
  expect_identical(got[["lti:test_Series_construction:1259"]], "ERR TypeError: All arguments must be of type <class 'sympy.physics.control.lti.SISOLinearTimeInvariant'>.")
  expect_identical(got[["lti:test_Series_construction:1260"]], "ERR TypeError: All arguments must be of type <class 'sympy.physics.control.lti.SISOLinearTimeInvariant'>.")
})

test_that("lti test_MIMOSeries_construction matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("a, x, b, c, s, g, d,    p, k, tau, zeta, wn, T, z"); a <- .py_getitem(.tmpu_, 0L); x <- .py_getitem(.tmpu_, 1L); b <- .py_getitem(.tmpu_, 2L); c <- .py_getitem(.tmpu_, 3L); s <- .py_getitem(.tmpu_, 4L); g <- .py_getitem(.tmpu_, 5L); d <- .py_getitem(.tmpu_, 6L); p <- .py_getitem(.tmpu_, 7L); k <- .py_getitem(.tmpu_, 8L); tau <- .py_getitem(.tmpu_, 9L); zeta <- .py_getitem(.tmpu_, 10L); wn <- .py_getitem(.tmpu_, 11L); T <- .py_getitem(.tmpu_, 12L); z <- .py_getitem(.tmpu_, 13L)
})
.py_try({
.tmpu_ <- .py_symbols("a0:4, b0:5, c0:4, d0:4"); a0 <- .py_getitem(.tmpu_, 0L); a1 <- .py_getitem(.tmpu_, 1L); a2 <- .py_getitem(.tmpu_, 2L); a3 <- .py_getitem(.tmpu_, 3L); b0 <- .py_getitem(.tmpu_, 4L); b1 <- .py_getitem(.tmpu_, 5L); b2 <- .py_getitem(.tmpu_, 6L); b3 <- .py_getitem(.tmpu_, 7L); b4 <- .py_getitem(.tmpu_, 8L); c0 <- .py_getitem(.tmpu_, 9L); c1 <- .py_getitem(.tmpu_, 10L); c2 <- .py_getitem(.tmpu_, 11L); c3 <- .py_getitem(.tmpu_, 12L); d0 <- .py_getitem(.tmpu_, 13L); d1 <- .py_getitem(.tmpu_, 14L); d2 <- .py_getitem(.tmpu_, 15L); d3 <- .py_getitem(.tmpu_, 16L)
})
.py_try({
TF1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
TF2 <- TransferFunction(k, 1L, s)
})
.py_try({
TF3 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
TF4 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(tau, s), 1L), s)
})
.py_try({
TF5 <- TransferFunction(k, p, s)
})
.py_try({
TF6 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(k, s), p), s)
})
p_are_floats_equals <- function(n1, n2) {
  "\n    Helper function to compare two floats.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return(.py_method(Float(n1), "is_same", Float(n2), .py_isclose))
}
p_are_complex_equals <- function(n1, n2) {
  "\n    Helper function to compare two complex numbers.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return((.py_truthy(p_are_floats_equals(re(n1), re(n2))) && .py_truthy(p_are_floats_equals(im(n1), im(n2)))))
}
.py_try({
tf_1 <- TransferFunction(.py_sub(.py_add(.py_mul(a0, .py_pow(s, 3L)), .py_mul(a1, .py_pow(s, 2L))), .py_mul(a2, s)), .py_sub(.py_add(.py_mul(b0, .py_pow(p, 4L)), .py_mul(b1, .py_pow(p, 3L))), .py_mul(.py_mul(b2, s), p)), s)
})
.py_try({
tf_2 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
tf_3 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
dtf_1 <- DiscreteTransferFunction(.py_sub(.py_add(.py_mul(a0, .py_pow(s, 3L)), .py_mul(a1, .py_pow(s, 2L))), .py_mul(a2, s)), .py_sub(.py_add(.py_mul(b0, .py_pow(p, 4L)), .py_mul(b1, .py_pow(p, 3L))), .py_mul(.py_mul(b2, s), p)), s, 2L)
})
.py_try({
dtf_2 <- DiscreteTransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s, 2L)
})
.py_try({
dtf_3 <- DiscreteTransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s, 2L)
})
.py_try({
dTF1 <- DiscreteTransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s, 2L)
})
.py_try({
dTF2 <- DiscreteTransferFunction(k, 1L, s, 2L)
})
.py_try({
dTF3 <- DiscreteTransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s, 2L)
})
.py_try({
tfm_1 <- TransferFunctionMatrix(list(list(tf_1, tf_2, tf_3), list(.py_neg(tf_3), .py_neg(tf_2), tf_1)))
})
.py_try({
tfm_2 <- TransferFunctionMatrix(list(list(.py_neg(tf_2)), list(.py_neg(tf_2)), list(.py_neg(tf_3))))
})
.py_try({
tfm_3 <- TransferFunctionMatrix(list(list(.py_neg(tf_3))))
})
.py_try({
tfm_4 <- TransferFunctionMatrix(list(list(TF3), list(TF2), list(.py_neg(TF1))))
})
.py_try({
tfm_5 <- .py_method(TransferFunctionMatrix, "from_Matrix", .py_Matrix(list(.py_div(1L, p))), p)
})
.py_try({
dtfm_1 <- TransferFunctionMatrix(list(list(dtf_1, dtf_2, dtf_3), list(.py_neg(dtf_3), .py_neg(dtf_2), dtf_1)))
})
.py_try({
dtfm_2 <- TransferFunctionMatrix(list(list(.py_neg(dtf_2)), list(.py_neg(dtf_2)), list(.py_neg(dtf_3))))
})
.py_try({
dtfm_3 <- TransferFunctionMatrix(list(list(.py_neg(dtf_3))))
})
.py_try({
dtfm_4 <- TransferFunctionMatrix(list(list(dTF3), list(dTF2), list(.py_neg(dTF1))))
})
.py_try({
s8 <- MIMOSeries(tfm_2, tfm_1)
})
.py_try({
.py_rec("lti:test_MIMOSeries_construction:1290", function() .py_attr(s8, "args"))
})
.py_try({
.py_rec("lti:test_MIMOSeries_construction:1291", function() .py_attr(s8, "var"))
})
.py_try({
.py_rec("lti:test_MIMOSeries_construction:1292", function() (.py_truthy(.py_eq(.py_attr(s8, "shape"), .py_tup(.py_attr(s8, "num_outputs"), .py_attr(s8, "num_inputs")))) && .py_truthy(.py_eq(.py_tup(.py_attr(s8, "num_outputs"), .py_attr(s8, "num_inputs")), .py_tup(2L, 1L)))))
})
.py_try({
.py_rec("lti:test_MIMOSeries_construction:1293", function() .py_attr(s8, "is_continuous"))
})
.py_try({
s9 <- MIMOSeries(tfm_3, tfm_2, tfm_1)
})
.py_try({
.py_rec("lti:test_MIMOSeries_construction:1296", function() .py_attr(s9, "args"))
})
.py_try({
.py_rec("lti:test_MIMOSeries_construction:1297", function() .py_attr(s9, "var"))
})
.py_try({
.py_rec("lti:test_MIMOSeries_construction:1298", function() (.py_truthy(.py_eq(.py_attr(s9, "shape"), .py_tup(.py_attr(s9, "num_outputs"), .py_attr(s9, "num_inputs")))) && .py_truthy(.py_eq(.py_tup(.py_attr(s9, "num_outputs"), .py_attr(s9, "num_inputs")), .py_tup(2L, 1L)))))
})
.py_try({
s11 <- MIMOSeries(tfm_3, MIMOParallel(.py_neg(tfm_2), .py_neg(tfm_4)), tfm_1)
})
.py_try({
.py_rec("lti:test_MIMOSeries_construction:1301", function() .py_attr(s11, "args"))
})
.py_try({
.py_rec("lti:test_MIMOSeries_construction:1302", function() (.py_truthy(.py_eq(.py_attr(s11, "shape"), .py_tup(.py_attr(s11, "num_outputs"), .py_attr(s11, "num_inputs")))) && .py_truthy(.py_eq(.py_tup(.py_attr(s11, "num_outputs"), .py_attr(s11, "num_inputs")), .py_tup(2L, 1L)))))
})
.py_try({
ds8 <- MIMOSeries(dtfm_2, dtfm_1)
})
.py_try({
.py_rec("lti:test_MIMOSeries_construction:1306", function() .py_attr(ds8, "args"))
})
.py_try({
.py_rec("lti:test_MIMOSeries_construction:1307", function() .py_attr(ds8, "var"))
})
.py_try({
.py_rec("lti:test_MIMOSeries_construction:1308", function() (.py_truthy(.py_eq(.py_attr(ds8, "shape"), .py_tup(.py_attr(ds8, "num_outputs"), .py_attr(ds8, "num_inputs")))) && .py_truthy(.py_eq(.py_tup(.py_attr(ds8, "num_outputs"), .py_attr(ds8, "num_inputs")), .py_tup(2L, 1L)))))
})
.py_try({
.py_rec("lti:test_MIMOSeries_construction:1309", function() .py_attr(ds8, "is_continuous"))
})
.py_try({
.py_rec("lti:test_MIMOSeries_construction:1310", function() .py_attr(ds8, "sampling_time"))
})
.py_try({
ds9 <- MIMOSeries(dtfm_3, dtfm_2, dtfm_1)
})
.py_try({
.py_rec("lti:test_MIMOSeries_construction:1313", function() .py_attr(ds9, "args"))
})
.py_try({
.py_rec("lti:test_MIMOSeries_construction:1314", function() .py_attr(ds9, "var"))
})
.py_try({
.py_rec("lti:test_MIMOSeries_construction:1315", function() (.py_truthy(.py_eq(.py_attr(ds9, "shape"), .py_tup(.py_attr(ds9, "num_outputs"), .py_attr(ds9, "num_inputs")))) && .py_truthy(.py_eq(.py_tup(.py_attr(ds9, "num_outputs"), .py_attr(ds9, "num_inputs")), .py_tup(2L, 1L)))))
})
.py_try({
ds11 <- MIMOSeries(dtfm_3, MIMOParallel(.py_neg(dtfm_2), .py_neg(dtfm_4)), dtfm_1)
})
.py_try({
.py_rec("lti:test_MIMOSeries_construction:1318", function() .py_attr(ds11, "args"))
})
.py_try({
.py_rec("lti:test_MIMOSeries_construction:1319", function() (.py_truthy(.py_eq(.py_attr(ds11, "shape"), .py_tup(.py_attr(ds11, "num_outputs"), .py_attr(ds11, "num_inputs")))) && .py_truthy(.py_eq(.py_tup(.py_attr(ds11, "num_outputs"), .py_attr(ds11, "num_inputs")), .py_tup(2L, 1L)))))
})
.py_try({
.py_rec("lti:test_MIMOSeries_construction:1320", function() .py_attr(ds11, "sampling_time"))
})
.py_try({
.py_rec("lti:test_MIMOSeries_construction:1323", function() MIMOSeries())
})
.py_try({
.py_rec("lti:test_MIMOSeries_construction:1326", function() MIMOSeries(tfm_1, tf_1))
})
.py_try({
.py_rec("lti:test_MIMOSeries_construction:1330", function() MIMOSeries(tfm_1, tfm_2, .py_neg(tfm_1)))
})
.py_try({
.py_rec("lti:test_MIMOSeries_construction:1333", function() MIMOSeries(tfm_3, tfm_5))
})
.py_try({
.py_rec("lti:test_MIMOSeries_construction:1336", function() MIMOSeries(2L, tfm_2, tfm_3))
})
.py_try({
.py_rec("lti:test_MIMOSeries_construction:1337", function() MIMOSeries(.py_add(.py_pow(s, 2L), .py_mul(p, s)), .py_neg(tfm_2), tfm_3))
})
.py_try({
.py_rec("lti:test_MIMOSeries_construction:1338", function() MIMOSeries(.py_Matrix(list(.py_div(1L, p))), tfm_3))
})
.py_try({
dtfm_6 <- .py_method(TransferFunctionMatrix, "from_Matrix", .py_mul(.py_Matrix(eye(3L)), s), s, 12L)
})
.py_try({
.py_rec("lti:test_MIMOSeries_construction:1342", function() MIMOSeries(dtfm_6, dtfm_1))
})
.py_try({
.py_rec("lti:test_MIMOSeries_construction:1343", function() MIMOSeries(dtfm_2, tfm_1))
}))---")
  expect_identical(got[["lti:test_MIMOSeries_construction:1290"]], "(TransferFunctionMatrix(((TransferFunction(-a2*p + s, a2*s + p, s),), (TransferFunction(-a2*p + s, a2*s + p, s),), (TransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s),))), TransferFunctionMatrix(((TransferFunction(a0*s^3 + a1*s^2 - a2*s, b0*p^4 + b1*p^3 - b2*p*s, s), TransferFunction(a2*p - s, a2*s + p, s), TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s)), (TransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s), TransferFunction(-a2*p + s, a2*s + p, s), TransferFunction(a0*s^3 + a1*s^2 - a2*s, b0*p^4 + b1*p^3 - b2*p*s, s)))))")
  expect_identical(got[["lti:test_MIMOSeries_construction:1291"]], "s")
  expect_identical(got[["lti:test_MIMOSeries_construction:1292"]], "True")
  expect_identical(got[["lti:test_MIMOSeries_construction:1293"]], "True")
  expect_identical(got[["lti:test_MIMOSeries_construction:1296"]], "(TransferFunctionMatrix(((TransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s),),)), TransferFunctionMatrix(((TransferFunction(-a2*p + s, a2*s + p, s),), (TransferFunction(-a2*p + s, a2*s + p, s),), (TransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s),))), TransferFunctionMatrix(((TransferFunction(a0*s^3 + a1*s^2 - a2*s, b0*p^4 + b1*p^3 - b2*p*s, s), TransferFunction(a2*p - s, a2*s + p, s), TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s)), (TransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s), TransferFunction(-a2*p + s, a2*s + p, s), TransferFunction(a0*s^3 + a1*s^2 - a2*s, b0*p^4 + b1*p^3 - b2*p*s, s)))))")
  expect_identical(got[["lti:test_MIMOSeries_construction:1297"]], "s")
  expect_identical(got[["lti:test_MIMOSeries_construction:1298"]], "True")
  expect_identical(got[["lti:test_MIMOSeries_construction:1301"]], "(TransferFunctionMatrix(((TransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s),),)), MIMOParallel(TransferFunctionMatrix(((TransferFunction(a2*p - s, a2*s + p, s),), (TransferFunction(a2*p - s, a2*s + p, s),), (TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s),))), TransferFunctionMatrix(((TransferFunction(-a2*p + s, a2*s + p, s),), (TransferFunction(-k, 1, s),), (TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s),)))), TransferFunctionMatrix(((TransferFunction(a0*s^3 + a1*s^2 - a2*s, b0*p^4 + b1*p^3 - b2*p*s, s), TransferFunction(a2*p - s, a2*s + p, s), TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s)), (TransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s), TransferFunction(-a2*p + s, a2*s + p, s), TransferFunction(a0*s^3 + a1*s^2 - a2*s, b0*p^4 + b1*p^3 - b2*p*s, s)))))")
  expect_identical(got[["lti:test_MIMOSeries_construction:1302"]], "True")
  expect_identical(got[["lti:test_MIMOSeries_construction:1306"]], "(TransferFunctionMatrix(((DiscreteTransferFunction(-a2*p + s, a2*s + p, s, 2),), (DiscreteTransferFunction(-a2*p + s, a2*s + p, s, 2),), (DiscreteTransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s, 2),))), TransferFunctionMatrix(((DiscreteTransferFunction(a0*s^3 + a1*s^2 - a2*s, b0*p^4 + b1*p^3 - b2*p*s, s, 2), DiscreteTransferFunction(a2*p - s, a2*s + p, s, 2), DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 2)), (DiscreteTransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s, 2), DiscreteTransferFunction(-a2*p + s, a2*s + p, s, 2), DiscreteTransferFunction(a0*s^3 + a1*s^2 - a2*s, b0*p^4 + b1*p^3 - b2*p*s, s, 2)))))")
  expect_identical(got[["lti:test_MIMOSeries_construction:1307"]], "s")
  expect_identical(got[["lti:test_MIMOSeries_construction:1308"]], "True")
  expect_identical(got[["lti:test_MIMOSeries_construction:1309"]], "False")
  expect_identical(got[["lti:test_MIMOSeries_construction:1310"]], "2")
  expect_identical(got[["lti:test_MIMOSeries_construction:1313"]], "(TransferFunctionMatrix(((DiscreteTransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s, 2),),)), TransferFunctionMatrix(((DiscreteTransferFunction(-a2*p + s, a2*s + p, s, 2),), (DiscreteTransferFunction(-a2*p + s, a2*s + p, s, 2),), (DiscreteTransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s, 2),))), TransferFunctionMatrix(((DiscreteTransferFunction(a0*s^3 + a1*s^2 - a2*s, b0*p^4 + b1*p^3 - b2*p*s, s, 2), DiscreteTransferFunction(a2*p - s, a2*s + p, s, 2), DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 2)), (DiscreteTransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s, 2), DiscreteTransferFunction(-a2*p + s, a2*s + p, s, 2), DiscreteTransferFunction(a0*s^3 + a1*s^2 - a2*s, b0*p^4 + b1*p^3 - b2*p*s, s, 2)))))")
  expect_identical(got[["lti:test_MIMOSeries_construction:1314"]], "s")
  expect_identical(got[["lti:test_MIMOSeries_construction:1315"]], "True")
  expect_identical(got[["lti:test_MIMOSeries_construction:1318"]], "(TransferFunctionMatrix(((DiscreteTransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s, 2),),)), MIMOParallel(TransferFunctionMatrix(((DiscreteTransferFunction(a2*p - s, a2*s + p, s, 2),), (DiscreteTransferFunction(a2*p - s, a2*s + p, s, 2),), (DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 2),))), TransferFunctionMatrix(((DiscreteTransferFunction(-a2*p + s, a2*s + p, s, 2),), (DiscreteTransferFunction(-k, 1, s, 2),), (DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 2),)))), TransferFunctionMatrix(((DiscreteTransferFunction(a0*s^3 + a1*s^2 - a2*s, b0*p^4 + b1*p^3 - b2*p*s, s, 2), DiscreteTransferFunction(a2*p - s, a2*s + p, s, 2), DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 2)), (DiscreteTransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s, 2), DiscreteTransferFunction(-a2*p + s, a2*s + p, s, 2), DiscreteTransferFunction(a0*s^3 + a1*s^2 - a2*s, b0*p^4 + b1*p^3 - b2*p*s, s, 2)))))")
  expect_identical(got[["lti:test_MIMOSeries_construction:1319"]], "True")
  expect_identical(got[["lti:test_MIMOSeries_construction:1320"]], "2")
  expect_identical(got[["lti:test_MIMOSeries_construction:1323"]], "ERR ValueError: At least 1 argument must be passed.")
  expect_identical(got[["lti:test_MIMOSeries_construction:1326"]], "ERR TypeError: All arguments must be of type <class 'sympy.physics.control.lti.MIMOLinearTimeInvariant'>.")
  expect_identical(got[["lti:test_MIMOSeries_construction:1330"]], "ERR ValueError: \\nNumber of input signals do not match the number of output signals of\\nadjacent systems for some args.")
  expect_identical(got[["lti:test_MIMOSeries_construction:1333"]], "ERR ValueError: \\nAll transfer functions should use the same complex variable of the\\nLaplace transform or z-transform. 2 different values found.")
  expect_identical(got[["lti:test_MIMOSeries_construction:1336"]], "ERR TypeError: All arguments must be of type <class 'sympy.physics.control.lti.MIMOLinearTimeInvariant'>.")
  expect_identical(got[["lti:test_MIMOSeries_construction:1337"]], "ERR TypeError: All arguments must be of type <class 'sympy.physics.control.lti.MIMOLinearTimeInvariant'>.")
  expect_identical(got[["lti:test_MIMOSeries_construction:1338"]], "ERR TypeError: All arguments must be of type <class 'sympy.physics.control.lti.MIMOLinearTimeInvariant'>.")
  expect_identical(got[["lti:test_MIMOSeries_construction:1342"]], "ERR TypeError: \\n                Incompatible sampling times detected.\\n                All systems should have the same sampling time.\\n                Found sampling time: 12 and 2")
  expect_identical(got[["lti:test_MIMOSeries_construction:1343"]], "ERR TypeError: \\n                Incompatible systems detected.\\n                All systems should be either continuous-time or discrete-time.\\n                Found: TransferFunctionMatrix(((DiscreteTransferFunction(-a2*p + s, a2*s + p, s, 2),), (DiscreteTransferFunction(-a2*p + s, a2*s + p, s, 2),), (DiscreteTransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s, 2),))) and TransferFunctionMatrix(((TransferFunction(a0*s^3 + a1*s^2 - a2*s, b0*p^4 + b1*p^3 - b2*p*s, s), TransferFunction(a2*p - s, a2*s + p, s), TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s)), (TransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s), TransferFunction(-a2*p + s, a2*s + p, s), TransferFunction(a0*s^3 + a1*s^2 - a2*s, b0*p^4 + b1*p^3 - b2*p*s, s))))")
})

test_that("lti test_Series_functions matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("a, x, b, c, s, g, d,    p, k, tau, zeta, wn, T, z"); a <- .py_getitem(.tmpu_, 0L); x <- .py_getitem(.tmpu_, 1L); b <- .py_getitem(.tmpu_, 2L); c <- .py_getitem(.tmpu_, 3L); s <- .py_getitem(.tmpu_, 4L); g <- .py_getitem(.tmpu_, 5L); d <- .py_getitem(.tmpu_, 6L); p <- .py_getitem(.tmpu_, 7L); k <- .py_getitem(.tmpu_, 8L); tau <- .py_getitem(.tmpu_, 9L); zeta <- .py_getitem(.tmpu_, 10L); wn <- .py_getitem(.tmpu_, 11L); T <- .py_getitem(.tmpu_, 12L); z <- .py_getitem(.tmpu_, 13L)
})
.py_try({
.tmpu_ <- .py_symbols("a0:4, b0:5, c0:4, d0:4"); a0 <- .py_getitem(.tmpu_, 0L); a1 <- .py_getitem(.tmpu_, 1L); a2 <- .py_getitem(.tmpu_, 2L); a3 <- .py_getitem(.tmpu_, 3L); b0 <- .py_getitem(.tmpu_, 4L); b1 <- .py_getitem(.tmpu_, 5L); b2 <- .py_getitem(.tmpu_, 6L); b3 <- .py_getitem(.tmpu_, 7L); b4 <- .py_getitem(.tmpu_, 8L); c0 <- .py_getitem(.tmpu_, 9L); c1 <- .py_getitem(.tmpu_, 10L); c2 <- .py_getitem(.tmpu_, 11L); c3 <- .py_getitem(.tmpu_, 12L); d0 <- .py_getitem(.tmpu_, 13L); d1 <- .py_getitem(.tmpu_, 14L); d2 <- .py_getitem(.tmpu_, 15L); d3 <- .py_getitem(.tmpu_, 16L)
})
.py_try({
TF1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
TF2 <- TransferFunction(k, 1L, s)
})
.py_try({
TF3 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
TF4 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(tau, s), 1L), s)
})
.py_try({
TF5 <- TransferFunction(k, p, s)
})
.py_try({
TF6 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(k, s), p), s)
})
p_are_floats_equals <- function(n1, n2) {
  "\n    Helper function to compare two floats.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return(.py_method(Float(n1), "is_same", Float(n2), .py_isclose))
}
p_are_complex_equals <- function(n1, n2) {
  "\n    Helper function to compare two complex numbers.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return((.py_truthy(p_are_floats_equals(re(n1), re(n2))) && .py_truthy(p_are_floats_equals(im(n1), im(n2)))))
}
.py_try({
tf1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
tf2 <- TransferFunction(k, 1L, s)
})
.py_try({
tf3 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
tf4 <- TransferFunction(.py_sub(.py_add(.py_mul(a0, p), .py_pow(p, a1)), s), p, p)
})
.py_try({
tf5 <- TransferFunction(.py_sub(.py_add(.py_mul(a1, .py_pow(s, 2L)), .py_mul(a2, s)), a0), .py_add(s, a0), s)
})
.py_try({
dtf1 <- DiscreteTransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s, .py_flt(0.2))
})
.py_try({
dtf2 <- DiscreteTransferFunction(k, 1L, s, .py_flt(0.2))
})
.py_try({
dtf3 <- DiscreteTransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s, .py_flt(0.2))
})
.py_try({
dtf4 <- DiscreteTransferFunction(.py_sub(.py_add(.py_mul(a0, p), .py_pow(p, a1)), s), p, p, .py_flt(0.2))
})
.py_try({
dtf5 <- DiscreteTransferFunction(.py_sub(.py_add(.py_mul(a1, .py_pow(s, 2L)), .py_mul(a2, s)), a0), .py_add(s, a0), s, .py_flt(0.2))
})
.py_try({
.py_rec("lti:test_Series_functions:1360", function() (.py_truthy(.py_eq(.py_mul(.py_mul(tf1, tf2), tf3), Series(tf1, tf2, tf3))) && .py_truthy(.py_eq(Series(tf1, tf2, tf3), Series(Series(tf1, tf2), tf3))) && .py_truthy(.py_eq(Series(Series(tf1, tf2), tf3), Series(tf1, Series(tf2, tf3))))))
})
.py_try({
.py_rec("lti:test_Series_functions:1361", function() .py_mul(tf1, .py_add(tf2, tf3)))
})
.py_try({
.py_rec("lti:test_Series_functions:1362", function() .py_add(.py_mul(tf1, tf2), tf5))
})
.py_try({
.py_rec("lti:test_Series_functions:1363", function() .py_sub(.py_mul(tf1, tf2), tf5))
})
.py_try({
.py_rec("lti:test_Series_functions:1364", function() .py_add(.py_add(.py_mul(tf1, tf2), tf3), tf5))
})
.py_try({
.py_rec("lti:test_Series_functions:1365", function() .py_sub(.py_sub(.py_mul(tf1, tf2), tf3), tf5))
})
.py_try({
.py_rec("lti:test_Series_functions:1366", function() .py_add(.py_sub(.py_mul(tf1, tf2), tf3), tf5))
})
.py_try({
.py_rec("lti:test_Series_functions:1367", function() .py_add(.py_mul(tf1, tf2), .py_mul(tf3, tf5)))
})
.py_try({
.py_rec("lti:test_Series_functions:1368", function() .py_sub(.py_mul(tf1, tf2), .py_mul(tf3, tf5)))
})
.py_try({
.py_rec("lti:test_Series_functions:1369", function() .py_mul(.py_mul(.py_mul(tf2, tf3), .py_sub(tf2, tf1)), tf3))
})
.py_try({
.py_rec("lti:test_Series_functions:1370", function() .py_mul(.py_neg(tf1), tf2))
})
.py_try({
.py_rec("lti:test_Series_functions:1371", function() .py_neg(.py_mul(tf1, tf2)))
})
.py_try({
.py_rec("lti:test_Series_functions:1372", function() .py_mul(.py_mul(tf1, tf2), tf4))
})
.py_try({
.py_rec("lti:test_Series_functions:1373", function() .py_mul(tf1, .py_sub(tf2, tf4)))
})
.py_try({
.py_rec("lti:test_Series_functions:1374", function() .py_mul(tf3, .py_Matrix(list(1L, 2L, 3L))))
})
.py_try({
.py_rec("lti:test_Series_functions:1377", function() (.py_truthy(.py_eq(Series(tf1, tf2, evaluate = TRUE), .py_method(Series(tf1, tf2), "doit"))) && .py_truthy(.py_eq(.py_method(Series(tf1, tf2), "doit"), TransferFunction(k, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L)), s)))))
})
.py_try({
.py_rec("lti:test_Series_functions:1379", function() (.py_truthy(.py_eq(Series(tf1, tf2, Parallel(tf1, .py_neg(tf3)), evaluate = TRUE), .py_method(Series(tf1, tf2, Parallel(tf1, .py_neg(tf3))), "doit"))) && .py_truthy(.py_eq(.py_method(Series(tf1, tf2, Parallel(tf1, .py_neg(tf3))), "doit"), TransferFunction(.py_mul(k, .py_add(.py_add(.py_mul(a2, s), p), .py_mul(.py_add(.py_mul(.py_neg(a2), p), s), .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L))))), .py_mul(.py_add(.py_mul(a2, s), p), .py_pow(.py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L)), 2L)), s)))))
})
.py_try({
.py_rec("lti:test_Series_functions:1381", function() (.py_truthy(.py_eq(Series(tf2, tf1, .py_neg(tf3), evaluate = TRUE), .py_method(Series(tf2, tf1, .py_neg(tf3)), "doit"))) && .py_truthy(.py_eq(.py_method(Series(tf2, tf1, .py_neg(tf3)), "doit"), TransferFunction(.py_mul(k, .py_add(.py_mul(.py_neg(a2), p), s)), .py_mul(.py_add(.py_mul(a2, s), p), .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L))), s)))))
})
.py_try({
.py_rec("lti:test_Series_functions:1383", function() Series(tf1, .py_neg(tf2), evaluate = FALSE))
})
.py_try({
.py_rec("lti:test_Series_functions:1385", function() .py_method(Series(Parallel(tf1, tf2), Parallel(tf2, .py_neg(tf3))), "doit"))
})
.py_try({
.py_rec("lti:test_Series_functions:1387", function() .py_method(Series(.py_neg(tf1), .py_neg(tf2), .py_neg(tf3)), "doit"))
})
.py_try({
.py_rec("lti:test_Series_functions:1389", function() .py_neg(.py_method(Series(tf1, tf2, tf3), "doit")))
})
.py_try({
.py_rec("lti:test_Series_functions:1391", function() .py_method(Series(tf2, tf3, Parallel(tf2, .py_neg(tf1)), tf3), "doit"))
})
.py_try({
.py_rec("lti:test_Series_functions:1394", function() .py_method(Series(tf1, tf2), "rewrite", TransferFunction))
})
.py_try({
.py_rec("lti:test_Series_functions:1395", function() .py_method(Series(tf2, tf1, .py_neg(tf3)), "rewrite", TransferFunction))
})
.py_try({
S1 <- Series(Parallel(tf1, tf2), Parallel(tf2, .py_neg(tf3)))
})
.py_try({
.py_rec("lti:test_Series_functions:1399", function() .py_attr(S1, "is_proper"))
})
.py_try({
.py_rec("lti:test_Series_functions:1400", function() .py_attr(S1, "is_strictly_proper"))
})
.py_try({
.py_rec("lti:test_Series_functions:1401", function() .py_attr(S1, "is_biproper"))
})
.py_try({
S2 <- Series(tf1, tf2, tf3)
})
.py_try({
.py_rec("lti:test_Series_functions:1404", function() .py_attr(S2, "is_proper"))
})
.py_try({
.py_rec("lti:test_Series_functions:1405", function() .py_attr(S2, "is_strictly_proper"))
})
.py_try({
.py_rec("lti:test_Series_functions:1406", function() .py_attr(S2, "is_biproper"))
})
.py_try({
S3 <- Series(tf1, .py_neg(tf2), Parallel(tf1, .py_neg(tf3)))
})
.py_try({
.py_rec("lti:test_Series_functions:1409", function() .py_attr(S3, "is_proper"))
})
.py_try({
.py_rec("lti:test_Series_functions:1410", function() .py_attr(S3, "is_strictly_proper"))
})
.py_try({
.py_rec("lti:test_Series_functions:1411", function() .py_attr(S3, "is_biproper"))
})
.py_try({
.py_rec("lti:test_Series_functions:1414", function() (.py_truthy(.py_eq(.py_mul(.py_mul(dtf1, dtf2), dtf3), Series(dtf1, dtf2, dtf3))) && .py_truthy(.py_eq(Series(dtf1, dtf2, dtf3), Series(Series(dtf1, dtf2), dtf3))) && .py_truthy(.py_eq(Series(Series(dtf1, dtf2), dtf3), Series(dtf1, Series(dtf2, dtf3))))))
})
.py_try({
.py_rec("lti:test_Series_functions:1416", function() .py_mul(dtf1, .py_add(dtf2, dtf3)))
})
.py_try({
.py_rec("lti:test_Series_functions:1417", function() .py_add(.py_mul(dtf1, dtf2), dtf5))
})
.py_try({
.py_rec("lti:test_Series_functions:1418", function() .py_sub(.py_mul(dtf1, dtf2), dtf5))
})
.py_try({
.py_rec("lti:test_Series_functions:1419", function() .py_add(.py_add(.py_mul(dtf1, dtf2), dtf3), dtf5))
})
.py_try({
.py_rec("lti:test_Series_functions:1420", function() .py_sub(.py_sub(.py_mul(dtf1, dtf2), dtf3), dtf5))
})
.py_try({
.py_rec("lti:test_Series_functions:1421", function() .py_add(.py_sub(.py_mul(dtf1, dtf2), dtf3), dtf5))
})
.py_try({
.py_rec("lti:test_Series_functions:1422", function() .py_add(.py_mul(dtf1, dtf2), .py_mul(dtf3, dtf5)))
})
.py_try({
.py_rec("lti:test_Series_functions:1424", function() .py_sub(.py_mul(dtf1, dtf2), .py_mul(dtf3, dtf5)))
})
.py_try({
.py_rec("lti:test_Series_functions:1427", function() .py_mul(.py_mul(.py_mul(dtf2, dtf3), .py_sub(dtf2, dtf1)), dtf3))
})
.py_try({
.py_rec("lti:test_Series_functions:1429", function() .py_mul(.py_neg(dtf1), dtf2))
})
.py_try({
.py_rec("lti:test_Series_functions:1430", function() .py_neg(.py_mul(dtf1, dtf2)))
})
.py_try({
.py_rec("lti:test_Series_functions:1433", function() .py_mul(.py_mul(dtf1, dtf2), dtf4))
})
.py_try({
.py_rec("lti:test_Series_functions:1434", function() .py_mul(dtf1, .py_sub(dtf2, dtf4)))
})
.py_try({
.py_rec("lti:test_Series_functions:1435", function() .py_mul(dtf3, .py_Matrix(list(1L, 2L, 3L))))
})
.py_try({
.py_rec("lti:test_Series_functions:1438", function() (.py_truthy(.py_eq(Series(dtf1, dtf2, evaluate = TRUE), .py_method(Series(dtf1, dtf2), "doit"))) && .py_truthy(.py_eq(.py_method(Series(dtf1, dtf2), "doit"), DiscreteTransferFunction(k, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L)), s, .py_flt(0.2))))))
})
.py_try({
.py_rec("lti:test_Series_functions:1440", function() (.py_truthy(.py_eq(Series(dtf1, dtf2, Parallel(dtf1, .py_neg(dtf3)), evaluate = TRUE), .py_method(Series(dtf1, dtf2, Parallel(dtf1, .py_neg(dtf3))), "doit"))) && .py_truthy(.py_eq(.py_method(Series(dtf1, dtf2, Parallel(dtf1, .py_neg(dtf3))), "doit"), DiscreteTransferFunction(.py_mul(k, .py_add(.py_add(.py_mul(a2, s), p), .py_mul(.py_add(.py_mul(.py_neg(a2), p), s), .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L))))), .py_mul(.py_add(.py_mul(a2, s), p), .py_pow(.py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L)), 2L)), s, .py_flt(0.2))))))
})
.py_try({
.py_rec("lti:test_Series_functions:1445", function() (.py_truthy(.py_eq(Series(dtf2, dtf1, .py_neg(dtf3), evaluate = TRUE), .py_method(Series(dtf2, dtf1, .py_neg(dtf3)), "doit"))) && .py_truthy(.py_eq(.py_method(Series(dtf2, dtf1, .py_neg(dtf3)), "doit"), DiscreteTransferFunction(.py_mul(k, .py_add(.py_mul(.py_neg(a2), p), s)), .py_mul(.py_add(.py_mul(a2, s), p), .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L))), s, .py_flt(0.2))))))
})
.py_try({
.py_rec("lti:test_Series_functions:1449", function() Series(dtf1, .py_neg(dtf2), evaluate = FALSE))
})
.py_try({
.py_rec("lti:test_Series_functions:1451", function() .py_method(Series(Parallel(dtf1, dtf2), Parallel(dtf2, .py_neg(dtf3))), "doit"))
})
.py_try({
.py_rec("lti:test_Series_functions:1455", function() .py_method(Series(.py_neg(dtf1), .py_neg(dtf2), .py_neg(dtf3)), "doit"))
})
.py_try({
.py_rec("lti:test_Series_functions:1458", function() .py_neg(.py_method(Series(dtf1, dtf2, dtf3), "doit")))
})
.py_try({
.py_rec("lti:test_Series_functions:1461", function() .py_method(Series(dtf2, dtf3, Parallel(dtf2, .py_neg(dtf1)), dtf3), "doit"))
})
.py_try({
.py_rec("lti:test_Series_functions:1465", function() .py_method(Series(dtf1, dtf2), "rewrite", DiscreteTransferFunction))
})
.py_try({
.py_rec("lti:test_Series_functions:1467", function() .py_method(Series(dtf2, dtf1, .py_neg(dtf3)), "rewrite", DiscreteTransferFunction))
})
.py_try({
S1 <- Series(Parallel(dtf1, dtf2), Parallel(dtf2, .py_neg(dtf3)))
})
.py_try({
.py_rec("lti:test_Series_functions:1473", function() .py_attr(S1, "is_proper"))
})
.py_try({
.py_rec("lti:test_Series_functions:1474", function() .py_attr(S1, "is_strictly_proper"))
})
.py_try({
.py_rec("lti:test_Series_functions:1475", function() .py_attr(S1, "is_biproper"))
})
.py_try({
S2 <- Series(dtf1, dtf2, dtf3)
})
.py_try({
.py_rec("lti:test_Series_functions:1478", function() .py_attr(S2, "is_proper"))
})
.py_try({
.py_rec("lti:test_Series_functions:1479", function() .py_attr(S2, "is_strictly_proper"))
})
.py_try({
.py_rec("lti:test_Series_functions:1480", function() .py_attr(S2, "is_biproper"))
})
.py_try({
S3 <- Series(dtf1, .py_neg(dtf2), Parallel(dtf1, .py_neg(dtf3)))
})
.py_try({
.py_rec("lti:test_Series_functions:1483", function() .py_attr(S3, "is_proper"))
})
.py_try({
.py_rec("lti:test_Series_functions:1484", function() .py_attr(S3, "is_strictly_proper"))
})
.py_try({
.py_rec("lti:test_Series_functions:1485", function() .py_attr(S3, "is_biproper"))
})
.py_try({
.py_rec("lti:test_Series_functions:1488", function() .py_method(Series(dtf1, dtf2), "rewrite", TransferFunction))
})
.py_try({
.py_rec("lti:test_Series_functions:1490", function() .py_method(Series(tf1, tf2), "rewrite", DiscreteTransferFunction))
}))---")
  expect_identical(got[["lti:test_Series_functions:1360"]], "True")
  expect_identical(got[["lti:test_Series_functions:1361"]], "Series(TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s), Parallel(TransferFunction(k, 1, s), TransferFunction(a2*p - s, a2*s + p, s)))")
  expect_identical(got[["lti:test_Series_functions:1362"]], "Parallel(Series(TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s), TransferFunction(k, 1, s)), TransferFunction(-a0 + a1*s^2 + a2*s, a0 + s, s))")
  expect_identical(got[["lti:test_Series_functions:1363"]], "Parallel(Series(TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s), TransferFunction(k, 1, s)), TransferFunction(a0 - a1*s^2 - a2*s, a0 + s, s))")
  expect_identical(got[["lti:test_Series_functions:1364"]], "Parallel(Series(TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s), TransferFunction(k, 1, s)), TransferFunction(a2*p - s, a2*s + p, s), TransferFunction(-a0 + a1*s^2 + a2*s, a0 + s, s))")
  expect_identical(got[["lti:test_Series_functions:1365"]], "Parallel(Series(TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s), TransferFunction(k, 1, s)), TransferFunction(-a2*p + s, a2*s + p, s), TransferFunction(a0 - a1*s^2 - a2*s, a0 + s, s))")
  expect_identical(got[["lti:test_Series_functions:1366"]], "Parallel(Series(TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s), TransferFunction(k, 1, s)), TransferFunction(-a2*p + s, a2*s + p, s), TransferFunction(-a0 + a1*s^2 + a2*s, a0 + s, s))")
  expect_identical(got[["lti:test_Series_functions:1367"]], "Parallel(Series(TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s), TransferFunction(k, 1, s)), Series(TransferFunction(a2*p - s, a2*s + p, s), TransferFunction(-a0 + a1*s^2 + a2*s, a0 + s, s)))")
  expect_identical(got[["lti:test_Series_functions:1368"]], "Parallel(Series(TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s), TransferFunction(k, 1, s)), Series(TransferFunction(-1, 1, s), TransferFunction(a2*p - s, a2*s + p, s), TransferFunction(-a0 + a1*s^2 + a2*s, a0 + s, s)))")
  expect_identical(got[["lti:test_Series_functions:1369"]], "Series(TransferFunction(k, 1, s), TransferFunction(a2*p - s, a2*s + p, s), Parallel(TransferFunction(k, 1, s), TransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s)), TransferFunction(a2*p - s, a2*s + p, s))")
  expect_identical(got[["lti:test_Series_functions:1370"]], "Series(TransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s), TransferFunction(k, 1, s))")
  expect_identical(got[["lti:test_Series_functions:1371"]], "Series(TransferFunction(-1, 1, s), TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s), TransferFunction(k, 1, s))")
  expect_identical(got[["lti:test_Series_functions:1372"]], "ERR ValueError: \\nAll transfer functions should use the same complex variable of the\\nLaplace transform or z-transform. 2 different values found.")
  expect_identical(got[["lti:test_Series_functions:1373"]], "ERR ValueError: \\nAll the transfer functions should use the same complex variable of the\\nLaplace domain or z-domain.")
  expect_identical(got[["lti:test_Series_functions:1374"]], "ERR ValueError: <class 'sympy.physics.control.lti.TransferFunction'> cannot be multiplied with <class 'sympy.matrices.dense.MutableDenseMatrix'>.")
  expect_identical(got[["lti:test_Series_functions:1377"]], "True")
  expect_identical(got[["lti:test_Series_functions:1379"]], "True")
  expect_identical(got[["lti:test_Series_functions:1381"]], "True")
  expect_identical(got[["lti:test_Series_functions:1383"]], "Series(TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s), TransferFunction(-k, 1, s))")
  expect_identical(got[["lti:test_Series_functions:1385"]], "TransferFunction((k*(s^2 + 2*s*wn*zeta + wn^2) + 1)*(-a2*p + k*(a2*s + p) + s), (a2*s + p)*(s^2 + 2*s*wn*zeta + wn^2), s)")
  expect_identical(got[["lti:test_Series_functions:1387"]], "TransferFunction(k*(-a2*p + s), (a2*s + p)*(s^2 + 2*s*wn*zeta + wn^2), s)")
  expect_identical(got[["lti:test_Series_functions:1389"]], "TransferFunction(-k*(a2*p - s), (a2*s + p)*(s^2 + 2*s*wn*zeta + wn^2), s)")
  expect_identical(got[["lti:test_Series_functions:1391"]], "TransferFunction(k*(a2*p - s)^2*(k*(s^2 + 2*s*wn*zeta + wn^2) - 1), (a2*s + p)^2*(s^2 + 2*s*wn*zeta + wn^2), s)")
  expect_identical(got[["lti:test_Series_functions:1394"]], "TransferFunction(k, s^2 + 2*s*wn*zeta + wn^2, s)")
  expect_identical(got[["lti:test_Series_functions:1395"]], "TransferFunction(k*(-a2*p + s), (a2*s + p)*(s^2 + 2*s*wn*zeta + wn^2), s)")
  expect_identical(got[["lti:test_Series_functions:1399"]], "True")
  expect_identical(got[["lti:test_Series_functions:1400"]], "False")
  expect_identical(got[["lti:test_Series_functions:1401"]], "True")
  expect_identical(got[["lti:test_Series_functions:1404"]], "True")
  expect_identical(got[["lti:test_Series_functions:1405"]], "True")
  expect_identical(got[["lti:test_Series_functions:1406"]], "False")
  expect_identical(got[["lti:test_Series_functions:1409"]], "True")
  expect_identical(got[["lti:test_Series_functions:1410"]], "True")
  expect_identical(got[["lti:test_Series_functions:1411"]], "False")
  expect_identical(got[["lti:test_Series_functions:1414"]], "True")
  expect_identical(got[["lti:test_Series_functions:1416"]], "Series(DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.2), Parallel(DiscreteTransferFunction(k, 1, s, 0.2), DiscreteTransferFunction(a2*p - s, a2*s + p, s, 0.2)))")
  expect_identical(got[["lti:test_Series_functions:1417"]], "Parallel(Series(DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.2), DiscreteTransferFunction(k, 1, s, 0.2)), DiscreteTransferFunction(-a0 + a1*s^2 + a2*s, a0 + s, s, 0.2))")
  expect_identical(got[["lti:test_Series_functions:1418"]], "Parallel(Series(DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.2), DiscreteTransferFunction(k, 1, s, 0.2)), DiscreteTransferFunction(a0 - a1*s^2 - a2*s, a0 + s, s, 0.2))")
  expect_identical(got[["lti:test_Series_functions:1419"]], "Parallel(Series(DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.2), DiscreteTransferFunction(k, 1, s, 0.2)), DiscreteTransferFunction(a2*p - s, a2*s + p, s, 0.2), DiscreteTransferFunction(-a0 + a1*s^2 + a2*s, a0 + s, s, 0.2))")
  expect_identical(got[["lti:test_Series_functions:1420"]], "Parallel(Series(DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.2), DiscreteTransferFunction(k, 1, s, 0.2)), DiscreteTransferFunction(-a2*p + s, a2*s + p, s, 0.2), DiscreteTransferFunction(a0 - a1*s^2 - a2*s, a0 + s, s, 0.2))")
  expect_identical(got[["lti:test_Series_functions:1421"]], "Parallel(Series(DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.2), DiscreteTransferFunction(k, 1, s, 0.2)), DiscreteTransferFunction(-a2*p + s, a2*s + p, s, 0.2), DiscreteTransferFunction(-a0 + a1*s^2 + a2*s, a0 + s, s, 0.2))")
  expect_identical(got[["lti:test_Series_functions:1422"]], "Parallel(Series(DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.2), DiscreteTransferFunction(k, 1, s, 0.2)), Series(DiscreteTransferFunction(a2*p - s, a2*s + p, s, 0.2), DiscreteTransferFunction(-a0 + a1*s^2 + a2*s, a0 + s, s, 0.2)))")
  expect_identical(got[["lti:test_Series_functions:1424"]], "Parallel(Series(DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.2), DiscreteTransferFunction(k, 1, s, 0.2)), Series(DiscreteTransferFunction(-1, 1, s, 0.2), DiscreteTransferFunction(a2*p - s, a2*s + p, s, 0.2), DiscreteTransferFunction(-a0 + a1*s^2 + a2*s, a0 + s, s, 0.2)))")
  expect_identical(got[["lti:test_Series_functions:1427"]], "Series(DiscreteTransferFunction(k, 1, s, 0.2), DiscreteTransferFunction(a2*p - s, a2*s + p, s, 0.2), Parallel(DiscreteTransferFunction(k, 1, s, 0.2), DiscreteTransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s, 0.2)), DiscreteTransferFunction(a2*p - s, a2*s + p, s, 0.2))")
  expect_identical(got[["lti:test_Series_functions:1429"]], "Series(DiscreteTransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s, 0.2), DiscreteTransferFunction(k, 1, s, 0.2))")
  expect_identical(got[["lti:test_Series_functions:1430"]], "Series(DiscreteTransferFunction(-1, 1, s, 0.2), DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.2), DiscreteTransferFunction(k, 1, s, 0.2))")
  expect_identical(got[["lti:test_Series_functions:1433"]], "ERR ValueError: \\nAll transfer functions should use the same complex variable of the\\nLaplace transform or z-transform. 2 different values found.")
  expect_identical(got[["lti:test_Series_functions:1434"]], "ERR ValueError: \\nAll the transfer functions should use the same complex variable of the\\nLaplace domain or z-domain.")
  expect_identical(got[["lti:test_Series_functions:1435"]], "ERR ValueError: <class 'sympy.physics.control.lti.DiscreteTransferFunction'> cannot be multiplied with <class 'sympy.matrices.dense.MutableDenseMatrix'>.")
  expect_identical(got[["lti:test_Series_functions:1438"]], "True")
  expect_identical(got[["lti:test_Series_functions:1440"]], "True")
  expect_identical(got[["lti:test_Series_functions:1445"]], "True")
  expect_identical(got[["lti:test_Series_functions:1449"]], "Series(DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.2), DiscreteTransferFunction(-k, 1, s, 0.2))")
  expect_identical(got[["lti:test_Series_functions:1451"]], "DiscreteTransferFunction((k*(s^2 + 2*s*wn*zeta + wn^2) + 1)*(-a2*p + k*(a2*s + p) + s), (a2*s + p)*(s^2 + 2*s*wn*zeta + wn^2), s, 0.2)")
  expect_identical(got[["lti:test_Series_functions:1455"]], "DiscreteTransferFunction(k*(-a2*p + s), (a2*s + p)*(s^2 + 2*s*wn*zeta + wn^2), s, 0.2)")
  expect_identical(got[["lti:test_Series_functions:1458"]], "DiscreteTransferFunction(-k*(a2*p - s), (a2*s + p)*(s^2 + 2*s*wn*zeta + wn^2), s, 0.2)")
  expect_identical(got[["lti:test_Series_functions:1461"]], "DiscreteTransferFunction(k*(a2*p - s)^2*(k*(s^2 + 2*s*wn*zeta + wn^2) - 1), (a2*s + p)^2*(s^2 + 2*s*wn*zeta + wn^2), s, 0.2)")
  expect_identical(got[["lti:test_Series_functions:1465"]], "DiscreteTransferFunction(k, s^2 + 2*s*wn*zeta + wn^2, s, 0.2)")
  expect_identical(got[["lti:test_Series_functions:1467"]], "DiscreteTransferFunction(k*(-a2*p + s), (a2*s + p)*(s^2 + 2*s*wn*zeta + wn^2), s, 0.2)")
  expect_identical(got[["lti:test_Series_functions:1473"]], "True")
  expect_identical(got[["lti:test_Series_functions:1474"]], "False")
  expect_identical(got[["lti:test_Series_functions:1475"]], "True")
  expect_identical(got[["lti:test_Series_functions:1478"]], "True")
  expect_identical(got[["lti:test_Series_functions:1479"]], "True")
  expect_identical(got[["lti:test_Series_functions:1480"]], "False")
  expect_identical(got[["lti:test_Series_functions:1483"]], "True")
  expect_identical(got[["lti:test_Series_functions:1484"]], "True")
  expect_identical(got[["lti:test_Series_functions:1485"]], "False")
  expect_identical(got[["lti:test_Series_functions:1488"]], "ERR TypeError: \\n                    Cannot rewrite a discrete-time Series object as a\\n                    TransferFunction.")
  expect_identical(got[["lti:test_Series_functions:1490"]], "ERR TypeError: \\n                    Cannot rewrite a continuous-time Series object as a\\n                    DiscreteTransferFunction.")
})

test_that("lti test_MIMOSeries_functions matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("a, x, b, c, s, g, d,    p, k, tau, zeta, wn, T, z"); a <- .py_getitem(.tmpu_, 0L); x <- .py_getitem(.tmpu_, 1L); b <- .py_getitem(.tmpu_, 2L); c <- .py_getitem(.tmpu_, 3L); s <- .py_getitem(.tmpu_, 4L); g <- .py_getitem(.tmpu_, 5L); d <- .py_getitem(.tmpu_, 6L); p <- .py_getitem(.tmpu_, 7L); k <- .py_getitem(.tmpu_, 8L); tau <- .py_getitem(.tmpu_, 9L); zeta <- .py_getitem(.tmpu_, 10L); wn <- .py_getitem(.tmpu_, 11L); T <- .py_getitem(.tmpu_, 12L); z <- .py_getitem(.tmpu_, 13L)
})
.py_try({
.tmpu_ <- .py_symbols("a0:4, b0:5, c0:4, d0:4"); a0 <- .py_getitem(.tmpu_, 0L); a1 <- .py_getitem(.tmpu_, 1L); a2 <- .py_getitem(.tmpu_, 2L); a3 <- .py_getitem(.tmpu_, 3L); b0 <- .py_getitem(.tmpu_, 4L); b1 <- .py_getitem(.tmpu_, 5L); b2 <- .py_getitem(.tmpu_, 6L); b3 <- .py_getitem(.tmpu_, 7L); b4 <- .py_getitem(.tmpu_, 8L); c0 <- .py_getitem(.tmpu_, 9L); c1 <- .py_getitem(.tmpu_, 10L); c2 <- .py_getitem(.tmpu_, 11L); c3 <- .py_getitem(.tmpu_, 12L); d0 <- .py_getitem(.tmpu_, 13L); d1 <- .py_getitem(.tmpu_, 14L); d2 <- .py_getitem(.tmpu_, 15L); d3 <- .py_getitem(.tmpu_, 16L)
})
.py_try({
TF1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
TF2 <- TransferFunction(k, 1L, s)
})
.py_try({
TF3 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
TF4 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(tau, s), 1L), s)
})
.py_try({
TF5 <- TransferFunction(k, p, s)
})
.py_try({
TF6 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(k, s), p), s)
})
p_are_floats_equals <- function(n1, n2) {
  "\n    Helper function to compare two floats.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return(.py_method(Float(n1), "is_same", Float(n2), .py_isclose))
}
p_are_complex_equals <- function(n1, n2) {
  "\n    Helper function to compare two complex numbers.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return((.py_truthy(p_are_floats_equals(re(n1), re(n2))) && .py_truthy(p_are_floats_equals(im(n1), im(n2)))))
}
.py_try({
tfm1 <- TransferFunctionMatrix(list(list(TF1, TF2, TF3), list(.py_neg(TF3), .py_neg(TF2), TF1)))
})
.py_try({
tfm2 <- TransferFunctionMatrix(list(list(.py_neg(TF1)), list(.py_neg(TF2)), list(.py_neg(TF3))))
})
.py_try({
tfm3 <- TransferFunctionMatrix(list(list(.py_neg(TF1))))
})
.py_try({
tfm4 <- TransferFunctionMatrix(list(list(.py_neg(TF2), .py_neg(TF3)), list(.py_neg(TF1), TF2)))
})
.py_try({
tfm5 <- TransferFunctionMatrix(list(list(TF2, .py_neg(TF2)), list(.py_neg(TF3), .py_neg(TF2))))
})
.py_try({
tfm6 <- TransferFunctionMatrix(list(list(.py_neg(TF3)), list(TF1)))
})
.py_try({
tfm7 <- TransferFunctionMatrix(list(list(TF1), list(.py_neg(TF2))))
})
.py_try({
dTF1 <- DiscreteTransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s, .py_flt(0.001))
})
.py_try({
dTF2 <- DiscreteTransferFunction(k, 1L, s, .py_flt(0.001))
})
.py_try({
dTF3 <- DiscreteTransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s, .py_flt(0.001))
})
.py_try({
dtfm1 <- TransferFunctionMatrix(list(list(dTF1, dTF2, dTF3), list(.py_neg(dTF3), .py_neg(dTF2), dTF1)))
})
.py_try({
dtfm2 <- TransferFunctionMatrix(list(list(.py_neg(dTF1)), list(.py_neg(dTF2)), list(.py_neg(dTF3))))
})
.py_try({
dtfm3 <- TransferFunctionMatrix(list(list(.py_neg(dTF1))))
})
.py_try({
dtfm4 <- TransferFunctionMatrix(list(list(.py_neg(dTF2), .py_neg(dTF3)), list(.py_neg(dTF1), dTF2)))
})
.py_try({
dtfm5 <- TransferFunctionMatrix(list(list(dTF2, .py_neg(dTF2)), list(.py_neg(dTF3), .py_neg(dTF2))))
})
.py_try({
dtfm6 <- TransferFunctionMatrix(list(list(.py_neg(dTF3)), list(dTF1)))
})
.py_try({
dtfm7 <- TransferFunctionMatrix(list(list(dTF1), list(.py_neg(dTF2))))
})
.py_try({
.py_rec("lti:test_MIMOSeries_functions:1515", function() .py_add(.py_mul(tfm1, tfm2), tfm6))
})
.py_try({
.py_rec("lti:test_MIMOSeries_functions:1516", function() .py_add(.py_add(.py_mul(tfm1, tfm2), tfm7), tfm6))
})
.py_try({
.py_rec("lti:test_MIMOSeries_functions:1517", function() .py_sub(.py_sub(.py_mul(tfm1, tfm2), tfm6), tfm7))
})
.py_try({
.py_rec("lti:test_MIMOSeries_functions:1518", function() .py_add(.py_mul(tfm4, tfm5), .py_sub(tfm4, tfm5)))
})
.py_try({
.py_rec("lti:test_MIMOSeries_functions:1519", function() .py_add(.py_mul(tfm4, .py_neg(tfm6)), .py_mul(.py_neg(tfm4), tfm6)))
})
.py_try({
.py_rec("lti:test_MIMOSeries_functions:1521", function() .py_add(.py_mul(tfm1, tfm2), TF1))
})
.py_try({
.py_rec("lti:test_MIMOSeries_functions:1522", function() .py_add(.py_mul(tfm1, tfm2), a0))
})
.py_try({
.py_rec("lti:test_MIMOSeries_functions:1523", function() .py_sub(.py_mul(tfm4, tfm6), .py_sub(s, 1L)))
})
.py_try({
.py_rec("lti:test_MIMOSeries_functions:1524", function() .py_sub(.py_mul(tfm4, .py_neg(tfm6)), 8L))
})
.py_try({
.py_rec("lti:test_MIMOSeries_functions:1525", function() .py_add(.py_add(-1L, .py_pow(p, 5L)), .py_mul(tfm1, tfm2)))
})
.py_try({
.py_rec("lti:test_MIMOSeries_functions:1529", function() .py_add(.py_mul(.py_neg(tfm1), tfm2), tfm4))
})
.py_try({
.py_rec("lti:test_MIMOSeries_functions:1530", function() .py_add(.py_sub(.py_mul(tfm1, tfm2), tfm4), tfm5))
})
.py_try({
.py_rec("lti:test_MIMOSeries_functions:1531", function() .py_sub(.py_mul(tfm1, tfm2), .py_mul(tfm4, tfm5)))
})
.py_try({
.py_rec("lti:test_MIMOSeries_functions:1533", function() .py_mul(.py_mul(tfm1, tfm2), .py_neg(tfm3)))
})
.py_try({
.py_rec("lti:test_MIMOSeries_functions:1534", function() .py_mul(.py_mul(tfm1, .py_neg(tfm2)), tfm3))
})
.py_try({
.py_rec("lti:test_MIMOSeries_functions:1538", function() .py_mul(.py_mul(tfm4, tfm5), TF1))
})
.py_try({
.py_rec("lti:test_MIMOSeries_functions:1539", function() .py_mul(.py_mul(tfm4, tfm5), a1))
})
.py_try({
.py_rec("lti:test_MIMOSeries_functions:1540", function() .py_mul(.py_mul(tfm4, .py_neg(tfm5)), .py_sub(s, 2L)))
})
.py_try({
.py_rec("lti:test_MIMOSeries_functions:1541", function() .py_mul(.py_mul(tfm5, tfm4), 9L))
})
.py_try({
.py_rec("lti:test_MIMOSeries_functions:1542", function() .py_mul(.py_mul(.py_add(.py_neg(.py_pow(p, 3L)), 1L), tfm5), tfm4))
})
.py_try({
.py_rec("lti:test_MIMOSeries_functions:1545", function() (.py_truthy(.py_eq(MIMOSeries(tfm2, tfm1, evaluate = TRUE), .py_method(MIMOSeries(tfm2, tfm1), "doit"))) && .py_truthy(.py_eq(.py_method(MIMOSeries(tfm2, tfm1), "doit"), TransferFunctionMatrix(.py_tup(.py_tup(TransferFunction(.py_sub(.py_add(.py_mul(.py_mul(.py_neg(.py_pow(k, 2L)), .py_pow(.py_add(.py_mul(a2, s), p), 2L)), .py_pow(.py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L)), 2L)), .py_mul(.py_mul(.py_add(.py_mul(.py_neg(a2), p), s), .py_sub(.py_mul(a2, p), s)), .py_pow(.py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L)), 2L))), .py_pow(.py_add(.py_mul(a2, s), p), 2L)), .py_mul(.py_pow(.py_add(.py_mul(a2, s), p), 2L), .py_pow(.py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L)), 2L)), s)), .py_tup(TransferFunction(.py_add(.py_add(.py_mul(.py_mul(.py_pow(k, 2L), .py_pow(.py_add(.py_mul(a2, s), p), 2L)), .py_pow(.py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L)), 2L)), .py_mul(.py_mul(.py_add(.py_mul(.py_neg(a2), p), s), .py_add(.py_mul(a2, s), p)), .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L)))), .py_mul(.py_mul(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p)), .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L)))), .py_mul(.py_pow(.py_add(.py_mul(a2, s), p), 2L), .py_pow(.py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L)), 2L)), s))))))))
})
.py_try({
mat_1 <- .py_Matrix(list(list(.py_div(1L, .py_add(1L, s)), .py_div(.py_add(1L, s), .py_pow(.py_add(.py_add(1L, .py_pow(s, 2L)), .py_mul(2L, s)), 3L)))))
})
.py_try({
mat_2 <- .py_Matrix(list(list(.py_add(1L, s)), list(.py_div(.py_pow(.py_add(.py_add(1L, .py_pow(s, 2L)), .py_mul(2L, s)), 3L), .py_add(1L, s)))))
})
.py_try({
.tmpu_ <- .py_tup(.py_method(TransferFunctionMatrix, "from_Matrix", mat_1, s), .py_method(TransferFunctionMatrix, "from_Matrix", mat_2, s)); tm_1 <- .py_getitem(.tmpu_, 0L); tm_2 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
.py_rec("lti:test_MIMOSeries_functions:1555", function() .py_method(MIMOSeries(tm_2, tm_1), "doit"))
})
.py_try({
.py_rec("lti:test_MIMOSeries_functions:1557", function() .py_method(.py_method(MIMOSeries(tm_2, tm_1), "doit"), "simplify"))
})
.py_try({
.py_rec("lti:test_MIMOSeries_functions:1560", function() (.py_truthy(.py_eq(MIMOSeries(.py_neg(tfm3), .py_neg(tfm2), tfm1, evaluate = TRUE), .py_method(MIMOSeries(.py_neg(tfm3), .py_neg(tfm2), tfm1), "doit"))) && .py_truthy(.py_eq(.py_method(MIMOSeries(.py_neg(tfm3), .py_neg(tfm2), tfm1), "doit"), TransferFunctionMatrix(.py_tup(.py_tup(TransferFunction(.py_add(.py_add(.py_mul(.py_mul(.py_pow(k, 2L), .py_pow(.py_add(.py_mul(a2, s), p), 2L)), .py_pow(.py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L)), 2L)), .py_mul(.py_pow(.py_sub(.py_mul(a2, p), s), 2L), .py_pow(.py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L)), 2L))), .py_pow(.py_add(.py_mul(a2, s), p), 2L)), .py_mul(.py_pow(.py_add(.py_mul(a2, s), p), 2L), .py_pow(.py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L)), 3L)), s)), .py_tup(TransferFunction(.py_add(.py_add(.py_mul(.py_mul(.py_neg(.py_pow(k, 2L)), .py_pow(.py_add(.py_mul(a2, s), p), 2L)), .py_pow(.py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L)), 2L)), .py_mul(.py_mul(.py_add(.py_mul(.py_neg(a2), p), s), .py_add(.py_mul(a2, s), p)), .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L)))), .py_mul(.py_mul(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p)), .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L)))), .py_mul(.py_pow(.py_add(.py_mul(a2, s), p), 2L), .py_pow(.py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L)), 3L)), s))))))))
})
.py_try({
.py_rec("lti:test_MIMOSeries_functions:1566", function() (.py_truthy(.py_eq(MIMOSeries(MIMOParallel(tfm4, tfm5), tfm5, evaluate = TRUE), .py_method(MIMOSeries(MIMOParallel(tfm4, tfm5), tfm5), "doit"))) && .py_truthy(.py_eq(.py_method(MIMOSeries(MIMOParallel(tfm4, tfm5), tfm5), "doit"), TransferFunctionMatrix(.py_tup(.py_tup(TransferFunction(.py_mul(.py_neg(k), .py_add(.py_sub(.py_mul(.py_neg(a2), s), p), .py_mul(.py_add(.py_mul(.py_neg(a2), p), s), .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L))))), .py_mul(.py_add(.py_mul(a2, s), p), .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L))), s), TransferFunction(.py_mul(k, .py_add(.py_sub(.py_mul(.py_neg(a2), p), .py_mul(k, .py_add(.py_mul(a2, s), p))), s)), .py_add(.py_mul(a2, s), p), s)), .py_tup(TransferFunction(.py_mul(.py_neg(k), .py_add(.py_sub(.py_mul(.py_neg(a2), s), p), .py_mul(.py_add(.py_mul(.py_neg(a2), p), s), .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L))))), .py_mul(.py_add(.py_mul(a2, s), p), .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L))), s), TransferFunction(.py_mul(.py_add(.py_mul(.py_neg(a2), p), s), .py_add(.py_sub(.py_mul(.py_neg(a2), p), .py_mul(k, .py_add(.py_mul(a2, s), p))), s)), .py_pow(.py_add(.py_mul(a2, s), p), 2L), s)))))) && .py_truthy(.py_eq(TransferFunctionMatrix(.py_tup(.py_tup(TransferFunction(.py_mul(.py_neg(k), .py_add(.py_sub(.py_mul(.py_neg(a2), s), p), .py_mul(.py_add(.py_mul(.py_neg(a2), p), s), .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L))))), .py_mul(.py_add(.py_mul(a2, s), p), .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L))), s), TransferFunction(.py_mul(k, .py_add(.py_sub(.py_mul(.py_neg(a2), p), .py_mul(k, .py_add(.py_mul(a2, s), p))), s)), .py_add(.py_mul(a2, s), p), s)), .py_tup(TransferFunction(.py_mul(.py_neg(k), .py_add(.py_sub(.py_mul(.py_neg(a2), s), p), .py_mul(.py_add(.py_mul(.py_neg(a2), p), s), .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L))))), .py_mul(.py_add(.py_mul(a2, s), p), .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L))), s), TransferFunction(.py_mul(.py_add(.py_mul(.py_neg(a2), p), s), .py_add(.py_sub(.py_mul(.py_neg(a2), p), .py_mul(k, .py_add(.py_mul(a2, s), p))), s)), .py_pow(.py_add(.py_mul(a2, s), p), 2L), s)))), .py_method(MIMOSeries(MIMOParallel(tfm4, tfm5), tfm5), "rewrite", TransferFunctionMatrix)))))
})
.py_try({
.py_rec("lti:test_MIMOSeries_functions:1583", function() .py_add(.py_mul(dtfm1, dtfm2), dtfm6))
})
.py_try({
.py_rec("lti:test_MIMOSeries_functions:1584", function() .py_add(.py_add(.py_mul(dtfm1, dtfm2), dtfm7), dtfm6))
})
.py_try({
.py_rec("lti:test_MIMOSeries_functions:1586", function() .py_sub(.py_sub(.py_mul(dtfm1, dtfm2), dtfm6), dtfm7))
})
.py_try({
.py_rec("lti:test_MIMOSeries_functions:1588", function() .py_add(.py_mul(dtfm4, dtfm5), .py_sub(dtfm4, dtfm5)))
})
.py_try({
.py_rec("lti:test_MIMOSeries_functions:1590", function() .py_add(.py_mul(dtfm4, .py_neg(dtfm6)), .py_mul(.py_neg(dtfm4), dtfm6)))
})
.py_try({
.py_rec("lti:test_MIMOSeries_functions:1593", function() .py_add(.py_mul(dtfm1, dtfm2), dTF1))
})
.py_try({
.py_rec("lti:test_MIMOSeries_functions:1594", function() .py_add(.py_mul(dtfm1, dtfm2), a0))
})
.py_try({
.py_rec("lti:test_MIMOSeries_functions:1595", function() .py_sub(.py_mul(dtfm4, dtfm6), .py_sub(s, 1L)))
})
.py_try({
.py_rec("lti:test_MIMOSeries_functions:1596", function() .py_sub(.py_mul(dtfm4, .py_neg(dtfm6)), 8L))
})
.py_try({
.py_rec("lti:test_MIMOSeries_functions:1597", function() .py_add(.py_add(-1L, .py_pow(p, 5L)), .py_mul(dtfm1, dtfm2)))
})
.py_try({
.py_rec("lti:test_MIMOSeries_functions:1601", function() .py_add(.py_mul(.py_neg(dtfm1), dtfm2), dtfm4))
})
.py_try({
.py_rec("lti:test_MIMOSeries_functions:1602", function() .py_add(.py_sub(.py_mul(dtfm1, dtfm2), dtfm4), dtfm5))
})
.py_try({
.py_rec("lti:test_MIMOSeries_functions:1603", function() .py_sub(.py_mul(dtfm1, dtfm2), .py_mul(dtfm4, dtfm5)))
})
.py_try({
.py_rec("lti:test_MIMOSeries_functions:1605", function() .py_mul(.py_mul(dtfm1, dtfm2), .py_neg(dtfm3)))
})
.py_try({
.py_rec("lti:test_MIMOSeries_functions:1606", function() .py_mul(.py_mul(dtfm1, .py_neg(dtfm2)), dtfm3))
})
.py_try({
.py_rec("lti:test_MIMOSeries_functions:1610", function() .py_mul(.py_mul(dtfm4, dtfm5), dTF1))
})
.py_try({
.py_rec("lti:test_MIMOSeries_functions:1611", function() .py_mul(.py_mul(dtfm4, dtfm5), a1))
})
.py_try({
.py_rec("lti:test_MIMOSeries_functions:1612", function() .py_mul(.py_mul(dtfm4, .py_neg(dtfm5)), .py_sub(s, 2L)))
})
.py_try({
.py_rec("lti:test_MIMOSeries_functions:1613", function() .py_mul(.py_mul(dtfm5, dtfm4), 9L))
})
.py_try({
.py_rec("lti:test_MIMOSeries_functions:1614", function() .py_mul(.py_mul(.py_add(.py_neg(.py_pow(p, 3L)), 1L), dtfm5), dtfm4))
})
.py_try({
.py_rec("lti:test_MIMOSeries_functions:1617", function() (.py_truthy(.py_eq(MIMOSeries(dtfm2, dtfm1, evaluate = TRUE), .py_method(MIMOSeries(dtfm2, dtfm1), "doit"))) && .py_truthy(.py_eq(.py_method(MIMOSeries(dtfm2, dtfm1), "doit"), TransferFunctionMatrix(.py_tup(.py_tup(DiscreteTransferFunction(.py_sub(.py_add(.py_mul(.py_mul(.py_neg(.py_pow(k, 2L)), .py_pow(.py_add(.py_mul(a2, s), p), 2L)), .py_pow(.py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L)), 2L)), .py_mul(.py_mul(.py_add(.py_mul(.py_neg(a2), p), s), .py_sub(.py_mul(a2, p), s)), .py_pow(.py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L)), 2L))), .py_pow(.py_add(.py_mul(a2, s), p), 2L)), .py_mul(.py_pow(.py_add(.py_mul(a2, s), p), 2L), .py_pow(.py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L)), 2L)), s, .py_flt(0.001))), .py_tup(DiscreteTransferFunction(.py_add(.py_add(.py_mul(.py_mul(.py_pow(k, 2L), .py_pow(.py_add(.py_mul(a2, s), p), 2L)), .py_pow(.py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L)), 2L)), .py_mul(.py_mul(.py_add(.py_mul(.py_neg(a2), p), s), .py_add(.py_mul(a2, s), p)), .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L)))), .py_mul(.py_mul(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p)), .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L)))), .py_mul(.py_pow(.py_add(.py_mul(a2, s), p), 2L), .py_pow(.py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L)), 2L)), s, .py_flt(0.001)))))))))
})
.py_try({
dmat_1 <- .py_Matrix(list(list(.py_div(1L, .py_add(1L, s)), .py_div(.py_add(1L, s), .py_pow(.py_add(.py_add(1L, .py_pow(s, 2L)), .py_mul(2L, s)), 3L)))))
})
.py_try({
dmat_2 <- .py_Matrix(list(list(.py_add(1L, s)), list(.py_div(.py_pow(.py_add(.py_add(1L, .py_pow(s, 2L)), .py_mul(2L, s)), 3L), .py_add(1L, s)))))
})
.py_try({
dtm_1 <- .py_method(TransferFunctionMatrix, "from_Matrix", dmat_1, s, .py_flt(0.001))
})
.py_try({
dtm_2 <- .py_method(TransferFunctionMatrix, "from_Matrix", dmat_2, s, .py_flt(0.001))
})
.py_try({
.py_rec("lti:test_MIMOSeries_functions:1638", function() .py_method(MIMOSeries(dtm_2, dtm_1), "doit"))
})
.py_try({
.py_rec("lti:test_MIMOSeries_functions:1642", function() .py_method(.py_method(MIMOSeries(dtm_2, dtm_1), "doit"), "simplify"))
})
.py_try({
.py_rec("lti:test_MIMOSeries_functions:1646", function() (.py_truthy(.py_eq(MIMOSeries(.py_neg(dtfm3), .py_neg(dtfm2), dtfm1, evaluate = TRUE), .py_method(MIMOSeries(.py_neg(dtfm3), .py_neg(dtfm2), dtfm1), "doit"))) && .py_truthy(.py_eq(.py_method(MIMOSeries(.py_neg(dtfm3), .py_neg(dtfm2), dtfm1), "doit"), TransferFunctionMatrix(.py_tup(.py_tup(DiscreteTransferFunction(.py_add(.py_add(.py_mul(.py_mul(.py_pow(k, 2L), .py_pow(.py_add(.py_mul(a2, s), p), 2L)), .py_pow(.py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L)), 2L)), .py_mul(.py_pow(.py_sub(.py_mul(a2, p), s), 2L), .py_pow(.py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L)), 2L))), .py_pow(.py_add(.py_mul(a2, s), p), 2L)), .py_mul(.py_pow(.py_add(.py_mul(a2, s), p), 2L), .py_pow(.py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L)), 3L)), s, .py_flt(0.001))), .py_tup(DiscreteTransferFunction(.py_add(.py_add(.py_mul(.py_mul(.py_neg(.py_pow(k, 2L)), .py_pow(.py_add(.py_mul(a2, s), p), 2L)), .py_pow(.py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L)), 2L)), .py_mul(.py_mul(.py_add(.py_mul(.py_neg(a2), p), s), .py_add(.py_mul(a2, s), p)), .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L)))), .py_mul(.py_mul(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p)), .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L)))), .py_mul(.py_pow(.py_add(.py_mul(a2, s), p), 2L), .py_pow(.py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L)), 3L)), s, .py_flt(0.001)))))))))
})
.py_try({
.py_rec("lti:test_MIMOSeries_functions:1661", function() (.py_truthy(.py_eq(MIMOSeries(MIMOParallel(dtfm4, dtfm5), dtfm5, evaluate = TRUE), .py_method(MIMOSeries(MIMOParallel(dtfm4, dtfm5), dtfm5), "doit"))) && .py_truthy(.py_eq(.py_method(MIMOSeries(MIMOParallel(dtfm4, dtfm5), dtfm5), "doit"), TransferFunctionMatrix(.py_tup(.py_tup(DiscreteTransferFunction(.py_mul(.py_neg(k), .py_add(.py_sub(.py_mul(.py_neg(a2), s), p), .py_mul(.py_add(.py_mul(.py_neg(a2), p), s), .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L))))), .py_mul(.py_add(.py_mul(a2, s), p), .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L))), s, .py_flt(0.001)), DiscreteTransferFunction(.py_mul(k, .py_add(.py_sub(.py_mul(.py_neg(a2), p), .py_mul(k, .py_add(.py_mul(a2, s), p))), s)), .py_add(.py_mul(a2, s), p), s, .py_flt(0.001))), .py_tup(DiscreteTransferFunction(.py_mul(.py_neg(k), .py_add(.py_sub(.py_mul(.py_neg(a2), s), p), .py_mul(.py_add(.py_mul(.py_neg(a2), p), s), .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L))))), .py_mul(.py_add(.py_mul(a2, s), p), .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L))), s, .py_flt(0.001)), DiscreteTransferFunction(.py_mul(.py_add(.py_mul(.py_neg(a2), p), s), .py_add(.py_sub(.py_mul(.py_neg(a2), p), .py_mul(k, .py_add(.py_mul(a2, s), p))), s)), .py_pow(.py_add(.py_mul(a2, s), p), 2L), s, .py_flt(0.001))))))) && .py_truthy(.py_eq(TransferFunctionMatrix(.py_tup(.py_tup(DiscreteTransferFunction(.py_mul(.py_neg(k), .py_add(.py_sub(.py_mul(.py_neg(a2), s), p), .py_mul(.py_add(.py_mul(.py_neg(a2), p), s), .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L))))), .py_mul(.py_add(.py_mul(a2, s), p), .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L))), s, .py_flt(0.001)), DiscreteTransferFunction(.py_mul(k, .py_add(.py_sub(.py_mul(.py_neg(a2), p), .py_mul(k, .py_add(.py_mul(a2, s), p))), s)), .py_add(.py_mul(a2, s), p), s, .py_flt(0.001))), .py_tup(DiscreteTransferFunction(.py_mul(.py_neg(k), .py_add(.py_sub(.py_mul(.py_neg(a2), s), p), .py_mul(.py_add(.py_mul(.py_neg(a2), p), s), .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L))))), .py_mul(.py_add(.py_mul(a2, s), p), .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L))), s, .py_flt(0.001)), DiscreteTransferFunction(.py_mul(.py_add(.py_mul(.py_neg(a2), p), s), .py_add(.py_sub(.py_mul(.py_neg(a2), p), .py_mul(k, .py_add(.py_mul(a2, s), p))), s)), .py_pow(.py_add(.py_mul(a2, s), p), 2L), s, .py_flt(0.001))))), .py_method(MIMOSeries(MIMOParallel(dtfm4, dtfm5), dtfm5), "rewrite", TransferFunctionMatrix)))))
}))---")
  expect_identical(got[["lti:test_MIMOSeries_functions:1515"]], "MIMOParallel(MIMOSeries(TransferFunctionMatrix(((TransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s),), (TransferFunction(-k, 1, s),), (TransferFunction(-a2*p + s, a2*s + p, s),))), TransferFunctionMatrix(((TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s), TransferFunction(k, 1, s), TransferFunction(a2*p - s, a2*s + p, s)), (TransferFunction(-a2*p + s, a2*s + p, s), TransferFunction(-k, 1, s), TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s))))), TransferFunctionMatrix(((TransferFunction(-a2*p + s, a2*s + p, s),), (TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s),))))")
  expect_identical(got[["lti:test_MIMOSeries_functions:1516"]], "MIMOParallel(MIMOSeries(TransferFunctionMatrix(((TransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s),), (TransferFunction(-k, 1, s),), (TransferFunction(-a2*p + s, a2*s + p, s),))), TransferFunctionMatrix(((TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s), TransferFunction(k, 1, s), TransferFunction(a2*p - s, a2*s + p, s)), (TransferFunction(-a2*p + s, a2*s + p, s), TransferFunction(-k, 1, s), TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s))))), TransferFunctionMatrix(((TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s),), (TransferFunction(-k, 1, s),))), TransferFunctionMatrix(((TransferFunction(-a2*p + s, a2*s + p, s),), (TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s),))))")
  expect_identical(got[["lti:test_MIMOSeries_functions:1517"]], "MIMOParallel(MIMOSeries(TransferFunctionMatrix(((TransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s),), (TransferFunction(-k, 1, s),), (TransferFunction(-a2*p + s, a2*s + p, s),))), TransferFunctionMatrix(((TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s), TransferFunction(k, 1, s), TransferFunction(a2*p - s, a2*s + p, s)), (TransferFunction(-a2*p + s, a2*s + p, s), TransferFunction(-k, 1, s), TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s))))), TransferFunctionMatrix(((TransferFunction(a2*p - s, a2*s + p, s),), (TransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s),))), TransferFunctionMatrix(((TransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s),), (TransferFunction(k, 1, s),))))")
  expect_identical(got[["lti:test_MIMOSeries_functions:1518"]], "MIMOParallel(MIMOSeries(TransferFunctionMatrix(((TransferFunction(k, 1, s), TransferFunction(-k, 1, s)), (TransferFunction(-a2*p + s, a2*s + p, s), TransferFunction(-k, 1, s)))), TransferFunctionMatrix(((TransferFunction(-k, 1, s), TransferFunction(-a2*p + s, a2*s + p, s)), (TransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s), TransferFunction(k, 1, s))))), TransferFunctionMatrix(((TransferFunction(-k, 1, s), TransferFunction(-a2*p + s, a2*s + p, s)), (TransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s), TransferFunction(k, 1, s)))), TransferFunctionMatrix(((TransferFunction(-k, 1, s), TransferFunction(k, 1, s)), (TransferFunction(a2*p - s, a2*s + p, s), TransferFunction(k, 1, s)))))")
  expect_identical(got[["lti:test_MIMOSeries_functions:1519"]], "MIMOParallel(MIMOSeries(TransferFunctionMatrix(((TransferFunction(a2*p - s, a2*s + p, s),), (TransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s),))), TransferFunctionMatrix(((TransferFunction(-k, 1, s), TransferFunction(-a2*p + s, a2*s + p, s)), (TransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s), TransferFunction(k, 1, s))))), MIMOSeries(TransferFunctionMatrix(((TransferFunction(-a2*p + s, a2*s + p, s),), (TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s),))), TransferFunctionMatrix(((TransferFunction(k, 1, s), TransferFunction(a2*p - s, a2*s + p, s)), (TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s), TransferFunction(-k, 1, s))))))")
  expect_identical(got[["lti:test_MIMOSeries_functions:1521"]], "ERR ValueError: <class 'sympy.physics.control.lti.TransferFunction'> cannot be added with <class 'sympy.physics.control.lti.MIMOSeries'>.")
  expect_identical(got[["lti:test_MIMOSeries_functions:1522"]], "ERR TypeError: unsupported operand type(s) for +: 'MIMOSeries' and 'Symbol'")
  expect_identical(got[["lti:test_MIMOSeries_functions:1523"]], "ERR TypeError: unsupported operand type(s) for -: 'MIMOSeries' and 'Add'")
  expect_identical(got[["lti:test_MIMOSeries_functions:1524"]], "ERR TypeError: unsupported operand type(s) for -: 'MIMOSeries' and 'int'")
  expect_identical(got[["lti:test_MIMOSeries_functions:1525"]], "ERR TypeError: unsupported operand type(s) for +: 'Add' and 'MIMOSeries'")
  expect_identical(got[["lti:test_MIMOSeries_functions:1529"]], "ERR TypeError: Shape of all the args is not equal.")
  expect_identical(got[["lti:test_MIMOSeries_functions:1530"]], "ERR TypeError: Shape of all the args is not equal.")
  expect_identical(got[["lti:test_MIMOSeries_functions:1531"]], "ERR TypeError: Shape of all the args is not equal.")
  expect_identical(got[["lti:test_MIMOSeries_functions:1533"]], "MIMOSeries(TransferFunctionMatrix(((TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s),),)), TransferFunctionMatrix(((TransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s),), (TransferFunction(-k, 1, s),), (TransferFunction(-a2*p + s, a2*s + p, s),))), TransferFunctionMatrix(((TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s), TransferFunction(k, 1, s), TransferFunction(a2*p - s, a2*s + p, s)), (TransferFunction(-a2*p + s, a2*s + p, s), TransferFunction(-k, 1, s), TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s)))))")
  expect_identical(got[["lti:test_MIMOSeries_functions:1534"]], "MIMOSeries(TransferFunctionMatrix(((TransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s),),)), TransferFunctionMatrix(((TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s),), (TransferFunction(k, 1, s),), (TransferFunction(a2*p - s, a2*s + p, s),))), TransferFunctionMatrix(((TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s), TransferFunction(k, 1, s), TransferFunction(a2*p - s, a2*s + p, s)), (TransferFunction(-a2*p + s, a2*s + p, s), TransferFunction(-k, 1, s), TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s)))))")
  expect_identical(got[["lti:test_MIMOSeries_functions:1538"]], "ERR ValueError: <class 'sympy.physics.control.lti.TransferFunction'> cannot be multiplied with <class 'sympy.physics.control.lti.MIMOSeries'>.")
  expect_identical(got[["lti:test_MIMOSeries_functions:1539"]], "ERR TypeError: unsupported operand type(s) for *: 'MIMOSeries' and 'Symbol'")
  expect_identical(got[["lti:test_MIMOSeries_functions:1540"]], "ERR TypeError: unsupported operand type(s) for *: 'MIMOSeries' and 'Add'")
  expect_identical(got[["lti:test_MIMOSeries_functions:1541"]], "ERR TypeError: unsupported operand type(s) for *: 'MIMOSeries' and 'int'")
  expect_identical(got[["lti:test_MIMOSeries_functions:1542"]], "ERR TypeError: unsupported operand type(s) for *: 'Add' and 'TransferFunctionMatrix'")
  expect_identical(got[["lti:test_MIMOSeries_functions:1545"]], "True")
  expect_identical(got[["lti:test_MIMOSeries_functions:1555"]], "TransferFunctionMatrix(((TransferFunction(2*(s + 1)^2*(s^2 + 2*s + 1)^3, (s + 1)^2*(s^2 + 2*s + 1)^3, s),),))")
  expect_identical(got[["lti:test_MIMOSeries_functions:1557"]], "TransferFunctionMatrix(((TransferFunction(2, 1, s),),))")
  expect_identical(got[["lti:test_MIMOSeries_functions:1560"]], "True")
  expect_identical(got[["lti:test_MIMOSeries_functions:1566"]], "True")
  expect_identical(got[["lti:test_MIMOSeries_functions:1583"]], "MIMOParallel(MIMOSeries(TransferFunctionMatrix(((DiscreteTransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s, 0.001),), (DiscreteTransferFunction(-k, 1, s, 0.001),), (DiscreteTransferFunction(-a2*p + s, a2*s + p, s, 0.001),))), TransferFunctionMatrix(((DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.001), DiscreteTransferFunction(k, 1, s, 0.001), DiscreteTransferFunction(a2*p - s, a2*s + p, s, 0.001)), (DiscreteTransferFunction(-a2*p + s, a2*s + p, s, 0.001), DiscreteTransferFunction(-k, 1, s, 0.001), DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.001))))), TransferFunctionMatrix(((DiscreteTransferFunction(-a2*p + s, a2*s + p, s, 0.001),), (DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.001),))))")
  expect_identical(got[["lti:test_MIMOSeries_functions:1584"]], "MIMOParallel(MIMOSeries(TransferFunctionMatrix(((DiscreteTransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s, 0.001),), (DiscreteTransferFunction(-k, 1, s, 0.001),), (DiscreteTransferFunction(-a2*p + s, a2*s + p, s, 0.001),))), TransferFunctionMatrix(((DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.001), DiscreteTransferFunction(k, 1, s, 0.001), DiscreteTransferFunction(a2*p - s, a2*s + p, s, 0.001)), (DiscreteTransferFunction(-a2*p + s, a2*s + p, s, 0.001), DiscreteTransferFunction(-k, 1, s, 0.001), DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.001))))), TransferFunctionMatrix(((DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.001),), (DiscreteTransferFunction(-k, 1, s, 0.001),))), TransferFunctionMatrix(((DiscreteTransferFunction(-a2*p + s, a2*s + p, s, 0.001),), (DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.001),))))")
  expect_identical(got[["lti:test_MIMOSeries_functions:1586"]], "MIMOParallel(MIMOSeries(TransferFunctionMatrix(((DiscreteTransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s, 0.001),), (DiscreteTransferFunction(-k, 1, s, 0.001),), (DiscreteTransferFunction(-a2*p + s, a2*s + p, s, 0.001),))), TransferFunctionMatrix(((DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.001), DiscreteTransferFunction(k, 1, s, 0.001), DiscreteTransferFunction(a2*p - s, a2*s + p, s, 0.001)), (DiscreteTransferFunction(-a2*p + s, a2*s + p, s, 0.001), DiscreteTransferFunction(-k, 1, s, 0.001), DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.001))))), TransferFunctionMatrix(((DiscreteTransferFunction(a2*p - s, a2*s + p, s, 0.001),), (DiscreteTransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s, 0.001),))), TransferFunctionMatrix(((DiscreteTransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s, 0.001),), (DiscreteTransferFunction(k, 1, s, 0.001),))))")
  expect_identical(got[["lti:test_MIMOSeries_functions:1588"]], "MIMOParallel(MIMOSeries(TransferFunctionMatrix(((DiscreteTransferFunction(k, 1, s, 0.001), DiscreteTransferFunction(-k, 1, s, 0.001)), (DiscreteTransferFunction(-a2*p + s, a2*s + p, s, 0.001), DiscreteTransferFunction(-k, 1, s, 0.001)))), TransferFunctionMatrix(((DiscreteTransferFunction(-k, 1, s, 0.001), DiscreteTransferFunction(-a2*p + s, a2*s + p, s, 0.001)), (DiscreteTransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s, 0.001), DiscreteTransferFunction(k, 1, s, 0.001))))), TransferFunctionMatrix(((DiscreteTransferFunction(-k, 1, s, 0.001), DiscreteTransferFunction(-a2*p + s, a2*s + p, s, 0.001)), (DiscreteTransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s, 0.001), DiscreteTransferFunction(k, 1, s, 0.001)))), TransferFunctionMatrix(((DiscreteTransferFunction(-k, 1, s, 0.001), DiscreteTransferFunction(k, 1, s, 0.001)), (DiscreteTransferFunction(a2*p - s, a2*s + p, s, 0.001), DiscreteTransferFunction(k, 1, s, 0.001)))))")
  expect_identical(got[["lti:test_MIMOSeries_functions:1590"]], "MIMOParallel(MIMOSeries(TransferFunctionMatrix(((DiscreteTransferFunction(a2*p - s, a2*s + p, s, 0.001),), (DiscreteTransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s, 0.001),))), TransferFunctionMatrix(((DiscreteTransferFunction(-k, 1, s, 0.001), DiscreteTransferFunction(-a2*p + s, a2*s + p, s, 0.001)), (DiscreteTransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s, 0.001), DiscreteTransferFunction(k, 1, s, 0.001))))), MIMOSeries(TransferFunctionMatrix(((DiscreteTransferFunction(-a2*p + s, a2*s + p, s, 0.001),), (DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.001),))), TransferFunctionMatrix(((DiscreteTransferFunction(k, 1, s, 0.001), DiscreteTransferFunction(a2*p - s, a2*s + p, s, 0.001)), (DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.001), DiscreteTransferFunction(-k, 1, s, 0.001))))))")
  expect_identical(got[["lti:test_MIMOSeries_functions:1593"]], "ERR ValueError: <class 'sympy.physics.control.lti.DiscreteTransferFunction'> cannot be added with <class 'sympy.physics.control.lti.MIMOSeries'>.")
  expect_identical(got[["lti:test_MIMOSeries_functions:1594"]], "ERR TypeError: unsupported operand type(s) for +: 'MIMOSeries' and 'Symbol'")
  expect_identical(got[["lti:test_MIMOSeries_functions:1595"]], "ERR TypeError: unsupported operand type(s) for -: 'MIMOSeries' and 'Add'")
  expect_identical(got[["lti:test_MIMOSeries_functions:1596"]], "ERR TypeError: unsupported operand type(s) for -: 'MIMOSeries' and 'int'")
  expect_identical(got[["lti:test_MIMOSeries_functions:1597"]], "ERR TypeError: unsupported operand type(s) for +: 'Add' and 'MIMOSeries'")
  expect_identical(got[["lti:test_MIMOSeries_functions:1601"]], "ERR TypeError: Shape of all the args is not equal.")
  expect_identical(got[["lti:test_MIMOSeries_functions:1602"]], "ERR TypeError: Shape of all the args is not equal.")
  expect_identical(got[["lti:test_MIMOSeries_functions:1603"]], "ERR TypeError: Shape of all the args is not equal.")
  expect_identical(got[["lti:test_MIMOSeries_functions:1605"]], "MIMOSeries(TransferFunctionMatrix(((DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.001),),)), TransferFunctionMatrix(((DiscreteTransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s, 0.001),), (DiscreteTransferFunction(-k, 1, s, 0.001),), (DiscreteTransferFunction(-a2*p + s, a2*s + p, s, 0.001),))), TransferFunctionMatrix(((DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.001), DiscreteTransferFunction(k, 1, s, 0.001), DiscreteTransferFunction(a2*p - s, a2*s + p, s, 0.001)), (DiscreteTransferFunction(-a2*p + s, a2*s + p, s, 0.001), DiscreteTransferFunction(-k, 1, s, 0.001), DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.001)))))")
  expect_identical(got[["lti:test_MIMOSeries_functions:1606"]], "MIMOSeries(TransferFunctionMatrix(((DiscreteTransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s, 0.001),),)), TransferFunctionMatrix(((DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.001),), (DiscreteTransferFunction(k, 1, s, 0.001),), (DiscreteTransferFunction(a2*p - s, a2*s + p, s, 0.001),))), TransferFunctionMatrix(((DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.001), DiscreteTransferFunction(k, 1, s, 0.001), DiscreteTransferFunction(a2*p - s, a2*s + p, s, 0.001)), (DiscreteTransferFunction(-a2*p + s, a2*s + p, s, 0.001), DiscreteTransferFunction(-k, 1, s, 0.001), DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.001)))))")
  expect_identical(got[["lti:test_MIMOSeries_functions:1610"]], "ERR ValueError: <class 'sympy.physics.control.lti.DiscreteTransferFunction'> cannot be multiplied with <class 'sympy.physics.control.lti.MIMOSeries'>.")
  expect_identical(got[["lti:test_MIMOSeries_functions:1611"]], "ERR TypeError: unsupported operand type(s) for *: 'MIMOSeries' and 'Symbol'")
  expect_identical(got[["lti:test_MIMOSeries_functions:1612"]], "ERR TypeError: unsupported operand type(s) for *: 'MIMOSeries' and 'Add'")
  expect_identical(got[["lti:test_MIMOSeries_functions:1613"]], "ERR TypeError: unsupported operand type(s) for *: 'MIMOSeries' and 'int'")
  expect_identical(got[["lti:test_MIMOSeries_functions:1614"]], "ERR TypeError: unsupported operand type(s) for *: 'Add' and 'TransferFunctionMatrix'")
  expect_identical(got[["lti:test_MIMOSeries_functions:1617"]], "True")
  expect_identical(got[["lti:test_MIMOSeries_functions:1638"]], "TransferFunctionMatrix(((DiscreteTransferFunction(2*(s + 1)^2*(s^2 + 2*s + 1)^3, (s + 1)^2*(s^2 + 2*s + 1)^3, s, 0.001),),))")
  expect_identical(got[["lti:test_MIMOSeries_functions:1642"]], "TransferFunctionMatrix(((DiscreteTransferFunction(2, 1, s, 0.001),),))")
  expect_identical(got[["lti:test_MIMOSeries_functions:1646"]], "True")
  expect_identical(got[["lti:test_MIMOSeries_functions:1661"]], "True")
})

test_that("lti test_Parallel_construction matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("a, x, b, c, s, g, d,    p, k, tau, zeta, wn, T, z"); a <- .py_getitem(.tmpu_, 0L); x <- .py_getitem(.tmpu_, 1L); b <- .py_getitem(.tmpu_, 2L); c <- .py_getitem(.tmpu_, 3L); s <- .py_getitem(.tmpu_, 4L); g <- .py_getitem(.tmpu_, 5L); d <- .py_getitem(.tmpu_, 6L); p <- .py_getitem(.tmpu_, 7L); k <- .py_getitem(.tmpu_, 8L); tau <- .py_getitem(.tmpu_, 9L); zeta <- .py_getitem(.tmpu_, 10L); wn <- .py_getitem(.tmpu_, 11L); T <- .py_getitem(.tmpu_, 12L); z <- .py_getitem(.tmpu_, 13L)
})
.py_try({
.tmpu_ <- .py_symbols("a0:4, b0:5, c0:4, d0:4"); a0 <- .py_getitem(.tmpu_, 0L); a1 <- .py_getitem(.tmpu_, 1L); a2 <- .py_getitem(.tmpu_, 2L); a3 <- .py_getitem(.tmpu_, 3L); b0 <- .py_getitem(.tmpu_, 4L); b1 <- .py_getitem(.tmpu_, 5L); b2 <- .py_getitem(.tmpu_, 6L); b3 <- .py_getitem(.tmpu_, 7L); b4 <- .py_getitem(.tmpu_, 8L); c0 <- .py_getitem(.tmpu_, 9L); c1 <- .py_getitem(.tmpu_, 10L); c2 <- .py_getitem(.tmpu_, 11L); c3 <- .py_getitem(.tmpu_, 12L); d0 <- .py_getitem(.tmpu_, 13L); d1 <- .py_getitem(.tmpu_, 14L); d2 <- .py_getitem(.tmpu_, 15L); d3 <- .py_getitem(.tmpu_, 16L)
})
.py_try({
TF1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
TF2 <- TransferFunction(k, 1L, s)
})
.py_try({
TF3 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
TF4 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(tau, s), 1L), s)
})
.py_try({
TF5 <- TransferFunction(k, p, s)
})
.py_try({
TF6 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(k, s), p), s)
})
p_are_floats_equals <- function(n1, n2) {
  "\n    Helper function to compare two floats.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return(.py_method(Float(n1), "is_same", Float(n2), .py_isclose))
}
p_are_complex_equals <- function(n1, n2) {
  "\n    Helper function to compare two complex numbers.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return((.py_truthy(p_are_floats_equals(re(n1), re(n2))) && .py_truthy(p_are_floats_equals(im(n1), im(n2)))))
}
.py_try({
tf <- TransferFunction(.py_sub(.py_add(.py_mul(a0, .py_pow(s, 3L)), .py_mul(a1, .py_pow(s, 2L))), .py_mul(a2, s)), .py_sub(.py_add(.py_mul(b0, .py_pow(p, 4L)), .py_mul(b1, .py_pow(p, 3L))), .py_mul(.py_mul(b2, s), p)), s)
})
.py_try({
tf2 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
tf3 <- TransferFunction(.py_sub(.py_add(.py_mul(a0, p), .py_pow(p, a1)), s), p, p)
})
.py_try({
tf4 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
inp <- Function("X_d")(s)
})
.py_try({
out <- Function("X")(s)
})
.py_try({
dtf <- DiscreteTransferFunction(.py_sub(.py_add(.py_mul(a0, .py_pow(s, 3L)), .py_mul(a1, .py_pow(s, 2L))), .py_mul(a2, s)), .py_sub(.py_add(.py_mul(b0, .py_pow(p, 4L)), .py_mul(b1, .py_pow(p, 3L))), .py_mul(.py_mul(b2, s), p)), s, .py_flt(0.5))
})
.py_try({
dtf2 <- DiscreteTransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s, .py_flt(0.5))
})
.py_try({
dtf3 <- DiscreteTransferFunction(.py_sub(.py_add(.py_mul(a0, p), .py_pow(p, a1)), s), p, p, .py_flt(0.5))
})
.py_try({
dtf4 <- DiscreteTransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s, .py_flt(0.5))
})
.py_try({
.py_rec("lti:test_Parallel_construction:1692", function() .py_attr(tf, "is_continuous"))
})
.py_try({
.py_rec("lti:test_Parallel_construction:1693", function() .py_attr(dtf, "is_continuous"))
})
.py_try({
p0 <- Parallel(tf, tf2)
})
.py_try({
.py_rec("lti:test_Parallel_construction:1697", function() .py_attr(p0, "args"))
})
.py_try({
.py_rec("lti:test_Parallel_construction:1698", function() .py_attr(p0, "var"))
})
.py_try({
p1 <- Parallel(Series(tf, .py_neg(tf2)), tf2)
})
.py_try({
.py_rec("lti:test_Parallel_construction:1701", function() .py_attr(p1, "args"))
})
.py_try({
.py_rec("lti:test_Parallel_construction:1702", function() .py_attr(p1, "var"))
})
.py_try({
tf3_ <- TransferFunction(inp, 1L, s)
})
.py_try({
tf4_ <- TransferFunction(.py_neg(out), 1L, s)
})
.py_try({
p2 <- Parallel(tf, Series(tf3_, .py_neg(tf4_)), tf2)
})
.py_try({
.py_rec("lti:test_Parallel_construction:1707", function() .py_attr(p2, "args"))
})
.py_try({
p3 <- Parallel(tf, tf2, tf4)
})
.py_try({
.py_rec("lti:test_Parallel_construction:1710", function() .py_attr(p3, "args"))
})
.py_try({
p4 <- Parallel(tf3_, tf4_)
})
.py_try({
.py_rec("lti:test_Parallel_construction:1713", function() .py_attr(p4, "args"))
})
.py_try({
.py_rec("lti:test_Parallel_construction:1714", function() .py_attr(p4, "var"))
})
.py_try({
p5 <- Parallel(tf, tf2)
})
.py_try({
.py_rec("lti:test_Parallel_construction:1717", function() p0)
})
.py_try({
.py_rec("lti:test_Parallel_construction:1718", function() p0)
})
.py_try({
p6 <- Parallel(tf2, tf4, Series(tf2, .py_neg(tf4)))
})
.py_try({
.py_rec("lti:test_Parallel_construction:1721", function() .py_attr(p6, "args"))
})
.py_try({
p7 <- Parallel(tf2, tf4, Series(tf2, .py_neg(tf)), tf4)
})
.py_try({
.py_rec("lti:test_Parallel_construction:1724", function() .py_attr(p7, "args"))
})
.py_try({
.py_rec("lti:test_Parallel_construction:1726", function() Parallel(tf, tf3))
})
.py_try({
.py_rec("lti:test_Parallel_construction:1727", function() Parallel(tf, tf2, tf3, tf4))
})
.py_try({
.py_rec("lti:test_Parallel_construction:1728", function() Parallel(.py_neg(tf3), tf4))
})
.py_try({
.py_rec("lti:test_Parallel_construction:1729", function() Parallel(2L, tf, tf4))
})
.py_try({
.py_rec("lti:test_Parallel_construction:1730", function() Parallel(.py_add(.py_pow(s, 2L), .py_mul(p, s)), tf3, tf2))
})
.py_try({
.py_rec("lti:test_Parallel_construction:1731", function() Parallel(tf3, .py_Matrix(list(1L, 2L, 3L, 4L))))
})
.py_try({
dp0 <- Parallel(dtf, dtf2)
})
.py_try({
.py_rec("lti:test_Parallel_construction:1735", function() .py_attr(dp0, "args"))
})
.py_try({
.py_rec("lti:test_Parallel_construction:1736", function() .py_attr(dp0, "var"))
})
.py_try({
.py_rec("lti:test_Parallel_construction:1737", function() .py_attr(dp0, "sampling_time"))
})
.py_try({
dp1 <- Parallel(Series(dtf, .py_neg(dtf2)), dtf2)
})
.py_try({
.py_rec("lti:test_Parallel_construction:1740", function() .py_attr(dp1, "args"))
})
.py_try({
.py_rec("lti:test_Parallel_construction:1741", function() .py_attr(dp1, "var"))
})
.py_try({
.py_rec("lti:test_Parallel_construction:1742", function() .py_attr(dp1, "sampling_time"))
})
.py_try({
dtf3_ <- DiscreteTransferFunction(inp, 1L, s, .py_flt(0.5))
})
.py_try({
dtf4_ <- DiscreteTransferFunction(.py_neg(out), 1L, s, .py_flt(0.5))
})
.py_try({
dp2 <- Parallel(dtf, Series(dtf3_, .py_neg(dtf4_)), dtf2)
})
.py_try({
.py_rec("lti:test_Parallel_construction:1748", function() .py_attr(dp2, "args"))
})
.py_try({
.py_rec("lti:test_Parallel_construction:1749", function() .py_attr(dp2, "sampling_time"))
})
.py_try({
dp3 <- Parallel(dtf, dtf2, dtf4)
})
.py_try({
.py_rec("lti:test_Parallel_construction:1752", function() .py_attr(dp3, "args"))
})
.py_try({
.py_rec("lti:test_Parallel_construction:1753", function() .py_attr(dp3, "sampling_time"))
})
.py_try({
dp4 <- Parallel(dtf3_, dtf4_)
})
.py_try({
.py_rec("lti:test_Parallel_construction:1756", function() .py_attr(dp4, "args"))
})
.py_try({
.py_rec("lti:test_Parallel_construction:1757", function() .py_attr(dp4, "var"))
})
.py_try({
.py_rec("lti:test_Parallel_construction:1758", function() .py_attr(dp4, "sampling_time"))
})
.py_try({
dp5 <- Parallel(dtf, dtf2)
})
.py_try({
.py_rec("lti:test_Parallel_construction:1761", function() dp0)
})
.py_try({
.py_rec("lti:test_Parallel_construction:1762", function() dp0)
})
.py_try({
.py_rec("lti:test_Parallel_construction:1763", function() .py_attr(dp5, "sampling_time"))
})
.py_try({
dp6 <- Parallel(dtf2, dtf4, Series(dtf2, .py_neg(dtf4)))
})
.py_try({
.py_rec("lti:test_Parallel_construction:1766", function() .py_attr(dp6, "args"))
})
.py_try({
.py_rec("lti:test_Parallel_construction:1767", function() .py_attr(dp6, "sampling_time"))
})
.py_try({
dp7 <- Parallel(dtf2, dtf4, Series(dtf2, .py_neg(dtf)), dtf4)
})
.py_try({
.py_rec("lti:test_Parallel_construction:1770", function() .py_attr(dp7, "args"))
})
.py_try({
.py_rec("lti:test_Parallel_construction:1771", function() .py_attr(dp7, "sampling_time"))
})
.py_try({
dtf5 <- DiscreteTransferFunction(s, .py_sub(s, 1L), s, 1L)
})
.py_try({
.py_rec("lti:test_Parallel_construction:1774", function() Parallel(dtf, dtf3))
})
.py_try({
.py_rec("lti:test_Parallel_construction:1775", function() Parallel(dtf, dtf2, dtf3, dtf4))
})
.py_try({
.py_rec("lti:test_Parallel_construction:1776", function() Parallel(.py_neg(dtf3), dtf4))
})
.py_try({
.py_rec("lti:test_Parallel_construction:1777", function() Parallel(dtf5, dtf))
})
.py_try({
.py_rec("lti:test_Parallel_construction:1778", function() Parallel(dtf5, tf))
})
.py_try({
.py_rec("lti:test_Parallel_construction:1779", function() Parallel(2L, dtf, dtf4))
})
.py_try({
.py_rec("lti:test_Parallel_construction:1780", function() Parallel(.py_add(.py_pow(s, 2L), .py_mul(p, s)), dtf3, dtf2))
})
.py_try({
.py_rec("lti:test_Parallel_construction:1781", function() Parallel(dtf3, .py_Matrix(list(1L, 2L, 3L, 4L))))
}))---")
  expect_identical(got[["lti:test_Parallel_construction:1692"]], "True")
  expect_identical(got[["lti:test_Parallel_construction:1693"]], "False")
  expect_identical(got[["lti:test_Parallel_construction:1697"]], "(TransferFunction(a0*s^3 + a1*s^2 - a2*s, b0*p^4 + b1*p^3 - b2*p*s, s), TransferFunction(a2*p - s, a2*s + p, s))")
  expect_identical(got[["lti:test_Parallel_construction:1698"]], "s")
  expect_identical(got[["lti:test_Parallel_construction:1701"]], "(Series(TransferFunction(a0*s^3 + a1*s^2 - a2*s, b0*p^4 + b1*p^3 - b2*p*s, s), TransferFunction(-a2*p + s, a2*s + p, s)), TransferFunction(a2*p - s, a2*s + p, s))")
  expect_identical(got[["lti:test_Parallel_construction:1702"]], "s")
  expect_identical(got[["lti:test_Parallel_construction:1707"]], "(TransferFunction(a0*s^3 + a1*s^2 - a2*s, b0*p^4 + b1*p^3 - b2*p*s, s), Series(TransferFunction(X_d(s), 1, s), TransferFunction(X(s), 1, s)), TransferFunction(a2*p - s, a2*s + p, s))")
  expect_identical(got[["lti:test_Parallel_construction:1710"]], "(TransferFunction(a0*s^3 + a1*s^2 - a2*s, b0*p^4 + b1*p^3 - b2*p*s, s), TransferFunction(a2*p - s, a2*s + p, s), TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s))")
  expect_identical(got[["lti:test_Parallel_construction:1713"]], "(TransferFunction(X_d(s), 1, s), TransferFunction(-X(s), 1, s))")
  expect_identical(got[["lti:test_Parallel_construction:1714"]], "s")
  expect_identical(got[["lti:test_Parallel_construction:1717"]], "Parallel(TransferFunction(a0*s^3 + a1*s^2 - a2*s, b0*p^4 + b1*p^3 - b2*p*s, s), TransferFunction(a2*p - s, a2*s + p, s))")
  expect_identical(got[["lti:test_Parallel_construction:1718"]], "Parallel(TransferFunction(a0*s^3 + a1*s^2 - a2*s, b0*p^4 + b1*p^3 - b2*p*s, s), TransferFunction(a2*p - s, a2*s + p, s))")
  expect_identical(got[["lti:test_Parallel_construction:1721"]], "(TransferFunction(a2*p - s, a2*s + p, s), TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s), Series(TransferFunction(a2*p - s, a2*s + p, s), TransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s)))")
  expect_identical(got[["lti:test_Parallel_construction:1724"]], "(TransferFunction(a2*p - s, a2*s + p, s), TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s), Series(TransferFunction(a2*p - s, a2*s + p, s), TransferFunction(-a0*s^3 - a1*s^2 + a2*s, b0*p^4 + b1*p^3 - b2*p*s, s)), TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s))")
  expect_identical(got[["lti:test_Parallel_construction:1726"]], "ERR ValueError: \\nAll transfer functions should use the same complex variable of the\\nLaplace transform or z-transform. 2 different values found.")
  expect_identical(got[["lti:test_Parallel_construction:1727"]], "ERR ValueError: \\nAll transfer functions should use the same complex variable of the\\nLaplace transform or z-transform. 2 different values found.")
  expect_identical(got[["lti:test_Parallel_construction:1728"]], "ERR ValueError: \\nAll transfer functions should use the same complex variable of the\\nLaplace transform or z-transform. 2 different values found.")
  expect_identical(got[["lti:test_Parallel_construction:1729"]], "ERR TypeError: All arguments must be of type <class 'sympy.physics.control.lti.SISOLinearTimeInvariant'>.")
  expect_identical(got[["lti:test_Parallel_construction:1730"]], "ERR TypeError: All arguments must be of type <class 'sympy.physics.control.lti.SISOLinearTimeInvariant'>.")
  expect_identical(got[["lti:test_Parallel_construction:1731"]], "ERR TypeError: All arguments must be of type <class 'sympy.physics.control.lti.SISOLinearTimeInvariant'>.")
  expect_identical(got[["lti:test_Parallel_construction:1735"]], "(DiscreteTransferFunction(a0*s^3 + a1*s^2 - a2*s, b0*p^4 + b1*p^3 - b2*p*s, s, 0.5), DiscreteTransferFunction(a2*p - s, a2*s + p, s, 0.5))")
  expect_identical(got[["lti:test_Parallel_construction:1736"]], "s")
  expect_identical(got[["lti:test_Parallel_construction:1737"]], "0.500000000000000")
  expect_identical(got[["lti:test_Parallel_construction:1740"]], "(Series(DiscreteTransferFunction(a0*s^3 + a1*s^2 - a2*s, b0*p^4 + b1*p^3 - b2*p*s, s, 0.5), DiscreteTransferFunction(-a2*p + s, a2*s + p, s, 0.5)), DiscreteTransferFunction(a2*p - s, a2*s + p, s, 0.5))")
  expect_identical(got[["lti:test_Parallel_construction:1741"]], "s")
  expect_identical(got[["lti:test_Parallel_construction:1742"]], "0.500000000000000")
  expect_identical(got[["lti:test_Parallel_construction:1748"]], "(DiscreteTransferFunction(a0*s^3 + a1*s^2 - a2*s, b0*p^4 + b1*p^3 - b2*p*s, s, 0.5), Series(DiscreteTransferFunction(X_d(s), 1, s, 0.5), DiscreteTransferFunction(X(s), 1, s, 0.5)), DiscreteTransferFunction(a2*p - s, a2*s + p, s, 0.5))")
  expect_identical(got[["lti:test_Parallel_construction:1749"]], "0.500000000000000")
  expect_identical(got[["lti:test_Parallel_construction:1752"]], "(DiscreteTransferFunction(a0*s^3 + a1*s^2 - a2*s, b0*p^4 + b1*p^3 - b2*p*s, s, 0.5), DiscreteTransferFunction(a2*p - s, a2*s + p, s, 0.5), DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.5))")
  expect_identical(got[["lti:test_Parallel_construction:1753"]], "0.500000000000000")
  expect_identical(got[["lti:test_Parallel_construction:1756"]], "(DiscreteTransferFunction(X_d(s), 1, s, 0.5), DiscreteTransferFunction(-X(s), 1, s, 0.5))")
  expect_identical(got[["lti:test_Parallel_construction:1757"]], "s")
  expect_identical(got[["lti:test_Parallel_construction:1758"]], "0.500000000000000")
  expect_identical(got[["lti:test_Parallel_construction:1761"]], "Parallel(DiscreteTransferFunction(a0*s^3 + a1*s^2 - a2*s, b0*p^4 + b1*p^3 - b2*p*s, s, 0.5), DiscreteTransferFunction(a2*p - s, a2*s + p, s, 0.5))")
  expect_identical(got[["lti:test_Parallel_construction:1762"]], "Parallel(DiscreteTransferFunction(a0*s^3 + a1*s^2 - a2*s, b0*p^4 + b1*p^3 - b2*p*s, s, 0.5), DiscreteTransferFunction(a2*p - s, a2*s + p, s, 0.5))")
  expect_identical(got[["lti:test_Parallel_construction:1763"]], "0.500000000000000")
  expect_identical(got[["lti:test_Parallel_construction:1766"]], "(DiscreteTransferFunction(a2*p - s, a2*s + p, s, 0.5), DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.5), Series(DiscreteTransferFunction(a2*p - s, a2*s + p, s, 0.5), DiscreteTransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s, 0.5)))")
  expect_identical(got[["lti:test_Parallel_construction:1767"]], "0.500000000000000")
  expect_identical(got[["lti:test_Parallel_construction:1770"]], "(DiscreteTransferFunction(a2*p - s, a2*s + p, s, 0.5), DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.5), Series(DiscreteTransferFunction(a2*p - s, a2*s + p, s, 0.5), DiscreteTransferFunction(-a0*s^3 - a1*s^2 + a2*s, b0*p^4 + b1*p^3 - b2*p*s, s, 0.5)), DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.5))")
  expect_identical(got[["lti:test_Parallel_construction:1771"]], "0.500000000000000")
  expect_identical(got[["lti:test_Parallel_construction:1774"]], "ERR ValueError: \\nAll transfer functions should use the same complex variable of the\\nLaplace transform or z-transform. 2 different values found.")
  expect_identical(got[["lti:test_Parallel_construction:1775"]], "ERR ValueError: \\nAll transfer functions should use the same complex variable of the\\nLaplace transform or z-transform. 2 different values found.")
  expect_identical(got[["lti:test_Parallel_construction:1776"]], "ERR ValueError: \\nAll transfer functions should use the same complex variable of the\\nLaplace transform or z-transform. 2 different values found.")
  expect_identical(got[["lti:test_Parallel_construction:1777"]], "ERR TypeError: \\n                Incompatible sampling times detected.\\n                All systems should have the same sampling time.\\n                Found sampling time: 1 and 0.500000000000000")
  expect_identical(got[["lti:test_Parallel_construction:1778"]], "ERR TypeError: \\n                Incompatible systems detected.\\n                All systems should be either continuous-time or discrete-time.\\n                Found: DiscreteTransferFunction(s, s - 1, s, 1) and TransferFunction(a0*s^3 + a1*s^2 - a2*s, b0*p^4 + b1*p^3 - b2*p*s, s)")
  expect_identical(got[["lti:test_Parallel_construction:1779"]], "ERR TypeError: All arguments must be of type <class 'sympy.physics.control.lti.SISOLinearTimeInvariant'>.")
  expect_identical(got[["lti:test_Parallel_construction:1780"]], "ERR TypeError: All arguments must be of type <class 'sympy.physics.control.lti.SISOLinearTimeInvariant'>.")
  expect_identical(got[["lti:test_Parallel_construction:1781"]], "ERR TypeError: All arguments must be of type <class 'sympy.physics.control.lti.SISOLinearTimeInvariant'>.")
})

test_that("lti test_MIMOParallel_construction matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("a, x, b, c, s, g, d,    p, k, tau, zeta, wn, T, z"); a <- .py_getitem(.tmpu_, 0L); x <- .py_getitem(.tmpu_, 1L); b <- .py_getitem(.tmpu_, 2L); c <- .py_getitem(.tmpu_, 3L); s <- .py_getitem(.tmpu_, 4L); g <- .py_getitem(.tmpu_, 5L); d <- .py_getitem(.tmpu_, 6L); p <- .py_getitem(.tmpu_, 7L); k <- .py_getitem(.tmpu_, 8L); tau <- .py_getitem(.tmpu_, 9L); zeta <- .py_getitem(.tmpu_, 10L); wn <- .py_getitem(.tmpu_, 11L); T <- .py_getitem(.tmpu_, 12L); z <- .py_getitem(.tmpu_, 13L)
})
.py_try({
.tmpu_ <- .py_symbols("a0:4, b0:5, c0:4, d0:4"); a0 <- .py_getitem(.tmpu_, 0L); a1 <- .py_getitem(.tmpu_, 1L); a2 <- .py_getitem(.tmpu_, 2L); a3 <- .py_getitem(.tmpu_, 3L); b0 <- .py_getitem(.tmpu_, 4L); b1 <- .py_getitem(.tmpu_, 5L); b2 <- .py_getitem(.tmpu_, 6L); b3 <- .py_getitem(.tmpu_, 7L); b4 <- .py_getitem(.tmpu_, 8L); c0 <- .py_getitem(.tmpu_, 9L); c1 <- .py_getitem(.tmpu_, 10L); c2 <- .py_getitem(.tmpu_, 11L); c3 <- .py_getitem(.tmpu_, 12L); d0 <- .py_getitem(.tmpu_, 13L); d1 <- .py_getitem(.tmpu_, 14L); d2 <- .py_getitem(.tmpu_, 15L); d3 <- .py_getitem(.tmpu_, 16L)
})
.py_try({
TF1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
TF2 <- TransferFunction(k, 1L, s)
})
.py_try({
TF3 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
TF4 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(tau, s), 1L), s)
})
.py_try({
TF5 <- TransferFunction(k, p, s)
})
.py_try({
TF6 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(k, s), p), s)
})
p_are_floats_equals <- function(n1, n2) {
  "\n    Helper function to compare two floats.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return(.py_method(Float(n1), "is_same", Float(n2), .py_isclose))
}
p_are_complex_equals <- function(n1, n2) {
  "\n    Helper function to compare two complex numbers.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return((.py_truthy(p_are_floats_equals(re(n1), re(n2))) && .py_truthy(p_are_floats_equals(im(n1), im(n2)))))
}
.py_try({
tfm1 <- TransferFunctionMatrix(list(list(TF1), list(TF2), list(TF3)))
})
.py_try({
tfm2 <- TransferFunctionMatrix(list(list(.py_neg(TF3)), list(TF2), list(TF1)))
})
.py_try({
tfm3 <- TransferFunctionMatrix(list(list(TF1)))
})
.py_try({
tfm4 <- TransferFunctionMatrix(list(list(TF2), list(TF1), list(TF3)))
})
.py_try({
tfm5 <- TransferFunctionMatrix(list(list(TF1, TF2), list(TF2, TF1)))
})
.py_try({
tfm6 <- TransferFunctionMatrix(list(list(TF2, TF1), list(TF1, TF2)))
})
.py_try({
tfm7 <- .py_method(TransferFunctionMatrix, "from_Matrix", .py_Matrix(list(list(.py_div(1L, p)))), p)
})
.py_try({
dTF1 <- DiscreteTransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s, .py_flt(0.5))
})
.py_try({
dTF2 <- DiscreteTransferFunction(k, 1L, s, .py_flt(0.5))
})
.py_try({
dTF3 <- DiscreteTransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s, .py_flt(0.5))
})
.py_try({
dtfm1 <- TransferFunctionMatrix(list(list(dTF1), list(dTF2), list(dTF3)))
})
.py_try({
dtfm2 <- TransferFunctionMatrix(list(list(.py_neg(dTF3)), list(dTF2), list(dTF1)))
})
.py_try({
dtfm3 <- TransferFunctionMatrix(list(list(dTF1)))
})
.py_try({
dtfm4 <- TransferFunctionMatrix(list(list(dTF2), list(dTF1), list(dTF3)))
})
.py_try({
dtfm5 <- TransferFunctionMatrix(list(list(dTF1, dTF2), list(dTF2, dTF1)))
})
.py_try({
dtfm6 <- TransferFunctionMatrix(list(list(dTF2, dTF1), list(dTF1, dTF2)))
})
.py_try({
p8 <- MIMOParallel(tfm1, tfm2)
})
.py_try({
.py_rec("lti:test_MIMOParallel_construction:1805", function() .py_attr(p8, "args"))
})
.py_try({
.py_rec("lti:test_MIMOParallel_construction:1806", function() .py_attr(p8, "var"))
})
.py_try({
.py_rec("lti:test_MIMOParallel_construction:1807", function() (.py_truthy(.py_eq(.py_attr(p8, "shape"), .py_tup(.py_attr(p8, "num_outputs"), .py_attr(p8, "num_inputs")))) && .py_truthy(.py_eq(.py_tup(.py_attr(p8, "num_outputs"), .py_attr(p8, "num_inputs")), .py_tup(3L, 1L)))))
})
.py_try({
.py_rec("lti:test_MIMOParallel_construction:1808", function() .py_attr(p8, "is_continuous"))
})
.py_try({
p9 <- MIMOParallel(MIMOSeries(tfm3, tfm1), tfm2)
})
.py_try({
.py_rec("lti:test_MIMOParallel_construction:1811", function() .py_attr(p9, "args"))
})
.py_try({
.py_rec("lti:test_MIMOParallel_construction:1812", function() .py_attr(p9, "var"))
})
.py_try({
.py_rec("lti:test_MIMOParallel_construction:1813", function() (.py_truthy(.py_eq(.py_attr(p9, "shape"), .py_tup(.py_attr(p9, "num_outputs"), .py_attr(p9, "num_inputs")))) && .py_truthy(.py_eq(.py_tup(.py_attr(p9, "num_outputs"), .py_attr(p9, "num_inputs")), .py_tup(3L, 1L)))))
})
.py_try({
p10 <- MIMOParallel(tfm1, MIMOSeries(tfm3, tfm4), tfm2)
})
.py_try({
.py_rec("lti:test_MIMOParallel_construction:1816", function() .py_attr(p10, "args"))
})
.py_try({
.py_rec("lti:test_MIMOParallel_construction:1817", function() .py_attr(p10, "var"))
})
.py_try({
.py_rec("lti:test_MIMOParallel_construction:1818", function() (.py_truthy(.py_eq(.py_attr(p10, "shape"), .py_tup(.py_attr(p10, "num_outputs"), .py_attr(p10, "num_inputs")))) && .py_truthy(.py_eq(.py_tup(.py_attr(p10, "num_outputs"), .py_attr(p10, "num_inputs")), .py_tup(3L, 1L)))))
})
.py_try({
p11 <- MIMOParallel(tfm2, tfm1, tfm4)
})
.py_try({
.py_rec("lti:test_MIMOParallel_construction:1821", function() .py_attr(p11, "args"))
})
.py_try({
.py_rec("lti:test_MIMOParallel_construction:1822", function() (.py_truthy(.py_eq(.py_attr(p11, "shape"), .py_tup(.py_attr(p11, "num_outputs"), .py_attr(p11, "num_inputs")))) && .py_truthy(.py_eq(.py_tup(.py_attr(p11, "num_outputs"), .py_attr(p11, "num_inputs")), .py_tup(3L, 1L)))))
})
.py_try({
p12 <- MIMOParallel(tfm6, tfm5)
})
.py_try({
.py_rec("lti:test_MIMOParallel_construction:1825", function() .py_attr(p12, "args"))
})
.py_try({
.py_rec("lti:test_MIMOParallel_construction:1826", function() (.py_truthy(.py_eq(.py_attr(p12, "shape"), .py_tup(.py_attr(p12, "num_outputs"), .py_attr(p12, "num_inputs")))) && .py_truthy(.py_eq(.py_tup(.py_attr(p12, "num_outputs"), .py_attr(p12, "num_inputs")), .py_tup(2L, 2L)))))
})
.py_try({
p13 <- MIMOParallel(tfm2, tfm4, MIMOSeries(.py_neg(tfm3), tfm4), .py_neg(tfm4))
})
.py_try({
.py_rec("lti:test_MIMOParallel_construction:1829", function() .py_attr(p13, "args"))
})
.py_try({
.py_rec("lti:test_MIMOParallel_construction:1830", function() (.py_truthy(.py_eq(.py_attr(p13, "shape"), .py_tup(.py_attr(p13, "num_outputs"), .py_attr(p13, "num_inputs")))) && .py_truthy(.py_eq(.py_tup(.py_attr(p13, "num_outputs"), .py_attr(p13, "num_inputs")), .py_tup(3L, 1L)))))
})
.py_try({
dp8 <- MIMOParallel(dtfm1, dtfm2)
})
.py_try({
.py_rec("lti:test_MIMOParallel_construction:1834", function() .py_attr(dp8, "args"))
})
.py_try({
.py_rec("lti:test_MIMOParallel_construction:1835", function() .py_attr(dp8, "var"))
})
.py_try({
.py_rec("lti:test_MIMOParallel_construction:1836", function() (.py_truthy(.py_eq(.py_attr(dp8, "shape"), .py_tup(.py_attr(dp8, "num_outputs"), .py_attr(dp8, "num_inputs")))) && .py_truthy(.py_eq(.py_tup(.py_attr(dp8, "num_outputs"), .py_attr(dp8, "num_inputs")), .py_tup(3L, 1L)))))
})
.py_try({
.py_rec("lti:test_MIMOParallel_construction:1837", function() .py_attr(dp8, "sampling_time"))
})
.py_try({
.py_rec("lti:test_MIMOParallel_construction:1838", function() .py_attr(dp8, "is_continuous"))
})
.py_try({
dp9 <- MIMOParallel(MIMOSeries(dtfm3, dtfm1), dtfm2)
})
.py_try({
.py_rec("lti:test_MIMOParallel_construction:1841", function() .py_attr(dp9, "args"))
})
.py_try({
.py_rec("lti:test_MIMOParallel_construction:1842", function() .py_attr(dp9, "var"))
})
.py_try({
.py_rec("lti:test_MIMOParallel_construction:1843", function() (.py_truthy(.py_eq(.py_attr(dp9, "shape"), .py_tup(.py_attr(dp9, "num_outputs"), .py_attr(dp9, "num_inputs")))) && .py_truthy(.py_eq(.py_tup(.py_attr(dp9, "num_outputs"), .py_attr(dp9, "num_inputs")), .py_tup(3L, 1L)))))
})
.py_try({
dp10 <- MIMOParallel(dtfm1, MIMOSeries(dtfm3, dtfm4), dtfm2)
})
.py_try({
.py_rec("lti:test_MIMOParallel_construction:1846", function() .py_attr(dp10, "args"))
})
.py_try({
.py_rec("lti:test_MIMOParallel_construction:1847", function() .py_attr(dp10, "var"))
})
.py_try({
.py_rec("lti:test_MIMOParallel_construction:1848", function() (.py_truthy(.py_eq(.py_attr(dp10, "shape"), .py_tup(.py_attr(dp10, "num_outputs"), .py_attr(dp10, "num_inputs")))) && .py_truthy(.py_eq(.py_tup(.py_attr(dp10, "num_outputs"), .py_attr(dp10, "num_inputs")), .py_tup(3L, 1L)))))
})
.py_try({
dp11 <- MIMOParallel(dtfm2, dtfm1, dtfm4)
})
.py_try({
.py_rec("lti:test_MIMOParallel_construction:1851", function() .py_attr(dp11, "args"))
})
.py_try({
.py_rec("lti:test_MIMOParallel_construction:1852", function() (.py_truthy(.py_eq(.py_attr(dp11, "shape"), .py_tup(.py_attr(dp11, "num_outputs"), .py_attr(dp11, "num_inputs")))) && .py_truthy(.py_eq(.py_tup(.py_attr(dp11, "num_outputs"), .py_attr(dp11, "num_inputs")), .py_tup(3L, 1L)))))
})
.py_try({
dp12 <- MIMOParallel(dtfm6, dtfm5)
})
.py_try({
.py_rec("lti:test_MIMOParallel_construction:1855", function() .py_attr(dp12, "args"))
})
.py_try({
.py_rec("lti:test_MIMOParallel_construction:1856", function() (.py_truthy(.py_eq(.py_attr(dp12, "shape"), .py_tup(.py_attr(dp12, "num_outputs"), .py_attr(dp12, "num_inputs")))) && .py_truthy(.py_eq(.py_tup(.py_attr(dp12, "num_outputs"), .py_attr(dp12, "num_inputs")), .py_tup(2L, 2L)))))
})
.py_try({
dp13 <- MIMOParallel(dtfm2, dtfm4, MIMOSeries(.py_neg(dtfm3), dtfm4), .py_neg(dtfm4))
})
.py_try({
.py_rec("lti:test_MIMOParallel_construction:1859", function() .py_attr(dp13, "args"))
})
.py_try({
.py_rec("lti:test_MIMOParallel_construction:1860", function() (.py_truthy(.py_eq(.py_attr(dp13, "shape"), .py_tup(.py_attr(dp13, "num_outputs"), .py_attr(dp13, "num_inputs")))) && .py_truthy(.py_eq(.py_tup(.py_attr(dp13, "num_outputs"), .py_attr(dp13, "num_inputs")), .py_tup(3L, 1L)))))
})
.py_try({
.py_rec("lti:test_MIMOParallel_construction:1863", function() MIMOParallel(.py_tup()))
})
.py_try({
.py_rec("lti:test_MIMOParallel_construction:1866", function() MIMOParallel(tfm1, tfm2, TF1))
})
.py_try({
.py_rec("lti:test_MIMOParallel_construction:1869", function() MIMOParallel(tfm1, tfm3, tfm4))
})
.py_try({
.py_rec("lti:test_MIMOParallel_construction:1872", function() MIMOParallel(tfm3, tfm7))
})
.py_try({
.py_rec("lti:test_MIMOParallel_construction:1875", function() MIMOParallel(2L, tfm1, tfm4))
})
.py_try({
.py_rec("lti:test_MIMOParallel_construction:1876", function() MIMOParallel(.py_add(.py_pow(s, 2L), .py_mul(p, s)), .py_neg(tfm4), tfm2))
})
.py_try({
dtfm_6 <- .py_method(TransferFunctionMatrix, "from_Matrix", .py_mul(.py_Matrix(eye(2L)), s), s, T)
})
.py_try({
.py_rec("lti:test_MIMOParallel_construction:1880", function() MIMOParallel(dtfm_6, dtfm1))
})
.py_try({
.py_rec("lti:test_MIMOParallel_construction:1881", function() MIMOParallel(dtfm1, tfm2))
}))---")
  expect_identical(got[["lti:test_MIMOParallel_construction:1805"]], "(TransferFunctionMatrix(((TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s),), (TransferFunction(k, 1, s),), (TransferFunction(a2*p - s, a2*s + p, s),))), TransferFunctionMatrix(((TransferFunction(-a2*p + s, a2*s + p, s),), (TransferFunction(k, 1, s),), (TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s),))))")
  expect_identical(got[["lti:test_MIMOParallel_construction:1806"]], "s")
  expect_identical(got[["lti:test_MIMOParallel_construction:1807"]], "True")
  expect_identical(got[["lti:test_MIMOParallel_construction:1808"]], "True")
  expect_identical(got[["lti:test_MIMOParallel_construction:1811"]], "(MIMOSeries(TransferFunctionMatrix(((TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s),),)), TransferFunctionMatrix(((TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s),), (TransferFunction(k, 1, s),), (TransferFunction(a2*p - s, a2*s + p, s),)))), TransferFunctionMatrix(((TransferFunction(-a2*p + s, a2*s + p, s),), (TransferFunction(k, 1, s),), (TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s),))))")
  expect_identical(got[["lti:test_MIMOParallel_construction:1812"]], "s")
  expect_identical(got[["lti:test_MIMOParallel_construction:1813"]], "True")
  expect_identical(got[["lti:test_MIMOParallel_construction:1816"]], "(TransferFunctionMatrix(((TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s),), (TransferFunction(k, 1, s),), (TransferFunction(a2*p - s, a2*s + p, s),))), MIMOSeries(TransferFunctionMatrix(((TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s),),)), TransferFunctionMatrix(((TransferFunction(k, 1, s),), (TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s),), (TransferFunction(a2*p - s, a2*s + p, s),)))), TransferFunctionMatrix(((TransferFunction(-a2*p + s, a2*s + p, s),), (TransferFunction(k, 1, s),), (TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s),))))")
  expect_identical(got[["lti:test_MIMOParallel_construction:1817"]], "s")
  expect_identical(got[["lti:test_MIMOParallel_construction:1818"]], "True")
  expect_identical(got[["lti:test_MIMOParallel_construction:1821"]], "(TransferFunctionMatrix(((TransferFunction(-a2*p + s, a2*s + p, s),), (TransferFunction(k, 1, s),), (TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s),))), TransferFunctionMatrix(((TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s),), (TransferFunction(k, 1, s),), (TransferFunction(a2*p - s, a2*s + p, s),))), TransferFunctionMatrix(((TransferFunction(k, 1, s),), (TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s),), (TransferFunction(a2*p - s, a2*s + p, s),))))")
  expect_identical(got[["lti:test_MIMOParallel_construction:1822"]], "True")
  expect_identical(got[["lti:test_MIMOParallel_construction:1825"]], "(TransferFunctionMatrix(((TransferFunction(k, 1, s), TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s)), (TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s), TransferFunction(k, 1, s)))), TransferFunctionMatrix(((TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s), TransferFunction(k, 1, s)), (TransferFunction(k, 1, s), TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s)))))")
  expect_identical(got[["lti:test_MIMOParallel_construction:1826"]], "True")
  expect_identical(got[["lti:test_MIMOParallel_construction:1829"]], "(TransferFunctionMatrix(((TransferFunction(-a2*p + s, a2*s + p, s),), (TransferFunction(k, 1, s),), (TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s),))), TransferFunctionMatrix(((TransferFunction(k, 1, s),), (TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s),), (TransferFunction(a2*p - s, a2*s + p, s),))), MIMOSeries(TransferFunctionMatrix(((TransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s),),)), TransferFunctionMatrix(((TransferFunction(k, 1, s),), (TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s),), (TransferFunction(a2*p - s, a2*s + p, s),)))), TransferFunctionMatrix(((TransferFunction(-k, 1, s),), (TransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s),), (TransferFunction(-a2*p + s, a2*s + p, s),))))")
  expect_identical(got[["lti:test_MIMOParallel_construction:1830"]], "True")
  expect_identical(got[["lti:test_MIMOParallel_construction:1834"]], "(TransferFunctionMatrix(((DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.5),), (DiscreteTransferFunction(k, 1, s, 0.5),), (DiscreteTransferFunction(a2*p - s, a2*s + p, s, 0.5),))), TransferFunctionMatrix(((DiscreteTransferFunction(-a2*p + s, a2*s + p, s, 0.5),), (DiscreteTransferFunction(k, 1, s, 0.5),), (DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.5),))))")
  expect_identical(got[["lti:test_MIMOParallel_construction:1835"]], "s")
  expect_identical(got[["lti:test_MIMOParallel_construction:1836"]], "True")
  expect_identical(got[["lti:test_MIMOParallel_construction:1837"]], "0.500000000000000")
  expect_identical(got[["lti:test_MIMOParallel_construction:1838"]], "False")
  expect_identical(got[["lti:test_MIMOParallel_construction:1841"]], "(MIMOSeries(TransferFunctionMatrix(((DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.5),),)), TransferFunctionMatrix(((DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.5),), (DiscreteTransferFunction(k, 1, s, 0.5),), (DiscreteTransferFunction(a2*p - s, a2*s + p, s, 0.5),)))), TransferFunctionMatrix(((DiscreteTransferFunction(-a2*p + s, a2*s + p, s, 0.5),), (DiscreteTransferFunction(k, 1, s, 0.5),), (DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.5),))))")
  expect_identical(got[["lti:test_MIMOParallel_construction:1842"]], "s")
  expect_identical(got[["lti:test_MIMOParallel_construction:1843"]], "True")
  expect_identical(got[["lti:test_MIMOParallel_construction:1846"]], "(TransferFunctionMatrix(((DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.5),), (DiscreteTransferFunction(k, 1, s, 0.5),), (DiscreteTransferFunction(a2*p - s, a2*s + p, s, 0.5),))), MIMOSeries(TransferFunctionMatrix(((DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.5),),)), TransferFunctionMatrix(((DiscreteTransferFunction(k, 1, s, 0.5),), (DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.5),), (DiscreteTransferFunction(a2*p - s, a2*s + p, s, 0.5),)))), TransferFunctionMatrix(((DiscreteTransferFunction(-a2*p + s, a2*s + p, s, 0.5),), (DiscreteTransferFunction(k, 1, s, 0.5),), (DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.5),))))")
  expect_identical(got[["lti:test_MIMOParallel_construction:1847"]], "s")
  expect_identical(got[["lti:test_MIMOParallel_construction:1848"]], "True")
  expect_identical(got[["lti:test_MIMOParallel_construction:1851"]], "(TransferFunctionMatrix(((DiscreteTransferFunction(-a2*p + s, a2*s + p, s, 0.5),), (DiscreteTransferFunction(k, 1, s, 0.5),), (DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.5),))), TransferFunctionMatrix(((DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.5),), (DiscreteTransferFunction(k, 1, s, 0.5),), (DiscreteTransferFunction(a2*p - s, a2*s + p, s, 0.5),))), TransferFunctionMatrix(((DiscreteTransferFunction(k, 1, s, 0.5),), (DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.5),), (DiscreteTransferFunction(a2*p - s, a2*s + p, s, 0.5),))))")
  expect_identical(got[["lti:test_MIMOParallel_construction:1852"]], "True")
  expect_identical(got[["lti:test_MIMOParallel_construction:1855"]], "(TransferFunctionMatrix(((DiscreteTransferFunction(k, 1, s, 0.5), DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.5)), (DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.5), DiscreteTransferFunction(k, 1, s, 0.5)))), TransferFunctionMatrix(((DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.5), DiscreteTransferFunction(k, 1, s, 0.5)), (DiscreteTransferFunction(k, 1, s, 0.5), DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.5)))))")
  expect_identical(got[["lti:test_MIMOParallel_construction:1856"]], "True")
  expect_identical(got[["lti:test_MIMOParallel_construction:1859"]], "(TransferFunctionMatrix(((DiscreteTransferFunction(-a2*p + s, a2*s + p, s, 0.5),), (DiscreteTransferFunction(k, 1, s, 0.5),), (DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.5),))), TransferFunctionMatrix(((DiscreteTransferFunction(k, 1, s, 0.5),), (DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.5),), (DiscreteTransferFunction(a2*p - s, a2*s + p, s, 0.5),))), MIMOSeries(TransferFunctionMatrix(((DiscreteTransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s, 0.5),),)), TransferFunctionMatrix(((DiscreteTransferFunction(k, 1, s, 0.5),), (DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.5),), (DiscreteTransferFunction(a2*p - s, a2*s + p, s, 0.5),)))), TransferFunctionMatrix(((DiscreteTransferFunction(-k, 1, s, 0.5),), (DiscreteTransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s, 0.5),), (DiscreteTransferFunction(-a2*p + s, a2*s + p, s, 0.5),))))")
  expect_identical(got[["lti:test_MIMOParallel_construction:1860"]], "True")
  expect_identical(got[["lti:test_MIMOParallel_construction:1863"]], "ERR TypeError: All arguments must be of type <class 'sympy.physics.control.lti.MIMOLinearTimeInvariant'>.")
  expect_identical(got[["lti:test_MIMOParallel_construction:1866"]], "ERR TypeError: All arguments must be of type <class 'sympy.physics.control.lti.MIMOLinearTimeInvariant'>.")
  expect_identical(got[["lti:test_MIMOParallel_construction:1869"]], "ERR TypeError: Shape of all the args is not equal.")
  expect_identical(got[["lti:test_MIMOParallel_construction:1872"]], "ERR ValueError: \\nAll transfer functions should use the same complex variable of the\\nLaplace transform or z-transform. 2 different values found.")
  expect_identical(got[["lti:test_MIMOParallel_construction:1875"]], "ERR TypeError: All arguments must be of type <class 'sympy.physics.control.lti.MIMOLinearTimeInvariant'>.")
  expect_identical(got[["lti:test_MIMOParallel_construction:1876"]], "ERR TypeError: All arguments must be of type <class 'sympy.physics.control.lti.MIMOLinearTimeInvariant'>.")
  expect_identical(got[["lti:test_MIMOParallel_construction:1880"]], "ERR TypeError: Shape of all the args is not equal.")
  expect_identical(got[["lti:test_MIMOParallel_construction:1881"]], "ERR TypeError: \\n                Incompatible systems detected.\\n                All systems should be either continuous-time or discrete-time.\\n                Found: TransferFunctionMatrix(((DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.5),), (DiscreteTransferFunction(k, 1, s, 0.5),), (DiscreteTransferFunction(a2*p - s, a2*s + p, s, 0.5),))) and TransferFunctionMatrix(((TransferFunction(-a2*p + s, a2*s + p, s),), (TransferFunction(k, 1, s),), (TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s),)))")
})

test_that("lti test_Parallel_functions matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("a, x, b, c, s, g, d,    p, k, tau, zeta, wn, T, z"); a <- .py_getitem(.tmpu_, 0L); x <- .py_getitem(.tmpu_, 1L); b <- .py_getitem(.tmpu_, 2L); c <- .py_getitem(.tmpu_, 3L); s <- .py_getitem(.tmpu_, 4L); g <- .py_getitem(.tmpu_, 5L); d <- .py_getitem(.tmpu_, 6L); p <- .py_getitem(.tmpu_, 7L); k <- .py_getitem(.tmpu_, 8L); tau <- .py_getitem(.tmpu_, 9L); zeta <- .py_getitem(.tmpu_, 10L); wn <- .py_getitem(.tmpu_, 11L); T <- .py_getitem(.tmpu_, 12L); z <- .py_getitem(.tmpu_, 13L)
})
.py_try({
.tmpu_ <- .py_symbols("a0:4, b0:5, c0:4, d0:4"); a0 <- .py_getitem(.tmpu_, 0L); a1 <- .py_getitem(.tmpu_, 1L); a2 <- .py_getitem(.tmpu_, 2L); a3 <- .py_getitem(.tmpu_, 3L); b0 <- .py_getitem(.tmpu_, 4L); b1 <- .py_getitem(.tmpu_, 5L); b2 <- .py_getitem(.tmpu_, 6L); b3 <- .py_getitem(.tmpu_, 7L); b4 <- .py_getitem(.tmpu_, 8L); c0 <- .py_getitem(.tmpu_, 9L); c1 <- .py_getitem(.tmpu_, 10L); c2 <- .py_getitem(.tmpu_, 11L); c3 <- .py_getitem(.tmpu_, 12L); d0 <- .py_getitem(.tmpu_, 13L); d1 <- .py_getitem(.tmpu_, 14L); d2 <- .py_getitem(.tmpu_, 15L); d3 <- .py_getitem(.tmpu_, 16L)
})
.py_try({
TF1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
TF2 <- TransferFunction(k, 1L, s)
})
.py_try({
TF3 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
TF4 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(tau, s), 1L), s)
})
.py_try({
TF5 <- TransferFunction(k, p, s)
})
.py_try({
TF6 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(k, s), p), s)
})
p_are_floats_equals <- function(n1, n2) {
  "\n    Helper function to compare two floats.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return(.py_method(Float(n1), "is_same", Float(n2), .py_isclose))
}
p_are_complex_equals <- function(n1, n2) {
  "\n    Helper function to compare two complex numbers.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return((.py_truthy(p_are_floats_equals(re(n1), re(n2))) && .py_truthy(p_are_floats_equals(im(n1), im(n2)))))
}
.py_try({
tf1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
tf2 <- TransferFunction(k, 1L, s)
})
.py_try({
tf3 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
tf4 <- TransferFunction(.py_sub(.py_add(.py_mul(a0, p), .py_pow(p, a1)), s), p, p)
})
.py_try({
tf5 <- TransferFunction(.py_sub(.py_add(.py_mul(a1, .py_pow(s, 2L)), .py_mul(a2, s)), a0), .py_add(s, a0), s)
})
.py_try({
dtf1 <- DiscreteTransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s, .py_flt(0.01))
})
.py_try({
dtf2 <- DiscreteTransferFunction(k, 1L, s, .py_flt(0.01))
})
.py_try({
dtf3 <- DiscreteTransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s, .py_flt(0.01))
})
.py_try({
dtf4 <- DiscreteTransferFunction(.py_sub(.py_add(.py_mul(a0, p), .py_pow(p, a1)), s), p, p, .py_flt(0.01))
})
.py_try({
dtf5 <- DiscreteTransferFunction(.py_sub(.py_add(.py_mul(a1, .py_pow(s, 2L)), .py_mul(a2, s)), a0), .py_add(s, a0), s, .py_flt(0.01))
})
.py_try({
.py_rec("lti:test_Parallel_functions:1897", function() .py_add(.py_add(tf1, tf2), tf3))
})
.py_try({
.py_rec("lti:test_Parallel_functions:1898", function() .py_add(.py_add(.py_add(tf1, tf2), tf3), tf5))
})
.py_try({
.py_rec("lti:test_Parallel_functions:1899", function() .py_sub(.py_sub(.py_add(tf1, tf2), tf3), tf5))
})
.py_try({
.py_rec("lti:test_Parallel_functions:1900", function() .py_add(tf1, .py_mul(tf2, tf3)))
})
.py_try({
.py_rec("lti:test_Parallel_functions:1901", function() .py_sub(tf1, .py_mul(tf2, tf3)))
})
.py_try({
.py_rec("lti:test_Parallel_functions:1902", function() .py_sub(.py_neg(tf1), tf2))
})
.py_try({
.py_rec("lti:test_Parallel_functions:1903", function() .py_neg(.py_add(tf1, tf2)))
})
.py_try({
.py_rec("lti:test_Parallel_functions:1904", function() .py_mul(.py_add(tf2, tf3), tf1))
})
.py_try({
.py_rec("lti:test_Parallel_functions:1905", function() .py_mul(.py_add(tf1, tf2), .py_mul(tf3, tf5)))
})
.py_try({
.py_rec("lti:test_Parallel_functions:1906", function() .py_mul(.py_neg(.py_add(tf2, tf3)), .py_neg(tf5)))
})
.py_try({
.py_rec("lti:test_Parallel_functions:1907", function() .py_add(.py_add(.py_add(tf2, tf3), .py_mul(tf2, tf1)), tf5))
})
.py_try({
.py_rec("lti:test_Parallel_functions:1908", function() .py_sub(.py_add(.py_add(tf2, tf3), .py_mul(tf2, tf1)), tf3))
})
.py_try({
.py_rec("lti:test_Parallel_functions:1909", function() .py_mul(.py_add(.py_add(tf1, tf2), tf5), .py_add(tf3, tf5)))
})
.py_try({
.py_rec("lti:test_Parallel_functions:1910", function() .py_add(.py_add(tf1, tf2), tf4))
})
.py_try({
.py_rec("lti:test_Parallel_functions:1911", function() .py_sub(tf1, .py_mul(tf2, tf4)))
})
.py_try({
.py_rec("lti:test_Parallel_functions:1912", function() .py_add(tf3, .py_Matrix(list(1L, 2L, 3L))))
})
.py_try({
.py_rec("lti:test_Parallel_functions:1915", function() (.py_truthy(.py_eq(Parallel(tf1, tf2, evaluate = TRUE), .py_method(Parallel(tf1, tf2), "doit"))) && .py_truthy(.py_eq(.py_method(Parallel(tf1, tf2), "doit"), TransferFunction(.py_add(.py_mul(k, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L))), 1L), .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L)), s)))))
})
.py_try({
.py_rec("lti:test_Parallel_functions:1917", function() (.py_truthy(.py_eq(Parallel(tf1, tf2, Series(.py_neg(tf1), tf3), evaluate = TRUE), .py_method(Parallel(tf1, tf2, Series(.py_neg(tf1), tf3)), "doit"))) && .py_truthy(.py_eq(.py_method(Parallel(tf1, tf2, Series(.py_neg(tf1), tf3)), "doit"), TransferFunction(.py_add(.py_add(.py_mul(.py_mul(k, .py_add(.py_mul(a2, s), p)), .py_pow(.py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L)), 2L)), .py_mul(.py_add(.py_mul(.py_neg(a2), p), s), .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L)))), .py_mul(.py_add(.py_mul(a2, s), p), .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L)))), .py_mul(.py_add(.py_mul(a2, s), p), .py_pow(.py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L)), 2L)), s)))))
})
.py_try({
.py_rec("lti:test_Parallel_functions:1921", function() (.py_truthy(.py_eq(Parallel(tf2, tf1, .py_neg(tf3), evaluate = TRUE), .py_method(Parallel(tf2, tf1, .py_neg(tf3)), "doit"))) && .py_truthy(.py_eq(.py_method(Parallel(tf2, tf1, .py_neg(tf3)), "doit"), TransferFunction(.py_add(.py_add(.py_add(.py_mul(a2, s), .py_mul(.py_mul(k, .py_add(.py_mul(a2, s), p)), .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L)))), p), .py_mul(.py_add(.py_mul(.py_neg(a2), p), s), .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L)))), .py_mul(.py_add(.py_mul(a2, s), p), .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L))), s)))))
})
.py_try({
.py_rec("lti:test_Parallel_functions:1924", function() Parallel(tf1, .py_neg(tf2), evaluate = FALSE))
})
.py_try({
.py_rec("lti:test_Parallel_functions:1926", function() .py_method(Parallel(Series(tf1, tf2), Series(tf2, tf3)), "doit"))
})
.py_try({
.py_rec("lti:test_Parallel_functions:1928", function() .py_method(Parallel(.py_neg(tf1), .py_neg(tf2), .py_neg(tf3)), "doit"))
})
.py_try({
.py_rec("lti:test_Parallel_functions:1931", function() .py_neg(.py_method(Parallel(tf1, tf2, tf3), "doit")))
})
.py_try({
.py_rec("lti:test_Parallel_functions:1934", function() .py_method(Parallel(tf2, tf3, Series(tf2, .py_neg(tf1)), tf3), "doit"))
})
.py_try({
.py_rec("lti:test_Parallel_functions:1938", function() .py_method(Parallel(tf1, tf2), "rewrite", TransferFunction))
})
.py_try({
.py_rec("lti:test_Parallel_functions:1940", function() .py_method(Parallel(tf2, tf1, .py_neg(tf3)), "rewrite", TransferFunction))
})
.py_try({
.py_rec("lti:test_Parallel_functions:1944", function() (.py_truthy(.py_eq(Parallel(tf1, Parallel(tf2, tf3)), Parallel(tf1, tf2, tf3))) && .py_truthy(.py_eq(Parallel(tf1, tf2, tf3), Parallel(Parallel(tf1, tf2), tf3)))))
})
.py_try({
P1 <- Parallel(Series(tf1, tf2), Series(tf2, tf3))
})
.py_try({
.py_rec("lti:test_Parallel_functions:1947", function() .py_attr(P1, "is_proper"))
})
.py_try({
.py_rec("lti:test_Parallel_functions:1948", function() .py_attr(P1, "is_strictly_proper"))
})
.py_try({
.py_rec("lti:test_Parallel_functions:1949", function() .py_attr(P1, "is_biproper"))
})
.py_try({
P2 <- Parallel(tf1, .py_neg(tf2), .py_neg(tf3))
})
.py_try({
.py_rec("lti:test_Parallel_functions:1952", function() .py_attr(P2, "is_proper"))
})
.py_try({
.py_rec("lti:test_Parallel_functions:1953", function() .py_attr(P2, "is_strictly_proper"))
})
.py_try({
.py_rec("lti:test_Parallel_functions:1954", function() .py_attr(P2, "is_biproper"))
})
.py_try({
P3 <- Parallel(tf1, .py_neg(tf2), Series(tf1, tf3))
})
.py_try({
.py_rec("lti:test_Parallel_functions:1957", function() .py_attr(P3, "is_proper"))
})
.py_try({
.py_rec("lti:test_Parallel_functions:1958", function() .py_attr(P3, "is_strictly_proper"))
})
.py_try({
.py_rec("lti:test_Parallel_functions:1959", function() .py_attr(P3, "is_biproper"))
})
.py_try({
.py_rec("lti:test_Parallel_functions:1962", function() .py_add(.py_add(dtf1, dtf2), dtf3))
})
.py_try({
.py_rec("lti:test_Parallel_functions:1963", function() .py_add(.py_add(.py_add(dtf1, dtf2), dtf3), dtf5))
})
.py_try({
.py_rec("lti:test_Parallel_functions:1964", function() .py_sub(.py_sub(.py_add(dtf1, dtf2), dtf3), dtf5))
})
.py_try({
.py_rec("lti:test_Parallel_functions:1965", function() .py_add(dtf1, .py_mul(dtf2, dtf3)))
})
.py_try({
.py_rec("lti:test_Parallel_functions:1966", function() .py_sub(dtf1, .py_mul(dtf2, dtf3)))
})
.py_try({
.py_rec("lti:test_Parallel_functions:1967", function() .py_sub(.py_neg(dtf1), dtf2))
})
.py_try({
.py_rec("lti:test_Parallel_functions:1968", function() .py_neg(.py_add(dtf1, dtf2)))
})
.py_try({
.py_rec("lti:test_Parallel_functions:1970", function() .py_mul(.py_add(dtf2, dtf3), dtf1))
})
.py_try({
.py_rec("lti:test_Parallel_functions:1971", function() .py_mul(.py_add(dtf1, dtf2), .py_mul(dtf3, dtf5)))
})
.py_try({
.py_rec("lti:test_Parallel_functions:1972", function() .py_mul(.py_neg(.py_add(dtf2, dtf3)), .py_neg(dtf5)))
})
.py_try({
.py_rec("lti:test_Parallel_functions:1974", function() .py_add(.py_add(.py_add(dtf2, dtf3), .py_mul(dtf2, dtf1)), dtf5))
})
.py_try({
.py_rec("lti:test_Parallel_functions:1976", function() .py_sub(.py_add(.py_add(dtf2, dtf3), .py_mul(dtf2, dtf1)), dtf3))
})
.py_try({
.py_rec("lti:test_Parallel_functions:1978", function() .py_mul(.py_add(.py_add(dtf1, dtf2), dtf5), .py_add(dtf3, dtf5)))
})
.py_try({
.py_rec("lti:test_Parallel_functions:1981", function() .py_add(.py_add(dtf1, dtf2), dtf4))
})
.py_try({
.py_rec("lti:test_Parallel_functions:1982", function() .py_sub(dtf1, .py_mul(dtf2, dtf4)))
})
.py_try({
.py_rec("lti:test_Parallel_functions:1983", function() .py_add(dtf3, .py_Matrix(list(1L, 2L, 3L))))
})
.py_try({
.py_rec("lti:test_Parallel_functions:1986", function() (.py_truthy(.py_eq(Parallel(dtf1, dtf2, evaluate = TRUE), .py_method(Parallel(dtf1, dtf2), "doit"))) && .py_truthy(.py_eq(.py_method(Parallel(dtf1, dtf2), "doit"), DiscreteTransferFunction(.py_add(.py_mul(k, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L))), 1L), .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L)), s, .py_flt(0.01))))))
})
.py_try({
.py_rec("lti:test_Parallel_functions:1990", function() (.py_truthy(.py_eq(Parallel(dtf1, dtf2, Series(.py_neg(dtf1), dtf3), evaluate = TRUE), .py_method(Parallel(dtf1, dtf2, Series(.py_neg(dtf1), dtf3)), "doit"))) && .py_truthy(.py_eq(.py_method(Parallel(dtf1, dtf2, Series(.py_neg(dtf1), dtf3)), "doit"), DiscreteTransferFunction(.py_add(.py_add(.py_mul(.py_mul(k, .py_add(.py_mul(a2, s), p)), .py_pow(.py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L)), 2L)), .py_mul(.py_add(.py_mul(.py_neg(a2), p), s), .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L)))), .py_mul(.py_add(.py_mul(a2, s), p), .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L)))), .py_mul(.py_add(.py_mul(a2, s), p), .py_pow(.py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L)), 2L)), s, .py_flt(0.01))))))
})
.py_try({
.py_rec("lti:test_Parallel_functions:1996", function() (.py_truthy(.py_eq(Parallel(dtf2, dtf1, .py_neg(dtf3), evaluate = TRUE), .py_method(Parallel(dtf2, dtf1, .py_neg(dtf3)), "doit"))) && .py_truthy(.py_eq(.py_method(Parallel(dtf2, dtf1, .py_neg(dtf3)), "doit"), DiscreteTransferFunction(.py_add(.py_add(.py_add(.py_mul(a2, s), .py_mul(.py_mul(k, .py_add(.py_mul(a2, s), p)), .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L)))), p), .py_mul(.py_add(.py_mul(.py_neg(a2), p), s), .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L)))), .py_mul(.py_add(.py_mul(a2, s), p), .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L))), s, .py_flt(0.01))))))
})
.py_try({
.py_rec("lti:test_Parallel_functions:2001", function() Parallel(dtf1, .py_neg(dtf2), evaluate = FALSE))
})
.py_try({
.py_rec("lti:test_Parallel_functions:2004", function() .py_method(Parallel(Series(dtf1, dtf2), Series(dtf2, dtf3)), "doit"))
})
.py_try({
.py_rec("lti:test_Parallel_functions:2008", function() .py_method(Parallel(.py_neg(dtf1), .py_neg(dtf2), .py_neg(dtf3)), "doit"))
})
.py_try({
.py_rec("lti:test_Parallel_functions:2012", function() .py_neg(.py_method(Parallel(dtf1, dtf2, dtf3), "doit")))
})
.py_try({
.py_rec("lti:test_Parallel_functions:2016", function() .py_method(Parallel(dtf2, dtf3, Series(dtf2, .py_neg(dtf1)), dtf3), "doit"))
})
.py_try({
.py_rec("lti:test_Parallel_functions:2021", function() .py_method(Parallel(dtf1, dtf2), "rewrite", DiscreteTransferFunction))
})
.py_try({
.py_rec("lti:test_Parallel_functions:2024", function() .py_method(Parallel(dtf2, dtf1, .py_neg(dtf3)), "rewrite", DiscreteTransferFunction))
})
.py_try({
.py_rec("lti:test_Parallel_functions:2029", function() (.py_truthy(.py_eq(Parallel(dtf1, Parallel(dtf2, dtf3)), Parallel(dtf1, dtf2, dtf3))) && .py_truthy(.py_eq(Parallel(dtf1, dtf2, dtf3), Parallel(Parallel(dtf1, dtf2), dtf3)))))
})
.py_try({
P1 <- Parallel(Series(dtf1, dtf2), Series(dtf2, dtf3))
})
.py_try({
.py_rec("lti:test_Parallel_functions:2033", function() .py_attr(P1, "is_proper"))
})
.py_try({
.py_rec("lti:test_Parallel_functions:2034", function() .py_attr(P1, "is_strictly_proper"))
})
.py_try({
.py_rec("lti:test_Parallel_functions:2035", function() .py_attr(P1, "is_biproper"))
})
.py_try({
P2 <- Parallel(dtf1, .py_neg(dtf2), .py_neg(dtf3))
})
.py_try({
.py_rec("lti:test_Parallel_functions:2038", function() .py_attr(P2, "is_proper"))
})
.py_try({
.py_rec("lti:test_Parallel_functions:2039", function() .py_attr(P2, "is_strictly_proper"))
})
.py_try({
.py_rec("lti:test_Parallel_functions:2040", function() .py_attr(P2, "is_biproper"))
})
.py_try({
P3 <- Parallel(dtf1, .py_neg(dtf2), Series(dtf1, dtf3))
})
.py_try({
.py_rec("lti:test_Parallel_functions:2043", function() .py_attr(P3, "is_proper"))
})
.py_try({
.py_rec("lti:test_Parallel_functions:2044", function() .py_attr(P3, "is_strictly_proper"))
})
.py_try({
.py_rec("lti:test_Parallel_functions:2045", function() .py_attr(P3, "is_biproper"))
})
.py_try({
.py_rec("lti:test_Parallel_functions:2048", function() .py_method(Parallel(dtf1, dtf2), "rewrite", TransferFunction))
})
.py_try({
.py_rec("lti:test_Parallel_functions:2050", function() .py_method(Parallel(tf1, tf2), "rewrite", DiscreteTransferFunction))
}))---")
  expect_identical(got[["lti:test_Parallel_functions:1897"]], "Parallel(TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s), TransferFunction(k, 1, s), TransferFunction(a2*p - s, a2*s + p, s))")
  expect_identical(got[["lti:test_Parallel_functions:1898"]], "Parallel(TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s), TransferFunction(k, 1, s), TransferFunction(a2*p - s, a2*s + p, s), TransferFunction(-a0 + a1*s^2 + a2*s, a0 + s, s))")
  expect_identical(got[["lti:test_Parallel_functions:1899"]], "Parallel(TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s), TransferFunction(k, 1, s), TransferFunction(-a2*p + s, a2*s + p, s), TransferFunction(a0 - a1*s^2 - a2*s, a0 + s, s))")
  expect_identical(got[["lti:test_Parallel_functions:1900"]], "Parallel(TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s), Series(TransferFunction(k, 1, s), TransferFunction(a2*p - s, a2*s + p, s)))")
  expect_identical(got[["lti:test_Parallel_functions:1901"]], "Parallel(TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s), Series(TransferFunction(-1, 1, s), TransferFunction(k, 1, s), TransferFunction(a2*p - s, a2*s + p, s)))")
  expect_identical(got[["lti:test_Parallel_functions:1902"]], "Parallel(TransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s), TransferFunction(-k, 1, s))")
  expect_identical(got[["lti:test_Parallel_functions:1903"]], "Series(TransferFunction(-1, 1, s), Parallel(TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s), TransferFunction(k, 1, s)))")
  expect_identical(got[["lti:test_Parallel_functions:1904"]], "Series(Parallel(TransferFunction(k, 1, s), TransferFunction(a2*p - s, a2*s + p, s)), TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s))")
  expect_identical(got[["lti:test_Parallel_functions:1905"]], "Series(Parallel(TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s), TransferFunction(k, 1, s)), TransferFunction(a2*p - s, a2*s + p, s), TransferFunction(-a0 + a1*s^2 + a2*s, a0 + s, s))")
  expect_identical(got[["lti:test_Parallel_functions:1906"]], "Series(TransferFunction(-1, 1, s), Parallel(TransferFunction(k, 1, s), TransferFunction(a2*p - s, a2*s + p, s)), TransferFunction(a0 - a1*s^2 - a2*s, a0 + s, s))")
  expect_identical(got[["lti:test_Parallel_functions:1907"]], "Parallel(TransferFunction(k, 1, s), TransferFunction(a2*p - s, a2*s + p, s), Series(TransferFunction(k, 1, s), TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s)), TransferFunction(-a0 + a1*s^2 + a2*s, a0 + s, s))")
  expect_identical(got[["lti:test_Parallel_functions:1908"]], "Parallel(TransferFunction(k, 1, s), TransferFunction(a2*p - s, a2*s + p, s), Series(TransferFunction(k, 1, s), TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s)), TransferFunction(-a2*p + s, a2*s + p, s))")
  expect_identical(got[["lti:test_Parallel_functions:1909"]], "Series(Parallel(TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s), TransferFunction(k, 1, s), TransferFunction(-a0 + a1*s^2 + a2*s, a0 + s, s)), Parallel(TransferFunction(a2*p - s, a2*s + p, s), TransferFunction(-a0 + a1*s^2 + a2*s, a0 + s, s)))")
  expect_identical(got[["lti:test_Parallel_functions:1910"]], "ERR ValueError: \\nAll transfer functions should use the same complex variable of the\\nLaplace transform or z-transform. 2 different values found.")
  expect_identical(got[["lti:test_Parallel_functions:1911"]], "ERR ValueError: \\nAll the transfer functions should use the same complex variable of the\\nLaplace domain or z-domain.")
  expect_identical(got[["lti:test_Parallel_functions:1912"]], "ERR ValueError: <class 'sympy.physics.control.lti.TransferFunction'> cannot be added with <class 'sympy.matrices.dense.MutableDenseMatrix'>.")
  expect_identical(got[["lti:test_Parallel_functions:1915"]], "True")
  expect_identical(got[["lti:test_Parallel_functions:1917"]], "True")
  expect_identical(got[["lti:test_Parallel_functions:1921"]], "True")
  expect_identical(got[["lti:test_Parallel_functions:1924"]], "Parallel(TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s), TransferFunction(-k, 1, s))")
  expect_identical(got[["lti:test_Parallel_functions:1926"]], "TransferFunction(k*(a2*p - s)*(s^2 + 2*s*wn*zeta + wn^2) + k*(a2*s + p), (a2*s + p)*(s^2 + 2*s*wn*zeta + wn^2), s)")
  expect_identical(got[["lti:test_Parallel_functions:1928"]], "TransferFunction(-a2*s - k*(a2*s + p)*(s^2 + 2*s*wn*zeta + wn^2) - p + (-a2*p + s)*(s^2 + 2*s*wn*zeta + wn^2), (a2*s + p)*(s^2 + 2*s*wn*zeta + wn^2), s)")
  expect_identical(got[["lti:test_Parallel_functions:1931"]], "TransferFunction(-a2*s - k*(a2*s + p)*(s^2 + 2*s*wn*zeta + wn^2) - p - (a2*p - s)*(s^2 + 2*s*wn*zeta + wn^2), (a2*s + p)*(s^2 + 2*s*wn*zeta + wn^2), s)")
  expect_identical(got[["lti:test_Parallel_functions:1934"]], "TransferFunction(k*(a2*s + p)*(s^2 + 2*s*wn*zeta + wn^2) - k*(a2*s + p) + (2*a2*p - 2*s)*(s^2 + 2*s*wn*zeta + wn^2), (a2*s + p)*(s^2 + 2*s*wn*zeta + wn^2), s)")
  expect_identical(got[["lti:test_Parallel_functions:1938"]], "TransferFunction(k*(s^2 + 2*s*wn*zeta + wn^2) + 1, s^2 + 2*s*wn*zeta + wn^2, s)")
  expect_identical(got[["lti:test_Parallel_functions:1940"]], "TransferFunction(a2*s + k*(a2*s + p)*(s^2 + 2*s*wn*zeta + wn^2) + p + (-a2*p + s)*(s^2 + 2*s*wn*zeta + wn^2), (a2*s + p)*(s^2 + 2*s*wn*zeta + wn^2), s)")
  expect_identical(got[["lti:test_Parallel_functions:1944"]], "True")
  expect_identical(got[["lti:test_Parallel_functions:1947"]], "True")
  expect_identical(got[["lti:test_Parallel_functions:1948"]], "False")
  expect_identical(got[["lti:test_Parallel_functions:1949"]], "True")
  expect_identical(got[["lti:test_Parallel_functions:1952"]], "True")
  expect_identical(got[["lti:test_Parallel_functions:1953"]], "False")
  expect_identical(got[["lti:test_Parallel_functions:1954"]], "True")
  expect_identical(got[["lti:test_Parallel_functions:1957"]], "True")
  expect_identical(got[["lti:test_Parallel_functions:1958"]], "False")
  expect_identical(got[["lti:test_Parallel_functions:1959"]], "True")
  expect_identical(got[["lti:test_Parallel_functions:1962"]], "Parallel(DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.01), DiscreteTransferFunction(k, 1, s, 0.01), DiscreteTransferFunction(a2*p - s, a2*s + p, s, 0.01))")
  expect_identical(got[["lti:test_Parallel_functions:1963"]], "Parallel(DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.01), DiscreteTransferFunction(k, 1, s, 0.01), DiscreteTransferFunction(a2*p - s, a2*s + p, s, 0.01), DiscreteTransferFunction(-a0 + a1*s^2 + a2*s, a0 + s, s, 0.01))")
  expect_identical(got[["lti:test_Parallel_functions:1964"]], "Parallel(DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.01), DiscreteTransferFunction(k, 1, s, 0.01), DiscreteTransferFunction(-a2*p + s, a2*s + p, s, 0.01), DiscreteTransferFunction(a0 - a1*s^2 - a2*s, a0 + s, s, 0.01))")
  expect_identical(got[["lti:test_Parallel_functions:1965"]], "Parallel(DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.01), Series(DiscreteTransferFunction(k, 1, s, 0.01), DiscreteTransferFunction(a2*p - s, a2*s + p, s, 0.01)))")
  expect_identical(got[["lti:test_Parallel_functions:1966"]], "Parallel(DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.01), Series(DiscreteTransferFunction(-1, 1, s, 0.01), DiscreteTransferFunction(k, 1, s, 0.01), DiscreteTransferFunction(a2*p - s, a2*s + p, s, 0.01)))")
  expect_identical(got[["lti:test_Parallel_functions:1967"]], "Parallel(DiscreteTransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s, 0.01), DiscreteTransferFunction(-k, 1, s, 0.01))")
  expect_identical(got[["lti:test_Parallel_functions:1968"]], "Series(DiscreteTransferFunction(-1, 1, s, 0.01), Parallel(DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.01), DiscreteTransferFunction(k, 1, s, 0.01)))")
  expect_identical(got[["lti:test_Parallel_functions:1970"]], "Series(Parallel(DiscreteTransferFunction(k, 1, s, 0.01), DiscreteTransferFunction(a2*p - s, a2*s + p, s, 0.01)), DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.01))")
  expect_identical(got[["lti:test_Parallel_functions:1971"]], "Series(Parallel(DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.01), DiscreteTransferFunction(k, 1, s, 0.01)), DiscreteTransferFunction(a2*p - s, a2*s + p, s, 0.01), DiscreteTransferFunction(-a0 + a1*s^2 + a2*s, a0 + s, s, 0.01))")
  expect_identical(got[["lti:test_Parallel_functions:1972"]], "Series(DiscreteTransferFunction(-1, 1, s, 0.01), Parallel(DiscreteTransferFunction(k, 1, s, 0.01), DiscreteTransferFunction(a2*p - s, a2*s + p, s, 0.01)), DiscreteTransferFunction(a0 - a1*s^2 - a2*s, a0 + s, s, 0.01))")
  expect_identical(got[["lti:test_Parallel_functions:1974"]], "Parallel(DiscreteTransferFunction(k, 1, s, 0.01), DiscreteTransferFunction(a2*p - s, a2*s + p, s, 0.01), Series(DiscreteTransferFunction(k, 1, s, 0.01), DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.01)), DiscreteTransferFunction(-a0 + a1*s^2 + a2*s, a0 + s, s, 0.01))")
  expect_identical(got[["lti:test_Parallel_functions:1976"]], "Parallel(DiscreteTransferFunction(k, 1, s, 0.01), DiscreteTransferFunction(a2*p - s, a2*s + p, s, 0.01), Series(DiscreteTransferFunction(k, 1, s, 0.01), DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.01)), DiscreteTransferFunction(-a2*p + s, a2*s + p, s, 0.01))")
  expect_identical(got[["lti:test_Parallel_functions:1978"]], "Series(Parallel(DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.01), DiscreteTransferFunction(k, 1, s, 0.01), DiscreteTransferFunction(-a0 + a1*s^2 + a2*s, a0 + s, s, 0.01)), Parallel(DiscreteTransferFunction(a2*p - s, a2*s + p, s, 0.01), DiscreteTransferFunction(-a0 + a1*s^2 + a2*s, a0 + s, s, 0.01)))")
  expect_identical(got[["lti:test_Parallel_functions:1981"]], "ERR ValueError: \\nAll transfer functions should use the same complex variable of the\\nLaplace transform or z-transform. 2 different values found.")
  expect_identical(got[["lti:test_Parallel_functions:1982"]], "ERR ValueError: \\nAll the transfer functions should use the same complex variable of the\\nLaplace domain or z-domain.")
  expect_identical(got[["lti:test_Parallel_functions:1983"]], "ERR ValueError: <class 'sympy.physics.control.lti.DiscreteTransferFunction'> cannot be added with <class 'sympy.matrices.dense.MutableDenseMatrix'>.")
  expect_identical(got[["lti:test_Parallel_functions:1986"]], "True")
  expect_identical(got[["lti:test_Parallel_functions:1990"]], "True")
  expect_identical(got[["lti:test_Parallel_functions:1996"]], "True")
  expect_identical(got[["lti:test_Parallel_functions:2001"]], "Parallel(DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 0.01), DiscreteTransferFunction(-k, 1, s, 0.01))")
  expect_identical(got[["lti:test_Parallel_functions:2004"]], "DiscreteTransferFunction(k*(a2*p - s)*(s^2 + 2*s*wn*zeta + wn^2) + k*(a2*s + p), (a2*s + p)*(s^2 + 2*s*wn*zeta + wn^2), s, 0.01)")
  expect_identical(got[["lti:test_Parallel_functions:2008"]], "DiscreteTransferFunction(-a2*s - k*(a2*s + p)*(s^2 + 2*s*wn*zeta + wn^2) - p + (-a2*p + s)*(s^2 + 2*s*wn*zeta + wn^2), (a2*s + p)*(s^2 + 2*s*wn*zeta + wn^2), s, 0.01)")
  expect_identical(got[["lti:test_Parallel_functions:2012"]], "DiscreteTransferFunction(-a2*s - k*(a2*s + p)*(s^2 + 2*s*wn*zeta + wn^2) - p - (a2*p - s)*(s^2 + 2*s*wn*zeta + wn^2), (a2*s + p)*(s^2 + 2*s*wn*zeta + wn^2), s, 0.01)")
  expect_identical(got[["lti:test_Parallel_functions:2016"]], "DiscreteTransferFunction(k*(a2*s + p)*(s^2 + 2*s*wn*zeta + wn^2) - k*(a2*s + p) + (2*a2*p - 2*s)*(s^2 + 2*s*wn*zeta + wn^2), (a2*s + p)*(s^2 + 2*s*wn*zeta + wn^2), s, 0.01)")
  expect_identical(got[["lti:test_Parallel_functions:2021"]], "DiscreteTransferFunction(k*(s^2 + 2*s*wn*zeta + wn^2) + 1, s^2 + 2*s*wn*zeta + wn^2, s, 0.01)")
  expect_identical(got[["lti:test_Parallel_functions:2024"]], "DiscreteTransferFunction(a2*s + k*(a2*s + p)*(s^2 + 2*s*wn*zeta + wn^2) + p + (-a2*p + s)*(s^2 + 2*s*wn*zeta + wn^2), (a2*s + p)*(s^2 + 2*s*wn*zeta + wn^2), s, 0.01)")
  expect_identical(got[["lti:test_Parallel_functions:2029"]], "True")
  expect_identical(got[["lti:test_Parallel_functions:2033"]], "True")
  expect_identical(got[["lti:test_Parallel_functions:2034"]], "False")
  expect_identical(got[["lti:test_Parallel_functions:2035"]], "True")
  expect_identical(got[["lti:test_Parallel_functions:2038"]], "True")
  expect_identical(got[["lti:test_Parallel_functions:2039"]], "False")
  expect_identical(got[["lti:test_Parallel_functions:2040"]], "True")
  expect_identical(got[["lti:test_Parallel_functions:2043"]], "True")
  expect_identical(got[["lti:test_Parallel_functions:2044"]], "False")
  expect_identical(got[["lti:test_Parallel_functions:2045"]], "True")
  expect_identical(got[["lti:test_Parallel_functions:2048"]], "ERR TypeError: \\n                    Cannot rewrite a discrete-time Parallel object as a\\n                    TransferFunction.")
  expect_identical(got[["lti:test_Parallel_functions:2050"]], "ERR TypeError: \\n                    Cannot rewrite a continuous-time Parallel object as a\\n                    DiscreteTransferFunction.")
})

test_that("lti test_MIMOParallel_functions matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("a, x, b, c, s, g, d,    p, k, tau, zeta, wn, T, z"); a <- .py_getitem(.tmpu_, 0L); x <- .py_getitem(.tmpu_, 1L); b <- .py_getitem(.tmpu_, 2L); c <- .py_getitem(.tmpu_, 3L); s <- .py_getitem(.tmpu_, 4L); g <- .py_getitem(.tmpu_, 5L); d <- .py_getitem(.tmpu_, 6L); p <- .py_getitem(.tmpu_, 7L); k <- .py_getitem(.tmpu_, 8L); tau <- .py_getitem(.tmpu_, 9L); zeta <- .py_getitem(.tmpu_, 10L); wn <- .py_getitem(.tmpu_, 11L); T <- .py_getitem(.tmpu_, 12L); z <- .py_getitem(.tmpu_, 13L)
})
.py_try({
.tmpu_ <- .py_symbols("a0:4, b0:5, c0:4, d0:4"); a0 <- .py_getitem(.tmpu_, 0L); a1 <- .py_getitem(.tmpu_, 1L); a2 <- .py_getitem(.tmpu_, 2L); a3 <- .py_getitem(.tmpu_, 3L); b0 <- .py_getitem(.tmpu_, 4L); b1 <- .py_getitem(.tmpu_, 5L); b2 <- .py_getitem(.tmpu_, 6L); b3 <- .py_getitem(.tmpu_, 7L); b4 <- .py_getitem(.tmpu_, 8L); c0 <- .py_getitem(.tmpu_, 9L); c1 <- .py_getitem(.tmpu_, 10L); c2 <- .py_getitem(.tmpu_, 11L); c3 <- .py_getitem(.tmpu_, 12L); d0 <- .py_getitem(.tmpu_, 13L); d1 <- .py_getitem(.tmpu_, 14L); d2 <- .py_getitem(.tmpu_, 15L); d3 <- .py_getitem(.tmpu_, 16L)
})
.py_try({
TF1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
TF2 <- TransferFunction(k, 1L, s)
})
.py_try({
TF3 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
TF4 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(tau, s), 1L), s)
})
.py_try({
TF5 <- TransferFunction(k, p, s)
})
.py_try({
TF6 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(k, s), p), s)
})
p_are_floats_equals <- function(n1, n2) {
  "\n    Helper function to compare two floats.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return(.py_method(Float(n1), "is_same", Float(n2), .py_isclose))
}
p_are_complex_equals <- function(n1, n2) {
  "\n    Helper function to compare two complex numbers.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return((.py_truthy(p_are_floats_equals(re(n1), re(n2))) && .py_truthy(p_are_floats_equals(im(n1), im(n2)))))
}
.py_try({
tf4 <- TransferFunction(.py_sub(.py_add(.py_mul(a0, p), .py_pow(p, a1)), s), p, p)
})
.py_try({
tf5 <- TransferFunction(.py_sub(.py_add(.py_mul(a1, .py_pow(s, 2L)), .py_mul(a2, s)), a0), .py_add(s, a0), s)
})
.py_try({
dtf4 <- DiscreteTransferFunction(.py_sub(.py_add(.py_mul(a0, p), .py_pow(p, a1)), s), p, p, Rational(1L, 2L))
})
.py_try({
dtf5 <- DiscreteTransferFunction(.py_sub(.py_add(.py_mul(a1, .py_pow(s, 2L)), .py_mul(a2, s)), a0), .py_add(s, a0), s, Rational(1L, 2L))
})
.py_try({
tfm1 <- TransferFunctionMatrix(list(list(TF1), list(TF2), list(TF3)))
})
.py_try({
tfm2 <- TransferFunctionMatrix(list(list(.py_neg(TF2)), list(tf5), list(.py_neg(TF1))))
})
.py_try({
tfm3 <- TransferFunctionMatrix(list(list(tf5), list(.py_neg(tf5)), list(TF2)))
})
.py_try({
tfm4 <- TransferFunctionMatrix(list(list(TF2, .py_neg(tf5)), list(TF1, tf5)))
})
.py_try({
tfm5 <- TransferFunctionMatrix(list(list(TF1, TF2), list(TF3, .py_neg(tf5))))
})
.py_try({
tfm6 <- TransferFunctionMatrix(list(list(.py_neg(TF2))))
})
.py_try({
tfm7 <- TransferFunctionMatrix(list(list(tf4), list(.py_neg(tf4)), list(tf4)))
})
.py_try({
dTF1 <- DiscreteTransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s, Rational(1L, 2L))
})
.py_try({
dTF2 <- DiscreteTransferFunction(k, 1L, s, Rational(1L, 2L))
})
.py_try({
dTF3 <- DiscreteTransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s, Rational(1L, 2L))
})
.py_try({
dtfm1 <- TransferFunctionMatrix(list(list(dTF1), list(dTF2), list(dTF3)))
})
.py_try({
dtfm2 <- TransferFunctionMatrix(list(list(.py_neg(dTF2)), list(dtf5), list(.py_neg(dTF1))))
})
.py_try({
dtfm3 <- TransferFunctionMatrix(list(list(dtf5), list(.py_neg(dtf5)), list(dTF2)))
})
.py_try({
dtfm4 <- TransferFunctionMatrix(list(list(dTF2, .py_neg(dtf5)), list(dTF1, dtf5)))
})
.py_try({
dtfm5 <- TransferFunctionMatrix(list(list(dTF1, dTF2), list(dTF3, .py_neg(dtf5))))
})
.py_try({
dtfm6 <- TransferFunctionMatrix(list(list(.py_neg(dTF2))))
})
.py_try({
dtfm7 <- TransferFunctionMatrix(list(list(dtf4), list(.py_neg(dtf4)), list(dtf4)))
})
.py_try({
.py_rec("lti:test_MIMOParallel_functions:2082", function() (.py_truthy(.py_eq(.py_add(.py_add(tfm1, tfm2), tfm3), MIMOParallel(tfm1, tfm2, tfm3))) && .py_truthy(.py_eq(MIMOParallel(tfm1, tfm2, tfm3), MIMOParallel(MIMOParallel(tfm1, tfm2), tfm3)))))
})
.py_try({
.py_rec("lti:test_MIMOParallel_functions:2083", function() .py_sub(.py_sub(tfm2, tfm1), tfm3))
})
.py_try({
.py_rec("lti:test_MIMOParallel_functions:2084", function() .py_add(.py_sub(tfm2, tfm3), .py_mul(.py_mul(.py_neg(tfm1), tfm6), .py_neg(tfm6))))
})
.py_try({
.py_rec("lti:test_MIMOParallel_functions:2085", function() .py_sub(.py_add(tfm1, tfm1), .py_mul(.py_neg(tfm1), tfm6)))
})
.py_try({
.py_rec("lti:test_MIMOParallel_functions:2086", function() .py_add(.py_sub(.py_sub(tfm2, tfm3), tfm1), tfm2))
})
.py_try({
.py_rec("lti:test_MIMOParallel_functions:2087", function() .py_sub(.py_sub(.py_add(tfm1, tfm2), tfm3), tfm1))
})
.py_try({
.py_rec("lti:test_MIMOParallel_functions:2088", function() .py_add(.py_add(tfm1, tfm2), TF2))
})
.py_try({
.py_rec("lti:test_MIMOParallel_functions:2089", function() .py_sub(.py_sub(tfm1, tfm2), a1))
})
.py_try({
.py_rec("lti:test_MIMOParallel_functions:2090", function() .py_sub(.py_sub(tfm2, tfm3), .py_sub(s, 1L)))
})
.py_try({
.py_rec("lti:test_MIMOParallel_functions:2091", function() .py_sub(.py_sub(.py_neg(tfm3), tfm2), 9L))
})
.py_try({
.py_rec("lti:test_MIMOParallel_functions:2092", function() .py_sub(.py_sub(.py_sub(1L, .py_pow(p, 3L)), tfm3), tfm2))
})
.py_try({
.py_rec("lti:test_MIMOParallel_functions:2094", function() .py_sub(.py_sub(tfm3, tfm2), tfm7))
})
.py_try({
.py_rec("lti:test_MIMOParallel_functions:2095", function() .py_add(.py_sub(tfm2, tfm1), tfm7))
})
.py_try({
.py_rec("lti:test_MIMOParallel_functions:2097", function() .py_add(.py_add(tfm1, tfm2), tfm4))
})
.py_try({
.py_rec("lti:test_MIMOParallel_functions:2098", function() .py_sub(.py_sub(tfm1, tfm2), tfm4))
})
.py_try({
.py_rec("lti:test_MIMOParallel_functions:2100", function() .py_mul(.py_add(tfm1, tfm2), tfm6))
})
.py_try({
.py_rec("lti:test_MIMOParallel_functions:2101", function() .py_mul(.py_mul(.py_sub(tfm2, tfm3), tfm6), .py_neg(tfm6)))
})
.py_try({
.py_rec("lti:test_MIMOParallel_functions:2102", function() .py_mul(.py_sub(.py_sub(tfm2, tfm1), tfm3), .py_add(tfm6, tfm6)))
})
.py_try({
.py_rec("lti:test_MIMOParallel_functions:2103", function() .py_mul(.py_add(tfm4, tfm5), TF1))
})
.py_try({
.py_rec("lti:test_MIMOParallel_functions:2104", function() .py_mul(.py_sub(tfm2, tfm3), a2))
})
.py_try({
.py_rec("lti:test_MIMOParallel_functions:2105", function() .py_mul(.py_add(tfm3, tfm2), .py_sub(s, 6L)))
})
.py_try({
.py_rec("lti:test_MIMOParallel_functions:2106", function() .py_mul(.py_add(.py_add(tfm1, tfm2), tfm3), 0L))
})
.py_try({
.py_rec("lti:test_MIMOParallel_functions:2107", function() .py_mul(.py_sub(1L, .py_pow(p, 3L)), .py_add(tfm1, tfm3)))
})
.py_try({
.py_rec("lti:test_MIMOParallel_functions:2110", function() .py_mul(.py_mul(.py_sub(tfm3, tfm2), tfm4), tfm5))
})
.py_try({
.py_rec("lti:test_MIMOParallel_functions:2112", function() .py_mul(.py_sub(tfm1, tfm2), tfm5))
})
.py_try({
.py_rec("lti:test_MIMOParallel_functions:2115", function() (.py_truthy(.py_eq(MIMOParallel(tfm1, tfm2, evaluate = TRUE), .py_method(MIMOParallel(tfm1, tfm2), "doit"))) && .py_truthy(.py_eq(.py_method(MIMOParallel(tfm1, tfm2), "doit"), .py_method(MIMOParallel(tfm1, tfm2), "rewrite", TransferFunctionMatrix))) && .py_truthy(.py_eq(.py_method(MIMOParallel(tfm1, tfm2), "rewrite", TransferFunctionMatrix), TransferFunctionMatrix(.py_tup(.py_tup(TransferFunction(.py_add(.py_mul(.py_neg(k), .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L))), 1L), .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L)), s)), .py_tup(TransferFunction(.py_add(.py_add(.py_add(.py_neg(a0), .py_mul(a1, .py_pow(s, 2L))), .py_mul(a2, s)), .py_mul(k, .py_add(a0, s))), .py_add(a0, s), s)), .py_tup(TransferFunction(.py_add(.py_sub(.py_mul(.py_neg(a2), s), p), .py_mul(.py_sub(.py_mul(a2, p), s), .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L)))), .py_mul(.py_add(.py_mul(a2, s), p), .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L))), s))))))))
})
.py_try({
.py_rec("lti:test_MIMOParallel_functions:2126", function() (.py_truthy(.py_eq(.py_add(.py_add(dtfm1, dtfm2), dtfm3), MIMOParallel(dtfm1, dtfm2, dtfm3))) && .py_truthy(.py_eq(MIMOParallel(dtfm1, dtfm2, dtfm3), MIMOParallel(MIMOParallel(dtfm1, dtfm2), dtfm3)))))
})
.py_try({
.py_rec("lti:test_MIMOParallel_functions:2128", function() .py_sub(.py_sub(dtfm2, dtfm1), dtfm3))
})
.py_try({
.py_rec("lti:test_MIMOParallel_functions:2129", function() .py_add(.py_sub(dtfm2, dtfm3), .py_mul(.py_mul(.py_neg(dtfm1), dtfm6), .py_neg(dtfm6))))
})
.py_try({
.py_rec("lti:test_MIMOParallel_functions:2131", function() .py_sub(.py_add(dtfm1, dtfm1), .py_mul(.py_neg(dtfm1), dtfm6)))
})
.py_try({
.py_rec("lti:test_MIMOParallel_functions:2133", function() .py_add(.py_sub(.py_sub(dtfm2, dtfm3), dtfm1), dtfm2))
})
.py_try({
.py_rec("lti:test_MIMOParallel_functions:2135", function() .py_sub(.py_sub(.py_add(dtfm1, dtfm2), dtfm3), dtfm1))
})
.py_try({
.py_rec("lti:test_MIMOParallel_functions:2137", function() .py_add(.py_add(dtfm1, dtfm2), dTF2))
})
.py_try({
.py_rec("lti:test_MIMOParallel_functions:2138", function() .py_sub(.py_sub(dtfm1, dtfm2), a1))
})
.py_try({
.py_rec("lti:test_MIMOParallel_functions:2139", function() .py_sub(.py_sub(dtfm2, dtfm3), .py_sub(s, 1L)))
})
.py_try({
.py_rec("lti:test_MIMOParallel_functions:2140", function() .py_sub(.py_sub(.py_neg(dtfm3), dtfm2), 9L))
})
.py_try({
.py_rec("lti:test_MIMOParallel_functions:2141", function() .py_sub(.py_sub(.py_sub(1L, .py_pow(p, 3L)), dtfm3), dtfm2))
})
.py_try({
.py_rec("lti:test_MIMOParallel_functions:2143", function() .py_sub(.py_sub(dtfm3, dtfm2), dtfm7))
})
.py_try({
.py_rec("lti:test_MIMOParallel_functions:2144", function() .py_add(.py_sub(dtfm2, dtfm1), dtfm7))
})
.py_try({
.py_rec("lti:test_MIMOParallel_functions:2146", function() .py_add(.py_add(dtfm1, dtfm2), dtfm4))
})
.py_try({
.py_rec("lti:test_MIMOParallel_functions:2147", function() .py_sub(.py_sub(dtfm1, dtfm2), dtfm4))
})
.py_try({
.py_rec("lti:test_MIMOParallel_functions:2149", function() .py_mul(.py_add(dtfm1, dtfm2), dtfm6))
})
.py_try({
.py_rec("lti:test_MIMOParallel_functions:2151", function() .py_mul(.py_mul(.py_sub(dtfm2, dtfm3), dtfm6), .py_neg(dtfm6)))
})
.py_try({
.py_rec("lti:test_MIMOParallel_functions:2153", function() .py_mul(.py_sub(.py_sub(dtfm2, dtfm1), dtfm3), .py_add(dtfm6, dtfm6)))
})
.py_try({
.py_rec("lti:test_MIMOParallel_functions:2156", function() .py_mul(.py_add(dtfm4, dtfm5), dTF1))
})
.py_try({
.py_rec("lti:test_MIMOParallel_functions:2157", function() .py_mul(.py_sub(dtfm2, dtfm3), a2))
})
.py_try({
.py_rec("lti:test_MIMOParallel_functions:2158", function() .py_mul(.py_add(dtfm3, dtfm2), .py_sub(s, 6L)))
})
.py_try({
.py_rec("lti:test_MIMOParallel_functions:2159", function() .py_mul(.py_add(.py_add(dtfm1, dtfm2), dtfm3), 0L))
})
.py_try({
.py_rec("lti:test_MIMOParallel_functions:2160", function() .py_mul(.py_sub(1L, .py_pow(p, 3L)), .py_add(dtfm1, dtfm3)))
})
.py_try({
.py_rec("lti:test_MIMOParallel_functions:2163", function() .py_mul(.py_mul(.py_sub(dtfm3, dtfm2), dtfm4), dtfm5))
})
.py_try({
.py_rec("lti:test_MIMOParallel_functions:2165", function() .py_mul(.py_sub(dtfm1, dtfm2), dtfm5))
})
.py_try({
.py_rec("lti:test_MIMOParallel_functions:2168", function() (.py_truthy(.py_eq(MIMOParallel(dtfm1, dtfm2, evaluate = TRUE), .py_method(MIMOParallel(dtfm1, dtfm2), "doit"))) && .py_truthy(.py_eq(.py_method(MIMOParallel(dtfm1, dtfm2), "doit"), .py_method(MIMOParallel(dtfm1, dtfm2), "rewrite", TransferFunctionMatrix))) && .py_truthy(.py_eq(.py_method(MIMOParallel(dtfm1, dtfm2), "rewrite", TransferFunctionMatrix), TransferFunctionMatrix(.py_tup(.py_tup(DiscreteTransferFunction(.py_add(.py_mul(.py_neg(k), .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L))), 1L), .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L)), s, Rational(1L, 2L))), .py_tup(DiscreteTransferFunction(.py_add(.py_add(.py_add(.py_neg(a0), .py_mul(a1, .py_pow(s, 2L))), .py_mul(a2, s)), .py_mul(k, .py_add(a0, s))), .py_add(a0, s), s, Rational(1L, 2L))), .py_tup(DiscreteTransferFunction(.py_add(.py_sub(.py_mul(.py_neg(a2), s), p), .py_mul(.py_sub(.py_mul(a2, p), s), .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L)))), .py_mul(.py_add(.py_mul(a2, s), p), .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, s), wn), zeta)), .py_pow(wn, 2L))), s, Rational(1L, 2L)))))))))
}))---")
  expect_identical(got[["lti:test_MIMOParallel_functions:2082"]], "True")
  expect_identical(got[["lti:test_MIMOParallel_functions:2083"]], "MIMOParallel(TransferFunctionMatrix(((TransferFunction(-k, 1, s),), (TransferFunction(-a0 + a1*s^2 + a2*s, a0 + s, s),), (TransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s),))), TransferFunctionMatrix(((TransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s),), (TransferFunction(-k, 1, s),), (TransferFunction(-a2*p + s, a2*s + p, s),))), TransferFunctionMatrix(((TransferFunction(a0 - a1*s^2 - a2*s, a0 + s, s),), (TransferFunction(-a0 + a1*s^2 + a2*s, a0 + s, s),), (TransferFunction(-k, 1, s),))))")
  expect_identical(got[["lti:test_MIMOParallel_functions:2084"]], "MIMOParallel(TransferFunctionMatrix(((TransferFunction(-k, 1, s),), (TransferFunction(-a0 + a1*s^2 + a2*s, a0 + s, s),), (TransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s),))), TransferFunctionMatrix(((TransferFunction(a0 - a1*s^2 - a2*s, a0 + s, s),), (TransferFunction(-a0 + a1*s^2 + a2*s, a0 + s, s),), (TransferFunction(-k, 1, s),))), MIMOSeries(TransferFunctionMatrix(((TransferFunction(k, 1, s),),)), TransferFunctionMatrix(((TransferFunction(-k, 1, s),),)), TransferFunctionMatrix(((TransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s),), (TransferFunction(-k, 1, s),), (TransferFunction(-a2*p + s, a2*s + p, s),)))))")
  expect_identical(got[["lti:test_MIMOParallel_functions:2085"]], "MIMOParallel(TransferFunctionMatrix(((TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s),), (TransferFunction(k, 1, s),), (TransferFunction(a2*p - s, a2*s + p, s),))), TransferFunctionMatrix(((TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s),), (TransferFunction(k, 1, s),), (TransferFunction(a2*p - s, a2*s + p, s),))), MIMOSeries(TransferFunctionMatrix(((TransferFunction(k, 1, s),),)), TransferFunctionMatrix(((TransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s),), (TransferFunction(-k, 1, s),), (TransferFunction(-a2*p + s, a2*s + p, s),)))))")
  expect_identical(got[["lti:test_MIMOParallel_functions:2086"]], "MIMOParallel(TransferFunctionMatrix(((TransferFunction(-k, 1, s),), (TransferFunction(-a0 + a1*s^2 + a2*s, a0 + s, s),), (TransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s),))), TransferFunctionMatrix(((TransferFunction(a0 - a1*s^2 - a2*s, a0 + s, s),), (TransferFunction(-a0 + a1*s^2 + a2*s, a0 + s, s),), (TransferFunction(-k, 1, s),))), TransferFunctionMatrix(((TransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s),), (TransferFunction(-k, 1, s),), (TransferFunction(-a2*p + s, a2*s + p, s),))), TransferFunctionMatrix(((TransferFunction(-k, 1, s),), (TransferFunction(-a0 + a1*s^2 + a2*s, a0 + s, s),), (TransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s),))))")
  expect_identical(got[["lti:test_MIMOParallel_functions:2087"]], "MIMOParallel(TransferFunctionMatrix(((TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s),), (TransferFunction(k, 1, s),), (TransferFunction(a2*p - s, a2*s + p, s),))), TransferFunctionMatrix(((TransferFunction(-k, 1, s),), (TransferFunction(-a0 + a1*s^2 + a2*s, a0 + s, s),), (TransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s),))), TransferFunctionMatrix(((TransferFunction(a0 - a1*s^2 - a2*s, a0 + s, s),), (TransferFunction(-a0 + a1*s^2 + a2*s, a0 + s, s),), (TransferFunction(-k, 1, s),))), TransferFunctionMatrix(((TransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s),), (TransferFunction(-k, 1, s),), (TransferFunction(-a2*p + s, a2*s + p, s),))))")
  expect_identical(got[["lti:test_MIMOParallel_functions:2088"]], "ERR ValueError: <class 'sympy.physics.control.lti.TransferFunction'> cannot be added with <class 'sympy.physics.control.lti.MIMOParallel'>.")
  expect_identical(got[["lti:test_MIMOParallel_functions:2089"]], "ERR TypeError: unsupported operand type(s) for -: 'MIMOParallel' and 'Symbol'")
  expect_identical(got[["lti:test_MIMOParallel_functions:2090"]], "ERR TypeError: unsupported operand type(s) for -: 'MIMOParallel' and 'Add'")
  expect_identical(got[["lti:test_MIMOParallel_functions:2091"]], "ERR TypeError: unsupported operand type(s) for -: 'MIMOParallel' and 'int'")
  expect_identical(got[["lti:test_MIMOParallel_functions:2092"]], "ERR TypeError: unsupported operand type(s) for -: 'Add' and 'TransferFunctionMatrix'")
  expect_identical(got[["lti:test_MIMOParallel_functions:2094"]], "ERR ValueError: \\nAll transfer functions should use the same complex variable of the\\nLaplace transform or z-transform. 2 different values found.")
  expect_identical(got[["lti:test_MIMOParallel_functions:2095"]], "ERR ValueError: \\nAll transfer functions should use the same complex variable of the\\nLaplace transform or z-transform. 2 different values found.")
  expect_identical(got[["lti:test_MIMOParallel_functions:2097"]], "ERR TypeError: Shape of all the args is not equal.")
  expect_identical(got[["lti:test_MIMOParallel_functions:2098"]], "ERR TypeError: Shape of all the args is not equal.")
  expect_identical(got[["lti:test_MIMOParallel_functions:2100"]], "MIMOSeries(TransferFunctionMatrix(((TransferFunction(-k, 1, s),),)), MIMOParallel(TransferFunctionMatrix(((TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s),), (TransferFunction(k, 1, s),), (TransferFunction(a2*p - s, a2*s + p, s),))), TransferFunctionMatrix(((TransferFunction(-k, 1, s),), (TransferFunction(-a0 + a1*s^2 + a2*s, a0 + s, s),), (TransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s),)))))")
  expect_identical(got[["lti:test_MIMOParallel_functions:2101"]], "MIMOSeries(TransferFunctionMatrix(((TransferFunction(k, 1, s),),)), TransferFunctionMatrix(((TransferFunction(-k, 1, s),),)), MIMOParallel(TransferFunctionMatrix(((TransferFunction(-k, 1, s),), (TransferFunction(-a0 + a1*s^2 + a2*s, a0 + s, s),), (TransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s),))), TransferFunctionMatrix(((TransferFunction(a0 - a1*s^2 - a2*s, a0 + s, s),), (TransferFunction(-a0 + a1*s^2 + a2*s, a0 + s, s),), (TransferFunction(-k, 1, s),)))))")
  expect_identical(got[["lti:test_MIMOParallel_functions:2102"]], "MIMOSeries(MIMOParallel(TransferFunctionMatrix(((TransferFunction(-k, 1, s),),)), TransferFunctionMatrix(((TransferFunction(-k, 1, s),),))), MIMOParallel(TransferFunctionMatrix(((TransferFunction(-k, 1, s),), (TransferFunction(-a0 + a1*s^2 + a2*s, a0 + s, s),), (TransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s),))), TransferFunctionMatrix(((TransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s),), (TransferFunction(-k, 1, s),), (TransferFunction(-a2*p + s, a2*s + p, s),))), TransferFunctionMatrix(((TransferFunction(a0 - a1*s^2 - a2*s, a0 + s, s),), (TransferFunction(-a0 + a1*s^2 + a2*s, a0 + s, s),), (TransferFunction(-k, 1, s),)))))")
  expect_identical(got[["lti:test_MIMOParallel_functions:2103"]], "ERR ValueError: <class 'sympy.physics.control.lti.TransferFunction'> cannot be multiplied with <class 'sympy.physics.control.lti.MIMOParallel'>.")
  expect_identical(got[["lti:test_MIMOParallel_functions:2104"]], "ERR TypeError: unsupported operand type(s) for *: 'MIMOParallel' and 'Symbol'")
  expect_identical(got[["lti:test_MIMOParallel_functions:2105"]], "ERR TypeError: unsupported operand type(s) for *: 'MIMOParallel' and 'Add'")
  expect_identical(got[["lti:test_MIMOParallel_functions:2106"]], "ERR TypeError: unsupported operand type(s) for *: 'MIMOParallel' and 'int'")
  expect_identical(got[["lti:test_MIMOParallel_functions:2107"]], "ERR TypeError: unsupported operand type(s) for *: 'Add' and 'MIMOParallel'")
  expect_identical(got[["lti:test_MIMOParallel_functions:2110"]], "ERR ValueError: \\nNumber of input signals do not match the number of output signals of\\nadjacent systems for some args.")
  expect_identical(got[["lti:test_MIMOParallel_functions:2112"]], "ERR ValueError: \\nNumber of input signals do not match the number of output signals of\\nadjacent systems for some args.")
  expect_identical(got[["lti:test_MIMOParallel_functions:2115"]], "True")
  expect_identical(got[["lti:test_MIMOParallel_functions:2126"]], "True")
  expect_identical(got[["lti:test_MIMOParallel_functions:2128"]], "MIMOParallel(TransferFunctionMatrix(((DiscreteTransferFunction(-k, 1, s, 1/2),), (DiscreteTransferFunction(-a0 + a1*s^2 + a2*s, a0 + s, s, 1/2),), (DiscreteTransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s, 1/2),))), TransferFunctionMatrix(((DiscreteTransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s, 1/2),), (DiscreteTransferFunction(-k, 1, s, 1/2),), (DiscreteTransferFunction(-a2*p + s, a2*s + p, s, 1/2),))), TransferFunctionMatrix(((DiscreteTransferFunction(a0 - a1*s^2 - a2*s, a0 + s, s, 1/2),), (DiscreteTransferFunction(-a0 + a1*s^2 + a2*s, a0 + s, s, 1/2),), (DiscreteTransferFunction(-k, 1, s, 1/2),))))")
  expect_identical(got[["lti:test_MIMOParallel_functions:2129"]], "MIMOParallel(TransferFunctionMatrix(((DiscreteTransferFunction(-k, 1, s, 1/2),), (DiscreteTransferFunction(-a0 + a1*s^2 + a2*s, a0 + s, s, 1/2),), (DiscreteTransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s, 1/2),))), TransferFunctionMatrix(((DiscreteTransferFunction(a0 - a1*s^2 - a2*s, a0 + s, s, 1/2),), (DiscreteTransferFunction(-a0 + a1*s^2 + a2*s, a0 + s, s, 1/2),), (DiscreteTransferFunction(-k, 1, s, 1/2),))), MIMOSeries(TransferFunctionMatrix(((DiscreteTransferFunction(k, 1, s, 1/2),),)), TransferFunctionMatrix(((DiscreteTransferFunction(-k, 1, s, 1/2),),)), TransferFunctionMatrix(((DiscreteTransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s, 1/2),), (DiscreteTransferFunction(-k, 1, s, 1/2),), (DiscreteTransferFunction(-a2*p + s, a2*s + p, s, 1/2),)))))")
  expect_identical(got[["lti:test_MIMOParallel_functions:2131"]], "MIMOParallel(TransferFunctionMatrix(((DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 1/2),), (DiscreteTransferFunction(k, 1, s, 1/2),), (DiscreteTransferFunction(a2*p - s, a2*s + p, s, 1/2),))), TransferFunctionMatrix(((DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 1/2),), (DiscreteTransferFunction(k, 1, s, 1/2),), (DiscreteTransferFunction(a2*p - s, a2*s + p, s, 1/2),))), MIMOSeries(TransferFunctionMatrix(((DiscreteTransferFunction(k, 1, s, 1/2),),)), TransferFunctionMatrix(((DiscreteTransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s, 1/2),), (DiscreteTransferFunction(-k, 1, s, 1/2),), (DiscreteTransferFunction(-a2*p + s, a2*s + p, s, 1/2),)))))")
  expect_identical(got[["lti:test_MIMOParallel_functions:2133"]], "MIMOParallel(TransferFunctionMatrix(((DiscreteTransferFunction(-k, 1, s, 1/2),), (DiscreteTransferFunction(-a0 + a1*s^2 + a2*s, a0 + s, s, 1/2),), (DiscreteTransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s, 1/2),))), TransferFunctionMatrix(((DiscreteTransferFunction(a0 - a1*s^2 - a2*s, a0 + s, s, 1/2),), (DiscreteTransferFunction(-a0 + a1*s^2 + a2*s, a0 + s, s, 1/2),), (DiscreteTransferFunction(-k, 1, s, 1/2),))), TransferFunctionMatrix(((DiscreteTransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s, 1/2),), (DiscreteTransferFunction(-k, 1, s, 1/2),), (DiscreteTransferFunction(-a2*p + s, a2*s + p, s, 1/2),))), TransferFunctionMatrix(((DiscreteTransferFunction(-k, 1, s, 1/2),), (DiscreteTransferFunction(-a0 + a1*s^2 + a2*s, a0 + s, s, 1/2),), (DiscreteTransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s, 1/2),))))")
  expect_identical(got[["lti:test_MIMOParallel_functions:2135"]], "MIMOParallel(TransferFunctionMatrix(((DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 1/2),), (DiscreteTransferFunction(k, 1, s, 1/2),), (DiscreteTransferFunction(a2*p - s, a2*s + p, s, 1/2),))), TransferFunctionMatrix(((DiscreteTransferFunction(-k, 1, s, 1/2),), (DiscreteTransferFunction(-a0 + a1*s^2 + a2*s, a0 + s, s, 1/2),), (DiscreteTransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s, 1/2),))), TransferFunctionMatrix(((DiscreteTransferFunction(a0 - a1*s^2 - a2*s, a0 + s, s, 1/2),), (DiscreteTransferFunction(-a0 + a1*s^2 + a2*s, a0 + s, s, 1/2),), (DiscreteTransferFunction(-k, 1, s, 1/2),))), TransferFunctionMatrix(((DiscreteTransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s, 1/2),), (DiscreteTransferFunction(-k, 1, s, 1/2),), (DiscreteTransferFunction(-a2*p + s, a2*s + p, s, 1/2),))))")
  expect_identical(got[["lti:test_MIMOParallel_functions:2137"]], "ERR ValueError: <class 'sympy.physics.control.lti.DiscreteTransferFunction'> cannot be added with <class 'sympy.physics.control.lti.MIMOParallel'>.")
  expect_identical(got[["lti:test_MIMOParallel_functions:2138"]], "ERR TypeError: unsupported operand type(s) for -: 'MIMOParallel' and 'Symbol'")
  expect_identical(got[["lti:test_MIMOParallel_functions:2139"]], "ERR TypeError: unsupported operand type(s) for -: 'MIMOParallel' and 'Add'")
  expect_identical(got[["lti:test_MIMOParallel_functions:2140"]], "ERR TypeError: unsupported operand type(s) for -: 'MIMOParallel' and 'int'")
  expect_identical(got[["lti:test_MIMOParallel_functions:2141"]], "ERR TypeError: unsupported operand type(s) for -: 'Add' and 'TransferFunctionMatrix'")
  expect_identical(got[["lti:test_MIMOParallel_functions:2143"]], "ERR ValueError: \\nAll transfer functions should use the same complex variable of the\\nLaplace transform or z-transform. 2 different values found.")
  expect_identical(got[["lti:test_MIMOParallel_functions:2144"]], "ERR ValueError: \\nAll transfer functions should use the same complex variable of the\\nLaplace transform or z-transform. 2 different values found.")
  expect_identical(got[["lti:test_MIMOParallel_functions:2146"]], "ERR TypeError: Shape of all the args is not equal.")
  expect_identical(got[["lti:test_MIMOParallel_functions:2147"]], "ERR TypeError: Shape of all the args is not equal.")
  expect_identical(got[["lti:test_MIMOParallel_functions:2149"]], "MIMOSeries(TransferFunctionMatrix(((DiscreteTransferFunction(-k, 1, s, 1/2),),)), MIMOParallel(TransferFunctionMatrix(((DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 1/2),), (DiscreteTransferFunction(k, 1, s, 1/2),), (DiscreteTransferFunction(a2*p - s, a2*s + p, s, 1/2),))), TransferFunctionMatrix(((DiscreteTransferFunction(-k, 1, s, 1/2),), (DiscreteTransferFunction(-a0 + a1*s^2 + a2*s, a0 + s, s, 1/2),), (DiscreteTransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s, 1/2),)))))")
  expect_identical(got[["lti:test_MIMOParallel_functions:2151"]], "MIMOSeries(TransferFunctionMatrix(((DiscreteTransferFunction(k, 1, s, 1/2),),)), TransferFunctionMatrix(((DiscreteTransferFunction(-k, 1, s, 1/2),),)), MIMOParallel(TransferFunctionMatrix(((DiscreteTransferFunction(-k, 1, s, 1/2),), (DiscreteTransferFunction(-a0 + a1*s^2 + a2*s, a0 + s, s, 1/2),), (DiscreteTransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s, 1/2),))), TransferFunctionMatrix(((DiscreteTransferFunction(a0 - a1*s^2 - a2*s, a0 + s, s, 1/2),), (DiscreteTransferFunction(-a0 + a1*s^2 + a2*s, a0 + s, s, 1/2),), (DiscreteTransferFunction(-k, 1, s, 1/2),)))))")
  expect_identical(got[["lti:test_MIMOParallel_functions:2153"]], "MIMOSeries(MIMOParallel(TransferFunctionMatrix(((DiscreteTransferFunction(-k, 1, s, 1/2),),)), TransferFunctionMatrix(((DiscreteTransferFunction(-k, 1, s, 1/2),),))), MIMOParallel(TransferFunctionMatrix(((DiscreteTransferFunction(-k, 1, s, 1/2),), (DiscreteTransferFunction(-a0 + a1*s^2 + a2*s, a0 + s, s, 1/2),), (DiscreteTransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s, 1/2),))), TransferFunctionMatrix(((DiscreteTransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s, 1/2),), (DiscreteTransferFunction(-k, 1, s, 1/2),), (DiscreteTransferFunction(-a2*p + s, a2*s + p, s, 1/2),))), TransferFunctionMatrix(((DiscreteTransferFunction(a0 - a1*s^2 - a2*s, a0 + s, s, 1/2),), (DiscreteTransferFunction(-a0 + a1*s^2 + a2*s, a0 + s, s, 1/2),), (DiscreteTransferFunction(-k, 1, s, 1/2),)))))")
  expect_identical(got[["lti:test_MIMOParallel_functions:2156"]], "ERR ValueError: <class 'sympy.physics.control.lti.DiscreteTransferFunction'> cannot be multiplied with <class 'sympy.physics.control.lti.MIMOParallel'>.")
  expect_identical(got[["lti:test_MIMOParallel_functions:2157"]], "ERR TypeError: unsupported operand type(s) for *: 'MIMOParallel' and 'Symbol'")
  expect_identical(got[["lti:test_MIMOParallel_functions:2158"]], "ERR TypeError: unsupported operand type(s) for *: 'MIMOParallel' and 'Add'")
  expect_identical(got[["lti:test_MIMOParallel_functions:2159"]], "ERR TypeError: unsupported operand type(s) for *: 'MIMOParallel' and 'int'")
  expect_identical(got[["lti:test_MIMOParallel_functions:2160"]], "ERR TypeError: unsupported operand type(s) for *: 'Add' and 'MIMOParallel'")
  expect_identical(got[["lti:test_MIMOParallel_functions:2163"]], "ERR ValueError: \\nNumber of input signals do not match the number of output signals of\\nadjacent systems for some args.")
  expect_identical(got[["lti:test_MIMOParallel_functions:2165"]], "ERR ValueError: \\nNumber of input signals do not match the number of output signals of\\nadjacent systems for some args.")
  expect_identical(got[["lti:test_MIMOParallel_functions:2168"]], "True")
})

test_that("lti test_Feedback_construction matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("a, x, b, c, s, g, d,    p, k, tau, zeta, wn, T, z"); a <- .py_getitem(.tmpu_, 0L); x <- .py_getitem(.tmpu_, 1L); b <- .py_getitem(.tmpu_, 2L); c <- .py_getitem(.tmpu_, 3L); s <- .py_getitem(.tmpu_, 4L); g <- .py_getitem(.tmpu_, 5L); d <- .py_getitem(.tmpu_, 6L); p <- .py_getitem(.tmpu_, 7L); k <- .py_getitem(.tmpu_, 8L); tau <- .py_getitem(.tmpu_, 9L); zeta <- .py_getitem(.tmpu_, 10L); wn <- .py_getitem(.tmpu_, 11L); T <- .py_getitem(.tmpu_, 12L); z <- .py_getitem(.tmpu_, 13L)
})
.py_try({
.tmpu_ <- .py_symbols("a0:4, b0:5, c0:4, d0:4"); a0 <- .py_getitem(.tmpu_, 0L); a1 <- .py_getitem(.tmpu_, 1L); a2 <- .py_getitem(.tmpu_, 2L); a3 <- .py_getitem(.tmpu_, 3L); b0 <- .py_getitem(.tmpu_, 4L); b1 <- .py_getitem(.tmpu_, 5L); b2 <- .py_getitem(.tmpu_, 6L); b3 <- .py_getitem(.tmpu_, 7L); b4 <- .py_getitem(.tmpu_, 8L); c0 <- .py_getitem(.tmpu_, 9L); c1 <- .py_getitem(.tmpu_, 10L); c2 <- .py_getitem(.tmpu_, 11L); c3 <- .py_getitem(.tmpu_, 12L); d0 <- .py_getitem(.tmpu_, 13L); d1 <- .py_getitem(.tmpu_, 14L); d2 <- .py_getitem(.tmpu_, 15L); d3 <- .py_getitem(.tmpu_, 16L)
})
.py_try({
TF1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
TF2 <- TransferFunction(k, 1L, s)
})
.py_try({
TF3 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
TF4 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(tau, s), 1L), s)
})
.py_try({
TF5 <- TransferFunction(k, p, s)
})
.py_try({
TF6 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(k, s), p), s)
})
p_are_floats_equals <- function(n1, n2) {
  "\n    Helper function to compare two floats.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return(.py_method(Float(n1), "is_same", Float(n2), .py_isclose))
}
p_are_complex_equals <- function(n1, n2) {
  "\n    Helper function to compare two complex numbers.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return((.py_truthy(p_are_floats_equals(re(n1), re(n2))) && .py_truthy(p_are_floats_equals(im(n1), im(n2)))))
}
.py_try({
tf1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
tf2 <- TransferFunction(k, 1L, s)
})
.py_try({
tf3 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
tf4 <- TransferFunction(.py_sub(.py_add(.py_mul(a0, p), .py_pow(p, a1)), s), p, p)
})
.py_try({
tf5 <- TransferFunction(.py_sub(.py_add(.py_mul(a1, .py_pow(s, 2L)), .py_mul(a2, s)), a0), .py_add(s, a0), s)
})
.py_try({
tf6 <- TransferFunction(.py_sub(s, p), .py_add(p, s), p)
})
.py_try({
dtf1 <- DiscreteTransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s, 2L)
})
.py_try({
dtf2 <- DiscreteTransferFunction(k, 1L, s, 2L)
})
.py_try({
dtf3 <- DiscreteTransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s, 2L)
})
.py_try({
dtf4 <- DiscreteTransferFunction(.py_sub(.py_add(.py_mul(a0, p), .py_pow(p, a1)), s), p, p, 2L)
})
.py_try({
dtf5 <- DiscreteTransferFunction(.py_sub(.py_add(.py_mul(a1, .py_pow(s, 2L)), .py_mul(a2, s)), a0), .py_add(s, a0), s, 2L)
})
.py_try({
dtf6 <- DiscreteTransferFunction(.py_sub(s, p), .py_add(p, s), p, 2L)
})
.py_try({
.py_rec("lti:test_Feedback_construction:2196", function() .py_attr(tf1, "is_continuous"))
})
.py_try({
.py_rec("lti:test_Feedback_construction:2197", function() .py_attr(dtf1, "is_continuous"))
})
.py_try({
f1 <- Feedback(TransferFunction(1L, 1L, s), .py_mul(.py_mul(tf1, tf2), tf3))
})
.py_try({
.py_rec("lti:test_Feedback_construction:2201", function() .py_attr(f1, "args"))
})
.py_try({
.py_rec("lti:test_Feedback_construction:2202", function() .py_attr(f1, "sys1"))
})
.py_try({
.py_rec("lti:test_Feedback_construction:2203", function() .py_attr(f1, "sys2"))
})
.py_try({
.py_rec("lti:test_Feedback_construction:2204", function() .py_attr(f1, "var"))
})
.py_try({
f2 <- Feedback(tf1, .py_mul(tf2, tf3))
})
.py_try({
.py_rec("lti:test_Feedback_construction:2207", function() .py_attr(f2, "args"))
})
.py_try({
.py_rec("lti:test_Feedback_construction:2208", function() .py_attr(f2, "sys1"))
})
.py_try({
.py_rec("lti:test_Feedback_construction:2209", function() .py_attr(f2, "sys2"))
})
.py_try({
.py_rec("lti:test_Feedback_construction:2210", function() .py_attr(f2, "var"))
})
.py_try({
f3 <- Feedback(.py_mul(tf1, tf2), tf5)
})
.py_try({
.py_rec("lti:test_Feedback_construction:2213", function() .py_attr(f3, "args"))
})
.py_try({
.py_rec("lti:test_Feedback_construction:2214", function() .py_attr(f3, "sys1"))
})
.py_try({
f4 <- Feedback(tf4, tf6)
})
.py_try({
.py_rec("lti:test_Feedback_construction:2217", function() .py_attr(f4, "args"))
})
.py_try({
.py_rec("lti:test_Feedback_construction:2218", function() .py_attr(f4, "sys1"))
})
.py_try({
.py_rec("lti:test_Feedback_construction:2219", function() .py_attr(f4, "var"))
})
.py_try({
f5 <- Feedback(tf5, TransferFunction(1L, 1L, s))
})
.py_try({
.py_rec("lti:test_Feedback_construction:2222", function() .py_attr(f5, "args"))
})
.py_try({
.py_rec("lti:test_Feedback_construction:2223", function() .py_attr(f5, "var"))
})
.py_try({
.py_rec("lti:test_Feedback_construction:2224", function() f5)
})
.py_try({
f6 <- Feedback(TransferFunction(1L, 1L, p), tf4)
})
.py_try({
.py_rec("lti:test_Feedback_construction:2227", function() .py_attr(f6, "args"))
})
.py_try({
.py_rec("lti:test_Feedback_construction:2228", function() .py_attr(f6, "var"))
})
.py_try({
f7 <- .py_neg(Feedback(.py_mul(tf4, tf6), TransferFunction(1L, 1L, p)))
})
.py_try({
.py_rec("lti:test_Feedback_construction:2231", function() .py_attr(f7, "args"))
})
.py_try({
.py_rec("lti:test_Feedback_construction:2232", function() .py_attr(f7, "sys1"))
})
.py_try({
.py_rec("lti:test_Feedback_construction:2235", function() Feedback(tf1, .py_add(tf2, tf3)))
})
.py_try({
.py_rec("lti:test_Feedback_construction:2236", function() Feedback(tf1, .py_Matrix(list(1L, 2L, 3L))))
})
.py_try({
.py_rec("lti:test_Feedback_construction:2237", function() Feedback(TransferFunction(1L, 1L, s), .py_sub(s, 1L)))
})
.py_try({
.py_rec("lti:test_Feedback_construction:2238", function() Feedback(1L, 1L))
})
.py_try({
.py_rec("lti:test_Feedback_construction:2240", function() Feedback(tf2, .py_mul(tf4, tf5)))
})
.py_try({
.py_rec("lti:test_Feedback_construction:2241", function() Feedback(tf2, tf1, .py_flt(1.5)))
})
.py_try({
.py_rec("lti:test_Feedback_construction:2242", function() Feedback(tf1, .py_neg(.py_pow(tf1, -1L))))
})
.py_try({
.py_rec("lti:test_Feedback_construction:2243", function() Feedback(tf4, tf5))
})
.py_try({
df1 <- Feedback(DiscreteTransferFunction(1L, 1L, s, 2L), .py_mul(.py_mul(dtf1, dtf2), dtf3))
})
.py_try({
.py_rec("lti:test_Feedback_construction:2247", function() .py_attr(df1, "args"))
})
.py_try({
.py_rec("lti:test_Feedback_construction:2249", function() .py_attr(df1, "sys1"))
})
.py_try({
.py_rec("lti:test_Feedback_construction:2250", function() .py_attr(df1, "sys2"))
})
.py_try({
.py_rec("lti:test_Feedback_construction:2251", function() .py_attr(df1, "var"))
})
.py_try({
df2 <- Feedback(dtf1, .py_mul(dtf2, dtf3))
})
.py_try({
.py_rec("lti:test_Feedback_construction:2254", function() .py_attr(df2, "args"))
})
.py_try({
.py_rec("lti:test_Feedback_construction:2255", function() .py_attr(df2, "sys1"))
})
.py_try({
.py_rec("lti:test_Feedback_construction:2256", function() .py_attr(df2, "sys2"))
})
.py_try({
.py_rec("lti:test_Feedback_construction:2257", function() .py_attr(df2, "var"))
})
.py_try({
.py_rec("lti:test_Feedback_construction:2258", function() .py_attr(df2, "sampling_time"))
})
.py_try({
df3 <- Feedback(.py_mul(dtf1, dtf2), dtf5)
})
.py_try({
.py_rec("lti:test_Feedback_construction:2261", function() .py_attr(df3, "args"))
})
.py_try({
.py_rec("lti:test_Feedback_construction:2262", function() .py_attr(df3, "sys1"))
})
.py_try({
.py_rec("lti:test_Feedback_construction:2263", function() .py_attr(df3, "sampling_time"))
})
.py_try({
df4 <- Feedback(dtf4, dtf6)
})
.py_try({
.py_rec("lti:test_Feedback_construction:2266", function() .py_attr(df4, "args"))
})
.py_try({
.py_rec("lti:test_Feedback_construction:2267", function() .py_attr(df4, "sys1"))
})
.py_try({
.py_rec("lti:test_Feedback_construction:2268", function() .py_attr(df4, "var"))
})
.py_try({
.py_rec("lti:test_Feedback_construction:2269", function() .py_attr(df4, "sampling_time"))
})
.py_try({
df5 <- Feedback(dtf5, DiscreteTransferFunction(1L, 1L, s, 2L))
})
.py_try({
.py_rec("lti:test_Feedback_construction:2272", function() .py_attr(df5, "args"))
})
.py_try({
.py_rec("lti:test_Feedback_construction:2273", function() .py_attr(df5, "var"))
})
.py_try({
.py_rec("lti:test_Feedback_construction:2274", function() df5)
})
.py_try({
.py_rec("lti:test_Feedback_construction:2275", function() .py_attr(df5, "sampling_time"))
})
.py_try({
df6 <- Feedback(DiscreteTransferFunction(1L, 1L, p, 2L), dtf4)
})
.py_try({
.py_rec("lti:test_Feedback_construction:2278", function() .py_attr(df6, "args"))
})
.py_try({
.py_rec("lti:test_Feedback_construction:2279", function() .py_attr(df6, "var"))
})
.py_try({
.py_rec("lti:test_Feedback_construction:2280", function() .py_attr(df6, "sampling_time"))
})
.py_try({
df7 <- .py_neg(Feedback(.py_mul(dtf4, dtf6), DiscreteTransferFunction(1L, 1L, p, 2L)))
})
.py_try({
.py_rec("lti:test_Feedback_construction:2283", function() .py_attr(df7, "args"))
})
.py_try({
.py_rec("lti:test_Feedback_construction:2286", function() .py_attr(df7, "sys1"))
})
.py_try({
.py_rec("lti:test_Feedback_construction:2288", function() .py_attr(df7, "sampling_time"))
})
.py_try({
.py_rec("lti:test_Feedback_construction:2291", function() Feedback(dtf1, .py_add(dtf2, dtf3)))
})
.py_try({
.py_rec("lti:test_Feedback_construction:2292", function() Feedback(dtf1, .py_Matrix(list(1L, 2L, 3L))))
})
.py_try({
.py_rec("lti:test_Feedback_construction:2293", function() Feedback(DiscreteTransferFunction(1L, 1L, s), .py_sub(s, 1L)))
})
.py_try({
.py_rec("lti:test_Feedback_construction:2294", function() Feedback(1L, 1L))
})
.py_try({
dtf_ <- DiscreteTransferFunction(s, .py_add(s, 1L), s, .py_flt(0.1))
})
.py_try({
.py_rec("lti:test_Feedback_construction:2296", function() Feedback(dtf_, dtf1))
})
.py_try({
.py_rec("lti:test_Feedback_construction:2297", function() Feedback(f1, df1))
})
.py_try({
.py_rec("lti:test_Feedback_construction:2298", function() Feedback(dtf2, .py_mul(dtf4, dtf5)))
})
.py_try({
.py_rec("lti:test_Feedback_construction:2299", function() Feedback(dtf2, dtf1, .py_flt(1.5)))
})
.py_try({
.py_rec("lti:test_Feedback_construction:2300", function() Feedback(dtf1, .py_neg(.py_pow(dtf1, -1L))))
})
.py_try({
.py_rec("lti:test_Feedback_construction:2301", function() Feedback(dtf4, dtf5))
}))---")
  expect_identical(got[["lti:test_Feedback_construction:2196"]], "True")
  expect_identical(got[["lti:test_Feedback_construction:2197"]], "False")
  expect_identical(got[["lti:test_Feedback_construction:2201"]], "(TransferFunction(1, 1, s), Series(TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s), TransferFunction(k, 1, s), TransferFunction(a2*p - s, a2*s + p, s)), -1)")
  expect_identical(got[["lti:test_Feedback_construction:2202"]], "TransferFunction(1, 1, s)")
  expect_identical(got[["lti:test_Feedback_construction:2203"]], "Series(TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s), TransferFunction(k, 1, s), TransferFunction(a2*p - s, a2*s + p, s))")
  expect_identical(got[["lti:test_Feedback_construction:2204"]], "s")
  expect_identical(got[["lti:test_Feedback_construction:2207"]], "(TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s), Series(TransferFunction(k, 1, s), TransferFunction(a2*p - s, a2*s + p, s)), -1)")
  expect_identical(got[["lti:test_Feedback_construction:2208"]], "TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s)")
  expect_identical(got[["lti:test_Feedback_construction:2209"]], "Series(TransferFunction(k, 1, s), TransferFunction(a2*p - s, a2*s + p, s))")
  expect_identical(got[["lti:test_Feedback_construction:2210"]], "s")
  expect_identical(got[["lti:test_Feedback_construction:2213"]], "(Series(TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s), TransferFunction(k, 1, s)), TransferFunction(-a0 + a1*s^2 + a2*s, a0 + s, s), -1)")
  expect_identical(got[["lti:test_Feedback_construction:2214"]], "Series(TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s), TransferFunction(k, 1, s))")
  expect_identical(got[["lti:test_Feedback_construction:2217"]], "(TransferFunction(a0*p + p^a1 - s, p, p), TransferFunction(-p + s, p + s, p), -1)")
  expect_identical(got[["lti:test_Feedback_construction:2218"]], "TransferFunction(a0*p + p^a1 - s, p, p)")
  expect_identical(got[["lti:test_Feedback_construction:2219"]], "p")
  expect_identical(got[["lti:test_Feedback_construction:2222"]], "(TransferFunction(-a0 + a1*s^2 + a2*s, a0 + s, s), TransferFunction(1, 1, s), -1)")
  expect_identical(got[["lti:test_Feedback_construction:2223"]], "s")
  expect_identical(got[["lti:test_Feedback_construction:2224"]], "Feedback(TransferFunction(-a0 + a1*s^2 + a2*s, a0 + s, s), TransferFunction(1, 1, s), -1)")
  expect_identical(got[["lti:test_Feedback_construction:2227"]], "(TransferFunction(1, 1, p), TransferFunction(a0*p + p^a1 - s, p, p), -1)")
  expect_identical(got[["lti:test_Feedback_construction:2228"]], "p")
  expect_identical(got[["lti:test_Feedback_construction:2231"]], "(Series(TransferFunction(-1, 1, p), TransferFunction(a0*p + p^a1 - s, p, p), TransferFunction(-p + s, p + s, p)), TransferFunction(-1, 1, p), -1)")
  expect_identical(got[["lti:test_Feedback_construction:2232"]], "Series(TransferFunction(-1, 1, p), TransferFunction(a0*p + p^a1 - s, p, p), TransferFunction(-p + s, p + s, p))")
  expect_identical(got[["lti:test_Feedback_construction:2235"]], "ERR TypeError: Unsupported type for `sys2` in Feedback.")
  expect_identical(got[["lti:test_Feedback_construction:2236"]], "ERR TypeError: Unsupported type for `sys2` in Feedback.")
  expect_identical(got[["lti:test_Feedback_construction:2237"]], "ERR TypeError: Unsupported type for `sys2` in Feedback.")
  expect_identical(got[["lti:test_Feedback_construction:2238"]], "ERR TypeError: Unsupported type for `sys1` in Feedback.")
  expect_identical(got[["lti:test_Feedback_construction:2240"]], "ERR ValueError: \\nAll the transfer functions should use the same complex variable of the\\nLaplace domain or z-domain.")
  expect_identical(got[["lti:test_Feedback_construction:2241"]], "ERR ValueError: \\nUnsupported type for feedback. `sign` arg should either be 1 (positive\\nfeedback loop) or -1 (negative feedback loop).")
  expect_identical(got[["lti:test_Feedback_construction:2242"]], "ERR ValueError: \\nThe equivalent system will have zero\\ndenominator.")
  expect_identical(got[["lti:test_Feedback_construction:2243"]], "ERR ValueError: \\nBoth `sys1` and `sys2` should be                     using the same\\ncomplex variable.")
  expect_identical(got[["lti:test_Feedback_construction:2247"]], "(DiscreteTransferFunction(1, 1, s, 2), Series(DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 2), DiscreteTransferFunction(k, 1, s, 2), DiscreteTransferFunction(a2*p - s, a2*s + p, s, 2)), -1)")
  expect_identical(got[["lti:test_Feedback_construction:2249"]], "DiscreteTransferFunction(1, 1, s, 2)")
  expect_identical(got[["lti:test_Feedback_construction:2250"]], "Series(DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 2), DiscreteTransferFunction(k, 1, s, 2), DiscreteTransferFunction(a2*p - s, a2*s + p, s, 2))")
  expect_identical(got[["lti:test_Feedback_construction:2251"]], "s")
  expect_identical(got[["lti:test_Feedback_construction:2254"]], "(DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 2), Series(DiscreteTransferFunction(k, 1, s, 2), DiscreteTransferFunction(a2*p - s, a2*s + p, s, 2)), -1)")
  expect_identical(got[["lti:test_Feedback_construction:2255"]], "DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 2)")
  expect_identical(got[["lti:test_Feedback_construction:2256"]], "Series(DiscreteTransferFunction(k, 1, s, 2), DiscreteTransferFunction(a2*p - s, a2*s + p, s, 2))")
  expect_identical(got[["lti:test_Feedback_construction:2257"]], "s")
  expect_identical(got[["lti:test_Feedback_construction:2258"]], "2")
  expect_identical(got[["lti:test_Feedback_construction:2261"]], "(Series(DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 2), DiscreteTransferFunction(k, 1, s, 2)), DiscreteTransferFunction(-a0 + a1*s^2 + a2*s, a0 + s, s, 2), -1)")
  expect_identical(got[["lti:test_Feedback_construction:2262"]], "Series(DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 2), DiscreteTransferFunction(k, 1, s, 2))")
  expect_identical(got[["lti:test_Feedback_construction:2263"]], "2")
  expect_identical(got[["lti:test_Feedback_construction:2266"]], "(DiscreteTransferFunction(a0*p + p^a1 - s, p, p, 2), DiscreteTransferFunction(-p + s, p + s, p, 2), -1)")
  expect_identical(got[["lti:test_Feedback_construction:2267"]], "DiscreteTransferFunction(a0*p + p^a1 - s, p, p, 2)")
  expect_identical(got[["lti:test_Feedback_construction:2268"]], "p")
  expect_identical(got[["lti:test_Feedback_construction:2269"]], "2")
  expect_identical(got[["lti:test_Feedback_construction:2272"]], "(DiscreteTransferFunction(-a0 + a1*s^2 + a2*s, a0 + s, s, 2), DiscreteTransferFunction(1, 1, s, 2), -1)")
  expect_identical(got[["lti:test_Feedback_construction:2273"]], "s")
  expect_identical(got[["lti:test_Feedback_construction:2274"]], "Feedback(DiscreteTransferFunction(-a0 + a1*s^2 + a2*s, a0 + s, s, 2), DiscreteTransferFunction(1, 1, s, 2), -1)")
  expect_identical(got[["lti:test_Feedback_construction:2275"]], "2")
  expect_identical(got[["lti:test_Feedback_construction:2278"]], "(DiscreteTransferFunction(1, 1, p, 2), DiscreteTransferFunction(a0*p + p^a1 - s, p, p, 2), -1)")
  expect_identical(got[["lti:test_Feedback_construction:2279"]], "p")
  expect_identical(got[["lti:test_Feedback_construction:2280"]], "2")
  expect_identical(got[["lti:test_Feedback_construction:2283"]], "(Series(DiscreteTransferFunction(-1, 1, p, 2), DiscreteTransferFunction(a0*p + p^a1 - s, p, p, 2), DiscreteTransferFunction(-p + s, p + s, p, 2)), DiscreteTransferFunction(-1, 1, p, 2), -1)")
  expect_identical(got[["lti:test_Feedback_construction:2286"]], "Series(DiscreteTransferFunction(-1, 1, p, 2), DiscreteTransferFunction(a0*p + p^a1 - s, p, p, 2), DiscreteTransferFunction(-p + s, p + s, p, 2))")
  expect_identical(got[["lti:test_Feedback_construction:2288"]], "2")
  expect_identical(got[["lti:test_Feedback_construction:2291"]], "ERR TypeError: Unsupported type for `sys2` in Feedback.")
  expect_identical(got[["lti:test_Feedback_construction:2292"]], "ERR TypeError: Unsupported type for `sys2` in Feedback.")
  expect_identical(got[["lti:test_Feedback_construction:2293"]], "ERR TypeError: Unsupported type for `sys2` in Feedback.")
  expect_identical(got[["lti:test_Feedback_construction:2294"]], "ERR TypeError: Unsupported type for `sys1` in Feedback.")
  expect_identical(got[["lti:test_Feedback_construction:2296"]], "ERR TypeError: \\n                Incompatible sampling times detected.\\n                All systems should have the same sampling time.\\n                Found sampling time: 0.100000000000000 and 2")
  expect_identical(got[["lti:test_Feedback_construction:2297"]], "ERR TypeError: \\n                Incompatible systems detected.\\n                All systems should be either continuous-time or discrete-time.\\n                Found: Feedback(TransferFunction(1, 1, s), Series(TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s), TransferFunction(k, 1, s), TransferFunction(a2*p - s, a2*s + p, s)), -1) and Feedback(DiscreteTransferFunction(1, 1, s, 2), Series(DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 2), DiscreteTransferFunction(k, 1, s, 2), DiscreteTransferFunction(a2*p - s, a2*s + p, s, 2)), -1)")
  expect_identical(got[["lti:test_Feedback_construction:2298"]], "ERR ValueError: \\nAll the transfer functions should use the same complex variable of the\\nLaplace domain or z-domain.")
  expect_identical(got[["lti:test_Feedback_construction:2299"]], "ERR ValueError: \\nUnsupported type for feedback. `sign` arg should either be 1 (positive\\nfeedback loop) or -1 (negative feedback loop).")
  expect_identical(got[["lti:test_Feedback_construction:2300"]], "ERR ValueError: \\nThe equivalent system will have zero\\ndenominator.")
  expect_identical(got[["lti:test_Feedback_construction:2301"]], "ERR ValueError: \\nBoth `sys1` and `sys2` should be                     using the same\\ncomplex variable.")
})

test_that("lti test_Feedback_functions matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("a, x, b, c, s, g, d,    p, k, tau, zeta, wn, T, z"); a <- .py_getitem(.tmpu_, 0L); x <- .py_getitem(.tmpu_, 1L); b <- .py_getitem(.tmpu_, 2L); c <- .py_getitem(.tmpu_, 3L); s <- .py_getitem(.tmpu_, 4L); g <- .py_getitem(.tmpu_, 5L); d <- .py_getitem(.tmpu_, 6L); p <- .py_getitem(.tmpu_, 7L); k <- .py_getitem(.tmpu_, 8L); tau <- .py_getitem(.tmpu_, 9L); zeta <- .py_getitem(.tmpu_, 10L); wn <- .py_getitem(.tmpu_, 11L); T <- .py_getitem(.tmpu_, 12L); z <- .py_getitem(.tmpu_, 13L)
})
.py_try({
.tmpu_ <- .py_symbols("a0:4, b0:5, c0:4, d0:4"); a0 <- .py_getitem(.tmpu_, 0L); a1 <- .py_getitem(.tmpu_, 1L); a2 <- .py_getitem(.tmpu_, 2L); a3 <- .py_getitem(.tmpu_, 3L); b0 <- .py_getitem(.tmpu_, 4L); b1 <- .py_getitem(.tmpu_, 5L); b2 <- .py_getitem(.tmpu_, 6L); b3 <- .py_getitem(.tmpu_, 7L); b4 <- .py_getitem(.tmpu_, 8L); c0 <- .py_getitem(.tmpu_, 9L); c1 <- .py_getitem(.tmpu_, 10L); c2 <- .py_getitem(.tmpu_, 11L); c3 <- .py_getitem(.tmpu_, 12L); d0 <- .py_getitem(.tmpu_, 13L); d1 <- .py_getitem(.tmpu_, 14L); d2 <- .py_getitem(.tmpu_, 15L); d3 <- .py_getitem(.tmpu_, 16L)
})
.py_try({
TF1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
TF2 <- TransferFunction(k, 1L, s)
})
.py_try({
TF3 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
TF4 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(tau, s), 1L), s)
})
.py_try({
TF5 <- TransferFunction(k, p, s)
})
.py_try({
TF6 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(k, s), p), s)
})
p_are_floats_equals <- function(n1, n2) {
  "\n    Helper function to compare two floats.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return(.py_method(Float(n1), "is_same", Float(n2), .py_isclose))
}
p_are_complex_equals <- function(n1, n2) {
  "\n    Helper function to compare two complex numbers.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return((.py_truthy(p_are_floats_equals(re(n1), re(n2))) && .py_truthy(p_are_floats_equals(im(n1), im(n2)))))
}
.py_try({
tf <- TransferFunction(1L, 1L, s)
})
.py_try({
tf1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
tf2 <- TransferFunction(k, 1L, s)
})
.py_try({
tf3 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
tf4 <- TransferFunction(.py_sub(.py_add(.py_mul(a0, p), .py_pow(p, a1)), s), p, p)
})
.py_try({
tf5 <- TransferFunction(.py_sub(.py_add(.py_mul(a1, .py_pow(s, 2L)), .py_mul(a2, s)), a0), .py_add(s, a0), s)
})
.py_try({
tf6 <- TransferFunction(.py_sub(s, p), .py_add(p, s), p)
})
.py_try({
dtf <- DiscreteTransferFunction(1L, 1L, s, 10L)
})
.py_try({
dtf1 <- DiscreteTransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s, 10L)
})
.py_try({
dtf2 <- DiscreteTransferFunction(k, 1L, s, 10L)
})
.py_try({
dtf3 <- DiscreteTransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s, 10L)
})
.py_try({
dtf4 <- DiscreteTransferFunction(.py_sub(.py_add(.py_mul(a0, p), .py_pow(p, a1)), s), p, p, 10L)
})
.py_try({
dtf5 <- DiscreteTransferFunction(.py_sub(.py_add(.py_mul(a1, .py_pow(s, 2L)), .py_mul(a2, s)), a0), .py_add(s, a0), s, 10L)
})
.py_try({
dtf6 <- DiscreteTransferFunction(.py_sub(s, p), .py_add(p, s), p, 10L)
})
.py_try({
.py_rec("lti:test_Feedback_functions:2322", function() .py_mul(.py_div(.py_mul(.py_mul(tf1, tf2), tf3), tf3), tf5))
})
.py_try({
.py_rec("lti:test_Feedback_functions:2323", function() .py_div(.py_mul(.py_mul(tf1, tf2), tf3), .py_mul(tf3, tf5)))
})
.py_try({
.py_rec("lti:test_Feedback_functions:2324", function() .py_div(tf, .py_add(tf, tf1)))
})
.py_try({
.py_rec("lti:test_Feedback_functions:2325", function() .py_div(tf, .py_add(tf, .py_mul(.py_mul(tf1, tf2), tf3))))
})
.py_try({
.py_rec("lti:test_Feedback_functions:2326", function() .py_div(tf1, .py_add(tf, .py_mul(.py_mul(tf1, tf2), tf3))))
})
.py_try({
.py_rec("lti:test_Feedback_functions:2327", function() .py_div(.py_mul(tf1, tf2), .py_add(tf, .py_mul(tf1, tf2))))
})
.py_try({
.py_rec("lti:test_Feedback_functions:2328", function() .py_div(.py_mul(tf1, tf2), .py_add(tf, .py_mul(.py_mul(tf1, tf2), tf5))))
})
.py_try({
.py_rec("lti:test_Feedback_functions:2329", function() .py_in(.py_div(.py_mul(tf1, tf2), .py_add(tf, .py_mul(.py_mul(.py_mul(tf1, tf2), tf5), tf3))), .py_tup(Feedback(.py_mul(tf1, tf2), .py_mul(tf5, tf3)), Feedback(.py_mul(tf1, tf2), .py_mul(tf3, tf5)))))
})
.py_try({
.py_rec("lti:test_Feedback_functions:2330", function() .py_div(tf4, .py_add(TransferFunction(1L, 1L, p), .py_mul(tf4, tf6))))
})
.py_try({
.py_rec("lti:test_Feedback_functions:2331", function() .py_div(tf5, .py_add(tf, tf5)))
})
.py_try({
.py_rec("lti:test_Feedback_functions:2333", function() .py_div(.py_mul(.py_mul(tf1, tf2), tf3), .py_add(1L, .py_mul(.py_mul(tf1, tf2), tf3))))
})
.py_try({
.py_rec("lti:test_Feedback_functions:2334", function() .py_div(.py_mul(tf2, tf3), .py_add(tf, .py_mul(.py_mul(tf2, tf3), tf4))))
})
.py_try({
.py_rec("lti:test_Feedback_functions:2336", function() .py_method(Feedback(tf, .py_mul(.py_mul(tf1, tf2), tf3)), "doit"))
})
.py_try({
.py_rec("lti:test_Feedback_functions:2339", function() .py_attr(Feedback(tf, .py_mul(.py_mul(tf1, tf2), tf3)), "sensitivity"))
})
.py_try({
.py_rec("lti:test_Feedback_functions:2341", function() .py_method(Feedback(tf1, .py_mul(tf2, tf3)), "doit"))
})
.py_try({
.py_rec("lti:test_Feedback_functions:2344", function() .py_attr(Feedback(tf1, .py_mul(tf2, tf3)), "sensitivity"))
})
.py_try({
.py_rec("lti:test_Feedback_functions:2346", function() .py_method(Feedback(.py_mul(tf1, tf2), tf5), "doit"))
})
.py_try({
.py_rec("lti:test_Feedback_functions:2349", function() .py_attr(Feedback(.py_mul(tf1, tf2), tf5, 1L), "sensitivity"))
})
.py_try({
.py_rec("lti:test_Feedback_functions:2351", function() .py_method(Feedback(tf4, tf6), "doit"))
})
.py_try({
.py_rec("lti:test_Feedback_functions:2353", function() .py_neg(.py_method(Feedback(.py_mul(tf4, tf6), TransferFunction(1L, 1L, p)), "doit")))
})
.py_try({
.py_rec("lti:test_Feedback_functions:2355", function() .py_method(Feedback(tf, tf), "doit"))
})
.py_try({
.py_rec("lti:test_Feedback_functions:2357", function() .py_method(Feedback(tf1, .py_mul(tf2, tf5)), "rewrite", TransferFunction))
})
.py_try({
.py_rec("lti:test_Feedback_functions:2362", function() .py_method(Feedback(TransferFunction(1L, 1L, p), tf4), "rewrite", TransferFunction))
})
.py_try({
.py_rec("lti:test_Feedback_functions:2367", function() .py_mul(.py_div(.py_mul(.py_mul(dtf1, dtf2), dtf3), dtf3), dtf5))
})
.py_try({
.py_rec("lti:test_Feedback_functions:2369", function() .py_div(.py_mul(.py_mul(dtf1, dtf2), dtf3), .py_mul(dtf3, dtf5)))
})
.py_try({
.py_rec("lti:test_Feedback_functions:2371", function() .py_div(dtf, .py_add(dtf, dtf1)))
})
.py_try({
.py_rec("lti:test_Feedback_functions:2372", function() .py_div(dtf, .py_add(dtf, .py_mul(.py_mul(dtf1, dtf2), dtf3))))
})
.py_try({
.py_rec("lti:test_Feedback_functions:2373", function() .py_div(dtf1, .py_add(dtf, .py_mul(.py_mul(dtf1, dtf2), dtf3))))
})
.py_try({
.py_rec("lti:test_Feedback_functions:2374", function() .py_div(.py_mul(dtf1, dtf2), .py_add(dtf, .py_mul(dtf1, dtf2))))
})
.py_try({
.py_rec("lti:test_Feedback_functions:2375", function() .py_div(.py_mul(dtf1, dtf2), .py_add(dtf, .py_mul(.py_mul(dtf1, dtf2), dtf5))))
})
.py_try({
.py_rec("lti:test_Feedback_functions:2376", function() .py_in(.py_div(.py_mul(dtf1, dtf2), .py_add(dtf, .py_mul(.py_mul(.py_mul(dtf1, dtf2), dtf5), dtf3))), .py_tup(Feedback(.py_mul(dtf1, dtf2), .py_mul(dtf5, dtf3)), Feedback(.py_mul(dtf1, dtf2), .py_mul(dtf3, dtf5)))))
})
.py_try({
.py_rec("lti:test_Feedback_functions:2378", function() .py_div(dtf4, .py_add(DiscreteTransferFunction(1L, 1L, p, 10L), .py_mul(dtf4, dtf6))))
})
.py_try({
.py_rec("lti:test_Feedback_functions:2380", function() .py_div(dtf5, .py_add(dtf, dtf5)))
})
.py_try({
.py_rec("lti:test_Feedback_functions:2382", function() .py_div(.py_mul(.py_mul(dtf1, dtf2), dtf3), .py_add(1L, .py_mul(.py_mul(dtf1, dtf2), dtf3))))
})
.py_try({
.py_rec("lti:test_Feedback_functions:2383", function() .py_div(.py_mul(dtf2, dtf3), .py_add(dtf, .py_mul(.py_mul(dtf2, dtf3), dtf4))))
})
.py_try({
.py_rec("lti:test_Feedback_functions:2385", function() .py_method(Feedback(dtf, .py_mul(.py_mul(dtf1, dtf2), dtf3)), "doit"))
})
.py_try({
.py_rec("lti:test_Feedback_functions:2389", function() .py_attr(Feedback(dtf, .py_mul(.py_mul(dtf1, dtf2), dtf3)), "sensitivity"))
})
.py_try({
.py_rec("lti:test_Feedback_functions:2391", function() .py_method(Feedback(dtf1, .py_mul(dtf2, dtf3)), "doit"))
})
.py_try({
.py_rec("lti:test_Feedback_functions:2396", function() .py_attr(Feedback(dtf1, .py_mul(dtf2, dtf3)), "sensitivity"))
})
.py_try({
.py_rec("lti:test_Feedback_functions:2398", function() .py_method(Feedback(.py_mul(dtf1, dtf2), dtf5), "doit"))
})
.py_try({
.py_rec("lti:test_Feedback_functions:2403", function() .py_attr(Feedback(.py_mul(dtf1, dtf2), dtf5, 1L), "sensitivity"))
})
.py_try({
.py_rec("lti:test_Feedback_functions:2406", function() .py_method(Feedback(dtf4, dtf6), "doit"))
})
.py_try({
.py_rec("lti:test_Feedback_functions:2409", function() .py_neg(.py_method(Feedback(.py_mul(dtf4, dtf6), DiscreteTransferFunction(1L, 1L, p, 10L)), "doit")))
})
.py_try({
.py_rec("lti:test_Feedback_functions:2413", function() .py_method(Feedback(dtf, dtf), "doit"))
})
.py_try({
.py_rec("lti:test_Feedback_functions:2415", function() .py_method(Feedback(dtf1, .py_mul(dtf2, dtf5)), "rewrite", DiscreteTransferFunction))
})
.py_try({
.py_rec("lti:test_Feedback_functions:2420", function() .py_method(Feedback(DiscreteTransferFunction(1L, 1L, p, 10L), dtf4), "rewrite", DiscreteTransferFunction))
})
.py_try({
.py_rec("lti:test_Feedback_functions:2425", function() .py_method(Feedback(dtf1, dtf2), "rewrite", TransferFunction))
})
.py_try({
.py_rec("lti:test_Feedback_functions:2427", function() .py_method(Feedback(tf1, tf2), "rewrite", DiscreteTransferFunction))
}))---")
  expect_identical(got[["lti:test_Feedback_functions:2322"]], "Series(TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s), TransferFunction(k, 1, s), TransferFunction(a2*p - s, a2*s + p, s), TransferFunction(a2*s + p, a2*p - s, s), TransferFunction(-a0 + a1*s^2 + a2*s, a0 + s, s))")
  expect_identical(got[["lti:test_Feedback_functions:2323"]], "Series(TransferFunction(k*(a2*p - s), (a2*s + p)*(s^2 + 2*s*wn*zeta + wn^2), s), TransferFunction((a0 + s)*(a2*s + p), (a2*p - s)*(-a0 + a1*s^2 + a2*s), s))")
  expect_identical(got[["lti:test_Feedback_functions:2324"]], "Feedback(TransferFunction(1, 1, s), TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s), -1)")
  expect_identical(got[["lti:test_Feedback_functions:2325"]], "Feedback(TransferFunction(1, 1, s), Series(TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s), TransferFunction(k, 1, s), TransferFunction(a2*p - s, a2*s + p, s)), -1)")
  expect_identical(got[["lti:test_Feedback_functions:2326"]], "Feedback(TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s), Series(TransferFunction(k, 1, s), TransferFunction(a2*p - s, a2*s + p, s)), -1)")
  expect_identical(got[["lti:test_Feedback_functions:2327"]], "Feedback(Series(TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s), TransferFunction(k, 1, s)), TransferFunction(1, 1, s), -1)")
  expect_identical(got[["lti:test_Feedback_functions:2328"]], "Feedback(Series(TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s), TransferFunction(k, 1, s)), TransferFunction(-a0 + a1*s^2 + a2*s, a0 + s, s), -1)")
  expect_identical(got[["lti:test_Feedback_functions:2329"]], "True")
  expect_identical(got[["lti:test_Feedback_functions:2330"]], "Feedback(TransferFunction(a0*p + p^a1 - s, p, p), TransferFunction(-p + s, p + s, p), -1)")
  expect_identical(got[["lti:test_Feedback_functions:2331"]], "Feedback(TransferFunction(-a0 + a1*s^2 + a2*s, a0 + s, s), TransferFunction(1, 1, s), -1)")
  expect_identical(got[["lti:test_Feedback_functions:2333"]], "ERR TypeError: unsupported operand type(s) for +: 'int' and 'Series'")
  expect_identical(got[["lti:test_Feedback_functions:2334"]], "ERR ValueError: \\nAll transfer functions should use the same complex variable of the\\nLaplace transform or z-transform. 2 different values found.")
  expect_identical(got[["lti:test_Feedback_functions:2336"]], "TransferFunction((a2*s + p)*(s^2 + 2*s*wn*zeta + wn^2), k*(a2*p - s) + (a2*s + p)*(s^2 + 2*s*wn*zeta + wn^2), s)")
  expect_identical(got[["lti:test_Feedback_functions:2339"]], "1/(k*(a2*p - s)/((a2*s + p)*(s^2 + 2*s*wn*zeta + wn^2)) + 1)")
  expect_identical(got[["lti:test_Feedback_functions:2341"]], "TransferFunction((a2*s + p)*(s^2 + 2*s*wn*zeta + wn^2), (k*(a2*p - s) + (a2*s + p)*(s^2 + 2*s*wn*zeta + wn^2))*(s^2 + 2*s*wn*zeta + wn^2), s)")
  expect_identical(got[["lti:test_Feedback_functions:2344"]], "1/(k*(a2*p - s)/((a2*s + p)*(s^2 + 2*s*wn*zeta + wn^2)) + 1)")
  expect_identical(got[["lti:test_Feedback_functions:2346"]], "TransferFunction(k*(a0 + s)*(s^2 + 2*s*wn*zeta + wn^2), (k*(-a0 + a1*s^2 + a2*s) + (a0 + s)*(s^2 + 2*s*wn*zeta + wn^2))*(s^2 + 2*s*wn*zeta + wn^2), s)")
  expect_identical(got[["lti:test_Feedback_functions:2349"]], "1/(-k*(-a0 + a1*s^2 + a2*s)/((a0 + s)*(s^2 + 2*s*wn*zeta + wn^2)) + 1)")
  expect_identical(got[["lti:test_Feedback_functions:2351"]], "TransferFunction(p*(p + s)*(a0*p + p^a1 - s), p*(p*(p + s) + (-p + s)*(a0*p + p^a1 - s)), p)")
  expect_identical(got[["lti:test_Feedback_functions:2353"]], "TransferFunction(-p*(-p + s)*(p + s)*(a0*p + p^a1 - s), p*(p + s)*(p*(p + s) + (-p + s)*(a0*p + p^a1 - s)), p)")
  expect_identical(got[["lti:test_Feedback_functions:2355"]], "TransferFunction(1, 2, s)")
  expect_identical(got[["lti:test_Feedback_functions:2357"]], "TransferFunction((a0 + s)*(s^2 + 2*s*wn*zeta + wn^2), (k*(-a0 + a1*s^2 + a2*s) + (a0 + s)*(s^2 + 2*s*wn*zeta + wn^2))*(s^2 + 2*s*wn*zeta + wn^2), s)")
  expect_identical(got[["lti:test_Feedback_functions:2362"]], "TransferFunction(p, a0*p + p + p^a1 - s, p)")
  expect_identical(got[["lti:test_Feedback_functions:2367"]], "Series(DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 10), DiscreteTransferFunction(k, 1, s, 10), DiscreteTransferFunction(a2*p - s, a2*s + p, s, 10), DiscreteTransferFunction(a2*s + p, a2*p - s, s, 10), DiscreteTransferFunction(-a0 + a1*s^2 + a2*s, a0 + s, s, 10))")
  expect_identical(got[["lti:test_Feedback_functions:2369"]], "Series(DiscreteTransferFunction(k*(a2*p - s), (a2*s + p)*(s^2 + 2*s*wn*zeta + wn^2), s, 10), DiscreteTransferFunction((a0 + s)*(a2*s + p), (a2*p - s)*(-a0 + a1*s^2 + a2*s), s, 10))")
  expect_identical(got[["lti:test_Feedback_functions:2371"]], "Feedback(DiscreteTransferFunction(1, 1, s, 10), DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 10), -1)")
  expect_identical(got[["lti:test_Feedback_functions:2372"]], "Feedback(DiscreteTransferFunction(1, 1, s, 10), Series(DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 10), DiscreteTransferFunction(k, 1, s, 10), DiscreteTransferFunction(a2*p - s, a2*s + p, s, 10)), -1)")
  expect_identical(got[["lti:test_Feedback_functions:2373"]], "Feedback(DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 10), Series(DiscreteTransferFunction(k, 1, s, 10), DiscreteTransferFunction(a2*p - s, a2*s + p, s, 10)), -1)")
  expect_identical(got[["lti:test_Feedback_functions:2374"]], "Feedback(Series(DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 10), DiscreteTransferFunction(k, 1, s, 10)), DiscreteTransferFunction(1, 1, s, 10), -1)")
  expect_identical(got[["lti:test_Feedback_functions:2375"]], "Feedback(Series(DiscreteTransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s, 10), DiscreteTransferFunction(k, 1, s, 10)), DiscreteTransferFunction(-a0 + a1*s^2 + a2*s, a0 + s, s, 10), -1)")
  expect_identical(got[["lti:test_Feedback_functions:2376"]], "True")
  expect_identical(got[["lti:test_Feedback_functions:2378"]], "Feedback(DiscreteTransferFunction(a0*p + p^a1 - s, p, p, 10), DiscreteTransferFunction(-p + s, p + s, p, 10), -1)")
  expect_identical(got[["lti:test_Feedback_functions:2380"]], "Feedback(DiscreteTransferFunction(-a0 + a1*s^2 + a2*s, a0 + s, s, 10), DiscreteTransferFunction(1, 1, s, 10), -1)")
  expect_identical(got[["lti:test_Feedback_functions:2382"]], "ERR TypeError: unsupported operand type(s) for +: 'int' and 'Series'")
  expect_identical(got[["lti:test_Feedback_functions:2383"]], "ERR ValueError: \\nAll transfer functions should use the same complex variable of the\\nLaplace transform or z-transform. 2 different values found.")
  expect_identical(got[["lti:test_Feedback_functions:2385"]], "DiscreteTransferFunction((a2*s + p)*(s^2 + 2*s*wn*zeta + wn^2), k*(a2*p - s) + (a2*s + p)*(s^2 + 2*s*wn*zeta + wn^2), s, 10)")
  expect_identical(got[["lti:test_Feedback_functions:2389"]], "1/(k*(a2*p - s)/((a2*s + p)*(s^2 + 2*s*wn*zeta + wn^2)) + 1)")
  expect_identical(got[["lti:test_Feedback_functions:2391"]], "DiscreteTransferFunction((a2*s + p)*(s^2 + 2*s*wn*zeta + wn^2), (k*(a2*p - s) + (a2*s + p)*(s^2 + 2*s*wn*zeta + wn^2))*(s^2 + 2*s*wn*zeta + wn^2), s, 10)")
  expect_identical(got[["lti:test_Feedback_functions:2396"]], "1/(k*(a2*p - s)/((a2*s + p)*(s^2 + 2*s*wn*zeta + wn^2)) + 1)")
  expect_identical(got[["lti:test_Feedback_functions:2398"]], "DiscreteTransferFunction(k*(a0 + s)*(s^2 + 2*s*wn*zeta + wn^2), (k*(-a0 + a1*s^2 + a2*s) + (a0 + s)*(s^2 + 2*s*wn*zeta + wn^2))*(s^2 + 2*s*wn*zeta + wn^2), s, 10)")
  expect_identical(got[["lti:test_Feedback_functions:2403"]], "1/(-k*(-a0 + a1*s^2 + a2*s)/((a0 + s)*(s^2 + 2*s*wn*zeta + wn^2)) + 1)")
  expect_identical(got[["lti:test_Feedback_functions:2406"]], "DiscreteTransferFunction(p*(p + s)*(a0*p + p^a1 - s), p*(p*(p + s) + (-p + s)*(a0*p + p^a1 - s)), p, 10)")
  expect_identical(got[["lti:test_Feedback_functions:2409"]], "DiscreteTransferFunction(-p*(-p + s)*(p + s)*(a0*p + p^a1 - s), p*(p + s)*(p*(p + s) + (-p + s)*(a0*p + p^a1 - s)), p, 10)")
  expect_identical(got[["lti:test_Feedback_functions:2413"]], "DiscreteTransferFunction(1, 2, s, 10)")
  expect_identical(got[["lti:test_Feedback_functions:2415"]], "DiscreteTransferFunction((a0 + s)*(s^2 + 2*s*wn*zeta + wn^2), (k*(-a0 + a1*s^2 + a2*s) + (a0 + s)*(s^2 + 2*s*wn*zeta + wn^2))*(s^2 + 2*s*wn*zeta + wn^2), s, 10)")
  expect_identical(got[["lti:test_Feedback_functions:2420"]], "DiscreteTransferFunction(p, a0*p + p + p^a1 - s, p, 10)")
  expect_identical(got[["lti:test_Feedback_functions:2425"]], "ERR TypeError: \\n                    Cannot rewrite a discrete-time Feedback object as a\\n                    TransferFunction.")
  expect_identical(got[["lti:test_Feedback_functions:2427"]], "ERR TypeError: \\n                    Cannot rewrite a continuous-time Feedback object as a\\n                    DiscreteTransferFunction.")
})

test_that("lti test_Feedback_with_Series matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("a, x, b, c, s, g, d,    p, k, tau, zeta, wn, T, z"); a <- .py_getitem(.tmpu_, 0L); x <- .py_getitem(.tmpu_, 1L); b <- .py_getitem(.tmpu_, 2L); c <- .py_getitem(.tmpu_, 3L); s <- .py_getitem(.tmpu_, 4L); g <- .py_getitem(.tmpu_, 5L); d <- .py_getitem(.tmpu_, 6L); p <- .py_getitem(.tmpu_, 7L); k <- .py_getitem(.tmpu_, 8L); tau <- .py_getitem(.tmpu_, 9L); zeta <- .py_getitem(.tmpu_, 10L); wn <- .py_getitem(.tmpu_, 11L); T <- .py_getitem(.tmpu_, 12L); z <- .py_getitem(.tmpu_, 13L)
})
.py_try({
.tmpu_ <- .py_symbols("a0:4, b0:5, c0:4, d0:4"); a0 <- .py_getitem(.tmpu_, 0L); a1 <- .py_getitem(.tmpu_, 1L); a2 <- .py_getitem(.tmpu_, 2L); a3 <- .py_getitem(.tmpu_, 3L); b0 <- .py_getitem(.tmpu_, 4L); b1 <- .py_getitem(.tmpu_, 5L); b2 <- .py_getitem(.tmpu_, 6L); b3 <- .py_getitem(.tmpu_, 7L); b4 <- .py_getitem(.tmpu_, 8L); c0 <- .py_getitem(.tmpu_, 9L); c1 <- .py_getitem(.tmpu_, 10L); c2 <- .py_getitem(.tmpu_, 11L); c3 <- .py_getitem(.tmpu_, 12L); d0 <- .py_getitem(.tmpu_, 13L); d1 <- .py_getitem(.tmpu_, 14L); d2 <- .py_getitem(.tmpu_, 15L); d3 <- .py_getitem(.tmpu_, 16L)
})
.py_try({
TF1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
TF2 <- TransferFunction(k, 1L, s)
})
.py_try({
TF3 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
TF4 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(tau, s), 1L), s)
})
.py_try({
TF5 <- TransferFunction(k, p, s)
})
.py_try({
TF6 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(k, s), p), s)
})
p_are_floats_equals <- function(n1, n2) {
  "\n    Helper function to compare two floats.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return(.py_method(Float(n1), "is_same", Float(n2), .py_isclose))
}
p_are_complex_equals <- function(n1, n2) {
  "\n    Helper function to compare two complex numbers.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return((.py_truthy(p_are_floats_equals(re(n1), re(n2))) && .py_truthy(p_are_floats_equals(im(n1), im(n2)))))
}
.py_try({
tf1 <- TransferFunction(.py_add(s, 1L), 1L, s)
})
.py_try({
tf2 <- TransferFunction(.py_add(s, 2L), 1L, s)
})
.py_try({
fd1 <- Feedback(tf1, tf2, -1L)
})
.py_try({
fd2 <- Feedback(tf1, tf2, 1L)
})
.py_try({
unit <- TransferFunction(1L, 1L, s)
})
.py_try({
dtf1 <- DiscreteTransferFunction(.py_add(s, 1L), 1L, s)
})
.py_try({
dtf2 <- DiscreteTransferFunction(.py_add(s, 2L), 1L, s)
})
.py_try({
dfd1 <- Feedback(dtf1, dtf2, -1L)
})
.py_try({
dfd2 <- Feedback(dtf1, dtf2, 1L)
})
.py_try({
dunit <- DiscreteTransferFunction(1L, 1L, s)
})
.py_try({
.py_rec("lti:test_Feedback_with_Series:2444", function() .py_attr(tf1, "is_continuous"))
})
.py_try({
.py_rec("lti:test_Feedback_with_Series:2445", function() .py_attr(dtf1, "is_continuous"))
})
.py_try({
.py_rec("lti:test_Feedback_with_Series:2449", function() .py_isinstance(fd1, SISOLinearTimeInvariant))
})
.py_try({
.py_rec("lti:test_Feedback_with_Series:2450", function() .py_isinstance(fd1, Feedback))
})
.py_try({
.py_rec("lti:test_Feedback_with_Series:2453", function() .py_attr(fd1, "num"))
})
.py_try({
.py_rec("lti:test_Feedback_with_Series:2454", function() .py_attr(fd2, "num"))
})
.py_try({
.py_rec("lti:test_Feedback_with_Series:2455", function() .py_attr(fd1, "den"))
})
.py_try({
.py_rec("lti:test_Feedback_with_Series:2456", function() .py_attr(fd2, "den"))
})
.py_try({
s1 <- Series(tf1, fd1)
})
.py_try({
p1 <- Parallel(tf1, fd1)
})
.py_try({
.py_rec("lti:test_Feedback_with_Series:2461", function() .py_mul(tf1, fd1))
})
.py_try({
.py_rec("lti:test_Feedback_with_Series:2462", function() .py_add(tf1, fd1))
})
.py_try({
.py_rec("lti:test_Feedback_with_Series:2463", function() .py_method(s1, "doit"))
})
.py_try({
.py_rec("lti:test_Feedback_with_Series:2464", function() .py_method(p1, "doit"))
})
.py_try({
fd3 <- Feedback(.py_mul(tf1, fd1), tf2, -1L)
})
.py_try({
.py_rec("lti:test_Feedback_with_Series:2468", function() fd3)
})
.py_try({
.py_rec("lti:test_Feedback_with_Series:2469", function() .py_attr(fd3, "num"))
})
.py_try({
.py_rec("lti:test_Feedback_with_Series:2470", function() .py_attr(fd3, "den"))
})
.py_try({
tf3 <- TransferFunction(.py_mul(tf1, fd1), tf2, s)
})
.py_try({
.py_rec("lti:test_Feedback_with_Series:2474", function() tf3)
})
.py_try({
.py_rec("lti:test_Feedback_with_Series:2475", function() .py_attr(tf3, "num"))
})
.py_try({
.py_rec("lti:test_Feedback_with_Series:2479", function() .py_isinstance(dfd1, SISOLinearTimeInvariant))
})
.py_try({
.py_rec("lti:test_Feedback_with_Series:2480", function() .py_isinstance(dfd1, Feedback))
})
.py_try({
.py_rec("lti:test_Feedback_with_Series:2483", function() .py_attr(dfd1, "num"))
})
.py_try({
.py_rec("lti:test_Feedback_with_Series:2484", function() .py_attr(dfd2, "num"))
})
.py_try({
.py_rec("lti:test_Feedback_with_Series:2485", function() .py_attr(dfd1, "den"))
})
.py_try({
.py_rec("lti:test_Feedback_with_Series:2486", function() .py_attr(dfd2, "den"))
})
.py_try({
ds1 <- Series(dtf1, dfd1)
})
.py_try({
dp1 <- Parallel(dtf1, dfd1)
})
.py_try({
.py_rec("lti:test_Feedback_with_Series:2491", function() .py_mul(dtf1, dfd1))
})
.py_try({
.py_rec("lti:test_Feedback_with_Series:2492", function() .py_add(dtf1, dfd1))
})
.py_try({
.py_rec("lti:test_Feedback_with_Series:2493", function() .py_method(ds1, "doit"))
})
.py_try({
.py_rec("lti:test_Feedback_with_Series:2494", function() .py_method(dp1, "doit"))
})
.py_try({
dfd3 <- Feedback(.py_mul(dtf1, dfd1), dtf2, -1L)
})
.py_try({
.py_rec("lti:test_Feedback_with_Series:2500", function() dfd3)
})
.py_try({
.py_rec("lti:test_Feedback_with_Series:2501", function() .py_attr(dfd3, "num"))
})
.py_try({
.py_rec("lti:test_Feedback_with_Series:2502", function() .py_attr(dfd3, "den"))
})
.py_try({
dtf3 <- DiscreteTransferFunction(.py_mul(dtf1, dfd1), dtf2, s)
})
.py_try({
.py_rec("lti:test_Feedback_with_Series:2506", function() dtf3)
})
.py_try({
.py_rec("lti:test_Feedback_with_Series:2507", function() .py_attr(dtf3, "num"))
}))---")
  expect_identical(got[["lti:test_Feedback_with_Series:2444"]], "True")
  expect_identical(got[["lti:test_Feedback_with_Series:2445"]], "False")
  expect_identical(got[["lti:test_Feedback_with_Series:2449"]], "True")
  expect_identical(got[["lti:test_Feedback_with_Series:2450"]], "True")
  expect_identical(got[["lti:test_Feedback_with_Series:2453"]], "TransferFunction(s + 1, 1, s)")
  expect_identical(got[["lti:test_Feedback_with_Series:2454"]], "TransferFunction(s + 1, 1, s)")
  expect_identical(got[["lti:test_Feedback_with_Series:2455"]], "Parallel(TransferFunction(1, 1, s), Series(TransferFunction(s + 2, 1, s), TransferFunction(s + 1, 1, s)))")
  expect_identical(got[["lti:test_Feedback_with_Series:2456"]], "Parallel(TransferFunction(1, 1, s), Series(TransferFunction(-1, 1, s), TransferFunction(s + 2, 1, s), TransferFunction(s + 1, 1, s)))")
  expect_identical(got[["lti:test_Feedback_with_Series:2461"]], "Series(TransferFunction(s + 1, 1, s), Feedback(TransferFunction(s + 1, 1, s), TransferFunction(s + 2, 1, s), -1))")
  expect_identical(got[["lti:test_Feedback_with_Series:2462"]], "Parallel(TransferFunction(s + 1, 1, s), Feedback(TransferFunction(s + 1, 1, s), TransferFunction(s + 2, 1, s), -1))")
  expect_identical(got[["lti:test_Feedback_with_Series:2463"]], "TransferFunction((s + 1)^2, (s + 1)*(s + 2) + 1, s)")
  expect_identical(got[["lti:test_Feedback_with_Series:2464"]], "TransferFunction(s + (s + 1)*((s + 1)*(s + 2) + 1) + 1, (s + 1)*(s + 2) + 1, s)")
  expect_identical(got[["lti:test_Feedback_with_Series:2468"]], "Feedback(Series(TransferFunction(s + 1, 1, s), Feedback(TransferFunction(s + 1, 1, s), TransferFunction(s + 2, 1, s), -1)), TransferFunction(s + 2, 1, s), -1)")
  expect_identical(got[["lti:test_Feedback_with_Series:2469"]], "Series(TransferFunction(s + 1, 1, s), Feedback(TransferFunction(s + 1, 1, s), TransferFunction(s + 2, 1, s), -1))")
  expect_identical(got[["lti:test_Feedback_with_Series:2470"]], "Parallel(TransferFunction(1, 1, s), Series(TransferFunction(s + 2, 1, s), TransferFunction(s + 1, 1, s), Feedback(TransferFunction(s + 1, 1, s), TransferFunction(s + 2, 1, s), -1)))")
  expect_identical(got[["lti:test_Feedback_with_Series:2474"]], "TransferFunction(Series(TransferFunction(s + 1, 1, s), Feedback(TransferFunction(s + 1, 1, s), TransferFunction(s + 2, 1, s), -1)), TransferFunction(s + 2, 1, s), s)")
  expect_identical(got[["lti:test_Feedback_with_Series:2475"]], "Series(TransferFunction(s + 1, 1, s), Feedback(TransferFunction(s + 1, 1, s), TransferFunction(s + 2, 1, s), -1))")
  expect_identical(got[["lti:test_Feedback_with_Series:2479"]], "True")
  expect_identical(got[["lti:test_Feedback_with_Series:2480"]], "True")
  expect_identical(got[["lti:test_Feedback_with_Series:2483"]], "DiscreteTransferFunction(s + 1, 1, s, 1)")
  expect_identical(got[["lti:test_Feedback_with_Series:2484"]], "DiscreteTransferFunction(s + 1, 1, s, 1)")
  expect_identical(got[["lti:test_Feedback_with_Series:2485"]], "Parallel(DiscreteTransferFunction(1, 1, s, 1), Series(DiscreteTransferFunction(s + 2, 1, s, 1), DiscreteTransferFunction(s + 1, 1, s, 1)))")
  expect_identical(got[["lti:test_Feedback_with_Series:2486"]], "Parallel(DiscreteTransferFunction(1, 1, s, 1), Series(DiscreteTransferFunction(-1, 1, s, 1), DiscreteTransferFunction(s + 2, 1, s, 1), DiscreteTransferFunction(s + 1, 1, s, 1)))")
  expect_identical(got[["lti:test_Feedback_with_Series:2491"]], "Series(DiscreteTransferFunction(s + 1, 1, s, 1), Feedback(DiscreteTransferFunction(s + 1, 1, s, 1), DiscreteTransferFunction(s + 2, 1, s, 1), -1))")
  expect_identical(got[["lti:test_Feedback_with_Series:2492"]], "Parallel(DiscreteTransferFunction(s + 1, 1, s, 1), Feedback(DiscreteTransferFunction(s + 1, 1, s, 1), DiscreteTransferFunction(s + 2, 1, s, 1), -1))")
  expect_identical(got[["lti:test_Feedback_with_Series:2493"]], "DiscreteTransferFunction((s + 1)^2, (s + 1)*(s + 2) + 1, s, 1)")
  expect_identical(got[["lti:test_Feedback_with_Series:2494"]], "DiscreteTransferFunction(s + (s + 1)*((s + 1)*(s + 2) + 1) + 1, (s + 1)*(s + 2) + 1, s, 1)")
  expect_identical(got[["lti:test_Feedback_with_Series:2500"]], "Feedback(Series(DiscreteTransferFunction(s + 1, 1, s, 1), Feedback(DiscreteTransferFunction(s + 1, 1, s, 1), DiscreteTransferFunction(s + 2, 1, s, 1), -1)), DiscreteTransferFunction(s + 2, 1, s, 1), -1)")
  expect_identical(got[["lti:test_Feedback_with_Series:2501"]], "Series(DiscreteTransferFunction(s + 1, 1, s, 1), Feedback(DiscreteTransferFunction(s + 1, 1, s, 1), DiscreteTransferFunction(s + 2, 1, s, 1), -1))")
  expect_identical(got[["lti:test_Feedback_with_Series:2502"]], "Parallel(DiscreteTransferFunction(1, 1, s, 1), Series(DiscreteTransferFunction(s + 2, 1, s, 1), DiscreteTransferFunction(s + 1, 1, s, 1), Feedback(DiscreteTransferFunction(s + 1, 1, s, 1), DiscreteTransferFunction(s + 2, 1, s, 1), -1)))")
  expect_identical(got[["lti:test_Feedback_with_Series:2506"]], "DiscreteTransferFunction(Series(DiscreteTransferFunction(s + 1, 1, s, 1), Feedback(DiscreteTransferFunction(s + 1, 1, s, 1), DiscreteTransferFunction(s + 2, 1, s, 1), -1)), DiscreteTransferFunction(s + 2, 1, s, 1), s, 1)")
  expect_identical(got[["lti:test_Feedback_with_Series:2507"]], "Series(DiscreteTransferFunction(s + 1, 1, s, 1), Feedback(DiscreteTransferFunction(s + 1, 1, s, 1), DiscreteTransferFunction(s + 2, 1, s, 1), -1))")
})

test_that("lti test_issue_26161 matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("a, x, b, c, s, g, d,    p, k, tau, zeta, wn, T, z"); a <- .py_getitem(.tmpu_, 0L); x <- .py_getitem(.tmpu_, 1L); b <- .py_getitem(.tmpu_, 2L); c <- .py_getitem(.tmpu_, 3L); s <- .py_getitem(.tmpu_, 4L); g <- .py_getitem(.tmpu_, 5L); d <- .py_getitem(.tmpu_, 6L); p <- .py_getitem(.tmpu_, 7L); k <- .py_getitem(.tmpu_, 8L); tau <- .py_getitem(.tmpu_, 9L); zeta <- .py_getitem(.tmpu_, 10L); wn <- .py_getitem(.tmpu_, 11L); T <- .py_getitem(.tmpu_, 12L); z <- .py_getitem(.tmpu_, 13L)
})
.py_try({
.tmpu_ <- .py_symbols("a0:4, b0:5, c0:4, d0:4"); a0 <- .py_getitem(.tmpu_, 0L); a1 <- .py_getitem(.tmpu_, 1L); a2 <- .py_getitem(.tmpu_, 2L); a3 <- .py_getitem(.tmpu_, 3L); b0 <- .py_getitem(.tmpu_, 4L); b1 <- .py_getitem(.tmpu_, 5L); b2 <- .py_getitem(.tmpu_, 6L); b3 <- .py_getitem(.tmpu_, 7L); b4 <- .py_getitem(.tmpu_, 8L); c0 <- .py_getitem(.tmpu_, 9L); c1 <- .py_getitem(.tmpu_, 10L); c2 <- .py_getitem(.tmpu_, 11L); c3 <- .py_getitem(.tmpu_, 12L); d0 <- .py_getitem(.tmpu_, 13L); d1 <- .py_getitem(.tmpu_, 14L); d2 <- .py_getitem(.tmpu_, 15L); d3 <- .py_getitem(.tmpu_, 16L)
})
.py_try({
TF1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
TF2 <- TransferFunction(k, 1L, s)
})
.py_try({
TF3 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
TF4 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(tau, s), 1L), s)
})
.py_try({
TF5 <- TransferFunction(k, p, s)
})
.py_try({
TF6 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(k, s), p), s)
})
p_are_floats_equals <- function(n1, n2) {
  "\n    Helper function to compare two floats.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return(.py_method(Float(n1), "is_same", Float(n2), .py_isclose))
}
p_are_complex_equals <- function(n1, n2) {
  "\n    Helper function to compare two complex numbers.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return((.py_truthy(p_are_floats_equals(re(n1), re(n2))) && .py_truthy(p_are_floats_equals(im(n1), im(n2)))))
}
.py_try({
.tmpu_ <- .py_symbols("I_b, I_s, m, h, l2, l1", real = TRUE, nonnegative = TRUE); Ib <- .py_getitem(.tmpu_, 0L); Is <- .py_getitem(.tmpu_, 1L); m <- .py_getitem(.tmpu_, 2L); h <- .py_getitem(.tmpu_, 3L); l2 <- .py_getitem(.tmpu_, 4L); l1 <- .py_getitem(.tmpu_, 5L)
})
.py_try({
.tmpu_ <- .py_symbols("K_D, K_P, v", real = TRUE); KD <- .py_getitem(.tmpu_, 0L); KP <- .py_getitem(.tmpu_, 1L); v <- .py_getitem(.tmpu_, 2L)
})
.py_try({
tau1_sq <- .py_div(.py_div(.py_div(.py_add(Ib, .py_mul(m, .py_pow(h, 2L))), m), g), h)
})
.py_try({
tau2 <- .py_div(l2, v)
})
.py_try({
tau3 <- .py_div(v, .py_add(l1, l2))
})
.py_try({
K <- .py_div(.py_div(.py_pow(v, 2L), g), .py_add(l1, l2))
})
.py_try({
Gtheta <- TransferFunction(.py_mul(.py_neg(K), .py_add(.py_mul(tau2, s), 1L)), .py_sub(.py_mul(tau1_sq, .py_pow(s, 2L)), 1L), s)
})
.py_try({
Gdelta <- TransferFunction(1L, .py_add(.py_mul(Is, .py_pow(s, 2L)), .py_mul(c, s)), s)
})
.py_try({
Gpsi <- TransferFunction(1L, .py_mul(tau3, s), s)
})
.py_try({
Dcont <- TransferFunction(.py_mul(KD, s), 1L, s)
})
.py_try({
PIcont <- TransferFunction(KP, s, s)
})
.py_try({
Gunity <- TransferFunction(1L, 1L, s)
})
.py_try({
Ginner <- Feedback(.py_mul(Dcont, Gdelta), Gtheta)
})
.py_try({
Gouter <- Feedback(.py_mul(.py_mul(PIcont, Ginner), Gpsi), Gunity)
})
.py_try({
.py_rec("lti:test_issue_26161:2530", function() Gouter)
})
.py_try({
.py_rec("lti:test_issue_26161:2531", function() .py_attr(Gouter, "num"))
})
.py_try({
.py_rec("lti:test_issue_26161:2532", function() .py_attr(Gouter, "den"))
})
.py_try({
expr <- .py_div(.py_mul(.py_mul(.py_mul(.py_mul(.py_mul(.py_mul(.py_mul(.py_mul(KD, KP), g), .py_pow(s, 3L)), .py_pow(v, 2L)), .py_add(l1, l2)), .py_pow(.py_add(.py_mul(Is, .py_pow(s, 2L)), .py_mul(c, s)), 2L)), .py_add(.py_mul(.py_mul(.py_neg(g), h), m), .py_mul(.py_pow(s, 2L), .py_add(Ib, .py_mul(.py_pow(h, 2L), m))))), .py_add(.py_mul(.py_mul(.py_mul(.py_mul(.py_mul(.py_mul(.py_neg(KD), g), h), m), s), .py_pow(v, 2L)), .py_add(.py_mul(l2, s), v)), .py_mul(.py_mul(.py_mul(.py_mul(g, v), .py_add(l1, l2)), .py_add(.py_mul(Is, .py_pow(s, 2L)), .py_mul(c, s))), .py_add(.py_mul(.py_mul(.py_neg(g), h), m), .py_mul(.py_pow(s, 2L), .py_add(Ib, .py_mul(.py_pow(h, 2L), m))))))), .py_div(.py_mul(.py_mul(.py_mul(.py_mul(.py_pow(s, 2L), v), .py_add(.py_mul(Is, .py_pow(s, 2L)), .py_mul(c, s))), .py_add(.py_mul(.py_mul(.py_mul(.py_mul(.py_mul(.py_mul(.py_neg(KD), g), h), m), s), .py_pow(v, 2L)), .py_add(.py_mul(l2, s), v)), .py_mul(.py_mul(.py_mul(.py_mul(g, v), .py_add(l1, l2)), .py_add(.py_mul(Is, .py_pow(s, 2L)), .py_mul(c, s))), .py_add(.py_mul(.py_mul(.py_neg(g), h), m), .py_mul(.py_pow(s, 2L), .py_add(Ib, .py_mul(.py_pow(h, 2L), m))))))), .py_add(.py_mul(.py_mul(.py_mul(.py_mul(.py_mul(.py_mul(.py_mul(KD, KP), g), s), v), .py_pow(.py_add(l1, l2), 2L)), .py_add(.py_mul(Is, .py_pow(s, 2L)), .py_mul(c, s))), .py_add(.py_mul(.py_mul(.py_neg(g), h), m), .py_mul(.py_pow(s, 2L), .py_add(Ib, .py_mul(.py_pow(h, 2L), m))))), .py_mul(.py_mul(.py_mul(.py_pow(s, 2L), v), .py_add(.py_mul(Is, .py_pow(s, 2L)), .py_mul(c, s))), .py_add(.py_mul(.py_mul(.py_mul(.py_mul(.py_mul(.py_mul(.py_neg(KD), g), h), m), s), .py_pow(v, 2L)), .py_add(.py_mul(l2, s), v)), .py_mul(.py_mul(.py_mul(.py_mul(g, v), .py_add(l1, l2)), .py_add(.py_mul(Is, .py_pow(s, 2L)), .py_mul(c, s))), .py_add(.py_mul(.py_mul(.py_neg(g), h), m), .py_mul(.py_pow(s, 2L), .py_add(Ib, .py_mul(.py_pow(h, 2L), m))))))))), .py_add(l1, l2)))
})
.py_try({
.py_rec("lti:test_issue_26161:2539", function() .py_method(.py_sub(.py_method(Gouter, "to_expr"), expr), "simplify"))
}))---")
  expect_identical(got[["lti:test_issue_26161:2530"]], "Feedback(Series(TransferFunction(K_P, s, s), Feedback(Series(TransferFunction(K_D*s, 1, s), TransferFunction(1, I_s*s^2 + c*s, s)), TransferFunction(-v^2*(l2*s/v + 1)/(g*(l1 + l2)), -1 + s^2*(I_b + h^2*m)/(g*h*m), s), -1), TransferFunction(1, s*v/(l1 + l2), s)), TransferFunction(1, 1, s), -1)")
  expect_identical(got[["lti:test_issue_26161:2531"]], "Series(TransferFunction(K_P, s, s), Feedback(Series(TransferFunction(K_D*s, 1, s), TransferFunction(1, I_s*s^2 + c*s, s)), TransferFunction(-v^2*(l2*s/v + 1)/(g*(l1 + l2)), -1 + s^2*(I_b + h^2*m)/(g*h*m), s), -1), TransferFunction(1, s*v/(l1 + l2), s))")
  expect_identical(got[["lti:test_issue_26161:2532"]], "Parallel(TransferFunction(1, 1, s), Series(TransferFunction(1, 1, s), TransferFunction(K_P, s, s), Feedback(Series(TransferFunction(K_D*s, 1, s), TransferFunction(1, I_s*s^2 + c*s, s)), TransferFunction(-v^2*(l2*s/v + 1)/(g*(l1 + l2)), -1 + s^2*(I_b + h^2*m)/(g*h*m), s), -1), TransferFunction(1, s*v/(l1 + l2), s)))")
  expect_identical(got[["lti:test_issue_26161:2539"]], "0")
})

test_that("lti test_MIMOFeedback_construction matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("a, x, b, c, s, g, d,    p, k, tau, zeta, wn, T, z"); a <- .py_getitem(.tmpu_, 0L); x <- .py_getitem(.tmpu_, 1L); b <- .py_getitem(.tmpu_, 2L); c <- .py_getitem(.tmpu_, 3L); s <- .py_getitem(.tmpu_, 4L); g <- .py_getitem(.tmpu_, 5L); d <- .py_getitem(.tmpu_, 6L); p <- .py_getitem(.tmpu_, 7L); k <- .py_getitem(.tmpu_, 8L); tau <- .py_getitem(.tmpu_, 9L); zeta <- .py_getitem(.tmpu_, 10L); wn <- .py_getitem(.tmpu_, 11L); T <- .py_getitem(.tmpu_, 12L); z <- .py_getitem(.tmpu_, 13L)
})
.py_try({
.tmpu_ <- .py_symbols("a0:4, b0:5, c0:4, d0:4"); a0 <- .py_getitem(.tmpu_, 0L); a1 <- .py_getitem(.tmpu_, 1L); a2 <- .py_getitem(.tmpu_, 2L); a3 <- .py_getitem(.tmpu_, 3L); b0 <- .py_getitem(.tmpu_, 4L); b1 <- .py_getitem(.tmpu_, 5L); b2 <- .py_getitem(.tmpu_, 6L); b3 <- .py_getitem(.tmpu_, 7L); b4 <- .py_getitem(.tmpu_, 8L); c0 <- .py_getitem(.tmpu_, 9L); c1 <- .py_getitem(.tmpu_, 10L); c2 <- .py_getitem(.tmpu_, 11L); c3 <- .py_getitem(.tmpu_, 12L); d0 <- .py_getitem(.tmpu_, 13L); d1 <- .py_getitem(.tmpu_, 14L); d2 <- .py_getitem(.tmpu_, 15L); d3 <- .py_getitem(.tmpu_, 16L)
})
.py_try({
TF1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
TF2 <- TransferFunction(k, 1L, s)
})
.py_try({
TF3 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
TF4 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(tau, s), 1L), s)
})
.py_try({
TF5 <- TransferFunction(k, p, s)
})
.py_try({
TF6 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(k, s), p), s)
})
p_are_floats_equals <- function(n1, n2) {
  "\n    Helper function to compare two floats.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return(.py_method(Float(n1), "is_same", Float(n2), .py_isclose))
}
p_are_complex_equals <- function(n1, n2) {
  "\n    Helper function to compare two complex numbers.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return((.py_truthy(p_are_floats_equals(re(n1), re(n2))) && .py_truthy(p_are_floats_equals(im(n1), im(n2)))))
}
.py_try({
tf1 <- TransferFunction(1L, s, s)
})
.py_try({
tf2 <- TransferFunction(s, .py_sub(.py_pow(s, 3L), 1L), s)
})
.py_try({
tf3 <- TransferFunction(s, .py_add(s, 1L), s)
})
.py_try({
tf4 <- TransferFunction(s, .py_add(.py_pow(s, 2L), 1L), s)
})
.py_try({
dtf1 <- DiscreteTransferFunction(1L, s, s, .py_flt(0.5))
})
.py_try({
dtf2 <- DiscreteTransferFunction(s, .py_sub(.py_pow(s, 3L), 1L), s, .py_flt(0.5))
})
.py_try({
dtf3 <- DiscreteTransferFunction(s, .py_add(s, 1L), s, .py_flt(0.5))
})
.py_try({
dtf4 <- DiscreteTransferFunction(s, .py_add(.py_pow(s, 2L), 1L), s, .py_flt(0.5))
})
.py_try({
tfm_1 <- TransferFunctionMatrix(list(list(tf1, tf2), list(tf3, tf4)))
})
.py_try({
tfm_2 <- TransferFunctionMatrix(list(list(tf2, tf3), list(tf4, tf1)))
})
.py_try({
tfm_3 <- TransferFunctionMatrix(list(list(tf3, tf4), list(tf1, tf2)))
})
.py_try({
dtfm_1 <- TransferFunctionMatrix(list(list(dtf1, dtf2), list(dtf3, dtf4)))
})
.py_try({
dtfm_2 <- TransferFunctionMatrix(list(list(dtf2, dtf3), list(dtf4, dtf1)))
})
.py_try({
dtfm_3 <- TransferFunctionMatrix(list(list(dtf3, dtf4), list(dtf1, dtf2)))
})
.py_try({
f1 <- MIMOFeedback(tfm_1, tfm_2)
})
.py_try({
.py_rec("lti:test_MIMOFeedback_construction:2563", function() .py_attr(f1, "args"))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_construction:2564", function() .py_attr(f1, "sys1"))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_construction:2565", function() .py_attr(f1, "sys2"))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_construction:2566", function() .py_attr(f1, "var"))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_construction:2567", function() .py_attr(f1, "sign"))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_construction:2568", function() .py_neg(.py_neg(f1)))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_construction:2569", function() .py_attr(f1, "is_continuous"))
})
.py_try({
f2 <- MIMOFeedback(tfm_2, tfm_1, 1L)
})
.py_try({
.py_rec("lti:test_MIMOFeedback_construction:2572", function() .py_attr(f2, "args"))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_construction:2573", function() .py_attr(f2, "sys1"))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_construction:2574", function() .py_attr(f2, "sys2"))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_construction:2575", function() .py_attr(f2, "var"))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_construction:2576", function() .py_attr(f2, "sign"))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_construction:2577", function() .py_attr(f2, "is_continuous"))
})
.py_try({
f3 <- MIMOFeedback(tfm_1, MIMOSeries(tfm_3, tfm_2))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_construction:2580", function() .py_attr(f3, "args"))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_construction:2581", function() .py_attr(f3, "sys1"))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_construction:2582", function() .py_attr(f3, "sys2"))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_construction:2583", function() .py_attr(f3, "var"))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_construction:2584", function() .py_attr(f3, "sign"))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_construction:2585", function() .py_attr(f3, "is_continuous"))
})
.py_try({
mat <- .py_Matrix(list(list(1L, .py_div(1L, s)), list(0L, 1L)))
})
.py_try({
.tmp_ <- .py_method(TransferFunctionMatrix, "from_Matrix", mat, s)
sys1 <- .tmp_
controller <- .tmp_
})
.py_try({
f4 <- MIMOFeedback(sys1, controller)
})
.py_try({
.py_rec("lti:test_MIMOFeedback_construction:2590", function() .py_attr(f4, "args"))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_construction:2591", function() (.py_truthy(.py_eq(.py_attr(f4, "sys1"), .py_attr(f4, "sys2"))) && .py_truthy(.py_eq(.py_attr(f4, "sys2"), sys1))))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_construction:2592", function() .py_attr(f4, "is_continuous"))
})
.py_try({
df1 <- MIMOFeedback(dtfm_1, dtfm_2)
})
.py_try({
.py_rec("lti:test_MIMOFeedback_construction:2596", function() .py_attr(df1, "args"))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_construction:2597", function() .py_attr(df1, "sys1"))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_construction:2598", function() .py_attr(df1, "sys2"))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_construction:2599", function() .py_attr(df1, "var"))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_construction:2600", function() .py_attr(df1, "sign"))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_construction:2601", function() .py_attr(df1, "sampling_time"))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_construction:2602", function() .py_neg(.py_neg(df1)))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_construction:2603", function() .py_attr(df1, "is_continuous"))
})
.py_try({
df2 <- MIMOFeedback(dtfm_2, dtfm_1, 1L)
})
.py_try({
.py_rec("lti:test_MIMOFeedback_construction:2606", function() .py_attr(df2, "args"))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_construction:2607", function() .py_attr(df2, "sys1"))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_construction:2608", function() .py_attr(df2, "sys2"))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_construction:2609", function() .py_attr(df2, "var"))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_construction:2610", function() .py_attr(df2, "sign"))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_construction:2611", function() .py_attr(df2, "sampling_time"))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_construction:2612", function() .py_attr(df2, "is_continuous"))
})
.py_try({
df3 <- MIMOFeedback(dtfm_1, MIMOSeries(dtfm_3, dtfm_2))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_construction:2615", function() .py_attr(df3, "args"))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_construction:2616", function() .py_attr(df3, "sys1"))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_construction:2617", function() .py_attr(df3, "sys2"))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_construction:2618", function() .py_attr(df3, "var"))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_construction:2619", function() .py_attr(df3, "sign"))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_construction:2620", function() .py_attr(df3, "sampling_time"))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_construction:2621", function() .py_attr(df3, "is_continuous"))
})
.py_try({
dmat <- .py_Matrix(list(list(1L, .py_div(1L, s)), list(0L, 1L)))
})
.py_try({
.tmp_ <- .py_method(TransferFunctionMatrix, "from_Matrix", dmat, s, .py_flt(0.5))
dsys1 <- .tmp_
dcontroller <- .tmp_
})
.py_try({
df4 <- MIMOFeedback(dsys1, dcontroller)
})
.py_try({
.py_rec("lti:test_MIMOFeedback_construction:2627", function() .py_attr(df4, "args"))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_construction:2628", function() (.py_truthy(.py_eq(.py_attr(df4, "sys1"), .py_attr(df4, "sys2"))) && .py_truthy(.py_eq(.py_attr(df4, "sys2"), dsys1))))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_construction:2629", function() .py_attr(df4, "sampling_time"))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_construction:2630", function() .py_attr(df4, "is_continuous"))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_construction:2633", function() MIMOFeedback(dtfm_1, tfm_2))
})
.py_try({
dtfm_4 <- .py_method(TransferFunctionMatrix, "from_Matrix", .py_mul(eye(2L), s), var = s, sampling_time = .py_flt(0.1))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_construction:2636", function() MIMOFeedback(dtfm_1, dtfm_4))
}))---")
  expect_identical(got[["lti:test_MIMOFeedback_construction:2563"]], "(TransferFunctionMatrix(((TransferFunction(1, s, s), TransferFunction(s, s^3 - 1, s)), (TransferFunction(s, s + 1, s), TransferFunction(s, s^2 + 1, s)))), TransferFunctionMatrix(((TransferFunction(s, s^3 - 1, s), TransferFunction(s, s + 1, s)), (TransferFunction(s, s^2 + 1, s), TransferFunction(1, s, s)))), -1)")
  expect_identical(got[["lti:test_MIMOFeedback_construction:2564"]], "TransferFunctionMatrix(((TransferFunction(1, s, s), TransferFunction(s, s^3 - 1, s)), (TransferFunction(s, s + 1, s), TransferFunction(s, s^2 + 1, s))))")
  expect_identical(got[["lti:test_MIMOFeedback_construction:2565"]], "TransferFunctionMatrix(((TransferFunction(s, s^3 - 1, s), TransferFunction(s, s + 1, s)), (TransferFunction(s, s^2 + 1, s), TransferFunction(1, s, s))))")
  expect_identical(got[["lti:test_MIMOFeedback_construction:2566"]], "s")
  expect_identical(got[["lti:test_MIMOFeedback_construction:2567"]], "-1")
  expect_identical(got[["lti:test_MIMOFeedback_construction:2568"]], "MIMOFeedback(TransferFunctionMatrix(((TransferFunction(1, s, s), TransferFunction(s, s^3 - 1, s)), (TransferFunction(s, s + 1, s), TransferFunction(s, s^2 + 1, s)))), TransferFunctionMatrix(((TransferFunction(s, s^3 - 1, s), TransferFunction(s, s + 1, s)), (TransferFunction(s, s^2 + 1, s), TransferFunction(1, s, s)))), -1)")
  expect_identical(got[["lti:test_MIMOFeedback_construction:2569"]], "True")
  expect_identical(got[["lti:test_MIMOFeedback_construction:2572"]], "(TransferFunctionMatrix(((TransferFunction(s, s^3 - 1, s), TransferFunction(s, s + 1, s)), (TransferFunction(s, s^2 + 1, s), TransferFunction(1, s, s)))), TransferFunctionMatrix(((TransferFunction(1, s, s), TransferFunction(s, s^3 - 1, s)), (TransferFunction(s, s + 1, s), TransferFunction(s, s^2 + 1, s)))), 1)")
  expect_identical(got[["lti:test_MIMOFeedback_construction:2573"]], "TransferFunctionMatrix(((TransferFunction(s, s^3 - 1, s), TransferFunction(s, s + 1, s)), (TransferFunction(s, s^2 + 1, s), TransferFunction(1, s, s))))")
  expect_identical(got[["lti:test_MIMOFeedback_construction:2574"]], "TransferFunctionMatrix(((TransferFunction(1, s, s), TransferFunction(s, s^3 - 1, s)), (TransferFunction(s, s + 1, s), TransferFunction(s, s^2 + 1, s))))")
  expect_identical(got[["lti:test_MIMOFeedback_construction:2575"]], "s")
  expect_identical(got[["lti:test_MIMOFeedback_construction:2576"]], "1")
  expect_identical(got[["lti:test_MIMOFeedback_construction:2577"]], "True")
  expect_identical(got[["lti:test_MIMOFeedback_construction:2580"]], "(TransferFunctionMatrix(((TransferFunction(1, s, s), TransferFunction(s, s^3 - 1, s)), (TransferFunction(s, s + 1, s), TransferFunction(s, s^2 + 1, s)))), MIMOSeries(TransferFunctionMatrix(((TransferFunction(s, s + 1, s), TransferFunction(s, s^2 + 1, s)), (TransferFunction(1, s, s), TransferFunction(s, s^3 - 1, s)))), TransferFunctionMatrix(((TransferFunction(s, s^3 - 1, s), TransferFunction(s, s + 1, s)), (TransferFunction(s, s^2 + 1, s), TransferFunction(1, s, s))))), -1)")
  expect_identical(got[["lti:test_MIMOFeedback_construction:2581"]], "TransferFunctionMatrix(((TransferFunction(1, s, s), TransferFunction(s, s^3 - 1, s)), (TransferFunction(s, s + 1, s), TransferFunction(s, s^2 + 1, s))))")
  expect_identical(got[["lti:test_MIMOFeedback_construction:2582"]], "MIMOSeries(TransferFunctionMatrix(((TransferFunction(s, s + 1, s), TransferFunction(s, s^2 + 1, s)), (TransferFunction(1, s, s), TransferFunction(s, s^3 - 1, s)))), TransferFunctionMatrix(((TransferFunction(s, s^3 - 1, s), TransferFunction(s, s + 1, s)), (TransferFunction(s, s^2 + 1, s), TransferFunction(1, s, s)))))")
  expect_identical(got[["lti:test_MIMOFeedback_construction:2583"]], "s")
  expect_identical(got[["lti:test_MIMOFeedback_construction:2584"]], "-1")
  expect_identical(got[["lti:test_MIMOFeedback_construction:2585"]], "True")
  expect_identical(got[["lti:test_MIMOFeedback_construction:2590"]], "(TransferFunctionMatrix(((TransferFunction(1, 1, s), TransferFunction(1, s, s)), (TransferFunction(0, 1, s), TransferFunction(1, 1, s)))), TransferFunctionMatrix(((TransferFunction(1, 1, s), TransferFunction(1, s, s)), (TransferFunction(0, 1, s), TransferFunction(1, 1, s)))), -1)")
  expect_identical(got[["lti:test_MIMOFeedback_construction:2591"]], "True")
  expect_identical(got[["lti:test_MIMOFeedback_construction:2592"]], "True")
  expect_identical(got[["lti:test_MIMOFeedback_construction:2596"]], "(TransferFunctionMatrix(((DiscreteTransferFunction(1, s, s, 0.5), DiscreteTransferFunction(s, s^3 - 1, s, 0.5)), (DiscreteTransferFunction(s, s + 1, s, 0.5), DiscreteTransferFunction(s, s^2 + 1, s, 0.5)))), TransferFunctionMatrix(((DiscreteTransferFunction(s, s^3 - 1, s, 0.5), DiscreteTransferFunction(s, s + 1, s, 0.5)), (DiscreteTransferFunction(s, s^2 + 1, s, 0.5), DiscreteTransferFunction(1, s, s, 0.5)))), -1)")
  expect_identical(got[["lti:test_MIMOFeedback_construction:2597"]], "TransferFunctionMatrix(((DiscreteTransferFunction(1, s, s, 0.5), DiscreteTransferFunction(s, s^3 - 1, s, 0.5)), (DiscreteTransferFunction(s, s + 1, s, 0.5), DiscreteTransferFunction(s, s^2 + 1, s, 0.5))))")
  expect_identical(got[["lti:test_MIMOFeedback_construction:2598"]], "TransferFunctionMatrix(((DiscreteTransferFunction(s, s^3 - 1, s, 0.5), DiscreteTransferFunction(s, s + 1, s, 0.5)), (DiscreteTransferFunction(s, s^2 + 1, s, 0.5), DiscreteTransferFunction(1, s, s, 0.5))))")
  expect_identical(got[["lti:test_MIMOFeedback_construction:2599"]], "s")
  expect_identical(got[["lti:test_MIMOFeedback_construction:2600"]], "-1")
  expect_identical(got[["lti:test_MIMOFeedback_construction:2601"]], "0.500000000000000")
  expect_identical(got[["lti:test_MIMOFeedback_construction:2602"]], "MIMOFeedback(TransferFunctionMatrix(((DiscreteTransferFunction(1, s, s, 0.5), DiscreteTransferFunction(s, s^3 - 1, s, 0.5)), (DiscreteTransferFunction(s, s + 1, s, 0.5), DiscreteTransferFunction(s, s^2 + 1, s, 0.5)))), TransferFunctionMatrix(((DiscreteTransferFunction(s, s^3 - 1, s, 0.5), DiscreteTransferFunction(s, s + 1, s, 0.5)), (DiscreteTransferFunction(s, s^2 + 1, s, 0.5), DiscreteTransferFunction(1, s, s, 0.5)))), -1)")
  expect_identical(got[["lti:test_MIMOFeedback_construction:2603"]], "False")
  expect_identical(got[["lti:test_MIMOFeedback_construction:2606"]], "(TransferFunctionMatrix(((DiscreteTransferFunction(s, s^3 - 1, s, 0.5), DiscreteTransferFunction(s, s + 1, s, 0.5)), (DiscreteTransferFunction(s, s^2 + 1, s, 0.5), DiscreteTransferFunction(1, s, s, 0.5)))), TransferFunctionMatrix(((DiscreteTransferFunction(1, s, s, 0.5), DiscreteTransferFunction(s, s^3 - 1, s, 0.5)), (DiscreteTransferFunction(s, s + 1, s, 0.5), DiscreteTransferFunction(s, s^2 + 1, s, 0.5)))), 1)")
  expect_identical(got[["lti:test_MIMOFeedback_construction:2607"]], "TransferFunctionMatrix(((DiscreteTransferFunction(s, s^3 - 1, s, 0.5), DiscreteTransferFunction(s, s + 1, s, 0.5)), (DiscreteTransferFunction(s, s^2 + 1, s, 0.5), DiscreteTransferFunction(1, s, s, 0.5))))")
  expect_identical(got[["lti:test_MIMOFeedback_construction:2608"]], "TransferFunctionMatrix(((DiscreteTransferFunction(1, s, s, 0.5), DiscreteTransferFunction(s, s^3 - 1, s, 0.5)), (DiscreteTransferFunction(s, s + 1, s, 0.5), DiscreteTransferFunction(s, s^2 + 1, s, 0.5))))")
  expect_identical(got[["lti:test_MIMOFeedback_construction:2609"]], "s")
  expect_identical(got[["lti:test_MIMOFeedback_construction:2610"]], "1")
  expect_identical(got[["lti:test_MIMOFeedback_construction:2611"]], "0.500000000000000")
  expect_identical(got[["lti:test_MIMOFeedback_construction:2612"]], "False")
  expect_identical(got[["lti:test_MIMOFeedback_construction:2615"]], "(TransferFunctionMatrix(((DiscreteTransferFunction(1, s, s, 0.5), DiscreteTransferFunction(s, s^3 - 1, s, 0.5)), (DiscreteTransferFunction(s, s + 1, s, 0.5), DiscreteTransferFunction(s, s^2 + 1, s, 0.5)))), MIMOSeries(TransferFunctionMatrix(((DiscreteTransferFunction(s, s + 1, s, 0.5), DiscreteTransferFunction(s, s^2 + 1, s, 0.5)), (DiscreteTransferFunction(1, s, s, 0.5), DiscreteTransferFunction(s, s^3 - 1, s, 0.5)))), TransferFunctionMatrix(((DiscreteTransferFunction(s, s^3 - 1, s, 0.5), DiscreteTransferFunction(s, s + 1, s, 0.5)), (DiscreteTransferFunction(s, s^2 + 1, s, 0.5), DiscreteTransferFunction(1, s, s, 0.5))))), -1)")
  expect_identical(got[["lti:test_MIMOFeedback_construction:2616"]], "TransferFunctionMatrix(((DiscreteTransferFunction(1, s, s, 0.5), DiscreteTransferFunction(s, s^3 - 1, s, 0.5)), (DiscreteTransferFunction(s, s + 1, s, 0.5), DiscreteTransferFunction(s, s^2 + 1, s, 0.5))))")
  expect_identical(got[["lti:test_MIMOFeedback_construction:2617"]], "MIMOSeries(TransferFunctionMatrix(((DiscreteTransferFunction(s, s + 1, s, 0.5), DiscreteTransferFunction(s, s^2 + 1, s, 0.5)), (DiscreteTransferFunction(1, s, s, 0.5), DiscreteTransferFunction(s, s^3 - 1, s, 0.5)))), TransferFunctionMatrix(((DiscreteTransferFunction(s, s^3 - 1, s, 0.5), DiscreteTransferFunction(s, s + 1, s, 0.5)), (DiscreteTransferFunction(s, s^2 + 1, s, 0.5), DiscreteTransferFunction(1, s, s, 0.5)))))")
  expect_identical(got[["lti:test_MIMOFeedback_construction:2618"]], "s")
  expect_identical(got[["lti:test_MIMOFeedback_construction:2619"]], "-1")
  expect_identical(got[["lti:test_MIMOFeedback_construction:2620"]], "0.500000000000000")
  expect_identical(got[["lti:test_MIMOFeedback_construction:2621"]], "False")
  expect_identical(got[["lti:test_MIMOFeedback_construction:2627"]], "(TransferFunctionMatrix(((DiscreteTransferFunction(1, 1, s, 0.5), DiscreteTransferFunction(1, s, s, 0.5)), (DiscreteTransferFunction(0, 1, s, 0.5), DiscreteTransferFunction(1, 1, s, 0.5)))), TransferFunctionMatrix(((DiscreteTransferFunction(1, 1, s, 0.5), DiscreteTransferFunction(1, s, s, 0.5)), (DiscreteTransferFunction(0, 1, s, 0.5), DiscreteTransferFunction(1, 1, s, 0.5)))), -1)")
  expect_identical(got[["lti:test_MIMOFeedback_construction:2628"]], "True")
  expect_identical(got[["lti:test_MIMOFeedback_construction:2629"]], "0.500000000000000")
  expect_identical(got[["lti:test_MIMOFeedback_construction:2630"]], "False")
  expect_identical(got[["lti:test_MIMOFeedback_construction:2633"]], "ERR TypeError: \\n                Incompatible systems detected.\\n                All systems should be either continuous-time or discrete-time.\\n                Found: TransferFunctionMatrix(((DiscreteTransferFunction(1, s, s, 0.5), DiscreteTransferFunction(s, s^3 - 1, s, 0.5)), (DiscreteTransferFunction(s, s + 1, s, 0.5), DiscreteTransferFunction(s, s^2 + 1, s, 0.5)))) and TransferFunctionMatrix(((TransferFunction(s, s^3 - 1, s), TransferFunction(s, s + 1, s)), (TransferFunction(s, s^2 + 1, s), TransferFunction(1, s, s))))")
  expect_identical(got[["lti:test_MIMOFeedback_construction:2636"]], "ERR TypeError: \\n                Incompatible sampling times detected.\\n                All systems should have the same sampling time.\\n                Found sampling time: 0.500000000000000 and 0.100000000000000")
})

test_that("lti test_MIMOFeedback_errors matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("a, x, b, c, s, g, d,    p, k, tau, zeta, wn, T, z"); a <- .py_getitem(.tmpu_, 0L); x <- .py_getitem(.tmpu_, 1L); b <- .py_getitem(.tmpu_, 2L); c <- .py_getitem(.tmpu_, 3L); s <- .py_getitem(.tmpu_, 4L); g <- .py_getitem(.tmpu_, 5L); d <- .py_getitem(.tmpu_, 6L); p <- .py_getitem(.tmpu_, 7L); k <- .py_getitem(.tmpu_, 8L); tau <- .py_getitem(.tmpu_, 9L); zeta <- .py_getitem(.tmpu_, 10L); wn <- .py_getitem(.tmpu_, 11L); T <- .py_getitem(.tmpu_, 12L); z <- .py_getitem(.tmpu_, 13L)
})
.py_try({
.tmpu_ <- .py_symbols("a0:4, b0:5, c0:4, d0:4"); a0 <- .py_getitem(.tmpu_, 0L); a1 <- .py_getitem(.tmpu_, 1L); a2 <- .py_getitem(.tmpu_, 2L); a3 <- .py_getitem(.tmpu_, 3L); b0 <- .py_getitem(.tmpu_, 4L); b1 <- .py_getitem(.tmpu_, 5L); b2 <- .py_getitem(.tmpu_, 6L); b3 <- .py_getitem(.tmpu_, 7L); b4 <- .py_getitem(.tmpu_, 8L); c0 <- .py_getitem(.tmpu_, 9L); c1 <- .py_getitem(.tmpu_, 10L); c2 <- .py_getitem(.tmpu_, 11L); c3 <- .py_getitem(.tmpu_, 12L); d0 <- .py_getitem(.tmpu_, 13L); d1 <- .py_getitem(.tmpu_, 14L); d2 <- .py_getitem(.tmpu_, 15L); d3 <- .py_getitem(.tmpu_, 16L)
})
.py_try({
TF1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
TF2 <- TransferFunction(k, 1L, s)
})
.py_try({
TF3 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
TF4 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(tau, s), 1L), s)
})
.py_try({
TF5 <- TransferFunction(k, p, s)
})
.py_try({
TF6 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(k, s), p), s)
})
p_are_floats_equals <- function(n1, n2) {
  "\n    Helper function to compare two floats.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return(.py_method(Float(n1), "is_same", Float(n2), .py_isclose))
}
p_are_complex_equals <- function(n1, n2) {
  "\n    Helper function to compare two complex numbers.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return((.py_truthy(p_are_floats_equals(re(n1), re(n2))) && .py_truthy(p_are_floats_equals(im(n1), im(n2)))))
}
.py_try({
tf1 <- TransferFunction(1L, s, s)
})
.py_try({
tf2 <- TransferFunction(s, .py_sub(.py_pow(s, 3L), 1L), s)
})
.py_try({
tf3 <- TransferFunction(s, .py_sub(s, 1L), s)
})
.py_try({
tf4 <- TransferFunction(s, .py_add(.py_pow(s, 2L), 1L), s)
})
.py_try({
tf5 <- TransferFunction(1L, 1L, s)
})
.py_try({
tf6 <- TransferFunction(-1L, .py_sub(s, 1L), s)
})
.py_try({
dtf1 <- DiscreteTransferFunction(1L, s, s, .py_flt(0.5))
})
.py_try({
dtf2 <- DiscreteTransferFunction(s, .py_sub(.py_pow(s, 3L), 1L), s, .py_flt(0.5))
})
.py_try({
dtf3 <- DiscreteTransferFunction(s, .py_sub(s, 1L), s, .py_flt(0.5))
})
.py_try({
dtf4 <- DiscreteTransferFunction(s, .py_add(.py_pow(s, 2L), 1L), s, .py_flt(0.5))
})
.py_try({
dtf5 <- DiscreteTransferFunction(1L, 1L, s, .py_flt(0.5))
})
.py_try({
dtf6 <- DiscreteTransferFunction(-1L, .py_sub(s, 1L), s, .py_flt(0.5))
})
.py_try({
tfm_1 <- TransferFunctionMatrix(list(list(tf1, tf2), list(tf3, tf4)))
})
.py_try({
tfm_2 <- TransferFunctionMatrix(list(list(tf2, tf3), list(tf4, tf1)))
})
.py_try({
tfm_3 <- .py_method(TransferFunctionMatrix, "from_Matrix", eye(2L), var = s)
})
.py_try({
tfm_4 <- TransferFunctionMatrix(list(list(tf1, tf5), list(tf5, tf5)))
})
.py_try({
tfm_5 <- TransferFunctionMatrix(list(list(.py_neg(tf3), tf3), list(tf3, tf6)))
})
.py_try({
tfm_6 <- TransferFunctionMatrix(list(list(.py_neg(tf3))))
})
.py_try({
tfm_7 <- TransferFunctionMatrix(list(list(tf3, tf4)))
})
.py_try({
dtfm_1 <- TransferFunctionMatrix(list(list(dtf1, dtf2), list(dtf3, dtf4)))
})
.py_try({
dtfm_2 <- TransferFunctionMatrix(list(list(dtf2, dtf3), list(dtf4, dtf1)))
})
.py_try({
dtfm_3 <- .py_method(TransferFunctionMatrix, "from_Matrix", eye(2L), var = s)
})
.py_try({
dtfm_4 <- TransferFunctionMatrix(list(list(dtf1, dtf5), list(dtf5, dtf5)))
})
.py_try({
dtfm_5 <- TransferFunctionMatrix(list(list(.py_neg(dtf3), dtf3), list(dtf3, dtf6)))
})
.py_try({
dtfm_6 <- TransferFunctionMatrix(list(list(.py_neg(dtf3))))
})
.py_try({
dtfm_7 <- TransferFunctionMatrix(list(list(dtf3, dtf4)))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_errors:2674", function() MIMOFeedback(tf1, tf2))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_errors:2675", function() MIMOFeedback(MIMOParallel(tfm_1, tfm_2), tfm_3))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_errors:2677", function() MIMOFeedback(tfm_1, tfm_6, 1L))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_errors:2678", function() MIMOFeedback(tfm_7, tfm_7))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_errors:2680", function() MIMOFeedback(tfm_1, tfm_2, -2L))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_errors:2682", function() MIMOFeedback(tfm_5, tfm_4, 1L))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_errors:2683", function() MIMOFeedback(tfm_4, .py_neg(tfm_5)))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_errors:2684", function() MIMOFeedback(tfm_3, tfm_3, 1L))
})
.py_try({
tfm_8 <- .py_method(TransferFunctionMatrix, "from_Matrix", eye(2L), var = p)
})
.py_try({
.py_rec("lti:test_MIMOFeedback_errors:2687", function() MIMOFeedback(tfm_1, tfm_8, 1L))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_errors:2691", function() MIMOFeedback(dtf1, dtf2))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_errors:2692", function() MIMOFeedback(MIMOParallel(dtfm_1, dtfm_2), dtfm_3))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_errors:2694", function() MIMOFeedback(dtfm_1, dtfm_6, 1L))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_errors:2695", function() MIMOFeedback(dtfm_7, dtfm_7))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_errors:2697", function() MIMOFeedback(dtfm_1, dtfm_2, -2L))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_errors:2699", function() MIMOFeedback(dtfm_5, dtfm_4, 1L))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_errors:2700", function() MIMOFeedback(dtfm_4, .py_neg(dtfm_5)))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_errors:2701", function() MIMOFeedback(dtfm_3, dtfm_3, 1L))
})
.py_try({
dtfm_8 <- .py_method(TransferFunctionMatrix, "from_Matrix", eye(2L), var = p)
})
.py_try({
.py_rec("lti:test_MIMOFeedback_errors:2704", function() MIMOFeedback(dtfm_1, dtfm_8, 1L))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_errors:2706", function() MIMOFeedback(tfm_1, dtfm_2))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_errors:2707", function() MIMOFeedback(dtfm_1, tfm_2))
}))---")
  expect_identical(got[["lti:test_MIMOFeedback_errors:2674"]], "ERR TypeError: Unsupported type for `sys1` in MIMO Feedback.")
  expect_identical(got[["lti:test_MIMOFeedback_errors:2675"]], "ERR TypeError: Unsupported type for `sys1` in MIMO Feedback.")
  expect_identical(got[["lti:test_MIMOFeedback_errors:2677"]], "ERR ValueError: \\nProduct of `sys1` and `sys2` must yield a square matrix.")
  expect_identical(got[["lti:test_MIMOFeedback_errors:2678"]], "ERR ValueError: \\nProduct of `sys1` and `sys2` must yield a square matrix.")
  expect_identical(got[["lti:test_MIMOFeedback_errors:2680"]], "ERR ValueError: \\nUnsupported type for feedback. `sign` arg should either be 1 (positive\\nfeedback loop) or -1 (negative feedback loop).")
  expect_identical(got[["lti:test_MIMOFeedback_errors:2682"]], "ERR ValueError: Non-Invertible system inputted.")
  expect_identical(got[["lti:test_MIMOFeedback_errors:2683"]], "ERR ValueError: Non-Invertible system inputted.")
  expect_identical(got[["lti:test_MIMOFeedback_errors:2684"]], "ERR ValueError: Non-Invertible system inputted.")
  expect_identical(got[["lti:test_MIMOFeedback_errors:2687"]], "ERR ValueError: \\nBoth `sys1` and `sys2` should be using the same complex variable.")
  expect_identical(got[["lti:test_MIMOFeedback_errors:2691"]], "ERR TypeError: Unsupported type for `sys1` in MIMO Feedback.")
  expect_identical(got[["lti:test_MIMOFeedback_errors:2692"]], "ERR TypeError: Unsupported type for `sys1` in MIMO Feedback.")
  expect_identical(got[["lti:test_MIMOFeedback_errors:2694"]], "ERR ValueError: \\nProduct of `sys1` and `sys2` must yield a square matrix.")
  expect_identical(got[["lti:test_MIMOFeedback_errors:2695"]], "ERR ValueError: \\nProduct of `sys1` and `sys2` must yield a square matrix.")
  expect_identical(got[["lti:test_MIMOFeedback_errors:2697"]], "ERR ValueError: \\nUnsupported type for feedback. `sign` arg should either be 1 (positive\\nfeedback loop) or -1 (negative feedback loop).")
  expect_identical(got[["lti:test_MIMOFeedback_errors:2699"]], "ERR ValueError: Non-Invertible system inputted.")
  expect_identical(got[["lti:test_MIMOFeedback_errors:2700"]], "ERR ValueError: Non-Invertible system inputted.")
  expect_identical(got[["lti:test_MIMOFeedback_errors:2701"]], "ERR ValueError: Non-Invertible system inputted.")
  expect_identical(got[["lti:test_MIMOFeedback_errors:2704"]], "ERR ValueError: \\nBoth `sys1` and `sys2` should be using the same complex variable.")
  expect_identical(got[["lti:test_MIMOFeedback_errors:2706"]], "ERR TypeError: \\n                Incompatible systems detected.\\n                All systems should be either continuous-time or discrete-time.\\n                Found: TransferFunctionMatrix(((TransferFunction(1, s, s), TransferFunction(s, s^3 - 1, s)), (TransferFunction(s, s - 1, s), TransferFunction(s, s^2 + 1, s)))) and TransferFunctionMatrix(((DiscreteTransferFunction(s, s^3 - 1, s, 0.5), DiscreteTransferFunction(s, s - 1, s, 0.5)), (DiscreteTransferFunction(s, s^2 + 1, s, 0.5), DiscreteTransferFunction(1, s, s, 0.5))))")
  expect_identical(got[["lti:test_MIMOFeedback_errors:2707"]], "ERR TypeError: \\n                Incompatible systems detected.\\n                All systems should be either continuous-time or discrete-time.\\n                Found: TransferFunctionMatrix(((DiscreteTransferFunction(1, s, s, 0.5), DiscreteTransferFunction(s, s^3 - 1, s, 0.5)), (DiscreteTransferFunction(s, s - 1, s, 0.5), DiscreteTransferFunction(s, s^2 + 1, s, 0.5)))) and TransferFunctionMatrix(((TransferFunction(s, s^3 - 1, s), TransferFunction(s, s - 1, s)), (TransferFunction(s, s^2 + 1, s), TransferFunction(1, s, s))))")
})

test_that("lti test_MIMOFeedback_functions matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("a, x, b, c, s, g, d,    p, k, tau, zeta, wn, T, z"); a <- .py_getitem(.tmpu_, 0L); x <- .py_getitem(.tmpu_, 1L); b <- .py_getitem(.tmpu_, 2L); c <- .py_getitem(.tmpu_, 3L); s <- .py_getitem(.tmpu_, 4L); g <- .py_getitem(.tmpu_, 5L); d <- .py_getitem(.tmpu_, 6L); p <- .py_getitem(.tmpu_, 7L); k <- .py_getitem(.tmpu_, 8L); tau <- .py_getitem(.tmpu_, 9L); zeta <- .py_getitem(.tmpu_, 10L); wn <- .py_getitem(.tmpu_, 11L); T <- .py_getitem(.tmpu_, 12L); z <- .py_getitem(.tmpu_, 13L)
})
.py_try({
.tmpu_ <- .py_symbols("a0:4, b0:5, c0:4, d0:4"); a0 <- .py_getitem(.tmpu_, 0L); a1 <- .py_getitem(.tmpu_, 1L); a2 <- .py_getitem(.tmpu_, 2L); a3 <- .py_getitem(.tmpu_, 3L); b0 <- .py_getitem(.tmpu_, 4L); b1 <- .py_getitem(.tmpu_, 5L); b2 <- .py_getitem(.tmpu_, 6L); b3 <- .py_getitem(.tmpu_, 7L); b4 <- .py_getitem(.tmpu_, 8L); c0 <- .py_getitem(.tmpu_, 9L); c1 <- .py_getitem(.tmpu_, 10L); c2 <- .py_getitem(.tmpu_, 11L); c3 <- .py_getitem(.tmpu_, 12L); d0 <- .py_getitem(.tmpu_, 13L); d1 <- .py_getitem(.tmpu_, 14L); d2 <- .py_getitem(.tmpu_, 15L); d3 <- .py_getitem(.tmpu_, 16L)
})
.py_try({
TF1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
TF2 <- TransferFunction(k, 1L, s)
})
.py_try({
TF3 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
TF4 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(tau, s), 1L), s)
})
.py_try({
TF5 <- TransferFunction(k, p, s)
})
.py_try({
TF6 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(k, s), p), s)
})
p_are_floats_equals <- function(n1, n2) {
  "\n    Helper function to compare two floats.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return(.py_method(Float(n1), "is_same", Float(n2), .py_isclose))
}
p_are_complex_equals <- function(n1, n2) {
  "\n    Helper function to compare two complex numbers.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return((.py_truthy(p_are_floats_equals(re(n1), re(n2))) && .py_truthy(p_are_floats_equals(im(n1), im(n2)))))
}
.py_try({
tf1 <- TransferFunction(1L, s, s)
})
.py_try({
tf2 <- TransferFunction(s, .py_sub(s, 1L), s)
})
.py_try({
tf3 <- TransferFunction(1L, 1L, s)
})
.py_try({
tf4 <- TransferFunction(-1L, .py_sub(s, 1L), s)
})
.py_try({
dtf1 <- DiscreteTransferFunction(1L, s, s, .py_flt(0.4))
})
.py_try({
dtf2 <- DiscreteTransferFunction(s, .py_sub(s, 1L), s, .py_flt(0.4))
})
.py_try({
dtf3 <- DiscreteTransferFunction(1L, 1L, s, .py_flt(0.4))
})
.py_try({
dtf4 <- DiscreteTransferFunction(-1L, .py_sub(s, 1L), s, .py_flt(0.4))
})
.py_try({
tfm_1 <- .py_method(TransferFunctionMatrix, "from_Matrix", eye(2L), var = s)
})
.py_try({
tfm_2 <- TransferFunctionMatrix(list(list(tf1, tf3), list(tf3, tf3)))
})
.py_try({
tfm_3 <- TransferFunctionMatrix(list(list(.py_neg(tf2), tf2), list(tf2, tf4)))
})
.py_try({
tfm_4 <- TransferFunctionMatrix(list(list(tf1, tf2), list(.py_neg(tf2), tf1)))
})
.py_try({
dtfm_1 <- .py_method(TransferFunctionMatrix, "from_Matrix", eye(2L), var = s, sampling_time = .py_flt(0.4))
})
.py_try({
dtfm_2 <- TransferFunctionMatrix(list(list(dtf1, dtf3), list(dtf3, dtf3)))
})
.py_try({
dtfm_3 <- TransferFunctionMatrix(list(list(.py_neg(dtf2), dtf2), list(dtf2, dtf4)))
})
.py_try({
dtfm_4 <- TransferFunctionMatrix(list(list(dtf1, dtf2), list(.py_neg(dtf2), dtf1)))
})
.py_try({
F_1 <- MIMOFeedback(tfm_2, tfm_3)
})
.py_try({
F_2 <- MIMOFeedback(tfm_2, MIMOSeries(tfm_4, .py_neg(tfm_1)), 1L)
})
.py_try({
.py_rec("lti:test_MIMOFeedback_functions:2736", function() .py_attr(F_1, "sensitivity"))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_functions:2737", function() .py_attr(F_2, "sensitivity"))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_functions:2740", function() (.py_truthy(.py_eq(.py_method(F_1, "doit"), TransferFunctionMatrix(.py_tup(.py_tup(TransferFunction(1L, .py_mul(2L, s), s), TransferFunction(1L, 2L, s)), .py_tup(TransferFunction(1L, 2L, s), TransferFunction(1L, 2L, s)))))) && .py_truthy(.py_eq(TransferFunctionMatrix(.py_tup(.py_tup(TransferFunction(1L, .py_mul(2L, s), s), TransferFunction(1L, 2L, s)), .py_tup(TransferFunction(1L, 2L, s), TransferFunction(1L, 2L, s)))), .py_method(F_1, "rewrite", TransferFunctionMatrix)))))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_functions:2744", function() .py_method(F_2, "doit", cancel = FALSE, expand = TRUE))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_functions:2747", function() .py_method(F_2, "doit", cancel = FALSE))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_functions:2751", function() .py_method(F_2, "doit"))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_functions:2754", function() .py_method(F_2, "doit", expand = TRUE))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_functions:2758", function() .py_neg(.py_method(F_1, "doit")))
})
.py_try({
dF_1 <- MIMOFeedback(dtfm_2, dtfm_3)
})
.py_try({
dF_2 <- MIMOFeedback(dtfm_2, MIMOSeries(dtfm_4, .py_neg(dtfm_1)), 1L)
})
.py_try({
.py_rec("lti:test_MIMOFeedback_functions:2764", function() .py_attr(dF_1, "sampling_time"))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_functions:2765", function() .py_attr(dF_2, "sampling_time"))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_functions:2767", function() .py_attr(dF_1, "sensitivity"))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_functions:2768", function() .py_attr(dF_2, "sensitivity"))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_functions:2771", function() (.py_truthy(.py_eq(.py_method(dF_1, "doit"), TransferFunctionMatrix(.py_tup(.py_tup(DiscreteTransferFunction(1L, .py_mul(2L, s), s, .py_flt(0.4)), DiscreteTransferFunction(1L, 2L, s, .py_flt(0.4))), .py_tup(DiscreteTransferFunction(1L, 2L, s, .py_flt(0.4)), DiscreteTransferFunction(1L, 2L, s, .py_flt(0.4))))))) && .py_truthy(.py_eq(TransferFunctionMatrix(.py_tup(.py_tup(DiscreteTransferFunction(1L, .py_mul(2L, s), s, .py_flt(0.4)), DiscreteTransferFunction(1L, 2L, s, .py_flt(0.4))), .py_tup(DiscreteTransferFunction(1L, 2L, s, .py_flt(0.4)), DiscreteTransferFunction(1L, 2L, s, .py_flt(0.4))))), .py_method(dF_1, "rewrite", TransferFunctionMatrix)))))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_functions:2778", function() .py_method(dF_2, "doit", cancel = FALSE, expand = TRUE))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_functions:2785", function() .py_method(dF_2, "doit", cancel = FALSE))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_functions:2793", function() .py_method(dF_2, "doit"))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_functions:2800", function() .py_method(dF_2, "doit", expand = TRUE))
})
.py_try({
.py_rec("lti:test_MIMOFeedback_functions:2807", function() .py_neg(.py_method(dF_1, "doit")))
}))---")
  expect_identical(got[["lti:test_MIMOFeedback_functions:2736"]], "Matrix([[1/2, 0], [0, 1/2]])")
  expect_identical(got[["lti:test_MIMOFeedback_functions:2737"]], "Matrix([[(-2*s^4 + s^2)/(s^2 - s + 1), (2*s^3 - s^2)/(s^2 - s + 1)], [-s^2, s]])")
  expect_identical(got[["lti:test_MIMOFeedback_functions:2740"]], "True")
  expect_identical(got[["lti:test_MIMOFeedback_functions:2744"]], "TransferFunctionMatrix(((TransferFunction(-s^5 + 2*s^4 - 2*s^3 + s^2, s^5 - 2*s^4 + 3*s^3 - 2*s^2 + s, s), TransferFunction(-2*s^4 + 2*s^3, s^2 - s + 1, s)), (TransferFunction(0, 1, s), TransferFunction(-s^2 + s, 1, s))))")
  expect_identical(got[["lti:test_MIMOFeedback_functions:2747"]], "TransferFunctionMatrix(((TransferFunction(s*(2*s^3 - s^2)*(s^2 - s + 1) + (-2*s^4 + s^2)*(s^2 - s + 1), s*(s^2 - s + 1)^2, s), TransferFunction(-2*s^4 + 2*s^3, s^2 - s + 1, s)), (TransferFunction(0, 1, s), TransferFunction(-s^2 + s, 1, s))))")
  expect_identical(got[["lti:test_MIMOFeedback_functions:2751"]], "TransferFunctionMatrix(((TransferFunction(s*(-2*s^2 + s*(2*s - 1) + 1), s^2 - s + 1, s), TransferFunction(-2*s^3*(s - 1), s^2 - s + 1, s)), (TransferFunction(0, 1, s), TransferFunction(s*(1 - s), 1, s))))")
  expect_identical(got[["lti:test_MIMOFeedback_functions:2754"]], "TransferFunctionMatrix(((TransferFunction(-s^2 + s, s^2 - s + 1, s), TransferFunction(-2*s^4 + 2*s^3, s^2 - s + 1, s)), (TransferFunction(0, 1, s), TransferFunction(-s^2 + s, 1, s))))")
  expect_identical(got[["lti:test_MIMOFeedback_functions:2758"]], "TransferFunctionMatrix(((TransferFunction(-1, 2*s, s), TransferFunction(-1, 2, s)), (TransferFunction(-1, 2, s), TransferFunction(-1, 2, s))))")
  expect_identical(got[["lti:test_MIMOFeedback_functions:2764"]], "0.400000000000000")
  expect_identical(got[["lti:test_MIMOFeedback_functions:2765"]], "0.400000000000000")
  expect_identical(got[["lti:test_MIMOFeedback_functions:2767"]], "Matrix([[1/2, 0], [0, 1/2]])")
  expect_identical(got[["lti:test_MIMOFeedback_functions:2768"]], "Matrix([[(-2*s^4 + s^2)/(s^2 - s + 1), (2*s^3 - s^2)/(s^2 - s + 1)], [-s^2, s]])")
  expect_identical(got[["lti:test_MIMOFeedback_functions:2771"]], "True")
  expect_identical(got[["lti:test_MIMOFeedback_functions:2778"]], "TransferFunctionMatrix(((DiscreteTransferFunction(-s^5 + 2*s^4 - 2*s^3 + s^2, s^5 - 2*s^4 + 3*s^3 - 2*s^2 + s, s, 0.4), DiscreteTransferFunction(-2*s^4 + 2*s^3, s^2 - s + 1, s, 0.4)), (DiscreteTransferFunction(0, 1, s, 0.4), DiscreteTransferFunction(-s^2 + s, 1, s, 0.4))))")
  expect_identical(got[["lti:test_MIMOFeedback_functions:2785"]], "TransferFunctionMatrix(((DiscreteTransferFunction(s*(2*s^3 - s^2)*(s^2 - s + 1) + (-2*s^4 + s^2)*(s^2 - s + 1), s*(s^2 - s + 1)^2, s, 0.4), DiscreteTransferFunction(-2*s^4 + 2*s^3, s^2 - s + 1, s, 0.4)), (DiscreteTransferFunction(0, 1, s, 0.4), DiscreteTransferFunction(-s^2 + s, 1, s, 0.4))))")
  expect_identical(got[["lti:test_MIMOFeedback_functions:2793"]], "TransferFunctionMatrix(((DiscreteTransferFunction(s*(-2*s^2 + s*(2*s - 1) + 1), s^2 - s + 1, s, 0.4), DiscreteTransferFunction(-2*s^3*(s - 1), s^2 - s + 1, s, 0.4)), (DiscreteTransferFunction(0, 1, s, 0.4), DiscreteTransferFunction(s*(1 - s), 1, s, 0.4))))")
  expect_identical(got[["lti:test_MIMOFeedback_functions:2800"]], "TransferFunctionMatrix(((DiscreteTransferFunction(-s^2 + s, s^2 - s + 1, s, 0.4), DiscreteTransferFunction(-2*s^4 + 2*s^3, s^2 - s + 1, s, 0.4)), (DiscreteTransferFunction(0, 1, s, 0.4), DiscreteTransferFunction(-s^2 + s, 1, s, 0.4))))")
  expect_identical(got[["lti:test_MIMOFeedback_functions:2807"]], "TransferFunctionMatrix(((DiscreteTransferFunction(-1, 2*s, s, 0.4), DiscreteTransferFunction(-1, 2, s, 0.4)), (DiscreteTransferFunction(-1, 2, s, 0.4), DiscreteTransferFunction(-1, 2, s, 0.4))))")
})

test_that("lti test_TransferFunctionMatrix_construction matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("a, x, b, c, s, g, d,    p, k, tau, zeta, wn, T, z"); a <- .py_getitem(.tmpu_, 0L); x <- .py_getitem(.tmpu_, 1L); b <- .py_getitem(.tmpu_, 2L); c <- .py_getitem(.tmpu_, 3L); s <- .py_getitem(.tmpu_, 4L); g <- .py_getitem(.tmpu_, 5L); d <- .py_getitem(.tmpu_, 6L); p <- .py_getitem(.tmpu_, 7L); k <- .py_getitem(.tmpu_, 8L); tau <- .py_getitem(.tmpu_, 9L); zeta <- .py_getitem(.tmpu_, 10L); wn <- .py_getitem(.tmpu_, 11L); T <- .py_getitem(.tmpu_, 12L); z <- .py_getitem(.tmpu_, 13L)
})
.py_try({
.tmpu_ <- .py_symbols("a0:4, b0:5, c0:4, d0:4"); a0 <- .py_getitem(.tmpu_, 0L); a1 <- .py_getitem(.tmpu_, 1L); a2 <- .py_getitem(.tmpu_, 2L); a3 <- .py_getitem(.tmpu_, 3L); b0 <- .py_getitem(.tmpu_, 4L); b1 <- .py_getitem(.tmpu_, 5L); b2 <- .py_getitem(.tmpu_, 6L); b3 <- .py_getitem(.tmpu_, 7L); b4 <- .py_getitem(.tmpu_, 8L); c0 <- .py_getitem(.tmpu_, 9L); c1 <- .py_getitem(.tmpu_, 10L); c2 <- .py_getitem(.tmpu_, 11L); c3 <- .py_getitem(.tmpu_, 12L); d0 <- .py_getitem(.tmpu_, 13L); d1 <- .py_getitem(.tmpu_, 14L); d2 <- .py_getitem(.tmpu_, 15L); d3 <- .py_getitem(.tmpu_, 16L)
})
.py_try({
TF1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
TF2 <- TransferFunction(k, 1L, s)
})
.py_try({
TF3 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
TF4 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(tau, s), 1L), s)
})
.py_try({
TF5 <- TransferFunction(k, p, s)
})
.py_try({
TF6 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(k, s), p), s)
})
p_are_floats_equals <- function(n1, n2) {
  "\n    Helper function to compare two floats.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return(.py_method(Float(n1), "is_same", Float(n2), .py_isclose))
}
p_are_complex_equals <- function(n1, n2) {
  "\n    Helper function to compare two complex numbers.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return((.py_truthy(p_are_floats_equals(re(n1), re(n2))) && .py_truthy(p_are_floats_equals(im(n1), im(n2)))))
}
.py_try({
tf5 <- TransferFunction(.py_sub(.py_add(.py_mul(a1, .py_pow(s, 2L)), .py_mul(a2, s)), a0), .py_add(s, a0), s)
})
.py_try({
tf4 <- TransferFunction(.py_sub(.py_add(.py_mul(a0, p), .py_pow(p, a1)), s), p, p)
})
.py_try({
dtf1 <- DiscreteTransferFunction(.py_sub(.py_add(.py_mul(a1, .py_pow(z, 2L)), .py_mul(a2, z)), a0), .py_add(z, a0), z, 2L)
})
.py_try({
dtf2 <- DiscreteTransferFunction(z, .py_sub(z, 1L), z, 2L)
})
.py_try({
dtf3 <- DiscreteTransferFunction(.py_sub(.py_add(.py_mul(a0, p), .py_pow(p, a1)), s), p, p, 3L)
})
.py_try({
tfm3_ <- TransferFunctionMatrix(list(list(.py_neg(TF3))))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_construction:2819", function() (.py_truthy(.py_eq(.py_attr(tfm3_, "shape"), .py_tup(.py_attr(tfm3_, "num_outputs"), .py_attr(tfm3_, "num_inputs")))) && .py_truthy(.py_eq(.py_tup(.py_attr(tfm3_, "num_outputs"), .py_attr(tfm3_, "num_inputs")), .py_tup(1L, 1L)))))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_construction:2820", function() .py_attr(tfm3_, "args"))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_construction:2821", function() .py_attr(tfm3_, "var"))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_construction:2822", function() .py_attr(tfm3_, "sampling_time"))
})
.py_try({
tfm5 <- TransferFunctionMatrix(list(list(TF1, .py_neg(TF2)), list(TF3, tf5)))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_construction:2825", function() (.py_truthy(.py_eq(.py_attr(tfm5, "shape"), .py_tup(.py_attr(tfm5, "num_outputs"), .py_attr(tfm5, "num_inputs")))) && .py_truthy(.py_eq(.py_tup(.py_attr(tfm5, "num_outputs"), .py_attr(tfm5, "num_inputs")), .py_tup(2L, 2L)))))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_construction:2826", function() .py_attr(tfm5, "args"))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_construction:2827", function() .py_attr(tfm5, "var"))
})
.py_try({
tfm7 <- TransferFunctionMatrix(list(list(TF1, TF2), list(TF3, .py_neg(tf5)), list(.py_neg(tf5), TF2)))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_construction:2830", function() (.py_truthy(.py_eq(.py_attr(tfm7, "shape"), .py_tup(.py_attr(tfm7, "num_outputs"), .py_attr(tfm7, "num_inputs")))) && .py_truthy(.py_eq(.py_tup(.py_attr(tfm7, "num_outputs"), .py_attr(tfm7, "num_inputs")), .py_tup(3L, 2L)))))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_construction:2831", function() .py_attr(tfm7, "args"))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_construction:2832", function() .py_attr(tfm7, "var"))
})
.py_try({
dtfm1 <- TransferFunctionMatrix(list(list(dtf1, dtf2), list(.py_add(.py_neg(dtf1), dtf2), .py_mul(dtf2, dtf1))))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_construction:2835", function() (.py_truthy(.py_eq(.py_attr(dtfm1, "shape"), .py_tup(.py_attr(dtfm1, "num_outputs"), .py_attr(dtfm1, "num_inputs")))) && .py_truthy(.py_eq(.py_tup(.py_attr(dtfm1, "num_outputs"), .py_attr(dtfm1, "num_inputs")), .py_tup(2L, 2L)))))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_construction:2836", function() .py_attr(dtfm1, "args"))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_construction:2839", function() .py_attr(dtfm1, "var"))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_construction:2840", function() .py_attr(dtfm1, "sampling_time"))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_construction:2842", function() .py_attr(tfm3_, "is_continuous"))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_construction:2843", function() .py_attr(dtfm1, "is_continuous"))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_construction:2846", function() TransferFunctionMatrix(list(list(TF1), list(TF2), list(tf4))))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_construction:2847", function() TransferFunctionMatrix(list(list(TF1, tf4), list(TF3, tf5))))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_construction:2850", function() TransferFunctionMatrix(list(list(TF1), list(TF3, tf5))))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_construction:2851", function() TransferFunctionMatrix(list(list(TF1, TF3), list(tf5))))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_construction:2854", function() TransferFunctionMatrix(list(list(TF1, TF2), list(TF3, .py_Matrix(list(1L, 2L))))))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_construction:2855", function() TransferFunctionMatrix(list(list(TF1, .py_Matrix(list(1L, 2L))), list(TF3, TF2))))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_construction:2858", function() TransferFunctionMatrix(list(TF1, TF2, tf5)))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_construction:2859", function() TransferFunctionMatrix(list(TF1)))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_construction:2862", function() TransferFunctionMatrix(list(list(TF1, dtf2), list(TF3, tf5))))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_construction:2865", function() TransferFunctionMatrix(list(list(dtf1, dtf2), list(dtf3, dtf1))))
}))---")
  expect_identical(got[["lti:test_TransferFunctionMatrix_construction:2819"]], "True")
  expect_identical(got[["lti:test_TransferFunctionMatrix_construction:2820"]], "(((TransferFunction(-a2*p + s, a2*s + p, s),),),)")
  expect_identical(got[["lti:test_TransferFunctionMatrix_construction:2821"]], "s")
  expect_identical(got[["lti:test_TransferFunctionMatrix_construction:2822"]], "0")
  expect_identical(got[["lti:test_TransferFunctionMatrix_construction:2825"]], "True")
  expect_identical(got[["lti:test_TransferFunctionMatrix_construction:2826"]], "(((TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s), TransferFunction(-k, 1, s)), (TransferFunction(a2*p - s, a2*s + p, s), TransferFunction(-a0 + a1*s^2 + a2*s, a0 + s, s))),)")
  expect_identical(got[["lti:test_TransferFunctionMatrix_construction:2827"]], "s")
  expect_identical(got[["lti:test_TransferFunctionMatrix_construction:2830"]], "True")
  expect_identical(got[["lti:test_TransferFunctionMatrix_construction:2831"]], "(((TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s), TransferFunction(k, 1, s)), (TransferFunction(a2*p - s, a2*s + p, s), TransferFunction(a0 - a1*s^2 - a2*s, a0 + s, s)), (TransferFunction(a0 - a1*s^2 - a2*s, a0 + s, s), TransferFunction(k, 1, s))),)")
  expect_identical(got[["lti:test_TransferFunctionMatrix_construction:2832"]], "s")
  expect_identical(got[["lti:test_TransferFunctionMatrix_construction:2835"]], "True")
  expect_identical(got[["lti:test_TransferFunctionMatrix_construction:2836"]], "(((DiscreteTransferFunction(-a0 + a1*z^2 + a2*z, a0 + z, z, 2), DiscreteTransferFunction(z, z - 1, z, 2)), (Parallel(DiscreteTransferFunction(a0 - a1*z^2 - a2*z, a0 + z, z, 2), DiscreteTransferFunction(z, z - 1, z, 2)), Series(DiscreteTransferFunction(z, z - 1, z, 2), DiscreteTransferFunction(-a0 + a1*z^2 + a2*z, a0 + z, z, 2)))),)")
  expect_identical(got[["lti:test_TransferFunctionMatrix_construction:2839"]], "z")
  expect_identical(got[["lti:test_TransferFunctionMatrix_construction:2840"]], "2")
  expect_identical(got[["lti:test_TransferFunctionMatrix_construction:2842"]], "True")
  expect_identical(got[["lti:test_TransferFunctionMatrix_construction:2843"]], "False")
  expect_identical(got[["lti:test_TransferFunctionMatrix_construction:2846"]], "ERR ValueError: \\nConflicting value(s) found for `var`. All TransferFunction instances\\nin TransferFunctionMatrix should use the same complex variable in\\nLaplace domain or z-domain.")
  expect_identical(got[["lti:test_TransferFunctionMatrix_construction:2847"]], "ERR ValueError: \\nConflicting value(s) found for `var`. All TransferFunction instances\\nin TransferFunctionMatrix should use the same complex variable in\\nLaplace domain or z-domain.")
  expect_identical(got[["lti:test_TransferFunctionMatrix_construction:2850"]], "ERR ValueError: mismatched dimensions")
  expect_identical(got[["lti:test_TransferFunctionMatrix_construction:2851"]], "ERR ValueError: mismatched dimensions")
  expect_identical(got[["lti:test_TransferFunctionMatrix_construction:2854"]], "ERR TypeError: \\nEach element is expected to be of type `SISOLinearTimeInvariant`.")
  expect_identical(got[["lti:test_TransferFunctionMatrix_construction:2855"]], "ERR TypeError: \\nEach element is expected to be of type `SISOLinearTimeInvariant`.")
  expect_identical(got[["lti:test_TransferFunctionMatrix_construction:2858"]], "ERR ValueError: \\n`arg` param in TransferFunctionMatrix should strictly be a nested list\\ncontaining TransferFunctionBase objects.")
  expect_identical(got[["lti:test_TransferFunctionMatrix_construction:2859"]], "ERR ValueError: \\n`arg` param in TransferFunctionMatrix should strictly be a nested list\\ncontaining TransferFunctionBase objects.")
  expect_identical(got[["lti:test_TransferFunctionMatrix_construction:2862"]], "ERR ValueError: \\nConflicting value(s) found for `var`. All TransferFunction instances\\nin TransferFunctionMatrix should use the same complex variable in\\nLaplace domain or z-domain.")
  expect_identical(got[["lti:test_TransferFunctionMatrix_construction:2865"]], "ERR ValueError: \\nConflicting value(s) found for `var`. All TransferFunction instances\\nin TransferFunctionMatrix should use the same complex variable in\\nLaplace domain or z-domain.")
})

test_that("lti test_TransferFunctionMatrix_functions matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("a, x, b, c, s, g, d,    p, k, tau, zeta, wn, T, z"); a <- .py_getitem(.tmpu_, 0L); x <- .py_getitem(.tmpu_, 1L); b <- .py_getitem(.tmpu_, 2L); c <- .py_getitem(.tmpu_, 3L); s <- .py_getitem(.tmpu_, 4L); g <- .py_getitem(.tmpu_, 5L); d <- .py_getitem(.tmpu_, 6L); p <- .py_getitem(.tmpu_, 7L); k <- .py_getitem(.tmpu_, 8L); tau <- .py_getitem(.tmpu_, 9L); zeta <- .py_getitem(.tmpu_, 10L); wn <- .py_getitem(.tmpu_, 11L); T <- .py_getitem(.tmpu_, 12L); z <- .py_getitem(.tmpu_, 13L)
})
.py_try({
.tmpu_ <- .py_symbols("a0:4, b0:5, c0:4, d0:4"); a0 <- .py_getitem(.tmpu_, 0L); a1 <- .py_getitem(.tmpu_, 1L); a2 <- .py_getitem(.tmpu_, 2L); a3 <- .py_getitem(.tmpu_, 3L); b0 <- .py_getitem(.tmpu_, 4L); b1 <- .py_getitem(.tmpu_, 5L); b2 <- .py_getitem(.tmpu_, 6L); b3 <- .py_getitem(.tmpu_, 7L); b4 <- .py_getitem(.tmpu_, 8L); c0 <- .py_getitem(.tmpu_, 9L); c1 <- .py_getitem(.tmpu_, 10L); c2 <- .py_getitem(.tmpu_, 11L); c3 <- .py_getitem(.tmpu_, 12L); d0 <- .py_getitem(.tmpu_, 13L); d1 <- .py_getitem(.tmpu_, 14L); d2 <- .py_getitem(.tmpu_, 15L); d3 <- .py_getitem(.tmpu_, 16L)
})
.py_try({
TF1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
TF2 <- TransferFunction(k, 1L, s)
})
.py_try({
TF3 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
TF4 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(tau, s), 1L), s)
})
.py_try({
TF5 <- TransferFunction(k, p, s)
})
.py_try({
TF6 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(k, s), p), s)
})
p_are_floats_equals <- function(n1, n2) {
  "\n    Helper function to compare two floats.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return(.py_method(Float(n1), "is_same", Float(n2), .py_isclose))
}
p_are_complex_equals <- function(n1, n2) {
  "\n    Helper function to compare two complex numbers.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return((.py_truthy(p_are_floats_equals(re(n1), re(n2))) && .py_truthy(p_are_floats_equals(im(n1), im(n2)))))
}
.py_try({
tf5 <- TransferFunction(.py_sub(.py_add(.py_mul(a1, .py_pow(s, 2L)), .py_mul(a2, s)), a0), .py_add(s, a0), s)
})
.py_try({
mat_1 <- .py_Matrix(list(list(.py_div(.py_mul(.py_mul(s, .py_add(s, 1L)), .py_sub(s, 3L)), .py_add(.py_pow(s, 4L), 1L)), 2L), list(p, .py_div(.py_mul(p, .py_add(s, 1L)), .py_mul(s, .py_add(.py_pow(s, 1L), 1L))))))
})
.py_try({
mat_2 <- .py_Matrix(list(list(.py_div(.py_add(.py_mul(2L, s), 1L), .py_sub(.py_pow(s, 2L), 9L)))))
})
.py_try({
mat_3 <- .py_Matrix(list(list(1L, 2L), list(3L, 4L)))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_functions:2879", function() .py_method(TransferFunctionMatrix, "from_Matrix", mat_1, s))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_functions:2882", function() .py_method(TransferFunctionMatrix, "from_Matrix", mat_2, s))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_functions:2884", function() .py_method(TransferFunctionMatrix, "from_Matrix", mat_3, p))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_functions:2888", function() .py_method(TransferFunctionMatrix, "from_Matrix", mat_1, s, 2L))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_functions:2894", function() .py_method(TransferFunctionMatrix, "from_Matrix", mat_2, s, .py_flt(0.1)))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_functions:2896", function() .py_method(TransferFunctionMatrix, "from_Matrix", mat_3, p, T))
})
.py_try({
tfm1 <- TransferFunctionMatrix(list(list(TF1), list(TF2)))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_functions:2904", function() .py_neg(tfm1))
})
.py_try({
tfm2 <- TransferFunctionMatrix(list(list(TF1, TF2, TF3), list(tf5, .py_neg(TF1), .py_neg(TF3))))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_functions:2907", function() .py_neg(tfm2))
})
.py_try({
H_1 <- .py_method(TransferFunctionMatrix, "from_Matrix", mat_1, s)
})
.py_try({
H_2 <- TransferFunctionMatrix(list(list(TransferFunction(.py_mul(.py_mul(a, p), s), .py_mul(k, .py_pow(s, 2L)), s), TransferFunction(.py_mul(p, s), .py_mul(k, .py_sub(.py_pow(s, 2L), a)), s))))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_functions:2913", function() .py_method(H_1, "subs", p, 1L))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_functions:2914", function() .py_method(H_1, "subs", .py_dict(list(p), list(1L))))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_functions:2915", function() .py_method(H_1, "subs", .py_dict(list(p, s), list(1L, 1L))))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_functions:2916", function() .py_method(H_2, "subs", p, 2L))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_functions:2917", function() .py_method(H_2, "subs", k, 1L))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_functions:2918", function() .py_method(H_2, "subs", a, 0L))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_functions:2919", function() .py_method(H_2, "subs", .py_dict(list(p, k, a), list(1L, 1L, a0))))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_functions:2922", function() .py_method(H_2, "eval_frequency", .py_add(.py_div(.py_S(1L), 2L), sym_I)))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_functions:2926", function() .py_method(H_1, "transpose"))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_functions:2927", function() .py_method(H_2, "transpose"))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_functions:2928", function() .py_method(.py_method(H_1, "transpose"), "transpose"))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_functions:2929", function() .py_method(.py_method(H_2, "transpose"), "transpose"))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_functions:2933", function() .py_method(H_1, "elem_poles"))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_functions:2935", function() .py_method(H_2, "elem_poles"))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_functions:2936", function() .py_method(tfm2, "elem_poles"))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_functions:2941", function() .py_method(H_1, "elem_zeros"))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_functions:2942", function() .py_method(H_2, "elem_zeros"))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_functions:2943", function() .py_method(tfm2, "elem_zeros"))
})
.py_try({
H_3 <- TransferFunctionMatrix(list(list(Series(TransferFunction(1L, .py_sub(.py_pow(s, 3L), 3L), s), TransferFunction(.py_add(.py_sub(.py_pow(s, 2L), .py_mul(2L, s)), 5L), 1L, s), TransferFunction(1L, s, s)))))
})
.py_try({
H_4 <- TransferFunctionMatrix(list(list(Parallel(TransferFunction(.py_sub(.py_pow(s, 3L), 3L), .py_add(.py_sub(.py_sub(.py_mul(4L, .py_pow(s, 4L)), .py_pow(s, 2L)), .py_mul(2L, s)), 5L), s), TransferFunction(.py_sub(4L, .py_pow(s, 3L)), .py_add(.py_sub(.py_sub(.py_mul(4L, .py_pow(s, 4L)), .py_pow(s, 2L)), .py_mul(2L, s)), 5L), s)))))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_functions:2951", function() .py_method(H_3, "doit"))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_functions:2952", function() .py_method(H_4, "doit"))
})
.py_try({
dH_3 <- TransferFunctionMatrix(list(list(Series(DiscreteTransferFunction(1L, .py_sub(.py_pow(s, 3L), 3L), s, .py_flt(0.1)), DiscreteTransferFunction(.py_add(.py_sub(.py_pow(s, 2L), .py_mul(2L, s)), 5L), 1L, s, .py_flt(0.1)), DiscreteTransferFunction(1L, s, s, .py_flt(0.1))))))
})
.py_try({
dH_4 <- TransferFunctionMatrix(list(list(Parallel(DiscreteTransferFunction(.py_sub(.py_pow(s, 3L), 3L), .py_add(.py_sub(.py_sub(.py_mul(4L, .py_pow(s, 4L)), .py_pow(s, 2L)), .py_mul(2L, s)), 5L), s, .py_flt(0.1)), DiscreteTransferFunction(.py_sub(4L, .py_pow(s, 3L)), .py_add(.py_sub(.py_sub(.py_mul(4L, .py_pow(s, 4L)), .py_pow(s, 2L)), .py_mul(2L, s)), 5L), s, .py_flt(0.1))))))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_functions:2963", function() .py_method(dH_3, "doit"))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_functions:2965", function() .py_method(dH_4, "doit"))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_functions:2970", function() .py_method(H_1, "_flat"))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_functions:2971", function() .py_method(H_2, "_flat"))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_functions:2972", function() .py_method(H_3, "_flat"))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_functions:2973", function() .py_method(H_4, "_flat"))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_functions:2977", function() .py_method(H_1, "evalf"))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_functions:2979", function() .py_method(.py_method(H_2, "subs", .py_dict(list(a, p, k), list(.py_flt(3.141), .py_flt(2.88), 2L))), "evalf"))
})
.py_try({
dH_1 <- .py_method(TransferFunctionMatrix, "from_Matrix", mat_1, s, T)
})
.py_try({
dH_2 <- TransferFunctionMatrix(list(list(DiscreteTransferFunction(.py_mul(.py_mul(a, p), s), .py_mul(k, .py_pow(s, 2L)), s, T), DiscreteTransferFunction(.py_mul(p, s), .py_mul(k, .py_sub(.py_pow(s, 2L), a)), s, T))))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_functions:2987", function() .py_method(dH_1, "evalf"))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_functions:2993", function() .py_method(.py_method(dH_2, "subs", .py_dict(list(a, p, k), list(.py_flt(3.141), .py_flt(2.88), 2L))), "evalf"))
})
.py_try({
H_5 <- TransferFunctionMatrix(list(list(TransferFunction(.py_add(.py_add(.py_pow(s, 5L), .py_pow(s, 3L)), s), .py_sub(s, .py_pow(s, 2L)), s), TransferFunction(.py_mul(.py_add(s, 3L), .py_sub(s, 1L)), .py_mul(.py_sub(s, 1L), .py_add(s, 5L)), s))))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_functions:3001", function() (.py_truthy(.py_eq(.py_method(H_5, "simplify"), simplify(H_5))) && .py_truthy(.py_eq(simplify(H_5), TransferFunctionMatrix(.py_tup(.py_tup(TransferFunction(.py_sub(.py_sub(.py_neg(.py_pow(s, 4L)), .py_pow(s, 2L)), 1L), .py_sub(s, 1L), s), TransferFunction(.py_add(s, 3L), .py_add(s, 5L), s))))))))
})
.py_try({
dH_5 <- TransferFunctionMatrix(list(list(DiscreteTransferFunction(.py_add(.py_add(.py_pow(s, 5L), .py_pow(s, 3L)), s), .py_sub(s, .py_pow(s, 2L)), s, 2L), DiscreteTransferFunction(.py_mul(.py_add(s, 3L), .py_sub(s, 1L)), .py_mul(.py_sub(s, 1L), .py_add(s, 5L)), s, 2L))))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_functions:3010", function() (.py_truthy(.py_eq(.py_method(dH_5, "simplify"), simplify(dH_5))) && .py_truthy(.py_eq(simplify(dH_5), TransferFunctionMatrix(.py_tup(.py_tup(DiscreteTransferFunction(.py_sub(.py_sub(.py_neg(.py_pow(s, 4L)), .py_pow(s, 2L)), 1L), .py_sub(s, 1L), s, 2L), DiscreteTransferFunction(.py_add(s, 3L), .py_add(s, 5L), s, 2L))))))))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_functions:3017", function() .py_method(H_1, "expand"))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_functions:3020", function() .py_method(H_5, "expand"))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_functions:3025", function() .py_method(dH_1, "expand"))
})
.py_try({
.py_rec("lti:test_TransferFunctionMatrix_functions:3031", function() .py_method(dH_5, "expand"))
}))---")
  expect_identical(got[["lti:test_TransferFunctionMatrix_functions:2879"]], "TransferFunctionMatrix(((TransferFunction(s*(s - 3)*(s + 1), s^4 + 1, s), TransferFunction(2, 1, s)), (TransferFunction(p, 1, s), TransferFunction(p, s, s))))")
  expect_identical(got[["lti:test_TransferFunctionMatrix_functions:2882"]], "TransferFunctionMatrix(((TransferFunction(2*s + 1, s^2 - 9, s),),))")
  expect_identical(got[["lti:test_TransferFunctionMatrix_functions:2884"]], "TransferFunctionMatrix(((TransferFunction(1, 1, p), TransferFunction(2, 1, p)), (TransferFunction(3, 1, p), TransferFunction(4, 1, p))))")
  expect_identical(got[["lti:test_TransferFunctionMatrix_functions:2888"]], "TransferFunctionMatrix(((DiscreteTransferFunction(s*(s - 3)*(s + 1), s^4 + 1, s, 2), DiscreteTransferFunction(2, 1, s, 2)), (DiscreteTransferFunction(p, 1, s, 2), DiscreteTransferFunction(p, s, s, 2))))")
  expect_identical(got[["lti:test_TransferFunctionMatrix_functions:2894"]], "TransferFunctionMatrix(((DiscreteTransferFunction(2*s + 1, s^2 - 9, s, 0.1),),))")
  expect_identical(got[["lti:test_TransferFunctionMatrix_functions:2896"]], "TransferFunctionMatrix(((DiscreteTransferFunction(1, 1, p, T), DiscreteTransferFunction(2, 1, p, T)), (DiscreteTransferFunction(3, 1, p, T), DiscreteTransferFunction(4, 1, p, T))))")
  expect_identical(got[["lti:test_TransferFunctionMatrix_functions:2904"]], "TransferFunctionMatrix(((TransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s),), (TransferFunction(-k, 1, s),)))")
  expect_identical(got[["lti:test_TransferFunctionMatrix_functions:2907"]], "TransferFunctionMatrix(((TransferFunction(-1, s^2 + 2*s*wn*zeta + wn^2, s), TransferFunction(-k, 1, s), TransferFunction(-a2*p + s, a2*s + p, s)), (TransferFunction(a0 - a1*s^2 - a2*s, a0 + s, s), TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s), TransferFunction(a2*p - s, a2*s + p, s))))")
  expect_identical(got[["lti:test_TransferFunctionMatrix_functions:2913"]], "TransferFunctionMatrix(((TransferFunction(s*(s - 3)*(s + 1), s^4 + 1, s), TransferFunction(2, 1, s)), (TransferFunction(1, 1, s), TransferFunction(1, s, s))))")
  expect_identical(got[["lti:test_TransferFunctionMatrix_functions:2914"]], "TransferFunctionMatrix(((TransferFunction(s*(s - 3)*(s + 1), s^4 + 1, s), TransferFunction(2, 1, s)), (TransferFunction(1, 1, s), TransferFunction(1, s, s))))")
  expect_identical(got[["lti:test_TransferFunctionMatrix_functions:2915"]], "TransferFunctionMatrix(((TransferFunction(s*(s - 3)*(s + 1), s^4 + 1, s), TransferFunction(2, 1, s)), (TransferFunction(1, 1, s), TransferFunction(1, s, s))))")
  expect_identical(got[["lti:test_TransferFunctionMatrix_functions:2916"]], "TransferFunctionMatrix(((TransferFunction(2*a*s, k*s^2, s), TransferFunction(2*s, k*(-a + s^2), s)),))")
  expect_identical(got[["lti:test_TransferFunctionMatrix_functions:2917"]], "TransferFunctionMatrix(((TransferFunction(a*p*s, s^2, s), TransferFunction(p*s, -a + s^2, s)),))")
  expect_identical(got[["lti:test_TransferFunctionMatrix_functions:2918"]], "TransferFunctionMatrix(((TransferFunction(0, k*s^2, s), TransferFunction(p*s, k*s^2, s)),))")
  expect_identical(got[["lti:test_TransferFunctionMatrix_functions:2919"]], "TransferFunctionMatrix(((TransferFunction(a0*s, s^2, s), TransferFunction(s, -a0 + s^2, s)),))")
  expect_identical(got[["lti:test_TransferFunctionMatrix_functions:2922"]], "Matrix([[2*a*p/(5*k) - 4*I*a*p/(5*k), I*p/(-a*k - 3*k/4 + I*k) + p/(-2*a*k - 3*k/2 + 2*I*k)]])")
  expect_identical(got[["lti:test_TransferFunctionMatrix_functions:2926"]], "TransferFunctionMatrix(((TransferFunction(s*(s - 3)*(s + 1), s^4 + 1, s), TransferFunction(p, 1, s)), (TransferFunction(2, 1, s), TransferFunction(p, s, s))))")
  expect_identical(got[["lti:test_TransferFunctionMatrix_functions:2927"]], "TransferFunctionMatrix(((TransferFunction(a*p*s, k*s^2, s),), (TransferFunction(p*s, k*(-a + s^2), s),)))")
  expect_identical(got[["lti:test_TransferFunctionMatrix_functions:2928"]], "TransferFunctionMatrix(((TransferFunction(s*(s - 3)*(s + 1), s^4 + 1, s), TransferFunction(2, 1, s)), (TransferFunction(p, 1, s), TransferFunction(p, s, s))))")
  expect_identical(got[["lti:test_TransferFunctionMatrix_functions:2929"]], "TransferFunctionMatrix(((TransferFunction(a*p*s, k*s^2, s), TransferFunction(p*s, k*(-a + s^2), s)),))")
  expect_identical(got[["lti:test_TransferFunctionMatrix_functions:2933"]], "[[[-sqrt(2)/2 - sqrt(2)*I/2, -sqrt(2)/2 + sqrt(2)*I/2, sqrt(2)/2 - sqrt(2)*I/2, sqrt(2)/2 + sqrt(2)*I/2], []], [[], [0]]]")
  expect_identical(got[["lti:test_TransferFunctionMatrix_functions:2935"]], "[[[0, 0], [sqrt(a), -sqrt(a)]]]")
  expect_identical(got[["lti:test_TransferFunctionMatrix_functions:2936"]], "[[[wn*(-zeta + sqrt((zeta - 1)*(zeta + 1))), wn*(-zeta - sqrt((zeta - 1)*(zeta + 1)))], [], [-p/a2]], [[-a0], [wn*(-zeta + sqrt((zeta - 1)*(zeta + 1))), wn*(-zeta - sqrt((zeta - 1)*(zeta + 1)))], [-p/a2]]]")
  expect_identical(got[["lti:test_TransferFunctionMatrix_functions:2941"]], "[[[-1, 0, 3], []], [[], []]]")
  expect_identical(got[["lti:test_TransferFunctionMatrix_functions:2942"]], "[[[0], [0]]]")
  expect_identical(got[["lti:test_TransferFunctionMatrix_functions:2943"]], "[[[], [], [a2*p]], [[-a2/(2*a1) - sqrt(4*a0*a1 + a2^2)/(2*a1), -a2/(2*a1) + sqrt(4*a0*a1 + a2^2)/(2*a1)], [], [a2*p]]]")
  expect_identical(got[["lti:test_TransferFunctionMatrix_functions:2951"]], "TransferFunctionMatrix(((TransferFunction(s^2 - 2*s + 5, s*(s^3 - 3), s),),))")
  expect_identical(got[["lti:test_TransferFunctionMatrix_functions:2952"]], "TransferFunctionMatrix(((TransferFunction(1, 4*s^4 - s^2 - 2*s + 5, s),),))")
  expect_identical(got[["lti:test_TransferFunctionMatrix_functions:2963"]], "TransferFunctionMatrix(((DiscreteTransferFunction(s^2 - 2*s + 5, s*(s^3 - 3), s, 0.1),),))")
  expect_identical(got[["lti:test_TransferFunctionMatrix_functions:2965"]], "TransferFunctionMatrix(((DiscreteTransferFunction(1, 4*s^4 - s^2 - 2*s + 5, s, 0.1),),))")
  expect_identical(got[["lti:test_TransferFunctionMatrix_functions:2970"]], "[TransferFunction(s*(s - 3)*(s + 1), s^4 + 1, s), TransferFunction(2, 1, s), TransferFunction(p, 1, s), TransferFunction(p, s, s)]")
  expect_identical(got[["lti:test_TransferFunctionMatrix_functions:2971"]], "[TransferFunction(a*p*s, k*s^2, s), TransferFunction(p*s, k*(-a + s^2), s)]")
  expect_identical(got[["lti:test_TransferFunctionMatrix_functions:2972"]], "[Series(TransferFunction(1, s^3 - 3, s), TransferFunction(s^2 - 2*s + 5, 1, s), TransferFunction(1, s, s))]")
  expect_identical(got[["lti:test_TransferFunctionMatrix_functions:2973"]], "[Parallel(TransferFunction(s^3 - 3, 4*s^4 - s^2 - 2*s + 5, s), TransferFunction(4 - s^3, 4*s^4 - s^2 - 2*s + 5, s))]")
  expect_identical(got[["lti:test_TransferFunctionMatrix_functions:2977"]], "TransferFunctionMatrix(((TransferFunction(s*(s - 3.0)*(s + 1.0), s^4 + 1.0, s), TransferFunction(2.0, 1, s)), (TransferFunction(1.0*p, 1, s), TransferFunction(p, s, s))))")
  expect_identical(got[["lti:test_TransferFunctionMatrix_functions:2979"]], "TransferFunctionMatrix(((TransferFunction(4.52304, s, s), TransferFunction(2.88*s, 2.0*s^2 - 6.282, s)),))")
  expect_identical(got[["lti:test_TransferFunctionMatrix_functions:2987"]], "TransferFunctionMatrix(((DiscreteTransferFunction(s*(s - 3.0)*(s + 1.0), s^4 + 1.0, s, T), DiscreteTransferFunction(2.0, 1, s, T)), (DiscreteTransferFunction(1.0*p, 1, s, T), DiscreteTransferFunction(p, s, s, T))))")
  expect_identical(got[["lti:test_TransferFunctionMatrix_functions:2993"]], "TransferFunctionMatrix(((DiscreteTransferFunction(4.52304, s, s, T), DiscreteTransferFunction(2.88*s, 2.0*s^2 - 6.282, s, T)),))")
  expect_identical(got[["lti:test_TransferFunctionMatrix_functions:3001"]], "True")
  expect_identical(got[["lti:test_TransferFunctionMatrix_functions:3010"]], "True")
  expect_identical(got[["lti:test_TransferFunctionMatrix_functions:3017"]], "TransferFunctionMatrix(((TransferFunction(s^3 - 2*s^2 - 3*s, s^4 + 1, s), TransferFunction(2, 1, s)), (TransferFunction(p, 1, s), TransferFunction(p, s, s))))")
  expect_identical(got[["lti:test_TransferFunctionMatrix_functions:3020"]], "TransferFunctionMatrix(((TransferFunction(s^5 + s^3 + s, -s^2 + s, s), TransferFunction(s^2 + 2*s - 3, s^2 + 4*s - 5, s)),))")
  expect_identical(got[["lti:test_TransferFunctionMatrix_functions:3025"]], "TransferFunctionMatrix(((DiscreteTransferFunction(s^3 - 2*s^2 - 3*s, s^4 + 1, s, T), DiscreteTransferFunction(2, 1, s, T)), (DiscreteTransferFunction(p, 1, s, T), DiscreteTransferFunction(p, s, s, T))))")
  expect_identical(got[["lti:test_TransferFunctionMatrix_functions:3031"]], "TransferFunctionMatrix(((DiscreteTransferFunction(s^5 + s^3 + s, -s^2 + s, s, 2), DiscreteTransferFunction(s^2 + 2*s - 3, s^2 + 4*s - 5, s, 2)),))")
})

test_that("lti test_TransferFunction_gbt matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("a, x, b, c, s, g, d,    p, k, tau, zeta, wn, T, z"); a <- .py_getitem(.tmpu_, 0L); x <- .py_getitem(.tmpu_, 1L); b <- .py_getitem(.tmpu_, 2L); c <- .py_getitem(.tmpu_, 3L); s <- .py_getitem(.tmpu_, 4L); g <- .py_getitem(.tmpu_, 5L); d <- .py_getitem(.tmpu_, 6L); p <- .py_getitem(.tmpu_, 7L); k <- .py_getitem(.tmpu_, 8L); tau <- .py_getitem(.tmpu_, 9L); zeta <- .py_getitem(.tmpu_, 10L); wn <- .py_getitem(.tmpu_, 11L); T <- .py_getitem(.tmpu_, 12L); z <- .py_getitem(.tmpu_, 13L)
})
.py_try({
.tmpu_ <- .py_symbols("a0:4, b0:5, c0:4, d0:4"); a0 <- .py_getitem(.tmpu_, 0L); a1 <- .py_getitem(.tmpu_, 1L); a2 <- .py_getitem(.tmpu_, 2L); a3 <- .py_getitem(.tmpu_, 3L); b0 <- .py_getitem(.tmpu_, 4L); b1 <- .py_getitem(.tmpu_, 5L); b2 <- .py_getitem(.tmpu_, 6L); b3 <- .py_getitem(.tmpu_, 7L); b4 <- .py_getitem(.tmpu_, 8L); c0 <- .py_getitem(.tmpu_, 9L); c1 <- .py_getitem(.tmpu_, 10L); c2 <- .py_getitem(.tmpu_, 11L); c3 <- .py_getitem(.tmpu_, 12L); d0 <- .py_getitem(.tmpu_, 13L); d1 <- .py_getitem(.tmpu_, 14L); d2 <- .py_getitem(.tmpu_, 15L); d3 <- .py_getitem(.tmpu_, 16L)
})
.py_try({
TF1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
TF2 <- TransferFunction(k, 1L, s)
})
.py_try({
TF3 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
TF4 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(tau, s), 1L), s)
})
.py_try({
TF5 <- TransferFunction(k, p, s)
})
.py_try({
TF6 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(k, s), p), s)
})
p_are_floats_equals <- function(n1, n2) {
  "\n    Helper function to compare two floats.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return(.py_method(Float(n1), "is_same", Float(n2), .py_isclose))
}
p_are_complex_equals <- function(n1, n2) {
  "\n    Helper function to compare two complex numbers.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return((.py_truthy(p_are_floats_equals(re(n1), re(n2))) && .py_truthy(p_are_floats_equals(im(n1), im(n2)))))
}
.py_try({
tf <- TransferFunction(1L, .py_add(.py_mul(a, s), b), s)
})
.py_try({
.tmpu_ <- gbt(tf, T, .py_flt(0.5)); numZ <- .py_getitem(.tmpu_, 0L); denZ <- .py_getitem(.tmpu_, 1L)
})
.py_try({
dtf_test_bilinear1 <- .py_method(.py_method(DiscreteTransferFunction, "from_coeff_lists", numZ, denZ, z, T), "expand")
})
.py_try({
dtf_test_bilinear2 <- .py_method(.py_method(DiscreteTransferFunction, "from_gbt", tf, T, .py_flt(0.5), z), "expand")
})
.py_try({
num_manual <- list(.py_div(T, .py_mul(2L, .py_add(a, .py_div(.py_mul(b, T), 2L)))), .py_div(T, .py_mul(2L, .py_add(a, .py_div(.py_mul(b, T), 2L)))))
})
.py_try({
den_manual <- list(1L, .py_div(.py_add(.py_neg(a), .py_div(.py_mul(b, T), 2L)), .py_add(a, .py_div(.py_mul(b, T), 2L))))
})
.py_try({
dtf_manual <- .py_method(.py_method(DiscreteTransferFunction, "from_coeff_lists", num_manual, den_manual, z, T), "expand")
})
.py_try({
num_manual <- list(.py_div(T, .py_mul(2L, .py_add(a, .py_div(.py_mul(b, T), 2L)))), .py_div(T, .py_mul(2L, .py_add(a, .py_div(.py_mul(b, T), 2L)))))
})
.py_try({
den_manual <- list(1L, .py_div(.py_add(.py_neg(a), .py_div(.py_mul(b, T), 2L)), .py_add(a, .py_div(.py_mul(b, T), 2L))))
})
.py_try({
dtf_manual <- .py_method(.py_method(DiscreteTransferFunction, "from_coeff_lists", num_manual, den_manual, z, T), "expand")
})
.py_try({
.py_rec("lti:test_TransferFunction_gbt:3056", function() (.py_truthy(.py_eq(dtf_test_bilinear1, dtf_test_bilinear2)) && .py_truthy(.py_eq(dtf_test_bilinear2, dtf_manual))))
})
.py_try({
.py_rec("lti:test_TransferFunction_gbt:3057", function() (.py_truthy(.py_eq(dtf_test_bilinear1, dtf_test_bilinear2)) && .py_truthy(.py_eq(dtf_test_bilinear2, dtf_manual))))
})
.py_try({
tf <- TransferFunction(1L, .py_add(.py_mul(a, s), b), s)
})
.py_try({
.tmpu_ <- gbt(tf, T, 0L); numZ <- .py_getitem(.tmpu_, 0L); denZ <- .py_getitem(.tmpu_, 1L)
})
.py_try({
dtf_test_forward1 <- .py_method(.py_method(DiscreteTransferFunction, "from_coeff_lists", numZ, denZ, z, T), "expand")
})
.py_try({
dtf_test_forward2 <- .py_method(.py_method(DiscreteTransferFunction, "from_gbt", tf, T, 0L, z), "expand")
})
.py_try({
num_manual <- list(.py_div(T, a))
})
.py_try({
den_manual <- list(1L, .py_div(.py_add(.py_neg(a), .py_mul(b, T)), a))
})
.py_try({
dtf_manual <- .py_method(.py_method(DiscreteTransferFunction, "from_coeff_lists", num_manual, den_manual, z, T), "expand")
})
.py_try({
.py_rec("lti:test_TransferFunction_gbt:3073", function() (.py_truthy(.py_eq(dtf_test_forward1, dtf_test_forward2)) && .py_truthy(.py_eq(dtf_test_forward2, dtf_manual))))
})
.py_try({
tf <- TransferFunction(1L, .py_add(.py_mul(a, s), b), s)
})
.py_try({
.tmpu_ <- gbt(tf, T, 1L); numZ <- .py_getitem(.tmpu_, 0L); denZ <- .py_getitem(.tmpu_, 1L)
})
.py_try({
dtf_test_backward1 <- .py_method(.py_method(DiscreteTransferFunction, "from_coeff_lists", numZ, denZ, z, T), "expand")
})
.py_try({
dtf_test_backward2 <- .py_method(.py_method(DiscreteTransferFunction, "from_gbt", tf, T, 1L, z), "expand")
})
.py_try({
num_manual <- list(.py_div(T, .py_add(a, .py_mul(b, T))), 0L)
})
.py_try({
den_manual <- list(1L, .py_div(.py_neg(a), .py_add(a, .py_mul(b, T))))
})
.py_try({
dtf_manual <- .py_method(.py_method(DiscreteTransferFunction, "from_coeff_lists", num_manual, den_manual, z, T), "expand")
})
.py_try({
.py_rec("lti:test_TransferFunction_gbt:3088", function() (.py_truthy(.py_eq(dtf_test_backward1, dtf_test_backward2)) && .py_truthy(.py_eq(dtf_test_backward2, dtf_manual))))
})
.py_try({
tf <- TransferFunction(1L, .py_add(.py_mul(a, s), b), s)
})
.py_try({
.tmpu_ <- gbt(tf, T, .py_flt(0.3)); numZ <- .py_getitem(.tmpu_, 0L); denZ <- .py_getitem(.tmpu_, 1L)
})
.py_try({
dtf_test_gbt1 <- .py_method(.py_method(DiscreteTransferFunction, "from_coeff_lists", numZ, denZ, z, T), "expand")
})
.py_try({
dtf_test_gbt2 <- .py_method(.py_method(DiscreteTransferFunction, "from_gbt", tf, T, .py_flt(0.3), z), "expand")
})
.py_try({
num_manual <- list(.py_div(.py_mul(3L, T), .py_mul(10L, .py_add(a, .py_div(.py_mul(.py_mul(3L, b), T), 10L)))), .py_div(.py_mul(7L, T), .py_mul(10L, .py_add(a, .py_div(.py_mul(.py_mul(3L, b), T), 10L)))))
})
.py_try({
den_manual <- list(1L, .py_div(.py_add(.py_neg(a), .py_div(.py_mul(.py_mul(7L, b), T), 10L)), .py_add(a, .py_div(.py_mul(.py_mul(3L, b), T), 10L))))
})
.py_try({
dtf_manual <- .py_method(.py_method(DiscreteTransferFunction, "from_coeff_lists", num_manual, den_manual, z, T), "expand")
})
.py_try({
.py_rec("lti:test_TransferFunction_gbt:3103", function() (.py_truthy(.py_eq(dtf_test_gbt1, dtf_test_gbt2)) && .py_truthy(.py_eq(dtf_test_gbt2, dtf_manual))))
}))---")
  expect_identical(got[["lti:test_TransferFunction_gbt:3056"]], "True")
  expect_identical(got[["lti:test_TransferFunction_gbt:3057"]], "True")
  expect_identical(got[["lti:test_TransferFunction_gbt:3073"]], "True")
  expect_identical(got[["lti:test_TransferFunction_gbt:3088"]], "True")
  expect_identical(got[["lti:test_TransferFunction_gbt:3103"]], "True")
})

test_that("lti test_TransferFunction_bilinear matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("a, x, b, c, s, g, d,    p, k, tau, zeta, wn, T, z"); a <- .py_getitem(.tmpu_, 0L); x <- .py_getitem(.tmpu_, 1L); b <- .py_getitem(.tmpu_, 2L); c <- .py_getitem(.tmpu_, 3L); s <- .py_getitem(.tmpu_, 4L); g <- .py_getitem(.tmpu_, 5L); d <- .py_getitem(.tmpu_, 6L); p <- .py_getitem(.tmpu_, 7L); k <- .py_getitem(.tmpu_, 8L); tau <- .py_getitem(.tmpu_, 9L); zeta <- .py_getitem(.tmpu_, 10L); wn <- .py_getitem(.tmpu_, 11L); T <- .py_getitem(.tmpu_, 12L); z <- .py_getitem(.tmpu_, 13L)
})
.py_try({
.tmpu_ <- .py_symbols("a0:4, b0:5, c0:4, d0:4"); a0 <- .py_getitem(.tmpu_, 0L); a1 <- .py_getitem(.tmpu_, 1L); a2 <- .py_getitem(.tmpu_, 2L); a3 <- .py_getitem(.tmpu_, 3L); b0 <- .py_getitem(.tmpu_, 4L); b1 <- .py_getitem(.tmpu_, 5L); b2 <- .py_getitem(.tmpu_, 6L); b3 <- .py_getitem(.tmpu_, 7L); b4 <- .py_getitem(.tmpu_, 8L); c0 <- .py_getitem(.tmpu_, 9L); c1 <- .py_getitem(.tmpu_, 10L); c2 <- .py_getitem(.tmpu_, 11L); c3 <- .py_getitem(.tmpu_, 12L); d0 <- .py_getitem(.tmpu_, 13L); d1 <- .py_getitem(.tmpu_, 14L); d2 <- .py_getitem(.tmpu_, 15L); d3 <- .py_getitem(.tmpu_, 16L)
})
.py_try({
TF1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
TF2 <- TransferFunction(k, 1L, s)
})
.py_try({
TF3 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
TF4 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(tau, s), 1L), s)
})
.py_try({
TF5 <- TransferFunction(k, p, s)
})
.py_try({
TF6 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(k, s), p), s)
})
p_are_floats_equals <- function(n1, n2) {
  "\n    Helper function to compare two floats.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return(.py_method(Float(n1), "is_same", Float(n2), .py_isclose))
}
p_are_complex_equals <- function(n1, n2) {
  "\n    Helper function to compare two complex numbers.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return((.py_truthy(p_are_floats_equals(re(n1), re(n2))) && .py_truthy(p_are_floats_equals(im(n1), im(n2)))))
}
.py_try({
tf <- TransferFunction(1L, .py_add(.py_mul(a, s), b), s)
})
.py_try({
.tmpu_ <- bilinear(tf, T); numZ <- .py_getitem(.tmpu_, 0L); denZ <- .py_getitem(.tmpu_, 1L)
})
.py_try({
dtf_test_bilinear1 <- .py_method(.py_method(DiscreteTransferFunction, "from_coeff_lists", numZ, denZ, z, T), "expand")
})
.py_try({
dtf_test_bilinear2 <- .py_method(.py_method(DiscreteTransferFunction, "from_bilinear", tf, T, z), "expand")
})
.py_try({
num_manual <- list(.py_div(T, .py_mul(2L, .py_add(a, .py_div(.py_mul(b, T), 2L)))), .py_div(T, .py_mul(2L, .py_add(a, .py_div(.py_mul(b, T), 2L)))))
})
.py_try({
den_manual <- list(1L, .py_div(.py_add(.py_neg(a), .py_div(.py_mul(b, T), 2L)), .py_add(a, .py_div(.py_mul(b, T), 2L))))
})
.py_try({
dtf_manual <- .py_method(.py_method(DiscreteTransferFunction, "from_coeff_lists", num_manual, den_manual, z, T), "expand")
})
.py_try({
.py_rec("lti:test_TransferFunction_bilinear:3121", function() (.py_truthy(.py_eq(dtf_test_bilinear1, dtf_test_bilinear2)) && .py_truthy(.py_eq(dtf_test_bilinear2, dtf_manual))))
}))---")
  expect_identical(got[["lti:test_TransferFunction_bilinear:3121"]], "True")
})

test_that("lti test_TransferFunction_forward_diff matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("a, x, b, c, s, g, d,    p, k, tau, zeta, wn, T, z"); a <- .py_getitem(.tmpu_, 0L); x <- .py_getitem(.tmpu_, 1L); b <- .py_getitem(.tmpu_, 2L); c <- .py_getitem(.tmpu_, 3L); s <- .py_getitem(.tmpu_, 4L); g <- .py_getitem(.tmpu_, 5L); d <- .py_getitem(.tmpu_, 6L); p <- .py_getitem(.tmpu_, 7L); k <- .py_getitem(.tmpu_, 8L); tau <- .py_getitem(.tmpu_, 9L); zeta <- .py_getitem(.tmpu_, 10L); wn <- .py_getitem(.tmpu_, 11L); T <- .py_getitem(.tmpu_, 12L); z <- .py_getitem(.tmpu_, 13L)
})
.py_try({
.tmpu_ <- .py_symbols("a0:4, b0:5, c0:4, d0:4"); a0 <- .py_getitem(.tmpu_, 0L); a1 <- .py_getitem(.tmpu_, 1L); a2 <- .py_getitem(.tmpu_, 2L); a3 <- .py_getitem(.tmpu_, 3L); b0 <- .py_getitem(.tmpu_, 4L); b1 <- .py_getitem(.tmpu_, 5L); b2 <- .py_getitem(.tmpu_, 6L); b3 <- .py_getitem(.tmpu_, 7L); b4 <- .py_getitem(.tmpu_, 8L); c0 <- .py_getitem(.tmpu_, 9L); c1 <- .py_getitem(.tmpu_, 10L); c2 <- .py_getitem(.tmpu_, 11L); c3 <- .py_getitem(.tmpu_, 12L); d0 <- .py_getitem(.tmpu_, 13L); d1 <- .py_getitem(.tmpu_, 14L); d2 <- .py_getitem(.tmpu_, 15L); d3 <- .py_getitem(.tmpu_, 16L)
})
.py_try({
TF1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
TF2 <- TransferFunction(k, 1L, s)
})
.py_try({
TF3 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
TF4 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(tau, s), 1L), s)
})
.py_try({
TF5 <- TransferFunction(k, p, s)
})
.py_try({
TF6 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(k, s), p), s)
})
p_are_floats_equals <- function(n1, n2) {
  "\n    Helper function to compare two floats.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return(.py_method(Float(n1), "is_same", Float(n2), .py_isclose))
}
p_are_complex_equals <- function(n1, n2) {
  "\n    Helper function to compare two complex numbers.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return((.py_truthy(p_are_floats_equals(re(n1), re(n2))) && .py_truthy(p_are_floats_equals(im(n1), im(n2)))))
}
.py_try({
tf <- TransferFunction(1L, .py_add(.py_mul(a, s), b), s)
})
.py_try({
.tmpu_ <- forward_diff(tf, T); numZ <- .py_getitem(.tmpu_, 0L); denZ <- .py_getitem(.tmpu_, 1L)
})
.py_try({
dtf_test_forward1 <- .py_method(.py_method(DiscreteTransferFunction, "from_coeff_lists", numZ, denZ, z, T), "expand")
})
.py_try({
dtf_test_forward2 <- .py_method(.py_method(DiscreteTransferFunction, "from_forward_diff", tf, T, z), "expand")
})
.py_try({
num_manual <- list(.py_div(T, a))
})
.py_try({
den_manual <- list(1L, .py_div(.py_add(.py_neg(a), .py_mul(b, T)), a))
})
.py_try({
dtf_manual <- .py_method(.py_method(DiscreteTransferFunction, "from_coeff_lists", num_manual, den_manual, z, T), "expand")
})
.py_try({
.py_rec("lti:test_TransferFunction_forward_diff:3139", function() (.py_truthy(.py_eq(dtf_test_forward1, dtf_test_forward2)) && .py_truthy(.py_eq(dtf_test_forward2, dtf_manual))))
}))---")
  expect_identical(got[["lti:test_TransferFunction_forward_diff:3139"]], "True")
})

test_that("lti test_TransferFunction_backward_diff matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("a, x, b, c, s, g, d,    p, k, tau, zeta, wn, T, z"); a <- .py_getitem(.tmpu_, 0L); x <- .py_getitem(.tmpu_, 1L); b <- .py_getitem(.tmpu_, 2L); c <- .py_getitem(.tmpu_, 3L); s <- .py_getitem(.tmpu_, 4L); g <- .py_getitem(.tmpu_, 5L); d <- .py_getitem(.tmpu_, 6L); p <- .py_getitem(.tmpu_, 7L); k <- .py_getitem(.tmpu_, 8L); tau <- .py_getitem(.tmpu_, 9L); zeta <- .py_getitem(.tmpu_, 10L); wn <- .py_getitem(.tmpu_, 11L); T <- .py_getitem(.tmpu_, 12L); z <- .py_getitem(.tmpu_, 13L)
})
.py_try({
.tmpu_ <- .py_symbols("a0:4, b0:5, c0:4, d0:4"); a0 <- .py_getitem(.tmpu_, 0L); a1 <- .py_getitem(.tmpu_, 1L); a2 <- .py_getitem(.tmpu_, 2L); a3 <- .py_getitem(.tmpu_, 3L); b0 <- .py_getitem(.tmpu_, 4L); b1 <- .py_getitem(.tmpu_, 5L); b2 <- .py_getitem(.tmpu_, 6L); b3 <- .py_getitem(.tmpu_, 7L); b4 <- .py_getitem(.tmpu_, 8L); c0 <- .py_getitem(.tmpu_, 9L); c1 <- .py_getitem(.tmpu_, 10L); c2 <- .py_getitem(.tmpu_, 11L); c3 <- .py_getitem(.tmpu_, 12L); d0 <- .py_getitem(.tmpu_, 13L); d1 <- .py_getitem(.tmpu_, 14L); d2 <- .py_getitem(.tmpu_, 15L); d3 <- .py_getitem(.tmpu_, 16L)
})
.py_try({
TF1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
TF2 <- TransferFunction(k, 1L, s)
})
.py_try({
TF3 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
TF4 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(tau, s), 1L), s)
})
.py_try({
TF5 <- TransferFunction(k, p, s)
})
.py_try({
TF6 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(k, s), p), s)
})
p_are_floats_equals <- function(n1, n2) {
  "\n    Helper function to compare two floats.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return(.py_method(Float(n1), "is_same", Float(n2), .py_isclose))
}
p_are_complex_equals <- function(n1, n2) {
  "\n    Helper function to compare two complex numbers.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return((.py_truthy(p_are_floats_equals(re(n1), re(n2))) && .py_truthy(p_are_floats_equals(im(n1), im(n2)))))
}
.py_try({
tf <- TransferFunction(1L, .py_add(.py_mul(a, s), b), s)
})
.py_try({
.tmpu_ <- backward_diff(tf, T); numZ <- .py_getitem(.tmpu_, 0L); denZ <- .py_getitem(.tmpu_, 1L)
})
.py_try({
dtf_test_backward1 <- .py_method(.py_method(DiscreteTransferFunction, "from_coeff_lists", numZ, denZ, z, T), "expand")
})
.py_try({
dtf_test_backward2 <- .py_method(.py_method(DiscreteTransferFunction, "from_backward_diff", tf, T, z), "expand")
})
.py_try({
num_manual <- list(.py_div(T, .py_add(a, .py_mul(b, T))), 0L)
})
.py_try({
den_manual <- list(1L, .py_div(.py_neg(a), .py_add(a, .py_mul(b, T))))
})
.py_try({
dtf_manual <- .py_method(.py_method(DiscreteTransferFunction, "from_coeff_lists", num_manual, den_manual, z, T), "expand")
})
.py_try({
.py_rec("lti:test_TransferFunction_backward_diff:3157", function() (.py_truthy(.py_eq(dtf_test_backward1, dtf_test_backward2)) && .py_truthy(.py_eq(dtf_test_backward2, dtf_manual))))
}))---")
  expect_identical(got[["lti:test_TransferFunction_backward_diff:3157"]], "True")
})

test_that("lti test_TransferFunction_phase_margin matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("a, x, b, c, s, g, d,    p, k, tau, zeta, wn, T, z"); a <- .py_getitem(.tmpu_, 0L); x <- .py_getitem(.tmpu_, 1L); b <- .py_getitem(.tmpu_, 2L); c <- .py_getitem(.tmpu_, 3L); s <- .py_getitem(.tmpu_, 4L); g <- .py_getitem(.tmpu_, 5L); d <- .py_getitem(.tmpu_, 6L); p <- .py_getitem(.tmpu_, 7L); k <- .py_getitem(.tmpu_, 8L); tau <- .py_getitem(.tmpu_, 9L); zeta <- .py_getitem(.tmpu_, 10L); wn <- .py_getitem(.tmpu_, 11L); T <- .py_getitem(.tmpu_, 12L); z <- .py_getitem(.tmpu_, 13L)
})
.py_try({
.tmpu_ <- .py_symbols("a0:4, b0:5, c0:4, d0:4"); a0 <- .py_getitem(.tmpu_, 0L); a1 <- .py_getitem(.tmpu_, 1L); a2 <- .py_getitem(.tmpu_, 2L); a3 <- .py_getitem(.tmpu_, 3L); b0 <- .py_getitem(.tmpu_, 4L); b1 <- .py_getitem(.tmpu_, 5L); b2 <- .py_getitem(.tmpu_, 6L); b3 <- .py_getitem(.tmpu_, 7L); b4 <- .py_getitem(.tmpu_, 8L); c0 <- .py_getitem(.tmpu_, 9L); c1 <- .py_getitem(.tmpu_, 10L); c2 <- .py_getitem(.tmpu_, 11L); c3 <- .py_getitem(.tmpu_, 12L); d0 <- .py_getitem(.tmpu_, 13L); d1 <- .py_getitem(.tmpu_, 14L); d2 <- .py_getitem(.tmpu_, 15L); d3 <- .py_getitem(.tmpu_, 16L)
})
.py_try({
TF1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
TF2 <- TransferFunction(k, 1L, s)
})
.py_try({
TF3 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
TF4 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(tau, s), 1L), s)
})
.py_try({
TF5 <- TransferFunction(k, p, s)
})
.py_try({
TF6 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(k, s), p), s)
})
p_are_floats_equals <- function(n1, n2) {
  "\n    Helper function to compare two floats.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return(.py_method(Float(n1), "is_same", Float(n2), .py_isclose))
}
p_are_complex_equals <- function(n1, n2) {
  "\n    Helper function to compare two complex numbers.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return((.py_truthy(p_are_floats_equals(re(n1), re(n2))) && .py_truthy(p_are_floats_equals(im(n1), im(n2)))))
}
.py_try({
tf1 <- TransferFunction(10L, .py_add(.py_pow(p, 3L), 1L), p)
})
.py_try({
tf2 <- TransferFunction(.py_pow(s, 2L), 10L, s)
})
.py_try({
tf3 <- TransferFunction(1L, .py_add(.py_mul(a, s), b), s)
})
.py_try({
tf4 <- TransferFunction(.py_mul(.py_add(s, 1L), .py_exp(.py_div(s, tau))), .py_add(.py_pow(s, 2L), 2L), s)
})
.py_try({
tf_m <- TransferFunctionMatrix(list(list(tf2), list(tf3)))
})
.py_try({
.py_rec("lti:test_TransferFunction_phase_margin:3168", function() phase_margin(tf2))
})
.py_try({
.py_rec("lti:test_TransferFunction_phase_margin:3170", function() phase_margin(tf4))
})
.py_try({
.py_rec("lti:test_TransferFunction_phase_margin:3171", function() phase_margin(tf3))
})
.py_try({
.py_rec("lti:test_TransferFunction_phase_margin:3172", function() phase_margin(MIMOSeries(tf_m)))
}))---")
  expect_identical(got[["lti:test_TransferFunction_phase_margin:3168"]], "0")
  expect_identical(got[["lti:test_TransferFunction_phase_margin:3170"]], "ERR NotImplementedError: Margins for systems with Time delay terms are not supported.")
  expect_identical(got[["lti:test_TransferFunction_phase_margin:3171"]], "ERR ValueError: Extra degree of freedom found. Make sure that there are no free symbols in the dynamical system other than the variable of Laplace transform.")
  expect_identical(got[["lti:test_TransferFunction_phase_margin:3172"]], "ERR ValueError: Margins are only applicable for SISO LTI systems.")
})

test_that("lti test_TransferFunction_gain_margin matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("a, x, b, c, s, g, d,    p, k, tau, zeta, wn, T, z"); a <- .py_getitem(.tmpu_, 0L); x <- .py_getitem(.tmpu_, 1L); b <- .py_getitem(.tmpu_, 2L); c <- .py_getitem(.tmpu_, 3L); s <- .py_getitem(.tmpu_, 4L); g <- .py_getitem(.tmpu_, 5L); d <- .py_getitem(.tmpu_, 6L); p <- .py_getitem(.tmpu_, 7L); k <- .py_getitem(.tmpu_, 8L); tau <- .py_getitem(.tmpu_, 9L); zeta <- .py_getitem(.tmpu_, 10L); wn <- .py_getitem(.tmpu_, 11L); T <- .py_getitem(.tmpu_, 12L); z <- .py_getitem(.tmpu_, 13L)
})
.py_try({
.tmpu_ <- .py_symbols("a0:4, b0:5, c0:4, d0:4"); a0 <- .py_getitem(.tmpu_, 0L); a1 <- .py_getitem(.tmpu_, 1L); a2 <- .py_getitem(.tmpu_, 2L); a3 <- .py_getitem(.tmpu_, 3L); b0 <- .py_getitem(.tmpu_, 4L); b1 <- .py_getitem(.tmpu_, 5L); b2 <- .py_getitem(.tmpu_, 6L); b3 <- .py_getitem(.tmpu_, 7L); b4 <- .py_getitem(.tmpu_, 8L); c0 <- .py_getitem(.tmpu_, 9L); c1 <- .py_getitem(.tmpu_, 10L); c2 <- .py_getitem(.tmpu_, 11L); c3 <- .py_getitem(.tmpu_, 12L); d0 <- .py_getitem(.tmpu_, 13L); d1 <- .py_getitem(.tmpu_, 14L); d2 <- .py_getitem(.tmpu_, 15L); d3 <- .py_getitem(.tmpu_, 16L)
})
.py_try({
TF1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
TF2 <- TransferFunction(k, 1L, s)
})
.py_try({
TF3 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
TF4 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(tau, s), 1L), s)
})
.py_try({
TF5 <- TransferFunction(k, p, s)
})
.py_try({
TF6 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(k, s), p), s)
})
p_are_floats_equals <- function(n1, n2) {
  "\n    Helper function to compare two floats.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return(.py_method(Float(n1), "is_same", Float(n2), .py_isclose))
}
p_are_complex_equals <- function(n1, n2) {
  "\n    Helper function to compare two complex numbers.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return((.py_truthy(p_are_floats_equals(re(n1), re(n2))) && .py_truthy(p_are_floats_equals(im(n1), im(n2)))))
}
.py_try({
tf1 <- TransferFunction(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(5L, .py_add(s, 1L)), .py_sub(s, 5L)), .py_sub(s, 10L)), s)
})
.py_try({
tf2 <- TransferFunction(.py_add(.py_add(.py_pow(s, 2L), .py_mul(2L, s)), 1L), 1L, s)
})
.py_try({
tf3 <- TransferFunction(1L, .py_add(.py_mul(a, s), b), s)
})
.py_try({
tf4 <- TransferFunction(.py_mul(.py_add(s, 1L), .py_exp(.py_div(s, tau))), .py_add(.py_pow(s, 2L), 2L), s)
})
.py_try({
tf_m <- TransferFunctionMatrix(list(list(tf2), list(tf3)))
})
.py_try({
.py_rec("lti:test_TransferFunction_gain_margin:3182", function() gain_margin(tf1))
})
.py_try({
.py_rec("lti:test_TransferFunction_gain_margin:3183", function() gain_margin(tf2))
})
.py_try({
.py_rec("lti:test_TransferFunction_gain_margin:3185", function() gain_margin(tf4))
})
.py_try({
.py_rec("lti:test_TransferFunction_gain_margin:3186", function() gain_margin(tf3))
})
.py_try({
.py_rec("lti:test_TransferFunction_gain_margin:3187", function() gain_margin(MIMOSeries(tf_m)))
}))---")
  expect_identical(got[["lti:test_TransferFunction_gain_margin:3182"]], "-20*log(7/540)/log(10)")
  expect_identical(got[["lti:test_TransferFunction_gain_margin:3183"]], "oo")
  expect_identical(got[["lti:test_TransferFunction_gain_margin:3185"]], "ERR NotImplementedError: Margins for systems with Time delay terms are not supported.")
  expect_identical(got[["lti:test_TransferFunction_gain_margin:3186"]], "ERR ValueError: Extra degree of freedom found. Make sure that there are no free symbols in the dynamical system other than the variable of Laplace transform.")
  expect_identical(got[["lti:test_TransferFunction_gain_margin:3187"]], "ERR ValueError: Margins are only applicable for SISO LTI systems.")
})

test_that("lti test_create_state_space matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("a, x, b, c, s, g, d,    p, k, tau, zeta, wn, T, z"); a <- .py_getitem(.tmpu_, 0L); x <- .py_getitem(.tmpu_, 1L); b <- .py_getitem(.tmpu_, 2L); c <- .py_getitem(.tmpu_, 3L); s <- .py_getitem(.tmpu_, 4L); g <- .py_getitem(.tmpu_, 5L); d <- .py_getitem(.tmpu_, 6L); p <- .py_getitem(.tmpu_, 7L); k <- .py_getitem(.tmpu_, 8L); tau <- .py_getitem(.tmpu_, 9L); zeta <- .py_getitem(.tmpu_, 10L); wn <- .py_getitem(.tmpu_, 11L); T <- .py_getitem(.tmpu_, 12L); z <- .py_getitem(.tmpu_, 13L)
})
.py_try({
.tmpu_ <- .py_symbols("a0:4, b0:5, c0:4, d0:4"); a0 <- .py_getitem(.tmpu_, 0L); a1 <- .py_getitem(.tmpu_, 1L); a2 <- .py_getitem(.tmpu_, 2L); a3 <- .py_getitem(.tmpu_, 3L); b0 <- .py_getitem(.tmpu_, 4L); b1 <- .py_getitem(.tmpu_, 5L); b2 <- .py_getitem(.tmpu_, 6L); b3 <- .py_getitem(.tmpu_, 7L); b4 <- .py_getitem(.tmpu_, 8L); c0 <- .py_getitem(.tmpu_, 9L); c1 <- .py_getitem(.tmpu_, 10L); c2 <- .py_getitem(.tmpu_, 11L); c3 <- .py_getitem(.tmpu_, 12L); d0 <- .py_getitem(.tmpu_, 13L); d1 <- .py_getitem(.tmpu_, 14L); d2 <- .py_getitem(.tmpu_, 15L); d3 <- .py_getitem(.tmpu_, 16L)
})
.py_try({
TF1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
TF2 <- TransferFunction(k, 1L, s)
})
.py_try({
TF3 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
TF4 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(tau, s), 1L), s)
})
.py_try({
TF5 <- TransferFunction(k, p, s)
})
.py_try({
TF6 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(k, s), p), s)
})
p_are_floats_equals <- function(n1, n2) {
  "\n    Helper function to compare two floats.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return(.py_method(Float(n1), "is_same", Float(n2), .py_isclose))
}
p_are_complex_equals <- function(n1, n2) {
  "\n    Helper function to compare two complex numbers.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return((.py_truthy(p_are_floats_equals(re(n1), re(n2))) && .py_truthy(p_are_floats_equals(im(n1), im(n2)))))
}
.py_try({
A <- .py_Matrix(list(list(0L, 1L), list(1L, 0L)))
})
.py_try({
B <- .py_Matrix(list(1L, 0L))
})
.py_try({
C <- .py_Matrix(list(list(0L, 1L)))
})
.py_try({
D <- .py_Matrix(list(0L))
})
.py_try({
cont_ss1 <- create_state_space(A, B, C, D)
})
.py_try({
.py_rec("lti:test_create_state_space:3196", function() .py_isinstance(cont_ss1, StateSpace))
})
.py_try({
cont_ss2 <- create_state_space(A, B, C, D, sampling_time = 0L)
})
.py_try({
.py_rec("lti:test_create_state_space:3199", function() .py_isinstance(cont_ss2, StateSpace))
})
.py_try({
disc_ss1 <- create_state_space(A, B, C, D, .py_flt(0.1))
})
.py_try({
.py_rec("lti:test_create_state_space:3202", function() .py_isinstance(disc_ss1, DiscreteStateSpace))
})
.py_try({
disc_ss2 <- create_state_space(A, B, C, D, T)
})
.py_try({
.py_rec("lti:test_create_state_space:3205", function() .py_isinstance(disc_ss2, DiscreteStateSpace))
}))---")
  expect_identical(got[["lti:test_create_state_space:3196"]], "True")
  expect_identical(got[["lti:test_create_state_space:3199"]], "True")
  expect_identical(got[["lti:test_create_state_space:3202"]], "True")
  expect_identical(got[["lti:test_create_state_space:3205"]], "True")
})

test_that("lti test_StateSpace_construction matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("a, x, b, c, s, g, d,    p, k, tau, zeta, wn, T, z"); a <- .py_getitem(.tmpu_, 0L); x <- .py_getitem(.tmpu_, 1L); b <- .py_getitem(.tmpu_, 2L); c <- .py_getitem(.tmpu_, 3L); s <- .py_getitem(.tmpu_, 4L); g <- .py_getitem(.tmpu_, 5L); d <- .py_getitem(.tmpu_, 6L); p <- .py_getitem(.tmpu_, 7L); k <- .py_getitem(.tmpu_, 8L); tau <- .py_getitem(.tmpu_, 9L); zeta <- .py_getitem(.tmpu_, 10L); wn <- .py_getitem(.tmpu_, 11L); T <- .py_getitem(.tmpu_, 12L); z <- .py_getitem(.tmpu_, 13L)
})
.py_try({
.tmpu_ <- .py_symbols("a0:4, b0:5, c0:4, d0:4"); a0 <- .py_getitem(.tmpu_, 0L); a1 <- .py_getitem(.tmpu_, 1L); a2 <- .py_getitem(.tmpu_, 2L); a3 <- .py_getitem(.tmpu_, 3L); b0 <- .py_getitem(.tmpu_, 4L); b1 <- .py_getitem(.tmpu_, 5L); b2 <- .py_getitem(.tmpu_, 6L); b3 <- .py_getitem(.tmpu_, 7L); b4 <- .py_getitem(.tmpu_, 8L); c0 <- .py_getitem(.tmpu_, 9L); c1 <- .py_getitem(.tmpu_, 10L); c2 <- .py_getitem(.tmpu_, 11L); c3 <- .py_getitem(.tmpu_, 12L); d0 <- .py_getitem(.tmpu_, 13L); d1 <- .py_getitem(.tmpu_, 14L); d2 <- .py_getitem(.tmpu_, 15L); d3 <- .py_getitem(.tmpu_, 16L)
})
.py_try({
TF1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
TF2 <- TransferFunction(k, 1L, s)
})
.py_try({
TF3 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
TF4 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(tau, s), 1L), s)
})
.py_try({
TF5 <- TransferFunction(k, p, s)
})
.py_try({
TF6 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(k, s), p), s)
})
p_are_floats_equals <- function(n1, n2) {
  "\n    Helper function to compare two floats.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return(.py_method(Float(n1), "is_same", Float(n2), .py_isclose))
}
p_are_complex_equals <- function(n1, n2) {
  "\n    Helper function to compare two complex numbers.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return((.py_truthy(p_are_floats_equals(re(n1), re(n2))) && .py_truthy(p_are_floats_equals(im(n1), im(n2)))))
}
.py_try({
A1 <- .py_Matrix(list(list(0L, 1L), list(1L, 0L)))
})
.py_try({
B1 <- .py_Matrix(list(1L, 0L))
})
.py_try({
C1 <- .py_Matrix(list(list(0L, 1L)))
})
.py_try({
D1 <- .py_Matrix(list(0L))
})
.py_try({
ss1 <- StateSpace(A1, B1, C1, D1)
})
.py_try({
.py_rec("lti:test_StateSpace_construction:3215", function() .py_attr(ss1, "state_matrix"))
})
.py_try({
.py_rec("lti:test_StateSpace_construction:3216", function() .py_attr(ss1, "input_matrix"))
})
.py_try({
.py_rec("lti:test_StateSpace_construction:3217", function() .py_attr(ss1, "output_matrix"))
})
.py_try({
.py_rec("lti:test_StateSpace_construction:3218", function() .py_attr(ss1, "feedforward_matrix"))
})
.py_try({
.py_rec("lti:test_StateSpace_construction:3219", function() .py_attr(ss1, "sampling_time"))
})
.py_try({
.py_rec("lti:test_StateSpace_construction:3220", function() .py_attr(ss1, "args"))
})
.py_try({
ss2 <- StateSpace(.py_Matrix(list(a0)), .py_Matrix(list(a1)), .py_Matrix(list(a2)), .py_Matrix(list(a3)))
})
.py_try({
.py_rec("lti:test_StateSpace_construction:3226", function() .py_attr(ss2, "state_matrix"))
})
.py_try({
.py_rec("lti:test_StateSpace_construction:3227", function() .py_attr(ss2, "input_matrix"))
})
.py_try({
.py_rec("lti:test_StateSpace_construction:3228", function() .py_attr(ss2, "output_matrix"))
})
.py_try({
.py_rec("lti:test_StateSpace_construction:3229", function() .py_attr(ss2, "feedforward_matrix"))
})
.py_try({
.py_rec("lti:test_StateSpace_construction:3230", function() .py_attr(ss2, "args"))
})
.py_try({
ss3 <- StateSpace(.py_Matrix(list(list(.py_flt(-1.5), -2L), list(1L, 0L))), .py_Matrix(list(list(.py_flt(0.5), 0L), list(0L, 1L))), .py_Matrix(list(list(0L, 1L), list(0L, 2L))), .py_Matrix(list(list(2L, 2L), list(1L, 1L))))
})
.py_try({
.py_rec("lti:test_StateSpace_construction:3238", function() .py_attr(ss3, "state_matrix"))
})
.py_try({
.py_rec("lti:test_StateSpace_construction:3239", function() .py_attr(ss3, "input_matrix"))
})
.py_try({
.py_rec("lti:test_StateSpace_construction:3240", function() .py_attr(ss3, "output_matrix"))
})
.py_try({
.py_rec("lti:test_StateSpace_construction:3241", function() .py_attr(ss3, "feedforward_matrix"))
})
.py_try({
.py_rec("lti:test_StateSpace_construction:3242", function() .py_attr(ss3, "args"))
})
.py_try({
A4 <- .py_Matrix(list(list(a0, a1), list(a2, a3)))
})
.py_try({
B4 <- .py_Matrix(list(list(b0, b1), list(b2, b3)))
})
.py_try({
C4 <- .py_Matrix(list(list(c0, c1), list(c2, c3)))
})
.py_try({
D4 <- .py_Matrix(list(list(d0, d1), list(d2, d3)))
})
.py_try({
ss4 <- StateSpace(A4, B4, C4, D4)
})
.py_try({
.py_rec("lti:test_StateSpace_construction:3258", function() .py_attr(ss4, "state_matrix"))
})
.py_try({
.py_rec("lti:test_StateSpace_construction:3259", function() .py_attr(ss4, "input_matrix"))
})
.py_try({
.py_rec("lti:test_StateSpace_construction:3260", function() .py_attr(ss4, "output_matrix"))
})
.py_try({
.py_rec("lti:test_StateSpace_construction:3261", function() .py_attr(ss4, "feedforward_matrix"))
})
.py_try({
.py_rec("lti:test_StateSpace_construction:3262", function() .py_attr(ss4, "args"))
})
.py_try({
ss5 <- StateSpace()
})
.py_try({
.py_rec("lti:test_StateSpace_construction:3273", function() .py_attr(ss5, "args"))
})
.py_try({
A6 <- .py_Matrix(list(list(0L, 1L), list(1L, 0L)))
})
.py_try({
B6 <- .py_Matrix(list(1L, 1L))
})
.py_try({
ss6 <- StateSpace(A6, B6)
})
.py_try({
.py_rec("lti:test_StateSpace_construction:3279", function() .py_attr(ss6, "state_matrix"))
})
.py_try({
.py_rec("lti:test_StateSpace_construction:3280", function() .py_attr(ss6, "input_matrix"))
})
.py_try({
.py_rec("lti:test_StateSpace_construction:3281", function() .py_attr(ss6, "output_matrix"))
})
.py_try({
.py_rec("lti:test_StateSpace_construction:3282", function() .py_attr(ss6, "feedforward_matrix"))
})
.py_try({
.py_rec("lti:test_StateSpace_construction:3283", function() .py_attr(ss6, "args"))
})
.py_try({
.py_rec("lti:test_StateSpace_construction:3293", function() .py_attr(ss1, "is_SISO"))
})
.py_try({
.py_rec("lti:test_StateSpace_construction:3294", function() .py_attr(ss2, "is_SISO"))
})
.py_try({
.py_rec("lti:test_StateSpace_construction:3295", function() .py_attr(ss3, "is_SISO"))
})
.py_try({
.py_rec("lti:test_StateSpace_construction:3296", function() .py_attr(ss4, "is_SISO"))
})
.py_try({
.py_rec("lti:test_StateSpace_construction:3297", function() .py_attr(ss5, "is_SISO"))
})
.py_try({
.py_rec("lti:test_StateSpace_construction:3298", function() .py_attr(ss6, "is_SISO"))
})
.py_try({
.py_rec("lti:test_StateSpace_construction:3301", function() StateSpace(.py_Matrix(list(s, .py_pow(.py_add(s, 1L), 2L))), .py_Matrix(list(.py_add(s, 1L))), .py_Matrix(list(.py_sub(.py_pow(s, 2L), 1L))), .py_Matrix(list(.py_mul(2L, s)))))
})
.py_try({
.py_rec("lti:test_StateSpace_construction:3303", function() StateSpace(.py_Matrix(list(s)), .py_Matrix(list(.py_add(s, 1L), .py_add(.py_pow(s, 3L), 1L))), .py_Matrix(list(.py_sub(.py_pow(s, 2L), 1L))), .py_Matrix(list(.py_mul(2L, s)))))
})
.py_try({
.py_rec("lti:test_StateSpace_construction:3305", function() StateSpace(.py_Matrix(list(s)), .py_Matrix(list(.py_add(s, 1L))), .py_Matrix(list(list(.py_sub(.py_pow(s, 2L), 1L)), list(.py_add(.py_add(.py_pow(s, 2L), .py_mul(2L, s)), 1L)))), .py_Matrix(list(.py_mul(2L, s)))))
})
.py_try({
.py_rec("lti:test_StateSpace_construction:3307", function() StateSpace(.py_Matrix(list(list(.py_neg(s), .py_neg(s)), list(s, 0L))), .py_Matrix(list(list(.py_div(s, 2L), 0L), list(0L, s))), .py_Matrix(list(list(0L, s))), .py_Matrix(list(list(.py_mul(2L, s), .py_mul(2L, s)), list(s, s)))))
})
.py_try({
.py_rec("lti:test_StateSpace_construction:3313", function() StateSpace(.py_pow(s, 2L), .py_add(s, 1L), .py_mul(2L, s), 1L))
})
.py_try({
.py_rec("lti:test_StateSpace_construction:3314", function() StateSpace(.py_Matrix(list(2L, .py_flt(0.5))), .py_Matrix(list(-1L)), .py_Matrix(list(1L)), 0L))
}))---")
  expect_identical(got[["lti:test_StateSpace_construction:3215"]], "Matrix([[0, 1], [1, 0]])")
  expect_identical(got[["lti:test_StateSpace_construction:3216"]], "Matrix([[1], [0]])")
  expect_identical(got[["lti:test_StateSpace_construction:3217"]], "Matrix([[0, 1]])")
  expect_identical(got[["lti:test_StateSpace_construction:3218"]], "Matrix([[0]])")
  expect_identical(got[["lti:test_StateSpace_construction:3219"]], "0")
  expect_identical(got[["lti:test_StateSpace_construction:3220"]], "(Matrix([[0, 1], [1, 0]]), Matrix([[1], [0]]), Matrix([[0, 1]]), Matrix([[0]]))")
  expect_identical(got[["lti:test_StateSpace_construction:3226"]], "Matrix([[a0]])")
  expect_identical(got[["lti:test_StateSpace_construction:3227"]], "Matrix([[a1]])")
  expect_identical(got[["lti:test_StateSpace_construction:3228"]], "Matrix([[a2]])")
  expect_identical(got[["lti:test_StateSpace_construction:3229"]], "Matrix([[a3]])")
  expect_identical(got[["lti:test_StateSpace_construction:3230"]], "(Matrix([[a0]]), Matrix([[a1]]), Matrix([[a2]]), Matrix([[a3]]))")
  expect_identical(got[["lti:test_StateSpace_construction:3238"]], "Matrix([[-1.50000000000000, -2], [1, 0]])")
  expect_identical(got[["lti:test_StateSpace_construction:3239"]], "Matrix([[0.500000000000000, 0], [0, 1]])")
  expect_identical(got[["lti:test_StateSpace_construction:3240"]], "Matrix([[0, 1], [0, 2]])")
  expect_identical(got[["lti:test_StateSpace_construction:3241"]], "Matrix([[2, 2], [1, 1]])")
  expect_identical(got[["lti:test_StateSpace_construction:3242"]], "(Matrix([[-1.50000000000000, -2], [1, 0]]), Matrix([[0.500000000000000, 0], [0, 1]]), Matrix([[0, 1], [0, 2]]), Matrix([[2, 2], [1, 1]]))")
  expect_identical(got[["lti:test_StateSpace_construction:3258"]], "Matrix([[a0, a1], [a2, a3]])")
  expect_identical(got[["lti:test_StateSpace_construction:3259"]], "Matrix([[b0, b1], [b2, b3]])")
  expect_identical(got[["lti:test_StateSpace_construction:3260"]], "Matrix([[c0, c1], [c2, c3]])")
  expect_identical(got[["lti:test_StateSpace_construction:3261"]], "Matrix([[d0, d1], [d2, d3]])")
  expect_identical(got[["lti:test_StateSpace_construction:3262"]], "(Matrix([[a0, a1], [a2, a3]]), Matrix([[b0, b1], [b2, b3]]), Matrix([[c0, c1], [c2, c3]]), Matrix([[d0, d1], [d2, d3]]))")
  expect_identical(got[["lti:test_StateSpace_construction:3273"]], "(Matrix([[0]]), Matrix([[0]]), Matrix([[0]]), Matrix([[0]]))")
  expect_identical(got[["lti:test_StateSpace_construction:3279"]], "Matrix([[0, 1], [1, 0]])")
  expect_identical(got[["lti:test_StateSpace_construction:3280"]], "Matrix([[1], [1]])")
  expect_identical(got[["lti:test_StateSpace_construction:3281"]], "Matrix([[0, 0]])")
  expect_identical(got[["lti:test_StateSpace_construction:3282"]], "Matrix([[0]])")
  expect_identical(got[["lti:test_StateSpace_construction:3283"]], "(Matrix([[0, 1], [1, 0]]), Matrix([[1], [1]]), Matrix([[0, 0]]), Matrix([[0]]))")
  expect_identical(got[["lti:test_StateSpace_construction:3293"]], "True")
  expect_identical(got[["lti:test_StateSpace_construction:3294"]], "True")
  expect_identical(got[["lti:test_StateSpace_construction:3295"]], "False")
  expect_identical(got[["lti:test_StateSpace_construction:3296"]], "False")
  expect_identical(got[["lti:test_StateSpace_construction:3297"]], "True")
  expect_identical(got[["lti:test_StateSpace_construction:3298"]], "True")
  expect_identical(got[["lti:test_StateSpace_construction:3301"]], "ERR ShapeError: Matrix A must be a square matrix.")
  expect_identical(got[["lti:test_StateSpace_construction:3303"]], "ERR ShapeError: Matrices A and B must have the same number of rows.")
  expect_identical(got[["lti:test_StateSpace_construction:3305"]], "ERR ShapeError: Matrices C and D must have the same number of rows.")
  expect_identical(got[["lti:test_StateSpace_construction:3307"]], "ERR ShapeError: Matrices C and D must have the same number of rows.")
  expect_identical(got[["lti:test_StateSpace_construction:3313"]], "ERR TypeError: A, B, C and D inputs must all be sympy Matrices.")
  expect_identical(got[["lti:test_StateSpace_construction:3314"]], "ERR TypeError: A, B, C and D inputs must all be sympy Matrices.")
})

test_that("lti test_StateSpace_add_mul matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("a, x, b, c, s, g, d,    p, k, tau, zeta, wn, T, z"); a <- .py_getitem(.tmpu_, 0L); x <- .py_getitem(.tmpu_, 1L); b <- .py_getitem(.tmpu_, 2L); c <- .py_getitem(.tmpu_, 3L); s <- .py_getitem(.tmpu_, 4L); g <- .py_getitem(.tmpu_, 5L); d <- .py_getitem(.tmpu_, 6L); p <- .py_getitem(.tmpu_, 7L); k <- .py_getitem(.tmpu_, 8L); tau <- .py_getitem(.tmpu_, 9L); zeta <- .py_getitem(.tmpu_, 10L); wn <- .py_getitem(.tmpu_, 11L); T <- .py_getitem(.tmpu_, 12L); z <- .py_getitem(.tmpu_, 13L)
})
.py_try({
.tmpu_ <- .py_symbols("a0:4, b0:5, c0:4, d0:4"); a0 <- .py_getitem(.tmpu_, 0L); a1 <- .py_getitem(.tmpu_, 1L); a2 <- .py_getitem(.tmpu_, 2L); a3 <- .py_getitem(.tmpu_, 3L); b0 <- .py_getitem(.tmpu_, 4L); b1 <- .py_getitem(.tmpu_, 5L); b2 <- .py_getitem(.tmpu_, 6L); b3 <- .py_getitem(.tmpu_, 7L); b4 <- .py_getitem(.tmpu_, 8L); c0 <- .py_getitem(.tmpu_, 9L); c1 <- .py_getitem(.tmpu_, 10L); c2 <- .py_getitem(.tmpu_, 11L); c3 <- .py_getitem(.tmpu_, 12L); d0 <- .py_getitem(.tmpu_, 13L); d1 <- .py_getitem(.tmpu_, 14L); d2 <- .py_getitem(.tmpu_, 15L); d3 <- .py_getitem(.tmpu_, 16L)
})
.py_try({
TF1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
TF2 <- TransferFunction(k, 1L, s)
})
.py_try({
TF3 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
TF4 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(tau, s), 1L), s)
})
.py_try({
TF5 <- TransferFunction(k, p, s)
})
.py_try({
TF6 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(k, s), p), s)
})
p_are_floats_equals <- function(n1, n2) {
  "\n    Helper function to compare two floats.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return(.py_method(Float(n1), "is_same", Float(n2), .py_isclose))
}
p_are_complex_equals <- function(n1, n2) {
  "\n    Helper function to compare two complex numbers.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return((.py_truthy(p_are_floats_equals(re(n1), re(n2))) && .py_truthy(p_are_floats_equals(im(n1), im(n2)))))
}
.py_try({
A1 <- .py_Matrix(list(list(4L, 1L), list(2L, -3L)))
})
.py_try({
B1 <- .py_Matrix(list(list(5L, 2L), list(-3L, -3L)))
})
.py_try({
C1 <- .py_Matrix(list(list(2L, -4L), list(0L, 1L)))
})
.py_try({
D1 <- .py_Matrix(list(list(3L, 2L), list(1L, -1L)))
})
.py_try({
ss1 <- StateSpace(A1, B1, C1, D1)
})
.py_try({
A2 <- .py_Matrix(list(list(-3L, 4L, 2L), list(-1L, -3L, 0L), list(2L, 5L, 3L)))
})
.py_try({
B2 <- .py_Matrix(list(list(1L, 4L), list(-3L, -3L), list(-2L, 1L)))
})
.py_try({
C2 <- .py_Matrix(list(list(4L, 2L, -3L), list(1L, 4L, 3L)))
})
.py_try({
D2 <- .py_Matrix(list(list(-2L, 4L), list(0L, 1L)))
})
.py_try({
ss2 <- StateSpace(A2, B2, C2, D2)
})
.py_try({
ss3 <- StateSpace()
})
.py_try({
ss4 <- StateSpace(.py_Matrix(list(1L)), .py_Matrix(list(2L)), .py_Matrix(list(3L)), .py_Matrix(list(4L)))
})
.py_try({
expected_add <- StateSpace(.py_Matrix(list(list(4L, 1L, 0L, 0L, 0L), list(2L, -3L, 0L, 0L, 0L), list(0L, 0L, -3L, 4L, 2L), list(0L, 0L, -1L, -3L, 0L), list(0L, 0L, 2L, 5L, 3L))), .py_Matrix(list(list(5L, 2L), list(-3L, -3L), list(1L, 4L), list(-3L, -3L), list(-2L, 1L))), .py_Matrix(list(list(2L, -4L, 4L, 2L, -3L), list(0L, 1L, 1L, 4L, 3L))), .py_Matrix(list(list(1L, 6L), list(1L, 0L))))
})
.py_try({
expected_mul <- StateSpace(.py_Matrix(list(list(-3L, 4L, 2L, 0L, 0L), list(-1L, -3L, 0L, 0L, 0L), list(2L, 5L, 3L, 0L, 0L), list(22L, 18L, -9L, 4L, 1L), list(-15L, -18L, 0L, 2L, -3L))), .py_Matrix(list(list(1L, 4L), list(-3L, -3L), list(-2L, 1L), list(-10L, 22L), list(6L, -15L))), .py_Matrix(list(list(14L, 14L, -3L, 2L, -4L), list(3L, -2L, -6L, 0L, 1L))), .py_Matrix(list(list(-6L, 14L), list(-2L, 3L))))
})
.py_try({
.py_rec("lti:test_StateSpace_add_mul:3374", function() .py_add(ss1, ss2))
})
.py_try({
.py_rec("lti:test_StateSpace_add_mul:3375", function() .py_mul(ss1, ss2))
})
.py_try({
.py_rec("lti:test_StateSpace_add_mul:3376", function() .py_add(ss3, .py_div(1L, 2L)))
})
.py_try({
.py_rec("lti:test_StateSpace_add_mul:3377", function() .py_mul(ss4, .py_flt(1.5)))
})
.py_try({
.py_rec("lti:test_StateSpace_add_mul:3378", function() .py_mul(.py_flt(1.5), ss4))
})
.py_try({
.py_rec("lti:test_StateSpace_add_mul:3379", function() .py_add(ss1, ss3))
})
.py_try({
.py_rec("lti:test_StateSpace_add_mul:3380", function() .py_mul(ss2, ss4))
}))---")
  expect_identical(got[["lti:test_StateSpace_add_mul:3374"]], "StateSpace(Matrix([\\n[4,  1,  0,  0, 0],\\n[2, -3,  0,  0, 0],\\n[0,  0, -3,  4, 2],\\n[0,  0, -1, -3, 0],\\n[0,  0,  2,  5, 3]]), Matrix([\\n[ 5,  2],\\n[-3, -3],\\n[ 1,  4],\\n[-3, -3],\\n[-2,  1]]), Matrix([\\n[2, -4, 4, 2, -3],\\n[0,  1, 1, 4,  3]]), Matrix([\\n[1, 6],\\n[1, 0]]))")
  expect_identical(got[["lti:test_StateSpace_add_mul:3375"]], "StateSpace(Matrix([\\n[ -3,   4,  2, 0,  0],\\n[ -1,  -3,  0, 0,  0],\\n[  2,   5,  3, 0,  0],\\n[ 22,  18, -9, 4,  1],\\n[-15, -18,  0, 2, -3]]), Matrix([\\n[  1,   4],\\n[ -3,  -3],\\n[ -2,   1],\\n[-10,  22],\\n[  6, -15]]), Matrix([\\n[14, 14, -3, 2, -4],\\n[ 3, -2, -6, 0,  1]]), Matrix([\\n[-6, 14],\\n[-2,  3]]))")
  expect_identical(got[["lti:test_StateSpace_add_mul:3376"]], "StateSpace(Matrix([[0]]), Matrix([[0]]), Matrix([[0]]), Matrix([[0.5]]))")
  expect_identical(got[["lti:test_StateSpace_add_mul:3377"]], "StateSpace(Matrix([[1]]), Matrix([[2]]), Matrix([[4.5]]), Matrix([[6.0]]))")
  expect_identical(got[["lti:test_StateSpace_add_mul:3378"]], "StateSpace(Matrix([[1]]), Matrix([[3.0]]), Matrix([[3]]), Matrix([[6.0]]))")
  expect_identical(got[["lti:test_StateSpace_add_mul:3379"]], "ERR ShapeError: Systems with incompatible inputs and outputs cannot be added.")
  expect_identical(got[["lti:test_StateSpace_add_mul:3380"]], "ERR ShapeError: Systems with incompatible inputs and outputs cannot be multiplied.")
})

test_that("lti test_StateSpace_negation matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("a, x, b, c, s, g, d,    p, k, tau, zeta, wn, T, z"); a <- .py_getitem(.tmpu_, 0L); x <- .py_getitem(.tmpu_, 1L); b <- .py_getitem(.tmpu_, 2L); c <- .py_getitem(.tmpu_, 3L); s <- .py_getitem(.tmpu_, 4L); g <- .py_getitem(.tmpu_, 5L); d <- .py_getitem(.tmpu_, 6L); p <- .py_getitem(.tmpu_, 7L); k <- .py_getitem(.tmpu_, 8L); tau <- .py_getitem(.tmpu_, 9L); zeta <- .py_getitem(.tmpu_, 10L); wn <- .py_getitem(.tmpu_, 11L); T <- .py_getitem(.tmpu_, 12L); z <- .py_getitem(.tmpu_, 13L)
})
.py_try({
.tmpu_ <- .py_symbols("a0:4, b0:5, c0:4, d0:4"); a0 <- .py_getitem(.tmpu_, 0L); a1 <- .py_getitem(.tmpu_, 1L); a2 <- .py_getitem(.tmpu_, 2L); a3 <- .py_getitem(.tmpu_, 3L); b0 <- .py_getitem(.tmpu_, 4L); b1 <- .py_getitem(.tmpu_, 5L); b2 <- .py_getitem(.tmpu_, 6L); b3 <- .py_getitem(.tmpu_, 7L); b4 <- .py_getitem(.tmpu_, 8L); c0 <- .py_getitem(.tmpu_, 9L); c1 <- .py_getitem(.tmpu_, 10L); c2 <- .py_getitem(.tmpu_, 11L); c3 <- .py_getitem(.tmpu_, 12L); d0 <- .py_getitem(.tmpu_, 13L); d1 <- .py_getitem(.tmpu_, 14L); d2 <- .py_getitem(.tmpu_, 15L); d3 <- .py_getitem(.tmpu_, 16L)
})
.py_try({
TF1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
TF2 <- TransferFunction(k, 1L, s)
})
.py_try({
TF3 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
TF4 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(tau, s), 1L), s)
})
.py_try({
TF5 <- TransferFunction(k, p, s)
})
.py_try({
TF6 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(k, s), p), s)
})
p_are_floats_equals <- function(n1, n2) {
  "\n    Helper function to compare two floats.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return(.py_method(Float(n1), "is_same", Float(n2), .py_isclose))
}
p_are_complex_equals <- function(n1, n2) {
  "\n    Helper function to compare two complex numbers.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return((.py_truthy(p_are_floats_equals(re(n1), re(n2))) && .py_truthy(p_are_floats_equals(im(n1), im(n2)))))
}
.py_try({
A <- .py_Matrix(list(list(a0, a1), list(a2, a3)))
})
.py_try({
B <- .py_Matrix(list(list(b0, b1), list(b2, b3)))
})
.py_try({
C <- .py_Matrix(list(list(c0, c1), list(c1, c2), list(c2, c3)))
})
.py_try({
D <- .py_Matrix(list(list(d0, d1), list(d1, d2), list(d2, d3)))
})
.py_try({
SS <- StateSpace(A, B, C, D)
})
.py_try({
SS_neg <- .py_neg(SS)
})
.py_try({
state_mat <- .py_Matrix(list(list(-1L, 1L), list(1L, -1L)))
})
.py_try({
input_mat <- .py_Matrix(list(1L, -1L))
})
.py_try({
output_mat <- .py_Matrix(list(list(-1L, 1L)))
})
.py_try({
feedforward_mat <- .py_Matrix(list(1L))
})
.py_try({
system <- StateSpace(state_mat, input_mat, output_mat, feedforward_mat)
})
.py_try({
.py_rec("lti:test_StateSpace_negation:3396", function() SS_neg)
})
.py_try({
.py_rec("lti:test_StateSpace_negation:3407", function() .py_neg(system))
})
.py_try({
.py_rec("lti:test_StateSpace_negation:3413", function() .py_neg(SS_neg))
})
.py_try({
.py_rec("lti:test_StateSpace_negation:3414", function() .py_neg(.py_neg(.py_neg(.py_neg(system)))))
}))---")
  expect_identical(got[["lti:test_StateSpace_negation:3396"]], "StateSpace(Matrix([\\n[a0, a1],\\n[a2, a3]]), Matrix([\\n[b0, b1],\\n[b2, b3]]), Matrix([\\n[-c0, -c1],\\n[-c1, -c2],\\n[-c2, -c3]]), Matrix([\\n[-d0, -d1],\\n[-d1, -d2],\\n[-d2, -d3]]))")
  expect_identical(got[["lti:test_StateSpace_negation:3407"]], "StateSpace(Matrix([\\n[-1,  1],\\n[ 1, -1]]), Matrix([\\n[ 1],\\n[-1]]), Matrix([[1, -1]]), Matrix([[-1]]))")
  expect_identical(got[["lti:test_StateSpace_negation:3413"]], "StateSpace(Matrix([\\n[a0, a1],\\n[a2, a3]]), Matrix([\\n[b0, b1],\\n[b2, b3]]), Matrix([\\n[c0, c1],\\n[c1, c2],\\n[c2, c3]]), Matrix([\\n[d0, d1],\\n[d1, d2],\\n[d2, d3]]))")
  expect_identical(got[["lti:test_StateSpace_negation:3414"]], "StateSpace(Matrix([\\n[-1,  1],\\n[ 1, -1]]), Matrix([\\n[ 1],\\n[-1]]), Matrix([[-1, 1]]), Matrix([[1]]))")
})

test_that("lti test_SymPy_substitution_functions matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("a, x, b, c, s, g, d,    p, k, tau, zeta, wn, T, z"); a <- .py_getitem(.tmpu_, 0L); x <- .py_getitem(.tmpu_, 1L); b <- .py_getitem(.tmpu_, 2L); c <- .py_getitem(.tmpu_, 3L); s <- .py_getitem(.tmpu_, 4L); g <- .py_getitem(.tmpu_, 5L); d <- .py_getitem(.tmpu_, 6L); p <- .py_getitem(.tmpu_, 7L); k <- .py_getitem(.tmpu_, 8L); tau <- .py_getitem(.tmpu_, 9L); zeta <- .py_getitem(.tmpu_, 10L); wn <- .py_getitem(.tmpu_, 11L); T <- .py_getitem(.tmpu_, 12L); z <- .py_getitem(.tmpu_, 13L)
})
.py_try({
.tmpu_ <- .py_symbols("a0:4, b0:5, c0:4, d0:4"); a0 <- .py_getitem(.tmpu_, 0L); a1 <- .py_getitem(.tmpu_, 1L); a2 <- .py_getitem(.tmpu_, 2L); a3 <- .py_getitem(.tmpu_, 3L); b0 <- .py_getitem(.tmpu_, 4L); b1 <- .py_getitem(.tmpu_, 5L); b2 <- .py_getitem(.tmpu_, 6L); b3 <- .py_getitem(.tmpu_, 7L); b4 <- .py_getitem(.tmpu_, 8L); c0 <- .py_getitem(.tmpu_, 9L); c1 <- .py_getitem(.tmpu_, 10L); c2 <- .py_getitem(.tmpu_, 11L); c3 <- .py_getitem(.tmpu_, 12L); d0 <- .py_getitem(.tmpu_, 13L); d1 <- .py_getitem(.tmpu_, 14L); d2 <- .py_getitem(.tmpu_, 15L); d3 <- .py_getitem(.tmpu_, 16L)
})
.py_try({
TF1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
TF2 <- TransferFunction(k, 1L, s)
})
.py_try({
TF3 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
TF4 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(tau, s), 1L), s)
})
.py_try({
TF5 <- TransferFunction(k, p, s)
})
.py_try({
TF6 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(k, s), p), s)
})
p_are_floats_equals <- function(n1, n2) {
  "\n    Helper function to compare two floats.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return(.py_method(Float(n1), "is_same", Float(n2), .py_isclose))
}
p_are_complex_equals <- function(n1, n2) {
  "\n    Helper function to compare two complex numbers.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return((.py_truthy(p_are_floats_equals(re(n1), re(n2))) && .py_truthy(p_are_floats_equals(im(n1), im(n2)))))
}
.py_try({
ss1 <- StateSpace(.py_Matrix(list(s)), .py_Matrix(list(.py_pow(.py_add(s, 1L), 2L))), .py_Matrix(list(.py_sub(.py_pow(s, 2L), 1L))), .py_Matrix(list(.py_mul(2L, s))))
})
.py_try({
ss2 <- StateSpace(.py_Matrix(list(.py_add(s, p))), .py_Matrix(list(.py_mul(.py_add(s, 1L), .py_sub(p, 1L)))), .py_Matrix(list(.py_sub(.py_pow(p, 3L), .py_pow(s, 3L)))), .py_Matrix(list(.py_sub(s, p))))
})
.py_try({
dss1 <- DiscreteStateSpace(.py_Matrix(list(s)), .py_Matrix(list(.py_pow(.py_add(s, 1L), 2L))), .py_Matrix(list(.py_sub(.py_pow(s, 2L), 1L))), .py_Matrix(list(.py_mul(2L, s))), .py_flt(0.1))
})
.py_try({
dss2 <- DiscreteStateSpace(.py_Matrix(list(.py_add(s, p))), .py_Matrix(list(.py_mul(.py_add(s, 1L), .py_sub(p, 1L)))), .py_Matrix(list(.py_sub(.py_pow(p, 3L), .py_pow(s, 3L)))), .py_Matrix(list(.py_sub(s, p))), .py_flt(0.1))
})
.py_try({
.py_rec("lti:test_SymPy_substitution_functions:3426", function() .py_method(ss1, "subs", .py_dict(list(s), list(5L))))
})
.py_try({
.py_rec("lti:test_SymPy_substitution_functions:3427", function() .py_method(ss2, "subs", .py_dict(list(p), list(1L))))
})
.py_try({
.py_rec("lti:test_SymPy_substitution_functions:3429", function() .py_method(dss1, "subs", .py_dict(list(s), list(5L))))
})
.py_try({
.py_rec("lti:test_SymPy_substitution_functions:3431", function() .py_method(dss2, "subs", .py_dict(list(p), list(1L))))
})
.py_try({
.py_rec("lti:test_SymPy_substitution_functions:3436", function() .py_method(ss1, "xreplace", .py_dict(list(s), list(p))))
})
.py_try({
.py_rec("lti:test_SymPy_substitution_functions:3438", function() .py_method(ss2, "xreplace", .py_dict(list(s, p), list(a, b))))
})
.py_try({
.py_rec("lti:test_SymPy_substitution_functions:3441", function() .py_method(dss1, "xreplace", .py_dict(list(s), list(p))))
})
.py_try({
.py_rec("lti:test_SymPy_substitution_functions:3444", function() .py_method(dss2, "xreplace", .py_dict(list(s, p), list(a, b))))
})
.py_try({
p1 <- .py_add(.py_mul(a1, s), a0)
})
.py_try({
p2 <- .py_add(.py_add(.py_mul(b2, .py_pow(s, 2L)), .py_mul(b1, s)), b0)
})
.py_try({
G <- StateSpace(.py_Matrix(list(p1)), .py_Matrix(list(p2)))
})
.py_try({
dG <- DiscreteStateSpace(.py_Matrix(list(p1)), .py_Matrix(list(p2)), sampling_time = .py_flt(0.3))
})
.py_try({
expect <- StateSpace(.py_Matrix(list(list(.py_add(.py_mul(2L, s), 1L)))), .py_Matrix(list(list(.py_add(.py_add(.py_mul(5L, .py_pow(s, 2L)), .py_mul(4L, s)), 3L)))), .py_Matrix(list(list(0L))), .py_Matrix(list(list(0L))))
})
.py_try({
expect_ <- StateSpace(.py_Matrix(list(list(.py_add(.py_mul(.py_flt(2.0), s), .py_flt(1.0))))), .py_Matrix(list(list(.py_add(.py_add(.py_mul(.py_flt(5.0), .py_pow(s, 2L)), .py_mul(.py_flt(4.0), s)), .py_flt(3.0))))), .py_Matrix(list(list(0L))), .py_Matrix(list(list(0L))))
})
.py_try({
dexpect <- DiscreteStateSpace(.py_Matrix(list(list(.py_add(.py_mul(2L, s), 1L)))), .py_Matrix(list(list(.py_add(.py_add(.py_mul(5L, .py_pow(s, 2L)), .py_mul(4L, s)), 3L)))), .py_Matrix(list(list(0L))), .py_Matrix(list(list(0L))), .py_flt(0.3))
})
.py_try({
dexpect_ <- DiscreteStateSpace(.py_Matrix(list(list(.py_add(.py_mul(.py_flt(2.0), s), .py_flt(1.0))))), .py_Matrix(list(list(.py_add(.py_add(.py_mul(.py_flt(5.0), .py_pow(s, 2L)), .py_mul(.py_flt(4.0), s)), .py_flt(3.0))))), .py_Matrix(list(list(0L))), .py_Matrix(list(list(0L))), .py_flt(0.3))
})
.py_try({
.py_rec("lti:test_SymPy_substitution_functions:3465", function() .py_method(G, "subs", .py_dict(list(a0, a1, b0, b1, b2), list(1L, 2L, 3L, 4L, 5L))))
})
.py_try({
.py_rec("lti:test_SymPy_substitution_functions:3466", function() .py_method(.py_method(G, "subs", .py_dict(list(a0, a1, b0, b1, b2), list(1L, 2L, 3L, 4L, 5L))), "evalf"))
})
.py_try({
.py_rec("lti:test_SymPy_substitution_functions:3467", function() .py_method(expect, "evalf"))
})
.py_try({
.py_rec("lti:test_SymPy_substitution_functions:3469", function() .py_method(dG, "subs", .py_dict(list(a0, a1, b0, b1, b2), list(1L, 2L, 3L, 4L, 5L))))
})
.py_try({
.py_rec("lti:test_SymPy_substitution_functions:3470", function() .py_method(.py_method(dG, "subs", .py_dict(list(a0, a1, b0, b1, b2), list(1L, 2L, 3L, 4L, 5L))), "evalf"))
})
.py_try({
.py_rec("lti:test_SymPy_substitution_functions:3471", function() .py_method(dexpect, "evalf"))
}))---")
  expect_identical(got[["lti:test_SymPy_substitution_functions:3426"]], "StateSpace(Matrix([[5]]), Matrix([[36]]), Matrix([[24]]), Matrix([[10]]))")
  expect_identical(got[["lti:test_SymPy_substitution_functions:3427"]], "StateSpace(Matrix([[s + 1]]), Matrix([[0]]), Matrix([[1 - s^3]]), Matrix([[s - 1]]))")
  expect_identical(got[["lti:test_SymPy_substitution_functions:3429"]], "DiscreteStateSpace(Matrix([[5]]), Matrix([[36]]), Matrix([[24]]), Matrix([[10]]), 0.1)")
  expect_identical(got[["lti:test_SymPy_substitution_functions:3431"]], "DiscreteStateSpace(Matrix([[s + 1]]), Matrix([[0]]), Matrix([[1 - s^3]]), Matrix([[s - 1]]), 0.1)")
  expect_identical(got[["lti:test_SymPy_substitution_functions:3436"]], "StateSpace(Matrix([[p]]), Matrix([[(p + 1)^2]]), Matrix([[p^2 - 1]]), Matrix([[2*p]]))")
  expect_identical(got[["lti:test_SymPy_substitution_functions:3438"]], "StateSpace(Matrix([[a + b]]), Matrix([[(a + 1)*(b - 1)]]), Matrix([[-a^3 + b^3]]), Matrix([[a - b]]))")
  expect_identical(got[["lti:test_SymPy_substitution_functions:3441"]], "DiscreteStateSpace(Matrix([[p]]), Matrix([[(p + 1)^2]]), Matrix([[p^2 - 1]]), Matrix([[2*p]]), 0.1)")
  expect_identical(got[["lti:test_SymPy_substitution_functions:3444"]], "DiscreteStateSpace(Matrix([[a + b]]), Matrix([[(a + 1)*(b - 1)]]), Matrix([[-a^3 + b^3]]), Matrix([[a - b]]), 0.1)")
  expect_identical(got[["lti:test_SymPy_substitution_functions:3465"]], "StateSpace(Matrix([[2*s + 1]]), Matrix([[5*s^2 + 4*s + 3]]), Matrix([[0]]), Matrix([[0]]))")
  expect_identical(got[["lti:test_SymPy_substitution_functions:3466"]], "StateSpace(Matrix([[2.0*s + 1.0]]), Matrix([[5.0*s^2 + 4.0*s + 3.0]]), Matrix([[0]]), Matrix([[0]]))")
  expect_identical(got[["lti:test_SymPy_substitution_functions:3467"]], "StateSpace(Matrix([[2.0*s + 1.0]]), Matrix([[5.0*s^2 + 4.0*s + 3.0]]), Matrix([[0]]), Matrix([[0]]))")
  expect_identical(got[["lti:test_SymPy_substitution_functions:3469"]], "DiscreteStateSpace(Matrix([[2*s + 1]]), Matrix([[5*s^2 + 4*s + 3]]), Matrix([[0]]), Matrix([[0]]), 0.3)")
  expect_identical(got[["lti:test_SymPy_substitution_functions:3470"]], "DiscreteStateSpace(Matrix([[2.0*s + 1.0]]), Matrix([[5.0*s^2 + 4.0*s + 3.0]]), Matrix([[0]]), Matrix([[0]]), 0.3)")
  expect_identical(got[["lti:test_SymPy_substitution_functions:3471"]], "DiscreteStateSpace(Matrix([[2.0*s + 1.0]]), Matrix([[5.0*s^2 + 4.0*s + 3.0]]), Matrix([[0]]), Matrix([[0]]), 0.3)")
})

test_that("lti test_conversion matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("a, x, b, c, s, g, d,    p, k, tau, zeta, wn, T, z"); a <- .py_getitem(.tmpu_, 0L); x <- .py_getitem(.tmpu_, 1L); b <- .py_getitem(.tmpu_, 2L); c <- .py_getitem(.tmpu_, 3L); s <- .py_getitem(.tmpu_, 4L); g <- .py_getitem(.tmpu_, 5L); d <- .py_getitem(.tmpu_, 6L); p <- .py_getitem(.tmpu_, 7L); k <- .py_getitem(.tmpu_, 8L); tau <- .py_getitem(.tmpu_, 9L); zeta <- .py_getitem(.tmpu_, 10L); wn <- .py_getitem(.tmpu_, 11L); T <- .py_getitem(.tmpu_, 12L); z <- .py_getitem(.tmpu_, 13L)
})
.py_try({
.tmpu_ <- .py_symbols("a0:4, b0:5, c0:4, d0:4"); a0 <- .py_getitem(.tmpu_, 0L); a1 <- .py_getitem(.tmpu_, 1L); a2 <- .py_getitem(.tmpu_, 2L); a3 <- .py_getitem(.tmpu_, 3L); b0 <- .py_getitem(.tmpu_, 4L); b1 <- .py_getitem(.tmpu_, 5L); b2 <- .py_getitem(.tmpu_, 6L); b3 <- .py_getitem(.tmpu_, 7L); b4 <- .py_getitem(.tmpu_, 8L); c0 <- .py_getitem(.tmpu_, 9L); c1 <- .py_getitem(.tmpu_, 10L); c2 <- .py_getitem(.tmpu_, 11L); c3 <- .py_getitem(.tmpu_, 12L); d0 <- .py_getitem(.tmpu_, 13L); d1 <- .py_getitem(.tmpu_, 14L); d2 <- .py_getitem(.tmpu_, 15L); d3 <- .py_getitem(.tmpu_, 16L)
})
.py_try({
TF1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
TF2 <- TransferFunction(k, 1L, s)
})
.py_try({
TF3 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
TF4 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(tau, s), 1L), s)
})
.py_try({
TF5 <- TransferFunction(k, p, s)
})
.py_try({
TF6 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(k, s), p), s)
})
p_are_floats_equals <- function(n1, n2) {
  "\n    Helper function to compare two floats.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return(.py_method(Float(n1), "is_same", Float(n2), .py_isclose))
}
p_are_complex_equals <- function(n1, n2) {
  "\n    Helper function to compare two complex numbers.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return((.py_truthy(p_are_floats_equals(re(n1), re(n2))) && .py_truthy(p_are_floats_equals(im(n1), im(n2)))))
}
.py_try({
A1 <- .py_Matrix(list(list(-5L, -1L), list(3L, -1L)))
})
.py_try({
B1 <- .py_Matrix(list(2L, 5L))
})
.py_try({
C1 <- .py_Matrix(list(list(1L, 2L)))
})
.py_try({
D1 <- .py_Matrix(list(0L))
})
.py_try({
H1 <- StateSpace(A1, B1, C1, D1)
})
.py_try({
H3 <- StateSpace(.py_Matrix(list(list(a0, a1), list(a2, a3))), B = .py_Matrix(list(list(b1), list(b2))), C = .py_Matrix(list(list(c1, c2))))
})
.py_try({
tm1 <- .py_method(H1, "rewrite", TransferFunction)
})
.py_try({
.py_rec("lti:test_conversion:3482", function() .py_method(H1, "rewrite", DiscreteTransferFunction))
})
.py_try({
tm2 <- .py_method(.py_neg(H1), "rewrite", TransferFunction)
})
.py_try({
tf1 <- .py_getitem(.py_getitem(tm1, 0L), 0L)
})
.py_try({
tf2 <- .py_getitem(.py_getitem(tm2, 0L), 0L)
})
.py_try({
.py_rec("lti:test_conversion:3488", function() tf1)
})
.py_try({
.py_rec("lti:test_conversion:3489", function() .py_attr(tf2, "num"))
})
.py_try({
.py_rec("lti:test_conversion:3490", function() .py_attr(tf2, "den"))
})
.py_try({
A2 <- .py_Matrix(list(list(.py_flt(-1.5), -2L, 3L), list(1L, 0L, 1L), list(2L, 1L, 1L)))
})
.py_try({
B2 <- .py_Matrix(list(list(.py_flt(0.5), 0L, 1L), list(0L, 1L, 2L), list(2L, 2L, 3L)))
})
.py_try({
C2 <- .py_Matrix(list(list(0L, 1L, 0L), list(0L, 2L, 1L), list(1L, 0L, 2L)))
})
.py_try({
D2 <- .py_Matrix(list(list(2L, 2L, 0L), list(1L, 1L, 1L), list(3L, 2L, 1L)))
})
.py_try({
H2 <- StateSpace(A2, B2, C2, D2)
})
.py_try({
tm3 <- .py_method(H2, "rewrite", TransferFunction)
})
.py_try({
.py_rec("lti:test_conversion:3504", function() .py_method(H3, "rewrite", TransferFunction))
})
.py_try({
SS <- .py_method(TF1, "rewrite", StateSpace)
})
.py_try({
.py_rec("lti:test_conversion:3515", function() .py_getitem(.py_getitem(.py_method(SS, "rewrite", TransferFunction), 0L), 0L))
})
.py_try({
SS1 <- .py_method(TF4, "rewrite", StateSpace)
})
.py_try({
.py_rec("lti:test_conversion:3518", function() SS1)
})
.py_try({
.py_rec("lti:test_conversion:3521", function() .py_getitem(.py_getitem(.py_method(SS1, "rewrite", TransferFunction), 0L), 0L))
})
.py_try({
SS2 <- .py_method(TF5, "rewrite", StateSpace)
})
.py_try({
.py_rec("lti:test_conversion:3524", function() SS2)
})
.py_try({
.py_rec("lti:test_conversion:3526", function() .py_getitem(.py_getitem(.py_method(SS2, "rewrite", TransferFunction), 0L), 0L))
})
.py_try({
SS3 <- .py_method(TF6, "rewrite", StateSpace)
})
.py_try({
.py_rec("lti:test_conversion:3529", function() SS3)
})
.py_try({
.py_rec("lti:test_conversion:3533", function() .py_method(TF1, "rewrite", DiscreteStateSpace))
})
.py_try({
.py_rec("lti:test_conversion:3536", function() .py_method(TransferFunction(.py_add(.py_sub(.py_add(.py_mul(b, .py_pow(s, 2L)), .py_pow(p, 2L)), .py_mul(a, p)), s), .py_sub(b, .py_pow(p, 2L)), s), "rewrite", StateSpace))
})
.py_try({
DH1 <- DiscreteStateSpace(A1, B1, C1, D1, .py_flt(0.1))
})
.py_try({
DH3 <- DiscreteStateSpace(.py_Matrix(list(list(a0, a1), list(a2, a3))), B = .py_Matrix(list(list(b1), list(b2))), C = .py_Matrix(list(list(c1, c2))), sampling_time = T)
})
.py_try({
dtm1 <- .py_method(DH1, "rewrite", DiscreteTransferFunction)
})
.py_try({
.py_rec("lti:test_conversion:3544", function() .py_method(DH1, "rewrite", TransferFunction))
})
.py_try({
dtm2 <- .py_method(.py_neg(DH1), "rewrite", DiscreteTransferFunction)
})
.py_try({
dtf1 <- .py_getitem(.py_getitem(dtm1, 0L), 0L)
})
.py_try({
dtf2 <- .py_getitem(.py_getitem(dtm2, 0L), 0L)
})
.py_try({
.py_rec("lti:test_conversion:3550", function() dtf1)
})
.py_try({
.py_rec("lti:test_conversion:3551", function() .py_attr(dtf2, "num"))
})
.py_try({
.py_rec("lti:test_conversion:3552", function() .py_attr(dtf2, "den"))
})
.py_try({
.py_rec("lti:test_conversion:3553", function() .py_attr(dtf1, "sampling_time"))
})
.py_try({
.py_rec("lti:test_conversion:3554", function() .py_attr(dtf2, "sampling_time"))
})
.py_try({
DH2 <- DiscreteStateSpace(A2, B2, C2, D2, T)
})
.py_try({
dtm3 <- .py_method(DH2, "rewrite", DiscreteTransferFunction)
})
.py_try({
.py_rec("lti:test_conversion:3564", function() .py_method(DH3, "rewrite", DiscreteTransferFunction))
})
.py_try({
dtf1 <- DiscreteTransferFunction(.py_add(z, 5L), .py_add(.py_add(.py_add(.py_pow(z, 3L), .py_mul(2L, .py_pow(z, 2L))), .py_mul(4L, z)), 3L), z, .py_flt(0.1))
})
.py_try({
DSS <- .py_method(dtf1, "rewrite", DiscreteStateSpace)
})
.py_try({
A <- .py_Matrix(list(list(0L, 1L, 0L), list(0L, 0L, 1L), list(-3L, -4L, -2L)))
})
.py_try({
B <- .py_Matrix(list(0L, 0L, 1L))
})
.py_try({
C <- .py_Matrix(list(list(5L, 1L, 0L)))
})
.py_try({
D <- .py_Matrix(list(list(0L)))
})
.py_try({
.py_rec("lti:test_conversion:3574", function() DSS)
})
.py_try({
.py_rec("lti:test_conversion:3576", function() .py_getitem(.py_getitem(.py_method(DSS, "rewrite", DiscreteTransferFunction), 0L), 0L))
})
.py_try({
.py_rec("lti:test_conversion:3579", function() .py_method(dtf1, "rewrite", StateSpace))
})
.py_try({
.py_rec("lti:test_conversion:3582", function() .py_method(DiscreteTransferFunction(.py_add(.py_add(.py_mul(b2, .py_pow(z, 2L)), .py_mul(b1, z)), b0), .py_add(z, a0), z), "rewrite", DiscreteStateSpace))
}))---")
  expect_identical(got[["lti:test_conversion:3482"]], "ERR TypeError: \\n            The continuous state space model cannot be rewritten as a\\n            discrete-time transfer function model.\\n            ")
  expect_identical(got[["lti:test_conversion:3488"]], "TransferFunction(12*s + 59, s^2 + 6*s + 8, s)")
  expect_identical(got[["lti:test_conversion:3489"]], "-12*s - 59")
  expect_identical(got[["lti:test_conversion:3490"]], "s^2 + 6*s + 8")
  expect_identical(got[["lti:test_conversion:3504"]], "[[TransferFunction(-c1*(a1*b2 - a3*b1 + b1*s) - c2*(-a0*b2 + a2*b1 + b2*s), -a0*a3 + a0*s + a1*a2 + a3*s - s^2, s)]]")
  expect_identical(got[["lti:test_conversion:3515"]], "TransferFunction(1, s^2 + 2*s*wn*zeta + wn^2, s)")
  expect_identical(got[["lti:test_conversion:3518"]], "StateSpace(Matrix([[-1/tau]]), Matrix([[1]]), Matrix([[-k/tau^2 + p/tau]]), Matrix([[k/tau]]))")
  expect_identical(got[["lti:test_conversion:3521"]], "TransferFunction(k*s + p, s*tau + 1, s)")
  expect_identical(got[["lti:test_conversion:3524"]], "StateSpace(Matrix([[0]]), Matrix([[0]]), Matrix([[0]]), Matrix([[k/p]]))")
  expect_identical(got[["lti:test_conversion:3526"]], "TransferFunction(k, p, s)")
  expect_identical(got[["lti:test_conversion:3529"]], "StateSpace(Matrix([[0]]), Matrix([[0]]), Matrix([[0]]), Matrix([[1]]))")
  expect_identical(got[["lti:test_conversion:3533"]], "ERR TypeError: \\n            The continuous transfer function model cannot be rewritten as a\\n            discrete-time state space model.\\n            ")
  expect_identical(got[["lti:test_conversion:3536"]], "ERR ValueError: Transfer Function must be proper.")
  expect_identical(got[["lti:test_conversion:3544"]], "ERR TypeError: \\n            The discrete state space model cannot be rewritten as a\\n            continuous-time transfer function model.\\n            ")
  expect_identical(got[["lti:test_conversion:3550"]], "DiscreteTransferFunction(12*z + 59, z^2 + 6*z + 8, z, 0.1)")
  expect_identical(got[["lti:test_conversion:3551"]], "-12*z - 59")
  expect_identical(got[["lti:test_conversion:3552"]], "z^2 + 6*z + 8")
  expect_identical(got[["lti:test_conversion:3553"]], "0.100000000000000")
  expect_identical(got[["lti:test_conversion:3554"]], "0.100000000000000")
  expect_identical(got[["lti:test_conversion:3564"]], "[[DiscreteTransferFunction(-c1*(a1*b2 - a3*b1 + b1*z) - c2*(-a0*b2 + a2*b1 + b2*z), -a0*a3 + a0*z + a1*a2 + a3*z - z^2, z, T)]]")
  expect_identical(got[["lti:test_conversion:3574"]], "DiscreteStateSpace(Matrix([\\n[ 0,  1,  0],\\n[ 0,  0,  1],\\n[-3, -4, -2]]), Matrix([\\n[0],\\n[0],\\n[1]]), Matrix([[5, 1, 0]]), Matrix([[0]]), 0.1)")
  expect_identical(got[["lti:test_conversion:3576"]], "DiscreteTransferFunction(z + 5, z^3 + 2*z^2 + 4*z + 3, z, 0.1)")
  expect_identical(got[["lti:test_conversion:3579"]], "ERR TypeError: \\n            The discrete transfer function model cannot be rewritten as a\\n            continuous-time state space model.\\n            ")
  expect_identical(got[["lti:test_conversion:3582"]], "ERR ValueError: Transfer Function must be proper.")
})

test_that("lti test_StateSpace_functions matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("a, x, b, c, s, g, d,    p, k, tau, zeta, wn, T, z"); a <- .py_getitem(.tmpu_, 0L); x <- .py_getitem(.tmpu_, 1L); b <- .py_getitem(.tmpu_, 2L); c <- .py_getitem(.tmpu_, 3L); s <- .py_getitem(.tmpu_, 4L); g <- .py_getitem(.tmpu_, 5L); d <- .py_getitem(.tmpu_, 6L); p <- .py_getitem(.tmpu_, 7L); k <- .py_getitem(.tmpu_, 8L); tau <- .py_getitem(.tmpu_, 9L); zeta <- .py_getitem(.tmpu_, 10L); wn <- .py_getitem(.tmpu_, 11L); T <- .py_getitem(.tmpu_, 12L); z <- .py_getitem(.tmpu_, 13L)
})
.py_try({
.tmpu_ <- .py_symbols("a0:4, b0:5, c0:4, d0:4"); a0 <- .py_getitem(.tmpu_, 0L); a1 <- .py_getitem(.tmpu_, 1L); a2 <- .py_getitem(.tmpu_, 2L); a3 <- .py_getitem(.tmpu_, 3L); b0 <- .py_getitem(.tmpu_, 4L); b1 <- .py_getitem(.tmpu_, 5L); b2 <- .py_getitem(.tmpu_, 6L); b3 <- .py_getitem(.tmpu_, 7L); b4 <- .py_getitem(.tmpu_, 8L); c0 <- .py_getitem(.tmpu_, 9L); c1 <- .py_getitem(.tmpu_, 10L); c2 <- .py_getitem(.tmpu_, 11L); c3 <- .py_getitem(.tmpu_, 12L); d0 <- .py_getitem(.tmpu_, 13L); d1 <- .py_getitem(.tmpu_, 14L); d2 <- .py_getitem(.tmpu_, 15L); d3 <- .py_getitem(.tmpu_, 16L)
})
.py_try({
TF1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
TF2 <- TransferFunction(k, 1L, s)
})
.py_try({
TF3 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
TF4 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(tau, s), 1L), s)
})
.py_try({
TF5 <- TransferFunction(k, p, s)
})
.py_try({
TF6 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(k, s), p), s)
})
p_are_floats_equals <- function(n1, n2) {
  "\n    Helper function to compare two floats.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return(.py_method(Float(n1), "is_same", Float(n2), .py_isclose))
}
p_are_complex_equals <- function(n1, n2) {
  "\n    Helper function to compare two complex numbers.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return((.py_truthy(p_are_floats_equals(re(n1), re(n2))) && .py_truthy(p_are_floats_equals(im(n1), im(n2)))))
}
.py_try({
A_mat <- .py_Matrix(list(list(.py_flt(-1.5), -2L), list(1L, 0L)))
})
.py_try({
B_mat <- .py_Matrix(list(.py_flt(0.5), 0L))
})
.py_try({
C_mat <- .py_Matrix(list(list(0L, 1L)))
})
.py_try({
D_mat <- .py_Matrix(list(1L))
})
.py_try({
SS1 <- StateSpace(A_mat, B_mat, C_mat, D_mat)
})
.py_try({
SS2 <- StateSpace(.py_Matrix(list(list(1L, 1L), list(4L, -2L))), .py_Matrix(list(list(0L, 1L), list(0L, 2L))), .py_Matrix(list(list(-1L, 1L), list(1L, -1L))))
})
.py_try({
SS3 <- StateSpace(.py_Matrix(list(list(1L, 1L), list(4L, -2L))), .py_Matrix(list(list(1L, -1L), list(1L, -1L))))
})
.py_try({
SS4 <- StateSpace(.py_Matrix(list(list(a0, a1), list(a2, a3))), .py_Matrix(list(list(b1), list(b2))), .py_Matrix(list(list(c1, c2))))
})
.py_try({
.py_rec("lti:test_StateSpace_functions:3660", function() .py_method(SS1, "is_observable"))
})
.py_try({
.py_rec("lti:test_StateSpace_functions:3661", function() .py_method(SS2, "is_observable"))
})
.py_try({
.py_rec("lti:test_StateSpace_functions:3662", function() .py_method(SS1, "observability_matrix"))
})
.py_try({
.py_rec("lti:test_StateSpace_functions:3663", function() .py_method(SS2, "observability_matrix"))
})
.py_try({
.py_rec("lti:test_StateSpace_functions:3664", function() .py_method(SS1, "observable_subspace"))
})
.py_try({
.py_rec("lti:test_StateSpace_functions:3665", function() .py_method(SS2, "observable_subspace"))
})
.py_try({
.py_rec("lti:test_StateSpace_functions:3666", function() .py_method(SS4, "observable_subspace"))
})
.py_try({
.py_rec("lti:test_StateSpace_functions:3667", function() .py_method(SS1, "unobservable_subspace"))
})
.py_try({
.py_rec("lti:test_StateSpace_functions:3668", function() .py_method(SS2, "unobservable_subspace"))
})
.py_try({
.py_rec("lti:test_StateSpace_functions:3669", function() .py_method(SS4, "unobservable_subspace"))
})
.py_try({
Qo <- .py_method(.py_method(SS4, "observability_matrix"), "subs", list(.py_tup(a0, 0L), .py_tup(a1, -6L), .py_tup(a2, 1L), .py_tup(a3, -5L), .py_tup(c1, 0L), .py_tup(c2, 1L)))
})
.py_try({
.py_rec("lti:test_StateSpace_functions:3672", function() Qo)
})
.py_try({
ss_obs <- .py_method(StateSpace(.py_Matrix(list(list(1L, 0L, 1L), list(0L, 0L, 0L), list(0L, 0L, -2L))), .py_Matrix(list(1L, 1L, 0L)), .py_attr(.py_Matrix(list(1L, 1L, Rational(1L, 3L))), "T")), "to_observable_form")
})
.py_try({
A_obs <- .py_getitem(.py_attr(ss_obs, "A"), .py_tup(.py_slice(NULL, 2L, NULL), .py_slice(NULL, 2L, NULL)))
})
.py_try({
A_nobs <- .py_getitem(.py_attr(ss_obs, "A"), .py_tup(.py_slice(2L, NULL, NULL), .py_slice(2L, NULL, NULL)))
})
.py_try({
.py_rec("lti:test_StateSpace_functions:3679", function() .py_method(A_nobs, "eigenvals"))
})
.py_try({
.py_rec("lti:test_StateSpace_functions:3682", function() .py_method(SS1, "is_controllable"))
})
.py_try({
.py_rec("lti:test_StateSpace_functions:3683", function() .py_method(SS3, "is_controllable"))
})
.py_try({
.py_rec("lti:test_StateSpace_functions:3684", function() .py_method(SS1, "controllability_matrix"))
})
.py_try({
.py_rec("lti:test_StateSpace_functions:3685", function() .py_method(SS3, "controllability_matrix"))
})
.py_try({
.py_rec("lti:test_StateSpace_functions:3686", function() .py_method(SS4, "controllability_matrix"))
})
.py_try({
.py_rec("lti:test_StateSpace_functions:3688", function() .py_method(SS1, "controllable_subspace"))
})
.py_try({
.py_rec("lti:test_StateSpace_functions:3689", function() .py_method(SS3, "controllable_subspace"))
})
.py_try({
.py_rec("lti:test_StateSpace_functions:3690", function() .py_method(SS1, "uncontrollable_subspace"))
})
.py_try({
.py_rec("lti:test_StateSpace_functions:3691", function() .py_method(SS3, "uncontrollable_subspace"))
})
.py_try({
.py_rec("lti:test_StateSpace_functions:3692", function() .py_method(SS4, "uncontrollable_subspace"))
})
.py_try({
Qc <- .py_method(.py_method(SS4, "controllability_matrix"), "subs", list(.py_tup(a0, 0L), .py_tup(a1, 1L), .py_tup(a2, -6L), .py_tup(a3, -5L), .py_tup(b1, 0L), .py_tup(b2, 1L)))
})
.py_try({
.py_rec("lti:test_StateSpace_functions:3695", function() Qc)
})
.py_try({
ss_contr <- .py_method(StateSpace(.py_Matrix(list(list(1L, 0L, 1L), list(0L, 0L, 0L), list(0L, 0L, -2L))), .py_Matrix(list(1L, 1L, 0L)), .py_attr(.py_Matrix(list(1L, 1L, 0L)), "T")), "to_controllable_form")
})
.py_try({
A_contr <- .py_getitem(.py_attr(ss_contr, "A"), .py_tup(.py_slice(NULL, 2L, NULL), .py_slice(NULL, 2L, NULL)))
})
.py_try({
A_ncontr <- .py_getitem(.py_attr(ss_contr, "A"), .py_tup(.py_slice(2L, NULL, NULL), .py_slice(2L, NULL, NULL)))
})
.py_try({
.py_rec("lti:test_StateSpace_functions:3700", function() .py_method(A_ncontr, "eigenvals"))
})
.py_try({
A1 <- .py_Matrix(list(list(0L, 1L), list(1L, 0L)))
})
.py_try({
B1 <- .py_Matrix(list(list(0L), list(1L)))
})
.py_try({
C1 <- .py_Matrix(list(list(0L, 1L)))
})
.py_try({
D1 <- .py_Matrix(list(list(0L)))
})
.py_try({
ss1 <- StateSpace(A1, B1, C1, D1)
})
.py_try({
ss2 <- StateSpace(.py_Matrix(list(list(1L, 0L), list(0L, 1L))), .py_Matrix(list(list(1L), list(0L))), .py_Matrix(list(list(1L, 0L))), .py_Matrix(list(list(1L))))
})
.py_try({
ss3 <- .py_method(ss1, "append", ss2)
})
.py_try({
ss4 <- .py_method(SS4, "append", ss1)
})
.py_try({
.py_rec("lti:test_StateSpace_functions:3712", function() .py_attr(ss3, "num_states"))
})
.py_try({
.py_rec("lti:test_StateSpace_functions:3713", function() .py_attr(ss3, "num_inputs"))
})
.py_try({
.py_rec("lti:test_StateSpace_functions:3714", function() .py_attr(ss3, "num_outputs"))
})
.py_try({
.py_rec("lti:test_StateSpace_functions:3715", function() .py_attr(ss3, "state_matrix"))
})
.py_try({
.py_rec("lti:test_StateSpace_functions:3716", function() .py_attr(ss3, "input_matrix"))
})
.py_try({
.py_rec("lti:test_StateSpace_functions:3717", function() .py_attr(ss3, "output_matrix"))
})
.py_try({
.py_rec("lti:test_StateSpace_functions:3718", function() .py_attr(ss3, "feedforward_matrix"))
})
.py_try({
.py_rec("lti:test_StateSpace_functions:3721", function() .py_attr(ss4, "num_states"))
})
.py_try({
.py_rec("lti:test_StateSpace_functions:3722", function() .py_attr(ss4, "num_inputs"))
})
.py_try({
.py_rec("lti:test_StateSpace_functions:3723", function() .py_attr(ss4, "num_outputs"))
})
.py_try({
.py_rec("lti:test_StateSpace_functions:3724", function() .py_attr(ss4, "state_matrix"))
})
.py_try({
.py_rec("lti:test_StateSpace_functions:3725", function() .py_attr(ss4, "input_matrix"))
})
.py_try({
.py_rec("lti:test_StateSpace_functions:3726", function() .py_attr(ss4, "output_matrix"))
})
.py_try({
.py_rec("lti:test_StateSpace_functions:3727", function() .py_attr(ss4, "feedforward_matrix"))
}))---")
  expect_identical(got[["lti:test_StateSpace_functions:3660"]], "True")
  expect_identical(got[["lti:test_StateSpace_functions:3661"]], "False")
  expect_identical(got[["lti:test_StateSpace_functions:3662"]], "Matrix([[0, 1], [1, 0]])")
  expect_identical(got[["lti:test_StateSpace_functions:3663"]], "Matrix([[-1, 1], [1, -1], [3, -3], [-3, 3]])")
  expect_identical(got[["lti:test_StateSpace_functions:3664"]], "[Matrix([[1], [0]]), Matrix([[0], [1]])]")
  expect_identical(got[["lti:test_StateSpace_functions:3665"]], "[Matrix([[-1], [1]])]")
  expect_identical(got[["lti:test_StateSpace_functions:3666"]], "ERR NotImplementedError: Unbservable subspace cannot be determined for symbolic matrices.")
  expect_identical(got[["lti:test_StateSpace_functions:3667"]], "[]")
  expect_identical(got[["lti:test_StateSpace_functions:3668"]], "[Matrix([[1], [1]])]")
  expect_identical(got[["lti:test_StateSpace_functions:3669"]], "ERR NotImplementedError: Unbservable subspace cannot be determined for symbolic matrices.")
  expect_identical(got[["lti:test_StateSpace_functions:3672"]], "Matrix([[0, 1], [1, -5]])")
  expect_identical(got[["lti:test_StateSpace_functions:3679"]], "{-2: 1}")
  expect_identical(got[["lti:test_StateSpace_functions:3682"]], "True")
  expect_identical(got[["lti:test_StateSpace_functions:3683"]], "False")
  expect_identical(got[["lti:test_StateSpace_functions:3684"]], "Matrix([[0.500000000000000, -0.750000000000000], [0, 0.500000000000000]])")
  expect_identical(got[["lti:test_StateSpace_functions:3685"]], "Matrix([[1, -1, 2, -2], [1, -1, 2, -2]])")
  expect_identical(got[["lti:test_StateSpace_functions:3686"]], "Matrix([[b1, a0*b1 + a1*b2], [b2, a2*b1 + a3*b2]])")
  expect_identical(got[["lti:test_StateSpace_functions:3688"]], "[Matrix([[0.500000000000000], [0]]), Matrix([[-0.750000000000000], [0.500000000000000]])]")
  expect_identical(got[["lti:test_StateSpace_functions:3689"]], "[Matrix([[1], [1]])]")
  expect_identical(got[["lti:test_StateSpace_functions:3690"]], "[]")
  expect_identical(got[["lti:test_StateSpace_functions:3691"]], "[Matrix([[-1], [1]])]")
  expect_identical(got[["lti:test_StateSpace_functions:3692"]], "ERR NotImplementedError: Uncontrollable subspace cannot be determined for symbolic matrices.")
  expect_identical(got[["lti:test_StateSpace_functions:3695"]], "Matrix([[0, 1], [1, -5]])")
  expect_identical(got[["lti:test_StateSpace_functions:3700"]], "{-2: 1}")
  expect_identical(got[["lti:test_StateSpace_functions:3712"]], "4")
  expect_identical(got[["lti:test_StateSpace_functions:3713"]], "2")
  expect_identical(got[["lti:test_StateSpace_functions:3714"]], "2")
  expect_identical(got[["lti:test_StateSpace_functions:3715"]], "Matrix([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])")
  expect_identical(got[["lti:test_StateSpace_functions:3716"]], "Matrix([[0, 0], [1, 0], [0, 1], [0, 0]])")
  expect_identical(got[["lti:test_StateSpace_functions:3717"]], "Matrix([[0, 1, 0, 0], [0, 0, 1, 0]])")
  expect_identical(got[["lti:test_StateSpace_functions:3718"]], "Matrix([[0, 0], [0, 1]])")
  expect_identical(got[["lti:test_StateSpace_functions:3721"]], "4")
  expect_identical(got[["lti:test_StateSpace_functions:3722"]], "2")
  expect_identical(got[["lti:test_StateSpace_functions:3723"]], "2")
  expect_identical(got[["lti:test_StateSpace_functions:3724"]], "Matrix([[a0, a1, 0, 0], [a2, a3, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]])")
  expect_identical(got[["lti:test_StateSpace_functions:3725"]], "Matrix([[b1, 0], [b2, 0], [0, 0], [0, 1]])")
  expect_identical(got[["lti:test_StateSpace_functions:3726"]], "Matrix([[c1, c2, 0, 0], [0, 0, 0, 1]])")
  expect_identical(got[["lti:test_StateSpace_functions:3727"]], "Matrix([[0, 0], [0, 0]])")
})

test_that("lti test_StateSpace_series matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("a, x, b, c, s, g, d,    p, k, tau, zeta, wn, T, z"); a <- .py_getitem(.tmpu_, 0L); x <- .py_getitem(.tmpu_, 1L); b <- .py_getitem(.tmpu_, 2L); c <- .py_getitem(.tmpu_, 3L); s <- .py_getitem(.tmpu_, 4L); g <- .py_getitem(.tmpu_, 5L); d <- .py_getitem(.tmpu_, 6L); p <- .py_getitem(.tmpu_, 7L); k <- .py_getitem(.tmpu_, 8L); tau <- .py_getitem(.tmpu_, 9L); zeta <- .py_getitem(.tmpu_, 10L); wn <- .py_getitem(.tmpu_, 11L); T <- .py_getitem(.tmpu_, 12L); z <- .py_getitem(.tmpu_, 13L)
})
.py_try({
.tmpu_ <- .py_symbols("a0:4, b0:5, c0:4, d0:4"); a0 <- .py_getitem(.tmpu_, 0L); a1 <- .py_getitem(.tmpu_, 1L); a2 <- .py_getitem(.tmpu_, 2L); a3 <- .py_getitem(.tmpu_, 3L); b0 <- .py_getitem(.tmpu_, 4L); b1 <- .py_getitem(.tmpu_, 5L); b2 <- .py_getitem(.tmpu_, 6L); b3 <- .py_getitem(.tmpu_, 7L); b4 <- .py_getitem(.tmpu_, 8L); c0 <- .py_getitem(.tmpu_, 9L); c1 <- .py_getitem(.tmpu_, 10L); c2 <- .py_getitem(.tmpu_, 11L); c3 <- .py_getitem(.tmpu_, 12L); d0 <- .py_getitem(.tmpu_, 13L); d1 <- .py_getitem(.tmpu_, 14L); d2 <- .py_getitem(.tmpu_, 15L); d3 <- .py_getitem(.tmpu_, 16L)
})
.py_try({
TF1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
TF2 <- TransferFunction(k, 1L, s)
})
.py_try({
TF3 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
TF4 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(tau, s), 1L), s)
})
.py_try({
TF5 <- TransferFunction(k, p, s)
})
.py_try({
TF6 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(k, s), p), s)
})
p_are_floats_equals <- function(n1, n2) {
  "\n    Helper function to compare two floats.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return(.py_method(Float(n1), "is_same", Float(n2), .py_isclose))
}
p_are_complex_equals <- function(n1, n2) {
  "\n    Helper function to compare two complex numbers.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return((.py_truthy(p_are_floats_equals(re(n1), re(n2))) && .py_truthy(p_are_floats_equals(im(n1), im(n2)))))
}
.py_try({
a1 <- .py_Matrix(list(list(0L, 1L), list(1L, 0L)))
})
.py_try({
b1 <- .py_Matrix(list(list(0L), list(1L)))
})
.py_try({
c1 <- .py_Matrix(list(list(0L, 1L)))
})
.py_try({
d1 <- .py_Matrix(list(list(0L)))
})
.py_try({
a2 <- .py_Matrix(list(list(1L, 0L), list(0L, 1L)))
})
.py_try({
b2 <- .py_Matrix(list(list(1L), list(0L)))
})
.py_try({
c2 <- .py_Matrix(list(list(1L, 0L)))
})
.py_try({
d2 <- .py_Matrix(list(list(1L)))
})
.py_try({
ss1 <- StateSpace(a1, b1, c1, d1)
})
.py_try({
ss2 <- StateSpace(a2, b2, c2, d2)
})
.py_try({
tf1 <- TransferFunction(s, .py_add(s, 1L), s)
})
.py_try({
ser1 <- Series(ss1, ss2)
})
.py_try({
.py_rec("lti:test_StateSpace_series:3744", function() ser1)
})
.py_try({
.py_rec("lti:test_StateSpace_series:3753", function() .py_method(ser1, "doit"))
})
.py_try({
.py_rec("lti:test_StateSpace_series:3767", function() .py_attr(ser1, "num_inputs"))
})
.py_try({
.py_rec("lti:test_StateSpace_series:3768", function() .py_attr(ser1, "num_outputs"))
})
.py_try({
.py_rec("lti:test_StateSpace_series:3769", function() .py_method(ser1, "rewrite", TransferFunction))
})
.py_try({
.py_rec("lti:test_StateSpace_series:3770", function() .py_method(ser1, "rewrite", DiscreteTransferFunction))
})
.py_try({
ser2 <- Series(ss1)
})
.py_try({
ser3 <- Series(ser2, ss2)
})
.py_try({
.py_rec("lti:test_StateSpace_series:3774", function() .py_method(ser3, "doit"))
})
.py_try({
ss_issue <- StateSpace(a2, .py_add(b1, b2), .py_add(c1, c2), d1)
})
.py_try({
.py_rec("lti:test_StateSpace_series:3778", function() .py_method(Series(ss_issue, ss_issue, ss_issue), "doit"))
})
.py_try({
ser_tf <- Series(tf1, ss1)
})
.py_try({
.py_rec("lti:test_StateSpace_series:3798", function() ser_tf)
})
.py_try({
.py_rec("lti:test_StateSpace_series:3803", function() .py_method(ser_tf, "doit"))
})
.py_try({
.py_rec("lti:test_StateSpace_series:3814", function() .py_method(ser_tf, "rewrite", TransferFunction))
})
.py_try({
a3 <- .py_Matrix(list(list(4L, 1L), list(2L, -3L)))
})
.py_try({
b3 <- .py_Matrix(list(list(5L, 2L), list(-3L, -3L)))
})
.py_try({
c3 <- .py_Matrix(list(list(2L, -4L), list(0L, 1L)))
})
.py_try({
d3 <- .py_Matrix(list(list(3L, 2L), list(1L, -1L)))
})
.py_try({
a4 <- .py_Matrix(list(list(-3L, 4L, 2L), list(-1L, -3L, 0L), list(2L, 5L, 3L)))
})
.py_try({
b4 <- .py_Matrix(list(list(1L, 4L), list(-3L, -3L), list(-2L, 1L)))
})
.py_try({
c4 <- .py_Matrix(list(list(4L, 2L, -3L), list(1L, 4L, 3L)))
})
.py_try({
d4 <- .py_Matrix(list(list(-2L, 4L), list(0L, 1L)))
})
.py_try({
ss3 <- StateSpace(a3, b3, c3, d3)
})
.py_try({
ss4 <- StateSpace(a4, b4, c4, d4)
})
.py_try({
ser4 <- MIMOSeries(ss3, ss4)
})
.py_try({
.py_rec("lti:test_StateSpace_series:3828", function() ser4)
})
.py_try({
.py_rec("lti:test_StateSpace_series:3847", function() .py_method(ser4, "doit"))
})
.py_try({
.py_rec("lti:test_StateSpace_series:3866", function() .py_attr(ser4, "num_inputs"))
})
.py_try({
.py_rec("lti:test_StateSpace_series:3867", function() .py_attr(ser4, "num_outputs"))
})
.py_try({
ser5 <- MIMOSeries(ss3)
})
.py_try({
ser6 <- MIMOSeries(ser5, ss4)
})
.py_try({
.py_rec("lti:test_StateSpace_series:3870", function() .py_method(ser6, "doit"))
})
.py_try({
.py_rec("lti:test_StateSpace_series:3871", function() .py_method(ser6, "rewrite", TransferFunctionMatrix))
})
.py_try({
tf2 <- TransferFunction(1L, s, s)
})
.py_try({
tf3 <- TransferFunction(1L, .py_add(s, 1L), s)
})
.py_try({
tf4 <- TransferFunction(s, .py_add(s, 2L), s)
})
.py_try({
tfm <- TransferFunctionMatrix(list(list(tf1, tf2), list(tf3, tf4)))
})
.py_try({
ser6 <- MIMOSeries(ss3, tfm)
})
.py_try({
.py_rec("lti:test_StateSpace_series:3877", function() ser6)
}))---")
  expect_identical(got[["lti:test_StateSpace_series:3744"]], "Series(StateSpace(Matrix([\\n[0, 1],\\n[1, 0]]), Matrix([\\n[0],\\n[1]]), Matrix([[0, 1]]), Matrix([[0]])), StateSpace(Matrix([\\n[1, 0],\\n[0, 1]]), Matrix([\\n[1],\\n[0]]), Matrix([[1, 0]]), Matrix([[1]])))")
  expect_identical(got[["lti:test_StateSpace_series:3753"]], "StateSpace(Matrix([\\n[0, 1, 0, 0],\\n[1, 0, 0, 0],\\n[0, 1, 1, 0],\\n[0, 0, 0, 1]]), Matrix([\\n[0],\\n[1],\\n[0],\\n[0]]), Matrix([[0, 1, 1, 0]]), Matrix([[0]]))")
  expect_identical(got[["lti:test_StateSpace_series:3767"]], "1")
  expect_identical(got[["lti:test_StateSpace_series:3768"]], "1")
  expect_identical(got[["lti:test_StateSpace_series:3769"]], "TransferFunction(s^2, s^3 - s^2 - s + 1, s)")
  expect_identical(got[["lti:test_StateSpace_series:3770"]], "ERR TypeError: \\n            The continuous state space model cannot be rewritten as a\\n            discrete-time transfer function model.\\n            ")
  expect_identical(got[["lti:test_StateSpace_series:3774"]], "StateSpace(Matrix([\\n[0, 1, 0, 0],\\n[1, 0, 0, 0],\\n[0, 1, 1, 0],\\n[0, 0, 0, 1]]), Matrix([\\n[0],\\n[1],\\n[0],\\n[0]]), Matrix([[0, 1, 1, 0]]), Matrix([[0]]))")
  expect_identical(got[["lti:test_StateSpace_series:3778"]], "StateSpace(Matrix([\\n[1, 0, 0, 0, 0, 0],\\n[0, 1, 0, 0, 0, 0],\\n[1, 1, 1, 0, 0, 0],\\n[1, 1, 0, 1, 0, 0],\\n[0, 0, 1, 1, 1, 0],\\n[0, 0, 1, 1, 0, 1]]), Matrix([\\n[1],\\n[1],\\n[0],\\n[0],\\n[0],\\n[0]]), Matrix([[0, 0, 0, 0, 1, 1]]), Matrix([[0]]))")
  expect_identical(got[["lti:test_StateSpace_series:3798"]], "Series(TransferFunction(s, s + 1, s), StateSpace(Matrix([\\n[0, 1],\\n[1, 0]]), Matrix([\\n[0],\\n[1]]), Matrix([[0, 1]]), Matrix([[0]])))")
  expect_identical(got[["lti:test_StateSpace_series:3803"]], "StateSpace(Matrix([\\n[-1, 0, 0],\\n[ 0, 0, 1],\\n[-1, 1, 0]]), Matrix([\\n[1],\\n[0],\\n[1]]), Matrix([[0, 0, 1]]), Matrix([[0]]))")
  expect_identical(got[["lti:test_StateSpace_series:3814"]], "TransferFunction(s^2, s^3 + s^2 - s - 1, s)")
  expect_identical(got[["lti:test_StateSpace_series:3828"]], "MIMOSeries(StateSpace(Matrix([\\n[4,  1],\\n[2, -3]]), Matrix([\\n[ 5,  2],\\n[-3, -3]]), Matrix([\\n[2, -4],\\n[0,  1]]), Matrix([\\n[3,  2],\\n[1, -1]])), StateSpace(Matrix([\\n[-3,  4, 2],\\n[-1, -3, 0],\\n[ 2,  5, 3]]), Matrix([\\n[ 1,  4],\\n[-3, -3],\\n[-2,  1]]), Matrix([\\n[4, 2, -3],\\n[1, 4,  3]]), Matrix([\\n[-2, 4],\\n[ 0, 1]])))")
  expect_identical(got[["lti:test_StateSpace_series:3847"]], "StateSpace(Matrix([\\n[ 4,  1,  0,  0, 0],\\n[ 2, -3,  0,  0, 0],\\n[ 2,  0, -3,  4, 2],\\n[-6,  9, -1, -3, 0],\\n[-4,  9,  2,  5, 3]]), Matrix([\\n[  5,  2],\\n[ -3, -3],\\n[  7, -2],\\n[-12, -3],\\n[ -5, -5]]), Matrix([\\n[-4, 12, 4, 2, -3],\\n[ 0,  1, 1, 4,  3]]), Matrix([\\n[-2, -8],\\n[ 1, -1]]))")
  expect_identical(got[["lti:test_StateSpace_series:3866"]], "2")
  expect_identical(got[["lti:test_StateSpace_series:3867"]], "2")
  expect_identical(got[["lti:test_StateSpace_series:3870"]], "StateSpace(Matrix([\\n[ 4,  1,  0,  0, 0],\\n[ 2, -3,  0,  0, 0],\\n[ 2,  0, -3,  4, 2],\\n[-6,  9, -1, -3, 0],\\n[-4,  9,  2,  5, 3]]), Matrix([\\n[  5,  2],\\n[ -3, -3],\\n[  7, -2],\\n[-12, -3],\\n[ -5, -5]]), Matrix([\\n[-4, 12, 4, 2, -3],\\n[ 0,  1, 1, 4,  3]]), Matrix([\\n[-2, -8],\\n[ 1, -1]]))")
  expect_identical(got[["lti:test_StateSpace_series:3871"]], "[[TransferFunction(-2*s^5 - 41*s^4 + 130*s^3 + 2383*s^2 + 4038*s - 42024, s^5 + 2*s^4 - 26*s^3 - 74*s^2 + 167*s + 574, s), TransferFunction(-8*s^5 - 59*s^4 + 361*s^3 + 2080*s^2 + 804*s - 34236, s^5 + 2*s^4 - 26*s^3 - 74*s^2 + 167*s + 574, s)], [TransferFunction(s^5 - 57*s^4 - 630*s^3 - 231*s^2 + 7719*s + 11778, s^5 + 2*s^4 - 26*s^3 - 74*s^2 + 167*s + 574, s), TransferFunction(-s^5 - 34*s^4 - 356*s^3 - 53*s^2 + 5498*s + 8968, s^5 + 2*s^4 - 26*s^3 - 74*s^2 + 167*s + 574, s)]]")
  expect_identical(got[["lti:test_StateSpace_series:3877"]], "MIMOSeries(StateSpace(Matrix([\\n[4,  1],\\n[2, -3]]), Matrix([\\n[ 5,  2],\\n[-3, -3]]), Matrix([\\n[2, -4],\\n[0,  1]]), Matrix([\\n[3,  2],\\n[1, -1]])), TransferFunctionMatrix(((TransferFunction(s, s + 1, s), TransferFunction(1, s, s)), (TransferFunction(1, s + 1, s), TransferFunction(s, s + 2, s)))))")
})

test_that("lti test_StateSpace_parallel matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("a, x, b, c, s, g, d,    p, k, tau, zeta, wn, T, z"); a <- .py_getitem(.tmpu_, 0L); x <- .py_getitem(.tmpu_, 1L); b <- .py_getitem(.tmpu_, 2L); c <- .py_getitem(.tmpu_, 3L); s <- .py_getitem(.tmpu_, 4L); g <- .py_getitem(.tmpu_, 5L); d <- .py_getitem(.tmpu_, 6L); p <- .py_getitem(.tmpu_, 7L); k <- .py_getitem(.tmpu_, 8L); tau <- .py_getitem(.tmpu_, 9L); zeta <- .py_getitem(.tmpu_, 10L); wn <- .py_getitem(.tmpu_, 11L); T <- .py_getitem(.tmpu_, 12L); z <- .py_getitem(.tmpu_, 13L)
})
.py_try({
.tmpu_ <- .py_symbols("a0:4, b0:5, c0:4, d0:4"); a0 <- .py_getitem(.tmpu_, 0L); a1 <- .py_getitem(.tmpu_, 1L); a2 <- .py_getitem(.tmpu_, 2L); a3 <- .py_getitem(.tmpu_, 3L); b0 <- .py_getitem(.tmpu_, 4L); b1 <- .py_getitem(.tmpu_, 5L); b2 <- .py_getitem(.tmpu_, 6L); b3 <- .py_getitem(.tmpu_, 7L); b4 <- .py_getitem(.tmpu_, 8L); c0 <- .py_getitem(.tmpu_, 9L); c1 <- .py_getitem(.tmpu_, 10L); c2 <- .py_getitem(.tmpu_, 11L); c3 <- .py_getitem(.tmpu_, 12L); d0 <- .py_getitem(.tmpu_, 13L); d1 <- .py_getitem(.tmpu_, 14L); d2 <- .py_getitem(.tmpu_, 15L); d3 <- .py_getitem(.tmpu_, 16L)
})
.py_try({
TF1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
TF2 <- TransferFunction(k, 1L, s)
})
.py_try({
TF3 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
TF4 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(tau, s), 1L), s)
})
.py_try({
TF5 <- TransferFunction(k, p, s)
})
.py_try({
TF6 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(k, s), p), s)
})
p_are_floats_equals <- function(n1, n2) {
  "\n    Helper function to compare two floats.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return(.py_method(Float(n1), "is_same", Float(n2), .py_isclose))
}
p_are_complex_equals <- function(n1, n2) {
  "\n    Helper function to compare two complex numbers.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return((.py_truthy(p_are_floats_equals(re(n1), re(n2))) && .py_truthy(p_are_floats_equals(im(n1), im(n2)))))
}
.py_try({
a1 <- .py_Matrix(list(list(0L, 1L), list(1L, 0L)))
})
.py_try({
b1 <- .py_Matrix(list(list(0L), list(1L)))
})
.py_try({
c1 <- .py_Matrix(list(list(0L, 1L)))
})
.py_try({
d1 <- .py_Matrix(list(list(0L)))
})
.py_try({
a2 <- .py_Matrix(list(list(1L, 0L), list(0L, 1L)))
})
.py_try({
b2 <- .py_Matrix(list(list(1L), list(0L)))
})
.py_try({
c2 <- .py_Matrix(list(list(1L, 0L)))
})
.py_try({
d2 <- .py_Matrix(list(list(1L)))
})
.py_try({
ss1 <- StateSpace(a1, b1, c1, d1)
})
.py_try({
ss2 <- StateSpace(a2, b2, c2, d2)
})
.py_try({
p1 <- Parallel(ss1, ss2)
})
.py_try({
.py_rec("lti:test_StateSpace_parallel:3902", function() p1)
})
.py_try({
.py_rec("lti:test_StateSpace_parallel:3904", function() .py_method(p1, "doit"))
})
.py_try({
.py_rec("lti:test_StateSpace_parallel:3916", function() .py_method(p1, "rewrite", TransferFunction))
})
.py_try({
.py_rec("lti:test_StateSpace_parallel:3917", function() .py_method(p1, "rewrite", DiscreteTransferFunction))
})
.py_try({
tf1 <- TransferFunction(s, .py_add(s, 1L), s)
})
.py_try({
p2 <- Parallel(ss1, tf1)
})
.py_try({
.py_rec("lti:test_StateSpace_parallel:3922", function() p2)
})
.py_try({
.py_rec("lti:test_StateSpace_parallel:3927", function() .py_method(p2, "doit"))
})
.py_try({
.py_rec("lti:test_StateSpace_parallel:3938", function() .py_method(p2, "rewrite", TransferFunction))
})
.py_try({
a3 <- .py_Matrix(list(list(4L, 1L), list(2L, -3L)))
})
.py_try({
b3 <- .py_Matrix(list(list(5L, 2L), list(-3L, -3L)))
})
.py_try({
c3 <- .py_Matrix(list(list(2L, -4L), list(0L, 1L)))
})
.py_try({
d3 <- .py_Matrix(list(list(3L, 2L), list(1L, -1L)))
})
.py_try({
a4 <- .py_Matrix(list(list(-3L, 4L, 2L), list(-1L, -3L, 0L), list(2L, 5L, 3L)))
})
.py_try({
b4 <- .py_Matrix(list(list(1L, 4L), list(-3L, -3L), list(-2L, 1L)))
})
.py_try({
c4 <- .py_Matrix(list(list(4L, 2L, -3L), list(1L, 4L, 3L)))
})
.py_try({
d4 <- .py_Matrix(list(list(-2L, 4L), list(0L, 1L)))
})
.py_try({
ss3 <- StateSpace(a3, b3, c3, d3)
})
.py_try({
ss4 <- StateSpace(a4, b4, c4, d4)
})
.py_try({
p3 <- MIMOParallel(ss3, ss4)
})
.py_try({
.py_rec("lti:test_StateSpace_parallel:3952", function() p3)
})
.py_try({
.py_rec("lti:test_StateSpace_parallel:3971", function() .py_method(p3, "doit"))
})
.py_try({
tf2 <- TransferFunction(1L, s, s)
})
.py_try({
tf3 <- TransferFunction(1L, .py_add(s, 1L), s)
})
.py_try({
tf4 <- TransferFunction(s, .py_add(s, 2L), s)
})
.py_try({
tfm <- TransferFunctionMatrix(list(list(tf1, tf2), list(tf3, tf4)))
})
.py_try({
p4 <- MIMOParallel(tfm, ss3)
})
.py_try({
.py_rec("lti:test_StateSpace_parallel:3996", function() p4)
}))---")
  expect_identical(got[["lti:test_StateSpace_parallel:3902"]], "Parallel(StateSpace(Matrix([\\n[0, 1],\\n[1, 0]]), Matrix([\\n[0],\\n[1]]), Matrix([[0, 1]]), Matrix([[0]])), StateSpace(Matrix([\\n[1, 0],\\n[0, 1]]), Matrix([\\n[1],\\n[0]]), Matrix([[1, 0]]), Matrix([[1]])))")
  expect_identical(got[["lti:test_StateSpace_parallel:3904"]], "StateSpace(Matrix([\\n[0, 1, 0, 0],\\n[1, 0, 0, 0],\\n[0, 0, 1, 0],\\n[0, 0, 0, 1]]), Matrix([\\n[0],\\n[1],\\n[1],\\n[0]]), Matrix([[0, 1, 1, 0]]), Matrix([[1]]))")
  expect_identical(got[["lti:test_StateSpace_parallel:3916"]], "TransferFunction(s*(s + 2), s^2 - 1, s)")
  expect_identical(got[["lti:test_StateSpace_parallel:3917"]], "ERR TypeError: \\n            The continuous state space model cannot be rewritten as a\\n            discrete-time transfer function model.\\n            ")
  expect_identical(got[["lti:test_StateSpace_parallel:3922"]], "Parallel(StateSpace(Matrix([\\n[0, 1],\\n[1, 0]]), Matrix([\\n[0],\\n[1]]), Matrix([[0, 1]]), Matrix([[0]])), TransferFunction(s, s + 1, s))")
  expect_identical(got[["lti:test_StateSpace_parallel:3927"]], "StateSpace(Matrix([\\n[0, 1,  0],\\n[1, 0,  0],\\n[0, 0, -1]]), Matrix([\\n[0],\\n[1],\\n[1]]), Matrix([[0, 1, -1]]), Matrix([[1]]))")
  expect_identical(got[["lti:test_StateSpace_parallel:3938"]], "TransferFunction(s^2, s^2 - 1, s)")
  expect_identical(got[["lti:test_StateSpace_parallel:3952"]], "MIMOParallel(StateSpace(Matrix([\\n[4,  1],\\n[2, -3]]), Matrix([\\n[ 5,  2],\\n[-3, -3]]), Matrix([\\n[2, -4],\\n[0,  1]]), Matrix([\\n[3,  2],\\n[1, -1]])), StateSpace(Matrix([\\n[-3,  4, 2],\\n[-1, -3, 0],\\n[ 2,  5, 3]]), Matrix([\\n[ 1,  4],\\n[-3, -3],\\n[-2,  1]]), Matrix([\\n[4, 2, -3],\\n[1, 4,  3]]), Matrix([\\n[-2, 4],\\n[ 0, 1]])))")
  expect_identical(got[["lti:test_StateSpace_parallel:3971"]], "StateSpace(Matrix([\\n[4,  1,  0,  0, 0],\\n[2, -3,  0,  0, 0],\\n[0,  0, -3,  4, 2],\\n[0,  0, -1, -3, 0],\\n[0,  0,  2,  5, 3]]), Matrix([\\n[ 5,  2],\\n[-3, -3],\\n[ 1,  4],\\n[-3, -3],\\n[-2,  1]]), Matrix([\\n[2, -4, 4, 2, -3],\\n[0,  1, 1, 4,  3]]), Matrix([\\n[1, 6],\\n[1, 0]]))")
  expect_identical(got[["lti:test_StateSpace_parallel:3996"]], "MIMOParallel(TransferFunctionMatrix(((TransferFunction(s, s + 1, s), TransferFunction(1, s, s)), (TransferFunction(1, s + 1, s), TransferFunction(s, s + 2, s)))), StateSpace(Matrix([\\n[4,  1],\\n[2, -3]]), Matrix([\\n[ 5,  2],\\n[-3, -3]]), Matrix([\\n[2, -4],\\n[0,  1]]), Matrix([\\n[3,  2],\\n[1, -1]])))")
})

test_that("lti test_StateSpace_feedback matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("a, x, b, c, s, g, d,    p, k, tau, zeta, wn, T, z"); a <- .py_getitem(.tmpu_, 0L); x <- .py_getitem(.tmpu_, 1L); b <- .py_getitem(.tmpu_, 2L); c <- .py_getitem(.tmpu_, 3L); s <- .py_getitem(.tmpu_, 4L); g <- .py_getitem(.tmpu_, 5L); d <- .py_getitem(.tmpu_, 6L); p <- .py_getitem(.tmpu_, 7L); k <- .py_getitem(.tmpu_, 8L); tau <- .py_getitem(.tmpu_, 9L); zeta <- .py_getitem(.tmpu_, 10L); wn <- .py_getitem(.tmpu_, 11L); T <- .py_getitem(.tmpu_, 12L); z <- .py_getitem(.tmpu_, 13L)
})
.py_try({
.tmpu_ <- .py_symbols("a0:4, b0:5, c0:4, d0:4"); a0 <- .py_getitem(.tmpu_, 0L); a1 <- .py_getitem(.tmpu_, 1L); a2 <- .py_getitem(.tmpu_, 2L); a3 <- .py_getitem(.tmpu_, 3L); b0 <- .py_getitem(.tmpu_, 4L); b1 <- .py_getitem(.tmpu_, 5L); b2 <- .py_getitem(.tmpu_, 6L); b3 <- .py_getitem(.tmpu_, 7L); b4 <- .py_getitem(.tmpu_, 8L); c0 <- .py_getitem(.tmpu_, 9L); c1 <- .py_getitem(.tmpu_, 10L); c2 <- .py_getitem(.tmpu_, 11L); c3 <- .py_getitem(.tmpu_, 12L); d0 <- .py_getitem(.tmpu_, 13L); d1 <- .py_getitem(.tmpu_, 14L); d2 <- .py_getitem(.tmpu_, 15L); d3 <- .py_getitem(.tmpu_, 16L)
})
.py_try({
TF1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
TF2 <- TransferFunction(k, 1L, s)
})
.py_try({
TF3 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
TF4 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(tau, s), 1L), s)
})
.py_try({
TF5 <- TransferFunction(k, p, s)
})
.py_try({
TF6 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(k, s), p), s)
})
p_are_floats_equals <- function(n1, n2) {
  "\n    Helper function to compare two floats.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return(.py_method(Float(n1), "is_same", Float(n2), .py_isclose))
}
p_are_complex_equals <- function(n1, n2) {
  "\n    Helper function to compare two complex numbers.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return((.py_truthy(p_are_floats_equals(re(n1), re(n2))) && .py_truthy(p_are_floats_equals(im(n1), im(n2)))))
}
.py_try({
a1 <- .py_Matrix(list(list(0L, 1L), list(1L, 0L)))
})
.py_try({
b1 <- .py_Matrix(list(list(0L), list(1L)))
})
.py_try({
c1 <- .py_Matrix(list(list(0L, 1L)))
})
.py_try({
d1 <- .py_Matrix(list(list(0L)))
})
.py_try({
a2 <- .py_Matrix(list(list(1L, 0L), list(0L, 1L)))
})
.py_try({
b2 <- .py_Matrix(list(list(1L), list(0L)))
})
.py_try({
c2 <- .py_Matrix(list(list(1L, 0L)))
})
.py_try({
d2 <- .py_Matrix(list(list(1L)))
})
.py_try({
ss1 <- StateSpace(a1, b1, c1, d1)
})
.py_try({
ss2 <- StateSpace(a2, b2, c2, d2)
})
.py_try({
fd1 <- Feedback(ss1, ss2)
})
.py_try({
.py_rec("lti:test_StateSpace_feedback:4024", function() fd1)
})
.py_try({
.py_rec("lti:test_StateSpace_feedback:4026", function() .py_method(fd1, "doit"))
})
.py_try({
.py_rec("lti:test_StateSpace_feedback:4037", function() .py_method(fd1, "rewrite", TransferFunction))
})
.py_try({
.py_rec("lti:test_StateSpace_feedback:4038", function() .py_method(fd1, "rewrite", DiscreteTransferFunction))
})
.py_try({
fd2 <- Feedback(ss1, ss2, 1L)
})
.py_try({
.py_rec("lti:test_StateSpace_feedback:4042", function() .py_method(fd2, "doit"))
})
.py_try({
.py_rec("lti:test_StateSpace_feedback:4053", function() .py_method(fd2, "rewrite", TransferFunction))
})
.py_try({
tf1 <- TransferFunction(s, .py_add(s, 1L), s)
})
.py_try({
fd3 <- Feedback(ss1, tf1)
})
.py_try({
.py_rec("lti:test_StateSpace_feedback:4058", function() fd3)
})
.py_try({
.py_rec("lti:test_StateSpace_feedback:4064", function() .py_method(fd3, "doit"))
})
.py_try({
a3 <- .py_Matrix(list(list(4L, 1L), list(2L, -3L)))
})
.py_try({
b3 <- .py_Matrix(list(list(5L, 2L), list(-3L, -3L)))
})
.py_try({
c3 <- .py_Matrix(list(list(2L, -4L), list(0L, 1L)))
})
.py_try({
d3 <- .py_Matrix(list(list(3L, 2L), list(1L, -1L)))
})
.py_try({
a4 <- .py_Matrix(list(list(-3L, 4L, 2L), list(-1L, -3L, 0L), list(2L, 5L, 3L)))
})
.py_try({
b4 <- .py_Matrix(list(list(1L, 4L), list(-3L, -3L), list(-2L, 1L)))
})
.py_try({
c4 <- .py_Matrix(list(list(4L, 2L, -3L), list(1L, 4L, 3L)))
})
.py_try({
d4 <- .py_Matrix(list(list(-2L, 4L), list(0L, 1L)))
})
.py_try({
ss3 <- StateSpace(a3, b3, c3, d3)
})
.py_try({
ss4 <- StateSpace(a4, b4, c4, d4)
})
.py_try({
fd4 <- MIMOFeedback(ss3, ss4)
})
.py_try({
.py_rec("lti:test_StateSpace_feedback:4088", function() fd4)
})
.py_try({
.py_rec("lti:test_StateSpace_feedback:4107", function() .py_method(fd4, "doit"))
})
.py_try({
fd5 <- MIMOFeedback(ss3, ss4, 1L)
})
.py_try({
.py_rec("lti:test_StateSpace_feedback:4125", function() .py_method(fd5, "doit"))
}))---")
  expect_identical(got[["lti:test_StateSpace_feedback:4024"]], "Feedback(StateSpace(Matrix([\\n[0, 1],\\n[1, 0]]), Matrix([\\n[0],\\n[1]]), Matrix([[0, 1]]), Matrix([[0]])), StateSpace(Matrix([\\n[1, 0],\\n[0, 1]]), Matrix([\\n[1],\\n[0]]), Matrix([[1, 0]]), Matrix([[1]])), -1)")
  expect_identical(got[["lti:test_StateSpace_feedback:4026"]], "StateSpace(Matrix([\\n[0,  1,  0, 0],\\n[1, -1, -1, 0],\\n[0,  1,  1, 0],\\n[0,  0,  0, 1]]), Matrix([\\n[0],\\n[1],\\n[0],\\n[0]]), Matrix([[0, 1, 0, 0]]), Matrix([[0]]))")
  expect_identical(got[["lti:test_StateSpace_feedback:4037"]], "TransferFunction(s*(s - 1), s^3 - s + 1, s)")
  expect_identical(got[["lti:test_StateSpace_feedback:4038"]], "ERR TypeError: \\n            The continuous state space model cannot be rewritten as a\\n            discrete-time transfer function model.\\n            ")
  expect_identical(got[["lti:test_StateSpace_feedback:4042"]], "StateSpace(Matrix([\\n[0, 1, 0, 0],\\n[1, 1, 1, 0],\\n[0, 1, 1, 0],\\n[0, 0, 0, 1]]), Matrix([\\n[0],\\n[1],\\n[0],\\n[0]]), Matrix([[0, 1, 0, 0]]), Matrix([[0]]))")
  expect_identical(got[["lti:test_StateSpace_feedback:4053"]], "TransferFunction(s*(s - 1), s^3 - 2*s^2 - s + 1, s)")
  expect_identical(got[["lti:test_StateSpace_feedback:4058"]], "Feedback(StateSpace(Matrix([\\n[0, 1],\\n[1, 0]]), Matrix([\\n[0],\\n[1]]), Matrix([[0, 1]]), Matrix([[0]])), TransferFunction(s, s + 1, s), -1)")
  expect_identical(got[["lti:test_StateSpace_feedback:4064"]], "StateSpace(Matrix([\\n[0,  1,  0],\\n[1, -1,  1],\\n[0,  1, -1]]), Matrix([\\n[0],\\n[1],\\n[0]]), Matrix([[0, 1, 0]]), Matrix([[0]]))")
  expect_identical(got[["lti:test_StateSpace_feedback:4088"]], "MIMOFeedback(StateSpace(Matrix([\\n[4,  1],\\n[2, -3]]), Matrix([\\n[ 5,  2],\\n[-3, -3]]), Matrix([\\n[2, -4],\\n[0,  1]]), Matrix([\\n[3,  2],\\n[1, -1]])), StateSpace(Matrix([\\n[-3,  4, 2],\\n[-1, -3, 0],\\n[ 2,  5, 3]]), Matrix([\\n[ 1,  4],\\n[-3, -3],\\n[-2,  1]]), Matrix([\\n[4, 2, -3],\\n[1, 4,  3]]), Matrix([\\n[-2, 4],\\n[ 0, 1]])), -1)")
  expect_identical(got[["lti:test_StateSpace_feedback:4107"]], "StateSpace(Matrix([\\n[   3,  -3/4, -15/4, -37/2, -15],\\n[ 7/2, -39/8,   9/8,  39/4,   9],\\n[   3, -41/4, -45/4, -51/2, -19],\\n[-9/2, 129/8,  73/8, 171/4,  36],\\n[-3/2,  47/8,  31/8,  85/4,  18]]), Matrix([\\n[-1/4,  19/4],\\n[ 3/8, -21/8],\\n[ 1/4,  29/4],\\n[ 3/8, -93/8],\\n[ 5/8, -35/8]]), Matrix([\\n[  1, -15/4,  -7/4, -21/2, -9],\\n[1/2, -13/8, -13/8, -19/4, -3]]), Matrix([\\n[-1/4, 11/4],\\n[ 1/8,  9/8]]))")
  expect_identical(got[["lti:test_StateSpace_feedback:4125"]], "StateSpace(Matrix([\\n[  4/7,    62/7,    1,  -8, -69/7],\\n[ 32/7, -135/14, -3/2,   3,  36/7],\\n[-10/7,    41/7,   -4, -12, -97/7],\\n[ 12/7, -111/14, -5/2,  18, 171/7],\\n[  2/7,  -29/14, -1/2,  10,  81/7]]), Matrix([\\n[   6/7, -17/7],\\n[ -9/14, 15/14],\\n[   6/7, -31/7],\\n[-27/14, 87/14],\\n[-15/14, 25/14]]), Matrix([\\n[-2/7,  11/7,    1, -4, -39/7],\\n[-2/7, 15/14, -1/2, -3, -18/7]]), Matrix([\\n[ 4/7,   -9/7],\\n[1/14, -11/14]]))")
})

test_that("lti test_StateSpace_stability matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("a, x, b, c, s, g, d,    p, k, tau, zeta, wn, T, z"); a <- .py_getitem(.tmpu_, 0L); x <- .py_getitem(.tmpu_, 1L); b <- .py_getitem(.tmpu_, 2L); c <- .py_getitem(.tmpu_, 3L); s <- .py_getitem(.tmpu_, 4L); g <- .py_getitem(.tmpu_, 5L); d <- .py_getitem(.tmpu_, 6L); p <- .py_getitem(.tmpu_, 7L); k <- .py_getitem(.tmpu_, 8L); tau <- .py_getitem(.tmpu_, 9L); zeta <- .py_getitem(.tmpu_, 10L); wn <- .py_getitem(.tmpu_, 11L); T <- .py_getitem(.tmpu_, 12L); z <- .py_getitem(.tmpu_, 13L)
})
.py_try({
.tmpu_ <- .py_symbols("a0:4, b0:5, c0:4, d0:4"); a0 <- .py_getitem(.tmpu_, 0L); a1 <- .py_getitem(.tmpu_, 1L); a2 <- .py_getitem(.tmpu_, 2L); a3 <- .py_getitem(.tmpu_, 3L); b0 <- .py_getitem(.tmpu_, 4L); b1 <- .py_getitem(.tmpu_, 5L); b2 <- .py_getitem(.tmpu_, 6L); b3 <- .py_getitem(.tmpu_, 7L); b4 <- .py_getitem(.tmpu_, 8L); c0 <- .py_getitem(.tmpu_, 9L); c1 <- .py_getitem(.tmpu_, 10L); c2 <- .py_getitem(.tmpu_, 11L); c3 <- .py_getitem(.tmpu_, 12L); d0 <- .py_getitem(.tmpu_, 13L); d1 <- .py_getitem(.tmpu_, 14L); d2 <- .py_getitem(.tmpu_, 15L); d3 <- .py_getitem(.tmpu_, 16L)
})
.py_try({
TF1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
TF2 <- TransferFunction(k, 1L, s)
})
.py_try({
TF3 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
TF4 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(tau, s), 1L), s)
})
.py_try({
TF5 <- TransferFunction(k, p, s)
})
.py_try({
TF6 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(k, s), p), s)
})
p_are_floats_equals <- function(n1, n2) {
  "\n    Helper function to compare two floats.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return(.py_method(Float(n1), "is_same", Float(n2), .py_isclose))
}
p_are_complex_equals <- function(n1, n2) {
  "\n    Helper function to compare two complex numbers.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return((.py_truthy(p_are_floats_equals(re(n1), re(n2))) && .py_truthy(p_are_floats_equals(im(n1), im(n2)))))
}
.py_try({
k <- .py_symbols("k")
})
.py_try({
B <- .py_Matrix(list(1L, 0L, 0L))
})
.py_try({
C <- .py_Matrix(list(list(0L, 1L, 0L)))
})
.py_try({
D <- .py_Matrix(list(0L))
})
.py_try({
A1 <- .py_Matrix(list(list(0L, 1L, 0L), list(0L, 0L, 1L), list(.py_sub(k, 1L), .py_mul(-2L, k), -1L)))
})
.py_try({
ss1 <- StateSpace(A1, B, C, D)
})
.py_try({
ineq <- .py_method(ss1, "get_asymptotic_stability_conditions")
})
.py_try({
.py_rec("lti:test_StateSpace_stability:4158", function() ineq)
})
.py_try({
A2 <- .py_Matrix(list(list(1L, 0L, 0L), list(0L, -1L, k), list(0L, 0L, -1L)))
})
.py_try({
ss2 <- StateSpace(A2, B, C, D)
})
.py_try({
ineq <- .py_method(ss2, "get_asymptotic_stability_conditions")
})
.py_try({
.py_rec("lti:test_StateSpace_stability:4163", function() ineq)
}))---")
  expect_identical(got[["lti:test_StateSpace_stability:4158"]], "[True, 3*k - 1 > 0, 1 - k > 0]")
  expect_identical(got[["lti:test_StateSpace_stability:4163"]], "[False]")
})

test_that("lti test_DiscreteStateSpace_construction matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("a, x, b, c, s, g, d,    p, k, tau, zeta, wn, T, z"); a <- .py_getitem(.tmpu_, 0L); x <- .py_getitem(.tmpu_, 1L); b <- .py_getitem(.tmpu_, 2L); c <- .py_getitem(.tmpu_, 3L); s <- .py_getitem(.tmpu_, 4L); g <- .py_getitem(.tmpu_, 5L); d <- .py_getitem(.tmpu_, 6L); p <- .py_getitem(.tmpu_, 7L); k <- .py_getitem(.tmpu_, 8L); tau <- .py_getitem(.tmpu_, 9L); zeta <- .py_getitem(.tmpu_, 10L); wn <- .py_getitem(.tmpu_, 11L); T <- .py_getitem(.tmpu_, 12L); z <- .py_getitem(.tmpu_, 13L)
})
.py_try({
.tmpu_ <- .py_symbols("a0:4, b0:5, c0:4, d0:4"); a0 <- .py_getitem(.tmpu_, 0L); a1 <- .py_getitem(.tmpu_, 1L); a2 <- .py_getitem(.tmpu_, 2L); a3 <- .py_getitem(.tmpu_, 3L); b0 <- .py_getitem(.tmpu_, 4L); b1 <- .py_getitem(.tmpu_, 5L); b2 <- .py_getitem(.tmpu_, 6L); b3 <- .py_getitem(.tmpu_, 7L); b4 <- .py_getitem(.tmpu_, 8L); c0 <- .py_getitem(.tmpu_, 9L); c1 <- .py_getitem(.tmpu_, 10L); c2 <- .py_getitem(.tmpu_, 11L); c3 <- .py_getitem(.tmpu_, 12L); d0 <- .py_getitem(.tmpu_, 13L); d1 <- .py_getitem(.tmpu_, 14L); d2 <- .py_getitem(.tmpu_, 15L); d3 <- .py_getitem(.tmpu_, 16L)
})
.py_try({
TF1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
TF2 <- TransferFunction(k, 1L, s)
})
.py_try({
TF3 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
TF4 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(tau, s), 1L), s)
})
.py_try({
TF5 <- TransferFunction(k, p, s)
})
.py_try({
TF6 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(k, s), p), s)
})
p_are_floats_equals <- function(n1, n2) {
  "\n    Helper function to compare two floats.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return(.py_method(Float(n1), "is_same", Float(n2), .py_isclose))
}
p_are_complex_equals <- function(n1, n2) {
  "\n    Helper function to compare two complex numbers.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return((.py_truthy(p_are_floats_equals(re(n1), re(n2))) && .py_truthy(p_are_floats_equals(im(n1), im(n2)))))
}
.py_try({
A1 <- .py_Matrix(list(list(0L, 1L), list(1L, 0L)))
})
.py_try({
B1 <- .py_Matrix(list(1L, 0L))
})
.py_try({
C1 <- .py_Matrix(list(list(0L, 1L)))
})
.py_try({
D1 <- .py_Matrix(list(0L))
})
.py_try({
ss1 <- DiscreteStateSpace(A1, B1, C1, D1)
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_construction:4173", function() .py_attr(ss1, "state_matrix"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_construction:4174", function() .py_attr(ss1, "input_matrix"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_construction:4175", function() .py_attr(ss1, "output_matrix"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_construction:4176", function() .py_attr(ss1, "feedforward_matrix"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_construction:4177", function() .py_attr(ss1, "sampling_time"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_construction:4178", function() .py_attr(ss1, "args"))
})
.py_try({
ss2 <- DiscreteStateSpace(.py_Matrix(list(a0)), .py_Matrix(list(a1)), .py_Matrix(list(a2)), .py_Matrix(list(a3)), .py_flt(0.1))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_construction:4185", function() .py_attr(ss2, "state_matrix"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_construction:4186", function() .py_attr(ss2, "input_matrix"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_construction:4187", function() .py_attr(ss2, "output_matrix"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_construction:4188", function() .py_attr(ss2, "feedforward_matrix"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_construction:4189", function() .py_attr(ss2, "sampling_time"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_construction:4190", function() .py_isinstance(.py_attr(ss2, "sampling_time"), "Number"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_construction:4191", function() .py_attr(ss2, "args"))
})
.py_try({
ss3 <- DiscreteStateSpace(.py_Matrix(list(list(.py_flt(-1.5), -2L), list(1L, 0L))), .py_Matrix(list(list(.py_flt(0.5), 0L), list(0L, 1L))), .py_Matrix(list(list(0L, 1L), list(0L, 2L))), .py_Matrix(list(list(2L, 2L), list(1L, 1L))), T)
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_construction:4200", function() .py_attr(ss3, "state_matrix"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_construction:4201", function() .py_attr(ss3, "input_matrix"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_construction:4202", function() .py_attr(ss3, "output_matrix"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_construction:4203", function() .py_attr(ss3, "feedforward_matrix"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_construction:4204", function() .py_attr(ss3, "sampling_time"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_construction:4205", function() .py_attr(ss3, "args"))
})
.py_try({
A4 <- .py_Matrix(list(list(a0, a1), list(a2, a3)))
})
.py_try({
B4 <- .py_Matrix(list(list(b0, b1), list(b2, b3)))
})
.py_try({
C4 <- .py_Matrix(list(list(c0, c1), list(c2, c3)))
})
.py_try({
D4 <- .py_Matrix(list(list(d0, d1), list(d2, d3)))
})
.py_try({
ss4 <- DiscreteStateSpace(A4, B4, C4, D4, 12L)
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_construction:4221", function() .py_attr(ss4, "state_matrix"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_construction:4222", function() .py_attr(ss4, "input_matrix"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_construction:4223", function() .py_attr(ss4, "output_matrix"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_construction:4224", function() .py_attr(ss4, "feedforward_matrix"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_construction:4225", function() .py_attr(ss4, "sampling_time"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_construction:4226", function() .py_attr(ss4, "args"))
})
.py_try({
ss5 <- DiscreteStateSpace()
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_construction:4237", function() .py_attr(ss5, "args"))
})
.py_try({
A6 <- .py_Matrix(list(list(0L, 1L), list(1L, 0L)))
})
.py_try({
B6 <- .py_Matrix(list(1L, 1L))
})
.py_try({
ss6 <- DiscreteStateSpace(A6, B6, sampling_time = 2L)
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_construction:4244", function() .py_attr(ss6, "state_matrix"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_construction:4245", function() .py_attr(ss6, "input_matrix"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_construction:4246", function() .py_attr(ss6, "output_matrix"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_construction:4247", function() .py_attr(ss6, "feedforward_matrix"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_construction:4248", function() .py_attr(ss6, "sampling_time"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_construction:4249", function() .py_attr(ss6, "args"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_construction:4259", function() .py_attr(ss1, "is_SISO"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_construction:4260", function() .py_attr(ss2, "is_SISO"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_construction:4261", function() .py_attr(ss3, "is_SISO"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_construction:4262", function() .py_attr(ss4, "is_SISO"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_construction:4263", function() .py_attr(ss5, "is_SISO"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_construction:4264", function() .py_attr(ss6, "is_SISO"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_construction:4267", function() StateSpace(.py_Matrix(list(s, .py_pow(.py_add(s, 1L), 2L))), .py_Matrix(list(.py_add(s, 1L))), .py_Matrix(list(.py_sub(.py_pow(s, 2L), 1L))), .py_Matrix(list(.py_mul(2L, s)))))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_construction:4269", function() StateSpace(.py_Matrix(list(s)), .py_Matrix(list(.py_add(s, 1L), .py_add(.py_pow(s, 3L), 1L))), .py_Matrix(list(.py_sub(.py_pow(s, 2L), 1L))), .py_Matrix(list(.py_mul(2L, s)))))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_construction:4271", function() StateSpace(.py_Matrix(list(s)), .py_Matrix(list(.py_add(s, 1L))), .py_Matrix(list(list(.py_sub(.py_pow(s, 2L), 1L)), list(.py_add(.py_add(.py_pow(s, 2L), .py_mul(2L, s)), 1L)))), .py_Matrix(list(.py_mul(2L, s)))))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_construction:4273", function() StateSpace(.py_Matrix(list(list(.py_neg(s), .py_neg(s)), list(s, 0L))), .py_Matrix(list(list(.py_div(s, 2L), 0L), list(0L, s))), .py_Matrix(list(list(0L, s))), .py_Matrix(list(list(.py_mul(2L, s), .py_mul(2L, s)), list(s, s)))))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_construction:4279", function() StateSpace(.py_pow(s, 2L), .py_add(s, 1L), .py_mul(2L, s), 1L))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_construction:4280", function() StateSpace(.py_Matrix(list(2L, .py_flt(0.5))), .py_Matrix(list(-1L)), .py_Matrix(list(1L)), 0L))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_construction:4283", function() DiscreteStateSpace(A1, B1, C1, D1, 0L))
}))---")
  expect_identical(got[["lti:test_DiscreteStateSpace_construction:4173"]], "Matrix([[0, 1], [1, 0]])")
  expect_identical(got[["lti:test_DiscreteStateSpace_construction:4174"]], "Matrix([[1], [0]])")
  expect_identical(got[["lti:test_DiscreteStateSpace_construction:4175"]], "Matrix([[0, 1]])")
  expect_identical(got[["lti:test_DiscreteStateSpace_construction:4176"]], "Matrix([[0]])")
  expect_identical(got[["lti:test_DiscreteStateSpace_construction:4177"]], "1")
  expect_identical(got[["lti:test_DiscreteStateSpace_construction:4178"]], "(Matrix([[0, 1], [1, 0]]), Matrix([[1], [0]]), Matrix([[0, 1]]), Matrix([[0]]), 1)")
  expect_identical(got[["lti:test_DiscreteStateSpace_construction:4185"]], "Matrix([[a0]])")
  expect_identical(got[["lti:test_DiscreteStateSpace_construction:4186"]], "Matrix([[a1]])")
  expect_identical(got[["lti:test_DiscreteStateSpace_construction:4187"]], "Matrix([[a2]])")
  expect_identical(got[["lti:test_DiscreteStateSpace_construction:4188"]], "Matrix([[a3]])")
  expect_identical(got[["lti:test_DiscreteStateSpace_construction:4189"]], "0.100000000000000")
  expect_identical(got[["lti:test_DiscreteStateSpace_construction:4190"]], "True")
  expect_identical(got[["lti:test_DiscreteStateSpace_construction:4191"]], "(Matrix([[a0]]), Matrix([[a1]]), Matrix([[a2]]), Matrix([[a3]]), 0.100000000000000)")
  expect_identical(got[["lti:test_DiscreteStateSpace_construction:4200"]], "Matrix([[-1.50000000000000, -2], [1, 0]])")
  expect_identical(got[["lti:test_DiscreteStateSpace_construction:4201"]], "Matrix([[0.500000000000000, 0], [0, 1]])")
  expect_identical(got[["lti:test_DiscreteStateSpace_construction:4202"]], "Matrix([[0, 1], [0, 2]])")
  expect_identical(got[["lti:test_DiscreteStateSpace_construction:4203"]], "Matrix([[2, 2], [1, 1]])")
  expect_identical(got[["lti:test_DiscreteStateSpace_construction:4204"]], "T")
  expect_identical(got[["lti:test_DiscreteStateSpace_construction:4205"]], "(Matrix([[-1.50000000000000, -2], [1, 0]]), Matrix([[0.500000000000000, 0], [0, 1]]), Matrix([[0, 1], [0, 2]]), Matrix([[2, 2], [1, 1]]), T)")
  expect_identical(got[["lti:test_DiscreteStateSpace_construction:4221"]], "Matrix([[a0, a1], [a2, a3]])")
  expect_identical(got[["lti:test_DiscreteStateSpace_construction:4222"]], "Matrix([[b0, b1], [b2, b3]])")
  expect_identical(got[["lti:test_DiscreteStateSpace_construction:4223"]], "Matrix([[c0, c1], [c2, c3]])")
  expect_identical(got[["lti:test_DiscreteStateSpace_construction:4224"]], "Matrix([[d0, d1], [d2, d3]])")
  expect_identical(got[["lti:test_DiscreteStateSpace_construction:4225"]], "12")
  expect_identical(got[["lti:test_DiscreteStateSpace_construction:4226"]], "(Matrix([[a0, a1], [a2, a3]]), Matrix([[b0, b1], [b2, b3]]), Matrix([[c0, c1], [c2, c3]]), Matrix([[d0, d1], [d2, d3]]), 12)")
  expect_identical(got[["lti:test_DiscreteStateSpace_construction:4237"]], "(Matrix([[0]]), Matrix([[0]]), Matrix([[0]]), Matrix([[0]]), 1)")
  expect_identical(got[["lti:test_DiscreteStateSpace_construction:4244"]], "Matrix([[0, 1], [1, 0]])")
  expect_identical(got[["lti:test_DiscreteStateSpace_construction:4245"]], "Matrix([[1], [1]])")
  expect_identical(got[["lti:test_DiscreteStateSpace_construction:4246"]], "Matrix([[0, 0]])")
  expect_identical(got[["lti:test_DiscreteStateSpace_construction:4247"]], "Matrix([[0]])")
  expect_identical(got[["lti:test_DiscreteStateSpace_construction:4248"]], "2")
  expect_identical(got[["lti:test_DiscreteStateSpace_construction:4249"]], "(Matrix([[0, 1], [1, 0]]), Matrix([[1], [1]]), Matrix([[0, 0]]), Matrix([[0]]), 2)")
  expect_identical(got[["lti:test_DiscreteStateSpace_construction:4259"]], "True")
  expect_identical(got[["lti:test_DiscreteStateSpace_construction:4260"]], "True")
  expect_identical(got[["lti:test_DiscreteStateSpace_construction:4261"]], "False")
  expect_identical(got[["lti:test_DiscreteStateSpace_construction:4262"]], "False")
  expect_identical(got[["lti:test_DiscreteStateSpace_construction:4263"]], "True")
  expect_identical(got[["lti:test_DiscreteStateSpace_construction:4264"]], "True")
  expect_identical(got[["lti:test_DiscreteStateSpace_construction:4267"]], "ERR ShapeError: Matrix A must be a square matrix.")
  expect_identical(got[["lti:test_DiscreteStateSpace_construction:4269"]], "ERR ShapeError: Matrices A and B must have the same number of rows.")
  expect_identical(got[["lti:test_DiscreteStateSpace_construction:4271"]], "ERR ShapeError: Matrices C and D must have the same number of rows.")
  expect_identical(got[["lti:test_DiscreteStateSpace_construction:4273"]], "ERR ShapeError: Matrices C and D must have the same number of rows.")
  expect_identical(got[["lti:test_DiscreteStateSpace_construction:4279"]], "ERR TypeError: A, B, C and D inputs must all be sympy Matrices.")
  expect_identical(got[["lti:test_DiscreteStateSpace_construction:4280"]], "ERR TypeError: A, B, C and D inputs must all be sympy Matrices.")
  expect_identical(got[["lti:test_DiscreteStateSpace_construction:4283"]], "ERR ValueError: \\nThe sampling time cannot be zero. If you want to create a continuous\\nstate space, use the StateSpace class instead.")
})

test_that("lti test_DiscreteStateSpace_add_mul matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("a, x, b, c, s, g, d,    p, k, tau, zeta, wn, T, z"); a <- .py_getitem(.tmpu_, 0L); x <- .py_getitem(.tmpu_, 1L); b <- .py_getitem(.tmpu_, 2L); c <- .py_getitem(.tmpu_, 3L); s <- .py_getitem(.tmpu_, 4L); g <- .py_getitem(.tmpu_, 5L); d <- .py_getitem(.tmpu_, 6L); p <- .py_getitem(.tmpu_, 7L); k <- .py_getitem(.tmpu_, 8L); tau <- .py_getitem(.tmpu_, 9L); zeta <- .py_getitem(.tmpu_, 10L); wn <- .py_getitem(.tmpu_, 11L); T <- .py_getitem(.tmpu_, 12L); z <- .py_getitem(.tmpu_, 13L)
})
.py_try({
.tmpu_ <- .py_symbols("a0:4, b0:5, c0:4, d0:4"); a0 <- .py_getitem(.tmpu_, 0L); a1 <- .py_getitem(.tmpu_, 1L); a2 <- .py_getitem(.tmpu_, 2L); a3 <- .py_getitem(.tmpu_, 3L); b0 <- .py_getitem(.tmpu_, 4L); b1 <- .py_getitem(.tmpu_, 5L); b2 <- .py_getitem(.tmpu_, 6L); b3 <- .py_getitem(.tmpu_, 7L); b4 <- .py_getitem(.tmpu_, 8L); c0 <- .py_getitem(.tmpu_, 9L); c1 <- .py_getitem(.tmpu_, 10L); c2 <- .py_getitem(.tmpu_, 11L); c3 <- .py_getitem(.tmpu_, 12L); d0 <- .py_getitem(.tmpu_, 13L); d1 <- .py_getitem(.tmpu_, 14L); d2 <- .py_getitem(.tmpu_, 15L); d3 <- .py_getitem(.tmpu_, 16L)
})
.py_try({
TF1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
TF2 <- TransferFunction(k, 1L, s)
})
.py_try({
TF3 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
TF4 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(tau, s), 1L), s)
})
.py_try({
TF5 <- TransferFunction(k, p, s)
})
.py_try({
TF6 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(k, s), p), s)
})
p_are_floats_equals <- function(n1, n2) {
  "\n    Helper function to compare two floats.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return(.py_method(Float(n1), "is_same", Float(n2), .py_isclose))
}
p_are_complex_equals <- function(n1, n2) {
  "\n    Helper function to compare two complex numbers.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return((.py_truthy(p_are_floats_equals(re(n1), re(n2))) && .py_truthy(p_are_floats_equals(im(n1), im(n2)))))
}
.py_try({
A1 <- .py_Matrix(list(list(4L, 1L), list(2L, -3L)))
})
.py_try({
B1 <- .py_Matrix(list(list(5L, 2L), list(-3L, -3L)))
})
.py_try({
C1 <- .py_Matrix(list(list(2L, -4L), list(0L, 1L)))
})
.py_try({
D1 <- .py_Matrix(list(list(3L, 2L), list(1L, -1L)))
})
.py_try({
ss1 <- DiscreteStateSpace(A1, B1, C1, D1, T)
})
.py_try({
A2 <- .py_Matrix(list(list(-3L, 4L, 2L), list(-1L, -3L, 0L), list(2L, 5L, 3L)))
})
.py_try({
B2 <- .py_Matrix(list(list(1L, 4L), list(-3L, -3L), list(-2L, 1L)))
})
.py_try({
C2 <- .py_Matrix(list(list(4L, 2L, -3L), list(1L, 4L, 3L)))
})
.py_try({
D2 <- .py_Matrix(list(list(-2L, 4L), list(0L, 1L)))
})
.py_try({
ss2 <- DiscreteStateSpace(A2, B2, C2, D2, T)
})
.py_try({
ss3 <- DiscreteStateSpace(sampling_time = T)
})
.py_try({
ss4 <- DiscreteStateSpace(.py_Matrix(list(1L)), .py_Matrix(list(2L)), .py_Matrix(list(3L)), .py_Matrix(list(4L)), T)
})
.py_try({
expected_add <- DiscreteStateSpace(.py_Matrix(list(list(4L, 1L, 0L, 0L, 0L), list(2L, -3L, 0L, 0L, 0L), list(0L, 0L, -3L, 4L, 2L), list(0L, 0L, -1L, -3L, 0L), list(0L, 0L, 2L, 5L, 3L))), .py_Matrix(list(list(5L, 2L), list(-3L, -3L), list(1L, 4L), list(-3L, -3L), list(-2L, 1L))), .py_Matrix(list(list(2L, -4L, 4L, 2L, -3L), list(0L, 1L, 1L, 4L, 3L))), .py_Matrix(list(list(1L, 6L), list(1L, 0L))), T)
})
.py_try({
expected_mul <- DiscreteStateSpace(.py_Matrix(list(list(-3L, 4L, 2L, 0L, 0L), list(-1L, -3L, 0L, 0L, 0L), list(2L, 5L, 3L, 0L, 0L), list(22L, 18L, -9L, 4L, 1L), list(-15L, -18L, 0L, 2L, -3L))), .py_Matrix(list(list(1L, 4L), list(-3L, -3L), list(-2L, 1L), list(-10L, 22L), list(6L, -15L))), .py_Matrix(list(list(14L, 14L, -3L, 2L, -4L), list(3L, -2L, -6L, 0L, 1L))), .py_Matrix(list(list(-6L, 14L), list(-2L, 3L))), T)
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_add_mul:4343", function() .py_add(ss1, ss2))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_add_mul:4344", function() .py_mul(ss1, ss2))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_add_mul:4345", function() .py_add(ss3, .py_div(1L, 2L)))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_add_mul:4347", function() .py_mul(ss4, .py_flt(1.5)))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_add_mul:4349", function() .py_mul(.py_flt(1.5), ss4))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_add_mul:4352", function() .py_add(ss1, ss3))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_add_mul:4353", function() .py_mul(ss2, ss4))
})
.py_try({
ss5 <- DiscreteStateSpace(A1, B1, C1, D1, .py_flt(0.1))
})
.py_try({
cont_ss <- StateSpace(A2, B2, C2, D2)
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_add_mul:4359", function() .py_add(ss5, ss2))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_add_mul:4360", function() .py_mul(ss5, ss2))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_add_mul:4361", function() .py_add(ss1, cont_ss))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_add_mul:4362", function() .py_mul(ss1, cont_ss))
}))---")
  expect_identical(got[["lti:test_DiscreteStateSpace_add_mul:4343"]], "DiscreteStateSpace(Matrix([\\n[4,  1,  0,  0, 0],\\n[2, -3,  0,  0, 0],\\n[0,  0, -3,  4, 2],\\n[0,  0, -1, -3, 0],\\n[0,  0,  2,  5, 3]]), Matrix([\\n[ 5,  2],\\n[-3, -3],\\n[ 1,  4],\\n[-3, -3],\\n[-2,  1]]), Matrix([\\n[2, -4, 4, 2, -3],\\n[0,  1, 1, 4,  3]]), Matrix([\\n[1, 6],\\n[1, 0]]), T)")
  expect_identical(got[["lti:test_DiscreteStateSpace_add_mul:4344"]], "DiscreteStateSpace(Matrix([\\n[ -3,   4,  2, 0,  0],\\n[ -1,  -3,  0, 0,  0],\\n[  2,   5,  3, 0,  0],\\n[ 22,  18, -9, 4,  1],\\n[-15, -18,  0, 2, -3]]), Matrix([\\n[  1,   4],\\n[ -3,  -3],\\n[ -2,   1],\\n[-10,  22],\\n[  6, -15]]), Matrix([\\n[14, 14, -3, 2, -4],\\n[ 3, -2, -6, 0,  1]]), Matrix([\\n[-6, 14],\\n[-2,  3]]), T)")
  expect_identical(got[["lti:test_DiscreteStateSpace_add_mul:4345"]], "DiscreteStateSpace(Matrix([[0]]), Matrix([[0]]), Matrix([[0]]), Matrix([[0.5]]), T)")
  expect_identical(got[["lti:test_DiscreteStateSpace_add_mul:4347"]], "DiscreteStateSpace(Matrix([[1]]), Matrix([[2]]), Matrix([[4.5]]), Matrix([[6.0]]), T)")
  expect_identical(got[["lti:test_DiscreteStateSpace_add_mul:4349"]], "DiscreteStateSpace(Matrix([[1]]), Matrix([[3.0]]), Matrix([[3]]), Matrix([[6.0]]), T)")
  expect_identical(got[["lti:test_DiscreteStateSpace_add_mul:4352"]], "ERR ShapeError: Systems with incompatible inputs and outputs cannot be added.")
  expect_identical(got[["lti:test_DiscreteStateSpace_add_mul:4353"]], "ERR ShapeError: Systems with incompatible inputs and outputs cannot be multiplied.")
  expect_identical(got[["lti:test_DiscreteStateSpace_add_mul:4359"]], "ERR TypeError: \\n                Incompatible sampling times detected.\\n                All systems should have the same sampling time.\\n                Found sampling time: 0.100000000000000 and T")
  expect_identical(got[["lti:test_DiscreteStateSpace_add_mul:4360"]], "ERR TypeError: \\n                Incompatible sampling times detected.\\n                All systems should have the same sampling time.\\n                Found sampling time: 0.100000000000000 and T")
  expect_identical(got[["lti:test_DiscreteStateSpace_add_mul:4361"]], "ERR TypeError: \\n                Incompatible systems detected.\\n                All systems should be either continuous-time or discrete-time.\\n                Found: DiscreteStateSpace(Matrix([\\n[4,  1],\\n[2, -3]]), Matrix([\\n[ 5,  2],\\n[-3, -3]]), Matrix([\\n[2, -4],\\n[0,  1]]), Matrix([\\n[3,  2],\\n[1, -1]]), T) and StateSpace(Matrix([\\n[-3,  4, 2],\\n[-1, -3, 0],\\n[ 2,  5, 3]]), Matrix([\\n[ 1,  4],\\n[-3, -3],\\n[-2,  1]]), Matrix([\\n[4, 2, -3],\\n[1, 4,  3]]), Matrix([\\n[-2, 4],\\n[ 0, 1]]))")
  expect_identical(got[["lti:test_DiscreteStateSpace_add_mul:4362"]], "ERR TypeError: \\n                Incompatible systems detected.\\n                All systems should be either continuous-time or discrete-time.\\n                Found: DiscreteStateSpace(Matrix([\\n[4,  1],\\n[2, -3]]), Matrix([\\n[ 5,  2],\\n[-3, -3]]), Matrix([\\n[2, -4],\\n[0,  1]]), Matrix([\\n[3,  2],\\n[1, -1]]), T) and StateSpace(Matrix([\\n[-3,  4, 2],\\n[-1, -3, 0],\\n[ 2,  5, 3]]), Matrix([\\n[ 1,  4],\\n[-3, -3],\\n[-2,  1]]), Matrix([\\n[4, 2, -3],\\n[1, 4,  3]]), Matrix([\\n[-2, 4],\\n[ 0, 1]]))")
})

test_that("lti test_DiscreteStateSpace_negation matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("a, x, b, c, s, g, d,    p, k, tau, zeta, wn, T, z"); a <- .py_getitem(.tmpu_, 0L); x <- .py_getitem(.tmpu_, 1L); b <- .py_getitem(.tmpu_, 2L); c <- .py_getitem(.tmpu_, 3L); s <- .py_getitem(.tmpu_, 4L); g <- .py_getitem(.tmpu_, 5L); d <- .py_getitem(.tmpu_, 6L); p <- .py_getitem(.tmpu_, 7L); k <- .py_getitem(.tmpu_, 8L); tau <- .py_getitem(.tmpu_, 9L); zeta <- .py_getitem(.tmpu_, 10L); wn <- .py_getitem(.tmpu_, 11L); T <- .py_getitem(.tmpu_, 12L); z <- .py_getitem(.tmpu_, 13L)
})
.py_try({
.tmpu_ <- .py_symbols("a0:4, b0:5, c0:4, d0:4"); a0 <- .py_getitem(.tmpu_, 0L); a1 <- .py_getitem(.tmpu_, 1L); a2 <- .py_getitem(.tmpu_, 2L); a3 <- .py_getitem(.tmpu_, 3L); b0 <- .py_getitem(.tmpu_, 4L); b1 <- .py_getitem(.tmpu_, 5L); b2 <- .py_getitem(.tmpu_, 6L); b3 <- .py_getitem(.tmpu_, 7L); b4 <- .py_getitem(.tmpu_, 8L); c0 <- .py_getitem(.tmpu_, 9L); c1 <- .py_getitem(.tmpu_, 10L); c2 <- .py_getitem(.tmpu_, 11L); c3 <- .py_getitem(.tmpu_, 12L); d0 <- .py_getitem(.tmpu_, 13L); d1 <- .py_getitem(.tmpu_, 14L); d2 <- .py_getitem(.tmpu_, 15L); d3 <- .py_getitem(.tmpu_, 16L)
})
.py_try({
TF1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
TF2 <- TransferFunction(k, 1L, s)
})
.py_try({
TF3 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
TF4 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(tau, s), 1L), s)
})
.py_try({
TF5 <- TransferFunction(k, p, s)
})
.py_try({
TF6 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(k, s), p), s)
})
p_are_floats_equals <- function(n1, n2) {
  "\n    Helper function to compare two floats.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return(.py_method(Float(n1), "is_same", Float(n2), .py_isclose))
}
p_are_complex_equals <- function(n1, n2) {
  "\n    Helper function to compare two complex numbers.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return((.py_truthy(p_are_floats_equals(re(n1), re(n2))) && .py_truthy(p_are_floats_equals(im(n1), im(n2)))))
}
.py_try({
A <- .py_Matrix(list(list(a0, a1), list(a2, a3)))
})
.py_try({
B <- .py_Matrix(list(list(b0, b1), list(b2, b3)))
})
.py_try({
C <- .py_Matrix(list(list(c0, c1), list(c1, c2), list(c2, c3)))
})
.py_try({
D <- .py_Matrix(list(list(d0, d1), list(d1, d2), list(d2, d3)))
})
.py_try({
SS <- DiscreteStateSpace(A, B, C, D, T)
})
.py_try({
SS_neg <- .py_neg(SS)
})
.py_try({
state_mat <- .py_Matrix(list(list(-1L, 1L), list(1L, -1L)))
})
.py_try({
input_mat <- .py_Matrix(list(1L, -1L))
})
.py_try({
output_mat <- .py_Matrix(list(list(-1L, 1L)))
})
.py_try({
feedforward_mat <- .py_Matrix(list(1L))
})
.py_try({
system <- DiscreteStateSpace(state_mat, input_mat, output_mat, feedforward_mat, T)
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_negation:4379", function() SS_neg)
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_negation:4390", function() .py_neg(system))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_negation:4396", function() .py_neg(SS_neg))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_negation:4397", function() .py_neg(.py_neg(.py_neg(.py_neg(system)))))
}))---")
  expect_identical(got[["lti:test_DiscreteStateSpace_negation:4379"]], "DiscreteStateSpace(Matrix([\\n[a0, a1],\\n[a2, a3]]), Matrix([\\n[b0, b1],\\n[b2, b3]]), Matrix([\\n[-c0, -c1],\\n[-c1, -c2],\\n[-c2, -c3]]), Matrix([\\n[-d0, -d1],\\n[-d1, -d2],\\n[-d2, -d3]]), T)")
  expect_identical(got[["lti:test_DiscreteStateSpace_negation:4390"]], "DiscreteStateSpace(Matrix([\\n[-1,  1],\\n[ 1, -1]]), Matrix([\\n[ 1],\\n[-1]]), Matrix([[1, -1]]), Matrix([[-1]]), T)")
  expect_identical(got[["lti:test_DiscreteStateSpace_negation:4396"]], "DiscreteStateSpace(Matrix([\\n[a0, a1],\\n[a2, a3]]), Matrix([\\n[b0, b1],\\n[b2, b3]]), Matrix([\\n[c0, c1],\\n[c1, c2],\\n[c2, c3]]), Matrix([\\n[d0, d1],\\n[d1, d2],\\n[d2, d3]]), T)")
  expect_identical(got[["lti:test_DiscreteStateSpace_negation:4397"]], "DiscreteStateSpace(Matrix([\\n[-1,  1],\\n[ 1, -1]]), Matrix([\\n[ 1],\\n[-1]]), Matrix([[-1, 1]]), Matrix([[1]]), T)")
})

test_that("lti test_DiscreteStateSpace_functions matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("a, x, b, c, s, g, d,    p, k, tau, zeta, wn, T, z"); a <- .py_getitem(.tmpu_, 0L); x <- .py_getitem(.tmpu_, 1L); b <- .py_getitem(.tmpu_, 2L); c <- .py_getitem(.tmpu_, 3L); s <- .py_getitem(.tmpu_, 4L); g <- .py_getitem(.tmpu_, 5L); d <- .py_getitem(.tmpu_, 6L); p <- .py_getitem(.tmpu_, 7L); k <- .py_getitem(.tmpu_, 8L); tau <- .py_getitem(.tmpu_, 9L); zeta <- .py_getitem(.tmpu_, 10L); wn <- .py_getitem(.tmpu_, 11L); T <- .py_getitem(.tmpu_, 12L); z <- .py_getitem(.tmpu_, 13L)
})
.py_try({
.tmpu_ <- .py_symbols("a0:4, b0:5, c0:4, d0:4"); a0 <- .py_getitem(.tmpu_, 0L); a1 <- .py_getitem(.tmpu_, 1L); a2 <- .py_getitem(.tmpu_, 2L); a3 <- .py_getitem(.tmpu_, 3L); b0 <- .py_getitem(.tmpu_, 4L); b1 <- .py_getitem(.tmpu_, 5L); b2 <- .py_getitem(.tmpu_, 6L); b3 <- .py_getitem(.tmpu_, 7L); b4 <- .py_getitem(.tmpu_, 8L); c0 <- .py_getitem(.tmpu_, 9L); c1 <- .py_getitem(.tmpu_, 10L); c2 <- .py_getitem(.tmpu_, 11L); c3 <- .py_getitem(.tmpu_, 12L); d0 <- .py_getitem(.tmpu_, 13L); d1 <- .py_getitem(.tmpu_, 14L); d2 <- .py_getitem(.tmpu_, 15L); d3 <- .py_getitem(.tmpu_, 16L)
})
.py_try({
TF1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
TF2 <- TransferFunction(k, 1L, s)
})
.py_try({
TF3 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
TF4 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(tau, s), 1L), s)
})
.py_try({
TF5 <- TransferFunction(k, p, s)
})
.py_try({
TF6 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(k, s), p), s)
})
p_are_floats_equals <- function(n1, n2) {
  "\n    Helper function to compare two floats.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return(.py_method(Float(n1), "is_same", Float(n2), .py_isclose))
}
p_are_complex_equals <- function(n1, n2) {
  "\n    Helper function to compare two complex numbers.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return((.py_truthy(p_are_floats_equals(re(n1), re(n2))) && .py_truthy(p_are_floats_equals(im(n1), im(n2)))))
}
.py_try({
A_mat <- .py_Matrix(list(list(.py_flt(-1.5), -2L), list(1L, 0L)))
})
.py_try({
B_mat <- .py_Matrix(list(.py_flt(0.5), 0L))
})
.py_try({
C_mat <- .py_Matrix(list(list(0L, 1L)))
})
.py_try({
D_mat <- .py_Matrix(list(1L))
})
.py_try({
SS1 <- DiscreteStateSpace(A_mat, B_mat, C_mat, D_mat)
})
.py_try({
SS2 <- DiscreteStateSpace(.py_Matrix(list(list(1L, 1L), list(4L, -2L))), .py_Matrix(list(list(0L, 1L), list(0L, 2L))), .py_Matrix(list(list(-1L, 1L), list(1L, -1L))))
})
.py_try({
SS3 <- DiscreteStateSpace(.py_Matrix(list(list(1L, 1L), list(4L, -2L))), .py_Matrix(list(list(1L, -1L), list(1L, -1L))))
})
.py_try({
SS4 <- DiscreteStateSpace(.py_Matrix(list(list(a0, a1), list(a2, a3))), .py_Matrix(list(list(b1), list(b2))), .py_Matrix(list(list(c1, c2))))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_functions:4415", function() .py_method(SS1, "is_observable"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_functions:4416", function() .py_method(SS2, "is_observable"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_functions:4417", function() .py_method(SS1, "observability_matrix"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_functions:4418", function() .py_method(SS2, "observability_matrix"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_functions:4419", function() .py_method(SS1, "observable_subspace"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_functions:4420", function() .py_method(SS2, "observable_subspace"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_functions:4421", function() .py_method(SS4, "observable_subspace"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_functions:4422", function() .py_method(SS1, "unobservable_subspace"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_functions:4423", function() .py_method(SS2, "unobservable_subspace"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_functions:4424", function() .py_method(SS4, "unobservable_subspace"))
})
.py_try({
Qo <- .py_method(.py_method(SS4, "observability_matrix"), "subs", list(.py_tup(a0, 0L), .py_tup(a1, -6L), .py_tup(a2, 1L), .py_tup(a3, -5L), .py_tup(c1, 0L), .py_tup(c2, 1L)))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_functions:4427", function() Qo)
})
.py_try({
ss_obs <- .py_method(StateSpace(.py_Matrix(list(list(1L, 0L, 1L), list(0L, 0L, 0L), list(0L, 0L, -2L))), .py_Matrix(list(1L, 1L, 0L)), .py_attr(.py_Matrix(list(1L, 1L, Rational(1L, 3L))), "T")), "to_observable_form")
})
.py_try({
A_obs <- .py_getitem(.py_attr(ss_obs, "A"), .py_tup(.py_slice(NULL, 2L, NULL), .py_slice(NULL, 2L, NULL)))
})
.py_try({
A_nobs <- .py_getitem(.py_attr(ss_obs, "A"), .py_tup(.py_slice(2L, NULL, NULL), .py_slice(2L, NULL, NULL)))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_functions:4434", function() .py_method(A_nobs, "eigenvals"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_functions:4437", function() .py_method(SS1, "is_controllable"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_functions:4438", function() .py_method(SS3, "is_controllable"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_functions:4439", function() .py_method(SS1, "controllability_matrix"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_functions:4440", function() .py_method(SS3, "controllability_matrix"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_functions:4441", function() .py_method(SS4, "controllability_matrix"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_functions:4443", function() .py_method(SS1, "controllable_subspace"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_functions:4444", function() .py_method(SS3, "controllable_subspace"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_functions:4445", function() .py_method(SS1, "uncontrollable_subspace"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_functions:4446", function() .py_method(SS3, "uncontrollable_subspace"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_functions:4447", function() .py_method(SS4, "uncontrollable_subspace"))
})
.py_try({
Qc <- .py_method(.py_method(SS4, "controllability_matrix"), "subs", list(.py_tup(a0, 0L), .py_tup(a1, 1L), .py_tup(a2, -6L), .py_tup(a3, -5L), .py_tup(b1, 0L), .py_tup(b2, 1L)))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_functions:4450", function() Qc)
})
.py_try({
ss_contr <- .py_method(StateSpace(.py_Matrix(list(list(1L, 0L, 1L), list(0L, 0L, 0L), list(0L, 0L, -2L))), .py_Matrix(list(1L, 1L, 0L)), .py_attr(.py_Matrix(list(1L, 1L, 0L)), "T")), "to_controllable_form")
})
.py_try({
A_contr <- .py_getitem(.py_attr(ss_contr, "A"), .py_tup(.py_slice(NULL, 2L, NULL), .py_slice(NULL, 2L, NULL)))
})
.py_try({
A_ncontr <- .py_getitem(.py_attr(ss_contr, "A"), .py_tup(.py_slice(2L, NULL, NULL), .py_slice(2L, NULL, NULL)))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_functions:4455", function() .py_method(A_ncontr, "eigenvals"))
})
.py_try({
A1 <- .py_Matrix(list(list(0L, 1L), list(1L, 0L)))
})
.py_try({
B1 <- .py_Matrix(list(list(0L), list(1L)))
})
.py_try({
C1 <- .py_Matrix(list(list(0L, 1L)))
})
.py_try({
D1 <- .py_Matrix(list(list(0L)))
})
.py_try({
ss1 <- DiscreteStateSpace(A1, B1, C1, D1)
})
.py_try({
ss2 <- DiscreteStateSpace(.py_Matrix(list(list(1L, 0L), list(0L, 1L))), .py_Matrix(list(list(1L), list(0L))), .py_Matrix(list(list(1L, 0L))), .py_Matrix(list(list(1L))))
})
.py_try({
ss3 <- .py_method(ss1, "append", ss2)
})
.py_try({
ss4 <- .py_method(SS4, "append", ss1)
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_functions:4468", function() .py_attr(ss3, "num_states"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_functions:4469", function() .py_attr(ss3, "num_inputs"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_functions:4470", function() .py_attr(ss3, "num_outputs"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_functions:4471", function() .py_attr(ss3, "state_matrix"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_functions:4473", function() .py_attr(ss3, "input_matrix"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_functions:4474", function() .py_attr(ss3, "output_matrix"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_functions:4475", function() .py_attr(ss3, "feedforward_matrix"))
})
.py_try({
cont_ss1 <- StateSpace(A1, B1, C1, D1)
})
.py_try({
ss1_ <- DiscreteStateSpace(A1, B1, C1, D1, sampling_time = 2L)
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_functions:4479", function() .py_method(ss1, "append", cont_ss1))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_functions:4480", function() .py_method(ss1, "append", ss1_))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_functions:4483", function() .py_attr(ss4, "num_states"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_functions:4484", function() .py_attr(ss4, "num_inputs"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_functions:4485", function() .py_attr(ss4, "num_outputs"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_functions:4486", function() .py_attr(ss4, "state_matrix"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_functions:4488", function() .py_attr(ss4, "input_matrix"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_functions:4489", function() .py_attr(ss4, "output_matrix"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_functions:4490", function() .py_attr(ss4, "feedforward_matrix"))
}))---")
  expect_identical(got[["lti:test_DiscreteStateSpace_functions:4415"]], "True")
  expect_identical(got[["lti:test_DiscreteStateSpace_functions:4416"]], "False")
  expect_identical(got[["lti:test_DiscreteStateSpace_functions:4417"]], "Matrix([[0, 1], [1, 0]])")
  expect_identical(got[["lti:test_DiscreteStateSpace_functions:4418"]], "Matrix([[-1, 1], [1, -1], [3, -3], [-3, 3]])")
  expect_identical(got[["lti:test_DiscreteStateSpace_functions:4419"]], "[Matrix([[1], [0]]), Matrix([[0], [1]])]")
  expect_identical(got[["lti:test_DiscreteStateSpace_functions:4420"]], "[Matrix([[-1], [1]])]")
  expect_identical(got[["lti:test_DiscreteStateSpace_functions:4421"]], "ERR NotImplementedError: Unbservable subspace cannot be determined for symbolic matrices.")
  expect_identical(got[["lti:test_DiscreteStateSpace_functions:4422"]], "[]")
  expect_identical(got[["lti:test_DiscreteStateSpace_functions:4423"]], "[Matrix([[1], [1]])]")
  expect_identical(got[["lti:test_DiscreteStateSpace_functions:4424"]], "ERR NotImplementedError: Unbservable subspace cannot be determined for symbolic matrices.")
  expect_identical(got[["lti:test_DiscreteStateSpace_functions:4427"]], "Matrix([[0, 1], [1, -5]])")
  expect_identical(got[["lti:test_DiscreteStateSpace_functions:4434"]], "{-2: 1}")
  expect_identical(got[["lti:test_DiscreteStateSpace_functions:4437"]], "True")
  expect_identical(got[["lti:test_DiscreteStateSpace_functions:4438"]], "False")
  expect_identical(got[["lti:test_DiscreteStateSpace_functions:4439"]], "Matrix([[0.500000000000000, -0.750000000000000], [0, 0.500000000000000]])")
  expect_identical(got[["lti:test_DiscreteStateSpace_functions:4440"]], "Matrix([[1, -1, 2, -2], [1, -1, 2, -2]])")
  expect_identical(got[["lti:test_DiscreteStateSpace_functions:4441"]], "Matrix([[b1, a0*b1 + a1*b2], [b2, a2*b1 + a3*b2]])")
  expect_identical(got[["lti:test_DiscreteStateSpace_functions:4443"]], "[Matrix([[0.500000000000000], [0]]), Matrix([[-0.750000000000000], [0.500000000000000]])]")
  expect_identical(got[["lti:test_DiscreteStateSpace_functions:4444"]], "[Matrix([[1], [1]])]")
  expect_identical(got[["lti:test_DiscreteStateSpace_functions:4445"]], "[]")
  expect_identical(got[["lti:test_DiscreteStateSpace_functions:4446"]], "[Matrix([[-1], [1]])]")
  expect_identical(got[["lti:test_DiscreteStateSpace_functions:4447"]], "ERR NotImplementedError: Uncontrollable subspace cannot be determined for symbolic matrices.")
  expect_identical(got[["lti:test_DiscreteStateSpace_functions:4450"]], "Matrix([[0, 1], [1, -5]])")
  expect_identical(got[["lti:test_DiscreteStateSpace_functions:4455"]], "{-2: 1}")
  expect_identical(got[["lti:test_DiscreteStateSpace_functions:4468"]], "4")
  expect_identical(got[["lti:test_DiscreteStateSpace_functions:4469"]], "2")
  expect_identical(got[["lti:test_DiscreteStateSpace_functions:4470"]], "2")
  expect_identical(got[["lti:test_DiscreteStateSpace_functions:4471"]], "Matrix([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])")
  expect_identical(got[["lti:test_DiscreteStateSpace_functions:4473"]], "Matrix([[0, 0], [1, 0], [0, 1], [0, 0]])")
  expect_identical(got[["lti:test_DiscreteStateSpace_functions:4474"]], "Matrix([[0, 1, 0, 0], [0, 0, 1, 0]])")
  expect_identical(got[["lti:test_DiscreteStateSpace_functions:4475"]], "Matrix([[0, 0], [0, 1]])")
  expect_identical(got[["lti:test_DiscreteStateSpace_functions:4479"]], "ERR TypeError: \\n                Incompatible systems detected.\\n                All systems should be either continuous-time or discrete-time.\\n                Found: DiscreteStateSpace(Matrix([\\n[0, 1],\\n[1, 0]]), Matrix([\\n[0],\\n[1]]), Matrix([[0, 1]]), Matrix([[0]]), 1) and StateSpace(Matrix([\\n[0, 1],\\n[1, 0]]), Matrix([\\n[0],\\n[1]]), Matrix([[0, 1]]), Matrix([[0]]))")
  expect_identical(got[["lti:test_DiscreteStateSpace_functions:4480"]], "ERR TypeError: \\n                Incompatible sampling times detected.\\n                All systems should have the same sampling time.\\n                Found sampling time: 1 and 2")
  expect_identical(got[["lti:test_DiscreteStateSpace_functions:4483"]], "4")
  expect_identical(got[["lti:test_DiscreteStateSpace_functions:4484"]], "2")
  expect_identical(got[["lti:test_DiscreteStateSpace_functions:4485"]], "2")
  expect_identical(got[["lti:test_DiscreteStateSpace_functions:4486"]], "Matrix([[a0, a1, 0, 0], [a2, a3, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]])")
  expect_identical(got[["lti:test_DiscreteStateSpace_functions:4488"]], "Matrix([[b1, 0], [b2, 0], [0, 0], [0, 1]])")
  expect_identical(got[["lti:test_DiscreteStateSpace_functions:4489"]], "Matrix([[c1, c2, 0, 0], [0, 0, 0, 1]])")
  expect_identical(got[["lti:test_DiscreteStateSpace_functions:4490"]], "Matrix([[0, 0], [0, 0]])")
})

test_that("lti test_DiscreteStateSpace_series matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("a, x, b, c, s, g, d,    p, k, tau, zeta, wn, T, z"); a <- .py_getitem(.tmpu_, 0L); x <- .py_getitem(.tmpu_, 1L); b <- .py_getitem(.tmpu_, 2L); c <- .py_getitem(.tmpu_, 3L); s <- .py_getitem(.tmpu_, 4L); g <- .py_getitem(.tmpu_, 5L); d <- .py_getitem(.tmpu_, 6L); p <- .py_getitem(.tmpu_, 7L); k <- .py_getitem(.tmpu_, 8L); tau <- .py_getitem(.tmpu_, 9L); zeta <- .py_getitem(.tmpu_, 10L); wn <- .py_getitem(.tmpu_, 11L); T <- .py_getitem(.tmpu_, 12L); z <- .py_getitem(.tmpu_, 13L)
})
.py_try({
.tmpu_ <- .py_symbols("a0:4, b0:5, c0:4, d0:4"); a0 <- .py_getitem(.tmpu_, 0L); a1 <- .py_getitem(.tmpu_, 1L); a2 <- .py_getitem(.tmpu_, 2L); a3 <- .py_getitem(.tmpu_, 3L); b0 <- .py_getitem(.tmpu_, 4L); b1 <- .py_getitem(.tmpu_, 5L); b2 <- .py_getitem(.tmpu_, 6L); b3 <- .py_getitem(.tmpu_, 7L); b4 <- .py_getitem(.tmpu_, 8L); c0 <- .py_getitem(.tmpu_, 9L); c1 <- .py_getitem(.tmpu_, 10L); c2 <- .py_getitem(.tmpu_, 11L); c3 <- .py_getitem(.tmpu_, 12L); d0 <- .py_getitem(.tmpu_, 13L); d1 <- .py_getitem(.tmpu_, 14L); d2 <- .py_getitem(.tmpu_, 15L); d3 <- .py_getitem(.tmpu_, 16L)
})
.py_try({
TF1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
TF2 <- TransferFunction(k, 1L, s)
})
.py_try({
TF3 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
TF4 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(tau, s), 1L), s)
})
.py_try({
TF5 <- TransferFunction(k, p, s)
})
.py_try({
TF6 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(k, s), p), s)
})
p_are_floats_equals <- function(n1, n2) {
  "\n    Helper function to compare two floats.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return(.py_method(Float(n1), "is_same", Float(n2), .py_isclose))
}
p_are_complex_equals <- function(n1, n2) {
  "\n    Helper function to compare two complex numbers.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return((.py_truthy(p_are_floats_equals(re(n1), re(n2))) && .py_truthy(p_are_floats_equals(im(n1), im(n2)))))
}
.py_try({
a1 <- .py_Matrix(list(list(0L, 1L), list(1L, 0L)))
})
.py_try({
b1 <- .py_Matrix(list(list(0L), list(1L)))
})
.py_try({
c1 <- .py_Matrix(list(list(0L, 1L)))
})
.py_try({
d1 <- .py_Matrix(list(list(0L)))
})
.py_try({
a2 <- .py_Matrix(list(list(1L, 0L), list(0L, 1L)))
})
.py_try({
b2 <- .py_Matrix(list(list(1L), list(0L)))
})
.py_try({
c2 <- .py_Matrix(list(list(1L, 0L)))
})
.py_try({
d2 <- .py_Matrix(list(list(1L)))
})
.py_try({
ss1 <- DiscreteStateSpace(a1, b1, c1, d1, T)
})
.py_try({
ss2 <- DiscreteStateSpace(a2, b2, c2, d2, T)
})
.py_try({
tf1 <- DiscreteTransferFunction(s, .py_add(s, 1L), s, T)
})
.py_try({
ser1 <- Series(ss1, ss2)
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_series:4507", function() ser1)
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_series:4514", function() .py_method(ser1, "doit"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_series:4518", function() .py_attr(ser1, "num_inputs"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_series:4519", function() .py_attr(ser1, "num_outputs"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_series:4520", function() .py_method(ser1, "rewrite", DiscreteTransferFunction))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_series:4522", function() .py_method(ser1, "rewrite", TransferFunction))
})
.py_try({
ser2 <- Series(ss1)
})
.py_try({
ser3 <- Series(ser2, ss2)
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_series:4525", function() .py_method(ser3, "doit"))
})
.py_try({
ser_tf <- Series(tf1, ss1)
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_series:4529", function() ser_tf)
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_series:4533", function() .py_method(ser_tf, "doit"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_series:4537", function() .py_method(ser_tf, "rewrite", DiscreteTransferFunction))
})
.py_try({
a3 <- .py_Matrix(list(list(4L, 1L), list(2L, -3L)))
})
.py_try({
b3 <- .py_Matrix(list(list(5L, 2L), list(-3L, -3L)))
})
.py_try({
c3 <- .py_Matrix(list(list(2L, -4L), list(0L, 1L)))
})
.py_try({
d3 <- .py_Matrix(list(list(3L, 2L), list(1L, -1L)))
})
.py_try({
a4 <- .py_Matrix(list(list(-3L, 4L, 2L), list(-1L, -3L, 0L), list(2L, 5L, 3L)))
})
.py_try({
b4 <- .py_Matrix(list(list(1L, 4L), list(-3L, -3L), list(-2L, 1L)))
})
.py_try({
c4 <- .py_Matrix(list(list(4L, 2L, -3L), list(1L, 4L, 3L)))
})
.py_try({
d4 <- .py_Matrix(list(list(-2L, 4L), list(0L, 1L)))
})
.py_try({
ss3 <- DiscreteStateSpace(a3, b3, c3, d3, T)
})
.py_try({
ss4 <- DiscreteStateSpace(a4, b4, c4, d4, T)
})
.py_try({
ser4 <- MIMOSeries(ss3, ss4)
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_series:4552", function() ser4)
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_series:4571", function() .py_method(ser4, "doit"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_series:4590", function() .py_attr(ser4, "num_inputs"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_series:4591", function() .py_attr(ser4, "num_outputs"))
})
.py_try({
ser5 <- MIMOSeries(ss3)
})
.py_try({
ser6 <- MIMOSeries(ser5, ss4)
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_series:4594", function() .py_method(ser6, "doit"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_series:4595", function() .py_method(ser6, "rewrite", TransferFunctionMatrix))
})
.py_try({
tf2 <- DiscreteTransferFunction(1L, s, s, T)
})
.py_try({
tf3 <- DiscreteTransferFunction(1L, .py_add(s, 1L), s, T)
})
.py_try({
tf4 <- DiscreteTransferFunction(s, .py_add(s, 2L), s, T)
})
.py_try({
tfm <- TransferFunctionMatrix(list(list(tf1, tf2), list(tf3, tf4)))
})
.py_try({
ser6 <- MIMOSeries(ss3, tfm)
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_series:4602", function() ser6)
})
.py_try({
ss5 <- DiscreteStateSpace(a1, b1, c1, d1, .py_flt(0.2))
})
.py_try({
cont_ss1 <- StateSpace(a2, b2, c2, d2)
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_series:4621", function() Series(ss5, ss2))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_series:4622", function() Series(ss1, cont_ss1))
})
.py_try({
ss6 <- DiscreteStateSpace(a3, b3, c3, d3)
})
.py_try({
cont_ss2 <- StateSpace(a4, b4, c4, d4)
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_series:4627", function() MIMOSeries(ss6, ss4))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_series:4628", function() MIMOSeries(ss3, cont_ss2))
}))---")
  expect_identical(got[["lti:test_DiscreteStateSpace_series:4507"]], "Series(DiscreteStateSpace(Matrix([\\n[0, 1],\\n[1, 0]]), Matrix([\\n[0],\\n[1]]), Matrix([[0, 1]]), Matrix([[0]]), T), DiscreteStateSpace(Matrix([\\n[1, 0],\\n[0, 1]]), Matrix([\\n[1],\\n[0]]), Matrix([[1, 0]]), Matrix([[1]]), T))")
  expect_identical(got[["lti:test_DiscreteStateSpace_series:4514"]], "DiscreteStateSpace(Matrix([\\n[0, 1, 0, 0],\\n[1, 0, 0, 0],\\n[0, 1, 1, 0],\\n[0, 0, 0, 1]]), Matrix([\\n[0],\\n[1],\\n[0],\\n[0]]), Matrix([[0, 1, 1, 0]]), Matrix([[0]]), T)")
  expect_identical(got[["lti:test_DiscreteStateSpace_series:4518"]], "1")
  expect_identical(got[["lti:test_DiscreteStateSpace_series:4519"]], "1")
  expect_identical(got[["lti:test_DiscreteStateSpace_series:4520"]], "DiscreteTransferFunction(z^2, z^3 - z^2 - z + 1, z, T)")
  expect_identical(got[["lti:test_DiscreteStateSpace_series:4522"]], "ERR TypeError: \\n            The discrete state space model cannot be rewritten as a\\n            continuous-time transfer function model.\\n            ")
  expect_identical(got[["lti:test_DiscreteStateSpace_series:4525"]], "DiscreteStateSpace(Matrix([\\n[0, 1, 0, 0],\\n[1, 0, 0, 0],\\n[0, 1, 1, 0],\\n[0, 0, 0, 1]]), Matrix([\\n[0],\\n[1],\\n[0],\\n[0]]), Matrix([[0, 1, 1, 0]]), Matrix([[0]]), T)")
  expect_identical(got[["lti:test_DiscreteStateSpace_series:4529"]], "Series(DiscreteTransferFunction(s, s + 1, s, T), DiscreteStateSpace(Matrix([\\n[0, 1],\\n[1, 0]]), Matrix([\\n[0],\\n[1]]), Matrix([[0, 1]]), Matrix([[0]]), T))")
  expect_identical(got[["lti:test_DiscreteStateSpace_series:4533"]], "DiscreteStateSpace(Matrix([\\n[-1, 0, 0],\\n[ 0, 0, 1],\\n[-1, 1, 0]]), Matrix([\\n[1],\\n[0],\\n[1]]), Matrix([[0, 0, 1]]), Matrix([[0]]), T)")
  expect_identical(got[["lti:test_DiscreteStateSpace_series:4537"]], "DiscreteTransferFunction(z^2, z^3 + z^2 - z - 1, z, T)")
  expect_identical(got[["lti:test_DiscreteStateSpace_series:4552"]], "MIMOSeries(DiscreteStateSpace(Matrix([\\n[4,  1],\\n[2, -3]]), Matrix([\\n[ 5,  2],\\n[-3, -3]]), Matrix([\\n[2, -4],\\n[0,  1]]), Matrix([\\n[3,  2],\\n[1, -1]]), T), DiscreteStateSpace(Matrix([\\n[-3,  4, 2],\\n[-1, -3, 0],\\n[ 2,  5, 3]]), Matrix([\\n[ 1,  4],\\n[-3, -3],\\n[-2,  1]]), Matrix([\\n[4, 2, -3],\\n[1, 4,  3]]), Matrix([\\n[-2, 4],\\n[ 0, 1]]), T))")
  expect_identical(got[["lti:test_DiscreteStateSpace_series:4571"]], "DiscreteStateSpace(Matrix([\\n[ 4,  1,  0,  0, 0],\\n[ 2, -3,  0,  0, 0],\\n[ 2,  0, -3,  4, 2],\\n[-6,  9, -1, -3, 0],\\n[-4,  9,  2,  5, 3]]), Matrix([\\n[  5,  2],\\n[ -3, -3],\\n[  7, -2],\\n[-12, -3],\\n[ -5, -5]]), Matrix([\\n[-4, 12, 4, 2, -3],\\n[ 0,  1, 1, 4,  3]]), Matrix([\\n[-2, -8],\\n[ 1, -1]]), T)")
  expect_identical(got[["lti:test_DiscreteStateSpace_series:4590"]], "2")
  expect_identical(got[["lti:test_DiscreteStateSpace_series:4591"]], "2")
  expect_identical(got[["lti:test_DiscreteStateSpace_series:4594"]], "DiscreteStateSpace(Matrix([\\n[ 4,  1,  0,  0, 0],\\n[ 2, -3,  0,  0, 0],\\n[ 2,  0, -3,  4, 2],\\n[-6,  9, -1, -3, 0],\\n[-4,  9,  2,  5, 3]]), Matrix([\\n[  5,  2],\\n[ -3, -3],\\n[  7, -2],\\n[-12, -3],\\n[ -5, -5]]), Matrix([\\n[-4, 12, 4, 2, -3],\\n[ 0,  1, 1, 4,  3]]), Matrix([\\n[-2, -8],\\n[ 1, -1]]), T)")
  expect_identical(got[["lti:test_DiscreteStateSpace_series:4595"]], "[[DiscreteTransferFunction(-2*z^5 - 41*z^4 + 130*z^3 + 2383*z^2 + 4038*z - 42024, z^5 + 2*z^4 - 26*z^3 - 74*z^2 + 167*z + 574, z, T), DiscreteTransferFunction(-8*z^5 - 59*z^4 + 361*z^3 + 2080*z^2 + 804*z - 34236, z^5 + 2*z^4 - 26*z^3 - 74*z^2 + 167*z + 574, z, T)], [DiscreteTransferFunction(z^5 - 57*z^4 - 630*z^3 - 231*z^2 + 7719*z + 11778, z^5 + 2*z^4 - 26*z^3 - 74*z^2 + 167*z + 574, z, T), DiscreteTransferFunction(-z^5 - 34*z^4 - 356*z^3 - 53*z^2 + 5498*z + 8968, z^5 + 2*z^4 - 26*z^3 - 74*z^2 + 167*z + 574, z, T)]]")
  expect_identical(got[["lti:test_DiscreteStateSpace_series:4602"]], "MIMOSeries(DiscreteStateSpace(Matrix([\\n[4,  1],\\n[2, -3]]), Matrix([\\n[ 5,  2],\\n[-3, -3]]), Matrix([\\n[2, -4],\\n[0,  1]]), Matrix([\\n[3,  2],\\n[1, -1]]), T), TransferFunctionMatrix(((DiscreteTransferFunction(s, s + 1, s, T), DiscreteTransferFunction(1, s, s, T)), (DiscreteTransferFunction(1, s + 1, s, T), DiscreteTransferFunction(s, s + 2, s, T)))))")
  expect_identical(got[["lti:test_DiscreteStateSpace_series:4621"]], "ERR TypeError: \\n                Incompatible sampling times detected.\\n                All systems should have the same sampling time.\\n                Found sampling time: 0.200000000000000 and T")
  expect_identical(got[["lti:test_DiscreteStateSpace_series:4622"]], "ERR TypeError: \\n                Incompatible systems detected.\\n                All systems should be either continuous-time or discrete-time.\\n                Found: DiscreteStateSpace(Matrix([\\n[0, 1],\\n[1, 0]]), Matrix([\\n[0],\\n[1]]), Matrix([[0, 1]]), Matrix([[0]]), T) and StateSpace(Matrix([\\n[1, 0],\\n[0, 1]]), Matrix([\\n[1],\\n[0]]), Matrix([[1, 0]]), Matrix([[1]]))")
  expect_identical(got[["lti:test_DiscreteStateSpace_series:4627"]], "ERR TypeError: \\n                Incompatible sampling times detected.\\n                All systems should have the same sampling time.\\n                Found sampling time: 1 and T")
  expect_identical(got[["lti:test_DiscreteStateSpace_series:4628"]], "ERR TypeError: \\n                Incompatible systems detected.\\n                All systems should be either continuous-time or discrete-time.\\n                Found: DiscreteStateSpace(Matrix([\\n[4,  1],\\n[2, -3]]), Matrix([\\n[ 5,  2],\\n[-3, -3]]), Matrix([\\n[2, -4],\\n[0,  1]]), Matrix([\\n[3,  2],\\n[1, -1]]), T) and StateSpace(Matrix([\\n[-3,  4, 2],\\n[-1, -3, 0],\\n[ 2,  5, 3]]), Matrix([\\n[ 1,  4],\\n[-3, -3],\\n[-2,  1]]), Matrix([\\n[4, 2, -3],\\n[1, 4,  3]]), Matrix([\\n[-2, 4],\\n[ 0, 1]]))")
})

test_that("lti test_DiscreteStateSpace_parallel matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("a, x, b, c, s, g, d,    p, k, tau, zeta, wn, T, z"); a <- .py_getitem(.tmpu_, 0L); x <- .py_getitem(.tmpu_, 1L); b <- .py_getitem(.tmpu_, 2L); c <- .py_getitem(.tmpu_, 3L); s <- .py_getitem(.tmpu_, 4L); g <- .py_getitem(.tmpu_, 5L); d <- .py_getitem(.tmpu_, 6L); p <- .py_getitem(.tmpu_, 7L); k <- .py_getitem(.tmpu_, 8L); tau <- .py_getitem(.tmpu_, 9L); zeta <- .py_getitem(.tmpu_, 10L); wn <- .py_getitem(.tmpu_, 11L); T <- .py_getitem(.tmpu_, 12L); z <- .py_getitem(.tmpu_, 13L)
})
.py_try({
.tmpu_ <- .py_symbols("a0:4, b0:5, c0:4, d0:4"); a0 <- .py_getitem(.tmpu_, 0L); a1 <- .py_getitem(.tmpu_, 1L); a2 <- .py_getitem(.tmpu_, 2L); a3 <- .py_getitem(.tmpu_, 3L); b0 <- .py_getitem(.tmpu_, 4L); b1 <- .py_getitem(.tmpu_, 5L); b2 <- .py_getitem(.tmpu_, 6L); b3 <- .py_getitem(.tmpu_, 7L); b4 <- .py_getitem(.tmpu_, 8L); c0 <- .py_getitem(.tmpu_, 9L); c1 <- .py_getitem(.tmpu_, 10L); c2 <- .py_getitem(.tmpu_, 11L); c3 <- .py_getitem(.tmpu_, 12L); d0 <- .py_getitem(.tmpu_, 13L); d1 <- .py_getitem(.tmpu_, 14L); d2 <- .py_getitem(.tmpu_, 15L); d3 <- .py_getitem(.tmpu_, 16L)
})
.py_try({
TF1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
TF2 <- TransferFunction(k, 1L, s)
})
.py_try({
TF3 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
TF4 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(tau, s), 1L), s)
})
.py_try({
TF5 <- TransferFunction(k, p, s)
})
.py_try({
TF6 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(k, s), p), s)
})
p_are_floats_equals <- function(n1, n2) {
  "\n    Helper function to compare two floats.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return(.py_method(Float(n1), "is_same", Float(n2), .py_isclose))
}
p_are_complex_equals <- function(n1, n2) {
  "\n    Helper function to compare two complex numbers.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return((.py_truthy(p_are_floats_equals(re(n1), re(n2))) && .py_truthy(p_are_floats_equals(im(n1), im(n2)))))
}
.py_try({
a1 <- .py_Matrix(list(list(0L, 1L), list(1L, 0L)))
})
.py_try({
b1 <- .py_Matrix(list(list(0L), list(1L)))
})
.py_try({
c1 <- .py_Matrix(list(list(0L, 1L)))
})
.py_try({
d1 <- .py_Matrix(list(list(0L)))
})
.py_try({
a2 <- .py_Matrix(list(list(1L, 0L), list(0L, 1L)))
})
.py_try({
b2 <- .py_Matrix(list(list(1L), list(0L)))
})
.py_try({
c2 <- .py_Matrix(list(list(1L, 0L)))
})
.py_try({
d2 <- .py_Matrix(list(list(1L)))
})
.py_try({
ss1 <- DiscreteStateSpace(a1, b1, c1, d1, .py_flt(0.12))
})
.py_try({
ss2 <- DiscreteStateSpace(a2, b2, c2, d2, .py_flt(0.12))
})
.py_try({
p1 <- Parallel(ss1, ss2)
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_parallel:4643", function() p1)
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_parallel:4649", function() .py_method(p1, "doit"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_parallel:4661", function() .py_method(p1, "rewrite", DiscreteTransferFunction))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_parallel:4663", function() .py_method(p1, "rewrite", TransferFunction))
})
.py_try({
tf1 <- DiscreteTransferFunction(s, .py_add(s, 1L), s, .py_flt(0.12))
})
.py_try({
p2 <- Parallel(ss1, tf1)
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_parallel:4668", function() p2)
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_parallel:4674", function() .py_method(p2, "doit"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_parallel:4685", function() .py_method(p2, "rewrite", DiscreteTransferFunction))
})
.py_try({
a3 <- .py_Matrix(list(list(4L, 1L), list(2L, -3L)))
})
.py_try({
b3 <- .py_Matrix(list(list(5L, 2L), list(-3L, -3L)))
})
.py_try({
c3 <- .py_Matrix(list(list(2L, -4L), list(0L, 1L)))
})
.py_try({
d3 <- .py_Matrix(list(list(3L, 2L), list(1L, -1L)))
})
.py_try({
a4 <- .py_Matrix(list(list(-3L, 4L, 2L), list(-1L, -3L, 0L), list(2L, 5L, 3L)))
})
.py_try({
b4 <- .py_Matrix(list(list(1L, 4L), list(-3L, -3L), list(-2L, 1L)))
})
.py_try({
c4 <- .py_Matrix(list(list(4L, 2L, -3L), list(1L, 4L, 3L)))
})
.py_try({
d4 <- .py_Matrix(list(list(-2L, 4L), list(0L, 1L)))
})
.py_try({
ss3 <- DiscreteStateSpace(a3, b3, c3, d3, .py_flt(0.12))
})
.py_try({
ss4 <- DiscreteStateSpace(a4, b4, c4, d4, .py_flt(0.12))
})
.py_try({
p3 <- MIMOParallel(ss3, ss4)
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_parallel:4700", function() p3)
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_parallel:4719", function() .py_method(p3, "doit"))
})
.py_try({
tf2 <- DiscreteTransferFunction(1L, s, s, .py_flt(0.12))
})
.py_try({
tf3 <- DiscreteTransferFunction(1L, .py_add(s, 1L), s, .py_flt(0.12))
})
.py_try({
tf4 <- DiscreteTransferFunction(s, .py_add(s, 2L), s, .py_flt(0.12))
})
.py_try({
tfm <- TransferFunctionMatrix(list(list(tf1, tf2), list(tf3, tf4)))
})
.py_try({
p4 <- MIMOParallel(tfm, ss3)
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_parallel:4744", function() p4)
})
.py_try({
ss5 <- DiscreteStateSpace(a1, b1, c1, d1, .py_flt(0.33))
})
.py_try({
cont_ss1 <- DiscreteStateSpace(a2, b2, c2, d2)
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_parallel:4764", function() Parallel(ss5, ss2))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_parallel:4765", function() Parallel(ss1, cont_ss1))
})
.py_try({
ss6 <- DiscreteStateSpace(a3, b3, c3, d3, .py_flt(0.33))
})
.py_try({
cont_ss2 <- DiscreteStateSpace(a4, b4, c4, d4)
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_parallel:4770", function() MIMOParallel(ss6, ss4))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_parallel:4771", function() MIMOParallel(ss3, cont_ss2))
}))---")
  expect_identical(got[["lti:test_DiscreteStateSpace_parallel:4643"]], "Parallel(DiscreteStateSpace(Matrix([\\n[0, 1],\\n[1, 0]]), Matrix([\\n[0],\\n[1]]), Matrix([[0, 1]]), Matrix([[0]]), 0.12), DiscreteStateSpace(Matrix([\\n[1, 0],\\n[0, 1]]), Matrix([\\n[1],\\n[0]]), Matrix([[1, 0]]), Matrix([[1]]), 0.12))")
  expect_identical(got[["lti:test_DiscreteStateSpace_parallel:4649"]], "DiscreteStateSpace(Matrix([\\n[0, 1, 0, 0],\\n[1, 0, 0, 0],\\n[0, 0, 1, 0],\\n[0, 0, 0, 1]]), Matrix([\\n[0],\\n[1],\\n[1],\\n[0]]), Matrix([[0, 1, 1, 0]]), Matrix([[1]]), 0.12)")
  expect_identical(got[["lti:test_DiscreteStateSpace_parallel:4661"]], "DiscreteTransferFunction(z*(z + 2), z^2 - 1, z, 0.12)")
  expect_identical(got[["lti:test_DiscreteStateSpace_parallel:4663"]], "ERR TypeError: \\n            The discrete state space model cannot be rewritten as a\\n            continuous-time transfer function model.\\n            ")
  expect_identical(got[["lti:test_DiscreteStateSpace_parallel:4668"]], "Parallel(DiscreteStateSpace(Matrix([\\n[0, 1],\\n[1, 0]]), Matrix([\\n[0],\\n[1]]), Matrix([[0, 1]]), Matrix([[0]]), 0.12), DiscreteTransferFunction(s, s + 1, s, 0.12))")
  expect_identical(got[["lti:test_DiscreteStateSpace_parallel:4674"]], "DiscreteStateSpace(Matrix([\\n[0, 1,  0],\\n[1, 0,  0],\\n[0, 0, -1]]), Matrix([\\n[0],\\n[1],\\n[1]]), Matrix([[0, 1, -1]]), Matrix([[1]]), 0.12)")
  expect_identical(got[["lti:test_DiscreteStateSpace_parallel:4685"]], "DiscreteTransferFunction(z^2, z^2 - 1, z, 0.12)")
  expect_identical(got[["lti:test_DiscreteStateSpace_parallel:4700"]], "MIMOParallel(DiscreteStateSpace(Matrix([\\n[4,  1],\\n[2, -3]]), Matrix([\\n[ 5,  2],\\n[-3, -3]]), Matrix([\\n[2, -4],\\n[0,  1]]), Matrix([\\n[3,  2],\\n[1, -1]]), 0.12), DiscreteStateSpace(Matrix([\\n[-3,  4, 2],\\n[-1, -3, 0],\\n[ 2,  5, 3]]), Matrix([\\n[ 1,  4],\\n[-3, -3],\\n[-2,  1]]), Matrix([\\n[4, 2, -3],\\n[1, 4,  3]]), Matrix([\\n[-2, 4],\\n[ 0, 1]]), 0.12))")
  expect_identical(got[["lti:test_DiscreteStateSpace_parallel:4719"]], "DiscreteStateSpace(Matrix([\\n[4,  1,  0,  0, 0],\\n[2, -3,  0,  0, 0],\\n[0,  0, -3,  4, 2],\\n[0,  0, -1, -3, 0],\\n[0,  0,  2,  5, 3]]), Matrix([\\n[ 5,  2],\\n[-3, -3],\\n[ 1,  4],\\n[-3, -3],\\n[-2,  1]]), Matrix([\\n[2, -4, 4, 2, -3],\\n[0,  1, 1, 4,  3]]), Matrix([\\n[1, 6],\\n[1, 0]]), 0.12)")
  expect_identical(got[["lti:test_DiscreteStateSpace_parallel:4744"]], "MIMOParallel(TransferFunctionMatrix(((DiscreteTransferFunction(s, s + 1, s, 0.12), DiscreteTransferFunction(1, s, s, 0.12)), (DiscreteTransferFunction(1, s + 1, s, 0.12), DiscreteTransferFunction(s, s + 2, s, 0.12)))), DiscreteStateSpace(Matrix([\\n[4,  1],\\n[2, -3]]), Matrix([\\n[ 5,  2],\\n[-3, -3]]), Matrix([\\n[2, -4],\\n[0,  1]]), Matrix([\\n[3,  2],\\n[1, -1]]), 0.12))")
  expect_identical(got[["lti:test_DiscreteStateSpace_parallel:4764"]], "ERR TypeError: \\n                Incompatible sampling times detected.\\n                All systems should have the same sampling time.\\n                Found sampling time: 0.330000000000000 and 0.120000000000000")
  expect_identical(got[["lti:test_DiscreteStateSpace_parallel:4765"]], "ERR TypeError: \\n                Incompatible sampling times detected.\\n                All systems should have the same sampling time.\\n                Found sampling time: 0.120000000000000 and 1")
  expect_identical(got[["lti:test_DiscreteStateSpace_parallel:4770"]], "ERR TypeError: \\n                Incompatible sampling times detected.\\n                All systems should have the same sampling time.\\n                Found sampling time: 0.330000000000000 and 0.120000000000000")
  expect_identical(got[["lti:test_DiscreteStateSpace_parallel:4771"]], "ERR TypeError: \\n                Incompatible sampling times detected.\\n                All systems should have the same sampling time.\\n                Found sampling time: 0.120000000000000 and 1")
})

test_that("lti test_DiscreteStateSpace_feedback matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("a, x, b, c, s, g, d,    p, k, tau, zeta, wn, T, z"); a <- .py_getitem(.tmpu_, 0L); x <- .py_getitem(.tmpu_, 1L); b <- .py_getitem(.tmpu_, 2L); c <- .py_getitem(.tmpu_, 3L); s <- .py_getitem(.tmpu_, 4L); g <- .py_getitem(.tmpu_, 5L); d <- .py_getitem(.tmpu_, 6L); p <- .py_getitem(.tmpu_, 7L); k <- .py_getitem(.tmpu_, 8L); tau <- .py_getitem(.tmpu_, 9L); zeta <- .py_getitem(.tmpu_, 10L); wn <- .py_getitem(.tmpu_, 11L); T <- .py_getitem(.tmpu_, 12L); z <- .py_getitem(.tmpu_, 13L)
})
.py_try({
.tmpu_ <- .py_symbols("a0:4, b0:5, c0:4, d0:4"); a0 <- .py_getitem(.tmpu_, 0L); a1 <- .py_getitem(.tmpu_, 1L); a2 <- .py_getitem(.tmpu_, 2L); a3 <- .py_getitem(.tmpu_, 3L); b0 <- .py_getitem(.tmpu_, 4L); b1 <- .py_getitem(.tmpu_, 5L); b2 <- .py_getitem(.tmpu_, 6L); b3 <- .py_getitem(.tmpu_, 7L); b4 <- .py_getitem(.tmpu_, 8L); c0 <- .py_getitem(.tmpu_, 9L); c1 <- .py_getitem(.tmpu_, 10L); c2 <- .py_getitem(.tmpu_, 11L); c3 <- .py_getitem(.tmpu_, 12L); d0 <- .py_getitem(.tmpu_, 13L); d1 <- .py_getitem(.tmpu_, 14L); d2 <- .py_getitem(.tmpu_, 15L); d3 <- .py_getitem(.tmpu_, 16L)
})
.py_try({
TF1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
TF2 <- TransferFunction(k, 1L, s)
})
.py_try({
TF3 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
TF4 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(tau, s), 1L), s)
})
.py_try({
TF5 <- TransferFunction(k, p, s)
})
.py_try({
TF6 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(k, s), p), s)
})
p_are_floats_equals <- function(n1, n2) {
  "\n    Helper function to compare two floats.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return(.py_method(Float(n1), "is_same", Float(n2), .py_isclose))
}
p_are_complex_equals <- function(n1, n2) {
  "\n    Helper function to compare two complex numbers.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return((.py_truthy(p_are_floats_equals(re(n1), re(n2))) && .py_truthy(p_are_floats_equals(im(n1), im(n2)))))
}
.py_try({
a1 <- .py_Matrix(list(list(0L, 1L), list(1L, 0L)))
})
.py_try({
b1 <- .py_Matrix(list(list(0L), list(1L)))
})
.py_try({
c1 <- .py_Matrix(list(list(0L, 1L)))
})
.py_try({
d1 <- .py_Matrix(list(list(0L)))
})
.py_try({
a2 <- .py_Matrix(list(list(1L, 0L), list(0L, 1L)))
})
.py_try({
b2 <- .py_Matrix(list(list(1L), list(0L)))
})
.py_try({
c2 <- .py_Matrix(list(list(1L, 0L)))
})
.py_try({
d2 <- .py_Matrix(list(list(1L)))
})
.py_try({
ss1 <- DiscreteStateSpace(a1, b1, c1, d1, T)
})
.py_try({
ss2 <- DiscreteStateSpace(a2, b2, c2, d2, T)
})
.py_try({
fd1 <- Feedback(ss1, ss2)
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_feedback:4788", function() fd1)
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_feedback:4795", function() .py_method(fd1, "doit"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_feedback:4806", function() .py_method(fd1, "rewrite", DiscreteTransferFunction))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_feedback:4808", function() .py_method(fd1, "rewrite", TransferFunction))
})
.py_try({
fd2 <- Feedback(ss1, ss2, 1L)
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_feedback:4812", function() .py_method(fd2, "doit"))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_feedback:4823", function() .py_method(fd2, "rewrite", DiscreteTransferFunction))
})
.py_try({
tf1 <- DiscreteTransferFunction(s, .py_add(s, 1L), s, T)
})
.py_try({
fd3 <- Feedback(ss1, tf1)
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_feedback:4829", function() fd3)
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_feedback:4835", function() .py_method(fd3, "doit"))
})
.py_try({
a3 <- .py_Matrix(list(list(4L, 1L), list(2L, -3L)))
})
.py_try({
b3 <- .py_Matrix(list(list(5L, 2L), list(-3L, -3L)))
})
.py_try({
c3 <- .py_Matrix(list(list(2L, -4L), list(0L, 1L)))
})
.py_try({
d3 <- .py_Matrix(list(list(3L, 2L), list(1L, -1L)))
})
.py_try({
a4 <- .py_Matrix(list(list(-3L, 4L, 2L), list(-1L, -3L, 0L), list(2L, 5L, 3L)))
})
.py_try({
b4 <- .py_Matrix(list(list(1L, 4L), list(-3L, -3L), list(-2L, 1L)))
})
.py_try({
c4 <- .py_Matrix(list(list(4L, 2L, -3L), list(1L, 4L, 3L)))
})
.py_try({
d4 <- .py_Matrix(list(list(-2L, 4L), list(0L, 1L)))
})
.py_try({
ss3 <- DiscreteStateSpace(a3, b3, c3, d3, T)
})
.py_try({
ss4 <- DiscreteStateSpace(a4, b4, c4, d4, T)
})
.py_try({
fd4 <- MIMOFeedback(ss3, ss4)
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_feedback:4859", function() fd4)
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_feedback:4878", function() .py_method(fd4, "doit"))
})
.py_try({
fd5 <- MIMOFeedback(ss3, ss4, 1L)
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_feedback:4903", function() .py_method(fd5, "doit"))
})
.py_try({
ss5 <- DiscreteStateSpace(a1, b1, c1, d1, .py_flt(0.12))
})
.py_try({
cont_ss1 <- DiscreteStateSpace(a2, b2, c2, d2)
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_feedback:4931", function() Feedback(ss5, ss2))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_feedback:4932", function() Feedback(ss1, cont_ss1))
})
.py_try({
ss6 <- DiscreteStateSpace(a3, b3, c3, d3, .py_flt(0.12))
})
.py_try({
cont_ss2 <- DiscreteStateSpace(a4, b4, c4, d4)
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_feedback:4937", function() MIMOFeedback(ss6, ss4))
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_feedback:4938", function() MIMOFeedback(ss3, cont_ss2))
}))---")
  expect_identical(got[["lti:test_DiscreteStateSpace_feedback:4788"]], "Feedback(DiscreteStateSpace(Matrix([\\n[0, 1],\\n[1, 0]]), Matrix([\\n[0],\\n[1]]), Matrix([[0, 1]]), Matrix([[0]]), T), DiscreteStateSpace(Matrix([\\n[1, 0],\\n[0, 1]]), Matrix([\\n[1],\\n[0]]), Matrix([[1, 0]]), Matrix([[1]]), T), -1)")
  expect_identical(got[["lti:test_DiscreteStateSpace_feedback:4795"]], "DiscreteStateSpace(Matrix([\\n[0,  1,  0, 0],\\n[1, -1, -1, 0],\\n[0,  1,  1, 0],\\n[0,  0,  0, 1]]), Matrix([\\n[0],\\n[1],\\n[0],\\n[0]]), Matrix([[0, 1, 0, 0]]), Matrix([[0]]), T)")
  expect_identical(got[["lti:test_DiscreteStateSpace_feedback:4806"]], "DiscreteTransferFunction(z*(z - 1), z^3 - z + 1, z, T)")
  expect_identical(got[["lti:test_DiscreteStateSpace_feedback:4808"]], "ERR TypeError: \\n            The discrete state space model cannot be rewritten as a\\n            continuous-time transfer function model.\\n            ")
  expect_identical(got[["lti:test_DiscreteStateSpace_feedback:4812"]], "DiscreteStateSpace(Matrix([\\n[0, 1, 0, 0],\\n[1, 1, 1, 0],\\n[0, 1, 1, 0],\\n[0, 0, 0, 1]]), Matrix([\\n[0],\\n[1],\\n[0],\\n[0]]), Matrix([[0, 1, 0, 0]]), Matrix([[0]]), T)")
  expect_identical(got[["lti:test_DiscreteStateSpace_feedback:4823"]], "DiscreteTransferFunction(z*(z - 1), z^3 - 2*z^2 - z + 1, z, T)")
  expect_identical(got[["lti:test_DiscreteStateSpace_feedback:4829"]], "Feedback(DiscreteStateSpace(Matrix([\\n[0, 1],\\n[1, 0]]), Matrix([\\n[0],\\n[1]]), Matrix([[0, 1]]), Matrix([[0]]), T), DiscreteTransferFunction(s, s + 1, s, T), -1)")
  expect_identical(got[["lti:test_DiscreteStateSpace_feedback:4835"]], "DiscreteStateSpace(Matrix([\\n[0,  1,  0],\\n[1, -1,  1],\\n[0,  1, -1]]), Matrix([\\n[0],\\n[1],\\n[0]]), Matrix([[0, 1, 0]]), Matrix([[0]]), T)")
  expect_identical(got[["lti:test_DiscreteStateSpace_feedback:4859"]], "MIMOFeedback(DiscreteStateSpace(Matrix([\\n[4,  1],\\n[2, -3]]), Matrix([\\n[ 5,  2],\\n[-3, -3]]), Matrix([\\n[2, -4],\\n[0,  1]]), Matrix([\\n[3,  2],\\n[1, -1]]), T), DiscreteStateSpace(Matrix([\\n[-3,  4, 2],\\n[-1, -3, 0],\\n[ 2,  5, 3]]), Matrix([\\n[ 1,  4],\\n[-3, -3],\\n[-2,  1]]), Matrix([\\n[4, 2, -3],\\n[1, 4,  3]]), Matrix([\\n[-2, 4],\\n[ 0, 1]]), T), -1)")
  expect_identical(got[["lti:test_DiscreteStateSpace_feedback:4878"]], "DiscreteStateSpace(Matrix([\\n[   3,  -3/4, -15/4, -37/2, -15],\\n[ 7/2, -39/8,   9/8,  39/4,   9],\\n[   3, -41/4, -45/4, -51/2, -19],\\n[-9/2, 129/8,  73/8, 171/4,  36],\\n[-3/2,  47/8,  31/8,  85/4,  18]]), Matrix([\\n[-1/4,  19/4],\\n[ 3/8, -21/8],\\n[ 1/4,  29/4],\\n[ 3/8, -93/8],\\n[ 5/8, -35/8]]), Matrix([\\n[  1, -15/4,  -7/4, -21/2, -9],\\n[1/2, -13/8, -13/8, -19/4, -3]]), Matrix([\\n[-1/4, 11/4],\\n[ 1/8,  9/8]]), T)")
  expect_identical(got[["lti:test_DiscreteStateSpace_feedback:4903"]], "DiscreteStateSpace(Matrix([\\n[  4/7,    62/7,    1,  -8, -69/7],\\n[ 32/7, -135/14, -3/2,   3,  36/7],\\n[-10/7,    41/7,   -4, -12, -97/7],\\n[ 12/7, -111/14, -5/2,  18, 171/7],\\n[  2/7,  -29/14, -1/2,  10,  81/7]]), Matrix([\\n[   6/7, -17/7],\\n[ -9/14, 15/14],\\n[   6/7, -31/7],\\n[-27/14, 87/14],\\n[-15/14, 25/14]]), Matrix([\\n[-2/7,  11/7,    1, -4, -39/7],\\n[-2/7, 15/14, -1/2, -3, -18/7]]), Matrix([\\n[ 4/7,   -9/7],\\n[1/14, -11/14]]), T)")
  expect_identical(got[["lti:test_DiscreteStateSpace_feedback:4931"]], "ERR TypeError: \\n                Incompatible sampling times detected.\\n                All systems should have the same sampling time.\\n                Found sampling time: 0.120000000000000 and T")
  expect_identical(got[["lti:test_DiscreteStateSpace_feedback:4932"]], "ERR TypeError: \\n                Incompatible sampling times detected.\\n                All systems should have the same sampling time.\\n                Found sampling time: T and 1")
  expect_identical(got[["lti:test_DiscreteStateSpace_feedback:4937"]], "ERR TypeError: \\n                Incompatible sampling times detected.\\n                All systems should have the same sampling time.\\n                Found sampling time: 0.120000000000000 and T")
  expect_identical(got[["lti:test_DiscreteStateSpace_feedback:4938"]], "ERR TypeError: \\n                Incompatible sampling times detected.\\n                All systems should have the same sampling time.\\n                Found sampling time: T and 1")
})

test_that("lti test_DiscreteStateSpace_stability matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("a, x, b, c, s, g, d,    p, k, tau, zeta, wn, T, z"); a <- .py_getitem(.tmpu_, 0L); x <- .py_getitem(.tmpu_, 1L); b <- .py_getitem(.tmpu_, 2L); c <- .py_getitem(.tmpu_, 3L); s <- .py_getitem(.tmpu_, 4L); g <- .py_getitem(.tmpu_, 5L); d <- .py_getitem(.tmpu_, 6L); p <- .py_getitem(.tmpu_, 7L); k <- .py_getitem(.tmpu_, 8L); tau <- .py_getitem(.tmpu_, 9L); zeta <- .py_getitem(.tmpu_, 10L); wn <- .py_getitem(.tmpu_, 11L); T <- .py_getitem(.tmpu_, 12L); z <- .py_getitem(.tmpu_, 13L)
})
.py_try({
.tmpu_ <- .py_symbols("a0:4, b0:5, c0:4, d0:4"); a0 <- .py_getitem(.tmpu_, 0L); a1 <- .py_getitem(.tmpu_, 1L); a2 <- .py_getitem(.tmpu_, 2L); a3 <- .py_getitem(.tmpu_, 3L); b0 <- .py_getitem(.tmpu_, 4L); b1 <- .py_getitem(.tmpu_, 5L); b2 <- .py_getitem(.tmpu_, 6L); b3 <- .py_getitem(.tmpu_, 7L); b4 <- .py_getitem(.tmpu_, 8L); c0 <- .py_getitem(.tmpu_, 9L); c1 <- .py_getitem(.tmpu_, 10L); c2 <- .py_getitem(.tmpu_, 11L); c3 <- .py_getitem(.tmpu_, 12L); d0 <- .py_getitem(.tmpu_, 13L); d1 <- .py_getitem(.tmpu_, 14L); d2 <- .py_getitem(.tmpu_, 15L); d3 <- .py_getitem(.tmpu_, 16L)
})
.py_try({
TF1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_mul(.py_mul(2L, zeta), wn), s)), .py_pow(wn, 2L)), s)
})
.py_try({
TF2 <- TransferFunction(k, 1L, s)
})
.py_try({
TF3 <- TransferFunction(.py_sub(.py_mul(a2, p), s), .py_add(.py_mul(a2, s), p), s)
})
.py_try({
TF4 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(tau, s), 1L), s)
})
.py_try({
TF5 <- TransferFunction(k, p, s)
})
.py_try({
TF6 <- TransferFunction(.py_add(.py_mul(k, s), p), .py_add(.py_mul(k, s), p), s)
})
p_are_floats_equals <- function(n1, n2) {
  "\n    Helper function to compare two floats.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return(.py_method(Float(n1), "is_same", Float(n2), .py_isclose))
}
p_are_complex_equals <- function(n1, n2) {
  "\n    Helper function to compare two complex numbers.\n    Returns True if they are approximately equal, otherwise False.\n\n    "
  return((.py_truthy(p_are_floats_equals(re(n1), re(n2))) && .py_truthy(p_are_floats_equals(im(n1), im(n2)))))
}
.py_try({
k <- .py_symbols("k")
})
.py_try({
B <- .py_Matrix(list(1L, 0L, 0L))
})
.py_try({
C <- .py_Matrix(list(list(0L, 1L, 0L)))
})
.py_try({
D <- .py_Matrix(list(0L))
})
.py_try({
A1 <- .py_Matrix(list(list(0L, 1L, 0L), list(0L, 0L, 1L), list(.py_sub(k, 1L), .py_mul(-2L, k), -1L)))
})
.py_try({
ss1 <- DiscreteStateSpace(A1, B, C, D)
})
.py_try({
ineq <- .py_method(ss1, "get_asymptotic_stability_conditions")
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_stability:4949", function() ineq)
})
.py_try({
A2 <- .py_Matrix(list(list(1L, 0L, 0L), list(0L, -1L, k), list(0L, 0L, -1L)))
})
.py_try({
ss2 <- DiscreteStateSpace(A2, B, C, D)
})
.py_try({
ineq <- .py_method(ss2, "get_asymptotic_stability_conditions")
})
.py_try({
.py_rec("lti:test_DiscreteStateSpace_stability:4955", function() ineq)
}))---")
  expect_identical(got[["lti:test_DiscreteStateSpace_stability:4949"]], "[-15*k^2 + 20*k - 5 > 0, -8*k^2 - 8*k + 8 > 0, 3*k^2 + 8*k - 3 > 0]")
  expect_identical(got[["lti:test_DiscreteStateSpace_stability:4955"]], "[False]")
})

test_that("routh test_table matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
s <- .py_symbols("s")
})
.py_try({
.tmpu_ <- .py_symbols("b_0 b_1 b_2 b_3 b_4"); b0 <- .py_getitem(.tmpu_, 0L); b1 <- .py_getitem(.tmpu_, 1L); b2 <- .py_getitem(.tmpu_, 2L); b3 <- .py_getitem(.tmpu_, 3L); b4 <- .py_getitem(.tmpu_, 4L)
})
.py_try({
p1 <- .py_add(.py_add(.py_add(.py_add(.py_mul(b4, .py_pow(s, 4L)), .py_mul(b3, .py_pow(s, 3L))), .py_mul(b2, .py_pow(s, 2L))), .py_mul(b1, s)), b0)
})
.py_try({
t1 <- RouthHurwitz(p1, s)
})
.py_try({
expected1 <- .py_Matrix(list(list(b4, b2, b0), list(b3, b1, 0L), list(.py_add(.py_div(.py_mul(.py_neg(b1), b4), b3), b2), b0, 0L), list(.py_div(.py_add(.py_mul(b0, .py_pow(b3, 2L)), .py_mul(b1, .py_sub(.py_mul(b1, b4), .py_mul(b2, b3)))), .py_sub(.py_mul(b1, b4), .py_mul(b2, b3))), 0L, 0L), list(b0, 0L, 0L)))
})
.py_try({
.py_rec("routh:test_table:21", function() .py_method(t1, "equals", expected1))
})
.py_try({
expected1_1 <- .py_Matrix(list(list(b4), list(b3), list(.py_add(.py_div(.py_mul(.py_neg(b1), b4), b3), b2)), list(.py_div(.py_add(.py_mul(b0, .py_pow(b3, 2L)), .py_mul(b1, .py_sub(.py_mul(b1, b4), .py_mul(b2, b3)))), .py_sub(.py_mul(b1, b4), .py_mul(b2, b3)))), list(b0)))
})
.py_try({
.py_rec("routh:test_table:28", function() .py_method(.py_getitem(t1, .py_tup(.py_slice(NULL, NULL, NULL), 0L)), "equals", expected1_1))
})
.py_try({
.py_rec("routh:test_table:29", function() .py_attr(t1, "zero_row_case"))
})
.py_try({
.py_rec("routh:test_table:30", function() .py_attr(t1, "auxiliary_polynomials"))
})
.py_try({
p2 <- .py_add(.py_add(.py_add(.py_add(.py_pow(s, 4L), .py_pow(s, 3L)), .py_mul(3L, .py_pow(s, 2L))), .py_mul(3L, s)), 3L)
})
.py_try({
t2 <- RouthHurwitz(p2, s)
})
.py_try({
expected2 <- .py_Matrix(list(list(1L, 3L, 3L), list(1L, 3L, 0L), list(-3L, 3L, 0L), list(4L, 0L, 0L), list(3L, 0L, 0L)))
})
.py_try({
.py_rec("routh:test_table:40", function() .py_method(t2, "equals", expected2))
})
.py_try({
.py_rec("routh:test_table:41", function() .py_attr(t2, "zero_row_case"))
})
.py_try({
p3 <- .py_add(.py_add(.py_add(.py_add(.py_add(.py_add(.py_pow(s, 6L), .py_mul(2L, .py_pow(s, 5L))), .py_mul(8L, .py_pow(s, 4L))), .py_mul(12L, .py_pow(s, 3L))), .py_mul(20L, .py_pow(s, 2L))), .py_mul(16L, s)), 16L)
})
.py_try({
t3 <- RouthHurwitz(p3, s)
})
.py_try({
expected3 <- .py_Matrix(list(list(1L, 8L, 20L, 16L), list(2L, 12L, 16L, 0L), list(2L, 12L, 16L, 0L), list(8L, 24L, 0L, 0L), list(6L, 16L, 0L, 0L), list(.py_div(8L, 3L), 0L, 0L, 0L), list(16L, 0L, 0L, 0L)))
})
.py_try({
.py_rec("routh:test_table:51", function() .py_method(t3, "equals", expected3))
})
.py_try({
.py_rec("routh:test_table:52", function() .py_attr(t3, "zero_row_case"))
})
.py_try({
.py_rec("routh:test_table:53", function() .py_attr(t3, "auxiliary_polynomials"))
}))---")
  expect_identical(got[["routh:test_table:21"]], "True")
  expect_identical(got[["routh:test_table:28"]], "True")
  expect_identical(got[["routh:test_table:29"]], "False")
  expect_identical(got[["routh:test_table:30"]], "None")
  expect_identical(got[["routh:test_table:40"]], "True")
  expect_identical(got[["routh:test_table:41"]], "False")
  expect_identical(got[["routh:test_table:51"]], "True")
  expect_identical(got[["routh:test_table:52"]], "True")
  expect_identical(got[["routh:test_table:53"]], "[Poly(2*s^4 + 12*s^2 + 16, s, domain='ZZ')]")
})

test_that("control_plots test_errors matches SymPy", {
  got <- .physmisc_battery_run(r"---(s <- Symbol("s")
p <- Symbol("p")
a <- Symbol("a")
.py_try({
matplotlib <- import_module("matplotlib", import_kwargs = .py_dict(list("fromlist"), list(list("pyplot"))), catch = .py_tup(RuntimeError))
})
.py_try({
numpy <- import_module("numpy")
})
.py_try({
tf1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(p, 2L), .py_mul(.py_flt(0.5), p)), 2L), p)
})
.py_try({
tf2 <- TransferFunction(p, .py_add(.py_add(.py_mul(6L, .py_pow(p, 2L)), .py_mul(3L, p)), 1L), p)
})
.py_try({
tf3 <- TransferFunction(p, .py_sub(.py_pow(p, 3L), 1L), p)
})
.py_try({
tf4 <- TransferFunction(10L, .py_pow(p, 3L), p)
})
.py_try({
tf5 <- TransferFunction(5L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(2L, s)), 10L), s)
})
.py_try({
tf6 <- TransferFunction(1L, 1L, s)
})
.py_try({
tf7 <- TransferFunction(.py_add(.py_add(.py_add(.py_mul(.py_mul(4L, s), 3L), .py_mul(9L, .py_pow(s, 2L))), .py_mul(.py_flt(0.1), s)), 11L), .py_add(.py_add(.py_mul(8L, .py_pow(s, 6L)), .py_mul(9L, .py_pow(s, 4L))), 11L), s)
})
.py_try({
tf8 <- TransferFunction(5L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_add(2L, sym_I), s)), 10L), s)
})
.py_try({
ser1 <- Series(tf4, TransferFunction(1L, .py_sub(p, 5L), p))
})
.py_try({
ser2 <- Series(tf3, TransferFunction(p, .py_add(p, 2L), p))
})
.py_try({
par1 <- Parallel(tf1, tf2)
})
p_to_tuple <- function(a, b) {
  return(.py_tup(.py_tuple(a), .py_tuple(b)))
}
p_trim_tuple <- function(a, b) {
  .tmpu_ <- p_to_tuple(a, b); a <- .py_getitem(.tmpu_, 0L); b <- .py_getitem(.tmpu_, 1L)
  return(.py_tup(.py_tuple(.py_add(.py_add(.py_getitem(a, .py_slice(0L, 2L, NULL)), .py_getitem(a, .py_slice(.py_floordiv(.py_len(a), 2L), .py_add(.py_floordiv(.py_len(a), 2L), 1L), NULL))), .py_getitem(a, .py_slice(-2L, NULL, NULL)))), .py_tuple(.py_add(.py_add(.py_getitem(b, .py_slice(0L, 2L, NULL)), .py_getitem(b, .py_slice(.py_floordiv(.py_len(b), 2L), .py_add(.py_floordiv(.py_len(b), 2L), 1L), NULL))), .py_getitem(b, .py_slice(-2L, NULL, NULL))))))
}
y_coordinate_equality <- function(plot_data_func, evalf_func, system) {
  "Checks whether the y-coordinate value of the plotted\n    data point is equal to the value of the function at a\n    particular x."
  .tmpu_ <- plot_data_func(system); x <- .py_getitem(.tmpu_, 0L); y <- .py_getitem(.tmpu_, 1L)
  .tmpu_ <- p_trim_tuple(x, y); x <- .py_getitem(.tmpu_, 0L); y <- .py_getitem(.tmpu_, 1L)
  y_exp <- .py_tuple(.py_comp(x, function(.it) { x_i <- .it; list(evalf_func(system, x_i)) }))
  return(.py_all(.py_comp(.py_zip(y_exp, y), function(.it) { y_exp_i <- .py_getitem(.it, 0L); y_i <- .py_getitem(.it, 1L); list(.py_lt(Abs(.py_sub(y_exp_i, y_i)), .py_flt(1e-08))) })))
}
check_point_accuracy <- function(a, b) {
  return(.py_all(.py_comp(.py_zip(a, b), function(.it) { p_ <- .it; list(.py_call(.py_isclose, .py_splat(p_), rel_tol = .py_flt(0.1), abs_tol = .py_flt(1e-06))) })))
}
.py_try({
if (.py_truthy((!.py_truthy(matplotlib)))) {
  .py_try({
  NULL
  })
}
})
.py_try({
tfm <- TransferFunctionMatrix(list(list(tf6, tf5), list(tf5, tf6)))
})
.py_try({
expr <- .py_div(1L, .py_sub(.py_pow(s, 2L), 1L))
})
.py_try({
.py_rec("control_plots:test_errors:69", function() pole_zero_plot(tfm))
})
.py_try({
.py_rec("control_plots:test_errors:70", function() pole_zero_numerical_data(expr))
})
.py_try({
.py_rec("control_plots:test_errors:71", function() impulse_response_plot(expr))
})
.py_try({
.py_rec("control_plots:test_errors:72", function() impulse_response_numerical_data(tfm))
})
.py_try({
.py_rec("control_plots:test_errors:73", function() step_response_plot(tfm))
})
.py_try({
.py_rec("control_plots:test_errors:74", function() step_response_numerical_data(expr))
})
.py_try({
.py_rec("control_plots:test_errors:75", function() ramp_response_plot(expr))
})
.py_try({
.py_rec("control_plots:test_errors:76", function() ramp_response_numerical_data(tfm))
})
.py_try({
.py_rec("control_plots:test_errors:77", function() bode_plot(tfm))
})
.py_try({
tf_a <- TransferFunction(a, .py_add(s, 1L), s)
})
.py_try({
.py_rec("control_plots:test_errors:81", function() pole_zero_plot(tf_a))
})
.py_try({
.py_rec("control_plots:test_errors:82", function() pole_zero_numerical_data(tf_a))
})
.py_try({
.py_rec("control_plots:test_errors:83", function() impulse_response_plot(tf_a))
})
.py_try({
.py_rec("control_plots:test_errors:84", function() impulse_response_numerical_data(tf_a))
})
.py_try({
.py_rec("control_plots:test_errors:85", function() step_response_plot(tf_a))
})
.py_try({
.py_rec("control_plots:test_errors:86", function() step_response_numerical_data(tf_a))
})
.py_try({
.py_rec("control_plots:test_errors:87", function() ramp_response_plot(tf_a))
})
.py_try({
.py_rec("control_plots:test_errors:88", function() ramp_response_numerical_data(tf_a))
})
.py_try({
.py_rec("control_plots:test_errors:89", function() bode_plot(tf_a))
})
.py_try({
.py_rec("control_plots:test_errors:92", function() impulse_response_plot(tf1, lower_limit = -1L))
})
.py_try({
.py_rec("control_plots:test_errors:93", function() step_response_plot(tf1, lower_limit = .py_flt(-0.1)))
})
.py_try({
.py_rec("control_plots:test_errors:94", function() ramp_response_plot(tf1, lower_limit = .py_div(-4L, 3L)))
})
.py_try({
.py_rec("control_plots:test_errors:97", function() ramp_response_plot(tf1, slope = .py_flt(-0.1)))
})
.py_try({
.py_rec("control_plots:test_errors:100", function() bode_plot(tf1, freq_unit = "hz"))
})
.py_try({
.py_rec("control_plots:test_errors:101", function() bode_plot(tf1, phase_unit = "degree"))
}))---")
  expect_identical(got[["control_plots:test_errors:69"]], "ERR NotImplementedError: Only SISO LTI systems are currently supported.")
  expect_identical(got[["control_plots:test_errors:70"]], "ERR NotImplementedError: Only SISO LTI systems are currently supported.")
  expect_identical(got[["control_plots:test_errors:71"]], "ERR NotImplementedError: Only SISO LTI systems are currently supported.")
  expect_identical(got[["control_plots:test_errors:72"]], "ERR NotImplementedError: Only SISO LTI systems are currently supported.")
  expect_identical(got[["control_plots:test_errors:73"]], "ERR NotImplementedError: Only SISO LTI systems are currently supported.")
  expect_identical(got[["control_plots:test_errors:74"]], "ERR NotImplementedError: Only SISO LTI systems are currently supported.")
  expect_identical(got[["control_plots:test_errors:75"]], "ERR NotImplementedError: Only SISO LTI systems are currently supported.")
  expect_identical(got[["control_plots:test_errors:76"]], "ERR NotImplementedError: Only SISO LTI systems are currently supported.")
  expect_identical(got[["control_plots:test_errors:77"]], "ERR NotImplementedError: Only SISO LTI systems are currently supported.")
  expect_identical(got[["control_plots:test_errors:81"]], "ERR ValueError: Extra degree of freedom found. Make sure that there are no free symbols in the dynamical system other than the variable of Laplace transform.")
  expect_identical(got[["control_plots:test_errors:82"]], "ERR ValueError: Extra degree of freedom found. Make sure that there are no free symbols in the dynamical system other than the variable of Laplace transform.")
  expect_identical(got[["control_plots:test_errors:83"]], "ERR ValueError: Extra degree of freedom found. Make sure that there are no free symbols in the dynamical system other than the variable of Laplace transform.")
  expect_identical(got[["control_plots:test_errors:84"]], "ERR ValueError: Extra degree of freedom found. Make sure that there are no free symbols in the dynamical system other than the variable of Laplace transform.")
  expect_identical(got[["control_plots:test_errors:85"]], "ERR ValueError: Extra degree of freedom found. Make sure that there are no free symbols in the dynamical system other than the variable of Laplace transform.")
  expect_identical(got[["control_plots:test_errors:86"]], "ERR ValueError: Extra degree of freedom found. Make sure that there are no free symbols in the dynamical system other than the variable of Laplace transform.")
  expect_identical(got[["control_plots:test_errors:87"]], "ERR ValueError: Extra degree of freedom found. Make sure that there are no free symbols in the dynamical system other than the variable of Laplace transform.")
  expect_identical(got[["control_plots:test_errors:88"]], "ERR ValueError: Extra degree of freedom found. Make sure that there are no free symbols in the dynamical system other than the variable of Laplace transform.")
  expect_identical(got[["control_plots:test_errors:89"]], "ERR ValueError: Extra degree of freedom found. Make sure that there are no free symbols in the dynamical system other than the variable of Laplace transform.")
  expect_identical(got[["control_plots:test_errors:92"]], "ERR ValueError: Lower limit of time must be greater than or equal to zero.")
  expect_identical(got[["control_plots:test_errors:93"]], "ERR ValueError: Lower limit of time must be greater than or equal to zero.")
  expect_identical(got[["control_plots:test_errors:94"]], "ERR ValueError: Lower limit of time must be greater than or equal to zero.")
  expect_identical(got[["control_plots:test_errors:97"]], "ERR ValueError: Slope must be greater than or equal to zero.")
  expect_identical(got[["control_plots:test_errors:100"]], "ERR ValueError: Only \"rad/sec\" and \"Hz\" are accepted frequency units.")
  expect_identical(got[["control_plots:test_errors:101"]], "ERR ValueError: Only \"rad\" and \"deg\" are accepted phase units.")
})

test_that("control_plots test_pole_zero matches SymPy", {
  got <- .physmisc_battery_run(r"---(s <- Symbol("s")
p <- Symbol("p")
a <- Symbol("a")
.py_try({
matplotlib <- import_module("matplotlib", import_kwargs = .py_dict(list("fromlist"), list(list("pyplot"))), catch = .py_tup(RuntimeError))
})
.py_try({
numpy <- import_module("numpy")
})
.py_try({
tf1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(p, 2L), .py_mul(.py_flt(0.5), p)), 2L), p)
})
.py_try({
tf2 <- TransferFunction(p, .py_add(.py_add(.py_mul(6L, .py_pow(p, 2L)), .py_mul(3L, p)), 1L), p)
})
.py_try({
tf3 <- TransferFunction(p, .py_sub(.py_pow(p, 3L), 1L), p)
})
.py_try({
tf4 <- TransferFunction(10L, .py_pow(p, 3L), p)
})
.py_try({
tf5 <- TransferFunction(5L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(2L, s)), 10L), s)
})
.py_try({
tf6 <- TransferFunction(1L, 1L, s)
})
.py_try({
tf7 <- TransferFunction(.py_add(.py_add(.py_add(.py_mul(.py_mul(4L, s), 3L), .py_mul(9L, .py_pow(s, 2L))), .py_mul(.py_flt(0.1), s)), 11L), .py_add(.py_add(.py_mul(8L, .py_pow(s, 6L)), .py_mul(9L, .py_pow(s, 4L))), 11L), s)
})
.py_try({
tf8 <- TransferFunction(5L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_add(2L, sym_I), s)), 10L), s)
})
.py_try({
ser1 <- Series(tf4, TransferFunction(1L, .py_sub(p, 5L), p))
})
.py_try({
ser2 <- Series(tf3, TransferFunction(p, .py_add(p, 2L), p))
})
.py_try({
par1 <- Parallel(tf1, tf2)
})
p_to_tuple <- function(a, b) {
  return(.py_tup(.py_tuple(a), .py_tuple(b)))
}
p_trim_tuple <- function(a, b) {
  .tmpu_ <- p_to_tuple(a, b); a <- .py_getitem(.tmpu_, 0L); b <- .py_getitem(.tmpu_, 1L)
  return(.py_tup(.py_tuple(.py_add(.py_add(.py_getitem(a, .py_slice(0L, 2L, NULL)), .py_getitem(a, .py_slice(.py_floordiv(.py_len(a), 2L), .py_add(.py_floordiv(.py_len(a), 2L), 1L), NULL))), .py_getitem(a, .py_slice(-2L, NULL, NULL)))), .py_tuple(.py_add(.py_add(.py_getitem(b, .py_slice(0L, 2L, NULL)), .py_getitem(b, .py_slice(.py_floordiv(.py_len(b), 2L), .py_add(.py_floordiv(.py_len(b), 2L), 1L), NULL))), .py_getitem(b, .py_slice(-2L, NULL, NULL))))))
}
y_coordinate_equality <- function(plot_data_func, evalf_func, system) {
  "Checks whether the y-coordinate value of the plotted\n    data point is equal to the value of the function at a\n    particular x."
  .tmpu_ <- plot_data_func(system); x <- .py_getitem(.tmpu_, 0L); y <- .py_getitem(.tmpu_, 1L)
  .tmpu_ <- p_trim_tuple(x, y); x <- .py_getitem(.tmpu_, 0L); y <- .py_getitem(.tmpu_, 1L)
  y_exp <- .py_tuple(.py_comp(x, function(.it) { x_i <- .it; list(evalf_func(system, x_i)) }))
  return(.py_all(.py_comp(.py_zip(y_exp, y), function(.it) { y_exp_i <- .py_getitem(.it, 0L); y_i <- .py_getitem(.it, 1L); list(.py_lt(Abs(.py_sub(y_exp_i, y_i)), .py_flt(1e-08))) })))
}
check_point_accuracy <- function(a, b) {
  return(.py_all(.py_comp(.py_zip(a, b), function(.it) { p_ <- .it; list(.py_call(.py_isclose, .py_splat(p_), rel_tol = .py_flt(0.1), abs_tol = .py_flt(1e-06))) })))
}
pz_tester <- function(sys, expected_value) {
  .tmpu_ <- pole_zero_numerical_data(sys); p_z <- .py_getitem(.tmpu_, 0L); p_p <- .py_getitem(.tmpu_, 1L)
  z_check <- .py_all_close(p_z, .py_getitem(expected_value, 0L))
  p_check <- .py_all_close(p_p, .py_getitem(expected_value, 1L))
  return((.py_truthy(p_check) && .py_truthy(z_check)))
}
.py_try({
exp1 <- list(list(), list(.py_sub(.py_flt(-0.24999999999999994), complex(real = 0.0, imaginary = 1.3919410907075054)), .py_add(.py_flt(-0.24999999999999994), complex(real = 0.0, imaginary = 1.3919410907075054))))
})
.py_try({
exp2 <- list(list(.py_flt(0.0)), list(.py_sub(.py_flt(-0.25), complex(real = 0.0, imaginary = 0.3227486121839514)), .py_add(.py_flt(-0.25), complex(real = 0.0, imaginary = 0.3227486121839514))))
})
.py_try({
exp3 <- list(list(.py_flt(0.0)), list(.py_add(.py_flt(0.9999999999999998), complex(real = 0.0, imaginary = 0.0)), .py_sub(.py_flt(-0.5000000000000004), complex(real = 0.0, imaginary = 0.8660254037844395)), .py_add(.py_flt(-0.5000000000000004), complex(real = 0.0, imaginary = 0.8660254037844395))))
})
.py_try({
exp4 <- list(list(), list(.py_flt(0.0), .py_flt(0.0), .py_flt(0.0), .py_flt(5.0)))
})
.py_try({
exp5 <- list(list(.py_flt(-5.645751311064592), .py_flt(-0.5000000000000008), .py_flt(-0.3542486889354093)), list(.py_sub(.py_flt(-0.24999999999999986), complex(real = 0.0, imaginary = 0.32274861218395134)), .py_add(.py_flt(-0.2499999999999998), complex(real = 0.0, imaginary = 0.32274861218395134)), .py_sub(.py_flt(-0.24999999999999986), complex(real = 0.0, imaginary = 1.3919410907075052)), .py_add(.py_flt(-0.2499999999999998), complex(real = 0.0, imaginary = 1.3919410907075052))))
})
.py_try({
exp6 <- list(list(), list(.py_sub(.py_flt(-1.1641600331447917), complex(real = 0.0, imaginary = 3.545808351896439)), .py_add(.py_flt(-0.8358399668552097), complex(real = 0.0, imaginary = 2.5458083518964383))))
})
.py_try({
.py_rec("control_plots:test_pole_zero:125", function() pz_tester(tf1, exp1))
})
.py_try({
.py_rec("control_plots:test_pole_zero:126", function() pz_tester(tf2, exp2))
})
.py_try({
.py_rec("control_plots:test_pole_zero:127", function() pz_tester(tf3, exp3))
})
.py_try({
.py_rec("control_plots:test_pole_zero:128", function() pz_tester(ser1, exp4))
})
.py_try({
.py_rec("control_plots:test_pole_zero:130", function() pz_tester(tf8, exp6))
}))---")
  expect_identical(got[["control_plots:test_pole_zero:125"]], "True")
  expect_identical(got[["control_plots:test_pole_zero:126"]], "True")
  expect_identical(got[["control_plots:test_pole_zero:127"]], "True")
  expect_identical(got[["control_plots:test_pole_zero:128"]], "True")
  expect_identical(got[["control_plots:test_pole_zero:130"]], "True")
})

test_that("control_plots test_bode matches SymPy", {
  got <- .physmisc_battery_run(r"---(s <- Symbol("s")
p <- Symbol("p")
a <- Symbol("a")
.py_try({
matplotlib <- import_module("matplotlib", import_kwargs = .py_dict(list("fromlist"), list(list("pyplot"))), catch = .py_tup(RuntimeError))
})
.py_try({
numpy <- import_module("numpy")
})
.py_try({
tf1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(p, 2L), .py_mul(.py_flt(0.5), p)), 2L), p)
})
.py_try({
tf2 <- TransferFunction(p, .py_add(.py_add(.py_mul(6L, .py_pow(p, 2L)), .py_mul(3L, p)), 1L), p)
})
.py_try({
tf3 <- TransferFunction(p, .py_sub(.py_pow(p, 3L), 1L), p)
})
.py_try({
tf4 <- TransferFunction(10L, .py_pow(p, 3L), p)
})
.py_try({
tf5 <- TransferFunction(5L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(2L, s)), 10L), s)
})
.py_try({
tf6 <- TransferFunction(1L, 1L, s)
})
.py_try({
tf7 <- TransferFunction(.py_add(.py_add(.py_add(.py_mul(.py_mul(4L, s), 3L), .py_mul(9L, .py_pow(s, 2L))), .py_mul(.py_flt(0.1), s)), 11L), .py_add(.py_add(.py_mul(8L, .py_pow(s, 6L)), .py_mul(9L, .py_pow(s, 4L))), 11L), s)
})
.py_try({
tf8 <- TransferFunction(5L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_add(2L, sym_I), s)), 10L), s)
})
.py_try({
ser1 <- Series(tf4, TransferFunction(1L, .py_sub(p, 5L), p))
})
.py_try({
ser2 <- Series(tf3, TransferFunction(p, .py_add(p, 2L), p))
})
.py_try({
par1 <- Parallel(tf1, tf2)
})
p_to_tuple <- function(a, b) {
  return(.py_tup(.py_tuple(a), .py_tuple(b)))
}
p_trim_tuple <- function(a, b) {
  .tmpu_ <- p_to_tuple(a, b); a <- .py_getitem(.tmpu_, 0L); b <- .py_getitem(.tmpu_, 1L)
  return(.py_tup(.py_tuple(.py_add(.py_add(.py_getitem(a, .py_slice(0L, 2L, NULL)), .py_getitem(a, .py_slice(.py_floordiv(.py_len(a), 2L), .py_add(.py_floordiv(.py_len(a), 2L), 1L), NULL))), .py_getitem(a, .py_slice(-2L, NULL, NULL)))), .py_tuple(.py_add(.py_add(.py_getitem(b, .py_slice(0L, 2L, NULL)), .py_getitem(b, .py_slice(.py_floordiv(.py_len(b), 2L), .py_add(.py_floordiv(.py_len(b), 2L), 1L), NULL))), .py_getitem(b, .py_slice(-2L, NULL, NULL))))))
}
y_coordinate_equality <- function(plot_data_func, evalf_func, system) {
  "Checks whether the y-coordinate value of the plotted\n    data point is equal to the value of the function at a\n    particular x."
  .tmpu_ <- plot_data_func(system); x <- .py_getitem(.tmpu_, 0L); y <- .py_getitem(.tmpu_, 1L)
  .tmpu_ <- p_trim_tuple(x, y); x <- .py_getitem(.tmpu_, 0L); y <- .py_getitem(.tmpu_, 1L)
  y_exp <- .py_tuple(.py_comp(x, function(.it) { x_i <- .it; list(evalf_func(system, x_i)) }))
  return(.py_all(.py_comp(.py_zip(y_exp, y), function(.it) { y_exp_i <- .py_getitem(.it, 0L); y_i <- .py_getitem(.it, 1L); list(.py_lt(Abs(.py_sub(y_exp_i, y_i)), .py_flt(1e-08))) })))
}
check_point_accuracy <- function(a, b) {
  return(.py_all(.py_comp(.py_zip(a, b), function(.it) { p_ <- .it; list(.py_call(.py_isclose, .py_splat(p_), rel_tol = .py_flt(0.1), abs_tol = .py_flt(1e-06))) })))
}
.py_try({
if (.py_truthy((!.py_truthy(numpy)))) {
  .py_try({
  NULL
  })
}
})
bode_phase_evalf <- function(system, point) {
  expr <- .py_method(system, "to_expr")
  p_w <- Dummy("w", real = TRUE)
  w_expr <- .py_method(expr, "subs", .py_dict(list(.py_attr(system, "var")), list(.py_mul(sym_I, p_w))))
  return(.py_method(.py_method(.py_arg(w_expr), "subs", .py_dict(list(p_w), list(point))), "evalf"))
}
bode_mag_evalf <- function(system, point) {
  expr <- .py_method(system, "to_expr")
  p_w <- Dummy("w", real = TRUE)
  w_expr <- .py_method(expr, "subs", .py_dict(list(.py_attr(system, "var")), list(.py_mul(sym_I, p_w))))
  return(.py_mul(20L, .py_method(.py_method(.py_log(Abs(w_expr), 10L), "subs", .py_dict(list(p_w), list(point))), "evalf")))
}
test_bode_data <- function(sys) {
  return((.py_truthy(y_coordinate_equality(bode_magnitude_numerical_data, bode_mag_evalf, sys)) && .py_truthy(y_coordinate_equality(bode_phase_numerical_data, bode_phase_evalf, sys))))
}
.py_try({
.py_rec("control_plots:test_bode:153", function() test_bode_data(tf1))
})
.py_try({
.py_rec("control_plots:test_bode:154", function() test_bode_data(tf2))
})
.py_try({
.py_rec("control_plots:test_bode:155", function() test_bode_data(tf3))
})
.py_try({
.py_rec("control_plots:test_bode:156", function() test_bode_data(tf4))
})
.py_try({
.py_rec("control_plots:test_bode:157", function() test_bode_data(tf5))
}))---")
  expect_identical(got[["control_plots:test_bode:153"]], "True")
  expect_identical(got[["control_plots:test_bode:154"]], "True")
  expect_identical(got[["control_plots:test_bode:155"]], "True")
  expect_identical(got[["control_plots:test_bode:156"]], "True")
  expect_identical(got[["control_plots:test_bode:157"]], "True")
})

test_that("control_plots test_impulse_response matches SymPy", {
  got <- .physmisc_battery_run(r"---(s <- Symbol("s")
p <- Symbol("p")
a <- Symbol("a")
.py_try({
matplotlib <- import_module("matplotlib", import_kwargs = .py_dict(list("fromlist"), list(list("pyplot"))), catch = .py_tup(RuntimeError))
})
.py_try({
numpy <- import_module("numpy")
})
.py_try({
tf1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(p, 2L), .py_mul(.py_flt(0.5), p)), 2L), p)
})
.py_try({
tf2 <- TransferFunction(p, .py_add(.py_add(.py_mul(6L, .py_pow(p, 2L)), .py_mul(3L, p)), 1L), p)
})
.py_try({
tf3 <- TransferFunction(p, .py_sub(.py_pow(p, 3L), 1L), p)
})
.py_try({
tf4 <- TransferFunction(10L, .py_pow(p, 3L), p)
})
.py_try({
tf5 <- TransferFunction(5L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(2L, s)), 10L), s)
})
.py_try({
tf6 <- TransferFunction(1L, 1L, s)
})
.py_try({
tf7 <- TransferFunction(.py_add(.py_add(.py_add(.py_mul(.py_mul(4L, s), 3L), .py_mul(9L, .py_pow(s, 2L))), .py_mul(.py_flt(0.1), s)), 11L), .py_add(.py_add(.py_mul(8L, .py_pow(s, 6L)), .py_mul(9L, .py_pow(s, 4L))), 11L), s)
})
.py_try({
tf8 <- TransferFunction(5L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_add(2L, sym_I), s)), 10L), s)
})
.py_try({
ser1 <- Series(tf4, TransferFunction(1L, .py_sub(p, 5L), p))
})
.py_try({
ser2 <- Series(tf3, TransferFunction(p, .py_add(p, 2L), p))
})
.py_try({
par1 <- Parallel(tf1, tf2)
})
p_to_tuple <- function(a, b) {
  return(.py_tup(.py_tuple(a), .py_tuple(b)))
}
p_trim_tuple <- function(a, b) {
  .tmpu_ <- p_to_tuple(a, b); a <- .py_getitem(.tmpu_, 0L); b <- .py_getitem(.tmpu_, 1L)
  return(.py_tup(.py_tuple(.py_add(.py_add(.py_getitem(a, .py_slice(0L, 2L, NULL)), .py_getitem(a, .py_slice(.py_floordiv(.py_len(a), 2L), .py_add(.py_floordiv(.py_len(a), 2L), 1L), NULL))), .py_getitem(a, .py_slice(-2L, NULL, NULL)))), .py_tuple(.py_add(.py_add(.py_getitem(b, .py_slice(0L, 2L, NULL)), .py_getitem(b, .py_slice(.py_floordiv(.py_len(b), 2L), .py_add(.py_floordiv(.py_len(b), 2L), 1L), NULL))), .py_getitem(b, .py_slice(-2L, NULL, NULL))))))
}
y_coordinate_equality <- function(plot_data_func, evalf_func, system) {
  "Checks whether the y-coordinate value of the plotted\n    data point is equal to the value of the function at a\n    particular x."
  .tmpu_ <- plot_data_func(system); x <- .py_getitem(.tmpu_, 0L); y <- .py_getitem(.tmpu_, 1L)
  .tmpu_ <- p_trim_tuple(x, y); x <- .py_getitem(.tmpu_, 0L); y <- .py_getitem(.tmpu_, 1L)
  y_exp <- .py_tuple(.py_comp(x, function(.it) { x_i <- .it; list(evalf_func(system, x_i)) }))
  return(.py_all(.py_comp(.py_zip(y_exp, y), function(.it) { y_exp_i <- .py_getitem(.it, 0L); y_i <- .py_getitem(.it, 1L); list(.py_lt(Abs(.py_sub(y_exp_i, y_i)), .py_flt(1e-08))) })))
}
check_point_accuracy <- function(a, b) {
  return(.py_all(.py_comp(.py_zip(a, b), function(.it) { p_ <- .it; list(.py_call(.py_isclose, .py_splat(p_), rel_tol = .py_flt(0.1), abs_tol = .py_flt(1e-06))) })))
}
.py_try({
if (.py_truthy((!.py_truthy(numpy)))) {
  .py_try({
  NULL
  })
}
})
impulse_res_tester <- function(sys, expected_value) {
  .tmpu_ <- .py_call(p_to_tuple, .py_splat(impulse_response_numerical_data(sys, adaptive = FALSE, n = 10L))); x <- .py_getitem(.tmpu_, 0L); y <- .py_getitem(.tmpu_, 1L)
  x_check <- check_point_accuracy(x, .py_getitem(expected_value, 0L))
  y_check <- check_point_accuracy(y, .py_getitem(expected_value, 1L))
  return((.py_truthy(x_check) && .py_truthy(y_check)))
}
.py_try({
exp1 <- .py_tup(.py_tup(.py_flt(0.0), .py_flt(1.1111111111111112), .py_flt(2.2222222222222223), .py_flt(3.3333333333333335), .py_flt(4.444444444444445), .py_flt(5.555555555555555), .py_flt(6.666666666666667), .py_flt(7.777777777777779), .py_flt(8.88888888888889), .py_flt(10.0)), .py_tup(.py_flt(0.0), .py_flt(0.544019738507865), .py_flt(0.01993849743234938), .py_flt(-0.31140243360893216), .py_flt(-0.022852779906491996), .py_flt(0.1778306498155759), .py_flt(0.01962941084328499), .py_flt(-0.1013115194573652), .py_flt(-0.014975541213105696), .py_flt(0.0575789724730714)))
})
.py_try({
exp2 <- .py_tup(.py_tup(.py_flt(0.0), .py_flt(1.1111111111111112), .py_flt(2.2222222222222223), .py_flt(3.3333333333333335), .py_flt(4.444444444444445), .py_flt(5.555555555555555), .py_flt(6.666666666666667), .py_flt(7.777777777777779), .py_flt(8.88888888888889), .py_flt(10.0)), .py_tup(.py_flt(0.1666666675), .py_flt(0.08389223412935855), .py_flt(0.02338051973475047), .py_flt(-0.014966807776379383), .py_flt(-0.034645954223054234), .py_flt(-0.040560075735512804), .py_flt(-0.037658628907103885), .py_flt(-0.030149507719590022), .py_flt(-0.021162090730736834), .py_flt(-0.012721292737437523)))
})
.py_try({
exp3 <- .py_tup(.py_tup(.py_flt(0.0), .py_flt(1.1111111111111112), .py_flt(2.2222222222222223), .py_flt(3.3333333333333335), .py_flt(4.444444444444445), .py_flt(5.555555555555555), .py_flt(6.666666666666667), .py_flt(7.777777777777779), .py_flt(8.88888888888889), .py_flt(10.0)), .py_tup(.py_flt(4.369893391586999e-09), .py_flt(1.1750333000630964), .py_flt(3.2922404058312473), .py_flt(9.432290008148343), .py_flt(28.37098083007151), .py_flt(86.18577464367974), .py_flt(261.90356653762115), .py_flt(795.6538758627842), .py_flt(2416.9920942096983), .py_flt(7342.159505206647)))
})
.py_try({
exp4 <- .py_tup(.py_tup(.py_flt(0.0), .py_flt(1.1111111111111112), .py_flt(2.2222222222222223), .py_flt(3.3333333333333335), .py_flt(4.444444444444445), .py_flt(5.555555555555555), .py_flt(6.666666666666667), .py_flt(7.777777777777779), .py_flt(8.88888888888889), .py_flt(10.0)), .py_tup(.py_flt(0.0), .py_flt(6.17283950617284), .py_flt(24.69135802469136), .py_flt(55.555555555555564), .py_flt(98.76543209876544), .py_flt(154.320987654321), .py_flt(222.22222222222226), .py_flt(302.46913580246917), .py_flt(395.0617283950618), .py_flt(500.0)))
})
.py_try({
exp5 <- .py_tup(.py_tup(.py_flt(0.0), .py_flt(1.1111111111111112), .py_flt(2.2222222222222223), .py_flt(3.3333333333333335), .py_flt(4.444444444444445), .py_flt(5.555555555555555), .py_flt(6.666666666666667), .py_flt(7.777777777777779), .py_flt(8.88888888888889), .py_flt(10.0)), .py_tup(.py_flt(0.0), .py_flt(-0.10455606138085417), .py_flt(0.06757671513476461), .py_flt(-0.03234567568833768), .py_flt(0.013582514927757873), .py_flt(-0.005273419510705473), .py_flt(0.0019364083003354075), .py_flt(-0.000680070134067832), .py_flt(0.00022969845960406913), .py_flt(-7.476094359583917e-05)))
})
.py_try({
exp6 <- .py_tup(.py_tup(.py_flt(0.0), .py_flt(1.1111111111111112), .py_flt(2.2222222222222223), .py_flt(3.3333333333333335), .py_flt(4.444444444444445), .py_flt(5.555555555555555), .py_flt(6.666666666666667), .py_flt(7.777777777777779), .py_flt(8.88888888888889), .py_flt(10.0)), .py_tup(.py_flt(-6.016699583000218e-09), .py_flt(0.35039802056107394), .py_flt(3.3728423827689884), .py_flt(12.119846079276684), .py_flt(25.86101014293389), .py_flt(29.352480635282088), .py_flt(-30.49475907497664), .py_flt(-273.8717189554019), .py_flt(-863.2381702029659), .py_flt(-1747.0262164682233)))
})
.py_try({
exp7 <- .py_tup(.py_tup(.py_flt(0.0), .py_flt(1.1111111111111112), .py_flt(2.2222222222222223), .py_flt(3.3333333333333335), .py_flt(4.444444444444445), .py_flt(5.555555555555555), .py_flt(6.666666666666667), .py_flt(7.777777777777779), .py_flt(8.88888888888889), .py_flt(10.0)), .py_tup(.py_flt(0.0), .py_flt(18.934638095560974), .py_flt(5346.93244680907), .py_flt(1384609.8718249386), .py_flt(358161126.65801865), .py_flt(92645770015.70108), .py_flt(23964739753087.42), .py_flt(6198974342083139.0), .py_flt(1.603492601616059e+18), .py_flt(4.147764422869658e+20)))
})
.py_try({
.py_rec("control_plots:test_impulse_response:207", function() impulse_res_tester(tf1, exp1))
})
.py_try({
.py_rec("control_plots:test_impulse_response:208", function() impulse_res_tester(tf2, exp2))
})
.py_try({
.py_rec("control_plots:test_impulse_response:209", function() impulse_res_tester(tf3, exp3))
})
.py_try({
.py_rec("control_plots:test_impulse_response:210", function() impulse_res_tester(tf4, exp4))
})
.py_try({
.py_rec("control_plots:test_impulse_response:211", function() impulse_res_tester(tf5, exp5))
})
.py_try({
.py_rec("control_plots:test_impulse_response:212", function() impulse_res_tester(tf7, exp6))
})
.py_try({
.py_rec("control_plots:test_impulse_response:213", function() impulse_res_tester(ser1, exp7))
}))---")
  expect_identical(got[["control_plots:test_impulse_response:207"]], "True")
  expect_identical(got[["control_plots:test_impulse_response:208"]], "True")
  expect_identical(got[["control_plots:test_impulse_response:209"]], "True")
  expect_identical(got[["control_plots:test_impulse_response:210"]], "True")
  expect_identical(got[["control_plots:test_impulse_response:211"]], "True")
  expect_identical(got[["control_plots:test_impulse_response:212"]], "True")
  expect_identical(got[["control_plots:test_impulse_response:213"]], "True")
})

test_that("control_plots test_step_response matches SymPy", {
  got <- .physmisc_battery_run(r"---(s <- Symbol("s")
p <- Symbol("p")
a <- Symbol("a")
.py_try({
matplotlib <- import_module("matplotlib", import_kwargs = .py_dict(list("fromlist"), list(list("pyplot"))), catch = .py_tup(RuntimeError))
})
.py_try({
numpy <- import_module("numpy")
})
.py_try({
tf1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(p, 2L), .py_mul(.py_flt(0.5), p)), 2L), p)
})
.py_try({
tf2 <- TransferFunction(p, .py_add(.py_add(.py_mul(6L, .py_pow(p, 2L)), .py_mul(3L, p)), 1L), p)
})
.py_try({
tf3 <- TransferFunction(p, .py_sub(.py_pow(p, 3L), 1L), p)
})
.py_try({
tf4 <- TransferFunction(10L, .py_pow(p, 3L), p)
})
.py_try({
tf5 <- TransferFunction(5L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(2L, s)), 10L), s)
})
.py_try({
tf6 <- TransferFunction(1L, 1L, s)
})
.py_try({
tf7 <- TransferFunction(.py_add(.py_add(.py_add(.py_mul(.py_mul(4L, s), 3L), .py_mul(9L, .py_pow(s, 2L))), .py_mul(.py_flt(0.1), s)), 11L), .py_add(.py_add(.py_mul(8L, .py_pow(s, 6L)), .py_mul(9L, .py_pow(s, 4L))), 11L), s)
})
.py_try({
tf8 <- TransferFunction(5L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_add(2L, sym_I), s)), 10L), s)
})
.py_try({
ser1 <- Series(tf4, TransferFunction(1L, .py_sub(p, 5L), p))
})
.py_try({
ser2 <- Series(tf3, TransferFunction(p, .py_add(p, 2L), p))
})
.py_try({
par1 <- Parallel(tf1, tf2)
})
p_to_tuple <- function(a, b) {
  return(.py_tup(.py_tuple(a), .py_tuple(b)))
}
p_trim_tuple <- function(a, b) {
  .tmpu_ <- p_to_tuple(a, b); a <- .py_getitem(.tmpu_, 0L); b <- .py_getitem(.tmpu_, 1L)
  return(.py_tup(.py_tuple(.py_add(.py_add(.py_getitem(a, .py_slice(0L, 2L, NULL)), .py_getitem(a, .py_slice(.py_floordiv(.py_len(a), 2L), .py_add(.py_floordiv(.py_len(a), 2L), 1L), NULL))), .py_getitem(a, .py_slice(-2L, NULL, NULL)))), .py_tuple(.py_add(.py_add(.py_getitem(b, .py_slice(0L, 2L, NULL)), .py_getitem(b, .py_slice(.py_floordiv(.py_len(b), 2L), .py_add(.py_floordiv(.py_len(b), 2L), 1L), NULL))), .py_getitem(b, .py_slice(-2L, NULL, NULL))))))
}
y_coordinate_equality <- function(plot_data_func, evalf_func, system) {
  "Checks whether the y-coordinate value of the plotted\n    data point is equal to the value of the function at a\n    particular x."
  .tmpu_ <- plot_data_func(system); x <- .py_getitem(.tmpu_, 0L); y <- .py_getitem(.tmpu_, 1L)
  .tmpu_ <- p_trim_tuple(x, y); x <- .py_getitem(.tmpu_, 0L); y <- .py_getitem(.tmpu_, 1L)
  y_exp <- .py_tuple(.py_comp(x, function(.it) { x_i <- .it; list(evalf_func(system, x_i)) }))
  return(.py_all(.py_comp(.py_zip(y_exp, y), function(.it) { y_exp_i <- .py_getitem(.it, 0L); y_i <- .py_getitem(.it, 1L); list(.py_lt(Abs(.py_sub(y_exp_i, y_i)), .py_flt(1e-08))) })))
}
check_point_accuracy <- function(a, b) {
  return(.py_all(.py_comp(.py_zip(a, b), function(.it) { p_ <- .it; list(.py_call(.py_isclose, .py_splat(p_), rel_tol = .py_flt(0.1), abs_tol = .py_flt(1e-06))) })))
}
.py_try({
if (.py_truthy((!.py_truthy(numpy)))) {
  .py_try({
  NULL
  })
}
})
step_res_tester <- function(sys, expected_value) {
  .tmpu_ <- .py_call(p_to_tuple, .py_splat(step_response_numerical_data(sys, adaptive = FALSE, n = 10L))); x <- .py_getitem(.tmpu_, 0L); y <- .py_getitem(.tmpu_, 1L)
  x_check <- check_point_accuracy(x, .py_getitem(expected_value, 0L))
  y_check <- check_point_accuracy(y, .py_getitem(expected_value, 1L))
  return((.py_truthy(x_check) && .py_truthy(y_check)))
}
.py_try({
exp1 <- .py_tup(.py_tup(.py_flt(0.0), .py_flt(1.1111111111111112), .py_flt(2.2222222222222223), .py_flt(3.3333333333333335), .py_flt(4.444444444444445), .py_flt(5.555555555555555), .py_flt(6.666666666666667), .py_flt(7.777777777777779), .py_flt(8.88888888888889), .py_flt(10.0)), .py_tup(.py_flt(-1.9193285738516863e-08), .py_flt(0.42283495488246126), .py_flt(0.7840485977945262), .py_flt(0.5546841805655717), .py_flt(0.33903033806932087), .py_flt(0.4627251747410237), .py_flt(0.5909907598988051), .py_flt(0.5247213989553071), .py_flt(0.4486997874319281), .py_flt(0.4839358435839171)))
})
.py_try({
exp2 <- .py_tup(.py_tup(.py_flt(0.0), .py_flt(1.1111111111111112), .py_flt(2.2222222222222223), .py_flt(3.3333333333333335), .py_flt(4.444444444444445), .py_flt(5.555555555555555), .py_flt(6.666666666666667), .py_flt(7.777777777777779), .py_flt(8.88888888888889), .py_flt(10.0)), .py_tup(.py_flt(0.0), .py_flt(0.13728409095645816), .py_flt(0.19474559355325086), .py_flt(0.1974909129243011), .py_flt(0.16841657696573073), .py_flt(0.12559777736159378), .py_flt(0.08153828016664713), .py_flt(0.04360471317348958), .py_flt(0.015072994568868221), .py_flt(-0.003636420058445484)))
})
.py_try({
exp3 <- .py_tup(.py_tup(.py_flt(0.0), .py_flt(1.1111111111111112), .py_flt(2.2222222222222223), .py_flt(3.3333333333333335), .py_flt(4.444444444444445), .py_flt(5.555555555555555), .py_flt(6.666666666666667), .py_flt(7.777777777777779), .py_flt(8.88888888888889), .py_flt(10.0)), .py_tup(.py_flt(0.0), .py_flt(0.6314542141914303), .py_flt(2.9356520038101035), .py_flt(9.37731009663807), .py_flt(28.452300356688376), .py_flt(86.25721933273988), .py_flt(261.9236645044672), .py_flt(795.6435410577224), .py_flt(2416.9786984578764), .py_flt(7342.154119725917)))
})
.py_try({
exp4 <- .py_tup(.py_tup(.py_flt(0.0), .py_flt(1.1111111111111112), .py_flt(2.2222222222222223), .py_flt(3.3333333333333335), .py_flt(4.444444444444445), .py_flt(5.555555555555555), .py_flt(6.666666666666667), .py_flt(7.777777777777779), .py_flt(8.88888888888889), .py_flt(10.0)), .py_tup(.py_flt(0.0), .py_flt(2.286236899862826), .py_flt(18.28989519890261), .py_flt(61.72839629629631), .py_flt(146.31916159122088), .py_flt(285.7796124828532), .py_flt(493.8271703703705), .py_flt(784.1792566529494), .py_flt(1170.553292729767), .py_flt(1666.6667)))
})
.py_try({
exp5 <- .py_tup(.py_tup(.py_flt(0.0), .py_flt(1.1111111111111112), .py_flt(2.2222222222222223), .py_flt(3.3333333333333335), .py_flt(4.444444444444445), .py_flt(5.555555555555555), .py_flt(6.666666666666667), .py_flt(7.777777777777779), .py_flt(8.88888888888889), .py_flt(10.0)), .py_tup(.py_flt(-3.999999997894577e-09), .py_flt(0.6720357068882895), .py_flt(0.4429938256137113), .py_flt(0.5182010838004518), .py_flt(0.4944139147159695), .py_flt(0.5016379853883338), .py_flt(0.4995466896527733), .py_flt(0.5001154784851325), .py_flt(0.49997448824584123), .py_flt(0.5000039745919259)))
})
.py_try({
exp6 <- .py_tup(.py_tup(.py_flt(0.0), .py_flt(1.1111111111111112), .py_flt(2.2222222222222223), .py_flt(3.3333333333333335), .py_flt(4.444444444444445), .py_flt(5.555555555555555), .py_flt(6.666666666666667), .py_flt(7.777777777777779), .py_flt(8.88888888888889), .py_flt(10.0)), .py_tup(.py_flt(-1.5433688493882158e-09), .py_flt(0.3428705539937336), .py_flt(1.1253619102202777), .py_flt(3.1849962651016517), .py_flt(9.47532757182671), .py_flt(28.727231099148135), .py_flt(87.29426924860557), .py_flt(265.2138681048606), .py_flt(805.6636260007757), .py_flt(2447.387582370878)))
})
.py_try({
.py_rec("control_plots:test_step_response:256", function() step_res_tester(tf1, exp1))
})
.py_try({
.py_rec("control_plots:test_step_response:257", function() step_res_tester(tf2, exp2))
})
.py_try({
.py_rec("control_plots:test_step_response:258", function() step_res_tester(tf3, exp3))
})
.py_try({
.py_rec("control_plots:test_step_response:259", function() step_res_tester(tf4, exp4))
})
.py_try({
.py_rec("control_plots:test_step_response:260", function() step_res_tester(tf5, exp5))
})
.py_try({
.py_rec("control_plots:test_step_response:261", function() step_res_tester(ser2, exp6))
}))---")
  expect_identical(got[["control_plots:test_step_response:256"]], "True")
  expect_identical(got[["control_plots:test_step_response:257"]], "True")
  expect_identical(got[["control_plots:test_step_response:258"]], "True")
  expect_identical(got[["control_plots:test_step_response:259"]], "True")
  expect_identical(got[["control_plots:test_step_response:260"]], "True")
  expect_identical(got[["control_plots:test_step_response:261"]], "True")
})

test_that("control_plots test_ramp_response matches SymPy", {
  got <- .physmisc_battery_run(r"---(s <- Symbol("s")
p <- Symbol("p")
a <- Symbol("a")
.py_try({
matplotlib <- import_module("matplotlib", import_kwargs = .py_dict(list("fromlist"), list(list("pyplot"))), catch = .py_tup(RuntimeError))
})
.py_try({
numpy <- import_module("numpy")
})
.py_try({
tf1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(p, 2L), .py_mul(.py_flt(0.5), p)), 2L), p)
})
.py_try({
tf2 <- TransferFunction(p, .py_add(.py_add(.py_mul(6L, .py_pow(p, 2L)), .py_mul(3L, p)), 1L), p)
})
.py_try({
tf3 <- TransferFunction(p, .py_sub(.py_pow(p, 3L), 1L), p)
})
.py_try({
tf4 <- TransferFunction(10L, .py_pow(p, 3L), p)
})
.py_try({
tf5 <- TransferFunction(5L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(2L, s)), 10L), s)
})
.py_try({
tf6 <- TransferFunction(1L, 1L, s)
})
.py_try({
tf7 <- TransferFunction(.py_add(.py_add(.py_add(.py_mul(.py_mul(4L, s), 3L), .py_mul(9L, .py_pow(s, 2L))), .py_mul(.py_flt(0.1), s)), 11L), .py_add(.py_add(.py_mul(8L, .py_pow(s, 6L)), .py_mul(9L, .py_pow(s, 4L))), 11L), s)
})
.py_try({
tf8 <- TransferFunction(5L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_add(2L, sym_I), s)), 10L), s)
})
.py_try({
ser1 <- Series(tf4, TransferFunction(1L, .py_sub(p, 5L), p))
})
.py_try({
ser2 <- Series(tf3, TransferFunction(p, .py_add(p, 2L), p))
})
.py_try({
par1 <- Parallel(tf1, tf2)
})
p_to_tuple <- function(a, b) {
  return(.py_tup(.py_tuple(a), .py_tuple(b)))
}
p_trim_tuple <- function(a, b) {
  .tmpu_ <- p_to_tuple(a, b); a <- .py_getitem(.tmpu_, 0L); b <- .py_getitem(.tmpu_, 1L)
  return(.py_tup(.py_tuple(.py_add(.py_add(.py_getitem(a, .py_slice(0L, 2L, NULL)), .py_getitem(a, .py_slice(.py_floordiv(.py_len(a), 2L), .py_add(.py_floordiv(.py_len(a), 2L), 1L), NULL))), .py_getitem(a, .py_slice(-2L, NULL, NULL)))), .py_tuple(.py_add(.py_add(.py_getitem(b, .py_slice(0L, 2L, NULL)), .py_getitem(b, .py_slice(.py_floordiv(.py_len(b), 2L), .py_add(.py_floordiv(.py_len(b), 2L), 1L), NULL))), .py_getitem(b, .py_slice(-2L, NULL, NULL))))))
}
y_coordinate_equality <- function(plot_data_func, evalf_func, system) {
  "Checks whether the y-coordinate value of the plotted\n    data point is equal to the value of the function at a\n    particular x."
  .tmpu_ <- plot_data_func(system); x <- .py_getitem(.tmpu_, 0L); y <- .py_getitem(.tmpu_, 1L)
  .tmpu_ <- p_trim_tuple(x, y); x <- .py_getitem(.tmpu_, 0L); y <- .py_getitem(.tmpu_, 1L)
  y_exp <- .py_tuple(.py_comp(x, function(.it) { x_i <- .it; list(evalf_func(system, x_i)) }))
  return(.py_all(.py_comp(.py_zip(y_exp, y), function(.it) { y_exp_i <- .py_getitem(.it, 0L); y_i <- .py_getitem(.it, 1L); list(.py_lt(Abs(.py_sub(y_exp_i, y_i)), .py_flt(1e-08))) })))
}
check_point_accuracy <- function(a, b) {
  return(.py_all(.py_comp(.py_zip(a, b), function(.it) { p_ <- .it; list(.py_call(.py_isclose, .py_splat(p_), rel_tol = .py_flt(0.1), abs_tol = .py_flt(1e-06))) })))
}
.py_try({
if (.py_truthy((!.py_truthy(numpy)))) {
  .py_try({
  NULL
  })
}
})
ramp_res_tester <- function(sys, num_points, expected_value, slope = 1L) {
  .tmpu_ <- .py_call(p_to_tuple, .py_splat(ramp_response_numerical_data(sys, slope = slope, adaptive = FALSE, n = num_points))); x <- .py_getitem(.tmpu_, 0L); y <- .py_getitem(.tmpu_, 1L)
  x_check <- check_point_accuracy(x, .py_getitem(expected_value, 0L))
  y_check <- check_point_accuracy(y, .py_getitem(expected_value, 1L))
  return((.py_truthy(x_check) && .py_truthy(y_check)))
}
.py_try({
exp1 <- .py_tup(.py_tup(.py_flt(0.0), .py_flt(2.0), .py_flt(4.0), .py_flt(6.0), .py_flt(8.0), .py_flt(10.0)), .py_tup(.py_flt(0.0), .py_flt(0.7324667795033895), .py_flt(1.9909720978650398), .py_flt(2.7956587704217783), .py_flt(3.9224897567931514), .py_flt(4.85022655284895)))
})
.py_try({
exp2 <- .py_tup(.py_tup(.py_flt(0.0), .py_flt(1.1111111111111112), .py_flt(2.2222222222222223), .py_flt(3.3333333333333335), .py_flt(4.444444444444445), .py_flt(5.555555555555555), .py_flt(6.666666666666667), .py_flt(7.777777777777779), .py_flt(8.88888888888889), .py_flt(10.0)), .py_tup(.py_flt(2.4360213402019326e-08), .py_flt(0.10175320182493253), .py_flt(0.33057612497658406), .py_flt(0.5967937263298935), .py_flt(0.8431511866718248), .py_flt(1.0398805391471613), .py_flt(1.1776043125035738), .py_flt(1.2600994825747305), .py_flt(1.2981042689274653), .py_flt(1.304684417610106)))
})
.py_try({
exp3 <- .py_tup(.py_tup(.py_flt(0.0), .py_flt(1.1111111111111112), .py_flt(2.2222222222222223), .py_flt(3.3333333333333335), .py_flt(4.444444444444445), .py_flt(5.555555555555555), .py_flt(6.666666666666667), .py_flt(7.777777777777779), .py_flt(8.88888888888889), .py_flt(10.0)), .py_tup(.py_flt(-3.9329040468771836e-08), .py_flt(0.34686634635794555), .py_flt(2.9998828170537903), .py_flt(12.33303690737476), .py_flt(40.993913948137795), .py_flt(127.84145222317912), .py_flt(391.41713691996), .py_flt(1192.0006858708389), .py_flt(3623.9808672503405), .py_flt(11011.728034546572)))
})
.py_try({
exp4 <- .py_tup(.py_tup(.py_flt(0.0), .py_flt(1.1111111111111112), .py_flt(2.2222222222222223), .py_flt(3.3333333333333335), .py_flt(4.444444444444445), .py_flt(5.555555555555555), .py_flt(6.666666666666667), .py_flt(7.777777777777779), .py_flt(8.88888888888889), .py_flt(10.0)), .py_tup(.py_flt(0.0), .py_flt(1.9051973784484078), .py_flt(30.483158055174524), .py_flt(154.32098765432104), .py_flt(487.7305288827924), .py_flt(1190.7483615302544), .py_flt(2469.1358024691367), .py_flt(4574.3789056546275), .py_flt(7803.688462124678), .py_flt(12500.0)))
})
.py_try({
exp5 <- .py_tup(.py_tup(.py_flt(0.0), .py_flt(1.1111111111111112), .py_flt(2.2222222222222223), .py_flt(3.3333333333333335), .py_flt(4.444444444444445), .py_flt(5.555555555555555), .py_flt(6.666666666666667), .py_flt(7.777777777777779), .py_flt(8.88888888888889), .py_flt(10.0)), .py_tup(.py_flt(0.0), .py_flt(3.8844361856975635), .py_flt(9.141792069209865), .py_flt(14.096349157657231), .py_flt(19.09783068994694), .py_flt(24.10179770390321), .py_flt(29.09907319114121), .py_flt(34.10040420185154), .py_flt(39.09983919254265), .py_flt(44.10006013058409)))
})
.py_try({
exp6 <- .py_tup(.py_tup(.py_flt(0.0), .py_flt(1.1111111111111112), .py_flt(2.2222222222222223), .py_flt(3.3333333333333335), .py_flt(4.444444444444445), .py_flt(5.555555555555555), .py_flt(6.666666666666667), .py_flt(7.777777777777779), .py_flt(8.88888888888889), .py_flt(10.0)), .py_tup(.py_flt(0.0), .py_flt(1.1111111111111112), .py_flt(2.2222222222222223), .py_flt(3.3333333333333335), .py_flt(4.444444444444445), .py_flt(5.555555555555555), .py_flt(6.666666666666667), .py_flt(7.777777777777779), .py_flt(8.88888888888889), .py_flt(10.0)))
})
.py_try({
.py_rec("control_plots:test_ramp_response:298", function() ramp_res_tester(tf1, 6L, exp1))
})
.py_try({
.py_rec("control_plots:test_ramp_response:299", function() ramp_res_tester(tf2, 10L, exp2, .py_flt(1.2)))
})
.py_try({
.py_rec("control_plots:test_ramp_response:300", function() ramp_res_tester(tf3, 10L, exp3, .py_flt(1.5)))
})
.py_try({
.py_rec("control_plots:test_ramp_response:301", function() ramp_res_tester(tf4, 10L, exp4, 3L))
})
.py_try({
.py_rec("control_plots:test_ramp_response:302", function() ramp_res_tester(tf5, 10L, exp5, 9L))
})
.py_try({
.py_rec("control_plots:test_ramp_response:303", function() ramp_res_tester(tf6, 10L, exp6))
}))---")
  expect_identical(got[["control_plots:test_ramp_response:298"]], "True")
  expect_identical(got[["control_plots:test_ramp_response:299"]], "True")
  expect_identical(got[["control_plots:test_ramp_response:300"]], "True")
  expect_identical(got[["control_plots:test_ramp_response:301"]], "True")
  expect_identical(got[["control_plots:test_ramp_response:302"]], "True")
  expect_identical(got[["control_plots:test_ramp_response:303"]], "True")
})

test_that("control_plots test_nyquist_plot_expr matches SymPy", {
  got <- .physmisc_battery_run(r"---(s <- Symbol("s")
p <- Symbol("p")
a <- Symbol("a")
.py_try({
matplotlib <- import_module("matplotlib", import_kwargs = .py_dict(list("fromlist"), list(list("pyplot"))), catch = .py_tup(RuntimeError))
})
.py_try({
numpy <- import_module("numpy")
})
.py_try({
tf1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(p, 2L), .py_mul(.py_flt(0.5), p)), 2L), p)
})
.py_try({
tf2 <- TransferFunction(p, .py_add(.py_add(.py_mul(6L, .py_pow(p, 2L)), .py_mul(3L, p)), 1L), p)
})
.py_try({
tf3 <- TransferFunction(p, .py_sub(.py_pow(p, 3L), 1L), p)
})
.py_try({
tf4 <- TransferFunction(10L, .py_pow(p, 3L), p)
})
.py_try({
tf5 <- TransferFunction(5L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(2L, s)), 10L), s)
})
.py_try({
tf6 <- TransferFunction(1L, 1L, s)
})
.py_try({
tf7 <- TransferFunction(.py_add(.py_add(.py_add(.py_mul(.py_mul(4L, s), 3L), .py_mul(9L, .py_pow(s, 2L))), .py_mul(.py_flt(0.1), s)), 11L), .py_add(.py_add(.py_mul(8L, .py_pow(s, 6L)), .py_mul(9L, .py_pow(s, 4L))), 11L), s)
})
.py_try({
tf8 <- TransferFunction(5L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_add(2L, sym_I), s)), 10L), s)
})
.py_try({
ser1 <- Series(tf4, TransferFunction(1L, .py_sub(p, 5L), p))
})
.py_try({
ser2 <- Series(tf3, TransferFunction(p, .py_add(p, 2L), p))
})
.py_try({
par1 <- Parallel(tf1, tf2)
})
p_to_tuple <- function(a, b) {
  return(.py_tup(.py_tuple(a), .py_tuple(b)))
}
p_trim_tuple <- function(a, b) {
  .tmpu_ <- p_to_tuple(a, b); a <- .py_getitem(.tmpu_, 0L); b <- .py_getitem(.tmpu_, 1L)
  return(.py_tup(.py_tuple(.py_add(.py_add(.py_getitem(a, .py_slice(0L, 2L, NULL)), .py_getitem(a, .py_slice(.py_floordiv(.py_len(a), 2L), .py_add(.py_floordiv(.py_len(a), 2L), 1L), NULL))), .py_getitem(a, .py_slice(-2L, NULL, NULL)))), .py_tuple(.py_add(.py_add(.py_getitem(b, .py_slice(0L, 2L, NULL)), .py_getitem(b, .py_slice(.py_floordiv(.py_len(b), 2L), .py_add(.py_floordiv(.py_len(b), 2L), 1L), NULL))), .py_getitem(b, .py_slice(-2L, NULL, NULL))))))
}
y_coordinate_equality <- function(plot_data_func, evalf_func, system) {
  "Checks whether the y-coordinate value of the plotted\n    data point is equal to the value of the function at a\n    particular x."
  .tmpu_ <- plot_data_func(system); x <- .py_getitem(.tmpu_, 0L); y <- .py_getitem(.tmpu_, 1L)
  .tmpu_ <- p_trim_tuple(x, y); x <- .py_getitem(.tmpu_, 0L); y <- .py_getitem(.tmpu_, 1L)
  y_exp <- .py_tuple(.py_comp(x, function(.it) { x_i <- .it; list(evalf_func(system, x_i)) }))
  return(.py_all(.py_comp(.py_zip(y_exp, y), function(.it) { y_exp_i <- .py_getitem(.it, 0L); y_i <- .py_getitem(.it, 1L); list(.py_lt(Abs(.py_sub(y_exp_i, y_i)), .py_flt(1e-08))) })))
}
check_point_accuracy <- function(a, b) {
  return(.py_all(.py_comp(.py_zip(a, b), function(.it) { p_ <- .it; list(.py_call(.py_isclose, .py_splat(p_), rel_tol = .py_flt(0.1), abs_tol = .py_flt(1e-06))) })))
}
.py_try({
.tmpu_ <- nyquist_plot_expr(tf1); r1 <- .py_getitem(.tmpu_, 0L); i1 <- .py_getitem(.tmpu_, 1L); w1 <- .py_getitem(.tmpu_, 2L)
})
.py_try({
.tmpu_ <- nyquist_plot_expr(tf2); r2 <- .py_getitem(.tmpu_, 0L); i2 <- .py_getitem(.tmpu_, 1L); w2 <- .py_getitem(.tmpu_, 2L)
})
.py_try({
.tmpu_ <- nyquist_plot_expr(tf3); r3 <- .py_getitem(.tmpu_, 0L); i3 <- .py_getitem(.tmpu_, 1L); w3 <- .py_getitem(.tmpu_, 2L)
})
.py_try({
.tmpu_ <- nyquist_plot_expr(tf4); r4 <- .py_getitem(.tmpu_, 0L); i4 <- .py_getitem(.tmpu_, 1L); w4 <- .py_getitem(.tmpu_, 2L)
})
.py_try({
.py_rec("control_plots:test_nyquist_plot_expr:311", function() r1)
})
.py_try({
.py_rec("control_plots:test_nyquist_plot_expr:312", function() i1)
})
.py_try({
.py_rec("control_plots:test_nyquist_plot_expr:313", function() r2)
})
.py_try({
.py_rec("control_plots:test_nyquist_plot_expr:314", function() i2)
})
.py_try({
.py_rec("control_plots:test_nyquist_plot_expr:315", function() r3)
})
.py_try({
.py_rec("control_plots:test_nyquist_plot_expr:316", function() i3)
})
.py_try({
.py_rec("control_plots:test_nyquist_plot_expr:317", function() r4)
})
.py_try({
.py_rec("control_plots:test_nyquist_plot_expr:318", function() i4)
}))---")
  expect_identical(got[["control_plots:test_nyquist_plot_expr:311"]], "(2 - _w^2)/(0.25*_w^2 + (2 - _w^2)^2)")
  expect_identical(got[["control_plots:test_nyquist_plot_expr:312"]], "-0.5*_w/(0.25*_w^2 + (2 - _w^2)^2)")
  expect_identical(got[["control_plots:test_nyquist_plot_expr:313"]], "3*_w^2/(9*_w^2 + (1 - 6*_w^2)^2)")
  expect_identical(got[["control_plots:test_nyquist_plot_expr:314"]], "_w*(1 - 6*_w^2)/(9*_w^2 + (1 - 6*_w^2)^2)")
  expect_identical(got[["control_plots:test_nyquist_plot_expr:315"]], "-_w^4/(_w^6 + 1)")
  expect_identical(got[["control_plots:test_nyquist_plot_expr:316"]], "-_w/(_w^6 + 1)")
  expect_identical(got[["control_plots:test_nyquist_plot_expr:317"]], "0")
  expect_identical(got[["control_plots:test_nyquist_plot_expr:318"]], "10/_w^3")
})

test_that("control_plots test_nichols_expr matches SymPy", {
  got <- .physmisc_battery_run(r"---(s <- Symbol("s")
p <- Symbol("p")
a <- Symbol("a")
.py_try({
matplotlib <- import_module("matplotlib", import_kwargs = .py_dict(list("fromlist"), list(list("pyplot"))), catch = .py_tup(RuntimeError))
})
.py_try({
numpy <- import_module("numpy")
})
.py_try({
tf1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(p, 2L), .py_mul(.py_flt(0.5), p)), 2L), p)
})
.py_try({
tf2 <- TransferFunction(p, .py_add(.py_add(.py_mul(6L, .py_pow(p, 2L)), .py_mul(3L, p)), 1L), p)
})
.py_try({
tf3 <- TransferFunction(p, .py_sub(.py_pow(p, 3L), 1L), p)
})
.py_try({
tf4 <- TransferFunction(10L, .py_pow(p, 3L), p)
})
.py_try({
tf5 <- TransferFunction(5L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(2L, s)), 10L), s)
})
.py_try({
tf6 <- TransferFunction(1L, 1L, s)
})
.py_try({
tf7 <- TransferFunction(.py_add(.py_add(.py_add(.py_mul(.py_mul(4L, s), 3L), .py_mul(9L, .py_pow(s, 2L))), .py_mul(.py_flt(0.1), s)), 11L), .py_add(.py_add(.py_mul(8L, .py_pow(s, 6L)), .py_mul(9L, .py_pow(s, 4L))), 11L), s)
})
.py_try({
tf8 <- TransferFunction(5L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_add(2L, sym_I), s)), 10L), s)
})
.py_try({
ser1 <- Series(tf4, TransferFunction(1L, .py_sub(p, 5L), p))
})
.py_try({
ser2 <- Series(tf3, TransferFunction(p, .py_add(p, 2L), p))
})
.py_try({
par1 <- Parallel(tf1, tf2)
})
p_to_tuple <- function(a, b) {
  return(.py_tup(.py_tuple(a), .py_tuple(b)))
}
p_trim_tuple <- function(a, b) {
  .tmpu_ <- p_to_tuple(a, b); a <- .py_getitem(.tmpu_, 0L); b <- .py_getitem(.tmpu_, 1L)
  return(.py_tup(.py_tuple(.py_add(.py_add(.py_getitem(a, .py_slice(0L, 2L, NULL)), .py_getitem(a, .py_slice(.py_floordiv(.py_len(a), 2L), .py_add(.py_floordiv(.py_len(a), 2L), 1L), NULL))), .py_getitem(a, .py_slice(-2L, NULL, NULL)))), .py_tuple(.py_add(.py_add(.py_getitem(b, .py_slice(0L, 2L, NULL)), .py_getitem(b, .py_slice(.py_floordiv(.py_len(b), 2L), .py_add(.py_floordiv(.py_len(b), 2L), 1L), NULL))), .py_getitem(b, .py_slice(-2L, NULL, NULL))))))
}
y_coordinate_equality <- function(plot_data_func, evalf_func, system) {
  "Checks whether the y-coordinate value of the plotted\n    data point is equal to the value of the function at a\n    particular x."
  .tmpu_ <- plot_data_func(system); x <- .py_getitem(.tmpu_, 0L); y <- .py_getitem(.tmpu_, 1L)
  .tmpu_ <- p_trim_tuple(x, y); x <- .py_getitem(.tmpu_, 0L); y <- .py_getitem(.tmpu_, 1L)
  y_exp <- .py_tuple(.py_comp(x, function(.it) { x_i <- .it; list(evalf_func(system, x_i)) }))
  return(.py_all(.py_comp(.py_zip(y_exp, y), function(.it) { y_exp_i <- .py_getitem(.it, 0L); y_i <- .py_getitem(.it, 1L); list(.py_lt(Abs(.py_sub(y_exp_i, y_i)), .py_flt(1e-08))) })))
}
check_point_accuracy <- function(a, b) {
  return(.py_all(.py_comp(.py_zip(a, b), function(.it) { p_ <- .it; list(.py_call(.py_isclose, .py_splat(p_), rel_tol = .py_flt(0.1), abs_tol = .py_flt(1e-06))) })))
}
.py_try({
.tmpu_ <- nichols_plot_expr(tf1); m1 <- .py_getitem(.tmpu_, 0L); p1 <- .py_getitem(.tmpu_, 1L); w1 <- .py_getitem(.tmpu_, 2L)
})
.py_try({
.tmpu_ <- nichols_plot_expr(tf2); m2 <- .py_getitem(.tmpu_, 0L); p2 <- .py_getitem(.tmpu_, 1L); w2 <- .py_getitem(.tmpu_, 2L)
})
.py_try({
.tmpu_ <- nichols_plot_expr(tf3); m3 <- .py_getitem(.tmpu_, 0L); p3 <- .py_getitem(.tmpu_, 1L); w3 <- .py_getitem(.tmpu_, 2L)
})
.py_try({
.tmpu_ <- nichols_plot_expr(tf4); m4 <- .py_getitem(.tmpu_, 0L); p4 <- .py_getitem(.tmpu_, 1L); w4 <- .py_getitem(.tmpu_, 2L)
})
.py_try({
.py_rec("control_plots:test_nichols_expr:326", function() m1)
})
.py_try({
.py_rec("control_plots:test_nichols_expr:327", function() p1)
})
.py_try({
.py_rec("control_plots:test_nichols_expr:328", function() m2)
})
.py_try({
.py_rec("control_plots:test_nichols_expr:329", function() p2)
})
.py_try({
.py_rec("control_plots:test_nichols_expr:330", function() m3)
}))---")
  expect_identical(got[["control_plots:test_nichols_expr:326"]], "20*log(1/sqrt(_w^4 - 3.75*_w^2 + 4))/log(10)")
  expect_identical(got[["control_plots:test_nichols_expr:327"]], "180*arg(1/(-_w^2 + 0.5*_w*I + 2))/pi")
  expect_identical(got[["control_plots:test_nichols_expr:328"]], "20*log(Abs(_w)/sqrt(36*_w^4 - 3*_w^2 + 1))/log(10)")
  expect_identical(got[["control_plots:test_nichols_expr:329"]], "180*arg(_w*I/(-6*_w^2 + 3*_w*I + 1))/pi")
  expect_identical(got[["control_plots:test_nichols_expr:330"]], "20*log(Abs(_w)/sqrt(_w^6 + 1))/log(10)")
})

test_that("control_plots test_issue_29530 matches SymPy", {
  got <- .physmisc_battery_run(r"---(s <- Symbol("s")
p <- Symbol("p")
a <- Symbol("a")
.py_try({
matplotlib <- import_module("matplotlib", import_kwargs = .py_dict(list("fromlist"), list(list("pyplot"))), catch = .py_tup(RuntimeError))
})
.py_try({
numpy <- import_module("numpy")
})
.py_try({
tf1 <- TransferFunction(1L, .py_add(.py_add(.py_pow(p, 2L), .py_mul(.py_flt(0.5), p)), 2L), p)
})
.py_try({
tf2 <- TransferFunction(p, .py_add(.py_add(.py_mul(6L, .py_pow(p, 2L)), .py_mul(3L, p)), 1L), p)
})
.py_try({
tf3 <- TransferFunction(p, .py_sub(.py_pow(p, 3L), 1L), p)
})
.py_try({
tf4 <- TransferFunction(10L, .py_pow(p, 3L), p)
})
.py_try({
tf5 <- TransferFunction(5L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(2L, s)), 10L), s)
})
.py_try({
tf6 <- TransferFunction(1L, 1L, s)
})
.py_try({
tf7 <- TransferFunction(.py_add(.py_add(.py_add(.py_mul(.py_mul(4L, s), 3L), .py_mul(9L, .py_pow(s, 2L))), .py_mul(.py_flt(0.1), s)), 11L), .py_add(.py_add(.py_mul(8L, .py_pow(s, 6L)), .py_mul(9L, .py_pow(s, 4L))), 11L), s)
})
.py_try({
tf8 <- TransferFunction(5L, .py_add(.py_add(.py_pow(s, 2L), .py_mul(.py_add(2L, sym_I), s)), 10L), s)
})
.py_try({
ser1 <- Series(tf4, TransferFunction(1L, .py_sub(p, 5L), p))
})
.py_try({
ser2 <- Series(tf3, TransferFunction(p, .py_add(p, 2L), p))
})
.py_try({
par1 <- Parallel(tf1, tf2)
})
p_to_tuple <- function(a, b) {
  return(.py_tup(.py_tuple(a), .py_tuple(b)))
}
p_trim_tuple <- function(a, b) {
  .tmpu_ <- p_to_tuple(a, b); a <- .py_getitem(.tmpu_, 0L); b <- .py_getitem(.tmpu_, 1L)
  return(.py_tup(.py_tuple(.py_add(.py_add(.py_getitem(a, .py_slice(0L, 2L, NULL)), .py_getitem(a, .py_slice(.py_floordiv(.py_len(a), 2L), .py_add(.py_floordiv(.py_len(a), 2L), 1L), NULL))), .py_getitem(a, .py_slice(-2L, NULL, NULL)))), .py_tuple(.py_add(.py_add(.py_getitem(b, .py_slice(0L, 2L, NULL)), .py_getitem(b, .py_slice(.py_floordiv(.py_len(b), 2L), .py_add(.py_floordiv(.py_len(b), 2L), 1L), NULL))), .py_getitem(b, .py_slice(-2L, NULL, NULL))))))
}
y_coordinate_equality <- function(plot_data_func, evalf_func, system) {
  "Checks whether the y-coordinate value of the plotted\n    data point is equal to the value of the function at a\n    particular x."
  .tmpu_ <- plot_data_func(system); x <- .py_getitem(.tmpu_, 0L); y <- .py_getitem(.tmpu_, 1L)
  .tmpu_ <- p_trim_tuple(x, y); x <- .py_getitem(.tmpu_, 0L); y <- .py_getitem(.tmpu_, 1L)
  y_exp <- .py_tuple(.py_comp(x, function(.it) { x_i <- .it; list(evalf_func(system, x_i)) }))
  return(.py_all(.py_comp(.py_zip(y_exp, y), function(.it) { y_exp_i <- .py_getitem(.it, 0L); y_i <- .py_getitem(.it, 1L); list(.py_lt(Abs(.py_sub(y_exp_i, y_i)), .py_flt(1e-08))) })))
}
check_point_accuracy <- function(a, b) {
  return(.py_all(.py_comp(.py_zip(a, b), function(.it) { p_ <- .it; list(.py_call(.py_isclose, .py_splat(p_), rel_tol = .py_flt(0.1), abs_tol = .py_flt(1e-06))) })))
}
.py_try({
if (.py_truthy((!.py_truthy(numpy)))) {
  .py_try({
  NULL
  })
}
})
.py_try({
num <- .py_np_array(list(.py_add(.py_flt(0.00348799648), complex(real = 0.0, imaginary = 0.0)), .py_sub(.py_flt(3.94714557), complex(real = 0.0, imaginary = 327.602726)), .py_sub(.py_flt(-12637021.6), complex(real = 0.0, imaginary = 862936.992)), .py_add(.py_flt(-28582292400.0), complex(real = 0.0, imaginary = 268995812000.0)), .py_add(.py_flt(3495380640000000.0), complex(real = 0.0, imaginary = 372323599000000.0)), .py_sub(.py_flt(2.17691576e+18), complex(real = 0.0, imaginary = 2.86131354e+19)), .py_sub(.py_flt(-1.47410187e+23), complex(real = 0.0, imaginary = 3.86894682e+21)), .py_add(.py_flt(1.76676107e+25), complex(real = 0.0, imaginary = 4.66860012e+26)), .py_sub(.py_flt(8.63781663e+29), complex(real = 0.0, imaginary = 1.02652261e+29)), .py_sub(.py_flt(-1.90093975e+32), complex(real = 0.0, imaginary = 8.42072953e+32)), .py_add(.py_flt(-3.29925768e+35), complex(real = 0.0, imaginary = 1.22585e+35))))
})
.py_try({
den <- .py_np_array(list(.py_sub(.py_flt(-9.19119158e+36), complex(real = 0.0, imaginary = 5.78774471e+36))))
})
.py_try({
num_poly <- .py_np_sum(.py_comp(.py_range(.py_len(num)), function(.it) { i <- .it; list(.py_mul(.py_getitem(num, i), .py_pow(s, .py_int(.py_sub(.py_sub(.py_len(num), 1L), i))))) }))
})
.py_try({
den_poly <- .py_np_sum(.py_comp(.py_range(.py_len(den)), function(.it) { i <- .it; list(.py_mul(.py_getitem(den, i), .py_pow(s, .py_int(.py_sub(.py_sub(.py_len(den), 1L), i))))) }))
})
.py_try({
sym_tf <- TransferFunction(num_poly, den_poly, s)
})
.py_try({
.tmpu_ <- bode_magnitude_numerical_data(sym_tf, freq_unit = "Hz", initial_exp = 0L, final_exp = 5L); freq <- .py_getitem(.tmpu_, 0L); mag <- .py_getitem(.tmpu_, 1L)
})
.py_try({
.py_rec("control_plots:test_issue_29530:353", function() .py_method(.py_np_isnan(mag), "any"))
}))---")
  expect_identical(got[["control_plots:test_issue_29530:353"]], "False")
})
