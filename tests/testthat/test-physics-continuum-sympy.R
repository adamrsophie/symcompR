skip_unless_full_tests()

# Expected values are the outputs of the same calls in SymPy (physics.continuum_mechanics);
# the cases come from its test files and run through the helpers in
# helper-physics-misc-battery.R.

test_that("beam test_Beam matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
x <- Symbol("x")
})
.py_try({
y <- Symbol("y")
})
.py_try({
.tmpu_ <- .py_symbols("R1, R2"); R1 <- .py_getitem(.tmpu_, 0L); R2 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
E <- Symbol("E")
})
.py_try({
E_1 <- Symbol("E_1")
})
.py_try({
I <- Symbol("I")
})
.py_try({
I_1 <- Symbol("I_1")
})
.py_try({
A <- Symbol("A")
})
.py_try({
b <- Beam(1L, E, I)
})
.py_try({
.py_rec("beam:test_Beam:29", function() .py_attr(b, "length"))
})
.py_try({
.py_rec("beam:test_Beam:30", function() .py_attr(b, "elastic_modulus"))
})
.py_try({
.py_rec("beam:test_Beam:31", function() .py_attr(b, "second_moment"))
})
.py_try({
.py_rec("beam:test_Beam:32", function() .py_attr(b, "variable"))
})
.py_try({
.py_setattr(b, "length", 4L)
})
.py_try({
.py_rec("beam:test_Beam:36", function() .py_attr(b, "length"))
})
.py_try({
.py_setattr(b, "elastic_modulus", E_1)
})
.py_try({
.py_rec("beam:test_Beam:40", function() .py_attr(b, "elastic_modulus"))
})
.py_try({
.py_setattr(b, "second_moment", I_1)
})
.py_try({
.py_rec("beam:test_Beam:44", function() .py_attr(b, "second_moment"))
})
.py_try({
.py_setattr(b, "variable", y)
})
.py_try({
.py_rec("beam:test_Beam:48", function() .py_attr(b, "variable"))
})
.py_try({
.py_setattr(b, "bc_deflection", list(.py_tup(0L, 2L)))
})
.py_try({
.py_setattr(b, "bc_slope", list(.py_tup(0L, 1L)))
})
.py_try({
.py_setattr(b, "bc_bending_moment", list(.py_tup(0L, 5L)))
})
.py_try({
.py_setattr(b, "bc_shear_force", list(.py_tup(2L, 1L)))
})
.py_try({
.py_rec("beam:test_Beam:55", function() .py_attr(b, "boundary_conditions"))
})
.py_try({
.py_setattr(b, "bc_shear_force", .py_list_mut(.py_attr(b, "bc_shear_force"), "extend", list(.py_tup(1L, 1L), .py_tup(2L, 3L))))
})
.py_try({
sf_bcs <- .py_attr(b, "bc_shear_force")
})
.py_try({
.py_rec("beam:test_Beam:61", function() sf_bcs)
})
.py_try({
.py_setattr(b, "bc_bending_moment", .py_list_mut(.py_attr(b, "bc_bending_moment"), "extend", list(.py_tup(1L, 3L), .py_tup(5L, 3L))))
})
.py_try({
bm_bcs <- .py_attr(b, "bc_bending_moment")
})
.py_try({
.py_rec("beam:test_Beam:66", function() bm_bcs)
})
.py_try({
.py_setattr(b, "bc_slope", .py_list_mut(.py_attr(b, "bc_slope"), "extend", list(.py_tup(4L, 3L), .py_tup(5L, 0L))))
})
.py_try({
s_bcs <- .py_attr(b, "bc_slope")
})
.py_try({
.py_rec("beam:test_Beam:71", function() s_bcs)
})
.py_try({
.py_setattr(b, "bc_deflection", .py_list_mut(.py_attr(b, "bc_deflection"), "extend", list(.py_tup(4L, 3L), .py_tup(5L, 0L))))
})
.py_try({
d_bcs <- .py_attr(b, "bc_deflection")
})
.py_try({
.py_rec("beam:test_Beam:76", function() d_bcs)
})
.py_try({
bcs_new <- .py_attr(b, "boundary_conditions")
})
.py_try({
.py_rec("beam:test_Beam:80", function() bcs_new)
})
.py_try({
b1 <- Beam(30L, E, I)
})
.py_try({
.py_method(b1, "apply_load", -8L, 0L, -1L)
})
.py_try({
.py_method(b1, "apply_load", R1, 10L, -1L)
})
.py_try({
.py_method(b1, "apply_load", R2, 30L, -1L)
})
.py_try({
.py_method(b1, "apply_load", 120L, 30L, -2L)
})
.py_try({
.py_setattr(b1, "bc_deflection", list(.py_tup(10L, 0L), .py_tup(30L, 0L)))
})
.py_try({
.py_method(b1, "solve_for_reaction_loads", R1, R2)
})
.py_try({
p <- .py_attr(b1, "reaction_loads")
})
.py_try({
q <- .py_dict(list(R1, R2), list(6L, 2L))
})
.py_try({
p <- .py_attr(b1, "load")
})
.py_try({
q <- .py_add(.py_add(.py_add(.py_mul(-8L, SingularityFunction(x, 0L, -1L)), .py_mul(6L, SingularityFunction(x, 10L, -1L))), .py_mul(120L, SingularityFunction(x, 30L, -2L))), .py_mul(2L, SingularityFunction(x, 30L, -1L)))
})
.py_try({
.py_rec("beam:test_Beam:103", function() p)
})
.py_try({
p <- .py_method(b1, "shear_force")
})
.py_try({
q <- .py_sub(.py_sub(.py_sub(.py_mul(8L, SingularityFunction(x, 0L, 0L)), .py_mul(6L, SingularityFunction(x, 10L, 0L))), .py_mul(120L, SingularityFunction(x, 30L, -1L))), .py_mul(2L, SingularityFunction(x, 30L, 0L)))
})
.py_try({
.py_rec("beam:test_Beam:109", function() p)
})
.py_try({
p <- .py_method(b1, "shear_stress")
})
.py_try({
q <- .py_div(.py_sub(.py_sub(.py_sub(.py_mul(8L, SingularityFunction(x, 0L, 0L)), .py_mul(6L, SingularityFunction(x, 10L, 0L))), .py_mul(120L, SingularityFunction(x, 30L, -1L))), .py_mul(2L, SingularityFunction(x, 30L, 0L))), A)
})
.py_try({
.py_rec("beam:test_Beam:116", function() p)
})
.py_try({
p <- .py_method(b1, "bending_moment")
})
.py_try({
q <- .py_sub(.py_sub(.py_sub(.py_mul(8L, SingularityFunction(x, 0L, 1L)), .py_mul(6L, SingularityFunction(x, 10L, 1L))), .py_mul(120L, SingularityFunction(x, 30L, 0L))), .py_mul(2L, SingularityFunction(x, 30L, 1L)))
})
.py_try({
.py_rec("beam:test_Beam:122", function() p)
})
.py_try({
p <- .py_method(b1, "slope")
})
.py_try({
q <- .py_add(.py_add(.py_add(.py_add(.py_mul(-4L, SingularityFunction(x, 0L, 2L)), .py_mul(3L, SingularityFunction(x, 10L, 2L))), .py_mul(120L, SingularityFunction(x, 30L, 1L))), SingularityFunction(x, 30L, 2L)), Rational(4000L, 3L))
})
.py_try({
.py_rec("beam:test_Beam:129", function() p)
})
.py_try({
p <- .py_method(b1, "deflection")
})
.py_try({
q <- .py_sub(.py_add(.py_add(.py_add(.py_sub(.py_mul(x, Rational(4000L, 3L)), .py_div(.py_mul(4L, SingularityFunction(x, 0L, 3L)), 3L)), SingularityFunction(x, 10L, 3L)), .py_mul(60L, SingularityFunction(x, 30L, 2L))), .py_div(SingularityFunction(x, 30L, 3L), 3L)), 12000L)
})
.py_try({
.py_rec("beam:test_Beam:136", function() p)
})
.py_try({
l <- Symbol("l")
})
.py_try({
w0 <- Symbol("w0")
})
.py_try({
w2 <- Symbol("w2")
})
.py_try({
a1 <- Symbol("a1")
})
.py_try({
c <- Symbol("c")
})
.py_try({
c1 <- Symbol("c1")
})
.py_try({
d <- Symbol("d")
})
.py_try({
e <- Symbol("e")
})
.py_try({
f <- Symbol("f")
})
.py_try({
b2 <- Beam(l, E, I)
})
.py_try({
.py_method(b2, "apply_load", w0, a1, 1L)
})
.py_try({
.py_method(b2, "apply_load", w2, c1, -1L)
})
.py_try({
.py_setattr(b2, "bc_deflection", list(.py_tup(c, d)))
})
.py_try({
.py_setattr(b2, "bc_slope", list(.py_tup(e, f)))
})
.py_try({
p <- .py_attr(b2, "load")
})
.py_try({
q <- .py_add(.py_mul(w0, SingularityFunction(x, a1, 1L)), .py_mul(w2, SingularityFunction(x, c1, -1L)))
})
.py_try({
.py_rec("beam:test_Beam:160", function() p)
})
.py_try({
p <- .py_method(b2, "shear_force")
})
.py_try({
q <- .py_sub(.py_div(.py_mul(.py_neg(w0), SingularityFunction(x, a1, 2L)), 2L), .py_mul(w2, SingularityFunction(x, c1, 0L)))
})
.py_try({
.py_rec("beam:test_Beam:166", function() p)
})
.py_try({
p <- .py_method(b2, "shear_stress")
})
.py_try({
q <- .py_div(.py_sub(.py_div(.py_mul(.py_neg(w0), SingularityFunction(x, a1, 2L)), 2L), .py_mul(w2, SingularityFunction(x, c1, 0L))), A)
})
.py_try({
.py_rec("beam:test_Beam:172", function() p)
})
.py_try({
p <- .py_method(b2, "bending_moment")
})
.py_try({
q <- .py_sub(.py_div(.py_mul(.py_neg(w0), SingularityFunction(x, a1, 3L)), 6L), .py_mul(w2, SingularityFunction(x, c1, 1L)))
})
.py_try({
.py_rec("beam:test_Beam:177", function() p)
})
.py_try({
p <- .py_method(b2, "slope")
})
.py_try({
q <- .py_add(.py_div(.py_add(.py_div(.py_mul(w0, SingularityFunction(x, a1, 4L)), 24L), .py_div(.py_mul(w2, SingularityFunction(x, c1, 2L)), 2L)), .py_mul(E, I)), .py_div(.py_sub(.py_sub(.py_mul(.py_mul(E, I), f), .py_div(.py_mul(w0, SingularityFunction(e, a1, 4L)), 24L)), .py_div(.py_mul(w2, SingularityFunction(e, c1, 2L)), 2L)), .py_mul(E, I)))
})
.py_try({
.py_rec("beam:test_Beam:182", function() expand(p))
})
.py_try({
p <- .py_method(b2, "deflection")
})
.py_try({
q <- .py_add(.py_add(.py_div(.py_mul(x, .py_sub(.py_sub(.py_mul(.py_mul(E, I), f), .py_div(.py_mul(w0, SingularityFunction(e, a1, 4L)), 24L)), .py_div(.py_mul(w2, SingularityFunction(e, c1, 2L)), 2L))), .py_mul(E, I)), .py_div(.py_add(.py_div(.py_mul(w0, SingularityFunction(x, a1, 5L)), 120L), .py_div(.py_mul(w2, SingularityFunction(x, c1, 3L)), 6L)), .py_mul(E, I))), .py_div(.py_sub(.py_sub(.py_add(.py_add(.py_mul(.py_mul(E, I), .py_add(.py_mul(.py_neg(c), f), d)), .py_div(.py_mul(.py_mul(c, w0), SingularityFunction(e, a1, 4L)), 24L)), .py_div(.py_mul(.py_mul(c, w2), SingularityFunction(e, c1, 2L)), 2L)), .py_div(.py_mul(w0, SingularityFunction(c, a1, 5L)), 120L)), .py_div(.py_mul(w2, SingularityFunction(c, c1, 3L)), 6L)), .py_mul(E, I)))
})
.py_try({
.py_rec("beam:test_Beam:194", function() simplify(.py_sub(p, q)))
})
.py_try({
b3 <- Beam(9L, E, I, 2L)
})
.py_try({
.py_method(b3, "apply_load", value = -2L, start = 2L, order = 2L, end = 3L)
})
.py_try({
.py_setattr(b3, "bc_slope", .py_list_mut(.py_attr(b3, "bc_slope"), "append", .py_tup(0L, 2L)))
})
.py_try({
C3 <- .py_symbols("C3")
})
.py_try({
C4 <- .py_symbols("C4")
})
.py_try({
p <- .py_attr(b3, "load")
})
.py_try({
q <- .py_add(.py_add(.py_add(.py_mul(-2L, SingularityFunction(x, 2L, 2L)), .py_mul(2L, SingularityFunction(x, 3L, 0L))), .py_mul(4L, SingularityFunction(x, 3L, 1L))), .py_mul(2L, SingularityFunction(x, 3L, 2L)))
})
.py_try({
.py_rec("beam:test_Beam:205", function() p)
})
.py_try({
p <- .py_method(b3, "shear_force")
})
.py_try({
q <- .py_sub(.py_sub(.py_sub(.py_div(.py_mul(2L, SingularityFunction(x, 2L, 3L)), 3L), .py_mul(2L, SingularityFunction(x, 3L, 1L))), .py_mul(2L, SingularityFunction(x, 3L, 2L))), .py_div(.py_mul(2L, SingularityFunction(x, 3L, 3L)), 3L))
})
.py_try({
.py_rec("beam:test_Beam:210", function() p)
})
.py_try({
p <- .py_method(b3, "shear_stress")
})
.py_try({
q <- .py_sub(.py_sub(.py_sub(.py_div(SingularityFunction(x, 2L, 3L), 3L), .py_mul(1L, SingularityFunction(x, 3L, 1L))), .py_mul(1L, SingularityFunction(x, 3L, 2L))), .py_div(.py_mul(1L, SingularityFunction(x, 3L, 3L)), 3L))
})
.py_try({
.py_rec("beam:test_Beam:215", function() p)
})
.py_try({
p <- .py_method(b3, "slope")
})
.py_try({
q <- .py_sub(2L, .py_div(.py_sub(.py_sub(.py_sub(.py_div(SingularityFunction(x, 2L, 5L), 30L), .py_div(SingularityFunction(x, 3L, 3L), 3L)), .py_div(SingularityFunction(x, 3L, 4L), 6L)), .py_div(SingularityFunction(x, 3L, 5L), 30L)), .py_mul(E, I)))
})
.py_try({
.py_rec("beam:test_Beam:220", function() p)
})
.py_try({
p <- .py_method(b3, "deflection")
})
.py_try({
q <- .py_sub(.py_mul(2L, x), .py_div(.py_sub(.py_sub(.py_sub(.py_div(SingularityFunction(x, 2L, 6L), 180L), .py_div(SingularityFunction(x, 3L, 4L), 12L)), .py_div(SingularityFunction(x, 3L, 5L), 30L)), .py_div(SingularityFunction(x, 3L, 6L), 180L)), .py_mul(E, I)))
})
.py_try({
.py_rec("beam:test_Beam:226", function() p)
})
.py_try({
b4 <- Beam(4L, E, I, 3L)
})
.py_try({
.py_method(b4, "apply_load", -3L, 0L, 0L, end = 3L)
})
.py_try({
p <- .py_attr(b4, "load")
})
.py_try({
q <- .py_add(.py_mul(-3L, SingularityFunction(x, 0L, 0L)), .py_mul(3L, SingularityFunction(x, 3L, 0L)))
})
.py_try({
.py_rec("beam:test_Beam:233", function() p)
})
.py_try({
p <- .py_method(b4, "shear_force")
})
.py_try({
q <- .py_sub(.py_mul(3L, SingularityFunction(x, 0L, 1L)), .py_mul(3L, SingularityFunction(x, 3L, 1L)))
})
.py_try({
.py_rec("beam:test_Beam:238", function() p)
})
.py_try({
p <- .py_method(b4, "shear_stress")
})
.py_try({
q <- .py_sub(SingularityFunction(x, 0L, 1L), SingularityFunction(x, 3L, 1L))
})
.py_try({
.py_rec("beam:test_Beam:242", function() p)
})
.py_try({
p <- .py_method(b4, "slope")
})
.py_try({
q <- .py_add(.py_div(.py_mul(-3L, SingularityFunction(x, 0L, 3L)), 6L), .py_div(.py_mul(3L, SingularityFunction(x, 3L, 3L)), 6L))
})
.py_try({
.py_rec("beam:test_Beam:246", function() p)
})
.py_try({
p <- .py_method(b4, "deflection")
})
.py_try({
q <- .py_add(.py_div(.py_mul(-3L, SingularityFunction(x, 0L, 4L)), 24L), .py_div(.py_mul(3L, SingularityFunction(x, 3L, 4L)), 24L))
})
.py_try({
.py_rec("beam:test_Beam:250", function() p)
})
.py_try({
.py_rec("beam:test_Beam:253", function() .py_method(b4, "apply_load", -3L, 0L, -1L, end = 3L))
})
.t254 <- function() {
  .py_setattr(b4, "variable", 1L)
}
.py_try({
.py_rec("beam:test_Beam:254", .t254)
}))---")
  expect_identical(got[["beam:test_Beam:29"]], "1")
  expect_identical(got[["beam:test_Beam:30"]], "E")
  expect_identical(got[["beam:test_Beam:31"]], "I")
  expect_identical(got[["beam:test_Beam:32"]], "x")
  expect_identical(got[["beam:test_Beam:36"]], "4")
  expect_identical(got[["beam:test_Beam:40"]], "E_1")
  expect_identical(got[["beam:test_Beam:44"]], "I_1")
  expect_identical(got[["beam:test_Beam:48"]], "y")
  expect_identical(got[["beam:test_Beam:55"]], "{'deflection': [(0, 2)], 'slope': [(0, 1)], 'bending_moment': [(0, 5)], 'shear_force': [(2, 1)]}")
  expect_identical(got[["beam:test_Beam:61"]], "[(2, 1), (1, 1), (2, 3)]")
  expect_identical(got[["beam:test_Beam:66"]], "[(0, 5), (1, 3), (5, 3)]")
  expect_identical(got[["beam:test_Beam:71"]], "[(0, 1), (4, 3), (5, 0)]")
  expect_identical(got[["beam:test_Beam:76"]], "[(0, 2), (4, 3), (5, 0)]")
  expect_identical(got[["beam:test_Beam:80"]], "{'deflection': [(0, 2), (4, 3), (5, 0)], 'slope': [(0, 1), (4, 3), (5, 0)], 'bending_moment': [(0, 5), (1, 3), (5, 3)], 'shear_force': [(2, 1), (1, 1), (2, 3)]}")
  expect_identical(got[["beam:test_Beam:103"]], "-8*SingularityFunction(x, 0, -1) + 6*SingularityFunction(x, 10, -1) + 120*SingularityFunction(x, 30, -2) + 2*SingularityFunction(x, 30, -1)")
  expect_identical(got[["beam:test_Beam:109"]], "8*SingularityFunction(x, 0, 0) - 6*SingularityFunction(x, 10, 0) - 120*SingularityFunction(x, 30, -1) - 2*SingularityFunction(x, 30, 0)")
  expect_identical(got[["beam:test_Beam:116"]], "(8*SingularityFunction(x, 0, 0) - 6*SingularityFunction(x, 10, 0) - 120*SingularityFunction(x, 30, -1) - 2*SingularityFunction(x, 30, 0))/A")
  expect_identical(got[["beam:test_Beam:122"]], "8*SingularityFunction(x, 0, 1) - 6*SingularityFunction(x, 10, 1) - 120*SingularityFunction(x, 30, 0) - 2*SingularityFunction(x, 30, 1)")
  expect_identical(got[["beam:test_Beam:129"]], "(-4*SingularityFunction(x, 0, 2) + 3*SingularityFunction(x, 10, 2) + 120*SingularityFunction(x, 30, 1) + SingularityFunction(x, 30, 2) + 4000/3)/(E*I)")
  expect_identical(got[["beam:test_Beam:136"]], "(4000*x/3 - 4*SingularityFunction(x, 0, 3)/3 + SingularityFunction(x, 10, 3) + 60*SingularityFunction(x, 30, 2) + SingularityFunction(x, 30, 3)/3 - 12000)/(E*I)")
  expect_identical(got[["beam:test_Beam:160"]], "w0*SingularityFunction(x, a1, 1) + w2*SingularityFunction(x, c1, -1)")
  expect_identical(got[["beam:test_Beam:166"]], "-w0*SingularityFunction(x, a1, 2)/2 - w2*SingularityFunction(x, c1, 0)")
  expect_identical(got[["beam:test_Beam:172"]], "(-w0*SingularityFunction(x, a1, 2)/2 - w2*SingularityFunction(x, c1, 0))/A")
  expect_identical(got[["beam:test_Beam:177"]], "-w0*SingularityFunction(x, a1, 3)/6 - w2*SingularityFunction(x, c1, 1)")
  expect_identical(got[["beam:test_Beam:182"]], "f - w0*SingularityFunction(e, a1, 4)/(24*E*I) + w0*SingularityFunction(x, a1, 4)/(24*E*I) - w2*SingularityFunction(e, c1, 2)/(2*E*I) + w2*SingularityFunction(x, c1, 2)/(2*E*I)")
  expect_identical(got[["beam:test_Beam:194"]], "0")
  expect_identical(got[["beam:test_Beam:205"]], "-2*SingularityFunction(x, 2, 2) + 2*SingularityFunction(x, 3, 0) + 4*SingularityFunction(x, 3, 1) + 2*SingularityFunction(x, 3, 2)")
  expect_identical(got[["beam:test_Beam:210"]], "2*SingularityFunction(x, 2, 3)/3 - 2*SingularityFunction(x, 3, 1) - 2*SingularityFunction(x, 3, 2) - 2*SingularityFunction(x, 3, 3)/3")
  expect_identical(got[["beam:test_Beam:215"]], "SingularityFunction(x, 2, 3)/3 - SingularityFunction(x, 3, 1) - SingularityFunction(x, 3, 2) - SingularityFunction(x, 3, 3)/3")
  expect_identical(got[["beam:test_Beam:220"]], "2 - (SingularityFunction(x, 2, 5)/30 - SingularityFunction(x, 3, 3)/3 - SingularityFunction(x, 3, 4)/6 - SingularityFunction(x, 3, 5)/30)/(E*I)")
  expect_identical(got[["beam:test_Beam:226"]], "C4 + 2*x - (SingularityFunction(x, 2, 6)/180 - SingularityFunction(x, 3, 4)/12 - SingularityFunction(x, 3, 5)/30 - SingularityFunction(x, 3, 6)/180)/(E*I)")
  expect_identical(got[["beam:test_Beam:233"]], "-3*SingularityFunction(x, 0, 0) + 3*SingularityFunction(x, 3, 0)")
  expect_identical(got[["beam:test_Beam:238"]], "3*SingularityFunction(x, 0, 1) - 3*SingularityFunction(x, 3, 1)")
  expect_identical(got[["beam:test_Beam:242"]], "SingularityFunction(x, 0, 1) - SingularityFunction(x, 3, 1)")
  expect_identical(got[["beam:test_Beam:246"]], "C3 + (-SingularityFunction(x, 0, 3)/2 + SingularityFunction(x, 3, 3)/2)/(E*I)")
  expect_identical(got[["beam:test_Beam:250"]], "C3*x + C4 + (-SingularityFunction(x, 0, 4)/8 + SingularityFunction(x, 3, 4)/8)/(E*I)")
  expect_identical(got[["beam:test_Beam:253"]], "ERR ValueError: If 'end' is provided the 'order' of the load cannot be negative, i.e. 'end' is only valid for distributed loads.")
  expect_identical(got[["beam:test_Beam:254"]], "ERR TypeError: The variable should be a Symbol object.")
})

test_that("beam test_insufficient_bconditions matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
x <- Symbol("x")
})
.py_try({
y <- Symbol("y")
})
.py_try({
.tmpu_ <- .py_symbols("R1, R2"); R1 <- .py_getitem(.tmpu_, 0L); R2 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
L <- .py_symbols("L", positive = TRUE)
})
.py_try({
.tmpu_ <- .py_symbols("E I P a3 a4"); E <- .py_getitem(.tmpu_, 0L); I <- .py_getitem(.tmpu_, 1L); P <- .py_getitem(.tmpu_, 2L); a3 <- .py_getitem(.tmpu_, 3L); a4 <- .py_getitem(.tmpu_, 4L)
})
.py_try({
b <- Beam(L, E, I, base_char = "a")
})
.py_try({
.py_method(b, "apply_load", R2, L, -1L)
})
.py_try({
.py_method(b, "apply_load", R1, 0L, -1L)
})
.py_try({
.py_method(b, "apply_load", .py_neg(P), .py_div(L, 2L), -1L)
})
.py_try({
.py_method(b, "solve_for_reaction_loads", R1, R2)
})
.py_try({
p <- .py_method(b, "slope")
})
.py_try({
q <- .py_add(.py_sub(.py_div(.py_mul(P, SingularityFunction(x, 0L, 2L)), 4L), .py_div(.py_mul(P, SingularityFunction(x, .py_div(L, 2L), 2L)), 2L)), .py_div(.py_mul(P, SingularityFunction(x, L, 2L)), 4L))
})
.py_try({
.py_rec("beam:test_insufficient_bconditions:272", function() p)
})
.py_try({
p <- .py_method(b, "deflection")
})
.py_try({
q <- .py_add(.py_sub(.py_div(.py_mul(P, SingularityFunction(x, 0L, 3L)), 12L), .py_div(.py_mul(P, SingularityFunction(x, .py_div(L, 2L), 3L)), 6L)), .py_div(.py_mul(P, SingularityFunction(x, L, 3L)), 12L))
})
.py_try({
.py_rec("beam:test_insufficient_bconditions:276", function() p)
})
.py_try({
.py_setattr(b, "bc_deflection", list(.py_tup(0L, 0L)))
})
.py_try({
p <- .py_method(b, "deflection")
})
.py_try({
q <- .py_add(.py_sub(.py_add(.py_mul(a3, x), .py_div(.py_mul(P, SingularityFunction(x, 0L, 3L)), 12L)), .py_div(.py_mul(P, SingularityFunction(x, .py_div(L, 2L), 3L)), 6L)), .py_div(.py_mul(P, SingularityFunction(x, L, 3L)), 12L))
})
.py_try({
.py_rec("beam:test_insufficient_bconditions:281", function() p)
})
.py_try({
.py_setattr(b, "bc_deflection", list(.py_tup(0L, 0L), .py_tup(L, 0L)))
})
.py_try({
p <- .py_method(b, "deflection")
})
.py_try({
q <- .py_add(.py_sub(.py_add(.py_div(.py_mul(.py_mul(.py_neg(.py_pow(L, 2L)), P), x), 16L), .py_div(.py_mul(P, SingularityFunction(x, 0L, 3L)), 12L)), .py_div(.py_mul(P, SingularityFunction(x, .py_div(L, 2L), 3L)), 6L)), .py_div(.py_mul(P, SingularityFunction(x, L, 3L)), 12L))
})
.py_try({
.py_rec("beam:test_insufficient_bconditions:286", function() p)
}))---")
  expect_identical(got[["beam:test_insufficient_bconditions:272"]], "a3 + (P*SingularityFunction(x, 0, 2)/4 - P*SingularityFunction(x, L/2, 2)/2 + P*SingularityFunction(x, L, 2)/4)/(E*I)")
  expect_identical(got[["beam:test_insufficient_bconditions:276"]], "a3*x + a4 + (P*SingularityFunction(x, 0, 3)/12 - P*SingularityFunction(x, L/2, 3)/6 + P*SingularityFunction(x, L, 3)/12)/(E*I)")
  expect_identical(got[["beam:test_insufficient_bconditions:281"]], "(P*SingularityFunction(x, 0, 3)/12 - P*SingularityFunction(x, L/2, 3)/6 + P*SingularityFunction(x, L, 3)/12 + a3*x)/(E*I)")
  expect_identical(got[["beam:test_insufficient_bconditions:286"]], "(-L^2*P*x/16 + P*SingularityFunction(x, 0, 3)/12 - P*SingularityFunction(x, L/2, 3)/6 + P*SingularityFunction(x, L, 3)/12)/(E*I)")
})

test_that("beam test_beam_units matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
x <- Symbol("x")
})
.py_try({
y <- Symbol("y")
})
.py_try({
.tmpu_ <- .py_symbols("R1, R2"); R1 <- .py_getitem(.tmpu_, 0L); R2 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
E <- Symbol("E")
})
.py_try({
I <- Symbol("I")
})
.py_try({
.tmpu_ <- .py_symbols("R1, R2"); R1 <- .py_getitem(.tmpu_, 0L); R2 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
kN <- .py_mul(.py_unit("kilo"), .py_unit("newton"))
})
.py_try({
gN <- .py_mul(.py_unit("giga"), .py_unit("newton"))
})
.py_try({
b <- Beam(.py_mul(8L, .py_unit("meter")), .py_div(.py_mul(200L, gN), .py_pow(.py_unit("meter"), 2L)), .py_mul(.py_mul(400L, 1000000L), .py_pow(.py_mul(.py_unit("milli"), .py_unit("meter")), 4L)))
})
.py_try({
.py_method(b, "apply_load", .py_mul(5L, kN), .py_mul(2L, .py_unit("meter")), -1L)
})
.py_try({
.py_method(b, "apply_load", R1, .py_mul(0L, .py_unit("meter")), -1L)
})
.py_try({
.py_method(b, "apply_load", R2, .py_mul(8L, .py_unit("meter")), -1L)
})
.py_try({
.py_method(b, "apply_load", .py_div(.py_mul(10L, kN), .py_unit("meter")), .py_mul(4L, .py_unit("meter")), 0L, end = .py_mul(8L, .py_unit("meter")))
})
.py_try({
.py_setattr(b, "bc_deflection", list(.py_tup(.py_mul(0L, .py_unit("meter")), .py_mul(0L, .py_unit("meter"))), .py_tup(.py_mul(8L, .py_unit("meter")), .py_mul(0L, .py_unit("meter")))))
})
.py_try({
.py_method(b, "solve_for_reaction_loads", R1, R2)
})
.py_try({
b <- Beam(.py_mul(3L, .py_unit("meter")), .py_div(.py_mul(E, .py_unit("newton")), .py_pow(.py_unit("meter"), 2L)), .py_mul(I, .py_pow(.py_unit("meter"), 4L)))
})
.py_try({
.py_method(b, "apply_load", .py_mul(8L, kN), .py_mul(1L, .py_unit("meter")), -1L)
})
.py_try({
.py_method(b, "apply_load", R1, .py_mul(0L, .py_unit("meter")), -1L)
})
.py_try({
.py_method(b, "apply_load", R2, .py_mul(3L, .py_unit("meter")), -1L)
})
.py_try({
.py_method(b, "apply_load", .py_mul(.py_mul(12L, kN), .py_unit("meter")), .py_mul(2L, .py_unit("meter")), -2L)
})
.py_try({
.py_setattr(b, "bc_deflection", list(.py_tup(.py_mul(0L, .py_unit("meter")), .py_mul(0L, .py_unit("meter"))), .py_tup(.py_mul(3L, .py_unit("meter")), .py_mul(0L, .py_unit("meter")))))
})
.py_try({
.py_method(b, "solve_for_reaction_loads", R1, R2)
})
.py_try({
.py_rec("beam:test_beam_units:337", function() .py_method(.py_method(b, "deflection"), "subs", x, .py_mul(1L, .py_unit("meter"))))
}))---")
  expect_identical(got[["beam:test_beam_units:337"]], "62000*meter/(9*E*I)")
})

test_that("beam test_composite_beam matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
x <- Symbol("x")
})
.py_try({
y <- Symbol("y")
})
.py_try({
.tmpu_ <- .py_symbols("R1, R2"); R1 <- .py_getitem(.tmpu_, 0L); R2 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
E <- Symbol("E")
})
.py_try({
I <- Symbol("I")
})
.py_try({
b1 <- Beam(2L, E, .py_mul(.py_flt(1.5), I))
})
.py_try({
b2 <- Beam(2L, E, I)
})
.py_try({
b <- .py_method(b1, "join", b2, "fixed")
})
.py_try({
.py_method(b, "apply_load", -20L, 0L, -1L)
})
.py_try({
.py_method(b, "apply_load", 80L, 0L, -2L)
})
.py_try({
.py_method(b, "apply_load", 20L, 4L, -1L)
})
.py_try({
.py_setattr(b, "bc_slope", list(.py_tup(0L, 0L)))
})
.py_try({
.py_setattr(b, "bc_deflection", list(.py_tup(0L, 0L)))
})
.py_try({
.py_rec("beam:test_composite_beam:384", function() .py_attr(b, "length"))
})
.py_try({
.py_rec("beam:test_composite_beam:385", function() .py_attr(b, "second_moment"))
})
.py_try({
.py_rec("beam:test_composite_beam:386", function() .py_method(.py_method(b, "slope"), "subs", x, 4L))
})
.py_try({
.py_rec("beam:test_composite_beam:387", function() .py_method(.py_method(b, "slope"), "subs", x, 2L))
})
.py_try({
.py_rec("beam:test_composite_beam:388", function() .py_int(.py_getitem(.py_attr(.py_method(.py_method(b, "deflection"), "subs", x, 4L), "args"), 0L)))
})
.py_try({
l <- .py_symbols("l", positive = TRUE)
})
.py_try({
.tmpu_ <- .py_symbols("R1 M1 R2 R3 P"); R1 <- .py_getitem(.tmpu_, 0L); M1 <- .py_getitem(.tmpu_, 1L); R2 <- .py_getitem(.tmpu_, 2L); R3 <- .py_getitem(.tmpu_, 3L); P <- .py_getitem(.tmpu_, 4L)
})
.py_try({
b1 <- Beam(.py_mul(2L, l), E, I)
})
.py_try({
b2 <- Beam(.py_mul(2L, l), E, I)
})
.py_try({
b <- .py_method(b1, "join", b2, "hinge")
})
.py_try({
.py_method(b, "apply_load", M1, 0L, -2L)
})
.py_try({
.py_method(b, "apply_load", R1, 0L, -1L)
})
.py_try({
.py_method(b, "apply_load", R2, l, -1L)
})
.py_try({
.py_method(b, "apply_load", R3, .py_mul(4L, l), -1L)
})
.py_try({
.py_method(b, "apply_load", P, .py_mul(3L, l), -1L)
})
.py_try({
.py_setattr(b, "bc_slope", list(.py_tup(0L, 0L)))
})
.py_try({
.py_setattr(b, "bc_deflection", list(.py_tup(0L, 0L), .py_tup(l, 0L), .py_tup(.py_mul(4L, l), 0L)))
})
.py_try({
.py_method(b, "solve_for_reaction_loads", M1, R1, R2, R3)
})
.py_try({
.py_rec("beam:test_composite_beam:404", function() .py_method(.py_method(b, "slope"), "subs", x, .py_mul(3L, l)))
})
.py_try({
.py_rec("beam:test_composite_beam:405", function() .py_method(.py_method(b, "deflection"), "subs", x, .py_mul(2L, l)))
})
.py_try({
.py_rec("beam:test_composite_beam:406", function() .py_method(.py_method(b, "deflection"), "subs", x, .py_mul(3L, l)))
})
.py_try({
b1 <- Beam(2L, 500L, 10L)
})
.py_try({
b2 <- Beam(2L, 500L, 10L)
})
.py_try({
b <- .py_method(b1, "join", b2, "fixed")
})
.py_try({
.py_method(b, "apply_load", M1, 0L, -2L)
})
.py_try({
.py_method(b, "apply_load", R1, 0L, -1L)
})
.py_try({
.py_method(b, "apply_load", R2, 1L, -1L)
})
.py_try({
.py_method(b, "apply_load", R3, 4L, -1L)
})
.py_try({
.py_method(b, "apply_load", 10L, 3L, -1L)
})
.py_try({
.py_setattr(b, "bc_slope", list(.py_tup(0L, 0L)))
})
.py_try({
.py_setattr(b, "bc_deflection", list(.py_tup(0L, 0L), .py_tup(1L, 0L), .py_tup(4L, 0L)))
})
.py_try({
.py_method(b, "solve_for_reaction_loads", M1, R1, R2, R3)
})
.py_try({
.py_rec("beam:test_composite_beam:420", function() .py_method(b, "slope"))
})
.py_try({
.py_rec("beam:test_composite_beam:423", function() .py_method(b, "deflection"))
})
.py_try({
E <- Symbol("E")
})
.py_try({
I <- Symbol("I")
})
.py_try({
b <- Beam(6L, E, I)
})
.py_try({
c <- Beam(6L, E, I)
})
.t431 <- function() {
  .py_method(b, "join", c, "hige")
}
.py_try({
.py_rec("beam:test_composite_beam:431", .t431)
}))---")
  expect_identical(got[["beam:test_composite_beam:384"]], "4")
  expect_identical(got[["beam:test_composite_beam:385"]], "Piecewise((1.5*I, x <= 2), (I, x <= 4))")
  expect_identical(got[["beam:test_composite_beam:386"]], "120.0/(E*I)")
  expect_identical(got[["beam:test_composite_beam:387"]], "80.0/(E*I)")
  expect_identical(got[["beam:test_composite_beam:388"]], "-302")
  expect_identical(got[["beam:test_composite_beam:404"]], "-7*P*l^2/(48*E*I)")
  expect_identical(got[["beam:test_composite_beam:405"]], "7*P*l^3/(24*E*I)")
  expect_identical(got[["beam:test_composite_beam:406"]], "5*P*l^3/(16*E*I)")
  expect_identical(got[["beam:test_composite_beam:420"]], "-2*SingularityFunction(x, 0, 1)/5625 + SingularityFunction(x, 0, 2)/1875 - 133*SingularityFunction(x, 1, 2)/135000 + SingularityFunction(x, 3, 2)/1000 - 37*SingularityFunction(x, 4, 2)/67500")
  expect_identical(got[["beam:test_composite_beam:423"]], "-SingularityFunction(x, 0, 2)/5625 + SingularityFunction(x, 0, 3)/5625 - 133*SingularityFunction(x, 1, 3)/405000 + SingularityFunction(x, 3, 3)/3000 - 37*SingularityFunction(x, 4, 3)/202500")
  expect_identical(got[["beam:test_composite_beam:431"]], "ERR ValueError: Invalid joining method. Choose from 'fixed' or 'hinge'.")
})

test_that("beam test_point_cflexure_single matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
x <- Symbol("x")
})
.py_try({
y <- Symbol("y")
})
.py_try({
.tmpu_ <- .py_symbols("R1, R2"); R1 <- .py_getitem(.tmpu_, 0L); R2 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
E <- Symbol("E")
})
.py_try({
I <- Symbol("I")
})
.py_try({
b <- Beam(10L, E, I)
})
.py_try({
.py_method(b, "apply_load", -4L, 0L, -1L)
})
.py_try({
.py_method(b, "apply_load", -46L, 6L, -1L)
})
.py_try({
.py_method(b, "apply_load", 10L, 2L, -1L)
})
.py_try({
.py_method(b, "apply_load", 20L, 4L, -1L)
})
.py_try({
.py_method(b, "apply_load", 3L, 6L, 0L)
})
.py_try({
.py_rec("beam:test_point_cflexure_single:444", function() .py_method(b, "point_cflexure"))
}))---")
  expect_identical(got[["beam:test_point_cflexure_single:444"]], "[10/3]")
})

test_that("beam test_point_cflexure_multiple matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
x <- Symbol("x")
})
.py_try({
y <- Symbol("y")
})
.py_try({
.tmpu_ <- .py_symbols("R1, R2"); R1 <- .py_getitem(.tmpu_, 0L); R2 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
E <- Symbol("E")
})
.py_try({
I <- Symbol("I")
})
.py_try({
b <- Beam(15L, E, I)
})
.py_try({
r0 <- .py_method(b, "apply_support", 0L, type = "pin")
})
.py_try({
r10 <- .py_method(b, "apply_support", 10L, type = "pin")
})
.py_try({
.tmpu_ <- .py_method(b, "apply_support", 15L, type = "fixed"); r15 <- .py_getitem(.tmpu_, 0L); m15 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
.py_method(b, "apply_rotation_hinge", 12L)
})
.py_try({
.py_method(b, "apply_load", -10L, 5L, -1L)
})
.py_try({
.py_method(b, "apply_load", -5L, 10L, 0L, 15L)
})
.py_try({
.py_method(b, "solve_for_reaction_loads", r0, r10, r15, m15)
})
.py_try({
.py_rec("beam:test_point_cflexure_multiple:459", function() .py_method(b, "point_cflexure"))
}))---")
  expect_identical(got[["beam:test_point_cflexure_multiple:459"]], "[1200/163, 12, 163/12]")
})

test_that("beam test_point_cflexure_slow_cases matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
x <- Symbol("x")
})
.py_try({
y <- Symbol("y")
})
.py_try({
.tmpu_ <- .py_symbols("R1, R2"); R1 <- .py_getitem(.tmpu_, 0L); R2 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
E <- Symbol("E")
})
.py_try({
I <- Symbol("I")
})
.py_try({
b <- Beam(10L, E, I)
})
.py_try({
.tmpu_ <- .py_method(b, "apply_support", 0L, type = "fixed"); r0 <- .py_getitem(.tmpu_, 0L); m0 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
r4 <- .py_method(b, "apply_support", 4L, type = "pin")
})
.py_try({
.tmpu_ <- .py_method(b, "apply_support", 7L, type = "fixed"); r7 <- .py_getitem(.tmpu_, 0L); m7 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
.py_method(b, "apply_rotation_hinge", 2L)
})
.py_try({
.py_method(b, "apply_rotation_hinge", 5L)
})
.py_try({
.py_method(b, "apply_load", -5L, 0L, 0L, 2L)
})
.py_try({
.py_method(b, "apply_load", 7L, 3L, -1L)
})
.py_try({
.py_method(b, "solve_for_reaction_loads", r0, m0, r4, r7, m7)
})
.py_try({
.py_rec("beam:test_point_cflexure_slow_cases:476", function() .py_method(b, "point_cflexure"))
})
.py_try({
E <- .py_symbols("E")
})
.py_try({
I <- .py_symbols("I")
})
.py_try({
P <- .py_symbols("P", positive = TRUE)
})
.py_try({
b <- Beam(15L, E, I)
})
.py_try({
r0 <- .py_method(b, "apply_support", 0L, type = "pin")
})
.py_try({
r10 <- .py_method(b, "apply_support", 10L, type = "pin")
})
.py_try({
.tmpu_ <- .py_method(b, "apply_support", 15L, type = "fixed"); r15 <- .py_getitem(.tmpu_, 0L); m15 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
.py_method(b, "apply_rotation_hinge", 12L)
})
.py_try({
.py_method(b, "apply_load", .py_neg(P), 5L, -1L)
})
.py_try({
.py_method(b, "solve_for_reaction_loads", r0, r10, r15, m15)
})
.py_try({
.py_rec("beam:test_point_cflexure_slow_cases:489", function() .py_method(b, "point_cflexure"))
})
.py_try({
E <- Symbol("E")
})
.py_try({
I <- Symbol("I")
})
.py_try({
.tmpu_ <- .py_symbols("M1, M2, R1,R2"); M1 <- .py_getitem(.tmpu_, 0L); M2 <- .py_getitem(.tmpu_, 1L); R1 <- .py_getitem(.tmpu_, 2L); R2 <- .py_getitem(.tmpu_, 3L)
})
.py_try({
F <- Symbol("F", positive = TRUE)
})
.py_try({
b5 <- Beam(10L, E, I)
})
.py_try({
.py_setattr(b5, "bc_deflection", list(.py_tup(0L, 0L), .py_tup(10L, 0L)))
})
.py_try({
.py_setattr(b5, "bc_slope", list(.py_tup(0L, 0L), .py_tup(10L, 0L)))
})
.py_try({
.py_method(b5, "apply_load", R1, 0L, -1L)
})
.py_try({
.py_method(b5, "apply_load", M1, 0L, -2L)
})
.py_try({
.py_method(b5, "apply_load", R2, 10L, -1L)
})
.py_try({
.py_method(b5, "apply_load", M2, 10L, -2L)
})
.py_try({
.py_method(b5, "apply_load", .py_neg(F), 5L, -1L)
})
.py_try({
.py_method(b5, "solve_for_reaction_loads", R1, R2, M1, M2)
})
.py_try({
.py_rec("beam:test_point_cflexure_slow_cases:505", function() .py_method(b5, "point_cflexure"))
})
.py_try({
E <- Symbol("E")
})
.py_try({
I <- Symbol("I")
})
.py_try({
b <- Beam(10L, E, I)
})
.py_try({
.tmpu_ <- .py_method(b, "apply_support", 0L, type = "fixed"); r0 <- .py_getitem(.tmpu_, 0L); m0 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
r4 <- .py_method(b, "apply_support", 4L, type = "pin")
})
.py_try({
.tmpu_ <- .py_method(b, "apply_support", 7L, type = "fixed"); r7 <- .py_getitem(.tmpu_, 0L); m7 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
.py_method(b, "apply_rotation_hinge", 2L)
})
.py_try({
.py_method(b, "apply_rotation_hinge", 5L)
})
.py_try({
.py_method(b, "apply_load", -5L, 0L, 0L, 2L)
})
.py_try({
.py_method(b, "apply_load", 7L, 3L, -1L)
})
.py_try({
.py_method(b, "solve_for_reaction_loads", r0, m0, r4, r7, m7)
})
.py_try({
.py_rec("beam:test_point_cflexure_slow_cases:519", function() .py_method(b, "point_cflexure"))
})
.py_try({
E <- .py_symbols("E")
})
.py_try({
I <- .py_symbols("I")
})
.py_try({
P <- .py_symbols("P", positive = TRUE)
})
.py_try({
b <- Beam(15L, E, I)
})
.py_try({
r0 <- .py_method(b, "apply_support", 0L, type = "pin")
})
.py_try({
r10 <- .py_method(b, "apply_support", 10L, type = "pin")
})
.py_try({
.tmpu_ <- .py_method(b, "apply_support", 15L, type = "fixed"); r15 <- .py_getitem(.tmpu_, 0L); m15 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
.py_method(b, "apply_rotation_hinge", 5L)
})
.py_try({
.py_method(b, "apply_rotation_hinge", 12L)
})
.py_try({
.py_method(b, "apply_load", -10L, 5L, -1L)
})
.py_try({
.py_method(b, "apply_load", -5L, 10L, 0L, 15L)
})
.py_try({
.py_method(b, "solve_for_reaction_loads", r0, r10, r15, m15)
})
.py_try({
.py_rec("beam:test_point_cflexure_slow_cases:534", function() .py_method(b, "point_cflexure"))
})
.py_try({
E <- Symbol("E")
})
.py_try({
E <- 200L
})
.py_try({
I <- Symbol("I")
})
.py_try({
I <- 1000L
})
.py_try({
b1 <- Beam(5L, E, I)
})
.py_try({
b2 <- Beam(10L, E, I)
})
.py_try({
b <- .py_method(b1, "join", b2, via = "hinge")
})
.py_try({
.py_method(b, "apply_support", 0L, type = "pin")
})
.py_try({
.py_method(b, "apply_support", 10L, type = "pin")
})
.py_try({
.py_method(b, "apply_support", 15L, type = "fixed")
})
.py_try({
.py_method(b, "apply_load", -10L, 5L, -1L)
})
.py_try({
.py_method(b, "apply_load", -5L, 0L, 0L, 5L)
})
.py_try({
.py_method(b, "apply_load", -5L, 10L, 0L, 15L)
})
.py_try({
.tmpu_ <- .py_symbols("R_0, R_10, R_15, M_15"); r0 <- .py_getitem(.tmpu_, 0L); r10 <- .py_getitem(.tmpu_, 1L); r15 <- .py_getitem(.tmpu_, 2L); m15 <- .py_getitem(.tmpu_, 3L)
})
.py_try({
.py_method(b, "solve_for_reaction_loads", r0, r10, r15, m15)
})
.py_try({
.py_rec("beam:test_point_cflexure_slow_cases:552", function() .py_method(b, "point_cflexure"))
})
.py_try({
E <- Symbol("E")
})
.py_try({
I <- Symbol("I")
})
.py_try({
b <- Beam(10L, E, I)
})
.py_try({
.py_method(b, "apply_support", 0L, type = "pin")
})
.py_try({
.py_method(b, "apply_support", 10L, type = "pin")
})
.py_try({
.py_method(b, "apply_load", -20L, 5L, -1L)
})
.py_try({
.tmpu_ <- .py_symbols("R_0, R_10"); R_0 <- .py_getitem(.tmpu_, 0L); R_10 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
.py_method(b, "solve_for_reaction_loads", R_0, R_10)
})
.py_try({
.py_rec("beam:test_point_cflexure_slow_cases:563", function() .py_method(b, "point_cflexure"))
})
.py_try({
E <- Symbol("E")
})
.py_try({
I <- Symbol("I")
})
.py_try({
P <- Symbol("P", positive = TRUE)
})
.py_try({
Z <- Symbol("Z", positive = TRUE)
})
.py_try({
b1 <- Beam(5L, E, I)
})
.py_try({
b2 <- Beam(3L, E, I)
})
.py_try({
b <- .py_method(b1, "join", b2, via = "hinge")
})
.py_try({
.py_method(b, "apply_support", 0L, type = "fixed")
})
.py_try({
.py_method(b, "apply_support", 8L, type = "pin")
})
.py_try({
.py_method(b, "apply_load", .py_neg(P), 2L, -1L)
})
.py_try({
.py_method(b, "apply_load", .py_neg(Z), 5L, 0L, 8L)
})
.py_try({
.tmpu_ <- .py_symbols("R_0, M_0, R_8"); r0 <- .py_getitem(.tmpu_, 0L); m0 <- .py_getitem(.tmpu_, 1L); r8 <- .py_getitem(.tmpu_, 2L)
})
.py_try({
.py_method(b, "solve_for_reaction_loads", r0, m0, r8)
})
.py_try({
.py_rec("beam:test_point_cflexure_slow_cases:579", function() .py_method(b, "point_cflexure"))
})
.py_try({
E <- Symbol("E")
})
.py_try({
I <- Symbol("I")
})
.py_try({
b1 <- Beam(5L, E, I)
})
.py_try({
b2 <- Beam(3L, E, I)
})
.py_try({
b <- .py_method(b1, "join", b2, via = "hinge")
})
.py_try({
.py_method(b, "apply_support", 0L, type = "fixed")
})
.py_try({
.py_method(b, "apply_support", 8L, type = "pin")
})
.py_try({
.py_method(b, "apply_load", -10L, 5L, -1L)
})
.py_try({
.py_method(b, "apply_load", 10L, 7L, -1L)
})
.py_try({
.tmpu_ <- .py_symbols("R_0, M_0, R_8"); r0 <- .py_getitem(.tmpu_, 0L); m0 <- .py_getitem(.tmpu_, 1L); r8 <- .py_getitem(.tmpu_, 2L)
})
.py_try({
.py_method(b, "solve_for_reaction_loads", r0, m0, r8)
})
.py_try({
.py_rec("beam:test_point_cflexure_slow_cases:593", function() .py_method(b, "point_cflexure"))
})
.py_try({
E <- Symbol("E")
})
.py_try({
I <- Symbol("I")
})
.py_try({
b1 <- Beam(5L, E, I)
})
.py_try({
b2 <- Beam(3L, E, I)
})
.py_try({
b <- .py_method(b1, "join", b2, via = "hinge")
})
.py_try({
.py_method(b, "apply_support", 0L, type = "fixed")
})
.py_try({
.py_method(b, "apply_support", 8L, type = "pin")
})
.py_try({
.py_method(b, "apply_load", -10L, 5L, -1L)
})
.py_try({
.tmpu_ <- .py_symbols("R_0, M_0, R_8"); r0 <- .py_getitem(.tmpu_, 0L); m0 <- .py_getitem(.tmpu_, 1L); r8 <- .py_getitem(.tmpu_, 2L)
})
.py_try({
.py_method(b, "solve_for_reaction_loads", r0, m0, r8)
})
.py_try({
.py_rec("beam:test_point_cflexure_slow_cases:606", function() .py_method(b, "point_cflexure"))
})
.py_try({
E <- Symbol("E")
})
.py_try({
I <- Symbol("I")
})
.py_try({
b <- Beam(7L, E, I)
})
.py_try({
.tmpu_ <- .py_method(b, "apply_support", 0L, type = "fixed"); r0 <- .py_getitem(.tmpu_, 0L); m0 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
r4 <- .py_method(b, "apply_support", 4L, type = "pin")
})
.py_try({
r7 <- .py_method(b, "apply_support", 7L, type = "pin")
})
.py_try({
.py_method(b, "apply_rotation_hinge", 2L)
})
.py_try({
.py_method(b, "apply_rotation_hinge", 5L)
})
.py_try({
.py_method(b, "apply_load", -5L, 0L, 0L, 2L)
})
.py_try({
.py_method(b, "apply_load", -10L, 6L, -1L)
})
.py_try({
.py_method(b, "solve_for_reaction_loads", r0, m0, r4, r7)
})
.py_try({
.py_rec("beam:test_point_cflexure_slow_cases:620", function() .py_method(b, "point_cflexure"))
})
.py_try({
E <- Symbol("E")
})
.py_try({
I <- Symbol("I")
})
.py_try({
b <- Beam(10L, E, I)
})
.py_try({
.tmpu_ <- .py_method(b, "apply_support", 0L, type = "fixed"); r0 <- .py_getitem(.tmpu_, 0L); m0 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
r4 <- .py_method(b, "apply_support", 4L, type = "pin")
})
.py_try({
.tmpu_ <- .py_method(b, "apply_support", 7L, type = "fixed"); r7 <- .py_getitem(.tmpu_, 0L); m7 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
.py_method(b, "apply_rotation_hinge", 2L)
})
.py_try({
.py_method(b, "apply_rotation_hinge", 5L)
})
.py_try({
.py_method(b, "apply_load", -5L, 0L, 0L, 2L)
})
.py_try({
.py_method(b, "apply_load", -10L, 6L, -1L)
})
.py_try({
.py_method(b, "apply_load", 7L, 3L, -1L)
})
.py_try({
.py_method(b, "solve_for_reaction_loads", r0, m0, r4, r7, m7)
})
.py_try({
.py_rec("beam:test_point_cflexure_slow_cases:634", function() .py_method(b, "point_cflexure"))
})
.py_try({
E <- Symbol("E")
})
.py_try({
I <- Symbol("I")
})
.py_try({
b <- Beam(10L, E, I)
})
.py_try({
.tmpu_ <- .py_method(b, "apply_support", 0L, type = "fixed"); r0 <- .py_getitem(.tmpu_, 0L); m0 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
.tmpu_ <- .py_method(b, "apply_support", 4L, type = "fixed"); r4 <- .py_getitem(.tmpu_, 0L); m4 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
.tmpu_ <- .py_method(b, "apply_support", 7L, type = "fixed"); r7 <- .py_getitem(.tmpu_, 0L); m7 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
.py_method(b, "apply_rotation_hinge", 2L)
})
.py_try({
.py_method(b, "apply_rotation_hinge", 5L)
})
.py_try({
.py_method(b, "apply_load", -5L, 0L, 0L, 2L)
})
.py_try({
.py_method(b, "apply_load", -10L, 6L, -1L)
})
.py_try({
.py_method(b, "apply_load", 7L, 3L, -1L)
})
.py_try({
.py_method(b, "solve_for_reaction_loads", r0, m0, r4, m4, r7, m7)
})
.py_try({
.py_rec("beam:test_point_cflexure_slow_cases:649", function() .py_method(b, "point_cflexure"))
}))---")
  expect_identical(got[["beam:test_point_cflexure_slow_cases:476"]], "[173/260, 2, 1490/381, 5]")
  expect_identical(got[["beam:test_point_cflexure_slow_cases:489"]], "[25/3, 12]")
  expect_identical(got[["beam:test_point_cflexure_slow_cases:505"]], "[5/2, 15/2]")
  expect_identical(got[["beam:test_point_cflexure_slow_cases:519"]], "[173/260, 2, 1490/381, 5]")
  expect_identical(got[["beam:test_point_cflexure_slow_cases:534"]], "[12]")
  expect_identical(got[["beam:test_point_cflexure_slow_cases:552"]], "[5, 149/8 - 3*sqrt(209)/8]")
  expect_identical(got[["beam:test_point_cflexure_slow_cases:563"]], "[]")
  expect_identical(got[["beam:test_point_cflexure_slow_cases:579"]], "[5]")
  expect_identical(got[["beam:test_point_cflexure_slow_cases:593"]], "[]")
  expect_identical(got[["beam:test_point_cflexure_slow_cases:606"]], "[]")
  expect_identical(got[["beam:test_point_cflexure_slow_cases:620"]], "[1, 2, 5]")
  expect_identical(got[["beam:test_point_cflexure_slow_cases:634"]], "[73/260, 2, 5, 2705/437]")
  expect_identical(got[["beam:test_point_cflexure_slow_cases:649"]], "[13/16, 2, 482/129, 4, 5, 83/13]")
})

test_that("beam test_remove_load matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
x <- Symbol("x")
})
.py_try({
y <- Symbol("y")
})
.py_try({
.tmpu_ <- .py_symbols("R1, R2"); R1 <- .py_getitem(.tmpu_, 0L); R2 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
E <- Symbol("E")
})
.py_try({
I <- Symbol("I")
})
.py_try({
b <- Beam(4L, E, I)
})
.py_try({
.ok_ <- tryCatch({ .py_method(b, "remove_load", 2L, 1L, -1L); FALSE }, error = function(e) if (identical(.err_class(e), "ValueError")) TRUE else stop(e))
if (.ok_) .py_rec("beam:test_remove_load:661", function() TRUE) else .py_rec("beam:test_remove_load:663", function() FALSE)
})
.py_try({
.py_method(b, "apply_load", -3L, 0L, -2L)
})
.py_try({
.py_method(b, "apply_load", 4L, 2L, -1L)
})
.py_try({
.py_method(b, "apply_load", -2L, 2L, 2L, end = 3L)
})
.py_try({
.py_method(b, "remove_load", -2L, 2L, 2L, end = 3L)
})
.py_try({
.py_rec("beam:test_remove_load:669", function() .py_attr(b, "load"))
})
.py_try({
.py_rec("beam:test_remove_load:670", function() .py_attr(b, "applied_loads"))
})
.py_try({
.ok_ <- tryCatch({ .py_method(b, "remove_load", 1L, 2L, -1L); FALSE }, error = function(e) if (identical(.err_class(e), "ValueError")) TRUE else stop(e))
if (.ok_) .py_rec("beam:test_remove_load:677", function() TRUE) else .py_rec("beam:test_remove_load:679", function() FALSE)
})
.py_try({
.py_method(b, "remove_load", -3L, 0L, -2L)
})
.py_try({
.py_method(b, "remove_load", 4L, 2L, -1L)
})
.py_try({
.py_rec("beam:test_remove_load:683", function() .py_attr(b, "load"))
})
.py_try({
.py_rec("beam:test_remove_load:684", function() .py_attr(b, "applied_loads"))
}))---")
  expect_identical(got[["beam:test_remove_load:661"]], "True")
  expect_identical(got[["beam:test_remove_load:669"]], "-3*SingularityFunction(x, 0, -2) + 4*SingularityFunction(x, 2, -1)")
  expect_identical(got[["beam:test_remove_load:670"]], "[(-3, 0, -2, None), (4, 2, -1, None)]")
  expect_identical(got[["beam:test_remove_load:677"]], "True")
  expect_identical(got[["beam:test_remove_load:683"]], "0")
  expect_identical(got[["beam:test_remove_load:684"]], "[]")
})

test_that("beam test_apply_support matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
x <- Symbol("x")
})
.py_try({
y <- Symbol("y")
})
.py_try({
.tmpu_ <- .py_symbols("R1, R2"); R1 <- .py_getitem(.tmpu_, 0L); R2 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
E <- Symbol("E")
})
.py_try({
I <- Symbol("I")
})
.py_try({
b <- Beam(4L, E, I)
})
.py_try({
.py_method(b, "apply_support", 0L, "fixed")
})
.py_try({
.py_method(b, "apply_load", 20L, 4L, -1L)
})
.py_try({
.tmpu_ <- .py_symbols("M_0, R_0"); M_0 <- .py_getitem(.tmpu_, 0L); R_0 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
.py_method(b, "solve_for_reaction_loads", R_0, M_0)
})
.py_try({
.py_rec("beam:test_apply_support:696", function() simplify(.py_method(b, "slope")))
})
.py_try({
.py_rec("beam:test_apply_support:698", function() simplify(.py_method(b, "deflection")))
})
.py_try({
b <- Beam(30L, E, I)
})
.py_try({
p0 <- .py_method(b, "apply_support", 10L, "pin")
})
.py_try({
p1 <- .py_method(b, "apply_support", 30L, "roller")
})
.py_try({
.py_method(b, "apply_load", -8L, 0L, -1L)
})
.py_try({
.py_method(b, "apply_load", 120L, 30L, -2L)
})
.py_try({
.py_method(b, "solve_for_reaction_loads", p0, p1)
})
.py_try({
.py_rec("beam:test_apply_support:707", function() .py_method(b, "slope"))
})
.py_try({
.py_rec("beam:test_apply_support:709", function() .py_method(b, "deflection"))
})
.py_try({
R_10 <- Symbol("R_10")
})
.py_try({
R_30 <- Symbol("R_30")
})
.py_try({
.py_rec("beam:test_apply_support:713", function() p0)
})
.py_try({
.py_rec("beam:test_apply_support:715", function() .py_getitem(.py_attr(b, "reaction_loads"), p0))
})
.py_try({
b <- Beam(8L, E, I)
})
.py_try({
.tmpu_ <- .py_method(b, "apply_support", 0L, "fixed"); p0 <- .py_getitem(.tmpu_, 0L); m0 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
p1 <- .py_method(b, "apply_support", 8L, "roller")
})
.py_try({
.py_method(b, "apply_load", -5L, 0L, 0L, 8L)
})
.py_try({
.py_method(b, "solve_for_reaction_loads", p0, m0, p1)
})
.py_try({
R_0 <- Symbol("R_0")
})
.py_try({
M_0 <- Symbol("M_0")
})
.py_try({
R_8 <- Symbol("R_8")
})
.py_try({
.py_rec("beam:test_apply_support:725", function() p0)
})
.py_try({
.py_rec("beam:test_apply_support:726", function() m0)
})
.py_try({
.py_rec("beam:test_apply_support:727", function() p1)
})
.py_try({
.py_rec("beam:test_apply_support:729", function() .py_getitem(.py_attr(b, "reaction_loads"), m0))
})
.py_try({
b <- Beam(8L, E, I)
})
.py_try({
.tmpu_ <- .py_method(b, "apply_support", 0L, "fixed"); p0 <- .py_getitem(.tmpu_, 0L); m0 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
p1 <- .py_method(b, "apply_support", 8L, "roller")
})
.py_try({
.py_method(b, "apply_load", -5L, 0L, 0L, 8L)
})
.py_try({
.py_method(b, "solve_for_reaction_loads")
})
.py_try({
R_0 <- Symbol("R_0")
})
.py_try({
M_0 <- Symbol("M_0")
})
.py_try({
R_8 <- Symbol("R_8")
})
.py_try({
.py_rec("beam:test_apply_support:739", function() p0)
})
.py_try({
.py_rec("beam:test_apply_support:740", function() m0)
})
.py_try({
.py_rec("beam:test_apply_support:741", function() p1)
})
.py_try({
.py_rec("beam:test_apply_support:743", function() .py_getitem(.py_attr(b, "reaction_loads"), m0))
})
.py_try({
P <- Symbol("P", positive = TRUE)
})
.py_try({
L <- Symbol("L", positive = TRUE)
})
.py_try({
b <- Beam(L, E, I)
})
.py_try({
.py_method(b, "apply_support", 0L, type = "fixed")
})
.py_try({
.py_method(b, "apply_support", L, type = "fixed")
})
.py_try({
.py_method(b, "apply_load", .py_neg(P), .py_div(L, 2L), -1L)
})
.py_try({
.tmpu_ <- .py_symbols("R_0, R_L, M_0, M_L"); R_0 <- .py_getitem(.tmpu_, 0L); R_L <- .py_getitem(.tmpu_, 1L); M_0 <- .py_getitem(.tmpu_, 2L); M_L <- .py_getitem(.tmpu_, 3L)
})
.py_try({
.py_method(b, "solve_for_reaction_loads", R_0, R_L, M_0, M_L)
})
.py_try({
P <- Symbol("P", positive = TRUE)
})
.py_try({
L <- Symbol("L", positive = TRUE)
})
.py_try({
b <- Beam(L, E, I)
})
.py_try({
.py_method(b, "apply_support", 0L, type = "fixed")
})
.py_try({
.py_method(b, "apply_support", L, type = "fixed")
})
.py_try({
.py_method(b, "apply_load", .py_neg(P), .py_div(L, 2L), -1L)
})
.py_try({
.tmpu_ <- .py_symbols("R_0, R_L, M_0, M_L"); R_0 <- .py_getitem(.tmpu_, 0L); R_L <- .py_getitem(.tmpu_, 1L); M_0 <- .py_getitem(.tmpu_, 2L); M_L <- .py_getitem(.tmpu_, 3L)
})
.py_try({
.py_method(b, "solve_for_reaction_loads")
})
.py_try({
P <- Symbol("P", positive = TRUE)
})
.py_try({
L <- Symbol("L", positive = TRUE)
})
.py_try({
c <- Beam(L, E, I)
})
.py_try({
.py_method(c, "apply_support", 0L, type = "fixed")
})
.py_try({
.py_method(c, "apply_support", L, type = "fixed")
})
.py_try({
.py_method(c, "apply_load", .py_neg(P), .py_div(L, 2L), -1L)
})
.py_try({
.tmpu_ <- .py_symbols("R_0, R_L, M_0, M_L"); R_0 <- .py_getitem(.tmpu_, 0L); R_L <- .py_getitem(.tmpu_, 1L); M_0 <- .py_getitem(.tmpu_, 2L); M_L <- .py_getitem(.tmpu_, 3L)
})
.py_try({
.py_method(c, "solve_for_reaction_loads")
})
.py_try({
E <- Symbol("E")
})
.py_try({
I <- Symbol("I")
})
.py_try({
b <- Beam(6L, E, I)
})
.t778 <- function() {
  .py_method(b, "apply_support", 0L, "pen")
}
.py_try({
.py_rec("beam:test_apply_support:778", .t778)
}))---")
  expect_identical(got[["beam:test_apply_support:696"]], "10*(8*SingularityFunction(x, 0, 1) - SingularityFunction(x, 0, 2) + SingularityFunction(x, 4, 2))/(E*I)")
  expect_identical(got[["beam:test_apply_support:698"]], "10*(12*SingularityFunction(x, 0, 2) - SingularityFunction(x, 0, 3) + SingularityFunction(x, 4, 3))/(3*E*I)")
  expect_identical(got[["beam:test_apply_support:707"]], "(-4*SingularityFunction(x, 0, 2) + 3*SingularityFunction(x, 10, 2) + 120*SingularityFunction(x, 30, 1) + SingularityFunction(x, 30, 2) + 4000/3)/(E*I)")
  expect_identical(got[["beam:test_apply_support:709"]], "(4000*x/3 - 4*SingularityFunction(x, 0, 3)/3 + SingularityFunction(x, 10, 3) + 60*SingularityFunction(x, 30, 2) + SingularityFunction(x, 30, 3)/3 - 12000)/(E*I)")
  expect_identical(got[["beam:test_apply_support:713"]], "R_10")
  expect_identical(got[["beam:test_apply_support:715"]], "6")
  expect_identical(got[["beam:test_apply_support:725"]], "R_0")
  expect_identical(got[["beam:test_apply_support:726"]], "M_0")
  expect_identical(got[["beam:test_apply_support:727"]], "R_8")
  expect_identical(got[["beam:test_apply_support:729"]], "-40")
  expect_identical(got[["beam:test_apply_support:739"]], "R_0")
  expect_identical(got[["beam:test_apply_support:740"]], "M_0")
  expect_identical(got[["beam:test_apply_support:741"]], "R_8")
  expect_identical(got[["beam:test_apply_support:743"]], "-40")
  expect_identical(got[["beam:test_apply_support:778"]], "ERR ValueError: Invalid support type. Choose from 'pin', 'roller', or 'fixed'.")
})

test_that("beam test_for_duplicate_symbols_to_solve_for_reaction_loads matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
x <- Symbol("x")
})
.py_try({
y <- Symbol("y")
})
.py_try({
.tmpu_ <- .py_symbols("R1, R2"); R1 <- .py_getitem(.tmpu_, 0L); R2 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
.tmpu_ <- .py_symbols("E I"); E <- .py_getitem(.tmpu_, 0L); I <- .py_getitem(.tmpu_, 1L)
})
.py_try({
b <- Beam(4L, E, I)
})
.py_try({
.py_method(b, "apply_load", 20L, 4L, -1L)
})
.py_try({
.tmpu_ <- .py_symbols("R, M"); R <- .py_getitem(.tmpu_, 0L); M <- .py_getitem(.tmpu_, 1L)
})
.py_try({
.py_method(b, "apply_load", R, 0L, -1L)
})
.py_try({
.py_method(b, "apply_load", M, 0L, -2L)
})
.py_try({
.py_setattr(b, "bc_deflection", list(.py_tup(0L, 0L)))
})
.py_try({
.py_setattr(b, "bc_slope", list(.py_tup(0L, 0L)))
})
.py_try({
.py_method(b, "solve_for_reaction_loads", R, M)
})
.t792 <- function() {
  .py_method(b, "solve_for_reaction_loads", R, R)
}
.py_try({
.py_rec("beam:test_for_duplicate_symbols_to_solve_for_reaction_loads:792", .t792)
}))---")
  expect_identical(got[["beam:test_for_duplicate_symbols_to_solve_for_reaction_loads:792"]], "ERR ValueError: Duplicate Symbols passed into solve_for_reaction_loads()")
})

test_that("beam test_apply_rotation_hinge matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
x <- Symbol("x")
})
.py_try({
y <- Symbol("y")
})
.py_try({
.tmpu_ <- .py_symbols("R1, R2"); R1 <- .py_getitem(.tmpu_, 0L); R2 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
b <- Beam(15L, 20L, 20L)
})
.py_try({
.tmpu_ <- .py_method(b, "apply_support", 0L, type = "fixed"); r0 <- .py_getitem(.tmpu_, 0L); m0 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
r10 <- .py_method(b, "apply_support", 10L, type = "pin")
})
.py_try({
r15 <- .py_method(b, "apply_support", 15L, type = "pin")
})
.py_try({
p7 <- .py_method(b, "apply_rotation_hinge", 7L)
})
.py_try({
p12 <- .py_method(b, "apply_rotation_hinge", 12L)
})
.py_try({
.py_method(b, "apply_load", -10L, 7L, -1L)
})
.py_try({
.py_method(b, "apply_load", -2L, 10L, 0L, 15L)
})
.py_try({
.py_method(b, "solve_for_reaction_loads", r0, m0, r10, r15)
})
.py_try({
.tmpu_ <- .py_symbols("R_0, M_0, R_10, R_15, P_7, P_12"); R_0 <- .py_getitem(.tmpu_, 0L); M_0 <- .py_getitem(.tmpu_, 1L); R_10 <- .py_getitem(.tmpu_, 2L); R_15 <- .py_getitem(.tmpu_, 3L); P_7 <- .py_getitem(.tmpu_, 4L); P_12 <- .py_getitem(.tmpu_, 5L)
})
.py_try({
expected_reactions <- .py_dict(list(R_0, M_0, R_10, R_15), list(.py_div(20L, 3L), .py_div(-140L, 3L), .py_div(31L, 3L), 3L))
})
.py_try({
expected_rotations <- .py_dict(list(P_7, P_12), list(.py_div(2281L, 2160L), .py_div(-5137L, 5184L)))
})
.py_try({
reaction_symbols <- list(r0, m0, r10, r15)
})
.py_try({
rotation_symbols <- list(p7, p12)
})
.py_try({
tolerance <- .py_flt(1e-06)
})
.py_try({
.py_rec("beam:test_apply_rotation_hinge:812", function() .py_all(.py_comp(reaction_symbols, function(.it) { r <- .it; list(.py_lt(.py_abs(.py_sub(.py_getitem(.py_attr(b, "reaction_loads"), r), .py_getitem(expected_reactions, r))), tolerance)) })))
})
.py_try({
.py_rec("beam:test_apply_rotation_hinge:813", function() .py_all(.py_comp(rotation_symbols, function(.it) { r <- .it; list(.py_lt(.py_abs(.py_sub(.py_getitem(.py_attr(b, "rotation_jumps"), r), .py_getitem(expected_rotations, r))), tolerance)) })))
})
.py_try({
expected_bending_moment <- .py_sub(.py_sub(.py_add(.py_add(.py_sub(.py_add(.py_sub(.py_sub(.py_div(.py_mul(140L, SingularityFunction(x, 0L, 0L)), 3L), .py_div(.py_mul(20L, SingularityFunction(x, 0L, 1L)), 3L)), .py_div(.py_mul(11405L, SingularityFunction(x, 7L, -1L)), 27L)), .py_mul(10L, SingularityFunction(x, 7L, 1L))), .py_div(.py_mul(31L, SingularityFunction(x, 10L, 1L)), 3L)), SingularityFunction(x, 10L, 2L)), .py_div(.py_mul(128425L, SingularityFunction(x, 12L, -1L)), 324L)), .py_mul(3L, SingularityFunction(x, 15L, 1L))), SingularityFunction(x, 15L, 2L))
})
.py_try({
.py_rec("beam:test_apply_rotation_hinge:819", function() .py_method(.py_method(b, "bending_moment"), "expand"))
})
.py_try({
expected_slope <- .py_add(.py_add(.py_sub(.py_sub(.py_add(.py_sub(.py_add(.py_add(.py_div(.py_mul(-7L, SingularityFunction(x, 0L, 1L)), 60L), .py_div(SingularityFunction(x, 0L, 2L), 120L)), .py_div(.py_mul(2281L, SingularityFunction(x, 7L, 0L)), 2160L)), .py_div(SingularityFunction(x, 7L, 2L), 80L)), .py_div(.py_mul(31L, SingularityFunction(x, 10L, 2L)), 2400L)), .py_div(SingularityFunction(x, 10L, 3L), 1200L)), .py_div(.py_mul(5137L, SingularityFunction(x, 12L, 0L)), 5184L)), .py_div(.py_mul(3L, SingularityFunction(x, 15L, 2L)), 800L)), .py_div(SingularityFunction(x, 15L, 3L), 1200L))
})
.py_try({
.py_rec("beam:test_apply_rotation_hinge:825", function() .py_method(.py_method(b, "slope"), "expand"))
})
.py_try({
expected_deflection <- .py_add(.py_add(.py_sub(.py_sub(.py_add(.py_sub(.py_add(.py_add(.py_div(.py_mul(-7L, SingularityFunction(x, 0L, 2L)), 120L), .py_div(SingularityFunction(x, 0L, 3L), 360L)), .py_div(.py_mul(2281L, SingularityFunction(x, 7L, 1L)), 2160L)), .py_div(SingularityFunction(x, 7L, 3L), 240L)), .py_div(.py_mul(31L, SingularityFunction(x, 10L, 3L)), 7200L)), .py_div(SingularityFunction(x, 10L, 4L), 4800L)), .py_div(.py_mul(5137L, SingularityFunction(x, 12L, 1L)), 5184L)), .py_div(SingularityFunction(x, 15L, 3L), 800L)), .py_div(SingularityFunction(x, 15L, 4L), 4800L))
})
.py_try({
.py_rec("beam:test_apply_rotation_hinge:831", function() .py_method(.py_method(b, "deflection"), "expand"))
})
.py_try({
E <- Symbol("E")
})
.py_try({
I <- Symbol("I")
})
.py_try({
F <- Symbol("F")
})
.py_try({
b <- Beam(10L, E, I)
})
.py_try({
.tmpu_ <- .py_method(b, "apply_support", 0L, type = "fixed"); r0 <- .py_getitem(.tmpu_, 0L); m0 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
r10 <- .py_method(b, "apply_support", 10L, type = "pin")
})
.py_try({
.py_method(b, "apply_rotation_hinge", 6L)
})
.py_try({
.py_method(b, "apply_load", F, 8L, -1L)
})
.py_try({
.py_method(b, "solve_for_reaction_loads", r0, m0, r10)
})
.py_try({
.py_rec("beam:test_apply_rotation_hinge:843", function() .py_method(b, "bending_moment"))
})
.py_try({
expected_deflection <- .py_div(.py_neg(.py_add(.py_sub(.py_add(.py_add(.py_div(.py_mul(.py_mul(-3L, F), SingularityFunction(x, 0L, 2L)), 2L), .py_div(.py_mul(F, SingularityFunction(x, 0L, 3L)), 12L)), .py_mul(.py_mul(17L, F), SingularityFunction(x, 6L, 1L))), .py_div(.py_mul(F, SingularityFunction(x, 8L, 3L)), 6L)), .py_div(.py_mul(F, SingularityFunction(x, 10L, 3L)), 12L))), .py_mul(E, I))
})
.py_try({
.py_rec("beam:test_apply_rotation_hinge:849", function() .py_method(.py_method(b, "deflection"), "expand"))
})
.py_try({
E <- Symbol("E")
})
.py_try({
I <- Symbol("I")
})
.py_try({
F <- Symbol("F")
})
.py_try({
l1 <- Symbol("l1", positive = TRUE)
})
.py_try({
l2 <- Symbol("l2", positive = TRUE)
})
.py_try({
l3 <- Symbol("l3", positive = TRUE)
})
.py_try({
L <- .py_add(.py_add(l1, l2), l3)
})
.py_try({
b <- Beam(L, E, I)
})
.py_try({
.tmpu_ <- .py_method(b, "apply_support", 0L, type = "fixed"); r0 <- .py_getitem(.tmpu_, 0L); m0 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
r1 <- .py_method(b, "apply_support", L, type = "pin")
})
.py_try({
.py_method(b, "apply_rotation_hinge", l1)
})
.py_try({
.py_method(b, "apply_load", F, .py_add(l1, l2), -1L)
})
.py_try({
.py_method(b, "solve_for_reaction_loads", r0, m0, r1)
})
.py_try({
.py_rec("beam:test_apply_rotation_hinge:864", function() .py_getitem(.py_attr(b, "reaction_loads"), r0))
})
.py_try({
.py_rec("beam:test_apply_rotation_hinge:865", function() .py_getitem(.py_attr(b, "reaction_loads"), m0))
})
.py_try({
.py_rec("beam:test_apply_rotation_hinge:866", function() .py_getitem(.py_attr(b, "reaction_loads"), r1))
})
.py_try({
expected_bending_moment <- .py_sub(.py_sub(.py_add(.py_add(.py_div(.py_mul(.py_mul(.py_mul(.py_neg(F), l1), l3), SingularityFunction(x, 0L, 0L)), .py_add(l2, l3)), .py_div(.py_mul(.py_mul(F, l2), SingularityFunction(x, .py_add(.py_add(l1, l2), l3), 1L)), .py_add(l2, l3))), .py_div(.py_mul(.py_mul(F, l3), SingularityFunction(x, 0L, 1L)), .py_add(l2, l3))), .py_mul(F, SingularityFunction(x, .py_add(l1, l2), 1L))), .py_div(.py_mul(.py_add(.py_add(.py_add(.py_sub(.py_sub(.py_mul(.py_mul(.py_mul(-2L, F), .py_pow(l1, 3L)), l3), .py_mul(.py_mul(.py_mul(.py_mul(3L, F), .py_pow(l1, 2L)), l2), l3)), .py_mul(.py_mul(.py_mul(3L, F), .py_pow(l1, 2L)), .py_pow(l3, 2L))), .py_mul(.py_mul(F, .py_pow(l2, 3L)), l3)), .py_mul(.py_mul(.py_mul(3L, F), .py_pow(l2, 2L)), .py_pow(l3, 2L))), .py_mul(.py_mul(.py_mul(2L, F), l2), .py_pow(l3, 3L))), SingularityFunction(x, l1, -1L)), .py_add(.py_add(.py_mul(6L, .py_pow(l2, 2L)), .py_mul(.py_mul(12L, l2), l3)), .py_mul(6L, .py_pow(l3, 2L)))))
})
.py_try({
.py_rec("beam:test_apply_rotation_hinge:872", function() simplify(.py_method(.py_method(b, "bending_moment"), "expand")))
}))---")
  expect_identical(got[["beam:test_apply_rotation_hinge:812"]], "True")
  expect_identical(got[["beam:test_apply_rotation_hinge:813"]], "True")
  expect_identical(got[["beam:test_apply_rotation_hinge:819"]], "140*SingularityFunction(x, 0, 0)/3 - 20*SingularityFunction(x, 0, 1)/3 - 11405*SingularityFunction(x, 7, -1)/27 + 10*SingularityFunction(x, 7, 1) - 31*SingularityFunction(x, 10, 1)/3 + SingularityFunction(x, 10, 2) + 128425*SingularityFunction(x, 12, -1)/324 - 3*SingularityFunction(x, 15, 1) - SingularityFunction(x, 15, 2)")
  expect_identical(got[["beam:test_apply_rotation_hinge:825"]], "-7*SingularityFunction(x, 0, 1)/60 + SingularityFunction(x, 0, 2)/120 + 2281*SingularityFunction(x, 7, 0)/2160 - SingularityFunction(x, 7, 2)/80 + 31*SingularityFunction(x, 10, 2)/2400 - SingularityFunction(x, 10, 3)/1200 - 5137*SingularityFunction(x, 12, 0)/5184 + 3*SingularityFunction(x, 15, 2)/800 + SingularityFunction(x, 15, 3)/1200")
  expect_identical(got[["beam:test_apply_rotation_hinge:831"]], "-7*SingularityFunction(x, 0, 2)/120 + SingularityFunction(x, 0, 3)/360 + 2281*SingularityFunction(x, 7, 1)/2160 - SingularityFunction(x, 7, 3)/240 + 31*SingularityFunction(x, 10, 3)/7200 - SingularityFunction(x, 10, 4)/4800 - 5137*SingularityFunction(x, 12, 1)/5184 + SingularityFunction(x, 15, 3)/800 + SingularityFunction(x, 15, 4)/4800")
  expect_identical(got[["beam:test_apply_rotation_hinge:843"]], "-3*F*SingularityFunction(x, 0, 0) + F*SingularityFunction(x, 0, 1)/2 + 17*F*SingularityFunction(x, 6, -1) - F*SingularityFunction(x, 8, 1) + F*SingularityFunction(x, 10, 1)/2")
  expect_identical(got[["beam:test_apply_rotation_hinge:849"]], "3*F*SingularityFunction(x, 0, 2)/(2*E*I) - F*SingularityFunction(x, 0, 3)/(12*E*I) - 17*F*SingularityFunction(x, 6, 1)/(E*I) + F*SingularityFunction(x, 8, 3)/(6*E*I) - F*SingularityFunction(x, 10, 3)/(12*E*I)")
  expect_identical(got[["beam:test_apply_rotation_hinge:864"]], "-F*l3/(l2 + l3)")
  expect_identical(got[["beam:test_apply_rotation_hinge:865"]], "F*l1*l3/(l2 + l3)")
  expect_identical(got[["beam:test_apply_rotation_hinge:866"]], "-F*l2/(l2 + l3)")
  expect_identical(got[["beam:test_apply_rotation_hinge:872"]], "F*(l3*(l2 + l3)*(2*l1^3 + 3*l1^2*l2 + 3*l1^2*l3 - l2^3 - 3*l2^2*l3 - 2*l2*l3^2)*SingularityFunction(x, l1, -1) - 6*(l2 + l3)*(l2^2 + 2*l2*l3 + l3^2)*SingularityFunction(x, l1 + l2, 1) + 6*(l2^2 + 2*l2*l3 + l3^2)*(-l1*l3*SingularityFunction(x, 0, 0) + l2*SingularityFunction(x, l1 + l2 + l3, 1) + l3*SingularityFunction(x, 0, 1)))/(6*(l2 + l3)*(l2^2 + 2*l2*l3 + l3^2))")
})

test_that("beam test_apply_sliding_hinge matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
x <- Symbol("x")
})
.py_try({
y <- Symbol("y")
})
.py_try({
.tmpu_ <- .py_symbols("R1, R2"); R1 <- .py_getitem(.tmpu_, 0L); R2 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
b <- Beam(13L, 20L, 20L)
})
.py_try({
.tmpu_ <- .py_method(b, "apply_support", 0L, type = "fixed"); r0 <- .py_getitem(.tmpu_, 0L); m0 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
w8 <- .py_method(b, "apply_sliding_hinge", 8L)
})
.py_try({
r13 <- .py_method(b, "apply_support", 13L, type = "pin")
})
.py_try({
.py_method(b, "apply_load", -10L, 5L, -1L)
})
.py_try({
.py_method(b, "solve_for_reaction_loads", r0, m0, r13)
})
.py_try({
.tmpu_ <- .py_symbols("R_0, M_0, R_13 W_8"); R_0 <- .py_getitem(.tmpu_, 0L); M_0 <- .py_getitem(.tmpu_, 1L); R_13 <- .py_getitem(.tmpu_, 2L); W_8 <- .py_getitem(.tmpu_, 3L)
})
.py_try({
tolerance <- .py_flt(1e-06)
})
.py_try({
.py_rec("beam:test_apply_sliding_hinge:884", function() .py_lt(.py_abs(.py_sub(.py_getitem(.py_attr(b, "deflection_jumps"), w8), .py_div(85L, 24L))), tolerance))
})
.py_try({
.py_rec("beam:test_apply_sliding_hinge:885", function() .py_method(b, "bending_moment"))
})
.py_try({
.py_rec("beam:test_apply_sliding_hinge:887", function() .py_method(b, "deflection"))
})
.py_try({
E <- Symbol("E")
})
.py_try({
I <- Symbol("I")
})
.py_try({
I2 <- Symbol("I2")
})
.py_try({
b1 <- Beam(5L, E, I)
})
.py_try({
b2 <- Beam(8L, E, I2)
})
.py_try({
b <- .py_method(b1, "join", b2)
})
.py_try({
.tmpu_ <- .py_method(b, "apply_support", 0L, type = "fixed"); r0 <- .py_getitem(.tmpu_, 0L); m0 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
.py_method(b, "apply_sliding_hinge", 8L)
})
.py_try({
r13 <- .py_method(b, "apply_support", 13L, type = "pin")
})
.py_try({
.py_method(b, "apply_load", -10L, 5L, -1L)
})
.py_try({
.py_method(b, "solve_for_reaction_loads", r0, m0, r13)
})
.py_try({
W_8 <- Symbol("W_8")
})
.py_try({
.py_rec("beam:test_apply_sliding_hinge:902", function() .py_attr(b, "deflection_jumps"))
})
.py_try({
E <- Symbol("E")
})
.py_try({
I <- Symbol("I")
})
.py_try({
q <- Symbol("q")
})
.py_try({
l1 <- Symbol("l1", positive = TRUE)
})
.py_try({
l2 <- Symbol("l2", positive = TRUE)
})
.py_try({
l3 <- Symbol("l3", positive = TRUE)
})
.py_try({
L <- .py_add(.py_add(l1, l2), l3)
})
.py_try({
b <- Beam(L, E, I)
})
.py_try({
r0 <- .py_method(b, "apply_support", 0L, type = "pin")
})
.py_try({
r3 <- .py_method(b, "apply_support", l1, type = "pin")
})
.py_try({
.py_method(b, "apply_sliding_hinge", .py_add(l1, l2))
})
.py_try({
r10 <- .py_method(b, "apply_support", L, type = "pin")
})
.py_try({
.py_method(b, "apply_load", q, 0L, 0L, l1)
})
.py_try({
.py_method(b, "solve_for_reaction_loads", r0, r3, r10)
})
.py_try({
.py_rec("beam:test_apply_sliding_hinge:918", function() .py_method(b, "bending_moment"))
})
.py_try({
.py_rec("beam:test_apply_sliding_hinge:921", function() .py_method(b, "deflection"))
}))---")
  expect_identical(got[["beam:test_apply_sliding_hinge:884"]], "True")
  expect_identical(got[["beam:test_apply_sliding_hinge:885"]], "50*SingularityFunction(x, 0, 0) - 10*SingularityFunction(x, 0, 1) + 10*SingularityFunction(x, 5, 1) - 4250*SingularityFunction(x, 8, -2)/3")
  expect_identical(got[["beam:test_apply_sliding_hinge:887"]], "-SingularityFunction(x, 0, 2)/16 + SingularityFunction(x, 0, 3)/240 - SingularityFunction(x, 5, 3)/240 + 85*SingularityFunction(x, 8, 0)/24")
  expect_identical(got[["beam:test_apply_sliding_hinge:902"]], "{W_8: 4250/(3*E*I2)}")
  expect_identical(got[["beam:test_apply_sliding_hinge:918"]], "l1*q*SingularityFunction(x, 0, 1)/2 + l1*q*SingularityFunction(x, l1, 1)/2 - q*SingularityFunction(x, 0, 2)/2 + q*SingularityFunction(x, l1, 2)/2 + (-l1^3*l2*q/24 - l1^3*l3*q/24)*SingularityFunction(x, l1 + l2, -2)")
  expect_identical(got[["beam:test_apply_sliding_hinge:921"]], "(l1^3*q*x/24 - l1*q*SingularityFunction(x, 0, 3)/12 - l1*q*SingularityFunction(x, l1, 3)/12 + q*SingularityFunction(x, 0, 4)/24 - q*SingularityFunction(x, l1, 4)/24 + (l1^3*l2*q/24 + l1^3*l3*q/24)*SingularityFunction(x, l1 + l2, 0))/(E*I)")
})

test_that("beam test_max_shear_force matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
x <- Symbol("x")
})
.py_try({
y <- Symbol("y")
})
.py_try({
.tmpu_ <- .py_symbols("R1, R2"); R1 <- .py_getitem(.tmpu_, 0L); R2 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
E <- Symbol("E")
})
.py_try({
I <- Symbol("I")
})
.py_try({
b <- Beam(3L, E, I)
})
.py_try({
.tmpu_ <- .py_symbols("R, M"); R <- .py_getitem(.tmpu_, 0L); M <- .py_getitem(.tmpu_, 1L)
})
.py_try({
.py_method(b, "apply_load", R, 0L, -1L)
})
.py_try({
.py_method(b, "apply_load", M, 0L, -2L)
})
.py_try({
.py_method(b, "apply_load", 2L, 3L, -1L)
})
.py_try({
.py_method(b, "apply_load", 4L, 2L, -1L)
})
.py_try({
.py_method(b, "apply_load", 2L, 2L, 0L, end = 3L)
})
.py_try({
.py_method(b, "solve_for_reaction_loads", R, M)
})
.py_try({
.py_rec("beam:test_max_shear_force:938", function() .py_method(b, "max_shear_force"))
}))---")
  expect_identical(got[["beam:test_max_shear_force:938"]], "(Interval(0, 2), 8)")
})

test_that("beam test_max_shear_force_symbolic_load matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
x <- Symbol("x")
})
.py_try({
y <- Symbol("y")
})
.py_try({
.tmpu_ <- .py_symbols("R1, R2"); R1 <- .py_getitem(.tmpu_, 0L); R2 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
E <- Symbol("E")
})
.py_try({
I <- Symbol("I")
})
.py_try({
l <- .py_symbols("l", positive = TRUE)
})
.py_try({
P <- Symbol("P")
})
.py_try({
b <- Beam(l, E, I)
})
.py_try({
.tmpu_ <- .py_symbols("R1, R2"); R1 <- .py_getitem(.tmpu_, 0L); R2 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
.py_method(b, "apply_load", R1, 0L, -1L)
})
.py_try({
.py_method(b, "apply_load", R2, l, -1L)
})
.py_try({
.py_method(b, "apply_load", P, 0L, 0L, end = l)
})
.py_try({
.py_method(b, "solve_for_reaction_loads", R1, R2)
})
.py_try({
max_shear <- .py_method(b, "max_shear_force")
})
.py_try({
.py_rec("beam:test_max_shear_force_symbolic_load:954", function() .py_getitem(max_shear, 0L))
})
.py_try({
.py_rec("beam:test_max_shear_force_symbolic_load:955", function() simplify(.py_sub(.py_getitem(max_shear, 1L), .py_div(.py_mul(l, Abs(P)), 2L))))
}))---")
  expect_identical(got[["beam:test_max_shear_force_symbolic_load:954"]], "0")
  expect_identical(got[["beam:test_max_shear_force_symbolic_load:955"]], "0")
})

test_that("beam test_max_bmoment matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
x <- Symbol("x")
})
.py_try({
y <- Symbol("y")
})
.py_try({
.tmpu_ <- .py_symbols("R1, R2"); R1 <- .py_getitem(.tmpu_, 0L); R2 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
E <- Symbol("E")
})
.py_try({
I <- Symbol("I")
})
.py_try({
.tmpu_ <- .py_symbols("l, P", positive = TRUE); l <- .py_getitem(.tmpu_, 0L); P <- .py_getitem(.tmpu_, 1L)
})
.py_try({
b <- Beam(l, E, I)
})
.py_try({
.tmpu_ <- .py_symbols("R1, R2"); R1 <- .py_getitem(.tmpu_, 0L); R2 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
.py_method(b, "apply_load", R1, 0L, -1L)
})
.py_try({
.py_method(b, "apply_load", R2, l, -1L)
})
.py_try({
.py_method(b, "apply_load", P, .py_div(l, 2L), -1L)
})
.py_try({
.py_method(b, "solve_for_reaction_loads", R1, R2)
})
.py_try({
.py_attr(b, "reaction_loads")
})
.py_try({
.py_rec("beam:test_max_bmoment:970", function() .py_method(b, "max_bmoment"))
})
.py_try({
b <- Beam(l, E, I)
})
.py_try({
.tmpu_ <- .py_symbols("R1, R2"); R1 <- .py_getitem(.tmpu_, 0L); R2 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
.py_method(b, "apply_load", R1, 0L, -1L)
})
.py_try({
.py_method(b, "apply_load", R2, l, -1L)
})
.py_try({
.py_method(b, "apply_load", P, 0L, 0L, end = l)
})
.py_try({
.py_method(b, "solve_for_reaction_loads", R1, R2)
})
.py_try({
.py_rec("beam:test_max_bmoment:978", function() .py_method(b, "max_bmoment"))
}))---")
  expect_identical(got[["beam:test_max_bmoment:970"]], "(l/2, P*l/4)")
  expect_identical(got[["beam:test_max_bmoment:978"]], "(l/2, P*l^2/8)")
})

test_that("beam test_max_deflection matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
x <- Symbol("x")
})
.py_try({
y <- Symbol("y")
})
.py_try({
.tmpu_ <- .py_symbols("R1, R2"); R1 <- .py_getitem(.tmpu_, 0L); R2 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
.tmpu_ <- .py_symbols("E, I, l, F", positive = TRUE); E <- .py_getitem(.tmpu_, 0L); I <- .py_getitem(.tmpu_, 1L); l <- .py_getitem(.tmpu_, 2L); F <- .py_getitem(.tmpu_, 3L)
})
.py_try({
b <- Beam(l, E, I)
})
.py_try({
.py_setattr(b, "bc_deflection", list(.py_tup(0L, 0L), .py_tup(l, 0L)))
})
.py_try({
.py_setattr(b, "bc_slope", list(.py_tup(0L, 0L), .py_tup(l, 0L)))
})
.py_try({
.py_method(b, "apply_load", .py_div(F, 2L), 0L, -1L)
})
.py_try({
.py_method(b, "apply_load", .py_div(.py_mul(.py_neg(F), l), 8L), 0L, -2L)
})
.py_try({
.py_method(b, "apply_load", .py_div(F, 2L), l, -1L)
})
.py_try({
.py_method(b, "apply_load", .py_div(.py_mul(F, l), 8L), l, -2L)
})
.py_try({
.py_method(b, "apply_load", .py_neg(F), .py_div(l, 2L), -1L)
})
.py_try({
.py_rec("beam:test_max_deflection:991", function() .py_method(b, "max_deflection"))
}))---")
  expect_identical(got[["beam:test_max_deflection:991"]], "(l/2, F*l^3/(192*E*I))")
})

test_that("beam test_solve_for_ild_reactions matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
x <- Symbol("x")
})
.py_try({
y <- Symbol("y")
})
.py_try({
.tmpu_ <- .py_symbols("R1, R2"); R1 <- .py_getitem(.tmpu_, 0L); R2 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
E <- Symbol("E")
})
.py_try({
I <- Symbol("I")
})
.py_try({
b <- Beam(10L, E, I)
})
.py_try({
.py_method(b, "apply_support", 0L, type = "pin")
})
.py_try({
.py_method(b, "apply_support", 10L, type = "pin")
})
.py_try({
.tmpu_ <- .py_symbols("R_0, R_10"); R_0 <- .py_getitem(.tmpu_, 0L); R_10 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
.py_method(b, "solve_for_ild_reactions", 1L, R_0, R_10)
})
.py_try({
a <- .py_attr(b, "ild_variable")
})
.py_try({
.py_rec("beam:test_solve_for_ild_reactions:1002", function() .py_attr(b, "ild_reactions"))
})
.py_try({
E <- Symbol("E")
})
.py_try({
I <- Symbol("I")
})
.py_try({
F <- Symbol("F")
})
.py_try({
L <- Symbol("L", positive = TRUE)
})
.py_try({
b <- Beam(L, E, I)
})
.py_try({
.py_method(b, "apply_support", L, type = "fixed")
})
.py_try({
.py_method(b, "apply_load", F, 0L, -1L)
})
.py_try({
.tmpu_ <- .py_symbols("R_L, M_L"); R_L <- .py_getitem(.tmpu_, 0L); M_L <- .py_getitem(.tmpu_, 1L)
})
.py_try({
.py_method(b, "solve_for_ild_reactions", F, R_L, M_L)
})
.py_try({
a <- .py_attr(b, "ild_variable")
})
.py_try({
.py_rec("beam:test_solve_for_ild_reactions:1017", function() .py_attr(b, "ild_reactions"))
})
.py_try({
E <- Symbol("E")
})
.py_try({
I <- Symbol("I")
})
.py_try({
b <- Beam(20L, E, I)
})
.py_try({
r0 <- .py_method(b, "apply_support", 0L, type = "pin")
})
.py_try({
r5 <- .py_method(b, "apply_support", 5L, type = "pin")
})
.py_try({
r10 <- .py_method(b, "apply_support", 10L, type = "pin")
})
.py_try({
.tmpu_ <- .py_method(b, "apply_support", 20L, type = "fixed"); r20 <- .py_getitem(.tmpu_, 0L); m20 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
.py_method(b, "solve_for_ild_reactions", 1L, r0, r5, r10, r20, m20)
})
.py_try({
a <- .py_attr(b, "ild_variable")
})
.py_try({
.py_rec("beam:test_solve_for_ild_reactions:1030", function() .py_method(.py_getitem(.py_attr(b, "ild_reactions"), r0), "subs", a, 4L))
})
.py_try({
.py_rec("beam:test_solve_for_ild_reactions:1031", function() .py_method(.py_getitem(.py_attr(b, "ild_reactions"), r5), "subs", a, 4L))
})
.py_try({
.py_rec("beam:test_solve_for_ild_reactions:1032", function() .py_method(.py_getitem(.py_attr(b, "ild_reactions"), r10), "subs", a, 4L))
})
.py_try({
.py_rec("beam:test_solve_for_ild_reactions:1033", function() .py_method(.py_getitem(.py_attr(b, "ild_reactions"), r20), "subs", a, 12L))
})
.py_try({
.py_rec("beam:test_solve_for_ild_reactions:1034", function() .py_method(.py_getitem(.py_attr(b, "ild_reactions"), m20), "subs", a, 12L))
}))---")
  expect_identical(got[["beam:test_solve_for_ild_reactions:1002"]], "{R_0: -SingularityFunction(a, 0, 0) + SingularityFunction(a, 0, 1)/10 - SingularityFunction(a, 10, 1)/10, R_10: -SingularityFunction(a, 0, 1)/10 + SingularityFunction(a, 10, 0) + SingularityFunction(a, 10, 1)/10}")
  expect_identical(got[["beam:test_solve_for_ild_reactions:1017"]], "{R_L: -F*SingularityFunction(a, 0, 0) + F*SingularityFunction(a, L, 0) - F, M_L: -F*L*SingularityFunction(a, 0, 0) - F*L + F*SingularityFunction(a, 0, 1) - F*SingularityFunction(a, L, 1)}")
  expect_identical(got[["beam:test_solve_for_ild_reactions:1030"]], "-59/475")
  expect_identical(got[["beam:test_solve_for_ild_reactions:1031"]], "-2296/2375")
  expect_identical(got[["beam:test_solve_for_ild_reactions:1032"]], "243/2375")
  expect_identical(got[["beam:test_solve_for_ild_reactions:1033"]], "-83/475")
  expect_identical(got[["beam:test_solve_for_ild_reactions:1034"]], "-264/475")
})

test_that("beam test_solve_for_ild_shear matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
x <- Symbol("x")
})
.py_try({
y <- Symbol("y")
})
.py_try({
.tmpu_ <- .py_symbols("R1, R2"); R1 <- .py_getitem(.tmpu_, 0L); R2 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
E <- Symbol("E")
})
.py_try({
I <- Symbol("I")
})
.py_try({
F <- Symbol("F")
})
.py_try({
L1 <- Symbol("L1", positive = TRUE)
})
.py_try({
L2 <- Symbol("L2", positive = TRUE)
})
.py_try({
b <- Beam(.py_add(L1, L2), E, I)
})
.py_try({
r0 <- .py_method(b, "apply_support", 0L, type = "pin")
})
.py_try({
rL <- .py_method(b, "apply_support", .py_add(L1, L2), type = "pin")
})
.py_try({
.py_method(b, "solve_for_ild_reactions", F, r0, rL)
})
.py_try({
.py_method(b, "solve_for_ild_shear", L1, F, r0, rL)
})
.py_try({
a <- .py_attr(b, "ild_variable")
})
.py_try({
expected_shear <- .py_sub(.py_sub(.py_add(.py_add(.py_add(.py_sub(.py_sub(.py_div(.py_mul(.py_mul(.py_neg(F), L1), SingularityFunction(a, 0L, 0L)), .py_add(L1, L2)), .py_div(.py_mul(.py_mul(F, L2), SingularityFunction(a, 0L, 0L)), .py_add(L1, L2))), .py_mul(F, SingularityFunction(.py_neg(a), 0L, 0L))), .py_mul(F, SingularityFunction(a, .py_add(L1, L2), 0L))), F), .py_div(.py_mul(F, SingularityFunction(a, 0L, 1L)), .py_add(L1, L2))), .py_div(.py_mul(F, SingularityFunction(a, .py_add(L1, L2), 1L)), .py_add(L1, L2))), .py_mul(.py_add(.py_add(.py_sub(.py_add(.py_div(.py_mul(.py_mul(.py_neg(F), L1), SingularityFunction(a, 0L, 0L)), .py_add(L1, L2)), .py_div(.py_mul(.py_mul(F, L1), SingularityFunction(a, .py_add(L1, L2), 0L)), .py_add(L1, L2))), .py_div(.py_mul(.py_mul(F, L2), SingularityFunction(a, 0L, 0L)), .py_add(L1, L2))), .py_div(.py_mul(.py_mul(F, L2), SingularityFunction(a, .py_add(L1, L2), 0L)), .py_add(L1, L2))), .py_mul(2L, F)), SingularityFunction(a, L1, 0L)))
})
.py_try({
.py_rec("beam:test_solve_for_ild_shear:1054", function() .py_method(.py_attr(b, "ild_shear"), "expand"))
})
.py_try({
E <- Symbol("E")
})
.py_try({
I <- Symbol("I")
})
.py_try({
b <- Beam(20L, E, I)
})
.py_try({
r0 <- .py_method(b, "apply_support", 0L, type = "pin")
})
.py_try({
r5 <- .py_method(b, "apply_support", 5L, type = "pin")
})
.py_try({
r10 <- .py_method(b, "apply_support", 10L, type = "pin")
})
.py_try({
.tmpu_ <- .py_method(b, "apply_support", 20L, type = "fixed"); r20 <- .py_getitem(.tmpu_, 0L); m20 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
.py_method(b, "solve_for_ild_reactions", 1L, r0, r5, r10, r20, m20)
})
.py_try({
.py_method(b, "solve_for_ild_shear", 6L, 1L, r0, r5, r10, r20, m20)
})
.py_try({
a <- .py_attr(b, "ild_variable")
})
.py_try({
.py_rec("beam:test_solve_for_ild_shear:1066", function() .py_method(.py_attr(b, "ild_shear"), "subs", a, 12L))
})
.py_try({
.py_rec("beam:test_solve_for_ild_shear:1067", function() .py_method(.py_attr(b, "ild_shear"), "subs", a, 4L))
}))---")
  expect_identical(got[["beam:test_solve_for_ild_shear:1054"]], "F*L1*SingularityFunction(a, 0, 0)*SingularityFunction(a, L1, 0)/(L1 + L2) - F*L1*SingularityFunction(a, 0, 0)/(L1 + L2) - F*L1*SingularityFunction(a, L1, 0)*SingularityFunction(a, L1 + L2, 0)/(L1 + L2) + F*L2*SingularityFunction(a, 0, 0)*SingularityFunction(a, L1, 0)/(L1 + L2) - F*L2*SingularityFunction(a, 0, 0)/(L1 + L2) - F*L2*SingularityFunction(a, L1, 0)*SingularityFunction(a, L1 + L2, 0)/(L1 + L2) - F*SingularityFunction(-a, 0, 0) - 2*F*SingularityFunction(a, L1, 0) + F*SingularityFunction(a, L1 + L2, 0) + F + F*SingularityFunction(a, 0, 1)/(L1 + L2) - F*SingularityFunction(a, L1 + L2, 1)/(L1 + L2)")
  expect_identical(got[["beam:test_solve_for_ild_shear:1066"]], "96/475")
  expect_identical(got[["beam:test_solve_for_ild_shear:1067"]], "-216/2375")
})

test_that("beam test_solve_for_ild_moment matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
x <- Symbol("x")
})
.py_try({
y <- Symbol("y")
})
.py_try({
.tmpu_ <- .py_symbols("R1, R2"); R1 <- .py_getitem(.tmpu_, 0L); R2 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
E <- Symbol("E")
})
.py_try({
I <- Symbol("I")
})
.py_try({
F <- Symbol("F")
})
.py_try({
L1 <- Symbol("L1", positive = TRUE)
})
.py_try({
L2 <- Symbol("L2", positive = TRUE)
})
.py_try({
b <- Beam(.py_add(L1, L2), E, I)
})
.py_try({
r0 <- .py_method(b, "apply_support", 0L, type = "pin")
})
.py_try({
rL <- .py_method(b, "apply_support", .py_add(L1, L2), type = "pin")
})
.py_try({
a <- .py_attr(b, "ild_variable")
})
.py_try({
.py_method(b, "solve_for_ild_reactions", F, r0, rL)
})
.py_try({
.py_method(b, "solve_for_ild_moment", L1, F, r0, rL)
})
.py_try({
.py_rec("beam:test_solve_for_ild_moment:1081", function() .py_method(.py_method(.py_method(.py_attr(b, "ild_moment"), "subs", a, 3L), "subs", L1, 5L), "subs", L2, 5L))
})
.py_try({
E <- Symbol("E")
})
.py_try({
I <- Symbol("I")
})
.py_try({
b <- Beam(20L, E, I)
})
.py_try({
r0 <- .py_method(b, "apply_support", 0L, type = "pin")
})
.py_try({
r5 <- .py_method(b, "apply_support", 5L, type = "pin")
})
.py_try({
r10 <- .py_method(b, "apply_support", 10L, type = "pin")
})
.py_try({
.tmpu_ <- .py_method(b, "apply_support", 20L, type = "fixed"); r20 <- .py_getitem(.tmpu_, 0L); m20 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
.py_method(b, "solve_for_ild_reactions", 1L, r0, r5, r10, r20, m20)
})
.py_try({
.py_method(b, "solve_for_ild_moment", 5L, 1L, r0, r5, r10, r20, m20)
})
.py_try({
.py_rec("beam:test_solve_for_ild_moment:1092", function() .py_method(.py_attr(b, "ild_moment"), "subs", a, 12L))
})
.py_try({
.py_rec("beam:test_solve_for_ild_moment:1093", function() .py_method(.py_attr(b, "ild_moment"), "subs", a, 4L))
}))---")
  expect_identical(got[["beam:test_solve_for_ild_moment:1081"]], "-3*F/2")
  expect_identical(got[["beam:test_solve_for_ild_moment:1092"]], "-96/475")
  expect_identical(got[["beam:test_solve_for_ild_moment:1093"]], "36/95")
})

test_that("beam test_ild_with_rotation_hinge matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
x <- Symbol("x")
})
.py_try({
y <- Symbol("y")
})
.py_try({
.tmpu_ <- .py_symbols("R1, R2"); R1 <- .py_getitem(.tmpu_, 0L); R2 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
E <- Symbol("E")
})
.py_try({
I <- Symbol("I")
})
.py_try({
F <- Symbol("F")
})
.py_try({
L1 <- Symbol("L1", positive = TRUE)
})
.py_try({
L2 <- Symbol("L2", positive = TRUE)
})
.py_try({
L3 <- Symbol("L3", positive = TRUE)
})
.py_try({
b <- Beam(.py_add(.py_add(L1, L2), L3), E, I)
})
.py_try({
r0 <- .py_method(b, "apply_support", 0L, type = "pin")
})
.py_try({
r1 <- .py_method(b, "apply_support", .py_add(L1, L2), type = "pin")
})
.py_try({
r2 <- .py_method(b, "apply_support", .py_add(.py_add(L1, L2), L3), type = "pin")
})
.py_try({
.py_method(b, "apply_rotation_hinge", .py_add(L1, L2))
})
.py_try({
.py_method(b, "solve_for_ild_reactions", F, r0, r1, r2)
})
.py_try({
a <- .py_attr(b, "ild_variable")
})
.py_try({
.py_rec("beam:test_ild_with_rotation_hinge:1109", function() .py_method(.py_method(.py_method(.py_method(.py_getitem(.py_attr(b, "ild_reactions"), r0), "subs", a, 4L), "subs", L1, 5L), "subs", L2, 5L), "subs", L3, 10L))
})
.py_try({
.py_rec("beam:test_ild_with_rotation_hinge:1110", function() .py_method(.py_method(.py_method(.py_method(.py_getitem(.py_attr(b, "ild_reactions"), r0), "subs", a, -10L), "subs", L1, 5L), "subs", L2, 5L), "subs", L3, 10L))
})
.py_try({
.py_rec("beam:test_ild_with_rotation_hinge:1111", function() .py_method(.py_method(.py_method(.py_method(.py_getitem(.py_attr(b, "ild_reactions"), r0), "subs", a, 25L), "subs", L1, 5L), "subs", L2, 5L), "subs", L3, 10L))
})
.py_try({
.py_rec("beam:test_ild_with_rotation_hinge:1112", function() .py_method(.py_method(.py_method(.py_method(.py_getitem(.py_attr(b, "ild_reactions"), r1), "subs", a, 4L), "subs", L1, 5L), "subs", L2, 5L), "subs", L3, 10L))
})
.py_try({
.py_rec("beam:test_ild_with_rotation_hinge:1113", function() .py_method(.py_method(.py_method(.py_method(.py_getitem(.py_attr(b, "ild_reactions"), r2), "subs", a, 18L), "subs", L1, 5L), "subs", L2, 5L), "subs", L3, 10L))
})
.py_try({
.py_method(b, "solve_for_ild_shear", L1, F, r0, r1, r2)
})
.py_try({
.py_rec("beam:test_ild_with_rotation_hinge:1115", function() .py_method(.py_method(.py_method(.py_method(.py_attr(b, "ild_shear"), "subs", a, 7L), "subs", L1, 5L), "subs", L2, 5L), "subs", L3, 10L))
})
.py_try({
.py_rec("beam:test_ild_with_rotation_hinge:1116", function() .py_method(.py_method(.py_method(.py_method(.py_attr(b, "ild_shear"), "subs", a, 70L), "subs", L1, 5L), "subs", L2, 5L), "subs", L3, 10L))
})
.py_try({
.py_method(b, "solve_for_ild_moment", L1, F, r0, r1, r2)
})
.py_try({
.py_rec("beam:test_ild_with_rotation_hinge:1118", function() .py_method(.py_method(.py_method(.py_method(.py_attr(b, "ild_moment"), "subs", a, 1L), "subs", L1, 5L), "subs", L2, 5L), "subs", L3, 10L))
})
.py_try({
.py_rec("beam:test_ild_with_rotation_hinge:1119", function() .py_method(.py_method(.py_method(.py_method(.py_attr(b, "ild_moment"), "subs", a, 8L), "subs", L1, 5L), "subs", L2, 5L), "subs", L3, 10L))
}))---")
  expect_identical(got[["beam:test_ild_with_rotation_hinge:1109"]], "-3*F/5")
  expect_identical(got[["beam:test_ild_with_rotation_hinge:1110"]], "0")
  expect_identical(got[["beam:test_ild_with_rotation_hinge:1111"]], "0")
  expect_identical(got[["beam:test_ild_with_rotation_hinge:1112"]], "-2*F/5")
  expect_identical(got[["beam:test_ild_with_rotation_hinge:1113"]], "-4*F/5")
  expect_identical(got[["beam:test_ild_with_rotation_hinge:1115"]], "-3*F/10")
  expect_identical(got[["beam:test_ild_with_rotation_hinge:1116"]], "0")
  expect_identical(got[["beam:test_ild_with_rotation_hinge:1118"]], "-F/2")
  expect_identical(got[["beam:test_ild_with_rotation_hinge:1119"]], "-F")
})

test_that("beam test_ild_with_sliding_hinge matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
x <- Symbol("x")
})
.py_try({
y <- Symbol("y")
})
.py_try({
.tmpu_ <- .py_symbols("R1, R2"); R1 <- .py_getitem(.tmpu_, 0L); R2 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
b <- Beam(13L, 200L, 200L)
})
.py_try({
r0 <- .py_method(b, "apply_support", 0L, type = "pin")
})
.py_try({
r6 <- .py_method(b, "apply_support", 6L, type = "pin")
})
.py_try({
.tmpu_ <- .py_method(b, "apply_support", 13L, type = "fixed"); r13 <- .py_getitem(.tmpu_, 0L); m13 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
w3 <- .py_method(b, "apply_sliding_hinge", 3L)
})
.py_try({
.py_method(b, "solve_for_ild_reactions", 1L, r0, r6, r13, m13)
})
.py_try({
a <- .py_attr(b, "ild_variable")
})
.py_try({
.py_rec("beam:test_ild_with_sliding_hinge:1129", function() .py_method(.py_getitem(.py_attr(b, "ild_reactions"), r0), "subs", a, 3L))
})
.py_try({
.py_rec("beam:test_ild_with_sliding_hinge:1130", function() .py_method(.py_getitem(.py_attr(b, "ild_reactions"), r6), "subs", a, 3L))
})
.py_try({
.py_rec("beam:test_ild_with_sliding_hinge:1131", function() .py_method(.py_getitem(.py_attr(b, "ild_reactions"), r13), "subs", a, 9L))
})
.py_try({
.py_rec("beam:test_ild_with_sliding_hinge:1132", function() .py_method(.py_getitem(.py_attr(b, "ild_reactions"), m13), "subs", a, 9L))
})
.py_try({
.py_rec("beam:test_ild_with_sliding_hinge:1133", function() .py_method(.py_getitem(.py_attr(b, "ild_reactions"), m13), "subs", a, 15L))
})
.py_try({
.py_rec("beam:test_ild_with_sliding_hinge:1134", function() .py_method(.py_getitem(.py_attr(b, "ild_reactions"), m13), "subs", a, -3L))
})
.py_try({
.py_rec("beam:test_ild_with_sliding_hinge:1135", function() .py_method(.py_getitem(.py_attr(b, "ild_deflection_jumps"), w3), "subs", a, 9L))
})
.py_try({
.py_method(b, "solve_for_ild_shear", 7L, 1L, r0, r6, r13, m13)
})
.py_try({
.py_rec("beam:test_ild_with_sliding_hinge:1137", function() .py_method(.py_attr(b, "ild_shear"), "subs", a, 8L))
})
.py_try({
.py_method(b, "solve_for_ild_moment", 8L, 1L, r0, r6, r13, m13)
})
.py_try({
.py_rec("beam:test_ild_with_sliding_hinge:1139", function() .py_method(.py_attr(b, "ild_moment"), "subs", a, 3L))
}))---")
  expect_identical(got[["beam:test_ild_with_sliding_hinge:1129"]], "-1")
  expect_identical(got[["beam:test_ild_with_sliding_hinge:1130"]], "9/14")
  expect_identical(got[["beam:test_ild_with_sliding_hinge:1131"]], "-207/343")
  expect_identical(got[["beam:test_ild_with_sliding_hinge:1132"]], "-60/49")
  expect_identical(got[["beam:test_ild_with_sliding_hinge:1133"]], "0")
  expect_identical(got[["beam:test_ild_with_sliding_hinge:1134"]], "0")
  expect_identical(got[["beam:test_ild_with_sliding_hinge:1135"]], "-9/35000")
  expect_identical(got[["beam:test_ild_with_sliding_hinge:1137"]], "-200/343")
  expect_identical(got[["beam:test_ild_with_sliding_hinge:1139"]], "-12/7")
})

test_that("beam test_apply_load_ramp_start_greater_than_end matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
x <- Symbol("x")
})
.py_try({
y <- Symbol("y")
})
.py_try({
.tmpu_ <- .py_symbols("R1, R2"); R1 <- .py_getitem(.tmpu_, 0L); R2 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
.tmpu_ <- .py_symbols("E I"); E <- .py_getitem(.tmpu_, 0L); I <- .py_getitem(.tmpu_, 1L)
})
.py_try({
b <- Beam(10L, E, I)
})
.py_try({
.py_method(b, "apply_load", 5L, 10L, 1L, end = 0L)
})
.py_try({
expected_load <- .py_add(.py_sub(.py_mul(50L, SingularityFunction(x, 0L, 0L)), .py_mul(5L, SingularityFunction(x, 0L, 1L))), .py_mul(5L, SingularityFunction(x, 10L, 1L)))
})
.py_try({
.py_rec("beam:test_apply_load_ramp_start_greater_than_end:1153", function() simplify(.py_sub(.py_attr(b, "load"), expected_load)))
})
.py_try({
b1 <- Beam(10L, E, I)
})
.py_try({
.py_method(b1, "apply_load", 5L, 10L, 3L, end = 0L)
})
.py_try({
expected_load <- .py_add(.py_sub(.py_add(.py_sub(.py_mul(5000L, SingularityFunction(x, 0L, 0L)), .py_mul(1500L, SingularityFunction(x, 0L, 1L))), .py_mul(150L, SingularityFunction(x, 0L, 2L))), .py_mul(5L, SingularityFunction(x, 0L, 3L))), .py_mul(5L, SingularityFunction(x, 10L, 3L)))
})
.py_try({
.py_rec("beam:test_apply_load_ramp_start_greater_than_end:1163", function() simplify(.py_sub(.py_attr(b1, "load"), expected_load)))
})
.py_try({
b2 <- Beam(10L, E, I)
})
.py_try({
.py_method(b2, "apply_load", 5L, 6L, 1L, end = 2L)
})
.py_try({
expected_load <- .py_add(.py_sub(.py_mul(20L, SingularityFunction(x, 2L, 0L)), .py_mul(5L, SingularityFunction(x, 2L, 1L))), .py_mul(5L, SingularityFunction(x, 6L, 1L)))
})
.py_try({
.py_rec("beam:test_apply_load_ramp_start_greater_than_end:1171", function() simplify(.py_sub(.py_attr(b2, "load"), expected_load)))
}))---")
  expect_identical(got[["beam:test_apply_load_ramp_start_greater_than_end:1153"]], "0")
  expect_identical(got[["beam:test_apply_load_ramp_start_greater_than_end:1163"]], "0")
  expect_identical(got[["beam:test_apply_load_ramp_start_greater_than_end:1171"]], "0")
})

test_that("beam test_apply_load_symbolic_start_end matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
x <- Symbol("x")
})
.py_try({
y <- Symbol("y")
})
.py_try({
.tmpu_ <- .py_symbols("R1, R2"); R1 <- .py_getitem(.tmpu_, 0L); R2 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
.tmpu_ <- .py_symbols("E I"); E <- .py_getitem(.tmpu_, 0L); I <- .py_getitem(.tmpu_, 1L)
})
.py_try({
.tmpu_ <- .py_symbols("a b", real = TRUE); a <- .py_getitem(.tmpu_, 0L); b_ <- .py_getitem(.tmpu_, 1L)
})
.py_try({
beam <- Beam(b_, E, I)
})
.py_try({
.py_method(beam, "apply_load", 5L, b_, 1L, end = a)
})
.py_try({
f <- .py_mul(5L, .py_pow(x, 1L))
})
.py_try({
expected_load <- .py_add(.py_sub(.py_mul(.py_method(f, "subs", x, .py_sub(b_, a)), SingularityFunction(x, a, 0L)), .py_div(.py_mul(.py_method(.py_method(f, "diff", x, 1L), "subs", x, .py_sub(b_, a)), SingularityFunction(x, a, 1L)), sym_factorial(1L))), .py_mul(5L, SingularityFunction(x, b_, 1L)))
})
.py_try({
.py_rec("beam:test_apply_load_symbolic_start_end:1189", function() simplify(.py_sub(.py_attr(beam, "load"), expected_load)))
}))---")
  expect_identical(got[["beam:test_apply_load_symbolic_start_end:1189"]], "0")
})

test_that("beam test_Beam3D matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
x <- Symbol("x")
})
.py_try({
y <- Symbol("y")
})
.py_try({
.tmpu_ <- .py_symbols("R1, R2"); R1 <- .py_getitem(.tmpu_, 0L); R2 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
.tmpu_ <- .py_symbols("l, E, G, I, A"); l <- .py_getitem(.tmpu_, 0L); E <- .py_getitem(.tmpu_, 1L); G <- .py_getitem(.tmpu_, 2L); I <- .py_getitem(.tmpu_, 3L); A <- .py_getitem(.tmpu_, 4L)
})
.py_try({
.tmpu_ <- .py_symbols("R1, R2, R3, R4"); R1 <- .py_getitem(.tmpu_, 0L); R2 <- .py_getitem(.tmpu_, 1L); R3 <- .py_getitem(.tmpu_, 2L); R4 <- .py_getitem(.tmpu_, 3L)
})
.py_try({
b <- Beam3D(l, E, G, I, A)
})
.py_try({
.tmpu_ <- .py_symbols("m, q"); m <- .py_getitem(.tmpu_, 0L); q <- .py_getitem(.tmpu_, 1L)
})
.py_try({
.py_method(b, "apply_load", q, 0L, 0L, dir = "y")
})
.py_try({
.py_method(b, "apply_moment_load", m, 0L, 0L, dir = "z")
})
.py_try({
.py_setattr(b, "bc_slope", list(.py_tup(0L, list(0L, 0L, 0L)), .py_tup(l, list(0L, 0L, 0L))))
})
.py_try({
.py_setattr(b, "bc_deflection", list(.py_tup(0L, list(0L, 0L, 0L)), .py_tup(l, list(0L, 0L, 0L))))
})
.py_try({
.py_method(b, "solve_slope_deflection")
})
.py_try({
.py_rec("beam:test_Beam3D:1204", function() .py_method(b, "polar_moment"))
})
.py_try({
.py_rec("beam:test_Beam3D:1205", function() .py_method(b, "shear_force"))
})
.py_try({
.py_rec("beam:test_Beam3D:1206", function() .py_method(b, "shear_stress"))
})
.py_try({
.py_rec("beam:test_Beam3D:1207", function() .py_method(b, "axial_stress"))
})
.py_try({
.py_rec("beam:test_Beam3D:1208", function() .py_method(b, "bending_moment"))
})
.py_try({
expected_deflection <- .py_div(.py_mul(x, .py_add(.py_add(.py_add(.py_div(.py_mul(.py_mul(.py_mul(A, G), q), .py_pow(x, 3L)), 4L), .py_mul(.py_mul(.py_mul(A, G), .py_pow(x, 2L)), .py_sub(.py_div(.py_div(.py_mul(.py_neg(l), .py_add(.py_mul(.py_mul(.py_mul(A, G), l), .py_sub(.py_mul(l, q), .py_mul(2L, m))), .py_mul(.py_mul(.py_mul(12L, E), I), q))), .py_add(.py_mul(.py_mul(A, G), .py_pow(l, 2L)), .py_mul(.py_mul(12L, E), I))), 2L), m))), .py_div(.py_mul(.py_mul(.py_mul(.py_mul(3L, E), I), l), .py_add(.py_mul(.py_mul(.py_mul(A, G), l), .py_sub(.py_mul(l, q), .py_mul(2L, m))), .py_mul(.py_mul(.py_mul(12L, E), I), q))), .py_add(.py_mul(.py_mul(A, G), .py_pow(l, 2L)), .py_mul(.py_mul(12L, E), I)))), .py_mul(x, .py_sub(.py_add(.py_add(.py_div(.py_mul(.py_mul(.py_mul(.py_neg(A), G), .py_pow(l, 2L)), q), 2L), .py_div(.py_div(.py_mul(.py_mul(.py_mul(.py_mul(3L, A), G), .py_pow(l, 2L)), .py_add(.py_mul(.py_mul(.py_mul(A, G), l), .py_sub(.py_mul(l, q), .py_mul(2L, m))), .py_mul(.py_mul(.py_mul(12L, E), I), q))), .py_add(.py_mul(.py_mul(A, G), .py_pow(l, 2L)), .py_mul(.py_mul(12L, E), I))), 4L)), .py_mul(.py_mul(.py_mul(.py_mul(A, G), l), m), Rational(3L, 2L))), .py_mul(.py_mul(.py_mul(3L, E), I), q))))), .py_mul(.py_mul(.py_mul(.py_mul(6L, A), E), G), I))
})
.py_try({
.tmpu_ <- .py_method(b, "deflection"); dx <- .py_getitem(.tmpu_, 0L); dy <- .py_getitem(.tmpu_, 1L); dz <- .py_getitem(.tmpu_, 2L)
})
.py_try({
.py_rec("beam:test_Beam3D:1215", function() (.py_truthy(.py_eq(dx, dz)) && .py_truthy(.py_eq(dz, 0L))))
})
.py_try({
.py_rec("beam:test_Beam3D:1216", function() simplify(.py_sub(dy, expected_deflection)))
})
.py_try({
b2 <- Beam3D(30L, E, G, I, A, x)
})
.py_try({
.py_method(b2, "apply_load", 50L, start = 0L, order = 0L, dir = "y")
})
.py_try({
.py_setattr(b2, "bc_deflection", list(.py_tup(0L, list(0L, 0L, 0L)), .py_tup(30L, list(0L, 0L, 0L))))
})
.py_try({
.py_method(b2, "apply_load", R1, start = 0L, order = -1L, dir = "y")
})
.py_try({
.py_method(b2, "apply_load", R2, start = 30L, order = -1L, dir = "y")
})
.py_try({
.py_method(b2, "solve_for_reaction_loads", R1, R2)
})
.py_try({
.py_rec("beam:test_Beam3D:1224", function() .py_attr(b2, "reaction_loads"))
})
.py_try({
.py_method(b2, "solve_slope_deflection")
})
.py_try({
.py_rec("beam:test_Beam3D:1227", function() .py_method(b2, "slope"))
})
.py_try({
expected_deflection <- .py_add(.py_sub(.py_add(.py_sub(.py_div(.py_mul(25L, .py_pow(x, 4L)), .py_mul(.py_mul(12L, E), I)), .py_div(.py_mul(125L, .py_pow(x, 3L)), .py_mul(E, I))), .py_div(.py_mul(1875L, .py_pow(x, 2L)), .py_mul(E, I))), .py_div(.py_mul(25L, .py_pow(x, 2L)), .py_mul(A, G))), .py_div(.py_mul(750L, x), .py_mul(A, G)))
})
.py_try({
.tmpu_ <- .py_method(b2, "deflection"); dx <- .py_getitem(.tmpu_, 0L); dy <- .py_getitem(.tmpu_, 1L); dz <- .py_getitem(.tmpu_, 2L)
})
.py_try({
.py_rec("beam:test_Beam3D:1231", function() (.py_truthy(.py_eq(dx, dz)) && .py_truthy(.py_eq(dz, 0L))))
})
.py_try({
.py_rec("beam:test_Beam3D:1232", function() dy)
})
.py_try({
b3 <- Beam3D(30L, E, G, I, A, x)
})
.py_try({
.py_method(b3, "apply_load", 8L, start = 0L, order = 0L, dir = "y")
})
.py_try({
.py_method(b3, "apply_load", .py_mul(9L, x), start = 0L, order = 0L, dir = "z")
})
.py_try({
.py_method(b3, "apply_load", R1, start = 0L, order = -1L, dir = "y")
})
.py_try({
.py_method(b3, "apply_load", R2, start = 30L, order = -1L, dir = "y")
})
.py_try({
.py_method(b3, "apply_load", R3, start = 0L, order = -1L, dir = "z")
})
.py_try({
.py_method(b3, "apply_load", R4, start = 30L, order = -1L, dir = "z")
})
.py_try({
.py_method(b3, "solve_for_reaction_loads", R1, R2, R3, R4)
})
.py_try({
.py_rec("beam:test_Beam3D:1243", function() .py_attr(b3, "reaction_loads"))
}))---")
  expect_identical(got[["beam:test_Beam3D:1204"]], "2*I")
  expect_identical(got[["beam:test_Beam3D:1205"]], "[0, -q*x, 0]")
  expect_identical(got[["beam:test_Beam3D:1206"]], "[0, -q*x/A, 0]")
  expect_identical(got[["beam:test_Beam3D:1207"]], "0")
  expect_identical(got[["beam:test_Beam3D:1208"]], "[0, 0, -m*x + q*x^2/2]")
  expect_identical(got[["beam:test_Beam3D:1215"]], "True")
  expect_identical(got[["beam:test_Beam3D:1216"]], "0")
  expect_identical(got[["beam:test_Beam3D:1224"]], "{R1: -750, R2: -750}")
  expect_identical(got[["beam:test_Beam3D:1227"]], "[0, 0, 25*x^3/(3*E*I) - 375*x^2/(E*I) + 3750*x/(E*I)]")
  expect_identical(got[["beam:test_Beam3D:1231"]], "True")
  expect_identical(got[["beam:test_Beam3D:1232"]], "25*x^4/(12*E*I) - 125*x^3/(E*I) + 1875*x^2/(E*I) - 25*x^2/(A*G) + 750*x/(A*G)")
  expect_identical(got[["beam:test_Beam3D:1243"]], "{R1: -120, R2: -120, R3: -1350, R4: -2700}")
})

test_that("beam test_polar_moment_Beam3D matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
x <- Symbol("x")
})
.py_try({
y <- Symbol("y")
})
.py_try({
.tmpu_ <- .py_symbols("R1, R2"); R1 <- .py_getitem(.tmpu_, 0L); R2 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
.tmpu_ <- .py_symbols("l, E, G, A, I1, I2"); l <- .py_getitem(.tmpu_, 0L); E <- .py_getitem(.tmpu_, 1L); G <- .py_getitem(.tmpu_, 2L); A <- .py_getitem(.tmpu_, 3L); I1 <- .py_getitem(.tmpu_, 4L); I2 <- .py_getitem(.tmpu_, 5L)
})
.py_try({
I <- list(I1, I2)
})
.py_try({
b <- Beam3D(l, E, G, I, A)
})
.py_try({
.py_rec("beam:test_polar_moment_Beam3D:1251", function() .py_method(b, "polar_moment"))
}))---")
  expect_identical(got[["beam:test_polar_moment_Beam3D:1251"]], "I1 + I2")
})

test_that("beam test_parabolic_loads matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
x <- Symbol("x")
})
.py_try({
y <- Symbol("y")
})
.py_try({
.tmpu_ <- .py_symbols("R1, R2"); R1 <- .py_getitem(.tmpu_, 0L); R2 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
.tmpu_ <- .py_symbols("E, I, L", positive = TRUE, real = TRUE); E <- .py_getitem(.tmpu_, 0L); I <- .py_getitem(.tmpu_, 1L); L <- .py_getitem(.tmpu_, 2L)
})
.py_try({
.tmpu_ <- .py_symbols("R, M, P", real = TRUE); R <- .py_getitem(.tmpu_, 0L); M <- .py_getitem(.tmpu_, 1L); P <- .py_getitem(.tmpu_, 2L)
})
.py_try({
beam <- Beam(L, E, I)
})
.py_try({
.py_setattr(beam, "bc_deflection", .py_list_mut(.py_attr(beam, "bc_deflection"), "append", .py_tup(0L, 0L)))
})
.py_try({
.py_setattr(beam, "bc_slope", .py_list_mut(.py_attr(beam, "bc_slope"), "append", .py_tup(0L, 0L)))
})
.py_try({
.py_method(beam, "apply_load", R, 0L, -1L)
})
.py_try({
.py_method(beam, "apply_load", M, 0L, -2L)
})
.py_try({
.py_method(beam, "apply_load", 1L, 0L, 2L)
})
.py_try({
.py_method(beam, "solve_for_reaction_loads", R, M)
})
.py_try({
.py_rec("beam:test_parabolic_loads:1273", function() .py_getitem(.py_attr(beam, "reaction_loads"), R))
})
.py_try({
beam <- Beam(.py_mul(2L, L), E, I)
})
.py_try({
.py_setattr(beam, "bc_deflection", .py_list_mut(.py_attr(beam, "bc_deflection"), "append", .py_tup(0L, 0L)))
})
.py_try({
.py_setattr(beam, "bc_slope", .py_list_mut(.py_attr(beam, "bc_slope"), "append", .py_tup(0L, 0L)))
})
.py_try({
.py_method(beam, "apply_load", R, 0L, -1L)
})
.py_try({
.py_method(beam, "apply_load", M, 0L, -2L)
})
.py_try({
.py_method(beam, "apply_load", 1L, 0L, 2L, end = L)
})
.py_try({
.py_method(beam, "solve_for_reaction_loads", R, M)
})
.py_try({
.py_rec("beam:test_parabolic_loads:1290", function() .py_getitem(.py_attr(beam, "reaction_loads"), R))
})
.py_try({
beam <- Beam(.py_mul(2L, L), E, I)
})
.py_try({
.py_method(beam, "apply_load", P, 0L, 0L, end = L)
})
.py_try({
loading <- .py_method(.py_attr(beam, "load"), "xreplace", .py_dict(list(L, E, I, P), list(10L, 20L, 30L, 40L)))
})
.py_try({
.py_rec("beam:test_parabolic_loads:1296", function() .py_method(loading, "xreplace", .py_dict(list(x), list(5L))))
})
.py_try({
.py_rec("beam:test_parabolic_loads:1297", function() .py_method(loading, "xreplace", .py_dict(list(x), list(15L))))
})
.py_try({
beam <- Beam(.py_mul(2L, L), E, I)
})
.py_try({
.py_method(beam, "apply_load", P, 0L, 1L, end = L)
})
.py_try({
.py_rec("beam:test_parabolic_loads:1302", function() .py_attr(beam, "load"))
})
.py_try({
beam <- Beam(.py_mul(2L, L), E, I)
})
.py_try({
.py_method(beam, "apply_load", P, 0L, 8L, end = L)
})
.py_try({
loading <- .py_method(.py_attr(beam, "load"), "xreplace", .py_dict(list(L, E, I, P), list(10L, 20L, 30L, 40L)))
})
.py_try({
.py_rec("beam:test_parabolic_loads:1310", function() .py_method(loading, "xreplace", .py_dict(list(x), list(5L))))
})
.py_try({
.py_rec("beam:test_parabolic_loads:1311", function() .py_method(loading, "xreplace", .py_dict(list(x), list(15L))))
}))---")
  expect_identical(got[["beam:test_parabolic_loads:1273"]], "-L^3/3")
  expect_identical(got[["beam:test_parabolic_loads:1290"]], "-L^3/3")
  expect_identical(got[["beam:test_parabolic_loads:1296"]], "40")
  expect_identical(got[["beam:test_parabolic_loads:1297"]], "0")
  expect_identical(got[["beam:test_parabolic_loads:1302"]], "-L*P*SingularityFunction(x, L, 0) + P*SingularityFunction(x, 0, 1) - P*SingularityFunction(x, L, 1)")
  expect_identical(got[["beam:test_parabolic_loads:1310"]], "15625000")
  expect_identical(got[["beam:test_parabolic_loads:1311"]], "0")
})

test_that("beam test_cross_section matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
x <- Symbol("x")
})
.py_try({
y <- Symbol("y")
})
.py_try({
.tmpu_ <- .py_symbols("R1, R2"); R1 <- .py_getitem(.tmpu_, 0L); R2 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
I <- Symbol("I")
})
.py_try({
l <- Symbol("l")
})
.py_try({
E <- Symbol("E")
})
.py_try({
.tmpu_ <- .py_symbols("C3, C4"); C3 <- .py_getitem(.tmpu_, 0L); C4 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
.tmpu_ <- .py_symbols("a, c, g, h, r, n"); a <- .py_getitem(.tmpu_, 0L); c <- .py_getitem(.tmpu_, 1L); g <- .py_getitem(.tmpu_, 2L); h <- .py_getitem(.tmpu_, 3L); r <- .py_getitem(.tmpu_, 4L); n <- .py_getitem(.tmpu_, 5L)
})
.py_try({
b0 <- Beam(l, E, I)
})
.py_try({
.py_rec("beam:test_cross_section:1323", function() .py_attr(b0, "second_moment"))
})
.py_try({
.py_rec("beam:test_cross_section:1324", function() .py_attr(b0, "cross_section"))
})
.py_try({
.py_setattr(b0, "cross_section", Circle(.py_tup(0L, 0L), 5L))
})
.py_try({
.py_rec("beam:test_cross_section:1326", function() .py_attr(b0, "second_moment"))
})
.py_try({
.py_rec("beam:test_cross_section:1327", function() .py_attr(b0, "cross_section"))
})
.py_try({
.py_setattr(b0, "second_moment", .py_sub(.py_mul(2L, n), 6L))
})
.py_try({
.py_rec("beam:test_cross_section:1329", function() .py_attr(b0, "second_moment"))
})
.py_try({
.py_rec("beam:test_cross_section:1330", function() .py_attr(b0, "cross_section"))
})
.t1331 <- function() {
  .py_setattr(b0, "second_moment", Circle(.py_tup(0L, 0L), 5L))
}
.py_try({
.py_rec("beam:test_cross_section:1331", .t1331)
})
.py_try({
b1 <- Beam(50L, E, Circle(.py_tup(0L, 0L), r))
})
.py_try({
.py_rec("beam:test_cross_section:1336", function() .py_attr(b1, "cross_section"))
})
.py_try({
.py_rec("beam:test_cross_section:1337", function() .py_attr(b1, "second_moment"))
})
.py_try({
.py_method(b1, "apply_load", -10L, 0L, -1L)
})
.py_try({
.py_method(b1, "apply_load", R1, 5L, -1L)
})
.py_try({
.py_method(b1, "apply_load", R2, 50L, -1L)
})
.py_try({
.py_method(b1, "apply_load", 90L, 45L, -2L)
})
.py_try({
.py_method(b1, "solve_for_reaction_loads", R1, R2)
})
.py_try({
.py_rec("beam:test_cross_section:1344", function() .py_attr(b1, "load"))
})
.py_try({
.py_rec("beam:test_cross_section:1346", function() .py_method(b1, "bending_moment"))
})
.py_try({
q <- .py_div(.py_add(.py_add(.py_add(.py_mul(-5L, SingularityFunction(x, 0L, 2L)), .py_div(.py_mul(41L, SingularityFunction(x, 5L, 2L)), .py_S(9L))), .py_mul(90L, SingularityFunction(x, 45L, 1L))), .py_div(.py_mul(4L, SingularityFunction(x, 50L, 2L)), .py_S(9L))), .py_mul(.py_mul(.py_mul(S_Pi, E), r), .py_pow(Abs(r), 3L)))
})
.py_try({
.py_rec("beam:test_cross_section:1350", function() .py_method(b1, "slope"))
})
.py_try({
q <- .py_div(.py_add(.py_add(.py_add(.py_div(.py_mul(-5L, SingularityFunction(x, 0L, 3L)), 3L), .py_div(.py_mul(41L, SingularityFunction(x, 5L, 3L)), 27L)), .py_mul(45L, SingularityFunction(x, 45L, 2L))), .py_div(.py_mul(4L, SingularityFunction(x, 50L, 3L)), 27L)), .py_mul(.py_mul(.py_mul(S_Pi, E), r), .py_pow(Abs(r), 3L)))
})
.py_try({
.py_rec("beam:test_cross_section:1353", function() .py_method(b1, "deflection"))
})
.py_try({
b2 <- Beam(20L, E, Polygon(.py_tup(0L, 0L), .py_tup(a, 0L), .py_tup(a, c), .py_tup(0L, c)))
})
.py_try({
.py_rec("beam:test_cross_section:1357", function() .py_attr(b2, "cross_section"))
})
.py_try({
.py_rec("beam:test_cross_section:1358", function() .py_attr(b2, "second_moment"))
})
.py_try({
b3 <- Beam(15L, E, Triangle(.py_tup(0L, 0L), .py_tup(g, 0L), .py_tup(.py_div(g, 2L), h)))
})
.py_try({
.py_rec("beam:test_cross_section:1361", function() .py_attr(b3, "cross_section"))
})
.py_try({
.py_rec("beam:test_cross_section:1362", function() .py_attr(b3, "second_moment"))
})
.py_try({
b <- .py_method(b2, "join", b3, "fixed")
})
.py_try({
.py_method(b, "apply_load", -30L, 0L, -1L)
})
.py_try({
.py_method(b, "apply_load", 65L, 0L, -2L)
})
.py_try({
.py_method(b, "apply_load", 40L, 0L, -1L)
})
.py_try({
.py_setattr(b, "bc_slope", list(.py_tup(0L, 0L)))
})
.py_try({
.py_setattr(b, "bc_deflection", list(.py_tup(0L, 0L)))
})
.py_try({
.py_rec("beam:test_cross_section:1372", function() .py_attr(b, "second_moment"))
})
.py_try({
.py_rec("beam:test_cross_section:1373", function() .py_attr(b, "cross_section"))
})
.py_try({
.py_rec("beam:test_cross_section:1374", function() .py_attr(b, "length"))
})
.py_try({
.py_rec("beam:test_cross_section:1375", function() .py_method(.py_method(b, "slope"), "subs", x, 7L))
})
.py_try({
.py_rec("beam:test_cross_section:1376", function() .py_method(.py_method(b, "slope"), "subs", x, 25L))
})
.py_try({
.py_rec("beam:test_cross_section:1377", function() .py_method(.py_method(b, "deflection"), "subs", x, 30L))
}))---")
  expect_identical(got[["beam:test_cross_section:1323"]], "I")
  expect_identical(got[["beam:test_cross_section:1324"]], "None")
  expect_identical(got[["beam:test_cross_section:1326"]], "625*pi/4")
  expect_identical(got[["beam:test_cross_section:1327"]], "Circle(Point2D(0, 0), 5)")
  expect_identical(got[["beam:test_cross_section:1329"]], "2*n - 6")
  expect_identical(got[["beam:test_cross_section:1330"]], "None")
  expect_identical(got[["beam:test_cross_section:1331"]], "ERR ValueError: To update cross-section geometry use `cross_section` attribute")
  expect_identical(got[["beam:test_cross_section:1336"]], "Circle(Point2D(0, 0), r)")
  expect_identical(got[["beam:test_cross_section:1337"]], "pi*r*Abs(r)^3/4")
  expect_identical(got[["beam:test_cross_section:1344"]], "-10*SingularityFunction(x, 0, -1) + 82*SingularityFunction(x, 5, -1)/9 + 90*SingularityFunction(x, 45, -2) + 8*SingularityFunction(x, 50, -1)/9")
  expect_identical(got[["beam:test_cross_section:1346"]], "10*SingularityFunction(x, 0, 1) - 82*SingularityFunction(x, 5, 1)/9 - 90*SingularityFunction(x, 45, 0) - 8*SingularityFunction(x, 50, 1)/9")
  expect_identical(got[["beam:test_cross_section:1350"]], "C3 + 4*(-5*SingularityFunction(x, 0, 2) + 41*SingularityFunction(x, 5, 2)/9 + 90*SingularityFunction(x, 45, 1) + 4*SingularityFunction(x, 50, 2)/9)/(pi*E*r*Abs(r)^3)")
  expect_identical(got[["beam:test_cross_section:1353"]], "C3*x + C4 + 4*(-5*SingularityFunction(x, 0, 3)/3 + 41*SingularityFunction(x, 5, 3)/27 + 45*SingularityFunction(x, 45, 2) + 4*SingularityFunction(x, 50, 3)/27)/(pi*E*r*Abs(r)^3)")
  expect_identical(got[["beam:test_cross_section:1357"]], "Polygon(Point2D(0, 0), Point2D(a, 0), Point2D(a, c), Point2D(0, c))")
  expect_identical(got[["beam:test_cross_section:1358"]], "a*c^3/12")
  expect_identical(got[["beam:test_cross_section:1361"]], "Triangle(Point2D(0, 0), Point2D(g, 0), Point2D(g/2, h))")
  expect_identical(got[["beam:test_cross_section:1362"]], "g*h^3/36")
  expect_identical(got[["beam:test_cross_section:1372"]], "Piecewise((a*c^3/12, x <= 20), (g*h^3/36, x <= 35))")
  expect_identical(got[["beam:test_cross_section:1373"]], "None")
  expect_identical(got[["beam:test_cross_section:1374"]], "35")
  expect_identical(got[["beam:test_cross_section:1375"]], "8400/(E*a*c^3)")
  expect_identical(got[["beam:test_cross_section:1376"]], "52200/(E*g*h^3) + 39600/(E*a*c^3)")
  expect_identical(got[["beam:test_cross_section:1377"]], "-537000/(E*g*h^3) - 712000/(E*a*c^3)")
})

test_that("beam test_max_shear_force_Beam3D matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
x <- Symbol("x")
})
.py_try({
y <- Symbol("y")
})
.py_try({
.tmpu_ <- .py_symbols("R1, R2"); R1 <- .py_getitem(.tmpu_, 0L); R2 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
x <- .py_symbols("x")
})
.py_try({
b <- Beam3D(20L, 40L, 21L, 100L, 25L)
})
.py_try({
.py_method(b, "apply_load", 15L, start = 0L, order = 0L, dir = "z")
})
.py_try({
.py_method(b, "apply_load", .py_mul(12L, x), start = 0L, order = 0L, dir = "y")
})
.py_try({
.py_setattr(b, "bc_deflection", list(.py_tup(0L, list(0L, 0L, 0L)), .py_tup(20L, list(0L, 0L, 0L))))
})
.py_try({
.py_rec("beam:test_max_shear_force_Beam3D:1385", function() .py_method(b, "max_shear_force"))
}))---")
  expect_identical(got[["beam:test_max_shear_force_Beam3D:1385"]], "[(0, 0), (20, 2400), (20, 300)]")
})

test_that("beam test_max_bending_moment_Beam3D matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
x <- Symbol("x")
})
.py_try({
y <- Symbol("y")
})
.py_try({
.tmpu_ <- .py_symbols("R1, R2"); R1 <- .py_getitem(.tmpu_, 0L); R2 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
x <- .py_symbols("x")
})
.py_try({
b <- Beam3D(20L, 40L, 21L, 100L, 25L)
})
.py_try({
.py_method(b, "apply_load", 15L, start = 0L, order = 0L, dir = "z")
})
.py_try({
.py_method(b, "apply_load", .py_mul(12L, x), start = 0L, order = 0L, dir = "y")
})
.py_try({
.py_setattr(b, "bc_deflection", list(.py_tup(0L, list(0L, 0L, 0L)), .py_tup(20L, list(0L, 0L, 0L))))
})
.py_try({
.py_rec("beam:test_max_bending_moment_Beam3D:1393", function() .py_method(b, "max_bmoment"))
}))---")
  expect_identical(got[["beam:test_max_bending_moment_Beam3D:1393"]], "[(0, 0), (20, 3000), (20, 16000)]")
})

test_that("beam test_torsion_Beam3D matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
x <- Symbol("x")
})
.py_try({
y <- Symbol("y")
})
.py_try({
.tmpu_ <- .py_symbols("R1, R2"); R1 <- .py_getitem(.tmpu_, 0L); R2 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
x <- .py_symbols("x")
})
.py_try({
b <- Beam3D(20L, 40L, 21L, 100L, 25L)
})
.py_try({
.py_method(b, "apply_moment_load", 15L, 5L, -2L, dir = "x")
})
.py_try({
.py_method(b, "apply_moment_load", 25L, 10L, -2L, dir = "x")
})
.py_try({
.py_method(b, "apply_moment_load", -5L, 20L, -2L, dir = "x")
})
.py_try({
.py_method(b, "solve_for_torsion")
})
.py_try({
.py_rec("beam:test_torsion_Beam3D:1414", function() .py_method(.py_method(b, "angular_deflection"), "subs", x, 3L))
})
.py_try({
.py_rec("beam:test_torsion_Beam3D:1415", function() .py_method(.py_method(b, "angular_deflection"), "subs", x, 9L))
})
.py_try({
.py_rec("beam:test_torsion_Beam3D:1416", function() .py_method(.py_method(b, "angular_deflection"), "subs", x, 12L))
})
.py_try({
.py_rec("beam:test_torsion_Beam3D:1417", function() .py_method(.py_method(b, "angular_deflection"), "subs", x, 17L))
})
.py_try({
.py_rec("beam:test_torsion_Beam3D:1418", function() .py_method(.py_method(b, "angular_deflection"), "subs", x, 20L))
}))---")
  expect_identical(got[["beam:test_torsion_Beam3D:1414"]], "1/40")
  expect_identical(got[["beam:test_torsion_Beam3D:1415"]], "17/280")
  expect_identical(got[["beam:test_torsion_Beam3D:1416"]], "53/840")
  expect_identical(got[["beam:test_torsion_Beam3D:1417"]], "2/35")
  expect_identical(got[["beam:test_torsion_Beam3D:1418"]], "3/56")
})

test_that("column test_column matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
x <- Symbol("x")
})
.py_try({
x_1 <- Symbol("x_1")
})
.py_try({
E <- Symbol("E")
})
.py_try({
E_1 <- Symbol("E_1")
})
.py_try({
A <- Symbol("A")
})
.py_try({
A_1 <- Symbol("A_1")
})
.py_try({
c <- Column(10L, E, A)
})
.py_try({
.py_rec("column:test_column:19", function() .py_attr(c, "length"))
})
.py_try({
.py_rec("column:test_column:20", function() .py_attr(c, "elastic_modulus"))
})
.py_try({
.py_rec("column:test_column:21", function() .py_attr(c, "area"))
})
.py_try({
.py_rec("column:test_column:22", function() .py_attr(c, "variable"))
})
.py_try({
.py_rec("column:test_column:25", function() .py_str(c))
})
.py_try({
.py_setattr(c, "length", 20L)
})
.py_try({
.py_rec("column:test_column:29", function() .py_attr(c, "length"))
})
.py_try({
.py_setattr(c, "elastic_modulus", E_1)
})
.py_try({
.py_rec("column:test_column:33", function() .py_attr(c, "elastic_modulus"))
})
.py_try({
.py_setattr(c, "area", A_1)
})
.py_try({
.py_rec("column:test_column:37", function() .py_attr(c, "area"))
})
.py_try({
.py_setattr(c, "variable", x_1)
})
.py_try({
.py_rec("column:test_column:41", function() .py_attr(c, "variable"))
})
.py_try({
c2 <- Column(10L, E, A)
})
.py_try({
.py_method(c2, "apply_support", 0L)
})
.py_try({
.py_rec("column:test_column:46", function() .py_attr(c2, "_bc_extension"))
})
.py_try({
.py_method(c2, "apply_support", 10L)
})
.py_try({
.py_rec("column:test_column:48", function() .py_attr(c2, "_bc_extension"))
})
.py_try({
.py_method(c2, "apply_support", 5L)
})
.py_try({
.py_rec("column:test_column:50", function() .py_attr(c2, "_bc_extension"))
})
.py_try({
c3 <- Column(10L, E, A)
})
.py_try({
.py_method(c3, "apply_load", -10L, 10L, -1L)
})
.py_try({
.py_rec("column:test_column:55", function() .py_attr(c3, "applied_loads"))
})
.py_try({
.py_method(c3, "apply_load", 10L, 5L, -1L)
})
.py_try({
.py_rec("column:test_column:57", function() .py_attr(c3, "applied_loads"))
})
.py_try({
p <- .py_attr(c3, "load")
})
.py_try({
q <- .py_add(.py_mul(-10L, SingularityFunction(x, 10L, -1L)), .py_mul(10L, SingularityFunction(x, 5L, -1L)))
})
.py_try({
.py_rec("column:test_column:62", function() p)
})
.py_try({
c4 <- Column(10L, E, A)
})
.py_try({
.tmpu_ <- .py_symbols("F G"); F <- .py_getitem(.tmpu_, 0L); G <- .py_getitem(.tmpu_, 1L)
})
.py_try({
.py_method(c4, "apply_load", F, 0L, -1L)
})
.py_try({
.py_method(c4, "apply_load", .py_neg(G), 10L, -1L)
})
.py_try({
.py_rec("column:test_column:69", function() .py_attr(c4, "applied_loads"))
})
.py_try({
p <- .py_attr(c4, "load")
})
.py_try({
q <- .py_sub(.py_mul(F, SingularityFunction(x, 0L, -1L)), .py_mul(G, SingularityFunction(x, 10L, -1L)))
})
.py_try({
.py_rec("column:test_column:74", function() p)
})
.py_try({
c5 <- Column(10L, E, A)
})
.py_try({
.py_method(c5, "apply_support", 0L)
})
.py_try({
.py_method(c5, "apply_load", 10L, 0L, -1L)
})
.py_try({
.py_method(c5, "apply_load", 10L, 5L, -1L)
})
.py_try({
p <- .py_attr(c5, "load")
})
.py_try({
R_0 <- Symbol("R_0")
})
.py_try({
q <- .py_add(.py_add(.py_mul(R_0, SingularityFunction(x, 0L, -1L)), .py_mul(10L, SingularityFunction(x, 0L, -1L))), .py_mul(10L, SingularityFunction(x, 5L, -1L)))
})
.py_try({
.py_rec("column:test_column:85", function() p)
}))---")
  expect_identical(got[["column:test_column:19"]], "10")
  expect_identical(got[["column:test_column:20"]], "E")
  expect_identical(got[["column:test_column:21"]], "A")
  expect_identical(got[["column:test_column:22"]], "x")
  expect_identical(got[["column:test_column:25"]], "'Column(10, E, A)'")
  expect_identical(got[["column:test_column:29"]], "20")
  expect_identical(got[["column:test_column:33"]], "E_1")
  expect_identical(got[["column:test_column:37"]], "A_1")
  expect_identical(got[["column:test_column:41"]], "x_1")
  expect_identical(got[["column:test_column:46"]], "[0]")
  expect_identical(got[["column:test_column:48"]], "[0, 10]")
  expect_identical(got[["column:test_column:50"]], "[0, 10, 5]")
  expect_identical(got[["column:test_column:55"]], "[(-10, 10, -1, None)]")
  expect_identical(got[["column:test_column:57"]], "[(-10, 10, -1, None), (10, 5, -1, None)]")
  expect_identical(got[["column:test_column:62"]], "10*SingularityFunction(x, 5, -1) - 10*SingularityFunction(x, 10, -1)")
  expect_identical(got[["column:test_column:69"]], "[(F, 0, -1, None), (-G, 10, -1, None)]")
  expect_identical(got[["column:test_column:74"]], "F*SingularityFunction(x, 0, -1) - G*SingularityFunction(x, 10, -1)")
  expect_identical(got[["column:test_column:85"]], "R_0*SingularityFunction(x, 0, -1) + 10*SingularityFunction(x, 0, -1) + 10*SingularityFunction(x, 5, -1)")
})

test_that("column test_supports_using_manual_method matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
x <- Symbol("x")
})
.py_try({
x_1 <- Symbol("x_1")
})
.py_try({
E <- Symbol("E")
})
.py_try({
A <- Symbol("A")
})
.py_try({
R_0 <- Symbol("R_0")
})
.py_try({
c5 <- Column(10L, E, A)
})
.py_try({
.py_method(c5, "apply_load", R_0, 0L, -1L)
})
.py_try({
.py_setattr(c5, "_bc_extension", .py_list_mut(.py_attr(c5, "_bc_extension"), "append", 0L))
})
.py_try({
.py_method(c5, "apply_load", 10L, 0L, -1L)
})
.py_try({
.py_method(c5, "apply_load", 10L, 5L, -1L)
})
.py_try({
p <- .py_attr(c5, "load")
})
.py_try({
q <- .py_add(.py_add(.py_mul(R_0, SingularityFunction(x, 0L, -1L)), .py_mul(10L, SingularityFunction(x, 0L, -1L))), .py_mul(10L, SingularityFunction(x, 5L, -1L)))
})
.py_try({
.py_rec("column:test_supports_using_manual_method:103", function() p)
})
.py_try({
.tmpu_ <- .py_symbols("E A"); E <- .py_getitem(.tmpu_, 0L); A <- .py_getitem(.tmpu_, 1L)
})
.py_try({
R_0 <- Symbol("R_0")
})
.py_try({
c <- Column(5L, E, A)
})
.py_try({
.py_method(c, "apply_load", R_0, 0L, -1L)
})
.py_try({
.py_setattr(c, "_bc_extension", .py_list_mut(.py_attr(c, "_bc_extension"), "append", 0L))
})
.py_try({
.py_method(c, "apply_load", -1L, .py_flt(2.5), -1L)
})
.py_try({
.py_method(c, "apply_load", 2L, 5L, -1L)
})
.py_try({
.py_method(c, "solve_for_reaction_loads", R_0)
})
.py_try({
p <- .py_attr(c, "reaction_loads")
})
.py_try({
q <- .py_dict(list(R_0), list(-1L))
})
.py_try({
.py_rec("column:test_supports_using_manual_method:118", function() p)
})
.py_try({
c1 <- Column(5L, E, A)
})
.py_try({
R_5 <- Symbol("R_5")
})
.py_try({
.py_method(c1, "apply_load", R_5, 5L, -1L)
})
.py_try({
.py_setattr(c1, "_bc_extension", .py_list_mut(.py_attr(c1, "_bc_extension"), "append", 5L))
})
.py_try({
.py_method(c1, "apply_load", -4L, .py_flt(2.5), -1L)
})
.py_try({
.py_method(c1, "apply_load", 2L, 5L, -1L)
})
.py_try({
.py_method(c1, "solve_for_reaction_loads", R_5)
})
.py_try({
p <- .py_attr(c1, "reaction_loads")
})
.py_try({
q <- .py_dict(list(R_5), list(2L))
})
.py_try({
.py_rec("column:test_supports_using_manual_method:131", function() p)
})
.py_try({
c2 <- Column(10L, E, A)
})
.py_try({
R_0 <- Symbol("R_0")
})
.py_try({
R_10 <- Symbol("R_10")
})
.py_try({
.py_method(c2, "apply_load", R_10, 10L, -1L)
})
.py_try({
.py_method(c2, "apply_load", R_0, 0L, -1L)
})
.py_try({
.py_setattr(c2, "_bc_extension", .py_list_mut(.py_attr(c2, "_bc_extension"), "append", 0L))
})
.py_try({
.py_setattr(c2, "_bc_extension", .py_list_mut(.py_attr(c2, "_bc_extension"), "append", 10L))
})
.py_try({
.py_method(c2, "apply_load", -1L, 5L, -1L)
})
.py_try({
.py_method(c2, "solve_for_reaction_loads", R_0, R_10)
})
.py_try({
p <- .py_attr(c2, "reaction_loads")
})
.py_try({
q <- .py_dict(list(R_0, R_10), list(Rational(1L, 2L), Rational(1L, 2L)))
})
.py_try({
.py_rec("column:test_supports_using_manual_method:146", function() p)
}))---")
  expect_identical(got[["column:test_supports_using_manual_method:103"]], "R_0*SingularityFunction(x, 0, -1) + 10*SingularityFunction(x, 0, -1) + 10*SingularityFunction(x, 5, -1)")
  expect_identical(got[["column:test_supports_using_manual_method:118"]], "{R_0: -1}")
  expect_identical(got[["column:test_supports_using_manual_method:131"]], "{R_5: 2}")
  expect_identical(got[["column:test_supports_using_manual_method:146"]], "{R_0: 1/2, R_10: 1/2}")
})

test_that("column test_distributed_loads matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
x <- Symbol("x")
})
.py_try({
x_1 <- Symbol("x_1")
})
.py_try({
.tmpu_ <- .py_symbols("E A"); E <- .py_getitem(.tmpu_, 0L); A <- .py_getitem(.tmpu_, 1L)
})
.py_try({
c <- Column(10L, E, A)
})
.py_try({
.py_method(c, "apply_load", 5L, 0L, -1L)
})
.py_try({
.py_method(c, "apply_load", 5L, 2L, 0L, end = 4L)
})
.py_try({
.py_method(c, "apply_load", 5L, 4L, 1L, end = 6L)
})
.py_try({
p <- .py_attr(c, "applied_loads")
})
.py_try({
q <- list(.py_tup(5L, 0L, -1L, NULL), .py_tup(5L, 2L, 0L, 4L), .py_tup(5L, 4L, 1L, 6L))
})
.py_try({
.py_rec("column:test_distributed_loads:160", function() p)
})
.py_try({
p <- .py_attr(c, "load")
})
.py_try({
q <- .py_sub(.py_sub(.py_add(.py_sub(.py_add(.py_mul(5L, SingularityFunction(x, 0L, -1L)), .py_mul(5L, SingularityFunction(x, 2L, 0L))), .py_mul(5L, SingularityFunction(x, 4L, 0L))), .py_mul(5L, SingularityFunction(x, 4L, 1L))), .py_mul(10L, SingularityFunction(x, 6L, 0L))), .py_mul(5L, SingularityFunction(x, 6L, 1L)))
})
.py_try({
.py_rec("column:test_distributed_loads:171", function() p)
})
.py_try({
.tmpu_ <- .py_symbols("E A F L"); E <- .py_getitem(.tmpu_, 0L); A <- .py_getitem(.tmpu_, 1L); F <- .py_getitem(.tmpu_, 2L); L <- .py_getitem(.tmpu_, 3L)
})
.py_try({
c1 <- Column(L, E, A)
})
.py_try({
.py_method(c1, "apply_load", F, 0L, -1L)
})
.py_try({
.py_method(c1, "apply_load", F, .py_div(L, 4L), 0L, end = .py_div(L, 2L))
})
.py_try({
.py_method(c1, "apply_load", F, .py_div(L, 2L), 1L, end = L)
})
.py_try({
p <- .py_attr(c1, "applied_loads")
})
.py_try({
q <- list(.py_tup(F, 0L, -1L, NULL), .py_tup(F, .py_div(L, 4L), 0L, .py_div(L, 2L)), .py_tup(F, .py_div(L, 2L), 1L, L))
})
.py_try({
.py_rec("column:test_distributed_loads:184", function() p)
})
.py_try({
p <- .py_attr(c1, "load")
})
.py_try({
q <- .py_sub(.py_sub(.py_add(.py_sub(.py_add(.py_mul(F, SingularityFunction(x, 0L, -1L)), .py_mul(F, SingularityFunction(x, .py_div(L, 4L), 0L))), .py_mul(F, SingularityFunction(x, .py_div(L, 2L), 0L))), .py_mul(F, SingularityFunction(x, .py_div(L, 2L), 1L))), .py_mul(.py_mul(.py_div(L, 2L), F), SingularityFunction(x, L, 0L))), .py_mul(F, SingularityFunction(x, L, 1L)))
})
.py_try({
.py_rec("column:test_distributed_loads:196", function() p)
}))---")
  expect_identical(got[["column:test_distributed_loads:160"]], "[(5, 0, -1, None), (5, 2, 0, 4), (5, 4, 1, 6)]")
  expect_identical(got[["column:test_distributed_loads:171"]], "5*SingularityFunction(x, 0, -1) + 5*SingularityFunction(x, 2, 0) - 5*SingularityFunction(x, 4, 0) + 5*SingularityFunction(x, 4, 1) - 10*SingularityFunction(x, 6, 0) - 5*SingularityFunction(x, 6, 1)")
  expect_identical(got[["column:test_distributed_loads:184"]], "[(F, 0, -1, None), (F, L/4, 0, L/2), (F, L/2, 1, L)]")
  expect_identical(got[["column:test_distributed_loads:196"]], "-F*L*SingularityFunction(x, L, 0)/2 + F*SingularityFunction(x, 0, -1) + F*SingularityFunction(x, L/4, 0) - F*SingularityFunction(x, L/2, 0) + F*SingularityFunction(x, L/2, 1) - F*SingularityFunction(x, L, 1)")
})

test_that("column test_remove_load matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
x <- Symbol("x")
})
.py_try({
x_1 <- Symbol("x_1")
})
.py_try({
.tmpu_ <- .py_symbols("E A"); E <- .py_getitem(.tmpu_, 0L); A <- .py_getitem(.tmpu_, 1L)
})
.py_try({
c <- Column(10L, E, A)
})
.py_try({
.py_method(c, "apply_load", 10L, 0L, -1L)
})
.py_try({
.py_method(c, "apply_load", -10L, 10L, -1L)
})
.py_try({
.py_method(c, "remove_load", 10L, 0L, -1L, NULL)
})
.py_try({
p <- .py_attr(c, "applied_loads")
})
.py_try({
q <- list(.py_tup(-10L, 10L, -1L, NULL))
})
.py_try({
.py_rec("column:test_remove_load:211", function() p)
})
.py_try({
.py_method(c, "apply_load", 5L, 0L, 0L, end = 10L)
})
.py_try({
.py_method(c, "remove_load", 5L, 0L, 0L, end = 10L)
})
.py_try({
p <- .py_attr(c, "load")
})
.py_try({
q <- .py_mul(-10L, SingularityFunction(x, 10L, -1L))
})
.py_try({
.py_rec("column:test_remove_load:218", function() p)
})
.py_try({
.py_method(c, "apply_load", 5L, 0L, 1L, end = 5L)
})
.py_try({
.py_method(c, "remove_load", 5L, 0L, 1L, 5L)
})
.py_try({
p <- .py_attr(c, "load")
})
.py_try({
q <- .py_mul(-10L, SingularityFunction(x, 10L, -1L))
})
.py_try({
.py_rec("column:test_remove_load:225", function() p)
})
.py_try({
.tmpu_ <- .py_symbols("E A F L"); E <- .py_getitem(.tmpu_, 0L); A <- .py_getitem(.tmpu_, 1L); F <- .py_getitem(.tmpu_, 2L); L <- .py_getitem(.tmpu_, 3L)
})
.py_try({
c <- Column(L, E, A)
})
.py_try({
.py_method(c, "apply_load", F, 0L, -1L)
})
.py_try({
.py_method(c, "apply_load", .py_neg(F), L, -1L)
})
.py_try({
.py_method(c, "remove_load", F, 0L, -1L, NULL)
})
.py_try({
p <- .py_attr(c, "applied_loads")
})
.py_try({
q <- list(.py_tup(.py_neg(F), L, -1L, NULL))
})
.py_try({
.py_rec("column:test_remove_load:237", function() p)
})
.py_try({
.py_method(c, "apply_load", F, 0L, 0L, end = L)
})
.py_try({
.py_method(c, "remove_load", F, 0L, 0L, end = L)
})
.py_try({
p <- .py_attr(c, "load")
})
.py_try({
q <- .py_mul(.py_neg(F), SingularityFunction(x, L, -1L))
})
.py_try({
.py_rec("column:test_remove_load:244", function() p)
})
.py_try({
.py_method(c, "apply_load", F, 0L, 1L, end = .py_div(L, 2L))
})
.py_try({
.py_method(c, "remove_load", F, 0L, 1L, .py_div(L, 2L))
})
.py_try({
p <- .py_attr(c, "load")
})
.py_try({
q <- .py_mul(.py_neg(F), SingularityFunction(x, L, -1L))
})
.py_try({
.py_rec("column:test_remove_load:251", function() p)
})
.py_try({
c1 <- Column(10L, E, A)
})
.py_try({
.py_method(c1, "apply_load", 1L, 5L, -1L)
})
.py_try({
.py_method(c1, "apply_load", 10L, 10L, 1L, end = 0L)
})
.py_try({
.py_method(c1, "remove_load", 10L, 10L, 1L, end = 0L)
})
.py_try({
p <- .py_attr(c1, "applied_loads")
})
.py_try({
q <- list(.py_tup(1L, 5L, -1L, NULL))
})
.py_try({
.py_rec("column:test_remove_load:261", function() p)
})
.py_try({
p <- .py_attr(c1, "load")
})
.py_try({
q <- SingularityFunction(x, 5L, -1L)
})
.py_try({
.py_rec("column:test_remove_load:265", function() p)
}))---")
  expect_identical(got[["column:test_remove_load:211"]], "[(-10, 10, -1, None)]")
  expect_identical(got[["column:test_remove_load:218"]], "-10*SingularityFunction(x, 10, -1)")
  expect_identical(got[["column:test_remove_load:225"]], "-10*SingularityFunction(x, 10, -1)")
  expect_identical(got[["column:test_remove_load:237"]], "[(-F, L, -1, None)]")
  expect_identical(got[["column:test_remove_load:244"]], "-F*SingularityFunction(x, L, -1)")
  expect_identical(got[["column:test_remove_load:251"]], "-F*SingularityFunction(x, L, -1)")
  expect_identical(got[["column:test_remove_load:261"]], "[(1, 5, -1, None)]")
  expect_identical(got[["column:test_remove_load:265"]], "SingularityFunction(x, 5, -1)")
})

test_that("column test_reactions_point_loads matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
x <- Symbol("x")
})
.py_try({
x_1 <- Symbol("x_1")
})
.py_try({
.tmpu_ <- .py_symbols("E A"); E <- .py_getitem(.tmpu_, 0L); A <- .py_getitem(.tmpu_, 1L)
})
.py_try({
c <- Column(5L, E, A)
})
.py_try({
.py_method(c, "apply_support", 0L)
})
.py_try({
.py_method(c, "apply_load", -1L, .py_flt(2.5), -1L)
})
.py_try({
.py_method(c, "apply_load", 2L, 5L, -1L)
})
.py_try({
.py_method(c, "solve_for_reaction_loads")
})
.py_try({
p <- .py_attr(c, "reaction_loads")
})
.py_try({
R_0 <- Symbol("R_0")
})
.py_try({
q <- .py_dict(list(R_0), list(-1L))
})
.py_try({
.py_rec("column:test_reactions_point_loads:283", function() p)
})
.py_try({
c1 <- Column(5L, E, A)
})
.py_try({
.py_method(c1, "apply_support", 5L)
})
.py_try({
.py_method(c1, "apply_load", -4L, .py_flt(2.5), -1L)
})
.py_try({
.py_method(c1, "apply_load", 2L, 5L, -1L)
})
.py_try({
.py_method(c1, "solve_for_reaction_loads")
})
.py_try({
p <- .py_attr(c1, "reaction_loads")
})
.py_try({
R_5 <- Symbol("R_5")
})
.py_try({
q <- .py_dict(list(R_5), list(2L))
})
.py_try({
.py_rec("column:test_reactions_point_loads:295", function() p)
})
.py_try({
c2 <- Column(10L, E, A)
})
.py_try({
.py_method(c2, "apply_support", 10L)
})
.py_try({
.py_method(c2, "apply_support", 0L)
})
.py_try({
.py_method(c2, "apply_load", -1L, 5L, -1L)
})
.py_try({
.py_method(c2, "solve_for_reaction_loads")
})
.py_try({
p <- .py_attr(c2, "reaction_loads")
})
.py_try({
R_10 <- Symbol("R_10")
})
.py_try({
q <- .py_dict(list(R_0, R_10), list(Rational(1L, 2L), Rational(1L, 2L)))
})
.py_try({
.py_rec("column:test_reactions_point_loads:307", function() p)
})
.py_try({
c3 <- Column(10L, E, A)
})
.py_try({
.py_method(c3, "apply_support", 2L)
})
.py_try({
.py_method(c3, "apply_support", 8L)
})
.py_try({
.py_method(c3, "apply_load", -1L, 5L, -1L)
})
.py_try({
.py_method(c3, "apply_load", -1L, 0L, -1L)
})
.py_try({
.py_method(c3, "apply_load", -1L, 10L, -1L)
})
.py_try({
.py_method(c3, "solve_for_reaction_loads")
})
.py_try({
p <- .py_attr(c3, "reaction_loads")
})
.py_try({
.tmpu_ <- .py_symbols("R_2 R_8"); R_2 <- .py_getitem(.tmpu_, 0L); R_8 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
q <- .py_dict(list(R_2, R_8), list(Rational(3L, 2L), Rational(3L, 2L)))
})
.py_try({
.py_rec("column:test_reactions_point_loads:321", function() p)
})
.py_try({
c4 <- Column(10L, E, A)
})
.py_try({
.py_method(c4, "apply_support", 0L)
})
.py_try({
.py_method(c4, "apply_support", 10L)
})
.py_try({
.py_method(c4, "apply_load", -1L, 6L, -1L)
})
.py_try({
.py_method(c4, "solve_for_reaction_loads")
})
.py_try({
p <- .py_attr(c4, "reaction_loads")
})
.py_try({
q <- .py_dict(list(R_0, R_10), list(Rational(2L, 5L), Rational(3L, 5L)))
})
.py_try({
.py_rec("column:test_reactions_point_loads:332", function() p)
})
.py_try({
c5 <- Column(10L, E, A)
})
.py_try({
.py_method(c5, "apply_support", 10L)
})
.py_try({
.py_method(c5, "apply_support", 2L)
})
.py_try({
.py_method(c5, "apply_load", -1L, 0L, -1L)
})
.py_try({
.py_method(c5, "apply_load", 2L, 5L, -1L)
})
.py_try({
.py_method(c5, "apply_load", -3L, 8L, -1L)
})
.py_try({
.py_method(c5, "solve_for_reaction_loads")
})
.py_try({
p <- .py_attr(c5, "reaction_loads")
})
.py_try({
q <- .py_dict(list(R_2, R_10), list(Rational(1L, 2L), Rational(3L, 2L)))
})
.py_try({
.py_rec("column:test_reactions_point_loads:345", function() p)
})
.py_try({
c6 <- Column(10L, E, A)
})
.py_try({
.py_method(c6, "apply_support", 0L)
})
.py_try({
.py_method(c6, "apply_support", 5L)
})
.py_try({
.py_method(c6, "apply_support", 10L)
})
.py_try({
.py_method(c6, "apply_load", -1L, 2L, -1L)
})
.py_try({
.py_method(c6, "apply_load", -1L, 8L, -1L)
})
.py_try({
.py_method(c6, "solve_for_reaction_loads")
})
.py_try({
p <- .py_attr(c6, "reaction_loads")
})
.py_try({
q <- .py_dict(list(R_0, R_5, R_10), list(Rational(3L, 5L), Rational(4L, 5L), Rational(3L, 5L)))
})
.py_try({
.py_rec("column:test_reactions_point_loads:358", function() p)
})
.py_try({
c6 <- Column(10L, E, A)
})
.py_try({
.py_method(c6, "apply_support", 0L)
})
.py_try({
.py_method(c6, "apply_support", 5L)
})
.py_try({
.py_method(c6, "apply_support", 10L)
})
.py_try({
.py_method(c6, "apply_load", -1L, 2L, -1L)
})
.py_try({
.py_method(c6, "solve_for_reaction_loads")
})
.py_try({
p <- .py_attr(c6, "reaction_loads")
})
.py_try({
q <- .py_dict(list(R_0, R_5, R_10), list(Rational(3L, 5L), Rational(2L, 5L), 0L))
})
.py_try({
.py_rec("column:test_reactions_point_loads:370", function() p)
})
.py_try({
.tmpu_ <- .py_symbols("L F", positive = TRUE); L <- .py_getitem(.tmpu_, 0L); F <- .py_getitem(.tmpu_, 1L)
})
.py_try({
c7 <- Column(L, E, A)
})
.py_try({
.py_method(c7, "apply_support", 0L)
})
.py_try({
.py_method(c7, "apply_support", L)
})
.py_try({
.py_method(c7, "apply_load", .py_neg(F), .py_div(L, 2L), -1L)
})
.py_try({
.py_method(c7, "solve_for_reaction_loads")
})
.py_try({
p <- .py_attr(c7, "reaction_loads")
})
.py_try({
R_L <- Symbol("R_L")
})
.py_try({
q <- .py_dict(list(R_0, R_L), list(.py_div(F, 2L), .py_div(F, 2L)))
})
.py_try({
.py_rec("column:test_reactions_point_loads:384", function() p)
})
.py_try({
c8 <- Column(10L, E, A)
})
.py_try({
.py_method(c8, "apply_support", 0L)
})
.py_try({
.py_method(c8, "apply_support", 10L)
})
.py_try({
.py_method(c8, "apply_load", -1L, 0L, -1L)
})
.py_try({
.py_method(c8, "apply_load", -1L, 5L, -1L)
})
.py_try({
.py_method(c8, "apply_load", -1L, 10L, -1L)
})
.py_try({
.py_method(c8, "solve_for_reaction_loads")
})
.py_try({
p <- .py_attr(c8, "reaction_loads")
})
.py_try({
q <- .py_dict(list(R_0, R_10), list(Rational(3L, 2L), Rational(3L, 2L)))
})
.py_try({
.py_rec("column:test_reactions_point_loads:397", function() p)
})
.py_try({
.py_method(c8, "remove_load", -1L, 10L, -1L)
})
.py_try({
.py_method(c8, "solve_for_reaction_loads")
})
.py_try({
p <- .py_attr(c8, "reaction_loads")
})
.py_try({
q <- .py_dict(list(R_0, R_10), list(Rational(3L, 2L), Rational(1L, 2L)))
})
.py_try({
.py_rec("column:test_reactions_point_loads:404", function() p)
}))---")
  expect_identical(got[["column:test_reactions_point_loads:283"]], "{R_0: -1}")
  expect_identical(got[["column:test_reactions_point_loads:295"]], "{R_5: 2}")
  expect_identical(got[["column:test_reactions_point_loads:307"]], "{R_10: 1/2, R_0: 1/2}")
  expect_identical(got[["column:test_reactions_point_loads:321"]], "{R_2: 3/2, R_8: 3/2}")
  expect_identical(got[["column:test_reactions_point_loads:332"]], "{R_0: 2/5, R_10: 3/5}")
  expect_identical(got[["column:test_reactions_point_loads:345"]], "{R_10: 3/2, R_2: 1/2}")
  expect_identical(got[["column:test_reactions_point_loads:358"]], "{R_0: 3/5, R_5: 4/5, R_10: 3/5}")
  expect_identical(got[["column:test_reactions_point_loads:370"]], "{R_0: 3/5, R_5: 2/5, R_10: 0}")
  expect_identical(got[["column:test_reactions_point_loads:384"]], "{R_0: F/2, R_L: F/2}")
  expect_identical(got[["column:test_reactions_point_loads:397"]], "{R_0: 3/2, R_10: 3/2}")
  expect_identical(got[["column:test_reactions_point_loads:404"]], "{R_0: 3/2, R_10: 1/2}")
})

test_that("column test_reactions_higher_orders matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
x <- Symbol("x")
})
.py_try({
x_1 <- Symbol("x_1")
})
.py_try({
.tmpu_ <- .py_symbols("E A"); E <- .py_getitem(.tmpu_, 0L); A <- .py_getitem(.tmpu_, 1L)
})
.py_try({
c <- Column(10L, E, A)
})
.py_try({
.py_method(c, "apply_support", 0L)
})
.py_try({
.py_method(c, "apply_load", -1L, 0L, 0L, end = 10L)
})
.py_try({
.py_method(c, "solve_for_reaction_loads")
})
.py_try({
p <- .py_attr(c, "reaction_loads")
})
.py_try({
R_0 <- Symbol("R_0")
})
.py_try({
q <- .py_dict(list(R_0), list(10L))
})
.py_try({
.py_rec("column:test_reactions_higher_orders:421", function() p)
})
.py_try({
c1 <- Column(10L, E, A)
})
.py_try({
.py_method(c1, "apply_support", 0L)
})
.py_try({
.py_method(c1, "apply_support", 10L)
})
.py_try({
.py_method(c1, "apply_load", 1L, 0L, 0L, end = 5L)
})
.py_try({
.py_method(c1, "apply_load", 2L, 5L, 0L, end = 10L)
})
.py_try({
.py_method(c1, "solve_for_reaction_loads")
})
.py_try({
p <- .py_attr(c1, "reaction_loads")
})
.py_try({
R_10 <- Symbol("R_10")
})
.py_try({
q <- .py_dict(list(R_0, R_10), list(.py_neg(Rational(25L, 4L)), .py_neg(Rational(35L, 4L))))
})
.py_try({
.py_rec("column:test_reactions_higher_orders:434", function() p)
})
.py_try({
c2 <- Column(10L, E, A)
})
.py_try({
.py_method(c2, "apply_support", 10L)
})
.py_try({
.py_method(c2, "apply_load", 1L, 0L, 1L, end = 5L)
})
.py_try({
.py_method(c2, "solve_for_reaction_loads")
})
.py_try({
p <- .py_attr(c2, "reaction_loads")
})
.py_try({
q <- .py_dict(list(R_10), list(.py_neg(Rational(25L, 2L))))
})
.py_try({
.py_rec("column:test_reactions_higher_orders:444", function() p)
})
.py_try({
c3 <- Column(10L, E, A)
})
.py_try({
.py_method(c3, "apply_support", 0L)
})
.py_try({
.py_method(c3, "apply_support", 10L)
})
.py_try({
.py_method(c3, "apply_load", 1L, 0L, 1L, end = 10L)
})
.py_try({
.py_method(c3, "solve_for_reaction_loads")
})
.py_try({
p <- .py_attr(c3, "reaction_loads")
})
.py_try({
q <- .py_dict(list(R_0, R_10), list(.py_neg(Rational(50L, 3L)), .py_neg(Rational(100L, 3L))))
})
.py_try({
.py_rec("column:test_reactions_higher_orders:455", function() p)
})
.py_try({
c4 <- Column(10L, E, A)
})
.py_try({
.py_method(c4, "apply_support", 10L)
})
.py_try({
.py_method(c4, "apply_load", 1L, 0L, 2L, end = 10L)
})
.py_try({
.py_method(c4, "solve_for_reaction_loads")
})
.py_try({
p <- .py_attr(c4, "reaction_loads")
})
.py_try({
q <- .py_dict(list(R_10), list(.py_neg(Rational(1000L, 3L))))
})
.py_try({
.py_rec("column:test_reactions_higher_orders:465", function() p)
})
.py_try({
c5 <- Column(10L, E, A)
})
.py_try({
.py_method(c5, "apply_support", 0L)
})
.py_try({
.py_method(c5, "apply_support", 10L)
})
.py_try({
.py_method(c5, "apply_load", 100L, 0L, -1L)
})
.py_try({
.py_method(c5, "apply_load", 20L, 0L, 0L, end = 4L)
})
.py_try({
.py_method(c5, "apply_load", 20L, 4L, 0L, end = 8L)
})
.py_try({
.py_method(c5, "apply_load", 5L, 4L, 1L, end = 8L)
})
.py_try({
.py_method(c5, "apply_load", 40L, 8L, 0L, end = 10L)
})
.py_try({
.py_method(c5, "apply_load", 2L, 8L, 2L, end = 10L)
})
.py_try({
.py_method(c5, "solve_for_reaction_loads")
})
.py_try({
p <- .py_attr(c5, "reaction_loads")
})
.py_try({
q <- .py_dict(list(R_0, R_10), list(.py_neg(Rational(1088L, 5L)), .py_neg(Rational(2516L, 15L))))
})
.py_try({
.py_rec("column:test_reactions_higher_orders:482", function() p)
})
.py_try({
c6 <- Column(8L, 20000L, .py_flt(0.75))
})
.py_try({
.py_method(c6, "apply_support", 0L)
})
.py_try({
.py_method(c6, "apply_support", 8L)
})
.py_try({
.py_method(c6, "apply_load", -100L, 8L, -1L)
})
.py_try({
.py_method(c6, "apply_load", -20L, 0L, 0L, end = 8L)
})
.py_try({
.py_method(c6, "apply_load", -10L, 4L, 1L, end = 0L)
})
.py_try({
.py_method(c6, "solve_for_reaction_loads")
})
.py_try({
p <- .py_attr(c6, "load")
})
.py_try({
R_8 <- Symbol("R_8")
})
.py_try({
q <- .py_add(.py_sub(.py_sub(.py_add(.py_sub(.py_add(.py_mul(R_0, SingularityFunction(x, 0L, -1L)), .py_mul(R_8, SingularityFunction(x, 8L, -1L))), .py_mul(60L, SingularityFunction(x, 0L, 0L))), .py_mul(10L, SingularityFunction(x, 0L, 1L))), .py_mul(10L, SingularityFunction(x, 4L, 1L))), .py_mul(100L, SingularityFunction(x, 8L, -1L))), .py_mul(20L, SingularityFunction(x, 8L, 0L)))
})
.py_try({
.py_rec("column:test_reactions_higher_orders:502", function() p)
})
.py_try({
p <- .py_attr(c6, "reaction_loads")
})
.py_try({
q <- .py_dict(list(R_0, R_10), list(Rational(440L, 3L), Rational(580L, 3L)))
}))---")
  expect_identical(got[["column:test_reactions_higher_orders:421"]], "{R_0: 10}")
  expect_identical(got[["column:test_reactions_higher_orders:434"]], "{R_0: -25/4, R_10: -35/4}")
  expect_identical(got[["column:test_reactions_higher_orders:444"]], "{R_10: -25/2}")
  expect_identical(got[["column:test_reactions_higher_orders:455"]], "{R_0: -50/3, R_10: -100/3}")
  expect_identical(got[["column:test_reactions_higher_orders:465"]], "{R_10: -1000/3}")
  expect_identical(got[["column:test_reactions_higher_orders:482"]], "{R_0: -1088/5, R_10: -2516/15}")
  expect_identical(got[["column:test_reactions_higher_orders:502"]], "R_0*SingularityFunction(x, 0, -1) + R_8*SingularityFunction(x, 8, -1) - 60*SingularityFunction(x, 0, 0) + 10*SingularityFunction(x, 0, 1) - 10*SingularityFunction(x, 4, 1) - 100*SingularityFunction(x, 8, -1) + 20*SingularityFunction(x, 8, 0)")
})

test_that("column test_telescope_hinge matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
x <- Symbol("x")
})
.py_try({
x_1 <- Symbol("x_1")
})
.py_try({
.tmpu_ <- .py_symbols("E A"); E <- .py_getitem(.tmpu_, 0L); A <- .py_getitem(.tmpu_, 1L)
})
.py_try({
c <- Column(10L, E, A)
})
.py_try({
.py_method(c, "apply_support", 0L)
})
.py_try({
.py_method(c, "apply_support", 10L)
})
.py_try({
.py_method(c, "apply_load", 10L, 5L, -1L)
})
.py_try({
.py_method(c, "apply_telescope_hinge", .py_flt(7.5))
})
.py_try({
p <- .py_attr(c, "_bc_hinge")
})
.py_try({
q <- list(.py_flt(7.5))
})
.py_try({
.py_rec("column:test_telescope_hinge:521", function() p)
})
.py_try({
p <- .py_attr(c, "_applied_hinges")
})
.py_try({
q <- list(Symbol("u_7.5"))
})
.py_try({
.py_rec("column:test_telescope_hinge:525", function() p)
})
.py_try({
p <- .py_attr(c, "load")
})
.py_try({
.tmpu_ <- .py_symbols("R_0, R_10"); R_0 <- .py_getitem(.tmpu_, 0L); R_10 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
q <- .py_add(.py_add(.py_add(.py_mul(R_0, SingularityFunction(x, 0L, -1L)), .py_mul(10L, SingularityFunction(x, 5L, -1L))), .py_mul(.py_mul(.py_mul(E, A), Symbol("u_7.5")), SingularityFunction(x, .py_flt(7.5), -2L))), .py_mul(R_10, SingularityFunction(x, 10L, -1L)))
})
.py_try({
.py_rec("column:test_telescope_hinge:536", function() p)
})
.py_try({
.py_method(c, "solve_for_reaction_loads")
})
.py_try({
p <- .py_attr(c, "reaction_loads")
})
.py_try({
q <- .py_dict(list(R_0, R_10), list(-10L, 0L))
})
.py_try({
.py_rec("column:test_telescope_hinge:542", function() p)
})
.py_try({
p <- .py_attr(c, "hinge_extensions")
})
.py_try({
q <- .py_dict(list(Symbol("u_7.5")), list(.py_div(50L, .py_mul(E, A))))
})
.py_try({
.py_rec("column:test_telescope_hinge:546", function() p)
})
.py_try({
c2 <- Column(10L, 20000L, .py_flt(0.5))
})
.py_try({
.py_method(c2, "apply_support", 0L)
})
.py_try({
.py_method(c2, "apply_support", 10L)
})
.py_try({
.py_method(c2, "apply_telescope_hinge", 4L)
})
.py_try({
.py_method(c2, "apply_load", -5L, 3L, -1L)
})
.py_try({
.py_method(c2, "apply_load", -10L, 8L, -1L)
})
.py_try({
.py_method(c2, "solve_for_reaction_loads")
})
.py_try({
p <- .py_attr(c2, "reaction_loads")
})
.py_try({
q <- .py_dict(list(R_0, R_10), list(5L, 10L))
})
.py_try({
.py_rec("column:test_telescope_hinge:560", function() p)
})
.py_try({
p <- .py_attr(c2, "hinge_extensions")
})
.py_try({
u_4 <- Symbol("u_4")
})
.py_try({
q <- .py_dict(list(u_4), list(Rational(1L, 2000L)))
})
.py_try({
.py_rec("column:test_telescope_hinge:565", function() p)
}))---")
  expect_identical(got[["column:test_telescope_hinge:521"]], "[7.50000000000000]")
  expect_identical(got[["column:test_telescope_hinge:525"]], "[u_7.5]")
  expect_identical(got[["column:test_telescope_hinge:536"]], "A*E*u_7.5*SingularityFunction(x, 7.5, -2) + R_0*SingularityFunction(x, 0, -1) + R_10*SingularityFunction(x, 10, -1) + 10*SingularityFunction(x, 5, -1)")
  expect_identical(got[["column:test_telescope_hinge:542"]], "{R_0: -10, R_10: 0}")
  expect_identical(got[["column:test_telescope_hinge:546"]], "{u_7.5: 50/(A*E)}")
  expect_identical(got[["column:test_telescope_hinge:560"]], "{R_0: 5, R_10: 10}")
  expect_identical(got[["column:test_telescope_hinge:565"]], "{u_4: 1/2000}")
})

test_that("column test_equations matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
x <- Symbol("x")
})
.py_try({
x_1 <- Symbol("x_1")
})
.py_try({
c <- Column(10L, 210000L, 1L)
})
.py_try({
.py_method(c, "apply_support", 0L)
})
.py_try({
.py_method(c, "apply_support", 10L)
})
.py_try({
.py_method(c, "apply_load", 5L, 8L, -1L)
})
.py_try({
.tmpu_ <- .py_symbols("R_0 R_10"); R_0 <- .py_getitem(.tmpu_, 0L); R_10 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
.tmpu_ <- .py_symbols("C_N C_u"); C_N <- .py_getitem(.tmpu_, 0L); C_u <- .py_getitem(.tmpu_, 1L)
})
.py_try({
p <- .py_method(c, "axial_force")
})
.py_try({
q <- .py_sub(.py_sub(.py_sub(C_N, .py_mul(R_0, SingularityFunction(x, 0L, 0L))), .py_mul(R_10, SingularityFunction(x, 10L, 0L))), .py_mul(5L, SingularityFunction(x, 8L, 0L)))
})
.py_try({
.py_rec("column:test_equations:586", function() p)
})
.py_try({
p <- .py_method(c, "extension")
})
.py_try({
q <- .py_sub(.py_sub(.py_sub(.py_add(.py_mul(C_N, x), C_u), .py_div(.py_mul(R_0, SingularityFunction(x, 0L, 1L)), 210000L)), .py_div(.py_mul(R_10, SingularityFunction(x, 10L, 1L)), 210000L)), .py_div(SingularityFunction(x, 8L, 1L), 42000L))
})
.py_try({
.py_rec("column:test_equations:595", function() p)
})
.py_try({
.py_method(c, "solve_for_reaction_loads")
})
.py_try({
p <- .py_method(c, "axial_force")
})
.py_try({
q <- .py_add(.py_sub(SingularityFunction(x, 0L, 0L), .py_mul(5L, SingularityFunction(x, 8L, 0L))), .py_mul(4L, SingularityFunction(x, 10L, 0L)))
})
.py_try({
.py_rec("column:test_equations:606", function() p)
})
.py_try({
p <- .py_method(c, "extension")
})
.py_try({
q <- .py_add(.py_sub(.py_div(SingularityFunction(x, 0L, 1L), 210000L), .py_div(SingularityFunction(x, 8L, 1L), 42000L)), .py_div(SingularityFunction(x, 10L, 1L), 52500L))
})
.py_try({
.py_rec("column:test_equations:614", function() p)
}))---")
  expect_identical(got[["column:test_equations:586"]], "C_N - R_0*SingularityFunction(x, 0, 0) - R_10*SingularityFunction(x, 10, 0) - 5*SingularityFunction(x, 8, 0)")
  expect_identical(got[["column:test_equations:595"]], "C_N*x + C_u - R_0*SingularityFunction(x, 0, 1)/210000 - R_10*SingularityFunction(x, 10, 1)/210000 - SingularityFunction(x, 8, 1)/42000")
  expect_identical(got[["column:test_equations:606"]], "SingularityFunction(x, 0, 0) - 5*SingularityFunction(x, 8, 0) + 4*SingularityFunction(x, 10, 0)")
  expect_identical(got[["column:test_equations:614"]], "SingularityFunction(x, 0, 1)/210000 - SingularityFunction(x, 8, 1)/42000 + SingularityFunction(x, 10, 1)/52500")
})

test_that("truss test_truss matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
A <- Symbol("A")
})
.py_try({
B <- Symbol("B")
})
.py_try({
C <- Symbol("C")
})
.py_try({
.tmpu_ <- .py_symbols("AB, BC, AC"); AB <- .py_getitem(.tmpu_, 0L); BC <- .py_getitem(.tmpu_, 1L); AC <- .py_getitem(.tmpu_, 2L)
})
.py_try({
P <- Symbol("P")
})
.py_try({
t <- Truss()
})
.py_try({
.py_rec("truss:test_truss:15", function() .py_attr(t, "nodes"))
})
.py_try({
.py_rec("truss:test_truss:16", function() .py_attr(t, "node_labels"))
})
.py_try({
.py_rec("truss:test_truss:17", function() .py_attr(t, "node_positions"))
})
.py_try({
.py_rec("truss:test_truss:18", function() .py_attr(t, "members"))
})
.py_try({
.py_rec("truss:test_truss:19", function() .py_attr(t, "loads"))
})
.py_try({
.py_rec("truss:test_truss:20", function() .py_attr(t, "supports"))
})
.py_try({
.py_rec("truss:test_truss:21", function() .py_attr(t, "reaction_loads"))
})
.py_try({
.py_rec("truss:test_truss:22", function() .py_attr(t, "internal_forces"))
})
.py_try({
.py_method(t, "add_node", .py_tup(A, 0L, 0L), .py_tup(B, 2L, 2L), .py_tup(C, 3L, 0L))
})
.py_try({
.py_rec("truss:test_truss:26", function() .py_attr(t, "nodes"))
})
.py_try({
.py_rec("truss:test_truss:27", function() .py_attr(t, "node_labels"))
})
.py_try({
.py_rec("truss:test_truss:28", function() .py_attr(t, "node_positions"))
})
.py_try({
.py_rec("truss:test_truss:29", function() .py_attr(t, "loads"))
})
.py_try({
.py_rec("truss:test_truss:30", function() .py_attr(t, "supports"))
})
.py_try({
.py_rec("truss:test_truss:31", function() .py_attr(t, "reaction_loads"))
})
.py_try({
.py_method(t, "remove_node", C)
})
.py_try({
.py_rec("truss:test_truss:35", function() .py_attr(t, "nodes"))
})
.py_try({
.py_rec("truss:test_truss:36", function() .py_attr(t, "node_labels"))
})
.py_try({
.py_rec("truss:test_truss:37", function() .py_attr(t, "node_positions"))
})
.py_try({
.py_rec("truss:test_truss:38", function() .py_attr(t, "loads"))
})
.py_try({
.py_rec("truss:test_truss:39", function() .py_attr(t, "supports"))
})
.py_try({
.py_method(t, "add_node", .py_tup(C, 3L, 0L))
})
.py_try({
.py_method(t, "add_member", .py_tup(AB, A, B), .py_tup(BC, B, C), .py_tup(AC, A, C))
})
.py_try({
.py_rec("truss:test_truss:45", function() .py_attr(t, "members"))
})
.py_try({
.py_rec("truss:test_truss:46", function() .py_attr(t, "internal_forces"))
})
.py_try({
.py_method(t, "remove_member", BC)
})
.py_try({
.py_rec("truss:test_truss:50", function() .py_attr(t, "members"))
})
.py_try({
.py_rec("truss:test_truss:51", function() .py_attr(t, "internal_forces"))
})
.py_try({
.py_method(t, "add_member", .py_tup(BC, B, C))
})
.py_try({
.tmpu_ <- .py_symbols("D, CD"); D <- .py_getitem(.tmpu_, 0L); CD <- .py_getitem(.tmpu_, 1L)
})
.py_try({
.py_method(t, "change_node_label", .py_tup(B, D))
})
.py_try({
.py_rec("truss:test_truss:59", function() .py_attr(t, "nodes"))
})
.py_try({
.py_rec("truss:test_truss:60", function() .py_attr(t, "node_labels"))
})
.py_try({
.py_rec("truss:test_truss:61", function() .py_attr(t, "loads"))
})
.py_try({
.py_rec("truss:test_truss:62", function() .py_attr(t, "supports"))
})
.py_try({
.py_rec("truss:test_truss:63", function() .py_attr(t, "members"))
})
.py_try({
.py_method(t, "change_member_label", .py_tup(BC, CD))
})
.py_try({
.py_rec("truss:test_truss:66", function() .py_attr(t, "members"))
})
.py_try({
.py_rec("truss:test_truss:67", function() .py_attr(t, "internal_forces"))
})
.py_try({
.py_method(t, "apply_load", .py_tup(A, P, 90L), .py_tup(A, .py_div(P, 4L), 90L), .py_tup(A, .py_mul(2L, P), 45L), .py_tup(D, .py_div(P, 2L), 90L))
})
.py_try({
.py_rec("truss:test_truss:72", function() .py_attr(t, "loads"))
})
.py_try({
.py_rec("truss:test_truss:73", function() .py_getitem(.py_attr(t, "loads"), A))
})
.py_try({
.py_method(t, "remove_load", .py_tup(A, .py_div(P, 4L), 90L))
})
.py_try({
.py_rec("truss:test_truss:77", function() .py_attr(t, "loads"))
})
.py_try({
.py_rec("truss:test_truss:78", function() .py_getitem(.py_attr(t, "loads"), A))
})
.py_try({
.py_method(t, "apply_support", .py_tup(A, "pinned"), .py_tup(D, "roller"))
})
.py_try({
.py_rec("truss:test_truss:82", function() .py_attr(t, "supports"))
})
.py_try({
.py_rec("truss:test_truss:83", function() .py_attr(t, "reaction_loads"))
})
.py_try({
.py_rec("truss:test_truss:84", function() .py_attr(t, "loads"))
})
.py_try({
.py_method(t, "remove_support", A)
})
.py_try({
.py_rec("truss:test_truss:88", function() .py_attr(t, "supports"))
})
.py_try({
.py_rec("truss:test_truss:89", function() .py_attr(t, "reaction_loads"))
})
.py_try({
.py_rec("truss:test_truss:90", function() .py_attr(t, "loads"))
})
.py_try({
.py_method(t, "apply_support", .py_tup(A, "pinned"))
})
.py_try({
.py_method(t, "solve")
})
.py_try({
.py_rec("truss:test_truss:96", function() .py_getitem(.py_attr(t, "reaction_loads"), "R_A_x"))
})
.py_try({
.py_rec("truss:test_truss:97", function() .py_getitem(.py_attr(t, "reaction_loads"), "R_A_y"))
})
.py_try({
.py_rec("truss:test_truss:98", function() .py_getitem(.py_attr(t, "reaction_loads"), "R_D_y"))
})
.py_try({
.py_rec("truss:test_truss:99", function() .py_div(.py_getitem(.py_attr(t, "internal_forces"), AB), P))
})
.py_try({
.py_rec("truss:test_truss:100", function() .py_getitem(.py_attr(t, "internal_forces"), CD))
})
.py_try({
.py_rec("truss:test_truss:101", function() .py_getitem(.py_attr(t, "internal_forces"), AC))
}))---")
  expect_identical(got[["truss:test_truss:15"]], "[]")
  expect_identical(got[["truss:test_truss:16"]], "[]")
  expect_identical(got[["truss:test_truss:17"]], "[]")
  expect_identical(got[["truss:test_truss:18"]], "{}")
  expect_identical(got[["truss:test_truss:19"]], "{}")
  expect_identical(got[["truss:test_truss:20"]], "{}")
  expect_identical(got[["truss:test_truss:21"]], "{}")
  expect_identical(got[["truss:test_truss:22"]], "{}")
  expect_identical(got[["truss:test_truss:26"]], "[(A, 0, 0), (B, 2, 2), (C, 3, 0)]")
  expect_identical(got[["truss:test_truss:27"]], "[A, B, C]")
  expect_identical(got[["truss:test_truss:28"]], "[(0, 0), (2, 2), (3, 0)]")
  expect_identical(got[["truss:test_truss:29"]], "{}")
  expect_identical(got[["truss:test_truss:30"]], "{}")
  expect_identical(got[["truss:test_truss:31"]], "{}")
  expect_identical(got[["truss:test_truss:35"]], "[(A, 0, 0), (B, 2, 2)]")
  expect_identical(got[["truss:test_truss:36"]], "[A, B]")
  expect_identical(got[["truss:test_truss:37"]], "[(0, 0), (2, 2)]")
  expect_identical(got[["truss:test_truss:38"]], "{}")
  expect_identical(got[["truss:test_truss:39"]], "{}")
  expect_identical(got[["truss:test_truss:45"]], "{AB: [A, B], BC: [B, C], AC: [A, C]}")
  expect_identical(got[["truss:test_truss:46"]], "{AB: 0, BC: 0, AC: 0}")
  expect_identical(got[["truss:test_truss:50"]], "{AB: [A, B], AC: [A, C]}")
  expect_identical(got[["truss:test_truss:51"]], "{AB: 0, AC: 0}")
  expect_identical(got[["truss:test_truss:59"]], "[(A, 0, 0), (D, 2, 2), (C, 3, 0)]")
  expect_identical(got[["truss:test_truss:60"]], "[A, D, C]")
  expect_identical(got[["truss:test_truss:61"]], "{}")
  expect_identical(got[["truss:test_truss:62"]], "{}")
  expect_identical(got[["truss:test_truss:63"]], "{AB: [A, D], AC: [A, C], BC: [D, C]}")
  expect_identical(got[["truss:test_truss:66"]], "{AB: [A, D], AC: [A, C], CD: [D, C]}")
  expect_identical(got[["truss:test_truss:67"]], "{AB: 0, AC: 0, CD: 0}")
  expect_identical(got[["truss:test_truss:72"]], "{A: [[P, 90], [P/4, 90], [2*P, 45]], D: [[P/2, 90]]}")
  expect_identical(got[["truss:test_truss:73"]], "[[P, 90], [P/4, 90], [2*P, 45]]")
  expect_identical(got[["truss:test_truss:77"]], "{A: [[P, 90], [2*P, 45]], D: [[P/2, 90]]}")
  expect_identical(got[["truss:test_truss:78"]], "[[P, 90], [2*P, 45]]")
  expect_identical(got[["truss:test_truss:82"]], "{A: 'pinned', D: 'roller'}")
  expect_identical(got[["truss:test_truss:83"]], "{}")
  expect_identical(got[["truss:test_truss:84"]], "{A: [[P, 90], [2*P, 45], [R_A_x, 0], [R_A_y, 90]], D: [[P/2, 90], [R_D_y, 90]]}")
  expect_identical(got[["truss:test_truss:88"]], "{D: 'roller'}")
  expect_identical(got[["truss:test_truss:89"]], "{}")
  expect_identical(got[["truss:test_truss:90"]], "{A: [[P, 90], [2*P, 45]], D: [[P/2, 90], [R_D_y, 90]]}")
  expect_identical(got[["truss:test_truss:96"]], "-sqrt(2)*P")
  expect_identical(got[["truss:test_truss:97"]], "-sqrt(2)*P - P")
  expect_identical(got[["truss:test_truss:98"]], "-P/2")
  expect_identical(got[["truss:test_truss:99"]], "0")
  expect_identical(got[["truss:test_truss:100"]], "0")
  expect_identical(got[["truss:test_truss:101"]], "0")
})

test_that("cable test_cable matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
c <- Cable(.py_tup("A", 0L, 10L), .py_tup("B", 10L, 10L))
})
.py_try({
.py_rec("cable:test_cable:8", function() .py_attr(c, "supports"))
})
.py_try({
.py_rec("cable:test_cable:9", function() .py_attr(c, "left_support"))
})
.py_try({
.py_rec("cable:test_cable:10", function() .py_attr(c, "right_support"))
})
.py_try({
.py_rec("cable:test_cable:11", function() .py_attr(c, "loads"))
})
.py_try({
.py_rec("cable:test_cable:12", function() .py_attr(c, "loads_position"))
})
.py_try({
.py_rec("cable:test_cable:13", function() .py_attr(c, "length"))
})
.py_try({
.py_rec("cable:test_cable:14", function() .py_attr(c, "reaction_loads"))
})
.py_try({
.py_method(c, "change_support", "A", .py_tup("C", 12L, 3L))
})
.py_try({
.py_rec("cable:test_cable:18", function() .py_attr(c, "supports"))
})
.py_try({
.py_rec("cable:test_cable:19", function() .py_attr(c, "left_support"))
})
.py_try({
.py_rec("cable:test_cable:20", function() .py_attr(c, "right_support"))
})
.py_try({
.py_rec("cable:test_cable:21", function() .py_attr(c, "reaction_loads"))
})
.py_try({
.py_method(c, "change_support", "C", .py_tup("A", 0L, 10L))
})
.py_try({
.py_method(c, "apply_load", -1L, .py_tup("X", 2L, 5L, 3L, 30L))
})
.py_try({
.py_method(c, "apply_load", -1L, .py_tup("Y", 5L, 8L, 5L, 60L))
})
.py_try({
.py_rec("cable:test_cable:28", function() .py_attr(c, "loads"))
})
.py_try({
.py_rec("cable:test_cable:29", function() .py_attr(c, "loads_position"))
})
.py_try({
.py_rec("cable:test_cable:30", function() .py_attr(c, "length"))
})
.py_try({
.py_rec("cable:test_cable:31", function() .py_attr(c, "reaction_loads"))
})
.py_try({
.py_method(c, "remove_loads", "X")
})
.py_try({
.py_rec("cable:test_cable:35", function() .py_attr(c, "loads"))
})
.py_try({
.py_rec("cable:test_cable:36", function() .py_attr(c, "loads_position"))
})
.py_try({
.py_rec("cable:test_cable:37", function() .py_attr(c, "length"))
})
.py_try({
.py_rec("cable:test_cable:38", function() .py_attr(c, "reaction_loads"))
})
.py_try({
.py_method(c, "remove_loads", "Y")
})
.py_try({
.py_method(c, "apply_load", 0L, .py_tup("Z", 9L))
})
.py_try({
.py_rec("cable:test_cable:44", function() .py_attr(c, "loads"))
})
.py_try({
.py_rec("cable:test_cable:45", function() .py_attr(c, "loads_position"))
})
.py_try({
.py_rec("cable:test_cable:46", function() .py_attr(c, "length"))
})
.py_try({
.py_rec("cable:test_cable:47", function() .py_attr(c, "reaction_loads"))
})
.py_try({
.py_method(c, "apply_length", 20L)
})
.py_try({
.py_rec("cable:test_cable:51", function() .py_attr(c, "length"))
})
.py_try({
NULL
})
.py_try({
c <- Cable(.py_tup("A", 0L, 10L), .py_tup("B", .py_flt(5.5), 8L))
})
.py_try({
.py_method(c, "apply_load", -1L, .py_tup("Z", 2L, .py_flt(7.26), 3L, 270L))
})
.py_try({
.py_method(c, "apply_load", -1L, .py_tup("X", 4L, 6L, 8L, 270L))
})
.py_try({
.py_method(c, "solve")
})
.py_try({
.py_rec("cable:test_cable:61", function() .py_lt(.py_abs(.py_sub(.py_getitem(.py_attr(c, "tension"), Symbol("A_Z")), .py_flt(6.89488895397307))), .py_flt(1e-11)))
})
.py_try({
.py_rec("cable:test_cable:62", function() .py_lt(.py_abs(.py_sub(.py_getitem(.py_attr(c, "tension"), Symbol("Z_X")), .py_flt(4.79150773600774))), .py_flt(1e-11)))
})
.py_try({
.py_rec("cable:test_cable:63", function() .py_lt(.py_abs(.py_sub(.py_getitem(.py_attr(c, "tension"), Symbol("X_B")), .py_flt(6.78571428571429))), .py_flt(1e-11)))
})
.py_try({
.py_rec("cable:test_cable:65", function() .py_lt(.py_abs(.py_add(.py_getitem(.py_attr(c, "reaction_loads"), Symbol("R_A_x")), .py_flt(4.0650406504065))), .py_flt(1e-11)))
})
.py_try({
.py_rec("cable:test_cable:66", function() .py_lt(.py_abs(.py_sub(.py_getitem(.py_attr(c, "reaction_loads"), Symbol("R_A_y")), .py_flt(5.56910569105691))), .py_flt(1e-11)))
})
.py_try({
.py_rec("cable:test_cable:67", function() .py_lt(.py_abs(.py_sub(.py_getitem(.py_attr(c, "reaction_loads"), Symbol("R_B_x")), .py_flt(4.0650406504065))), .py_flt(1e-11)))
})
.py_try({
.py_rec("cable:test_cable:68", function() .py_lt(.py_abs(.py_sub(.py_getitem(.py_attr(c, "reaction_loads"), Symbol("R_B_y")), .py_flt(5.43089430894309))), .py_flt(1e-11)))
})
.py_try({
.py_rec("cable:test_cable:69", function() .py_lt(.py_abs(.py_sub(.py_attr(c, "length"), .py_flt(8.2560958484519))), .py_flt(1e-11)))
})
.py_try({
NULL
})
.py_try({
c <- Cable(.py_tup("A", 0L, 40L), .py_tup("B", 100L, 20L))
})
.py_try({
.py_method(c, "apply_load", 0L, .py_tup("X", 850L))
})
.py_try({
.py_method(c, "solve", .py_flt(58.58), 0L)
})
.py_try({
.py_rec("cable:test_cable:79", function() .py_lt(.py_abs(.py_sub(.py_method(c, "tension_at", 0L), .py_flt(61717.4130533677))), .py_flt(1e-10)))
})
.py_try({
.py_rec("cable:test_cable:80", function() .py_lt(.py_abs(.py_sub(.py_method(c, "tension_at", 40L), .py_flt(39738.0809048449))), .py_flt(1e-10)))
})
.py_try({
.py_rec("cable:test_cable:81", function() .py_lt(.py_abs(.py_sub(.py_getitem(.py_attr(c, "reaction_loads"), Symbol("R_A_x")), .py_flt(36465.0))), .py_flt(1e-10)))
})
.py_try({
.py_rec("cable:test_cable:82", function() .py_lt(.py_abs(.py_add(.py_getitem(.py_attr(c, "reaction_loads"), Symbol("R_A_y")), .py_flt(49793.0))), .py_flt(1e-10)))
})
.py_try({
.py_rec("cable:test_cable:83", function() .py_lt(.py_abs(.py_sub(.py_getitem(.py_attr(c, "reaction_loads"), Symbol("R_B_x")), .py_flt(44399.9537590861))), .py_flt(1e-10)))
})
.py_try({
.py_rec("cable:test_cable:84", function() .py_lt(.py_abs(.py_sub(.py_getitem(.py_attr(c, "reaction_loads"), Symbol("R_B_y")), .py_flt(42868.2071025955))), .py_flt(1e-10)))
}))---")
  expect_identical(got[["cable:test_cable:8"]], "{'A': [0, 10], 'B': [10, 10]}")
  expect_identical(got[["cable:test_cable:9"]], "[0, 10]")
  expect_identical(got[["cable:test_cable:10"]], "[10, 10]")
  expect_identical(got[["cable:test_cable:11"]], "{'distributed': {}, 'point_load': {}}")
  expect_identical(got[["cable:test_cable:12"]], "{}")
  expect_identical(got[["cable:test_cable:13"]], "0")
  expect_identical(got[["cable:test_cable:14"]], "{R_A_x: 0, R_A_y: 0, R_B_x: 0, R_B_y: 0}")
  expect_identical(got[["cable:test_cable:18"]], "{'B': [10, 10], 'C': [12, 3]}")
  expect_identical(got[["cable:test_cable:19"]], "[10, 10]")
  expect_identical(got[["cable:test_cable:20"]], "[12, 3]")
  expect_identical(got[["cable:test_cable:21"]], "{R_B_x: 0, R_B_y: 0, R_C_x: 0, R_C_y: 0}")
  expect_identical(got[["cable:test_cable:28"]], "{'distributed': {}, 'point_load': {'X': [3, 30], 'Y': [5, 60]}}")
  expect_identical(got[["cable:test_cable:29"]], "{'X': [2, 5], 'Y': [5, 8]}")
  expect_identical(got[["cable:test_cable:30"]], "0")
  expect_identical(got[["cable:test_cable:31"]], "{R_A_x: 0, R_A_y: 0, R_B_x: 0, R_B_y: 0}")
  expect_identical(got[["cable:test_cable:35"]], "{'distributed': {}, 'point_load': {'Y': [5, 60]}}")
  expect_identical(got[["cable:test_cable:36"]], "{'Y': [5, 8]}")
  expect_identical(got[["cable:test_cable:37"]], "0")
  expect_identical(got[["cable:test_cable:38"]], "{R_A_x: 0, R_A_y: 0, R_B_x: 0, R_B_y: 0}")
  expect_identical(got[["cable:test_cable:44"]], "{'distributed': {'Z': 9}, 'point_load': {}}")
  expect_identical(got[["cable:test_cable:45"]], "{}")
  expect_identical(got[["cable:test_cable:46"]], "0")
  expect_identical(got[["cable:test_cable:47"]], "{R_A_x: 0, R_A_y: 0, R_B_x: 0, R_B_y: 0}")
  expect_identical(got[["cable:test_cable:51"]], "20")
  expect_identical(got[["cable:test_cable:61"]], "True")
  expect_identical(got[["cable:test_cable:62"]], "True")
  expect_identical(got[["cable:test_cable:63"]], "True")
  expect_identical(got[["cable:test_cable:65"]], "True")
  expect_identical(got[["cable:test_cable:66"]], "True")
  expect_identical(got[["cable:test_cable:67"]], "True")
  expect_identical(got[["cable:test_cable:68"]], "True")
  expect_identical(got[["cable:test_cable:69"]], "True")
  expect_identical(got[["cable:test_cable:79"]], "True")
  expect_identical(got[["cable:test_cable:80"]], "True")
  expect_identical(got[["cable:test_cable:81"]], "True")
  expect_identical(got[["cable:test_cable:82"]], "True")
  expect_identical(got[["cable:test_cable:83"]], "True")
  expect_identical(got[["cable:test_cable:84"]], "True")
})

test_that("arch test_arch_init matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
x <- Symbol("x")
})
.py_try({
t <- Symbol("t")
})
.py_try({
a <- Arch(.py_tup(0L, 0L), .py_tup(10L, 0L), crown_x = 5L, crown_y = 5L)
})
.py_try({
.py_rec("arch:test_arch_init:10", function() .py_attr(a, "get_loads"))
})
.py_try({
.py_rec("arch:test_arch_init:11", function() .py_attr(a, "reaction_force"))
})
.py_try({
.py_rec("arch:test_arch_init:12", function() .py_attr(a, "supports"))
})
.py_try({
.py_rec("arch:test_arch_init:13", function() .py_attr(a, "left_support"))
})
.py_try({
.py_rec("arch:test_arch_init:14", function() .py_attr(a, "right_support"))
})
.py_try({
.py_rec("arch:test_arch_init:15", function() .py_attr(a, "get_shape_eqn"))
})
.py_try({
a <- Arch(.py_tup(0L, 0L), .py_tup(10L, 1L), crown_x = 6L)
})
.py_try({
.py_method(a, "change_support_type", left_support = "roller")
})
.py_try({
.py_method(a, "add_member", .py_flt(0.5))
})
.py_try({
.py_rec("arch:test_arch_init:20", function() .py_attr(a, "supports"))
})
.py_try({
.py_rec("arch:test_arch_init:21", function() simplify(.py_attr(a, "get_shape_eqn")))
}))---")
  expect_identical(got[["arch:test_arch_init:10"]], "{'concentrated': {}, 'distributed': {}}")
  expect_identical(got[["arch:test_arch_init:11"]], "{R_A_x: 0, R_A_y: 0, R_B_x: 0, R_B_y: 0}")
  expect_identical(got[["arch:test_arch_init:12"]], "{'left': 'hinge', 'right': 'hinge'}")
  expect_identical(got[["arch:test_arch_init:13"]], "(0, 0)")
  expect_identical(got[["arch:test_arch_init:14"]], "(10, 0)")
  expect_identical(got[["arch:test_arch_init:15"]], "5 - (x - 5)^2/5")
  expect_identical(got[["arch:test_arch_init:20"]], "{'left': 'roller', 'right': 'hinge'}")
  expect_identical(got[["arch:test_arch_init:21"]], "x*(12 - x)/20")
})

test_that("arch test_arch_support matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
x <- Symbol("x")
})
.py_try({
t <- Symbol("t")
})
.py_try({
a <- Arch(.py_tup(0L, 0L), .py_tup(40L, 0L), crown_x = 20L, crown_y = 12L)
})
.py_try({
.py_method(a, "apply_load", -1L, "C", 8L, 150L, angle = 270L)
})
.py_try({
.py_method(a, "apply_load", 0L, "D", start = 20L, end = 40L, mag = -4L)
})
.py_try({
.py_method(a, "solve")
})
.py_try({
.py_rec("arch:test_arch_support:28", function() .py_lt(.py_abs(.py_sub(.py_getitem(.py_attr(a, "reaction_force"), Symbol("R_A_x")), .py_flt(83.33333333333333))), .py_flt(1e-11)))
})
.py_try({
.py_rec("arch:test_arch_support:29", function() .py_lt(.py_abs(.py_sub(.py_getitem(.py_attr(a, "reaction_force"), Symbol("R_B_y")), .py_flt(90.0))), .py_flt(1e-11)))
})
.py_try({
.py_rec("arch:test_arch_support:30", function() .py_lt(.py_abs(.py_add(.py_getitem(.py_attr(a, "reaction_force"), Symbol("R_B_x")), .py_flt(83.33333333333333))), .py_flt(1e-11)))
})
.py_try({
.py_rec("arch:test_arch_support:31", function() .py_lt(.py_abs(.py_sub(.py_getitem(.py_attr(a, "reaction_force"), Symbol("R_A_y")), .py_flt(140.0))), .py_flt(1e-11)))
}))---")
  expect_identical(got[["arch:test_arch_support:28"]], "True")
  expect_identical(got[["arch:test_arch_support:29"]], "True")
  expect_identical(got[["arch:test_arch_support:30"]], "True")
  expect_identical(got[["arch:test_arch_support:31"]], "True")
})

test_that("arch test_arch_member matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
x <- Symbol("x")
})
.py_try({
t <- Symbol("t")
})
.py_try({
a <- Arch(.py_tup(0L, 0L), .py_tup(40L, 0L), crown_x = 20L, crown_y = 15L)
})
.py_try({
.py_method(a, "change_support_type", right_support = "roller")
})
.py_try({
.py_method(a, "add_member", 0L)
})
.py_try({
.py_method(a, "apply_load", -1L, "D", start = 12L, mag = 3L, angle = 270L)
})
.py_try({
.py_method(a, "apply_load", -1L, "E", start = 6L, mag = 4L, angle = 270L)
})
.py_try({
.py_method(a, "apply_load", -1L, "C", start = 30L, mag = 5L, angle = 270L)
})
.py_try({
.py_method(a, "solve")
})
.py_try({
.py_rec("arch:test_arch_member:41", function() .py_getitem(.py_attr(a, "reaction_force"), Symbol("R_A_x")))
})
.py_try({
.py_rec("arch:test_arch_member:42", function() .py_lt(.py_abs(.py_sub(.py_getitem(.py_attr(a, "reaction_force"), Symbol("R_A_y")), .py_flt(6.75))), .py_flt(1e-11)))
})
.py_try({
.py_rec("arch:test_arch_member:43", function() .py_getitem(.py_attr(a, "reaction_force"), Symbol("R_B_x")))
})
.py_try({
.py_rec("arch:test_arch_member:44", function() .py_lt(.py_abs(.py_sub(.py_getitem(.py_attr(a, "reaction_force"), Symbol("R_B_y")), .py_flt(5.25))), .py_flt(1e-11)))
}))---")
  expect_identical(got[["arch:test_arch_member:41"]], "0")
  expect_identical(got[["arch:test_arch_member:42"]], "True")
  expect_identical(got[["arch:test_arch_member:43"]], "0")
  expect_identical(got[["arch:test_arch_member:44"]], "True")
})

test_that("arch test_symbol_magnitude matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
x <- Symbol("x")
})
.py_try({
t <- Symbol("t")
})
.py_try({
a <- Arch(.py_tup(0L, 0L), .py_tup(16L, 0L), crown_x = 8L, crown_y = 5L)
})
.py_try({
.py_method(a, "apply_load", 0L, "C", start = 3L, end = 5L, mag = t)
})
.py_try({
.py_method(a, "solve")
})
.py_try({
.py_rec("arch:test_symbol_magnitude:50", function() .py_getitem(.py_attr(a, "reaction_force"), Symbol("R_A_x")))
})
.py_try({
.py_rec("arch:test_symbol_magnitude:51", function() .py_getitem(.py_attr(a, "reaction_force"), Symbol("R_A_y")))
})
.py_try({
.py_rec("arch:test_symbol_magnitude:52", function() .py_getitem(.py_attr(a, "reaction_force"), Symbol("R_B_x")))
})
.py_try({
.py_rec("arch:test_symbol_magnitude:53", function() .py_getitem(.py_attr(a, "reaction_force"), Symbol("R_B_y")))
})
.py_try({
.py_rec("arch:test_symbol_magnitude:54", function() .py_method(a, "bending_moment_at", 4L))
}))---")
  expect_identical(got[["arch:test_symbol_magnitude:50"]], "-4*t/5")
  expect_identical(got[["arch:test_symbol_magnitude:51"]], "-3*t/2")
  expect_identical(got[["arch:test_symbol_magnitude:52"]], "4*t/5")
  expect_identical(got[["arch:test_symbol_magnitude:53"]], "-t/2")
  expect_identical(got[["arch:test_symbol_magnitude:54"]], "-5*t/2")
})

test_that("arch test_forces matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
x <- Symbol("x")
})
.py_try({
t <- Symbol("t")
})
.py_try({
a <- Arch(.py_tup(0L, 0L), .py_tup(40L, 0L), crown_x = 20L, crown_y = 12L)
})
.py_try({
.py_method(a, "apply_load", -1L, "C", 8L, 150L, angle = 270L)
})
.py_try({
.py_method(a, "apply_load", 0L, "D", start = 20L, end = 40L, mag = -4L)
})
.py_try({
.py_method(a, "solve")
})
.py_try({
.py_rec("arch:test_forces:61", function() .py_lt(.py_abs(.py_sub(.py_method(a, "axial_force_at", .py_flt(7.999999999999999)), .py_flt(149.430523405935))), .py_flt(1e-12)))
})
.py_try({
.py_rec("arch:test_forces:62", function() .py_lt(.py_abs(.py_sub(.py_method(a, "shear_force_at", .py_flt(7.999999999999999)), .py_flt(64.9227473161196))), .py_flt(1e-12)))
}))---")
  expect_identical(got[["arch:test_forces:61"]], "True")
  expect_identical(got[["arch:test_forces:62"]], "True")
})
