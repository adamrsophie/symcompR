skip_unless_full_tests()

# Expected values are the outputs of the same calls in SymPy (physics.optics);
# the cases come from its test files and run through the helpers in
# helper-physics-misc-battery.R.

test_that("gaussopt test_gauss_opt matches SymPy", {
  got <- .physmisc_battery_run(r"---(streq <- function(a, b) {
  return(.py_eq(.py_str(a), .py_str(b)))
}
.py_try({
mat <- RayTransferMatrix(1L, 2L, 3L, 4L)
})
.py_try({
.py_rec("gaussopt:test_gauss_opt:23", function() mat)
})
.py_try({
.py_rec("gaussopt:test_gauss_opt:24", function() mat)
})
.py_try({
.py_rec("gaussopt:test_gauss_opt:25", function() list(.py_attr(mat, "A"), .py_attr(mat, "B"), .py_attr(mat, "C"), .py_attr(mat, "D")))
})
.py_try({
.tmpu_ <- .py_symbols("d f h n1 n2 R"); d <- .py_getitem(.tmpu_, 0L); f <- .py_getitem(.tmpu_, 1L); h <- .py_getitem(.tmpu_, 2L); n1 <- .py_getitem(.tmpu_, 3L); n2 <- .py_getitem(.tmpu_, 4L); R <- .py_getitem(.tmpu_, 5L)
})
.py_try({
lens <- ThinLens(f)
})
.py_try({
.py_rec("gaussopt:test_gauss_opt:29", function() lens)
})
.py_try({
.py_rec("gaussopt:test_gauss_opt:30", function() .py_attr(lens, "C"))
})
.py_try({
.py_rec("gaussopt:test_gauss_opt:31", function() FreeSpace(d))
})
.py_try({
.py_rec("gaussopt:test_gauss_opt:32", function() FlatRefraction(n1, n2))
})
.py_try({
.py_rec("gaussopt:test_gauss_opt:33", function() CurvedRefraction(R, n1, n2))
})
.py_try({
.py_rec("gaussopt:test_gauss_opt:35", function() FlatMirror())
})
.py_try({
.py_rec("gaussopt:test_gauss_opt:36", function() CurvedMirror(R))
})
.py_try({
.py_rec("gaussopt:test_gauss_opt:37", function() ThinLens(f))
})
.py_try({
mul <- .py_mul(CurvedMirror(R), FreeSpace(d))
})
.py_try({
mul_mat <- .py_mul(.py_Matrix(list(list(1L, 0L), list(.py_div(-2L, R), 1L))), .py_Matrix(list(list(1L, d), list(0L, 1L))))
})
.py_try({
.py_rec("gaussopt:test_gauss_opt:41", function() .py_attr(mul, "A"))
})
.py_try({
.py_rec("gaussopt:test_gauss_opt:42", function() .py_attr(mul, "B"))
})
.py_try({
.py_rec("gaussopt:test_gauss_opt:43", function() .py_attr(mul, "C"))
})
.py_try({
.py_rec("gaussopt:test_gauss_opt:44", function() .py_attr(mul, "D"))
})
.py_try({
angle <- .py_symbols("angle")
})
.py_try({
.py_rec("gaussopt:test_gauss_opt:47", function() GeometricRay(h, angle))
})
.py_try({
.py_rec("gaussopt:test_gauss_opt:48", function() .py_mul(FreeSpace(d), GeometricRay(h, angle)))
})
.py_try({
.py_rec("gaussopt:test_gauss_opt:50", function() GeometricRay(.py_Matrix(.py_tup(.py_tup(h), .py_tup(angle)))))
})
.py_try({
.py_rec("gaussopt:test_gauss_opt:51", function() .py_attr(.py_mul(FreeSpace(d), GeometricRay(h, angle)), "height"))
})
.py_try({
.py_rec("gaussopt:test_gauss_opt:52", function() .py_attr(.py_mul(FreeSpace(d), GeometricRay(h, angle)), "angle"))
})
.py_try({
p <- BeamParameter(.py_flt(5.3e-07), 1L, w = .py_flt(0.001))
})
.py_try({
.py_rec("gaussopt:test_gauss_opt:55", function() streq(.py_attr(p, "q"), .py_add(1L, .py_mul(.py_mul(.py_flt(1.88679245283019), sym_I), S_Pi))))
})
.py_try({
.py_rec("gaussopt:test_gauss_opt:56", function() streq(.py_N(.py_attr(p, "q")), .py_add(.py_flt(1.0), .py_mul(.py_flt(5.92753330865999), sym_I))))
})
.py_try({
.py_rec("gaussopt:test_gauss_opt:57", function() streq(.py_N(.py_attr(p, "w_0")), Float(.py_flt(0.001))))
})
.py_try({
.py_rec("gaussopt:test_gauss_opt:58", function() streq(.py_N(.py_attr(p, "z_r")), Float(.py_flt(5.92753330865999))))
})
.py_try({
fs <- FreeSpace(10L)
})
.py_try({
p1 <- .py_mul(fs, p)
})
.py_try({
.py_rec("gaussopt:test_gauss_opt:61", function() streq(.py_N(.py_attr(p, "w")), Float(.py_flt(0.00101413072159615))))
})
.py_try({
.py_rec("gaussopt:test_gauss_opt:62", function() streq(.py_N(.py_attr(p1, "w")), Float(.py_flt(0.00210803120913829))))
})
.py_try({
.tmpu_ <- .py_symbols("w wavelen"); w <- .py_getitem(.tmpu_, 0L); wavelen <- .py_getitem(.tmpu_, 1L)
})
.py_try({
.py_rec("gaussopt:test_gauss_opt:65", function() waist2rayleigh(w, wavelen))
})
.py_try({
.tmpu_ <- .py_symbols("z_r wavelen"); z_r <- .py_getitem(.tmpu_, 0L); wavelen <- .py_getitem(.tmpu_, 1L)
})
.py_try({
.py_rec("gaussopt:test_gauss_opt:67", function() rayleigh2waist(z_r, wavelen))
})
.py_try({
.tmpu_ <- .py_symbols("a b f"); a <- .py_getitem(.tmpu_, 0L); b <- .py_getitem(.tmpu_, 1L); f <- .py_getitem(.tmpu_, 2L)
})
.py_try({
.py_rec("gaussopt:test_gauss_opt:70", function() geometric_conj_ab(a, b))
})
.py_try({
.py_rec("gaussopt:test_gauss_opt:71", function() geometric_conj_af(a, f))
})
.py_try({
.py_rec("gaussopt:test_gauss_opt:72", function() geometric_conj_bf(b, f))
})
.py_try({
.py_rec("gaussopt:test_gauss_opt:73", function() geometric_conj_ab(oo, b))
})
.py_try({
.py_rec("gaussopt:test_gauss_opt:74", function() geometric_conj_ab(a, oo))
})
.py_try({
.tmpu_ <- .py_symbols("s_in z_r_in f"); s_in <- .py_getitem(.tmpu_, 0L); z_r_in <- .py_getitem(.tmpu_, 1L); f <- .py_getitem(.tmpu_, 2L)
})
.py_try({
.py_rec("gaussopt:test_gauss_opt:77", function() .py_getitem(gaussian_conj(s_in, z_r_in, f), 0L))
})
.py_try({
.py_rec("gaussopt:test_gauss_opt:79", function() .py_getitem(gaussian_conj(s_in, z_r_in, f), 1L))
})
.py_try({
.py_rec("gaussopt:test_gauss_opt:81", function() .py_getitem(gaussian_conj(s_in, z_r_in, f), 2L))
})
.py_try({
.tmpu_ <- .py_symbols("l w_i w_o f"); l <- .py_getitem(.tmpu_, 0L); w_i <- .py_getitem(.tmpu_, 1L); w_o <- .py_getitem(.tmpu_, 2L); f <- .py_getitem(.tmpu_, 3L)
})
.py_try({
.py_rec("gaussopt:test_gauss_opt:85", function() .py_getitem(conjugate_gauss_beams(l, w_i, w_o, f = f), 0L))
})
.py_try({
.py_rec("gaussopt:test_gauss_opt:87", function() sym_factor(.py_getitem(conjugate_gauss_beams(l, w_i, w_o, f = f), 1L)))
})
.py_try({
.py_rec("gaussopt:test_gauss_opt:89", function() .py_getitem(conjugate_gauss_beams(l, w_i, w_o, f = f), 2L))
})
.py_try({
.tmpu_ <- .py_symbols("z l w_0", positive = TRUE); z <- .py_getitem(.tmpu_, 0L); l <- .py_getitem(.tmpu_, 1L); w_0 <- .py_getitem(.tmpu_, 2L)
})
.py_try({
p <- BeamParameter(l, z, w = w_0)
})
.py_try({
.py_rec("gaussopt:test_gauss_opt:93", function() .py_attr(p, "radius"))
})
.py_try({
.py_rec("gaussopt:test_gauss_opt:94", function() .py_attr(p, "w"))
})
.py_try({
.py_rec("gaussopt:test_gauss_opt:95", function() .py_attr(p, "w_0"))
})
.py_try({
.py_rec("gaussopt:test_gauss_opt:96", function() .py_attr(p, "divergence"))
})
.py_try({
.py_rec("gaussopt:test_gauss_opt:97", function() .py_attr(p, "gouy"))
})
.py_try({
.py_rec("gaussopt:test_gauss_opt:98", function() .py_attr(p, "waist_approximation_limit"))
})
.py_try({
p <- BeamParameter(.py_flt(5.3e-07), 1L, w = .py_flt(0.001), n = 2L)
})
.py_try({
.py_rec("gaussopt:test_gauss_opt:101", function() streq(.py_attr(p, "q"), .py_add(1L, .py_mul(.py_mul(.py_flt(3.77358490566038), sym_I), S_Pi))))
})
.py_try({
.py_rec("gaussopt:test_gauss_opt:102", function() streq(.py_N(.py_attr(p, "z_r")), Float(.py_flt(11.85506661732))))
})
.py_try({
.py_rec("gaussopt:test_gauss_opt:103", function() streq(.py_N(.py_attr(p, "w_0")), Float(.py_flt(0.001))))
}))---")
  expect_identical(got[["gaussopt:test_gauss_opt:23"]], "Matrix([[1, 2], [3, 4]])")
  expect_identical(got[["gaussopt:test_gauss_opt:24"]], "Matrix([[1, 2], [3, 4]])")
  expect_identical(got[["gaussopt:test_gauss_opt:25"]], "[1, 2, 3, 4]")
  expect_identical(got[["gaussopt:test_gauss_opt:29"]], "Matrix([[1, 0], [-1/f, 1]])")
  expect_identical(got[["gaussopt:test_gauss_opt:30"]], "-1/f")
  expect_identical(got[["gaussopt:test_gauss_opt:31"]], "Matrix([[1, d], [0, 1]])")
  expect_identical(got[["gaussopt:test_gauss_opt:32"]], "Matrix([[1, 0], [0, n1/n2]])")
  expect_identical(got[["gaussopt:test_gauss_opt:33"]], "Matrix([[1, 0], [(n1 - n2)/(R*n2), n1/n2]])")
  expect_identical(got[["gaussopt:test_gauss_opt:35"]], "Matrix([[1, 0], [0, 1]])")
  expect_identical(got[["gaussopt:test_gauss_opt:36"]], "Matrix([[1, 0], [-2/R, 1]])")
  expect_identical(got[["gaussopt:test_gauss_opt:37"]], "Matrix([[1, 0], [-1/f, 1]])")
  expect_identical(got[["gaussopt:test_gauss_opt:41"]], "1")
  expect_identical(got[["gaussopt:test_gauss_opt:42"]], "d")
  expect_identical(got[["gaussopt:test_gauss_opt:43"]], "-2/R")
  expect_identical(got[["gaussopt:test_gauss_opt:44"]], "1 - 2*d/R")
  expect_identical(got[["gaussopt:test_gauss_opt:47"]], "Matrix([[h], [angle]])")
  expect_identical(got[["gaussopt:test_gauss_opt:48"]], "Matrix([[angle*d + h], [angle]])")
  expect_identical(got[["gaussopt:test_gauss_opt:50"]], "Matrix([[h], [angle]])")
  expect_identical(got[["gaussopt:test_gauss_opt:51"]], "angle*d + h")
  expect_identical(got[["gaussopt:test_gauss_opt:52"]], "angle")
  expect_identical(got[["gaussopt:test_gauss_opt:55"]], "True")
  expect_identical(got[["gaussopt:test_gauss_opt:56"]], "True")
  expect_identical(got[["gaussopt:test_gauss_opt:57"]], "True")
  expect_identical(got[["gaussopt:test_gauss_opt:58"]], "True")
  expect_identical(got[["gaussopt:test_gauss_opt:61"]], "True")
  expect_identical(got[["gaussopt:test_gauss_opt:62"]], "True")
  expect_identical(got[["gaussopt:test_gauss_opt:65"]], "pi*w^2/wavelen")
  expect_identical(got[["gaussopt:test_gauss_opt:67"]], "sqrt(wavelen*z_r)/sqrt(pi)")
  expect_identical(got[["gaussopt:test_gauss_opt:70"]], "a*b/(a + b)")
  expect_identical(got[["gaussopt:test_gauss_opt:71"]], "a*f/(a - f)")
  expect_identical(got[["gaussopt:test_gauss_opt:72"]], "b*f/(b - f)")
  expect_identical(got[["gaussopt:test_gauss_opt:73"]], "b")
  expect_identical(got[["gaussopt:test_gauss_opt:74"]], "a")
  expect_identical(got[["gaussopt:test_gauss_opt:77"]], "1/(-1/(s_in + z_r_in^2/(-f + s_in)) + 1/f)")
  expect_identical(got[["gaussopt:test_gauss_opt:79"]], "z_r_in/(1 - s_in^2/f^2 + z_r_in^2/f^2)")
  expect_identical(got[["gaussopt:test_gauss_opt:81"]], "1/sqrt(1 - s_in^2/f^2 + z_r_in^2/f^2)")
  expect_identical(got[["gaussopt:test_gauss_opt:85"]], "f*(1 - sqrt(w_i^2/w_o^2 - pi^2*w_i^4/(f^2*l^2)))")
  expect_identical(got[["gaussopt:test_gauss_opt:87"]], "f*w_o^2*(w_i^2/w_o^2 - sqrt(w_i^2/w_o^2 - pi^2*w_i^4/(f^2*l^2)))/w_i^2")
  expect_identical(got[["gaussopt:test_gauss_opt:89"]], "f")
  expect_identical(got[["gaussopt:test_gauss_opt:93"]], "z*(1 + pi^2*w_0^4/(l^2*z^2))")
  expect_identical(got[["gaussopt:test_gauss_opt:94"]], "w_0*sqrt(l^2*z^2/(pi^2*w_0^4) + 1)")
  expect_identical(got[["gaussopt:test_gauss_opt:95"]], "w_0")
  expect_identical(got[["gaussopt:test_gauss_opt:96"]], "l/(pi*w_0)")
  expect_identical(got[["gaussopt:test_gauss_opt:97"]], "atan(l*z/(pi*w_0^2))")
  expect_identical(got[["gaussopt:test_gauss_opt:98"]], "2*l/pi")
  expect_identical(got[["gaussopt:test_gauss_opt:101"]], "True")
  expect_identical(got[["gaussopt:test_gauss_opt:102"]], "True")
  expect_identical(got[["gaussopt:test_gauss_opt:103"]], "True")
})

test_that("medium test_medium matches SymPy", {
  got <- .physmisc_battery_run(r"---(epsilon <- Symbol("epsilon")
mu <- Symbol("mu")
n <- Symbol("n")
.py_try({
c <- .py_method(.py_unit("speed_of_light"), "convert_to", .py_div(.py_unit("m"), .py_unit("s")))
})
.py_try({
e0 <- .py_method(.py_unit("e0"), "convert_to", .py_div(.py_mul(.py_pow(.py_unit("A"), 2L), .py_pow(.py_unit("s"), 4L)), .py_mul(.py_unit("kg"), .py_pow(.py_unit("m"), 3L))))
})
.py_try({
u0 <- .py_method(.py_unit("u0"), "convert_to", .py_div(.py_mul(.py_unit("m"), .py_unit("kg")), .py_mul(.py_pow(.py_unit("A"), 2L), .py_pow(.py_unit("s"), 2L))))
})
.py_try({
m1 <- Medium("m1")
})
.py_try({
.py_rec("medium:test_medium:16", function() .py_attr(m1, "intrinsic_impedance"))
})
.py_try({
.py_rec("medium:test_medium:17", function() .py_attr(m1, "speed"))
})
.py_try({
.py_rec("medium:test_medium:18", function() .py_attr(m1, "refractive_index"))
})
.py_try({
.py_rec("medium:test_medium:19", function() .py_attr(m1, "permittivity"))
})
.py_try({
.py_rec("medium:test_medium:20", function() .py_attr(m1, "permeability"))
})
.py_try({
m2 <- Medium("m2", epsilon, mu)
})
.py_try({
.py_rec("medium:test_medium:22", function() .py_attr(m2, "intrinsic_impedance"))
})
.py_try({
.py_rec("medium:test_medium:23", function() .py_attr(m2, "speed"))
})
.py_try({
.py_rec("medium:test_medium:24", function() .py_attr(m2, "refractive_index"))
})
.py_try({
.py_rec("medium:test_medium:25", function() .py_attr(m2, "permittivity"))
})
.py_try({
.py_rec("medium:test_medium:26", function() .py_attr(m2, "permeability"))
})
.py_try({
m3 <- Medium("m3", .py_div(.py_mul(.py_mul(.py_mul(.py_flt(9.0), .py_pow(10L, -12L)), .py_pow(.py_unit("s"), 4L)), .py_pow(.py_unit("A"), 2L)), .py_mul(.py_pow(.py_unit("m"), 3L), .py_unit("kg"))), .py_div(.py_mul(.py_mul(.py_mul(.py_flt(1.45), .py_pow(10L, -6L)), .py_unit("kg")), .py_unit("m")), .py_mul(.py_pow(.py_unit("A"), 2L), .py_pow(.py_unit("s"), 2L))))
})
.py_try({
.py_rec("medium:test_medium:30", function() .py_gt(.py_attr(m3, "refractive_index"), .py_attr(m1, "refractive_index")))
})
.py_try({
.py_rec("medium:test_medium:31", function() m3)
})
.py_try({
m4 <- Medium("m4", .py_div(.py_mul(.py_mul(.py_mul(.py_flt(7.0), .py_pow(10L, -12L)), .py_pow(.py_unit("s"), 4L)), .py_pow(.py_unit("A"), 2L)), .py_mul(.py_pow(.py_unit("m"), 3L), .py_unit("kg"))), .py_div(.py_mul(.py_mul(.py_mul(.py_flt(1.15), .py_pow(10L, -6L)), .py_unit("kg")), .py_unit("m")), .py_mul(.py_pow(.py_unit("A"), 2L), .py_pow(.py_unit("s"), 2L))))
})
.py_try({
.py_rec("medium:test_medium:35", function() .py_lt(.py_attr(m4, "refractive_index"), .py_attr(m1, "refractive_index")))
})
.py_try({
m5 <- Medium("m5", permittivity = .py_div(.py_mul(.py_mul(.py_mul(710L, .py_pow(10L, -12L)), .py_pow(.py_unit("s"), 4L)), .py_pow(.py_unit("A"), 2L)), .py_mul(.py_pow(.py_unit("m"), 3L), .py_unit("kg"))), n = .py_flt(1.33))
})
.py_try({
.py_rec("medium:test_medium:37", function() .py_lt(.py_abs(.py_sub(.py_attr(m5, "intrinsic_impedance"), .py_div(.py_mul(.py_mul(.py_flt(6.24845417765552), .py_unit("kg")), .py_pow(.py_unit("m"), 2L)), .py_mul(.py_pow(.py_unit("A"), 2L), .py_pow(.py_unit("s"), 3L))))), .py_div(.py_mul(.py_mul(.py_flt(1e-12), .py_unit("kg")), .py_pow(.py_unit("m"), 2L)), .py_mul(.py_pow(.py_unit("A"), 2L), .py_pow(.py_unit("s"), 3L)))))
})
.py_try({
.py_rec("medium:test_medium:39", function() .py_lt(.py_abs(.py_sub(.py_attr(m5, "speed"), .py_div(.py_mul(.py_flt(225407863.157895), .py_unit("m")), .py_unit("s")))), .py_div(.py_mul(.py_flt(1e-06), .py_unit("m")), .py_unit("s"))))
})
.py_try({
.py_rec("medium:test_medium:40", function() .py_lt(.py_abs(.py_sub(.py_attr(m5, "refractive_index"), .py_flt(1.33))), .py_flt(1e-12)))
})
.py_try({
.py_rec("medium:test_medium:41", function() .py_lt(.py_abs(.py_sub(.py_attr(m5, "permittivity"), .py_div(.py_mul(.py_mul(.py_flt(7.1e-10), .py_pow(.py_unit("A"), 2L)), .py_pow(.py_unit("s"), 4L)), .py_mul(.py_unit("kg"), .py_pow(.py_unit("m"), 3L))))), .py_div(.py_mul(.py_mul(.py_flt(1e-20), .py_pow(.py_unit("A"), 2L)), .py_pow(.py_unit("s"), 4L)), .py_mul(.py_unit("kg"), .py_pow(.py_unit("m"), 3L)))))
})
.py_try({
.py_rec("medium:test_medium:43", function() .py_lt(.py_abs(.py_sub(.py_attr(m5, "permeability"), .py_div(.py_mul(.py_mul(.py_flt(2.77206575232851e-08), .py_unit("kg")), .py_unit("m")), .py_mul(.py_pow(.py_unit("A"), 2L), .py_pow(.py_unit("s"), 2L))))), .py_div(.py_mul(.py_mul(.py_flt(1e-20), .py_unit("kg")), .py_unit("m")), .py_mul(.py_pow(.py_unit("A"), 2L), .py_pow(.py_unit("s"), 2L)))))
})
.py_try({
m6 <- Medium("m6", NULL, mu, n)
})
.py_try({
.py_rec("medium:test_medium:46", function() .py_attr(m6, "permittivity"))
})
.py_try({
.py_rec("medium:test_medium:48", function() .py_attr(Medium("m7"), "refractive_index"))
})
.py_try({
.py_rec("medium:test_medium:49", function() Medium("m9", e0, u0, 2L))
}))---")
  expect_identical(got[["medium:test_medium:16"]], "149896229*pi*kilogram*meter^2/(1250000*ampere^2*second^3)")
  expect_identical(got[["medium:test_medium:17"]], "299792458*meter/second")
  expect_identical(got[["medium:test_medium:18"]], "1")
  expect_identical(got[["medium:test_medium:19"]], "625000*ampere^2*second^4/(22468879468420441*pi*kilogram*meter^3)")
  expect_identical(got[["medium:test_medium:20"]], "pi*kilogram*meter/(2500000*ampere^2*second^2)")
  expect_identical(got[["medium:test_medium:22"]], "sqrt(mu/epsilon)")
  expect_identical(got[["medium:test_medium:23"]], "1/sqrt(epsilon*mu)")
  expect_identical(got[["medium:test_medium:24"]], "299792458*meter*sqrt(epsilon*mu)/second")
  expect_identical(got[["medium:test_medium:25"]], "epsilon")
  expect_identical(got[["medium:test_medium:26"]], "mu")
  expect_identical(got[["medium:test_medium:30"]], "True")
  expect_identical(got[["medium:test_medium:31"]], "MediumPP(m3, 9.0e-12*ampere^2*second^4/(kilogram*meter^3), 1.45e-6*kilogram*meter/(ampere^2*second^2))")
  expect_identical(got[["medium:test_medium:35"]], "True")
  expect_identical(got[["medium:test_medium:37"]], "True")
  expect_identical(got[["medium:test_medium:39"]], "True")
  expect_identical(got[["medium:test_medium:40"]], "True")
  expect_identical(got[["medium:test_medium:41"]], "True")
  expect_identical(got[["medium:test_medium:43"]], "True")
  expect_identical(got[["medium:test_medium:46"]], "second^2*n^2/(89875517873681764*meter^2*mu)")
  expect_identical(got[["medium:test_medium:48"]], "1")
  expect_identical(got[["medium:test_medium:49"]], "ERR ValueError: Specifying all of permittivity, permeability, and n is not allowed")
})

test_that("polarization test_polarization matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.py_rec("polarization:test_polarization:14", function() jones_vector(0L, 0L))
})
.py_try({
.py_rec("polarization:test_polarization:15", function() jones_vector(.py_div(S_Pi, 2L), 0L))
})
.py_try({
.py_rec("polarization:test_polarization:17", function() stokes_vector(0L, 0L))
})
.py_try({
.py_rec("polarization:test_polarization:18", function() stokes_vector(.py_div(S_Pi, 2L), 0L))
})
.py_try({
H <- jones_vector(0L, 0L)
})
.py_try({
V <- jones_vector(.py_div(S_Pi, 2L), 0L)
})
.py_try({
D <- jones_vector(.py_div(S_Pi, 4L), 0L)
})
.py_try({
A <- jones_vector(.py_div(.py_neg(S_Pi), 4L), 0L)
})
.py_try({
R <- jones_vector(0L, .py_div(S_Pi, 4L))
})
.py_try({
L <- jones_vector(0L, .py_div(.py_neg(S_Pi), 4L))
})
.py_try({
res <- list(.py_Matrix(list(1L, 1L, 0L, 0L)), .py_Matrix(list(1L, -1L, 0L, 0L)), .py_Matrix(list(1L, 0L, 1L, 0L)), .py_Matrix(list(1L, 0L, -1L, 0L)), .py_Matrix(list(1L, 0L, 0L, 1L)), .py_Matrix(list(1L, 0L, 0L, -1L)))
})
.py_try({
.py_rec("polarization:test_polarization:34", function() .py_comp(list(H, V, D, A, R, L), function(.it) { e <- .it; list(jones_2_stokes(e)) }))
})
.py_try({
.py_rec("polarization:test_polarization:36", function() linear_polarizer(0L))
})
.py_try({
delta <- .py_symbols("delta", real = TRUE)
})
.py_try({
res <- .py_Matrix(list(list(.py_exp(.py_div(.py_mul(.py_neg(sym_I), delta), 2L)), 0L), list(0L, .py_exp(.py_div(.py_mul(sym_I, delta), 2L)))))
})
.py_try({
.py_rec("polarization:test_polarization:40", function() phase_retarder(0L, delta))
})
.py_try({
.py_rec("polarization:test_polarization:42", function() half_wave_retarder(0L))
})
.py_try({
res <- .py_Matrix(list(list(.py_exp(.py_div(.py_mul(.py_neg(sym_I), S_Pi), 4L)), 0L), list(0L, .py_mul(sym_I, .py_exp(.py_div(.py_mul(.py_neg(sym_I), S_Pi), 4L))))))
})
.py_try({
.py_rec("polarization:test_polarization:45", function() quarter_wave_retarder(0L))
})
.py_try({
.py_rec("polarization:test_polarization:47", function() transmissive_filter(1L))
})
.py_try({
.py_rec("polarization:test_polarization:49", function() reflective_filter(1L))
})
.py_try({
res <- .py_Matrix(list(list(.py_div(.py_S(1L), 2L), .py_div(.py_S(1L), 2L), 0L, 0L), list(.py_div(.py_S(1L), 2L), .py_div(.py_S(1L), 2L), 0L, 0L), list(0L, 0L, 0L, 0L), list(0L, 0L, 0L, 0L)))
})
.py_try({
.py_rec("polarization:test_polarization:55", function() mueller_matrix(linear_polarizer(0L)))
})
.py_try({
res <- .py_Matrix(list(list(1L, 0L, 0L, 0L), list(0L, 0L, 0L, .py_neg(sym_I)), list(0L, 0L, 1L, 0L), list(0L, .py_neg(sym_I), 0L, 0L)))
})
.py_try({
.py_rec("polarization:test_polarization:58", function() polarizing_beam_splitter())
}))---")
  expect_identical(got[["polarization:test_polarization:14"]], "Matrix([[1], [0]])")
  expect_identical(got[["polarization:test_polarization:15"]], "Matrix([[0], [1]])")
  expect_identical(got[["polarization:test_polarization:17"]], "Matrix([[1], [1], [0], [0]])")
  expect_identical(got[["polarization:test_polarization:18"]], "Matrix([[1], [-1], [0], [0]])")
  expect_identical(got[["polarization:test_polarization:34"]], "[Matrix([[1], [1], [0], [0]]), Matrix([[1], [-1], [0], [0]]), Matrix([[1], [0], [1], [0]]), Matrix([[1], [0], [-1], [0]]), Matrix([[1], [0], [0], [1]]), Matrix([[1], [0], [0], [-1]])]")
  expect_identical(got[["polarization:test_polarization:36"]], "Matrix([[1, 0], [0, 0]])")
  expect_identical(got[["polarization:test_polarization:40"]], "Matrix([[exp(-I*delta/2), 0], [0, exp(I*delta/2)]])")
  expect_identical(got[["polarization:test_polarization:42"]], "Matrix([[-I, 0], [0, I]])")
  expect_identical(got[["polarization:test_polarization:45"]], "Matrix([[exp(-I*pi/4), 0], [0, I*exp(-I*pi/4)]])")
  expect_identical(got[["polarization:test_polarization:47"]], "Matrix([[1, 0], [0, 1]])")
  expect_identical(got[["polarization:test_polarization:49"]], "Matrix([[1, 0], [0, -1]])")
  expect_identical(got[["polarization:test_polarization:55"]], "Matrix([[1/2, 1/2, 0, 0], [1/2, 1/2, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])")
  expect_identical(got[["polarization:test_polarization:58"]], "Matrix([[1, 0, 0, 0], [0, 0, 0, -I], [0, 0, 1, 0], [0, -I, 0, 0]])")
})

test_that("optics_utils test_refraction_angle matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
ae <- function(a, b, n) comp(a, b, .py_pow(10L, .py_neg(n)))
})
.py_try({
.tmpu_ <- .py_symbols("n1, n2"); n1 <- .py_getitem(.tmpu_, 0L); n2 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
m1 <- Medium("m1")
})
.py_try({
m2 <- Medium("m2")
})
.py_try({
r1 <- Ray3D(Point3D(-1L, -1L, 1L), Point3D(0L, 0L, 0L))
})
.py_try({
i <- .py_Matrix(list(1L, 1L, 1L))
})
.py_try({
n <- .py_Matrix(list(0L, 0L, 1L))
})
.py_try({
normal_ray <- Ray3D(Point3D(0L, 0L, 0L), Point3D(0L, 0L, 1L))
})
.py_try({
P <- Plane(Point3D(0L, 0L, 0L), normal_vector = list(0L, 0L, 1L))
})
.py_try({
.py_rec("optics_utils:test_refraction_angle:33", function() refraction_angle(r1, 1L, 1L, n))
})
.py_try({
.py_rec("optics_utils:test_refraction_angle:37", function() refraction_angle(list(1L, 1L, 1L), 1L, 1L, n))
})
.py_try({
.py_rec("optics_utils:test_refraction_angle:41", function() refraction_angle(.py_tup(1L, 1L, 1L), 1L, 1L, n))
})
.py_try({
.py_rec("optics_utils:test_refraction_angle:45", function() refraction_angle(i, 1L, 1L, list(0L, 0L, 1L)))
})
.py_try({
.py_rec("optics_utils:test_refraction_angle:49", function() refraction_angle(i, 1L, 1L, .py_tup(0L, 0L, 1L)))
})
.py_try({
.py_rec("optics_utils:test_refraction_angle:53", function() refraction_angle(i, 1L, 1L, normal_ray))
})
.py_try({
.py_rec("optics_utils:test_refraction_angle:57", function() refraction_angle(i, 1L, 1L, plane = P))
})
.py_try({
.py_rec("optics_utils:test_refraction_angle:61", function() refraction_angle(r1, 1L, 1L, plane = P))
})
.py_try({
.py_rec("optics_utils:test_refraction_angle:63", function() refraction_angle(r1, m1, .py_flt(1.33), plane = P))
})
.py_try({
.py_rec("optics_utils:test_refraction_angle:65", function() refraction_angle(r1, 1L, m2, plane = P))
})
.py_try({
.py_rec("optics_utils:test_refraction_angle:67", function() refraction_angle(r1, n1, n2, plane = P))
})
.py_try({
.py_rec("optics_utils:test_refraction_angle:69", function() refraction_angle(r1, .py_flt(1.33), 1L, plane = P))
})
.py_try({
.py_rec("optics_utils:test_refraction_angle:70", function() refraction_angle(r1, 1L, 1L, normal_ray))
})
.py_try({
.py_rec("optics_utils:test_refraction_angle:72", function() ae(refraction_angle(.py_flt(0.5), 1L, 2L), .py_flt(0.24207), 5L))
})
.py_try({
.py_rec("optics_utils:test_refraction_angle:73", function() ae(refraction_angle(.py_flt(0.5), 2L, 1L), .py_flt(1.28293), 5L))
})
.py_try({
.py_rec("optics_utils:test_refraction_angle:74", function() refraction_angle(r1, m1, m2, normal_ray, P))
})
.py_try({
.py_rec("optics_utils:test_refraction_angle:75", function() refraction_angle(m1, m1, m2))
})
.py_try({
.py_rec("optics_utils:test_refraction_angle:76", function() refraction_angle(r1, m1, m2, NULL, i))
})
.py_try({
.py_rec("optics_utils:test_refraction_angle:77", function() refraction_angle(r1, m1, m2, m2))
}))---")
  expect_identical(got[["optics_utils:test_refraction_angle:33"]], "Matrix([[1], [1], [-1]])")
  expect_identical(got[["optics_utils:test_refraction_angle:37"]], "Matrix([[1], [1], [-1]])")
  expect_identical(got[["optics_utils:test_refraction_angle:41"]], "Matrix([[1], [1], [-1]])")
  expect_identical(got[["optics_utils:test_refraction_angle:45"]], "Matrix([[1], [1], [-1]])")
  expect_identical(got[["optics_utils:test_refraction_angle:49"]], "Matrix([[1], [1], [-1]])")
  expect_identical(got[["optics_utils:test_refraction_angle:53"]], "Matrix([[1], [1], [-1]])")
  expect_identical(got[["optics_utils:test_refraction_angle:57"]], "Matrix([[1], [1], [-1]])")
  expect_identical(got[["optics_utils:test_refraction_angle:61"]], "Ray3D(Point3D(0, 0, 0), Point3D(1, 1, -1))")
  expect_identical(got[["optics_utils:test_refraction_angle:63"]], "Ray3D(Point3D(0, 0, 0), Point3D(100/133, 100/133, -789378201649271*sqrt(3)/1000000000000000))")
  expect_identical(got[["optics_utils:test_refraction_angle:65"]], "Ray3D(Point3D(0, 0, 0), Point3D(1, 1, -1))")
  expect_identical(got[["optics_utils:test_refraction_angle:67"]], "Ray3D(Point3D(0, 0, 0), Point3D(n1/n2, n1/n2, -sqrt(3)*sqrt(-2*n1^2/(3*n2^2) + 1)))")
  expect_identical(got[["optics_utils:test_refraction_angle:69"]], "0")
  expect_identical(got[["optics_utils:test_refraction_angle:70"]], "Ray3D(Point3D(0, 0, 0), Point3D(1, 1, -1))")
  expect_identical(got[["optics_utils:test_refraction_angle:72"]], "True")
  expect_identical(got[["optics_utils:test_refraction_angle:73"]], "True")
  expect_identical(got[["optics_utils:test_refraction_angle:74"]], "ERR ValueError: Either plane or normal is acceptable.")
  expect_identical(got[["optics_utils:test_refraction_angle:75"]], "ERR TypeError: incident should be a Matrix, Ray3D, or sequence")
  expect_identical(got[["optics_utils:test_refraction_angle:76"]], "ERR TypeError: plane should be an instance of geometry.plane.Plane")
  expect_identical(got[["optics_utils:test_refraction_angle:77"]], "ERR TypeError: Normal should be a Matrix, Ray3D, or sequence")
})

test_that("optics_utils test_fresnel_coefficients matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
ae <- function(a, b, n) comp(a, b, .py_pow(10L, .py_neg(n)))
})
.py_try({
.py_rec("optics_utils:test_fresnel_coefficients:81", function() .py_all(.py_comp(.py_zip(fresnel_coefficients(.py_flt(0.5), 1L, .py_flt(1.33)), list(.py_flt(0.11163), .py_flt(-0.17138), .py_flt(0.83581), .py_flt(0.82862))), function(.it) { i <- .py_getitem(.it, 0L); j <- .py_getitem(.it, 1L); list(ae(i, j, 5L)) })))
})
.py_try({
.py_rec("optics_utils:test_fresnel_coefficients:84", function() .py_all(.py_comp(.py_zip(fresnel_coefficients(.py_flt(0.5), .py_flt(1.33), 1L), list(.py_flt(-0.07726), .py_flt(0.20482), .py_flt(1.22724), .py_flt(1.20482))), function(.it) { i <- .py_getitem(.it, 0L); j <- .py_getitem(.it, 1L); list(ae(i, j, 5L)) })))
})
.py_try({
m1 <- Medium("m1")
})
.py_try({
m2 <- Medium("m2", n = 2L)
})
.py_try({
.py_rec("optics_utils:test_fresnel_coefficients:89", function() .py_all(.py_comp(.py_zip(fresnel_coefficients(.py_flt(0.3), m1, m2), list(.py_flt(0.31784), .py_flt(-0.34865), .py_flt(0.65892), .py_flt(0.65135))), function(.it) { i <- .py_getitem(.it, 0L); j <- .py_getitem(.it, 1L); list(ae(i, j, 5L)) })))
})
.py_try({
ans <- list(list(.py_flt(-0.23563), .py_flt(-0.97184)), list(.py_flt(0.81648), .py_flt(-0.57738)))
})
.py_try({
got <- fresnel_coefficients(.py_flt(0.6), m2, m1)
})
.py_try({
for (.it_ in .py_iter(.py_cap(.py_zip(got, ans), 60L))) {
  i <- .py_getitem(.it_, 0L); j <- .py_getitem(.it_, 1L)
  .py_try({
  for (.it_ in .py_iter(.py_cap(.py_zip(.py_method(i, "as_real_imag"), j), 60L))) {
    a <- .py_getitem(.it_, 0L); b <- .py_getitem(.it_, 1L)
    .py_try({
    .py_rec("optics_utils:test_fresnel_coefficients:96", function() ae(a, b, 5L))
    })
  }
  })
}
}))---")
  expect_identical(got[["optics_utils:test_fresnel_coefficients:81"]], "True")
  expect_identical(got[["optics_utils:test_fresnel_coefficients:84"]], "True")
  expect_identical(got[["optics_utils:test_fresnel_coefficients:89"]], "True")
  expect_identical(got[["optics_utils:test_fresnel_coefficients:96"]], "True")
  expect_identical(got[["optics_utils:test_fresnel_coefficients:96#2"]], "True")
  expect_identical(got[["optics_utils:test_fresnel_coefficients:96#3"]], "True")
  expect_identical(got[["optics_utils:test_fresnel_coefficients:96#4"]], "True")
})

test_that("optics_utils test_deviation matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
ae <- function(a, b, n) comp(a, b, .py_pow(10L, .py_neg(n)))
})
.py_try({
.tmpu_ <- .py_symbols("n1, n2"); n1 <- .py_getitem(.tmpu_, 0L); n2 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
r1 <- Ray3D(Point3D(-1L, -1L, 1L), Point3D(0L, 0L, 0L))
})
.py_try({
n <- .py_Matrix(list(0L, 0L, 1L))
})
.py_try({
i <- .py_Matrix(list(-1L, -1L, -1L))
})
.py_try({
normal_ray <- Ray3D(Point3D(0L, 0L, 0L), Point3D(0L, 0L, 1L))
})
.py_try({
P <- Plane(Point3D(0L, 0L, 0L), normal_vector = list(0L, 0L, 1L))
})
.py_try({
.py_rec("optics_utils:test_deviation:106", function() deviation(r1, 1L, 1L, normal = n))
})
.py_try({
.py_rec("optics_utils:test_deviation:107", function() deviation(r1, 1L, 1L, plane = P))
})
.py_try({
.py_rec("optics_utils:test_deviation:108", function() .py_lt(.py_add(.py_method(deviation(r1, 1L, .py_flt(1.1), plane = P), "evalf", 3L), .py_flt(0.119)), .py_flt(0.001)))
})
.py_try({
.py_rec("optics_utils:test_deviation:109", function() .py_lt(.py_add(.py_method(deviation(i, 1L, .py_flt(1.1), normal = normal_ray), "evalf", 3L), .py_flt(0.119)), .py_flt(0.001)))
})
.py_try({
.py_rec("optics_utils:test_deviation:110", function() deviation(r1, .py_flt(1.33), 1L, plane = P))
})
.py_try({
.py_rec("optics_utils:test_deviation:111", function() deviation(r1, 1L, 1L, normal = list(0L, 0L, 1L)))
})
.py_try({
.py_rec("optics_utils:test_deviation:112", function() deviation(list(-1L, -1L, -1L), 1L, 1L, normal = list(0L, 0L, 1L)))
})
.py_try({
.py_rec("optics_utils:test_deviation:113", function() ae(deviation(.py_flt(0.5), 1L, 2L), .py_flt(-0.25793), 5L))
})
.py_try({
.py_rec("optics_utils:test_deviation:114", function() ae(deviation(.py_flt(0.5), 2L, 1L), .py_flt(0.78293), 5L))
}))---")
  expect_identical(got[["optics_utils:test_deviation:106"]], "0")
  expect_identical(got[["optics_utils:test_deviation:107"]], "0")
  expect_identical(got[["optics_utils:test_deviation:108"]], "True")
  expect_identical(got[["optics_utils:test_deviation:109"]], "True")
  expect_identical(got[["optics_utils:test_deviation:110"]], "None")
  expect_identical(got[["optics_utils:test_deviation:111"]], "0")
  expect_identical(got[["optics_utils:test_deviation:112"]], "0")
  expect_identical(got[["optics_utils:test_deviation:113"]], "True")
  expect_identical(got[["optics_utils:test_deviation:114"]], "True")
})

test_that("optics_utils test_brewster_angle matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
ae <- function(a, b, n) comp(a, b, .py_pow(10L, .py_neg(n)))
})
.py_try({
m1 <- Medium("m1", n = 1L)
})
.py_try({
m2 <- Medium("m2", n = .py_flt(1.33))
})
.py_try({
.py_rec("optics_utils:test_brewster_angle:120", function() ae(brewster_angle(m1, m2), .py_flt(0.93), 2L))
})
.py_try({
m1 <- Medium("m1", permittivity = .py_unit("e0"), n = 1L)
})
.py_try({
m2 <- Medium("m2", permittivity = .py_unit("e0"), n = .py_flt(1.33))
})
.py_try({
.py_rec("optics_utils:test_brewster_angle:123", function() ae(brewster_angle(m1, m2), .py_flt(0.93), 2L))
})
.py_try({
.py_rec("optics_utils:test_brewster_angle:124", function() ae(brewster_angle(1L, .py_flt(1.33)), .py_flt(0.93), 2L))
}))---")
  expect_identical(got[["optics_utils:test_brewster_angle:120"]], "True")
  expect_identical(got[["optics_utils:test_brewster_angle:123"]], "True")
  expect_identical(got[["optics_utils:test_brewster_angle:124"]], "True")
})

test_that("optics_utils test_critical_angle matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
ae <- function(a, b, n) comp(a, b, .py_pow(10L, .py_neg(n)))
})
.py_try({
m1 <- Medium("m1", n = 1L)
})
.py_try({
m2 <- Medium("m2", n = .py_flt(1.33))
})
.py_try({
.py_rec("optics_utils:test_critical_angle:130", function() ae(critical_angle(m2, m1), .py_flt(0.85), 2L))
}))---")
  expect_identical(got[["optics_utils:test_critical_angle:130"]], "True")
})

test_that("optics_utils test_lens_makers_formula matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
ae <- function(a, b, n) comp(a, b, .py_pow(10L, .py_neg(n)))
})
.py_try({
.tmpu_ <- .py_symbols("n1, n2"); n1 <- .py_getitem(.tmpu_, 0L); n2 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
m1 <- Medium("m1", permittivity = .py_unit("e0"), n = 1L)
})
.py_try({
m2 <- Medium("m2", permittivity = .py_unit("e0"), n = .py_flt(1.33))
})
.py_try({
.py_rec("optics_utils:test_lens_makers_formula:137", function() lens_makers_formula(n1, n2, 10L, -10L))
})
.py_try({
.py_rec("optics_utils:test_lens_makers_formula:138", function() ae(lens_makers_formula(m1, m2, 10L, -10L), .py_flt(-20.15), 2L))
})
.py_try({
.py_rec("optics_utils:test_lens_makers_formula:139", function() ae(lens_makers_formula(.py_flt(1.33), 1L, 10L, -10L), .py_flt(15.15), 2L))
}))---")
  expect_identical(got[["optics_utils:test_lens_makers_formula:137"]], "5.0*n2/(n1 - n2)")
  expect_identical(got[["optics_utils:test_lens_makers_formula:138"]], "True")
  expect_identical(got[["optics_utils:test_lens_makers_formula:139"]], "True")
})

test_that("optics_utils test_mirror_formula matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
ae <- function(a, b, n) comp(a, b, .py_pow(10L, .py_neg(n)))
})
.py_try({
.tmpu_ <- .py_symbols("u, v, f"); u <- .py_getitem(.tmpu_, 0L); v <- .py_getitem(.tmpu_, 1L); f <- .py_getitem(.tmpu_, 2L)
})
.py_try({
.py_rec("optics_utils:test_mirror_formula:144", function() mirror_formula(focal_length = f, u = u))
})
.py_try({
.py_rec("optics_utils:test_mirror_formula:145", function() mirror_formula(focal_length = f, v = v))
})
.py_try({
.py_rec("optics_utils:test_mirror_formula:146", function() mirror_formula(u = u, v = v))
})
.py_try({
.py_rec("optics_utils:test_mirror_formula:147", function() mirror_formula(u = oo, v = v))
})
.py_try({
.py_rec("optics_utils:test_mirror_formula:148", function() mirror_formula(u = oo, v = oo))
})
.py_try({
.py_rec("optics_utils:test_mirror_formula:149", function() mirror_formula(focal_length = oo, u = u))
})
.py_try({
.py_rec("optics_utils:test_mirror_formula:150", function() mirror_formula(u = u, v = oo))
})
.py_try({
.py_rec("optics_utils:test_mirror_formula:151", function() mirror_formula(focal_length = oo, v = oo))
})
.py_try({
.py_rec("optics_utils:test_mirror_formula:152", function() mirror_formula(focal_length = f, v = oo))
})
.py_try({
.py_rec("optics_utils:test_mirror_formula:153", function() mirror_formula(focal_length = oo, v = v))
})
.py_try({
.py_rec("optics_utils:test_mirror_formula:154", function() mirror_formula(focal_length = oo, u = oo))
})
.py_try({
.py_rec("optics_utils:test_mirror_formula:155", function() mirror_formula(focal_length = f, u = oo))
})
.py_try({
.py_rec("optics_utils:test_mirror_formula:156", function() mirror_formula(focal_length = oo, u = u))
})
.py_try({
.py_rec("optics_utils:test_mirror_formula:157", function() mirror_formula(focal_length = f, u = u, v = v))
}))---")
  expect_identical(got[["optics_utils:test_mirror_formula:144"]], "f*u/(-f + u)")
  expect_identical(got[["optics_utils:test_mirror_formula:145"]], "f*v/(-f + v)")
  expect_identical(got[["optics_utils:test_mirror_formula:146"]], "u*v/(u + v)")
  expect_identical(got[["optics_utils:test_mirror_formula:147"]], "v")
  expect_identical(got[["optics_utils:test_mirror_formula:148"]], "oo")
  expect_identical(got[["optics_utils:test_mirror_formula:149"]], "-u")
  expect_identical(got[["optics_utils:test_mirror_formula:150"]], "u")
  expect_identical(got[["optics_utils:test_mirror_formula:151"]], "oo")
  expect_identical(got[["optics_utils:test_mirror_formula:152"]], "f")
  expect_identical(got[["optics_utils:test_mirror_formula:153"]], "-v")
  expect_identical(got[["optics_utils:test_mirror_formula:154"]], "oo")
  expect_identical(got[["optics_utils:test_mirror_formula:155"]], "f")
  expect_identical(got[["optics_utils:test_mirror_formula:156"]], "-u")
  expect_identical(got[["optics_utils:test_mirror_formula:157"]], "ERR ValueError: Please provide only two parameters")
})

test_that("optics_utils test_lens_formula matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
ae <- function(a, b, n) comp(a, b, .py_pow(10L, .py_neg(n)))
})
.py_try({
.tmpu_ <- .py_symbols("u, v, f"); u <- .py_getitem(.tmpu_, 0L); v <- .py_getitem(.tmpu_, 1L); f <- .py_getitem(.tmpu_, 2L)
})
.py_try({
.py_rec("optics_utils:test_lens_formula:162", function() lens_formula(focal_length = f, u = u))
})
.py_try({
.py_rec("optics_utils:test_lens_formula:163", function() lens_formula(focal_length = f, v = v))
})
.py_try({
.py_rec("optics_utils:test_lens_formula:164", function() lens_formula(u = u, v = v))
})
.py_try({
.py_rec("optics_utils:test_lens_formula:165", function() lens_formula(u = oo, v = v))
})
.py_try({
.py_rec("optics_utils:test_lens_formula:166", function() lens_formula(u = oo, v = oo))
})
.py_try({
.py_rec("optics_utils:test_lens_formula:167", function() lens_formula(focal_length = oo, u = u))
})
.py_try({
.py_rec("optics_utils:test_lens_formula:168", function() lens_formula(u = u, v = oo))
})
.py_try({
.py_rec("optics_utils:test_lens_formula:169", function() lens_formula(focal_length = oo, v = oo))
})
.py_try({
.py_rec("optics_utils:test_lens_formula:170", function() lens_formula(focal_length = oo, v = v))
})
.py_try({
.py_rec("optics_utils:test_lens_formula:171", function() lens_formula(focal_length = f, v = oo))
})
.py_try({
.py_rec("optics_utils:test_lens_formula:172", function() lens_formula(focal_length = oo, u = oo))
})
.py_try({
.py_rec("optics_utils:test_lens_formula:173", function() lens_formula(focal_length = oo, u = u))
})
.py_try({
.py_rec("optics_utils:test_lens_formula:174", function() lens_formula(focal_length = f, u = oo))
})
.py_try({
.py_rec("optics_utils:test_lens_formula:175", function() lens_formula(focal_length = f, u = u, v = v))
}))---")
  expect_identical(got[["optics_utils:test_lens_formula:162"]], "f*u/(f + u)")
  expect_identical(got[["optics_utils:test_lens_formula:163"]], "f*v/(f - v)")
  expect_identical(got[["optics_utils:test_lens_formula:164"]], "u*v/(u - v)")
  expect_identical(got[["optics_utils:test_lens_formula:165"]], "v")
  expect_identical(got[["optics_utils:test_lens_formula:166"]], "oo")
  expect_identical(got[["optics_utils:test_lens_formula:167"]], "u")
  expect_identical(got[["optics_utils:test_lens_formula:168"]], "-u")
  expect_identical(got[["optics_utils:test_lens_formula:169"]], "-oo")
  expect_identical(got[["optics_utils:test_lens_formula:170"]], "v")
  expect_identical(got[["optics_utils:test_lens_formula:171"]], "-f")
  expect_identical(got[["optics_utils:test_lens_formula:172"]], "oo")
  expect_identical(got[["optics_utils:test_lens_formula:173"]], "u")
  expect_identical(got[["optics_utils:test_lens_formula:174"]], "f")
  expect_identical(got[["optics_utils:test_lens_formula:175"]], "ERR ValueError: Please provide only two parameters")
})

test_that("optics_utils test_hyperfocal_distance matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
ae <- function(a, b, n) comp(a, b, .py_pow(10L, .py_neg(n)))
})
.py_try({
.tmpu_ <- .py_symbols("f, N, c"); f <- .py_getitem(.tmpu_, 0L); N <- .py_getitem(.tmpu_, 1L); c <- .py_getitem(.tmpu_, 2L)
})
.py_try({
.py_rec("optics_utils:test_hyperfocal_distance:180", function() hyperfocal_distance(f = f, N = N, c = c))
})
.py_try({
.py_rec("optics_utils:test_hyperfocal_distance:181", function() ae(hyperfocal_distance(f = .py_flt(0.5), N = 8L, c = .py_flt(0.0033)), .py_flt(9.47), 2L))
}))---")
  expect_identical(got[["optics_utils:test_hyperfocal_distance:180"]], "f^2/(N*c)")
  expect_identical(got[["optics_utils:test_hyperfocal_distance:181"]], "True")
})

test_that("optics_utils test_transverse_magnification matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
ae <- function(a, b, n) comp(a, b, .py_pow(10L, .py_neg(n)))
})
.py_try({
.tmpu_ <- .py_symbols("si, so"); si <- .py_getitem(.tmpu_, 0L); so <- .py_getitem(.tmpu_, 1L)
})
.py_try({
.py_rec("optics_utils:test_transverse_magnification:186", function() transverse_magnification(si, so))
})
.py_try({
.py_rec("optics_utils:test_transverse_magnification:187", function() transverse_magnification(30L, 15L))
}))---")
  expect_identical(got[["optics_utils:test_transverse_magnification:186"]], "-si/so")
  expect_identical(got[["optics_utils:test_transverse_magnification:187"]], "-2")
})

test_that("optics_utils test_lens_makers_formula_thick_lens matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
ae <- function(a, b, n) comp(a, b, .py_pow(10L, .py_neg(n)))
})
.py_try({
.tmpu_ <- .py_symbols("n1, n2"); n1 <- .py_getitem(.tmpu_, 0L); n2 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
m1 <- Medium("m1", permittivity = .py_unit("e0"), n = 1L)
})
.py_try({
m2 <- Medium("m2", permittivity = .py_unit("e0"), n = .py_flt(1.33))
})
.py_try({
.py_rec("optics_utils:test_lens_makers_formula_thick_lens:194", function() ae(lens_makers_formula(m1, m2, 10L, -10L, d = 1L), .py_flt(-19.82), 2L))
})
.py_try({
.py_rec("optics_utils:test_lens_makers_formula_thick_lens:195", function() lens_makers_formula(n1, n2, 1L, -1L, d = .py_flt(0.1)))
}))---")
  expect_identical(got[["optics_utils:test_lens_makers_formula_thick_lens:194"]], "True")
  expect_identical(got[["optics_utils:test_lens_makers_formula_thick_lens:195"]], "n2/((2.0 - (0.1*n1 - 0.1*n2)/n1)*(n1 - n2))")
})

test_that("optics_utils test_lens_makers_formula_plano_lens matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
ae <- function(a, b, n) comp(a, b, .py_pow(10L, .py_neg(n)))
})
.py_try({
.tmpu_ <- .py_symbols("n1, n2"); n1 <- .py_getitem(.tmpu_, 0L); n2 <- .py_getitem(.tmpu_, 1L)
})
.py_try({
m1 <- Medium("m1", permittivity = .py_unit("e0"), n = 1L)
})
.py_try({
m2 <- Medium("m2", permittivity = .py_unit("e0"), n = .py_flt(1.33))
})
.py_try({
.py_rec("optics_utils:test_lens_makers_formula_plano_lens:202", function() ae(lens_makers_formula(m1, m2, 10L, oo), .py_flt(-40.3), 2L))
})
.py_try({
.py_rec("optics_utils:test_lens_makers_formula_plano_lens:203", function() lens_makers_formula(n1, n2, 10L, oo))
}))---")
  expect_identical(got[["optics_utils:test_lens_makers_formula_plano_lens:202"]], "True")
  expect_identical(got[["optics_utils:test_lens_makers_formula_plano_lens:203"]], "10.0*n2/(n1 - n2)")
})

test_that("waves test_twave matches SymPy", {
  got <- .physmisc_battery_run(r"---(epsilon <- Symbol("epsilon")
mu <- Symbol("mu")
.py_try({
c <- .py_method(.py_unit("speed_of_light"), "convert_to", .py_div(.py_unit("m"), .py_unit("s")))
})
.py_try({
.tmpu_ <- .py_symbols("A1, phi1, A2, phi2, f"); A1 <- .py_getitem(.tmpu_, 0L); phi1 <- .py_getitem(.tmpu_, 1L); A2 <- .py_getitem(.tmpu_, 2L); phi2 <- .py_getitem(.tmpu_, 3L); f <- .py_getitem(.tmpu_, 4L)
})
.py_try({
n <- Symbol("n")
})
.py_try({
t <- Symbol("t")
})
.py_try({
x <- Symbol("x")
})
.py_try({
E <- Function("E")
})
.py_try({
w1 <- TWave(A1, f, phi1)
})
.py_try({
w2 <- TWave(A2, f, phi2)
})
.py_try({
.py_rec("waves:test_twave:25", function() .py_attr(w1, "amplitude"))
})
.py_try({
.py_rec("waves:test_twave:26", function() .py_attr(w1, "frequency"))
})
.py_try({
.py_rec("waves:test_twave:27", function() .py_attr(w1, "phase"))
})
.py_try({
.py_rec("waves:test_twave:28", function() .py_attr(w1, "wavelength"))
})
.py_try({
.py_rec("waves:test_twave:29", function() .py_attr(w1, "time_period"))
})
.py_try({
.py_rec("waves:test_twave:30", function() .py_attr(w1, "angular_velocity"))
})
.py_try({
.py_rec("waves:test_twave:31", function() .py_attr(w1, "wavenumber"))
})
.py_try({
.py_rec("waves:test_twave:32", function() .py_attr(w1, "speed"))
})
.py_try({
w3 <- .py_add(w1, w2)
})
.py_try({
.py_rec("waves:test_twave:35", function() .py_attr(w3, "amplitude"))
})
.py_try({
.py_rec("waves:test_twave:36", function() .py_attr(w3, "frequency"))
})
.py_try({
.py_rec("waves:test_twave:37", function() .py_attr(w3, "phase"))
})
.py_try({
.py_rec("waves:test_twave:38", function() .py_attr(w3, "wavelength"))
})
.py_try({
.py_rec("waves:test_twave:39", function() .py_attr(w3, "time_period"))
})
.py_try({
.py_rec("waves:test_twave:40", function() .py_attr(w3, "angular_velocity"))
})
.py_try({
.py_rec("waves:test_twave:41", function() .py_attr(w3, "wavenumber"))
})
.py_try({
.py_rec("waves:test_twave:42", function() .py_attr(w3, "speed"))
})
.py_try({
.py_rec("waves:test_twave:44", function() .py_method(w3, "rewrite", "pde"))
})
.py_try({
.py_rec("waves:test_twave:45", function() .py_method(w3, "rewrite", .py_cos))
})
.py_try({
.py_rec("waves:test_twave:48", function() .py_method(w3, "rewrite", .py_exp))
})
.py_try({
w4 <- TWave(A1, NULL, 0L, .py_div(1L, f))
})
.py_try({
.py_rec("waves:test_twave:53", function() .py_attr(w4, "frequency"))
})
.py_try({
w5 <- .py_sub(w1, w2)
})
.py_try({
.py_rec("waves:test_twave:56", function() .py_attr(w5, "amplitude"))
})
.py_try({
.py_rec("waves:test_twave:57", function() .py_attr(w5, "frequency"))
})
.py_try({
.py_rec("waves:test_twave:58", function() .py_attr(w5, "phase"))
})
.py_try({
.py_rec("waves:test_twave:59", function() .py_attr(w5, "wavelength"))
})
.py_try({
.py_rec("waves:test_twave:60", function() .py_attr(w5, "time_period"))
})
.py_try({
.py_rec("waves:test_twave:61", function() .py_attr(w5, "angular_velocity"))
})
.py_try({
.py_rec("waves:test_twave:62", function() .py_attr(w5, "wavenumber"))
})
.py_try({
.py_rec("waves:test_twave:63", function() .py_attr(w5, "speed"))
})
.py_try({
.py_rec("waves:test_twave:65", function() .py_method(w5, "rewrite", "pde"))
})
.py_try({
.py_rec("waves:test_twave:66", function() .py_method(w5, "rewrite", .py_cos))
})
.py_try({
.py_rec("waves:test_twave:69", function() .py_method(w5, "rewrite", .py_exp))
})
.py_try({
w6 <- .py_mul(2L, w1)
})
.py_try({
.py_rec("waves:test_twave:74", function() .py_attr(w6, "amplitude"))
})
.py_try({
.py_rec("waves:test_twave:75", function() .py_attr(w6, "frequency"))
})
.py_try({
.py_rec("waves:test_twave:76", function() .py_attr(w6, "phase"))
})
.py_try({
w7 <- .py_neg(w6)
})
.py_try({
.py_rec("waves:test_twave:78", function() .py_attr(w7, "amplitude"))
})
.py_try({
.py_rec("waves:test_twave:79", function() .py_attr(w7, "frequency"))
})
.py_try({
.py_rec("waves:test_twave:80", function() .py_attr(w7, "phase"))
})
.py_try({
.py_rec("waves:test_twave:82", function() TWave(A1))
})
.py_try({
.py_rec("waves:test_twave:83", function() TWave(A1, f, phi1, t))
}))---")
  expect_identical(got[["waves:test_twave:25"]], "A1")
  expect_identical(got[["waves:test_twave:26"]], "f")
  expect_identical(got[["waves:test_twave:27"]], "phi1")
  expect_identical(got[["waves:test_twave:28"]], "299792458*meter/(second*f*n)")
  expect_identical(got[["waves:test_twave:29"]], "1/f")
  expect_identical(got[["waves:test_twave:30"]], "2*pi*f")
  expect_identical(got[["waves:test_twave:31"]], "pi*second*f*n/(149896229*meter)")
  expect_identical(got[["waves:test_twave:32"]], "299792458*meter/(second*n)")
  expect_identical(got[["waves:test_twave:35"]], "sqrt(A1^2 + 2*A1*A2*cos(phi1 - phi2) + A2^2)")
  expect_identical(got[["waves:test_twave:36"]], "f")
  expect_identical(got[["waves:test_twave:37"]], "atan2(A1*sin(phi1) + A2*sin(phi2), A1*cos(phi1) + A2*cos(phi2))")
  expect_identical(got[["waves:test_twave:38"]], "299792458*meter/(second*f*n)")
  expect_identical(got[["waves:test_twave:39"]], "1/f")
  expect_identical(got[["waves:test_twave:40"]], "2*pi*f")
  expect_identical(got[["waves:test_twave:41"]], "pi*second*f*n/(149896229*meter)")
  expect_identical(got[["waves:test_twave:42"]], "299792458*meter/(second*n)")
  expect_identical(got[["waves:test_twave:44"]], "epsilon*mu*Derivative(E(x, t), (t, 2)) + Derivative(E(x, t), (x, 2))")
  expect_identical(got[["waves:test_twave:45"]], "sqrt(A1^2 + 2*A1*A2*cos(phi1 - phi2) + A2^2)*cos(-2*pi*f*t + atan2(A1*sin(phi1) + A2*sin(phi2), A1*cos(phi1) + A2*cos(phi2)) + pi*second*f*n*x/(149896229*meter))")
  expect_identical(got[["waves:test_twave:48"]], "sqrt(A1^2 + 2*A1*A2*cos(phi1 - phi2) + A2^2)*exp(I*(-2*pi*f*t + atan2(A1*sin(phi1) + A2*sin(phi2), A1*cos(phi1) + A2*cos(phi2)) + pi*second*f*n*x/(149896229*meter)))")
  expect_identical(got[["waves:test_twave:53"]], "f")
  expect_identical(got[["waves:test_twave:56"]], "sqrt(A1^2 - 2*A1*A2*cos(phi1 - phi2) + A2^2)")
  expect_identical(got[["waves:test_twave:57"]], "f")
  expect_identical(got[["waves:test_twave:58"]], "atan2(A1*sin(phi1) - A2*sin(phi2), A1*cos(phi1) - A2*cos(phi2))")
  expect_identical(got[["waves:test_twave:59"]], "299792458*meter/(second*f*n)")
  expect_identical(got[["waves:test_twave:60"]], "1/f")
  expect_identical(got[["waves:test_twave:61"]], "2*pi*f")
  expect_identical(got[["waves:test_twave:62"]], "pi*second*f*n/(149896229*meter)")
  expect_identical(got[["waves:test_twave:63"]], "299792458*meter/(second*n)")
  expect_identical(got[["waves:test_twave:65"]], "epsilon*mu*Derivative(E(x, t), (t, 2)) + Derivative(E(x, t), (x, 2))")
  expect_identical(got[["waves:test_twave:66"]], "sqrt(A1^2 - 2*A1*A2*cos(phi1 - phi2) + A2^2)*cos(-2*pi*f*t + atan2(A1*sin(phi1) - A2*sin(phi2), A1*cos(phi1) - A2*cos(phi2)) + pi*second*f*n*x/(149896229*meter))")
  expect_identical(got[["waves:test_twave:69"]], "sqrt(A1^2 - 2*A1*A2*cos(phi1 - phi2) + A2^2)*exp(I*(-2*pi*f*t + atan2(A1*sin(phi1) - A2*sin(phi2), A1*cos(phi1) - A2*cos(phi2)) + pi*second*f*n*x/(149896229*meter)))")
  expect_identical(got[["waves:test_twave:74"]], "2*A1")
  expect_identical(got[["waves:test_twave:75"]], "f")
  expect_identical(got[["waves:test_twave:76"]], "phi1")
  expect_identical(got[["waves:test_twave:78"]], "-2*A1")
  expect_identical(got[["waves:test_twave:79"]], "f")
  expect_identical(got[["waves:test_twave:80"]], "phi1")
  expect_identical(got[["waves:test_twave:82"]], "ERR ValueError: Either frequency or time period is needed.")
  expect_identical(got[["waves:test_twave:83"]], "ERR ValueError: frequency and time_period should be consistent.")
})
