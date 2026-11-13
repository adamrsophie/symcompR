skip_unless_full_tests()

# Expected values are the outputs of the same calls in SymPy (hydrogen, wigner, sho, qho_1d, pring and physics.matrices);
# the cases come from its test files and run through the helpers in
# helper-physics-misc-battery.R.

test_that("hydrogen test_wavefunction matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("n r Z"); n <- .py_getitem(.tmpu_, 0L); r <- .py_getitem(.tmpu_, 1L); Z <- .py_getitem(.tmpu_, 2L)
})
feq <- function(a, b, max_relative_error = .py_flt(1e-12), max_absolute_error = .py_flt(1e-12)) {
  a <- .py_float(a)
  b <- .py_float(b)
  if (.py_truthy(.py_lt(.py_abs(.py_sub(a, b)), max_absolute_error))) {
    return(TRUE)
  }
  if (.py_truthy(.py_gt(.py_abs(b), .py_abs(a)))) {
    relative_error <- .py_abs(.py_div(.py_sub(a, b), b))
  } else {
    relative_error <- .py_abs(.py_div(.py_sub(a, b), a))
  }
  return(.py_le(relative_error, max_relative_error))
}
.py_try({
a <- .py_div(1L, Z)
})
.py_try({
R <- .py_dict(list(.py_tup(1L, 0L), .py_tup(2L, 0L), .py_tup(2L, 1L), .py_tup(3L, 0L), .py_tup(3L, 1L), .py_tup(3L, 2L), .py_tup(4L, 0L), .py_tup(4L, 1L), .py_tup(4L, 2L), .py_tup(4L, 3L)), list(.py_mul(.py_mul(2L, .py_sqrt(.py_div(1L, .py_pow(a, 3L)))), .py_exp(.py_div(.py_neg(r), a))), .py_mul(.py_mul(.py_sqrt(.py_div(1L, .py_mul(2L, .py_pow(a, 3L)))), .py_exp(.py_div(.py_neg(r), .py_mul(2L, a)))), .py_sub(1L, .py_div(r, .py_mul(2L, a)))), .py_div(.py_mul(.py_mul(.py_mul(S_Half, .py_sqrt(.py_div(1L, .py_mul(6L, .py_pow(a, 3L))))), .py_exp(.py_div(.py_neg(r), .py_mul(2L, a)))), r), a), .py_mul(.py_mul(.py_mul(Rational(2L, 3L), .py_sqrt(.py_div(1L, .py_mul(3L, .py_pow(a, 3L))))), .py_exp(.py_div(.py_neg(r), .py_mul(3L, a)))), .py_add(.py_sub(1L, .py_div(.py_mul(2L, r), .py_mul(3L, a))), .py_mul(Rational(2L, 27L), .py_pow(.py_div(r, a), 2L)))), .py_div(.py_mul(.py_mul(.py_mul(.py_mul(Rational(4L, 27L), .py_sqrt(.py_div(2L, .py_mul(3L, .py_pow(a, 3L))))), .py_exp(.py_div(.py_neg(r), .py_mul(3L, a)))), .py_sub(1L, .py_div(r, .py_mul(6L, a)))), r), a), .py_mul(.py_mul(.py_mul(Rational(2L, 81L), .py_sqrt(.py_div(2L, .py_mul(15L, .py_pow(a, 3L))))), .py_exp(.py_div(.py_neg(r), .py_mul(3L, a)))), .py_pow(.py_div(r, a), 2L)), .py_mul(.py_mul(.py_mul(Rational(1L, 4L), .py_sqrt(.py_div(1L, .py_pow(a, 3L)))), .py_exp(.py_div(.py_neg(r), .py_mul(4L, a)))), .py_sub(.py_add(.py_sub(1L, .py_div(.py_mul(3L, r), .py_mul(4L, a))), .py_mul(Rational(1L, 8L), .py_pow(.py_div(r, a), 2L))), .py_mul(Rational(1L, 192L), .py_pow(.py_div(r, a), 3L)))), .py_mul(.py_mul(.py_mul(.py_mul(Rational(1L, 16L), .py_sqrt(.py_div(5L, .py_mul(3L, .py_pow(a, 3L))))), .py_exp(.py_div(.py_neg(r), .py_mul(4L, a)))), .py_add(.py_sub(1L, .py_div(r, .py_mul(4L, a))), .py_mul(Rational(1L, 80L), .py_pow(.py_div(r, a), 2L)))), .py_div(r, a)), .py_mul(.py_mul(.py_mul(.py_mul(Rational(1L, 64L), .py_sqrt(.py_div(1L, .py_mul(5L, .py_pow(a, 3L))))), .py_exp(.py_div(.py_neg(r), .py_mul(4L, a)))), .py_sub(1L, .py_div(r, .py_mul(12L, a)))), .py_pow(.py_div(r, a), 2L)), .py_mul(.py_mul(.py_mul(Rational(1L, 768L), .py_sqrt(.py_div(1L, .py_mul(35L, .py_pow(a, 3L))))), .py_exp(.py_div(.py_neg(r), .py_mul(4L, a)))), .py_pow(.py_div(r, a), 3L))))
})
.py_try({
for (.it_ in .py_iter(.py_cap(R, 60L))) {
  n <- .py_getitem(.it_, 0L); l <- .py_getitem(.it_, 1L)
  .py_try({
  .py_rec("hydrogen:test_wavefunction:50", function() simplify(.py_sub(hydrogen_R_nl(n, l, r, Z), .py_getitem(R, .py_tup(n, l)))))
  })
}
}))---")
  expect_identical(got[["hydrogen:test_wavefunction:50"]], "0")
  expect_identical(got[["hydrogen:test_wavefunction:50#2"]], "0")
  expect_identical(got[["hydrogen:test_wavefunction:50#3"]], "0")
  expect_identical(got[["hydrogen:test_wavefunction:50#4"]], "0")
  expect_identical(got[["hydrogen:test_wavefunction:50#5"]], "0")
  expect_identical(got[["hydrogen:test_wavefunction:50#6"]], "0")
  expect_identical(got[["hydrogen:test_wavefunction:50#7"]], "0")
  expect_identical(got[["hydrogen:test_wavefunction:50#8"]], "0")
  expect_identical(got[["hydrogen:test_wavefunction:50#9"]], "0")
  expect_identical(got[["hydrogen:test_wavefunction:50#10"]], "0")
})

test_that("hydrogen test_norm matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("n r Z"); n <- .py_getitem(.tmpu_, 0L); r <- .py_getitem(.tmpu_, 1L); Z <- .py_getitem(.tmpu_, 2L)
})
feq <- function(a, b, max_relative_error = .py_flt(1e-12), max_absolute_error = .py_flt(1e-12)) {
  a <- .py_float(a)
  b <- .py_float(b)
  if (.py_truthy(.py_lt(.py_abs(.py_sub(a, b)), max_absolute_error))) {
    return(TRUE)
  }
  if (.py_truthy(.py_gt(.py_abs(b), .py_abs(a)))) {
    relative_error <- .py_abs(.py_div(.py_sub(a, b), b))
  } else {
    relative_error <- .py_abs(.py_div(.py_sub(a, b), a))
  }
  return(.py_le(relative_error, max_relative_error))
}
.py_try({
n_max <- 2L
})
.py_try({
for (.it_ in .py_iter(.py_cap(.py_range(.py_add(n_max, 1L)), 60L))) {
  n <- .it_
  .py_try({
  for (.it_ in .py_iter(.py_cap(.py_range(n), 60L))) {
    l <- .it_
    .py_try({
    .py_rec("hydrogen:test_norm:58", function() sym_integrate(.py_mul(.py_pow(hydrogen_R_nl(n, l, r), 2L), .py_pow(r, 2L)), .py_tup(r, 0L, oo)))
    })
  }
  })
}
}))---")
  expect_identical(got[["hydrogen:test_norm:58"]], "1")
  expect_identical(got[["hydrogen:test_norm:58#2"]], "1")
  expect_identical(got[["hydrogen:test_norm:58#3"]], "1")
})

test_that("hydrogen test_psi_nlm matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("n r Z"); n <- .py_getitem(.tmpu_, 0L); r <- .py_getitem(.tmpu_, 1L); Z <- .py_getitem(.tmpu_, 2L)
})
feq <- function(a, b, max_relative_error = .py_flt(1e-12), max_absolute_error = .py_flt(1e-12)) {
  a <- .py_float(a)
  b <- .py_float(b)
  if (.py_truthy(.py_lt(.py_abs(.py_sub(a, b)), max_absolute_error))) {
    return(TRUE)
  }
  if (.py_truthy(.py_gt(.py_abs(b), .py_abs(a)))) {
    relative_error <- .py_abs(.py_div(.py_sub(a, b), b))
  } else {
    relative_error <- .py_abs(.py_div(.py_sub(a, b), a))
  }
  return(.py_le(relative_error, max_relative_error))
}
.py_try({
r <- .py_S("r")
})
.py_try({
phi <- .py_S("phi")
})
.py_try({
theta <- .py_S("theta")
})
.py_try({
.py_rec("hydrogen:test_psi_nlm:64", function() hydrogen_Psi_nlm(1L, 0L, 0L, r, phi, theta))
})
.py_try({
.py_rec("hydrogen:test_psi_nlm:65", function() hydrogen_Psi_nlm(2L, 1L, -1L, r, phi, theta))
})
.py_try({
.py_rec("hydrogen:test_psi_nlm:67", function() hydrogen_Psi_nlm(3L, 2L, 1L, r, phi, theta, 2L))
}))---")
  expect_identical(got[["hydrogen:test_psi_nlm:64"]], "exp(-r)/sqrt(pi)")
  expect_identical(got[["hydrogen:test_psi_nlm:65"]], "r*exp(-r/2)*exp(-I*phi)*sin(theta)/(8*sqrt(pi))")
  expect_identical(got[["hydrogen:test_psi_nlm:67"]], "-8*sqrt(2)*r^2*exp(-2*r/3)*exp(I*phi)*sin(theta)*cos(theta)/(81*sqrt(pi))")
})

test_that("hydrogen test_hydrogen_energies matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("n r Z"); n <- .py_getitem(.tmpu_, 0L); r <- .py_getitem(.tmpu_, 1L); Z <- .py_getitem(.tmpu_, 2L)
})
feq <- function(a, b, max_relative_error = .py_flt(1e-12), max_absolute_error = .py_flt(1e-12)) {
  a <- .py_float(a)
  b <- .py_float(b)
  if (.py_truthy(.py_lt(.py_abs(.py_sub(a, b)), max_absolute_error))) {
    return(TRUE)
  }
  if (.py_truthy(.py_gt(.py_abs(b), .py_abs(a)))) {
    relative_error <- .py_abs(.py_div(.py_sub(a, b), b))
  } else {
    relative_error <- .py_abs(.py_div(.py_sub(a, b), a))
  }
  return(.py_le(relative_error, max_relative_error))
}
.py_try({
.py_rec("hydrogen:test_hydrogen_energies:72", function() hydrogen_E_nl(n, Z))
})
.py_try({
.py_rec("hydrogen:test_hydrogen_energies:73", function() hydrogen_E_nl(n))
})
.py_try({
.py_rec("hydrogen:test_hydrogen_energies:75", function() hydrogen_E_nl(1L, 47L))
})
.py_try({
.py_rec("hydrogen:test_hydrogen_energies:76", function() hydrogen_E_nl(2L, 47L))
})
.py_try({
.py_rec("hydrogen:test_hydrogen_energies:78", function() hydrogen_E_nl(1L))
})
.py_try({
.py_rec("hydrogen:test_hydrogen_energies:79", function() hydrogen_E_nl(2L))
})
.py_try({
.py_rec("hydrogen:test_hydrogen_energies:80", function() hydrogen_E_nl(3L))
})
.py_try({
.py_rec("hydrogen:test_hydrogen_energies:81", function() hydrogen_E_nl(4L))
})
.py_try({
.py_rec("hydrogen:test_hydrogen_energies:82", function() hydrogen_E_nl(100L))
})
.py_try({
.py_rec("hydrogen:test_hydrogen_energies:84", function() hydrogen_E_nl(0L))
}))---")
  expect_identical(got[["hydrogen:test_hydrogen_energies:72"]], "-Z^2/(2*n^2)")
  expect_identical(got[["hydrogen:test_hydrogen_energies:73"]], "-1/(2*n^2)")
  expect_identical(got[["hydrogen:test_hydrogen_energies:75"]], "-2209/2")
  expect_identical(got[["hydrogen:test_hydrogen_energies:76"]], "-2209/8")
  expect_identical(got[["hydrogen:test_hydrogen_energies:78"]], "-1/2")
  expect_identical(got[["hydrogen:test_hydrogen_energies:79"]], "-1/8")
  expect_identical(got[["hydrogen:test_hydrogen_energies:80"]], "-1/18")
  expect_identical(got[["hydrogen:test_hydrogen_energies:81"]], "-1/32")
  expect_identical(got[["hydrogen:test_hydrogen_energies:82"]], "-1/20000")
  expect_identical(got[["hydrogen:test_hydrogen_energies:84"]], "ERR ValueError: 'n' must be positive integer")
})

test_that("hydrogen test_hydrogen_energies_relat matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("n r Z"); n <- .py_getitem(.tmpu_, 0L); r <- .py_getitem(.tmpu_, 1L); Z <- .py_getitem(.tmpu_, 2L)
})
feq <- function(a, b, max_relative_error = .py_flt(1e-12), max_absolute_error = .py_flt(1e-12)) {
  a <- .py_float(a)
  b <- .py_float(b)
  if (.py_truthy(.py_lt(.py_abs(.py_sub(a, b)), max_absolute_error))) {
    return(TRUE)
  }
  if (.py_truthy(.py_gt(.py_abs(b), .py_abs(a)))) {
    relative_error <- .py_abs(.py_div(.py_sub(a, b), b))
  } else {
    relative_error <- .py_abs(.py_div(.py_sub(a, b), a))
  }
  return(.py_le(relative_error, max_relative_error))
}
.py_try({
.py_rec("hydrogen:test_hydrogen_energies_relat:89", function() hydrogen_E_nl_dirac(2L, 0L, Z = 1L, c = 1L))
})
.py_try({
.py_rec("hydrogen:test_hydrogen_energies_relat:90", function() simplify(.py_sub(hydrogen_E_nl_dirac(2L, 0L, Z = 1L, c = 2L), .py_sub(.py_div(.py_add(.py_mul(8L, .py_sqrt(3L)), 16L), .py_sqrt(.py_add(.py_mul(16L, .py_sqrt(3L)), 32L))), 4L))))
})
.py_try({
.py_rec("hydrogen:test_hydrogen_energies_relat:92", function() simplify(.py_sub(hydrogen_E_nl_dirac(2L, 0L, Z = 1L, c = 3L), .py_sub(.py_div(.py_add(.py_mul(54L, .py_sqrt(2L)), 81L), .py_sqrt(.py_add(.py_mul(108L, .py_sqrt(2L)), 162L))), 9L))))
})
.py_try({
.py_rec("hydrogen:test_hydrogen_energies_relat:97", function() simplify(.py_sub(hydrogen_E_nl_dirac(2L, 0L, Z = 1L, c = 137L), .py_sub(.py_div(.py_add(352275361L, .py_mul(10285412L, .py_sqrt(1173L))), .py_sqrt(.py_add(704550722L, .py_mul(20570824L, .py_sqrt(1173L))))), 18769L))))
})
.py_try({
.py_rec("hydrogen:test_hydrogen_energies_relat:99", function() simplify(.py_sub(hydrogen_E_nl_dirac(2L, 0L, Z = 82L, c = 137L), .py_sub(.py_div(.py_add(352275361L, .py_mul(2571353L, .py_sqrt(12045L))), .py_sqrt(.py_add(704550722L, .py_mul(5142706L, .py_sqrt(12045L))))), 18769L))))
})
.py_try({
for (.it_ in .py_iter(.py_cap(.py_range(1L, 5L), 60L))) {
  n <- .it_
  .py_try({
  for (.it_ in .py_iter(.py_cap(.py_range(n), 60L))) {
    l <- .it_
    .py_try({
    .py_rec("hydrogen:test_hydrogen_energies_relat:106", function() feq(hydrogen_E_nl_dirac(n, l), hydrogen_E_nl(n), .py_flt(1e-05), .py_flt(1e-05)))
    })
    .py_try({
    if (.py_truthy(.py_gt(l, 0L))) {
      .py_try({
      .py_rec("hydrogen:test_hydrogen_energies_relat:108", function() feq(hydrogen_E_nl_dirac(n, l, FALSE), hydrogen_E_nl(n), .py_flt(1e-05), .py_flt(1e-05)))
      })
    }
    })
  }
  })
}
})
.py_try({
Z <- 2L
})
.py_try({
for (.it_ in .py_iter(.py_cap(.py_range(1L, 5L), 60L))) {
  n <- .it_
  .py_try({
  for (.it_ in .py_iter(.py_cap(.py_range(n), 60L))) {
    l <- .it_
    .py_try({
    .py_rec("hydrogen:test_hydrogen_energies_relat:113", function() feq(hydrogen_E_nl_dirac(n, l, Z = Z), hydrogen_E_nl(n, Z), .py_flt(0.0001), .py_flt(0.0001)))
    })
    .py_try({
    if (.py_truthy(.py_gt(l, 0L))) {
      .py_try({
      .py_rec("hydrogen:test_hydrogen_energies_relat:115", function() feq(hydrogen_E_nl_dirac(n, l, FALSE, Z), hydrogen_E_nl(n, Z), .py_flt(0.0001), .py_flt(0.0001)))
      })
    }
    })
  }
  })
}
})
.py_try({
Z <- 3L
})
.py_try({
for (.it_ in .py_iter(.py_cap(.py_range(1L, 5L), 60L))) {
  n <- .it_
  .py_try({
  for (.it_ in .py_iter(.py_cap(.py_range(n), 60L))) {
    l <- .it_
    .py_try({
    .py_rec("hydrogen:test_hydrogen_energies_relat:120", function() feq(hydrogen_E_nl_dirac(n, l, Z = Z), hydrogen_E_nl(n, Z), .py_flt(0.001), .py_flt(0.001)))
    })
    .py_try({
    if (.py_truthy(.py_gt(l, 0L))) {
      .py_try({
      .py_rec("hydrogen:test_hydrogen_energies_relat:122", function() feq(hydrogen_E_nl_dirac(n, l, FALSE, Z), hydrogen_E_nl(n, Z), .py_flt(0.001), .py_flt(0.001)))
      })
    }
    })
  }
  })
}
})
.py_try({
.py_rec("hydrogen:test_hydrogen_energies_relat:125", function() hydrogen_E_nl_dirac(0L, 0L))
})
.py_try({
.py_rec("hydrogen:test_hydrogen_energies_relat:126", function() hydrogen_E_nl_dirac(1L, -1L))
})
.py_try({
.py_rec("hydrogen:test_hydrogen_energies_relat:127", function() hydrogen_E_nl_dirac(1L, 0L, FALSE))
}))---")
  expect_identical(got[["hydrogen:test_hydrogen_energies_relat:89"]], "-1 + sqrt(2)/2")
  expect_identical(got[["hydrogen:test_hydrogen_energies_relat:90"]], "0")
  expect_identical(got[["hydrogen:test_hydrogen_energies_relat:92"]], "0")
  expect_identical(got[["hydrogen:test_hydrogen_energies_relat:97"]], "0")
  expect_identical(got[["hydrogen:test_hydrogen_energies_relat:99"]], "0")
  expect_identical(got[["hydrogen:test_hydrogen_energies_relat:106"]], "True")
  expect_identical(got[["hydrogen:test_hydrogen_energies_relat:106#2"]], "True")
  expect_identical(got[["hydrogen:test_hydrogen_energies_relat:106#3"]], "True")
  expect_identical(got[["hydrogen:test_hydrogen_energies_relat:108"]], "True")
  expect_identical(got[["hydrogen:test_hydrogen_energies_relat:106#4"]], "True")
  expect_identical(got[["hydrogen:test_hydrogen_energies_relat:106#5"]], "True")
  expect_identical(got[["hydrogen:test_hydrogen_energies_relat:108#2"]], "True")
  expect_identical(got[["hydrogen:test_hydrogen_energies_relat:106#6"]], "True")
  expect_identical(got[["hydrogen:test_hydrogen_energies_relat:108#3"]], "True")
  expect_identical(got[["hydrogen:test_hydrogen_energies_relat:106#7"]], "True")
  expect_identical(got[["hydrogen:test_hydrogen_energies_relat:106#8"]], "True")
  expect_identical(got[["hydrogen:test_hydrogen_energies_relat:108#4"]], "True")
  expect_identical(got[["hydrogen:test_hydrogen_energies_relat:106#9"]], "True")
  expect_identical(got[["hydrogen:test_hydrogen_energies_relat:108#5"]], "True")
  expect_identical(got[["hydrogen:test_hydrogen_energies_relat:106#10"]], "True")
  expect_identical(got[["hydrogen:test_hydrogen_energies_relat:108#6"]], "True")
  expect_identical(got[["hydrogen:test_hydrogen_energies_relat:113"]], "True")
  expect_identical(got[["hydrogen:test_hydrogen_energies_relat:113#2"]], "True")
  expect_identical(got[["hydrogen:test_hydrogen_energies_relat:113#3"]], "True")
  expect_identical(got[["hydrogen:test_hydrogen_energies_relat:115"]], "True")
  expect_identical(got[["hydrogen:test_hydrogen_energies_relat:113#4"]], "True")
  expect_identical(got[["hydrogen:test_hydrogen_energies_relat:113#5"]], "True")
  expect_identical(got[["hydrogen:test_hydrogen_energies_relat:115#2"]], "True")
  expect_identical(got[["hydrogen:test_hydrogen_energies_relat:113#6"]], "True")
  expect_identical(got[["hydrogen:test_hydrogen_energies_relat:115#3"]], "True")
  expect_identical(got[["hydrogen:test_hydrogen_energies_relat:113#7"]], "True")
  expect_identical(got[["hydrogen:test_hydrogen_energies_relat:113#8"]], "True")
  expect_identical(got[["hydrogen:test_hydrogen_energies_relat:115#4"]], "True")
  expect_identical(got[["hydrogen:test_hydrogen_energies_relat:113#9"]], "True")
  expect_identical(got[["hydrogen:test_hydrogen_energies_relat:115#5"]], "True")
  expect_identical(got[["hydrogen:test_hydrogen_energies_relat:113#10"]], "True")
  expect_identical(got[["hydrogen:test_hydrogen_energies_relat:115#6"]], "True")
  expect_identical(got[["hydrogen:test_hydrogen_energies_relat:120"]], "True")
  expect_identical(got[["hydrogen:test_hydrogen_energies_relat:120#2"]], "True")
  expect_identical(got[["hydrogen:test_hydrogen_energies_relat:120#3"]], "True")
  expect_identical(got[["hydrogen:test_hydrogen_energies_relat:122"]], "True")
  expect_identical(got[["hydrogen:test_hydrogen_energies_relat:120#4"]], "True")
  expect_identical(got[["hydrogen:test_hydrogen_energies_relat:120#5"]], "True")
  expect_identical(got[["hydrogen:test_hydrogen_energies_relat:122#2"]], "True")
  expect_identical(got[["hydrogen:test_hydrogen_energies_relat:120#6"]], "True")
  expect_identical(got[["hydrogen:test_hydrogen_energies_relat:122#3"]], "True")
  expect_identical(got[["hydrogen:test_hydrogen_energies_relat:120#7"]], "True")
  expect_identical(got[["hydrogen:test_hydrogen_energies_relat:120#8"]], "True")
  expect_identical(got[["hydrogen:test_hydrogen_energies_relat:122#4"]], "True")
  expect_identical(got[["hydrogen:test_hydrogen_energies_relat:120#9"]], "True")
  expect_identical(got[["hydrogen:test_hydrogen_energies_relat:122#5"]], "True")
  expect_identical(got[["hydrogen:test_hydrogen_energies_relat:120#10"]], "True")
  expect_identical(got[["hydrogen:test_hydrogen_energies_relat:122#6"]], "True")
  expect_identical(got[["hydrogen:test_hydrogen_energies_relat:125"]], "ERR ValueError: 'n' must be greater than 'l'")
  expect_identical(got[["hydrogen:test_hydrogen_energies_relat:126"]], "ERR ValueError: 'l' must be positive or zero")
  expect_identical(got[["hydrogen:test_hydrogen_energies_relat:127"]], "ERR ValueError: Spin must be up for l==0.")
})

test_that("hydrogen test_y_lm matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("n r Z"); n <- .py_getitem(.tmpu_, 0L); r <- .py_getitem(.tmpu_, 1L); Z <- .py_getitem(.tmpu_, 2L)
})
feq <- function(a, b, max_relative_error = .py_flt(1e-12), max_absolute_error = .py_flt(1e-12)) {
  a <- .py_float(a)
  b <- .py_float(b)
  if (.py_truthy(.py_lt(.py_abs(.py_sub(a, b)), max_absolute_error))) {
    return(TRUE)
  }
  if (.py_truthy(.py_gt(.py_abs(b), .py_abs(a)))) {
    relative_error <- .py_abs(.py_div(.py_sub(a, b), b))
  } else {
    relative_error <- .py_abs(.py_div(.py_sub(a, b), a))
  }
  return(.py_le(relative_error, max_relative_error))
}
.py_try({
.tmpu_ <- .py_symbols("l, m, phi, theta", real = TRUE); l <- .py_getitem(.tmpu_, 0L); m <- .py_getitem(.tmpu_, 1L); phi <- .py_getitem(.tmpu_, 2L); theta <- .py_getitem(.tmpu_, 3L)
})
.py_try({
.py_rec("hydrogen:test_y_lm:133", function() hydrogen_Y_lm(1L, 0L, 0L, S_Pi))
})
.py_try({
.py_rec("hydrogen:test_y_lm:134", function() hydrogen_Y_lm(0L, 0L, 0L, 0L))
})
.py_try({
.py_rec("hydrogen:test_y_lm:135", function() hydrogen_Y_lm(0L, 0L, phi, theta))
})
.py_try({
.py_rec("hydrogen:test_y_lm:136", function() hydrogen_Y_lm(1L, 0L, phi, theta))
})
.py_try({
.py_rec("hydrogen:test_y_lm:137", function() hydrogen_Y_lm(1L, 1L, phi, theta))
})
.py_try({
.py_rec("hydrogen:test_y_lm:138", function() hydrogen_Y_lm(1L, -1L, phi, theta))
}))---")
  expect_identical(got[["hydrogen:test_y_lm:133"]], "-sqrt(3)/(2*sqrt(pi))")
  expect_identical(got[["hydrogen:test_y_lm:134"]], "1/(2*sqrt(pi))")
  expect_identical(got[["hydrogen:test_y_lm:135"]], "1/(2*sqrt(pi))")
  expect_identical(got[["hydrogen:test_y_lm:136"]], "sqrt(3)*cos(theta)/(2*sqrt(pi))")
  expect_identical(got[["hydrogen:test_y_lm:137"]], "-sqrt(6)*exp(I*phi)*sin(theta)/(4*sqrt(pi))")
  expect_identical(got[["hydrogen:test_y_lm:138"]], "sqrt(6)*exp(-I*phi)*sin(theta)/(4*sqrt(pi))")
})

test_that("hydrogen test_z_lm matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("n r Z"); n <- .py_getitem(.tmpu_, 0L); r <- .py_getitem(.tmpu_, 1L); Z <- .py_getitem(.tmpu_, 2L)
})
feq <- function(a, b, max_relative_error = .py_flt(1e-12), max_absolute_error = .py_flt(1e-12)) {
  a <- .py_float(a)
  b <- .py_float(b)
  if (.py_truthy(.py_lt(.py_abs(.py_sub(a, b)), max_absolute_error))) {
    return(TRUE)
  }
  if (.py_truthy(.py_gt(.py_abs(b), .py_abs(a)))) {
    relative_error <- .py_abs(.py_div(.py_sub(a, b), b))
  } else {
    relative_error <- .py_abs(.py_div(.py_sub(a, b), a))
  }
  return(.py_le(relative_error, max_relative_error))
}
.py_try({
.tmpu_ <- .py_symbols("l, m, phi, theta", real = TRUE); l <- .py_getitem(.tmpu_, 0L); m <- .py_getitem(.tmpu_, 1L); phi <- .py_getitem(.tmpu_, 2L); theta <- .py_getitem(.tmpu_, 3L)
})
.py_try({
.py_rec("hydrogen:test_z_lm:144", function() hydrogen_Z_lm(2L, 0L, 0L, 0L))
})
.py_try({
.py_rec("hydrogen:test_z_lm:145", function() hydrogen_Z_lm(1L, -1L, .py_div(S_Pi, 2L), .py_div(S_Pi, 2L)))
})
.py_try({
.py_rec("hydrogen:test_z_lm:146", function() .py_method(hydrogen_Z_lm(2L, 2L, phi, theta), "simplify"))
})
.py_try({
.py_rec("hydrogen:test_z_lm:148", function() .py_method(hydrogen_Z_lm(2L, -2L, phi, theta), "simplify"))
}))---")
  expect_identical(got[["hydrogen:test_z_lm:144"]], "sqrt(5)/(2*sqrt(pi))")
  expect_identical(got[["hydrogen:test_z_lm:145"]], "-sqrt(3)/(2*sqrt(pi))")
  expect_identical(got[["hydrogen:test_z_lm:146"]], "sqrt(15)*sin(theta)^2*cos(2*phi)/(4*sqrt(pi))")
  expect_identical(got[["hydrogen:test_z_lm:148"]], "-sqrt(15)*sin(2*phi)*sin(theta)^2/(4*sqrt(pi))")
})

test_that("matrices test_parallel_axis_theorem matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
mat1 <- .py_Matrix(.py_tup(.py_tup(2L, -1L, -1L), .py_tup(-1L, 2L, -1L), .py_tup(-1L, -1L, 2L)))
})
.py_try({
.py_rec("matrices:test_parallel_axis_theorem:16", function() pat_matrix(1L, 1L, 1L, 1L))
})
.py_try({
.py_rec("matrices:test_parallel_axis_theorem:17", function() pat_matrix(2L, 1L, 1L, 1L))
})
.py_try({
mat2 <- .py_Matrix(.py_tup(.py_tup(0L, 0L, 0L), .py_tup(0L, 1L, 0L), .py_tup(0L, 0L, 1L)))
})
.py_try({
.py_rec("matrices:test_parallel_axis_theorem:21", function() pat_matrix(1L, 1L, 0L, 0L))
})
.py_try({
.py_rec("matrices:test_parallel_axis_theorem:22", function() pat_matrix(2L, 1L, 0L, 0L))
})
.py_try({
mat3 <- .py_Matrix(.py_tup(.py_tup(1L, 0L, 0L), .py_tup(0L, 0L, 0L), .py_tup(0L, 0L, 1L)))
})
.py_try({
.py_rec("matrices:test_parallel_axis_theorem:26", function() pat_matrix(1L, 0L, 1L, 0L))
})
.py_try({
.py_rec("matrices:test_parallel_axis_theorem:27", function() pat_matrix(2L, 0L, 1L, 0L))
})
.py_try({
mat4 <- .py_Matrix(.py_tup(.py_tup(1L, 0L, 0L), .py_tup(0L, 1L, 0L), .py_tup(0L, 0L, 0L)))
})
.py_try({
.py_rec("matrices:test_parallel_axis_theorem:31", function() pat_matrix(1L, 0L, 0L, 1L))
})
.py_try({
.py_rec("matrices:test_parallel_axis_theorem:32", function() pat_matrix(2L, 0L, 0L, 1L))
}))---")
  expect_identical(got[["matrices:test_parallel_axis_theorem:16"]], "Matrix([[2, -1, -1], [-1, 2, -1], [-1, -1, 2]])")
  expect_identical(got[["matrices:test_parallel_axis_theorem:17"]], "Matrix([[4, -2, -2], [-2, 4, -2], [-2, -2, 4]])")
  expect_identical(got[["matrices:test_parallel_axis_theorem:21"]], "Matrix([[0, 0, 0], [0, 1, 0], [0, 0, 1]])")
  expect_identical(got[["matrices:test_parallel_axis_theorem:22"]], "Matrix([[0, 0, 0], [0, 2, 0], [0, 0, 2]])")
  expect_identical(got[["matrices:test_parallel_axis_theorem:26"]], "Matrix([[1, 0, 0], [0, 0, 0], [0, 0, 1]])")
  expect_identical(got[["matrices:test_parallel_axis_theorem:27"]], "Matrix([[2, 0, 0], [0, 0, 0], [0, 0, 2]])")
  expect_identical(got[["matrices:test_parallel_axis_theorem:31"]], "Matrix([[1, 0, 0], [0, 1, 0], [0, 0, 0]])")
  expect_identical(got[["matrices:test_parallel_axis_theorem:32"]], "Matrix([[2, 0, 0], [0, 2, 0], [0, 0, 0]])")
})

test_that("matrices test_Pauli matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
sigma1 <- msigma(1L)
})
.py_try({
sigma2 <- msigma(2L)
})
.py_try({
sigma3 <- msigma(3L)
})
.py_try({
.py_rec("matrices:test_Pauli:43", function() sigma1)
})
.py_try({
.py_rec("matrices:test_Pauli:44", function() sigma1)
})
.py_try({
.py_rec("matrices:test_Pauli:47", function() .py_mul(sigma1, sigma2))
})
.py_try({
.py_rec("matrices:test_Pauli:48", function() .py_mul(sigma3, sigma1))
})
.py_try({
.py_rec("matrices:test_Pauli:49", function() .py_mul(sigma2, sigma3))
})
.py_try({
.py_rec("matrices:test_Pauli:51", function() .py_mul(sigma1, sigma1))
})
.py_try({
.py_rec("matrices:test_Pauli:52", function() .py_mul(sigma2, sigma2))
})
.py_try({
.py_rec("matrices:test_Pauli:53", function() .py_mul(sigma3, sigma3))
})
.py_try({
.py_rec("matrices:test_Pauli:55", function() .py_mul(.py_mul(sigma1, 2L), sigma1))
})
.py_try({
.py_rec("matrices:test_Pauli:56", function() .py_mul(.py_mul(sigma1, sigma3), sigma1))
}))---")
  expect_identical(got[["matrices:test_Pauli:43"]], "Matrix([[0, 1], [1, 0]])")
  expect_identical(got[["matrices:test_Pauli:44"]], "Matrix([[0, 1], [1, 0]])")
  expect_identical(got[["matrices:test_Pauli:47"]], "Matrix([[I, 0], [0, -I]])")
  expect_identical(got[["matrices:test_Pauli:48"]], "Matrix([[0, 1], [-1, 0]])")
  expect_identical(got[["matrices:test_Pauli:49"]], "Matrix([[0, I], [I, 0]])")
  expect_identical(got[["matrices:test_Pauli:51"]], "Matrix([[1, 0], [0, 1]])")
  expect_identical(got[["matrices:test_Pauli:52"]], "Matrix([[1, 0], [0, 1]])")
  expect_identical(got[["matrices:test_Pauli:53"]], "Matrix([[1, 0], [0, 1]])")
  expect_identical(got[["matrices:test_Pauli:55"]], "Matrix([[2, 0], [0, 2]])")
  expect_identical(got[["matrices:test_Pauli:56"]], "Matrix([[-1, 0], [0, 1]])")
})

test_that("matrices test_Dirac matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
gamma0 <- mgamma(0L)
})
.py_try({
gamma1 <- mgamma(1L)
})
.py_try({
gamma2 <- mgamma(2L)
})
.py_try({
gamma3 <- mgamma(3L)
})
.py_try({
gamma5 <- mgamma(5L)
})
.py_try({
.py_rec("matrices:test_Dirac:67", function() gamma5)
})
.py_try({
.py_rec("matrices:test_Dirac:68", function() .py_add(.py_mul(gamma1, gamma2), .py_mul(gamma2, gamma1)))
})
.py_try({
.py_rec("matrices:test_Dirac:69", function() .py_mul(gamma0, gamma0))
})
.py_try({
.py_rec("matrices:test_Dirac:70", function() .py_mul(gamma2, gamma2))
})
.py_try({
.py_rec("matrices:test_Dirac:71", function() .py_mul(gamma2, gamma2))
})
.py_try({
.py_rec("matrices:test_Dirac:73", function() mgamma(5L, TRUE))
}))---")
  expect_identical(got[["matrices:test_Dirac:67"]], "Matrix([[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]])")
  expect_identical(got[["matrices:test_Dirac:68"]], "Matrix([[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])")
  expect_identical(got[["matrices:test_Dirac:69"]], "Matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])")
  expect_identical(got[["matrices:test_Dirac:70"]], "Matrix([[-1, 0, 0, 0], [0, -1, 0, 0], [0, 0, -1, 0], [0, 0, 0, -1]])")
  expect_identical(got[["matrices:test_Dirac:71"]], "Matrix([[-1, 0, 0, 0], [0, -1, 0, 0], [0, 0, -1, 0], [0, 0, 0, -1]])")
  expect_identical(got[["matrices:test_Dirac:73"]], "Matrix([[0, 0, -1, 0], [0, 0, 0, -1], [-1, 0, 0, 0], [0, -1, 0, 0]])")
})

test_that("matrices test_mdft matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.py_rec("matrices:test_mdft:78", function() mdft(1L))
})
.py_try({
.py_rec("matrices:test_mdft:80", function() mdft(2L))
})
.py_try({
.py_rec("matrices:test_mdft:82", function() mdft(4L))
}))---")
  expect_identical(got[["matrices:test_mdft:78"]], "Matrix([[1]])")
  expect_identical(got[["matrices:test_mdft:80"]], "Matrix([[sqrt(2)/2, sqrt(2)/2], [sqrt(2)/2, -sqrt(2)/2]])")
  expect_identical(got[["matrices:test_mdft:82"]], "Matrix([[1/2, 1/2, 1/2, 1/2], [1/2, -I/2, -1/2, I/2], [1/2, -1/2, 1/2, -1/2], [1/2, I/2, -1/2, -I/2]])")
})

test_that("pring test_wavefunction matches SymPy", {
  got <- .physmisc_battery_run(r"---(m <- Symbol("m")
x <- Symbol("x")
r <- Symbol("r")
.py_try({
Psi <- .py_dict(list(0L, 1L, 2L, 3L), list(.py_div(1L, .py_sqrt(.py_mul(2L, S_Pi))), .py_mul(.py_div(1L, .py_sqrt(.py_mul(2L, S_Pi))), .py_exp(.py_mul(sym_I, x))), .py_mul(.py_div(1L, .py_sqrt(.py_mul(2L, S_Pi))), .py_exp(.py_mul(.py_mul(2L, sym_I), x))), .py_mul(.py_div(1L, .py_sqrt(.py_mul(2L, S_Pi))), .py_exp(.py_mul(.py_mul(3L, sym_I), x)))))
})
.py_try({
for (.it_ in .py_iter(.py_cap(Psi, 60L))) {
  n <- .it_
  .py_try({
  .py_rec("pring:test_wavefunction:20", function() simplify(.py_sub(pring_wavefunction(n, x), .py_getitem(Psi, n))))
  })
}
}))---")
  expect_identical(got[["pring:test_wavefunction:20"]], "0")
  expect_identical(got[["pring:test_wavefunction:20#2"]], "0")
  expect_identical(got[["pring:test_wavefunction:20#3"]], "0")
  expect_identical(got[["pring:test_wavefunction:20#4"]], "0")
})

test_that("pring test_norm matches SymPy", {
  got <- .physmisc_battery_run(r"---(m <- Symbol("m")
x <- Symbol("x")
r <- Symbol("r")
.py_try({
n <- 1L
})
.py_try({
for (.it_ in .py_iter(.py_cap(.py_range(.py_add(n, 1L)), 60L))) {
  i <- .it_
  .py_try({
  .py_rec("pring:test_norm:26", function() sym_integrate(.py_mul(pring_wavefunction(i, x), pring_wavefunction(.py_neg(i), x)), .py_tup(x, 0L, .py_mul(2L, S_Pi))))
  })
}
}))---")
  expect_identical(got[["pring:test_norm:26"]], "1")
  expect_identical(got[["pring:test_norm:26#2"]], "1")
})

test_that("pring test_orthogonality matches SymPy", {
  got <- .physmisc_battery_run(r"---(m <- Symbol("m")
x <- Symbol("x")
r <- Symbol("r")
.py_try({
n <- 1L
})
.py_try({
for (.it_ in .py_iter(.py_cap(.py_range(.py_add(n, 1L)), 60L))) {
  i <- .it_
  .py_try({
  for (.it_ in .py_iter(.py_cap(.py_range(.py_add(i, 1L), .py_add(n, 1L)), 60L))) {
    j <- .it_
    .py_try({
    .py_rec("pring:test_orthogonality:34", function() sym_integrate(.py_mul(pring_wavefunction(i, x), pring_wavefunction(j, x)), .py_tup(x, 0L, .py_mul(2L, S_Pi))))
    })
  }
  })
}
}))---")
  expect_identical(got[["pring:test_orthogonality:34"]], "0")
})

test_that("pring test_energy matches SymPy", {
  got <- .physmisc_battery_run(r"---(m <- Symbol("m")
x <- Symbol("x")
r <- Symbol("r")
.py_try({
n <- 1L
})
.py_try({
for (.it_ in .py_iter(.py_cap(.py_range(.py_add(n, 1L)), 60L))) {
  i <- .it_
  .py_try({
  .py_rec("pring:test_energy:41", function() simplify(.py_sub(pring_energy(i, m, r), .py_div(.py_mul(.py_pow(i, 2L), .py_pow(.py_hbar, 2L)), .py_mul(.py_mul(2L, m), .py_pow(r, 2L))))))
  })
}
}))---")
  expect_identical(got[["pring:test_energy:41"]], "0")
  expect_identical(got[["pring:test_energy:41#2"]], "0")
})

test_that("qho_1d test_wavefunction matches SymPy", {
  got <- .physmisc_battery_run(r"---(omega <- Symbol("omega")
m <- Symbol("m")
x <- Symbol("x")
.py_try({
nu <- .py_div(.py_mul(m, omega), .py_hbar)
})
.py_try({
Psi <- .py_dict(list(0L, 1L, 2L, 3L), list(.py_mul(.py_pow(.py_div(nu, S_Pi), Rational(1L, 4L)), .py_exp(.py_div(.py_mul(.py_neg(nu), .py_pow(x, 2L)), 2L))), .py_mul(.py_mul(.py_mul(.py_pow(.py_div(nu, S_Pi), Rational(1L, 4L)), .py_sqrt(.py_mul(2L, nu))), x), .py_exp(.py_div(.py_mul(.py_neg(nu), .py_pow(x, 2L)), 2L))), .py_mul(.py_div(.py_mul(.py_pow(.py_div(nu, S_Pi), Rational(1L, 4L)), .py_sub(.py_mul(.py_mul(2L, nu), .py_pow(x, 2L)), 1L)), .py_sqrt(2L)), .py_exp(.py_div(.py_mul(.py_neg(nu), .py_pow(x, 2L)), 2L))), .py_mul(.py_mul(.py_mul(.py_pow(.py_div(nu, S_Pi), Rational(1L, 4L)), .py_sqrt(.py_div(nu, 3L))), .py_sub(.py_mul(.py_mul(2L, nu), .py_pow(x, 3L)), .py_mul(3L, x))), .py_exp(.py_div(.py_mul(.py_neg(nu), .py_pow(x, 2L)), 2L)))))
})
.py_try({
for (.it_ in .py_iter(.py_cap(Psi, 60L))) {
  n <- .it_
  .py_try({
  .py_rec("qho_1d:test_wavefunction:24", function() simplify(.py_sub(qho_psi_n(n, x, m, omega), .py_getitem(Psi, n))))
  })
}
}))---")
  expect_identical(got[["qho_1d:test_wavefunction:24"]], "0")
  expect_identical(got[["qho_1d:test_wavefunction:24#2"]], "0")
  expect_identical(got[["qho_1d:test_wavefunction:24#3"]], "0")
  expect_identical(got[["qho_1d:test_wavefunction:24#4"]], "0")
})

test_that("qho_1d test_norm matches SymPy", {
  got <- .physmisc_battery_run(r"---(omega <- Symbol("omega")
m <- Symbol("m")
x <- Symbol("x")
.py_try({
nu <- .py_div(.py_mul(m, omega), .py_hbar)
})
.py_try({
n <- 1L
})
.py_try({
for (.it_ in .py_iter(.py_cap(.py_range(.py_add(n, 1L)), 60L))) {
  i <- .it_
  .py_try({
  .py_rec("qho_1d:test_norm:30", function() sym_integrate(.py_pow(qho_psi_n(i, x, 1L, 1L), 2L), .py_tup(x, .py_neg(oo), oo)))
  })
}
}))---")
  expect_identical(got[["qho_1d:test_norm:30"]], "1")
  expect_identical(got[["qho_1d:test_norm:30#2"]], "1")
})

test_that("qho_1d test_orthogonality matches SymPy", {
  got <- .physmisc_battery_run(r"---(omega <- Symbol("omega")
m <- Symbol("m")
x <- Symbol("x")
.py_try({
nu <- .py_div(.py_mul(m, omega), .py_hbar)
})
.py_try({
n <- 1L
})
.py_try({
for (.it_ in .py_iter(.py_cap(.py_range(.py_add(n, 1L)), 60L))) {
  i <- .it_
  .py_try({
  for (.it_ in .py_iter(.py_cap(.py_range(.py_add(i, 1L), .py_add(n, 1L)), 60L))) {
    j <- .it_
    .py_try({
    .py_rec("qho_1d:test_orthogonality:37", function() sym_integrate(.py_mul(qho_psi_n(i, x, 1L, 1L), qho_psi_n(j, x, 1L, 1L)), .py_tup(x, .py_neg(oo), oo)))
    })
  }
  })
}
}))---")
  expect_identical(got[["qho_1d:test_orthogonality:37"]], "0")
})

test_that("qho_1d test_energies matches SymPy", {
  got <- .physmisc_battery_run(r"---(omega <- Symbol("omega")
m <- Symbol("m")
x <- Symbol("x")
.py_try({
nu <- .py_div(.py_mul(m, omega), .py_hbar)
})
.py_try({
n <- 1L
})
.py_try({
for (.it_ in .py_iter(.py_cap(.py_range(.py_add(n, 1L)), 60L))) {
  i <- .it_
  .py_try({
  .py_rec("qho_1d:test_energies:44", function() qho_E_n(i, omega))
  })
}
}))---")
  expect_identical(got[["qho_1d:test_energies:44"]], "hbar*omega/2")
  expect_identical(got[["qho_1d:test_energies:44#2"]], "3*hbar*omega/2")
})

test_that("qho_1d test_coherent_state matches SymPy", {
  got <- .physmisc_battery_run(r"---(omega <- Symbol("omega")
m <- Symbol("m")
x <- Symbol("x")
.py_try({
nu <- .py_div(.py_mul(m, omega), .py_hbar)
})
.py_try({
n <- 10L
})
.py_try({
alpha <- Symbol("alpha")
})
.py_try({
for (.it_ in .py_iter(.py_cap(.py_range(.py_add(n, 1L)), 60L))) {
  i <- .it_
  .py_try({
  .py_rec("qho_1d:test_coherent_state:51", function() simplify(.py_mul(.py_sqrt(.py_add(n, 1L)), qho_coherent_state(.py_add(n, 1L), alpha))))
  })
}
}))---")
  expect_identical(got[["qho_1d:test_coherent_state:51"]], "sqrt(7)*alpha^11*exp(-Abs(alpha)^2/2)/5040")
  expect_identical(got[["qho_1d:test_coherent_state:51#2"]], "sqrt(7)*alpha^11*exp(-Abs(alpha)^2/2)/5040")
  expect_identical(got[["qho_1d:test_coherent_state:51#3"]], "sqrt(7)*alpha^11*exp(-Abs(alpha)^2/2)/5040")
  expect_identical(got[["qho_1d:test_coherent_state:51#4"]], "sqrt(7)*alpha^11*exp(-Abs(alpha)^2/2)/5040")
  expect_identical(got[["qho_1d:test_coherent_state:51#5"]], "sqrt(7)*alpha^11*exp(-Abs(alpha)^2/2)/5040")
  expect_identical(got[["qho_1d:test_coherent_state:51#6"]], "sqrt(7)*alpha^11*exp(-Abs(alpha)^2/2)/5040")
  expect_identical(got[["qho_1d:test_coherent_state:51#7"]], "sqrt(7)*alpha^11*exp(-Abs(alpha)^2/2)/5040")
  expect_identical(got[["qho_1d:test_coherent_state:51#8"]], "sqrt(7)*alpha^11*exp(-Abs(alpha)^2/2)/5040")
  expect_identical(got[["qho_1d:test_coherent_state:51#9"]], "sqrt(7)*alpha^11*exp(-Abs(alpha)^2/2)/5040")
  expect_identical(got[["qho_1d:test_coherent_state:51#10"]], "sqrt(7)*alpha^11*exp(-Abs(alpha)^2/2)/5040")
  expect_identical(got[["qho_1d:test_coherent_state:51#11"]], "sqrt(7)*alpha^11*exp(-Abs(alpha)^2/2)/5040")
})

test_that("sho test_energy matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("n l hw"); n <- .py_getitem(.tmpu_, 0L); l <- .py_getitem(.tmpu_, 1L); hw <- .py_getitem(.tmpu_, 2L)
})
.py_try({
.py_rec("sho:test_energy:22", function() simplify(.py_sub(sho_E_nl(n, l, hw), .py_mul(.py_add(.py_add(.py_mul(2L, n), l), Rational(3L, 2L)), hw))))
}))---")
  expect_identical(got[["sho:test_energy:22"]], "0")
})

test_that("wigner test_clebsch_gordan_docs matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.py_rec("wigner:test_clebsch_gordan_docs:17", function() clebsch_gordan(Rational(3L, 2L), S_Half, 2L, Rational(3L, 2L), S_Half, 2L))
})
.py_try({
.py_rec("wigner:test_clebsch_gordan_docs:18", function() clebsch_gordan(Rational(3L, 2L), S_Half, 1L, Rational(3L, 2L), Rational(-1L, 2L), 1L))
})
.py_try({
.py_rec("wigner:test_clebsch_gordan_docs:19", function() clebsch_gordan(Rational(3L, 2L), S_Half, 1L, Rational(-1L, 2L), S_Half, 0L))
}))---")
  expect_identical(got[["wigner:test_clebsch_gordan_docs:17"]], "1")
  expect_identical(got[["wigner:test_clebsch_gordan_docs:18"]], "sqrt(3)/2")
  expect_identical(got[["wigner:test_clebsch_gordan_docs:19"]], "-sqrt(2)/2")
})

test_that("wigner test_trival_zero matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
x <- .py_flt(0.5)
})
.py_try({
.py_rec("wigner:test_trival_zero:24", function() clebsch_gordan(x, x, 1L, x, x, -1L))
})
.py_try({
.py_rec("wigner:test_trival_zero:25", function() clebsch_gordan(x, x, x, x, x, x))
})
.py_try({
.py_rec("wigner:test_trival_zero:26", function() clebsch_gordan(1L, 2L, 10L, 0L, 0L, 0L))
})
.py_try({
.py_rec("wigner:test_trival_zero:27", function() clebsch_gordan(-1L, 2L, 3L, 0L, 0L, 0L))
})
.py_try({
.py_rec("wigner:test_trival_zero:28", function() clebsch_gordan(x, x, 1L, .py_flt(1.5), x, 2L))
})
.py_try({
.py_rec("wigner:test_trival_zero:29", function() clebsch_gordan(1L, 2L, 3L, x, x, 1L))
})
.py_try({
.py_rec("wigner:test_trival_zero:30", function() clebsch_gordan(2L, 3L, 4L, 0L, 0L, 0L))
})
.py_try({
.py_rec("wigner:test_trival_zero:32", function() wigner_6j(x, x, x, x, x, x))
})
.py_try({
.py_rec("wigner:test_trival_zero:33", function() wigner_6j(x, x, 1L, x, x, x))
})
.py_try({
.py_rec("wigner:test_trival_zero:34", function() wigner_6j(x, x, 1L, -1L, x, 1L))
})
.py_try({
.py_rec("wigner:test_trival_zero:35", function() wigner_6j(1L, 1L, 2L, x, x, x))
})
.py_try({
.py_rec("wigner:test_trival_zero:37", function() wigner_9j(x, x, x, x, x, x, x, x, x))
})
.py_try({
.py_rec("wigner:test_trival_zero:38", function() wigner_9j(x, x, 1L, x, x, x, x, x, x))
})
.py_try({
.py_rec("wigner:test_trival_zero:39", function() wigner_9j(x, x, 1L, x, x, 1L, x, x, x))
})
.py_try({
.py_rec("wigner:test_trival_zero:40", function() wigner_9j(x, x, 1L, x, x, 1L, x, x, 1L))
})
.py_try({
.py_rec("wigner:test_trival_zero:41", function() wigner_9j(x, x, 1L, x, x, 1L, 1L, x, x))
})
.py_try({
.py_rec("wigner:test_trival_zero:42", function() wigner_9j(x, x, 1L, x, x, 0L, 1L, 1L, 2L))
}))---")
  expect_identical(got[["wigner:test_trival_zero:24"]], "0")
  expect_identical(got[["wigner:test_trival_zero:25"]], "0")
  expect_identical(got[["wigner:test_trival_zero:26"]], "0")
  expect_identical(got[["wigner:test_trival_zero:27"]], "0")
  expect_identical(got[["wigner:test_trival_zero:28"]], "0")
  expect_identical(got[["wigner:test_trival_zero:29"]], "0")
  expect_identical(got[["wigner:test_trival_zero:30"]], "0")
  expect_identical(got[["wigner:test_trival_zero:32"]], "0")
  expect_identical(got[["wigner:test_trival_zero:33"]], "0")
  expect_identical(got[["wigner:test_trival_zero:34"]], "0")
  expect_identical(got[["wigner:test_trival_zero:35"]], "0")
  expect_identical(got[["wigner:test_trival_zero:37"]], "0")
  expect_identical(got[["wigner:test_trival_zero:38"]], "0")
  expect_identical(got[["wigner:test_trival_zero:39"]], "0")
  expect_identical(got[["wigner:test_trival_zero:40"]], "0")
  expect_identical(got[["wigner:test_trival_zero:41"]], "0")
  expect_identical(got[["wigner:test_trival_zero:42"]], "0")
})

test_that("wigner test_clebsch_gordan matches SymPy", {
  got <- .physmisc_battery_run(r"---(tn <- function(a, b) {
  return(.py_lt(.py_method(.py_sub(a, b), "n", 64L), .py_S("1e-64")))
}
.py_try({
h <- S_One
})
.py_try({
k <- S_Half
})
.py_try({
l <- Rational(3L, 2L)
})
.py_try({
i <- Rational(-1L, 2L)
})
.py_try({
n <- Rational(7L, 2L)
})
.py_try({
p <- Rational(5L, 2L)
})
.py_try({
.py_rec("wigner:test_clebsch_gordan:55", function() clebsch_gordan(k, k, 1L, k, k, 1L))
})
.py_try({
.py_rec("wigner:test_clebsch_gordan:56", function() clebsch_gordan(k, k, 1L, k, k, 0L))
})
.py_try({
.py_rec("wigner:test_clebsch_gordan:57", function() clebsch_gordan(k, k, 1L, i, i, -1L))
})
.py_try({
.py_rec("wigner:test_clebsch_gordan:58", function() clebsch_gordan(k, k, 1L, k, i, 0L))
})
.py_try({
.py_rec("wigner:test_clebsch_gordan:59", function() clebsch_gordan(k, k, 0L, k, i, 0L))
})
.py_try({
.py_rec("wigner:test_clebsch_gordan:60", function() clebsch_gordan(k, k, 1L, i, k, 0L))
})
.py_try({
.py_rec("wigner:test_clebsch_gordan:61", function() clebsch_gordan(k, k, 0L, i, k, 0L))
})
.py_try({
.py_rec("wigner:test_clebsch_gordan:62", function() clebsch_gordan(h, k, l, 1L, k, l))
})
.py_try({
.py_rec("wigner:test_clebsch_gordan:63", function() clebsch_gordan(h, k, l, 1L, i, k))
})
.py_try({
.py_rec("wigner:test_clebsch_gordan:64", function() clebsch_gordan(h, k, k, 1L, i, k))
})
.py_try({
.py_rec("wigner:test_clebsch_gordan:65", function() clebsch_gordan(h, k, k, 0L, k, k))
})
.py_try({
.py_rec("wigner:test_clebsch_gordan:66", function() clebsch_gordan(h, k, l, 0L, k, k))
})
.py_try({
.py_rec("wigner:test_clebsch_gordan:67", function() clebsch_gordan(h, h, .py_S(2L), 1L, 1L, .py_S(2L)))
})
.py_try({
.py_rec("wigner:test_clebsch_gordan:68", function() clebsch_gordan(h, h, .py_S(2L), 1L, 0L, 1L))
})
.py_try({
.py_rec("wigner:test_clebsch_gordan:69", function() clebsch_gordan(h, h, .py_S(2L), 0L, 1L, 1L))
})
.py_try({
.py_rec("wigner:test_clebsch_gordan:70", function() clebsch_gordan(h, h, 1L, 1L, 0L, 1L))
})
.py_try({
.py_rec("wigner:test_clebsch_gordan:71", function() clebsch_gordan(h, h, 1L, 0L, 1L, 1L))
})
.py_try({
.py_rec("wigner:test_clebsch_gordan:72", function() clebsch_gordan(l, l, .py_S(3L), l, l, .py_S(3L)))
})
.py_try({
.py_rec("wigner:test_clebsch_gordan:73", function() clebsch_gordan(l, l, .py_S(2L), l, k, .py_S(2L)))
})
.py_try({
.py_rec("wigner:test_clebsch_gordan:74", function() clebsch_gordan(l, l, .py_S(3L), l, k, .py_S(2L)))
})
.py_try({
.py_rec("wigner:test_clebsch_gordan:75", function() clebsch_gordan(.py_S(2L), .py_S(2L), .py_S(4L), .py_S(2L), .py_S(2L), .py_S(4L)))
})
.py_try({
.py_rec("wigner:test_clebsch_gordan:76", function() clebsch_gordan(.py_S(2L), .py_S(2L), .py_S(3L), .py_S(2L), 1L, .py_S(3L)))
})
.py_try({
.py_rec("wigner:test_clebsch_gordan:77", function() clebsch_gordan(.py_S(2L), .py_S(2L), .py_S(3L), 1L, 1L, .py_S(2L)))
})
.py_try({
.py_rec("wigner:test_clebsch_gordan:78", function() clebsch_gordan(p, h, n, p, 1L, n))
})
.py_try({
.py_rec("wigner:test_clebsch_gordan:79", function() clebsch_gordan(p, h, p, p, 0L, p))
})
.py_try({
.py_rec("wigner:test_clebsch_gordan:80", function() clebsch_gordan(p, h, l, k, 1L, l))
})
.py_try({
.py_rec("wigner:test_clebsch_gordan:81", function() tn(clebsch_gordan(5L, 5L, 5L, -1L, 3L, 2L, prec = 64L), .py_div(.py_sqrt(195L), 39L)))
}))---")
  expect_identical(got[["wigner:test_clebsch_gordan:55"]], "1")
  expect_identical(got[["wigner:test_clebsch_gordan:56"]], "0")
  expect_identical(got[["wigner:test_clebsch_gordan:57"]], "1")
  expect_identical(got[["wigner:test_clebsch_gordan:58"]], "sqrt(2)/2")
  expect_identical(got[["wigner:test_clebsch_gordan:59"]], "sqrt(2)/2")
  expect_identical(got[["wigner:test_clebsch_gordan:60"]], "sqrt(2)/2")
  expect_identical(got[["wigner:test_clebsch_gordan:61"]], "-sqrt(2)/2")
  expect_identical(got[["wigner:test_clebsch_gordan:62"]], "1")
  expect_identical(got[["wigner:test_clebsch_gordan:63"]], "sqrt(3)/3")
  expect_identical(got[["wigner:test_clebsch_gordan:64"]], "sqrt(6)/3")
  expect_identical(got[["wigner:test_clebsch_gordan:65"]], "-sqrt(3)/3")
  expect_identical(got[["wigner:test_clebsch_gordan:66"]], "sqrt(6)/3")
  expect_identical(got[["wigner:test_clebsch_gordan:67"]], "1")
  expect_identical(got[["wigner:test_clebsch_gordan:68"]], "sqrt(2)/2")
  expect_identical(got[["wigner:test_clebsch_gordan:69"]], "sqrt(2)/2")
  expect_identical(got[["wigner:test_clebsch_gordan:70"]], "sqrt(2)/2")
  expect_identical(got[["wigner:test_clebsch_gordan:71"]], "-sqrt(2)/2")
  expect_identical(got[["wigner:test_clebsch_gordan:72"]], "1")
  expect_identical(got[["wigner:test_clebsch_gordan:73"]], "sqrt(2)/2")
  expect_identical(got[["wigner:test_clebsch_gordan:74"]], "sqrt(2)/2")
  expect_identical(got[["wigner:test_clebsch_gordan:75"]], "1")
  expect_identical(got[["wigner:test_clebsch_gordan:76"]], "sqrt(2)/2")
  expect_identical(got[["wigner:test_clebsch_gordan:77"]], "0")
  expect_identical(got[["wigner:test_clebsch_gordan:78"]], "1")
  expect_identical(got[["wigner:test_clebsch_gordan:79"]], "sqrt(35)/7")
  expect_identical(got[["wigner:test_clebsch_gordan:80"]], "sqrt(15)/15")
  expect_identical(got[["wigner:test_clebsch_gordan:81"]], "True")
})

test_that("wigner test_clebsch_gordan_numpy matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
NULL
})
.py_try({
.py_rec("wigner:test_clebsch_gordan_numpy:89", function() .py_call(clebsch_gordan, .py_splat(.py_method(.py_np_zeros(6L), "astype", .py_np_int64))))
})
.py_try({
.py_rec("wigner:test_clebsch_gordan_numpy:90", function() wigner_3j(2L, .py_np_float64(.py_flt(6.0)), .py_flt(4.0), 0L, 0L, 0L))
})
.py_try({
.py_rec("wigner:test_clebsch_gordan_numpy:91", function() wigner_3j(0L, .py_flt(0.5), .py_flt(0.5), 0L, .py_flt(0.5), .py_flt(-0.5)))
})
.py_try({
.py_rec("wigner:test_clebsch_gordan_numpy:92", function() wigner_3j(.py_flt(2.1), 6L, 4L, 0L, 0L, 0L))
}))---")
  expect_identical(got[["wigner:test_clebsch_gordan_numpy:89"]], "1")
  expect_identical(got[["wigner:test_clebsch_gordan_numpy:90"]], "sqrt(715)/143")
  expect_identical(got[["wigner:test_clebsch_gordan_numpy:91"]], "sqrt(2)/2")
  expect_identical(got[["wigner:test_clebsch_gordan_numpy:92"]], "ERR ValueError: expecting integer or half-integer, got 2.1")
})

test_that("wigner test_wigner matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
NULL
})
tn <- function(a, b) {
  return(.py_lt(.py_method(.py_sub(a, b), "n", 64L), .py_S("1e-64")))
}
.py_try({
.py_rec("wigner:test_wigner:102", function() tn(wigner_9j(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 0L, prec = 64L), Rational(1L, 18L)))
})
.py_try({
.py_rec("wigner:test_wigner:103", function() wigner_9j(3L, 3L, 2L, 3L, 3L, 2L, 3L, 3L, 2L))
})
.py_try({
.py_rec("wigner:test_wigner:105", function() wigner_6j(5L, 5L, 5L, 5L, 5L, 5L))
})
.py_try({
.py_rec("wigner:test_wigner:106", function() wigner_6j(2L, 2L, 2L, .py_flt(1.5), .py_flt(1.5), .py_flt(1.5)))
})
.py_try({
.py_rec("wigner:test_wigner:107", function() tn(wigner_6j(8L, 8L, 8L, 8L, 8L, 8L, prec = 64L), Rational(-12219L, 965770L)))
})
.py_try({
.py_rec("wigner:test_wigner:108", function() wigner_6j(1L, 1L, 1L, .py_flt(1.0), .py_np_float64(.py_flt(1.0)), 1L))
})
.py_try({
.py_rec("wigner:test_wigner:109", function() wigner_6j(.py_flt(3.0), .py_np_float32(3L), .py_flt(3.0), 3L, 3L, 3L))
})
.py_try({
half <- S_Half
})
.py_try({
.py_rec("wigner:test_wigner:112", function() wigner_9j(0L, 0L, 0L, 0L, half, half, 0L, half, half))
})
.py_try({
.py_rec("wigner:test_wigner:113", function() wigner_9j(3L, 5L, 4L, .py_mul(7L, half), .py_mul(5L, half), 4L, .py_mul(9L, half), .py_mul(9L, half), 0L))
})
.py_try({
.py_rec("wigner:test_wigner:117", function() wigner_9j(1L, 4L, 3L, .py_mul(5L, half), 4L, .py_mul(5L, half), .py_mul(5L, half), 2L, .py_mul(7L, half)))
})
.py_try({
.py_rec("wigner:test_wigner:121", function() wigner_9j(4L, .py_mul(9L, half), .py_mul(5L, half), 2L, 4L, 4L, 5L, .py_mul(7L, half), .py_mul(7L, half)))
})
.py_try({
.py_rec("wigner:test_wigner:125", function() wigner_9j(5L, 5L, .py_flt(5.0), .py_np_float64(.py_flt(5.0)), 5L, 5L, 5L, 5L, 5L))
})
.py_try({
.py_rec("wigner:test_wigner:129", function() wigner_9j(.py_flt(1.0), .py_flt(2.0), .py_flt(3.0), 3L, 2L, 1L, 2L, 1L, 3L))
}))---")
  expect_identical(got[["wigner:test_wigner:102"]], "True")
  expect_identical(got[["wigner:test_wigner:103"]], "17*sqrt(6)/4410")
  expect_identical(got[["wigner:test_wigner:105"]], "1/52")
  expect_identical(got[["wigner:test_wigner:106"]], "0")
  expect_identical(got[["wigner:test_wigner:107"]], "True")
  expect_identical(got[["wigner:test_wigner:108"]], "1/6")
  expect_identical(got[["wigner:test_wigner:109"]], "-1/14")
  expect_identical(got[["wigner:test_wigner:112"]], "1/2")
  expect_identical(got[["wigner:test_wigner:113"]], "-19*sqrt(210)/207900")
  expect_identical(got[["wigner:test_wigner:117"]], "-19*sqrt(330)/105840")
  expect_identical(got[["wigner:test_wigner:121"]], "-59*sqrt(105)/727650")
  expect_identical(got[["wigner:test_wigner:125"]], "0")
  expect_identical(got[["wigner:test_wigner:129"]], "-4*sqrt(70)/11025")
})

test_that("wigner test_gaunt matches SymPy", {
  got <- .physmisc_battery_run(r"---(tn <- function(a, b) {
  return(.py_lt(.py_method(.py_sub(a, b), "n", 64L), .py_S("1e-64")))
}
.py_try({
.py_rec("wigner:test_gaunt:139", function() gaunt(1L, 0L, 1L, 1L, 0L, -1L))
})
.py_try({
.py_rec("wigner:test_gaunt:140", function() .py_isinstance(.py_getitem(.py_attr(gaunt(1L, 1L, 0L, -1L, 1L, 0L), "args"), 0L), Rational))
})
.py_try({
.py_rec("wigner:test_gaunt:141", function() .py_isinstance(.py_getitem(.py_attr(gaunt(0L, 1L, 1L, 0L, -1L, 1L), "args"), 0L), Rational))
})
.py_try({
.py_rec("wigner:test_gaunt:143", function() tn(gaunt(10L, 10L, 12L, 9L, 3L, -12L, prec = 64L), .py_div(.py_mul(Rational(-98L, 62031L), .py_sqrt(6279L)), .py_sqrt(S_Pi))))
})
gaunt_ref <- function(l1, l2, l3, m1, m2, m3) {
  return(.py_mul(.py_mul(.py_sqrt(.py_div(.py_mul(.py_mul(.py_add(.py_mul(2L, l1), 1L), .py_add(.py_mul(2L, l2), 1L)), .py_add(.py_mul(2L, l3), 1L)), .py_mul(4L, S_Pi))), wigner_3j(l1, l2, l3, 0L, 0L, 0L)), wigner_3j(l1, l2, l3, m1, m2, m3)))
}
.py_try({
threshold <- .py_flt(1e-10)
})
.py_try({
l_max <- 3L
})
.py_try({
l3_max <- 24L
})
.py_try({
for (.it_ in .py_iter(.py_cap(.py_range(.py_add(l_max, 1L)), 2L))) {
  l1 <- .it_
  .py_try({
  for (.it_ in .py_iter(.py_cap(.py_range(.py_add(l_max, 1L)), 2L))) {
    l2 <- .it_
    .py_try({
    for (.it_ in .py_iter(.py_cap(.py_range(.py_add(l3_max, 1L)), 2L))) {
      l3 <- .it_
      .py_try({
      for (.it_ in .py_iter(.py_cap(.py_range(.py_neg(l1), .py_add(l1, 1L)), 2L))) {
        m1 <- .it_
        .py_try({
        for (.it_ in .py_iter(.py_cap(.py_range(.py_neg(l2), .py_add(l2, 1L)), 2L))) {
          m2 <- .it_
          .py_try({
          for (.it_ in .py_iter(.py_cap(.py_range(.py_neg(l3), .py_add(l3, 1L)), 2L))) {
            m3 <- .it_
            .py_try({
            args <- .py_tup(l1, l2, l3, m1, m2, m3)
            })
            .py_try({
            g <- .py_call(gaunt, .py_splat(args))
            })
            .py_try({
            g0 <- .py_call(gaunt_ref, .py_splat(args))
            })
            .py_try({
            .py_rec("wigner:test_gaunt:163", function() .py_lt(.py_abs(.py_sub(g, g0)), threshold))
            })
            .py_try({
            if (.py_truthy(.py_ne(.py_add(.py_add(m1, m2), m3), 0L))) {
              .py_try({
              .py_rec("wigner:test_gaunt:165", function() .py_lt(.py_abs(g), threshold))
              })
            }
            })
            .py_try({
            if (.py_truthy(.py_mod(.py_add(.py_add(l1, l2), l3), 2L))) {
              .py_try({
              .py_rec("wigner:test_gaunt:167", function() .py_lt(.py_abs(g), threshold))
              })
            }
            })
          }
          })
        }
        })
      }
      })
    }
    })
  }
  })
}
})
.py_try({
.py_rec("wigner:test_gaunt:168", function() gaunt(1L, 1L, 0L, 0L, 2L, -2L))
}))---")
  expect_identical(got[["wigner:test_gaunt:139"]], "-1/(2*sqrt(pi))")
  expect_identical(got[["wigner:test_gaunt:140"]], "True")
  expect_identical(got[["wigner:test_gaunt:141"]], "True")
  expect_identical(got[["wigner:test_gaunt:143"]], "True")
  expect_identical(got[["wigner:test_gaunt:163"]], "True")
  expect_identical(got[["wigner:test_gaunt:163#2"]], "True")
  expect_identical(got[["wigner:test_gaunt:165"]], "True")
  expect_identical(got[["wigner:test_gaunt:167"]], "True")
  expect_identical(got[["wigner:test_gaunt:163#3"]], "True")
  expect_identical(got[["wigner:test_gaunt:167#2"]], "True")
  expect_identical(got[["wigner:test_gaunt:163#4"]], "True")
  expect_identical(got[["wigner:test_gaunt:165#2"]], "True")
  expect_identical(got[["wigner:test_gaunt:167#3"]], "True")
  expect_identical(got[["wigner:test_gaunt:163#5"]], "True")
  expect_identical(got[["wigner:test_gaunt:167#4"]], "True")
  expect_identical(got[["wigner:test_gaunt:163#6"]], "True")
  expect_identical(got[["wigner:test_gaunt:165#3"]], "True")
  expect_identical(got[["wigner:test_gaunt:163#7"]], "True")
  expect_identical(got[["wigner:test_gaunt:165#4"]], "True")
  expect_identical(got[["wigner:test_gaunt:163#8"]], "True")
  expect_identical(got[["wigner:test_gaunt:165#5"]], "True")
  expect_identical(got[["wigner:test_gaunt:163#9"]], "True")
  expect_identical(got[["wigner:test_gaunt:163#10"]], "True")
  expect_identical(got[["wigner:test_gaunt:165#6"]], "True")
  expect_identical(got[["wigner:test_gaunt:167#5"]], "True")
  expect_identical(got[["wigner:test_gaunt:163#11"]], "True")
  expect_identical(got[["wigner:test_gaunt:167#6"]], "True")
  expect_identical(got[["wigner:test_gaunt:163#12"]], "True")
  expect_identical(got[["wigner:test_gaunt:165#7"]], "True")
  expect_identical(got[["wigner:test_gaunt:163#13"]], "True")
  expect_identical(got[["wigner:test_gaunt:165#8"]], "True")
  expect_identical(got[["wigner:test_gaunt:163#14"]], "True")
  expect_identical(got[["wigner:test_gaunt:165#9"]], "True")
  expect_identical(got[["wigner:test_gaunt:163#15"]], "True")
  expect_identical(got[["wigner:test_gaunt:163#16"]], "True")
  expect_identical(got[["wigner:test_gaunt:165#10"]], "True")
  expect_identical(got[["wigner:test_gaunt:163#17"]], "True")
  expect_identical(got[["wigner:test_gaunt:165#11"]], "True")
  expect_identical(got[["wigner:test_gaunt:163#18"]], "True")
  expect_identical(got[["wigner:test_gaunt:165#12"]], "True")
  expect_identical(got[["wigner:test_gaunt:163#19"]], "True")
  expect_identical(got[["wigner:test_gaunt:163#20"]], "True")
  expect_identical(got[["wigner:test_gaunt:165#13"]], "True")
  expect_identical(got[["wigner:test_gaunt:167#7"]], "True")
  expect_identical(got[["wigner:test_gaunt:163#21"]], "True")
  expect_identical(got[["wigner:test_gaunt:165#14"]], "True")
  expect_identical(got[["wigner:test_gaunt:167#8"]], "True")
  expect_identical(got[["wigner:test_gaunt:163#22"]], "True")
  expect_identical(got[["wigner:test_gaunt:165#15"]], "True")
  expect_identical(got[["wigner:test_gaunt:167#9"]], "True")
  expect_identical(got[["wigner:test_gaunt:163#23"]], "True")
  expect_identical(got[["wigner:test_gaunt:165#16"]], "True")
  expect_identical(got[["wigner:test_gaunt:167#10"]], "True")
  expect_identical(got[["wigner:test_gaunt:163#24"]], "True")
  expect_identical(got[["wigner:test_gaunt:165#17"]], "True")
  expect_identical(got[["wigner:test_gaunt:167#11"]], "True")
  expect_identical(got[["wigner:test_gaunt:163#25"]], "True")
  expect_identical(got[["wigner:test_gaunt:165#18"]], "True")
  expect_identical(got[["wigner:test_gaunt:167#12"]], "True")
  expect_identical(got[["wigner:test_gaunt:163#26"]], "True")
  expect_identical(got[["wigner:test_gaunt:165#19"]], "True")
  expect_identical(got[["wigner:test_gaunt:167#13"]], "True")
  expect_identical(got[["wigner:test_gaunt:163#27"]], "True")
  expect_identical(got[["wigner:test_gaunt:167#14"]], "True")
  expect_identical(got[["wigner:test_gaunt:168"]], "0")
})

test_that("wigner test_realgaunt matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
for (.it_ in .py_iter(.py_cap(.py_range(3L), 60L))) {
  l <- .it_
  .py_try({
  for (.it_ in .py_iter(.py_cap(.py_range(.py_neg(l), .py_add(l, 1L)), 60L))) {
    m <- .it_
    .py_try({
    .py_rec("wigner:test_realgaunt:175", function() real_gaunt(0L, l, l, 0L, m, m))
    })
  }
  })
}
})
.py_try({
.py_rec("wigner:test_realgaunt:176", function() real_gaunt(1L, 1L, 2L, 0L, 0L, 0L))
})
.py_try({
.py_rec("wigner:test_realgaunt:177", function() real_gaunt(1L, 1L, 2L, 1L, 1L, 0L))
})
.py_try({
.py_rec("wigner:test_realgaunt:178", function() real_gaunt(2L, 2L, 2L, 0L, 0L, 0L))
})
.py_try({
.py_rec("wigner:test_realgaunt:179", function() real_gaunt(2L, 2L, 2L, 0L, 2L, 2L))
})
.py_try({
.py_rec("wigner:test_realgaunt:180", function() real_gaunt(2L, 2L, 2L, -2L, -2L, 0L))
})
.py_try({
.py_rec("wigner:test_realgaunt:181", function() real_gaunt(1L, 1L, 2L, -1L, 0L, -1L))
})
.py_try({
.py_rec("wigner:test_realgaunt:182", function() real_gaunt(1L, 1L, 2L, 0L, 1L, 1L))
})
.py_try({
.py_rec("wigner:test_realgaunt:183", function() real_gaunt(1L, 1L, 2L, 1L, 1L, 2L))
})
.py_try({
.py_rec("wigner:test_realgaunt:184", function() real_gaunt(1L, 1L, 2L, -1L, 1L, -2L))
})
.py_try({
.py_rec("wigner:test_realgaunt:185", function() real_gaunt(1L, 1L, 2L, -1L, -1L, 2L))
})
.py_try({
.py_rec("wigner:test_realgaunt:186", function() real_gaunt(2L, 2L, 2L, 0L, 1L, 1L))
})
.py_try({
.py_rec("wigner:test_realgaunt:187", function() real_gaunt(2L, 2L, 2L, 1L, 1L, 2L))
})
.py_try({
.py_rec("wigner:test_realgaunt:188", function() real_gaunt(2L, 2L, 2L, -1L, -1L, 2L))
})
.py_try({
.py_rec("wigner:test_realgaunt:190", function() real_gaunt(-2L, -2L, -2L, -2L, -2L, 0L))
})
.py_try({
.py_rec("wigner:test_realgaunt:191", function() real_gaunt(-2L, 1L, 0L, 1L, 1L, 1L))
})
.py_try({
.py_rec("wigner:test_realgaunt:192", function() real_gaunt(-2L, -1L, -2L, -1L, -1L, 0L))
})
.py_try({
.py_rec("wigner:test_realgaunt:193", function() real_gaunt(-2L, -2L, -2L, -2L, -2L, -2L))
})
.py_try({
.py_rec("wigner:test_realgaunt:194", function() real_gaunt(-2L, -1L, -2L, -1L, -1L, -1L))
})
.py_try({
x <- .py_symbols("x", integer = TRUE)
})
.py_try({
v <- .py_mul(list(0L), 6L)
})
.py_try({
for (.it_ in .py_iter(.py_cap(.py_range(.py_len(v)), 60L))) {
  i <- .it_
  .py_try({
  v <- .py_setitem(v, i, x)
  })
  .py_try({
  .py_rec("wigner:test_realgaunt:200", function() .py_call(real_gaunt, .py_splat(v)))
  })
  .py_try({
  v <- .py_setitem(v, i, 0L)
  })
}
}))---")
  expect_identical(got[["wigner:test_realgaunt:175"]], "1/(2*sqrt(pi))")
  expect_identical(got[["wigner:test_realgaunt:175#2"]], "1/(2*sqrt(pi))")
  expect_identical(got[["wigner:test_realgaunt:175#3"]], "1/(2*sqrt(pi))")
  expect_identical(got[["wigner:test_realgaunt:175#4"]], "1/(2*sqrt(pi))")
  expect_identical(got[["wigner:test_realgaunt:175#5"]], "1/(2*sqrt(pi))")
  expect_identical(got[["wigner:test_realgaunt:175#6"]], "1/(2*sqrt(pi))")
  expect_identical(got[["wigner:test_realgaunt:175#7"]], "1/(2*sqrt(pi))")
  expect_identical(got[["wigner:test_realgaunt:175#8"]], "1/(2*sqrt(pi))")
  expect_identical(got[["wigner:test_realgaunt:175#9"]], "1/(2*sqrt(pi))")
  expect_identical(got[["wigner:test_realgaunt:176"]], "sqrt(5)/(5*sqrt(pi))")
  expect_identical(got[["wigner:test_realgaunt:177"]], "-sqrt(5)/(10*sqrt(pi))")
  expect_identical(got[["wigner:test_realgaunt:178"]], "sqrt(5)/(7*sqrt(pi))")
  expect_identical(got[["wigner:test_realgaunt:179"]], "-sqrt(5)/(7*sqrt(pi))")
  expect_identical(got[["wigner:test_realgaunt:180"]], "-sqrt(5)/(7*sqrt(pi))")
  expect_identical(got[["wigner:test_realgaunt:181"]], "sqrt(15)/(10*sqrt(pi))")
  expect_identical(got[["wigner:test_realgaunt:182"]], "sqrt(15)/(10*sqrt(pi))")
  expect_identical(got[["wigner:test_realgaunt:183"]], "sqrt(15)/(10*sqrt(pi))")
  expect_identical(got[["wigner:test_realgaunt:184"]], "sqrt(15)/(10*sqrt(pi))")
  expect_identical(got[["wigner:test_realgaunt:185"]], "-sqrt(15)/(10*sqrt(pi))")
  expect_identical(got[["wigner:test_realgaunt:186"]], "sqrt(5)/(14*sqrt(pi))")
  expect_identical(got[["wigner:test_realgaunt:187"]], "sqrt(15)/(14*sqrt(pi))")
  expect_identical(got[["wigner:test_realgaunt:188"]], "-sqrt(15)/(14*sqrt(pi))")
  expect_identical(got[["wigner:test_realgaunt:190"]], "0")
  expect_identical(got[["wigner:test_realgaunt:191"]], "0")
  expect_identical(got[["wigner:test_realgaunt:192"]], "0")
  expect_identical(got[["wigner:test_realgaunt:193"]], "0")
  expect_identical(got[["wigner:test_realgaunt:194"]], "0")
  expect_identical(got[["wigner:test_realgaunt:200"]], "ERR ValueError: x is not an integer")
  expect_identical(got[["wigner:test_realgaunt:200#2"]], "ERR ValueError: x is not an integer")
  expect_identical(got[["wigner:test_realgaunt:200#3"]], "ERR ValueError: x is not an integer")
  expect_identical(got[["wigner:test_realgaunt:200#4"]], "ERR ValueError: x is not an integer")
  expect_identical(got[["wigner:test_realgaunt:200#5"]], "ERR ValueError: x is not an integer")
  expect_identical(got[["wigner:test_realgaunt:200#6"]], "ERR ValueError: x is not an integer")
})

test_that("wigner test_racah matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.py_rec("wigner:test_racah:205", function() racah(2L, 2L, 2L, 2L, 1L, 3L))
})
.py_try({
.py_rec("wigner:test_racah:206", function() racah(3L, 3L, 3L, 3L, 3L, 3L))
})
.py_try({
.py_rec("wigner:test_racah:207", function() racah(2L, 2L, 2L, 2L, 2L, 2L))
})
.py_try({
.py_rec("wigner:test_racah:208", function() .py_attr(racah(7L, 8L, 7L, 1L, 7L, 7L, prec = 4L), "is_Float"))
})
.py_try({
.py_rec("wigner:test_racah:209", function() racah(.py_flt(5.5), .py_flt(7.5), .py_flt(9.5), .py_flt(6.5), 8L, 9L))
})
.py_try({
.py_rec("wigner:test_racah:210", function() .py_lt(.py_abs(.py_sub(racah(.py_flt(5.5), .py_flt(7.5), .py_flt(9.5), .py_flt(6.5), 8L, 9L, prec = 4L), .py_flt(-0.01517))), .py_S("1e-4")))
}))---")
  expect_identical(got[["wigner:test_racah:205"]], "0")
  expect_identical(got[["wigner:test_racah:206"]], "-1/14")
  expect_identical(got[["wigner:test_racah:207"]], "-3/70")
  expect_identical(got[["wigner:test_racah:208"]], "True")
  expect_identical(got[["wigner:test_racah:209"]], "-719*sqrt(598)/1158924")
  expect_identical(got[["wigner:test_racah:210"]], "True")
})

test_that("wigner test_dot_rota_grad_SH matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
.tmpu_ <- .py_symbols("theta phi"); theta <- .py_getitem(.tmpu_, 0L); phi <- .py_getitem(.tmpu_, 1L)
})
.py_try({
.py_rec("wigner:test_dot_rota_grad_SH:215", function() dot_rot_grad_Ynm(1L, 1L, 1L, 1L, 1L, 0L))
})
.py_try({
.py_rec("wigner:test_dot_rota_grad_SH:217", function() .py_method(dot_rot_grad_Ynm(1L, 1L, 1L, 1L, 1L, 0L), "doit"))
})
.py_try({
.py_rec("wigner:test_dot_rota_grad_SH:219", function() dot_rot_grad_Ynm(1L, 5L, 1L, 1L, 1L, 2L))
})
.py_try({
.py_rec("wigner:test_dot_rota_grad_SH:221", function() .py_method(dot_rot_grad_Ynm(1L, 5L, 1L, 1L, 1L, 2L), "doit"))
})
.py_try({
.py_rec("wigner:test_dot_rota_grad_SH:223", function() .py_method(dot_rot_grad_Ynm(3L, 3L, 3L, 3L, theta, phi), "doit"))
})
.py_try({
.py_rec("wigner:test_dot_rota_grad_SH:225", function() .py_method(dot_rot_grad_Ynm(3L, 3L, 1L, 1L, theta, phi), "doit"))
})
.py_try({
.py_rec("wigner:test_dot_rota_grad_SH:227", function() .py_method(dot_rot_grad_Ynm(3L, 2L, 2L, 0L, theta, phi), "doit"))
})
.py_try({
.py_rec("wigner:test_dot_rota_grad_SH:229", function() .py_method(.py_method(dot_rot_grad_Ynm(3L, 2L, 3L, 2L, theta, phi), "doit"), "expand"))
}))---")
  expect_identical(got[["wigner:test_dot_rota_grad_SH:215"]], "Sum(sqrt(18*_k + 9)*(_k^2 + _k - 4)*Wigner3j(1, 1, _k, 0, 0, 0)*Wigner3j(1, 1, _k, 1, 1, -2)*Ynm(_k, 2, 1, 0)/(4*sqrt(pi)), (_k, 0, 2))")
  expect_identical(got[["wigner:test_dot_rota_grad_SH:217"]], "sqrt(30)*Ynm(2, 2, 1, 0)/(10*sqrt(pi))")
  expect_identical(got[["wigner:test_dot_rota_grad_SH:219"]], "Sum(sqrt(18*_k + 9)*(_k^2 + _k - 4)*Wigner3j(1, 1, _k, 0, 0, 0)*Wigner3j(1, 1, _k, 5, 1, -6)*Ynm(_k, 6, 1, 2)/(4*sqrt(pi)), (_k, 0, 2))")
  expect_identical(got[["wigner:test_dot_rota_grad_SH:221"]], "0")
  expect_identical(got[["wigner:test_dot_rota_grad_SH:223"]], "15*sqrt(3003)*Ynm(6, 6, theta, phi)/(143*sqrt(pi))")
  expect_identical(got[["wigner:test_dot_rota_grad_SH:225"]], "sqrt(3)*Ynm(4, 4, theta, phi)/sqrt(pi)")
  expect_identical(got[["wigner:test_dot_rota_grad_SH:227"]], "3*sqrt(55)*Ynm(5, 2, theta, phi)/(11*sqrt(pi))")
  expect_identical(got[["wigner:test_dot_rota_grad_SH:229"]], "-sqrt(70)*Ynm(4, 4, theta, phi)/(11*sqrt(pi)) + 45*sqrt(182)*Ynm(6, 4, theta, phi)/(143*sqrt(pi))")
})

test_that("wigner test_wigner_d matches SymPy", {
  got <- .physmisc_battery_run(r"---(.py_try({
half <- .py_div(.py_S(1L), 2L)
})
.py_try({
.py_rec("wigner:test_wigner_d:236", function() wigner_d_small(half, 0L))
})
.py_try({
.py_rec("wigner:test_wigner_d:237", function() wigner_d_small(half, .py_div(S_Pi, 2L)))
})
.py_try({
.py_rec("wigner:test_wigner_d:238", function() wigner_d_small(half, S_Pi))
})
.py_try({
.tmpu_ <- .py_symbols("alpha, beta, gamma", real = TRUE); alpha <- .py_getitem(.tmpu_, 0L); beta <- .py_getitem(.tmpu_, 1L); gamma <- .py_getitem(.tmpu_, 2L)
})
.py_try({
D <- wigner_d(half, alpha, beta, gamma)
})
.py_try({
.py_rec("wigner:test_wigner_d:242", function() .py_getitem(D, .py_tup(0L, 0L)))
})
.py_try({
.py_rec("wigner:test_wigner_d:243", function() .py_getitem(D, .py_tup(0L, 1L)))
})
.py_try({
.py_rec("wigner:test_wigner_d:244", function() .py_getitem(D, .py_tup(1L, 0L)))
})
.py_try({
.py_rec("wigner:test_wigner_d:245", function() .py_getitem(D, .py_tup(1L, 1L)))
})
.py_try({
.tmpu_ <- .py_symbols("theta phi", real = TRUE); theta <- .py_getitem(.tmpu_, 0L); phi <- .py_getitem(.tmpu_, 1L)
})
.py_try({
v <- .py_Matrix(.py_comp(.py_range(1L, -2L, -1L), function(.it) { mj <- .it; list(Ynm(1L, mj, theta, phi)) }))
})
.py_try({
w <- .py_matmul(wigner_d(1L, .py_div(.py_neg(S_Pi), 2L), .py_div(S_Pi, 2L), .py_div(.py_neg(S_Pi), 2L)), .py_method(v, "subs", .py_dict(list(theta, phi), list(.py_div(S_Pi, 4L), S_Pi))))
})
.py_try({
w_ <- .py_method(v, "subs", .py_dict(list(theta, phi), list(.py_div(S_Pi, 2L), .py_div(S_Pi, 4L))))
})
.py_try({
.py_rec("wigner:test_wigner_d:252", function() .py_method(.py_method(w, "expand", func = TRUE), "as_real_imag"))
}))---")
  expect_identical(got[["wigner:test_wigner_d:236"]], "Matrix([[1, 0], [0, 1]])")
  expect_identical(got[["wigner:test_wigner_d:237"]], "Matrix([[sqrt(2)/2, sqrt(2)/2], [-sqrt(2)/2, sqrt(2)/2]])")
  expect_identical(got[["wigner:test_wigner_d:238"]], "Matrix([[0, 1], [-1, 0]])")
  expect_identical(got[["wigner:test_wigner_d:242"]], "exp(I*alpha/2)*exp(I*gamma/2)*cos(beta/2)")
  expect_identical(got[["wigner:test_wigner_d:243"]], "exp(I*alpha/2)*exp(-I*gamma/2)*sin(beta/2)")
  expect_identical(got[["wigner:test_wigner_d:244"]], "-exp(-I*alpha/2)*exp(I*gamma/2)*sin(beta/2)")
  expect_identical(got[["wigner:test_wigner_d:245"]], "exp(-I*alpha/2)*exp(-I*gamma/2)*cos(beta/2)")
  expect_identical(got[["wigner:test_wigner_d:252"]], "(Matrix([[-sqrt(3)/(4*sqrt(pi))], [0], [sqrt(3)/(4*sqrt(pi))]]), Matrix([[-sqrt(3)/(4*sqrt(pi))], [0], [-sqrt(3)/(4*sqrt(pi))]]))")
})
