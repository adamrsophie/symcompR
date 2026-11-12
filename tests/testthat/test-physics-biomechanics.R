test_that("characteristic curves print, release and differentiate as in SymPy", {
  l <- Symbol("l_T_tilde")
  cs <- lapply(0:3, function(i) Symbol(paste0("c_", i)))
  f <- do.call(TendonForceLengthDeGroote2016, c(list(l), cs))
  expect_equal(sstr(f), "TendonForceLengthDeGroote2016(l_T_tilde, c_0, c_1, c_2, c_3)")
  expect_equal(sstr(doit(f)), "c_0*exp(c_3*(-c_1 + l_T_tilde)) - c_2")
  expect_equal(sstr(diff(f, cs[[4L]])), "c_0*(-c_1 + l_T_tilde)*exp(c_3*(-c_1 + l_T_tilde))")
  expect_identical(inverse(f), TendonForceLengthInverseDeGroote2016)
  expect_equal(sstr(TendonForceLengthDeGroote2016_with_defaults(l)),
               "TendonForceLengthDeGroote2016(l_T_tilde, 0.2, 0.995, 0.25, 33.93669377311689)")
  lm <- Symbol("l_M_tilde")
  p <- FiberForceLengthPassiveDeGroote2016(lm, cs[[1L]], cs[[2L]])
  expect_equal(sstr(doit(p, evaluate = FALSE)), "(-1 + exp(c_1*(l_M_tilde - 1)/c_0))/(exp(c_1) - 1)")
  expect_error(CharacteristicCurveFunction(), "Cannot directly instantiate")
})

test_that("activation dynamics follow SymPy", {
  act <- FirstOrderActivationDeGroote2016("name", Symbol("tau_a"), Symbol("tau_d"), Symbol("b"))
  expect_equal(format(act), paste0("FirstOrderActivationDeGroote2016('name', activation_time_constant=tau_a, ",
                                   "deactivation_time_constant=tau_d, smoothing_rate=b)"))
  expect_equal(sstr(act$p), "Matrix([[tau_a], [tau_d], [b]])")
  expect_equal(sstr(act$x), "Matrix([[a_name(t)]])")
  z <- ZerothOrderActivation("name")
  expect_equal(sstr(z$M), "Matrix(0, 0, [])")
  expect_equal(sstr(z$a), "e_name(t)")
  expect_error(z$e <- 1, "has no setter")
  expect_error(ActivationBase("x"), "Can't instantiate abstract class ActivationBase")
})

test_that("a rigid tendon musculotendon collects its states and constants", {
  N <- ReferenceFrame("N")
  O <- PhysPoint("pO"); I <- PhysPoint("pI")
  I$set_pos(O, dynamicsymbols("q") * N$x)
  m <- MusculotendonDeGroote2016("name", LinearPathway(O, I), FirstOrderActivationDeGroote2016("name"))
  expect_equal(sstr(m$x), "Matrix([[a_name(t)]])")
  expect_equal(m$constants$nrow, 31L)
  expect_equal(format(MusculotendonFormulation("TENDON_FORCE_EXPLICIT")), "2")
  expect_error(MusculotendonFormulation(5), "5 is not a valid MusculotendonFormulation")
})
