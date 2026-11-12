# Expected values from sympy/physics/units/tests.
s <- function(e) if (is_sym(e)) sstr(e) else format(e)

test_that("quantities print, convert and simplify as in SymPy", {
  u <- physics_units()
  x <- Symbol("x")
  expect_equal(s(5 * u$kilometer), "5*kilometer")
  expect_equal(s(u$speed_of_light * u$meter * x), "speed_of_light*meter*x")
  expect_equal(s(convert_to(u$mile, u$kilometer)), "25146*kilometer/15625")
  expect_equal(s(convert_to(3 * u$newton, list(u$centimeter, u$gram, u$second))), "300000*centimeter*gram/second^2")
  expect_equal(s(convert_to(u$joule * u$second, u$joule)), "joule*second")
  expect_equal(s(convert_to(u$statcoulomb, u$coulomb, u$cgs_gauss)), "coulomb/2997924580")
  expect_equal(s(quantity_simplify(u$kilo * u$foot * u$inch)), "250*foot^2/3")
  expect_equal(s(quantity_simplify(5 * u$joule / u$coulomb, across_dimensions = TRUE, unit_system = "SI")), "5*volt")
  expect_equal(s(simplify(10 * x * u$K * u$km^2 / u$m / u$cm)), "1000000000*kelvin*x")
  expect_equal(find_unit("charge"), c("C", "coulomb", "coulombs", "planck_charge", "elementary_charge"))
})

test_that("dimensions, prefixes and unit systems follow SymPy", {
  u <- physics_units()
  expect_equal(s(u$length / u$time), "Dimension(length/time)")
  expect_equal(s(2 * u$length), "Dimension(length, L)")
  expect_equal(s(get_dimensional_dependencies(u$dimsys_SI, u$length / u$time)),
               "{Dimension(length, L): 1, Dimension(time, T): -1}")
  expect_error(u$length + u$foot, "cannot sum dimension and quantity", class = "TypeError")
  expect_equal(s(u$kilo * u$kilo), "M")
  expect_equal(s(u$kilo * u$meter), "1000*meter")
  r <- collect_factor_and_dimension(u$SI, u$meter / u$second + 36 * u$km / (10 * u$hour))
  expect_equal(c(s(r[[1L]]), s(r[[2L]])), c("2", "Dimension(length/time)"))
  expect_error(check_dimensions(u$meter + u$second), class = "ValueError")
})
