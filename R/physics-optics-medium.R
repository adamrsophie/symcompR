# Optical media and transverse waves, after sympy/physics/optics/medium.py
# and waves.py.

.optics_units <- function() {
  if (is.null(.optics_state$c)) {
    U <- physics_units()
    .optics_state$c <- convert_to(U$speed_of_light, U$meter / U$second)
    .optics_state$e0 <- convert_to(U$e0, U$ampere^2 * U$second^4 / (U$kilogram * U$meter^3))
    .optics_state$u0 <- convert_to(U$u0, U$meter * U$kilogram / (U$ampere^2 * U$second^2))
  }
  .optics_state
}
.optics_state <- new.env()
.optics_c <- function() .optics_units()$c

#' An optical medium
#'
#' `Medium(name, permittivity, permeability, n)` is a medium given by any two
#' of its permittivity, permeability and refractive index `n`, or the vacuum
#' when none is given. Its properties are `$speed`, `$refractive_index` and
#' `$n`, and for a medium given by its permittivity and permeability also
#' `$permittivity`, `$permeability` and `$intrinsic_impedance`. After
#' `sympy.physics.optics.medium`.
#'
#' @param name The name of the medium.
#' @param permittivity,permeability The electric permittivity and magnetic
#'   permeability, with units.
#' @param n The refractive index.
#' @return A `Medium` object (of class `MediumPP` or `MediumN`).
#' @examples
#' m <- Medium("m", n = 2)
#' m$speed
#' Medium("vacuum")$refractive_index
#' @export
Medium <- function(name, permittivity = NULL, permeability = NULL, n = NULL) {
  if (is_sym(name) && identical(name$type, "Str")) name <- name$name
  if (!is.null(permittivity)) permittivity <- as_sym(permittivity)
  if (!is.null(permeability)) permeability <- as_sym(permeability)
  if (!is.null(n)) n <- as_sym(n)
  c <- .optics_c()
  if (!is.null(n)) {
    if (!is.null(permittivity) && is.null(permeability)) {
      permeability <- n^2 / (c^2 * permittivity)
      return(.medium_pp(name, permittivity, permeability))
    } else if (!is.null(permeability) && is.null(permittivity)) {
      permittivity <- n^2 / (c^2 * permeability)
      return(.medium_pp(name, permittivity, permeability))
    } else if (!is.null(permittivity)) {
      .physics_error("ValueError", "Specifying all of permittivity, permeability, and n is not allowed")
    }
    return(.pyobj_new("MediumN", list(.Str(name), n), name = name))
  }
  if (!is.null(permittivity) && !is.null(permeability)) return(.medium_pp(name, permittivity, permeability))
  if (is.null(permittivity) && is.null(permeability)) {
    u <- .optics_units()
    return(.medium_pp(name, u$e0, u$u0))
  }
  .physics_error("ValueError", "Arguments are underspecified. Either specify n or any two of permittivity, permeability, and n")
}

.medium_pp <- function(name, permittivity, permeability)
  .pyobj_new("MediumPP", list(.Str(name), permittivity, permeability), name = name)

.pyclass_def("Medium",
  props = list(
    speed = function(self) .optics_c() / self$n,
    refractive_index = function(self) .optics_c() / self$speed))
.pyclass_def("MediumN", parent = "Medium",
  ctor = function(name, n) .pyobj_new("MediumN", list(name, n), name = name$name),
  props = list(n = function(self) self$args[[2L]]))
.pyclass_def("MediumPP", parent = "Medium",
  ctor = function(name, permittivity, permeability) .medium_pp(name$name, permittivity, permeability),
  props = list(
    intrinsic_impedance = function(self) sqrt(self$permeability / self$permittivity),
    permittivity = function(self) self$args[[2L]],
    permeability = function(self) self$args[[3L]],
    n = function(self) .optics_c() * sqrt(self$permittivity * self$permeability)))

#' A transverse wave
#'
#' `TWave(amplitude, frequency, phase, time_period, n)` is a simple
#' transverse sinusoidal wave in a medium of refractive index `n`, given by
#' its frequency or its time period. Its properties are `$amplitude`,
#' `$frequency`, `$phase`, `$time_period`, `$n`, `$wavelength`, `$speed`,
#' `$angular_velocity` and `$wavenumber`. Waves of equal frequency add by
#' interference; a wave times a number scales its amplitude.
#' `w$rewrite(target)` writes the wave with `"sin"`, `"cos"` or `"exp"`, or
#' gives the wave equation with `"pde"`. After `sympy.physics.optics.waves`.
#'
#' @param amplitude The amplitude.
#' @param frequency The frequency.
#' @param phase The phase angle.
#' @param time_period The time period.
#' @param n The refractive index of the medium.
#' @return A `TWave` object.
#' @examples
#' f <- Symbol("f")
#' w1 <- TWave(Symbol("A1"), f, Symbol("phi1"))
#' w2 <- TWave(Symbol("A2"), f, Symbol("phi2"))
#' (w1 + w2)$amplitude
#' w1$rewrite("cos")
#' @export
TWave <- function(amplitude, frequency = NULL, phase = S_Zero, time_period = NULL, n = Symbol("n")) {
  if (!is.null(time_period)) {
    time_period <- as_sym(time_period)
    .frequency <- S_One / time_period
  }
  if (!is.null(frequency)) {
    frequency <- as_sym(frequency)
    .time_period <- S_One / frequency
    if (!is.null(time_period) && !identical(frequency$key, (S_One / time_period)$key))
      .physics_error("ValueError", "frequency and time_period should be consistent.")
  }
  if (is.null(frequency) && is.null(time_period)) .physics_error("ValueError", "Either frequency or time period is needed.")
  if (is.null(frequency)) frequency <- .frequency
  if (is.null(time_period)) time_period <- .time_period
  .pyobj_new("TWave", list(as_sym(amplitude), frequency, as_sym(phase), time_period, as_sym(n)))
}

local({
  ta <- function(i) function(self) self$args[[i]]
  tw_add <- function(self, other) {
    if (!inherits(other, "TWave"))
      .physics_error("TypeError", paste0(.py_typename(other), " and TWave objects cannot be added."))
    if (identical(self$frequency$key, other$frequency$key) && identical(self$wavelength$key, other$wavelength$key)) {
      a1 <- self$amplitude; a2 <- other$amplitude; p1 <- self$phase; p2 <- other$phase
      return(TWave(sqrt(a1^2 + a2^2 + 2 * a1 * a2 * cos(p1 - p2)), self$frequency,
                   sym_atan2(a1 * sin(p1) + a2 * sin(p2), a1 * cos(p1) + a2 * cos(p2))))
    }
    .physics_error("NotImplementedError", "Interference of waves with different frequencies has not been implemented.")
  }
  tw_mul <- function(self, other) {
    other <- as_sym(other)
    if (other$type %in% c("Integer", "Rational", "Float"))
      return(do.call(TWave, c(list(self$amplitude * other), self$args[-1L])))
    .physics_error("TypeError", paste0(.py_typename(other), " and TWave objects cannot be multiplied."))
  }
  tw_neg <- function(self) tw_mul(self, -1L)
  phase_arg <- function(self) self$wavenumber * Symbol("x") - self$angular_velocity * Symbol("t") + self$phase
  .pyclass_def("TWave", expr = TRUE,
    ctor = function(amplitude, frequency, phase, time_period, n) TWave(amplitude, frequency, phase, time_period, n),
    props = list(
      amplitude = ta(1L), frequency = ta(2L), phase = ta(3L), time_period = ta(4L), n = ta(5L),
      wavelength = function(self) .optics_c() / (self$frequency * self$n),
      speed = function(self) self$wavelength * self$frequency,
      angular_velocity = function(self) 2 * S_Pi * self$frequency,
      wavenumber = function(self) 2 * S_Pi / self$wavelength),
    methods = list(
      rewrite = function(self, target, ...) {
        switch(target,
          sin = self$amplitude * new_sym("Function", list(.add(list(phase_arg(self), S_Pi / 2))), name = "sin"),
          cos = self$amplitude * cos(phase_arg(self)),
          exp = self$amplitude * exp(S_I * phase_arg(self)),
          pde = {
            x <- Symbol("x"); t <- Symbol("t")
            E <- Function("E")
            Derivative(E(x, t), x, 2) + Symbol("mu") * Symbol("epsilon") * Derivative(E(x, t), t, 2)
          },
          self)
      }),
    ops = list(
      `__add__` = tw_add, `__radd__` = tw_add,
      `__sub__` = function(self, other) tw_add(self, if (inherits(other, "TWave")) tw_neg(other) else -1L * as_sym(other)),
      `__rsub__` = function(self, other) tw_add(tw_neg(self), other),
      `__mul__` = tw_mul, `__rmul__` = tw_mul,
      `__neg__` = tw_neg))
})
