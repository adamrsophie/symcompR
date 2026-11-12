# The predefined dimensions, prefixes, units, constants and unit systems,
# after sympy/physics/units/definitions/ and sympy/physics/units/systems/.
#
# They are built on first use, in SymPy's import order: the scale factors
# are resolved eagerly when they are set, so the order of the definitions
# matters. physics_units() returns the namespace, with every name of
# sympy.physics.units (aliases included) plus the systems and dictionaries.

#' Predefined physical units
#'
#' Returns the environment that holds the dimensions, prefixes, units,
#' physical constants and unit systems of `sympy.physics.units`, under their
#' SymPy names and aliases (`meter`, `m`, `kilometer`, `speed_of_light`,
#' `c`, `length`, `kilo`, ...), together with the systems `SI`, `MKS`,
#' `MKSA`, `natural` and `cgs_gauss`, the dimension systems `dimsys_SI`,
#' `dimsys_MKSA`, `dimsys_length_weight_time`, `dimsys_cgs` and
#' `dimsys_default`, and the prefix lists `PREFIXES` and `BIN_PREFIXES`.
#' Since several names (`length`, `time`, `c`, `t`, `T`, `F`, `S`, `N`)
#' would mask R or package functions, they are not exported one by one; use
#' `u <- physics_units()` and `u$meter`, or `with(physics_units(), ...)`.
#'
#' @return An environment.
#' @examples
#' u <- physics_units()
#' 5 * u$kilometer
#' with(physics_units(), convert_to(mile, kilometer))
#' @export
physics_units <- function() {
  .units_ensure()
  .units_state$ns
}

.units_ensure <- function() {
  if (!.units_state$built) {
    .units_state$built <- TRUE
    ok <- FALSE
    on.exit(if (!ok) .units_state$built <- FALSE)
    .units_build()
    ok <- TRUE
  }
  invisible(NULL)
}

.units_build <- function() {
  ns <- new.env()
  .units_state$ns <- ns
  Q <- function(name, abbrev = NULL, latex_repr = NULL) Quantity(name, abbrev = abbrev, latex_repr = latex_repr)
  PC <- function(name, abbrev = NULL, latex_repr = NULL) PhysicalConstant(name, abbrev = abbrev, latex_repr = latex_repr)
  def <- function(value, ...) {
    for (nm in c(...)) assign(nm, value, envir = ns)
    value
  }
  rsf <- set_global_relative_scale_factor
  gdim <- set_global_dimension
  One <- S_One
  R <- function(p, q) Rational(p, q)

  # definitions/dimension_definitions.py
  angle <- def(Dimension("angle"), "angle")
  length <- def(Dimension("length", "L"), "length")
  mass <- def(Dimension("mass", "M"), "mass")
  time <- def(Dimension("time", "T"), "time")
  current <- def(Dimension("current", "I"), "current")
  temperature <- def(Dimension("temperature", "T"), "temperature")
  amount_of_substance <- def(Dimension("amount_of_substance"), "amount_of_substance", "amount")
  luminous_intensity <- def(Dimension("luminous_intensity"), "luminous_intensity", "luminosity")
  velocity <- def(Dimension("velocity"), "velocity", "speed")
  acceleration <- def(Dimension("acceleration"), "acceleration")
  momentum <- def(Dimension("momentum"), "momentum")
  force <- def(Dimension("force", "F"), "force")
  energy <- def(Dimension("energy", "E"), "energy")
  power <- def(Dimension("power"), "power")
  pressure <- def(Dimension("pressure"), "pressure")
  frequency <- def(Dimension("frequency", "f"), "frequency")
  action <- def(Dimension("action", "A"), "action")
  area <- def(Dimension("area"), "area")
  volume <- def(Dimension("volume"), "volume")
  voltage <- def(Dimension("voltage", "U"), "voltage")
  impedance <- def(Dimension("impedance", "Z"), "impedance")
  conductance <- def(Dimension("conductance", "G"), "conductance")
  capacitance <- def(Dimension("capacitance"), "capacitance")
  inductance <- def(Dimension("inductance"), "inductance")
  charge <- def(Dimension("charge", "Q"), "charge")
  magnetic_density <- def(Dimension("magnetic_density", "B"), "magnetic_density", "magnetic_flux_density")
  magnetic_flux <- def(Dimension("magnetic_flux"), "magnetic_flux")
  information <- def(Dimension("information"), "information")

  # prefixes.py
  pfx <- list(c("yotta", "Y", 24), c("zetta", "Z", 21), c("exa", "E", 18), c("peta", "P", 15),
              c("tera", "T", 12), c("giga", "G", 9), c("mega", "M", 6), c("kilo", "k", 3),
              c("hecto", "h", 2), c("deca", "da", 1), c("deci", "d", -1), c("centi", "c", -2),
              c("milli", "m", -3), c("micro", "mu", -6), c("nano", "n", -9), c("pico", "p", -12),
              c("femto", "f", -15), c("atto", "a", -18), c("zepto", "z", -21), c("yocto", "y", -24))
  PREFIXES <- list()
  for (p in pfx) {
    pr <- Prefix(p[1L], p[2L], as.integer(p[3L]), latex_repr = if (p[1L] == "micro") "\\mu" else NULL)
    assign(p[1L], pr, envir = ns)
    PREFIXES[[p[2L]]] <- pr
  }
  bpfx <- list(c("kibi", "Ki", 10), c("mebi", "Mi", 20), c("gibi", "Gi", 30), c("tebi", "Ti", 40),
               c("pebi", "Pi", 50), c("exbi", "Ei", 60))
  BIN_PREFIXES <- list()
  for (p in bpfx) {
    pr <- Prefix(p[1L], p[2L], as.integer(p[3L]), 2)
    assign(p[1L], pr, envir = ns)
    BIN_PREFIXES[[p[2L]]] <- pr
  }
  ns$PREFIXES <- PREFIXES
  ns$BIN_PREFIXES <- BIN_PREFIXES
  kilo <- ns$kilo; mega <- ns$mega; milli <- ns$milli; micro <- ns$micro; deci <- ns$deci
  centi <- ns$centi; nano <- ns$nano; pico <- ns$pico
  kibi <- ns$kibi; mebi <- ns$mebi; gibi <- ns$gibi; tebi <- ns$tebi; pebi <- ns$pebi; exbi <- ns$exbi

  # definitions/unit_definitions.py
  percent <- def(Q("percent", latex_repr = "\\%"), "percent", "percents")
  rsf(percent, R(1, 100), One)
  permille <- def(Q("permille"), "permille")
  rsf(permille, R(1, 1000), One)
  radian <- def(Q("radian", abbrev = "rad"), "rad", "radian", "radians")
  gdim(radian, angle)
  degree <- def(Q("degree", abbrev = "deg", latex_repr = "^\\circ"), "deg", "degree", "degrees")
  rsf(degree, S_Pi / 180, radian)
  steradian <- def(Q("steradian", abbrev = "sr"), "sr", "steradian", "steradians")
  angular_mil <- def(Q("angular_mil", abbrev = "mil"), "mil", "angular_mil", "angular_mils")
  meter <- def(Q("meter", abbrev = "m"), "m", "meter", "meters")
  gram <- def(Q("gram", abbrev = "g"), "g", "gram", "grams")
  kilogram <- def(Q("kilogram", abbrev = "kg"), "kg", "kilogram", "kilograms")
  rsf(kilogram, kilo, gram)
  second <- def(Q("second", abbrev = "s"), "s", "second", "seconds")
  ampere <- def(Q("ampere", abbrev = "A"), "A", "ampere", "amperes")
  gdim(ampere, current)
  kelvin <- def(Q("kelvin", abbrev = "K"), "K", "kelvin", "kelvins")
  gdim(kelvin, temperature)
  mole <- def(Q("mole", abbrev = "mol"), "mol", "mole", "moles")
  gdim(mole, amount_of_substance)
  candela <- def(Q("candela", abbrev = "cd"), "cd", "candela", "candelas")
  gdim(candela, luminous_intensity)
  newton <- def(Q("newton", abbrev = "N"), "newton", "newtons", "N")
  kilonewton <- def(Q("kilonewton", abbrev = "kN"), "kilonewton", "kilonewtons", "kN")
  rsf(kilonewton, kilo, newton)
  meganewton <- def(Q("meganewton", abbrev = "MN"), "meganewton", "meganewtons", "MN")
  rsf(meganewton, mega, newton)
  joule <- def(Q("joule", abbrev = "J"), "joule", "joules", "J")
  watt <- def(Q("watt", abbrev = "W"), "watt", "watts", "W")
  pascal <- def(Q("pascal", abbrev = "Pa"), "pascal", "pascals", "Pa", "pa")
  hertz <- def(Q("hertz", abbrev = "Hz"), "hertz", "hz", "Hz")
  dyne <- def(Q("dyne"), "dyne")
  rsf(dyne, One / Integer(10)^5, newton)
  erg <- def(Q("erg"), "erg")
  rsf(erg, One / Integer(10)^7, joule)
  coulomb <- def(Q("coulomb", abbrev = "C"), "coulomb", "coulombs", "C")
  gdim(coulomb, charge)
  volt <- def(Q("volt", abbrev = "V"), "volt", "volts", "v", "V")
  gdim(volt, voltage)
  ohm <- def(Q("ohm", abbrev = "ohm", latex_repr = "\\Omega"), "ohm", "ohms")
  gdim(ohm, impedance)
  siemens <- def(Q("siemens", abbrev = "S"), "siemens", "S", "mho", "mhos")
  gdim(siemens, conductance)
  farad <- def(Q("farad", abbrev = "F"), "farad", "farads", "F")
  gdim(farad, capacitance)
  henry <- def(Q("henry", abbrev = "H"), "henry", "henrys", "H")
  gdim(henry, inductance)
  tesla <- def(Q("tesla", abbrev = "T"), "tesla", "teslas", "T")
  gdim(tesla, magnetic_density)
  weber <- def(Q("weber", abbrev = "Wb"), "weber", "webers", "Wb", "wb")
  gdim(weber, magnetic_flux)
  statampere <- def(Q("statampere"), "statampere")
  statcoulomb <- def(Q("statcoulomb", abbrev = "statC"), "statcoulomb", "statC", "franklin")
  statvolt <- def(Q("statvolt"), "statvolt")
  gauss <- def(Q("gauss"), "gauss")
  maxwell <- def(Q("maxwell"), "maxwell")
  debye <- def(Q("debye"), "debye")
  oersted <- def(Q("oersted"), "oersted")
  dioptre <- def(Q("dioptre"), "optical_power", "dioptre", "diopter", "D")
  lux <- def(Q("lux", abbrev = "lx"), "lux", "lx")
  katal <- def(Q("katal", abbrev = "kat"), "katal", "kat")
  gray <- def(Q("gray"), "gray", "Gy")
  becquerel <- def(Q("becquerel", abbrev = "Bq"), "becquerel", "Bq")
  mg <- def(Q("milligram", abbrev = "mg"), "mg", "milligram", "milligrams")
  rsf(mg, milli, gram)
  ug <- def(Q("microgram", abbrev = "ug", latex_repr = "\\mu\\text{g}"), "ug", "microgram", "micrograms")
  rsf(ug, micro, gram)
  atomic_mass_constant <- def(PC("atomic_mass_constant"), "Da", "dalton", "amu", "amus", "atomic_mass_unit",
                              "atomic_mass_constant")
  tonne <- def(Q("tonne", abbrev = "t"), "t", "metric_ton", "tonne")
  rsf(tonne, mega, gram)
  electron_rest_mass <- def(Q("electron_rest_mass", abbrev = "me"), "me", "electron_rest_mass")
  km <- def(Q("kilometer", abbrev = "km"), "km", "kilometer", "kilometers")
  rsf(km, kilo, meter)
  dm <- def(Q("decimeter", abbrev = "dm"), "dm", "decimeter", "decimeters")
  rsf(dm, deci, meter)
  cm <- def(Q("centimeter", abbrev = "cm"), "cm", "centimeter", "centimeters")
  rsf(cm, centi, meter)
  mm <- def(Q("millimeter", abbrev = "mm"), "mm", "millimeter", "millimeters")
  rsf(mm, milli, meter)
  um <- def(Q("micrometer", abbrev = "um", latex_repr = "\\mu\\text{m}"), "um", "micrometer", "micrometers",
            "micron", "microns")
  rsf(um, micro, meter)
  nm <- def(Q("nanometer", abbrev = "nm"), "nm", "nanometer", "nanometers")
  rsf(nm, nano, meter)
  pm <- def(Q("picometer", abbrev = "pm"), "pm", "picometer", "picometers")
  rsf(pm, pico, meter)
  ft <- def(Q("foot", abbrev = "ft"), "ft", "foot", "feet")
  rsf(ft, R(3048, 10000), meter)
  inch <- def(Q("inch"), "inch", "inches")
  rsf(inch, R(1, 12), ft)
  yd <- def(Q("yard", abbrev = "yd"), "yd", "yard", "yards")
  rsf(yd, 3, ft)
  mi <- def(Q("mile"), "mi", "mile", "miles")
  rsf(mi, 5280, ft)
  nmi <- def(Q("nautical_mile"), "nmi", "nautical_mile", "nautical_miles")
  rsf(nmi, 6076, ft)
  angstrom <- def(Q("angstrom", latex_repr = "\\r{A}"), "angstrom", "angstroms")
  rsf(angstrom, One / Integer(10)^10, meter)
  hectare <- def(Q("hectare", abbrev = "ha"), "ha", "hectare")
  liter <- def(Q("liter", abbrev = "l"), "l", "L", "liter", "liters")
  dl <- def(Q("deciliter", abbrev = "dl"), "dl", "dL", "deciliter", "deciliters")
  rsf(dl, R(1, 10), liter)
  cl <- def(Q("centiliter", abbrev = "cl"), "cl", "cL", "centiliter", "centiliters")
  rsf(cl, R(1, 100), liter)
  ml <- def(Q("milliliter", abbrev = "ml"), "ml", "mL", "milliliter", "milliliters")
  rsf(ml, R(1, 1000), liter)
  ms <- def(Q("millisecond", abbrev = "ms"), "ms", "millisecond", "milliseconds")
  rsf(ms, milli, second)
  us <- def(Q("microsecond", abbrev = "us", latex_repr = "\\mu\\text{s}"), "us", "microsecond", "microseconds")
  rsf(us, micro, second)
  ns_ <- def(Q("nanosecond", abbrev = "ns"), "ns", "nanosecond", "nanoseconds")
  rsf(ns_, nano, second)
  ps <- def(Q("picosecond", abbrev = "ps"), "ps", "picosecond", "picoseconds")
  rsf(ps, pico, second)
  minute <- def(Q("minute"), "minute", "minutes")
  rsf(minute, 60, second)
  hour <- def(Q("hour"), "h", "hour", "hours")
  rsf(hour, 60, minute)
  day <- def(Q("day"), "day", "days")
  rsf(day, 24, hour)
  y1 <- def(Q("anomalistic_year"), "anomalistic_year", "anomalistic_years")
  rsf(y1, 365.259636, day)
  y2 <- def(Q("sidereal_year"), "sidereal_year", "sidereal_years")
  rsf(y2, 31558149.540, second)
  tropical_year <- def(Q("tropical_year"), "tropical_year", "tropical_years", "year", "years")
  rsf(tropical_year, 365.24219, day)
  y3 <- def(Q("common_year"), "common_year", "common_years")
  rsf(y3, 365, day)
  julian_year <- def(Q("julian_year"), "julian_year", "julian_years")
  rsf(julian_year, 365 + One / 4, day)
  y4 <- def(Q("draconic_year"), "draconic_year", "draconic_years")
  rsf(y4, 346.62, day)
  y5 <- def(Q("gaussian_year"), "gaussian_year", "gaussian_years")
  rsf(y5, 365.2568983, day)
  y6 <- def(Q("full_moon_cycle"), "full_moon_cycle", "full_moon_cycles")
  rsf(y6, 411.78443029, day)

  G <- def(PC("gravitational_constant", abbrev = "G"), "G", "gravitational_constant")
  c <- def(PC("speed_of_light", abbrev = "c"), "c", "speed_of_light")
  elementary_charge <- def(PC("elementary_charge", abbrev = "e"), "elementary_charge")
  planck <- def(PC("planck", abbrev = "h"), "planck")
  hbar <- def(PC("hbar", abbrev = "hbar"), "hbar")
  eV <- def(PC("electronvolt", abbrev = "eV"), "eV", "electronvolt", "electronvolts")
  avogadro_number <- def(PC("avogadro_number"), "avogadro_number")
  avogadro_constant <- def(PC("avogadro_constant"), "avogadro", "avogadro_constant")
  boltzmann <- def(PC("boltzmann_constant"), "boltzmann", "boltzmann_constant")
  stefan <- def(PC("stefan_boltzmann_constant"), "stefan", "stefan_boltzmann_constant")
  molar_gas_constant <- def(PC("molar_gas_constant", abbrev = "R"), "R", "molar_gas_constant")
  faraday_constant <- def(PC("faraday_constant"), "faraday_constant")
  josephson_constant <- def(PC("josephson_constant", abbrev = "K_j"), "josephson_constant")
  von_klitzing_constant <- def(PC("von_klitzing_constant", abbrev = "R_k"), "von_klitzing_constant")
  gee <- def(PC("acceleration_due_to_gravity", abbrev = "g"), "gee", "gees", "acceleration_due_to_gravity")
  u0 <- def(PC("magnetic_constant"), "u0", "magnetic_constant", "vacuum_permeability")
  e0 <- def(PC("vacuum_permittivity"), "e0", "electric_constant", "vacuum_permittivity")
  Z0 <- def(PC("vacuum_impedance", abbrev = "Z_0", latex_repr = "Z_{0}"), "Z0", "vacuum_impedance")
  coulomb_constant <- def(PC("coulomb_constant", abbrev = "k_e"), "coulomb_constant", "coulombs_constant",
                          "electric_force_constant")
  atmosphere <- def(Q("atmosphere", abbrev = "atm"), "atmosphere", "atmospheres", "atm")
  kPa <- def(Q("kilopascal", abbrev = "kPa"), "kPa", "kilopascal")
  rsf(kPa, kilo, pascal)
  bar <- def(Q("bar", abbrev = "bar"), "bar", "bars")
  pound <- def(Q("pound"), "pound", "pounds")
  psi <- def(Q("psi"), "psi")
  dHg0 <- def(Float(13.5951), "dHg0")
  mmHg <- def(Q("mmHg"), "mmHg", "torr")
  rsf(atmosphere, 101325, pascal)
  rsf(bar, 100, kPa)
  rsf(pound, R(45359237, 100000000), kilogram)
  mmu <- def(Q("milli_mass_unit"), "mmu", "mmus", "milli_mass_unit")
  quart <- def(Q("quart"), "quart", "quarts")
  ly <- def(Q("lightyear", abbrev = "ly"), "ly", "lightyear", "lightyears")
  au <- def(Q("astronomical_unit", abbrev = "AU"), "au", "astronomical_unit", "astronomical_units")
  planck_mass <- def(Q("planck_mass", abbrev = "m_P", latex_repr = "m_\\text{P}"), "planck_mass")
  planck_time <- def(Q("planck_time", abbrev = "t_P", latex_repr = "t_\\text{P}"), "planck_time")
  planck_temperature <- def(Q("planck_temperature", abbrev = "T_P", latex_repr = "T_\\text{P}"), "planck_temperature")
  planck_length <- def(Q("planck_length", abbrev = "l_P", latex_repr = "l_\\text{P}"), "planck_length")
  planck_charge <- def(Q("planck_charge", abbrev = "q_P", latex_repr = "q_\\text{P}"), "planck_charge")
  planck_area <- def(Q("planck_area"), "planck_area")
  planck_volume <- def(Q("planck_volume"), "planck_volume")
  planck_momentum <- def(Q("planck_momentum"), "planck_momentum")
  planck_energy <- def(Q("planck_energy", abbrev = "E_P", latex_repr = "E_\\text{P}"), "planck_energy")
  planck_force <- def(Q("planck_force", abbrev = "F_P", latex_repr = "F_\\text{P}"), "planck_force")
  planck_power <- def(Q("planck_power", abbrev = "P_P", latex_repr = "P_\\text{P}"), "planck_power")
  planck_density <- def(Q("planck_density", abbrev = "rho_P", latex_repr = "\\rho_\\text{P}"), "planck_density")
  planck_energy_density <- def(Q("planck_energy_density", abbrev = "rho^E_P"), "planck_energy_density")
  planck_intensity <- def(Q("planck_intensity", abbrev = "I_P", latex_repr = "I_\\text{P}"), "planck_intensity")
  planck_angular_frequency <- def(Q("planck_angular_frequency", abbrev = "omega_P", latex_repr = "\\omega_\\text{P}"),
                                  "planck_angular_frequency")
  planck_pressure <- def(Q("planck_pressure", abbrev = "p_P", latex_repr = "p_\\text{P}"), "planck_pressure")
  planck_current <- def(Q("planck_current", abbrev = "I_P", latex_repr = "I_\\text{P}"), "planck_current")
  planck_voltage <- def(Q("planck_voltage", abbrev = "V_P", latex_repr = "V_\\text{P}"), "planck_voltage")
  planck_impedance <- def(Q("planck_impedance", abbrev = "Z_P", latex_repr = "Z_\\text{P}"), "planck_impedance")
  planck_acceleration <- def(Q("planck_acceleration", abbrev = "a_P", latex_repr = "a_\\text{P}"), "planck_acceleration")
  bit <- def(Q("bit"), "bit", "bits")
  gdim(bit, information)
  byte <- def(Q("byte"), "byte", "bytes")
  kibibyte <- def(Q("kibibyte"), "kibibyte", "kibibytes")
  mebibyte <- def(Q("mebibyte"), "mebibyte", "mebibytes")
  gibibyte <- def(Q("gibibyte"), "gibibyte", "gibibytes")
  tebibyte <- def(Q("tebibyte"), "tebibyte", "tebibytes")
  pebibyte <- def(Q("pebibyte"), "pebibyte", "pebibytes")
  exbibyte <- def(Q("exbibyte"), "exbibyte", "exbibytes")
  rsf(byte, 8, bit)
  rsf(kibibyte, kibi, byte)
  rsf(mebibyte, mebi, byte)
  rsf(gibibyte, gibi, byte)
  rsf(tebibyte, tebi, byte)
  rsf(pebibyte, pebi, byte)
  rsf(exbibyte, exbi, byte)
  curie <- def(Q("curie", abbrev = "Ci"), "curie", "Ci")
  rutherford <- def(Q("rutherford", abbrev = "Rd"), "rutherford", "Rd")

  m <- meter; s <- second; kg <- kilogram

  # systems/length_weight_time.py
  lwt <- DimensionSystem(list(length, mass, time), dimensional_dependencies = list(
    velocity = list(length = 1, time = -1),
    acceleration = list(length = 1, time = -2),
    momentum = list(mass = 1, length = 1, time = -1),
    force = list(mass = 1, length = 1, time = -2),
    energy = list(mass = 1, length = 2, time = -2),
    power = list(length = 2, mass = 1, time = -3),
    pressure = list(mass = 1, length = -1, time = -2),
    frequency = list(time = -1),
    action = list(length = 2, mass = 1, time = -1),
    area = list(length = 2),
    volume = list(length = 3)))
  ns$dimsys_length_weight_time <- lwt
  sqd <- set_quantity_dimension
  sqs <- set_quantity_scale_factor
  sqd(lwt, meter, length); sqs(lwt, meter, One)
  sqd(lwt, gram, mass); sqs(lwt, gram, One)
  sqd(lwt, second, time); sqs(lwt, second, One)
  sqd(lwt, newton, force); sqs(lwt, newton, kilogram * meter / second^2)
  sqd(lwt, joule, energy); sqs(lwt, joule, newton * meter)
  sqd(lwt, watt, power); sqs(lwt, watt, joule / second)
  sqd(lwt, pascal, pressure); sqs(lwt, pascal, newton / meter^2)
  sqd(lwt, hertz, frequency); sqs(lwt, hertz, One)
  sqd(lwt, dioptre, 1 / length); sqs(lwt, dioptre, 1 / meter)
  sqd(lwt, hectare, length^2); sqs(lwt, hectare, (meter^2) * 10000)
  sqd(lwt, liter, length^3); sqs(lwt, liter, meter^3 / 1000)
  sqd(lwt, G, length^3 * mass^-1 * time^-2); sqs(lwt, G, 6.67430e-11 * m^3 / (kg * s^2))
  sqd(lwt, c, velocity); sqs(lwt, c, 299792458 * meter / second)
  sqd(lwt, planck, action); sqs(lwt, planck, 6.62607015e-34 * joule * second)
  sqd(lwt, hbar, action); sqs(lwt, hbar, planck / (2 * S_Pi))

  # systems/mks.py
  all_units <- list(gram, joule, newton, watt, pascal, hertz)
  for (u in list(meter, gram, second, joule, newton, watt, pascal, hertz))
    all_units <- c(all_units, prefix_unit(u, PREFIXES))
  all_units <- c(all_units, list(G, c))
  MKS <- UnitSystem(list(meter, kilogram, second), all_units, name = "MKS", dimension_system = lwt,
                    derived_units = list(list(power, watt), list(time, second), list(pressure, pascal),
                                         list(length, meter), list(frequency, hertz), list(mass, kilogram),
                                         list(force, newton), list(energy, joule), list(velocity, meter / second),
                                         list(acceleration, meter / (second^2))))
  ns$MKS <- MKS

  # systems/mksa.py
  mksa_units <- list(ampere, volt, ohm, siemens, farad, henry, coulomb, tesla, weber)
  all_units <- list()
  for (u in mksa_units) all_units <- c(all_units, prefix_unit(u, PREFIXES))
  all_units <- c(all_units, mksa_units, list(Z0))
  dimsys_MKSA <- extend(lwt, list(current), new_dim_deps = list(
    voltage = list(mass = 1, length = 2, current = -1, time = -3),
    impedance = list(mass = 1, length = 2, current = -2, time = -3),
    conductance = list(mass = -1, length = -2, current = 2, time = 3),
    capacitance = list(mass = -1, length = -2, current = 2, time = 4),
    inductance = list(mass = 1, length = 2, current = -2, time = -2),
    charge = list(current = 1, time = 1),
    magnetic_density = list(mass = 1, current = -1, time = -2),
    magnetic_flux = list(length = 2, mass = 1, current = -1, time = -2)))
  ns$dimsys_MKSA <- dimsys_MKSA
  MKSA <- extend(MKS, list(ampere), all_units, name = "MKSA", dimension_system = dimsys_MKSA,
                 derived_units = list(list(magnetic_flux, weber), list(impedance, ohm), list(current, ampere),
                                      list(voltage, volt), list(inductance, henry), list(conductance, siemens),
                                      list(magnetic_density, tesla), list(charge, coulomb), list(capacitance, farad)))
  ns$MKSA <- MKSA

  # systems/natural.py
  ns$natural <- UnitSystem(list(hbar, eV, c), prefix_unit(eV, PREFIXES), name = "Natural system")

  # systems/si.py
  si_units <- list(mole, candela, kelvin, lux, hertz, newton, pascal, joule, watt, coulomb, volt, farad, ohm,
                   siemens, weber, tesla, henry, candela, lux, becquerel, gray, katal)
  all_units <- list()
  for (u in si_units) all_units <- c(all_units, prefix_unit(u, PREFIXES))
  all_units <- c(all_units, si_units, list(mole, candela, kelvin, lux))
  dimsys_SI <- extend(dimsys_MKSA, list(temperature, amount_of_substance, luminous_intensity))
  ns$dimsys_SI <- dimsys_SI
  ns$dimsys_default <- extend(dimsys_SI, list(information))
  SI <- extend(MKSA, list(mole, candela, kelvin), all_units, name = "SI", dimension_system = dimsys_SI,
               derived_units = list(list(power, watt), list(magnetic_flux, weber), list(time, second),
                                    list(impedance, ohm), list(pressure, pascal), list(current, ampere),
                                    list(voltage, volt), list(length, meter), list(frequency, hertz),
                                    list(inductance, henry), list(temperature, kelvin),
                                    list(amount_of_substance, mole), list(luminous_intensity, candela),
                                    list(conductance, siemens), list(mass, kilogram), list(magnetic_density, tesla),
                                    list(charge, coulomb), list(force, newton), list(capacitance, farad),
                                    list(energy, joule), list(velocity, meter / second)))
  ns$SI <- SI
  sqd(SI, radian, One)
  sqs(SI, ampere, One)
  sqs(SI, kelvin, One)
  sqs(SI, mole, One)
  sqs(SI, candela, One)
  sqs(SI, coulomb, One)
  sqs(SI, volt, joule / coulomb)
  sqs(SI, ohm, volt / ampere)
  sqs(SI, siemens, ampere / volt)
  sqs(SI, farad, coulomb / volt)
  sqs(SI, henry, volt * second / ampere)
  sqs(SI, tesla, volt * second / meter^2)
  sqs(SI, weber, joule / ampere)
  sqd(SI, lux, luminous_intensity / length^2); sqs(SI, lux, steradian * candela / meter^2)
  sqd(SI, katal, amount_of_substance / time); sqs(SI, katal, mole / second)
  sqd(SI, gray, energy / mass); sqs(SI, gray, meter^2 / second^2)
  sqd(SI, becquerel, 1 / time); sqs(SI, becquerel, 1 / second)
  sqd(SI, elementary_charge, charge); sqs(SI, elementary_charge, 1.602176634e-19 * coulomb)
  sqd(SI, eV, energy); sqs(SI, eV, 1.602176634e-19 * joule)
  sqd(SI, avogadro_number, One); sqs(SI, avogadro_number, Float(6.02214076e23))
  sqd(SI, avogadro_constant, amount_of_substance^-1); sqs(SI, avogadro_constant, avogadro_number / mole)
  sqd(SI, boltzmann, energy / temperature); sqs(SI, boltzmann, 1.380649e-23 * joule / kelvin)
  sqd(SI, stefan, energy * time^-1 * length^-2 * temperature^-4)
  sqs(SI, stefan, S_Pi^2 * boltzmann^4 / (60 * hbar^3 * c^2))
  sqd(SI, atomic_mass_constant, mass); sqs(SI, atomic_mass_constant, 1.66053906660e-24 * gram)
  sqd(SI, molar_gas_constant, energy / (temperature * amount_of_substance))
  sqs(SI, molar_gas_constant, boltzmann * avogadro_constant)
  sqd(SI, faraday_constant, charge / amount_of_substance)
  sqs(SI, faraday_constant, elementary_charge * avogadro_constant)
  sqd(SI, josephson_constant, frequency / voltage); sqs(SI, josephson_constant, 2 * elementary_charge / planck)
  sqd(SI, von_klitzing_constant, voltage / current)
  sqs(SI, von_klitzing_constant, planck / elementary_charge^2)
  sqd(SI, gee, acceleration); sqs(SI, gee, 9.80665 * meter / second^2)
  sqd(SI, u0, force / current^2); sqs(SI, u0, 4 * S_Pi / Integer(10)^7 * newton / ampere^2)
  sqd(SI, e0, capacitance / length); sqs(SI, e0, 1 / (u0 * c^2))
  sqd(SI, Z0, impedance); sqs(SI, Z0, u0 * c)
  sqd(SI, electron_rest_mass, mass); sqs(SI, electron_rest_mass, 9.1093837015e-31 * kilogram)
  sqd(SI, coulomb_constant, force * length^2 / charge^2); sqs(SI, coulomb_constant, 1 / (4 * S_Pi * e0))
  sqd(SI, psi, pressure); sqs(SI, psi, pound * gee / inch^2)
  sqd(SI, mmHg, pressure); sqs(SI, mmHg, dHg0 * gee * kilogram / meter^2)
  sqd(SI, mmu, mass); sqs(SI, mmu, atomic_mass_constant / 1000)
  sqd(SI, quart, length^3); sqs(SI, quart, R(231, 4) * inch^3)
  sqd(SI, ly, length); sqs(SI, ly, c * julian_year)
  sqd(SI, au, length); sqs(SI, au, 149597870691 * meter)
  sqd(SI, planck_mass, mass); sqs(SI, planck_mass, sqrt(hbar * c / G))
  sqd(SI, planck_time, time); sqs(SI, planck_time, sqrt(hbar * G / c^5))
  sqd(SI, planck_temperature, temperature); sqs(SI, planck_temperature, sqrt(hbar * c^5 / G / boltzmann^2))
  sqd(SI, planck_length, length); sqs(SI, planck_length, sqrt(hbar * G / c^3))
  sqd(SI, planck_charge, charge); sqs(SI, planck_charge, sqrt(4 * S_Pi * e0 * hbar * c))
  sqd(SI, planck_area, length^2); sqs(SI, planck_area, planck_length^2)
  sqd(SI, planck_volume, length^3); sqs(SI, planck_volume, planck_length^3)
  sqd(SI, planck_momentum, mass * velocity); sqs(SI, planck_momentum, planck_mass * c)
  sqd(SI, planck_energy, energy); sqs(SI, planck_energy, planck_mass * c^2)
  sqd(SI, planck_force, force); sqs(SI, planck_force, planck_energy / planck_length)
  sqd(SI, planck_power, power); sqs(SI, planck_power, planck_energy / planck_time)
  sqd(SI, planck_density, mass / length^3); sqs(SI, planck_density, planck_mass / planck_length^3)
  sqd(SI, planck_energy_density, energy / length^3)
  sqs(SI, planck_energy_density, planck_energy / planck_length^3)
  sqd(SI, planck_intensity, mass * time^(-3)); sqs(SI, planck_intensity, planck_energy_density * c)
  sqd(SI, planck_angular_frequency, 1 / time); sqs(SI, planck_angular_frequency, 1 / planck_time)
  sqd(SI, planck_pressure, pressure); sqs(SI, planck_pressure, planck_force / planck_length^2)
  sqd(SI, planck_current, current); sqs(SI, planck_current, planck_charge / planck_time)
  sqd(SI, planck_voltage, voltage); sqs(SI, planck_voltage, planck_energy / planck_charge)
  sqd(SI, planck_impedance, impedance); sqs(SI, planck_impedance, planck_voltage / planck_current)
  sqd(SI, planck_acceleration, acceleration); sqs(SI, planck_acceleration, c / planck_time)
  sqd(SI, curie, 1 / time); sqs(SI, curie, 37000000000 * becquerel)
  sqd(SI, rutherford, 1 / time); sqs(SI, rutherford, 1000000 * becquerel)

  # systems/cgs.py
  dimsys_cgs <- extend(lwt, list(), new_dim_deps = list(
    impedance = list(time = 1, length = -1),
    conductance = list(time = -1, length = 1),
    capacitance = list(length = 1),
    inductance = list(time = 2, length = -1),
    charge = list(mass = S_Half, length = R(3, 2), time = -1),
    current = list(mass = One / 2, length = 3 * One / 2, time = -2),
    voltage = list(length = -One / 2, mass = One / 2, time = -1),
    magnetic_density = list(length = -One / 2, mass = One / 2, time = -1),
    magnetic_flux = list(length = 3 * One / 2, mass = One / 2, time = -1)))
  ns$dimsys_cgs <- dimsys_cgs
  cgs <- UnitSystem(list(cm, gram, second), list(), name = "cgs_gauss", dimension_system = dimsys_cgs)
  ns$cgs_gauss <- cgs
  sqs(cgs, coulomb_constant, 1)
  sqd(cgs, statcoulomb, charge)
  sqs(cgs, statcoulomb, cm^R(3, 2) * gram^S_Half / second)
  sqd(cgs, coulomb, charge)
  sqd(cgs, statampere, current)
  sqs(cgs, statampere, statcoulomb / second)
  sqd(cgs, statvolt, voltage)
  sqs(cgs, statvolt, erg / statcoulomb)
  sqd(cgs, volt, voltage)
  sqd(cgs, gauss, magnetic_density)
  sqs(cgs, gauss, sqrt(gram / cm) / second)
  sqd(cgs, tesla, magnetic_density)
  sqd(cgs, maxwell, magnetic_flux)
  sqs(cgs, maxwell, sqrt(cm^3 * gram) / second)
  sqs(cgs, coulomb, 10 * c * statcoulomb)
  sqs(cgs, ampere, 10 * c * statcoulomb / second)
  sqs(cgs, volt, Integer(10)^6 / c * statvolt)
  sqs(cgs, weber, Integer(10)^8 * maxwell)
  sqs(cgs, tesla, Integer(10)^4 * gauss)
  sqs(cgs, debye, One / Integer(10)^18 * statcoulomb * cm)
  sqs(cgs, oersted, sqrt(gram / cm) / second)
  sqs(cgs, ohm, Integer(10)^5 / c^2 * second / cm)
  sqs(cgs, farad, One / Integer(10)^5 * c^2 * cm)
  sqs(cgs, henry, Integer(10)^5 / c^2 / cm * second^2)
  sqd(cgs, coulomb_constant, 1)
  sqs(cgs, coulomb_constant, 1)
  invisible(NULL)
}
