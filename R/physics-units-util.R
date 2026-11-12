# Conversion and simplification of expressions with units, after
# sympy/physics/units/util.py, and find_unit() from
# sympy/physics/units/__init__.py.

# the target units of convert_to(): one unit, an R list or a sym Tuple
.u_targets <- function(target_units) {
  if (is_sym(target_units)) return(if (target_units$type == "Tuple") target_units$args else list(target_units))
  if (is.list(target_units)) return(lapply(target_units, as_sym))
  if (length(target_units) > 1L) return(lapply(as.list(target_units), as_sym))
  list(as_sym(target_units))
}

.get_conversion_matrix_for_expr <- function(expr, target_units, unit_system) {
  expr_dim <- Dimension(get_dimensional_expr(unit_system, expr))
  ds <- unit_system$dimension_system
  if (is.null(ds)) .u_stop("AttributeError", "'NoneType' object has no attribute 'get_dimensional_dependencies'")
  dim_dependencies <- get_dimensional_dependencies(ds, expr_dim, mark_dimensionless = TRUE)
  target_dims <- lapply(target_units, function(x) Dimension(get_dimensional_expr(unit_system, x)))
  tdeps <- lapply(target_dims, function(x) get_dimensional_dependencies(ds, x, mark_dimensionless = TRUE))
  canon_dim_units <- do.call(c, lapply(tdeps, `[[`, "k"))
  if (!all(vapply(dim_dependencies$k, function(k) .u_in(k, canon_dim_units), TRUE))) return(NULL)
  canon_dim_units <- canon_dim_units[!duplicated(vapply(canon_dim_units, .k, ""))]
  rows <- lapply(canon_dim_units, function(j)
    c(lapply(tdeps, function(d) .ud_get(d, j, S_Zero)), list(.ud_get(dim_dependencies, j, S_Zero))))
  sol <- .gauss_solve(rows, length(target_dims))
  # Matrix.solve raises NonInvertibleMatrixError unless the solution is unique
  if (is.null(sol) || length(attr(sol, "pivots")) < length(target_dims)) return(NULL)
  attr(sol, "pivots") <- NULL
  sol
}

#' Convert and simplify expressions with units
#'
#' `convert_to(expr, target_units, unit_system)` rewrites the units and
#' quantities of `expr` as products of powers of `target_units` (one unit, or
#' a list of units) whenever the dimensions allow it, as SymPy's
#' `convert_to()`. `quantity_simplify()` replaces prefixes by their values and
#' units of the same dimension by one of them; with `across_dimensions = TRUE`
#' it also converts to the preferred unit of the resulting dimension in
#' `unit_system`. `check_dimensions()` returns `expr` when the terms of every
#' sum agree in dimension, and raises a `ValueError` otherwise.
#' `find_unit()` lists the names of the predefined units that contain a
#' string, or that share the dimension of a unit, dimension or expression.
#'
#' @param expr An expression.
#' @param target_units A unit, a list of units, or a `Tuple`.
#' @param unit_system A `UnitSystem` or its name.
#' @param across_dimensions Combine units of different dimensions.
#' @param quantity A string, quantity, dimension or expression.
#' @return An expression, or a character vector for `find_unit()`.
#' @examples
#' u <- physics_units()
#' convert_to(u$mile, u$kilometer)
#' convert_to(u$speed_of_light, list(u$meter, u$second))
#' quantity_simplify(u$kilo * u$foot * u$inch)
#' find_unit("volt")
#' @export
convert_to <- function(expr, target_units, unit_system = "SI") {
  .units_ensure()
  unit_system <- get_unit_system(unit_system)
  targets <- .u_targets(target_units)
  handle_adds <- function(e) .add(lapply(.sympy_args(e), function(i) convert_to(i, targets, unit_system)))
  if (is_sym(expr) && expr$type == "Add") return(handle_adds(expr))
  if (is_sym(expr) && expr$type == "Pow" && expr$args[[1L]]$type == "Add")
    return(.upy_pow(handle_adds(expr$args[[1L]]), expr$args[[2L]]))
  expr <- as_sym(expr)
  if (.u_is_function(expr)) expr <- together(expr)
  if (!.is_quantity(expr) && .has_type(expr, "Quantity"))
    expr <- .map_bottom_up(expr, function(z) if (.is_quantity(z)) convert_to(z, targets, unit_system) else z)
  total_sf <- function(e) {
    if (e$type == "Mul") return(Reduce(.upy_mul, lapply(.sympy_args(e), total_sf)))
    if (e$type == "Pow") return(.upy_pow(total_sf(e$args[[1L]]), e$args[[2L]]))
    if (.is_quantity(e)) return(.qm_get_qsf(unit_system, e))
    e
  }
  depmat <- .get_conversion_matrix_for_expr(expr, targets, unit_system)
  if (is.null(depmat)) return(expr)
  expr_scale_factor <- total_sf(expr)
  .upy_mul(expr_scale_factor, .mul(Map(function(u, p) .upy_pow(.upy_mul(.upy_div(S_One, total_sf(u)), u), p),
                                      targets, depmat)))
}

#' @rdname convert_to
#' @export
quantity_simplify <- function(expr, across_dimensions = FALSE, unit_system = NULL) {
  .units_ensure()
  expr <- as_sym(expr)
  if (!length(expr$args) || !.has_type(expr, c("Prefix", "Quantity"))) return(expr)
  ps <- atoms(expr, "Prefix")
  if (length(ps)) expr <- xreplace(expr, lapply(ps, function(p) list(p, p$extra$scale_factor)))
  qs <- atoms(expr, "Quantity")
  if (length(qs)) {
    dkeys <- vapply(qs, function(q) quantity_dimension(q)$key, "")
    for (k in unique(dkeys)) {
      grp <- qs[dkeys == k]
      if (length(grp) == 1L) next
      v <- .bool_print_order(grp)
      ref <- .div(v[[1L]], scale_factor(v[[1L]]))
      expr <- xreplace(expr, lapply(v[-1L], function(vi) list(vi, .mul(list(ref, scale_factor(vi))))))
    }
  }
  if (across_dimensions) {
    if (is.null(unit_system)) .u_stop("ValueError", "unit_system must be specified if across_dimensions is True")
    unit_system <- get_unit_system(unit_system)
    ds <- unit_system$dimension_system
    dim_expr <- get_dimensional_expr(unit_system, expr)
    dim_deps <- get_dimensional_dependencies(ds, dim_expr, mark_dimensionless = TRUE)
    target <- NULL
    for (i in seq_along(ds$deps$k)) {
      if (.ud_eq(ds$deps$v[[i]], dim_deps)) {
        target <- ds$deps$k[[i]]
        break
      }
    }
    if (is.null(target)) return(expr)
    target_unit <- .ud_get(unit_system$derived_units, target)
    if (!is.null(target_unit)) expr <- convert_to(expr, target_unit, unit_system)
  }
  expr
}

# default_sort_key of a Python tuple of (Dimension, power) pairs
.u_tuple_key <- function(items) {
  list(list(10, 0, "tuple"), list(length(items), lapply(items, function(it)
    if (is.list(it) && !is_sym(it)) .u_tuple_key(it) else .sort_key(it))), .ONE_SORT_KEY, 1)
}
.u_tuple_nodes <- function(items) 1 + sum(vapply(items, function(it)
  if (is.list(it) && !is_sym(it)) .u_tuple_nodes(it) else .nodes(it), 0))

#' @rdname convert_to
#' @export
check_dimensions <- function(expr, unit_system = "SI") {
  .units_ensure()
  expr <- as_sym(expr)
  unit_system <- get_unit_system(unit_system)
  add_dict <- function(d1, d2) {
    d3 <- d1
    for (i in seq_along(d2$k)) {
      k <- d2$k[[i]]
      v <- d2$v[[i]]
      v1 <- .ud_get(d1, k)
      d3 <- .ud_set(d3, k, if (is.null(v1)) v else .upy_add(v, v1))
    }
    keep <- vapply(d3$v, function(v) !.u_val_eq(v, S_Zero), TRUE)
    .ud(d3$k[keep], d3$v[keep])
  }
  ds <- .us_dimsys(unit_system, "get_dimensional_dependencies")
  for (a in atoms(expr, "Add")) {
    deset <- list()
    dkeys <- character()
    add_tup <- function(tup) {
      tk <- paste(vapply(tup, function(p) paste0(p[[1L]]$key, "=", p[[2L]]$key), ""), collapse = ";")
      if (!(tk %in% dkeys)) {
        dkeys <<- c(dkeys, tk)
        deset[[length(deset) + 1L]] <<- tup
      }
    }
    for (ai in .sympy_args(a)) {
      if (.u_is_number(ai)) {
        add_tup(list())
        next
      }
      skip <- FALSE
      dimdict <- .ud()
      for (i in (if (ai$type == "Mul") .sympy_args(ai) else list(ai))) {
        if (.has_type(i, "Quantity")) i <- Dimension(get_dimensional_expr(unit_system, i))
        if (.has_type(i, "Dimension")) {
          dimdict <- add_dict(dimdict, get_dimensional_dependencies(ds, i))
        } else if (length(free_symbols(i))) {
          skip <- TRUE
          break
        }
      }
      if (skip) next
      tup <- Map(function(k, v) list(k, v), dimdict$k, dimdict$v)
      if (length(tup) > 1L) tup <- .key_sort(lapply(tup, .u_tuple_key), tup)
      add_tup(tup)
      if (length(deset) > 1L) {
        srt <- .key_sort(lapply(deset, function(t) list(.u_tuple_nodes(t), .u_tuple_key(t))), deset)
        .u_stop("ValueError", paste0("addends have incompatible dimensions: ", .upy_repr(srt)))
      }
    }
  }
  reps <- list()
  for (m in atoms(expr, "Mul")) {
    if (any(vapply(m$args, .is_dim, TRUE)))
      reps[[length(reps) + 1L]] <- list(m, .mul(Filter(function(i) !.u_is_number(i), m$args)))
  }
  if (!length(reps)) return(expr)
  xreplace(expr, reps)
}

# dir(sympy.physics.units) without the dunder names
.U_MODULE_DIR <- strsplit(paste(c(
  "A Bq C D Da Dimension DimensionSystem F G Gy H Hz J K L N Pa Quantity R S T Unit UnitSystem V W Wb Z0",
  "acceleration acceleration_due_to_gravity action amount amount_of_substance ampere amperes amu amus angstrom",
  "angstroms angular_mil angular_mils anomalistic_year anomalistic_years area astronomical_unit",
  "astronomical_units atm atmosphere atmospheres atomic_mass_constant atomic_mass_unit atto au avogadro",
  "avogadro_constant avogadro_number bar bars becquerel bit bits boltzmann boltzmann_constant byte c cL candela",
  "candelas capacitance cd centi centiliter centiliters centimeter centimeters charge cl cm common_year",
  "common_years conductance convert_to coulomb coulomb_constant coulombs current dHg0 dL dalton day days deca",
  "deci deciliter deciliters decimeter decimeters definitions deg degree degrees dimensions dioptre dl dm",
  "draconic_year draconic_years e0 eV electric_constant electric_force_constant electron_rest_mass electronvolt",
  "electronvolts elementary_charge energy exa exbi exbibyte exbibytes farad faraday_constant farads feet femto",
  "find_unit foot force frequency ft full_moon_cycle full_moon_cycles g gaussian_year gaussian_years gee gees",
  "gibi gibibyte gibibytes giga gram grams gravitational_constant gray h ha hbar hectare hecto henry henrys hertz",
  "hour hours hz impedance inch inches inductance josephson_constant joule joules julian_year julian_years kPa",
  "kat katal kelvin kelvins kg kibi kibibyte kibibytes kilo kilogram kilograms kilometer kilometers km l length",
  "lightyear lightyears liter liters luminosity luminous_intensity lux lx ly m mL magnetic_constant",
  "magnetic_density magnetic_flux magnetic_flux_density mass me mebi mebibyte mebibytes mega meter meters",
  "metric_ton mg mho mhos mi micro microgram micrograms micrometer micrometers micron microns microsecond",
  "microseconds mil mile miles milli milli_mass_unit milligram milligrams milliliter milliliters millimeter",
  "millimeters millisecond milliseconds minute minutes mks mksa ml mm mmHg mmu mmus mol molar_gas_constant mole",
  "moles momentum ms nano nanometer nanometers nanosecond nanoseconds nautical_mile nautical_miles newton newtons",
  "nm nmi ns ohm ohms optical_power pa pascal pascals pebi pebibyte pebibytes percent percents permille peta pico",
  "picometer picometers picosecond picoseconds planck planck_acceleration planck_angular_frequency planck_area",
  "planck_charge planck_current planck_density planck_energy planck_energy_density planck_force planck_impedance",
  "planck_intensity planck_length planck_mass planck_momentum planck_power planck_pressure planck_temperature",
  "planck_time planck_voltage planck_volume pm pound pounds power prefixes pressure ps psi quantities quart",
  "quarts rad radian radians s second seconds si sidereal_year sidereal_years siemens speed speed_of_light sr",
  "stefan stefan_boltzmann_constant steradian steradians systems t tebi tebibyte tebibytes temperature tera tesla",
  "teslas time tonne torr tropical_year tropical_years u0 ug um unitsystem us util v vacuum_impedance",
  "vacuum_permeability vacuum_permittivity velocity volt voltage volts volume von_klitzing_constant watt watts wb",
  "weber webers yard yards yd year years yocto yotta zepto zetta"
), collapse = " "), " ", fixed = TRUE)[[1L]]

#' @rdname convert_to
#' @export
find_unit <- function(quantity, unit_system = "SI") {
  .units_ensure()
  unit_system <- get_unit_system(unit_system)
  ns <- .units_state$ns
  getq <- function(nm) if (nm %in% .U_MODULE_DIR) ns[[nm]] else NULL
  rv <- character()
  if (is.character(quantity)) {
    rv <- Filter(function(i) grepl(quantity, i, fixed = TRUE) && .is_quantity(getq(i)), .U_MODULE_DIR)
    if (!(quantity %in% .U_MODULE_DIR))
      .u_stop("AttributeError", paste0("module 'sympy.physics.units' has no attribute '", quantity, "'"))
    dim <- getq(quantity)
    if (.is_dim(dim)) rv <- c(rv, find_unit(dim))
  } else {
    quantity <- as_sym(quantity)
    for (i in sort(.U_MODULE_DIR, method = "radix")) {
      other <- getq(i)
      if (!.is_quantity(other)) next
      od <- quantity_dimension(other)
      if (.is_quantity(quantity)) {
        if (identical(quantity_dimension(quantity)$key, od$key)) rv <- c(rv, i)
      } else if (.is_dim(quantity)) {
        if (identical(od$key, quantity$key)) rv <- c(rv, i)
      } else if (identical(od$key, Dimension(get_dimensional_expr(unit_system, quantity))$key)) {
        rv <- c(rv, i)
      }
    }
  }
  rv <- unique(rv)
  rv[order(nchar(rv), rv, method = "radix")]
}
