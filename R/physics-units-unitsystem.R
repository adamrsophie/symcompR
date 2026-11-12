# Unit systems, after sympy/physics/units/unitsystem.py.

#' Unit systems
#'
#' `UnitSystem(base_units, units, name, descr, dimension_system,
#' derived_units)` is SymPy's coherent set of units on top of a
#' [DimensionSystem()]. Creating one registers it under its name for
#' `get_unit_system()`; the predefined systems are `"SI"`, `"MKS"`, `"MKSA"`,
#' `"Natural system"` and `"cgs_gauss"`. `get_dimensional_expr()` gives the
#' dimension of an expression as an expression in the dimension names, and
#' `collect_factor_and_dimension()` (SymPy's `_collect_factor_and_dimension`)
#' returns `list(scale_factor, dimension)`, checking that added terms and
#' function arguments agree in dimension. `get_units_non_prefixed()` lists
#' the units without a prefix that are not physical constants.
#'
#' @param base_units A list of base quantities.
#' @param units,base A list of further quantities.
#' @param name The name of the system.
#' @param descr,description A description.
#' @param dimension_system A `DimensionSystem`.
#' @param derived_units Preferred units of dimensions, as a list of
#'   `list(dimension, unit)` pairs or a `units_dict`.
#' @param unit_system,us,x A `UnitSystem` or its name.
#' @param expr An expression.
#' @param ... Unused.
#' @return A `UnitSystem`, an expression, a list or a `units_dict`.
#' @examples
#' u <- physics_units()
#' SI <- get_unit_system("SI")
#' get_dimensional_expr(SI, u$newton / u$meter)
#' collect_factor_and_dimension(SI, u$meter + u$kilometer)
#' @export
UnitSystem <- function(base_units, units = list(), name = "", descr = "", dimension_system = NULL,
                       derived_units = list()) {
  .units_ensure()
  base_units <- lapply(as.list(base_units), as_sym)
  all <- c(base_units, lapply(as.list(units), as_sym))
  obj <- new.env()
  obj$name <- name
  obj$descr <- descr
  obj$base_units <- base_units
  obj$units <- all[!duplicated(vapply(all, .k, ""))]
  obj$dimension_system <- dimension_system
  obj$derived_units <- .as_od(derived_units)
  class(obj) <- "UnitSystem"
  obj <- .u_mapper_init(obj)
  nms <- vapply(.units_state$systems, `[[`, "", 1L)
  i <- match(name, nms)
  if (is.na(i)) i <- length(nms) + 1L
  .units_state$systems[[i]] <- list(name, obj)
  obj
}

#' @export
format.UnitSystem <- function(x, ...) {
  if (!identical(x$name, "")) return(x$name)
  paste0("UnitSystem((", paste(vapply(x$base_units, .str, ""), collapse = ", "), "))")
}

#' @export
print.UnitSystem <- function(x, ...) {
  cat(format(x), "\n", sep = "")
  invisible(x)
}

#' @rdname UnitSystem
#' @export
extend.UnitSystem <- function(x, base, units = list(), name = "", description = "", dimension_system = NULL,
                              derived_units = list(), ...) {
  du <- x$derived_units
  nd <- .as_od(derived_units)
  for (i in seq_along(nd$k)) du <- .ud_set(du, nd$k[[i]], nd$v[[i]])
  UnitSystem(c(x$base_units, as.list(base)), c(x$units, as.list(units)), name, description, dimension_system, du)
}

#' @rdname UnitSystem
#' @export
get_dimension_system <- function(us) us$dimension_system

#' @rdname UnitSystem
#' @export
get_unit_system <- function(unit_system) {
  .units_ensure()
  if (inherits(unit_system, "UnitSystem")) return(unit_system)
  nms <- vapply(.units_state$systems, `[[`, "", 1L)
  i <- if (is.character(unit_system) && length(unit_system) == 1L) match(unit_system, nms) else NA_integer_
  us <- if (is.na(i)) NULL else .units_state$systems[[i]][[2L]]
  if (is.null(us)) {
    .u_stop("ValueError", paste0("Unit system is not supported. Currentlysupported unit systems are ",
                                 paste(nms[order(nms, method = "radix")], collapse = ", ")))
  }
  us
}

#' @rdname UnitSystem
#' @export
get_default_unit_system <- function() get_unit_system("SI")

#' @rdname UnitSystem
#' @export
derived_units <- function(us) us$derived_units

# the SymPy Function classes; Min and Max are applications only
.u_is_function <- function(e) e$type == "Function" && !(e$name %in% .BASIC_KEY_FUNS)

#' @rdname UnitSystem
#' @export
get_dimensional_expr <- function(us, expr) {
  .units_ensure()
  expr <- as_sym(expr)
  t <- expr$type
  if (t == "Mul") return(.mul(lapply(.sympy_args(expr), function(i) get_dimensional_expr(us, i))))
  if (t == "Pow") return(.upy_pow(get_dimensional_expr(us, expr$args[[1L]]), expr$args[[2L]]))
  if (t == "Add") return(get_dimensional_expr(us, .sympy_args(expr)[[1L]]))
  if (t == "Derivative") {
    dim <- get_dimensional_expr(us, expr$args[[1L]])
    for (vc in expr$args[-1L])
      dim <- .upy_div(dim, .upy_pow(get_dimensional_expr(us, vc$args[[1L]]), vc$args[[2L]]))
    return(dim)
  }
  if (.u_is_function(expr)) {
    args <- lapply(expr$args, function(a) get_dimensional_expr(us, a))
    if (all(vapply(args, function(i) .u_val_eq(i, S_One), TRUE))) return(S_One)
    return(fcall(expr$name, args))
  }
  if (t == "Quantity") return(.qm_get_qdim(us, expr)$args[[1L]])
  S_One
}

#' @rdname UnitSystem
#' @export
collect_factor_and_dimension <- function(us, expr) {
  .units_ensure()
  .u_cfd(us, as_sym(expr))
}

.u_cfd <- function(us, expr) {
  t <- expr$type
  ds <- function() .us_dimsys(us, "is_dimensionless")
  if (t == "Quantity") return(list(scale_factor(expr), quantity_dimension(expr)))
  if (t == "Mul") {
    factor <- S_One
    dimension <- Dimension(S_One)
    for (arg in .sympy_args(expr)) {
      fd <- .u_cfd(us, arg)
      factor <- .upy_mul(factor, fd[[1L]])
      dimension <- .upy_mul(dimension, fd[[2L]])
    }
    return(list(factor, dimension))
  }
  if (t == "Pow") {
    fd <- .u_cfd(us, expr$args[[1L]])
    ed <- .u_cfd(us, expr$args[[2L]])
    exp_dim <- if (is_dimensionless(ds(), ed[[2L]])) S_One else ed[[2L]]
    return(list(.upy_pow(fd[[1L]], ed[[1L]]), .upy_pow(fd[[2L]], .upy_mul(ed[[1L]], exp_dim))))
  }
  if (t == "Add") {
    args <- .sympy_args(expr)
    fd <- .u_cfd(us, args[[1L]])
    factor <- fd[[1L]]
    dim <- fd[[2L]]
    for (addend in args[-1L]) {
      ad <- .u_cfd(us, addend)
      if (!equivalent_dims(ds(), dim, ad[[2L]]))
        .u_stop("ValueError", paste0("Dimension of \"", .str(addend), "\" is ", .str(ad[[2L]]), ", but it should be ", .str(dim)))
      factor <- .upy_add(factor, ad[[1L]])
    }
    return(list(factor, dim))
  }
  if (t == "Derivative") {
    fd <- .u_cfd(us, expr$args[[1L]])
    factor <- fd[[1L]]
    dim <- fd[[2L]]
    for (vc in expr$args[-1L]) {
      id <- .u_cfd(us, vc$args[[1L]])
      factor <- .upy_div(factor, .upy_pow(id[[1L]], vc$args[[2L]]))
      dim <- .upy_div(dim, .upy_pow(id[[2L]], vc$args[[2L]]))
    }
    return(list(factor, dim))
  }
  if (t == "Function") {
    fds <- lapply(expr$args, function(a) .u_cfd(us, a))
    dims <- lapply(fds, function(d) if (is_dimensionless(ds(), d[[2L]])) Dimension(S_One) else d[[2L]])
    dim <- if (length(dims)) dims[[1L]] else Dimension(S_One)
    for (d in dims[-1L]) {
      if (!equivalent_dims(ds(), dim, d))
        .u_stop("ValueError", paste0("Dimension of \"", .str(expr), "\" is ", .str(d), ", but it should be ", .str(dim)))
    }
    return(list(fcall(expr$name, lapply(fds, `[[`, 1L)), dim))
  }
  if (t == "Dimension") return(list(S_One, expr))
  list(expr, Dimension(S_One))
}

#' @rdname UnitSystem
#' @export
get_units_non_prefixed <- function(us) {
  Filter(function(u) !is_prefixed(u) && !is_physical_constant(u), us$units)
}
