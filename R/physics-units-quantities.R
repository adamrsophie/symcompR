# Physical quantities, after sympy/physics/units/quantities.py.
#
# A Quantity is an atom of type "Quantity" (SymPy's AtomicExpr) that prints
# as its name. Its identity is its class, name and abbreviation, so the
# registries of scale factors and dimensions are keyed by the sym key.
# PhysicalConstant shares the type and carries extra$pc, so that
# .has_type(e, "Quantity") finds both, as isinstance() does in SymPy.

.new_quantity <- function(name, abbrev, latex_repr, is_prefixed, pc) {
  if (is_sym(name)) name <- name$name
  abbrev <- if (is.null(abbrev)) Symbol(name) else if (is.character(abbrev)) Symbol(abbrev) else as_sym(abbrev)
  new_sym("Quantity", name = name,
          extra = list(abbrev = abbrev, pc = pc, latex = latex_repr, is_prefixed = isTRUE(is_prefixed)),
          key = paste0(if (pc) "PQ" else "Q", "(", Symbol(name)$key, ",", abbrev$key, ")"))
}

#' Physical quantities
#'
#' `Quantity(name, abbrev)` is SymPy's physical quantity: a unit of measure,
#' a constant or a generic quantity. It is a positive real atom that prints
#' as its name and takes part in arithmetic like a symbol (`Unit` is an
#' alias, as in `sympy.physics.units`). `PhysicalConstant()` marks a physical
#' constant such as `speed_of_light`. A quantity has no scale factor or
#' dimension of its own: set them globally with
#' `set_global_relative_scale_factor()` and `set_global_dimension()`, or per
#' system with [set_quantity_scale_factor()] and [set_quantity_dimension()].
#' `quantity_dimension()` and [scale_factor()] read them in the SI system.
#' The predefined units and constants are in [physics_units()].
#'
#' @param name The name, a string or symbol.
#' @param abbrev The abbreviation; the name by default.
#' @param latex_repr,pretty_unicode_repr,pretty_ascii_repr,mathml_presentation_repr
#'   Optional printed forms.
#' @param is_prefixed Whether the quantity carries a prefix.
#' @param ... Ignored assumptions.
#' @param q A quantity.
#' @param scale_factor A number or expression, possibly with prefixes.
#' @param reference_quantity The quantity the scale factor refers to.
#' @param dimension A `Dimension`.
#' @return A `sym` Quantity, an expression, `TRUE`/`FALSE`, or `NULL`
#'   invisibly for the setters.
#' @examples
#' u <- physics_units()
#' q <- Quantity("q1")
#' set_global_relative_scale_factor(q, 5000, u$meter)
#' convert_to(q, u$meter)
#' scale_factor(u$kilometer)
#' quantity_dimension(u$newton)
#' @export
Quantity <- function(name, abbrev = NULL, latex_repr = NULL, pretty_unicode_repr = NULL,
                     pretty_ascii_repr = NULL, mathml_presentation_repr = NULL, is_prefixed = FALSE, ...) {
  .new_quantity(name, abbrev, latex_repr, is_prefixed, pc = FALSE)
}

#' @rdname Quantity
#' @export
PhysicalConstant <- function(name, abbrev = NULL, latex_repr = NULL, pretty_unicode_repr = NULL,
                             pretty_ascii_repr = NULL, mathml_presentation_repr = NULL, is_prefixed = FALSE, ...) {
  .new_quantity(name, abbrev, latex_repr, is_prefixed, pc = TRUE)
}

#' @rdname Quantity
#' @export
Unit <- Quantity

register_str_printer("Quantity", function(e) e$name)
register_latex_printer("Quantity", function(e) {
  if (!is.null(e$extra$latex)) e$extra$latex else paste0("\\text{", .str(e$extra$abbrev), "}")
})
.facts_hooks$Quantity <- function(e) c(commutative = TRUE, real = TRUE, nonzero = TRUE, positive = TRUE)

#' @rdname Quantity
#' @export
set_global_dimension <- function(q, dimension) {
  .units_ensure()
  assign(q$key, dimension, envir = .units_state$dim_global)
  invisible(NULL)
}

#' @rdname Quantity
#' @export
set_global_relative_scale_factor <- function(q, scale_factor, reference_quantity) {
  .units_ensure()
  sf <- as_sym(scale_factor)
  if (.is_prefix(sf)) assign(q$key, TRUE, envir = .units_state$prefixed)
  sf <- .map_bottom_up(sf, function(z) if (.is_prefix(z)) z$extra$scale_factor else z)
  reference_quantity <- as_sym(reference_quantity)
  assign(q$key, list(sf, reference_quantity), envir = .units_state$sf_global)
  assign(q$key, reference_quantity, envir = .units_state$equiv_global)
  invisible(NULL)
}

#' @rdname Quantity
#' @export
quantity_name <- function(q) Symbol(q$name)

#' @rdname Quantity
#' @export
quantity_abbrev <- function(q) q$extra$abbrev

#' @rdname Quantity
#' @export
quantity_dimension <- function(q) get_quantity_dimension(get_default_unit_system(), q)

#' @rdname Quantity
#' @export
is_prefixed <- function(q) isTRUE(q$extra$is_prefixed) || isTRUE(.units_state$prefixed[[q$key]])

#' @rdname Quantity
#' @export
is_physical_constant <- function(q) .is_quantity(q) && isTRUE(q$extra$pc)

# evalf keeps quantities, dimensions and prefixes as they are, like
# symbols: N() evaluates with placeholders in their place
.units_evalf <- function(e, digits, subs) {
  nodes <- atoms(e, c("Quantity", "Dimension", "Prefix"))
  ph <- lapply(nodes, function(z) if (z$type == "Prefix") Dummy("p") else Dummy("q", positive = TRUE))
  r <- N(xreplace(e, Map(list, nodes, ph)), digits, subs)
  xreplace(r, Map(list, ph, nodes))
}
