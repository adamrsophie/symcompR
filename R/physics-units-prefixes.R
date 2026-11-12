# Unit prefixes, after sympy/physics/units/prefixes.py.

#' Unit prefixes
#'
#' `Prefix(name, abbrev, exponent, base)` is SymPy's prefix with scale factor
#' `base^exponent`, such as `kilo` (`10^3`) or `kibi` (`2^10`). A prefix times
#' a quantity gives the scaled quantity, and products and quotients of
#' prefixes simplify to a prefix or a number, as in `sympy.physics.units`.
#' The predefined prefixes are in [physics_units()], with the dictionaries
#' `PREFIXES` and `BIN_PREFIXES`. `prefix_unit(unit, prefixes)` returns the
#' prefixed units of `unit` for a named list of prefixes. `scale_factor()`
#' gives the scale factor of a prefix or of a quantity in the SI system.
#'
#' @param name The name.
#' @param abbrev The abbreviation.
#' @param exponent The exponent of the base.
#' @param base The base, 10 by default.
#' @param latex_repr Optional LaTeX form.
#' @param unit A quantity.
#' @param prefixes A named list of prefixes.
#' @param x A prefix or quantity.
#' @return A `sym` Prefix, a list of quantities, or an expression.
#' @examples
#' u <- physics_units()
#' u$kilo * u$meter
#' u$kilo * u$kilo
#' prefix_unit(u$meter, u$PREFIXES[c("m", "c")])
#' @export
Prefix <- function(name, abbrev, exponent, base = 10, latex_repr = NULL) {
  name <- if (is.character(name)) sympify(name) else as_sym(name)
  abbrev <- if (is.character(abbrev)) Symbol(abbrev) else as_sym(abbrev)
  exponent <- as_sym(exponent)
  base <- as_sym(base)
  new_sym("Prefix", list(name, abbrev, exponent, base),
          extra = list(scale_factor = .pow(base, exponent), latex = latex_repr))
}

register_str_printer("Prefix", function(e) .str(e$args[[1L]]))
register_rebuilder("Prefix", function(e, args) Prefix(args[[1L]], args[[2L]], args[[3L]], args[[4L]]))

#' @rdname Prefix
#' @export
scale_factor <- function(x) {
  x <- as_sym(x)
  if (.is_prefix(x)) return(x$extra$scale_factor)
  if (.is_quantity(x)) return(get_quantity_scale_factor(get_default_unit_system(), x))
  .u_stop("AttributeError", paste0("'", .sympy_class(x), "' object has no attribute 'scale_factor'"))
}

# Prefix._prefix_with_scale_factor
.prefix_with_scale_factor <- function(fact) {
  .units_ensure()
  for (prefixes in list(.units_state$ns$PREFIXES, .units_state$ns$BIN_PREFIXES))
    for (p in prefixes) if (.u_val_eq(p$extra$scale_factor, fact)) return(p)
  NULL
}

.u_has_scale_factor <- function(x) .is_prefix(x) || .is_quantity(x)

# Prefix.__mul__
.prefix_mul <- function(self, other) {
  if (!.u_has_scale_factor(other)) return(.mul(list(self, other)))
  fact <- .upy_mul(self$extra$scale_factor, scale_factor(other))
  if (.is_prefix(other)) {
    if (.u_val_eq(fact, S_One)) return(S_One)
    p <- .prefix_with_scale_factor(fact)
    if (!is.null(p)) return(p)
    return(fact)
  }
  .upy_mul(self$extra$scale_factor, other)
}

# Prefix.__truediv__
.prefix_truediv <- function(self, other) {
  if (!.u_has_scale_factor(other)) return(.mul(list(self, .pow(other, S_NegativeOne))))
  fact <- .upy_div(self$extra$scale_factor, scale_factor(other))
  if (.u_val_eq(fact, S_One)) return(S_One)
  if (.is_prefix(other)) {
    p <- .prefix_with_scale_factor(fact)
    if (!is.null(p)) return(p)
    return(fact)
  }
  .upy_div(self$extra$scale_factor, other)
}

# Prefix.__rtruediv__
.prefix_rtruediv <- function(self, other) {
  if (.u_val_eq(other, S_One)) {
    p <- .prefix_with_scale_factor(.div(S_One, self$extra$scale_factor))
    if (!is.null(p)) return(p)
  }
  .upy_div(other, self$extra$scale_factor)
}

#' @rdname Prefix
#' @export
prefix_unit <- function(unit, prefixes) {
  .units_ensure()
  out <- list()
  for (p in prefixes) {
    q <- Quantity(paste0(.str(p$args[[1L]]), unit$name),
                  abbrev = paste0(.str(p$args[[2L]]), .str(unit$extra$abbrev)), is_prefixed = TRUE)
    assign(q$key, unit, envir = .units_state$equiv_global)
    assign(q$key, list(p$extra$scale_factor, unit), envir = .units_state$sf_global)
    out[[length(out) + 1L]] <- q
  }
  out
}
