# Physical dimensions and dimension systems, after
# sympy/physics/units/dimensions.py.
#
# A Dimension is a sym of type "Dimension" whose single argument is its name;
# the optional symbol sits in `extra` and, as in SymPy, takes no part in
# equality. Python's operators on dimensions and prefixes (`__mul__`,
# `__add__`, ...) differ from the Mul and Add constructors, so the R
# operators reach them through .ops_hooks and the ported algorithms call the
# .py_* helpers wherever SymPy uses an operator.

# ---------------------------------------------------------------------------
# Shared state and small helpers

.units_state <- new.env()
.units_state$built <- FALSE
.units_state$sf_global <- new.env()     # _quantity_scale_factors_global
.units_state$equiv_global <- new.env()  # _quantity_dimensional_equivalence_map_global
.units_state$dim_global <- new.env()    # _quantity_dimension_global
.units_state$prefixed <- new.env()      # quantities marked prefixed after creation
.units_state$systems <- list()          # UnitSystem._unit_systems

.u_stop <- function(cls, msg) stop(structure(class = c(cls, "error", "condition"), list(message = msg, call = NULL)))

.is_dim <- function(x) is_sym(x) && x$type == "Dimension"
.is_quantity <- function(x) is_sym(x) && x$type == "Quantity"
.is_prefix <- function(x) is_sym(x) && x$type == "Prefix"

# Python dicts keep insertion order; a key compares by its sym key (or by
# its string, for the string keys DimensionSystem accepts), and a repeated
# key keeps the first key object with the last value.
.ud <- function(keys = list(), values = list()) {
  d <- structure(list(k = list(), v = list()), class = "units_dict")
  for (i in seq_along(keys)) d <- .ud_set(d, keys[[i]], values[[i]])
  d
}
.ud_id <- function(key) if (is.character(key)) paste0("\001", key) else key$key
.ud_idx <- function(d, key) {
  if (!length(d$k)) return(NA_integer_)
  match(.ud_id(key), vapply(d$k, .ud_id, ""))
}
.ud_has <- function(d, key) !is.na(.ud_idx(d, key))
.ud_get <- function(d, key, default = NULL) {
  i <- .ud_idx(d, key)
  if (is.na(i)) default else d$v[[i]]
}
.ud_set <- function(d, key, value) {
  i <- .ud_idx(d, key)
  if (is.na(i)) {
    d$k[[length(d$k) + 1L]] <- key
    d$v[[length(d$v) + 1L]] <- value
  } else {
    d$v[[i]] <- value
  }
  d
}
.ud_len <- function(d) length(d$k)
# Python ==: numbers compare by value, everything else structurally
.u_val_eq <- function(a, b) {
  if (inherits(a, "units_dict") || inherits(b, "units_dict"))
    return(inherits(a, "units_dict") && inherits(b, "units_dict") && .ud_eq(a, b))
  if (identical(a$key, b$key)) return(TRUE)
  if (is_number_type(a) && is_number_type(b)) return(num_cmp(a, b) == 0)
  FALSE
}
.ud_eq <- function(a, b) {
  if (.ud_len(a) != .ud_len(b)) return(FALSE)
  for (i in seq_along(a$k)) {
    j <- .ud_idx(b, a$k[[i]])
    if (is.na(j) || !.u_val_eq(a$v[[i]], b$v[[j]])) return(FALSE)
  }
  TRUE
}

# Python's repr of the pieces of a result
.upy_repr <- function(x) {
  if (is.null(x)) return("None")
  if (is.logical(x)) return(if (isTRUE(x)) "True" else "False")
  if (is.character(x)) return(paste0("'", x, "'"))
  if (inherits(x, "units_dict")) return(format(x))
  # repr() is sstr(), which prints a prefix by its name
  if (is_sym(x)) return(if (x$type == "Prefix") .str(x) else .str_top(x))
  if (is.numeric(x)) return(format(x))
  if (is.list(x)) {
    inner <- vapply(x, .upy_repr, "")
    if (length(x) == 1L) return(paste0("(", inner, ",)"))
    return(paste0("(", paste(inner, collapse = ", "), ")"))
  }
  format(x)
}

# A SymPy Dict (marked "sorted") prints its keys in default_sort_key order,
# a Python dict in insertion order.
#' @export
format.units_dict <- function(x, ...) {
  if (!length(x$k)) return("{}")
  o <- seq_along(x$k)
  if (isTRUE(attr(x, "sorted"))) o <- .key_sort(lapply(x$k, .sort_key), o)
  paste0("{", paste(vapply(o, function(i) paste0(.upy_repr(x$k[[i]]), ": ", .upy_repr(x$v[[i]])), ""),
                    collapse = ", "), "}")
}

#' @export
print.units_dict <- function(x, ...) {
  cat(format(x), "\n", sep = "")
  invisible(x)
}

# Expr.is_number: quantities and dimensions are never numbers
.u_is_number <- function(e) .expr_is_number(e)

# The classes of SymPy that are not Expr, so cannot name a dimension
.U_NON_EXPR <- c("Tuple", "Dict", "BooleanAtom", "And", "Or", "Not", "Xor", "Implies", "Equivalent",
                 "FiniteSet", "Interval", "Union", "Intersection", "Complement", "EmptySet", "SpecialSet",
                 "Equality", "Unequality", "StrictLessThan", "LessThan", "StrictGreaterThan", "GreaterThan")

# ---------------------------------------------------------------------------
# Dimension

#' Physical dimensions
#'
#' `Dimension(name, symbol)` is SymPy's `Dimension`: a name such as `"length"`
#' (or any expression) and an optional symbol such as `"L"`. Dimensions
#' multiply, divide and take powers to form new dimensions, as in
#' `sympy.physics.units`; adding equal dimensions gives the dimension again
#' and adding a quantity is an error. Equality ignores the symbol.
#' `dimension_name()` and `dimension_symbol()` return the parts, and
#' `has_integer_powers()` tells whether every exponent of the dependencies in
#' a dimension system is an integer.
#'
#' @param name A string or expression.
#' @param symbol `NULL`, a string or a symbol.
#' @param d A `Dimension`.
#' @param dim_sys A `DimensionSystem`.
#' @return A `sym` Dimension, an expression, or `TRUE`/`FALSE`.
#' @examples
#' u <- physics_units()
#' L <- Dimension("length", "L")
#' L / Dimension("time")
#' has_integer_powers(sqrt(u$length), u$dimsys_SI)
#' @export
Dimension <- function(name, symbol = NULL) {
  if (is.character(name) && length(name) == 1L) {
    name <- Symbol(name)
  } else if (is.list(name) && !is_sym(name)) {
    .u_stop("TypeError", "Dimension name needs to be a valid math expression")
  } else {
    name <- as_sym(name)
  }
  if (name$type %in% .U_NON_EXPR) .u_stop("TypeError", "Dimension name needs to be a valid math expression")
  if (is.character(symbol)) symbol <- Symbol(symbol)
  else if (!is.null(symbol) && !(is_sym(symbol) && symbol$type == "Symbol")) .u_stop("AssertionError", "")
  new_sym("Dimension", list(name), extra = list(symbol = symbol))
}

#' @rdname Dimension
#' @export
dimension_name <- function(d) d$args[[1L]]

#' @rdname Dimension
#' @export
dimension_symbol <- function(d) d$extra$symbol

#' @rdname Dimension
#' @export
has_integer_powers <- function(d, dim_sys) {
  deps <- get_dimensional_dependencies(dim_sys, d)
  all(vapply(deps$v, function(p) p$type == "Integer", TRUE))
}

register_str_printer("Dimension", function(e) {
  s <- e$extra$symbol
  if (is.null(s)) paste0("Dimension(", .str(e$args[[1L]]), ")")
  else paste0("Dimension(", .str(e$args[[1L]]), ", ", .str(s), ")")
})
# func(*args) builds a fresh Dimension, without the symbol
register_rebuilder("Dimension", function(e, args) Dimension(args[[1L]]))
.facts_hooks$Dimension <- function(e) c(commutative = TRUE, positive = TRUE, real = TRUE)

# Dimension._eval_power, reached from the Pow constructor
.dim_eval_power <- function(b, e) Dimension(.upy_pow(b$args[[1L]], e))

# Dimension.__add__
.dim_add <- function(self, other) {
  if (.has_type(other, "Quantity")) .u_stop("TypeError", "cannot sum dimension and quantity")
  if (.is_dim(other) && identical(self$key, other$key)) return(self)
  .add(list(self, other))
}

# Dimension.__mul__
.dim_mul <- function(self, other) {
  if (.has_type(other, "Quantity")) .u_stop("TypeError", "cannot sum dimension and quantity")
  if (.is_dim(other)) return(Dimension(.upy_mul(self$args[[1L]], other$args[[1L]])))
  if (!length(free_symbols(other))) return(self)
  .mul(list(self, other))
}

# Python's binary operators with SymPy's _op_priority dispatch: Dimension
# and Prefix (priority 13) win over every other expression.
.upy_add <- function(a, b) {
  if (.is_dim(a)) return(.dim_add(a, b))
  if (.is_dim(b)) return(.dim_add(b, a))
  .add(list(a, b))
}
.upy_sub <- function(a, b) {
  if (.is_dim(a)) return(.upy_add(a, b))
  if (.is_dim(b)) return(.upy_add(b, a))
  .sub(a, b)
}
.upy_mul <- function(a, b) {
  if (.is_dim(a)) return(.dim_mul(a, b))
  if (.is_prefix(a)) return(.prefix_mul(a, b))
  if (.is_dim(b)) return(.dim_mul(b, a))
  .mul(list(a, b))
}
.upy_div <- function(a, b) {
  if (.is_dim(a)) return(.upy_mul(a, .pow(b, S_NegativeOne)))
  if (.is_prefix(a)) return(.prefix_truediv(a, b))
  if (.is_dim(b)) return(.upy_mul(a, .upy_pow(b, S_NegativeOne)))
  if (.is_prefix(b)) return(.prefix_rtruediv(b, a))
  .div(a, b)
}
.upy_pow <- function(a, b) {
  if (.is_dim(a)) return(.dim_eval_power(a, b))
  .pow(a, b)
}
.upy_neg <- function(a) if (.is_dim(a)) a else .neg(a)

.upy_binop <- function(op, a, b) {
  switch(op,
         "+" = .upy_add(a, b), "-" = .upy_sub(a, b), "*" = .upy_mul(a, b),
         "/" = .upy_div(a, b), "^" = .upy_pow(a, b), NULL)
}

local({
  old <- .ops_hooks$dispatch
  .ops_hooks$dispatch <- function(gen, e1, e2) {
    if ((is_sym(e1) && e1$type %in% c("Dimension", "Prefix")) || (is_sym(e2) && e2$type %in% c("Dimension", "Prefix"))) {
      r <- .upy_binop(gen, as_sym(e1), as_sym(e2))
      if (!is.null(r)) return(r)
    }
    if (is.null(old)) NULL else old(gen, e1, e2)
  }
})

# Dimension._from_dimensional_dependencies
.dim_from_deps <- function(deps) {
  r <- S_One
  for (i in seq_along(deps$k)) r <- .upy_mul(r, .upy_pow(deps$k[[i]], deps$v[[i]]))
  r
}

# ---------------------------------------------------------------------------
# _QuantityMapper: the quantity maps of a dimension or unit system

.u_mapper_init <- function(obj) {
  obj$qdm <- new.env()
  obj$qsf <- new.env()
  obj
}

#' Dimensions and scale factors of quantities in a system
#'
#' The `_QuantityMapper` methods shared by [DimensionSystem()] and
#' [UnitSystem()]. `set_quantity_dimension()` and
#' `set_quantity_scale_factor()` record a quantity's dimension and scale
#' factor in one system (prefixes and quantities inside the scale factor are
#' replaced by their scale factors at once); `get_quantity_dimension()` and
#' `get_quantity_scale_factor()` look them up, falling back to the global
#' definitions made with [set_global_dimension()] and
#' [set_global_relative_scale_factor()].
#'
#' @param system A `DimensionSystem` or `UnitSystem`.
#' @param quantity,unit A quantity (or, for the scale factor, a product or
#'   power of quantities).
#' @param dimension A `Dimension` or `1`.
#' @param scale_factor An expression.
#' @return `NULL` invisibly for the setters; a `Dimension` or expression.
#' @examples
#' u <- physics_units()
#' get_quantity_dimension(u$SI, u$newton)
#' get_quantity_scale_factor(u$SI, u$kilometer * u$gram)
#' @export
set_quantity_dimension <- function(system, quantity, dimension) {
  .units_ensure()
  dimension <- as_sym(dimension)
  if (!.is_dim(dimension)) {
    if (.u_val_eq(dimension, S_One)) dimension <- Dimension(S_One)
    else .u_stop("ValueError", "expected dimension or 1")
  }
  assign(quantity$key, dimension, envir = system$qdm)
  invisible(NULL)
}

#' @rdname set_quantity_dimension
#' @export
set_quantity_scale_factor <- function(system, quantity, scale_factor) {
  .units_ensure()
  sf <- as_sym(scale_factor)
  sf <- .map_bottom_up(sf, function(z) if (.is_prefix(z)) z$extra$scale_factor else z)
  sf <- .map_bottom_up(sf, function(z) if (.is_quantity(z)) .qm_get_qsf(system, z) else z)
  assign(quantity$key, sf, envir = system$qsf)
  invisible(NULL)
}

#' @rdname set_quantity_dimension
#' @export
get_quantity_dimension <- function(system, unit) {
  .units_ensure()
  .qm_get_qdim(system, as_sym(unit))
}

#' @rdname set_quantity_dimension
#' @export
get_quantity_scale_factor <- function(system, unit) {
  .units_ensure()
  .qm_get_qsf(system, as_sym(unit))
}

.us_dimsys <- function(us, field) {
  ds <- us$dimension_system
  if (is.null(ds)) .u_stop("AttributeError", paste0("'NoneType' object has no attribute '", field, "'"))
  ds
}

.qm_get_qdim <- function(obj, unit) {
  k <- unit$key
  if (inherits(obj, "UnitSystem")) {
    v <- .us_dimsys(obj, "_quantity_dimension_map")$qdm[[k]]
    if (!is.null(v)) return(v)
  }
  v <- obj$qdm[[k]]
  if (!is.null(v)) return(v)
  v <- .units_state$dim_global[[k]]
  if (!is.null(v)) return(v)
  dep <- .units_state$equiv_global[[k]]
  if (!is.null(dep)) {
    if (.is_quantity(dep)) return(.qm_get_qdim(obj, dep))
    if (!inherits(obj, "UnitSystem")) .u_stop("AttributeError", "'DimensionSystem' object has no attribute 'get_dimensional_expr'")
    return(Dimension(get_dimensional_expr(obj, dep)))
  }
  if (.is_quantity(unit)) return(Dimension(Symbol(unit$name)))
  Dimension(S_One)
}

.qm_get_qsf <- function(obj, unit) {
  k <- unit$key
  if (inherits(obj, "UnitSystem")) {
    v <- .us_dimsys(obj, "_quantity_scale_factors")$qsf[[k]]
    if (!is.null(v)) return(v)
  }
  if (unit$type == "Mul") return(.mul(lapply(.sympy_args(unit), function(a) .qm_get_qsf(obj, a))))
  if (unit$type == "Pow") return(.upy_pow(.qm_get_qsf(obj, unit$args[[1L]]), unit$args[[2L]]))
  v <- obj$qsf[[k]]
  if (!is.null(v)) return(v)
  g <- .units_state$sf_global[[k]]
  if (!is.null(g)) return(.upy_mul(g[[1L]], .qm_get_qsf(obj, g[[2L]])))
  S_One
}

# ---------------------------------------------------------------------------
# DimensionSystem

# a dict given as a units_dict, a named list or vector (string keys), or a
# list of list(key, value) pairs
.as_od <- function(x) {
  if (inherits(x, "units_dict")) return(x)
  if (is.null(x) || !length(x)) return(.ud())
  nm <- names(x)
  if (!is.null(nm) && all(nzchar(nm))) return(.ud(as.list(nm), as.list(x)))
  .ud(lapply(x, `[[`, 1L), lapply(x, `[[`, 2L))
}

.ds_parse_dim <- function(dim) {
  if (is.character(dim)) return(Dimension(Symbol(dim)))
  if (.is_dim(dim)) return(dim)
  if (is_sym(dim) && dim$type == "Symbol") return(Dimension(dim))
  .u_stop("TypeError", paste0(.upy_repr(dim), " wrong type"))
}

.ds_parse_dim_name <- function(dim) {
  if (.is_dim(dim)) return(dim)
  if (is.character(dim)) return(Dimension(Symbol(dim)))
  if (is_sym(dim) && dim$type == "Symbol") return(Dimension(dim))
  .u_stop("TypeError", paste0("unrecognized type ", class(dim)[1L], " for ", .upy_repr(dim)))
}

.u_in <- function(x, lst) any(vapply(lst, function(y) identical(x$key, y$key), TRUE))

#' Dimension systems
#'
#' `DimensionSystem(base_dims, derived_dims, dimensional_dependencies)` is
#' SymPy's coherent set of dimensions: the base dimensions, the derived ones
#' and how the derived ones depend on the base ones. Dimensions may be given
#' as `Dimension` objects or names, and `dimensional_dependencies` as a named
#' list such as `list(velocity = list(length = 1, time = -1))`, a list of
#' `list(dimension, dependencies)` pairs, or a `units_dict`.
#' `get_dimensional_dependencies()` expresses a dimension through the base
#' dimensions as a `units_dict` (fields `k` and `v` hold the keys and the
#' powers, in Python's dict order); `equivalent_dims()` compares two
#' dimensions that way and `is_dimensionless()` tests for no dependencies.
#' `extend()` adds base and derived dimensions. The other functions are the
#' legacy helpers of SymPy (`list_can_dims`, `dim_can_vector`, ...).
#'
#' @param base_dims,new_base_dims Base dimensions.
#' @param derived_dims,new_derived_dims Derived dimensions.
#' @param dimensional_dependencies,new_dim_deps Dependencies of the derived
#'   dimensions.
#' @param dim_sys,x A `DimensionSystem` (or, for `extend()`, `system_dim()`
#'   and `is_consistent()`, a `UnitSystem`).
#' @param name,dim,dim1,dim2,dimension Dimensions or dimension expressions.
#' @param mark_dimensionless Return `{Dimension(1): 1}` for a dimensionless
#'   input.
#' @param ... Further arguments of the `UnitSystem` method.
#' @return A `DimensionSystem`, `units_dict`, list, matrix, expression,
#'   integer or `TRUE`/`FALSE`.
#' @examples
#' u <- physics_units()
#' get_dimensional_dependencies(u$dimsys_SI, u$force)
#' equivalent_dims(u$dimsys_SI, u$energy, u$force * u$length)
#' ds <- DimensionSystem(list(u$length, u$time), list(u$velocity))
#' base_dims(ds)
#' @export
DimensionSystem <- function(base_dims, derived_dims = list(), dimensional_dependencies = list()) {
  .units_ensure()
  deps <- .as_od(dimensional_dependencies)
  base_dims <- lapply(as.list(base_dims), .ds_parse_dim)
  derived_dims <- lapply(as.list(derived_dims), .ds_parse_dim)
  for (dim in base_dims) {
    cur <- .ud_get(deps, dim)
    if (!is.null(cur)) {
      cur <- .as_od(cur)
      one <- .ud_get(cur, dim)
      if (.ud_len(cur) != 1L || is.null(one) || !.u_val_eq(as_sym(one), S_One))
        .u_stop("IndexError", "Repeated value in base dimensions")
    }
    deps <- .ud_set(deps, dim, .ud(list(dim), list(S_One)))
  }
  for (dim in deps$k) {
    dim <- .ds_parse_dim(dim)
    if (!.u_in(dim, derived_dims) && !.u_in(dim, base_dims)) derived_dims[[length(derived_dims) + 1L]] <- dim
  }
  parse_dict <- function(d) {
    d <- .as_od(d)
    .ud(lapply(d$k, .ds_parse_dim_name), lapply(d$v, as_sym))
  }
  nd <- .ud()
  for (i in seq_along(deps$k)) nd <- .ud_set(nd, .ds_parse_dim_name(deps$k[[i]]), parse_dict(deps$v[[i]]))
  for (dim in derived_dims) {
    if (.u_in(dim, base_dims)) .u_stop("ValueError", paste0("Dimension ", .str(dim), " both in base and derived"))
    if (!.ud_has(nd, dim)) nd <- .ud_set(nd, dim, .ud(list(dim), list(S_One)))
  }
  obj <- new.env()
  obj$base_dims <- .key_sort(lapply(base_dims, .sort_key), base_dims)
  obj$derived_dims <- .key_sort(lapply(derived_dims, .sort_key), derived_dims)
  obj$deps <- nd
  class(obj) <- "DimensionSystem"
  .u_mapper_init(obj)
}

#' @export
format.DimensionSystem <- function(x, ...) {
  tup <- function(l) if (length(l) == 1L) paste0("(", .str(l[[1L]]), ",)") else paste0("(", paste(vapply(l, .str, ""), collapse = ", "), ")")
  srt <- function(d) {
    o <- .key_sort(lapply(d$k, .sort_key), seq_along(d$k))
    paste0("{", paste(vapply(o, function(i) paste0(.str(d$k[[i]]), ": ",
                                                    if (inherits(d$v[[i]], "units_dict")) srt(d$v[[i]]) else .str(d$v[[i]])), ""),
                      collapse = ", "), "}")
  }
  paste0("DimensionSystem(", tup(x$base_dims), ", ", tup(x$derived_dims), ", ", srt(x$deps), ")")
}

#' @export
print.DimensionSystem <- function(x, ...) {
  cat(format(x), "\n", sep = "")
  invisible(x)
}

#' @rdname DimensionSystem
#' @export
base_dims <- function(dim_sys) dim_sys$base_dims

#' @rdname DimensionSystem
#' @export
derived_dims <- function(dim_sys) dim_sys$derived_dims

#' @rdname DimensionSystem
#' @export
dimensional_dependencies <- function(dim_sys) {
  d <- dim_sys$deps
  d$v <- lapply(d$v, function(v) structure(v, sorted = TRUE))
  structure(d, sorted = TRUE)
}

.ds_deps_for_name <- function(ds, dimension) {
  if (is.character(dimension)) dimension <- Dimension(Symbol(dimension))
  else if (!.is_dim(dimension)) dimension <- Dimension(dimension)
  nm <- dimension$args[[1L]]
  if (nm$type == "Symbol") {
    v <- .ud_get(ds$deps, dimension)
    if (is.null(v)) return(.ud(list(dimension), list(S_One)))
    return(v)
  }
  if (.u_is_number(nm) || (nm$type == "Constant" && !(nm$name %in% c("I", "oo", "-oo", "zoo", "nan")))) return(.ud())
  get_for_name <- function(z) .ds_deps_for_name(ds, z)
  if (nm$type == "Mul") {
    ret <- .ud()
    for (a in .sympy_args(nm)) {
      d <- get_for_name(a)
      for (i in seq_along(d$k)) ret <- .ud_set(ret, d$k[[i]], .upy_add(.ud_get(ret, d$k[[i]], S_Zero), d$v[[i]]))
    }
    keep <- vapply(ret$v, function(v) !.u_val_eq(v, S_Zero), TRUE)
    return(.ud(ret$k[keep], ret$v[keep]))
  }
  if (nm$type == "Add") {
    dicts <- lapply(.sympy_args(nm), get_for_name)
    if (all(vapply(dicts[-1L], function(d) .ud_eq(d, dicts[[1L]]), TRUE))) return(dicts[[1L]])
    .u_stop("TypeError", "Only equivalent dimensions can be added or subtracted.")
  }
  if (nm$type == "Pow") {
    ex <- nm$args[[2L]]
    dim_base <- get_for_name(nm$args[[1L]])
    dim_exp <- get_for_name(ex)
    if (!.ud_len(dim_exp) || ex$type == "Symbol")
      return(.ud(dim_base$k, lapply(dim_base$v, function(v) .upy_mul(v, ex))))
    .u_stop("TypeError", "The exponent for the power operator must be a Symbol or dimensionless.")
  }
  if (nm$type == "Function") {
    args <- lapply(nm$args, function(a) .dim_from_deps(get_for_name(a)))
    result <- fcall(nm$name, args)
    dicts <- lapply(nm$args, get_for_name)
    if (.is_dim(result)) return(get_dimensional_dependencies(ds, result))
    if (result$type == "Function" && result$name == nm$name) {
      if (nm$name %in% c("sin", "cos", "tan", "cot", "sec", "csc", "sinc")) {
        d0 <- dicts[[1L]]
        if (!.ud_len(d0) || .ud_eq(d0, .ud(list(Dimension("angle")), list(S_One)))) return(.ud())
        .u_stop("TypeError", "The input argument for the function <class 'sympy.physics.units.dimensions.Dimension'> must be dimensionless or have dimensions of angle.")
      }
      if (all(vapply(dicts, function(d) !.ud_len(d), TRUE))) return(.ud())
      .u_stop("TypeError", "The input arguments for the function <class 'sympy.physics.units.dimensions.Dimension'> must be dimensionless.")
    }
    return(get_for_name(result))
  }
  .u_stop("TypeError", paste0("Type ", .u_py_class(nm), " not implemented for get_dimensional_dependencies"))
}

# the repr of type(e) for an expression, as in SymPy's error messages
.u_py_class <- function(e) {
  cls <- .sympy_class(e)
  mod <- switch(cls, Integer = , Rational = , Float = , Zero = , One = , NegativeOne = , Half = "sympy.core.numbers",
                Symbol = "sympy.core.symbol", Dimension = "sympy.physics.units.dimensions",
                Derivative = "sympy.core.function", Integral = "sympy.integrals.integrals", "sympy.core.basic")
  paste0("<class '", mod, ".", cls, "'>")
}

#' @rdname DimensionSystem
#' @export
get_dimensional_dependencies <- function(dim_sys, name, mark_dimensionless = FALSE) {
  .units_ensure()
  dimdep <- .ds_deps_for_name(dim_sys, name)
  if (mark_dimensionless && !.ud_len(dimdep)) return(.ud(list(Dimension(S_One)), list(S_One)))
  dimdep
}

#' @rdname DimensionSystem
#' @export
equivalent_dims <- function(dim_sys, dim1, dim2) {
  if (!inherits(dim_sys, "DimensionSystem"))
    .u_stop("TypeError", "DimensionSystem.equivalent_dims() missing 1 required positional argument: 'dim2'")
  .ud_eq(get_dimensional_dependencies(dim_sys, dim1), get_dimensional_dependencies(dim_sys, dim2))
}

#' @rdname DimensionSystem
#' @export
extend <- function(x, ...) UseMethod("extend")

#' @rdname DimensionSystem
#' @export
extend.DimensionSystem <- function(x, new_base_dims, new_derived_dims = list(), new_dim_deps = NULL, ...) {
  deps <- x$deps
  if (length(new_dim_deps)) {
    nd <- .as_od(new_dim_deps)
    for (i in seq_along(nd$k)) deps <- .ud_set(deps, nd$k[[i]], nd$v[[i]])
  }
  new <- DimensionSystem(c(x$base_dims, as.list(new_base_dims)), c(x$derived_dims, as.list(new_derived_dims)), deps)
  for (k in ls(x$qdm, all.names = TRUE)) assign(k, x$qdm[[k]], envir = new$qdm)
  for (k in ls(x$qsf, all.names = TRUE)) assign(k, x$qsf[[k]], envir = new$qsf)
  new
}

#' @rdname DimensionSystem
#' @export
is_dimensionless <- function(dim_sys, dimension) {
  if (.u_val_eq(dimension$args[[1L]], S_One)) return(TRUE)
  !.ud_len(get_dimensional_dependencies(dim_sys, dimension))
}

#' @rdname DimensionSystem
#' @export
list_can_dims <- function(dim_sys) {
  dimset <- list()
  for (i in dim_sys$base_dims) {
    for (k in get_dimensional_dependencies(dim_sys, i)$k) if (!.u_in(k, dimset)) dimset[[length(dimset) + 1L]] <- k
  }
  if (!length(dimset)) return(dimset)
  dimset[order(vapply(dimset, .str, ""), method = "radix")]
}

#' @rdname DimensionSystem
#' @export
dim_can_vector <- function(dim_sys, dim) {
  deps <- get_dimensional_dependencies(dim_sys, dim)
  Matrix(lapply(list_can_dims(dim_sys), function(d) list(.ud_get(deps, d, S_Zero))))
}

.u_row_join <- function(ms) Reduce(function(a, b) cbind(a, b), ms)

#' @rdname DimensionSystem
#' @export
inv_can_transf_matrix <- function(dim_sys) .u_row_join(lapply(dim_sys$base_dims, function(d) dim_can_vector(dim_sys, d)))

#' @rdname DimensionSystem
#' @export
can_transf_matrix <- function(dim_sys) {
  bd <- dim_sys$base_dims
  bd <- bd[order(vapply(bd, .str, ""), method = "radix")]
  inverse(.u_row_join(lapply(bd, function(d) dim_can_vector(dim_sys, d))))
}

#' @rdname DimensionSystem
#' @export
dim_vector <- function(dim_sys, dim) can_transf_matrix(dim_sys) %*% dim_can_vector(dim_sys, dim)

#' @rdname DimensionSystem
#' @export
print_dim_base <- function(dim_sys, dim) {
  dims <- dim_vector(dim_sys, dim)
  syms <- lapply(dim_sys$base_dims, function(i) if (!is.null(i$extra$symbol)) i$extra$symbol else i$args[[1L]])
  res <- S_One
  for (i in seq_len(min(length(syms), dims$nrow))) res <- .upy_mul(res, .upy_pow(syms[[i]], .mget(dims, i, 1L)))
  res
}

#' @rdname DimensionSystem
#' @export
system_dim <- function(x) length(if (inherits(x, "UnitSystem")) x$base_units else x$base_dims)

#' @rdname DimensionSystem
#' @export
is_consistent <- function(x) {
  if (inherits(x, "UnitSystem")) x <- .us_dimsys(x, "is_consistent")
  M <- inv_can_transf_matrix(x)
  M$nrow == M$ncol
}
