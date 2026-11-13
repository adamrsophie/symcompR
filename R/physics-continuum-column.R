# Axially loaded columns, after
# sympy/physics/continuum_mechanics/column.py.

# shared helpers of the continuum mechanics module

.cm_sympify <- function(v) if (is.null(v)) NULL else as_sym(v)

.cm_same <- function(a, b) {
  if (is.null(a) || is.null(b)) return(is.null(a) && is.null(b))
  identical(as_sym(a)$key, as_sym(b)$key)
}

.cm_tuple_eq <- function(t1, t2) length(t1) == length(t2) && all(vapply(seq_along(t1), function(i) .cm_same(t1[[i]], t2[[i]]), TRUE))

# a support or hinge symbol named after its location, as f'{float(loc):g}'
.cm_loc_name <- function(prefix, loc) {
  if (.expr_is_number(loc)) return(Symbol(paste0(prefix, sprintf("%g", as.numeric(N(loc))))))
  Symbol(paste0(prefix, .str(loc)))
}

# linsolve(eqs, syms).args[0]: the solution values in the order of syms; a
# free unknown stays itself. linsolve() solves over a fraction field, so a
# value prints as its expanded numerator over its expanded denominator.
.cm_linsolve <- function(eqs, syms) {
  sol <- linsolve(eqs, syms)
  if (is.null(sol)) .physics_error("IndexError", "tuple index out of range")
  lapply(syms, function(s) {
    v <- sol[[s$name]]
    if (is.null(v)) return(s)
    nd <- fraction(cancel(v))
    expand(nd[[1L]]) / expand(nd[[2L]])
  })
}

# the field of an object without going through its properties
.cm_get <- function(self, name) get(name, envir = self, inherits = FALSE)
.cm_set <- function(self, name, value) assign(name, value, envir = self)
.cm_push <- function(self, name, value) assign(name, c(.cm_get(self, name), list(value)), envir = self)

# the end load of a distributed load cut at `end`
.cm_taylor <- function(x, value, start, order, end) {
  f <- value * (x - start)^order
  o <- as.integer(as.character(order))
  .add(lapply(0:o, function(i) {
    .div(.mul(list(subs(diff(f, x, i), x, end), SingularityFunction(x, end, i))), factorial(as_sym(i)))
  }))
}

# a base R plot of an expression over [lower, upper]
.cm_plot <- function(e, x, upper, main, ylab, col, n = 200L, lower = 0) {
  xs <- seq(as.numeric(N(as_sym(lower))), as.numeric(N(as_sym(upper))), length.out = n)
  ys <- vapply(xs, function(v) suppressWarnings(as.numeric(N(subs(e, x, v)))), 0)
  graphics::plot(xs, ys, type = "l", col = col, main = main, xlab = "x", ylab = ylab)
  invisible(list(xs, ys))
}

#' Columns under axial loads
#'
#' `Column(length, elastic_modulus, area)` is a column loaded along its
#' axis. `$apply_support(loc)` adds a support, `$apply_load(value, start,
#' order, end)` a load in singularity function form and
#' `$apply_telescope_hinge(loc)` a telescope hinge;
#' `$solve_for_reaction_loads()` finds the reactions, after which
#' `$axial_force()` and `$extension()` give the axial force and the
#' extension along the column. The properties are `$length`,
#' `$elastic_modulus`, `$area`, `$variable`, `$load`, `$applied_loads`,
#' `$reaction_loads` and `$hinge_extensions`. `$plot_axial_force()` and
#' `$plot_extension()` draw base R plots. After
#' `sympy.physics.continuum_mechanics.column`.
#'
#' @param length,elastic_modulus,area The length, the elastic modulus and
#'   the cross section area.
#' @param variable The coordinate along the column.
#' @param base_char The base name of the integration constants.
#' @return A `Column` object.
#' @examples
#' E <- Symbol("E"); A <- Symbol("A")
#' c <- Column(10, E, A)
#' c$apply_support(0)
#' c$apply_load(-1, 5, -1)
#' c$solve_for_reaction_loads()
#' c$reaction_loads
#' @export
Column <- function(length, elastic_modulus, area, variable = Symbol("x"), base_char = "C") {
  self <- .pyenv_new("Column")
  self$length <- length
  self$elastic_modulus <- elastic_modulus
  self$area <- area
  self$variable <- variable
  for (f in c("_bc_extension", "_bc_hinge", "_applied_supports", "_applied_hinges", "_applied_loads"))
    .cm_set(self, f, list())
  .cm_set(self, "_base_char", base_char)
  .cm_set(self, "_reaction_loads", .pydict(list(), list()))
  .cm_set(self, "_hinge_extensions", .pydict(list(), list()))
  .cm_set(self, "_integration_constants", .pydict(list(), list()))
  .cm_set(self, "_load", S_Zero)
  .cm_set(self, "_is_solved", FALSE)
  self
}

.col_apply_load <- function(self, value, start, order, end = NULL) {
  value <- as_sym(value); start <- as_sym(start); order <- as_sym(order)
  x <- self$variable
  if (!is.null(end)) {
    end <- as_sym(end)
    .cm_push(self, "_applied_loads", .pytuple(list(value, start, order, end)))
    ld <- .cm_get(self, "_load") + value * SingularityFunction(x, start, order)
    if (isTRUE(ask_fact(simplify(start - end), "positive"))) ld <- ld + .col_taylor(x, -value, start, order, end)
    else ld <- ld - .col_taylor(x, value, start, order, end)
    .cm_set(self, "_load", ld)
  } else {
    .cm_push(self, "_applied_loads", .pytuple(list(value, start, order, NULL)))
    .cm_set(self, "_load", .cm_get(self, "_load") + value * SingularityFunction(x, start, order))
  }
  invisible(NULL)
}

.col_taylor <- function(x, value, start, order, end) {
  if (isTRUE(ask_fact(order, "negative")))
    .physics_error("ValueError", "If 'end' is provided the 'order' of the load cannot be negative.")
  .cm_taylor(x, value, start, order, end)
}

.col_remove_load <- function(self, value, start, order, end = NULL) {
  value <- as_sym(value); start <- as_sym(start); order <- as_sym(order)
  x <- self$variable
  loads <- .cm_get(self, "_applied_loads")
  key <- list(value, start, order, end)
  i <- which(vapply(loads, function(t) .cm_tuple_eq(unclass(t), key), TRUE))
  if (!length(i)) .physics_error("ValueError", "No such load distribution exists on the column object.")
  .cm_set(self, "_load", .cm_get(self, "_load") - value * SingularityFunction(x, start, order))
  .cm_set(self, "_applied_loads", loads[-i[[1L]]])
  if (!is.null(end)) {
    end <- as_sym(end)
    ld <- .cm_get(self, "_load")
    if (isTRUE(ask_fact(simplify(start - end), "positive"))) ld <- ld - .col_taylor(x, -value, start, order, end)
    else ld <- ld + .col_taylor(x, value, start, order, end)
    .cm_set(self, "_load", ld)
  }
  invisible(NULL)
}

.col_solve <- function(self, ...) {
  x <- self$variable; qx <- .cm_get(self, "_load"); L <- self$length
  A <- self$area; E <- self$elastic_modulus
  applied_supports <- c(list(...), .cm_get(self, "_applied_supports"))
  C_N <- Symbol("C_N"); C_u <- Symbol("C_u")
  axial_force <- -integrate(qx, x) + C_N
  extension <- -(integrate(integrate(qx, x), x)) / (E * A) + C_N * x + C_u
  eqs <- list(limit(axial_force, x, 0L, dir = "-"), limit(axial_force, x, L, dir = "+"))
  for (loc in .cm_get(self, "_bc_extension")) eqs[[length(eqs) + 1L]] <- Eq(subs(extension, x, loc), 0L)
  for (loc in .cm_get(self, "_bc_hinge")) eqs[[length(eqs) + 1L]] <- Eq(limit(axial_force, x, loc, dir = "+"), 0L)
  hinges <- .cm_get(self, "_applied_hinges")
  unknowns <- c(applied_supports, hinges, list(C_N, C_u))
  solution <- lapply(.cm_linsolve(eqs, unknowns), function(s) nsimplify(s, rational = TRUE))
  ns <- length(applied_supports); nh <- length(hinges)
  .cm_set(self, "_reaction_loads", .pydict(applied_supports, solution[seq_len(ns)]))
  .cm_set(self, "_hinge_extensions", .pydict(hinges, solution[ns + seq_len(nh)]))
  .cm_set(self, "_integration_constants", .pydict(list(C_N, C_u), solution[ns + nh + 1:2]))
  .cm_set(self, "_is_solved", TRUE)
  invisible(NULL)
}

.cm_dict_subs <- function(e, d) {
  for (i in seq_along(d$keys)) e <- subs(e, d$keys[[i]], d$values[[i]])
  e
}

.col_solved_load <- function(self)
  .cm_dict_subs(.cm_dict_subs(.cm_get(self, "_load"), .cm_get(self, "_reaction_loads")), .cm_get(self, "_hinge_extensions"))

.col_constant <- function(self, name) {
  if (!isTRUE(.cm_get(self, "_is_solved"))) return(Symbol(name))
  d <- .cm_get(self, "_integration_constants")
  d$values[[match(name, vapply(d$keys, function(k) k$name, ""))]]
}

local({
  sym_setter <- function(field) function(self, v) .cm_set(self, field, as_sym(v))
  .pyclass_def("Column",
    props = list(
      length = function(self) .cm_get(self, "_length"),
      elastic_modulus = function(self) .cm_get(self, "_elastic_modulus"),
      area = function(self) .cm_get(self, "_area"),
      variable = function(self) .cm_get(self, "_variable"),
      reaction_loads = function(self) .cm_get(self, "_reaction_loads"),
      hinge_extensions = function(self) .cm_get(self, "_hinge_extensions"),
      applied_loads = function(self) .cm_get(self, "_applied_loads"),
      load = function(self) .cm_get(self, "_load")),
    setters = list(
      length = sym_setter("_length"), elastic_modulus = sym_setter("_elastic_modulus"), area = sym_setter("_area"),
      variable = function(self, v) {
        if (!(is_sym(v) && v$type == "Symbol")) .physics_error("TypeError", "The variable should be a Symbol object.")
        .cm_set(self, "_variable", v)
      }),
    methods = list(
      `__str__` = function(self) sprintf("Column(%s, %s, %s)", sstr(self$length), sstr(self$elastic_modulus), sstr(self$area)),
      apply_support = function(self, loc) {
        loc <- as_sym(loc)
        if (any(vapply(.cm_get(self, "_applied_supports"), function(s) .cm_same(s, loc), TRUE)))
          .physics_error("ValueError", "There is already a support at this location.")
        r <- .cm_loc_name("R_", loc)
        .cm_push(self, "_applied_supports", r)
        .cm_set(self, "_load", .cm_get(self, "_load") + r * SingularityFunction(self$variable, loc, -1L))
        .cm_push(self, "_bc_extension", loc)
        invisible(NULL)
      },
      apply_load = .col_apply_load,
      remove_load = .col_remove_load,
      apply_telescope_hinge = function(self, loc) {
        loc <- as_sym(loc)
        .cm_push(self, "_bc_hinge", loc)
        u <- .cm_loc_name("u_", loc)
        .cm_push(self, "_applied_hinges", u)
        .cm_set(self, "_load", .cm_get(self, "_load") +
                  (self$elastic_modulus * self$area * u) * SingularityFunction(self$variable, loc, -2L))
        invisible(NULL)
      },
      solve_for_reaction_loads = .col_solve,
      axial_force = function(self) -integrate(.col_solved_load(self), self$variable) + .col_constant(self, "C_N"),
      extension = function(self) {
        x <- self$variable
        -(integrate(integrate(.col_solved_load(self), x), x)) / (self$elastic_modulus * self$area) +
          .col_constant(self, "C_N") * x + .col_constant(self, "C_u")
      },
      plot_axial_force = function(self) .cm_plot(self$axial_force(), self$variable, self$length, "Axial Force", "N(x)", "red"),
      plot_extension = function(self) .cm_plot(self$extension(), self$variable, self$length, "Extension", "u(x)", "blue")))
})
